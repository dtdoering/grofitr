
options(stringsAsFactors = FALSE)

library("grofit")
library("doBy")
library("RColorBrewer")
library("ggplot2")
library("gdata")
library("magrittr")
library("dplyr")

# Inputs: Enter the proper locations and files for analysis ===================
# Need to set 'ExperimentNames' and 'path' objects before running script

# for (i in seq_along(ExperimentNames)){
  # PlateInfos[i] <- paste(path, "PlateInfo_", ExperimentNames[i], ".csv", sep = "")
# }

truncTime = 40 # hours

# Set output destination of GroFit_df and PDF of plots
out_dest <- "/Users/dtdoering/1_Research/Lab/DATA/phenotyping/plate\ reader/Output/"

# setwd(path)

# Load Additional Experimental data ===========================================
# This data will comprise the first 3 columns of the data tables for grofit

# Create file names for the saved plate reader data
DataFiles <- list.files(path) %>% grep("(-E\\d-P\\d|TRno\\d*)\\.(CSV|csv)", ., value = T)
# Create an array of object names using barcodes found in csv file
PlateNames <- character()
for (i in seq_along(DataFiles)) {
    PlateNames[i] <- grep("ID1: ", readLines(paste(path, DataFiles[i], sep = "")), value = T) %>% substr(6,12) %>% substr(regexpr("[^0]", .), nchar(.))
}

PlateInfos <- list()
for (i in seq_along(PlateNames)) {
    PlateInfos[[i]] <- paste(path, "PlateInfo_", PlateNames[i], ".csv", sep = "")
    PlateInfos[[i]] <- read.csv(PlateInfos[[i]], header = TRUE, check.names = T)

}

# Load in plate data using file names in Experiment_list
for (i in seq_along(PlateNames)) {
  assign(PlateNames[i], # Creates an object name for each plate
         read.csv(paste(path, DataFiles[i], sep = ""),
                  skip = which(
                    grepl("Raw Data \\(600\\)", readLines(paste(path, DataFiles[i], sep = "")))
                  )[1] - 1,
                  header = TRUE,
                  check.names = T,
                  row.names = NULL
         )
  )
}

# Creates a list of time points used in each experiment
# This list is based on the the numeric values not the long names (e.g. 1hr - 2hr)
timePoints_list <- list()
for (i in PlateNames) {
  timePoints_list[[i]] <- list()
  timePoints <- data.frame(get(i)[1,4:ncol(get(i))])
  timePoints <- apply(timePoints,2,function(x) gsub('\\.*(\\d*) h ?(\\d*).*', '\\1.0\\2',x))
  for (j in seq_along(timePoints)) {
    timePoints[j] <- as.numeric(unlist(strsplit(timePoints[j], '\\.'))[1]) + as.numeric(unlist(strsplit(timePoints[j], '\\.'))[2])/60
  }
  timePoints_list[[i]]$timePoints <- as.numeric(timePoints[as.numeric(timePoints) <= truncTime & is.na(timePoints) == F])
}

# Create a time matrix for each plate
timeMatrix_list <- list()
for (i in names(timePoints_list)) {
  timeMatrix_list[[i]] <- timePoints_list[[i]]$timePoints %>%
    replicate(nrow(get(i)), .) %>%
    t() %>% as.matrix()
}

# Modify the input raw data ===================================================
# Remove the timepoint row in the data table
for (i in seq_along(PlateNames)) {
  filename <- PlateNames[[i]]
  assign(filename, get(filename)[-1,])
}

# Merge Plate Info with growth data based on well
for (i in seq_along(PlateNames)) {
  filename <- PlateNames[i]
  assign(filename,
    merge(
      PlateInfos[[i]][ , c("Well.Row", "Well.Col", "Homolog", "feConc", "cuConc")],
      get(filename),
      by = c("Well.Row", "Well.Col")
    )
  )
}

# Drop columns no longer needed from the plate data objects
DropColumns <- c("Well.Row",
                "Well.Col",
                "Content")
for (i in seq_along(PlateNames)) {
  filename <- PlateNames[[i]]
  assign(filename, get(filename)[,!(colnames(get(filename)) %in% DropColumns)])
}

# Run GroFit for all plates ==================================================
if (!exists("GroFitResults")) {
  cat(noquote("GroFitResults not found - creating now"), '\n')
  GroFitResults <- list()
} else {
  cat(noquote("GroFitResults exists - Results will be added"), '\n')
}
cat(noquote("Running GroFit on all plates..."), '\n')
for (i in PlateNames) {
  GroFitResults[[i]] <- list()
  for (j in 1:nrow(get(i))) {
    GroFitResults[[i]][[j]] <- list()
    GroFitResults[[i]][[j]]$Plate <- i

    # Use the 3 columns taken from Plate Info to name corresponding data in GroFitResults
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[1]]] <- get(i)[j,1]
                                            # Should  ^ this be i?
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[2]]] <- get(i)[j,2]
                                            # Should  ^ this be i?
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[3]]] <- get(i)[j,3]
                                            # Should  ^ this be i?

    TimeData <- t(data.frame(timeMatrix_list[[i]][j,]))
    TimeData <- TimeData[, 1:ncol(TimeData), drop = F]

    GrowthData <- t(data.frame(as.numeric(get(i)[j,])))
    GrowthData <- GrowthData[, 1:(ncol(TimeData) + 3), drop = F]

    GroFitResults[[i]][[j]]$GroFitResults = grofit(TimeData, GrowthData,
      control = grofit.control(
        suppress.messages = TRUE,
        fit.opt = "b", # what kind of fit? m=model, s=spline, *b=both
        interactive = FALSE,
        nboot.gc = 100
        # smooth.gc  = 5
      )
    )
  }
  cat(noquote(paste("  Added results of plate ", i, ".", sep = "")), '\n')
}
cat(noquote("  Results complete!"), '\n')

# Creates a data table of lag, growth rate, and saturation ====================
cat(noquote("Summarizing results..."), '\n')
GroFit_df <- data.frame(
  Plate = character(),
  A = character(),
  B = character(),
  C = character(),
  Lag = numeric(),
  GrowthRate = numeric(),
  Saturation = numeric()
  )

colnames(GroFit_df)[2:4] <- c(colnames(get(PlateNames[[1]]))[1], # "Homolog"
                              colnames(get(PlateNames[[1]]))[2], # "feConc"
                              colnames(get(PlateNames[[1]]))[3]) # "cuConc"

for (i in seq_along(GroFitResults)) {
  for (j in seq_along(GroFitResults[[i]])) {
    GroFit_df[j,1] <- GroFitResults[[i]][[j]]$Plate
    GroFit_df[j,2] <- GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[1]]]
    GroFit_df[j,3] <- GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[2]]]
    GroFit_df[j,4] <- GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[3]]]
    GroFit_df[j,5] <- GroFitResults[[i]][[j]]$GroFitResults$gcFit$gcTable$lambda.model
    GroFit_df[j,6] <- GroFitResults[[i]][[j]]$GroFitResults$gcFit$gcTable$mu.model
    GroFit_df[j,7] <- GroFitResults[[i]][[j]]$GroFitResults$gcFit$gcTable$A.model
  }
  assign(paste(names(GroFitResults)[[i]], "_GroFit_df", sep = ""), GroFit_df)
  write.csv(GroFit_df, file = paste(out_dest, names(GroFitResults[i]), "_GroFit_df.csv", sep = ""))
}

# Plot growth curves ======================================================
makeplots <- function(platename) {
  pdf(
    paste(out_dest,
          Sys.time() %>% format("%Y-%m-%d-%H%M"),
          "_",
          platename,
          ".pdf",
          sep = ""
    )
  )
  for (well in seq_along(GroFitResults[[platename]])) {
    times <- GroFitResults[[platename]][[well]]$GroFitResults$gcFit$gcFittedModels[[1]]$raw.time
    od <- GroFitResults[[platename]][[well]]$GroFitResults$gcFit$gcFittedModels[[1]]$raw.data
    temp_df <- data.frame(cbind(times, od))
    colnames(temp_df) = c("Time", "Absorbance")

    Aobs <- summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$A.model
    A.upCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.A.model.up
    A.loCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.A.model.lo

    muobs <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$mu.model
    mu.upCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.mu.model.up
    mu.loCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.mu.model.lo

    lambdaobs <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$lambda.model
    lambda.upCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.lambda.model.up
    lambda.loCI <-  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.lambda.model.lo

    if (is.null(lambdaobs)) {
      lambdaobs <- 0
    }

    # ggplot statements =======================================================
    # growth data points ------------------------------------------------------
    curve <- ggplot(temp_df, aes(x = Time, y = Absorbance))
    curve <- curve + geom_point(pch = 19)

    # saturation point with confidence intervals ------------------------------
    if (!is.null(Aobs) & !is.na(Aobs) & !is.nan(A.upCI) & !is.nan(A.loCI)) {
      curve <- curve +
      geom_hline(
        yintercept = Aobs,
        color = "blue"
      )
    }
    if (!is.nan(A.upCI) & !is.na(A.upCI)) {
      curve <- curve +
      geom_hline(
        yintercept = A.upCI,
        col = "cyan",
        lty = 2
      )
    }
    if (!is.nan(A.loCI) & !is.na(A.loCI)) {
      curve <- curve +
      geom_hline(
        yintercept = A.loCI,
        col = "cyan",
        lty = 2
      )
    }

    # lag time with confidence intervals --------------------------------------
    if (!is.null(lambdaobs) & !is.na(lambdaobs) & !is.nan(lambda.upCI) & !is.nan(lambda.loCI) & 0 < lambdaobs & lambdaobs < truncTime) {
      curve <- curve +
      geom_vline(
        xintercept = lambdaobs,
        color = "green4"
      )
    }
    if (!is.nan(lambda.upCI) & !is.na(lambda.upCI) & 0 < lambda.upCI & lambda.upCI < truncTime) {
      curve <- curve +
      geom_vline(
        xintercept = lambda.upCI,
        col = "green",
        lty = 2
      )
    }
    if (!is.nan(lambda.loCI) & !is.na(lambda.loCI) & 0 < lambda.loCI & lambda.loCI < truncTime) {
      curve <- curve +
      geom_vline(
        xintercept = lambda.loCI,
        col = "green",
        lty = 2
      )
    }

    # max growth rate with confidence intervals -------------------------------
    if (!is.null(muobs) & !is.na(muobs) & !is.nan(mu.upCI) & !is.nan(mu.loCI) & 0 < lambdaobs & lambdaobs < truncTime) {
      curve <- curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = muobs,
        col = "red"
      )
    }
    if (!is.nan(mu.upCI) & !is.na(mu.upCI)) {
      curve <- curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = mu.upCI,
        col = "orange",
        lty = 2
      )
    }
    if (!is.nan(mu.loCI) & !is.na(mu.loCI) & 0 < lambdaobs & lambdaobs < truncTime) {
      curve <- curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = mu.loCI,
        col = "orange",
        lty = 2
      )
    }

    # format plots - axes, max/min, gridlines, etc. ---------------------------
    curve <- curve +
    ylim(0, 2) +
    scale_x_continuous(breaks = seq(0,150,10)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(color = "black")
    ) +
    xlab("Time (h)") +
    ylab("Absorbance") +

    ggtitle(paste(GroFitResults[[platename]][[well]]$Plate, ": ",
      GroFitResults[[platename]][[well]][[colnames(get(PlateNames[[1]]))[1]]],
      " LICM(",
      GroFitResults[[platename]][[well]][[colnames(get(PlateNames[[1]]))[2]]],
      ",",
      GroFitResults[[platename]][[well]][[colnames(get(PlateNames[[1]]))[3]]],
      ")",
      sep = "")
    )
    print(curve)
  }
  dev.off()
}

for (i in names(GroFitResults)) {
  cat(noquote(paste('  Plotting plate ', i, '...', sep = "")))
  makeplots(i)
  cat(noquote('Done.'), '\n')
}

# Cleanup - remove intermediate variables that aren't part of final output ====
rm(list = c("A.loCI", "A.upCI", "Aobs", "curve", "out_dest", "DropColumns", "Experiment_list", "ExperimentName", "ExperimentPlates", "filename", "GroFit_df", "GrowthData", "i", "j", "k", "lambda.loCI", "lambda.upCI", "lambdaobs", "mu.loCI", "mu.upCI", "muobs", "od", "PlateName_list", "PlateNames", "temp_df", "TimeData", "timeMatrix_list", "timePoints", "timePoints_list", "timePoints_m", "TIMES", "times", "truncTime", "usage"))

cat(noquote("GroFit script complete."), '\n')
