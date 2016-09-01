usage <- "This script is designed to run with data saved in the following format:
   - Saved as CSV
   - The filename in the format of ExperimentalID-ReplicateNumber-PlateNumber
     for example:
                CGR-E1-P1
              - CGR represents the experimental identifier: Carbon Growth Rate
              - E1 represents the replicate of the experiment - replicate 1
              - P1 represents the plate number - plate 1
This should be run on a single replicate at a time, replicates will be combined after the fact."

options(stringsAsFactors = FALSE)

library("grofit")
library("doBy")
library("RColorBrewer")
library("ggplot2")
library("gdata")
library("magrittr")
library("dplyr")

# Inputs: Enter the proper locations and files for analysis ===================
DataFile =
  "~/1_Research/Lab/DATA/phenotyping/plate\ reader/2016-08-05-Stacker4/R5CZ-E2-P1.csv"
PlateInfo =
  "~/1_Research/Lab/DATA/phenotyping/plate\ reader/2016-08-05-Stacker4/2016-08-01-PlateInfo_R5CZ.csv"

truncTime = 40 # hours

ExperimentNumber = 2 # Add replicate number here
NumberOfPlates= 1 # Add number of plates here

# Set output destination of GroFit_df and PDF of plots
df_dest =
  "/Users/dtdoering/1_Research/Lab/DATA/phenotyping/plate\ reader/Output/R5CZ_GroFit_df.csv"
plot_dest =
  "/Users/dtdoering/1_Research/Lab/DATA/phenotyping/plate\ reader/Output/"

DataLoc = dirname(DataFile)
setwd(DataLoc)

# Load Additional Experimental data ===========================================
# This data will comprise the first 3 columns of the data tables for grofit
PlateInfo = read.csv(PlateInfo, header = TRUE, check.names = T)

# Pull the plate name out of the file name
ExperimentName <- DataFile %>% strsplit("[/.-]") %>% unlist %>% tail(4) %>% head(1)

# Create file names for the saved plate reader data
Experiment_list = list()
for(i in 1:1){
  Experiment_list[[i]] = list()
  for(j in 1:NumberOfPlates){
    Experiment_list[[i]][[j]] = list()
    Experiment_list[[i]][[j]] = paste(ExperimentName,"-E", ExperimentNumber,"-P",j,".csv", sep = "")
  }
}
ExperimentPlates=unlist(Experiment_list)

# Load in plate data using file names in Experiment_list
for(i in seq_along(ExperimentPlates)){
  filename = paste(ExperimentName) # Creates an object name for each plate
  assign(filename,
         read.csv(ExperimentPlates[i],
                  skip = which(
                    grepl("Raw", readLines(DataFile)) |
                    grepl("M\\d+", readLines(DataFile))
                  )[1]-1,
                  header = TRUE,
                  check.names = T,
                  row.names = NULL
         )
  )
}

# Create a list of object names assigned above
# This name will consist of the experiment name name plus a number to represent each plate
PlateName_list = list()
for(i in 1:NumberOfPlates){
  PlateName_list[[i]] = list()
  PlateName_list[[i]] = paste(ExperimentName)#, i, sep = "")
}
PlateNames = unlist(PlateName_list)

# Creates a list of time points used in each experiment
# This list is based on the the numeric values not the long names (e.g. 1hr - 2hr)
timePoints_list = list()
for(i in PlateNames){
  timePoints_list[[i]] = list()
  timePoints = data.frame(get(i)[1,4:ncol(get(i))])
  timePoints = apply(timePoints,2,function(x) gsub('\\.*(\\d*) h ?(\\d*).*', '\\1.0\\2',x))
  for(j in seq_along(timePoints)){
    timePoints[j] <- as.numeric(unlist(strsplit(timePoints[j], '\\.'))[1]) + as.numeric(unlist(strsplit(timePoints[j], '\\.'))[2])/60
  }
  timePoints_list[[i]]$timePoints = as.numeric(timePoints[as.numeric(timePoints) <= truncTime & is.na(timePoints) == F])
}

# Create a time matrix for each plate
timeMatrix_list = list()
for(i in names(timePoints_list)){
  timeMatrix_list[[i]] = list()
  TIMES = rbind(timePoints_list[[i]]$timePoints)
  timePoints_m = data.frame(TIMES)
  for(j in 2:nrow(get(i))){
    timePoints_m[j,] = timePoints_m
  }
  timeMatrix_list[[i]] = as.matrix(timePoints_m)
}

# Modify the input raw data ===================================================
# Remove the timepoint row in the data table
for(i in seq_along(PlateNames)){
  filename = PlateNames[[i]]
  assign(filename, get(filename)[-1,])
}

# Merge PlateInfo with growth data based on well
for(i in seq_along(PlateNames)){
  filename = PlateNames[[i]]
  assign(filename, merge(PlateInfo[ , c("Well.Row", "Well.Col", "Homolog", "feConc", "cuConc")],get(filename), by = c("Well.Row", "Well.Col")))
}

# Drop columns no longer needed
DropColumns = c("Well.Row",
                "Well.Col",
                "Content")
for(i in seq_along(PlateNames)){
  filename = PlateNames[[i]]
  assign(filename, get(filename)[,!(colnames(get(filename)) %in% DropColumns)])
}

# Run GroFit for all plates ==================================================
if(!exists("GroFitResults")){
  cat(noquote("GroFitResults not found - creating now"), '\n')
  GroFitResults = list()
} else {
  cat(noquote("GroFitResults exists - Results will be added"), '\n')
}
cat(noquote("Running GroFit on all plates..."), '\n')
for(i in PlateNames){
  GroFitResults[[i]] = list()
  for(j in 1:nrow(get(i))){
    GroFitResults[[i]][[j]] = list()
    GroFitResults[[i]][[j]]$Plate = i

    # Use the 3 columns taken from PlateInfo to name corresponding data in GroFitResults
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[1]]] = get(i)[j,1]
                                            # Should  ^ this be i?
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[2]]] = get(i)[j,2]
                                            # Should  ^ this be i?
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[3]]] = get(i)[j,3]
                                            # Should  ^ this be i?

    TimeData = t(data.frame(timeMatrix_list[[i]][j,]))
    TimeData = TimeData[, 1:ncol(TimeData), drop = F]

    GrowthData = t(data.frame(as.numeric(get(i)[j,])))
    GrowthData = GrowthData[, 1:(ncol(TimeData)+3), drop = F]

    GroFitResults[[i]][[j]]$GroFitResults = grofit(TimeData, GrowthData,
      control = grofit.control(
        suppress.messages = TRUE,
        fit.opt = "m",
        interactive = FALSE,
        # model.type = c("logistic"),
        nboot.gc = 0,
        # smooth.gc  = 5
        )
    )
  }
  cat(noquote(paste("  Added results of plate ", filename, ".", sep = "")), '\n')
}
cat(noquote("  Results complete!"), '\n')

# Creates a data table of lag, growth rate, and saturation ====================
cat(noquote("Summarizing results..."), '\n')
GroFit_df = data.frame(
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
k = 1
for(i in seq_along(GroFitResults)){
  for(j in seq_along(GroFitResults[[i]])){
    GroFit_df[k,1] = GroFitResults[[i]][[j]]$Plate
    GroFit_df[k,2] = GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[1]]]
    GroFit_df[k,3] = GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[2]]]
    GroFit_df[k,4] = GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[3]]]
    GroFit_df[k,5] = GroFitResults[[i]][[j]]$GroFitResults[["gcFit"]][["gcTable"]][["lambda.model"]]
    GroFit_df[k,6] = GroFitResults[[i]][[j]]$GroFitResults[["gcFit"]][["gcTable"]][["mu.model"]]
    GroFit_df[k,7] = GroFitResults[[i]][[j]]$GroFitResults[["gcFit"]][["gcTable"]][["A.model"]]
    k = 1 + k
  }
  assign(paste(names(GroFitResults)[[i]], "_GroFit_df", sep = ""), GroFit_df)
}
# # Add a column to include the replicate number for the experiment
# GroFit_df$Replicate = rep("Rep2", nrow(GroFit_df))

# Output GroFit_df to CSV
write.csv(GroFit_df, file = df_dest) # unique(GroFit_df$Replicate), "_", Sys.Date(), ".csv", sep = ""), row.names = FALSE)

# Plot growth curves ======================================================

# Function to plot curves for an input plate name -------------------------
# use with for(i in GroFitResults){makeplots(i)}


makeplots <- function(platename){
# -------------------------------------------------------------------------

  pdf(
    paste(plot_dest,
          Sys.time() %>% format("%Y-%m-%d-%H%M"),
          "_",
          platename,
          ".pdf",
          sep = ""
    )
  )

  for(well in seq_along(GroFitResults[[platename]])){
    times = GroFitResults[[platename]][[well]]$GroFitResults$gcFit$gcFittedModels[[1]]$raw.time
    od = GroFitResults[[platename]][[well]]$GroFitResults$gcFit$gcFittedModels[[1]]$raw.data
    temp_df = data.frame(cbind(times, od))
    colnames(temp_df) = c("Time", "Absorbance")

    Aobs = summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$A.model
    A.upCI =  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.A.model.up
    A.loCI =  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.A.model.lo

    muobs =  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$mu.model
    mu.upCI =  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.mu.model.up
    mu.loCI =  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.mu.model.lo

    lambdaobs =  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$lambda.model
    lambda.upCI =  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.lambda.model.up
    lambda.loCI =  summary(GroFitResults[[platename]][[well]]$GroFitResults$gcFit)$ci95.lambda.model.lo

    if(is.null(lambdaobs)){
      lambdaobs = 0
    }

    # ggplot statements ===========================================================
    # growth data points ----------------------------------------------------------
    curve = ggplot(temp_df, aes(x = Time, y = Absorbance))
    curve = curve + geom_point(pch = 19)

    # saturation point with confidence intervals ----------------------------------
    if(!is.null(Aobs) & !is.na(Aobs) & !is.nan(A.upCI) & !is.nan(A.loCI)){
      curve = curve +
      geom_hline(
        yintercept = Aobs,
        color = "blue"
        )
    }
    if(!is.nan(A.upCI) & !is.na(A.upCI)){
      curve = curve +
      geom_hline(
        yintercept = A.upCI,
        col = "cyan",
        lty = 2
        )
    }
    if(!is.nan(A.loCI) & !is.na(A.loCI)){
      curve = curve +
      geom_hline(
        yintercept = A.loCI,
        col = "cyan",
        lty = 2
        )
    }

    # lag time with confidence intervals ------------------------------------------
    if(!is.null(lambdaobs) & !is.na(lambdaobs) & !is.nan(lambda.upCI) & !is.nan(lambda.loCI) & 0 < lambdaobs & lambdaobs < truncTime){
     curve = curve +
     geom_vline(
       xintercept = lambdaobs,
       color = "green4"
       )
    }
    if(!is.nan(lambda.upCI) & !is.na(lambda.upCI) & 0 < lambda.upCI & lambda.upCI < truncTime){
      curve = curve +
      geom_vline(
        xintercept = lambda.upCI,
        col = "green",
        lty = 2
        )
    }
    if(!is.nan(lambda.loCI) & !is.na(lambda.loCI) & 0 < lambda.loCI & lambda.loCI < truncTime){
      curve = curve +
      geom_vline(
        xintercept = lambda.loCI,
        col = "green",
        lty = 2
        )
    }

    # max growth rate with confidence intervals -----------------------------------
    if(!is.null(muobs) & !is.na(muobs) & !is.nan(mu.upCI) & !is.nan(mu.loCI) & 0 < lambdaobs & lambdaobs < truncTime){
      curve = curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = muobs,
        col = "red"
        )
    }
    if(!is.nan(mu.upCI) & !is.na(mu.upCI)){
      curve = curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = mu.upCI,
        col = "orange",
        lty = 2
        )
    }
    if(!is.nan(mu.loCI) & !is.na(mu.loCI) & 0 < lambdaobs & lambdaobs < truncTime){
      curve = curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = mu.loCI,
        col = "orange",
        lty = 2
        )
    }

    # format plots - axes, max/min, gridlines, etc. -------------------------------
    curve = curve +
    ylim(0, 2) +
    scale_x_continuous(breaks = seq(0,150,10)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(color = "black")
      )
    curve = curve +
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

for(i in names(GroFitResults)){
  cat(noquote(paste('  Plotting plate ', i, '...', sep = "")))
  makeplots(i)
  cat(noquote('Done.'), '\n')
}

# Cleanup - remove intermediate variables that aren't part of final output ==
rm(list = c("A.loCI", "A.upCI", "Aobs", "curve", "DataFile", "DataLoc", "df_dest", "DropColumns", "Experiment_list", "ExperimentName", "ExperimentNumber", "ExperimentPlates", "filename", "GroFit_df", "GrowthData", "i", "j", "k", "lambda.loCI", "lambda.upCI", "lambdaobs", "mu.loCI", "mu.upCI", "muobs", "NumberOfPlates", "od", "PlateInfo", "PlateName_list", "PlateNames", "plot_dest", "temp_df", "TimeData", "timeMatrix_list", "timePoints", "timePoints_list", "timePoints_m", "TIMES", "times", "truncTime", "usage"))

cat(noquote("GroFit script complete."), '\n')
