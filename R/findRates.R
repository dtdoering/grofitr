#'
#'
#'
#'

options(stringsAsFactors = FALSE, warn = 0)

findRates <- function(path,
                      trunctime = 1000,
                      out_dest = "/Users/dtdoering/1_Research/DATA/phenotyping/plate\ reader/Output/") {
  if (missing(path)) {
    stop("'path' is undefined.")
  }

  # Load Additional Experimental data =========================================
  # Create file names for the saved plate reader data
  if (path %>% substr(nchar(.) - 2, nchar(.)) %>% tolower() == "csv") {
    DataFiles <- basename(path)
    path %<>% dirname() %>% paste("/", sep = "")
  } else {
    DataFiles <- list.files(path) %>% grep("(-E\\d-P\\d|TRno\\d*)\\.(CSV|csv)", ., value = T)
  }

  PlateInfos <- list() # List of data frames indexed by number
  for (i in seq_along(PlateNames)) {
      PlateInfos[[i]] <- paste(path, "PlateInfo_", PlateNames[i], ".csv", sep = "")
      PlateInfos[[i]] <- read.csv(PlateInfos[[i]], header = TRUE, check.names = T)

    # Load in plate data using file names
    assign(PlateNames[i], # Makes an object for each plate
           read.csv(paste(path, DataFiles[i], sep = ""),
                    skip = which(
                      grepl("Raw Data \\(\\d{3}\\)", readLines(paste(path, DataFiles[i], sep = "")))
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
    timePoints_list[[i]]$timePoints <- as.numeric(timePoints[as.numeric(timePoints) <= trunctime & is.na(timePoints) == F])
  }

  # Create a time matrix for each plate
  timeMatrix_list <- list()
  for (i in names(timePoints_list)) {
    timeMatrix_list[[i]] <- timePoints_list[[i]]$timePoints %>%
      replicate(nrow(get(i)), .) %>%
      t() %>% as.matrix()
  }

  # Modify the input raw data =================================================
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
        PlateInfos[[i]][ , c("Well.Row", "Well.Col", "Homolog")],
        get(filename),
        by = c("Well.Row", "Well.Col")
      )
    )
  }

  # Run GroFit for all plates =================================================
  if (!exists("GroFitResults")) {
    cat(noquote("GroFitResults not found - creating now"), '\n')
    GroFitResults <- list()
  } else {
    cat(noquote("GroFitResults exists - Results will be added"), '\n')
  }
  cat(noquote(paste("Running GroFit on", length(PlateNames), "plates...")), '\n')
  for (i in seq_along(PlateNames)) {
    GroFitResults[[PlateNames[i]]] <- list()
    for (j in 1:nrow(get(PlateNames[i]))) {
      GroFitResults[[PlateNames[i]]][[j]] <- list()
      GroFitResults[[PlateNames[i]]][[j]]$Plate <- PlateNames[i]

      # Use the 3 columns taken from Plate Info to name corresponding data in GroFitResults
      GroFitResults[[PlateNames[i]]][[j]][[colnames(get(PlateNames[[i]]))[1]]] <- get(PlateNames[i])[j,1]
                                              # Should  ^ this be i?
      GroFitResults[[PlateNames[i]]][[j]][[colnames(get(PlateNames[[i]]))[2]]] <- get(PlateNames[i])[j,2]
                                              # Should  ^ this be i?
      GroFitResults[[PlateNames[i]]][[j]][[colnames(get(PlateNames[[i]]))[3]]] <- get(PlateNames[i])[j,3]
                                              # Should  ^ this be i?

      TimeData <- t(data.frame(timeMatrix_list[[PlateNames[i]]][j,]))
      TimeData <- TimeData[, 1:ncol(TimeData), drop = F]

      GrowthData <- t(data.frame(as.numeric(get(PlateNames[i])[j,])))
      GrowthData <- GrowthData[, 2:(ncol(TimeData) + 4), drop = F]
      GroFitResults[[PlateNames[i]]][[j]]$GroFitResults = grofit(TimeData, GrowthData,
        control = grofit.control(
          suppress.messages = TRUE,
          fit.opt = "b", # what kind of fit? m=model, s=spline, *b=both
          interactive = FALSE,
          nboot.gc = 100
          # smooth.gc  = 5
        )
      )
      cat(j, " ")
    }
    cat(noquote(paste("  Added results of plate ", PlateNames[i], " (", i, "/", length(PlateNames), ")", sep = "")), '\n')
  }
  cat(noquote("  Results complete!"), '\n')

  # Creates a data table of lag, growth rate, and saturation ==================
  cat(noquote("Summarizing results..."))
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

  # Store growth parameters from GroFitResults in X_GroFit_df
  for (plate in PlateNames) {
    for (well in seq_along(GroFitResults[[plate]])) {
      GroFit_df[well,1] <- GroFitResults[[plate]][[well]]$Plate
      GroFit_df[well,2] <- GroFitResults[[plate]][[well]][[colnames(get(PlateNames[[1]]))[1]]]
      GroFit_df[well,3] <- GroFitResults[[plate]][[well]][[colnames(get(PlateNames[[1]]))[2]]]
      GroFit_df[well,4] <- GroFitResults[[plate]][[well]][[colnames(get(PlateNames[[1]]))[3]]]
      GroFit_df[well,5] <- GroFitResults[[plate]][[well]]$GroFitResults$gcFit$gcTable$lambda.model
      GroFit_df[well,6] <- GroFitResults[[plate]][[well]]$GroFitResults$gcFit$gcTable$mu.model
      GroFit_df[well,7] <- GroFitResults[[plate]][[well]]$GroFitResults$gcFit$gcTable$A.model
    }
    assign(paste(plate, "_GroFit_df", sep = ""), GroFit_df, envir = .GlobalEnv)
    write.csv(GroFit_df, file = paste(out_dest, plate, "_GroFit_df.csv", sep = ""))
  }
  assign("GroFitResults", GroFitResults, envir = .GlobalEnv)
  cat(noquote('Done.'), '\n')

}

cat(noquote("GroFit script sourced."), '\n')
cat(noquote("Functions available:"))
cat(noquote("
1. getPlateNames(path)"), '\n')
cat(noquote("
2. findRates(path,
             trunctime = 1000,
             out_dest = \"/Users/dtdoering/1_Research/Lab/DATA/phenotyping/plate\ reader/Output/\")"), '\n')
cat(noquote("
3. makePlots(platename,
             trunctime = 1000,
             out_dest = \"/Users/dtdoering/1_Research/Lab/DATA/phenotyping/plate\ reader/Output/\")"), '\n')
