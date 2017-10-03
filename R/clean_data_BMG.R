#'
#' clean_data_BMG <- function(path) {
#'
#'   # Load Additional Experimental data =========================================
#'   # Create file names for the saved plate reader data
#'
#'
#'   if (path %>% substr(nchar(.) - 2, nchar(.)) %>% tolower() == "csv") { # a single *.csv provided
#'     DataFiles <- basename(path)
#'     path %<>% dirname() %>% paste("/", sep = "")
#'   } else {# a path containing *.csv files of the form "-Ex-Px.csv" or "TRnoX.csv"
#'     DataFiles <- list.files(path) %>% grep("(-E\\d-P\\d|TRno\\d*)\\.(CSV|csv)", ., value = T)
#'   }
#'
#'   PlateInfos <- list() # List of data frames indexed by number
#'   for (i in seq_along(PlateNames)) {
#'     PlateInfos[[i]] <- paste(path, "PlateInfo_", PlateNames[i], ".csv", sep = "")
#'     PlateInfos[[i]] <- read.csv(PlateInfos[[i]], header = TRUE, check.names = T)
#'
#'     # Load in plate data using file names
#'     assign(PlateNames[i], # Makes an object for each plate
#'            read.csv(paste(path, DataFiles[i], sep = ""),
#'                     skip = which(
#'                       grepl("Raw Data \\(\\d{3}\\)", readLines(paste(path, DataFiles[i], sep = "")))
#'                     )[1] - 1,
#'                     header = TRUE,
#'                     check.names = T,
#'                     row.names = NULL
#'            )
#'     )
#'   }
#'
#'   # Creates a list of time points used in each experiment
#'   # This list is based on the the numeric values not the long names (e.g. 1hr - 2hr)
#'   timePoints_list <- list()
#'   for (i in PlateNames) {
#'     timePoints_list[[i]] <- list()
#'     timePoints <- data.frame(get(i)[1,4:ncol(get(i))])
#'     timePoints <- apply(timePoints,2,function(x) gsub('\\.*(\\d*) h ?(\\d*).*', '\\1.0\\2',x))
#'     for (j in seq_along(timePoints)) {
#'       timePoints[j] <- as.numeric(unlist(strsplit(timePoints[j], '\\.'))[1]) + as.numeric(unlist(strsplit(timePoints[j], '\\.'))[2])/60
#'     }
#'     timePoints_list[[i]]$timePoints <- as.numeric(timePoints[as.numeric(timePoints) <= trunctime & is.na(timePoints) == F])
#'   }
#'
#'   # Create a time matrix for each plate
#'   timeMatrix_list <- list()
#'   for (i in names(timePoints_list)) {
#'     timeMatrix_list[[i]] <- timePoints_list[[i]]$timePoints %>%
#'       replicate(nrow(get(i)), .) %>%
#'       t() %>% as.matrix()
#'   }
#'
#'   # Modify the input raw data =================================================
#'   # Remove the timepoint row in the data table
#'   for (i in seq_along(PlateNames)) {
#'     filename <- PlateNames[[i]]
#'     assign(filename, get(filename)[-1,])
#'   }
#'
#'   # Merge Plate Info with growth data based on well
#'   for (i in seq_along(PlateNames)) {
#'     filename <- PlateNames[i]
#'     assign(filename,
#'            merge(
#'              PlateInfos[[i]][ , c("Well.Row", "Well.Col", "Homolog")],
#'              get(filename),
#'              by = c("Well.Row", "Well.Col")
#'            )
#'     )
#'   }
#'   }
