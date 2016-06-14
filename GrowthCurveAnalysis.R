###################################### NOTES ######################################
#                                                                                 #
# This script is designed to run with data saved in the following format:         #
#    - Saved as CSV                                                               #
#    - The filename in the format of ExperimentalID-ReplicateNumber-PlateNumber   #
#      for example:                                                               #
#                 CGR-E1-P1                                                       #
#                                                                                 #
#               - CGR represents the experimental identifier - Carbon Growth Rate #
#               - E1 represents the replicate of the experiment - replicate 1     #
#               - P1 represents the plate number - plate 1                        #
# This should be run on a single replicate at a time, replicates will be combined #
# after the fact.                                                                 #
#                                                                                 #
###################################################################################

options(stringsAsFactors = FALSE)

library("grofit")
library("doBy")
library("RColorBrewer")
library("ggplot2")
library("gdata")
library("magrittr")
library("dplyr")

### INPUTS: Enter the proper locations and files for analysis ###

DataFile = "~/1_Research/Lab/DATA/Plate_Reader/2016-05-25-Stacker/O7ED-E2-P1.csv"
PlateInfo = "~/1_Research/Lab/DATA/Plate_Reader/2016-05-25-Stacker/2016-05-25-PlateInfo_O7ED.csv"
# TreatmentInfo =           "~/1_Research/Lab/DATA/Plate_Reader/TreatmentInfo/20160128_MFC4-E2_TreatmentInfo.csv"

truncTime = 50 # hours

ExperimentNumber = 2 # ADD THE REPLICATE NUMBER HERE
NumberOfPlates= 1 # ADD NUMBER OF PLATES HERE

df_dest = "/Users/dtdoering/1_Research/Lab/DATA/Plate_Reader/Output/O7ED_GroFit_df.csv"
plot_dest = "/Users/dtdoering/1_Research/Lab/DATA/Plate_Reader/Output/"

# Uses parent directory of DataFile as working directory
DataLoc = paste(paste(head(unlist(strsplit(DataFile, "/")), n = length(unlist(strsplit(DataFile, "/"))) - 1), sep = '/', collapse = '/'), "/", sep = '')
setwd(DataLoc)

######## Load Additional Experimental data ########
######## NOTE: This data should be stored in a separate directory from the plate reader data

# This is the data that should be included in the first three columns of the data tables for grofit
# Here: PlateInfo contains the species and lab identifier that is in each well and TreatmentInfo includes the carbon sources in each plate

# In order to do the merge Plate Info should contain a column called "Well.Row" and "Well.Col"
# this will represent whatever is in the that row and column of the plate for example Row A Col 1 contain the yeast yHDO1
ExperimentName <- DataFile %>% strsplit("[/.-]") %>% unlist %>% tail(4) %>% head(1)

PlateInfo = read.csv(PlateInfo, header = TRUE, check.names = T)
# TreatmentInfo = data.frame(Plate = paste(ExperimentName, "", sep = ""), Treatment = "LICM")
# TreatmentInfo = read.csv(TreatmentInfo, header = TRUE, check.names = T)


######## Pulls the beginning of the plate name out of the name ########
# For example, the Experiment name pulled is "CGR" the data saved as CRG-E#-P#
# ExperimentName = dir()[1]   # Dana's method, requires data to be only/first file in working directory

ExperimentName <- DataFile %>% strsplit("[/.-]") %>% unlist %>% tail(4) %>% head(1)

######## Creates the file names for the saved plate reader data ########
# NOTE: the "Experiment Name" was already pulled earlier when looking in the file
# Here you can change the other parts of the saved name if you do not use "E" to denote replicate number or "P" to represent plate number
Experiment_list = list()
for(i in 1:1){
  Experiment_list[[i]] = list()
  for(j in 1:NumberOfPlates){
    Experiment_list[[i]][[j]] = list()
    Experiment_list[[i]][[j]] = paste(ExperimentName,"-E", ExperimentNumber,"-P",j,".csv", sep = "") #Creates the file names to read in the letters "E" and "P" can be changed here
  }
}

ExperimentPlates=unlist(Experiment_list)

######## Loads in plate data based on the names created above ########
for(i in 1:length(ExperimentPlates)){
  filename = paste(ExperimentName)#, i, sep = "") # Creates an object name for each plate
  assign(filename, read.csv(ExperimentPlates[i], skip = 6, header = TRUE, check.names = T, row.names = NULL)) # Reads in files from the directory based on the file names created above

}

######## Creates a list of object names assigned above ########
# This name will consist of the experiment name name plus a number to represent each plate
PlateName_list = list()
for(i in 1:NumberOfPlates){
  PlateName_list[[i]] = list()
  PlateName_list[[i]] = paste(ExperimentName)#, i, sep = "")
}
PlateNames = unlist(PlateName_list)

########  Creates a list of time points used in each experiment ########
# This list is based on the the numeric values not the long names (e.g. 1hr - 2hr)
timePoints_list = list()
for(i in 1:length(PlateNames)){
  timePoints_list[[i]] = list()
  timePoints = data.frame(get(PlateNames[[i]])[1,4:ncol(get(PlateNames[[i]]))])
  timePoints = apply(timePoints,2,function(x) gsub('\\.*(\\d*) h ?(\\d*).*', '\\1.0\\2',x))
  for(j in 1:length(timePoints)){
    timePoints[j] <- as.numeric(unlist(strsplit(timePoints[j], '\\.'))[1]) + as.numeric(unlist(strsplit(timePoints[j], '\\.'))[2])/60
  }
  timePoints_list[[i]]$timePoints = as.numeric(timePoints[as.numeric(timePoints) <= truncTime & is.na(timePoints) == F])
}

######## Removes the time-point row in the data table
for(i in 1:length(PlateNames)){
  filename = PlateNames[[i]]
  assign(filename, get(filename)[-1,])
}

######## Merges the information stored in Plate info with the data tables based on well and column information ########
######## This also drops the well and column information in the original data table ########
for(i in 1:length(PlateNames)){
  filename = PlateNames[[i]]
  assign(filename, merge(PlateInfo[ , c("Well.Row", "Well.Col", "Homolog", "feConc", "cuConc")],get(filename), by = c("Well.Row", "Well.Col")))
}

######## Adds a column of additional information based on Plate number from TreatmentInfo ########
# For example, TreatmentInfo contains carbon source used for each plate
# If you are using different concentrations of a treatment this should contain that information
######## This also deletes the column called "Content"  this is the column that should include concentrations ########
# for(i in 1:length(PlateNames)){
#   filename = PlateNames[[i]]
#   Treat = TreatmentInfo[which(TreatmentInfo$Plate == filename),"Treatment"]
#   Treatment = rep(Treat, nrow(get(filename)))
#   assign(filename, cbind(Treatment, get(filename)))
# }

######## Drops columns no longer needed ########
DropColumns = c("Well.Row",
                "Well.Col",
                "Content")
for(i in 1:length(PlateNames)){
  filename = PlateNames[[i]]
  assign(filename, get(filename)[,-which(colnames(get(filename)) %in% DropColumns)])
}

########  Creates a time matrix for each plate ########
timeMatrix_list = list()
for(i in 1:length(timePoints_list)){
  timeMatrix_list[[i]] = list()
  TIMES = rbind(timePoints_list[[i]]$timePoints)
  timePoints_m = data.frame(TIMES)
  for(j in 2:nrow(get(PlateNames[[i]]))){
    timePoints_m[j,] = timePoints_m
  }
  timeMatrix_list[[i]] = as.matrix(timePoints_m)
}

######## Runs GroFit for all plates ########
# All Names following a dollar sign or in quotes can be modified (see specific lines)
GroFitResults = list()
print(noquote("Generating GroFit results..."))
for(i in 1:length(PlateNames)){
  GroFitResults[[i]] = list()
  for(j in 1:nrow(get(PlateNames[[i]]))){
    GroFitResults[[i]][[j]] = list()
    GroFitResults[[i]][[j]]$Plate = PlateNames[[i]]

    # Use the 3 columns taken from PlateInfo to name corresponding data in GroFitResults
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[1]]] = get(PlateNames[[i]])[j,1]
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[2]]] = get(PlateNames[[i]])[j,2]
    GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[3]]] = get(PlateNames[[i]])[j,3]

    TimeData = t(data.frame(timeMatrix_list[[i]][j,]))
    TimeData = TimeData[, 1:ncol(TimeData), drop = F]

    GrowthData = t(data.frame(as.numeric(get(PlateNames[[i]])[j,])))
    GrowthData = GrowthData[, 1:(ncol(TimeData)+3), drop = F]

    GroFitResults[[i]][[j]]$GroFitResults = grofit(TimeData, GrowthData,
      control = grofit.control(
        suppress.messages = TRUE,
        fit.opt = "m",
        interactive = FALSE,
        model.type = c("logistic"),
        nboot.gc = 0,
        smooth.gc  = 5)
        )
  }
  assign(paste(PlateNames[[i]], "_GroFit_df", sep = ""), GroFit_df)
  print(noquote(paste("  Finished with plate ", i, ".", sep = "")))
}
print(noquote("  Results complete!"))

######## Creates a data table of lag, growth rate, and saturation ########
# Column names can be modified for relevant plate information - "Strain", "Media", "Treatment"
print(noquote("Generating data frame of growth parameters..."))
GroFit_df = data.frame(
  Plate = character(),
  A = character(),
  B = character(),
  C = character(),
  Lag = numeric(),
  GrowthRate = numeric(),
  Saturation = numeric()
  )

colnames(GroFit_df)[2:4] <- c(colnames(get(PlateNames[[1]]))[1],
                              colnames(get(PlateNames[[1]]))[2],
                              colnames(get(PlateNames[[1]]))[3])
k = 1
for(i in 1:length(GroFitResults)){
  for(j in 1:length(GroFitResults[[i]])){
    GroFit_df[k,1] = GroFitResults[[i]][[j]]$Plate
    GroFit_df[k,2] = GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[1]]]
    GroFit_df[k,3] = GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[2]]]
    GroFit_df[k,4] = GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[3]]]
    GroFit_df[k,5] = GroFitResults[[i]][[j]]$GroFitResults[["gcFit"]][["gcTable"]][["lambda.model"]]
    GroFit_df[k,6] = GroFitResults[[i]][[j]]$GroFitResults[["gcFit"]][["gcTable"]][["mu.model"]]
    GroFit_df[k,7] = GroFitResults[[i]][[j]]$GroFitResults[["gcFit"]][["gcTable"]][["A.model"]]
    k = 1 + k
  }
}
######## Adds a column to include the replicate number for the experiment this should be changed for each replicate ########
# GroFit_df$Replicate = rep("Rep2", nrow(GroFit_df)) ### CHANGE THE NAME FOR REPLICATE NUMBER

# SET DIRECTORY to SAVE GroFit_df TO
write.csv(GroFit_df, file = df_dest) # unique(GroFit_df$Replicate), "_", Sys.Date(), ".csv", sep = ""), row.names = FALSE)

######## Creates a growth curve ########
# SET DIRECTORY TO SAVE PDF OUTPUT TO #
pdf(paste(plot_dest, format(Sys.time(), "%Y-%m-%d-%H%M"), "_", noquote(unique(GroFit_df$Plate)),".pdf", sep = "", collapse = "+"))
print(noquote("Plotting curves..."))

for(i in 1:length(GroFitResults)){
  for(j in 1:length(GroFitResults[[i]])){
    times = GroFitResults[[i]][[j]]$GroFitResults$gcFit$gcFittedModels[[1]]$raw.time
    od = GroFitResults[[i]][[j]]$GroFitResults$gcFit$gcFittedModels[[1]]$raw.data
    temp_df = data.frame(cbind(times, od))
    colnames(temp_df) = c("Time", "Absorbance")

    Aobs = summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$A.model
    A.upCI =  summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$ci95.A.model.up
    A.loCI =  summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$ci95.A.model.lo

    muobs =  summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$mu.model
    mu.upCI =  summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$ci95.mu.model.up
    mu.loCI =  summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$ci95.mu.model.lo

    lambdaobs =  summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$lambda.model
    lambda.upCI =  summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$ci95.lambda.model.up
    lambda.loCI =  summary(GroFitResults[[i]][[j]]$GroFitResults$gcFit)$ci95.lambda.model.lo

    if(is.null(lambdaobs)){
      lambdaobs = 0
    }
    # NOTE: If you changed the list names in GroFitResults they should also be changed here
    curve = ggplot(temp_df, aes(x = Time, y = Absorbance))
    curve = curve + geom_point(pch = 19)

	# Saturation point - horizontal line with confidence intervals
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

	# Lag time - vertical line with confidence intervals
    if(!is.null(lambdaobs) & !is.na(lambdaobs) & !is.nan(lambda.upCI) & !is.nan(lambda.loCI)){
     curve = curve +
     geom_vline(
       xintercept = lambdaobs,
       color = "green4"
       )
    }
    if(!is.nan(lambda.upCI) & !is.na(lambda.upCI)){
      curve = curve +
      geom_vline(
        xintercept = lambda.upCI,
        col = "green",
        lty = 2
        )
    }
    if(!is.nan(lambda.loCI) & !is.na(lambda.loCI)){
      curve = curve +
      geom_vline(
        xintercept = lambda.loCI,
        col = "green",
        lty = 2
        )
    }

	# Max growth rate - sloped line with confidence intervals
    if(!is.null(muobs) & !is.na(muobs) & !is.nan(mu.upCI) & !is.nan(mu.loCI)){
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
    if(!is.nan(mu.loCI) & !is.na(mu.loCI)){
      curve = curve +
      geom_abline(
        intercept = -(muobs*lambdaobs),
        slope = mu.loCI,
        col = "orange",
        lty = 2
        )
    }

    # Plot formatting - axes, max/min, gridlines, etc.
    curve = curve +
    ylim(0,2) +
    scale_x_continuous(breaks = seq(0,150,10)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
      )
    curve = curve +
    xlab("Time") +
    ylab("Absorbance") +

    ggtitle(paste(GroFitResults[[i]][[j]]$Plate, ": ",
      GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[1]]],
      " LICM(",
      GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[2]]],
      ",",
      GroFitResults[[i]][[j]][[colnames(get(PlateNames[[1]]))[3]]],
      ")",
      sep = "")
      )
    print(curve)
    # print(noquote(paste("  Plate ", i, ", well ", j, " plotted.", sep = "")))

  #------------Operations on Well-by-well basis must be done above here----------#
  }

  print(noquote(paste("Plate ", i, " finished.", sep = "")))

  #------------Operations on plate-by-plate basis must be done above here--------#
}
dev.off()

print(noquote("GroFit script complete."))
