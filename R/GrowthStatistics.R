options(stringsAsFactors = FALSE)

library("grofit")
library("doBy")
library("RColorBrewer")
library("ggplot2")
library("gdata")
library("plotrix")

setwd("C:/Users/Dana/OneDrive/Documents/Phenotypic Variation/Analyses/GalactoseAnalyses") # SET TO DIRECTORY WHERE DATA IS

AnalysisName = dir()
#AnalysisName = AnalysisName[-1] # DROPPING REP1 CONTAMINATION
RepNumber = 4 #### ADD THE REPLICATE NUMBER HERE ####

######## Creates the file names for the saved plate reader data ########
# NOTE: the "Experiment Name" was already pulled earlier when looking in the file
# Here you can change the other parts of the saved name if you do not use "E" to denote replicate number or "P" to represent plate number

######## Loads in plate data based on the names created above ########
for(i in 1:length(AnalysisName)){
  filename = paste(AnalysisName[[i]]) # Creates an object name for each plate
  assign(filename, read.csv(AnalysisName[[i]], header = TRUE, check.names = T, row.names = NULL)) # Reads in files from the directory based on the file names created above

}

GrowthData_df = get(AnalysisName[[1]])
for(i in 2:length(AnalysisName)){
  GrowthData_df = rbind(GrowthData_df, get(AnalysisName[[i]]))
}

write.csv(GrowthData_df, file = "C:/Users/Dana/OneDrive/Documents/Phenotypic Variation/Analyses/CombinedGrowthData_df.csv")

Species = unique(GrowthData_df$Species)
Strain = unique(GrowthData_df$Strain)

SummaryStats= data.frame(Species = character(), Strain = character(), Treatment = character(), AvgLag = numeric(), AvgGrowth = numeric(), AvgSaturation = numeric(),
                         VarLag = numeric(), VarGrowth = numeric(), VarSaturation = numeric(),
                         StErLag = numeric(), StErGrowth = numeric(), StErSaturation = numeric(),
                         LagNAs = numeric(), GrowthNAs = numeric(), SaturationNAs = numeric())

LAG = which(colnames(GrowthData_df) == "Lag")
GRATE = which(colnames(GrowthData_df) == "GrowthRate")
SAT = which(colnames(GrowthData_df) == "Saturation")

r = 1

for(s in 1:length(Species)){
  for(i in 1:length(Strain)){
    tempDF=GrowthData_df[which(GrowthData_df$Species == Species[s] & GrowthData_df$Strain == Strain[i]),]
    if(nrow(tempDF) > 0){
      SummaryStats[r, "Species"] = unique(tempDF$Species)
      SummaryStats[r, "Strain"] = unique(tempDF$Strain)
      SummaryStats[r, "Treatment"] = unique(tempDF$Treatment)
      SummaryStats[r, "AvgLag"] = mean(tempDF[,LAG], na.rm = TRUE)
      SummaryStats[r, "AvgGrowth"] = mean(tempDF[,GRATE], na.rm = TRUE)
      SummaryStats[r, "AvgSaturation"] = mean(tempDF[,SAT], na.rm = TRUE)
      SummaryStats[r, "VarLag"] = var(tempDF[,LAG], na.rm = TRUE)
      SummaryStats[r, "VarGrowth"] = var(tempDF[,GRATE], na.rm = TRUE)
      SummaryStats[r, "VarSaturation"] = var(tempDF[,SAT], na.rm = TRUE)
      SummaryStats[r, "StErLag"] = std.error(tempDF[,LAG], na.rm = TRUE)
      SummaryStats[r, "StErGrowth"] = std.error(tempDF[,GRATE], na.rm = TRUE)
      SummaryStats[r, "StErSaturation"] = std.error(tempDF[,SAT], na.rm = TRUE)
      SummaryStats[r, "LagNAs"] = length(which(is.na(tempDF[,LAG]) == TRUE))
      SummaryStats[r, "GrowthNAs"] = length(which(is.na(tempDF[,GRATE]) == TRUE))
      SummaryStats[r, "SaturationNAs"] = length(which(is.na(tempDF[,SAT]) == TRUE))
      r = r+1
    }
  }
}
write.csv(SummaryStats, file = "C:/Users/Dana/OneDrive/Documents/Phenotypic Variation/Analyses/SummaryStats_df.csv", row.names = FALSE)

Growth_AOV = aov(GrowthRate ~ Strain, data = GrowthData_df)
summary(Growth_AOV)
Growth_Tukey = TukeyHSD(Growth_AOV)
write.csv(Growth_Tukey[[1]], file = "GrowthRate_TukeyHSD.csv")



a = ggplot(SummaryStats, aes(y = AvgGrowth, x = Strain)) +
geom_bar(stat = "identity", fill = "white") +
geom_errorbar(aes(ymin=AvgGrowth-StErGrowth, ymax=AvgGrowth+StErGrowth)) +
coord_flip()

Saturation_AOV = aov(Saturation ~ Strain, data = GrowthData_df)
summary(Saturation_AOV)
Saturation_Tukey = TukeyHSD(Saturation_AOV)

Lag_AOV = aov(Lag ~ Strain, data = GrowthData_df)
summary(Lag_AOV)
lag_Tukey = TukeyHSD(Lag_AOV)

#summary(lm(GrowthData_df$Saturation ~ GrowthData_df$GrowthRate))
#plot(GrowthData_df$Saturation ~ GrowthData_df$GrowthRate)
