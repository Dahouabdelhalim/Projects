rm(list = ls())
source("src/Rfunctions.R")


data <- read.csv("sister_data_compact.csv", header = TRUE, na.strings = ".")

## correlations (hybridization and diversification)
cor.test(data$HybPropdiff, data$mom_diff, method = "pearson")
cor.test(data$HybPropdiff, data$medusa_diff, method = "pearson") 
cor.test(data$HybPropdiff, data$bamm_diff, method = "pearson") 

cor.test(data$HybRatiodiff, data$mom_diff, method = "pearson")
cor.test(data$HybRatiodiff, data$medusa_diff, method = "pearson")
cor.test(data$HybRatiodiff, data$bamm_diff, method = "pearson") #


## correlations (hybridization and richness)
cor.test(data$HybPropdiff, log(data$RichnessDiff), method = "pearson") 
cor.test(data$HybRatiodiff, log(data$RichnessDiff), method = "pearson") 