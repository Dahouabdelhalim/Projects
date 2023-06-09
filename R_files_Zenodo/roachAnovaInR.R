#Basic code for anova analyses performed in R using compiled data across all experimental rounds

setwd('C:/Users/OEB131-G/Dropbox/oeb131analysis/');

data <- read.csv('roachMasterDataForR.csv')
data$idNumber <- as.character(data$idNumber)
#Run separate anovas on each round
r1data <- subset(data, experimentalRound == 1)
r1model <- aov(trackingPerf~idNumber, data=r1data)
summary(r1model)

r2data <- subset(data, experimentalRound == 2)
r2model <- aov(trackingPerf~idNumber, data=r2data)
summary(r2model)

r3data <- subset(data, experimentalRound == 3)
r3model <- aov(trackingPerf~idNumber, data=r3data)
summary(r3model)

r4data <- subset(data, experimentalRound == 4)
r4model <- aov(trackingPerf~idNumber, data=r4data)
summary(r4model)
