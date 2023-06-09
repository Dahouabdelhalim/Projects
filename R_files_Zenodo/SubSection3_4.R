#This source code reproduces the results shown in subsection
#3.4 "Testing the risk factors"

#author: Andrea Falcon-Cortes

#Functions to compute the linear models
data_lm_rad = function(data, datasub){
  
  d13 = length(which(data$RAD >= 0.5 & data$Year == 2013))/length(which(data$Year == 2013))
  d14 = length(which(data$RAD >= 0.5 & data$Year == 2014))/length(which(data$Year == 2014))
  d15 = length(which(data$RAD >= 0.5 & data$Year == 2015))/length(which(data$Year == 2015))
  d16 = length(which(data$RAD >= 0.5 & data$Year == 2016))/length(which(data$Year == 2016))
  d17 = length(which(data$RAD >= 0.5 & data$Year == 2017))/length(which(data$Year == 2017))
  d18 = length(which(data$RAD >= 0.5 & data$Year == 2018))/length(which(data$Year == 2018))
  d19 = length(which(data$RAD >= 0.5 & data$Year == 2019))/length(which(data$Year == 2019))
  d20 = length(which(data$RAD >= 0.5 & data$Year == 2020))/length(which(data$Year == 2020))
  
  d13.rep = rep(d13, length(which(datasub$Year == 2013)))
  d14.rep = rep(d14, length(which(datasub$Year == 2014)))
  d15.rep = rep(d15, length(which(datasub$Year == 2015)))
  d16.rep = rep(d16, length(which(datasub$Year == 2016)))
  d17.rep = rep(d17, length(which(datasub$Year == 2017)))
  d18.rep = rep(d18, length(which(datasub$Year == 2018)))
  d19.rep = rep(d19, length(which(datasub$Year == 2019)))
  d20.rep = rep(d20, length(which(datasub$Year == 2020)))
  
  data05 = c(d13, d14, d15, d16, d17, d18, d19, d20)
  data05.rep = c(d13.rep, d14.rep, d15.rep, d16.rep,
                 d17.rep, d18.rep, d19.rep, d20.rep)
  
  return(data05.rep)
  
}

data_lm_ptad = function(data, datasub){
  
  d13 = length(which(data$PT.AD == 1 & data$Year == 2013))/length(which(data$Year == 2013))
  d14 = length(which(data$PT.AD == 1 & data$Year == 2014))/length(which(data$Year == 2014))
  d15 = length(which(data$PT.AD == 1 & data$Year == 2015))/length(which(data$Year == 2015))
  d16 = length(which(data$PT.AD == 1 & data$Year == 2016))/length(which(data$Year == 2016))
  d17 = length(which(data$PT.AD == 1 & data$Year == 2017))/length(which(data$Year == 2017))
  d18 = length(which(data$PT.AD == 1 & data$Year == 2018))/length(which(data$Year == 2018))
  d19 = length(which(data$PT.AD == 1 & data$Year == 2019))/length(which(data$Year == 2019))
  d20 = length(which(data$PT.AD == 1 & data$Year == 2020))/length(which(data$Year == 2020))
  
  d13.rep = rep(d13, length(which(datasub$Year == 2013)))
  d14.rep = rep(d14, length(which(datasub$Year == 2014)))
  d15.rep = rep(d15, length(which(datasub$Year == 2015)))
  d16.rep = rep(d16, length(which(datasub$Year == 2016)))
  d17.rep = rep(d17, length(which(datasub$Year == 2017)))
  d18.rep = rep(d18, length(which(datasub$Year == 2018)))
  d19.rep = rep(d19, length(which(datasub$Year == 2019)))
  d20.rep = rep(d20, length(which(datasub$Year == 2020)))
  
  data05 = c(d13, d14, d15, d16, d17, d18, d19, d20)
  data05.rep = c(d13.rep, d14.rep, d15.rep, d16.rep,
                 d17.rep, d18.rep, d19.rep, d20.rep)
  
  return(data05.rep)
  
}

#Necessary packages
library(dplyr)
library(tidyverse)
library(lmtest)
library(latex2exp)

#Data relative to contracts made from 2013 to 2020

#Upload the source list
data <- read.csv("~/RPS_Cont_2013-2020.csv", 
                 stringsAsFactors=FALSE)

#Here we only show an example of the first and second experiments, and the linear model for
#EFOS class and RAD risk factor. The same code works for PCS and NC classes, and the remain
#risk factors.

#The source list only take the EFOS class
efos = data %>% filter(Status == "EFOS")

#First test
#Compute the value in y axis from which the CDF exceeds the threshold RAD=0.5
#i.e., P(RAD>=0.5 | EFOS)
#If this probability is > 0.5, then the first experiment was successful
threshold.y = environment(ecdf(efos$RAD))$y[which(environment(ecdf(efos$RAD))$x >= 0.5)[1]]
threshold.y

#Generate the graph to show the result of test 1
png("Fig5-Left.png", width = 5, height = 5, units = "in", res = 300)
plot(ecdf(efos$RAD), col = "red", verticals = T,col.01line = NULL,
     do.points = T, las =1, pch = 1, xlab = "RAD", ylab = "CDF", main = "EFOS",
     cex.axis = 1.0, cex.lab  =1.25, cex = 0.75)
abline(v = 0.5, col = "grey", lwd = 2, lty = 2)
abline(h = threshold.y, col = "grey", lwd = 2, lty = 2)
dev.off()

#Second test
#Compute the fraction of EFOS in the source list
frac.efos = nrow(efos)/nrow(data)

#Define a set of thresholds
rad05 = data %>% filter (RAD >= 0.5)
rad04 = data %>% filter (RAD >= 0.4)
rad06 = data %>% filter (RAD >= 0.6)

#Compute the fraction of EFOS that exceeded these thresholds, i.e.,
#P(EFOS|RAD >= threshold)
#If this probability is larger than the probability to randomly find an EFOS
#in the source list, then the 2nd test was successful
rad05.efos = length(which(rad05$Status == "EFOS"))/nrow(rad05)
rad04.efos = length(which(rad04$Status == "EFOS"))/nrow(rad04)
rad06.efos = length(which(rad06$Status == "EFOS"))/nrow(rad06)

#Plot the results of the 2nd test
rad = c(0.4, 0.5, 0.6)
rad.efos = c(rad04.efos, rad05.efos, rad06.efos)

png("Fig6-Up.png", units="in", width=5, height=5, res=300)
plot(rad, rad.efos, col = "red", type = "b", pch = 1, lty = 1, log = "y",
     ylim = c(0.001, 0.002), xlab = TeX("$r_i$"),
     ylab = TeX("$P(EFOS \\\\;|\\\\; RAD \\\\geq r_i)$"),
     cex.axis = 1.0, cex.lab  =1.0, cex = 0.75,
     las = 1)
abline(h=frac.efos, col = "red", lty = 2)
dev.off()

#Linear model

#From the source list consider only contracts in PCS and NC classes
fdata.wefos  = data %>% filter(Status == "NC" | Status == "PCS" )

#Consider only weighted data with RAD>= 0.5 (for data without EFOS) and
#PT.AD =1 for EFOS data
rad05.wefos = data_lm_rad(fdata.wefos, efos)
rad05.efos = data_lm_ptad(efos, efos)

#Compute the linear model
ml.efos = lm(rad05.efos ~ rad05.wefos)
summary(ml.efos)

#Define a set of colors for generate the plot of the linear model
col.efos = c(rep("green", length(which(efos$Year == 2013 |
                                         efos$Year == 2014 |
                                         efos$Year == 2015 |
                                         efos$Year == 2016 |
                                         efos$Year == 2017 |
                                         efos$Year == 2018))),
             rep("purple", length(which(efos$Year == 2019 |
                                          efos$Year == 2020))))

#Plot the linear model
png("Fig7.png", width = 5, height = 5, units = "in", res = 300)
plot(rad05.wefos, rad05.efos, col = col.efos, 
     ylab = TeX("$F_{EFOS}(PT.AD=1)$"),
     xlab = TeX("$F_{NC_{U}PCS}(RAD\\\\geq 0.5)$"),cex.axis = 1.0, cex.lab  =1.0, cex = 0.75,
     las = 1)
abline(ml.efos, col = "blue")
legend("bottomright", legend = c("1st Period", "2nd Period"),
       col = c("green", "purple"), pch = c(1,1), bty = "n", cex = 1.0)
dev.off()

#To generate Fig8 the user can use the main idea to plot CI showed in script
#SubSection3_3-NDV.R without separate by classes. 
