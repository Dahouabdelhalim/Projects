#This source code reproduces the results shown in subsection
#3.3 "Comparison between government periods"
#for dummy variables

#author: Andrea Falcon-Cortes

#Necessary packages
library(dplyr)
library(tidyverse)

#Data relative to contracts made from 2013 to 2020

#Upload the source list
data <- read.csv("~/RPS_Cont_2013-2020.csv", 
                 stringsAsFactors=FALSE)

#From the source list take only the dummy variables
fdata = subset(data, select = c(Status, DepID, ProvContID,
                                GO.APF, GO.GE, GO.GM,
                                PC.N, PC.I, PC.ITLC,
                                CT.OP, CT.S, CT.ADQ,
                                CT.AR, CT.SLAOP,
                                PT.AD, PT.I3P, PT.LP,
                                S.NOM, S.MED, S.PEQ,
                                S.MIC, S.NA, Year))

#To make the boxplots related to the comparison between
#government periods inside each class, it is necessary 
#to separate for years and classes. Here we only show and example
#for EFOS class. The same codes works for PCS and NC classes.

dataEPN_13.efos=fdata %>% filter(Year==2013 & Status == "EFOS")
dataEPN_14.efos=fdata %>% filter(Year==2014 & Status == "EFOS")
dataEPN_15.efos=fdata %>% filter(Year==2015 & Status == "EFOS")
dataEPN_16.efos=fdata %>% filter(Year==2016 & Status == "EFOS")
dataEPN_17.efos=fdata %>% filter(Year==2017 & Status == "EFOS")
dataEPN_18.efos=fdata %>% filter(Year==2018 & Status == "EFOS")

dataEPN_13.efos = subset(dataEPN_13.efos,select=-c(Year))
dataEPN_14.efos = subset(dataEPN_14.efos,select=-c(Year))
dataEPN_15.efos = subset(dataEPN_15.efos,select=-c(Year))
dataEPN_16.efos = subset(dataEPN_16.efos,select=-c(Year))
dataEPN_17.efos = subset(dataEPN_17.efos,select=-c(Year))
dataEPN_18.efos = subset(dataEPN_18.efos,select=-c(Year))

dataAMLO_19.efos=fdata %>% filter(Year==2019 & Status == "EFOS")
dataAMLO_20.efos=fdata %>% filter(Year==2020 & Status == "EFOS")

dataAMLO_19.efos = subset(dataAMLO_19.efos,select=-c(Year))
dataAMLO_20.efos = subset(dataAMLO_20.efos,select=-c(Year))


#Generate the boxplots
for(i in 4:22){
  
  l.main = paste("",sep = "")
  
  #Compute the fraction of contracts with dummy variable = 1
  fc1_13 = length(which(dataEPN_13.efos[,i]==1))/length(dataEPN_13.efos[,i])
  fc1_14 = length(which(dataEPN_14.efos[,i]==1))/length(dataEPN_14.efos[,i])
  fc1_15 = length(which(dataEPN_15.efos[,i]==1))/length(dataEPN_15.efos[,i])
  fc1_16 = length(which(dataEPN_16.efos[,i]==1))/length(dataEPN_16.efos[,i])
  fc1_17 = length(which(dataEPN_17.efos[,i]==1))/length(dataEPN_17.efos[,i])
  fc1_18 = length(which(dataEPN_18.efos[,i]==1))/length(dataEPN_18.efos[,i])
  fc1_19 = length(which(dataAMLO_19.efos[,i]==1))/length(dataAMLO_19.efos[,i])
  fc1_20 = length(which(dataAMLO_20.efos[,i]==1))/length(dataAMLO_20.efos[,i])
  
  fc1 = c(fc1_13, fc1_14, fc1_15, fc1_16, fc1_17, fc1_18)
  
  #range for y axis
  y.max = max(fc1, fc1_19, fc1_20)
  y.min = min(fc1, fc1_19, fc1_20)

  #define colors
  c1 <- rgb(0.4,0.8,0)
  c2 <- rgb(0.4,0.8,0, alpha = 0.2)
  c3 <- rgb(0.4,0.8,0, alpha = 0.6)
  
  name.g = paste("EFOS-",names(dataEPN_13.efos)[i],"-BP.png", sep = "")
  png(name.g, width = 5, height = 5, units = "in", res = 300)
  
  par(xpd = T)
  
  boxplot(fc1, names=names(dataEPN_13.efos)[i], ylab = "Fraction",
          col = c2, whiskcol = c1, staplecol = c3, medcol = "orange", 
          boxcol = c3, outcol = c3,
          pch = 0, ylim = c(y.min,y.max),
          xlab = names(dataEPN_13.efos)[i], main = l.main,
          cex.main = 1.25, cex.lab = 1.25, cex.axis = 1.0, cex.names = 1.25, las = 1)
  
  points(1,fc1_19, pch = 1, col = adjustcolor("purple", alpha.f = 0.4), cex = 1.25, lwd = 2)
  points(1,fc1_20, pch = 2, col = "purple", cex = 1.25, lwd = 2)
  
  dev.off()
  
}
