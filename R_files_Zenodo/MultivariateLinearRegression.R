#This source code reproduces the results 
#of multivariate (MV) linear models.

#author: Andrea Falcon-Cortes

#Functions to compute the MV linear models

frac.conf.dummy = function(data, variable, value, datasub){
  d13 = length(which(data[,variable] == value & data$Year == 2013))/length(which(data$Year == 2013))
  d14 = length(which(data[,variable] == value & data$Year == 2014))/length(which(data$Year == 2014))
  d15 = length(which(data[,variable] == value & data$Year == 2015))/length(which(data$Year == 2015))
  d16 = length(which(data[,variable] == value & data$Year == 2016))/length(which(data$Year == 2016))
  d17 = length(which(data[,variable] == value & data$Year == 2017))/length(which(data$Year == 2017))
  d18 = length(which(data[,variable] == value & data$Year == 2018))/length(which(data$Year == 2018))
  d19 = length(which(data[,variable] == value & data$Year == 2019))/length(which(data$Year == 2019))
  d20 = length(which(data[,variable] == value & data$Year == 2020))/length(which(data$Year == 2020))
  
  d13.rep = rep(d13, length(which(datasub$Year == 2013)))
  d14.rep = rep(d14, length(which(datasub$Year == 2014)))
  d15.rep = rep(d15, length(which(datasub$Year == 2015)))
  d16.rep = rep(d16, length(which(datasub$Year == 2016)))
  d17.rep = rep(d17, length(which(datasub$Year == 2017)))
  d18.rep = rep(d18, length(which(datasub$Year == 2018)))
  d19.rep = rep(d19, length(which(datasub$Year == 2019)))
  d20.rep = rep(d20, length(which(datasub$Year == 2020)))
  
  d = c(d13.rep, d14.rep, d15.rep, d16.rep,
        d18.rep, d17.rep, d19.rep, d20.rep)
  
  #The independent and dependent variables are standarized 
  #to compare the coefficients of independent variables in
  #MV model. 
  
  d = (d-mean(d))/sd(d)
  
  return(d)
}

frac.conf.rf = function(data, variable, value, datasub){
  d13 = length(which(data[,variable] >= value & data$Year == 2013))/length(which(data$Year == 2013))
  d14 = length(which(data[,variable] >= value & data$Year == 2014))/length(which(data$Year == 2014))
  d15 = length(which(data[,variable] >= value & data$Year == 2015))/length(which(data$Year == 2015))
  d16 = length(which(data[,variable] >= value & data$Year == 2016))/length(which(data$Year == 2016))
  d17 = length(which(data[,variable] >= value & data$Year == 2017))/length(which(data$Year == 2017))
  d18 = length(which(data[,variable] >= value & data$Year == 2018))/length(which(data$Year == 2018))
  d19 = length(which(data[,variable] >= value & data$Year == 2019))/length(which(data$Year == 2019))
  d20 = length(which(data[,variable] >= value & data$Year == 2020))/length(which(data$Year == 2020))
  
  d13.rep = rep(d13, length(which(datasub$Year == 2013)))
  d14.rep = rep(d14, length(which(datasub$Year == 2014)))
  d15.rep = rep(d15, length(which(datasub$Year == 2015)))
  d16.rep = rep(d16, length(which(datasub$Year == 2016)))
  d17.rep = rep(d17, length(which(datasub$Year == 2017)))
  d18.rep = rep(d18, length(which(datasub$Year == 2018)))
  d19.rep = rep(d19, length(which(datasub$Year == 2019)))
  d20.rep = rep(d20, length(which(datasub$Year == 2020)))
  
  d = c(d13.rep, d14.rep, d15.rep, d16.rep,
        d18.rep, d17.rep, d19.rep, d20.rep)
  
  d = (d-mean(d))/sd(d)
  
  return(d)
}

frac.conf.wc = function(data, variable1, value1, variable2, value2, datasub){
  d13 = length(which(data[,variable1] >= value1 & data[,variable2] >= value2 & data$Year == 2013))/length(which(data$Year == 2013))
  d14 = length(which(data[,variable1] >= value1 & data[,variable2] >= value2 & data$Year == 2014))/length(which(data$Year == 2014))
  d15 = length(which(data[,variable1] >= value1 & data[,variable2] >= value2 & data$Year == 2015))/length(which(data$Year == 2015))
  d16 = length(which(data[,variable1] >= value1 & data[,variable2] >= value2 & data$Year == 2016))/length(which(data$Year == 2016))
  d17 = length(which(data[,variable1] >= value1 & data[,variable2] >= value2 & data$Year == 2017))/length(which(data$Year == 2017))
  d18 = length(which(data[,variable1] >= value1 & data[,variable2] >= value2 & data$Year == 2018))/length(which(data$Year == 2018))
  d19 = length(which(data[,variable1] >= value1 & data[,variable2] >= value2 & data$Year == 2019))/length(which(data$Year == 2019))
  d20 = length(which(data[,variable1] >= value1 & data$Year == 2020))/length(which(data$Year == 2020))
  
  d13.rep = rep(d13, length(which(datasub$Year == 2013)))
  d14.rep = rep(d14, length(which(datasub$Year == 2014)))
  d15.rep = rep(d15, length(which(datasub$Year == 2015)))
  d16.rep = rep(d16, length(which(datasub$Year == 2016)))
  d17.rep = rep(d17, length(which(datasub$Year == 2017)))
  d18.rep = rep(d18, length(which(datasub$Year == 2018)))
  d19.rep = rep(d19, length(which(datasub$Year == 2019)))
  d20.rep = rep(d20, length(which(datasub$Year == 2020)))
  
  d = c(d13.rep, d14.rep, d15.rep, d16.rep,
        d18.rep, d17.rep, d19.rep, d20.rep)
  
  d = (d-mean(d))/sd(d)
  
  return(d)
}

#Necessary packages
library(dplyr)
library(tidyverse)
library(car)

#Data relative to contracts made from 2013 to 2020

#Upload the source list
data <- read.csv("~/RPS_Cont_2013-2020.csv", 
                 stringsAsFactors=FALSE)

#Here we show an example of MV models from EFOS class. 

#Obtain EFOS class from the source list
#This will be used to compute the dependent variable and to weight the model
efos = data %>% filter(Status == "EFOS")

#Take the source list except EFOS class (to avoid autocorrelation)
#This will be used to compute the independent variables
fdata.wefos  = data %>% filter(Status == "NC" | Status == "PCS" )

#Compute of dependent variable (ptad.efos - fraction of single-bidder contracts
#in EFOS class). And independet variable (rad05.wefos - fraction of contracts
#which fraction of single-bidder between buyer and supplier is larger than 0.5).
#This independent variable will be in all MV models.
rad05.wefos = frac.conf.rf(fdata.wefos, "RAD", 0.5, efos)
ptad.efos = frac.conf.dummy(efos, "PT.AD",1, efos)


############ MV model ########################

#Define the remain independent variables
name.iv = c("GO.APF","GO.GE","GO.GM",
            "PC.N", "PC.I", "PC.ITLC",
            "CT.OP", "CT.SLAOP", "CT.S", "CT.ADQ", "CT.AR",
            "PT.AD","PT.I3P","PT.LP",
            "S.NOM","S.MED","S.PEQ","S.MIC","S.NA",
            "CPW","SPW","Fav")

#Lists to save the results of MV models
multi.ml.efos = vector("list", length = length(name.iv))
names(multi.ml.efos) = name.iv
vf.efos = vector("list", length = length(name.iv))
names(vf.efos) = name.iv
extra.iv.efos = vector("list", length = length(name.iv))
names(extra.iv.efos) = name.iv

#Compute the MV models
for(i in 1:length(name.iv)){
  if(i>=1 & i<=19){
    #Compute the 2nd independent variable
    extra.iv.efos[[i]] = frac.conf.dummy(fdata.wefos, name.iv[i], 1, efos)
    #Compute the MV model
    ml = lm(ptad.efos ~rad05.wefos + extra.iv.efos[[i]])
    #Save the results
    multi.ml.efos[[i]] = summary(ml)
    vf.efos[[i]] = car::vif(ml)
  } else if (i==20){
    extra.iv.efos[[i]] = frac.conf.wc(fdata.wefos, name.iv[i], 5, name.iv[i+1],10000, efos)
    ml = lm(ptad.efos ~rad05.wefos + extra.iv.efos[[i]])
    multi.ml.efos[[i]] = summary(ml)
    vf.efos[[i]] = car::vif(ml)
  } else if (i == 21){
    extra.iv.efos[[i]] =extra.iv.efos[[i-1]]
    multi.ml.efos[[i]] = multi.ml.efos[[i-1]]
    vf.efos[[i]] = vf.efos[[i-1]]
  } else{
    extra.iv.efos[[i]] = frac.conf.rf(fdata.wefos, name.iv[i], 0.9, efos)
    ml = lm(ptad.efos ~rad05.wefos + extra.iv.efos[[i]])
    multi.ml.efos[[i]] = summary(ml)
    vf.efos[[i]] = car::vif(ml)
  }
}

