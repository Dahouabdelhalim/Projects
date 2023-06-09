#This source code reproduces the results shown in subsection
#3.2 "Comparison between contract classes"
#for Binomial Test section

#author: Andrea Falcon-Cortes

#Function to compute the z-score and p-value 
#for binomial test
z_score_test = function(s1,s2){
  
  x = length(which(s1 == 1))
  y = length(which(s2 == 1))
  
  n1 = length(s1)
  n2 = length(s2)
  
  p1 = x/n1
  p2 = y/n2
  
  p = (x+y)/(n1+n2)
  q = 1-p
  
  z1 = p1-p2
  z2 = p*q
  z3 = (1/n1) + (1/n2)
  
  z = z1/(sqrt(z2*z3))
  
  p_val = 2*pnorm(q=abs(z), lower.tail = F)
  
  zscore = list(z,p_val)
  
  return(zscore)
  
}

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

#To compute all the results shown in section 3.2 (BTest) it
#is necessary to separate the fdata list by periods and
#classes. The user needs to remember that the test is for
#different classes (by pairs) from the same period. 

#Here we show an example for compare EFOS vs PCS in the
#1st period (2013-2018). The same code works for the 
#remaining comparisons. 

#Separate fdata by periods (Only 1st period is shown)
dataEPN=fdata %>% filter(Year==2013 | Year==2014 | 
                         Year==2015 | Year==2016 | 
                         Year==2017 | Year==2018 )
dataEPN = subset(dataEPN,select=-c(Year))

#Separate each period by classes 
#(Only EFOS and PCS classes are shown)
dataEPN.efos = dataEPN %>% filter(Status == "EFOS")
dataEPN.pcs = dataEPN %>% filter(Status == "PCS")

#Run the Binomial test for the period and classes selected

#Array to save the z-score for each dummy variable
ts.z = vector("numeric", length = ncol(dataEPN)-3)
#Array to save the name of each dummy variable
ts.names = vector("character", length = ncol(dataEPN)-3)
#Array to save the p-vale for each dummy variable
ts.pv = vector("numeric", length = ncol(dataEPN)-3)
#Array to save the fraction of contracts with dummy variable=1
#EFOS
ts.f.efos = vector("numeric", length = ncol(dataEPN)-3)
#Array to save the fraction of contracts with dummy variable=1
#PCS
ts.f.pcs = vector("numeric", length = ncol(dataEPN)-3)
#Array to save the difference between fractions
ts.f = vector("numeric", length = ncol(dataEPN)-3)

#Binomial test for each dummy variable in different classes
for(j in 4:ncol(dataEPN)){
  ts.z[j-3] = z_score_test(dataEPN.efos[,j],
                           dataEPN.pcs[,j])[[1]]
  ts.names[j-3] = names(dataEPN)[j]
  ts.pv[j-3] = z_score_test(dataEPN.efos[,j],
                            dataEPN.pcs[,j])[[2]]
  ts.f.efos[j-3] = length(which(dataEPN.efos[,j] == 1))/nrow(dataEPN.efos)
  ts.f.pcs[j-3] = length(which(dataEPN.pcs[,j] == 1))/nrow(dataEPN.pcs)
  ts.f[j-3] = abs(ts.f.efos[j-3] - ts.f.pcs[j-3])
}

#Make a data frame with the obtained results
ts.EPN.e_p = data.frame(ts.names, ts.f.efos, ts.f.pcs,
                        ts.f, abs(ts.z), ts.pv, 
                        stringsAsFactors = F) 
#Arranged from the biggest ts.f to the smaller
ts.EPN.e_p = ts.EPN.e_p %>% arrange(desc(ts.f))
#Take only those dummy variables from which the 
#difference between classes is larger than 10%
#and p_value < 0.05
ts.EPN.ep.diff  = which(ts.EPN.e_p$ts.f >= .1 &
                          ts.EPN.e_p$ts.pv < 0.05)
#Shown the variables that passed the test
ts.EPN.e_p[ts.EPN.ep.diff,]

