#1.Beta regression model fitting and generating predictions of volumetric excess ice content
#Written by Ariane Castagner, finalized in June 2022

#Import packages

library(betareg)
library(stats)
library(Metrics)

#Set working directory (please change to your path to the project folder)

setwd("...//ExcessIceBetaRegression")

#Import data

sentinel <- read.csv("data//SentinelBHdata.csv") # Import Sentinel data
cryostrat <- read.csv("data//CryostratigraphyBHdata.csv") # Import cryostratigraphic data

#Prepare data (remove glaciofluvial samples before fitting model)

sentinel <- sentinel[which(sentinel$Deposit.type != "Glaciofluvial"),]

#Fit the beta regression model using the Sentinel data

fit <- betareg(EITrans ~ Depth + Material.simplified + Deposit.type + Deposit.type.conf + Mid_VI, data = sentinel)

#Make predictions using the Sentinel and cryostratigraphic data and add prediction column to dataset

cryostrat$predicted_EI_beta <- predict(fit, cryostrat)
sentinel$predicted_EI_beta <- predict(fit, sentinel)

#Generate model performance parameters

AIC(fit) #Get AIC value
rmse(sentinel$EITrans, sentinel$predicted_EI_beta) #Get RMSE value
summary <- summary(fit)
summary$pseudo.r.squared #Get R2 value

#Export new dataset with predictions

write.csv(cryostrat, "output//data//CryostratPred.csv")