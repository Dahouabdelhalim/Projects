rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(ggplot2)

#############################
###Nitrogen Concentrations###
#############################

#import soil nitrogen concentration data
soilConc <- read.csv("nohe141clean.csv")

#rename nitrate/nitrite and ammonium variables
soilConc <- soilConc %>% rename("NOX" ="NO2.NO3.mg.kg.",
                        "NH4"="NH4.mg.kg.")

#Note there are two observations with missing values in the NOX column
#any(is.na(soilConc)) 
#which(is.na(soilConc$NOX))

#Restrict to ambient CO2 and N treatment and to zero-species plots
Amb0 <- soilConc %>% filter(CO2.Treatment == "Camb" & Nitrogen.Treatment == "Namb" & CountOfSpecies == 0)

#Take mean concentrations of nitrate/nitrite (NOX) and ammonium (NH4) by depth, excluding observations with measurement of 0
NOXmean020 <- mean(Amb0$NOX[Amb0$NOX>0 & Amb0$Depth == "0-20"], na.rm = TRUE)
NOXmean4060 <- mean(Amb0$NOX[Amb0$NOX>0 & Amb0$Depth == "40-60"], na.rm = TRUE)
NH4mean020 <- mean(Amb0$NH4[Amb0$NH4>0 & Amb0$Depth == "0-20"], na.rm = TRUE)
NH4mean4060 <- mean(Amb0$NH4[Amb0$NH4>0 & Amb0$Depth == "40-60"], na.rm = TRUE)

#Interpolate and extrapolate soil concentrations linearly by depth in increments of 20cm until negative values of reached
#NOX (nitrate and nitrite)
dnox <- (NOXmean4060-NOXmean020)/2 #the increment dnox is half the difference between measured 0-20cm and 40-60cm values
nox <- NOXmean020 #initialize at top layer
NOXconcList <- c(nox) #stores nitrate concentrations (mg / kg soil) by depth
while (nox > -dnox){
  nox <- nox + dnox
  NOXconcList[(length(NOXconcList)+1)] <- nox
}

#NH4 (ammonium)
dnh4 <- (NH4mean4060-NH4mean020)/2 #the increment da is half the difference between measured 0-20cm and 40-60cm values
nh4 <- NH4mean020 #initialize at top layer
NH4concList <- c(nh4) #stores ammonium concentrations (mg / kg soil) by depth
while (nh4 > -dnh4){
  nh4 <- nh4 + dnh4
  NH4concList[(length(NH4concList)+1)] <- nh4
}

#Add two entries of zero to give NH4concList the same length as NOXconcList
NH4concList[(length(NH4concList)+1)] <- 0
NH4concList[(length(NH4concList)+1)] <- 0

#convert from molecular to elemental nitrogen (N) concentrations
##conversion ratios
NOX_to_N_high <- (14.007)/(14.007+2*15.999) #NO2
NOX_to_N_low <- (14.007)/(14.007+3*15.999) #NO3
NOX_to_N_avg <- 0.5*(NOX_to_N_high + NOX_to_N_low) #average NOX (simple assumption)
NH4_to_N <- (14.007)/(14.007+4*1.0078) #NH4
##elemental nitrogen concentrations (mg/kg)
NconcList <- NOX_to_N_avg * NOXconcList + NH4_to_N * NH4concList

#############################
######### Soil Mass #########
#############################

#import soil bulk density data
soilDens <- read.csv("bde141clean.csv")

#rename density variable (units are g/cm^3)
soilDens <- soilDens %>% rename("Dens" ="Bulk.Density..g.cm.3.")

#Restrict to ambient CO2 and N treatment and to zero-species plots
Amb0_D <- soilDens%>% filter(CO2.Treatment == "Camb" & Nitrogen.Treatment == "Namb" & CountOfSpecies == 0)

#visualize trend in density measurements by depth
#boxplot(Amb0_D$Dens ~ Amb0_D$Depth)

#take mean density across plots and depths (g/cm^3)
meanD <-mean(Amb0_D$Dens)

#calculate volume of one square meter by 20 cm in cm^3
layerVol <- 100 * 100 * 20

#convert soil bulk density to soil mass
layerSoilMass <- meanD * layerVol #units of *(g/cm^3)*(cm^3)=g

#convert soil mass from units of g to kg
layerSoilMass <- 0.001 * layerSoilMass

#############################
####### Nitrogen Mass #######
#############################

#compute N mass in each layer as 
layerNmass <- NconcList * layerSoilMass #units of (mg/kg)*(kg) = mg

#sum across masses to estimate elemental N mass (in grams) below one square meter of surface
Nmass <- sum(layerNmass)

#convert from mg N to g N
Nmass <- 0.001 * Nmass

#Nmass is grams of elemental nitrogen below one square meter of surface.
#Code returns ~0.633 g N per square meter.
print(Nmass)