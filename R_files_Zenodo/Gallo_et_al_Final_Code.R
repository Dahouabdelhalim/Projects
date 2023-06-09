# Code to accompany Gallo et al. Characterizing deep-water oxygen variability and seafloor 
    # community responses using a novel autonomous lander. Biogeosciences.

# Running R version 3.6.0 (2019-04-26) "Planting of a Tree" on platform x86_64-apple-darwin15.6.0
# Version numbers for packages are indicated when packages are loaded

#_____________________________________________________________________________

#clear workspace
rm(list = ls())

#Download all datafiles and set to your working directory
setwd("")

#### Nanolander oceanographic data and description of deployments ####

# Open datasheets
# The datasheet Gallo_et_al_BEEBE_Final.csv contains all of the 7 Nanolander deployments

All = read.csv("Gallo_et_al_BEEBE_Final.csv", header=TRUE) 

# Data were downloaded from the SBE MicroCAT, then data were inspected to remove samples
  #from before the deployment and after the recovery (data now only include seafloor measurements)
# Data were also inspected for outliers. Outliers were removed using 3*IQR and in rare cases, 
    #additional outliers were removed visually
#Information on deployment notation and dates of deployment
#D3 = D200-LJ-1 #Deployed Aug17, Recovered Sep1, Location: Scripps Reserve
#D4 = D200-LJ-2 #Deployed Sep7, Recovered Sep25, Location: Scripps Reserve
#D5 = D100-DM-Fall #Deployed Sep29, Recovered Nov3, Location: Del Mar Steeples
#D6 = D200-DM #Deployed Nov9, Recovered Nov29, Location: Del Mar Steeples
#D7 = D300-DM #Deployed Dec12, Recovered Jan5, Location: Del Mar Steeples
#D8 = D400-DM #Deployed Jan22, Recovered Feb8, Location: Del Mar Steeples
#D10 = D100-DM-Spr #Deployed Mar8, Recovered Apr27, Location: Del Mar Steeples, then Baja

#Notes on outliers removed
  #D3 15 NAs (outliers removed) for Ox and Ox_Perc_Sat
  #D4 no outliers removed
  #D5 no outliers removed
  #D6 no outliers removed
  #D7 no outliers removed
  #D8 2 NAs (outliers removed) for density and salinity, first oxygen and oxygen 
    #saturation measurement replaced with second measurement, due to obvious
    #oxygen sensor lag
  #D10 1 NA (outliers removed) for Ox and Ox_Perc_Sat

#"Spice" calculated using "oce" R package
  #library("oce")
  #All$Spice = swSpice(All$Salinity, temperature = All$Temp, pressure = All$Depth)
#spiciness unit kg/m^3

#"pO2_matm" calculating using Hofmann et al. (2011) Hypoxia by degrees:
  #Establishing definitions for a changing ocean. Supplementary R code. 
    #All$pO2_matm = pO2(All$Ox_umol, All$Salinity, All$Temp, All$Depth, lat = All$Latitude)
    #All$pO2_matm = as.numeric(All$pO2_matm)

#"pO2_kPa" calculated by converting matm to kPa
  #All$pO2_kPa = ((All$pO2_matm)*(1/1000)*(101.325))

#"pHest" Estimate pH calculated using Alin et al. (2012) Robust empirical relationships for estimating the 
  #carbonate system in the southern California Current System and application to CalCOFI 
  #hydrographic cruise data (2005–2011)
      #a0 = 7.758
      #a1 = 1.42*10^(-2)
      #a2 = 1.62*10^(-3)
      #a3 = 4.24*10^(-5)
      #Tr = 10.28
      #O2r = 138.46
  #All$pHest = a0 + (a1*(All$Temp - Tr)) + (a2*(All$Ox_umol - O2r)) + 
      #(a3*((All$Temp - Tr)*(All$Ox_umol - O2r)))

#"AragSSest" Aragonite saturation state calculated based on Alin et al. (2012)
  #Coefficients for aragonite saturation state:
      #a0 = 1.112
      #a1 = 9.59*10^(-2)
      #a2 = 3.54*10^(-3)
      #a3 = 5.91*10^(-4)
      #Tr = 10.28
      #O2r = 138.46
  #All$AragSSest = a0 + (a1*(All$Temp - Tr)) + (a2*(All$Ox_umol - O2r)) + 
      #(a3*((All$Temp - Tr)*(All$Ox_umol - O2r)))

#"CalcSSest" Calcite saturation state calculated based on Alin et al. (2012)
  ##Coefficients for calcite saturation state:
      #a0 = 1.749
      #a1 = 0.147
      #a2 = 5.61*10^(-3)
      #a3 = 8.02*10^(-4)
      #Tr = 10.28
      #O2r = 138.46
  #All$CalcSSest = a0 + (a1*(All$Temp - Tr)) + (a2*(All$Ox_umol - O2r)) + 
      #(a3*((All$Temp - Tr)*(All$Ox_umol - O2r)))

head(All)
str(All)
summary(All)

All$Date = as.character(All$Date) #change from factor to character
All$Time = as.character(All$Time) #change from factor to character
All$DateTime <- paste(All$Date, All$Time, sep=" ") #Makes a column that contains the two
All$Date = as.Date(All$Date, format = "%m/%d/%y") #Make R recognize it as a date class
All$DateTime <- as.POSIXct(All$DateTime, format="%m/%d/%y %H:%M:%S") #format time

str(All)

require("ggplot2")      #packageVersion("ggplot2"): ‘3.3.0’
require("dplyr")        #packageVersion("dplyr"): ‘0.8.1’
require("tidyverse")    #packageVersion("tidyverse"): ‘1.2.1’

#Split into individual deployments (see above for designation); CD simply indicates that it is the
  #clean data
#Calculate %time hypoxic, %severely hypoxic, and %time undersaturated with respect to 
  #aragonite or calcite

D3_CD = subset(All, Deployment == "D3")
Hypoxic_D3 = subset(D3_CD, Ox <= 60) 
UndersatArg_D3 = subset(D3_CD, AragSSest < 1) 
UndersatCalc_D3 = subset(D3_CD, CalcSSest < 1)
Perc_Hypoxic = (544/4303)*100 #D3
Perc_Undersaturated_Arag = (4288/4303)*100 #D3
Perc_Undersaturated_Calc = (0/4303)*100 #D3
#D3 hypoxic conditions present 12.64% of the time
#D3 undersaturated with respect to aragonite 99.65% of the time, but
#never undersaturated with respect to calcite

D4_CD = subset(All, Deployment == "D4")
Hypoxic_D4 = subset(D4_CD, Ox <= 60) 
UndersatArg_D4 = subset(D4_CD, AragSSest < 1) 
UndersatCalc_D4 = subset(D4_CD, CalcSSest < 1)
Perc_Hypoxic = (98/5212)*100 #D4
Perc_Undersaturated_Arag = (5200/5212)*100 #D4
Perc_Undersaturated_Calc = (0/5212)*100 #D4
#D4 hypoxic conditions present 1.88% of the time
#D4 undersaturated with respect to aragonite 99.77% of the time, but
#never undersaturated with respect to calcite

D5_CD = subset(All, Deployment == "D5")
Hypoxic_D5 = subset(D5_CD, Ox <= 60) 
UndersatArg_D5 = subset(D5_CD, AragSSest < 1) 
UndersatCalc_D5 = subset(D5_CD, CalcSSest < 1)
Perc_Undersaturated_Arag = (0/10061)*100 #D5
Perc_Undersaturated_Calc = (0/10061)*100 #D5
#D5 conditions never hypoxic
#D5 never undersaturated with respect to aragonite or calcite

D6_CD = subset(All, Deployment == "D6")
Hypoxic_D6 = subset(D6_CD, Ox <= 60) 
UndersatArg_D6 = subset(D6_CD, AragSSest < 1) 
UndersatCalc_D6 = subset(D6_CD, CalcSSest < 1)
Perc_Undersaturated_Arag = (5757/5757)*100 #D6
Perc_Undersaturated_Calc = (0/5757)*100 #D6
#D6 conditions never hypoxic
#D6 undersaturated with respect to aragonite 100% of the time, but 
#never undersaturated with respect to calcite

D7_CD = subset(All, Deployment == "D7")
Hypoxic_D7 = subset(D7_CD, Ox <= 60) 
OMZ_D7 = subset(D7_CD, Ox <=22.5)
UndersatArg_D7 = subset(D7_CD, AragSSest < 1) 
UndersatCalc_D7 = subset(D7_CD, CalcSSest < 1)
Perc_Hypoxic = (6912/6912)*100 #D7
Perc_Undersaturated_Arag = (6912/6912)*100 #D7
Perc_Undersaturated_Calc = (0/6912)*100 #D7
#D7 hypoxic conditions present 100% of time (never severely hypoxic)
#D7 undersaturated with respect to aragonite 100% of the time, but 
#never undersaturated with respect to calcite

D8_CD = subset(All, Deployment == "D8")
Hypoxic_D8 = subset(D8_CD, Ox <= 60) 
OMZ_D8 = subset(D8_CD, Ox <=22.5)
UndersatArg_D8 = subset(D8_CD, AragSSest < 1) 
UndersatCalc_D8 = subset(D8_CD, CalcSSest < 1)
Perc_Hypoxic = (4663/4663)*100 #D8
Perc_Severely_Hypoxic = (52/4663)*100 #D8
Perc_Undersaturated_Arag = (4663/4663)*100 #D8
Perc_Undersaturated_Calc = (4537/4663)*100 #D8
#D8 hypoxic conditions present 100% of time (severely hypoxic 1.12% of time)
#D8 undersaturated with respect to aragonite 100% of the time, and
#undersaturated with respect to calcite 97.30% of the time

D10_CD = subset(All, Deployment == "D10")
Hypoxic_D10 = subset(D10_CD, Ox <= 60) 
UndersatArg_D10 = subset(D10_CD, AragSSest < 1) 
UndersatCalc_D10 = subset(D10_CD, CalcSSest < 1)
Perc_Undersaturated_Arag = (5682/6123)*100 #D10
Perc_Undersaturated_Calc = (0/6123)*100 #D10
#D10 conditions never hypoxic
#D8 undersaturated with respect to aragonite 92.80% of the time, and
#never undersaturated with respect to calcite

#Determine means and ranges for all variables
summary(D3_CD)
summary(D4_CD)
summary(D5_CD)
summary(D6_CD)
summary(D7_CD)
summary(D8_CD)
summary(D10_CD)

#Calculate the coefficient of dispersion as percentage (sd/mean) for each deployment
CV_Ox_D3 = sd(D3_CD$Ox_umol, na.rm=TRUE)/mean(D3_CD$Ox_umol, na.rm=TRUE)*100
CV_Temp_D3 = sd(D3_CD$Temp, na.rm=TRUE)/mean(D3_CD$Temp, na.rm=TRUE)*100
CV_pH_D3 = sd(D3_CD$pHest, na.rm=TRUE)/mean(D3_CD$pHest, na.rm=TRUE)*100

CV_Ox_D4 = sd(D4_CD$Ox_umol, na.rm=TRUE)/mean(D4_CD$Ox_umol, na.rm=TRUE)*100
CV_Temp_D4 = sd(D4_CD$Temp, na.rm=TRUE)/mean(D4_CD$Temp, na.rm=TRUE)*100
CV_pH_D4 = sd(D4_CD$pHest, na.rm=TRUE)/mean(D4_CD$pHest, na.rm=TRUE)*100

CV_Ox_D5 = sd(D5_CD$Ox_umol, na.rm=TRUE)/mean(D5_CD$Ox_umol, na.rm=TRUE)*100
CV_Temp_D5 = sd(D5_CD$Temp, na.rm=TRUE)/mean(D5_CD$Temp, na.rm=TRUE)*100
CV_pH_D5 = sd(D5_CD$pHest, na.rm=TRUE)/mean(D5_CD$pHest, na.rm=TRUE)*100

CV_Ox_D6 = sd(D6_CD$Ox_umol, na.rm=TRUE)/mean(D6_CD$Ox_umol, na.rm=TRUE)*100
CV_Temp_D6 = sd(D6_CD$Temp, na.rm=TRUE)/mean(D6_CD$Temp, na.rm=TRUE)*100
CV_pH_D6 = sd(D6_CD$pHest, na.rm=TRUE)/mean(D6_CD$pHest, na.rm=TRUE)*100

CV_Ox_D7 = sd(D7_CD$Ox_umol, na.rm=TRUE)/mean(D7_CD$Ox_umol, na.rm=TRUE)*100
CV_Temp_D7 = sd(D7_CD$Temp, na.rm=TRUE)/mean(D7_CD$Temp, na.rm=TRUE)*100
CV_pH_D7 = sd(D7_CD$pHest, na.rm=TRUE)/mean(D7_CD$pHest, na.rm=TRUE)*100

CV_Ox_D8 = sd(D8_CD$Ox_umol, na.rm=TRUE)/mean(D8_CD$Ox_umol, na.rm=TRUE)*100
CV_Temp_D8 = sd(D8_CD$Temp, na.rm=TRUE)/mean(D8_CD$Temp, na.rm=TRUE)*100
CV_pH_D8 = sd(D8_CD$pHest, na.rm=TRUE)/mean(D8_CD$pHest, na.rm=TRUE)*100

CV_Ox_D10 = sd(D10_CD$Ox_umol, na.rm=TRUE)/mean(D10_CD$Ox_umol, na.rm=TRUE)*100
CV_Temp_D10 = sd(D10_CD$Temp, na.rm=TRUE)/mean(D10_CD$Temp, na.rm=TRUE)*100
CV_pH_D10 = sd(D10_CD$pHest, na.rm=TRUE)/mean(D10_CD$pHest, na.rm=TRUE)*100

#Anova to test difference between 100 m Del Mar deployments during fall and spring
Seasonal_Comp = subset(All, Deployment == "D10" | Deployment == "D5")
Seasonal_aov <- aov(Ox_umol ~ Deployment, data=Seasonal_Comp)
summary(Seasonal_aov)
#Df  Sum Sq Mean Sq F value Pr(>F)    
#Deployment      1 2989663 2989663   72477 <2e-16 ***
#  Residuals   16182  667507      41                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### Reproduce Figure 4 ####

Means = All %>% 
  group_by(Deployment) %>%
  summarise_all(list(mean), na.rm = TRUE)

print(levels(All$Deployment))

#Look at probability density of observations for multiple environmental variables for each
  #deployment
Ox_umol.plot = ggplot(All, aes(x=Ox_umol, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=Ox_umol, colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Oxygen concentration (umol/kg)") + 
  theme(legend.position = "none") 
Ox_umol.plot

Ox_sat.plot = ggplot(All, aes(x=Ox_Perc_Sat, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=Ox_Perc_Sat, colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Oxygen saturation (%)") +
  theme(legend.position = "none") 
Ox_sat.plot

pO2_matm.plot = ggplot(All, aes(x=pO2_matm, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=pO2_matm,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("pO2 (matm)") +
  theme(legend.position = "none") 
pO2_matm.plot

pO2_kPa.plot = ggplot(All, aes(x=pO2_kPa, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=pO2_kPa,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("pO2 (kPa)") +
  theme(legend.position = "none") 
pO2_kPa.plot

pHest.plot = ggplot(All, aes(x=pHest, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=pHest,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("pHest") +
  theme(legend.position = "none") 
pHest.plot

AragSSest.plot = ggplot(All, aes(x=AragSSest, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=AragSSest,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Aragonite Saturation State (estimated)") +
  theme(legend.position = "none") 
AragSSest.plot

CalcSSest.plot = ggplot(All, aes(x=CalcSSest, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=CalcSSest,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Calcite Saturation State (estimated)") +
  theme(legend.position = "none") 
CalcSSest.plot

Temp.plot = ggplot(All, aes(x=Temp, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=Temp,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Temperature (°C)") +
  theme(legend.position = "none") 
Temp.plot

Sal.plot = ggplot(All, aes(x=Salinity, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=Salinity,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Salinity (psu)") +
  theme(legend.position = "none") 
Sal.plot

Density.plot = ggplot(All, aes(x=Density, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=Density,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Density (kg/m^3)") +
  theme(legend.position = "none") 
Density.plot

Spice.plot = ggplot(All, aes(x=Spice, fill=Deployment)) +
  geom_density(alpha=.3) +
  geom_vline(data=Means, aes(xintercept=Spice,  colour=Deployment),
             linetype="dashed", size=1) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Spiciness (kg/m^3)") +
  theme(legend.position = "none") 
Spice.plot

#These probability density panels were used in Figure 4
#Figure 4: Mean and variance of near-seafloor temperature, oxygen concentration, 
  #oxygen partial pressure, pHest, salinity, and spiciness. The probability density of data 
  #collected for each deployment is shown, with the color of the data distributions
  #corresponding to each deployment (as indicated in the color legend). The mean is 
  #indicated with a dotted line in the same color and exact values are given in Table 1. 
  #pHest is estimated pH, calculated using empirical relationships from Alin et al. (2012). 
  #Sampling dates for each deployment are given in Table 1.
#Read in Multiplot_Function.R
multiplot(Temp.plot, Ox_umol.plot, pO2_kPa.plot, pHest.plot, Sal.plot, Spice.plot, cols = 1)

#### Reproduce Figure 5 ####

## Figure 5: 2017-2018 near-bottom dissolved oxygen concentration in the Southern 
  #California Bight shown in relation to temperature (A) and spiciness (B). Data points 
  #represent samples taken every five minutes with the SBE MicroCAT-ODO sensor
  #during the seven deployments. Deployments are distinguished by color, as indicated in the 
  #color legend. Sampling dates for each deployment are given in Table 1.
Temp_v_Ox = ggplot(All, aes(x=Temp, y = Ox_umol, col=Deployment)) +
  geom_point() +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Temperature (°C)") +
  ylab("Oxygen concentration (umol/kg)")
Temp_v_Ox

#Test correlation between temperature and O2 with linear regressions
D3_lm_Temp = lm(Ox_umol~Temp,data=subset(All,Deployment=="D3"))
summary(D3_lm_Temp) 
#Multiple R-squared:  0.8033,  Adjusted R-squared:  0.8033 
#F-statistic: 1.751e+04 on 1 and 4287 DF,  p-value: < 2.2e-16
D4_lm_Temp = lm(Ox_umol~Temp,data=subset(All,Deployment=="D4"))
summary(D4_lm_Temp)
#Multiple R-squared:  0.7753,  Adjusted R-squared:  0.7753 
#F-statistic: 1.798e+04 on 1 and 5210 DF,  p-value: < 2.2e-16
D5_lm_Temp = lm(Ox_umol~Temp,data=subset(All,Deployment=="D5"))
summary(D5_lm_Temp)
#Multiple R-squared:  0.8443,  Adjusted R-squared:  0.8443 
#F-statistic: 5.454e+04 on 1 and 10059 DF,  p-value: < 2.2e-16
D6_lm_Temp = lm(Ox_umol~Temp,data=subset(All,Deployment=="D6"))
summary(D6_lm_Temp)
#Multiple R-squared:  0.4063,  Adjusted R-squared:  0.4062 
#F-statistic:  3939 on 1 and 5755 DF,  p-value: < 2.2e-16
D7_lm_Temp = lm(Ox_umol~Temp,data=subset(All,Deployment=="D7"))
summary(D7_lm_Temp)
#Multiple R-squared:  0.8512,  Adjusted R-squared:  0.8511 
#F-statistic: 3.952e+04 on 1 and 6910 DF,  p-value: < 2.2e-16
D8_lm_Temp = lm(Ox_umol~Temp,data=subset(All,Deployment=="D8"))
summary(D8_lm_Temp)
#Multiple R-squared:  0.9027,  Adjusted R-squared:  0.9027 
#F-statistic: 4.323e+04 on 1 and 4661 DF,  p-value: < 2.2e-16
D10_lm_Temp = lm(Ox_umol~Temp,data=subset(All,Deployment=="D10"))
summary(D10_lm_Temp)
#Multiple R-squared:  0.8076,  Adjusted R-squared:  0.8076 
#F-statistic: 2.569e+04 on 1 and 6121 DF,  p-value: < 2.2e-16

Spice_v_Ox = ggplot(All, aes(x=Spice, y = Ox_umol, col=Deployment)) +
  geom_point() +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Spiciness (kg/m^3)") +
  ylab("Oxygen concentration (umol/kg)")
Spice_v_Ox

#Test correlation between spiciness and O2 with linear regressions
D3_lm_Spice = lm(Ox_umol~Spice,data=subset(All,Deployment=="D3"))
summary(D3_lm_Spice)
#Multiple R-squared:  0.9764,  Adjusted R-squared:  0.9764 
#F-statistic: 1.772e+05 on 1 and 4287 DF,  p-value: < 2.2e-16
D4_lm_Spice = lm(Ox_umol~Spice,data=subset(All,Deployment=="D4"))
summary(D4_lm_Spice)
#Multiple R-squared:  0.9191,  Adjusted R-squared:  0.9191 
#F-statistic: 5.921e+04 on 1 and 5210 DF,  p-value: < 2.2e-16
D5_lm_Spice = lm(Ox_umol~Spice,data=subset(All,Deployment=="D5"))
summary(D5_lm_Spice)
#Multiple R-squared:  0.3109,  Adjusted R-squared:  0.3108 
#F-statistic:  4538 on 1 and 10059 DF,  p-value: < 2.2e-16
D6_lm_Spice = lm(Ox_umol~Spice,data=subset(All,Deployment=="D6"))
summary(D6_lm_Spice)
#Multiple R-squared:  0.6142,  Adjusted R-squared:  0.6142 
#F-statistic:  9163 on 1 and 5755 DF,  p-value: < 2.2e-16
D7_lm_Spice = lm(Ox_umol~Spice,data=subset(All,Deployment=="D7"))
summary(D7_lm_Spice)
#Multiple R-squared:  0.6111,  Adjusted R-squared:  0.611 
#F-statistic: 1.086e+04 on 1 and 6910 DF,  p-value: < 2.2e-16
D8_lm_Spice = lm(Ox_umol~Spice,data=subset(All,Deployment=="D8"))
summary(D8_lm_Spice)
#Multiple R-squared:  0.6768,  Adjusted R-squared:  0.6767 
#F-statistic:  9756 on 1 and 4659 DF,  p-value: < 2.2e-16
D10_lm_Spice = lm(Ox_umol~Spice,data=subset(All,Deployment=="D10"))
summary(D10_lm_Spice)
#Multiple R-squared:  0.8139,  Adjusted R-squared:  0.8139 
#F-statistic: 2.677e+04 on 1 and 6121 DF,  p-value: < 2.2e-16

multiplot(Temp_v_Ox, Spice_v_Ox, cols = 1)
#saved as a 7 by 8 pdf titled Ox_Corr.pdf

#### Reproduce Supplement 1E ####

#Plot the time series of oxygen for each deployment for Supplement 1E
Ox_time_D3 = ggplot(D3_CD, aes(x=Days, y = Ox_umol)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  xlab("Days") +
  ylab("Oxygen concentration (umol/kg)") +
  theme(legend.position = "none") 
Ox_time_D3

Ox_time_D4 = ggplot(D4_CD, aes(x=Days, y = Ox_umol)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  xlab("Days") +
  ylab("Oxygen concentration (umol/kg)") +
  theme(legend.position = "none") 
Ox_time_D4

Ox_time_D5 = ggplot(D5_CD, aes(x=Days, y = Ox_umol)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  xlab("Days") +
  ylab("Oxygen concentration (umol/kg)") +
  theme(legend.position = "none") 
Ox_time_D5

Ox_time_D6 = ggplot(D6_CD, aes(x=Days, y = Ox_umol)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  xlab("Days") +
  ylab("Oxygen concentration (umol/kg)") +
  theme(legend.position = "none") 
Ox_time_D6

Ox_time_D7 = ggplot(D7_CD, aes(x=Days, y = Ox_umol)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  xlab("Days") +
  ylab("Oxygen concentration (umol/kg)") +
  theme(legend.position = "none") 
Ox_time_D7

Ox_time_D8 = ggplot(D8_CD, aes(x=Days, y = Ox_umol)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  xlab("Days") +
  ylab("Oxygen concentration (umol/kg)") +
  theme(legend.position = "none") 
Ox_time_D8

Ox_time_D10 = ggplot(D10_CD, aes(x=Days, y = Ox_umol)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  xlab("Days") +
  ylab("Oxygen concentration (umol/kg)") +
  theme(legend.position = "none") 
Ox_time_D10

#Oxygen time series were plotted together and included as Supplement 1E
multiplot(Ox_time_D5, Ox_time_D10, Ox_time_D3, Ox_time_D4, Ox_time_D6,
          Ox_time_D7, Ox_time_D8, cols = 2)

#### Reproduce Supplement 1F ####

#Plot oxygen concentration relative to tidal variability for the deployment time series
#These figures were included as Supplement 1F
Ox_Depth_D3 = ggplot(D3_CD, aes(x=Hours, y = Depth, col = Ox_umol)) +
  geom_line() + theme_classic() + xlim(c(250,350)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Hours") +
  ylab("Depth (m)") +
  ggtitle("D3 Oxygen and the Tide (250-350 hours)") +
  scale_color_gradient(low = "red", high = "blue")
Ox_Depth_D3

Ox_Depth_D4 = ggplot(D4_CD, aes(x=Hours, y = Depth, col = Ox_umol)) +
  geom_line() + theme_classic() + xlim(c(0,150)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Hours") +
  ylab("Depth (m)") +
  ggtitle("D4 Oxygen and the Tide (0-150 hours)") +
  scale_color_gradient(low = "red", high = "blue")
Ox_Depth_D4

Ox_Depth_D5 = ggplot(D5_CD, aes(x=Hours, y = Depth, col = Ox_umol)) +
  geom_line() + theme_classic() + xlim(c(700,800)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Hours") +
  ylab("Depth (m)") +
  ggtitle("D5 Oxygen and the Tide (700-800 hours)") +
  scale_color_gradient(low = "red", high = "blue")
Ox_Depth_D5

Ox_Depth_D6 = ggplot(D6_CD, aes(x=Hours, y = Depth, col = Ox_umol)) +
  geom_line() + theme_classic() + xlim(c(0,150)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Hours") +
  ylab("Depth (m)") +
  ggtitle("D6 Oxygen and the Tide (0-150 hours)") +
  scale_color_gradient(low = "red", high = "blue")
Ox_Depth_D6

Ox_Depth_D7 = ggplot(D7_CD, aes(x=Hours, y = Depth, col = Ox_umol)) +
  geom_line() + theme_classic() + xlim(c(400,600)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Hours") +
  ylab("Depth (m)") +
  ggtitle("D7 Oxygen and the Tide (400-600 hours)") +
  scale_color_gradient(low = "red", high = "blue")
Ox_Depth_D7

Ox_Depth_D8 = ggplot(D8_CD, aes(x=Hours, y = Depth, col = Ox_umol)) +
  geom_line() + theme_classic() + xlim(c(150,300)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Hours") +
  ylab("Depth (m)") +
  ggtitle("D8 Oxygen and the Tide (150-300 hours)") +
  scale_color_gradient(low = "red", high = "blue")
Ox_Depth_D8

Ox_Depth_D10 = ggplot(D10_CD, aes(x=Hours, y = Depth, col = Ox_umol)) +
  geom_line() + theme_classic() + xlim(c(0,150)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=13)) +
  theme(axis.text.y = element_text(size=13)) +
  xlab("Hours") +
  ylab("Depth (m)") +
  ggtitle("D10 Oxygen and the Tide (0-150 hours)") +
  scale_color_gradient(low = "red", high = "blue")
Ox_Depth_D10

## Also plot oxygen time-series relative to tidal cycle
Tidal_D3 = ggplot(D3_CD, aes(x=Days, y = Depth)) +
  geom_line() + theme_classic() + theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) + theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) + xlab("Days") + ylab("Depth (m)") +
  theme(legend.position = "none") 
Tidal_D3

Tidal_D4 = ggplot(D4_CD, aes(x=Days, y = Depth)) +
  geom_line() + theme_classic() + theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) + theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) + xlab("Days") + ylab("Depth (m)") +
  theme(legend.position = "none") 
Tidal_D4

Tidal_D5 = ggplot(D5_CD, aes(x=Days, y = Depth)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) + theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) + theme(axis.text.y = element_text(size=12)) +
  xlab("Days") + ylab("Depth (m)") + theme(legend.position = "none") 
Tidal_D5

Tidal_D6 = ggplot(D6_CD, aes(x=Days, y = Depth)) +
  geom_line() + theme_classic() +
  theme(axis.title.x = element_text(size=12)) + theme(axis.title.y = element_text(size=12)) +
  theme(axis.text.x = element_text(size=12)) + theme(axis.text.y = element_text(size=12)) +
  xlab("Days") + ylab("Depth (m)") + theme(legend.position = "none") 
Tidal_D6

Tidal_D7 = ggplot(D7_CD, aes(x=Days, y = Depth)) +
  geom_line() + theme_classic() + theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) + theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) + xlab("Days") + ylab("Depth (m)") +
  theme(legend.position = "none") 
Tidal_D7

Tidal_D8 = ggplot(D8_CD, aes(x=Days, y = Depth)) +
  geom_line() + theme_classic() + theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) + theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) + xlab("Days") + ylab("Depth (m)") +
  theme(legend.position = "none") 
Tidal_D8

Tidal_D10 = ggplot(D10_CD, aes(x=Days, y = Depth)) +
  geom_line() + theme_classic() + theme(axis.title.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=12)) + theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) + xlab("Days") + ylab("Depth (m)") +
  theme(legend.position = "none") 
Tidal_D10

multiplot(Ox_time_D5, Tidal_D5, Ox_time_D10, Tidal_D10, cols = 2)
multiplot(Ox_time_D3, Tidal_D3, Ox_time_D4, Tidal_D4, cols = 2)
multiplot(Ox_time_D6, Tidal_D6, Ox_time_D7, Tidal_D7, cols = 2)
multiplot(Ox_time_D8, Tidal_D8, cols = 2)

## Now look for any linear relationships between amount of tidal variability
  # and amount of oxygen variability in one day; compare across deployments
CV_D3 = D3_CD %>%
  group_by(Date) %>%
  summarize(meanOx = mean(Ox_umol, na.rm = TRUE),
            meanDepth = mean(Depth, na.rm = TRUE),
            stdevOx = sd(Ox_umol, na.rm = TRUE),
            stdevDepth = sd(Depth, na.rm = TRUE),
            minOx = min(Ox_umol, na.rm = TRUE),
            maxOx = max(Ox_umol, na.rm = TRUE),
            minDepth = min(Depth, na.rm = TRUE),
            maxDepth = max(Depth, na.rm = TRUE),
            n=n()) %>%
  mutate(CV_Ox = (stdevOx/meanOx)*100, 
         CV_Depth = (stdevDepth/meanDepth)*100,
         Range_Ox = maxOx-minOx,
         Range_Depth = maxDepth-minDepth,
         Deployment = "D3",
         Depth_CAT = "200 m") %>%
  as.data.frame()

CV_D4 = D4_CD %>%
  group_by(Date) %>%
  summarize(meanOx = mean(Ox_umol, na.rm = TRUE),
            meanDepth = mean(Depth, na.rm = TRUE),
            stdevOx = sd(Ox_umol, na.rm = TRUE),
            stdevDepth = sd(Depth, na.rm = TRUE),
            minOx = min(Ox_umol, na.rm = TRUE),
            maxOx = max(Ox_umol, na.rm = TRUE),
            minDepth = min(Depth, na.rm = TRUE),
            maxDepth = max(Depth, na.rm = TRUE),
            n=n()) %>%
  mutate(CV_Ox = (stdevOx/meanOx)*100, 
         CV_Depth = (stdevDepth/meanDepth)*100,
         Range_Ox = maxOx-minOx,
         Range_Depth = maxDepth-minDepth,
         Deployment = "D4",
         Depth_CAT = "200 m") %>%
  as.data.frame()

CV_D5 = D5_CD %>%
  group_by(Date) %>%
  summarize(meanOx = mean(Ox_umol, na.rm = TRUE),
            meanDepth = mean(Depth, na.rm = TRUE),
            stdevOx = sd(Ox_umol, na.rm = TRUE),
            stdevDepth = sd(Depth, na.rm = TRUE),
            minOx = min(Ox_umol, na.rm = TRUE),
            maxOx = max(Ox_umol, na.rm = TRUE),
            minDepth = min(Depth, na.rm = TRUE),
            maxDepth = max(Depth, na.rm = TRUE),
            n=n()) %>%
  mutate(CV_Ox = (stdevOx/meanOx)*100, 
         CV_Depth = (stdevDepth/meanDepth)*100,
         Range_Ox = maxOx-minOx,
         Range_Depth = maxDepth-minDepth,
         Deployment = "D5",
         Depth_CAT = "100 m") %>%
  as.data.frame()

CV_D6 = D6_CD %>%
  group_by(Date) %>%
  summarize(meanOx = mean(Ox_umol, na.rm = TRUE),
            meanDepth = mean(Depth, na.rm = TRUE),
            stdevOx = sd(Ox_umol, na.rm = TRUE),
            stdevDepth = sd(Depth, na.rm = TRUE),
            minOx = min(Ox_umol, na.rm = TRUE),
            maxOx = max(Ox_umol, na.rm = TRUE),
            minDepth = min(Depth, na.rm = TRUE),
            maxDepth = max(Depth, na.rm = TRUE),
            n=n()) %>%
  mutate(CV_Ox = (stdevOx/meanOx)*100, 
         CV_Depth = (stdevDepth/meanDepth)*100,
         Range_Ox = maxOx-minOx,
         Range_Depth = maxDepth-minDepth,
         Deployment = "D6",
         Depth_CAT = "200 m") %>%
  as.data.frame()

CV_D7 = D7_CD %>%
  group_by(Date) %>%
  summarize(meanOx = mean(Ox_umol, na.rm = TRUE),
            meanDepth = mean(Depth, na.rm = TRUE),
            stdevOx = sd(Ox_umol, na.rm = TRUE),
            stdevDepth = sd(Depth, na.rm = TRUE),
            minOx = min(Ox_umol, na.rm = TRUE),
            maxOx = max(Ox_umol, na.rm = TRUE),
            minDepth = min(Depth, na.rm = TRUE),
            maxDepth = max(Depth, na.rm = TRUE),
            n=n()) %>%
  mutate(CV_Ox = (stdevOx/meanOx)*100, 
         CV_Depth = (stdevDepth/meanDepth)*100,
         Range_Ox = maxOx-minOx,
         Range_Depth = maxDepth-minDepth,
         Deployment = "D7",
         Depth_CAT = "300 m") %>%
  as.data.frame()

CV_D8 = D8_CD %>%
  group_by(Date) %>%
  summarize(meanOx = mean(Ox_umol, na.rm = TRUE),
            meanDepth = mean(Depth, na.rm = TRUE),
            stdevOx = sd(Ox_umol, na.rm = TRUE),
            stdevDepth = sd(Depth, na.rm = TRUE),
            minOx = min(Ox_umol, na.rm = TRUE),
            maxOx = max(Ox_umol, na.rm = TRUE),
            minDepth = min(Depth, na.rm = TRUE),
            maxDepth = max(Depth, na.rm = TRUE),
            n=n()) %>%
  mutate(CV_Ox = (stdevOx/meanOx)*100, 
         CV_Depth = (stdevDepth/meanDepth)*100,
         Range_Ox = maxOx-minOx,
         Range_Depth = maxDepth-minDepth,
         Deployment = "D8",
         Depth_CAT = "400 m") %>%
  as.data.frame()

CV_D10 = D10_CD %>%
  group_by(Date) %>%
  summarize(meanOx = mean(Ox_umol, na.rm = TRUE),
            meanDepth = mean(Depth, na.rm = TRUE),
            stdevOx = sd(Ox_umol, na.rm = TRUE),
            stdevDepth = sd(Depth, na.rm = TRUE),
            minOx = min(Ox_umol, na.rm = TRUE),
            maxOx = max(Ox_umol, na.rm = TRUE),
            minDepth = min(Depth, na.rm = TRUE),
            maxDepth = max(Depth, na.rm = TRUE),
            n=n()) %>%
  mutate(CV_Ox = (stdevOx/meanOx)*100, 
         CV_Depth = (stdevDepth/meanDepth)*100,
         Range_Ox = maxOx-minOx,
         Range_Depth = maxDepth-minDepth,
         Deployment = "D10",
         Depth_CAT = "100 m") %>%
  as.data.frame()

#100 m deployments
mean(CV_D5$CV_Ox) #3.421773
mean(CV_D10$CV_Ox) #3.246505

#200 m deployments
mean(CV_D3$CV_Ox) #11.90748
mean(CV_D4$CV_Ox) #11.44809
mean(CV_D6$CV_Ox) #6.759599

#300 m deployment
mean(CV_D7$CV_Ox) #5.257815

#400 m deployment
mean(CV_D8$CV_Ox) #6.617346

#100 m deployments
mean(CV_D5$Range_Ox) #20.16778
mean(CV_D10$Range_Ox) #14.19

#200 m deployments
mean(CV_D3$Range_Ox) #31.35813
mean(CV_D4$Range_Ox) #34.17632
mean(CV_D6$Range_Ox) #23.86429

#300 m deployment
mean(CV_D7$Range_Ox) #11.2112

#400 m deployment
mean(CV_D8$Range_Ox) #8.256471

Full_Updated <- rbind(CV_D3, CV_D4, CV_D5, CV_D6, CV_D7, CV_D8, CV_D10)
head(Full_Updated)
str(Full_Updated)
Full_Updated$Deployment = as.factor(Full_Updated$Deployment)
Full_Updated$Depth_CAT = as.factor(Full_Updated$Depth_CAT)

#Remove any days where there is not a full day of samples (< 288 points)
Full_Updated = Full_Updated %>%
  filter(n == 288) %>%
  as.data.table()

ggplot(data = Full_Updated, aes(x = Range_Depth, y = Range_Ox, col = Deployment)) +
  geom_point() + geom_smooth(method = lm) + xlab("Daily Range in Depth (m)") +
  ylab("Daily Range in Oxygen Concentration (umol/kg)") +
  theme(legend.position = "none") 

ggplot(data = Full_Updated, aes(x = CV_Depth, y = CV_Ox, col = Deployment)) +
  geom_point() + geom_smooth(method = lm) #+ facet_grid(Deployment ~ .)

head(Full_Updated)

Full_Updated_100m = Full_Updated %>%
filter(Depth_CAT == "100 m") %>%
  as.data.table()
head(Full_Updated_100m)
min(Full_Updated_100m$Range_Ox) #6.86 umol/kg
max(Full_Updated_100m$Range_Ox) #34.04 umol/kg

Full_Updated_200m = Full_Updated %>%
  filter(Depth_CAT == "200 m") %>%
  as.data.table()
head(Full_Updated_200m)
min(Full_Updated_200m$Range_Ox) #14.68 umol/kg
max(Full_Updated_200m$Range_Ox) #46.23 umol/kg

Full_Updated_300m = Full_Updated %>%
  filter(Depth_CAT == "300 m") %>%
  as.data.table()
head(Full_Updated_300m)
min(Full_Updated_300m$Range_Ox) #8.01 umol/kg
max(Full_Updated_300m$Range_Ox) #14.85 umol/kg

Full_Updated_400m = Full_Updated %>%
  filter(Depth_CAT == "400 m") %>%
  as.data.table()
head(Full_Updated_400m)
min(Full_Updated_400m$Range_Ox) #3.45 umol/kg
max(Full_Updated_400m$Range_Ox) #11.71 umol/kg

#### Reproduce Spectral analysis of oxygen concentration for Supplement 1D and 1E ####

#Simple spectral analysis notes from https://stat.ethz.ch/pipermail/r-help/2007-January/123320.html
#https://anomaly.io/seasonal-trend-decomposition-in-r/index.html

library(zoo)    #packageVersion("zoo"): ‘1.8.6’

All$Mins = as.numeric(All$Mins)
ggplot(data = All, aes(x = DateTime, y = Ox_umol, col = Deployment)) +
  geom_point()

## For Deployment 3 (D3)
Ox_D3 <- D3_CD$Ox_umol
Time_D3 <- D3_CD$DateTime
plot(Ox_D3 ~ Time_D3)
Ox_D3 <- na.locf(Ox_D3) #fills NAs using zoo package, replaces NA with most recent, non-NA value prior to it
Ox_D3_ts = ts(Ox_D3, freq = 288) #288 samples a day
Deco_Ox_D3 = decompose(Ox_D3_ts, type ="additive")
plot(Deco_Ox_D3)
N <- length(Time_D3) #4303

#detrend data
trend = stl(Ox_D3_ts, s.window = "periodic")$time.series[,2]
detrend_Ox_D3 = Ox_D3_ts - (trend - trend[1])
transform <- fft(detrend_Ox_D3) # Using fft
dc <- Mod(transform[1])/N # Extract DC component from transform (mean of time-series)
periodogram <- round( Mod(transform)^2/N, 3)
periodogram <- periodogram[-1] # Drop first element, which is the mean
periodogram <- periodogram[1:(N/2)] # keep first half up to Nyquist limit
# Approximate number of data points in single cycle:
print( N / which(max(periodogram) == periodogram) ) 
#148.3793 data-points in single cycle, or 12.36494 hr

# plot spectrum against Fourier Frequency index for D3
plot(periodogram, col="red", type="o",
     xlab="Fourier Frequency Index", xlim=c(0,100),
     ylab="Periodogram",
     main="Periodogram derived from 'fft' for D3")

## For Deployment 4 (D4)
Ox_D4 <- D4_CD$Ox_umol
Time_D4 <- D4_CD$DateTime
plot(Ox_D4 ~ Time_D4)
Ox_D4_ts = ts(Ox_D4, freq = 288) #288 samples a day
Deco_Ox_D4 = decompose(Ox_D4_ts, type ="additive")
plot(Deco_Ox_D4)

N <- length(Time_D4) #5212
#detrend data
trend = stl(Ox_D4_ts, s.window = "periodic")$time.series[,2]
detrend_Ox_D4 = Ox_D4_ts - (trend - trend[1])
transform <- fft(detrend_Ox_D4) # Using fft
dc <- Mod(transform[1])/N # Extract DC component from transform (mean of time-series)
periodogram <- round( Mod(transform)^2/N, 3)
periodogram <- periodogram[-1] # Drop first element, which is the mean
periodogram <- periodogram[1:(N/2)] # keep first half up to Nyquist limit
# Approximate number of data points in single cycle:
print( N / which(max(periodogram) == periodogram) ) 
#148.9143 data-points in single cycle, or 12.40953 hr

# plot spectrum against Fourier Frequency index
plot(periodogram, col="red", type="o",
     xlab="Fourier Frequency Index", xlim=c(0,100),
     ylab="Periodogram",
     main="Periodogram derived from 'fft' for D4")

## For Deployment 5 (D5)
Ox_D5 <- D5_CD$Ox_umol
Time_D5 <- D5_CD$DateTime
plot(Ox_D5 ~ Time_D5)
Ox_D5_ts = ts(Ox_D5, freq = 288) #288 samples a day
Deco_Ox_D5 = decompose(Ox_D5_ts, type ="additive")
plot(Deco_Ox_D5)

N <- length(Time_D5) #10061
#detrend data
trend = stl(Ox_D5_ts, s.window = "periodic")$time.series[,2]
detrend_Ox_D5 = Ox_D5_ts - (trend - trend[1])
transform <- fft(detrend_Ox_D5) # Using fft
dc <- Mod(transform[1])/N # Extract DC component from transform (mean of time-series)
periodogram <- round( Mod(transform)^2/N, 3)
periodogram <- periodogram[-1] # Drop first element, which is the mean
periodogram <- periodogram[1:(N/2)] # keep first half up to Nyquist limit
# Approximate number of data points in single cycle:
print( N / which(max(periodogram) == periodogram) ) 
#143.7286 data-points in single cycle, or 11.97738 hr

# plot spectrum against Fourier Frequency index
plot(periodogram, col="red", type="o",
     xlab="Fourier Frequency Index", xlim=c(0,150),
     ylab="Periodogram",
     main="Periodogram derived from 'fft' for D5")

## For Deployment 6 (D6)
Ox_D6 <- D6_CD$Ox_umol
Time_D6 <- D6_CD$DateTime
plot(Ox_D6 ~ Time_D6)
Ox_D6_ts = ts(Ox_D6, freq = 288) #288 samples a day
Deco_Ox_D6 = decompose(Ox_D6_ts, type ="additive")
plot(Deco_Ox_D6)

N <- length(Time_D6) #5757
#detrend data
trend = stl(Ox_D6_ts, s.window = "periodic")$time.series[,2]
detrend_Ox_D6 = Ox_D6_ts - (trend - trend[1])
transform <- fft(detrend_Ox_D6) # Using fft
dc <- Mod(transform[1])/N # Extract DC component from transform (mean of time-series)
periodogram <- round( Mod(transform)^2/N, 3)
periodogram <- periodogram[-1] # Drop first element, which is the mean
periodogram <- periodogram[1:(N/2)] # keep first half up to Nyquist limit
# Approximate number of data points in single cycle:
print( N / which(max(periodogram) == periodogram) ) 
#147.6154 data-points in single cycle, or 12.30128 hr

# plot spectrum against Fourier Frequency index
plot(periodogram, col="red", type="o",
     xlab="Fourier Frequency Index", xlim=c(0,150),
     ylab="Periodogram",
     main="Periodogram derived from 'fft' for D6")

## For Deployment 7 (D7)
Ox_D7 <- D7_CD$Ox_umol
Time_D7 <- D7_CD$DateTime
plot(Ox_D7 ~ Time_D7)
Ox_D7_ts = ts(Ox_D7, freq = 288) #288 samples a day
Deco_Ox_D7 = decompose(Ox_D7_ts, type ="additive")
plot(Deco_Ox_D7)

N <- length(Time_D7) #6912
#detrend data
trend = stl(Ox_D7_ts, s.window = "periodic")$time.series[,2]
detrend_Ox_D7 = Ox_D7_ts - (trend - trend[1])
transform <- fft(detrend_Ox_D7) # Using fft
dc <- Mod(transform[1])/N # Extract DC component from transform (mean of time-series)
periodogram <- round( Mod(transform)^2/N, 3)
periodogram <- periodogram[-1] # Drop first element, which is the mean
periodogram <- periodogram[1:(N/2)] # keep first half up to Nyquist limit
# Approximate number of data points in single cycle:
print( N / which(max(periodogram) == periodogram) ) 
#160.7442 data-points in single cycle, or 13.39535 hr

# plot spectrum against Fourier Frequency index
plot(periodogram, col="red", type="o",
     xlab="Fourier Frequency Index", xlim=c(0,300),
     ylab="Periodogram",
     main="Periodogram derived from 'fft' for D7")

## For Deployment 8 (D8)
Ox_D8 <- D8_CD$Ox_umol
Time_D8 <- D8_CD$DateTime
plot(Ox_D8 ~ Time_D8)
Ox_D8_ts = ts(Ox_D8, freq = 288) #288 samples a day
Deco_Ox_D8 = decompose(Ox_D8_ts, type ="additive")
plot(Deco_Ox_D8)

N <- length(Time_D8) #4663
#detrend data
trend = stl(Ox_D8_ts, s.window = "periodic")$time.series[,2]
detrend_Ox_D8 = Ox_D8_ts - (trend - trend[1])
transform <- fft(detrend_Ox_D8) # Using fft
dc <- Mod(transform[1])/N # Extract DC component from transform (mean of time-series)
periodogram <- round( Mod(transform)^2/N, 3)
periodogram <- periodogram[-1] # Drop first element, which is the mean
periodogram <- periodogram[1:(N/2)] # keep first half up to Nyquist limit
# Approximate number of data points in single cycle:
print( N / which(max(periodogram) == periodogram) ) 
#145.7188 data-points in single cycle, or 12.14323 hr

# plot spectrum against Fourier Frequency index
plot(periodogram, col="red", type="o",
     xlab="Fourier Frequency Index", xlim=c(0,150),
     ylab="Periodogram",
     main="Periodogram derived from 'fft' for D8")

## For Deployment 10 (D10)
Ox_D10 <- D10_CD$Ox_umol
Time_D10 <- D10_CD$DateTime
plot(Ox_D10 ~ Time_D10)
Ox_D10_ts = ts(Ox_D10, freq = 288) #288 samples a day
Deco_Ox_D10 = decompose(Ox_D10_ts, type ="additive")
plot(Deco_Ox_D10)

N <- length(Time_D10) #6123
#detrend data
trend = stl(Ox_D10_ts, s.window = "periodic")$time.series[,2]
detrend_Ox_D10 = Ox_D10_ts - (trend - trend[1])
transform <- fft(detrend_Ox_D10) # Using fft
dc <- Mod(transform[1])/N # Extract DC component from transform (mean of time-series)
periodogram <- round( Mod(transform)^2/N, 3)
periodogram <- periodogram[-1] # Drop first element, which is the mean
periodogram <- periodogram[1:(N/2)] # keep first half up to Nyquist limit
# Approximate number of data points in single cycle:
print( N / which(max(periodogram) == periodogram) ) 
#149.3415 data-points in single cycle, or 12.44512 hr

# plot spectrum against Fourier Frequency index
plot(periodogram, col="red", type="o",
     xlab="Fourier Frequency Index", xlim=c(0,150),
     ylab="Periodogram",
     main="Periodogram derived from 'fft' for D10")

#### Reproduce Community Analyses for Figure 6 ####

#Load packages
library(tidyverse)    #packageVersion("tidyverse"): ‘1.2.1’
library(vegan)        #packageVersion("vegan"): ‘2.5.5’
library(ggplot2)      #packageVersion("ggplot2"): ‘3.3.0’

# Open datasheets

#The datasheet "Comm_Dominants_Morethan1.csv" contains nanolander community data. Rows represent
  #20-second video samples, columns represent unique species, counts represent number of individuals
  #of each unique species observed during each 20 second video clip. 
#Only samples with good visibility were retained (Visibility_CAT == "A")
#All samples with no animal observations were removed (n = 49 samples)
#All samples with only one animal observation were removed (n = 105 samples)
#Species that had fewer than 8 observations across all samples were removed. Removed species were:
    #Scorpaena_gutatta, Sebastes_caurinus, Sebastes_ovalis, Oxylebius_pictus, Seriphus_politus,
    #Lyconema_barbatum, Anarrhichthys_ocellatus, Unidentified_shark, Squalus_suckleyi, Sebastes_diploproa,
    #Apristurus_brunneus, Lyopsetta_exilis, Raja_rhina, Skate, Cephaloscyllium_ventriosum,
    #Algae_covered_crab, Lopholithodes_foraminatus, Paralithodes_californiensis, Octopus_californicus,
    #Caridean_shrimp

Community_Dominants_Morethan1 = read.csv("Comm_Dominants_Morethan1.csv", header=TRUE, row.names = 1)
head(Community_Dominants_Morethan1)
Community_Env_Matrix = read.csv("Community_Env_Matrix.csv", header=TRUE, row.names = 1)
head(Community_Env_Matrix)

#All blanks automatically got replaced with NAs, now change all NAs to zeroes 
Community_Dominants_Morethan1[is.na(Community_Dominants_Morethan1)] <- 0 
head(Community_Dominants_Morethan1)
str(Community_Dominants_Morethan1)
#3357 samples; 43 species

#Can either run this analysis or (since it's time consuming), just open
  #the .RData file
load("nMDS_output.RData")

NMDS_AllDeployments_Dominants_Morethan1=metaMDS(Community_Dominants_Morethan1,k=2,trymax=10)
NMDS_AllDeployments_Dominants_Morethan1 #no convergent solutions, stress is 0.1089617 
stressplot(NMDS_AllDeployments_Dominants_Morethan1)
plot(NMDS_AllDeployments_Dominants_Morethan1, scaling = 3)

##Now create environmental matrices to go along with the community matrices for the nMDS plots
plot(NMDS_AllDeployments_Dominants_Morethan1, scaling = 3, xlim = c(-1,1))
ordihull(NMDS_AllDeployments_Dominants_Morethan1, 
         groups=Community_Env_Matrix$Deployment,draw="polygon",col="grey90",label=T)
orditorp(NMDS_AllDeployments_Dominants_Morethan1,display="species",col="red",air=0.01) #all species look good
orditorp(NMDS_AllDeployments_Dominants_Morethan1,display="sites",col="black",air=0.01)

All_colvec = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
print(levels(Community_Env_Matrix$Deployment)) 
#"D10" "D4"  "D5"  "D6"  "D7"  "D8" 
#"#999999" = D10, 
#"#E69F00" = D4, 
#"#56B4E9" = D5, 
#"#009E73" = D6, 
#"#F0E442" = D7, 
#"#0072B2" = D8, 
#"#D55E00" = EXTRA, 
#"#CC79A7" = EXTRA)

## Figure 6A
plot(NMDS_AllDeployments_Dominants_Morethan1, type = "n", xlim = c(-1.5,1.5))
with(Community_Env_Matrix, points(NMDS_AllDeployments_Dominants_Morethan1, 
                display = "sites", col = All_colvec[Deployment], pch = 19, 
                bg = All_colvec[Deployment]))
ordihull(NMDS_AllDeployments_Dominants_Morethan1, 
         groups=Community_Env_Matrix$Deployment,draw="polygon",
         col=All_colvec,label=F)

#### Figure 6B ####
plot(NMDS_AllDeployments_Dominants_Morethan1, type = "n", xlim = c(-1.5,1.5))
with(Community_Env_Matrix, points(NMDS_AllDeployments_Dominants_Morethan1, 
    display = "sites", col = All_colvec[FISH_INVERT_CAT], pch = 19, bg = All_colvec[FISH_INVERT_CAT]))
ordisurf(NMDS_AllDeployments_Dominants_Morethan1, 
         Community_Env_Matrix$Ox_umol, main="",col="blue", add = TRUE)

## Now look only at community similarity across the 2 200 m deployments

### Only for D4
## Any species that had 2 or fewer observations were removed 
    #these were: Sebastes_semicinctus, Zaniolepis_spp, Pacific_sanddab
## Euphausiacea removed
## Any samples with < 2 observations removed
## D4_REC_0729, D4_REC_0262, D4_REC_0205, D4_REC_0256 removed

D4_Only = read.csv("D4_Only_Community.csv", header=TRUE, row.names = 1)
head(D4_Only)
#All blanks automatically got replaced with NAs, now change all NAs to zeroes 
D4_Only[is.na(D4_Only)] <- 0 
head(D4_Only)
str(D4_Only)
#844 samples; 19 species

D4_Only_nMDS=metaMDS(D4_Only,k=2,trymax=10)
D4_Only_nMDS #no convergent solutions, stress is 0.1962331
stressplot(D4_Only_nMDS)
plot(D4_Only_nMDS, scaling = 3)
orditorp(D4_Only_nMDS,display="sites",col="black",air=0.01)

D4_Only_Env = read.csv("D4_Only_Env.csv", header=TRUE)
D4_Only_Env2 = read.csv("D4_Only_Env2.csv", header=TRUE)
head(D4_Only_Env)
str(D4_Only_Env2)

#And you need to remerge it so the order is preserved as for the community matrix
D4_Only_Env3 = left_join(D4_Only_Env, D4_Only_Env2, by = "Sample_ID") 
head(D4_Only_Env$Sample_ID)
head(D4_Only_Env3$Sample_ID)
#check to make sure the sample_ID order is the same #YES

#### Figure 6D ####

#Now we can plot the NMDS
# NMDS output with: 
#colors by Day_Night #goldenrod2 = day, darkslateblue = night
#with ellipses (50% SD)
Diurnal_colvec = c("goldenrod2", "darkslateblue")
plot(D4_Only_nMDS, type = "n", xlim = c(-1.5,2))
with(D4_Only_Env3, points(D4_Only_nMDS, display = "sites", col = Diurnal_colvec[Day_Night], 
                          pch = 19, bg = Diurnal_colvec[Day_Night]))
#orditorp(D4_Only_nMDS,display="species",col="red",air=0.01)
ordiellipse(D4_Only_nMDS, groups = D4_Only_Env3$Day_Night, kind = "sd", conf = 0.50,
            draw = "polygon", col=Diurnal_colvec, 
            label = T)

#### Figure 6E ####

# NMDS output with: 
#colors by Ox_CAT_Extremes_umol #cyan3 = High, lightsteelblue4 = Intermediate, red1 = Low
#with ellipses (50% SD)
Ox_colvec = c("cyan3", "lightsteelblue4", "red1")
print(levels(D4_Only_Env3$Ox_CAT_Extremes_umol)) 

plot(D4_Only_nMDS, type = "n", xlim = c(-1.5,2))
with(D4_Only_Env3, points(D4_Only_nMDS, display = "sites", col = Ox_colvec[Ox_CAT_Extremes_umol], 
                          pch = 19, bg = Ox_colvec[Ox_CAT_Extremes_umol]))
ordiellipse(D4_Only_nMDS, groups = D4_Only_Env3$Ox_CAT_Extremes_umol, kind = "sd", conf = 0.25,
            draw = "polygon", col= Ox_colvec, label = F)
#orditorp(D4_Only_nMDS,display="species",col="red",air=0.01)

### Only for D6
## Any species that have 2 or fewer observations were removed
    #These were: Lycodes_sp., Hydrolagus_colliei, Mysids
## Any samples with < 1 observation removed (none)

D6_Only = read.csv("D6_Only.csv", header=TRUE, row.names = 1)
head(D6_Only)
#All blanks automatically got replaced with NAs, now change all NAs to zeroes 
D6_Only[is.na(D6_Only)] <- 0 
head(D6_Only)
str(D6_Only)
#645 samples; 17 species

D6_Only_nMDS=metaMDS(D6_Only,k=2,trymax=10)
D6_Only_nMDS #no convergent solutions, stress is 0.1858216
stressplot(D6_Only_nMDS)
plot(D6_Only_nMDS, scaling = 3)
orditorp(D6_Only_nMDS,display="species",col="red",air=0.01)

D6_Only_Env = read.csv("D6_Only_Env.csv", header=TRUE)
D6_Only_Env2 = read.csv("D6_Only_Env2.csv", header=TRUE)
head(D6_Only_Env2)
str(D6_Only_Env2)

#And you need to remerge it so the order is preserved as for the community matrix
D6_Only_Env3 = left_join(D6_Only_Env, D6_Only_Env2, by = "Sample_ID") 
head(D6_Only_Env$Sample_ID)
head(D6_Only_Env3$Sample_ID)
#check to make sure the sample_ID order is the same #YES

#### Figure 6F ####

#Now we can plot the NMDS
# NMDS output with: 
#colors by Day_Night #goldenrod2 = day, darkslateblue = night
#with ellipses (50% SD)
Diurnal_colvec = c("goldenrod2", "darkslateblue")
plot(D6_Only_nMDS, type = "n", xlim = c(-1,1))
with(D6_Only_Env3, points(D6_Only_nMDS, display = "sites", col = Diurnal_colvec[Day_Night], 
                          pch = 19, bg = Diurnal_colvec[Day_Night]))
#orditorp(D6_Only_nMDS,display="species",col="red",air=0.01)
ordiellipse(D6_Only_nMDS, groups = D6_Only_Env3$Day_Night, kind = "sd", conf = 0.50,
            draw = "polygon", col=Diurnal_colvec, 
            label = T)

#### Figure 6G ####

# NMDS output with: 
#colors by Ox_CAT_Extremes_umol #cyan3 = High, lightsteelblue4 = Intermediate, red1 = Low
#with ellipses (50% SD)
Ox_colvec = c("cyan3", "lightsteelblue4", "red1")
plot(D6_Only_nMDS, type = "n", xlim = c(-1,1))
with(D6_Only_Env3, points(D6_Only_nMDS, display = "sites", col = Ox_colvec[Ox_CAT_Extremes_umol], 
                          pch = 19, bg = Ox_colvec[Ox_CAT_Extremes_umol]))
#orditorp(D6_Only_nMDS,display="species",col="red",air=0.01)
ordiellipse(D6_Only_nMDS, groups = D6_Only_Env3$Ox_CAT_Extremes_umol, kind = "sd", conf = 0.50,
            draw = "polygon", col = Ox_colvec, 
            label = F)

#### Species Accumulation Code for Figure 6C ####

## All counts from the community matrix were converted to 0's or 1's (presence
# v. absence) for this analysis. 

# Load required packages
library(vegan)      #packageVersion("vegan"): ‘2.5.5’
library(ggplot2)    #packageVersion("ggplot2"): ‘3.3.0’
library(scales)     #packageVersion("scales"): ‘1.1.1’

#can load the final workspace in using this or run the full analysis yourself
load("Fish_Species_Accumulation_Curves.RData")

#Load data file
Fish_Comm <- read.csv("FishOnlyNonZeroOnly.csv")
Fish_Comm$X = NULL

# Turn all species amounts into presence / absence (0 and 1). 
Fish_Comm[, 5:(ncol(Fish_Comm))] <- ifelse(is.na(Fish_Comm)[ , 5:(ncol(Fish_Comm))], 0, 1)
head(Fish_Comm)

# Make a vector of all the unique deployment codes. This will ensure the code makes an individual curve
# for each deployment. There should be 6 deployments. 
Deploys <- unique(Fish_Comm$Deployment)
Deploys <- Deploys[-c(7)]

# Build the data frame
fish.df <- data.frame()

# Getting the data for species accumulation curves. 
# Going through each deployment -- makes an individual curve for each different location. 
for(h in Deploys){
  df.deploy <- Fish_Comm[Fish_Comm$Deploy == h, ]
  column.sums <- colSums(df.deploy[5:(ncol(Fish_Comm))])
  # Keep it presence absence - all numbers should be 0 or 1
  column.sums2 <- ifelse(column.sums == 0, 0, 1)
  # Goes through every species for each deployment. 
  for(i in 1:nrow(df.deploy)){
    # Pulls a random amount of videos from 1-100 according to what number j is. And counts
    # the number of species observed. 
    for(j in 1:100){
      video.clips <- df.deploy[sample(1:nrow(df.deploy),i), ]
      column.sums3 <- colSums(video.clips[5:(ncol(Fish_Comm))])
      # Keeps everything 0 and 1
      column.sums4 <- ifelse(column.sums3 == 0, 0, 1)
      column.sums4 <- ifelse(column.sums4 == 0, 0, 1)
      num.fish.video<- sum(column.sums4)
      # put everything into the data frame made in the beginning. 
      fish.unique <- num.fish.video
      fish.vec <- cbind(h, i, fish.unique)
      fish.df <- rbind(fish.df, fish.vec)
      
    }
  }
}

head(fish.df)
tail(fish.df)
str(fish.df)

# Everything is filed as a factor so you must coerce it. 

# Make a data.frame with the first column as a character vector of deployment codes
# second column is numeric with how many video clips it used, 
# third column is the proportion of fish species captured in the iteration. 
data.unique <- data.frame("Deploy" = as.character(fish.df$h), 
    "Num.video.clip" = as.numeric(as.character(fish.df$i)), 
    "Unique.fish" = as.numeric(as.character(fish.df$fish.unique)))
# Check that coercion worked. 
head(data.unique)
str(data.unique)

# Aggregate it = take the mean of the repeated trials so that the data.frame is back
# to the same size you started with . 
agg_mean.unique  <-  aggregate(data.unique$Unique.fish,
    by = list(data.unique$Deploy, data.unique$Num.video.clip),
    FUN = mean)

# Rename the columns 
colnames(agg_mean.unique) <- c("Deployment", "Num.video.clip", "Unique.fish")
# Make a number of days column - This is in the event you'd like to switch the x-axis from number of 
# samples to number of days the lander was deployed. 
Number.days <- agg_mean.unique$Num.video.clip
Number.days <- Number.days/72

agg_mean.unique <- cbind(agg_mean.unique, Number.days)

agg_mean.unique
str(agg_mean.unique)

# Aggregate again but using the standard error so we can make 95% confidence
# intervals. 
agg_sd.unique  <-  aggregate(data.unique$Unique.fish,
        by = list(data.unique$Deploy, data.unique$Num.video.clip),
        FUN = sd)  
colnames(agg_sd.unique) <- c("Deployment", "Num.video.clip", "Unique.fish")
agg_sd.unique
str(agg_sd.unique)

# Graphing the species accumulation curves using ggplot2
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
rare.plot.Fish <- ggplot(data = agg_mean.unique, mapping = aes(x = Num.video.clip, y = Unique.fish, 
  group = Deployment)) +
  geom_line(aes(color = Deployment)) +
  labs(x = "Number of Samples", y = "Number of Unique Fish Species") +
  theme_bw() + scale_colour_manual(values = cbPalette) 

rare.plot.Fish

#### CalCOFI Comparison with Nanolander Data ####

library(data.table)     #packageVersion("data.table"): ‘1.12.8’
library(tidyverse)      #packageVersion("tidyverse"): ‘1.2.1’
library(cowplot)        #packageVersion("cowplot"): ‘0.9.4’
library(ggplot2)        #packageVersion("ggplot2"): ‘3.3.0’
library(dplyr)          #packageVersion("dplyr"): ‘0.8.1’

CalSta_93.3_28 = read.csv("CalSta_93.3_28.csv", header = TRUE)
CalSta_93.3_28$X = NULL

CalSta_93.3_28$Date = as.factor(CalSta_93.3_28$Date)
CalSta_93.3_28$DateTime <- as.POSIXct(CalSta_93.3_28$Date_Time_PST, format="%d-%b-%Y %H:%M:%S") #format time
CalSta_93.3_28$Identifier = paste(CalSta_93.3_28$Sta_ID, CalSta_93.3_28$Depth, sep = '_')

head(CalSta_93.3_28)
summary(CalSta_93.3_28) 
str(CalSta_93.3_28)

#Check for outliers or weird data and set outliers to NAs
layout(matrix(1:2, ncol = 2))

CalSta_93.3_28$Temp1[CalSta_93.3_28$Temp1 > 30] <- NA 
CalSta_93.3_28$Temp1[CalSta_93.3_28$Temp1 < 5] <- NA
CalSta_93.3_28$Temp2[CalSta_93.3_28$Temp2 < 5] <- NA

CalSta_93.3_28$Salt1_Corr[CalSta_93.3_28$Salt1_Corr < 25] <- NA
CalSta_93.3_28$Salt1_Corr[CalSta_93.3_28$Salt1_Corr > 40] <- NA
CalSta_93.3_28$Salt2_Corr[CalSta_93.3_28$Salt2_Corr < 25] <- NA
CalSta_93.3_28$Salt2_Corr[CalSta_93.3_28$Salt2_Corr > 40] <- NA

#Look at which columns have no NAs (TempAve, SaltAve_Corr, but Ox has many NAs)
#Need to create an additional column that takes the mean of the two oxygen measurements and removes NAs
CalSta_93.3_28 = CalSta_93.3_28 %>%
  rowwise() %>% mutate(AvgOx=mean(c(Ox1uM_StaCorr, Ox2uM_StaCorr), na.rm=TRUE)) %>%
  as.data.table()

#Even though TempAve exists, sometimes (Jan-2015, for example), the calculation is wrong
CalSta_93.3_28 = CalSta_93.3_28 %>%
  rowwise() %>% mutate(AvgTemp=mean(c(Temp1, Temp2), na.rm=TRUE)) %>%
  as.data.table()

#Do it for salinity too, just in case
CalSta_93.3_28 = CalSta_93.3_28 %>%
  rowwise() %>% mutate(AvgSal=mean(c(Salt1_Corr, Salt2_Corr), na.rm=TRUE)) %>%
  as.data.table()

#look for any outliers
layout(matrix(1:4, ncol = 2))
hist(CalSta_93.3_28$AvgSal)
hist(CalSta_93.3_28$AvgTemp)
hist(CalSta_93.3_28$AvgOx)
head(CalSta_93.3_28)
str(CalSta_93.3_28)

#Calculate the coefficient of variation for each depth 0-500 m for T and O2
CV_93.3_28 = CalSta_93.3_28%>%
  group_by(Depth) %>%
  summarize(meanT = mean(AvgTemp, na.rm = TRUE),
            meanOx = mean(AvgOx, na.rm = TRUE),
            stdevT = sd(AvgTemp, na.rm = TRUE),
            stdevOx = sd(AvgOx, na.rm = TRUE),
            minT = min(AvgTemp, na.rm = TRUE),
            maxT = max(AvgTemp, na.rm = TRUE),
            minOx = min(AvgOx, na.rm = TRUE),
            maxOx = max(AvgOx, na.rm = TRUE),
            n=n()) %>%
  mutate(CV_Temp = (stdevT/meanT)*100, 
         CV_Ox = (stdevOx/meanOx)*100,
         Range_Temp = maxT-minT,
         Range_Ox = maxOx-minOx) %>%
  as.data.frame()

layout(matrix(1:1, ncol = 1))
ggplot(CV_93.3_28, aes(x=n, y = Depth*-1)) + geom_point()

#### Figure 7A-D ####

Temp_Plot_93.3_28_2 = qplot(meanT, Depth, data=CV_93.3_28, na.rm = TRUE) + 
  scale_y_reverse() + geom_errorbarh(aes(xmax = meanT+(2*stdevT), xmin = meanT-(2*stdevT)), height = 0, color = "#CC0000") +
  geom_errorbarh(aes(xmax = meanT+stdevT, xmin = meanT-stdevT), height = 0, col = "#FF6666") +
  geom_hline(yintercept=c(100, 200, 300, 400), linetype="dashed", color = c("#7bccc4", "#4292c6", "#6a51a3", "#08519c")) +
  geom_point(aes(meanT, Depth), size = 1) + ylab("Depth (m)") + 
  xlab("Temperature (°C)") +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
Ox_Plot_93.3_28_2 = qplot(meanOx, Depth, data=CV_93.3_28, na.rm = TRUE) + 
  scale_y_reverse() + geom_errorbarh(aes(xmax = meanOx+(2*stdevOx), xmin = meanOx-(2*stdevOx)), height = 0, color = "#0072B2") + 
  geom_errorbarh(aes(xmax = meanOx+stdevOx, xmin = meanOx-stdevOx), height = 0, color = "#56B4E9") + 
  geom_hline(yintercept=c(100, 200, 300, 400), linetype="dashed", color = c("#7bccc4", "#4292c6", "#6a51a3", "#08519c")) +
  geom_point(aes(meanOx, Depth), size = 1) + ylab("Depth (m)") + 
  xlab("Oxygen Concentration (umol/kg)") +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
CV_Temp_93.3_28 = qplot(CV_Temp, Depth, data = CV_93.3_28, na.rm = TRUE, size=I(1)) + scale_y_reverse() + 
  xlab("CV Temp") + ylab("Depth (m)") +
  geom_hline(yintercept=c(100, 200, 300, 400), linetype="dashed", color = c("#7bccc4", "#4292c6", "#6a51a3", "#08519c")) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
CV_Ox_93.3_28 = qplot(CV_Ox, Depth, data = CV_93.3_28, na.rm = TRUE, size=I(1)) + scale_y_reverse() + 
  xlab("CV Ox") + ylab("Depth (m)") +
  geom_hline(yintercept=c(100, 200, 300, 400), linetype="dashed", color = c("#7bccc4", "#4292c6", "#6a51a3", "#08519c")) +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))
plot_grid(Temp_Plot_93.3_28_2, CV_Temp_93.3_28, Ox_Plot_93.3_28_2, CV_Ox_93.3_28, ncol=4)

#check how many unique dates
unique(CalSta_93.3_28$Study) #says 61 unique levels (first: 20-Oct-2003; 
    #last: 04-Nov-2019)
head(CalSta_93.3_28)
tail(CalSta_93.3_28)

#### Figure 7 E-H #### 

# Extracted all data for depths 100, 200, 300, and 400 m from CalSta_93.3_28 and 
  #saved them into a datasheet with the Nanolander datasets 

CalCOFI_Comp=read.csv("CalCOFI_Comp_Updated.csv", header=TRUE)

summary(CalCOFI_Comp)
str(CalCOFI_Comp)
Comp_100m = subset(CalCOFI_Comp, Deployment == "D5" | Deployment == "D10" | Deployment == "CalCOFI_100m")
Comp_200m = subset(CalCOFI_Comp, Deployment == "D4" | Deployment == "D3" | Deployment == "D6" | Deployment == "CalCOFI_200m")
Comp_300m = subset(CalCOFI_Comp, Deployment == "D7" | Deployment == "CalCOFI_300m")
Comp_400m = subset(CalCOFI_Comp, Deployment == "D8" | Deployment == "CalCOFI_400m")

Means_100 = Comp_100m %>% 
  group_by(Deployment) %>%
  summarise_all(list(mean), na.rm = TRUE)

Means_200 = Comp_200m %>% 
  group_by(Deployment) %>%
  summarise_all(list(mean), na.rm = TRUE)

Means_300 = Comp_300m %>% 
  group_by(Deployment) %>%
  summarise_all(list(mean), na.rm = TRUE)

Means_400 = Comp_400m %>% 
  group_by(Deployment) %>%
  summarise_all(list(mean), na.rm = TRUE)

Violinplot_100 <- ggplot(Comp_100m, aes(x=Deployment, y=Ox_umol, fill = Deployment)) + 
  geom_violin(trim=TRUE) +
  scale_fill_manual(values=c("#7bccc4", "#7bccc4", "#7bccc4")) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 2), geom = "pointrange", color = "white") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange") +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x=element_blank(), axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) + 
  ylab("Oxygen concentration (umol/kg)") +guides(fill=FALSE)

Violinplot_200 <- ggplot(Comp_200m, aes(x=Deployment, y=Ox_umol, fill = Deployment)) + 
  geom_violin(trim=TRUE) +
  scale_fill_manual(values=c("#4292c6", "#4292c6", "#4292c6", "#4292c6", "#4292c6")) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 2), geom = "pointrange", color = "white") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange") +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x=element_blank(), axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  ylab("Oxygen concentration (umol/kg)") + guides(fill=FALSE)

Violinplot_300 <- ggplot(Comp_300m, aes(x=Deployment, y=Ox_umol, fill = Deployment)) + 
  geom_violin(trim=TRUE) + scale_fill_manual(values=c("#6a51a3", "#6a51a3")) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 2), geom = "pointrange", color = "white") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange") +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x=element_blank(), axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  ylab("Oxygen concentration (umol/kg)") + guides(fill=FALSE)

Violinplot_400 <- ggplot(Comp_400m, aes(x=Deployment, y=Ox_umol, fill = Deployment)) + 
  geom_violin(trim=TRUE) + scale_fill_manual(values=c("#08519c", "#08519c")) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 2), geom = "pointrange", color = "white") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "pointrange") +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x=element_blank(), axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) +
  ylab("Oxygen concentration (umol/kg)") + guides(fill=FALSE)

plot_grid(Violinplot_100, Violinplot_200, Violinplot_300, Violinplot_400, cols = 2)

#### Figure 7 I-J #### 

CalSta_93.3_28_Depths_Combined = CalSta_93.3_28 %>%
  filter(Depth == 100 | Depth == 200 | Depth == 300 | Depth == 400) %>%
  as.data.frame()

CalSta_93.3_28_Depths_Combined$Depth_As_Factor = CalSta_93.3_28_Depths_Combined$Depth
CalSta_93.3_28_Depths_Combined$Depth_As_Factor = as.factor(CalSta_93.3_28_Depths_Combined$Depth_As_Factor)

Cal_D100m =read.csv("Cal_D100m.csv", header=TRUE)
Cal_D100m$DateTime <- as.POSIXct(Cal_D100m$Date_Time_PST, format="%d-%b-%Y %H:%M:%S") #format time
Cal_D200m =read.csv("Cal_D200m.csv", header=TRUE)
Cal_D200m$DateTime <- as.POSIXct(Cal_D200m$Date_Time_PST, format="%d-%b-%Y %H:%M:%S") #format time
Cal_D300m =read.csv("Cal_D300m.csv", header=TRUE)
Cal_D300m$DateTime <- as.POSIXct(Cal_D300m$Date_Time_PST, format="%d-%b-%Y %H:%M:%S") #format time
Cal_D400m =read.csv("Cal_D400m.csv", header=TRUE)
Cal_D400m$DateTime <- as.POSIXct(Cal_D400m$Date_Time_PST, format="%d-%b-%Y %H:%M:%S") #format time

#The timeseries for station 93.3 28 is:
#From and including: Thursday, October 20, 2003
#To, but not including Tuesday, November 4, 2019
#Result: 5859 days
#It is 5859 days from the start date to the end date, but not including the end date.
#Or 16 years, 15 days excluding the end date.
#Or 192 months, 15 days excluding the end date.
#Alternative time units
#5859 days can be converted to one of these units:
#506,217,600 seconds
#8,436,960 minutes
#140,616 hours
#837 weeks
#1605.21% of a common year (365 days)

Model_Ox = CalSta_93.3_28_Depths_Combined %>% 
  group_by(Depth_As_Factor) %>% 
  do({
    mod = lm(AvgOx ~ DateTime, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  }) %>%
  mutate(Slope_Day = Slope*(60*60*24)) %>%
  mutate(Slope_Year = Slope*(60*60*24*365)) %>%
  mutate(Overall_Change_Years = Slope_Year*16) %>%
  mutate(Overall_Change_Days = Slope_Day*5859) %>%
  as.data.frame()
Model_Ox

str(Cal_D100m)
Mod_Ox_100 = lm(AvgOx ~ DateTime, data = Cal_D100m)
summary(Mod_Ox_100) #p = 0.973
Mod_Ox_200 = lm(AvgOx ~ DateTime, data = Cal_D200m)
summary(Mod_Ox_200) #p = 0.314
Mod_Ox_300 = lm(AvgOx ~ DateTime, data = Cal_D300m)
summary(Mod_Ox_300) #p = 0.0004619***
Mod_Ox_400 = lm(AvgOx ~ DateTime, data = Cal_D400m)
summary(Mod_Ox_400) #p = 1.026e-07***

Ox_TS = ggplot(CalSta_93.3_28_Depths_Combined, aes(x=DateTime, y = AvgOx, group = Depth_As_Factor)) +
  geom_point(aes(color=Depth_As_Factor)) + 
  scale_color_manual(values=c("#7bccc4", "#4292c6", "#6a51a3", "#08519c")) +
  geom_abline(slope = Mod_Ox_100$coefficients[2], intercept = Mod_Ox_100$coefficients[1], color = c("#7bccc4"), linetype="dashed", size = 0.8) +
  geom_abline(slope = Mod_Ox_200$coefficients[2], intercept = Mod_Ox_200$coefficients[1], color = c("#4292c6"), linetype="dashed", size = 0.8) +
  geom_abline(slope = Mod_Ox_300$coefficients[2], intercept = Mod_Ox_300$coefficients[1], color = c("#6a51a3"), size = 0.8) +
  geom_abline(slope = Mod_Ox_400$coefficients[2], intercept = Mod_Ox_400$coefficients[1], color = c("#08519c"), size = 0.8) + 
  ylab("Oxygen concentration (umol/kg)") + xlab("Date") +
  #ggtitle("Oxygen Timeseries") + 
  theme(legend.position = "none") 
Ox_TS

Model_Temp = CalSta_93.3_28_Depths_Combined %>% 
  group_by(Depth_As_Factor) %>% 
  do({
    mod = lm(AvgTemp ~ DateTime, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  }) %>%
  mutate(Slope_Day = Slope*(60*60*24)) %>%
  mutate(Slope_Year = Slope*(60*60*24*365)) %>%
  mutate(Overall_Change_Years = Slope_Year*16) %>%
  mutate(Overall_Change_Days = Slope_Day*5859) %>%
  as.data.frame()
Model_Temp

Mod_Temp_100 = lm(AvgTemp ~ DateTime, data = Cal_D100m)
summary(Mod_Temp_100) #p = 0.404
Mod_Temp_200 = lm(AvgTemp ~ DateTime, data = Cal_D200m)
summary(Mod_Temp_200) #p = 0.7049
Mod_Temp_300 = lm(AvgTemp ~ DateTime, data = Cal_D300m)
summary(Mod_Temp_300) #p = 0.9174
Mod_Temp_400 = lm(AvgTemp ~ DateTime, data = Cal_D400m)
summary(Mod_Temp_400) #p = 0.4387

Temp_TS = ggplot(CalSta_93.3_28_Depths_Combined, aes(x=DateTime, y = AvgTemp, group=Depth_As_Factor)) +
  geom_point(aes(color=Depth_As_Factor)) +
  scale_color_manual(values=c("#7bccc4", "#4292c6", "#6a51a3", "#08519c")) +
  geom_abline(slope = Mod_Temp_100$coefficients[2], intercept = Mod_Temp_100$coefficients[1], color = c("#7bccc4"), linetype="dashed", size = 0.8) +
  geom_abline(slope = Mod_Temp_200$coefficients[2], intercept = Mod_Temp_200$coefficients[1], color = c("#4292c6"), linetype="dashed", size = 0.8) +
  geom_abline(slope = Mod_Temp_300$coefficients[2], intercept = Mod_Temp_300$coefficients[1], color = c("#6a51a3"), linetype="dashed", size = 0.8) +
  geom_abline(slope = Mod_Temp_400$coefficients[2], intercept = Mod_Temp_400$coefficients[1], color = c("#08519c"), linetype="dashed", size = 0.8) + 
  ylab("Temperature (°C)") + xlab("Date") +
  #ggtitle("Temperature Timeseries") + 
  theme(legend.position = "none") 
Temp_TS

plot_grid(Temp_TS, Ox_TS, ncol = 2)

#### Reproduce Supplement 1G Cross Shelf Structure #### 

# Only selected data from stations 093.3 028.0, 093.3 026.7, 093.3 030.0, 
    #093.3 035.0, 093.3 040.0 (the closest line to the nanolander deployments)
  #Only selected CalCOFI cruises that were during the same times as the 
    #Nanolander deployments: 1708SR, 1711SR, 1802SH, 1804SH

CalSta_SelectDates=read.csv("CalSta_SelectDates.csv", header=TRUE)
CalSta_SelectDates$X = NULL

CalSta_SelectDates$Date = as.factor(CalSta_SelectDates$Date)
CalSta_SelectDates$DateTime <- as.POSIXct(CalSta_SelectDates$Date_Time_PST, format="%d-%b-%Y %H:%M:%S") #format time
CalSta_SelectDates$Identifier = paste(CalSta_SelectDates$Sta_ID, CalSta_SelectDates$Depth, sep = '_')

head(CalSta_SelectDates)
summary(CalSta_SelectDates) 
str(CalSta_SelectDates)

#Check for outliers or weird data and set outliers to NAs
layout(matrix(1:2, ncol = 2))

CalSta_SelectDates$Temp1[CalSta_SelectDates$Temp1 > 30] <- NA 
CalSta_SelectDates$Temp1[CalSta_SelectDates$Temp1 < 5] <- NA
CalSta_SelectDates$Temp2[CalSta_SelectDates$Temp2 < 5] <- NA

CalSta_SelectDates$Salt1_Corr[CalSta_SelectDates$Salt1_Corr < 25] <- NA
CalSta_SelectDates$Salt1_Corr[CalSta_SelectDates$Salt1_Corr > 40] <- NA
CalSta_SelectDates$Salt2_Corr[CalSta_SelectDates$Salt2_Corr < 25] <- NA
CalSta_SelectDates$Salt2_Corr[CalSta_SelectDates$Salt2_Corr > 40] <- NA

#Look at which columns have no NAs (TempAve, SaltAve_Corr, but Ox has many NAs)
#Need to create an additional column that takes the mean of the two oxygen measurements and removes NAs
CalSta_SelectDates = CalSta_SelectDates %>%
  rowwise() %>% mutate(AvgOx=mean(c(Ox1uM_StaCorr, Ox2uM_StaCorr), na.rm=TRUE)) %>%
  as.data.table()

#Even though TempAve exists, sometimes (Jan-2015, for example), the calculation is wrong
CalSta_SelectDates = CalSta_SelectDates %>%
  rowwise() %>% mutate(AvgTemp=mean(c(Temp1, Temp2), na.rm=TRUE)) %>%
  as.data.table()

#Do it for salinity too, just in case
CalSta_SelectDates = CalSta_SelectDates %>%
  rowwise() %>% mutate(AvgSal=mean(c(Salt1_Corr, Salt2_Corr), na.rm=TRUE)) %>%
  as.data.table()

#look for any outliers
layout(matrix(1:4, ncol = 2))
hist(CalSta_SelectDates$AvgSal)
hist(CalSta_SelectDates$AvgTemp)
hist(CalSta_SelectDates$AvgOx)
head(CalSta_SelectDates)
str(CalSta_SelectDates)

CalSta_SelectDates_1708SR = CalSta_SelectDates %>%
  filter(Sta_ID != "093.3 026.7") %>%
  filter(Study == "1708SR") %>% as.data.frame() 

CalSta_SelectDates_1711SR = CalSta_SelectDates %>%
  filter(Sta_ID != "093.3 026.7") %>%
  filter(Study == "1711SR") %>% as.data.frame() 

CalSta_SelectDates_1802SH = CalSta_SelectDates %>%
  filter(Sta_ID != "093.3 026.7") %>%
  filter(Study == "1802SH") %>% as.data.frame() 

CalSta_SelectDates_1804SH = CalSta_SelectDates %>%
  filter(Sta_ID != "093.3 026.7") %>%
  filter(Study == "1804SH") %>% as.data.frame() 

library(calecopal)      #packageVersion("calecopal"): ‘0.1.0’

Temp_1708SR = ggplot(CalSta_SelectDates_1708SR, aes(x=AvgTemp, y = Depth, color=Sta_ID)) + scale_y_reverse() +
  geom_point() + scale_color_manual(values = cal_palette("sierra2")) + 
  theme(legend.position = "none") + xlab("") + ylab("") + xlim(5,23) 
Ox_1708SR = ggplot(CalSta_SelectDates_1708SR, aes(x=AvgOx, y = Depth, color=Sta_ID)) + scale_y_reverse() +
  geom_point() + scale_color_manual(values = cal_palette("sierra2")) +
  theme(legend.position = "none") + xlab("") + ylab("") + xlim(0,400)
plot_grid(Temp_1708SR, Ox_1708SR, ncol = 2)

Temp_1711SR = ggplot(CalSta_SelectDates_1711SR, aes(x=AvgTemp, y = Depth, color=Sta_ID)) + scale_y_reverse() +
  geom_point() + scale_color_manual(values = cal_palette("sierra2")) +
  theme(legend.position = "none") + xlab("") + ylab("") + xlim(5,23)
Ox_1711SR = ggplot(CalSta_SelectDates_1711SR, aes(x=AvgOx, y = Depth, color=Sta_ID)) + scale_y_reverse() +
  geom_point() + scale_color_manual(values = cal_palette("sierra2")) +
  theme(legend.position = "none") + xlab("") + ylab("") + xlim(0,400)
plot_grid(Temp_1711SR, Ox_1711SR, ncol = 2)

Temp_1804SH = ggplot(CalSta_SelectDates_1804SH, aes(x=AvgTemp, y = Depth, color=Sta_ID)) + scale_y_reverse() +
  geom_point() + scale_color_manual(values = cal_palette("sierra2")) +
  theme(legend.position = "none") + xlab("") + ylab("") + xlim(5,23)
Ox_1804SH = ggplot(CalSta_SelectDates_1804SH, aes(x=AvgOx, y = Depth, color=Sta_ID)) + scale_y_reverse() +
  geom_point() + scale_color_manual(values = cal_palette("sierra2")) +
theme(legend.position = "none") + xlab("") + ylab("") + xlim(0,400)
plot_grid(Temp_1804SH, Ox_1804SH, ncol = 2)

plot_grid(Temp_1708SR, Temp_1711SR, Temp_1804SH, Ox_1708SR, Ox_1711SR, Ox_1804SH, ncol = 3)