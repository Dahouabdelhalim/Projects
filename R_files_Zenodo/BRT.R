### R Code to run a boosted regression tree (BRT) model
### Created by: Hannah Clipp


### LIBRARIES

library(mgcv)
library(dismo)
library(gbm)
library(foreign)
library(stringr)



### LOAD DATA

OUT<-"[Insert file path to desired data file]"
setwd(OUT)
winddata <- readRDS("winddata_wk.rds")

#for monthly and weekly scales
winddata <- subset(winddata, winddata$rosariwmea < 1000)
winddata$logVIR <- log(winddata$rosariwmea + 0.001)

#for daily scales
winddata <- subset(winddata, winddata$ros < 1000)
winddata$logVIR <- log(winddata$ros + 0.001)

winddata$ndvi[is.nan(winddata$ndvi)] <- NA

#check distribution
plot(winddata$Longitude,winddata$Latitude,pch='.')



### BRT MODEL TIME

brt_model  <- gbm.step(data=winddata, 
                       gbm.x = c("Longitude","Dist_coast","ha_5","ag_5","ur_5","Rel_elevat","Dist_Radar","distDMSP12","temp","NDVI","u_ful_925","v_ful_925","u_Car_925","v_Car_925","u_Atl_925","v_Atl_925","ex_uwind","ex_vwind","year","week"),
                       gbm.y = "logVIR",
                       family="gaussian",
                       tree.complexity = 2, #keep at 2
                       learning.rate = 0.5,
                       bag.fraction = 0.5,n.trees=100) #keep n.trees at 100; adjust learning.rate if needed

saveRDS(brt_model,"brt_wind_wk.rds")
