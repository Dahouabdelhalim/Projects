# Created by MRU on May 2022
# To do the biomass predictions with uncertainty quantification

# Libraries: --------------------------------------------------------------
library(raster)
library(tidyr)
library(ncdf4)
library(dplyr)
library(stringr)
library(tidyverse)
library(Hmisc)
library(rgdal)
library(rgeos)
library(boot)
library(bayestestR)
library(parallel)
library(ranger)

# Directories: ------------------------------------------------------------
biomassdir = "/data/homezvol1/uribedim/biomass/"

# Area to test: -----------------------------------------------------------

# extent to crop data to:
myarea <- extent(-113, -30, -23.5, 0)    #neotrops: -113, -30; africa: -20, 65; asia: 65, 180 
myareaasspoly <- as(myarea, 'SpatialPolygons')
crs(myareaasspoly) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Load Biomass data: -----------------------------------------------------------

biomassAvg_Bacc <- read_csv(paste0(biomassdir, "biomassavgBacc_allaimcwd_19802009_nohumanwater.csv"), 
                            col_names = FALSE) #%>%
colnames(biomassAvg_Bacc) = c("filename", "biomass_avg", "biomass_low", "biomass_upp", "biomass_med", 
                              "biomass_lowq", "biomass_upq", "biomass_sd", "pixelcount")
biomassAvg_Bacc <- biomassAvg_Bacc %>%
  mutate(ai = as.factor(unlist(lapply(biomassAvg_Bacc$filename,function(x) str_split(x, "_")[[1]][4]))))

biomassAvg_Xu <- read_csv(paste0(biomassdir, "biomassavgXu_allaimcwd_19802009_nohumanwater.csv"), 
                          col_names = FALSE)
colnames(biomassAvg_Xu) = c("filename", "biomass_avg", "biomass_low", "biomass_upp", "biomass_med", 
                            "biomass_lowq", "biomass_upq", "biomass_sd", "pixelcount")
biomassAvg_Xu <- biomassAvg_Xu %>%
  mutate(ai = as.factor(unlist(lapply(biomassAvg_Xu$filename,function(x) str_split(x, "_")[[1]][4]))))

biomassAvg_ESA <- read_csv(paste0(biomassdir, "biomassavgESA_allaimcwd_19802009_nohumanwater.csv"), 
                           col_names = FALSE) 
colnames(biomassAvg_ESA) = c("filename", "biomass_avg", "biomass_low", "biomass_upp", "biomass_med", 
                             "biomass_lowq", "biomass_upq", "biomass_sd", "pixelcount")
biomassAvg_ESA <- biomassAvg_ESA %>%
  mutate(ai = as.factor(unlist(lapply(biomassAvg_ESA$filename,function(x) str_split(x, "_")[[1]][4]))))

biomassAvg_all = list(biomassAvg_Bacc, biomassAvg_Xu, biomassAvg_ESA)

# Load Models: -----------------------------------------------------------------

rfrg_Bacc <- readRDS(paste0(biomassdir, "rfrg_Bacc2.rds"))
rfrg_Xu <- readRDS(paste0(biomassdir, "rfrg_Xu2.rds"))
rfrg_ESA <- readRDS(paste0(biomassdir, "rfrg_ESA2.rds"))
rfrg_all = list(rfrg_Bacc, rfrg_Xu, rfrg_ESA)

# Load Climate data: ----------------------------------------------

## GCMs:
### Precipitation data:
i.gcmP = stack(paste0(biomassdir, "Global_Mgcm_P.30yr.nc"),
               varname = "rcp85mean")
### MCWD data:
i.gcmMCWD.pm = stack(paste0(biomassdir, "Global_Mgcm_MCWDpm.30yr.nc"),
                     varname = "rcp45mean")

### Inclunding CRU lower and upper Q analysis:
i.gcmP.low = raster::stack(paste0(biomassdir, "Global_Mgcm_P.30yr_rcp45cru_lowQ.nc"))
i.gcmMCWD.low = raster::stack(paste0(biomassdir, "Global_Mgcm_MCWDpm.30yr_rcp45cru_lowQ.nc"))

i.gcmP.up = raster::stack(paste0(biomassdir, "Global_Mgcm_P.30yr_rcp45cru_upQ.nc"))
i.gcmMCWD.up = raster::stack(paste0(biomassdir, "Global_Mgcm_MCWDpm.30yr_rcp45cru_upQ.nc"))

# Estimate biomass per climate zone for each timeframe, using biomass and climate uncertainty (distribution): -------------------------------------------

# Function to find the corresponding AI zone from P and MCWD data:
findai <- function(x){
  P = x[1]
  MCWD = x[2]
  theai <- P/MCWD 
  theai2 <- ifelse(theai < -3.8 | theai == Inf, 6, 
                   ifelse(theai >= -3.8  & theai < -1.8, 5, 
                          ifelse(theai >= -1.8  & theai < -1, 4,
                                 ifelse(theai >= -1  & theai < -0.25, 3,
                                        ifelse(theai >= -0.25  & theai < -1/19, 2,
                                               ifelse(theai >= -1/19  & theai <= 0, 1, NA))))))
  theai2 <- ifelse(theai2 == 6 & P >= 1700, 7, theai2)
  return(theai2)
}

# Function to predict biomass quantiles by pixel using RF quantile regression models and sample from that distribution:
preditcbio_RF <- function(climraster, themodel){
  tempdf <- as.data.frame(climraster, xy = TRUE, na.rm = TRUE)
  theprediction <- predict(themodel,
                           data = data.frame(P = tempdf[,3], MCWD = tempdf[, 4]),
                           predict.all = FALSE, type = "quantiles", quantiles = c(0.25, 0.5, 0.75))
  preddf <- data.frame(x = tempdf[,1], y = tempdf[,2],
                       avg = theprediction$predictions[,2], 
                       lowq = theprediction$predictions[,1],
                       upq = theprediction$predictions[,3])
  sd1RF <- abs((preddf$lowq - preddf$avg)/0.675)
  sd2RF <- abs((preddf$upq - preddf$avg)/0.675)
  preddf$RFsd <- rowMeans(cbind(sd1RF, sd2RF)) #or sd.final=max(c(sd1,sd2))
  preddf$pred <- apply(preddf, 1, 
                       function(x) rnorm(1, mean =  x[3], sd = x[6]))
  preddf$pred <- ifelse(preddf$pred < 0, 0, preddf$pred)
  return(preddf)
}

# Crop climate data to desired extent:
i.gcmP.sub <- crop(i.gcmP, myareaasspoly)
i.gcmPlow.sub <- crop(i.gcmP.low, myareaasspoly)
i.gcmPup.sub <- crop(i.gcmP.up, myareaasspoly)
i.gcmMCWD.sub <- crop(i.gcmMCWD.pm, myareaasspoly)
i.gcmMCWDlow.sub <- crop(i.gcmMCWD.low, myareaasspoly)
i.gcmMCWDup.sub <- crop(i.gcmMCWD.up, myareaasspoly)

## Function to predict biomass with random sampling of climate, dataset and biomass density:
clim_unc_wbio <- function(precdata, mcwddata, 
                          precup, preclow, mcwdlow, mcwdup, 
                          biomassdatasets, biomassmodels, 
                          method){

  # Climate distribution:
  sd1P <- abs((precup - precdata)/0.675)
  sd2P <-  abs((preclow - precdata)/0.675)
  sdP.final <- calc(stack(sd1P, sd2P), mean, na.rm = TRUE) #or sd.final=max(c(sd1,sd2))
  Pdata <- stack(precdata, sdP.final)
  sd1MCWD <- abs((mcwdup - mcwddata)/0.675)
  sd2MCWD <-  abs((mcwdlow - mcwddata)/0.675)
  sdMCWD.final <- calc(stack(sd1MCWD, sd2MCWD), mean, na.rm = TRUE) #or sd.final=max(c(sd1,sd2))
  MCWDdata <- stack(abs(mcwddata), sdMCWD.final)
  
  # Climate sampling:
  newPdata <- calc(Pdata, fun = function(w){rnorm(1, mean = w[1], sd = w[2])})
  newPdata[newPdata < 0] = 0
  newMCWDdata <- calc(MCWDdata, fun = function(z){rnorm(1, mean = z[1], sd = z[2])}) * -1
  newMCWDdata[newMCWDdata > 0] = 0
  
  # Find AI for the sampled climate:
  newclimdata <- stack(newPdata, newMCWDdata)
  theai <- calc(newclimdata, fun = findai)
  
  # Choose one biomass dataset:
  thedatasetindex = sample(c(1:length(biomassdatasets)), 1)
  biomassdata <- biomassdatasets[[thedatasetindex]]
  #biomassmodel <- biomassmodels[[thedatasetindex]]
  
  # Biomass sampling:
  if(method == "Avg"){
    tempras <- theai
    thesds <- lapply(c(1:7), function(ai){
      print(ai)
      sd1ai <- abs((biomassdata$biomass_upp[ai] - biomassdata$biomass_avg[ai])/1.96)
      sd2ai <- abs((biomassdata$biomass_low[ai] - biomassdata$biomass_avg[ai])/1.96)
      sdai.final <- mean(c(sd1ai,sd2ai))
      print(length(tempras[tempras == ai]))
      #if(length(tempras[tempras == ai]) > 0){
        #print(tempras)
        #tempras[tempras == ai] <- rnorm(length(tempras[tempras == ai]), mean = biomassdata$biomass_avg[ai], 
        #                                sd = sdai.final)
        #print(tempras)
      return(sdai.final)
    })
    tempras[tempras == 1] <- rnorm(length(tempras[tempras == 1]), mean = biomassdata$biomass_avg[1], 
                                   sd = thesds[[1]])
    tempras[tempras == 2] <- rnorm(length(tempras[tempras == 2]), mean = biomassdata$biomass_avg[2], 
                                   sd = thesds[[2]])
    tempras[tempras == 3] <- rnorm(length(tempras[tempras == 3]), mean = biomassdata$biomass_avg[3], 
                                   sd = thesds[[3]])
    tempras[tempras == 4] <- rnorm(length(tempras[tempras == 4]), mean = biomassdata$biomass_avg[4], 
                                   sd = thesds[[4]])
    tempras[tempras == 5] <- rnorm(length(tempras[tempras == 5]), mean = biomassdata$biomass_avg[5], 
                                   sd = thesds[[5]])
    tempras[tempras == 6] <- rnorm(length(tempras[tempras == 6]), mean = biomassdata$biomass_avg[6], 
                                   sd = thesds[[6]])
    tempras[tempras == 7] <- rnorm(length(tempras[tempras == 7]), mean = biomassdata$biomass_avg[7], 
                                   sd = thesds[[7]])
    print(tempras)
    }
  
  if(method == "RF"){
    RFpred <- preditcbio_RF(newclimdata, 
                            themodel = biomassmodels[[thedatasetindex]])
    print("predictRF ok")
    print(RFpred)
    tempras <- rasterFromXYZ(cbind(RFpred$x, RFpred$y, RFpred$pred))
    print(tempras)
  }
  
  return(tempras) # output is a raster with biomass predictions for each pixel
}

# Run function (with 1000 replicates/simulations) for each of the five time frames:
start_time <- Sys.time()
clim_unc_alltimes <- mclapply(c(1:5), function(x){
  print(x)
  reps = 1000
  done = 0
  themethod = "Avg"
  clim_unc_rep <- replicate(reps,  clim_unc_wbio(precdata = i.gcmP.sub[[x]], 
                                                mcwddata = i.gcmMCWD.sub[[x]], 
                                                preclow = i.gcmPlow.sub[[x]], 
                                                precup = i.gcmPup.sub[[x]], 
                                                mcwdlow = i.gcmMCWDlow.sub[[x]],
                                                mcwdup = i.gcmMCWDup.sub[[x]],
                                                biomassdatasets = biomassAvg_all,
                                                biomassmodels = rfrg_all,
                                                method = themethod))
  # Write the result of each replicate to a raster file so that I can use it later:
  y = done + c(1:reps)
  print(y)
  lapply(c(1:reps), function(z) writeRaster(clim_unc_rep[[z]], 
                                            filename = paste0("./allreps/rcp45/avg/neosouth2/biomass_estwunc_rcp45_", themethod, "_neosouth_alldata_", x, "_rep_", y[z],".tif")))
  return(clim_unc_rep)
}, mc.cores = 5)
end_time <- Sys.time()
print(end_time - start_time)

