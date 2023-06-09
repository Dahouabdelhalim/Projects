#*********************************************************************************
# Cesar D. Jimenez-Rodriguez
# Luxembourg Institute of Science and Technology (LIST)
# Environmental Research and Innovation Department (ERIN)
# Environmental Sensing and Modelling Unit
# Agro-Environmental Systems Group
# E-mail: cdjimenezcr@gmail.com
#         cesar.jimenez@list.lu 
#*********************************************************************************
#
#-Script Description------------------------------------------------------------
# This script read the sapflow data from SapFluxNet (Poyatos et al. 2021) [doi: 10.5194/essd-13-2607-2021].
# It uses the files in csv format provided by Poyatos et al. (2021) and follows 
# the recommendations of Nelson et al. (2020) [doi: 10.1111/gcb.15314] for upscaling 
# sapflow measurements to stand level. These recommendations were applied to hourly 
# sampling intervals instead of sub-hourly sampling intervals. The main output consist of
# a csv file with daily transpiration for a selected species. The output file is composed
# by the following columns:
#
# "time"    : date format "YYYY-MM-DD"
# "stand"   : mean daily transpiration in mm/d
# "sd_stand": standard deviation of the mean daily transpiration in mm/d
# "n_trees" : number of trees used for the estimation of daily transpiration
#
# The files are saved as sap_SITE_NAME.csv
#
# NOTE: Sap_csv_v6_xsps:202211#

#-Cleaning R environment-------------------------------------------------------- 
rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
in_path  = "./SAPFLUXNET/" # Folder path where the input files from SapFluxNet are stored.
out_path = "./results/"                     # Folder path to store the output file.

#-Select ID of the SapFluxNet site --------------------------------------------#
#site_id = 'ESP_ALT_TRI'; stid <- 'ES-Alt' ; pftn=6 # NOTE: for DE-Hin, you have to switch line 115 for 116
#site_id = 'FRA_PUE'; stid <- 'FR-Pue' ; pftn=6 # NOTE: for DE-Hin, you have to switch line 115 for 116
#sps ="Quercus ilex"

#site_id = 'FRA_HES_HE2_NON'; stid <- 'FR-Hes' ; pftn=8 # NOTE: for DE-Hin, you have to switch line 115 for 116
site_id = 'DEU_HIN_TER'; stid <- 'DE-Hin' ; pftn=8 # NOTE: for DE-Hin, you have to switch line 115 for 116
sps ="Fagus sylvatica"

#### DO NOT CHANGE FROM THIS POINT ONWARDS #####################################
#-Required Packages-------------------------------------------------------------
library(dplyr)
library(ncdf4)
library(lubridate)

# Extra Functions --------------------------------------------------------------
VPD <- function(ta,rh){
  es   <- 0.6108*exp((17.27*ta)/(ta+237.3))
  ea   <- es*rh/100
  vpd  <- es-ea
  return(vpd)
}

#-Data Retrieval ----------------------------------------------------------------------------------
setwd(in_path)
env <- read.csv(paste(site_id,"_env_md.csv",sep = ""))
site <- read.csv(paste(site_id,"_site_md.csv",sep = ""))
sap <- read.csv(paste(site_id,"_sapf_data.csv",sep = ""))
trees <- read.csv(paste(site_id,"_plant_md.csv",sep = ""))
stand <- read.csv(paste(site_id,"_stand_md.csv",sep = ""))

cnms <- colnames(sap)
cnms[2] <- "solar_TIMESTAMP"
colnames(sap) <- cnms

rm(cnms)

# MODEL RESULTS ----------------------------------------------------------------

#-Sampling Times & Hourly Summary----------------------------------------------------------------------------
tz <- as.numeric(sub(".*UTC *(.*?) *:.*", "\\\\1", env$env_time_zone))
tznc <- sub(".*UTC *(.*?) *,.*", "\\\\1", env$env_time_zone)

sap$TIMESTAMP <- gsub("T"," ", sap$TIMESTAMP)
sap$TIMESTAMP <- gsub("Z"," ", sap$TIMESTAMP)
sap$time <- paste(substr(sap$TIMESTAMP,0,13),":00:00",sep="")
sap <- subset(sap, select = -c(solar_TIMESTAMP,TIMESTAMP))
sap <- sap %>% group_by(sap$time) %>% summarise_all(list(mean),na.rm = TRUE)
sap$time <- NULL
colnames(sap)[1] <- "time"

trs <- unique(trees$pl_code)

treex <- trees[c("pl_code","pl_dbh","pl_species")]
treex$g <- pi*(treex$pl_dbh/2)^2 #m2
treex$pl_dbh <-NULL

# Done to remove the species 
sps_lst <- NA
for (k in treex$pl_code) {
 if(treex$pl_species[treex$pl_code==k] == sps){
   next
 }else{
   sps_lst <- append(sps_lst, k)
 }
}
sps_lst <- na.omit(sps_lst)

sap[sps_lst] <- NULL
sap[sap == "NaN"] <- NA

spt <- data.frame(trees$pl_code)
spt <- data.frame(spt[!spt$trees.pl_code %in% sps_lst,]);colnames(spt)<-"sps"
n_trees <- as.numeric(length(spt$sps))

#Transpiration per site in mm hr
sap_st <- sap
for (i in spt$sps) {
  #sap_st[i] <- (sap_st[i]/treex$g[treex$pl_code==i])*(stand$st_basal_area/10000)*10
  sap_st[i] <- (sap_st[i]/treex$g[treex$pl_code==i])*(37.33/10000)*10 # Data based on the paper https://rmets.onlinelibrary.wiley.com/doi/10.1002/gdj3.45
} #mm/tree/d
sap_st$et_mm <- apply(sap_st[,2:n_trees], 1, mean,na.rm = TRUE) # hourly average
sap_st <- sap_st[c("time","et_mm")]
sap_st[sap_st == "NaN"] <- NA

sap <- sap_st

sap <- sap[!(format(strptime(sap$time,format = "%Y-%m-%d %H:%M:%S"),"%m") == "02" & format(strptime(sap$time,format = "%Y-%m-%d %H:%M:%S"), "%d") == "29"), , drop = FALSE] # With this linee we remove the leap day of leap year
sap$date <- date(strptime(sap$time,format = "%Y-%m-%d %H:%M:%S"))

sap$min <- NA
for (t in unique(sap$date)) {
  sapmin <- na.omit(sap$et_mm[sap$date==t])
  sapmin <- min(sapmin)
  sap$min[sap$date==t] <- sapmin
}

sap$min <- replace(sap$min, is.infinite(sap$min),NA)
sap$et_mm <- sap$et_mm-sap$min

sap[c("min","date")]<-NULL

sap_n <- sap[c("time","et_mm")]
sap_n$time <- substr(sap$time,0,10)

sap_n <- sap_n %>% group_by(sap_n$time) %>% summarise_at(vars(-time),sum,na.rm=TRUE) #mm/d

####################
sap <- sap_n; colnames(sap) <- c("time","stand")
sap$doy <- yday(as.Date(paste(1993,substr(sap$time,6,7),substr(sap$time,9,10),sep="-")))
sap$YEAR <- year(strptime(sap$time,format = "%Y-%m-%d"))

# Path to folder ---------------------------------------------------------------
write.csv(sap,file = paste(out_path,"sap_",stid ,".csv",sep=""))

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#