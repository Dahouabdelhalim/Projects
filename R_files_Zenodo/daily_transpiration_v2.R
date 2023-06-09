#*******************************************************************************
# Cesar D. Jimenez-Rodriguez
# Luxembourg Institute of Science and Technology (LIST)
# Environmental Research and Innovation Department (ERIN)
# Environmental Sensing and Modelling Unit
# Agro-Environmental Systems Group
# E-mail: cesar.jimenez@list.lu
#*******************************************************************************

#-Script Description-----------------------------------------------------------#
# This script read the sapflow data from SapFluxNet (Poyatos et al. 2021) [doi: 10.5194/essd-13-2607-2021].
# It uses the files in csv format provided by Poyatos et al. (2021) and follows 
# the recommendations of Nelson et al. (2020) [doi: 10.1111/gcb.15314] for upscaling 
# sapflow measurements to stand level. These recommendations were applied to hourly 
# sampling intervals instead of subhourly sampling intervals. The main output consist of
# a csv file with daily transpiration at stand level. The output file is composed
# by the following columns:
#
# "time"    : date format "YYYY-MM-DD"
# "stand"   : mean daily transpiration in mm/d
# "sd_stand": standard deviation of the mean daily transpiration in mm/d
# "n_trees" : number of trees used for the estimation of daily transpiration
#
# The files are saved as sap_SITE_NAME.csv
#
# NOTE: Sap_csv_v3:202108#

#-Cleaning R environment-------------------------------------------------------- 
rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
in_path  = "/0.1.5 SapFlux/csv/plant"    # Folder path where the input files from SapFluxNet are stored.
out_path =' '                     # Folder path to store the output file.

#-Select ID of the SapFluxNet site --------------------------------------------#
site_id = 'RUS_FYO'        # SAPFLUXNET code used for each experimental location | Examples: 'FRA_HES_HE2_NON' ; 'FRA_PUE' ; "ESP_ALT_TRI" ; RUS_Fyo'

#### DO NOT CHANGE FROM THIS POINT ONWARDS #####################################
#-Required Packages-------------------------------------------------------------
library(dplyr)

#-Data Retrieval ---------------------------------------------------------------
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

#-Sampling Times & Hourly Summary-----------------------------------------------
tz <- as.numeric(sub(".*UTC *(.*?) *:.*", "\\\\1", env$env_time_zone))
tznc <- sub(".*UTC *(.*?) *,.*", "\\\\1", env$env_time_zone)

sap$TIMESTAMP <- gsub("T"," ", sap$TIMESTAMP)
sap$TIMESTAMP <- gsub("Z"," ", sap$TIMESTAMP)
sap$time <- paste(substr(sap$TIMESTAMP,0,13),":00:00",sep="")
sap <- subset(sap, select = -c(solar_TIMESTAMP,TIMESTAMP))
sap <- sap %>% group_by(sap$time) %>% summarise_all(list(mean),na.rm = TRUE)
sap$time <- NULL
colnames(sap)[1] <- "time"
sap$time <- substr(sap$time,0,10)

sap_t <- sap %>% group_by(sap$time) %>% summarise_at(vars(-time),sum,na.rm=TRUE)
colnames(sap_t)[1] <- "time"
sap_n <- sap
sap_n[c(2:ncol(sap_n))] <- ifelse(is.na(sap_n[c(2:ncol(sap_n))]),1,0)
sap_n <- sap_n %>% group_by(sap_n$time) %>% summarise_at(vars(-time),sum,na.rm=TRUE) #cm3/tree/d
colnames(sap_n)[1] <- "time"

trs <- unique(trees$pl_code)
dts <- unique(sap_n$time) 
for (i in trs) {
  x <- sap_n[c("time",as.character(i))]
  y <- sap_t[c("time",as.character(i))]
  z <- merge(x,y,by="time");colnames(z) <- c("time","n","t")
  for (j in dts){
    z$t[z$time==j] <- ifelse(z$n[z$time==j] > 4, -1000,  z$t[z$time==j])
  }
  sap_t[as.character(i)] <- z$t
  }

sap_t[sap_t == -1000] <- NA

treex <- trees[c("pl_code","pl_dbh")]
treex$g <- pi*(treex$pl_dbh/2)^2; treex$pl_dbh <-NULL
cnms <- colnames(sap_t);cnms <- cnms[c(-1)]
for (i in cnms) {
  sap_t[i] <- (sap_t[i]/treex$g[treex$pl_code==i])*(stand$st_basal_area/10000)*10
} #mm/tree/d

sap_t$tree <- apply(sap_t[,2:(nrow(trees))], 1, mean,na.rm = TRUE) # daily average
sap_t$tree_sd <- apply(sap_t[,2:(nrow(trees))], 1, sd,na.rm = TRUE) # Hourly average
sap_t$n <- rowSums(!is.na(sap_t[c(2:length(cnms))]))

sap <- sap_t[c("time","tree","tree_sd","n")]
colnames(sap) <- c("time","stand","sd_stand","n_trees")

# Path to folder ---------------------------------------------------------------
dir.create(paste(out_path,"/",site_id,sep=""))
write.csv(sap,file = paste(out_path,"/",site_id,"/","sap_",site_id,".csv",sep=""))

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#
