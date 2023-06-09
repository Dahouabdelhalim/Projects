#*******************************************************************************
#
# Cesar D. Jimenez-Rodriguez
# Luxembourg Institute of Science and Technology (LIST)
# Environmental Research and Innovation Department (ERIN)
# Environmental Sensing and Modelling Unit
# Agro-Environmental Systems Group
# E-mail: cesar.jimenez@list.lu | cdjimenezcr@gmail.com
# 
#*******************************************************************************
#
#-Script Description------------------------------------------------------------
# This script prepares the atmospheric forcing required by CLM 5.0. based on the
# available data from SapFluxNet (Poyatos et al. 2020) [doi: 10.5281/zenodo.2530797].
# All Data is summarized into hourly time steps to create the ncdf file. 
# The previously identified data gaps are filled with the hourly data available 
# in "COSMO Regional Reanalysis" (ftp://opendata.dwd.de/climate_environment/REA/COSMO_REA6/).
# All NetCDF files are saved as yyyy-mm.nc
# NOTE: SapFLuxNCDF_csv_v6:202103
#------------------------------------------------------------------------------#
#-Cleaning the environment------------------------------------------------------
rm(list=ls()) # Clean the R environment
dev.off()
graphics.off()
#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
in_path  = "~/Documents/CAPACITY/DataSets/0.1.5 SapFlux/csv/plant" # Location of Sapfluxnet data set.
out_path = '~/Documents/CAPACITY/Outputs/prep' # Location of NetCDF files.
gap_path = "~/Documents/CAPACITY/Outputs/gaps/" # Location of data retrieved from COSMO REA for filling the gaps.

#-Select ID of the SapFluxNet site --------------------------------------------#
site_id = 'FRA_HES_HE2_NON'; tw = c("2001-01-01","2005-12-31")

#-Missing variables to be filled with the COSMO REA product.
mv = c("ps","rf") #full missing variables on the data set: ps:atmospheric pressure, rf: precipitation, ta: air temperature

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#### DO NOT CHANGE FROM THIS POINT ONWARDS #####################################
#-Required Packages-------------------------------------------------------------
library(ncdf4)
library(sapfluxnetr)
library(dplyr)
library(lubridate)
library(tidyverse)
library(rgdal)
library(sp)
library(ggplot2)
#-Data Retrieval ---------------------------------------------------------------
setwd(in_path)
atm <- read.csv(paste(site_id,"_env_data.csv",sep = ""))
env <- read.csv(paste(site_id,"_env_md.csv",sep = ""))
site <- read.csv(paste(site_id,"_site_md.csv",sep = ""))
stand <- read.csv(paste(site_id,"_stand_md.csv",sep = ""))

cnms <- colnames(atm)
cnms[2] <- "solar_TIMESTAMP"
colnames(atm) <- cnms

rm(cnms)

#-Sampling Times & Hourly Summary-----------------------------------------------
tw <- as.POSIXct(tw)

tz <- as.numeric(sub(".*UTC *(.*?) *:.*", "\\\\1", env$env_time_zone))
tznc <- sub(".*UTC *(.*?) *,.*", "\\\\1", env$env_time_zone)

atm$TIMESTAMP <- gsub("T"," ", atm$TIMESTAMP)
atm$TIMESTAMP <- gsub("Z",paste(" :",tz,sep = ""), atm$TIMESTAMP)
atm$time <- paste(substr(atm$TIMESTAMP,0,13),":00:00",sep="")

rf <- atm[c("time","precip")]
rf <- rf %>% group_by(rf$time) %>% summarise_at(vars(-time), sum)
colnames(rf)<-c("tm","precip")
atm <- subset(atm, select = -c(solar_TIMESTAMP,TIMESTAMP,precip))
atm <- atm %>% group_by(atm$time) %>% summarise_all(list(mean),na.rm = TRUE)
colnames(atm)[1] <- c("tm")
atm <- merge(rf,atm,by="tm")
atm$date <- paste(substr(atm$tm,0,10),sep="")
rm(rf)

#-Computing missing atmospheric forcing variables-------------------------------
# Incomming longwave Radiation
lwin <- function(ta,rh){
  e   <- 0.6108*exp((17.27*ta)/(ta+237.3))
  ecs <- 0.7+(e*(5.95*10^(-4))*exp(1500/(ta+273.15)))
  lw  <- (5.67*10^(-8))*ecs*(ta+273.15)^4
  return(lw)
} # This formula retrieve the incomming long wave radiation with a clear sky emisivity based on Idso (1981) [doi: 10.1029/WR017i002p00295].

#-Path to folder storing the data in nc format ---------------------------------
dir.create(paste(out_path,paste("nc_",site_id,sep = ""),sep="/"))
ncpath <- paste(out_path,paste("nc_",site_id,sep = ""),sep="/")
ncname <- site_id

#- The file name is assigned to each group of variables within the file.
atm$fname <- paste(substr(atm$tm,0,7),sep="")

#-Creating final ncfiles -------------------------------------------------------
#-Retrieving all data required from SapFLUXNET Files
atmd <- atm
atmd <- atmd[(between(as.numeric(strptime(paste(atmd$fname,"-01",sep=""),format = "%Y-%m-%d")),
                     as.numeric(tw[1]),
                     as.numeric(tw[2]))),]
atmd <- atmd[!is.na(atmd$fname),]
atmd$rn <- as.numeric(row_number(atmd$tm))-1
atmd$sec <- 0
atmd$sec <- ifelse(atmd$rn == 0,
                   as.numeric(strptime(atmd$tm,format = "%Y-%m-%d %H:%M:%S")),
                   as.numeric(strptime(head(atmd$tm,1),format = "%Y-%m-%d %H:%M:%S"))+(3600*atmd$rn))

x<-atmd
x %>% tally()
x <- x%>% group_by(fname) %>% tally()

#-Retrieveing data for the missing variables from COSMO-REA6 -------------------

for(i in mv){
  u <- read.csv(paste(gap_path,"csv_",site_id,"/",site_id,"_gap_",i,".csv",sep = ""));
  u <- u[!duplicated(u$sec),]
  u[c("X.1","X")] <- NULL
  mn <- min(atmd$sec);mx <- max(atmd$sec)
  u <- u[between(u$sec,mn,mx),]
  if(i == "ps"){
    atmd$psrf <- NULL
    colnames(u) <- c("psrf","sec")
    atmd <- merge(atmd,u,all=T)
  }else if(i == "rf"){
    atmd$precip <- NULL 
    colnames(u) <- c("precip","sec")
    atmd <- merge(atmd,u,all=T)
  }else if(i == "rh"){
    atmd$rh <- NULL 
    colnames(u) <- c("rh","sec")
    atmd <- merge(atmd,u,all=T)
  }else{xxx <- u$V1}
}

x<-atmd
x %>% tally()
x <- x%>% group_by(fname) %>% tally()

atmd$date<- NULL

atmd$day <- as.numeric(substr(atmd$tm,9,10))+(as.numeric(substr(atmd$tm,12,13))/24)-1

#-NaNs Check -------------------------------------------------------------------
#-Atmospheric Pressure
atmd[is.nan(atmd$psrf)] <- NA
atmd$psrf <- ifelse(is.na(atmd$psrf)==TRUE,mean(atmd$psrf,na.rm=T),atmd$psrf)
na_pa <- atmd[is.na(atmd$psrf),] %>%  print()
hist(atmd$psrf)

#-Precipitation
atmd$precip[is.nan(atmd$precip)] <- NA
u <- read.csv(paste(gap_path,"csv_",site_id,"/",site_id,"_gap_rf.csv",sep = ""))
u[c("X.1","X")] <- NULL
mn <- min(atmd$sec);mx <- max(atmd$sec)
u <- u[between(u$sec,mn,mx),]
atmd$precip <- ifelse(is.na(atmd$precip) == T, u$V1[u$sec==atmd$sec],atmd$precip)
atmd$precip <- ifelse(is.na(atmd$precip)==TRUE,0,atmd$precip)
na_precip <- atmd[is.na(atmd$precip),] %>%  print()
hist(atmd$precip)

#-Temperature
atmd$ta[is.nan(atmd$ta)] <- NA
u <- read.csv(paste(gap_path,"csv_",site_id,"/",site_id,"_gap_ta.csv",sep = ""));
u[c("X.1","X")] <- NULL
mn <- min(atmd$sec);mx <- max(atmd$sec)
u <- u[between(u$sec,mn,mx),]
atmd$ta <- ifelse(is.na(atmd$ta)==T,u$V1[u$sec==atmd$sec],atmd$ta)
na_tbot <- atmd[is.na(atmd$ta),] %>%  print()
hist(atmd$ta)

#-Relative Humidity
atmd$rh[is.nan(atmd$rh)] <- NA
u <- read.csv(paste(gap_path,"csv_",site_id,"/",site_id,"_gap_rh.csv",sep = ""));
u[c("X.1","X")] <- NULL
mn <- min(atmd$sec);mx <- max(atmd$sec)
u <- u[between(u$sec,mn,mx),]
atmd$rh <- ifelse(is.na(atmd$rh)==T,u$V1[u$sec==atmd$sec],atmd$rh)
na_rh <- atmd[is.na(atmd$rh),] %>%  print()
hist(atmd$rh)

#-Wind Speed
atmd$ws[is.nan(atmd$ws)] <- NA
u <- read.csv(paste(gap_path,"csv_",site_id,"/",site_id,"_gap_ws.csv",sep = ""));
u[c("X.1","X")] <- NULL
mn <- min(atmd$sec);mx <- max(atmd$sec)
u <- u[between(u$sec,mn,mx),]
atmd$ws <- ifelse(is.na(atmd$ws)==TRUE,mean(atmd$ws,na.rm=T),atmd$ws)
atmd$ws <- ifelse(atmd$ws < 0,0,atmd$ws)
na_ws <- atmd[is.na(atmd$ws),] %>%  print()
hist(atmd$ws)

#-Short Wave Radiation
atmd$sw_in[is.nan(atmd$sw_in)] <- NA
u <- read.csv(paste(gap_path,"csv_",site_id,"/",site_id,"_gap_swr.csv",sep = ""));
u[c("X.1","X")] <- NULL
mn <- min(atmd$sec);mx <- max(atmd$sec)
u <- u[between(u$sec,mn,mx),]
atmd$sw_in <- ifelse(is.na(atmd$sw_in)==T,u$V1[u$sec==atmd$sec],atmd$sw_in)
na_swr <- atmd[is.na(atmd$sw_in),] %>%  print()
hist(atmd$sw_in)
  
atmd$lw_in <- with(atmd,lwin(ta,rh))
na_lw_in <- atmd[is.na(atmd$lw_in),] %>%  print()
hist(atmd$lw_in)

atmd <- atmd[!is.na(atmd$fname),]
atmd <- atmd[!(format(strptime(atmd$tm,format = "%Y-%m-%d %H:%M:%S"),"%m") == "02" & format(strptime(atmd$tm,format = "%Y-%m-%d %H:%M:%S"), "%d") == "29"), , drop = FALSE] # With this linee we remove the leap day of leap year
atmd <- atmd[order(atmd$fname,atmd$day),]

x<-atmd
x %>% tally()
x <- x%>% group_by(fname) %>% tally()

atmd$time <- NULL
write.csv(atmd,file = paste(ncpath,"/",site_id,"atmfor.csv",sep = ""))

#-Creating NetCDF files---------------------------------------------------------
for (i in unique(atmd$fname)) {
  TM <- atmd$day[atmd$fname==i]
  
scalar <- ncdim_def(name = "scalar",units="",vals = 1)
time <- ncdim_def(name = "time", units = paste('days since ',i,'-01 00:00:00',sep=""),vals = TM,unlim = TRUE,calendar = "gregorian") ## Addition after first iteration with the CLM
lat <- ncdim_def(name = "lat", units = "degrees",vals = 1,unlim = FALSE)
lon <- ncdim_def(name = "lon", units = "degrees",vals = 1,unlim = FALSE)

#-Defining Variables
edgew <- ncvar_def(name="EDGEW", units='degrees E', dim=list(scalar), missval=-1, longname='western edge in atmospheric data', prec="float")
edgee <- ncvar_def(name="EDGEE", units='degrees E', dim=list(scalar), missval=-1, longname='eastern edge in atmospheric data', prec="float")
edges <- ncvar_def(name="EDGES", units='degrees N', dim=list(scalar), missval=-1, longname='southern edge in atmospheric data', prec="float")
edgen <- ncvar_def(name="EDGEN", units='degrees N', dim=list(scalar), missval=-1, longname='northern edge in atmospheric data', prec="float")

longxy <- ncvar_def(name="LONGXY", units='degrees E', dim=list(lon,lat), missval=-1, longname='longitude', prec="float")
latixy <- ncvar_def(name="LATIXY", units='degrees N', dim=list(lon,lat), missval=-1, longname='latitude', prec="float")

flds <- ncvar_def(name="FLDS", units='W/m2', dim=list(lon,lat,time), missval=-1, longname='incident longwave (FLDS)', prec="float")
fsds <- ncvar_def(name="FSDS", units='W/m2', dim=list(lon,lat,time), missval=-1, longname='incident shortwave (FSDS)', prec="float")
prectmms <- ncvar_def(name="PRECTmms", units='mm/s', dim=list(lon,lat,time), missval=-1, longname='precipitation (PRECTmms)', prec="float")
psrf <- ncvar_def(name="PSRF", units='Pa', dim=list(lon,lat,time), missval=-1, longname='pressure at the lowest atm level (PSRF)', prec="float")
tbot <- ncvar_def(name="TBOT", units='K', dim=list(lon,lat,time), missval=-1, longname='temperature at the lowest atm level (TBOT)', prec="float")
rh <- ncvar_def(name="RH", units='%', dim=list(lon,lat,time), missval=-1, longname='relative humidity at the lowest atm level (RH)', prec="float")
wind <- ncvar_def(name="WIND", units='m/s', dim=list(lon,lat,time), missval=-1, longname='wind at the lowest atm level (WIND)', prec="float")
zbot <- ncvar_def(name="ZBOT", units='m', dim=list(lon,lat,time), missval=-1, longname='observational height', prec="float")

vars <- list(edgew,edgee,edges,edgen,longxy,latixy,flds,prectmms,psrf,fsds,tbot,rh,wind,zbot)

  ta <- atmd$ta[atmd$fname==i]+273.15
  relh <- atmd$rh[atmd$fname==i]
  ws <- atmd$ws[atmd$fname==i]
  swd <- atmd$sw_in[atmd$fname==i]
  lwd <- atmd$lw_in[atmd$fname==i]
  prec <- atmd$precip[atmd$fname==i]/3600
  a_prsr <- atmd$psrf[atmd$fname==i]
  z <- rep(as.numeric(stand$st_height)+1,length(TM))
  
# Create netCDF file and put arrays
ncfname <- paste(ncpath,"/",i,".nc", sep="")
xx <- nc_create(filename = ncfname,vars = vars,force_v4 = TRUE,verbose = FALSE)
ncvar_put(xx,edgew,floor(site$si_long))
ncvar_put(xx,edgee,ceiling(site$si_long))
ncvar_put(xx,edges,floor(site$si_lat))
ncvar_put(xx,edgen,ceiling(site$si_lat))
ncvar_put(xx,longxy,site$si_long)
ncvar_put(xx,latixy,site$si_lat)
ncvar_put(xx,flds,lwd)
ncvar_put(xx,fsds,swd)
ncvar_put(xx,prectmms,prec)
ncvar_put(xx,psrf,a_prsr)
ncvar_put(xx,tbot,ta)
ncvar_put(xx,rh,relh)
ncvar_put(xx,wind,ws)
ncvar_put(xx,zbot,z)

nc_close(xx)
}

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#
