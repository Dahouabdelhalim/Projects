#*******************************************************************************
# Cesar D. Jimenez-Rodriguez
# Luxembourg Institute of Science and Technology (LIST)
# Environmental Research and Innovation Department (ERIN)
# Environmental Sensing and Modelling Unit
# Agro-Environmental Systems Group
# E-mail: cesar.jimenez@list.lu 
#*******************************************************************************

#-Script Description--------------------------------------------------------------------------------
# This script reads the results file from CLM 5 as ".nc" file and summarizes the 
# main data in daily averages.
# The files are saved as SITE_NAME_EXPERIMENT.csv
#
# NOTE: Post_Proce_v14:202108#

#-Required Packages---------------------------------------------------------------------------------
library(ncdf4)
library(dplyr)
library(lubridate)
#library(tidyverse)

#-User settings-----------------------------------------------------------------

rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
res_path  = "/model_results/"  # Folder path where CLM 5 results are store in ".nc" format.
par_path = '/parameters'               # Folder path where the experiment parameters are store in ".nc" format.
out_path =''          # Folder path to store the output files.
expr = "std"                                             # Experiment ID use to identify the model results.
srfc ='std'                                                # Name used for the surface data set to select from the parameter file.

#-Select ID of the SapFluxNet site --------------------------------------------#
site_id = 'RUS_FYO'; stid <- 'RU-Fyo' ;pftn=2

#### DO NOT CHANGE FROM THIS POINT ONWARDS #####################################
# Extra Functions --------------------------------------------------------------
VPD <- function(ta,rh){
  es   <- 0.6108*exp((17.27*ta)/(ta+237.3))
  ea   <- es*rh/100
  vpd  <- es-ea
  return(vpd)
}

MPa <- (9.80655/1000000) # Used to convert from MPa into mmH20

# Parameters Setup -------------------------------------------------------------
def_exp = c("std","brk91","brknn") # Experiments using only the default parameter file of CLM 5.0

if (expr %in% def_exp) {
  nc <- nc_open(paste(par_path,"clm5_params.c171117.nc",sep="/"))
  ck <- ncvar_get(nc,"ck");ck<-ck[pftn]
  kmax <- ncvar_get(nc,"kmax"); kmax <- kmax[pftn,]
  p50 <- ncvar_get(nc,"psi50");p50 <- p50[pftn,]
  rt_cond <- ncvar_get(nc,"krmax");rt_cond<-rt_cond[pftn]
  p50_sun <- p50[1]*MPa#mm H20
  p50_shade <- p50[2]*MPa#mm H20
  p50_xyl <- p50[3]*MPa#mm H20
  p50_root <- p50[4]*MPa#mm H20
}else{
  nc <- nc_open(paste(par_path,"/",stid,"/","clm5_params.c171117_",stid,"_",expr,".nc",sep=""))
  ck <- ncvar_get(nc,"ck");ck<-ck[pftn]
  kmax <- ncvar_get(nc,"kmax"); kmax <- kmax[pftn,]
  p50 <- ncvar_get(nc,"psi50");p50 <- p50[pftn,]
  rt_cond <- ncvar_get(nc,"krmax");rt_cond<-rt_cond[pftn]
  p50_sun <- p50[1]*MPa#mm H20
  p50_shade <- p50[2]*MPa#mm H20
  p50_xyl <- p50[3]*MPa#mm H20
  p50_root <- p50[4]*MPa#mm H20
}

# MODEL RESULTS ----------------------------------------------------------------
mdl_rslts <- ifelse(expr %in% def_exp,paste(expr,sep = ""),paste(srfc,expr,sep = "_"))

lst <- list.files(paste(res_path,"/",site_id,"/",expr,sep=""))
lst <- lst[-length(lst)]
setwd(paste(res_path,"/",site_id,"/",expr,sep=""))

TIME <- NA; ET <- NA; EI <- NA; ES <- NA
YEAR <- NA
RF <- NA
TA <- NA
RH <- NA
SWR <- NA
LWR <- NA
WND <- NA
VEGWP <- matrix(NA,nrow=1,ncol=4)
SMP <- matrix(NA,nrow=1,ncol=25)
VSWC <- matrix(NA,nrow=1,ncol=20)

SHG <- NA
SHV <- NA
LAI <- NA
TSOI <- matrix(NA,nrow=1,ncol=25)

GSSHA <- NA
GSSUN <- NA

FPSN  <- NA
BT <- NA
FSDS <- NA
FIRA <- NA
GHF <- NA #Ground heat flux
P <- NA

ZWT <- NA

ZBOT <- NA

QOVER <- NA
QDRAI <- NA
QCHARGE <- NA

for (i in unique(lst)) {
  nc <- nc_open(i)
  qover <- ncvar_get(nc,varid = "QOVER",raw_datavals = F)
  qdrai <- ncvar_get(nc,varid = "QDRAI",raw_datavals = F)
  qcharge <- ncvar_get(nc,varid = "QCHARGE",raw_datavals = F)
  et <- ncvar_get(nc,varid = "QVEGT",raw_datavals = F)
  ei <- ncvar_get(nc,varid = "QVEGE",raw_datavals = F)
  es <- ncvar_get(nc,varid = "QSOIL",raw_datavals = F)
  rf <- ncvar_get(nc,varid = "RAIN",raw_datavals = F)
  ta <- ncvar_get(nc,varid = "TBOT",raw_datavals = F)
  rh <- ncvar_get(nc,varid = "RH2M",raw_datavals = F)
  swr <- ncvar_get(nc,varid = "FSDS",raw_datavals = F)
  lwr <- ncvar_get(nc,varid = "FLDS",raw_datavals = F)
  wnd <- ncvar_get(nc,varid = "WIND",raw_datavals = F)
  vegwp <- t(ncvar_get(nc,varid = "VEGWP",raw_datavals = F))
  smp <- t(ncvar_get(nc,varid = "SMP",raw_datavals = F))
  time <- nc[["dim"]][["time"]][["vals"]]
  yr <- rep(as.numeric(substr(i,16,19)),length(time))
  
  shg <- ncvar_get(nc,varid = "FSH_G",raw_datavals = F)
  shv <- ncvar_get(nc,varid = "FSH_V",raw_datavals = F)
  lai <- ncvar_get(nc,varid = "TLAI",raw_datavals = F)
  tsoi <- t(ncvar_get(nc,varid = "TSOI",raw_datavals = F))
  vswc <- t(ncvar_get(nc,varid = "H2OSOI",raw_datavals = F))
  
  gssha <- ncvar_get(nc,varid = "GSSHA",raw_datavals = F)
  gssun <- ncvar_get(nc,varid = "GSSUN",raw_datavals = F)
  
  fpsn <- ncvar_get(nc,varid="FPSN",raw_datavals=F)
  
  bt <- ncvar_get(nc, varid = "BTRANMN", raw_datavals = F)
  fsds <- ncvar_get(nc, varid = "FSDS", raw_datavals = F)
  fira <- ncvar_get(nc, varid = "FIRA", raw_datavals = F)
  ghf <- ncvar_get(nc, varid = "FGR12", raw_datavals = F)
  zwt <- ncvar_get(nc, varid = "ZWT", raw_datavals = F)
  zbot <- ncvar_get(nc,varid = "ZBOT", raw_datavals = F)
  p <- ncvar_get(nc, varid = "PBOT", raw_datavals = F)
  
  TIME <- c(TIME,time)
  ET <- c(ET,et)
  EI <- c(EI,ei)
  ES <- c(ES,es)
  RF <- c(RF,rf)
  TA <- c(TA,ta)
  RH <- c(RH,rh)
  SWR <- c(SWR,swr)
  LWR <- c(LWR,lwr)
  WND <- c(WND,wnd)
  YEAR <- c(YEAR,yr)
  VEGWP <- rbind(VEGWP,vegwp)
  SMP <- rbind(SMP,smp)
  VSWC <- rbind(VSWC,vswc)

  QOVER <- c(QOVER,qover)
  QDRAI <- c(QDRAI,qdrai)
  QCHARGE <- c(QCHARGE,qcharge)
  
  BT <- c(BT,bt)
  FSDS <- c(FSDS,fsds)
  FIRA <- c(FIRA,fira)
  GHF <- c(GHF,ghf)
  P <- c(P,p)
  
  SHG <- c(SHG,shg)
  SHV <- c(SHV,shv)
  LAI <- c(LAI,lai)
  TSOI <- rbind(TSOI,tsoi)
  
  GSSHA <- c(GSSHA,gssha)
  GSSUN <- c(GSSUN,gssun)
  FPSN <- c(FPSN,fpsn)
  
  ZWT <- c(ZWT,zwt)
  
  ZBOT <- c(ZBOT,zbot)
}

ZBOT <- mean(na.omit(ZBOT))

TIME <- na.omit(TIME)
YEAR <- na.omit(YEAR)

QOVER <- na.omit(ET)*3600
QDRAI <- na.omit(QDRAI)*3600
QCHARGE <- na.omit(QCHARGE)*3600

ET <- na.omit(ET)*3600
ES <- na.omit(ES)*3600
EI <- na.omit(EI)*3600
TA <- na.omit(TA)
RF <- na.omit(RF)*3600
RH <- na.omit(RH)
SWR <- na.omit(SWR)
LWR <- na.omit(LWR)
WND <- na.omit(WND)

BT <- BT[-1]
FSDS <- FSDS[-1]
FIRA <- FIRA[-1]
P <- P[-1]
GHF <- GHF[-1]

VEGWP <- na.omit(VEGWP)
VEGWP <- VEGWP*MPa
colnames(VEGWP) <- c("sun","shade","xyl","root")
SMP <- na.omit(SMP)
SMP <- SMP*MPa
smp <- as.data.frame(SMP[,c(1:7)]) # The first 7 layers represent 68 cm depth
colnames(smp) <- c("sl1","sl2","sl3","sl4","sl5","sl6","sl7")

FPSN <- na.omit(FPSN)

ZWT <- na.omit(ZWT)

GSSHA <- GSSHA[-1]
GSSUN <- GSSUN[-1]

SHG <- na.omit(SHG)
SHV <- na.omit(SHV)
LAI <- na.omit(LAI)
TSOI <- na.omit(TSOI)
VSWC <- na.omit(VSWC); VSWC <- data.frame(VSWC)

mdr <- data.frame(TIME,YEAR,ET,EI,ES,TA,RH,RF,SWR,LWR,WND,VEGWP,smp,SHG,SHV,LAI,VSWC,FPSN,ZWT,GSSUN,GSSHA,QDRAI,QOVER,QCHARGE,FSDS,FIRA,P,GHF)

mdr$doy <- floor(mdr$TIME)

colnames(mdr)[1:2] <- c("time","yr")

mdr$TA <- mdr$TA-273.15
mdr$vpd <- ifelse(VPD(mdr$TA,mdr$RH) < 0,0,VPD(mdr$TA,mdr$RH))

dval <- mdr[c("doy","yr","RF","ET","EI","ES","ZWT","root","xyl","sun","shade","GSSHA","GSSUN","vpd","QDRAI","QOVER","QCHARGE","FSDS","FIRA","WND","TA","RH","P","GHF","sl1","sl2","sl3","sl4","sl5","sl6","sl7")]

dval$RF <- ifelse(dval$RF < 0,0,dval$RF)
dval$ET <- ifelse(dval$ET < 0,0,dval$ET)
dval$EI <- ifelse(dval$EI < 0,0,dval$EI)
dval$ES <- ifelse(dval$ES < 0,0,dval$ES)
dval$QCHARGE <- ifelse(dval$QCHARGE < 0,0,dval$QCHARGE)
dval$QDRAI <- ifelse(dval$QDRAI < 0,0,dval$QDRAI)
dval$QOVER <- ifelse(dval$QOVER < 0,0,dval$QOVER)

dval$ZWT <- dval$ZWT/24
dval$root <- dval$root/24
dval$xyl <- dval$xyl/24
dval$sun <- dval$sun/24
dval$shade <- dval$shade/24
dval$GSSHA <- dval$GSSHA/24
dval$GSSUN <- dval$GSSUN/24
dval$vpd <- dval$vpd/24
dval$yd <- paste(dval$yr,dval$doy,sep = "-")
dval$TD <- NA
for (i in unique(dval$yd)) {
  dval$TD[dval$yd==i] <- max(dval$TA[dval$yd==i])-min(dval$TA[dval$yd==i])
};dval$yd <-NULL
dval$TA <- dval$TA/24
dval$TD <- dval$TD/24
dval$FSDS <- dval$FSDS/24
dval$FIRA <- dval$FIRA/24
dval$GHF <- dval$GHF/24
dval$WND <- dval$WND/24
dval$RH <- dval$RH/24
dval$P <- dval$P/24

dval[,c("sl1","sl2","sl3","sl4","sl5","sl6","sl7")] <- dval[,c("sl1","sl2","sl3","sl4","sl5","sl6","sl7")]/24

dval$doy <- paste(dval$yr,dval$doy,sep="-")
dval <- dval %>% group_by(dval$doy) %>% summarise_at(vars(-doy), sum)
colnames(dval)[1] <- "doy"
dval$yr <- as.numeric(substr(dval$doy,1,4))
dval$doy <- as.numeric(sub('.*-', '', dval$doy))

#-Daily Reference Evapotranspiration--------------------------------------------
dval$D <- (4098*(0.6108*exp((17.27*dval$TA)/(dval$TA+237.3))))/((dval$TA+237.3)^2) # Slope of vapour pressure curve
dval$g <- (0.001013*(dval$P/1000))/(2.45*0.622) #psycometric constant
dval$NR <- 0.0864*((dval$FSDS*(1-0.23))-dval$FIRA)

dval$Eref <- (1/(dval$D+dval$g*(1+0.34*dval$WND)))*((0.408*dval$D*dval$NR)+(dval$g*dval$WND*dval$vpd*(900/(dval$TA+273))))
dval$Eref <- ifelse(dval$Eref<=0,0,dval$Eref)
dval$kleaf <- (kmax[1])*(2^(-((dval$xyl/p50_sun)^ck)))
dval$kxyl <- (kmax[3])*(2^(-((dval$root/p50_xyl)^ck)))

dval[c("D","g","NR","FSDS","FIRA","TD")]<-NULL

df <- dval

df$root_sd <- df$root/p50_root 
df$xyl_sd <- df$xyl/p50_xyl
df$sun_sd <- df$sun/p50_sun

# Path to folder ---------------------------------------------------------------
dir.create(paste(out_path,"/",site_id,sep=""))
write.csv(df,file = paste(out_path,"/",site_id,"/",site_id,"_",expr,".csv",sep=""))

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#
