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
# This script reads the summary results obtained with script 1 for each experimental
# site save as .csv file. It creates figure 6 and uses the following model outputs:
# kolr: it has the best fitted kmax and the vulnerable configuration of p50 and ck
# kostd: it has the best fitted kmax and the default configuration
# kohr: it has the best fitted kmax and the resistant configuraiton of p50 and ck
#
# NOTE: paper2_exp_koptimum_v5.R:202108#
#-Required Packages---------------------------------------------------------------------------------
library(ncdf4)
library(dplyr)
library(lubridate)
library(vioplot)
library(tidyverse)

#-User settings-----------------------------------------------------------------
rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
res_path  = "./results/"  # Folder path where are stored the csv files from script 1
out_path = "./figures/" # Folder path to store the figure files.

#-Select ID of the SapFluxNet site --------------------------------------------#
stid = c("FR-Hes","DE-Hin","FR-Pue","ES-Alt")
kmax = c(2.0e-8,1.1e-8,1.1e-8,6.5e-9)
yri = c(2001,2012,2001,2012)
yrf = c(2005,2014,2005,2014)
sites = data.frame(stid,kmax,yri,yrf)
  
#### DO NOT CHANGE FROM THIS POINT ONWARDS #####################################
# Extra Functions --------------------------------------------------------------
VPD <- function(ta,rh){
    es   <- 0.6108*exp((17.27*ta)/(ta+237.3))
    ea   <- es*rh/100
    vpd  <- es-ea
    return(vpd)
}
  
MPa <- (9.80655/1000000) # Used to convert from MPa into mmH20
  
# Extra Functions --------------------------------------------------------------
  xpr <- c("kolr","kostd","kohr")
  
  setwd(res_path)
  
  pdf(paste(res_path,"figure6.pdf",sep=""),width = 9,height = (8))
  layout(matrix(c(1,2,3,4,9,10,11,12,5,6,7,8,13,14,15,16,17,18,19,20,25,26,27,28,21,22,23,24,29,30,31,32,33,33,33,33,33,33,33,33),5,8, byrow = T), widths=c(0.25,0.5,0.5,0.5,0.25,0.5,0.5,0.5), heights=c(1,1,1,1,0.2))
  
  for (stid in unique(sites$stid)) {
  kmax = sites$kmax[sites$stid==stid]
  yrp = c(sites$yri[sites$stid==stid]:sites$yrf[sites$stid==stid])
  j=stid
  sap <- read.csv(paste(res_path,"sap_",stid,".csv",sep=""))
  sap$X <-NULL
  colnames(sap) <-c("time","sap","doy","yr")
  
  for (i in xpr) {
  df <- read.csv(paste(res_path,stid,"/",stid,"_",i,"_hr.csv",sep="")) # Read the experiment data
  df$beta <- ((df$BSHA*df$SHALAI)+(df$BSUN*df$SUNLAI))/(df$SHALAI+df$SUNLAI)
  
  et <- aggregate(df$ET,list(df$yr,df$doy),FUN=sum);colnames(et)<-c("yr","doy","et")
  beta <- aggregate(df$beta,list(df$yr,df$doy),FUN=mean);colnames(beta)<-c("yr","doy","beta")
  plc_x <- aggregate(df$kxyl,list(df$yr,df$doy),FUN=mean);colnames(plc_x)<-c("yr","doy","plc_x")
  plc_x$plc_x <- 100*(1-plc_x$plc_x/kmax)
  plc_l <- aggregate(df$kleaf,list(df$yr,df$doy),FUN=mean);colnames(plc_l)<-c("yr","doy","plc_l")
  plc_l$plc_l <- 100*(1-(plc_l$plc_l/kmax))
  
  df <- merge(et,beta,by=c("yr","doy"))
  df <- merge(df,plc_x,by=c("yr","doy"))
  df <- merge(df,plc_l,by=c("yr","doy"))
  df$exp <- i
  
  assign(i,df)
  rm(df,et,beta,plc_x,plc_l)
  }
  
  df <- rbind(kolr,kostd,kohr)
  df <- df[df$yr %in% yrp,]
  
  df$season <-NA
  df$season[df$doy %in% c(60:151)] <- "Spring"
  df$season[df$doy %in% c(152:243)] <- "Summer"
  df$season[df$doy %in% c(243:334)] <- "Autumn"
  df$season[is.na(df$season)] <- "Winter"
  df<-merge(df,sap,by=c("yr","doy"))
  df$exp <- factor(df$exp,levels = c("kolr","kostd","kohr"))
  df <- df[!(df$season=="Winter"),]
  
ylm = ceiling(max(df$et,df$sap)+2)
  
par(mar=c(0.25,4.5,2,0.25),xpd=F,font.axis = 1) #(DLUR)
plot(df$et,ylim=c(0,ylm),type="n",xaxt="n",bty="n",xlab = "\\n",ylab=expression(italic("E")["T"]*" [mm d"^{-1}*"]"),cex.lab=1.5,cex.axis=1.15)
par(mar=c(0.25,0,2,0.25),xpd=F,font.axis = 1) #(DLUR)  
  for (i in unique(df$season)) {
    vioplot(df$sap[df$season == i]~df$exp[df$season == i],col = "forestgreen", plotCentre = "line", side = "right",ylim=c(0,ylm),axes = T,yaxt="n")
    par(new=T)
    vioplot(df$et[df$season == i]~df$exp[df$season == i],col=c("red","grey","blue"), plotCentre = "line", side = "left",ylim=c(0,ylm),axes = T,yaxt="n")
    mtext(i,side = 3,line = 0.5)
    #legend("topleft",legend = c(j,expression(italic("E")["T-measured"])),cex=1.25,bty="n",pch=c(NA,15),col=c(NA,"forestgreen"))
  }

par(mar=c(1,4.5,2,0.25),xpd=F,font.axis = 1) #(DLUR)
plot(df$plc_l,ylim=c(100,0),type="n",xaxt="n",bty="n",xlab = "\\n",ylab=expression("PLC [%]"),cex.lab=1.5,cex.axis=1.15)
par(mar=c(1,0,2,0.25),xpd=F,font.axis = 1) #(DLUR)  
  for (i in unique(df$season)) {
    vioplot(df$plc_l[df$season == i]~df$exp[df$season == i],col = "lawngreen", plotCentre = "line", side = "right",ylim=c(110,0),axes = T,yaxt="n")
    par(new=T)
    vioplot(df$plc_x[df$season == i]~df$exp[df$season == i],col="peru", plotCentre = "line", side = "left",ylim=c(110,0),axes = T,yaxt="n")
    axis(3,at=c(1,2,3),labels = c("VC","DC","RC"))
    #legend("bottomright",legend = c("Root-Stem","Stem-Leaf"),pch = 15,col = c("peru","lawngreen"),bty="n",cex=1.25,horiz = T)
    if(i == "Spring"){
      legend("bottomleft",legend = c(stid),bty="n",cex=2,horiz = T)
    }else{next}
  }
  }
  par(mar=c(0,0,0,0),xpd=F,font.axis = 1) #(DLUR)  
  plot.new()
  legend("bottom",legend = c(expression(italic("E")["T"]),expression(italic("E")["T-VC"]),expression(italic("E")["T-DC"]),expression(italic("E")["T-RC"]),"R-S","S-L"),pch = 15,col = c("forestgreen","red","grey","blue","peru","lawngreen"),bty="n",cex=1.75,horiz = T)
  dev.off() 

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#