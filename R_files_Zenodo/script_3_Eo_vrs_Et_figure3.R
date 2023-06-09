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
# Non-Parametric comparison among models

#-Required Packages-------------------------------------------------------------
library(dplyr)
library(readr)  
library(lubridate)
library(tidyverse)
library(data.table)
library(tibble) # Required for a further step
library(vioplot)

#-User settings-----------------------------------------------------------------
rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
res_path  = "./results/"  # Folder path where are stored the csv files from script 1 and 2

#-Select ID of the SapFluxNet site --------------------------------------------#
stid <- c("FR-Hes","DE-Hin","FR-Pue","ES-Alt")
pftn <- c(8,8,6,6)
yri <- c(2001,2012,2001,2012) # Initial year selected for analysis per site
yrf <- c(2005,2014,2005,2014) # Final year selected for analysis per site (inclusive)
sps <- c("Fagus sylvatica","Fagus sylvatica","Quercus ilex","Quercus ilex")

#-Extra Functions --------------------------------------------------------------
  
#-Data Retrieval----------------------------------------------------------------
sites <- data.frame(stid,pftn,yri,yrf,sps);rm(stid,yri,yrf,pftn,sps)
u <- as.numeric(length(sites$stid))
dfa <- NA

pdf(paste(res_path,"/figure3.pdf",sep=""),width = 10,height = (8))
layout(matrix(c(1:16),4,4, byrow = FALSE), widths=c(1), heights=c(1,1,1,0.25))

for (i in sites$stid) {
  wd <- paste(res_path,"/",i,sep="")
  setwd(wd)

sap <- read.csv(paste("./","sap_",i,".csv",sep = ""))
sap <- sap[!(format(strptime(sap$time,format = "%Y-%m-%d"),"%m") == "02" & format(strptime(sap$time,format = "%Y-%m-%d"), "%d") == "29"), , drop = FALSE] # With this one we remove the leap day of leap year
sap$site <- i
sap$date <- sap$time
sap$doy <- as.numeric(strftime(strptime(substr(sap$date,6,10),format = "%m-%d"), format = "%j"))
sap$time <- paste(substr(sap$date,0,4),"-",sap$doy,sep="");sap[c("X","doy")]<-NULL

s <- list.files(pattern="*.csv$")
s <- s[lapply(s,function(x) length(grep("sap",x,value=FALSE))) == 0]

df_input_list <- lapply(s, fread)
names(df_input_list) <- gsub(s, pattern=paste(i,"_",sep=""), replacement="")

df <- bind_rows(df_input_list, .id = "id")
df$id <- gsub(df$id, pattern=".csv", replacement="")
df$doy <- df$doy+1
df$time <- paste(df$yr,"-",df$doy,sep="")
df <- merge(sap,df,by="time",all=T)
df$yr <- as.integer(substr(df$time,0,4))

yrk <- c(sites$yri[sites$stid==i]:sites$yrf[sites$stid==i])

df <- df[df$yr %in% yrk,]
df$stand[is.na(df$stand)] <- 0
df$dates <- strptime(df$date,format = "%Y-%m-%d")
df$mnth <- month(df$dates)
df$pmdif <- df$shade-df$root #Plant matirc differential from shaded/roots

sd <- c(0.02,0.04,0.06,0.08,0.12,0.16,0.2)
df$smp <- df$sl1*sd[1]+df$sl2*sd[2]+df$sl3*sd[3]+df$sl4*sd[4]+df$sl5*sd[5]+df$sl6*sd[6]+df$sl7*sd[7]
df$dif_s_r <- df$smp-df$root
df$mnth <- as.numeric(df$mnth)

df <- df[order(df$dates),]

dt <- df[df$id=="std",]
dt <- dt[order(dt$dates),]

df$season <-NA
for (j in unique(df$doy)) {
  if(df$doy[df$doy==j] %in% c(60:151)){
    df$season[df$doy==j] <- "Spring"
  }else if(df$doy[df$doy==j] %in% c(152:243)){
    df$season[df$doy==j] <- "Summer"
  }else if(df$doy[df$doy==j] %in% c(243:334)){
    df$season[df$doy==j] <- "Autumn"
  }else{
    df$season[df$doy==j] <- "Winter" 
  }
}

df <- df[!df$season=="Winter",]

df$season <- factor(df$season,levels = c("Spring","Summer","Autumn"))

spring = rgb(0.2,0.8,0.2,0.25); spring_l=rgb(0.2,0.8,0.2)
summer = rgb(1,0.6,0.2,0.25); summer_l = rgb(1,0.6,0.2)
fall = rgb(0.4,0.4,0,0.5); fall_l = rgb(0.4,0.4,0)
k= NA

for (j in unique(df$season)) {
  if(j == "Spring"){
    k = spring_l
  }else if(j == "Summer"){
    k = summer_l
  }else{k = fall_l}
    
  par(mar=c(0.25,4.5,0.25,0.25),xpd=F,font.axis = 1) #(DLUR)
  vioplot(df$stand[df$season==j],col = rgb(0.14,0.46,0.14,0.5), plotCentre = "line", side = "left",ylim=c(0,10),ylab = "\\n",axes = F,names = c(""))
  xm <- round(median(df$stand[df$season==j]),2)
  text(0.525,xm,labels = xm,font=2)
  par(new=T)
  vioplot(df$Eref[df$season==j],col = "grey", plotCentre = "line", side = "right",ylim=c(0,10),ylab = expression("mm d"^{-1}),axes=F,names = c(""))
  xm <- round(median(df$Eref[df$season==j]),2)
  text(1.475,xm,labels = xm,font = 2)
  legend("topleft",legend = c(j,expression(italic(E)["T"]),expression(italic(E)["o"])),pch=c(NA,15,15),col = c(NA,"forestgreen","grey"),bty="n",cex=1.5)#,text.col =c(k,"black","black"))
}
par(mar=c(0.25,0.25,0.25,0.25),xpd=F,font.axis = 1) #(DLUR)
plot.new()
legend("top",legend=i,bty="n",cex=1.5,horiz = F)
legend("bottom",legend=paste("[",sites$yri[sites$stid==i],"-",sites$yrf[sites$stid==i],"]",sep=""),bty="n",cex=1.25,horiz = F)
}
dev.off()

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#