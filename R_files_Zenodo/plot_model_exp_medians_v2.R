#*******************************************************************************
# Cesar D. Jimenez-Rodriguez
# Luxembourg Institute of Science and Technology (LIST)
# Environmental Research and Innovation Department (ERIN)
# Environmental Sensing and Modelling Unit
# Agro-Environmental Systems Group
# E-mail: cesar.jimenez@list.lu 
#*******************************************************************************

#-Script Description------------------------------------------------------------
# This script produces a set plots per experimental site comparing the medians 
# per day of the year for the three model configurations.

#-Required Packages---------------------------------------------------------------------------------
library(dplyr)
library(readr)  
library(lubridate)
library(tidyverse)
library(data.table)
library(tibble)

#-User settings-----------------------------------------------------------------
rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
res_path  = "" # Directory path where site results are located
par_path = "" # Directory path where the parameter file is located
#-Select ID of the SapFluxNet site --------------------------------------------#
stid <- c("FRA_HES_HE2_NON","FRA_PUE","ESP_ALT_TRI","RUS_FYO")
pftn <- c(8,6,6,2)
yri <- c(2001,2001,2012,2001) # Initial year selected for analysis per site
yrf <- c(2005,2011,2014,2004) # Final year selected for analysis per site (inclusive)
st <- c("FR-Hes","FR-Pue","ES-Alt","RU-Fyo")
p50 <- c(-2.648,-2.648,-2.648,-5.197)

exp <- c("std","brknn","brk91")
exp_ppr <- c("DMC","DBC","FBC")

#-Data Retrieval and Plotting --------------------------------------------------
sites <- data.frame(stid,pftn,yri,yrf,st,p50)
u <- as.numeric(length(sites$stid))

sites$stid <- factor(sites$stid, levels=c("ESP_ALT_TRI","FRA_PUE","FRA_HES_HE2_NON","RUS_FYO"))

pdf(paste(res_path,"/models.pdf",sep=""),width = 7,height = (10))
layout(matrix(c(1:4),4,1, byrow = TRUE), widths=c(1), heights=c(1))
for (i in unique(sites$stid)) {
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

rf <- df[df$id=="std",]; rf <- rf[c("dates","RF","yr","doy")]
rf <- rf[order(rf$dates),]
rf$accu <- NA
for (j in yrk) {
rf$accu[rf$yr==j] <- cumsum(rf$RF[rf$yr==j])   
}

df$mnth <- month(df$dates)
df$pmdif <- df$shade-df$root

sd <- c(0.02,0.04,0.06,0.08,0.12,0.16,0.2)
df$smp <- df$sl1*sd[1]+df$sl2*sd[2]+df$sl3*sd[3]+df$sl4*sd[4]+df$sl5*sd[5]+df$sl6*sd[6]+df$sl7*sd[7]
df$dif_s_r <- df$smp-df$root
df$mnth <- as.numeric(df$mnth)

df$bias <- df$ET-df$stand

y <-df[c("sd_stand","doy")]; y$sd_stand[is.na(y$sd_stand)] <- 0

y <- aggregate(y$sd_stand,list(y$doy),FUN=sd)

ss <- as.character(sites$st[sites$stid==i]);ss

df$smp_n <- df$smp/sites$p50[sites$stid==i]
df$smp_col <- NA

e <- c("sl1","sl2","sl3","sl4","sl5","sl6","sl7")

dt <- df[df$id=="std",]

par(mar=c(4,5,0.5,0.5),xpd=F,font.axis = 1) #(DLUR)

x <- aggregate(dt$stand,list(dt$doy),FUN=median);colnames(x)<-c("doy","stand")
plot(x$doy,x$stand,ylim=c(0,max(df$ET,df$stand)),pch=16,cex=0.5,lwd=1.5,type="l",col="forestgreen",xaxt="n",xlab = "Day of the Year",ylab="\\n")
axis(1,at=c(1,60,120,180,240,300,365))
mtext(expression(italic("E")["T"]*" [mm d"^{-1}*"]"),side = 2,line = 3,lwd=1.5)
legend("topleft",legend=c(as.character(sites$st[sites$stid==i]),paste(sites$yri[sites$stid==i],"-",sites$yrf[sites$stid==i],sep=""),expression(italic(E)["T"]),"DMC","DBC","FBC"),bty="n",
       lty=c(NA,NA,1,1,2,2),lwd=2,col=c(NA,NA,"forestgreen","honeydew4","red","purple"),cex=1)

x <- aggregate(df$ET[df$id=="std"],list(df$doy[df$id=="std"]),FUN=median);colnames(x)<-c("doy","et")
points(x$doy,x$et,col="honeydew4",pch=16,cex=0.5,type = "l",lwd=1.75,lty=1)

x <- aggregate(df$ET[df$id=="brknn"],list(df$doy[df$id=="brknn"]),FUN=median);colnames(x)<-c("doy","et")
points(x$doy,x$et,col="red",pch=16,cex=0.5,type = "l",lwd=1.5,lty=2)

x <- aggregate(df$ET[df$id=="brk91"],list(df$doy[df$id=="brk91"]),FUN=median);colnames(x)<-c("doy","et")
points(x$doy,x$et,col="purple",pch=16,cex=0.5,type = "l",lwd=2.5,lty=3)
}
dev.off()

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#