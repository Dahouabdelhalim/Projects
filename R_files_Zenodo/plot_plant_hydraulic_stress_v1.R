#*******************************************************************************
# Cesar D. Jimenez-Rodriguez
# Luxembourg Institute of Science and Technology (LIST)
# Environmental Research and Innovation Department (ERIN)
# Environmental Sensing and Modelling Unit
# Agro-Environmental Systems Group
# E-mail: cesar.jimenez@list.lu 
#*******************************************************************************

#-Script Description--------------------------------------------------------------------------------
# This script produces a set of violin plots comparing the stress in plant water
# transport for the months June, July, August, September, and October. The 
# comparison is based on the xylem-stem and xylem-leaf.

#-Required Packages---------------------------------------------------------------------------------
library(vioplot)

#-User settings-----------------------------------------------------------------
rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
site_id = c('ESP_ALT_TRI','FRA_PUE','FRA_HES_HE2_NON','RUS_FYO') # Site code from SapFluxNet.
stid <- c('ES-Alt','FR-Pue','FR-Hes','RU-Fyo')    # Site code from Fluxnet.
pftn= c(6,8,6,8)              # Plant functional type from the CLM list.
yri = c(2012,2001,2001,2001)
yrf = c(2014,2011,2005,2005)

mp <- data.frame(site_id,stid,pftn,yri,yrf)
mp$site_id <- as.character(mp$site_id)
mp$stid <- as.character(mp$stid)

sns <- data.frame(c(0:364));colnames(sns)<-"doy"
for (j in unique(sns$doy)) {
  if(sns$doy[sns$doy==j] %in% c(60:151)){
    sns$season[sns$doy==j] <- "Spring"
  }else if(sns$doy[sns$doy==j] %in% c(152:243)){
    sns$season[sns$doy==j] <- "Summer"
  }else if(sns$doy[sns$doy==j] %in% c(243:334)){
    sns$season[sns$doy==j] <- "Autumn"
  }else{
    sns$season[sns$doy==j] <- "Winter" 
  }
}
sns$season <- factor(sns$season,levels = c("Winter","Spring","Summer","Autumn"))

for (j in unique(sns$doy)) {
  if(sns$doy[sns$doy==j] %in% c(0:31)){
    sns$mnt[sns$doy==j] <- 1
  }else if(sns$doy[sns$doy==j] %in% c(32:59)){
    sns$mnt[sns$doy==j] <- 2
  }else if(sns$doy[sns$doy==j] %in% c(60:90)){
    sns$mnt[sns$doy==j] <- 3
  }else if(sns$doy[sns$doy==j] %in% c(91:120)){
    sns$mnt[sns$doy==j] <- 4
  }else if(sns$doy[sns$doy==j] %in% c(121:151)){
    sns$mnt[sns$doy==j] <- 5
  }else if(sns$doy[sns$doy==j] %in% c(152:181)){
    sns$mnt[sns$doy==j] <- 6
  }else if(sns$doy[sns$doy==j] %in% c(182:212)){
    sns$mnt[sns$doy==j] <- 7
  }else if(sns$doy[sns$doy==j] %in% c(212:243)){
    sns$mnt[sns$doy==j] <- 8
  }else if(sns$doy[sns$doy==j] %in% c(244:273)){
    sns$mnt[sns$doy==j] <- 9
  }else if(sns$doy[sns$doy==j] %in% c(274:304)){
    sns$mnt[sns$doy==j] <- 10
  }else if(sns$doy[sns$doy==j] %in% c(305:334)){
    sns$mnt[sns$doy==j] <- 11
  }else{
    sns$mnt[sns$doy==j] <- 12 
  }
}

pdf("~/Documents/CAPACITY/analysis/1_paper/PLC_3config.pdf",width = 12,height = (8))
layout(matrix(c(1:12),3,4, byrow = FALSE), widths=c(1), heights=c(1))
par(mar=c(2.5,4.5,0.35,0.25),xpd=F,font.axis = 1) #(DLUR)

for (i in mp$site_id) {

 x <- mp[mp$site_id==i,]
std <- read.csv(paste("~/",i,"/",i,"_std.csv",sep=""))
fbc <- read.csv(paste("~/",i,"/",i,"_brk91.csv",sep=""))
dbc <- read.csv(paste("~/",i,"/",i,"_brknn.csv",sep=""))

std <- std[c("yr","doy","kxyl","kleaf","ET")];std<-merge(std,sns,bby="doy");std<-std[std$yr %in% c(x$yri:x$yrf),];std<-std[std$mnt %in% c(6:10),]
fbc <- fbc[c("yr","doy","kxyl","kleaf","ET")];fbc<-merge(fbc,sns,bby="doy");fbc<-fbc[fbc$yr %in% c(x$yri:x$yrf),];fbc<-fbc[fbc$mnt %in% c(6:10),]
dbc <- dbc[c("yr","doy","kxyl","kleaf","ET")];dbc<-merge(dbc,sns,bby="doy");dbc<-dbc[dbc$yr %in% c(x$yri:x$yrf),];dbc<-dbc[dbc$mnt %in% c(6:10),]

vioplot(100*((std$kxyl/(2e-8)))~std$mnt, plotCentre = "line", side = "left",axes = F,col="lemonchiffon3",ylim=c(-10,100),ylab=expression(Xi*" [%]"),names = c("Jun","Jul","Aug","Sep","Oct"))
par(new=T)
vioplot(100*((std$kleaf/(2e-8)))~std$mnt, plotCentre = "line", side = "right",axes = F,col="greenyellow",ylim=c(-10,100),ylab=expression(Xi*" [%]"),names = c("Jun","Jul","Aug","Sep","Oct"))
legend("bottomleft",legend = c(x$stid,"DMC"),col=c(NA,NA),pch=15,bty="n",cex=1)
legend("bottomright",legend = c("Stem","Leaf"),col=c("lemonchiffon3","yellowgreen"),pch=15,bty="n",cex=1)

vioplot(100*((dbc$kxyl/(2e-8)))~std$mnt, plotCentre = "line", side = "left",axes = F,col="lemonchiffon3",ylim=c(-10,100),ylab=expression(Xi*" [%]"),names = c("Jun","Jul","Aug","Sep","Oct"))
par(new=T)
vioplot(100*((dbc$kleaf/(2e-8)))~std$mnt, plotCentre = "line", side = "right",axes = F,col="greenyellow",ylim=c(-10,100),ylab=expression(Xi*" [%]"),names = c("Jun","Jul","Aug","Sep","Oct"))
legend("bottomleft",legend = c(x$stid,"DBC"),col=c(NA,NA),pch=15,bty="n",cex=1)
legend("bottomright",legend = c("Stem","Leaf"),col=c("lemonchiffon3","yellowgreen"),pch=15,bty="n",cex=1)

vioplot(100*((fbc$kxyl/(2e-8)))~fbc$mnt, plotCentre = "line", side = "left",axes = F,col="lemonchiffon3",ylim=c(-10,100),ylab=expression(Xi*" [%]"),names = c("Jun","Jul","Aug","Sep","Oct"))
par(new=T)
vioplot(100*((fbc$kleaf/(2e-8)))~fbc$mnt, plotCentre = "line", side = "right",axes = F,col="greenyellow",ylim=c(-10,100),ylab=expression(Xi*" [%]"),names = c("Jun","Jul","Aug","Sep","Oct"))
legend("bottomleft",legend = c(x$stid,"FBC"),col=c(NA,NA),pch=15,bty="n",cex=1)
legend("bottomright",legend = c("Stem","Leaf"),col=c("lemonchiffon3","yellowgreen"),pch=15,bty="n",cex=1)
}
dev.off()

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#