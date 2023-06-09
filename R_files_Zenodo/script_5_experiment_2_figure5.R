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
# site save as .csv file. It creates figure 5 and uses the following model outputs:
# klr: it has the high kmax configuration
# kalr: it has the intermediate high kmax configuration
# std: it has the default configuration
# kahr: it has the intermediate low kmax configuration
# khr: it has the low kmax configuration
#
# NOTE: paper2_exp_koptimum_v5.R:202108#
#-Required Packages-------------------------------------------------------------
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
yri = c(2001,2012,2001,2012)
yrf = c(2005,2014,2005,2014)
sites = data.frame(stid,yri,yrf)
  
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
xpr <- c("klr","kalr","std","kahr","khr")
  
setwd(res_path)
  
  pdf(paste(res_path,"figure5.pdf",sep=""),width = 12,height = (13))
  layout(matrix(c(1:6,31,7:12,31,13:18,31,19:24,31,25:30,31),7,5, byrow = F), widths=c(0.25,1,1,1,1), heights=c(1,1,1,1,1,0.15,0.20))
  par(mar=c(2,5,0.25,0),xpd=T,font.axis = 1) #(DLUR)
  plot(c(50,330),c(-0.5,6),xlim = c(50,330),type = "n",xaxt="n",ylim = c(-0.5,6),ylab = expression(italic("E")*" [mm d"^{-1}*"]"),bty="n",xlab = "\\n",cex.lab=2,cex.axis=1.5)
  plot(c(50,330),c(-0.5,6),xlim = c(50,330),type = "n",xaxt="n",ylim = c(-0.5,6),ylab = expression(italic("E")*" [mm d"^{-1}*"]"),bty="n",xlab = "\\n",cex.lab=2,cex.axis=1.5)
  plot(c(50,330),c(-0.5,6),xlim = c(50,330),type = "n",xaxt="n",ylim = c(-0.5,6),ylab = expression(italic("E")*" [mm d"^{-1}*"]"),bty="n",xlab = "\\n",cex.lab=2,cex.axis=1.5)
  plot(c(50,330),c(-0.5,6),xlim = c(50,330),type = "n",xaxt="n",ylim = c(-0.5,6),ylab = expression(italic("E")*" [mm d"^{-1}*"]"),bty="n",xlab = "\\n",cex.lab=2,cex.axis=1.5)
  plot(c(50,330),c(-0.5,6),xlim = c(50,330),type = "n",xaxt="n",ylim = c(-0.5,6),ylab = expression(italic("E")*" [mm d"^{-1}*"]"),bty="n",xlab = "\\n",cex.lab=2,cex.axis=1.5)
  plot.new()
  
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
      es <- aggregate(df$ES,list(df$yr,df$doy),FUN=sum);colnames(es)<-c("yr","doy","es")
      et <- aggregate(df$ET,list(df$yr,df$doy),FUN=sum);colnames(et)<-c("yr","doy","et")
      beta <- aggregate(df$beta,list(df$yr,df$doy),FUN=mean);colnames(beta)<-c("yr","doy","beta")
      
      df <- merge(et,beta,by=c("yr","doy"))
      df <- merge(df,es,by=c("yr","doy"))
      df$exp <- i
      
      assign(i,df)
      rm(df,et,beta)
    }
    
    df <- rbind(klr,kalr,std,kahr,khr)
    df <- df[df$yr %in% yrp,]
    
    df$season <-NA
    df$season[df$doy %in% c(60:151)] <- "Spring"
    df$season[df$doy %in% c(152:243)] <- "Summer"
    df$season[df$doy %in% c(243:334)] <- "Autumn"
    df$season[is.na(df$season)] <- "Winter"
    df<-merge(df,sap,by=c("yr","doy"))
    df$exp <- factor(df$exp,levels = c("klr","kalr","std","kahr","khr"))
    df <- df[!(df$season=="Winter"),]
    
    ylm = round(max(df$et,df$sap)+1,0)
    
    color.gradient <- function(x, colors=c("firebrick3","antiquewhite","cornflowerblue"), colsteps=100) {
      return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(0,1, length.out=colsteps)) ] )
    }
  
    bt <- c(0,0.5,1)
    bt<-data.frame(bt);colnames(bt)<-"bt"
    bt$col<- color.gradient(bt$bt)

    for (k in unique(xpr)) {
    nm <- c("H","IH","DC","IL","L")
    ix <- nm[as.numeric(which(xpr == k))]
      
    da <- df[df$exp==k,] 
    sp<- aggregate(da$sap,list(da$doy),FUN=mean);colnames(sp)<-c("doy","sap")
    x<- aggregate(da$et,list(da$doy),FUN=mean);colnames(x)<-c("doy","et")
    y<- aggregate(da$beta,list(da$doy),FUN=mean);colnames(y)<-c("doy","beta")
    y$col <- color.gradient(y$beta);y$y <--10
    z<- aggregate(da$es,list(da$doy),FUN=mean);colnames(z)<-c("doy","es")
    x<-merge(x,y,by="doy")
    
    par(mar=c(2,0,0.25,0.25),xpd=F,font.axis = 1) #(DLUR)
    plot(x$doy,x$et,xlim = c(50,330),type = "h",ylim = c(-0.5,6),yaxt="n",ylab="\\n",xlab="\\n",lwd=2,col=x$col,cex.lab=1.5,cex.axis=1.25)
    points(x$doy,x$et,type = "l",lwd=1.5,col="darkgrey")
    points(y$doy,y$y,col=y$col,type = "h",lwd=2,lend=3)
    points(sp$doy,sp$sap,col="forestgreen",type = "l",lwd=1.5)
    points(z$doy,z$es,col="purple",type = "l",lwd=1.5)
  
    da$mae <- abs(da$sap-da$et)
    mae <- round(mean(da$mae),3)
    r2 <- round(cor(da$sap,da$et,method = "pearson"),4)
    mylabel1 = bquote(italic(R)^2 == .(r2))
    mylabel2 = bquote(.(ix)*italic(k)["max"])
    mylabel3 = bquote(MAE == .(mae))
    legend("topright",legend = c(mylabel2,mylabel3),bty = "n",cex=1.75)
    }
    par(mar=c(0,0,0,0),xpd=F,font.axis = 1) #(DLUR)
    plot.new()
    legend("top",legend = c(stid),bty = "n",cex=2)
  }
  plot.new()
  legend("bottom",legend = c(expression(italic("E")["T mea"]),expression(italic("E")["T mod"]),expression(italic("E")["S"]),expression(beta*"=0.0"),expression(beta*"=0.5"),expression(beta*"=1.0")),col=c("forestgreen","grey","purple",bt$col),pch=c(NA,NA,NA,15,15,15),lwd=c(2,2,2,NA,NA,NA),lty=c(1,1,1,NA,NA,NA),bty = "n",cex=2,horiz = T)
  dev.off()
  
#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#