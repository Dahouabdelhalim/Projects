#*******************************************************************************
# Cesar D. Jimenez-Rodriguez
# Luxembourg Institute of Science and Technology (LIST)
# Environmental Research and Innovation Department (ERIN)
# Environmental Sensing and Modelling Unit
# Agro-Environmental Systems Group
# E-mail: cesar.jimenez@list.lu 
#*******************************************************************************

#-Script Description------------------------------------------------------------
# This script plot the plant vulnerability curve per site, experiment, and 
# season.

#-Required Packages-------------------------------------------------------------
library(dplyr)
library(readr)  
library(lubridate)
library(tidyverse)
library(data.table)
library(tibble) # Required for a further step

#-User settings-----------------------------------------------------------------
rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
res_path  = ""
par_path = ''
#-Select ID of the SapFluxNet site --------------------------------------------#
stid <- c("ESP_ALT_TRI","FRA_PUE","FRA_HES_HE2_NON","RUS_FYO")
pftn <- c(6,6,8,2)
yri <- c(2012,2001,2001,2001) # Initial year selected for analysis per site
yrf <- c(2014,2011,2005,2004) # Final year selected for analysis per site (inclusive)
st <- c("ES-Alt","FR-Pue","FR-Hes","RU-Fyo")
p50 <- c(-2.648,-2.648,-2.648,-5.197)

exp <- c("std","brknn","brk91")
exp_ppr <- c("DMC","DBC","FBC")

#-Data Retrieval----------------------------------------------------------------
sites <- data.frame(stid,pftn,yri,yrf,st,p50)#;rm(stid,yri,yrf,pftn,st)
u <- as.numeric(length(sites$stid))

sites$stid <- factor(sites$stid, levels=c("ESP_ALT_TRI","FRA_PUE","FRA_HES_HE2_NON","RUS_FYO"))

for (i in unique(sites$stid)) {
  wd <- paste(res_path,"/",i,sep="")
  setwd(wd)

pdf(paste(res_path,"/",i,"_models_vulcurv.pdf",sep=""),width = 15,height = (10))
layout(matrix(c(1:36),4,9, byrow = T), widths=c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.2), heights=c(1,1,1,0.3))
  
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
df$pmdif <- df$shade-df$root #Plant matirc differential from shaded/roots

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

df$season <-NA
for (j in unique(df$mnth)) {
  if(df$mnth[df$mnth==j] %in% c(3,4,5)){
    df$season[df$mnth==j] <- "Spring"
  }else if(df$mnth[df$mnth==j] %in% c(6,7,8)){
    df$season[df$mnth==j] <- "Summer"
  }else if(df$mnth[df$mnth==j] %in% c(9,10,11)){
    df$season[df$mnth==j] <- "Autumn"
  }else{
    df$season[df$mnth==j] <- "Winter" 
  }
}

EXP <- data.frame(exp,exp_ppr)
EXP$exp <- factor(EXP$exp, levels=c("std","brknn","brk91"))
df$id <- factor(df$id, levels=c("std","brknn","brk91"))

for (k in unique(EXP$exp)) {
  x <- df[df$id==k,]
  par(mar=c(0.25,5,0.25,0),xpd=T,font.axis = 1) #(DLUR)
  plot(c(1,2),c(0,100),type="n",xaxt="n",yaxt="n",bty="n",ylab = "\\n",xlab = "\\n")
  mtext(expression(italic(Xi)*" [%]"),side = 2,line = 2.5, cex = 1.5)
  axis(2, at=c(0,20,50,80,100),cex.axis=1.5, line=-0.25)
  for (j in unique(x$season)) {
    par(mar=c(0.25,0,0.25,0),xpd=T,font.axis = 1) #(DLUR)
    d=as.character(sites$st[sites$stid==i])
    g=as.character(EXP$exp_ppr[EXP$exp==k])
    plot(x$root[x$season==j],(100*x$kxyl[x$season==j]/2e-08),xlim=c(-10,0),pch=19,col=rgb(0.5,0.5,0.5,0.01),cex.lab=1,cex.axis=1,cex=3,ylim=c(0,100),xlab="\\n",ylab="\\n",xaxt="n",yaxt="n") #ylab="PLC [ % ]",xlab=expression(italic(psi)["stem"]*" [MPa]"))
    points(x$xyl[x$season==j],(100*x$kleaf[x$season==j]/2e-08),pch=15,col=rgb(0.13,0.67,0.13,0.01),cex=1.5)
    legend("topleft",legend = c("Leaves [L]","Xylem [X]","Frequency",5,25,50,">100"),bty="n",cex=1.75,pch=c(15,15,NA,16,16,16,16),col=c("forestgreen","grey",NA,rgb(0.5,0.5,0.5,0.05),rgb(0.5,0.5,0.5,0.25),rgb(0.5,0.5,0.5,0.5),rgb(0.5,0.5,0.5,1)))
    legend("bottomleft",legend = c(d,g),cex=1.75)
    legend("bottomright",legend = j,bty="n",cex=1.75)    
    par(mar=c(0.25,0,0.25,0),xpd=T,font.axis = 1) #(DLUR)
    boxplot(100*x$kxyl[x$season==j]/2e-08,100*x$kleaf[x$season==j]/2e-08,ylim=c(0,100),border = c("grey50","forestgreen"),col="grey90",lwd=2,yaxt="n",cex.lab=2,cex.axis=2,cex=1,pch=16)
  }
}
plot.new()
par(mar=c(6,0,0,0),xpd=T,font.axis = 1) #(DLUR)
plot(c(-10,0),c(1,2),type="n",yaxt="n",bty="n",ylab = "\\n",xlab = "\\n",xaxt="n")
mtext(expression(italic(psi)*" [MPa]"),side = 1,line = 3, cex = 1.5)
axis(1, at=c(-10,-5,0),cex.axis=1.5, line=-1)
plot.new()
par(mar=c(6,0,0,0),xpd=T,font.axis = 1) #(DLUR)
plot(c(-10,0),c(1,2),type="n",yaxt="n",bty="n",ylab = "\\n",xlab = "\\n",xaxt="n")
mtext(expression(italic(psi)*" [MPa]"),side = 1,line = 3, cex = 1.5)
axis(1, at=c(-10,-5,0),cex.axis=1.5, line=-1)
plot.new()
par(mar=c(6,0,0,0),xpd=T,font.axis = 1) #(DLUR)
plot(c(-10,0),c(1,2),type="n",yaxt="n",bty="n",ylab = "\\n",xlab = "\\n",xaxt="n")
mtext(expression(italic(psi)*" [MPa]"),side = 1,line = 3, cex = 1.5)
axis(1, at=c(-10,-5,0),cex.axis=1.5, line=-1)
plot.new()
par(mar=c(6,0,0,0),xpd=T,font.axis = 1) #(DLUR)
plot(c(-10,0),c(1,2),type="n",yaxt="n",bty="n",ylab = "\\n",xlab = "\\n",xaxt="n")
mtext(expression(italic(psi)*" [MPa]"),side = 1,line = 3, cex = 1.5)
axis(1,at=c(-10,-5,0),cex.axis=1.5, line=-1)
plot.new()

dev.off()
}

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#