#*******************************************************************************
# Cesar D. Jimenez-Rodriguez
# Luxembourg Institute of Science and Technology (LIST)
# Environmental Research and Innovation Department (ERIN)
# Environmental Sensing and Modelling Unit
# Agro-Environmental Systems Group
# E-mail: cesar.jimenez@list.lu 
#*******************************************************************************

#-Script Description------------------------------------------------------------
# The following script plots the index of agreement (IOA) proposed by Duveiller 
# et al. (2016). [doi =  https://doi.org/10.1038/srep19401].

#-Required Packages---------------------------------------------------------------------------------
library("readr")  
library(ncdf4)
library(sapfluxnetr)
library(dplyr)
library(lubridate)
library(tidyverse)
library(raster)
library(rgdal)
library(sp)
library(htmlwidgets)
library(RCurl)
library(XML)
library(rdwd)
library(qdapRegex)
library(plyr)
library(hexbin)
library(RColorBrewer)
library(ggplot2)
library(gganimate)
library(plotly)
library(ggpubr)
library(data.table)
library(plyr)
library(openair)
library(plotrix)
library(network)

#-User settings-----------------------------------------------------------------
rm(list = ls())

#-User settings-----------------------------------------------------------------
# Define absolute path of input and output files
res_path  = "~/Documents/CAPACITY/analysis/1_paper"
par_path = '~/Documents/CAPACITY/parameters'
#-Select ID of the SapFluxNet site --------------------------------------------#

site_id = c('ESP_ALT_TRI','FRA_HES_HE2_NON','FRA_PUE','RUS_FYO') # Site code from SapFluxNet.
stid <- c('ES-Alt','FR-Hes','FR-Pue','RU-Fyo')    # Site code from Fluxnet.
pftn= c(6,8,6,8)              # Plant functional type from the CLM list.
yri = c(2012,2001,2001,2001)
yrf = c(2014,2005,2011,2005)

mp <- data.frame(site_id,stid,pftn,yri,yrf)

expr <- c("std","brknn","brk91") # List of experiments to be evaluated according to the individual project code asisgned by the user.

#-Extra Functions --------------------------------------------------------------
# IOA corresponds to the Index of Agreement
IOA <- function(Tme,Tmo){
  r <- cor(Tme,Tmo,method = "pearson")
  alpha <- ifelse(r<0,0,2/((sd(Tme)/sd(Tmo))+(sd(Tmo)/sd(Tme))+(((mean(Tme)-mean(Tmo))^2)/(sd(Tme)*sd(Tmo)))))
  ioa <- r*alpha
return(c(r,ioa,alpha))
}

#-Data Retrieval----------------------------------------------------------------
pdf(paste(res_path,"/IOA_v2.pdf",sep=""),width = 10,height = 8)
layout(matrix(c(1:4),2,2, byrow = T), widths=c(1), heights=c(1),)
for (i in mp$stid) {
  j <- mp[mp$stid==i,]
yrk <- c(j$yri:j$yrf) # This is the range of years where the analysis is focused on.

wd <- paste(res_path,"/",as.character(j$site_id),sep="")

setwd(wd)

sap <- read.csv(paste(res_path,"/",as.character(j$site_id),"/","sap_",as.character(j$site_id),".csv",sep = ""))
sap$date <- sap$time
sap$doy <- as.numeric(strftime(strptime(substr(sap$date,6,10),format = "%m-%d"), format = "%j"))
sap$time <- paste(substr(sap$date,0,4),"-",sap$doy,sep="");sap[c("X","doy")]<-NULL

s <- list.files(pattern="*.csv$")
s <- s[lapply(s,function(x) length(grep("sap",x,value=FALSE))) == 0]

df_input_list <- lapply(s, fread)
names(df_input_list) <- gsub(s, pattern=paste(as.character(j$site_id),"_",sep=""), replacement="")

df <- bind_rows(df_input_list, .id = "id")
df$id <- gsub(df$id, pattern=".csv", replacement="")
df$doy <- df$doy+1
df$time <- paste(df$yr,"-",df$doy,sep="")
df <- merge(sap,df,by="time",all=T)
df$yr <- as.integer(substr(df$time,0,4))

df <- df[df$yr %in% yrk,]
df$stand[is.na(df$stand)] <- 0

df$id <- factor(df$id, levels = c("std","brknn","brk91"))
bsc <- df[df$id=="std",]

df$month <- as.numeric(substr(df$date,6,7))

x <-na.omit(df[c("id","month","stand","ET")])
w <- c("id","month","r","ioa","alpha")

for (e in unique(x$id)) {
  y <- x[x$id==e,]
  for (f in c(1:12)) {
   z <- y[y$month==f,]
   a <- c(e,f,round(IOA(z$stand,z$ET),2))
   w <- rbind(w,a)
  }
}
nms <- w[1,]
w<- as.data.frame(w)
colnames(w)<- nms
w <- w[-1,]

w$month <- as.numeric(as.character(w$month))
w$alpha <- as.numeric(as.character(w$alpha))
w$ioa <- as.numeric(as.character(w$ioa))
w$r <- as.numeric(as.character(w$r))
w$mdl <- as.character(w$id)
for (l in expr) {
  w$mdl[w$mdl==l] <- which(expr == l)
}

pch = c(18,15,16)
clr = c("azure4","darkcyan","firebrick2")
nm <- c("DMC","DBC","FBC")
stngs <- data.frame(expr,nm,pch,clr)

par(mar=c(3,4,0.25,0.25),xpd=F,font.axis = 1) #(DLUR)
plot(w$month,w$ioa,ylim=c(0,1.1),type="n",xaxt="n",ylab=expression(Gamma*" or r"))
axis(1,at = c(1:12),labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
legend("topleft",legend = j$stid,bty="n",cex=1.5)
legend("topright",legend=c(nm,expression(Gamma),"r"),col=c(as.character(clr),"black","black"),lty=c(1,1,1,1,3),lwd=2,bty="n",cex=1.15)
for (k in w$id) {
  t <- w[w$id==k,]
  points(t$month,t$r,pch=stngs$pch[stngs$expr==k],col=as.character(stngs$clr[stngs$expr==k]),type="l",lty=3,lwd=3)
  points(t$month,t$ioa,pch=stngs$pch[stngs$expr==k],col=as.character(stngs$clr[stngs$expr==k]),type="l",lwd=3,cex=2)
}
}
dev.off()

#-Cleaning RStudio -------------------------------------------------------------
rm(list=ls()) 
#-End of Script ---------------------------------------------------------------#