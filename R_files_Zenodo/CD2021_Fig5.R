#Author: Joeri Reinders
#Purpose: Script which creates correlation analyses of precipitation over the Brazos, Trinity and Neches rivers, and geopotentialheight and wind field reanalysis data for days of flooding [given by an csv document] - figure 5 in paper.
# Also constructs the Great Plains Low Level Jet (GPLLJ) index and the Western Bermuda High (WBH) index and correlates it with the above reanalysis data.
#Data: contact reinders.j@northeastern.edu for orginal data files

rm(list=ls())
x11()

library(ncdf4)
library(ggplot2)
library(gridExtra)
library(raster)
library(sf)
library(tmap)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(R.matlab)
library(lubridate)
library(metR)
library(cowplot)
library(RColorBrewer)
library(dplR)
library(extRemes)
library(multiApply)
library(patchwork)
library(boot)

rm(list=ls())

drive_e <- 'D:/'
drive_l <- 'C:/'
region <- "Texas"
NASH_index <- data.frame(x=c(-75,-92), y=c(30,30)) #30°N and 75°W and 30°N and 92 °W 

myPallette <- c(rev(brewer.pal(11, "RdBu")))
plot_list = list()

## Precipitation load ##
########################
directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
precipdata <- read.csv(paste(drive_l,directory,"Discharg_data/Precip_RiversWatersheds.csv",sep=""))

#create a column with the sum of rain per watershed - (standardized by rastersize!)
for(X in 1:length(precipdata$Neches)){precipdata$sum[X] <- sum((precipdata$Neches[X]*2),(precipdata$Trinity[X]*2),(precipdata$Brazos[X]*53))/(53+2+2)}

#extract spring days
dates <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "days") 
months <- format(dates, "%m")
dates_lente <- which((months=="04")|(months=="05")|(months=="06"))
spring_precipdata <- precipdata[dates_lente,]

#7 day total rainfall and maximum
day7sum_precip <- 1:(length(spring_precipdata$sum)/7)
for(p7 in 1:length(day7sum_precip)){
  if(p7==1){
    A <- p7
    B <- p7+6
  }
  day7sum_precip[p7] <- sum(spring_precipdata$sum[A:B])
  A <- A+7
  B <- B+7
}

day7max_precip <- 1:(length(spring_precipdata$sum)/7)
for(p7 in 1:length(day7max_precip)){
  if(p7==1){
    A <- p7
    B <- p7+6
  }
  day7max_precip[p7] <- max(spring_precipdata$sum[A:B])
  A <- A+7
  B <- B+7
}

#monthly total rainfall
days.in.month <- rep(c(30,31,30), 73)
monthsum_precip <- 1:(73*3)
for(pm in 1:length(monthsum_precip)){ #length(monthsum_precip)){
  if(pm==1){
    A <- 1
    B <- 30
  }
  monthsum_precip[pm] <- sum(spring_precipdata$sum[A:B])
  A <- B+1
  B <- A+days.in.month[pm+1]-1
}

monthmax_precip <- 1:(73*3)
for(pm in 1:length(monthsum_precip)){ #length(monthsum_precip)){
  if(pm==1){
    A <- 1
    B <- 30
  }
  monthmax_precip[pm] <- max(spring_precipdata$sum[A:B])
  A <- B+1
  B <- A+days.in.month[pm+1]-1
}

## Geopotential Height load ##
##############################
level <- 850
directory <- "Data/Reanalysis/NCEP-NCAR/GeopotentialHeight/daily_data/"  #change this into the folder you need it to be
file_name <- paste('All.',level,'.nc',sep="")
var <- "hgt"

nc <- nc_open(paste(drive_e,directory,file_name, sep=""))
attributes(nc$var)

lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat", verbose = F)
time <- ncvar_get(nc, "time")

ndvi.array <- ncvar_get(nc, var) # store the data in a 3-dimensional array
fillvalue <- ncatt_get(nc, var, "_FillValue")
ndvi.array[ndvi.array == fillvalue$value] <- NA

ndvi.array <- ndvi.array[,,dates_lente]

## V-wind load ##
#################
w.level <- 925
directory_vwind <- "Data/Reanalysis/NCEP-NCAR/Vwinds/daily_data/"  #change this into the folder you need it to be
file_name_vwind <- paste('All.',w.level,'.nc',sep="")
var_vwind <- "vwnd"

nc_vwind <- nc_open(paste(drive_e,directory_vwind,file_name_vwind, sep=""))

vlon <- ncvar_get(nc_vwind, "lon")
vlat <- ncvar_get(nc_vwind, "lat", verbose = F)
vtime <- ncvar_get(nc_vwind, "time")

v.ndvi.array <- ncvar_get(nc_vwind, var_vwind) # store the data in a 3-dimensional array
fillvalue <- ncatt_get(nc_vwind, var_vwind, "_FillValue")
v.ndvi.array[v.ndvi.array == fillvalue$value] <- NA

v.ndvi.array <- v.ndvi.array[,,dates_lente]

#7 day averages
mean_fun <- function(data){
  resultaat <- mean(data,na.rm=T)
  return(resultaat)
} #

day7averageGPH <- array(data=NA, dim=c(144,73,(length(spring_precipdata$Neches)/7)))
for(g7 in 1:length(day7averageGPH[1,1,])){
  if(g7==1){
    A <- g7
    B <- g7+6
  }
  day7averageGPH[,,g7] <- apply(ndvi.array[,,A:B],1:2, mean_fun)
  A <- A+7
  B <- B+7
}

day7averageVW <- array(data=NA, dim=c(144,73,(length(spring_precipdata$Neches)/7)))
for(v7 in 1:length(day7averageVW[1,1,])){
  if(v7==1){
    A <- v7
    B <- v7+6
  }
  day7averageVW[,,v7] <- apply(v.ndvi.array[,,A:B],1:2, mean_fun)
  A <- A+7
  B <- B+7
}

#Montlhy averages
monthsum_GPH <- array(data=NA, dim=c(144,73,(73*3)))
for(gm in 1:length(monthsum_GPH[1,1,])){
  if(gm==1){
    A <- 1
    B <- 30
  }
  monthsum_GPH[,,gm] <- apply(ndvi.array[,,A:B],1:2, mean_fun)
  
  if(gm<length(monthsum_GPH[1,1,])){
    A <- B+1
    B <- A+days.in.month[gm+1]-1}
}

monthsum_VW <- array(data=NA, dim=c(144,73,(73*3)))
for(vm in 1:length(monthsum_VW[1,1,])){ 
  if(vm==1){
    A <- 1
    B <- 30
  }
  monthsum_VW[,,vm] <- apply(v.ndvi.array[,,A:B],1:2, mean_fun)
  
  if(vm<length(monthsum_VW[1,1,])){
    A <- B+1
    B <- A+days.in.month[vm+1]-1}
}

correlation.data.function = list(day7sum_precip, day7sum_precip)
correlation.data.input = list(day7averageGPH, day7averageVW)
titelplot <- c("(a)  Weekly Total Precipitation ~ Weekly 850hPA Geopotential Height", "(b)  Weekly Total Precipitation ~ Weekly 925hPA Meridional Wind")

# GP LLJ index coordinates are:longitude-latitude box (25°-35°N, 102°-97°W)
#PLOT LOOP
for(plot in 1:2){
  correlatie <- function(x){
    corra1 <- cor.test(x,correlation.data.function[[plot]], method="pearson")
    corra2 <- as.numeric(corra1$estimate)
    return(corra2)} #here you need to change what to corrolate with what! 
  correlatie.p <- function(x){
    x <- x
    y <- correlation.data.function[[plot]]
    dat <- data.frame(x, y)
    
    bootCorTest <- function(data, i){
      d <- data[i, ]
      cor.test(d$x, d$y)$p.value
    }
    
    b <- boot(dat, bootCorTest, R = 10)
    
    corra2 <- mean(b$t)
    return(corra2)} #here you need to change what to corrolate with what! 
  
  corr_map <- data.matrix(as.data.frame(Apply(correlation.data.input[[plot]], 3, correlatie)))
  corr_map.p <- data.matrix(as.data.frame(Apply(correlation.data.input[[plot]], 3, correlatie.p)))

  r.corr <- raster(t(corr_map), xmn=lon[1], xmx=lon[144], ymn=lat[1], ymx=lat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  #r.corr <- flip(r.corr,"y")
  r.corr.p <- raster(t(corr_map.p), xmn=lon[1], xmx=lon[144], ymn=lat[1], ymx=lat[73],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  df.corr <- raster::as.data.frame(r.corr, xy = TRUE)
  df.corr.p <- raster::as.data.frame(r.corr.p, xy = TRUE)
  df.corr$layerp <- df.corr.p$layer
  cells <- length(df.corr$x)
  for(i in 1:cells){ 
    if(df.corr$x[i]<180){df.corr$x[i]<-df.corr$x[i]-((360/144)/2)}
    if(df.corr$x[i]>180){df.corr$x[i]<-df.corr$x[i]-360+((360/144)/2)}
  }
  
  df.corr_point <- subset(df.corr, (layer>0.1 | layer < -0.1))
  df.corr_point <- subset(df.corr_point, layerp<0.01)
  
  mybreak <- c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)
  
  q <- ggplot()+
    metR::geom_contour_fill(data=df.corr, aes(x=x, y=y, z=layer),binwidth = 0.05)+
    geom_point(data=df.corr_point, aes(x=x, y=y), color="black", size=0.5)+
    borders('world', xlim=range(df.corr$x), ylim=range(df.corr$y), colour='black', size=.2)+
    borders('state', colour="black", size=0.1)+
    theme_minimal()+ 
    theme(panel.ontop=TRUE, panel.background=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    {if(plot==1)geom_line(data=NASH_index, aes(x=x, y=y), color="red", size=1.5)}+
    {if(plot==2)geom_rect(aes(xmin = -102, xmax = -97, ymin = 25, ymax = 35), fill = "transparent", color = "red", size = 1.25)}+
    {if(plot==2)geom_rect(aes(xmin = -92.5, xmax = -87.5, ymin = 28, ymax = 33), fill = "transparent", color = "red", size = 1.25)}+
    scale_fill_gradientn(colors =myPallette, limits=c(-0.50,0.50), breaks=mybreak,
                         guide = guide_colorsteps(frame.colour = "black", ticks.colour = "black", barwidth=32, barheight=0.8, title.size=4, title.position="top",title=expression(paste("Pearson's ",italic("r"))), 
                                                  title.hjust=0.5, title.theme=element_text(size=16),label.theme=element_text(size=13)))+
    coord_quickmap(xlim=c(-140,-10), ylim=c(0,70))+
    labs(x="Longitude", y="Latitude",fill=paste("", sep=""))+
    theme(legend.position = "bottom")+
    theme(legend.margin=margin(20, 0, 0, 0))+
    theme(plot.title = element_text(hjust = 0, face="bold", size=13))+
    theme(axis.title.x = element_text(size = 15),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 15),axis.text.y = element_text(size = 14))+
    labs(title=paste(titelplot[plot], sep=""))
  
  plot_list[[plot]] <- q
  
}


# Monthly Anomolies Rain, GPLLJ, WBHI #
#######################################

index_wijzer <- 1:219

#Monthly anomolies precip:
monthsum_precip_anom <- 1:length(monthsum_precip)
for(a in 1:3){
  values <- monthsum_precip[seq(a, length(monthsum_precip), 3)]
  index.values <- index_wijzer[seq(a, length(index_wijzer), 3)]
  precip.normal.m <- mean(values)
  precip.normal.sd <- sd(values)
  
  monthsum_precip_anom[index.values] <- (monthsum_precip[index.values]-precip.normal.m)/precip.normal.sd
}

names_plotb <- c("Montly Maximum Precipitation ~ Great Plain Low Level Jet Index","Montly Maximum Precipitation ~ Western Bermuda High Index" )

for(plotb in 1:2){
  
  if(plotb==1){
    #To facilitate analysis of GPLLJ variability, an index
    #is constructed from areal averaging of the meridional
    #wind in a 5° 10° longitude-latitude box (25°-35°N,
    #                                          102°-97°W). 
    #The GPLLJ index is
    #finally defined as the box-averaged 900-hPa meridional
    #wind. (Scott and Weaver 2008)
    
    #Monthly anomolies GPLLJ:
    GPLLJ.index <- apply(monthsum_VW[104:106,23:27,],3, mean_fun)
    index.variable <- GPLLJ.index
  }

  if(plotb==2){
    #The BHI is an index of the difference between normalized sea-level pressure over Bermuda and New Orleans,
    #Louisiana (Stahle and Cleve- land, 1992); thus, monthly sea-level pressure values were acquired for grid cells 
    #corresponding to 32.5 °N and 64 °W and 30°N and 90 °W. The centre of the Bermuda High is typically located approximately 
    #2500 km east of Bermuda during summer; therefore, the BHI shows the standardized pressure gradient across the western side 
    #of the Bermuda High. The
    
    #The WBHI is a newly developed index based on the BHI, with the major difference being that it uses 850-hPa heights and is
    #centred over the southeastern United States. Monthly 850-hPa geopoten- tial heights were acquired for grid cells corresponding
    #to 30°N and 75°W and 30°N and 92 °W (Figure 1).  (Diem,2013)
    
    #Monthly anomolies WBHI:
    lon
    WBHI_A <- monthsum_GPH[115,25,]
    WBHI_B <- monthsum_GPH[108,25,]
    WBHI <- (WBHI_A-WBHI_B)
    WBHI.norm <- 1:219
    for(a in 1:3){
      values <- WBHI[seq(a, length(WBHI), 3)]
      index.values <- index_wijzer[seq(a, length(index_wijzer), 3)]
      WBHI.m <- mean(values)
      WBHI.sd <- sd(values)
      WBHI.norm[index.values] <- (WBHI[index.values]-WBHI.m)/WBHI.sd
    }
    index.variable <- WBHI.norm
  }
  
  dates_m <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "month") 
  months_m <- format(dates_m, "%m")
  dates_lente_m <- which((months_m=="04")|(months_m=="05")|(months_m=="06"))
  
  low.pass_filter <- function(data){
    resultaat <- pass.filt(data, W=0.1, type="low", method="Butterworth" )
    return(resultaat)
  }
  filtered.index <- low.pass_filter(index.variable)
  filtered.Precipm <- low.pass_filter(monthsum_precip_anom)
  
  df_gpllj <- data.frame(monthsum_precip_anom,index.variable,filtered.Precipm,filtered.index)
  df_gpllj$date <- dates_m[dates_lente_m]
  
  
  cor <- cor.test(monthsum_precip,index.variable,method="pearson")
  
  j <- ggplot()+
    geom_line(data=df_gpllj, aes(x=date,y= index.variable),color="black",size=0.5)+
    geom_line(data=df_gpllj, aes(x=date,y= monthsum_precip_anom), color="red",size=0.5)+
    annotate(geom="text",x=df_gpllj$date[18], y=2.5, label=paste("R= ", round(as.numeric(cor$estimate),3),sep=""), size= 5)+
    theme_classic()+
    labs(x="Date", y="Standardized Anomoly")+
    theme(plot.title = element_text(hjust = 0, face="bold", size=15))+
    theme(axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 13),axis.title.y = element_text(size = 14),axis.text.y = element_text(size = 13))+
    labs(title=paste(names_plotb[plotb], sep=""))
  
  t <- ggplot()+
    geom_point(data=df_gpllj, aes(x=monthsum_precip_anom,y=index.variable), color="black",size=1.0)+
    geom_smooth(data=df_gpllj, aes(x=monthsum_precip_anom,y=index.variable), method=lm)+
    annotate(geom="text",x=-1.5, y=2.4, label=paste("R= ", round(as.numeric(cor$estimate),3),sep=""), size= 5)+
    annotate(geom="text",x=-0.5, y=2.4, label=paste("(P-value = ", round(as.numeric(cor$p.value),3), ")",sep=""), size= 4.5)+
    theme_classic()+
    xlim(-2, 2.7)+
    ylim(-2, 2.7)+
    labs(x="Precipitation Anomily", y="Index Anomoly")+
    theme(plot.title = element_text(hjust = 0, face="bold", size=13))+
    theme(axis.title.x = element_text(size = 14),axis.text.x = element_text(size = 13),axis.title.y = element_text(size = 14),axis.text.y = element_text(size = 13))+
    labs(title=paste(names_plotb[plotb], sep=""))
  
  plot_list[[plotb+2]] <- t
}

combinedA <- (plot_list[[1]] + plot_list[[2]]) & theme(legend.position = "bottom")
combinedB <- (plot_list[[4]] + plot_list[[3]]) + plot_layout(widths = unit(c(15), c('cm'))) & theme(legend.position = "bottom")
CA <- combinedA + plot_layout(guides = "collect")
CB <- combinedB + plot_layout(guides = "collect") 
CA

