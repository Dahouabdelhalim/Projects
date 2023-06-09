#Author: Joeri Reinders
#Purpose: Script which creates composite maps of soil moisture data for days with precipitation events do and do not result in flood events [given by an csv document] - figure 6 in paper.
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
library(patchwork)


rm(list=ls())

drive_e <- 'D:/'
drive_l <- 'C:/'
region <- "Texas"

directory <- "Users/joeri/OneDrive/Bureaublad/Northeastern University/Research Chapters/Chapter 4/"
data <- read.csv(paste(drive_l,directory,"Discharg_data/precip_peak_floods.csv",sep=""))

myPallette <- c(rev(brewer.pal(11, "RdBu")))
plot_list = list()
mean_fun <- function(data){
  resultaat <- mean(data,na.rm=T)
  return(resultaat)
} #
sd_fun <- function(data){
  resultaat <- sd(data,na.rm=T)
  return(resultaat)
} #
houston <- data.frame(x=c(-95.391), y=(c=29.5787))
level <- 500

# Soil Moisture load #
######################
directory <- "Data/CPC/"  #change this into the folder you need it to be
file_name <- paste('soilw.mon.mean.v2.nc',sep="")
var <- "soilw"

nc <- nc_open(paste(drive_e,directory,file_name, sep=""))
attributes(nc$var)

lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat", verbose = F)
time <- ncvar_get(nc, "time")

ndvi.array <- ncvar_get(nc, var) # store the data in a 3-dimensional array
dim(ndvi.array) 

fillvalue <- ncatt_get(nc, var, "_FillValue")
ndvi.array[ndvi.array == fillvalue$value] <- NA

####Monthly Anomoly
ndvi.normal.m <- array(NA, dim=c(720,360,3))
ndvi.normal.sd <- array(NA, dim=c(720,360,3))

dates <- seq(as.Date("1948-01-01"), as.Date("2020-12-31"), by = "months") 
index.min <- which(dates=="2020-12-01")
months_dates <- month(dates)
indexs <- which((months_dates==4)|(months_dates==5)|(months_dates==6))

ndvi.array <- ndvi.array[,,1:index.min]
ndvi.array <- ndvi.array[,,indexs]

dates <- dates[indexs]
date_dates <- format(dates, "%m")

dates_anom <- seq(as.Date("1971-01-01"), as.Date("2010-12-31"), by = "months")
index_anom.min <- which((dates==as.Date("1971-04-01")))-1
month_dates_anom <- month(dates_anom)
indexs_anom <- which((month_dates_anom==4)|(month_dates_anom==5)|(month_dates_anom==6))

dates_anom <- dates_anom[indexs_anom]
date_dates_anom <- format(dates_anom, "%m")

dates_ref <- format((seq(as.Date("1948-04-01"), as.Date("1948-06-30"), by = "months")),"%m") 

index_anom.min
length(date_dates_anom)
length(date_dates)

for(day in 1:3){
  date_index <- (which(dates_ref[day]==date_dates_anom)+index_anom.min)
  ndvi.normal <- ndvi.array[,,date_index]
  ndvi.normal.m[,,day] <- apply(ndvi.normal,1:2, mean_fun)
  ndvi.normal.sd[,,day] <- apply(ndvi.normal,1:2, sd_fun)
  
}

for(days in 1:3){
  date_index <- which(dates_ref[days]==date_dates)
  ndvi.array[,,date_index] <- sweep(ndvi.array[,,date_index],1:2,ndvi.normal.m[,,days])/replicate(length(date_index), ndvi.normal.sd[,,days], simplify="array")
}

ndvi.array.subset <- ndvi.array[500:600,90:130,]  
dim(ndvi.array.subset)

mybreak <- seq(-1,1,0.2)



directory <- "QGIS/HoustonRiver_wsheddel/"
watershed_brazos <- st_read(paste(drive_e,directory,"Brazos_Richmond_Watersheduf.shp", sep=""))
watershed_trinity <- st_read(paste(drive_e,directory,"Trinity_Romayor_Watershed.shp", sep=""))
watershed_neches <- st_read(paste(drive_e,directory,"Neches_Evadale_Watershed.shp", sep=""))

directory <- "QGIS/Figure_map_RM_2021s/"
houston_area  <- st_read(paste(drive_e,directory,"Houston_urbanarea.shp", sep=""))
houston_stad  <- st_read(paste(drive_e,directory,"houstonstad.shp", sep=""))
river_Trinity  <- st_read(paste(drive_e,directory,"Trinity_det.shp", sep=""))
river_Brazos  <- st_read(paste(drive_e,directory,"Brazos_det.shp", sep=""))
river_Neches  <- st_read(paste(drive_e,directory,"Neches_det.shp", sep=""))
river_Neches2  <- st_read(paste(drive_e,directory,"Neches2.shp", sep=""))
lakes  <- st_read(paste(drive_e,directory,"Lakes_map.shp", sep=""))

directory <- "QGIS/packages/Natural_Earth_quick_start/50m_cultural/"
coastline  <- st_read(paste(drive_e,directory,"ne_50m_admin_0_countries.shp", sep=""))
directory <- "QGIS/packages/Natural_Earth_quick_start/50m_physical/"
ocean  <- st_read(paste(drive_e,directory,"ne_50m_ocean.shp", sep=""))

gages <- data.frame("x"=c(-95.75, -94.84, -94.09, -98),"y"= c(29.58,30.45,30.35, 29))

blauw <- "dodgerblue4"

plot_title <- c("(a)   No Flood", "(b)   Flood")
floods <- 1:30
flood_data <- c(as.Date(data$NechesFlood), as.Date(data$TrinityFlood), as.Date(data$BrazosFlood))
no_flood_data <- c(as.Date(data$NechesNoFlood), as.Date(data$TrinityNoFlood), as.Date(data$BrazosNoFlood))

bootstrap_function <- function(B,n,data,ref,alpha){
  
  results.bootstrap <- array(data=NA, dim=c(101,41,B))
  
  for(bt in 1:B){
    dates_sampled <- sample(1:219,n,replace=TRUE)
    ndvi.slices.bootstrap <- data[,,dates_sampled]
    results.bootstrap[,,bt] <- apply(ndvi.slices.bootstrap,1:2, mean_fun)
  }
  
  
  results.bootstrap2 <- array(data=NA, dim=c(101,41))
  for(lon.b in 1:101){
    for(lat.b in 1:41){
      bina <- sum(abs(results.bootstrap[lon.b,lat.b,])>abs(ref[lon.b,lat.b]))/B  
      if(bina<alpha){results.bootstrap2[lon.b,lat.b] <- 1}
      else{results.bootstrap2[lon.b,lat.b] <- 0}
    }
  }
  
  return(results.bootstrap2)
}
alpha <- 0.01

for(plot in 1:2){

  if(plot==2){  
    for(index in 1:30){
    floods[index] <- which(format(dates,"%Y-%m")==format(as.Date(flood_data[index]),"%Y-%m"))
    }
  }
  
  if(plot==1){  
    for(index in 1:30){
      floods[index] <- which(format(dates,"%Y-%m")==format(as.Date(no_flood_data[index]),"%Y-%m"))
    }
  }
  
  ndvi.slices <- ndvi.array.subset[,,floods] #netcdfs of the picked days
  ndvi.slice <- apply(ndvi.slices,1:2, mean_fun)
  
  B <- 10
  results.bootstrap <- array(data=NA, dim=c(101,41,B))
  
  for(bt in 1:B){
    dates_sampled <- sample(1:219,30,replace=TRUE)
    ndvi.slices.bootstrap <- ndvi.array.subset[,,dates_sampled]
    results.bootstrap[,,bt] <- apply(ndvi.slices.bootstrap,1:2, mean_fun)
  }
  
  
  results.bootstrap2 <- array(data=NA, dim=c(101,41))
  for(lon.b in 1:101){
    for(lat.b in 1:41){
      if(is.na(ndvi.slice[lon.b,lat.b])==TRUE){results.bootstrap2[lon.b,lat.b] <- 0}else{
      bina <- sum(abs(results.bootstrap[lon.b,lat.b,])>abs(ndvi.slice[lon.b,lat.b]))/B  
      if(bina<alpha){results.bootstrap2[lon.b,lat.b] <- 1}
      else{results.bootstrap2[lon.b,lat.b] <- 0}}
    }
  }
  
  
  #results.bootstrap2 <- bootstrap_function(10,30,ndvi.array.subset,ndvi.slice,alpha)
  
  r <- raster(t(ndvi.slice), xmn=lon[500], xmx=lon[600], ymn=lat[130], ymx=lat[90],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # first the raster is molded into a dataframe which can be plotted in ggplot and the lot is adjusted for the world map
  flood <- raster::as.data.frame(r, xy = TRUE)
  cells <- length(flood$x)
  
  for(i in 1:cells){ 
    if(flood$x[i]<180){flood$x[i]<-flood$x[i]}
    if(flood$x[i]>180){flood$x[i]<-(flood$x[i]-359.75)}
  }
  flood$y <- flood$y-0.25
  
  #Prepares the bootstrap results in a raster and dataframe format for ggplot
  r.b <- raster(t(results.bootstrap2), xmn=lon[500], xmx=lon[600], ymn=lat[130], ymx=lat[90],crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  flood.b <- raster::as.data.frame(r.b, xy = TRUE)
  cells <- length(flood.b$x)
  for(i in 1:cells){ 
    if(flood.b$x[i]<180){flood.b$x[i]<-flood.b$x[i]}
    if(flood.b$x[i]>180){flood.b$x[i]<-(flood.b$x[i]-359.75)}
  }
  flood.b$y <- flood.b$y-0.25
  flood.b <- subset(flood.b, layer==1)
  
  j <- ggplot()+
        metR::geom_contour_fill(data=flood, aes(x=x, y=y, z=layer),binwidth = 0.01)+
        geom_point(data=flood.b, aes(x=x, y=y), color="black", size=0.8)+
        geom_sf(data= houston_area, colour="gray", fill="gray")+
        geom_sf(data= houston_stad, colour="black", fill="black")+
        geom_sf(data= river_Neches, colour=blauw, fill=blauw, size =0.8)+
        geom_sf(data= river_Neches2, colour=blauw, fill=blauw, size =0.8)+
        geom_sf(data= river_Brazos, colour=blauw, fill=blauw, size =0.8)+
        geom_sf(data= river_Trinity, colour=blauw, fill=blauw, size =0.8)+
        geom_sf(data= lakes, colour=blauw, fill=blauw)+
        geom_sf(data= coastline, colour="black", fill=NA, size =1)+
        geom_sf(data = watershed_brazos,colour = "black", fill = NA, alpha = 1, size =1)+
        geom_sf(data = watershed_trinity,colour = "black", fill = NA, alpha = 1, size =1)+
        geom_sf(data = watershed_neches,colour = "black", fill = NA, alpha = 1, size =1)+
        geom_sf(data= ocean, colour="black", fill="#CFDAED", size =1)+
        geom_point(data=gages,aes(x=x, y=y), size=3.5, shape=21, col="black", fill="#FDBF6F")+
        theme_minimal()+
        scale_fill_gradientn(colors =rev(myPallette), limits=c(-1,1), breaks=mybreak,
                             guide = guide_colorsteps(frame.colour = "black", ticks.colour = "black", barwidth=40, barheight=0.8, title.size=4, title.position="top",title="Standardized soil moisture anomoly",
                                                      title.hjust=0.5, title.theme=element_text(size=16),label.theme=element_text(size=13)))+
        coord_sf(xlim=c(-98,-94), ylim=c(29,32))+
        labs(x="Longitude", y="Latitude",fill=paste("", sep=""))+
        labs(title=plot_title[plot])+ 
        theme(legend.position = "bottom")+
        theme(plot.title = element_text(hjust = 0, face="bold", size=19))+
        theme(legend.margin=margin(20, 0, 0, 0))+
        theme(axis.title.x = element_text(size = 15),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 15),axis.text.y = element_text(size = 14))+
        annotate(geom="text", x=-95.47, y=30.64, label="   Trinity\\nStudy Area",color="black", size=4.5)+
        annotate(geom="text", x=-94.49, y=30.8, label="  Neches\\nStudy Area",color="black", size=4.5)+
        annotate(geom="text", x=-97.6, y=30.25, label=" Brazos\\nStudy Area",color="black", size=4.5)+
        annotate(geom="text", x=-95.55, y=29.9, label="Houston",color="black", size=4.5)+
        annotate(geom="text", x=-94.34, y=29.15, label="Gulf of Mexico",color="black", size=6)+
        annotate(geom="text", x=-97.4, y=29, label="Gauging Station",color="black", size=5.5)
  
  
  
  plot_list[[plot]] <- j

}


#CFDAED
combineda <- plot_list[[1]] + plot_list[[2]] & theme(legend.position = "bottom")
combineda + plot_layout(guides = "collect")
