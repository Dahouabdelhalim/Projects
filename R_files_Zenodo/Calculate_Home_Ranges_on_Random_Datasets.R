##### This is the latest script for computing home ranges for Cape gannets
### This script used a single csv file with all the tracked individuals and then sperates them according to year and sex
### Additional features include functions to calculate area of home ranges and compute overlap indices

# JA Botha 06 April 2017

library(adehabitatHR)
library(sp)
require(rgdal)
library(plyr)


setwd("C:/Users/Jonathan/Desktop/Ale-Resource Partitioning/HomeRange/RandomDatasets")
rm(list=ls()) # clear data list
# Get everything loaded
files<-dir(pattern="HR",full.names=T)
list(files)

for(i in 2:length(files)){
  
   locs <-read.csv(files[i],sep=",",header=T) # Full data set containing all the individuals for each sex in each year
  head(locs)
  
  locs<- rename(locs,c("dta.Date"="Date","dta.Species"="Species","dta.Lat"="Lat", "dta.Long"="Long")) 
  
  #Basic plot
  #plot(x=locs$Long,y=locs$Lat,pch=3, col="blue",ylab="Lat",xlab="Long")
  #plot(locs$Long,locs$Lat)
  
  # Create a SpatialPointsDataFrame
  sptracks <- SpatialPointsDataFrame(coords = as.matrix(locs[,c(4,3)]), 
                                     data = data.frame(locs[,c(5,1,4,3)]), 
                                     coords.nrs = numeric(0), 
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
str(sptracks)
names(sptracks)
#plot(sptracks)

#Plotting spatial data to check it is correct.
#plot(sptracks, col="blue", pch=19)
#to make plotting extent for khr (grid argument i.e whole map dimensions)
ext<-extent(c(22.60, 24.4, -34.32,-33.76)) #xmin, xmax, ymin, ymax

#making ext a raster called r
r<-raster(nrows=(ext[4]-ext[3])*1000, ncol=(ext[2]-ext[1])*1000, ext=ext,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#making r a spatial pixels class called sp
sp<-as(r,"SpatialPixels")

# Calculate the UD using the kernel based method 
kud <- kernelUD(sptracks[,1], h="LSCV",grid=sp,hlim=c(0.01,10), kern="bivnorm", extent=5)
image(kud)

kud[[1]]@h
kud[[2]]@h

#image(kudl)

ii <- estUDm2spixdf(kud)
class(ii)

vud <- getvolumeUD(kud)
vud

ii <- kernel.area(kud, percent=seq(50, 95, by=5))
ii

# Compute the overlap for various density contours
overlap95 <- kerneloverlaphr(kud, method = c("UDOI"), percent = 95, conditional = FALSE)          
overlap95
OV95 <- overlap95[2,"BND"]
OV95

overlap50 <- kerneloverlaphr(kud, method = c("UDOI"), percent = 50, conditional = FALSE)          
overlap50
OV50 <- overlap50[2,"BND"]
OV50

# Get the 95% density contours
#ALL95 <- getverticeshr(kud,percent=95)
BND95 <- getverticeshr(kud[[1]],percent=95)
HBD95  <- getverticeshr(kud[[2]],percent=95)

# Get the area for each 95% density contour
Area95BND <- as.data.frame(BND95)
Area95BND
BND.area95 <- Area95BND[1,"area"]
BND.area95
Area95HBD <- as.data.frame(HBD95)
Area95HBD
HBD.area95 <- Area95HBD[1,"area"]
HBD.area95

# Get the 50% density contours
BND50 <- getverticeshr(kud[[1]],percent=50)
HBD50  <- getverticeshr(kud[[2]],percent=50)

# Get the area for each 95% density contour
Area50BND <- as.data.frame(BND50)
Area50BND
BND.area50 <- Area50BND[1,"area"]
BND.area50
Area50HBD <- as.data.frame(HBD50)
Area50HBD
HBD.area50 <- Area50HBD[1,"area"]
HBD.area50


addparameters<-data.frame(OV95,OV50,BND.area95,HBD.area95,BND.area50,HBD.area50)

parameters<-rbind(parameters,addparameters)

}

write.csv(parameters, "RandomOvelaps.csv", row.names=T)

