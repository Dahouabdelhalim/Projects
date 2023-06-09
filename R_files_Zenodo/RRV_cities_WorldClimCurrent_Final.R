##### Sadie Ryan, University of Florida
##### Updated April 2018

##### Purpose: Overlay RRV City lat/long coordinates onto Current WorldClim for 12 months of min, max, mean temps

library(raster)
library(sp)
library(rgdal)
library(dismo)

setwd("E:/Dropbox/RRV R0")
#Read in RRV points
RRVCities<-read.csv("CitiesForWorldClim.csv")
coordinates(RRVCities) <- c("Long", "Lat")
projection(RRVCities) <- CRS("+proj=lonlat +ellps=WGS84")

australia.mrc <- gmap("Australia")

setwd("R:/Ryan_Lab/RRV")

##CURRENT temp data (WorldClim)

#maximum temperature Data
w <- getData('worldclim', var='tmax', res=5)
maxTCur_Stack10<-w*0.1

#Minimum temperature data
x <- getData('worldclim', var='tmin', res=5)
minTCur_Stack10<-x*0.1

#Mean temperature data
x <- getData('worldclim', var='tmean', res=5)
meanTCur_Stack10<-x*0.1

plot(meanTCur_Stack10[[1]])
plot(RRVCities, col="red", add=TRUE)

#Overlay points on rasters (bricks)
RRVCitiesMax<-extract(maxTCur_Stack10, RRVCities)
RRVCitiesMin<-extract(minTCur_Stack10, RRVCities)
RRVCitiesMean<-extract(meanTCur_Stack10, RRVCities)

RRVMax<-data.frame(RRVCities$City,RRVCitiesMax)
RRVMin<-data.frame(RRVCities$City,RRVCitiesMin)
RRVMean<-data.frame(RRVCities$City,RRVCitiesMean)

setwd("E:/Dropbox/RRV R0")
#Write out CSV of points with 12 months x 3
write.csv(RRVMax, "RRVCitiesmax.csv")
write.csv(RRVMin, "RRVCitiesmin.csv")
write.csv(RRVMean, "RRVCitiesmean.csv")
