# Calculate home ranges for bottlenose and humpback dolphins
# Data for both species is in one csv file
library(adehabitatHR)
library(sp)
require(rgdal)

setwd("C:/Users/Jonathan/Desktop/Ale-Resource Partitioning/HomeRange")
rm(list=ls()) # clear data list
# Get everything loaded
files<-dir(pattern="DolphinLocs",full.names=T)
list(files)

locs <-read.csv(files[1],sep=",",header=T) # Full data set containing all the individuals for each sex in each year
head(locs)

#Basic plot
plot(x=locs$Long,y=locs$Lat,pch=3, col="blue",ylab="Lat",xlab="Long")
plot(locs$Long,locs$Lat)

#filename <- as.factor(substr(basename(files[1]),1,8))   # Creates a filename which is needed when exporting cs
plot(locs$Long,locs$Lat)

# Create a SpatialPointsDataFrame
sptracks <- SpatialPointsDataFrame(coords = as.matrix(locs[,c(4,3)]), 
                                   data = data.frame(locs[,c(2,1,4,3)]), 
                                   coords.nrs = numeric(0), 
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"), match.ID = TRUE)

str(sptracks)
names(sptracks)
plot(sptracks)


#Plotting spatial data to check it is correct.
plot(sptracks, col="blue", pch=19)
#to make plotting extent for khr (grid argument i.e whole map dimensions)
ext<-extent(c(22.60, 24.4, -34.32,-33.76)) #xmin, xmax, ymin, ymax

#making ext a raster called r
r<-raster(nrows=(ext[4]-ext[3])*1000, ncol=(ext[2]-ext[1])*1000, ext=ext,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#making r a spatial pixels class called sp
sp<-as(r,"SpatialPixels")

# Calculate the UD using the kernel based method 
#khr<-kernelUD(xy=spatpt, h="LSCV", grid=sp, hlim=c(0.01,10), kern="bivnorm", extent=5)

kud <- kernelUD(sptracks[,1], h="LSCV",grid=sp,hlim=c(0.01,10), kern="bivnorm", extent=5)

image(kud)

kud[[1]]@h
kud[[2]]@h

#kudl <- kernelUD(sptracks[,1], h="LSCV")
image(kud)

ii <- estUDm2spixdf(kud)
class(ii)

vud <- getvolumeUD(kud)
vud

ii <- kernel.area(kud, percent=seq(50, 95, by=5))
ii

# Compute the overlap for various density contours
overlap95 <- kerneloverlaphr(kud, method = c("UDOI"), percent = 95, conditional = FALSE)          
overlap95

overlap50 <- kerneloverlaphr(kud, method = c("UDOI"), percent = 50, conditional = FALSE)          
overlap50

# Export the overlap information
write.csv(overlap95, "Overlap_Indices_95.csv", row.names=T)
write.csv(overlap50, "Overlap_Indices_50.csv", row.names=T)

############### Creating the shapefiles ##############

# Get the 95% density contours
ALL95 <- getverticeshr(kud,percent=95)
UDBND95 <- getverticeshr(kud[[1]],percent=95)
UDHBD95  <- getverticeshr(kud[[2]],percent=95)

# Plot the 95% density contours
plot(ALL95, col=1:4)
plot(UDBND95 , col=3)
plot(UDHBD95, col=4)


# Get the area for each 95% density contour
Area95 <- as.data.frame(ALL95)
Area95

write.csv(Area95, "Contour_Area_95.csv", row.names=T)

# Export the 95% density contours as a shapefile
writeOGR(UDBND95,"shapes","testShape",dsn = "C:/Users/Jonathan/Desktop/Ale-Resource Partitioning/HomeRange/Shapefiles", 
         layer = "BND-HR95",driver="ESRI Shapefile")

writeOGR(UDHBD95,"shapes","testShape",dsn = "C:/Users/Jonathan/Desktop/Ale-Resource Partitioning/HomeRange/Shapefiles", 
         layer = "HBD-HR95",driver="ESRI Shapefile")

# Get the 50% density contours
ALL50 <- getverticeshr(kud,percent=50)
UDBND50 <- getverticeshr(kud[[1]],percent=50)
UDHBD50  <- getverticeshr(kud[[2]],percent=50)

# Plot the 95% density contours
plot(ALL50, col=1:2)
plot(UDBND50, col=3)
plot(UDHBD50, col=4)

# Get the area for each 95% density contour
Area50 <- as.data.frame(ALL50)
Area50
write.csv(Area50, "Contour_Area_50.csv", row.names=T)

writeOGR(UDBND50,"shapes","testShape",dsn = "C:/Users/Jonathan/Desktop/Ale-Resource Partitioning/HomeRange/Shapefiles", 
         layer = "BND-HR50",driver="ESRI Shapefile")

writeOGR(UDHBD50,"shapes","testShape",dsn = "C:/Users/Jonathan/Desktop/Ale-Resource Partitioning/HomeRange/Shapefiles", 
         layer = "HBD-HR50",driver="ESRI Shapefile")

