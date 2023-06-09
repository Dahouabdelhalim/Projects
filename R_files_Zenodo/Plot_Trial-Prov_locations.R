#This script produces a map showing all the provenances and trials compiled by the BeechCOSTe52

#Required libraries
library(rgdal)
library(rworldmap)
library(ggmap)
library(ggplot2)
library(lattice)
data(coastsCoarse)
data(countriesLow)
library(plyr)

#set working directory 
setwd("/WorkingDirectory/")

#Read databases
Fsyl <- read.csv("Fsylvatica.csv")
Trials <- read.csv("Trial_coords.csv")
Prov <- read.csv("Prov_coords.csv")

#Rename coordinates 
colnames(Trials)[14] <- "lonTrial"
colnames(Trials)[15] <- "latTrial"
colnames(Prov)[13] <- "lonProv"
colnames(Prov)[14] <- "latProv"

#Add tree age for each measure
Fsyl$H95cmAge <- (1995 - Fsyl$YearPlantation)
Fsyl$H96cmAge <- (1996 - Fsyl$YearPlantation)
Fsyl$H97cmAge <- (1997 - Fsyl$YearPlantation)
Fsyl$H98cmAge <- (1998 - Fsyl$YearPlantation)
Fsyl$H99cmAge <- (1999 - Fsyl$YearPlantation)
Fsyl$H00cmAge <- (2000 - Fsyl$YearPlantation)
Fsyl$H01cmAge <- (2001 - Fsyl$YearPlantation)
Fsyl$H02cmAge <- (2002 - Fsyl$YearPlantation)
Fsyl$H03cmAge <- (2003 - Fsyl$YearPlantation)
Fsyl$H04cmAge <- (2004 - Fsyl$YearPlantation)
Fsyl$H05cmAge <- (2005 - Fsyl$YearPlantation)
Fsyl$H06cmAge <- (2006 - Fsyl$YearPlantation)
Fsyl$H07cmAge <- (2007 - Fsyl$YearPlantation)
Fsyl$H08cmAge <- (2008 - Fsyl$YearPlantation)


#Merge geographical coordinates of trials and provenances with phenotypic data
new = merge(Fsyl,Trials, "Trial")
data = merge(new, Prov, by="ID_ProvCode")

#Add Europe layer
wmap_laea<- spTransform(countriesLow, CRS("+proj=laea")) 

names(data)
dataBU19 <- subset(data, (SERIE=="BU19"))
dataBU20 <- subset(data, (SERIE=="BU20"))

dataprovBU19 <-cbind(dataBU19$latProv,dataBU19$lonProv)
dataprovBU20 <-cbind(dataBU20$latProv,dataBU20$lonProv)

datatrialsBU19 <-cbind(dataBU19$latTrial,dataBU19$lonTrial)
datatrialsBU20 <-cbind(dataBU20$latTrial,dataBU20$lonTrial)

dataprovBU19 <-as.data.frame(dataprovBU19)
dataprovBU20 <-as.data.frame(dataprovBU20) 
 
datatrialsBU19 <-as.data.frame(datatrialsBU19)
datatrialsBU20 <-as.data.frame(datatrialsBU20)

names(datatrialsBU19) <-c("latitude","longitude")
names(datatrialsBU20) <-c("latitude","longitude")

names(dataprovBU19) <-c("latitude","longitude")
names(dataprovBU20) <-c("latitude","longitude")

#Reproject trials 
coordinates(datatrialsBU19)<-c("longitude","latitude")
proj4string(datatrialsBU19) <- CRS("+proj=longlat")
datatrialsBU19_laea<-spTransform(datatrialsBU19, CRS("+proj=laea"))

coordinates(datatrialsBU20)<-c("longitude","latitude")
proj4string(datatrialsBU20) <- CRS("+proj=longlat")
datatrialsBU20_laea<-spTransform(datatrialsBU20, CRS("+proj=laea"))

countries_df<-fortify(countriesLow)
datatrialsBU19_df<-data.frame(datatrialsBU19)
datatrialsBU20_df<-data.frame(datatrialsBU20)

#Reproject Europe map
wmap_laea_df<-fortify(wmap_laea)
datatrialsBU19_laea_df<-data.frame(datatrialsBU19_laea)
datatrialsBU20_laea_df<-data.frame(datatrialsBU20_laea)

#Extract georaphical tile from trial coordinates
xmin<-min(datatrialsBU20_laea_df$longitude)
xmax<-max(datatrialsBU20_laea_df$longitude)
ymin<-min(datatrialsBU20_laea_df$latitude)
ymax<-max(datatrialsBU20_laea_df$latitude)
buff<-700000


#Reproject trials 
dataprovBU19<-  subset(dataprovBU19,(latitude != "NA" ))
dataprovBU20<-  subset(dataprovBU20,(latitude != "NA" ))

coordinates(dataprovBU19)<-c("longitude","latitude")
proj4string(dataprovBU19) <- CRS("+proj=longlat")
dataprovBU19_laea<-spTransform(dataprovBU19, CRS("+proj=laea"))

coordinates(dataprovBU20)<-c("longitude","latitude")
proj4string(dataprovBU20) <- CRS("+proj=longlat")
dataprovBU20_laea<-spTransform(dataprovBU20, CRS("+proj=laea"))

countries_df<-fortify(countriesLow)
dataprovBU19_df<-data.frame(dataprovBU19)
dataprovBU20_df<-data.frame(dataprovBU20)


#reprojected data
wmap_laea_df<-fortify(wmap_laea)
dataprovBU19_laea_df<-data.frame(dataprovBU19_laea)
dataprovBU20_laea_df<-data.frame(dataprovBU20_laea)

#Extract geographical tile from provenance coordinates
xmin<-min(dataprovBU20_laea_df$longitude)
xmax<-max(dataprovBU20_laea_df$longitude)
ymin<-min(dataprovBU20_laea_df$latitude)
ymax<-max(dataprovBU20_laea_df$latitude)
buff<-700000


#Read and reproject distribution range from Euforgen map
shapefile <- readOGR(dsn="Rcode", "Fagus_sylvatica_EUFORGEN")
shapefile_df <- fortify(shapefile)
shapefile_laea<- spTransform(shapefile, CRS("+proj=laea")) 
shapefile_df_laea <- fortify(shapefile_laea)


xmin<-min(dataprovBU20_laea_df$longitude)
xmax<-max(dataprovBU20_laea_df$longitude)
ymin<-min(dataprovBU20_laea_df$latitude)
ymax<-max(dataprovBU20_laea_df$latitude)
buff<-700000


#Save map with trials and provenance locations 
tiff("Fsyl_prov_trial.tif", width = 2780, height = 1680, units = "px", pointsize = 12, res = 300,  compression = 'lzw')

ggplot() + 

  theme_minimal()+  
  geom_polygon(data=wmap_laea_df, aes(long,lat,group=group), fill="white", color="black")+
  geom_polygon(data = shapefile_df_laea, aes(x = long, y = lat, group = group), color = 'dark grey', fill = 'dark grey', size = .2, alpha= 0.5)+	
  geom_point(data=datatrialsBU19_laea_df, aes(longitude, latitude, group=NULL),
            color="green3", size=2.5, alpha=0.5) +
  geom_point(data=datatrialsBU20_laea_df, aes(longitude, latitude, group=NULL),
             color="blue3", size=2.5, alpha=0.5) +	
  geom_point(data=dataprovBU19_laea_df, aes(longitude, latitude, group=NULL),
             color="green3", size=2, alpha=0.5, shape=1) +
  geom_point(data=dataprovBU20_laea_df, aes(longitude, latitude, group=NULL),
             color="blue3", size=2, alpha=0.5, shape=1) +
  coord_cartesian(xlim = c((xmin-buff),(xmax+buff)), ylim = c((ymin-buff),(ymax+buff))) +
  theme(aspect.ratio=1)

dev.off()



