#This script produces maps of tree height averaged by trial or provenance

#required libraries
library(rgdal)
library(rworldmap)
library(ggmap)
library(ggplot2)
library(lattice)
data(coastsCoarse)
data(countriesLow)
library(plyr)
library(sp) 

#set working directory 
setwd("/WorkingDirectory/")


#Read files
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


#save individual beech data with trial and provenance coordinates: 

write.csv(Fsyl, file = "Fsylvatica_coord.csv", sep="\\t", col.names=TRUE)


############################################################


##Example of plot of tree height for a given age####

#Select trees of a given age (year 7 from plantation time)
T_Y <-  subset(data, ( H02cm != "NA" & H02cmAge == "7"))
unique(T_Y$H02cmAge )
T_Y2 <-  subset(data,(H05cm != "NA" & H05cmAge == "7"))
unique(T_Y2$H05cmAge)

T_Y$H7y <- T_Y$H02cm
T_Y2$H7y <- T_Y2$H05cm

year7<- merge(T_Y,T_Y2, all=T)


#Height at 7 year old trees averaged by trial
aggYear7 <- aggregate(H7y~Trial, data=year7, FUN=function(x) c(mean=mean(x)))
colnames(aggYear7) <- c("Trial","H7ymean_trial")
meanH7_by_trial<- merge(aggYear7,year7,by="Trial")
meanY7<-  meanH7_by_trial[!duplicated(meanH7_by_trial$H7ymean_trial),]
meanY7$H7ymean_trial #18 trials 


#Reproject to lambert azimutal equal-area projection (laea)

coordinates(meanY7)<-c("lonTrial","latTrial")
proj4string(meanY7) <- CRS("+proj=longlat")
year7_laea<-spTransform(meanY7, CRS("+proj=laea"))
year7_laea_df<-data.frame(year7_laea)

countries_df<-fortify(countriesLow)
year7_df<-data.frame(meanY7)

#Reproject to laea 
wmap_laea<- spTransform(countriesLow, CRS("+proj=laea")) 
wmap_laea_df<-fortify(wmap_laea)

#Read distribution range from Euforgen
shapefile <- readOGR(dsn="./", "Fagus_sylvatica_EUFORGEN")
shapefile_df <- fortify(shapefile)
shapefile_laea<- spTransform(shapefile, CRS("+proj=laea")) 
shapefile_df_laea <- fortify(shapefile_laea)


#Extract georaphical tile from trial coordinates
xmin<-min(year7_laea_df$lonTrial)
xmax<-max(year7_laea_df$lonTrial)
ymin<-min(year7_laea_df$latTrial)
ymax<-max(year7_laea_df$latTrial)
xmin<-min(year7_laea_df$lonTrial)
xmax<-max(year7_laea_df$lonTrial)
ymin<-min(year7_laea_df$latTrial)
ymax<-max(year7_laea_df$latTrial)
buff<-900000 

#plot tree height by trial 
tiff("Height_7yr_by_trial_colors.tif", width = 2180, height = 1480, units = "px", pointsize = 12, res = 300,  compression = 'lzw')

ggplot() + 
  
  theme_minimal()+  	
  geom_polygon(data=wmap_laea_df, aes(long,lat,group=group), fill="white", color="black")+
  geom_polygon(data = shapefile_df_laea, aes(x = long, y = lat, group = group), color = 'Dark Grey', fill = 'Dark Grey', size = .2, alpha= 0.5)+		
  geom_point(data=year7_laea_df, aes(lonTrial, latTrial, group=NULL,fill=NULL, size=H7ymean_trial,  colour=H7ymean_trial), alpha=0.8) +
  scale_color_gradient(low="red", high='blue', limits=c(0,400),  breaks=c(100,200,300,400)) +
  guides(color= guide_legend(), size=guide_legend()) +
  scale_size_continuous(limits=c(0,400), breaks=c(100,200,300,400)) + 
  coord_cartesian(xlim = c((xmin-buff),(xmax+buff)), ylim = c((ymin-buff),(ymax+buff))) +
  theme(aspect.ratio=1)

dev.off()


#Height at 7 year old trees averaged by provenance
aggYear7_P <- aggregate(H7y~ID_ProvCode, data=year7, FUN=function(x) c(mean=mean(x)))
colnames(aggYear7_P) <- c("ID_ProvCode","H7ymean_prov")
meanH7_by_prov<- merge(aggYear7_P,year7,by="ID_ProvCode")
meanY7_P<-  meanH7_by_prov[!duplicated(meanH7_by_prov$H7ymean_prov),]
meanY7_P$H7ymean_prov #188 provenances with H measures at 7 years old 
meanY7_P<-  subset(meanY7_P,(latProv != "NA" ))


#Reproject to laea  
coordinates(meanY7_P)<-c("lonProv","latProv")
coordinates(meanY7_P)<-c("lonProv","latProv")
proj4string(meanY7_P) <- CRS("+proj=longlat")
year7_P_laea<-spTransform(meanY7_P, CRS("+proj=laea"))

countries_df<-fortify(countriesLow)
year7_P_df<-data.frame(meanY7_P)

wmap_laea_df<-fortify(wmap_laea)
year7_P_laea_df<-data.frame(year7_P_laea)

#plot tree height by provenance 
tiff("Height_7yr_by_prov_color.tif", width = 2180, height = 1480, units = "px", pointsize = 12, res = 300,  compression = 'lzw')

ggplot() + 
 
  theme_minimal()+
  geom_polygon(data=wmap_laea_df, aes(long,lat,group=group), fill="white", color="black")+
  geom_polygon(data = shapefile_df_laea, aes(x = long, y = lat, group = group), color = 'Dark Grey', fill = 'Dark Grey', size = .2, alpha= 0.5)+	
  geom_point(data=year7_P_laea_df, aes(lonProv, latProv,     group=NULL,size=H7ymean_prov,colour=H7ymean_prov), alpha=I(6/10)) +
  scale_color_gradient(low="red", high='blue', limits=c(0,400), breaks=c(100,200,300,400)) + 	
  guides(color= guide_legend(), size=guide_legend()) +
  scale_size_continuous(limits=c(0,400), breaks=c(100,200,300,400)) + 
  coord_cartesian(xlim = c((xmin-buff),(xmax+buff)), ylim = c((ymin-buff),(ymax+buff))) +
  theme(aspect.ratio=1)

dev.off()

