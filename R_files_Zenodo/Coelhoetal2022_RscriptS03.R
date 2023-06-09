#Install packages
install.packages("ENMeval")
install.packages("spocc")
install.packages("raster")
install.packages("rgdal")
install.packages("maptools")
install.packages("dplyr")
install.packages("ecospat")

##Load packages 
library(ENMeval)
library(spocc)
library(raster)
library(rgdal)
library(maptools)
library(dplyr)
library(ecospat)

#Read files
files<-stack(list.files(path = "PathTo/var",pattern='bil', full.names=T)) ### import Bioclim variables
e<-extent(-78.352, -34.099, -20.812, 12.839) #Delimit space to South America
files<-crop(files,e) #Crop to South America
plot(files[[1]])
dev.off()

projection(files)<- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # coodinate reference system for WGS84
plot(files[[1]])

####remove highly correlated variables 
files.sub<- dropLayer(files, c(1,2,5,8,9,10,11,13,16,17,18,19)) #### remove selected layers

names(files.sub)
plot(files.sub[[2]])

predictors<-files.sub


#### AMF ####
#Import occurrence points from Amazonia
amf<-read.csv("PathTo/occ_thin/amf.csv",header = T)

bg_bio_ls<-list()
occ_bio_ls<-list()

  #Import background shapefile (ecoregions where occurrences can be found)
  sp_shape_amf <- readOGR(paste('PathTo/background/amf.shp',sep=""))
  climatelayers <- predictors
  #Extract background climatic data
  sp_Bg_amf<- mask(climatelayers, sp_shape_amf)
  sp_Bg_xy_amf<-rasterToPoints(sp_Bg_amf)
  sp_Bg_xy_amf<-sp_Bg_xy_amf[complete.cases(sp_Bg_xy_amf),]


  species_amf<-data.frame(rep("amf",nrow(sp_Bg_xy_amf)))
  colnames(species_amf)<-"sp"
  sp_Bg_xy_amf<-cbind(species_amf,sp_Bg_xy_amf)
  sp_Bg_xy_amf<-rename(sp_Bg_xy_amf,lon=x,lat=y)
  bg_bio_ls_amf<-sp_Bg_xy_amf

  #Extract occurrence points climatic data
  occ.file_amf <- amf[which(amf[,1]==paste(unique(amf$sp))),]
  occ_amf <- occ.file_amf[,2:3]
  sp_occ_bg_amf<-extract(climatelayers, occ_amf)

  bioclim_sp_amf<- cbind(occ.file_amf,sp_occ_bg_amf)
  bioclim_sp_amf<-bioclim_sp_amf[complete.cases(sp_occ_bg_amf),]
  occ_bio_ls_amf<-bioclim_sp_amf

#Combine both datasets
bioclim_bg_amf<-rbind(bg_bio_ls_amf)
bioclim_sp_amf<-rbind(occ_bio_ls_amf)

bioclimdata_amf<-rbind(bioclim_bg_amf,bioclim_sp_amf)
sp<-bioclimdata_amf[,1]

bioclimdata_amf$obs_type<-rep(c("background","occ"),t=c(nrow(bioclim_bg_amf),nrow(bioclim_sp_amf)))
dim(bioclimdata_amf)

bioclimdata_amf<-bioclimdata_amf[complete.cases(bioclimdata_amf),]

data.frame(table(bioclimdata_amf[bioclimdata_amf$obs_type == "background",1]))
data.frame(table(bioclimdata_amf[bioclimdata_amf$obs_type == "occ",1]))




#### ATF ####
#Import occurrence points from Atlantic forest
atf<-read.csv("08ENMeval/occ_thin/atf.csv",header = T)

bg_bio_ls_atf<-list()
occ_bio_ls_atf<-list()

  #Import background shapefile (ecoregions where occurrences can be found)
  sp_shape_atf <- readOGR(paste('PathTo/background/atf.shp',sep=""))
  climatelayers <- predictors
  sp_Bg_atf<- mask(climatelayers, sp_shape_atf)
  sp_Bg_xy_atf<-rasterToPoints(sp_Bg_atf)
  sp_Bg_xy_atf<-sp_Bg_xy_atf[complete.cases(sp_Bg_xy_atf),]

  species_atf<-data.frame(rep("atf",nrow(sp_Bg_xy_atf)))
  colnames(species_atf)<-"sp"
  sp_Bg_xy_atf<-cbind(species_atf,sp_Bg_xy_atf)
  sp_Bg_xy_atf<-rename(sp_Bg_xy_atf,lon=x,lat=y)
  bg_bio_ls_atf<-sp_Bg_xy_atf

  #Extract occurrence points climatic data
  occ.file_atf <- atf[which(atf[,1]==paste(unique(atf$sp))),]
  occ_atf <- occ.file_atf[,2:3]
  sp_occ_bg_atf<-extract(climatelayers, occ_atf)

  bioclim_sp_atf<- cbind(occ.file_atf,sp_occ_bg_atf)
  bioclim_sp_atf<-bioclim_sp_atf[complete.cases(sp_occ_bg_atf),]
  occ_bio_ls_atf<-bioclim_sp_atf

#Combine both datasets
bioclim_bg_atf<-rbind(bg_bio_ls_atf)
bioclim_sp_atf<-rbind(occ_bio_ls_atf)

bioclimdata_atf<-rbind(bioclim_bg_atf,bioclim_sp_atf)
sp<-bioclimdata_atf[,1]

bioclimdata_atf$obs_type<-rep(c("background","occ"),t=c(nrow(bioclim_bg_atf),nrow(bioclim_sp_atf)))
dim(bioclimdata_atf)

bioclimdata_atf<-bioclimdata_atf[complete.cases(bioclimdata_atf),]

data.frame(table(bioclimdata_atf[bioclimdata_atf$obs_type == "background",1]))
data.frame(table(bioclimdata_atf[bioclimdata_atf$obs_type == "occ",1]))


#### Merge Amazonia and Atlantic forest climatic information into one csv file ####
bioclimdata_atf
bioclimdata_amf

bioclimdata_all<-rbind(bioclimdata_atf,bioclimdata_amf)

write.csv(bioclimdata_all,"08ENMeval/broennimann/bioclimdata_all.csv",row.names = F)