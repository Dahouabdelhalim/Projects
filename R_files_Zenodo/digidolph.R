#Digidolph- create spatially-explicit null association model from dolphin data
#V. Foroughirad
#vjf2@duke.edu
#Created August 8th, 2016
#Modified June 28th, 2017 

library(adehabitatHR)
library(maptools)
library(rgdal)
library(spatstat)
library(Digiroo2)
library(coda)
library(spdep)
library(raster)
library(PBSmapping)
library(gdata)
library(rgeos)

options(stringsAsFactors = FALSE)

source("digidolph_helper_functions.R")

#optional parameters to set for digidolph

#minimum number of sightings per animal

min_sightings<-45

#number of simulations to run

num_sim<-2

#add in optimum gprox to test, or derive from model

gprox_opt<-1000

#filter raw sightings to include only those which fell within an area with yearly coverage

raw_sightings<-read.csv("dryad_dolphin_data.csv")

raw_sightings$Date<-as.Date(raw_sightings$Date, format=c("%d-%b-%Y"))

raw_sightings$year<-format(raw_sightings$Date,"%Y") 

xydata<-cbind(raw_sightings$gps_east,raw_sightings$gps_south)
xydata2<-as.data.frame(project(xydata, "+proj=tmerc +lat_0=-25 +lon_0=113 +k=0.99999 +x_0=50000 +y_0=100000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
xydata3<-cbind(raw_sightings[,c("Date", "observation_id", "dolphin_id", "year")],trunc(xydata2,0))
colnames(xydata3)<-c("Date", "observation_id", "dolphin_id", "year","X","Y")

yearly_xydata<-SpatialPointsDataFrame(xydata3[,c("X","Y")],xydata3["year"])

mcps<-mcp(yearly_xydata[,1], percent=100, unin=c("m"), unout=c("m2"))

#plot(mcps, col=c(1:10), axes=TRUE)

mcpolyset<-SpatialPolygons2PolySet(mcps)
mcp_allyears<-joinPolys(mcpolyset, operation="INT")
mcp_allyears<-PolySet2SpatialPolygons(mcp_allyears)
#plot(mcp_allyears, add=TRUE, col="green")

#add buffer to intersection

buffered<-raster::buffer(mcp_allyears, width=2500)
#plot(buffered, axes=TRUE, col="yellow")

#intersect buffered with hrxydata to get surveys that fall in that range

sighting_xydata<-SpatialPointsDataFrame(xydata3[,c("X","Y")],xydata3["observation_id"])

buff_surveys<-intersect(sighting_xydata, buffered)
#plot(buff_surveys, add=TRUE)

buff_surveys<-as.data.frame(buff_surveys)

fs<-subset(xydata3, xydata3$observation_id %in% buff_surveys$observation_id)

#calculate number of sightings for each individual in set

fs$sightings<-sapply(fs$dolphin_id, function(i) length(fs$dolphin_id[which(fs$dolphin_id==i)]))

#filter based on number of sightings

fs45<-subset(fs, fs$sightings>=min_sightings)

#Calculate the real HWI for this data set

availability <- read.csv("dryad_dolphin_availability.csv")
availability$entry<-as.Date(availability$entry, format=c("%d-%b-%Y"))
availability$depart<-as.Date(availability$depart, format=c("%d-%b-%Y"))

realHWI<-hwi_filtered(sightings=fs45, group_variable="observation_id", dates="Date", IDs="dolphin_id", symmetric=FALSE, availability=availability)

#Construct home ranges for each animal

#Create grid with 5km buffer

grid_buffer=5000

x <- seq(min(fs45[,"X"])-grid_buffer,max(fs45[,"X"])+grid_buffer,by=100) # where resolution is the pixel size you desire
y <- seq(min(fs45[,"Y"])-grid_buffer,max(fs45[,"Y"])+grid_buffer,by=100)

xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

#Create UDs for each animal and extract h values (need to manually select h for boundary method)

hrxydata<-SpatialPointsDataFrame(fs45[,c("X","Y")],fs45["dolphin_id"])

uds_href<-kernelUD(hrxydata[,1],grid=xy)

hvalues<-list()
for (i in 1:length(uds_href)) {
  h<-uds_href[[i]]@h$h
  id<-names(uds_href)[[i]]
  hvalues[[i]]<-c(h, id)
}

h<-as.data.frame(do.call("rbind", hvalues), stringsAsFactors=FALSE)

names(h)<-c("h_opt", "dolphin_id")

h$h_opt<-as.numeric(h$h_opt)

#Create simplified coastline, length of segments must be greater than 3*h
bound <- structure(list(x = c(122000,122000,116500,110000,108000), y = c(1000,10500,14500,20800,31280)), .Names = c("x", "y"))
bound <- do.call("cbind",bound)
Slo1 <- Line(bound)
Sli1 <- Lines(list(Slo1), ID="frontier1")
barrier <- SpatialLines(list(Sli1))

optud<-list()

#Filter out unnecesary parts from loop, extract polygons by ID

for (i in 1:dim(h)[1]){

  cdol<-hrxydata[hrxydata$dolphin_id==h$dolphin_id[i],]
  hopt<-h$h_opt[i]
  uds_man<-kernelUD(cdol,h=hopt,grid=xy, boundary=barrier)
  optud[[i]]<-uds_man
  cat(i)
}

uddf<-unlist(optud)
class(uddf)<-"estUDm"


#remove land from new estimates

coast_polygon<-readOGR("coastpolygon", "coastpolygon")

coast_polygon<-spTransform(coast_polygon, CRS("+proj=tmerc +lat_0=-25 +lon_0=113 +k=0.99999 +x_0=50000 +y_0=100000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

#plot(coast, axes=TRUE, lwd=2);plot(kernelcontours90, add=TRUE, col="blue");plot(coast_polygon, add=TRUE, lwd=2)

udsgdf <- as(estUDm2spixdf(uddf),"SpatialGridDataFrame")

#use coast polygon as mask for spatial grid

rgrid <- raster(udsgdf)
#plot(rgrid)
rgrid_msk <- mask(rgrid,coast_polygon, inverse=TRUE)
#plot(rgrid_msk)

#convert back to sgdf

grid_ae <- as(rgrid_msk, 'SpatialGridDataFrame')
gridded(grid_ae) <- TRUE
grid_ae[[1]] <- as.numeric(!is.na(grid_ae[[1]])) 

## Then, you just have to multiply each column of udsgdf by the mask 
resu <- lapply(1:ncol(udsgdf), function(i) { udsgdf[[i]] * grid_ae[[1]] }) 
#Alternatively, you can re-standardize

resu <- lapply(1:ncol(udsgdf), function(i) {udsgdf[[i]] * grid_ae[[1]] / sum(udsgdf[[i]] * grid_ae[[1]]) })

resu <- as.data.frame(resu) 
names(resu) <- names(udsgdf@data) 

## and define it as data slot for udsgdf 
udsgdf@data <- resu 

#Will also need to do this for contours using the masked grid we just creasted 

masked_grid<-udsgdf

fullgrid(masked_grid) <- FALSE 
re <- lapply(1:ncol(masked_grid), function(i) { 
  so <- new("estUD", masked_grid[,i]) 
  so@h <- list(h=0, meth="specified") # fake value 
  so@vol <- FALSE  
  return(so) 
}) 
names(re) <- names(masked_grid) 
class(re) <- "estUDm" 

#Create daily minimum convex polygons with 1km buffers, always including launch point as one of the vertices.

groupings<-unique(fs45[,c("Date", "observation_id", "X", "Y")])

days<-split(groupings, groupings$Date)

#Add a few points near to launch area to create mcp

launch<-c("2001-01-01","launch",122241,11966)

launch<-as.data.frame(t(launch))

names(launch)<-names(days[[1]])    

launch[2,]<-c("2001-01-01","launch",122241,11965)
launch[3,]<-c("2001-01-01","launch",122241,11964)
launch[4,]<-c("2001-01-01","launch",122241,11963)

launch[,1]<-as.Date(launch[,1])

days1<-lapply(days, function(x) add_launch(x))

survey_days<-as.data.frame(do.call("rbind", days1))

survey_days[,c("X","Y")]<-apply(survey_days[,c("X","Y")],2, as.numeric)

daily_xydata<-SpatialPointsDataFrame(survey_days[,c("X","Y")],survey_days["Date"])

mcps<-mcp(daily_xydata[,1], percent=100, unin=c("m"), unout=c("m2"))

#Add buffer, make sure whole area is covered

buff_days<-gBuffer(mcps, byid=TRUE,width=1000)

#Number of animals in study
n<-length(unique(fs45$dolphin_id))

#Number of survey days
d<-length(unique(fs45$Date))
dates<-unique(fs45$Date)

#Get availability matrix

dolphins<-sort(unique(fs45$dolphin_id))

matnames<-list(dates,dolphins)

alive<-Vectorize(FUN=function(r,c) isTRUE(r>=availability$entry[which(availability$dolphin_id==c)] & r<=availability$depart[which(availability$dolphin_id==c)]))

schedule<-outer(dates, dolphins, FUN=alive)
dimnames(schedule)<-matnames

dolphin_density_per_km<-dim(fs45)[1]/(sum(area(buff_days))/1000000)

#################Ok now set up gprox testing?
Gprox_set<-sort(rep(seq(600, 1400, by=100),10))

digidolph_gprox<-list()

for (i in 1:length(Gprox_set)){
daily_sim<- one_sim(d=d,
            buff_days=buff_days, 
            udsgdf=udsgdf, 
            dolphin_density_per_km=dolphin_density_per_km, 
            schedule=schedule, 
            gprox=Gprox_set[i])

Assoctable<-do.call("rbind",daily_sim)

Assoctable<-digu(Assoctable)

digidolph_gprox[[i]]<-Assoctable
}

rand_mats<-lapply(digidolph_gprox,function(x) hwi_filtered(sightings=x,group_variable="Group", dates="Permutation", IDs="IDs", symmetric=TRUE, availability=availability))
rm<-mergeMatrices(rand_mats)

#check output to determine optimum gprox

#set optimum gprox

gprox<-gprox_opt

digidolph<-list()

for (k in 1:num_sim){
  daily_sim<- one_sim(d=d,
                      buff_days=buff_days, 
                      udsgdf=udsgdf, 
                      dolphin_density_per_km=dolphin_density_per_km, 
                      schedule=schedule, 
                      gprox=gprox)
  
  Assoctable<-do.call("rbind",daily_sim)
  
  Assoctable<-digu(Assoctable)
  
  digidolph[[k]]<-Assoctable
}

rand_mats<-lapply(digidolph,function(x) hwi_filtered(sightings=x,group_variable="Group", dates="Permutation", IDs="IDs", symmetric=TRUE, availability=availability))
rm<-mergeMatrices(rand_mats)

realh<-as.data.frame(na.omit(unmatrix(realHWI)))
names(realh)<-"realHWI"

#Merge real association index with random values

output<-merge(rm, realh, by="row.names") 




