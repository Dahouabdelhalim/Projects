##creating a spatially explicit null model for eastern water dragon data 
#K.Strickland

library(adehabitatHR)
library(maptools)
library(rgdal)
library(spatstat)
#library(devtools) # only need for install
#install_git("git://github.com/gsk3/taRifx.geo.git")
library(taRifx.geo)
library(Digiroo2)
library(coda)
library(spdep)
library(gtools)
library(gdata)
library(raster)

source("hwi_functions.R")

#### read in data and create home ranges; subset to minimum sightings required to create stable home range ##

HRdata<- read.csv("dragon_data2.csv",stringsAsFactors = FALSE)
fs25<-subset(HRdata,HRdata$Sightings>=25)
xydata<-cbind(fs25$X,fs25$Y)
xydata2<-as.data.frame(project(xydata, "+proj=tmerc +lat_0=-28 +lon_0=153 +k=0.99999 +x_0=50000 +y_0=100000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
xydata3<-cbind(fs25$Name,xydata2)
colnames(xydata3)<-c("Name","X","Y")

x <- seq(min(xydata3[,"X"])-20,max(xydata3[,"X"])+20,by=1) # where resolution is the pixel size you desire
y <- seq(min(xydata3[,"Y"])-20,max(xydata3[,"Y"])+20,by=1)
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE
class(xy)

hrxydata<-SpatialPointsDataFrame(xydata3[,2:3],xydata3["Name"])
uds7<-kernelUD(hrxydata[,1],h= 7,grid=xy) 

udsgdf <- as(estUDm2spixdf(uds7),"SpatialGridDataFrame")

## real half-weights

realHWI<-hwi(fs25, group_variable="observation_id", dates="Date", IDs="Name")
rHWI<-unmatrix(realHWI)

## create probability tables for use in simulations

#Number of animals in study
n<-length(unique(fs25$Name))

#Number of survey days
d<-length(unique(fs25$Date))

ID_counts <- tapply(rep(1,length(fs25$Name)),fs25$Name,sum)
No_SurveyDays <- d
pObs <- as.vector(ID_counts)/No_SurveyDays # proportion of time individuals are sighted in surveys
nameObs <- names(ID_counts)
TotalAnimals <- n # Total number of animals in study
samplesize <- round(dim(fs25)/d) # Number of IDs to include in simulations (with replacement) - average number of ids per survey in the observed data
bootlength <- d # Number of permutations (equal to number of survey days)
ExpProb <- data.frame(Kangaroo=nameObs,Probability=pObs) # Table of Expected Probabilities

########### RUN model to assign the gprox value

# set distances you wish to test

Gprox_full<-sort(rep(seq(1,10,by=1),10))

fullgprox<-list()

for (k in 1:length(Gprox_full)) {
  
  Gprox<-Gprox_full[[k]]
  pID <- sapply(1:bootlength,function(i) sample(x=TotalAnimals,size=samplesize,replace=FALSE,prob=pObs))
  result_full<- fAssocmatrix(sPerm=1:bootlength,Gprox=Gprox,iextract=udsgdf,iID=pID)
  fullgprox[[k]]<-result_full
  
  #verbose
  cat(k)
  
}

full<-lapply(fullgprox, digu)
rand_mats<-lapply(full,function(x) hwi(sightings=x,group_variable="Group", dates="Permutation", IDs="IDs",symmetric = TRUE))
ran_m<-mergeMatrices(rand_mats)
rm<-lapply(ran_m,function(x) mean(x))

rHWI<-na.exclude(rHWI)
ob_mean<-mean(rHWI)

##check output to determine optimum gprox and use the gprox value closest to observed mean (ob_mean)

## set optimum gprox (e.g. 5m for dragons)

Gprox <- 5

########### RUN FULL PERMUTATIONS


full<-list()

 for (k in 1:1000) {
  
  pID <- sapply(1:bootlength,function(i) sample(x=TotalAnimals,size=samplesize,replace=FALSE,prob=pObs))
  result_full<- fAssocmatrix(sPerm=1:bootlength,Gprox=Gprox,iextract=udsgdf,iID=pID)
  full[[k]]<-result_full
  #verbose
  cat(k)
  
}

full<-lapply(full, digu)
rand_mats<-lapply(full,function(x) hwi(sightings=x,group_variable="Group", dates="Permutation", IDs="IDs",symmetric = TRUE))
full_mats<-list()

for (i in 1:length(rand_mats)){
  
  rand_mats[[i]]->mat
  unmat<-data.frame(unmatrix(mat))
  IDs<-row.names(unmat)
  unmat<-data.frame(cbind(IDs,unmat$unmatrix.mat.))
  colnames(unmat)<-c("ID","HWI")
  mat_f<-unmat[order(unmat$ID),]
  full_mats[[i]]<-mat_f
  
  #verbose
  
  cat(i)
  
}

final<-do.call("cbind",full_mats)
finalf<-final[ , -which(names(final) %in% c("ID"))]
final_f<-cbind(final[,1],finalf)
rownames(final_f)<-final_f$`final[, 1]`
#Merge real association index with random values

output<-merge(final_f, rHWI, by="row.names") 
