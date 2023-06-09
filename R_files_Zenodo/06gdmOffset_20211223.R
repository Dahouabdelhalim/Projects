
setwd("F:/data_analysis/GradientForest20210517/environment_m5/14GDM_model_each_env/env1_outlier")

require(raster)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(fields)



require(geosphere)
#function to remove intercept from gdm predictions
removeIntercept <- function(mod,pred){
  adjust <- 0 - log(1-pred) - mod$intercept
  adjustDissim <- 1-exp(0-adjust)
  return(adjustDissim)
}

##############
#Load and prep FST data
##############
adptMat <- read.csv("Qa_22pop_lfmm_each_env_total.csv")

#Read in population locations
pops <- read.csv("Qa_22pop_info.csv")

#make sure the two are in the same order
all(adptMat$pop == pops$code)

##############
#Load and prep shapefile and climate data
##############
#load shapefile
shp <- shapefile("RasterT_Int_Rec1.shp")

#choose predictors
predNames <- c("bio1","bio2","bio3","bio4","bio8","bio13","bio15","bio19")

raster1 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10current_climate/bio1.tif")
raster2 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10current_climate/bio2.tif")
raster3 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10current_climate/bio3.tif")
raster4 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10current_climate/bio4.tif")
raster5 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10current_climate/bio8.tif")
raster6 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10current_climate/bio13.tif")
raster7 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10current_climate/bio15.tif")
raster8 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10current_climate/bio19.tif")

presClim <- stack(raster1,raster2,raster3,raster4,raster5,raster6,raster7,raster8)
presClim <- presClim[[predNames]]

#Creates pred data for gdm (cols = population name, long, lat, climate data)
pred <- data.frame(pop=pops$code,long=pops$long, lat=pops$lat, extract(presClim, y=pops[,c("long","lat")]),
                   stringsAsFactors=FALSE)

######################
#GDM model
######################
#Create site pair table
sitePair <- formatsitepair(bioDat = adptMat, bioFormat=3, siteColumn="pop", XColumn="long",YColumn="lat",predData=pred)

#Create and plot gdm
mod <- gdm(na.omit(sitePair), geo=FALSE)

#load future climate data
raster1 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10future_cilmate_RCP85/bio1.tif")
raster2 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10future_cilmate_RCP85/bio2.tif")
raster3 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10future_cilmate_RCP85/bio3.tif")
raster4 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10future_cilmate_RCP85/bio4.tif")
raster5 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10future_cilmate_RCP85/bio8.tif")
raster6 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10future_cilmate_RCP85/bio13.tif")
raster7 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10future_cilmate_RCP85/bio15.tif")
raster8 <- raster("F:/data_analysis/GradientForest20210517/environment_m5/06GDM_model_PC1_BF25/10future_cilmate_RCP85/bio19.tif")

futClims <- stack(raster1,raster2,raster3,raster4,raster5,raster6,raster7,raster8)
futClims <- futClims[[predNames]]
futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)

#Getting all coordinates in range map
popDat <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
popDat <- data.frame(distance=1, weight=1, popDat)
popDat <- split(popDat, seq(nrow(popDat)))


###############
#Forward offset calculation
##############
cl <- makeCluster(4) #ideally should be run in parallel to reduce computing time
registerDoParallel(cl)
forwardOffsetGDM <- foreach(i = 1:length(popDat), .packages=c("fields","gdm","geosphere")) %dopar%{
  
  #get the focal population
  onePop <- popDat[[i]]
  
  #set up a dataframe where the first site is the focal population, and the second population
  #are sites across North America
  setUp <- cbind(onePop,futClimDat)
  colnames(setUp) <- c("distance","weights",
                       "s1.xCoord", "s1.yCoord",paste("s1.", predNames, sep=""), 
                       "s2.xCoord", "s2.yCoord",paste("s2.", predNames, sep=""))
  
  #rearrange the colums for the gdm prediction
  dat <- setUp[,c("distance","weights","s1.xCoord", "s1.yCoord","s2.xCoord", "s2.yCoord",
                  paste("s1.", predNames, sep=""), 
                  paste("s2.", predNames, sep=""))]
  
  #do the prediction and set up a dataframe with second sites x/y and predicted Fst
  combinedDat <- predict(object=mod, dat, time=FALSE)
  combinedDat <- data.frame(dat[,c("s2.xCoord","s2.yCoord")], predFst=removeIntercept(mod, combinedDat))
  
  ##Get metrics for the focal population
  #coordinate of focal population
  coord <- onePop[,c("x","y")]
  
  #choose the pixels with the minimum fst
  minCoords <- combinedDat[which(combinedDat$predFst == min(combinedDat$predFst)),]
  
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoords["dists"] <- distGeo(p1=coord, p2=minCoords[,1:2])
  minCoords <- minCoords[which(minCoords$dists == min(minCoords$dists)),]
  
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoords <- minCoords[sample(1:nrow(minCoords),1),]
  
  #get local offset
  offset <- combinedDat[which(combinedDat$s2.xCoord == coord$x & combinedDat$s2.yCoord == coord$y),"predFst"]
  
  #get the minimum predicted fst - forward offset in this case
  minVal <- minCoords$predFst
  
  #get distance and coordinates of site that minimizes fst
  toGo <- minCoords$dists
  minPt <- minCoords[,c("s2.xCoord", "s2.yCoord")]
  
  #get bearing to the site that minimizes fst
  bear <- bearing(coord, minPt)
  
  #write out
  out <- c(x1=coord[[1]], y1=coord[[2]],local=offset,forwardFst=minVal, predDist=toGo, bearing=bear,x2=minPt[[1]],y2=minPt[[2]])
  
}

stopCluster(cl)

#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#forwardFst: forward offset
#predDist: distance to site of forward offset
#bearing: bearing to site of forward offset
#x2/y2: coordinate of site of forward offset
forwardOffsetGDM <- do.call(rbind, forwardOffsetGDM)

write.csv(forwardOffsetGDM,paste0("./Qa_each_env_lfmm_total_rcp85_forwardOffsetGDM.csv"), row.names=FALSE)


###############
#Reverse offset calculation
##############
#Getting all coordinates in the range in current climate
popDat <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
popDat <- data.frame(popDat)

#Gets climate data from the range in future climate
futClimMask <- mask(x=futClims, mask=shp)
futClimDat <- as.data.frame(futClimMask, xy=TRUE, na.rm=TRUE)

#set up for prediction
futClimDat <- data.frame(distance=1, weight=1, futClimDat)

###############
#Reverse offset calculation
##############
cl <- makeCluster(4)
registerDoParallel(cl)
reverseOffsetGDM <- foreach(i = 1:nrow(futClimDat), .packages=c("fields","gdm","geosphere")) %dopar%{
  
  #get the focal population
  onePop <- futClimDat[i,]
  
  #set up a dataframe where the first site is the focal population, and the second population
  #are sites across the range
  setUp <- cbind(onePop,popDat)
  colnames(setUp) <- c("distance","weights",
                       "s1.xCoord", "s1.yCoord",paste("s1.", predNames, sep=""), 
                       "s2.xCoord", "s2.yCoord",paste("s2.", predNames, sep=""))
  
  #rearrange the colums for the gdm prediction
  dat <- setUp[,c("distance","weights","s1.xCoord", "s1.yCoord","s2.xCoord", "s2.yCoord",
                  paste("s1.", predNames, sep=""), 
                  paste("s2.", predNames, sep=""))]
  
  #do the prediction and set up a dataframe with second sites x/y and predicted Fst
  combinedDat <- predict(object=mod, dat, time=FALSE)
  combinedDat <- data.frame(dat[,c("s2.xCoord","s2.yCoord")], predFst=removeIntercept(mod, combinedDat))
  
  ##Get metrics for the focal population
  #coordinate of focal population
  coord <- onePop[,c("x","y")]
  
  #choose the pixels with the minimum fst
  minCoords <- combinedDat[which(combinedDat$predFst == min(combinedDat$predFst)),]
  
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoords["dists"] <- distGeo(p1=coord, p2=minCoords[,1:2])
  minCoords <- minCoords[which(minCoords$dists == min(minCoords$dists)),]
  
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoords <- minCoords[sample(1:nrow(minCoords),1),]
  
  #get local offset
  offset <- combinedDat[which(combinedDat$s2.xCoord == coord$x & combinedDat$s2.yCoord == coord$y),"predFst"]
  
  #get the minimum predicted fst - reverse offset in this case
  minVal <- minCoords$predFst
  
  #get distance and coordinates of site that minimizes fst
  toGo <- minCoords$dists
  minPt <- minCoords[,c("s2.xCoord", "s2.yCoord")]
  
  #get bearing to the site that minimizes fst
  bear <- bearing(coord, minPt)
  
  #write out
  out <- c(x1=coord[[1]], y1=coord[[2]],local=offset,reverseFst=minVal, predDist=toGo, bearing=bear,x2=minPt[[1]],y2=minPt[[2]])
  
}

stopCluster(cl)

#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset - included as sanity check - should be identical to 'offset' from the calculation of forward offset above
#reverseFst: reverse offset
#predDist: distance to site of reverse offset
#bearing: bearing to site of reverse offset
#x2/y2: coordinate of site of reverse offset
reverseOffsetGDM <- do.call(rbind, reverseOffsetGDM)

write.csv(reverseOffsetGDM,paste0("./Qa_each_env_lfmm_total_rcp85_reverseOffsetGDM.csv"), row.names=FALSE)
