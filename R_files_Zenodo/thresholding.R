
library(maptools)
library(raster)
WI <- readShapeSpatial("~/Wisconsin_tz.shp", verbose=T,proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
names <- unique(species_names)
thresholds <- as.data.frame(matrix(nrow=length(names),ncol=2))
colnames(thresholds) <- c("species","threshold")


now <- paste("~/",names[1],"_bioclim.asc",sep="")
nowraster <- raster(now)
resampled <- raster(nrow=225,ncol=726/2) 
extent(resampled) <- extent(nowraster)
nowraster <- resample(nowraster,resampled,method="bilinear")
for(j in 1:10) {
  quantile_test[j,3] <- quantile(nowraster,quantile_test[j,2])/quantile(nowraster,quantile_test[j,1])
}
print(quantile_test)
min <- quantile_test[which.min(quantile_test$V3),]
rname <- as.numeric(rownames(min))
if(rname<2){
  threshold <- quantile(nowraster,0.3)
  replace_with_one <-  which(values(nowraster)>=threshold)
  replace_with_zero <- which(values(nowraster)<threshold)
  values(nowraster)[replace_with_one]<-1
  values(nowraster)[replace_with_zero]<-0
} else{
  threshold <- quantile(nowraster,min[,2])
  replace_with_one <-  which(values(nowraster)>=threshold)
  replace_with_zero <- which(values(nowraster)<threshold)
  values(nowraster)[replace_with_one]<-1
  values(nowraster)[replace_with_zero]<-0
}
thresholds[1,1] <- names[1]
thresholds[1,2] <- threshold
now_brick <- stack(nowraster)


for (i in 11:length(names)){
  print(i)
  now <- paste("~/",names[i],"_bioclim.asc",sep="")
  nowraster <- raster(now)
  nowraster <- resample(nowraster,resampled,method="bilinear")
  pointsfile <- paste("~/",names[i],"_samplePredictions.csv",sep="")
  data <- read.csv(pointsfile, header=T)
  SpatialPoints(data[,1:2],proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs")) -> point_shp
  SpatialPointsDataFrame(coords=point_shp,data=data)->pointsSPDF
  sub <- pointsSPDF[WI,]
  features <- nrow(sub@data)
  if(features<1){
    replace_with_zero <- which(values(nowraster)<1)
    values(nowraster)[replace_with_zero] <- 0
    thresholds_climate[i,1] <- names[i]
    thresholds_climate[i,2] <- 0
  } 
  if(features>=0 & features<11){
    threshold <- quantile(nowraster,0.9)
    replace_with_one <-  which(values(nowraster)>=threshold)
    replace_with_zero <- which(values(nowraster)<threshold)
    values(nowraster)[replace_with_one]<-1
    values(nowraster)[replace_with_zero]<-0
  }
  if(features>=11){
    for(j in 1:10) {
      quantile_test[j,3] <- quantile(nowraster,quantile_test[j,2])/quantile(nowraster,quantile_test[j,1])
    }
    min <- quantile_test[which.min(quantile_test$V3),]
    rname <- as.numeric(rownames(min))
    if(rname<7){
      threshold <- mean(na.omit(values(nowraster)))-0.3*mean(na.omit(values(nowraster)))
      replace_with_one <-  which(values(nowraster)>=threshold)
      replace_with_zero <- which(values(nowraster)<threshold)
      values(nowraster)[replace_with_one]<-1
      values(nowraster)[replace_with_zero]<-0
    }
    if(rname>-2 & rname<10){
      threshold <- quantile(nowraster,min[,1])
      replace_with_one <-  which(values(nowraster)>=threshold)
      replace_with_zero <- which(values(nowraster)<threshold)
      values(nowraster)[replace_with_one]<-1
      values(nowraster)[replace_with_zero]<-0
    }
    if(rname>=7){
      threshold <- quantile(nowraster,min[,2])
      replace_with_one <-  which(values(nowraster)>=threshold)
      replace_with_zero <- which(values(nowraster)<threshold)
      values(nowraster)[replace_with_one]<-1
      values(nowraster)[replace_with_zero]<-0
    }
  }
  thresholds[i,1] <- names[i]
  thresholds[i,2] <- threshold
  now_brick <- stack(now_brick,nowraster, quick=T)
}
