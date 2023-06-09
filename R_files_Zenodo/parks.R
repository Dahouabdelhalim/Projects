library(sp)
library(spdep)
library(INLA)
library(raster)
library(ggplot2)
library(SpatialEpi)
library(maptools)
library(spatialEco)
library(rgdal)
library(reshape2)
library(plyr)
library(ROCR)
library(OptimalCutpoints)
library(viridis)
library(ggplot2)
library(ggthemes)

#requires the file "fires.Rdata" in the same folder in additon to the shapefiles below

#read in the shapefiles with fires data
fires2017<-readOGR("fires_year_2017.shp", layer="fires_year_2017")
fires2018<-readOGR("fires_year_2018.shp", layer="fires_year_2018")
parks<-readOGR("MacTinPica.shp", layer="MacTinPica")

#generate a raster for downstream distance calculation
rast <- raster(ncol = 516, nrow = 280)
extent(rast) <- extent(parks)

#rasterize the parks shapefile
pras<-rasterize(parks, rast, "hectareas_")
crs(pras)<-crs(fires2017)

#expand the parks to make a clip of the fire data
e <- extent(pras)+2

#generate a clip raster
clipr <- raster(e,nrow=280,ncol=516)
clipr[]<-runif(280*516)
clipp <- as(clipr, 'SpatialPolygons')
crs(clipp)<-crs(fires2017)

#crop the fire points to speed calculations
cfires2017 <- crop(fires2017,  clipp)
cfires2018 <- crop(fires2018,  clipp)

#calculate distance (this takes a while)
d2f2017<-distanceFromPoints(pras, cfires2017)
d2f2018<-distanceFromPoints(pras, cfires2018)

#bring data to same projection
parks<-spTransform(parks, crs(d2f2017))

#import the raster with forest cover data for 2016
cover <- raster("3PA_Forest2016.tif")

#get rid of noise
cover[cover == 128] <- NA

#resample raster to the resolution of other data
cover <- resample(cover, d2f2017, method='ngb')

#evaluate only the parks
d2f2017<-mask(d2f2017, parks)
d2f2017<-mask(d2f2017, parks)
cover<-mask(cover, parks)

#convert to data frame
d2017<-as.data.frame(d2f2017, xy=T)
d2018<-as.data.frame(d2f2018, xy=T)
dcover<-as.data.frame(cover, xy=T)

#take the Guaviare modeling data
load("fires.Rdata")

#import the raster with forest cover data for 2017
loss17 <- raster("3PA_Hansen_v1.5_loss2017.tif")

#eliminate the issue of high values introducing noise
loss17[loss17 == 128] <- NA

#eliminate previously converted pixels
loss17[loss17 == 1] <- NA
loss17[loss17 == 17] <- 1

#resample raster to the resolution of other data
ploss17 <- resample(loss17, d2f2017, method='bilinear')

#evaluate only the parks
ploss17<-mask(ploss17, parks)

#convert to data frame
dloss17<-as.data.frame(ploss17, xy=T)

#generate a dataframe for prediction
npar<-d2017
npar$y17<-d2017$layer
npar$y18<-d2018$layer
npar$y17<-log10(npar$y17/1000+1)
npar$y18<-log10(npar$y18/1000+1)
npar$loss17<-dloss17$layer
npar$cover<-dcover$layer

#keep only forest pixels
npar<-subset(npar, cover == 2, select = 1:6)

#there are no matches in geographic units in the training set
#here we randomly assign the labels from the largest veredas to the new data 
lab<-unique(defo$lab)
npar$lab<-sample(lab, dim(npar)[1], replace=T)
npar$layer<-NULL
npar<-na.omit(npar)
nare<-melt(npar[,c(1:4,6:6)], id=c("x","y","lab"), value.name="difi", variable.name="year")

#number of tries per pixel is 1
nare$dum<-1

#make columns match exactly
defo$dist<-NULL
defo$label<-NULL

#make a data set that includes observations and predictions
new<-as.data.frame(rbind.fill(defo, nare))

#this link tells inla which rows to predict as opposed to estimate
link <- rep(NA, dim(new)[1])
link[which(is.na(new$defo))] <- 1

#run the model (this is how predictions are obtained in inla)
mnew = inla(f1, data=new, family="binomial", Ntrials=dum, control.inla = list(strategy = "gaussian", int.strategy="ccd"), control.predictor=list(link = link), control.compute=list(waic=TRUE))

#write prediction for data without observations
nare$mnewp<-mnew$summary.fitted.values$mean[(dim(defo)[1]+1):dim(new)[1]]

#plot predictions
box<-ggplot(aes(y = mnewp, x = year, fill = year), data = nare) + geom_boxplot()+labs(x="Year", y="Expected probability of forest loss")+ guides(fill=guide_legend(title = "Year"))+theme_pander()+scale_fill_manual(values=c("goldenrod1", "firebrick"))
ggsave("predicted_p.png", h=4, w=4)

#to determine cutpoint and auc, pull out the estimates that match observations
new<-defo
t<-unlist(strsplit(as.character(defo$year) , "y"))
t<-t[c(FALSE, TRUE)]
t<-as.numeric(t)
t<-t+2000
new$year<-t

#Guaviare is training data
#axis is test data
#make a data frame of observation/prediction data for validation (Guaviare) data
pred<-as.data.frame(cbind(new[,c(1:5)], mnew$summary.fitted.values$mean[1:dim(defo)[1]]))
colnames(pred)[1:6]<-c("x", "y","year", "defo","difi", "pnew")
prednew.rocr<- prediction(pred$pnew,pred$defo)
prednew.auc<-performance(prednew.rocr, measure = "auc")

#sample data set for optimal cutpoint
train<-pred[sample(dim(pred)[1],dim(pred)[1]/2),]
ocnew <- optimal.cutpoints(X = pnew ~ defo, tag.healthy = 0, methods = "MaxSp", data = train, pop.prev = NULL, control = control.cutpoints(valueDLR.Positive=15), ci.fit = FALSE, conf.level = 0.95, trace = F)

#make a data frame of observation/prediction data for Hansen axis data
#note this is fraction as a result of resampling to lower resolution
pHan<-as.data.frame(cbind(nare[1:dim(npar)[1],c(1:2, 5:5, 7:7)], npar$loss17))
colnames(pHan)[1:5]<-c("x", "y","difi", "pnew","defo")

#use arbitrary cutpoint to cutoff the axis Hansen observations (not a probabilistic model)
#using the Guaviare cutpoint leads to overestimate of size of deforestation but same AUC
pHan$def<-as.ordered(ifelse(pHan$defo >.5,1,0))

#get AUC for test data
pHan.rocr<- prediction(pHan$pnew,pHan$def)
pHan.auc<-performance(pHan.rocr, measure = "auc")

#separate the years into columns
dnare<-dcast(nare, x+y ~ year , value.var="mnewp")

#calculate change in probability
dnare$change<-(dnare$y18-dnare$y17)*100/dnare$y17

#apply a cutoff for prediction each year
dnare$bin17<-as.ordered(ifelse(dnare$y17>summary(ocnew)[[5]]$Global$MaxSp[[1]][1], 1, 0))
dnare$bin18<-as.ordered(ifelse(dnare$y18 > summary(ocnew)[[5]]$Global$MaxSp[[1]][1], 1, 0))

#print out the summary of predictions
sink("summary_prediction.txt")
print(summary(dnare))
sink()

#print out the summary of Hansen
sink("Hansen.txt")
print(summary(pHan))
sink()

#print out the AUC of the test data
sink("AUC_Meta_data2017.txt")
print(pHan.auc)
sink()

#make the data frame geographic
coordinates(dnare) = ~ x+y

#make raster to plot predictions
rast <- raster(ncol = length(unique(coordinates(dnare)[,1])), nrow = length(unique(coordinates(dnare)[,2])))

#make raster same size as spatial data frame
extent(rast) <- extent(dnare)

#rasterize predictions 
delta1718<-rasterize(dnare, rast, "change", fun = median)

tiff(file = "change.tiff", width = 2400, height = 1600, units = "px", res = 300)
plot(delta1718, col=rev(inferno(256)))
dev.off()

#plot the Hansen data set
tiff(file = "Hansen2017.tiff", width = 2400, height = 1600, units = "px", res = 300)
plot(loss17, col=rev(viridis(256)))
dev.off()

rm(list=setdiff(ls(), c("cfires2017","cfires2018", "clip", "clipp","clipr","d2017", "d2018", "d2f2017", "d2f2018", "defo", "delta1718", "dnare", "lab", "link", "fires2017", "fires2018" , "mnew", "nare","npar", "ocnew", "parks","pras","pred", "prednew.auc", "prednew.rocr", "rast", "nare", "t", "train", "g", "f1", "pHan", "pHan.rocr", "pHan.auc", "ploss17", "dloss17", "cover", "dcover", "test")))
save.image("parks.Rdata")