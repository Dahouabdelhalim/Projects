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

#requires the file "parks.Rdata" in the same folder in additon to the shapefiles below
#take the previous modeling data
load("parks.Rdata")

#import the raster with forest cover data for 2018
def18 <- raster("defdeg20181.tif")

#change extent
extent(def18) <- extent(parks)

#plot the Dolors data set
tiff(file = "Dolors2018.tiff", width = 2400, height = 1600, units = "px", res = 300)
plot(def18, col=rev(viridis(256)))
dev.off()

#replace values
values(def18) <- ifelse(values(def18) ==3, 1, ifelse(values(def18) ==1,NA, ifelse(values(def18) ==2, 1,0))) 

#resample raster to the resolution of other data
ploss18 <- resample(def18, d2f2017, method='bilinear')

#evaluate only the parks
ploss18<-mask(ploss18, parks)

#convert to data frame
dloss18<-as.data.frame(ploss18, xy=T)

#make a data frame of observation/prediction data for Hansen axis data
#note this is fraction as a result of resampling to lower resolution
pD18<-nare[1:dim(npar)[1],c(1:2, 5:5, 7:7)]
dloss18<-na.omit(dloss18)
pD18<-merge(pD18, dloss18, by.x = c("x", "y"), by.y = c("x","y"))
colnames(pD18)[1:5]<-c("x", "y","difi", "pnew","defo")

#use arbitrary cutpoint to cutoff the axis Hansen observations (not a probabilistic model)
#using the Guaviare cutpoint leads to overestimate of size of deforestation but same AUC
pD18$def<-as.ordered(ifelse(pD18$defo >.5,1,0))

#get AUC for test data and plot
pD18.rocr<- prediction(pD18$pnew,pD18$def)
pD18.auc<-performance(pD18.rocr, measure = "auc")
perf18 <- performance(pD18.rocr, measure = "tpr", x.measure = "fpr") 
pdf("AUC2018.pdf")
plot(perf18)
dev.off()

#print out the summary of Google Earth
sink("def2018lowres.txt")
print(summary(pD18))
sink()

#print out the AUC of the test data
sink("AUC_Meta_data2018.txt")
print(pD18.auc)
sink()

rm(list=setdiff(ls(), c("cfires2017","cfires2018", "clip", "clipp","clipr","d2017", "d2018", "d2f2017", "d2f2018", "defo", "delta1718", "dnare", "lab", "link", "fires2017", "fires2018" , "mnew", "ocnew", "parks","pras","pred", "prednew.auc", "prednew.rocr", "rast", "nare", "npar","t", "train", "g", "f1", "pHan", "pHan.rocr", "pHan.auc", "pD18", "pD18.rocr", "pD18.auc", "ploss17", "dloss17", "cover", "dcover", "def18")))
save.image("parks.Rdata")