library(PresenceAbsence)
library(virtualspecies)
library(SDMTools)
library(raster)
library(downloader)
library(maptools)
library(rgdal)
library(sp)
library(dismo)
library(XML)
library(foreign)
library(jsonlite)
library(GISTools)
library(rJava)
library(ape)
library(randomForest)
library(gam)
library(mcgv)
library(ade4)
library(rworldmap)

jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

setwd("Your_directory/here")

## Study area in Scandinavia
Extent=shapefile("StudyArea.shp")
e = extent(Extent)

Extent=shapefile("StudyArea.shp")
e = extent(Extent)

## Load in elevation DEM
Elev <- raster("Elevation.tif")
Elev <- crop(Elev, e)

## Load in aspect
Aspect <- raster("Aspect.tif")
Aspect <- crop(Aspect, e)
## Load in slope

## Load in forest
Forest <- raster("Percent_Forest.tif")
Forest=crop(Forest, e)

## Line up predictors
Forest <- resample(Forest, Elev)

# Stack predictors
PredictorsVS <- stack(Forest, Elev, Aspect)

plot(PredictorsVS)
legend(lty = c(2,-1,1))
HetHomLS <- randomPoints(Forest, 1000000)
ForVals <- extract(Forest, HetHomLS)
AspVals <- extract(Aspect, HetHomLS)
EleVals <- extract(Elev, HetHomLS)

par(mfrow=c(1,3))
hist(ForVals)
hist(AspVals)
hist(EleVals)

## Stdvs
sd(ForVals)
sd(AspVals)
sd(EleVals)

## Check Finland landscape
sd(ForValsF)
sd(AspValsF)
sd(EleValsF)

## Parameters
my.parameters <- formatFunctions(Percent_Forest = c(fun = 'dnorm', mean = 75, sd = 45),
                                 Elevation = c(fun = 'dnorm', mean = 750, sd = 350),
                                 Aspect = c(fun='dnorm', mean=180, sd=100))


### Make a species suitability out of these rasters
DireWolf <- generateSpFromFun(PredictorsVS, my.parameters, plot = TRUE)
#DireWolf2 <- generateSpFromFun(PredictorsVS, formula = "2*forestfilled2 + ElevationS + AspectS", my.parameters, plot = TRUE)
#DireWolf2 
## Save stuff
#save(DireWolf2, file = "DireWolf")
#load(file="MyVirtualSpecies")
#load(file="DireWolf")
#DireWolf2 <- generateSpFromFun(Predictors200m, my.parameters, plot = TRUE,rescale=FALSE)
## Convert to Presence/Absence
#paDW <- convertToPA(DireWolf2, beta=0.8, alpha = -0.05, plot = TRUE)
paDW <- convertToPA(DireWolf, beta=0.5, alpha = -0.05, plot = TRUE)
#plot(paDW, breaks=seq(min(0),max(1),length.out=200),legend=F,col=rev(terrain.colors(255)))
## See response plots
#plotResponse(DireWolf2)
plotResponse(DireWolf, ylab = " ")

plotResponse(DirewolfF)
str(DireWolf)
## Sample presence Points
#presence.points <- sampleOccurrences(paDW,   n = 200, # The number of points to sample
# type = "presence only")
#presence.points2 <- sampleOccurrences(paDW2, n=200, type = "presence only")
#DWP2 <- presence.points
#DWP3 <- presence.points2
#title(main="Simulated Presence Points")

#DWP2 <- DWP2$sample.points[,1:2]
#DWP3 <- DWP3$sample.points[,1:2]
PA.points <- sampleOccurrences(paDW, n=2000, type="presence-absence")
write.table(data.frame(c(1:length(PA.points$sample.points)),coordinates(PA.points$sample.points)),
            row.names=F,col.names=F,quote=F,sep="\\t",file="DWPresent_Absent")
#DWPresent <- PA.points$sample.points   #NOTE: ALL THIS POUNDED OUT BECAUSE NOW I JUST NEED TO READ TABLE IN EACH TIME
#DWPresent<-subset(DWPresent, Real==1)
#DWPresent <- subset(DWPresent, select=c("x","y"))
#DWPVals <- extract(PredictorsVS, DWPresent)
DWPresent <- read.table("DWPresent_Absent")
DWPresent <- subset(DWPresent, V4==1)
DWPresent <- subset(DWPresent, select=c("V2","V3"))
DWPVals <- extract(PredictorsVS, DWPresent)
### Same for absences
DWAbsent <- read.table("DWPresent_Absent")
DWAbsent <- subset(DWAbsent, V4==0)
DWAbsent <- subset(DWAbsent, select=c("V2","V3"))
DWAVals <- extract(PredictorsVS, DWAbsent)


DWP2 <- DWPresent


### MaxEnt Models at different Grain Sizes
#### aggregate to lower resolutions
# Predictors @ 50m
Predictors50m <- aggregate(PredictorsVS, fact=2, fun=mean)

# Predictors @ 100m
Predictors100m <- aggregate(Predictors50m, fact=2, fun=mean)

# Predictors @ 200m
Predictors200m <- aggregate(Predictors100m, fact=2, fun=mean)

# Predictors @ 400m
Predictors400m <- aggregate(Predictors200m, fact=2, fun=mean)

# Predictors @ 800m
Predictors800m <- aggregate(Predictors400m, fact=2, fun=mean)

## PRedictors @ 1600m
Predictors1600m <- aggregate(Predictors800m, fact=2, fun=mean)

## fix envs error
options( java.parameters = "-Xmx4g" )


bg.DW <- randomPoints(PredictorsVS, 10000, p=DWP2)

### MaxEnt Models of species

### 25m model
eDW <- list()
r.DW <- list()
resultsDW=list()
bg.DW <- randomPoints(PredictorsVS, 10000, p=DWP2)
fold.DW <- kfold(DWP2, k=5)
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DW <- DWP2[fold.DW == i, ]
  occtrain.DW <- DWP2[fold.DW != i, ]
  me.DW <- maxent(PredictorsVS, occtrain.DW,  args=('jackknife=true'))
  plot(me.DW)
  r.DW[[i]] <- predict(me.DW, PredictorsVS)
  eDW[[i]]<- evaluate(me.DW, p=occtest.DW, a=bg.DW, x=PredictorsVS)
  resultsDW[[i]] <- me.DW@results
}
aucDW1 <- sapply( eDW, function(x){slot(x, 'auc')} )
mean(aucDW1)
pDWME <- stack(r.DW)
pDWME <- calc(pDWME, fun=mean)
kappa1 <- sapply( eDW, function(x){slot(x, 'kappa')} )
par(mfrow=c(1,1)) 
maxTPRTNR1 <-  eDW[[1@t(which.max(eDW[[1@TPR]] + eDW[[1@TNR]])]]# threshold at maximum of the sum of the sensitivity (true positive rate) and specificity (true negative rate)

resultsDWdf <- sapply(resultsDW,"[",1:56) 
resultsDWdf <- as.data.frame(resultsDWdf)
names4 <- dimnames(resultsDW[[1]])
names4 <- unlist(names4)
resultsDWdf$variables <- names4
meanresultsDW <- data.frame(variables=resultsDWdf[,6], Means=rowMeans(resultsDWdf[,-6]))
write.csv(meanresultsDW, file="meanValuesDW.csv")

DWTSS <- sapply(eDW, threshold)
DWTSS <- unlist(DWTSS)
DWTSS <- DWTSS[seq(2, length(DWTSS), 6)]
DWTSS <- mean(DWTSS)

par(mar=c(8,4,4,2))
plot(pDWME,main="MaxEnt 25m Resolution Model of Finer Species Distribution",breaks=seq(min(0),max(1),length.out=200),legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)


eDW2 <- list()
r.DW2 <- list()
resultsDW2 <- list()

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DW2 <- DWP2[fold.DW == i, ]
  occtrain.DW2 <- DWP2[fold.DW != i, ]
  me.DW2 <- maxent(Predictors50m, occtrain.DW2, args=('jackknife=true'))
  #plot(me.DW2)
  r.DW2[[i]] <- dismo::predict(me.DW2, Predictors50m)
  eDW2[[i]] <- evaluate(me.DW2, p=occtest.DW2, a=bg.DW, x=Predictors50m)
  resultsDW2[[i]] <- me.DW2@results
}
auc2 <- sapply( eDW2, function(x){slot(x, 'auc')} )
mean(auc2)
pDWME2 <- stack(r.DW2)
pDWME2 <- calc(pDWME2, fun=mean)

## Get TSS out of evaluations
DW2TSS <- sapply(eDW2, threshold)
DW2TSS <- unlist(DW2TSS)
DW2TSS <- DW2TSS[seq(2, length(DW2TSS), 6)]
DW2TSS <- mean(DW2TSS)

resultsDW2df <- sapply(resultsDW2,"[",1:56) 
resultsDW2df <- as.data.frame(resultsDW2df)
names4 <- dimnames(resultsDW2[[1]])
names4 <- unlist(names4)
resultsDW2df$variables <- names4
meanresultsDW2 <- data.frame(variables=resultsDW2df[,6], Means=rowMeans(resultsDW2df[,-6]))
write.csv(meanresultsDW2, file="meanValuesDW2.csv")
par(mfrow=c(1,1)) 
plot(pDWME2,main="MaxEnt 50m Resolution Model of Finer Species Distribution",breaks=seq(min(0),max(1),length.out=200),legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

## MaxEnt @ 100m
eDW3 <- list()
fold.DW3 <- kfold(DWP2, k=5)
r.DW3 <- list()
resultsDW3 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DW3 <- DWP2[fold.DW == i, ]
  occtrain.DW3 <- DWP2[fold.DW != i, ]
  me.DW3 <- maxent(Predictors100m, occtrain.DW3, args=('jackknife=true'))
  plot(me.DW3)
  r.DW3[[i]] <- dismo::predict(me.DW3, Predictors100m)
  eDW3[[i]] <- evaluate(me.DW3, p=occtest.DW3, a=bg.DW, x=Predictors100m)
  resultsDW3[[i]] <- me.DW3@results
}
auc3 <- sapply( eDW3, function(x){slot(x, 'auc')} )
pDWME3 <- stack(r.DW3)
pDWME3 <- calc(pDWME3, fun=mean)

DW3TSS <- sapply(eDW3, threshold)
DW3TSS <- unlist(DW3TSS)
DW3TSS <- DW3TSS[seq(2, length(DW3TSS), 6)]
DW3TSS <- mean(DW3TSS)

resultsDW3df <- sapply(resultsDW3,"[",1:56) 
resultsDW3df <- as.data.frame(resultsDW3df)
names4 <- dimnames(resultsDW3[[1]])
names4 <- unlist(names4)
resultsDW3df$variables <- names4
meanresultsDW3 <- data.frame(variables=resultsDW3df[,6], Means=rowMeans(resultsDW3df[,-6]))
write.csv(meanresultsDW3, file="meanValuesDW3.csv")
par(mfrow=c(1,1)) 
plot(pDWME3,main="MaxEnt 100m Resolution Model of Finer Species Distribution",breaks=seq(min(0),max(1),length.out=200),legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)


### MaxEnt @200m
eDW4 <- list()
fold.DW4 <- kfold(DWP2, k=5)
r.DW4 <- list()
resultsDW4 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DW4 <- DWP2[fold.DW == i, ]
  occtrain.DW4 <- DWP2[fold.DW != i, ]
  me.DW4 <- maxent(Predictors200m, occtrain.DW4, args=('jackknife=true'))
  plot(me.DW4)
  r.DW4[[i]] <- dismo::predict(me.DW4, Predictors200m)
  eDW4[[i]] <- evaluate(me.DW4, p=occtest.DW4, a=bg.DW, x=Predictors200m)
  resultsDW4[[i]] <- me.DW4@results
}
auc4 <- sapply( eDW4, function(x){slot(x, 'auc')} )
pDWME4 <- stack(r.DW4)
pDWME4 <- calc(pDWME4, fun=mean)
resultsDW4df <- sapply(resultsDW4,"[",1:56) 
resultsDW4df <- as.data.frame(resultsDW4df)

DW4TSS <- sapply(eDW4, threshold)
DW4TSS <- unlist(DW4TSS)
DW4TSS <- DW4TSS[seq(2, length(DW4TSS), 6)]
DW4TSS <- mean(DW4TSS)

names4 <- dimnames(resultsDW4[[1]])
names4 <- unlist(names4)
resultsDW4df$variables <- names4
meanresultsDW4 <- data.frame(variables=resultsDW4df[,6], Means=rowMeans(resultsDW4df[,-6]))
write.csv(meanresultsDW4, file="meanValuesDW4.csv")
par(mfrow=c(1,1)) 
plot(pDWME4,main="MaxEnt 200m Resolution Model of Finer Species Distribution",breaks=seq(min(0),max(1),length.out=200),legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

## MaxEnt 400m Model
eDW400 <- list()
fold.DW400 <- kfold(DWP2, k=5)
r.DW400 <- list()
resultsDW400 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DW400 <- DWP2[fold.DW == i, ]
  occtrain.DW400 <- DWP2[fold.DW != i, ]
  me.DW400 <- maxent(Predictors400m, occtrain.DW400, args=('jackknife=true'))
  plot(me.DW400)
  r.DW400[[i]] <- dismo::predict(me.DW400, Predictors400m)
  eDW400[[i]] <- evaluate(me.DW400, p=occtest.DW400, a=bg.DW, x=Predictors400m)
  resultsDW400[[i]] <- me.DW400@results
}
auc400 <- sapply( eDW400, function(x){slot(x, 'auc')} )
mean(auc400)
pDWME400 <- stack(r.DW400)
pDWME400 <- calc(pDWME400, fun=mean)

DW400TSS <- sapply(eDW400, threshold)
DW400TSS <- unlist(DW400TSS)
DW400TSS <- DW400TSS[seq(2, length(DW400TSS), 6)]
DW400TSS <- mean(DW400TSS)

resultsDW400df <- sapply(resultsDW400,"[",1:56) 
resultsDW400df <- as.data.frame(resultsDW400df)
names4 <- dimnames(resultsDW400[[1]])
names4 <- unlist(names4)
resultsDW400df$variables <- names4
meanresultsDW400 <- data.frame(variables=resultsDW400df[,6], Means=rowMeans(resultsDW400df[,-6]))
write.csv(meanresultsDW400, file="meanValuesDW400.csv")
par(mfrow=c(1,1)) 
plot(pDWME400,main="MaxEnt 400m Resolution Model of Finer Species Distribution",breaks=seq(min(0),max(1),length.out=200),legend=F,col=rev(terrain.colors(255)))

y<- c(0:100)*0.005
colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)


## MaxEnt 800m Model
eDW5 <- list()
fold.DW5 <- kfold(DWP2, k=5)
r.DW5 <- list()
resultsDW5 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DW5 <- DWP2[fold.DW == i, ]
  occtrain.DW5 <- DWP2[fold.DW != i, ]
  me.DW5 <- maxent(Predictors800m, occtrain.DW5, args=('jackknife=true'))
  plot(me.DW5)
  r.DW5[[i]] <- dismo::predict(me.DW5, Predictors800m)
  eDW5[[i]] <- evaluate(me.DW5, p=occtest.DW5, a=bg.DW, x=Predictors800m)
  resultsDW5[[i]] <- me.DW5@results
}
auc5 <- sapply( eDW5, function(x){slot(x, 'auc')} )
pDWME5 <- stack(r.DW5)
pDWME5 <- calc(pDWME5, fun=mean)

DW5TSS <- sapply(eDW5, threshold)
DW5TSS <- unlist(DW5TSS)
DW5TSS <- DW5TSS[seq(2, length(DW5TSS), 6)]
DW5TSS <- mean(DW5TSS)

resultsDW5df <- sapply(resultsDW5,"[",1:56) 
resultsDW5df <- as.data.frame(resultsDW5df)
resultsDW5df$variables <- names4
meanresultsDW5 <- data.frame(variables=resultsDW5df[,6], Means=rowMeans(resultsDW5df[,-6]))
write.csv(meanresultsDW5, file="meanValuesDW5.csv")
par(mfrow=c(1,1)) 
plot(pDWME5,main="MaxEnt 800m Resolution Model of Finer Species Distribution",breaks=seq(min(0),max(1),length.out=200),legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

### MaxEnt @1600m
eDW6 = list()
fold.DW6 <- kfold(DWP2, k=5)
r.DW6 <- list()
resultsDW6 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DW6 <- DWP2[fold.DW == i, ]
  occtrain.DW6 <- DWP2[fold.DW != i, ]
  me.DW6 <- maxent(Predictors1600m, occtrain.DW6, args=('jackknife=true'))
  me.DW6
  #plot(me.DW6)
  r.DW6[[i]] <- dismo::predict(me.DW6, Predictors1600m)
  eDW6[[i]] <- evaluate(me.DW6, p=occtest.DW6, a=bg.DW, x=Predictors1600m)
  resultsDW6[[i]] <- me.DW6@results
}
auc6 <- sapply( eDW6, function(x){slot(x, 'auc')} )

resultsDW6df <- sapply(resultsDW6,"[",1:56) 
resultsDW6df <- as.data.frame(resultsDW6df)

DW6TSS <- sapply(eDW6, threshold)
DW6TSS <- unlist(DW6TSS)
DW6TSS <- DW6TSS[seq(2, length(DW6TSS), 6)]
DW6TSS <- mean(DW6TSS)

names4 <- dimnames(resultsDW6[[1]])
names4 <- unlist(names4)
resultsDW6df$variables <- names4
meanresultsDW6 <- data.frame(variables=resultsDW6df[,6], Means=rowMeans(resultsDW6df[,-6]))
write.csv(meanresultsDW6, file="meanValuesDW6.csv")
pDWME6 <- stack(r.DW6)
pDWME6 <- calc(pDWME6, fun=mean)
confm6 <- sapply( eDW6, function(x){slot(x, 'confusion')} )
# threshold at maximum of the sum of the sensitivity (true positive rate) and specificity (true negative rate)
eDW6@t[which.max(eDW4E@TPR + eDW4E@TNR)]
plot(eDW4E@t, eDW4E@TPR+eDW4E@TNR, type="l", lwd=2, cex.lab=1.5)

EDW400t <- sapply( eDW400, function(x){slot(x, 't')} )
write.table(EDW400t[[3]], "threshold_vector.txt")
par(mfrow=c(1,1)) 
plot(pDWME6,main="MaxEnt 1600m Resolution Model of Finer Species Distribution", breaks=seq(min(0),max(1),length.out=200),
     legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)


## Probably run models on all points as well, see what happens. 
#Compute TSS/kappa with confusion matrix from evaluation
mean(auc1)
mean(auc2)
mean(auc3)
mean(auc4)
mean(auc400)
mean(auc5)
mean(auc6)

## Response plots
response(me.DW, main="25 m")
response(me.DW2,main="50 m")
response(me.DW3,main="100 m")
response(me.DW4,main="200 m")
response(me.DW400,main="400 m")
response(me.DW5,main="800 m")
response(me.DW6,main="1600 m")

### Extract points across all predictions (Maybe average prediction rasters somehow in loop?)
## Or perhaps I can just run model on all points, since want to see how well model can do?)
#SuitPoints <- randomPoints(pDWME, 10000)
#write.table(SuitPoints, "SuitPoints.txt")
SuitPoints <- read.table("SuitPoints.txt")
PredictionsDWME25m <- extract(pDWME, SuitPoints)
PredictionsDWME50m <- extract(pDWME2, SuitPoints)
PredictionsDWME100m <- extract(pDWME3, SuitPoints)
PredictionsDWME200m <- extract(pDWME4, SuitPoints)
PredictionsDWME400m <- extract(pDWME400, SuitPoints)
PredictionsDWME800m <- extract(pDWME5, SuitPoints)
PredictionsDWME1600m <- extract(pDWME6, SuitPoints)
# TRUTH!
PredictionsDWtrue <- extract(DireWolf$suitab.raster, SuitPoints)

## Correlation between predictions and truth
par(mfrow=c(4,2), pty="s", mai=c(.5,.2,.3,.3))

cor(PredictionsDWME25m, PredictionsDWtrue)
plot(PredictionsDWME25m,PredictionsDWtrue, xlim=c(0,0.8), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("25 m (Pearson's r= ", formatC(cor(PredictionsDWME25m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2),cex.lab=1.5,cex.main=2, cex.axis=1.5)
abline(a=0,b=1, col="green3")
cor(PredictionsDWME50m, PredictionsDWtrue)
plot(PredictionsDWME50m, PredictionsDWtrue, xlim=c(0,0.8), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("50 m (Pearson's r= ", formatC(cor(PredictionsDWME50m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2),cex.lab=1.5,cex.main=2, cex.axis=1.5)
abline(a=0,b=1, col="green3")
cor(PredictionsDWME100m, PredictionsDWtrue)
plot(PredictionsDWME100m, PredictionsDWtrue, xlim=c(0,0.8), ylim=c(0,1),xlab="Prediction", ylab="Truth",main=paste0("100 m (Pearson's r= ", formatC(cor(PredictionsDWME100m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2),cex.lab=1.5,cex.main=2, cex.axis=1.5)
abline(a=0,b=1, col="green3")
(PredictionsDWME200m, PredictionsDWtrue)
plot(PredictionsDWME200m, PredictionsDWtrue, xlim=c(0,0.8), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("200 m (Pearson's r= ", formatC(cor(PredictionsDWME200m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2),cex.lab=1.5,cex.main=2, cex.axis=1.5)
abline(a=0,b=1, col="green3")
cor(PredictionsDWME400m, PredictionsDWtrue)
plot(PredictionsDWME400m, PredictionsDWtrue,xlim=c(0,0.8), ylim=c(0,1), xlab="Prediction", ylab="Truth", main=paste0("400 m (Pearson's r= ", formatC(cor(PredictionsDWME400m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2),cex.lab=1.5,cex.main=2, cex.axis=1.5)
abline(a=0,b=1, col="green3")
cor(PredictionsDWME800m, PredictionsDWtrue)
plot(PredictionsDWME800m, PredictionsDWtrue, xlim=c(0,0.8), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("800 m (Pearson's r= ", formatC(cor(PredictionsDWME800m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2),cex.lab=1.5,cex.main=2, cex.axis=1.5)
abline(a=0,b=1, col="green3")
cor(PredictionsDWME1600m, PredictionsDWtrue)
plot(PredictionsDWME1600m,PredictionsDWtrue, xlim=c(0,0.8), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("1600 m (Pearson's r= ", formatC(cor(PredictionsDWME1600m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2),cex.lab=1.5,cex.main=2, cex.axis=1.5)
abline(a=0,b=1, col="green3")

# Mean Predictions
mean(PredictionsDWME25m)
mean(PredictionsDWME50m)
mean(PredictionsDWME100m)
mean(PredictionsDWME200m)
mean(PredictionsDWME400m)
mean(PredictionsDWME800m)
mean(PredictionsDWME1600m)
# TRUTH!
mean(PredictionsDWtrue)

sd(PredictionsDWME25m)
sd(PredictionsDWME50m)
sd(PredictionsDWME100m)
sd(PredictionsDWME200m)
sd(PredictionsDWME400m)
sd(PredictionsDWME800m)
sd(PredictionsDWME1600m)
# TRUTH!
sd(PredictionsDWtrue)
DWdev <- calc.deviance(PredictionsDWtrue, PredictionsDWME25m,weights = rep(1,length(PredictionsDWtrue)), 
                       family="poisson", calc.mean = TRUE) 
DWdev2 <- calc.deviance(PredictionsDWtrue, PredictionsDWME50m,weights = rep(1,length(PredictionsDWtrue)), 
                        family="poisson", calc.mean = TRUE)
DWdev3 <- calc.deviance(PredictionsDWtrue, PredictionsDWME100m,weights = rep(1,length(PredictionsDWtrue)), 
                        family="poisson", calc.mean = TRUE)
DWdev4 <- calc.deviance(PredictionsDWtrue, PredictionsDWME200m,weights = rep(1,length(PredictionsDWtrue)), 
                        family="poisson", calc.mean = TRUE)
DWdev5 <- calc.deviance(PredictionsDWtrue, PredictionsDWME400m,weights = rep(1,length(PredictionsDWtrue)), 
                        family="poisson", calc.mean = TRUE)
DWdev6 <- calc.deviance(PredictionsDWtrue, PredictionsDWME800m,weights = rep(1,length(PredictionsDWtrue)), 
                        family="poisson", calc.mean = TRUE)
DWdev7 <- calc.deviance(PredictionsDWtrue, PredictionsDWME1600m,weights = rep(1,length(PredictionsDWtrue)), 
                        family="poisson", calc.mean = TRUE)



## Build the same species, but at a coarser resolution
#DirewolfF3 <-  generateSpFromFun(Predictors250m, formula = "2*forestfilled2 + ElevationS + AspectS", my.parameters, plot = TRUE)
DirewolfF <- generateSpFromFun(Predictors200m, my.parameters, plot = TRUE)
## Convert to Presence/Absence
paDWF <- convertToPA(DirewolfF, beta=0.5, alpha = -0.05, plot = TRUE)
#plot(paDWF, breaks=seq(min(0),max(1),length.out=200),
# legend=F,col=rev(terrain.colors(255)))

## Make Points, presences
PA.pointsF <- sampleOccurrences(paDWF, n=2000, type="presence-absence")


write.table(data.frame(c(1:length(PA.pointsF$sample.points)),coordinates(PA.pointsF$sample.points)),
            row.names=F,col.names=F,quote=F,sep="\\t",file="DWFPresent_Absent")

DWFPresent <- read.table("DWFPresent_Absent")
DWFPresent <- subset(DWFPresent, V4==1)
DWFPresent <- subset(DWFPresent, select=c("V2","V3"))
DWFPVals <- extract(PredictorsVS, DWFPresent)
### Same for absences
DWFAbsent <- read.table("DWFPresent_Absent")
DWFAbsent <- subset(DWFAbsent, V4==0)
DWFAbsent <- subset(DWFAbsent, select=c("V2","V3"))
DWFAVals <- extract(PredictorsVS, DWFAbsent)

DWPF <- DWFPresent

## See response plots
plotResponse(DirewolfF)

## Save stuff
save(DirewolfF, file = "FatDirewolf")
load(file ="FatDirewolf")

#presence.pointsF <- sampleOccurrences(paDWF,   n = 200, # The number of points to sample
# type = "presence only")
#DWPF <- presence.pointsF
#title(main="Simulated Presence Points")

#DWPF <- DWPF$sample.points[,1:2]

## MaxEnt 25m for Direwolf 2
eDWF <- list()
r.DWF <- list()
resultsDWF <- list()
#bg.DW <- randomPoints(PredictorsVS, 10000, p=DWPF)
fold.DWF <- kfold(DWPF, k=5)

## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DWF <- DWPF[fold.DWF == i, ]
  occtrain.DWF <- DWPF[fold.DWF != i, ]
  me.DWF <- maxent(PredictorsVS, occtrain.DWF, args=('jackknife=true'))
  #plot(me.DW)
  r.DWF[[i]] <- predict(me.DWF, PredictorsVS)
  eDWF[[i]]<- evaluate(me.DWF, p=occtest.DWF, a=bg.DW, x=PredictorsVS)
  resultsDWF[[i]] <- me.DWF@results
}
aucF <- sapply( eDWF, function(x){slot(x, 'auc')} )
mean(aucF)
kappaF <- sapply( eDWF, function(x){slot(x, 'kappa')} )
pDWMEF <- stack(r.DWF)
pDWMEF <- calc(pDWMEF, fun=mean)

DWFTSS <- sapply(eDWF, threshold)
DWFTSS <- unlist(DWFTSS)
DWFTSS <- DWFTSS[seq(2, length(DWFTSS), 6)]
DWFTSS <- mean(DWFTSS)

resultsDWFdf <- sapply(resultsDWF,"[",1:56) 
resultsDWFdf <- as.data.frame(resultsDWFdf)
names4 <- dimnames(resultsDWF[[1]])
names4 <- unlist(names4)
resultsDWFdf$variables <- names4
meanresultsDWF <- data.frame(variables=resultsDWFdf[,6], Means=rowMeans(resultsDWFdf[,-6]))
write.csv(meanresultsDWF, file="meanValuesDWF.csv")
par(mfrow=c(1,1)) 
plot(pDWMEF,main="MaxEnt 25m Resolution Model of Coarser Species Distribution" ,breaks=seq(min(0),max(1),length.out=200),
     legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)
### MaxEnt @ 50m
eDWF2 <- list()
fold.DWF2 <- kfold(DWPF, k=5)
r.DWF2 <- list()
resultsDWF2 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DWF2 <- DWPF[fold.DWF == i, ]
  occtrain.DWF2 <- DWPF[fold.DWF != i, ]
  me.DWF2 <- maxent(Predictors50m, occtrain.DWF2, args=('jackknife=true'))
  #plot(me.DW2)
  r.DWF2[[i]] <- dismo::predict(me.DWF2, Predictors50m)
  eDWF2[[i]] <- evaluate(me.DWF2, p=occtest.DWF2, a=bg.DW, x=Predictors50m)
  resultsDWF2[[i]] <- me.DWF2@results
}
aucF2 <- sapply( eDWF2, function(x){slot(x, 'auc')} )
mean(aucF2)
pDWMEF2 <- stack(r.DWF2)
pDWMEF2 <- calc(pDWMEF2, fun=mean)
resultsDWF2df <- sapply(resultsDWF2,"[",1:56) 
resultsDWF2df <- as.data.frame(resultsDWF2df)

DWF2TSS <- sapply(eDWF2, threshold)
DWF2TSS <- unlist(DWF2TSS)
DWF2TSS <- DWF2TSS[seq(2, length(DWF2TSS), 6)]
DWF2TSS <- mean(DWF2TSS)

names4 <- dimnames(resultsDWF2[[1]])
names4 <- unlist(names4)
resultsDWF2df$variables <- names4
meanresultsDWF2 <- data.frame(variables=resultsDWF2df[,6], Means=rowMeans(resultsDWF2df[,-6]))
write.csv(meanresultsDWF2, file="meanValuesDWF2.csv")
par(mfrow=c(1,1)) 
plot(pDWMEF2,main="MaxEnt 50m Resolution Model of Coarser Species Distribution",  breaks=seq(min(0),max(1),length.out=200),
     legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)
## MaxEnt @ 100m
eDWF3 <- list()
fold.DWF3 <- kfold(DWPF, k=5)
r.DWF3 <- list()
resultsDWF3 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DWF3 <- DWPF[fold.DWF == i, ]
  occtrain.DWF3 <- DWPF[fold.DWF != i, ]
  me.DWF3 <- maxent(Predictors100m, occtrain.DWF3, args=('jackknife=true'))
  #plot(me.DW2)
  r.DWF3[[i]] <- dismo::predict(me.DWF3, Predictors100m)
  eDWF3[[i]] <- evaluate(me.DWF3, p=occtest.DWF3, a=bg.DW, x=Predictors100m)
  resultsDWF3[[i]] <- me.DWF3@results
}
aucF3 <- sapply( eDWF3, function(x){slot(x, 'auc')} )
mean(aucF3)
pDWMEF3 <- stack(r.DWF3)
pDWMEF3 <- calc(pDWMEF3, fun=mean)
resultsDWF3df <- sapply(resultsDWF3,"[",1:56) 
resultsDWF3df <- as.data.frame(resultsDWF3df)

DWF3TSS <- sapply(eDWF3, threshold)
DWF3TSS <- unlist(DWF3TSS)
DWF3TSS <- DWF3TSS[seq(2, length(DWF3TSS), 6)]
DWF3TSS <- mean(DWF3TSS)

names4 <- dimnames(resultsDWF3[[1]])
names4 <- unlist(names4)
resultsDWF3df$variables <- names4
meanresultsDWF3 <- data.frame(variables=resultsDWF3df[,6], Means=rowMeans(resultsDWF3df[,-6]))
write.csv(meanresultsDWF3, file="meanValuesDWF3.csv")
par(mfrow=c(1,1)) 
plot(pDWMEF3,main="MaxEnt 100m Resolution Model of Coarser Species Distribution", breaks=seq(min(0),max(1),length.out=200),
     legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

## MaxEnt @ 200m
eDWF4 <- list()
fold.DWF4 <- kfold(DWPF, k=5)
r.DWF4 <- list()
resultsDWF4 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DWF4 <- DWPF[fold.DWF == i, ]
  occtrain.DWF4 <- DWPF[fold.DWF != i, ]
  me.DWF4 <- maxent(Predictors200m, occtrain.DWF4, args=('jackknife=true'))
  #plot(me.DW2)
  r.DWF4[[i]] <- dismo::predict(me.DWF4, Predictors200m)
  eDWF4[[i]] <- evaluate(me.DWF4, p=occtest.DWF4, a=bg.DW, x=Predictors200m)
  resultsDWF4[[i]] <- me.DWF4@results
}
aucF4 <- sapply( eDWF4, function(x){slot(x, 'auc')} )
mean(aucF4)
pDWMEF4 <- stack(r.DWF4)
pDWMEF4 <- calc(pDWMEF4, fun=mean)
resultsDWF4df <- sapply(resultsDWF4,"[",1:56) 
resultsDWF4df <- as.data.frame(resultsDWF4df)

DWF4TSS <- sapply(eDWF4, threshold)
DWF4TSS <- unlist(DWF4TSS)
DWF4TSS <- DWF4TSS[seq(2, length(DWF4TSS), 6)]
DWF4TSS <- mean(DWF4TSS)

names4 <- dimnames(resultsDWF4[[1]])
names4 <- unlist(names4)
resultsDWF4df$variables <- names4
meanresultsDWF4 <- data.frame(variables=resultsDWF4df[,6], Means=rowMeans(resultsDWF4df[,-6]))
write.csv(meanresultsDWF4, file="meanValuesDWF4.csv")
par(mfrow=c(1,1)) 
plot(pDWMEF4,main="MaxEnt 200m Resolution Model of Species Distribution", breaks=seq(min(0),max(1),length.out=200),
     legend=F,col=rev(terrain.colors(255)))

y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

## MaxEnt @ 400m
eDWF400 <- list()
fold.DWF400 <- kfold(DWPF, k=5)
r.DWF400 <- list()
resultsDWF400 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DWF400 <- DWPF[fold.DWF == i, ]
  occtrain.DWF400 <- DWPF[fold.DWF != i, ]
  me.DWF400 <- maxent(Predictors400m, occtrain.DWF400, args=('jackknife=true'))
  #plot(me.DW2)
  r.DWF400[[i]] <- dismo::predict(me.DWF400, Predictors400m)
  eDWF400[[i]] <- evaluate(me.DWF400, p=occtest.DWF400, a=bg.DW, x=Predictors400m)
  resultsDWF400[[i]] <- me.DWF400@results
}
aucF400 <- sapply( eDWF400, function(x){slot(x, 'auc')} )
mean(aucF400)
pDWMEF400 <- stack(r.DWF400)
pDWMEF400 <- calc(pDWMEF400, fun=mean)
resultsDWF400df <- sapply(resultsDWF400,"[",1:56) 
resultsDWF400df <- as.data.frame(resultsDWF400df)

DWF400TSS <- sapply(eDWF400, threshold)
DWF400TSS <- unlist(DWF400TSS)
DWF400TSS <- DWF400TSS[seq(2, length(DWF400TSS), 6)]
DWF400TSS <- mean(DWF400TSS)

names4 <- dimnames(resultsDWF400[[1]])
names4 <- unlist(names4)
resultsDWF400df$variables <- names4
meanresultsDWF400 <- data.frame(variables=resultsDWF400df[,6], Means=rowMeans(resultsDWF400df[,-6]))
write.csv(meanresultsDWF400, file="meanValuesDWF400.csv")
par(mfrow=c(1,1)) 
plot(pDWMEF400,main="MaxEnt 400m Resolution Model of Coarser Species Distribution", breaks=seq(min(0),max(1),length.out=200),
     legend=F,col=rev(terrain.colors(255)))

## MaxEnt @ 800m
eDWF5 <- list()
fold.DWF5 <- kfold(DWPF, k=5)
r.DWF5 <- list()
resultsDWF5 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DWF5 <- DWPF[fold.DWF == i, ]
  occtrain.DWF5 <- DWPF[fold.DWF != i, ]
  me.DWF5 <- maxent(Predictors800m, occtrain.DWF5, args=('jackknife=true'))
  #plot(me.DW2)
  r.DWF5[[i]] <- dismo::predict(me.DWF5, Predictors800m)
  eDWF5[[i]] <- evaluate(me.DWF5, p=occtest.DWF5, a=bg.DW, x=Predictors800m)
  resultsDWF5[[i]] <- me.DWF5@results
}
aucF5 <- sapply( eDWF5, function(x){slot(x, 'auc')} )
mean(aucF5)
pDWMEF5 <- stack(r.DWF5)
pDWMEF5 <- calc(pDWMEF5, fun=mean)
resultsDWF5df <- sapply(resultsDWF5,"[",1:56) 
resultsDWF5df <- as.data.frame(resultsDWF5df)

DWF5TSS <- sapply(eDWF5, threshold)
DWF5TSS <- unlist(DWF5TSS)
DWF5TSS <- DWF5TSS[seq(2, length(DWF5TSS), 6)]
DWF5TSS <- mean(DWF5TSS)

names4 <- dimnames(resultsDWF5[[1]])
names4 <- unlist(names4)
resultsDWF5df$variables <- names4
meanresultsDWF5 <- data.frame(variables=resultsDWF5df[,6], Means=rowMeans(resultsDWF5df[,-6]))
write.csv(meanresultsDWF5, file="meanValuesDWF5.csv")
par(mfrow=c(1,1)) 
plot(pDWMEF5,main="MaxEnt 800m Resolution Model of Coarser Species Distribution", breaks=seq(min(0),max(1),length.out=200),
     legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

## MaxEnt @ 1600m
eDWF6 <- list()
fold.DWF6 <- kfold(DWPF, k=5)
r.DWF6 <- list()
resultsDWF6 <- list()
## Create training and testing data (take out 20% for training, test with that 20%)

par(mfrow=c(1,3)) 
for(i in 1:5){
  occtest.DWF6 <- DWPF[fold.DWF == i, ]
  occtrain.DWF6 <- DWPF[fold.DWF != i, ]
  me.DWF6 <- maxent(Predictors1600m, occtrain.DWF6, args=('jackknife=true'))
  #plot(me.DW2)
  r.DWF6[[i]]<- dismo::predict(me.DWF6, Predictors1600m)
  eDWF6[[i]] <- evaluate(me.DWF6, p=occtest.DWF6, a=bg.DW, x=Predictors1600m)
  resultsDWF6[[i]] <- me.DWF6@results
}
aucF6 <- sapply( eDWF6, function(x){slot(x, 'auc')} )
mean(aucF6)
pDWMEF6 <- stack(r.DWF6)
pDWMEF6 <- calc(pDWMEF6, fun=mean)
resultsDWF6df <- sapply(resultsDWF6,"[",1:56) 
resultsDWF6df <- as.data.frame(resultsDWF6df)

DWF6TSS <- sapply(eDWF6, threshold)
DWF6TSS <- unlist(DWF6TSS)
DWF6TSS <- DWF6TSS[seq(2, length(DWF6TSS), 6)]
DWF6TSS <- mean(DWF6TSS)

names4 <- dimnames(resultsDWF6[[1]])
names4 <- unlist(names4)
resultsDWF6df$variables <- names4
meanresultsDWF6 <- data.frame(variables=resultsDWF6df[,6], Means=rowMeans(resultsDWF6df[,-6]))
write.csv(meanresultsDWF6, file="meanValuesDWF6.csv")
evalDWF6df <- sapply()
par(mfrow=c(1,1)) 
plot(pDWMEF6,main="MaxEnt 1600m Resolution Model of Coarser Species Distribution", breaks=seq(min(0),max(1),length.out=200),
     legend=F,col=rev(terrain.colors(255)))
y<- c(0:100)*0.005

colorbar.plot( 4450000, 4120000, y, col=rev(terrain.colors(255)), strip.width=0.02, strip.length=(.5), horizontal=FALSE )
text( 4450000, 4030000, labels = 0, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4110000, labels = 0.5, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

text( 4450000, 4185000, labels = 1, adj = NULL,
      pos = 4, offset = 0.5, vfont = NULL,
      cex = 1, col = NULL, font = NULL)

## Comparisons
mean(aucF)
mean(aucF2)
mean(aucF3)
mean(aucF4)
mean(aucF400)
mean(aucF5)
mean(aucF6)

## Response plots
response(me.DWF, main="25 m")
response(me.DWF2,main="50 m")
response(me.DWF3,main="100 m")
response(me.DWF4,main="200 m")
response(me.DWF400,main="400 m")
response(me.DWF5,main="800 m")
response(me.DWF6,main="1600 m")


PredictionsDWFME25m <- extract(pDWMEF, SuitPoints)
PredictionsDWFME50m <- extract(pDWMEF2, SuitPoints)
PredictionsDWFME100m <- extract(pDWMEF3, SuitPoints)
PredictionsDWFME200m <- extract(pDWMEF4, SuitPoints)
PredictionsDWFME400m <- extract(pDWMEF400, SuitPoints)
PredictionsDWFME800m <- extract(pDWMEF5, SuitPoints)
PredictionsDWFME1600m <- extract(pDWMEF6, SuitPoints)
# TRUTH!
PredictionsDWFtrue <- extract(DirewolfF$suitab.raster, SuitPoints)

## Correlation between predictions and truth
par(mfrow=c(4,2), pty="s", mai=c(.3,.15,.15,.2))

cor(PredictionsDWFME25m, PredictionsDWFtrue)
plot(PredictionsDWFME25m,PredictionsDWFtrue, xlim=c(0,1), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("25 m (Pearson's r= ", formatC(cor(PredictionsDWFME25m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFME50m, PredictionsDWFtrue)
plot(PredictionsDWFME50m, PredictionsDWFtrue, xlim=c(0,1), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("50 m (Pearson's r= ", formatC(cor(PredictionsDWFME50m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFME100m, PredictionsDWFtrue)
plot(PredictionsDWFME100m, PredictionsDWFtrue, xlim=c(0,1), ylim=c(0,1),xlab="Prediction", ylab="Truth",main=paste0("100 m (Pearson's r= ", formatC(cor(PredictionsDWFME100m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFME200m, PredictionsDWFtrue)
plot(PredictionsDWFME200m, PredictionsDWFtrue, xlim=c(0,1), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("200 m (Pearson's r= ", formatC(cor(PredictionsDWFME200m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFME400m, PredictionsDWFtrue)
plot(PredictionsDWFME400m, PredictionsDWFtrue, xlim=c(0,1), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("400 m (Pearson's r= ", formatC(cor(PredictionsDWFME400m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFME800m, PredictionsDWFtrue)
plot(PredictionsDWFME800m, PredictionsDWFtrue, xlim=c(0,1), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("800 m (Pearson's r= ", formatC(cor(PredictionsDWFME800m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFME1600m, PredictionsDWFtrue)
plot(PredictionsDWFME1600m,PredictionsDWFtrue, xlim=c(0,1), ylim=c(0,1),xlab="Prediction", ylab="Truth", main=paste0("1600 m (Pearson's r= ", formatC(cor(PredictionsDWFME1600m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")

# Mean Predictions
mean(PredictionsDWFME25m)
mean(PredictionsDWFME50m)
mean(PredictionsDWFME100m)
mean(PredictionsDWFME200m)
mean(PredictionsDWFME400m)
mean(PredictionsDWFME800m)
mean(PredictionsDWFME1600m)
# TRUTH!
mean(PredictionsDWFtrue)

sd(PredictionsDWFME25m)
sd(PredictionsDWFME50m)
sd(PredictionsDWFME100m)
sd(PredictionsDWFME200m)
sd(PredictionsDWFME400m)
sd(PredictionsDWFME800m)
sd(PredictionsDWFME1600m)
# TRUTH!
sd(PredictionsDWFtrue)

### GLM Models of virtual species
# Fine Grain Species



### Format data
DWPVals <- as.data.frame(DWPVals)
DWAVals <- as.data.frame(DWAVals)

DWp <- c(rep(1, nrow(DWPVals)))
DWp <- data.frame(cbind(DWp, rbind(DWPVals)))
DWa <- c(rep(0, nrow(DWAVals)))
DWa <- data.frame(cbind(DWa, rbind(DWAVals)))

pafoldDWp <- DWPVals
pafoldDWa <- DWAVals
k <- 5
groupp <- kfold(pafoldDWp, k)
groupa <- kfold(pafoldDWa, k)

eDWGLM <- list()
pDWGLM <- list()
for (i in 1:k) {
  ptrainDW <- DWp[groupp !=i,]
  atrainDW <- DWa[groupa !=i,]
  ptestDW <- DWp[groupp == i,]
  atestDW <- DWa[groupa == i,]
  pa <- c(rep(1, nrow(ptrainDW)), rep(0, nrow(atrainDW)))
  sdmdatatrainDW <- data.frame(cbind(pa, rbind(ptrainDW[2:4], atrainDW[2:4])))
  pmDW <- glm(pa ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
              * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDW, 
              family=binomial)
  pmDW <- step(pmDW)
  pDWGLM[[i]] <- predict(PredictorsVS, pmDW, type="response")
  eDWGLM[[i]] <- evaluate(pmDW, p=ptestDW, a=atestDW)
}
eDWGLM
pDWGLMm <- stack(pDWGLM)
PredictionsDWGLM <- calc(pDWGLMm, fun=mean)
aucDWGLM <- sapply( eDWGLM, function(x){slot(x, 'auc')} )
mean(aucDWGLM)
plot(PredictionsDWGLM)

DWGLMTSS <- sapply(eDWGLM, threshold)
DWGLMTSS <- unlist(DWGLMTSS)
DWGLMTSS <- DWGLMTSS[seq(2, length(DWGLMTSS), 6)]
DWGLMTSS <- mean(DWGLMTSS)

### GLM 50m
DWPVals2 <- extract(Predictors50m, DWPresent)
DWAVals2 <- extract(Predictors50m, DWAbsent)

### Format data
DWPVals2 <- as.data.frame(DWPVals2)
DWAVals2 <- as.data.frame(DWAVals2)
DWp2 <- c(rep(1, nrow(DWPVals2)))
DWp2 <- data.frame(cbind(DWp2, rbind(DWPVals2)))
DWa2 <- c(rep(0, nrow(DWAVals2)))
DWa2 <- data.frame(cbind(DWa2, rbind(DWAVals2)))

pafoldDWp2 <- DWPVals2
pafoldDWa2 <- DWAVals2
k <- 5
group2p <- kfold(pafoldDWp2, k)
group2a <- kfold(pafoldDWa2, k)

eDWGLM2 <- list()
pDWGLM2 <- list()
for (i in 1:k) {
  ptrainDW2 <- DWp2[group2p !=i,]
  atrainDW2 <- DWa2[group2a !=i,]
  ptestDW2 <- DWp2[group2p == i,]
  atestDW2 <- DWa2[group2a == i,]
  pa2 <- c(rep(1, nrow(ptrainDW2)), rep(0, nrow(atrainDW2)))
  sdmdatatrainDW2 <- data.frame(cbind(pa2, rbind(ptrainDW2[2:4], atrainDW2[2:4])))
  pmDW2 <- glm(pa2 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
               * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDW2, 
               family=binomial)
  pmDW2 <- step(pmDW2)
  pDWGLM2[[i]] <- predict(Predictors50m, pmDW2, type="response")
  eDWGLM2[[i]] <- evaluate(pmDW2, p=ptestDW2, a=atestDW2)
}
pDWGLMm2 <- stack(pDWGLM2)
PredictionsDWGLM2 <- calc(pDWGLMm2, fun=mean)
aucDWGLM2 <- sapply( eDWGLM2, function(x){slot(x, 'auc')} )
mean(aucDWGLM2)
response(pmDW2)
plot(PredictionsDWGLM2)

DWGLM2TSS <- sapply(eDWGLM2, threshold)
DWGLM2TSS <- unlist(DWGLM2TSS)
DWGLM2TSS <- DWGLM2TSS[seq(2, length(DWGLM2TSS), 6)]
DWGLM2TSS <- mean(DWGLM2TSS)
### GLM 100m
DWPVals3 <- extract(Predictors100m, DWPresent)
DWAVals3 <- extract(Predictors100m, DWAbsent)

### Format data
DWPVals3 <- as.data.frame(DWPVals3)
DWAVals3 <- as.data.frame(DWAVals3)
DWp3 <- c(rep(1, nrow(DWPVals3)))
DWp3 <- data.frame(cbind(DWp3, rbind(DWPVals3)))
DWa3 <- c(rep(0, nrow(DWAVals3)))
DWa3 <- data.frame(cbind(DWa3, rbind(DWAVals3)))

pafoldDWp3 <- DWPVals3
pafoldDWa3 <- DWAVals3
k <- 5
group3p <- kfold(pafoldDWp3, k)
group3a <- kfold(pafoldDWa3, k)

eDWGLM3 <- list()
pDWGLM3 <- list()
for (i in 1:k) {
  ptrainDW3 <- DWp3[group3p !=i,]
  atrainDW3 <- DWa3[group3a !=i,]
  ptestDW3 <- DWp3[group3p == i,]
  atestDW3 <- DWa3[group3a == i,]
  pa3 <- c(rep(1, nrow(ptrainDW3)), rep(0, nrow(atrainDW3)))
  sdmdatatrainDW3 <- data.frame(cbind(pa3, rbind(ptrainDW3[2:4], atrainDW3[2:4])))
  pmDW3 <- glm(pa3 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
               * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDW3, 
               family=binomial)
  pmDW3 <- step(pmDW3)
  pDWGLM3[[i]] <- predict(Predictors100m, pmDW3, type="response")
  eDWGLM3[[i]] <- evaluate(pmDW3, p=ptestDW3, a=atestDW3)
}
pDWGLMm3 <- stack(pDWGLM3)
PredictionsDWGLM3 <- calc(pDWGLMm3, fun=mean)
aucDWGLM3 <- sapply( eDWGLM3, function(x){slot(x, 'auc')} )
mean(aucDWGLM3)
response(pmDW3)
plot(PredictionsDWGLM3)
DWGLM3TSS <- sapply(eDWGLM3, threshold)
DWGLM3TSS <- unlist(DWGLM3TSS)
DWGLM3TSS <- DWGLM3TSS[seq(2, length(DWGLM3TSS), 6)]
DWGLM3TSS <- mean(DWGLM3TSS)

### GLM 200m
DWPVals4 <- extract(Predictors200m, DWPresent)
DWAVals4 <- extract(Predictors200m, DWAbsent)

### Format data
DWPVals4 <- as.data.frame(DWPVals4)
DWAVals4 <- as.data.frame(DWAVals4)
DWp4 <- c(rep(1, nrow(DWPVals4)))
DWp4 <- data.frame(cbind(DWp4, rbind(DWPVals4)))
DWa4 <- c(rep(0, nrow(DWAVals4)))
DWa4 <- data.frame(cbind(DWa4, rbind(DWAVals4)))

pafoldDWp4 <- DWPVals4
pafoldDWa4 <- DWAVals4
k <- 5
group4p <- kfold(pafoldDWp4, k)
group4a <- kfold(pafoldDWa4, k)

eDWGLM4 <- list()
pDWGLM4 <- list()
for (i in 1:k) {
  ptrainDW4 <- DWp4[group4p !=i,]
  atrainDW4 <- DWa4[group4a !=i,]
  ptestDW4 <- DWp4[group4p == i,]
  atestDW4 <- DWa4[group4a == i,]
  pa4 <- c(rep(1, nrow(ptrainDW4)), rep(0, nrow(atrainDW4)))
  sdmdatatrainDW4 <- data.frame(cbind(pa4, rbind(ptrainDW4[2:4], atrainDW4[2:4])))
  pmDW4 <- glm(pa4 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
               * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDW4, 
               family=binomial)
  pmDW4 <- step(pmDW4)
  pDWGLM4[[i]] <- predict(Predictors200m, pmDW4, type="response")
  eDWGLM4[[i]] <- evaluate(pmDW4, p=ptestDW4, a=atestDW4)
}
pDWGLMm4 <- stack(pDWGLM4)
PredictionsDWGLM4 <- calc(pDWGLMm4, fun=mean)
aucDWGLM4 <- sapply( eDWGLM4, function(x){slot(x, 'auc')} )
mean(aucDWGLM4)
response(pmDW4)
plot(PredictionsDWGLM4)

DWGLM4TSS <- sapply(eDWGLM4, threshold)
DWGLM4TSS <- unlist(DWGLM4TSS)
DWGLM4TSS <- DWGLM4TSS[seq(2, length(DWGLM4TSS), 6)]
DWGLM4TSS <- mean(DWGLM4TSS)

### GLM 400m
DWPVals5 <- extract(Predictors400m, DWPresent)
DWAVals5 <- extract(Predictors400m, DWAbsent)

### Format data
DWPVals5 <- as.data.frame(DWPVals5)
DWAVals5 <- as.data.frame(DWAVals5)
DWp5 <- c(rep(1, nrow(DWPVals5)))
DWp5 <- data.frame(cbind(DWp5, rbind(DWPVals5)))
DWa5 <- c(rep(0, nrow(DWAVals5)))
DWa5 <- data.frame(cbind(DWa5, rbind(DWAVals5)))

pafoldDWp5 <- DWPVals5
pafoldDWa5 <- DWAVals5
k <- 5
group5p <- kfold(pafoldDWp5, k)
group5a <- kfold(pafoldDWa5, k)

eDWGLM5 <- list()
pDWGLM5 <- list()
for (i in 1:k) {
  ptrainDW5 <- DWp5[group5p !=i,]
  atrainDW5 <- DWa5[group5a !=i,]
  ptestDW5 <- DWp5[group5p == i,]
  atestDW5 <- DWa5[group5a == i,]
  pa5 <- c(rep(1, nrow(ptrainDW5)), rep(0, nrow(atrainDW5)))
  sdmdatatrainDW5 <- data.frame(cbind(pa5, rbind(ptrainDW5[2:4], atrainDW5[2:4])))
  pmDW5 <- glm(pa5 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
               * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDW5, 
               family=binomial)
  pmDW5 <- step(pmDW5)
  pDWGLM5[[i]] <- predict(Predictors400m, pmDW5, type="response")
  eDWGLM5[[i]] <- evaluate(pmDW5, p=ptestDW5, a=atestDW5)
}
pDWGLMm5 <- stack(pDWGLM5)
PredictionsDWGLM5 <- calc(pDWGLMm5, fun=mean)
aucDWGLM5 <- sapply( eDWGLM5, function(x){slot(x, 'auc')} )
mean(aucDWGLM5)
response(pmDW5)
plot(PredictionsDWGLM5)

DWGLM5TSS <- sapply(eDWGLM5, threshold)
DWGLM5TSS <- unlist(DWGLM5TSS)
DWGLM5TSS <- DWGLM5TSS[seq(2, length(DWGLM5TSS), 6)]
DWGLM5TSS <- mean(DWGLM5TSS)

### GLM 800
DWPVals6 <- extract(Predictors800m, DWPresent)
DWAVals6 <- extract(Predictors800m, DWAbsent)

### Format data
DWPVals6 <- as.data.frame(DWPVals6)
DWAVals6 <- as.data.frame(DWAVals6)
DWp6 <- c(rep(1, nrow(DWPVals6)))
DWp6 <- data.frame(cbind(DWp6, rbind(DWPVals6)))
DWa6 <- c(rep(0, nrow(DWAVals6)))
DWa6 <- data.frame(cbind(DWa6, rbind(DWAVals6)))

pafoldDWp6 <- DWPVals6
pafoldDWa6 <- DWAVals6
k <- 5
group6p <- kfold(pafoldDWp6, k)
group6a <- kfold(pafoldDWa6, k)

eDWGLM6 <- list()
pDWGLM6 <- list()
for (i in 1:k) {
  ptrainDW6 <- DWp6[group6p !=i,]
  atrainDW6 <- DWa6[group6a !=i,]
  ptestDW6 <- DWp6[group6p == i,]
  atestDW6 <- DWa6[group6a == i,]
  pa6 <- c(rep(1, nrow(ptrainDW6)), rep(0, nrow(atrainDW6)))
  sdmdatatrainDW6 <- data.frame(cbind(pa6, rbind(ptrainDW6[2:4], atrainDW6[2:4])))
  pmDW6 <- glm(pa6 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
               * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDW6, 
               family=binomial)
  pmDW6 <- step(pmDW6)
  pDWGLM6[[i]] <- predict(Predictors800m, pmDW6, type="response")
  eDWGLM6[[i]] <- evaluate(pmDW6, p=ptestDW6, a=atestDW6)
}
pDWGLMm6 <- stack(pDWGLM6)
PredictionsDWGLM6 <- calc(pDWGLMm6, fun=mean)
aucDWGLM6 <- sapply( eDWGLM6, function(x){slot(x, 'auc')} )
mean(aucDWGLM6)
response(pmDW6)
plot(PredictionsDWGLM6)
DWGLM6TSS <- sapply(eDWGLM6, threshold)
DWGLM6TSS <- unlist(DWGLM6TSS)
DWGLM6TSS <- DWGLM6TSS[seq(2, length(DWGLM6TSS), 6)]
DWGLM6TSS <- mean(DWGLM6TSS)

### GLM 1600
DWPVals7 <- extract(Predictors1600m, DWPresent)
DWAVals7 <- extract(Predictors1600m, DWAbsent)

### Format data
DWPVals7 <- as.data.frame(DWPVals7)
DWAVals7 <- as.data.frame(DWAVals7)
DWp7 <- c(rep(1, nrow(DWPVals7)))
DWp7 <- data.frame(cbind(DWp7, rbind(DWPVals7)))
DWa7 <- c(rep(0, nrow(DWAVals7)))
DWa7 <- data.frame(cbind(DWa7, rbind(DWAVals7)))

pafoldDWp7 <- DWPVals7
pafoldDWa7 <- DWAVals7
k <- 5
group7p <- kfold(pafoldDWp7, k)
group7a <- kfold(pafoldDWa7, k)

eDWGLM7 <- list()
pDWGLM7 <- list()
for (i in 1:k) {
  ptrainDW7 <- DWp7[group7p !=i,]
  atrainDW7 <- DWa7[group7a !=i,]
  ptestDW7 <- DWp7[group7p == i,]
  atestDW7 <- DWa7[group7a == i,]
  pa7 <- c(rep(1, nrow(ptrainDW7)), rep(0, nrow(atrainDW7)))
  sdmdatatrainDW7 <- data.frame(cbind(pa7, rbind(ptrainDW7[2:4], atrainDW7[2:4])))
  pmDW7 <- glm(pa7 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
               * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDW7, 
               family=binomial)
  pmDW7 <- step(pmDW7)
  pDWGLM7[[i]] <- predict(Predictors1600m, pmDW7, type="response")
  eDWGLM7[[i]] <- evaluate(pmDW7, p=ptestDW7, a=atestDW7)
}
pDWGLMm7 <- stack(pDWGLM7)
PredictionsDWGLM7 <- calc(pDWGLMm7, fun=mean)
aucDWGLM7 <- sapply( eDWGLM7, function(x){slot(x, 'auc')} )
mean(aucDWGLM7)
response(pmDW7)
plot(PredictionsDWGLM7)
DWGLM7TSS <- sapply(eDWGLM7, threshold)
DWGLM7TSS <- unlist(DWGLM7TSS)
DWGLM7TSS <- DWGLM7TSS[seq(2, length(DWGLM7TSS), 6)]
DWGLM7TSS <- mean(DWGLM7TSS)

mult.test <- glm(pa7 ~ Percent_Forest * Elevation * Aspect, data = sdmdatatrainDW7 )
## Average AUC
mean(aucDWGLM)
mean(aucDWGLM2)
mean(aucDWGLM3)
mean(aucDWGLM4)
mean(aucDWGLM5)
mean(aucDWGLM6)
mean(aucDWGLM7)

## Variable responses
response(pmDW, main="25 m")
response(pmDW2, main="50 m")
response(pmDW3, main="100 m")
response(pmDW4, main="200 m")
response(pmDW5, main="400 m")
response(pmDW6, main="800 m")
response(pmDW7, main="1600 m")


## Average predictions
PredictionsDWGLM25m <- extract(PredictionsDWGLM, SuitPoints)
PredictionsDWGLM50m <- extract(PredictionsDWGLM2, SuitPoints)
PredictionsDWGLM100m <- extract(PredictionsDWGLM3, SuitPoints)
PredictionsDWGLM200m <- extract(PredictionsDWGLM4, SuitPoints)
PredictionsDWGLM400m <- extract(PredictionsDWGLM5, SuitPoints)
PredictionsDWGLM800m <- extract(PredictionsDWGLM6, SuitPoints)
PredictionsDWGLM1600m <- extract(PredictionsDWGLM7, SuitPoints)
# TRUTH!
PredictionsDWGLMtrue <- extract(HabSpec$suitab.raster, SuitPoints)

par(mfrow=c(4,2), pty="s", mai=c(.4,.2,.2,.2))

cor(PredictionsDWGLM25m, PredictionsDWtrue)
plot(PredictionsDWGLM25m,PredictionsDWtrue, xlab="Prediction", ylab="Truth", main=paste0("25 m (Pearson's r= ", formatC(cor(PredictionsDWGLM25m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWGLM50m, PredictionsDWtrue)
plot(PredictionsDWGLM50m, PredictionsDWtrue, xlab="Prediction", ylab="Truth", main=paste0("50 m (Pearson's r= ", formatC(cor(PredictionsDWGLM50m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWGLM100m, PredictionsDWtrue)
plot(PredictionsDWGLM100m, PredictionsDWtrue, xlab="Prediction", ylab="Truth",main=paste0("100 m (Pearson's r= ", formatC(cor(PredictionsDWGLM100m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
(PredictionsDWGLM200m, PredictionsDWtrue)
plot(PredictionsDWGLM200m, PredictionsDWtrue, xlab="Prediction", ylab="Truth", main=paste0("200 m (Pearson's r= ", formatC(cor(PredictionsDWGLM200m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWGLM400m, PredictionsDWtrue)
plot(PredictionsDWGLM400m, PredictionsDWtrue, xlab="Prediction", ylab="Truth", main=paste0("400 m (Pearson's r= ", formatC(cor(PredictionsDWGLM400m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWGLM800m, PredictionsDWtrue)
plot(PredictionsDWGLM800m, PredictionsDWtrue, xlab="Prediction", ylab="Truth", main=paste0("800 m (Pearson's r= ", formatC(cor(PredictionsDWGLM800m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWGLM1600m, PredictionsDWtrue)
plot(PredictionsDWGLM1600m,PredictionsDWtrue, xlab="Prediction", ylab="Truth", main=paste0("1600 m (Pearson's r= ", formatC(cor(PredictionsDWGLM1600m, PredictionsDWtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")

# Mean Predictions
mean(PredictionsDWGLM25m)
mean(PredictionsDWGLM50m)
mean(PredictionsDWGLM100m)
mean(PredictionsDWGLM200m)
mean(PredictionsDWGLM400m)
mean(PredictionsDWGLM800m)
mean(PredictionsDWGLM1600m)
# TRUTH!
mean(PredictionsDWGLMtrue)

sd(PredictionsDWGLM25m)
sd(PredictionsDWGLM50m)
sd(PredictionsDWGLM100m)
sd(PredictionsDWGLM200m)
sd(PredictionsDWGLM400m)
sd(PredictionsDWGLM800m)
sd(PredictionsDWGLM1600m)
# TRUTH!
sd(PredictionsDWGLMtrue)

## Coarse Grain Species

## 25m
DWFPVals <- extract(PredictorsVS, DWFPresent)
DWFAVals <- extract(PredictorsVS, DWFAbsent)

### Format data
DWFPVals <- as.data.frame(DWFPVals)
DWFAVals <- as.data.frame(DWFAVals)
DWFp <- c(rep(1, nrow(DWFPVals)))
DWFp <- data.frame(cbind(DWFp, rbind(DWFPVals)))
DWFa <- c(rep(0, nrow(DWFAVals)))
DWFa <- data.frame(cbind(DWFa, rbind(DWFAVals)))

pafoldDWFp <- DWFPVals
pafoldDWFa <- DWFAVals
k <- 5
groupp <- kfold(pafoldDWFp, k)
groupa <- kfold(pafoldDWFa, k)

eDWFGLM <- list()
pDWFGLM <- list()
for (i in 1:k) {
  ptrainDWF <- DWFp[groupp !=i,]
  atrainDWF <- DWFa[groupa !=i,]
  ptestDWF <- DWFp[groupp == i,]
  atestDWF <- DWFa[groupa == i,]
  pa <- c(rep(1, nrow(ptrainDWF)), rep(0, nrow(atrainDWF)))
  sdmdatatrainDWF <- data.frame(cbind(pa, rbind(ptrainDWF[2:4], atrainDWF[2:4])))
  pmDWF <- glm(pa ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
               * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDWF, 
               family=binomial)
  pmDWF <- step(pmDWF)
  pDWFGLM[[i]] <- predict(PredictorsVS, pmDWF, type="response")
  eDWFGLM[[i]] <- evaluate(pmDWF, p=ptestDWF, a=atestDWF)
}
pDWFGLMm <- stack(pDWFGLM)
PredictionsDWFGLM <- calc(pDWFGLMm, fun=mean)
aucDWFGLM <- sapply( eDWFGLM, function(x){slot(x, 'auc')} )
mean(aucDWFGLM)
response(pmDWF)
plot(PredictionsDWFGLM)

DWFGLMTSS <- sapply(eDWFGLM, threshold)
DWFGLMTSS <- unlist(DWFGLMTSS)
DWFGLMTSS <- DWFGLMTSS[seq(2, length(DWFGLMTSS), 6)]
DWFGLMTSS <- mean(DWFGLMTSS)

### GLM 50m
DWFPVals2 <- extract(Predictors50m, DWFPresent)
DWFAVals2 <- extract(Predictors50m, DWFAbsent)

### Format data
DWFPVals2 <- as.data.frame(DWFPVals2)
DWFAVals2 <- as.data.frame(DWFAVals2)
DWFp2 <- c(rep(1, nrow(DWFPVals2)))
DWFp2 <- data.frame(cbind(DWFp2, rbind(DWFPVals2)))
DWFa2 <- c(rep(0, nrow(DWFAVals2)))
DWFa2 <- data.frame(cbind(DWFa2, rbind(DWFAVals2)))

pafoldDWFp2 <- DWFPVals2
pafoldDWFa2 <- DWFAVals2
k <- 5
group2p <- kfold(pafoldDWFp2, k)
group2a <- kfold(pafoldDWFa2, k)

eDWFGLM2 <- list()
pDWFGLM2 <- list()
for (i in 1:k) {
  ptrainDWF2 <- DWFp2[group2p !=i,]
  atrainDWF2 <- DWFa2[group2a !=i,]
  ptestDWF2 <- DWFp2[group2p == i,]
  atestDWF2 <- DWFa2[group2a == i,]
  pa2 <- c(rep(1, nrow(ptrainDWF2)), rep(0, nrow(atrainDWF2)))
  sdmdatatrainDWF2 <- data.frame(cbind(pa2, rbind(ptrainDWF2[2:4], atrainDWF2[2:4])))
  pmDWF2 <- glm(pa2 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
                * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDWF2, 
                family=binomial)
  pmDWF2 <- step(pmDWF2)
  pDWFGLM2[[i]] <- predict(Predictors50m, pmDWF2, type="response")
  eDWFGLM2[[i]] <- evaluate(pmDWF2, p=ptestDWF2, a=atestDWF2)
}
pDWFGLMm2 <- stack(pDWFGLM2)
PredictionsDWFGLM2 <- calc(pDWFGLMm2, fun=mean)
aucDWFGLM2 <- sapply( eDWFGLM2, function(x){slot(x, 'auc')} )
mean(aucDWFGLM2)
response(pmDWF2)
plot(PredictionsDWFGLM2)
DWFGLM2TSS <- sapply(eDWFGLM2, threshold)
DWFGLM2TSS <- unlist(DWFGLM2TSS)
DWFGLM2TSS <- DWFGLM2TSS[seq(2, length(DWFGLM2TSS), 6)]
DWFGLM2TSS <- mean(DWFGLM2TSS)

### GLM 100m
DWFPVals3 <- extract(Predictors100m, DWFPresent)
DWFAVals3 <- extract(Predictors100m, DWFAbsent)

### Format data
DWFPVals3 <- as.data.frame(DWFPVals3)
DWFAVals3 <- as.data.frame(DWFAVals3)
DWFp3 <- c(rep(1, nrow(DWFPVals3)))
DWFp3 <- data.frame(cbind(DWFp3, rbind(DWFPVals3)))
DWFa3 <- c(rep(0, nrow(DWFAVals3)))
DWFa3 <- data.frame(cbind(DWFa3, rbind(DWFAVals3)))

pafoldDWFp3 <- DWFPVals3
pafoldDWFa3 <- DWFAVals3
k <- 5
group3p <- kfold(pafoldDWFp3, k)
group3a <- kfold(pafoldDWFa3, k)

eDWFGLM3 <- list()
pDWFGLM3 <- list()
for (i in 1:k) {
  ptrainDWF3 <- DWFp3[group3p !=i,]
  atrainDWF3 <- DWFa3[group3a !=i,]
  ptestDWF3 <- DWFp3[group3p == i,]
  atestDWF3 <- DWFa3[group3a == i,]
  pa3 <- c(rep(1, nrow(ptrainDWF3)), rep(0, nrow(atrainDWF3)))
  sdmdatatrainDWF3 <- data.frame(cbind(pa3, rbind(ptrainDWF3[2:4], atrainDWF3[2:4])))
  pmDWF3 <- glm(pa3 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
                * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDWF3, 
                family=binomial)
  pmDWF3 <- step(pmDWF3)  
  pDWFGLM3[[i]] <- predict(Predictors100m, pmDWF3, type="response")
  eDWFGLM3[[i]] <- evaluate(pmDWF3, p=ptestDWF3, a=atestDWF3)
}
pDWFGLMm3 <- stack(pDWFGLM3)
PredictionsDWFGLM3 <- calc(pDWFGLMm3, fun=mean)
aucDWFGLM3 <- sapply( eDWFGLM3, function(x){slot(x, 'auc')} )
mean(aucDWFGLM3)
response(pmDWF3)
plot(PredictionsDWFGLM3)
DWFGLM3TSS <- sapply(eDWFGLM3, threshold)
DWFGLM3TSS <- unlist(DWFGLM3TSS)
DWFGLM3TSS <- DWFGLM3TSS[seq(2, length(DWFGLM3TSS), 6)]
DWFGLM3TSS <- mean(DWFGLM3TSS)

### GLM 200m
DWFPVals4 <- extract(Predictors200m, DWFPresent)
DWFAVals4 <- extract(Predictors200m, DWFAbsent)

### Format data
DWFPVals4 <- as.data.frame(DWFPVals4)
DWFAVals4 <- as.data.frame(DWFAVals4)
DWFp4 <- c(rep(1, nrow(DWFPVals4)))
DWFp4 <- data.frame(cbind(DWFp4, rbind(DWFPVals4)))
DWFa4 <- c(rep(0, nrow(DWFAVals4)))
DWFa4 <- data.frame(cbind(DWFa4, rbind(DWFAVals4)))

pafoldDWFp4 <- DWFPVals4
pafoldDWFa4 <- DWFAVals4
k <- 5
group4p <- kfold(pafoldDWFp4, k)
group4a <- kfold(pafoldDWFa4, k)

eDWFGLM4 <- list()
pDWFGLM4 <- list()
for (i in 1:k) {
  ptrainDWF4 <- DWFp4[group4p !=i,]
  atrainDWF4 <- DWFa4[group4a !=i,]
  ptestDWF4 <- DWFp4[group4p == i,]
  atestDWF4 <- DWFa4[group4a == i,]
  pa4 <- c(rep(1, nrow(ptrainDWF4)), rep(0, nrow(atrainDWF4)))
  sdmdatatrainDWF4 <- data.frame(cbind(pa4, rbind(ptrainDWF4[2:4], atrainDWF4[2:4])))
  pmDWF4 <- glm(pa4 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
                * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDWF4, 
                family=binomial)
  pmDWF4 <- step(pmDWF4) 
  pDWFGLM4[[i]] <- predict(Predictors200m, pmDWF4, type="response")
  eDWFGLM4[[i]] <- evaluate(pmDWF4, p=ptestDWF4, a=atestDWF4)
}
pDWFGLMm4 <- stack(pDWFGLM4)
PredictionsDWFGLM4 <- calc(pDWFGLMm4, fun=mean)
aucDWFGLM4 <- sapply( eDWFGLM4, function(x){slot(x, 'auc')} )
mean(aucDWFGLM4)
response(pmDWF4)
plot(PredictionsDWFGLM4)

DWFGLM4TSS <- sapply(eDWFGLM4, threshold)
DWFGLM4TSS <- unlist(DWFGLM4TSS)
DWFGLM4TSS <- DWFGLM4TSS[seq(2, length(DWFGLM4TSS), 6)]
DWFGLM4TSS <- mean(DWFGLM4TSS)

### GLM 400m
DWFPVals5 <- extract(Predictors400m, DWFPresent)
DWFAVals5 <- extract(Predictors400m, DWFAbsent)

### Format data
DWFPVals5 <- as.data.frame(DWFPVals5)
DWFAVals5 <- as.data.frame(DWFAVals5)
DWFp5 <- c(rep(1, nrow(DWFPVals5)))
DWFp5 <- data.frame(cbind(DWFp5, rbind(DWFPVals5)))
DWFa5 <- c(rep(0, nrow(DWFAVals5)))
DWFa5 <- data.frame(cbind(DWFa5, rbind(DWFAVals5)))

pafoldDWFp5 <- DWFPVals5
pafoldDWFa5 <- DWFAVals5
k <- 5
group5p <- kfold(pafoldDWFp5, k)
group5a <- kfold(pafoldDWFa5, k)

eDWFGLM5 <- list()
pDWFGLM5 <- list()
for (i in 1:k) {
  ptrainDWF5 <- DWFp5[group5p !=i,]
  atrainDWF5 <- DWFa5[group5a !=i,]
  ptestDWF5 <- DWFp5[group5p == i,]
  atestDWF5 <- DWFa5[group5a == i,]
  pa5 <- c(rep(1, nrow(ptrainDWF5)), rep(0, nrow(atrainDWF5)))
  sdmdatatrainDWF5 <- data.frame(cbind(pa5, rbind(ptrainDWF5[2:4], atrainDWF5[2:4])))
  pmDWF5 <- glm(pa5 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
                * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDWF5, 
                family=binomial)
  pmDWF5 <- step(pmDWF5)
  pDWFGLM5[[i]] <- predict(Predictors400m, pmDWF5, type="response")
  eDWFGLM5[[i]] <- evaluate(pmDWF5, p=ptestDWF5, a=atestDWF5)
}
pDWFGLMm5 <- stack(pDWFGLM5)
PredictionsDWFGLM5 <- calc(pDWFGLMm5, fun=mean)
aucDWFGLM5 <- sapply( eDWFGLM5, function(x){slot(x, 'auc')} )
mean(aucDWFGLM5)
response(pmDWF5)
plot(PredictionsDWFGLM5)
DWFGLM5TSS <- sapply(eDWFGLM5, threshold)
DWFGLM5TSS <- unlist(DWFGLM5TSS)
DWFGLM5TSS <- DWFGLM5TSS[seq(2, length(DWFGLM5TSS), 6)]
DWFGLM5TSS <- mean(DWFGLM5TSS)

### GLM 800
DWFPVals6 <- extract(Predictors800m, DWFPresent)
DWFAVals6 <- extract(Predictors800m, DWFAbsent)

### Format data
DWFPVals6 <- as.data.frame(DWFPVals6)
DWFAVals6 <- as.data.frame(DWFAVals6)
DWFp6 <- c(rep(1, nrow(DWFPVals6)))
DWFp6 <- data.frame(cbind(DWFp6, rbind(DWFPVals6)))
DWFa6 <- c(rep(0, nrow(DWFAVals6)))
DWFa6 <- data.frame(cbind(DWFa6, rbind(DWFAVals6)))

pafoldDWFp6 <- DWFPVals6
pafoldDWFa6 <- DWFAVals6
k <- 5
group6p <- kfold(pafoldDWFp6, k)
group6a <- kfold(pafoldDWFa6, k)

eDWFGLM6 <- list()
pDWFGLM6 <- list()
for (i in 1:k) {
  ptrainDWF6 <- DWFp6[group6p !=i,]
  atrainDWF6 <- DWFa6[group6a !=i,]
  ptestDWF6 <- DWFp6[group6p == i,]
  atestDWF6 <- DWFa6[group6a == i,]
  pa6 <- c(rep(1, nrow(ptrainDWF6)), rep(0, nrow(atrainDWF6)))
  sdmdatatrainDWF6 <- data.frame(cbind(pa6, rbind(ptrainDWF6[2:4], atrainDWF6[2:4])))
  pmDWF6 <- glm(pa6 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
                * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDWF6, 
                family=binomial)
  pmDWF6 <- step(pmDWF6)
  pDWFGLM6[[i]] <- predict(Predictors800m, pmDWF6, type="response")
  eDWFGLM6[[i]] <- evaluate(pmDWF6, p=ptestDWF6, a=atestDWF6)
}
pDWFGLMm6 <- stack(pDWFGLM6)
PredictionsDWFGLM6 <- calc(pDWFGLMm6, fun=mean)
aucDWFGLM6 <- sapply( eDWFGLM6, function(x){slot(x, 'auc')} )
mean(aucDWFGLM6)
response(pmDWF6)
plot(PredictionsDWFGLM6)
DWFGLM6TSS <- sapply(eDWFGLM6, threshold)
DWFGLM6TSS <- unlist(DWFGLM6TSS)
DWFGLM6TSS <- DWFGLM6TSS[seq(2, length(DWFGLM6TSS), 6)]
DWFGLM6TSS <- mean(DWFGLM6TSS)

### GLM 1600
DWFPVals7 <- extract(Predictors1600m, DWFPresent)
DWFAVals7 <- extract(Predictors1600m, DWFAbsent)

### Format data
DWFPVals7 <- as.data.frame(DWFPVals7)
DWFAVals7 <- as.data.frame(DWFAVals7)
DWFp7 <- c(rep(1, nrow(DWFPVals7)))
DWFp7 <- data.frame(cbind(DWFp7, rbind(DWFPVals7)))
DWFa7 <- c(rep(0, nrow(DWFAVals7)))
DWFa7 <- data.frame(cbind(DWFa7, rbind(DWFAVals7)))

pafoldDWFp7 <- DWFPVals7
pafoldDWFa7 <- DWFAVals7
k <- 5
group7p <- kfold(pafoldDWFp7, k)
group7a <- kfold(pafoldDWFa7, k)

eDWFGLM7 <- list()
pDWFGLM7 <- list()
for (i in 1:k) {
  ptrainDWF7 <- DWFp7[group7p !=i,]
  atrainDWF7 <- DWFa7[group7a !=i,]
  ptestDWF7 <- DWFp7[group7p == i,]
  atestDWF7 <- DWFa7[group7a == i,]
  pa7 <- c(rep(1, nrow(ptrainDWF7)), rep(0, nrow(atrainDWF7)))
  sdmdatatrainDWF7 <- data.frame(cbind(pa7, rbind(ptrainDWF7[2:4], atrainDWF7[2:4])))
  pmDWF7 <- glm(pa7 ~ Percent_Forest^3 * Percent_Forest^2 * Percent_Forest
                * Elevation^3 * Elevation^2 * Elevation * Aspect^3 * Aspect^2 * Aspect, data=sdmdatatrainDWF7, 
                family=binomial)
  pmDWF7 <- step(pmDWF7)
  pDWFGLM7[[i]] <- predict(Predictors1600m, pmDWF7, type="response")
  eDWFGLM7[[i]] <- evaluate(pmDWF7, p=ptestDWF7, a=atestDWF7)
}
pDWFGLMm7 <- stack(pDWFGLM7)
PredictionsDWFGLM7 <- calc(pDWFGLMm7, fun=mean)
aucDWFGLM7 <- sapply( eDWFGLM7, function(x){slot(x, 'auc')} )
mean(aucDWFGLM7)
response(pmDWF7)
plot(PredictionsDWFGLM7)

DWFGLM7TSS <- sapply(eDWFGLM7, threshold)
DWFGLM7TSS <- unlist(DWFGLM7TSS)
DWFGLM7TSS <- DWFGLM7TSS[seq(2, length(DWFGLM7TSS), 6)]
DWFGLM7TSS <- mean(DWFGLM7TSS)

### Variable responses
response(pmDWF, main="25 m")
response(pmDWF2, main="50 m")
response(pmDWF3, main="100 m")
response(pmDWF4, main="200 m")
response(pmDWF5, main="400 m")
response(pmDWF6, main="800 m")
response(pmDWF7, main="1600 m")
## Average AUC
mean(aucDWFGLM)
mean(aucDWFGLM2)
mean(aucDWFGLM3)
mean(aucDWFGLM4)
mean(aucDWFGLM5)
mean(aucDWFGLM6)
mean(aucDWFGLM7)

## Correlations
par(mfrow=c(4,2), pty="s", mai=c(.4,.2,.2,.2))

cor(PredictionsDWFGLM25m, PredictionsDWFtrue)
plot(PredictionsDWFGLM25m,PredictionsDWFtrue, xlab="Prediction", ylab="Truth", main=paste0("25 m (Pearson's r= ", formatC(cor(PredictionsDWFGLM25m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFGLM50m, PredictionsDWFtrue)
plot(PredictionsDWFGLM50m, PredictionsDWFtrue, xlab="Prediction", ylab="Truth", main=paste0("50 m (Pearson's r= ", formatC(cor(PredictionsDWFGLM50m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFGLM100m, PredictionsDWFtrue)
plot(PredictionsDWFGLM100m, PredictionsDWFtrue, xlab="Prediction", ylab="Truth",main=paste0("100 m (Pearson's r= ", formatC(cor(PredictionsDWFGLM100m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
plot(PredictionsDWFGLM200m, PredictionsDWFtrue, xlab="Prediction", ylab="Truth", main=paste0("200 m (Pearson's r= ", formatC(cor(PredictionsDWFGLM200m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFGLM400m, PredictionsDWFtrue)
plot(PredictionsDWFGLM400m, PredictionsDWFtrue, xlab="Prediction", ylab="Truth", main=paste0("400 m (Pearson's r= ", formatC(cor(PredictionsDWFGLM400m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFGLM800m, PredictionsDWFtrue)
plot(PredictionsDWFGLM800m, PredictionsDWFtrue, xlab="Prediction", ylab="Truth", main=paste0("800 m (Pearson's r= ", formatC(cor(PredictionsDWFGLM800m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")
cor(PredictionsDWFGLM1600m, PredictionsDWFtrue)
plot(PredictionsDWFGLM1600m,PredictionsDWFtrue, xlab="Prediction", ylab="Truth", main=paste0("1600 m (Pearson's r= ", formatC(cor(PredictionsDWFGLM1600m, PredictionsDWFtrue), 3, format="f"),")"),mgp=c(2,.5,.2))
abline(a=0,b=1, col="green3")

## Average predictions
PredictionsDWFGLM25m <- extract(PredictionsDWFGLM, SuitPoints)
PredictionsDWFGLM50m <- extract(PredictionsDWFGLM2, SuitPoints)
PredictionsDWFGLM100m <- extract(PredictionsDWFGLM3, SuitPoints)
PredictionsDWFGLM200m <- extract(PredictionsDWFGLM4, SuitPoints)
PredictionsDWFGLM400m <- extract(PredictionsDWFGLM5, SuitPoints)
PredictionsDWFGLM800m <- extract(PredictionsDWFGLM6, SuitPoints)
PredictionsDWFGLM1600m <- extract(PredictionsDWFGLM7, SuitPoints)
# TRUTH!
PredictionsDWFGLMtrue <- extract(DirewolfF$suitab.raster, SuitPoints)

# Mean Predictions
mean(PredictionsDWFGLM25m)
mean(PredictionsDWFGLM50m)
mean(PredictionsDWFGLM100m)
mean(PredictionsDWFGLM200m)
mean(PredictionsDWFGLM400m)
mean(PredictionsDWFGLM800m)
mean(PredictionsDWFGLM1600m)
# TRUTH!
mean(PredictionsDWFGLMtrue)

sd(PredictionsDWFGLM25m)
sd(PredictionsDWFGLM50m)
sd(PredictionsDWFGLM100m)
sd(PredictionsDWFGLM200m)
sd(PredictionsDWFGLM400m)
sd(PredictionsDWFGLM800m)
sd(PredictionsDWFGLM1600m)
# TRUTH!
sd(PredictionsDWFGLMtrue)

### Try to do a multiple comparison of model AUCs
library(PMCMR)
TestResultsMaxEnt <- c(auc1, auc2, auc3, auc4, auc400, auc5, auc6)

###Methods <- c("A","A","A","A","A","B","B","B","B","B","C","C","C","C","C","D","D","D",
"D","D","E","E","E","E","E","F","F","F","F","F","G","G","G","G","G")
##Methods <- as.factor(Methods)
Methodsalt <- c(25,25,25,25,25,50,50,50,50,50,100,100,100,100,100,
                200,200,200,200,200,400,400,400,400,400,800,800,800,800,800,1600,1600,1600,1600,1600)
vanWaerden.test(x=TestResultsMaxEnt, g=Methodsalt)
posthoc.vanWaerden.test(x=TestResultsMaxEnt, g=Methodsalt, p.adjust.method="none")

### MaxEnt Coarse Grain
#TestResultsMaxEntF <- c(aucF, aucF2, aucF3, aucF4, aucF400, aucF5, aucF6)

TestResultsMaxEntF <- c(aucF,aucF2, aucF3, aucF4, aucF400, aucF5, aucF6)
vanWaerden.test(x=TestResultsMaxEntF, g=Methodsalt)
posthoc.vanWaerden.test(x=TestResultsMaxEntF, g=Methodsalt, p.adjust.method="none")

## GLM Fine Grain
TestResultsGLM <- c(aucDWGLM,aucDWGLM2,aucDWGLM3,aucDWGLM4,aucDWGLM5,aucDWGLM6,aucDWGLM7)

vanWaerden.test(x=TestResultsGLM, g=Methodsalt)
posthoc.vanWaerden.test(x=TestResultsGLM, g=Methodsalt, p.adjust.method="none")


TestResultsGLMF <- c(aucDWFGLM,aucDWFGLM2,aucDWFGLM3,aucDWFGLM4,aucDWFGLM5,aucDWFGLM6,aucDWFGLM7)
#AUC_comp <- data.frame(x=TestResults, g=Methods)
vanWaerden.test(x=TestResultsGLMF, g=Methodsalt)
posthoc.vanWaerden.test(x=TestResultsGLMF, g=Methodsalt, p.adjust.method="none")


#kruskal.test(x ~ g, data=AUC_comp)
#posthoc.kruskal.nemenyi.test(x=TestResults, g=Methods, dist="Tukey")

## Get all the thresholds outa the evaluates
## MaxEnt Fine
eDW
eDW2
eDW3
eDW4
eDW400
eDW5
eDW6

## MaxEnt Coarse
eDWF
eDWF2
eDWF3
eDWF4
eDWF400
eDWF5
eDWF6

## GLM Fine
eDWGLM
eDWGLM2
eDWGLM3
eDWGLM4
eDWGLM5
eDWGLM6
eDWGLM7

## GLM Coarse
eDWFGLM
eDWFGLM2
eDWFGLM3
eDWFGLM4
eDWFGLM5
eDWFGLM6
eDWFGLM7

### Rescale suitability to fall between 0 and 1 via: (val - min) / (max - min)
##RescalePDWFGLM7 <- calc(PredictionsDWFGLM7, fun=function(x){(x-minValue(PredictionsDWFGLM7))/(maxValue(PredictionsDWFGLM7)-minValue(PredictionsDWFGLM7))})
##PredictionsDWFGLM7RE <- calc(RescalePDWFGLM7, fun=mean)
##PredictionsDWFGLM1600mRE <- extract(PredictionsDWFGLM7RE, SuitPoints)
##mean(PredictionsDWFGLM1600mRE)
## Check full models:
me.DW4E<- maxent(Predictors200m, DWP2)
#plot(me.DW2)
r.DW4E <- dismo::predict(me.DWF, Predictors200m)
eDW4E <- evaluate(me.DW4E, p=DWP2, a=bg.DW, x=Predictors200m)



# threshold at maximum of the sum of the sensitivity (true positive rate) and specificity (true negative rate)
eDW4E@t[which.max(eDW4E@TPR + eDW4E@TNR)]
plot(eDW4E, 'kappa')
ET <- threshold(eDW4E)
plot(eDW4E@t, eDW4E@TPR+eDW4E@TNR, type="l", lwd=2, cex.lab=1.5)

EDW400t <- sapply( eDW400, function(x){slot(x, 't')} )
write.table(EDW400t[[3]], "threshold_vector.txt")

## Reclassify the raster to fit max TPR+ TNR
## 25m
matre <- c(0, DWTSS, NA, DWTSS, 1, 1)
rclmat <- matrix(matre, ncol=3, byrow=TRUE)
rcME <- reclassify(pDWME, rclmat)
plot(rcME, legend=F, main="25 m")

f <- freq(rcME) 
PixelArea <- 25*25

HabAreaDW <- f*PixelArea/1000000
HabAreaDW

## 50m
matre2 <- c(0, DW2TSS, NA, DW2TSS, 1, 1)
rclmat2 <- matrix(matre2, ncol=3, byrow=TRUE)
rcME2 <- reclassify(pDWME2, rclmat2)
plot(rcME2, legend=F, main="50 m")

f2 <- freq(rcME2) 
PixelArea2 <- 50*50

HabAreaDW2 <- f2*PixelArea2/1000000
HabAreaDW2

## 100m
matre3 <- c(0, DW3TSS, NA, DW3TSS, 1, 1)
rclmat3 <- matrix(matre3, ncol=3, byrow=TRUE)
rcME3 <- reclassify(pDWME3, rclmat3)
plot(rcME3, legend=F, main="100 m")

f3 <- freq(rcME3) 
PixelArea3 <- 100*100

HabAreaDW3 <- f3*PixelArea3/1000000
HabAreaDW3



## 200m
matre4 <- c(0, DW4TSS, NA, DW4TSS, 1, 1)
rclmat4 <- matrix(matre4, ncol=3, byrow=TRUE)
rcME4 <- reclassify(pDWME4, rclmat4)
plot(rcME4, legend=F, main="200 m")

f4 <- freq(rcME4) 
PixelArea4 <- 200*200

HabAreaDW4 <- f4*PixelArea4/1000000
HabAreaDW4

## 400m
matre5 <- c(0, DW400TSS, NA, DW400TSS, 1, 1)
rclmat5 <- matrix(matre5, ncol=3, byrow=TRUE)
rcME5 <- reclassify(pDWME400, rclmat5)
plot(rcME5, legend=F, main="400 m")

f5 <- freq(rcME5) 
PixelArea5 <- 400*400

HabAreaDW5 <- f5*PixelArea5/1000000
HabAreaDW5


## 800m
matre6 <- c(0, DW5TSS, NA, DW5TSS, 1, 1)
rclmat6 <- matrix(matre6, ncol=3, byrow=TRUE)
rcME6 <- reclassify(pDWME5, rclmat6)
plot(rcME6, legend=F, main="800 m")

f6 <- freq(rcME6) 
PixelArea6 <- 800*800

HabAreaDW6 <- f6*PixelArea6/1000000
HabAreaDW6

## 1600m
matre7 <- c(0, DW6TSS, NA, DW6TSS, 1, 1)
rclmat7 <- matrix(matre7, ncol=3, byrow=TRUE)
rcME7 <- reclassify(pDWME6, rclmat7)
plot(rcME7, legend=F, main="1600 m")

f7 <- freq(rcME7) 
PixelArea7 <- 1600*1600

HabAreaDW7 <- f7*PixelArea7/1000000
HabAreaDW7


### Coarse Species MaxEnt 
## 25m
matreMEF <- c(0, DWFTSS, NA, DWFTSS, 1, 1)
rclmatMEF <- matrix(matreMEF, ncol=3, byrow=TRUE)
rcMEF <- reclassify(pDWMEF, rclmatMEF)
plot(rcMEF, legend=F, main="25 m")

fF <- freq(rcMEF) 
PixelArea <- 25*25

HabAreaDWF <- fF*PixelArea/1000000
HabAreaDWF

## 50m
matreMEF2 <- c(0, DWF2TSS, NA, DWF2TSS, 1, 1)
rclmatMEF2 <- matrix(matreMEF2, ncol=3, byrow=TRUE)
rcMEF2 <- reclassify(pDWMEF2, rclmatMEF2)
plot(rcMEF2, legend=F, main="50 m")

fF2 <- freq(rcMEF2) 

HabAreaDWF2 <- fF2*PixelArea2/1000000
HabAreaDWF2

## 100m
matreMEF3 <- c(0, DWF3TSS, NA, DWF3TSS, 1, 1)
rclmatMEF3 <- matrix(matreMEF3, ncol=3, byrow=TRUE)
rcMEF3 <- reclassify(pDWMEF3, rclmatMEF3)
plot(rcMEF3, legend=F, main="100 m")

fF3 <- freq(rcMEF3) 
HabAreaDWF3 <- fF3*PixelArea3/1000000
HabAreaDWF3

## 200 m
matreMEF4 <- c(0, DWF4TSS, NA, DWF4TSS, 1, 1)
rclmatMEF4 <- matrix(matreMEF4, ncol=3, byrow=TRUE)
rcMEF4 <- reclassify(pDWMEF4, rclmatMEF4)
plot(rcMEF4, legend=F, main="200 m")

fF4 <- freq(rcMEF4) 

HabAreaDWF4 <- fF4*PixelArea4/1000000
HabAreaDWF4

## 400 m
matreMEF5 <- c(0, DWF400TSS, NA, DWF400TSS, 1, 1)
rclmatMEF5 <- matrix(matreMEF5, ncol=3, byrow=TRUE)
rcMEF5 <- reclassify(pDWMEF400, rclmatMEF5)
plot(rcMEF5, legend=F, main="400 m")

fF5 <- freq(rcMEF5) 

HabAreaDWF5 <- fF5*PixelArea5/1000000
HabAreaDWF5

## 800 m
matreMEF6 <- c(0, DWF5TSS, NA, DWF5TSS, 1, 1)
rclmatMEF6 <- matrix(matreMEF6, ncol=3, byrow=TRUE)
rcMEF6 <- reclassify(pDWMEF5, rclmatMEF6)
plot(rcMEF6, legend=F, main="800 m")

fF6 <- freq(rcMEF6) 

HabAreaDWF6 <- fF6*PixelArea6/1000000
HabAreaDWF6

## 1600 m
matreMEF7 <- c(0, DWF6TSS, NA, DWF6TSS, 1, 1)
rclmatMEF7 <- matrix(matreMEF7, ncol=3, byrow=TRUE)
rcMEF7 <- reclassify(pDWMEF6, rclmatMEF7)
plot(rcMEF7, legend=F, main="1600 m")

fF7 <- freq(rcMEF7) 

HabAreaDWF7 <- fF7*PixelArea7/1000000
HabAreaDWF7

## GLM Fine grain
## 25m
matreGLM <- c(-0.5, DWGLMTSS, NA, DWGLMTSS, 1, 1)
rclmatGLM <- matrix(matreGLM, ncol=3, byrow=TRUE)
rcGLM <- reclassify(PredictionsDWGLM, rclmatGLM)
plot(rcGLM, legend=F, main="25 m")

fGLM <- freq(rcGLM) 

HabAreaDWGLM <- fGLM*PixelArea/1000000
HabAreaDWGLM

## 50m
matreGLM2 <- c(-0.5, DWGLM2TSS, NA, DWGLM2TSS, 1, 1)
rclmatGLM2 <- matrix(matreGLM2, ncol=3, byrow=TRUE)
rcGLM2 <- reclassify(PredictionsDWGLM2, rclmatGLM2)
plot(rcGLM2, legend=F, main="50 m")

fGLM2 <- freq(rcGLM2) 

HabAreaDWGLM2 <- fGLM2*PixelArea2/1000000
HabAreaDWGLM2

## 100m
matreGLM3 <- c(-0.5, DWGLM3TSS, NA, DWGLM3TSS, 1, 1)
rclmatGLM3 <- matrix(matreGLM3, ncol=3, byrow=TRUE)
rcGLM3 <- reclassify(PredictionsDWGLM3, rclmatGLM3)
plot(rcGLM3, legend=F, main="100 m")

fGLM3 <- freq(rcGLM3) 

HabAreaDWGLM3 <- fGLM3*PixelArea3/1000000
HabAreaDWGLM3

## 200m
matreGLM4 <- c(-0.5, DWGLM4TSS, NA, DWGLM4TSS, 1, 1)
rclmatGLM4 <- matrix(matreGLM4, ncol=3, byrow=TRUE)
rcGLM4 <- reclassify(PredictionsDWGLM4, rclmatGLM4)
plot(rcGLM4, legend=F, main="200 m")

fGLM4 <- freq(rcGLM4) 

HabAreaDWGLM4 <- fGLM4*PixelArea4/1000000
HabAreaDWGLM4
## 400m
matreGLM5 <- c(-0.5, DWGLM5TSS, NA, DWGLM5TSS, 1, 1)
rclmatGLM5 <- matrix(matreGLM5, ncol=3, byrow=TRUE)
rcGLM5 <- reclassify(PredictionsDWGLM5, rclmatGLM5)
plot(rcGLM5, legend=F, main="400 m")

fGLM5 <- freq(rcGLM5) 

HabAreaDWGLM5 <- fGLM5*PixelArea5/1000000
HabAreaDWGLM5

## 800m
matreGLM6 <- c(-0.5, DWGLM6TSS, NA, DWGLM6TSS, 1, 1)
rclmatGLM6 <- matrix(matreGLM6, ncol=3, byrow=TRUE)
rcGLM6 <- reclassify(PredictionsDWGLM6, rclmatGLM6)
plot(rcGLM6, legend=F, main="800 m")

fGLM6 <- freq(rcGLM6) 

HabAreaDWGLM6 <- fGLM6*PixelArea6/1000000
HabAreaDWGLM6

## 1600m
matreGLM7 <- c(-0.5, DWGLM7TSS, NA, DWGLM7TSS, 1, 1)
rclmatGLM7 <- matrix(matreGLM7, ncol=3, byrow=TRUE)
rcGLM7 <- reclassify(PredictionsDWGLM7, rclmatGLM7)
plot(rcGLM7, legend=F, main="1600 m")

fGLM7 <- freq(rcGLM7) 

HabAreaDWGLM7 <- fGLM7*PixelArea7/1000000
HabAreaDWGLM7

### Coarse Species GLM

## 25m
matreGLMF <- c(-0.5, DWFGLMTSS, NA, DWFGLMTSS, 2, 1)
rclmatGLMF <- matrix(matreGLMF, ncol=3, byrow=TRUE)
rcGLMF <- reclassify(PredictionsDWFGLM, rclmatGLMF)
plot(rcGLMF, legend=T, main="25 m")

fGLMF <- freq(rcGLMF) 

HabAreaDWGLMF <- fGLMF*PixelArea/1000000
HabAreaDWGLMF

## 50m
matreGLMF2 <- c(-0.5, DWFGLM2TSS, NA, DWFGLM2TSS, 2, 1)
rclmatGLMF2 <- matrix(matreGLMF2, ncol=3, byrow=TRUE)
rcGLMF2 <- reclassify(PredictionsDWFGLM2, rclmatGLMF2)
plot(rcGLMF2, legend=F, main="50 m")

fGLMF2 <- freq(rcGLMF2) 

HabAreaDWGLMF2 <- fGLMF2*PixelArea2/1000000
HabAreaDWGLMF2

## 100m
matreGLMF3 <- c(-0.5, DWFGLM3TSS, NA, DWFGLM3TSS, 2, 1)
rclmatGLMF3 <- matrix(matreGLMF3, ncol=3, byrow=TRUE)
rcGLMF3 <- reclassify(PredictionsDWFGLM3, rclmatGLMF3)
plot(rcGLMF3, legend=F, main="100 m")

fGLMF3 <- freq(rcGLMF3) 

HabAreaDWGLMF3 <- fGLMF3*PixelArea3/1000000
HabAreaDWGLMF3

## 200m
matreGLMF4 <- c(-0.5, DWFGLM4TSS, NA, DWFGLM4TSS, 2, 1)
rclmatGLMF4 <- matrix(matreGLMF4, ncol=3, byrow=TRUE)
rcGLMF4 <- reclassify(PredictionsDWFGLM4, rclmatGLMF4)
plot(rcGLMF4, legend=F, main="200 m")

fGLMF4 <- freq(rcGLMF4) 

HabAreaDWGLMF4 <- fGLMF4*PixelArea4/1000000
HabAreaDWGLMF4

## 400m
matreGLMF5 <- c(-0.5, DWFGLM5TSS, NA, DWFGLM5TSS, 2, 1)
rclmatGLMF5 <- matrix(matreGLMF5, ncol=3, byrow=TRUE)
rcGLMF5 <- reclassify(PredictionsDWFGLM5, rclmatGLMF5)
plot(rcGLMF5, legend=F, main="400 m")

fGLMF5 <- freq(rcGLMF5) 

HabAreaDWGLMF5 <- fGLMF5*PixelArea5/1000000
HabAreaDWGLMF5

## 800m
matreGLMF6 <- c(-0.5, DWFGLM6TSS, NA, DWFGLM6TSS, 2, 1)
rclmatGLMF6 <- matrix(matreGLMF6, ncol=3, byrow=TRUE)
rcGLMF6 <- reclassify(PredictionsDWFGLM6, rclmatGLMF6)
plot(rcGLMF6, legend=F, main="800 m")

fGLMF6 <- freq(rcGLMF6) 

HabAreaDWGLMF6 <- fGLMF6*PixelArea6/1000000
HabAreaDWGLMF6

## 1600m
matreGLMF7 <- c(-0.5, DWFGLM7TSS, NA, DWFGLM7TSS, 2, 1)
rclmatGLMF7 <- matrix(matreGLMF7, ncol=3, byrow=TRUE)
rcGLMF7 <- reclassify(PredictionsDWFGLM7, rclmatGLMF7)
plot(rcGLMF7, legend=F, main="1600 m")

fGLMF7 <- freq(rcGLMF7) 

HabAreaDWGLMF7 <- fGLMF7*PixelArea7/1000000
HabAreaDWGLMF7

## True Habitat Areas
## 25m Species
fTruthDW <- freq(paDW$pa.raster)

HabAreaDWTruth <- fTruthDW*PixelArea/1000000
HabAreaDWTruth

### 200m species
fTruthDWF <- freq(paDWF$pa.raster)

HabAreaDWFTruth <- fTruthDWF*PixelArea4/1000000
HabAreaDWFTruth

## Areas in a row
HabAreaDW
HabAreaDW2
HabAreaDW3
HabAreaDW4
HabAreaDW5
HabAreaDW6
HabAreaDW7

HabAreaDWF
HabAreaDWF2
HabAreaDWF3
HabAreaDWF4
HabAreaDWF5
HabAreaDWF6
HabAreaDWF7

HabAreaDWGLM
HabAreaDWGLM2
HabAreaDWGLM3
HabAreaDWGLM4
HabAreaDWGLM5
HabAreaDWGLM6
HabAreaDWGLM7

HabAreaDWGLMF
HabAreaDWGLMF2
HabAreaDWGLMF3
HabAreaDWGLMF4
HabAreaDWGLMF5
HabAreaDWGLMF6
HabAreaDWGLMF7

mai = c(.1, 0.1, 0.1, 0.1))
## Plot it all
par(mar = c(0.01, 0.01, 2, 2))

par(mfrow=c(4,2))
,mai = c(.05, .05, .05, .05))
par(mfrow = c(4,2), mar = c(0,0,.5,.5) + 0.1)
par(mfrow = c(4, 2),     # 2x2 layout
    # two rows of text at the outer left and bottom margin
    mar = c(1, 1, 0, 0)) # space for one row of text at ticks and to separate plot    # axis label at 2 rows distance, tick labels at 1 row

plot(paDW$pa.raster, main="Truth", legend=FALSE,  xaxt='n', yaxt='n', cex.main=2)
plot(rcME, legend=F, main="25 m",  xaxt='n', yaxt='n', cex.main=2)
plot(rcME2, legend=F, main="50 m",  xaxt='n', yaxt='n', cex.main=2)
plot(rcME3, legend=F, main="100 m",  xaxt='n', yaxt='n', cex.main=2)
plot(rcME4, legend=F, main="200 m",  xaxt='n', yaxt='n', cex.main=2)
plot(rcME5, legend=F, main="400 m",  xaxt='n', yaxt='n', cex.main=2)
plot(rcME6, legend=F, main="800 m",  xaxt='n', yaxt='n', cex.main=2)
plot(rcME7, legend=F, main="1600 m",  xaxt='n', yaxt='n', cex.main=2)
