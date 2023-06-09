#setwd("H:/BigLake/Paper/FigsSI1and2")
library(rgdal)
library(raster)
library(knitr)
library(rspatial)

region <- readOGR(dsn="H:/BigLake/Paper/FigsSI1and2", layer="StudyRegionSmallDissolve")
lakes <- readOGR(dsn="H:/BigLake/Paper/FigsSI1and2", layer="AllLakes")
secchilakes <- readOGR(dsn="H:/BigLake/Paper/FigsSI1and2", layer="SecchiLakes")
maxdepthlakes <- readOGR(dsn="H:/BigLake/Paper/FigsSI1and2", layer="MaxDepthLakes")

#### Plot lakes #########################
# Plot all lakes
plot(region, col='light blue')
points(lakes, col='red', cex=.5, pch='+')

# Create raster of region 
r <- raster(region)
res(r) <- 10000
r <- rasterize(region, r)
plot(r)
quads <- as(r, 'SpatialPolygons')

# Plot all lakes on raster of region
plot(quads, add=TRUE)
points(lakes, col='red', cex=.5)

# Plot just Secchi lakes
plot(quads, add=TRUE)
points(secchilakes, col='red', cex=.5)

# Plot just Secchi lakes
plot(quads, add=TRUE)
points(maxdepthlakes, col='red', cex=.25)

#################################################################
# Plot raster count for all lakes
nlk <- rasterize(coordinates(lakes), r, fun='count', background=0)
plot(nlk)
plot(region, add=TRUE)

# Plot raster count for all lakes exclusively in region
nlakes <- mask(nlk, r)
plot(nlakes)
plot(region, add=TRUE)

## Plot frequency of lakes per raster cell
fall <- freq(nlakes, useNA='no')
plot(fall, pch=20)

## Number of raster cells
quadratsall <- sum(fall[,2])
## Number of lakes in region
casesall <- sum(fall[,1]*fall[,2])
## Average number of lakes per cell
muall <- casesall / quadratsall
muall


ffall <- data.frame(fall)
colnames(ffall) <- c('K', 'X')
ffall$Kmu <- ffall$K - muall
ffall$Kmu2 <- ffall$Kmu^2
ffall$XKmu2 <- ffall$Kmu2 * ffall$X
## Variance in the number of lakes per cell
s2all <- sum(ffall$XKmu2)/(sum(ffall$X)-1)
s2all
## Coefficient of variation in the number of lakes per cell
VMRall <- s2all / muall
VMRall

#################################################################
# Plot raster count for Secchi lakes
nslk <- rasterize(coordinates(secchilakes), r, fun='count', background=0)
plot(nslk)
plot(region, add=TRUE)

# Plot raster count for all Secchi lakes exclusively in region
nsecchilakes <- mask(nslk, r)
plot(nsecchilakes)
plot(region, add=TRUE)

## Plot frequency of Secchi lakes per raster cell
fsecchi <- freq(nsecchilakes, useNA='no')
plot(fsecchi, pch=20)

## Number of raster cells
quadratsecchi <- sum(fsecchi[,2])
## Number of Secchi lakes in region
casessecchi <- sum(fsecchi[,1]*fsecchi[,2])
## Average number of Secchi lakes per cell
musecchi <- casessecchi / quadratsecchi
musecchi


ffsecchi <- data.frame(fsecchi)
colnames(ffsecchi) <- c('K', 'X')
ffsecchi$Kmu <- ffsecchi$K - musecchi
ffsecchi$Kmu2 <- ffsecchi$Kmu^2
ffsecchi$XKmu2 <- ffsecchi$Kmu2 * ffsecchi$X
## Variance in the number of Secchi lakes per cell
s2secchi <- sum(ffsecchi$XKmu2)/(sum(ffsecchi$X)-1)
s2secchi
## Coefficient of variation in the number of Secchi lakes per cell
VMRsecchi <- s2secchi / musecchi
VMRsecchi


#################################################################
# Plot raster count for Max Depth lakes
ndlk <- rasterize(coordinates(maxdepthlakes), r, fun='count', background=0)
plot(ndlk)
plot(region, add=TRUE)

# Plot raster count for all depth lakes exclusively in region
ndepthlakes <- mask(ndlk, r)
plot(ndepthlakes)
plot(region, add=TRUE)

## Plot frequency of depth lakes per raster cell
fdepth <- freq(ndepthlakes, useNA='no')
plot(fdepth, pch=20)

## Number of raster cells
quadratdepth <- sum(fdepth[,2])
## Number of depth lakes in region
casesdepth <- sum(fdepth[,1]*fdepth[,2])
## Average number of depth lakes per cell
mudepth <- casesdepth / quadratdepth
mudepth

ffdepth <- data.frame(fdepth)
colnames(ffdepth) <- c('K', 'X')
ffdepth$Kmu <- ffdepth$K - mudepth
ffdepth$Kmu2 <- ffdepth$Kmu^2
ffdepth$XKmu2 <- ffdepth$Kmu2 * ffdepth$X
## Variance in the number of depth lakes per cell
s2depth <- sum(ffdepth$XKmu2)/(sum(ffdepth$X)-1)
s2depth
## Coefficient of variation in the number of depth lakes per cell
VMRdepth <- s2depth / mudepth
VMRdepth

################################################

xyall <- coordinates(lakes)
xysecchi <- coordinates(secchilakes)
xydepth <- coordinates(maxdepthlakes)

StudyArea <- raster::area(region)
denssecchi <- nrow(xysecchi) / StudyArea # Estimate of lambda of Secchi lakes
densdepth <- nrow(xydepth) / StudyArea # Estimate of lambda of depth lakes

### Calculate Secchi lakes distance matrix and min distances ###
dsecchi <- dist(xysecchi)
dmsecchi <- as.matrix(dsecchi)
diag(dmsecchi) <- NA
dsecchimin <- apply(dmsecchi, 1, min, na.rm=TRUE)
meandsecchimin <- mean(dsecchimin)

### Calculate depth lakes distance matrix and min distances ###
ddepth <- dist(xydepth)
dmdepth <- as.matrix(ddepth)
diag(dmdepth) <- NA
ddepthmin <- apply(dmdepth, 1, min, na.rm=TRUE)
meanddepthmin <- mean(ddepthmin)


### Calculate K and L functions for Secchi lakes #######
distanceksecchi <- seq(1, 70000, 100)
Kdsecchi <- sapply(distanceksecchi, function(x) sum(dsecchi < x))
Kdsecchi <- Kdsecchi / (length(Kdsecchi)*denssecchi)
Lsecchi <- sqrt(Kdsecchi/pi)-distanceksecchi  


### Calculate K and L functions for depth lakes #######
distancekdepth <- seq(1, 70000, 100)
Kddepth <- sapply(distancekdepth, function(x) sum(ddepth < x))
Kddepth <- Kddepth / (length(Kddepth)*densdepth)
Ldepth <- sqrt(Kddepth/pi)-distancekdepth 


### Start loop here for Secchi mean distances...
meandsampminsecchi <- matrix(0,1000,1)
distanceksamp <- seq(1, 70000, 100)
Lsampsecchi <- matrix(0,length(distanceksamp),1000)

for (i in 1:1000) 
{

xyallsampindex <- sample(1:nrow(xyall),nrow(xysecchi),replace=F)
xyallsamp <- xyall[xyallsampindex,]

dsamp <- dist(xyallsamp)
dmsamp <- as.matrix(dsamp)
diag(dmsamp) <- NA
dsampmin <- apply(dmsamp, 1, min, na.rm=TRUE)
meandsampminsecchi[i,1] <- mean(dsampmin)

message('Processing ', i)

}

hist(meandsampminsecchi[,1], breaks=seq(4400,5400,by=10),
     border="black", col="grey",xlab="Mean nearest lake distance over a subset of 8,733 lakes (meters)",
     main=NA)
abline(v = meandsecchimin, col="red", lwd=3, lty=2)

meanmeandsampminsecchi <- mean(meandsampminsecchi[,2])
temp <- meandsampminsecchi[,2]-meanmeandsampminsecchi
tempsq <- temp^2
varsecchi <- sum(tempsq/999)
sdsecchi <- var^0.5

#############################

### Start loop here for depth mean distances...
meandsampmindepth <- matrix(0,1000,1)
distanceksamp <- seq(1, 70000, 100)
Lsampdepth <- matrix(0,length(distanceksamp),1000)

for (i in 1:1000) 
{
  
  xyallsampindex <- sample(1:nrow(xyall),nrow(xydepth),replace=F)
  xyallsamp <- xyall[xyallsampindex,]
  
  dsamp <- dist(xyallsamp)
  dmsamp <- as.matrix(dsamp)
  diag(dmsamp) <- NA
  dsampmin <- apply(dmsamp, 1, min, na.rm=TRUE)
  meandsampmindepth[i,1] <- mean(dsampmin)
  
  message('Processing ', i)
  
}

hist(meandsampmindepth[,1], breaks=seq(4200,5400,by=10),
     border="black", col="grey",xlab="Mean nearest lake distance over a subset of 9,371 lakes (meters)",
     main=NA)
abline(v = meanddepthmin, col="red", lwd=3, lty=2)

meanmeandsampmindepth <- mean(meandsampmindepth[,1])
temp <- meandsampmin[,1]-meanmeandsampmin
tempsq <- temp^2
vardepth <- sum(tempsq/999)
sddepth <- var^0.5


####### Calculate K and L functions for sample #######
for (i in 1:1000) 
{
  
xyallsampindex <- sample(1:nrow(xyall),nrow(xysecchi),replace=F)
xyallsamp <- xyall[xyallsampindex,]

dsamp <- dist(xyallsamp)

Kdsamp <- sapply(distanceksamp, function(x) sum(dsamp < x))
Kdsamp <- Kdsamp / (length(Kdsamp)*denssecchi)
Lsampsecchi[1:length(distanceksamp),i] <- sqrt(Kdsamp/pi)-distanceksamp  

### End loop here...

}

Lsampimportsecchi <- read.csv("Lsampsecchi.csv", header=FALSE)
Lsampimportsecchi <- Lsampimportsecchi[,1:500]
plot(distanceksecchi, Lsecchi, type='l', lwd=2)
Lsampimportsecchi <- t(Lsampimportsecchi)

for (i in 3:500) {
lines(distanceksamp, Lsampimportsecchi[i,], lwd=2, col='blue')
}


####### Calculate K and L functions for depth sample #######
for (i in 1:1000) 
{
  
  xyallsampindex <- sample(1:nrow(xyall),nrow(xydepth),replace=F)
  xyallsamp <- xyall[xyallsampindex,]
  
  dsamp <- dist(xyallsamp)
  
  Kdsamp <- sapply(distanceksamp, function(x) sum(dsamp < x))
  Kdsamp <- Kdsamp / (length(Kdsamp)*densdepth)
  Lsampdepth[1:length(distanceksamp),i] <- sqrt(Kdsamp/pi)-distanceksamp  

  message('Processing ', i)  
  ### End loop here...
  
}

Lsampimportdepth <- read.csv("Lsampdepth.csv", header=FALSE)
Lsampimportdepth <- Lsampimportdepth[,1:500]
plot(distancekdepth, Ldepth, type='l', lwd=2)
Lsampimportdepth <- t(Lsampimportdepth)

for (i in 3:500) {
  lines(distanceksamp, Lsampimportdepth[i,], lwd=2, col='blue')
}

