

rm(list=ls())

library(raster)
library(wsyn)
library(rgdal)
library(maptools)
library(RColorBrewer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load datasets
ndvi.raw <- stack("../Data/compositeMODIS_ALL_20220713.tif")
guan <- readOGR("../Data/guanacaste.shp")
groundtruth <- read.csv("../Data/forestPlots_ACG.csv")
image.dates <- read.csv("../Data/ImageMetadata_20220714.csv")[,1]

image.years <- as.numeric(substr(image.dates, 1, 4))
years <- unique(image.years)

image.dates <- as.Date(image.dates, format=c("%Y_%m_%d"))


#ndvi.raw <- projectRaster(ndvi.raw, crs=CRS("+init=EPSG:4326"))

## assign and match projections
prj <- "+proj=lcc +datum=NAD27 +ellps=clrk66 +lat_0=10.4667 +lat_1=9.9333 +lat_2=11.0000 +lon_0=-84.3333 +x_0=500000.000 + y_0=271820.522 +units=m"
crs(guan) <- CRS(prj) 
#guan.prj <- spTransform(guan, CRS("+init=EPSG:4326"))
guan.prj <-spTransform(guan, CRS(proj4string(ndvi.raw)), use_ob_tran=TRUE)

#writeOGR(guan.prj, "../Data/", "guan_prj", driver="ESRI Shapefile")

# quartz()
# plot(ndvi.raw, 1, xlim=c(-86,-85), ylim=c(10.5, 11.5))
# plot(guan.prj, add=T)

ndvi <- crop(ndvi.raw, guan.prj) #subset to square bounding study area
ndvi <- mask(ndvi, guan.prj) #mask out areas outside border

ndvi.array <- as.array(ndvi) #convert raster into y-by-x-by-time array

grid.inds <- expand.grid(rev(1:(dim(ndvi.array)[1])), 1:dim(ndvi.array)[2]) 

ndvi.matrix <- NULL #turn into locations-by-time matrix
for(ii in 1:dim(ndvi.array)[3]){
  ndvi.matrix <- cbind(ndvi.matrix, c(ndvi.array[,,ii]))
}
#ndvi.matrix[ndvi.matrix==0] <- NA
ndvi.matrix[ndvi.matrix<=2000] <- NA

keep<- apply(ndvi.matrix, 1, function(x){!all(is.na(x))}) #find pixels outside guanacaste border for removal

ndvi.matrix <- ndvi.matrix[keep,]
grid.coords <- grid.inds[keep,]
ndvi.matrix <- ndvi.matrix[,image.years < 2022]


## Synchrony over all years -----------------------------------------------------------------------

cormat <- cor(t(ndvi.matrix), use="pairwise.complete.obs")
diag(cormat) <- 0
ndvi.clust <- cluseigen(cormat)
ndvi.clust <- ndvi.clust[[2]] #perform clustering to find groups
ndvi.modularity <- modularity(cormat, membership=ndvi.clust, decomp=TRUE) #modularity provides uncertainty

nodestr.scale <- ndvi.modularity$nodeQ
for(ii in 1:max(ndvi.clust)){
  nodestr.scale[ndvi.clust==ii] <- nodestr.scale[ndvi.clust==ii]/max(nodestr.scale[ndvi.clust==ii])
}

groundtruth$pch <- NA
groundtruth$pch[groundtruth$forestType=="dry"] <- 1
groundtruth$pch[groundtruth$forestType=="transition"] <- 2
groundtruth$pch[groundtruth$forestType=="wet"] <- 3
groundtruth$pch[groundtruth$forestType=="cloud"] <- 4


plot(grid.coords[,2], grid.coords[,1], pch=16, col=ndvi.clust, cex=nodestr.scale)
#points(grid.coords[nodestr.scale < 0.2,2], grid.coords[nodestr.scale < 0.2,1], pch=16, col="blue", cex=0.5)
text(grid.coords[ndvi.clust==2 & nodestr.scale==1,2], grid.coords[ndvi.clust==2 & nodestr.scale==1,1], "1", col="grey")
text(grid.coords[ndvi.clust==1 & nodestr.scale==1,2], grid.coords[ndvi.clust==1 & nodestr.scale==1,1], "2", col="grey")
text(grid.coords[which.min(nodestr.scale),2], grid.coords[which.min(nodestr.scale),1], "3", col="grey")

## rasterize the outputs
group.matrix <- matrix(NA, nrow=dim(ndvi.array)[1], ncol=dim(ndvi.array)[2])
certainty.matrix <- matrix(NA, dim(ndvi.array)[1], dim(ndvi.array)[2])

for(ii in 1:nrow(grid.coords)){
  x1 <- grid.coords[ii,1]
  y1 <- grid.coords[ii,2]
  group.matrix[x1,y1] <- ndvi.clust[ii]
  certainty.matrix[x1,y1] <- nodestr.scale[ii]
}

group.raster <- raster(group.matrix[rev(1:dim(group.matrix)[1]),], template=ndvi)
certainty.raster <- raster(certainty.matrix[rev(1:dim(group.matrix)[1]),], template=ndvi)

# writeRaster(group.raster, filename=paste0("../Outputs/Intermediate/groupRaster","allYrs",".tif"), format="GTiff", overwrite=TRUE)
# writeRaster(certainty.raster, filename=paste0("../Outputs/Intermediate/certaintyRaster","allYrs",".tif"), format="GTiff", overwrite=TRUE)


# Look at time series of particular points

# dry <- ndvi.matrix[ndvi.clust==2 & nodestr.scale==1,]
# wet <- ndvi.matrix[ndvi.clust==1 & nodestr.scale==1,]
# ecotone <- ndvi.matrix[which.min(nodestr.scale),]

dry <- colMeans(ndvi.matrix[ndvi.clust==2,], na.rm=T)
wet <- colMeans(ndvi.matrix[ndvi.clust==1,], na.rm=T)
ecotone <- colMeans(ndvi.matrix[nodestr.scale < 0.2,], na.rm=T)

cols = brewer.pal(3, "Set2")

png("../Outputs/fig1_delineation.png", width=6.5, height=4, units="in", res=300)

layout(matrix(c(1,1,1,2:4), ncol=2, nrow=3), widths=c(0.6,0.4))
par(mar=c(2.1,1.1,2.1,0.6))

plot(grid.coords[,2], grid.coords[,1], pch=16, col=cols[ndvi.clust], cex=nodestr.scale, asp=1,
     xaxt="n", yaxt="n")
points(grid.coords[nodestr.scale < 0.2,2], grid.coords[nodestr.scale < 0.2,1], pch=1, col=cols[3], cex=0.5, lwd=0.5)
pp <- par("usr")
text(pp[1]+0.05*abs(diff(pp[1:2])), pp[4]-0.05*abs(diff(pp[3:4])), "a)")
#text(grid.coords[ndvi.clust==2 & nodestr.scale==1,2], grid.coords[ndvi.clust==2 & nodestr.scale==1,1], "1", col="black")
#text(grid.coords[ndvi.clust==1 & nodestr.scale==1,2], grid.coords[ndvi.clust==1 & nodestr.scale==1,1], "2", col="black")
#text(grid.coords[which.min(nodestr.scale),2], grid.coords[which.min(nodestr.scale),1], "3", col="black")
legend("top", pch=c(19,19,1), col=cols[c(2,1,3)], legend=c("Dry forest", "Rainforest", "Ecotone"), 
       ncol=3, bty="n", inset=-0.08, xpd=TRUE)

par(mar=c(2.1,3.1,2.1,1.1), mgp=c(2,0.8,0))

plot(image.dates[image.years < 2022], dry/10000, type="l", ylab="NDVI", ylim=c(0.2,1))
abline(v=as.Date(paste(2000:2022,"01","01", sep="-")), col="lightgrey", lty=3)
mtext("Dry forest", cex=2/3, line=0.1)
pp <- par("usr")
text(pp[1]+0.05*abs(diff(pp[1:2])), pp[4]-0.1*abs(diff(pp[3:4])), "b)")
plot(image.dates[image.years < 2022], ecotone/10000, type="l", ylab="NDVI", ylim=c(0.2,1))
abline(v=as.Date(paste(2000:2022,"01","01", sep="-")), col="lightgrey", lty=3)
mtext("Ecotone", cex=2/3, line=0.1)
text(pp[1]+0.05*abs(diff(pp[1:2])), pp[4]-0.1*abs(diff(pp[3:4])), "c)")
plot(image.dates[image.years < 2022], wet/10000, type="l", ylab="NDVI", ylim=c(0.2,1))
abline(v=as.Date(paste(2000:2022,"01","01", sep="-")), col="lightgrey", lty=3)
mtext("Rainforest", cex=2/3, line=0.1)
text(pp[1]+0.05*abs(diff(pp[1:2])), pp[4]-0.1*abs(diff(pp[3:4])), "d)")

dev.off()


## Loop over individual years ---------------------------------------------------------------------

pdf("../Outputs/RoughFigs/clusterGroups_allYrs_force2clusters_removeLowNDVI.pdf", onefile=TRUE)

## loop over years and do analysis
for(yy in years){
  
  if(yy ==  2022){break}

  ndvi.yy <- ndvi.matrix[ , image.years==yy]
  dmin <- ceiling(ncol(ndvi.yy))*(2/3)

  keep.yy <- apply(ndvi.yy, 1, function(x){sum(!is.na(x)) > dmin})

  if(sum(keep.yy) == 0){next}

  ndvi.yy <- ndvi.yy[keep.yy, ]
  grid.coords.yy <- grid.coords[keep.yy, ]


  cormat <- cor(t(ndvi.yy), use="pairwise.complete.obs") #create correlation matrix
  diag(cormat) <- 0

  ndvi.clust <- cluseigen(cormat)
  ndvi.clust <- ndvi.clust[[2]] #perform clustering to find groups
  #print(paste("n clusters = ", length(ndvi.clust), sep=""))
  ndvi.modularity <- modularity(cormat, membership=ndvi.clust, decomp=TRUE) #modularity provides uncertainty

  nodestr.scale <- ndvi.modularity$nodeQ
  for(ii in 1:max(ndvi.clust)){
    nodestr.scale[ndvi.clust==ii] <- nodestr.scale[ndvi.clust==ii]/max(nodestr.scale[ndvi.clust==ii])
  }

  plot(grid.coords.yy[,2], grid.coords.yy[,1], pch=16, col=ndvi.clust, cex=nodestr.scale, main=yy)
  #points(grid.inds.sub[nodestr.scale < 0.15,2], grid.inds.sub[nodestr.scale < 0.15,1], pch=16, col="blue", cex=0.5)

  ## rasterize the outputs
  group.matrix <- matrix(NA, nrow=dim(ndvi.array)[1], ncol=dim(ndvi.array)[2])
  certainty.matrix <- matrix(NA, dim(ndvi.array)[1], dim(ndvi.array)[2])
  
  for(ii in 1:nrow(grid.coords.yy)){
    x1 <- grid.coords.yy[ii,1]
    y1 <- grid.coords.yy[ii,2]
    group.matrix[x1,y1] <- ndvi.clust[ii]
    certainty.matrix[x1,y1] <- nodestr.scale[ii]
  }
  
  # ecotone.thresh <- 0.1 #threshold level of group assignment uncertainty; below this is considered ecotone
  # ecotone.matrix <- group.matrix
  # ecotone.matrix[certainty.matrix < ecotone.thresh] <- 3
  
  group.raster <- raster(group.matrix[rev(1:dim(group.matrix)[1]),], template=ndvi)
  certainty.raster <- raster(certainty.matrix[rev(1:dim(group.matrix)[1]),], template=ndvi)
  #ecotone.raster <- raster(ecotone.matrix[rev(1:dim(group.matrix)[1]),], template=ndvi)
  #plot(ecotone.raster)
  
  writeRaster(group.raster, filename=paste0("../Outputs/Intermediate/groupRaster",yy,".tif"), format="GTiff", overwrite=TRUE)
  writeRaster(certainty.raster, filename=paste0("../Outputs/Intermediate/certaintyRaster",yy,".tif"), format="GTiff", overwrite=TRUE)
  
}

dev.off()


## make DEM map

dem <- raster("../Data/10n090w_20101117_gmted_med075.tif")
dem.guan <- crop(dem, guan.prj)


png("../Outputs/fig_dem.png", units="in", width=6.5, height=6.5, res=300)

plot(dem.guan)
plot(guan.prj, add=TRUE)  

dev.off()


## 