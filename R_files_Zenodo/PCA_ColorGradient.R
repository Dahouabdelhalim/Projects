#PCA with color gradient
# PCA function created by C.A. Buerkle, rest of script created by M.E.F. LaCava

####Create colors grid ####
## Test color mixing
#install.packages("colorspace")
library(colorspace)
mixcolor(0.5, RGB(1,0,0), RGB(0,0,1))
barplot(1:3,col=c(rgb(1, 0, 0),rgb(0.5,0,0.5), rgb(0,0,1)))

#blues <- c(rgb(0,0,1),rgb(.24,.24,1),rgb(.48,.48,1),rgb(.73,.73,1))
blues <- c(rgb(0,0,1),rgb(.3,.3,1),rgb(.61,.61,1),rgb(.91,.91,1))
barplot(1:4,col=blues)
#reds <- c(rgb(1,0,0),rgb(1,.18,.18),rgb(1,.37,.37),rgb(1,.55,.55))
reds <- c(rgb(1,0,0),rgb(1,.3,.3),rgb(1,.61,.61),rgb(1,.91,.91))
barplot(1:4,col=reds)

## Make matrix of 16 colors for grid
# used website to pick initial 4 colors in blue and 4 in red
#   https://meyerweb.com/eric/tools/color-blend/#:::rgbd
col.grid <- matrix(nrow=4,ncol=4)
# two gradients of color
#blues <- list(c(.73,.73,1),c(.48,.48,1),c(.24,.24,1),c(0,0,1))
blues <- list(c(.91,.91,1),c(.61,.61,1),c(.3,.3,1),c(0,0,1))
#reds <- list(c(1,.55,.55),c(1,.37,.37),c(1,.18,.18),c(1,0,0))
reds <- list(c(1,.91,.91),c(1,.61,.61),c(1,.3,.3),c(1,0,0))
# mix gradients to create grid
for (i in 1:nrow(col.grid)){
  for (j in 1:ncol(col.grid)){
    r <- (blues[[i]][1]+reds[[j]][1])/2
    g <- (blues[[i]][2]+reds[[j]][2])/2
    b <- (blues[[i]][3]+reds[[j]][3])/2
    col.grid[i,j] <- rgb(r,g,b)
  }
}
col.grid
barplot(1:16,col=col.grid)
col.vect <- as.vector(t(col.grid))

## Import state polygon and sample locations
state <- readOGR(dsn=paste(getwd(),"/LandscapeFeatures",sep=""), layer="state")
proj4string(state) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
latlong <- read.csv("LatLong_398ind.csv", header=F, stringsAsFactors=F)
#Or microsat samples
#WD <- "/Users/melanielacava/Data/PHmsats/"
#setwd(WD)
#latlong <- read.csv("LatLong_274ind.csv", header=F, stringsAsFactors=F)

names(latlong) <-c("id","x","y")
coordinates(latlong) <- ~ x + y
proj4string(latlong) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj4string(latlong)==proj4string(state)
plot(state)
plot(latlong,add=T)
#add grid to state
bb <- bbox(state)
cs <- c((bb[1,2]-bb[1,1])/4,(bb[2,2]-bb[2,1])/4) #cell size
cc <- bb[,1]+(cs/2)
cd <- ceiling(diff(t(bb))/cs) #number of cells per direction
grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
grd
sp_grd <- SpatialGridDataFrame(grd,
                               data=data.frame(id=1:prod(cd)),
                               proj4string=CRS(proj4string(state)))
#import interstate to overlay on grid
highway <- readOGR(dsn=paste(getwd(),"/LandscapeFeatures",sep=""), layer="InterstateOnly")
#check it worked
library("lattice")
spplot(sp_grd, colorkey=F, region=T, col.regions=as.vector(t(col.grid)),
       panel = function(...) {
         panel.gridplot(..., border="black")
         #         sp.polygons(state)
         sp.polygons(highway,lwd=2) #add highways to help orient map
         sp.points(latlong, cex=1,col="black",pch=16)
         #         panel.text(...)
       })
#export as pdf
pdf(file="samples_gradientmap.pdf",width=5,height=5,useDingbats=FALSE)
spplot(sp_grd, colorkey=F, region=T, col.regions=as.vector(t(col.grid)),
       panel = function(...) {
         panel.gridplot(..., border="black")
         #         sp.polygons(state)
         sp.polygons(highway,lwd=2) #add highways to help orient map
         sp.points(latlong, cex=0.8,col="black",pch=16)
         #         panel.text(...)
       })
dev.off()

## Assign colors to points
over(latlong,sp_grd) #which grid cell is each sample in
head(latlong)
latlong$cell <- over(latlong,sp_grd)[,1] 
latlong$color <- NA
for (i in 1:nrow(latlong)){
  latlong$color[i] <- col.vect[latlong$cell[i]]
}
head(latlong)
#check it looks right on map
plot(state)
plot(latlong,add=T,col=latlong$color)


#### SNP PCA ####
#PCA function
do.pca<-function(gmat, write.gcov=FALSE, inds=""){
  gmn<-apply(gmat,1,mean, na.rm=T) #takes a mean across individuals for each locus
  gmnmat<-matrix(gmn,nrow=nrow(gmat),ncol=ncol(gmat)) #creates matrix with mean filled in for entire row
  gprime<-gmat-gmnmat ## remove mean 
  
  gcovarmat<-matrix(NA,nrow=ncol(gmat),ncol=ncol(gmat))
  for(i in 1:ncol(gmat)){
    for(j in i:ncol(gmat)){
      if (i==j){
          gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs") #only use loci that both samples in pair have a genotype for
      }
      else{
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs") #only use loci that both samples in pair have a genotype for
        gcovarmat[j,i]<-gcovarmat[i,j]
      }
    }
  }
  if(write.gcov==TRUE){
    inds<-ifelse(inds == "", paste("i", 1:ncol(gmat)), inds)
    write.table(round(gcovarmat,5),file="gcovarmat.txt",
                quote=F,row.names=F,col.names=inds)
  }
  prcomp(x=gcovarmat,center=TRUE,scale=FALSE)
}

#import data
geno <- read.table("pntest_mean_398_sorted.txt", header=F, stringsAsFactors=F)
#convert genotypes to matrix
pcatest <- as.matrix(geno)
pca.snp <- do.pca(pcatest)
pca.summary.snp<-summary(pca.snp)
#plot with color grid
plot(pca.snp$x[,'PC2']~pca.snp$x[,'PC1'], cex=0.9, col=latlong$color, pch = 16, ann=F, mgp=c(3,0.5,0),cex.axis=0.8)
mtext(side=1,text=paste("PC1 (", round(pca.summary.snp$importance[2,1]*100, 1), "%)", sep=""),line=1.6)
mtext(side=2,text=paste("PC2 (", round(pca.summary.snp$importance[2,3]*100, 1), "%)", sep=""),line=1.6)
#export as pdf
pdf(file="SNP_gradientPCA2.pdf",width=5,height=5,useDingbats=FALSE)
plot(pca.snp$x[,'PC2']~pca.snp$x[,'PC1'], cex=0.9, col=latlong$color, pch = 16, ann=F, mgp=c(3,0.5,0),cex.axis=0.8)
mtext(side=1,text=paste("PC1 (", round(pca.summary.snp$importance[2,1]*100, 1), "%)", sep=""),line=1.6)
mtext(side=2,text=paste("PC2 (", round(pca.summary.snp$importance[2,3]*100, 1), "%)", sep=""),line=1.6)
dev.off()

#### Msat PCA ####
#import data
msat.colors<- read.csv("ColorGrid_274ind.csv",header=T,stringsAsFactors=F)
library(adegenet)
ph.msat <- read.structure("PHmsats_11loci_HuntRegions_sorted_GenAlEx.stru", n.ind=274, n.loc=11, onerowperind=T, col.lab=1, col.pop=2, row.marknames=1)
# create scaling agent to replace missing data with mean allele frequency (ok b/c almost no missing data)
scl <- scaleGen(ph.msat, NA.method="mean")
pca.msat <- dudi.pca(scl,cent=F,scale=F,scannf=F,nf=20)
pca.msat$eig #the eigenvalues of the analysis, indicating the amount of variance represented by each principal component (PC)
#plot
plot(pca.msat$li$Axis1,pca.msat$li$Axis2, col=msat.colors$color, cex=0.9, pch=16, ann=F, mgp=c(3,0.5,0),cex.axis=0.8)
mtext(side=1,text=paste("PC1 (", round(pca.msat$eig[1], 1), "%)", sep=""),line=1.6)
mtext(side=2,text=paste("PC2 (", round(pca.msat$eig[2], 1), "%)", sep=""),line=1.6)
#export as pdf
pdf(file="msat_gradientPCA.pdf",width=5,height=5,useDingbats=FALSE)
plot(pca.msat$li$Axis1,pca.msat$li$Axis2, col=msat.colors$color, cex=0.9, pch=16, ann=F, mgp=c(3,0.5,0),cex.axis=0.8)
mtext(side=1,text=paste("PC1 (", round(pca.msat$eig[1], 1), "%)", sep=""),line=1.6)
mtext(side=2,text=paste("PC2 (", round(pca.msat$eig[2], 1), "%)", sep=""),line=1.6)
dev.off()
