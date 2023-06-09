library(geosphere)
library(ggplot2)
library(raster)
library(rgdal)

setwd("/Users/Jeremiah/Documents")
psi.values <- read.csv("psi_df_20.csv",header=T,stringsAsFactors = F)

#for western analysis
setwd("/Users/Jeremiah/Dropbox/Campanula mating system/Papers/Range Expansion Routes/paper/MolEcol Revision/Analyses/TDoA/")
psi.values <- read.csv("psi_df_20_west.csv",header=T,stringsAsFactors = F)
setwd("/Users/Jeremiah/Documents")
coords <- read.csv("Rad_Pop_LATLON_CP.csv",header=T,stringsAsFactors = F)
pops <- unique(c(psi.values$pop1,psi.values$pop2))

#for western analysis, use next 3 lines
#removed "MI127"
west <- c("AR125", "AR56", "IA10", "IN77", "KS60", "MN117", "MN118", "MI126", "MI127", "MO115", "MO116", "MO49", "MO57", "OK61", "WI128") 
psi.values <- psi.values[psi.values$pop1 %in% west,]
psi.values <- psi.values[psi.values$pop2 %in% west,]

#always use this
coords <- coords[coords$pop %in%pops,]
sites <- coords

#subsetting for west
psi.values <- psi.values[psi.values$pop1 %in% west,]
psi.values <- psi.values[psi.values$pop2 %in% west,]
pops <- pops[pops %in% west]
coords <- coords[coords$pop %in% west,]
sites <- coords

output <- c()
for(i in pops) {
  for(j in pops) {
  output = c(output,sum(0,psi.values[psi.values$pop1==i & psi.values$pop2==j,]$psi))
  }
}

psi <- matrix(unlist(output), ncol = length(pops), byrow = T)
psi[lower.tri(psi)] <- t(psi)[lower.tri(psi)]

model.1d <- function(xy, data, pct=0.01, f.dist="haversine"){
  if (f.dist=="euclidean"){
    f.dist <- function(i, j){
      sqrt((ix -jx)^2 + (iy-jy)^2 )
    }
  }else{ if(f.dist=="haversine"){
    f.dist <- distHaversine
  }}
  
  y = xy[2] 
  x = xy[1]
  ixy = data[,1:2]
  jxy = data[,3:4]
  psi = data[,5]
  
  
  d = f.dist(ixy, c(x,y)) - f.dist(jxy, c(x,y))
  l = lm( psi ~ d ) 
  l$f.e = .5 * pct/ l$coefficients[2]
  l$rsq = summary(l)$r.squared
  
  return (l)
}

prep.tdoa.data <- function(coords, psi){
  locs <- coords[,c('long', 'lat')]
  n.locs <- nrow(locs)
  
  tdoa.data <- c()
  for(i in 1:n.locs){
    for(j in (i+1):n.locs){
      if( i>=j | j >n.locs) break
      tdoa.data <- rbind( tdoa.data, c(locs[i,], locs[j,], psi[i,j]))
    }
  }
  
  tdoa.data <- matrix(unlist(tdoa.data), ncol=5)
  
  return(tdoa.data)
}

summary.origin.results <- function(res){
  bbox <- res$bbox
  orig <- which(t(res[[2]])==max((res[[1]]>0)*res[[2]],
                                 na.rm=T), arr.ind=T)        
  s1<-seq(bbox[1,1],bbox[1,2],length.out=res$xlen) 
  s2<-seq(bbox[2,1],bbox[2,2],length.out=res$ylen) 
  
  coords_orig <- c(s1[orig[1,2]], s2[orig[1,1]])
  max.reg <- res[[3]][[orig[1,2]]][[orig[1,1]]]
  
  smry <- summary(max.reg)
  q <- smry$coefficients[2,1]
  r1 = 1/(1+2*1000*q)
  r2 = 1/(1+2*10000*q)
  r3 = 1/(1+2*100000*q)
  r <- 0.99
  d1pc <- (1-r)/(2*q*r)
  
  res.data <- ( c( coords_orig, q*1000, r1, r2, r3, d1pc/1000, smry$adj.r.squared,
                   smry$coefficients[2,4] * 10^4) )
  res.data <- as.data.frame(t(res.data))
  names(res.data) <- c("long", "lat", "q", 
                       "r1","r10","r100",
                       "d1","rsq","pval")
  return(res.data)
}

summary.origin.result.list <- function(res){
  a <- sapply(res$tbl, summary)
  a <- data.frame(a)
  na <- sapply(res$regions, paste, collapse="+")
  na[na==""] <- "ALL"
  names(a) <- na
  return(a)
}

coords2country = function(points){  
  countriesSP <- getMap(resolution='high')
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  indices$ADMIN  
}

single.origin <- function(tdoa.data, bbox,  pop.coords,
                          pct=0.01,
                          xlen=100, ylen=100, 
                          exclude.ocean=T,
                          exclude.land=F,
                          ...){
  
  
  #define locs for estimate
  s1<-seq(bbox[1,1],bbox[1,2], length.out=xlen)
  s2<-seq(bbox[2,1],bbox[2,2], length.out=ylen)
  coords <- expand.grid(s1,s2)
  ij <- expand.grid(1:length(s1), 1:length(s2))
  
  if(exclude.ocean){
    cc <- coords2country(coords)
    to.keep <- !is.na(cc)
  } else {
    to.keep <- rep(T, nrow(coords))
  }
  if(exclude.land){
    cc <- coords2country(coords)
    to.keep <- is.na(cc)
  }
  
  # init output
  d0 <- matrix(NA, ncol=ylen, nrow=xlen)
  rsq <- matrix(NA, ncol=ylen, nrow=xlen)
  mdlq <- list()
  for(i in 1:xlen){
    mdlq[[i]] <- list()
  }
  
  for(r in 1:nrow(coords)){
    i <- ij[r,1]
    j <- ij[r,2]
    x <- coords[r,1]
    y <- coords[r,2]
    
    if(to.keep[r]){
      mdl <- model.1d( xy=c(x,y), data=tdoa.data, pct=pct)
      d0[i, j] <- mdl$f.e
      rsq[i, j] <- mdl$rsq
      mdlq[[i]][[j]] <- mdl
    }
  }
  res <- list( d0=d0, rsq=rsq, mdlq=mdlq, bbox=bbox, xlen=xlen, 
               ylen=ylen, coords=pop.coords)
  class(res) <- 'origin.results'
  return(res)
}

#analysis +-5, plotting +-10
get.sample.bbox <- function(samples) {
  
  s <- c('long', 'lat')
  mins <- apply(samples[,s],2, function(x) min(x)-10)
  maxs <- apply(samples[,s],2, function(x) max(x)+10)
  return(cbind(mins, maxs))
  
}

find.origin <- function(coords, psi, 
                        xlen=50,  ylen=50, doPlot=F, doPdf=F, 
                        ...){
  tdoa.data <- prep.tdoa.data(coords, psi)
  bbox <- get.sample.bbox(coords)
  
  res <- single.origin(tdoa.data,
                        bbox = bbox,
                        coords,
                        xlen=xlen,
                        ylen=ylen, ...)
   
   
  return( res )
}

origin <- find.origin(coords,psi)

plot.origin.results <- function(x, n.levels=100, color.function=heat.colors,
                                color.negative='transparent',
                                add.map=T,
                                add.samples=T,
                                add.sample.het=F,
                                add.likely.origin=T,
                                asp=1,
                                ...){
  plot.default(NA, xlim=x$bbox[1,], ylim=x$bbox[2,], 
               xlab="", ylab="",
               xaxt="n", yaxt="n",
               xaxs='i', yaxs='i',
               asp=asp, ...)
  
  s1<-seq(x$bbox[1,1],x$bbox[1,2],length.out=x$xlen)
  s2<-seq(x$bbox[2,1],x$bbox[2,2],length.out=x$ylen)
  
  rel <- (x[[1]]>0) * (x[[2]]-min(x[[2]],na.rm=T)) /
    (max(x[[2]],na.rm=T)-min(x[[2]],na.rm=T))+0.001
  
  levels <- c(0,quantile(rel[rel>0.001], 0:n.levels/n.levels, na.rm=T) + 
                1e-6 * 0:n.levels/n.levels)
  cols <- c(color.negative, color.function(n.levels-1))
  
  
  .filled.contour(s1, s2, rel, levels, cols)
  
  
  # rect(x$bbox[1,1], x$bbox[2,1], x$bbox[1,2], x$bbox[2,2], border='black',
  #      lwd=2, col=NULL)
  
  if(add.likely.origin){
    points(summary(x)[,1:2], col='black', pch='x', cex=2)
  }
  if(add.map){
    require(rworldmap)
    m <- getMap("high")
    plot(m, add=T, lwd=1.3)
  }
  
  if(add.sample.het){
    samples <- x$coords
    hets <- (samples$hets - min(samples$hets) )/(
      max(samples$hets) - min(samples$hets))
    points( samples$long, samples$lat,
            pch=16, cex=3, col=grey(hets) )
    points( samples$long, samples$lat,
            pch=1, cex=3, col="black",lwd=2 )
  }
  else if(add.samples){
    samples <- x$coords
    points( samples$long, samples$lat,
            pch=16, cex=1, col="black",lwd=1 )
  }
}

require(rworldmap)
plot.origin.results(origin, add.sample.het=F)

#get likely origin
summary(origin)[,1:2]

#making nice countour plots
n.levels <-30
rel <- (origin[[1]]>0) * (origin[[2]]-min(origin[[2]],na.rm=T)) /
  (max(origin[[2]],na.rm=T)-min(origin[[2]],na.rm=T))+0.001

levels <- c(0,quantile(rel[rel>0.001], 0:n.levels/n.levels, na.rm=T) +
              1e-6 * 0:n.levels/n.levels)

s1<-seq(origin$bbox[1,1],origin$bbox[1,2],length.out=origin$xlen)
s2<-seq(origin$bbox[2,1],origin$bbox[2,2],length.out=origin$ylen)

plt.data <- cbind(expand.grid(long=s1,lat=s2),z=as.vector(rel))
plt.data[is.na(plt.data$z),] <- 0 
plt.data[plt.data$z<0.4,]$z <- NA

#From Nate May 31st
library(devtools)
install_github("clauswilke/ggisoband")
library(ggplot2)
library(ggisoband)
plt.data <- cbind(expand.grid(long=s1,lat=s2),z=as.vector(rel))
plt.data[is.na(plt.data$z),] <- 0 
plt.data[plt.data$z<.02,]$z <- NA

ggplot(plt.data, aes(x=long, y=lat, z=z, fill=..level.., alpha=0.2)) +
  scale_fill_gradientn(colors=c("transparent","yellow","red"), values=c(0,0.02,1), na.value="transparent") +
  geom_isobands(color="orange")+theme(legend.position="none") 

setwd("/Users/Jeremiah/Dropbox/Campanula mating system/Papers/Range Expansion Routes/analyses/pca")
shapeData <- shapefile("ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp")
shapeData@data$id <- rownames(shapeData@data)

watershedPoints <- fortify(shapeData, region = "id")
watershedDF <- merge(watershedPoints, shapeData@data, by = "id")
watershedDF.1 <- subset(watershedDF,name %in% c("Mississippi") & long >-99.65)

setwd("/Users/Jeremiah/Dropbox/Campanula mating system/Papers/Range Expansion Routes/analyses/map")
extent <- read.csv("lgm_extent.csv", header=T)
sites <- read.csv("Rad_Pop_LATLON_CP.csv", stringsAsFactors=F)
sites$long <- sites$lon
pops <- c("AL2012", "AL79", "ALBG", "AR125", "AR56", "IA10", "IN46", "IN77", "KS60", "KY51", "MI126", "MI127", "MN117", "MN118", "MO115", "MO116", "MO49", "MO57", "OH119", "OH64", "OK61", "PA27", "TN19", "WI128")
sites <- sites[sites$pop %in% west,]
states <- map_data("state")
apple <- read.csv("apple_poly.csv")
ggplot(data = subset(states,long > -100 & long < -75 & lat > 25 & lat < 50)) + theme_bw() + ylim(c(25,50)) + xlim(c(-100,-75))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_polygon(data=apple, aes(long,lat), fill="dimgrey", alpha=1) +
  geom_raster(data=plt.data, aes(long,lat,fill=z,alpha=0.6),interpolate = T) +
  scale_fill_gradientn(colors = c("transparent","green", "red"), values=c(0.05,0.5,1), na.value="transparent") +
  geom_contour(data=plt.data,aes(long,lat,z=z,fill=z),color="white") +
  #geom_polygon(data=extent, aes(long,lat), fill="deepskyblue2", alpha=0.2) +
  geom_polygon(aes(x = long, y = lat, group = group), fill=NA,color = "black",size=0.25) + 
  geom_path(data=watershedDF.1, aes(x = long, y = lat, group = group), color = 'blue',size=0.75) +
  geom_point(aes(x=summary(origin)[,1],y=summary(origin)[,2]), shape='x', size=9,color="black") +
  geom_point(data=sites, aes(long,lat), size=3, color="grey27") +
  coord_fixed(1.3) +
  xlab("longitude") +
  ylab("latitude") +
  theme(legend.position="none") +
  guides(fill=FALSE)  # do this to leave off the color legend

