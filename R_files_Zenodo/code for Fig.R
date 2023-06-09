library(gtools)
library(ncdf4)
library(MuMIn)
library(zoo)
library(scales) 
library(nlme)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(pracma)
library(FactoMineR)
library(lmtest)
library(MuMIn)
library(broom)
library(reshape2)

# reload ERSST data
nc <- nc_open("/Users/MikeLitzow 1/Documents/R/NSF-GOA/sst.mnmean.v4.nc") 

# assign dates (Days since January 1, 1800)
d <- dates(ncvar_get(nc, "time"), origin=c(1,15,1800))

# Extract SST data for desired period and locations:
# Pick start and end dates (January 1950-December 2012):
d1 <- d[1153:1908]

# Select latitude and longitude: 20-66 deg. N, 132-250 deg. E:
x <- ncvar_get(nc, "lon", start=67, count=60)
y <- ncvar_get(nc, "lat", start=12, count=24)

SST1 <- ncvar_get(nc, "sst", start=c(67,12,1153), count=c(60,24,length(d1)))

# process
SST1 <- aperm(SST1, 3:1)  # First, reverse order of dimensions - "transpose" array
SST1 <- SST1[,24:1,]  # Reverse order of latitudes to be increasing for convenience in later plotting
y <- rev(y)  # Also reverse corresponding vector of lattidues
SST1 <- matrix(SST1, nrow=dim(SST1)[1], ncol=prod(dim(SST1)[2:3]))  # Change to matrix with column for each grid point, rows for monthly means

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes
dimnames(SST1)  <- list(as.character(d1), paste("N", lat, "E", lon, sep=""))

# identify columns in SST matrix corresponding to land
land <- is.na(colMeans(SST1)) 

# For analysis, we only use the columns of the matrix with non-missing values:
X1 <- SST1[,!land] 

# To remove seasonal trend we compute long-term means for each month and substract them:
m1 <- months(d1)  # Extracts months from the date vector
y1 <- years(d1)
f <- function(x) tapply(x, m1, mean)  # function to compute monthly means for a single time series
mu1 <- apply(X1, 2, f)	# compute monthly means for each time series (cell)
mu1 <- mu1[rep(1:12, length(d1)/12),]  # replicate means matrix for each year at each location
X1.anom <- X1 - mu1   # compute matrix of anomalies

# now detrend
X1.anom.detr <- detrend(X1.anom)

# get a vector of weights (square root of the cosine of latitude)
lat.weights <- lat[!land]
weight <- sqrt(cos(lat.weights*pi/180))

# EOF by era
EOF.early <- svd.triplet(cov(X1.anom.detr[y1 <= 1988,]), col.w=weight) #weighting the columns
EOF.late <- svd.triplet(cov(X1.anom.detr[y1 > 1988,]), col.w=weight)

# get loadings for PC1/2 by era
eig.1.early <- EOF.early$U[,1]
eig.2.early <- EOF.early$U[,2]

eig.1.late <- EOF.late$U[,1]
eig.2.late <- EOF.late$U[,2]

# get % variance explained by era
var.early <- 100*round(prop.table(EOF.early$vs),3)
            
var.late <- 100*round(prop.table(EOF.late$vs),3)

# set colors
new.col <- my.col <- tim.colors(64)
grays <- c("gray98", "gray97", "gray96", "gray95", "gray94", "gray93", "gray92", "gray91", "gray90", "gray89", "gray88")

my.col[22:43] <- c(grays[11:1], grays)
new.col[27:36] <- c(grays[5:1], grays[1:5])

# and plot

png("Fig S2.png", 1.2*11.4/2.54, 11.4/2.54, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,2), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# set the limit for plotting 
lim <- range(eig.1.early, eig.1.late, eig.2.early, eig.2.late, na.rm=T)

z <- rep(NA, ncol(SST1))
z[!land] <- eig.1.early
z <- t(matrix(z, length(y))) 

image.plot(x,y,z, col=my.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("a", adj=0.05, line=-1.3, cex=mt.cex)
mtext(paste("PC1 1950-1988 (", var.early[1], "%)", sep=""), cex=0.8)

z <- rep(NA, ncol(SST1))
z[!land] <- -eig.1.late # reversing the loadings to match 1950-1988
z <- t(matrix(z, length(y)))

image.plot(x,y,z, col=my.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("b", adj=0.05, line=-1.3, cex=mt.cex)
mtext(paste("PC1 1989-2012 (", var.late[1], "%)", sep=""), cex=0.8)

z <- rep(NA, ncol(SST1))
z[!land] <- eig.2.early
z <- t(matrix(z, length(y))) 

image.plot(x,y,z, col=my.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("c", adj=0.05, line=-1.3, cex=mt.cex)
mtext(paste("PC2 1950-1988 (", var.early[2], "%)", sep=""), cex=0.8)

z <- rep(NA, ncol(SST1))
z[!land] <- eig.2.late
z <- t(matrix(z, length(y)))

image.plot(x,y,z, col=my.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("d", adj=0.05, line=-1.3, cex=mt.cex)
mtext(paste("PC2 1989-2012 (", var.late[2], "%)", sep=""), cex=0.8)

dev.off()
