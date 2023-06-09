library(raster)
library(abc)

pdf("Figure5.pdf", paper="special", family = "Helvetica", onefile = FALSE, width=6.83, height=3.27, colormodel="rgb",  pointsize=18)

op <- par(mar = c(2, 2, 0.5, 2), mgp=c(1, 0.25, 0), las=1, ps=12)
L = c(0.03,0.36,0.1,0.95) 	#Top L
M = c(0.45,0.675,0.3,0.85)	#Bottom L
R = c(0.765,0.99,0.3,0.85)	#Top R

plot.new()
par(new = "TRUE",              
    plt = L,
    las = 1,                      
    ps=12, pty="m",tck = -0.03) 	



#### Fig A
nrows = 73
ncols = 138
df1 <- read.table("Rcustom/Map_100x100_Nean_veg2.asc", skip = 6, header = FALSE, sep = "\\t", nrows=nrows)
r.mat <- as.matrix(df1)
rmat2<-ifelse(r.mat==-9999,0, r.mat)
rmat2<-ifelse(rmat2>9 & rmat2 < 13,3, rmat2)

#0 sea
#1 H. sapiens
#2 Neanderthal and H. sapiens
#3 Himalaya

breakpoints <- c(-0.5,0.5,1.5, 2.5, 3.5)
colors <- c(grey(0.9),grey(0.4), "darkred", grey(0.7))

r <- raster(rmat2, xmn=0, xmx=ncols-1, ymn=0, ymx=nrows-1)
plot(r, breaks=breakpoints, col=colors, legend=F, las=1, axes=F, box=F)

points<-read.table("Rcustom/GeneSamples.sam", skip = 1, header = FALSE, sep = "\\t")[,-c(3,6)]
colnames(points)<-c("Country", "Samples", "Coor1", "Coor2")

points(points[1,"Coor2"]/100, points[1,"Coor1"]/100, pch=16, cex= 1, col="goldenrod2" )
points(points[2,"Coor2"]/100, points[2,"Coor1"]/100, pch=16, cex= 1, col="goldenrod2" )
points(2700/100, 2300/100, pch=16, cex= 1, col="darkslategray3" )
#####	

#### Fig B
par(new = "TRUE",              
    plt = M,   
    las = 1,                      
    ps=12, pty="s", tck = -0.03) 	
    
gamma<-read.csv("Results/neanderthal.csv", header=T)[,1]
intro<-1-read.csv("Results/neanderthal.csv", header=T)[,-1]
dat<-intro[, c("France", "China")]

obs<-dat[1,]
obs[1,]<-c(0.02, 0.02)

mpar2<-abc(obs, gamma, dat,tol=0.02,method="neuralnet")
summary(mpar2)

dens<-density(mpar2$adj.values[,1], weights= mpar2$weights/sum(mpar2$weights), adjust=3)

plot(dens, xlim=c(0,0.04), ylim=c(0,150), type="n", lty=1, lwd=3, bty="l", xlab=expression(paste("Interbreeding success rate ", (italic(gamma)))), ylab="Density", cex.lab=0.6, cex.axis=0.6, xaxs="i",yaxs="i", las=1, main="")

xx <- cbind(dens$x,dens$x)
yy <- cbind(dens$y,dens$y)
polygon(xx,yy, col=grey(0.9),border = "black", lty= 1, lwd=1)

dens2<-cbind(x=dens$x, y=dens$y)
dens2<-subset(dens2, dens2[,"x"]>=0.004)
dens2<-subset(dens2, dens2[,"x"]<=0.013 )

xx2 <- cbind(dens2[,"x"],rev(dens2[,"x"]))
yy2 <- cbind(dens2[,"y"],0)
polygon(xx2,yy2, col=grey(0.4),border = NA, lty= 1, lwd=1)

lines(c(0,0.04), dunif(c(0,0.04), min=0, max=0.04), lty=2)


#### Fig C
par(new = "TRUE",              
    plt = R,   
    las = 1,                      
    ps=12, pty="s", tck = -0.03) 	

superpose.eb <- function (x, y, ebl, ebu = ebl, length = 0.08, ...)
arrows(x, y + ebu, x, y - ebl, angle = 90, code = 1, lwd=0.8, length = length, ...)

dat2<-cbind(gamma=gamma, round(dat,2))
dat2<-subset(dat2, (dat2$France>=0.01 & dat2$France<=0.03) & (dat2$China>=0.01 & dat2$China<=0.03) )

dat2$enrich<-(dat2[,"China"]-dat2[,"France"])/dat2[,"France"]
datM<-aggregate(dat2, list(dat2[,1]), mean, na.rm=T)[,-1]

datSD<-aggregate(dat2, list(dat2[,1]), sd, na.rm=T)[,-1]
datCount<-aggregate(dat2, list(dat2[,1]), length)[,-1]


dat3<-cbind(datM[1,"enrich"], datM[5,"enrich"], datM[10,"enrich"])
colnames(dat3)<-c(0.004, 0.008, 0.013)

ic<-c(1.96*datSD[1,"enrich"]/sqrt(datCount[1,"enrich"]), 1.96*datSD[5,"enrich"]/sqrt(datCount[5,"enrich"]), 1.96*datSD[10,"enrich"]/sqrt(datCount[10,"enrich"]))

bar1<-barplot(100*dat3, space=c(0.1,0.25),width=1,axis.lty=1,beside=TRUE, col=grey(0.8), mgp=c(1, 0.25, 0.1), ylim=c(0,40), border=NA, las=1, ylab="Augmentation (%)", cex.lab=0.6, cex.axis=0.6, cex.names=0.6,  xlab=expression(paste("Interbreeding success rate ", (italic(gamma)))))

superpose.eb(bar1, 100*dat3, ebl=0, ebu=100*ic)
superpose.eb(bar1, 100*dat3, ebl=0, ebu=-100*ic)

lines(c(-6, 2, 3.75), c(12, 12, 12), lty=2)
lines(c(-6, 2, 3.75), c(20, 20, 20), lty=2)


par(xpd=NA)
text(-0.85, 40, "C", font=2)
text(-6.15, 40, "B", font=2)
text(-12.5, 40, "A", font=2)


par(op)             
dev.off()
