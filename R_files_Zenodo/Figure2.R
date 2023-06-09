library("shape")
source("Rcustom/colorlegend2.R") #modified from library 'shape'
source("Rcustom/filled.contour4.R") #modified from http://wiki.cbr.washington.edu/qerm/sites/qerm/images/1/16/filled.contour3.R
mycolors<-colorRampPalette(c("black", "cyan", "aquamarine4", "darkolivegreen","burlywood", "gold1","darkred"))
rotate <- function(x) t(apply(x, 2, rev))

pdf("Figure2.pdf", paper="special", family = "Helvetica", onefile = FALSE, width= 6.83, height= 6.83, colormodel="rgb",  pointsize=12)

op <- par(mar = c(0.5, 0, 0.5, 0), oma=c(4,4.8,2,5.2), mfcol=c(6,5), pty="s", las=1, ps=12)

#########################
#####    NC1    ######### 
#########################

#NC1 | Scen 1 | 3000gen | gamma = 0.05
inv<- read.csv("Results/Scen1_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen1_NCbis_mean.csv", header=T)

NC=1

propI <- sapply(0:24, function(i){ round(1-inv[22,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]


filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(2, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=2)

par(xpd=NA)
text(x=0.5, y=1.17, "NC1", cex=1, font=2)
text(x=0.5, y=1.05, expression(paste(italic('K'), " = 50, ", italic('Km'), " = 10", " | ", italic('K'), " = 50, ", italic('Km'), " = 10")), cex=0.5, font=2)
text(x=-0.3, y=0.5, "Invasive", cex=1, font=2, srt=90)
text(x=-0.48, y=0.5, "Whole area", cex=1, font=2, srt=90)
par(xpd=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(2, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=2)

par(xpd=NA)
text(x=-0.3, y=0.50, "Local", cex=1, font=2, srt=90)
text(x=-0.48, y=0.50, "Whole area", cex=1, font=2, srt=90)
par(xpd=F)


#######################################
#NC1 | Scen 2 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen2_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen2_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(2, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=2)

par(xpd=NA)
text(x=-0.3, y=0.50, "Invasive", cex=1, font=2, srt=90)
text(x=-0.48, y=0.50, "Restricted area", cex=1, font=2, srt=90)
par(xpd=F)

filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
rect(0.5,0, 1, 1, col="grey", border=NA)
rect(0,0.5, 1, 1, col="grey", border=NA)
axis(2, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=2)

par(xpd=NA)
text(x=-0.3, y=0.50, "Local", cex=1, font=2, srt=90)
text(x=-0.48, y=0.50, "Restricted area", cex=1, font=2, srt=90)
par(xpd=F)


#######################################
#NC1 | Scen 3 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen3_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen3_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(2, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=2)

par(xpd=NA)
text(x=-0.3, y=0.50, "Invasive 1", cex=1, font=2, srt=90)
text(x=-0.48, y=0.50, "Two invasions", cex=1, font=2, srt=90)
par(xpd=F)

filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(2, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=2)
axis(1, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=1)

par(xpd=NA)
text(x=-0.3, y=0.50, "Invasive 2", cex=1, font=2, srt=90)
text(x=-0.48, y=0.50, "Two invasions", cex=1, font=2, srt=90)
par(xpd=F)


#########################
#####    NC2    ######### 
#########################

#NC2 | Scen 1 | 3000gen | gamma = 0.05

NC=2
k=(24*(NC-1))
inv<- read.csv("Results/Scen1_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen1_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)

par(xpd=NA)
text(x=0.50, y=1.17, "NC2", cex=1, font=2)
text(x= 0.50, y=1.05, expression(paste(italic('K'), " = 500, ", italic('Km'), " = 10", " | ", italic('K'), " = 500, ", italic('Km'), " = 10")), cex=0.5, font=2)
par(xpd=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)

#######################################
#NC2 | Scen 2 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen2_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen2_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
rect(0.5,0, 1, 1, col="grey", border=NA)
rect(0,0.5, 1, 1, col="grey", border=NA)


#######################################
#NC2 | Scen 3 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen3_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen3_NCbis_mean.csv", header=T)

propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(1, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=1)


#########################
#####    NC3    ######### 
#########################

#NC3 | Scen 1 | 3000gen | gamma = 0.05
NC=3
k=(24*(NC-1))
inv<- read.csv("Results/Scen1_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen1_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
par(xpd=NA)
text(x=0.50, y=1.17, "NC3", cex=1, font=2)
text(x= 0.50, y=1.05, expression(paste(italic('K'), " = 500, ", italic('Km'), " = 100", " | ", italic('K'), " = 500, ", italic('Km'), " = 100")), cex=0.5, font=2)

par(xpd=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
#######################################
#NC3 | Scen 2 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen2_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen2_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
rect(0.5,0, 1, 1, col="grey", border=NA)
rect(0,0.5, 1, 1, col="grey", border=NA)


#######################################
#NC3 | Scen 3 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen3_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen3_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(1, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=1)



#########################
#####    NC4    ######### 
#########################

#NC4 | Scen 1 | 3000gen | gamma = 0.05
NC=4
k=(24*(NC-1))
inv<- read.csv("Results/Scen1_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen1_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)

par(xpd=NA)
text(x=0.50, y=1.17, "NC4", cex=1, font=2)
text(x= 0.50, y=1.05, expression(paste(italic('K'), " = 50, ", italic('Km'), " = 10", " | ", italic('K'), " = 500, ", italic('Km'), " = 100")), cex=0.5, font=2)
par(xpd=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
#######################################
#NC4 | Scen 2 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen2_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen2_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
rect(0.5,0, 1, 1, col="grey", border=NA)
rect(0,0.5, 1, 1, col="grey", border=NA)


#######################################
#NC4 | Scen 3 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen3_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen3_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(1, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=1)


#########################
#####    NC5    ######### 
#########################

#NC5 | Scen 1 | 3000gen | gamma = 0.05
NC=5
k=(24*(NC-1))
inv<- read.csv("Results/Scen1_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen1_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)

par(xpd=NA)
text(x=0.50, y=1.17, "NC5", cex=1, font=2)
text(x= 0.50, y=1.05, expression(paste(italic('K'), " = 500, ", italic('Km'), " = 100", " | ", italic('K'), " = 50, ", italic('Km'), " = 10")), cex=0.5, font=2)
par(xpd=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)

#######################################
#NC5 | Scen 2 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen2_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen2_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
rect(0.5,0, 1, 1, col="grey", border=NA)
rect(0,0.5, 1, 1, col="grey", border=NA)


#######################################
#NC5 | Scen 3 | 3000gen | gamma = 0.05#
#######################################
inv<- read.csv("Results/Scen3_NC_mean.csv", header=T)
loc<- read.csv("Results/Scen3_NCbis_mean.csv", header=T)


propI <- sapply(0:24, function(i){ round(1-inv[22+k,5+i], 3 )})
propL <- sapply(0:24, function(i){ round(loc[22+k,5+i], 3 )})

matI<-matrix(propI, nrow=5, ncol=5, byrow=T)
matI<-matI[5:1,]
matL<-matrix(propL, nrow=5, ncol=5, byrow=T)
matL<-matL[5:1,]

filled.contour4(rotate(matI), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
filled.contour4(rotate(matL), col=mycolors(100), levels=seq(0, 1, length=101), axes=F)
axis(1, at=seq(0,1, length=5), labels=seq(0,1, length=5)*100, tck=-0.01,  padj=0.5, cex.axis=0.8, mgp=c(3, 0.5, 0), lwd.ticks=1, las=1)

#####################
par(xpd=NA)
text(x=-2.02, y=-0.40, "Distance (x)", cex=1.5, font=2)
text(x=-5.73, y=3.4, "Distance (y)", cex=1.5, srt=90, font=2)
text(x=1.44, y=5.4, "Introgression\\n(%)", cex=1, font=2)

par(xpd=NA)
colorlegend2(posy = c(1.16, 5.16), posx = c(1.15, 1.45), col = mycolors(100), zlim = c(0, 100), zlevels = 11)


par(op)
dev.off()
