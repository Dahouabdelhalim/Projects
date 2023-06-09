setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology")

library(abind)
library(ncdf4)
#library(ncdf4.helpers)
library(maptools)
library(sp)
library(raster)
library(rgdal)


load("CobiaShelfMonthHistsShelfDepthWeightsBayHists.RData")
load("ShelfMonthHistsUSShelf.RData")


#########
#Calculate Ratios for whole month
#########
#to remove outlier data points that may inaccurately impact the data
#in some instances, there were very few observations at extreme temps which occurred more often than in the oceangraphic observations
#this would in turn cause untrue massive ratio values at these extreme temps
#so this now basically if temp did not represent more than 1% of observations, remove those observations
JanAll[which(JanAll<0.01)]<-0
FebAll[which(FebAll<0.01)]<-0
MarAll[which(MarAll<0.01)]<-0
AprAll[which(AprAll<0.01)]<-0
MayAll[which(MayAll<0.01)]<-0
JunAll[which(JunAll<0.01)]<-0
JulAll[which(JulAll<0.01)]<-0
AugAll[which(AugAll<0.01)]<-0
SeptAll[which(SeptAll<0.01)]<-0
OctAll[which(OctAll<0.01)]<-0
NovAll[which(NovAll<0.01)]<-0
DecAll[which(DecAll<0.01)]<-0


JanRatio<-JanAll/JanAllshelfd
JanRatio[is.na(JanRatio)]<-0
JanRatio[is.infinite(JanRatio)]<-0 #do this as well because at least for Sept there was a time when fish was in water (>30°C) but that temp did not occur in hycom habitat on US shelf during that month
FebRatio<-FebAll/FebAllshelfd
FebRatio[is.na(FebRatio)]<-0
FebRatio[is.infinite(FebRatio)]<-0
MarRatio<-MarAll/MarAllshelfd
MarRatio[is.na(MarRatio)]<-0
MarRatio[is.infinite(MarRatio)]<-0
AprRatio<-AprAll/AprAllshelfd
AprRatio[is.na(AprRatio)]<-0
AprRatio[is.infinite(AprRatio)]<-0
MayRatio<-MayAll/MayAllshelfd
MayRatio[is.na(MayRatio)]<-0
MayRatio[is.infinite(MayRatio)]<-0
JunRatio<-JunAll/JunAllshelfd
JunRatio[is.na(JunRatio)]<-0
JunRatio[is.infinite(JunRatio)]<-0
JulRatio<-JulAll/JulAllshelfd
JulRatio[is.na(JulRatio)]<-0
JulRatio[is.infinite(JulRatio)]<-0
AugRatio<-AugAll/AugAllshelfd
AugRatio[is.na(AugRatio)]<-0
AugRatio[is.infinite(AugRatio)]<-0
SeptRatio<-SeptAll/SeptAllshelfd
SeptRatio[is.na(SeptRatio)]<-0
SeptRatio[is.infinite(SeptRatio)]<-0
OctRatio<-OctAll/OctAllshelfd
OctRatio[is.na(OctRatio)]<-0
OctRatio[is.infinite(OctRatio)]<-0
NovRatio<-NovAll/NovAllshelfd
NovRatio[is.na(NovRatio)]<-0
NovRatio[is.infinite(NovRatio)]<-0
DecRatio<-DecAll/DecAllshelfd
DecRatio[is.na(DecRatio)]<-0
DecRatio[is.infinite(DecRatio)]<-0

############
#Plot Ratios
############
#first plot habitat use, habitat availability, and ratio stacked for each month
temps<-seq(from=1.5,to=33.5,by=0.5)#go to 33 not 33.5 (like I did for breaks in hists in earlier scripts cause 33.5 is really the farthest right break in hist, aka right side of last bin which is between 33 and 33.5)
par(mfrow=c(3,1))
#Jan
barplot(height=JanAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Jan Habitat Use")
barplot(height=JanAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Jan Habitat Availability")
barplot(height=JanRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Jan Ratio")
abline(h=1,lty=2)
#Feb
barplot(height=FebAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Feb Habitat Use")
barplot(height=FebAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Feb Habitat Availability")
barplot(height=FebRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Feb Ratio")
abline(h=1,lty=2)
#Mar
barplot(height=MarAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Mar Habitat Use")
barplot(height=MarAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Mar Habitat Availability")
barplot(height=MarRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Mar Ratio")
abline(h=1,lty=2)
#Apr
barplot(height=AprAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Apr Habitat Use")
barplot(height=AprAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Apr Habitat Availability")
barplot(height=AprRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Apr Ratio")
abline(h=1,lty=2)
#May
barplot(height=MayAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="May Habitat Use")
barplot(height=MayAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="May Habitat Availability")
barplot(height=MayRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="May Ratio")
abline(h=1,lty=2)
#Jun
barplot(height=JunAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Jun Habitat Use")
barplot(height=JunAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Jun Habitat Availability")
barplot(height=JunRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Jun Ratio")
abline(h=1,lty=2)
#Jul
barplot(height=JulAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Jul Habitat Use")
barplot(height=JulAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Jul Habitat Availability")
barplot(height=JulRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Jul Ratio")
abline(h=1,lty=2)
#Aug
barplot(height=AugAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Aug Habitat Use")
barplot(height=AugAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Aug Habitat Availability")
barplot(height=AugRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Aug Ratio")
abline(h=1,lty=2)
#Sept
barplot(height=SeptAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Sept Habitat Use")
barplot(height=SeptAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Sept Habitat Availability")
barplot(height=SeptRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Sept Ratio")
abline(h=1,lty=2)
#Oct
barplot(height=OctAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Oct Habitat Use")
barplot(height=OctAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Oct Habitat Availability")
barplot(height=OctRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Oct Ratio")
abline(h=1,lty=2)
#Nov
barplot(height=NovAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Nov Habitat Use")
barplot(height=NovAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Nov Habitat Availability")
barplot(height=NovRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Nov Ratio")
abline(h=1,lty=2)
#Dec
barplot(height=DecAll[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Dec Habitat Use")
barplot(height=DecAllshelfd[18:62], names.arg=temps[18:62],ylim=c(0,0.5),ylab="Densities",main="Dec Habitat Availability")
barplot(height=DecRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="Densities",main="Dec Ratio")
abline(h=1,lty=2)



#Plot just Ratios of each month
par(mfrow = c(3,4),oma = c(5,3,0,3) + 0.1,mar = c(0,2,3.5,2.5) + 0.1)
barplot(height=JanRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(JanAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n",lwd=1.5, axes=F, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=F)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(JanAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=2,"Ratios",line=3,cex=1.25)
mtext(side=3,"January",cex=1.25)

barplot(height=FebRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(FebAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n",axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=F)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(FebAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=3,"February",cex=1.25)

barplot(height=MarRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(MarAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=F)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(MarAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=3,"March",cex=1.25)

barplot(height=AprRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(AprAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=F)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(AprAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=4,"Habitat Densities",line=3.5,cex=1.25)
mtext(side=3,"April",cex=1.25)

barplot(height=MayRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(MayAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=F)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(MayAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=2,"Ratios",line=3,cex=1.25)
mtext(side=3,"May",cex=1.25)

barplot(height=JunRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(JunAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=F)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(JunAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=3,"June",cex=1.25)

barplot(height=JulRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(JulAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=F)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(JulAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=3,"July",cex=1.25)

barplot(height=AugRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(AugAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=F)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(AugAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=4,"Habitat Densities",line=3.5,cex=1.25)
mtext(side=3,"August",cex=1.25)

barplot(height=SeptRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(SeptAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=T)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(SeptAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=2,"Ratios",line=3,cex=1.25)
mtext(side=1,"Temperature (°C)",line=3,cex=1.25)
mtext(side=3,"September",cex=1.25)

barplot(height=OctRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(OctAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=T)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(OctAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=1,"Temperature (°C)",line=3,cex=1.25)
mtext(side=3,"October",cex=1.25)

barplot(height=NovRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(NovAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=T)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(NovAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=1,"Temperature (°C)",line=3,cex=1.25)
mtext(side=3,"November",cex=1.25)

barplot(height=DecRatio[18:62], names.arg=temps[18:62],ylim=c(0,18),ylab="",xaxt="n")
abline(h=1,lty=2)
par(new=TRUE)
plot(DecAll[18:62]~temps[18:62],type="l",ylim=c(0,0.35),xaxt="n", axes=F,lwd=1.5, xlab="",ylab="",yaxs="i",col="coral1") #yaxs="i" makes yaxis 0 start at bottom like barplot
axis(side=4)
axis(side=1,labels=T)
#mtext(side=4,line=2.5,"Habitat Densities")
lines(DecAllshelfd[18:62]~temps[18:62],type="l",yaxs="i",col="light blue",lwd=1.5)
mtext(side=4,"Habitat Densities",line=3.5,cex=1.25)
mtext(side=1,"Temperature (°C)",line=3,cex=1.25)
mtext(side=3,"December",cex=1.25)


#############
#Plot Depth Weighting Factors
#############
depthnames<-c("0-20","20-40","40-60","60-80","80-100","100-250")
par(mfrow = c(3,4),oma = c(5,3,0,3) + 0.1,mar = c(1.5,4,2.5,1) + 0.1)
barplot(height=rev(JanDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="January",xaxt="n")
axis(side=1,labels=F)
mtext("Depth bins (m)",side=2,line=4.5)
barplot(height=rev(FebDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="February",xaxt="n")
axis(side=1,labels=F)
barplot(height=rev(MarDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="March",xaxt="n")
axis(side=1,labels=F)
barplot(height=rev(AprDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="April",xaxt="n")
axis(side=1,labels=F)
barplot(height=rev(MayDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="May",xaxt="n")
mtext("Depth bins (m)",side=2,line=4.5)
axis(side=1,labels=F)
barplot(height=rev(JunDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="June",xaxt="n")
axis(side=1,labels=F)
barplot(height=rev(JulDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="July",xaxt="n")
axis(side=1,labels=F)
barplot(height=rev(AugDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="August",xaxt="n")
axis(side=1,labels=F)
barplot(height=rev(SeptDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="September")
mtext("Proportion used",side=1,line=3)
mtext("Depth bins (m)",side=2,line=4.5)
barplot(height=rev(OctDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="October")
mtext("Proportion used",side=1,line=3)
barplot(height=rev(NovDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="November")
mtext("Proportion used",side=1,line=3)
barplot(height=rev(DecDepthWeight),names.arg=rev(depthnames),horiz=TRUE,xlim=c(0,1),las=2,main="December")
mtext("Proportion used",side=1,line=3)



##############
#Compare temps of climatology and extreme years
##############
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology")
load("ShelfMonthClimatology.RData") #noticed that anywhere where its 10p should be 100p
load("ShelfMonthExtremeYears.RData")
#climatology
Jan_stack<-stack((brick(c(Jan_array_depth_sum$temp20_clim$raster,Jan_array_depth_sum$temp40_clim$raster,
                          Jan_array_depth_sum$temp60_clim$raster,Jan_array_depth_sum$temp80_clim$raster,
                          Jan_array_depth_sum$temp100_clim$raster,Jan_array_depth_sum$temp10p_clim$raster))))
Feb_stack<-stack((brick(c(Feb_array_depth_sum$temp20_clim$raster,Feb_array_depth_sum$temp40_clim$raster,
                          Feb_array_depth_sum$temp60_clim$raster,Feb_array_depth_sum$temp80_clim$raster,
                          Feb_array_depth_sum$temp100_clim$raster,Feb_array_depth_sum$temp10p_clim$raster))))
Mar_stack<-stack((brick(c(Mar_array_depth_sum$temp20_clim$raster,Mar_array_depth_sum$temp40_clim$raster,
                          Mar_array_depth_sum$temp60_clim$raster,Mar_array_depth_sum$temp80_clim$raster,
                          Mar_array_depth_sum$temp100_clim$raster,Mar_array_depth_sum$temp10p_clim$raster))))
Apr_stack<-stack((brick(c(Apr_array_depth_sum$temp20_clim$raster,Apr_array_depth_sum$temp40_clim$raster,
                          Apr_array_depth_sum$temp60_clim$raster,Apr_array_depth_sum$temp80_clim$raster,
                          Apr_array_depth_sum$temp100_clim$raster,Apr_array_depth_sum$temp10p_clim$raster))))
May_stack<-stack((brick(c(May_array_depth_sum$temp20_clim$raster,May_array_depth_sum$temp40_clim$raster,
                          May_array_depth_sum$temp60_clim$raster,May_array_depth_sum$temp80_clim$raster,
                          May_array_depth_sum$temp100_clim$raster,May_array_depth_sum$temp10p_clim$raster))))
Jun_stack<-stack((brick(c(Jun_array_depth_sum$temp20_clim$raster,Jun_array_depth_sum$temp40_clim$raster,
                          Jun_array_depth_sum$temp60_clim$raster,Jun_array_depth_sum$temp80_clim$raster,
                          Jun_array_depth_sum$temp100_clim$raster,Jun_array_depth_sum$temp10p_clim$raster))))
Jul_stack<-stack((brick(c(Jul_array_depth_sum$temp20_clim$raster,Jul_array_depth_sum$temp40_clim$raster,
                          Jul_array_depth_sum$temp60_clim$raster,Jul_array_depth_sum$temp80_clim$raster,
                          Jul_array_depth_sum$temp100_clim$raster,Jul_array_depth_sum$temp10p_clim$raster))))
Aug_stack<-stack((brick(c(Aug_array_depth_sum$temp20_clim$raster,Aug_array_depth_sum$temp40_clim$raster,
                          Aug_array_depth_sum$temp60_clim$raster,Aug_array_depth_sum$temp80_clim$raster,
                          Aug_array_depth_sum$temp100_clim$raster,Aug_array_depth_sum$temp10p_clim$raster))))
Sept_stack<-stack((brick(c(Sept_array_depth_sum$temp20_clim$raster,Sept_array_depth_sum$temp40_clim$raster,
                           Sept_array_depth_sum$temp60_clim$raster,Sept_array_depth_sum$temp80_clim$raster,
                           Sept_array_depth_sum$temp100_clim$raster,Sept_array_depth_sum$temp10p_clim$raster))))
Oct_stack<-stack((brick(c(Oct_array_depth_sum$temp20_clim$raster,Oct_array_depth_sum$temp40_clim$raster,
                          Oct_array_depth_sum$temp60_clim$raster,Oct_array_depth_sum$temp80_clim$raster,
                          Oct_array_depth_sum$temp100_clim$raster,Oct_array_depth_sum$temp10p_clim$raster))))
Nov_stack<-stack((brick(c(Nov_array_depth_sum$temp20_clim$raster,Nov_array_depth_sum$temp40_clim$raster,
                          Nov_array_depth_sum$temp60_clim$raster,Nov_array_depth_sum$temp80_clim$raster,
                          Nov_array_depth_sum$temp100_clim$raster,Nov_array_depth_sum$temp10p_clim$raster))))
Dec_stack<-stack((brick(c(Dec_array_depth_sum$temp20_clim$raster,Dec_array_depth_sum$temp40_clim$raster,
                          Dec_array_depth_sum$temp60_clim$raster,Dec_array_depth_sum$temp80_clim$raster,
                          Dec_array_depth_sum$temp100_clim$raster,Dec_array_depth_sum$temp10p_clim$raster))))

#extreme years
ex_array201201_stack<-stack((brick(c(ex_array201201_depth_sum$temp20_clim$raster,ex_array201201_depth_sum$temp40_clim$raster,
                                     ex_array201201_depth_sum$temp60_clim$raster,ex_array201201_depth_sum$temp80_clim$raster,
                                     ex_array201201_depth_sum$temp100_clim$raster,ex_array201201_depth_sum$temp100p_clim$raster))))
ex_array201202_stack<-stack((brick(c(ex_array201202_depth_sum$temp20_clim$raster,ex_array201202_depth_sum$temp40_clim$raster,
                                     ex_array201202_depth_sum$temp60_clim$raster,ex_array201202_depth_sum$temp80_clim$raster,
                                     ex_array201202_depth_sum$temp100_clim$raster,ex_array201202_depth_sum$temp100p_clim$raster))))
ex_array201203_stack<-stack((brick(c(ex_array201203_depth_sum$temp20_clim$raster,ex_array201203_depth_sum$temp40_clim$raster,
                                     ex_array201203_depth_sum$temp60_clim$raster,ex_array201203_depth_sum$temp80_clim$raster,
                                     ex_array201203_depth_sum$temp100_clim$raster,ex_array201203_depth_sum$temp100p_clim$raster))))
ex_array201204_stack<-stack((brick(c(ex_array201204_depth_sum$temp20_clim$raster,ex_array201204_depth_sum$temp40_clim$raster,
                                     ex_array201204_depth_sum$temp60_clim$raster,ex_array201204_depth_sum$temp80_clim$raster,
                                     ex_array201204_depth_sum$temp100_clim$raster,ex_array201204_depth_sum$temp100p_clim$raster))))
ex_array201205_stack<-stack((brick(c(ex_array201205_depth_sum$temp20_clim$raster,ex_array201205_depth_sum$temp40_clim$raster,
                                     ex_array201205_depth_sum$temp60_clim$raster,ex_array201205_depth_sum$temp80_clim$raster,
                                     ex_array201205_depth_sum$temp100_clim$raster,ex_array201205_depth_sum$temp100p_clim$raster))))
ex_array201206_stack<-stack((brick(c(ex_array201206_depth_sum$temp20_clim$raster,ex_array201206_depth_sum$temp40_clim$raster,
                                     ex_array201206_depth_sum$temp60_clim$raster,ex_array201206_depth_sum$temp80_clim$raster,
                                     ex_array201206_depth_sum$temp100_clim$raster,ex_array201206_depth_sum$temp100p_clim$raster))))
ex_array201207_stack<-stack((brick(c(ex_array201207_depth_sum$temp20_clim$raster,ex_array201207_depth_sum$temp40_clim$raster,
                                     ex_array201207_depth_sum$temp60_clim$raster,ex_array201207_depth_sum$temp80_clim$raster,
                                     ex_array201207_depth_sum$temp100_clim$raster,ex_array201207_depth_sum$temp100p_clim$raster))))
ex_array201208_stack<-stack((brick(c(ex_array201208_depth_sum$temp20_clim$raster,ex_array201208_depth_sum$temp40_clim$raster,
                                     ex_array201208_depth_sum$temp60_clim$raster,ex_array201208_depth_sum$temp80_clim$raster,
                                     ex_array201208_depth_sum$temp100_clim$raster,ex_array201208_depth_sum$temp100p_clim$raster))))
ex_array201209_stack<-stack((brick(c(ex_array201209_depth_sum$temp20_clim$raster,ex_array201209_depth_sum$temp40_clim$raster,
                                     ex_array201209_depth_sum$temp60_clim$raster,ex_array201209_depth_sum$temp80_clim$raster,
                                     ex_array201209_depth_sum$temp100_clim$raster,ex_array201209_depth_sum$temp100p_clim$raster))))
ex_array201210_stack<-stack((brick(c(ex_array201210_depth_sum$temp20_clim$raster,ex_array201210_depth_sum$temp40_clim$raster,
                                     ex_array201210_depth_sum$temp60_clim$raster,ex_array201210_depth_sum$temp80_clim$raster,
                                     ex_array201210_depth_sum$temp100_clim$raster,ex_array201210_depth_sum$temp100p_clim$raster))))
ex_array201211_stack<-stack((brick(c(ex_array201211_depth_sum$temp20_clim$raster,ex_array201211_depth_sum$temp40_clim$raster,
                                     ex_array201211_depth_sum$temp60_clim$raster,ex_array201211_depth_sum$temp80_clim$raster,
                                     ex_array201211_depth_sum$temp100_clim$raster,ex_array201211_depth_sum$temp100p_clim$raster))))
ex_array201212_stack<-stack((brick(c(ex_array201212_depth_sum$temp20_clim$raster,ex_array201212_depth_sum$temp40_clim$raster,
                                     ex_array201212_depth_sum$temp60_clim$raster,ex_array201212_depth_sum$temp80_clim$raster,
                                     ex_array201212_depth_sum$temp100_clim$raster,ex_array201212_depth_sum$temp100p_clim$raster))))

ex_array199601_stack<-stack((brick(c(ex_array199601_depth_sum$temp20_clim$raster,ex_array199601_depth_sum$temp40_clim$raster,
                                     ex_array199601_depth_sum$temp60_clim$raster,ex_array199601_depth_sum$temp80_clim$raster,
                                     ex_array199601_depth_sum$temp100_clim$raster,ex_array199601_depth_sum$temp100p_clim$raster))))
ex_array199602_stack<-stack((brick(c(ex_array199602_depth_sum$temp20_clim$raster,ex_array199602_depth_sum$temp40_clim$raster,
                                     ex_array199602_depth_sum$temp60_clim$raster,ex_array199602_depth_sum$temp80_clim$raster,
                                     ex_array199602_depth_sum$temp100_clim$raster,ex_array199602_depth_sum$temp100p_clim$raster))))
ex_array199603_stack<-stack((brick(c(ex_array199603_depth_sum$temp20_clim$raster,ex_array199603_depth_sum$temp40_clim$raster,
                                     ex_array199603_depth_sum$temp60_clim$raster,ex_array199603_depth_sum$temp80_clim$raster,
                                     ex_array199603_depth_sum$temp100_clim$raster,ex_array199603_depth_sum$temp100p_clim$raster))))
ex_array199604_stack<-stack((brick(c(ex_array199604_depth_sum$temp20_clim$raster,ex_array199604_depth_sum$temp40_clim$raster,
                                     ex_array199604_depth_sum$temp60_clim$raster,ex_array199604_depth_sum$temp80_clim$raster,
                                     ex_array199604_depth_sum$temp100_clim$raster,ex_array199604_depth_sum$temp100p_clim$raster))))
ex_array199605_stack<-stack((brick(c(ex_array199605_depth_sum$temp20_clim$raster,ex_array199605_depth_sum$temp40_clim$raster,
                                     ex_array199605_depth_sum$temp60_clim$raster,ex_array199605_depth_sum$temp80_clim$raster,
                                     ex_array199605_depth_sum$temp100_clim$raster,ex_array199605_depth_sum$temp100p_clim$raster))))
ex_array199606_stack<-stack((brick(c(ex_array199606_depth_sum$temp20_clim$raster,ex_array199606_depth_sum$temp40_clim$raster,
                                     ex_array199606_depth_sum$temp60_clim$raster,ex_array199606_depth_sum$temp80_clim$raster,
                                     ex_array199606_depth_sum$temp100_clim$raster,ex_array199606_depth_sum$temp100p_clim$raster))))
ex_array199607_stack<-stack((brick(c(ex_array199607_depth_sum$temp20_clim$raster,ex_array199607_depth_sum$temp40_clim$raster,
                                     ex_array199607_depth_sum$temp60_clim$raster,ex_array199607_depth_sum$temp80_clim$raster,
                                     ex_array199607_depth_sum$temp100_clim$raster,ex_array199607_depth_sum$temp100p_clim$raster))))
ex_array199608_stack<-stack((brick(c(ex_array199608_depth_sum$temp20_clim$raster,ex_array199608_depth_sum$temp40_clim$raster,
                                     ex_array199608_depth_sum$temp60_clim$raster,ex_array199608_depth_sum$temp80_clim$raster,
                                     ex_array199608_depth_sum$temp100_clim$raster,ex_array199608_depth_sum$temp100p_clim$raster))))
ex_array199609_stack<-stack((brick(c(ex_array199609_depth_sum$temp20_clim$raster,ex_array199609_depth_sum$temp40_clim$raster,
                                     ex_array199609_depth_sum$temp60_clim$raster,ex_array199609_depth_sum$temp80_clim$raster,
                                     ex_array199609_depth_sum$temp100_clim$raster,ex_array199609_depth_sum$temp100p_clim$raster))))
ex_array199610_stack<-stack((brick(c(ex_array199610_depth_sum$temp20_clim$raster,ex_array199610_depth_sum$temp40_clim$raster,
                                     ex_array199610_depth_sum$temp60_clim$raster,ex_array199610_depth_sum$temp80_clim$raster,
                                     ex_array199610_depth_sum$temp100_clim$raster,ex_array199610_depth_sum$temp100p_clim$raster))))
ex_array199611_stack<-stack((brick(c(ex_array199611_depth_sum$temp20_clim$raster,ex_array199611_depth_sum$temp40_clim$raster,
                                     ex_array199611_depth_sum$temp60_clim$raster,ex_array199611_depth_sum$temp80_clim$raster,
                                     ex_array199611_depth_sum$temp100_clim$raster,ex_array199611_depth_sum$temp100p_clim$raster))))
ex_array199612_stack<-stack((brick(c(ex_array199612_depth_sum$temp20_clim$raster,ex_array199612_depth_sum$temp40_clim$raster,
                                     ex_array199612_depth_sum$temp60_clim$raster,ex_array199612_depth_sum$temp80_clim$raster,
                                     ex_array199612_depth_sum$temp100_clim$raster,ex_array199612_depth_sum$temp100p_clim$raster))))

#mask for US waters
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/FederalAndStateWaters/GIS_shp")
uswaters<-readOGR(".","US_State_Fed_waters")
uswaters<-spTransform(uswaters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#climatology
Jan_stack<-mask(Jan_stack,mask=uswaters)
Feb_stack<-mask(Feb_stack,mask=uswaters)
Mar_stack<-mask(Mar_stack,mask=uswaters)
Apr_stack<-mask(Apr_stack,mask=uswaters)
May_stack<-mask(May_stack,mask=uswaters)
Jun_stack<-mask(Jun_stack,mask=uswaters)
Jul_stack<-mask(Jul_stack,mask=uswaters)
Aug_stack<-mask(Aug_stack,mask=uswaters)
Sept_stack<-mask(Sept_stack,mask=uswaters)
Oct_stack<-mask(Oct_stack,mask=uswaters)
Nov_stack<-mask(Nov_stack,mask=uswaters)
Dec_stack<-mask(Dec_stack,mask=uswaters)

#extreme years
ex_array199601_stack<-mask(ex_array199601_stack,mask=uswaters)
ex_array199602_stack<-mask(ex_array199602_stack,mask=uswaters)
ex_array199603_stack<-mask(ex_array199603_stack,mask=uswaters)
ex_array199604_stack<-mask(ex_array199604_stack,mask=uswaters)
ex_array199605_stack<-mask(ex_array199605_stack,mask=uswaters)
ex_array199606_stack<-mask(ex_array199606_stack,mask=uswaters)
ex_array199607_stack<-mask(ex_array199607_stack,mask=uswaters)
ex_array199608_stack<-mask(ex_array199608_stack,mask=uswaters)
ex_array199609_stack<-mask(ex_array199609_stack,mask=uswaters)
ex_array199610_stack<-mask(ex_array199610_stack,mask=uswaters)
ex_array199611_stack<-mask(ex_array199611_stack,mask=uswaters)
ex_array199612_stack<-mask(ex_array199612_stack,mask=uswaters)
ex_array201201_stack<-mask(ex_array201201_stack,mask=uswaters)
ex_array201202_stack<-mask(ex_array201202_stack,mask=uswaters)
ex_array201203_stack<-mask(ex_array201203_stack,mask=uswaters)
ex_array201204_stack<-mask(ex_array201204_stack,mask=uswaters)
ex_array201205_stack<-mask(ex_array201205_stack,mask=uswaters)
ex_array201206_stack<-mask(ex_array201206_stack,mask=uswaters)
ex_array201207_stack<-mask(ex_array201207_stack,mask=uswaters)
ex_array201208_stack<-mask(ex_array201208_stack,mask=uswaters)
ex_array201209_stack<-mask(ex_array201209_stack,mask=uswaters)
ex_array201210_stack<-mask(ex_array201210_stack,mask=uswaters)
ex_array201211_stack<-mask(ex_array201211_stack,mask=uswaters)
ex_array201212_stack<-mask(ex_array201212_stack,mask=uswaters)

#read in future years temps
setwd("/Users/danielcrear/Documents/PhD Project/SDM/SDM/Regional_climatology/Climate_Deltas")
load("ShelfMonthClimatologyFuture.RData")
#need to mask with us waters before can take means
Jan_0_20_future_mask<-mask(Jan_0_20_future,mask=uswaters)
Jan_20_40_future_mask<-mask(Jan_20_40_future,mask=uswaters)
Jan_40_60_future_mask<-mask(Jan_40_60_future,mask=uswaters)
Jan_60_80_future_mask<-mask(Jan_60_80_future,mask=uswaters)
Feb_0_20_future_mask<-mask(Feb_0_20_future,mask=uswaters)
Feb_20_40_future_mask<-mask(Feb_20_40_future,mask=uswaters)
Feb_40_60_future_mask<-mask(Feb_40_60_future,mask=uswaters)
Feb_60_80_future_mask<-mask(Feb_60_80_future,mask=uswaters)
Mar_0_20_future_mask<-mask(Mar_0_20_future,mask=uswaters)
Mar_20_40_future_mask<-mask(Mar_20_40_future,mask=uswaters)
Mar_40_60_future_mask<-mask(Mar_40_60_future,mask=uswaters)
Mar_60_80_future_mask<-mask(Mar_60_80_future,mask=uswaters)
Apr_0_20_future_mask<-mask(Apr_0_20_future,mask=uswaters)
Apr_20_40_future_mask<-mask(Apr_20_40_future,mask=uswaters)
Apr_40_60_future_mask<-mask(Apr_40_60_future,mask=uswaters)
Apr_60_80_future_mask<-mask(Apr_60_80_future,mask=uswaters)
May_0_20_future_mask<-mask(May_0_20_future,mask=uswaters)
May_20_40_future_mask<-mask(May_20_40_future,mask=uswaters)
May_40_60_future_mask<-mask(May_40_60_future,mask=uswaters)
May_60_80_future_mask<-mask(May_60_80_future,mask=uswaters)
Jun_0_20_future_mask<-mask(Jun_0_20_future,mask=uswaters)
Jun_20_40_future_mask<-mask(Jun_20_40_future,mask=uswaters)
Jun_40_60_future_mask<-mask(Jun_40_60_future,mask=uswaters)
Jun_60_80_future_mask<-mask(Jun_60_80_future,mask=uswaters)
Jul_0_20_future_mask<-mask(Jul_0_20_future,mask=uswaters)
Jul_20_40_future_mask<-mask(Jul_20_40_future,mask=uswaters)
Jul_40_60_future_mask<-mask(Jul_40_60_future,mask=uswaters)
Jul_60_80_future_mask<-mask(Jul_60_80_future,mask=uswaters)
Aug_0_20_future_mask<-mask(Aug_0_20_future,mask=uswaters)
Aug_20_40_future_mask<-mask(Aug_20_40_future,mask=uswaters)
Aug_40_60_future_mask<-mask(Aug_40_60_future,mask=uswaters)
Aug_60_80_future_mask<-mask(Aug_60_80_future,mask=uswaters)
Sept_0_20_future_mask<-mask(Sept_0_20_future,mask=uswaters)
Sept_20_40_future_mask<-mask(Sept_20_40_future,mask=uswaters)
Sept_40_60_future_mask<-mask(Sept_40_60_future,mask=uswaters)
Sept_60_80_future_mask<-mask(Sept_60_80_future,mask=uswaters)
Oct_0_20_future_mask<-mask(Oct_0_20_future,mask=uswaters)
Oct_20_40_future_mask<-mask(Oct_20_40_future,mask=uswaters)
Oct_40_60_future_mask<-mask(Oct_40_60_future,mask=uswaters)
Oct_60_80_future_mask<-mask(Oct_60_80_future,mask=uswaters)
Nov_0_20_future_mask<-mask(Nov_0_20_future,mask=uswaters)
Nov_20_40_future_mask<-mask(Nov_20_40_future,mask=uswaters)
Nov_40_60_future_mask<-mask(Nov_40_60_future,mask=uswaters)
Nov_60_80_future_mask<-mask(Nov_60_80_future,mask=uswaters)
Dec_0_20_future_mask<-mask(Dec_0_20_future,mask=uswaters)
Dec_20_40_future_mask<-mask(Dec_20_40_future,mask=uswaters)
Dec_40_60_future_mask<-mask(Dec_40_60_future,mask=uswaters)
Dec_60_80_future_mask<-mask(Dec_60_80_future,mask=uswaters)

setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/Climate_Deltas/ExtremeFutureTemps")
load("future70extremetemps.RData")
Jan_70_future_mask<-mask(future70extremetempsJan, mask=uswaters)
Feb_70_future_mask<-mask(future70extremetempsFeb, mask=uswaters)
Mar_70_future_mask<-mask(future70extremetempsMar, mask=uswaters)
Apr_70_future_mask<-mask(future70extremetempsApr, mask=uswaters)
May_70_future_mask<-mask(future70extremetempsMay, mask=uswaters)
Jun_70_future_mask<-mask(future70extremetempsJun, mask=uswaters)
Jul_70_future_mask<-mask(future70extremetempsJul, mask=uswaters)
Aug_70_future_mask<-mask(future70extremetempsAug, mask=uswaters)
Sept_70_future_mask<-mask(future70extremetempsSept, mask=uswaters)
Oct_70_future_mask<-mask(future70extremetempsOct, mask=uswaters)
Nov_70_future_mask<-mask(future70extremetempsNov, mask=uswaters)
Dec_70_future_mask<-mask(future70extremetempsDec, mask=uswaters)


#take mean of temps
clim_temps20<-c(mean(values(Jan_stack$Jan_20),na.rm=T),mean(values(Feb_stack$Feb_20),na.rm=T),mean(values(Mar_stack$Mar_20),na.rm=T),
                mean(values(Apr_stack$Apr_20),na.rm=T),mean(values(May_stack$May_20),na.rm=T),mean(values(Jun_stack$Jun_20),na.rm=T),
                mean(values(Jul_stack$Jul_20),na.rm=T),mean(values(Aug_stack$Aug_20),na.rm=T),mean(values(Sept_stack$Sept_20),na.rm=T),
                mean(values(Oct_stack$Oct_20),na.rm=T),mean(values(Nov_stack$Nov_20),na.rm=T),mean(values(Dec_stack$Dec_20),na.rm=T))
clim_temps40<-c(mean(values(Jan_stack$Jan_40),na.rm=T),mean(values(Feb_stack$Feb_40),na.rm=T),mean(values(Mar_stack$Mar_40),na.rm=T),
                mean(values(Apr_stack$Apr_40),na.rm=T),mean(values(May_stack$May_40),na.rm=T),mean(values(Jun_stack$Jun_40),na.rm=T),
                mean(values(Jul_stack$Jul_40),na.rm=T),mean(values(Aug_stack$Aug_40),na.rm=T),mean(values(Sept_stack$Sept_40),na.rm=T),
                mean(values(Oct_stack$Oct_40),na.rm=T),mean(values(Nov_stack$Nov_40),na.rm=T),mean(values(Dec_stack$Dec_40),na.rm=T))
clim_temps60<-c(mean(values(Jan_stack$Jan_60),na.rm=T),mean(values(Feb_stack$Feb_60),na.rm=T),mean(values(Mar_stack$Mar_60),na.rm=T),
                mean(values(Apr_stack$Apr_60),na.rm=T),mean(values(May_stack$May_60),na.rm=T),mean(values(Jun_stack$Jun_60),na.rm=T),
                mean(values(Jul_stack$Jul_60),na.rm=T),mean(values(Aug_stack$Aug_60),na.rm=T),mean(values(Sept_stack$Sept_60),na.rm=T),
                mean(values(Oct_stack$Oct_60),na.rm=T),mean(values(Nov_stack$Nov_60),na.rm=T),mean(values(Dec_stack$Dec_60),na.rm=T))
clim_temps80<-c(mean(values(Jan_stack$Jan_80),na.rm=T),mean(values(Feb_stack$Feb_80),na.rm=T),mean(values(Mar_stack$Mar_80),na.rm=T),
                mean(values(Apr_stack$Apr_80),na.rm=T),mean(values(May_stack$May_80),na.rm=T),mean(values(Jun_stack$Jun_80),na.rm=T),
                mean(values(Jul_stack$Jul_80),na.rm=T),mean(values(Aug_stack$Aug_80),na.rm=T),mean(values(Sept_stack$Sept_80),na.rm=T),
                mean(values(Oct_stack$Oct_80),na.rm=T),mean(values(Nov_stack$Nov_80),na.rm=T),mean(values(Dec_stack$Dec_80),na.rm=T))
clim_temps100<-c(mean(values(Jan_stack$Jan_100),na.rm=T),mean(values(Feb_stack$Feb_100),na.rm=T),mean(values(Mar_stack$Mar_100),na.rm=T),
                 mean(values(Apr_stack$Apr_100),na.rm=T),mean(values(May_stack$May_100),na.rm=T),mean(values(Jun_stack$Jun_100),na.rm=T),
                 mean(values(Jul_stack$Jul_100),na.rm=T),mean(values(Aug_stack$Aug_100),na.rm=T),mean(values(Sept_stack$Sept_100),na.rm=T),
                 mean(values(Oct_stack$Oct_100),na.rm=T),mean(values(Nov_stack$Nov_100),na.rm=T),mean(values(Dec_stack$Dec_100),na.rm=T))
clim_temps250<-c(mean(values(Jan_stack$Jan_250),na.rm=T),mean(values(Feb_stack$Feb_250),na.rm=T),mean(values(Mar_stack$Mar_250),na.rm=T),
                 mean(values(Apr_stack$Apr_250),na.rm=T),mean(values(May_stack$May_250),na.rm=T),mean(values(Jun_stack$Jun_250),na.rm=T),
                 mean(values(Jul_stack$Jul_250),na.rm=T),mean(values(Aug_stack$Aug_250),na.rm=T),mean(values(Sept_stack$Sept_250),na.rm=T),
                 mean(values(Oct_stack$Oct_250),na.rm=T),mean(values(Nov_stack$Nov_250),na.rm=T),mean(values(Dec_stack$Dec_250),na.rm=T))

t2012_temps20<-c(mean(values(ex_array201201_stack$X201201_20),na.rm=T),mean(values(ex_array201202_stack$X201202_20),na.rm=T),mean(values(ex_array201203_stack$X201203_20),na.rm=T),
                 mean(values(ex_array201204_stack$X201204_20),na.rm=T),mean(values(ex_array201205_stack$X201205_20),na.rm=T),mean(values(ex_array201206_stack$X201206_20),na.rm=T),
                 mean(values(ex_array201207_stack$X201207_20),na.rm=T),mean(values(ex_array201208_stack$X201208_20),na.rm=T),mean(values(ex_array201209_stack$X201209_20),na.rm=T),
                 mean(values(ex_array201210_stack$X201210_20),na.rm=T),mean(values(ex_array201211_stack$X201211_20),na.rm=T),mean(values(ex_array201212_stack$X201212_20),na.rm=T))
t2012_temps40<-c(mean(values(ex_array201201_stack$X201201_40),na.rm=T),mean(values(ex_array201202_stack$X201202_40),na.rm=T),mean(values(ex_array201203_stack$X201203_40),na.rm=T),
                 mean(values(ex_array201204_stack$X201204_40),na.rm=T),mean(values(ex_array201205_stack$X201205_40),na.rm=T),mean(values(ex_array201206_stack$X201206_40),na.rm=T),
                 mean(values(ex_array201207_stack$X201207_40),na.rm=T),mean(values(ex_array201208_stack$X201208_40),na.rm=T),mean(values(ex_array201209_stack$X201209_40),na.rm=T),
                 mean(values(ex_array201210_stack$X201210_40),na.rm=T),mean(values(ex_array201211_stack$X201211_40),na.rm=T),mean(values(ex_array201212_stack$X201212_40),na.rm=T))
t2012_temps60<-c(mean(values(ex_array201201_stack$X201201_60),na.rm=T),mean(values(ex_array201202_stack$X201202_60),na.rm=T),mean(values(ex_array201203_stack$X201203_60),na.rm=T),
                 mean(values(ex_array201204_stack$X201204_60),na.rm=T),mean(values(ex_array201205_stack$X201205_60),na.rm=T),mean(values(ex_array201206_stack$X201206_60),na.rm=T),
                 mean(values(ex_array201207_stack$X201207_60),na.rm=T),mean(values(ex_array201208_stack$X201208_60),na.rm=T),mean(values(ex_array201209_stack$X201209_60),na.rm=T),
                 mean(values(ex_array201210_stack$X201210_60),na.rm=T),mean(values(ex_array201211_stack$X201211_60),na.rm=T),mean(values(ex_array201212_stack$X201212_60),na.rm=T))
t2012_temps80<-c(mean(values(ex_array201201_stack$X201201_80),na.rm=T),mean(values(ex_array201202_stack$X201202_80),na.rm=T),mean(values(ex_array201203_stack$X201203_80),na.rm=T),
                 mean(values(ex_array201204_stack$X201204_80),na.rm=T),mean(values(ex_array201205_stack$X201205_80),na.rm=T),mean(values(ex_array201206_stack$X201206_80),na.rm=T),
                 mean(values(ex_array201207_stack$X201207_80),na.rm=T),mean(values(ex_array201208_stack$X201208_80),na.rm=T),mean(values(ex_array201209_stack$X201209_80),na.rm=T),
                 mean(values(ex_array201210_stack$X201210_80),na.rm=T),mean(values(ex_array201211_stack$X201211_80),na.rm=T),mean(values(ex_array201212_stack$X201212_80),na.rm=T))
t2012_temps100<-c(mean(values(ex_array201201_stack$X201201_100),na.rm=T),mean(values(ex_array201202_stack$X201202_100),na.rm=T),mean(values(ex_array201203_stack$X201203_100),na.rm=T),
                  mean(values(ex_array201204_stack$X201204_100),na.rm=T),mean(values(ex_array201205_stack$X201205_100),na.rm=T),mean(values(ex_array201206_stack$X201206_100),na.rm=T),
                  mean(values(ex_array201207_stack$X201207_100),na.rm=T),mean(values(ex_array201208_stack$X201208_100),na.rm=T),mean(values(ex_array201209_stack$X201209_100),na.rm=T),
                  mean(values(ex_array201210_stack$X201210_100),na.rm=T),mean(values(ex_array201211_stack$X201211_100),na.rm=T),mean(values(ex_array201212_stack$X201212_100),na.rm=T))
t2012_temps250<-c(mean(values(ex_array201201_stack$X201201_250),na.rm=T),mean(values(ex_array201202_stack$X201202_250),na.rm=T),mean(values(ex_array201203_stack$X201203_250),na.rm=T),
                  mean(values(ex_array201204_stack$X201204_250),na.rm=T),mean(values(ex_array201205_stack$X201205_250),na.rm=T),mean(values(ex_array201206_stack$X201206_250),na.rm=T),
                  mean(values(ex_array201207_stack$X201207_250),na.rm=T),mean(values(ex_array201208_stack$X201208_250),na.rm=T),mean(values(ex_array201209_stack$X201209_250),na.rm=T),
                  mean(values(ex_array201210_stack$X201210_250),na.rm=T),mean(values(ex_array201211_stack$X201211_250),na.rm=T),mean(values(ex_array201212_stack$X201212_250),na.rm=T))

t1996_temps20<-c(mean(values(ex_array199601_stack$X199601_20),na.rm=T),mean(values(ex_array199602_stack$X199602_20),na.rm=T),mean(values(ex_array199603_stack$X199603_20),na.rm=T),
                 mean(values(ex_array199604_stack$X199604_20),na.rm=T),mean(values(ex_array199605_stack$X199605_20),na.rm=T),mean(values(ex_array199606_stack$X199606_20),na.rm=T),
                 mean(values(ex_array199607_stack$X199607_20),na.rm=T),mean(values(ex_array199608_stack$X199608_20),na.rm=T),mean(values(ex_array199609_stack$X199609_20),na.rm=T),
                 mean(values(ex_array199610_stack$X199610_20),na.rm=T),mean(values(ex_array199611_stack$X199611_20),na.rm=T),mean(values(ex_array199612_stack$X199612_20),na.rm=T))
t1996_temps40<-c(mean(values(ex_array199601_stack$X199601_40),na.rm=T),mean(values(ex_array199602_stack$X199602_40),na.rm=T),mean(values(ex_array199603_stack$X199603_40),na.rm=T),
                 mean(values(ex_array199604_stack$X199604_40),na.rm=T),mean(values(ex_array199605_stack$X199605_40),na.rm=T),mean(values(ex_array199606_stack$X199606_40),na.rm=T),
                 mean(values(ex_array199607_stack$X199607_40),na.rm=T),mean(values(ex_array199608_stack$X199608_40),na.rm=T),mean(values(ex_array199609_stack$X199609_40),na.rm=T),
                 mean(values(ex_array199610_stack$X199610_40),na.rm=T),mean(values(ex_array199611_stack$X199611_40),na.rm=T),mean(values(ex_array199612_stack$X199612_40),na.rm=T))
t1996_temps60<-c(mean(values(ex_array199601_stack$X199601_60),na.rm=T),mean(values(ex_array199602_stack$X199602_60),na.rm=T),mean(values(ex_array199603_stack$X199603_60),na.rm=T),
                 mean(values(ex_array199604_stack$X199604_60),na.rm=T),mean(values(ex_array199605_stack$X199605_60),na.rm=T),mean(values(ex_array199606_stack$X199606_60),na.rm=T),
                 mean(values(ex_array199607_stack$X199607_60),na.rm=T),mean(values(ex_array199608_stack$X199608_60),na.rm=T),mean(values(ex_array199609_stack$X199609_60),na.rm=T),
                 mean(values(ex_array199610_stack$X199610_60),na.rm=T),mean(values(ex_array199611_stack$X199611_60),na.rm=T),mean(values(ex_array199612_stack$X199612_60),na.rm=T))
t1996_temps80<-c(mean(values(ex_array199601_stack$X199601_80),na.rm=T),mean(values(ex_array199602_stack$X199602_80),na.rm=T),mean(values(ex_array199603_stack$X199603_80),na.rm=T),
                 mean(values(ex_array199604_stack$X199604_80),na.rm=T),mean(values(ex_array199605_stack$X199605_80),na.rm=T),mean(values(ex_array199606_stack$X199606_80),na.rm=T),
                 mean(values(ex_array199607_stack$X199607_80),na.rm=T),mean(values(ex_array199608_stack$X199608_80),na.rm=T),mean(values(ex_array199609_stack$X199609_80),na.rm=T),
                 mean(values(ex_array199610_stack$X199610_80),na.rm=T),mean(values(ex_array199611_stack$X199611_80),na.rm=T),mean(values(ex_array199612_stack$X199612_80),na.rm=T))
t1996_temps100<-c(mean(values(ex_array199601_stack$X199601_100),na.rm=T),mean(values(ex_array199602_stack$X199602_100),na.rm=T),mean(values(ex_array199603_stack$X199603_100),na.rm=T),
                  mean(values(ex_array199604_stack$X199604_100),na.rm=T),mean(values(ex_array199605_stack$X199605_100),na.rm=T),mean(values(ex_array199606_stack$X199606_100),na.rm=T),
                  mean(values(ex_array199607_stack$X199607_100),na.rm=T),mean(values(ex_array199608_stack$X199608_100),na.rm=T),mean(values(ex_array199609_stack$X199609_100),na.rm=T),
                  mean(values(ex_array199610_stack$X199610_100),na.rm=T),mean(values(ex_array199611_stack$X199611_100),na.rm=T),mean(values(ex_array199612_stack$X199612_100),na.rm=T))
t1996_temps250<-c(mean(values(ex_array199601_stack$X199601_250),na.rm=T),mean(values(ex_array199602_stack$X199602_250),na.rm=T),mean(values(ex_array199603_stack$X199603_250),na.rm=T),
                  mean(values(ex_array199604_stack$X199604_250),na.rm=T),mean(values(ex_array199605_stack$X199605_250),na.rm=T),mean(values(ex_array199606_stack$X199606_250),na.rm=T),
                  mean(values(ex_array199607_stack$X199607_250),na.rm=T),mean(values(ex_array199608_stack$X199608_250),na.rm=T),mean(values(ex_array199609_stack$X199609_250),na.rm=T),
                  mean(values(ex_array199610_stack$X199610_250),na.rm=T),mean(values(ex_array199611_stack$X199611_250),na.rm=T),mean(values(ex_array199612_stack$X199612_250),na.rm=T))

#future years
#0-20
future0_20_temps20<-c(mean(values(Jan_0_20_future_mask$layer.1),na.rm=T),mean(values(Feb_0_20_future_mask$layer.1),na.rm=T),mean(values(Mar_0_20_future_mask$layer.1),na.rm=T),
                      mean(values(Apr_0_20_future_mask$layer.1),na.rm=T),mean(values(May_0_20_future_mask$layer.1),na.rm=T),mean(values(Jun_0_20_future_mask$layer.1),na.rm=T),
                      mean(values(Jul_0_20_future_mask$layer.1),na.rm=T),mean(values(Aug_0_20_future_mask$layer.1),na.rm=T),mean(values(Sept_0_20_future_mask$layer.1),na.rm=T),
                      mean(values(Oct_0_20_future_mask$layer.1),na.rm=T),mean(values(Nov_0_20_future_mask$layer.1),na.rm=T),mean(values(Dec_0_20_future_mask$layer.1),na.rm=T))
future0_20_temps40<-c(mean(values(Jan_0_20_future_mask$layer.2),na.rm=T),mean(values(Feb_0_20_future_mask$layer.2),na.rm=T),mean(values(Mar_0_20_future_mask$layer.2),na.rm=T),
                      mean(values(Apr_0_20_future_mask$layer.2),na.rm=T),mean(values(May_0_20_future_mask$layer.2),na.rm=T),mean(values(Jun_0_20_future_mask$layer.2),na.rm=T),
                      mean(values(Jul_0_20_future_mask$layer.2),na.rm=T),mean(values(Aug_0_20_future_mask$layer.2),na.rm=T),mean(values(Sept_0_20_future_mask$layer.2),na.rm=T),
                      mean(values(Oct_0_20_future_mask$layer.2),na.rm=T),mean(values(Nov_0_20_future_mask$layer.2),na.rm=T),mean(values(Dec_0_20_future_mask$layer.2),na.rm=T))
future0_20_temps60<-c(mean(values(Jan_0_20_future_mask$layer.3),na.rm=T),mean(values(Feb_0_20_future_mask$layer.3),na.rm=T),mean(values(Mar_0_20_future_mask$layer.3),na.rm=T),
                      mean(values(Apr_0_20_future_mask$layer.3),na.rm=T),mean(values(May_0_20_future_mask$layer.3),na.rm=T),mean(values(Jun_0_20_future_mask$layer.3),na.rm=T),
                      mean(values(Jul_0_20_future_mask$layer.3),na.rm=T),mean(values(Aug_0_20_future_mask$layer.3),na.rm=T),mean(values(Sept_0_20_future_mask$layer.3),na.rm=T),
                      mean(values(Oct_0_20_future_mask$layer.3),na.rm=T),mean(values(Nov_0_20_future_mask$layer.3),na.rm=T),mean(values(Dec_0_20_future_mask$layer.3),na.rm=T))
future0_20_temps80<-c(mean(values(Jan_0_20_future_mask$layer.4),na.rm=T),mean(values(Feb_0_20_future_mask$layer.4),na.rm=T),mean(values(Mar_0_20_future_mask$layer.4),na.rm=T),
                      mean(values(Apr_0_20_future_mask$layer.4),na.rm=T),mean(values(May_0_20_future_mask$layer.4),na.rm=T),mean(values(Jun_0_20_future_mask$layer.4),na.rm=T),
                      mean(values(Jul_0_20_future_mask$layer.4),na.rm=T),mean(values(Aug_0_20_future_mask$layer.4),na.rm=T),mean(values(Sept_0_20_future_mask$layer.4),na.rm=T),
                      mean(values(Oct_0_20_future_mask$layer.4),na.rm=T),mean(values(Nov_0_20_future_mask$layer.4),na.rm=T),mean(values(Dec_0_20_future_mask$layer.4),na.rm=T))
future0_20_temps100<-c(mean(values(Jan_0_20_future_mask$layer.5),na.rm=T),mean(values(Feb_0_20_future_mask$layer.5),na.rm=T),mean(values(Mar_0_20_future_mask$layer.5),na.rm=T),
                       mean(values(Apr_0_20_future_mask$layer.5),na.rm=T),mean(values(May_0_20_future_mask$layer.5),na.rm=T),mean(values(Jun_0_20_future_mask$layer.5),na.rm=T),
                       mean(values(Jul_0_20_future_mask$layer.5),na.rm=T),mean(values(Aug_0_20_future_mask$layer.5),na.rm=T),mean(values(Sept_0_20_future_mask$layer.5),na.rm=T),
                       mean(values(Oct_0_20_future_mask$layer.5),na.rm=T),mean(values(Nov_0_20_future_mask$layer.5),na.rm=T),mean(values(Dec_0_20_future_mask$layer.5),na.rm=T))
future0_20_temps250<-c(mean(values(Jan_0_20_future_mask$layer.6),na.rm=T),mean(values(Feb_0_20_future_mask$layer.6),na.rm=T),mean(values(Mar_0_20_future_mask$layer.6),na.rm=T),
                       mean(values(Apr_0_20_future_mask$layer.6),na.rm=T),mean(values(May_0_20_future_mask$layer.6),na.rm=T),mean(values(Jun_0_20_future_mask$layer.6),na.rm=T),
                       mean(values(Jul_0_20_future_mask$layer.6),na.rm=T),mean(values(Aug_0_20_future_mask$layer.6),na.rm=T),mean(values(Sept_0_20_future_mask$layer.6),na.rm=T),
                       mean(values(Oct_0_20_future_mask$layer.6),na.rm=T),mean(values(Nov_0_20_future_mask$layer.6),na.rm=T),mean(values(Dec_0_20_future_mask$layer.6),na.rm=T))
#20-40
future20_40_temps20<-c(mean(values(Jan_20_40_future_mask$layer.1),na.rm=T),mean(values(Feb_20_40_future_mask$layer.1),na.rm=T),mean(values(Mar_20_40_future_mask$layer.1),na.rm=T),
                       mean(values(Apr_20_40_future_mask$layer.1),na.rm=T),mean(values(May_20_40_future_mask$layer.1),na.rm=T),mean(values(Jun_20_40_future_mask$layer.1),na.rm=T),
                       mean(values(Jul_20_40_future_mask$layer.1),na.rm=T),mean(values(Aug_20_40_future_mask$layer.1),na.rm=T),mean(values(Sept_20_40_future_mask$layer.1),na.rm=T),
                       mean(values(Oct_20_40_future_mask$layer.1),na.rm=T),mean(values(Nov_20_40_future_mask$layer.1),na.rm=T),mean(values(Dec_20_40_future_mask$layer.1),na.rm=T))
future20_40_temps40<-c(mean(values(Jan_20_40_future_mask$layer.2),na.rm=T),mean(values(Feb_20_40_future_mask$layer.2),na.rm=T),mean(values(Mar_20_40_future_mask$layer.2),na.rm=T),
                       mean(values(Apr_20_40_future_mask$layer.2),na.rm=T),mean(values(May_20_40_future_mask$layer.2),na.rm=T),mean(values(Jun_20_40_future_mask$layer.2),na.rm=T),
                       mean(values(Jul_20_40_future_mask$layer.2),na.rm=T),mean(values(Aug_20_40_future_mask$layer.2),na.rm=T),mean(values(Sept_20_40_future_mask$layer.2),na.rm=T),
                       mean(values(Oct_20_40_future_mask$layer.2),na.rm=T),mean(values(Nov_20_40_future_mask$layer.2),na.rm=T),mean(values(Dec_20_40_future_mask$layer.2),na.rm=T))
future20_40_temps60<-c(mean(values(Jan_20_40_future_mask$layer.3),na.rm=T),mean(values(Feb_20_40_future_mask$layer.3),na.rm=T),mean(values(Mar_20_40_future_mask$layer.3),na.rm=T),
                       mean(values(Apr_20_40_future_mask$layer.3),na.rm=T),mean(values(May_20_40_future_mask$layer.3),na.rm=T),mean(values(Jun_20_40_future_mask$layer.3),na.rm=T),
                       mean(values(Jul_20_40_future_mask$layer.3),na.rm=T),mean(values(Aug_20_40_future_mask$layer.3),na.rm=T),mean(values(Sept_20_40_future_mask$layer.3),na.rm=T),
                       mean(values(Oct_20_40_future_mask$layer.3),na.rm=T),mean(values(Nov_20_40_future_mask$layer.3),na.rm=T),mean(values(Dec_20_40_future_mask$layer.3),na.rm=T))
future20_40_temps80<-c(mean(values(Jan_20_40_future_mask$layer.4),na.rm=T),mean(values(Feb_20_40_future_mask$layer.4),na.rm=T),mean(values(Mar_20_40_future_mask$layer.4),na.rm=T),
                       mean(values(Apr_20_40_future_mask$layer.4),na.rm=T),mean(values(May_20_40_future_mask$layer.4),na.rm=T),mean(values(Jun_20_40_future_mask$layer.4),na.rm=T),
                       mean(values(Jul_20_40_future_mask$layer.4),na.rm=T),mean(values(Aug_20_40_future_mask$layer.4),na.rm=T),mean(values(Sept_20_40_future_mask$layer.4),na.rm=T),
                       mean(values(Oct_20_40_future_mask$layer.4),na.rm=T),mean(values(Nov_20_40_future_mask$layer.4),na.rm=T),mean(values(Dec_20_40_future_mask$layer.4),na.rm=T))
future20_40_temps100<-c(mean(values(Jan_20_40_future_mask$layer.5),na.rm=T),mean(values(Feb_20_40_future_mask$layer.5),na.rm=T),mean(values(Mar_20_40_future_mask$layer.5),na.rm=T),
                        mean(values(Apr_20_40_future_mask$layer.5),na.rm=T),mean(values(May_20_40_future_mask$layer.5),na.rm=T),mean(values(Jun_20_40_future_mask$layer.5),na.rm=T),
                        mean(values(Jul_20_40_future_mask$layer.5),na.rm=T),mean(values(Aug_20_40_future_mask$layer.5),na.rm=T),mean(values(Sept_20_40_future_mask$layer.5),na.rm=T),
                        mean(values(Oct_20_40_future_mask$layer.5),na.rm=T),mean(values(Nov_20_40_future_mask$layer.5),na.rm=T),mean(values(Dec_20_40_future_mask$layer.5),na.rm=T))
future20_40_temps250<-c(mean(values(Jan_20_40_future_mask$layer.6),na.rm=T),mean(values(Feb_20_40_future_mask$layer.6),na.rm=T),mean(values(Mar_20_40_future_mask$layer.6),na.rm=T),
                        mean(values(Apr_20_40_future_mask$layer.6),na.rm=T),mean(values(May_20_40_future_mask$layer.6),na.rm=T),mean(values(Jun_20_40_future_mask$layer.6),na.rm=T),
                        mean(values(Jul_20_40_future_mask$layer.6),na.rm=T),mean(values(Aug_20_40_future_mask$layer.6),na.rm=T),mean(values(Sept_20_40_future_mask$layer.6),na.rm=T),
                        mean(values(Oct_20_40_future_mask$layer.6),na.rm=T),mean(values(Nov_20_40_future_mask$layer.6),na.rm=T),mean(values(Dec_20_40_future_mask$layer.6),na.rm=T))
#40-60
future40_60_temps20<-c(mean(values(Jan_40_60_future_mask$layer.1),na.rm=T),mean(values(Feb_40_60_future_mask$layer.1),na.rm=T),mean(values(Mar_40_60_future_mask$layer.1),na.rm=T),
                       mean(values(Apr_40_60_future_mask$layer.1),na.rm=T),mean(values(May_40_60_future_mask$layer.1),na.rm=T),mean(values(Jun_40_60_future_mask$layer.1),na.rm=T),
                       mean(values(Jul_40_60_future_mask$layer.1),na.rm=T),mean(values(Aug_40_60_future_mask$layer.1),na.rm=T),mean(values(Sept_40_60_future_mask$layer.1),na.rm=T),
                       mean(values(Oct_40_60_future_mask$layer.1),na.rm=T),mean(values(Nov_40_60_future_mask$layer.1),na.rm=T),mean(values(Dec_40_60_future_mask$layer.1),na.rm=T))
future40_60_temps40<-c(mean(values(Jan_40_60_future_mask$layer.2),na.rm=T),mean(values(Feb_40_60_future_mask$layer.2),na.rm=T),mean(values(Mar_40_60_future_mask$layer.2),na.rm=T),
                       mean(values(Apr_40_60_future_mask$layer.2),na.rm=T),mean(values(May_40_60_future_mask$layer.2),na.rm=T),mean(values(Jun_40_60_future_mask$layer.2),na.rm=T),
                       mean(values(Jul_40_60_future_mask$layer.2),na.rm=T),mean(values(Aug_40_60_future_mask$layer.2),na.rm=T),mean(values(Sept_40_60_future_mask$layer.2),na.rm=T),
                       mean(values(Oct_40_60_future_mask$layer.2),na.rm=T),mean(values(Nov_40_60_future_mask$layer.2),na.rm=T),mean(values(Dec_40_60_future_mask$layer.2),na.rm=T))
future40_60_temps60<-c(mean(values(Jan_40_60_future_mask$layer.3),na.rm=T),mean(values(Feb_40_60_future_mask$layer.3),na.rm=T),mean(values(Mar_40_60_future_mask$layer.3),na.rm=T),
                       mean(values(Apr_40_60_future_mask$layer.3),na.rm=T),mean(values(May_40_60_future_mask$layer.3),na.rm=T),mean(values(Jun_40_60_future_mask$layer.3),na.rm=T),
                       mean(values(Jul_40_60_future_mask$layer.3),na.rm=T),mean(values(Aug_40_60_future_mask$layer.3),na.rm=T),mean(values(Sept_40_60_future_mask$layer.3),na.rm=T),
                       mean(values(Oct_40_60_future_mask$layer.3),na.rm=T),mean(values(Nov_40_60_future_mask$layer.3),na.rm=T),mean(values(Dec_40_60_future_mask$layer.3),na.rm=T))
future40_60_temps80<-c(mean(values(Jan_40_60_future_mask$layer.4),na.rm=T),mean(values(Feb_40_60_future_mask$layer.4),na.rm=T),mean(values(Mar_40_60_future_mask$layer.4),na.rm=T),
                       mean(values(Apr_40_60_future_mask$layer.4),na.rm=T),mean(values(May_40_60_future_mask$layer.4),na.rm=T),mean(values(Jun_40_60_future_mask$layer.4),na.rm=T),
                       mean(values(Jul_40_60_future_mask$layer.4),na.rm=T),mean(values(Aug_40_60_future_mask$layer.4),na.rm=T),mean(values(Sept_40_60_future_mask$layer.4),na.rm=T),
                       mean(values(Oct_40_60_future_mask$layer.4),na.rm=T),mean(values(Nov_40_60_future_mask$layer.4),na.rm=T),mean(values(Dec_40_60_future_mask$layer.4),na.rm=T))
future40_60_temps100<-c(mean(values(Jan_40_60_future_mask$layer.5),na.rm=T),mean(values(Feb_40_60_future_mask$layer.5),na.rm=T),mean(values(Mar_40_60_future_mask$layer.5),na.rm=T),
                        mean(values(Apr_40_60_future_mask$layer.5),na.rm=T),mean(values(May_40_60_future_mask$layer.5),na.rm=T),mean(values(Jun_40_60_future_mask$layer.5),na.rm=T),
                        mean(values(Jul_40_60_future_mask$layer.5),na.rm=T),mean(values(Aug_40_60_future_mask$layer.5),na.rm=T),mean(values(Sept_40_60_future_mask$layer.5),na.rm=T),
                        mean(values(Oct_40_60_future_mask$layer.5),na.rm=T),mean(values(Nov_40_60_future_mask$layer.5),na.rm=T),mean(values(Dec_40_60_future_mask$layer.5),na.rm=T))
future40_60_temps250<-c(mean(values(Jan_40_60_future_mask$layer.6),na.rm=T),mean(values(Feb_40_60_future_mask$layer.6),na.rm=T),mean(values(Mar_40_60_future_mask$layer.6),na.rm=T),
                        mean(values(Apr_40_60_future_mask$layer.6),na.rm=T),mean(values(May_40_60_future_mask$layer.6),na.rm=T),mean(values(Jun_40_60_future_mask$layer.6),na.rm=T),
                        mean(values(Jul_40_60_future_mask$layer.6),na.rm=T),mean(values(Aug_40_60_future_mask$layer.6),na.rm=T),mean(values(Sept_40_60_future_mask$layer.6),na.rm=T),
                        mean(values(Oct_40_60_future_mask$layer.6),na.rm=T),mean(values(Nov_40_60_future_mask$layer.6),na.rm=T),mean(values(Dec_40_60_future_mask$layer.6),na.rm=T))
#60-80
future60_80_temps20<-c(mean(values(Jan_60_80_future_mask$layer.1),na.rm=T),mean(values(Feb_60_80_future_mask$layer.1),na.rm=T),mean(values(Mar_60_80_future_mask$layer.1),na.rm=T),
                       mean(values(Apr_60_80_future_mask$layer.1),na.rm=T),mean(values(May_60_80_future_mask$layer.1),na.rm=T),mean(values(Jun_60_80_future_mask$layer.1),na.rm=T),
                       mean(values(Jul_60_80_future_mask$layer.1),na.rm=T),mean(values(Aug_60_80_future_mask$layer.1),na.rm=T),mean(values(Sept_60_80_future_mask$layer.1),na.rm=T),
                       mean(values(Oct_60_80_future_mask$layer.1),na.rm=T),mean(values(Nov_60_80_future_mask$layer.1),na.rm=T),mean(values(Dec_60_80_future_mask$layer.1),na.rm=T))
future60_80_temps40<-c(mean(values(Jan_60_80_future_mask$layer.2),na.rm=T),mean(values(Feb_60_80_future_mask$layer.2),na.rm=T),mean(values(Mar_60_80_future_mask$layer.2),na.rm=T),
                       mean(values(Apr_60_80_future_mask$layer.2),na.rm=T),mean(values(May_60_80_future_mask$layer.2),na.rm=T),mean(values(Jun_60_80_future_mask$layer.2),na.rm=T),
                       mean(values(Jul_60_80_future_mask$layer.2),na.rm=T),mean(values(Aug_60_80_future_mask$layer.2),na.rm=T),mean(values(Sept_60_80_future_mask$layer.2),na.rm=T),
                       mean(values(Oct_60_80_future_mask$layer.2),na.rm=T),mean(values(Nov_60_80_future_mask$layer.2),na.rm=T),mean(values(Dec_60_80_future_mask$layer.2),na.rm=T))
future60_80_temps60<-c(mean(values(Jan_60_80_future_mask$layer.3),na.rm=T),mean(values(Feb_60_80_future_mask$layer.3),na.rm=T),mean(values(Mar_60_80_future_mask$layer.3),na.rm=T),
                       mean(values(Apr_60_80_future_mask$layer.3),na.rm=T),mean(values(May_60_80_future_mask$layer.3),na.rm=T),mean(values(Jun_60_80_future_mask$layer.3),na.rm=T),
                       mean(values(Jul_60_80_future_mask$layer.3),na.rm=T),mean(values(Aug_60_80_future_mask$layer.3),na.rm=T),mean(values(Sept_60_80_future_mask$layer.3),na.rm=T),
                       mean(values(Oct_60_80_future_mask$layer.3),na.rm=T),mean(values(Nov_60_80_future_mask$layer.3),na.rm=T),mean(values(Dec_60_80_future_mask$layer.3),na.rm=T))
future60_80_temps80<-c(mean(values(Jan_60_80_future_mask$layer.4),na.rm=T),mean(values(Feb_60_80_future_mask$layer.4),na.rm=T),mean(values(Mar_60_80_future_mask$layer.4),na.rm=T),
                       mean(values(Apr_60_80_future_mask$layer.4),na.rm=T),mean(values(May_60_80_future_mask$layer.4),na.rm=T),mean(values(Jun_60_80_future_mask$layer.4),na.rm=T),
                       mean(values(Jul_60_80_future_mask$layer.4),na.rm=T),mean(values(Aug_60_80_future_mask$layer.4),na.rm=T),mean(values(Sept_60_80_future_mask$layer.4),na.rm=T),
                       mean(values(Oct_60_80_future_mask$layer.4),na.rm=T),mean(values(Nov_60_80_future_mask$layer.4),na.rm=T),mean(values(Dec_60_80_future_mask$layer.4),na.rm=T))
future60_80_temps100<-c(mean(values(Jan_60_80_future_mask$layer.5),na.rm=T),mean(values(Feb_60_80_future_mask$layer.5),na.rm=T),mean(values(Mar_60_80_future_mask$layer.5),na.rm=T),
                        mean(values(Apr_60_80_future_mask$layer.5),na.rm=T),mean(values(May_60_80_future_mask$layer.5),na.rm=T),mean(values(Jun_60_80_future_mask$layer.5),na.rm=T),
                        mean(values(Jul_60_80_future_mask$layer.5),na.rm=T),mean(values(Aug_60_80_future_mask$layer.5),na.rm=T),mean(values(Sept_60_80_future_mask$layer.5),na.rm=T),
                        mean(values(Oct_60_80_future_mask$layer.5),na.rm=T),mean(values(Nov_60_80_future_mask$layer.5),na.rm=T),mean(values(Dec_60_80_future_mask$layer.5),na.rm=T))
future60_80_temps250<-c(mean(values(Jan_60_80_future_mask$layer.6),na.rm=T),mean(values(Feb_60_80_future_mask$layer.6),na.rm=T),mean(values(Mar_60_80_future_mask$layer.6),na.rm=T),
                        mean(values(Apr_60_80_future_mask$layer.6),na.rm=T),mean(values(May_60_80_future_mask$layer.6),na.rm=T),mean(values(Jun_60_80_future_mask$layer.6),na.rm=T),
                        mean(values(Jul_60_80_future_mask$layer.6),na.rm=T),mean(values(Aug_60_80_future_mask$layer.6),na.rm=T),mean(values(Sept_60_80_future_mask$layer.6),na.rm=T),
                        mean(values(Oct_60_80_future_mask$layer.6),na.rm=T),mean(values(Nov_60_80_future_mask$layer.6),na.rm=T),mean(values(Dec_60_80_future_mask$layer.6),na.rm=T))

#Future Extreme Year (70)
future70_temps20<-c(mean(values(Jan_70_future_mask$layer.70.1),na.rm=T),mean(values(Feb_70_future_mask$layer.70.1),na.rm=T),mean(values(Mar_70_future_mask$layer.70.1),na.rm=T),
                    mean(values(Apr_70_future_mask$layer.70.1),na.rm=T),mean(values(May_70_future_mask$layer.70.1),na.rm=T),mean(values(Jun_70_future_mask$layer.70.1),na.rm=T),
                    mean(values(Jul_70_future_mask$layer.70.1),na.rm=T),mean(values(Aug_70_future_mask$layer.70.1),na.rm=T),mean(values(Sept_70_future_mask$layer.70.1),na.rm=T),
                    mean(values(Oct_70_future_mask$layer.70.1),na.rm=T),mean(values(Nov_70_future_mask$layer.70.1),na.rm=T),mean(values(Dec_70_future_mask$layer.70.1),na.rm=T))
future70_temps40<-c(mean(values(Jan_70_future_mask$layer.70.2),na.rm=T),mean(values(Feb_70_future_mask$layer.70.2),na.rm=T),mean(values(Mar_70_future_mask$layer.70.2),na.rm=T),
                    mean(values(Apr_70_future_mask$layer.70.2),na.rm=T),mean(values(May_70_future_mask$layer.70.2),na.rm=T),mean(values(Jun_70_future_mask$layer.70.2),na.rm=T),
                    mean(values(Jul_70_future_mask$layer.70.2),na.rm=T),mean(values(Aug_70_future_mask$layer.70.2),na.rm=T),mean(values(Sept_70_future_mask$layer.70.2),na.rm=T),
                    mean(values(Oct_70_future_mask$layer.70.2),na.rm=T),mean(values(Nov_70_future_mask$layer.70.2),na.rm=T),mean(values(Dec_70_future_mask$layer.70.2),na.rm=T))
future70_temps60<-c(mean(values(Jan_70_future_mask$layer.70.3),na.rm=T),mean(values(Feb_70_future_mask$layer.70.3),na.rm=T),mean(values(Mar_70_future_mask$layer.70.3),na.rm=T),
                    mean(values(Apr_70_future_mask$layer.70.3),na.rm=T),mean(values(May_70_future_mask$layer.70.3),na.rm=T),mean(values(Jun_70_future_mask$layer.70.3),na.rm=T),
                    mean(values(Jul_70_future_mask$layer.70.3),na.rm=T),mean(values(Aug_70_future_mask$layer.70.3),na.rm=T),mean(values(Sept_70_future_mask$layer.70.3),na.rm=T),
                    mean(values(Oct_70_future_mask$layer.70.3),na.rm=T),mean(values(Nov_70_future_mask$layer.70.3),na.rm=T),mean(values(Dec_70_future_mask$layer.70.3),na.rm=T))
future70_temps80<-c(mean(values(Jan_70_future_mask$layer.70.4),na.rm=T),mean(values(Feb_70_future_mask$layer.70.4),na.rm=T),mean(values(Mar_70_future_mask$layer.70.4),na.rm=T),
                    mean(values(Apr_70_future_mask$layer.70.4),na.rm=T),mean(values(May_70_future_mask$layer.70.4),na.rm=T),mean(values(Jun_70_future_mask$layer.70.4),na.rm=T),
                    mean(values(Jul_70_future_mask$layer.70.4),na.rm=T),mean(values(Aug_70_future_mask$layer.70.4),na.rm=T),mean(values(Sept_70_future_mask$layer.70.4),na.rm=T),
                    mean(values(Oct_70_future_mask$layer.70.4),na.rm=T),mean(values(Nov_70_future_mask$layer.70.4),na.rm=T),mean(values(Dec_70_future_mask$layer.70.4),na.rm=T))
future70_temps100<-c(mean(values(Jan_70_future_mask$layer.70.5),na.rm=T),mean(values(Feb_70_future_mask$layer.70.5),na.rm=T),mean(values(Mar_70_future_mask$layer.70.5),na.rm=T),
                     mean(values(Apr_70_future_mask$layer.70.5),na.rm=T),mean(values(May_70_future_mask$layer.70.5),na.rm=T),mean(values(Jun_70_future_mask$layer.70.5),na.rm=T),
                     mean(values(Jul_70_future_mask$layer.70.5),na.rm=T),mean(values(Aug_70_future_mask$layer.70.5),na.rm=T),mean(values(Sept_70_future_mask$layer.70.5),na.rm=T),
                     mean(values(Oct_70_future_mask$layer.70.5),na.rm=T),mean(values(Nov_70_future_mask$layer.70.5),na.rm=T),mean(values(Dec_70_future_mask$layer.70.5),na.rm=T))
future70_temps250<-c(mean(values(Jan_70_future_mask$layer.70.6),na.rm=T),mean(values(Feb_70_future_mask$layer.70.6),na.rm=T),mean(values(Mar_70_future_mask$layer.70.6),na.rm=T),
                     mean(values(Apr_70_future_mask$layer.70.6),na.rm=T),mean(values(May_70_future_mask$layer.70.6),na.rm=T),mean(values(Jun_70_future_mask$layer.70.6),na.rm=T),
                     mean(values(Jul_70_future_mask$layer.70.6),na.rm=T),mean(values(Aug_70_future_mask$layer.70.6),na.rm=T),mean(values(Sept_70_future_mask$layer.70.6),na.rm=T),
                     mean(values(Oct_70_future_mask$layer.70.6),na.rm=T),mean(values(Nov_70_future_mask$layer.70.6),na.rm=T),mean(values(Dec_70_future_mask$layer.70.6),na.rm=T))

#Histograms of available temps
#Jan
par(mfrow=c(2,4))
hist201201<-hist(c(values(ex_array201201_stack$X201201_20),values(ex_array201201_stack$X201201_40),values(ex_array201201_stack$X201201_60),
                   values(ex_array201201_stack$X201201_80),values(ex_array201201_stack$X201201_100),values(ex_array201201_stack$X201201_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_01")
hist199601<-hist(c(values(ex_array199601_stack$X199601_20),values(ex_array199601_stack$X199601_40),values(ex_array199601_stack$X199601_60),
                   values(ex_array199601_stack$X199601_80),values(ex_array199601_stack$X199601_100),values(ex_array199601_stack$X199601_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_01")
histclim01<-hist(c(values(Jan_stack$Jan_20),values(Jan_stack$Jan_40),values(Jan_stack$Jan_60),
                   values(Jan_stack$Jan_80),values(Jan_stack$Jan_100),values(Jan_stack$Jan_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_01")
hist0_20f01<-hist(c(values(Jan_0_20_future_mask$layer.1),values(Jan_0_20_future_mask$layer.2),values(Jan_0_20_future_mask$layer.3),
                    values(Jan_0_20_future_mask$layer.4),values(Jan_0_20_future_mask$layer.5),values(Jan_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_01")
hist20_40f01<-hist(c(values(Jan_20_40_future_mask$layer.1),values(Jan_20_40_future_mask$layer.2),values(Jan_20_40_future_mask$layer.3),
                     values(Jan_20_40_future_mask$layer.4),values(Jan_20_40_future_mask$layer.5),values(Jan_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_01")
hist40_60f01<-hist(c(values(Jan_40_60_future_mask$layer.1),values(Jan_40_60_future_mask$layer.2),values(Jan_40_60_future_mask$layer.3),
                     values(Jan_40_60_future_mask$layer.4),values(Jan_40_60_future_mask$layer.5),values(Jan_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_01")
hist60_80f01<-hist(c(values(Jan_60_80_future_mask$layer.1),values(Jan_60_80_future_mask$layer.2),values(Jan_60_80_future_mask$layer.3),
                     values(Jan_60_80_future_mask$layer.4),values(Jan_60_80_future_mask$layer.5),values(Jan_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_01")
hist70f01<-hist(c(values(Jan_70_future_mask$layer.70.1),values(Jan_70_future_mask$layer.70.2),values(Jan_70_future_mask$layer.70.3),
                  values(Jan_70_future_mask$layer.70.4),values(Jan_70_future_mask$layer.70.5),values(Jan_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_01")
#Feb
par(mfrow=c(2,4))
hist201202<-hist(c(values(ex_array201202_stack$X201202_20),values(ex_array201202_stack$X201202_40),values(ex_array201202_stack$X201202_60),
                   values(ex_array201202_stack$X201202_80),values(ex_array201202_stack$X201202_100),values(ex_array201202_stack$X201202_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_02")
hist199602<-hist(c(values(ex_array199602_stack$X199602_20),values(ex_array199602_stack$X199602_40),values(ex_array199602_stack$X199602_60),
                   values(ex_array199602_stack$X199602_80),values(ex_array199602_stack$X199602_100),values(ex_array199602_stack$X199602_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_02")
histclim02<-hist(c(values(Feb_stack$Feb_20),values(Feb_stack$Feb_40),values(Feb_stack$Feb_60),
                   values(Feb_stack$Feb_80),values(Feb_stack$Feb_100),values(Feb_stack$Feb_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_02")
hist0_20f02<-hist(c(values(Feb_0_20_future_mask$layer.1),values(Feb_0_20_future_mask$layer.2),values(Feb_0_20_future_mask$layer.3),
                    values(Feb_0_20_future_mask$layer.4),values(Feb_0_20_future_mask$layer.5),values(Feb_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_02")
hist20_40f02<-hist(c(values(Feb_20_40_future_mask$layer.1),values(Feb_20_40_future_mask$layer.2),values(Feb_20_40_future_mask$layer.3),
                     values(Feb_20_40_future_mask$layer.4),values(Feb_20_40_future_mask$layer.5),values(Feb_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_02")
hist40_60f02<-hist(c(values(Feb_40_60_future_mask$layer.1),values(Feb_40_60_future_mask$layer.2),values(Feb_40_60_future_mask$layer.3),
                     values(Feb_40_60_future_mask$layer.4),values(Feb_40_60_future_mask$layer.5),values(Feb_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_02")
hist60_80f02<-hist(c(values(Feb_60_80_future_mask$layer.1),values(Feb_60_80_future_mask$layer.2),values(Feb_60_80_future_mask$layer.3),
                     values(Feb_60_80_future_mask$layer.4),values(Feb_60_80_future_mask$layer.5),values(Feb_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_02")
hist70f02<-hist(c(values(Feb_70_future_mask$layer.70.1),values(Feb_70_future_mask$layer.70.2),values(Feb_70_future_mask$layer.70.3),
                  values(Feb_70_future_mask$layer.70.4),values(Feb_70_future_mask$layer.70.5),values(Feb_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_02")
#Mar
par(mfrow=c(2,4))
hist201203<-hist(c(values(ex_array201203_stack$X201203_20),values(ex_array201203_stack$X201203_40),values(ex_array201203_stack$X201203_60),
                   values(ex_array201203_stack$X201203_80),values(ex_array201203_stack$X201203_100),values(ex_array201203_stack$X201203_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_03")
hist199603<-hist(c(values(ex_array199603_stack$X199603_20),values(ex_array199603_stack$X199603_40),values(ex_array199603_stack$X199603_60),
                   values(ex_array199603_stack$X199603_80),values(ex_array199603_stack$X199603_100),values(ex_array199603_stack$X199603_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_03")
histclim03<-hist(c(values(Mar_stack$Mar_20),values(Mar_stack$Mar_40),values(Mar_stack$Mar_60),
                   values(Mar_stack$Mar_80),values(Mar_stack$Mar_100),values(Mar_stack$Mar_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_03")
hist0_20f03<-hist(c(values(Mar_0_20_future_mask$layer.1),values(Mar_0_20_future_mask$layer.2),values(Mar_0_20_future_mask$layer.3),
                    values(Mar_0_20_future_mask$layer.4),values(Mar_0_20_future_mask$layer.5),values(Mar_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_03")
hist20_40f03<-hist(c(values(Mar_20_40_future_mask$layer.1),values(Mar_20_40_future_mask$layer.2),values(Mar_20_40_future_mask$layer.3),
                     values(Mar_20_40_future_mask$layer.4),values(Mar_20_40_future_mask$layer.5),values(Mar_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_03")
hist40_60f03<-hist(c(values(Mar_40_60_future_mask$layer.1),values(Mar_40_60_future_mask$layer.2),values(Mar_40_60_future_mask$layer.3),
                     values(Mar_40_60_future_mask$layer.4),values(Mar_40_60_future_mask$layer.5),values(Mar_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_03")
hist60_80f03<-hist(c(values(Mar_60_80_future_mask$layer.1),values(Mar_60_80_future_mask$layer.2),values(Mar_60_80_future_mask$layer.3),
                     values(Mar_60_80_future_mask$layer.4),values(Mar_60_80_future_mask$layer.5),values(Mar_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_03")
hist70f03<-hist(c(values(Mar_70_future_mask$layer.70.1),values(Mar_70_future_mask$layer.70.2),values(Mar_70_future_mask$layer.70.3),
                  values(Mar_70_future_mask$layer.70.4),values(Mar_70_future_mask$layer.70.5),values(Mar_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_03")
#Apr
par(mfrow=c(2,4))
hist201204<-hist(c(values(ex_array201204_stack$X201204_20),values(ex_array201204_stack$X201204_40),values(ex_array201204_stack$X201204_60),
                   values(ex_array201204_stack$X201204_80),values(ex_array201204_stack$X201204_100),values(ex_array201204_stack$X201204_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_04")
hist199604<-hist(c(values(ex_array199604_stack$X199604_20),values(ex_array199604_stack$X199604_40),values(ex_array199604_stack$X199604_60),
                   values(ex_array199604_stack$X199604_80),values(ex_array199604_stack$X199604_100),values(ex_array199604_stack$X199604_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_04")
histclim04<-hist(c(values(Apr_stack$Apr_20),values(Apr_stack$Apr_40),values(Apr_stack$Apr_60),
                   values(Apr_stack$Apr_80),values(Apr_stack$Apr_100),values(Apr_stack$Apr_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_04")
hist0_20f04<-hist(c(values(Apr_0_20_future_mask$layer.1),values(Apr_0_20_future_mask$layer.2),values(Apr_0_20_future_mask$layer.3),
                    values(Apr_0_20_future_mask$layer.4),values(Apr_0_20_future_mask$layer.5),values(Apr_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_04")
hist20_40f04<-hist(c(values(Apr_20_40_future_mask$layer.1),values(Apr_20_40_future_mask$layer.2),values(Apr_20_40_future_mask$layer.3),
                     values(Apr_20_40_future_mask$layer.4),values(Apr_20_40_future_mask$layer.5),values(Apr_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_04")
hist40_60f04<-hist(c(values(Apr_40_60_future_mask$layer.1),values(Apr_40_60_future_mask$layer.2),values(Apr_40_60_future_mask$layer.3),
                     values(Apr_40_60_future_mask$layer.4),values(Apr_40_60_future_mask$layer.5),values(Apr_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_04")
hist60_80f04<-hist(c(values(Apr_60_80_future_mask$layer.1),values(Apr_60_80_future_mask$layer.2),values(Apr_60_80_future_mask$layer.3),
                     values(Apr_60_80_future_mask$layer.4),values(Apr_60_80_future_mask$layer.5),values(Apr_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_04")
hist70f04<-hist(c(values(Apr_70_future_mask$layer.70.1),values(Apr_70_future_mask$layer.70.2),values(Apr_70_future_mask$layer.70.3),
                  values(Apr_70_future_mask$layer.70.4),values(Apr_70_future_mask$layer.70.5),values(Apr_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_04")
#May
par(mfrow=c(2,4))
hist201205<-hist(c(values(ex_array201205_stack$X201205_20),values(ex_array201205_stack$X201205_40),values(ex_array201205_stack$X201205_60),
                   values(ex_array201205_stack$X201205_80),values(ex_array201205_stack$X201205_100),values(ex_array201205_stack$X201205_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_05")
hist199605<-hist(c(values(ex_array199605_stack$X199605_20),values(ex_array199605_stack$X199605_40),values(ex_array199605_stack$X199605_60),
                   values(ex_array199605_stack$X199605_80),values(ex_array199605_stack$X199605_100),values(ex_array199605_stack$X199605_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_05")
histclim05<-hist(c(values(May_stack$May_20),values(May_stack$May_40),values(May_stack$May_60),
                   values(May_stack$May_80),values(May_stack$May_100),values(May_stack$May_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_05")
hist0_20f05<-hist(c(values(May_0_20_future_mask$layer.1),values(May_0_20_future_mask$layer.2),values(May_0_20_future_mask$layer.3),
                    values(May_0_20_future_mask$layer.4),values(May_0_20_future_mask$layer.5),values(May_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_05")
hist20_40f05<-hist(c(values(May_20_40_future_mask$layer.1),values(May_20_40_future_mask$layer.2),values(May_20_40_future_mask$layer.3),
                     values(May_20_40_future_mask$layer.4),values(May_20_40_future_mask$layer.5),values(May_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_05")
hist40_60f05<-hist(c(values(May_40_60_future_mask$layer.1),values(May_40_60_future_mask$layer.2),values(May_40_60_future_mask$layer.3),
                     values(May_40_60_future_mask$layer.4),values(May_40_60_future_mask$layer.5),values(May_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_05")
hist60_80f05<-hist(c(values(May_60_80_future_mask$layer.1),values(May_60_80_future_mask$layer.2),values(May_60_80_future_mask$layer.3),
                     values(May_60_80_future_mask$layer.4),values(May_60_80_future_mask$layer.5),values(May_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_05")
hist70f05<-hist(c(values(May_70_future_mask$layer.70.1),values(May_70_future_mask$layer.70.2),values(May_70_future_mask$layer.70.3),
                  values(May_70_future_mask$layer.70.4),values(May_70_future_mask$layer.70.5),values(May_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_05")
#Jun
par(mfrow=c(2,4))
hist201206<-hist(c(values(ex_array201206_stack$X201206_20),values(ex_array201206_stack$X201206_40),values(ex_array201206_stack$X201206_60),
                   values(ex_array201206_stack$X201206_80),values(ex_array201206_stack$X201206_100),values(ex_array201206_stack$X201206_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_06")
hist199606<-hist(c(values(ex_array199606_stack$X199606_20),values(ex_array199606_stack$X199606_40),values(ex_array199606_stack$X199606_60),
                   values(ex_array199606_stack$X199606_80),values(ex_array199606_stack$X199606_100),values(ex_array199606_stack$X199606_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_06")
histclim06<-hist(c(values(Jun_stack$Jun_20),values(Jun_stack$Jun_40),values(Jun_stack$Jun_60),
                   values(Jun_stack$Jun_80),values(Jun_stack$Jun_100),values(Jun_stack$Jun_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_06")
hist0_20f06<-hist(c(values(Jun_0_20_future_mask$layer.1),values(Jun_0_20_future_mask$layer.2),values(Jun_0_20_future_mask$layer.3),
                    values(Jun_0_20_future_mask$layer.4),values(Jun_0_20_future_mask$layer.5),values(Jun_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_06")
hist20_40f06<-hist(c(values(Jun_20_40_future_mask$layer.1),values(Jun_20_40_future_mask$layer.2),values(Jun_20_40_future_mask$layer.3),
                     values(Jun_20_40_future_mask$layer.4),values(Jun_20_40_future_mask$layer.5),values(Jun_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_06")
hist40_60f06<-hist(c(values(Jun_40_60_future_mask$layer.1),values(Jun_40_60_future_mask$layer.2),values(Jun_40_60_future_mask$layer.3),
                     values(Jun_40_60_future_mask$layer.4),values(Jun_40_60_future_mask$layer.5),values(Jun_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_06")
hist60_80f06<-hist(c(values(Jun_60_80_future_mask$layer.1),values(Jun_60_80_future_mask$layer.2),values(Jun_60_80_future_mask$layer.3),
                     values(Jun_60_80_future_mask$layer.4),values(Jun_60_80_future_mask$layer.5),values(Jun_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_06")
hist70f06<-hist(c(values(Jun_70_future_mask$layer.70.1),values(Jun_70_future_mask$layer.70.2),values(Jun_70_future_mask$layer.70.3),
                  values(Jun_70_future_mask$layer.70.4),values(Jun_70_future_mask$layer.70.5),values(Jun_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_06")
#Jul
par(mfrow=c(2,4))
hist201207<-hist(c(values(ex_array201207_stack$X201207_20),values(ex_array201207_stack$X201207_40),values(ex_array201207_stack$X201207_60),
                   values(ex_array201207_stack$X201207_80),values(ex_array201207_stack$X201207_100),values(ex_array201207_stack$X201207_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_07")
hist199607<-hist(c(values(ex_array199607_stack$X199607_20),values(ex_array199607_stack$X199607_40),values(ex_array199607_stack$X199607_60),
                   values(ex_array199607_stack$X199607_80),values(ex_array199607_stack$X199607_100),values(ex_array199607_stack$X199607_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_07")
histclim07<-hist(c(values(Jul_stack$Jul_20),values(Jul_stack$Jul_40),values(Jul_stack$Jul_60),
                   values(Jul_stack$Jul_80),values(Jul_stack$Jul_100),values(Jul_stack$Jul_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_07")
hist0_20f07<-hist(c(values(Jul_0_20_future_mask$layer.1),values(Jul_0_20_future_mask$layer.2),values(Jul_0_20_future_mask$layer.3),
                    values(Jul_0_20_future_mask$layer.4),values(Jul_0_20_future_mask$layer.5),values(Jul_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_07")
hist20_40f07<-hist(c(values(Jul_20_40_future_mask$layer.1),values(Jul_20_40_future_mask$layer.2),values(Jul_20_40_future_mask$layer.3),
                     values(Jul_20_40_future_mask$layer.4),values(Jul_20_40_future_mask$layer.5),values(Jul_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_07")
hist40_60f07<-hist(c(values(Jul_40_60_future_mask$layer.1),values(Jul_40_60_future_mask$layer.2),values(Jul_40_60_future_mask$layer.3),
                     values(Jul_40_60_future_mask$layer.4),values(Jul_40_60_future_mask$layer.5),values(Jul_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_07")
hist60_80f07<-hist(c(values(Jul_60_80_future_mask$layer.1),values(Jul_60_80_future_mask$layer.2),values(Jul_60_80_future_mask$layer.3),
                     values(Jul_60_80_future_mask$layer.4),values(Jul_60_80_future_mask$layer.5),values(Jul_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_07")
hist70f07<-hist(c(values(Jul_70_future_mask$layer.70.1),values(Jul_70_future_mask$layer.70.2),values(Jul_70_future_mask$layer.70.3),
                  values(Jul_70_future_mask$layer.70.4),values(Jul_70_future_mask$layer.70.5),values(Jul_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_07")
#Aug
par(mfrow=c(2,4))
hist201208<-hist(c(values(ex_array201208_stack$X201208_20),values(ex_array201208_stack$X201208_40),values(ex_array201208_stack$X201208_60),
                   values(ex_array201208_stack$X201208_80),values(ex_array201208_stack$X201208_100),values(ex_array201208_stack$X201208_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_08")
hist199608<-hist(c(values(ex_array199608_stack$X199608_20),values(ex_array199608_stack$X199608_40),values(ex_array199608_stack$X199608_60),
                   values(ex_array199608_stack$X199608_80),values(ex_array199608_stack$X199608_100),values(ex_array199608_stack$X199608_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_08")
histclim08<-hist(c(values(Aug_stack$Aug_20),values(Aug_stack$Aug_40),values(Aug_stack$Aug_60),
                   values(Aug_stack$Aug_80),values(Aug_stack$Aug_100),values(Aug_stack$Aug_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_08")
hist0_20f08<-hist(c(values(Aug_0_20_future_mask$layer.1),values(Aug_0_20_future_mask$layer.2),values(Aug_0_20_future_mask$layer.3),
                    values(Aug_0_20_future_mask$layer.4),values(Aug_0_20_future_mask$layer.5),values(Aug_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_08")
hist20_40f08<-hist(c(values(Aug_20_40_future_mask$layer.1),values(Aug_20_40_future_mask$layer.2),values(Aug_20_40_future_mask$layer.3),
                     values(Aug_20_40_future_mask$layer.4),values(Aug_20_40_future_mask$layer.5),values(Aug_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_08")
hist40_60f08<-hist(c(values(Aug_40_60_future_mask$layer.1),values(Aug_40_60_future_mask$layer.2),values(Aug_40_60_future_mask$layer.3),
                     values(Aug_40_60_future_mask$layer.4),values(Aug_40_60_future_mask$layer.5),values(Aug_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_08")
hist60_80f08<-hist(c(values(Aug_60_80_future_mask$layer.1),values(Aug_60_80_future_mask$layer.2),values(Aug_60_80_future_mask$layer.3),
                     values(Aug_60_80_future_mask$layer.4),values(Aug_60_80_future_mask$layer.5),values(Aug_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_08")
hist70f08<-hist(c(values(Aug_70_future_mask$layer.70.1),values(Aug_70_future_mask$layer.70.2),values(Aug_70_future_mask$layer.70.3),
                  values(Aug_70_future_mask$layer.70.4),values(Aug_70_future_mask$layer.70.5),values(Aug_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_08")
#Sept
par(mfrow=c(2,4))
hist201209<-hist(c(values(ex_array201209_stack$X201209_20),values(ex_array201209_stack$X201209_40),values(ex_array201209_stack$X201209_60),
                   values(ex_array201209_stack$X201209_80),values(ex_array201209_stack$X201209_100),values(ex_array201209_stack$X201209_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_09")
hist199609<-hist(c(values(ex_array199609_stack$X199609_20),values(ex_array199609_stack$X199609_40),values(ex_array199609_stack$X199609_60),
                   values(ex_array199609_stack$X199609_80),values(ex_array199609_stack$X199609_100),values(ex_array199609_stack$X199609_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_09")
histclim09<-hist(c(values(Sept_stack$Sept_20),values(Sept_stack$Sept_40),values(Sept_stack$Sept_60),
                   values(Sept_stack$Sept_80),values(Sept_stack$Sept_100),values(Sept_stack$Sept_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_09")
hist0_20f09<-hist(c(values(Sept_0_20_future_mask$layer.1),values(Sept_0_20_future_mask$layer.2),values(Sept_0_20_future_mask$layer.3),
                    values(Sept_0_20_future_mask$layer.4),values(Sept_0_20_future_mask$layer.5),values(Sept_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_09")
hist20_40f09<-hist(c(values(Sept_20_40_future_mask$layer.1),values(Sept_20_40_future_mask$layer.2),values(Sept_20_40_future_mask$layer.3),
                     values(Sept_20_40_future_mask$layer.4),values(Sept_20_40_future_mask$layer.5),values(Sept_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_09")
hist40_60f09<-hist(c(values(Sept_40_60_future_mask$layer.1),values(Sept_40_60_future_mask$layer.2),values(Sept_40_60_future_mask$layer.3),
                     values(Sept_40_60_future_mask$layer.4),values(Sept_40_60_future_mask$layer.5),values(Sept_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_09")
hist60_80f09<-hist(c(values(Sept_60_80_future_mask$layer.1),values(Sept_60_80_future_mask$layer.2),values(Sept_60_80_future_mask$layer.3),
                     values(Sept_60_80_future_mask$layer.4),values(Sept_60_80_future_mask$layer.5),values(Sept_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_09")
hist70f09<-hist(c(values(Sept_70_future_mask$layer.70.1),values(Sept_70_future_mask$layer.70.2),values(Sept_70_future_mask$layer.70.3),
                  values(Sept_70_future_mask$layer.70.4),values(Sept_70_future_mask$layer.70.5),values(Sept_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_09")
#Oct
par(mfrow=c(2,4))
hist201210<-hist(c(values(ex_array201210_stack$X201210_20),values(ex_array201210_stack$X201210_40),values(ex_array201210_stack$X201210_60),
                   values(ex_array201210_stack$X201210_80),values(ex_array201210_stack$X201210_100),values(ex_array201210_stack$X201210_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_10")
hist199610<-hist(c(values(ex_array199610_stack$X199610_20),values(ex_array199610_stack$X199610_40),values(ex_array199610_stack$X199610_60),
                   values(ex_array199610_stack$X199610_80),values(ex_array199610_stack$X199610_100),values(ex_array199610_stack$X199610_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_10")
histclim10<-hist(c(values(Oct_stack$Oct_20),values(Oct_stack$Oct_40),values(Oct_stack$Oct_60),
                   values(Oct_stack$Oct_80),values(Oct_stack$Oct_100),values(Oct_stack$Oct_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_10")
hist0_20f10<-hist(c(values(Oct_0_20_future_mask$layer.1),values(Oct_0_20_future_mask$layer.2),values(Oct_0_20_future_mask$layer.3),
                    values(Oct_0_20_future_mask$layer.4),values(Oct_0_20_future_mask$layer.5),values(Oct_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_10")
hist20_40f10<-hist(c(values(Oct_20_40_future_mask$layer.1),values(Oct_20_40_future_mask$layer.2),values(Oct_20_40_future_mask$layer.3),
                     values(Oct_20_40_future_mask$layer.4),values(Oct_20_40_future_mask$layer.5),values(Oct_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_10")
hist40_60f10<-hist(c(values(Oct_40_60_future_mask$layer.1),values(Oct_40_60_future_mask$layer.2),values(Oct_40_60_future_mask$layer.3),
                     values(Oct_40_60_future_mask$layer.4),values(Oct_40_60_future_mask$layer.5),values(Oct_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_10")
hist60_80f10<-hist(c(values(Oct_60_80_future_mask$layer.1),values(Oct_60_80_future_mask$layer.2),values(Oct_60_80_future_mask$layer.3),
                     values(Oct_60_80_future_mask$layer.4),values(Oct_60_80_future_mask$layer.5),values(Oct_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_10")
hist70f10<-hist(c(values(Oct_70_future_mask$layer.70.1),values(Oct_70_future_mask$layer.70.2),values(Oct_70_future_mask$layer.70.3),
                  values(Oct_70_future_mask$layer.70.4),values(Oct_70_future_mask$layer.70.5),values(Oct_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_10")
#Nov
par(mfrow=c(2,4))
hist201211<-hist(c(values(ex_array201211_stack$X201211_20),values(ex_array201211_stack$X201211_40),values(ex_array201211_stack$X201211_60),
                   values(ex_array201211_stack$X201211_80),values(ex_array201211_stack$X201211_100),values(ex_array201211_stack$X201211_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_11")
hist199611<-hist(c(values(ex_array199611_stack$X199611_20),values(ex_array199611_stack$X199611_40),values(ex_array199611_stack$X199611_60),
                   values(ex_array199611_stack$X199611_80),values(ex_array199611_stack$X199611_100),values(ex_array199611_stack$X199611_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_11")
histclim11<-hist(c(values(Mar_stack$Mar_20),values(Mar_stack$Mar_40),values(Mar_stack$Mar_60),
                   values(Mar_stack$Mar_80),values(Mar_stack$Mar_100),values(Mar_stack$Mar_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_11")
hist0_20f11<-hist(c(values(Mar_0_20_future_mask$layer.1),values(Mar_0_20_future_mask$layer.2),values(Mar_0_20_future_mask$layer.3),
                    values(Mar_0_20_future_mask$layer.4),values(Mar_0_20_future_mask$layer.5),values(Mar_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_11")
hist20_40f11<-hist(c(values(Mar_20_40_future_mask$layer.1),values(Mar_20_40_future_mask$layer.2),values(Mar_20_40_future_mask$layer.3),
                     values(Mar_20_40_future_mask$layer.4),values(Mar_20_40_future_mask$layer.5),values(Mar_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_11")
hist40_60f11<-hist(c(values(Mar_40_60_future_mask$layer.1),values(Mar_40_60_future_mask$layer.2),values(Mar_40_60_future_mask$layer.3),
                     values(Mar_40_60_future_mask$layer.4),values(Mar_40_60_future_mask$layer.5),values(Mar_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_11")
hist60_80f11<-hist(c(values(Mar_60_80_future_mask$layer.1),values(Mar_60_80_future_mask$layer.2),values(Mar_60_80_future_mask$layer.3),
                     values(Mar_60_80_future_mask$layer.4),values(Mar_60_80_future_mask$layer.5),values(Mar_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_11")
hist70f11<-hist(c(values(Mar_70_future_mask$layer.70.1),values(Mar_70_future_mask$layer.70.2),values(Mar_70_future_mask$layer.70.3),
                  values(Mar_70_future_mask$layer.70.4),values(Mar_70_future_mask$layer.70.5),values(Mar_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_11")
#Dec
par(mfrow=c(2,4))
hist201212<-hist(c(values(ex_array201212_stack$X201212_20),values(ex_array201212_stack$X201212_40),values(ex_array201212_stack$X201212_60),
                   values(ex_array201212_stack$X201212_80),values(ex_array201212_stack$X201212_100),values(ex_array201212_stack$X201212_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="2012_12")
hist199612<-hist(c(values(ex_array199612_stack$X199612_20),values(ex_array199612_stack$X199612_40),values(ex_array199612_stack$X199612_60),
                   values(ex_array199612_stack$X199612_80),values(ex_array199612_stack$X199612_100),values(ex_array199612_stack$X199612_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="1996_12")
histclim12<-hist(c(values(Dec_stack$Dec_20),values(Dec_stack$Dec_40),values(Dec_stack$Dec_60),
                   values(Dec_stack$Dec_80),values(Dec_stack$Dec_100),values(Dec_stack$Dec_250)),breaks=seq(-2,38,1),ylim=c(0,2000),main="clim_12")
hist0_20f12<-hist(c(values(Dec_0_20_future_mask$layer.1),values(Dec_0_20_future_mask$layer.2),values(Dec_0_20_future_mask$layer.3),
                    values(Dec_0_20_future_mask$layer.4),values(Dec_0_20_future_mask$layer.5),values(Dec_0_20_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="0-20_12")
hist20_40f12<-hist(c(values(Dec_20_40_future_mask$layer.1),values(Dec_20_40_future_mask$layer.2),values(Dec_20_40_future_mask$layer.3),
                     values(Dec_20_40_future_mask$layer.4),values(Dec_20_40_future_mask$layer.5),values(Dec_20_40_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="20-40_12")
hist40_60f12<-hist(c(values(Dec_40_60_future_mask$layer.1),values(Dec_40_60_future_mask$layer.2),values(Dec_40_60_future_mask$layer.3),
                     values(Dec_40_60_future_mask$layer.4),values(Dec_40_60_future_mask$layer.5),values(Dec_40_60_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="40-60_12")
hist60_80f12<-hist(c(values(Dec_60_80_future_mask$layer.1),values(Dec_60_80_future_mask$layer.2),values(Dec_60_80_future_mask$layer.3),
                     values(Dec_60_80_future_mask$layer.4),values(Dec_60_80_future_mask$layer.5),values(Dec_60_80_future_mask$layer.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="60-80_12")
hist70f12<-hist(c(values(Dec_70_future_mask$layer.70.1),values(Dec_70_future_mask$layer.70.2),values(Dec_70_future_mask$layer.70.3),
                  values(Dec_70_future_mask$layer.70.4),values(Dec_70_future_mask$layer.70.5),values(Dec_70_future_mask$layer.70.6)),breaks=seq(-2,38,1),ylim=c(0,2000),main="70_12")


#Mean temp over month at each depth bin
par(mfrow = c(2,3),oma = c(5,4,0,0) + 0.1,mar = c(0,1,2,1) + 0.1)
plot(clim_temps20~c(1:12), type="b",ylim=c(9,28),ylab="Temperature °C)",xlab="",cex.lab=1.5)
text(2.5,27,"0-20 m",cex=2)
mtext(side=2,"Temperature (°C)",line=2.75)
lines(t2012_temps20~c(1:12), type="b",col="darkred")
lines(t1996_temps20~c(1:12), type="b",col="purple")
lines(future0_20_temps20~c(1:12), type="b",col="blue")
lines(future20_40_temps20~c(1:12), type="b",col="green")
lines(future40_60_temps20~c(1:12), type="b",col="orange")
lines(future60_80_temps20~c(1:12), type="b",col="red")
lines(future70_temps20~c(1:12), type="b",col="grey")

plot(clim_temps40~c(1:12), type="b",ylim=c(9,28),ylab="",xlab="")
text(3,27,"20-40 m",cex=2)
lines(t2012_temps40~c(1:12), type="b",col="darkred")
lines(t1996_temps40~c(1:12), type="b",col="purple")
lines(future0_20_temps40~c(1:12), type="b",col="blue")
lines(future20_40_temps40~c(1:12), type="b",col="green")
lines(future40_60_temps40~c(1:12), type="b",col="orange")
lines(future60_80_temps40~c(1:12), type="b",col="red")
lines(future70_temps40~c(1:12), type="b",col="grey")

plot(clim_temps60~c(1:12), type="b",ylim=c(9,28),ylab="",xlab="")
text(3,27,"40-60 m",cex=2)
lines(t2012_temps60~c(1:12), type="b",col="darkred")
lines(t1996_temps60~c(1:12), type="b",col="purple")
lines(future0_20_temps60~c(1:12), type="b",col="blue")
lines(future20_40_temps60~c(1:12), type="b",col="green")
lines(future40_60_temps60~c(1:12), type="b",col="orange")
lines(future60_80_temps60~c(1:12), type="b",col="red")
lines(future70_temps60~c(1:12), type="b",col="grey")

plot(clim_temps80~c(1:12), type="b",ylim=c(9,28),ylab="Temperature °C)",xlab="Month",cex.lab=1.5)
text(3,27,"60-80 m",cex=2)
mtext(side=2,"Temperature (°C)",line=2.75)
mtext(side=1,"Month",line=2.5)
lines(t2012_temps80~c(1:12), type="b",col="darkred")
lines(t1996_temps80~c(1:12), type="b",col="purple")
lines(future0_20_temps80~c(1:12), type="b",col="blue")
lines(future20_40_temps80~c(1:12), type="b",col="green")
lines(future40_60_temps80~c(1:12), type="b",col="orange")
lines(future60_80_temps80~c(1:12), type="b",col="red")
lines(future70_temps80~c(1:12), type="b",col="grey")

plot(clim_temps100~c(1:12), type="b",ylim=c(9,28),xlab="Month",ylab="",cex.lab=1.5)
text(3.5,27,"80-100 m",cex=2)
mtext(side=1,"Month",line=2.5)
lines(t2012_temps100~c(1:12), type="b",col="darkred")
lines(t1996_temps100~c(1:12), type="b",col="purple")
lines(future0_20_temps100~c(1:12), type="b",col="blue")
lines(future20_40_temps100~c(1:12), type="b",col="green")
lines(future40_60_temps100~c(1:12), type="b",col="orange")
lines(future60_80_temps100~c(1:12), type="b",col="red")
lines(future70_temps100~c(1:12), type="b",col="grey")

plot(clim_temps250~c(1:12), type="b",ylim=c(9,28),xlab="Month",ylab="",cex.lab=1.5)
text(3.5,27,"100-250 m",cex=2)
mtext(side=1,"Month",line=2.5)
lines(t2012_temps250~c(1:12), type="b",col="darkred")
lines(t1996_temps250~c(1:12), type="b",col="purple")
lines(future0_20_temps250~c(1:12), type="b",col="blue")
lines(future20_40_temps250~c(1:12), type="b",col="green")
lines(future40_60_temps250~c(1:12), type="b",col="orange")
lines(future60_80_temps250~c(1:12), type="b",col="red")
lines(future70_temps250~c(1:12), type="b",col="grey")
legend(1,24, legend=c("Future Warm","60-80","40-60","20-40","0-20","Warm Year","Cool Year","Climatology"),
       fill=c("grey","red","orange","green","blue","darkred","purple","black"),ncol=2)


###############
#Determine years that are similar to 2012 and warmer than 2012
###############
#Calculate difference between climatology and 2012 for the 0-20, 20-40, 40-60 depth bins
t2012_20difs<-abs(t2012_temps20-clim_temps20)
t2012_40difs<-abs(t2012_temps40-clim_temps40)
t2012_60difs<-abs(t2012_temps60-clim_temps60)

t2012_temps20upper<-t2012_temps20+t2012_20difs
t2012_temps40upper<-t2012_temps40+t2012_40difs
t2012_temps60upper<-t2012_temps60+t2012_60difs

t2012_temps20lower<-t2012_temps20-t2012_20difs
t2012_temps40lower<-t2012_temps40-t2012_40difs
t2012_temps60lower<-t2012_temps60-t2012_60difs

setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/Climate_Deltas")
load("tempmeans20_60.RData")

similarpercentwarmtempsallyears<-NULL
for(i in 1:80){
  year20<-allmeans20[i,]
  year40<-allmeans40[i,]
  year60<-allmeans60[i,]
  
  numtemp20<-length(year20[year20<t2012_temps20upper & year20>t2012_temps20lower])
  numtemp40<-length(year40[year40<t2012_temps40upper & year40>t2012_temps40lower])
  numtemp60<-length(year60[year60<t2012_temps60upper & year60>t2012_temps60lower])
  
  totalwarmtemps<-numtemp20+numtemp40+numtemp60
  percentwarmtemps<-(totalwarmtemps/36)*100 #36 becauase there are 3 years of 12 months, want to get percentage of years that are "similar" to 2012
  similarpercentwarmtempsallyears<-c(similarpercentwarmtempsallyears,percentwarmtemps)
}
similarpercentwarmtempsallyears#years that are similar to 2012 will have percentage greater than 66%
#number of years in each time period that have years similar to 2012
length(similarpercentwarmtempsallyears[1:20][similarpercentwarmtempsallyears[1:20]>66]) #3 years
length(similarpercentwarmtempsallyears[21:40][similarpercentwarmtempsallyears[21:40]>66]) #11 years
length(similarpercentwarmtempsallyears[41:60][similarpercentwarmtempsallyears[41:60]>66]) #1 year
length(similarpercentwarmtempsallyears[61:80][similarpercentwarmtempsallyears[61:80]>66]) #0 years


hotterpercentwarmtempsallyears<-NULL
for(i in 1:80){
  year20<-allmeans20[i,]
  year40<-allmeans40[i,]
  year60<-allmeans60[i,]
  
  numtemp20<-length(year20[year20>t2012_temps20upper])
  numtemp40<-length(year40[year40>t2012_temps40upper])
  numtemp60<-length(year60[year60>t2012_temps60upper])
  
  totalwarmtemps<-numtemp20+numtemp40+numtemp60
  percentwarmtemps<-(totalwarmtemps/36)*100 #36 becauase there are 3 years of 12 months, want to get percentage of years that are "similar" to 2012
  hotterpercentwarmtempsallyears<-c(hotterpercentwarmtempsallyears,percentwarmtemps)
}

hotterpercentwarmtempsallyears#years that are warmer than 2012 will have percentage greater than 66%
#number of years in each time period that have years warmer than 2012
length(hotterpercentwarmtempsallyears[1:20][hotterpercentwarmtempsallyears[1:20]>66]) #0 years
length(hotterpercentwarmtempsallyears[21:40][hotterpercentwarmtempsallyears[21:40]>66]) #1 year
length(hotterpercentwarmtempsallyears[41:60][hotterpercentwarmtempsallyears[41:60]>66]) #15 years
length(hotterpercentwarmtempsallyears[61:80][hotterpercentwarmtempsallyears[61:80]>66]) #20 years

warm_warmer<-data.frame(Similar_to_2012=c(length(similarpercentwarmtempsallyears[1:20][similarpercentwarmtempsallyears[1:20]>66]),
                                          length(similarpercentwarmtempsallyears[21:40][similarpercentwarmtempsallyears[21:40]>66]),
                                          length(similarpercentwarmtempsallyears[41:60][similarpercentwarmtempsallyears[41:60]>66]),
                                          length(similarpercentwarmtempsallyears[61:80][similarpercentwarmtempsallyears[61:80]>66])),
                        Warmer_than_2012=c(length(hotterpercentwarmtempsallyears[1:20][hotterpercentwarmtempsallyears[1:20]>66]),
                                           length(hotterpercentwarmtempsallyears[21:40][hotterpercentwarmtempsallyears[21:40]>66]),
                                           length(hotterpercentwarmtempsallyears[41:60][hotterpercentwarmtempsallyears[41:60]>66]),
                                           length(hotterpercentwarmtempsallyears[61:80][hotterpercentwarmtempsallyears[61:80]>66])))

#Pull out one year that was very warm (year 70)


##############
#Assign Ratio to climatology & extreme years
#for extreme years (2012-warmest, 1996-coldest)
#https://www.ncdc.noaa.gov/cag/global/time-series/0.0,0.0/land_ocean/ann/8/1994-2015
#https://www.nefsc.noaa.gov/ecosys/current-conditions/#longterm
#Mills et al. 2013
##############
#bind depth arrays by month
#climatology (noticed that anywhere where its 10p should be 100p, but leave it for now I know what it means)
Jan_arrays<-abind(Jan_array_depth_sum$temp20_clim$array,Jan_array_depth_sum$temp40_clim$array,
                  Jan_array_depth_sum$temp60_clim$array,Jan_array_depth_sum$temp80_clim$array,
                  Jan_array_depth_sum$temp100_clim$array,Jan_array_depth_sum$temp10p_clim$array,along=3)
Feb_arrays<-abind(Feb_array_depth_sum$temp20_clim$array,Feb_array_depth_sum$temp40_clim$array,
                  Feb_array_depth_sum$temp60_clim$array,Feb_array_depth_sum$temp80_clim$array,
                  Feb_array_depth_sum$temp100_clim$array,Feb_array_depth_sum$temp10p_clim$array,along=3)
Mar_arrays<-abind(Mar_array_depth_sum$temp20_clim$array,Mar_array_depth_sum$temp40_clim$array,
                  Mar_array_depth_sum$temp60_clim$array,Mar_array_depth_sum$temp80_clim$array,
                  Mar_array_depth_sum$temp100_clim$array,Mar_array_depth_sum$temp10p_clim$array,along=3)
Apr_arrays<-abind(Apr_array_depth_sum$temp20_clim$array,Apr_array_depth_sum$temp40_clim$array,
                  Apr_array_depth_sum$temp60_clim$array,Apr_array_depth_sum$temp80_clim$array,
                  Apr_array_depth_sum$temp100_clim$array,Apr_array_depth_sum$temp10p_clim$array,along=3)
May_arrays<-abind(May_array_depth_sum$temp20_clim$array,May_array_depth_sum$temp40_clim$array,
                  May_array_depth_sum$temp60_clim$array,May_array_depth_sum$temp80_clim$array,
                  May_array_depth_sum$temp100_clim$array,May_array_depth_sum$temp10p_clim$array,along=3)
Jun_arrays<-abind(Jun_array_depth_sum$temp20_clim$array,Jun_array_depth_sum$temp40_clim$array,
                  Jun_array_depth_sum$temp60_clim$array,Jun_array_depth_sum$temp80_clim$array,
                  Jun_array_depth_sum$temp100_clim$array,Jun_array_depth_sum$temp10p_clim$array,along=3)
Jul_arrays<-abind(Jul_array_depth_sum$temp20_clim$array,Jul_array_depth_sum$temp40_clim$array,
                  Jul_array_depth_sum$temp60_clim$array,Jul_array_depth_sum$temp80_clim$array,
                  Jul_array_depth_sum$temp100_clim$array,Jul_array_depth_sum$temp10p_clim$array,along=3)
Aug_arrays<-abind(Aug_array_depth_sum$temp20_clim$array,Aug_array_depth_sum$temp40_clim$array,
                  Aug_array_depth_sum$temp60_clim$array,Aug_array_depth_sum$temp80_clim$array,
                  Aug_array_depth_sum$temp100_clim$array,Aug_array_depth_sum$temp10p_clim$array,along=3)
Sept_arrays<-abind(Sept_array_depth_sum$temp20_clim$array,Sept_array_depth_sum$temp40_clim$array,
                   Sept_array_depth_sum$temp60_clim$array,Sept_array_depth_sum$temp80_clim$array,
                   Sept_array_depth_sum$temp100_clim$array,Sept_array_depth_sum$temp10p_clim$array,along=3)
Oct_arrays<-abind(Oct_array_depth_sum$temp20_clim$array,Oct_array_depth_sum$temp40_clim$array,
                  Oct_array_depth_sum$temp60_clim$array,Oct_array_depth_sum$temp80_clim$array,
                  Oct_array_depth_sum$temp100_clim$array,Oct_array_depth_sum$temp10p_clim$array,along=3)
Nov_arrays<-abind(Nov_array_depth_sum$temp20_clim$array,Nov_array_depth_sum$temp40_clim$array,
                  Nov_array_depth_sum$temp60_clim$array,Nov_array_depth_sum$temp80_clim$array,
                  Nov_array_depth_sum$temp100_clim$array,Nov_array_depth_sum$temp10p_clim$array,along=3)
Dec_arrays<-abind(Dec_array_depth_sum$temp20_clim$array,Dec_array_depth_sum$temp40_clim$array,
                  Dec_array_depth_sum$temp60_clim$array,Dec_array_depth_sum$temp80_clim$array,
                  Dec_array_depth_sum$temp100_clim$array,Dec_array_depth_sum$temp10p_clim$array,along=3)

#extreme years
ex_arrays201201_full<-abind(ex_array201201_depth_sum$temp20_clim$array,ex_array201201_depth_sum$temp40_clim$array,
                            ex_array201201_depth_sum$temp60_clim$array,ex_array201201_depth_sum$temp80_clim$array,
                            ex_array201201_depth_sum$temp100_clim$array,ex_array201201_depth_sum$temp100p_clim$array,along=3)
ex_arrays201202_full<-abind(ex_array201202_depth_sum$temp20_clim$array,ex_array201202_depth_sum$temp40_clim$array,
                            ex_array201202_depth_sum$temp60_clim$array,ex_array201202_depth_sum$temp80_clim$array,
                            ex_array201202_depth_sum$temp100_clim$array,ex_array201202_depth_sum$temp100p_clim$array,along=3)
ex_arrays201203_full<-abind(ex_array201203_depth_sum$temp20_clim$array,ex_array201203_depth_sum$temp40_clim$array,
                            ex_array201203_depth_sum$temp60_clim$array,ex_array201203_depth_sum$temp80_clim$array,
                            ex_array201203_depth_sum$temp100_clim$array,ex_array201203_depth_sum$temp100p_clim$array,along=3)
ex_arrays201204_full<-abind(ex_array201204_depth_sum$temp20_clim$array,ex_array201204_depth_sum$temp40_clim$array,
                            ex_array201204_depth_sum$temp60_clim$array,ex_array201204_depth_sum$temp80_clim$array,
                            ex_array201204_depth_sum$temp100_clim$array,ex_array201204_depth_sum$temp100p_clim$array,along=3)
ex_arrays201205_full<-abind(ex_array201205_depth_sum$temp20_clim$array,ex_array201205_depth_sum$temp40_clim$array,
                            ex_array201205_depth_sum$temp60_clim$array,ex_array201205_depth_sum$temp80_clim$array,
                            ex_array201205_depth_sum$temp100_clim$array,ex_array201205_depth_sum$temp100p_clim$array,along=3)
ex_arrays201206_full<-abind(ex_array201206_depth_sum$temp20_clim$array,ex_array201206_depth_sum$temp40_clim$array,
                            ex_array201206_depth_sum$temp60_clim$array,ex_array201206_depth_sum$temp80_clim$array,
                            ex_array201206_depth_sum$temp100_clim$array,ex_array201206_depth_sum$temp100p_clim$array,along=3)
ex_arrays201207_full<-abind(ex_array201207_depth_sum$temp20_clim$array,ex_array201207_depth_sum$temp40_clim$array,
                            ex_array201207_depth_sum$temp60_clim$array,ex_array201207_depth_sum$temp80_clim$array,
                            ex_array201207_depth_sum$temp100_clim$array,ex_array201207_depth_sum$temp100p_clim$array,along=3)
ex_arrays201208_full<-abind(ex_array201208_depth_sum$temp20_clim$array,ex_array201208_depth_sum$temp40_clim$array,
                            ex_array201208_depth_sum$temp60_clim$array,ex_array201208_depth_sum$temp80_clim$array,
                            ex_array201208_depth_sum$temp100_clim$array,ex_array201208_depth_sum$temp100p_clim$array,along=3)
ex_arrays201209_full<-abind(ex_array201209_depth_sum$temp20_clim$array,ex_array201209_depth_sum$temp40_clim$array,
                            ex_array201209_depth_sum$temp60_clim$array,ex_array201209_depth_sum$temp80_clim$array,
                            ex_array201209_depth_sum$temp100_clim$array,ex_array201209_depth_sum$temp100p_clim$array,along=3)
ex_arrays201210_full<-abind(ex_array201210_depth_sum$temp20_clim$array,ex_array201210_depth_sum$temp40_clim$array,
                            ex_array201210_depth_sum$temp60_clim$array,ex_array201210_depth_sum$temp80_clim$array,
                            ex_array201210_depth_sum$temp100_clim$array,ex_array201210_depth_sum$temp100p_clim$array,along=3)
ex_arrays201211_full<-abind(ex_array201211_depth_sum$temp20_clim$array,ex_array201211_depth_sum$temp40_clim$array,
                            ex_array201211_depth_sum$temp60_clim$array,ex_array201211_depth_sum$temp80_clim$array,
                            ex_array201211_depth_sum$temp100_clim$array,ex_array201211_depth_sum$temp100p_clim$array,along=3)
ex_arrays201212_full<-abind(ex_array201212_depth_sum$temp20_clim$array,ex_array201212_depth_sum$temp40_clim$array,
                            ex_array201212_depth_sum$temp60_clim$array,ex_array201212_depth_sum$temp80_clim$array,
                            ex_array201212_depth_sum$temp100_clim$array,ex_array201212_depth_sum$temp100p_clim$array,along=3)
ex_arrays199601_full<-abind(ex_array199601_depth_sum$temp20_clim$array,ex_array199601_depth_sum$temp40_clim$array,
                            ex_array199601_depth_sum$temp60_clim$array,ex_array199601_depth_sum$temp80_clim$array,
                            ex_array199601_depth_sum$temp100_clim$array,ex_array199601_depth_sum$temp100p_clim$array,along=3)
ex_arrays199602_full<-abind(ex_array199602_depth_sum$temp20_clim$array,ex_array199602_depth_sum$temp40_clim$array,
                            ex_array199602_depth_sum$temp60_clim$array,ex_array199602_depth_sum$temp80_clim$array,
                            ex_array199602_depth_sum$temp100_clim$array,ex_array199602_depth_sum$temp100p_clim$array,along=3)
ex_arrays199603_full<-abind(ex_array199603_depth_sum$temp20_clim$array,ex_array199603_depth_sum$temp40_clim$array,
                            ex_array199603_depth_sum$temp60_clim$array,ex_array199603_depth_sum$temp80_clim$array,
                            ex_array199603_depth_sum$temp100_clim$array,ex_array199603_depth_sum$temp100p_clim$array,along=3)
ex_arrays199604_full<-abind(ex_array199604_depth_sum$temp20_clim$array,ex_array199604_depth_sum$temp40_clim$array,
                            ex_array199604_depth_sum$temp60_clim$array,ex_array199604_depth_sum$temp80_clim$array,
                            ex_array199604_depth_sum$temp100_clim$array,ex_array199604_depth_sum$temp100p_clim$array,along=3)
ex_arrays199605_full<-abind(ex_array199605_depth_sum$temp20_clim$array,ex_array199605_depth_sum$temp40_clim$array,
                            ex_array199605_depth_sum$temp60_clim$array,ex_array199605_depth_sum$temp80_clim$array,
                            ex_array199605_depth_sum$temp100_clim$array,ex_array199605_depth_sum$temp100p_clim$array,along=3)
ex_arrays199606_full<-abind(ex_array199606_depth_sum$temp20_clim$array,ex_array199606_depth_sum$temp40_clim$array,
                            ex_array199606_depth_sum$temp60_clim$array,ex_array199606_depth_sum$temp80_clim$array,
                            ex_array199606_depth_sum$temp100_clim$array,ex_array199606_depth_sum$temp100p_clim$array,along=3)
ex_arrays199607_full<-abind(ex_array199607_depth_sum$temp20_clim$array,ex_array199607_depth_sum$temp40_clim$array,
                            ex_array199607_depth_sum$temp60_clim$array,ex_array199607_depth_sum$temp80_clim$array,
                            ex_array199607_depth_sum$temp100_clim$array,ex_array199607_depth_sum$temp100p_clim$array,along=3)
ex_arrays199608_full<-abind(ex_array199608_depth_sum$temp20_clim$array,ex_array199608_depth_sum$temp40_clim$array,
                            ex_array199608_depth_sum$temp60_clim$array,ex_array199608_depth_sum$temp80_clim$array,
                            ex_array199608_depth_sum$temp100_clim$array,ex_array199608_depth_sum$temp100p_clim$array,along=3)
ex_arrays199609_full<-abind(ex_array199609_depth_sum$temp20_clim$array,ex_array199609_depth_sum$temp40_clim$array,
                            ex_array199609_depth_sum$temp60_clim$array,ex_array199609_depth_sum$temp80_clim$array,
                            ex_array199609_depth_sum$temp100_clim$array,ex_array199609_depth_sum$temp100p_clim$array,along=3)
ex_arrays199610_full<-abind(ex_array199610_depth_sum$temp20_clim$array,ex_array199610_depth_sum$temp40_clim$array,
                            ex_array199610_depth_sum$temp60_clim$array,ex_array199610_depth_sum$temp80_clim$array,
                            ex_array199610_depth_sum$temp100_clim$array,ex_array199610_depth_sum$temp100p_clim$array,along=3)
ex_arrays199611_full<-abind(ex_array199611_depth_sum$temp20_clim$array,ex_array199611_depth_sum$temp40_clim$array,
                            ex_array199611_depth_sum$temp60_clim$array,ex_array199611_depth_sum$temp80_clim$array,
                            ex_array199611_depth_sum$temp100_clim$array,ex_array199611_depth_sum$temp100p_clim$array,along=3)
ex_arrays199612_full<-abind(ex_array199612_depth_sum$temp20_clim$array,ex_array199612_depth_sum$temp40_clim$array,
                            ex_array199612_depth_sum$temp60_clim$array,ex_array199612_depth_sum$temp80_clim$array,
                            ex_array199612_depth_sum$temp100_clim$array,ex_array199612_depth_sum$temp100p_clim$array,along=3)


#function assigns a ratio to the climatology or extreme years for each depth bin and each month
temps<-seq(from=1.5,to=33.5,by=0.5)
ratioassign<-function(month_array=NA,ratio=NA,weight=NA){
  month_array[which(month_array<1.5)]<-0 #do this because temp values less than 1.5 end up not be changed to ratio value because I don't have a ratio value less than 1.5 (assumed to be 0 because no cobia found that cold)
  for(i in 1:(length(temps)-1)){
    play_array<-month_array
    month_array[which(play_array>=temps[i] & play_array<temps[i+1])]<-ratio[i]
  }
  month_array[,,1]<-month_array[,,1]*weight[1]
  month_array[,,2]<-month_array[,,2]*weight[2]
  month_array[,,3]<-month_array[,,3]*weight[3]
  month_array[,,4]<-month_array[,,4]*weight[4]
  month_array[,,5]<-month_array[,,5]*weight[5]
  month_array[,,6]<-month_array[,,6]*weight[6]
  
  return(month_array=month_array)
}

#climatology
Jan_ratio_weighted<-ratioassign(month_array=Jan_arrays,ratio=JanRatio,weight=JanDepthWeight)
Feb_ratio_weighted<-ratioassign(month_array=Feb_arrays,ratio=FebRatio,weight=FebDepthWeight)
Mar_ratio_weighted<-ratioassign(month_array=Mar_arrays,ratio=MarRatio,weight=MarDepthWeight)
Apr_ratio_weighted<-ratioassign(month_array=Apr_arrays,ratio=AprRatio,weight=AprDepthWeight)
May_ratio_weighted<-ratioassign(month_array=May_arrays,ratio=MayRatio,weight=MayDepthWeight)
Jun_ratio_weighted<-ratioassign(month_array=Jun_arrays,ratio=JunRatio,weight=JunDepthWeight)
Jul_ratio_weighted<-ratioassign(month_array=Jul_arrays,ratio=JulRatio,weight=JulDepthWeight)
Aug_ratio_weighted<-ratioassign(month_array=Aug_arrays,ratio=AugRatio,weight=AugDepthWeight)
Sept_ratio_weighted<-ratioassign(month_array=Sept_arrays,ratio=SeptRatio,weight=SeptDepthWeight)
Oct_ratio_weighted<-ratioassign(month_array=Oct_arrays,ratio=OctRatio,weight=OctDepthWeight)
Nov_ratio_weighted<-ratioassign(month_array=Nov_arrays,ratio=NovRatio,weight=NovDepthWeight)
Dec_ratio_weighted<-ratioassign(month_array=Dec_arrays,ratio=DecRatio,weight=DecDepthWeight)

#extreme years
ex_arrays201201_weighted<-ratioassign(month_array=ex_arrays201201_full,ratio=JanRatio,weight=JanDepthWeight)
ex_arrays201202_weighted<-ratioassign(month_array=ex_arrays201202_full,ratio=FebRatio,weight=FebDepthWeight)
ex_arrays201203_weighted<-ratioassign(month_array=ex_arrays201203_full,ratio=MarRatio,weight=MarDepthWeight)
ex_arrays201204_weighted<-ratioassign(month_array=ex_arrays201204_full,ratio=AprRatio,weight=AprDepthWeight)
ex_arrays201205_weighted<-ratioassign(month_array=ex_arrays201205_full,ratio=MayRatio,weight=MayDepthWeight)
ex_arrays201206_weighted<-ratioassign(month_array=ex_arrays201206_full,ratio=JunRatio,weight=JunDepthWeight)
ex_arrays201207_weighted<-ratioassign(month_array=ex_arrays201207_full,ratio=JulRatio,weight=JulDepthWeight)
ex_arrays201208_weighted<-ratioassign(month_array=ex_arrays201208_full,ratio=AugRatio,weight=AugDepthWeight)
ex_arrays201209_weighted<-ratioassign(month_array=ex_arrays201209_full,ratio=SeptRatio,weight=SeptDepthWeight)
ex_arrays201210_weighted<-ratioassign(month_array=ex_arrays201210_full,ratio=OctRatio,weight=OctDepthWeight)
ex_arrays201211_weighted<-ratioassign(month_array=ex_arrays201211_full,ratio=NovRatio,weight=NovDepthWeight)
ex_arrays201212_weighted<-ratioassign(month_array=ex_arrays201212_full,ratio=DecRatio,weight=DecDepthWeight)
ex_arrays199601_weighted<-ratioassign(month_array=ex_arrays199601_full,ratio=JanRatio,weight=JanDepthWeight)
ex_arrays199602_weighted<-ratioassign(month_array=ex_arrays199602_full,ratio=FebRatio,weight=FebDepthWeight)
ex_arrays199603_weighted<-ratioassign(month_array=ex_arrays199603_full,ratio=MarRatio,weight=MarDepthWeight)
ex_arrays199604_weighted<-ratioassign(month_array=ex_arrays199604_full,ratio=AprRatio,weight=AprDepthWeight)
ex_arrays199605_weighted<-ratioassign(month_array=ex_arrays199605_full,ratio=MayRatio,weight=MayDepthWeight)
ex_arrays199606_weighted<-ratioassign(month_array=ex_arrays199606_full,ratio=JunRatio,weight=JunDepthWeight)
ex_arrays199607_weighted<-ratioassign(month_array=ex_arrays199607_full,ratio=JulRatio,weight=JulDepthWeight)
ex_arrays199608_weighted<-ratioassign(month_array=ex_arrays199608_full,ratio=AugRatio,weight=AugDepthWeight)
ex_arrays199609_weighted<-ratioassign(month_array=ex_arrays199609_full,ratio=SeptRatio,weight=SeptDepthWeight)
ex_arrays199610_weighted<-ratioassign(month_array=ex_arrays199610_full,ratio=OctRatio,weight=OctDepthWeight)
ex_arrays199611_weighted<-ratioassign(month_array=ex_arrays199611_full,ratio=NovRatio,weight=NovDepthWeight)
ex_arrays199612_weighted<-ratioassign(month_array=ex_arrays199612_full,ratio=DecRatio,weight=DecDepthWeight)



#then add up ratios though water column, value greater than 1 within a cell
#means that area habitat is prefered
#climatology
Jan_total_sum<-apply(Jan_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Feb_total_sum<-apply(Feb_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Mar_total_sum<-apply(Mar_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Apr_total_sum<-apply(Apr_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
May_total_sum<-apply(May_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Jun_total_sum<-apply(Jun_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Jul_total_sum<-apply(Jul_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Aug_total_sum<-apply(Aug_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Sept_total_sum<-apply(Sept_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Oct_total_sum<-apply(Oct_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Nov_total_sum<-apply(Nov_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Dec_total_sum<-apply(Dec_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)

#extreme years
ex_array201201_total_sum<-apply(ex_arrays201201_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201202_total_sum<-apply(ex_arrays201202_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201203_total_sum<-apply(ex_arrays201203_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201204_total_sum<-apply(ex_arrays201204_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201205_total_sum<-apply(ex_arrays201205_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201206_total_sum<-apply(ex_arrays201206_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201207_total_sum<-apply(ex_arrays201207_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201208_total_sum<-apply(ex_arrays201208_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201209_total_sum<-apply(ex_arrays201209_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201210_total_sum<-apply(ex_arrays201210_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201211_total_sum<-apply(ex_arrays201211_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array201212_total_sum<-apply(ex_arrays201212_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199601_total_sum<-apply(ex_arrays199601_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199602_total_sum<-apply(ex_arrays199602_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199603_total_sum<-apply(ex_arrays199603_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199604_total_sum<-apply(ex_arrays199604_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199605_total_sum<-apply(ex_arrays199605_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199606_total_sum<-apply(ex_arrays199606_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199607_total_sum<-apply(ex_arrays199607_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199608_total_sum<-apply(ex_arrays199608_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199609_total_sum<-apply(ex_arrays199609_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199610_total_sum<-apply(ex_arrays199610_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199611_total_sum<-apply(ex_arrays199611_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
ex_array199612_total_sum<-apply(ex_arrays199612_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)


############
#Rasterize
############
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology")
exmp_hycom<-nc_open("expt_53.X_19940101_example.nc")
lon <- ncvar_get(exmp_hycom, varid = "lon") #starts most western point
lat <- ncvar_get(exmp_hycom, varid = "lat") #starts most southern point
temp_exmp<-ncvar_get(exmp_hycom, varid = "water_temp")
nc_close(exmp_hycom)

#for arrayTrim function
lon_lat<-expand.grid(lon=lon, lat=lat)
coordinates(lon_lat)<- ~lon + lat 
proj4string(lon_lat)<-CRS("+init=epsg:4326")#sets it to lat-long
e<-extent(c(min(lon),max(lon),min(lat),max(lat)))+5
r<-raster(e,nrow=338,ncol=207, crs="+proj=longlat +datum=WGS84 +no_defs") # hycom data is 338 lats by 207 longs 

#this function rasterizes ratio array
#when summing further up it turned all NAs to 0 so now need to turn the areas where NAs should be back to NAs
arrayTrim<-function(month_array=NA){
  
  #need to turn all land back to NAs from 0
  #this takes an example array (from 53.X) and finds all locations where all depths have NAs (aka land) 
  #and saves them in a new array (really a matrix) as NAs and all other grid cell locations as 1
  #multiply month_array by NAmat so that all times when there are actually values it will keep those same values, and times when there are only NAs those grid cells will be given NAs
  NAmat<-matrix(0,nrow=length(lon),ncol=length(lat))
  for(i in 1:length(lon)){
    for(j in 1:length(lat)){
      NAmat[i,j]<-ifelse(all(is.na(temp_exmp[i,j,]))==TRUE,NA,1)
    }
    print(i)
  }
  
  new_month_array<-month_array*NAmat
  
  #make values on the west coast of FL NAs
  #cut off is south of 27.7 and west of -80.5
  new_month_array[1:19,1:57]<-NA
  
  #Turn all cells in water east of shelf to NA (use >350m depth as shelf eastern cutoff), NE science center survey uses 366m (100fathoms) as there cutoff for sampling the shelf
  #Need to find first lon where there is a temp value at 400m and make all cells at and beyond that lon NA
  for(l in 1:length(lat)){
    #point of if statements are to find most west position where there are temp values >350m (meaning depth is greater than 350m, which we want to turn to NAs)
    locnonNA<-which(!is.na(temp_exmp[,l,27])) #which lons have NAs >350m
    b<-diff(locnonNA,lag=1) #if any lags are greater than 1 then means there are shallow areas between deep areas (mostly around bahamas)
    if(length(b)>0){ #if b has no length, means that all lon are in shallower water than 400m (gulf of maine)
      if(any(b>1)){ #if b has length greater 0 then need to see if any lags greater than 1, if yes go to next lines
        c<-max(which(b>1))+1 #finds furthest east location where there is a temp deeper than 350m (end of shallow bit)
        minlocnonNA<-locnonNA[c] #actual locates index value for lon (should be first value after the last NA in locnonNA)
      }else{
        minlocnonNA<-min(locnonNA) #if all lags equal 1 then can just select the furthest west location (no shallow areas betwween deep areas)
      }
      new_month_array[minlocnonNA:207,l]<-NA #once find most west position >350m (minlocnonNA) then make all lons east of that position NAs for all depths and days at the given lat
    }else{
      new_month_array[,l]<-new_month_array[,l] #again if all shelf is in shallower than 400m then no need to make anything NAs
    }
  }
  
  #Summarize temp data by depth bins
  month_rasterlayer<-rasterize(lon_lat, r,new_month_array, fun=mean) #units of raster is in degrees
  
  return(month_rasterlayer=month_rasterlayer)
}

#climatology
JanWRationrasterlayer<-arrayTrim(month_array = Jan_total_sum)
FebWRationrasterlayer<-arrayTrim(month_array = Feb_total_sum)
MarWRationrasterlayer<-arrayTrim(month_array = Mar_total_sum)
AprWRationrasterlayer<-arrayTrim(month_array = Apr_total_sum)
MayWRationrasterlayer<-arrayTrim(month_array = May_total_sum)
JunWRationrasterlayer<-arrayTrim(month_array = Jun_total_sum)
JulWRationrasterlayer<-arrayTrim(month_array = Jul_total_sum)
AugWRationrasterlayer<-arrayTrim(month_array = Aug_total_sum)
SeptWRationrasterlayer<-arrayTrim(month_array = Sept_total_sum)
OctWRationrasterlayer<-arrayTrim(month_array = Oct_total_sum)
NovWRationrasterlayer<-arrayTrim(month_array = Nov_total_sum)
DecWRationrasterlayer<-arrayTrim(month_array = Dec_total_sum)

#extreme years
ex201201WRationrasterlayer<-arrayTrim(month_array = ex_array201201_total_sum)
ex201202WRationrasterlayer<-arrayTrim(month_array = ex_array201202_total_sum)
ex201203WRationrasterlayer<-arrayTrim(month_array = ex_array201203_total_sum)
ex201204WRationrasterlayer<-arrayTrim(month_array = ex_array201204_total_sum)
ex201205WRationrasterlayer<-arrayTrim(month_array = ex_array201205_total_sum)
ex201206WRationrasterlayer<-arrayTrim(month_array = ex_array201206_total_sum)
ex201207WRationrasterlayer<-arrayTrim(month_array = ex_array201207_total_sum)
ex201208WRationrasterlayer<-arrayTrim(month_array = ex_array201208_total_sum)
ex201209WRationrasterlayer<-arrayTrim(month_array = ex_array201209_total_sum)
ex201210WRationrasterlayer<-arrayTrim(month_array = ex_array201210_total_sum)
ex201211WRationrasterlayer<-arrayTrim(month_array = ex_array201211_total_sum)
ex201212WRationrasterlayer<-arrayTrim(month_array = ex_array201212_total_sum)
ex199601WRationrasterlayer<-arrayTrim(month_array = ex_array199601_total_sum)
ex199602WRationrasterlayer<-arrayTrim(month_array = ex_array199602_total_sum)
ex199603WRationrasterlayer<-arrayTrim(month_array = ex_array199603_total_sum)
ex199604WRationrasterlayer<-arrayTrim(month_array = ex_array199604_total_sum)
ex199605WRationrasterlayer<-arrayTrim(month_array = ex_array199605_total_sum)
ex199606WRationrasterlayer<-arrayTrim(month_array = ex_array199606_total_sum)
ex199607WRationrasterlayer<-arrayTrim(month_array = ex_array199607_total_sum)
ex199608WRationrasterlayer<-arrayTrim(month_array = ex_array199608_total_sum)
ex199609WRationrasterlayer<-arrayTrim(month_array = ex_array199609_total_sum)
ex199610WRationrasterlayer<-arrayTrim(month_array = ex_array199610_total_sum)
ex199611WRationrasterlayer<-arrayTrim(month_array = ex_array199611_total_sum)
ex199612WRationrasterlayer<-arrayTrim(month_array = ex_array199612_total_sum)




###############
#Model Validation using Acoustic Telemetry
###############
############
#Assign Ratio to Model Validating data and Validate Model
###########
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/Hycom_habitat")
load("Validation_Data_Model.RData")
#also make sure to read in ratio_assign (1233-1235) and arrayTrim funtions (line 1355)
#also read in monthly ratios at the beginning of script (lines 1-72)

#first need to average each month over days, dim: 207,338, #of days
mean_temp20_201707<-apply(temp20_expt_57.7_201707.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201707<-apply(temp40_expt_57.7_201707.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201707<-apply(temp60_expt_57.7_201707.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201707<-apply(temp80_expt_57.7_201707.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201707<-apply(temp100_expt_57.7_201707.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201707<-apply(temp100p_expt_57.7_201707.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201708<-apply(temp20_expt_57.7_201708.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201708<-apply(temp40_expt_57.7_201708.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201708<-apply(temp60_expt_57.7_201708.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201708<-apply(temp80_expt_57.7_201708.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201708<-apply(temp100_expt_57.7_201708.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201708<-apply(temp100p_expt_57.7_201708.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201709<-apply(temp20_expt_57.7_201709.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201709<-apply(temp40_expt_57.7_201709.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201709<-apply(temp60_expt_57.7_201709.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201709<-apply(temp80_expt_57.7_201709.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201709<-apply(temp100_expt_57.7_201709.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201709<-apply(temp100p_expt_57.7_201709.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201710<-apply(temp20_expt_92.9_201710.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201710<-apply(temp40_expt_92.9_201710.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201710<-apply(temp60_expt_92.9_201710.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201710<-apply(temp80_expt_92.9_201710.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201710<-apply(temp100_expt_92.9_201710.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201710<-apply(temp100p_expt_92.9_201710.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201711<-apply(temp20_expt_92.9_201711.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201711<-apply(temp40_expt_92.9_201711.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201711<-apply(temp60_expt_92.9_201711.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201711<-apply(temp80_expt_92.9_201711.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201711<-apply(temp100_expt_92.9_201711.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201711<-apply(temp100p_expt_92.9_201711.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201712<-apply(temp20_expt_92.9_201712.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201712<-apply(temp40_expt_92.9_201712.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201712<-apply(temp60_expt_92.9_201712.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201712<-apply(temp80_expt_92.9_201712.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201712<-apply(temp100_expt_92.9_201712.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201712<-apply(temp100p_expt_92.9_201712.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201801<-apply(temp20_expt_93.0_201801.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201801<-apply(temp40_expt_93.0_201801.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201801<-apply(temp60_expt_93.0_201801.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201801<-apply(temp80_expt_93.0_201801.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201801<-apply(temp100_expt_93.0_201801.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201801<-apply(temp100p_expt_93.0_201801.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201802<-apply(temp20_expt_93.0_201802.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201802<-apply(temp40_expt_93.0_201802.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201802<-apply(temp60_expt_93.0_201802.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201802<-apply(temp80_expt_93.0_201802.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201802<-apply(temp100_expt_93.0_201802.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201802<-apply(temp100p_expt_93.0_201802.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201803<-apply(temp20_expt_93.0_201803.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201803<-apply(temp40_expt_93.0_201803.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201803<-apply(temp60_expt_93.0_201803.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201803<-apply(temp80_expt_93.0_201803.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201803<-apply(temp100_expt_93.0_201803.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201803<-apply(temp100p_expt_93.0_201803.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201804<-apply(temp20_expt_93.0_201804.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201804<-apply(temp40_expt_93.0_201804.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201804<-apply(temp60_expt_93.0_201804.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201804<-apply(temp80_expt_93.0_201804.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201804<-apply(temp100_expt_93.0_201804.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201804<-apply(temp100p_expt_93.0_201804.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201805<-apply(temp20_expt_93.0_201805.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201805<-apply(temp40_expt_93.0_201805.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201805<-apply(temp60_expt_93.0_201805.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201805<-apply(temp80_expt_93.0_201805.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201805<-apply(temp100_expt_93.0_201805.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201805<-apply(temp100p_expt_93.0_201805.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201806<-apply(temp20_expt_93.0_201806.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201806<-apply(temp40_expt_93.0_201806.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201806<-apply(temp60_expt_93.0_201806.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201806<-apply(temp80_expt_93.0_201806.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201806<-apply(temp100_expt_93.0_201806.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201806<-apply(temp100p_expt_93.0_201806.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201807<-apply(temp20_expt_93.0_201807.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201807<-apply(temp40_expt_93.0_201807.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201807<-apply(temp60_expt_93.0_201807.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201807<-apply(temp80_expt_93.0_201807.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201807<-apply(temp100_expt_93.0_201807.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201807<-apply(temp100p_expt_93.0_201807.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201808<-apply(temp20_expt_93.0_201808.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201808<-apply(temp40_expt_93.0_201808.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201808<-apply(temp60_expt_93.0_201808.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201808<-apply(temp80_expt_93.0_201808.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201808<-apply(temp100_expt_93.0_201808.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201808<-apply(temp100p_expt_93.0_201808.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201809<-apply(temp20_expt_93.0_201809.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201809<-apply(temp40_expt_93.0_201809.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201809<-apply(temp60_expt_93.0_201809.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201809<-apply(temp80_expt_93.0_201809.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201809<-apply(temp100_expt_93.0_201809.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201809<-apply(temp100p_expt_93.0_201809.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201810<-apply(temp20_expt_93.0_201810.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201810<-apply(temp40_expt_93.0_201810.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201810<-apply(temp60_expt_93.0_201810.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201810<-apply(temp80_expt_93.0_201810.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201810<-apply(temp100_expt_93.0_201810.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201810<-apply(temp100p_expt_93.0_201810.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201811<-apply(temp20_expt_93.0_201811.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201811<-apply(temp40_expt_93.0_201811.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201811<-apply(temp60_expt_93.0_201811.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201811<-apply(temp80_expt_93.0_201811.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201811<-apply(temp100_expt_93.0_201811.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201811<-apply(temp100p_expt_93.0_201811.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201812<-apply(temp20_expt_93.0_201812.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201812<-apply(temp40_expt_93.0_201812.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201812<-apply(temp60_expt_93.0_201812.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201812<-apply(temp80_expt_93.0_201812.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201812<-apply(temp100_expt_93.0_201812.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201812<-apply(temp100p_expt_93.0_201812.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201901<-apply(temp20_expt_93.0_201901.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201901<-apply(temp40_expt_93.0_201901.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201901<-apply(temp60_expt_93.0_201901.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201901<-apply(temp80_expt_93.0_201901.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201901<-apply(temp100_expt_93.0_201901.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201901<-apply(temp100p_expt_93.0_201901.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201902<-apply(temp20_expt_93.0_201902.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201902<-apply(temp40_expt_93.0_201902.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201902<-apply(temp60_expt_93.0_201902.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201902<-apply(temp80_expt_93.0_201902.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201902<-apply(temp100_expt_93.0_201902.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201902<-apply(temp100p_expt_93.0_201902.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201903<-apply(temp20_expt_93.0_201903.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201903<-apply(temp40_expt_93.0_201903.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201903<-apply(temp60_expt_93.0_201903.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201903<-apply(temp80_expt_93.0_201903.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201903<-apply(temp100_expt_93.0_201903.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201903<-apply(temp100p_expt_93.0_201903.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201904<-apply(temp20_expt_93.0_201904.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201904<-apply(temp40_expt_93.0_201904.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201904<-apply(temp60_expt_93.0_201904.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201904<-apply(temp80_expt_93.0_201904.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201904<-apply(temp100_expt_93.0_201904.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201904<-apply(temp100p_expt_93.0_201904.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201905<-apply(temp20_expt_93.0_201905.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201905<-apply(temp40_expt_93.0_201905.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201905<-apply(temp60_expt_93.0_201905.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201905<-apply(temp80_expt_93.0_201905.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201905<-apply(temp100_expt_93.0_201905.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201905<-apply(temp100p_expt_93.0_201905.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201906<-apply(temp20_expt_93.0_201906.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201906<-apply(temp40_expt_93.0_201906.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201906<-apply(temp60_expt_93.0_201906.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201906<-apply(temp80_expt_93.0_201906.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201906<-apply(temp100_expt_93.0_201906.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201906<-apply(temp100p_expt_93.0_201906.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201907<-apply(temp20_expt_93.0_201907.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201907<-apply(temp40_expt_93.0_201907.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201907<-apply(temp60_expt_93.0_201907.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201907<-apply(temp80_expt_93.0_201907.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201907<-apply(temp100_expt_93.0_201907.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201907<-apply(temp100p_expt_93.0_201907.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp20_201908<-apply(temp20_expt_93.0_201908.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp40_201908<-apply(temp40_expt_93.0_201908.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp60_201908<-apply(temp60_expt_93.0_201908.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp80_201908<-apply(temp80_expt_93.0_201908.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100_201908<-apply(temp100_expt_93.0_201908.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
mean_temp100p_201908<-apply(temp100p_expt_93.0_201908.nc$array[1:207,1:338,],MARGIN=c(1,2),mean,na.rm=TRUE) #takes mean of temps over days for that month
#dim: 207, 338

#bind depth arrays by month
Jul17_arrays<-abind(mean_temp20_201707,mean_temp40_201707,mean_temp60_201707,
                    mean_temp80_201707,mean_temp100_201707,mean_temp100p_201707,along=3)
Aug17_arrays<-abind(mean_temp20_201708,mean_temp40_201708,mean_temp60_201708,
                    mean_temp80_201708,mean_temp100_201708,mean_temp100p_201708,along=3)
Sept17_arrays<-abind(mean_temp20_201709,mean_temp40_201709,mean_temp60_201709,
                     mean_temp80_201709,mean_temp100_201709,mean_temp100p_201709,along=3)
Oct17_arrays<-abind(mean_temp20_201710,mean_temp40_201710,mean_temp60_201710,
                    mean_temp80_201710,mean_temp100_201710,mean_temp100p_201710,along=3)
Nov17_arrays<-abind(mean_temp20_201711,mean_temp40_201711,mean_temp60_201711,
                    mean_temp80_201711,mean_temp100_201711,mean_temp100p_201711,along=3)
Dec17_arrays<-abind(mean_temp20_201712,mean_temp40_201712,mean_temp60_201712,
                    mean_temp80_201712,mean_temp100_201712,mean_temp100p_201712,along=3)
Jan18_arrays<-abind(mean_temp20_201801,mean_temp40_201801,mean_temp60_201801,
                    mean_temp80_201801,mean_temp100_201801,mean_temp100p_201801,along=3)
Feb18_arrays<-abind(mean_temp20_201802,mean_temp40_201802,mean_temp60_201802,
                    mean_temp80_201802,mean_temp100_201802,mean_temp100p_201802,along=3)
Mar18_arrays<-abind(mean_temp20_201803,mean_temp40_201803,mean_temp60_201803,
                    mean_temp80_201803,mean_temp100_201803,mean_temp100p_201803,along=3)
Apr18_arrays<-abind(mean_temp20_201804,mean_temp40_201804,mean_temp60_201804,
                    mean_temp80_201804,mean_temp100_201804,mean_temp100p_201804,along=3)
May18_arrays<-abind(mean_temp20_201805,mean_temp40_201805,mean_temp60_201805,
                    mean_temp80_201805,mean_temp100_201805,mean_temp100p_201805,along=3)
Jun18_arrays<-abind(mean_temp20_201806,mean_temp40_201806,mean_temp60_201806,
                    mean_temp80_201806,mean_temp100_201806,mean_temp100p_201806,along=3)
Jul18_arrays<-abind(mean_temp20_201807,mean_temp40_201807,mean_temp60_201807,
                    mean_temp80_201807,mean_temp100_201807,mean_temp100p_201807,along=3)
Aug18_arrays<-abind(mean_temp20_201808,mean_temp40_201808,mean_temp60_201808,
                    mean_temp80_201808,mean_temp100_201808,mean_temp100p_201808,along=3)
Sept18_arrays<-abind(mean_temp20_201809,mean_temp40_201809,mean_temp60_201809,
                     mean_temp80_201809,mean_temp100_201809,mean_temp100p_201809,along=3)
Oct18_arrays<-abind(mean_temp20_201810,mean_temp40_201810,mean_temp60_201810,
                    mean_temp80_201810,mean_temp100_201810,mean_temp100p_201810,along=3)
Nov18_arrays<-abind(mean_temp20_201811,mean_temp40_201811,mean_temp60_201811,
                    mean_temp80_201811,mean_temp100_201811,mean_temp100p_201811,along=3)
Dec18_arrays<-abind(mean_temp20_201812,mean_temp40_201812,mean_temp60_201812,
                    mean_temp80_201812,mean_temp100_201812,mean_temp100p_201812,along=3)
Jan19_arrays<-abind(mean_temp20_201901,mean_temp40_201901,mean_temp60_201901,
                    mean_temp80_201901,mean_temp100_201901,mean_temp100p_201901,along=3)
Feb19_arrays<-abind(mean_temp20_201902,mean_temp40_201902,mean_temp60_201902,
                    mean_temp80_201902,mean_temp100_201902,mean_temp100p_201902,along=3)
Mar19_arrays<-abind(mean_temp20_201903,mean_temp40_201903,mean_temp60_201903,
                    mean_temp80_201903,mean_temp100_201903,mean_temp100p_201903,along=3)
Apr19_arrays<-abind(mean_temp20_201904,mean_temp40_201904,mean_temp60_201904,
                    mean_temp80_201904,mean_temp100_201904,mean_temp100p_201904,along=3)
May19_arrays<-abind(mean_temp20_201905,mean_temp40_201905,mean_temp60_201905,
                    mean_temp80_201905,mean_temp100_201905,mean_temp100p_201905,along=3)
Jun19_arrays<-abind(mean_temp20_201906,mean_temp40_201906,mean_temp60_201906,
                    mean_temp80_201906,mean_temp100_201906,mean_temp100p_201906,along=3)
Jul19_arrays<-abind(mean_temp20_201907,mean_temp40_201907,mean_temp60_201907,
                    mean_temp80_201907,mean_temp100_201907,mean_temp100p_201907,along=3)
Aug19_arrays<-abind(mean_temp20_201908,mean_temp40_201908,mean_temp60_201908,
                    mean_temp80_201908,mean_temp100_201908,mean_temp100p_201908,along=3)
#dim:207,338,6


#assigns a ratio to the validation data for each depth bin and each month
Jul17_ratio_weighted<-ratioassign(month_array=Jul17_arrays,ratio=JulRatio,weight=JulDepthWeight)
Aug17_ratio_weighted<-ratioassign(month_array=Aug17_arrays,ratio=AugRatio,weight=AugDepthWeight)
Sept17_ratio_weighted<-ratioassign(month_array=Sept17_arrays,ratio=SeptRatio,weight=SeptDepthWeight)
Oct17_ratio_weighted<-ratioassign(month_array=Oct17_arrays,ratio=OctRatio,weight=OctDepthWeight)
Nov17_ratio_weighted<-ratioassign(month_array=Nov17_arrays,ratio=NovRatio,weight=NovDepthWeight)
Dec17_ratio_weighted<-ratioassign(month_array=Dec17_arrays,ratio=DecRatio,weight=DecDepthWeight)
Jan18_ratio_weighted<-ratioassign(month_array=Jan18_arrays,ratio=JanRatio,weight=JanDepthWeight)
Feb18_ratio_weighted<-ratioassign(month_array=Feb18_arrays,ratio=FebRatio,weight=FebDepthWeight)
Mar18_ratio_weighted<-ratioassign(month_array=Mar18_arrays,ratio=MarRatio,weight=MarDepthWeight)
Apr18_ratio_weighted<-ratioassign(month_array=Apr18_arrays,ratio=AprRatio,weight=AprDepthWeight)
May18_ratio_weighted<-ratioassign(month_array=May18_arrays,ratio=MayRatio,weight=MayDepthWeight)
Jun18_ratio_weighted<-ratioassign(month_array=Jun18_arrays,ratio=JunRatio,weight=JunDepthWeight)
Jul18_ratio_weighted<-ratioassign(month_array=Jul18_arrays,ratio=JulRatio,weight=JulDepthWeight)
Aug18_ratio_weighted<-ratioassign(month_array=Aug18_arrays,ratio=AugRatio,weight=AugDepthWeight)
Sept18_ratio_weighted<-ratioassign(month_array=Sept18_arrays,ratio=SeptRatio,weight=SeptDepthWeight)
Oct18_ratio_weighted<-ratioassign(month_array=Oct18_arrays,ratio=OctRatio,weight=OctDepthWeight)
Nov18_ratio_weighted<-ratioassign(month_array=Nov18_arrays,ratio=NovRatio,weight=NovDepthWeight)
Dec18_ratio_weighted<-ratioassign(month_array=Dec18_arrays,ratio=DecRatio,weight=DecDepthWeight)
Jan19_ratio_weighted<-ratioassign(month_array=Jan19_arrays,ratio=JanRatio,weight=JanDepthWeight)
Feb19_ratio_weighted<-ratioassign(month_array=Feb19_arrays,ratio=FebRatio,weight=FebDepthWeight)
Mar19_ratio_weighted<-ratioassign(month_array=Mar19_arrays,ratio=MarRatio,weight=MarDepthWeight)
Apr19_ratio_weighted<-ratioassign(month_array=Apr19_arrays,ratio=AprRatio,weight=AprDepthWeight)
May19_ratio_weighted<-ratioassign(month_array=May19_arrays,ratio=MayRatio,weight=MayDepthWeight)
Jun19_ratio_weighted<-ratioassign(month_array=Jun19_arrays,ratio=JunRatio,weight=JunDepthWeight)
Jul19_ratio_weighted<-ratioassign(month_array=Jul19_arrays,ratio=JulRatio,weight=JulDepthWeight)
Aug19_ratio_weighted<-ratioassign(month_array=Aug19_arrays,ratio=AugRatio,weight=AugDepthWeight)
#dim:207,338,6

#then add up ratios though water column, value greater than 1 within a cell
#means that area habitat is prefered
Jul17_total_sum<-apply(Jul17_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Aug17_total_sum<-apply(Aug17_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Sept17_total_sum<-apply(Sept17_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Oct17_total_sum<-apply(Oct17_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Nov17_total_sum<-apply(Nov17_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Dec17_total_sum<-apply(Dec17_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Jan18_total_sum<-apply(Jan18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Feb18_total_sum<-apply(Feb18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Mar18_total_sum<-apply(Mar18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Apr18_total_sum<-apply(Apr18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
May18_total_sum<-apply(May18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Jun18_total_sum<-apply(Jun18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Jul18_total_sum<-apply(Jul18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Aug18_total_sum<-apply(Aug18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Sept18_total_sum<-apply(Sept18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Oct18_total_sum<-apply(Oct18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Nov18_total_sum<-apply(Nov18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Dec18_total_sum<-apply(Dec18_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Jan19_total_sum<-apply(Jan19_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Feb19_total_sum<-apply(Feb19_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Mar19_total_sum<-apply(Mar19_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Apr19_total_sum<-apply(Apr19_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
May19_total_sum<-apply(May19_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Jun19_total_sum<-apply(Jun19_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Jul19_total_sum<-apply(Jul19_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
Aug19_total_sum<-apply(Aug19_ratio_weighted[1:207,1:338,],MARGIN = c(1,2),sum,na.rm=TRUE)
#dim: 207, 338

#when summing further up it turned all NAs to 0 so now need to turn the areas where NAs should be back to NAs
#this is for arrayTrim function if not already read in
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology")
exmp_hycom<-nc_open("expt_53.X_19940101_example.nc")
lon <- ncvar_get(exmp_hycom, varid = "lon") #starts most western point
lat <- ncvar_get(exmp_hycom, varid = "lat") #starts most southern point
temp_exmp<-ncvar_get(exmp_hycom, varid = "water_temp")
nc_close(exmp_hycom)

#for arrayTrim function
lon_lat<-expand.grid(lon=lon, lat=lat)
coordinates(lon_lat)<- ~lon + lat 
proj4string(lon_lat)<-CRS("+init=epsg:4326")#sets it to lat-long
e<-extent(c(min(lon),max(lon),min(lat),max(lat)))+5
r<-raster(e,nrow=338,ncol=207, crs="+proj=longlat +datum=WGS84 +no_defs") # hycom data is 338 lats by 207 longs 

Jul17WRationrasterlayer<-arrayTrim(month_array = Jul17_total_sum)
Aug17WRationrasterlayer<-arrayTrim(month_array = Aug17_total_sum)
Sept17WRationrasterlayer<-arrayTrim(month_array = Sept17_total_sum)
Oct17WRationrasterlayer<-arrayTrim(month_array = Oct17_total_sum)
Nov17WRationrasterlayer<-arrayTrim(month_array = Nov17_total_sum)
Dec17WRationrasterlayer<-arrayTrim(month_array = Dec17_total_sum)
Jan18WRationrasterlayer<-arrayTrim(month_array = Jan18_total_sum)
Feb18WRationrasterlayer<-arrayTrim(month_array = Feb18_total_sum)
Mar18WRationrasterlayer<-arrayTrim(month_array = Mar18_total_sum)
Apr18WRationrasterlayer<-arrayTrim(month_array = Apr18_total_sum)
May18WRationrasterlayer<-arrayTrim(month_array = May18_total_sum)
Jun18WRationrasterlayer<-arrayTrim(month_array = Jun18_total_sum)
Jul18WRationrasterlayer<-arrayTrim(month_array = Jul18_total_sum)
Aug18WRationrasterlayer<-arrayTrim(month_array = Aug18_total_sum)
Sept18WRationrasterlayer<-arrayTrim(month_array = Sept18_total_sum)
Oct18WRationrasterlayer<-arrayTrim(month_array = Oct18_total_sum)
Nov18WRationrasterlayer<-arrayTrim(month_array = Nov18_total_sum)
Dec18WRationrasterlayer<-arrayTrim(month_array = Dec18_total_sum)
Jan19WRationrasterlayer<-arrayTrim(month_array = Jan19_total_sum)
Feb19WRationrasterlayer<-arrayTrim(month_array = Feb19_total_sum)
Mar19WRationrasterlayer<-arrayTrim(month_array = Mar19_total_sum)
Apr19WRationrasterlayer<-arrayTrim(month_array = Apr19_total_sum)
May19WRationrasterlayer<-arrayTrim(month_array = May19_total_sum)
Jun19WRationrasterlayer<-arrayTrim(month_array = Jun19_total_sum)
Jul19WRationrasterlayer<-arrayTrim(month_array = Jul19_total_sum)
Aug19WRationrasterlayer<-arrayTrim(month_array = Aug19_total_sum)


setwd("~/Documents/Miscellaneous_R/Miscellaneous/Cobia_Movements/VIMS_cobia_detections")
detects<-read.csv("20200110_ALL_ARRAYS_VIMS_cobia_detections.csv")
detects$Date_Time_Local<-as.POSIXct(detects$Date_Time_Local, format="%Y-%m-%d %H:%M:%S")
detects$Month_Year_Local<-strftime(detects$Date_Time_Local,format="%Y-%m")
detects$Year<-strftime(detects$Date_Time_Local,format="%Y")


detects$ID<-sapply(strsplit(as.character(detects$Transmitter),split= "-",fixed=TRUE),'[',3)
#get rid of detections from fish that had archival tags in them 18448, 18433, 18437, 18441, 18430

detects_no_archive<-detects[-which(detects$ID=="18448" | detects$ID=="18433" | detects$ID=="18437" | detects$ID=="18441" | detects$ID=="18430"),]

#weird missing lat/lon for station B5 for May and June 2019, so just copy and paste B5 lat/lon from other rows
detects_no_archive$Latitude[which(detects_no_archive$Station.Name=="B5")]<-36.97968
detects_no_archive$Longitude[which(detects_no_archive$Station.Name=="B5")]<- -76.00367

validFunc<-function(month_year=NA,raster_month_year=NA){
  d<-which(detects_no_archive$Month_Year_Local==month_year)
  d1<-detects_no_archive[d,]
  yearS<-d1$Year[1]
  #alpha calculation (prop of detections within preferred area)
  detects_coords<-subset(d1, select=c(Longitude,Latitude))
  detects_sp_coords<-SpatialPoints(detects_coords,proj4string = CRS("+init=epsg:4326"))
  detects_raster<-extract(raster_month_year,detects_sp_coords) #some NAs, (e.g. detects at mouth of James River)
  if(any(is.na(detects_raster))){
    detects_raster<-detects_raster[-which(is.na(detects_raster)==TRUE)]
  }else{
    detects_raster<-detects_raster
  }
  num_preferred<-length(which(detects_raster>1))
  num_total<-length(detects_raster)
  alpha<-num_preferred/num_total
  #beta calculation (prop of receivers [that received detections] within preferred area)
  stations_month_year<-unique(d1$Station.Name)
  stat_latRec<-NULL
  stat_lonRec<-NULL
  for(j in stations_month_year){
    s<-which(d1$Station.Name==j)
    s1<-d1[s,]
    stat_lat<-s1$Latitude[1]
    stat_lon<-s1$Longitude[1]
    stat_latRec<-c(stat_latRec,stat_lat)
    stat_lonRec<-c(stat_lonRec,stat_lon)
  }
  stat_lon_lat<-data.frame(lon=stat_lonRec,lat=stat_latRec)
  station_sp_coords<-SpatialPoints(stat_lon_lat,proj4string = CRS("+init=epsg:4326"))
  station_raster<-extract(raster_month_year,station_sp_coords) #some NAs, (e.g. stations at mouth of James River)
  if(any(is.na(station_raster))){
    station_raster<-station_raster[-which(is.na(station_raster)==TRUE)]
  }else{
    station_raster<-station_raster
  }
  num_stats_preferred<-length(which(station_raster>1)) #number of stations in preferred habitat that month and year that received a detection that month
  dY<-which(detects$Year==yearS)
  dY1<-detects[dY,]
  num_stats_total<-length(unique(dY1$Station.Name))#total number of stations that recieved detections that year
  beta<-num_stats_preferred/num_stats_total #if beta equals 0 then that means no receivers that received a detection was in preferred habitat
  validate_metric<-alpha/beta #validation metric, >1 means model is informative
  
  return(c(alpha=alpha, beta=beta, validate_metric=validate_metric))
}

Jul17_validation_metric<-validFunc(month_year="2017-07",raster_month_year=Jul17WRationrasterlayer)
Aug17_validation_metric<-validFunc(month_year="2017-08",raster_month_year=Aug17WRationrasterlayer)
Sept17_validation_metric<-validFunc(month_year="2017-09",raster_month_year=Sept17WRationrasterlayer)
Oct17_validation_metric<-validFunc(month_year="2017-10",raster_month_year=Oct17WRationrasterlayer)
Nov17_validation_metric<-validFunc(month_year="2017-11",raster_month_year=Nov17WRationrasterlayer)
Dec17_validation_metric<-validFunc(month_year="2017-12",raster_month_year=Dec17WRationrasterlayer)
#Jan18_validation_metric<-validFunc(month_year="2018-01",raster_month_year=Jan18WRationrasterlayer) #no detections
#Feb18_validation_metric<-validFunc(month_year="2018-02",raster_month_year=Feb18WRationrasterlayer) #no detections
#Mar18_validation_metric<-validFunc(month_year="2018-03",raster_month_year=Mar18WRationrasterlayer) #no detections
#Apr18_validation_metric<-validFunc(month_year="2018-04",raster_month_year=Apr18WRationrasterlayer) #no detections
May18_validation_metric<-validFunc(month_year="2018-05",raster_month_year=May18WRationrasterlayer)
Jun18_validation_metric<-validFunc(month_year="2018-06",raster_month_year=Jun18WRationrasterlayer)
Jul18_validation_metric<-validFunc(month_year="2018-07",raster_month_year=Jul18WRationrasterlayer)
Aug18_validation_metric<-validFunc(month_year="2018-08",raster_month_year=Aug18WRationrasterlayer)
Sept18_validation_metric<-validFunc(month_year="2018-09",raster_month_year=Sept18WRationrasterlayer)
Oct18_validation_metric<-validFunc(month_year="2018-10",raster_month_year=Oct18WRationrasterlayer)
Nov18_validation_metric<-validFunc(month_year="2018-11",raster_month_year=Nov18WRationrasterlayer)
Dec18_validation_metric<-validFunc(month_year="2018-12",raster_month_year=Dec18WRationrasterlayer)
Jan19_validation_metric<-validFunc(month_year="2019-01",raster_month_year=Jan19WRationrasterlayer)
Feb19_validation_metric<-validFunc(month_year="2019-02",raster_month_year=Feb19WRationrasterlayer)
Mar19_validation_metric<-validFunc(month_year="2019-03",raster_month_year=Mar19WRationrasterlayer)
Apr19_validation_metric<-validFunc(month_year="2019-04",raster_month_year=Apr19WRationrasterlayer)
May19_validation_metric<-validFunc(month_year="2019-05",raster_month_year=May19WRationrasterlayer)
Jun19_validation_metric<-validFunc(month_year="2019-06",raster_month_year=Jun19WRationrasterlayer)

#make table
val2017<-c("-","-","-","-","-","-",round(Jul17_validation_metric[3],1),round(Aug17_validation_metric[3],1),
           round(Sept17_validation_metric[3],1),round(Oct17_validation_metric[3],1),round(Nov17_validation_metric[3],1),
           round(Dec17_validation_metric[3],1))
val2018<-c("-","-","-","-",round(May18_validation_metric[3],1),round(Jun18_validation_metric[3],1),round(Jul18_validation_metric[3],1),
           round(Aug18_validation_metric[3],1),round(Sept18_validation_metric[3],1),round(Oct18_validation_metric[3],1),
           round(Nov18_validation_metric[3],1),round(Dec18_validation_metric[3],1))
val2019<-c(round(Jan19_validation_metric[3],1),round(Feb19_validation_metric[3],1),round(Mar19_validation_metric[3],1),round(Apr19_validation_metric[3],1),
           round(May19_validation_metric[3],1),round(Jun19_validation_metric[3],1),"-","-","-","-","-","-")
valmonth<-c("January","February","March","April","May","June","July","August","September","October","November","December")
valid_dat<-data.frame(Month=valmonth,"2017"=val2017,"2018"=val2018,"2019"=val2019)
colnames(valid_dat)<-c("Month","2017","2018","2019")


#for values for the paper
forarchivestudy<-detects_no_archive[-which(detects_no_archive$Month_Year_Local=="2019-07"|detects_no_archive$Month_Year_Local=="2019-08"|detects_no_archive$Month_Year_Local=="2019-09"|detects_no_archive$Month_Year_Local=="2019-10"|
                                             detects_no_archive$Month_Year_Local=="2019-11"),]
#number of indivudals
length(unique(forarchivestudy$Transmitter))
#number of receiver stations
length(unique(forarchivestudy$Station.Name))


setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology")
write.csv(valid_dat, file="validation_table.csv")


#############
#Mask for only US waters
#############
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/FederalAndStateWaters/GIS_shp")
uswaters<-readOGR(".","US_State_Fed_waters")
uswaters<-spTransform(uswaters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ME_waters<-readOGR(".","ME_waters")
ME_waters<-spTransform(ME_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
NH_waters<-readOGR(".","NH_waters")
NH_waters<-spTransform(NH_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
MA_waters<-readOGR(".","MA_waters")
MA_waters<-spTransform(MA_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
RI_waters<-readOGR(".","RI_waters")
RI_waters<-spTransform(RI_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
CT_waters<-readOGR(".","CT_waters")
CT_waters<-spTransform(CT_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
NY_waters<-readOGR(".","NY_waters")
NY_waters<-spTransform(NY_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
NJ_waters<-readOGR(".","NJ_waters")
NJ_waters<-spTransform(NJ_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
DE_waters<-readOGR(".","DE_waters")
DE_waters<-spTransform(DE_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
MD_waters<-readOGR(".","MD_waters")
MD_waters<-spTransform(MD_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
VA_waters<-readOGR(".","VA_waters")
VA_waters<-spTransform(VA_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
NC_waters<-readOGR(".","NC_waters")
NC_waters<-spTransform(NC_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
SC_waters<-readOGR(".","SC_waters")
SC_waters<-spTransform(SC_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
GA_waters<-readOGR(".","GA_waters")
GA_waters<-spTransform(GA_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
FL_waters<-readOGR(".","FL_waters")
FL_waters<-spTransform(FL_waters,"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#climatology
JanWRationrasterlayer<-mask(JanWRationrasterlayer,mask=uswaters)
FebWRationrasterlayer<-mask(FebWRationrasterlayer,mask=uswaters)
MarWRationrasterlayer<-mask(MarWRationrasterlayer,mask=uswaters)
AprWRationrasterlayer<-mask(AprWRationrasterlayer,mask=uswaters)
MayWRationrasterlayer<-mask(MayWRationrasterlayer,mask=uswaters)
JunWRationrasterlayer<-mask(JunWRationrasterlayer,mask=uswaters)
JulWRationrasterlayer<-mask(JulWRationrasterlayer,mask=uswaters)
AugWRationrasterlayer<-mask(AugWRationrasterlayer,mask=uswaters)
SeptWRationrasterlayer<-mask(SeptWRationrasterlayer,mask=uswaters)
OctWRationrasterlayer<-mask(OctWRationrasterlayer,mask=uswaters)
NovWRationrasterlayer<-mask(NovWRationrasterlayer,mask=uswaters)
DecWRationrasterlayer<-mask(DecWRationrasterlayer,mask=uswaters)

#extreme years
ex201201WRationrasterlayer<-mask(ex201201WRationrasterlayer,mask=uswaters)
ex201202WRationrasterlayer<-mask(ex201202WRationrasterlayer,mask=uswaters)
ex201203WRationrasterlayer<-mask(ex201203WRationrasterlayer,mask=uswaters)
ex201204WRationrasterlayer<-mask(ex201204WRationrasterlayer,mask=uswaters)
ex201205WRationrasterlayer<-mask(ex201205WRationrasterlayer,mask=uswaters)
ex201206WRationrasterlayer<-mask(ex201206WRationrasterlayer,mask=uswaters)
ex201207WRationrasterlayer<-mask(ex201207WRationrasterlayer,mask=uswaters)
ex201208WRationrasterlayer<-mask(ex201208WRationrasterlayer,mask=uswaters)
ex201209WRationrasterlayer<-mask(ex201209WRationrasterlayer,mask=uswaters)
ex201210WRationrasterlayer<-mask(ex201210WRationrasterlayer,mask=uswaters)
ex201211WRationrasterlayer<-mask(ex201211WRationrasterlayer,mask=uswaters)
ex201212WRationrasterlayer<-mask(ex201212WRationrasterlayer,mask=uswaters)
ex199601WRationrasterlayer<-mask(ex199601WRationrasterlayer,mask=uswaters)
ex199602WRationrasterlayer<-mask(ex199602WRationrasterlayer,mask=uswaters)
ex199603WRationrasterlayer<-mask(ex199603WRationrasterlayer,mask=uswaters)
ex199604WRationrasterlayer<-mask(ex199604WRationrasterlayer,mask=uswaters)
ex199605WRationrasterlayer<-mask(ex199605WRationrasterlayer,mask=uswaters)
ex199606WRationrasterlayer<-mask(ex199606WRationrasterlayer,mask=uswaters)
ex199607WRationrasterlayer<-mask(ex199607WRationrasterlayer,mask=uswaters)
ex199608WRationrasterlayer<-mask(ex199608WRationrasterlayer,mask=uswaters)
ex199609WRationrasterlayer<-mask(ex199609WRationrasterlayer,mask=uswaters)
ex199610WRationrasterlayer<-mask(ex199610WRationrasterlayer,mask=uswaters)
ex199611WRationrasterlayer<-mask(ex199611WRationrasterlayer,mask=uswaters)
ex199612WRationrasterlayer<-mask(ex199612WRationrasterlayer,mask=uswaters)


#write rasters
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/Raster_outputs")
writeRaster(JanWRationrasterlayer,file="JanWRationrasterlayer.tiff",overwrite=T)
writeRaster(FebWRationrasterlayer,file="FebWRationrasterlayer.tiff",overwrite=T)
writeRaster(MarWRationrasterlayer,file="MarWRationrasterlayer.tiff",overwrite=T)
writeRaster(AprWRationrasterlayer,file="AprWRationrasterlayer.tiff",overwrite=T)
writeRaster(MayWRationrasterlayer,file="MayWRationrasterlayer.tiff",overwrite=T)
writeRaster(JunWRationrasterlayer,file="JunWRationrasterlayer.tiff",overwrite=T)
writeRaster(JulWRationrasterlayer,file="JulWRationrasterlayer.tiff",overwrite=T)
writeRaster(AugWRationrasterlayer,file="AugWRationrasterlayer.tiff",overwrite=T)
writeRaster(SeptWRationrasterlayer,file="SeptWRationrasterlayer.tiff",overwrite=T)
writeRaster(OctWRationrasterlayer,file="OctWRationrasterlayer.tiff",overwrite=T)
writeRaster(NovWRationrasterlayer,file="NovWRationrasterlayer.tiff",overwrite=T)
writeRaster(DecWRationrasterlayer,file="DecWRationrasterlayer.tiff",overwrite=T)

writeRaster(ex201201WRationrasterlayer,file="ex201201WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201202WRationrasterlayer,file="ex201202WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201203WRationrasterlayer,file="ex201203WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201204WRationrasterlayer,file="ex201204WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201205WRationrasterlayer,file="ex201205WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201206WRationrasterlayer,file="ex201206WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201207WRationrasterlayer,file="ex201207WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201208WRationrasterlayer,file="ex201208WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201209WRationrasterlayer,file="ex201209WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201210WRationrasterlayer,file="ex201210WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201211WRationrasterlayer,file="ex201211WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex201212WRationrasterlayer,file="ex201212WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199601WRationrasterlayer,file="ex199601WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199602WRationrasterlayer,file="ex199602WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199603WRationrasterlayer,file="ex199603WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199604WRationrasterlayer,file="ex199604WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199605WRationrasterlayer,file="ex199605WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199606WRationrasterlayer,file="ex199606WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199607WRationrasterlayer,file="ex199607WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199608WRationrasterlayer,file="ex199608WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199609WRationrasterlayer,file="ex199609WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199610WRationrasterlayer,file="ex199610WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199611WRationrasterlayer,file="ex199611WRationrasterlayer.tiff",overwrite=T)
writeRaster(ex199612WRationrasterlayer,file="ex199612WRationrasterlayer.tiff",overwrite=T)



##################
#Assign Ratio to Future climatology & extreme year
##################
#make sure habitat ratios and depth weights are read in from above (lines 13-73)
#also make sure to read in uswaters shapefile (lines 1559-1591)
setwd("/Users/danielcrear/Documents/PhD Project/SDM/SDM/Regional_climatology/Climate_Deltas")
load("ShelfMonthClimatologyFuture.RData")
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/Climate_Deltas/ExtremeFutureTemps")
load("future70extremetemps.RData")


#need to read in a raster that has already been masked because when summing further down it turns all NAs to 0
#so need to mask it afterwards with a raster that already has the correct boundaries (including the uswaters mask)
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/Raster_outputs")
JanWRationrasterlayer<-raster("JanWRationrasterlayer.tif")


#function assigns a ratio to the climatology or extreme years for each depth bin and each month
#then sums up ratios for all layers together to get one surface
#and only includes US waters by masking with US waters shapefile
temps<-seq(from=1.5,to=33.5,by=0.5)
ratioassignRaster<-function(month_raster=NA,ratio=NA,weight=NA){
  month_raster[month_raster<1.5]<-0 #do this because temp values less than 1.5 end up not be changed to ratio value because I don't have a ratio value less than 1.5 (assumed to be 0 because no cobia found that cold)
  month_raster[month_raster>33.5]<-0 #do this because temp values greater than 33.5 end up not be changed to ratio value because I don't have a ratio value greater than 1.5 (assumed to be 0 because no cobia found that warm)
  for(i in 1:(length(temps)-1)){
    play_raster<-month_raster
    month_raster[play_raster>=temps[i] & play_raster<temps[i+1]]<-ratio[i]
  }
  month_raster$layer.1<-month_raster$layer.1*weight[1]
  month_raster$layer.2<-month_raster$layer.2*weight[2]
  month_raster$layer.3<-month_raster$layer.3*weight[3]
  month_raster$layer.4<-month_raster$layer.4*weight[4]
  month_raster$layer.5<-month_raster$layer.5*weight[5]
  month_raster$layer.6<-month_raster$layer.6*weight[6]
  
  #add ratios for all layers together to get one surface
  month_raster_weighted<-calc(month_raster, fun=sum, na.rm=T)
  
  #only include US waters by masking with US waters shapefile
  #instead of US water shapefile use a raster that has everything correct because when add up all layers above, all NAs turn to 0s and boundarys disappear
  month_raster_weighted<-mask(month_raster_weighted,mask=JanWRationrasterlayer)
  
  return(month_raster_weighted=month_raster_weighted)
}

#Future climatology
Jan0_20WRationrasterlayer<-ratioassignRaster(month_raster=Jan_0_20_future,ratio=JanRatio,weight=JanDepthWeight)
Jan20_40WRationrasterlayer<-ratioassignRaster(month_raster=Jan_20_40_future,ratio=JanRatio,weight=JanDepthWeight)
Jan40_60WRationrasterlayer<-ratioassignRaster(month_raster=Jan_40_60_future,ratio=JanRatio,weight=JanDepthWeight)
Jan60_80WRationrasterlayer<-ratioassignRaster(month_raster=Jan_60_80_future,ratio=JanRatio,weight=JanDepthWeight)
Feb0_20WRationrasterlayer<-ratioassignRaster(month_raster=Feb_0_20_future,ratio=FebRatio,weight=FebDepthWeight)
Feb20_40WRationrasterlayer<-ratioassignRaster(month_raster=Feb_20_40_future,ratio=FebRatio,weight=FebDepthWeight)
Feb40_60WRationrasterlayer<-ratioassignRaster(month_raster=Feb_40_60_future,ratio=FebRatio,weight=FebDepthWeight)
Feb60_80WRationrasterlayer<-ratioassignRaster(month_raster=Feb_60_80_future,ratio=FebRatio,weight=FebDepthWeight)
Mar0_20WRationrasterlayer<-ratioassignRaster(month_raster=Mar_0_20_future,ratio=MarRatio,weight=MarDepthWeight)
Mar20_40WRationrasterlayer<-ratioassignRaster(month_raster=Mar_20_40_future,ratio=MarRatio,weight=MarDepthWeight)
Mar40_60WRationrasterlayer<-ratioassignRaster(month_raster=Mar_40_60_future,ratio=MarRatio,weight=MarDepthWeight)
Mar60_80WRationrasterlayer<-ratioassignRaster(month_raster=Mar_60_80_future,ratio=MarRatio,weight=MarDepthWeight)
Apr0_20WRationrasterlayer<-ratioassignRaster(month_raster=Apr_0_20_future,ratio=AprRatio,weight=AprDepthWeight)
Apr20_40WRationrasterlayer<-ratioassignRaster(month_raster=Apr_20_40_future,ratio=AprRatio,weight=AprDepthWeight)
Apr40_60WRationrasterlayer<-ratioassignRaster(month_raster=Apr_40_60_future,ratio=AprRatio,weight=AprDepthWeight)
Apr60_80WRationrasterlayer<-ratioassignRaster(month_raster=Apr_60_80_future,ratio=AprRatio,weight=AprDepthWeight)
May0_20WRationrasterlayer<-ratioassignRaster(month_raster=May_0_20_future,ratio=MayRatio,weight=MayDepthWeight)
May20_40WRationrasterlayer<-ratioassignRaster(month_raster=May_20_40_future,ratio=MayRatio,weight=MayDepthWeight)
May40_60WRationrasterlayer<-ratioassignRaster(month_raster=May_40_60_future,ratio=MayRatio,weight=MayDepthWeight)
May60_80WRationrasterlayer<-ratioassignRaster(month_raster=May_60_80_future,ratio=MayRatio,weight=MayDepthWeight)
Jun0_20WRationrasterlayer<-ratioassignRaster(month_raster=Jun_0_20_future,ratio=JunRatio,weight=JunDepthWeight)
Jun20_40WRationrasterlayer<-ratioassignRaster(month_raster=Jun_20_40_future,ratio=JunRatio,weight=JunDepthWeight)
Jun40_60WRationrasterlayer<-ratioassignRaster(month_raster=Jun_40_60_future,ratio=JunRatio,weight=JunDepthWeight)
Jun60_80WRationrasterlayer<-ratioassignRaster(month_raster=Jun_60_80_future,ratio=JunRatio,weight=JunDepthWeight)
Jul0_20WRationrasterlayer<-ratioassignRaster(month_raster=Jul_0_20_future,ratio=JulRatio,weight=JulDepthWeight)
Jul20_40WRationrasterlayer<-ratioassignRaster(month_raster=Jul_20_40_future,ratio=JulRatio,weight=JulDepthWeight)
Jul40_60WRationrasterlayer<-ratioassignRaster(month_raster=Jul_40_60_future,ratio=JulRatio,weight=JulDepthWeight)
Jul60_80WRationrasterlayer<-ratioassignRaster(month_raster=Jul_60_80_future,ratio=JulRatio,weight=JulDepthWeight)
Aug0_20WRationrasterlayer<-ratioassignRaster(month_raster=Aug_0_20_future,ratio=AugRatio,weight=AugDepthWeight)
Aug20_40WRationrasterlayer<-ratioassignRaster(month_raster=Aug_20_40_future,ratio=AugRatio,weight=AugDepthWeight)
Aug40_60WRationrasterlayer<-ratioassignRaster(month_raster=Aug_40_60_future,ratio=AugRatio,weight=AugDepthWeight)
Aug60_80WRationrasterlayer<-ratioassignRaster(month_raster=Aug_60_80_future,ratio=AugRatio,weight=AugDepthWeight)
Sept0_20WRationrasterlayer<-ratioassignRaster(month_raster=Sept_0_20_future,ratio=SeptRatio,weight=SeptDepthWeight)
Sept20_40WRationrasterlayer<-ratioassignRaster(month_raster=Sept_20_40_future,ratio=SeptRatio,weight=SeptDepthWeight)
Sept40_60WRationrasterlayer<-ratioassignRaster(month_raster=Sept_40_60_future,ratio=SeptRatio,weight=SeptDepthWeight)
Sept60_80WRationrasterlayer<-ratioassignRaster(month_raster=Sept_60_80_future,ratio=SeptRatio,weight=SeptDepthWeight)
Oct0_20WRationrasterlayer<-ratioassignRaster(month_raster=Oct_0_20_future,ratio=OctRatio,weight=OctDepthWeight)
Oct20_40WRationrasterlayer<-ratioassignRaster(month_raster=Oct_20_40_future,ratio=OctRatio,weight=OctDepthWeight)
Oct40_60WRationrasterlayer<-ratioassignRaster(month_raster=Oct_40_60_future,ratio=OctRatio,weight=OctDepthWeight)
Oct60_80WRationrasterlayer<-ratioassignRaster(month_raster=Oct_60_80_future,ratio=OctRatio,weight=OctDepthWeight)
Nov0_20WRationrasterlayer<-ratioassignRaster(month_raster=Nov_0_20_future,ratio=NovRatio,weight=NovDepthWeight)
Nov20_40WRationrasterlayer<-ratioassignRaster(month_raster=Nov_20_40_future,ratio=NovRatio,weight=NovDepthWeight)
Nov40_60WRationrasterlayer<-ratioassignRaster(month_raster=Nov_40_60_future,ratio=NovRatio,weight=NovDepthWeight)
Nov60_80WRationrasterlayer<-ratioassignRaster(month_raster=Nov_60_80_future,ratio=NovRatio,weight=NovDepthWeight)
Dec0_20WRationrasterlayer<-ratioassignRaster(month_raster=Dec_0_20_future,ratio=DecRatio,weight=DecDepthWeight)
Dec20_40WRationrasterlayer<-ratioassignRaster(month_raster=Dec_20_40_future,ratio=DecRatio,weight=DecDepthWeight)
Dec40_60WRationrasterlayer<-ratioassignRaster(month_raster=Dec_40_60_future,ratio=DecRatio,weight=DecDepthWeight)
Dec60_80WRationrasterlayer<-ratioassignRaster(month_raster=Dec_60_80_future,ratio=DecRatio,weight=DecDepthWeight)

#Future Extreme year (year 70)
names(future70extremetempsJan)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsFeb)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsMar)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsApr)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsMay)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsJun)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsJul)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsAug)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsSept)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsOct)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsNov)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")
names(future70extremetempsDec)<-c("layer.1","layer.2","layer.3","layer.4","layer.5","layer.6")

Jan70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsJan,ratio=JanRatio,weight=JanDepthWeight)
Feb70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsFeb,ratio=FebRatio,weight=FebDepthWeight)
Mar70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsMar,ratio=MarRatio,weight=MarDepthWeight)
Apr70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsApr,ratio=AprRatio,weight=AprDepthWeight)
May70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsMay,ratio=MayRatio,weight=MayDepthWeight)
Jun70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsJun,ratio=JunRatio,weight=JunDepthWeight)
Jul70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsJul,ratio=JulRatio,weight=JulDepthWeight)
Aug70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsAug,ratio=AugRatio,weight=AugDepthWeight)
Sept70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsSept,ratio=SeptRatio,weight=SeptDepthWeight)
Oct70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsOct,ratio=OctRatio,weight=OctDepthWeight)
Nov70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsNov,ratio=NovRatio,weight=NovDepthWeight)
Dec70WRationrasterlayer<-ratioassignRaster(month_raster=future70extremetempsDec,ratio=DecRatio,weight=DecDepthWeight)

#write Rasters
setwd("/Users/danielcrear/Documents/PhD Project/SDM/SDM/Regional_climatology/Raster_outputs_future")
writeRaster(Jan0_20WRationrasterlayer,file="Jan0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Jan20_40WRationrasterlayer,file="Jan20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Jan40_60WRationrasterlayer,file="Jan40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Jan60_80WRationrasterlayer,file="Jan60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Feb0_20WRationrasterlayer,file="Feb0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Feb20_40WRationrasterlayer,file="Feb20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Feb40_60WRationrasterlayer,file="Feb40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Feb60_80WRationrasterlayer,file="Feb60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Mar0_20WRationrasterlayer,file="Mar0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Mar20_40WRationrasterlayer,file="Mar20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Mar40_60WRationrasterlayer,file="Mar40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Mar60_80WRationrasterlayer,file="Mar60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Apr0_20WRationrasterlayer,file="Apr0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Apr20_40WRationrasterlayer,file="Apr20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Apr40_60WRationrasterlayer,file="Apr40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Apr60_80WRationrasterlayer,file="Apr60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(May0_20WRationrasterlayer,file="May0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(May20_40WRationrasterlayer,file="May20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(May40_60WRationrasterlayer,file="May40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(May60_80WRationrasterlayer,file="May60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Jun0_20WRationrasterlayer,file="Jun0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Jun20_40WRationrasterlayer,file="Jun20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Jun40_60WRationrasterlayer,file="Jun40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Jun60_80WRationrasterlayer,file="Jun60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Jul0_20WRationrasterlayer,file="Jul0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Jul20_40WRationrasterlayer,file="Jul20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Jul40_60WRationrasterlayer,file="Jul40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Jul60_80WRationrasterlayer,file="Jul60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Aug0_20WRationrasterlayer,file="Aug0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Aug20_40WRationrasterlayer,file="Aug20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Aug40_60WRationrasterlayer,file="Aug40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Aug60_80WRationrasterlayer,file="Aug60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Sept0_20WRationrasterlayer,file="Sept0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Sept20_40WRationrasterlayer,file="Sept20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Sept40_60WRationrasterlayer,file="Sept40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Sept60_80WRationrasterlayer,file="Sept60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Oct0_20WRationrasterlayer,file="Oct0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Oct20_40WRationrasterlayer,file="Oct20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Oct40_60WRationrasterlayer,file="Oct40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Oct60_80WRationrasterlayer,file="Oct60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Nov0_20WRationrasterlayer,file="Nov0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Nov20_40WRationrasterlayer,file="Nov20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Nov40_60WRationrasterlayer,file="Nov40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Nov60_80WRationrasterlayer,file="Nov60_80WRationrasterlayer.tif",overwrite=T)
writeRaster(Dec0_20WRationrasterlayer,file="Dec0_20WRationrasterlayer.tif",overwrite=T)
writeRaster(Dec20_40WRationrasterlayer,file="Dec20_40WRationrasterlayer.tif",overwrite=T)
writeRaster(Dec40_60WRationrasterlayer,file="Dec40_60WRationrasterlayer.tif",overwrite=T)
writeRaster(Dec60_80WRationrasterlayer,file="Dec60_80WRationrasterlayer.tif",overwrite=T)


#write Future extreme rasters
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/Raster_outputs_future")
writeRaster(Jan70WRationrasterlayer,file="Jan70WRationrasterlayer.tif",overwrite=T)
writeRaster(Feb70WRationrasterlayer,file="Feb70WRationrasterlayer.tif",overwrite=T)
writeRaster(Mar70WRationrasterlayer,file="Mar70WRationrasterlayer.tif",overwrite=T)
writeRaster(Apr70WRationrasterlayer,file="Apr70WRationrasterlayer.tif",overwrite=T)
writeRaster(May70WRationrasterlayer,file="May70WRationrasterlayer.tif",overwrite=T)
writeRaster(Jun70WRationrasterlayer,file="Jun70WRationrasterlayer.tif",overwrite=T)
writeRaster(Jul70WRationrasterlayer,file="Jul70WRationrasterlayer.tif",overwrite=T)
writeRaster(Aug70WRationrasterlayer,file="Aug70WRationrasterlayer.tif",overwrite=T)
writeRaster(Sept70WRationrasterlayer,file="Sept70WRationrasterlayer.tif",overwrite=T)
writeRaster(Oct70WRationrasterlayer,file="Oct70WRationrasterlayer.tif",overwrite=T)
writeRaster(Nov70WRationrasterlayer,file="Nov70WRationrasterlayer.tif",overwrite=T)
writeRaster(Dec70WRationrasterlayer,file="Dec70WRationrasterlayer.tif",overwrite=T)


#################
#Calculate Metrics
#Read rasters in if not already in
################
#make sure have states read in from above (lines 1889)
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology")
exmp_hycom<-nc_open("expt_53.X_19940101_example.nc")
lon <- ncvar_get(exmp_hycom, varid = "lon") #starts most western point
lat <- ncvar_get(exmp_hycom, varid = "lat") #starts most southern point
temp_exmp<-ncvar_get(exmp_hycom, varid = "water_temp")
nc_close(exmp_hycom)
setwd("~/Documents/PhD Project/SDM/SDM/Regional_climatology/Raster_outputs")
#climatology
JanWRationrasterlayer<-raster("JanWRationrasterlayer.tif")
FebWRationrasterlayer<-raster("FebWRationrasterlayer.tif")
MarWRationrasterlayer<-raster("MarWRationrasterlayer.tif")
AprWRationrasterlayer<-raster("AprWRationrasterlayer.tif")
MayWRationrasterlayer<-raster("MayWRationrasterlayer.tif")
JunWRationrasterlayer<-raster("JunWRationrasterlayer.tif")
JulWRationrasterlayer<-raster("JulWRationrasterlayer.tif")
AugWRationrasterlayer<-raster("AugWRationrasterlayer.tif")
SeptWRationrasterlayer<-raster("SeptWRationrasterlayer.tif")
OctWRationrasterlayer<-raster("OctWRationrasterlayer.tif")
NovWRationrasterlayer<-raster("NovWRationrasterlayer.tif")
DecWRationrasterlayer<-raster("DecWRationrasterlayer.tif")
#extreme years
ex201201WRationrasterlayer<-raster("ex201201WRationrasterlayer.tif")
ex201202WRationrasterlayer<-raster("ex201202WRationrasterlayer.tif")
ex201203WRationrasterlayer<-raster("ex201203WRationrasterlayer.tif")
ex201204WRationrasterlayer<-raster("ex201204WRationrasterlayer.tif")
ex201205WRationrasterlayer<-raster("ex201205WRationrasterlayer.tif")
ex201206WRationrasterlayer<-raster("ex201206WRationrasterlayer.tif")
ex201207WRationrasterlayer<-raster("ex201207WRationrasterlayer.tif")
ex201208WRationrasterlayer<-raster("ex201208WRationrasterlayer.tif")
ex201209WRationrasterlayer<-raster("ex201209WRationrasterlayer.tif")
ex201210WRationrasterlayer<-raster("ex201210WRationrasterlayer.tif")
ex201211WRationrasterlayer<-raster("ex201211WRationrasterlayer.tif")
ex201212WRationrasterlayer<-raster("ex201212WRationrasterlayer.tif")
ex199601WRationrasterlayer<-raster("ex199601WRationrasterlayer.tif")
ex199602WRationrasterlayer<-raster("ex199602WRationrasterlayer.tif")
ex199603WRationrasterlayer<-raster("ex199603WRationrasterlayer.tif")
ex199604WRationrasterlayer<-raster("ex199604WRationrasterlayer.tif")
ex199605WRationrasterlayer<-raster("ex199605WRationrasterlayer.tif")
ex199606WRationrasterlayer<-raster("ex199606WRationrasterlayer.tif")
ex199607WRationrasterlayer<-raster("ex199607WRationrasterlayer.tif")
ex199608WRationrasterlayer<-raster("ex199608WRationrasterlayer.tif")
ex199609WRationrasterlayer<-raster("ex199609WRationrasterlayer.tif")
ex199610WRationrasterlayer<-raster("ex199610WRationrasterlayer.tif")
ex199611WRationrasterlayer<-raster("ex199611WRationrasterlayer.tif")
ex199612WRationrasterlayer<-raster("ex199612WRationrasterlayer.tif")
#future climatologies
setwd("/Users/danielcrear/Documents/PhD Project/SDM/SDM/Regional_climatology/Raster_outputs_future")
Jan0_20WRationrasterlayer<-raster("Jan0_20WRationrasterlayer.tif")
Jan20_40WRationrasterlayer<-raster("Jan20_40WRationrasterlayer.tif")
Jan40_60WRationrasterlayer<-raster("Jan40_60WRationrasterlayer.tif")
Jan60_80WRationrasterlayer<-raster("Jan60_80WRationrasterlayer.tif")
Feb0_20WRationrasterlayer<-raster("Feb0_20WRationrasterlayer.tif")
Feb20_40WRationrasterlayer<-raster("Feb20_40WRationrasterlayer.tif")
Feb40_60WRationrasterlayer<-raster("Feb40_60WRationrasterlayer.tif")
Feb60_80WRationrasterlayer<-raster("Feb60_80WRationrasterlayer.tif")
Mar0_20WRationrasterlayer<-raster("Mar0_20WRationrasterlayer.tif")
Mar20_40WRationrasterlayer<-raster("Mar20_40WRationrasterlayer.tif")
Mar40_60WRationrasterlayer<-raster("Mar40_60WRationrasterlayer.tif")
Mar60_80WRationrasterlayer<-raster("Mar60_80WRationrasterlayer.tif")
Apr0_20WRationrasterlayer<-raster("Apr0_20WRationrasterlayer.tif")
Apr20_40WRationrasterlayer<-raster("Apr20_40WRationrasterlayer.tif")
Apr40_60WRationrasterlayer<-raster("Apr40_60WRationrasterlayer.tif")
Apr60_80WRationrasterlayer<-raster("Apr60_80WRationrasterlayer.tif")
May0_20WRationrasterlayer<-raster("May0_20WRationrasterlayer.tif")
May20_40WRationrasterlayer<-raster("May20_40WRationrasterlayer.tif")
May40_60WRationrasterlayer<-raster("May40_60WRationrasterlayer.tif")
May60_80WRationrasterlayer<-raster("May60_80WRationrasterlayer.tif")
Jun0_20WRationrasterlayer<-raster("Jun0_20WRationrasterlayer.tif")
Jun20_40WRationrasterlayer<-raster("Jun20_40WRationrasterlayer.tif")
Jun40_60WRationrasterlayer<-raster("Jun40_60WRationrasterlayer.tif")
Jun60_80WRationrasterlayer<-raster("Jun60_80WRationrasterlayer.tif")
Jul0_20WRationrasterlayer<-raster("Jul0_20WRationrasterlayer.tif")
Jul20_40WRationrasterlayer<-raster("Jul20_40WRationrasterlayer.tif")
Jul40_60WRationrasterlayer<-raster("Jul40_60WRationrasterlayer.tif")
Jul60_80WRationrasterlayer<-raster("Jul60_80WRationrasterlayer.tif")
Aug0_20WRationrasterlayer<-raster("Aug0_20WRationrasterlayer.tif")
Aug20_40WRationrasterlayer<-raster("Aug20_40WRationrasterlayer.tif")
Aug40_60WRationrasterlayer<-raster("Aug40_60WRationrasterlayer.tif")
Aug60_80WRationrasterlayer<-raster("Aug60_80WRationrasterlayer.tif")
Sept0_20WRationrasterlayer<-raster("Sept0_20WRationrasterlayer.tif")
Sept20_40WRationrasterlayer<-raster("Sept20_40WRationrasterlayer.tif")
Sept40_60WRationrasterlayer<-raster("Sept40_60WRationrasterlayer.tif")
Sept60_80WRationrasterlayer<-raster("Sept60_80WRationrasterlayer.tif")
Oct0_20WRationrasterlayer<-raster("Oct0_20WRationrasterlayer.tif")
Oct20_40WRationrasterlayer<-raster("Oct20_40WRationrasterlayer.tif")
Oct40_60WRationrasterlayer<-raster("Oct40_60WRationrasterlayer.tif")
Oct60_80WRationrasterlayer<-raster("Oct60_80WRationrasterlayer.tif")
Nov0_20WRationrasterlayer<-raster("Nov0_20WRationrasterlayer.tif")
Nov20_40WRationrasterlayer<-raster("Nov20_40WRationrasterlayer.tif")
Nov40_60WRationrasterlayer<-raster("Nov40_60WRationrasterlayer.tif")
Nov60_80WRationrasterlayer<-raster("Nov60_80WRationrasterlayer.tif")
Dec0_20WRationrasterlayer<-raster("Dec0_20WRationrasterlayer.tif")
Dec20_40WRationrasterlayer<-raster("Dec20_40WRationrasterlayer.tif")
Dec40_60WRationrasterlayer<-raster("Dec40_60WRationrasterlayer.tif")
Dec60_80WRationrasterlayer<-raster("Dec60_80WRationrasterlayer.tif")

Jan70WRationrasterlayer<-raster("Jan70WRationrasterlayer.tif")
Feb70WRationrasterlayer<-raster("Feb70WRationrasterlayer.tif")
Mar70WRationrasterlayer<-raster("Mar70WRationrasterlayer.tif")
Apr70WRationrasterlayer<-raster("Apr70WRationrasterlayer.tif")
May70WRationrasterlayer<-raster("May70WRationrasterlayer.tif")
Jun70WRationrasterlayer<-raster("Jun70WRationrasterlayer.tif")
Jul70WRationrasterlayer<-raster("Jul70WRationrasterlayer.tif")
Aug70WRationrasterlayer<-raster("Aug70WRationrasterlayer.tif")
Sept70WRationrasterlayer<-raster("Sept70WRationrasterlayer.tif")
Oct70WRationrasterlayer<-raster("Oct70WRationrasterlayer.tif")
Nov70WRationrasterlayer<-raster("Nov70WRationrasterlayer.tif")
Dec70WRationrasterlayer<-raster("Dec70WRationrasterlayer.tif")


##################
#State water use
##################
#function calculates habitat value for each state
sumvalFunc<-function(raster_month=NA,waters=NA){
  state_raster<-mask(raster_month,mask=waters) #masks cell only if the cell's centroid is in the mask (state shapefile)
  statearearaster<-area(state_raster,na.rm=TRUE) #area is in km2
  statearearasteriso<-statearearaster[!is.na(statearearaster)]#only isolates (gives vector of) cell area values without NAs
  staterastervals<-values(state_raster) #values of all cells
  staterastervalsiso<-staterastervals[!is.na(staterastervals)] #only isolates (gives vector of) cell values (ratio value) without NAs
  statetotalval<-statearearasteriso*staterastervalsiso #gives area value x ratio value
  statesumval<-sum(statetotalval) #total sum of all "values" on shelf
  
  return(statesumval=statesumval) #sum within that state's waters
}

#this function calculate overall habitat value along US shelf and divides it by habitat value for each state to get a percentage
stateHabitatUse<-function(raster_month=NA){
  raster_month[raster_month <= 1] <- NA
  arearaster<-area(raster_month,na.rm=TRUE) #area is in km2
  arearasteriso<-arearaster[!is.na(arearaster)]#only isolates (gives vector of) cell area values without NAs
  rastervals<-values(raster_month) #values of all cells
  rastervalsiso<-rastervals[!is.na(rastervals)] #only isolates (gives vector of) cell values (ratio value) without NAs
  totalval<-arearasteriso*rastervalsiso #gives area value x ratio value
  sumval<-sum(totalval) #total sum of all "values" on shelf
  
  MEstatesumval<-sumvalFunc(waters=ME_waters,raster_month=raster_month)
  NHstatesumval<-sumvalFunc(waters=NH_waters,raster_month=raster_month)
  MAstatesumval<-sumvalFunc(waters=MA_waters,raster_month=raster_month)
  RIstatesumval<-sumvalFunc(waters=RI_waters,raster_month=raster_month)
  CTstatesumval<-sumvalFunc(waters=CT_waters,raster_month=raster_month)
  NYstatesumval<-sumvalFunc(waters=NY_waters,raster_month=raster_month)
  NJstatesumval<-sumvalFunc(waters=NJ_waters,raster_month=raster_month)
  DEstatesumval<-sumvalFunc(waters=DE_waters,raster_month=raster_month)
  MDstatesumval<-sumvalFunc(waters=MD_waters,raster_month=raster_month)
  VAstatesumval<-sumvalFunc(waters=VA_waters,raster_month=raster_month)
  NCstatesumval<-sumvalFunc(waters=NC_waters,raster_month=raster_month)
  SCstatesumval<-sumvalFunc(waters=SC_waters,raster_month=raster_month)
  GAstatesumval<-sumvalFunc(waters=GA_waters,raster_month=raster_month)
  FLstatesumval<-sumvalFunc(waters=FL_waters,raster_month=raster_month)
  
  MEstate_percent<-MEstatesumval/sumval
  NHstate_percent<-NHstatesumval/sumval
  MAstate_percent<-MAstatesumval/sumval
  RIstate_percent<-RIstatesumval/sumval
  CTstate_percent<-CTstatesumval/sumval
  NYstate_percent<-NYstatesumval/sumval
  NJstate_percent<-NJstatesumval/sumval
  DEstate_percent<-DEstatesumval/sumval
  MDstate_percent<-MDstatesumval/sumval
  VAstate_percent<-VAstatesumval/sumval
  NCstate_percent<-NCstatesumval/sumval
  SCstate_percent<-SCstatesumval/sumval
  GAstate_percent<-GAstatesumval/sumval
  FLstate_percent<-FLstatesumval/sumval
  
  states_percents<-c(MEstate_percent,NHstate_percent,MAstate_percent,RIstate_percent,CTstate_percent,
                     NYstate_percent,NJstate_percent,DEstate_percent,MDstate_percent,VAstate_percent,
                     NCstate_percent,SCstate_percent,GAstate_percent,FLstate_percent,sumval)
  return(states_percents=states_percents)
}

States<-c("ME","NH","MA","RI","CT","NY","NJ","DE","MD","VA","NC","SC","GA","FL")
#climatology
Jan_state_percents<-stateHabitatUse(raster_month=JanWRationrasterlayer)
Feb_state_percents<-stateHabitatUse(raster_month=FebWRationrasterlayer)
Mar_state_percents<-stateHabitatUse(raster_month=MarWRationrasterlayer)
Apr_state_percents<-stateHabitatUse(raster_month=AprWRationrasterlayer)
May_state_percents<-stateHabitatUse(raster_month=MayWRationrasterlayer)
Jun_state_percents<-stateHabitatUse(raster_month=JunWRationrasterlayer)
Jul_state_percents<-stateHabitatUse(raster_month=JulWRationrasterlayer)
Aug_state_percents<-stateHabitatUse(raster_month=AugWRationrasterlayer)
Sept_state_percents<-stateHabitatUse(raster_month=SeptWRationrasterlayer)
Oct_state_percents<-stateHabitatUse(raster_month=OctWRationrasterlayer)
Nov_state_percents<-stateHabitatUse(raster_month=NovWRationrasterlayer)
Dec_state_percents<-stateHabitatUse(raster_month=DecWRationrasterlayer)

State_water_percents<-data.frame(States=States,Jan=Jan_state_percents[1:14]*100,Feb=Feb_state_percents[1:14]*100,Mar=Mar_state_percents[1:14]*100,
                                 Apr=Apr_state_percents[1:14]*100,May=May_state_percents[1:14]*100,Jun=Jun_state_percents[1:14]*100,Jul=Jul_state_percents[1:14]*100,
                                 Aug=Aug_state_percents[1:14]*100,Sept=Sept_state_percents[1:14]*100,Oct=Oct_state_percents[1:14]*100,Nov=Nov_state_percents[1:14]*100,Dec=Dec_state_percents[1:14]*100)
revState_water_percents<-State_water_percents[nrow(State_water_percents):1,]

par(mfrow=c(2,6))
barplot(height=revState_water_percents$Jan,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jan Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Feb,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Feb Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Mar,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Mar Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Apr,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Apr Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$May,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="May Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Jun,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jun Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Jul,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jul Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Aug,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Aug Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Sept,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Sept Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Oct,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Oct Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Nov,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Nov Climatology",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents$Dec,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,15),las=1,main="Dec Climatology",xlab="% Suitable Habitat in State Waters")

sumjuststatewaters<-c(sum(State_water_percents$Jan),sum(State_water_percents$Feb),sum(State_water_percents$Mar),sum(State_water_percents$Apr),
                      sum(State_water_percents$May),sum(State_water_percents$Jun),sum(State_water_percents$Jul),sum(State_water_percents$Aug),
                      sum(State_water_percents$Sept),sum(State_water_percents$Oct),sum(State_water_percents$Nov),sum(State_water_percents$Dec))

#extreme years
ex201201_state_percents<-stateHabitatUse(raster_month=ex201201WRationrasterlayer)
ex201202_state_percents<-stateHabitatUse(raster_month=ex201202WRationrasterlayer)
ex201203_state_percents<-stateHabitatUse(raster_month=ex201203WRationrasterlayer)
ex201204_state_percents<-stateHabitatUse(raster_month=ex201204WRationrasterlayer)
ex201205_state_percents<-stateHabitatUse(raster_month=ex201205WRationrasterlayer)
ex201206_state_percents<-stateHabitatUse(raster_month=ex201206WRationrasterlayer)
ex201207_state_percents<-stateHabitatUse(raster_month=ex201207WRationrasterlayer)
ex201208_state_percents<-stateHabitatUse(raster_month=ex201208WRationrasterlayer)
ex201209_state_percents<-stateHabitatUse(raster_month=ex201209WRationrasterlayer)
ex201210_state_percents<-stateHabitatUse(raster_month=ex201210WRationrasterlayer)
ex201211_state_percents<-stateHabitatUse(raster_month=ex201211WRationrasterlayer)
ex201212_state_percents<-stateHabitatUse(raster_month=ex201212WRationrasterlayer)
ex199601_state_percents<-stateHabitatUse(raster_month=ex199601WRationrasterlayer)
ex199602_state_percents<-stateHabitatUse(raster_month=ex199602WRationrasterlayer)
ex199603_state_percents<-stateHabitatUse(raster_month=ex199603WRationrasterlayer)
ex199604_state_percents<-stateHabitatUse(raster_month=ex199604WRationrasterlayer)
ex199605_state_percents<-stateHabitatUse(raster_month=ex199605WRationrasterlayer)
ex199606_state_percents<-stateHabitatUse(raster_month=ex199606WRationrasterlayer)
ex199607_state_percents<-stateHabitatUse(raster_month=ex199607WRationrasterlayer)
ex199608_state_percents<-stateHabitatUse(raster_month=ex199608WRationrasterlayer)
ex199609_state_percents<-stateHabitatUse(raster_month=ex199609WRationrasterlayer)
ex199610_state_percents<-stateHabitatUse(raster_month=ex199610WRationrasterlayer)
ex199611_state_percents<-stateHabitatUse(raster_month=ex199611WRationrasterlayer)
ex199612_state_percents<-stateHabitatUse(raster_month=ex199612WRationrasterlayer)


State_water_percents2012ex<-data.frame(States=States,Jan=ex201201_state_percents[1:14]*100,Feb=ex201202_state_percents[1:14]*100,Mar=ex201203_state_percents[1:14]*100,
                                       Apr=ex201204_state_percents[1:14]*100,May=ex201205_state_percents[1:14]*100,Jun=ex201206_state_percents[1:14]*100,Jul=ex201207_state_percents[1:14]*100,
                                       Aug=ex201208_state_percents[1:14]*100,Sept=ex201209_state_percents[1:14]*100,Oct=ex201210_state_percents[1:14]*100,Nov=ex201211_state_percents[1:14]*100,Dec=ex201212_state_percents[1:14]*100)
revState_water_percents2012ex<-State_water_percents2012ex[nrow(State_water_percents2012ex):1,]

State_water_percents1996ex<-data.frame(States=States,Jan=ex199601_state_percents[1:14]*100,Feb=ex199602_state_percents[1:14]*100,Mar=ex199603_state_percents[1:14]*100,
                                       Apr=ex199604_state_percents[1:14]*100,May=ex199605_state_percents[1:14]*100,Jun=ex199606_state_percents[1:14]*100,Jul=ex199607_state_percents[1:14]*100,
                                       Aug=ex199608_state_percents[1:14]*100,Sept=ex199609_state_percents[1:14]*100,Oct=ex199610_state_percents[1:14]*100,Nov=ex199611_state_percents[1:14]*100,Dec=ex199612_state_percents[1:14]*100)
revState_water_percents1996ex<-State_water_percents1996ex[nrow(State_water_percents1996ex):1,]

par(mfrow=c(2,6))
barplot(height=revState_water_percents2012ex$Jan,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jan EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Feb,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Feb EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Mar,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Mar EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Apr,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Apr EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$May,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="May EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Jun,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jun EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Jul,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jul EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Aug,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Aug EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Sept,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Sept EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Oct,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Oct EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Nov,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Nov EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents2012ex$Dec,names.arg=revState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Dec EX 2012 (Warm)",xlab="% Suitable Habitat in State Waters")

par(mfrow=c(2,6))
barplot(height=revState_water_percents1996ex$Jan,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jan EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Feb,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Feb EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Mar,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Mar EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Apr,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Apr EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$May,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="May EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Jun,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jun EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Jul,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Jul EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Aug,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Aug EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Sept,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Sept EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Oct,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Oct EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Nov,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Nov EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")
barplot(height=revState_water_percents1996ex$Dec,names.arg=revState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,15),las=1,main="Dec EX 1996 (Cool)",xlab="% Suitable Habitat in State Waters")

sumjuststatewaters2012ex<-c(sum(State_water_percents2012ex$Jan),sum(State_water_percents2012ex$Feb),sum(State_water_percents2012ex$Mar),sum(State_water_percents2012ex$Apr),
                            sum(State_water_percents2012ex$May),sum(State_water_percents2012ex$Jun),sum(State_water_percents2012ex$Jul),sum(State_water_percents2012ex$Aug),
                            sum(State_water_percents2012ex$Sept),sum(State_water_percents2012ex$Oct),sum(State_water_percents2012ex$Nov),sum(State_water_percents2012ex$Dec))

sumjuststatewaters1996ex<-c(sum(State_water_percents1996ex$Jan),sum(State_water_percents1996ex$Feb),sum(State_water_percents1996ex$Mar),sum(State_water_percents1996ex$Apr),
                            sum(State_water_percents1996ex$May),sum(State_water_percents1996ex$Jun),sum(State_water_percents1996ex$Jul),sum(State_water_percents1996ex$Aug),
                            sum(State_water_percents1996ex$Sept),sum(State_water_percents1996ex$Oct),sum(State_water_percents1996ex$Nov),sum(State_water_percents1996ex$Dec))

#future climatologies
#0-20
Jan0_20_state_percents<-stateHabitatUse(raster_month=Jan0_20WRationrasterlayer)
Feb0_20_state_percents<-stateHabitatUse(raster_month=Feb0_20WRationrasterlayer)
Mar0_20_state_percents<-stateHabitatUse(raster_month=Mar0_20WRationrasterlayer)
Apr0_20_state_percents<-stateHabitatUse(raster_month=Apr0_20WRationrasterlayer)
May0_20_state_percents<-stateHabitatUse(raster_month=May0_20WRationrasterlayer)
Jun0_20_state_percents<-stateHabitatUse(raster_month=Jun0_20WRationrasterlayer)
Jul0_20_state_percents<-stateHabitatUse(raster_month=Jul0_20WRationrasterlayer)
Aug0_20_state_percents<-stateHabitatUse(raster_month=Aug0_20WRationrasterlayer)
Sept0_20_state_percents<-stateHabitatUse(raster_month=Sept0_20WRationrasterlayer)
Oct0_20_state_percents<-stateHabitatUse(raster_month=Oct0_20WRationrasterlayer)
Nov0_20_state_percents<-stateHabitatUse(raster_month=Nov0_20WRationrasterlayer)
Dec0_20_state_percents<-stateHabitatUse(raster_month=Dec0_20WRationrasterlayer)

State0_20_water_percents<-data.frame(States=States,Jan=Jan0_20_state_percents[1:14]*100,Feb=Feb0_20_state_percents[1:14]*100,Mar=Mar0_20_state_percents[1:14]*100,
                                     Apr=Apr0_20_state_percents[1:14]*100,May=May0_20_state_percents[1:14]*100,Jun=Jun0_20_state_percents[1:14]*100,Jul=Jul0_20_state_percents[1:14]*100,
                                     Aug=Aug0_20_state_percents[1:14]*100,Sept=Sept0_20_state_percents[1:14]*100,Oct=Oct0_20_state_percents[1:14]*100,Nov=Nov0_20_state_percents[1:14]*100,Dec=Dec0_20_state_percents[1:14]*100)
revState0_20_water_percents<-State0_20_water_percents[nrow(State0_20_water_percents):1,]
#20_40
Jan20_40_state_percents<-stateHabitatUse(raster_month=Jan20_40WRationrasterlayer)
Feb20_40_state_percents<-stateHabitatUse(raster_month=Feb20_40WRationrasterlayer)
Mar20_40_state_percents<-stateHabitatUse(raster_month=Mar20_40WRationrasterlayer)
Apr20_40_state_percents<-stateHabitatUse(raster_month=Apr20_40WRationrasterlayer)
May20_40_state_percents<-stateHabitatUse(raster_month=May20_40WRationrasterlayer)
Jun20_40_state_percents<-stateHabitatUse(raster_month=Jun20_40WRationrasterlayer)
Jul20_40_state_percents<-stateHabitatUse(raster_month=Jul20_40WRationrasterlayer)
Aug20_40_state_percents<-stateHabitatUse(raster_month=Aug20_40WRationrasterlayer)
Sept20_40_state_percents<-stateHabitatUse(raster_month=Sept20_40WRationrasterlayer)
Oct20_40_state_percents<-stateHabitatUse(raster_month=Oct20_40WRationrasterlayer)
Nov20_40_state_percents<-stateHabitatUse(raster_month=Nov20_40WRationrasterlayer)
Dec20_40_state_percents<-stateHabitatUse(raster_month=Dec20_40WRationrasterlayer)

State20_40_water_percents<-data.frame(States=States,Jan=Jan20_40_state_percents[1:14]*100,Feb=Feb20_40_state_percents[1:14]*100,Mar=Mar20_40_state_percents[1:14]*100,
                                      Apr=Apr20_40_state_percents[1:14]*100,May=May20_40_state_percents[1:14]*100,Jun=Jun20_40_state_percents[1:14]*100,Jul=Jul20_40_state_percents[1:14]*100,
                                      Aug=Aug20_40_state_percents[1:14]*100,Sept=Sept20_40_state_percents[1:14]*100,Oct=Oct20_40_state_percents[1:14]*100,Nov=Nov20_40_state_percents[1:14]*100,Dec=Dec20_40_state_percents[1:14]*100)
revState20_40_water_percents<-State20_40_water_percents[nrow(State20_40_water_percents):1,]
#40_60
Jan40_60_state_percents<-stateHabitatUse(raster_month=Jan40_60WRationrasterlayer)
Feb40_60_state_percents<-stateHabitatUse(raster_month=Feb40_60WRationrasterlayer)
Mar40_60_state_percents<-stateHabitatUse(raster_month=Mar40_60WRationrasterlayer)
Apr40_60_state_percents<-stateHabitatUse(raster_month=Apr40_60WRationrasterlayer)
May40_60_state_percents<-stateHabitatUse(raster_month=May40_60WRationrasterlayer)
Jun40_60_state_percents<-stateHabitatUse(raster_month=Jun40_60WRationrasterlayer)
Jul40_60_state_percents<-stateHabitatUse(raster_month=Jul40_60WRationrasterlayer)
Aug40_60_state_percents<-stateHabitatUse(raster_month=Aug40_60WRationrasterlayer)
Sept40_60_state_percents<-stateHabitatUse(raster_month=Sept40_60WRationrasterlayer)
Oct40_60_state_percents<-stateHabitatUse(raster_month=Oct40_60WRationrasterlayer)
Nov40_60_state_percents<-stateHabitatUse(raster_month=Nov40_60WRationrasterlayer)
Dec40_60_state_percents<-stateHabitatUse(raster_month=Dec40_60WRationrasterlayer)

State40_60_water_percents<-data.frame(States=States,Jan=Jan40_60_state_percents[1:14]*100,Feb=Feb40_60_state_percents[1:14]*100,Mar=Mar40_60_state_percents[1:14]*100,
                                      Apr=Apr40_60_state_percents[1:14]*100,May=May40_60_state_percents[1:14]*100,Jun=Jun40_60_state_percents[1:14]*100,Jul=Jul40_60_state_percents[1:14]*100,
                                      Aug=Aug40_60_state_percents[1:14]*100,Sept=Sept40_60_state_percents[1:14]*100,Oct=Oct40_60_state_percents[1:14]*100,Nov=Nov40_60_state_percents[1:14]*100,Dec=Dec40_60_state_percents[1:14]*100)
revState40_60_water_percents<-State40_60_water_percents[nrow(State40_60_water_percents):1,]
#60_80
Jan60_80_state_percents<-stateHabitatUse(raster_month=Jan60_80WRationrasterlayer)
Feb60_80_state_percents<-stateHabitatUse(raster_month=Feb60_80WRationrasterlayer)
Mar60_80_state_percents<-stateHabitatUse(raster_month=Mar60_80WRationrasterlayer)
Apr60_80_state_percents<-stateHabitatUse(raster_month=Apr60_80WRationrasterlayer)
May60_80_state_percents<-stateHabitatUse(raster_month=May60_80WRationrasterlayer)
Jun60_80_state_percents<-stateHabitatUse(raster_month=Jun60_80WRationrasterlayer)
Jul60_80_state_percents<-stateHabitatUse(raster_month=Jul60_80WRationrasterlayer)
Aug60_80_state_percents<-stateHabitatUse(raster_month=Aug60_80WRationrasterlayer)
Sept60_80_state_percents<-stateHabitatUse(raster_month=Sept60_80WRationrasterlayer)
Oct60_80_state_percents<-stateHabitatUse(raster_month=Oct60_80WRationrasterlayer)
Nov60_80_state_percents<-stateHabitatUse(raster_month=Nov60_80WRationrasterlayer)
Dec60_80_state_percents<-stateHabitatUse(raster_month=Dec60_80WRationrasterlayer)

State60_80_water_percents<-data.frame(States=States,Jan=Jan60_80_state_percents[1:14]*100,Feb=Feb60_80_state_percents[1:14]*100,Mar=Mar60_80_state_percents[1:14]*100,
                                      Apr=Apr60_80_state_percents[1:14]*100,May=May60_80_state_percents[1:14]*100,Jun=Jun60_80_state_percents[1:14]*100,Jul=Jul60_80_state_percents[1:14]*100,
                                      Aug=Aug60_80_state_percents[1:14]*100,Sept=Sept60_80_state_percents[1:14]*100,Oct=Oct60_80_state_percents[1:14]*100,Nov=Nov60_80_state_percents[1:14]*100,Dec=Dec60_80_state_percents[1:14]*100)
revState60_80_water_percents<-State60_80_water_percents[nrow(State60_80_water_percents):1,]

#ExtremeFutureYear
Jan70_state_percents<-stateHabitatUse(raster_month=Jan70WRationrasterlayer)
Feb70_state_percents<-stateHabitatUse(raster_month=Feb70WRationrasterlayer)
Mar70_state_percents<-stateHabitatUse(raster_month=Mar70WRationrasterlayer)
Apr70_state_percents<-stateHabitatUse(raster_month=Apr70WRationrasterlayer)
May70_state_percents<-stateHabitatUse(raster_month=May70WRationrasterlayer)
Jun70_state_percents<-stateHabitatUse(raster_month=Jun70WRationrasterlayer)
Jul70_state_percents<-stateHabitatUse(raster_month=Jul70WRationrasterlayer)
Aug70_state_percents<-stateHabitatUse(raster_month=Aug70WRationrasterlayer)
Sept70_state_percents<-stateHabitatUse(raster_month=Sept70WRationrasterlayer)
Oct70_state_percents<-stateHabitatUse(raster_month=Oct70WRationrasterlayer)
Nov70_state_percents<-stateHabitatUse(raster_month=Nov70WRationrasterlayer)
Dec70_state_percents<-stateHabitatUse(raster_month=Dec70WRationrasterlayer)

State70_water_percents<-data.frame(States=States,Jan=Jan70_state_percents[1:14]*100,Feb=Feb70_state_percents[1:14]*100,Mar=Mar70_state_percents[1:14]*100,
                                   Apr=Apr70_state_percents[1:14]*100,May=May70_state_percents[1:14]*100,Jun=Jun70_state_percents[1:14]*100,Jul=Jul70_state_percents[1:14]*100,
                                   Aug=Aug70_state_percents[1:14]*100,Sept=Sept70_state_percents[1:14]*100,Oct=Oct70_state_percents[1:14]*100,Nov=Nov70_state_percents[1:14]*100,Dec=Dec70_state_percents[1:14]*100)
revState70_water_percents<-State70_water_percents[nrow(State70_water_percents):1,]


par(mfrow=c(2,6))
par(mar=c(3,3,4,1))
white=rgb(1,1,1,alpha=0)
barplot(height=revState0_20_water_percents$Jan,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="January",xlab="",col=white,border="blue")
barplot(height=revState20_40_water_percents$Jan,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Jan,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Jan,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Jan,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
legend(2,16,legend=c("Climatology","0-20","20-40","40-60","60-80"),pch=c(22,22,22,22,22),col=c("black","blue","green","orange","red"),cex=1.5)
barplot(height=revState0_20_water_percents$Feb,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="February",xlab="",col=white,border="blue")
barplot(height=revState20_40_water_percents$Feb,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Feb,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Feb,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Feb,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$Mar,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="March",xlab="",col=white,border="blue")
barplot(height=revState20_40_water_percents$Mar,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Mar,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Mar,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Mar,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$Apr,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="April",xlab="",col=white,border="blue")
barplot(height=revState20_40_water_percents$Apr,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Apr,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Apr,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Apr,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$May,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="May",xlab="",col=white,border="blue")
barplot(height=revState20_40_water_percents$May,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$May,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$May,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$May,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$Jun,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="June",xlab="",col=white,border="blue")
barplot(height=revState20_40_water_percents$Jun,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Jun,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Jun,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Jun,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
par(mar=c(5,3,2,1))
barplot(height=revState0_20_water_percents$Jul,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="July",xlab="% SH (In State Waters)",col=white,border="blue")
barplot(height=revState20_40_water_percents$Jul,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Jul,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Jul,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Jul,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$Aug,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="August",xlab="% SH (In State Waters)",col=white,border="blue")
barplot(height=revState20_40_water_percents$Aug,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Aug,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Aug,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Aug,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$Sept,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="September",xlab="% SH (In State Waters)",col=white,border="blue")
barplot(height=revState20_40_water_percents$Sept,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Sept,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Sept,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Sept,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$Oct,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="October",xlab="% SH (In State Waters)",col=white,border="blue")
barplot(height=revState20_40_water_percents$Oct,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Oct,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Oct,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Oct,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$Nov,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="November",xlab="% SH (In State Waters)",col=white,border="blue")
barplot(height=revState20_40_water_percents$Nov,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Nov,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Nov,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Nov,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revState0_20_water_percents$Dec,names.arg=revState_water_percents$States,horiz=TRUE,xlim=c(0,18),las=1,main="December",xlab="% SH (In State Waters)",col=white,border="blue")
barplot(height=revState20_40_water_percents$Dec,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revState40_60_water_percents$Dec,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revState60_80_water_percents$Dec,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revState_water_percents$Dec,names.arg=revState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")



##################
#All water offshore use
##################
ME_NH_line<-42.8723
MA_RI_line<-41.3043
CT_NY_line<-40.4897
NJ_line<-38.8801
DE_line<-38.4512
MD_line_Ocean<-38.0266
MD_line_Bay<-37.8869
MD_line<-mean(c(MD_line_Ocean,MD_line_Bay))
VA_line<-36.5508
NC_line<-33.8106
SC_line<-32.0333
GA_line<-30.7122

#this function crops offshore state waters from the overall large raster of cobia habitat use over shelf
#then calculates habitat value for each offshore state waters
offshoreraster<-function(state_line1=NA,state_line2=NA,raster_month=NA){
  new_lats<-lat[which(lat<state_line1 & lat>=state_line2)]
  lon_newlats<-expand.grid(lon=lon, lat=new_lats)
  coordinates(lon_newlats)<- ~lon + lat 
  proj4string(lon_newlats)<-CRS("+init=epsg:4326")#sets it to lat-long
  e_newlats<-extent(c(min(lon),max(lon),min(new_lats),max(new_lats)))
  
  offshorecrop<-crop(raster_month,e_newlats,snap="out") 
  offshorearearaster<-area(offshorecrop,na.rm=TRUE) #area is in km2
  offshorearearasteriso<-offshorearearaster[!is.na(offshorearearaster)]#only isolates (gives vector of) cell area values without NAs
  offshorerastervals<-values(offshorecrop) #values of all cells
  offshorerastervalsiso<-offshorerastervals[!is.na(offshorerastervals)] #only isolates (gives vector of) cell values (ratio value) without NAs
  offshoretotalval<-offshorearearasteriso*offshorerastervalsiso #gives area value x ratio value
  offshoresumval<-sum(offshoretotalval) #total sum of all "values" on shelf
  
  return(offshoresumval=offshoresumval)
}

#this function calculate overall habitat value along US shelf and divides it by habitat value for each offshore state waters to get a percentage
stateOffshoreHabitatUse<-function(raster_month=NA){
  raster_month[raster_month <= 1] <- NA
  arearaster<-area(raster_month,na.rm=TRUE) #area is in km2
  arearasteriso<-arearaster[!is.na(arearaster)]#only isolates (gives vector of) cell area values without NAs
  rastervals<-values(raster_month) #values of all cells
  rastervalsiso<-rastervals[!is.na(rastervals)] #only isolates (gives vector of) cell values (ratio value) without NAs
  totalval<-arearasteriso*rastervalsiso #gives area value x ratio value
  sumval<-sum(totalval) #total sum of all "values" on shelf
  
  ME_NHoffshoreumval<-offshoreraster(state_line1=max(lat),state_line2=ME_NH_line,raster_month=raster_month)
  MA_RIoffshoreumval<-offshoreraster(state_line1=ME_NH_line,state_line2=MA_RI_line,raster_month=raster_month)
  CT_NYoffshoreumval<-offshoreraster(state_line1=MA_RI_line,state_line2=CT_NY_line,raster_month=raster_month)
  NJoffshoreumval<-offshoreraster(state_line1=CT_NY_line,state_line2=NJ_line,raster_month=raster_month)
  DEoffshoreumval<-offshoreraster(state_line1=NJ_line,state_line2=DE_line,raster_month=raster_month)
  MDoffshoreumval<-offshoreraster(state_line1=DE_line,state_line2=MD_line,raster_month=raster_month)
  VAoffshoreumval<-offshoreraster(state_line1=MD_line,state_line2=VA_line,raster_month=raster_month)
  NCoffshoreumval<-offshoreraster(state_line1=VA_line,state_line2=NC_line,raster_month=raster_month)
  SCoffshoreumval<-offshoreraster(state_line1=NC_line,state_line2=SC_line,raster_month=raster_month)
  GAoffshoreumval<-offshoreraster(state_line1=SC_line,state_line2=GA_line,raster_month=raster_month)
  FLoffshoreumval<-offshoreraster(state_line1=GA_line,state_line2=min(lat),raster_month=raster_month)
  
  ME_NHoffshore_percent<-ME_NHoffshoreumval/sumval
  MA_RIoffshore_percent<-MA_RIoffshoreumval/sumval
  CT_NYoffshore_percent<-CT_NYoffshoreumval/sumval
  NJoffshore_percent<-NJoffshoreumval/sumval
  DEoffshore_percent<-DEoffshoreumval/sumval
  MDoffshore_percent<-MDoffshoreumval/sumval
  VAoffshore_percent<-VAoffshoreumval/sumval
  NCoffshore_percent<-NCoffshoreumval/sumval
  SCoffshore_percent<-SCoffshoreumval/sumval
  GAoffshore_percent<-GAoffshoreumval/sumval
  FLoffshore_percent<-FLoffshoreumval/sumval
  
  
  offshore_percents<-c(ME_NHoffshore_percent,MA_RIoffshore_percent,CT_NYoffshore_percent,NJoffshore_percent,
                       DEoffshore_percent,MDoffshore_percent,VAoffshore_percent,NCoffshore_percent,SCoffshore_percent,
                       GAoffshore_percent,FLoffshore_percent,sumval)
  return(offshore_percents=offshore_percents)
  
}

Offshore_States<-c("ME_NH","MA_RI","CT_NY","NJ","DE","MD","VA","NC","SC","GA","FL")
#climatology
Jan_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=JanWRationrasterlayer)
Feb_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=FebWRationrasterlayer)
Mar_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=MarWRationrasterlayer)
Apr_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=AprWRationrasterlayer)
May_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=MayWRationrasterlayer)
Jun_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=JunWRationrasterlayer)
Jul_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=JulWRationrasterlayer)
Aug_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=AugWRationrasterlayer)
Sept_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=SeptWRationrasterlayer)
Oct_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=OctWRationrasterlayer)
Nov_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=NovWRationrasterlayer)
Dec_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=DecWRationrasterlayer)


OffshoreState_water_percents<-data.frame(States=Offshore_States,Jan=Jan_offshorestate_percents[1:11]*100,Feb=Feb_offshorestate_percents[1:11]*100,Mar=Mar_offshorestate_percents[1:11]*100,
                                         Apr=Apr_offshorestate_percents[1:11]*100,May=May_offshorestate_percents[1:11]*100,Jun=Jun_offshorestate_percents[1:11]*100,Jul=Jul_offshorestate_percents[1:11]*100,
                                         Aug=Aug_offshorestate_percents[1:11]*100,Sept=Sept_offshorestate_percents[1:11]*100,Oct=Oct_offshorestate_percents[1:11]*100,Nov=Nov_offshorestate_percents[1:11]*100,Dec=Dec_offshorestate_percents[1:11]*100)
revOffshoreState_water_percents<-OffshoreState_water_percents[nrow(OffshoreState_water_percents):1,]

par(mfrow=c(2,6))
barplot(height=revOffshoreState_water_percents$Jan,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Jan Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Feb,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Feb Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Mar,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Mar Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Apr,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Apr Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$May,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="May Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Jun,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Jun Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Jul,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Jul Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Aug,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Aug Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Sept,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Sept Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Oct,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Oct Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Nov,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Nov Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents$Dec,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,50),main="Dec Climatology",las=1,xlab="% Suit Hab in Offshore State Waters")

#sum should add to 100
sumjustoffshorestatewaters<-c(sum(OffshoreState_water_percents$Jan),sum(OffshoreState_water_percents$Feb),sum(OffshoreState_water_percents$Mar),sum(OffshoreState_water_percents$Apr),
                              sum(OffshoreState_water_percents$May),sum(OffshoreState_water_percents$Jun),sum(OffshoreState_water_percents$Jul),sum(OffshoreState_water_percents$Aug),
                              sum(OffshoreState_water_percents$Sept),sum(OffshoreState_water_percents$Oct),sum(OffshoreState_water_percents$Nov),sum(OffshoreState_water_percents$Dec))

#extreme years
ex201201_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201201WRationrasterlayer)
ex201202_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201202WRationrasterlayer)
ex201203_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201203WRationrasterlayer)
ex201204_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201204WRationrasterlayer)
ex201205_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201205WRationrasterlayer)
ex201206_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201206WRationrasterlayer)
ex201207_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201207WRationrasterlayer)
ex201208_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201208WRationrasterlayer)
ex201209_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201209WRationrasterlayer)
ex201210_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201210WRationrasterlayer)
ex201211_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201211WRationrasterlayer)
ex201212_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex201212WRationrasterlayer)
ex199601_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199601WRationrasterlayer)
ex199602_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199602WRationrasterlayer)
ex199603_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199603WRationrasterlayer)
ex199604_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199604WRationrasterlayer)
ex199605_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199605WRationrasterlayer)
ex199606_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199606WRationrasterlayer)
ex199607_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199607WRationrasterlayer)
ex199608_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199608WRationrasterlayer)
ex199609_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199609WRationrasterlayer)
ex199610_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199610WRationrasterlayer)
ex199611_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199611WRationrasterlayer)
ex199612_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=ex199612WRationrasterlayer)


OffshoreState_water_percents2012ex<-data.frame(States=Offshore_States,Jan=ex201201_offshorestate_percents[1:11]*100,Feb=ex201202_offshorestate_percents[1:11]*100,Mar=ex201203_offshorestate_percents[1:11]*100,
                                               Apr=ex201204_offshorestate_percents[1:11]*100,May=ex201205_offshorestate_percents[1:11]*100,Jun=ex201206_offshorestate_percents[1:11]*100,Jul=ex201207_offshorestate_percents[1:11]*100,
                                               Aug=ex201208_offshorestate_percents[1:11]*100,Sept=ex201209_offshorestate_percents[1:11]*100,Oct=ex201210_offshorestate_percents[1:11]*100,Nov=ex201211_offshorestate_percents[1:11]*100,Dec=ex201212_offshorestate_percents[1:11]*100)
revOffshoreState_water_percents2012ex<-OffshoreState_water_percents2012ex[nrow(OffshoreState_water_percents2012ex):1,]

OffshoreState_water_percents1996ex<-data.frame(States=Offshore_States,Jan=ex199601_offshorestate_percents[1:11]*100,Feb=ex199602_offshorestate_percents[1:11]*100,Mar=ex199603_offshorestate_percents[1:11]*100,
                                               Apr=ex199604_offshorestate_percents[1:11]*100,May=ex199605_offshorestate_percents[1:11]*100,Jun=ex199606_offshorestate_percents[1:11]*100,Jul=ex199607_offshorestate_percents[1:11]*100,
                                               Aug=ex199608_offshorestate_percents[1:11]*100,Sept=ex199609_offshorestate_percents[1:11]*100,Oct=ex199610_offshorestate_percents[1:11]*100,Nov=ex199611_offshorestate_percents[1:11]*100,Dec=ex199612_offshorestate_percents[1:11]*100)
revOffshoreState_water_percents1996ex<-OffshoreState_water_percents1996ex[nrow(OffshoreState_water_percents1996ex):1,]

par(mfrow=c(2,6))
barplot(height=revOffshoreState_water_percents2012ex$Jan,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Jan EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Feb,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Feb EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Mar,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Mar EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Apr,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Apr EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$May,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="May EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Jun,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Jun EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Jul,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Jul EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Aug,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Aug EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Sept,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Sept EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Oct,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Oct EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Nov,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Nov EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents2012ex$Dec,names.arg=revOffshoreState_water_percents2012ex$States,horiz=TRUE,xlim=c(0,50),main="Dec EX 2012 (Warm)",las=1,xlab="% Suit Hab in Offshore State Waters")

barplot(height=revOffshoreState_water_percents1996ex$Jan,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Jan EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Feb,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Feb EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Mar,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Mar EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Apr,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Apr EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$May,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="May EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Jun,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Jun EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Jul,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Jul EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Aug,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Aug EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Sept,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Sept EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Oct,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Oct EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Nov,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Nov EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")
barplot(height=revOffshoreState_water_percents1996ex$Dec,names.arg=revOffshoreState_water_percents1996ex$States,horiz=TRUE,xlim=c(0,50),main="Dec EX 1996 (Cool)",las=1,xlab="% Suit Hab in Offshore State Waters")

#sum should add to 100
sumjustoffshorestatewaters2012ex<-c(sum(OffshoreState_water_percents2012ex$Jan),sum(OffshoreState_water_percents2012ex$Feb),sum(OffshoreState_water_percents2012ex$Mar),sum(OffshoreState_water_percents2012ex$Apr),
                                    sum(OffshoreState_water_percents2012ex$May),sum(OffshoreState_water_percents2012ex$Jun),sum(OffshoreState_water_percents2012ex$Jul),sum(OffshoreState_water_percents2012ex$Aug),
                                    sum(OffshoreState_water_percents2012ex$Sept),sum(OffshoreState_water_percents2012ex$Oct),sum(OffshoreState_water_percents2012ex$Nov),sum(OffshoreState_water_percents2012ex$Dec))
#sum should add to 100
sumjustoffshorestatewaters1996ex<-c(sum(OffshoreState_water_percents1996ex$Jan),sum(OffshoreState_water_percents1996ex$Feb),sum(OffshoreState_water_percents1996ex$Mar),sum(OffshoreState_water_percents1996ex$Apr),
                                    sum(OffshoreState_water_percents1996ex$May),sum(OffshoreState_water_percents1996ex$Jun),sum(OffshoreState_water_percents1996ex$Jul),sum(OffshoreState_water_percents1996ex$Aug),
                                    sum(OffshoreState_water_percents1996ex$Sept),sum(OffshoreState_water_percents1996ex$Oct),sum(OffshoreState_water_percents1996ex$Nov),sum(OffshoreState_water_percents1996ex$Dec))


#future climatologies
#0_20
Jan0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jan0_20WRationrasterlayer)
Feb0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Feb0_20WRationrasterlayer)
Mar0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Mar0_20WRationrasterlayer)
Apr0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Apr0_20WRationrasterlayer)
May0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=May0_20WRationrasterlayer)
Jun0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jun0_20WRationrasterlayer)
Jul0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jul0_20WRationrasterlayer)
Aug0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Aug0_20WRationrasterlayer)
Sept0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Sept0_20WRationrasterlayer)
Oct0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Oct0_20WRationrasterlayer)
Nov0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Nov0_20WRationrasterlayer)
Dec0_20_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Dec0_20WRationrasterlayer)

OffshoreState0_20_water_percents<-data.frame(States=Offshore_States,Jan=Jan0_20_offshorestate_percents[1:11]*100,Feb=Feb0_20_offshorestate_percents[1:11]*100,Mar=Mar0_20_offshorestate_percents[1:11]*100,
                                             Apr=Apr0_20_offshorestate_percents[1:11]*100,May=May0_20_offshorestate_percents[1:11]*100,Jun=Jun0_20_offshorestate_percents[1:11]*100,Jul=Jul0_20_offshorestate_percents[1:11]*100,
                                             Aug=Aug0_20_offshorestate_percents[1:11]*100,Sept=Sept0_20_offshorestate_percents[1:11]*100,Oct=Oct0_20_offshorestate_percents[1:11]*100,Nov=Nov0_20_offshorestate_percents[1:11]*100,Dec=Dec0_20_offshorestate_percents[1:11]*100)
revOffshoreState0_20_water_percents<-OffshoreState0_20_water_percents[nrow(OffshoreState0_20_water_percents):1,]
#sum should add to 100
sumjustoffshorestate0_20waters<-c(sum(OffshoreState0_20_water_percents$Jan),sum(OffshoreState0_20_water_percents$Feb),sum(OffshoreState0_20_water_percents$Mar),sum(OffshoreState0_20_water_percents$Apr),
                                  sum(OffshoreState0_20_water_percents$May),sum(OffshoreState0_20_water_percents$Jun),sum(OffshoreState0_20_water_percents$Jul),sum(OffshoreState0_20_water_percents$Aug),
                                  sum(OffshoreState0_20_water_percents$Sept),sum(OffshoreState0_20_water_percents$Oct),sum(OffshoreState0_20_water_percents$Nov),sum(OffshoreState0_20_water_percents$Dec))
#20_40
Jan20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jan20_40WRationrasterlayer)
Feb20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Feb20_40WRationrasterlayer)
Mar20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Mar20_40WRationrasterlayer)
Apr20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Apr20_40WRationrasterlayer)
May20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=May20_40WRationrasterlayer)
Jun20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jun20_40WRationrasterlayer)
Jul20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jul20_40WRationrasterlayer)
Aug20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Aug20_40WRationrasterlayer)
Sept20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Sept20_40WRationrasterlayer)
Oct20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Oct20_40WRationrasterlayer)
Nov20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Nov20_40WRationrasterlayer)
Dec20_40_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Dec20_40WRationrasterlayer)

OffshoreState20_40_water_percents<-data.frame(States=Offshore_States,Jan=Jan20_40_offshorestate_percents[1:11]*100,Feb=Feb20_40_offshorestate_percents[1:11]*100,Mar=Mar20_40_offshorestate_percents[1:11]*100,
                                              Apr=Apr20_40_offshorestate_percents[1:11]*100,May=May20_40_offshorestate_percents[1:11]*100,Jun=Jun20_40_offshorestate_percents[1:11]*100,Jul=Jul20_40_offshorestate_percents[1:11]*100,
                                              Aug=Aug20_40_offshorestate_percents[1:11]*100,Sept=Sept20_40_offshorestate_percents[1:11]*100,Oct=Oct20_40_offshorestate_percents[1:11]*100,Nov=Nov20_40_offshorestate_percents[1:11]*100,Dec=Dec20_40_offshorestate_percents[1:11]*100)
revOffshoreState20_40_water_percents<-OffshoreState20_40_water_percents[nrow(OffshoreState20_40_water_percents):1,]
#sum should add to 100
sumjustoffshorestate20_40waters<-c(sum(OffshoreState20_40_water_percents$Jan),sum(OffshoreState20_40_water_percents$Feb),sum(OffshoreState20_40_water_percents$Mar),sum(OffshoreState20_40_water_percents$Apr),
                                   sum(OffshoreState20_40_water_percents$May),sum(OffshoreState20_40_water_percents$Jun),sum(OffshoreState20_40_water_percents$Jul),sum(OffshoreState20_40_water_percents$Aug),
                                   sum(OffshoreState20_40_water_percents$Sept),sum(OffshoreState20_40_water_percents$Oct),sum(OffshoreState20_40_water_percents$Nov),sum(OffshoreState20_40_water_percents$Dec))
#40_60
Jan40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jan40_60WRationrasterlayer)
Feb40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Feb40_60WRationrasterlayer)
Mar40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Mar40_60WRationrasterlayer)
Apr40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Apr40_60WRationrasterlayer)
May40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=May40_60WRationrasterlayer)
Jun40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jun40_60WRationrasterlayer)
Jul40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jul40_60WRationrasterlayer)
Aug40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Aug40_60WRationrasterlayer)
Sept40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Sept40_60WRationrasterlayer)
Oct40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Oct40_60WRationrasterlayer)
Nov40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Nov40_60WRationrasterlayer)
Dec40_60_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Dec40_60WRationrasterlayer)

OffshoreState40_60_water_percents<-data.frame(States=Offshore_States,Jan=Jan40_60_offshorestate_percents[1:11]*100,Feb=Feb40_60_offshorestate_percents[1:11]*100,Mar=Mar40_60_offshorestate_percents[1:11]*100,
                                              Apr=Apr40_60_offshorestate_percents[1:11]*100,May=May40_60_offshorestate_percents[1:11]*100,Jun=Jun40_60_offshorestate_percents[1:11]*100,Jul=Jul40_60_offshorestate_percents[1:11]*100,
                                              Aug=Aug40_60_offshorestate_percents[1:11]*100,Sept=Sept40_60_offshorestate_percents[1:11]*100,Oct=Oct40_60_offshorestate_percents[1:11]*100,Nov=Nov40_60_offshorestate_percents[1:11]*100,Dec=Dec40_60_offshorestate_percents[1:11]*100)
revOffshoreState40_60_water_percents<-OffshoreState40_60_water_percents[nrow(OffshoreState40_60_water_percents):1,]
#sum should add to 100
sumjustoffshorestate40_60waters<-c(sum(OffshoreState40_60_water_percents$Jan),sum(OffshoreState40_60_water_percents$Feb),sum(OffshoreState40_60_water_percents$Mar),sum(OffshoreState40_60_water_percents$Apr),
                                   sum(OffshoreState40_60_water_percents$May),sum(OffshoreState40_60_water_percents$Jun),sum(OffshoreState40_60_water_percents$Jul),sum(OffshoreState40_60_water_percents$Aug),
                                   sum(OffshoreState40_60_water_percents$Sept),sum(OffshoreState40_60_water_percents$Oct),sum(OffshoreState40_60_water_percents$Nov),sum(OffshoreState40_60_water_percents$Dec))
#60_80
Jan60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jan60_80WRationrasterlayer)
Feb60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Feb60_80WRationrasterlayer)
Mar60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Mar60_80WRationrasterlayer)
Apr60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Apr60_80WRationrasterlayer)
May60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=May60_80WRationrasterlayer)
Jun60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jun60_80WRationrasterlayer)
Jul60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jul60_80WRationrasterlayer)
Aug60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Aug60_80WRationrasterlayer)
Sept60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Sept60_80WRationrasterlayer)
Oct60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Oct60_80WRationrasterlayer)
Nov60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Nov60_80WRationrasterlayer)
Dec60_80_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Dec60_80WRationrasterlayer)

OffshoreState60_80_water_percents<-data.frame(States=Offshore_States,Jan=Jan60_80_offshorestate_percents[1:11]*100,Feb=Feb60_80_offshorestate_percents[1:11]*100,Mar=Mar60_80_offshorestate_percents[1:11]*100,
                                              Apr=Apr60_80_offshorestate_percents[1:11]*100,May=May60_80_offshorestate_percents[1:11]*100,Jun=Jun60_80_offshorestate_percents[1:11]*100,Jul=Jul60_80_offshorestate_percents[1:11]*100,
                                              Aug=Aug60_80_offshorestate_percents[1:11]*100,Sept=Sept60_80_offshorestate_percents[1:11]*100,Oct=Oct60_80_offshorestate_percents[1:11]*100,Nov=Nov60_80_offshorestate_percents[1:11]*100,Dec=Dec60_80_offshorestate_percents[1:11]*100)
revOffshoreState60_80_water_percents<-OffshoreState60_80_water_percents[nrow(OffshoreState60_80_water_percents):1,]
#sum should add to 100
sumjustoffshorestate60_80waters<-c(sum(OffshoreState60_80_water_percents$Jan),sum(OffshoreState60_80_water_percents$Feb),sum(OffshoreState60_80_water_percents$Mar),sum(OffshoreState60_80_water_percents$Apr),
                                   sum(OffshoreState60_80_water_percents$May),sum(OffshoreState60_80_water_percents$Jun),sum(OffshoreState60_80_water_percents$Jul),sum(OffshoreState60_80_water_percents$Aug),
                                   sum(OffshoreState60_80_water_percents$Sept),sum(OffshoreState60_80_water_percents$Oct),sum(OffshoreState60_80_water_percents$Nov),sum(OffshoreState60_80_water_percents$Dec))


#FutureExtremeYears
Jan70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jan70WRationrasterlayer)
Feb70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Feb70WRationrasterlayer)
Mar70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Mar70WRationrasterlayer)
Apr70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Apr70WRationrasterlayer)
May70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=May70WRationrasterlayer)
Jun70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jun70WRationrasterlayer)
Jul70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Jul70WRationrasterlayer)
Aug70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Aug70WRationrasterlayer)
Sept70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Sept70WRationrasterlayer)
Oct70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Oct70WRationrasterlayer)
Nov70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Nov70WRationrasterlayer)
Dec70_offshorestate_percents<-stateOffshoreHabitatUse(raster_month=Dec70WRationrasterlayer)

OffshoreState70_water_percents<-data.frame(States=Offshore_States,Jan=Jan70_offshorestate_percents[1:11]*100,Feb=Feb70_offshorestate_percents[1:11]*100,Mar=Mar70_offshorestate_percents[1:11]*100,
                                           Apr=Apr70_offshorestate_percents[1:11]*100,May=May70_offshorestate_percents[1:11]*100,Jun=Jun70_offshorestate_percents[1:11]*100,Jul=Jul70_offshorestate_percents[1:11]*100,
                                           Aug=Aug70_offshorestate_percents[1:11]*100,Sept=Sept70_offshorestate_percents[1:11]*100,Oct=Oct70_offshorestate_percents[1:11]*100,Nov=Nov70_offshorestate_percents[1:11]*100,Dec=Dec70_offshorestate_percents[1:11]*100)
revOffshoreState70_water_percents<-OffshoreState70_water_percents[nrow(OffshoreState70_water_percents):1,]
#sum should add to 100
sumjustoffshorestate70waters<-c(sum(OffshoreState70_water_percents$Jan),sum(OffshoreState70_water_percents$Feb),sum(OffshoreState70_water_percents$Mar),sum(OffshoreState70_water_percents$Apr),
                                sum(OffshoreState70_water_percents$May),sum(OffshoreState70_water_percents$Jun),sum(OffshoreState70_water_percents$Jul),sum(OffshoreState70_water_percents$Aug),
                                sum(OffshoreState70_water_percents$Sept),sum(OffshoreState70_water_percents$Oct),sum(OffshoreState70_water_percents$Nov),sum(OffshoreState70_water_percents$Dec))


#keep in mind all of these plots are showing relative distributions, so % suitable habitat within that given year, doesn't take in account actual amount of suitable habitat that may have changed from one year to another
par(mfrow=c(2,6))
par(mar=c(3,4,4,0.5))
white=rgb(1,1,1,alpha=0)
barplot(height=revOffshoreState0_20_water_percents$Jan,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="January",xlab="",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Jan,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Jan,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Jan,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Jan,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
legend(5,13,legend=c("Climatology","0-20","20-40","40-60","60-80"),pch=c(22,22,22,22,22),col=c("black","blue","green","orange","red"),cex=1.5)
barplot(height=revOffshoreState0_20_water_percents$Feb,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="February",xlab="",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Feb,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Feb,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Feb,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Feb,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$Mar,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="March",xlab="",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Mar,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Mar,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Mar,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Mar,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$Apr,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="April",xlab="",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Apr,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Apr,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Apr,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Apr,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$May,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="May",xlab="",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$May,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$May,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$May,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$May,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$Jun,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="June",xlab="",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Jun,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Jun,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Jun,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Jun,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
par(mar=c(5,4,2,0.5))
barplot(height=revOffshoreState0_20_water_percents$Jul,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="July",xlab="% SH (All Waters)",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Jul,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Jul,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Jul,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Jul,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$Aug,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="August",xlab="% SH (All Waters)",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Aug,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Aug,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Aug,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Aug,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$Sept,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="September",xlab="% SH (All Waters)",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Sept,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Sept,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Sept,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Sept,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$Oct,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="October",xlab="% SH (All Waters)",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Oct,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Oct,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Oct,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Oct,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$Nov,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="November",xlab="% SH (All Waters)",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Nov,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Nov,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Nov,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Nov,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")
barplot(height=revOffshoreState0_20_water_percents$Dec,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlim=c(0,63),las=1,main="December",xlab="% SH (All Waters)",col=white,border="blue")
barplot(height=revOffshoreState20_40_water_percents$Dec,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="green",yaxt="n",xaxt="n")
barplot(height=revOffshoreState40_60_water_percents$Dec,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="orange",yaxt="n",xaxt="n")
barplot(height=revOffshoreState60_80_water_percents$Dec,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="red",yaxt="n",xaxt="n")
barplot(height=revOffshoreState_water_percents$Dec,names.arg=revOffshoreState_water_percents$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="black",yaxt="n",xaxt="n")





##############
#Calculate Change of state and offshore use between climatology (avg) and the extremes
#Calculate sum habitat suitability
##############
#inshore state waters
inshoredifchange2012<-data.frame(States=States)
for(i in 2:ncol(State_water_percents)){
  dif<-State_water_percents2012ex[,i]-State_water_percents[,i]
  inshoredifchange2012[,i]<-dif
}
colnames(inshoredifchange2012)<-c("States","Jan","Feb","Mar","Apr","May","Jun",
                                  "Jul","Aug","Sept","Oct","Nov","Dec")

inshoredifchange1996<-data.frame(States=States)
for(i in 2:ncol(State_water_percents)){
  dif<-State_water_percents1996ex[,i]-State_water_percents[,i]
  inshoredifchange1996[,i]<-dif
}
colnames(inshoredifchange1996)<-c("States","Jan","Feb","Mar","Apr","May","Jun",
                                  "Jul","Aug","Sept","Oct","Nov","Dec")

inshoredifchange70<-data.frame(States=States)
for(i in 2:ncol(State_water_percents)){
  dif<-State70_water_percents[,i]-State_water_percents[,i]
  inshoredifchange70[,i]<-dif
}
colnames(inshoredifchange70)<-c("States","Jan","Feb","Mar","Apr","May","Jun",
                                "Jul","Aug","Sept","Oct","Nov","Dec")

#offshore state waters
offshoredifchange2012<-data.frame(States=Offshore_States)
for(i in 2:ncol(OffshoreState_water_percents)){
  dif<-OffshoreState_water_percents2012ex[,i]-OffshoreState_water_percents[,i]
  offshoredifchange2012[,i]<-dif
}
colnames(offshoredifchange2012)<-c("States","Jan","Feb","Mar","Apr","May","Jun",
                                   "Jul","Aug","Sept","Oct","Nov","Dec")

offshoredifchange1996<-data.frame(States=Offshore_States)
for(i in 2:ncol(State_water_percents)){
  dif<-OffshoreState_water_percents1996ex[,i]-OffshoreState_water_percents[,i]
  offshoredifchange1996[,i]<-dif
}
colnames(offshoredifchange1996)<-c("States","Jan","Feb","Mar","Apr","May","Jun",
                                   "Jul","Aug","Sept","Oct","Nov","Dec")

offshoredifchange70<-data.frame(States=Offshore_States)
for(i in 2:ncol(State_water_percents)){
  dif<-OffshoreState70_water_percents[,i]-OffshoreState_water_percents[,i]
  offshoredifchange70[,i]<-dif
}
colnames(offshoredifchange70)<-c("States","Jan","Feb","Mar","Apr","May","Jun",
                                 "Jul","Aug","Sept","Oct","Nov","Dec")

revinshoredifchange2012<-inshoredifchange2012[nrow(inshoredifchange2012):1,]
revinshoredifchange1996<-inshoredifchange1996[nrow(inshoredifchange1996):1,]
revinshoredifchange70<-inshoredifchange70[nrow(inshoredifchange70):1,]
revoffshoredifchange2012<-offshoredifchange2012[nrow(offshoredifchange2012):1,]
revoffshoredifchange1996<-offshoredifchange1996[nrow(offshoredifchange1996):1,]
revoffshoredifchange70<-offshoredifchange70[nrow(offshoredifchange70):1,]

#inshore 2012 diff
par(mfrow=c(2,6))
barplot(height=revinshoredifchange2012$Jan,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Jan diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Feb,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Feb diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Mar,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Mar diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Apr,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Apr diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$May,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="May diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Jun,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Jun diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Jul,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Jul diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Aug,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Aug diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Sept,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Sept diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Oct,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Oct diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Nov,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Nov diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange2012$Dec,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),main="Dec diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in Inshore State Waters")

#inshore 1996 diff
par(mfrow=c(2,6))
barplot(height=revinshoredifchange1996$Jan,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Jan diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Feb,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Feb diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Mar,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Mar diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Apr,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Apr diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$May,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="May diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Jun,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Jun diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Jul,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Jul diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Aug,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Aug diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Sept,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Sept diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Oct,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Oct diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Nov,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Nov diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")
barplot(height=revinshoredifchange1996$Dec,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlim=c(-5,5),main="Dec diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in Inshore State Waters")

#inshore diff 2012, 1996, future year 70
par(mfrow=c(2,6))
par(mar=c(3,4,4,0.5))
white=rgb(1,1,1,alpha=0)
barplot(height=revinshoredifchange2012$Jan,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="January",xlab="",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Jan,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
legend(-5,17,legend=c("Warm Year", "Cool Year"),pch=c(22,22),col=c("darkred","purple"),cex=1.5)
#barplot(height=revinshoredifchange70$Jan,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Feb,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="February",xlab="",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Feb,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Feb,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Mar,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="March",xlab="",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Mar,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Mar,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Apr,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="April",xlab="",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Apr,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Apr,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$May,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="May",xlab="",col=white,border="darkred")
barplot(height=revinshoredifchange1996$May,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$May,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Jun,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="June",xlab="",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Jun,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Jun,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
par(mar=c(5,4,2,0.5))
barplot(height=revinshoredifchange2012$Jul,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="July",xlab="SH Difference (In State Waters)",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Jul,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Jul,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Aug,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="August",xlab="SH Difference (In State Waters)",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Aug,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Aug,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Sept,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="September",xlab="SH Difference (In State Waters)",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Sept,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Sept,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Oct,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="October",xlab="SH Difference (In State Waters)",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Oct,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Oct,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Nov,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="November",xlab="SH Difference (In State Waters)",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Nov,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Nov,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revinshoredifchange2012$Dec,names.arg=revinshoredifchange2012$States,horiz=TRUE,xlim=c(-5,5),las=1,main="December",xlab="SH Difference (In State Waters)",col=white,border="darkred")
barplot(height=revinshoredifchange1996$Dec,names.arg=revinshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revinshoredifchange70$Dec,names.arg=revinshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")



#offshore 2012 diff
par(mfrow=c(2,6))
barplot(height=revoffshoredifchange2012$Jan,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Jan diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Feb,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Feb diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Mar,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Mar diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Apr,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Apr diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$May,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="May diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Jun,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Jun diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Jul,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Jul diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Aug,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Aug diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Sept,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Sept diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Oct,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Oct diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Nov,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Nov diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange2012$Dec,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-15,15),main="Dec diff EX 2012 (Warm)",las=1,xlab="% Suit Hab in offshore State Waters")

#offshore 1996 diff
par(mfrow=c(2,6))
barplot(height=revoffshoredifchange1996$Jan,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Jan diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Feb,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Feb diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Mar,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Mar diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Apr,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Apr diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$May,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="May diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Jun,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Jun diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Jul,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Jul diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Aug,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Aug diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Sept,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Sept diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Oct,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Oct diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Nov,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Nov diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")
barplot(height=revoffshoredifchange1996$Dec,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlim=c(-15,15),main="Dec diff EX 1996 (Cool)",las=1,xlab="% Suit Hab in offshore State Waters")


#offshore diff 2012, 1996, future year 70
par(mfrow=c(2,6))
par(mar=c(3,4,4,0.5))
white=rgb(1,1,1,alpha=0)
barplot(height=revoffshoredifchange2012$Jan,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="January",xlab="",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Jan,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
legend(-12,13.3,legend=c("Warm Year", "Cool Year"),pch=c(22,22),col=c("darkred","purple"),cex=1.5)
#barplot(height=revoffshoredifchange70$Jan,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Feb,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="February",xlab="",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Feb,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Feb,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Mar,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="March",xlab="",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Mar,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Mar,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Apr,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="April",xlab="",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Apr,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Apr,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$May,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="May",xlab="",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$May,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$May,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Jun,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="June",xlab="",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Jun,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Jun,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
par(mar=c(5,4,2,0.5))
barplot(height=revoffshoredifchange2012$Jul,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="July",xlab="SH Difference (All Waters)",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Jul,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Jul,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Aug,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="August",xlab="SH Difference (All Waters)",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Aug,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Aug,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Sept,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="September",xlab="SH Difference (All Waters)",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Sept,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Sept,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Oct,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="October",xlab="SH Difference (All Waters)",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Oct,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Oct,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Nov,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="November",xlab="SH Difference (All Waters)",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Nov,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Nov,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")
barplot(height=revoffshoredifchange2012$Dec,names.arg=revoffshoredifchange2012$States,horiz=TRUE,xlim=c(-18,18),las=1,main="December",xlab="SH Difference (All Waters)",col=white,border="darkred")
barplot(height=revoffshoredifchange1996$Dec,names.arg=revoffshoredifchange1996$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="purple",yaxt="n",xaxt="n")
#barplot(height=revoffshoredifchange70$Dec,names.arg=revoffshoredifchange70$States,horiz=TRUE,xlab="",ylab="",add=TRUE,col=white,border="grey",yaxt="n",xaxt="n")


#Sum habitat suitability, 12 values because it is the sum of all states

sum_suitableClim<-c(Jan_offshorestate_percents[12],Feb_offshorestate_percents[12],Mar_offshorestate_percents[12],
                    Apr_offshorestate_percents[12],May_offshorestate_percents[12],Jun_offshorestate_percents[12],Jul_offshorestate_percents[12],
                    Aug_offshorestate_percents[12],Sept_offshorestate_percents[12],Oct_offshorestate_percents[12],Nov_offshorestate_percents[12],Dec_offshorestate_percents[12])

sum_suitable2012<-c(ex201201_offshorestate_percents[12],ex201202_offshorestate_percents[12],ex201203_offshorestate_percents[12],
                    ex201204_offshorestate_percents[12],ex201205_offshorestate_percents[12],ex201206_offshorestate_percents[12],ex201207_offshorestate_percents[12],
                    ex201208_offshorestate_percents[12],ex201209_offshorestate_percents[12],ex201210_offshorestate_percents[12],ex201211_offshorestate_percents[12],ex201212_offshorestate_percents[12])

sum_suitable1996<-c(ex199601_offshorestate_percents[12],ex199602_offshorestate_percents[12],ex199603_offshorestate_percents[12],
                    ex199604_offshorestate_percents[12],ex199605_offshorestate_percents[12],ex199606_offshorestate_percents[12],ex199607_offshorestate_percents[12],
                    ex199608_offshorestate_percents[12],ex199609_offshorestate_percents[12],ex199610_offshorestate_percents[12],ex199611_offshorestate_percents[12],ex199612_offshorestate_percents[12])

sum_suitable0_20<-c(Jan0_20_offshorestate_percents[12],Feb0_20_offshorestate_percents[12],Mar0_20_offshorestate_percents[12],
                    Apr0_20_offshorestate_percents[12],May0_20_offshorestate_percents[12],Jun0_20_offshorestate_percents[12],Jul0_20_offshorestate_percents[12],
                    Aug0_20_offshorestate_percents[12],Sept0_20_offshorestate_percents[12],Oct0_20_offshorestate_percents[12],Nov0_20_offshorestate_percents[12],Dec0_20_offshorestate_percents[12])

sum_suitable20_40<-c(Jan20_40_offshorestate_percents[12],Feb20_40_offshorestate_percents[12],Mar20_40_offshorestate_percents[12],
                     Apr20_40_offshorestate_percents[12],May20_40_offshorestate_percents[12],Jun20_40_offshorestate_percents[12],Jul20_40_offshorestate_percents[12],
                     Aug20_40_offshorestate_percents[12],Sept20_40_offshorestate_percents[12],Oct20_40_offshorestate_percents[12],Nov20_40_offshorestate_percents[12],Dec20_40_offshorestate_percents[12])

sum_suitable40_60<-c(Jan40_60_offshorestate_percents[12],Feb40_60_offshorestate_percents[12],Mar40_60_offshorestate_percents[12],
                     Apr40_60_offshorestate_percents[12],May40_60_offshorestate_percents[12],Jun40_60_offshorestate_percents[12],Jul40_60_offshorestate_percents[12],
                     Aug40_60_offshorestate_percents[12],Sept40_60_offshorestate_percents[12],Oct40_60_offshorestate_percents[12],Nov40_60_offshorestate_percents[12],Dec40_60_offshorestate_percents[12])

sum_suitable60_80<-c(Jan60_80_offshorestate_percents[12],Feb60_80_offshorestate_percents[12],Mar60_80_offshorestate_percents[12],
                     Apr60_80_offshorestate_percents[12],May60_80_offshorestate_percents[12],Jun60_80_offshorestate_percents[12],Jul60_80_offshorestate_percents[12],
                     Aug60_80_offshorestate_percents[12],Sept60_80_offshorestate_percents[12],Oct60_80_offshorestate_percents[12],Nov60_80_offshorestate_percents[12],Dec60_80_offshorestate_percents[12])

sum_suitable70<-c(Jan70_offshorestate_percents[12],Feb70_offshorestate_percents[12],Mar70_offshorestate_percents[12],
                  Apr70_offshorestate_percents[12],May70_offshorestate_percents[12],Jun70_offshorestate_percents[12],Jul70_offshorestate_percents[12],
                  Aug70_offshorestate_percents[12],Sept70_offshorestate_percents[12],Oct70_offshorestate_percents[12],Nov70_offshorestate_percents[12],Dec70_offshorestate_percents[12])


#plot % suitable habitat over year or climatology for each month, so % suitable habitat each month over a year or climatology
#not as interesting
par(mfrow=c(1,1))
plot((sum_suitableClim/sum(sum_suitableClim))*100~c(1:12),type="b",ylim=c(0,20),xlab="Month",ylab="% Suitable Habitat Over Year",lwd=2)
lines((sum_suitable2012/sum(sum_suitable2012))*100~c(1:12),type="b",col="darksalmon",lwd=2)
lines((sum_suitable1996/sum(sum_suitable1996))*100~c(1:12),type="b",col="purple",lwd=2)
lines((sum_suitable0_20/sum(sum_suitable0_20))*100~c(1:12),type="b",col="blue",lwd=2)
lines((sum_suitable20_40/sum(sum_suitable20_40))*100~c(1:12),type="b",col="green",lwd=2)
lines((sum_suitable40_60/sum(sum_suitable40_60))*100~c(1:12),type="b",col="orange",lwd=2)
lines((sum_suitable60_80/sum(sum_suitable60_80))*100~c(1:12),type="b",col="red",lwd=2)
lines((sum_suitable70/sum(sum_suitable70))*100~c(1:12),type="b",col="grey",lwd=2)

legend("topleft",legend=c("Climatology","Warm Year", "Cool Year","0-20","20-40","40-60","60-80","Future Warm"), 
       fill=c("black","darksalmon","purple","blue","green","orange","red","grey"))


#summed habitat suitability (area of grid cellxratio values added up) for each time period for each month
par(mfrow=c(1,1))
plot(sum_suitableClim~c(1:12),type="b",ylim=c(50000,1000000),xlab="Month",ylab="Total Habitat Suitability Index",lwd=2,xaxt="n",cex.lab=1.25)
axis(side=1,at=c(1:12),c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),cex.axis=1.25)
lines(sum_suitable2012~c(1:12),type="b",col="darkred",lwd=2)
lines(sum_suitable1996~c(1:12),type="b",col="purple",lwd=2)
lines(sum_suitable0_20~c(1:12),type="b",col="blue",lwd=2)
lines(sum_suitable20_40~c(1:12),type="b",col="green",lwd=2)
lines(sum_suitable40_60~c(1:12),type="b",col="orange",lwd=2)
lines(sum_suitable60_80~c(1:12),type="b",col="red",lwd=2)
lines(sum_suitable70~c(1:12),type="b",col="grey",lwd=2)

legend("topleft",legend=c("Climatology","Warm Year", "Cool Year","0-20","20-40","40-60","60-80","Future Warm"), 
       lty=c(1,1,1,1,1,1,1,1),lwd=c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5),col=c("black","darkred","purple","blue","green","orange","red","grey"))

#for storymap
par(mfrow=c(1,1),mar=c(5,5,2,2))
plot(sum_suitableClim~c(1:12),type="b",ylim=c(100000,800000),xlab="Month",ylab="Total Available Habitat on U.S. Shelf",lwd=2,yaxt="n",cex.axis=1.5,cex.lab=1.5,xaxt="n")
axis(side=1,at=c(1:12),c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),cex.axis=1.25)
#lines(sum_suitable2012~c(1:12),type="b",col="darkred",lwd=2)
#lines(sum_suitable1996~c(1:12),type="b",col="purple",lwd=2)
lines(sum_suitable0_20~c(1:12),type="b",col="blue",lwd=2)
lines(sum_suitable20_40~c(1:12),type="b",col="green",lwd=2)
lines(sum_suitable40_60~c(1:12),type="b",col="orange",lwd=2)
lines(sum_suitable60_80~c(1:12),type="b",col="red",lwd=2)
#lines(sum_suitable70~c(1:12),type="b",col="grey",lwd=2)

legend("topleft",legend=c("Current","0-20","20-40","40-60","60-80"), 
       lty=c(1,1,1,1,1),lwd=c(1.5,1.5,1.5,1.5,1.5),col=c("black","blue","green","orange","red"),cex=1.25)


































