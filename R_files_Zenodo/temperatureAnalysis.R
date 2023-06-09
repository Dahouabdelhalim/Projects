#comparison of temperatures from early 1970s to 2018 for California
#using antioch pump station (closest to contemporary site with sufficient data)

setwd("~/RDocs/ColiasPostdoc/LabExperiments/Final/DryadSubmission")
tempData<-read.table("./CA(Antioch)_Weather_Data_NOAA.csv",sep=",",header=T)

tempData$DATE<-as.Date(tempData$DATE,"%m/%d/%Y")
tempData$year<-as.integer(format(tempData$DATE,"%Y"))
tempData$month<-as.integer(format(tempData$DATE,"%m"))

#sets up focal periods for our study (10 years prededing each), past=1962-1971, present=2009-2018
tempData$period<-ifelse(tempData$year>=1962&tempData$year<=1971,"past",ifelse(tempData$year>=2009&tempData$year<=2018,"present",""))

#comparing tmax and tmin between 1962-1971 and 2009-2018 in March, July and October
#_________________________________________________________________________________________________________________________________

focalData<-subset(tempData,(period=="past"|period=="present")) #limits data to just focal periods and Antioch *closest to contemporary site
focalData<-subset(focalData,month==3|month==7|month==10) #limits to march, July, October
focalData$period<-factor(focalData$period)

#max temp results
maxRes<-lm(TMAX~period*factor(month),data=focalData)

#min temp results
minRes<-lm(TMIN~period*months(DATE),data=focalData)

#plot for paper; first calculate means and SE
#tempSummary<-data.frame(period=c("past","present","past","present","past","present"),month=c("July","July","March","March","October","October"))
tempSummary<-aggregate(focalData$TMAX,list(period=focalData$period,month=focalData$month),mean)
tempSummary<-aggregate(TMAX~period+month,data=focalData,mean)
tempSummary$mean<-tempSummary$TMAX
tempSummary$sd<-aggregate(TMAX~period+month,data=focalData,sd)$TMAX
tempSummary$n<-aggregate(TMAX~period+month,data=focalData,length)$TMAX
tempSummary$se<-tempSummary$sd/sqrt(tempSummary$n)

#plot minimum temperatures
minSummary<-aggregate(TMIN~period+month,data=focalData,mean)
minSummary$mean<-minSummary$TMIN
minSummary$sd<-aggregate(TMIN~period+month,data=focalData,sd)$TMIN
minSummary$n<-aggregate(TMIN~period+month,data=focalData,length)$TMIN
minSummary$se<-minSummary$sd/sqrt(minSummary$n)

#for paper; both min and max; currently 600x900
par(fig=c(0,1,0.45,1))
plot(tempSummary$mean~tempSummary$month, col=c("black","blue")[tempSummary$period], ylab="Mean Daily\\nHigh Temperature (\\u00B0C)",xlab="", xaxt="n",ylim=c(15,35),xlim=c(2,11),bty="l",pch=c(1,2)[tempSummary$period],cex=1.2, cex.axis=1,cex.lab=1.2,  mgp=c(2.0,0.8,0))
arrows(tempSummary$month, tempSummary$mean-tempSummary$se*1.96, tempSummary$month, tempSummary$mean+tempSummary$se*1.96, length=0.05, angle=90, code=3,lwd=2, col=c("black","blue")[tempSummary$period])
lines(subset(tempSummary,period=="past")$mean~subset(tempSummary,period=="past")$month,lwd=2, lty=2)
lines(subset(tempSummary,period=="present")$mean~subset(tempSummary,period=="present")$month,lwd=2, col="blue")
text(x=2.5, y=33, labels='a', cex=1.5)

par(fig=c(0,1,0,.55),new=T) #minumum plot
plot(minSummary$mean~minSummary$month, col=c("black","blue")[minSummary$period], ylab="Mean Daily\\nLow Temperature (\\u00B0C)",xlab="Month",xaxt="n",ylim=c(0,20),xlim=c(2,11),bty="l",pch=c(1,2)[minSummary$period],cex=1.2, cex.axis=1,cex.lab=1.2,mgp=c(2.0,0.8,0))
axis(1,at=c(3,7,10),labels=c("March","July", "October"),cex.axis=1,mgp=c(2.0,0.8,0))
arrows(minSummary$month, minSummary$mean-minSummary$se*1.96, minSummary$month, minSummary$mean+minSummary$se*1.96, length=0.05, angle=90, code=3,lwd=2, col=c("black","blue")[minSummary$period])
lines(subset(minSummary,period=="past")$mean~subset(minSummary,period=="past")$month,lwd=2, lty=2)
lines(subset(minSummary,period=="present")$mean~subset(minSummary,period=="present")$month,lwd=2, col="blue")
text(x=2.5, y=18, labels='b', cex=1.5)


#for presentation (just march/July; maxima); scale 450x500
tempSummaryShow<-subset(tempSummary,month!=10)
plot(tempSummaryShow$mean~tempSummaryShow$month, col=c("black","blue")[tempSummaryShow$period], ylab="Mean Daily High Temperature (\\u00B0C)",xlab="Month",xaxt="n",ylim=c(15,35),xlim=c(1.5,8.5),bty="l",pch=c(1,2)[tempSummaryShow$period],cex=1.5, cex.axis=1.2,cex.lab=1.3)
axis(1,at=c(3,7,10),labels=c("March","July", "October"),cex.axis=1.2)
arrows(tempSummaryShow$month, tempSummaryShow$mean-tempSummaryShow$se*1.96, tempSummaryShow$month, tempSummaryShow$mean+tempSummaryShow$se*1.96, length=0.05, angle=90, code=3,lwd=2, col=c("black","blue")[tempSummaryShow$period])
lines(subset(tempSummaryShow,period=="past")$mean~subset(tempSummaryShow,period=="past")$month,lwd=2, lty=2)
lines(subset(tempSummaryShow,period=="present")$mean~subset(tempSummaryShow,period=="present")$month,lwd=2, col="blue")
legend(c(1.75,4.5),c(35,31.5),c("1962-1971","2009-2018"),col=c("black","blue"),lwd=2,lty=c(2,1))


#for presentation (just march/July; minima); scale 450x500
minSummaryShow<-subset(minSummary,month!=10)
plot(minSummaryShow$mean~minSummaryShow$month, col=c("black","blue")[minSummaryShow$period], ylab="Mean Daily Low Temperature (\\u00B0C)",xlab="Month",xaxt="n",ylim=c(0,20),xlim=c(1.5,8.5),bty="l",pch=c(1,2)[minSummaryShow$period],cex=1.5, cex.axis=1.2,cex.lab=1.3)
axis(1,at=c(3,7,10),labels=c("March","July", "October"),cex.axis=1.2)
arrows(minSummaryShow$month, minSummaryShow$mean-minSummaryShow$se*1.96, minSummaryShow$month, minSummaryShow$mean+minSummaryShow$se*1.96, length=0.05, angle=90, code=3,lwd=2, col=c("black","blue")[minSummaryShow$period])
lines(subset(minSummaryShow,period=="past")$mean~subset(minSummaryShow,period=="past")$month,lwd=2, lty=2)
lines(subset(minSummaryShow,period=="present")$mean~subset(minSummaryShow,period=="present")$month,lwd=2, col="blue")
legend(c(1.75,4.5),c(20,16.5),c("1962-1971","2009-2018"),col=c("black","blue"),lwd=2,lty=c(2,1))