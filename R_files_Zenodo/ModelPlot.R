library(Hmisc)
library(lme4)
library(segmented)
library(car)
#library(LMERConvenienceFunctions)

modelplot.ContinentalN<-function()
  
{
sink(file = "Output.ContinentalN.txt", append = TRUE, type = c("output"),split = FALSE)  
  
print("-----------------------------------")
print(date())  

##This is the master data file###
datatemp<-data.frame(read.csv("InputData15N.csv"))
dataKeep<-(subset(datatemp,datatemp$Exclude=="Keep"))


###This has the site data
datatempSite<-data.frame(read.csv("InputDataSite.csv"))

###This has the ring width data
dataBAI<-data.frame(read.csv("InputBAI01.csv"))

###merge in site data
data2<-merge(dataKeep, datatempSite, by="GroupSite", all=TRUE)

###modern definition
modern=1970
print(paste("modern = ",modern))

premodern=1930
print(paste("pre-modern = ",premodern))

prepremodern=1890
print(paste("pre-pre-modern = ",prepremodern))


###calculate number of observations
print(paste("number of observations = ",length(dataKeep$N15)))

### subset to modern date range
data3<-subset(data2, data2$AvgYear>=modern)

####standardize 15N by core ######
#standardize data so each core has same mean
###standardization is not necessary any more, but keeping for panel graphs###
data4<-aggregate(data3, by=list(data3$Core), FUN=mean)
dataAllRing<-merge(data2, data4[,c("Core","N15")], by="Core", all=TRUE)

#also can standardize by GroupSite, better for panel graphs to standardize by cores
#data4<-aggregate(data3, by=list(data3$GroupSite), FUN=mean)
#dataAllRing<-merge(data2, data4[,c("GroupSite","N15")], by="GroupSite", all=TRUE)

names(dataAllRing)[names(dataAllRing)=="N15.y"] <- "N15.mean"
names(dataAllRing)[names(dataAllRing)=="N15.x"] <- "N15"


####calculate standardized 15N
dataAllRing[,"N15.std"]<-dataAllRing$N15-dataAllRing$N15.mean

###order data by Year for later plotting
dataAllRing<-dataAllRing[order(dataAllRing$AvgYear),]
write.csv(dataAllRing,"OutputDataAllRing.csv")

##summary file for functional groups
dataAllRingFxnlSumm<-aggregate(dataAllRing$AvgYear, by=list(dataAllRing$Core, dataAllRing$FxnlGrp1), FUN=max)
names(dataAllRingFxnlSumm)[1] <- "Core"
names(dataAllRingFxnlSumm)[2] <- "FxnlGrp1"
dataAllRingFxnlSumm<-dataAllRingFxnlSumm[c("Core", "FxnlGrp1")]

##summarize by decade
dataAllRing[,"Decade"]<-floor(dataAllRing$AvgYear/10)*10+5
#sums by groupsite and decade
dataAllRingByDecade<-aggregate(cbind(dataAllRing$AvgYear,dataAllRing$N15.std)~GroupSite+Decade, data=dataAllRing, mean, na.rm=TRUE)
dataAllRingByDecade<-aggregate(dataAllRingByDecade, by=list(dataAllRingByDecade$Decade), FUN=mean)
dataAllRingByDecade<-dataAllRingByDecade[-c(1:3)]
names(dataAllRingByDecade)[1] <- "AvgYearDecade"
names(dataAllRingByDecade)[2] <- "AvgN15.std"

write.csv(dataAllRingByDecade,"OutputDataAllRingByDecade.csv")


###summarize by decade for just pre-pre-modern trees
dataAllRingPPM<-dataAllRing[which(dataAllRing$MinYear<=prepremodern),]

dataAllRingPPMByDecade<-aggregate(cbind(dataAllRingPPM$AvgYear,dataAllRingPPM$N15.std)~GroupSite+Decade, data=dataAllRingPPM, mean, na.rm=TRUE)
dataAllRingPPMByDecade<-aggregate(dataAllRingPPMByDecade, by=list(dataAllRingPPMByDecade$Decade), FUN=mean)
dataAllRingPPMByDecade<-dataAllRingPPMByDecade[-c(1:3)]
names(dataAllRingPPMByDecade)[1] <- "AvgYearDecade"
names(dataAllRingPPMByDecade)[2] <- "AvgN15.std"


######################
###run regressions####
######################
###core diagnostics
dataCoreDiag<-aggregate(data2$AvgYear, by=list(data2$Core), FUN=max)
names(dataCoreDiag)[1] <- "Core"
names(dataCoreDiag)[2] <- "YearMax"
dataCoreDiag$YearMin<-aggregate(data2$AvgYear, by=list(data2$Core), FUN=min)[,2]
dataCoreDiag$CountModern<-aggregate(data3$AvgYear, by=list(data3$Core), FUN=length)[,2]

write.csv(dataCoreDiag, "OutputDataCoreDiag.csv")

##stats on youngest cores, and number of cores overall
print(paste("Youngest Core = ", max(dataCoreDiag$YearMin)))
print(paste("# cores overall = ",nrow(dataCoreDiag)))

### slopes for each core using 15N data that is not standardized
dataModern<-subset(dataAllRing, dataAllRing$AvgYear>=modern)
fit<-lmList(N15~AvgYear|Core,  data=dataModern)
CoreCoef<-(coef(fit))
CoreCoef$Core<-rownames(CoreCoef)
write.csv(CoreCoef,"CoreCoef.csv")

#Pre modern slopes (40 y before)
dataPreModern<-subset(dataAllRing, dataAllRing$AvgYear>=premodern & dataAllRing$AvgYear<=modern)
fit<-lmList(N15.std~AvgYear|Core,  data=dataPreModern)
CoreOldCoef<-(coef(fit))
CoreOldCoef$Core<-rownames(CoreOldCoef)

#Pre Pre modern slopes (80 y before)
dataPPM<-subset(dataAllRing, dataAllRing$AvgYear>=prepremodern & dataAllRing$AvgYear<=premodern)
fit<-lmList(N15.std~AvgYear|Core,  data=dataPPM)
CoreOld2Coef<-(coef(fit))
CoreOld2Coef$Core<-rownames(CoreOld2Coef)

  
### BAI summary statistics
###BAI slopes
dataBAIModern<-subset(dataBAI, dataBAI$Year>=modern)
fit<-lmList(BAI~Year|Core,  data=dataBAIModern)
BAICoreCoef<-(coef(fit))
BAICoreCoef$Core<-rownames(BAICoreCoef)

###mean BAI
BAICoreMeanModern<-aggregate(dataBAIModern, by=list(dataBAIModern$Core), FUN=mean)

###merge BAI mean and slopes
BAICoreCoef2<-merge(BAICoreMeanModern, BAICoreCoef,by="Core", all=TRUE)
BAICoreCoef2<-BAICoreCoef2[c("Core","BAI","(Intercept)","Year.y")]

###merge BAI and core data
dataCore1<-merge(CoreCoef,BAICoreCoef2 , by="Core", all.x=TRUE)
dataCore1<-merge(dataCore1,CoreOldCoef , by="Core", all.x=TRUE)
dataCore1<-merge(dataCore1,CoreOld2Coef , by="Core", all.x=TRUE)
colnames(dataCore1)<-c("Core","Intercept.15NYear","Slope.15NYear","MeanBAI","Intercept.BAIYear","Slope.BAI.Year", "Intercept.15NYear.PM","Slope.15NYear.PM", "Intercept.15NYear.PPM","Slope.15NYear.PPM")


### slopes for each group site starting from raw
fit<-lmList(N15.std~AvgYear|GroupSite,  data=dataModern)
GroupSiteCoef<-(coef(fit))
GroupSiteCoef$GroupSite<-rownames(GroupSiteCoef)
dataGroupSite1<-merge(dataModern, GroupSiteCoef, by="GroupSite", all=TRUE)

##merge BAI data into Group Data
##this gets df that links Cores and Group Sites
dataCoreGroup<-aggregate(datatemp, by=list(datatemp$Core), FUN=mean)
dataCoreGroup<-dataCoreGroup[c("GroupSite", "Core")]
dataCore2<-merge(dataCoreGroup, dataCore1, by="Core", all=TRUE)
write.csv(dataCore2,"OutputDataCore2.csv")


##aggregate core level data to Group level
###note in aggregate, na.rm=TRUE allows calculating mean with NA data present in category
dataCore2GS<-aggregate(dataCore2, by=list(dataCore2$GroupSite), FUN=mean,na.action=na.pass, na.rm=TRUE)
dataCore2GS<-dataCore2GS[c(-1,-2)]

data7GS<-aggregate(dataGroupSite1, by=list(dataGroupSite1$GroupSite), FUN=mean)
write.csv(data7GS,"data7GS.csv")
drops <- c("Species","Lab", "Group1","Core","FxnlGrp1","StartYear","EndYear","Exclude","State", "Fxnl")
dataAllGS<-data7GS[ , !(names(data7GS) %in% drops)]

dataAllGS<-merge(dataCore2GS, dataAllGS, by="GroupSite", all=TRUE)
drops <- c("Group.1")
dataAllGS<-dataAllGS[ , !(names(dataAllGS) %in% drops)]
dataAllGS<-as.data.frame(dataAllGS)

names(dataAllGS)[names(dataAllGS) == 'Intercept.15NYear'] <- 'Intercept.15N.Year.Core'
names(dataAllGS)[names(dataAllGS) == 'Slope.15NYear'] <- 'Slope.15N.Year.Core'
names(dataAllGS)[names(dataAllGS) == '(Intercept)'] <- 'Intercept.15N.Year.GS'
names(dataAllGS)[names(dataAllGS) == 'AvgYear.y'] <- 'Slope.15N.Year.GS'

write.csv(dataAllGS, "OutputDataAllGS.csv")

print(paste("MAP range = ", min(dataAllGS$MAP.mm, na.rm=TRUE), "to", max(dataAllGS$MAP.mm, na.rm=TRUE), "...range = ", max(dataAllGS$MAP.mm, na.rm=TRUE)-min(dataAllGS$MAP.mm, na.rm=TRUE)))
print(paste("MAT range = ", min(dataAllGS$MAT, na.rm=TRUE), "to", max(dataAllGS$MAT, na.rm=TRUE), "...range = ", max(dataAllGS$MAT, na.rm=TRUE)-min(dataAllGS$MAT, na.rm=TRUE)))
print(paste("Ndep range = ", min(dataAllGS$Ndep, na.rm=TRUE), "to", max(dataAllGS$Ndep, na.rm=TRUE), "...range = ", max(dataAllGS$Ndep, na.rm=TRUE)-min(dataAllGS$Ndep, na.rm=TRUE)))

#calculate # sites each year
###make a dataframe that runs from 1850-2015
##repeat counting how many sites have a minimum age greater than the year-1, e.g. how many
##sites have a minimum age older than 1799 for 1800 data point
sitecount<-data.frame(Year=1800:2015)
for(i in 1:216){
sitecount$count[i]<-sum(dataAllGS$MinYear<i+1798, na.rm=TRUE)
}

##now merge core data with Group level data
dataCore3<-merge(dataCore2, dataAllGS, by="GroupSite", all=TRUE)

###now merge it with functional group data
dataCore3<-merge(dataCore3, dataAllRingFxnlSumm, by="Core", all=TRUE)

write.csv(dataCore3,"OutputDataCore3.csv")

### regress slopes vs. site metrics

print("regression results using slope data derived from averaging core data")
fit<-lm(Slope.15N.Year.Core~log.MAP+MAT  + Ndep  , data=dataAllGS)
print(summary(fit)$coefficients)
print("Outlier Test")
print(outlierTest(fit))

dataAllRingNoOutlier <- dataAllRing[which(dataAllRing$GroupSite!=240),]
dataAllGSNoOutlier <- dataAllGS[which(dataAllGS$GroupSite!=240),]
dataCore3NoOutlier <- dataCore3[which(dataCore3$GroupSite!=240),]


###Test BAI against residuals of model
print("BAI ranges")
print(quantile(dataCore1$MeanBAI,probs=0.025, na.rm=TRUE))
print(quantile(dataCore1$MeanBAI,probs=0.975, na.rm=TRUE))

print("BAI slope ranges")
print(quantile(dataCore1$Slope.BAI.Year,probs=0.025, na.rm=TRUE))
print(quantile(dataCore1$Slope.BAI.Year,probs=0.975, na.rm=TRUE))

dataAllGS$resid<-residuals(fit)

###Test age of trees against residuals of model

print("Test residuals against age of oldest tree")
dataCoreDiag2<-merge(dataCoreDiag, dataCore2,by="Core", all=TRUE)
dataCoreDiag3<-(aggregate(dataCoreDiag2$YearMin, by=list(dataCoreDiag2$GroupSite), FUN=min, na.rm=TRUE))
names(dataCoreDiag3)[1] <- "GroupSite"
names(dataCoreDiag3)[2] <- "YearMin"
dataAllGStemp<-merge(dataAllGS, dataCoreDiag3, by="GroupSite", all=TRUE)
write.csv(dataAllGStemp,"dataAllGStemp.csv")

fitresid<-lm(resid~YearMin, data=dataAllGStemp, na.action=na.omit)
print(summary(fitresid)$coefficients)
print(summary(fitresid))
print(outlierTest(fitresid))

####random effects regression 
fitnested<-lmer(Slope.15NYear~log.MAP+MAT + Ndep  + (1|GroupSite) , data=dataCore3)
print("*************************")
print("###random effects regression without functional group")
print(summary(fitnested))
coefs <- data.frame(coef(summary(fitnested)))
coefs$p.value<-(2 * (1 - pnorm(abs(coefs$t.value))))
print(coef(fitnested)$FxnlGrp1)
print(coefs)

####nested regression 
print("###random effects regression with functional group")
fitnested<-lmer(Slope.15NYear~log.MAP+MAT + Ndep + FxnlGrp1 + (1|GroupSite) , data=dataCore3)
print(summary(fitnested))
coefs <- data.frame(coef(summary(fitnested)))
coefs$p.value<-(2 * (1 - pnorm(abs(coefs$t.value))))
print(coef(fitnested)$FxnlGrp1)
print(coefs)

###estimate declines pre-modern
dataCore3PM <- dataCore3[which(dataCore3$MinYear<=premodern),]
write.csv(dataCore3PM,"OutputDataCore3PM.csv")
dataCore3PMGS<-aggregate(dataCore3PM, by=list(dataCore3PM$GroupSite), FUN=mean)
slopePMmean<-round(mean(dataCore3PMGS$Slope.15NYear.PM.y),4)
slopePMse<-sd(dataCore3PMGS$Slope.15NYear.PM.y)/sqrt(length(dataCore3PMGS$Slope.15NYear.PM.y))
print(paste("mean slope 15N from ", premodern, "to ", modern, "= ", slopePMmean, " ± ", (round(slopePMse,3)), "n = ",length(dataCore3PMGS$Slope.15NYear.PPM.y), "total change over 40 years = ", slopePMmean*40))
print(paste("# cores older than 1930 =", nrow(dataCore3PMGS)))

###estimate declines before pre-modern
dataCore3PPM <- dataCore3[which(dataCore3$MinYear<=prepremodern),]
write.csv(dataCore3PM,"OutputDataCore3PPM.csv")
dataCore3PPMGS<-aggregate(dataCore3PPM, by=list(dataCore3PPM$GroupSite), FUN=mean)

slopePPMmean<-round(mean(dataCore3PPMGS$Slope.15NYear.PPM.y),6)
slopePPMse<-sd(dataCore3PPMGS$Slope.15NYear.PPM.y)/sqrt(length(dataCore3PPMGS$Slope.15NYear.PPM.y))
print(paste("mean slope 15N from ", prepremodern, " to ", premodern, " = ", slopePPMmean, " ± ", (round(slopePPMse,3)), "n = ",length(dataCore3PPMGS$Slope.15NYear.PPM.y)))
print(paste("# cores older than 1890 =", nrow(dataCore3PPMGS)))

###calculate slopes continuously
DateStart=1850
DateEnd=2015
DateRange=40

dataSlopeRange<-dataCoreGroup[c("GroupSite", "Core")]
for (i in 1:(DateEnd-DateStart+1)){
dataAllRing2<-dataAllRing[c("GroupSite","Core","AvgYear","N15")]
###calculate minimum year for each core from dataAllRing
###Then merge it back into dataAllRing2
dataAllRingMin<-aggregate(dataAllRing$AvgYear, by=list(dataAllRing$Core), FUN=min)
names(dataAllRingMin)[names(dataAllRingMin)=="Group.1"] <- "Core"
dataAllRing2<-merge(dataAllRing2,dataAllRingMin,by="Core", all=TRUE)

dataRange<-subset(dataAllRing2, dataAllRing2$AvgYear>=DateStart-DateRange/2+i-1 & dataAllRing2$AvgYear<DateStart+DateRange/2+i-1)

# if want to standardize the data by group site, substitute this code for next lines
#dataRange2<-aggregate(dataRange, by=list(dataRange$GroupSite), FUN=mean)
#dataRange3<-merge(dataRange, dataRange2[,c("GroupSite","N15")], by="GroupSite", all=TRUE)

#standardizing by core
dataRange2<-aggregate(dataRange, by=list(dataRange$Core), FUN=mean)
dataRange3<-merge(dataRange, dataRange2[,c("Core","N15")], by="Core", all=TRUE)

names(dataRange3)[names(dataRange3)=="N15.y"] <- "N15.mean"
names(dataRange3)[names(dataRange3)=="N15.x"] <- "N15"
names(dataRange3)[names(dataRange3)=="x"] <- "YearMin"

####calculate standardized 15N
dataRange3[,"N15.std"]<-dataRange3$N15-dataRange3$N15.mean
dataRange3<-dataRange3[which(dataRange3$YearMin<DateStart-DateRange/2+i-1),]

fit<-lmList(N15.std~AvgYear|Core,  data=dataRange3)
CoreCoef<-(coef(fit))
CoreCoef$Core<-rownames(CoreCoef)
CoreCoef<-CoreCoef[c("Core", "AvgYear")]

dataSlopeRange<-merge(dataSlopeRange,CoreCoef, by="Core", all=TRUE)
names(dataSlopeRange)[i+2] <- DateStart-1+i

}
write.csv(dataSlopeRange,"dataSlopeRange.csv")

print("dataSlopeRangeCount")
dataSlopeRangeMatrix = as.matrix(dataSlopeRange)
dataSlopeRangeCount<-colSums(dataSlopeRangeMatrix != 0, na.rm = TRUE)

dataSlopeRangeCount<-dataSlopeRangeCount[-c(1:3)]
write.csv(dataSlopeRangeCount,"dataSlopeRangeCount.csv")

dataSlopeRangeMean<-aggregate(dataSlopeRange, by=list(dataSlopeRange$GroupSite), FUN=mean, na.action=na.pass, na.rm=TRUE)

write.csv(dataSlopeRangeMean, "dataSlopeRangeMean.csv")

dataSlopeRangeMean2<-colMeans(dataSlopeRangeMean[,4:(DateEnd-DateStart+4)], na.rm=TRUE)

write.csv(dataSlopeRangeMean2, "dataSlopeRangeMean2.csv")

###now stack the core-level slope data so that can plot it out. 
dataSlopeRangeMean3<-na.omit(stack(dataSlopeRangeMean, select=4:(DateEnd-DateStart+4)))
write.csv(dataSlopeRangeMean3, "dataSlopeRangeMean3.csv")



sink()
################
#### Plots  #####
################

##### plot multipanel figure
pdf(file="15N.Panels.pdf",11,8, encoding="WinAnsi")
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), 5,4, byrow = FALSE), heights=c(5,5,5,5,7), widths=c(6,5,5,5))

Sites=c("200","201","202","203","204","205","208","209","210","212","213","214","216","217","218","219","220","221","222","223","224","226","227","228","229","230","231","233","234","235","236","238","239","240","241","242","243","244","245","246","247","249","250","251","252","253","254","255","257")
para1=c(0,0,0,0,5,0,0,0,0,5,0,0,0,0,5,0,0,0,0,5)
para2=c(5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
ax1=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE)
ax2=c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)

for(i in 1:20){
datapanel<-subset(dataAllRing, dataAllRing$GroupSite==Sites[i])

par(mar=c(para1[i],para2[i],1,1))

plot(datapanel$AvgYear, datapanel$N15.std, type="p", pch=16, cex.axis=1, cex.lab=1.5, cex=1, col=rgb(100,100,100,75,maxColorValue=255), xlab="Year", ylab=expression(paste('Wood ',delta^{15},'N',' (‰)')),xlim=c(1750,2015),ylim=c(-4,6),axes=FALSE, frame.plot=TRUE)

lines(smooth.spline(datapanel$AvgYear, datapanel$N15.std, df=10), col="red", lwd=1)

Axis(side=1, labels=ax1[i], cex.axis=1, cex.lab=1.5)
Axis(side=2, labels=ax2[i], cex.axis=1, cex.lab=1.5)
minor.tick(nx=5, ny=1, tick.ratio=0.5)
mtext(Sites[i], side = 3, line=-1.2,adj=0.02,cex=0.75, col = "black")
}

dev.off()

pdf(file="15N.Panelsb.pdf",11,8, encoding="WinAnsi")
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), 5,4, byrow = FALSE), heights=c(5,5,5,5,7), widths=c(6,5,5,5))

for(i in 1:20){
  datapanel<-subset(dataAllRing, dataAllRing$GroupSite==Sites[i+20])
  par(mar=c(para1[i],para2[i],1,1))
  
  plot(datapanel$AvgYear, datapanel$N15.std, type="p", pch=16, cex.axis=1, cex.lab=1.5, cex=1, col=rgb(100,100,100,75,maxColorValue=255), xlab="Year", ylab=expression(paste('Wood ',delta^{15},'N',' (‰)')),xlim=c(1750,2015),ylim=c(-4,6),axes=FALSE, frame.plot=TRUE)
  lines(smooth.spline(datapanel$AvgYear, datapanel$N15.std, df=10), col="red", lwd=1)
  
  Axis(side=1, labels=ax1[i], cex.axis=1, cex.lab=1.5)
  Axis(side=2, labels=ax2[i], cex.axis=1, cex.lab=1.5)
  minor.tick(nx=5, ny=1, tick.ratio=0.5)
  mtext(Sites[i+20], side = 3, line=-1.2,adj=0.02,cex=0.75, col = "black")
  
  
}

dev.off()

pdf(file="15N.Panelsc.pdf",11,8, encoding="WinAnsi")
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), 5,4, byrow = FALSE), heights=c(5,5,5,5,7), widths=c(6,5,5,5))

para1=c(0,0,0,0,5,0,0,0,0,5,0,0,0,0,5,0,0,0,0,5)
para2=c(5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
ax1=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE)
ax2=c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)


for(i in 1:9){
  datapanel<-subset(dataAllRing, dataAllRing$GroupSite==Sites[i+40])
  par(mar=c(para1[i],para2[i],1,1))
  
  plot(datapanel$AvgYear, datapanel$N15.std, type="p", pch=16, cex.axis=1, cex.lab=1.5, cex=1, col=rgb(100,100,100,75,maxColorValue=255), xlab="Year", ylab=expression(paste('Wood ',delta^{15},'N',' (‰)')),xlim=c(1750,2015),ylim=c(-4,6),axes=FALSE, frame.plot=TRUE)
  lines(smooth.spline(datapanel$AvgYear, datapanel$N15.std, df=10), col="red", lwd=1)
  
  Axis(side=1, labels=ax1[i], cex.axis=1, cex.lab=1)
  Axis(side=2, labels=ax2[i], cex.axis=1, cex.lab=1)
  minor.tick(nx=5, ny=1, tick.ratio=0.5)
  mtext(Sites[i+40], side = 3, line=-1.2,adj=0.02,cex=0.75, col = "black")

}

dev.off()

###custom partial residual plots
###Custom partial residual plot for slope data averaged from cores 
pdf(file="15NLeverageCore.pdf",10,3.5, encoding="WinAnsi")
layout(matrix(c(1,2,3),1,3, byrow = FALSE), widths=c(6,5,5))

lm.results.Core = lm(Slope.15N.Year.Core ~ log.MAP + MAT  + Ndep , data=dataAllGS)

dataAllGSCent <- transform(dataAllGS, log.MAPCent = scale(dataAllGS$log.MAP,center=TRUE, scale=FALSE), MATCent = scale(dataAllGS$MAT,center=TRUE, scale=FALSE), NdepCent = scale(dataAllGS$Ndep,center=TRUE, scale=FALSE))
lm.terms.results<-lm(Slope.15N.Year.Core ~ log.MAPCent + MATCent  + NdepCent , data=dataAllGSCent)
pterms <- predict(lm.results.Core, type="terms")
partial.residuals <- apply(pterms,2,function(x)x+resid(lm.results.Core)+mean(dataAllGSCent$Slope.15N.Year.Core))

xlabel=c("log MAP", "MAT (°C)", expression(paste("N deposition (g m"^{-2}," d"^-1,")")))
ylabelTF=c(TRUE, FALSE,FALSE)
mtextlabel=c("a","b","c")
MAPlabel=c(2.6, 2.8, 3, 3.2, 3.4)

pary=c(5,1,1)
# the model in lm.results includes the response in first column, so we index with i + 1
for(i in 1:3){
  par(mar=c(5,pary[i],1,1))
  plot(x=lm.results.Core$model[,(i+1)],y=partial.residuals[,i], type="p", pch=16, cex.axis=1.5, cex.lab=1.5, cex=1, col="black", xlab=xlabel[i], ylab=expression(paste(Delta,' Wood ',delta^{15},'N',' (‰ y'^"-1",")")), axes=FALSE, frame.plot=TRUE)
  modelx<-(lm.results.Core$model[,(i+1)])
  mylm<-lm(partial.residuals[,i] ~ modelx)
  ###!!!!!!make panel 3 dashed!!!!!!!!!
  ltypanel<-ifelse(i==3,2,1)
  abline(mylm, col="red", lty=ltypanel)
  
  xseq<-lm.results.Core$model[,(i+1)]
  xseq<-xseq[order(xseq)]
  
  prd <- predict(mylm, newdata=data.frame(modelx=xseq), interval="confidence") 
  lines(xseq, prd[,2], lty=2)
  lines(xseq, prd[,3], lty=2) 
  write.csv(prd, paste("OutputPartialResidPred",i,".csv"))
  
  Axis(side=1, labels=TRUE, cex.axis=1.5, cex.lab=1)
  Axis(side=2, labels=ylabelTF[i], cex.axis=1.5, cex.lab=1)
  minor.tick(nx=5, ny=5, tick.ratio=0.5)

}
dev.off()

###plot of running means of slopes

pdf(file="15NSlopesYear.pdf",7,7, encoding="WinAnsi")
layout(matrix(c(1,1,2),3,1, byrow = FALSE), heights=c(5,5,5))

#1
par(mar=c(1,5,1,1))
dataSlopeRangeMean3$ind<-as.numeric(as.character(dataSlopeRangeMean3$ind))
plot(dataSlopeRangeMean3$ind, dataSlopeRangeMean3$values, type="p", pch=16, cex.axis=1.5, cex.lab=1.5, cex=0.5, col=rgb(100,100,100,50,maxColorValue=255), xlab="Year", ylab=expression(paste(Delta,' Wood ',delta^{15},'N',' (‰ y'^"-1",")")), ylim=c(-0.05, 0.05), xlim=c(DateStart,2015), axes=FALSE, frame.plot=TRUE)

plx<-predict(loess(dataSlopeRangeMean3$values ~ dataSlopeRangeMean3$ind, span=0.2, family="gaussian", evaluation=100), se=T)

dataSlopeRangeMean3$fit<-plx$fit

dataSlopeRangeMean3pos<-dataSlopeRangeMean3[which(dataSlopeRangeMean3$fit>0),]
dataSlopeRangeMean3neg<-dataSlopeRangeMean3[which(dataSlopeRangeMean3$fit<0),]
points(dataSlopeRangeMean3pos$ind,dataSlopeRangeMean3pos$fit, lwd=2, pch=16, col="red", cex=1)
points(dataSlopeRangeMean3neg$ind,dataSlopeRangeMean3neg$fit, lwd=2, pch=16, col="blue", cex=1)

lines(dataSlopeRangeMean3$ind,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(dataSlopeRangeMean3$ind,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)


Axis(side=1, labels=FALSE, cex.axis=1, cex.lab=1)
Axis(side=2, labels=TRUE, cex.axis=1, cex.lab=1)
minor.tick(nx=5, ny=5, tick.ratio=0.5)
mtext("a", side = 3, line=-2,adj=0.02,cex=1, col = "black")

#2
par(mar=c(5,5,1,1))

##take predicted column and put it into working df again. 
##then replace values for predicted slopes with mean modern slope
##then aggregate to one value per year
##then generate a cumulative sum to generate wood d15N data
dataSlopeRangeMean3$fit2<-plx$fit
dataSlopeRangeMean3$fit2[which(dataSlopeRangeMean3$ind>=1990)]=mean(dataAllGSCent$Slope.15N.Year.Core)
dataSlopeRangeMean4<-aggregate(dataSlopeRangeMean3, by=list(dataSlopeRangeMean3$ind), FUN=mean)
dataSlopeRangeMean4$Cum<-cumsum(dataSlopeRangeMean4$fit2)

dataAllRing18301870<-subset(dataAllRing, dataAllRing$AvgYear>=1830 & dataAllRing$AvgYear<=1870)
startdel1850<-mean(dataAllRing18301870$N15)

sink(file = "Output.ContinentalN.txt", append = TRUE, type = c("output"),split = FALSE)    
print("initial del 15N 1830-1870")
print(startdel1850)
sink()
write.csv(dataSlopeRangeMean4, "dataSlopeRangeMean4.csv")

plot(dataAllRing$AvgYear, dataAllRing$N15.std, type="p", pch=16, cex.axis=1.5, cex.lab=1.5, cex=0.5, col="white", xlab="Year", ylab=expression(paste('Calc Wood ',delta^{15},'N',' (‰)')),ylim=c(-1.5,0.5),xlim=c(DateStart,2015), axes=FALSE, frame.plot=TRUE)


dataSlopeRangeMean4pos<-dataSlopeRangeMean4[which(dataSlopeRangeMean4$fit>0),]
dataSlopeRangeMean4neg<-dataSlopeRangeMean4[which(dataSlopeRangeMean4$fit<0),]
points(dataSlopeRangeMean4pos$ind,dataSlopeRangeMean4pos$Cum+startdel1850, lwd=2, pch=16, col="red", cex=1)
points(dataSlopeRangeMean4neg$ind,dataSlopeRangeMean4neg$Cum+startdel1850, lwd=2, pch=16, col="blue", cex=1)

Axis(side=1, labels=TRUE, cex.axis=1, cex.lab=1)
Axis(side=2, labels=TRUE, cex.axis=1, cex.lab=1)
minor.tick(nx=5, ny=5, tick.ratio=0.5)
mtext("b", side = 3, line=-2,adj=0.02,cex=1, col = "black")

dev.off()


###plot num sites x year

pdf(file="15NCountSitesYear.pdf",6,6, encoding="WinAnsi")

#1
par(mar=c(5,5,1,1))

plot(sitecount[,1], sitecount[,2], type="p", pch=16, cex.axis=1.5, cex.lab=1.5, cex=0.5, col="black", xlab="Year", ylab="# Sites", xlim=c(1800,2015), axes=FALSE, frame.plot=TRUE)

Axis(side=1, labels=TRUE, cex.axis=1, cex.lab=1)
Axis(side=2, labels=TRUE, cex.axis=1, cex.lab=1)
minor.tick(nx=5, ny=5, tick.ratio=0.5)

dev.off()

###plot binned centered 15N data by decade

pdf(file="15NstdDecade1850.pdf",6,6, encoding="WinAnsi")

par(mar=c(5,5,1,1))

plot(dataAllRingByDecade$AvgYearDecade, dataAllRingByDecade$AvgN15.std, type="p", pch=16, cex.axis=1.5, cex.lab=1.5, cex=0.5, col="black", xlab="Year", ylab=expression(paste('Wood ',delta^{15},'N Anomaly',' (‰)')), xlim=c(1850,2015), ylim=c(-2,3), axes=FALSE, frame.plot=TRUE)
points(dataAllRing$AvgYear,dataAllRing$N15.std, pch=16, col="gray", cex=0.4)

points(dataAllRingByDecade$AvgYearDecade, dataAllRingByDecade$AvgN15.std, pch=16, col="black", cex=1)

lines(smooth.spline(dataAllRingByDecade$AvgYearDecade, dataAllRingByDecade$AvgN15.std, df=20), col="black", lwd=2)



plx<-predict(loess(dataAllRingByDecade$AvgN15.std ~ dataAllRingByDecade$AvgYearDecade, span=0.2, family="gaussian", evaluation=100), se=T)

dataAllRingByDecade$fit<-plx$fit

xseq<-seq(from=1, to=length(dataAllRingByDecade[,1]), by =1)

polygon(x=c(dataAllRingByDecade[,1], rev(dataAllRingByDecade[,1])), y=c(plx$fit[xseq]- qt(0.975,plx$df)*plx$se[xseq], plx$fit[rev(xseq)]+ qt(0.975,plx$df)*plx$se[rev(xseq)]), col=rgb(0,0,255,100,maxColorValue=255))


lines(dataAllRingByDecade$AvgYearDecade,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(dataAllRingByDecade$AvgYearDecade,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)

Axis(side=1, labels=TRUE, cex.axis=1, cex.lab=1)
Axis(side=2, labels=TRUE, cex.axis=1, cex.lab=1)
minor.tick(nx=5, ny=5, tick.ratio=0.5)

dev.off()

sink(file = "Output.ContinentalN.txt", append = TRUE, type = c("output"),split = FALSE)  

print(paste("For all trees: difference 2010-1850", round(dataAllRingByDecade[42,2]-dataAllRingByDecade[26,2],2), " ± ", (round((plx$se[26]^2+plx$se[42]^2)^0.5,2))))
print(paste("For all trees: difference 2010-1930", round(dataAllRingByDecade[42,2]-dataAllRingByDecade[34,2],2), " ± ", (round((plx$se[34]^2+plx$se[42]^2)^0.5,2))))

sink()



###plot binned centered 15N data by decade for PrePreModern

pdf(file="15NstdPPMDecade1850.pdf",6,6, encoding="WinAnsi")

par(mar=c(5,5,1,1))

plot(dataAllRingPPMByDecade$AvgYearDecade, dataAllRingPPMByDecade$AvgN15.std, type="p", pch=16, cex.axis=1.5, cex.lab=1.5, cex=0.5, col="black", xlab="Year", ylab=expression(paste('Wood ',delta^{15},'N Anomaly',' (‰)')), xlim=c(1850,2015), ylim=c(-2,3), axes=FALSE, frame.plot=TRUE)
points(dataAllRingPPM$AvgYear,dataAllRingPPM$N15.std, pch=16, col="gray", cex=0.4)

points(dataAllRingPPMByDecade$AvgYearDecade, dataAllRingPPMByDecade$AvgN15.std, pch=16, col="black", cex=1)

lines(smooth.spline(dataAllRingPPMByDecade$AvgYearDecade, dataAllRingPPMByDecade$AvgN15.std, df=20), col="black", lwd=2)



plx<-predict(loess(dataAllRingPPMByDecade$AvgN15.std ~ dataAllRingPPMByDecade$AvgYearDecade, span=0.2, family="gaussian", evaluation=100), se=T)

dataAllRingPPMByDecade$fit<-plx$fit

xseq<-seq(from=1, to=length(dataAllRingPPMByDecade[,1]), by =1)

polygon(x=c(dataAllRingPPMByDecade[,1], rev(dataAllRingPPMByDecade[,1])), y=c(plx$fit[xseq]- qt(0.975,plx$df)*plx$se[xseq], plx$fit[rev(xseq)]+ qt(0.975,plx$df)*plx$se[rev(xseq)]), col=rgb(0,0,255,100,maxColorValue=255))


lines(dataAllRingPPMByDecade$AvgYearDecade,plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
lines(dataAllRingPPMByDecade$AvgYearDecade,plx$fit + qt(0.975,plx$df)*plx$se, lty=2)


Axis(side=1, labels=TRUE, cex.axis=1, cex.lab=1)
Axis(side=2, labels=TRUE, cex.axis=1, cex.lab=1)
minor.tick(nx=5, ny=5, tick.ratio=0.5)

dev.off()

sink(file = "Output.ContinentalN.txt", append = TRUE, type = c("output"),split = FALSE)  

print(paste("For PPM trees: difference 2010-1850", round(dataAllRingPPMByDecade[42,2]-dataAllRingPPMByDecade[26,2],2), " ± ", (round((plx$se[26]^2+plx$se[42]^2)^0.5,2))))


sink()

#################
###print out means
##################
sink(file = "Output.ContinentalN.txt", append = TRUE, type = c("output"),split = FALSE)  

print("means of slopes and parameters")

slopemean<-round(mean(dataAllGSCent$Slope.15N.Year.Core),4)
slopese<-sd(dataAllGSCent$Slope.15N.Year.Core)/sqrt(length(dataAllGSCent$Slope.15N.Year.Core))
print(paste("mean slope 15N modern = ", slopemean, " ± ", (round(slopese,3))))


###print means of predictor variables.
logMAPmean<-round(mean(dataAllGSCent$log.MAP, na.rm=TRUE),4)
logMAPse<-sd(dataAllGSCent$log.MAP, na.rm=TRUE)/sqrt(length(dataAllGSCent$log.MAP))
print(paste("mean log MAP  = ", logMAPmean, " ± ", (round(logMAPse,3))))

MATmean<-round(mean(dataAllGSCent$MAT, na.rm=TRUE),4)
MATse<-sd(dataAllGSCent$MAT, na.rm=TRUE)/sqrt(length(dataAllGSCent$MAT))
print(paste("mean MAT  = ", MATmean, " ± ", (round(MATse,3))))

Ndepmean<-round(mean(dataAllGSCent$Ndep, na.rm=TRUE),4)
Ndepse<-sd(dataAllGSCent$Ndep, na.rm=TRUE)/sqrt(length(dataAllGSCent$Ndep))
print(paste("mean Ndep  = ", Ndepmean, " ± ", (round(Ndepse,3))))

sink()
}