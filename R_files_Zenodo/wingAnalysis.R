#analysis of lab experiments recreating Hoffmann 1973 raising Colias eurytheme at different photoperiods; plus some additional observations

#load data
setwd("~/RDocs/ColiasPostdoc/LabExperiments/Final/DryadSubmission")
oldData<-read.table("./Hoffmann_1973_extracted_data.csv",sep=",",header=T) #data from historic experiment; extracted from Hoffmann 1973, fig 3
spectraData<-read.table("./CA_spectra_gray.csv",sep=",",header=T) #spectra data to relate image data to historic data
imageData<-read.table("./photo_measurement_data.csv",sep=",",header=T) #data from images of butterfly wings; mean value of each color channel for each wing
specimenData<-read.table("./photo_idx.csv",sep=",",header=T) #data about source and treatment of each specimen

#condense spectra data to averages (originally two measurements of each individual)
spectraData<-aggregate(spectraData,list(spectraData$ID),mean)
spectraData$ID<-spectraData$Group.1

#filters out bad data/measures, setting to NA
imageData[which(imageData$colorQualityUL!=""),c(which(colnames(imageData)=="meanULb"):which(colnames(imageData)=="stdDevULr"))]<-NA
imageData[which(imageData$colorQualityUR!=""),c(which(colnames(imageData)=="meanURb"):which(colnames(imageData)=="stdDevURr"))]<-NA
imageData[which(imageData$colorQualityLL!=""),c(which(colnames(imageData)=="meanLLb"):which(colnames(imageData)=="stdDevLLr"))]<-NA
imageData[which(imageData$colorQualityLR!=""),c(which(colnames(imageData)=="meanLRb"):which(colnames(imageData)=="stdDevLRr"))]<-NA
imageData[which(imageData$sizeQualityUL!=""),c(which(colnames(imageData)=="areaUL"):which(colnames(imageData)=="lengthUL"))]<-NA
imageData[which(imageData$sizeQualityUR!=""),c(which(colnames(imageData)=="areaUR"):which(colnames(imageData)=="lengthUR"))]<-NA
imageData[which(imageData$sizeQualityLL!=""),c(which(colnames(imageData)=="areaLL"):which(colnames(imageData)=="lengthLL"))]<-NA
imageData[which(imageData$sizeQualityLR!=""),c(which(colnames(imageData)=="areaLR"):which(colnames(imageData)=="lengthLR"))]<-NA


#rescales length measurements to mm rather than pixels
imageData$scaledLengthLL<-imageData$scale*imageData$lengthLL
imageData$scaledLengthLR<-imageData$scale*imageData$lengthLR
imageData$scaledLengthUL<-imageData$scale*imageData$lengthUL
imageData$scaledLengthUR<-imageData$scale*imageData$lengthUR

#rescales reflectances to proportions; cm to mm
oldData$Reflectance<-oldData$Reflectance/100
oldData$RefMinusSE<-oldData$RefMinusSE/100
oldData$RefPlusSE<-oldData$RefPlusSE/100
oldData$UpperLength<-oldData$UpperLength*10
oldData$ULenMinusSE<-oldData$ULenMinusSE*10
oldData$ULenPlusSE<-oldData$ULenPlusSE*10
oldData$LowerLength<-oldData$LowerLength*10
oldData$LLenMinusSE<-oldData$LLenMinusSE*10
oldData$LLenPlusSE<-oldData$LLenPlusSE*10

#set up additonal variables
specimenData$trueFam<-paste(specimenData$Population,specimenData$Family,sep="") #creates unique family identifier
imageData$meanRG<-rowMeans(subset(imageData,select=c("meanLLg","meanLLr","meanLRg","meanLRr")),na.rm=T) #current melanization metric, averages red and green channels for both hindwings
imageData$meanLengthUp<-rowMeans(subset(imageData,select=c("scaledLengthUL","scaledLengthUR")),na.rm=T) #current size metric, averages length of both forewings

#calculates full SEs for Hoffman's figures
oldData$RefNegSE<-oldData$Reflectance-oldData$RefMinusSE
oldData$RefPosSE<-oldData$RefPlusSE-oldData$Reflectance
oldData$RefSE<-(oldData$RefNegSE+oldData$RefPosSE)/2 #average standard error of old measurements
oldData$ULenNegSE<-oldData$UpperLength-oldData$ULenMinusSE
oldData$ULenPosSE<-oldData$ULenPlusSE-oldData$UpperLength
oldData$ULenSE<-(oldData$ULenNegSE+oldData$ULenPosSE)/2
oldData$LLenNegSE<-oldData$LowerLength-oldData$LLenMinusSE
oldData$LLenPosSE<-oldData$LLenPlusSE-oldData$LowerLength
oldData$LLenSE<-(oldData$LLenNegSE+oldData$LLenPosSE)/2

#associates specimen data with other datasets, using ID
imageData<-merge(imageData,specimenData,by="ID")
spectraData<-merge(spectraData,specimenData,by="ID")

#Relation meanRG to reflectance and estimation of corresponding values; only for males currently (independent effect of sex on regression); using temperature treatment data in correlation
#_____________________________________________________________________________________________________________________________________
ventMaleCAImg<-subset(imageData,side=="V"&Sex=="Male")

ventMaleCAFull<-merge(ventMaleCAImg,spectraData,by="ID")

cor(ventMaleCAFull$X650,ventMaleCAFull$meanRG) #correlation between measures; report in paper

conversion<-lm(X650~meanRG,data=ventMaleCAFull) #male only conversion equation

#generates estimates reflectance at 650 form all images
imageData$X650est<-conversion$coefficients[2]*imageData$meanRG+conversion$coefficients[1]

#Comparing reflectance with white and gray background, to demonstrate that gray is the conservative choice
#____________________________________________________________________________________________________________________________________
whiteSpectraData<-read.table("./CA_spectra_white.csv",sep=",",header=T) #spectra data from measurements on white background

comboSpectra<-merge(spectraData,whiteSpectraData,by="ID")

cor(comboSpectra$X650.x,comboSpectra$X650.y)
range(comboSpectra$X650.y)
range(comboSpectra$X650.x)
plot(comboSpectra$X650.y~comboSpectra$Larval_Photoperiod)

mean(comboSpectra$X650.x)
mean(comboSpectra$X650.y)

reflectSummary<-data.frame(Larval_Photoperiod=c(10,11,12,13,14))
reflectSummary$grayMean<-aggregate(comboSpectra$X650.x,list(comboSpectra$Larval_Photoperiod),mean)$x
reflectSummary$whiteMean<-aggregate(comboSpectra$X650.y,list(comboSpectra$Larval_Photoperiod),mean)$x

#Comparison of old and current data for color
#_____________________________________________________________________________________________________________________________________
ventMaleCAPhotoImg<-subset(imageData,side=="V"&Sex=="Male"&Pupal_Temp==25&!is.nan(meanRG))
ventMaleCAPhotoImgFam<-subset(ventMaleCAPhotoImg, Family!="") #alternative dataset excluding individual without known family

#analysis of just current data
library(lme4)
contemporaryRes<-lmer(X650est~factor(Larval_Photoperiod)+(1|Family), REML=F, data=ventMaleCAPhotoImgFam)
contemporaryResRed<-lmer(X650est~(1|Family), REML=F, data=ventMaleCAPhotoImgFam)
anova(contemporaryRes,contemporaryResRed)

#calcuating means and SE for contemporary data
contempSummary<-data.frame(Larval_Photoperiod=c(10,11,12,13,14))
contempSummary$means<-aggregate(ventMaleCAPhotoImg$X650est,list(ventMaleCAPhotoImg$Larval_Photoperiod),mean)$x
contempSummary$sd<-aggregate(ventMaleCAPhotoImg$X650est,list(ventMaleCAPhotoImg$Larval_Photoperiod),sd)$x
contempSummary$n<-aggregate(ventMaleCAPhotoImg$X650est,list(ventMaleCAPhotoImg$Larval_Photoperiod),length)$x
contempSummary$se<-contempSummary$sd/sqrt(contempSummary$n)

#plot export dim (paper): 750x500 pixels
plot(oldData$Photoperiod,oldData$Reflectance,xlim=c(9,17),ylim=c(0.2,0.5),cex=1.5, cex.axis=1.2,cex.lab=1.3,ylab="Reflectance (at 650 nm)",xlab="Photoperiod (hrs)", bty="l")
arrows(oldData$Photoperiod, oldData$Reflectance-oldData$RefSE*1.96, oldData$Photoperiod, oldData$Reflectance+oldData$RefSE*1.96, length=0.05, angle=90, code=3,lwd=2)
lines(oldData$Photoperiod[5:10],oldData$Reflectance[5:10],lwd=2, lty=2) #range limits range of photoperiods graphed
points(contempSummary$Larval_Photoperiod,contempSummary$means,col="blue",cex=1.5, pch=2)
arrows(contempSummary$Larval_Photoperiod, contempSummary$means-contempSummary$se*1.96, contempSummary$Larval_Photoperiod, contempSummary$means+contempSummary$se*1.96, length=0.05, angle=90, code=3,lwd=2, col="blue")
lines(contempSummary$Larval_Photoperiod,contempSummary$means,col="blue",lwd=2)
legend(15,.3,c(1971,2018),col=c("black","blue"),lwd=2,lty=c(2,1)) #presentation only

#range of reflectances
contempSummary$means[5]-contempSummary$means[1]
oldData$Reflectance[9]-oldData$Reflectance[5]

#stats to compare to hoffman by simiulating data with n=10;
#____________________________________________________________________________________________________________________________________
oldPhotoSimData<-data.frame("Photoperiod"=c(),"Reflectance"=c())

for (i in seq(5,9)) { #loop to make simulated data (based on Larson 1992) with an assummed n=10
  refs<-rep(oldData$Reflectance[i]+oldData$RefSE[i],9)
  refs<-append(refs,10*oldData$Reflectance[i]-9*refs[1])
  oldPhotoSimData<-rbind(oldPhotoSimData,data.frame(Photoperiod=oldData$Photoperiod[i], Reflectance=refs))
}

oldPhotoSimData$Time<-"Past"

fullPhotoData<-rbind(oldPhotoSimData,data.frame(Photoperiod=ventMaleCAPhotoImg$Larval_Photoperiod,Reflectance=ventMaleCAPhotoImg$X650est, Time="Present"))

refResult<-lm(Reflectance~factor(Photoperiod)*Time,data=fullPhotoData)

refResultIndep<-lm(Reflectance~factor(Photoperiod)+Time,data=fullPhotoData)


#Table of stats for comparison with different values of n; uses oldData and ventMaleCAPhotoImg
#___________________________________________________________________________________________________________________________________
simulationTable<-data.frame("n"=c(),"interactF"=c(),"interactP"=c(),"indepF"=c(),"indepP"=c())

for (j in c(2,4,6,8,10,12,16,20,25,30)){
  simData<-data.frame("Photoperiod"=c(),"Reflectance"=c())
  
  for (i in seq(5,9)) { #loop to make simulated data (based on Larson 1992) with an assummed n=10
    refs<-rep(oldData$Reflectance[i]+oldData$RefSE[i],j-1)
    refs<-append(refs,j*oldData$Reflectance[i]-(j-1)*refs[1])
    simData<-rbind(simData,data.frame(Photoperiod=oldData$Photoperiod[i], Reflectance=refs))
  }
  
  simData$Time<-"Past"
  
  simData<-rbind(simData,data.frame(Photoperiod=ventMaleCAPhotoImg$Larval_Photoperiod,Reflectance=ventMaleCAPhotoImg$X650est, Time="Present"))
  
  simResultInteract<-anova(lm(Reflectance~factor(Photoperiod)*Time,data=simData))
  
  simResultIndep<-anova(lm(Reflectance~factor(Photoperiod)+Time,data=simData))
  
  #adds a new line with F values and P values from each model for that iteration of n
  simulationTable<-rbind(simulationTable,data.frame("n"=j,"interactF"=simResultInteract$'F value'[3],"interactP"=simResultInteract$'Pr(>F)'[3],"indepF"=simResultIndep$'F value'[2],"indepP"=simResultIndep$'Pr(>F)'[2]))
}


#alternative method for artificial data set (based on response to Larson 1992); gives same results
#____________________________________________________________________________________________________________________________________
altoldPhotoSimData<-data.frame("Photoperiod"=c(),"Reflectance"=c())

for (i in seq(5,9)) { #loop to make simulated data (based on first response to Larson 1992) with an assummed n=10
  refs<-rep(oldData$Reflectance[i],10-2)
  refs<-append(refs,oldData$Reflectance[i]+oldData$RefSE[i]*sqrt(10)*sqrt((10-1)/2))
  refs<-append(refs,oldData$Reflectance[i]-oldData$RefSE[i]*sqrt(10)*sqrt((10-1)/2))
  altoldPhotoSimData<-rbind(altoldPhotoSimData,data.frame(Photoperiod=oldData$Photoperiod[i], Reflectance=refs))
}

altoldPhotoSimData$Time<-"Past"

altfullPhotoData<-rbind(altoldPhotoSimData,data.frame(Photoperiod=ventMaleCAPhotoImg$Larval_Photoperiod,Reflectance=ventMaleCAPhotoImg$X650est, Time="Present"))

altRefResult<-lm(Reflectance~factor(Photoperiod)*Time,data=altfullPhotoData)

sd(oldPhotoSimData$Reflectance)
sd(altoldPhotoSimData$Reflectance)

#Comparison of old and current temperature data; for use in table
#_____________________________________________________________________________________________________________________________________
ventMaleCATempImg<-subset(imageData,side=="V"&Population=="CA"&Sex=="Male"&Larval_Photoperiod==14&!is.nan(meanRG))

#manual means, s.e./CI
mean20<-mean(subset(ventMaleCATempImg,Pupal_Temp==20)$X650est)
sd20<-sd(subset(ventMaleCATempImg,Pupal_Temp==20)$X650est)
n20<-length(subset(ventMaleCATempImg,Pupal_Temp==20)$X650est) #manual s.e. disagrees slightly with se.fit
mean25<-mean(subset(ventMaleCATempImg,Pupal_Temp==25)$X650est)
sd25<-sd(subset(ventMaleCATempImg,Pupal_Temp==25)$X650est)
n25<-length(subset(ventMaleCATempImg,Pupal_Temp==25)$X650est) #manual s.e. disagrees slightly with se.fit; may need to redo previous section

#comparisons to hoffman data
mean20old<-.440
n20old<-12
sd20old<-.022/1.96*sqrt(n20old)
mean25old<-.455
n25old<-8
sd25old<-.018/1.96*sqrt(n20old)

#simulates dataset from summary statistics (Larson 1992)
refs<-rep(mean20old+sd20old/sqrt(n20old),n20old-1)
refs<-append(refs,n20old*mean20old-(n20old-1)*refs[1])
oldTempSimData<-data.frame(Temp=20,Reflectance=refs)
refs<-rep(mean25old+sd25old/sqrt(n25old),n25old-1)
refs<-append(refs,n25old*mean25old-(n25old-1)*refs[1])
oldTempSimData<-rbind(oldTempSimData,data.frame(Temp=25, Reflectance=refs))

oldTempSimData$Time<-"Past"

fullTempData<-rbind(oldTempSimData,data.frame(Temp=ventMaleCATempImg$Pupal_Temp,Reflectance=ventMaleCATempImg$X650est, Time="Present"))

#anovas comparing current and simulated past data
tempResult<-lm(Reflectance~factor(Temp)*Time,data=fullTempData) #tests two way interaction
tempResult<-lm(Reflectance~factor(Temp)+Time,data=fullTempData) #both orders for type II sums of squares
tempResult<-lm(Reflectance~factor(Temp),data=fullTempData)

#wing length comparison; just upper wings
#___________________________________________________________________________________________________________________________________
ventMaleCAPhotoImg<-subset(imageData,side=="V"&Sex=="Male"&Pupal_Temp==25&!is.nan(meanLengthUp))
ventMaleCAPhotoImgFam<-subset(ventMaleCAPhotoImg, Family!="") #excludes individual without known family

library(lme4)
contemporaryRes<-lmer(meanLengthUp~factor(Larval_Photoperiod)+(1|Family), REML=F, data=ventMaleCAPhotoImgFam)
contemporaryResRed<-lmer(meanLengthUp~(1|Family), REML=F, data=ventMaleCAPhotoImgFam)
anova(contemporaryRes,contemporaryResRed)

#calcuating means and SE for contemporary data
contempSummary<-data.frame(Larval_Photoperiod=c(10,11,12,13,14))
contempSummary$means<-aggregate(ventMaleCAPhotoImg$meanLengthUp,list(ventMaleCAPhotoImg$Larval_Photoperiod),mean)$x
contempSummary$sd<-aggregate(ventMaleCAPhotoImg$meanLengthUp,list(ventMaleCAPhotoImg$Larval_Photoperiod),sd)$x
contempSummary$n<-aggregate(ventMaleCAPhotoImg$meanLengthUp,list(ventMaleCAPhotoImg$Larval_Photoperiod),length)$x
contempSummary$se<-contempSummary$sd/sqrt(contempSummary$n)

#dims (for paper): 750 x 500
plot(oldData$Photoperiod,oldData$UpperLength,xlim=c(9,17),ylim=c(18,24),cex=1.5, cex.axis=1.2,cex.lab=1.3,ylab="Forewing Length (mm)",xlab="Photoperiod (hrs)", bty="l")
arrows(oldData$Photoperiod, oldData$UpperLength-oldData$ULenSE*1.96, oldData$Photoperiod, oldData$UpperLength+oldData$ULenSE*1.96, length=0.05, angle=90, code=3,lwd=2)
lines(oldData$Photoperiod[5:10],oldData$UpperLength[5:10],lwd=2, lty=2) #range limits range of photoperiods graphed
points(contempSummary$Larval_Photoperiod,contempSummary$means,col="blue",cex=1.5,pch=2)
arrows(contempSummary$Larval_Photoperiod, contempSummary$means-contempSummary$se*1.96, contempSummary$Larval_Photoperiod, contempSummary$means+contempSummary$se*1.96, length=0.05, angle=90, code=3,lwd=2, col="blue")
lines(contempSummary$Larval_Photoperiod,contempSummary$means,col="blue",lwd=2)
legend(15,20,c(1971,2018),col=c("black","blue"),lwd=2,lty=c(2,1)) #presentation only

#stats to compare to hoffman by simiulating data with n=10
#____________________________________________________________________________________________________________________________________
oldLenSimData<-data.frame("Photoperiod"=c(),"UpperLength"=c())

for (i in seq(5,9)) { #loop to make simulated data (based on Larson 1992) with an assummed n=10
  refs<-rep(oldData$UpperLength[i]+oldData$ULenSE[i],9)
  refs<-append(refs,10*oldData$UpperLength[i]-9*refs[1])
  oldLenSimData<-rbind(oldLenSimData,data.frame(Photoperiod=oldData$Photoperiod[i], UpperLength=refs))
}

oldLenSimData$Time<-"Past"

fullLenData<-rbind(oldLenSimData,data.frame(Photoperiod=ventMaleCAPhotoImg$Larval_Photoperiod,UpperLength=ventMaleCAPhotoImg$meanLengthUp, Time="Present"))

lenResult<-lm(UpperLength~factor(Photoperiod)*Time,data=fullLenData)
lenResult<-lm(UpperLength~factor(Photoperiod)+Time,data=fullLenData) #without interactions
lenResult<-lm(UpperLength~Time+factor(Photoperiod),data=fullLenData) #without interactions

#tests of family effects (supplemental materials)
#___________________________________________________________________________________________________________________________________
ventMaleCAPhotoImg<-subset(imageData,side=="V"&Sex=="Male"&Pupal_Temp==25&!is.nan(meanRG))
ventMaleCAPhotoImgFam<-subset(ventMaleCAPhotoImg, Family!="") #excludes individual without known family

library(lme4)
contemporaryRes<-lmer(X650est~factor(Larval_Photoperiod)+(1|Family), REML=F, data=ventMaleCAPhotoImgFam)
contemporaryResInteract<-lmer(X650est~factor(Larval_Photoperiod)+(factor(Larval_Photoperiod)|Family), REML=F, data=ventMaleCAPhotoImgFam)
anova(contemporaryRes,contemporaryResInteract)


#plotting family by photoperiod with lines
#___________________________________________________________________________________________________________________________________

#generates family level summary statistics (only mean currently used)
famSummary<-aggregate(ventMaleCAPhotoImgFam$X650est,list(ventMaleCAPhotoImgFam$Family,ventMaleCAPhotoImgFam$Larval_Photoperiod),mean)
famSummary$sd<-aggregate(ventMaleCAPhotoImgFam$X650est,list(ventMaleCAPhotoImgFam$Family,ventMaleCAPhotoImgFam$Larval_Photoperiod),sd)$x
famSummary$n<-aggregate(ventMaleCAPhotoImgFam$X650est,list(ventMaleCAPhotoImgFam$Family,ventMaleCAPhotoImgFam$Larval_Photoperiod),length)$x
famSummary$se<-famSummary$sd/sqrt(famSummary$n)

#plots family means : export 750 x 500 pixels
plot(famSummary$Group.2, famSummary$x,col=famSummary$Group.1, xlim=c(9,15),ylim=c(0.2,0.5),cex=1.5, cex.axis=1.2,cex.lab=1.3,ylab="Reflectance (at 650 nm)",xlab="Photoperiod (hrs)", bty="l")
#arrows(famSummary$Group.2, famSummary$x-famSummary$sd, famSummary$Group.2, famSummary$x+famSummary$sd, length=0.05, angle=90, code=3,lwd=1, col=famSummary$Group.1) #not currently including family level sd; not always meaningful (somtimes only one sample)
lines(famSummary$Group.2[famSummary$Group.1=="A"],famSummary$x[famSummary$Group.1=="A"],col=2,lwd=1.5,lty=2)
lines(famSummary$Group.2[famSummary$Group.1=="B"],famSummary$x[famSummary$Group.1=="B"],col=3,lwd=1.5,lty=2)
lines(famSummary$Group.2[famSummary$Group.1=="C"],famSummary$x[famSummary$Group.1=="C"],col=4,lwd=1.5,lty=2)
lines(famSummary$Group.2[famSummary$Group.1=="D"],famSummary$x[famSummary$Group.1=="D"],col=5,lwd=1.5,lty=2)
lines(famSummary$Group.2[famSummary$Group.1=="E"],famSummary$x[famSummary$Group.1=="E"],col=6,lwd=1.5,lty=2)
lines(famSummary$Group.2[famSummary$Group.1=="F"],famSummary$x[famSummary$Group.1=="F"],col=7,lwd=1.5,lty=2)
lines(famSummary$Group.2[famSummary$Group.1=="G"],famSummary$x[famSummary$Group.1=="G"],col=8,lwd=1.5,lty=2)

#adds contemporary mean and sd to plot
#note that this redoes contempSummary without individual of unknown family
contempSummary<-data.frame(Larval_Photoperiod=c(10,11,12,13,14))
contempSummary$means<-aggregate(ventMaleCAPhotoImgFam$X650est,list(ventMaleCAPhotoImgFam$Larval_Photoperiod),mean)$x 
contempSummary$sd<-aggregate(ventMaleCAPhotoImgFam$X650est,list(ventMaleCAPhotoImgFam$Larval_Photoperiod),sd)$x

points(contempSummary$Larval_Photoperiod,contempSummary$means, cex=2, pch=2)
arrows(contempSummary$Larval_Photoperiod, contempSummary$means-contempSummary$sd, contempSummary$Larval_Photoperiod, contempSummary$means+contempSummary$sd, length=0.05, angle=90, code=3,lwd=2, col="black")
lines(contempSummary$Larval_Photoperiod,contempSummary$means,col="black",lwd=2)

#investigation of sex effect in current population
#____________________________________________________________________________________________________________________________________
ventCAImg<-subset(imageData,side=="V") #uses all CA data
ventCAFull<-merge(ventCAImg,spectraData,by="ID")

#conversion to reflectance at 650 nm
fullConversion<-lm(X650~meanRG+Sex.x, data=ventCAFull)
imageData$X650est2<-fullConversion$coefficients[2]*imageData$meanRG+fullConversion$coefficients[1]+fullConversion$coefficients[3]*(imageData$Sex=="Male")

library(lme4)

ventCAPhotoImg<-subset(imageData,side=="V"&Population=="CA"&Pupal_Temp==25)
ventCAPhotoImgFam<-subset(ventCAPhotoImg, Family!="") #drops individual without family

fullRes<-lmer(X650est2~factor(Larval_Photoperiod)*Sex+(1|Family), REML=F, data=ventCAPhotoImgFam)

lessRes<-lmer(X650est2~factor(Larval_Photoperiod)+Sex+(1|Family), REML=F, data=ventCAPhotoImgFam)
anova(fullRes,lessRes)

fullSizeRes<-lmer(meanLengthUp~factor(Larval_Photoperiod)*Sex+(1|Family), REML=F, data=ventCAPhotoImgFam)

lessSizeRes<-lmer(meanLengthUp~factor(Larval_Photoperiod)+Sex+(1|Family), REML=F, data=ventCAPhotoImgFam)
anova(fullSizeRes,lessSizeRes)

lesslessSizeRes<-lmer(meanLengthUp~factor(Larval_Photoperiod)+(1|Family), REML=F, data=ventCAPhotoImgFam)
anova(lesslessSizeRes,lessSizeRes)

otherlesslessSizeRes<-lmer(meanLengthUp~Sex+(1|Family), REML=F, data=ventCAPhotoImgFam)
anova(otherlesslessSizeRes,lessSizeRes)

#plots for supplemental materials
#____________________________________________________________________________________________________________________________________

#calculates summary data
sexSummary<-data.frame(Larval_Photoperiod=c(10,11,12,13,14), Sex=c(rep("Female",5),rep("Male",5)))

ventCAPhotoImg<-subset(imageData,side=="V"&Population=="CA"&Pupal_Temp==25&!is.nan(X650est2)) #excludes NANs for reflectance
ventCAPhotoImgFam<-subset(ventCAPhotoImg, Family!="") #drops individual without family

sexSummary$refMeans<-aggregate(ventCAPhotoImgFam$X650est2,list(ventCAPhotoImgFam$Larval_Photoperiod,ventCAPhotoImgFam$Sex),mean)$x
sexSummary$refsd<-aggregate(ventCAPhotoImgFam$X650est2,list(ventCAPhotoImgFam$Larval_Photoperiod,ventCAPhotoImgFam$Sex),sd)$x
sexSummary$refn<-aggregate(ventCAPhotoImgFam$X650est2,list(ventCAPhotoImgFam$Larval_Photoperiod,ventCAPhotoImgFam$Sex),length)$x
sexSummary$refse<-sexSummary$refsd/sqrt(sexSummary$refn)

ventCAPhotoImg<-subset(imageData,side=="V"&Population=="CA"&Pupal_Temp==25&!is.nan(meanLengthUp)) #excludes NANs for length
ventCAPhotoImgFam<-subset(ventCAPhotoImg, Family!="") #drops individual without family
sexSummary$lenMeans<-aggregate(ventCAPhotoImgFam$meanLengthUp,list(ventCAPhotoImgFam$Larval_Photoperiod,ventCAPhotoImgFam$Sex),mean)$x
sexSummary$lensd<-aggregate(ventCAPhotoImgFam$meanLengthUp,list(ventCAPhotoImgFam$Larval_Photoperiod,ventCAPhotoImgFam$Sex),sd)$x
sexSummary$lenn<-aggregate(ventCAPhotoImgFam$meanLengthUp,list(ventCAPhotoImgFam$Larval_Photoperiod,ventCAPhotoImgFam$Sex),length)$x
sexSummary$lense<-sexSummary$lensd/sqrt(sexSummary$lenn)

femaleSummary<-subset(sexSummary,Sex=="Female")
maleSummary<-subset(sexSummary,Sex=="Male")

#reflectance plot : 750x500 pixels
plot(femaleSummary$Larval_Photoperiod,femaleSummary$refMeans,xlim=c(9,17),ylim=c(0.2,0.5),cex=1.5, cex.axis=1.2,cex.lab=1.3,ylab="Reflectance (at 650 nm)",xlab="Photoperiod (hrs)", bty="l")
arrows(femaleSummary$Larval_Photoperiod, femaleSummary$refMeans-femaleSummary$refse*1.96, femaleSummary$Larval_Photoperiod, femaleSummary$refMeans+femaleSummary$refse*1.96, length=0.05, angle=90, code=3,lwd=2)
lines(femaleSummary$Larval_Photoperiod,femaleSummary$refMeans,lwd=2, lty=2) #range limits range of photoperiods graphed
points(maleSummary$Larval_Photoperiod,maleSummary$refMeans,col="blue",cex=1.5, pch=2)
arrows(maleSummary$Larval_Photoperiod, maleSummary$refMeans-maleSummary$refse*1.96, maleSummary$Larval_Photoperiod, maleSummary$refMeans+maleSummary$refse*1.96, length=0.05, angle=90, code=3,lwd=2, col="blue")
lines(maleSummary$Larval_Photoperiod,maleSummary$refMeans,col="blue",lwd=2)

#size plot : 750x500 pixels
plot(femaleSummary$Larval_Photoperiod,femaleSummary$lenMeans,xlim=c(9,17),ylim=c(18,26),cex=1.5, cex.axis=1.2,cex.lab=1.3,ylab="Forewing Length (mm)",xlab="Photoperiod (hrs)", bty="l")
arrows(femaleSummary$Larval_Photoperiod, femaleSummary$lenMeans-femaleSummary$lense*1.96, femaleSummary$Larval_Photoperiod, femaleSummary$lenMeans+femaleSummary$lense*1.96, length=0.05, angle=90, code=3,lwd=2)
lines(femaleSummary$Larval_Photoperiod,femaleSummary$lenMeans,lwd=2, lty=2) #range limits range of photoperiods graphed
points(maleSummary$Larval_Photoperiod,maleSummary$lenMeans,col="blue",cex=1.5, pch=2)
arrows(maleSummary$Larval_Photoperiod, maleSummary$lenMeans-maleSummary$lense*1.96, maleSummary$Larval_Photoperiod, maleSummary$lenMeans+maleSummary$lense*1.96, length=0.05, angle=90, code=3,lwd=2, col="blue")
lines(maleSummary$Larval_Photoperiod,maleSummary$lenMeans,col="blue",lwd=2)
