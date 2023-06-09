setwd("~/Documents/PhD Project/Respirometry/Respirometry/Cobia Corrected Data")
library(StreamMetabolism)

#read in all csvs ending with "data_corrected.csv" and putting csvs in list
mycsv = dir(pattern="data_corrected.csv")
n <- length(mycsv)
mylist <- vector("list", n)
for(i in 1:n) mylist[[i]] <- read.csv(mycsv[i])


sumdat<-data.frame()

for(q in 1:length(mylist)){
  dat<-mylist[q]
  dat<-data.frame(dat)
  dat$Date_Time<-paste(dat$Date, dat$time, seq=" ")
  dat$Date_Time<-as.POSIXct(dat$Date_Time, format="%m/%d/%y %H:%M:%S")
  dat<-dat[which(dat$r2>=0.80),] #gets rid of points where r2 is lower than .80
  dat$Latitude<- +37.2501
  dat$Longitude<- -76.4983
  dat$Date <- format(dat$Date_Time, format="%Y/%m/%d")
  
  for(i in 1:nrow(dat)){
    #incorporates DST!
    sun<-sunrise.set(dat$Latitude[i],dat$Longitude[i], dat$Date[i], timezone = "America/New_York",num.days=1)
    dat$Sunrise[i]<-as.character(sun$sunrise)
    dat$Sunset[i]<-as.character(sun$sunset)
  }
  dat$Sunrise<-as.POSIXct(as.factor(dat$Sunrise), format="%Y-%m-%d %H:%M:%S")
  dat$Sunset<-as.POSIXct(dat$Sunset, format="%Y-%m-%d %H:%M:%S")
  #compares date/time of measurement and sunrise and sunset of measurement
  dat$Sunrise1<-ifelse(as.numeric(dat$Sunrise)<as.numeric(dat$Date_Time),"earlier",dat$Sunrise) 
  dat$Sunset1<-ifelse(as.numeric(dat$Sunset)>as.numeric(dat$Date_Time),"later",dat$Sunset) 
  #assigns day/night for each measurement
  dat$Day_Night<-ifelse(dat$Sunrise1=="earlier" & dat$Sunset1=="later","Day","Night")
  
  #some trials, it goes over two nights so might need to specify night 1 or 2
  dates<-unique(dat$Date)
  dat$DayNum<-ifelse(dat$Date==dates[1],1,ifelse(dat$Date==dates[2],2,3))
  dat$Hour<-format(dat$Date_Time, "%H")
  for(i in 1:nrow(dat)){
    blah<-ifelse(dat$Day_Night[i]=="Night" & dat$Date[i]!=dat$Date[1], ifelse(as.numeric(dat$Hour[i])>12 | dat$DayNum[i]>2 ,"Night2",dat$Day_Night[i]),dat$Day_Night[i])
    dat$blah[i]<-blah
  }
  dat$Day_Night<-dat$blah
  
  #change O2 saturation based on temp and salinity (look at O2 saturation table)
  dat$Percent_Max_O2<-round(dat$max..O2./dat$X100_Perc_O2, digit=2)
  dat$Percent_Min_O2<-round(dat$min..O2./dat$X100_Perc_O2, digit=2) #just to have percent of min O2
  #if O2 percent is less than 80% and after start of hypoxia trial is called it hypoxic
  dat$hypoxia<-NA
  hypoxic_loc<-which(dat$Comments=="hypoxia")
  hypoxic_loc<-ifelse(length(hypoxic_loc)==0,0,hypoxic_loc)
  dat$hypoxia<-ifelse(dat$Percent_Max_O2<0.80 & dat$X>=hypoxic_loc,"yes","no")
  badloc<-which((dat$hypoxia=="no" & dat$Percent_Max_O2<0.80) | dat$finalMO2<0) #gets rid of rows where water became hypoxic during normoxic part of trial
  
  if(length(badloc)==0){
    dat=dat
  }else (dat=dat[-badloc,])  
  
  numdat<-nrow(dat)
  
  norm<-dat[-which(dat$hypoxia=="yes"),]
  hypox<-dat[which(dat$hypoxia=="yes"),]
  
  #plot entire trial, remember not all trials have two nights so will see warnings for some trials
  setwd("~/Documents/PhD Project/Respirometry/Respirometry/Trials")
  png(file=paste(dat$Animal_ID[1], dat$Test_Temp[1],"png", sep="."),width=650,height=450)
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,4))
  plot(finalMO2~Date_Timen, main=paste(dat$Animal_ID[1], dat$Test_Temp[1],sep="_") ,data=dat, xaxt="n", type="b", ylab="Metabolic rate (mg/kg/h)", xlab="")
  r<-as.POSIXct(round(range(dat$Date_Time), "hours"))
  axis.POSIXct(1, at = seq(r[1], r[2], by = 7200), format = "%H:%M",las=2)
  mtext(side=1,line=3.75,"Time of day")
  polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10), #for hypoxia trial
          x=c(dat$Date_Timen[min(which(dat$hypoxia=="yes"))],dat$Date_Timen[min(which(dat$hypoxia=="yes"))],
              dat$Date_Timen[max(which(dat$hypoxia=="yes"))],dat$Date_Timen[max(which(dat$hypoxia=="yes"))]),col=rgb(1,0,0,.5))
  polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10),  #for night shading
          x=c(dat$Sunset1[min(which(dat$Day_Night=="Night"))],dat$Sunset1[min(which(dat$Day_Night=="Night"))],
              dat$Sunrise1[max(which(dat$Day_Night=="Night"))],dat$Sunrise1[max(which(dat$Day_Night=="Night"))]),col=rgb(.5,.5,.5,.5))
  polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10),  #for night shading
          x=c(dat$Sunset1[min(which(dat$Day_Night=="Night2"))],dat$Sunset1[min(which(dat$Day_Night=="Night2"))],
              dat$Date_Timen[max(which(dat$Day_Night=="Night2"))],dat$Date_Timen[max(which(dat$Day_Night=="Night2"))]),col=rgb(.5,.5,.5,.5))
  par(new=T)
  plot(Percent_Max_O2*100~Date_Timen, data=dat, type="l", xaxt="n", axes=F, xlab="",ylab="",lwd=3,col="red")
  axis(side=4)
  mtext(side=4,line=2.5,"Oxygen saturation (%)")
  dev.off()
  
  datAS<-subset(dat, select=-c(Sunrise,Sunset,Sunrise1,Sunset1, Latitude, Longitude,blah,Hour))
  
  datAS<-datAS[which(datAS$hypoxia=="no"),]
  
  d<-hist(datAS$finalMO2, breaks=5)
  
  
  ###############################################################
  # Best Method
  # MMR: highest point
  # RMR: mean lowest 10% during normoxia
  ###############################################################
  
  #MMR, highest MO2 during first 3 hours
  chase_start<-which(datAS$Comments=="chase")
  datAS$chase<-ifelse(datAS$Date_Time<datAS$Date_Time[chase_start]+10800 & datAS$Date_Time>=datAS$Date_Time[chase_start],"yes","no")
  
  MMR<-max(datAS$finalMO2[which(datAS$chase=="yes")])
  MMR
  
  #SMR, mean of 10% of lowest points during normoxia
  #pull out lowest 10% points
  n<-10
  low10<-datAS[datAS$finalMO2 < quantile(datAS$finalMO2,prob=n/100),]
  nrowlow10<-nrow(low10)
  #mean of these points
  SMR<-mean(low10$finalMO2)
  SMRsd<-sd(low10$finalMO2)
  SMRse<-SMRsd/sqrt(nrow(low10))
  SMR
  
  AS<-MMR-SMR
  AS
  
  FS<-MMR/SMR
  FS
  
  
  ##########
  #Scrit
  #finds first value where MO2 drops below SMR and the remaining points are below SMR
  brec<-NULL
  hypox<-hypox[order(hypox$min..O2., decreasing = TRUE), ] #needed to do this because C23_32 decreased during Scrit strangely, this did not change any other Scrit or Ccrit values
  for(i in 1:nrow(hypox)){
    b<-ifelse(any(hypox$finalMO2[i:nrow(hypox)]>SMR)=="TRUE","bad","good") 
    brec[i]<-b
  }
  min<-hypox$min..O2.[min(which(brec=="good"))]
  
  #finds first value where MO2 drops below SMR and the remaining points drop below the previous
  #brec<-NULL
  #for(i in 1:nrow(hypox)){
  #  b<-ifelse(any(hypox$finalMO2[i:nrow(hypox)]>SMR)=="TRUE","bad",ifelse(all(hypox$finalMO2[(i+1):nrow(hypox)]<hypox$finalMO2[i]),"good","bad")) 
  #  brec[i]<-b
  #}
  #brec[length(brec)]<-"good"
  #min<-hypox$min..O2.[min(which(brec=="good"))]
  
  #make regression line through points where MO2 drops below SMR and the remaining points are below SMR
  regpts<-hypox[which(hypox$min..O2.<=min),]
  if(is.na(min)==FALSE){ #if statement for when animal doesn't make it to hypoxia trial
    reglm<-lm(finalMO2~min..O2., data=regpts)
    slope<-coef(reglm)[2]
    intercept<-coef(reglm)[1]
    xrange<-seq(from = min(dat$min..O2.), to= max(dat$min..O2.),by = 5)
    ypts<-(slope*xrange)+intercept
    lines(ypts~xrange)
    ccrit<-as.numeric((SMR-intercept)/slope) #solve x for when y = SMR
  }else (ccrit<-NA)
  
  #if statement because there are times when the slope is neg not pos or ccrit is NA, this messes up Scrit and Ccrit
  if(is.na(ccrit)==TRUE){
    hypoxxrange<-NA
  }else if(min(hypox$min..O2.)<ccrit){
    hypoxxrange<-seq(from = min(hypox$min..O2.), to= ccrit,by = 0.1)
  }else(hypoxxrange<-NA)
  
  hypoxypts<-(slope*hypoxxrange)+intercept
  
  #if statement because there are times when the ccrit is NA
  if(is.na(ccrit)==TRUE){
    SMRxrange<-NA
  }else if(max(dat$min..O2.)>ccrit){
    SMRxrange<-seq(from = ccrit, to= max(dat$min..O2.),by = 1)
  }else(SMRxrange<-NA)
  
  SMRypts<-rep(SMR,length(SMRxrange))
  
  #plot hypoxia trial
  setwd("~/Documents/PhD Project/Respirometry/Respirometry/Hypoxia Plots")
  png(file=paste(dat$Animal_ID[1], dat$Test_Temp[1],"png", sep="."),width=650,height=450)
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2))
  plot(finalMO2~min..O2., data=dat, pch=16, xlab="Min O2 (uM)", ylab="Metabolic rate (mg/kg/h)",main=paste(dat$Animal_ID[1], dat$Test_Temp[1], sep="_"))
  lines(hypoxypts~hypoxxrange, col="red")
  lines(SMRypts~SMRxrange, col="red")
  scrit<-round(ccrit/dat$X100_Perc_O2[1],3)*100
  scrit<-ifelse(hypoxxrange[1]=="NA",NA,scrit)
  scrit
  ccrit<-ifelse(hypoxxrange[1]=="NA",NA,ccrit)
  new_ccrit<-ccrit*32/1000
  new_ccrit
  text(x=(.7*max(dat$min..O2.)), y=(0.93*max(dat$finalMO2)),paste(scrit,"%",sep=""), cex=2, col="red")
  dev.off()
  
  #maybe if I want to determine scrit based on lower end of CI later
  #SMRnegCI<-SMR-(sd(norm$finalMO2[1:10])*1.96)
  #SMRposCI<-SMR+(sd(norm$finalMO2[1:10])*1.96)
  
  #pred<-predict(reglm, se.fit=T)
  #regpts$fit<-pred$fit
  #regpts$se<-pred$se.fit
  #regpts$lowerCI<-regpts$fit-(regpts$se*1.96)
  #regpts$upperCI<-regpts$fit+(regpts$se*1.96)
  
  #lines(fit~min..O2., data=regpts)

  
  animal<-dat$Animal_ID[1]
  test_temp<-dat$Test_Temp[1]
  animal_mass<-dat$animal.mass[1]
  
  sumdat[q,1]<-levels(animal)
  sumdat[q,2]<-test_temp
  sumdat[q,3]<-MMR
  sumdat[q,4]<-SMR
  sumdat[q,5]<-AS
  sumdat[q,6]<-FS
  sumdat[q,7]<-scrit
  sumdat[q,8]<-new_ccrit
  sumdat[q,9]<-animal_mass
  sumdat[q,10]<-nrowlow15
  sumdat[q,11]<-numdat
  
  
  print(q)
}

colnames(sumdat)<-c("Animal_ID","Test_Temp","MMR","SMR","AS","FS","Scrit","Ccrit","Animal_Mass","NumSMR","Numdat")


#write csv for sumdat











#### OTHER METHODS


###############################################################
# Method 1: 3 high 10 low approach
# MMR: average top 3 within first 8 hours
# RMR: average of 10 lowest over entire normoxic trial
# suffix is 3_10
###############################################################
#get rid of outliers (MO2 < 2 sd from mean)
#range<-c((avg+(2*sd)), (avg-(2*sd)))
#dat<-dat[-which(dat$finalMO2<range[2]),]
datAS<-datAS[with(datAS,order(finalMO2), decreasing=TRUE),] 


#SMR - took mean of lowest 10 values
SMR<-mean(datAS$finalMO2[1:10])
SMRsd<-sd(datAS$finalMO2[1:10])
SMRse<-SMRsd/sqrt(length(datAS$finalMO2[1:10]))
SMR3_10<-SMR

#MMR, mean of 3 hightest MO2 within 8 hours of chase
datAS<-datAS[with(datAS,order(Date_Time)),]
chase_start<-which(datAS$Comments=="chase")
datAS$chase<-ifelse(datAS$Date_Time<datAS$Date_Time[chase_start]+28800 & datAS$Date_Time>=datAS$Date_Time[chase_start],"yes","no")
chase_pts<-datAS[which(datAS$Date_Time<datAS$Date_Time[chase_start]+28800 & datAS$Date_Time>=datAS$Date_Time[chase_start]),]
chase_pts<-chase_pts[with(chase_pts,order(finalMO2)),] 
MMR<-mean(chase_pts$finalMO2[nrow(chase_pts):(nrow(chase_pts)-2)])
MMRsd<-sd(chase_pts$finalMO2[nrow(chase_pts):(nrow(chase_pts)-2)])
MMRse<-MMRsd/sqrt(length(chase_pts$finalMO2[nrow(chase_pts):(nrow(chase_pts)-2)]))
MMR3_10<-MMR

AS3_10<-MMR3_10-SMR3_10
AS3_10

##########
#Scrit
#finds first value where MO2 drops below SMR and the remaining points are below SMR
brec<-NULL
for(i in 1:nrow(hypox)){
  b<-ifelse(any(hypox$finalMO2[i:nrow(hypox)]>SMR)=="TRUE","bad","good") 
  brec[i]<-b
}
min<-hypox$min..O2.[min(which(brec=="good"))]

#make regression line through points where MO2 drops below SMR and the remaining points are below SMR
regpts<-hypox[which(hypox$min..O2.<=min),]
if(is.na(min)==FALSE){ #if statement for when animal doesn't make it to hypoxia trial
  reglm<-lm(finalMO2~min..O2., data=regpts)
  slope<-coef(reglm)[2]
  intercept<-coef(reglm)[1]
  xrange<-seq(from = min(dat$min..O2.), to= max(dat$min..O2.),by = 5)
  ypts<-(slope*xrange)+intercept
  lines(ypts~xrange)
  ccrit<-as.numeric((SMR-intercept)/slope) #solve x for when y = SMR
}else (ccrit<-NA)

#if statement because there are times when the slope is neg not pos or ccrit is NA, this messes up Scrit and Ccrit
if(is.na(ccrit)==TRUE){
  hypoxxrange<-NA
}else if(min(hypox$min..O2.)<ccrit){
  hypoxxrange<-seq(from = min(hypox$min..O2.), to= ccrit,by = 0.1)
}else(hypoxxrange<-NA)

hypoxypts<-(slope*hypoxxrange)+intercept

#if statement because there are times when the ccrit is NA
if(is.na(ccrit)==TRUE){
  SMRxrange<-NA
}else if(max(dat$min..O2.)>ccrit){
  SMRxrange<-seq(from = ccrit, to= max(dat$min..O2.),by = 1)
}else(SMRxrange<-NA)

SMRypts<-rep(SMR,length(SMRxrange))

plot(finalMO2~min..O2., data=dat, pch=16, xlab="Min O2 (uM)", ylab="Metabolic rate (mg/kg/h)",main=paste(dat$Animal_ID[1], dat$Test_Temp[1],"3_10" , sep="_"))
lines(hypoxypts~hypoxxrange, col="red")
lines(SMRypts~SMRxrange, col="red")
text(x=(.7*max(dat$min..O2.)), y=(0.93*max(dat$finalMO2)),round(ccrit/dat$X100_Perc_O2[1],3)*100, cex=2, col="red")

scrit<-round(ccrit/dat$X100_Perc_O2[1],3)*100
scrit<-ifelse(hypoxxrange[1]=="NA",NA,scrit)
scrit_3_10<-scrit
ccrit<-ifelse(hypoxxrange[1]=="NA",NA,ccrit)
new_ccrit<-ccrit*32/1000
new_ccrit_3_10<-new_ccrit

#maybe if I want to determine scrit based on lower end of CI later
#SMRnegCI<-SMR-(sd(norm$finalMO2[1:10])*1.96)
#SMRposCI<-SMR+(sd(norm$finalMO2[1:10])*1.96)

#pred<-predict(reglm, se.fit=T)
#regpts$fit<-pred$fit
#regpts$se<-pred$se.fit
#regpts$lowerCI<-regpts$fit-(regpts$se*1.96)
#regpts$upperCI<-regpts$fit+(regpts$se*1.96)

#lines(fit~min..O2., data=regpts)



###############################################################
# Method 2: 3 high 3 low approach
# MMR: average top 3 within first 8 hours
# RMR: average of 3 lowest over entire normoxic trial
# suffix is 3_3
###############################################################
#get rid of outliers (MO2 < 2 sd from mean)
#range<-c((avg+(2*sd)), (avg-(2*sd)))
#dat<-dat[-which(dat$finalMO2<range[2]),]
datAS<-datAS[with(datAS,order(finalMO2), decreasing=TRUE),] 


#SMR - took mean of lowest 3 values
SMR<-mean(datAS$finalMO2[1:3])
SMRsd<-sd(datAS$finalMO2[1:3])
SMRse<-SMRsd/sqrt(length(datAS$finalMO2[1:3]))
SMR3_3<-SMR

#MMR, mean of 3 hightest MO2 within 8 hours of chase
datAS<-datAS[with(datAS,order(Date_Time)),]
chase_start<-which(datAS$Comments=="chase")
datAS$chase<-ifelse(datAS$Date_Time<datAS$Date_Time[chase_start]+28800 & datAS$Date_Time>=datAS$Date_Time[chase_start],"yes","no")
chase_pts<-datAS[which(datAS$Date_Time<datAS$Date_Time[chase_start]+28800 & datAS$Date_Time>=datAS$Date_Time[chase_start]),]
chase_pts<-chase_pts[with(chase_pts,order(finalMO2)),] 
MMR<-mean(chase_pts$finalMO2[nrow(chase_pts):(nrow(chase_pts)-2)])
MMRsd<-sd(chase_pts$finalMO2[nrow(chase_pts):(nrow(chase_pts)-2)])
MMRse<-MMRsd/sqrt(length(chase_pts$finalMO2[nrow(chase_pts):(nrow(chase_pts)-2)]))
MMR3_3<-MMR

AS3_3<-MMR3_3-SMR3_3
AS3_3

##########
#Scrit
#finds first value where MO2 drops below SMR and the remaining points are below SMR
brec<-NULL
for(i in 1:nrow(hypox)){
  b<-ifelse(any(hypox$finalMO2[i:nrow(hypox)]>SMR)=="TRUE","bad","good") 
  brec[i]<-b
}
min<-hypox$min..O2.[min(which(brec=="good"))]

#make regression line through points where MO2 drops below SMR and the remaining points are below SMR
regpts<-hypox[which(hypox$min..O2.<=min),]
if(is.na(min)==FALSE){ #if statement for when animal doesn't make it to hypoxia trial
  reglm<-lm(finalMO2~min..O2., data=regpts)
  slope<-coef(reglm)[2]
  intercept<-coef(reglm)[1]
  xrange<-seq(from = min(dat$min..O2.), to= max(dat$min..O2.),by = 5)
  ypts<-(slope*xrange)+intercept
  lines(ypts~xrange)
  ccrit<-as.numeric((SMR-intercept)/slope) #solve x for when y = SMR
}else (ccrit<-NA)

#if statement because there are times when the slope is neg not pos or ccrit is NA, this messes up Scrit and Ccrit
if(is.na(ccrit)==TRUE){
  hypoxxrange<-NA
}else if(min(hypox$min..O2.)<ccrit){
  hypoxxrange<-seq(from = min(hypox$min..O2.), to= ccrit,by = 0.1)
}else(hypoxxrange<-NA)

hypoxypts<-(slope*hypoxxrange)+intercept

#if statement because there are times when the ccrit is NA
if(is.na(ccrit)==TRUE){
  SMRxrange<-NA
}else if(max(dat$min..O2.)>ccrit){
  SMRxrange<-seq(from = ccrit, to= max(dat$min..O2.),by = 1)
}else(SMRxrange<-NA)

SMRypts<-rep(SMR,length(SMRxrange))

plot(finalMO2~min..O2., data=dat, pch=16, xlab="Min O2 (uM)", ylab="Metabolic rate (mg/kg/h)",main=paste(dat$Animal_ID[1], dat$Test_Temp[1],"3_3" , sep="_"))
lines(hypoxypts~hypoxxrange)
lines(SMRypts~SMRxrange)
text(x=(.7*max(dat$min..O2.)), y=(0.93*max(dat$finalMO2)),round(ccrit/dat$X100_Perc_O2[1],3)*100, cex=2)

scrit<-round(ccrit/dat$X100_Perc_O2[1],3)*100
scrit<-ifelse(hypoxxrange[1]=="NA",NA,scrit)
scrit_3_3<-scrit
ccrit<-ifelse(hypoxxrange[1]=="NA",NA,ccrit)
new_ccrit<-ccrit*32/1000
new_ccrit_3_3<-new_ccrit

###############################################################
# Method 3: GAM
# MMR: use highest prediction from GAM model
# RMR: use lowest prediction from GAM model
# suffix is GAM
###############################################################

datGAM<-datAS
datGAM$Day_Night<-as.factor(datGAM$Day_Night)
datGAM$chase<-as.factor(datGAM$chase)

library(mgcv)


#predGAM<-predict(g1, newdata = newdat,type="response",se.fit=TRUE)
#d<-cbind(newdat,predGAM$fit[1:nrow(predGAM$fit)],predGAM$se[1:nrow(predGAM$se)])
#names(d)[4:5]<-c("Averages","Variance")
#results<-summaryBy(Averages~Date_Timen, data=d, FUN=c(mean))# calculate groupwise summary statistics (mean) across area
#vari<-summaryBy(Averages~Date_Timen, data=d, FUN=c(var))
#sd<-sqrt(vari$Averages.var) #gets sd from variances 

g1<-gam(finalMO2~s(Date_Timen), data=datGAM, family=gaussian)
summary(g1)
gam.check(g1)
newdat<-data.frame(Date_Timen=seq(from=min(datGAM$Date_Timen), to=max(datGAM$Date_Timen), by=1000))
newdat$predGAM<-predict(g1, newdata = newdat,type="response")


par(mar=c(5,4,4,4))
plot(finalMO2~Date_Timen, main=paste(dat$Animal_ID[1], dat$Test_Temp[1],"GAM" , sep="_"),data=dat, xaxt="n", type="b", ylab="Metabolic rate (mg/kg/h)", xlab="")
r<-as.POSIXct(round(range(dat$Date_Time), "hours"))
axis.POSIXct(1, at = seq(r[1], r[2], by = 7200), format = "%H:%M",las=2)
mtext(side=1,line=3.75,"Time of day")
polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10), #for hypoxia trial
        x=c(dat$Date_Timen[min(which(dat$hypoxia=="yes"))],dat$Date_Timen[min(which(dat$hypoxia=="yes"))],
            dat$Date_Timen[max(which(dat$hypoxia=="yes"))],dat$Date_Timen[max(which(dat$hypoxia=="yes"))]),col=rgb(1,0,0,.5))
polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10),  #for night shading
        x=c(dat$Sunset1[min(which(dat$Day_Night=="Night"))],dat$Sunset1[min(which(dat$Day_Night=="Night"))],
            dat$Sunrise1[max(which(dat$Day_Night=="Night"))],dat$Sunrise1[max(which(dat$Day_Night=="Night"))]),col=rgb(.5,.5,.5,.5))
polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10),  #for night shading
        x=c(dat$Sunset1[min(which(dat$Day_Night=="Night2"))],dat$Sunset1[min(which(dat$Day_Night=="Night2"))],
            dat$Date_Timen[max(which(dat$Day_Night=="Night2"))],dat$Date_Timen[max(which(dat$Day_Night=="Night2"))]),col=rgb(.5,.5,.5,.5))
lines(predGAM~Date_Timen,dat=newdat, type="b", col="blue")
par(new=T)
plot(Percent_Max_O2*100~Date_Timen, data=dat, type="l", xaxt="n", axes=F, xlab="",ylab="",lwd=3,col="red")
axis(side=4)
mtext(side=4,line=2.5,"Oxygen saturation (%)")

max_time<-newdat$Date_Timen[which.max(newdat$predGAM)]
chase_max<-max(chase_pts$Date_Timen)
chase_min<-min(chase_pts$Date_Timen)
MMR_GAM<-ifelse(max_time<=chase_max & max_time>=chase_min,max(newdat$predGAM, na.rm=T),NA)
MMR_GAM
SMR_GAM<-min(newdat$predGAM)
SMR_GAM
AS_GAM<-MMR_GAM-SMR_GAM
AS_GAM

##########
#Scrit
#finds first value where MO2 drops below SMR and the remaining points are below SMR
SMR<-SMR_GAM
brec<-NULL
for(i in 1:nrow(hypox)){
  b<-ifelse(any(hypox$finalMO2[i:nrow(hypox)]>SMR)=="TRUE","bad","good") 
  brec[i]<-b
}
min<-hypox$min..O2.[min(which(brec=="good"))]

#make regression line through points where MO2 drops below SMR and the remaining points are below SMR
regpts<-hypox[which(hypox$min..O2.<=min),]
if(is.na(min)==FALSE){ #if statement for when animal doesn't make it to hypoxia trial
  reglm<-lm(finalMO2~min..O2., data=regpts)
  slope<-coef(reglm)[2]
  intercept<-coef(reglm)[1]
  xrange<-seq(from = min(dat$min..O2.), to= max(dat$min..O2.),by = 5)
  ypts<-(slope*xrange)+intercept
  lines(ypts~xrange)
  ccrit<-as.numeric((SMR-intercept)/slope) #solve x for when y = SMR
}else (ccrit<-NA)

#if statement because there are times when the slope is neg not pos or ccrit is NA, this messes up Scrit and Ccrit
if(is.na(ccrit)==TRUE){
  hypoxxrange<-NA
}else if(min(hypox$min..O2.)<ccrit){
  hypoxxrange<-seq(from = min(hypox$min..O2.), to= ccrit,by = 0.1)
}else(hypoxxrange<-NA)

hypoxypts<-(slope*hypoxxrange)+intercept

#if statement because there are times when the ccrit is NA
if(is.na(ccrit)==TRUE){
  SMRxrange<-NA
}else if(max(dat$min..O2.)>ccrit){
  SMRxrange<-seq(from = ccrit, to= max(dat$min..O2.),by = 1)
}else(SMRxrange<-NA)

SMRypts<-rep(SMR,length(SMRxrange))

plot(finalMO2~min..O2., data=dat, pch=16, xlab="Min O2 (uM)", ylab="Metabolic rate (mg/kg/h)",main=paste(dat$Animal_ID[1], dat$Test_Temp[1],"GAM" , sep="_"))
lines(hypoxypts~hypoxxrange, col="blue")
lines(SMRypts~SMRxrange, col="blue")
text(x=(.7*max(dat$min..O2.)), y=(0.93*max(dat$finalMO2)),round(ccrit/dat$X100_Perc_O2[1],3)*100, cex=2, col="blue")

scrit<-round(ccrit/dat$X100_Perc_O2[1],3)*100
scrit<-ifelse(hypoxxrange[1]=="NA",NA,scrit)
scrit_GAM<-scrit
ccrit<-ifelse(hypoxxrange[1]=="NA",NA,ccrit)
new_ccrit<-ccrit*32/1000
new_ccrit_GAM<-new_ccrit



###############################################################
# Method 4: Simple moving average (SMA)
# MMR: use max value from SMA
# RMR: use low value from SMA
# suffix is SMA
###############################################################

datSMA<-datAS

#SMA taking the mean of 3 points (point before, actual point, point after) for every point
roll<-rollmean(datSMA$finalMO2, 3, fill = c(mean(c(datSMA$finalMO2[1],datSMA$finalMO2[2])), 
                                            align = "center", mean(c(datSMA$finalMO2[nrow(datSMA)],datSMA$finalMO2[nrow(datSMA)-1]))))

datSMA$RollMean<-roll

par(mar=c(5,4,4,4))
plot(finalMO2~Date_Timen, main=paste(dat$Animal_ID[1], dat$Test_Temp[1],"SMA", sep="_"),data=dat, xaxt="n", type="b", ylab="Metabolic rate (mg/kg/h)", xlab="")
r<-as.POSIXct(round(range(dat$Date_Time), "hours"))
axis.POSIXct(1, at = seq(r[1], r[2], by = 7200), format = "%H:%M",las=2)
mtext(side=1,line=3.75,"Time of day")
polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10), #for hypoxia trial
        x=c(dat$Date_Timen[min(which(dat$hypoxia=="yes"))],dat$Date_Timen[min(which(dat$hypoxia=="yes"))],
            dat$Date_Timen[max(which(dat$hypoxia=="yes"))],dat$Date_Timen[max(which(dat$hypoxia=="yes"))]),col=rgb(1,0,0,.5))
polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10),  #for night shading
        x=c(dat$Sunset1[min(which(dat$Day_Night=="Night"))],dat$Sunset1[min(which(dat$Day_Night=="Night"))],
            dat$Sunrise1[max(which(dat$Day_Night=="Night"))],dat$Sunrise1[max(which(dat$Day_Night=="Night"))]),col=rgb(.5,.5,.5,.5))
polygon(y=c(min(dat$finalMO2)/10,max(dat$finalMO2)*10,max(dat$finalMO2)*10,min(dat$finalMO2)/10),  #for night shading
        x=c(dat$Sunset1[min(which(dat$Day_Night=="Night2"))],dat$Sunset1[min(which(dat$Day_Night=="Night2"))],
            dat$Date_Timen[max(which(dat$Day_Night=="Night2"))],dat$Date_Timen[max(which(dat$Day_Night=="Night2"))]),col=rgb(.5,.5,.5,.5))
lines(roll~Date_Timen, data=datSMA, col="green", lwd=2)
par(new=T)
plot(Percent_Max_O2*100~Date_Timen, data=dat, type="l", xaxt="n", axes=F, xlab="",ylab="",lwd=3,col="red")
axis(side=4)
mtext(side=4,line=2.5,"Oxygen saturation (%)")

max_time<-datSMA$Date_Timen[which.max(datSMA$RollMean)]
chase_max<-max(chase_pts$Date_Timen)
chase_min<-min(chase_pts$Date_Timen)
MMR_SMA<-ifelse(max_time<=chase_max & max_time>=chase_min,max(datSMA$RollMean,na.rm=T),NA)
MMR_SMA
SMR_SMA<-min(datSMA$RollMean, na.rm=T)
SMR_SMA
AS_SMA<-MMR_SMA-SMR_SMA
AS_SMA

##########
#Scrit
#finds first value where MO2 drops below SMR and the remaining points are below SMR
SMR<-SMR_SMA
brec<-NULL
for(i in 1:nrow(hypox)){
  b<-ifelse(any(hypox$finalMO2[i:nrow(hypox)]>SMR)=="TRUE","bad","good") 
  brec[i]<-b
}
min<-hypox$min..O2.[min(which(brec=="good"))]

#make regression line through points where MO2 drops below SMR and the remaining points are below SMR
regpts<-hypox[which(hypox$min..O2.<=min),]
if(is.na(min)==FALSE){ #if statement for when animal doesn't make it to hypoxia trial
  reglm<-lm(finalMO2~min..O2., data=regpts)
  slope<-coef(reglm)[2]
  intercept<-coef(reglm)[1]
  xrange<-seq(from = min(dat$min..O2.), to= max(dat$min..O2.),by = 5)
  ypts<-(slope*xrange)+intercept
  lines(ypts~xrange)
  ccrit<-as.numeric((SMR-intercept)/slope) #solve x for when y = SMR
}else (ccrit<-NA)

if(is.na(ccrit)==TRUE){
  hypoxxrange<-NA
}else if(min(hypox$min..O2.)<ccrit){
  hypoxxrange<-seq(from = min(hypox$min..O2.), to= ccrit,by = 0.1)
}else(hypoxxrange<-NA)

hypoxypts<-(slope*hypoxxrange)+intercept

#if statement because there are times when the ccrit is NA
if(is.na(ccrit)==TRUE){
  SMRxrange<-NA
}else if(max(dat$min..O2.)>ccrit){
  SMRxrange<-seq(from = ccrit, to= max(dat$min..O2.),by = 1)
}else(SMRxrange<-NA)

SMRypts<-rep(SMR,length(SMRxrange))

plot(finalMO2~min..O2., data=dat, pch=16, xlab="Min O2 (uM)", ylab="Metabolic rate (mg/kg/h)",main=paste(dat$Animal_ID[1], dat$Test_Temp[1],"SMA" , sep="_"))
lines(hypoxypts~hypoxxrange, col="green")
lines(SMRypts~SMRxrange, col="green")
text(x=(.7*max(dat$min..O2.)), y=(0.93*max(dat$finalMO2)),round(ccrit/dat$X100_Perc_O2[1],3)*100, cex=2, col="green")

scrit<-round(ccrit/dat$X100_Perc_O2[1],3)*100
scrit<-ifelse(hypoxxrange[1]=="NA",NA,scrit)
scrit_SMA<-scrit
ccrit<-ifelse(hypoxxrange[1]=="NA",NA,ccrit)
new_ccrit<-ccrit*32/1000
new_ccrit_SMA<-new_ccrit


###############################################################
# Method 5: Frequency distribution
# MMR: median of second hump
# RMR: median of first hump
# suffix is freq
###############################################################
library(diptest)
hist(datAS$finalMO2, breaks=15)
plot(density(datAS$finalMO2))
dip.p<-dip.test(datAS$finalMO2)$"p.value"
dip.p

find_modes<- function(x) {
  modes <- NULL
  for ( i in 2:(length(x)-1) ){
    if ( (x[i] > x[i-1]) & (x[i] > x[i+1]) ) {
      modes <- c(modes,i)
    }
  }
  if ( length(modes) == 0 ) {
    modes = 'This is a monotonic distribution'
  }
  return(modes)
}

mymodes_indices <- find_modes(density(datAS$finalMO2)$y)
density(datAS$finalMO2)$y[mymodes_indices] 
peaks<-density(datAS$finalMO2)$x[mymodes_indices]

SMR_freq<-ifelse(dip.p<0.05,peaks[1],NA)
SMR_freq
MMR_freq<-ifelse(dip.p<0.05,peaks[2],NA)
MMR_freq

AS_freq<-MMR_freq-SMR_freq
AS_freq

##########
#Scrit
#finds first value where MO2 drops below SMR and the remaining points are below SMR
SMR<-SMR_freq
brec<-NULL
for(i in 1:nrow(hypox)){
  b<-ifelse(any(hypox$finalMO2[i:nrow(hypox)]>SMR)=="TRUE","bad","good") 
  brec[i]<-b
}
min<-hypox$min..O2.[min(which(brec=="good"))]

#make regression line through points where MO2 drops below SMR and the remaining points are below SMR
regpts<-hypox[which(hypox$min..O2.<=min),]
if(is.na(min)==FALSE){ #if statement for when animal doesn't make it to hypoxia trial
  reglm<-lm(finalMO2~min..O2., data=regpts)
  slope<-coef(reglm)[2]
  intercept<-coef(reglm)[1]
  xrange<-seq(from = min(dat$min..O2.), to= max(dat$min..O2.),by = 5)
  ypts<-(slope*xrange)+intercept
  lines(ypts~xrange)
  ccrit<-as.numeric((SMR-intercept)/slope) #solve x for when y = SMR
}else (ccrit<-NA)

if(is.na(ccrit)==TRUE){
  hypoxxrange<-NA
}else if(min(hypox$min..O2.)<ccrit){
  hypoxxrange<-seq(from = min(hypox$min..O2.), to= ccrit,by = 0.1)
}else(hypoxxrange<-NA)

hypoxypts<-(slope*hypoxxrange)+intercept

#if statement because there are times when the ccrit is NA
if(is.na(ccrit)==TRUE){
  SMRxrange<-NA
}else if(max(dat$min..O2.)>ccrit){
  SMRxrange<-seq(from = ccrit, to= max(dat$min..O2.),by = 1)
}else(SMRxrange<-NA)

SMRypts<-rep(SMR,length(SMRxrange))

plot(finalMO2~min..O2., data=dat, pch=16, xlab="Min O2 (uM)", ylab="Metabolic rate (mg/kg/h)",main=paste(dat$Animal_ID[1], dat$Test_Temp[1],"SMA" , sep="_"))
lines(hypoxypts~hypoxxrange, col="green")
lines(SMRypts~SMRxrange, col="green")
text(x=(.7*max(dat$min..O2.)), y=(0.93*max(dat$finalMO2)),round(ccrit/dat$X100_Perc_O2[1],3)*100, cex=2, col="green")

scrit<-round(ccrit/dat$X100_Perc_O2[1],3)*100
scrit<-ifelse(hypoxxrange[1]=="NA",NA,scrit)
scrit_freq<-scrit
ccrit<-ifelse(hypoxxrange[1]=="NA",NA,ccrit)
new_ccrit<-ccrit*32/1000
new_ccrit_freq<-new_ccrit




plot(sumdat$MMR3_10[which(sumdat$Animal_ID=="SB02")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB02")],lty=2, type="b",ylim=c(150,300), ylab="Metabolic Rate", xlab="Temperature", main="SB02") 
lines(sumdat$MMR3_3[which(sumdat$Animal_ID=="SB02")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB02")], lty=2,type="b",col="blue")
lines(sumdat$MMR_GAM[which(sumdat$Animal_ID=="SB02")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB02")],lty=2,type="b", col="red")
lines(sumdat$MMR_SMA[which(sumdat$Animal_ID=="SB02")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB02")],lty=2,type="b", col="green")
lines(sumdat$RMR3_10[which(sumdat$Animal_ID=="SB02")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB02")],type="b")
lines(sumdat$RMR3_3[which(sumdat$Animal_ID=="SB02")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB02")],type="b",col="blue")
lines(sumdat$RMR_GAM[which(sumdat$Animal_ID=="SB02")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB02")], type="b",col="red")
lines(sumdat$RMR_SMA[which(sumdat$Animal_ID=="SB02")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB02")], type="b",col="green")

plot(sumdat$MMR3_10[which(sumdat$Animal_ID=="SB03")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB03")], lty=2,type="b",ylim=c(100,400), ylab="Metabolic Rate", xlab="Temperature", main="SB03") 
lines(sumdat$MMR3_3[which(sumdat$Animal_ID=="SB03")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB03")],lty=2,type="b", col="blue")
lines(sumdat$MMR_GAM[which(sumdat$Animal_ID=="SB03")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB03")],lty=2,type="b", col="red")
lines(sumdat$RMR3_10[which(sumdat$Animal_ID=="SB03")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB03")],type="b")
lines(sumdat$RMR3_3[which(sumdat$Animal_ID=="SB03")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB03")],type="b",col="blue")
lines(sumdat$RMR_GAM[which(sumdat$Animal_ID=="SB03")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB03")],type="b", col="red")

plot(sumdat$MMR3_10[which(sumdat$Animal_ID=="SB05")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB05")], lty=2,type="b",ylim=c(75,300), ylab="Metabolic Rate", xlab="Temperature", main="SB05") 
lines(sumdat$MMR3_3[which(sumdat$Animal_ID=="SB05")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB05")],lty=2,type="b", col="blue")
lines(sumdat$MMR_GAM[which(sumdat$Animal_ID=="SB05")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB05")],lty=2,type="b", col="red")
lines(sumdat$MMR_SMA[which(sumdat$Animal_ID=="SB05")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB05")],lty=2,type="b", col="green")
lines(sumdat$RMR3_10[which(sumdat$Animal_ID=="SB05")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB05")],type="b")
lines(sumdat$RMR3_3[which(sumdat$Animal_ID=="SB05")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB05")],type="b",col="blue")
lines(sumdat$RMR_GAM[which(sumdat$Animal_ID=="SB05")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB05")],type="b", col="red")
lines(sumdat$RMR_SMA[which(sumdat$Animal_ID=="SB05")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB05")],type="b", col="green")

plot(sumdat$MMR3_10[which(sumdat$Animal_ID=="SB04")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB04")], lty=2,type="b",ylim=c(150,300), ylab="Metabolic Rate", xlab="Temperature", main="SB04") 
lines(sumdat$MMR3_3[which(sumdat$Animal_ID=="SB04")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB04")],lty=2,type="b", col="blue")
lines(sumdat$MMR_GAM[which(sumdat$Animal_ID=="SB04")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB04")],lty=2,type="b", col="red")
lines(sumdat$MMR_SMA[which(sumdat$Animal_ID=="SB04")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB04")],lty=2,type="b", col="green")
lines(sumdat$RMR3_10[which(sumdat$Animal_ID=="SB04")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB04")],type="b")
lines(sumdat$RMR3_3[which(sumdat$Animal_ID=="SB04")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB04")],type="b",col="blue")
lines(sumdat$RMR_GAM[which(sumdat$Animal_ID=="SB04")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB04")],type="b", col="red")
lines(sumdat$RMR_SMA[which(sumdat$Animal_ID=="SB04")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB04")],type="b", col="green")

plot(sumdat$MMR3_10[which(sumdat$Animal_ID=="SB12")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB12")], lty=2,type="b",ylim=c(150,350), ylab="Metabolic Rate", xlab="Temperature", main="SB12") 
lines(sumdat$MMR3_3[which(sumdat$Animal_ID=="SB12")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB12")],lty=2,type="b", col="blue")
lines(sumdat$MMR_GAM[which(sumdat$Animal_ID=="SB12")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB12")],lty=2,type="b", col="red")
lines(sumdat$MMR_SMA[which(sumdat$Animal_ID=="SB12")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB12")],lty=2,type="b", col="green")
lines(sumdat$RMR3_10[which(sumdat$Animal_ID=="SB12")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB12")],type="b")
lines(sumdat$RMR3_3[which(sumdat$Animal_ID=="SB12")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB12")],type="b",col="blue")
lines(sumdat$RMR_GAM[which(sumdat$Animal_ID=="SB12")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB12")],type="b", col="red")
lines(sumdat$RMR_SMA[which(sumdat$Animal_ID=="SB12")]~sumdat$Test_Temp[which(sumdat$Animal_ID=="SB12")],type="b", col="green")

temps<-unique(sumdat$Test_Temp)
meansAS<-data.frame()
x<-0
for(i in temps){
  x=x+1
  temp<-which(sumdat$Test_Temp==i)
  temp1<-sumdat[temp,]
  avgMMR3_10<-mean(temp1$MMR3_10, na.rm=T)
  sdMMR3_10<-sd(temp1$MMR3_10, na.rm=T)
  seMMR3_10<-sdMMR3_10/sqrt(nrow(temp1))
  avgRMR3_10<-mean(temp1$RMR3_10,na.rm=T)
  sdRMR3_10<-sd(temp1$RMR3_10, na.rm=T)
  seRMR3_10<-sdRMR3_10/sqrt(nrow(temp1))
  avgAS3_10<-mean(temp1$AS3_10,na.rm=T)
  sdAS3_10<-sd(temp1$AS3_10, na.rm=T)
  seAS3_10<-sdAS3_10/sqrt(nrow(temp1))
  avgMMR3_3<-mean(temp1$MMR3_3,na.rm=T)
  sdMMR3_3<-sd(temp1$MMR3_3, na.rm=T)
  seMMR3_3<-sdMMR3_3/sqrt(nrow(temp1))
  avgRMR3_3<-mean(temp1$RMR3_3,na.rm=T)
  sdRMR3_3<-sd(temp1$RMR3_3, na.rm=T)
  seRMR3_3<-sdRMR3_3/sqrt(nrow(temp1))
  avgAS3_3<-mean(temp1$AS3_3,na.rm=T)
  sdAS3_3<-sd(temp1$AS3_3, na.rm=T)
  seAS3_3<-sdAS3_3/sqrt(nrow(temp1))
  avgMMR_GAM<-mean(temp1$MMR_GAM,na.rm=T)
  sdMMR_GAM<-sd(temp1$MMR_GAM, na.rm=T)
  seMMR_GAM<-sdMMR_GAM/sqrt(nrow(temp1))
  avgRMR_GAM<-mean(temp1$RMR_GAM,na.rm=T)
  sdRMR_GAM<-sd(temp1$RMR_GAM, na.rm=T)
  seRMR_GAM<-sdRMR_GAM/sqrt(nrow(temp1))
  avgAS_GAM<-mean(temp1$AS_GAM,na.rm=T)
  sdAS_GAM<-sd(temp1$AS_GAM, na.rm=T)
  seAS_GAM<-sdAS_GAM/sqrt(nrow(temp1))
  avgMMR_SMA<-mean(temp1$MMR_SMA,na.rm=T)
  sdMMR_SMA<-sd(temp1$MMR_SMA, na.rm=T)
  seMMR_SMA<-sdMMR_SMA/sqrt(nrow(temp1))
  avgRMR_SMA<-mean(temp1$RMR_SMA,na.rm=T)
  sdRMR_SMA<-sd(temp1$RMR_SMA, na.rm=T)
  seRMR_SMA<-sdRMR_SMA/sqrt(nrow(temp1))
  avgAS_SMA<-mean(temp1$AS_SMA,na.rm=T)
  sdAS_SMA<-sd(temp1$AS_SMA, na.rm=T)
  seAS_SMA<-sdAS_SMA/sqrt(nrow(temp1))
  meansAS[x,1]<-i
  meansAS[x,2]<-avgMMR3_10
  meansAS[x,3]<-seMMR3_10
  meansAS[x,4]<-avgRMR3_10
  meansAS[x,5]<-seRMR3_10
  meansAS[x,6]<-avgAS3_10
  meansAS[x,7]<-seAS3_10
  meansAS[x,8]<-avgMMR3_3
  meansAS[x,9]<-seMMR3_3
  meansAS[x,10]<-avgRMR3_3
  meansAS[x,11]<-seRMR3_3
  meansAS[x,12]<-avgAS3_3
  meansAS[x,13]<-seAS3_3
  meansAS[x,14]<-avgMMR_GAM
  meansAS[x,15]<-seMMR_GAM
  meansAS[x,16]<-avgRMR_GAM
  meansAS[x,17]<-seRMR_GAM
  meansAS[x,18]<-avgAS_GAM
  meansAS[x,19]<-seAS_GAM
  meansAS[x,20]<-avgMMR_SMA
  meansAS[x,21]<-seMMR_SMA
  meansAS[x,22]<-avgRMR_SMA
  meansAS[x,23]<-seRMR_SMA
  meansAS[x,24]<-avgAS_SMA
  meansAS[x,25]<-seAS_SMA
}

colnames(meansAS)<-c("Temp","MMR3_10","SE_MMR3_10","RMR3_10","SE_RMR3_10","AS3_10","SE_AS3_10","MMR3_3","SE_MMR3_3",
                     "RMR3_3","SE_RMR3_3","AS3_3","SE_AS3_3","MMR_GAM","SE_MMR_GAM","RMR_GAM","SE_RMR_GAM","AS_GAM","SE_AS_GAM",
                     "MMR_SMA","SE_MMR_SMA","RMR_SMA","SE_RMR_SMA","AS_SMA","SE_AS_SMA")



#attempted hockey stick method for determining Ccrit, doesn't estimate well
#hockey stick Ccrit function
hs=function(p){
  slope=exp(p[1])
  Ccrit=exp(p[2])
  intercept<-exp(p[3])
  
  predCcrit=ifelse(dat$min..O2.<=Ccrit, slope*dat$min..O2.+intercept,SMR)
  
  residHS=(log(dat$finalMO2)-log(predCcrit))^2
  SS=sum(residHS)
  return(SS)
}

p=log(c(1.5,110,50)) #slope, Ccrit, intercept starting values
predCcrit=numeric(length(dat$finalMO2))

outHS=optim(p,hs)

RSS=outHS$value
n=nrow(dat)
AIC.HS=n*(1+log(2*pi*(RSS/n)))+2*(length(p)+1)

#extract out predicted values
slope=exp(outHS$par[1])
Ccrit=exp(outHS$par[2])
intercept=exp(outHS$par[3])

x1=subset(dat$min..O2.,dat$min..O2.<Ccrit)
x2=subset(dat$min..O2.,dat$min..O2.>=Ccrit)

y1=x1*slope+intercept
y2=rep(SMR,length(x2))
xnew=c(x1,x2)
ynew=c(y1,y2)
hockey=cbind(ynew,xnew)
lines(xnew,ynew,col="black",lwd=2)






temps<-c(24,28,32)
plot(MMR3_10~Temp, data=meansAS,type="b",ylab="Metabolic rate (mg/kg/h)",xlab="Temperature (°C)",ylim=c(125,350))
lines(RMR3_10~Temp,data=meansAS, type="b")
arrows(y0=(meansAS$MMR3_10-meansAS$SE_MMR3_10),y1=(meansAS$MMR3_10+meansAS$SE_MMR3_10), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1)
arrows(y0=(meansAS$RMR3_10-meansAS$SE_RMR3_10),y1=(meansAS$RMR3_10+meansAS$SE_RMR3_10), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1)
lines(MMR3_3~Temp,data=meansAS, type="b", col="blue")
lines(RMR3_3~Temp,data=meansAS, type="b", col="blue")
arrows(y0=(meansAS$MMR3_3-meansAS$SE_MMR3_3),y1=(meansAS$MMR3_3+meansAS$SE_MMR3_3), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="blue")
arrows(y0=(meansAS$RMR3_3-meansAS$SE_RMR3_3),y1=(meansAS$RMR3_3+meansAS$SE_RMR3_3), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="blue")
lines(RMR_GAM~Temp,data=meansAS, type="b", col="red")
lines(MMR_GAM~Temp,data=meansAS, type="b", col="red")
arrows(y0=(meansAS$MMR_GAM-meansAS$SE_MMR_GAM),y1=(meansAS$MMR_GAM+meansAS$SE_MMR_GAM), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="red")
arrows(y0=(meansAS$RMR_GAM-meansAS$SE_RMR_GAM),y1=(meansAS$RMR_GAM+meansAS$SE_RMR_GAM), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="red")
lines(RMR_SMA~Temp,data=meansAS, type="b", col="green")
lines(MMR_SMA~Temp,data=meansAS, type="b", col="green")
arrows(y0=(meansAS$MMR_SMA-meansAS$SE_MMR_SMA),y1=(meansAS$MMR_SMA+meansAS$SE_MMR_SMA), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="green")
arrows(y0=(meansAS$RMR_SMA-meansAS$SE_RMR_SMA),y1=(meansAS$RMR_SMA+meansAS$SE_RMR_SMA), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="green")


temps<-c(24,28,32)
plot(AS3_10~Temp, data=meansAS,type="b",ylab="Aerobic Scope (mg/kg/h)",xlab="Temperature (°C)",ylim=c(25,175))
arrows(y0=(meansAS$AS3_10-meansAS$SE_AS3_10),y1=(meansAS$AS3_10+meansAS$SE_AS3_10), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1)
lines(AS3_3~Temp,data=meansAS, type="b", col="blue")
arrows(y0=(meansAS$AS3_3-meansAS$SE_AS3_3),y1=(meansAS$AS3_3+meansAS$SE_AS3_3), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="blue")
lines(AS_GAM~Temp,data=meansAS, type="b", col="red")
arrows(y0=(meansAS$AS_GAM-meansAS$SE_AS_GAM),y1=(meansAS$AS_GAM+meansAS$SE_AS_GAM), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="red")
lines(AS_SMA~Temp,data=meansAS, type="b", col="green")
arrows(y0=(meansAS$AS_SMA-meansAS$SE_AS_SMA),y1=(meansAS$AS_SMA+meansAS$SE_AS_SMA), 
       x0=as.numeric(temps),x1=as.numeric(temps),code=3, angle=90,length=.1,col="green")




par(mfrow=c(2,2))
boxplot(AS3_10~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Aerobic scope", main="3_10 Aerobic Scope")
boxplot(AS3_3~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Aerobic scope", main="3_3 Aerobic Scope")
boxplot(AS_GAM~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Aerobic scope", main="GAM Aerobic Scope")
boxplot(AS_SMA~Test_Temp, data=sumdat, xlab="Temperature (°C)", ylab="Aerobic scope", main="SMA Aerobic Scope")





