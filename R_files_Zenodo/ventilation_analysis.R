setwd("~/Documents/PhD Project/Respirometry/Respirometry")
sumdat<-read.csv("MO2_methods10_cobia_20180914.csv")

setwd("~/Documents/PhD Project/Ventilation")

vent<-read.csv("cobia_ventilation.csv")
vent$Date_Time<-paste(vent$Date,vent$time, sep=" ")
vent$Date_Time<-as.POSIXct(vent$Date_Time, format="%m/%d/%y %H:%M:%S")

names(vent)[5]<-"Ventilation_Rate"

library(StreamMetabolism)

#read in all csvs ending with "data_corrected.csv" and putting csvs in list
mycsv = dir(pattern="data_corrected.csv")
n <- length(mycsv)
mylist <- vector("list", n)
for(i in 1:n) mylist[[i]] <- read.csv(mycsv[i])


datRec<-NULL
for(q in 1:length(mylist)){
  dat<-mylist[q]
  dat<-data.frame(dat)
  dat$Date_Time<-paste(dat$Date, dat$time, seq=" ")
  dat$Date_Time<-as.POSIXct(dat$Date_Time, format="%m/%d/%y %H:%M:%S")
  #dat<-dat[which(dat$r2>=0.80),] #gets rid of points where r2 is lower than .80
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
  
  vent_dat<-vent[which(as.character(vent$Animal_ID)==as.character(dat$Animal_ID[1]) & vent$Test_Temp==dat$Test_Temp[1]),]
  
  merged<-merge(dat,vent_dat,by="Date_Time",all.x=T)
  
  dat<-cbind(dat,Time.of.Gape=merged$Time.of.Gape,Ventilation_Rate=merged$Ventilation_Rate,Mouth_Gape=merged$Mouth.Gape)
  
  #want to divide mouth gape by minimum mouth gape during trial so that when mouth gape and vent rate are multiplied
  #the VRMGM is equal to just the ventilation rate at the smallest mouth gape, then as mouth gape increases so will its impact on ventilation rate and this VRMGM
  dat$Mouth_Gape_to_min<-dat$Mouth_Gape/min(dat$Mouth_Gape,na.rm=T)
  
  dat$VentMouthGape_Metric<-dat$Ventilation_Rate*dat$Mouth_Gape_to_min
  
  #plot entire trial, remember not all trials have two nights so will see warnings for some trials
  setwd("~/Documents/PhD Project/Ventilation/Trials")
  #png(file=paste(dat$Animal_ID[1], dat$Test_Temp[1],"ventilation","png", sep="."),width=650,height=450)
  tiff(file=paste(dat$Animal_ID[1], dat$Test_Temp[1],"ventilation","tiff", sep="."), width=30, height=15, units="cm",compression="lzw",res=150)
  par(mfrow=c(1,1))
  par(mar=c(6,7.3,2,10))
  plot(finalMO2~Date_Timen, main=paste(dat$Animal_ID[1], dat$Test_Temp[1],sep="_") ,data=dat, xaxt="n", type="b",lwd=2, ylab="", xlab="")
  r<-as.POSIXct(round(range(dat$Date_Time), "hours"))
  axis.POSIXct(1, at = seq(r[1], r[2], by = 7200), format = "%H:%M",las=2)
  mtext(side=1,line=4,"Time of day",cex=1.25)
  mtext(side=2,line=2,expression(paste("Metabolic Rate (mg kg"^"-1","h"^"-1",")")),cex=1.25)
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
  par(new=T)
  plot(Ventilation_Rate~Date_Timen, data=dat, type="l", xaxt="n", axes=F, xlab="",ylab="",lwd=1,lty=2,col="blue")
  axis(side=4, line=3.75)
  mtext(side=4,line=5.75,"Ventilation Rate")
  par(new=T)
  plot(VentMouthGape_Metric~Date_Timen, data=dat, type="b", xaxt="n", axes=F, xlab="",ylab="",lwd=1,lty=2,col="deeppink2")
  axis(side=2, line=3.75)
  mtext(side=2,line=6,"VRMG") #for Vent Rate Mouth Gape Metric
  par(new=T)
  plot(Mouth_Gape~Date_Timen, data=dat, type="l", xaxt="n", axes=F, xlab="",ylab="",lwd=1,lty=2,col="dark green")
  axis(side=4, line=7)
  mtext(side=4,line=9,"Mouth Gape (cm)")
  
  dev.off()
  
  SMR<-sumdat$SMR[which(as.character(sumdat$Animal_ID)==as.character(dat$Animal_ID[1]) & as.character(sumdat$Test_Temp)==as.character(dat$Test_Temp[1]))]
  SMRfirst<-min(which(dat$finalMO2<SMR))
  dat$before_after_SMR<-NA
  dat$before_after_SMR[1:(SMRfirst-1)]<-"before"
  dat$before_after_SMR[SMRfirst:nrow(dat)]<-"after"
  
  dat<-dat[which(dat$r2>=0.80),]
  
  datRec<-rbind(datRec,dat)
  print(q)
}

#add TL to dataframe
datRec$TL<-NA
datRec$TL[which(datRec$Animal_ID=="C08")]<-95
datRec$TL[which(datRec$Animal_ID=="C09")]<-88
datRec$TL[which(datRec$Animal_ID=="C11")]<-86
datRec$TL[which(datRec$Animal_ID=="C12")]<-105
datRec$TL[which(datRec$Animal_ID=="C13")]<-108
datRec$TL[which(datRec$Animal_ID=="C15")]<-84
datRec$TL[which(datRec$Animal_ID=="C16")]<-87
datRec$TL[which(datRec$Animal_ID=="C17")]<-88
datRec$TL[which(datRec$Animal_ID=="C18")]<-92
datRec$TL[which(datRec$Animal_ID=="C19")]<-91.5
datRec$TL[which(datRec$Animal_ID=="C20")]<-87.5
datRec$TL[which(datRec$Animal_ID=="C21")]<-82
datRec$TL[which(datRec$Animal_ID=="C22")]<-96
datRec$TL[which(datRec$Animal_ID=="C23")]<-83
datRec$TL[which(datRec$Animal_ID=="C24")]<-81
datRec$TL[which(datRec$Animal_ID=="C25")]<-88
datRec$TL[which(datRec$Animal_ID=="C26")]<-108
datRec$TL[which(datRec$Animal_ID=="C27")]<-88.5
datRec$TL[which(datRec$Animal_ID=="C28")]<-102
datRec$TL[which(datRec$Animal_ID=="C30")]<-81
datRec$TL[which(datRec$Animal_ID=="C31")]<-85
datRec$TL[which(datRec$Animal_ID=="C32")]<-81
datRec$TL[which(datRec$Animal_ID=="C33")]<-96
datRec$TL[which(datRec$Animal_ID=="C34")]<-89
datRec$TL[which(datRec$Animal_ID=="C35")]<-102.5

#remember not including inds that died (C31_28, C24_32, C32_32, C28_32)
datRec$Alive_Dead<-ifelse(datRec$Animal_ID=="C31" & datRec$Test_Temp=="28", "dead",
                         ifelse(datRec$Animal_ID=="C24" & datRec$Test_Temp=="32", "dead",
                                ifelse(datRec$Animal_ID=="C32" & datRec$Test_Temp=="32", "dead",
                                       ifelse(datRec$Animal_ID=="C28" & datRec$Test_Temp=="32","dead","alive"))))

datRec_dead<-datRec[which(datRec$Alive_Dead=="dead"),]
datRec<-datRec[which(datRec$Alive_Dead=="alive"),]

#####add Crit and Scrit to cobia dataframe
#need to pull it from cobia_parts_sas
setwd("~/Documents/PhD Project/Respirometry/Respirometry")
ccritdat<-read.csv("MO2_methods10_cobia_20180914.csv")

#assign NAs to Ccrit and Scrit for 30_28 and 34_28 because the automation method to calculate ccrit from SMR was not accurate for those two fish
ccritdat$Ccrit[which(ccritdat$Animal_ID=="C30"&ccritdat$Test_Temp=="28")]<-NA
ccritdat$Scrit[which(ccritdat$Animal_ID=="C30"&ccritdat$Test_Temp=="28")]<-NA
ccritdat$Ccrit[which(ccritdat$Animal_ID=="C34"&ccritdat$Test_Temp=="28")]<-NA
ccritdat$Scrit[which(ccritdat$Animal_ID=="C34"&ccritdat$Test_Temp=="28")]<-NA

ccritdat<-subset(ccritdat, select=c(Animal_ID, Test_Temp, Ccrit))

ccritdat$Ccritnew<-ccritdat$Ccrit*1000/32#puts ccrit into same units as O2 min units (umol/l)

#Test_Temp should be factor
datRec$Test_Temp<-as.factor(datRec$Test_Temp)

#going to use O2 middway between min and max
datRec$midO2<-(datRec$max..O2.+datRec$min..O2.)/2

#add observation column
animal<-unique(datRec$Animal_ID)
temp<-unique(datRec$Test_Temp)
obs_list<-NULL
test_list<-NULL
crash<-NULL
for(i in animal){
  a<-which(datRec$Animal_ID==i)
  a1<-datRec[a,]
  s<-which(ccritdat$Animal_ID==i)
  s1<-ccritdat[s,]
  for(j in temp){
    t<-which(a1$Test_Temp==j)
    t1<-a1[t,]
    p<-which(s1$Test_Temp==j)
    p1<-s1[p,]
    obs<-seq(from=1, to=nrow(t1))
    obs_list<-c(obs_list,obs)
    ccrit<-p1$Ccritnew[1]
    b_a<-ifelse(t1$min..O2.>p1$Ccritnew,"before","after")
    crash<-c(crash, b_a)
  }
  test<-seq(from=1, to=nrow(a1))
  test_list<-c(test_list,test)
}

badlocs<-which(obs_list==0)
badlocs1<-badlocs+1
locs<-c(badlocs,badlocs1)
obs_list<-obs_list[-locs]

datRec$Obs<-obs_list
datRec$test<-test_list
datRec$Crash<-crash


datRec<-subset(datRec, select=c(Date, time, animal.mass,midO2,Animal_ID,Test_Temp,X100_Perc_O2,
                                Date_Time,Time_Diff,finalMO2,Day_Night, hypoxia,Ventilation_Rate,Mouth_Gape,
                                VentMouthGape_Metric, TL, Obs, test, Crash,before_after_SMR))

#use Time_Diff as proxy for time after chase
datRec$Time_Diff<-datRec$Time_Diff+60 #added 60 seconds because it is roughly 60s from chase to start of trial, will use time_diff as proxy for time after chase
datRec$Crash<-as.factor(datRec$Crash)

#need to get rid of NA rows
cobia<-datRec[complete.cases(datRec$VentMouthGape_Metric),]#will have a lot of eliminated for when fish wasn't facing camera

#######################
###Normoxia analysis###
######################

cobia_norm<-cobia[which(cobia$hypoxia=="no"),]
cobia_norm$Animal_ID<-as.factor(as.character(cobia_norm$Animal_ID)) #do this to get rid of 2 cobia that died that were only tested at one temp to get rid of that level
cobia_norm$before_after_SMR<-as.factor(cobia_norm$before_after_SMR)

####Check for correlations in covariates
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(cbind(cobia_norm$Test_Temp,cobia_norm$Ventilation_Rate,cobia_norm$Mouth_Gape,cobia_norm$VentMouthGape_Metric,
            cobia_norm$Time_Diff,cobia_norm$animal.mass,cobia_norm$TL,cobia_norm$before_after_SMR),lower.panel = panel.cor)
#animal.mass and TL are correlated so I took animal.mass out
#time_diff and before_after_SMR are correlated, use time_diff
#everything else seems okay!


#######modeling
library(nlme)
library(HH)
mod0<-lm(finalMO2~Ventilation_Rate+Mouth_Gape+VentMouthGape_Metric+Test_Temp+TL+Time_Diff, data=cobia_norm,na.action = na.omit,x=TRUE)
vif(mod0) #covariates are not collinear, below 5 is okay

#no correlation structure and not modeling variance structure (heterogenity)
mod1<-lme(finalMO2~Ventilation_Rate+Mouth_Gape+VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit, method = "ML")
summary(mod1)
#diagnostics
plot(mod1)
plot(resid(mod1)~cobia_norm$Test_Temp)
plot(resid(mod1)~cobia_norm$TL)
plot(resid(mod1)~cobia_norm$VentMouthGape_Metric)
#equal variance looks slightly skewed in VentMouthGape_Metric, but can't seem to fix it but messing with weights, will ignore for now
qqnorm(mod1) #eh not too good
mod1.2<-lme(log(finalMO2)~Ventilation_Rate+Mouth_Gape+VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit, method = "ML")
plot(mod1.2) #much better now
plot(resid(mod1.2)~cobia_norm$Test_Temp)
plot(resid(mod1.2)~cobia_norm$TL)
plot(resid(mod1.2)~cobia_norm$VentMouthGape_Metric)
library(car)
vif(mod1.2) #looks like VRMG and Mouth Gape are colinear, will get rid of Mouth Gape and Ventilation rate from model
mod1.3<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit, method = "ML")
plot(mod1.3) #much better now
plot(resid(mod1.3)~cobia_norm$Test_Temp)
plot(resid(mod1.3)~cobia_norm$TL)
plot(resid(mod1.3)~cobia_norm$VentMouthGape_Metric)
qqnorm(mod1.3)#better
vif(mod1.3) #better
#have this just so have model not log transformed, doesn't seem as good as log transformed one
mod1.4<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit, method = "ML")

#check for temporal autocorrelation (ie. independence vs non independence)
E<-residuals(mod1.3, type="normalized")
I1<-!is.na(cobia_norm$VentMouthGape_Metric)
Efull<-vector(length=length(cobia_norm$Test_Temp))
Efull<-NA
Efull[I1]<-E
acf(Efull,na.action = na.pass) #definitely independence violated

#no correlation structure, model variance structure (heterogenity)
mod2<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,weights = varIdent(form=~1|Test_Temp), method="ML")
BIC(mod1.3,mod2) #model with variance structure is better

#test different correlation structures
mod3<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod4<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corCompSymm(form=~Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
#mod4 not converging anyway
BIC(mod1.3,mod2,mod3,mod4)#model with AR1 correlation structure is best

#Now adjust fixed covariates
library(mgcv)
mod3<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod5<-lme(log(finalMO2)~VentMouthGape_Metric*TL+Test_Temp+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod6<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod7<-lme(log(finalMO2)~VentMouthGape_Metric*Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod8<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod9<-gam(log(finalMO2)~s(VentMouthGape_Metric)+Test_Temp+s(Time_Diff)+s(Animal_ID,bs="re"), data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),method="ML",family=gaussian(link = "identity"))
mod10<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")

BIC(mod3,mod5,mod6,mod7,mod8,mod9,mod10) #best model is mod3 (BIC is -1504.1309)


#switch to reml
mod3a<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="REML")
summary(mod3a)
#VRMG is sig
#TL is sig
#Time_Diff is sig
anova(mod3a)
#temp is sig

library(multcomp)
library(lsmeans)
#look at interaction between TBF and Temp
lsmeans(mod3a, pairwise~Test_Temp, adjust="tukey")
#all are sig different from each other

#rsquared calc 
#conditional R2: describes the proportion of variance explained by both the fixed and random factors
#marginal R2: describes the proportion of variance explained by the fixed factor(s) alone
#either package will do
library(piecewiseSEM) #John Lefcheck's package
rsquared(mod3a)
library(gabtool)
r.squared(mod3a)
#conditional R2 is 0.71 which is very good

library(doBy)

disp.lme=function(mod=NA){
  E=resid(mod,type='pearson')
  d=sum(E^2)/mod$fixDF$X[1]
  return(d)
}

#try to sort out log space and random effects...couldn't really find a solution so just left data in log space for figures


######
#Bootstrapping to get errors for plots
#these are for the covariates in best model
sumfun <- function(x, ...){ #this function fixes the problem of when NAs are introduced before summaryBy, get rid of NAs before mean is taken, NAs are made when during bootstrap no points from an ind is selected so would get an NA when trying to predict over that ind
  c(mean=mean(x, na.rm=TRUE, ...), var=var(x, na.rm=TRUE, ...), length=length(x))
}
recVRMG<-matrix(nrow=100, ncol=1000, NA)
recTD<-matrix(nrow=100, ncol=1000, NA)
recTL<-matrix(nrow=100, ncol=1000, NA)
recTemp<-matrix(nrow=3, ncol=1000, NA)
set.seed(2)
for(i in 1:1000){
  df=cobia_norm[sample(nrow(cobia_norm),nrow(cobia_norm),replace=TRUE),]
  mod3a<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=df,na.action = na.omit,
             correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="REML")
  #VRMG
  newdatVRMG.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),Time_Diff=mean(df$Time_Diff),TL=mean(df$TL),
                               VentMouthGape_Metric=seq(from=min(df$VentMouthGape_Metric,na.rm = T),to=max(df$VentMouthGape_Metric,na.rm = T),length=100), 
                               Animal_ID=levels(df$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
  pVRMG.norm<-cbind(newdatVRMG.norm,pred=predict(mod3a, newdata = newdatVRMG.norm,se.fit=TRUE, type="response"))
  #pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred+(0.5*disp.lme(mod3a)))
  pVRMG.normV<-summaryBy(pred~VentMouthGape_Metric, data=pVRMG.norm)
  recVRMG[1:100,i]<-pVRMG.normV$pred.mean
  #Time_Diff
  newdatTD.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_Metric=mean(df$VentMouthGape_Metric),TL=mean(df$TL),
                             Time_Diff=seq(from=min(df$Time_Diff,na.rm = T),to=max(df$Time_Diff,na.rm = T),length=100), 
                             Animal_ID=levels(df$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
  pTD.norm<-cbind(newdatTD.norm,pred=predict(mod3a, newdata = newdatTD.norm,se.fit=TRUE, type="response"))
  #pTD.norm$pred.corrected<-exp(pTD.norm$pred+(0.5*disp.lme(mod3a)))
  pTD.normTD<-summaryBy(pred~Time_Diff, data=pTD.norm)
  recTD[1:100,i]<-pTD.normTD$pred.mean
  #TL
  newdatTL.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_Metric=mean(df$VentMouthGape_Metric),Time_Diff=mean(df$Time_Diff),
                             TL=seq(from=min(df$TL,na.rm = T),to=max(df$TL,na.rm = T),length=100), 
                             Animal_ID=levels(as.factor(as.character(df$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
  pTL.norm<-cbind(newdatTL.norm,pred=predict(mod3a, newdata = newdatTL.norm,se.fit=TRUE, type="response"))
  #pTL.norm$pred.corrected<-exp(pTL.norm$pred+(0.5*disp.lme(mod3a)))
  pTL.normTL<-summaryBy(pred~TL, data=pTL.norm)
  recTL[1:100,i]<-pTL.normTL$pred.mean
  #Temp
  newdatTemp.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_Metric=mean(df$VentMouthGape_Metric),Time_Diff=mean(df$Time_Diff),
                               TL=mean(df$TL), 
                               Animal_ID=levels(df$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
  pTemp.norm<-cbind(newdatTemp.norm,pred=predict(mod3a, newdata = newdatTemp.norm,se.fit=TRUE, type="response"))
  #pTemp.norm$pred.corrected<-exp(pTemp.norm$pred+(0.5*disp.lme(mod3a)))
  pTemp.normT<-summaryBy(pred~Test_Temp, data=pTemp.norm)
  recTemp[1:3,i]<-pTemp.normT$pred.mean
  
  print(i)
}

#calculating confidence intervals from bootstrapping
CIrec<-NULL
CIfunc<-function(data=NA){
  
  for(i in 1:nrow(data)){
    CI<-1.96*sd(data[i,],na.rm=T)
    CIrec<-c(CIrec,CI)
  }
  return(CIrec)
}

outVRMG<-CIfunc(recVRMG)
outTD<-CIfunc(recTD)
outTL<-CIfunc(recTL)
outTemp<-CIfunc(recTemp)


####predicting to get actual model output (non bootstrapping)
#rerun model so mod3a can be used to get actual predictions again (model was slightly different above due to boostrapping)
mod3a<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="REML")

#for VRMG
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Time_Diff=mean(cobia_norm$Time_Diff),TL=mean(cobia_norm$TL),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100), 
                             Animal_ID=levels(cobia_norm$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
pVRMG.norm<-cbind(newdatVRMG.norm,pred=predict(mod3a, newdata = newdatVRMG.norm,se.fit=TRUE, type="response"))
#pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred+(0.5*disp.lme(mod3a)))
pVRMG.normV<-summaryBy(pred~VentMouthGape_Metric, data=pVRMG.norm)
pVRMG.normV$UpperCI<-pVRMG.normV$pred.mean+outVRMG #add CI to mean
pVRMG.normV$LowerCI<-pVRMG.normV$pred.mean-outVRMG #subtract CI from mean

newdatTD.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),VentMouthGape_Metric=mean(cobia_norm$VentMouthGape_Metric),TL=mean(cobia_norm$TL),
                             Time_Diff=seq(from=min(cobia_norm$Time_Diff,na.rm = T),to=max(cobia_norm$Time_Diff,na.rm = T),length=100), 
                             Animal_ID=levels(as.factor(as.character(cobia_norm$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTD.norm<-cbind(newdatTD.norm,pred=predict(mod3a, newdata = newdatTD.norm,se.fit=TRUE, type="response"))
#pTD.norm$pred.corrected<-exp(pTD.norm$pred+(0.5*disp.lme(mod3a)))
pTD.normTD<-summaryBy(pred~Time_Diff, data=pTD.norm)
pTD.normTD$UpperCI<-pTD.normTD$pred.mean+outTD #add CI to mean
pTD.normTD$LowerCI<-pTD.normTD$pred.mean-outTD #subtract CI from mean

newdatTL.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),VentMouthGape_Metric=mean(cobia_norm$VentMouthGape_Metric),Time_Diff=mean(cobia_norm$Time_Diff),
                           TL=seq(from=min(cobia_norm$TL,na.rm = T),to=max(cobia_norm$TL,na.rm = T),length=100), 
                           Animal_ID=levels(as.factor(as.character(cobia_norm$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTL.norm<-cbind(newdatTL.norm,pred=predict(mod3a, newdata = newdatTL.norm,se.fit=TRUE, type="response"))
#pTL.norm$pred.corrected<-exp(pTL.norm$pred+(0.5*disp.lme(mod3a)))
pTL.normTL<-summaryBy(pred~TL, data=pTL.norm)
pTL.normTL$UpperCI<-pTL.normTL$pred.mean+outTL #add CI to mean
pTL.normTL$LowerCI<-pTL.normTL$pred.mean-outTL #subtract CI from mean

newdatTemp.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),VentMouthGape_Metric=mean(cobia_norm$VentMouthGape_Metric),Time_Diff=mean(cobia_norm$Time_Diff),
                           TL=mean(cobia_norm$TL), 
                           Animal_ID=levels(as.factor(as.character(cobia_norm$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTemp.norm<-cbind(newdatTemp.norm,pred=predict(mod3a, newdata = newdatTemp.norm,se.fit=TRUE, type="response"))
#pTemp.norm$pred.corrected<-exp(pTemp.norm$pred+(0.5*disp.lme(mod3a)))
pTemp.normT<-summaryBy(pred~Test_Temp, data=pTemp.norm)
pTemp.normT$UpperCI<-pTemp.normT$pred.mean+outTemp #add CI to mean
pTemp.normT$LowerCI<-pTemp.normT$pred.mean-outTemp #subtract CI from mean


#plot covariates
par(mar=(c(5,4.5,4,1)))
plot(log(finalMO2)~VentMouthGape_Metric, data=cobia_norm,pch=21, bg=c("blue","dark grey","red")[cobia_norm$Test_Temp], ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)
legend("topright",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pVRMG.normV$VentMouthGape_Metric,rev(pVRMG.normV$VentMouthGape_Metric)), y=c(pVRMG.normV$LowerCI,rev(pVRMG.normV$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~VentMouthGape_Metric,data=pVRMG.normV, lwd=2, pch=16)

plot(log(finalMO2)~Time_Diff, data=cobia_norm,pch=21, bg=c("blue","dark grey","red")[cobia_norm$Test_Temp], ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="Time since chase (s)",cex.axis=1.25,cex.lab=1.25)
legend("topright",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pTD.normTD$Time_Diff,rev(pTD.normTD$Time_Diff)), y=c(pTD.normTD$LowerCI,rev(pTD.normTD$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~Time_Diff,data=pTD.normTD, lwd=2, pch=16)

plot(log(finalMO2)~TL, data=cobia_norm,pch=21, bg=c("blue","dark grey","red")[cobia_norm$Test_Temp], ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)
legend("topright",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pTL.normTL$TL,rev(pTL.normTL$TL)), y=c(pTL.normTL$LowerCI,rev(pTL.normTL$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~TL,data=pTL.normTL, lwd=2, pch=16)

barplot<-barplot(height=pTemp.normT$pred.mean, names.arg=pTemp.normT$Test_Temp,ylab="log[Metabolic Rate (mg/kg/h)]",xlab="Temperature (°C)", ylim=c(0,5.5))
arrows(y0=pTemp.normT$LowerCI,y1=pTemp.normT$UpperCI, 
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)


cobia_norm$chase_time<-as.factor(ifelse(cobia_norm$Time_Diff<=10800,"yes","no"))


par(mar=(c(5,4.5,4,1)))
plot(log(finalMO2)~VentMouthGape_Metric, data=cobia_norm,pch=21, bg=c("blue","dark grey")[cobia_norm$chase_time], ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)
legend("topright",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pVRMG.normV$VentMouthGape_Metric,rev(pVRMG.normV$VentMouthGape_Metric)), y=c(pVRMG.normV$LowerCI,rev(pVRMG.normV$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~VentMouthGape_Metric,data=pVRMG.normV, lwd=2, pch=16)





#################
#Hypoxia analysis
#################

cobia_hypox<-cobia[which(cobia$hypoxia=="yes"),]#for model

pairs(cbind(cobia_hypox$Test_Temp,cobia_hypox$Ventilation_Rate,cobia_hypox$Mouth_Gape,cobia_hypox$VentMouthGape_Metric,
            cobia_hypox$animal.mass,cobia_hypox$TL,cobia_hypox$midO2,cobia_hypox$Crash),lower.panel = panel.cor)
#get NAs for crash metric because no ccrit for C34_28 or C30_28, so have NAs there and throws NAs, its ok though cause crash metric is correlated with midO2
#animal.mass and TL are correlated so I took animal.mass out
#seems like midO2 and crash are correlated, get rid of Crash metric
#everything else seems okay!


#######modeling
library(nlme)
library(HH)
h0<-lm(finalMO2~Ventilation_Rate+Mouth_Gape+VentMouthGape_Metric+Test_Temp+TL+midO2, data=cobia_hypox,na.action = na.omit,x=TRUE)
vif(h0) #seems good

#no correlation structure and not modeling variance structure (heterogenity)
h1<-lme(finalMO2~Ventilation_Rate+Mouth_Gape+VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,method="ML")
BIC(h0,h1)
vif(h1) #now mouth gape, and VRMG are colinear, will get rid of ventilation rate and mouth gape from model

h1.1<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,method="ML")
vif(h1.1) #all better

summary(h1.1)
#diagnostics
plot(h1.1)
plot(resid(h1.1)~cobia_hypox$Test_Temp)
plot(resid(h1.1)~cobia_hypox$TL)
plot(resid(h1.1)~cobia_hypox$VentMouthGape_Metric)
plot(resid(h1.1)~cobia_hypox$midO2)
#looks pretty homogenous but lets model variance really quick and see if it improves the model

#no correlation structure, model variance structure (heterogenity)
h2<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,weights = varIdent(form=~1|Test_Temp),method="ML")
h2.1<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,weights = varExp(form=~VentMouthGape_Metric),method="ML")
h2.2<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,weights = varFixed(~VentMouthGape_Metric),method="ML")
h2.3<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,weights = varPower(form=~VentMouthGape_Metric),method="ML")
BIC(h1.1,h2,h2.1,h2.2,h2.3) #model with temp in weights is best

#check for temporal autocorrelation (ie. independence vs non independence)
E<-residuals(h2, type="normalized")
I1<-!is.na(cobia_hypox$VentMouthGape_Metric)
Efull<-vector(length=length(cobia_hypox$Test_Temp))
Efull<-NA
Efull[I1]<-E
acf(Efull,na.action = na.pass) #independence violated slightly

#test different correlation structures
h3<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h4<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corCompSymm(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
BIC(h2,h3,h4)#model with AR1 correlation structure is best


#Now adjust fixed covariates
h3<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h5<-lme(finalMO2~VentMouthGape_Metric+Test_Temp+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h6<-lme(finalMO2~VentMouthGape_Metric*TL+Test_Temp+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h7<-lme(finalMO2~VentMouthGape_Metric+Test_Temp*midO2+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h8<-lme(finalMO2~VentMouthGape_Metric*Test_Temp+midO2+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h9<-lme(finalMO2~VentMouthGape_Metric*midO2+Test_Temp+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h10<-gam(finalMO2~s(VentMouthGape_Metric)+Test_Temp+s(midO2)+TL+s(Animal_ID,bs="re"), data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),method="ML",family=gaussian(link = "identity"))

BIC(h3,h5,h6,h7,h8,h9,h10) #best model is h9 (BIC was 2710.874)

#switch to reml
h9a<-lme(finalMO2~VentMouthGape_Metric*midO2+Test_Temp+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="REML")

summary(h9a)
#TL is sig
#interaction between VRMG and midO2 are sig
anova(h9a)
#temp is sig

#look at pairwise comparision
lsmeans(h9a, pairwise~Test_Temp, adjust="tukey")
#all are sig different from each other

#rsquared calc 
#either package will do
library(piecewiseSEM) #John Lefcheck's package
rsquared(h9a)
library(gabtool)
r.squared(h9a)
#conditional R2 is 0.80 which is very good


######
#Bootstrapping to get errors for plots
#these are for the covariates in best model
recTLh<-matrix(nrow=100, ncol=1000, NA)
recVRMGO2h<-matrix(nrow=10000, ncol=1000, NA)
recTemph<-matrix(nrow=3, ncol=1000, NA)
set.seed(3)
for(i in 1:1000){
  df=cobia_hypox[sample(nrow(cobia_hypox),nrow(cobia_hypox),replace=TRUE),]
  h9a<-lme(finalMO2~VentMouthGape_Metric*midO2+Test_Temp+TL,random=~1|Animal_ID, data=df,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="REML")
  #for TL
  newdatTL.hypox<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_Metric=mean(df$VentMouthGape_Metric),midO2=mean(df$midO2),
                              TL=seq(from=min(df$TL,na.rm = T),to=max(df$TL,na.rm = T),length=100), 
                              Animal_ID=levels(as.factor(as.character(df$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
  pTL.hypox<-cbind(newdatTL.hypox,pred=predict(h9a, newdata = newdatTL.hypox,se.fit=TRUE, type="response"))
  pTL.hypoxTL<-summaryBy(pred~TL, data=pTL.hypox,FUN=sumfun) #sometimes NAs fill these values because if the random dataset that is generated doesn't have all individuals randomly selected than when trying to predict over that ind the model will give an NA (happens sometimes for inds with little number of points)
  recTLh[1:100,i]<-pTL.hypoxTL$pred.mean
  #for temp
  newdatTemp.hypox<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_Metric=mean(df$VentMouthGape_Metric),midO2=mean(df$midO2),
                                TL=mean(df$TL), 
                                Animal_ID=levels(as.factor(as.character(df$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
  pTemp.hypox<-cbind(newdatTemp.hypox,pred=predict(h9a, newdata = newdatTemp.hypox,se.fit=TRUE, type="response"))
  pTemp.hypoxT<-summaryBy(pred~Test_Temp, data=pTemp.hypox,FUN=sumfun)
  recTemph[1:3,i]<-pTemp.hypoxT$pred.mean
  #for VRMG
  newdatVRMG.hypox<-expand.grid(Test_Temp=levels(df$Test_Temp),TL=mean(df$TL),Animal_ID=levels(as.factor(as.character(df$Animal_ID))),
                                VentMouthGape_Metric=seq(from=min(df$VentMouthGape_Metric,na.rm = T),to=max(df$VentMouthGape_Metric,na.rm = T),length=100), 
                                midO2=seq(from=min(df$midO2,na.rm = T),to=max(df$midO2,na.rm = T),length=100))
  pVRMG.hypox<-cbind(newdatVRMG.hypox,pred=predict(h9a, newdata = newdatVRMG.hypox,se.fit=TRUE, type="response"))
  pVRMG.hypoxV<-summaryBy(pred~VentMouthGape_Metric+midO2, data=pVRMG.hypox,FUN=sumfun)
  recVRMGO2h[1:10000,i]<-pVRMG.hypoxV$pred.mean
  
  print(i)
}

#calculating confidence intervals from bootstrapping
CIrec<-NULL

outVRMGO2h<-CIfunc(recVRMGO2h)
outTemph<-CIfunc(recTemph)
outTLh<-CIfunc(recTLh)


####predicting to get actual model output (non bootstrapping)
#rerun model so mod3a can be used to get actual predictions again (model was slightly different above due to boostrapping)
h9a<-lme(finalMO2~VentMouthGape_Metric*midO2+Test_Temp+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
         correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="REML")

#for TL
newdatTL.hypox<-expand.grid(Test_Temp=levels(cobia_hypox$Test_Temp),VentMouthGape_Metric=mean(cobia_hypox$VentMouthGape_Metric),midO2=mean(cobia_hypox$midO2),
                           TL=seq(from=min(cobia_hypox$TL,na.rm = T),to=max(cobia_hypox$TL,na.rm = T),length=100), 
                           Animal_ID=levels(as.factor(as.character(cobia_hypox$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTL.hypox<-cbind(newdatTL.hypox,pred=predict(h9a, newdata = newdatTL.hypox,se.fit=TRUE, type="response"))
pTL.hypoxTL<-summaryBy(pred~TL, data=pTL.hypox)
pTL.hypoxTL$UpperCI<-pTL.hypoxTL$pred.mean+outTLh #add CI to mean
pTL.hypoxTL$LowerCI<-pTL.hypoxTL$pred.mean-outTLh #subtract CI from mean

#for temp
newdatTemp.hypox<-expand.grid(Test_Temp=levels(cobia_hypox$Test_Temp),VentMouthGape_Metric=mean(cobia_hypox$VentMouthGape_Metric),midO2=mean(cobia_hypox$midO2),
                             TL=mean(cobia_hypox$TL), 
                             Animal_ID=levels(as.factor(as.character(cobia_hypox$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTemp.hypox<-cbind(newdatTemp.hypox,pred=predict(h9a, newdata = newdatTemp.hypox,se.fit=TRUE, type="response"))
pTemp.hypoxT<-summaryBy(pred~Test_Temp, data=pTemp.hypox)
pTemp.hypoxT$UpperCI<-pTemp.hypoxT$pred.mean+outTemph #add CI to mean
pTemp.hypoxT$LowerCI<-pTemp.hypoxT$pred.mean-outTemph #subtract CI from mean

#for VRMG
newdatVRMG.hypox<-expand.grid(Test_Temp=levels(cobia_hypox$Test_Temp),TL=mean(cobia_hypox$TL),Animal_ID=levels(as.factor(as.character(cobia_hypox$Animal_ID))),
                             VentMouthGape_Metric=seq(from=min(cobia_hypox$VentMouthGape_Metric,na.rm = T),to=max(cobia_hypox$VentMouthGape_Metric,na.rm = T),length=100), 
                             midO2=seq(from=min(cobia_hypox$midO2,na.rm = T),to=max(cobia_hypox$midO2,na.rm = T),length=100))
pVRMG.hypox<-cbind(newdatVRMG.hypox,pred=predict(h9a, newdata = newdatVRMG.hypox,se.fit=TRUE, type="response"))
pVRMGO2.hypoxV<-summaryBy(pred~VentMouthGape_Metric+midO2, data=pVRMG.hypox)


#####
#Plotting
plot(finalMO2~TL, data=cobia_hypox,pch=21, bg=c("blue","dark grey","red")[cobia_hypox$Test_Temp], ylim=c(min(cobia_hypox$finalMO2),max(cobia_hypox$finalMO2)),
     ylab=expression(paste("Metabolic Rate (mg kg"^"-1","h"^"-1",")")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)
#legend("topleft",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pTL.hypoxTL$TL,rev(pTL.hypoxTL$TL)), y=c(pTL.hypoxTL$LowerCI,rev(pTL.hypoxTL$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~TL,data=pTL.hypoxTL, lwd=2, pch=16)

barplot<-barplot(height=pTemp.hypoxT$pred.mean, names.arg=pTemp.hypoxT$Test_Temp,ylab="Metabolic Rate (mg/kg/h)",xlab="Temperature (°C)", ylim=c(0,200))
arrows(y0=pTemp.hypoxT$LowerCI,y1=pTemp.hypoxT$UpperCI, 
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)


#need to get data in matrix
mat<-matrix(NA,nrow=100, ncol=100)
VRMG_loop<-unique(pVRMGO2.hypoxV$VentMouthGape_Metric)
VRMG_loop<-sort(VRMG_loop, decreasing = FALSE)
O2_loops<-unique(pVRMGO2.hypoxV$midO2)
O2_loops<-sort(O2_loops, decreasing = FALSE)
for(i in 1:length(VRMG_loop)){
  o<-which(pVRMGO2.hypoxV$VentMouthGape_Metric==VRMG_loop[i])
  o1<-pVRMGO2.hypoxV[o,]
  for(j in 1:length(O2_loops)){
    ox<-which(o1$midO2==O2_loops[j])
    ox1<-o1[ox,]
    mat[i,j]<-ox1$pred.mean
  }
}

matT<-t(mat)

O2_loops1<-O2_loops*32/1000
cobia_hypox$midO21<-cobia_hypox$midO2*32/1000

#plotting (HIGH RES FIG)
library(fields)
#tiff("test1.tiff", width=25, height=15, units="cm",compression="lzw",res=300)
par(mar=c(6,5,2,3))
image.plot(O2_loops1,VRMG_loop,matT,zlim=c(45, 400), xlab=expression(paste("Oxygen Concentration (mg L"^"-1",")")),ylab="VRMG", axes=T,col=tim.colors(390),cex.axis=1.25,cex.lab=1.25,
           legend.lab=expression(paste("Metabolic Rate (mg kg"^"-1","h"^"-1",")")),legend.line=3,legend.mar=8.5,legend.cex=1.25)
points(VentMouthGape_Metric~midO21,data=cobia_hypox,pch=c(21,25)[cobia_hypox$Crash],bg=color.scale(cobia_hypox$finalMO2, col= tim.colors(390),zlim=c(45, 400)))
legend("bottomleft",legend=c("Pre Ccrit","Post Ccrit"),pch=c(25,21),col=c("black","black"))
#dev.off()






















##############
#Normalized VRMG (to max)
##############



setwd("~/Documents/PhD Project/Respirometry/Respirometry")
sumdat<-read.csv("MO2_methods10_cobia_20180914.csv")

setwd("~/Documents/PhD Project/Ventilation")

vent<-read.csv("cobia_ventilation.csv")
vent$Date_Time<-paste(vent$Date,vent$time, sep=" ")
vent$Date_Time<-as.POSIXct(vent$Date_Time, format="%m/%d/%y %H:%M:%S")

names(vent)[5]<-"Ventilation_Rate"

library(StreamMetabolism)

#read in all csvs ending with "data_corrected.csv" and putting csvs in list
mycsv = dir(pattern="data_corrected.csv")
n <- length(mycsv)
mylist <- vector("list", n)
for(i in 1:n) mylist[[i]] <- read.csv(mycsv[i])


datRec<-NULL
for(q in 1:length(mylist)){
  dat<-mylist[q]
  dat<-data.frame(dat)
  dat$Date_Time<-paste(dat$Date, dat$time, seq=" ")
  dat$Date_Time<-as.POSIXct(dat$Date_Time, format="%m/%d/%y %H:%M:%S")
  #dat<-dat[which(dat$r2>=0.80),] #gets rid of points where r2 is lower than .80
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
  
  vent_dat<-vent[which(as.character(vent$Animal_ID)==as.character(dat$Animal_ID[1]) & vent$Test_Temp==dat$Test_Temp[1]),]
  
  merged<-merge(dat,vent_dat,by="Date_Time",all.x=T)
  
  dat<-cbind(dat,Time.of.Gape=merged$Time.of.Gape,Ventilation_Rate=merged$Ventilation_Rate,Mouth_Gape=merged$Mouth.Gape)
  
  #want to divide mouth gape by minimum mouth gape during trial so that when mouth gape and vent rate are multiplied
  #the VRMGM is equal to just the ventilation rate at the smallest mouth gape, then as mouth gape increases so will its impact on ventilation rate and this VRMGM
  dat$Mouth_Gape_to_min<-dat$Mouth_Gape/min(dat$Mouth_Gape,na.rm=T)
  
  dat$VentMouthGape_Metric<-dat$Ventilation_Rate*dat$Mouth_Gape_to_min
  dat$VentMouthGape_MetricMax<-dat$VentMouthGape_Metric/max(dat$VentMouthGape_Metric,na.rm=T)
  
  #plot entire trial, remember not all trials have two nights so will see warnings for some trials
  setwd("~/Documents/PhD Project/Ventilation/Trials/Ventilation_to_max_figures")
  #png(file=paste(dat$Animal_ID[1], dat$Test_Temp[1],"ventilation","png", sep="."),width=650,height=450)
  tiff(file=paste(dat$Animal_ID[1], dat$Test_Temp[1],"ventilation.to.max","tiff", sep="."), width=30, height=15, units="cm",compression="lzw",res=150)
  par(mfrow=c(1,1))
  par(mar=c(6,7.3,2,10))
  plot(finalMO2~Date_Timen, main=paste(dat$Animal_ID[1], dat$Test_Temp[1],sep="_") ,data=dat, xaxt="n", type="b",lwd=2, ylab="", xlab="")
  r<-as.POSIXct(round(range(dat$Date_Time), "hours"))
  axis.POSIXct(1, at = seq(r[1], r[2], by = 7200), format = "%H:%M",las=2)
  mtext(side=1,line=4,"Time of day",cex=1.25)
  mtext(side=2,line=2,expression(paste("Metabolic Rate (mg kg"^"-1","h"^"-1",")")),cex=1.25)
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
  par(new=T)
  plot(Ventilation_Rate~Date_Timen, data=dat, type="l", xaxt="n", axes=F, xlab="",ylab="",lwd=1,lty=2,col="blue")
  axis(side=4, line=3.75)
  mtext(side=4,line=5.75,"Ventilation Rate")
  par(new=T)
  plot(VentMouthGape_MetricMax~Date_Timen, data=dat, type="b", xaxt="n", axes=F, xlab="",ylab="",lwd=1,lty=2,col="deeppink2")
  axis(side=2, line=3.75)
  mtext(side=2,line=6,"VRMG") #for Vent Rate Mouth Gape Metric
  par(new=T)
  plot(Mouth_Gape~Date_Timen, data=dat, type="l", xaxt="n", axes=F, xlab="",ylab="",lwd=1,lty=2,col="dark green")
  axis(side=4, line=7)
  mtext(side=4,line=9,"Mouth Gape (cm)")
  
  dev.off()
  
  SMR<-sumdat$SMR[which(as.character(sumdat$Animal_ID)==as.character(dat$Animal_ID[1]) & as.character(sumdat$Test_Temp)==as.character(dat$Test_Temp[1]))]
  SMRfirst<-min(which(dat$finalMO2<SMR))
  dat$before_after_SMR<-NA
  dat$before_after_SMR[1:(SMRfirst-1)]<-"before"
  dat$before_after_SMR[SMRfirst:nrow(dat)]<-"after"
  
  dat<-dat[which(dat$r2>=0.80),]
  
  datRec<-rbind(datRec,dat)
  print(q)
}

#add TL to dataframe
datRec$TL<-NA
datRec$TL[which(datRec$Animal_ID=="C08")]<-95
datRec$TL[which(datRec$Animal_ID=="C09")]<-88
datRec$TL[which(datRec$Animal_ID=="C11")]<-86
datRec$TL[which(datRec$Animal_ID=="C12")]<-105
datRec$TL[which(datRec$Animal_ID=="C13")]<-108
datRec$TL[which(datRec$Animal_ID=="C15")]<-84
datRec$TL[which(datRec$Animal_ID=="C16")]<-87
datRec$TL[which(datRec$Animal_ID=="C17")]<-88
datRec$TL[which(datRec$Animal_ID=="C18")]<-92
datRec$TL[which(datRec$Animal_ID=="C19")]<-91.5
datRec$TL[which(datRec$Animal_ID=="C20")]<-87.5
datRec$TL[which(datRec$Animal_ID=="C21")]<-82
datRec$TL[which(datRec$Animal_ID=="C22")]<-96
datRec$TL[which(datRec$Animal_ID=="C23")]<-83
datRec$TL[which(datRec$Animal_ID=="C24")]<-81
datRec$TL[which(datRec$Animal_ID=="C25")]<-88
datRec$TL[which(datRec$Animal_ID=="C26")]<-108
datRec$TL[which(datRec$Animal_ID=="C27")]<-88.5
datRec$TL[which(datRec$Animal_ID=="C28")]<-102
datRec$TL[which(datRec$Animal_ID=="C30")]<-81
datRec$TL[which(datRec$Animal_ID=="C31")]<-85
datRec$TL[which(datRec$Animal_ID=="C32")]<-81
datRec$TL[which(datRec$Animal_ID=="C33")]<-96
datRec$TL[which(datRec$Animal_ID=="C34")]<-89
datRec$TL[which(datRec$Animal_ID=="C35")]<-102.5

#remember not including inds that died (C31_28, C24_32, C32_32, C28_32)
datRec$Alive_Dead<-ifelse(datRec$Animal_ID=="C31" & datRec$Test_Temp=="28", "dead",
                          ifelse(datRec$Animal_ID=="C24" & datRec$Test_Temp=="32", "dead",
                                 ifelse(datRec$Animal_ID=="C32" & datRec$Test_Temp=="32", "dead",
                                        ifelse(datRec$Animal_ID=="C28" & datRec$Test_Temp=="32","dead","alive"))))

datRec_dead<-datRec[which(datRec$Alive_Dead=="dead"),]
datRec<-datRec[which(datRec$Alive_Dead=="alive"),]

#####add Crit and Scrit to cobia dataframe
#need to pull it from cobia_parts_sas
setwd("~/Documents/PhD Project/Respirometry/Respirometry")
ccritdat<-read.csv("MO2_methods10_cobia_20180914.csv")

#assign NAs to Ccrit and Scrit for 30_28 and 34_28 because the automation method to calculate ccrit from SMR was not accurate for those two fish
ccritdat$Ccrit[which(ccritdat$Animal_ID=="C30"&ccritdat$Test_Temp=="28")]<-NA
ccritdat$Scrit[which(ccritdat$Animal_ID=="C30"&ccritdat$Test_Temp=="28")]<-NA
ccritdat$Ccrit[which(ccritdat$Animal_ID=="C34"&ccritdat$Test_Temp=="28")]<-NA
ccritdat$Scrit[which(ccritdat$Animal_ID=="C34"&ccritdat$Test_Temp=="28")]<-NA

ccritdat<-subset(ccritdat, select=c(Animal_ID, Test_Temp, Ccrit))

ccritdat$Ccritnew<-ccritdat$Ccrit*1000/32#puts ccrit into same units as O2 min units (umol/l)

#Test_Temp should be factor
datRec$Test_Temp<-as.factor(datRec$Test_Temp)

#going to use O2 middway between min and max
datRec$midO2<-(datRec$max..O2.+datRec$min..O2.)/2

#add observation column
animal<-unique(datRec$Animal_ID)
temp<-unique(datRec$Test_Temp)
obs_list<-NULL
test_list<-NULL
crash<-NULL
for(i in animal){
  a<-which(datRec$Animal_ID==i)
  a1<-datRec[a,]
  s<-which(ccritdat$Animal_ID==i)
  s1<-ccritdat[s,]
  for(j in temp){
    t<-which(a1$Test_Temp==j)
    t1<-a1[t,]
    p<-which(s1$Test_Temp==j)
    p1<-s1[p,]
    obs<-seq(from=1, to=nrow(t1))
    obs_list<-c(obs_list,obs)
    ccrit<-p1$Ccritnew[1]
    b_a<-ifelse(t1$min..O2.>p1$Ccritnew,"before","after")
    crash<-c(crash, b_a)
  }
  test<-seq(from=1, to=nrow(a1))
  test_list<-c(test_list,test)
}

badlocs<-which(obs_list==0)
badlocs1<-badlocs+1
locs<-c(badlocs,badlocs1)
obs_list<-obs_list[-locs]

datRec$Obs<-obs_list
datRec$test<-test_list
datRec$Crash<-crash


datRec<-subset(datRec, select=c(Date, time, animal.mass,midO2,Animal_ID,Test_Temp,X100_Perc_O2,
                                Date_Time,Time_Diff,finalMO2,Day_Night, hypoxia,Ventilation_Rate,Mouth_Gape,
                                VentMouthGape_Metric,VentMouthGape_MetricMax, TL, Obs, test, Crash,before_after_SMR))

#use Time_Diff as proxy for time after chase
datRec$Time_Diff<-datRec$Time_Diff+60 #added 60 seconds because it is roughly 60s from chase to start of trial, will use time_diff as proxy for time after chase
datRec$Crash<-as.factor(datRec$Crash)

#need to get rid of NA rows
cobia<-datRec[complete.cases(datRec$VentMouthGape_Metric),]#will have a lot of eliminated for when fish wasn't facing camera

#######################
###Normoxia analysis (VRMG normalized)
######################

cobia_norm<-cobia[which(cobia$hypoxia=="no"),]
cobia_norm$Animal_ID<-as.factor(as.character(cobia_norm$Animal_ID)) #do this to get rid of 2 cobia that died that were only tested at one temp to get rid of that level
cobia_norm$before_after_SMR<-as.factor(cobia_norm$before_after_SMR)

####Check for correlations in covariates
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(cbind(cobia_norm$Test_Temp,cobia_norm$Ventilation_Rate,cobia_norm$Mouth_Gape,cobia_norm$VentMouthGape_MetricMax,
            cobia_norm$Time_Diff,cobia_norm$animal.mass,cobia_norm$TL,cobia_norm$before_after_SMR),lower.panel = panel.cor)
#animal.mass and TL are correlated so I took animal.mass out
#time_diff and before_after_SMR are correlated, use time_diff
#everything else seems okay!


#######modeling
library(nlme)
library(HH)
mod0<-lm(finalMO2~Ventilation_Rate+Mouth_Gape+VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff, data=cobia_norm,na.action = na.omit,x=TRUE)
vif(mod0) #covariates are not collinear, below 5 is okay

#no correlation structure and not modeling variance structure (heterogenity)
mod1<-lme(finalMO2~Ventilation_Rate+Mouth_Gape+VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit, method = "ML")
summary(mod1)
#diagnostics
plot(mod1)
plot(resid(mod1)~cobia_norm$Test_Temp)
plot(resid(mod1)~cobia_norm$TL)
plot(resid(mod1)~cobia_norm$VentMouthGape_MetricMax)
#equal variance looks slightly skewed in VentMouthGape_MetricMax, but can't seem to fix it but messing with weights, will ignore for now
qqnorm(mod1) #eh not too good
mod1.2<-lme(log(finalMO2)~Ventilation_Rate+Mouth_Gape+VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit, method = "ML")
plot(mod1.2) #much better now
plot(resid(mod1.2)~cobia_norm$Test_Temp)
plot(resid(mod1.2)~cobia_norm$TL)
plot(resid(mod1.2)~cobia_norm$VentMouthGape_MetricMax)
library(car)
vif(mod1.2) #looks like VRMG and Mouth Gape are colinear, will get rid of Mouth Gape and Ventilation rate from model
mod1.3<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit, method = "ML")
plot(mod1.3) #much better now
plot(resid(mod1.3)~cobia_norm$Test_Temp)
plot(resid(mod1.3)~cobia_norm$TL)
plot(resid(mod1.3)~cobia_norm$VentMouthGape_MetricMax)
qqnorm(mod1.3)#better
vif(mod1.3) #better
#have this just so have model not log transformed, doesn't seem as good as log transformed one
mod1.4<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit, method = "ML")

#check for temporal autocorrelation (ie. independence vs non independence)
E<-residuals(mod1.3, type="normalized")
I1<-!is.na(cobia_norm$VentMouthGape_MetricMax)
Efull<-vector(length=length(cobia_norm$Test_Temp))
Efull<-NA
Efull[I1]<-E
acf(Efull,na.action = na.pass) #definitely independence violated

#no correlation structure, model variance structure (heterogenity)
mod2<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,weights = varIdent(form=~1|Test_Temp), method="ML")
BIC(mod1.3,mod2) #model with variance structure is better

#test different correlation structures
mod3<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod4<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corCompSymm(form=~Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
#mod4 not converging anyway
BIC(mod1.3,mod2,mod3,mod4)#model with AR1 correlation structure is best

#Now adjust fixed covariates
library(mgcv)
mod3<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod5<-lme(log(finalMO2)~VentMouthGape_MetricMax*TL+Test_Temp+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod6<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod7<-lme(log(finalMO2)~VentMouthGape_MetricMax*Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod8<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")
mod9<-gam(log(finalMO2)~s(VentMouthGape_MetricMax)+Test_Temp+s(Time_Diff)+s(Animal_ID,bs="re"), data=cobia_norm,na.action = na.omit,
          correlation = corAR1(form=~1|Animal_ID),method="ML",family=gaussian(link = "identity"))
mod10<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="ML")

BIC(mod3,mod5,mod6,mod7,mod8,mod9,mod10) #best model is mod3 (BIC is -1498.6584)


#switch to reml
mod3a<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="REML")
summary(mod3a)
#VRMG is sig
#TL is sig
#Time_Diff is sig
anova(mod3a)
#temp is sig

library(multcomp)
library(lsmeans)
#look at interaction between TBF and Temp
lsmeans(mod3a, pairwise~Test_Temp, adjust="tukey")
#all are sig different from each other

#rsquared calc 
#conditional R2: describes the proportion of variance explained by both the fixed and random factors
#marginal R2: describes the proportion of variance explained by the fixed factor(s) alone
#either package will do
library(piecewiseSEM) #John Lefcheck's package
rsquared(mod3a)
library(gabtool)
r.squared(mod3a)
#conditional R2 is 0.71 which is very good

library(doBy)

disp.lme=function(mod=NA){
  E=resid(mod,type='pearson')
  d=sum(E^2)/mod$fixDF$X[1]
  return(d)
}

#try to sort out log space and random effects...couldn't really find a solution so just left data in log space for figures


######
#Bootstrapping to get errors for plots
#these are for the covariates in best model
sumfun <- function(x, ...){ #this function fixes the problem of when NAs are introduced before summaryBy, get rid of NAs before mean is taken, NAs are made when during bootstrap no points from an ind is selected so would get an NA when trying to predict over that ind
  c(mean=mean(x, na.rm=TRUE, ...), var=var(x, na.rm=TRUE, ...), length=length(x))
}
recVRMG<-matrix(nrow=100, ncol=1000, NA)
recTD<-matrix(nrow=100, ncol=1000, NA)
recTL<-matrix(nrow=100, ncol=1000, NA)
recTemp<-matrix(nrow=3, ncol=1000, NA)
set.seed(2)
for(i in 1:1000){
  df=cobia_norm[sample(nrow(cobia_norm),nrow(cobia_norm),replace=TRUE),]
  mod3a<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=df,na.action = na.omit,
             correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="REML")
  #VRMG
  newdatVRMG.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),Time_Diff=mean(df$Time_Diff),TL=mean(df$TL),
                               VentMouthGape_MetricMax=seq(from=min(df$VentMouthGape_MetricMax,na.rm = T),to=max(df$VentMouthGape_MetricMax,na.rm = T),length=100), 
                               Animal_ID=levels(df$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
  pVRMG.norm<-cbind(newdatVRMG.norm,pred=predict(mod3a, newdata = newdatVRMG.norm,se.fit=TRUE, type="response"))
  #pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred+(0.5*disp.lme(mod3a)))
  pVRMG.normV<-summaryBy(pred~VentMouthGape_MetricMax, data=pVRMG.norm)
  recVRMG[1:100,i]<-pVRMG.normV$pred.mean
  #Time_Diff
  newdatTD.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_MetricMax=mean(df$VentMouthGape_MetricMax),TL=mean(df$TL),
                             Time_Diff=seq(from=min(df$Time_Diff,na.rm = T),to=max(df$Time_Diff,na.rm = T),length=100), 
                             Animal_ID=levels(df$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
  pTD.norm<-cbind(newdatTD.norm,pred=predict(mod3a, newdata = newdatTD.norm,se.fit=TRUE, type="response"))
  #pTD.norm$pred.corrected<-exp(pTD.norm$pred+(0.5*disp.lme(mod3a)))
  pTD.normTD<-summaryBy(pred~Time_Diff, data=pTD.norm)
  recTD[1:100,i]<-pTD.normTD$pred.mean
  #TL
  newdatTL.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_MetricMax=mean(df$VentMouthGape_MetricMax),Time_Diff=mean(df$Time_Diff),
                             TL=seq(from=min(df$TL,na.rm = T),to=max(df$TL,na.rm = T),length=100), 
                             Animal_ID=levels(as.factor(as.character(df$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
  pTL.norm<-cbind(newdatTL.norm,pred=predict(mod3a, newdata = newdatTL.norm,se.fit=TRUE, type="response"))
  #pTL.norm$pred.corrected<-exp(pTL.norm$pred+(0.5*disp.lme(mod3a)))
  pTL.normTL<-summaryBy(pred~TL, data=pTL.norm)
  recTL[1:100,i]<-pTL.normTL$pred.mean
  #Temp
  newdatTemp.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_MetricMax=mean(df$VentMouthGape_MetricMax),Time_Diff=mean(df$Time_Diff),
                               TL=mean(df$TL), 
                               Animal_ID=levels(df$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
  pTemp.norm<-cbind(newdatTemp.norm,pred=predict(mod3a, newdata = newdatTemp.norm,se.fit=TRUE, type="response"))
  #pTemp.norm$pred.corrected<-exp(pTemp.norm$pred+(0.5*disp.lme(mod3a)))
  pTemp.normT<-summaryBy(pred~Test_Temp, data=pTemp.norm)
  recTemp[1:3,i]<-pTemp.normT$pred.mean
  
  print(i)
}

#calculating confidence intervals from bootstrapping
CIrec<-NULL
CIfunc<-function(data=NA){
  
  for(i in 1:nrow(data)){
    CI<-1.96*sd(data[i,],na.rm=T)
    CIrec<-c(CIrec,CI)
  }
  return(CIrec)
}

outVRMG<-CIfunc(recVRMG)
outTD<-CIfunc(recTD)
outTL<-CIfunc(recTL)
outTemp<-CIfunc(recTemp)


####predicting to get actual model output (non bootstrapping)
#rerun model so mod3a can be used to get actual predictions again (model was slightly different above due to boostrapping)
mod3a<-lme(log(finalMO2)~VentMouthGape_MetricMax+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="REML")

#for VRMG
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Time_Diff=mean(cobia_norm$Time_Diff),TL=mean(cobia_norm$TL),
                             VentMouthGape_MetricMax=seq(from=min(cobia_norm$VentMouthGape_MetricMax,na.rm = T),to=max(cobia_norm$VentMouthGape_MetricMax,na.rm = T),length=100), 
                             Animal_ID=levels(cobia_norm$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
pVRMG.norm<-cbind(newdatVRMG.norm,pred=predict(mod3a, newdata = newdatVRMG.norm,se.fit=TRUE, type="response"))
#pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred+(0.5*disp.lme(mod3a)))
pVRMG.normV<-summaryBy(pred~VentMouthGape_MetricMax, data=pVRMG.norm)
pVRMG.normV$UpperCI<-pVRMG.normV$pred.mean+outVRMG #add CI to mean
pVRMG.normV$LowerCI<-pVRMG.normV$pred.mean-outVRMG #subtract CI from mean

newdatTD.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),VentMouthGape_MetricMax=mean(cobia_norm$VentMouthGape_MetricMax),TL=mean(cobia_norm$TL),
                           Time_Diff=seq(from=min(cobia_norm$Time_Diff,na.rm = T),to=max(cobia_norm$Time_Diff,na.rm = T),length=100), 
                           Animal_ID=levels(as.factor(as.character(cobia_norm$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTD.norm<-cbind(newdatTD.norm,pred=predict(mod3a, newdata = newdatTD.norm,se.fit=TRUE, type="response"))
#pTD.norm$pred.corrected<-exp(pTD.norm$pred+(0.5*disp.lme(mod3a)))
pTD.normTD<-summaryBy(pred~Time_Diff, data=pTD.norm)
pTD.normTD$UpperCI<-pTD.normTD$pred.mean+outTD #add CI to mean
pTD.normTD$LowerCI<-pTD.normTD$pred.mean-outTD #subtract CI from mean

newdatTL.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),VentMouthGape_MetricMax=mean(cobia_norm$VentMouthGape_MetricMax),Time_Diff=mean(cobia_norm$Time_Diff),
                           TL=seq(from=min(cobia_norm$TL,na.rm = T),to=max(cobia_norm$TL,na.rm = T),length=100), 
                           Animal_ID=levels(as.factor(as.character(cobia_norm$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTL.norm<-cbind(newdatTL.norm,pred=predict(mod3a, newdata = newdatTL.norm,se.fit=TRUE, type="response"))
#pTL.norm$pred.corrected<-exp(pTL.norm$pred+(0.5*disp.lme(mod3a)))
pTL.normTL<-summaryBy(pred~TL, data=pTL.norm)
pTL.normTL$UpperCI<-pTL.normTL$pred.mean+outTL #add CI to mean
pTL.normTL$LowerCI<-pTL.normTL$pred.mean-outTL #subtract CI from mean

newdatTemp.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),VentMouthGape_MetricMax=mean(cobia_norm$VentMouthGape_MetricMax),Time_Diff=mean(cobia_norm$Time_Diff),
                             TL=mean(cobia_norm$TL), 
                             Animal_ID=levels(as.factor(as.character(cobia_norm$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTemp.norm<-cbind(newdatTemp.norm,pred=predict(mod3a, newdata = newdatTemp.norm,se.fit=TRUE, type="response"))
#pTemp.norm$pred.corrected<-exp(pTemp.norm$pred+(0.5*disp.lme(mod3a)))
pTemp.normT<-summaryBy(pred~Test_Temp, data=pTemp.norm)
pTemp.normT$UpperCI<-pTemp.normT$pred.mean+outTemp #add CI to mean
pTemp.normT$LowerCI<-pTemp.normT$pred.mean-outTemp #subtract CI from mean


#plot covariates
#VRMG
par(mar=(c(5,4.5,4,1)))
plot(log(finalMO2)~VentMouthGape_MetricMax, data=cobia_norm,pch=21, bg=c("blue","dark grey","red")[cobia_norm$Test_Temp], ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)
legend("topleft",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pVRMG.normV$VentMouthGape_MetricMax,rev(pVRMG.normV$VentMouthGape_MetricMax)), y=c(pVRMG.normV$LowerCI,rev(pVRMG.normV$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~VentMouthGape_MetricMax,data=pVRMG.normV, lwd=2, pch=16)

#VRMG: plot with plateaus from segmented by temp
library(segmented)
poly24<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="24"),],na.action = na.omit)
my.seg24 <- segmented(poly24, 
                    seg.Z = ~ VentMouthGape_MetricMax, 
                    psi = 0.8)
poly28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"),],na.action = na.omit)
my.seg28 <- segmented(poly28, 
                      seg.Z = ~ VentMouthGape_MetricMax, 
                      psi = 0.6)
poly32<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="32"),],na.action = na.omit)
my.seg32 <- segmented(poly32, 
                      seg.Z = ~ VentMouthGape_MetricMax, 
                      psi = 0.8)

par(mar=(c(5,4.5,4,1)))
plot(log(finalMO2)~VentMouthGape_MetricMax, data=cobia_norm,pch=21, bg=c("blue","dark grey","red")[cobia_norm$Test_Temp], ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)
legend("topleft",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
plot(my.seg24,add=TRUE,link=FALSE,lwd=2,col="blue")
plot(my.seg28,add=TRUE,link=FALSE,lwd=2,col="grey44")
plot(my.seg32,add=TRUE,link=FALSE,lwd=2,col="red")
polygon(x=c(pVRMG.normV$VentMouthGape_MetricMax,rev(pVRMG.normV$VentMouthGape_MetricMax)), y=c(pVRMG.normV$LowerCI,rev(pVRMG.normV$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~VentMouthGape_MetricMax,data=pVRMG.normV, lwd=3, pch=16)

#Time Diff
plot(log(finalMO2)~Time_Diff, data=cobia_norm,pch=21, bg=c("blue","dark grey","red")[cobia_norm$Test_Temp], ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="Time since chase (s)",cex.axis=1.25,cex.lab=1.25)
legend("topright",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pTD.normTD$Time_Diff,rev(pTD.normTD$Time_Diff)), y=c(pTD.normTD$LowerCI,rev(pTD.normTD$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~Time_Diff,data=pTD.normTD, lwd=2, pch=16)

#Total Length
plot(log(finalMO2)~TL, data=cobia_norm,pch=21, bg=c("blue","dark grey","red")[cobia_norm$Test_Temp], ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)
legend("topright",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pTL.normTL$TL,rev(pTL.normTL$TL)), y=c(pTL.normTL$LowerCI,rev(pTL.normTL$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~TL,data=pTL.normTL, lwd=2, pch=16)

#Temp
barplot<-barplot(height=pTemp.normT$pred.mean, names.arg=pTemp.normT$Test_Temp,ylab="log[Metabolic Rate (mg/kg/h)]",xlab="Temperature (°C)", ylim=c(0,5.5))
arrows(y0=pTemp.normT$LowerCI,y1=pTemp.normT$UpperCI, 
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)

#for showing that raw data based on chase, non chase
cobia_norm$chase_time<-as.factor(ifelse(cobia_norm$Time_Diff<=10800,"yes","no"))



#VRMG: plotting log(MO2) and VRMG for each trial
#can't use all trials because don't have enough VRMG values for all trials (threshold will be 30% of total number of points wthin trial)
propnrowRec<-NULL
cobia_norm$Trial<-paste(cobia_norm$Animal_ID,cobia_norm$Test_Temp,sep="_") #dataset with just VRMG
datRec$Trial<-paste(datRec$Animal_ID,datRec$Test_Temp,sep="_") #full dataset
trial<-unique(cobia_norm$Trial)
for(i in trial){
  tr<-which(cobia_norm$Trial==i)
  tr1<-cobia_norm[tr,]
  trdat<-which(datRec$Trial==i)
  trdat1<-datRec[trdat,]
  nrowtr1<-nrow(tr1)
  nrowdattr1<-nrow(trdat1)
  propnrow<-nrowtr1/nrowdattr1
  propnrowRec<-c(propnrowRec,propnrow)
}
cbind(trial,propnrowRec)


library(segmented)
modC12_24<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="24"&cobia_norm$Animal_ID=="C12"),],na.action = na.omit)
my.segC12_24 <- segmented(modC12_24, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC15_24<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="24"&cobia_norm$Animal_ID=="C15"),],na.action = na.omit)
my.segC15_24 <- segmented(modC15_24, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC21_24<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="24"&cobia_norm$Animal_ID=="C21"),],na.action = na.omit)
my.segC21_24 <- segmented(modC21_24, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC22_24<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="24"&cobia_norm$Animal_ID=="C22"),],na.action = na.omit)
my.segC22_24 <- segmented(modC22_24, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC23_28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"&cobia_norm$Animal_ID=="C23"),],na.action = na.omit)
my.segC23_28 <- segmented(modC23_28, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.6)
modC24_28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"&cobia_norm$Animal_ID=="C24"),],na.action = na.omit)
my.segC24_28 <- segmented(modC24_28, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC25_28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"&cobia_norm$Animal_ID=="C25"),],na.action = na.omit)
my.segC25_28 <- segmented(modC25_28, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.6)
modC26_28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"&cobia_norm$Animal_ID=="C26"),],na.action = na.omit)
my.segC26_28 <- segmented(modC26_28, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.4)
modC28_28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"&cobia_norm$Animal_ID=="C28"),],na.action = na.omit)
my.segC28_28 <- segmented(modC28_28, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.6)
modC30_28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"&cobia_norm$Animal_ID=="C30"),],na.action = na.omit)
my.segC30_28 <- segmented(modC30_28, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.6)
modC33_28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"&cobia_norm$Animal_ID=="C33"),],na.action = na.omit)
my.segC33_28 <- segmented(modC33_28, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.6)
modC34_28<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="28"&cobia_norm$Animal_ID=="C34"),],na.action = na.omit)
my.segC34_28 <- segmented(modC34_28, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC23_32<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="32"&cobia_norm$Animal_ID=="C23"),],na.action = na.omit)
my.segC23_32 <- segmented(modC23_32, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.9)
modC25_32<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="32"&cobia_norm$Animal_ID=="C25"),],na.action = na.omit)
my.segC25_32 <- segmented(modC25_32, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC26_32<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="32"&cobia_norm$Animal_ID=="C26"),],na.action = na.omit)
my.segC26_32 <- segmented(modC26_32, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC30_32<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="32"&cobia_norm$Animal_ID=="C30"),],na.action = na.omit)
my.segC30_32 <- segmented(modC30_32, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC33_32<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="32"&cobia_norm$Animal_ID=="C33"),],na.action = na.omit)
my.segC33_32 <- segmented(modC33_32, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.95)
modC34_32<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="32"&cobia_norm$Animal_ID=="C34"),],na.action = na.omit)
my.segC34_32 <- segmented(modC34_32, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)
modC35_32<-lm(log(finalMO2) ~ VentMouthGape_MetricMax,data=cobia_norm[which(cobia_norm$Test_Temp=="32"&cobia_norm$Animal_ID=="C35"),],na.action = na.omit)
my.segC35_32 <- segmented(modC35_32, seg.Z = ~ VentMouthGape_MetricMax, psi = 0.8)

plot(log(finalMO2)~VentMouthGape_MetricMax, data=cobia_norm[which(cobia_norm$Test_Temp=="32"&cobia_norm$Animal_ID=="C33"),],pch=21, bg="black", ylim=c(min(log(cobia_norm$finalMO2)),max(log(cobia_norm$finalMO2))),
     ylab=expression(paste("log[Metabolic Rate (mg kg"^"-1","h"^"-1",")]")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25,type="n",xlim=c(0.2,1))
legend("topleft",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
plot(my.segC12_24,add=TRUE,link=FALSE,lwd=1,col="blue")
plot(my.segC15_24,add=TRUE,link=FALSE,lwd=1,col="blue")
plot(my.segC21_24,add=TRUE,link=FALSE,lwd=1,col="blue")
plot(my.segC22_24,add=TRUE,link=FALSE,lwd=1,col="blue")
plot(my.segC23_28,add=TRUE,link=FALSE,lwd=1,col="grey")
plot(my.segC24_28,add=TRUE,link=FALSE,lwd=1,col="grey")
plot(my.segC25_28,add=TRUE,link=FALSE,lwd=1,col="grey")#
plot(my.segC26_28,add=TRUE,link=FALSE,lwd=1,col="grey")
plot(my.segC28_28,add=TRUE,link=FALSE,lwd=1,col="grey")
plot(my.segC30_28,add=TRUE,link=FALSE,lwd=1,col="grey")
plot(my.segC33_28,add=TRUE,link=FALSE,lwd=1,col="grey")
plot(my.segC34_28,add=TRUE,link=FALSE,lwd=1,col="grey")
plot(my.segC23_32,add=TRUE,link=FALSE,lwd=1,col="red")#
plot(my.segC25_32,add=TRUE,link=FALSE,lwd=1,col="red")
plot(my.segC26_32,add=TRUE,link=FALSE,lwd=1,col="red")#
plot(my.segC30_32,add=TRUE,link=FALSE,lwd=1,col="red")
plot(my.segC33_32,add=TRUE,link=FALSE,lwd=1,col="red")#
plot(my.segC34_32,add=TRUE,link=FALSE,lwd=1,col="red")
plot(my.segC35_32,add=TRUE,link=FALSE,lwd=1,col="red")
polygon(x=c(pVRMG.normV$VentMouthGape_MetricMax,rev(pVRMG.normV$VentMouthGape_MetricMax)), y=c(pVRMG.normV$LowerCI,rev(pVRMG.normV$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~VentMouthGape_MetricMax,data=pVRMG.normV, lwd=3, pch=16)



#################
#Hypoxia analysis (VRMG normalized)
#################

cobia_hypox<-cobia[which(cobia$hypoxia=="yes"),]


#######modeling
library(nlme)
library(HH)
h0<-lm(finalMO2~Ventilation_Rate+Mouth_Gape+VentMouthGape_MetricMax+Test_Temp+TL+midO2, data=cobia_hypox,na.action = na.omit,x=TRUE)
vif(h0) #seems good

#no correlation structure and not modeling variance structure (heterogenity)
h1<-lme(finalMO2~Ventilation_Rate+Mouth_Gape+VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,method="ML")
BIC(h0,h1)
vif(h1) #now mouth gape, and VRMG are colinear, will get rid of ventilation rate and mouth gape from model

h1.1<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,method="ML")
vif(h1.1) #all better

summary(h1.1)
#diagnostics
plot(h1.1)
plot(resid(h1.1)~cobia_hypox$Test_Temp)
plot(resid(h1.1)~cobia_hypox$TL)
plot(resid(h1.1)~cobia_hypox$VentMouthGape_MetricMax)
plot(resid(h1.1)~cobia_hypox$midO2)
#looks pretty homogenous but lets model variance really quick and see if it improves the model

#no correlation structure, model variance structure (heterogenity)
h2<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,weights = varIdent(form=~1|Test_Temp),method="ML")
h2.1<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,weights = varExp(form=~VentMouthGape_MetricMax),method="ML")
h2.2<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,weights = varFixed(~VentMouthGape_MetricMax),method="ML")
h2.3<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,weights = varPower(form=~VentMouthGape_MetricMax),method="ML")
BIC(h1.1,h2,h2.1,h2.2,h2.3) #model with temp in weights is best

#check for temporal autocorrelation (ie. independence vs non independence)
E<-residuals(h2, type="normalized")
I1<-!is.na(cobia_hypox$VentMouthGape_MetricMax)
Efull<-vector(length=length(cobia_hypox$Test_Temp))
Efull<-NA
Efull[I1]<-E
acf(Efull,na.action = na.pass) #independence violated slightly

#test different correlation structures
h3<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h4<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corCompSymm(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
BIC(h2,h3,h4)#model with AR1 correlation structure is best


#Now adjust fixed covariates
h3<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+TL+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h5<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h6<-lme(finalMO2~VentMouthGape_MetricMax*TL+Test_Temp+midO2,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h7<-lme(finalMO2~VentMouthGape_MetricMax+Test_Temp*midO2+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h8<-lme(finalMO2~VentMouthGape_MetricMax*Test_Temp+midO2+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h9<-lme(finalMO2~VentMouthGape_MetricMax*midO2+Test_Temp+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
        correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="ML")
h10<-gam(finalMO2~s(VentMouthGape_MetricMax)+Test_Temp+s(midO2)+TL+s(Animal_ID,bs="re"), data=cobia_hypox,na.action = na.omit,
         correlation = corAR1(form=~1|Animal_ID),method="ML",family=gaussian(link = "identity"))

BIC(h3,h5,h6,h7,h8,h9,h10) #best model is h9 (BIC was 2721.247)

#switch to reml
h9a<-lme(finalMO2~VentMouthGape_MetricMax*midO2+Test_Temp+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
         correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="REML")

summary(h9a)
#TL is sig
#interaction between VRMG and midO2 are sig
anova(h9a)
#temp is sig

#look at pairwise comparision
lsmeans(h9a, pairwise~Test_Temp, adjust="tukey")
#all are sig different from each other

#rsquared calc 
#either package will do
library(piecewiseSEM) #John Lefcheck's package
rsquared(h9a)
library(gabtool)
r.squared(h9a)
#conditional R2 is 0.80 which is very good


######
#Bootstrapping to get errors for plots
#these are for the covariates in best model
recTLh<-matrix(nrow=100, ncol=1000, NA)
recVRMGO2h<-matrix(nrow=10000, ncol=1000, NA)
recTemph<-matrix(nrow=3, ncol=1000, NA)
set.seed(3)
for(i in 1:1000){
  df=cobia_hypox[sample(nrow(cobia_hypox),nrow(cobia_hypox),replace=TRUE),]
  h9a<-lme(finalMO2~VentMouthGape_MetricMax*midO2+Test_Temp+TL,random=~1|Animal_ID, data=df,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="REML")
  #for TL
  newdatTL.hypox<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_MetricMax=mean(df$VentMouthGape_MetricMax),midO2=mean(df$midO2),
                              TL=seq(from=min(df$TL,na.rm = T),to=max(df$TL,na.rm = T),length=100), 
                              Animal_ID=levels(as.factor(as.character(df$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
  pTL.hypox<-cbind(newdatTL.hypox,pred=predict(h9a, newdata = newdatTL.hypox,se.fit=TRUE, type="response"))
  pTL.hypoxTL<-summaryBy(pred~TL, data=pTL.hypox,FUN=sumfun) #sometimes NAs fill these values because if the random dataset that is generated doesn't have all individuals randomly selected than when trying to predict over that ind the model will give an NA (happens sometimes for inds with little number of points)
  recTLh[1:100,i]<-pTL.hypoxTL$pred.mean
  #for temp
  newdatTemp.hypox<-expand.grid(Test_Temp=levels(df$Test_Temp),VentMouthGape_MetricMax=mean(df$VentMouthGape_MetricMax),midO2=mean(df$midO2),
                                TL=mean(df$TL), 
                                Animal_ID=levels(as.factor(as.character(df$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
  pTemp.hypox<-cbind(newdatTemp.hypox,pred=predict(h9a, newdata = newdatTemp.hypox,se.fit=TRUE, type="response"))
  pTemp.hypoxT<-summaryBy(pred~Test_Temp, data=pTemp.hypox,FUN=sumfun)
  recTemph[1:3,i]<-pTemp.hypoxT$pred.mean
  #for VRMG
  newdatVRMG.hypox<-expand.grid(Test_Temp=levels(df$Test_Temp),TL=mean(df$TL),Animal_ID=levels(as.factor(as.character(df$Animal_ID))),
                                VentMouthGape_MetricMax=seq(from=min(df$VentMouthGape_MetricMax,na.rm = T),to=max(df$VentMouthGape_MetricMax,na.rm = T),length=100), 
                                midO2=seq(from=min(df$midO2,na.rm = T),to=max(df$midO2,na.rm = T),length=100))
  pVRMG.hypox<-cbind(newdatVRMG.hypox,pred=predict(h9a, newdata = newdatVRMG.hypox,se.fit=TRUE, type="response"))
  pVRMG.hypoxV<-summaryBy(pred~VentMouthGape_MetricMax+midO2, data=pVRMG.hypox,FUN=sumfun)
  recVRMGO2h[1:10000,i]<-pVRMG.hypoxV$pred.mean
  
  print(i)
}

#calculating confidence intervals from bootstrapping
CIrec<-NULL

outVRMGO2h<-CIfunc(recVRMGO2h)
outTemph<-CIfunc(recTemph)
outTLh<-CIfunc(recTLh)


####predicting to get actual model output (non bootstrapping)
#rerun model so mod3a can be used to get actual predictions again (model was slightly different above due to boostrapping)
h9a<-lme(finalMO2~VentMouthGape_MetricMax*midO2+Test_Temp+TL,random=~1|Animal_ID, data=cobia_hypox,na.action = na.omit,
         correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp),method="REML")

#for TL
newdatTL.hypox<-expand.grid(Test_Temp=levels(cobia_hypox$Test_Temp),VentMouthGape_MetricMax=mean(cobia_hypox$VentMouthGape_MetricMax),midO2=mean(cobia_hypox$midO2),
                            TL=seq(from=min(cobia_hypox$TL,na.rm = T),to=max(cobia_hypox$TL,na.rm = T),length=100), 
                            Animal_ID=levels(as.factor(as.character(cobia_hypox$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTL.hypox<-cbind(newdatTL.hypox,pred=predict(h9a, newdata = newdatTL.hypox,se.fit=TRUE, type="response"))
pTL.hypoxTL<-summaryBy(pred~TL, data=pTL.hypox)
pTL.hypoxTL$UpperCI<-pTL.hypoxTL$pred.mean+outTLh #add CI to mean
pTL.hypoxTL$LowerCI<-pTL.hypoxTL$pred.mean-outTLh #subtract CI from mean

#for temp
newdatTemp.hypox<-expand.grid(Test_Temp=levels(cobia_hypox$Test_Temp),VentMouthGape_MetricMax=mean(cobia_hypox$VentMouthGape_MetricMax),midO2=mean(cobia_hypox$midO2),
                              TL=mean(cobia_hypox$TL), 
                              Animal_ID=levels(as.factor(as.character(cobia_hypox$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
pTemp.hypox<-cbind(newdatTemp.hypox,pred=predict(h9a, newdata = newdatTemp.hypox,se.fit=TRUE, type="response"))
pTemp.hypoxT<-summaryBy(pred~Test_Temp, data=pTemp.hypox)
pTemp.hypoxT$UpperCI<-pTemp.hypoxT$pred.mean+outTemph #add CI to mean
pTemp.hypoxT$LowerCI<-pTemp.hypoxT$pred.mean-outTemph #subtract CI from mean

#for VRMG
newdatVRMG.hypox<-expand.grid(Test_Temp=levels(cobia_hypox$Test_Temp),TL=mean(cobia_hypox$TL),Animal_ID=levels(as.factor(as.character(cobia_hypox$Animal_ID))),
                              VentMouthGape_MetricMax=seq(from=min(cobia_hypox$VentMouthGape_MetricMax,na.rm = T),to=max(cobia_hypox$VentMouthGape_MetricMax,na.rm = T),length=100), 
                              midO2=seq(from=min(cobia_hypox$midO2,na.rm = T),to=max(cobia_hypox$midO2,na.rm = T),length=100))
pVRMG.hypox<-cbind(newdatVRMG.hypox,pred=predict(h9a, newdata = newdatVRMG.hypox,se.fit=TRUE, type="response"))
pVRMGO2.hypoxV<-summaryBy(pred~VentMouthGape_MetricMax+midO2, data=pVRMG.hypox)


#####
#Plotting
plot(finalMO2~TL, data=cobia_hypox,pch=21, bg=c("blue","dark grey","red")[cobia_hypox$Test_Temp], ylim=c(min(cobia_hypox$finalMO2),max(cobia_hypox$finalMO2)),
     ylab=expression(paste("Metabolic Rate (mg kg"^"-1","h"^"-1",")")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)
#legend("topleft",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
polygon(x=c(pTL.hypoxTL$TL,rev(pTL.hypoxTL$TL)), y=c(pTL.hypoxTL$LowerCI,rev(pTL.hypoxTL$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~TL,data=pTL.hypoxTL, lwd=2, pch=16)

barplot<-barplot(height=pTemp.hypoxT$pred.mean, names.arg=pTemp.hypoxT$Test_Temp,ylab="Metabolic Rate (mg/kg/h)",xlab="Temperature (°C)", ylim=c(0,200))
arrows(y0=pTemp.hypoxT$LowerCI,y1=pTemp.hypoxT$UpperCI, 
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)


#need to get data in matrix
mat<-matrix(NA,nrow=100, ncol=100)
VRMG_loop<-unique(pVRMGO2.hypoxV$VentMouthGape_MetricMax)
VRMG_loop<-sort(VRMG_loop, decreasing = FALSE)
O2_loops<-unique(pVRMGO2.hypoxV$midO2)
O2_loops<-sort(O2_loops, decreasing = FALSE)
for(i in 1:length(VRMG_loop)){
  o<-which(pVRMGO2.hypoxV$VentMouthGape_MetricMax==VRMG_loop[i])
  o1<-pVRMGO2.hypoxV[o,]
  for(j in 1:length(O2_loops)){
    ox<-which(o1$midO2==O2_loops[j])
    ox1<-o1[ox,]
    mat[i,j]<-ox1$pred.mean
  }
}

matT<-t(mat)

O2_loops1<-O2_loops*32/1000
cobia_hypox$midO21<-cobia_hypox$midO2*32/1000

#plotting (HIGH RES FIG)
library(fields)
#tiff("test1.tiff", width=25, height=15, units="cm",compression="lzw",res=300)
par(mar=c(6,5,2,3)) #265 values comes from taking the 75th quartile of all MO2 vales in normoxia and hypoxia together
image.plot(O2_loops1,VRMG_loop,matT,zlim=c(45, 265), xlab=expression(paste("Oxygen Concentration (mg L"^"-1",")")),ylab="VRMG", axes=T,col=tim.colors(390),cex.axis=1.25,cex.lab=1.25,
           legend.lab=expression(paste("Metabolic Rate (mg kg"^"-1","h"^"-1",")")),legend.line=3,legend.mar=8.5,legend.cex=1.25)
points(VentMouthGape_MetricMax~midO21,data=cobia_hypox,pch=c(21,25)[cobia_hypox$Crash],bg=color.scale(cobia_hypox$finalMO2, col= tim.colors(265),zlim=c(45, 265)))
legend("bottomleft",legend=c("Pre Ccrit","Post Ccrit"),pch=c(25,21),col=c("black","black"))
#dev.off()














#trying to sort out log space and random effects...couldn't really find a solution so just left data in log space for figures
cobia_norm24<-cobia_norm[which(cobia_norm$Test_Temp=="24"),]
cobia_norm28<-cobia_norm[which(cobia_norm$Test_Temp=="28"),]
cobia_norm32<-cobia_norm[which(cobia_norm$Test_Temp=="32"),]

mod3a<-lme(log(finalMO2)~VentMouthGape_Metric,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID), method="REML")
newdatVRMG.norm<-expand.grid(VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100), 
                             Animal_ID=levels(cobia_norm$Animal_ID)) #needed to do that so that cobia that were dead were not included in levels anymore
pVRMG.norm<-cbind(newdatVRMG.norm,pred=predict(mod3a, newdata = newdatVRMG.norm,se.fit=TRUE, type="response"))
pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred+(0.5*disp.lme(mod3a)))
#pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred)
pVRMG.normV<-summaryBy(pred.corrected~VentMouthGape_Metric, data=pVRMG.norm)

plot(finalMO2~VentMouthGape_Metric, data=cobia_norm,pch=21, bg=c("blue","dark grey","red")[cobia_norm$Test_Temp], ylim=c(min(cobia_norm$finalMO2),max(cobia_norm$finalMO2)),
     ylab=expression(paste("Metabolic Rate (mg kg"^"-1","h"^"-1",")")),xlab="VRMG",cex.axis=1.25,cex.lab=1.25)

blah<-lm(finalMO2~VentMouthGape_Metric, data=cobia_norm)
newdatVRMG.norm<-expand.grid(VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blahS<-summaryBy(pred.fit~VentMouthGape_Metric, data=pred.blah)
lines(pred.fit.mean~VentMouthGape_Metric, data=pred.blahS,col="red")

blah<-lm(log(finalMO2)~VentMouthGape_Metric, data=cobia_norm)
newdatVRMG.norm<-expand.grid(VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blah$pred.corrected<-exp(pred.blah$pred.fit+(0.5*dispfun(blah)[3]))
pred.blahS<-summaryBy(pred.corrected~VentMouthGape_Metric, data=pred.blah)
lines(pred.corrected.mean~VentMouthGape_Metric, data=pred.blahS,col="pink")

blah<-lm(finalMO2~VentMouthGape_Metric+Test_Temp, data=cobia_norm)
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blahS<-summaryBy(pred.fit~VentMouthGape_Metric, data=pred.blah)
lines(pred.fit.mean~VentMouthGape_Metric, data=pred.blahS,col="blue")

blah<-lme(finalMO2~VentMouthGape_Metric+Test_Temp,random=~1|Animal_ID, data=cobia_norm)
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Animal_ID=levels(cobia_norm$Animal_ID),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blahS<-summaryBy(pred~VentMouthGape_Metric, data=pred.blah)
lines(pred.mean~VentMouthGape_Metric, data=pred.blahS,col="grey")

blah<-lme(finalMO2~VentMouthGape_Metric+TL+Test_Temp,random=~1|Animal_ID, data=cobia_norm)
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Animal_ID=levels(cobia_norm$Animal_ID),TL=mean(cobia_norm$TL),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blahS<-summaryBy(pred~VentMouthGape_Metric, data=pred.blah)
lines(pred.mean~VentMouthGape_Metric, data=pred.blahS,col="black")

blah<-lme(finalMO2~VentMouthGape_Metric+TL+Time_Diff+Test_Temp,random=~1|Animal_ID, data=cobia_norm)
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Animal_ID=levels(cobia_norm$Animal_ID),TL=mean(cobia_norm$TL),Time_Diff=mean(cobia_norm$Time_Diff),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blahS<-summaryBy(pred~VentMouthGape_Metric, data=pred.blah)
lines(pred.mean~VentMouthGape_Metric, data=pred.blahS,col="orange")

blah<-lme(log(finalMO2)~VentMouthGape_Metric+TL+Test_Temp,random=~1|Animal_ID, data=cobia_norm)
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Animal_ID=levels(cobia_norm$Animal_ID),TL=mean(cobia_norm$TL),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blah$pred.corrected<-exp(pred.blah$pred+(0.5*disp.lme(blah)))
pred.blahS<-summaryBy(pred.corrected~VentMouthGape_Metric, data=pred.blah)
lines(pred.corrected.mean~VentMouthGape_Metric, data=pred.blahS,col="green")

blah<-lme(log(finalMO2)~VentMouthGape_Metric+TL+Test_Temp+before_after_SMR,random=~1|Animal_ID, data=cobia_norm)
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Animal_ID=levels(cobia_norm$Animal_ID),TL=mean(cobia_norm$TL),before_after_SMR=levels(cobia_norm$before_after_SMR),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blah$pred.corrected<-exp(pred.blah$pred+(0.5*disp.lme(blah)))
pred.blahS<-summaryBy(pred.corrected~VentMouthGape_Metric, data=pred.blah)
lines(pred.corrected.mean~VentMouthGape_Metric, data=pred.blahS,col="purple")

blah<-lme(log(finalMO2)~VentMouthGape_Metric+before_after_SMR,random=~1|Animal_ID, data=cobia_norm)
newdatVRMG.norm<-expand.grid(Animal_ID=levels(cobia_norm$Animal_ID),before_after_SMR=levels(cobia_norm$before_after_SMR),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pred.blah<-cbind(newdatVRMG.norm,pred=predict(blah, newdata=newdatVRMG.norm,se.fit=TRUE, type="response"))
pred.blah$pred.corrected<-exp(pred.blah$pred+(0.5*disp.lme(blah)))
pred.blahS<-summaryBy(pred.corrected~VentMouthGape_Metric, data=pred.blah)
lines(pred.corrected.mean~VentMouthGape_Metric, data=pred.blahS,col="dark green")

#need to do this so that we have a time factor for AR1 when using glmmTMB
cobia_norm$Animal_ID_Temp<-as.factor(paste(cobia_norm$Animal_ID, cobia_norm$Test_Temp, sep="_"))
animal_id_temp<-unique(cobia_norm$Animal_ID_Temp)
timefactRec<-NULL
for(i in animal_id_temp){
  a<-which(cobia_norm$Animal_ID_Temp==i)
  a1<-cobia_norm[a,]
  timefact<-as.factor(1:nrow(a1))
  timefactRec<-c(timefactRec,timefact)
}

cobia_norm$TimeFactor<-timefactRec

library(glmmTMB)
mod3a<-glmmTMB(log(finalMO2)~VentMouthGape_Metric+Test_Temp+before_after_SMR+TL+(1|Animal_ID) + ar1(TimeFactor +0 |Animal_ID),data=cobia_norm,na.action = na.omit,
               family=gaussian(link = "identity"))

summary(mod3a)
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),before_after_SMR=levels(cobia_norm$before_after_SMR),TL=mean(cobia_norm$TL),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100), 
                             Animal_ID=levels(cobia_norm$Animal_ID),TimeFactor="1") #needed to do that so that cobia that were dead were not included in levels anymore
pVRMG.norm<-cbind(newdatVRMG.norm,pred=predict(mod3a, newdata = newdatVRMG.norm,se.fit=TRUE, type="response",re.form=NA))
pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred+(0.5*(sqrt(nobs(mod3a)/(1+df.residual(mod3a))))))
pVRMG.normV<-summaryBy(pred.corrected~VentMouthGape_Metric, data=pVRMG.norm)
lines(pred.corrected.mean~VentMouthGape_Metric, data=pVRMG.normV,col="dark green")


sqrt(nobs(mod3a)/(1+df.residual(mod3a)))


mod3a<-gls(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff, data=cobia_norm,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="REML")
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Time_Diff=mean(cobia_norm$Time_Diff),TL=mean(cobia_norm$TL),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100)) #needed to do that so that cobia that were dead were not included in levels anymore
pVRMG.norm<-cbind(newdatVRMG.norm,pred=predict(mod3a, newdata = newdatVRMG.norm,se.fit=TRUE, type="response",re.form=NA))
pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred+(0.5*(sqrt(nobs(mod3a)/(1+df.residual(mod3a))))))
pVRMG.normV<-summaryBy(pred.corrected~VentMouthGape_Metric, data=pVRMG.norm)
lines(pred.corrected.mean~VentMouthGape_Metric, data=pVRMG.normV,col="dark green")


mod3a<-lme(log(finalMO2)~VentMouthGape_Metric+Test_Temp+TL+Time_Diff,random=~1|Animal_ID, data=cobia_norm,na.action = na.omit,
           correlation = corAR1(form=~1|Animal_ID),weights = varIdent(form=~1|Test_Temp), method="REML")
newdatVRMG.norm<-expand.grid(Test_Temp=levels(cobia_norm$Test_Temp),Time_Diff=mean(cobia_norm$Time_Diff),TL=mean(cobia_norm$TL),
                             VentMouthGape_Metric=seq(from=min(cobia_norm$VentMouthGape_Metric,na.rm = T),to=max(cobia_norm$VentMouthGape_Metric,na.rm = T),length=100), 
                             Animal_ID=levels(cobia_norm$Animal_ID))
pVRMG.norm<-cbind(newdatVRMG.norm,pred=predict(mod3a, newdata = newdatVRMG.norm,se.fit=TRUE, type="response"))
pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred+(0.5*disp.lme(mod3a)))
#pVRMG.norm$pred.corrected<-exp(pVRMG.norm$pred)
pVRMG.normV<-summaryBy(pred.corrected~VentMouthGape_Metric, data=pVRMG.norm)
lines(pred.corrected.mean~VentMouthGape_Metric, data=pVRMG.normV,col="black")


