setwd("~/Documents/PhD Project/Respirometry/Respirometry")
sumdat<-read.csv("MO2_methods10_cobia_20180914.csv")

library(StreamMetabolism)

setwd("~/Documents/PhD Project/Respirometry/Respirometry/Cobia Corrected Data")
#read in all csvs ending with "data_corrected.csv" and putting csvs in list
mycsv = dir(pattern="data_corrected.csv")
n <- length(mycsv)
mylist <- vector("list", n)
for(i in 1:n) mylist[[i]] <- read.csv(mycsv[i])


timetoCalmRec<-NULL
timeofCalmRec<-NULL
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
  
  dat$Time_Diff<-dat$Time_Diff+60 #added 60 seconds because it is roughly 60s from end of chase to start of trial (aka time it takes to put top on and start trial on PC), will use time_diff as proxy for time after chase
  
  dat_norm<-dat[which(dat$hypoxia=="no"),]
  
  #MMR, highest MO2 during first 3 hours
  chase_start<-which(dat_norm$Comments=="chase")
  dat_norm$chase<-ifelse(dat_norm$Date_Time<dat_norm$Date_Time[chase_start]+10800 & dat_norm$Date_Time>=dat_norm$Date_Time[chase_start],"yes","no")
  
  MMR<-max(dat_norm$finalMO2[which(dat_norm$chase=="yes")])
  MMRloc<-which.max(dat_norm$finalMO2[which(dat_norm$chase=="yes")])
  
  SMR<-sumdat$SMR[which(as.character(sumdat$Animal_ID)==as.character(dat$Animal_ID[1]) & as.character(sumdat$Test_Temp)==as.character(dat$Test_Temp[1]))]
  
  n<-10
  SMR1<-mean(dat_norm$finalMO2[dat_norm$finalMO2 < quantile(dat_norm$finalMO2,prob=n/100)])
  
  if(all.equal(SMR,SMR1)){
    #SMRsd<-sd(dat_norm$finalMO2[dat_norm$finalMO2 < quantile(dat_norm$finalMO2,prob=n/100)])
    #SMR_high<-SMR1+SMRsd
    #SMR_low<-SMR1-SMRsd
    #time where the median of the all values between 1 SD of SMR
    #firstnearSMRpoint<-min(which(dat_norm$finalMO2<SMR_high & dat_norm$finalMO2>SMR_low))
    #lastnearSMRpoint<-max(which(dat_norm$finalMO2<SMR_high & dat_norm$finalMO2>SMR_low))
    #timetoCalm<-median(dat_norm$Time_Diff[firstnearSMRpoint:lastnearSMRpoint])
    #mintimefrommedian<-min(abs(dat_norm$Time_Diff-timetoCalm)) #if theres an even number of values that taken the median of, then need to make sure time of calm midway between the two points, if odd number of values then this number is 0
    #timeofCalm<-dat_norm$Date_Timen[which(dat_norm$Time_Diff==(timetoCalm-mintimefrommedian))] + mintimefrommedian
    timetoCalm<-dat_norm$Time_Diff[min(which(dat_norm$finalMO2<=SMR1 & dat_norm$Time_Diff>dat_norm$Time_Diff[MMRloc]))]
    timeofCalm<-dat_norm$Date_Timen[min(which(dat_norm$finalMO2<=SMR1 & dat_norm$Time_Diff>dat_norm$Time_Diff[MMRloc]))]
    }else{
    timetoCalm<-NA
    timeofCalm<-NA
  }
  
  timetoCalmRec<-c(timetoCalmRec,timetoCalm)
  timeofCalmRec<-c(timeofCalmRec,timeofCalm)
  print(q)
}
#NAs are caused from one of the cobia that died
sumdat$TimetoCalm<-timetoCalmRec
sumdat$TimeofCalm<-timeofCalmRec

sumdat$TimeofCalmP<-as.POSIXct(sumdat$TimeofCalm, origin="1970-01-01 00:00:00 UTC")

sumdat$Alive_Dead<-ifelse(sumdat$Animal_ID=="C31" & sumdat$Test_Temp=="28", "dead",
                          ifelse(sumdat$Animal_ID=="C24" & sumdat$Test_Temp=="32", "dead",
                                 ifelse(sumdat$Animal_ID=="C32" & sumdat$Test_Temp=="32", "dead",
                                        ifelse(sumdat$Animal_ID=="C28" & sumdat$Test_Temp=="32","dead","alive"))))

#add TL to dataframe
sumdat$TL<-NA
sumdat$TL[which(sumdat$Animal_ID=="C08")]<-95
sumdat$TL[which(sumdat$Animal_ID=="C09")]<-88
sumdat$TL[which(sumdat$Animal_ID=="C11")]<-86
sumdat$TL[which(sumdat$Animal_ID=="C12")]<-105
sumdat$TL[which(sumdat$Animal_ID=="C13")]<-108
sumdat$TL[which(sumdat$Animal_ID=="C15")]<-84
sumdat$TL[which(sumdat$Animal_ID=="C16")]<-87
sumdat$TL[which(sumdat$Animal_ID=="C17")]<-88
sumdat$TL[which(sumdat$Animal_ID=="C18")]<-92
sumdat$TL[which(sumdat$Animal_ID=="C19")]<-91.5
sumdat$TL[which(sumdat$Animal_ID=="C20")]<-87.5
sumdat$TL[which(sumdat$Animal_ID=="C21")]<-82
sumdat$TL[which(sumdat$Animal_ID=="C22")]<-96
sumdat$TL[which(sumdat$Animal_ID=="C23")]<-83
sumdat$TL[which(sumdat$Animal_ID=="C24")]<-81
sumdat$TL[which(sumdat$Animal_ID=="C25")]<-88
sumdat$TL[which(sumdat$Animal_ID=="C26")]<-108
sumdat$TL[which(sumdat$Animal_ID=="C27")]<-88.5
sumdat$TL[which(sumdat$Animal_ID=="C28")]<-102
sumdat$TL[which(sumdat$Animal_ID=="C30")]<-81
sumdat$TL[which(sumdat$Animal_ID=="C31")]<-85
sumdat$TL[which(sumdat$Animal_ID=="C32")]<-81
sumdat$TL[which(sumdat$Animal_ID=="C33")]<-96
sumdat$TL[which(sumdat$Animal_ID=="C34")]<-89
sumdat$TL[which(sumdat$Animal_ID=="C35")]<-102.5

#assign NAs to Ccrit and Scrit for 30_28 and 34_28 because the automation method to calculate ccrit from SMR was not accurate for those two fish
sumdat$Ccrit[which(sumdat$Animal_ID=="C30"&sumdat$Test_Temp=="28")]<-NA
sumdat$Scrit[which(sumdat$Animal_ID=="C30"&sumdat$Test_Temp=="28")]<-NA
sumdat$Ccrit[which(sumdat$Animal_ID=="C34"&sumdat$Test_Temp=="28")]<-NA
sumdat$Scrit[which(sumdat$Animal_ID=="C34"&sumdat$Test_Temp=="28")]<-NA



##########
#Run models
##########

sumdat$Test_Temp<-as.factor(sumdat$Test_Temp)
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
sumdatnoNA<-sumdat[complete.cases(sumdat),]
pairs(cbind(sumdatnoNA$Test_Temp,sumdatnoNA$MMR,sumdatnoNA$SMR,sumdatnoNA$AS,
            sumdatnoNA$Scrit,sumdatnoNA$Ccrit,sumdatnoNA$TL),lower.panel = panel.cor)
#seems like MMR, SMR are correlated with temp which is to be expected, get rid of MMR and SMR
#and ccrit and scrit are correlated, get rid of scrit
pairs(cbind(sumdatnoNA$Test_Temp,sumdatnoNA$AS,
            sumdatnoNA$Ccrit,sumdatnoNA$TL),lower.panel = panel.cor)

sumdatnodead<-sumdat[which(sumdat$Alive_Dead=="alive"),]

moddat<-sumdatnodead

#Ccrit cause problems with NAs and doesn't seem important so leave out of models
#######modeling
library(nlme)
library(HH)
mod0<-lm(TimetoCalm~Test_Temp+AS+TL,data=moddat,na.action = na.omit,x=TRUE)
vif(mod0) #covariates are not collinear, below 5 is okay

#no correlation structure and not modeling variance structure (heterogenity)
mod1<-lme(TimetoCalm~Test_Temp+AS+TL,random=~1|Animal_ID, data=moddat,na.action = na.omit, method = "ML")
BIC(mod0,mod1)#mod0 is better
summary(mod0)
#diagnostics
plot(mod0)
plot(resid(mod0)~moddat$Test_Temp)
plot(resid(mod0)~moddat$TL)
plot(resid(mod0)~moddat$AS)
plot(resid(mod0)~fitted(mod0))

#check for temporal autocorrelation (ie. independence vs non independence)
E<-resid(mod0)
I1<-!is.na(moddat$AS)
Efull<-vector(length=length(moddat$Test_Temp))
Efull<-NA
Efull[I1]<-E
acf(Efull,na.action = na.pass) #definitely independence violated

#no correlation structure, model variance structure (heterogenity) because not a mixed model

#test different correlation structures, nope cause not a mixed model


#different combos of covariates
mod0<-lm(TimetoCalm~Test_Temp+AS+TL,data=moddat,na.action = na.omit)
mod5<-lm(TimetoCalm~Test_Temp+AS,data=moddat,na.action = na.omit)
mod6<-lm(TimetoCalm~Test_Temp+TL, data=moddat,na.action = na.omit)
mod7<-lm(TimetoCalm~TL,data=moddat,na.action = na.omit)
mod8<-lm(TimetoCalm~Test_Temp*AS+TL,data=moddat,na.action = na.omit)
mod9<-lm(TimetoCalm~Test_Temp,data=moddat,na.action = na.omit)
mod10<-lm(TimetoCalm~Test_Temp*AS, data=moddat,na.action = na.omit)

BIC(mod0,mod5,mod6,mod7,mod8,mod9,mod10)
#even though mod10 is barely better than mod5, plots are very interpretable
#go with mod5

#switch to reml
mod10a<-lm(TimetoCalm~Test_Temp*AS, data=moddat,na.action = na.omit)
summary(mod10a)
anova(mod10a)
#temp is sig

library(multcomp)
library(lsmeans)
#look at interaction between AST and Temp
lsmeans(mod10a, pairwise~Test_Temp*AS, adjust="tukey")
#dif between 24 and 28 and 24 and 32
#Test_Temp lsmean   SE df lower.CL upper.CL
#24         63661 3326  4    54427    72896
#28         42220 3758  4    31785    52654
#32         32051 4926  4    18375    45728
slopeAS<-lstrends(mod10a, ~Test_Temp, var="AS")
summary(slopeAS, infer=TRUE)
#MO2 is not sig influenced by AS at 24°
#MO2 sig increases with AS at 28°
#MO2 sig increases with AS at at 32°
summary(pairs(slopeAS), infer=TRUE)
#slope for TBF at 24 is sig different than TBF at 28
#slope for TBF at 24 is sig different than TBF at 32
#slope for TBF at 28 is sig different than TBF at 32


#rsquared calc 
#conditional R2: describes the proportion of variance explained by both the fixed and random factors
#marginal R2: describes the proportion of variance explained by the fixed factor(s) alone
#either package will do
library(piecewiseSEM) #John Lefcheck's package
rsquared(mod10a)
library(gabtool)
r.squared(mod10a)
#conditional R2 is 0.55 which is very good


#from lsmeans (only alive cobia)
barplot<-barplot(height=c(63661/3600,42220/3600,32051/3600), names.arg=c("24","28","32"),ylab="Time to Recovery (h)",xlab="Temperature (°C)",ylim=c(0,22))
arrows(y0=c(54427/3600,31785/3600,18375/3600),y1=c(72896/3600,52654/3600,45728/3600), 
       x0=barplot,x1=barplot,code=3, angle=90,length=.1)


library(doBy)
recAS24<-matrix(nrow=100, ncol=1000, NA)
recAS28<-matrix(nrow=100, ncol=1000, NA)
recAS32<-matrix(nrow=100, ncol=1000, NA)
for(i in 1:1000){
  tryCatch({
  df=moddat[sample(nrow(moddat),nrow(moddat),replace=TRUE),]
  mod10a<-lm(TimetoCalm~Test_Temp*AS, data=df,na.action = na.omit)
  #AS*temp
  newdatAS.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),
                             AS=seq(from=min(df$AS,na.rm = T),to=max(df$AS,na.rm = T),length=100))
  pAS.norm<-cbind(newdatAS.norm,pred=predict(mod10a, newdata = newdatAS.norm, type="response"))
  pAS.normT<-summaryBy(pred~Test_Temp+AS, data=pAS.norm)
  pAS.normT24<-pAS.normT[which(pAS.normT$Test_Temp=="24"),]
  pAS.normT28<-pAS.normT[which(pAS.normT$Test_Temp=="28"),]
  pAS.normT32<-pAS.normT[which(pAS.normT$Test_Temp=="32"),]
  recAS24[1:100,i]<-pAS.normT24$pred.mean
  recAS28[1:100,i]<-pAS.normT28$pred.mean
  recAS32[1:100,i]<-pAS.normT32$pred.mean
    print(i)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\\n")})
}

recAS24<-recAS24[ , apply(recAS24, 2, function(x) !any(is.na(x)))]#gets rid of iterations where model didn't converge
recAS28<-recAS28[ , apply(recAS28, 2, function(x) !any(is.na(x)))]#gets rid of iterations where model didn't converge
recAS32<-recAS32[ , apply(recAS32, 2, function(x) !any(is.na(x)))]#gets rid of iterations where model didn't converge

#calculating confidence intervals from bootstrapping
CIrec<-NULL
CIfunc<-function(data=NA){
  
  for(i in 1:nrow(data)){
    CI<-1.96*sd(data[i,])
    CIrec<-c(CIrec,CI)
  }
  return(CIrec)
}
outAS24<-CIfunc(recAS24)
outAS28<-CIfunc(recAS28)
outAS32<-CIfunc(recAS32)


mod10a<-lm(TimetoCalm~Test_Temp*AS, data=moddat,na.action = na.omit)
newdatAS.norm<-expand.grid(Test_Temp=levels(moddat$Test_Temp),
                            AS=seq(from=min(moddat$AS,na.rm = T),to=max(moddat$AS,na.rm = T),length=100), Animal_ID=levels(as.factor(as.character(moddat$Animal_ID))))
pAS.norm<-cbind(newdatAS.norm,pred=predict(mod10a, newdata = newdatAS.norm, type="response"))
pAS.normT<-summaryBy(pred~Test_Temp+AS, data=pAS.norm)
pAS.normT24<-pAS.normT[which(pAS.normT$Test_Temp=="24"),]
pAS.normT28<-pAS.normT[which(pAS.normT$Test_Temp=="28"),]
pAS.normT32<-pAS.normT[which(pAS.normT$Test_Temp=="32"),]
pAS.normT24$UpperCI<-pAS.normT24$pred.mean+outAS24 #add CI to mean
pAS.normT24$LowerCI<-pAS.normT24$pred.mean-outAS24 #subtract CI from mean
pAS.normT28$UpperCI<-pAS.normT28$pred.mean+outAS28 #add CI to mean
pAS.normT28$LowerCI<-pAS.normT28$pred.mean-outAS28 #subtract CI from mean
pAS.normT32$UpperCI<-pAS.normT32$pred.mean+outAS32 #add CI to mean
pAS.normT32$LowerCI<-pAS.normT32$pred.mean-outAS32 #subtract CI from mean

#before plotting divide to get units in hours
pAS.normT24$pred.mean<-pAS.normT24$pred.mean/3600
pAS.normT24$UpperCI<-pAS.normT24$UpperCI/3600
pAS.normT24$LowerCI<-pAS.normT24$LowerCI/3600
pAS.normT28$pred.mean<-pAS.normT28$pred.mean/3600
pAS.normT28$UpperCI<-pAS.normT28$UpperCI/3600
pAS.normT28$LowerCI<-pAS.normT28$LowerCI/3600
pAS.normT32$pred.mean<-pAS.normT32$pred.mean/3600
pAS.normT32$UpperCI<-pAS.normT32$UpperCI/3600
pAS.normT32$LowerCI<-pAS.normT32$LowerCI/3600

plot(TimetoCalm/3600~AS, data=moddat,pch=21, bg=c("blue","dark grey","red")[moddat$Test_Temp], ylim=c(min(moddat$TimetoCalm/3600),max(moddat$TimetoCalm/3600)),
     ylab="Recovery Time (h)",xlab="AS",cex.axis=1.25,cex.lab=1.25)
#legend("topleft",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
#polygon(x=c(pAS.normT24$AS,rev(pAS.normT24$AS)), y=c(pAS.normT24$LowerCI,rev(pAS.normT24$UpperCI)), col=rgb(0,0,1,alpha=0.4) , border=NA)
lines(pred.mean~AS,data=pAS.normT24, lwd=2, pch=16, col="blue")
#polygon(x=c(pAS.normT28$AS,rev(pAS.normT28$AS)), y=c(pAS.normT28$LowerCI,rev(pAS.normT28$UpperCI)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~AS,data=pAS.normT28, lwd=2, pch=16, col="black")
#polygon(x=c(pAS.normT32$AS,rev(pAS.normT32$AS)), y=c(pAS.normT32$LowerCI,rev(pAS.normT32$UpperCI)), col=rgb(1,0,0,alpha=0.4) , border=NA)
lines(pred.mean~AS,data=pAS.normT32, lwd=2, pch=16, col="red")
#not enough data at 32 to get accurate 95% CI




#Try with SMR in model
#switch to reml
mod12a<-lme(TimetoCalm~SMR,random=~1|Animal_ID, data=moddat,na.action = na.omit, method = "REML")
summary(mod12a)
#SMR is sig

rsquared(mod12a)
sumfun <- function(x, ...){ #this function fixes the problem of when NAs are introduced before summaryBy, get rid of NAs before mean is taken, NAs are made when during bootstrap no points from an ind is selected so would get an NA when trying to predict over that ind
  c(mean=mean(x, na.rm=TRUE, ...), var=var(x, na.rm=TRUE, ...), length=length(x))
}
newdatSMR<-expand.grid(SMR=seq(from=min(moddat$SMR,na.rm = T),to=max(moddat$SMR,na.rm = T),length=100), 
                        Animal_ID=levels(as.factor(as.character(moddat$Animal_ID)))) #needed to do that so that cobia that were dead were not included in levels anymore
newdatSMR<-cbind(newdatSMR,pred=predict(mod12a, newdata = newdatSMR,se.fit=TRUE, type="response"))
newdatSMRsum<-summaryBy(pred~SMR, data=newdatSMR,sumfun=TRUE)

par(mar=(c(5,4.5,4,1)))
plot(TimetoCalm~SMR, data=moddat,pch=21, bg=c("blue","dark grey","red")[moddat$Test_Temp],
     ylab="Time to Calm",xlab="SMR",cex.axis=1.25,cex.lab=1.25)
legend("topright",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"))
lines(pred.mean~SMR,data=newdatSMRsum, lwd=2, pch=16)















######how to deal with interaction between AS and temp, not a good model though
library(doBy)
recAST24<-matrix(nrow=100, ncol=5000, NA)
recAST28<-matrix(nrow=100, ncol=5000, NA)
recAST32<-matrix(nrow=100, ncol=5000, NA)
for(i in 1:5000){
  tryCatch({
    df=moddat[sample(nrow(moddat),nrow(moddat),replace=TRUE),]
    mod10a<-lme(TimetoCalm~Test_Temp*AS,random=~1|Animal_ID, data=df,na.action = na.omit, method = "REML")
    #AS*temp
    newdatAST.norm<-expand.grid(Test_Temp=levels(df$Test_Temp),
                                AS=seq(from=min(df$AS,na.rm = T),to=max(df$AS,na.rm = T),length=100), Animal_ID=levels(as.factor(as.character(df$Animal_ID))))
    pAST.norm<-cbind(newdatAST.norm,pred=predict(mod10a, newdata = newdatAST.norm,se.fit=TRUE, type="response"))
    pAST.normT<-summaryBy(pred~Test_Temp+AS, data=pAST.norm)
    pAST.norm24<-pAST.normT[which(pAST.normT$Test_Temp=="24"),]
    pAST.norm28<-pAST.normT[which(pAST.normT$Test_Temp=="28"),]
    pAST.norm32<-pAST.normT[which(pAST.normT$Test_Temp=="32"),]
    recAST24[1:100,i]<-pAST.norm24$pred.mean
    recAST28[1:100,i]<-pAST.norm28$pred.mean
    recAST32[1:100,i]<-pAST.norm32$pred.mean
    
    print(i)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\\n")})
}

recAST24noNA<-recAST24[ , apply(recAST24, 2, function(x) !any(is.na(x)))]
recAST28noNA<-recAST28[ , apply(recAST28, 2, function(x) !any(is.na(x)))]
recAST32noNA<-recAST32[ , apply(recAST32, 2, function(x) !any(is.na(x)))]


#calculating confidence intervals from bootstrapping
CIrec<-NULL
CIfunc<-function(data=NA){
  
  for(i in 1:nrow(data)){
    CI<-1.96*sd(data[i,])
    CIrec<-c(CIrec,CI)
  }
  return(CIrec)
}

outAS24<-CIfunc(recAST24noNA)
outAS28<-CIfunc(recAST28noNA)
outAS32<-CIfunc(recAST32noNA)

#for AS and temp
mod10a<-lme(TimetoCalm~Test_Temp*AS,random=~1|Animal_ID, data=moddat,na.action = na.omit, method = "REML")
newdatAST.norm<-expand.grid(Test_Temp=levels(moddat$Test_Temp),
                            AS=seq(from=min(moddat$AS,na.rm = T),to=max(moddat$AS,na.rm = T),length=100), Animal_ID=levels(as.factor(as.character(moddat$Animal_ID))))
pAST.norm<-cbind(newdatAST.norm,pred=predict(mod10a, newdata = newdatAST.norm,se.fit=TRUE, type="response"))
pAST.normT<-summaryBy(pred~Test_Temp+AS, data=pAST.norm)
pAST.norm24<-pAST.normT[which(pAST.normT$Test_Temp=="24"),]
pAST.norm28<-pAST.normT[which(pAST.normT$Test_Temp=="28"),]
pAST.norm32<-pAST.normT[which(pAST.normT$Test_Temp=="32"),]
pAST.norm24$UpperCI<-pAST.norm24$pred.mean+outAS24 #add CI to mean
pAST.norm24$LowerCI<-pAST.norm24$pred.mean-outAS24 #subtract CI from mean
pAST.norm28$UpperCI<-pAST.norm28$pred.mean+outAS28 #add CI to mean
pAST.norm28$LowerCI<-pAST.norm28$pred.mean-outAS28 #subtract CI from mean
pAST.norm32$UpperCI<-pAST.norm32$pred.mean+outAS32 #add CI to mean
pAST.norm32$LowerCI<-pAST.norm32$pred.mean-outAS32 #subtract CI from mean

plot(TimetoCalm~AS, data=moddat,pch=21, bg=c("blue","dark grey","red")[moddat$Test_Temp], ylim=c(min(moddat$TimetoCalm),max(moddat$TimetoCalm)),
     ylab="",xlab="AS",cex.axis=1.25,cex.lab=1.25)
legend("bottomright",legend=c("24°C","28°C","32°C"),fill=c("blue","dark grey","red"),cex=1.5)
polygon(x=c(pAST.norm24$AS,rev(pAST.norm24$AS)), y=c(pAST.norm24$Lower,rev(pAST.norm24$Upper)), col=rgb(0,0,1,alpha=0.4) , border=NA)
lines(pred.mean~AS,data=pAST.norm24, lwd=2, pch=16,col="blue")
polygon(x=c(pAST.norm28$AS,rev(pAST.norm28$AS)), y=c(pAST.norm28$Lower,rev(pAST.norm28$Upper)), col=rgb(0,0,0,alpha=0.4) , border=NA)
lines(pred.mean~AS,data=pAST.norm28, lwd=2, pch=16)
polygon(x=c(pAST.norm32$AS,rev(pAST.norm32$AS)), y=c(pAST.norm32$Lower,rev(pAST.norm32$Upper)), col=rgb(1,0,0,alpha=0.4) , border=NA)
lines(pred.mean~AS,data=pAST.norm32, lwd=2, pch=16, col="red")
text(1.6,480,"C",cex=2)


