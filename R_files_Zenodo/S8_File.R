#1# dDTR across latitude and elevation 

#Fig. 1 b-d
library(raster)
library(s2dverification)

create_DTR<-function(my.raster){
  find_DTR<-function(x){max(x,na.rm=T)-min(x,na.rm=T)}
  out_raster<-calc(my.raster,find_DTR)
  return(out_raster)
} # function to generate a raster file of DTR from a raster file of daily temperature profile

T7_shaded<-brick("TA1cm_soil_100_7.nc")
## raster available from microclim database
DTR7_shaded<-create_DTR(T7_shaded) # create raster for DTR in shaded environment (forest)

T7_exposed<-brick("TA1cm_soil_0_7.nc")
## raster available from microclim database
DTR7_exposed<-create_DTR(T7_exposed) # create raster for DTR in exposed environment (deforested habitat)

brks <- seq(-36, 36, by=4)
nb <- length(brks)-1 
my.col <- s2dverification::clim.colors(nb)

plot(DTR7_exposed-DTR7_shaded,breaks=brks,col=my.col,las=1,ylim=c(0,66.5),xlim=c(-180,180),xaxt="n",yaxt="n")
axis(1,las=1,yaxp=c(-180,180,6))
axis(2,las=1,yaxp=c(0,70,7)) #Fig. 1b

elev<-raster("wc2.1_30s_elev.tif") #read global elevation data 
## raster available from WorldClim 2.1
Elevation_res<-resample(elev, DTR7_shaded, method="bilinear")
xy<-xyFromCell(Elevation_res,1:length(Elevation_res))

library(msir)
data_dDTR7<-data.frame(DTR_exposed=values(DTR7_exposed),DTR_shaded=values(DTR7_shaded),Elevation=values(Elevation_res),Lon=xy[,1],Lat=xy[,2])
colnames(data_dDTR7)[1:2]<-c("DTR_exposed","DTR7_shaded")
data_dDTR7<-data_dDTR7[complete.cases(data_dDTR7), ]
data_dDTR7$dDTR<-data_dDTR7$DTR_exposed-data_dDTR7$DTR7_shaded

#Temperate #Fig. 1c
#Elevation was log transformed before loess regression since samples from higher elevation are relatively scarce, 
#which makes itdifficult to dierctly estimate the trend and confidience interval at above XXX m. a. s. l. 
data_dDTR7_M<-data_dDTR7[(data_dDTR7$Lat)>=23.5&(data_dDTR7$Lat)<66.5,]
md_dDTR_M<-loess.sd(y=data_dDTR7_M$dDTR[data_dDTR7_M$Elevation>0],x=log(data_dDTR7_M$Elevation[data_dDTR7_M$Elevation>0]+1), nsigma = 1.96)
plot(x=data_dDTR7_M$Elevation,y=data_dDTR7_M$dDTR,xlab="Elevation",ylab="dDTR",xlim=c(0,6000),ylim=c(0,30),pch=".",col="lightblue",las=1)
subX<-c(seq(1,0.9*length(md_dDTR_M$x),100),seq(ceiling(0.9*length(md_dDTR_M$x)),length(md_dDTR_M$x)))
polygon(x=exp(c(md_dDTR_M$x[subX],rev(md_dDTR_M$x[subX])))-1,y=c(md_dDTR_M$upper[subX],rev(md_dDTR_M$lower[subX])),col=rgb(0.1,0,0,0.2),border=NA)
lines(x=exp(md_dDTR_M$x)-1,y=md_dDTR_M$y,col="darkblue",lwd=2)

#Tropics #Fig. 1d
data_dDTR7_L<-data_dDTR7[(data_dDTR7$Lat)>=0&(data_dDTR7$Lat)<23.5,]
md_dDTR_L<-loess.sd(y=data_dDTR7_L$dDTR[data_dDTR7_L$Elevation>0],x=data_dDTR7_L$Elevation[data_dDTR7_L$Elevation>0], nsigma = 1.96)
plot(x=data_dDTR7_L$Elevation,y=data_dDTR7_L$dDTR,xlab="Elevation",ylab="dDTR",xlim=c(0,6000),ylim=c(0,30),pch=".",col="lightblue",las=1)
polygon(x=c(md_dDTR_L$x,rev(md_dDTR_L$x)),y=c(md_dDTR_L$upper,rev(md_dDTR_L$lower)),col=rgb(0.1,0,0,0.2),border=NA)
lines(x=md_dDTR_L$x,y=md_dDTR_L$y,col="darkblue",lwd=2)


#2# Field temperature monitoring 
#Fig. 2 c-d; Table S1
library(lme4)
library(lmerTest)
data1<-read.csv("Data S1.csv")
data1$Elev_std<-(data1$Elevation-mean(data1$Elevation,na.rm=T))/(2*sd(data1$Elevation,na.rm=T))
model_Tmean<-lmer(Tmean~Elev_std*Type+(1|Plot),data=data1)   # test for table S1A, Fig 2c
model_DTR<-lmer(DTR~Elev_std*Type+(1|Plot),data=data1)       # test for table S1B, Fig 2d
model_Tmax<-lmer(Tmax~Elev_std*Type+(1|Plot),data=data1)     # test for table S1C
model_Tmin<-lmer(Tmin~Elev_std*Type+(1|Plot),data=data1)     # test for table S1D


#3# Field experiments
library(lme4)
library(lmerTest)
library(ordinal)
data2_1<-read.csv("Data S2.csv")
data2_1$Elev_std<-(data2_1$Elev-mean(data2_1$Elev,na.rm=T))/(2*sd(data2_1$Elev,na.rm=T))
data2_1$Tmean_std<-(data2_1$Tmean-mean(data2_1$Tmean,na.rm=T))/(2*sd(data2_1$Tmean,na.rm=T))
data2_1$DTR_std<-(data2_1$DTR-mean(data2_1$DTR,na.rm=T))/(2*sd(data2_1$DTR,na.rm=T))
data2_1$JD_std<-(data2_1$JulianDate-mean(data2_1$JulianDate,na.rm=T))/(2*sd(data2_1$JulianDate,na.rm=T))

model_canopy<-lmer(Canopy_Cover~Type+(1|Plot_Elev/Plot),data=data2_1) # test for Fig. 2b
model_TS2<-glmer(Bury~I(Elev_std^2)+Elev_std*Type+JD_std+(1|Plot_Elev/Plot),family=binomial,data=data2_1,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) # test for Table S2, Fig 2e
model_TS3<-glmer(Bury~I(Tmean_std^2)+Tmean_std*DTR_std+JD_std+(1|Plot_Elev/Plot),family=binomial,data=data2_1) # test for table S3, Fig 3

data2_2<-data2_1[is.na(data2_1$Day_Presence)==F,]
data2_2$Elev_std<-(data2_2$Elev-mean(data2_2$Elev,na.rm=T))/(2*sd(data2_2$Elev,na.rm=T))
data2_2$Tmean_std<-(data2_2$Tmean-mean(data2_2$Tmean,na.rm=T))/(2*sd(data2_2$Tmean,na.rm=T))
data2_2$DTR_std<-(data2_2$DTR-mean(data2_2$DTR,na.rm=T))/(2*sd(data2_2$DTR,na.rm=T))
data2_2$JD_std<-(data2_2$JulianDate-mean(data2_2$JulianDate,na.rm=T))/(2*sd(data2_2$JulianDate,na.rm=T))

model_TS4a<-clmm(ordered(Day_Presence)~Tmean_std+DTR_std+JD_std+(1|Plot),data=data2_2) #test for Fig.4 a & b; Table S4a

data2_3<-data2_1[is.na(data2_1$Day_Bury)==F,]
data2_3$Elev_std<-(data2_3$Elev-mean(data2_3$Elev,na.rm=T))/(2*sd(data2_3$Elev,na.rm=T))
data2_3$Tmean_std<-(data2_3$Tmean-mean(data2_3$Tmean,na.rm=T))/(2*sd(data2_3$Tmean,na.rm=T))
data2_3$DTR_std<-(data2_3$DTR-mean(data2_3$DTR,na.rm=T))/(2*sd(data2_3$DTR,na.rm=T))
data2_3$Tmax_std<-(data2_3$Tmax-mean(data2_3$Tmax,na.rm=T))/(2*sd(data2_3$Tmax,na.rm=T))
data2_3$Tmin_std<-(data2_3$Tmin-mean(data2_3$Tmin,na.rm=T))/(2*sd(data2_3$Tmin,na.rm=T))
data2_3$JD_std<-(data2_3$JulianDate-mean(data2_3$JulianDate,na.rm=T))/(2*sd(data2_3$JulianDate,na.rm=T))

model_TS4b<-clmm(ordered(Day_Bury)~Tmean_std+DTR_std+JD_std+(1|Plot),data=data2_3,nAGQ=13) #test for Fig.4 c & d; Table S4b
model_TS5<-clmm(ordered(Day_Bury)~Tmax_std+Tmin_std+JD_std+(1|Plot),data=data2_3,nAGQ=20) #test for Fig. S2; Table S5


#4# beetle arrival time 

#(results from re-analysing video taken in Chan et al (2019))
#Fig. S1
beetle_arrival_time<-read.csv("Data S3.csv")
beetle_arrival_time$Present_time<-strptime(as.character(beetle_arrival_time$Present_time), format="%H")
hist(beetle_arrival_time$Present_time,breaks="hour",freq=T,las=1,col="gray",border=NA,ylim=c(0,150),xlab="Time of day",main=NA)


#5# Lab exp. for maggots 
#Fig. S3

data.maggot.day<-read.csv("Data S4.csv")
wilcox.test(day~DTR,data.maggot.day)
boxplot(day~DTR,data=data.maggot.day,col=c("dodgerblue4","firebrick1"),las=1,ylim=c(2.5,8.5)) #test for Fig. S3a

data_maggot_length<-read.csv("Data S5.csv")
mmg<-lmer(LENGTH~TYPE+(1|SN),data=data_maggot_length) #test for Fig. S3b


#6# TPC for beetles 
#Fig S4a
beetle_TPC_data<-read.csv("Data S6.csv")
beetle_TPC_data$b_rate<-1/(beetle_TPC_data$buring.day)
beetle_TPC_data$b_rate[is.na(beetle_TPC_data$b_rate)==T]<-0

test_temperature<-data.frame(temperature=seq(8,22,0.1))
beetle_performance<-predict(loess(asin(sqrt(b_rate))~temperature,data=beetle_TPC_data),test_temperature$temperature)

topt_beetle=test_temperature$temperature[which(beetle_performance==max(beetle_performance))]
TPC_beetle<-nls(b_rate ~ ifelse(topt_beetle >= temperature, exp(-((temperature-topt_beetle)/2/sigma2)^2), 1-((temperature-topt_beetle)/(topt_beetle-tmax))^2),
                 start = list(tmax = 23, sigma2=4),
                 data = beetle_TPC_data)
summary(TPC_beetle) # Table S6


#7# TPC for maggots 
#Fig S4b
library(Metrics)
f_pinguis_1<-function(x){sin((-0.366+0.8539*x-0.1521*x^2)/180*pi)^2}# hatching rate
f_pinguis_2<-function(x){sin((-0.4480+0.5977*x-0.0943*x^2)/180*pi)^2}# maggot survival
f_pinguis_max<-optimize(function(x){(f_pinguis_1(x)*f_pinguis_2(x))},lower=15,upper=38,maximum=T)$objective
f_pinguis<-function(x){sqrt(f_pinguis_1(x)*f_pinguis_2(x)/f_pinguis_max)} # performance function of maggot based on Yang & Shiao (2014)
Topt_maggot<-optimize(f_pinguis,lower=15,upper=38,maximum=T)$maximum # Topt of maggot

TPC<-function(t,topt,tmax,sigma){
  if(t<=topt){p=exp(-((t-topt)/(2*sigma))^2)}else{
    p=1-((t-topt)/(topt-tmax))^2
  }
  if(p<0){p<-0}
  return(p)
} # general form of TPC

find_rmse<-function(tmax,sigma,topt=Topt_maggot){
  actual<-unlist(lapply(seq(15,37.3,0.1),f_pinguis)) #expectation from function f_pinguis
  TPC_maggot<-function(t){TPC(topt=topt,tmax=tmax,sigma=sigma,t=t)} #TPC based on given tmax and sigma
  predicted<-unlist(lapply(seq(15,37.3,0.1),TPC_maggot)) #expectation from function TPC_maggot
  return(rmse(actual,predicted)) # find rmse for two sets of expectations
} 
# function to calculate rmse between a TPC with given tmax and sigma, 
# and the expectation of the performance function of maggot based on Yang & Shiao (2014)

test_rmse<-data.frame(expand.grid(tmax=seq(20,40,0.01),sigma=seq(0.1,7,0.01))) # all tested combinations of tmax and sigma

test_rmse$rmse<-NA
for(i in 1:nrow(test_rmse)){
  test_rmse$rmse[i]<-find_rmse(tmax=test_rmse$tmax[i],sigma=test_rmse$sigma[i])
} # rmse for each combination of tmax and sigma 

test_rmse[test_rmse$rmse==min(test_rmse$rmse),] #best combination of tmax and sigma


#8# modeling realized TPC #Fig. 5 d-g & Fig. S5

Tprofile<-read.csv("Data S7.csv") #Daily temperature profiles obtained from the fields

#Data re-arangement
CT<-colnames(Tprofile)[6:53]
data.rearangement<-function(CT){
  dd<-Tprofile[,c("Plot","Type","ID","Elevation","Tmean","DTR")]
  dd$Time<-unlist(strsplit(CT,"X"))[2]
  dd$temperature<-Tprofile[,c(CT)]
  return(dd)
}
Tprofile<-do.call(rbind,lapply(CT,data.rearangement))
Tprofile$Time2<-factor(Tprofile$Time)
Tprofile$Elev_std<-(Tprofile$Elevation-mean(Tprofile$Elevation,na.rm=T))/(2*sd(Tprofile$Elevation,na.rm=T))
Tprofile$Tmean_std<-(Tprofile$Tmean-mean(Tprofile$Tmean,na.rm=T))/(2*sd(Tprofile$Tmean,na.rm=T))
Tprofile$DTR_std<-(Tprofile$DTR-mean(Tprofile$DTR,na.rm=T))/(2*sd(Tprofile$DTR,na.rm=T))

model_1<-lmer(temperature~Elev_std*Type*Time2+(1|Plot),data=Tprofile) 
#GLMM model for predicting temperature for any given combination of elevation,habitat type, and time 

model_2<-lmer(temperature~Tmean_std*DTR_std*Time2+(1|Plot),data=Tprofile)
#GLMM model for predicting temperature for any given combination of Tmean, DTR, and time 

library(emmeans)
emm_options(lmerTest.limit = 3300)
emm_options(pbkrtest.limit = 3300)

# Fig. S5
find_realized_TPC_1<-function(Elevation,Type,day=c("d","n","b"),bw=1.5,topt,tmax,sigma){
  require(kdensity)
  Elev_std<-(Elevation-mean(Tprofile$Elevation,na.rm=T))/(2*sd(Tprofile$Elevation,na.rm=T))
  dd<-data.frame(emmeans::emmeans(model_1,~Time2|Elev_std|Type,at=list(Elev_std=Elev_std,Type=Type)))
  dd$Time2<-as.numeric(as.character(dd$Time2))
  dd<-dd[order(dd$Time2),]
  dd$Day<-0
  dd$Day[12:38]<-1
  if(day=="d"){kde<-kdensity(dd$emmean[12:38],bw=bw)}
  if(day=="n"){kde<-kdensity(dd$emmean[c(1:11,39:48)],bw=bw)}
  if(day=="b"){kde<-kdensity(dd$emmean,bw=bw)}
  f1<-function(t){exp(-((t-topt)/(2*sigma))^2)*kde(t)}
  f2<-function(t){(1-((t-topt)/(topt-tmax))^2)*kde(t)}
  p1<-integrate(f1,lower=0,upper=topt)$value
  p2<-integrate(f2,lower=topt,upper=tmax)$value
  p=p1+p2
  if(day=="d"){return(27*p/48)}
  if(day=="n"){return(21*p/48)}
  if(day=="b"){return(p)}
}

test_elev=seq(1000,3000,1)
fitness_beetle_A<-unlist(lapply(test_elev,function(x){find_realized_TPC_1(x,"A",day="n",topt=16.7,tmax=23.05,sigma=2.04)}))
fitness_beetle_N<-unlist(lapply(test_elev,function(x){find_realized_TPC_1(x,"N",day="n",topt=16.7,tmax=23.05,sigma=2.04)}))

fitness_maggot_A<-unlist(lapply(test_elev,function(x){find_realized_TPC_1(x,"A",day="b",topt=29.45,tmax=37.09,sigma=4.58)}))
fitness_maggot_N<-unlist(lapply(test_elev,function(x){find_realized_TPC_1(x,"N",day="b",topt=29.45,tmax=37.09,sigma=4.58)}))

test_elev[which(fitness_beetle_N/fitness_maggot_N==max(fitness_beetle_N/fitness_maggot_N))]
test_elev[which(fitness_beetle_A/fitness_maggot_A==max(fitness_beetle_A/fitness_maggot_A))]

plot(x=test_elev,y=(fitness_beetle_A/fitness_maggot_A),type="l",col="orange", xlim=c(1000,3000),ylim=c(0,3),las=1,xlab="Elevation",ylab="W_beetle/W_maggot",lwd=2,xaxt="n")
lines(x=test_elev,y=(fitness_beetle_N/fitness_maggot_N),type="l",col="darkblue",lwd=2)
axis(1,at=seq(1000,3000,400),labels=seq(1000,3000,400))  

#Fig. 5 d-g
find_realized_TPC_2<-function(Tmean,DTR,day=c("d","n","b"),bw=1.5,topt,tmax,sigma){
  require(kdensity)
  Tmean_std<-(Tmean-mean(Tprofile$Tmean,na.rm=T))/(2*sd(Tprofile$Tmean,na.rm=T))
  DTR_std<-(DTR-mean(Tprofile$DTR,na.rm=T))/(2*sd(Tprofile$DTR,na.rm=T))
  dd<-data.frame(emmeans::emmeans(model_2,~Time2|Tmean_std|DTR_std,at=list(Tmean_std=Tmean_std,DTR_std=DTR_std)))
  dd$Tmean<-Tmean
  dd$DTR<-DTR
  dd$Time2<-as.numeric(as.character(dd$Time2))
  dd<-dd[order(dd$Time2),]
  dd$Day<-0
  dd$Day[12:38]<-1
  if(day=="d"){kde<-kdensity(dd$emmean[12:38],bw=bw)}
  if(day=="n"){kde<-kdensity(dd$emmean[c(1:11,39:48)],bw=bw)}
  if(day=="b"){kde<-kdensity(dd$emmean,bw=bw)}
  f1<-function(t){exp(-((t-topt)/(2*sigma))^2)*kde(t)}
  f2<-function(t){(1-((t-topt)/(topt-tmax))^2)*kde(t)}
  p1<-integrate(f1,lower=0,upper=topt)$value
  p2<-integrate(f2,lower=topt,upper=tmax)$value
  p=p1+p2
  if(day=="d"){return(27*p/48)}
  if(day=="n"){return(21*p/48)}
  if(day=="b"){return(p)}
}

test_Tmean=seq(12,21,0.1)
test_DTR=c(7,14)

fitness_beetle_7<-unlist(lapply(test_Tmean,function(x){find_realized_TPC_2(x,DTR=7,day="n",topt=16.7,tmax=23.05,sigma=2.04)}))
fitness_beetle_14<-unlist(lapply(test_Tmean,function(x){find_realized_TPC_2(x,DTR=14,day="n",topt=16.7,tmax=23.05,sigma=2.04)}))

fitness_maggot_7<-unlist(lapply(test_Tmean,function(x){find_realized_TPC_2(x,DTR=7,day="b",topt=29.45,tmax=37.09,sigma=4.58)}))
fitness_maggot_14<-unlist(lapply(test_Tmean,function(x){find_realized_TPC_2(x,DTR=14,day="b",topt=29.45,tmax=37.09,sigma=4.58)}))

#Fig. 5d
plot(x=test_Tmean,y=fitness_beetle_7,type="l",col="red", xlim=c(12,21),ylim=c(0,0.5),las=1,xlab="Mean daily temperature",ylab="W_beetle",lwd=2,xaxt="n")
lines(x=test_Tmean,y=fitness_maggot_7,type="l",col="black",lwd=2)
axis(1,at=seq(12,21,3),labels=seq(12,21,3)) 

#Fig. 5e
plot(x=test_Tmean,y=fitness_beetle_14,type="l",col="red", xlim=c(12,21),ylim=c(0,0.5),las=1,xlab="Mean daily temperature",ylab="W_beetle",lwd=2,xaxt="n")
lines(x=test_Tmean,y=fitness_maggot_14,type="l",col="black",lwd=2)
axis(1,at=seq(12,21,3),labels=seq(12,21,3))

#Fig. 5f
plot(x=test_Tmean,y=fitness_beetle_7/fitness_maggot_7,type="l",col="darkblue", xlim=c(12,21),ylim=c(0,2.5),las=1,xlab="Mean daily temperature",ylab="W_beetle/W_maggot",lwd=2,xaxt="n")
axis(1,at=seq(12,21,3),labels=seq(12,21,3))

#Fig. 5g
plot(x=test_Tmean,y=fitness_beetle_14/fitness_maggot_14,type="l",col="orange", xlim=c(12,21),ylim=c(0,2.5),las=1,xlab="Mean daily temperature",ylab="W_beetle/W_maggot",lwd=2,xaxt="n")
axis(1,at=seq(12,21,3),labels=seq(12,21,3))

