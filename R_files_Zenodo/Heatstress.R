#####-----Heat stress of Takydromus septentrionalis-----#####
#####-----Main calculation-----#####
#####-----Aug 27th 2020-----#####
####----load packages----####
library(data.table)
library(raster)
library(mgcv)
library(lmerTest)

####----import present soil temp data----####
T_2016<-fread("China_soilT5cm_summer_2016.csv") #from 2016 April to August
T_2016<-T_2016[,c(2:10)]
T_2016$SoilT5cm[which(T_2016$SoilT5cm>100)]<-NA
NA_2016<-T_2016[is.na(T_2016$SoilT5cm),]
SiteNA16<-unique(NA_2016$ID)
T_2016<-T_2016[!T_2016$ID %in% SiteNA16,] #remove all sites which contains NA


####----import change in mean air temp----####
#(Import current climate data and CMIP5 future climate data respectively)
data_current<-raster("current-bio1.bil") #get current(1960-1990, WorldClim 1.4,Generic grid format bio10min) annual mean temperature 
T_current<-as.data.frame(rasterToPoints(data_current))
colnames(T_current)<-c("Longitude","Latitude","T_current")
data_2050_rcp26<-raster("bc26bi501.tif")
data_2050_rcp45<-raster("bc45bi501.tif")
data_2050_rcp60<-raster("bc60bi501.tif")
data_2050_rcp85<-raster("bc85bi501.tif")
data_2070_rcp26<-raster("bc26bi701.tif")
data_2070_rcp45<-raster("bc45bi701.tif")
data_2070_rcp60<-raster("bc60bi701.tif")
data_2070_rcp85<-raster("bc85bi701.tif")
T_2050_rcp26<-as.data.frame(rasterToPoints(data_2050_rcp26))#get future(2050, RCP26, WorldClim 1.4,BCC-CSM1-1) annual mean temperature 
T_2050_rcp45<-as.data.frame(rasterToPoints(data_2050_rcp45))
T_2050_rcp60<-as.data.frame(rasterToPoints(data_2050_rcp60))
T_2050_rcp85<-as.data.frame(rasterToPoints(data_2050_rcp85))
T_2070_rcp26<-as.data.frame(rasterToPoints(data_2070_rcp26))#get future(2070, RCP26, WorldClim 1.4,BCC-CSM1-1) annual mean temperature 
T_2070_rcp45<-as.data.frame(rasterToPoints(data_2070_rcp45))
T_2070_rcp60<-as.data.frame(rasterToPoints(data_2070_rcp60))
T_2070_rcp85<-as.data.frame(rasterToPoints(data_2070_rcp85))
colnames(T_2070_rcp26)<-c("Longitude","Latitude","T_2070_rcp26")
colnames(T_2070_rcp45)<-c("Longitude","Latitude","T_2070_rcp45")
colnames(T_2070_rcp60)<-c("Longitude","Latitude","T_2070_rcp60")
colnames(T_2070_rcp85)<-c("Longitude","Latitude","T_2070_rcp85")
colnames(T_2050_rcp26)<-c("Longitude","Latitude","T_2050_rcp26")
colnames(T_2050_rcp45)<-c("Longitude","Latitude","T_2050_rcp45")
colnames(T_2050_rcp60)<-c("Longitude","Latitude","T_2050_rcp60")
colnames(T_2050_rcp85)<-c("Longitude","Latitude","T_2050_rcp85")
T_change<-cbind(T_current,T_2050_rcp26[,3],T_2050_rcp45[,3],T_2050_rcp60[,3],T_2050_rcp85[,3],T_2070_rcp26[,3],T_2070_rcp45[,3],T_2070_rcp60[,3],T_2070_rcp85[,3])
colnames(T_change)<-c("Longitude","Latitude","T_current","T_2050_rcp26","T_2050_rcp45","T_2050_rcp60","T_2050_rcp85","T_2070_rcp26","T_2070_rcp45","T_2070_rcp60","T_2070_rcp85")
T_change$Change_2050_rcp26<-T_change$T_2050_rcp26-T_change$T_current
T_change$Change_2050_rcp45<-T_change$T_2050_rcp45-T_change$T_current
T_change$Change_2050_rcp60<-T_change$T_2050_rcp60-T_change$T_current
T_change$Change_2050_rcp85<-T_change$T_2050_rcp85-T_change$T_current
T_change$Change_2070_rcp26<-T_change$T_2070_rcp26-T_change$T_current
T_change$Change_2070_rcp45<-T_change$T_2070_rcp45-T_change$T_current
T_change$Change_2070_rcp60<-T_change$T_2070_rcp60-T_change$T_current
T_change$Change_2070_rcp85<-T_change$T_2070_rcp85-T_change$T_current
T_change[,3:19]<-T_change[,3:19]/10
T_change<-T_change[T_change$Longitude>75,]
T_change<-T_change[T_change$Longitude<135,]
T_change<-T_change[T_change$Latitude>18,]
T_change<-T_change[T_change$Latitude<54,]  ##remove sites outside China



####----calculate soil temperature in 2050 and 2070----####
#creat a reference dataframe to relate coords and ID
coor_refer<-T_2016[,c(1:3)]
coor_refer<-coor_refer[!duplicated(coor_refer$ID),]
#Match T_change with coor_refer by lon and lat
coor_refer_change<-matrix(NA,nrow=0,ncol=8) #blank matrix for storing nearest T_change points for coor_refer
colnames(coor_refer_change)<-c("2050_rcp26","2050_rcp45","2050_rcp60","2050_rcp85","2070_rcp26","2070_rcp45","2070_rcp60","2070_rcp85")
for (i in 1:nrow(coor_refer)) {
  distant<-function(lon,lat){
    (lon-coor_refer$Longitude[i])^2+(lat-coor_refer$Latitude[i])^2
  }
  index<-which.min(mapply(distant,T_change$Longitude,T_change$Latitude))
  coor_refer_change<-rbind(coor_refer_change,T_change[index,12:19])
  # progress bar
  all <- nrow(coor_refer)
  # create progress bar
  pb <- txtProgressBar(min = 0, max = all, style = 3)
  # update progress bar
  setTxtProgressBar(pb, i)  
}
coor_refer<-cbind(coor_refer,coor_refer_change)
T_future<-merge(T_2016,coor_refer,by="ID")
T_future$SoilT5cm_2050_rcp26<-T_future$SoilT5cm+T_future$Change_2050_rcp26
T_future$SoilT5cm_2050_rcp45<-T_future$SoilT5cm+T_future$Change_2050_rcp45
T_future$SoilT5cm_2050_rcp60<-T_future$SoilT5cm+T_future$Change_2050_rcp60
T_future$SoilT5cm_2050_rcp85<-T_future$SoilT5cm+T_future$Change_2050_rcp85
T_future$SoilT5cm_2070_rcp26<-T_future$SoilT5cm+T_future$Change_2070_rcp26
T_future$SoilT5cm_2070_rcp45<-T_future$SoilT5cm+T_future$Change_2070_rcp45
T_future$SoilT5cm_2070_rcp60<-T_future$SoilT5cm+T_future$Change_2070_rcp60
T_future$SoilT5cm_2070_rcp85<-T_future$SoilT5cm+T_future$Change_2070_rcp85
T_future<-T_future[,c(1:8,20:27)]
colnames(T_future)<-c(colnames(T_2016)[1:8],"SoilT5cm_2050_rcp26","SoilT5cm_2050_rcp45","SoilT5cm_2050_rcp60","SoilT5cm_2050_rcp85","SoilT5cm_2070_rcp26","SoilT5cm_2070_rcp45","SoilT5cm_2070_rcp60","SoilT5cm_2070_rcp85")
write.csv(T_future,"T_future.csv",row.names=F)


####---fitting functions for CTmax and developmental rate----####
T_future<-fread("T_future.csv")
CTmax_YT<-42.16
CTmax_NJ<-41.30
CTmax_ND<-39.93
###fitting functions for CTmax
data_CTmax1<-read.csv("EAHT.csv")

data_fit1<-lmer(data=data_CTmax1,EAHT~Tincubation*Location+(1|Mother))
lmerControl(optimizer = "Nelder_Mead")
anova(data_fit1)
shapiro.test(residuals(data_fit1)) #w=0.99303, p=0.7453

data_CTmax2<-data_CTmax1
data_CTmax2$Location<-sub("NJ","aNJ",data_CTmax2$Location)
data_fit2<-lmer(data=data_CTmax2,EAHT~Tincubation*Location+(1|Mother))
summary(data_fit2)
anova(data_fit2)

data_CTmax3<-data_CTmax1
data_CTmax3$Location<-sub("YT","aYT",data_CTmax1$Location)
data_fit3<-lmer(data=data_CTmax3,EAHT~Tincubation*Location+(1|Mother))
summary(data_fit3)
anova(data_fit3)


data_NJ<-data_CTmax[data_CTmax$Location=="NJ",]
data_YT<-data_CTmax[data_CTmax$Location=="YT",]
data_ND<-data_CTmax[data_CTmax$Location=="ND",]

NJ_fit1<-lm(data=data_NJ,EAHT~Tincubation) #R2=0.00271 AIC=159.6604
NJ_fit2<-lm(data=data_NJ,EAHT~poly(Tincubation,2,raw=TRUE)) #Adjusted R2=0.2959 AIC=145.6277 !!!
summary(NJ_fit2)
NJ_fit3<-lm(data=data_NJ,EAHT~poly(Tincubation,3,raw=TRUE)) #Adjusted R2=0.2838 AIC=147.2739
NJ_fit4<-lm(data=data_NJ,EAHT~poly(Tincubation,4,raw=TRUE)) #Adjusted R2=0.2782 AIC=148.489
xx<-seq(22,34,length=50)
plot(x=data_NJ$Tincubation,y=data_NJ$EAHT,pch=19,ylim=c(35,47),xlim=c(23,33),xlab="Incubation temp (C)",ylab="CTmax of embryos (C)",col="green")
#lines(xx,predict(NJ_fit1,data.frame(Tincubation=xx)),col="red")
lines(xx,predict(NJ_fit2,data.frame(Tincubation=xx)),col="green")
#lines(xx,predict(NJ_fit3,data.frame(Tincubation=xx)),col="blue")
#lines(xx,predict(NJ_fit4,data.frame(Tincubation=xx)),col="purple")
YT_fit1<-lm(data=data_YT,EAHT~Tincubation) #Adjusted R2=-0.02251 AIC=127.1953 !!!
summary(YT_fit1)
YT_fit2<-lm(data=data_YT,EAHT~poly(Tincubation,2,raw=TRUE)) #Adjusted R2=-0.03281 AIC=128.5874 
YT_fit3<-lm(data=data_YT,EAHT~poly(Tincubation,3,raw=TRUE)) #Adjusted R2=-0.0575 AIC=130.5663
YT_fit4<-lm(data=data_YT,EAHT~poly(Tincubation,4,raw=TRUE)) #Adjusted R2=-0.0828 AIC=132.5187
xx<-seq(22,34,length=50)
#plot(x=data_YT$Tincubation,y=data_YT$EAHT,pch=19,ylim=c(32,50),main="Yantai",xlab="Incubation temp (C)",ylab="CTmax of embryos (C)")
points(x=data_YT$Tincubation,y=data_YT$EAHT,pch=19,col="blue")
lines(xx,predict(YT_fit1,data.frame(Tincubation=xx)),col="blue")
#lines(xx,predict(YT_fit2,data.frame(Tincubation=xx)),col="green")
#lines(xx,predict(YT_fit3,data.frame(Tincubation=xx)),col="blue")
#lines(xx,predict(YT_fit4,data.frame(Tincubation=xx)),col="purple")
ND_fit1<-lm(data=data_ND,EAHT~Tincubation) #Adjusted R2=0.3109 AIC=195.793
ND_fit2<-lm(data=data_ND,EAHT~poly(Tincubation,2,raw=TRUE)) #Adjusted R2=0.3257 AIC=195.6959 !!!
summary(ND_fit2)
ND_fit3<-lm(data=data_ND,EAHT~poly(Tincubation,3,raw=TRUE)) #Adjusted R2=0.3239 AIC=196.7443
ND_fit4<-lm(data=data_ND,EAHT~poly(Tincubation,4,raw=TRUE)) #Adjusted R2=0.3081 AIC=198.7443
xx<-seq(22,34,length=50)
#plot(x=data_ND$Tincubation,y=data_ND$EAHT,pch=19,ylim=c(32,50),main="Ningde",xlab="Incubation temp (C)",ylab="CTmax of embryos (C)")
points(x=data_ND$Tincubation,y=data_ND$EAHT,pch=19,col="red")
#lines(xx,predict(ND_fit1,data.frame(Tincubation=xx)),col="red")
lines(xx,predict(ND_fit2,data.frame(Tincubation=xx)),col="red")
#lines(xx,predict(ND_fit3,data.frame(Tincubation=xx)),col="blue")
#lines(xx,predict(ND_fit4,data.frame(Tincubation=xx)),col="purple")
#put all the points and lines in one figure
text("Yantai",x=30.9,y=46.5,type="bold",col="blue")
text("Nanjing",x=31,y=45.8,type="bold",col="green")
text("Ningde",x=30.99,y=45.1,type="bold",col="red")

###creat functions for CTmax
##--Yantai--##
FunCTm_YT<-function(T_mean){
  if (T_mean>=20 & T_mean<=36) CTm<-0.008634*T_mean+41.917876
  if (T_mean<20) CTm<-0.008634*20+41.917876
  if (T_mean>36) CTm<-0.008634*36+41.917876
  return(CTm)
}

##--Nanjing--##
FunCTm_NJ<-function(T_mean){
  if (T_mean>=20 & T_mean<=36) CTm<--0.1152*T_mean^2+6.3565*T_mean-45.3656
  if (T_mean<20) CTm<--0.1152*20^2+6.3565*20-45.3656
  if (T_mean>36) CTm<--0.1152*36^2+6.3565*36-45.3656
  return(CTm)
}

##--Ningde--##
FunCTm_ND<-function(T_mean){
  if (T_mean>=20 & T_mean<=36) CTm<--0.05227*T_mean^2+2.51287*T_mean+10.98166
  if (T_mean<20) CTm<--0.05227*20^2+2.51287*20+10.98166
  if (T_mean>36) CTm<--0.05227*36^2+2.51287*36+10.98166
  return(CTm)
}

###fitting function for developmental rate
dat_incuP<-read.csv("Incubation period.csv")
dat_incuP$D<-1/dat_incuP$Incubation_period
mod_YT<-lm(dat_incuP[dat_incuP$Pop=="Yantai",]$D~dat_incuP[dat_incuP$Pop=="Yantai",]$Tincubation)
summary(mod_YT)
mod_NJ<-lm(dat_incuP[dat_incuP$Pop=="Nanjing",]$D~dat_incuP[dat_incuP$Pop=="Nanjing",]$Tincubation)
summary(mod_NJ)
mod_ND<-lm(dat_incuP[dat_incuP$Pop=="Ningde",]$D~dat_incuP[dat_incuP$Pop=="Ningde",]$Tincubation)
summary(mod_ND)

###create functions for developmental rate
FunD_YT<-function(T){
  if (T<=10.56772) D=0
  else D=1/24*(0.0018930*T-0.0200047)
  return(D)
}

FunD_NJ<-function(T){
  if (T<=11.94933) D=0
  else D=1/24*((1.934e-03)*T-2.311e-02)
  return(D)
}

FunD_ND<-function(T){
  if (T<=11.12239) D=0
  else D=1/24*((1.871e-03)*T-2.081e-02)
  return(D)
}

####----Calculate heat stress----####
###Calculate developmental rate for current climate
T_2016$Julian_d<-do.call(paste, list(T_2016$Month, T_2016$Day, T_2016$Year))
T_2016$Julian_d<-as.Date(T_2016$Julian_d, format=c("%m %d %Y"))
T_2016$Julian_d<-as.numeric(format(T_2016$Julian_d, "%j"))
T_2016$Julian_h<-T_2016$Julian_d*24+T_2016$Hour
T_2016<-as.data.frame(T_2016)

T_2016$D_YT<-lapply(T_2016[,9],FunD_YT)
T_2016$D_NJ<-lapply(T_2016[,9],FunD_NJ)
T_2016$D_ND<-lapply(T_2016[,9],FunD_ND)

###Calculate maximun temperature for each sites
T16_sub<-T_2016[,c(1,9,11,12:14)] #reduced size of the matrix for "forloop"
###Calculate developmental rate for 2070 climate
T_future$Julian_d <- do.call(paste, list(T_future$Month, T_future$Day, T_future$Year))
T_future$Julian_d <- as.Date(T_future$Julian_d, format=c("%m %d %Y"))
T_future$Julian_d <- as.numeric(format(T_future$Julian_d, "%j"))
T_future$Julian_h <- T_future$Julian_d*24+T_future$Hour
T_future<-as.data.frame(T_future)
T_future$D_2050_YT<-lapply(T_future[,11],FunD_YT)
T_future$D_2050_NJ<-lapply(T_future[,11],FunD_NJ)
T_future$D_2050_ND<-lapply(T_future[,11],FunD_ND)
T_future$D_2070_YT<-lapply(T_future[,15],FunD_YT)
T_future$D_2070_NJ<-lapply(T_future[,15],FunD_NJ)
T_future$D_2070_ND<-lapply(T_future[,15],FunD_ND)
T_future$Julian_d <- do.call(paste, list(T_2016$Month, T_2016$Day, T_2016$Year))
T_future_sub<-T_future[,c(1,9:16,18:24)] #reduced size of the matrix for "forloop"

###-Calculate heat stress
Site<-unique(T16_sub$ID)
result16<-as.data.frame(Site)
result16$Stress16_YT<-NA
result16$Stress16_NJ<-NA
result16$Stress16_ND<-NA
result16$Stress16_YT_noplas<-NA
result16$Stress16_NJ_noplas<-NA
result16$Stress16_ND_noplas<-NA
result16$EUTT16_YT<-NA #mean EUTT in reproductive season
result16$EUTT16_NJ<-NA
result16$EUTT16_ND<-NA
result50<-as.data.frame(Site)
result50$Stress50_YT<-NA
result50$Stress50_NJ<-NA
result50$Stress50_ND<-NA
result50$Stress50_YT_noplas<-NA
result50$Stress50_NJ_noplas<-NA
result50$Stress50_ND_noplas<-NA
result50$EUTT50_YT<-NA #mean EUTT in reproductive season
result50$EUTT50_NJ<-NA
result50$EUTT50_ND<-NA
result70<-as.data.frame(Site)
result70$Stress70_YT<-NA
result70$Stress70_NJ<-NA
result70$Stress70_ND<-NA
result70$Stress70_YT_noplas<-NA
result70$Stress70_NJ_noplas<-NA
result70$Stress70_ND_noplas<-NA
result70$EUTT70_YT<-NA #mean EUTT in reproductive season
result70$EUTT70_NJ<-NA
result70$EUTT70_ND<-NA
change50<-as.data.frame(Site)
change50$change_YT<-NA
change50$change_NJ<-NA
change50$change_ND<-NA
change50$change_YT_noplas<-NA
change50$change_NJ_noplas<-NA
change50$change_ND_noplas<-NA
change70<-as.data.frame(Site)
change70$change_YT<-NA
change70$change_NJ<-NA
change70$change_ND<-NA
change70$change_YT_noplas<-NA
change70$change_NJ_noplas<-NA
change70$change_ND_noplas<-NA

for (i in 1:length(Site)) {
  dat16_sub<-subset(T16_sub,T16_sub$ID==Site[i]) #extract temp data for each site
  dat16_sub<-dat16_sub[order(dat16_sub$Julian_h),]
  Stress16_YT<-list()
  Stress16_NJ<-list()
  Stress16_ND<-list()
  Stress16_YT_noplas<-list()
  Stress16_NJ_noplas<-list()
  Stress16_ND_noplas<-list()
  EUTT16_YT<-list()
  EUTT16_NJ<-list()
  EUTT16_ND<-list()
  dat_future_sub<-subset(T_future_sub,T_future_sub$ID==Site[i]) #extract temp data for each site
  dat_future_sub<-dat_future_sub[order(dat_future_sub$Julian_h),]
  Stress50_YT<-list()
  Stress50_NJ<-list()
  Stress50_ND<-list()
  Stress50_YT_noplas<-list()
  Stress50_NJ_noplas<-list()
  Stress50_ND_noplas<-list()
  EUTT50_YT<-list()
  EUTT50_NJ<-list()
  EUTT50_ND<-list()
  Stress70_YT<-list()
  Stress70_NJ<-list()
  Stress70_ND<-list()
  Stress70_YT_noplas<-list()
  Stress70_NJ_noplas<-list()
  Stress70_ND_noplas<-list()
  EUTT70_YT<-list()
  EUTT70_NJ<-list()
  EUTT70_ND<-list()
  for (p in 721:nrow(dat16_sub)) {
    Dev16_list_YT<-dat16_sub$D_YT[(p-1):(p-720)] #a list of accomplished development in each hour, for 30days before certain hour
    Dev16_list_NJ<-dat16_sub$D_NJ[(p-1):(p-720)]
    Dev16_list_ND<-dat16_sub$D_ND[(p-1):(p-720)]
    accumDev16_list_YT<-cumsum(Dev16_list_YT)
    accumDev16_list_NJ<-cumsum(Dev16_list_NJ)
    accumDev16_list_ND<-cumsum(Dev16_list_ND)
    Half16_incuP_YT<-findInterval(0.5,accumDev16_list_YT) #count how many hours to develop to 50%
    Half16_incuP_NJ<-findInterval(0.5,accumDev16_list_NJ) 
    Half16_incuP_ND<-findInterval(0.5,accumDev16_list_ND)
    T16_mean_YT<-mean(dat16_sub$SoilT5cm[(p-1):(p-Half16_incuP_YT)])
    T16_mean_NJ<-mean(dat16_sub$SoilT5cm[(p-1):(p-Half16_incuP_NJ)])
    T16_mean_ND<-mean(dat16_sub$SoilT5cm[(p-1):(p-Half16_incuP_ND)])
    EUTT16_YT[p-720]<-FunCTm_YT(T16_mean_YT) #CTmax function of Yantai population
    EUTT16_NJ[p-720]<-FunCTm_NJ(T16_mean_NJ) #CTmax function of Nanjing population
    EUTT16_ND[p-720]<-FunCTm_ND(T16_mean_ND) #CTmax function of Ningde population
    
    Dev50_list_YT<-dat_future_sub$D_2050_YT[(p-1):(p-720)] #a list of accomplished development in each hour, for 30days before certain hour
    Dev50_list_NJ<-dat_future_sub$D_2050_NJ[(p-1):(p-720)]
    Dev50_list_ND<-dat_future_sub$D_2050_ND[(p-1):(p-720)]
    accumDev50_list_YT<-cumsum(Dev50_list_YT)
    accumDev50_list_NJ<-cumsum(Dev50_list_NJ)
    accumDev50_list_ND<-cumsum(Dev50_list_ND)
    Half50_incuP_YT<-findInterval(0.5,accumDev50_list_YT) #count how many hours to develop to 50%
    Half50_incuP_NJ<-findInterval(0.5,accumDev50_list_NJ) 
    Half50_incuP_ND<-findInterval(0.5,accumDev50_list_ND)
    T50_mean_YT<-mean(dat_future_sub$SoilT5cm_2050_rcp60[(p-1):(p-Half50_incuP_YT)])
    T50_mean_NJ<-mean(dat_future_sub$SoilT5cm_2050_rcp60[(p-1):(p-Half50_incuP_NJ)])
    T50_mean_ND<-mean(dat_future_sub$SoilT5cm_2050_rcp60[(p-1):(p-Half50_incuP_ND)])
    EUTT50_YT[p-720]<-FunCTm_YT(T50_mean_YT) #CTmax function of Yantai population
    EUTT50_NJ[p-720]<-FunCTm_NJ(T50_mean_NJ) #CTmax function of Nanjing population
    EUTT50_ND[p-720]<-FunCTm_ND(T50_mean_ND) #CTmax function of Ningde population
    
    Dev70_list_YT<-dat_future_sub$D_2070_YT[(p-1):(p-720)] #a list of accomplished development in each hour, for 30days before certain hour
    Dev70_list_NJ<-dat_future_sub$D_2070_NJ[(p-1):(p-720)]
    Dev70_list_ND<-dat_future_sub$D_2070_ND[(p-1):(p-720)]
    accumDev70_list_YT<-cumsum(Dev70_list_YT)
    accumDev70_list_NJ<-cumsum(Dev70_list_NJ)
    accumDev70_list_ND<-cumsum(Dev70_list_ND)
    Half70_incuP_YT<-findInterval(0.5,accumDev70_list_YT) #count how many hours to develop to 50%
    Half70_incuP_NJ<-findInterval(0.5,accumDev70_list_NJ) 
    Half70_incuP_ND<-findInterval(0.5,accumDev70_list_ND)
    T70_mean_YT<-mean(dat_future_sub$SoilT5cm_2070_rcp60[(p-1):(p-Half70_incuP_YT)])
    T70_mean_NJ<-mean(dat_future_sub$SoilT5cm_2070_rcp60[(p-1):(p-Half70_incuP_NJ)])
    T70_mean_ND<-mean(dat_future_sub$SoilT5cm_2070_rcp60[(p-1):(p-Half70_incuP_ND)])
    EUTT70_YT[p-720]<-FunCTm_YT(T70_mean_YT) #CTmax function of Yantai population
    EUTT70_NJ[p-720]<-FunCTm_NJ(T70_mean_NJ) #CTmax function of Nanjing population
    EUTT70_ND[p-720]<-FunCTm_ND(T70_mean_ND) #CTmax function of Ningde population
    
    if (dat16_sub$SoilT5cm[p]>EUTT16_YT[p-720] && Half16_incuP_YT!=720) Stress16_YT[length(Stress16_YT)+1]=1 else Stress16_YT[length(Stress16_YT)+1]=0
    if (dat16_sub$SoilT5cm[p]>EUTT16_NJ[p-720] && Half16_incuP_NJ!=720) Stress16_NJ[length(Stress16_NJ)+1]=1 else Stress16_NJ[length(Stress16_NJ)+1]=0
    if (dat16_sub$SoilT5cm[p]>EUTT16_ND[p-720] && Half16_incuP_ND!=720) Stress16_ND[length(Stress16_ND)+1]=1 else Stress16_ND[length(Stress16_ND)+1]=0
    if (dat_future_sub$SoilT5cm_2050_rcp60[p]>EUTT50_YT[p-720] && Half50_incuP_YT!=720) Stress50_YT[length(Stress50_YT)+1]=1 else Stress50_YT[length(Stress50_YT)+1]=0
    if (dat_future_sub$SoilT5cm_2050_rcp60[p]>EUTT50_NJ[p-720] && Half50_incuP_NJ!=720) Stress50_NJ[length(Stress50_NJ)+1]=1 else Stress50_NJ[length(Stress50_NJ)+1]=0
    if (dat_future_sub$SoilT5cm_2050_rcp60[p]>EUTT50_ND[p-720] && Half50_incuP_ND!=720) Stress50_ND[length(Stress50_ND)+1]=1 else Stress50_ND[length(Stress50_ND)+1]=0
    if (dat_future_sub$SoilT5cm_2070_rcp60[p]>EUTT70_YT[p-720] && Half70_incuP_YT!=720) Stress70_YT[length(Stress70_YT)+1]=1 else Stress70_YT[length(Stress70_YT)+1]=0
    if (dat_future_sub$SoilT5cm_2070_rcp60[p]>EUTT70_NJ[p-720] && Half70_incuP_NJ!=720) Stress70_NJ[length(Stress70_NJ)+1]=1 else Stress70_NJ[length(Stress70_NJ)+1]=0
    if (dat_future_sub$SoilT5cm_2070_rcp60[p]>EUTT70_ND[p-720] && Half70_incuP_ND!=720) Stress70_ND[length(Stress70_ND)+1]=1 else Stress70_ND[length(Stress70_ND)+1]=0
    if (dat16_sub$SoilT5cm[p]>CTmax_YT && Half16_incuP_YT!=720) Stress16_YT_noplas[length(Stress16_YT_noplas)+1]=1 else Stress16_YT_noplas[length(Stress16_YT_noplas)+1]=0
    if (dat16_sub$SoilT5cm[p]>CTmax_NJ && Half16_incuP_NJ!=720) Stress16_NJ_noplas[length(Stress16_NJ_noplas)+1]=1 else Stress16_NJ_noplas[length(Stress16_NJ_noplas)+1]=0
    if (dat16_sub$SoilT5cm[p]>CTmax_ND && Half16_incuP_ND!=720) Stress16_ND_noplas[length(Stress16_ND_noplas)+1]=1 else Stress16_ND_noplas[length(Stress16_ND_noplas)+1]=0
    if (dat_future_sub$SoilT5cm_2050_rcp60[p]>CTmax_YT && Half50_incuP_YT!=720) Stress50_YT_noplas[length(Stress50_YT_noplas)+1]=1 else Stress50_YT_noplas[length(Stress50_YT_noplas)+1]=0
    if (dat_future_sub$SoilT5cm_2050_rcp60[p]>CTmax_NJ && Half50_incuP_NJ!=720) Stress50_NJ_noplas[length(Stress50_NJ_noplas)+1]=1 else Stress50_NJ_noplas[length(Stress50_NJ_noplas)+1]=0
    if (dat_future_sub$SoilT5cm_2050_rcp60[p]>CTmax_ND && Half50_incuP_ND!=720) Stress50_ND_noplas[length(Stress50_ND_noplas)+1]=1 else Stress50_ND_noplas[length(Stress50_ND_noplas)+1]=0
    if (dat_future_sub$SoilT5cm_2070_rcp60[p]>CTmax_YT && Half70_incuP_YT!=720) Stress70_YT_noplas[length(Stress70_YT_noplas)+1]=1 else Stress70_YT_noplas[length(Stress70_YT_noplas)+1]=0
    if (dat_future_sub$SoilT5cm_2070_rcp60[p]>CTmax_NJ && Half70_incuP_NJ!=720) Stress70_NJ_noplas[length(Stress70_NJ_noplas)+1]=1 else Stress70_NJ_noplas[length(Stress70_NJ_noplas)+1]=0
    if (dat_future_sub$SoilT5cm_2070_rcp60[p]>CTmax_ND && Half70_incuP_ND!=720) Stress70_ND_noplas[length(Stress70_ND_noplas)+1]=1 else Stress70_ND_noplas[length(Stress70_ND_noplas)+1]=0

  }
  
  cal.HSday<-function(x){
    group<-rep(1:floor(length(x)/24),each=24)
    daily_HS<-tapply(x[1:length(group)], group, FUN = sum)
    return(length(which(daily_HS>0)))
  }
  
  result16$Stress16_YT[i]=cal.HSday(unlist(Stress16_YT))
  result16$Stress16_NJ[i]=cal.HSday(unlist(Stress16_NJ))
  result16$Stress16_ND[i]=cal.HSday(unlist(Stress16_ND))
  
  result16$Stress16_YT_noplas[i]=cal.HSday(unlist(Stress16_YT_noplas))
  result16$Stress16_NJ_noplas[i]=cal.HSday(unlist(Stress16_NJ_noplas))
  result16$Stress16_ND_noplas[i]=cal.HSday(unlist(Stress16_ND_noplas))
  
  result16$EUTT16_YT[i]=do.call(mean,EUTT16_YT) #mean EUTT in reproductive season
  result16$EUTT16_NJ[i]=do.call(mean,EUTT16_NJ)
  result16$EUTT16_ND[i]=do.call(mean,EUTT16_ND)
  
  result50$Stress50_YT[i]=cal.HSday(unlist(Stress50_YT))
  result50$Stress50_NJ[i]=cal.HSday(unlist(Stress50_NJ))
  result50$Stress50_ND[i]=cal.HSday(unlist(Stress50_ND))
  
  result50$Stress50_YT_noplas[i]=cal.HSday(unlist(Stress50_YT_noplas))
  result50$Stress50_NJ_noplas[i]=cal.HSday(unlist(Stress50_NJ_noplas))
  result50$Stress50_ND_noplas[i]=cal.HSday(unlist(Stress50_ND_noplas))
  
  result50$EUTT50_YT[i]=do.call(mean,EUTT50_YT) #mean EUTT in reproductive season
  result50$EUTT50_NJ[i]=do.call(mean,EUTT50_NJ)
  result50$EUTT50_ND[i]=do.call(mean,EUTT50_ND)
  
  result70$Stress70_YT[i]=cal.HSday(unlist(Stress70_YT))
  result70$Stress70_NJ[i]=cal.HSday(unlist(Stress70_NJ))
  result70$Stress70_ND[i]=cal.HSday(unlist(Stress70_ND))
  
  result70$Stress70_YT_noplas[i]=cal.HSday(unlist(Stress70_YT_noplas))
  result70$Stress70_NJ_noplas[i]=cal.HSday(unlist(Stress70_NJ_noplas))
  result70$Stress70_ND_noplas[i]=cal.HSday(unlist(Stress70_ND_noplas))
  
  result70$EUTT70_YT[i]=do.call(mean,EUTT70_YT) #mean EUTT in reproductive season
  result70$EUTT70_NJ[i]=do.call(mean,EUTT70_NJ)
  result70$EUTT70_ND[i]=do.call(mean,EUTT70_ND)
  
  change50$change_YT[i]=result50$Stress50_YT[i]-result16$Stress16_YT[i]
  change50$change_NJ[i]=result50$Stress50_NJ[i]-result16$Stress16_NJ[i]
  change50$change_ND[i]=result50$Stress50_ND[i]-result16$Stress16_ND[i]
  
  change50$change_YT_noplas[i]=result50$Stress50_YT_noplas[i]-result16$Stress16_YT_noplas[i]
  change50$change_NJ_noplas[i]=result50$Stress50_NJ_noplas[i]-result16$Stress16_NJ_noplas[i]
  change50$change_ND_noplas[i]=result50$Stress50_ND_noplas[i]-result16$Stress16_ND_noplas[i]
  
  change70$change_YT[i]=result70$Stress70_YT[i]-result16$Stress16_YT[i]
  change70$change_NJ[i]=result70$Stress70_NJ[i]-result16$Stress16_NJ[i]
  change70$change_ND[i]=result70$Stress70_ND[i]-result16$Stress16_ND[i]
  
  change70$change_YT_noplas[i]=result70$Stress70_YT_noplas[i]-result16$Stress16_YT_noplas[i]
  change70$change_NJ_noplas[i]=result70$Stress70_NJ_noplas[i]-result16$Stress16_NJ_noplas[i]
  change70$change_ND_noplas[i]=result70$Stress70_ND_noplas[i]-result16$Stress16_ND_noplas[i]
  
  # progress bar
  total <- length(Site)
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  # update progress bar
  setTxtProgressBar(pb, i)
}

#write.csv(change,"change_stress.csv")
result16<-merge(result16,coor_refer[,1:3],by.x="Site",by.y="ID",all.x=F,all.y=F)
result50<-merge(result50,coor_refer[,1:3],by.x="Site",by.y="ID",all.x=F,all.y=F)
result70<-merge(result70,coor_refer[,1:3],by.x="Site",by.y="ID",all.x=F,all.y=F)
change50<-merge(change50,coor_refer[,1:3],by.x="Site",by.y="ID",all.x=F,all.y=F)
change70<-merge(change70,coor_refer[,1:3],by.x="Site",by.y="ID",all.x=F,all.y=F)
daymean_2016<-aggregate(data=T_2016,SoilT5cm~ID+Year+Month+Day,FUN=mean)
daymean_2050<-aggregate(data=T_future,SoilT5cm_2050_rcp60~ID+Year+Month+Day,FUN=mean)
daymean_2070<-aggregate(data=T_future,SoilT5cm_2070_rcp60~ID+Year+Month+Day,FUN=mean)
daymean_seasonmean2016<-aggregate(data=T_2016,SoilT5cm~ID,FUN=mean)
daymean_seasonmean2050<-aggregate(data=T_future,SoilT5cm_2050_rcp60~ID,FUN=mean)
daymean_seasonmean2070<-aggregate(data=T_future,SoilT5cm_2070_rcp60~ID,FUN=mean)
daymax_2016<-aggregate(data=T_2016,SoilT5cm~ID+Year+Month+Day,FUN=max)
daymax_2050<-aggregate(data=T_future,SoilT5cm_2050_rcp60~ID+Year+Month+Day,FUN=max)
daymax_2070<-aggregate(data=T_future,SoilT5cm_2070_rcp60~ID+Year+Month+Day,FUN=max)
daymax_seasonmean2016<-aggregate(data=daymax_2016,SoilT5cm~ID,FUN=mean)
daymax_seasonmean2050<-aggregate(data=daymax_2050,SoilT5cm_2050_rcp60~ID,FUN=mean)
daymax_seasonmean2070<-aggregate(data=daymax_2070,SoilT5cm_2070_rcp60~ID,FUN=mean)
daymax_seasonmax2016<-aggregate(data=daymax_2016,SoilT5cm~ID,FUN=max)
daymax_seasonmax2050<-aggregate(data=daymax_2050,SoilT5cm_2050_rcp60~ID,FUN=max)
daymax_seasonmax2070<-aggregate(data=daymax_2070,SoilT5cm_2070_rcp60~ID,FUN=max)
season2016<-cbind(daymean_seasonmean2016,daymax_seasonmean2016[2],daymax_seasonmax2016[,2])
season2050<-cbind(daymean_seasonmean2050,daymax_seasonmean2050[2],daymax_seasonmax2050[,2])
season2070<-cbind(daymean_seasonmean2070,daymax_seasonmean2070[2],daymax_seasonmax2070[,2])
#mean daily maximum
result16<-merge(result16,season2016,by.x="Site",by.y="ID",all.x=F,all.y=F)
result50<-merge(result50,season2050,by.x="Site",by.y="ID",all.x=F,all.y=F)
result70<-merge(result70,season2070,by.x="Site",by.y="ID",all.x=F,all.y=F)
names(result16)[13]<-"Meandailymean"
names(result50)[13]<-"Meandailymean"
names(result70)[13]<-"Meandailymean"
names(result16)[14]<-"Meandailymax"
names(result50)[14]<-"Meandailymax"
names(result70)[14]<-"Meandailymax"
names(result16)[15]<-"Maxdailymax"
names(result50)[15]<-"Maxdailymax"
names(result70)[15]<-"Maxdailymax"
write.csv(result16,"result16.csv",row.names=F)
write.csv(result50,"result50.csv",row.names=F)
write.csv(result70,"result70.csv",row.names=F)
write.csv(change50,"change50.csv",row.names=F)
write.csv(change70,"change70.csv",row.names=F)

