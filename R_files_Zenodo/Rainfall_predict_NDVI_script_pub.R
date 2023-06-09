#Rainfall predict NDVI script

#Using this script to determine if rainfall predicts NDVI

library(here)

#Load data
rainall = read.csv(here::here("Output files","rainallyears.csv")) #Rain has slightly fewer rows bc only goes through Nov 2018
rainall$Date = as.Date(rainall$Date,"%m/%d/%y")

L8ndvi = read.csv(here::here("Output files","Landsat 8 NDVI yearly polygons noclouds.csv"))
L8ndvi$Date = as.Date(L8ndvi$Date,"%m/%d/%y")

#Get rainfall in cm
rainall$amountcm = rainall$amount/10



####Does rainfall predict NDVI during the field season?####

#Get field season NDVI
L8ndvi.16 = L8ndvi[which(L8ndvi$Date>="2015-06-04" & L8ndvi$Date<="2015-09-13"),]
L8ndvi.17 = L8ndvi[which(L8ndvi$Date>="2016-06-04" & L8ndvi$Date<="2016-09-13"),]
L8ndvi.18 = L8ndvi[which(L8ndvi$Date>="2017-06-04" & L8ndvi$Date<="2017-09-13"),]
L8ndvi.19 = L8ndvi[which(L8ndvi$Date>="2018-06-04" & L8ndvi$Date<="2018-09-13"),]

L8ndvi.fs = rbind(L8ndvi.16,L8ndvi.17,L8ndvi.18,L8ndvi.19)

#Check distribution
library(car)
library(MASS)
hist(L8ndvi.fs$NDVI)
qqp(L8ndvi.fs$NDVI,"norm")


#Sliding window model
library(climwin)
# ndvi.m = slidingwin(xvar=list(rain = rainall$amountcm),
#                     cdate = rainall$Date,
#                     bdate = L8ndvi.fs$Date,
#                     baseline = lm(NDVI~1,data=L8ndvi.fs),
#                     type="relative",
#                     cinterval="day",
#                     range=c(365,1),
#                     stat="sum",
#                     fun="lin")

#Save model file
#saveRDS(ndvi.m,here::here("RDS model outputs","rain_predicts_NDVI_season.rds"))

#Read in model file
ndvi.m = readRDS(here::here("Output files","Climwin model outputs","rain_predicts_NDVI_season.rds"))

#Get dataset
head(ndvi.m[[1]]$Dataset,n=10)
ndvi.m.dataset = ndvi.m[[1]]$Dataset

#Look at best model
ndvi.m[[1]]$BestModel
summary(ndvi.m[[1]]$BestModel)

#Plot relationship
plotbest(dataset = ndvi.m.dataset,
         bestmodel = ndvi.m[[1]]$BestModel, 
         bestmodeldata = ndvi.m[[1]]$BestModelData)

#Plot delta plot
plotdelta(dataset=ndvi.m.dataset)





####Does rainfall predict NDVI during the whole year?####


#Check distribution
hist(L8ndvi$NDVI)
qqp(L8ndvi$NDVI,"norm")


#Sliding window model
# ndvi.m.year = slidingwin(xvar=list(rain = rainall$amountcm),
#                     cdate = rainall$Date,
#                     bdate = L8ndvi$Date,
#                     baseline = lm(NDVI~1,data=L8ndvi),
#                     type="relative",
#                     cinterval="day",
#                     range=c(365,1),
#                     stat="sum",
#                     fun="lin")

#Save model file
#saveRDS(ndvi.m.year,here::here("RDS model outputs","rain_predicts_NDVI_year.rds"))

#Read in model file
ndvi.m.year = readRDS(here::here("Output files","Climwin model outputs","rain_predicts_NDVI_year.rds"))

#Get dataset
head(ndvi.m.year[[1]]$Dataset,n=10)
ndvi.m.year.dataset = ndvi.m.year[[1]]$Dataset

#Look at best model
ndvi.m.year[[1]]$BestModel
summary(ndvi.m.year[[1]]$BestModel)

#Plot relationship
plotbest(dataset = ndvi.m.year.dataset,
         bestmodel = ndvi.m.year[[1]]$BestModel, 
         bestmodeldata = ndvi.m.year[[1]]$BestModelData)

#Plot delta plot
plotdelta(dataset=ndvi.m.year.dataset)





####Compare observed to random####

#Run random window 100 times
# rand = randwin(repeats = 100,
#                window = "sliding",
#                xvar=list(rain = rainall$amountcm),
#                cdate = rainall$Date,
#                bdate = L8ndvi$Date,
#                baseline = lm(NDVI~1,data=L8ndvi),
#                type="relative",
#                cinterval="day",
#                range=c(365,1),
#                stat="sum",
#                fun="lin")



#Look at results
head(rand)
rand.dataset = rand[[1]]

#Save rand.dataset
#write.csv(rand.dataset,here::here("Output files","rain_predicts_NDVI_year_randomizations.csv"))

#Read in randomizations file
rand.dataset = read.csv(here::here("Output files","Climwin randomizations","rain_predicts_NDVI_year_randomizations.csv"))



#Plot histogram
library(cowplot)
hist(rand.dataset$deltaAICc,xlim = c(-50,0),xlab = "delta-AICc", main="")
abline(v=ndvi.m.year.dataset$deltaAICc[1],col="red",lwd=3)



