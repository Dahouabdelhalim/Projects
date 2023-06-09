
##A diet rich in C3 plants reveals the sensitivity of an alpine mammal to climate change



#Sabuj Bhattacharyya1, 2, Deborah A. Dawson2, Helen Hipperson2, Farah Ishtiaq1


#1Centre for Ecological Sciences, Indian Institute of Science, Bangalore, India
#2Department of Animal and Plant Sciences, Western Bank, Sheffield, S10 2TN, UK


#*Corresponding author:  Sabuj Bhattacharyya, bhattacharyyasabuj@gmail.com 


##R code for GLM models 


setwd("/Volumes/SABUJ/Barcode new analysis2/SR_GLMM")
data <-read.csv("SR_data.csv",header=TRUE)
attach(data)

names(data)
## to check how response variable distribution is
plot(data$Genus.richness.in.diet)
hist(data$Genus.richness.in.diet)
### histogram is left skewed
### to check the normality
shapiro.test(data$Genus.richness.in.diet)
##Shapiro-Wilk normality test

data:  data$Genus.richness.in.diet
##W = 0.92849, p-value = 7.758e-06
## results is significant or  the p-value is lower than 0.05, we can conclude that the sample deviates from normality.
## hence we should use GLM/glmm
## Z tranformation
ELET<-scale(data$ELE, scale = T,center = T)
data$ELET <- ELET
ROT<-scale(data$RO, scale = T,center = T)
data$ROT <- ROT
BAT<-scale(data$BA, scale = T,center = T)
data$BAT <- BAT

GRT<-scale(data$GR, scale = T,center = T)
data$GRT <- GRT

FOT<-scale(data$FO, scale = T,center = T)
data$FOT <- FOT

SHT<-scale(data$SH, scale = T,center = T)
data$SHT <- SHT

TRT<-scale(data$TR, scale = T,center = T)
data$TRT <- TRT
GFT<-scale(data$GF, scale = T,center = T)
data$GFT <- GFT

TVCT<-scale(data$TVC, scale = T,center = T)
data$TVCT <- TVCT

TAT<-scale(data$TA, scale = T,center = T)
data$TAT <- TAT

NTT<-scale(data$NT, scale = T,center = T)
data$NTT <- NTT

SFT<-scale(data$SF, scale = T,center = T)
data$SFT <- SFT
CDT<-scale(data$CD, scale = T,center = T)
data$CDT <- CDT
names(data)
library(lme4)
library(MuMIn)
library(sjPlot)

### SR= species richness, 
##ROCK VARIABELS: TAT=total area of talus, ROT=rock cover,NTT= distance to nearest talus,CDT= depth of crevices,
### VEGETATION VARAIBLES:FOT= forbes cover,TVCT=total vegetation cover, SHT=SHRUB COVER. TRT=TREE COVER,GRT=GRASS COVER
### TOPOGRAPHIC VARAIABLES:ASP=aspect,SLO=slope of plot,ELET=elevation
###### 
## 
## CDT AND TAT IS COREALATED
##TVC IS CORELATED WITH OTHER VEGETATION


## FOOD SELECTION IS GOVEREN BY PREDATION RISK(E.G.TALUS CHARECTERISTICS), FOOD AVILABILITY (E.G. FOBES COVER) AND TOPOHRAPHY
# Global Model WHICH HAS ALL NON CORELATED VAROABLES
gm1<-glm(Genus.richness.in.diet~TAT*CDT+NTT+ROT+FOT+SHT+TRT+GRT+ELET+ASP+SLO,data=data, family= poisson(link = "log")) 

m1<-glm(Genus.richness.in.diet~TAT+NTT+ROT+FOT+SHT+TRT+GRT+ELET+ASP+SLO,data=data, family= poisson(link = "log")) 

m2<-glm(Genus.richness.in.diet~TAT*FOT+TAT*SHT+ TAT*TRT+TAT*GRT+NTT+ROT+ELET+ASP+SLO,data=data, family= poisson(link = "log")) 

m3<-glm(Genus.richness.in.diet~TAT+NTT+ROT+ELET*FOT+ELET*SHT+ ELET*TRT+GRT*ELET+ASP+SLO,data=data, family= poisson(link = "log")) 

m4<-glm(Genus.richness.in.diet~TAT+NTT+ROT+ASP*FOT+ASP*SHT+ASP*TRT+ASP*GRT+ELET+SLO,data=data, family= poisson(link = "log")) 

m5<-glm(Genus.richness.in.diet~TAT+NTT+ROT+SLO*FOT+SLO*SHT+SLO*TRT+SLO*GRT+ELET+ASP,data=data, family= poisson(link = "log")) 

m6<-glm(Genus.richness.in.diet~TAT*CDT+NTT+ROT+FOT+SHT+TRT+GRT+ELET+ASP,data=data, family= poisson(link = "log")) 

m7<-glm(Genus.richness.in.diet~TAT*CDT+NTT+ROT+FOT+SHT+TRT+GRT+ELET,data=data, family= poisson(link = "log")) 

m8<-glm(Genus.richness.in.diet~TAT*CDT+NTT+ROT+FOT+SHT+GRT+,data=data, family= poisson(link = "log")) 

m9<-glm(Genus.richness.in.diet~TAT*CDT+NTT+ROT+FOT+SHT, data=data, family= poisson(link = "log"))

m10<-glm(Genus.richness.in.diet~TAT*CDT+NTT+ROT+FOT,data=data, family= poisson(link = "log"))  

m11<-glm(Genus.richness.in.diet~TAT*CDT+NTT+ROT,data=data, family= poisson(link = "log")) 

m12<-glm(Genus.richness.in.diet~TAT*CDT+NTT,data=data, family= poisson(link = "log")) 

m13<-glm(Genus.richness.in.diet~TAT*CDT,data=data, family= poisson(link = "log")) 

m14<-glm(Genus.richness.in.diet~NTT+ROT+FOT+SHT+GRT+TRT+ELET+ASP+SLO,data=data, family= poisson(link = "log")) 

m15<-glm(Genus.richness.in.diet~ROT+FOT+SHT+GRT+ELET+ASP+SLO,data=data, family= poisson(link 
                                                                                        = "log")) 

m16<-glm(Genus.richness.in.diet~FOT+SHT+GRT+ELET+ASP+SLO,data=data, family= poisson(link = 
                                                                                      "log")) 

m17<-glm(Genus.richness.in.diet~SHT+GRT+ELET+ASP+SLO,data=data, family= poisson(link = "log")) 

m18<-glm(Genus.richness.in.diet~GRT+ELET+ASP+SLO,data=data, family= poisson(link = "log")) 

m19<-glm(Genus.richness.in.diet~ELET+ASP+SLO,data=data, family= poisson(link = "log")) 

m20<-glm(Genus.richness.in.diet~ASP+SLO,data=data, family= poisson(link = "log")) 

m21<-glm(Genus.richness.in.diet~SLO,data=data, family= poisson(link = "log")) 

m22<-glm(Genus.richness.in.diet~NTT,data=data, family= poisson(link = "log")) 

m23<-glm(Genus.richness.in.diet~ROT, data=data, family= poisson(link = "log")) 

m24<-glm(Genus.richness.in.diet~FOT,data=data, family= poisson(link = "log")) 

m25<-glm(Genus.richness.in.diet~SHT,data=data, family= poisson(link = "log")) 

m26<-glm(Genus.richness.in.diet~TRT,data=data, family= poisson(link = "log")) 

m27<-glm(Genus.richness.in.diet~GRT,data=data, family= poisson(link = "log")) 

m28<-glm(Genus.richness.in.diet~ELET,data=data, family= poisson(link = "log"))

m29<-glm(Genus.richness.in.diet~ASP,data=data, family= poisson(link = "log"))  

m30<-glm(Genus.richness.in.diet~SLO,data=data, family= poisson(link = "log")) 

m31< glm(Genus.richness.in.diet~TAT,data=data, family= poisson(link = "log")) 

mnull<-glm(Genus.richness.in.diet~1,data=data,family=poisson(link="log"))

out.put<-model.sel(mnull,m1,m2,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31)
out.put



## Model Average
fm<-model.avg(m12,m13)
summary(fm)
## confidence interval
confint(fm)

library(effects)
plot(allEffects(m13))
plot(allEffects(m12))

