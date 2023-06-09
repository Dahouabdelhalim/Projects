library(ggplot2)
require(dplyr)

library(lmerTest)
library(MuMIn)
library(emmeans)

## Test1.  #############################
## The difference of otolith increment widths 
## among JP, CA, SA south-east and SA west subpopulations.

All<-read.csv('IW_10-100_2.csv',header = T)# Choose "IW_10-100_2.csv"
All$Fish.ID<-as.factor(All$Fish_ID)#set Fish.ID as factors
All$Region<-as.factor(All$Region)#set Region as factors
All$Age<-as.factor(All$Age)#set Age as factors

LMM1<-lmer(IW ~ Age*Region +(1|Fish.ID),All)
ems1=emmeans(LMM1, pairwise~Region|Age) #pairwise comparison 
pair1=pairs(ems1)
EF1=eff_size(ems1, sigma = sigma(LMM1), edf = Inf)
write.csv(ems1,'ems1.csv')
write.csv(pair1,'pair1.csv')
write.csv(EF1,'EF1.csv')
#Plot Supplementary Figure 6
par(mfrow=c(2,2))
qqnorm(resid(LMM1))
ypred = predict(LMM1)
res = residuals(LMM1, type = 'deviance')
plot(ypred,res)
hist(res)

## Test1 Ends. #############################