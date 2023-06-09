library(ggplot2)
require(dplyr)

library(lmerTest)
library(MuMIn)
library(emmeans)

## Test2.  #############################
## The difference of Moto  
## between JP and CA subpopulations.


All<-read.csv('Moto_re_stage_3_R3.csv',header = T) # Choose "Moto_re_stage_3_R3.csv"
All$Fish.ID<-as.factor(All$Fish.ID)#set Fish.ID as factors
All$Region<-as.factor(All$Region)#set Region as factors
All$Age.range2<-All$Age.range # Adjust for the different definitions for age range between JP and CA
All$Age.range2[(All$Region=='JP')&(All$Age.range=='C')]='B'
All$Age.range2[(All$Region=='JP')&(All$Age.range=='D')]='C'
All$Age.range2[(All$Region=='JP')&(All$Age.range=='E')]='C'
All$Age.range2[(All$Region=='JP')&(All$Age.range=='F')]='D'
All$Age.range2[(All$Region=='JP')&(All$Age.range=='G')]='D'

All$Age<-as.factor(All$Age.range2)#set Age.range2 as factors
All=subset(All,All$Age.range2!='E')#Remove 121-150 dph data of CA as there is no data for JP

hist(All$Moto.mean)
# Remove outliers prior to the test
boxplot(All$Moto.mean, plot=FALSE)$out
outliers <- boxplot(All$Moto.mean, plot=FALSE)$out
AllM<- All[-which(All$Moto.mean %in% outliers),]
hist(AllM$Moto.mean)

LMM2<-lmer(Moto.mean ~ Age*Region +(1|Fish.ID), AllM)
summary(LMM2)


ems2=emmeans(LMM2, pairwise~Region|Age)#pairwise comparison 
pair2=pairs(ems2)
EF2=eff_size(ems2, sigma = sigma(LMM2), edf = Inf)
write.csv(ems2,'ems2.csv')
write.csv(pair2,'pair2.csv')
write.csv(EF2,'EF2.csv')
#evaluate the normality of residuals (Supplementary Figure 7)
par(mfrow=c(2,2))
qqnorm(resid(LMM2))
ypred = predict(LMM2)
res = residuals(LMM2, type = 'deviance')
plot(ypred,res)
hist(res)


## Test2 Ends.  #############################


## Test3.  #############################
## The difference of experienced water temperature  
## between JP and CA subpopulations.

hist(All$Estimated.Temperature)
# Remove outliers prior to the test
boxplot(All$Estimated.Temperature, plot=FALSE)$out
outliers <- boxplot(All$Estimated.Temperature, plot=FALSE)$out
AllT<- All[-which(All$Estimated.Temperature %in% outliers),]
hist(AllT$Estimated.Temperature)

LMM3<-lmer(Estimated.Temperature ~ Age*Region +(1|Fish.ID), AllT)
summary(LMM3)

ems3=emmeans(LMM3, pairwise~Region|Age)#pairwise comparison 
pair3=pairs(ems3)
EF3=eff_size(ems3, sigma = sigma(LMM3), edf = Inf)
write.csv(ems3,'ems3.csv')
write.csv(pair3,'pair3.csv')
write.csv(EF3,'EF3.csv')
#evaluate the normality of residuals (Supplementary Figure 8)
par(mfrow=c(2,2))
qqnorm(resid(LMM3))
ypred = predict(LMM3)
res = residuals(LMM3, type = 'deviance')
plot(ypred,res)
hist(res)


## Test3 Ends.  #############################