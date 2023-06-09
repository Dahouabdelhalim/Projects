rm(list=ls())
library (nlme)

##  Daten einlesen
setwd("D:/ConsumerDemandPaper/VisitDuration")
VD=read.delim2("VisitDuration.txt")

str(VD)
attach(VD)
names(VD)
VD$RequiredNosepokeNumber<-as.factor(VD$RequiredNosepokeNumber)

## Model
M3 <- lme (log(VisitDuration) ~ RequiredNosepokeNumber, random= ~ 1 | Animal/RequiredNosepokeNumber, data= VD, method= 'ML')
summary (M3)
anova (M3, type= 'marginal')

qqnorm(residuals(M3)) 
plot(fitted(M3),residuals(M3))