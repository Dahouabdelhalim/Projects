rm(list=ls())
library(nlme)

##  Daten einlesen
setwd("D:/ConsumerDemandPaper/Consumer DemandAuswertung")
Curves=read.delim2("data_CosnumerDemandCurves.txt")
str(Curves)
attach(Curves)
names(Curves)

## Modell
cd1<-lme(log(DrinkingEvents+0.5)~log(RequiredNosepokeNumber)*Run, random = ~1|Animal/Run, data=Curves)
cd1
summary(cd1)

qqnorm(residuals(cd1)) 
plot(fitted(cd1),residuals(cd1)) 