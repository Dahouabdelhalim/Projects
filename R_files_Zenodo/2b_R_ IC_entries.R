rm(list=ls())
library (nlme)

##  Daten einlesen
setwd("D:/ConsumerDemandPaper/IC_entries")
ICE=read.delim2("data_ IC_entries.txt")

names(ICE)

ICE$Tag<-as.factor(ICE$Tag)
ICE$Phase<-as.factor(ICE$Phase)
str(ICE)

contrasts (ICE [, 'Phase']) <- contr.sum (2)


ICE [, 'DaySinceCleaningN'] <- (ICE [, 'DaySinceCleaning'] - mean (ICE [, 'DaySinceCleaning'])) / sd (ICE [, 'DaySinceCleaning'])

ICE [, 'DayN'] <- (ICE [, 'Day'] - mean (ICE [, 'Day'])) / sd (ICE [, 'Day'])

summary(ICE)

## Model
M2 <- lme(log(sum_entries) ~ DayN * Phase + DaySinceCleaningN, random = ~1|Tag/Phase, data=ICE, method='ML')
summary(M2)
anova (M2, type= 'marginal')

Entries.res = resid(M2)
plot(ICE$DaySinceCleaningN, Entries.res, ylab="Residuals", xlab="DaySinceCleaningN") 

scatter.smooth (resid (M2), ICE [, 'DaySinceCleaningN'])
with(ICE, scatter.smooth(DaySinceCleaningN, resid(M2)))


qqnorm(residuals(M2)) 
plot(fitted(M2),residuals(M2))