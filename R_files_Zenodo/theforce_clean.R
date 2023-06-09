# FORCE SCRIPT IN R
remove(collarstats)
rm(list=ls())

library(car)
library(lme4)
library(MuMIn)
library(lmodel2)

theforce<-read.table("/Users/Kayleigharose/Desktop/theforce.txt",header=T)
attach(theforce)
names(theforce)
plot(theforce)
str(theforce)

# does body mass influence relative force experienced by carrier ?
modpercentage<-lmer((percentagebodymass)~(bodymass)*category*gait+(1|name), REML=TRUE)
car::Anova(modpercentage, test.statistic="F")
modpercentagea=update(modpercentage,~.-bodymass:category:gait)
car::Anova(modpercentagea, test.statistic="F")


summary(modpercentagea)
coef(modpercentagea)
r.squaredGLMM(modpercentagea)
AIC(modpercentagea)
plot(modpercentagea)
qqnorm(resid(modpercentagea))
qqline(resid(modpercentagea))
hist(residuals(modpercentagea))

