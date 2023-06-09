# MAC WILD ANIMAL R SCRIPT

rm(list=ls())
wild<-read.table("/Users/Kayleigharose/Desktop/wildanimals.txt",header=T)
attach(wild)
names(wild)

# packages
library(car)
library(lme4)
library(MuMIn)

str(wild)
# amplitude with body mass
modsize<-lmer(amp~bodymass*gait+(1|id), REML=TRUE)
car::Anova(modsize, test.statistic="F")
modsize2=update(modsize,~.-bodymass:gait)
car::Anova(modsize2, test.statistic="F")
coef(modsize2)
r.squaredGLMM(modsize2)
AIC(modsize2)
plot(modsize2)
qqnorm(resid(modsize2))
qqline(resid(modsize2))
hist(residuals(modsize2))
# logged cause residuals not normally distributed 
modsizel<-lmer(log(amp)~log(bodymass)*gait+(1|id), REML=TRUE)
car::Anova(modsizel, test.statistic="F")
modsizel2=update(modsizel,~.-log(bodymass):gait)
car::Anova(modsizel2, test.statistic="F")
coef(modsizel2)
r.squaredGLMM(modsizel2)
AIC(modsizel2)
plot(modsizel2)
qqnorm(resid(modsizel2))
qqline(resid(modsizel2))
hist(residuals(modsizel2))

# period - nice 
modsizeperiod<-lmer(period~bodymass*gait+(1|id), REML=TRUE)
car::Anova(modsizeperiod, test.statistic="F")
coef(modsizeperiod)
r.squaredGLMM(modsizeperiod)
AIC(modsizeperiod)
plot(modsizeperiod)
qqnorm(resid(modsizeperiod))
qqline(resid(modsizeperiod))
hist(residuals(modsizeperiod))


# everything in one model
modeverything<-lmer(amp~period*gait+bodymass+(1|id), REML=TRUE)
car::Anova(modeverything, test.statistic="F")
coef(modeverything)
r.squaredGLMM(modeverything)
AIC(modeverything)
plot(modeverything)
qqnorm(resid(modeverything))
qqline(resid(modeverything))
hist(residuals(modeverything))
# needs logging 
modeverythingloggged<-lmer(log(amp)~log(period)*gait+bodymass+(1|id), REML=TRUE)
car::Anova(modeverythingloggged, test.statistic="F")
modeverythingloggged2=update(modeverythingloggged,~.-log(period):gait)
car::Anova(modeverythingloggged2, test.statistic="F")
coef(modeverythingloggged2)
r.squaredGLMM(modeverythingloggged2)
AIC(modeverythingloggged2)
plot(modeverythingloggged2)
qqnorm(resid(modeverythingloggged2))
qqline(resid(modeverythingloggged2))
hist(residuals(modeverythingloggged2))



modperiod<-lmer(period~bodymass*gait+(1|id), REML=TRUE)
car::Anova(modperiod, test.statistic="F")
coef(modperiod)
r.squaredGLMM(modperiod)
AIC(modperiod)
plot(modperiod)
qqnorm(resid(modperiod))
qqline(resid(modperiod))
hist(residuals(modperiod))
