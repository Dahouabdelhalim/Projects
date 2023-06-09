# MAC
rm(list=ls())
collarstats<-read.table("/Users/Kayleigharose/Desktop/collarstats.txt",header=T)
attach(collarstats)
names(collarstats)
plot(collarstats)

# packages
library(car)
library(lme4)
library(MuMIn)

# check variables
str(collarstats)

### A - does body mass influence period between peaks?
modb<-lmer((averageperiod)~(bodymass)*gait+(1|name), REML=TRUE)
car::Anova(modb, test.statistic="F")
plot(modb)
summary(modb)
plot(modb)
qqnorm(resid(modb))
qqline(resid(modb))
plot(modb, resid(., scaled=TRUE) ~ fitted(.) | name, abline = 0)
plot(modb, name ~ resid(., scaled=TRUE))
plot(modb, averagespeed ~ fitted(.) | name, abline = c(0,1))
plot(modb, resid(., scaled=TRUE) ~ bodymass | name, abline = 0)
require("lattice")
qqmath(modb, id=0.05)
hist(residuals(modb))
coef(modb)
summary(modb)
r.squaredGLMM(modb)



########### B Was gait-specific travel speed influenced by dog body mass or tag % body mass? #######
modspeed<-lmer((averagespeed)~(bodymass)*as.numeric(category_num)*gait+(1|name), REML=TRUE)
car::Anova(modspeed, test.statistic="F")
modspeeda=update(modspeed,~.-(bodymass):as.numeric(category_num):gait)
car::Anova(modspeeda, test.statistic="F")
modspeedb=update(modspeeda,~.-(bodymass):as.numeric(category_num))
car::Anova(modspeedb, test.statistic="F")
modspeedc=update(modspeedb,~.-as.numeric(category_num):gait)
car::Anova(modspeedc, test.statistic="F")
plot(modspeedc)
qqnorm(resid(modspeedc))
qqline(resid(modspeedc))
plot(modspeedc, resid(., scaled=TRUE) ~ fitted(.) | name, abline = 0)
plot(modspeedc, name ~ resid(., scaled=TRUE))
plot(modspeedc, averagespeed ~ fitted(.) | name, abline = c(0,1))
plot(modspeedc, resid(., scaled=TRUE) ~ bodymass | name, abline = 0)
require("lattice")
qqmath(modspeedc, id=0.05)
hist(residuals(modspeedc))
coef(modspeedc)
summary(modspeedc)
r.squaredGLMM(modspeedc)


########### C Does acceleration increase with speed?#######
mod1<-lmer((maxpeakvectorialsum)~0+(averagespeed)*category*gait+bodymass+(1|name), REML=TRUE)
car::Anova(mod1, test.statistic="F")
mod2=update(mod1,~.- averagespeed:category:gait)
car::Anova(mod2, test.statistic="F")
r.squaredGLMM(mod2)
summary(mod2)
coef(mod2)
r.squaredGLMM(mod2)
plot(mod2)
qqnorm(resid(mod2))
qqline(resid(mod2))
plot(mod2, resid(., scaled=TRUE) ~ fitted(.) | name, abline = 0)
plot(mod2, name ~ resid(., scaled=TRUE))
plot(mod2, averagespeed ~ fitted(.) | name, abline = c(0,1))
plot(mod2, resid(., scaled=TRUE) ~ averagespeed | name, abline = 0)
require("lattice")
qqmath(mod2, id=0.05)
hist(residuals(mod2))


######## FORCE - different script 



