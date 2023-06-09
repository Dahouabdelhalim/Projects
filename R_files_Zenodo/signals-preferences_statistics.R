rm(list=ls())
library(lme4)
library(nlme)
library("emmeans")
library(multcomp)
library(RLRsim)

setwd('/Users/cdesjonq/Documents/Rédactionarticles/1_in-prep/2021_Desjonqu_socialhybrids/submission/211014_PNAS/data/')

xdata <- read.csv(file = 'signals-preferences.csv', dec=',')

#full model fitting ----
res <- lmer(peak_pref~sp*sex*treat+year+z.temp+(1|aggregation), data=xdata,REML=FALSE)

# Assumptions
#Normality of residuals and independance
source('/Users/cdesjonq/Documents/stats_glmm/functions/diagnostic_fcns.r')
diagnostics.plot(res) # seems fine

## Colinearity
library(car)
res.lm <- lm(peak_pref~sp+treat+year+sex+z.temp, data=xdata)
vif(res.lm)# Ok, <4

## Stability
### excluding individual cases and levels of random effect, one at a time
source('/Users/cdesjonq/Documents/stats_glmm/functions/glmm_stability.r')
m.stab.ind <- glmm.model.stab(model.res=res, ind.cases=T)
warnings()#38 warnings... : In min(x, na.rm = na.rm) : aucun argument trouvé pour min ; Inf est renvoyé
round(m.stab.ind$summary[,-1], digits=2)

#Inference
drop1(res, test="Chisq")

#sp*sex
null1 <- lmer(peak_pref~sp*sex+sex*treat+sp*treat+year+z.temp+(1|aggregation), data=xdata)
null <- lmer(peak_pref~sex*treat+sp*treat+year+z.temp+(1|aggregation), data=xdata)
anova(null,null1, test="Chisq")
#treat*sex
null <- lmer(peak_pref~sp*sex+sp*treat+year+z.temp+(1|aggregation), data=xdata)
anova(null,null1, test="Chisq")
#sp*treat
null <- lmer(peak_pref~sp*sex+sex*treat+year+z.temp+(1|aggregation), data=xdata)
anova(null,null1, test="Chisq")
#sex
null1 <- lmer(peak_pref~sex+sp*treat+year+z.temp+(1|aggregation), data=xdata)
null <- lmer(peak_pref~sp*treat+year+z.temp+(1|aggregation), data=xdata)
anova(null,null1, test="Chisq")
#treat
null1 <- lmer(peak_pref~treat+sp*sex+year+z.temp+(1|aggregation), data=xdata)
null <- lmer(peak_pref~sp*sex+year+z.temp+(1|aggregation), data=xdata)
anova(null,null1, test="Chisq")
#sp
null1 <- lmer(peak_pref~sp+sex*treat+year+z.temp+(1|aggregation), data=xdata)
null <- lmer(peak_pref~sex*treat+year+z.temp+(1|aggregation), data=xdata)
anova(null,null1, test="Chisq")
#aggregation
null <- lm(peak_pref~sp*sex*treat+year+z.temp, data=xdata)
anova(res,null, test="Chisq")

## compare m0 and m1
exactLRT(res,null)

summary(res)

#LF ----
xdataLF <- xdata[xdata$sp=='LF',]
res <- lmer(peak_pref~sex*treat+treat*location+year+z.temp+(1|aggregation), data=xdataLF)
# Assumptions
#Normality of residuals and independance
source('/Users/cdesjonq/Documents/stats_glmm/functions/diagnostic_fcns.r')
diagnostics.plot(res) # seems fine

## Colinearity
library(car)
res.lm <- lm(peak_pref~treat+location+sex+year+z.temp, data=xdata)
vif(res.lm)# Ok, <4

## Stability
### excluding individual cases and levels of random effect, one at a time
source('/Users/cdesjonq/Documents/stats_glmm/functions/glmm_stability.r')
m.stab.ind <- glmm.model.stab(model.res=res, ind.cases=T)
warnings()#38 warnings... : In min(x, na.rm = na.rm) : aucun argument trouvé pour min ; Inf est renvoyé
round(m.stab.ind$summary[,-1], digits=2)


## Distribution of the random effects
ranef.diagn.plot(res) #Looks ok... 

std.pred <- scale(fitted(res))
std.res <- scale(residuals(res))
plot(std.pred, std.res)
text(std.pred, std.res, row.names(std.res), cex=0.6, pos=4, col="red") 

ranef(res)

#Inference
drop1(res, test="Chisq")

#sex
null1 <- lmer(peak_pref~sex+location*treat+year+z.temp+(1|aggregation), data=xdataLF)
null <- lmer(peak_pref~location*treat+year+z.temp+(1|aggregation), data=xdataLF)
anova(null,null1, test="Chisq")

#treat
null1 <- lmer(peak_pref~sex+location+treat+year+z.temp+(1|aggregation), data=xdataLF)
null <- lmer(peak_pref~sex+location+year+z.temp+(1|aggregation), data=xdataLF)
anova(null,null1, test="Chisq")

#location
null1 <- lmer(peak_pref~sex*treat+location+year+z.temp+(1|aggregation), data=xdataLF)
null <- lmer(peak_pref~sex*treat+year+z.temp+(1|aggregation), data=xdataLF)
anova(null,null1, test="Chisq")

#random effects
##aggreg
null <- lm(peak_pref~sex*treat+treat*location+year+z.temp, data=xdataLF)
anova(res,null, test="Chisq")

summary(res)
