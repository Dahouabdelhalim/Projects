rm(list=ls())
library(lme4)
library(nlme)
library("emmeans")
library(multcomp)
library(ggplot2)
setwd('/Users/cdesjonq/Documents/Rédactionarticles/1_in-prep/2021_Desjonqu_socialhybrids/submission/211014_PNAS/data/')

females <- read.csv("signals-preferences.csv", dec=',')      

head(females)
str(females)

#remove males
females <- females[!is.na(females$strength),]

#full model fitting ----
females$log.strength <- log(females$strength)
res <- lmer(log.strength~sp*treat+year+z.temp+(1|aggregation), data=females)


# Assumptions
#Normality of residuals and independance
source('/Users/cdesjonq/Documents/stats_glmm/functions/diagnostic_fcns.r')
diagnostics.plot(res) # seems fine

## Colinearity
library(car)
res.lm <- lm(strength~sp+treat+year+z.temp, data=females)
vif(res.lm)# Ok, <4

## Stability
### excluding individual cases and levels of random effect, one at a time
source('/Users/cdesjonq/Documents/stats_glmm/functions/glmm_stability.r')
m.stab.ind <- glmm.model.stab(model.res=res, ind.cases=T)
warnings()#38 warnings... : In min(x, na.rm = na.rm) : aucun argument trouvé pour min ; Inf est renvoyé
round(m.stab.ind$summary[,-1], digits=2)
#a <- m.stab.ind$summary[,-1]#stable except for within temp
#Inference
drop1(res, test="Chisq")

#sp
null1 <- lmer(log.strength~sp+treat+year+z.temp+(1|aggregation), data=females)
null <- lmer(log.strength~treat+year+z.temp+(1|aggregation), data=females)
anova(null1,null, test="Chisq")

#treat
null <- lmer(log.strength~sp+year+z.temp+(1|aggregation), data=females)
anova(null1,null, test="Chisq")

#random effects

##aggreg
null <- lm(log.strength~sp*treat+year+z.temp, data=females)
anova(res,null, test="Chisq")

#Figure ----

coefs <- fixef(res)
randeff <- ranef(res)
females$ranef <- NA
for (i in 1:length(females$aggregation)){
  females$ranef[i] <- randeff$aggregation[which(row.names(randeff$aggregation)==females$aggregation[i]), 1]
}
females$pred.data <- females$log.strength-(females$z.temp*coefs[6]+coefs[4]*as.numeric(females$year=='2019')+coefs[5]*as.numeric(females$year=='2020')+females$ranef)


means_m <- tapply(females$pred.data, INDEX=paste(females$sp,females$treat, females$loc), FUN=mean)
sd_m <- tapply(females$pred.data, INDEX=paste(females$sp,females$treat, females$loc), FUN=sd)
n_m <- tapply(females$pred.data, INDEX=paste(females$sp,females$treat, females$loc), FUN=length)
factors <- names(means_m)

se_m <- sd_m/sqrt(n_m)

ci_inf_m <- means_m-se_m
ci_sup_m <- means_m+se_m

bg <- ifelse(substr(factors, start=4, stop=5)=='ho', yes='blue', no='orange')

col <- ifelse(females$treat=='homo', yes='aquamarine2', no='orange')
col[females$loc%in%c('BOG', 'OLT', 'PNV')&females$treat=='homo'] <- 'blue'
col[females$loc%in%c('BOG', 'OLT', 'PNV')&females$treat=='hete'] <- 'brown'

bg <- ifelse(substr(factors, start=4, stop=5)=='ho', yes='aquamarine2', no='orange')
bg[substr(factors, start=9, stop=9)!='F'&substr(factors, start=4, stop=5)=='ho'] <- 'blue'
bg[substr(factors, start=9, stop=9)!='F'&substr(factors, start=4, stop=5)=='he'] <- 'brown'

index <- ifelse(substr(factors, start=1, stop=2)=='LF', yes=1, no=3)+ifelse(substr(factors, start=4, stop=5)=='ho', yes=0, no=1)


females$index <- 0
females$index[females$treat=='homo'&females$sp=='LF'] <- 1
females$index[females$treat=='hete'&females$sp=='LF'] <- 2
females$index[females$treat=='homo'&females$sp=='HF'] <- 3
females$index[females$treat=='hete'&females$sp=='HF'] <- 4

#ylim <- range(c(females$pred.data,means_m, ci_inf_m, ci_sup_m))
ylim <- range(c(means_m, ci_inf_m, ci_sup_m))
bitmap("females_socialhybrids_strength.jpg", width = 2400, height = 2400, units = 'px', res=900)
par(mar=c(4.1, 4.1, 1.1, 1.1))
plot(index,means_m, xaxt='n', ylab='preference strength (log)', xlab='', col=bg, pch=16, ylim=ylim, xlim=c(0.9, 4.3), frame=FALSE)
#points(jitter(females$index, 0.4), females$pred.data, col=alpha(col, 0.5), pch=16, cex=0.6)
arrows(x0 = index, y0 = ci_sup_m, x1 =index, y1 = ci_inf_m, length = 0, angle=90, code=3, col=bg)

het <- which(substr(factors,start= 4,stop=5)=='he')
hom <- which(substr(factors,start= 4,stop=5)=='ho')

arrows(x0=index[hom], x1=index[het],y0=means_m[hom], y1=means_m[het], length = 0, angle=90, code=3, col='grey', lty=3)
axis(side = 1, at = c(1, 1.5, 3.5, 4), labels = c('',expression(sp[low]), expression(sp[high]), ''))
mtext(side = 1, at = 2.5, text = c('species'), line = 2, font=1.5)
legend('bottomright', legend = c('mixed allopatric', 'mixed sympatric', 'own allopatric', 'own sympatric'), col=c('brown', 'orange', 'blue', 'aquamarine2'), pch=16, bty='n')
dev.off()


##LF ----
females_LF <- females[females$sp=='LF',]
res <- lmer(log.strength~treat*location+year+z.temp+(1|aggregation), data=females_LF)

#Asumptions
##Normality of residuals and independance
source('/Users/cdesjonq/Documents/stats_glmm/functions/diagnostic_fcns.r')
diagnostics.plot(res) # seems fine

## Colinearity
library(car)
res.lm <- lm(strength~location+temp+treat+sin+cos+daysintreat, data=females_LF)
vif(res.lm)# Ok, <4

## Stability
### excluding individual cases and levels of random effect, one at a time
source('/Users/cdesjonq/Documents/stats_glmm/functions/glmm_stability.r')
m.stab.ind <- glmm.model.stab(model.res=res, ind.cases=T)
warnings()#38 warnings... : In min(x, na.rm = na.rm) : aucun argument trouvé pour min ; Inf est renvoyé
m.stab.ind$summary[,-1]
a <- m.stab.ind$summary[,-1]#stable except for within temp


#Inferences
summary(res)
drop1(res, test="Chisq")

#location
null1 <- lmer(log.strength~treat+location+year+z.temp+(1|aggregation), data=females_LF)
null <- lmer(log.strength~treat+year+z.temp+(1|aggregation), data=females_LF)
anova(null, null1, test='Chisq')

#treat
null <- lmer(log.strength~location+year+z.temp+(1|aggregation), data=females_LF)
anova(null, null1, test='Chisq')

#random effects
#aggreg
null <- lm(log.strength~treat*location+year+z.temp, data=females_LF)
anova(res, null, test='Chisq')


