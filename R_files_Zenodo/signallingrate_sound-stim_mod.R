library(nlme)
library(emmeans)
xdata <- read.table(file="data_onto_playback.csv", sep=',', header=TRUE)

########################
######## Format ########
########################

xdata$z.age <- scale(xdata$age)
xdata$z.temp <- scale(xdata$temperature)
xdata$stimulus <- factor(xdata$stimulus, levels=levels(xdata$stimulus)[c(5,1:4)])
xdata$treatment <- factor(xdata$treatment, levels=levels(xdata$treat)[c(2,1)])

xdata$treatment <- factor(xdata$treatment, levels=c('TS','PB'))
levels(xdata$treatment) <- c('silent', 'sound')
xdata$stimulus <- factor(xdata$stimulus, levels=c('si','co','al', 'no', 'br'))
levels(xdata$stimulus) <- c('none', 'short', 'long', 'noise', 'brush')

xdata$id <- as.factor(xdata$id)

xdata <- xdata[xdata$treat=='sound',]
str(xdata)

####################################
########## Model fitting ###########
####################################

res <- lme(log.mod~time.in.rec+stimulus+sex*(z.age+I(z.age^2))+z.temp+I(z.temp^2), random=~1|id, data=xdata)#, family='poisson'

####################################
############Assumptions#############
####################################

#Normality of residuals and independance
source('diagnostic_fcns.r')
diagnostics.plot(res) # seems fine

## Colinearity
library(car)
res.lm <- lm(log.mod~z.age+sex+z.temp+stimulus+time.in.rec, data=xdata)
vif(res.lm)# Ok, <4

## Stability
stab <- influence(res, maxfun=20)
stability <- t(rbind(stab$fixed.effects, apply(stab$`fixed.effects[-case]`, FUN=quantile, probs=c(0.05, 0.95), MARGIN=2)))
#which are the unstable estimates
unstable <- which(sign(stability[,1])!=sign(stability[,2])|sign(stability[,1])!=sign(stability[,3]))

####################################
############ Inference #############
####################################

anova(res)

#contrasts

round(summary(res)$tTable[3:6,], digits=3)

xdata$stimulus <- factor(xdata$stimulus, levels=levels(xdata$stimulus)[c(2:5,1)])
res <- lme(log.mod~time.in.rec+stimulus+sex*(z.age+I(z.age^2))+z.temp+I(z.temp^2), random=~1|id, data=xdata)
round(summary(res)$tTable[3:5,], digits=3)

xdata$stimulus <- factor(xdata$stimulus, levels=levels(xdata$stimulus)[c(2:5,1)])
res <- lme(log.mod~time.in.rec+stimulus+sex*(z.age+I(z.age^2))+z.temp+I(z.temp^2), random=~1|id, data=xdata)
round(summary(res)$tTable[3:4,], digits=3)

xdata$stimulus <- factor(xdata$stimulus, levels=levels(xdata$stimulus)[c(2:5,1)])
res <- lme(log.mod~time.in.rec+stimulus+sex*(z.age+I(z.age^2))+z.temp+I(z.temp^2), random=~1|id, data=xdata)
round(summary(res)$tTable[3,], digits=3)

