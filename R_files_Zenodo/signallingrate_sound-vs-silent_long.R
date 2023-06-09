library(nlme)
library(emmeans)
xdata <- read.csv(file = "data_onto_playback_center.csv")

########################
#########Format#########
########################

xdata$z.age <- scale(xdata$age)
xdata$z.dose <- scale(xdata$dose)
xdata$z.temp <- scale(xdata$temperature)
levels(xdata$treatment) <- c("sound", "silent")
xdata$id <- as.factor(xdata$id)

str(xdata)

####################################
###########Model fitting############
####################################

res <- lme(log.long~z.temp+I(z.temp^2)+sex*(z.age+I(z.age^2))+treatment*z.dose, random=~1|id, data=xdata)

####################################
############Assumptions#############
####################################

#Normality of residuals and independance
source('diagnostic_fcns.r')
diagnostics.plot(res) # seems fine

## Colinearity
library(car)
res.lm <- lm(log.long~z.dose+sex+z.age+temperature+treatment, data=xdata)
vif(res.lm)# Ok, <4

## Stability
stab <- influence(res, maxfun=20)
stability <- t(rbind(stab$fixed.effects, apply(stab$`fixed.effects[-case]`, FUN=quantile, probs=c(0.05, 0.95), MARGIN=2)))
#which are the unstable estimates
unstable <- which(sign(stability[,1])!=sign(stability[,2])|sign(stability[,1])!=sign(stability[,3]))


anova(res)
