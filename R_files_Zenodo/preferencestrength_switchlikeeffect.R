library(lme4)

xdata <- read.csv(file="femalepref_dose.csv", sep=',')

####################################
########## Model fitting ###########
####################################

res <- lm(strength~exp.switch+pract.switch+z.temp, data=xdata)

####################################
########### Assumptions ############
####################################

#Normality of residuals and independance
source('/Users/cdesjonq/Documents/stats_glmm/functions/diagnostic_fcns.r')
diagnostics.plot(res) # okish

## Colinearity
library(car)
res.lm <- lm(strength~exp.switch+pract.switch+z.temp, data=xdata)
vif(res.lm)# Ok, <4

## Stability
### excluding individual cases and levels of random effect, one at a time

max(abs(dffits(res)))#Must be <2
hist(dffits(res))# no extremely influencial cases
round(cbind(coefficients(res), coefficients(res)+t(apply(X=dfbeta(res), MARGIN=2, FUN=range))), 5)#temp and z.cos not stable


####################################
########### Inferences #############
####################################

summary(res)
