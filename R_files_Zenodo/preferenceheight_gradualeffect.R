library(lme4)

xdata <- read.csv(file="femalepref_dose.csv", sep=',')

####################################
########## Model fitting ###########
####################################

res <- lm(peak_height~log.dose*log.pract+z.temp, data=xdata)

####################################
########### Assumptions ############
####################################

#Normality of residuals and independance
source('diagnostic_fcns.r')
diagnostics.plot(res)

## Colinearity
library(car)
res.lm <- lm(peak_height~dose+log.pract+z.temp+z.times.diff, data=xdata)
vif(res.lm)# Ok, <4

## Stability
### excluding individual cases and levels of random effect, one at a time
max(abs(dffits(res)))#Must be <2
hist(dffits(res))# no extremely influencial cases
round(cbind(coefficients(res), coefficients(res)+t(apply(X=dfbeta(res), MARGIN=2, FUN=range))), 5)

####################################
########### Inferences #############
####################################
anova(res)
