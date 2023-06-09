######################################################################################################################################
## Script by Arianne F. Messerman
## 13 July, 2021
## This script is adapted from the code of James T. Stroud, and calculates selection coefficients for mean residual respiratory surface
##   area water loss (RSAWL), log-body mass, and residual standard metabolic rate (SMR) values for individual juvenile spotted and 
##   marbled salamanders (Ambystoma maculatum and A. opacum, respectively).
## Please see the associated manuscript for full study details: 
##   Messerman, A.F. and M. Leal. The contributions of individual traits to survival among terrestrial juvenile pond-breeding 
##   salamanders. Functional Ecology.
######################################################################################################################################
setwd()
load() #This script requires the environment from Resp_Survival_2017-2018_Sept2021.R

#Find selection coefficients using ordinary least squares regression
amma$rel.surv.oct <- amma$surv.oct/mean(amma$surv.oct)
amma$rel.surv.apr <- amma$surv.may/mean(amma$surv.may)
amop$rel.surv.oct <- amop$surv.oct/mean(amop$surv.oct)
amop$rel.surv.apr <- amop$surv.may/mean(amop$surv.may)

#Linear selection coefficients
sc.amma.oct <- lm(rel.surv.oct ~ mass.std + resid.std + resid.mr.std, data=amma)
#Take the estimates for each linear term (presented as Beta) and standard error as variance
summary(sc.amma.oct)

sc.amop.oct <- lm(rel.surv.oct ~ mass.std + resid.std + resid.mr.std, data=amop)
summary(sc.amop.oct)

sc.amma.apr <- lm(rel.surv.apr ~ mass.std + resid.std + resid.mr.std, data=amma)
summary(sc.amma.apr)

sc.amop.apr <- lm(rel.surv.apr ~ mass.std + resid.std + resid.mr.std, data=amop)
summary(sc.amop.apr)

#Nonlinear selection coefficients
sc.amma.oct1 <- lm(rel.surv.oct ~ mass.std + I(mass.std^2) + resid.std + I(resid.std^2)+ resid.mr.std + I(resid.mr.std^2), data=amma)
#Take the estimates for each nonlinear term (presented as Gamma) and standard error as variance, take linear estimates from models above
summary(sc.amma.oct1)

sc.amop.oct1 <- lm(rel.surv.oct ~ mass.std + I(mass.std^2) + resid.std + I(resid.std^2)+ resid.mr.std + I(resid.mr.std^2), data=amop)
summary(sc.amop.oct1)

sc.amma.apr1 <- lm(rel.surv.apr ~ mass.std + I(mass.std^2) + resid.std + I(resid.std^2)+ resid.mr.std + I(resid.mr.std^2), data=amma)
summary(sc.amma.apr1)

sc.amop.apr1 <- lm(rel.surv.apr ~ mass.std + I(mass.std^2) + resid.std + I(resid.std^2)+ resid.mr.std + I(resid.mr.std^2), data=amop)
summary(sc.amop.apr1)
#Double your nonlinear coefficients and standard errors when you report them (you don't do this for linear estimates)