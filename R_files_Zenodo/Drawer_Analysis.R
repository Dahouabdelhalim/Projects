# packages
library("plyr")
library("lme4")
library("car")
library("multcomp")
library("MuMIn")
library(parallel)
# functions (provided by Roger Mundry)
source("/.../diagnostic_fcns.r")
source("/.../glmm_stability.r")
source("/.../boot_glmm.r")

# getting the data
setwd("/...")
test.data=read.csv(file="drawer_analysis_data.csv", header=T, sep=";", dec=",")

########################################################################################################## ##
## Q1: Do the two partners and/or their performance have different effects? (CA=co-action, CO=competition) ##
########################################################################################################## ##
##############
### for CA ###
##############
# make subset with only drawerset=CA and manipulation=NORM or SLOW
d=droplevels(subset(test.data, drawerset=="CA"))
d=droplevels(subset(d, manipulation=="normal"|manipulation == "slow"))

# log-transform responses to approach normality
d$log.lat=log(d$subj.lat)
hist(d$log.lat)

# z-transform predictors
d$total_trial_number=as.numeric(d$total_trial_number)
d$z.trial=as.vector(scale(d$total_trial_number))
d$z.partn.lat=as.vector(scale(d$partn.lat))
d$z.threat=as.vector(scale(d$threat_duration)) 

# random slope structure (function by Roger Mundry)
xx.fe.re=fe.re.tab(fe.model="log.lat ~ partn.lat + partner +
                   z.trial", re="(1|subject)",
                   other.vars=NULL, data=d
)
xx.fe.re$summary # include all random slopes

# dummy code factors
d$partner.code= as.numeric(d$partner== levels(d$partner)[2])
d$partner.code= d$partner.code - mean(d$partner.code)

# run the random slopes model
res=lmer(log.lat~ z.partn.lat + partner + z.trial +
           (1 + z.partn.lat + partner.code + z.trial ||subject),
         data = d, REML=F)

# check assumptions, residuals (function by Roger Mundry)
diagnostics.plot(res) # a bit off on the right end but still  acceptable 
# assumptions, random effects
ranef.diagn.plot(res) # ok

# check for model stability (function by Roger Mundry)
# check if ranges are very high, this indicates that the respective effects are unstable
m.stab=glmm.model.stab(model.res=res, contr=NULL)
table(m.stab$detailed$warnings) # No warning
round(m.stab$summary[,-1],3)
#                           orig    min    max
# (Intercept)              0.900  0.874  0.929
# z.partn.lat             -0.026 -0.035 -0.002
# partnersnickers         -0.088 -0.119 -0.026
# z.trial                  0.036 -0.006  0.086
# subject@(Intercept)@NA   0.042  0.000  0.091
# subject@z.partn.lat@NA   0.000  0.000  0.000
# subject@partner.code@NA  0.080  0.000  0.119
# subject@z.trial@NA       0.155  0.067  0.167
# Residual                 0.333  0.286  0.348
m.stab.plot(m.stab$summary [,-1]) 

# check for collinearity
# use fixed effects model to use variance inflation factor function
xres=lm(log.lat~ z.partn.lat + partner + z.trial + z.threat , data=d)
vif(xres) # max vif=1.3

# building the null model
null=lmer(log.lat~  z.trial +
            (1+ z.partn.lat + partner.code + z.trial ||subject),
          data = d, REML=F)

anova(null, res, test="Chisq") # the models are not different
# thus it looks like neither the partner nor their performance affect subjects' latencies
#      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# null  7 136.42 157.90 -61.209   122.42                         
# res   9 137.06 164.68 -59.529   119.06 3.3586      2     0.1865

# effect size
r.squaredGLMM(x=res)
# R2m        R2c 
# 0.03100209 0.22149630 

##############
### for CO ###
##############
# make subset with only drawerset=CO and manipulation=NORM or SLOW
d_co=droplevels(subset(test.data, drawerset=="CO"))
d_co=droplevels(subset(d_co, manipulation=="normal"|manipulation == "slow"))

# log-transform response
d_co$log.lat=log(d_co$subj.lat)
hist(d_co$log.lat)

# z-transform predictors
d_co$total_trial_number=as.numeric(d_co$total_trial_number)
d_co$z.trial=as.vector(scale(d_co$total_trial_number))
d_co$z.partn.lat=as.vector(scale(d_co$partn.lat))
d_co$z.threat=as.vector(scale(d_co$threat_duration)) 

# random slope structure
xx.fe.re=fe.re.tab(fe.model="log.lat ~ partn.lat + partner +
                   z.trial", re="(1|subject)",
                   other.vars=NULL, data=d_co
)
xx.fe.re$summary # include all random slopes

# dummy code factors
d_co$partner.code= as.numeric(d_co$partner== levels(d_co$partner)[2])
d_co$partner.code= d_co$partner.code - mean(d_co$partner.code)

# run the random slopes model
res=lmer(log.lat~ z.partn.lat + partner + z.trial +
           (1 + z.partn.lat + partner.code + z.trial ||subject),
         data = d_co, REML=F)

# check assumptions, residuals
diagnostics.plot(res) # ok
# assumptions, random effects
ranef.diagn.plot(res) #  ok

# check for model stability
# check if ranges are very high, this indicates that the respective effects are unstable
m.stab=glmm.model.stab(model.res=res, contr=NULL)
table(m.stab$detailed$warnings) # No warning
round(m.stab$summary[,-1],3)
#                           orig    min    max
# (Intercept)              0.889  0.823  0.946
# z.partn.lat              0.070  0.048  0.092
# partnersnickers         -0.154 -0.200 -0.074
# z.trial                 -0.249 -0.300 -0.193
# subject@(Intercept)@NA   0.224  0.180  0.260
# subject@z.partn.lat@NA   0.000  0.000  0.000
# subject@partner.code@NA  0.121  0.000  0.153
# subject@z.trial@NA       0.148  0.100  0.186
# Residual                 0.323  0.292  0.338
m.stab.plot(m.stab$summary [,-1]) 

# check for collinearity
# use fixed effects model to use variance inflation factor function
xres=lm(log.lat~ z.partn.lat + partner + z.trial + z.threat , data=d_co)
vif(xres) # max vif=1.3

# building the null model
null=lmer(log.lat~  z.trial +
            (1 + z.partn.lat + partner.code + z.trial ||subject),
          data = d_co, REML=F)

anova(null, res, test="Chisq") # the models are not different
# thus it looks like neither the partner nor their performance affect subjects' latencies
#      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# null  7 128.60 149.29 -57.301   114.60                         
# res   9 128.17 154.78 -55.087   110.17 4.4285      2     0.1092

# effect size
r.squaredGLMM(x=res)
# R2m       R2c 
# 0.3129429 0.6019132


## ############################################################# ##
## Main analysis: effects of social presence and drawer movement ##
## ############################################################# ##

# merge social (normal+slow) and SC (SC+SC2) conditions
d2=test.data
d2=droplevels(subset(test.data, manipulation != "BL"))

# aggregate normal and slow into "social"
d2$manip.help <- ifelse(d2$manipulation == "normal", "social",
                        ifelse(d2$manipulation == "slow", "social",
                               ifelse(d2$manipulation == "ghost", "ghost",
                                      ifelse(d2$manipulation == "SC", "SC",
                                             ifelse(d2$manipulation=="SC2", "SC",0)))))
d2$manip.help = as.factor(d2$manip.help)

# exclude subject Moritz because we do not have data for CO social for him
d2=subset(d2, subject != "moritz")
d2=droplevels(d2)

### ##################### ###
### Analysis CA condition ###
### ##################### ###
# make subset for Co-Action 
d2_CA=droplevels(subset(d2, drawerset=="CA")) 

# log-transfrom responses 
d2_CA$log.lat=log(d2_CA$subj.lat)
hist(d2_CA$log.lat)

# z-transform predictors
d2_CA$total_trial_number=as.numeric(d2_CA$total_trial_number)
d2_CA$z.trial=as.vector(scale(d2_CA$total_trial_number))

# design structure (random slopes)
xx.fe.re=fe.re.tab(fe.model="log.lat ~ manip.help + 
                   z.trial", re="(1|subject)",
                   other.vars=NULL, data=d2_CA)
xx.fe.re$summary # both random slopes can be included

# preparations for random slope model
# manually dummy-code factors
d2_CA$manipulation.code2= as.numeric(d2_CA$manip.help== levels(d2_CA$manip.help)[2])
d2_CA$manipulation.code2= d2_CA$manipulation.code2 - mean(d2_CA$manipulation.code2)  
d2_CA$manipulation.code3= as.numeric(d2_CA$manip.help== levels(d2_CA$manip.help)[3])
d2_CA$manipulation.code3= d2_CA$manipulation.code3 - mean(d2_CA$manipulation.code3)  

# run the random slopes model
res=lmer(log.lat ~ manip.help + z.trial + 
           (1+ manipulation.code2 + manipulation.code3 + z.trial ||subject),
         data = d2_CA, REML=F)

# check assumptions, residuals
diagnostics.plot(res) # ok
# assumptions, random effects
ranef.diagn.plot(res) # ok

# check for model stability
m.stab=glmm.model.stab(model.res=res, contr=NULL)
table(m.stab$detailed$warnings) # no warnings -> stable
round(m.stab$summary[,-1],3)
#                                 orig    min    max
# (Intercept)                    0.964  0.900  1.029
# manip.helpghost                0.014 -0.061  0.053
# manip.helpsocial              -0.103 -0.147 -0.069
# z.trial                       -0.098 -0.129 -0.079
# subject@(Intercept)@NA         0.157  0.118  0.176
# subject@manipulation.code2@NA  0.121  0.032  0.142
# subject@manipulation.code3@NA  0.000  0.000  0.071
# subject@z.trial@NA             0.040  0.000  0.062
# Residual                       0.326  0.287  0.341
m.stab.plot(m.stab$summary [,-1]) 

# check for collinearity
# use fixed effects model to use variance inflation factor function
xres=lm(log.lat~manip.help + z.trial , data=d2_CA)
vif(xres) 
# GVIF Df GVIF^(1/(2*Df))
# manip.help 1.005697  2        1.001421
# z.trial    1.005697  1        1.002844

# building the null model
null=lmer(log.lat~ z.trial +
            (1+ manipulation.code2 + manipulation.code3 + z.trial ||subject),
          data = d2_CA, REML=F)

# compare models
anova(null, res, test="Chisq") # models are not different
#      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# null  7 192.94 217.62 -89.473   178.94                         
# res   9 192.83 224.56 -87.417   174.83 4.1117      2      0.128

round(summary(res)$coefficients, 3)

as.data.frame(drop1(res,test="Chisq"))
#            Df      AIC      LRT     Pr(Chi)
# <none>     NA 192.8335       NA          NA
# manip.help  2 192.9452 4.111658 0.127986661
# z.trial     1 198.1105 7.277044 0.006984138

# effect size
r.squaredGLMM(res)
# R2m       R2c 
# 0.0809194 0.2797592 

# Confidence intervals (function by Roger Mundry)
boot.res=boot.glmm.pred(model.res = res, excl.warnings = T, nboots = 1000, para = T)
table(unlist(lapply(boot.res$all.warns, length)))
round(boot.res$ci.estimates, 3)
#                    orig  X2.5. X97.5.
# (Intercept)       0.964  0.823  1.106
# manip.helpghost   0.014 -0.142  0.162
# manip.helpsocial -0.103 -0.222  0.013
# z.trial          -0.098 -0.155 -0.038

## Post-hoc pairwise comparisons
aux <- glht(model = res, linfct= c("manip.helpghost = 0",
                                   "manip.helpsocial = 0",
                                   "manip.helpghost - manip.helpsocial = 0")) 

confint(aux)
summary(aux)

par(mar=c(10, 12, 3, 0.5)) # define graphics window to avoid cutting y-axis labels
plot(aux, cex.lab = 1,
     cex.axis = 1)

### ##################### ###
### Analysis CO condition ###
### ##################### ###
# make subset for Competition 
d2_CO=droplevels(subset(d2, drawerset=="CO")) 

# check responses 
hist(d2_CO$subj.lat)
min(d2_CO$subj.lat, na.rm = T)
d2_CO$log.lat=log(d2_CO$subj.lat)
hist(d2_CO$log.lat)

# z-transform predictors
d2_CO$total_trial_number=as.numeric(d2_CO$total_trial_number)
d2_CO$z.trial=as.vector(scale(d2_CO$total_trial_number))

# design structure (random slopes)
xx.fe.re=fe.re.tab(fe.model="log.lat ~ manip.help + 
                   z.trial", re="(1|subject)",
                   other.vars=NULL, data=d2_CO)
xx.fe.re$summary #both can be included

# preparations for random slope model
# manually dummy-code factors
d2_CO$manipulation.code2= as.numeric(d2_CO$manip.help== levels(d2_CO$manip.help)[2])
d2_CO$manipulation.code2= d2_CO$manipulation.code2 - mean(d2_CO$manipulation.code2)  
d2_CO$manipulation.code3= as.numeric(d2_CO$manip.help== levels(d2_CO$manip.help)[3])
d2_CO$manipulation.code3= d2_CO$manipulation.code3 - mean(d2_CO$manipulation.code3)  

# run the random slopes model
res=lmer(log.lat ~ manip.help + z.trial + 
           (1+ manipulation.code2 + manipulation.code3 + z.trial ||subject),
         data = d2_CO, REML=F)

# check assumptions, residuals
diagnostics.plot(res) # looks fine
# assumptions, random effects
ranef.diagn.plot(res) 

# check for model stability
m.stab=glmm.model.stab(model.res=res, contr=NULL)
table(m.stab$detailed$warnings) # no warnings 
round(m.stab$summary[,-1],3)
#                                 orig    min    max
# (Intercept)                    1.050  0.992  1.123
# manip.helpghost                0.119  0.073  0.160
# manip.helpsocial              -0.225 -0.266 -0.190
# z.trial                       -0.081 -0.090 -0.076
# subject@(Intercept)@NA         0.226  0.170  0.242
# subject@manipulation.code2@NA  0.192  0.136  0.221
# subject@manipulation.code3@NA  0.085  0.040  0.104
# subject@z.trial@NA             0.000  0.000  0.000
# Residual                       0.318  0.304  0.329
m.stab.plot(m.stab$summary [,-1]) 

# check for collinearity
# use fixed effects model to use variance inflation factor function
xres=lm(log.lat~manip.help + z.trial , data=d2_CO)
vif(xres) 
# GVIF Df GVIF^(1/(2*Df))
# manip.help 1.007774  2        1.001938
# z.trial    1.007774  1        1.003880

# building the null model
null=lmer(log.lat~ z.trial +
            (1+manipulation.code2 + manipulation.code3 + z.trial ||subject),
          data = d2_CO, REML=F)

# compare models
anova(null, res, test="Chisq") # models are different
#      Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# null  7 198.80 223.42 -92.401   184.80                            
# res   9 190.75 222.41 -86.378   172.75 12.047      2   0.002421 **

round(summary(res)$coefficients, 3)
#                   Estimate Std. Error t value
# (Intercept)         1.050      0.095  11.028
# manip.helpghost     0.119      0.091   1.305
# manip.helpsocial   -0.225      0.066  -3.416
# z.trial            -0.081      0.027  -3.048

as.data.frame(drop1(res,test="Chisq"))
#            Df      AIC       LRT     Pr(Chi)
# <none>     NA 190.7554        NA          NA
# manip.help  2 198.8023 12.046871 0.002421337
# z.trial     1 196.7852  8.029811 0.004601365

ddply(d2_CO, c("drawerset", "manip.help"), summarise,
      meanLatency = mean(subj.lat, na.rm = T))
#   drawerset manip.help meanLatency
# 1        CO         SC    3.057466
# 2        CO      ghost    3.456375
# 3        CO     social    2.594903

# Confidence intervals
boot.res=boot.glmm.pred(model.res = res, excl.warnings = T, nboots = 1000, para = T)
table(unlist(lapply(boot.res$all.warns, length)))
round(boot.res$ci.estimates, 3)
#                    orig  X2.5. X97.5.
# (Intercept)       1.050  0.866  1.226
# manip.helpghost   0.119 -0.054  0.292
# manip.helpsocial -0.225 -0.352 -0.088
# z.trial          -0.081 -0.135 -0.029

# effect size
r.squaredGLMM(x=res)
# R2m       R2c 
# 0.1519094 0.4688086 

## Post-hoc pairwise comparisons
aux <- glht(model = res, linfct= c("manip.helpghost = 0",
                                   "manip.helpsocial = 0",
                                   "manip.helpghost - manip.helpsocial = 0")) 
                                   

confint(aux)
# Linear Hypotheses:
#                                          Estimate lwr      upr     
# manip.helpghost == 0                     0.11936 -0.09405  0.33277
# manip.helpsocial == 0                   -0.22454 -0.37792 -0.07117
# manip.helpghost - manip.helpsocial == 0  0.34390  0.14781  0.53999
summary(aux)
#Linear Hypotheses:
# Linear Hypotheses:
#                                         Estimate Std. Error z value Pr(>|z|)    
# manip.helpghost == 0                     0.11936    0.09145   1.305 0.387576    
# manip.helpsocial == 0                   -0.22454    0.06573  -3.416 0.001657 ** 
# manip.helpghost - manip.helpsocial == 0  0.34390    0.08403   4.093 0.000101 ***


## ########################### ##
## Comparison to Baseline (BL) ##
## ########################### ##
# merge social and SC conditions
d3=test.data

# aggregate normal and slow into "social"
d3$manip.help <- ifelse(d3$manipulation == "normal", "social",
                        ifelse(d3$manipulation == "slow", "social",
                               ifelse(d3$manipulation == "ghost", "ghost",
                                      ifelse(d3$manipulation == "SC", "SC",
                                             ifelse(d3$manipulation=="SC2", "SC",
                                                    ifelse(d3$manipulation=="BL", "BL",
                                                           0))))))
d3$manip.help = as.factor(d3$manip.help)

# exclude Moritz because we do not have data for CO social for him
d3=subset(d3, subject != "moritz")
d3=droplevels(d3)

### ##################### ###
### Analysis CA condition ###
### ##################### ###
# make subset for Co-Action 
d3_CA=droplevels(subset(d3, drawerset=="CA")) 

# response -> log-transform
d3_CA$log.lat=log(d3_CA$subj.lat)
# z-transform predictors
d3_CA$total_trial_number=as.numeric(d3_CA$total_trial_number)
d3_CA$z.trial=as.vector(scale(d3_CA$total_trial_number))

# preparations for random slope model
# manually dummy-code factors
d3_CA$manipulation.code2= as.numeric(d3_CA$manip.help== levels(d3_CA$manip.help)[2])
d3_CA$manipulation.code2= d3_CA$manipulation.code2 - mean(d3_CA$manipulation.code2)  
d3_CA$manipulation.code3= as.numeric(d3_CA$manip.help== levels(d3_CA$manip.help)[3])
d3_CA$manipulation.code3= d3_CA$manipulation.code3 - mean(d3_CA$manipulation.code3)  

# run the random slopes model
res=lmer(log.lat ~ manip.help + z.trial + 
           (1+ manipulation.code2 + manipulation.code3 + z.trial ||subject),
         data = d3_CA, REML=F)

##  compare each condition with BL
aux <- glht(model = res, linfct= c("manip.helpghost = 0",
                                   "manip.helpSC = 0",
                                   "manip.helpsocial = 0"))
confint(aux)
# Linear Hypotheses:
#                        Estimate lwr      upr     
# manip.helpghost == 0   0.09784 -0.09140  0.28707
# manip.helpSC == 0      0.09870 -0.10439  0.30180
# manip.helpsocial == 0 -0.01009 -0.17243  0.15226

summary(aux)
# Linear Hypotheses:
#                        Estimate Std. Error z value Pr(>|z|)
# manip.helpghost == 0   0.09784    0.08241   1.187    0.435
# manip.helpSC == 0      0.09870    0.08845   1.116    0.480
# manip.helpsocial == 0 -0.01009    0.07070  -0.143    0.997

## none of the comparisons with BL is significant

### ##################### ###
### Analysis CO condition ###
### ##################### ###
# make subset for Co-Action 
d3_CO=droplevels(subset(d3, drawerset=="CO")) 

# response -> log-transform
d3_CO$log.lat=log(d3_CO$subj.lat)
# z-transform predictors
d3_CO$total_trial_number=as.numeric(d3_CO$total_trial_number)
d3_CO$z.trial=as.vector(scale(d3_CO$total_trial_number))

# preparations for random slope model
# manually dummy-code factors
d3_CO$manipulation.code2= as.numeric(d3_CO$manip.help== levels(d3_CO$manip.help)[2])
d3_CO$manipulation.code2= d3_CO$manipulation.code2 - mean(d3_CO$manipulation.code2)  
d3_CO$manipulation.code3= as.numeric(d3_CO$manip.help== levels(d3_CO$manip.help)[3])
d3_CO$manipulation.code3= d3_CO$manipulation.code3 - mean(d3_CO$manipulation.code3)  

# run the random slopes model
res=lmer(log.lat ~ manip.help + z.trial + 
           (1+ manipulation.code2 + manipulation.code3 + z.trial ||subject),
         data = d3_CO, REML=F)

##  compare each condition with BL
aux <- glht(model = res, linfct= c("manip.helpghost = 0",
                                   "manip.helpSC = 0",
                                   "manip.helpsocial = 0"))
confint(aux)
# Linear Hypotheses:
#                        Estimate lwr      upr     
# manip.helpghost == 0   0.12877 -0.10468  0.36222
# manip.helpSC == 0      0.01512 -0.17887  0.20911
# manip.helpsocial == 0 -0.20739 -0.37814 -0.03665

summary(aux)
# Linear Hypotheses:
#                        Estimate Std. Error z value Pr(>|z|)  
# manip.helpghost == 0   0.12877    0.10084   1.277   0.3952  
# manip.helpSC == 0      0.01512    0.08380   0.180   0.9943  
# manip.helpsocial == 0 -0.20739    0.07376  -2.812   0.0126 *

## only the comparison BL-competition is significant
