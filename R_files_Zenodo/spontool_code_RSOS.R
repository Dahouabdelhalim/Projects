########################spontaneous tool use in children
##### Authors: K Neldner, E Reindl, C Tennie, J Grant, K Tomaselli, M Nielsen
##### 
###DVs: tool pickup, correct use, correct success and incorrect success
## presented in that order

#####
####LOADING DATA AND CODING VARIABLES - GENERAL PREPARATION


rm(list=ls())

library(tidyverse)
spontool <- read.csv("spont_tool_data_RSOS.csv")
glimpse(spontool)

#install.packages("lme4")
library(lme4)
#install.packages("psych")
library(psych)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("car")
library(car)
#install.packages("glmmTMB")
library(glmmTMB)


glimpse(spontool)
head(spontool)
str(spontool)

spontool$Gender <- factor(spontool$Gender, levels = c(0,1), labels = c("Male","Female"))
summary(spontool$Gender)
spontool$Culture <- factor(spontool$Culture, levels = c(0,1), labels = c("SA", "BNE"))
summary(spontool$Culture)
spontool$Frequency <- factor(spontool$Frequency, levels = c(0,1), labels = c("low", "high"))
summary(spontool$Frequency)

spontool$Task <- factor(spontool$Task, levels = c(1,2,3,4,5,6,7,8,9,10,11,12), labels = c("MP", "FD", "LO", "TF", "IP", "P", "NH", "AS", "SE", "TFLMR", "GP", "ADW"))
summary(spontool$Task)


str(spontool)
View(spontool)

spontool$Frequency <-as.numeric (spontool$Frequency)
head(spontool$Frequency)
##see that low is 1 and high is 2, need to recode to low being 0 and high being 1
spontool$Frequency[spontool$Frequency==1] <- 0
head(spontool$Frequency)
spontool$Frequency[spontool$Frequency==2] <- 1
head(spontool$Frequency)

str(spontool)
View(spontool)
(table(spontool$Task, spontool$Frequency))

Pick_upnoADW = spontool %>%
  mutate(Task = as.character(Task))%>%
  filter(Task != 'ADW')

table(Pick_upnoADW$Task)
table(spontool$Task)


spontool$Age_cent = spontool$Age - mean(spontool$Age)
str(spontool)
Pick_upnoADW$Age_cent = Pick_upnoADW$Age - mean(Pick_upnoADW$Age)
str(Pick_upnoADW)



table(spontool$Pick_up, spontool$Frequency)
table(spontool$Pick_up, spontool$Task)
table(spontool$Correct_success, spontool$Task)
table(spontool$Correct_use, spontool$Task)


###############
######CHECK LOCATIONS IN SA ARE EQUAL BEFORE COLLAPSING TO "SA" cultural group########
###############

###here we use only the spontool_PLAWIT data frame and the Location variable (not Culture)
spontool$Location <- factor(spontool$Location, levels = c(1,2,3), labels = c("PLA", "WIT", "BNE"))
summary(spontool$Location)
head(spontool)

spontool_PLAWIT <- spontool[which(spontool$Location == "PLA" | spontool$Location == "WIT"),]
#or in tidyverse 
table(spontool_PLAWIT$Location)

View(spontool_PLAWIT)

############################
#do for CORRECT USE
############################

Correct_use.Africaonly <- glmer(Correct_use ~ Gender + Age_cent + Frequency + Location + (1 + Frequency|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_use.Africaonly)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Correct_use.Africaonlyiter <- glmer(Correct_use ~ Gender + Age_cent + Frequency + Location + (1 + Frequency|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_use.Africaonlyiter)

Correct_use.Africaonlynocorrel <- glmer(Correct_use ~ Gender + Age_cent + Frequency + Location + (1|P_No) + (0 + Frequency|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_use.Africaonlynocorrel)


Correct_use.Africaonlynoslope <- glmer(Correct_use ~ Gender + Age_cent + Frequency + Location + (1|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_use.Africaonlynoslope)

Correct_use.Africaonlynoslope <- glmmTMB (Correct_use ~ Gender + Age_cent + Frequency + Location + (1|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_use.Africaonlynoslope)

Correct_use.Africaonlynullnoslope <- glmmTMB (Correct_use ~ Gender + (1|P_No), data=spontool_PLAWIT, family=binomial)
summary (Correct_use.Africaonlynullnoslope)

anova (Correct_use.Africaonlynullnoslope, Correct_use.Africaonlynoslope, test="Chisq")

###########################
######do for CORRECT SUCCESS
##########################

Correct_success.Africaonly <- glmer(Correct_success ~ Gender + Age_cent + Frequency + Location + (1 + Frequency|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_success.Africaonly)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Correct_success.Africaonlyiter <- glmer(Correct_success ~ Gender + Age_cent + Frequency + Location + (1 + Frequency|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_success.Africaonlyiter)

Correct_success.Africaonlynocorrel <- glmer(Correct_success ~ Gender + Age_cent + Frequency + Location + (1|P_No) + (0 + Frequency|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_success.Africaonlynocorrel)

Correct_success.Africaonlynoslope <- glmer(Correct_success ~ Gender + Age_cent + Frequency + Location + (1|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_success.Africaonlynoslope)

Correct_success.Africaonlynoslopenull <- glmer(Correct_success ~ Gender + (1|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_success.Africaonlynoslopenull)

anova (Correct_success.Africaonlynoslopenull, Correct_success.Africaonlynoslope, test="Chisq")

#Age
Correct_success.Africaonlynoslopenoage <- glmer(Correct_success ~ Gender + Frequency + Location + (1|P_No), data=spontool_PLAWIT, family=binomial)
anova (Correct_success.Africaonlynoslopenoage, Correct_success.Africaonlynoslope, test="Chisq")
#df = 1 chi sq = 8.592 p = .003  

#Location
Correct_success.Africaonlynoslopenolocat <- glmer(Correct_success ~ Gender + Frequency + Age_cent + (1|P_No), data=spontool_PLAWIT, family=binomial)
anova (Correct_success.Africaonlynoslopenolocat, Correct_success.Africaonlynoslope, test="Chisq")
#non-sig df = 1 chi sq 0.81 p = 0.369 *very slightly diff to last time but overall the same

#Frequency
Correct_success.Africaonlynoslopenofreq <- glmer(Correct_success ~ Gender + Location + Age_cent + (1|P_No), data=spontool_PLAWIT, family=binomial)
anova (Correct_success.Africaonlynoslopenofreq, Correct_success.Africaonlynoslope, test="Chisq")
# nonsig df = 1 p = 0.241 chi sq 1.37 

#Gender
Correct_success.Africaonlynoslopenogender <- glmer(Correct_success ~ Frequency + Location + Age_cent + (1|P_No), data=spontool_PLAWIT, family=binomial)
anova (Correct_success.Africaonlynoslopenogender, Correct_success.Africaonlynoslope, test="Chisq")
#non sig chi sq 0.231 p = 0.631 df = 1

#effect sizes/CIs
se <- sqrt(diag(vcov(Correct_success.Africaonlynoslope)))
# table of estimates with 95% CI
(tab <- cbind(Est = fixef(Correct_success.Africaonlynoslope), 
              LL = fixef(Correct_success.Africaonlynoslope) - 1.96 * se, 
              UL = fixef(Correct_success.Africaonlynoslope) + 1.96 * se))

#Odds ratios and their confidence intervals
exp(tab)

Correct_success.Africaonlynoslopeintx <- glmer(Correct_success ~ Gender + Age_cent*Frequency + Frequency*Location + Location*Age_cent + (1|P_No), data=spontool_PLAWIT, family=binomial)
summary(Correct_success.Africaonlynoslopeintx)

anova(Correct_success.Africaonlynoslopeintx, Correct_success.Africaonlynoslope)

#check assumptions of model
source("diagnostic_fcns.r")
overdisp.test(Correct_success.Africaonlynoslope)

#check assumptions of model - absence of collinearity
Correct_success.collAfricaonlynoslope <- lm(Correct_success~Gender + Age_cent + Frequency + Location, data=spontool_PLAWIT)
vif(Correct_success.collAfricaonlynoslope)

##check assumptions of model -distribution of random effect
ranef.diagn.plot(Correct_success.Africaonlynoslope)

#model stability
source("glmm_stability.r")
bin.stab=glmm.model.stab(model.res=Correct_success.Africaonlynoslope)
bin.stab$detailed$warnings
bin.stab$summary

library(DHARMa)
citation("DHARMa")

#scaled residuals
simulationOutput_Africaonlynoslope <- simulateResiduals(fittedModel = Correct_success.Africaonlynoslope, n = 250)
hist(simulationOutput_Africaonlynoslope$scaledResiduals)
plot(simulationOutput_Africaonlynoslope)

#uniformity, outliers, dispersion
testResiduals(simulationOutput_Africaonlynoslope)
testDispersion(simulationOutput_Africaonlynoslope)

###########################
#### TOOL PICKUP ##########
##########################
### here we use Culture instead of Location, and the Pick_upnoADW data frame

Pick_up.full <- glmer(Pick_up ~ Gender + Age_cent + Frequency + Culture + (1 + Frequency|P_No), data=Pick_upnoADW, family=binomial)
summary(Pick_up.full)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Pick_up.full <- glmer(Pick_up ~ Gender + Age_cent + Frequency + Culture + (1 + Frequency|P_No), data=Pick_upnoADW, family=binomial)
summary(Pick_up.full)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Pick_up.fullnocorrel <- glmer(Pick_up ~ Gender + Age_cent + Frequency + Culture + (1|P_No) + (0 + Frequency|P_No), data=Pick_upnoADW, family=binomial)
summary(Pick_up.fullnocorrel)

Pick_up.fullnoslope <- glmer(Pick_up ~ Gender + Age_cent + Frequency + Culture + (1|P_No), data=Pick_upnoADW, family=binomial)
summary(Pick_up.fullnoslope)

Pick_up.fullnullnoslope <- glmer(Pick_up ~ Gender + (1|P_No), data=Pick_upnoADW, family=binomial)

Pick_up.fullnoslopeTMB <- glmmTMB (Pick_up ~ Gender + Age_cent + Frequency + Culture + (1|P_No), data=Pick_upnoADW, family=binomial)
summary(Pick_up.fullnoslopeTMB)

Pick_up.fullnullnoslopeTMB <- glmmTMB(Pick_up ~ Gender + (1|P_No), data=Pick_upnoADW, family=binomial)
summary(Pick_up.fullnullnoslopeTMB)

anova(Pick_up.fullnullnoslopeTMB, Pick_up.fullnoslopeTMB)

#age
Pick_up.fullnoslopeTMBnoage <- glmmTMB(Pick_up ~ Gender + Frequency + Culture + (1|P_No), data=Pick_upnoADW, family=binomial)
anova(Pick_up.fullnoslopeTMBnoage, Pick_up.fullnoslopeTMB)

#culture
Pick_up.fullnoslopeTMBnocult <- glmmTMB(Pick_up ~ Gender + Frequency + Age_cent + (1|P_No), data=Pick_upnoADW, family=binomial)
anova(Pick_up.fullnoslopeTMBnocult, Pick_up.fullnoslopeTMB)
summary(Pick_up.fullnoslopeTMBnocult)

#frequency
Pick_up.fullnoslopeTMBnofreq <- glmmTMB(Pick_up ~ Gender + Culture + Age_cent + (1|P_No), data=Pick_upnoADW, family=binomial)
anova(Pick_up.fullnoslopeTMBnofreq, Pick_up.fullnoslopeTMB)
#converges chi sq 4.64, p 0.031  #same as glmer

#gender
Pick_up.fullnoslopeTMBnogender <- glmmTMB(Pick_up ~ Frequency + Culture + Age_cent + (1|P_No), data=Pick_upnoADW, family=binomial)
anova(Pick_up.fullnoslopeTMBnogender, Pick_up.fullnoslopeTMB)
#converges chi sq 2.44, p 0.118

#Confidence intervals
confint(Pick_up.fullnoslopeTMB)
exp(confint(Pick_up.fullnoslopeTMB))

Pick_up.fullnoslopeTMB.coll <- lm(Pick_up ~ Age_cent + Gender + Frequency + Culture, 
                                    data=Pick_upnoADW)
vif(Pick_up.fullnoslopeTMB.coll)

source("diagnostic_fcns.r")
overdisp.test(Pick_up.fullnoslopeTMB)

#DHARMa
simulationOutput_Pick_up.fullnoslopeTMB <- simulateResiduals(fittedModel = Pick_up.fullnoslopeTMB, n = 250)
hist(simulationOutput_Pick_up.fullnoslopeTMB$scaledResiduals) # is more or less uniform

#plot the scaled residuals
plot(simulationOutput_Pick_up.fullnoslopeTMB)

#uniformity, outliers, dispersion
testResiduals(simulationOutput_Pick_up.fullnoslopeTMB)
testDispersion(simulationOutput_Pick_up.fullnoslopeTMB)

###########################
#### CORRECT USE ##########
##########################
#we move onto spontool dataframe and use Culture again

Correct_use.full <- glmer(Correct_use ~ Gender + Age_cent + Frequency + Culture + (1 + Frequency|P_No), data=spontool, family=binomial)
summary(Correct_use.full)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Correct_use.fulliter <- glmer(Correct_use ~ Gender + Age_cent + Frequency + Culture + (1 + Frequency|P_No), data=spontool, family=binomial)
summary(Correct_use.fulliter)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Correct_use.fullnocorrel <- glmer(Correct_use ~ Gender + Age_cent + Frequency + Culture + (1|P_No) + (0 + Frequency|P_No), data=spontool, family=binomial)
summary(Correct_use.fullnocorrel)

Correct_use.fullnoslope <- glmer(Correct_use ~ Gender + Age_cent + Frequency + Culture + (1|P_No), data=spontool, family=binomial)
summary(Correct_use.fullnoslope)

Correct_use.fullnoslopenull <- glmer(Correct_use ~ Gender + (1|P_No), data=spontool, family=binomial)
summary(Correct_use.fullnoslopenull)
anova (Correct_use.fullnoslopenull, Correct_use.fullnoslope)

#Age
Correct_use.fullnoslopenoage <- glmer(Correct_use ~ Gender + Frequency + Culture + (1|P_No), data=spontool, family=binomial)
anova (Correct_use.fullnoslopenoage, Correct_use.fullnoslope)
#chi sq 3.19 df 1 p .074 

#Culture
Correct_use.fullnoslopenocult <- glmer(Correct_use ~ Gender + Frequency + Age_cent + (1|P_No), data=spontool, family=binomial)
anova (Correct_use.fullnoslopenocult, Correct_use.fullnoslope)
#chi sq 57.92 df 1 p < .001

#Frequency
Correct_use.fullnoslopenofreq <- glmer(Correct_use ~ Gender + Culture + Age_cent + (1|P_No), data=spontool, family=binomial)
anova (Correct_use.fullnoslopenofreq, Correct_use.fullnoslope)
#df 1 chi sq 6.00 .014

#Gender
Correct_use.fullnoslopenogend <- glmer(Correct_use ~ Frequency + Culture + Age_cent + (1|P_No), data=spontool, family=binomial)
anova (Correct_use.fullnoslopenogend, Correct_use.fullnoslope)
# chi sq 0.57 df 1 .451

#effect sizes/CIs
se <- sqrt(diag(vcov(Correct_use.fullnoslope)))
(tab <- cbind(Est = fixef(Correct_use.fullnoslope), 
              LL = fixef(Correct_use.fullnoslope) - 1.96 * se, 
              UL = fixef(Correct_use.fullnoslope) + 1.96 * se))

#Odds ratios
exp(tab)

Correct_use.fullnoslopeintx <- glmer(Correct_use ~ Gender + Age_cent*Frequency + Frequency*Culture + Culture*Age_cent + (1|P_No), data=spontool, family=binomial)
summary(Correct_use.fullnoslopeintx)

source("diagnostic_fcns.r")
overdisp.test(Correct_use.fullnoslope)

Correct_use.collfullnoslope <- lm(Correct_use ~ Gender + Age_cent + Frequency + Culture, data=spontool)
vif(Correct_use.collfullnoslope)

ranef.diagn.plot(Correct_use.fullnoslope)

source("glmm_stability.r")
bin.stab=glmm.model.stab(model.res=Correct_use.fullnoslope)
bin.stab$detailed$warnings
bin.stab$summary

simulationOutput_fullnoslope <- simulateResiduals(fittedModel = Correct_use.fullnoslope, n = 250)
hist(simulationOutput_fullnoslope$scaledResiduals) 
plot(simulationOutput_fullnoslope)

#uniformity, outliers, dispersion
testResiduals(simulationOutput_fullnoslope)
testDispersion(simulationOutput_fullnoslope)

###########################
#### CORRECT SUCCESS ######
##########################
#use spontool data frame and Culture fixed effects

Correct_success.full <- glmer(Correct_success ~ Gender + Age_cent + Frequency + Culture + (1 + Frequency|P_No), data=spontool, family=binomial)
summary(Correct_success.full)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Correct_success.fulliter <- glmer(Correct_success ~ Gender + Age_cent + Frequency + Culture + (1 + Frequency|P_No), data=spontool, family=binomial)
summary(Correct_success.fulliter)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Correct_success.fullnocorrel <- glmer(Correct_success ~ Gender + Age_cent + Frequency + Culture + (1|P_No) + (0 + Frequency|P_No), data=spontool, family=binomial)
summary(Correct_success.fullnocorrel)

Correct_success.fullnoslope <- glmer(Correct_success ~ Gender + Age_cent + Frequency + Culture + (1|P_No), data=spontool, family=binomial)
summary(Correct_success.fullnoslope)

Correct_success.noslopeTMB <- glmmTMB(Correct_success ~ Gender + Age_cent + Frequency + Culture + (1|P_No), data=spontool, family=binomial)
summary(Correct_success.noslopeTMB)

Correct_success.noslopeTMBnull <- glmmTMB(Correct_success ~ Gender + (1 |P_No), data=spontool, family=binomial)
summary(Correct_success.noslopeTMBnull)
anova(Correct_success.noslopeTMBnull, Correct_success.noslopeTMB)

#Age
Correct_success.noslopeTMBnoage <- glmmTMB(Correct_success ~ Gender + Frequency + Culture + (1|P_No), data=spontool, family=binomial)
anova(Correct_success.noslopeTMBnoage, Correct_success.noslopeTMB)
#chi sq 20.46 df 1 <.001

#Frequency
Correct_success.noslopeTMBnofreq <- glmmTMB(Correct_success ~ Gender + Age_cent + Culture + (1|P_No), data=spontool, family=binomial)
anova(Correct_success.noslopeTMBnofreq, Correct_success.noslopeTMB)
#chi sq 10.69 df 1 <.001

#Culture
Correct_success.noslopeTMBnocult <- glmmTMB(Correct_success ~ Gender + Frequency + Age_cent + (1|P_No), data=spontool, family=binomial)
anova(Correct_success.noslopeTMBnocult, Correct_success.noslopeTMB)
#chi sq 37.94 df 1 <.001

#Gender
Correct_success.noslopeTMBnogend <- glmmTMB(Correct_success ~ Age_cent + Frequency + Culture + (1|P_No), data=spontool, family=binomial)
anova(Correct_success.noslopeTMBnogend, Correct_success.noslopeTMB)
#chi sq .82 df 1 .364

#confidence intervals
confint(Correct_success.noslopeTMB)

#Odds ratios and their confidence intervals
exp(confint(Correct_success.noslopeTMB))

#add intx terms to full model
Correct_success.noslopeTMBintx <- glmmTMB(Correct_success ~ Gender + Age_cent*Frequency + Culture*Frequency + Culture*Age_cent + (1|P_No), data=spontool, family=binomial)
summary(Correct_success.noslopeTMBintx)
anova(Correct_success.noslopeTMBintx, Correct_success.noslopeTMB)

#check assumptions and stability
Correct_success.noslopeTMB.coll <- lm(Correct_success ~ Gender + Age_cent + Frequency + Culture, data=spontool)
vif(Correct_success.noslopeTMB.coll)

source("diagnostic_fcns.r")
overdisp.test(Correct_success.noslopeTMB)

#DHARMa
simulationOutput_Correct_success.noslopeTMB <- simulateResiduals(fittedModel = Correct_success.noslopeTMB, n = 250)
#for a correctly specified model we would expect a uniform (flat) distribution of the overall residuals
hist(simulationOutput_Correct_success.noslopeTMB$scaledResiduals) # is more or less uniform

#plot the scaled residuals
plot(simulationOutput_Correct_success.noslopeTMB)

#uniformity, outliers, dispersion
testResiduals(simulationOutput_Correct_success.noslopeTMB)
testDispersion(simulationOutput_Correct_success.noslopeTMB)

###########################
#### INCORRECT SUCCESS ####
##########################
#use spont tool dataframe and Culture fixed effect

Incorrect_success.full <- glmer(Incorrect_success  ~ Gender + Age_cent + Frequency + Culture + (1+ Frequency|P_No), 
                           data=spontool, family=binomial)
summary(Incorrect_success.full)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Incorrect_success.fulliter <- glmer(Incorrect_success  ~ Gender + Age_cent + Frequency + Culture + (1+ Frequency|P_No), 
                                data=spontool, family=binomial)
summary(Incorrect_success.fulliter)

contr=glmerControl(optCtrl=list(maxfun=10000000))
Incorrect_success.nocorrel <- glmer(Incorrect_success  ~ Gender + Age_cent + Frequency + Culture + (1|P_No) + (0 + Frequency|P_No), 
                                    data=spontool, family=binomial)
summary(Incorrect_success.nocorrel)

Incorrect_success.noslope <- glmer(Incorrect_success  ~ Gender + Age_cent + Frequency + Culture + (1|P_No), 
                                    data=spontool, family=binomial)
summary(Incorrect_success.noslope)

#check better than a null model
Incorrect_success.noslopenull <- glmer(Incorrect_success ~ Gender + (1|P_No), 
                                data=spontool, family=binomial)
summary(Incorrect_success.noslopenull)
anova(Incorrect_success.noslopenull,Incorrect_success.noslope)
#chi sq 8.5 df 3 p .037

#age
Incorrect_success.noslopenoage <- glmer(Incorrect_success  ~ Gender + Frequency + Culture + (1|P_No), 
                                   data=spontool, family=binomial)
summary(Incorrect_success.noslopenoage)
anova(Incorrect_success.noslopenoage, Incorrect_success.noslope)
#chi sq 0.68 df 1 p 0.41

#frequency
Incorrect_success.noslopenofreq <- glmer(Incorrect_success  ~ Gender + Age_cent + Culture + (1|P_No), 
                                        data=spontool, family=binomial)
summary(Incorrect_success.noslopenofreq)
anova(Incorrect_success.noslopenofreq, Incorrect_success.noslope)
#chi sq 7.22 df 1 .007

#culture
Incorrect_success.noslopenocult <- glmer(Incorrect_success  ~ Gender + Frequency + Age_cent + (1|P_No), 
                                        data=spontool, family=binomial)
summary(Incorrect_success.noslopenocult)
anova(Incorrect_success.noslopenocult, Incorrect_success.noslope)
#chi sq 0.53, df 1 p .466

#gender
Incorrect_success.noslopenogend <- glmer(Incorrect_success  ~ Age_cent + Frequency + Culture + (1|P_No), 
                                        data=spontool, family=binomial)
summary(Incorrect_success.noslopenogend)
anova(Incorrect_success.noslopenogend, Incorrect_success.noslope)
#chi sq 1.11 df 1 .292

#effect sizes/CIs
se <- sqrt(diag(vcov(Incorrect_success.noslope)))
# table of estimates with 95% CI
(tab <- cbind(Est = fixef(Incorrect_success.noslope), 
              LL = fixef(Incorrect_success.noslope) - 1.96 * se, 
              UL = fixef(Incorrect_success.noslope) + 1.96 * se))

#Odds ratios and their confidence intervals
exp(tab)

Incorrect_success.noslopeintx <- glmer(Incorrect_success  ~ Gender + Age_cent*Frequency + Culture*Frequency + Culture*Age_cent + (1|P_No), 
                                   data=spontool, family=binomial)
summary(Incorrect_success.noslopeintx)
anova(Incorrect_success.noslopeintx, Incorrect_success.noslope)

#check assumptions and stability
Incorrect_success.noslope.coll <- lm(Incorrect_success ~ Gender + Age_cent + Frequency + Culture, data=spontool)
vif(Incorrect_success.noslope.coll)
source("Diagnostic_fcns.r")
ranef.diagn.plot(Incorrect_success.noslope)
source("diagnostic_fcns.r")
overdisp.test(Incorrect_success.noslope)
source("glmm_stability.r")
bin.stab=glmm.model.stab(model.res=Incorrect_success.noslope)
bin.stab$detailed$warnings
bin.stab$summary

#DHARMa
simulationOutput_Incorrect_success.noslope <- simulateResiduals(fittedModel = Incorrect_success.noslope, n = 250)
hist(simulationOutput_Incorrect_success.noslope$scaledResiduals) # is more or less uniform
plot(simulationOutput_Incorrect_success.noslope)

#uniformity, outliers, dispersion
testResiduals(simulationOutput_Incorrect_success.noslope)
testDispersion(simulationOutput_Incorrect_success.noslope)


