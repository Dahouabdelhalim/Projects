####Winter feeding and breeding modeling script

#These are the conservative analyses presented in the supplemental materials in table S6.
#It compares breeding performance of birds that never acessed feeders to those that completed at least 20 trials
#during cognitive testing

#Title: Long-term winter food supplementation shows no significant impact on reproductive performance in 
#Mountain Chickadees (Poecile gambeli) in the Sierra Nevada Mountains
#Journal: Ornithology

#Contact Joseph Welklin and Ben Sonnenberg with questions:
#jwelklin@gmail.com
#benjamin.r.sonnenberg@gmail.com

#Script structure:
#Basic plots then models for each variable - First egg date, Clutch size, Brood size, Mean mass, and Mass CV.
#Fancy raincloud plots at the end of the script for each variable. 


library(here)
library(dplyr)
library(reshape2)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(climwin)
library(car)
library(ggplot2)
library(cowplot)
library(emmeans)
library(ggdist)
library(tibble)
library(plotrix)
library(performance)

####Load data####

##Raw data
bf = read.csv(here::here("Breeding_Feeding_raw_040122_pub.csv"))

##Filter for birds that did not visit at all, or tested
bf = bf %>% filter(M_F.detected==M_F.tested)
#Down to 278 pairs from 356

##Get factors for modeling
bf = bf %>% mutate(YEAR = as.factor(YEAR)) %>% mutate(SEASON = as.factor(SEASON)) %>%
  mutate(ELEVATION = as.factor(ELEVATION)) %>%
  mutate(NEST.TYPE = as.factor(NEST.TYPE)) %>% mutate(NEST.ID = as.factor(NEST.ID)) %>%
  mutate(M.ID = as.factor(M.ID)) %>%
  mutate(F.ID = as.factor(F.ID)) %>%
  mutate(M_F.detected = as.factor(M_F.detected)) %>%
  mutate(M_F.tested = as.factor(M_F.tested))


#____________________________________________



####Set up to model first egg dates####

##Set up full dataframe - both high and low - select only initial nests
bf.fe = bf %>% filter(!is.na(J.FIRST.EGG)) %>% filter(NEST.TYPE=="INITIAL") #Both

###Remove outlier super late nests
bf.fe = bf.fe %>% filter(!NEST.ID %in% c("N_392","N_635"))

##High
bf.feH = bf.fe %>% filter(ELEVATION=="H") #High

##Low
bf.feL = bf.fe %>% filter(ELEVATION=="L") #Low


###Plot data

##By year and elevation
ggplot(data=bf.fe,aes(x=YEAR,y=J.FIRST.EGG)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() + 
  facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By age and elevation
ggplot(data=bf.fe,aes(x=M_F.AGE,y=J.FIRST.EGG)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() + 
  facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By tested and elevation
ggplot(data=bf.fe,aes(x=M_F.tested,y=J.FIRST.EGG)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() + 
  facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))



####Model first egg dates - parents detected####

###High
bf.feH.dm1 = lmer(J.FIRST.EGG ~ M_F.tested + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feH)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.feH.dm1, plot = T)
#Look good

##Model output
Anova(bf.feH.dm1, test.statistic = "Chisq")
#Year is important. No effect of tested or age.
summary(bf.feH.dm1)


##R2 values
#Full model 
r2(bf.feH.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#59.39

#Detected
bf.feH.dm1.nod = lmer(J.FIRST.EGG ~ YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feH)
-(r2(bf.feH.dm1)$R2_marginal - r2(bf.feH.dm1.nod)$R2_marginal) * 100
#-0.83

#Year
bf.feH.dm1.noy = lmer(J.FIRST.EGG ~ M_F.detected + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feH)
-(r2(bf.feH.dm1)$R2_marginal - r2(bf.feH.dm1.noy)$R2_marginal) * 100
#-55.97

#Age
bf.feH.dm1.noa = lmer(J.FIRST.EGG ~ M_F.detected + YEAR + (1|F.ID) + (1|M.ID),data=bf.feH)
-(r2(bf.feH.dm1)$R2_marginal - r2(bf.feH.dm1.noa)$R2_marginal) * 100
#-0.06



###Low
bf.feL.dm1 = lmer(J.FIRST.EGG ~ M_F.tested + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feL)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.feL.dm1, plot = T)
#Look decent

##Model output
Anova(bf.feL.dm1, test.statistic = "Chisq")
#Year and tested are important. Age not important
summary(bf.feL.dm1)
feL.dm1p = emmeans(bf.feH.dm1,specs="M_F.tested",adjust="tukey")
pairs(feL.dm1p)
#No differences when run comparisons

##R2 values
#Full model
r2(bf.feL.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#18.85

#Detected
bf.feL.dm1.nod = lmer(J.FIRST.EGG ~ YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feL)
-(r2(bf.feL.dm1)$R2_marginal - r2(bf.feL.dm1.nod)$R2_marginal) * 100
#-3.65

#Year
bf.feL.dm1.noy = lmer(J.FIRST.EGG ~ M_F.detected + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feL)
-(r2(bf.feL.dm1)$R2_marginal - r2(bf.feL.dm1.noy)$R2_marginal) * 100
#-7.81

#Age
bf.feL.dm1.noa = lmer(J.FIRST.EGG ~ M_F.detected + YEAR + (1|F.ID) + (1|M.ID),data=bf.feL)
-(r2(bf.feL.dm1)$R2_marginal - r2(bf.feL.dm1.noa)$R2_marginal) * 100
#-1.61



#____________________________________________



####Set up to model clutch size####

##Set up dataframes - high, low, both - Select only initial nests
bf.cs = bf %>% filter(!is.na(CLUTCH)) %>% filter(NEST.TYPE=="INITIAL") #Both

#High
bf.csH = bf.cs %>% filter(ELEVATION=="H") #High

##Low
bf.csL = bf.cs %>% filter(ELEVATION=="L") #Low


###Plot data
##By year
ggplot(data=bf.cs,aes(x=YEAR,y=CLUTCH)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By age and elevation
ggplot(data=bf.cs,aes(x=M_F.AGE,y=CLUTCH)) + geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By tested and elevation
ggplot(data=bf.cs,aes(x=M_F.tested,y=CLUTCH)) + geom_boxplot(outlier.alpha = 0) +
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By first egg date and elevation
ggplot(data=bf.cs,aes(x=J.FIRST.EGG,y=CLUTCH)) + geom_smooth(method="lm",color="black") +
  geom_point(aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))



####Model clutch size - parents detected####

###High 
bf.csH.dm1 = glmmTMB(CLUTCH ~ M_F.tested + YEAR + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")

##Check residuals of model
simulationOutput <- simulateResiduals(bf.csH.dm1, plot = T)
testDispersion(simulationOutput)
#Looks good

##Model output
Anova(bf.csH.dm1,test.statistic = "Chisq")
#Year, and first egg date are important. No effect of age or tested
summary(bf.csH.dm1)


##R2 values
#Full model
r2(bf.csH.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#34.72

#Detected
bf.csH.dm1.nod = glmmTMB(CLUTCH ~ YEAR + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")
-(r2(bf.csH.dm1)$R2_marginal - r2(bf.csH.dm1.nod)$R2_marginal) * 100
#-4.00

#Year
bf.csH.dm1.noy = glmmTMB(CLUTCH ~ M_F.detected + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")
-(r2(bf.csH.dm1)$R2_marginal - r2(bf.csH.dm1.noy)$R2_marginal) * 100
#-4.33

#Age
bf.csH.dm1.noa = glmmTMB(CLUTCH ~ M_F.detected + YEAR + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")
-(r2(bf.csH.dm1)$R2_marginal - r2(bf.csH.dm1.noa)$R2_marginal) * 100
#-2.18

#First egg
bf.csH.dm1.noj = glmmTMB(CLUTCH ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")
-(r2(bf.csH.dm1)$R2_marginal - r2(bf.csH.dm1.noj)$R2_marginal) * 100
#-3.21



###Low
bf.csL.dm1 = glmmTMB(CLUTCH ~ M_F.tested + YEAR + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")

##Check residuals of model
simulationOutput <- simulateResiduals(bf.csL.dm1, plot = T)
testDispersion(simulationOutput)
#Looks good

##Model output
Anova(bf.csL.dm1,test.statistic = "Chisq")
#Year, and tested are important.  No effect of age or first egg date.
summary(bf.csL.dm1)
csL.dm1p = emmeans(bf.csL.dm1,specs="M_F.tested",adjust="tukey")
pairs(csL.dm1p)
#M_F - no_F p=0.005
#M-no - no_F p=0.008
#no_F - no_no p=0.006
#Similar results as to full dataset at high - extremely small sample size for no_F


##R2 values
#Full model
r2(bf.csL.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#29.62

#Detected
bf.csL.dm1.nod = glmmTMB(CLUTCH ~ YEAR + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")
-(r2(bf.csL.dm1)$R2_marginal - r2(bf.csL.dm1.nod)$R2_marginal) * 100
#-11.72

#Year
bf.csL.dm1.noy = glmmTMB(CLUTCH ~ M_F.detected + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")
-(r2(bf.csL.dm1)$R2_marginal - r2(bf.csL.dm1.noy)$R2_marginal) * 100
#-2.98

#Age
bf.csL.dm1.noa = glmmTMB(CLUTCH ~ M_F.detected + YEAR + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")
-(r2(bf.csL.dm1)$R2_marginal - r2(bf.csL.dm1.noa)$R2_marginal) * 100
#-0.32

#First egg
bf.csL.dm1.noj = glmmTMB(CLUTCH ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")
-(r2(bf.csL.dm1)$R2_marginal - r2(bf.csL.dm1.noj)$R2_marginal) * 100
#-1.06



#____________________________________________



####Set up to model brood size####

##Set up dataframes - high, low, both - Select only initial nests
bf.bs = bf %>% filter(!is.na(BROOD)) %>% filter(NEST.TYPE=="INITIAL") #Both

#High
bf.bsH = bf.bs %>% filter(ELEVATION=="H") #High

##Low
bf.bsL = bf.bs %>% filter(ELEVATION=="L") #Low


###Plot data
##By year
ggplot(data=bf.bs,aes(x=YEAR,y=BROOD)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By age and elevation
ggplot(data=bf.bs,aes(x=M_F.AGE,y=BROOD)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By detected and elevation
ggplot(data=bf.bs,aes(x=M_F.tested,y=BROOD)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))



####Model brood size - parents detected####

###High 
bf.bsH.dm1 = glmmTMB(BROOD ~ M_F.tested + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsH,family="genpois")

##Check residuals of model
simulationOutput <- simulateResiduals(bf.bsH.dm1, plot = T)
testDispersion(simulationOutput)
#Looks good

##Model output
Anova(bf.bsH.dm1,test.statistic = "Chisq")
#Only year is important. No effect of detected or age.
summary(bf.bsH.dm1)


##R2 values
#Full model
r2(bf.bsH.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#31.72

#Detected
bf.bsH.dm1.nod = glmmTMB(BROOD ~ YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsH,family="genpois")
-(r2(bf.bsH.dm1)$R2_marginal - r2(bf.bsH.dm1.nod)$R2_marginal) * 100
#-6.48

#Year
bf.bsH.dm1.noy = glmmTMB(BROOD ~ M_F.detected + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsH,family="genpois")
-(r2(bf.bsH.dm1)$R2_marginal - r2(bf.bsH.dm1.noy)$R2_marginal) * 100
#-24.42

#Age
bf.bsH.dm1.noa = glmmTMB(BROOD ~ M_F.detected + YEAR + (1|F.ID) + (1|M.ID),data=bf.bsH,family="genpois")
-(r2(bf.bsH.dm1)$R2_marginal - r2(bf.bsH.dm1.noa)$R2_marginal) * 100
#-1.78



###Low
bf.bsL.dm1 = glmmTMB(BROOD ~ M_F.tested + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsL,family="genpois")

##Check residuals of model
simulationOutput <- simulateResiduals(bf.bsL.dm1, plot = T)
testDispersion(simulationOutput)
#Looks good

##Model output
Anova(bf.bsL.dm1,test.statistic = "Chisq")
#Nothing significant
summary(bf.bsL.dm1)


##R2 values
#Full model
r2(bf.bsL.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#7.92

#Detected
bf.bsL.dm1.nod = glmmTMB(BROOD ~ YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsL,family="genpois")
-(r2(bf.bsL.dm1)$R2_marginal - r2(bf.bsL.dm1.nod)$R2_marginal) * 100
#-2.98

#Year
bf.bsL.dm1.noy = glmmTMB(BROOD ~ M_F.detected + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsL,family="genpois")
-(r2(bf.bsL.dm1)$R2_marginal - r2(bf.bsL.dm1.noy)$R2_marginal) * 100
#-4.35

#Age
bf.bsL.dm1.noa = glmmTMB(BROOD ~ M_F.detected + YEAR + (1|F.ID) + (1|M.ID),data=bf.bsL,family="genpois")
-(r2(bf.bsL.dm1)$R2_marginal - r2(bf.bsL.dm1.noa)$R2_marginal) * 100
#-1.69



#____________________________________________



####Set up to model mean mass####

##Set up dataframes - high, low, both - Select only initial nests
bf.mm = bf %>% filter(!is.na(MEANMASS)) %>% filter(NEST.TYPE=="INITIAL")

##Remove nests banded on day 18 - normally banded day 16
bf.mm = bf.mm %>% filter(is.na(NOTES))

#High
bf.mmH = bf.mm %>% filter(ELEVATION=="H")

#Low
bf.mmL = bf.mm %>% filter(ELEVATION=="L")


###Plot data
##By year
ggplot(data=bf.mm,aes(x=YEAR,y=MEANMASS)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(8,9,10,11,12,13,14,15,16),limits=c(8,16)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By age and elevation
ggplot(data=bf.mm,aes(x=M_F.AGE,y=MEANMASS)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(8,9,10,11,12,13,14,15,16),limits=c(8,16)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By detected and elevation
ggplot(data=bf.mm,aes(x=M_F.tested,y=MEANMASS)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(8,9,10,11,12,13,14,15,16),limits=c(8,16)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))

##By brood size and elevation
ggplot(data=bf.mm,aes(x=BROOD,y=MEANMASS)) + geom_smooth(method="lm",color="black") +
  geom_point(aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(8,9,10,11,12,13,14,15,16),limits=c(8,16)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1")) 



####Model mean mass - parents detected####

###High
bf.mmH.dm1 = lmer(MEANMASS ~ M_F.tested + YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmH)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.mmH.dm1, plot = T)
#Look good

##Model output
Anova(bf.mmH.dm1,test.statistic = "Chisq")
#Year and is important. No effect of detected or brood size or age.
summary(bf.mmH.dm1)


##R2 values
#Full model
r2(bf.mmH.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#18.21

#Detected
bf.mmH.dm1.nod = lmer(MEANMASS ~ YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmH)
-(r2(bf.mmH.dm1)$R2_marginal - r2(bf.mmH.dm1.nod)$R2_marginal) * 100
#-1.03

#Year
bf.mmH.dm1.noy = lmer(MEANMASS ~ M_F.detected + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmH)
-(r2(bf.mmH.dm1)$R2_marginal - r2(bf.mmH.dm1.noy)$R2_marginal) * 100
#-12.93

#Age
bf.mmH.dm1.noa = lmer(MEANMASS ~ M_F.detected + YEAR + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmH)
-(r2(bf.mmH.dm1)$R2_marginal - r2(bf.mmH.dm1.noa)$R2_marginal) * 100
#-2.61

#Brood
bf.mmH.dm1.nob = lmer(MEANMASS ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.mmH)
-(r2(bf.mmH.dm1)$R2_marginal - r2(bf.mmH.dm1.nob)$R2_marginal) * 100
#-1.10



###Low only
bf.mmL.dm1 = lmer(MEANMASS ~ M_F.tested + YEAR + M_F.AGE + BROOD +  (1|F.ID) + (1|M.ID),data=bf.mmL)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.mmL.dm1, plot = T)
#Looks ok...

##Model output
Anova(bf.mmL.dm1,test.statistic = "Chisq")
#Year and brood size are important. No effect of detected or age. 
summary(bf.mmL.dm1)


##R2 values
#Full model
r2(bf.mmL.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#39.22

#Detected
bf.mmL.dm1.nod = lmer(MEANMASS ~ YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmL)
-(r2(bf.mmL.dm1)$R2_marginal - r2(bf.mmL.dm1.nod)$R2_marginal) * 100
#-2.10

#Year
bf.mmL.dm1.noy = lmer(MEANMASS ~ M_F.detected + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmL)
-(r2(bf.mmL.dm1)$R2_marginal - r2(bf.mmL.dm1.noy)$R2_marginal) * 100
#-28.92

#Age
bf.mmL.dm1.noa = lmer(MEANMASS ~ M_F.detected + YEAR + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmL)
-(r2(bf.mmL.dm1)$R2_marginal - r2(bf.mmL.dm1.noa)$R2_marginal) * 100
#-0.54

#Brood
bf.mmL.dm1.nob = lmer(MEANMASS ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.mmL)
-(r2(bf.mmL.dm1)$R2_marginal - r2(bf.mmL.dm1.nob)$R2_marginal) * 100
#-2.88



#____________________________________________



####Set up to model CV of mass####

##Set up dataframes - high, low, both - Select only initial nests
bf.cv = bf %>% filter(!is.na(CV)) %>% filter(NEST.TYPE=="INITIAL")

##Remove nests banded on day 18 - normally banded day 16
bf.cv = bf.cv %>% filter(is.na(NOTES))

#High
bf.cvH = bf.cv %>% filter(ELEVATION=="H")

#Low
bf.cvL = bf.cv %>% filter(ELEVATION=="L")


###Plot data
##By year
ggplot(data=bf.cv,aes(x=YEAR,y=CV)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  facet_wrap(facets="ELEVATION") + scale_color_manual(values=c("steelblue3","darkorange1"))

##By age and elevation
ggplot(data=bf.cv,aes(x=M_F.AGE,y=CV)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  facet_wrap(facets="ELEVATION") + scale_color_manual(values=c("steelblue3","darkorange1"))

##By detected and elevation
ggplot(data=bf.cv,aes(x=M_F.tested,y=CV)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  facet_wrap(facets="ELEVATION") + scale_color_manual(values=c("steelblue3","darkorange1"))

##By brood size and elevation
ggplot(data=bf.cv,aes(x=BROOD,y=CV)) + geom_smooth(method="lm",color="black") +
  geom_point(aes(color=ELEVATION)) + theme_cowplot() +
  facet_wrap(facets="ELEVATION") + scale_color_manual(values=c("steelblue3","darkorange1")) 



####Model CV of mass - parents detected####

###High
bf.cvH.dm1 = lmer(log(CV) ~ M_F.tested + YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvH)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.cvH.dm1, plot = T)
#Look good when CV is log transformed

##Model output
Anova(bf.cvH.dm1,test.statistic = "Chisq")
#Nothing is important
summary(bf.cvH.dm1)


##R2 values
#Full model
r2(bf.cvH.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#6.44

#Detected
bf.cvH.dm1.nod = lmer(log(CV) ~ YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvH)
-(r2(bf.cvH.dm1)$R2_marginal - r2(bf.cvH.dm1.nod)$R2_marginal) * 100
#-0.38

#Year
bf.cvH.dm1.noy = lmer(log(CV) ~ M_F.detected + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvH)
-(r2(bf.cvH.dm1)$R2_marginal - r2(bf.cvH.dm1.noy)$R2_marginal) * 100
#-1.99

#Age
bf.cvH.dm1.noa = lmer(log(CV) ~ M_F.detected + YEAR + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvH)
-(r2(bf.cvH.dm1)$R2_marginal - r2(bf.cvH.dm1.noa)$R2_marginal) * 100
#-2.19

#Brood
bf.cvH.dm1.nob = lmer(log(CV) ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.cvH)
-(r2(bf.cvH.dm1)$R2_marginal - r2(bf.cvH.dm1.nob)$R2_marginal) * 100
#-0.16



###Low only
bf.cvL.dm1 = lmer(log(CV) ~ M_F.tested + YEAR + M_F.AGE + BROOD +  (1|F.ID) + (1|M.ID),data=bf.cvL)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.cvL.dm1, plot = T)
#Look good when CV is log transformed

##Model output
Anova(bf.cvL.dm1,test.statistic = "Chisq")
#Year and brood size are important. No effect of detected or age. 
summary(bf.cvL.dm1)


##R2 values
#Full model
r2(bf.cvL.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#12.67

#Detected
bf.cvL.dm1.nod = lmer(log(CV) ~ YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvL)
-(r2(bf.cvL.dm1)$R2_marginal - r2(bf.cvL.dm1.nod)$R2_marginal) * 100
#-2.32

#Year
bf.cvL.dm1.noy = lmer(log(CV) ~ M_F.detected + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvL)
-(r2(bf.cvL.dm1)$R2_marginal - r2(bf.cvL.dm1.noy)$R2_marginal) * 100
#-7.93

#Age
bf.cvL.dm1.noa = lmer(log(CV) ~ M_F.detected + YEAR + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvL)
-(r2(bf.cvL.dm1)$R2_marginal - r2(bf.cvL.dm1.noa)$R2_marginal) * 100
#-0.42

#Brood
bf.cvL.dm1.nob = lmer(log(CV) ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.cvL)
-(r2(bf.cvL.dm1)$R2_marginal - r2(bf.cvL.dm1.nob)$R2_marginal) * 100
#-3.29




