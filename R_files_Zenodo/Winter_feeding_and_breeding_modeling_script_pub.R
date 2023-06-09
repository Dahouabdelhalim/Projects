####Winter feeding and breeding modeling script - conservative

##These are the main analyses presented in the manuscript

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

##Get factors for modeling
bf = bf %>% mutate(YEAR = as.factor(YEAR)) %>% mutate(SEASON = as.factor(SEASON)) %>%
  mutate(ELEVATION = as.factor(ELEVATION)) %>%
  mutate(NEST.TYPE = as.factor(NEST.TYPE)) %>% mutate(NEST.ID = as.factor(NEST.ID)) %>%
  mutate(M.ID = as.factor(M.ID)) %>%
  mutate(F.ID = as.factor(F.ID)) %>% mutate(M_F.detected = as.factor(M_F.detected)) %>%
  mutate(M_F.tested = as.factor(M_F.tested))


#____________________________________________



####Set up to model first egg dates####

##Set up full dataframe - both high and low - select only initial nests
bf.fe = bf %>% filter(!is.na(J.FIRST.EGG)) %>% filter(NEST.TYPE=="INITIAL") #Both

###Remove outlier super late nests - these are likely re-nests
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

##By detected and elevation
ggplot(data=bf.fe,aes(x=M_F.detected,y=J.FIRST.EGG)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() + 
  facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))



####Model first egg dates - parents detected####

###High
bf.feH.dm1 = lmer(J.FIRST.EGG ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feH)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.feH.dm1, plot = T)
#Look good

##Model output
Anova(bf.feH.dm1, test.statistic = "Chisq")
#Year is important. No effect of detected or age.
summary(bf.feH.dm1)


##R2 values
#Full model 
r2(bf.feH.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#59.28

#Detected
bf.feH.dm1.nod = lmer(J.FIRST.EGG ~ YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feH)
-(r2(bf.feH.dm1)$R2_marginal - r2(bf.feH.dm1.nod)$R2_marginal) * 100
#-0.28

#Year
bf.feH.dm1.noy = lmer(J.FIRST.EGG ~ M_F.detected + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feH)
-(r2(bf.feH.dm1)$R2_marginal - r2(bf.feH.dm1.noy)$R2_marginal) * 100
#-57.28

#Age
bf.feH.dm1.noa = lmer(J.FIRST.EGG ~ M_F.detected + YEAR + (1|F.ID) + (1|M.ID),data=bf.feH)
-(r2(bf.feH.dm1)$R2_marginal - r2(bf.feH.dm1.noa)$R2_marginal) * 100
#-0.48



###Low
bf.feL.dm1 = lmer(J.FIRST.EGG ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feL)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.feL.dm1, plot = T)
#Look decent

##Model output
Anova(bf.feL.dm1, test.statistic = "Chisq")
#Year and age are important. No effect of detected.
summary(bf.feL.dm1)


##R2 values
#Full model
r2(bf.feL.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#17.87

#Detected
bf.feL.dm1.nod = lmer(J.FIRST.EGG ~ YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feL)
-(r2(bf.feL.dm1)$R2_marginal - r2(bf.feL.dm1.nod)$R2_marginal) * 100
#-0.44

#Year
bf.feL.dm1.noy = lmer(J.FIRST.EGG ~ M_F.detected + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.feL)
-(r2(bf.feL.dm1)$R2_marginal - r2(bf.feL.dm1.noy)$R2_marginal) * 100
#-9.10

#Age
bf.feL.dm1.noa = lmer(J.FIRST.EGG ~ M_F.detected + YEAR + (1|F.ID) + (1|M.ID),data=bf.feL)
-(r2(bf.feL.dm1)$R2_marginal - r2(bf.feL.dm1.noa)$R2_marginal) * 100
#-1.27



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

##By detected and elevation
ggplot(data=bf.cs,aes(x=M_F.detected,y=CLUTCH)) + geom_boxplot(outlier.alpha = 0) +
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
bf.csH.dm1 = glmmTMB(CLUTCH ~ M_F.detected + YEAR + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")

##Check residuals of model
simulationOutput <- simulateResiduals(bf.csH.dm1, plot = T)
testDispersion(simulationOutput)
#Looks good

##Model output
Anova(bf.csH.dm1,test.statistic = "Chisq")
#Detected, year, and first egg date are important. No effect of age. 
summary(bf.csH.dm1)
#Some differences, check with emmeans and tukey test
csH.dm1p = emmeans(bf.csH.dm1,specs="M_F.detected",adjust="tukey")
pairs(csH.dm1p)
#Only difference is M_F - no_F p=0.0351
#no_no vs no_F almost different p=0.056, with no_no being greater clutch size


##R2 values
#Full model
r2(bf.csH.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#32.63
  
#Detected
bf.csH.dm1.nod = glmmTMB(CLUTCH ~ YEAR + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")
-(r2(bf.csH.dm1)$R2_marginal - r2(bf.csH.dm1.nod)$R2_marginal) * 100
#-4.17

#Year
bf.csH.dm1.noy = glmmTMB(CLUTCH ~ M_F.detected + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")
-(r2(bf.csH.dm1)$R2_marginal - r2(bf.csH.dm1.noy)$R2_marginal) * 100
#-10.62

#Age
bf.csH.dm1.noa = glmmTMB(CLUTCH ~ M_F.detected + YEAR + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")
-(r2(bf.csH.dm1)$R2_marginal - r2(bf.csH.dm1.noa)$R2_marginal) * 100
#-1.52

#First egg
bf.csH.dm1.noj = glmmTMB(CLUTCH ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.csH,family="genpois")
-(r2(bf.csH.dm1)$R2_marginal - r2(bf.csH.dm1.noj)$R2_marginal) * 100
#-4.26



###Low
bf.csL.dm1 = glmmTMB(CLUTCH ~ M_F.detected + YEAR + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")

##Check residuals of model
simulationOutput <- simulateResiduals(bf.csL.dm1, plot = T)
testDispersion(simulationOutput)
#Looks good

##Model output
Anova(bf.csL.dm1,test.statistic = "Chisq")
#Year and first egg date are important. Detected = 0.055. No effect of age.
summary(bf.csL.dm1)


##R2 values
#Full model
r2(bf.csL.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#25.67

#Detected
bf.csL.dm1.nod = glmmTMB(CLUTCH ~ YEAR + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")
-(r2(bf.csL.dm1)$R2_marginal - r2(bf.csL.dm1.nod)$R2_marginal) * 100
#-4.61

#Year
bf.csL.dm1.noy = glmmTMB(CLUTCH ~ M_F.detected + M_F.AGE + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")
-(r2(bf.csL.dm1)$R2_marginal - r2(bf.csL.dm1.noy)$R2_marginal) * 100
#-5.18

#Age
bf.csL.dm1.noa = glmmTMB(CLUTCH ~ M_F.detected + YEAR + scale(J.FIRST.EGG) + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")
-(r2(bf.csL.dm1)$R2_marginal - r2(bf.csL.dm1.noa)$R2_marginal) * 100
#-1.21

#First egg
bf.csL.dm1.noj = glmmTMB(CLUTCH ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.csL,family="genpois")
-(r2(bf.csL.dm1)$R2_marginal - r2(bf.csL.dm1.noj)$R2_marginal) * 100
#-3.86



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
ggplot(data=bf.bs,aes(x=M_F.detected,y=BROOD)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="ELEVATION") +
  scale_color_manual(values=c("steelblue3","darkorange1"))



####Model brood size - parents detected####

###High 
bf.bsH.dm1 = glmmTMB(BROOD ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsH,family="genpois")

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
#22.83

#Detected
bf.bsH.dm1.nod = glmmTMB(BROOD ~ YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsH,family="genpois")
-(r2(bf.bsH.dm1)$R2_marginal - r2(bf.bsH.dm1.nod)$R2_marginal) * 100
#-2.51

#Year
bf.bsH.dm1.noy = glmmTMB(BROOD ~ M_F.detected + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsH,family="genpois")
-(r2(bf.bsH.dm1)$R2_marginal - r2(bf.bsH.dm1.noy)$R2_marginal) * 100
#-16.35

#Age
bf.bsH.dm1.noa = glmmTMB(BROOD ~ M_F.detected + YEAR + (1|F.ID) + (1|M.ID),data=bf.bsH,family="genpois")
-(r2(bf.bsH.dm1)$R2_marginal - r2(bf.bsH.dm1.noa)$R2_marginal) * 100
#-4.27



###Low
bf.bsL.dm1 = glmmTMB(BROOD ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsL,family="genpois")

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
#9.37

#Detected
bf.bsL.dm1.nod = glmmTMB(BROOD ~ YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsL,family="genpois")
-(r2(bf.bsL.dm1)$R2_marginal - r2(bf.bsL.dm1.nod)$R2_marginal) * 100
#-2.46

#Year
bf.bsL.dm1.noy = glmmTMB(BROOD ~ M_F.detected + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.bsL,family="genpois")
-(r2(bf.bsL.dm1)$R2_marginal - r2(bf.bsL.dm1.noy)$R2_marginal) * 100
#-4.66

#Age
bf.bsL.dm1.noa = glmmTMB(BROOD ~ M_F.detected + YEAR + (1|F.ID) + (1|M.ID),data=bf.bsL,family="genpois")
-(r2(bf.bsL.dm1)$R2_marginal - r2(bf.bsL.dm1.noa)$R2_marginal) * 100
#-2.76



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
ggplot(data=bf.mm,aes(x=M_F.detected,y=MEANMASS)) + geom_boxplot(outlier.alpha = 0) + 
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
bf.mmH.dm1 = lmer(MEANMASS ~ M_F.detected + YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmH)

##Check residuals of model
simulationOutput <- simulateResiduals(bf.mmH.dm1, plot = T)
#Look ok...

##Model output
Anova(bf.mmH.dm1,test.statistic = "Chisq")
#Year and age are important. No effect of detected or brood size.
summary(bf.mmH.dm1)


##R2 values
#Full model
r2(bf.mmH.dm1)$R2_marginal * 100 #Marginal is variance explained by the fixed effects, conditional is the full model
#21.78

#Detected
bf.mmH.dm1.nod = lmer(MEANMASS ~ YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmH)
-(r2(bf.mmH.dm1)$R2_marginal - r2(bf.mmH.dm1.nod)$R2_marginal) * 100
#-2.38

#Year
bf.mmH.dm1.noy = lmer(MEANMASS ~ M_F.detected + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmH)
-(r2(bf.mmH.dm1)$R2_marginal - r2(bf.mmH.dm1.noy)$R2_marginal) * 100
#-13.08

#Age
bf.mmH.dm1.noa = lmer(MEANMASS ~ M_F.detected + YEAR + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmH)
-(r2(bf.mmH.dm1)$R2_marginal - r2(bf.mmH.dm1.noa)$R2_marginal) * 100
#-4.39

#Brood
bf.mmH.dm1.nob = lmer(MEANMASS ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.mmH)
-(r2(bf.mmH.dm1)$R2_marginal - r2(bf.mmH.dm1.nob)$R2_marginal) * 100
#-0.63



###Low only
bf.mmL.dm1 = lmer(MEANMASS ~ M_F.detected + YEAR + M_F.AGE + BROOD +  (1|F.ID) + (1|M.ID),data=bf.mmL)

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
#36.19

#Detected
bf.mmL.dm1.nod = lmer(MEANMASS ~ YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmL)
-(r2(bf.mmL.dm1)$R2_marginal - r2(bf.mmL.dm1.nod)$R2_marginal) * 100
#-1.20

#Year
bf.mmL.dm1.noy = lmer(MEANMASS ~ M_F.detected + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmL)
-(r2(bf.mmL.dm1)$R2_marginal - r2(bf.mmL.dm1.noy)$R2_marginal) * 100
#-28.03

#Age
bf.mmL.dm1.noa = lmer(MEANMASS ~ M_F.detected + YEAR + BROOD + (1|F.ID) + (1|M.ID),data=bf.mmL)
-(r2(bf.mmL.dm1)$R2_marginal - r2(bf.mmL.dm1.noa)$R2_marginal) * 100
#-0.70

#Brood
bf.mmL.dm1.nob = lmer(MEANMASS ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.mmL)
-(r2(bf.mmL.dm1)$R2_marginal - r2(bf.mmL.dm1.nob)$R2_marginal) * 100
#-1.99


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
ggplot(data=bf.cv,aes(x=M_F.detected,y=CV)) + geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = position_jitter(width=0.2,height=0),aes(color=ELEVATION)) + theme_cowplot() +
  facet_wrap(facets="ELEVATION") + scale_color_manual(values=c("steelblue3","darkorange1"))

##By brood size and elevation
ggplot(data=bf.cv,aes(x=BROOD,y=CV)) + geom_smooth(method="lm",color="black") +
  geom_point(aes(color=ELEVATION)) + theme_cowplot() +
  facet_wrap(facets="ELEVATION") + scale_color_manual(values=c("steelblue3","darkorange1")) 



####Model CV of mass - parents detected####

###High
bf.cvH.dm1 = lmer(log(CV) ~ M_F.detected + YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvH)

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
#6.70

#Detected
bf.cvH.dm1.nod = lmer(log(CV) ~ YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvH)
-(r2(bf.cvH.dm1)$R2_marginal - r2(bf.cvH.dm1.nod)$R2_marginal) * 100
#-0.25

#Year
bf.cvH.dm1.noy = lmer(log(CV) ~ M_F.detected + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvH)
-(r2(bf.cvH.dm1)$R2_marginal - r2(bf.cvH.dm1.noy)$R2_marginal) * 100
#-2.35

#Age
bf.cvH.dm1.noa = lmer(log(CV) ~ M_F.detected + YEAR + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvH)
-(r2(bf.cvH.dm1)$R2_marginal - r2(bf.cvH.dm1.noa)$R2_marginal) * 100
#-1.07

#Brood
bf.cvH.dm1.nob = lmer(log(CV) ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.cvH)
-(r2(bf.cvH.dm1)$R2_marginal - r2(bf.cvH.dm1.nob)$R2_marginal) * 100
#-0.58



###Low only
bf.cvL.dm1 = lmer(log(CV) ~ M_F.detected + YEAR + M_F.AGE + BROOD +  (1|F.ID) + (1|M.ID),data=bf.cvL)

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
#14.71

#Detected
bf.cvL.dm1.nod = lmer(log(CV) ~ YEAR + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvL)
-(r2(bf.cvL.dm1)$R2_marginal - r2(bf.cvL.dm1.nod)$R2_marginal) * 100
#-0.28

#Year
bf.cvL.dm1.noy = lmer(log(CV) ~ M_F.detected + M_F.AGE + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvL)
-(r2(bf.cvL.dm1)$R2_marginal - r2(bf.cvL.dm1.noy)$R2_marginal) * 100
#-11.18

#Age
bf.cvL.dm1.noa = lmer(log(CV) ~ M_F.detected + YEAR + BROOD + (1|F.ID) + (1|M.ID),data=bf.cvL)
-(r2(bf.cvL.dm1)$R2_marginal - r2(bf.cvL.dm1.noa)$R2_marginal) * 100
#-1.09

#Brood
bf.cvL.dm1.nob = lmer(log(CV) ~ M_F.detected + YEAR + M_F.AGE + (1|F.ID) + (1|M.ID),data=bf.cvL)
-(r2(bf.cvL.dm1)$R2_marginal - r2(bf.cvL.dm1.nob)$R2_marginal) * 100
#-2.17



#____________________________________________



####Plot first egg date fancy####

#Note: in stat_dots() binwidth=unit(0.015,"npc") means binwidths (width of y-axis bins) are 1.5% of visual area on graph
#dotsize=unit(0.8,"npc") means dots are 80% of the binwidth. See ?stat_dots() for more info. 

##Get facet titles for plotting
bf.fe = bf.fe %>% mutate(elevation.title=ifelse(ELEVATION=="H","High Elevation","Low Elevation"))


###Year and elevation
##Get groups into list for plotting
bf.fe.ylist = bf.fe %>% group_by(ELEVATION,YEAR) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
bf.fe.ysize = bf.fe %>% group_by(ELEVATION,YEAR,J.FIRST.EGG) %>% summarise(count=n()) %>% 
  group_by(ELEVATION,YEAR) %>% summarise(max=max(count)/10)

##Plot year
ggplot(data=bf.fe,aes(x=YEAR,y=J.FIRST.EGG)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.fe.ylist[[1]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[1]) +
  stat_halfeye(data=bf.fe.ylist[[2]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[2]) +
  stat_halfeye(data=bf.fe.ylist[[3]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[3]) +
  stat_halfeye(data=bf.fe.ylist[[4]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[4]) +
  stat_halfeye(data=bf.fe.ylist[[5]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[5]) +
  stat_halfeye(data=bf.fe.ylist[[6]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[6]) +
  stat_halfeye(data=bf.fe.ylist[[7]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[7]) +
  stat_halfeye(data=bf.fe.ylist[[8]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[8]) +
  stat_halfeye(data=bf.fe.ylist[[9]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[9]) +
  stat_halfeye(data=bf.fe.ylist[[10]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[10]) +
  stat_halfeye(data=bf.fe.ylist[[11]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[11]) +
  stat_halfeye(data=bf.fe.ylist[[12]],aes(x=YEAR,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.ysize$max[12]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + xlab("Year") + ylab("First egg date")


###Age and elevation
##Get groups into list for plotting
bf.fe.alist = bf.fe %>% group_by(ELEVATION,M_F.AGE) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
bf.fe.asize = bf.fe %>% group_by(ELEVATION,M_F.AGE,J.FIRST.EGG) %>% summarise(count=n()) %>% 
  group_by(ELEVATION,M_F.AGE) %>% summarise(max=max(count)/18)

##Plot age 
ggplot(data=bf.fe,aes(x=M_F.AGE,y=J.FIRST.EGG)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.fe.alist[[1]],aes(x=M_F.AGE,y=J.FIRST.EGG,fill=ELEVATION),adjust=1,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.asize$max[1]) +
  stat_halfeye(data=bf.fe.alist[[2]],aes(x=M_F.AGE,y=J.FIRST.EGG,fill=ELEVATION),adjust=1,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.asize$max[2]) +
  stat_halfeye(data=bf.fe.alist[[3]],aes(x=M_F.AGE,y=J.FIRST.EGG,fill=ELEVATION),adjust=1,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.asize$max[3]) +
  stat_halfeye(data=bf.fe.alist[[4]],aes(x=M_F.AGE,y=J.FIRST.EGG,fill=ELEVATION),adjust=1,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.asize$max[4]) +
  stat_halfeye(data=bf.fe.alist[[5]],aes(x=M_F.AGE,y=J.FIRST.EGG,fill=ELEVATION),adjust=1,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.asize$max[5]) +
  stat_halfeye(data=bf.fe.alist[[6]],aes(x=M_F.AGE,y=J.FIRST.EGG,fill=ELEVATION),adjust=1,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.asize$max[6]) +
  stat_halfeye(data=bf.fe.alist[[7]],aes(x=M_F.AGE,y=J.FIRST.EGG,fill=ELEVATION),adjust=1,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.asize$max[7]) +
  stat_halfeye(data=bf.fe.alist[[8]],aes(x=M_F.AGE,y=J.FIRST.EGG,fill=ELEVATION),adjust=1,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.asize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + xlab("Parents age") + ylab("First egg date") +
  scale_x_discrete(labels=c("Both\\nAdult","Adult M\\nJuv F","Juv M.\\nAdult F","Both\\nJuv"))


###Detected and elevation
##Get groups into list for plotting
bf.fe.dlist = bf.fe %>% group_by(ELEVATION,M_F.detected) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
bf.fe.dsize = bf.fe %>% group_by(ELEVATION,M_F.detected,J.FIRST.EGG) %>% summarise(count=n()) %>% 
  group_by(ELEVATION,M_F.detected) %>% summarise(max=max(count)/12)

##Plot detected 
ggplot(data=bf.fe,aes(x=M_F.detected,y=J.FIRST.EGG)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.3,aes(fill=ELEVATION),position = position_nudge(x=-0.2)) + 
  stat_halfeye(data=bf.fe.dlist[[1]],aes(x=M_F.detected,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.dsize$max[1]) +
  stat_halfeye(data=bf.fe.dlist[[2]],aes(x=M_F.detected,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.dsize$max[2]) +
  stat_halfeye(data=bf.fe.dlist[[3]],aes(x=M_F.detected,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.dsize$max[3]) +
  stat_halfeye(data=bf.fe.dlist[[4]],aes(x=M_F.detected,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.dsize$max[4]) +
  stat_halfeye(data=bf.fe.dlist[[5]],aes(x=M_F.detected,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.dsize$max[5]) +
  stat_halfeye(data=bf.fe.dlist[[6]],aes(x=M_F.detected,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.dsize$max[6]) +
  stat_halfeye(data=bf.fe.dlist[[7]],aes(x=M_F.detected,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.dsize$max[7]) +
  stat_halfeye(data=bf.fe.dlist[[8]],aes(x=M_F.detected,y=J.FIRST.EGG,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.fe.dsize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.02,"npc")) + xlab("Parents detected") + ylab("First egg date") +
  scale_x_discrete(labels=c("Both","Male\\nonly","Female\\nonly","Neither"))


#____________________________________________



####Plot clutch size fancy####

##Get facet titles for plotting
bf.cs = bf.cs %>% mutate(elevation.title=ifelse(ELEVATION=="H","High Elevation","Low Elevation"))

##Manually jitter clutch size data for plotting
bf.cs = bf.cs %>% group_by(ELEVATION,YEAR,CLUTCH) %>% mutate(n=n()) %>% mutate(count=1:n()) %>%
  mutate(CLUTCH.plot.year=ifelse(count<=(n/2),CLUTCH-0.08,CLUTCH+0.08))
bf.cs = bf.cs %>% group_by(ELEVATION,M_F.AGE,CLUTCH) %>% mutate(n=n()) %>% mutate(count=1:n()) %>%
  mutate(CLUTCH.plot.age = ifelse(count<=(n/3),CLUTCH-0.15,ifelse(count>(n-(n/3)),CLUTCH+0.15,CLUTCH)))
bf.cs = bf.cs %>% group_by(ELEVATION,M_F.detected,CLUTCH) %>% mutate(n=n()) %>% mutate(count=1:n()) %>%
  mutate(CLUTCH.plot.det = ifelse(count<=(n/3),CLUTCH-0.2,ifelse(count>(n-(n/3)),CLUTCH+0.2,CLUTCH)))



###Year and elevation
##Get groups into list for plotting
bf.cs.ylist = bf.cs %>% group_by(ELEVATION,YEAR) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.cs.ysize = bf.cs %>% group_by(ELEVATION,YEAR) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(3,4,7,5,6,9,3,6,7,9,8,8)) %>% mutate(max=max1/11)

##Plot year - custom adjustment for 2017 High
ggplot(data=bf.cs,aes(x=YEAR,y=CLUTCH.plot.year)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.cs.ylist[[1]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[1]) +
  stat_halfeye(data=bf.cs.ylist[[2]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=20,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[2]) +
  stat_halfeye(data=bf.cs.ylist[[3]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[3]) +
  stat_halfeye(data=bf.cs.ylist[[4]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[4]) +
  stat_halfeye(data=bf.cs.ylist[[5]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[5]) +
  stat_halfeye(data=bf.cs.ylist[[6]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[6]) +
  stat_halfeye(data=bf.cs.ylist[[7]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[7]) +
  stat_halfeye(data=bf.cs.ylist[[8]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[8]) +
  stat_halfeye(data=bf.cs.ylist[[9]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[9]) +
  stat_halfeye(data=bf.cs.ylist[[10]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[10]) +
  stat_halfeye(data=bf.cs.ylist[[11]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[11]) +
  stat_halfeye(data=bf.cs.ylist[[12]],aes(x=YEAR,y=CLUTCH.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.ysize$max[12]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + xlab("Year") + ylab("Clutch size") +
  scale_y_continuous(breaks=c(2,4,6,8,10,12),limits = c(1.9,11.1))



###Age and elevation
##Get groups into list for plotting
bf.cs.alist = bf.cs %>% group_by(ELEVATION,M_F.AGE) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.cs.asize = bf.cs %>% group_by(ELEVATION,M_F.AGE) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(15,3,2,1,19,2,3,1)) %>% mutate(max=max1/18)

##Plot age
ggplot(data=bf.cs,aes(x=M_F.AGE,y=CLUTCH.plot.age)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.cs.alist[[1]],aes(x=M_F.AGE,y=CLUTCH.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.asize$max[1]) +
  stat_halfeye(data=bf.cs.alist[[2]],aes(x=M_F.AGE,y=CLUTCH.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.asize$max[2]) +
  stat_halfeye(data=bf.cs.alist[[3]],aes(x=M_F.AGE,y=CLUTCH.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.asize$max[3]) +
  stat_halfeye(data=bf.cs.alist[[4]],aes(x=M_F.AGE,y=CLUTCH.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.asize$max[4]) +
  stat_halfeye(data=bf.cs.alist[[5]],aes(x=M_F.AGE,y=CLUTCH.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.asize$max[5]) +
  stat_halfeye(data=bf.cs.alist[[6]],aes(x=M_F.AGE,y=CLUTCH.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.asize$max[6]) +
  stat_halfeye(data=bf.cs.alist[[7]],aes(x=M_F.AGE,y=CLUTCH.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.asize$max[7]) +
  stat_halfeye(data=bf.cs.alist[[8]],aes(x=M_F.AGE,y=CLUTCH.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.asize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + xlab("Parents age") + ylab("Clutch size") +
  scale_x_discrete(labels=c("Both\\nAdult","Adult M\\nJuv F","Juv M.\\nAdult F","Both\\nJuv"))



###Detected and elevation
##Get groups into list for plotting
bf.cs.dlist = bf.cs %>% group_by(ELEVATION,M_F.detected) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.cs.dsize = bf.cs %>% group_by(ELEVATION,M_F.detected) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(9,4,3,2,8,4,1,12)) %>% mutate(max=max1/13)

##Plot detected 
ggplot(data=bf.cs,aes(x=M_F.detected,y=CLUTCH.plot.det)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.3,aes(fill=ELEVATION),position = position_nudge(x=-0.2)) + 
  stat_halfeye(data=bf.cs.dlist[[1]],aes(x=M_F.detected,y=CLUTCH.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.dsize$max[1]) +
  stat_halfeye(data=bf.cs.dlist[[2]],aes(x=M_F.detected,y=CLUTCH.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.dsize$max[2]) +
  stat_halfeye(data=bf.cs.dlist[[3]],aes(x=M_F.detected,y=CLUTCH.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.dsize$max[3]) +
  stat_halfeye(data=bf.cs.dlist[[4]],aes(x=M_F.detected,y=CLUTCH.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.dsize$max[4]) +
  stat_halfeye(data=bf.cs.dlist[[5]],aes(x=M_F.detected,y=CLUTCH.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.dsize$max[5]) +
  stat_halfeye(data=bf.cs.dlist[[6]],aes(x=M_F.detected,y=CLUTCH.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.dsize$max[6]) +
  stat_halfeye(data=bf.cs.dlist[[7]],aes(x=M_F.detected,y=CLUTCH.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.dsize$max[7]) +
  stat_halfeye(data=bf.cs.dlist[[8]],aes(x=M_F.detected,y=CLUTCH.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cs.dsize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.02,"npc")) + xlab("Parents detected") + ylab("Clutch size") +
  scale_x_discrete(labels=c("Both","Male\\nonly","Female\\nonly","Neither"))



###First egg date and elevation
ggplot(data=bf.cs,aes(x=J.FIRST.EGG,y=CLUTCH)) + geom_smooth(method="lm",color="black") +
  geom_point(aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12),limits = c(0,12)) + facet_wrap(facets="elevation.title") +
  scale_color_manual(values=c("steelblue3","darkorange1")) + xlab("First egg date") + ylab("Clutch size") +
  theme(legend.position = "")


#____________________________________________



####Plot brood size fancy####

##Get facet titles for plotting
bf.bs = bf.bs %>% mutate(elevation.title=ifelse(ELEVATION=="H","High Elevation","Low Elevation"))

##Manually jitter brood size data for plotting
bf.bs = bf.bs %>% group_by(ELEVATION,YEAR,BROOD) %>% mutate(n=n()) %>% mutate(count=1:n()) %>%
  mutate(BROOD.plot.year=ifelse(count<=(n/2),BROOD-0.07,BROOD+0.07))
bf.bs = bf.bs %>% group_by(ELEVATION,M_F.AGE,BROOD) %>% mutate(n=n()) %>% mutate(count=1:n()) %>%
  mutate(BROOD.plot.age = ifelse(count<=(n/3),BROOD-0.15,ifelse(count>(n-(n/3)),BROOD+0.15,BROOD)))
bf.bs = bf.bs %>% group_by(ELEVATION,M_F.detected,BROOD) %>% mutate(n=n()) %>% mutate(count=1:n()) %>%
  mutate(BROOD.plot.det = ifelse(count<=(n/3),BROOD-0.2,ifelse(count>(n-(n/3)),BROOD+0.2,BROOD)))



###Year and elevation
##Get groups into list for plotting
bf.bs.ylist = bf.bs %>% group_by(ELEVATION,YEAR) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.bs.ysize = bf.bs %>% group_by(ELEVATION,YEAR) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(3,3,4,3,6,7,3,4,5,8,6,7)) %>% mutate(max=max1/12)

##Plot year
ggplot(data=bf.bs,aes(x=YEAR,y=BROOD.plot.year)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.bs.ylist[[1]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[1]) +
  stat_halfeye(data=bf.bs.ylist[[2]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[2]) +
  stat_halfeye(data=bf.bs.ylist[[3]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[3]) +
  stat_halfeye(data=bf.bs.ylist[[4]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[4]) +
  stat_halfeye(data=bf.bs.ylist[[5]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[5]) +
  stat_halfeye(data=bf.bs.ylist[[6]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[6]) +
  stat_halfeye(data=bf.bs.ylist[[7]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[7]) +
  stat_halfeye(data=bf.bs.ylist[[8]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[8]) +
  stat_halfeye(data=bf.bs.ylist[[9]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[9]) +
  stat_halfeye(data=bf.bs.ylist[[10]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[10]) +
  stat_halfeye(data=bf.bs.ylist[[11]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[11]) +
  stat_halfeye(data=bf.bs.ylist[[12]],aes(x=YEAR,y=BROOD.plot.year,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.ysize$max[12]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + xlab("Year") + ylab("Brood size") +
  scale_y_continuous(breaks=c(2,4,6,8),limits = c(0.8,9.2))



###Age and elevation
##Get groups into list for plotting
bf.bs.alist = bf.bs %>% group_by(ELEVATION,M_F.AGE) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.bs.asize = bf.bs %>% group_by(ELEVATION,M_F.AGE) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(11,2,2,1,15,2,3,1)) %>% mutate(max=max1/18)


##Plot age
ggplot(data=bf.bs,aes(x=M_F.AGE,y=BROOD.plot.age)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.bs.alist[[1]],aes(x=M_F.AGE,y=BROOD.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.asize$max[1]) +
  stat_halfeye(data=bf.bs.alist[[2]],aes(x=M_F.AGE,y=BROOD.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.asize$max[2]) +
  stat_halfeye(data=bf.bs.alist[[3]],aes(x=M_F.AGE,y=BROOD.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.asize$max[3]) +
  stat_halfeye(data=bf.bs.alist[[4]],aes(x=M_F.AGE,y=BROOD.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.asize$max[4]) +
  stat_halfeye(data=bf.bs.alist[[5]],aes(x=M_F.AGE,y=BROOD.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.asize$max[5]) +
  stat_halfeye(data=bf.bs.alist[[6]],aes(x=M_F.AGE,y=BROOD.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.asize$max[6]) +
  stat_halfeye(data=bf.bs.alist[[7]],aes(x=M_F.AGE,y=BROOD.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.asize$max[7]) +
  stat_halfeye(data=bf.bs.alist[[8]],aes(x=M_F.AGE,y=BROOD.plot.age,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.asize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + ylab("Brood size") +
  scale_y_continuous(breaks=c(2,4,6,8),limits = c(0.8,9.2)) + xlab("Parents age") +
  scale_x_discrete(labels=c("Both\\nAdult","Adult M\\nJuv F","Juv M.\\nAdult F","Both\\nJuv"))



###Detected and elevation
##Get groups into list for plotting
bf.bs.dlist = bf.bs %>% group_by(ELEVATION,M_F.detected) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.bs.dsize = bf.bs %>% group_by(ELEVATION,M_F.detected) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(8,3,2,2,7,2,1,10)) %>% mutate(max=max1/12)

##Plot detected 
ggplot(data=bf.bs,aes(x=M_F.detected,y=BROOD.plot.det)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.3,aes(fill=ELEVATION),position = position_nudge(x=-0.2)) + 
  stat_halfeye(data=bf.bs.dlist[[1]],aes(x=M_F.detected,y=BROOD.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.dsize$max[1]) +
  stat_halfeye(data=bf.bs.dlist[[2]],aes(x=M_F.detected,y=BROOD.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.dsize$max[2]) +
  stat_halfeye(data=bf.bs.dlist[[3]],aes(x=M_F.detected,y=BROOD.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.dsize$max[3]) +
  stat_halfeye(data=bf.bs.dlist[[4]],aes(x=M_F.detected,y=BROOD.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.dsize$max[4]) +
  stat_halfeye(data=bf.bs.dlist[[5]],aes(x=M_F.detected,y=BROOD.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.dsize$max[5]) +
  stat_halfeye(data=bf.bs.dlist[[6]],aes(x=M_F.detected,y=BROOD.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.dsize$max[6]) +
  stat_halfeye(data=bf.bs.dlist[[7]],aes(x=M_F.detected,y=BROOD.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.dsize$max[7]) +
  stat_halfeye(data=bf.bs.dlist[[8]],aes(x=M_F.detected,y=BROOD.plot.det,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.bs.dsize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.02,"npc")) + xlab("Parents detected") + ylab("Brood size") +
  scale_x_discrete(labels=c("Both","Male\\nonly","Female\\nonly","Neither")) +
  scale_y_continuous(breaks=c(2,4,6,8),limits = c(0.8,9.2)) 



#____________________________________________



####Plot mean mass fancy####

##Get facet titles for plotting
bf.mm = bf.mm %>% mutate(elevation.title=ifelse(ELEVATION=="H","High Elevation","Low Elevation"))


###Year and elevation
##Get groups into list for plotting
bf.mm.ylist = bf.mm %>% group_by(ELEVATION,YEAR) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.mm.ysize = bf.mm %>% group_by(ELEVATION,YEAR) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(3,2,5,4,4,6,3,3,4,6,3,5)) %>% mutate(max=max1/10)

##Plot year
ggplot(data=bf.mm,aes(x=YEAR,y=MEANMASS)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.mm.ylist[[1]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[1]) +
  stat_halfeye(data=bf.mm.ylist[[2]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[2]) +
  stat_halfeye(data=bf.mm.ylist[[3]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[3]) +
  stat_halfeye(data=bf.mm.ylist[[4]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[4]) +
  stat_halfeye(data=bf.mm.ylist[[5]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[5]) +
  stat_halfeye(data=bf.mm.ylist[[6]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[6]) +
  stat_halfeye(data=bf.mm.ylist[[7]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[7]) +
  stat_halfeye(data=bf.mm.ylist[[8]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[8]) +
  stat_halfeye(data=bf.mm.ylist[[9]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[9]) +
  stat_halfeye(data=bf.mm.ylist[[10]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[10]) +
  stat_halfeye(data=bf.mm.ylist[[11]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[11]) +
  stat_halfeye(data=bf.mm.ylist[[12]],aes(x=YEAR,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.ysize$max[12]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + xlab("Year") + ylab("Mean nestling mass")



###Age and elevation
##Get groups into list for plotting
bf.mm.alist = bf.mm %>% group_by(ELEVATION,M_F.AGE) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.mm.asize = bf.mm %>% group_by(ELEVATION,M_F.AGE) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(10,2,2,1,12,2,2,2)) %>% mutate(max=max1/17)

##Plot age
ggplot(data=bf.mm,aes(x=M_F.AGE,y=MEANMASS)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.mm.alist[[1]],aes(x=M_F.AGE,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.asize$max[1]) +
  stat_halfeye(data=bf.mm.alist[[2]],aes(x=M_F.AGE,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.asize$max[2]) +
  stat_halfeye(data=bf.mm.alist[[3]],aes(x=M_F.AGE,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.asize$max[3]) +
  stat_halfeye(data=bf.mm.alist[[4]],aes(x=M_F.AGE,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.asize$max[4]) +
  stat_halfeye(data=bf.mm.alist[[5]],aes(x=M_F.AGE,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.asize$max[5]) +
  stat_halfeye(data=bf.mm.alist[[6]],aes(x=M_F.AGE,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.asize$max[6]) +
  stat_halfeye(data=bf.mm.alist[[7]],aes(x=M_F.AGE,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.asize$max[7]) +
  stat_halfeye(data=bf.mm.alist[[8]],aes(x=M_F.AGE,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.asize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + 
  xlab("Male_Female Age") + ylab("Mean nestling mass") + xlab("Parents age") +
  scale_x_discrete(labels=c("Both\\nAdult","Adult M\\nJuv F","Juv M.\\nAdult F","Both\\nJuv"))



###Detected and elevation
##Get groups into list for plotting
bf.mm.dlist = bf.mm %>% group_by(ELEVATION,M_F.detected) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.mm.dsize = bf.mm %>% group_by(ELEVATION,M_F.AGE) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(8,5,3,3,8,3,3,11)) %>% mutate(max=max1/10)

##Plot detected 
ggplot(data=bf.mm,aes(x=M_F.detected,y=MEANMASS)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  stat_halfeye(data=bf.mm.dlist[[1]],aes(x=M_F.detected,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.dsize$max[1]) +
  stat_halfeye(data=bf.mm.dlist[[2]],aes(x=M_F.detected,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.dsize$max[2]) +
  stat_halfeye(data=bf.mm.dlist[[3]],aes(x=M_F.detected,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.dsize$max[3]) +
  stat_halfeye(data=bf.mm.dlist[[4]],aes(x=M_F.detected,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.dsize$max[4]) +
  stat_halfeye(data=bf.mm.dlist[[5]],aes(x=M_F.detected,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.dsize$max[5]) +
  stat_halfeye(data=bf.mm.dlist[[6]],aes(x=M_F.detected,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.dsize$max[6]) +
  stat_halfeye(data=bf.mm.dlist[[7]],aes(x=M_F.detected,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.dsize$max[7]) +
  stat_halfeye(data=bf.mm.dlist[[8]],aes(x=M_F.detected,y=MEANMASS,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.mm.dsize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.02,"npc")) + 
  geom_boxplot(outlier.alpha = 0,width=0.3,aes(fill=ELEVATION),position = position_nudge(x=-0.2)) + 
  xlab("Parents detected") + ylab("Mean nestling mass") +
  scale_x_discrete(labels=c("Both","Male\\nonly","Female\\nonly","Neither"))



###Brood size and elevation
ggplot(data=bf.mm,aes(x=BROOD,y=MEANMASS)) + geom_smooth(method="lm",color="black") +
  geom_point(aes(color=ELEVATION)) + theme_cowplot() +
  scale_y_continuous(breaks=c(8,9,10,11,12,13,14,15,16),limits=c(8,16)) + facet_wrap(facets="elevation.title") +
  scale_color_manual(values=c("steelblue3","darkorange1")) + theme(legend.position = "") + xlab("Brood size") + ylab("Mean nestling mass")



#____________________________________________



####Plot CV fancy####

##Get facet titles for plotting
bf.cv = bf.cv %>% mutate(elevation.title=ifelse(ELEVATION=="H","High Elevation","Low Elevation"))


###Year and elevation
##Get groups into list for plotting
bf.cv.ylist = bf.cv %>% group_by(ELEVATION,YEAR) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.cv.ysize = bf.cv %>% group_by(ELEVATION,YEAR) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(3,4,6,5,7,5,2,3,6,6,6,5)) %>% mutate(max=max1/10)

##Plot year
ggplot(data=bf.cv,aes(x=YEAR,y=CV)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.cv.ylist[[1]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[1]) +
  stat_halfeye(data=bf.cv.ylist[[2]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[2]) +
  stat_halfeye(data=bf.cv.ylist[[3]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[3]) +
  stat_halfeye(data=bf.cv.ylist[[4]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[4]) +
  stat_halfeye(data=bf.cv.ylist[[5]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[5]) +
  stat_halfeye(data=bf.cv.ylist[[6]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[6]) +
  stat_halfeye(data=bf.cv.ylist[[7]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[7]) +
  stat_halfeye(data=bf.cv.ylist[[8]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[8]) +
  stat_halfeye(data=bf.cv.ylist[[9]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[9]) +
  stat_halfeye(data=bf.cv.ylist[[10]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[10]) +
  stat_halfeye(data=bf.cv.ylist[[11]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[11]) +
  stat_halfeye(data=bf.cv.ylist[[12]],aes(x=YEAR,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.ysize$max[12]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + xlab("Year") + ylab("CV of nestling mass")



###Age and elevation
##Get groups into list for plotting
bf.cv.alist = bf.cv %>% group_by(ELEVATION,M_F.AGE) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.cv.asize = bf.cv %>% group_by(ELEVATION,M_F.AGE) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(14,4,3,1,13,3,3,1)) %>% mutate(max=max1/18)

##Plot age
ggplot(data=bf.cv,aes(x=M_F.AGE,y=CV)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  geom_boxplot(outlier.alpha = 0,width=0.2,aes(fill=ELEVATION),position = position_nudge(x=-0.15)) + 
  stat_halfeye(data=bf.cv.alist[[1]],aes(x=M_F.AGE,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.asize$max[1]) +
  stat_halfeye(data=bf.cv.alist[[2]],aes(x=M_F.AGE,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.asize$max[2]) +
  stat_halfeye(data=bf.cv.alist[[3]],aes(x=M_F.AGE,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.asize$max[3]) +
  stat_halfeye(data=bf.cv.alist[[4]],aes(x=M_F.AGE,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.asize$max[4]) +
  stat_halfeye(data=bf.cv.alist[[5]],aes(x=M_F.AGE,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.asize$max[5]) +
  stat_halfeye(data=bf.cv.alist[[6]],aes(x=M_F.AGE,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.asize$max[6]) +
  stat_halfeye(data=bf.cv.alist[[7]],aes(x=M_F.AGE,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.asize$max[7]) +
  stat_halfeye(data=bf.cv.alist[[8]],aes(x=M_F.AGE,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.asize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.015,"npc")) + 
  xlab("Male_Female Age") + ylab("CV of nestling mass") + xlab("Parents age") +
  scale_x_discrete(labels=c("Both\\nAdult","Adult M\\nJuv F","Juv M.\\nAdult F","Both\\nJuv"))



###Detected and elevation
##Get groups into list for plotting
bf.cv.dlist = bf.cv %>% group_by(ELEVATION,M_F.detected) %>% group_split()

##Calculate sizes for rainclouds based on largest number of points in each group at each y-value
#Then divide to get a more reasonable number for plotting
#Did this manually based on max number of points per row in plot - not easy to determine ahead of time by rounding. 
#Bin width in plot for points is set by plot size. 
bf.cv.dsize = bf.cv %>% group_by(ELEVATION,M_F.AGE) %>%  summarise(count=n()) %>% 
  add_column(max1 = c(11,5,4,2,5,3,2,11)) %>% mutate(max=max1/10)

##Plot detected 
ggplot(data=bf.cv,aes(x=M_F.detected,y=CV)) + 
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  scale_fill_manual(values=c("steelblue3","darkorange1")) + theme_cowplot() + 
  theme(legend.position = "") +
  stat_halfeye(data=bf.cv.dlist[[1]],aes(x=M_F.detected,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.dsize$max[1]) +
  stat_halfeye(data=bf.cv.dlist[[2]],aes(x=M_F.detected,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.dsize$max[2]) +
  stat_halfeye(data=bf.cv.dlist[[3]],aes(x=M_F.detected,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.dsize$max[3]) +
  stat_halfeye(data=bf.cv.dlist[[4]],aes(x=M_F.detected,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.dsize$max[4]) +
  stat_halfeye(data=bf.cv.dlist[[5]],aes(x=M_F.detected,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.dsize$max[5]) +
  stat_halfeye(data=bf.cv.dlist[[6]],aes(x=M_F.detected,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.dsize$max[6]) +
  stat_halfeye(data=bf.cv.dlist[[7]],aes(x=M_F.detected,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.dsize$max[7]) +
  stat_halfeye(data=bf.cv.dlist[[8]],aes(x=M_F.detected,y=CV,fill=ELEVATION),adjust=2,.width=0,point_colour=NA,alpha=0.5,width=bf.cv.dsize$max[8]) +
  stat_dots(side="right",aes(fill=ELEVATION,color=ELEVATION),dotsize=unit(0.8,"npc"),binwidth=unit(0.02,"npc")) + 
  geom_boxplot(outlier.alpha = 0,width=0.3,aes(fill=ELEVATION),position = position_nudge(x=-0.2)) + 
  xlab("Parents detected") + ylab("CV of nestling mass") +
  scale_x_discrete(labels=c("Both","Male\\nonly","Female\\nonly","Neither"),expand=expansion(mult=c(0.2,0.1))) 



###Brood size and elevation
ggplot(data=bf.cv,aes(x=BROOD,y=CV)) + geom_smooth(method="lm",color="black") +
  geom_point(aes(color=ELEVATION)) + theme_cowplot() +
  facet_wrap(facets="elevation.title") + scale_color_manual(values=c("steelblue3","darkorange1")) +
  theme(legend.position = "") + xlab("Brood size") + ylab("CV of nestling mass") 



#____________________________________________



####Sample sizes graph####

###Sample size
bf.fe.table = bf.fe %>% group_by(elevation.title,YEAR,M_F.detected, .drop=F) %>% summarise(n=n()) %>% mutate(Variable="First egg date")
bf.cs.table = bf.cs %>% group_by(elevation.title,YEAR,M_F.detected, .drop=F) %>% summarise(n=n()) %>% mutate(Variable="Clutch size")
bf.bs.table = bf.bs %>% group_by(elevation.title,YEAR,M_F.detected, .drop=F) %>% summarise(n=n()) %>% mutate(Variable="Brood size")
bf.mm.table = bf.mm %>% group_by(elevation.title,YEAR,M_F.detected, .drop=F) %>% summarise(n=n()) %>% mutate(Variable="Mean nestling mass")
bf.cv.table = bf.cv %>% group_by(elevation.title,YEAR,M_F.detected, .drop=F) %>% summarise(n=n()) %>% mutate(Variable="CV of nestling mass")
bf.all.table = rbind(bf.fe.table,bf.cs.table,bf.bs.table,bf.mm.table,bf.cv.table)

##Reorder factors for plotting
bf.all.table$Variable = factor(bf.all.table$Variable,levels=c("First egg date","Clutch size","Brood size","Mean nestling mass","CV of nestling mass"))

##Plot sample sizes
ggplot(data=bf.all.table,aes(x=M_F.detected,y=n,group=Variable,fill=Variable)) + geom_bar(stat="identity",position=position_dodge()) + 
  theme_cowplot() + facet_wrap(facets=c("elevation.title","YEAR")) + scale_fill_grey(start=0,end=0.8) +
  xlab("Parents fed") + ylab("Sample size") + scale_x_discrete(labels=c("Both","Male\\nonly","Female\\nonly","Neither"))




