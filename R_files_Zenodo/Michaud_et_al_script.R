# Analyses for Michaud et al., "Pre-infection effects of nectar secondary compounds on a bumble bee gut pathogen"



library(tidyverse)

library(lme4)

library(survival)

library(effects)



####################################
### Thymine-Anabasine experiment ###
####################################

### 
Get data loaded and cleaned up - make numbers be numbers ###

thy.ana.inoc <- read.csv("thymol anabasine inoc counts.csv", na.strings=".", header = TRUE, dec = ",")

thy.ana.inoc$BeeID <- as.factor(thy.ana.inoc$BeeID) #make BeeID a factor

thy.ana.inoc$CallowMass <- as.numeric(as.character(thy.ana.inoc$CallowMass)) #make mass a number 


thy.ana.wings <- read.csv("thymol anabasine wings.csv", na.strings=".", header=TRUE, dec = ",")

thy.ana.wings$BeeID <- as.factor(thy.ana.wings$BeeID) #make BeeID a factor

thy.ana.wings$WingLength <- as.numeric(as.character(thy.ana.wings$ocular_units)) #make wing length a number 



#merge inoc/count and wing dataframes

thy.ana <- merge(thy.ana.inoc, thy.ana.wings, by=c("BeeID"))



#center and scale wing cell lengths

thy.ana$WingLC<-scale(thy.ana$WingLength, center = TRUE, scale = TRUE)



#relevel to have Control as reference level

thy.ana$Treatment <- relevel(thy.ana$Treatment, ref="Control")



### Fit Crithidia count model ###


#create observation-level random factor

thy.ana$obs<-c(1:length(thy.ana$Count))



#initially include and test inoculation date x colony random factor

thy.ana.mod1b <- glmer(Count ~ Treatment + WingLC + Colony + (1|InocDate) + (1|InocDate:Colony) + (1|obs), data=thy.ana, family=poisson)

thy.ana.mod2b <- glmer(Count ~ Treatment + WingLC + Colony + (1|InocDate) + (1|obs), data=thy.ana, family=poisson)

anova(thy.ana.mod1b,thy.ana.mod2b) 

#drop date x colony random factor


#test Treatment

thy.ana.mod3b <- glmer(Count ~ WingLC + Colony + (1|InocDate) + (1|obs), data=thy.ana, family=poisson)

anova(thy.ana.mod2b,thy.ana.mod3b) #treatment not significant



### probability of death ###

#load data

thy.ana.death <- read.csv("thymol anabasine counts death.csv", na.strings=".")

thy.ana.death$Colony <- as.factor(thy.ana.death$Colony)



#binomial model comparing death rates for each treatment

thy.ana.deathmod1 <- glmer(Dead ~ Treatment + (1|Colony), data=thy.ana.death, family=binomial)

thy.ana.deathmod2 <- glmer(Dead ~  (1|Colony), data=thy.ana.death, family=binomial)

anova(thy.ana.deathmod1,thy.ana.deathmod2)



#time to death using cox proportional hazards model:
thyanacox <- coxph(Surv(DaystoDeath, Dead) ~ Treatment, data = thy.ana.death)

summary(thyanacox)

thyanacox2 <- coxph(Surv(DaystoDeath, Dead) ~ Treatment + Colony, data = thy.ana.death)

summary(thyanacox2)



######################################
### Aucubin-Citric acid experiment ###
######################################


#get data loaded and cleaned up - make numbers be numbers

auc.ca.counts <- read.csv("aucubin citric acid counts.csv", na.strings=".", header = TRUE, dec = ",")

auc.ca.counts$BeeID <- as.factor(auc.ca.counts$BeeID) 

#make BeeID a factor

auc.ca.counts$Mass <- as.numeric(as.character(auc.ca.counts$Mass)) #weird that it's reading this as a factor, but I don't see any typos or characters inserted.  



#load wing data as size covariate

auc.ca.wings <- read.csv("aucubin citric acid wings.csv", na.strings=".", header = TRUE, dec = ",")

auc.ca.wings$BeeID <-as.factor(auc.ca.wings$BeeID)

auc.ca.wings$WingLength<-as.numeric(as.character(auc.ca.wings$WingLength))



#merge wing data and main dataframe

auc.ca <- merge(auc.ca.counts,auc.ca.wings,by=c("BeeID"))



#center and scale wing cell length

auc.ca$WingLC<-scale(auc.ca$WingLength, center = TRUE, scale = TRUE)



#relevel to have Control as reference 
level
auc.ca$Treatment <- relevel(auc.ca$Treatment, ref="Control")



### Fit Crithidia count model ###


#create observation-level random factor

auc.ca$obs<-c(1:112)



#initially include and test inoculation date x colony random factor

auc.ca.mod1b <- glmer(Count ~ Treatment + WingLC + Colony + (1|InocDate) + (1|InocDate:Colony) + (1|obs), data=auc.ca, family=poisson)

auc.ca.mod2b <- glmer(Count ~ Treatment + WingLC + Colony + (1|InocDate) + (1|obs), data=auc.ca, family=poisson)

anova(auc.ca.mod1b,auc.ca.mod2b) #can remove date x colony random factor



#test Treatment

auc.ca.mod3b <- glmer(Count ~ WingLC + Colony + (1|InocDate) + (1|obs), data=auc.ca, family=poisson)

anova(auc.ca.mod2b,auc.ca.mod3b)



### probability of death ###

#load data

auc.cit.death <- read.csv("aucubin citric acid counts death.csv", na.strings=".")

auc.cit.death$Colony <- as.factor(auc.cit.death$Colony)



#binomial model comparing death rates for each treatment

auc.cit.deathmod1 <- glmer(Dead ~ Treatment + (1|Colony), data=auc.cit.death, family=binomial)

auc.cit.deathmod2 <- glmer(Dead ~  (1|Colony), data=auc.cit.death, family=binomial)

anova(auc.cit.deathmod1,auc.cit.deathmod2)



#time to death using cox proportional hazards model:

auccitcox <- coxph(Surv(Days_til_death, Dead) ~ Treatment, data = auc.cit.death)

summary(auccitcox)


auccitcox2 <- coxph(Surv(Days_til_death, Dead) ~ Treatment + Colony, data = auc.cit.death)

summary(auccitcox2)




####################################
### Nicotine-Catalpol experiment ###
####################################


#get data loaded and cleaned up - make numbers be numbers

nic.cat.inoc <- read.csv("nicotine catalpol inoculations.csv", na.strings=".", header = TRUE, dec = ",")


nic.cat.inoc$BeeID <- as.factor(nic.cat.inoc$BeeID) #make BeeID a factor

nic.cat.inoc$CallowMass_g <- as.numeric(as.character(nic.cat.inoc$CallowMass_g)) #make mass a number 

nic.cat.counts <- read.csv("nicotine catalpol counts.csv", na.strings=".", header = TRUE, dec = ",")

nic.cat.counts$BeeID <- as.factor(nic.cat.counts$BeeID) #make BeeID a factor



#load wing data as size covariate

nic.cat.wings <- read.csv("nicotine catalpol wing lengths.csv", na.strings=".", header=TRUE, dec = ",")

nic.cat.wings$BeeID <- as.factor(nic.cat.wings$BeeID) #make BeeID a factor

nic.cat.wings$WingLength <- as.numeric(as.character(nic.cat.wings$Ocularunits)) #make wing length a number 



#merge inoc, count, and wing dataframes

nic.cat1 <- merge(nic.cat.inoc, nic.cat.counts, by=c("BeeID"))

nic.cat <- merge(nic.cat1,nic.cat.wings, by=c("BeeID"))



#center and scale wing cell lengths

nic.cat$WingLC<-scale(nic.cat$WingLength, center = TRUE, scale = TRUE)

 

#relevel to have Control as reference level

nic.cat$Treatment.x <- relevel(nic.cat$Treatment.x, ref="Control")



### Fit Crithidia count model ###


#create observation-level random factor

nic.cat$obs<-c(1:length(nic.cat$Count))



#initially include and test inoculation date x colony random factor

nic.cat.mod1b <- glmer(Count ~ Treatment.x + WingLC + ColonyID + (1|InocDate) + (1|InocDate:ColonyID) + (1|obs), data=subset(nic.cat, WingLC!='NA'), family=poisson)
ss <- getME(nic.cat.mod1b,c("theta","fixef"))

nic.cat.mod1b2 <- update(nic.cat.mod1b,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

nic.cat.mod2b <- glmer(Count ~ Treatment.x + WingLC + ColonyID + (1|InocDate) + (1|obs), data=subset(nic.cat, WingLC!='NA'), family=poisson)

ss <- getME(nic.cat.mod2b,c("theta","fixef"))

nic.cat.mod2b2 <- update(nic.cat.mod2b,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))


anova(nic.cat.mod1b2,nic.cat.mod2b2) #can remove date x colony random factor



#test Treatment

nic.cat.mod3b <- glmer(Count ~ WingLC + ColonyID + (1|InocDate) + (1|obs), data=subset(nic.cat, WingLC!='NA'), family=poisson)

ss <- getME(nic.cat.mod3b,c("theta","fixef"))

nic.cat.mod3b2 <- update(nic.cat.mod3b,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))

anova(nic.cat.mod2b2,nic.cat.mod3b2) #Treatment not significant



### probability of death ###


#load data
nicotine.catalpol.counts.death <- read.csv("nicotine catalpol counts death.csv", na.strings=".")

nic.cat.death <- nicotine.catalpol.counts.death



#binomial model comparing death rates for each treatment

nic.cat.deathmod1 <- glmer(Dead ~ Treatment + (1|ColonyID), data=nic.cat.death, family=binomial)

nic.cat.deathmod2 <- glmer(Dead ~  (1|ColonyID), data=nic.cat.death, family=binomial)

anova(nic.cat.deathmod1,nic.cat.deathmod2)



#time to death using cox proportional hazards model

niccatcox <- coxph(Surv(Days_til_death, Dead) ~ Treatment, data = nic.cat.death)

summary(niccatcox)

niccatcox2 <- coxph(Surv(Days_til_death, Dead) ~ Treatment + ColonyID, data = nic.cat.death)

summary(niccatcox2)
