## Hedgerows Project, Cleaned up data and analysis June-July 2019

##############################
# Load needed packages
library(tidyverse)
library(lme4)
library(emmeans)
##############################

### Analyses for "Plant species variation in pathogen transmission shapes colony-level bee 
### pathogens and performance" (or something similarly titled)
### Lynn S. Adler, Nicholas A. Barber, Olivia M. Biller and Rebecca E. Irwin



###############################
## Pathogen Infection analysis - only using bees in infected treatment
###############################

# load and wrangle data
cprop<-read.table("crithidia_infection_data.txt",header=T,na.strings=".")
treats<-read.table('flower_strip_treatments.txt',header=T,na.strings=".")

propjoin<-left_join(cprop,treats,by=c("tent"))
propjoin$round<-as.factor(propjoin$round)
propjoin$block<-as.factor(propjoin$block)
propjoin$tent<-as.factor(propjoin$tent)
propjoinI<-subset(propjoin,infect=='I')
nect2<-read.table("nectarresources.txt",header=T,na.strings=".")
nect2$tent<-as.factor(nect2$tent)
cjoinI<-left_join(propjoinI,nect2)
cjoinI$round<-as.factor(cjoinI$round)
cjoinI$block<-as.factor(cjoinI$block)
cjoinI$tent<-as.factor(cjoinI$tent)
cjoinI$nectar1<-scale(cjoinI$allnectar)
cjoinI$nectar2<-scale(cjoinI$noncanola)
cjoinI$nectar3<-scale(cjoinI$allnogold)
cjoinI$nectar4<-scale(cjoinI$noncanolanogold)

############################
# Proportion bees infected #
############################

prop1.1<-glmer(cbind(infectedbees,(totalbees-infectedbees))~treatment+(1|round/block),data=propjoinI,na.action=na.omit,family=binomial)
prop1.2<-glmer(cbind(infectedbees,(totalbees-infectedbees))~(1|round/block),data=propjoinI,na.action=na.omit,family=binomial)
anova(prop1.1,prop1.2) #no significant effect of flower strips treatment

########################
# Mean infection level #
########################

#without nectar
mall<-lmer(meancrithall2~treatment+(1|round/block),data=propjoinI,na.action=na.omit)
mall2<-lmer(meancrithall2~(1|round/block),data=propjoinI,na.action=na.omit)
anova(mall,mall2)
lsmeans(mall,pairwise~treatment)

#with nectar
malln<-lmer(meancrithall2~treatment+nectar1+(1|round/block),data=subset(cjoinI, nectar1!='NA'),na.action=na.omit)
malln2<-lmer(meancrithall2~nectar1+(1|round/block),data=subset(cjoinI, nectar1!='NA'),na.action=na.omit)
anova(malln,malln2)
malln3<-lmer(meancrithall2~treatment+(1|round/block),data=subset(cjoinI, nectar1!='NA'),na.action=na.omit)
anova(malln,malln3)



###############################
## Microcolony performance variables, without and with nectar as a covariate
## variables investigated:  egg number, mean larvae wt, number of adults, larvae number, mean egg 
## wt.
###############################

# Load and wrangle data
perfdata<-read.table("microcolony_performance_data.txt",header=T,na.strings=".")
nect2<-read.table("nectarresources.txt",header=T,na.strings=".")
perfjoin<-left_join(perfdata,nect2)
perfdata<-perfjoin
perfdata$block<-as.factor(perfdata$block)
perfdata$round<-as.factor(perfdata$round)

covariates<-read.table("microcolony_data.txt",header=T,na.strings=".")
covjoin<-left_join(perfdata,covariates)
perfdata<-covjoin

perfdata$nectar1<-scale(perfdata$allnectar)
perfdata$nectar2<-scale(perfdata$noncanola)
perfdata$nectar3<-scale(perfdata$allnogold)
perfdata$nectar4<-scale(perfdata$noncanolanogold)
perfdata$Obs<-c(1:88)

perfdata2<-subset(perfdata,nectar1!='NA' & egg_lay_age_at_deployment!='NA')
perfdata2$eggslaid<-perfdata2$egg_lay_age_at_deployment
perfdata2$numworkers<-perfdata2$num_workers_at_deployment

##################
### Egg number ###
##################

#egg number without nectar:
egg1b<-glmer(eggs~treatment*infect+(1|round/block)+(1|Obs),data=perfdata,family=poisson)
sum(residuals(egg1b,type="pearson")^2)/77 #not overdispersed
#not converging; change optimizer
ss <- getME(egg1b,c("theta","fixef"))
egg2 <- update(egg1b,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
sum(residuals(egg2,type="pearson")^2)/77 #not overdispersed
#converges well
egg2b<-glmer(eggs~treatment+infect+(1|round/block)+(1|Obs),data=perfdata,family=poisson)
ss <- getME(egg2b,c("theta","fixef"))
egg3 <- update(egg2b,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(egg2,egg3) #no interaction effect
egg3b<-glmer(eggs~infect+(1|round/block)+(1|Obs),data=perfdata,family=poisson)
ss <- getME(egg3b,c("theta","fixef"))
egg3b<-update(egg3b,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(egg3,egg3b) #no treatment effect
egg3c<-glmer(eggs~1+(1|round/block)+(1|Obs),data=perfdata,family=poisson)
anova(egg3b,egg3c) #marg sig effect, infection decreases egg production

#egg number, including nectar as covariate:
egg3d<-glmer(eggs~nectar1+treatment*infect+(1|round/block)+(1|Obs),data=subset(perfdata,nectar1!='NA') ,family=poisson)
ss <- getME(egg3d,c("theta","fixef"))
egg3e<-update(egg3d,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
egg3f<-glmer(eggs~nectar1+treatment+infect+(1|round/block)+(1|Obs),data=subset(perfdata,nectar1!='NA') ,family=poisson)
ss <- getME(egg3f,c("theta","fixef"))
egg3g<-update(egg3f,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(egg3e,egg3g) #no interaction
egg3hb<-glmer(eggs~nectar1+infect+(1|round/block)+(1|Obs),data=subset(perfdata,nectar1!='NA') ,family=poisson) #or drop treatment instead, no effect
ss <- getME(egg3hb,c("theta","fixef"))
egg3hb<-update(egg3hb,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(egg3g,egg3hb) #no effect of treatment
egg3i<-glmer(eggs~nectar1+(1|round/block)+(1|Obs),data=subset(perfdata,nectar1!='NA') ,family=poisson) 
anova(egg3hb,egg3i) #marginal effect of infect
egg3j<-glmer(eggs~(1|round/block)+(1|Obs),data=subset(perfdata,nectar1!='NA') ,family=poisson) 
anova(egg3j,egg3hb) #significant positive effect of nectar

##########################
### Mean larval weight ###
##########################

#mean larval weight, without nectar covariate:
larvwt1<-lmer(log(mnlarvwt)~ treatment*infect+(1|round/block), data=subset(perfdata,nectar1!='NA'))
larvwt2<-lmer(log(mnlarvwt)~ treatment+infect+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(larvwt1,larvwt2) #marg sig interaction effect; driven by high average weight of larvae in uninfected high-transmission tents
larvwt3<-lmer(log(mnlarvwt)~ infect+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(larvwt2,larvwt3)
larvwt4<-lmer(log(mnlarvwt)~ (1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(larvwt4,larvwt3)

#mean larval weight, including nectar as a covariate:
larvwt1<-lmer(log(mnlarvwt)~ nectar1+treatment*infect+(1|round/block), data=subset(perfdata,nectar1!='NA'))
larvwt2<-lmer(log(mnlarvwt)~ nectar1+treatment+infect+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(larvwt1,larvwt2) #marg sig interaction effect, again due to high average weight of larvae in uninfected high-transmission tents
larvwt3<-lmer(log(mnlarvwt)~ nectar1+infect+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(larvwt2,larvwt3)
larvwt4<-lmer(log(mnlarvwt)~ nectar1+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(larvwt3,larvwt4)
larvwt5<-lmer(log(mnlarvwt)~ (1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(larvwt4,larvwt5)

##################################
### Number of surviving adults ###
##################################

# includes starting number of adults at beginning of round as a covariate
# without nectar as a covariate:
bees1<-glmer(beescounted~treatment*infect+scale(workers_at_start)+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees1,c("theta","fixef"))
bees1b <- update(bees1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
bees2<-glmer(beescounted~treatment+infect+scale(workers_at_start)+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees2,c("theta","fixef"))
bees2b <- update(bees2,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
bees3<-glmer(beescounted~infect+scale(workers_at_start)+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees3,c("theta","fixef"))
bees3b <- update(bees3,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
bees4<-glmer(beescounted~treatment+scale(workers_at_start)+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees4,c("theta","fixef"))
bees4b <- update(bees4,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
bees5<-glmer(beescounted~scale(workers_at_start)+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees5,c("theta","fixef"))
bees5b <- update(bees5,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(bees1b,bees2b) #no treatment x infection interaction
anova(bees2b,bees3b) #significant treatment effect
anova(bees2b,bees4b) #no infection effect

#with nectar1 as covariate
bees1<-glmer(beescounted~nectar1+treatment*infect+workers_at_start+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees1,c("theta","fixef"))
bees11 <- update(bees1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
bees1b<-glmer(beescounted~nectar1+treatment+infect+workers_at_start+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees1b,c("theta","fixef"))
bees11b <- update(bees1b,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
bees1c<-glmer(beescounted~nectar1+infect+workers_at_start+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
bees1e<-glmer(beescounted~nectar1+workers_at_start+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
bees1f<-glmer(beescounted~workers_at_start+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees1f,c("theta","fixef"))
bees1f <- update(bees1f,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
bees1g<-glmer(beescounted~treatment+workers_at_start+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(bees1g,c("theta","fixef"))
bees1g <- update(bees1g,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(bees11,bees11b) #no treatment x infection interaction
anova(bees11b,bees1c) #no treatment effect
anova(bees1c,bees1e) #no infection effect
anova(bees1e,bees1f) #nectar highly significant, positively correlated with surviving adults
summary(bees1e)

########################
### Number of larvae ###
########################

# without nectar:
larv1<-glmer(larvae~treatment*infect+(1|round/block),data=perfdata,family=poisson)
sum(residuals(larv1,type="pearson")^2)/78 #overdispersed
Obs<-c(1:88)
larv1<-glmer(larvae~treatment*infect+(1|round/block)+(1|Obs),data=perfdata,family=poisson)
larv2<-glmer(larvae~treatment+infect+(1|round/block)+(1|Obs),data=perfdata,family=poisson)
anova(larv1,larv2) #no treatment x infection interaction
larv3<-glmer(larvae~infect+(1|round/block)+(1|Obs),data=perfdata,family=poisson)
anova(larv2,larv3) #treatment significant
larv4<-glmer(larvae~treatment+(1|round/block)+(1|Obs),data=perfdata,family=poisson)
anova(larv2,larv4) #no infection effect
emmeans(larv4, pairwise~treatment, transform="response")

# with nectar
larv1<-glmer(larvae~nectar1+treatment*infect+(1|round/block), data=subset(perfdata,nectar1!='NA'),family=poisson)
sum(residuals(larv1,type="pearson")^2)/76 #overdispersed
perfdata$Obs<-c(1:88)
larv1<-glmer(larvae~nectar1+treatment*infect+(1|round/block)+(1|Obs), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(larv1,c("theta","fixef"))
larv1b <- update(larv1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
larv2<-glmer(larvae~nectar1+treatment+infect+(1|round/block)+(1|Obs), data=subset(perfdata,nectar1!='NA'),family=poisson)
ss <- getME(larv2,c("theta","fixef"))
larv2b <- update(larv2,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(larv1b,larv2b) #no treatment x infection interaction
larv3<-glmer(larvae~nectar1+infect+(1|round/block)+(1|Obs),data=subset(perfdata,nectar1!='NA'),family=poisson)
anova(larv2b,larv3) #marginally significant treatment
larv3b<-glmer(larvae~treatment+nectar1+(1|round/block)+(1|Obs), data=subset(perfdata,nectar1!='NA'),family=poisson)
anova(larv2b,larv3b) #no effect of infection
larv4b<-glmer(larvae~treatment+(1|round/block)+(1|Obs), data=subset(perfdata,nectar1!='NA'),family=poisson)
anova(larv3b,larv4b) #significant positive effect of nectar

#######################
### Mean egg weight ###
#######################

# without nectar
eggwt1<-lmer(log(mneggwt)~ treatment*infect+(1|round/block),data=subset(perfdata,nectar1!='NA'))
eggwt2<-lmer(log(mneggwt)~ treatment+infect+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(eggwt1,eggwt2) # no treatment x infection interaction
eggwt3<-lmer(log(mneggwt)~ infect+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(eggwt2,eggwt3) #treatment significant
eggwt4<-lmer(log(mneggwt)~ treatment+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(eggwt2,eggwt4) #infection not significant

#with nectar
eggwt1<-lmer(log(mneggwt)~ treatment*infect+nectar1+(1|round/block),data=subset(perfdata,nectar1!='NA'))
eggwt2<-lmer(log(mneggwt)~ treatment+infect+nectar1+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(eggwt1,eggwt2) # no treatment x infection interaction
eggwt3<-lmer(log(mneggwt)~ infect+nectar1+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(eggwt2,eggwt3) #treatment not significant
eggwt4<-lmer(log(mneggwt)~ nectar1+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(eggwt3,eggwt4) #infection not significant
eggwt5<-lmer(log(mneggwt)~1+(1|round/block),data=subset(perfdata,nectar1!='NA'))
anova(eggwt4,eggwt5) #significant positive effect of nectar



###############################
## Pollination services and pollinator behavior from quick observations
## Response variables: number of foragers on canola, number of foragers on all plants, pollen
## grains deposited, fruit weight, fruit set, seed weight
###############################

# load and wrangle data
qdata<-read.table("pollinator_observation_data.txt",header=T,na.strings=".")
qdata$tent<-factor(qdata$tent)
qdata$block<-factor(qdata$block)
qdata$round<-factor(qdata$round)

# add in average nectar resources for potential covariate
nect2<-read.table("nectarresources.txt",header=T,na.strings=".")
nect2$tent<-as.factor(nect2$tent)
qjoin<-left_join(qdata,nect2)
qjoin$nectar1<-scale(qjoin$allnectar)
qjoin$nectar2<-scale(qjoin$noncanola)
qjoin$nectar3<-scale(qjoin$allnogold)
qjoin$nectar4<-scale(qjoin$noncanolanogold)
qjoin$tent<-as.factor(qjoin$tent)
#note that flower strip treatment called "hedge" in this code

####################################
### Number of foragers on canola ###
####################################

canmod1<-glmer(can~hedge*infect + (1|round/block/tent), na.action=na.omit, data=qdata, family=poisson)
canmod2<-glmer(can~hedge+infect + (1|round/block/tent), na.action=na.omit, data=qdata, family=poisson)
anova(canmod1,canmod2) #no significant treatment x infection interaction
canmod3<-glmer(can~infect + (1|round/block/tent), na.action=na.omit, data=qdata, family=poisson)
anova(canmod2,canmod3) #significant treatment effect
canmod4<-glmer(can~hedge + (1|round/block/tent), na.action=na.omit, data=qdata, family=poisson)
anova(canmod2,canmod4) #significant infection effect
lsmeans(canmod1, pairwise ~ infect | hedge, type="response")

########################################
### Number of foragers on all plants ###
########################################

summod1<-glmer(sum~hedge*infect + (1|round/block/tent), na.action=na.omit, data=qdata, family=poisson)
sum(residuals(summod1,type="pearson")^2)/455 #not overdispersed
summod2<-glmer(sum~hedge+infect + (1|round/block/tent), na.action=na.omit, data=qdata, family=poisson)
anova(summod1,summod2) #significant treatment x infection interaction
lsmeans(summod1, pairwise ~ infect | hedge, type="response")

###############################
### Pollen grains deposited ###
###############################

#load data and join treatments to emasculated data (with pollen counts)
emasc<-read.table("emasculated_flowers_data.txt",header=T,na.strings=".")
treats<-read.table("flower_strip_treatments.txt",header=T,na.strings=".")
ejoin<-left_join(emasc,treats,by=c("tent"))

depmod1.1<-glmer(canpollen~treatment*infect + (1|round/block/tent/polflower), na.action=na.omit, data=ejoin, family=poisson)
depmod1.2<-glmer(canpollen~treatment+infect + (1|round/block/tent/polflower), na.action=na.omit, data=ejoin, family=poisson)
anova(depmod1.1,depmod1.2) #treatment x infection interaction not significant
depmod1.3<-glmer(canpollen~infect + (1|round/block/tent/polflower), na.action=na.omit, data=ejoin, family=poisson)
anova(depmod1.2,depmod1.3) #flower strips treatment significant
depmod1.4<-glmer(canpollen~treatment + (1|round/block/tent/polflower), na.action=na.omit, data=ejoin, family=poisson)
anova(depmod1.2,depmod1.4) #marginally significant effect of infection
lsmeans(depmod1.4, pairwise ~ treatment)

#################
### Fruit set ###
#################

# load and organize fruit set data
hdata<-read.table('fruitset_data.txt',header=T,na.strings=".")
hdata$tent<-factor(hdata$tent)
hdata$block<-factor(hdata$block)
hdata$round<-factor(hdata$round)
hdata$plantdate<-factor(hdata$plantdate)

fsmod1<-glmer(cbind(frt_y,frt_n)~hedge*infect*handpol + (1|round/block/tent/plantdate), na.action=na.omit, data=hdata, family=binomial)
ss <- getME(fsmod1,c("theta","fixef"))
fsmod1b <- update(fsmod1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

fsmod2 <- glmer(cbind(frt_y,frt_n)~(hedge+infect+handpol)^2 + (1|round/block/tent/plantdate), na.action=na.omit, data=hdata, family=binomial)
ss <- getME(fsmod2,c("theta","fixef"))
fsmod2b <- update(fsmod2,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(fsmod1b,fsmod2b) #marginally significant 3-way treatment x infection x hand pollination interaction, but contrasts show no HP-OP differences for any flower composition/infection treatment combination:
lsmeans(fsmod1b, pairwise ~ handpol | hedge * infect)

fsmod3 <- glmer(cbind(frt_y,frt_n)~ handpol*hedge +handpol*infect + (1|round/block/tent/plantdate), na.action=na.omit, data=hdata, family=binomial)
ss <- getME(fsmod3,c("theta","fixef"))
fsmod3b <- update(fsmod3,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(fsmod2b,fsmod3b) #treatment x infection interaction not significant

fsmod4 <- glmer(cbind(frt_y,frt_n)~ handpol*hedge + infect + (1|round/block/tent/plantdate), na.action=na.omit, data=hdata, family=binomial)
ss <- getME(fsmod4,c("theta","fixef"))
fsmod4b <- update(fsmod4,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(fsmod3b,fsmod4b) #treatment x hand pollination interaction not significant

fsmod5 <- glmer(cbind(frt_y,frt_n)~ handpol+hedge + infect + (1|round/block/tent/plantdate), na.action=na.omit, data=hdata, family=binomial)
ss <- getME(fsmod5,c("theta","fixef"))
fsmod5b <- update(fsmod5,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(fsmod4b,fsmod5b) #treatment x hand pollination interaction not significant

fsmod6 <- glmer(cbind(frt_y,frt_n)~ handpol + infect + (1|round/block/tent/plantdate), na.action=na.omit, data=hdata, family=binomial)
anova(fsmod5b,fsmod6) #treatment not significant

fsmod7 <- glmer(cbind(frt_y,frt_n)~ handpol + (1|round/block/tent/plantdate), na.action=na.omit, data=hdata, family=binomial)
anova(fsmod7,fsmod6) #no effect of infection

fsmod8 <- glmer(cbind(frt_y,frt_n)~ (1|round/block/tent/plantdate), na.action=na.omit, data=hdata, family=binomial)
anova(fsmod7,fsmod8) #no effect of hand pollination

####################
### Fruit weight ###
####################

#load and organize data
sdata<-read.table('seedset_data.txt',header=T,na.strings=".")
sdata$tent<-factor(sdata$tent)
sdata$block<-factor(sdata$block)
sdata$round<-factor(sdata$round)
sdata$plant<-factor(sdata$plant)
sdata$seedwt<-sdata$devsdswt_g/sdata$devsds

#drop seed weights that were misrecorded
sdata$seedwt[25]<-NA
sdata$seedwt[50]<-NA
sdata$seedwt[86]<-NA
sdata$seedwt[299]<-NA
sdata$seedwt[346]<-NA

fwmod<-lmer(fruitwt_g~hedge*infect*hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata) #residuals heteroscedastic; log-transformation reduces this
fwmod1<-lmer(log(fruitwt_g)~hedge*infect*hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
fwmod2<-lmer(log(fruitwt_g)~(hedge+infect+hp_op)^2 +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(fwmod1,fwmod2) #no sigificant 3-way interaction
fwmod3<-lmer(log(fruitwt_g)~hp_op*hedge+infect*hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(fwmod2,fwmod3) #no significant treatment x infection interaction
fwmod4<-lmer(log(fruitwt_g)~hedge+infect*hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(fwmod3,fwmod4) #no significant treatment x hand pollination interaction
fwmod5<-lmer(log(fruitwt_g)~hedge+infect+hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(fwmod4,fwmod5) #no significant infection x hand pollination interaction
fwmod6<-lmer(log(fruitwt_g)~infect+hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(fwmod6,fwmod5) #no significant treatment effect
fwmod7<-lmer(log(fruitwt_g)~hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(fwmod6,fwmod7) #no significant infection effect
fwmod8<-lmer(log(fruitwt_g)~ (1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(fwmod8,fwmod7) #no significant hand pollination effect

###################
### Seed weight ###
###################

swmod1<-lmer(seedwt~hedge*infect*hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
swmod2<-lmer(seedwt~(hedge+infect+hp_op)^2 +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(swmod1,swmod2) #no significan 3-way interaction
swmod3<-lmer(seedwt~hp_op*hedge+infect*hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(swmod2,swmod3) #no significant treatment x infection interaction
swmod4<-lmer(seedwt~hedge+infect*hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(swmod3,swmod4) #no significant treatment x hand pollination interaction
swmod5<-lmer(seedwt~hedge+infect+hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(swmod4,swmod5) #no significant infection x hand pollination interaction
swmod6<-lmer(seedwt~infect+hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(swmod6,swmod5) #no significant treatment effect
swmod7<-lmer(seedwt~hp_op +(1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(swmod6,swmod7) #no significant infection effect
swmod8<-lmer(seedwt~ (1|round/block/tent/plant), na.action=na.omit, data=sdata)
anova(swmod8,swmod7) #no significant hand pollination effect

################
### Seed set ###
################

ssmod1<-glmer(cbind(devsds,undevsds)~hedge*infect*hp_op + (1|round/block/tent/plant), na.action=na.omit, family=binomial, data=sdata)
ss <- getME(ssmod1,c("theta","fixef"))
ssmod1b <- update(ssmod1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
ssmod2 <- glmer(cbind(devsds,undevsds)~(hedge+infect+hp_op)^2 + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=binomial)
ss <- getME(ssmod2,c("theta","fixef"))
ssmod2b <- update(ssmod2,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(ssmod1b,ssmod2b) #no significant 3-way interaction

ssmod3 <- glmer(cbind(devsds,undevsds)~hp_op*hedge +hp_op*infect + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=binomial)
ss <- getME(ssmod3,c("theta","fixef"))
ssmod3b <- update(ssmod3,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(ssmod2b,ssmod3b) #no significant treatment x infection interaction

ssmod4 <- glmer(cbind(devsds,undevsds)~hedge +hp_op*infect + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=binomial)
ss <- getME(ssmod4,c("theta","fixef"))
ssmod4b <- update(ssmod4,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(ssmod4b,ssmod3b) #no significant treatment x hand pollinatinon interaction

ssmod5 <- glmer(cbind(devsds,undevsds)~hedge + hp_op + infect + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=binomial)
ss <- getME(ssmod5,c("theta","fixef"))
ssmod5b <- update(ssmod5,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(ssmod4b,ssmod5b) #no significant infection x hand pollination interaction

ssmod6 <- glmer(cbind(devsds,undevsds)~ hp_op + infect + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=binomial)
anova(ssmod5b, ssmod6) #no significant treatment effect

ssmod7 <- glmer(cbind(devsds,undevsds)~ hp_op + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=binomial)
anova(ssmod7, ssmod6) #no significant infection effect

ssmod8 <- glmer(cbind(devsds,undevsds)~ (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=binomial)
anova(ssmod7, ssmod8) #significant hand pollination effect
lsmeans(ssmod7, pairwise ~ hp_op, type="response")

###################
### Seed number ###
###################

sds1<-glmer(devsds~infect*hedge*hp_op + (1|round/block/tent/plant), na.action=na.omit, family=poisson, data=sdata)
ss <- getME(sds1,c("theta","fixef"))
sds1b <- update(sds1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
sds2 <- glmer(devsds~(hedge+infect+hp_op)^2 + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=poisson)
ss <- getME(sds2,c("theta","fixef"))
sds2b <- update(sds2,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(sds1b,sds2b) #no significant 3-way interaction

sds3 <- glmer(devsds~hp_op*hedge +hp_op*infect + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=poisson)
ss <- getME(sds3,c("theta","fixef"))
sds3b <- update(sds3,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(sds2b,sds3b) #no significant treatment x infection interaction

sds4 <- glmer(devsds~hedge +hp_op*infect + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=poisson)
anova(sds4,sds3b) #treatment x hand pollination interaction significant, retained

sds5 <- glmer(devsds~hp_op*hedge +infect + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=poisson)
ss <- getME(sds5,c("theta","fixef"))
sds5b <- update(sds5,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(sds3b,sds5b) #infection x hand pollination interaction not significant

sds6 <- glmer(devsds~hp_op*hedge + (1|round/block/tent/plant), na.action=na.omit, data=sdata, family=poisson)
ss <- getME(sds6,c("theta","fixef"))
sds6b <- update(sds6,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(sds5b,sds6) #no significant infection effect
lsmeans(sds6, pairwise ~ hp_op | hedge)



###############################
## Pollinator Behavior based on detailed observations
## Response variables: total switches per minute between individuals plants regardless of species,
## proportion of switches that are to a new plant species, proportion of plants and proportion of 
## flowers visited that were low- vs. high-transmission.
###############################

#load and organize data
tentswitches<-read.table("pollinator_switches_data.txt",header=T,na.strings=".")
tenttreatments<-read.table("flower_strip_treatments.txt",header=T,na.strings=".")
tentswitches$tent<-factor(tentswitches$tent)
tenttreatments$tent<-factor(tenttreatments$tent)
tenttreatments$block<-factor(tenttreatments$block)
tenttreatments$round<-factor(tenttreatments$round)
switchdata<-left_join(tenttreatments,tentswitches)
switchdata2<-subset(switchdata,treatment!="Canola")
switchdata2$propdifspp <- switchdata2$Sum_of_plant_switches_observed_dif_sp/switchdata2$Sum_of_plant_switches_observed2
#this last variable is the proportion of plant switches observed that were switches to a different plant species

#add nectar data
nect2<-read.table("nectarresources.txt",header=T,na.strings=".")
nect2$tent<-as.factor(nect2$tent)
switchdata2<-left_join(switchdata2,nect2)
switchdata2$nectar1<-scale(switchdata2$allnectar)
switchdata2$obs <- c(1:72)

###########################
### Switches per minute ###
###########################

spmo1 <- lmer(switches_per_min_obs ~ treatment * infect + nectar1 + (1|round/block), data=switchdata2, na.action=na.omit)
spmo2 <- lmer(switches_per_min_obs ~ treatment + infect + nectar1 + (1|round/block), data=switchdata2, na.action=na.omit)
spmo3 <- lmer(switches_per_min_obs ~ treatment + nectar1 + (1|round/block), data=switchdata2, na.action=na.omit)
spmo4 <- lmer(switches_per_min_obs ~ nectar1 + (1|round/block), data=switchdata2, na.action=na.omit)
spmo5 <- lmer(switches_per_min_obs ~ (1|round/block), data=switchdata2, na.action=na.omit)
anova(spmo1,spmo2) #treatment x infection interaction not significant
anova(spmo2,spmo3) #infection not significant
anova(spmo3,spmo4) #treatment not significant
anova(spmo4,spmo5) #significant nectar effect; fewer switches with more nectar

##############################################################
### Proportion of switches that are to a new plant species ###
##############################################################

switchvar <- cbind(switchdata2$Sum_of_plant_switches_observed_dif_sp,(switchdata2$Sum_of_plant_switches_observed2-switchdata2$Sum_of_plant_switches_observed_dif_sp))
switchmod1 <- glmer(switchvar ~ treatment * infect + nectar1 + (1|round/block) + (1|obs), data=switchdata2, family="binomial", na.action=na.omit)
switchmod2 <- glmer(switchvar ~ treatment + infect + nectar1 + (1|round/block) + (1|obs), data=switchdata2, family="binomial", na.action=na.omit)
switchmod3 <- glmer(switchvar ~ treatment + nectar1 + (1|round/block) + (1|obs), data=switchdata2, family="binomial", na.action=na.omit)
switchmod4 <- glmer(switchvar ~ nectar1 + (1|round/block) + (1|obs), data=switchdata2, family="binomial", na.action=na.omit)
switchmod5 <- glmer(switchvar ~ (1|round/block) + (1|obs), data=switchdata2, family="binomial", na.action=na.omit)

anova(switchmod1,switchmod2) #treatment x infection interaction not signifixant
anova(switchmod2,switchmod3) #infection not significant
anova(switchmod3,switchmod4) #treatment not significant
anova(switchmod4,switchmod5) #nectar not significant

##############################################
### Proportion of plant visits low vs high ###
##############################################

propvis1 <- glmer(cbind(low_plants,high_plants) ~ treatment * infect + nectar1 + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
propvis2 <- glmer(cbind(low_plants,high_plants) ~ treatment + infect + nectar1 + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
propvis3 <- glmer(cbind(low_plants,high_plants) ~ treatment + nectar1 + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
propvis4 <- glmer(cbind(low_plants,high_plants) ~ nectar1 + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
propvis5 <- glmer(cbind(low_plants,high_plants) ~ treatment + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
ss <- getME(propvis5,c("theta","fixef"))
propvis5b <- update(propvis5,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(propvis1,propvis2) #no treatment x infection interaction
anova(propvis2,propvis3) #no infection effect
anova(propvis3,propvis4) #treatment highly significant, driven by more low visits in low treatment
anova(propvis3,propvis5b) #significant nectar effect; fewer low visits as nectar availability increases

################################################
### Proportion of flowers visits low vs high ###
################################################

#adding total nectar as covariate does not change results qualitatively
propfvis1 <- glmer(cbind(low_flowers,high_flowers) ~ treatment * infect + nectar1 + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
propfvis2 <- glmer(cbind(low_flowers,high_flowers) ~ treatment + infect + nectar1 + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
propfvis3 <- glmer(cbind(low_flowers,high_flowers) ~ treatment + nectar1 + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
propfvis4 <- glmer(cbind(low_flowers,high_flowers) ~ nectar1 + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
propfvis5 <- glmer(cbind(low_flowers,high_flowers) ~ treatment + (1|round/block) + (1|obs), data = switchdata2, family=binomial, na.action=na.omit)
anova(propfvis1,propfvis2) #no treatment x infection interaction
anova(propfvis2,propfvis3) #infection does not affect proportion visited
anova(propfvis3,propfvis4) #treatment highly significant; more low visits in low treatment
anova(propfvis3,propfvis5) #nectar significant; fewer low visits as nectar availability increases

##################################
### Thyme and Sunflower visits ###
##################################

#load data
sunthydata <- read.table("sun_thy_visit_data.txt",header=T,na.strings=".")
tenttreatments<-read.table("flower_strip_treatments.txt",header=T,na.strings=".")
sunthydata$tent<-factor(sunthydata$tent)
tenttreatments$tent<-factor(tenttreatments$tent)
tenttreatments$block<-factor(tenttreatments$block)
tenttreatments$round<-factor(tenttreatments$round)
sunthydata2<-left_join(sunthydata,tenttreatments)
sunthydata2$num_sun_flowers<-scale(sunthydata2$num_sun_flowers)
sunthydata2$num_thy_flowers<-scale(sunthydata2$num_thy_flowers)
sunthydata2$obs <- c(1:212)
#includes "num_sun_flowers" and "num_thy_flowers", which are daily estimates of flower number for those species, as covariates

#sunflower flowers - total visits, as Poisson models
sun1<-glmer(Sun_flowers_visited ~ treatment * infect + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=poisson)
ss <- getME(sun1,c("theta","fixef"))
sun1b <- update(sun1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
sun2<-glmer(Sun_flowers_visited ~ treatment + infect + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=poisson)
sun3<-glmer(Sun_flowers_visited ~ treatment + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=poisson)
sun4<-glmer(Sun_flowers_visited ~  num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=poisson)
anova(sun1,sun2) #no treatment x infection interaction
anova(sun2,sun3) #no infection effect
anova(sun2,sun3) #no treatment effect

#sunflower flowers - proportion of all visits, as binomial models
sunthydata2$sunprop <- cbind(sunthydata2$Sun_flowers_visited, (sunthydata2$Total_flowers_visited-sunthydata2$Sun_flowers_visited))
sun1<-glmer(sunprop ~ treatment * infect + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
ss <- getME(sun1,c("theta","fixef"))
sun1b <- update(sun1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
sun2<-glmer(sunprop ~ treatment + infect + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
sun3<-glmer(sunprop ~ treatment + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
sun4<-glmer(sunprop ~  num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
anova(sun1b,sun2) #no treatment x infection interaction
anova(sun2,sun3) #no infection effect
anova(sun3,sun4) #no treatment effect

#sunflower plants - total visits, as Poisson models
sun1<-glmer(Sun_plants_visited ~ treatment * infect + num_sun_flowers + (1|round/block/tent) , data=sunthydata2, na.action=na.omit, family=poisson)
sun2<-glmer(Sun_plants_visited ~ treatment + infect + num_sun_flowers + (1|round/block/tent) , data=sunthydata2, na.action=na.omit, family=poisson)
ss <- getME(sun2,c("theta","fixef"))
sun2b <- update(sun2,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
sun3<-glmer(Sun_plants_visited ~ treatment + num_sun_flowers + (1|round/block/tent) , data=sunthydata2, na.action=na.omit, family=poisson)
sun4<-glmer(Sun_plants_visited ~  num_sun_flowers + (1|round/block/tent) , data=sunthydata2, na.action=na.omit, family=poisson)
anova(sun1,sun2b) #no treatment x infection interaction
anova(sun2b,sun3) #no infection effect
anova(sun3,sun4) #treatment significant - more sunflower visits in low treatments (where there are more sunflowers)

#sunflower plants - proportion of all visits, as binomial models
sunthydata2$sunprop2 <- cbind(sunthydata2$Sun_plants_visited, (sunthydata2$Total_plants_visited-sunthydata2$Sun_plants_visited))
sun1<-glmer(sunprop2 ~ treatment * infect + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
sun2<-glmer(sunprop2 ~ treatment + infect + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
sun3<-glmer(sunprop2 ~ treatment + num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
sun4<-glmer(sunprop2 ~ num_sun_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
anova(sun1,sun2) #no treatment x infection interaction
anova(sun2,sun3) #no infection effect
anova(sun3,sun4) #treatment significant - higher proportion of sunflower visits in low treatments (where there are more sunflowers)

#thy flowers - total visits, as Poisson models
thy1<-glmer(Thy_flowers_visited ~ treatment * infect + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=poisson)
ss <- getME(thy1,c("theta","fixef"))
thy1b <- update(thy1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
thy2<-glmer(Thy_flowers_visited ~ treatment + infect + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=poisson)
thy3<-glmer(Thy_flowers_visited ~ treatment + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=poisson)
thy4<-glmer(Thy_flowers_visited ~  num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=poisson)
anova(thy1,thy2) #no treatment x infection interaction
anova(thy2,thy3) #no infection effect
anova(thy3,thy4) #no treatment effect

#thy flowers - proportion of all visits, as binomial models
thyprop <- cbind(sunthydata2$Thy_flowers_visited, (sunthydata2$Total_flowers_visited-sunthydata2$Thy_flowers_visited))
thy1<-glmer(thyprop ~ treatment * infect + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
ss <- getME(thy1,c("theta","fixef"))
thy1b <- update(thy1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
thy2<-glmer(thyprop ~ treatment + infect + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
thy3<-glmer(thyprop ~ treatment + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
thy4<-glmer(thyprop ~  num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
anova(thy1,thy2) #no treatment x infection interaction
anova(thy2,thy3) #no infection effect
anova(thy3,thy4) #no treatment effect

#thy plants - total visits, as Poisson models
thy1<-glmer(Thy_plants_visited ~ treatment * infect + num_thy_flowers + (1|round/block/tent), data=sunthydata2, na.action=na.omit, family=poisson)
ss <- getME(thy1,c("theta","fixef"))
thy1b <- update(thy1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
thy2<-glmer(Thy_plants_visited ~ treatment + infect + num_thy_flowers + (1|round/block/tent), data=sunthydata2, na.action=na.omit, family=poisson)
thy3<-glmer(Thy_plants_visited ~ treatment + num_thy_flowers + (1|round/block/tent), data=sunthydata2, na.action=na.omit, family=poisson)
thy4<-glmer(Thy_plants_visited ~  num_thy_flowers + (1|round/block/tent), data=sunthydata2, na.action=na.omit, family=poisson)
anova(thy1b,thy2) #no treatment x infection interaction
anova(thy2,thy3) #no infection effect
anova(thy3,thy4) #no treatment effect

#thy flowers - proportion of all visits, as binomial models
thyprop2 <- cbind(sunthydata2$Thy_plants_visited, (sunthydata2$Total_plants_visited-sunthydata2$Thy_plants_visited))
thy1<-glmer(thyprop2 ~ treatment * infect + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
thy2<-glmer(thyprop2 ~ treatment + infect + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
thy3<-glmer(thyprop2 ~ treatment + num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
thy4<-glmer(thyprop2 ~  num_thy_flowers + (1|round/block/tent) + (1|obs), data=sunthydata2, na.action=na.omit, family=binomial)
anova(thy1,thy2) #no treatment x infection interaction
anova(thy2,thy3) #no infection effect
anova(thy3,thy4) #no treatment effect
