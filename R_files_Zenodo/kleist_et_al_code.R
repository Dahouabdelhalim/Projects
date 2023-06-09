library(psych)
library(effects)
library(lme4)
library(MuMIn)

#################################################
#########ALL species nestling baseline cort#################
#################################################
#datasheet 1 kleist et al


all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)
all_nestling_cort_model <- lmer(dl_base_cort ~ scale(noise)+species+ scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+scale(time)+
                                  (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model)
#remove all random effects but year
all_nestling_cort_model_1<-lmer(dl_base_cort ~ scale(noise)+species+ scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+scale(time)+
                                  (1|year), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)

#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)

###best model
best<-lmer(dl_base_cort~noise +species+(1|year), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best)

null<-lmer(dl_base_cort~1+(1|year), data=all_dat)
anova(best,null, test="F")

############## Confidence intervals for models < 2 AIC #################
confint(best, level=0.85)
r.squaredGLMM(best)
###plot individual plots######
plot(allEffects(best, confidence.level = 0.85),ylab="Baseline cort (ng/ml)", main="Nestling, all species")
plot(effect('noise', best,confidence.level=0.85),ylim=(0:7),ylab=list(label="Baseline cort (ng/ml)",cex=0.8), xlab=list(label="", cex=0.8), main=list(label="Effect of noise nestling baseline stress hormones", cex=0.8))


#################################################
#########ALL species Mom baseline cort#################
#################################################
#datasheet 4 kleist et al
all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)

all_nestling_cort_model <- lmer(dl_base_cort ~ scale(noise)+species+ life_stage+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+scale(time)+
                                  (1|pair)+(1|site)+(1|nest_box)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model)

#keep nest ID and site

all_nestling_cort_model_1<-lmer(dl_base_cort ~ scale(noise)+species+ life_stage+ scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+scale(time)+
                                  (1|site)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)


summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)

#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)

###best model
best<-lmer(dl_base_cort~noise +life_stage+species+(1|site)+(1|nest_ID),
           REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(best)

null<-lmer(dl_base_cort~1+(1|site)+(1|nest_ID),data=all_dat)
anova(best,null, test="F")

r.squaredGLMM(best)
############## Confidence intervals for models < 2 AIC #################
confint(best, level=0.85)
###plot individual plots######
plot(allEffects(best, confidence.level = 0.85),ylab="Baseline cort (ng/ml)", main="Adult female, all species")

#################################################
#########Principal Components Analysis#################
#################################################

##this code was used for PCA on nestling wing, tail, tarsus and mass. However these are already included in the datasheetS

#datasheet  5 kleist et al

### intra-individual variation-create principal component rotated
attach(all_dat)
intra <- matrix(c(weight, tarsus, wing, tail),ncol=4)
fit <- principal(intra,nfactors=2,rotate="varimax",scores=TRUE)
fit
all_pc1 <- fit$scores[,1]
all_pc2 <- fit$scores[,2]

#based on loadings it looks like there is a wing/tail axis (pc1) and a weight/tarsus axis (pc2)


#################################################
#########ALL SPECIES PC1 #################
#################################################
#datasheet  5 kleist et al

all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)

all_nestling_cort_model_quad <- lmer(pc1~ species+poly(noise,2)+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+
                                       (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model_quad)

#remove site, pair
all_nestling_cort_model_1<-lmer(pc1~ species+poly(noise,2)+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+
                                  (1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)

#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)


best<-lmer(pc1~ species+poly(noise,2)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best)

confint(best, level=0.85)
r.squaredGLMM(best)

null<-lmer(pc1~1+(1|nestbox)+(1|year)+(1|nest_ID), data=all_dat)
anova(best,null, test="Chisq")

plot(allEffects(best, confidence.level = 0.95),ylab="Feather growth (PC1)", main="All species")

#################################################
#########All SPECIES PC2 #################
#################################################
#datasheet  5 kleist et al

all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)
#str(all_dat)

all_nestling_cort_model <- lmer(pc2~ species+poly(noise,2)+scale(laydate)+scale(distance)+scale(time)+scale(tree_cover_site)+ scale(chicks)+
                                       (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model)

#keep nestID

all_nestling_cort_model_1 <- lmer(pc2~ species+poly(noise,2)+scale(laydate)+scale(distance)+scale(time)+scale(tree_cover_site)+ scale(chicks)+
                                       (1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)

#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)

#best model
best<-lmer(pc2~poly(noise,2)+scale(chicks)+scale(distance)+scale(laydate)+species+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best)

confint(best, level=0.85)
r.squaredGLMM(best)

null<-lmer(pc2~1+(1|nest_ID), data=all_dat)
anova(best,null, test="F")

plot(allEffects(best, confidence.level = 0.85),ylab="Body size (PC2)", main="All species")

#################################################
#########ALL species nestling acute cort#################
#################################################
#datasheet 1 kleist et al, row 161 has no data for acute cort, removed from analysis
#all_dat
all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)
all_nestling_cort_model <- lmer(acute_cort ~ scale(noise)+species+ scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+scale(time)+
                                  (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model)

#remove everything but nest id

all_nestling_cort_model_1<-lmer(acute_cort ~ scale(noise)+species+ scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+scale(time)+
                                 (1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)


summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)

#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)


###best model
best<-lmer(acute_cort ~ scale(noise)+species+scale(tree_cover_site)+ scale(chicks)+scale(time)+
             (1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(best)
r.squaredGLMM(best)

null<-lmer(acute_cort~1+(1|nest_ID), data=all_dat)
anova(best,null, test="F")

############## Confidence intervals for models < 2 AIC #################
confint(best, level=0.85)
###plot individual plots######
plot(allEffects(best, confidence.level = 0.85),ylab="Acute cort (ng/ml)", main="Nestling, all species")

#################################################
#######################Hatching success###################
#################################################
#datasheet 2 kleist et al for this and the following all species models

rm(list = ls())
all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
##hatching success model
hatch_model <- glmer(hatched~ species+scale(noise)+ scale(laydate)+ scale(distance)+scale(tree_cover_site)+scale(eggs)+
                       (1|site) +(1|pair)+(1|year)+(1|nest_box),data=all_dat, family="binomial")

summary(hatch_model)

#remove all run as glm

hatch_model_1 <- glm(hatched~ species+ scale(noise)+ scale(laydate)+ scale(distance)+scale(tree_cover_site)+scale(eggs)
                     ,data=all_dat, family="binomial")


summary(hatch_model_1)
####dredge on model with random effects removed except species

options(na.action = "na.fail")
all_dredge<- dredge(hatch_model_1)
subset(all_dredge, delta < 2)
#subset of models within 2 AIC of best model

#Intercept is top model
#Run species models

##################################################
##########HATCHED_ATFL############################
##################################################

atfl_hatch <- subset(all_dat, species=="ATFL")
##hatching success model
hatch_model <-  glmer(hatched~ scale(noise)+ scale(laydate)+ scale(distance)+scale(tree_cover_site)+scale(eggs)+
                        (1|site) +(1|pair)+(1|year)+(1|nest_box),data=atfl_hatch, family="binomial")
summary(hatch_model)
#keep year
hatch_model_1 <- glmer(hatched~ scale(noise)+ scale(laydate)+ scale(distance)+scale(tree_cover_site)+scale(eggs)+
                         (1|year),data=atfl_hatch, family="binomial")

summary(hatch_model_1)
####dredge on model with random effects removed except species
options(na.action = "na.fail")
all_dredge<- dredge(hatch_model_1)
#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)

###Top model indistinguishable from intercept

###best model
best <- glmer(hatched~ noise+ (1|year),data=atfl_hatch, family="binomial")

summary(best)
r.squaredGLMM(best)

null<-glmer(hatched~1+(1|year), family="binomial", data=atfl_hatch)
anova(best,null, test="Chisq")

############## Confidence intervals for models < 2 AIC #################
confint(best, level=0.85)

##################################################
##########HATCHED_WEBL############################
##################################################


webl_hatch <- subset(all_dat, species=="WEBL")
##hatching success model
hatch_model <- glmer(hatched~ scale(noise)+ scale(laydate)+ scale(distance)+scale(tree_cover_site)+scale(eggs)+
                       (1|site) +(1|pair)+(1|year)+(1|nest_box),data=webl_hatch, family="binomial")
summary(hatch_model)
#remove all randoms, run as glm
hatch_model_1 <- glm(hatched~ scale(noise)+ scale(laydate)+ scale(distance)+scale(tree_cover_site)+scale(eggs),
                     data=webl_hatch, family="binomial")

summary(hatch_model_1)
####dredge on model with random effects removed except species
options(na.action = "na.fail")
hatch_dredge<- dredge(hatch_model_1)
subset(hatch_dredge, delta < 2)

###best model
best<-glm(hatched~noise+ tree_cover_site,data=webl_hatch, family="binomial")

summary(best)
r.squaredGLMM(best)

null<-glm(hatched~1, family="binomial",data=webl_hatch)

anova(best,null, test="Chisq")

############## Confidence intervals for models < 2 AIC #################
confint(best, level=0.85)
###plot individual plots######
plot(allEffects(best, confidence.level = 0.85),ylab="Body size (PC2)", main="All species")

eff.trans <- effect("noise", best, confidence.level=.85)
plot(eff.trans, type="response", ylab=list(label="Hatch success",cex=0.8), ylim=c(0:1),main=list(label="Effect of noise on WEBL hatching", cex=0.8))

##################################################
##########HATCHED_MOBL############################
##################################################

mobl_hatch <- subset(all_dat, species=="MOBL")
##hatching success model
hatch_model <-  glmer(hatched~ scale(noise)+ scale(laydate)+ scale(distance)+scale(tree_cover_site)+scale(eggs)+
                        (1|site) +(1|pair)+(1|year)+(1|nest_box),data=mobl_hatch, family="binomial")
summary(hatch_model)
#remove all randoms run as glm
hatch_model_1 <- glm(hatched~ scale(noise)+ scale(laydate)+ scale(distance)+scale(tree_cover_site)+scale(eggs)
                     , data=mobl_hatch, family="binomial")

summary(hatch_model_1)
####dredge on model with random effects removed except species
options(na.action = "na.fail")
hatch_dredge<- dredge(hatch_model_1)
subset(hatch_dredge, delta < 2)


#intercept only

#######################Maternal provisioning cort development###################
#######################Maternal provisioning cort development###################
#######################Maternal provisioning cort development###################

#################################################
#########All SPECIES PC1 #################
#################################################
#datasheet 3 kleist et al

##this analysis explores the effects of maternal baseline cort measured during
#provisioning on nestling feather length at day 12

all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)
##subset

all_nestling_cort_model <- 
  lmer(pc1~ scale(mom_base_cort)*species+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+
         (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model)

#keep pair and year

all_nestling_cort_model_1<-
  lmer(pc1~ scale(mom_base_cort)*species+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+
         (1|year)+(1|pair), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)

#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)


##best cort 1
best<-lmer(pc1~ scale(chicks)+scale(distance)+scale(mom_base_cort)*species+scale(tree_cover_site)
           +(1|pair)+(1|year), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best)

confint(best, level=0.85)
r.squaredGLMM(best)

null<-lmer(pc1~1+(1|pair)+(1|year), data=all_dat)
anova(best,null, test="F")

plot(allEffects(best, confidence.level = 0.85),ylab="PCfeather", main="Nestling, all species")

#################################################
#########All SPECIES PC2 #################
#################################################
#datasheet 3 kleist et al


##this analysis explores the effects of maternal baseline cort measured during
#provisioning on nestling body size (tarsus and mass PC2) at day 12
all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)
#str(all_dat)

all_nestling_cort_model <- lmer(pc2~ species*scale(mom_base_cort)+scale(laydate)+scale(distance)+scale(time)+scale(tree_cover_site)+ scale(chicks)+
                                  (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model)

#keep site

all_nestling_cort_model_1 <- lmer(pc2~ species*scale(mom_base_cort)+scale(laydate)+scale(distance)+scale(time)+scale(tree_cover_site)+ scale(chicks)+
                                    (1|site), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)

#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)


#best model
best<-lmer(pc2~scale(distance)+scale(laydate)+scale(mom_base_cort)*species+(1|site), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best)

confint(best, level=0.85)
r.squaredGLMM(best)

null<-lmer(pc2~1+(1|site),data=all_dat)
anova(best,null, test="F")

plot(allEffects(best, confidence.level = 0.85),ylab="Body size (PC2)", main="All species")


#######################ALL PC CORT NESTLING CORT###################
#######################ALL PC CORT NESTLING CORT###################
#######################ALL PC CORT NESTLING CORT###################

#################################################
#########All SPECIES PC1 NESTLING CORT #################
#################################################

#datasheet 1_kleist_et_al

##this analysis explores the effects of nestling baseline cort on day 12
#on nestling feather length (PC1) at day 12

all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)

all_nestling_cort_model <- 
  lmer(pc1~ scale(dl_base_cort)*species+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+
         (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model)

#keep nest ID and Site
all_nestling_cort_model_1 <- 
  lmer(pc1~ scale(dl_base_cort)*species+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+
         (1|site)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)
#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)

##best cort 1
best_cort<-lmer(pc1~ scale(dl_base_cort)+species
                +(1|site)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best_cort)

null<-lmer(pc1~1+(1|site)+(1|nest_ID),data=all_dat)
anova(best_cort,null, test="F")

confint(best_cort, level=0.85)
r.squaredGLMM(best_cort)

plot(allEffects(best_cort, confidence.level = 0.85),ylab="PCfeather", main="Nestling, all species")
#################################################
#########All SPECIES PC2 NESTLING CORT #################
#################################################
#datasheet 1_kleist_et_al

##this analysis explores the effects of nestlign baseline cort on day 12
#on nestling body size (tarsus and mass PC2) at day 12

all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)

all_nestling_cort_model <- lmer(pc2~ scale(dl_base_cort)*species+scale(laydate)+scale(distance)+scale(time)+scale(tree_cover_site)+ scale(chicks)+
                                  (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model)

#keep nestID

all_nestling_cort_model_1 <-  lmer(pc2~ scale(dl_base_cort)*species+scale(laydate)+scale(distance)+scale(time)+scale(tree_cover_site)+ scale(chicks)+
                                     (1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)
#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)

#best model
#no models with cort within 2 AIC
best<-lmer(pc2~scale(laydate)+species+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best)

null<-lmer(pc2~1+(1|nest_ID), data=all_dat)
anova(best,null, test="F")

confint(best, level=0.85)
r.squaredGLMM(best)

plot(allEffects(best, confidence.level = 0.85),ylab="Body size (PC2)", main="All species")


#######################ALL PC ACUTE CORT NESTLING ###################
#######################ALL PC ACUTE CORT NESTLING ###################
#######################ALL PC ACUTE CORT NESTLING ###################

#################################################
#########All SPECIES PC1 NESTLING ACUTE CORT #################
#################################################

#data sheet 1 kleist et all  (no value for acute cort on row 161, deleted for analysis)


all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)
##subset
atfl_dat <- subset(all_dat, species=="ATFL")
webl_dat <- subset(all_dat, species=="WEBL")
mobl_dat <- subset(all_dat, species=="MOBL")

plot(pc1~acute_cort, data=atfl_dat)
plot(pc1~acute_cort, data=webl_dat)
plot(pc1~acute_cort, data=mobl_dat)

##no evidence for nonlinear fit

all_nestling_cort_model <- 
  lmer(pc1~ scale(acute_cort)*species+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+
         (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model)

#keep nest ID
all_nestling_cort_model_1 <- 
  lmer(pc1~ scale(acute_cort)*species+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(chicks)+
         (1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)

summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)
#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)

##best cort 
best_cort<-lmer(pc1~ scale(acute_cort)+species
                +(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best_cort)

null<-lmer(pc1~1+(1|nest_ID), data=all_dat)
anova(best_cort,null, test="F")

confint(best_cort, level=0.85)
r.squaredGLMM(best_cort)

plot(allEffects(best_cort, confidence.level = 0.85),ylab="PCfeather", main="Nestling, all species")
#################################################
#########All SPECIES PC2 NESTLING ACUTE CORT #################
#################################################

#data sheet 1 kleist et al  (no value for acute cort on row 161, deleted for analysis)

all_dat$year<-as.character(all_dat$year)
all_dat$pair<-as.factor(all_dat$pair)
all_dat$site<-as.factor(all_dat$site)
all_dat$time<-as.numeric(all_dat$time)


all_nestling_cort_model <- lmer(pc2~ scale(acute_cort)*species+scale(laydate)+scale(distance)+scale(time)+scale(tree_cover_site)+ scale(chicks)+
                                  (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model)

#keep nestID

all_nestling_cort_model_1 <-  lmer(pc2~ scale(acute_cort)*species+scale(laydate)+scale(distance)+scale(time)+scale(tree_cover_site)+ scale(chicks)+
                                     (1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(all_nestling_cort_model_1)
####dredge on model with random effects removed
options(na.action = "na.fail")
all_dredge<- dredge(all_nestling_cort_model_1)
#subset of models within 2 AIC of best model
subset(all_dredge, delta < 2)

#best model with cort
best_cort<-lmer(pc2~scale(acute_cort)+species+(1|nest_ID), REML=FALSE, control=lmerControl(optimizer="bobyqa"), data=all_dat)
summary(best_cort)

confint(best_cort, level=0.85)
r.squaredGLMM(best_cort)

null<-lmer(pc2~1+(1|nest_ID),data=all_dat)
anova(best_cort,null, test="F")

plot(allEffects(best, confidence.level = 0.85),ylab="Body size (PC2)", main="All species")

############Incubation mom cort hatching success######################
############Incubation mom cort hatching success######################
############Incubation mom cort hatching success######################

##datasheet 2 kleist et al, must organize data set to remove rows
#without data for baseline cort measured during incubation lifestage 

hatch_dat$year<-as.character(hatch_dat$year)
hatch_dat$pair<-as.factor(hatch_dat$pair)
hatch_dat$site<-as.factor(hatch_dat$site)

##hatching success model
hatch_model <- 
  glmer(hatched~ scale(mom_base_cort)*species+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(eggs)+
          (1|pair)+(1|site)+(1|nestbox)+(1|year)+(1|nest_ID), data=hatch_dat, 
        family="binomial")

summary(hatch_model)

#remove all run glm

hatch_model_1 <- glm(hatched~ scale(mom_base_cort)*species+scale(laydate)+scale(distance)+scale(tree_cover_site)+ scale(eggs)
                     , data=hatch_dat, family="binomial")
summary(hatch_model_1)
####dredge on model with random effects removed except species

options(na.action = "na.fail")
hatch_dredge<- dredge(hatch_model_1)

#subset of models within 2 AIC of best model
subset(hatch_dredge, delta < 2)

best<-glm(hatched~ mom_base_cort+species,
          data=hatch_dat, family="binomial")
summary(best)
confint(best, level=0.85)
r.squaredGLMM(best)

null<-glm(hatched~1, data=hatch_dat, family="binomial")

anova(best,null, test="Chisq")

###plot ######

plot(allEffects(best, confidence.level = .85, xlevels=1))