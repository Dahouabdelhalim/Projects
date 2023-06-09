library("lme4")
library("emmeans")

dataAll <- read.csv("Fischer_et_al_Data_masks_lme_analysis.txt")
dataAll.frame <- data.frame(dataAll)
##########
# take only "mask on" data since we look at the changes!
dataUsed <- droplevels(subset(dataAll.frame, protocol == "mask" )) 

# d denotes delta changes, r denotes relative changes

## CBF
# LME model
mdl1.final<-lmer(rCBF ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## StO2
# LME model
mdl1.final<-lmer(dStO2 ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## tHb (THC)
# LME model
mdl1.final<-lmer(dTHC ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## OEF
# LME model
mdl1.final<-lmer(rOEF ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## CMRO2
# LME model
mdl1.final<-lmer(rCMRO2 ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## HR
# LME model
mdl1.final<-lmer(dHR_SPO2 ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## RR
# LME model
mdl1.final<-lmer(dRR ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## TcCO2 (TCO2)
# LME model
mdl1.final<-lmer(dTCO2 ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## PETCO2 (EtCO2)
# LME model
mdl1.final<-lmer(dEtCO2 ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## SpO2
# LME model
mdl1.final<-lmer(dSPO2 ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)

## MAP
# LME model
mdl1.final<-lmer(dMAPfina ~ 1  +  masktype  + (1|subjID) ,   data=dataUsed, REML=TRUE)
summary(mdl1.final)
# Post-hoc
emmdl.means  <- emmeans(mdl1.final, ~ masktype)

emmdl.means.contrast <- contrast(emmdl.means,adjust="fdr",method= list(
  "FFP2 "=c(1,0)-c(0,0),
  "surgical "=c(0,1)-c(0,0),
  "surgical-FFP2 "=c(0,1)-c(1,0)))

summary(emmdl.means.contrast)
confint(emmdl.means.contrast)
