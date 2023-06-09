#R code for the analyses on individual social network metrics in the manuscript "Resource distribution alters individual social network positions but not group network structures in experimental populations of forked fungus beetles"
#This script measures the effect of resource distribution on individual social network metrics (strength, betweenness, and clustering coefficient). This script first runs models with the observed dataset. Then, this script performs a node-level permutation to create 5000 randomizations of the observed dataset. This script then generates a null distribution of coefficients from models that use these randomized datasets. After generating a null distribution of coefficients, this script then compares the coefficients from the models that use the observed dataset and calculates p-values.
#Data are available at **add dryad link**

#set working directory to where you saved the data file
setwd("/Users/robincostello/Dropbox/Manuscripts/SN/For Resubmission")

#load libraries
library(lme4)
library(car)
library(glmmTMB)
library(DHARMa)
library(effects)
library(dplyr)
library(tidyr)
library(ggplot2)
library(foreach)
library(doParallel)
library(ggpubr)
library(ggeffects)
library(emmeans)
library(ggsignif)

#change R's default to contrasts 
options(contrasts=c("contr.sum", "contr.poly")) 

#upload data
IndividSN <- read.csv("IndividSN.csv")

#observed data models####

##function to standardize a variable 
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

##standardize variables
IndividSN %>%
  mutate(StnScansSeen = scale_this(ScansSeen)) %>%
  mutate(StnElytra = scale_this(Elytra)) -> IndividSN

##strength####
obs.strength.sex.0inflation <- glmmTMB(alpha~Treatment+Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
summary(obs.strength.sex.0inflation)
Anova(obs.strength.sex.0inflation, type=3)
plot(allEffects(obs.strength.sex.0inflation))

###test assumptions with DHARMa
obs.strength.sex.0inflation.resid<-simulateResiduals(obs.strength.sex.0inflation, n=250, integerResponse=F)
plot(obs.strength.sex.0inflation.resid)
testDispersion(obs.strength.sex.0inflation.resid)
testZeroInflation(obs.strength.sex.0inflation.resid)

###extract estimates
obs.strength.treatment.sex.coef.trt<-summary(obs.strength.sex.0inflation)$coefficients$cond[2,1]
obs.strength.treatment.sex.coef.sex<-summary(obs.strength.sex.0inflation)$coefficients$cond[3,1]
obs.strength.treatment.sex.coef.scansseen<-summary(obs.strength.sex.0inflation)$coefficients$cond[4,1]
obs.strength.treatment.sex.coef.elytra<-summary(obs.strength.sex.0inflation)$coefficients$cond[5,1]
obs.strength.treatment.sex.coef.period<-summary(obs.strength.sex.0inflation)$coefficients$cond[6,1]

###calculate marginal means
emmeans(obs.strength.sex.0inflation, "Treatment")
emmeans(obs.strength.sex.0inflation, "Period")
emmeans(obs.strength.sex.0inflation, "Survey_Sex")

###calculate effect size (Cohen's F2)
####r2full
r2full <- r2_nakagawa(obs.strength.sex.0inflation)
####r2treatment
obs.strength.sex.0inflation.notrt <- glmmTMB(alpha~Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2treatment <- r2_nakagawa(obs.strength.sex.0inflation.notrt)
####f2treament
f2treatment <- (r2full$R2_marginal - r2treatment$R2_marginal) / (1 - r2full$R2_marginal)
####r2sex
obs.strength.sex.0inflation.nosex <- glmmTMB(alpha~Treatment+StnScansSeen+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2sex <- r2_nakagawa(obs.strength.sex.0inflation.nosex)
####f2sex
f2sex <- (r2full$R2_marginal - r2sex$R2_marginal) / (1 - r2full$R2_marginal)
####r2elytra
obs.strength.sex.0inflation.noelytra <- glmmTMB(alpha~Treatment+Survey_Sex+StnScansSeen+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2elytra <- r2_nakagawa(obs.strength.sex.0inflation.noelytra)
####f2elytra
f2elytra <- (r2full$R2_marginal - r2elytra$R2_marginal) / (1 - r2full$R2_marginal)
####r2period
obs.strength.sex.0inflation.noperiod <- glmmTMB(alpha~Treatment+Survey_Sex+StnScansSeen+StnElytra+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2period <- r2_nakagawa(obs.strength.sex.0inflation.noperiod)
####f2period
f2period <- (r2full$R2_marginal - r2period$R2_marginal) / (1 - r2full$R2_marginal)
####r2scans
obs.strength.sex.0inflation.noscans <- glmmTMB(alpha~Treatment+Survey_Sex+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2scans <- r2_nakagawa(obs.strength.sex.0inflation.noscans)
####f2scans
f2scans <- (r2full$R2_marginal - r2scans$R2_marginal) / (1 - r2full$R2_marginal)

##betweenness####
obs.betweenness.sex.0inflation <- glmmTMB(betweenness~Treatment+Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), data=IndividSN, ziformula=~1)
summary(obs.betweenness.sex.0inflation)
Anova(obs.betweenness.sex.0inflation, type=3)
plot(allEffects(obs.betweenness.sex.0inflation))

###test for assumptions with DHARMa
obs.betweenness.sex.0inflation.resid<-simulateResiduals(obs.betweenness.sex.0inflation, n=250)
plot(obs.betweenness.sex.0inflation.resid)
testDispersion(obs.betweenness.sex.0inflation.resid)
testZeroInflation(obs.betweenness.sex.0inflation.resid)

###extract estimates
obs.betweenness.treatment.sex.coef.trt<-summary(obs.betweenness.sex.0inflation)$coefficients$cond[2,1]
obs.betweenness.treatment.sex.coef.sex<-summary(obs.betweenness.sex.0inflation)$coefficients$cond[3,1]
obs.betweenness.treatment.sex.coef.scansseen<-summary(obs.betweenness.sex.0inflation)$coefficients$cond[4,1]
obs.betweenness.treatment.sex.coef.elytra<-summary(obs.betweenness.sex.0inflation)$coefficients$cond[5,1]
obs.betweenness.treatment.sex.coef.period<-summary(obs.betweenness.sex.0inflation)$coefficients$cond[6,1]

###calculate marginal means
emmeans(obs.betweenness.sex.0inflation, "Treatment")
emmeans(obs.betweenness.sex.0inflation, "Period")
emmeans(obs.betweenness.sex.0inflation, "Survey_Sex")

###calculate effect size (Cohen's F2)
####r2full
r2full <- r2_nakagawa(obs.betweenness.sex.0inflation, tolerance = 1e-06)
####r2treatment
obs.betweenness.sex.0inflation.notrt <- glmmTMB(betweenness~Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2treatment <- r2_nakagawa(obs.betweenness.sex.0inflation.notrt, tolerance = 1e-06)
####f2treament
f2treatment <- (r2full$R2_marginal - r2treatment$R2_marginal) / (1 - r2full$R2_marginal)
####r2sex
obs.betweenness.sex.0inflation.nosex <- glmmTMB(betweenness~Treatment+StnScansSeen+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2sex <- r2_nakagawa(obs.betweenness.sex.0inflation.nosex, tolerance = 1e-06)
####f2sex
f2sex <- (r2full$R2_marginal - r2sex$R2_marginal) / (1 - r2full$R2_marginal)
####r2elytra
obs.betweenness.sex.0inflation.noelytra <- glmmTMB(betweenness~Treatment+Survey_Sex+StnScansSeen+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2elytra <- r2_nakagawa(obs.betweenness.sex.0inflation.noelytra, tolerance = 1e-06)
####f2elytra
f2elytra <- (r2full$R2_marginal - r2elytra$R2_marginal) / (1 - r2full$R2_marginal)
####r2period
obs.betweenness.sex.0inflation.noperiod <- glmmTMB(betweenness~Treatment+Survey_Sex+StnScansSeen+StnElytra+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2period <- r2_nakagawa(obs.betweenness.sex.0inflation.noperiod, tolerance = 1e-08)
####f2period
f2period <- (r2full$R2_marginal - r2period$R2_marginal) / (1 - r2full$R2_marginal)
####r2scans
obs.betweenness.sex.0inflation.noscans <- glmmTMB(betweenness~Treatment+Survey_Sex+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=IndividSN)
r2scans <- r2_nakagawa(obs.betweenness.sex.0inflation.noscans, tolerance = 1e-36)
####f2scans
f2scans <- (r2full$R2_marginal - r2scans$R2_marginal) / (1 - r2full$R2_marginal)

##clustering coefficient####
obs.cc.sex.0inflation<-glmmTMB(am~Treatment+Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), data=IndividSN, na.action=na.omit, ziformula=~1)
summary(obs.cc.sex.0inflation)
Anova(obs.cc.sex.0inflation, type=3)
plot(allEffects(obs.cc.sex.0inflation))

###test for assumptions with DHARMa
obs.cc.sex.0inflation.resid<-simulateResiduals(obs.cc.sex.0inflation, n=250, integerResponse=FALSE)
plot(obs.cc.sex.0inflation.resid)
testDispersion(obs.cc.sex.0inflation.resid)
testZeroInflation(obs.cc.sex.0inflation.resid)

###extract estimates
obs.cc.treatment.sex.coef.trt<-summary(obs.cc.sex.0inflation)$coefficients$cond[2,1]
obs.cc.treatment.sex.coef.sex<-summary(obs.cc.sex.0inflation)$coefficients$cond[3,1]
obs.cc.treatment.sex.coef.scansseen<-summary(obs.cc.sex.0inflation)$coefficients$cond[4,1]
obs.cc.treatment.sex.coef.elytra<-summary(obs.cc.sex.0inflation)$coefficients$cond[5,1]
obs.cc.treatment.sex.coef.period<-summary(obs.cc.sex.0inflation)$coefficients$cond[6,1]

###calculate marginal means
emmeans(obs.cc.sex.0inflation, "Treatment")
emmeans(obs.cc.sex.0inflation, "Period")
emmeans(obs.cc.sex.0inflation, "Survey_Sex")

###calculate effect size (Cohen's F2)
####r2full
r2full <- r2_nakagawa(obs.cc.sex.0inflation, tolerance = 1e-12)
####r2treatment
obs.cc.sex.0inflation.notrt <- glmmTMB(am~Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=IndividSN, na.action=na.omit)
r2treatment <- r2_nakagawa(obs.cc.sex.0inflation.notrt, tolerance = 1e-12)
####f2treament
f2treatment <- (r2full$R2_marginal - r2treatment$R2_marginal) / (1 - r2full$R2_marginal)
####r2sex
obs.cc.sex.0inflation.nosex <- glmmTMB(am~Treatment+StnScansSeen+StnElytra+Period+(1|Condo/ID), ziformula=~1, na.action=na.omit, data=IndividSN)
r2sex <- r2_nakagawa(obs.cc.sex.0inflation.nosex, tolerance = 1e-11)
####f2sex
f2sex <- (r2full$R2_marginal - r2sex$R2_marginal) / (1 - r2full$R2_marginal)
####r2elytra
obs.cc.sex.0inflation.noelytra <- glmmTMB(am~Treatment+Survey_Sex+StnScansSeen+Period+(1|Condo/ID), ziformula=~1, na.action=na.omit, data=IndividSN)
r2elytra <- r2_nakagawa(obs.cc.sex.0inflation.noelytra, tolerance = 1e-11)
####f2elytra
f2elytra <- (r2full$R2_marginal - r2elytra$R2_marginal) / (1 - r2full$R2_marginal)
####r2period
obs.cc.sex.0inflation.noperiod <- glmmTMB(am~Treatment+Survey_Sex+StnScansSeen+StnElytra+(1|Condo/ID), ziformula=~1, na.action=na.omit, data=IndividSN)
r2period <- r2_nakagawa(obs.cc.sex.0inflation.noperiod, tolerance = 1e-12)
####f2period
f2period <- (r2full$R2_marginal - r2period$R2_marginal) / (1 - r2full$R2_marginal)
####r2scans
obs.cc.sex.0inflation.noscans <- glmmTMB(am~Treatment+Survey_Sex+StnElytra+Period+(1|Condo/ID), ziformula=~1, na.action=na.omit, data=IndividSN)
r2scans <- r2_nakagawa(obs.cc.sex.0inflation.noscans, tolerance = 1e-11)
####f2scans
f2scans <- (r2full$R2_marginal - r2scans$R2_marginal) / (1 - r2full$R2_marginal)

#create permuted datasets####

##set number of permutations
permutations=10000

##shuffle observed data
for(i in 1:permutations) {
  PermIndividSN[[i]] <- IndividSN %>% mutate(StnElytra=sample(StnElytra)) %>% mutate(StnScansSeen=sample(StnScansSeen)) %>% mutate(Treatment=sample(Treatment)) %>% mutate(Survey_Sex=sample(Survey_Sex)) %>% mutate(Period=sample(Period))
}

##name each permutation sequentially
names(PermIndividSN) <- paste("perm", 1:length(PermIndividSN), sep="_")

#permuted data models####
##function that extracts the model summary
Model.Summary=function(x){
  setTimeLimit(cpu = Inf,  transient = FALSE)
  return(tryCatch(summary(x), error=function(e) NULL))
}

##strength####
###function that runs the model
perm.strength = function(PermIndividSN){
  perm.strength.mod<-glmmTMB(alpha~Treatment+Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), ziformula=~1, data=PermIndividSN)
  return(perm.strength.mod)
}

###run the model on each permutation
perm.strength.results<-lapply(PermIndividSN, FUN=perm.strength)

###extract model summary for each permutation
perm.strength.summary<-lapply(perm.strength.results, FUN=Model.Summary)

###remove models that do not converge (NA AIC values and NA P values)
####first extract AIC values from models
perm.strength.AIC<-foreach(i = 1:length(perm.strength.summary)) %dopar% {
  perm.strength.summary[[i]][["AICtab"]][["AIC"]]
}
names(perm.strength.AIC) <- paste("perm", 1:length(perm.strength.AIC), sep="_")
####make list of permutations with NA AIC values
perm.strength.AIC.NA<-perm.strength.AIC[is.na(perm.strength.AIC)]
####now, remove models that do not have AIC values
perm.strength.summary.converged<-perm.strength.summary[names(perm.strength.summary) %sans% names(perm.strength.AIC.NA)]
####because this doesn't account for all the models that did not converge, also need to remove models with NA P values. extra P values from models
perm.strength.p<-foreach(i = 1:length(perm.strength.summary.converged)) %dopar% {
  perm.strength.summary.converged[[i]][["coefficients"]][["cond"]][[11]]
}
names(perm.strength.p) <- names(perm.strength.summary.converged)
####make list of permutations with NA P values
perm.strength.p.NA<-perm.strength.p[is.na(perm.strength.p)]
####remove models that have NA P values
perm.strength.summary.converged<-perm.strength.summary.converged[names(perm.strength.summary.converged) %sans% names(perm.strength.p.NA)]

##select only first 5000 permutations
perm.strength.summary.converged %>%
  head(5000) -> perm.strength.summary.converged

##extract model estimate for treatment from model summary
perm.strength.treatment.coef<-foreach(i = 1:length(perm.strength.summary.converged)) %dopar% {
  perm.strength.summary.converged[[i]][["coefficients"]][["cond"]][[2]]
}
##flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.strength.treatment.coef.df<-do.call(rbind, perm.strength.treatment.coef)

##extract model estimate for sex from model summary
perm.strength.treatment.coef.sex<-foreach(i = 1:length(perm.strength.summary.converged)) %dopar% {
  perm.strength.summary.converged[[i]][["coefficients"]][["cond"]][[3]]
}
##flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.strength.treatment.coef.sex.df<-do.call(rbind, perm.strength.treatment.coef.sex)

##extract model estimate for scans seen from model summary
perm.strength.treatment.coef.scansseen<-foreach(i = 1:length(perm.strength.summary.converged)) %dopar% {
  perm.strength.summary.converged[[i]][["coefficients"]][["cond"]][[4]]
}
##flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.strength.treatment.coef.scansseen.df<-do.call(rbind, perm.strength.treatment.coef.scansseen)

##extract model estimate for elytra from model summary
perm.strength.treatment.coef.elytra<-foreach(i = 1:length(perm.strength.summary.converged)) %dopar% {
  perm.strength.summary.converged[[i]][["coefficients"]][["cond"]][[5]]
}
##flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.strength.treatment.coef.elytra.df<-do.call(rbind, perm.strength.treatment.coef.elytra)

##extract model estimate for period from model summary
perm.strength.treatment.coef.period<-foreach(i = 1:length(perm.strength.summary.converged)) %dopar% {
  perm.strength.summary.converged[[i]][["coefficients"]][["cond"]][[6]]
}
##flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.strength.treatment.coef.period.df<-do.call(rbind, perm.strength.treatment.coef.period)

##betweenness####
###function that runs the model
perm.betweenness = function(PermIndividSN){
  perm.betweenness.mod<-glmmTMB(betweenness~Treatment+Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), data=PermIndividSN, ziformula=~1)
  return(perm.betweenness.mod)
}

####run the model on each permutation
perm.betweenness.results<-lapply(PermIndividSN, FUN=perm.betweenness)

###extract model summary for each permutation
perm.betweenness.summary<-lapply(perm.betweenness.results, FUN=Model.Summary)

###remove models that do not converge (NA AIC values and NA P values)
####first extract AIC values from models
perm.betweenness.AIC<-foreach(i = 1:length(perm.betweenness.summary)) %dopar% {
  perm.betweenness.summary[[i]][["AICtab"]][["AIC"]]
}
names(perm.betweenness.AIC) <- paste("perm", 1:length(perm.betweenness.AIC), sep="_")
####make list of permutations with NA AIC values
perm.betweenness.AIC.NA<-perm.betweenness.AIC[is.na(perm.betweenness.AIC)]
####now, remove models that do not have AIC values
perm.betweenness.summary.converged<-perm.betweenness.summary[names(perm.betweenness.summary) %sans% names(perm.betweenness.AIC.NA)]
####because this doesn't account for all the models that did not converge, also need to remove models with NA P values. extra P values from models
perm.betweenness.p<-foreach(i = 1:length(perm.betweenness.summary.converged)) %dopar% {
  perm.betweenness.summary.converged[[i]][["coefficients"]][["cond"]][[11]]
}
names(perm.betweenness.p) <- names(perm.betweenness.summary.converged)
####make list of permutations with NA P values
perm.betweenness.p.NA<-perm.betweenness.p[is.na(perm.betweenness.p)]
####remove models that have NA P values
perm.betweenness.summary.converged<-perm.betweenness.summary.converged[names(perm.betweenness.summary.converged) %sans% names(perm.betweenness.p.NA)]

###select only first 5000 permutations
perm.betweenness.summary.converged %>%
  head(5000) -> perm.betweenness.summary.converged

###extract model estimate for treatment from model summary
perm.betweenness.treatment.coef<-foreach(i = 1:length(perm.betweenness.summary.converged)) %dopar% {
  perm.betweenness.summary[[i]][["coefficients"]][["cond"]][[2]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.betweenness.treatment.coef.df<-do.call(rbind, perm.betweenness.treatment.coef)

###extract model estimate for sex from model summary
perm.betweenness.treatment.coef.sex<-foreach(i = 1:length(perm.betweenness.summary.converged)) %dopar% {
  perm.betweenness.summary.converged[[i]][["coefficients"]][["cond"]][[3]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.betweenness.treatment.coef.sex.df<-do.call(rbind, perm.betweenness.treatment.coef.sex)

###extract model estimate for scans seen from model summary
perm.betweenness.treatment.coef.scansseen<-foreach(i = 1:length(perm.betweenness.summary.converged)) %dopar% {
  perm.betweenness.summary.converged[[i]][["coefficients"]][["cond"]][[4]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.betweenness.treatment.coef.scansseen.df<-do.call(rbind, perm.betweenness.treatment.coef.scansseen)

###extract model estimate for elytra from model summary
perm.betweenness.treatment.coef.elytra<-foreach(i = 1:length(perm.betweenness.summary.converged)) %dopar% {
  perm.betweenness.summary.converged[[i]][["coefficients"]][["cond"]][[5]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.betweenness.treatment.coef.elytra.df<-do.call(rbind, perm.betweenness.treatment.coef.elytra)

###extract model estimate for period from model summary
perm.betweenness.treatment.coef.period<-foreach(i = 1:length(perm.betweenness.summary.converged)) %dopar% {
  perm.betweenness.summary.converged[[i]][["coefficients"]][["cond"]][[6]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.betweenness.treatment.coef.period.df<-do.call(rbind, perm.betweenness.treatment.coef.period)

##clustering coefficient####
###function that runs the model
perm.cc = function(PermIndividSN){
  perm.cc.mod<-glmmTMB(am~Treatment+Survey_Sex+StnScansSeen+StnElytra+Period+(1|Condo/ID), data=PermIndividSN, na.action=na.omit, ziformula=~1)
  return(perm.cc.mod)
}

###run the model on each permutation
perm.cc.results<-lapply(PermIndividSN, FUN=perm.cc)

###extract model summary for each permutation
perm.cc.summary<-lapply(perm.cc.results, FUN=Model.Summary)

###remove models that do not converge (NA AIC values and NA P values)
####first extract AIC values from models
perm.cc.AIC<-foreach(i = 1:length(perm.cc.summary)) %dopar% {
  perm.cc.summary[[i]][["AICtab"]][["AIC"]]
}
names(perm.cc.AIC) <- paste("perm", 1:length(perm.cc.AIC), sep="_")
####make list of permutations with NA AIC values
perm.cc.AIC.NA<-perm.cc.AIC[is.na(perm.cc.AIC)]
####now, remove models that do not have AIC values
perm.cc.summary.converged<-perm.cc.summary[names(perm.cc.summary) %sans% names(perm.cc.AIC.NA)]
####because this doesn't account for all the models that did not converge, also need to remove models with NA P values. extra P values from models
perm.cc.p<-foreach(i = 1:length(perm.cc.summary.converged)) %dopar% {
  perm.cc.summary.converged[[i]][["coefficients"]][["cond"]][[11]]
}
names(perm.cc.p) <- names(perm.cc.summary.converged)
####make list of permutations with NA P values
perm.cc.p.NA<-perm.cc.p[is.na(perm.cc.p)]
####remove models that have NA P values
perm.cc.summary.converged<-perm.cc.summary.converged[names(perm.cc.summary.converged) %sans% names(perm.cc.p.NA)]

###select only first 5000 permutations
perm.cc.summary.converged %>%
  head(5000) -> perm.cc.summary.converged

###extract model estimate for treatment from model summary
perm.cc.treatment.coef<-foreach(i = 1:length(perm.cc.summary.converged)) %dopar% {
  perm.cc.summary.converged[[i]][["coefficients"]][["cond"]][[2]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.cc.treatment.coef.df<-do.call(rbind, perm.cc.treatment.coef)

###extract model estimate for sex from model summary
perm.cc.treatment.coef.sex<-foreach(i = 1:length(perm.cc.summary.converged)) %dopar% {
  perm.cc.summary.converged[[i]][["coefficients"]][["cond"]][[3]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.cc.treatment.coef.sex.df<-do.call(rbind, perm.cc.treatment.coef.sex)

###extract model estimate for scans seen from model summary
perm.cc.treatment.coef.scansseen<-foreach(i = 1:length(perm.cc.summary.converged)) %dopar% {
  perm.cc.summary.converged[[i]][["coefficients"]][["cond"]][[4]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.cc.treatment.coef.scansseen.df<-do.call(rbind, perm.cc.treatment.coef.scansseen)

###extract model estimate for elytra from model summary
perm.cc.treatment.coef.elytra<-foreach(i = 1:length(perm.cc.summary.converged)) %dopar% {
  perm.cc.summary.converged[[i]][["coefficients"]][["cond"]][[5]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.cc.treatment.coef.elytra.df<-do.call(rbind, perm.cc.treatment.coef.elytra)

###extract model estimate for period from model summary
perm.cc.treatment.coef.period<-foreach(i = 1:length(perm.cc.summary.converged)) %dopar% {
  perm.cc.summary.converged[[i]][["coefficients"]][["cond"]][[6]]
}
###flatten the list of permutations into a dataframe with one column that contains all the estimates from every permutation
perm.cc.treatment.coef.period.df<-do.call(rbind, perm.cc.treatment.coef.period)

#calculate p-values####
# FUNCTION: CWW (2015) Calculates one- and two-tailed p-values from a permuted distribution  
# Author Corlett Wolfe Wood

# One-tailed: The proportion of permuted values that are GREATER THAN the observed value
# Two-tailed: The proportion of permuted values that are MORE EXTREME than the observed

# NB: Calculation of two-tailed p-values on asymmetric null distributions
# On symmetric distributions, one-tailed p-values can be doubled to obtained two-tailed ones
# This doesn't work with asymmetric distributions because of skewness
# --> SOLUTION: This function below calculates two-tailed p-values as the proportion of
#     permuted.values that have a probability <= the probability of the observed value
# --> SOURCE: Pratt, J. W. and J. D. Gibbons. 1981. Concepts in Nonparametric Theory. 
#     Springer-Verlag, New York, NY pp. 29-32

# NB: p-values should NEVER be zero in a permutation test
# --> the observed data is one possible permutation, so at least one value
#     in the permuted data can be equal to the observed data (P>0)
# --> SOURCE: North, B. V., D. Curtis, and P. C. Sham. 2002. A note on the calculation of empirical 
#     p values from Monte Carlo procedures. Am. J. Hum. Genet. 71:439-441.
# --> SOURCE: Davison, A. C. and D. V. Hinkley. 1997. Bootstrap methods and their application. 
#     Cambridge University Press, Cambridge, United Kingdom
# --> SOLUTION: add 1 to the numerator and the denominator

pvals = function(observed.value, permuted.values){  
  # ONE-TAILED
  # Proportion of permuted values >= (if in upper tail) OR <= (if in lower tail) the observed value
  k.greater = length(which(permuted.values>=observed.value))
  k.less = length(which(permuted.values<=observed.value))
  n = length(permuted.values)
  # One-tailed p-value
  pval.onetail = (min(k.greater, k.less)+1)/(n+1) # See NB above for explanation of +1  
  
  # TWO-TAILED
  # Proportion of permuted values with PROBABILITIES <= the observed value
  df <- data.frame(permuted.values=permuted.values, prob=NA)
  names(df) <- c("permuted.values", "prob")
  for(p in 1:nrow(df)){
    # Calculate probability of each permuted value
    # Proportion of permuted values >= (if in upper tail) OR <= (if in lower tail) each permuted value
    k.greater = length(which(df$permuted.values>=df$permuted.values[p]))
    k.less = length(which(df$permuted.values<=df$permuted.values[p]))
    n = nrow(df)
    df$prob[p] = (min(k.greater, k.less)+1)/(n+1) # See NB above for explanation of +1
  }
  
  # Proportion of permuted.values with a have a probability <= the probability of the observed value
  prob.less <- length(which(df$prob<=pval.onetail))
  pval.twotail <- (prob.less+1)/(n+1) # See NB above for explanation of +1
  
  return(c(one.tailed=pval.onetail, 
           two.tailed=pval.twotail))
}

##strength####
strength.pvals<-pvals(obs.strength.treatment.sex.coef.trt, perm.strength.treatment.coef.df)
strength.pvals.sex<-pvals(obs.strength.treatment.sex.coef.sex, perm.strength.treatment.coef.sex.df)
strength.pvals.scansseen<-pvals(obs.strength.treatment.sex.coef.scansseen, perm.strength.treatment.coef.scansseen.df)
strength.pvals.elytra<-pvals(obs.strength.treatment.sex.coef.elytra, perm.strength.treatment.coef.elytra.df)
strength.pvals.period<-pvals(obs.strength.treatment.sex.coef.period, perm.strength.treatment.coef.period.df)

##betweenness####
betweenness.pvals<-pvals(obs.betweenness.treatment.sex.coef.trt, perm.betweenness.treatment.coef.df)
betweenness.pvals.sex<-pvals(obs.betweenness.treatment.sex.coef.sex, perm.betweenness.treatment.coef.sex.df)
betweenness.pvals.scansseen<-pvals(obs.betweenness.treatment.sex.coef.scansseen, perm.betweenness.treatment.coef.scansseen.df)
betweenness.pvals.elytra<-pvals(obs.betweenness.treatment.sex.coef.elytra, perm.betweenness.treatment.coef.elytra.df)
betweenness.pvals.period<-pvals(obs.betweenness.treatment.sex.coef.period, perm.betweenness.treatment.coef.period.df)

##clustering coefficient####
cc.pvals<-pvals(obs.cc.treatment.sex.coef.trt, perm.cc.treatment.coef.df)
cc.pvals.sex<-pvals(obs.cc.treatment.sex.coef.sex, perm.cc.treatment.coef.sex.df)
cc.pvals.scansseen<-pvals(obs.cc.treatment.sex.coef.scansseen, perm.cc.treatment.coef.scansseen.df)
cc.pvals.elytra<-pvals(obs.cc.treatment.sex.coef.elytra, perm.cc.treatment.coef.elytra.df)
cc.pvals.period<-pvals(obs.cc.treatment.sex.coef.period, perm.cc.treatment.coef.period.df)