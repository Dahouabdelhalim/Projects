#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# GLMMs assessing pairwise social and genetic predictors of microbial dissimilarity 

rm(list=ls())
graphics.off()

library(ggplot2)
library(lme4)
library(lmtest)
library(glmmADMB)
library(MuMIn)
library(MCMCglmm)

focal_dir <- "~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/"
setwd(focal_dir)
getwd()

load("data/beta_groom_same_12mo.RData") #betareg_df_groom_same
head(betareg_df_groom_same)
nrow(betareg_df_groom_same) #167
betareg_df_groom_same$Group <- as.factor(betareg_df_groom_same$Group)
#####################################################
# mcmcglmm
#####################################################

# make sure Ind1 and Ind2 have same number of levels

l <- list(betareg_df_groom_same$Ind1, betareg_df_groom_same$Ind2)
l <- unlist(l)
l <- unique(l)
Ind1 <- l
Ind2 <- l
mult.memb(~Ind1+Ind2)

betareg_df_groom_same$Ind1 <- factor(betareg_df_groom_same$Ind1, levels=l)
betareg_df_groom_same$Ind2 <- factor(betareg_df_groom_same$Ind2, levels=l)
setdiff(levels(betareg_df_groom_same$Ind1),levels(betareg_df_groom_same$Ind2))
setdiff(levels(betareg_df_groom_same$Ind2),levels(betareg_df_groom_same$Ind1))

model.update <- updateable(MCMCglmm)
set.seed(123)

#####################################################
## Bray-Curtis
#####################################################
# group only 
global.model.BC.group <- model.update(BC ~ Group + Related , random=~mm(Ind1 + Ind2), data=betareg_df_groom_same, nitt=300000, burnin=25000, thin=15, pr=TRUE)
print(summary(global.model.BC.group), digits=4)

# PL only
global.model.BC.PL <- model.update(BC ~ PL + Related, random=~mm(Ind1 + Ind2), data=betareg_df_groom_same, nitt=300000, burnin=25000, thin=15, pr=TRUE)
print(summary(global.model.BC.PL), digits=3)

# group and PL
global.model.BC.PL <- model.update(BC ~ Group + PL + Related, random=~mm(Ind1 + Ind2), data=betareg_df_groom_same, nitt=300000, burnin=25000, pr=TRUE)
print(summary(global.model.BC.PL), digits=3)

# Model comparison
dredge.MCMCglmm.BC.both <- dredge(global.model.BC.both, rank="DIC", subset= "Group"|"PL") 
dredge.BC <- as.data.frame(dredge.MCMCglmm.BC.both)

# only pairs in the same group
global.model.BC.PL.within <- model.update(BC ~ PL + Related , random=~mm(Ind1 + Ind2), data=betareg_df_groom_within, nitt=300000, burnin=25000, pr=TRUE)
summary(global.model.BC.PL.within)

# pairs in different groups
global.model.BC.PL.different <- model.update(BC ~ PL + Related, random=~mm(Ind1 + Ind2), data=betareg_df_groom_different, nitt=300000, burnin=25000, pr=TRUE)
summary(global.model.BC.PL.different)

#####################################################
## weighted unifrac
#####################################################
# group only 
global.model.WU.group <- model.update(WU ~ Group + Related , random=~mm(Ind1 + Ind2), data=betareg_df_groom_same, nitt=300000, burnin=25000, thin=15, pr=TRUE)
summary(global.model.WU.group)

# path length only
global.model.WU.PL <- model.update(WU ~ PL + Related, random=~mm(Ind1 + Ind2), data=betareg_df_groom_same, nitt=300000, burnin=25000, thin=15, pr=TRUE)
summary(global.model.WU.PL)

## both PL and group
global.model.WU.both <- model.update(WU ~ Group + PL + Related, random=~mm(Ind1 + Ind2), data=betareg_df_groom_same, nitt=300000, burnin=25000, thin=15, pr=TRUE)
summary(global.model.WU.both)

#### Model compairson
dredge.MCMCglmm.WU.both <- dredge(global.model.WU.both, rank="DIC", subset= "Group"|"PL") 

best <- get.models(dredge.MCMCglmm.WU.both, 1)[[1]] #summary of best model
summary(best)
confset.d4 <- get.models(dredge.MCMCglmm.WU.both, subset = DIC<4)
summary(confset.d4)

dredge.WU <- as.data.frame(dredge.MCMCglmm.WU.both)

## different groups only 
betareg_df_groom_different <- betareg_df_groom_same[betareg_df_groom_same$Group=="Different",]
betareg_df_groom_different <- droplevels(betareg_df_groom_different)
betareg_df_groom_different$Ind1 <- factor(betareg_df_groom_different$Ind1, levels=l)
betareg_df_groom_different$Ind2 <- factor(betareg_df_groom_different$Ind2, levels=l)

global.model.WU.PL.within <- model.update(WU ~ PL, random=~mm(Ind1 + Ind2), data=betareg_df_groom_within, nitt=300000, burnin=25000, thin=15,pr=TRUE)
print(summary(global.model.WU.PL.within), digits=4)

## within groups only
betareg_df_groom_within <- betareg_df_groom_same[betareg_df_groom_same$Group=="Same",]
betareg_df_groom_within <- droplevels(betareg_df_groom_within)
betareg_df_groom_within$Ind1 <- factor(betareg_df_groom_within$Ind1, levels=l)
betareg_df_groom_within$Ind2 <- factor(betareg_df_groom_within$Ind2, levels=l)

global.model.WU.PL.different <- model.update(WU ~ PL, random=~mm(Ind1 + Ind2), data=betareg_df_groom_different, nitt=300000, burnin=25000, thin=15,pr=TRUE)
print(summary(global.model.WU.PL.different), digits=3)

#####################################################
#glmmADMB (allows for Beta-distributed response variable)

fit_beta_groom_group_BC <- glmmadmb(BC ~ Group + Related + (1|Ind1) + (1|Ind2), data=betareg_df_groom_same, zeroInflation=FALSE, family="beta", link="logit")
summary(fit_beta_groom_group_BC)

fit_beta_groom_group_WU <- glmmadmb(WU ~ Group + Related + PL + (1|Ind1) + (1|Ind2), data=betareg_df_groom_same, zeroInflation=FALSE, family="beta", link="logit")
summary(fit_beta_groom_group_WU)

fit_beta_groom_PL_BC <- glmmadmb(BC ~ PL + Related + (1|Ind1) + (1|Ind2), data=betareg_df_groom_same, zeroInflation=FALSE, family="beta", link="logit")
summary(fit_beta_groom_PL_BC)

fit_beta_groom_PL_WU <- glmmadmb(WU ~ PL + Related + (1|Ind1) + (1|Ind2), data=betareg_df_groom_same, zeroInflation=FALSE, family="beta", link="logit")
summary(fit_beta_groom_PL_WU)

fit_beta_groom_PL_group_BC <- glmmadmb(BC ~ PL + Group + Related + (1|Ind1) + (1|Ind2), data=betareg_df_groom_same, zeroInflation=FALSE, family="beta", link="logit")
summary(fit_beta_groom_PL_group_BC)

fit_beta_groom_PL_group_WU <- glmmadmb(WU ~ PL + Group + Related + (1|Ind1) + (1|Ind2), data=betareg_df_groom_same, zeroInflation=FALSE, family="beta", link="logit")
summary(fit_beta_groom_PL_group_WU)
