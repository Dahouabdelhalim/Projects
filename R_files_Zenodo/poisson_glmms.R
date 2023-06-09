#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# Poisson GLMMs: predictors of individual microbiome richness 
rm(list=ls())
graphics.off()

library(MCMCglmm)
library(lme4)

setwd("~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/")

######################################################
#data
######################################################
##df for foraging data (N=16 individuals)
load("data/full_table_for_poisson_glmm.RData") #poisson_glmm

##df for social behavioral data (N=29 individuals)
load("data/diet_diversity_table_for_poisson_glmm.RData") #diet_richness_glmm

######################################################
# Diet GLMM
######################################################

part_tree_even <- glmer(Rarefied_OTUs ~  part_evenness + tree_evenness + (1 | Group), data = diet_richness_glmm, family = poisson(link = "log"), nAGQ = 100)
print(summary(part_tree_even)) 

######################################################
# Centrality and Age
######################################################

m1 <- glmer(Rarefied_OTUs ~  Weighted_Degree + Age + (1 | Group), data = poisson_glmm, family = poisson(link = "log"), nAGQ = 100)
print(summary(m1))

######################################################
# Centrality (Grooming Initiated vs Received)
######################################################

m2 <- glmer(Rarefied_OTUs ~  Weighted_In_Degree + Weighted_Out_Degree + (1 | Group), data = poisson_glmm, family = poisson(link = "log"), nAGQ = 100)
print(summary(m2))

######################################################
# Centrality (Adults Only) GLMM
######################################################

poisson_glmm_adult <- poisson_glmm[poisson_glmm$Age=="Adult"& !is.na(poisson_glmm$Centrality),]
poisson_glmm_adult <- droplevels(poisson_glmm_adult)

m3 <- glmer(Rarefied_OTUs ~ Weighted_Degree + (1 | Group), data = poisson_glmm_adult, family = poisson(link = "log"), nAGQ = 100)
print(summary(m3))

######################################################
# Scentmarking Rate GLMM
######################################################

poisson_glmm2 <- poisson_glmm[!is.na(poisson_glmm$Scentmark_Rate),]
poisson_glmm2$Individual <- droplevels(poisson_glmm2$Individual)
length(poisson_glmm2$Individual)
m4 <- glmer(Rarefied_OTUs ~  Scentmark_Rate + (1 | Group), data = poisson_glmm2, family = poisson(link = "log"), nAGQ = 100)
summary(m4)