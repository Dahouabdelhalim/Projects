# NSF FX Soil pH analysis

# This script runs the statistics for soil pH for the soils under the trees 
# Waiakea and Volcano from the Menge, Funk, Perakis, Wolf NSF grant 
# from 2015-2020.

rm(list=ls())

library(nlme)
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Load the soil pH data

V <- read.csv("HI_FX_soil_pH_Jul2019.csv")
W <- V[V$Species=="CAEQ" | V$Species=="GLSE" | V$Species=="PSCA",]
V <- V[V$Species=="ACKO" | V$Species=="MOFA" | V$Species=="DOVI",]

# The design is close to a basic anova

# Only use ones that are being used for the growth analysis (healthy).

###################################################################
########################### Soil pH ###############################
###################################################################

##### Waiakea 2019

# Get dataset for soil pH and for treatment and for tree
pH_HI_W_9 <- W[!is.na(W$Soil_pH) & W$Use_growth == 1,]$Soil_pH
Treatment_HI_W_9 <- W[!is.na(W$Soil_pH) & W$Use_growth == 1,]$Treatment
Species_HI_W_9 <- as.factor(W[!is.na(W$Soil_pH) & W$Use_growth == 1,]$Species)

summary(lm_soilpH_HI_W_3_06 <- lm(pH_HI_W_9 ~ 0 + Species_HI_W_9*Treatment_HI_W_9))
print(em_lm_soilpH_HI_W_3_06 <- emmeans(lm_soilpH_HI_W_3_06, list(pairwise ~ Species_HI_W_9*Treatment_HI_W_9), adjust = "tukey"))
# emmeans CH<GL,GH<GL,CM<GL,PM<GL
# b ab a ab; ab a a ab; ab a ab ab
tab <- as.data.frame(summary(em_lm_soilpH_HI_W_3_06)$emmean)
print(round(soilpH_GLSE_3_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4]),2))
print(round(soilpH_GLSE_3_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4]),2))
print(round(soilpH_GLSE_3_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4]),2))
print(round(soilpH_GLSE_3_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4]),2))
print(round(soilpH_CAEQ_3_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4]),2))
print(round(soilpH_CAEQ_3_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4]),2))
print(round(soilpH_CAEQ_3_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4]),2))
print(round(soilpH_CAEQ_3_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4]),2))
print(round(soilpH_PSCA_3_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4]),2))
print(round(soilpH_PSCA_3_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4]),2))
print(round(soilpH_PSCA_3_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4]),2))
print(round(soilpH_PSCA_3_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4]),2))




##### Volcano 2019

# Get dataset for soil KCl N and for treatment and for tree
pH_HI_V_9 <- V[!is.na(V$Soil_pH) & V$Use_growth == 1,]$Soil_pH
Treatment_HI_V_9 <- V[!is.na(V$Soil_pH) & V$Use_growth == 1,]$Treatment
Species_HI_V_9 <- as.factor(V[!is.na(V$Soil_pH) & V$Use_growth == 1,]$Species)

summary(lm_soilpH_HI_V_3_06 <- lm(pH_HI_V_9 ~ 0 + Species_HI_V_9*Treatment_HI_V_9))
print(em_lm_soilpH_HI_V_3_06 <- emmeans(lm_soilpH_HI_V_3_06, list(pairwise ~ Species_HI_V_9*Treatment_HI_V_9), adjust = "tukey"))
# emmeans none
# a a a a; a a a a; a a a a
tab <- as.data.frame(summary(em_lm_soilpH_HI_V_3_06)$emmean)
print(round(soilpH_ACKO_3_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4]),2))
print(round(soilpH_ACKO_3_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4]),2))
print(round(soilpH_ACKO_3_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4]),2))
print(round(soilpH_ACKO_3_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4]),2))
print(round(soilpH_MOFA_3_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4]),2))
print(round(soilpH_MOFA_3_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4]),2))
print(round(soilpH_MOFA_3_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4]),2))
print(round(soilpH_MOFA_3_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4]),2))
print(round(soilpH_DOVI_3_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4]),2))
print(round(soilpH_DOVI_3_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4]),2))
print(round(soilpH_DOVI_3_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4]),2))
print(round(soilpH_DOVI_3_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4]),2))






##########################################################################
##########################################################################
##########################################################################