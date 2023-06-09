# NSF FX Soil extractable N analysis

# This script runs the statistics for soil extractable (KCl) nitrate and 
# ammonium for the soils under the trees at Black Rock, Oregon, Waiakea, 
# and Volcano from the Menge, Funk, Perakis, Wolf NSF grant from 2015-2020.

rm(list=ls())

library(nlme)
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Load the soil %N data

V <- read.csv("HI_V_FX_Soil_KCl_N.csv")
W <- read.csv("HI_W_FX_Soil_KCl_N.csv")
N <- read.csv("NY_FX_Soil_KCl_N.csv")
O <- read.csv("OR_FX_Soil_KCl_N.csv")

# The design is close to a basic anova with repeated measures, so we used
# LMER across years with tree as a random effect.

# Only use ones that are being used for the growth analysis.

###################################################################
########################### Soil KCl N ############################
###################################################################

##### New York 2019

# Get dataset for soil KCl N and for treatment and for tree
NH4_NY_9 <- N[!is.na(N$NH4_ug_N.g_dry_soil) & N$Use_growth == 1,]$NH4_ug_N.g_dry_soil
NO3_NY_9 <- N[!is.na(N$NO3_ug_N.g_dry_soil) & N$Use_growth == 1,]$NO3_ug_N.g_dry_soil
NH4NO3_NY_9 <- N[!is.na(N$NH4NO3_ug_N.g_dry_soil) & N$Use_growth == 1,]$NH4NO3_ug_N.g_dry_soil
H2O_NY_9 <- N[!is.na(N$moisture_g_H2O.g_dry_soil) & N$Use_growth == 1,]$moisture_g_H2O.g_dry_soil
lnNH4_NY_9 <- log(NH4_NY_9)
# Convert 0s to 0.1s for log analysis
pNO3_NY_9 <- NO3_NY_9
pNO3_NY_9[pNO3_NY_9==0] <- 0.1
lnNO3_NY_9 <- log(pNO3_NY_9)
lnNH4NO3_NY_9 <- log(NH4NO3_NY_9)
Treatment_NY_9 <- N[!is.na(N$NH4NO3_ug_N.g_dry_soil) & N$Use_growth == 1,]$Treatment
Species_NY_9 <- as.factor(N[!is.na(N$NH4NO3_ug_N.g_dry_soil) & N$Use_growth == 1,]$Species)

# Soil NH4, NO3, and NH4+NO3 are highly heteroscedastic when not 
# log-transformed, so log-transforming for statistical analysis

# NH4
summary(lm_lnsoilNH4_NY_3_08 <- lm(lnNH4_NY_9 ~ 0 + Species_NY_9*Treatment_NY_9))
print(em_lm_lnsoilNH4_NY_3_08 <- emmeans(lm_lnsoilNH4_NY_3_08, list(pairwise ~ Species_NY_9*Treatment_NY_9), adjust = "tukey"))
# emmeans BL<BH, BL<BP, BL<RH, BL<RP, RL<RP
# ab abc bc c; a abc bc bc
tab <- as.data.frame(summary(em_lm_lnsoilNH4_NY_3_08)$emmean)
print(round(exp(soilNH4_ROPS_3_08_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNH4_ROPS_3_08_MN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNH4_ROPS_3_08_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNH4_ROPS_3_08_PHN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNH4_BENI_3_08_LN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNH4_BENI_3_08_MN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNH4_BENI_3_08_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNH4_BENI_3_08_PHN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))

# NO3

summary(lm_lnsoilNO3_NY_3_08 <- lm(lnNO3_NY_9 ~ 0 + Species_NY_9*Treatment_NY_9))
print(em_lm_lnsoilNO3_NY_3_08 <- emmeans(lm_lnsoilNO3_NY_3_08, list(pairwise ~ Species_NY_9*Treatment_NY_9), adjust = "tukey"))
# emmeans BL<BM, BL<BH, BL<BP, BL<RM, BL<RH, BL<RP, RL<BM, RL<BH, RL<BP, RL<RM, RL<RH, RL<RP
# a b b b; a b b b 
tab <- as.data.frame(summary(em_lm_lnsoilNO3_NY_3_08)$emmean)
print(round(exp(soilNO3_ROPS_3_08_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),2))
print(round(exp(soilNO3_ROPS_3_08_MN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNO3_ROPS_3_08_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNO3_ROPS_3_08_PHN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNO3_BENI_3_08_LN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),2))
print(round(exp(soilNO3_BENI_3_08_MN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNO3_BENI_3_08_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNO3_BENI_3_08_PHN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))

# NH4 + NO3
summary(lm_lnsoilNH4NO3_NY_3_08 <- lm(lnNH4NO3_NY_9 ~ 0 + Species_NY_9*Treatment_NY_9))
print(em_lm_lnsoilNH4NO3_NY_3_08 <- emmeans(lm_lnsoilNH4NO3_NY_3_08, list(pairwise ~ Species_NY_9*Treatment_NY_9), adjust = "tukey"))
# emmeans BL<BH, BL<BP, BL<RM, BL<RH, BL<RP, RL<RH, RL<RP
# ab bc c c; a abc bc bc
tab <- as.data.frame(summary(em_lm_lnsoilNH4NO3_NY_3_08)$emmean)
print(round(exp(soilNH4NO3_ROPS_3_08_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNH4NO3_ROPS_3_08_MN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNH4NO3_ROPS_3_08_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNH4NO3_ROPS_3_08_PHN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNH4NO3_BENI_3_08_LN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNH4NO3_BENI_3_08_MN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNH4NO3_BENI_3_08_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNH4NO3_BENI_3_08_PHN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))





##### Oregon 2018-2019

# Get dataset for soil KCl N and for treatment and for tree
NH4_OR_89 <- O[!is.na(O$NH4_ug_N.g_dry_soil) & O$Use_growth == 1,]$NH4_ug_N.g_dry_soil
NO3_OR_89 <- O[!is.na(O$NO3_ug_N.g_dry_soil) & O$Use_growth == 1,]$NO3_ug_N.g_dry_soil
NH4NO3_OR_89 <- O[!is.na(O$NH4NO3_ug_N.g_dry_soil) & O$Use_growth == 1,]$NH4NO3_ug_N.g_dry_soil
H2O_OR_89 <- O[!is.na(O$moisture_g_H2O.g_dry_soil) & O$Use_growth == 1,]$moisture_g_H2O.g_dry_soil
lnNH4_OR_89 <- log(NH4_OR_89)
# Convert 0s to 0.1s for log analysis
pNO3_OR_89 <- NO3_OR_89
pNO3_OR_89[pNO3_OR_89==0] <- 0.1
lnNO3_OR_89 <- log(pNO3_OR_89)
lnNH4NO3_OR_89 <- log(NH4NO3_OR_89)
Treatment_OR_89 <- O[!is.na(O$NH4NO3_ug_N.g_dry_soil) & O$Use_growth == 1,]$Treatment
Tree_OR_89 <- as.factor(O[!is.na(O$NH4NO3_ug_N.g_dry_soil) & O$Use_growth == 1,]$PID)
Species_OR_89 <- as.factor(O[!is.na(O$NH4NO3_ug_N.g_dry_soil) & O$Use_growth == 1,]$Species)

# Soil NH4, NO3, and NH4+NO3 are highly heteroscedastic when not 
# log-transformed, so log-transforming for statistical analysis

summary(lme_lnsoilNH4_OR_3_04_05 <- lme(lnNH4_OR_89 ~ 0 + Species_OR_89*Treatment_OR_89,
	random=~1 | Tree_OR_89))
print(em_lme_lnsoilNH4_OR_3_04_05 <- emmeans(lme_lnsoilNH4_OR_3_04_05, list(pairwise ~ Species_OR_89*Treatment_OR_89), adjust = "tukey"))
# emmeans PL<PM,PL<PH,PL<PP, PL<AM,PL<AH,PL<AP, AL<AM,AL<AH,AL<AP, AL<PM,AL<PH,AL<PP
# a b b b; a b b b 
tab <- as.data.frame(summary(em_lme_lnsoilNH4_OR_3_04_05)$emmean)
print(round(exp(soilNH4_ALRU_3_04_05_LN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNH4_ALRU_3_04_05_MN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNH4_ALRU_3_04_05_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNH4_ALRU_3_04_05_PHN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))
print(round(exp(soilNH4_PSME_3_04_05_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNH4_PSME_3_04_05_MN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNH4_PSME_3_04_05_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNH4_PSME_3_04_05_PHN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))

summary(lme_lnsoilNO3_OR_3_04_05 <- lme(lnNO3_OR_89 ~ 0 + Species_OR_89*Treatment_OR_89,
	random=~1 | Tree_OR_89))
print(em_lme_lnsoilNO3_OR_3_04_05 <- emmeans(lme_lnsoilNO3_OR_3_04_05, list(pairwise ~ Species_OR_89*Treatment_OR_89), adjust = "tukey"))
# emmeans PL<PM,PL<PH,PL<PP, PL<AM,PL<AH,PL<AP, AL<AM,AL<AH,AL<AP, AL<PM,AL<PH,AL<PP
# a b b b; a b b b 
tab <- as.data.frame(summary(em_lme_lnsoilNO3_OR_3_04_05)$emmean)
print(round(exp(soilNO3_ALRU_3_04_05_LN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNO3_ALRU_3_04_05_MN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNO3_ALRU_3_04_05_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNO3_ALRU_3_04_05_PHN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))
print(round(exp(soilNO3_PSME_3_04_05_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNO3_PSME_3_04_05_MN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNO3_PSME_3_04_05_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNO3_PSME_3_04_05_PHN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))

summary(lme_lnsoilNH4NO3_OR_3_04_05 <- lme(lnNH4NO3_OR_89 ~ 0 + Species_OR_89*Treatment_OR_89,
	random=~1 | Tree_OR_89))
print(em_lme_lnsoilNH4NO3_OR_3_04_05 <- emmeans(lme_lnsoilNH4NO3_OR_3_04_05, list(pairwise ~ Species_OR_89*Treatment_OR_89), adjust = "tukey"))
# emmeans PL<PM,PL<PH,PL<PP, PL<AM,PL<AH,PL<AP, AL<AM,AL<AH,AL<AP, AL<PM,AL<PH,AL<PP
# a b b b; a b b b 
tab <- as.data.frame(summary(em_lme_lnsoilNH4NO3_OR_3_04_05)$emmean)
print(round(exp(soilNH4NO3_ALRU_3_04_05_LN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNH4NO3_ALRU_3_04_05_MN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNH4NO3_ALRU_3_04_05_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNH4NO3_ALRU_3_04_05_PHN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))
print(round(exp(soilNH4NO3_PSME_3_04_05_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNH4NO3_PSME_3_04_05_MN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNH4NO3_PSME_3_04_05_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNH4NO3_PSME_3_04_05_PHN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))



##### Waiakea 2019

# Get dataset for soil KCl N and for treatment and for tree
NH4_HI_W_9 <- W[!is.na(W$NH4_ug_N.g_dry_soil) & W$Use_growth == 1,]$NH4_ug_N.g_dry_soil
NO3_HI_W_9 <- W[!is.na(W$NO3_ug_N.g_dry_soil) & W$Use_growth == 1,]$NO3_ug_N.g_dry_soil
NH4NO3_HI_W_9 <- W[!is.na(W$NH4NO3_ug_N.g_dry_soil) & W$Use_growth == 1,]$NH4NO3_ug_N.g_dry_soil
H2O_HI_W_9 <- W[!is.na(W$moisture_g_H2O.g_dry_soil) & W$Use_growth == 1,]$moisture_g_H2O.g_dry_soil
lnNH4_HI_W_9 <- log(NH4_HI_W_9)
# Convert 0s to 0.1s for log analysis
pNO3_HI_W_9 <- NO3_HI_W_9
pNO3_HI_W_9[pNO3_HI_W_9==0] <- 0.1
lnNO3_HI_W_9 <- log(pNO3_HI_W_9)
lnNH4NO3_HI_W_9 <- log(NH4NO3_HI_W_9)
Treatment_HI_W_9 <- W[!is.na(W$NH4NO3_ug_N.g_dry_soil) & W$Use_growth == 1,]$Treatment
Species_HI_W_9 <- as.factor(W[!is.na(W$NH4NO3_ug_N.g_dry_soil) & W$Use_growth == 1,]$Species)

# Soil NH4, NO3, and NH4+NO3 are highly heteroscedastic when not 
# log-transformed, so log-transforming for statistical analysis

summary(lm_lnsoilNH4_HI_W_3_06 <- lm(lnNH4_HI_W_9 ~ 0 + Species_HI_W_9*Treatment_HI_W_9))
print(em_lm_lnsoilNH4_HI_W_3_06 <- emmeans(lm_lnsoilNH4_HI_W_3_06, list(pairwise ~ Species_HI_W_9*Treatment_HI_W_9), adjust = "tukey"))
# emmeans none
# a a a a; a a a a; a a a a
tab <- as.data.frame(summary(em_lm_lnsoilNH4_HI_W_3_06)$emmean)
print(round(exp(soilNH4_GLSE_3_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNH4_GLSE_3_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNH4_GLSE_3_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNH4_GLSE_3_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4])),1))
print(round(exp(soilNH4_CAEQ_3_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNH4_CAEQ_3_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))
print(round(exp(soilNH4_CAEQ_3_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNH4_CAEQ_3_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4])),1))
print(round(exp(soilNH4_PSCA_3_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNH4_PSCA_3_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4])),1))
print(round(exp(soilNH4_PSCA_3_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNH4_PSCA_3_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4])),1))

summary(lm_lnsoilNO3_HI_W_3_06 <- lm(lnNO3_HI_W_9 ~ 0 + Species_HI_W_9*Treatment_HI_W_9))
print(em_lm_lnsoilNO3_HI_W_3_06 <- emmeans(lm_lnsoilNO3_HI_W_3_06, list(pairwise ~ Species_HI_W_9*Treatment_HI_W_9), adjust = "tukey"))
# emmeans none
# a a a a; a a a a; a a a a
tab <- as.data.frame(summary(em_lm_lnsoilNO3_HI_W_3_06)$emmean)
print(round(exp(soilNO3_GLSE_3_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),2))
print(round(exp(soilNO3_GLSE_3_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNO3_GLSE_3_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNO3_GLSE_3_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4])),1))
print(round(exp(soilNO3_CAEQ_3_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),2))
print(round(exp(soilNO3_CAEQ_3_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),2))
print(round(exp(soilNO3_CAEQ_3_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),2))
print(round(exp(soilNO3_CAEQ_3_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4])),1))
print(round(exp(soilNO3_PSCA_3_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),2))
print(round(exp(soilNO3_PSCA_3_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4])),1))
print(round(exp(soilNO3_PSCA_3_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),2))
print(round(exp(soilNO3_PSCA_3_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4])),1))

summary(lm_lnsoilNH4NO3_HI_W_3_06 <- lm(lnNH4NO3_HI_W_9 ~ 0 + Species_HI_W_9*Treatment_HI_W_9))
print(em_lm_lnsoilNH4NO3_HI_W_3_06 <- emmeans(lm_lnsoilNH4NO3_HI_W_3_06, list(pairwise ~ Species_HI_W_9*Treatment_HI_W_9), adjust = "tukey"))
# emmeans none
# a a a a; a a a a; a a a a
tab <- as.data.frame(summary(em_lm_lnsoilNH4NO3_HI_W_3_06)$emmean)
print(round(exp(soilNH4NO3_GLSE_3_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNH4NO3_GLSE_3_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNH4NO3_GLSE_3_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNH4NO3_GLSE_3_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4])),1))
print(round(exp(soilNH4NO3_CAEQ_3_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNH4NO3_CAEQ_3_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))
print(round(exp(soilNH4NO3_CAEQ_3_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNH4NO3_CAEQ_3_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4])),1))
print(round(exp(soilNH4NO3_PSCA_3_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNH4NO3_PSCA_3_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4])),1))
print(round(exp(soilNH4NO3_PSCA_3_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNH4NO3_PSCA_3_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4])),1))




##### Volcano 2019

# Get dataset for soil KCl N and for treatment and for tree
NH4_HI_V_9 <- V[!is.na(V$NH4_ug_N.g_dry_soil) & V$Use_growth == 1,]$NH4_ug_N.g_dry_soil
NO3_HI_V_9 <- V[!is.na(V$NO3_ug_N.g_dry_soil) & V$Use_growth == 1,]$NO3_ug_N.g_dry_soil
NH4NO3_HI_V_9 <- V[!is.na(V$NH4NO3_ug_N.g_dry_soil) & V$Use_growth == 1,]$NH4NO3_ug_N.g_dry_soil
H2O_HI_V_9 <- V[!is.na(V$moisture_g_H2O.g_dry_soil) & V$Use_growth == 1,]$moisture_g_H2O.g_dry_soil
lnNH4_HI_V_9 <- log(NH4_HI_V_9)
# Convert 0s to 0.1s for log analysis
pNO3_HI_V_9 <- NO3_HI_V_9
pNO3_HI_V_9[pNO3_HI_V_9==0] <- 0.1
lnNO3_HI_V_9 <- log(pNO3_HI_V_9)
lnNH4NO3_HI_V_9 <- log(NH4NO3_HI_V_9)
Treatment_HI_V_9 <- V[!is.na(V$NH4NO3_ug_N.g_dry_soil) & V$Use_growth == 1,]$Treatment
Species_HI_V_9 <- as.factor(V[!is.na(V$NH4NO3_ug_N.g_dry_soil) & V$Use_growth == 1,]$Species)

# Soil NH4, NO3, and NH4+NO3 are highly heteroscedastic when not 
# log-transformed, so log-transforming for statistical analysis

summary(lm_lnsoilNH4_HI_V_3_06 <- lm(lnNH4_HI_V_9 ~ 0 + Species_HI_V_9*Treatment_HI_V_9))
print(em_lm_lnsoilNH4_HI_V_3_06 <- emmeans(lm_lnsoilNH4_HI_V_3_06, list(pairwise ~ Species_HI_V_9*Treatment_HI_V_9), adjust = "tukey"))
# emmeans none
# a a a a; a a a a; a a a a
tab <- as.data.frame(summary(em_lm_lnsoilNH4_HI_V_3_06)$emmean)
print(round(exp(soilNH4_ACKO_3_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNH4_ACKO_3_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))
print(round(exp(soilNH4_ACKO_3_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNH4_ACKO_3_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4])),1))
print(round(exp(soilNH4_MOFA_3_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNH4_MOFA_3_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4])),1))
print(round(exp(soilNH4_MOFA_3_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNH4_MOFA_3_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4])),1))
print(round(exp(soilNH4_DOVI_3_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNH4_DOVI_3_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNH4_DOVI_3_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNH4_DOVI_3_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4])),1))

summary(lm_lnsoilNO3_HI_V_3_06 <- lm(lnNO3_HI_V_9 ~ 0 + Species_HI_V_9*Treatment_HI_V_9))
print(em_lm_lnsoilNO3_HI_V_3_06 <- emmeans(lm_lnsoilNO3_HI_V_3_06, list(pairwise ~ Species_HI_V_9*Treatment_HI_V_9), adjust = "tukey"))
# emmeans AL<AM,AL<AH,AL<AP, AL<ML,AL<MH,AL<MP, AL<DH,AL<DP
# a b b b; b ab b b; ab ab b b
tab <- as.data.frame(summary(em_lm_lnsoilNO3_HI_V_3_06)$emmean)
print(round(exp(soilNO3_ACKO_3_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),2))
print(round(exp(soilNO3_ACKO_3_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))
print(round(exp(soilNO3_ACKO_3_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNO3_ACKO_3_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4])),1))
print(round(exp(soilNO3_MOFA_3_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNO3_MOFA_3_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4])),1))
print(round(exp(soilNO3_MOFA_3_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNO3_MOFA_3_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4])),1))
print(round(exp(soilNO3_DOVI_3_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNO3_DOVI_3_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNO3_DOVI_3_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNO3_DOVI_3_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4])),1))

summary(lm_lnsoilNH4NO3_HI_V_3_06 <- lm(lnNH4NO3_HI_V_9 ~ 0 + Species_HI_V_9*Treatment_HI_V_9))
print(em_lm_lnsoilNH4NO3_HI_V_3_06 <- emmeans(lm_lnsoilNH4NO3_HI_V_3_06, list(pairwise ~ Species_HI_V_9*Treatment_HI_V_9), adjust = "tukey"))
# emmeans AL<MH
# a ab ab ab; ab ab b ab; ab ab ab ab 
tab <- as.data.frame(summary(em_lm_lnsoilNH4NO3_HI_V_3_06)$emmean)
print(round(exp(soilNH4NO3_ACKO_3_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4])),1))
print(round(exp(soilNH4NO3_ACKO_3_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4])),1))
print(round(exp(soilNH4NO3_ACKO_3_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4])),1))
print(round(exp(soilNH4NO3_ACKO_3_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4])),1))
print(round(exp(soilNH4NO3_MOFA_3_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4])),1))
print(round(exp(soilNH4NO3_MOFA_3_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4])),1))
print(round(exp(soilNH4NO3_MOFA_3_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4])),1))
print(round(exp(soilNH4NO3_MOFA_3_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4])),1))
print(round(exp(soilNH4NO3_DOVI_3_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4])),1))
print(round(exp(soilNH4NO3_DOVI_3_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4])),1))
print(round(exp(soilNH4NO3_DOVI_3_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4])),1))
print(round(exp(soilNH4NO3_DOVI_3_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4])),1))




##########################################################################
##########################################################################
##########################################################################