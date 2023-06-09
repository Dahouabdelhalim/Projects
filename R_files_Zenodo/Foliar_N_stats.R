# NSF FX Foliar N analysis

# This script runs the statistics for foliar %N for the trees at 
# Black Rock, Oregon, Waiakea, and Volcano from the Menge, Funk, Perakis, 
# Wolf NSF grant from 2015-2020.

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Load the foliar %N code

V <- read.csv("HI_V_FX_Size_FoliarCNIsotope_Data.csv")[1:96,1:111]
W <- read.csv("HI_W_FX_Size_FoliarCNIsotope_Data.csv")[1:108,1:111]
N <- read.csv("NY_FX_Size_FoliarCNIsotope_Data.csv")[1:66,1:166]
O <- read.csv("OR_FX_Size_FoliarCNIsotope_Data.csv")[1:64,1:101]

# The design is close to a basic anova with repeated measures, so we used
# LMER across years with tree as a random effect.

# Only use ones that are being used for the growth analysis (healthy)

#################################################################
########################### Foliar N ############################
#################################################################

##### New York 2016-2019

# Get dataset for foliar N and for treatment and for tree
foliarN_NY_6789 <- 
	c(
	N[!is.na(N$foliar_N_mg_g_02) & N$Use_growth_02 == 1,]$foliar_N_mg_g_02,
	N[!is.na(N$foliar_N_mg_g_04) & N$Use_growth_04 == 1,]$foliar_N_mg_g_04,
	N[!is.na(N$foliar_N_mg_g_07) & N$Use_growth_07 == 1,]$foliar_N_mg_g_07,
	N[!is.na(N$foliar_N_mg_g_08) & N$Use_growth_08 == 1,]$foliar_N_mg_g_08
	)
Treatment_NY_6789 <- 
	c(
	N[!is.na(N$foliar_N_mg_g_02) & N$Use_growth_02 == 1,]$Treatment,
	N[!is.na(N$foliar_N_mg_g_04) & N$Use_growth_04 == 1,]$Treatment,
	N[!is.na(N$foliar_N_mg_g_07) & N$Use_growth_07 == 1,]$Treatment,
	N[!is.na(N$foliar_N_mg_g_08) & N$Use_growth_08 == 1,]$Treatment
	)
Tree_NY_6789 <- 
	as.factor(c(
	N[!is.na(N$foliar_N_mg_g_02) & N$Use_growth_02 == 1,]$PID,
	N[!is.na(N$foliar_N_mg_g_04) & N$Use_growth_04 == 1,]$PID,
	N[!is.na(N$foliar_N_mg_g_07) & N$Use_growth_07 == 1,]$PID,
	N[!is.na(N$foliar_N_mg_g_08) & N$Use_growth_08 == 1,]$PID
	))
Species_NY_6789 <- 
	as.factor(c(
	N[!is.na(N$foliar_N_mg_g_02) & N$Use_growth_02 == 1,]$Species,
	N[!is.na(N$foliar_N_mg_g_04) & N$Use_growth_04 == 1,]$Species,
	N[!is.na(N$foliar_N_mg_g_07) & N$Use_growth_07 == 1,]$Species,
	N[!is.na(N$foliar_N_mg_g_08) & N$Use_growth_08 == 1,]$Species
	))
summary(lme_foliarN_NY_3_02_08 <- lme(foliarN_NY_6789 ~ 0 + Species_NY_6789*Treatment_NY_6789,
	random=~1 | Tree_NY_6789))
print(em_lme_foliarN_NY_3_02_08 <- emmeans(lme_foliarN_NY_3_02_08, list(pairwise ~ Species_NY_6789*Treatment_NY_6789), adjust = "tukey"))
# emmeans 
# All BENI < all ROPS
# Within BENI, L < H, L < P
# Within ROPS, none
# so c c c c, a ab b b
tab <- as.data.frame(summary(em_lme_foliarN_NY_3_02_08)$emmean)
print(round(foliarN_ROPS_3_02_08_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4]),1))
print(round(foliarN_ROPS_3_02_08_MN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4]),1))
print(round(foliarN_ROPS_3_02_08_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4]),1))
print(round(foliarN_ROPS_3_02_08_PHN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4]),1))
print(round(foliarN_BENI_3_02_08_LN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4]),1))
print(round(foliarN_BENI_3_02_08_MN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4]),1))
print(round(foliarN_BENI_3_02_08_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4]),1))
print(round(foliarN_BENI_3_02_08_PHN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4]),1))

##### Oregon 2018-2020

# Get dataset for foliar N and for treatment and for tree
foliarN_OR_890 <- 
	c(
	O[!is.na(O$foliar_N_mg_g_04) & O$Use_growth_04 == 1,]$foliar_N_mg_g_04,
	O[!is.na(O$foliar_N_mg_g_05) & O$Use_growth_05 == 1,]$foliar_N_mg_g_05,
	O[!is.na(O$foliar_N_mg_g_06) & O$Use_growth_06 == 1,]$foliar_N_mg_g_06
	)
Treatment_OR_890 <- 
	c(
	O[!is.na(O$foliar_N_mg_g_04) & O$Use_growth_04 == 1,]$Treatment,
	O[!is.na(O$foliar_N_mg_g_05) & O$Use_growth_05 == 1,]$Treatment,
	O[!is.na(O$foliar_N_mg_g_06) & O$Use_growth_06 == 1,]$Treatment
	)
Tree_OR_890 <- 
	as.factor(c(
	O[!is.na(O$foliar_N_mg_g_04) & O$Use_growth_04 == 1,]$PID,
	O[!is.na(O$foliar_N_mg_g_05) & O$Use_growth_05 == 1,]$PID,
	O[!is.na(O$foliar_N_mg_g_06) & O$Use_growth_06 == 1,]$PID
	))
Species_OR_890 <- 
	as.factor(c(
	O[!is.na(O$foliar_N_mg_g_04) & O$Use_growth_04 == 1,]$Species,
	O[!is.na(O$foliar_N_mg_g_05) & O$Use_growth_05 == 1,]$Species,
	O[!is.na(O$foliar_N_mg_g_06) & O$Use_growth_06 == 1,]$Species
	))
summary(lme_foliarN_OR_3_04_05_06 <- lme(foliarN_OR_890 ~ 0 + Species_OR_890*Treatment_OR_890,
	random=~1 | Tree_OR_890))
print(em_lme_foliarN_OR_3_04_05_06 <- emmeans(lme_foliarN_OR_3_04_05_06, list(pairwise ~ Species_OR_890*Treatment_OR_890), adjust = "tukey"))
# emmeans 
# All PSME < all ALRU
# Within PSME, L < P
# Within ALRU, none
# so c c c c, a ab ab b
tab <- as.data.frame(summary(em_lme_foliarN_OR_3_04_05_06)$emmean)
print(round(foliarN_ALRU_3_04_05_06_LN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4]),1))
print(round(foliarN_ALRU_3_04_05_06_MN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4]),1))
print(round(foliarN_ALRU_3_04_05_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4]),1))
print(round(foliarN_ALRU_3_04_05_06_PHN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4]),1))
print(round(foliarN_PSME_3_04_05_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4]),1))
print(round(foliarN_PSME_3_04_05_06_MN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4]),1))
print(round(foliarN_PSME_3_04_05_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4]),1))
print(round(foliarN_PSME_3_04_05_06_PHN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4]),1))

##### Waiakea 2018-2019

# Get dataset for foliar N and for treatment and for tree
foliarN_HI_W_89 <- 
	c(
	W[!is.na(W$foliar_N_mg_g_04) & W$Use_growth_04 == 1,]$foliar_N_mg_g_04,
	W[!is.na(W$foliar_N_mg_g_06) & W$Use_growth_06 == 1,]$foliar_N_mg_g_06
	)
Treatment_HI_W_89 <- 
	c(
	W[!is.na(W$foliar_N_mg_g_04) & W$Use_growth_04 == 1,]$Treatment,
	W[!is.na(W$foliar_N_mg_g_06) & W$Use_growth_06 == 1,]$Treatment
	)
Tree_HI_W_89 <- 
	as.factor(c(
	W[!is.na(W$foliar_N_mg_g_04) & W$Use_growth_04 == 1,]$PID,
	W[!is.na(W$foliar_N_mg_g_06) & W$Use_growth_06 == 1,]$PID
	))
Species_HI_W_89 <- 
	as.factor(c(
	W[!is.na(W$foliar_N_mg_g_04) & W$Use_growth_04 == 1,]$Species,
	W[!is.na(W$foliar_N_mg_g_06) & W$Use_growth_06 == 1,]$Species
	))
summary(lme_foliarN_HI_W_3_04_06 <- lme(foliarN_HI_W_89 ~ 0 + Species_HI_W_89*Treatment_HI_W_89,
	random=~1 | Tree_HI_W_89))
print(em_lme_foliarN_HI_W_3_04_06 <- emmeans(lme_foliarN_HI_W_3_04_06, list(pairwise ~ Species_HI_W_89*Treatment_HI_W_89), adjust = "tukey"))
# emmeans 
# All PSCA, CAEQ < all GLSE
# Within GLSE, L < P, M < P
# Within CAEQ and PSCA, PL<PP, PL<CL, PL<CM, PL<CH, PL<CP, PM<CP, PH<CP 
# so d d de e, bc bc bc c, a ab ab bc
tab <- as.data.frame(summary(em_lme_foliarN_HI_W_3_04_06)$emmean)
print(round(foliarN_GLSE_3_04_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4]),1))
print(round(foliarN_GLSE_3_04_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4]),1))
print(round(foliarN_GLSE_3_04_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4]),1))
print(round(foliarN_GLSE_3_04_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4]),1))
print(round(foliarN_CAEQ_3_04_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4]),1))
print(round(foliarN_CAEQ_3_04_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4]),1))
print(round(foliarN_CAEQ_3_04_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4]),1))
print(round(foliarN_CAEQ_3_04_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4]),1))
print(round(foliarN_PSCA_3_04_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4]),1))
print(round(foliarN_PSCA_3_04_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4]),1))
print(round(foliarN_PSCA_3_04_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4]),1))
print(round(foliarN_PSCA_3_04_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4]),1))

##### Volcano 2018-2019

# Get dataset for foliar N and for treatment and for tree
foliarN_HI_V_89 <- 
	c(
	V[!is.na(V$foliar_N_mg_g_04) & V$Use_growth_04 == 1,]$foliar_N_mg_g_04,
	V[!is.na(V$foliar_N_mg_g_06) & V$Use_growth_06 == 1,]$foliar_N_mg_g_06
	)
Treatment_HI_V_89 <- 
	c(
	V[!is.na(V$foliar_N_mg_g_04) & V$Use_growth_04 == 1,]$Treatment,
	V[!is.na(V$foliar_N_mg_g_06) & V$Use_growth_06 == 1,]$Treatment
	)
Tree_HI_V_89 <- 
	as.factor(c(
	V[!is.na(V$foliar_N_mg_g_04) & V$Use_growth_04 == 1,]$PID,
	V[!is.na(V$foliar_N_mg_g_06) & V$Use_growth_06 == 1,]$PID
	))
Species_HI_V_89 <- 
	as.factor(c(
	V[!is.na(V$foliar_N_mg_g_04) & V$Use_growth_04 == 1,]$Species,
	V[!is.na(V$foliar_N_mg_g_06) & V$Use_growth_06 == 1,]$Species
	))
summary(lme_foliarN_HI_V_3_04_06 <- lme(foliarN_HI_V_89 ~ 0 + Species_HI_V_89*Treatment_HI_V_89,
	random=~1 | Tree_HI_V_89))
print(em_lme_foliarN_HI_V_3_04_06 <- emmeans(lme_foliarN_HI_V_3_04_06, list(pairwise ~ Species_HI_V_89*Treatment_HI_V_89), adjust = "tukey"))
# emmeans 
# none different
tab <- as.data.frame(summary(em_lme_foliarN_HI_V_3_04_06)$emmean)
print(round(foliarN_ACKO_3_04_06_LN <- c(tab[4,3],tab[4,3]-tab[4,4],tab[4,3]+tab[4,4]),1))
print(round(foliarN_ACKO_3_04_06_MN <- c(tab[7,3],tab[7,3]-tab[7,4],tab[7,3]+tab[7,4]),1))
print(round(foliarN_ACKO_3_04_06_HN <- c(tab[1,3],tab[1,3]-tab[1,4],tab[1,3]+tab[1,4]),1))
print(round(foliarN_ACKO_3_04_06_PHN <- c(tab[10,3],tab[10,3]-tab[10,4],tab[10,3]+tab[10,4]),1))
print(round(foliarN_MOFA_3_04_06_LN <- c(tab[6,3],tab[6,3]-tab[6,4],tab[6,3]+tab[6,4]),1))
print(round(foliarN_MOFA_3_04_06_MN <- c(tab[9,3],tab[9,3]-tab[9,4],tab[9,3]+tab[9,4]),1))
print(round(foliarN_MOFA_3_04_06_HN <- c(tab[3,3],tab[3,3]-tab[3,4],tab[3,3]+tab[3,4]),1))
print(round(foliarN_MOFA_3_04_06_PHN <- c(tab[12,3],tab[12,3]-tab[12,4],tab[12,3]+tab[12,4]),1))
print(round(foliarN_DOVI_3_04_06_LN <- c(tab[5,3],tab[5,3]-tab[5,4],tab[5,3]+tab[5,4]),1))
print(round(foliarN_DOVI_3_04_06_MN <- c(tab[8,3],tab[8,3]-tab[8,4],tab[8,3]+tab[8,4]),1))
print(round(foliarN_DOVI_3_04_06_HN <- c(tab[2,3],tab[2,3]-tab[2,4],tab[2,3]+tab[2,4]),1))
print(round(foliarN_DOVI_3_04_06_PHN <- c(tab[11,3],tab[11,3]-tab[11,4],tab[11,3]+tab[11,4]),1))

##########################################################################
##########################################################################
##########################################################################