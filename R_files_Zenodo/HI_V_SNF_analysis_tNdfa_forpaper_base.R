# NSF FX N fixation analysis for the Volcano 15N data: total N fixation

# This script calculates N fixation estimates (total Ndfa; %Ndfa is in another
# file) for the trees at Volcano from the Menge, Funk, Perakis, Wolf NSF
# grant from 2015-2020. In other files we do a robustness check for %Ndfa. 
# Here we use the base case for %Ndfa.

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Run the %Ndfa analysis code to load the dataframe and calculate %Ndfa

source("HI_V_SNF_analysis_pNdfa_forpaper_base.R")

bmcols <- read.csv("HI_V_FX_AllBiomassComponents.csv")
# This file contains the biomass component (leaf, twig, secondary stem, main 
# stem) data, which were estimated by running the size measurements for each
# tree through the allometric equations developed for these trees (see
# Carreras Pereira et al.).

tissuechem <- read.csv("HI_V_FX_TissueCN_FilledData.csv")
# This file contains the twig, secondary stem, main stem, coarse root, and 
# fine root tissue C and N data. Leaf data are in the dat dataframe that
# was read in for the pNdfa file. For the tissuechem data, some of the 
# trees either didn't have a tissue type at harvest or were not harvested,
# so they were estimated from treatment means for each species as described
# in the methods text.

####################################################################
########################## Total N fixed ###########################
####################################################################

# Total Ndfa, not just %Ndfa
# To get total N fixed, need to multiply %Ndfa * total N
# Total N is sum of tissue%N*tissuebiomass for all tissues

# Calculate total N in each tree at different time steps

# Leaves
dat$Leaves_g_N_04 <- bmcols$Leaves_est_kg_04*dat$foliar_N_mg_g_04
dat$Leaves_g_N_06 <- bmcols$Leaves_est_kg_06*dat$foliar_N_mg_g_06

# Twigs (this assumes the same twig % N at each time step)
dat$Twigs_g_N_04 <- bmcols$Twigs_est_kg_04*tissuechem$twig_N_mg_g_06
dat$Twigs_g_N_06 <- bmcols$Twigs_est_kg_06*tissuechem$twig_N_mg_g_06

# SecondaryStem (this assumes the same secondary stem % N at each time step)
dat$SecondaryStem_g_N_04 <- bmcols$SecondaryStem_est_kg_04*tissuechem$secondarystem_N_mg_g_06
dat$SecondaryStem_g_N_06 <- bmcols$SecondaryStem_est_kg_06*tissuechem$secondarystem_N_mg_g_06
dat[is.na(dat$SecondaryStem_g_N_04),]$SecondaryStem_g_N_04 <- 0
dat[is.na(dat$SecondaryStem_g_N_06),]$SecondaryStem_g_N_06 <- 0

# MainStem (this assumes the same main stem % N at each time step)
dat$MainStem_g_N_04 <- bmcols$MainStem_est_kg_04*tissuechem$mainstem_N_mg_g_06
dat$MainStem_g_N_06 <- bmcols$MainStem_est_kg_06*tissuechem$mainstem_N_mg_g_06

### Add up biomass N components

dat$AGB_g_N_04 <- dat$Leaves_g_N_04 + dat$Twigs_g_N_04 + 
	dat$MainStem_g_N_04 + dat$SecondaryStem_g_N_04
dat$AGB_g_N_06 <- dat$Leaves_g_N_06 + dat$Twigs_g_N_06 + 
	dat$MainStem_g_N_06 + dat$SecondaryStem_g_N_06

### Calculate N fixed

dat$Nfix_g_N_u_04 <- dat$AGB_g_N_04*dat$Ndfa_u_04/100
dat$Nfix_g_N_u_06 <- dat$AGB_g_N_06*dat$Ndfa_u_06/100

### Calculate N fixed per year

# Average per year (divide by the number of years)
dat$Nfix_g_N_yr_u_04 <- dat$Nfix_g_N_u_04/3 # 3 years: 2016-2018
dat$Nfix_g_N_yr_u_06 <- dat$Nfix_g_N_u_06/4 # 4 years: 2016-2019

# This year
dat$Nfix_g_N_yr_u_04_06 <- (dat$AGB_g_N_06 - dat$AGB_g_N_04)*dat$Ndfa_u_06/100 # 2018 to 2019

### Calculate N fixed per year per biomass

# Average per year (divide by the number of years)
dat$Nfix_g_N_kg_biomass_yr_u_04 <- dat$Nfix_g_N_yr_u_04/dat$AGB_est_kg_04 # 3 years: 2016 and 2018
dat$Nfix_g_N_kg_biomass_yr_u_06 <- dat$Nfix_g_N_yr_u_06/dat$AGB_est_kg_06 # 4 years: 2016 and 2019

# This year
dat$Nfix_g_N_kg_biomass_yr_u_04_06 <- dat$Nfix_g_N_yr_u_04_06/dat$AGB_est_kg_06 # 4 years: 2016 and 2019

### Calculate N fixed per year per biomass N

# Average per year (divide by the number of years)
dat$Nfix_g_N_g_biomass_N_yr_u_04 <- dat$Nfix_g_N_yr_u_04/dat$AGB_g_N_04 # 3 years: 2016 and 2018
dat$Nfix_g_N_g_biomass_N_yr_u_06 <- dat$Nfix_g_N_yr_u_06/dat$AGB_g_N_06 # 4 years: 2016 and 2019

# This year
dat$Nfix_g_N_g_biomass_N_yr_u_04_06 <- dat$Nfix_g_N_yr_u_04_06/dat$AGB_g_N_06 # 4 years: 2016 and 2019

#####
### Stats and print out for table
#####

# ACKO

# 2018-2019
# Get dataset for tNdfa_kg_y and for treatment and for tree
tNdfa_kg_y_u_ACKO_89 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Nfix_g_N_kg_biomass_yr_u_04,
	dat[dat$Species=="ACKO" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Nfix_g_N_kg_biomass_yr_u_06
	)
Treatment_ACKO_89 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ACKO" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Tree_ACKO_89 <- 
	as.factor(c(
	dat[dat$Species=="ACKO" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ACKO" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$PID
	))

summary(lme_tNdfa_kg_y_base_ACKO <- lme(tNdfa_kg_y_u_ACKO_89 ~ 0 + Treatment_ACKO_89,
	random=~1 | Tree_ACKO_89))
print(em_lme_tNdfa_kg_y_base_ACKO <- emmeans(lme_tNdfa_kg_y_base_ACKO, list(pairwise ~ Treatment_ACKO_89), adjust = "tukey"))
# emmeans a a a a, all > 0 (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_tNdfa_kg_y_base_ACKO)$emmean)
print(round(tNdfa_kg_y_base_ACKO_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),2))
print(round(tNdfa_kg_y_base_ACKO_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),2))
print(round(tNdfa_kg_y_base_ACKO_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),2))
print(round(tNdfa_kg_y_base_ACKO_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),2))

# Get dataset for tNdfa_y and for treatment and for tree
tNdfa_y_u_ACKO_89 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$Nfix_g_N_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Nfix_g_N_yr_u_04,
	dat[dat$Species=="ACKO" & !is.na(dat$Nfix_g_N_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Nfix_g_N_yr_u_06
	)
summary(lme_tNdfa_y_base_ACKO <- lme(tNdfa_y_u_ACKO_89 ~ 0 + Treatment_ACKO_89,
	random=~1 | Tree_ACKO_89))
print(em_lme_tNdfa_y_base_ACKO <- emmeans(lme_tNdfa_y_base_ACKO, list(pairwise ~ Treatment_ACKO_89), adjust = "tukey"))
# emmeans a a a a, PHN > 0
tab <- as.data.frame(summary(em_lme_tNdfa_y_base_ACKO)$emmean)
print(round(tNdfa_y_base_ACKO_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(tNdfa_y_base_ACKO_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(tNdfa_y_base_ACKO_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(tNdfa_y_base_ACKO_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

# MOFA

# 2018-2019
# Get dataset for tNdfa_kg_y and for treatment and for tree
tNdfa_kg_y_u_MOFA_89 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Nfix_g_N_kg_biomass_yr_u_04,
	dat[dat$Species=="MOFA" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Nfix_g_N_kg_biomass_yr_u_06
	)
Treatment_MOFA_89 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="MOFA" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Tree_MOFA_89 <- 
	as.factor(c(
	dat[dat$Species=="MOFA" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="MOFA" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$PID
	))

summary(lme_tNdfa_kg_y_base_MOFA <- lme(tNdfa_kg_y_u_MOFA_89 ~ 0 + Treatment_MOFA_89,
	random=~1 | Tree_MOFA_89))
print(em_lme_tNdfa_kg_y_base_MOFA <- emmeans(lme_tNdfa_kg_y_base_MOFA, list(pairwise ~ Treatment_MOFA_89), adjust = "tukey"))
# emmeans a a a a, all > 0 except HN (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_tNdfa_kg_y_base_MOFA)$emmean)
print(round(tNdfa_kg_y_base_MOFA_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),2))
print(round(tNdfa_kg_y_base_MOFA_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),2))
print(round(tNdfa_kg_y_base_MOFA_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),2))
print(round(tNdfa_kg_y_base_MOFA_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),2))

# Get dataset for tNdfa_kg_y and for treatment and for tree
tNdfa_y_u_MOFA_89 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$Nfix_g_N_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Nfix_g_N_yr_u_04,
	dat[dat$Species=="MOFA" & !is.na(dat$Nfix_g_N_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Nfix_g_N_yr_u_06
	)
summary(lme_tNdfa_y_base_MOFA <- lme(tNdfa_y_u_MOFA_89 ~ 0 + Treatment_MOFA_89,
	random=~1 | Tree_MOFA_89))
print(em_lme_tNdfa_y_base_MOFA <- emmeans(lme_tNdfa_y_base_MOFA, list(pairwise ~ Treatment_MOFA_89), adjust = "tukey"))
# emmeans a a a a, HN > 0 (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_tNdfa_y_base_MOFA)$emmean)
print(round(tNdfa_y_base_MOFA_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(tNdfa_y_base_MOFA_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(tNdfa_y_base_MOFA_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(tNdfa_y_base_MOFA_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

# Uncomment this to write out a spreadsheet with the N fixation data
#write.csv(dat[,c(2,5,8,9,49,59,60,77,87,88,101,109,111,124,125,134:146)],
#	"HI_V_FX_SNF.csv")

##########################################################################
##########################################################################
##########################################################################