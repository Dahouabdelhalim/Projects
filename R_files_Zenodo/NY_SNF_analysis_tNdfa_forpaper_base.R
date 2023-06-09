# NSF FX N fixation analysis for the New York 15N data: total N fixation

# This script calculates N fixation estimates (total Ndfa; %Ndfa is in another
# file) for the trees at Black Rock from the Menge, Funk, Perakis, Wolf NSF
# grant from 2015-2020. In other files we do a robustness check for %Ndfa. 
# Here we use the base case for %Ndfa.

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Run the %Ndfa analysis code to load the dataframe and calculate %Ndfa

source("NY_SNF_analysis_pNdfa_forpaper_base.R")

bmcols <- read.csv("NY_FX_AllBiomassComponents.csv")
# This file contains the biomass component (leaf, twig, secondary stem, main 
# stem) data, which were estimated by running the size measurements for each
# tree through the allometric equations developed for these trees (see
# Carreras Pereira et al.).

tissuechem <- read.csv("NY_FX_TissueCN_FilledData.csv")
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
dat$Leaves_g_N_02 <- bmcols$Leaves_est_kg_02*dat$foliar_N_mg_g_02
dat$Leaves_g_N_04 <- bmcols$Leaves_est_kg_04*dat$foliar_N_mg_g_04
dat$Leaves_g_N_07 <- bmcols$Leaves_est_kg_07*dat$foliar_N_mg_g_07
dat$Leaves_g_N_08 <- bmcols$Leaves_est_kg_08*dat$foliar_N_mg_g_08

# Twigs (this assumes the same twig % N at each time step)
dat$Twigs_g_N_02 <- bmcols$Twigs_est_kg_02*tissuechem$twig_N_mg_g_08
dat$Twigs_g_N_04 <- bmcols$Twigs_est_kg_04*tissuechem$twig_N_mg_g_08
dat$Twigs_g_N_07 <- bmcols$Twigs_est_kg_07*tissuechem$twig_N_mg_g_08
dat$Twigs_g_N_08 <- bmcols$Twigs_est_kg_08*tissuechem$twig_N_mg_g_08

# SecondaryStem (this assumes the same secondary stem % N at each time step)
dat$SecondaryStem_g_N_02 <- bmcols$SecondaryStem_est_kg_02*tissuechem$secondarystem_N_mg_g_08
dat$SecondaryStem_g_N_04 <- bmcols$SecondaryStem_est_kg_04*tissuechem$secondarystem_N_mg_g_08
dat$SecondaryStem_g_N_07 <- bmcols$SecondaryStem_est_kg_07*tissuechem$secondarystem_N_mg_g_08
dat$SecondaryStem_g_N_08 <- bmcols$SecondaryStem_est_kg_08*tissuechem$secondarystem_N_mg_g_08
dat[is.na(dat$SecondaryStem_g_N_02),]$SecondaryStem_g_N_02 <- 0
dat[is.na(dat$SecondaryStem_g_N_04),]$SecondaryStem_g_N_04 <- 0
dat[is.na(dat$SecondaryStem_g_N_07),]$SecondaryStem_g_N_07 <- 0
dat[is.na(dat$SecondaryStem_g_N_08),]$SecondaryStem_g_N_08 <- 0

# MainStem (this assumes the same main stem % N at each time step)
dat$MainStem_g_N_02 <- bmcols$MainStem_est_kg_02*tissuechem$mainstem_N_mg_g_08
dat$MainStem_g_N_04 <- bmcols$MainStem_est_kg_04*tissuechem$mainstem_N_mg_g_08
dat$MainStem_g_N_07 <- bmcols$MainStem_est_kg_07*tissuechem$mainstem_N_mg_g_08
dat$MainStem_g_N_08 <- bmcols$MainStem_est_kg_08*tissuechem$mainstem_N_mg_g_08

# Roots (this assumes the same root % N at each time step)
dat$BGB_g_N_02 <- dat$BGB_est_kg_02*tissuechem$coarseroot_N_mg_g_08
dat$BGB_g_N_04 <- dat$BGB_est_kg_04*tissuechem$coarseroot_N_mg_g_08
dat$BGB_g_N_07 <- dat$BGB_est_kg_07*tissuechem$coarseroot_N_mg_g_08
dat$BGB_g_N_08 <- dat$BGB_est_kg_08*tissuechem$coarseroot_N_mg_g_08

### Add up biomass N components

dat$Biomass_g_N_02 <- dat$Leaves_g_N_02 + dat$Twigs_g_N_02 + 
	dat$MainStem_g_N_02 + dat$SecondaryStem_g_N_02 + dat$BGB_g_N_02
dat$Biomass_g_N_04 <- dat$Leaves_g_N_04 + dat$Twigs_g_N_04 + 
	dat$MainStem_g_N_04 + dat$SecondaryStem_g_N_04 + dat$BGB_g_N_04
dat$Biomass_g_N_07 <- dat$Leaves_g_N_07 + dat$Twigs_g_N_07 + 
	dat$MainStem_g_N_07 + dat$SecondaryStem_g_N_07 + dat$BGB_g_N_07
dat$Biomass_g_N_08 <- dat$Leaves_g_N_08 + dat$Twigs_g_N_08 + 
	dat$MainStem_g_N_08 + dat$SecondaryStem_g_N_08 + dat$BGB_g_N_08

### Calculate N fixed

dat$Nfix_g_N_u_04 <- dat$Biomass_g_N_04*dat$Ndfa_u_04/100
dat$Nfix_g_N_u_07 <- dat$Biomass_g_N_07*dat$Ndfa_u_07/100
dat$Nfix_g_N_u_08 <- dat$Biomass_g_N_08*dat$Ndfa_u_08/100

### Calculate N fixed per year

# Average per year (divide by the number of years)
dat$Nfix_g_N_yr_u_04 <- dat$Nfix_g_N_u_04/3 # 3 years: 2015-2017
dat$Nfix_g_N_yr_u_07 <- dat$Nfix_g_N_u_07/4 # 4 years: 2015-2018
dat$Nfix_g_N_yr_u_08 <- dat$Nfix_g_N_u_08/5 # 5 years: 2015-2019

# This year
dat$Nfix_g_N_yr_u_04_07 <- (dat$Biomass_g_N_07 - dat$Biomass_g_N_04)*dat$Ndfa_u_07/100 # 2017 to 2018
dat$Nfix_g_N_yr_u_07_08 <- (dat$Biomass_g_N_08 - dat$Biomass_g_N_07)*dat$Ndfa_u_08/100 # 2018 to 2019

### Calculate N fixed per year per biomass

# Average per year (divide by the number of years)
dat$Nfix_g_N_kg_biomass_yr_u_04 <- dat$Nfix_g_N_yr_u_04/dat$Biomass_est_kg_04 # 3 years: 2015 and 2017
dat$Nfix_g_N_kg_biomass_yr_u_07 <- dat$Nfix_g_N_yr_u_07/dat$Biomass_est_kg_07 # 4 years: 2015 and 2018
dat$Nfix_g_N_kg_biomass_yr_u_08 <- dat$Nfix_g_N_yr_u_08/dat$Biomass_est_kg_08 # 5 years: 2015 and 2019

# This year
dat$Nfix_g_N_kg_biomass_yr_u_04_07 <- dat$Nfix_g_N_yr_u_04_07/dat$Biomass_est_kg_07 # 4 years: 2015 and 2018
dat$Nfix_g_N_kg_biomass_yr_u_07_08 <- dat$Nfix_g_N_yr_u_07_08/dat$Biomass_est_kg_08 # 5 years: 2015 and 2019

### Calculate N fixed per year per biomass N

# Average per year (divide by the number of years)
dat$Nfix_g_N_g_biomass_N_yr_u_04 <- dat$Nfix_g_N_yr_u_04/dat$Biomass_g_N_04 # 3 years: 2015 and 2017
dat$Nfix_g_N_g_biomass_N_yr_u_07 <- dat$Nfix_g_N_yr_u_07/dat$Biomass_g_N_07 # 4 years: 2015 and 2018
dat$Nfix_g_N_g_biomass_N_yr_u_08 <- dat$Nfix_g_N_yr_u_08/dat$Biomass_g_N_08 # 5 years: 2015 and 2019

# This year
dat$Nfix_g_N_g_biomass_N_yr_u_04_07 <- dat$Nfix_g_N_yr_u_04_07/dat$Biomass_g_N_07 # 4 years: 2015 and 2018
dat$Nfix_g_N_g_biomass_N_yr_u_07_08 <- dat$Nfix_g_N_yr_u_07_08/dat$Biomass_g_N_08 # 5 years: 2015 and 2019

#####
### Stats and print out for table
#####

# 2017-2019
# Get dataset for tNdfa_kg_y and for treatment and for tree
tNdfa_kg_y_u_789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Nfix_g_N_kg_biomass_yr_u_04,
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_07) & 
		dat$Use_growth_07 == 1,]$Nfix_g_N_kg_biomass_yr_u_07,
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_08) & 
		dat$Use_growth_08 == 1,]$Nfix_g_N_kg_biomass_yr_u_08
	)
Treatment_ROPS_789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_07) & 
		dat$Use_growth_07 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_08) & 
		dat$Use_growth_08 == 1,]$Treatment
	)
Tree_ROPS_789 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_07) & 
		dat$Use_growth_07 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_08) & 
		dat$Use_growth_08 == 1,]$PID
	))

summary(lme_tNdfa_kg_y_base <- lme(tNdfa_kg_y_u_789 ~ 0 + Treatment_ROPS_789,
	random=~1 | Tree_ROPS_789))
print(em_lme_tNdfa_kg_y_base <- emmeans(lme_tNdfa_kg_y_base, list(pairwise ~ Treatment_ROPS_789), adjust = "tukey"))
# emmeans a b b b, all > 0 except PHN (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_tNdfa_kg_y_base)$emmean)
print(round(tNdfa_kg_y_base_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),2))
print(round(tNdfa_kg_y_base_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),2))
print(round(tNdfa_kg_y_base_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),2))
print(round(tNdfa_kg_y_base_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),2))

# 2017-2019
# Get dataset for tNdfa_y and for treatment and for tree
tNdfa_y_u_789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Nfix_g_N_yr_u_04,
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_yr_u_07) & 
		dat$Use_growth_07 == 1,]$Nfix_g_N_yr_u_07,
	dat[dat$Species=="ROPS" & !is.na(dat$Nfix_g_N_yr_u_08) & 
		dat$Use_growth_08 == 1,]$Nfix_g_N_yr_u_08
	)
summary(lme_tNdfa_y_base <- lme(tNdfa_y_u_789 ~ 0 + Treatment_ROPS_789,
	random=~1 | Tree_ROPS_789))
print(em_lme_tNdfa_y_base <- emmeans(lme_tNdfa_y_base, list(pairwise ~ Treatment_ROPS_789), adjust = "tukey"))
# emmeans a a a a, all > 0 except PHN (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_tNdfa_y_base)$emmean)
print(round(tNdfa_y_base_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(tNdfa_y_base_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(tNdfa_y_base_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(tNdfa_y_base_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

# Uncomment this to write out a spreadsheet with the N fixation data
#write.csv(dat[,c(2,5,7,8,54,64,65,89,99,100,111,123,124,138,145,149,151,
#	152,163,165,166,214:216,238:258)],"NY_FX_SNF.csv")

##########################################################################
##########################################################################
##########################################################################