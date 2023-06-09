# NSF FX N fixation analysis for the Oregon 15N data: total N fixation

# This script calculates N fixation estimates (total Ndfa; %Ndfa is in another
# file) for the trees at the Oregon site from the Menge, Funk, Perakis, Wolf NSF
# grant from 2015-2020. In other files we do a robustness check for %Ndfa. 
# Here we use the base case for %Ndfa.

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Run the %Ndfa analysis code to load the dataframe and calculate %Ndfa

source("OR_SNF_analysis_pNdfa_forpaper_base.R")

DFstemN_mg_g <- 4.0175 # mg N/g main stem for unfertilized Pseudotsuga
DFrootN_mg_g <- 5.4975 # mg N/g root for unfertilized Pseudotsuga

ARstemrootN_mg_g <- mean(c(5.1,7.4,7.0,3.6,4.2,3.9)) # mg N/g
# These are the six values from 3 twigs and 3 stem cookies from Alnus
# rubra that Steve measured right outside our plots. We don't have 
# root %N numbers, so we use the same values for roots and stems.

####################################################################
########################## Total N fixed ###########################
####################################################################

# Total Ndfa, not just %Ndfa
# To get total N fixed, need to multiply %Ndfa * total N
# Total N is sum of tissue%N*tissuebiomass for all tissues

###
###
### Scenario: tissue %N does not respond to fertilizer (it doesn't in Alnus foliage)
###
###

# Calculate total N in each tree at different time steps

# For leaves, use the measured N content
dat$Leaves_g_N_04 <- dat$Leaves_est_kg_04*dat$foliar_N_mg_g_04
dat$Leaves_g_N_05 <- dat$Leaves_est_kg_05*dat$foliar_N_mg_g_05
dat$Leaves_g_N_06 <- dat$Leaves_est_kg_06*dat$foliar_N_mg_g_06

dat$MainStem_g_N_04 <- dat$MainStem_est_kg_04*ARstemrootN_mg_g
dat[dat$Species=="PSME",]$MainStem_g_N_04 <-
	dat[dat$Species=="PSME",]$MainStem_est_kg_04*DFstemN_mg_g
dat$MainStem_g_N_05 <- dat$MainStem_est_kg_05*ARstemrootN_mg_g
dat[dat$Species=="PSME",]$MainStem_g_N_05 <-
	dat[dat$Species=="PSME",]$MainStem_est_kg_05*DFstemN_mg_g
dat$MainStem_g_N_06 <- dat$MainStem_est_kg_06*ARstemrootN_mg_g
dat[dat$Species=="PSME",]$MainStem_g_N_06 <-
	dat[dat$Species=="PSME",]$MainStem_est_kg_06*DFstemN_mg_g

dat$BGB_g_N_04 <- dat$BGB_est_kg_04*ARstemrootN_mg_g
dat[dat$Species=="PSME",]$BGB_g_N_04 <- 
	dat[dat$Species=="PSME",]$BGB_est_kg_04*DFrootN_mg_g
dat$BGB_g_N_05 <- dat$BGB_est_kg_05*ARstemrootN_mg_g
dat[dat$Species=="PSME",]$BGB_g_N_05 <- 
	dat[dat$Species=="PSME",]$BGB_est_kg_05*DFrootN_mg_g
dat$BGB_g_N_06 <- dat$BGB_est_kg_06*ARstemrootN_mg_g
dat[dat$Species=="PSME",]$BGB_g_N_06 <- 
	dat[dat$Species=="PSME",]$BGB_est_kg_06*DFrootN_mg_g

### Add up biomass N components

dat$Biomass_g_N_04 <- dat$Leaves_g_N_04 + 
	dat$MainStem_g_N_04 + dat$BGB_g_N_04
dat$Biomass_g_N_05 <- dat$Leaves_g_N_05 + 
	dat$MainStem_g_N_05 + dat$BGB_g_N_05
dat$Biomass_g_N_06 <- dat$Leaves_g_N_06 + 
	dat$MainStem_g_N_06 + dat$BGB_g_N_06

### Calculate N fixed

dat$Nfix_g_N_u_04 <- dat$Biomass_g_N_04*dat$Ndfa_u_04/100
dat$Nfix_g_N_u_05 <- dat$Biomass_g_N_05*dat$Ndfa_u_05/100
dat$Nfix_g_N_u_06 <- dat$Biomass_g_N_06*dat$Ndfa_u_06/100

### Calculate N fixed per year

# Average per year (divide by the number of years)
dat$Nfix_g_N_yr_u_04 <- dat$Nfix_g_N_u_04/3 # 3 years: 2016-2018
dat$Nfix_g_N_yr_u_05 <- dat$Nfix_g_N_u_05/4 # 4 years: 2016-2019
dat$Nfix_g_N_yr_u_06 <- dat$Nfix_g_N_u_06/5 # 5 years: 2016-2020

# This year
dat$Nfix_g_N_yr_u_04_05 <- (dat$Biomass_g_N_05 - dat$Biomass_g_N_04)*dat$Ndfa_u_05/100 # 2018 to 2019
dat$Nfix_g_N_yr_u_05_06 <- (dat$Biomass_g_N_06 - dat$Biomass_g_N_05)*dat$Ndfa_u_06/100 # 2019 to 2020

### Calculate N fixed per year per biomass

# Average per year (divide by the number of years)
dat$Nfix_g_N_kg_biomass_yr_u_04 <- dat$Nfix_g_N_yr_u_04/dat$Biomass_est_kg_04 # 3 years: 2016 and 2018
dat$Nfix_g_N_kg_biomass_yr_u_05 <- dat$Nfix_g_N_yr_u_05/dat$Biomass_est_kg_05 # 4 years: 2016 and 2019
dat$Nfix_g_N_kg_biomass_yr_u_06 <- dat$Nfix_g_N_yr_u_06/dat$Biomass_est_kg_06 # 5 years: 2016 and 2020

# This year
dat$Nfix_g_N_kg_biomass_yr_u_04_05 <- dat$Nfix_g_N_yr_u_04_05/dat$Biomass_est_kg_05 # 4 years: 2016 and 2019
dat$Nfix_g_N_kg_biomass_yr_u_05_06 <- dat$Nfix_g_N_yr_u_05_06/dat$Biomass_est_kg_06 # 4 years: 2016 and 2020

### Calculate N fixed per year per biomass N

# Average per year (divide by the number of years)
dat$Nfix_g_N_g_biomass_N_yr_u_04 <- dat$Nfix_g_N_yr_u_04/dat$Biomass_g_N_04 # 3 years: 2016 and 2018
dat$Nfix_g_N_g_biomass_N_yr_u_05 <- dat$Nfix_g_N_yr_u_05/dat$Biomass_g_N_05 # 4 years: 2016 and 2019
dat$Nfix_g_N_g_biomass_N_yr_u_06 <- dat$Nfix_g_N_yr_u_06/dat$Biomass_g_N_06 # 5 years: 2016 and 2020

# This year
dat$Nfix_g_N_g_biomass_N_yr_u_04_05 <- dat$Nfix_g_N_yr_u_04_05/dat$Biomass_g_N_05 # 4 years: 2016 and 2019
dat$Nfix_g_N_g_biomass_N_yr_u_05_06 <- dat$Nfix_g_N_yr_u_05_06/dat$Biomass_g_N_06 # 5 years: 2016 and 2020

#####
### Stats and print out for table
#####

### LME across 2016-2020 with tree as random, no error propagation
# 2018-2020
# Get dataset for tNdfa_kg_y and for treatment and for tree
tNdfa_kg_y_u_890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Nfix_g_N_kg_biomass_yr_u_04,
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_05) & 
		dat$Use_growth_05 == 1,]$Nfix_g_N_kg_biomass_yr_u_05,
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Nfix_g_N_kg_biomass_yr_u_06
	)
Treatment_ALRU_890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_05) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Tree_ALRU_890 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_05) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_kg_biomass_yr_u_06) & 
		dat$Use_growth_06 == 1,]$PID
	))

summary(lme_tNdfa_kg_y_base <- lme(tNdfa_kg_y_u_890 ~ 0 + Treatment_ALRU_890,
	random=~1 | Tree_ALRU_890))
print(em_lme_tNdfa_kg_y_base <- emmeans(lme_tNdfa_kg_y_base, list(pairwise ~ Treatment_ALRU_890), adjust = "tukey"))
# emmeans a b b b, all > 0 (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_tNdfa_kg_y_base)$emmean)
print(round(tNdfa_kg_y_base_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),2))
print(round(tNdfa_kg_y_base_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),2))
print(round(tNdfa_kg_y_base_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),2))
print(round(tNdfa_kg_y_base_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),2))

# Get dataset for tNdfa_y and for treatment and for tree
tNdfa_y_u_890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_yr_u_04) & 
		dat$Use_growth_04 == 1,]$Nfix_g_N_yr_u_04,
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_yr_u_05) & 
		dat$Use_growth_05 == 1,]$Nfix_g_N_yr_u_05,
	dat[dat$Species=="ALRU" & !is.na(dat$Nfix_g_N_yr_u_06) & 
		dat$Use_growth_06 == 1,]$Nfix_g_N_yr_u_06
	)
summary(lme_tNdfa_y_base <- lme(tNdfa_y_u_890 ~ 0 + Treatment_ALRU_890,
	random=~1 | Tree_ALRU_890))
print(em_lme_tNdfa_y_base <- emmeans(lme_tNdfa_y_base, list(pairwise ~ Treatment_ALRU_890), adjust = "tukey"))
# emmeans a a a a, all > 0 (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_tNdfa_y_base)$emmean)
print(round(tNdfa_y_base_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(tNdfa_y_base_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(tNdfa_y_base_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(tNdfa_y_base_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

# Uncomment this to write out a spreadsheet with the N fixation data
#write.csv(dat[,c(2,5,7,8,44,54,55,61,71,72,77,87,88,96:98,108:110,168:170,
#	180:200)],"OR_FX_SNF.csv")

##########################################################################
##########################################################################
##########################################################################