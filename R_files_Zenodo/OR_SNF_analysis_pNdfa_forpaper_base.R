# NSF FX N fixation analysis for the Oregon 15N data

# This script calculates N fixation estimates (%Ndfa; total Ndfa is in 
# another file) for the trees at the Oregon site from the Menge, Funk, Perakis, 
# Wolf NSF grant from 2015-2020. This focuses on the calculations and 
# statistics used in the paper.

# We want to do the %Ndfa calculation in a number of ways to assess how
# robust the results are to our assumptions. Here are details of these ways:
# - Which reference trees to use?
	# - Base case (this file): 
		# Paired trees when possible, 
		# means of the treatment when not possible
	# - Reference tree sensitivity check (different file): 
		# Means of treatments for all trees
# - Years (all in this file)
	# - 2018
	# - 2019
	# - 2020
# - Assumptions about retranslocation
	# - Base case (this file): 
		# Retranslocation is similar for Pseudotsuga and Alnus
	# - Retranslocation sensitivity check (different file):
		# Pseudotsuga uses more retranslocated N than Alnus
	# - Means of treatments for all trees

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Run the growth analysis code to load the dataframe and calculate growth rates

source("OR_growth_analysis.R")

#################################################################################
####################### N fixation estimates (%Ndfa) ############################
#################################################################################

purefix_15N_AP <- 0.3663

##################
### Base case: ###
##################

# Paired reference trees where possible, otherwise means of the treatment
# Retranslocation similar for Pseudotsuga and Alnus
# 2018, 2019, 2020

# Set up columns for all years and set the reference to the paired tree

dat$ref_15N_AP_u_04 <- NA
dat[dat$Species=="ALRU",]$ref_15N_AP_u_04 <- dat[dat$Species=="PSME",]$foliar_15N_AP_04

dat$ref_15N_AP_u_05 <- NA
dat[dat$Species=="ALRU",]$ref_15N_AP_u_05 <- dat[dat$Species=="PSME",]$foliar_15N_AP_05

dat$ref_15N_AP_u_06 <- NA
dat[dat$Species=="ALRU",]$ref_15N_AP_u_06 <- dat[dat$Species=="PSME",]$foliar_15N_AP_06

# Where there isn't a paired reference tree, set the reference to 
# the mean of all PSMEs paired to healthy ALRUes in the same treatment.
# Keeping this note here, but not applicable in Oregon, as all PSMEs were healthy.

#####
### Calculate %Ndfa
#####

# 2018

dat$Ndfa_u_04 <- NA
dat[dat$Species=="ALRU",]$Ndfa_u_04 <- 
	100*(dat[dat$Species=="ALRU",]$ref_15N_AP_u_04 - 
	dat[dat$Species=="ALRU",]$foliar_15N_AP_04)/
	(dat[dat$Species=="ALRU",]$ref_15N_AP_u_04 - purefix_15N_AP)

# 2019

dat$Ndfa_u_05 <- NA
dat[dat$Species=="ALRU",]$Ndfa_u_05 <- 
	100*(dat[dat$Species=="ALRU",]$ref_15N_AP_u_05 - 
	dat[dat$Species=="ALRU",]$foliar_15N_AP_05)/
	(dat[dat$Species=="ALRU",]$ref_15N_AP_u_05 - purefix_15N_AP)

# 2020

dat$Ndfa_u_06 <- NA
dat[dat$Species=="ALRU",]$Ndfa_u_06 <- 
	100*(dat[dat$Species=="ALRU",]$ref_15N_AP_u_06 - 
	dat[dat$Species=="ALRU",]$foliar_15N_AP_06)/
	(dat[dat$Species=="ALRU",]$ref_15N_AP_u_06 - purefix_15N_AP)

#####
### Stats and print out for table
#####

# 2018-2020
# Get dataset for pNdfa and for treatment and for tree
pNdfa_u_890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$Ndfa_u_04,
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_05) & 
		dat$Use_growth_05 == 1,]$Ndfa_u_05,
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$Ndfa_u_06
	)
Treatment_ALRU_890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_05) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Tree_ALRU_890 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_05) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_pNdfa_base <- lme(pNdfa_u_890 ~ 0 + Treatment_ALRU_890,
	random=~1 | Tree_ALRU_890))
print(em_lme_pNdfa_base <- emmeans(lme_pNdfa_base, list(pairwise ~ Treatment_ALRU_890), adjust = "tukey"))
# emmeans: a b b b, all > 0 (based on 95% CI)
tab <- as.data.frame(summary(em_lme_pNdfa_base)$emmean)
print(round(pNdfa_base_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(pNdfa_base_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(pNdfa_base_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(pNdfa_base_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

####################################################################
####################################################################
####################################################################