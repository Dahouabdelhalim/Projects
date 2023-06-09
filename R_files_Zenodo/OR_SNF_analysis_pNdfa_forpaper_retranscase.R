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
	# - Base case (different file): 
		# Retranslocation is similar for Pseudotsuga and Alnus
	# - Retranslocation sensitivity check (this file):
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

########################################################################
### Retranslocation Case: Alnus uses 10% less old N than Pseudotsuga ###
########################################################################

# Paired reference trees where possible, otherwise means of the treatment
# Retranslocation 10% less for Alnus than for Pseudotsuga
# 2018, 2019, 2020

# Reset columns for all years

dat$ref_15N_AP_u_04 <- NA
dat[dat$Species=="ALRU",]$ref_15N_AP_u_04 <- dat[dat$Species=="PSME",]$foliar_15N_AP_04

dat$ref_15N_AP_u_05 <- NA
dat[dat$Species=="ALRU",]$ref_15N_AP_u_05 <- dat[dat$Species=="PSME",]$foliar_15N_AP_05

dat$ref_15N_AP_u_06 <- NA
dat[dat$Species=="ALRU",]$ref_15N_AP_u_06 <- dat[dat$Species=="PSME",]$foliar_15N_AP_06

#####
### Calculate %Ndfa
#####

pnon <- 0.2 # Assume that the non-fixer uses 20% retranslocated N each new year
pfix <- 0.1 # Assume that the fixer uses 10% retranslocated N each new year

# Use this new equation for %Ndfa:
# %Ndfa,t = ((%15Nref,t ? pnon*%15Nref,t-1)/(1-pnon) ? 
#	(%15Nsample,t ? pfix*%15Nsample,t-1)/(1-pfix))/
#	((%15Nref,t ? pnon*%15Nref,t-1)/(1-pnon) ? %15Nfixation)

# %15Nref,t is ref_15N_AP_u_THISYEAR
# %15Nsample,t is foliar_15N_AP_THISYEAR
# %15Nref,t-1 is ref_15N_AP_u_LASTYEAR
# %15Nsample,t-1 is foliar_15N_AP_LASTYEAR
	# For _04, this is purefix_15N_AP
	# For _05, this is _04 
	# For _06, this is _05
# %15Nfixation is purefix_15N_AP

# 2018

dat$Ndfa_u_04 <- NA
dat[dat$Species=="ALRU",]$Ndfa_u_04 <- 
	100*((dat[dat$Species=="ALRU",]$ref_15N_AP_u_04 - pnon*purefix_15N_AP)/(1-pnon) -
	(dat[dat$Species=="ALRU",]$foliar_15N_AP_04 - pfix*purefix_15N_AP)/(1-pfix))/
	((dat[dat$Species=="ALRU",]$ref_15N_AP_u_04 - pnon*purefix_15N_AP)/(1-pnon) -
	purefix_15N_AP)

# 2019

dat$Ndfa_u_05 <- NA
dat[dat$Species=="ALRU",]$Ndfa_u_05 <- 
	100*((dat[dat$Species=="ALRU",]$ref_15N_AP_u_05 - 
	pnon*dat[dat$Species=="ALRU",]$ref_15N_AP_u_04)/(1-pnon) -
	(dat[dat$Species=="ALRU",]$foliar_15N_AP_05 - 
	pfix*dat[dat$Species=="ALRU",]$foliar_15N_AP_04)/(1-pfix))/
	((dat[dat$Species=="ALRU",]$ref_15N_AP_u_05 - 
	pnon*dat[dat$Species=="ALRU",]$ref_15N_AP_u_04)/(1-pnon) -
	purefix_15N_AP)

# 2020

dat$Ndfa_u_06 <- NA
dat[dat$Species=="ALRU",]$Ndfa_u_06 <- 
	100*((dat[dat$Species=="ALRU",]$ref_15N_AP_u_06 - 
	pnon*dat[dat$Species=="ALRU",]$ref_15N_AP_u_05)/(1-pnon) -
	(dat[dat$Species=="ALRU",]$foliar_15N_AP_06 - 
	pfix*dat[dat$Species=="ALRU",]$foliar_15N_AP_05)/(1-pfix))/
	((dat[dat$Species=="ALRU",]$ref_15N_AP_u_06 - 
	pnon*dat[dat$Species=="ALRU",]$ref_15N_AP_u_05)/(1-pnon) -
	purefix_15N_AP)

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
summary(lme_pNdfa_retranscase <- lme(pNdfa_u_890 ~ 0 + Treatment_ALRU_890,
	random=~1 | Tree_ALRU_890))
print(em_lme_pNdfa_retranscase <- emmeans(lme_pNdfa_retranscase, list(pairwise ~ Treatment_ALRU_890), adjust = "tukey"))
# emmeans a b b b, all > 0 (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_pNdfa_retranscase)$emmean)
print(round(pNdfa_retranscase_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(pNdfa_retranscase_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(pNdfa_retranscase_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(pNdfa_retranscase_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

####################################################################
####################################################################
####################################################################