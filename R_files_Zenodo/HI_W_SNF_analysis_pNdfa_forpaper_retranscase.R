# NSF FX N fixation analysis for the Waiakea 15N data

# This script calculates N fixation estimates (%Ndfa; total Ndfa is in 
# another file) for the trees at the Waiakea site from the Menge, Funk, Perakis, 
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
# - Years
	# - 2018
	# - 2019
# - Assumptions about retranslocation
	# - Base case (different file): 
		# Retranslocation is similar for Psidium and Gliricidia and Casuarina
	# - Retranslocation sensitivity check (this file):
		# Psidium uses more retranslocated N than Gliricidia and Casuarina

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Run the growth analysis code to load the dataframe and calculate growth rates

source("HI_W_growth_analysis.R")

#################################################################################
####################### N fixation estimates (%Ndfa) ############################
#################################################################################

purefix_15N_AP <- 0.3663

##################
##################
### Gliricidia ###
##################
##################

#########################################################################
### Retranslocation Case: Gliricidia uses 10% less old N than Psidium ###
#########################################################################

# Paired reference trees where possible, otherwise means of the treatment
# Retranslocation 10% less for Gliricidia than for Psidium
# 2018, 2019

# Set unhealthy PSCAs foliar d15N to NA
dat[dat$Species=="PSCA" & dat$Use_growth_04==0,]$foliar_15N_AP_04 <- NA
dat[dat$Species=="PSCA" & dat$Use_growth_06==0,]$foliar_15N_AP_06 <- NA

# Reset columns for all years

dat$ref_15N_AP_u_04 <- NA
dat[dat$Species=="GLSE",]$ref_15N_AP_u_04 <- dat[dat$Species=="PSCA",]$foliar_15N_AP_04

dat$ref_15N_AP_u_06 <- NA
dat[dat$Species=="GLSE",]$ref_15N_AP_u_06 <- dat[dat$Species=="PSCA",]$foliar_15N_AP_06

# Waiakea only: Filter out the ones that weren't labeled
dat[dat$Species=="GLSE",]$ref_15N_AP_u_04[c(9,18,27)] <- NA
dat[dat$Species=="GLSE",]$ref_15N_AP_u_06[c(9,18,27)] <- NA

# Where there isn't a paired reference tree, set the reference to 
# the mean of all PSCAs paired to healthy GLSEes in the same treatment.
# Don't use healthy PSCAs that are paired to unhealthy GLSEes for consistency
# with New York.

# LN
dat[dat$Species=="GLSE" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="LN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="GLSE" & dat$Use_growth_04==1 & dat$Treatment=="LN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="GLSE" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="LN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="GLSE" & dat$Use_growth_06>=0.5 & dat$Treatment=="LN",]$ref_15N_AP_u_06,na.rm=TRUE)
# MN
dat[dat$Species=="GLSE" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="MN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="GLSE" & dat$Use_growth_04==1 & dat$Treatment=="MN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="GLSE" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="MN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="GLSE" & dat$Use_growth_06>=0.5 & dat$Treatment=="MN",]$ref_15N_AP_u_06,na.rm=TRUE)
# HN
dat[dat$Species=="GLSE" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="HN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="GLSE" & dat$Use_growth_04==1 & dat$Treatment=="HN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="GLSE" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="HN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="GLSE" & dat$Use_growth_06>=0.5 & dat$Treatment=="HN",]$ref_15N_AP_u_06,na.rm=TRUE)
# PHN
dat[dat$Species=="GLSE" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="PHN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="GLSE" & dat$Use_growth_04==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="GLSE" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="PHN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="GLSE" & dat$Use_growth_06>=0.5 & dat$Treatment=="PHN",]$ref_15N_AP_u_06,na.rm=TRUE)

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
	# For _06, this is _04
# %15Nfixation is purefix_15N_AP

# 2018

dat$Ndfa_u_04 <- NA
dat[dat$Species=="GLSE",]$Ndfa_u_04 <- 
	100*((dat[dat$Species=="GLSE",]$ref_15N_AP_u_04 - pnon*purefix_15N_AP)/(1-pnon) -
	(dat[dat$Species=="GLSE",]$foliar_15N_AP_04 - pfix*purefix_15N_AP)/(1-pfix))/
	((dat[dat$Species=="GLSE",]$ref_15N_AP_u_04 - pnon*purefix_15N_AP)/(1-pnon) -
	purefix_15N_AP)

# 2019

dat$Ndfa_u_06 <- NA
dat[dat$Species=="GLSE",]$Ndfa_u_06 <- 
	100*((dat[dat$Species=="GLSE",]$ref_15N_AP_u_06 - 
	pnon*dat[dat$Species=="GLSE",]$ref_15N_AP_u_04)/(1-pnon) -
	(dat[dat$Species=="GLSE",]$foliar_15N_AP_06 - 
	pfix*dat[dat$Species=="GLSE",]$foliar_15N_AP_04)/(1-pfix))/
	((dat[dat$Species=="GLSE",]$ref_15N_AP_u_06 - 
	pnon*dat[dat$Species=="GLSE",]$ref_15N_AP_u_04)/(1-pnon) -
	purefix_15N_AP)

# Waiakea only: Filter out the ones that weren't labeled
dat[dat$Species=="GLSE",]$Ndfa_u_04[c(9,18,27)] <- NA
dat[dat$Species=="GLSE",]$Ndfa_u_06[c(9,18,27)] <- NA

#####
### Stats and print out for table
#####

# 2018-2019
# Get dataset for pNdfa and for treatment and for tree
pNdfa_u_GLSE_89 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04==1,]$Ndfa_u_04,
	dat[dat$Species=="GLSE" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06>=0.5,]$Ndfa_u_06
	)
Treatment_GLSE_89 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04==1,]$Treatment,
	dat[dat$Species=="GLSE" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06>=0.5,]$Treatment
	)
Tree_GLSE_89 <- 
	as.factor(c(
	dat[dat$Species=="GLSE" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04==1,]$PID,
	dat[dat$Species=="GLSE" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06>=0.5,]$PID
	))
summary(lme_pNdfa_retranscase <- lme(pNdfa_u_GLSE_89 ~ 0 + Treatment_GLSE_89,
	random=~1 | Tree_GLSE_89))
print(em_lme_pNdfa_retranscase <- emmeans(lme_pNdfa_retranscase, list(pairwise ~ Treatment_GLSE_89), adjust = "tukey"))
# emmeans a a a a, all > 0 (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_pNdfa_retranscase)$emmean)
print(round(pNdfa_retranscase_LN_GLSE <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(pNdfa_retranscase_MN_GLSE <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(pNdfa_retranscase_HN_GLSE <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(pNdfa_retranscase_PHN_GLSE <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

#################
#################
### Casuarina ###
#################
#################

#########################################################################
### Retranslocation Case: Casuarina uses 10% less old N than Psidium ###
#########################################################################

# Paired reference trees where possible, otherwise means of the treatment
# Retranslocation 10% less for Casuarina than for Psidium
# 2018, 2019

# Set unhealthy PSCAs foliar d15N to NA
dat[dat$Species=="PSCA" & dat$Use_growth_04==0,]$foliar_15N_AP_04 <- NA
dat[dat$Species=="PSCA" & dat$Use_growth_06==0,]$foliar_15N_AP_06 <- NA

# Reset columns for all years

dat$ref_15N_AP_u_04 <- NA
dat[dat$Species=="CAEQ",]$ref_15N_AP_u_04 <- dat[dat$Species=="PSCA",]$foliar_15N_AP_04

dat$ref_15N_AP_u_06 <- NA
dat[dat$Species=="CAEQ",]$ref_15N_AP_u_06 <- dat[dat$Species=="PSCA",]$foliar_15N_AP_06

# Waiakea only: Filter out the ones that weren't labeled
dat[dat$Species=="CAEQ",]$ref_15N_AP_u_04[c(9,18,27)] <- NA
dat[dat$Species=="CAEQ",]$ref_15N_AP_u_06[c(9,18,27)] <- NA

# Where there isn't a paired reference tree, set the reference to 
# the mean of all PSCAs paired to healthy CAEQes in the same treatment.
# Don't use healthy PSCAs that are paired to unhealthy CAEQes for consistency
# with New York.

# LN
dat[dat$Species=="CAEQ" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="LN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="CAEQ" & dat$Use_growth_04==1 & dat$Treatment=="LN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="CAEQ" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="LN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="CAEQ" & dat$Use_growth_06>=0.5 & dat$Treatment=="LN",]$ref_15N_AP_u_06,na.rm=TRUE)
# MN
dat[dat$Species=="CAEQ" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="MN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="CAEQ" & dat$Use_growth_04==1 & dat$Treatment=="MN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="CAEQ" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="MN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="CAEQ" & dat$Use_growth_06>=0.5 & dat$Treatment=="MN",]$ref_15N_AP_u_06,na.rm=TRUE)
# HN
dat[dat$Species=="CAEQ" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="HN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="CAEQ" & dat$Use_growth_04==1 & dat$Treatment=="HN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="CAEQ" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="HN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="CAEQ" & dat$Use_growth_06>=0.5 & dat$Treatment=="HN",]$ref_15N_AP_u_06,na.rm=TRUE)
# PHN
dat[dat$Species=="CAEQ" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="PHN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="CAEQ" & dat$Use_growth_04==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="CAEQ" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="PHN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="CAEQ" & dat$Use_growth_06>=0.5 & dat$Treatment=="PHN",]$ref_15N_AP_u_06,na.rm=TRUE)

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
	# For _06, this is _04
# %15Nfixation is purefix_15N_AP

# 2018

dat$Ndfa_u_04 <- NA
dat[dat$Species=="CAEQ",]$Ndfa_u_04 <- 
	100*((dat[dat$Species=="CAEQ",]$ref_15N_AP_u_04 - pnon*purefix_15N_AP)/(1-pnon) -
	(dat[dat$Species=="CAEQ",]$foliar_15N_AP_04 - pfix*purefix_15N_AP)/(1-pfix))/
	((dat[dat$Species=="CAEQ",]$ref_15N_AP_u_04 - pnon*purefix_15N_AP)/(1-pnon) -
	purefix_15N_AP)

# 2019

dat$Ndfa_u_06 <- NA
dat[dat$Species=="CAEQ",]$Ndfa_u_06 <- 
	100*((dat[dat$Species=="CAEQ",]$ref_15N_AP_u_06 - 
	pnon*dat[dat$Species=="CAEQ",]$ref_15N_AP_u_04)/(1-pnon) -
	(dat[dat$Species=="CAEQ",]$foliar_15N_AP_06 - 
	pfix*dat[dat$Species=="CAEQ",]$foliar_15N_AP_04)/(1-pfix))/
	((dat[dat$Species=="CAEQ",]$ref_15N_AP_u_06 - 
	pnon*dat[dat$Species=="CAEQ",]$ref_15N_AP_u_04)/(1-pnon) -
	purefix_15N_AP)

# Waiakea only: Filter out the ones that weren't labeled
dat[dat$Species=="CAEQ",]$Ndfa_u_04[c(9,18,27)] <- NA
dat[dat$Species=="CAEQ",]$Ndfa_u_06[c(9,18,27)] <- NA

#####
### Stats and print out for table
#####

# 2018-2019
# Get dataset for pNdfa and for treatment and for tree
pNdfa_u_CAEQ_89 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04==1,]$Ndfa_u_04,
	dat[dat$Species=="CAEQ" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06>=0.5,]$Ndfa_u_06
	)
Treatment_CAEQ_89 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04==1,]$Treatment,
	dat[dat$Species=="CAEQ" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06>=0.5,]$Treatment
	)
Tree_CAEQ_89 <- 
	as.factor(c(
	dat[dat$Species=="CAEQ" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04==1,]$PID,
	dat[dat$Species=="CAEQ" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06>=0.5,]$PID
	))
summary(lme_pNdfa_retranscase <- lme(pNdfa_u_CAEQ_89 ~ 0 + Treatment_CAEQ_89,
	random=~1 | Tree_CAEQ_89))
print(em_lme_pNdfa_retranscase <- emmeans(lme_pNdfa_retranscase, list(pairwise ~ Treatment_CAEQ_89), adjust = "tukey"))
# emmeans a ab ab b, all > 0 (based on 95% CIs)
tab <- as.data.frame(summary(em_lme_pNdfa_retranscase)$emmean)
print(round(pNdfa_retranscase_LN_CAEQ <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(pNdfa_retranscase_MN_CAEQ <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(pNdfa_retranscase_HN_CAEQ <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(pNdfa_retranscase_PHN_CAEQ <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

####################################################################
####################################################################
####################################################################