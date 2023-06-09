# NSF FX N fixation analysis for the Volcano 15N data

# This script calculates N fixation estimates (%Ndfa; total Ndfa is in 
# another file) for the trees at the Volcano site from the Menge, Funk, Perakis, 
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
# - Assumptions about retranslocation
	# - Base case (this file): 
		# Retranslocation is similar for Dodonaea and Acacia and Morella
	# - Retranslocation sensitivity check (different file):
		# Dodonaea uses more retranslocated N than Acacia and Morella

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Run the growth analysis code to load the dataframe and calculate growth rates

source("HI_V_growth_analysis.R")

#################################################################################
####################### N fixation estimates (%Ndfa) ############################
#################################################################################

purefix_15N_AP <- 0.3663

##################
##################
##### Acacia #####
##################
##################

##################
### Base case: ###
##################

# Paired reference trees where possible, otherwise means of the treatment
# Retranslocation similar for Dodonaea and Acacia
# 2018, 2019

# Set unhealthy DOVIs foliar d15N to NA
dat[dat$Species=="DOVI" & dat$Use_growth_04==0,]$foliar_15N_AP_04 <- NA
dat[dat$Species=="DOVI" & dat$Use_growth_06==0,]$foliar_15N_AP_06 <- NA

# Set up columns for all years and set the reference to the paired tree

dat$ref_15N_AP_u_04 <- NA
dat[dat$Species=="ACKO",]$ref_15N_AP_u_04 <- dat[dat$Species=="DOVI",]$foliar_15N_AP_04

dat$ref_15N_AP_u_06 <- NA
dat[dat$Species=="ACKO",]$ref_15N_AP_u_06 <- dat[dat$Species=="DOVI",]$foliar_15N_AP_06

# Where there isn't a paired reference tree, set the reference to 
# the mean of all DOVIs paired to healthy ACKOs in the same treatment.
# Don't use healthy DOVIs that are paired to unhealthy ACKOs for consistency
# with New York.

# LN
dat[dat$Species=="ACKO" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="LN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="ACKO" & dat$Use_growth_04==1 & dat$Treatment=="LN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="ACKO" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="LN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="ACKO" & dat$Use_growth_06==1 & dat$Treatment=="LN",]$ref_15N_AP_u_06,na.rm=TRUE)
# MN
dat[dat$Species=="ACKO" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="MN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="ACKO" & dat$Use_growth_04==1 & dat$Treatment=="MN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="ACKO" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="MN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="ACKO" & dat$Use_growth_06==1 & dat$Treatment=="MN",]$ref_15N_AP_u_06,na.rm=TRUE)
# HN
dat[dat$Species=="ACKO" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="HN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="ACKO" & dat$Use_growth_04==1 & dat$Treatment=="HN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="ACKO" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="HN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="ACKO" & dat$Use_growth_06==1 & dat$Treatment=="HN",]$ref_15N_AP_u_06,na.rm=TRUE)
# PHN
dat[dat$Species=="ACKO" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="PHN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="ACKO" & dat$Use_growth_04==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="ACKO" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="PHN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="ACKO" & dat$Use_growth_06==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_06,na.rm=TRUE)

#####
### Calculate %Ndfa
#####

# 2018

dat$Ndfa_u_04 <- NA
dat[dat$Species=="ACKO",]$Ndfa_u_04 <- 
	100*(dat[dat$Species=="ACKO",]$ref_15N_AP_u_04 - 
	dat[dat$Species=="ACKO",]$foliar_15N_AP_04)/
	(dat[dat$Species=="ACKO",]$ref_15N_AP_u_04 - purefix_15N_AP)

# 2019

dat$Ndfa_u_06 <- NA
dat[dat$Species=="ACKO",]$Ndfa_u_06 <- 
	100*(dat[dat$Species=="ACKO",]$ref_15N_AP_u_06 - 
	dat[dat$Species=="ACKO",]$foliar_15N_AP_06)/
	(dat[dat$Species=="ACKO",]$ref_15N_AP_u_06 - purefix_15N_AP)

#####
### Stats and print out for table
#####

# 2018-2019
# Get dataset for pNdfa and for treatment and for tree
pNdfa_u_ACKO_89 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$Ndfa_u_04,
	dat[dat$Species=="ACKO" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$Ndfa_u_06
	)
Treatment_ACKO_89 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ACKO" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Tree_ACKO_89 <- 
	as.factor(c(
	dat[dat$Species=="ACKO" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ACKO" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_pNdfa_base <- lme(pNdfa_u_ACKO_89 ~ 0 + Treatment_ACKO_89,
	random=~1 | Tree_ACKO_89))
print(em_lme_pNdfa_base <- emmeans(lme_pNdfa_base, list(pairwise ~ Treatment_ACKO_89), adjust = "tukey"))
# emmeans a b ab ab, all > 0 (based on 95% CI)
tab <- as.data.frame(summary(em_lme_pNdfa_base)$emmean)
print(round(pNdfa_base_LN_ACKO <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(pNdfa_base_MN_ACKO <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(pNdfa_base_HN_ACKO <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(pNdfa_base_PHN_ACKO <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))



#################
#################
### Morella ###
#################
#################

##################
### Base case: ###
##################

# Paired reference trees where possible, otherwise means of the treatment
# Retranslocation similar for Dodonaea and Morella
# 2018, 2019

# Set unhealthy DOVIs foliar d15N to NA
dat[dat$Species=="DOVI" & dat$Use_growth_04==0,]$foliar_15N_AP_04 <- NA
dat[dat$Species=="DOVI" & dat$Use_growth_06==0,]$foliar_15N_AP_06 <- NA

# Set up columns for all years and set the reference to the paired tree

dat[dat$Species=="MOFA",]$ref_15N_AP_u_04 <- dat[dat$Species=="DOVI",]$foliar_15N_AP_04

dat[dat$Species=="MOFA",]$ref_15N_AP_u_06 <- dat[dat$Species=="DOVI",]$foliar_15N_AP_06

# Where there isn't a paired reference tree, set the reference to 
# the mean of all DOVIs paired to healthy MOFAes in the same treatment.
# Don't use healthy DOVIs that are paired to unhealthy MOFAes for consistency
# with New York.

# LN
dat[dat$Species=="MOFA" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="LN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="MOFA" & dat$Use_growth_04==1 & dat$Treatment=="LN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="MOFA" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="LN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="MOFA" & dat$Use_growth_06==1 & dat$Treatment=="LN",]$ref_15N_AP_u_06,na.rm=TRUE)
# MN
dat[dat$Species=="MOFA" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="MN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="MOFA" & dat$Use_growth_04==1 & dat$Treatment=="MN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="MOFA" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="MN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="MOFA" & dat$Use_growth_06==1 & dat$Treatment=="MN",]$ref_15N_AP_u_06,na.rm=TRUE)
# HN
dat[dat$Species=="MOFA" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="HN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="MOFA" & dat$Use_growth_04==1 & dat$Treatment=="HN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="MOFA" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="HN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="MOFA" & dat$Use_growth_06==1 & dat$Treatment=="HN",]$ref_15N_AP_u_06,na.rm=TRUE)
# PHN
dat[dat$Species=="MOFA" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="PHN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="MOFA" & dat$Use_growth_04==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="MOFA" & is.na(dat$ref_15N_AP_u_06) & dat$Treatment=="PHN",]$ref_15N_AP_u_06 <-
	mean(dat[dat$Species=="MOFA" & dat$Use_growth_06==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_06,na.rm=TRUE)

#####
### Calculate %Ndfa
#####

# 2018

#dat$Ndfa_u_04 <- NA
dat[dat$Species=="MOFA",]$Ndfa_u_04 <- 
	100*(dat[dat$Species=="MOFA",]$ref_15N_AP_u_04 - 
	dat[dat$Species=="MOFA",]$foliar_15N_AP_04)/
	(dat[dat$Species=="MOFA",]$ref_15N_AP_u_04 - purefix_15N_AP)

# 2019

#dat$Ndfa_u_06 <- NA
dat[dat$Species=="MOFA",]$Ndfa_u_06 <- 
	100*(dat[dat$Species=="MOFA",]$ref_15N_AP_u_06 - 
	dat[dat$Species=="MOFA",]$foliar_15N_AP_06)/
	(dat[dat$Species=="MOFA",]$ref_15N_AP_u_06 - purefix_15N_AP)

#####
### Stats and print out for table
#####

# 2018-2019
# Get dataset for pNdfa and for treatment and for tree
pNdfa_u_MOFA_89 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$Ndfa_u_04,
	dat[dat$Species=="MOFA" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$Ndfa_u_06
	)
Treatment_MOFA_89 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="MOFA" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Tree_MOFA_89 <- 
	as.factor(c(
	dat[dat$Species=="MOFA" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="MOFA" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_pNdfa_base <- lme(pNdfa_u_MOFA_89 ~ 0 + Treatment_MOFA_89,
	random=~1 | Tree_MOFA_89))
print(em_lme_pNdfa_base <- emmeans(lme_pNdfa_base, list(pairwise ~ Treatment_MOFA_89), adjust = "tukey"))
# emmeans a a a a, all > 0 except HN (based on 95% CI)
tab <- as.data.frame(summary(em_lme_pNdfa_base)$emmean)
print(round(pNdfa_base_LN_MOFA <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(pNdfa_base_MN_MOFA <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(pNdfa_base_HN_MOFA <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(pNdfa_base_PHN_MOFA <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

####################################################################
####################################################################
####################################################################