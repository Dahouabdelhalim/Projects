# NSF FX N fixation analysis for the New York 15N data

# This script calculates N fixation estimates (%Ndfa; total Ndfa is in 
# another file) for the trees at Black Rock from the Menge, Funk, Perakis, 
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
	# - 2017
	# - 2018
	# - 2019
# - Assumptions about retranslocation
	# - Base case (this file): 
		# Retranslocation is similar for Betula and Robinia
	# - Retranslocation sensitivity check (different file):
		# Betula uses more retranslocated N than Robinia

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

# Run the growth analysis code to load the dataframe and calculate growth rates

source("NY_growth_analysis.R")

#################################################################################
####################### N fixation estimates (%Ndfa) ############################
#################################################################################

purefix_15N_AP <- 0.3663

##################
### Base case: ###
##################

# Paired reference trees where possible, otherwise means of the treatment
# Retranslocation similar for Betula and Robinia
# 2016, 2017, 2018, 2019

# Set unhealthy BENIs foliar d15N to NA
dat[dat$Species=="BENI" & dat$Use_growth_02==0,]$foliar_15N_AP_02 <- NA
dat[dat$Species=="BENI" & dat$Use_growth_04==0,]$foliar_15N_AP_04 <- NA
dat[dat$Species=="BENI" & dat$Use_growth_07==0,]$foliar_15N_AP_07 <- NA
dat[dat$Species=="BENI" & dat$Use_growth_08==0,]$foliar_15N_AP_08 <- NA

# Set up columns for all years and set the reference to the paired tree

dat$ref_15N_AP_u_02 <- NA
dat[dat$Species=="ROPS",]$ref_15N_AP_u_02 <- dat[dat$Species=="BENI",]$foliar_15N_AP_02

dat$ref_15N_AP_u_04 <- NA
dat[dat$Species=="ROPS",]$ref_15N_AP_u_04 <- dat[dat$Species=="BENI",]$foliar_15N_AP_04

dat$ref_15N_AP_u_07 <- NA
dat[dat$Species=="ROPS",]$ref_15N_AP_u_07 <- dat[dat$Species=="BENI",]$foliar_15N_AP_07

dat$ref_15N_AP_u_08 <- NA
dat[dat$Species=="ROPS",]$ref_15N_AP_u_08 <- dat[dat$Species=="BENI",]$foliar_15N_AP_08

# Where there isn't a paired reference tree, set the reference to 
# the mean of all BENIs paired to healthy ROPSes in the same treatment.
# Don't use healthy BENIs that are paired to unhealthy ROPSes because
# we stopped labeling them.

# LN
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_02) & dat$Treatment=="LN",]$ref_15N_AP_u_02 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_02==1 & dat$Treatment=="LN",]$ref_15N_AP_u_02,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="LN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_04==1 & dat$Treatment=="LN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_07) & dat$Treatment=="LN",]$ref_15N_AP_u_07 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_07==1 & dat$Treatment=="LN",]$ref_15N_AP_u_07,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_08) & dat$Treatment=="LN",]$ref_15N_AP_u_08 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_08==1 & dat$Treatment=="LN",]$ref_15N_AP_u_08,na.rm=TRUE)
# MN
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_02) & dat$Treatment=="MN",]$ref_15N_AP_u_02 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_02==1 & dat$Treatment=="MN",]$ref_15N_AP_u_02,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="MN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_04==1 & dat$Treatment=="MN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_07) & dat$Treatment=="MN",]$ref_15N_AP_u_07 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_07==1 & dat$Treatment=="MN",]$ref_15N_AP_u_07,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_08) & dat$Treatment=="MN",]$ref_15N_AP_u_08 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_08==1 & dat$Treatment=="MN",]$ref_15N_AP_u_08,na.rm=TRUE)
# HN
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_02) & dat$Treatment=="HN",]$ref_15N_AP_u_02 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_02==1 & dat$Treatment=="HN",]$ref_15N_AP_u_02,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="HN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_04==1 & dat$Treatment=="HN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_07) & dat$Treatment=="HN",]$ref_15N_AP_u_07 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_07==1 & dat$Treatment=="HN",]$ref_15N_AP_u_07,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_08) & dat$Treatment=="HN",]$ref_15N_AP_u_08 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_08==1 & dat$Treatment=="HN",]$ref_15N_AP_u_08,na.rm=TRUE)
# PHN
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_02) & dat$Treatment=="PHN",]$ref_15N_AP_u_02 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_02==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_02,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_04) & dat$Treatment=="PHN",]$ref_15N_AP_u_04 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_04==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_04,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_07) & dat$Treatment=="PHN",]$ref_15N_AP_u_07 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_07==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_07,na.rm=TRUE)
dat[dat$Species=="ROPS" & is.na(dat$ref_15N_AP_u_08) & dat$Treatment=="PHN",]$ref_15N_AP_u_08 <-
	mean(dat[dat$Species=="ROPS" & dat$Use_growth_08==1 & dat$Treatment=="PHN",]$ref_15N_AP_u_08,na.rm=TRUE)

#####
### Calculate %Ndfa
#####

# 2016

dat$Ndfa_u_02 <- NA
dat[dat$Species=="ROPS",]$Ndfa_u_02 <- 
	100*(dat[dat$Species=="ROPS",]$ref_15N_AP_u_02 - 
	dat[dat$Species=="ROPS",]$foliar_15N_AP_02)/
	(dat[dat$Species=="ROPS",]$ref_15N_AP_u_02 - purefix_15N_AP)

# 2017 

dat$Ndfa_u_04 <- NA
dat[dat$Species=="ROPS",]$Ndfa_u_04 <- 
	100*(dat[dat$Species=="ROPS",]$ref_15N_AP_u_04 - 
	dat[dat$Species=="ROPS",]$foliar_15N_AP_04)/
	(dat[dat$Species=="ROPS",]$ref_15N_AP_u_04 - purefix_15N_AP)

# 2018

dat$Ndfa_u_07 <- NA
dat[dat$Species=="ROPS",]$Ndfa_u_07 <- 
	100*(dat[dat$Species=="ROPS",]$ref_15N_AP_u_07 - 
	dat[dat$Species=="ROPS",]$foliar_15N_AP_07)/
	(dat[dat$Species=="ROPS",]$ref_15N_AP_u_07 - purefix_15N_AP)

# 2019

dat$Ndfa_u_08 <- NA
dat[dat$Species=="ROPS",]$Ndfa_u_08 <- 
	100*(dat[dat$Species=="ROPS",]$ref_15N_AP_u_08 - 
	dat[dat$Species=="ROPS",]$foliar_15N_AP_08)/
	(dat[dat$Species=="ROPS",]$ref_15N_AP_u_08 - purefix_15N_AP)

#####
### Stats and print out for table
#####

# 2017-2019
# Get dataset for pNdfa and for treatment and for tree
pNdfa_u_789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$Ndfa_u_04,
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_07) & 
		dat$Use_growth_07 == 1,]$Ndfa_u_07,
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_08) & 
		dat$Use_growth_08 == 1,]$Ndfa_u_08
	)
Treatment_ROPS_789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_07) & 
		dat$Use_growth_07 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_08) & 
		dat$Use_growth_08 == 1,]$Treatment
	)
Tree_ROPS_789 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_07) & 
		dat$Use_growth_07 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_08) & 
		dat$Use_growth_08 == 1,]$PID
	))
summary(lme_pNdfa_base <- lme(pNdfa_u_789 ~ 0 + Treatment_ROPS_789,
	random=~1 | Tree_ROPS_789))
print(em_lme_pNdfa_base <- emmeans(lme_pNdfa_base, list(pairwise ~ Treatment_ROPS_789), adjust = "tukey"))
# emmeans: a b b b, all > 0 (based on 95% CI)
tab <- as.data.frame(summary(em_lme_pNdfa_base)$emmean)
print(round(pNdfa_base_LN <- c(tab[2,2],tab[2,2]-tab[2,3],tab[2,2]+tab[2,3]),1))
print(round(pNdfa_base_MN <- c(tab[3,2],tab[3,2]-tab[3,3],tab[3,2]+tab[3,3]),1))
print(round(pNdfa_base_HN <- c(tab[1,2],tab[1,2]-tab[1,3],tab[1,2]+tab[1,3]),1))
print(round(pNdfa_base_PHN <- c(tab[4,2],tab[4,2]-tab[4,3],tab[4,2]+tab[4,3]),1))

####################################################################
####################################################################
####################################################################