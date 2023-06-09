# NSF FX growth analysis for the New York biomass data

# There are a number of ways we could analyze biomass growth rates
# for these trees. Each of these will be separate for each species.
# - Compute absolute and relative growth rates for each time period
# - Compute absolute and relative growth rates for each year
# - Compute absolute and relative growth rates from t0 to each tf
# For each of these, we want to analyze growth as a function of 
# treatment, which answers our initial scientific question:
# Is the species limited by N and/or P at each treatment level?
# We also want to account for the fact that growth rate varies with
# tree size, possibly individual trees, and possibly with year.
# Also want to do these for each year: Up to 2016, 2017, 2018, 2019.

# The "Use_growth_0X" columns in "NY_FX_Size_FoliarCNIsotope_Data.csv" 
# are the codes for whether or not to use data.
# 1 means use.
# 0 means dead or sufficiently ill or damaged that it shouldn't be used.

rm(list=ls())

library(lme4)
library(nlme)
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

dat <- read.csv("NY_FX_Size_FoliarCNIsotope_Data.csv")[1:66,1:166]

# Option to print out how many are in each category
# (uncomment the following lines)

#print(dat[dat$Species=="BENI" & dat$Treatment=="LN",c(32,64,99,123)])
#print(dat[dat$Species=="BENI" & dat$Treatment=="MN",c(32,64,99,123)])
#print(dat[dat$Species=="BENI" & dat$Treatment=="HN",c(32,64,99,123)])
#print(dat[dat$Species=="BENI" & dat$Treatment=="PHN",c(32,64,99,123)])

#print(dat[dat$Species=="ROPS" & dat$Treatment=="LN",c(32,64,99,123)])
#print(dat[dat$Species=="ROPS" & dat$Treatment=="MN",c(32,64,99,123)])
#print(dat[dat$Species=="ROPS" & dat$Treatment=="HN",c(32,64,99,123)])
#print(dat[dat$Species=="ROPS" & dat$Treatment=="PHN",c(32,64,99,123)])

#############################################################
################ Calculate growth rates #####################
#############################################################

# "AGR" is absolute growth rate and "RGR" is relative growth rate. 
# Each is calculated for aboveground, belowground, and total biomass, 
# for the following time periods: 

# Annual increments:

# Biomass

# June 2015 to October 2016
dat$AGR_Biomass_2015_2016 <- (dat$Biomass_est_kg_02 - dat$Biomass_est_kg_01)/((dat$Days_02 - dat$Days_01)/365.25)
dat$RGR_Biomass_2015_2016 <- (log(dat$Biomass_est_kg_02) - log(dat$Biomass_est_kg_01))/((dat$Days_02 - dat$Days_01)/365.25)

# October 2016 to October 2017
dat$AGR_Biomass_2016_2017 <- (dat$Biomass_est_kg_04 - dat$Biomass_est_kg_02)/((dat$Days_04 - dat$Days_02)/365.25)
dat$RGR_Biomass_2016_2017 <- (log(dat$Biomass_est_kg_04) - log(dat$Biomass_est_kg_02))/((dat$Days_04 - dat$Days_02)/365.25)

# October 2017 to September 2018
dat$AGR_Biomass_2017_2018 <- (dat$Biomass_est_kg_07 - dat$Biomass_est_kg_04)/((dat$Days_07 - dat$Days_04)/365.25)
dat$RGR_Biomass_2017_2018 <- (log(dat$Biomass_est_kg_07) - log(dat$Biomass_est_kg_04))/((dat$Days_07 - dat$Days_04)/365.25)

# September 2018 to October 2019
dat$AGR_Biomass_2018_2019 <- (dat$Biomass_est_kg_08 - dat$Biomass_est_kg_07)/((dat$Days_08 - dat$Days_07)/365.25)
dat$RGR_Biomass_2018_2019 <- (log(dat$Biomass_est_kg_08) - log(dat$Biomass_est_kg_07))/((dat$Days_08 - dat$Days_07)/365.25)

# AGB

# June 2015 to October 2016
dat$AGR_AGB_2015_2016 <- (dat$AGB_est_kg_02 - dat$AGB_est_kg_01)/((dat$Days_02 - dat$Days_01)/365.25)
dat$RGR_AGB_2015_2016 <- (log(dat$AGB_est_kg_02) - log(dat$AGB_est_kg_01))/((dat$Days_02 - dat$Days_01)/365.25)

# October 2016 to October 2017
dat$AGR_AGB_2016_2017 <- (dat$AGB_est_kg_04 - dat$AGB_est_kg_02)/((dat$Days_04 - dat$Days_02)/365.25)
dat$RGR_AGB_2016_2017 <- (log(dat$AGB_est_kg_04) - log(dat$AGB_est_kg_02))/((dat$Days_04 - dat$Days_02)/365.25)

# October 2017 to September 2018
dat$AGR_AGB_2017_2018 <- (dat$AGB_est_kg_07 - dat$AGB_est_kg_04)/((dat$Days_07 - dat$Days_04)/365.25)
dat$RGR_AGB_2017_2018 <- (log(dat$AGB_est_kg_07) - log(dat$AGB_est_kg_04))/((dat$Days_07 - dat$Days_04)/365.25)

# September 2018 to October 2019
dat$AGR_AGB_2018_2019 <- (dat$AGB_est_kg_08 - dat$AGB_est_kg_07)/((dat$Days_08 - dat$Days_07)/365.25)
dat$RGR_AGB_2018_2019 <- (log(dat$AGB_est_kg_08) - log(dat$AGB_est_kg_07))/((dat$Days_08 - dat$Days_07)/365.25)

# BGB

# June 2015 to October 2016
dat$AGR_BGB_2015_2016 <- (dat$BGB_est_kg_02 - dat$BGB_est_kg_01)/((dat$Days_02 - dat$Days_01)/365.25)
dat$RGR_BGB_2015_2016 <- (log(dat$BGB_est_kg_02) - log(dat$BGB_est_kg_01))/((dat$Days_02 - dat$Days_01)/365.25)

# October 2016 to October 2017
dat$AGR_BGB_2016_2017 <- (dat$BGB_est_kg_04 - dat$BGB_est_kg_02)/((dat$Days_04 - dat$Days_02)/365.25)
dat$RGR_BGB_2016_2017 <- (log(dat$BGB_est_kg_04) - log(dat$BGB_est_kg_02))/((dat$Days_04 - dat$Days_02)/365.25)

# October 2017 to September 2018
dat$AGR_BGB_2017_2018 <- (dat$BGB_est_kg_07 - dat$BGB_est_kg_04)/((dat$Days_07 - dat$Days_04)/365.25)
dat$RGR_BGB_2017_2018 <- (log(dat$BGB_est_kg_07) - log(dat$BGB_est_kg_04))/((dat$Days_07 - dat$Days_04)/365.25)

# September 2018 to October 2019
dat$AGR_BGB_2018_2019 <- (dat$BGB_est_kg_08 - dat$BGB_est_kg_07)/((dat$Days_08 - dat$Days_07)/365.25)
dat$RGR_BGB_2018_2019 <- (log(dat$BGB_est_kg_08) - log(dat$BGB_est_kg_07))/((dat$Days_08 - dat$Days_07)/365.25)



# 2015 to tf increments:

# Biomass

# June 2015 to October 2017
dat$AGR_Biomass_2015_2017 <- (dat$Biomass_est_kg_04 - dat$Biomass_est_kg_01)/((dat$Days_04 - dat$Days_01)/365.25)
dat$RGR_Biomass_2015_2017 <- (log(dat$Biomass_est_kg_04) - log(dat$Biomass_est_kg_01))/((dat$Days_04 - dat$Days_01)/365.25)

# June 2015 to September 2018
dat$AGR_Biomass_2015_2018 <- (dat$Biomass_est_kg_07 - dat$Biomass_est_kg_01)/((dat$Days_07 - dat$Days_01)/365.25)
dat$RGR_Biomass_2015_2018 <- (log(dat$Biomass_est_kg_07) - log(dat$Biomass_est_kg_01))/((dat$Days_07 - dat$Days_01)/365.25)

# June 2015 to October 2019
dat$AGR_Biomass_2015_2019 <- (dat$Biomass_est_kg_08 - dat$Biomass_est_kg_01)/((dat$Days_08 - dat$Days_01)/365.25)
dat$RGR_Biomass_2015_2019 <- (log(dat$Biomass_est_kg_08) - log(dat$Biomass_est_kg_01))/((dat$Days_08 - dat$Days_01)/365.25)

# AGB

# June 2015 to October 2017
dat$AGR_AGB_2015_2017 <- (dat$AGB_est_kg_04 - dat$AGB_est_kg_01)/((dat$Days_04 - dat$Days_01)/365.25)
dat$RGR_AGB_2015_2017 <- (log(dat$AGB_est_kg_04) - log(dat$AGB_est_kg_01))/((dat$Days_04 - dat$Days_01)/365.25)

# June 2015 to September 2018
dat$AGR_AGB_2015_2018 <- (dat$AGB_est_kg_07 - dat$AGB_est_kg_01)/((dat$Days_07 - dat$Days_01)/365.25)
dat$RGR_AGB_2015_2018 <- (log(dat$AGB_est_kg_07) - log(dat$AGB_est_kg_01))/((dat$Days_07 - dat$Days_01)/365.25)

# June 2015 to October 2019
dat$AGR_AGB_2015_2019 <- (dat$AGB_est_kg_08 - dat$AGB_est_kg_01)/((dat$Days_08 - dat$Days_01)/365.25)
dat$RGR_AGB_2015_2019 <- (log(dat$AGB_est_kg_08) - log(dat$AGB_est_kg_01))/((dat$Days_08 - dat$Days_01)/365.25)

# BGB

# June 2015 to October 2017
dat$AGR_BGB_2015_2017 <- (dat$BGB_est_kg_04 - dat$BGB_est_kg_01)/((dat$Days_04 - dat$Days_01)/365.25)
dat$RGR_BGB_2015_2017 <- (log(dat$BGB_est_kg_04) - log(dat$BGB_est_kg_01))/((dat$Days_04 - dat$Days_01)/365.25)

# June 2015 to September 2018
dat$AGR_BGB_2015_2018 <- (dat$BGB_est_kg_07 - dat$BGB_est_kg_01)/((dat$Days_07 - dat$Days_01)/365.25)
dat$RGR_BGB_2015_2018 <- (log(dat$BGB_est_kg_07) - log(dat$BGB_est_kg_01))/((dat$Days_07 - dat$Days_01)/365.25)

# June 2015 to October 2019
dat$AGR_BGB_2015_2019 <- (dat$BGB_est_kg_08 - dat$BGB_est_kg_01)/((dat$Days_08 - dat$Days_01)/365.25)
dat$RGR_BGB_2015_2019 <- (log(dat$BGB_est_kg_08) - log(dat$BGB_est_kg_01))/((dat$Days_08 - dat$Days_01)/365.25)


###############################################################
########################### Stats #############################
###############################################################

# Only use trees/measurements with Use_growth_0X == 1
# Response variables:
# BENI and ROPS ...
# Total biomass, AGB, and BGB for ...
# 2016: 
# - AGR 2015-2016
# - RGR 2015-2016
# 2017: 
# - AGR 2015-2016, 2016-2017
# - RGR 2015-2016, 2016-2017
# - AGR 2015-2017
# - RGR 2015-2017
# 2018: 
# - AGR 2015-2016, 2016-2017, 2017-2018
# - RGR 2015-2016, 2016-2017, 2017-2018
# - AGR 2015-2018
# - RGR 2015-2018
# 2019: 
# - AGR 2015-2016, 2016-2017, 2017-2018, 2018-2019
# - RGR 2015-2016, 2016-2017, 2017-2018, 2018-2019
# - AGR 2015-2019
# - RGR 2015-2019
#
# Driver: Treatment
# Covariates: For multi-year, biomass/AGB/BGB (fixed), tree (random)

########
# BENI #
########

### 
# BENI Total Biomass
###

##### RGR

### 2015-2016

summary(lm_BENI_Biomass_RGR_2015_2016 <- lm(RGR_Biomass_2015_2016 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_BENI_Biomass_RGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

### 2015-2017

# Cumulative RGR
summary(lm_BENI_Biomass_RGR_2015_2017 <- lm(RGR_Biomass_2015_2017 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_BENI_Biomass_RGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

# Annual RGR
RGR_BENI_Biomass_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_Biomass_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2016_2017
	)
Treatment_BENI_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
Biomass_BENI_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_02
	)
Tree_BENI_567 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_BENI_Biomass_RGR_201567 <- lme(RGR_BENI_Biomass_567 ~ 
	Treatment_BENI_567 + Biomass_BENI_567,
	random=~1 | Tree_BENI_567))
anova(lme_BENI_Biomass_RGR_201567)
emmeans(lme_BENI_Biomass_RGR_201567, list(pairwise ~ Treatment_BENI_567), adjust = "tukey")
# a b b ab

### 2015-2018

# Cumulative RGR
summary(lm_BENI_Biomass_RGR_2015_2018 <- lm(RGR_Biomass_2015_2018 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_BENI_Biomass_RGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

# Annual RGR
RGR_BENI_Biomass_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_Biomass_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_Biomass_2017_2018
	)
Treatment_BENI_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment
	)
Biomass_BENI_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_02,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$Biomass_est_kg_04
	)
Tree_BENI_5678 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID
	))
summary(lme_BENI_Biomass_RGR_2015678 <- lme(RGR_BENI_Biomass_5678 ~ 
	Treatment_BENI_5678 + Biomass_BENI_5678,
	random=~1 | Tree_BENI_5678))
anova(lme_BENI_Biomass_RGR_2015678)
emmeans(lme_BENI_Biomass_RGR_2015678, list(pairwise ~ Treatment_BENI_5678), adjust = "tukey")
# a b b ab

### 2015-2019

# Cumulative RGR
summary(lm_BENI_Biomass_RGR_2015_2019 <- lm(RGR_Biomass_2015_2019 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_BENI_Biomass_RGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_BENI_Biomass_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_Biomass_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_Biomass_2017_2018,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$RGR_Biomass_2018_2019
	)
Treatment_BENI_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$Treatment
	)
Biomass_BENI_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_02,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$Biomass_est_kg_04,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$Biomass_est_kg_07
	)
Tree_BENI_56789 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$PID
	))
summary(lme_BENI_Biomass_RGR_20156789 <- lme(RGR_BENI_Biomass_56789 ~ 
	Treatment_BENI_56789 + Biomass_BENI_56789,
	random=~1 | Tree_BENI_56789))
anova(lme_BENI_Biomass_RGR_20156789)
emmeans(lme_BENI_Biomass_RGR_20156789, list(pairwise ~ Treatment_BENI_56789), adjust = "tukey")
# a b b b

##### AGR

### 2015-2016

summary(lm_BENI_Biomass_AGR_2015_2016 <- lm(log(AGR_Biomass_2015_2016) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_BENI_Biomass_AGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2015-2017

# Cumulative AGR
summary(lm_BENI_Biomass_AGR_2015_2017 <- lm(log(AGR_Biomass_2015_2017) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_BENI_Biomass_AGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a ab b ab

# Annual AGR
AGR_BENI_Biomass_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_Biomass_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2016_2017
	)
summary(lme_BENI_Biomass_AGR_201567 <- lme(log(AGR_BENI_Biomass_567) ~ 
	Treatment_BENI_567 + Biomass_BENI_567,
	random=~1 | Tree_BENI_567))
anova(lme_BENI_Biomass_AGR_201567)
emmeans(lme_BENI_Biomass_AGR_201567, list(pairwise ~ Treatment_BENI_567), adjust = "tukey")
# a b b b

### 2015-2018

# Cumulative AGR
summary(lm_BENI_Biomass_AGR_2015_2018 <- lm(log(AGR_Biomass_2015_2018) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_BENI_Biomass_AGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a ab b ab

# Annual AGR
AGR_BENI_Biomass_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_Biomass_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_Biomass_2017_2018
	)
summary(lme_BENI_Biomass_AGR_2015678 <- lme(log(AGR_BENI_Biomass_5678) ~ 
	Treatment_BENI_5678 + Biomass_BENI_5678,
	random=~1 | Tree_BENI_5678))
anova(lme_BENI_Biomass_AGR_2015678)
emmeans(lme_BENI_Biomass_AGR_2015678, list(pairwise ~ Treatment_BENI_5678), adjust = "tukey")
# a b b b

### 2015-2019

# Cumulative AGR
summary(lm_BENI_Biomass_AGR_2015_2019 <- lm(log(AGR_Biomass_2015_2019) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_BENI_Biomass_AGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_BENI_Biomass_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_Biomass_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_Biomass_2017_2018,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$AGR_Biomass_2018_2019
	)
summary(lme_BENI_Biomass_AGR_20156789 <- lme(log(AGR_BENI_Biomass_56789) ~ 
	Treatment_BENI_56789 + Biomass_BENI_56789,
	random=~1 | Tree_BENI_56789))
anova(lme_BENI_Biomass_AGR_20156789)
emmeans(lme_BENI_Biomass_AGR_20156789, list(pairwise ~ Treatment_BENI_56789), adjust = "tukey")
# a b b ab






### 
# BENI Aboveground Biomass
###

##### RGR

### 2015-2016

summary(lm_BENI_AGB_RGR_2015_2016 <- lm(RGR_AGB_2015_2016 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_BENI_AGB_RGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

### 2015-2017

# Cumulative RGR
summary(lm_BENI_AGB_RGR_2015_2017 <- lm(RGR_AGB_2015_2017 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_BENI_AGB_RGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

# Annual RGR
RGR_BENI_AGB_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_AGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2016_2017
	)
Treatment_BENI_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
AGB_BENI_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_02
	)
Tree_BENI_567 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_BENI_AGB_RGR_201567 <- lme(RGR_BENI_AGB_567 ~ 
	Treatment_BENI_567 + AGB_BENI_567,
	random=~1 | Tree_BENI_567))
anova(lme_BENI_AGB_RGR_201567)
emmeans(lme_BENI_AGB_RGR_201567, list(pairwise ~ Treatment_BENI_567), adjust = "tukey")
# a b b ab

### 2015-2018

# Cumulative RGR
summary(lm_BENI_AGB_RGR_2015_2018 <- lm(RGR_AGB_2015_2018 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_BENI_AGB_RGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

# Annual RGR
RGR_BENI_AGB_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_AGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_AGB_2017_2018
	)
Treatment_BENI_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment
	)
AGB_BENI_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_02,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGB_est_kg_04
	)
Tree_BENI_5678 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID
	))
summary(lme_BENI_AGB_RGR_2015678 <- lme(RGR_BENI_AGB_5678 ~ 
	Treatment_BENI_5678 + AGB_BENI_5678,
	random=~1 | Tree_BENI_5678))
anova(lme_BENI_AGB_RGR_2015678)
emmeans(lme_BENI_AGB_RGR_2015678, list(pairwise ~ Treatment_BENI_5678), adjust = "tukey")
# a b b ab

### 2015-2019

# Cumulative RGR
summary(lm_BENI_AGB_RGR_2015_2019 <- lm(RGR_AGB_2015_2019 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_BENI_AGB_RGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_BENI_AGB_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_AGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_AGB_2017_2018,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$RGR_AGB_2018_2019
	)
Treatment_BENI_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$Treatment
	)
AGB_BENI_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_02,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGB_est_kg_04,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$AGB_est_kg_07
	)
Tree_BENI_56789 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$PID
	))
summary(lme_BENI_AGB_RGR_20156789 <- lme(RGR_BENI_AGB_56789 ~ 
	Treatment_BENI_56789 + AGB_BENI_56789,
	random=~1 | Tree_BENI_56789))
anova(lme_BENI_AGB_RGR_20156789)
emmeans(lme_BENI_AGB_RGR_20156789, list(pairwise ~ Treatment_BENI_56789), adjust = "tukey")
# a b b b

##### AGR

### 2015-2016

summary(lm_BENI_AGB_AGR_2015_2016 <- lm(log(AGR_AGB_2015_2016) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_BENI_AGB_AGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2015-2017

# Cumulative AGR
summary(lm_BENI_AGB_AGR_2015_2017 <- lm(log(AGR_AGB_2015_2017) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_BENI_AGB_AGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a ab b ab

# Annual AGR
AGR_BENI_AGB_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_AGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2016_2017
	)
summary(lme_BENI_AGB_AGR_201567 <- lme(log(AGR_BENI_AGB_567) ~ 
	Treatment_BENI_567 + AGB_BENI_567,
	random=~1 | Tree_BENI_567))
anova(lme_BENI_AGB_AGR_201567)
emmeans(lme_BENI_AGB_AGR_201567, list(pairwise ~ Treatment_BENI_567), adjust = "tukey")
# a b b ab

### 2015-2018

# Cumulative AGR
summary(lm_BENI_AGB_AGR_2015_2018 <- lm(log(AGR_AGB_2015_2018) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_BENI_AGB_AGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a ab b ab

# Annual AGR
AGR_BENI_AGB_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_AGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_AGB_2017_2018
	)
summary(lme_BENI_AGB_AGR_2015678 <- lme(log(AGR_BENI_AGB_5678) ~ 
	Treatment_BENI_5678 + AGB_BENI_5678,
	random=~1 | Tree_BENI_5678))
anova(lme_BENI_AGB_AGR_2015678)
emmeans(lme_BENI_AGB_AGR_2015678, list(pairwise ~ Treatment_BENI_5678), adjust = "tukey")
# a b b ab

### 2015-2019

# Cumulative AGR
summary(lm_BENI_AGB_AGR_2015_2019 <- lm(log(AGR_AGB_2015_2019) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_BENI_AGB_AGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_BENI_AGB_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_AGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_AGB_2017_2018,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$AGR_AGB_2018_2019
	)
summary(lme_BENI_AGB_AGR_20156789 <- lme(log(AGR_BENI_AGB_56789) ~ 
	Treatment_BENI_56789 + AGB_BENI_56789,
	random=~1 | Tree_BENI_56789))
anova(lme_BENI_AGB_AGR_20156789)
emmeans(lme_BENI_AGB_AGR_20156789, list(pairwise ~ Treatment_BENI_56789), adjust = "tukey")
# a b b ab



### 
# BENI Belowground Biomass
###

##### RGR

### 2015-2016

summary(lm_BENI_BGB_RGR_2015_2016 <- lm(RGR_BGB_2015_2016 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_BENI_BGB_RGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

### 2015-2017

# Cumulative RGR
summary(lm_BENI_BGB_RGR_2015_2017 <- lm(RGR_BGB_2015_2017 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_BENI_BGB_RGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

# Annual RGR
RGR_BENI_BGB_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_BGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2016_2017
	)
Treatment_BENI_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
BGB_BENI_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_02
	)
Tree_BENI_567 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_BENI_BGB_RGR_201567 <- lme(RGR_BENI_BGB_567 ~ 
	Treatment_BENI_567 + BGB_BENI_567,
	random=~1 | Tree_BENI_567))
anova(lme_BENI_BGB_RGR_201567)
emmeans(lme_BENI_BGB_RGR_201567, list(pairwise ~ Treatment_BENI_567), adjust = "tukey")
# a b b ab

### 2015-2018

# Cumulative RGR
summary(lm_BENI_BGB_RGR_2015_2018 <- lm(RGR_BGB_2015_2018 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_BENI_BGB_RGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

# Annual RGR
RGR_BENI_BGB_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_BGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_BGB_2017_2018
	)
Treatment_BENI_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment
	)
BGB_BENI_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_02,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$BGB_est_kg_04
	)
Tree_BENI_5678 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID
	))
summary(lme_BENI_BGB_RGR_2015678 <- lme(RGR_BENI_BGB_5678 ~ 
	Treatment_BENI_5678 + BGB_BENI_5678,
	random=~1 | Tree_BENI_5678))
anova(lme_BENI_BGB_RGR_2015678)
emmeans(lme_BENI_BGB_RGR_2015678, list(pairwise ~ Treatment_BENI_5678), adjust = "tukey")
# a b b ab

### 2015-2019

# Cumulative RGR
summary(lm_BENI_BGB_RGR_2015_2019 <- lm(RGR_BGB_2015_2019 ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_BENI_BGB_RGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_BENI_BGB_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_BGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_BGB_2017_2018,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$RGR_BGB_2018_2019
	)
Treatment_BENI_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$Treatment
	)
BGB_BENI_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_02,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$BGB_est_kg_04,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$BGB_est_kg_07
	)
Tree_BENI_56789 <- 
	as.factor(c(
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID,
	dat[dat$Species=="BENI" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$PID
	))
summary(lme_BENI_BGB_RGR_20156789 <- lme(RGR_BENI_BGB_56789 ~ 
	Treatment_BENI_56789 + BGB_BENI_56789,
	random=~1 | Tree_BENI_56789))
anova(lme_BENI_BGB_RGR_20156789)
emmeans(lme_BENI_BGB_RGR_20156789, list(pairwise ~ Treatment_BENI_56789), adjust = "tukey")
# a b b ab

##### AGR

### 2015-2016

summary(lm_BENI_BGB_AGR_2015_2016 <- lm(log(AGR_BGB_2015_2016) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_BENI_BGB_AGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a b ab ab

### 2015-2017

# Cumulative AGR
summary(lm_BENI_BGB_AGR_2015_2017 <- lm(log(AGR_BGB_2015_2017) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_BENI_BGB_AGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a ab b ab

# Annual AGR
AGR_BENI_BGB_567 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_BGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2016_2017
	)
summary(lme_BENI_BGB_AGR_201567 <- lme(log(AGR_BENI_BGB_567) ~ 
	Treatment_BENI_567 + BGB_BENI_567,
	random=~1 | Tree_BENI_567))
anova(lme_BENI_BGB_AGR_201567)
emmeans(lme_BENI_BGB_AGR_201567, list(pairwise ~ Treatment_BENI_567), adjust = "tukey")
# a b b ab

### 2015-2018

# Cumulative AGR
summary(lm_BENI_BGB_AGR_2015_2018 <- lm(log(AGR_BGB_2015_2018) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_BENI_BGB_AGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a b b ab

# Annual AGR
AGR_BENI_BGB_5678 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_BGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_BGB_2017_2018
	)
summary(lme_BENI_BGB_AGR_2015678 <- lme(log(AGR_BENI_BGB_5678) ~ 
	Treatment_BENI_5678 + BGB_BENI_5678,
	random=~1 | Tree_BENI_5678))
anova(lme_BENI_BGB_AGR_2015678)
emmeans(lme_BENI_BGB_AGR_2015678, list(pairwise ~ Treatment_BENI_5678), adjust = "tukey")
# a b b b

### 2015-2019

# Cumulative AGR
summary(lm_BENI_BGB_AGR_2015_2019 <- lm(log(AGR_BGB_2015_2019) ~ 
	Treatment,
	data=dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_BENI_BGB_AGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_BENI_BGB_56789 <- 
	c(
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_BGB_2015_2016,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_BGB_2017_2018,
	dat[dat$Species=="BENI" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$AGR_BGB_2018_2019
	)
summary(lme_BENI_BGB_AGR_20156789 <- lme(log(AGR_BENI_BGB_56789) ~ 
	Treatment_BENI_56789 + BGB_BENI_56789,
	random=~1 | Tree_BENI_56789))
anova(lme_BENI_BGB_AGR_20156789)
emmeans(lme_BENI_BGB_AGR_20156789, list(pairwise ~ Treatment_BENI_56789), adjust = "tukey")
# a b b b





########
# ROPS #
########

### 
# ROPS Total Biomass
###

##### RGR

### 2015-2016

summary(lm_ROPS_Biomass_RGR_2015_2016 <- lm(RGR_Biomass_2015_2016 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_ROPS_Biomass_RGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# b b ab a

### 2015-2017

# Cumulative RGR
summary(lm_ROPS_Biomass_RGR_2015_2017 <- lm(RGR_Biomass_2015_2017 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ROPS_Biomass_RGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# c bc ab a

# Annual RGR
RGR_ROPS_Biomass_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_Biomass_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2016_2017
	)
Treatment_ROPS_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
Biomass_ROPS_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_02
	)
Tree_ROPS_567 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ROPS_Biomass_RGR_201567 <- lme(RGR_ROPS_Biomass_567 ~ 
	Treatment_ROPS_567 + Biomass_ROPS_567,
	random=~1 | Tree_ROPS_567))
anova(lme_ROPS_Biomass_RGR_201567)
emmeans(lme_ROPS_Biomass_RGR_201567, list(pairwise ~ Treatment_ROPS_567), adjust = "tukey")
# a a a a

### 2015-2018

# Cumulative RGR
summary(lm_ROPS_Biomass_RGR_2015_2018 <- lm(RGR_Biomass_2015_2018 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_ROPS_Biomass_RGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# c b ab a

# Annual RGR
RGR_ROPS_Biomass_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_Biomass_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_Biomass_2017_2018
	)
Treatment_ROPS_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment
	)
Biomass_ROPS_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_02,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$Biomass_est_kg_04
	)
Tree_ROPS_5678 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID
	))
summary(lme_ROPS_Biomass_RGR_2015678 <- lme(RGR_ROPS_Biomass_5678 ~ 
	Treatment_ROPS_5678 + Biomass_ROPS_5678,
	random=~1 | Tree_ROPS_5678))
anova(lme_ROPS_Biomass_RGR_2015678)
emmeans(lme_ROPS_Biomass_RGR_2015678, list(pairwise ~ Treatment_ROPS_5678), adjust = "tukey")
# a a a a

### 2015-2019

# Cumulative RGR
summary(lm_ROPS_Biomass_RGR_2015_2019 <- lm(RGR_Biomass_2015_2019 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_ROPS_Biomass_RGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# c b b a

# Annual RGR
RGR_ROPS_Biomass_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_Biomass_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_Biomass_2017_2018,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$RGR_Biomass_2018_2019
	)
Treatment_ROPS_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$Treatment
	)
Biomass_ROPS_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_02,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$Biomass_est_kg_04,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$Biomass_est_kg_07
	)
Tree_ROPS_56789 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$PID
	))
summary(lme_ROPS_Biomass_RGR_20156789 <- lme(RGR_ROPS_Biomass_56789 ~ 
	Treatment_ROPS_56789 + Biomass_ROPS_56789,
	random=~1 | Tree_ROPS_56789))
anova(lme_ROPS_Biomass_RGR_20156789)
emmeans(lme_ROPS_Biomass_RGR_20156789, list(pairwise ~ Treatment_ROPS_56789), adjust = "tukey")
# a a a a

##### AGR

### 2015-2016

summary(lm_ROPS_Biomass_AGR_2015_2016 <- lm(log(AGR_Biomass_2015_2016) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_ROPS_Biomass_AGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a b b b

### 2015-2017

# Cumulative AGR
summary(lm_ROPS_Biomass_AGR_2015_2017 <- lm(log(AGR_Biomass_2015_2017) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ROPS_Biomass_AGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a ab ab b

# Annual AGR
AGR_ROPS_Biomass_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_Biomass_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2016_2017
	)
summary(lme_ROPS_Biomass_AGR_201567 <- lme(log(AGR_ROPS_Biomass_567) ~ 
	Treatment_ROPS_567 + Biomass_ROPS_567,
	random=~1 | Tree_ROPS_567))
anova(lme_ROPS_Biomass_AGR_201567)
emmeans(lme_ROPS_Biomass_AGR_201567, list(pairwise ~ Treatment_ROPS_567), adjust = "tukey")
# a a a a

### 2015-2018

# Cumulative AGR
summary(lm_ROPS_Biomass_AGR_2015_2018 <- lm(log(AGR_Biomass_2015_2018) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_ROPS_Biomass_AGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ROPS_Biomass_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_Biomass_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_Biomass_2017_2018
	)
summary(lme_ROPS_Biomass_AGR_2015678 <- lme(log(AGR_ROPS_Biomass_5678) ~ 
	Treatment_ROPS_5678 + Biomass_ROPS_5678,
	random=~1 | Tree_ROPS_5678))
anova(lme_ROPS_Biomass_AGR_2015678)
emmeans(lme_ROPS_Biomass_AGR_2015678, list(pairwise ~ Treatment_ROPS_5678), adjust = "tukey")
# a a a a

### 2015-2019

# Cumulative AGR
summary(lm_ROPS_Biomass_AGR_2015_2019 <- lm(log(AGR_Biomass_2015_2019) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_ROPS_Biomass_AGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ROPS_Biomass_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_Biomass_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_Biomass_2017_2018,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_08 == 1,]$AGR_Biomass_2018_2019
	)
summary(lme_ROPS_Biomass_AGR_20156789 <- lme(log(AGR_ROPS_Biomass_56789) ~ 
	Treatment_ROPS_56789 + Biomass_ROPS_56789,
	random=~1 | Tree_ROPS_56789))
anova(lme_ROPS_Biomass_AGR_20156789)
emmeans(lme_ROPS_Biomass_AGR_20156789, list(pairwise ~ Treatment_ROPS_56789), adjust = "tukey")
# a a a a






### 
# ROPS Aboveground Biomass
###

##### RGR

### 2015-2016

summary(lm_ROPS_AGB_RGR_2015_2016 <- lm(RGR_AGB_2015_2016 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_ROPS_AGB_RGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a ab ab b

### 2015-2017

# Cumulative RGR
summary(lm_ROPS_AGB_RGR_2015_2017 <- lm(RGR_AGB_2015_2017 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ROPS_AGB_RGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a b ab b

# Annual RGR
RGR_ROPS_AGB_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_AGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2016_2017
	)
Treatment_ROPS_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
AGB_ROPS_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_02
	)
Tree_ROPS_567 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ROPS_AGB_RGR_201567 <- lme(RGR_ROPS_AGB_567 ~ 
	Treatment_ROPS_567 + AGB_ROPS_567,
	random=~1 | Tree_ROPS_567))
anova(lme_ROPS_AGB_RGR_201567)
emmeans(lme_ROPS_AGB_RGR_201567, list(pairwise ~ Treatment_ROPS_567), adjust = "tukey")
# a a a a

### 2015-2018

# Cumulative RGR
summary(lm_ROPS_AGB_RGR_2015_2018 <- lm(RGR_AGB_2015_2018 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_ROPS_AGB_RGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a b b b

# Annual RGR
RGR_ROPS_AGB_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_AGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_AGB_2017_2018
	)
Treatment_ROPS_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment
	)
AGB_ROPS_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_02,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGB_est_kg_04
	)
Tree_ROPS_5678 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID
	))
summary(lme_ROPS_AGB_RGR_2015678 <- lme(RGR_ROPS_AGB_5678 ~ 
	Treatment_ROPS_5678 + AGB_ROPS_5678,
	random=~1 | Tree_ROPS_5678))
anova(lme_ROPS_AGB_RGR_2015678)
emmeans(lme_ROPS_AGB_RGR_2015678, list(pairwise ~ Treatment_ROPS_5678), adjust = "tukey")
# a a a a

### 2015-2019

# Cumulative RGR
summary(lm_ROPS_AGB_RGR_2015_2019 <- lm(RGR_AGB_2015_2019 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_ROPS_AGB_RGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# c b b a

# Annual RGR
RGR_ROPS_AGB_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_AGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_AGB_2017_2018,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$RGR_AGB_2018_2019
	)
Treatment_ROPS_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$Treatment
	)
AGB_ROPS_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_02,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGB_est_kg_04,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$AGB_est_kg_07
	)
Tree_ROPS_56789 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$PID
	))
summary(lme_ROPS_AGB_RGR_20156789 <- lme(RGR_ROPS_AGB_56789 ~ 
	Treatment_ROPS_56789 + AGB_ROPS_56789,
	random=~1 | Tree_ROPS_56789))
anova(lme_ROPS_AGB_RGR_20156789)
emmeans(lme_ROPS_AGB_RGR_20156789, list(pairwise ~ Treatment_ROPS_56789), adjust = "tukey")
# a a a a

##### AGR

### 2015-2016

summary(lm_ROPS_AGB_AGR_2015_2016 <- lm(log(AGR_AGB_2015_2016) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_ROPS_AGB_AGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a b b b

### 2015-2017

# Cumulative AGR
summary(lm_ROPS_AGB_AGR_2015_2017 <- lm(log(AGR_AGB_2015_2017) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ROPS_AGB_AGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ROPS_AGB_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_AGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2016_2017
	)
summary(lme_ROPS_AGB_AGR_201567 <- lme(log(AGR_ROPS_AGB_567) ~ 
	Treatment_ROPS_567 + AGB_ROPS_567,
	random=~1 | Tree_ROPS_567))
anova(lme_ROPS_AGB_AGR_201567)
emmeans(lme_ROPS_AGB_AGR_201567, list(pairwise ~ Treatment_ROPS_567), adjust = "tukey")
# a a a a

### 2015-2018

# Cumulative AGR
summary(lm_ROPS_AGB_AGR_2015_2018 <- lm(log(AGR_AGB_2015_2018) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_ROPS_AGB_AGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ROPS_AGB_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_AGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_AGB_2017_2018
	)
summary(lme_ROPS_AGB_AGR_2015678 <- lme(log(AGR_ROPS_AGB_5678) ~ 
	Treatment_ROPS_5678 + AGB_ROPS_5678,
	random=~1 | Tree_ROPS_5678))
anova(lme_ROPS_AGB_AGR_2015678)
emmeans(lme_ROPS_AGB_AGR_2015678, list(pairwise ~ Treatment_ROPS_5678), adjust = "tukey")
# a a a a

### 2015-2019

# Cumulative AGR
summary(lm_ROPS_AGB_AGR_2015_2019 <- lm(log(AGR_AGB_2015_2019) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_ROPS_AGB_AGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ROPS_AGB_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_AGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_AGB_2017_2018,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$AGR_AGB_2018_2019
	)
summary(lme_ROPS_AGB_AGR_20156789 <- lme(log(AGR_ROPS_AGB_56789) ~ 
	Treatment_ROPS_56789 + AGB_ROPS_56789,
	random=~1 | Tree_ROPS_56789))
anova(lme_ROPS_AGB_AGR_20156789)
emmeans(lme_ROPS_AGB_AGR_20156789, list(pairwise ~ Treatment_ROPS_56789), adjust = "tukey")
# a a a a



### 
# ROPS Belowground Biomass
###

##### RGR

### 2015-2016

summary(lm_ROPS_BGB_RGR_2015_2016 <- lm(RGR_BGB_2015_2016 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_ROPS_BGB_RGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# ab b ab a

### 2015-2017

# Cumulative RGR
summary(lm_ROPS_BGB_RGR_2015_2017 <- lm(RGR_BGB_2015_2017 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ROPS_BGB_RGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# ab b ab a

# Annual RGR
RGR_ROPS_BGB_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_BGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2016_2017
	)
Treatment_ROPS_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
BGB_ROPS_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_02
	)
Tree_ROPS_567 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ROPS_BGB_RGR_201567 <- lme(RGR_ROPS_BGB_567 ~ 
	Treatment_ROPS_567 + BGB_ROPS_567,
	random=~1 | Tree_ROPS_567))
anova(lme_ROPS_BGB_RGR_201567)
emmeans(lme_ROPS_BGB_RGR_201567, list(pairwise ~ Treatment_ROPS_567), adjust = "tukey")
# a a a a

### 2015-2018

# Cumulative RGR
summary(lm_ROPS_BGB_RGR_2015_2018 <- lm(RGR_BGB_2015_2018 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_ROPS_BGB_RGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# ab b ab a

# Annual RGR
RGR_ROPS_BGB_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_BGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_BGB_2017_2018
	)
Treatment_ROPS_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment
	)
BGB_ROPS_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_02,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$BGB_est_kg_04
	)
Tree_ROPS_5678 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID
	))
summary(lme_ROPS_BGB_RGR_2015678 <- lme(RGR_ROPS_BGB_5678 ~ 
	Treatment_ROPS_5678 + BGB_ROPS_5678,
	random=~1 | Tree_ROPS_5678))
anova(lme_ROPS_BGB_RGR_2015678)
emmeans(lme_ROPS_BGB_RGR_2015678, list(pairwise ~ Treatment_ROPS_5678), adjust = "tukey")
# a a a a

### 2015-2019

# Cumulative RGR
summary(lm_ROPS_BGB_RGR_2015_2019 <- lm(RGR_BGB_2015_2019 ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_ROPS_BGB_RGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# b c ab a

# Annual RGR
RGR_ROPS_BGB_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$RGR_BGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$RGR_BGB_2017_2018,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$RGR_BGB_2018_2019
	)
Treatment_ROPS_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$Treatment,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$Treatment
	)
BGB_ROPS_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_02,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$BGB_est_kg_04,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$BGB_est_kg_07
	)
Tree_ROPS_56789 <- 
	as.factor(c(
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$PID,
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$PID
	))
summary(lme_ROPS_BGB_RGR_20156789 <- lme(RGR_ROPS_BGB_56789 ~ 
	Treatment_ROPS_56789 + BGB_ROPS_56789,
	random=~1 | Tree_ROPS_56789))
anova(lme_ROPS_BGB_RGR_20156789)
emmeans(lme_ROPS_BGB_RGR_20156789, list(pairwise ~ Treatment_ROPS_56789), adjust = "tukey")
# a a a a

##### AGR

### 2015-2016

summary(lm_ROPS_BGB_AGR_2015_2016 <- lm(log(AGR_BGB_2015_2016) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2015_2016) & 
	dat$Use_growth_02 == 1,]))
emmeans(lm_ROPS_BGB_AGR_2015_2016, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2015-2017

# Cumulative AGR
summary(lm_ROPS_BGB_AGR_2015_2017 <- lm(log(AGR_BGB_2015_2017) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2015_2017) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ROPS_BGB_AGR_2015_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ROPS_BGB_567 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_BGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2016_2017
	)
summary(lme_ROPS_BGB_AGR_201567 <- lme(log(AGR_ROPS_BGB_567) ~ 
	Treatment_ROPS_567 + BGB_ROPS_567,
	random=~1 | Tree_ROPS_567))
anova(lme_ROPS_BGB_AGR_201567)
emmeans(lme_ROPS_BGB_AGR_201567, list(pairwise ~ Treatment_ROPS_567), adjust = "tukey")
# a a a a

### 2015-2018

# Cumulative AGR
summary(lm_ROPS_BGB_AGR_2015_2018 <- lm(log(AGR_BGB_2015_2018) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2015_2018) & 
	dat$Use_growth_07 == 1,]))
emmeans(lm_ROPS_BGB_AGR_2015_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ROPS_BGB_5678 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_BGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_BGB_2017_2018
	)
summary(lme_ROPS_BGB_AGR_2015678 <- lme(log(AGR_ROPS_BGB_5678) ~ 
	Treatment_ROPS_5678 + BGB_ROPS_5678,
	random=~1 | Tree_ROPS_5678))
anova(lme_ROPS_BGB_AGR_2015678)
emmeans(lme_ROPS_BGB_AGR_2015678, list(pairwise ~ Treatment_ROPS_5678), adjust = "tukey")
# a a a a

### 2015-2019

# Cumulative AGR
summary(lm_ROPS_BGB_AGR_2015_2019 <- lm(log(AGR_BGB_2015_2019) ~ 
	Treatment,
	data=dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2015_2019) & 
	dat$Use_growth_08 == 1,]))
emmeans(lm_ROPS_BGB_AGR_2015_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ROPS_BGB_56789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2015_2016) & 
		dat$Use_growth_02 == 1,]$AGR_BGB_2015_2016,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_07 == 1,]$AGR_BGB_2017_2018,
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_08 == 1,]$AGR_BGB_2018_2019
	)
summary(lme_ROPS_BGB_AGR_20156789 <- lme(log(AGR_ROPS_BGB_56789) ~ 
	Treatment_ROPS_56789 + BGB_ROPS_56789,
	random=~1 | Tree_ROPS_56789))
anova(lme_ROPS_BGB_AGR_20156789)
emmeans(lme_ROPS_BGB_AGR_20156789, list(pairwise ~ Treatment_ROPS_56789), adjust = "tukey")
# a a a a

#############################################################################
#############################################################################
#############################################################################