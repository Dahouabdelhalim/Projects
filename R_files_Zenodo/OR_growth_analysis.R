# NSF FX growth analysis for the Oregon biomass data

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
# Also want to do these for each year: Up to 2017, 2018, 2019, 2020.

# The "Use_growth_0X" columns in "OR_FX_Size_FoliarCNIsotope_Data.csv" 
# are the codes for whether or not to use data.
# 1 means use.
# 0 means dead or sufficiently ill or damaged that it shouldn't be used.
# Only one tree (a MN ALRU) in OR is odd.

rm(list=ls())

library(lme4)
library(nlme)
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

dat <- read.csv("OR_FX_Size_FoliarCNIsotope_Data.csv")[1:64,1:122]

# Option to print out how many are in each category
# (uncomment the following lines)

#print(dat[dat$Species=="PSME" & dat$Treatment=="LN",c(43,54,71,87)])
#print(dat[dat$Species=="PSME" & dat$Treatment=="MN",c(43,54,71,87)])
#print(dat[dat$Species=="PSME" & dat$Treatment=="HN",c(43,54,71,87)])
#print(dat[dat$Species=="PSME" & dat$Treatment=="PHN",c(43,54,71,87)])

#print(dat[dat$Species=="ALRU" & dat$Treatment=="LN",c(43,54,71,87)])
#print(dat[dat$Species=="ALRU" & dat$Treatment=="MN",c(43,54,71,87)])
#print(dat[dat$Species=="ALRU" & dat$Treatment=="HN",c(43,54,71,87)])
#print(dat[dat$Species=="ALRU" & dat$Treatment=="PHN",c(43,54,71,87)])

#############################################################
################ Calculate growth rates #####################
#############################################################

# "AGR" is absolute growth rate and "RGR" is relative growth rate. 
# Each is calculated for aboveground, belowground, and total biomass, 
# for the following time periods: 

# Annual increments:

# Biomass

# June 2016 to October 2017
dat$AGR_Biomass_2016_2017 <- (dat$Biomass_est_kg_03 - dat$Biomass_est_kg_01)/((dat$Days_03 - dat$Days_01)/365.25)
dat$RGR_Biomass_2016_2017 <- (log(dat$Biomass_est_kg_03) - log(dat$Biomass_est_kg_01))/((dat$Days_03 - dat$Days_01)/365.25)

# October 2017 to October 2018
dat$AGR_Biomass_2017_2018 <- (dat$Biomass_est_kg_04 - dat$Biomass_est_kg_03)/((dat$Days_04 - dat$Days_03)/365.25)
dat$RGR_Biomass_2017_2018 <- (log(dat$Biomass_est_kg_04) - log(dat$Biomass_est_kg_03))/((dat$Days_04 - dat$Days_03)/365.25)

# October 2018 to September 2019
dat$AGR_Biomass_2018_2019 <- (dat$Biomass_est_kg_05 - dat$Biomass_est_kg_04)/((dat$Days_05 - dat$Days_04)/365.25)
dat$RGR_Biomass_2018_2019 <- (log(dat$Biomass_est_kg_05) - log(dat$Biomass_est_kg_04))/((dat$Days_05 - dat$Days_04)/365.25)

# September 2019 to September 2020
dat$AGR_Biomass_2019_2020 <- (dat$Biomass_est_kg_06 - dat$Biomass_est_kg_05)/((dat$Days_06 - dat$Days_05)/365.25)
dat$RGR_Biomass_2019_2020 <- (log(dat$Biomass_est_kg_06) - log(dat$Biomass_est_kg_05))/((dat$Days_06 - dat$Days_05)/365.25)

# AGB

# June 2016 to October 2017
dat$AGR_AGB_2016_2017 <- (dat$AGB_est_kg_03 - dat$AGB_est_kg_01)/((dat$Days_03 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2017 <- (log(dat$AGB_est_kg_03) - log(dat$AGB_est_kg_01))/((dat$Days_03 - dat$Days_01)/365.25)

# October 2017 to October 2018
dat$AGR_AGB_2017_2018 <- (dat$AGB_est_kg_04 - dat$AGB_est_kg_03)/((dat$Days_04 - dat$Days_03)/365.25)
dat$RGR_AGB_2017_2018 <- (log(dat$AGB_est_kg_04) - log(dat$AGB_est_kg_03))/((dat$Days_04 - dat$Days_03)/365.25)

# October 2018 to September 2019
dat$AGR_AGB_2018_2019 <- (dat$AGB_est_kg_05 - dat$AGB_est_kg_04)/((dat$Days_05 - dat$Days_04)/365.25)
dat$RGR_AGB_2018_2019 <- (log(dat$AGB_est_kg_05) - log(dat$AGB_est_kg_04))/((dat$Days_05 - dat$Days_04)/365.25)

# September 2019 to September 2020
dat$AGR_AGB_2019_2020 <- (dat$AGB_est_kg_06 - dat$AGB_est_kg_05)/((dat$Days_06 - dat$Days_05)/365.25)
dat$RGR_AGB_2019_2020 <- (log(dat$AGB_est_kg_06) - log(dat$AGB_est_kg_05))/((dat$Days_06 - dat$Days_05)/365.25)

# BGB

# June 2016 to October 2017
dat$AGR_BGB_2016_2017 <- (dat$BGB_est_kg_03 - dat$BGB_est_kg_01)/((dat$Days_03 - dat$Days_01)/365.25)
dat$RGR_BGB_2016_2017 <- (log(dat$BGB_est_kg_03) - log(dat$BGB_est_kg_01))/((dat$Days_03 - dat$Days_01)/365.25)

# October 2017 to October 2018
dat$AGR_BGB_2017_2018 <- (dat$BGB_est_kg_04 - dat$BGB_est_kg_03)/((dat$Days_04 - dat$Days_03)/365.25)
dat$RGR_BGB_2017_2018 <- (log(dat$BGB_est_kg_04) - log(dat$BGB_est_kg_03))/((dat$Days_04 - dat$Days_03)/365.25)

# October 2018 to September 2019
dat$AGR_BGB_2018_2019 <- (dat$BGB_est_kg_05 - dat$BGB_est_kg_04)/((dat$Days_05 - dat$Days_04)/365.25)
dat$RGR_BGB_2018_2019 <- (log(dat$BGB_est_kg_05) - log(dat$BGB_est_kg_04))/((dat$Days_05 - dat$Days_04)/365.25)

# September 2019 to September 2020
dat$AGR_BGB_2019_2020 <- (dat$BGB_est_kg_06 - dat$BGB_est_kg_05)/((dat$Days_06 - dat$Days_05)/365.25)
dat$RGR_BGB_2019_2020 <- (log(dat$BGB_est_kg_06) - log(dat$BGB_est_kg_05))/((dat$Days_06 - dat$Days_05)/365.25)



# 2016 to tf increments:

# Biomass

# June 2016 to October 2018
dat$AGR_Biomass_2016_2018 <- (dat$Biomass_est_kg_04 - dat$Biomass_est_kg_01)/((dat$Days_04 - dat$Days_01)/365.25)
dat$RGR_Biomass_2016_2018 <- (log(dat$Biomass_est_kg_04) - log(dat$Biomass_est_kg_01))/((dat$Days_04 - dat$Days_01)/365.25)

# June 2016 to September 2019
dat$AGR_Biomass_2016_2019 <- (dat$Biomass_est_kg_05 - dat$Biomass_est_kg_01)/((dat$Days_05 - dat$Days_01)/365.25)
dat$RGR_Biomass_2016_2019 <- (log(dat$Biomass_est_kg_05) - log(dat$Biomass_est_kg_01))/((dat$Days_05 - dat$Days_01)/365.25)

# June 2016 to September 2020
dat$AGR_Biomass_2016_2020 <- (dat$Biomass_est_kg_06 - dat$Biomass_est_kg_01)/((dat$Days_06 - dat$Days_01)/365.25)
dat$RGR_Biomass_2016_2020 <- (log(dat$Biomass_est_kg_06) - log(dat$Biomass_est_kg_01))/((dat$Days_06 - dat$Days_01)/365.25)

# AGB

# June 2016 to October 2018
dat$AGR_AGB_2016_2018 <- (dat$AGB_est_kg_04 - dat$AGB_est_kg_01)/((dat$Days_04 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2018 <- (log(dat$AGB_est_kg_04) - log(dat$AGB_est_kg_01))/((dat$Days_04 - dat$Days_01)/365.25)

# June 2016 to September 2019
dat$AGR_AGB_2016_2019 <- (dat$AGB_est_kg_05 - dat$AGB_est_kg_01)/((dat$Days_05 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2019 <- (log(dat$AGB_est_kg_05) - log(dat$AGB_est_kg_01))/((dat$Days_05 - dat$Days_01)/365.25)

# June 2016 to September 2020
dat$AGR_AGB_2016_2020 <- (dat$AGB_est_kg_06 - dat$AGB_est_kg_01)/((dat$Days_06 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2020 <- (log(dat$AGB_est_kg_06) - log(dat$AGB_est_kg_01))/((dat$Days_06 - dat$Days_01)/365.25)

# BGB

# June 2016 to October 2018
dat$AGR_BGB_2016_2018 <- (dat$BGB_est_kg_04 - dat$BGB_est_kg_01)/((dat$Days_04 - dat$Days_01)/365.25)
dat$RGR_BGB_2016_2018 <- (log(dat$BGB_est_kg_04) - log(dat$BGB_est_kg_01))/((dat$Days_04 - dat$Days_01)/365.25)

# June 2016 to September 2019
dat$AGR_BGB_2016_2019 <- (dat$BGB_est_kg_05 - dat$BGB_est_kg_01)/((dat$Days_05 - dat$Days_01)/365.25)
dat$RGR_BGB_2016_2019 <- (log(dat$BGB_est_kg_05) - log(dat$BGB_est_kg_01))/((dat$Days_05 - dat$Days_01)/365.25)

# June 2016 to September 2020
dat$AGR_BGB_2016_2020 <- (dat$BGB_est_kg_06 - dat$BGB_est_kg_01)/((dat$Days_06 - dat$Days_01)/365.25)
dat$RGR_BGB_2016_2020 <- (log(dat$BGB_est_kg_06) - log(dat$BGB_est_kg_01))/((dat$Days_06 - dat$Days_01)/365.25)


###############################################################
########################### Stats #############################
###############################################################

# Only use trees/measurements with Use_growth_0X == 1
# Response variables:
# PSME and ALRU ...
# Total biomass, AGB, and BGB for ...
# 2017: 
# - AGR 2016-2017
# - RGR 2016-2017
# 2018: 
# - AGR 2016-2017, 2017-2018
# - RGR 2016-2017, 2017-2018
# - AGR 2016-2018
# - RGR 2016-2018
# 2019: 
# - AGR 2016-2017, 2017-2018, 2018-2019
# - RGR 2016-2017, 2017-2018, 2018-2019
# - AGR 2016-2019
# - RGR 2016-2019
# 2020: 
# - AGR 2016-2017, 2017-2018, 2018-2019, 2019-2020
# - RGR 2016-2017, 2017-2018, 2018-2019, 2019-2020
# - AGR 2016-2020
# - RGR 2016-2020
#
# Driver: Treatment
# Covariates: For multi-year, biomass/AGB/BGB (fixed), tree (random)

########
# PSME #
########

### 
# PSME Total Biomass
###

##### RGR

### 2016-2017

summary(lm_PSME_Biomass_RGR_2016_2017 <- lm(RGR_Biomass_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_PSME_Biomass_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_PSME_Biomass_RGR_2016_2018 <- lm(RGR_Biomass_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_PSME_Biomass_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_Biomass_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2017_2018
	)
Treatment_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
Biomass_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03
	)
Tree_PSME_678 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_PSME_Biomass_RGR_201678 <- lme(RGR_PSME_Biomass_678 ~ 
	Treatment_PSME_678 + Biomass_PSME_678,
	random=~1 | Tree_PSME_678))
anova(lme_PSME_Biomass_RGR_201678)
emmeans(lme_PSME_Biomass_RGR_201678, list(pairwise ~ Treatment_PSME_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_PSME_Biomass_RGR_2016_2019 <- lm(RGR_Biomass_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_PSME_Biomass_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_Biomass_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_Biomass_2018_2019
	)
Treatment_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
Biomass_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Biomass_est_kg_04
	)
Tree_PSME_6789 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_PSME_Biomass_RGR_2016789 <- lme(RGR_PSME_Biomass_6789 ~ 
	Treatment_PSME_6789 + Biomass_PSME_6789,
	random=~1 | Tree_PSME_6789))
anova(lme_PSME_Biomass_RGR_2016789)
emmeans(lme_PSME_Biomass_RGR_2016789, list(pairwise ~ Treatment_PSME_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative RGR
summary(lm_PSME_Biomass_RGR_2016_2020 <- lm(RGR_Biomass_2016_2020 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_PSME_Biomass_RGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_Biomass_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_Biomass_2018_2019,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$RGR_Biomass_2019_2020
	)
Treatment_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Biomass_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Biomass_est_kg_04,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$Biomass_est_kg_05
	)
Tree_PSME_67890 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_PSME_Biomass_RGR_20167890 <- lme(RGR_PSME_Biomass_67890 ~ 
	Treatment_PSME_67890 + Biomass_PSME_67890,
	random=~1 | Tree_PSME_67890))
anova(lme_PSME_Biomass_RGR_20167890)
emmeans(lme_PSME_Biomass_RGR_20167890, list(pairwise ~ Treatment_PSME_67890), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_PSME_Biomass_AGR_2016_2017 <- lm(log(AGR_Biomass_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_PSME_Biomass_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_PSME_Biomass_AGR_2016_2018 <- lm(log(AGR_Biomass_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_PSME_Biomass_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_Biomass_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2017_2018
	)
Treatment_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
Biomass_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03
	)
Tree_PSME_678 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_PSME_Biomass_AGR_201678 <- lme(log(AGR_PSME_Biomass_678) ~ 
	Treatment_PSME_678 + Biomass_PSME_678,
	random=~1 | Tree_PSME_678))
anova(lme_PSME_Biomass_AGR_201678)
emmeans(lme_PSME_Biomass_AGR_201678, list(pairwise ~ Treatment_PSME_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_PSME_Biomass_AGR_2016_2019 <- lm(log(AGR_Biomass_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_PSME_Biomass_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_Biomass_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_Biomass_2018_2019
	)
Treatment_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
Biomass_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Biomass_est_kg_04
	)
Tree_PSME_6789 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_PSME_Biomass_AGR_2016789 <- lme(log(AGR_PSME_Biomass_6789) ~ 
	Treatment_PSME_6789 + Biomass_PSME_6789,
	random=~1 | Tree_PSME_6789))
anova(lme_PSME_Biomass_AGR_2016789)
emmeans(lme_PSME_Biomass_AGR_2016789, list(pairwise ~ Treatment_PSME_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative AGR
summary(lm_PSME_Biomass_AGR_2016_2020 <- lm(log(AGR_Biomass_2016_2020) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_PSME_Biomass_AGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_Biomass_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_Biomass_2018_2019,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGR_Biomass_2019_2020
	)
Treatment_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Biomass_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Biomass_est_kg_04,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$Biomass_est_kg_05
	)
Tree_PSME_67890 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_PSME_Biomass_AGR_20167890 <- lme(log(AGR_PSME_Biomass_67890) ~ 
	Treatment_PSME_67890 + Biomass_PSME_67890,
	random=~1 | Tree_PSME_67890))
anova(lme_PSME_Biomass_AGR_20167890)
emmeans(lme_PSME_Biomass_AGR_20167890, list(pairwise ~ Treatment_PSME_67890), adjust = "tukey")
# a a a a





### 
# PSME Aboveground Biomass
###

##### RGR

### 2016-2017

summary(lm_PSME_AGB_RGR_2016_2017 <- lm(RGR_AGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_PSME_AGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_PSME_AGB_RGR_2016_2018 <- lm(RGR_AGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_PSME_AGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_AGB_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2017_2018
	)
Treatment_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
AGB_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03
	)
Tree_PSME_678 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_PSME_AGB_RGR_201678 <- lme(RGR_PSME_AGB_678 ~ 
	Treatment_PSME_678 + AGB_PSME_678,
	random=~1 | Tree_PSME_678))
anova(lme_PSME_AGB_RGR_201678)
emmeans(lme_PSME_AGB_RGR_201678, list(pairwise ~ Treatment_PSME_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_PSME_AGB_RGR_2016_2019 <- lm(RGR_AGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_PSME_AGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_AGB_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_AGB_2018_2019
	)
Treatment_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
AGB_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGB_est_kg_04
	)
Tree_PSME_6789 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_PSME_AGB_RGR_2016789 <- lme(RGR_PSME_AGB_6789 ~ 
	Treatment_PSME_6789 + AGB_PSME_6789,
	random=~1 | Tree_PSME_6789))
anova(lme_PSME_AGB_RGR_2016789)
emmeans(lme_PSME_AGB_RGR_2016789, list(pairwise ~ Treatment_PSME_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative RGR
summary(lm_PSME_AGB_RGR_2016_2020 <- lm(RGR_AGB_2016_2020 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_PSME_AGB_RGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_AGB_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_AGB_2018_2019,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$RGR_AGB_2019_2020
	)
Treatment_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
AGB_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGB_est_kg_04,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGB_est_kg_05
	)
Tree_PSME_67890 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_PSME_AGB_RGR_20167890 <- lme(RGR_PSME_AGB_67890 ~ 
	Treatment_PSME_67890 + AGB_PSME_67890,
	random=~1 | Tree_PSME_67890))
anova(lme_PSME_AGB_RGR_20167890)
emmeans(lme_PSME_AGB_RGR_20167890, list(pairwise ~ Treatment_PSME_67890), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_PSME_AGB_AGR_2016_2017 <- lm(log(AGR_AGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_PSME_AGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_PSME_AGB_AGR_2016_2018 <- lm(log(AGR_AGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_PSME_AGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_AGB_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2017_2018
	)
Treatment_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
AGB_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03
	)
Tree_PSME_678 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_PSME_AGB_AGR_201678 <- lme(log(AGR_PSME_AGB_678) ~ 
	Treatment_PSME_678 + AGB_PSME_678,
	random=~1 | Tree_PSME_678))
anova(lme_PSME_AGB_AGR_201678)
emmeans(lme_PSME_AGB_AGR_201678, list(pairwise ~ Treatment_PSME_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_PSME_AGB_AGR_2016_2019 <- lm(log(AGR_AGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_PSME_AGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_AGB_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_AGB_2018_2019
	)
Treatment_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
AGB_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGB_est_kg_04
	)
Tree_PSME_6789 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_PSME_AGB_AGR_2016789 <- lme(log(AGR_PSME_AGB_6789) ~ 
	Treatment_PSME_6789 + AGB_PSME_6789,
	random=~1 | Tree_PSME_6789))
anova(lme_PSME_AGB_AGR_2016789)
emmeans(lme_PSME_AGB_AGR_2016789, list(pairwise ~ Treatment_PSME_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative AGR
summary(lm_PSME_AGB_AGR_2016_2020 <- lm(log(AGR_AGB_2016_2020) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_PSME_AGB_AGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_AGB_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_AGB_2018_2019,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGR_AGB_2019_2020
	)
Treatment_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
AGB_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGB_est_kg_04,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGB_est_kg_05
	)
Tree_PSME_67890 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_PSME_AGB_AGR_20167890 <- lme(log(AGR_PSME_AGB_67890) ~ 
	Treatment_PSME_67890 + AGB_PSME_67890,
	random=~1 | Tree_PSME_67890))
anova(lme_PSME_AGB_AGR_20167890)
emmeans(lme_PSME_AGB_AGR_20167890, list(pairwise ~ Treatment_PSME_67890), adjust = "tukey")
# a a a a




### 
# PSME Belowground Biomass
###

##### RGR

### 2016-2017

summary(lm_PSME_BGB_RGR_2016_2017 <- lm(RGR_BGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_PSME_BGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_PSME_BGB_RGR_2016_2018 <- lm(RGR_BGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_PSME_BGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_BGB_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2017_2018
	)
Treatment_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
BGB_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03
	)
Tree_PSME_678 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_PSME_BGB_RGR_201678 <- lme(RGR_PSME_BGB_678 ~ 
	Treatment_PSME_678 + BGB_PSME_678,
	random=~1 | Tree_PSME_678))
anova(lme_PSME_BGB_RGR_201678)
emmeans(lme_PSME_BGB_RGR_201678, list(pairwise ~ Treatment_PSME_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_PSME_BGB_RGR_2016_2019 <- lm(RGR_BGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_PSME_BGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_BGB_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_BGB_2018_2019
	)
Treatment_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
BGB_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$BGB_est_kg_04
	)
Tree_PSME_6789 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_PSME_BGB_RGR_2016789 <- lme(RGR_PSME_BGB_6789 ~ 
	Treatment_PSME_6789 + BGB_PSME_6789,
	random=~1 | Tree_PSME_6789))
anova(lme_PSME_BGB_RGR_2016789)
emmeans(lme_PSME_BGB_RGR_2016789, list(pairwise ~ Treatment_PSME_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative RGR
summary(lm_PSME_BGB_RGR_2016_2020 <- lm(RGR_BGB_2016_2020 ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_PSME_BGB_RGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_PSME_BGB_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_BGB_2018_2019,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$RGR_BGB_2019_2020
	)
Treatment_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
BGB_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$BGB_est_kg_04,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$BGB_est_kg_05
	)
Tree_PSME_67890 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$RGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_PSME_BGB_RGR_20167890 <- lme(RGR_PSME_BGB_67890 ~ 
	Treatment_PSME_67890 + BGB_PSME_67890,
	random=~1 | Tree_PSME_67890))
anova(lme_PSME_BGB_RGR_20167890)
emmeans(lme_PSME_BGB_RGR_20167890, list(pairwise ~ Treatment_PSME_67890), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_PSME_BGB_AGR_2016_2017 <- lm(log(AGR_BGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_PSME_BGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_PSME_BGB_AGR_2016_2018 <- lm(log(AGR_BGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_PSME_BGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_BGB_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2017_2018
	)
Treatment_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
BGB_PSME_678 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03
	)
Tree_PSME_678 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_PSME_BGB_AGR_201678 <- lme(log(AGR_PSME_BGB_678) ~ 
	Treatment_PSME_678 + BGB_PSME_678,
	random=~1 | Tree_PSME_678))
anova(lme_PSME_BGB_AGR_201678)
emmeans(lme_PSME_BGB_AGR_201678, list(pairwise ~ Treatment_PSME_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_PSME_BGB_AGR_2016_2019 <- lm(log(AGR_BGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_PSME_BGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_BGB_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_BGB_2018_2019
	)
Treatment_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
BGB_PSME_6789 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$BGB_est_kg_04
	)
Tree_PSME_6789 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_PSME_BGB_AGR_2016789 <- lme(log(AGR_PSME_BGB_6789) ~ 
	Treatment_PSME_6789 + BGB_PSME_6789,
	random=~1 | Tree_PSME_6789))
anova(lme_PSME_BGB_AGR_2016789)
emmeans(lme_PSME_BGB_AGR_2016789, list(pairwise ~ Treatment_PSME_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative AGR
summary(lm_PSME_BGB_AGR_2016_2020 <- lm(log(AGR_BGB_2016_2020) ~ 
	Treatment,
	data=dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_PSME_BGB_AGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSME_BGB_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2017_2018,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_BGB_2018_2019,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGR_BGB_2019_2020
	)
Treatment_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
BGB_PSME_67890 <- 
	c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$BGB_est_kg_04,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$BGB_est_kg_05
	)
Tree_PSME_67890 <- 
	as.factor(c(
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="PSME" & !is.na(dat$AGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_PSME_BGB_AGR_20167890 <- lme(log(AGR_PSME_BGB_67890) ~ 
	Treatment_PSME_67890 + BGB_PSME_67890,
	random=~1 | Tree_PSME_67890))
anova(lme_PSME_BGB_AGR_20167890)
emmeans(lme_PSME_BGB_AGR_20167890, list(pairwise ~ Treatment_PSME_67890), adjust = "tukey")
# a a a a





########
# ALRU #
########

### 
# ALRU Total Biomass
###

##### RGR

### 2016-2017

summary(lm_ALRU_Biomass_RGR_2016_2017 <- lm(RGR_Biomass_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_ALRU_Biomass_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_ALRU_Biomass_RGR_2016_2018 <- lm(RGR_Biomass_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ALRU_Biomass_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_Biomass_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2017_2018
	)
Treatment_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
Biomass_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03
	)
Tree_ALRU_678 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ALRU_Biomass_RGR_201678 <- lme(RGR_ALRU_Biomass_678 ~ 
	Treatment_ALRU_678 + Biomass_ALRU_678,
	random=~1 | Tree_ALRU_678))
anova(lme_ALRU_Biomass_RGR_201678)
emmeans(lme_ALRU_Biomass_RGR_201678, list(pairwise ~ Treatment_ALRU_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_ALRU_Biomass_RGR_2016_2019 <- lm(RGR_Biomass_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_ALRU_Biomass_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_Biomass_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_Biomass_2018_2019
	)
Treatment_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
Biomass_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Biomass_est_kg_04
	)
Tree_ALRU_6789 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_ALRU_Biomass_RGR_2016789 <- lme(RGR_ALRU_Biomass_6789 ~ 
	Treatment_ALRU_6789 + Biomass_ALRU_6789,
	random=~1 | Tree_ALRU_6789))
anova(lme_ALRU_Biomass_RGR_2016789)
emmeans(lme_ALRU_Biomass_RGR_2016789, list(pairwise ~ Treatment_ALRU_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative RGR
summary(lm_ALRU_Biomass_RGR_2016_2020 <- lm(RGR_Biomass_2016_2020 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_ALRU_Biomass_RGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_Biomass_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_Biomass_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_Biomass_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_Biomass_2018_2019,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$RGR_Biomass_2019_2020
	)
Treatment_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Biomass_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Biomass_est_kg_04,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$Biomass_est_kg_05
	)
Tree_ALRU_67890 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_ALRU_Biomass_RGR_20167890 <- lme(RGR_ALRU_Biomass_67890 ~ 
	Treatment_ALRU_67890 + Biomass_ALRU_67890,
	random=~1 | Tree_ALRU_67890))
anova(lme_ALRU_Biomass_RGR_20167890)
emmeans(lme_ALRU_Biomass_RGR_20167890, list(pairwise ~ Treatment_ALRU_67890), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_ALRU_Biomass_AGR_2016_2017 <- lm(log(AGR_Biomass_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_ALRU_Biomass_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_ALRU_Biomass_AGR_2016_2018 <- lm(log(AGR_Biomass_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ALRU_Biomass_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_Biomass_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2017_2018
	)
Treatment_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
Biomass_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03
	)
Tree_ALRU_678 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ALRU_Biomass_AGR_201678 <- lme(log(AGR_ALRU_Biomass_678) ~ 
	Treatment_ALRU_678 + Biomass_ALRU_678,
	random=~1 | Tree_ALRU_678))
anova(lme_ALRU_Biomass_AGR_201678)
emmeans(lme_ALRU_Biomass_AGR_201678, list(pairwise ~ Treatment_ALRU_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_ALRU_Biomass_AGR_2016_2019 <- lm(log(AGR_Biomass_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_ALRU_Biomass_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_Biomass_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_Biomass_2018_2019
	)
Treatment_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
Biomass_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Biomass_est_kg_04
	)
Tree_ALRU_6789 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_ALRU_Biomass_AGR_2016789 <- lme(log(AGR_ALRU_Biomass_6789) ~ 
	Treatment_ALRU_6789 + Biomass_ALRU_6789,
	random=~1 | Tree_ALRU_6789))
anova(lme_ALRU_Biomass_AGR_2016789)
emmeans(lme_ALRU_Biomass_AGR_2016789, list(pairwise ~ Treatment_ALRU_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative AGR
summary(lm_ALRU_Biomass_AGR_2016_2020 <- lm(log(AGR_Biomass_2016_2020) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_ALRU_Biomass_AGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_Biomass_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_Biomass_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_Biomass_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_Biomass_2018_2019,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGR_Biomass_2019_2020
	)
Treatment_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
Biomass_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$Biomass_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$Biomass_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$Biomass_est_kg_04,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$Biomass_est_kg_05
	)
Tree_ALRU_67890 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_ALRU_Biomass_AGR_20167890 <- lme(log(AGR_ALRU_Biomass_67890) ~ 
	Treatment_ALRU_67890 + Biomass_ALRU_67890,
	random=~1 | Tree_ALRU_67890))
anova(lme_ALRU_Biomass_AGR_20167890)
emmeans(lme_ALRU_Biomass_AGR_20167890, list(pairwise ~ Treatment_ALRU_67890), adjust = "tukey")
# a a a a





### 
# ALRU Aboveground Biomass
###

##### RGR

### 2016-2017

summary(lm_ALRU_AGB_RGR_2016_2017 <- lm(RGR_AGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_ALRU_AGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_ALRU_AGB_RGR_2016_2018 <- lm(RGR_AGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ALRU_AGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_AGB_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2017_2018
	)
Treatment_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
AGB_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03
	)
Tree_ALRU_678 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ALRU_AGB_RGR_201678 <- lme(RGR_ALRU_AGB_678 ~ 
	Treatment_ALRU_678 + AGB_ALRU_678,
	random=~1 | Tree_ALRU_678))
anova(lme_ALRU_AGB_RGR_201678)
emmeans(lme_ALRU_AGB_RGR_201678, list(pairwise ~ Treatment_ALRU_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_ALRU_AGB_RGR_2016_2019 <- lm(RGR_AGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_ALRU_AGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_AGB_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_AGB_2018_2019
	)
Treatment_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
AGB_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGB_est_kg_04
	)
Tree_ALRU_6789 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_ALRU_AGB_RGR_2016789 <- lme(RGR_ALRU_AGB_6789 ~ 
	Treatment_ALRU_6789 + AGB_ALRU_6789,
	random=~1 | Tree_ALRU_6789))
anova(lme_ALRU_AGB_RGR_2016789)
emmeans(lme_ALRU_AGB_RGR_2016789, list(pairwise ~ Treatment_ALRU_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative RGR
summary(lm_ALRU_AGB_RGR_2016_2020 <- lm(RGR_AGB_2016_2020 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_ALRU_AGB_RGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_AGB_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_AGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_AGB_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_AGB_2018_2019,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$RGR_AGB_2019_2020
	)
Treatment_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
AGB_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGB_est_kg_04,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGB_est_kg_05
	)
Tree_ALRU_67890 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_ALRU_AGB_RGR_20167890 <- lme(RGR_ALRU_AGB_67890 ~ 
	Treatment_ALRU_67890 + AGB_ALRU_67890,
	random=~1 | Tree_ALRU_67890))
anova(lme_ALRU_AGB_RGR_20167890)
emmeans(lme_ALRU_AGB_RGR_20167890, list(pairwise ~ Treatment_ALRU_67890), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_ALRU_AGB_AGR_2016_2017 <- lm(log(AGR_AGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_ALRU_AGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_ALRU_AGB_AGR_2016_2018 <- lm(log(AGR_AGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ALRU_AGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_AGB_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2017_2018
	)
Treatment_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
AGB_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03
	)
Tree_ALRU_678 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ALRU_AGB_AGR_201678 <- lme(log(AGR_ALRU_AGB_678) ~ 
	Treatment_ALRU_678 + AGB_ALRU_678,
	random=~1 | Tree_ALRU_678))
anova(lme_ALRU_AGB_AGR_201678)
emmeans(lme_ALRU_AGB_AGR_201678, list(pairwise ~ Treatment_ALRU_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_ALRU_AGB_AGR_2016_2019 <- lm(log(AGR_AGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_ALRU_AGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_AGB_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_AGB_2018_2019
	)
Treatment_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
AGB_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGB_est_kg_04
	)
Tree_ALRU_6789 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_ALRU_AGB_AGR_2016789 <- lme(log(AGR_ALRU_AGB_6789) ~ 
	Treatment_ALRU_6789 + AGB_ALRU_6789,
	random=~1 | Tree_ALRU_6789))
anova(lme_ALRU_AGB_AGR_2016789)
emmeans(lme_ALRU_AGB_AGR_2016789, list(pairwise ~ Treatment_ALRU_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative AGR
summary(lm_ALRU_AGB_AGR_2016_2020 <- lm(log(AGR_AGB_2016_2020) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_ALRU_AGB_AGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_AGB_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_AGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_AGB_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_AGB_2018_2019,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGR_AGB_2019_2020
	)
Treatment_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
AGB_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGB_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGB_est_kg_04,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGB_est_kg_05
	)
Tree_ALRU_67890 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_ALRU_AGB_AGR_20167890 <- lme(log(AGR_ALRU_AGB_67890) ~ 
	Treatment_ALRU_67890 + AGB_ALRU_67890,
	random=~1 | Tree_ALRU_67890))
anova(lme_ALRU_AGB_AGR_20167890)
emmeans(lme_ALRU_AGB_AGR_20167890, list(pairwise ~ Treatment_ALRU_67890), adjust = "tukey")
# a a a a




### 
# ALRU Belowground Biomass
###

##### RGR

### 2016-2017

summary(lm_ALRU_BGB_RGR_2016_2017 <- lm(RGR_BGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_ALRU_BGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_ALRU_BGB_RGR_2016_2018 <- lm(RGR_BGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ALRU_BGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_BGB_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2017_2018
	)
Treatment_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
BGB_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03
	)
Tree_ALRU_678 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ALRU_BGB_RGR_201678 <- lme(RGR_ALRU_BGB_678 ~ 
	Treatment_ALRU_678 + BGB_ALRU_678,
	random=~1 | Tree_ALRU_678))
anova(lme_ALRU_BGB_RGR_201678)
emmeans(lme_ALRU_BGB_RGR_201678, list(pairwise ~ Treatment_ALRU_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_ALRU_BGB_RGR_2016_2019 <- lm(RGR_BGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_ALRU_BGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_BGB_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_BGB_2018_2019
	)
Treatment_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
BGB_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$BGB_est_kg_04
	)
Tree_ALRU_6789 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_ALRU_BGB_RGR_2016789 <- lme(RGR_ALRU_BGB_6789 ~ 
	Treatment_ALRU_6789 + BGB_ALRU_6789,
	random=~1 | Tree_ALRU_6789))
anova(lme_ALRU_BGB_RGR_2016789)
emmeans(lme_ALRU_BGB_RGR_2016789, list(pairwise ~ Treatment_ALRU_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative RGR
summary(lm_ALRU_BGB_RGR_2016_2020 <- lm(RGR_BGB_2016_2020 ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_ALRU_BGB_RGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ALRU_BGB_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$RGR_BGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$RGR_BGB_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$RGR_BGB_2018_2019,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$RGR_BGB_2019_2020
	)
Treatment_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
BGB_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$BGB_est_kg_04,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$BGB_est_kg_05
	)
Tree_ALRU_67890 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_ALRU_BGB_RGR_20167890 <- lme(RGR_ALRU_BGB_67890 ~ 
	Treatment_ALRU_67890 + BGB_ALRU_67890,
	random=~1 | Tree_ALRU_67890))
anova(lme_ALRU_BGB_RGR_20167890)
emmeans(lme_ALRU_BGB_RGR_20167890, list(pairwise ~ Treatment_ALRU_67890), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_ALRU_BGB_AGR_2016_2017 <- lm(log(AGR_BGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
	dat$Use_growth_03 == 1,]))
emmeans(lm_ALRU_BGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_ALRU_BGB_AGR_2016_2018 <- lm(log(AGR_BGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2018) & 
	dat$Use_growth_04 == 1,]))
emmeans(lm_ALRU_BGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_BGB_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2017_2018
	)
Treatment_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment
	)
BGB_ALRU_678 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03
	)
Tree_ALRU_678 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID
	))
summary(lme_ALRU_BGB_AGR_201678 <- lme(log(AGR_ALRU_BGB_678) ~ 
	Treatment_ALRU_678 + BGB_ALRU_678,
	random=~1 | Tree_ALRU_678))
anova(lme_ALRU_BGB_AGR_201678)
emmeans(lme_ALRU_BGB_AGR_201678, list(pairwise ~ Treatment_ALRU_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_ALRU_BGB_AGR_2016_2019 <- lm(log(AGR_BGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2019) & 
	dat$Use_growth_05 == 1,]))
emmeans(lm_ALRU_BGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_BGB_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_BGB_2018_2019
	)
Treatment_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment
	)
BGB_ALRU_6789 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$BGB_est_kg_04
	)
Tree_ALRU_6789 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID
	))
summary(lme_ALRU_BGB_AGR_2016789 <- lme(log(AGR_ALRU_BGB_6789) ~ 
	Treatment_ALRU_6789 + BGB_ALRU_6789,
	random=~1 | Tree_ALRU_6789))
anova(lme_ALRU_BGB_AGR_2016789)
emmeans(lme_ALRU_BGB_AGR_2016789, list(pairwise ~ Treatment_ALRU_6789), adjust = "tukey")
# a a a a

### 2016-2020

# Cumulative AGR
summary(lm_ALRU_BGB_AGR_2016_2020 <- lm(log(AGR_BGB_2016_2020) ~ 
	Treatment,
	data=dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2020) & 
	dat$Use_growth_06 == 1,]))
emmeans(lm_ALRU_BGB_AGR_2016_2020, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ALRU_BGB_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$AGR_BGB_2016_2017,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$AGR_BGB_2017_2018,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$AGR_BGB_2018_2019,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$AGR_BGB_2019_2020
	)
Treatment_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$Treatment,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$Treatment
	)
BGB_ALRU_67890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$BGB_est_kg_01,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$BGB_est_kg_03,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$BGB_est_kg_04,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$BGB_est_kg_05
	)
Tree_ALRU_67890 <- 
	as.factor(c(
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2016_2017) & 
		dat$Use_growth_03 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2017_2018) & 
		dat$Use_growth_04 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2018_2019) & 
		dat$Use_growth_05 == 1,]$PID,
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_BGB_2019_2020) & 
		dat$Use_growth_06 == 1,]$PID
	))
summary(lme_ALRU_BGB_AGR_20167890 <- lme(log(AGR_ALRU_BGB_67890) ~ 
	Treatment_ALRU_67890 + BGB_ALRU_67890,
	random=~1 | Tree_ALRU_67890))
anova(lme_ALRU_BGB_AGR_20167890)
emmeans(lme_ALRU_BGB_AGR_20167890, list(pairwise ~ Treatment_ALRU_67890), adjust = "tukey")
# a a a a





#############################################################################
#############################################################################
#############################################################################