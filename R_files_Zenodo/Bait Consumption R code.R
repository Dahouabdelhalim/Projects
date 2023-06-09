###################################################
## Nelson et al. 2018 Oikos R Code
## Purpose: To assess whether ant carbohydrate consumption as indicated by recruitment to baits 
##          depended on elevation or mean summer temperature
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
# 2. Did ant bait consumption rate depend on elevation and sugar concentration?
#    2.1 Import and clean data
#    2.2 Statistical analysis
# 3. Did ant forager abundance in baits depend on elevation and sugar concentration? 
#    2.1 Import and clean data
#    2.2 Statistical analysis
# 4. Did bait consumption per observed ant depend on elevation or mean summer temperature and sugar 
#    concentration?
#    2.1 Import and clean data
#    2.2 Statistical analysis

#####################################################
## 1. Preliminaries ----
#####################################################
## 1.1 Load required packages ----
library(dplyr)
library(lme4)
library(car)
library(glmmADMB)
library(bbmle)

## 1.2 Set working directory ----
setwd()

#####################################################
## 2. Did ant bait consumption rate depend on elevation and sugar concentration? ----
#####################################################
## 2.1 Import and clean data ----
# Import data, removing the water controls and observations with missing data
bait_data <- read.csv("BaitConsumption.csv", header = T) %>%
  filter(corr.cons.rate != "NA" & num.fpodz != "NA" & sugar.concentration != "0")

# When evaporation rates of water controls exceed consumption rates of sugar baits (producing negative adjusted consumption rates), set sugar consumption rates to 0.
bait_data <- bait_data %>% 
  mutate(corr.cons.rate = ifelse(as.numeric(as.character(corr.cons.rate)) <= 0, 0, corr.cons.rate))

# Define variables
bait_data$sugar.concentration <- as.factor(as.character(bait_data$sugar.concentration))
bait_data$date <- as.factor(as.character(bait_data$date))

## 2.2 Statistical analysis
# Linear mixed model testing whether ant bait consumption rates depended on elevation (high vs. low), sugar concentration, and their interaction
cons_rate_model <- lmer(log(corr.cons.rate + 1) ~ elevation*sugar.concentration + date + (1|elevation:valley) + (1|valley:mound.id), data = bait_data)
# Significance test
Anova(cons_rate_model, type = 3)

#####################################################
## 3. Did ant forager abundance in baits depend on elevation and sugar concentration? ----
#####################################################
## 3.1 Import and clean data
# Import data, removing the water controls and observations with missing data
ant_abund <- read.csv("BaitConsumption.csv", header = T) %>%
  filter(num.fpodz != "NA" & sugar.concentration != "0")

# Define variables
ant_abund$sugar.concentration <- as.factor(as.character(ant_abund$sugar.concentration))
ant_abund$date <- as.factor(as.character(ant_abund$date))

## 3.2 Statistical analysis
# Alternate possible generalized linear mixed models testing whether ant forager abundance in baits depended on elevation (high vs. low), sugar concentration, and their interaction
ant_abund_model1 <- glmmadmb(num.fpodz ~ elevation*sugar.concentration + date + (1|elevation:valley) + (1|valley:mound.id), data = ant_abund, family = "gaussian") # normal distribution
ant_abund_model2 <- glmmadmb(num.fpodz ~ elevation*sugar.concentration + date + (1|elevation:valley) + (1|valley:mound.id), data = ant_abund, family = "nbinom") # negative binomial distribution
ant_abund_model3 <- glmmadmb(num.fpodz ~ elevation*sugar.concentration + date + (1|elevation:valley) + (1|valley:mound.id), data = ant_abund, family = "nbinom1") # negative binomial distribution with a quasi-Poisson scale parameter
# Compare AIC values to select the best-fitting model
AICtab(ant_abund_model1, ant_abund_model2, ant_abund_model3)
# Significance test of the best-fitting model
Anova(ant_abund_model3, type = 3)

# Generalized linear mixed model with the non-significant interaction term removed to test for the significance of the main effects
ant_abund_model4 <- glmmadmb(num.fpodz ~ elevation + sugar.concentration + date + (1|elevation:valley) + (1|valley:mound.id), data = ant_abund, family = "nbinom1")
# Significance test
Anova(ant_abund_model4, type = 3)

#####################################################
## 4. Did bait consumption per observed ant depend on elevation or mean summer temperature and sugar 
##    concentration? ----
#####################################################
## 4.1 Import and clean data
# Add bait consumption rate per ant as a variable in the dataframe
cons_per_ant_data <- bait_data %>% mutate(cons.per.ant = corr.cons.rate/(num.fpodz + 1))

## 4.2 Statistical analysis
# Linear mixed model testing whether consumption rates per ant depended on elevation (high vs. low), sugar concentration, and their interaction
cons.ant.model <- lmer(log(cons.per.ant + 1) ~ elevation*sugar.concentration + date + (1|elevation:valley) + (1|valley:mound.id), data = cons_per_ant_data)
Anova(cons.ant.model, type = 3)

# Linear mixed model testing whether consumption rates per ant depended on mean summer temperature, sugar concentration, and their interaction
cons.ant.model <- lmer(log(cons.per.ant + 1) ~ summer.temp*sugar.concentration + date + (1|valley) + (1|valley:mound.id), data = cons_per_ant_data)
Anova(cons.ant.model, type = 3)
