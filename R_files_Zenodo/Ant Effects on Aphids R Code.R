###################################################
## Nelson et al. 2018 Oikos R Code
## Purpose: To assess whether the effects of ants on aphids depended on elevation or mean summer
##          temperature
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
# 2. Does the presence of ants tending aphids depend on elevation or mean summer temperature?
#    2.1 Import and clean data
#    2.2 Statistical analysis
# 3. Does ant per capita tending rate of aphids depend on elevation or mean summer temperature?
#    3.1 Import and clean data
#    3.2 Statistical analysis
# 4. Do the effects of ants on aphid colony survival depend on elevation or mean summer temperature?
#    4.1 Import and clean data
#    4.2 Statistical analysis
# 5. Do the effects of ants on aphid colony growth rate depend on elevation or mean summer temperature?
#    5.1 Import and clean data
#    5.2 Statistical analysis

#####################################################
## 1. Preliminaries ----
#####################################################
## 1.1 Load required packages ----
library(dplyr)
library(lme4)
library(car)

## 1.2 Set working directory ----
setwd()

#####################################################
## 2. Does the presence of ants tending aphids depend on elevation or mean summer temperature? ----
#####################################################
## 2.1 Import and clean data
# Import data for aphid colonies in the ant presence treatment only, removing missing observations
presence_data <- read.csv("AntEffectsOnAphids.csv", header = T) %>% 
  filter(num.fpodz != "NA" & num.aphids != "NA" & ant.treatment == "ants")

# Create variables indicating the mean number of ants and aphids across all observations as well as whether tending ants were ever present (vs. absent) at an aphid colony
presence_data <- presence_data %>%
  group_by(elevation, valley, mound.id, ant.treatment, summer.temp) %>%
  summarize(mean.aphids = mean(as.numeric(as.character(num.aphids))),
            mean.ants = mean(as.numeric(as.character(num.fpodz))),
            ants.present = ifelse(sum(as.numeric(as.character(num.fpodz))) > 0, "1", "0"))

# Define variables
presence_data$ants.present <- as.factor(as.character(presence_data$ants.present))

## 2.2 Statistical analysis
# Generalized linear model with a binomial distribution testing whether ant presence (vs. absence) depended on elevation and the mean number of aphids per colony
present_model <- glmmadmb(ants.present ~ elevation + mean.aphids + (1|elevation:valley), data = presence_data, family = "binomial")
# Significance test
Anova(present_model, type = 3)

# Generalized linear model with a binomial distribution testing whether ant presence (vs. absence) depended on mean summer temperature and the mean number of aphids per colony
present_temp_model <- glmmadmb(ants.present ~ summer.temp + mean.aphids + (1|elevation:valley), data = presence_data, family = "binomial")
# Significance test
Anova(present_temp_model, type = 3)

#####################################################
## 3. Does ant per capita tending rate of aphids depend on elevation or mean summer temperature? ----
#####################################################
## 3.1 Import and clean data
# Create data frame that only includes aphid colonies that were ant tended
tending_rate_data <- presence_data %>% filter(ants.present != "0")

## 3.2 Statistical analysis
# Linear mixed model testing whether the mean number of tending ants across all observations depended on elevation, the mean number of aphids per colony, and their interaction
tending_rate_model <- lmer(log(mean.ants) ~ mean.aphids*elevation + (1|elevation:valley), data = tending_rate_data)
# Significance test
Anova(tending_rate_model, type = 3)

# Linear mixed model with the non-significant interaction removed to test for the significance of the main effects
tending_rate_model_2 <- lmer(log(mean.ants) ~ mean.aphids + elevation + (1|elevation:valley), data = tending_rate_data)
# Significance test
Anova(tending_rate_model_2, type = 3)

# Linear mixed model testing whether the mean number of tending ants across all observations depended on mean summer temperature, the mean number of aphids per colony, and their interaction
tending_rate_temp_model <- lmer(log(mean.ants) ~ mean.aphids*summer.temp + (1|valley), data = tending_rate_data)
# Significance test
Anova(tending_rate_temp_model, type = 3)

#####################################################
## 4. Do the effects of ants on aphid colony survival depend on elevation or mean summer temperature? ----
#####################################################
## 4.1 Import and clean data
# Import data, removing missing observations and aphid colonies where multiple (>1) ants were observed to have breached the exclusions
survival_data <- read.csv("AntEffectsOnAphids.csv", header = T) %>% filter(survived != "NA") %>%
  mutate(ant.sneaker = ifelse(ant.treatment == "no ants" & ant.tended == "1", "1", "0")) %>%
  filter(ant.sneaker != "1")

# Define variables
survival_data$survived <- as.factor(as.character(survival_data$survived))

## 4.2 Statistical analysis
# Generalized linear model with a binomial distribution testing whether aphid colony survival (vs. extinction) depended on elevation, ant tending treatment, and their interaction
survival_model <- glmmadmb(survived ~ elevation*ant.treatment + (1|elevation:valley) + (1|valley:mound.id), family = "binomial", data = survival_data)
# Significance test
Anova(survival_model, type = 3)

# Generalized linear model with a binomial distribution testing whether ants affected aphid colony survival at high elevations
survival_model_high <- glmmadmb(survived ~ ant.treatment + (1|elevation:valley) + (1|valley:mound.id), family = "binomial", data = subset(survival_data, elevation == "High"))
Anova(survival_model_high, type = 3)

# Generalized linear model with a binomial distribution testing whether ants affected aphid colony survival at low elevations
survival_model_low <- glmmadmb(survived ~ ant.treatment + (1|elevation:valley) + (1|valley:mound.id), family = "binomial", data = subset(survival_data, elevation == "Low"))
Anova(survival_model_low, type = 3)

# Generalized linear model with a binomial distribution testing whether aphid colony survival (vs. extinction) depended on mean summer temperature, ant tending treatment, and their interaction
survival_temp_model <- glmmadmb(survived ~ summer.temp*ant.treatment + (1|valley) + (1|valley:mound.id), family = "binomial", data = survival_data)
# Significance test
Anova(survival_temp_model, type = 3)

#####################################################
## 5. Do the effects of ants on aphid colony growth rate depend on elevation or mean summer 
##    temperature? ----
#####################################################
## 5.1 Import and clean data
# Create data frame that only includes aphid colonies that survived
growth_data <- survival_data %>% filter(survived == "1")

## 5.2 Statistical analysis
# Linear mixed model testing whether aphid colony growth (r) depended on elevation (high vs. low), ant tending treatment, and their interaction
growth_model <- lmer(r ~ elevation*ant.treatment + (1|elevation:valley) + (1|valley:mound.id), data = growth_data)
# Significance test
Anova(growth_model, type = 3)

# Linear mixed model with the non-significant interaction removed to test for the significance of the main effects
growth_model_2 <- lmer(r ~ elevation + ant.treatment + (1|elevation:valley) + (1|valley:mound.id), data = growth_data)
# Significance test
Anova(growth_model_2, type = 3)

# Linear mixed model testing whether aphid colony growth (r) depended on mean summer temperature, ant tending treatment, and their interaction
growth_temp_model <- lmer(r ~ summer.temp*ant.treatment + (1|elevation:valley) + (1|valley:mound.id), data = growth_data)
# Significance test
Anova(growth_temp_model, type = 3)

# Linear mixed model with the non-significant interaction removed to test for the significance of the main effects
growth_temp_model_2 <- lmer(r ~ summer.temp + ant.treatment + (1|elevation:valley) + (1|valley:mound.id), data = growth_data)
# Significance test
Anova(growth_temp_model_2, type = 3)
