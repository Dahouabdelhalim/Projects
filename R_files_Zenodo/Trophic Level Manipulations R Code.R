###################################################
## Nelson et al. 2019, Journal of Animal Ecology R Code
## Purpose: To assess whether elevational changes in aridity mediated aphid performance directly and 
##          indirectly by altering aphid multi-trophic species interactions
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
# 2. Did aphid colony growth depend on aridity, trophic level treatment, and their interaction?
#    2.1 Import and clean data
#    2.2 Statistical analysis
# 3. Did natural enemy presence depend on aridity?
#    3.1 Import and clean data
#    3.2 Statistical analysis
# 4. Did natural enemy number when present depend on aridity?
#    4.1 Import and clean data
#    4.2 Statistical analysis
# 5. Did tending ant presence depend on aridity?
#    5.1 Import and clean data
#    5.2 Statistical analysis
# 6. Did tending ant number when present depend on aridity??
#    6.1 Import and clean data
#    6.2 Statistical analysis

###################################################
## 1. Preliminaries ----
###################################################
## 1.1 Load required packages ----
library(lme4)
library(car)
library(glmmADMB)
library(bbmle)
library(dplyr)
library(multcomp)

## 1.2 Set working directory ----
setwd()

###################################################
## 2. Did aphid colony growth depend on aridity, trophic level treatment, and their interaction? ----
###################################################
## 2.1 Import and clean data ----
# Import full data
r_data <- read.csv("TrophicLevelManipulations.csv", header = T)

# Create data frame for aphid colonies in the two trophic levels treatment only, to test for direct effects of aridity and variation in host plant effects
r_plant_eff <- r_data %>% filter(trophic.treatment == "Two")

# Create data frame for aphid colonies in the two and three trophic levels treatments only, to test for variation in natural enemy effects
r_pred_eff <- r_data %>% filter(trophic.treatment != "Four") # Create a data frame for the 

# Create data frame for aphid colonies in the three and four trophic levels treatments only, to test for variation in mutualist ant effects
r_ant_eff <- r_data %>% filter(trophic.treatment != "Two")

## 2.2 Statistical analysis ----
# Linear mixed effects model testing whether aphid colony per capita growth rates depended on aridity, trophic level treatment, and their interaction
r_model <- lmer(r ~ PC1*trophic.treatment + (1|valley) + (1|plot), contrasts = list(trophic.treatment = contr.sum), data = r_data)
Anova(r_model, type = 3, test.statistic = "F") # Significance test

# Linear mixed effects model testing for the direct effects of aridity and host plant effects on aphid colony per capita growth rates
r_model_plant_eff <- lmer(r ~ PC1 + (1|valley), data = r_plant_eff)
Anova(r_model_plant_eff, test.statistic = "F") # Significance test

# Linear mixed effects model testing for variation in natural enemy effects with aridity
r_model_pred_eff <- lmer(r ~ PC1*trophic.treatment + (1|valley) + (1|valley:plot), contrasts = list(trophic.treatment = contr.sum), data = r_pred_eff)
Anova(r_model_pred_eff, type = 3, test.statistic = "F") # Significance test

# Linear mixed effects model testing for variation in mutualist ant effects with aridity
r_model_ant_eff <- lmer(r ~ PC1*trophic.treatment + (1|valley) + (1|valley:plot), contrasts = list(trophic.treatment = contr.sum), data = r_ant_eff)
Anova(r_model_ant_eff, type = 3, test.statistic = "F") # Significance test

###################################################
## 3. Did natural enemy presence depend on aridity? ----
###################################################
## 3.1 Import and clean data ----
# Import data, removing aphid colonies where predators were excluded and creating a variable listing whether natural enemies were present (vs. absent) at an aphid colony
enemy_pres_data <- read.csv("TrophicLevelManipulations.csv", header = T) %>% 
  filter(trophic.treatment != "Two") %>%
  mutate(enemy.pres = ifelse(num.enemies == "0", "0", "1"))

# Define variables
enemy_pres_data$enemy.pres <- as.factor(as.character(enemy_pres_data$enemy.pres))

## 3.2 Statistical analysis ----
# Generalized linear mixed model testing whether the presence of natural enemies depended on aridity, the mean number of aphids, trophic level treatment (three vs. four), and sampling effort
enemy_pres_model <- glmmadmb(enemy.pres ~ PC1 + mean.aphids + trophic.treatment + num.obs + (1|valley) + (1|valley:plot), family = "binomial", data = enemy_pres_data)

# Significance test
summary(enemy_pres_model)$coefficients

###################################################
## 4. Did natural enemy number when present depend on aridity? ----
###################################################
## 4.1 Import and clean data ----
# Import data, removing aphid colonies with no observed natural enemies
enemy_num_data <- enemy_pres_data %>% filter(num.enemies > 0)

## 4.2 Statistical analysis ----
# Alternative generalized linear mixed models testing whether the number of natural enemies when present depended on aridity, the mean number of aphids, and their interaction, as well as the main effects of trophic level treatment and sampling effort
enemy_num_model_1 <- glmmadmb(num.enemies ~ PC1*mean.aphids + trophic.treatment + num.obs + (1|valley) + (1|valley:plot), family = "truncpoiss", data = enemy_num_data)

enemy_num_model_2 <- glmmadmb(num.enemies ~ PC1*mean.aphids + trophic.treatment + num.obs + (1|valley) + (1|valley:plot), family = "truncnbinom1", data = enemy_num_data)

enemy_num_model_3 <- glmmadmb(num.enemies ~ PC1*mean.aphids + trophic.treatment + num.obs + (1|valley) + (1|valley:plot), family = "truncnbinom", data = enemy_num_data)

# Select best-fitting model by comparing AIC values
AICtab(enemy_num_model_1, enemy_num_model_2, enemy_num_model_3)

# Significance test of best-fitting model
summary(enemy_num_model_2)$coefficients

# Remove the interaction to test for the significance of the main effects
enemy_num_model_4 <- glmmadmb(num.enemies ~ PC1 + mean.aphids + trophic.treatment + num.obs + (1|valley) + (1|valley:plot), family = "truncnbinom1", data = enemy_num_data)

# Significance test
summary(enemy_num_model_4)$coefficients

###################################################
## 5. Did tending ant presence depend on aridity? ----
###################################################
## 5.1 Import and clean data ----
# Import data, only including aphid colonies in the four trophic level treatment
ant_pres_data <- read.csv("TrophicLevelManipulations.csv", header = T) %>% 
  filter(trophic.treatment == "Four") %>%
  mutate(ant.pres = ifelse(num.ants == "0", "0", "1"))

# Define variables
ant_pres_data$ant.pres <- as.factor(as.character(ant_pres_data$ant.pres))

## 5.2 Statistical analysis ----
# Generalized linear mixed model testing whether the presence of tending ants depended on aridity and the mean number of aphids per colony
ant_pres_model <- glmmadmb(ant.pres ~ PC1 + mean.aphids + (1|valley), family = "binomial", data = ant_pres_data)

# Significance test
summary(ant_pres_model)$coefficients

###################################################
## 6. Did tending ant number when present depend on aridity? ----
###################################################
## 6.1 Import and clean data ----
# Import data, removing aphid colonies without ants present
ant_num_data <- ant_pres_data %>% filter(ant.pres != "0")

## 6.2 Statistical analysis ----
# Alternative generalized linear mixed models testing whether the number of tending ants depended on aridity, the maximum number of aphids per colony, and their interaction
ant_num_model_1 <- glmmadmb(num.ants ~ PC1*max.tended.aphids + (1|valley), data = ant_num_data, family = "truncpoiss")

ant_num_model_2 <- glmmadmb(num.ants ~ PC1*max.tended.aphids + (1|valley), data = ant_num_data, family = "truncnbinom1")

ant_num_model_3 <- glmmadmb(num.ants ~ PC1*max.tended.aphids + (1|valley), data = ant_num_data, family = "truncnbinom")

# Select the best-fitting model by comparing AIC values
AICtab(ant_num_model_1, ant_num_model_2, ant_num_model_3)

# Significance test of the best-fitting model
summary(ant_num_model_2)$coefficients
