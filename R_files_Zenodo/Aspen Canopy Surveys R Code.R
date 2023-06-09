###################################################
## Nelson et al. 2019, Journal of Animal Ecology R Code
## Purpose: To assess whether aphid and ant abundance in aspen canopy surveys varied along an elevational
##          gradient in aridity
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
# 2. Did aphid presence depend on aridity?
#    2.1 Import and clean data
#    2.2 Statistical analysis
# 3. Did aphid colony number when present depend on aridity?
#    3.1 Import and clean data
#    3.2 Statistical analysis
# 4. Did aphid colony size when present depend on aridity?
#    4.1 Import and clean data
#    4.2 Statistical analysis
# 5. Did tending ant presence depend on aridity?
#    5.1 Import and clean data
#    5.2 Statistical analysis
# 6. Did tending ant number when present depend on aridity?
#    6.1 Import and clean data
#    6.2 Statistical analysis
# 7. Did roaming ant presence depend on aridity?
#    7.1 Import and clean data
#    7.2 Statistical analysis
# 8. Did roaming ant number when present depend on aridity?
#    8.1 Import and clean data
#    8.2 Statistical analysis

###################################################
## 1. Preliminaries ----
###################################################
## 1.1 Load required packages ----
library(lme4)
library(car)
library(glmmADMB)
library(bbmle)
library(dplyr)

## 1.2 Set working directory ----
setwd()

###################################################
## 2. Did aphid presence depend on aridity? ----
###################################################
## 2.1 Import and clean data ----
# Import data and add a variable describing whether aphid colonies were present (vs. absent) in a tree
colony_pres_data <- read.csv("AspenCanopySurveys.csv", header = T) %>%
  filter(tree.id != "310") %>% # Exclude one tree with an unusually large number of aphid colonies
  mutate(colony.presence = ifelse(num.aphids == "0", "0", "1")) %>%
  group_by(elevation, valley, plot, elevation.m, PC1, tree.id, search.time) %>%
  summarize(colony.presence = max(as.numeric(as.character(colony.presence))))

# Define variables
colony_pres_data$colony.presence <- as.factor(as.character(colony_pres_data$colony.presence))
  
# 2.2 Statistical analysis ----
# Generalized linear mixed model testing whether aphid colony presence depended on aridity after accounting for search time
presence_model <- glmmadmb(colony.presence ~ PC1 + search.time + (1|valley) + (1|plot), data = colony_pres_data, family = "binomial") 

# Significance test
summary(presence_model)$coefficients

###################################################
## 3. Did aphid colony number when present depend on aridity? ----
###################################################
## 3.1 Import and clean data ----
# Import data and add a variable describing the number of aphid colonies per tree
colony_num_data <- read.csv("AspenCanopySurveys.csv", header = T) %>%
  mutate(colony.num = ifelse(num.aphids == "0", "0", "1")) %>%
  group_by(elevation, valley, plot, elevation.m, PC1, tree.id, search.time) %>%
  summarize(colony.num = sum(as.numeric(as.character(colony.num)))) %>%
  filter(colony.num > 0 & tree.id != "310") # Exclude trees with no aphids and one tree with an unusually large number of aphid colonies

## 3.2 Statistical analysis ----
# Alternative generalized linear mixed models testing whether aphid colony number per tree depended on aridity after accounting for search time
num_model_1 <- glmmadmb(colony.num ~ PC1 + search.time + (1|valley) + (1|plot), data = colony_num_data, family = "truncpoiss") # truncated-at-zero Poisson

num_model_2 <- glmmadmb(colony.num ~ PC1 + search.time + (1|valley) + (1|plot), data = colony_num_data, family = "truncnbinom1") # truncated-at-zero negative binomial with a quasi-Poisson scale parameter

num_model_3 <- glmmadmb(colony.num ~ PC1 + search.time + (1|valley) + (1|plot), data = colony_num_data, family = "truncnbinom") # truncated-at-zero negative binomial

# Select best-fitting model by comparing AIC values
AICtab(num_model_1, num_model_2, num_model_3)

#Significance test of the best-fitting model
summary(num_model_3)$coefficients

###################################################
## 4. Did aphid colony size when present depend on aridity? ----
###################################################
## 4.1 Import and clean data ----
# Import data, removing trees without aphids and trees with unusually high numbers of aphids per colony or aphid colonies per tree
colony_size_data <- read.csv("AspenCanopySurveys.csv", header = T) %>% 
  filter(num.aphids > 0 & num.aphids < 170 & tree.id != "310") 

# Define variables
colony_size_data$tree.id <- as.factor(as.character(colony_size_data$tree.id))

## 4.2 Statistical analysis ----
# Alternative generalized linear mixed models testing whether aphid colony size depended on aridity after accounting for search time
size_model_1 <- glmmadmb(num.aphids ~ PC1 + search.time + (1|valley) + (1|plot) + (1|tree.id), data = colony_size_data, family = "truncpoiss") # truncated-at-zero Poisson

size_model_2 <- glmmadmb(num.aphids ~ PC1 + search.time + (1|valley) + (1|plot) + (1|tree.id), data = colony_size_data, family = "truncnbinom1") # truncated-at-zero negative binomial with a quasi-Poisson scale parameter

size_model_3 <- glmmadmb(num.aphids ~ PC1 + search.time + (1|valley) + (1|plot) + (1|tree.id), data = colony_size_data, family = "truncnbinom") # truncated-at-zero negative binomial; model doesn't run

# Select best-fitting model by comparing AIC values
AICtab(size_model_1, size_model_2, size_model_3)

# Significance test of the best-fitting model
summary(size_model_3)$coefficients

###################################################
## 5. Did tending ant presence depend on aridity? ----
###################################################
## 5.1 Import and clean data ----
# Import data, creating a variable describing whether tending ants were present (vs. absent) and excluding trees without aphids and with three unusually large aphid colonies
tending_pres_data <- read.csv("AspenCanopySurveys.csv", header = T) %>%
  mutate(ant.tended = ifelse(num.ants == "0", "0", "1")) %>% 
  filter(num.aphids > 0 & num.aphids < 170)

# Define variables
tending_pres_data$ant.tended <- as.factor(as.character(tending_pres_data$ant.tended))
tending_pres_data$tree.id <- as.factor(as.character(tending_pres_data$tree.id))

## 5.2 Statistical analysis ----
# Generalized linear mixed model testing whether tending ant presence depended on aridity and the number of aphids per colony
tending_pres_model <- glmmadmb(ant.tended ~ PC1 + num.aphids + (1|valley) + (1|plot) + (1|tree.id), family = "binomial", data = tending_pres_data)

# Significance test
summary(tending_pres_model)$coefficients

###################################################
## 6. Did tending ant number when present depend on aridity? ----
###################################################
## 6.1 Import and clean data ----
# Import data, removing excluding trees without ants or aphids as well as with three unusually large aphid colonies
tending_num_data <- read.csv("AspenCanopySurveys.csv", header = T) %>%
  filter(num.aphids > 0 & num.aphids < 170 & num.ants > 0)

# Define variables
tending_num_data$tree.id <- as.factor(as.character(tending_num_data$tree.id))

## 6.2 Statistical analysis ----
# Alternative generalized linear mixed models testing whether tending ant number depended on aridity, aphid number per colony, and their interaction
tend_model_1 <- glmmadmb(num.ants ~ PC1*num.aphids + (1|valley) + (1|plot) + (1|tree.id), data = tending_num_data, family = "truncpoiss") # truncated-at-zero Poisson

tend_model_2 <- glmmadmb(num.ants ~ PC1*num.aphids + (1|valley) + (1|plot) + (1|tree.id), data = tending_num_data, family = "truncnbinom1") # truncated-at-zero negative binomial with a quasi-Poisson scale parameter

tend_model_3 <- glmmadmb(num.ants ~ PC1*num.aphids + (1|valley) + (1|plot) + (1|tree.id), data = tending_num_data, family = "truncnbinom") # truncated-at-zero negative binomial

# Select best-fitting model by comparing AIC values
AICtab(tend_model_1, tend_model_2, tend_model_3)

#Significance test of the best-fitting model
summary(tend_model_2)$coefficients

# Generalized linear mixed model with the non-significant interaction removed to test for the significance of the main effects
tend_model_4 <- glmmadmb(num.ants ~ PC1 + num.aphids + (1|valley) + (1|plot) + (1|tree.id), data = tending_num_data, family = "truncnbinom1")

# Significance test
summary(tend_model_4)$coefficients

###################################################
## 7. Did roaming ant presence depend on aridity? ----
###################################################
## 7.1 Import and clean data ----
# Import data, removing trees with aphids
roam_pres_data <- read.csv("AspenCanopySurveys.csv", header = T) %>%
  mutate(roam.pres = ifelse(num.ants == "0", "0", "1")) %>% 
  filter(num.aphids == "0" & other.aphid.sp == "N")

# Define variables
roam_pres_data$roam.pres <- as.factor(as.character(roam_pres_data$roam.pres))

## 7.2 Statistical analysis ----
# Generalized linear mixed model testing whether the presence of roaming ants depended on aridity after accounting for search time
roam_pres_model <- glmmadmb(roam.pres ~ PC1 + search.time + (1|valley) + (1|plot), data = roam_pres_data, family = "binomial")

# Significance test
summary(roam_pres_model)$coefficients

###################################################
## 8. Did roaming ant number when present depend on aridity? ----
###################################################
## 8.1 Import and clean data ----
# Import data, removing trees with aphids
roam_num_data <- read.csv("AspenCanopySurveys.csv", header = T) %>%
  filter(num.aphids == "0" & other.aphid.sp == "N" & num.ants > 0)

## 8.2 Statistical analysis ----
# Alternative generalized linear mixed models testing whether roaming ant number depended on aridity, after accounting for search time
roam_model_1 <- glmmadmb(num.ants ~ PC1 + search.time + (1|elevation:valley) + (1|valley:plot), data = roam_num_data, family = "truncpoiss") # truncated-at-zero Poisson

roam_model_2 <- glmmadmb(num.ants ~ PC1 + search.time + (1|elevation:valley) + (1|valley:plot), data = roam_num_data, family = "truncnbinom1") # truncated-at-zero negative binomial with a quasi-Poisson scale parameter

roam_model_3 <- glmmadmb(num.ants ~ PC1 + search.time + (1|elevation:valley) + (1|valley:plot), data = roam_num_data, family = "truncnbinom") # truncated-at-zero negative binomial with a quasi-Poisson scale parameter

# Select best-fitting model by comparing AIC values
AICtab(roam_model_1, roam_model_2, roam_model_3)

# Significance test of the best-fitting model
summary(roam_model_2)$coefficients
