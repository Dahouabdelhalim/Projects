###################################################
## Nelson et al. 2018 Oikos R Code
## Purpose: To assess whether aphid natural enemy abundance varied between one of the low (Spring Creek)
##          and one of the high elevation valleys (East River Valley)
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
# 2. Does natural enemy abundance differ between the high and low elevation valleys?
#    2.1 Import and clean data
#    2.2 Statistical analysis

###################################################
## 1. Preliminaries ----
###################################################
## 1.1 Load required packages ----
library(dplyr)
library(lme4)
library(car)

## 1.2 Set working directory ----
setwd()

#####################################################
## 2. Does natural enemy abundance differ between the high and low elevation valleys? ----
#####################################################
## 2.1 Import and clean data ----
# Import data and calculate the mean number of aphids and natural enemies per plant across all observations
enemy_data <- read.csv("NaturalEnemyAbundance.csv", header = T) %>%
  group_by(valley, plant.id, block) %>%
  summarize(mean.aphids = mean(num.aphids), mean.enemies = mean(num.enemies))

# Define variables
enemy_data$block <- as.factor(enemy_data$block)

# Remove plants where aphids were never present
enemy_data <- enemy_data %>% filter(mean.aphids > 0)

## 2.2 Statistical analysis ----
# Linear mixed model testing whether the mean number of natural enemies differed between valleys and with mean aphid number
enemy_model <- lmer(log(mean.enemies + 1) ~ valley + mean.aphids + (1|block), data = enemy_data)
Anova(enemy_model, type = 3)
