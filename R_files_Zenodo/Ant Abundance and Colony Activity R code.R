###################################################
## Nelson et al. 2018 Oikos R Code
## Purpose: To assess whether ant abundance and colony activity depended on elevation (high vs. low) or 
##          mean summer temperature
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
# 2. Did Formica podzolica ant abundance in pitfall traps differ with elevation or mean summer 
#    temperature?
#    2.1 Import and clean data
#    2.2 Statistical analysis
# 3. Did Formica podzolica ant nest mound surface area differ with elevation or mean summer 
#    temperature?
#    3.1 Import and clean data
#    3.2 Statistical analysis
# 4. Did Formica podzolica ant activity on the mound surface differ with elevation or mean summer 
#    temperature?
#    4.1 Import and clean data
#    4.2 Statistical analysis

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
## 2. Did Formica podzolica ant abundance in pitfall traps differ with elevation or mean summer 
##    temperature? ----
#####################################################
## 2.1 Import and clean data ----
trap_data <- read.csv("AntAbundance.csv", header = T)
# Calculate the number of ants per trap-day
trap_data$ants.per.trapday <- as.numeric(trap_data$num.fpodz/trap_data$trap.days)

# Create separate data frame for data collected in 2012
traps_2012 <- trap_data %>% filter(year == "2012")
# Create separate data frame for data collected in 2015
traps_2015 <- trap_data %>% filter(year == "2015")

## 2.2 Statistical analysis ----
# Linear mixed model testing whether ant abundance depended on elevation (high vs. low) in 2012
traps_2012_model <- lmer((ants.per.trapday)^(1/3) ~ elevation + (1|elevation:valley), data = traps_2012)
# Significance test
Anova(traps_2012_model, type = 3)

# Linear mixed model testing whether ant abundance depended on elevation (high vs. low) in 2015
traps_2015_model <- lmer((ants.per.trapday)^(1/3) ~ elevation + (1|elevation:valley), data = traps_2015)
# Significance test
Anova(traps_2015_model, type = 3)

# Linear mixed model testing whether ant abundance depended on mean summer temperature in 2015
traps_temp_model <- lmer((ants.per.trapday)^(1/3) ~ summer.temp + (1|valley), data = traps_2015)
# Significance test
Anova(traps_temp_model, type = 3)

#####################################################
## 3. Did Formica podzolica ant nest mound surface area differ with elevation or mean summer 
##    temperature? ----
#####################################################
## 3.1 Import and clean data ----
# Import data and summarize so that the surface area of each ant nest mound is listed in only one row
mound_area_data <- read.csv("MoundSizeActivity.csv", header = T) %>%
  group_by(elevation, summer.temp, valley, mound.id) %>%
  summarize(mound.area = mean(mound.area))
# Remove two unusually large ant nest mounds from the dataframe
mounds_no_outliers <- mound_area_data %>% filter(mound.area < 1)

## 3.2 Statistical analysis ----
# Linear mixed model testing whether ant nest mound surface area depended on elevation (high vs. low)
area_model <- lmer(log(mound.area) ~ elevation + (1|elevation:valley), data = mounds_no_outliers)
# Significance test
Anova(area_model, type = 3)

# Linear mixed model testing whether ant nest mound surface area depended on mean summer temperature
area_temp_model <- lmer(log(mound.area) ~ summer.temp + (1|valley), data = mounds_no_outliers)
# Significance test
Anova(area_temp_model, type = 3)

#####################################################
## 4. Did Formica podzolica ant activity on the mound surface differ with elevation or mean 
##    summer temperature? ----
#####################################################
## 4.1 Import and clean data ----
# Create data frame with the two unusually large ant nest mounds removed, and calculate the mean number of ants observed on each mound surface across all observations
activity_data <- read.csv("MoundSizeActivity.csv", header = T) %>% 
  filter(num.fpodz != "NA" & mound.area < 1) %>%
  group_by(elevation, summer.temp, valley, mound.id) %>%
  summarize(mean.fpodz = mean(num.fpodz))

## 4.2 Statistical analysis ----
# Linear mixed model testing whether the mean number of ants on the mound surface depended on elevation (high vs. low)
activity_model <- lmer(log(mean.fpodz) ~ elevation + (1|elevation:valley), data = activity_data)
# Significance test
Anova(activity_model, type = 3)

# Linear mixed model testing whether the mean number of ants on the mound surface depended on mean summer temperature
activity_temp_model <- lmer(log(mean.fpodz) ~ summer.temp + (1|valley), data = activity_data)
# Significance test
Anova(activity_temp_model, type = 3)
