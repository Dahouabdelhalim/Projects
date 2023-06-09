###################################################
## Nelson et al. 2017 R Code
## Purpose: Determine the mechanisms accounting for day/night variation in the critical thermal maxima of Ectatomma ruidum ant foragers
## Corresponding author:
## Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#   1.1 Load required packages
#   1.2 Set working directory
# 2. Did diurnally foraging ants have a higher CTmax than those foraging nocturnally, and were 
#    there differences in CTmax among ant colonies?
#   2.1 Import & clean data
#   2.2 Statistical analysis
# 3. Was there among-colony day/night partitioning of foraging?
#   3.1 Import & clean data
#   3.2 Statistical analysis
# 4. Was there within-colony day/night partitioning of foraging? 
#   4.1 Statistical analysis
#      4.1.1 Test whether the number of colonies where morning-marked foragers were recollected
#      depended on the time of day of collection
#      4.1.2 Test whether the number of colonies where afternoon-marked foragers were recollected
#      depended on the time of day of collection
# 5. Did ant foragers acclimate to experimental manipulations of temperature after 3 hours?
#   5.1 Import & clean data
#   5.2 Statistical analysis

# ###################################################
## 1. Preliminaries ----
#####################################################
## 1.1 Load required packages ----
library(dplyr)
library(lme4) 
library(car)

## 1.2 Set working directory ----
setwd()

###################################################
## 2. Did diurnally foraging ants have a higher CTmax than those foraging nocturnally, and were there differences in CTmax among ant colonies? ----
###################################################
## 2.1 Import & clean data ----
# Import the data
diel_data <- read.csv("DayNightCTmax.csv", header = T) %>% 
  filter(ant.colony.id != "6") # exclude Colony #6, where only two ants were collected
# Define Colony ID as a factor
diel_data$ant.colony.id <- as.factor(diel_data$ant.colony.id) 

## 2.2 Statistical analysis ----
# Generalized linear model testing whether CTmax depended on the time when ants were collected, ant source colony, and their interaction
diel_model <- glm(CT.max ~ time.collected*ant.colony.id, data = diel_data, family = quasipoisson(link = "log"))
# Significance test
Anova(diel_model, test.statistic = "F")

# Generalized linear model, with the interaction removed to test for the significance of the main effects
diel_model <- glm(CT.max ~ time.collected + ant.colony.id, data = diel_data, family = quasipoisson(link = "log")) 
# Significance test
Anova(diel_model, test.statistic = "F")


###################################################
## 3. Was there among-colony day/night partitioning of foraging? ----
###################################################
## 3.1 Import & clean data ---- 
partitioning_data <- read.csv("AntForagerNumber.csv", header = T)

## 3.2 Statistical analysis ---- 
# Linear model testing whether the number of ants collected diurnally depended on the number of ants collected nocturnally
model <- lm(diurnal.ant.num ~ nocturnal.ant.num, data = partitioning_data)
# Significance test
Anova(model, type = 3, test.statistic = "F")

# Test for normality of residuals
shapiro.test(resid(model)) 
hist(resid(model))


###################################################
## 4. Was there within-colony day/night partitioning of foraging? ----
###################################################
## 4.1 Statistical analysis ----
## 4.1.1 Test whether the number of colonies where morning-marked foragers were recollected depended on the time of day of collection
# These data are the number of colonies where morning-marked foragers were collected or not collected, both diurnally and nocturnally, respectively
n_morning <- c(5, 5, 5, 5)
# Format data for Fisher's exact test
n_colonies <- matrix(n_morning,nrow=2) 
# Fisher's exact test
fisher.test(n_colonies) 

## 4.1.2 Test whether the number of colonies where afternoon-marked foragers were recollected depended on the time of day of collection ----
# These data are the number of colonies where afternoon-marked foragers were collected or not collected, both diurnally and nocturnally, respectively
n_afternoon <- c(3, 7, 2, 8)
# Format data for Fisher's exact test
n_colonies <- matrix(n_afternoon,nrow=2) 
# Fisher's exact test
fisher.test(n_colonies)


###################################################
## 5. Did ant foragers acclimate to experimental manipulations of temperature after 3 hours? ----
###################################################
## 5.1 Import & clean data ----
# Load data
acclimation_data <- read.csv("AcclimationTest.csv", header = T) %>%
  filter(acclimation.time == "3" & alive == "1")  # exclude when ants acclimated to temperatures for 6 and 12 hours (leaving only 3 hours), and remove ants that died before CTmax measurements occurred

# Specify the structure of the variables
acclimation_data$CT.max <- as.numeric(as.character(acclimation_data$CT.max))
acclimation_data$temp.treatment <- as.factor(acclimation_data$temp.treatment)
acclimation_data$ant.colony.id <- as.factor(acclimation_data$ant.colony.id)

## 5.2 Statistical analysis ----
# General linear model testing whether CTmax depended on acclimation temperature, ant source colony, and their interaction
acclimation_model <- glm(CT.max ~ temp.treatment*ant.colony.id, data = acclimation_data, family = quasipoisson(link = "log"))
# Significance test
Anova(acclimation_model, test.statistic = "F")

# Generalized linear model, with the interaction removed to test for the significance of the main effects
acclimation_model <- glm(CT.max ~ temp.treatment + ant.colony.id, data = acclimation_data, family = quasipoisson(link = "log"))
# Significance test
Anova(acclimation_model, test.statistic = "F")
