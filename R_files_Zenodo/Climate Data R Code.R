###################################################
## Nelson et al. 2019, Journal of Animal Ecology R Code
## Purpose: To examine the relationship between elevation and multiple climatic variables and to use
##          principal component analysis to create a linear combination of climatic variables for each site
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
# 2. Are climatic variables correlated with elevation?
#    2.1 Import data
#    2.2 Statistical analysis
# 3. Principal component analysis (PCA) of climatic variables
#    3.1 Import and clean data
#    3.2 Statistical analysis
# 4. How are climatic variables correlated with PC1?
#    4.1 Import and clean data
#    4.2 Statistical analysis

###################################################
## 1. Preliminaries ----
###################################################
## 1.1 Load required packages ----
library(car)
library(dplyr)

## 1.2 Set working directory ----
setwd()

###################################################
## 2. Are climatic variables correlated with elevation? ----
###################################################
## 2.1 Import data ----
clim_data <- read.csv("ClimData.csv")

## 2.2 Statistical analysis ----
# Linear model testing whether the average mean summer temperature is correlated with elevation
mean_model <- lm(summer.temp ~ elevation.m, data = clim_data)
Anova(mean_model, test.statistic = "F")
summary(mean_model)

# Linear model testing whether the average maximum summer temperature is correlated with elevation
max_model <- lm(max.sum.temp ~ elevation.m, data = clim_data)
Anova(max_model, test.statistic = "F")
summary(max_model)

# Linear model testing whether the average minimum summer temperature is correlated with elevation
min_model <- lm(min.sum.temp ~ elevation.m, data = clim_data)
Anova(min_model, test.statistic = "F")
summary(min_model)

# Linear model testing whether the average range in summer temperatures is correlated with elevation
range_model <- lm(sum.temp.range ~ elevation.m, data = clim_data)
Anova(range_model, test.statistic = "F")
summary(range_model)

# Linear model testing whether mean summer precipitation is correlated with elevation
precip_model <- lm(summer.precip ~ elevation.m, data = clim_data)
Anova(precip_model, test.statistic = "F")
summary(precip_model)

###################################################
## 3. Principal component analysis (PCA) of climatic variables ----
###################################################
## 3.1 Import and clean data ----
# Clean data so that each independent estimate of the climatic variables is represented once
cleaned_data <- clim_data %>% group_by(summer.temp) %>%
  summarize(summer.precip = mean(summer.precip), 
            max.sum.temp = mean(max.sum.temp),
            min.sum.temp = mean(min.sum.temp),
            sum.temp.range = mean(sum.temp.range))

## 3.2 Statistical analysis ----
# Conduct the PCA
pca <- prcomp(cleaned_data, center = T, scale. = T)
summary(pca) # Summarize the PCA, get the proportion of variance explained by each component
pca$x # Get the principal component values for each site

###################################################
## 4. How are climatic variables correlated with PC1? ----
###################################################
## 4.1 Import and clean data ----
PC1_data <- clim_data %>% group_by(summer.temp) %>%
  summarize(summer.precip = mean(summer.precip), 
            max.sum.temp = mean(max.sum.temp),
            min.sum.temp = mean(min.sum.temp),
            sum.temp.range = mean(sum.temp.range), 
            PC1 = mean(PC1))

## 4.2 Statistical analysis ----
# Linear model testing whether average mean summer temperature is correlated with PC1
mean_model <- lm(summer.temp ~ PC1, data = PC1_data)
Anova(mean_model, test.statistic = "F")

# Linear model testing whether average maximum summer temperature is correlated with PC1
max_model <- lm(max.sum.temp ~ PC1, data = PC1_data)
Anova(max_model, test.statistic = "F")

# Linear model testing whether average minimum summer temperature is correlated with PC1
min_model <- lm(min.sum.temp ~ PC1, data = PC1_data)
Anova(min_model, test.statistic = "F")

# Linear model testing whether average range in summer temperatures is correlated with PC1
range_model <- lm(sum.temp.range ~ PC1, data = PC1_data)
Anova(range_model, test.statistic = "F")

# Linear model testing whether mean summer precipitation is correlated with PC1
precip_model <- lm(summer.precip ~ PC1, data = PC1_data)
Anova(precip_model, test.statistic = "F")
