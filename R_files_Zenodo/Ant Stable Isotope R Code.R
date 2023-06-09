###################################################
## Nelson et al. 2018 Oikos R Code
## Purpose: To assess whether ant diet as indicated by stable isotope analysis varied with elevation
##          or mean summer temperature
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
#    1.3 Import data
# 2. Does ant delta-N-15 depend on elevation or mean summer temperature?
# 3. Does ant C:N ratio depend on elevation or mean summer temperature?
# 4. Does ant delta-C-13 depend on elevation or mean summer temperature?
# 5. Does ant percent N depend on elevation or mean summer temperature?
# 6. Does ant percent C depend on elevation or mean summer temperature?

###################################################
## 1. Preliminaries ----
###################################################
## 1.1 Load required packages ----
library(lme4)
library(car)

## 1.2 Set working directory ----
setwd()

## 1.3 Import data
isotope_data <- read.csv("AntStableIsotopes.csv", header = T)

###################################################
## 2. Does ant delta-N-15 depend on elevation or mean summer temperature? ----
###################################################
# Linear mixed model testing whether ant delta-N-15 depended on elevation (high vs. low) and ant life stage
dN15_model <- lmer(dN15 ~ elevation + life.stage + (1|elevation:valley) + (1|valley:mound.id), data = isotope_data)
Anova(dN15_model, type = 3)

# Linear mixed model testing whether ant delta-N-15 depended on mean summer temperature and ant life stage
dN15_model_temp <- lmer(dN15 ~ summer.temp + life.stage + (1|valley) + (1|valley:mound.id), data = isotope_data)
Anova(dN15_model_temp, type = 3)

###################################################
## 3. Does ant C:N ratio depend on elevation or mean summer temperature? ----
###################################################
# Linear mixed model testing whether ant C:N ratio depended on elevation (high vs. low) and ant life stage
C.N_model <- lmer(C.N ~ elevation + life.stage + (1|elevation:valley) + (1|valley:mound.id), data = isotope_data)
Anova(C.N_model, type = 3)

# Linear mixed model testing whether ant C:N ratio depended on mean summer temperature and ant life stage
C.N_model_temp <- lmer(C.N ~ summer.temp + life.stage + (1|valley) + (1|valley:mound.id), data = isotope_data)
Anova(C.N_model_temp, type = 3)

###################################################
## 4. Does ant delta-C-13 depend on elevation or mean summer temperature? ----
###################################################
# Linear mixed model testing whether ant delta-C-13 depended on elevation (high vs. low) and ant life stage
dC13_model <- lmer(dC13 ~ elevation + life.stage + (1|elevation:valley) + (1|valley:mound.id), data = isotope_data)
Anova(dC13_model, type = 3)

# Linear mixed model testing whether ant delta-C-13 depended on mean summer temperature and ant life stage
dC13_model_temp <- lmer(dC13 ~ summer.temp + life.stage + (1|valley) + (1|valley:mound.id), data = isotope_data)
Anova(dC13_model_temp, type = 3)

###################################################
## 5. Does ant percent N depend on elevation or mean summer temperature? ----
###################################################
# Linear mixed model testing whether ant percent N depended on elevation (high vs. low) and ant life stage
perc_N_model <- lmer(perc.N ~ elevation + life.stage + (1|elevation:valley) + (1|valley:mound.id), data = isotope_data)
Anova(perc_N_model, type = 3)

# Linear mixed model testing whether ant percent N depended on mean summer temperature and ant life stage
perc_N_model_temp <- lmer(perc.N ~ summer.temp + life.stage + (1|valley) + (1|valley:mound.id), data = isotope_data)
Anova(perc_N_model_temp, type = 3)

###################################################
## 6. Does ant percent C depend on elevation or mean summer temperature? ----
###################################################
# Linear mixed model testing whether ant percent C depended on elevation (high vs. low) and ant life stage
perc_C_model <- lmer(perc.C ~ elevation + life.stage + (1|elevation:valley) + (1|valley:mound.id), data = isotope_data)
Anova(perc_C_model, type = 3)

# Linear mixed model testing whether ant percent C depended on mean summer temperature and ant life stage
perc_C_model_temp <- lmer(perc.C ~ summer.temp + life.stage + (1|elevation:valley) + (1|valley:mound.id), data = isotope_data)
Anova(perc_C_model_temp, type = 3)
