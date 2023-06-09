###################################################
## Nelson et al. 2019, Journal of Animal Ecology R Code
## Purpose: To test for differences in ant abundance and taxonomic composition in pitfall traps as well as
##          bait consumption along an elevational gradient in aridity
## Corresponding author: Annika S. Nelson (University of California, Irvine, annika.nelson@uci.edu)
###################################################
## README ----
# This script is divided into the following sections:
# 
# 1. Preliminaries
#    1.1 Load required packages
#    1.2 Set working directory
# 2. Did ant abundance in pitfall traps depend on aridity?
#    2.1 Import and clean data
#    2.2 Statistical analysis
# 3. Did dissimilarity in ant taxonomic composition in pitfall traps increase with differences in aridity?
#    3.1 Import and clean data
#    3.2 Statistical analysis
# 4. Did ant discovery of baits depend on aridity?
#    4.1 Import and clean data
#    4.2 Statistical analysis
# 5. Did ant consumption of baits depend on aridity?
#    5.1 Import and clean data
#    5.2 Statistical analysis

###################################################
## 1. Preliminaries ----
###################################################
## 1.1 Load required packages ----
library(lme4)
library(car)
library(glmmADMB)
library(vegan)
library(mixtools)
library(dplyr)

## 1.2 Set working directory ----
setwd()

###################################################
## 2. Did ant abundance in pitfall traps depend on aridity? ----
###################################################
## 2.1 Import and clean data ----
# Import data, creating a variable indicating the number of ants collected per trap per day of sampling in each plot
abund_data <- read.csv("AntAbundanceResourceUse.csv", header = T) %>%
  mutate(ants.per.trapday = num.ants/trap.days)

# Create separate data frames for 2012 and 2015 data
abund_2012_data <- abund_data %>% filter(year == "2012")
abund_2015_data <- abund_data %>% filter(year == "2015")

## 2.2 Statistical analysis ----
# Linear mixed model testing whether the number of ants collected per trap per day in each site in 2012 depended on aridity
abund_2012_model <- lmer(log(ants.per.trapday) ~ PC1 + (1|valley), data = abund_2012_data)
# Significance test
Anova(abund_2012_model, test.statistic = "F")

# Linear mixed model testing whether the number of ants collected per trap per day in each site in 2015 depended on PC1
abund_2015_model <- lmer(log(ants.per.trapday) ~ PC1 + (1|valley), data = abund_2015_data)
# Significance test
Anova(abund_2015_model, test.statistic = "F")

###################################################
## 3. Did dissimilarity in ant taxonomic composition in pitfall traps increase with differences in aridity? ----
###################################################
## 3.1 Import and clean data ----
# Import 2012 multivariate data
multivar_data_2012 <- read.csv("2012MultivariateData.csv", header = T, row.names = 1)

# Import 2015 multivariate data
multivar_data_2015 <- read.csv("2015MultivariateData.csv", header = T, row.names = 1)

# Create distance matrices for PC1 values in 2012 and 2015
clim_dist_2012 <- as.matrix(dist(multivar_data_2012$PC1))
clim_dist_2015 <- as.matrix(dist(multivar_data_2015$PC1))

# Create distance matrices describing differences in ant community composition among sites in 2012 and 2015
ant_dist_2012 <- as.matrix(dist(cbind(multivar_data_2012$rufa.per.trapday, multivar_data_2012$podz.per.trapday, multivar_data_2012$tap.per.trapday, multivar_data_2012$camp.per.trapday, multivar_data_2012$myr.per.trapday)))

ant_dist_2015 <- as.matrix(dist(cbind(multivar_data_2015$rufa.per.trapday, multivar_data_2015$podz.per.trapday, multivar_data_2015$tap.per.trapday, multivar_data_2015$camp.per.trapday, multivar_data_2015$myr.per.trapday)))

## 3.2 Statistical analysis ----
# Mantel test for 2012
mantel(clim_dist_2012, ant_dist_2012, permutations = 999)

# Mantel test for 2015
mantel(clim_dist_2015, ant_dist_2015, permutations = 999)

###################################################
## 4. Did ant discovery of baits depend on aridity? ----
###################################################
## 4.1 Import and clean data ----
# Import data, removing times when there were no baits in the field
bait_data <- abund_data %>% filter(bait.consumption != "NA")

# Define variables
bait_data$bait.consumption <- as.numeric(as.character(bait_data$bait.consumption))

# Create a mixture model to classify baits as discovered (vs. not discovered) based on bait consumption rates
mixture_model <- normalmixEM(bait_data$bait.consumption, mean.constr = c(0, NA))
# Summary of the two separate distributions
summary(mixture_model)
# Determine which distribution the points fall into based on 50% posterior probabilities
posterior_prob <- as.data.frame(cbind(bait.consumption = mixture_model$x, mixture_model$posterior)) 
posterior_prob <- posterior_prob %>% mutate(distribution = ifelse(comp.1 > 0.5, "0", "1"))

# Add factor to the data frame describing whether baits were discovered (discovered = 1) or not (discovered = 0) based on mixture model results
bait_data <- bait_data %>% mutate(discovered = as.factor(posterior_prob$distribution))

## 4.2 Statistical analysis ----
# Generalized linear mixed model testing whether bait discovery depended on aridity
disc_model <- glmmadmb(discovered ~ PC1 + (1|valley), family = "binomial", data = bait_data)

# Significance test
summary(disc_model)$coefficients

###################################################
## 5. Did ant consumption of baits depend on aridity? ----
###################################################
## 5.1 Import and clean data ----
# Create data frame that only includes baits that were discovered
consumption_data <- bait_data %>% filter(discovered == "1")

## 5.2 Statistical analysis ----
# Linear mixed model testing whether ant bait consumption rates depended on aridity
consumption_model <- lmer(bait.consumption ~ PC1 + (1|valley), data = consumption_data)

# Significance test
Anova(consumption_model, test.statistic = "F")
