############################################################################
# Project: Assort_Mating_Dunnocks
# Script name: Analyses_behavioural_effects_on_mating_system
# Script purpose: Test whether individual behaviour predicts the mating system
# Author: Benedikt Holtmann
# Date Created: 25/08/2020
# R version 4.0.2
############################################################################

rm(list = ls())

# Load packages ------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(MCMCglmm)

# Set working directory -----------------------------------------------------

setwd("~ set the working directory")

# Read data ----------------------------------------------------------------


FID_avg_year_first_broods <- read.csv("03_mating_system_FID_first_broods.csv", na.strings = c("", "NA"))


# Analyses for individuals using first broods ----------------------------------------------------------

### 1. Male model  ======================

# Make male dataset
FID_avg_year_first_broods_m <- FID_avg_year_first_broods %>% filter(Sex == "Male")
FID_avg_year_first_broods_m <- droplevels(FID_avg_year_first_broods_m)

# Remove NAs in fixed effects, necessary for MCMCglmm!
# FID_avg_year_first_broods_m <- FID_avg_year_first_broods_m %>% drop_na(c(mean_Hops_min, mean_FID, mean_Pecks_min, Age, SeptDay))
FID_avg_year_first_broods_m <- FID_avg_year_first_broods_m %>% drop_na(mean_FID)

FID_avg_year_first_broods_m <- droplevels(FID_avg_year_first_broods_m)


# Number of males
FID_avg_year_first_broods_m %>%
  distinct(BirdID) %>%
  tally()
# 57

# Observations per Season
FID_avg_year_first_broods_m %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    42
# 2           2013    49
# 3           2014    59

# Get the number of males for each Season
FID_avg_year_first_broods_m %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    30
# 2           2013    37
# 3           2014    38

# Repeated samples per season
FID_avg_year_first_broods_m %>%
  count(BreedingSeason, BirdID, name = "repeats") %>%
  count(BreedingSeason, repeats)


# Mean measurments per male
mean(tapply(FID_avg_year_first_broods_m$BirdID, FID_avg_year_first_broods_m$BirdID, length))
# [1] 2.631579

hist(scale(FID_avg_year_first_broods_m$mean_FID))

FID_avg_year_first_broods_m$MatingSystem_2 <- as.factor(FID_avg_year_first_broods_m$MatingSystem_2)


# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior_bin <- list(R = list(V = 1, fix = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

# MCMCglmm
MCMCglmm_assoc_FID_mating_first_broods_males <- MCMCglmm(MatingSystem_2 ~ 1 +
  scale(mean_FID) +
  scale(Age) +
  ObsPeriod,
random = ~BirdID, prior = prior_bin,
data = FID_avg_year_first_broods_m, singular.ok = TRUE,
nitt = 13000 * 100, burnin = 3000 * 100, thin = 10 * 50,
family = "categorical", pl = TRUE, saveX = TRUE
)

save(MCMCglmm_assoc_FID_mating_first_broods_males, file="~/.../.../MCMCglmm_assoc_FID_mating_first_broods_males.Rdata")

# Model diagnostics
plot(MCMCglmm_assoc_FID_mating_first_broods_males$VCV)
plot(MCMCglmm_assoc_FID_mating_first_broods_males$Sol)
autocorr(MCMCglmm_assoc_FID_mating_first_broods_males$VCV)
autocorr(MCMCglmm_assoc_FID_mating_first_broods_males$Sol)

summary(MCMCglmm_assoc_FID_mating_first_broods_males)


##### 1.2 Extract results male model  ======================
Intercept <- mean(MCMCglmm_assoc_FID_mating_first_broods_males$Sol[, 1])
CI_Intercept <- HPDinterval(MCMCglmm_assoc_FID_mating_first_broods_males$Sol[, 1])

Coef_z_mean_FID <- mean(MCMCglmm_assoc_FID_mating_first_broods_males$Sol[, 2])
CI_z_mean_FID <- HPDinterval(MCMCglmm_assoc_FID_mating_first_broods_males$Sol[, 2])

Coef_z_Age <- mean(MCMCglmm_assoc_FID_mating_first_broods_males$Sol[, 3])
CI_z_Age <- HPDinterval(MCMCglmm_assoc_FID_mating_first_broods_males$Sol[, 3])

Coef_BreedingSeason <- mean(MCMCglmm_assoc_FID_mating_first_broods_males$Sol[, 4])
CI_BreedingSeason <- HPDinterval(MCMCglmm_assoc_FID_mating_first_broods_males$Sol[, 4])

# Combine everything in one table

PosteriorMean <- rbind(Intercept, Coef_z_mean_FID, Coef_z_Age, Coef_BreedingSeason)
CredibleIntervals <- rbind(CI_Intercept, CI_z_mean_FID, CI_z_Age, CI_BreedingSeason)
Predictors <- rbind("Intercept", "z_mean_FID", "z_Age", "BreedingSeason")

results <- data.frame(cbind(Predictors, PosteriorMean, CredibleIntervals))
colnames(results) <- c("Predictors", "PosteriorMean", "CI_low", "CI_high")

# Save table
write.table(results, file = "Assoc_FID_mating_first_broods_males_results.csv", sep = ",", row.names = FALSE) 


### 2. Female model  ======================

# Make female dataset
FID_avg_year_first_broods_f <- FID_avg_year_first_broods %>% filter(Sex == "Female")
FID_avg_year_first_broods_f <- droplevels(FID_avg_year_first_broods_f)

# Remove NAs in fixed effects, necessary for MCMCglmm models
FID_avg_year_first_broods_f <- FID_avg_year_first_broods_f %>% drop_na(mean_FID)

FID_avg_year_first_broods_f <- droplevels(FID_avg_year_first_broods_f)


# Number of females
FID_avg_year_first_broods_f %>%
  distinct(BirdID) %>%
  tally()
# 41

# Observations per season
FID_avg_year_first_broods_f %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    20
# 2           2013    29
# 3           2014    35

# Get the number of males for each season
FID_avg_year_first_broods_f %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    17
# 2           2013    25
# 3           2014    29

# Mean measurments per female
mean(tapply(FID_avg_year_first_broods_f$BirdID, FID_avg_year_first_broods_f$BirdID, length))
# [1] 2.04878

# Repeated samples per season
FID_avg_year_first_broods_f %>%
  count(BreedingSeason, BirdID, name = "repeats") %>%
  count(BreedingSeason, repeats)

hist(scale(FID_avg_year_first_broods_f$mean_FID))

FID_avg_year_first_broods_f$MatingSystem_2 <- as.factor(FID_avg_year_first_broods_f$MatingSystem_2)


# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior_bin <- list(R = list(V = 1, fix = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

# MCMCglmm
MCMCglmm_assoc_FID_mating_first_broods_females <- MCMCglmm(MatingSystem_2 ~ 1 +
  scale(mean_FID) +
  scale(Age) +
  ObsPeriod,
random = ~BirdID, prior = prior_bin,
data = FID_avg_year_first_broods_f, singular.ok = TRUE,
nitt = 13000 * 200, burnin = 3000 * 200, thin = 10 * 100,
family = "categorical", pl = TRUE, saveX = TRUE
)

save(MCMCglmm_assoc_FID_mating_first_broods_females, file="~/.../.../MCMCglmm_assoc_FID_mating_first_broods_females.Rdata")

# Model diagnostics
plot(MCMCglmm_assoc_FID_mating_first_broods_females$VCV)
plot(MCMCglmm_assoc_FID_mating_first_broods_females$Sol)
autocorr(MCMCglmm_assoc_FID_mating_first_broods_females$VCV)
autocorr(MCMCglmm_assoc_FID_mating_first_broods_females$Sol)

summary(MCMCglmm_assoc_FID_mating_first_broods_females)


##### 2.2 Extract results female model  ======================
Intercept <- mean(MCMCglmm_assoc_FID_mating_first_broods_females$Sol[, 1])
CI_Intercept <- HPDinterval(MCMCglmm_assoc_FID_mating_first_broods_females$Sol[, 1])

Coef_z_mean_FID <- mean(MCMCglmm_assoc_FID_mating_first_broods_females$Sol[, 2])
CI_z_mean_FID <- HPDinterval(MCMCglmm_assoc_FID_mating_first_broods_females$Sol[, 2])

Coef_z_Age <- mean(MCMCglmm_assoc_FID_mating_first_broods_females$Sol[, 3])
CI_z_Age <- HPDinterval(MCMCglmm_assoc_FID_mating_first_broods_females$Sol[, 3])

Coef_BreedingSeason <- mean(MCMCglmm_assoc_FID_mating_first_broods_females$Sol[, 4])
CI_BreedingSeason <- HPDinterval(MCMCglmm_assoc_FID_mating_first_broods_females$Sol[, 4])

# Combine everything in one table
PosteriorMean <- rbind(Intercept, Coef_z_mean_FID, Coef_z_Age, Coef_BreedingSeason)
CredibleIntervals <- rbind(CI_Intercept, CI_z_mean_FID, CI_z_Age, CI_BreedingSeason)
Predictors <- rbind("Intercept", "z_mean_FID", "z_Age", "BreedingSeason")

results <- data.frame(cbind(Predictors, PosteriorMean, CredibleIntervals))
colnames(results) <- c("Predictors", "PosteriorMean", "CI_low", "CI_high")

# Save table
write.table(results, file = "Assoc_FID_mating_first_broods_females_results.csv", sep = ",", row.names = FALSE) 



### 3. Model testing whether FID of a male predicts its status in trios  ======================

# Get only polygamous groups
FID_avg_year_first_broods_m_poly <- FID_avg_year_first_broods_m %>% 
  filter(MatingSystem_2 != "Monogamous")

FID_avg_year_first_broods_m_poly <- droplevels(FID_avg_year_first_broods_m_poly)

# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior_bin = list(R = list(V = 1, fix = 1, nu=0.002), G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V= 1000)))

# MCMCglmm
MCMCglmm_assoc_status_mating_first_broods <- MCMCglmm(Status ~ 1 +
  scale(mean_FID) +
  scale(Age) +
  ObsPeriod,
random = ~BirdID, prior = prior_bin,
data = FID_avg_year_first_broods_m_poly, singular.ok = TRUE,
nitt = 13000 * 100, burnin = 3000 * 100, thin = 10 * 50,
family = "categorical", pl = TRUE, saveX = TRUE
)

save(MCMCglmm_assoc_status_mating_first_broods, file="~/.../.../MCMCglmm_assoc_status_mating_first_broods.Rdata")

# Model diagnostics
plot(MCMCglmm_assoc_status_mating_first_broods$VCV)
plot(MCMCglmm_assoc_status_mating_first_broods$Sol)
autocorr(MCMCglmm_assoc_status_mating_first_broods$VCV)
autocorr(MCMCglmm_assoc_status_mating_first_broods$Sol)

summary(MCMCglmm_assoc_status_mating_first_broods)


##### 3.2 Extract results male model  ======================
Intercept <- mean(MCMCglmm_assoc_status_mating_first_broods$Sol[, 1])
CI_Intercept <- HPDinterval(MCMCglmm_assoc_status_mating_first_broods$Sol[, 1])

Coef_z_mean_FID <- mean(MCMCglmm_assoc_status_mating_first_broods$Sol[, 2])
CI_z_mean_FID <- HPDinterval(MCMCglmm_assoc_status_mating_first_broods$Sol[, 2])

Coef_z_Age <- mean(MCMCglmm_assoc_status_mating_first_broods$Sol[, 3])
CI_z_Age <- HPDinterval(MCMCglmm_assoc_status_mating_first_broods$Sol[, 3])

Coef_BreedingSeason <- mean(MCMCglmm_assoc_status_mating_first_broods$Sol[, 4])
CI_BreedingSeason <- HPDinterval(MCMCglmm_assoc_status_mating_first_broods$Sol[, 4])

# Combine everything in one table
PosteriorMean <- rbind(Intercept, Coef_z_mean_FID, Coef_z_Age, Coef_BreedingSeason)
CredibleIntervals <- rbind(CI_Intercept, CI_z_mean_FID, CI_z_Age, CI_BreedingSeason)
Predictors <- rbind("Intercept", "z_mean_FID", "z_Age", "BreedingSeason")

results <- data.frame(cbind(Predictors, PosteriorMean, CredibleIntervals))
colnames(results) <- c("Predictors", "PosteriorMean", "CI_low", "CI_high")

# Save table
write.table(results, file = "Assoc_status_mating_first_broods_males_results.csv", sep = ",", row.names = FALSE) 
