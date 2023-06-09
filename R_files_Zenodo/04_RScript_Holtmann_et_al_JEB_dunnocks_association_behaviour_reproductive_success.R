############################################################################
# Project: Assort_Mating_Dunnocks
# Script name: Analyses_behavioural_effects_on_rep_success
# Script purpose: Test whether behaviours predict the mating system
# Author: Benedikt Holtmann
# Date Created: 25/08/2020
# R version 4.0.2
############################################################################

rm(list = ls())

# Load packages ------------------------------------------------------------

library(tidyverse)
library(MCMCglmm)
library(plotrix) # To use std.error() function


# Set working directory -----------------------------------------------------

setwd("~ set the working directory")


# 1. Association between FID and reproductive success ----------------------------------------------------------

# Read data ----------------------------------------------------------------

FID_avg_year_rep_success_males <- read.csv("0.4.1_assoc_rep_success_FID__males.csv", na.strings = c("", "NA"))
FID_avg_year_rep_success_females <- read.csv("0.4.2_assoc_rep_success_FID_females.csv", na.strings = c("", "NA"))

### 1.1 Male model ======================

# Change status so that monogamous males are the intercept
FID_avg_year_rep_success_males <- FID_avg_year_rep_success_males %>% mutate(Status_Male = ifelse(FID_avg_year_rep_success_males$Status_Male == "mono", "alpha_mono",
  ifelse(FID_avg_year_rep_success_males$Status_Male == "alpha", "alpha_poly", "beta_poly"))
)

# Number of males
FID_avg_year_rep_success_males %>%
  distinct(MaleID) %>%
  tally()
# 53

# Observations per season
FID_avg_year_rep_success_males %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    32
# 2           2013    45
# 3           2014    63

# Get the number of males for each season
FID_avg_year_rep_success_males %>%
  group_by(BreedingSeason) %>%
  distinct(MaleID) %>%
  tally()
# 1           2012    22
# 2           2013    33
# 3           2014    37

# Mean measurments per male
mean(tapply(FID_avg_year_rep_success_males$MaleID, FID_avg_year_rep_success_males$MaleID, length))
# [1] 2.641509

hist(FID_avg_year_rep_success_males$Fledged_Male)

mean(FID_avg_year_rep_success_males$Fledged_Male) # 0.4214286
sd(FID_avg_year_rep_success_males$Fledged_Male)
se<- sd(FID_avg_year_rep_success_males$Fledged_Male)/sqrt(length(FID_avg_year_rep_success_males$Fledged_Male)) # 0.06729837
std.error(FID_avg_year_rep_success_males$Fledged_Male)


# MCMCglmm ----------------------------------------------------------------
# Specify priors
prior_uni<- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

MCMCglmm_assoc_FID_RepSuccess_males_poisson <- MCMCglmm(Fledged_Male ~ 1 +
  scale(mean_FID_Male) * scale(mean_FID_Female) +
  Status_Male + 
  ObsPeriod,
random = ~MaleID, prior = prior_uni,
data = FID_avg_year_rep_success_males, singular.ok = TRUE,
nitt = 13000 * 100,  burnin = 3000 * 100, thin = 10 * 50,
family = "poisson", pl = TRUE, saveX = TRUE
)

save(MCMCglmm_assoc_FID_RepSuccess_males_poisson, file="~/.../.../MCMCglmm_assoc_FID_RepSuccess_males_poisson.Rdata")

# Model diagnostics
plot(MCMCglmm_assoc_FID_RepSuccess_males_poisson$VCV)
plot(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol)
autocorr(MCMCglmm_assoc_FID_RepSuccess_males_poisson$VCV)
autocorr(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol)

summary(MCMCglmm_assoc_FID_RepSuccess_males_poisson)


##### 1.2 Extract results male model  ======================
Intercept <- mean(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 1])
CI_Intercept <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 1])

Coef_z_mean_FID_Male <- mean(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 2])
CI_z_mean_FID_Male <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 2])

Coef_z_mean_FID_Female <- mean(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 3])
CI_z_mean_FID_Female <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 3])

Coef_Polygamous_alpha <- mean(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 4])
CI_Polygamous_alpha <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 4])

Coef_Polygamous_beta <- mean(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 5])
CI_Polygamous_beta <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 5])

Coef_BreedingSeason <- mean(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 6])
CI_BreedingSeason <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 6])

Coef_Interaction <- mean(MCMC_univariate_FID_repSuccess_males_groups_poisson$Sol[, 7])
CI_Coef_Interaction <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_males_poisson$Sol[, 7])

# Combine everything in one table
PosteriorMean <- rbind(Intercept, Coef_z_mean_FID_Male, Coef_z_mean_FID_Female, Coef_Polygamous_alpha, Coef_Polygamous_beta, Coef_BreedingSeason, Coef_Interaction)
CredibleIntervals <- rbind(CI_Intercept, CI_z_mean_FID_Male, CI_z_mean_FID_Female, CI_Polygamous_alpha, CI_Polygamous_beta, CI_BreedingSeason, CI_Coef_Interaction)
Predictors <- rbind("Intercept", "z_mean_FID_Male", "z_mean_FID_Female", "Polygamous_alpha", "Polygamous_beta",  "BreedingSeason", "Interaction")

results <- data.frame(cbind(Predictors, PosteriorMean, CredibleIntervals))
colnames(results) <- c("Predictors", "PosteriorMean", "CI_low", "CI_high")

# Save table
write.table(results, file = "MCMCglmm_assoc_FID_RepSuccess_males_poisson_results.csv", sep = ",", row.names = FALSE)


### 1.3 Female model ======================

# Number of females
FID_avg_year_rep_success_females %>%
  distinct(FemaleID) %>%
  tally()
# 40

# Observations per season
FID_avg_year_rep_success_females %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    21
# 2           2013    31
# 3           2014    40

# Get the number of males for each season
FID_avg_year_rep_success_females %>%
  group_by(BreedingSeason) %>%
  distinct(FemaleID) %>%
  tally()
# 1           2012    17
# 2           2013    24
# 3           2014    29

# Mean measurments per male
mean(tapply(FID_avg_year_rep_success_females$FemaleID, FID_avg_year_rep_success_females$FemaleID, length))
# [1] 2.3

hist(FID_avg_year_rep_success_females$Fledged_Female)

mean(FID_avg_year_rep_success_males$Fledged_Female) # 0.8285714
std.error(FID_avg_year_rep_success_males$Fledged_Female) # 0.1008407


# MCMCglmm ----------------------------------------------------------------
# Specify prior
prior_uni<- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

MCMCglmm_assoc_FID_RepSuccess_females_poisson <- MCMCglmm(Fledged_Female ~ 1 +
  scale(mean_FID_Female) * scale(mean_FID_Males) +
  Status_Female +
  ObsPeriod,
random = ~FemaleID, prior = prior_uni,
data = FID_avg_year_rep_success_females, singular.ok = TRUE,
nitt = 13000 * 300, burnin = 3000 * 300, thin = 10 * 150,
family = "poisson", pl = TRUE, saveX = TRUE
)

save(MCMCglmm_assoc_FID_RepSuccess_females_poisson, file="~/.../.../MCMCglmm_assoc_FID_RepSuccess_females_poissonn.Rdata")

# Model diagnostics
plot(MCMCglmm_assoc_FID_RepSuccess_females_poisson$VCV)
plot(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol)
autocorr(MCMCglmm_assoc_FID_RepSuccess_females_poisson$VCV)
autocorr(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol)

summary(MCMCglmm_assoc_FID_RepSuccess_females_poisson)


##### 1.4 Extract results male model  ======================
Intercept <- mean(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 1])
CI_Intercept <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 1])

Coef_z_mean_FID_Female <- mean(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 2])
CI_z_mean_FID_Female <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 2])

Coef_z_mean_FID_Male <- mean(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 3])
CI_z_mean_FID_Male <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 3])

Coef_Polygamous<- mean(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 4])
CI_Polygamous <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 4])

Coef_BreedingSeason <- mean(MCMC_univariate_FID_repSuccess_females_groups_poisson$Sol[, 5])
CI_BreedingSeason <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 5])

Coef_Interaction <- mean(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 6])
CI_Coef_Interaction <- HPDinterval(MCMCglmm_assoc_FID_RepSuccess_females_poisson$Sol[, 6])

# Combine everything in one table
PosteriorMean <- rbind(Intercept, Coef_z_mean_FID_Female, Coef_z_mean_FID_Male, Coef_Polygamous, Coef_BreedingSeason, Coef_Interaction)
CredibleIntervals <- rbind(CI_Intercept, CI_z_mean_FID_Female, CI_z_mean_FID_Male, CI_Polygamous, CI_BreedingSeason, CI_Coef_Interaction)
Predictors <- rbind("Intercept", "z_mean_FID_Female", "z_mean_FID_Male", "Polygamous", "BreedingSeason", "Interaction")

results <- data.frame(cbind(Predictors, PosteriorMean, CredibleIntervals))
colnames(results) <- c("Predictors", "PosteriorMean", "CI_low", "CI_high")

# Save table
write.table(results, file = "MCMCglmm_assoc_FID_RepSuccess_females_poisson_results.csv", sep = ",", row.names = FALSE) # save results in one table



# 2. Association between provisioning and reproductive success ----------------------------------------------------------

# Read data ----------------------------------------------------------------

provisioning_avg_year_rep_success_males <- read.csv("0.4.3_assoc_rep_success_provisioning_males.csv", na.strings = c("", "NA"))
provisioning_avg_year_rep_success_females <- read.csv("0.4.4_assoc_rep_success_provisioning_females.csv", na.strings = c("", "NA"))

### 2.1 Male model ======================

# Change status so that monogamous males are the intercept
provisioning_avg_year_rep_success_males <- provisioning_avg_year_rep_success_males %>% mutate(Status_Male = ifelse(provisioning_avg_year_rep_success_males$Status_Male == "mono", "alpha_mono",
  ifelse(provisioning_avg_year_rep_success_males$Status_Male == "alpha", "alpha_poly", "beta_poly")
))

# Number of males
provisioning_avg_year_rep_success_males %>%
  distinct(MaleID) %>%
  tally()
# 52

# Observations per season
provisioning_avg_year_rep_success_males %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    25
# 2           2013    38
# 3           2014    39

# Get the number of males for each season
provisioning_avg_year_rep_success_males %>%
  group_by(BreedingSeason) %>%
  distinct(MaleID) %>%
  tally()
# 1           2012    22
# 2           2013    33
# 3           2014    31

# Mean measurments per male
mean(tapply(provisioning_avg_year_rep_success_males$MaleID, provisioning_avg_year_rep_success_males$MaleID, length))
# [1] 1.961538

hist(provisioning_avg_year_rep_success_males$Fledged_Male)

mean(provisioning_avg_year_rep_success_males$Fledged_Male) # 0.8235294
std.error(provisioning_avg_year_rep_success_males$Fledged_Male) # 0.08747162


# MCMCglmm ----------------------------------------------------------------
# Specify priors
prior_uni<- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

MCMCglmm_assoc_provisioning_RepSuccess_males_poisson <- MCMCglmm(Fledged_Male ~ 1 +
  scale(mean_Visits_h_Male) * scale(mean_Visits_h_Female) +
  Status_Male +
  ObsPeriod,
random = ~MaleID, prior = prior_uni,
data = provisioning_avg_year_rep_success_males, singular.ok = TRUE,
nitt = 13000 * 200, burnin = 3000 * 200, thin = 10 * 100,
family = "poisson", pl = TRUE, saveX = TRUE
)

save(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson, file="~/.../.../MCMCglmm_assoc_provisioning_RepSuccess_males_poisson.Rdata")

# Model diagnostics
plot(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$VCV)
plot(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol)
autocorr(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$VCV)
autocorr(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol)

summary(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson)


##### 2.2 Extract results male model  ======================
Intercept <- mean(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 1])
CI_Intercept <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 1])

Coef_z_mean_Visit_Male <- mean(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 2])
CI_z_mean_Visit_Male <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 2])

Coef_z_mean_Visit_Female <- mean(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 3])
CI_z_mean_Visit_Female <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 3])

Coef_Polygamous_alpha <- mean(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 4])
CI_Polygamous_alpha <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 4])

Coef_Polygamous_beta <- mean(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 5])
CI_Polygamous_beta <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 5])

Coef_BreedingSeason <- mean(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 6])
CI_BreedingSeason <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 6])

Coef_Interaction <- mean(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 7])
CI_Coef_Interaction <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_males_poisson$Sol[, 7])

# Combine everything in one table
PosteriorMean <- rbind(Intercept, Coef_z_mean_Visit_Male, Coef_z_mean_Visit_Female, Coef_Polygamous_alpha, Coef_Polygamous_beta, Coef_BreedingSeason, Coef_Interaction)
CredibleIntervals <- rbind(CI_Intercept, CI_z_mean_Visit_Male, CI_z_mean_Visit_Female, CI_Polygamous_alpha, CI_Polygamous_beta, CI_BreedingSeason, CI_Coef_Interaction)
Predictors <- rbind("Intercept", "z_mean_FID_Male", "z_mean_Visit_Female", "Polygamous_alpha", "Polygamous_beta",  "BreedingSeason", "Interaction")

results <- data.frame(cbind(Predictors, PosteriorMean, CredibleIntervals))
colnames(results) <- c("Predictors", "PosteriorMean", "CI_low", "CI_high")

# Save table
write.table(results, file = "MCMCglmm_assoc_provisioning_RepSuccess_males_poisson_results.csv", sep = ",", row.names = FALSE)


### 2.3 Female model ======================

# Number of females
provisioning_avg_year_rep_success_females %>%
  distinct(FemaleID) %>%
  tally()
# 38

# Observations per season
provisioning_avg_year_rep_success_females %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    16
# 2           2013    26
# 3           2014    25

# Get the number of males for each season
provisioning_avg_year_rep_success_females %>%
  group_by(BreedingSeason) %>%
  distinct(FemaleID) %>%
  tally()
# 1           2012    16
# 2           2013    24
# 3           2014    21

# Mean measurments per female
mean(tapply(provisioning_avg_year_rep_success_females$FemaleID, provisioning_avg_year_rep_success_females$FemaleID, length))
# [1] 1.763158

hist(provisioning_avg_year_rep_success_females$Fledged_Female)

mean(provisioning_avg_year_rep_success_females$Fledged_Female) # 1.597015
std.error(provisioning_avg_year_rep_success_females$Fledged_Female) # 0.1443332

# MCMCglmm ----------------------------------------------------------------
# Specify priors
prior_uni<- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000)))

MCMCglmm_assoc_provisioning_RepSuccess_females_poisson <- MCMCglmm(Fledged_Female ~ 1 +
  scale(mean_Visits_h_Female) * scale(mean_Visits_h_Males) +
  Status_Female +
  ObsPeriod,
random = ~FemaleID, prior = prior_uni,
data = provisioning_avg_year_rep_success_females, singular.ok = TRUE,
nitt = 13000 * 400, burnin = 3000 * 400, thin = 10 * 200,
family = "poisson", pl = TRUE, saveX = TRUE
)

save(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson, file="~/.../.../MCMCglmm_assoc_provisioning_RepSuccess_females_poisson.Rdata")

# Model diagnostics
plot(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$VCV)
plot(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol)
autocorr(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$VCV)
autocorr(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol)

summary(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson)


##### 2.4 Extract results male model  ======================
Intercept <- mean(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 1])
CI_Intercept <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 1])

Coef_z_mean_Visit_Female <- mean(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 2])
CI_z_mean_Visit_Female <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 2])

Coef_z_mean_Visit_Male <- mean(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 3])
CI_z_mean_Visit_Male <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 3])

Coef_Polygamous <- mean(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 4])
CI_Polygamous <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 4])

Coef_BreedingSeason <- mean(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 5])
CI_BreedingSeason <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 5])

Coef_Interaction <- mean(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 6])
CI_Coef_Interaction <- HPDinterval(MCMCglmm_assoc_provisioning_RepSuccess_females_poisson$Sol[, 6])

# Combine everything in one table
PosteriorMean <- rbind(Intercept, Coef_z_mean_Visit_Female, Coef_z_mean_Visit_Male, Coef_Polygamous, Coef_BreedingSeason, Coef_Interaction)
CredibleIntervals <- rbind(CI_Intercept, CI_z_mean_Visit_Female, CI_z_mean_Visit_Male, CI_Polygamous, CI_BreedingSeason, CI_Coef_Interaction)
Predictors <- rbind("Intercept", "z_mean_Visit_Male", "z_mean_Visit_Female", "Polygamous", "BreedingSeason", "Interaction")

results <- data.frame(cbind(Predictors, PosteriorMean, CredibleIntervals))
colnames(results) <- c("Predictors", "PosteriorMean", "CI_low", "CI_high")

# Save table
write.table(results, file = "MCMCglmm_assoc_provisioning_RepSuccess_females_poisson_results.csv", sep = ",", row.names = FALSE)
