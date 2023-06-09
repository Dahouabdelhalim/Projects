############################################################################
# Project: Assortative_mating_GT_exploration_JEB
# Script name: Estimating_indirect_effects_and_assortative_mating_exploration
# Script purpose: Estimate indirect effects and assortative mating for behavioural and morphometric traits
# Author: Benedikt Holtmann
# Date Created: 15/06/2020
# R version 4.0.2
############################################################################
# Below we provide the R-code used for the analyses of exploratory behaviour
# All other traits have been analysed in same way.

# Please note that breathig rates have only been measured from 2015 onward.
# Thus, for the models on breathing rates we filterd the data by:
# exploration_data_2015 <- exploration_data %>% filter(Year >= 2015)

# Please also note that observer identities are called differently for each trait:
# FieldObserver: exploration
# MorphObserver: body mass, tarsus, wing
# TestObserver: breathing rate
############################################################################

rm(list = ls())

# Load packages ------------------------------------------------------------
library(tidyverse) # Formatting and manipulating data
library(MCMCglmm)


# Set working directory -----------------------------------------------------

setwd("~ set the working directory")

# 1. Read data ----------------------------------------------------------------
assort_mating_data <- read.csv("Data_Holtmann_&_Dingemanse_JEB_Parus_major_assort_mating.csv", na.strings = c("", "-99", "NA"), stringsAsFactors = FALSE)

# Get range of explorartion score
range(assort_mating_data$ExpScoreMale, na.rm = TRUE)
range(assort_mating_data$ExpScoreFemale, na.rm = TRUE)

hist(assort_mating_data$ExpScoreMale, breaks = 50)
hist(assort_mating_data$ExpScoreFemale, breaks = 50)

assort_mating_data <- as.data.frame(assort_mating_data)
class(assort_mating_data)

# Get sample sizes
n_distinct(assort_mating_data$FemaleID) # 990 Females
n_distinct(assort_mating_data$MaleID) # 979 Males
n_distinct(assort_mating_data$PairID) # 1260 Pairs

n_distinct(assort_mating_data$Plot) # 12 Plots
n_distinct(assort_mating_data$Year) # 10 Years
n_distinct(assort_mating_data$PlotYear) # 118 Plot-year
n_distinct(assort_mating_data$NestBox) # 531 Nestboxes
n_distinct(assort_mating_data$TestsObserverFemale) # 58 Observers
n_distinct(assort_mating_data$TestsObserverMale) # 58 Observers

assort_mating_data$TestsObserverFemale <- as.factor(assort_mating_data$TestsObserverFemale)
levels(assort_mating_data$TestsObserverFemale)
assort_mating_data$TestsObserverMale <- as.factor(assort_mating_data$TestsObserverMale)
levels(assort_mating_data$TestsObserverMale)


# Mean meaurements per pair
mean(tapply(assort_mating_data$PairID, assort_mating_data$PairID, length)) 
# [1] 1.235714


# 2. Analyses for exploratory behaviour ----------------------------------------------------------

### 2.1 Multivariate model including Female and Male ID as random effects (Model A)  ======================
# MCMCglmm ----

# Specify priors
prior_multi_1 <- list(R = list(V = diag(2), nu = 0.002), G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G2 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G3 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G4 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G5 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G6 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G7 = list(V = diag(4), nu = 4, alpha.mu = rep(0, 4), alpha.V = diag(4) * 1000)
))

# Bivariate model with trait specific fixed effects
Biv_MCMC_exploration_IndIDs <- MCMCglmm(cbind(scale(ExpScoreFemale), scale(ExpScoreMale)) ~ (trait - 1) + # trait-1 allows us to estimate the mean of each trait
  at.level(trait, 1):scale(AgeFemale) +
  at.level(trait, 1):I(scale(AgeFemale)^2) +
  at.level(trait, 1):scale(DecTimeFem) + # as decimal
  at.level(trait, 1):SequenceFemale +
  at.level(trait, 2):scale(AgeMale) +
  at.level(trait, 2):I(scale(AgeMale)^2) +
  at.level(trait, 2):scale(DecTimeMale) + # as decimal
  at.level(trait, 2):SequenceMale,
random = ~ us(trait):FemaleID + us(trait):MaleID + us(trait):Plot + us(trait):Year +
  us(trait):PlotYear + us(trait):NestBox + us(trait):str(FieldObserverFemale + FieldObserverMale), # structure of the variance-covariance matrix for the random effects,
# us allows different variances across male and female phenotypes and allows covariances to exist between them
rcov = ~ us(trait):units, # variance-covariance matrix for the residual variances
prior = prior_multi_1, data = assort_mating_data, nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
pr = TRUE, saveX = TRUE, saveZ = TRUE,
family = c("gaussian", "gaussian")
)

# Save model
# save(Biv_MCMC_exploration_IndIDs, file="~/.../.../Biv_MCMC_exploration_IndIDs.Rdata") # change to path/folder were you wnat to save your model
load("~/.../.../Biv_MCMC_exploration_IndIDs.Rdata")

summary(Biv_MCMC_exploration_IndIDs)

# Model diagnostics
plot(Biv_MCMC_exploration_IndIDs$Sol)
plot(Biv_MCMC_exploration_IndIDs$VCV)

autocorr.plot(Biv_MCMC_exploration_IndIDs$VCV)
autocorr(Biv_MCMC_exploration_IndIDs$VCV)

heidel.diag(Biv_MCMC_exploration_IndIDs$VCV)
heidel.diag(Biv_MCMC_exploration_IndIDs$Sol)


# 2.1.1 Calculate adjusted repeatabilities (direct and indirect effects) ----

# For Female exploration ====
# Females (direct effect)
FemRep_Expl_Fem <- (Biv_MCMC_exploration_IndIDs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.FemaleID"] /
  (Biv_MCMC_exploration_IndIDs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.FemaleID"] +
    Biv_MCMC_exploration_IndIDs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.MaleID"] +
    Biv_MCMC_exploration_IndIDs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_IndIDs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_IndIDs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_IndIDs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_IndIDs$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_IndIDs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"]))

# Repeatability
mean(FemRep_Expl_Fem)
# 0.34625

# 95% Credible intervals
HPDinterval(FemRep_Expl_Fem) 
# lower     upper
# var1 0.2777766 0.4172573
# attr(,"Probability")
# [1] 0.95

# Males (indirect effect on female exploration)
MaleRep_Expl_Fem <- (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.MaleID"] /
  (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.MaleID"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.FemaleID"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"]))

# Repeatability
mean(MaleRep_Expl_Fem)
# 0.02993924

# 95% Credible intervals
HPDinterval(MaleRep_Expl_Fem)
# lower      upper
# var1 3.259979e-09 0.08156807
# attr(,"Probability")
# [1] 0.95

# For male exploration ====
# Males (direct effect)
MaleRep_expl_Male <- (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.MaleID"] /
  (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.MaleID"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.FemaleID"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.units"]))

# Repeatability
mean(MaleRep_expl_Male) 
# 0.3027592

# 95% Credible intervals
HPDinterval(MaleRep_expl_Male) 
# lower     upper
# var1 0.2283531 0.3666248
# attr(,"Probability")
# [1] 0.95

# Females (indirect effect on male exploration)
FemRep_expl_Male <- (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.FemaleID"] /
  (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.FemaleID"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.MaleID"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.units"]))

# Repeatability
mean(FemRep_expl_Male)
# 0.01604696

# 95% Credible intervals
HPDinterval(FemRep_expl_Male)
# lower      upper
# var1 3.514759e-10 0.05735077
# attr(,"Probability")
# [1] 0.95


# Calculate average individual-level repeatability for exploration
Rep.MeanExp <- sqrt((Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.FemaleID"] /
  (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.FemaleID"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.MaleID"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"])) *
  (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.MaleID"] /
    (Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.MaleID"] +
      Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.FemaleID"] +
      Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
      Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
      Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
      Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
      Biv_MCMC_exploration_IndIDs_all_obs$VCV[, "traitExpScoreMale:traitExpScoreMale.units"])))

# Repeatability
mean(Rep.MeanExp)
# 0.3227857

# 95% Credible intervals
HPDinterval(Rep.MeanExp)
# lower     upper
# var1 0.269448 0.3676142
# attr(,"Probability")
# [1] 0.95


### 2.2 Bivariate model including PairID as random effect (Model B)  ======================
# MCMCglmm ----

# Specify priors
prior_multi_2 <- list(R = list(V = diag(2), nu = 0.002), G = list(
  G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G2 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G3 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G4 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G5 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2) * 1000),
  G6 = list(V = diag(4), nu = 4, alpha.mu = rep(0, 4), alpha.V = diag(4) * 1000)
))


# Specify model
# Model with trait specific fixed effects
Biv_MCMC_exploration_pairID <- MCMCglmm(cbind(scale(ExpScoreFemale), scale(ExpScoreMale)) ~ (trait - 1) + # trait-1 allows us to estimate the mean of each trait
  at.level(trait, 1):scale(AgeFemale) +
  at.level(trait, 1):I(scale(AgeFemale)^2) +
  at.level(trait, 1):scale(DecTimeFem) +
  at.level(trait, 1):SequenceFemale +
  at.level(trait, 2):scale(AgeMale) +
  at.level(trait, 2):I(scale(AgeMale)^2) +
  at.level(trait, 2):scale(DecTimeMale) +
  at.level(trait, 2):SequenceMale,
random = ~ us(trait):PairID + us(trait):Plot + us(trait):Year + us(trait):PlotYear + us(trait):NestBox +
  us(trait):str(FieldObserverFemale + FieldObserverMale), 
rcov = ~ us(trait):units, 
prior = prior_multi_2, data = assort_mating_data, nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
pr = TRUE, saveX = TRUE, saveZ = TRUE,
family = c("gaussian", "gaussian")
)

# save(Biv_MCMC_exploration_pairID, file="~/.../.../Biv_MCMC_exploration_pairID.Rdata")
load("~/.../.../Biv_MCMC_exploration_pairID.Rdata")

summary(Biv_MCMC_exploration_pairID)

# Model diagnostics
plot(Biv_MCMC_exploration_pairID$Sol)
plot(Biv_MCMC_exploration_pairID$VCV)

autocorr.plot(Biv_MCMC_exploration_pairID$VCV)

heidel.diag(Biv_MCMC_exploration_pairID$VCV)
heidel.diag(Biv_MCMC_exploration_pairID$Sol)


### 2.3 Calculate correlations ======================

### 2.3.1 Between-pair correlation (i.e. assortative mating) ----
pair.correlation_expl <- Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.PairID"] /
  sqrt(Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] *
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"])

# Posterior mode
Pair.PosteriorMode_expl <- posterior.mode(pair.correlation_expl)
# var1 
# 0.02289836  

# Posterior mean
Pair.PosteriorMean_expl <- mean(pair.correlation_expl)
# 0.05813303

# 95% Credible intervals
Pair.CI_Mean_expl <- HPDinterval(pair.correlation_expl)
# lower     upper
# var1 -0.1293717 0.2502615
# attr(,"Probability")
# [1] 0.95

# Calculate geometric mean repeatability for pair
Pair.GeomMeanExp <- sqrt((Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] /
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"])) *
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] /
    (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
      Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"])))

# Geometric mean repeatability
mean(Pair.GeomMeanExp)
# 0.3362674

# 95% credible intervals
HPDinterval(Pair.GeomMeanExp)
# lower     upper
# var1 0.2716115 0.3937992
# attr(,"Probability")
# [1] 0.95


### 2.3.2 Plot-level correlation ----
plot.correlation_expl <- Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.Plot"] /
  sqrt(Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] *
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"])

# Posterior mode
Plot.PosteriorMode_expl <- posterior.mode(plot.correlation_expl)
# var1
# 0.8283393

# Posterior mean
Plot.PosteriorMean_expl <- mean(plot.correlation_expl)
#0.5655275

# 95% Credible intervals
Plot.CI_Mean_expl <- HPDinterval(plot.correlation_expl)
# lower     upper
# var1 -0.007319196 0.9877738
# attr(,"Probability")
# [1] 0.95

# Calculate geometric mean repeatability for plot
Plot.GeomMeanExp <- sqrt((Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] /
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"])) *
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] /
    (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
      Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:FieldObserverMale:traitExpScoreMale"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"])))

# Geometric mean repeatability
mean(Plot.GeomMeanExp)
# 0.02241077

# 95% credible interval
HPDinterval(Plot.GeomMeanExp)
# lower      upper
# var1 0.003782382 0.04923053
# attr(,"Probability")
# [1] 0.95


### 2.3.3 Year-level correlation ----
year.correlation_expl <- Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.Year"] /
  sqrt(Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] *
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"])

# Posterior mode
Year.PosteriorMode_expl <- posterior.mode(year.correlation_expl)
# var1 
# 0.6658497 

# Posterior mean
Year.PosteriorMean_expl <- mean(year.correlation_expl)
# 0.4388143

# 95% Credible intervals
Year.CI_Mean_expl <- HPDinterval(year.correlation_expl)
# lower     upper
# var1 -0.1717393 0.9288838
# attr(,"Probability")
# [1] 0.95

# Calculate geometric mean repeatability for year
Year.GeomMeanExp <- sqrt((Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] /
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"])) *
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] /
    (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
      Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"])))

# Geometric mean repeatability
mean(Year.GeomMeanExp)
# 0.03450941

# 95% credible intervals
HPDinterval(Year.GeomMeanExp)
# lower      upper
# var1 0.006930181 0.06858512
# attr(,"Probability")
# [1] 0.95


### 2.3.4 PlotYear-level correlation ----
plotyear.correlation_expl <- Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.PlotYear"] /
  sqrt(Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] *
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"])

# Posterior mode
PlotYear.PosteriorMode_expl <- posterior.mode(plotyear.correlation_expl)
# var1
# 0.3912256 

# Posterior mean
PlotYear.PosteriorMean_expl <- mean(plotyear.correlation_expl)
# 0.06791574

# 95% Credible intervals
PlotYear.CI_Mean_expl <- HPDinterval(plotyear.correlation_expl)
# lower     upper
# var1 -0.7951757 0.9127234
# attr(,"Probability")
# [1] 0.95

# Calculate geometric mean repeatability for PlotYear
PlotYear.GeomMeanExp <- sqrt((Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] /
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"])) *
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] /
    (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
      Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"])))

# Geometric mean repeatability
mean(PlotYear.GeomMeanExp)
# 0.004168678

# 95% credible intervals
HPDinterval(PlotYear.GeomMeanExp)
# lower      upper
# var1 8.33947e-07 0.01243394
# attr(,"Probability")
# [1] 0.95


### 2.3.5 NestBox-level correlation ----
nestbox.correlation_expl <- Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.NestBox"] /
  sqrt(Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] *
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"])

# Posterioe mode
NestBox.PosteriorMode_expl <- posterior.mode(nestbox.correlation_expl)
# var1
# -0.3148954 

# Posterior mean
NestBox.PosteriorMean_expl <- mean(nestbox.correlation_expl)
# -0.067203156

# 95% Credible intervals
NestBox.CI_Mean_expl <- HPDinterval(nestbox.correlation_expl)
# lower     upper
# var1 -0.9579492 0.7558148
# attr(,"Probability")
# [1] 0.95


# Calculate geometric mean repeatability for NestBox
NestBox.GeomMeanExp <- sqrt((Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] /
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"])) *
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] /
    (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"])))

# Geometric mean repeatability
mean(NestBox.GeomMeanExp)
# 0.009475218

# 95% credible intervals 
HPDinterval(NestBox.GeomMeanExp)
# lower      upper
# var1 7.31432e-06 0.02722537
# attr(,"Probability")
# [1] 0.95


### 2.3.6 Observer-level correlation ----
obs.correlation_expl <- Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverMale:traitExpScoreMale"] /
  sqrt(Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] *
    Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"])

# Posterior mode
Obs.PosteriorMode_expl <- posterior.mode(nestbox.correlation_expl)
# var1 
# -0.3148954  

# Posterior mean
Obs.PosteriorMean_expl <- mean(nestbox.correlation_expl)
# -0.06720315

# 95% Credible intervals
Obs.CI_Mean_expl <- HPDinterval(nestbox.correlation_expl)
# lower     upper
# var1 -0.9579492 0.7558148
# attr(,"Probability")
# [1] 0.95


# Calculate geometric mean repeatability for observer
Observer.GeomMeanExp <- sqrt((Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] /
  (Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"])) *
  (Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] /
    (Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"])))

# Geometric mean repeatability
mean(Observer.GeomMeanExp)
# 0.00471408

# 95% Credible intervals
HPDinterval(Observer.GeomMeanExp)
# lower    upper
# var1 6.919399e-08 0.014336
# attr(,"Probability")
# [1] 0.95


### 2.3.7 Within-pair (residual) correlation ----
resid.correlation_expl <- Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.units"] /
  sqrt(Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"] *
         Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"])

# Posterior mode
Resid.PosteriorMode_expl <- posterior.mode(resid.correlation_expl)
# var1
# 0.1192411  

# Posterior mean
Resid.PosteriorMean_expl <- mean(resid.correlation_expl)
# 0.1364356

# 95% Credible intervals
Resid.CI_Mean_expl <- HPDinterval(resid.correlation_expl)
# lower     upper
# var1 0.03302035 0.2358914
# attr(,"Probability")
# [1] 0.95


# # Calculate geometric mean repeatability for residuals
resid.GeomMeanExp <- sqrt((Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"] /
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"] +
     Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] +
     Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] +
     Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
     Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
     Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] +
     Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"])) *
  (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"] /
    (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"] +
       Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] +
       Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] +
       Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
       Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
       Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] +
       Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"])))

# Geometric mean repeatability
mean(resid.GeomMeanExp)
# 0.5642698

# 95% Credible intervals
HPDinterval(resid.GeomMeanExp)
# lower    upper
# var1 0.5015277 0.627414
# attr(,"Probability")
# [1] 0.95


### 2.3.8 Phenotpypic correlation ----
phen.correlation_expl <- (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.PairID"] + Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.Plot"] +
  Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] + Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.PlotYear"] +
  Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.NestBox"] + Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverMale:traitExpScoreMale"] +
  Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreMale.units"]) /
  sqrt((Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PairID"] + Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Plot"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.Year"] + Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.PlotYear"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.NestBox"] + Biv_MCMC_exploration_pairID$VCV[, "FieldObserverFemale:traitExpScoreFemale.FieldObserverFemale:traitExpScoreFemale"] +
    Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreFemale:traitExpScoreFemale.units"]) *
    (Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PairID"] + Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Plot"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.Year"] + Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.PlotYear"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.NestBox"] + Biv_MCMC_exploration_pairID$VCV[, "FieldObserverMale:traitExpScoreMale.FieldObserverMale:traitExpScoreMale"] +
      Biv_MCMC_exploration_pairID$VCV[, "traitExpScoreMale:traitExpScoreMale.units"]))

# Posterior mode
Phen.PosteriorMode_expl <- posterior.mode(phen.correlation_expl)
# var1 
# 0.1301037 

# Posterior mean
Phen.PosteriorMean_expl <- mean(phen.correlation_expl)
# 0.1386944

# 95% Credible intervals
Phen.CI_Mean_expl <- HPDinterval(phen.correlation_expl)
# lower     upper
# var1 0.06929815 0.205547
# attr(,"Probability")
# [1] 0.95


# 3. Extract results and save in one table ----------------------------------------------------------
# Combine everything in one table
Correlation <- rbind(Pair.PosteriorMean_expl, Plot.PosteriorMean_expl, Year.PosteriorMean_expl, PlotYear.PosteriorMean_expl, NestBox.PosteriorMean_expl, Obs.PosteriorMean_expl, Resid.PosteriorMean_expl, Phen.PosteriorMean_expl)
CredibleIntervals <- rbind(Pair.CI_Mean_expl, Plot.CI_Mean_expl, Year.CI_Mean_expl, PlotYear.CI_Mean_expl, NestBox.CI_Mean_expl, Obs.CI_Mean_expl, Resid.CI_Mean_expl, Phen.CI_Mean_expl)
Parameter <- rbind("Pair", "Plot", "Year", "Plot-Year", "Nestbox",  "Observer", "Residual", "Phenotypic")

results <- data.frame(cbind(Parameter, Correlation, CredibleIntervals))
colnames(results) <- c("Parameter", "Correlation", "CI_low", "CI_high")

# Save table
write.table(results, file = "Biv_MCMC_exploration_pairID_results.csv", sep = ",", row.names = FALSE)
