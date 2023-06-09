############################################################################
# Project: Assort_Mating_Dunnocks
# Script name: Analyses_assortative_mating
# Script purpose: Test for assortative mating using one measurement per group per year
# Author: Benedikt Holtmann   
# Date Created: 14/08/2020
# R version 4.0.2
############################################################################
# To test for assortative mating, we follow the methods described in:
# Class et al. (2017), A statistical methodology for estimating assortative mating for phenotypic traits that are labile or measured with error. 
# Methods Ecol. Evol. 8, 1910-1919.
# 
# Specifically, we used multivariate mixed models fitting homologous male and female traits to compute between-pair correlations.
# For FID, we fitted seasonal mean FID scores as the response variables as we did not measure FID for each breeding attempt but rather continuously throughout a breeding season.
# In contrast to FID, we recorded provisioning rates for each breeding attempt and pair/group directly on 3-5 consecutive days allowing us to include repeated measurements of provisioning for each unique pair/group.
# Moreover, we only considered first broods for each unique pair/group, since these are less likely to be influenced by the previous breeding experience.
############################################################################

rm(list = ls())

# Load packages ------------------------------------------------------------

library(tidyverse)
library(MCMCglmm)

# Set working directory -----------------------------------------------------

setwd("~ set the working directory")

# Read data ----------------------------------------------------------------

FID_data <- read.csv("02.1_assort_mating_FID_first_broods.csv", na.strings = c("", "NA"))
provisioning_data <- read.csv("02.2_assort_mating_provisioning_first_broods.csv", na.strings = c("", "NA"))

### 1. Assortative mating for FID  ======================

##### 1.1 Bivariate model for monogamous pairs  ======================

FID_data_mono<- filter(FID_data, MatingSystem == "Monogamous")
n_distinct(FID_data_mono$PairID) # 33 Pairs
n_distinct(FID_data_mono$FemaleID) # 27 Females
n_distinct(FID_data_mono$MaleID) # 28 Males

# Remove NAs
FID_data_mono_male_number<- FID_data_mono %>% drop_na(mean_FID_Male)
n_distinct(FID_data_mono_male_number$Male)
# 25
FID_data_mono_female_number<- FID_data_mono %>% drop_na(mean_FID_Female)
n_distinct(FID_data_mono_female_number$Female)
# 22

FID_data_mono <- droplevels(FID_data_mono)

# Observations per season
FID_data_mono %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   14
# 2           2013   12
# 3           2014   12

# For how many pairs do we have repeated samples
FID_data_mono %>% 
  count(PairID, name = "repeats") %>%
  count(repeats) 

# repeats  n
# 1       1 29
# 2       2  3
# 3       3  1

# Mean measurments per unique pair
mean(tapply(FID_data_mono$PairID, FID_data_mono$PairID, length))
# [1] 1.151515

# Check distribution of male and female traits
hist(FID_data_mono$mean_FID_Male, breaks = 50)
hist(FID_data_mono$mean_FID_Female, breaks = 50)

# Remove NAs in fixed effects, necessary for MCMCglmm models.
FID_data_mono_2 <- FID_data_mono %>% drop_na(c(Male_Age, Female_Age))
n_distinct(FID_data_mono_2$PairID) # 32 Pairs
n_distinct(FID_data_mono_2$FemaleID) # 26 Females
n_distinct(FID_data_mono_2$MaleID) # 28 Males
mean(tapply(FID_data_mono_2$PairID, FID_data_mono_2$PairID, length)) # 1.15625

# MCMCglmm ----------------------------------------------------------------

# Specify priors
priorMono <- list(R=list(V=diag(2), nu=0.002), G=list(G1=list(V=diag(2), nu=2, alpha.mu = rep(0, 2), alpha.V= diag(2)*1000)))

Biv_FID_monogamous <- MCMCglmm(cbind(mean_FID_Male, mean_FID_Female) ~ (trait - 1) + # trait-1 allows us to estimate the mean of each trait
  at.level(trait, 1):scale(Male_Age) + at.level(trait, 2):scale(Female_Age) +
  trait:ObsPeriod,
random = ~ us(trait):PairID, # structure of the variance-covariance matrix for the random effects,
# us allows different variances across male and female phenotypes and allows covariances to exist between them.
rcov = ~ us(trait):units, # variance-covariance matrix for the residual variances
prior = priorMono, data = FID_data_mono_2, nitt = 13000 * 100, burnin = 3000 * 100, thin = 10 * 50,
family = c("gaussian", "gaussian")
)

save(Biv_FID_monogamous, file="~/.../.../Biv_FID_monogamous.Rdata")

summary(Biv_FID_monogamous)

# Model diagnostics
plot(Biv_FID_monogamous$VCV)
plot(Biv_FID_monogamous$Sol)
autocorr(Biv_FID_monogamous$VCV)


### Calculate among-pair correlation ----
pair.correlation_FID_mono <- Biv_FID_monogamous$VCV[,"traitmean_FID_Male:traitmean_FID_Female.PairID"]/
  sqrt(Biv_FID_monogamous$VCV[,"traitmean_FID_Male:traitmean_FID_Male.PairID"]*
         Biv_FID_monogamous$VCV[,"traitmean_FID_Female:traitmean_FID_Female.PairID"])

posterior.mode_FID_mono <- posterior.mode(pair.correlation_FID_mono)
# var1 
# 0.3335457  
posterior.mean_FID_mono <- mean(pair.correlation_FID_mono)
# 0.05446708

# 95% Credible intervals
CI_mean_FID_mono <- HPDinterval(pair.correlation_FID_mono)
# lower    upper
# var1 -0.7774429 0.7039008
# attr(,"Probability")
# [1] 0.95

### Calculate within-pair/residual correlation ----
resid.correlation_FID_mono <- Biv_FID_monogamous$VCV[,"traitmean_FID_Male:traitmean_FID_Female.units"]/
  sqrt(Biv_FID_monogamous$VCV[,"traitmean_FID_Male:traitmean_FID_Male.units"]*
         Biv_FID_monogamous$VCV[,"traitmean_FID_Female:traitmean_FID_Female.units"])

posterior.mode_FID_mono_resid <- posterior.mode(resid.correlation_FID_mono)
# var1 
# 0.9946926
posterior.mean_FID_mono_resid <- mean(resid.correlation_FID_mono)
# 0.9743684

# 95% Credible intervals
CI_mean_FID_mono_resid <- HPDinterval(resid.correlation_FID_mono)
# lower    upper
# var1 0.9041754 0.9991643
# attr(,"Probability")
# [1] 0.95


##### 1.2 Trivariate model for polygamnous groups  ======================

FID_data_poly <- filter(FID_data, MatingSystem == "Polygamous")
n_distinct(FID_data_poly$PairID) # 99 unique Pairs
n_distinct(FID_data_poly$GroupID) # 56 unique Groups


# Get separate columns for alpha and beta males 
# Data alpha males
FID_data_poly_alpha<- filter(FID_data_poly, StatusMale == "alpha" & MatingSystem == "Polygamous")
# Rename alpha male columns
FID_data_poly_alpha <- FID_data_poly_alpha %>%
  select(-PairID)  %>% 
  dplyr::rename(
    Alpha = Male, AlphaID = MaleID, AlphaBand = MaleBand, Alpha_Age = Male_Age, StatusA = StatusMale, mean_FID_alpha = mean_FID_Male
  )

FID_data_poly_alpha_number <- FID_data_poly_alpha %>% drop_na(mean_FID_alpha)
n_distinct(FID_data_poly_alpha_number$Alpha)
# 31

# Data beta males
FID_data_poly_beta<- filter(FID_data_poly, StatusMale == "beta" & MatingSystem == "Polygamous")
# Rename alpha male columns
FID_data_poly_beta<- FID_data_poly_beta %>% 
  select(-PairID)  %>% 
  dplyr::rename(
    Beta = Male, BetaID = MaleID, BetaBand = MaleBand, Beta_Age = Male_Age, StatusB = StatusMale, mean_FID_beta = mean_FID_Male
  )

FID_data_poly_beta_number <- FID_data_poly_beta %>% drop_na(mean_FID_beta)
n_distinct(FID_data_poly_beta_number$Beta)
# 36

# Combine alpha and beta males to have them in the same row
FID_data_poly_wide <- left_join(FID_data_poly_alpha, FID_data_poly_beta,
  by = c(
    "NestID", "GroupID", "Female", "FemaleID", "FemaleBand", "Female_Age", "StatusFemale",
    "BreedingSeason", "ObsPeriod", "mean_FID_Female", "MatingSystem"
  )
)


FID_data_poly_wide <- droplevels(FID_data_poly_wide)

FID_data_poly_female_number <- FID_data_poly_wide %>% drop_na(mean_FID_Female)
n_distinct(FID_data_poly_female_number$Female)
# 31

# Observations per season
FID_data_poly_wide %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   17
# 2           2013   21
# 3           2014   26

# For how many pairs do we have repeated samples
FID_data_poly_wide_2 %>% 
  count(GroupID, name = "repeats") %>%
  count(repeats) 

# repeats  n
# 1       1 43
# 2       2  7
# 3       3  1

# Mean measurments per unique group
mean(tapply(FID_data_poly_wide$GroupID, FID_data_poly_wide$GroupID, length))
# [1] 1.163636

# Check distributions
hist(FID_data_poly_wide$mean_FID_alpha, breaks = 50)
hist(FID_data_poly_wide$mean_FID_beta, breaks = 50)
hist(FID_data_poly_wide$mean_FID_Female, breaks = 50)

# Remove NAs in fixed effects, necessary for MCMCglmm models.
FID_data_poly_wide_2 <- FID_data_poly_wide %>% drop_na(c(Alpha_Age, Beta_Age, Female_Age))
n_distinct(FID_data_poly_wide_2$Female) # 37
n_distinct(FID_data_poly_wide_2$Alpha) # 33
n_distinct(FID_data_poly_wide_2$Beta) # 36
n_distinct(FID_data_poly_wide_2$GroupID) # 51
mean(tapply(FID_data_poly_wide_2$GroupID, FID_data_poly_wide_2$GroupID, length)) # 1.176471

# MCMCglmm ----------------------------------------------------------------

# Specify priors
priorPoly <- list(R=list(V=diag(3), nu=0.002), G=list(G1=list(V=diag(3), nu=3, alpha.mu = rep(0, 3), alpha.V= diag(3)*1000)))

Tri_FID_polygamous <- MCMCglmm(cbind(mean_FID_alpha, mean_FID_beta, mean_FID_Female) ~ (trait - 1) + # trait-1 allows us to estimate the mean of each trait
  at.level(trait, 1):scale(Alpha_Age) + at.level(trait, 2):scale(Beta_Age) + at.level(trait, 3):scale(Female_Age) +
  trait:ObsPeriod,
random = ~ us(trait):GroupID, # structure of the variance-covariance matrix for the random effects,
# us allows different variances across male and female phenotypes and allows covariances to exist between them
rcov = ~ us(trait):units, # variance-covariance matrix for the residual variances
prior = priorPoly, data = FID_data_poly_wide_2, nitt = 13000 * 100, burnin = 3000 * 100, thin = 10 * 50,
family = c("gaussian", "gaussian", "gaussian")
)

save(Tri_FID_polygamous, file="~/.../.../Tri_FID_polygamous.Rdata")

summary(Tri_FID_polygamous)

# Model diagnostics
plot(Tri_FID_polygamous$VCV)
plot(Tri_FID_polygamous$Sol)
autocorr(Tri_FID_polygamous$VCV)


### Calculate among-pair correlations ----
### Alpha males and females
pair.correlation_FID_poly_alpha_fem <- Tri_FID_polygamous$VCV[,"traitmean_FID_alpha:traitmean_FID_Female.GroupID"]/
  sqrt(Tri_FID_polygamous$VCV[,"traitmean_FID_alpha:traitmean_FID_alpha.GroupID"]*
         Tri_FID_polygamous$VCV[,"traitmean_FID_Female:traitmean_FID_Female.GroupID"])

posterior.mode_FID_poly_alpha_fem <- posterior.mode(pair.correlation_FID_poly_alpha_fem)
# var1 
# -0.004851628  
posterior.mean_FID_poly_alpha_fem <- mean(pair.correlation_FID_poly_alpha_fem)
# 0.1289425

# 95% Credible intervals
CI_mean_FID_poly_alpha_fem <- HPDinterval(pair.correlation_FID_poly_alpha_fem)
# lower    upper
# var1 -0.411192 0.7582392
# attr(,"Probability")
# [1] 0.95


### Beta males and females
pair.correlation_FID_poly_beta_fem <- Tri_FID_polygamous$VCV[,"traitmean_FID_beta:traitmean_FID_Female.GroupID"]/
  sqrt(Tri_FID_polygamous$VCV[,"traitmean_FID_beta:traitmean_FID_beta.GroupID"]*
         Tri_FID_polygamous$VCV[,"traitmean_FID_Female:traitmean_FID_Female.GroupID"])

posterior.mode_FID_poly_beta_fem <- posterior.mode(pair.correlation_FID_poly_beta_fem)
# var1 
# 0.2527185 
posterior.mean_FID_poly_beta_fem <- mean(pair.correlation_FID_poly_beta_fem)
# 0.186276

# 95% Credible intervals
CI_mean_FID_poly_beta_fem <- HPDinterval(pair.correlation_FID_poly_beta_fem)
# lower    upper
# var1 -0.2635876 0.6525775
# attr(,"Probability")
# [1] 0.95


### Alpha and beta males
pair.correlation_FID_poly_alpha_beta <- Tri_FID_polygamous$VCV[,"traitmean_FID_alpha:traitmean_FID_beta.GroupID"]/
  sqrt(Tri_FID_polygamous$VCV[,"traitmean_FID_alpha:traitmean_FID_alpha.GroupID"]*
         Tri_FID_polygamous$VCV[,"traitmean_FID_beta:traitmean_FID_beta.GroupID"])

posterior.mode_FID_poly_alpha_beta <- posterior.mode(pair.correlation_FID_poly_alpha_beta)
# var1 
# 0.3461024    
posterior.mean_FID_poly_alpha_beta <- mean(pair.correlation_FID_poly_alpha_beta)
# 0.3193115

# 95% Credible intervals
CI_mean_FID_poly_alpha_beta <- HPDinterval(pair.correlation_FID_poly_alpha_beta)
# lower    upper
# var1 -0.2058288 0.765769
# attr(,"Probability")
# [1] 0.95


### Calculate within-pair/residual correlations ----
### Alpha males and females
resid.correlation_FID_poly_alpha_fem <- Tri_FID_polygamous$VCV[,"traitmean_FID_alpha:traitmean_FID_Female.units"]/
  sqrt(Tri_FID_polygamous$VCV[,"traitmean_FID_alpha:traitmean_FID_alpha.units"]*
         Tri_FID_polygamous$VCV[,"traitmean_FID_Female:traitmean_FID_Female.units"])

posterior.mode_FID_poly_alpha_fem_resid <- posterior.mode(resid.correlation_FID_poly_alpha_fem)
# var1 
# -0.1820567   
posterior.mean_FID_poly_alpha_fem_resid <- mean(resid.correlation_FID_poly_alpha_fem)
# -0.09754572

# 95% Credible intervals
CI_mean_FID_poly_alpha_fem_resid <- HPDinterval(resid.correlation_FID_poly_alpha_fem)
# lower    upper
# var1 -0.7399224 0.5625051
# attr(,"Probability")
# [1] 0.95


### Beta males and females
resid.correlation_FID_poly_beta_fem <- Tri_FID_polygamous$VCV[,"traitmean_FID_beta:traitmean_FID_Female.units"]/
  sqrt(Tri_FID_polygamous$VCV[,"traitmean_FID_beta:traitmean_FID_beta.units"]*
         Tri_FID_polygamous$VCV[,"traitmean_FID_Female:traitmean_FID_Female.units"])

posterior.mode_FID_poly_beta_fem_resid <- posterior.mode(resid.correlation_FID_poly_beta_fem)
# var1 
# -0.6124251
posterior.mean_FID_poly_beta_fem_resid <- mean(resid.correlation_FID_poly_beta_fem)
# -0.2818196

# 95% Credible intervals
CI_mean_FID_poly_beta_fem_resid <- HPDinterval(resid.correlation_FID_poly_beta_fem)
# lower    upper
# var1 -0.9018819 0.5352381
# attr(,"Probability")
# [1] 0.95


### Alpha and beta males
resid.correlation_FID_poly_alpha_beta <- Tri_FID_polygamous$VCV[,"traitmean_FID_alpha:traitmean_FID_beta.units"]/
  sqrt(Tri_FID_polygamous$VCV[,"traitmean_FID_alpha:traitmean_FID_alpha.units"]*
         Tri_FID_polygamous$VCV[,"traitmean_FID_beta:traitmean_FID_beta.units"])

posterior.mode_FID_poly_alpha_beta_resid <- posterior.mode(resid.correlation_FID_poly_alpha_beta)
# var1 
# -0.3095147 
posterior.mean_FID_poly_alpha_beta_resid <- mean(resid.correlation_FID_poly_alpha_beta)
# -0.06650597

# 95% Credible intervals
CI_mean_FID_poly_alpha_beta_resid <- HPDinterval(resid.correlation_FID_poly_alpha_beta)
# lower    upper
# var1 -0.7500781 0.6632815
# attr(,"Probability")
# [1] 0.95


##### 1.3 Extract and save results  ======================

Pair.combo <- c("Monogamous", "PolyAlpha", "PolyBeta", "PolyAlphaBeta")

Posterior.mode_cor <- rbind(
  posterior.mode_FID_mono, posterior.mode_FID_poly_alpha_fem,
  posterior.mode_FID_poly_beta_fem, posterior.mode_FID_poly_alpha_beta
)

Posterior.mean_cor <- rbind(
  posterior.mean_FID_mono, posterior.mean_FID_poly_alpha_fem,
  posterior.mean_FID_poly_beta_fem, posterior.mean_FID_poly_alpha_beta
)

CI_mean.cor <- rbind(
  CI_mean_FID_mono, CI_mean_FID_poly_alpha_fem,
  CI_mean_FID_poly_beta_fem, CI_mean_FID_poly_alpha_beta
)

Posterior.mode_resid <- rbind(
  posterior.mode_FID_mono_resid, posterior.mode_FID_poly_alpha_fem_resid,
  posterior.mode_FID_poly_beta_fem_resid, posterior.mode_FID_poly_alpha_beta_resid
)

Posterior.mean_resid <- rbind(
  posterior.mean_FID_mono_resid, posterior.mean_FID_poly_alpha_fem_resid,
  posterior.mean_FID_poly_beta_fem_resid, posterior.mean_FID_poly_alpha_beta_resid
)

CI_mean.resid <- rbind(
  CI_mean_FID_mono_resid, CI_mean_FID_poly_alpha_fem_resid,
  CI_mean_FID_poly_beta_fem_resid, CI_mean_FID_poly_alpha_beta_resid
)

Model <- rbind("Bivariate", "Trivariate", "Trivariate", "Trivariate")

Results <- data.frame(cbind(
  Pair.combo, Posterior.mode_cor, Posterior.mean_cor, CI_mean.cor,
  Posterior.mode_resid, Posterior.mean_resid, CI_mean.resid, Model
))

colnames(Results) <- c(
  "PairCombo", "PosteriorMode", "CorPair", "CI_Cor_Pair_low", "CI_Cor_Pair_high",
 "PosteriorModeResid", "CorResid", "CI_Cor_Resid_low", "CI_Cor_Resid_high", "Model"
)

write.table(Results, file = "Assortative_mating_FID_results.csv", sep = ",", row.names = FALSE)



### 2. Assortative mating for provisioning  ======================

##### 2.1 Bivariate model for monogamous pairs  ======================

provisioning_data_mono<- filter(provisioning_data, MatingSystem == "Monogamous")
n_distinct(provisioning_data_mono$PairID) # 9 Pairs
n_distinct(provisioning_data_mono$FemaleID) # 9 Females
n_distinct(provisioning_data_mono$MaleID) # 9 Males

provisioning_data_mono <- droplevels(provisioning_data_mono)

# Observations per season
provisioning_data_mono %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   10
# 2           2013   18
# 3           2014   8

# Mean measurments per unique pair
mean(tapply(provisioning_data_mono$PairID, provisioning_data_mono$PairID, length))
# [1] 4

# Add sequence per year for each unique pair
provisioning_data_mono <- provisioning_data_mono %>% 
  arrange(NestID, BreedingSeason) %>% 
  group_by(BreedingSeason, PairID) %>%
  mutate(Sequence = row_number()) 
  
provisioning_data_mono <- ungroup(provisioning_data_mono)

provisioning_data_mono <- as.data.frame(provisioning_data_mono)

# Check distributions
# Visits per h
hist(provisioning_data_mono$Visits_h_Male, breaks = 50)
hist(provisioning_data_mono$Visits_h_Female, breaks = 50)

mean(provisioning_data_mono$Visits_h_Male)
mean(provisioning_data_mono$Visits_h_Female)


# Remove NAs in fixed effects, necessary for MCMCglmm models.
provisioning_data_mono_2 <- provisioning_data_mono %>% drop_na(c(Male_Age, Female_Age))

n_distinct(provisioning_data_mono_2$PairID) # 9 Pairs
n_distinct(provisioning_data_mono_2$FemaleID) # 9 Females
n_distinct(provisioning_data_mono_2$MaleID) # 9 Males
mean(tapply(provisioning_data_mono$PairID, provisioning_data_mono$PairID, length)) # 4

# MCMCglmm ----------------------------------------------------------------

# Specify priors
priorMono <- list(R=list(V=diag(2), nu=0.002), G=list(G1=list(V=diag(2), nu=2, alpha.mu = rep(0, 2), alpha.V= diag(2)*1000)))

Biv_Prov_monogamous <- MCMCglmm(cbind(Visits_h_Male, Visits_h_Female) ~ (trait - 1) + # trait-1 allows us to estimate the mean of each trait
  trait:scale(AgeChicks) +
  trait:scale(BroodSize) +
  at.level(trait, 1):scale(Male_Age) + at.level(trait, 2):scale(Female_Age) +
  trait:ObsPeriod,
random = ~ us(trait):PairID, # structure of the variance-covariance matrix for the random effects,
# us allows different variances across male and female phenotypes and allows covariances to exist between them
rcov = ~ us(trait):units, # variance-covariance matrix for the residual variances
prior = priorMono, data = provisioning_data_mono_2, nitt = 13000 * 100, burnin = 3000 * 100, thin = 10 * 50,
family = c("gaussian", "gaussian")
)

save(Biv_Prov_monogamous, file="~/.../...//Biv_Prov_monogamous.Rdata")

summary(Biv_Prov_monogamous)

# Model diagnostics
plot(Biv_Prov_monogamous$VCV)
plot(Biv_Prov_monogamous$Sol)
autocorr(Biv_Prov_monogamous$VCV)


### Calculate among-pair correlations ----
pair.correlation_prov_mono <- Biv_Prov_monogamous$VCV[,"traitVisits_h_Male:traitVisits_h_Female.PairID"]/
  sqrt(Biv_Prov_monogamous$VCV[,"traitVisits_h_Male:traitVisits_h_Male.PairID"]*
         Biv_Prov_monogamous$VCV[,"traitVisits_h_Female:traitVisits_h_Female.PairID"])

posterior.mode_prov_mono <- posterior.mode(pair.correlation_prov_mono)
# var1 
# -0.1929775 
posterior.mean_prov_mono <- mean(pair.correlation_prov_mono)
# 0.01306593

# 95% Credible intervals
CI_mean_prov_mono <- HPDinterval(pair.correlation_prov_mono)
# lower    upper
# var1 -0.8110319 0.9000718
# attr(,"Probability")
# [1] 0.95


### Calculate within-pair/residual correlations ----
resid.correlation_prov_mono <- Biv_Prov_monogamous$VCV[,"traitVisits_h_Male:traitVisits_h_Female.units"]/
  sqrt(Biv_Prov_monogamous$VCV[,"traitVisits_h_Male:traitVisits_h_Male.units"]*
         Biv_Prov_monogamous$VCV[,"traitVisits_h_Female:traitVisits_h_Female.units"])

posterior.mode_prov_mono_resid <- posterior.mode(resid.correlation_prov_mono)
# var1 
# 0.01589673   
posterior.mean_prov_mono_resid <- mean(resid.correlation_prov_mono)
# -0.02813214

# 95% Credible intervals
CI_mean_prov_mono_resid <- HPDinterval(resid.correlation_prov_mono)
# lower    upper
# var1 -0.3951188 0.3647676
# attr(,"Probability")
# [1] 0.95


##### 2.2 Trivariate model for polygamnous groups  ======================

provisioning_data_poly <- filter(provisioning_data, MatingSystem == "Polygamous")
n_distinct(provisioning_data_poly$PairID) # 44 unique Pairs
n_distinct(provisioning_data_poly$GroupID) # 25 unique Groups

provisioning_data_poly <- droplevels(provisioning_data_poly)


# Combine alpha and beta males of the same group in the same line
# Data alpha males
provisioning_data_poly_alpha<- filter(provisioning_data_poly, StatusMale == "alpha" & MatingSystem == "Polygamous")
# Rename alpha male columns
provisioning_data_poly_alpha <- provisioning_data_poly_alpha %>%
  dplyr::rename(
    Alpha = Male, AlphaID = MaleID, AlphaBand = MaleBand, Alpha_Age = Male_Age, StatusA = StatusMale, 
    Visits_Alpha = Visits_Male, DecimalH_Alpha = DecimalH_Male, Visits_min_Alpha = Visits_min_Male, 
    Visits_h_Alpha = Visits_h_Male, Visits_Chick_h_Alpha = Visits_Chick_h_Male
  )

provisioning_data_poly_alpha_number <- provisioning_data_poly_alpha %>% drop_na(Visits_Alpha)
n_distinct(provisioning_data_poly_alpha_number$Alpha)
# 19

# Data beta males
provisioning_data_poly_beta<- filter(provisioning_data_poly, StatusMale == "beta" & MatingSystem == "Polygamous")
# Rename alpha male columns
provisioning_data_poly_beta<- provisioning_data_poly_beta %>% 
  dplyr::rename(
    Beta = Male, BetaID = MaleID, BetaBand = MaleBand, Beta_Age = Male_Age, StatusB = StatusMale, 
    Visits_Beta = Visits_Male, DecimalH_Beta = DecimalH_Male, Visits_min_Beta = Visits_min_Male, 
    Visits_h_Beta = Visits_h_Male, Visits_Chick_h_Beta = Visits_Chick_h_Male
  )

provisioning_data_poly_beta_number <- provisioning_data_poly_beta %>% drop_na(Visits_Beta)
n_distinct(provisioning_data_poly_beta_number$Beta)
# 21

# Combine alpha and beta males to have them in the same row
provisioning_data_poly_wide <- left_join(provisioning_data_poly_alpha, provisioning_data_poly_beta,
  by = c(
    "NestID", "GroupID", "Female", "FemaleID", "FemaleBand", "Female_Age", "StatusFemale", "BreedingSeason",
    "ObsPeriod","Min", "Seconds", "BroodSize", "AgeChicks", "Visits_Female", "DecimalH_Female", "Visits_min_Female", 
    "Visits_h_Female", "Visits_Chick_h_Female"
  )
)

provisioning_data_poly_wide <- droplevels(provisioning_data_poly_wide)

provisioning_data_poly_female_number <- provisioning_data_poly_wide %>% drop_na(Visits_Female)
n_distinct(provisioning_data_poly_female_number$Female)
# 17

# Observations per season
provisioning_data_poly_wide %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   14
# 2           2013   58
# 3           2014   44

# Mean measurments per unique group
mean(tapply(provisioning_data_poly_wide$GroupID, provisioning_data_poly_wide$GroupID, length))
# [1] 4.64

# Add sequence per year for each unique pair
provisioning_data_poly_wide <- provisioning_data_poly_wide %>% 
  arrange(NestID, BreedingSeason) %>% 
  group_by(BreedingSeason, GroupID) %>%
  mutate(Sequence = row_number()) 

provisioning_data_poly_wide <-  ungroup(provisioning_data_poly_wide)
provisioning_data_poly_wide <- as.data.frame(provisioning_data_poly_wide) # make sure that it is a dataframe

# Check distributions
# Visits per h
hist(provisioning_data_poly_wide$Visits_h_Alpha, breaks = 50)
hist(provisioning_data_poly_wide$Visits_h_Beta, breaks = 50)
hist(provisioning_data_poly_wide$Visits_h_Fem, breaks = 50)

mean(provisioning_data_poly_wide$Visits_h_Alpha, na.rm = TRUE)
mean(provisioning_data_poly_wide$Visits_h_Beta, na.rm = TRUE)
mean(provisioning_data_poly_wide$Visits_h_Fem, na.rm = TRUE)

# Remove NAs in fixed effects, necessary for models!
provisioning_data_poly_wide_2 <- provisioning_data_poly_wide %>% drop_na(c(Alpha_Age, Beta_Age, Female_Age))

n_distinct(provisioning_data_poly_wide_2$GroupID) # 21 Trios
n_distinct(provisioning_data_poly_wide_2$FemaleID) # 17 Females
n_distinct(provisioning_data_poly_wide_2$Alpha) # 18 Males
n_distinct(provisioning_data_poly_wide_2$Beta) # 20 Males
mean(tapply(provisioning_data_poly_wide_2$GroupID, provisioning_data_poly_wide_2$GroupID, length)) # 4.71

# MCMCglmm ----------------------------------------------------------------

# Specify priors
priorPoly <- list(R=list(V=diag(3), nu=0.002), G=list(G1=list(V=diag(3), nu=3, alpha.mu = rep(0, 3), alpha.V= diag(3)*1000)))

Tri_Prov_polygamous <- MCMCglmm(cbind(Visits_h_Alpha, Visits_h_Beta, Visits_h_Female) ~ (trait - 1) + # trait-1 allows us to estimate the mean of each trait
  trait:scale(AgeChicks) +
  trait:scale(BroodSize) +
  at.level(trait, 1):scale(Alpha_Age) + at.level(trait, 2):scale(Beta_Age) + at.level(trait, 3):scale(Female_Age) +
  trait:ObsPeriod,
random = ~ us(trait):GroupID, # structure of the variance-covariance matrix for the random effects,
# us allows different variances across male and female phenotypes and allows covariances to exist between them
rcov = ~ us(trait):units, # variance-covariance matrix for the residual variances
prior = priorPoly, data = provisioning_data_poly_wide_2, nitt = 13000 * 100, burnin = 3000 * 100, thin = 10 * 50,
family = c("gaussian", "gaussian", "gaussian")
)

save(Tri_Prov_polygamous, file="~/.../.../Tri_Prov_polygamous.Rdata")

summary(Tri_Prov_polygamous)

# Model diagnostics
plot(Tri_Prov_polygamous$VCV)
plot(Tri_Prov_polygamous$Sol)
autocorr(Tri_Prov_polygamous$VCV)

plot(Tri_Prov_polygamous_poisson$VCV)
plot(Tri_Prov_polygamous_poisson$Sol)
autocorr(Tri_Prov_polygamous_poisson$VCV)


### Calculate among-pair correlations ----
### Alpha males and females
pair.correlation_prov_poly_alpha_fem <- Tri_Prov_polygamous$VCV[,"traitVisits_h_Alpha:traitVisits_h_Female.GroupID"]/
  sqrt(Tri_Prov_polygamous$VCV[,"traitVisits_h_Alpha:traitVisits_h_Alpha.GroupID"]*
         Tri_Prov_polygamous$VCV[,"traitVisits_h_Female:traitVisits_h_Female.GroupID"])

posterior.mode_prov_poly_alpha_fem <- posterior.mode(pair.correlation_prov_poly_alpha_fem)
# var1 
# 0.5730634     
posterior.mean_prov_poly_alpha_fem <- mean(pair.correlation_prov_poly_alpha_fem)
# 0.4673754

# 95% Credible intervals
CI_mean_prov_poly_alpha_fem <- HPDinterval(pair.correlation_prov_poly_alpha_fem)
# lower    upper
# var1 0.0124055 0.9061103
# attr(,"Probability")
# [1] 0.95


### Beta males and females
pair.correlation_prov_poly_beta_fem <- Tri_Prov_polygamous$VCV[,"traitVisits_h_Beta:traitVisits_h_Female.GroupID"]/
  sqrt(Tri_Prov_polygamous$VCV[,"traitVisits_h_Beta:traitVisits_h_Beta.GroupID"]*
         Tri_Prov_polygamous$VCV[,"traitVisits_h_Female:traitVisits_h_Female.GroupID"])

posterior.mode_prov_poly_beta_fem <- posterior.mode(pair.correlation_prov_poly_beta_fem)
# var1 
# -0.3793797
posterior.mean_prov_poly_beta_fem <- mean(pair.correlation_prov_poly_beta_fem)
# -0.2520335

# 95% Credible intervals
CI_mean_prov_poly_beta_fem <- HPDinterval(pair.correlation_prov_poly_beta_fem)
# lower    upper
# var1 -0.7168975 0.2700184
# attr(,"Probability")
# [1] 0.95


### Alpha and beta males
pair.correlation_prov_poly_alpha_beta <- Tri_Prov_polygamous$VCV[,"traitVisits_h_Alpha:traitVisits_h_Beta.GroupID"]/
  sqrt(Tri_Prov_polygamous$VCV[,"traitVisits_h_Alpha:traitVisits_h_Alpha.GroupID"]*
         Tri_Prov_polygamous$VCV[,"traitVisits_h_Beta:traitVisits_h_Beta.GroupID"])

posterior.mode_prov_poly_alpha_beta <- posterior.mode(pair.correlation_prov_poly_alpha_beta)
# var1 
# -0.5024891  
posterior.mean_prov_poly_alpha_beta <- mean(pair.correlation_prov_poly_alpha_beta)
# -0.3854232

# 95% Credible intervals
CI_mean_prov_poly_alpha_beta <- HPDinterval(pair.correlation_prov_poly_alpha_beta)
# lower    upper
# var1 -0.7766571 0.04172941
# attr(,"Probability")
# [1] 0.95


### Calculate within-pair/residual correlations ----
### Alpha males and females
resid.correlation_prov_poly_alpha_fem <- Tri_Prov_polygamous$VCV[,"traitVisits_h_Alpha:traitVisits_h_Female.units"]/
  sqrt(Tri_Prov_polygamous$VCV[,"traitVisits_h_Alpha:traitVisits_h_Alpha.units"]*
         Tri_Prov_polygamous$VCV[,"traitVisits_h_Female:traitVisits_h_Female.units"])

posterior.mode_prov_poly_alpha_fem_resid <- posterior.mode(resid.correlation_prov_poly_alpha_fem)
# var1 
# 0.01065771  
posterior.mean_prov_poly_alpha_fem_resid <- mean(resid.correlation_prov_poly_alpha_fem)
# 0.05999654

# 95% Credible intervals
CI_mean_prov_poly_alpha_fem_resid <- HPDinterval(resid.correlation_prov_poly_alpha_fem)
# lower    upper
# var1 -0.2017646 0.2801807
# attr(,"Probability")
# [1] 0.95


### Beta males and females
resid.correlation_prov_poly_beta_fem <- Tri_Prov_polygamous$VCV[,"traitVisits_h_Beta:traitVisits_h_Female.units"]/
  sqrt(Tri_Prov_polygamous$VCV[,"traitVisits_h_Beta:traitVisits_h_Beta.units"]*
         Tri_Prov_polygamous$VCV[,"traitVisits_h_Female:traitVisits_h_Female.units"])

posterior.mode_prov_poly_beta_fem_resid <- posterior.mode(resid.correlation_prov_poly_beta_fem)
# var1 
# 0.1125322 
posterior.mean_prov_poly_beta_fem_resid <- mean(resid.correlation_prov_poly_beta_fem)
# 0.09062651

# 95% Credible intervals
CI_mean_prov_poly_beta_fem_resid <- HPDinterval(resid.correlation_prov_poly_beta_fem)
# lower    upper
# var1 -0.1294914 0.3285948
# attr(,"Probability")
# [1] 0.95


### Alpha and beta males
resid.correlation_prov_poly_alpha_beta <- Tri_Prov_polygamous$VCV[,"traitVisits_h_Alpha:traitVisits_h_Beta.units"]/
  sqrt(Tri_Prov_polygamous$VCV[,"traitVisits_h_Alpha:traitVisits_h_Alpha.units"]*
         Tri_Prov_polygamous$VCV[,"traitVisits_h_Beta:traitVisits_h_Beta.units"])

posterior.mode_prov_poly_alpha_beta_resid <- posterior.mode(resid.correlation_prov_poly_alpha_beta)
# var1 
# -0.07645712 
posterior.mean_prov_poly_alpha_beta_resid <- mean(resid.correlation_prov_poly_alpha_beta)
# -0.07686552

# 95% Credible intervals
CI_mean_prov_poly_alpha_beta_resid <- HPDinterval(resid.correlation_prov_poly_alpha_beta)
# lower    upper
# var1 -0.2993149 0.1427212
# attr(,"Probability")
# [1] 0.95


##### 2.3 Extract and save results  ======================

Pair.combo <- c("Monogamous", "PolyAlpha", "PolyBeta", "PolyAlphaBeta")

Posterior.mode_cor <- rbind(
  posterior.mode_prov_mono, posterior.mode_prov_poly_alpha_fem,
  posterior.mode_prov_poly_beta_fem, posterior.mode_prov_poly_alpha_beta
)

Posterior.mean_cor <- rbind(
  posterior.mean_prov_mono, posterior.mean_prov_poly_alpha_fem,
  posterior.mean_prov_poly_beta_fem, posterior.mean_prov_poly_alpha_beta
)

CI_mean.cor <- rbind(
  CI_mean_prov_mono, CI_mean_prov_poly_alpha_fem,
  CI_mean_prov_poly_beta_fem, CI_mean_prov_poly_alpha_beta
)

Posterior.mode_resid <- rbind(
  posterior.mode_prov_mono_resid, posterior.mode_prov_poly_alpha_fem_resid,
  posterior.mode_prov_poly_beta_fem_resid, posterior.mode_prov_poly_alpha_beta_resid
)

Posterior.mean_resid <- rbind(
  posterior.mean_prov_mono_resid, posterior.mean_prov_poly_alpha_fem_resid,
  posterior.mean_prov_poly_beta_fem_resid, posterior.mean_prov_poly_alpha_beta_resid
)

CI_mean.resid <- rbind(
  CI_mean_prov_mono_resid, CI_mean_prov_poly_alpha_fem_resid,
  CI_mean_prov_poly_beta_fem_resid, CI_mean_prov_poly_alpha_beta_resid
)

Model <- rbind("Bivariate", "Trivariate", "Trivariate", "Trivariate")


Results <- data.frame(cbind(
  Pair.combo, Posterior.mode_cor, Posterior.mean_cor, CI_mean.cor,
  Posterior.mode_resid, Posterior.mean_resid, CI_mean.resid, Model
))

colnames(Results) <- c(
  "PairCombo", "PosteriorMode", "CorPair", "CI_Cor_Pair_low", "CI_Cor_Pair_high",
  "PosteriorModeResid", "CorResid", "CI_Cor_Resid_low", "CI_Cor_Resid_high", "Model"
)

write.table(Results, file = "Assortative_mating_Prov_results.csv", sep = ",", row.names = FALSE)
