############################################################################
# Project: Assort_Mating_Dunnocks
# Script name: Repeatability_of_behavioural_traits
# Script purpose: Estimate between- and within-individual variance and calculate repeatabilities
# Author: Benedikt Holtmann
# Date Created: 24/07/2020
# R version 4.0.2
############################################################################

rm(list = ls())

# Load packages ------------------------------------------------------------

library(tidyverse) # Formatting and manipulating data
library(MCMCglmm) # Bayesian framework using Markov chain Monte Carlo (MCMC) methods
library(lubridate) # Manipulate dates
library(plotrix)

# Set working directory -----------------------------------------------------

setwd("~ set the working directory ")

# Read data ----------------------------------------------------------------

activity_data <- read.csv("01.1_repeat_activity_data.csv", na.strings = c("", "NA"))
FID_data <- read.csv("01.2_repeat_FID_data.csv", na.strings = c("", "NA"))
vigilance_data <- read.csv("01.3_repeat_vigilance_data.csv", na.strings = c("", "NA"))
provisioning_data <- read.csv("01.4_repeat_provisioning_data.csv", na.strings = c("", "NA"))


### 1. Repeatability of activity  ======================

# Overview of activity measurements ----
# Total number of measurments for males and females
summary(activity_data$Sex)
# Female   Male   NA's
#    228    416      7

# Mean seconds per observations
mean(activity_data$Seconds) # 83.59601
sd(activity_data$Seconds) # 30.93275

# Mean hops per observation
mean(activity_data$Hops) # 85.08449
sd(activity_data$Hops) # 38.70763
se<- sd(activity_data$Hops)/length(activity_data$Hops)
std.error(activity_data$Hops)

# Observations per Season
activity_data %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   142
# 2           2013   249
# 3           2014   260

# Get the number of individuals for each Season
activity_data %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    48
# 2           2013    77
# 3           2014    78

# Mean measurments per bird
mean(tapply(activity_data$BirdID, activity_data$BirdID, length))
# [1] 5.425

# Change date into number of days starting from 1st September
# Specify that Date is a date in R
activity_data$Date <- as.character(as.Date(activity_data$Date, format = "%d/%m/%Y"))

# Reference dates for each season
RefDate2012 <- as.Date("01/09/12", format = "%d/%m/%y")
RefDate2013 <- as.Date("01/09/13", format = "%d/%m/%y")
RefDate2014 <- as.Date("01/09/14", format = "%d/%m/%y")

# Calculate SeptDay for each season
activity_data <- activity_data %>%
  mutate(SeptDay = ifelse(BreedingSeason == 2012, difftime(activity_data$Date, RefDate2012, units = c("days")),
    ifelse(BreedingSeason == 2013, difftime(activity_data$Date, RefDate2013, units = c("days")),
      difftime(activity_data$Date, RefDate2014, units = c("days"))
    )
  ))


# Remove NAs in fixed effects, necessary for MCMCglmm models.
activity_data <- activity_data %>% drop_na(c(Seconds, Time, SeptDay, Age))


### 1.1 Male model  ======================

activity_data_m <- activity_data %>% filter(Sex == 'Male')
activity_data_m <- droplevels(activity_data_m)

hist(activity_data_m$Hops, breaks = 50)

# Number of males 
activity_data_m %>%
  distinct(BirdID) %>%
  tally()
# 73

# Observations per Season
activity_data_m %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    91
# 2           2013   161
# 3           2014   163

# Get the number of males for each Season
activity_data_m %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    29
# 2           2013    47
# 3           2014    47

# Mean measurments per male
mean(tapply(activity_data_m$BirdID, activity_data_m$BirdID, length))
# [1] 5.684932


# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) # parameter expanded prior

# Model for activity all
RepAct_males <- MCMCglmm(Hops ~ 1 +
  scale(DecimalMin) + scale(Time) + scale(SeptDay) + scale(Age),
random = ~BirdID,
prior = prior1,
data = activity_data_m,
nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
family = "gaussian"
)

save(RepAct_males, file="~/.../.../RepAct_males.Rdata")

summary(RepAct_males)

# Model diagnostics
plot(RepAct_males$VCV)
plot(RepAct_males$Sol)

# Extract among-individual variance and credible intervals
VarAmong_Act_m <- mean(RepAct_males$VCV[, 1])
VarAmong_Act_CI_m <- HPDinterval(RepAct_males$VCV[, 1])

# Extract within-individual variance and credible intervals
VarWithin_Act_m <- mean(RepAct_males$VCV[, 2])
VarWithin_Act_CI_m <- HPDinterval(RepAct_males$VCV[, 2])

# Calculate repeatability for activity
RepeatAct_m <- mean(RepAct_males$VCV[, 1] / (RepAct_males$VCV[, 1] + RepAct_males$VCV[, 2]))
# [1] 0.1216713

# 95% Credible intervals
RepeatAct_CI_m <- HPDinterval(RepAct_males$VCV[, 1] / (RepAct_males$VCV[, 1] + RepAct_males$VCV[, 2]))
# lower     upper
# var1 0.02968203 0.2291111
# attr(,"Probability")
# [1] 0.95


### 1.2 Female model  ======================

activity_data_f <- activity_data %>% filter(Sex == 'Female')
activity_data_f <- droplevels(activity_data_f)

hist(activity_data_f$Hops, breaks = 50)

# Number of females 
activity_data_f %>%
  distinct(BirdID) %>%
  tally()
# 43

# Observations per Season
activity_data_f %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    50
# 2           2013    82
# 3           2014    96

# Get the number of males for each Season
activity_data_f %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    19
# 2           2013    27
# 3           2014    30

# Mean measurments per male
mean(tapply(activity_data_f$BirdID, activity_data_f$BirdID, length))
# [1] 5.302326


# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) # parameter expanded prior

# Model for activity all
RepAct_females <- MCMCglmm(Hops ~ 1 +
  scale(DecimalMin) + scale(Time) + scale(SeptDay) + scale(Age),
random = ~BirdID,
prior = prior1,
data = activity_data_f,
nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
family = "gaussian"
)

save(RepAct_females, file="~/.../.../RepAct_females.Rdata")

summary(RepAct_females)

# Model diagnostics
plot(RepAct_females$VCV)
plot(RepAct_females$Sol)

# Extract among-individual variance and credible intervals
VarAmong_Act_f <- mean(RepAct_females$VCV[, 1])
VarAmong_Act_CI_f <- HPDinterval(RepAct_females$VCV[, 1])

# Extract within-individual variance and credible intervals
VarWithin_Act_f <- mean(RepAct_females$VCV[, 2])
VarWithin_Act_CI_f <- HPDinterval(RepAct_females$VCV[, 2])

# Calculate repeatability for activity
RepeatAct_f <- mean(RepAct_females$VCV[, 1] / (RepAct_females$VCV[, 1] + RepAct_females$VCV[, 2]))
# [1] 0.06811202

# 95% Credible intervals
RepeatAct_CI_f <- HPDinterval(RepAct_females$VCV[, 1] / (RepAct_females$VCV[, 1] + RepAct_females$VCV[, 2]))
# lower     upper
# var1 3.384091e-08 0.1763171
# attr(,"Probability")
# [1] 0.95

# Calculate coefficient of variation (CVb)
CVb_Act_f <- mean(sqrt(RepAct_females$VCV[, 1]) / RepAct_females$Sol[, 1]) # among-individual variance divided by intercept
# [1] 0.07742327



### 2. Repeatability of FID  ======================

# Overview of FID measurements ----
# Total number of measurments for males and females
summary(FID_data$Sex)
# Female   Male   NA's
#    248    439      3

# Observations per Season
FID_data %>% group_by(BreedingSeason) %>% tally()
# 1           2012   174
# 2           2013   252
# 3           2014   264

# Get the number of individuals for each Season
FID_data %>% group_by(BreedingSeason) %>% distinct(BirdID) %>% tally()
# 1           2012    57
# 2           2013    75
# 3           2014    75

# Mean measurments per bird
mean(tapply(FID_data$BirdID, FID_data$BirdID, length))
# [1] 5.847458


# Change date into number of days starting from 1st September
# Specify that Date is a date in R
FID_data$Date <- as.character(as.Date(FID_data$Date, format = "%d/%m/%Y"))

# Reference dates for each season
RefDate2012 <- as.Date("01/09/12", format = "%d/%m/%y")
RefDate2013 <- as.Date("01/09/13", format = "%d/%m/%y")
RefDate2014 <- as.Date("01/09/14", format = "%d/%m/%y")

# Calculate SeptDay for each season
FID_data <- FID_data %>%
  mutate(SeptDay = ifelse(BreedingSeason == 2012, difftime(FID_data$Date, RefDate2012, units = c("days")),
    ifelse(BreedingSeason == 2013, difftime(FID_data$Date, RefDate2013, units = c("days")),
      difftime(FID_data$Date, RefDate2014, units = c("days"))
    )
  ))


# Remove NAs in fixed effects, necessary for models!
FID_data <- FID_data %>% drop_na(c(Time, Age))


### 2.1 Male model  ======================

FID_data_m <- FID_data  %>% filter(Sex == 'Male')
FID_data_m <- droplevels(FID_data_m)

hist(FID_data_m$FID, breaks = 50)

# Number of males
FID_data_m %>%
  distinct(BirdID) %>%
  tally()
# 72

# Observations per Season
FID_data_m %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   115
# 2           2013   165
# 3           2014   159

# Get the number of males for each Season
FID_data_m %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    36
# 2           2013    46
# 3           2014    44

# Mean measurments per male
mean(tapply(FID_data_m$BirdID, FID_data_m$BirdID, length))
# [1] 6.097222


# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) # parameter expanded prior

# Model for activity
RepFID_male <- MCMCglmm(FID ~ 1 +
  scale(Time) + scale(SeptDay) + scale(Age),
random = ~BirdID,
prior = prior1,
data = FID_data_m,
nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
family = "gaussian"
)

save(RepFID_male, file="~/.../.../RepFID_male.Rdata")

summary(RepFID_male)

# Model diagnostics
plot(RepFID_male$VCV)
plot(RepFID_male$Sol)

# Extract among-individual variance and credible intervals
VarAmong_FID_m <- mean(RepFID_male$VCV[, 1])
VarAmong_FID_CI_m <- HPDinterval(RepFID_male$VCV[, 1])

# Extract within-individual variance and credible intervals
VarWithin_FID_m  <- mean(RepFID_male$VCV[, 2])
VarWithin_FID_CI_m <- HPDinterval(RepFID_male$VCV[, 2])

# Calculate repeatability for FID
RepeatFID_m  <- mean(RepFID_male$VCV[, 1] / (RepFID_male$VCV[, 1] + RepFID_male$VCV[, 2]))
# [1] 0.6018936

# 95% Credible intervals
RepeatFID_CI_m <- HPDinterval(RepFID_male$VCV[, 1] / (RepFID_male$VCV[, 1] + RepFID_male$VCV[, 2]))
# lower     upper
# var1 0.4986549 0.7139631
# attr(,"Probability")
# [1] 0.95


### 2.2 Female model  ======================

FID_data_f <- FID_data  %>% filter(Sex == 'Female')
FID_data_f <- droplevels(FID_data_f)

hist(FID_data_f$FID, breaks = 50)

# Number of females 
FID_data_f %>%
  distinct(BirdID) %>%
  tally()
# 45

# Observations per Season
FID_data_f %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012    59
# 2           2013    84
# 3           2014    105

# Get the number of males for each Season
FID_data_f %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    21
# 2           2013    28
# 3           2014    31

# Mean measurments per male
mean(tapply(FID_data_f$BirdID, FID_data_f$BirdID, length))
# [1] 5.511111

# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) # parameter expanded prior

# Model for activity
RepFID_female <- MCMCglmm(FID ~ 1 +
  scale(Time) + scale(SeptDay) + scale(Age),
random = ~BirdID,
prior = prior1,
data = FID_data_f,
nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
family = "gaussian"
)

save(RepFID_female, file="~/.../.../RepFID_female.Rdata")

summary(RepFID_female)

# Model doagnistics
plot(RepFID_female$VCV)
plot(RepFID_female$Sol)

# Extract among-individual variance and credible intervals
VarAmong_FID_f <- mean(RepFID_female$VCV[, 1])
VarAmong_FID_m_CI <- HPDinterval(RepFID_female$VCV[, 1])

# Extract within-individual variance and credible intervals
VarWithin_FID_f  <- mean(RepFID_female$VCV[, 2])
VarWithin_FID_CI_f <- HPDinterval(RepFID_female$VCV[, 2])

# Calculate repeatability for FID
RepeatFID_f  <- mean(RepFID_female$VCV[, 1] / (RepFID_female$VCV[, 1] + RepFID_female$VCV[, 2]))
# [1] 0.6284872

# 95% Credible intervals
RepeatFID_CI_f <- HPDinterval(RepFID_female$VCV[, 1] / (RepFID_female$VCV[, 1] + RepFID_female$VCV[, 2]))
# lower     upper
# var1 0.5128581 0.7608908
# attr(,"Probability")
# [1] 0.95

# Calculate coefficient of variation (CVb)
CVb_FID_f <- mean(sqrt(RepFID_female$VCV[, 1]) / RepFID_female$Sol[, 1]) # among-individual variance divided by intercept
# [1] 0.3658056



### 3. Repeatability of vigilance  ======================

# Overview of vigilance measurements ----
# Total number of measurments for males and females
summary(vigilance_data$Sex)
# Female   Male   NA's 
#    219    381      6 

# Mean seconds per observations
mean(vigilance_data$Seconds) # 87.11881
sd(vigilance_data$Seconds) # 31.80934

# Mean hops per observation
mean(vigilance_data$Pecks) # 38.26238
sd(vigilance_data$Pecks) # 24.5387

# Observations per Season
vigilance_data %>% group_by(BreedingSeason) %>% tally()
# 1           2012   127
# 2           2013   210
# 3           2014   269

# Get the number of individuals for each Season
vigilance_data %>% group_by(BreedingSeason) %>% distinct(BirdID) %>% tally()
# 1           2012    43
# 2           2013    72
# 3           2014    79

# Mean measurments per bird
mean(tapply(vigilance_data$BirdID, vigilance_data$BirdID, length))
# [1] 5.135593


# Change date into number of days starting from 1st September
# Specify that Date is a date in R
vigilance_data$Date <- as.character(as.Date(vigilance_data$Date, format = "%d/%m/%Y"))

# Reference dates for each season
RefDate2012 <- as.Date("01/09/12", format = "%d/%m/%y")
RefDate2013 <- as.Date("01/09/13", format = "%d/%m/%y")
RefDate2014 <- as.Date("01/09/14", format = "%d/%m/%y")

# Calculate SeptDay for each season
vigilance_data <- vigilance_data %>%
  mutate(SeptDay = ifelse(BreedingSeason == 2012, difftime(vigilance_data$Date, RefDate2012, units = c("days")),
    ifelse(BreedingSeason == 2013, difftime(vigilance_data$Date, RefDate2013, units = c("days")),
      difftime(vigilance_data$Date, RefDate2014, units = c("days"))
    )
  ))


# Remove NAs in fixed effects, necessary for models!
vigilance_data <- vigilance_data %>% drop_na(c(Seconds, Time, SeptDay, Age))


### 3.1 Male model  ======================

vigilance_data_m <- vigilance_data %>% filter(Sex == "Male")
vigilance_data_m <- droplevels(vigilance_data_m)

hist(vigilance_data_m$Pecks, breaks = 50)

# Number of males
vigilance_data_m %>%
  distinct(BirdID) %>%
  tally()
# 73

# Observations per Season
vigilance_data_m %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   92
# 2           2013   131
# 3           2014   158

# Get the number of males for each Season
vigilance_data_m %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    31
# 2           2013    44
# 3           2014    48

# Mean measurments per male
mean(tapply(vigilance_data_m$BirdID, vigilance_data_m$BirdID, length))
# [1] 5.219178

# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) # parameter expanded prior

# Model for activity
RepVig_male <- MCMCglmm(Pecks ~ 1 +
  scale(DecimalMin) + scale(Time) + scale(SeptDay) + scale(Age),
random = ~BirdID,
prior = prior1,
data = vigilance_data_m,
nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
family = "gaussian"
)

save(RepVig_male, file="~/.../.../RepVig_male.Rdata")

summary(RepVig_male)

# Model diagnistics
plot(RepVig_male$VCV)
plot(RepVig_male$Sol)

# Extract among-individual variance and credible intervals
VarAmong_Vig_m <- mean(RepVig_male$VCV[, 1])
VarAmong_Vig_CI_m <- HPDinterval(RepVig_male$VCV[, 1])

# Extract within-individual variance and credible intervals
VarWithin_Vig_m <- mean(RepVig_male$VCV[, 2])
VarWithin_Vig_CI_m <- HPDinterval(RepVig_male$VCV[, 2])

# Calculate repeatability for vigilance
RepeatVig_m <- mean(RepVig_male$VCV[, 1] / (RepVig_male$VCV[, 1] + RepVig_male$VCV[, 2]))
# [1] 0.03535229

# 95% Credible intervals
RepeatVig_CI_m <- HPDinterval(RepVig_male$VCV[, 1] / (RepVig_male$VCV[, 1] + RepVig_male$VCV[, 2]))
# lower     upper
# var1 4.190274e-08 0.09695781
# attr(,"Probability")
# [1] 0.95


### 3.2 Female model  ======================

vigilance_data_f <- vigilance_data %>% filter(Sex == "Female")
vigilance_data_f <- droplevels(vigilance_data_f)

hist(vigilance_data_f$Pecks, breaks = 50)

# Number of males
vigilance_data_f %>%
  distinct(BirdID) %>%
  tally()
# 42

# Observations per Season
vigilance_data_f %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   35
# 2           2013   73
# 3           2014   111

# Get the number of males for each Season
vigilance_data_f %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    12
# 2           2013    25
# 3           2014    31

# Mean measurments per male
mean(tapply(vigilance_data_f$BirdID, vigilance_data_f$BirdID, length))
# [1] 5.214286


# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) # parameter expanded prior

# Model for activity
RepVig_female <- MCMCglmm(Pecks ~ 1 +
  scale(DecimalMin) + scale(Time) + scale(SeptDay) + scale(Age),
random = ~BirdID,
prior = prior1,
data = vigilance_data_f,
nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
family = "gaussian"
)

save(RepVig_female, file="~/.../.../RepVig_female.Rdata")

summary(RepVig_female)

# Model diagnistics
plot(RepVig_female$VCV)
plot(RepVig_female$Sol)

# Extract among-individual variance and credible intervals
VarAmong_Vig_f <- mean(RepVig_female$VCV[, 1])
VarAmong_Vig_CI_f <- HPDinterval(RepVig_female$VCV[, 1])

# Extract within-individual variance and credible intervals
VarWithin_Vig_f <- mean(RepVig_female$VCV[, 2])
VarWithin_Vig_CI_f <- HPDinterval(RepVig_female$VCV[, 2])

# Calculate repeatability for vigilance
RepeatVig_f <- mean(RepVig_female$VCV[, 1] / (RepVig_female$VCV[, 1] + RepVig_female$VCV[, 2]))
# [1] 0.06856051

# 95% Credible intervals
RepeatVig_CI_f <- HPDinterval(RepVig_female$VCV[, 1] / (RepVig_female$VCV[, 1] + RepVig_female$VCV[, 2]))
# lower     upper
# var1 1.481632e-07 0.1864457
# attr(,"Probability")
# [1] 0.95


### 4. Repeatability of provisioning  ======================

# Overview of provisioning measurements ----
# Total number of measurments for males and females
summary(provisioning_data$Sex)
# Female   Male 
# 345    531 

# Observations per Season
provisioning_data %>% group_by(BreedingSeason) %>% tally()
# 1           2012   171
# 2           2013   377
# 3           2014   328

# Get the number of individuals for each Season
provisioning_data %>% group_by(BreedingSeason) %>% distinct(BirdID) %>% tally()
# 1           2012    38
# 2           2013    57
# 3           2014    52

# Get the number of nests for each Season
provisioning_data %>% group_by(BreedingSeason) %>% distinct(NestID) %>% tally()
# 1           2012    20
# 2           2013    33
# 3           2014    29

# Get the number of videos for each Season
provisioning_data %>%
  arrange(NestID) %>%
  group_by(NestID, Date, BreedingSeason) %>%
  mutate(SequenceNest = row_number()) %>%
  filter(SequenceNest == 1) %>%
  ungroup(mating_systems_groups) %>%
  tally()
# 1   345

# Mean measurments per bird
mean(tapply(provisioning_data$BirdID, provisioning_data$BirdID, length))
# [1] 9.733333

# Change date into number of days starting from 1st September
# Specify that Date is a date in R
provisioning_data$Date <- as.character(as.Date(provisioning_data$Date, format = "%d/%m/%Y"))

# Reference dates for each season
RefDate2012 <- as.Date("01/09/12", format = "%d/%m/%y")
RefDate2013 <- as.Date("01/09/13", format = "%d/%m/%y")
RefDate2014 <- as.Date("01/09/14", format = "%d/%m/%y")

# Calculate SeptDay for each season
provisioning_data <- provisioning_data %>%
  mutate(SeptDay = ifelse(BreedingSeason == 2012, difftime(provisioning_data$Date, RefDate2012, units = c("days")),
    ifelse(BreedingSeason == 2013, difftime(provisioning_data$Date, RefDate2013, units = c("days")),
      difftime(provisioning_data$Date, RefDate2014, units = c("days"))
    )
  ))

# Remove NAs in fixed effects, necessary for models!
provisioning_data <- provisioning_data %>% drop_na(c(Min, SeptDay, Age))


### 4.1 Male model  ======================

provisioning_data_m <- provisioning_data %>% filter(Sex == "Male")
provisioning_data_m <- droplevels(provisioning_data_m)

hist(provisioning_data_m$Visits, breaks = 50) 

# Number of males
provisioning_data_m %>%
  distinct(BirdID) %>%
  tally()
# 52

# Observations per Season
provisioning_data_m %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   103
# 2           2013   224
# 3           2014   204

# Get the number of males for each Season
provisioning_data_m %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    22
# 2           2013    33
# 3           2014    31

# Mean measurments per male
mean(tapply(provisioning_data_m$BirdID, provisioning_data_m$BirdID, length))
# [1] 10.21154


# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) # parameter expanded prior

# Model for activity
RepProv_male <- MCMCglmm(Visits ~ 1 +
  scale(Min) + scale(SeptDay) + scale(Age),
random = ~BirdID,
prior = prior1,
data = provisioning_data_m,
nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
family = "gaussian"
)

save(RepProv_male, file="~/.../.../RepProv_male.Rdata")

summary(RepProv_male)

# Model diagnistics
plot(RepProv_male$VCV)
plot(RepProv_male$Sol)

# Extract among-individual variance and credible intervals
VarAmong_Prov_m <- mean(RepProv_male$VCV[, 1])
VarAmong_Prov_CI_m <- HPDinterval(RepProv_male$VCV[, 1])

# Extract within-individual variance and credible intervals
VarWithin_Prov_m <- mean(RepProv_male$VCV[, 2])
VarWithin_Prov_CI_m <- HPDinterval(RepProv_male$VCV[, 2])

# Calculate repeatability for provisioning
RepeatProv_m <- mean(RepProv_male$VCV[, 1] / (RepProv_male$VCV[, 1] + RepProv_male$VCV[, 2]))
# [1] 0.3969362

# 95% Credible intervals
RepeatProv_CI_m <- HPDinterval(RepProv_male$VCV[, 1] / (RepProv_male$VCV[, 1] + RepProv_male$VCV[, 2]))
# lower     upper
# var1 0.27087 0.5213708
# attr(,"Probability")
# [1] 0.95


### 4.2 Female model  ======================

provisioning_data_f <- provisioning_data %>% filter(Sex == "Female")
provisioning_data_f <- droplevels(provisioning_data_f)

hist(provisioning_data_f$Visits, breaks = 50

# Number of males
provisioning_data_f %>%
  distinct(BirdID) %>%
  tally()
# 38

# Observations per Season
provisioning_data_f %>%
  group_by(BreedingSeason) %>%
  tally()
# 1           2012   68
# 2           2013   153
# 3           2014   124

# Get the number of males for each Season
provisioning_data_f %>%
  group_by(BreedingSeason) %>%
  distinct(BirdID) %>%
  tally()
# 1           2012    16
# 2           2013    24
# 3           2014    21

# Mean measurments per male
mean(tapply(provisioning_data_f$BirdID, provisioning_data_f$BirdID, length))
# [1] 9.078947


# MCMCglmm ----------------------------------------------------------------

# Specify priors
prior1 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000))) # parameter expanded prior

# Model for activity
RepProv_female <- MCMCglmm(Visits ~ 1 +
  scale(Min) + scale(SeptDay) + scale(Age), 
random = ~BirdID,
prior = prior1,
data = provisioning_data_f,
nitt = 13000 * 50, burnin = 3000 * 50, thin = 10 * 25,
family = "gaussian"
)

save(RepProv_female, file="~/.../.../RepProv_female.Rdata")

summary(RepProv_female)

# Model diagnistics
plot(RepProv_female$VCV)
plot(RepProv_female$Sol)

# Extract among-individual variance and credible intervals
VarAmong_Prov_f <- mean(RepProv_female$VCV[, 1])
VarAmong_Prov_CI_f <- HPDinterval(RepProv_female$VCV[, 1])

# Extract within-individual variance and credible intervals
VarWithin_Prov_f <- mean(RepProv_female$VCV[, 2])
VarWithin_Prov_CI_f <- HPDinterval(RepProv_female$VCV[, 2])

# Calculate repeatability for provisioning
RepeatProv_f <- mean(RepProv_female$VCV[, 1] / (RepProv_female$VCV[, 1] + RepProv_female$VCV[, 2]))
# [1] 0.3424436

# 95% Credible intervals
RepeatProv_CI_f <- HPDinterval(RepProv_female$VCV[, 1] / (RepProv_female$VCV[, 1] + RepProv_female$VCV[, 2]))
# lower     upper
# var1 0.1958816 0.494044
# attr(,"Probability")
# [1] 0.95
