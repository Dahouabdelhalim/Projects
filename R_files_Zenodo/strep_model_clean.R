### Statistical analyses for the Ammophila-Paraxenos project (A Double-Edged Sword: Parental care increases risk of offspring infection by a maternally-vectored parasite)
## Authors: Rebecca Jean Millena, Jay Rosenheim

## Load in the whole dataset
library(tidyverse)
parasitism_lm <- read_csv("ammophila_strep_project.csv")

parasitism_lm

## Set each of the categorical effects of the model as factors
parasitism_lm$species_id <- as.factor(parasitism_lm$species_id)
parasitism_lm$county <- as.factor(parasitism_lm$county)
parasitism_lm$month <- as.factor(parasitism_lm$month)

## Running a generalized linear model with the base stats package in R; parasitism rate is the response, prey provisioning, month, species identification, and specimen size are all included as fixed effects
parasitism_log <- glm(parasitism ~ mean_prey_provisioned + month + species_id + size, data = parasitism_lm, 
                      family = binomial())

print(parasitism_log)
summary(parasitism_log)

## Running a generalized linear mixed-effects model with package "lme4"
library(lme4)

## This model iteration includes spatial variation (in the form of the "county" entries) as a random effect
parasitism_lme <- glmer(parasitism ~ mean_prey_provisioned + month + size + (1 | county) + (1 | species_id), data = parasitism_lm,
                        family = binomial())
## With random effects, we don't know what effect they might have on the parasitism rates

summary(parasitism_lme)

## Running ANOVA with the "car" package to get overall effects for all variables in the model
library(car)
ano_para_lme <- Anova(parasitism_lme)
ano_para_lme
summary(ano_para_lme)

## Checking for collinearity between variables
library(faraway)
pairs(parasitism ~ mean_prey_provisioned + month + size + county + species_id, data = parasitism_lm)
cor(parasitism_lm$size, parasitism_lm$mean_prey_provisioned)
