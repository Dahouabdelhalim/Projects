# Manuscript: "Anthropogenic noise impairs cooperation in bottlenose dolphins"
# by SÃ¸rensen, Haddock, Guarino, Jaakkola, McMullen, Jensen, Tyack and King

# Last updated: 14-11-2022

###################### Model estimates - behavioural data ######################
# Modelling on behaviours marked by PS with at least 'good agreement' between all 
# coders in the IRR analysis: 
# (1) Looking towards their partner
# (2) Visiting the other half of the lagoon 
# (3) Swimming together
################################################################################
## Load data and packages
################################################################################
rm(list = ls())

# load packages
library(glmmTMB) 
library(DHARMa) 
library(ggdist)
library(lmerTest)
library(nlme)
library(lme4)

## Load data
# The behavioural data are all the marked data done by PS that show good agreement between all three 
# observers in a IRR analysis: 
behav_data <- read.csv("Data_behaviour.csv", header=T)
attach(behav_data)
head(behav_data)
nrow(behav_data)

############################# Behaviour: "Looking" #############################
# Looking at their partner (binary) ~ EFD_combined (noise level per trial quantified as energy flux density, dB re 1 uPa^2*s) 
# with session number (1-20) and day (1-8) included as random effects 

full_look_occurrence <- glmer(Look_combined_binary ~ EFD_combined + (1|session) + (1|day), family = binomial, data = behav_data)
null_look_occurrence <- glmer(Look_combined_binary ~ 1 + (1|session) + (1|day), family = binomial, data = behav_data)

simulationOutput<-simulateResiduals(fittedModel = full_look_occurrence, plot=T)

anova(full_look_occurrence, null_look_occurrence) 
summary(full_look_occurrence)  
confint(full_look_occurrence)

############################ Other half ~treatment ######################################
# Visiting the other half of the lagoon  (binary) ~ EFD_combined (noise level per trial quantified as energy flux density, dB re 1 uPa^2*s) 
# with session number (1-20) and day (1-8) included as random effects 

full_otherhalf_occurrence <- glmer(Other_half_combined_binary ~ EFD_combined + (1|session) + (1|day), family = binomial, data = behav_data, nAGQ = 0) 
null_otherhalf_occurrence <- glmer(Other_half_combined_binary ~ 1 + (1|session) + (1|day), family = binomial, data = behav_data, nAGQ = 0)

simulationOutput<-simulateResiduals(fittedModel = full_otherhalf_occurrence, plot=T)

anova(full_otherhalf_occurrence, null_otherhalf_occurrence)
summary(full_otherhalf_occurrence)
confint(full_otherhalf_occurrence, method = "Wald") 

########################### Behaviour: Swim together ###########################
# swimming together (binary) ~ EFD_combined (noise level per trial quantified as energy flux density, dB re 1 uPa^2*s) 
# with session number (1-20) and day (1-8) included as random effects 

full_swim_occurrence <- glmer(Swim_together_binary ~ EFD_combined + (1|session) + (1|day), family = binomial, data = behav_data)
null_swim_occurrence <- glmer(Swim_together_binary~1 + (1|session) + (1|day), family = binomial, data = behav_data)

simulationOutput<-simulateResiduals(fittedModel = full_swim_occurrence, plot=T) 

anova(full_swim_occurrence, null_swim_occurrence)
summary(full_swim_occurrence) 
confint(full_swim_occurrence)

################################ END OF SCRIPT #################################s