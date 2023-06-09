# Manuscript: "Anthropogenic noise impairs cooperation in bottlenose dolphins"
# by SÃ¸rensen, Haddock, Guarino, Jaakkola, McMullen, Jensen, Tyack and King
#
# Script for investigating how the use of echolocation changes depending on the 
# noise levels (energy flux density, dB re 1 uPa^2*s) the dolphins were exposed
# to.
#
# Last updated: 14-11-2022
################################################################################
# clear environment
rm(list = ls())

# Load packages
library(DHARMa) 
library(lme4)
library(glmmTMB)  
library(lmerTest)

# Load data 
clickdata <- read.csv("Data_echolocation.csv", header=T) 
attach(clickdata)
str(clickdata)

# Subset data as per dolphin ID
Reese <- clickdata[which (clickdata$animal == 1), ]
Delta <- clickdata[which (clickdata$animal == 2), ]

################################################################################
#### Reese ####
# Run LMM with session number (1-20) and day (1-8) as random effects:
Reese$proportion_new <- Reese$proportion_clicking/100 
model1_full <- lmer(proportion_new ~ EFD + (1|session) + (1|day), data = Reese) 
model1_null <- lmer(proportion_new ~ 1 + (1|session) + (1|day), data = Reese)

# Model results
anova(model1_null, model1_full)
summary(model1_full) 
confint(model1_full)

# Modelfit using DHARMa package: 
simulationOutput<-simulateResiduals(fittedModel = model1_full, plot=T)

# Diagnostic plots:
par(mfrow=c(1,3))
hist(residuals(model1_full)) ; 
qqnorm(residuals(model1_full)) ; qqline(residuals(model1_full)); 
plot(residuals(model1_full) ~ fitted(model1_full))

################################################################################
#### Delta ####
# Run LMM with session number (1-20) and day (1-8) as random effects:
Delta$proportion_new <- Delta$proportion_clicking/100
model2_full <- lmer(proportion_new ~ EFD + (1|session) + (1|day), data = Delta)
model2_null <- lmer(proportion_new~ 1 + (1|session) + (1|day), data = Delta)

# Model results
anova(model2_null, model2_full)
summary(model2_full) 
confint(model2_full)

# Modelfit using DHARMa package: 
simulationOutput<-simulateResiduals(fittedModel = model2_full, plot=T) 

# Diagnostic plots:
par(mfrow=c(1,3))
hist(residuals(model2_full)) ; 
qqnorm(residuals(model2_full)) ; qqline(residuals(model2_full)); 
plot(residuals(model2_full) ~ fitted(model2_full))

################################# END-OF-SCRIPT ################################

