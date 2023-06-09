#### Small mammals reduce distance-dependence and increase seed predation risk in tropical rainforest fragments ####
#### Biotropica
#### Authors: Aparna KRISHNAN*, Anand M OSURI, Meghna KRISHNADAS

#### Description: This file contains the R codes for Table 3-- the statistical analysis of  seed survival using survivial analysis 

# Setting up the libraries and data sheet ####
library(dplyr)
library(survival)
library(coxme)
library(ggplot2)

#### TABLE 3 ####
seeds <- read.csv("seed_survival_data.csv", header = T, stringsAsFactors = F)

# Assigning the control treatments: Contiguous forests and Near
seeds$Plot_Location <- factor(seeds$Plot_Location, levels = c("N","F"))
seeds$Landscape_Type <- factor(seeds$Landscape_Type, levels = c("Cont","Frag"))

# Unique ID for seeds at each Plot and each Tree within the experimental block
seeds$ID <- paste(seeds$Species_Name, seeds$Location_Name, seeds$Block_No, seeds$Plot_Location)
seeds$DistID <-  paste(seeds$Species_Name, seeds$Location_Name, seeds$Block_No)

#### Survival Analysis ####
# Cox Proportional Hazard Mixed Effects Model 

## Acronychia pedunculata ####
ap.seeds <- seeds[seeds$Species_Name == "AP",]

# Open Plots (No Exclosure)
ap.seeds.nex <- ap.seeds[ap.seeds$Treatment == "No_Exclosure",]
# Plotting survival as a function of landscpe and plot location
# Nested Random effects- All seeds at a focal tree : All seeds within a plot
mod <- coxme(Surv(Weeks_to_Event, Event) ~ Landscape_Type * Plot_Location + (1|DistID/ID), data = ap.seeds.nex)
# Model Results
mod
# Exponent(confidence intervals)
exp(1)^confint(mod)

# Exclosure Plots
ap.seeds.ex <- ap.seeds[ap.seeds$Treatment == "Exclosure",]
mod <- coxme(Surv(Weeks_to_Event, Event) ~ Landscape_Type * Plot_Location + (1|DistID/ID), data = ap.seeds.ex)
mod
exp(1)^confint(mod)

## Cullenia exarillata ####
ce.seeds <- seeds[seeds$Species_Name == "CE",]

# Open Plots (No Exclosure)
ce.seeds.nex <- ce.seeds[ce.seeds$Treatment == "No_Exclosure",]
mod <- coxme(Surv(Weeks_to_Event, Event) ~ Landscape_Type * Plot_Location + (1|DistID/ID), data = ce.seeds.nex)
mod
exp(1)^confint(mod)

# Exclosure Plots
ce.seeds.ex <- ce.seeds[ce.seeds$Treatment == "Exclosure",]
mod <- coxme(Surv(Weeks_to_Event, Event) ~ Landscape_Type * Plot_Location + (1|DistID/ID), data = ce.seeds.ex)
mod
exp(1)^confint(mod)

## Syzygium rubicundum ####
sr.seeds <- seeds[seeds$Species_Name == "SR",]

# Open Plots (No Exclosure)
sr.seeds.nex <- sr.seeds[sr.seeds$Treatment == "No_Exclosure",]
mod <- coxme(Surv(Weeks_to_Event, Event) ~ Landscape_Type * Plot_Location + (1|DistID/ID), data = sr.seeds.nex)
mod
exp(1)^confint(mod)

# Exclosure Plots
sr.seeds.ex <- sr.seeds[sr.seeds$Treatment == "Exclosure",]
mod <- coxme(Surv(Weeks_to_Event, Event) ~ Landscape_Type * Plot_Location + (1|DistID/ID), data = sr.seeds.ex)
mod
exp(1)^confint(mod)

## Ormosia travancorica ####
ot.seeds <- seeds[seeds$Species_Name == "OT",]

# Open Plots (No Exclosure)
mod <- coxme(Surv(Weeks_to_Event, Event) ~ Landscape_Type * Plot_Location + (1|DistID/ID), data = ot.seeds)
mod
exp(1)^confint(mod)
