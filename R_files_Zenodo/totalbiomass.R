rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(ggplot2)
library(tidyr)

#################################
#### Initial Data Processing ####
#################################

#import biomass data
biomass <- read.csv("e141biomassAboveBelowLitterRoot.csv")

#shorten variable names
biomass <- biomass %>% rename("total"="Total.Biomass")

#restrict to ambient CO2 treatment and remove rows with N/A values in the ingrowth column
ambCO2 <- biomass %>% filter(CO2.Treatment == "Camb")

#############################
### Exotic Species Group ####
#############################

#generate "exotics" sub-dataset from monospecies plots of Elymus (Agropyron) repens and Poa pratensis
exotics <- ambCO2 %>% filter(monospecies == "Agropyron repens" | monospecies == "Poa pratensis")

#The following boxplots illustrate order of magnitude of biomass density (g/m^2) for exotics Elymus (Agropyron) repens and Poa pratensis
ggplot(data = exotics, aes(y = total, x = monospecies)) + geom_boxplot() 

#############################
### Native Species Group ####
#############################

#generate "natives" sub-dataset from monospecies plots of Schizachyrium scoparium and 16-species plots
#(the 16-species plots may contain some non-native species, but feature biodiversity, a linked attribute in our model)
natives <- ambCO2 %>% filter(monospecies == "Schizachyrium scoparium" | CountOfSpecies == "16")

#The following boxplots illustrate order of magnitude of biomass density (g/m^2) for 16 species plots (no label) and native Schizachyrium scoparium
ggplot(data = natives, aes(y = total, x = monospecies)) + geom_boxplot() 