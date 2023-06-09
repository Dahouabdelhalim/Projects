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
biomass <- biomass %>% rename("ingrowth"="Annual.Total.Root.Ingrowth..g.m.2.",
                              "total"="Total.Biomass",)

#restrict to ambient CO2 treatment and remove rows with N/A values in the ingrowth column
ambCO2 <- biomass %>% filter(CO2.Treatment == "Camb")
ambCO2 <- ambCO2 %>% drop_na(ingrowth)

#############################
### Exotic Species Group ####
#############################

#generate "exotics" sub-dataset from monospecies plots of Elymus (Agropyron) repens and Poa pratensis
exotics <- ambCO2 %>% filter(monospecies == "Agropyron repens" | monospecies == "Poa pratensis")

#estimate u_e 
#use root ingrowth as proxy for root senesence, assuming the two are in balance
#measurements of both root ingrowth and root mass are in the top 0-20cm of soil 
exotics <- exotics %>% mutate (ue = ingrowth / total) #divide yearly root ingrowth by total plant biomass, yielding value with units of of g g^(-1) yr^ = yr^(-1)

#Based on the following boxplots, we select the value ue=0.2 yr^(-1) 
#as reasonable for the exotic belowground tissue senesence parameter
ggplot(data = exotics, aes(y = ue, x = monospecies)) + geom_boxplot() 


#############################
### Native Species Group ####
#############################

#generate "natives" sub-dataset from monospecies plots of Schizachyrium scoparium and 16-species plots
#(the 16-species plots may contain some non-native species, but feature biodiversity, a linked attribute in our model)
natives <- ambCO2 %>% filter(monospecies == "Schizachyrium scoparium" | CountOfSpecies == "16")

#estimate u_n 
#use root ingrowth as proxy for root senesence, assuming the two are in balance
#measurements of both root ingrowth and root mass are in the top 0-20cm of soil 
natives <- natives %>% mutate (un = ingrowth / total) #divide yearly root ingrowth by total plant biomass, yielding value with units of of g g^(-1) yr^ = yr^(-1)

#Based on the following boxplots, we select the value un=0.2 yr^(-1) 
#as reasonable for the native belowground tissue senesence parameter
ggplot(data = natives, aes(y = un, x = monospecies)) + geom_boxplot() 


