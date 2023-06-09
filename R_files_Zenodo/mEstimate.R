rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(ggplot2)

#################################
#### Initial Data Processing ####
#################################

#import biomass data
biomass <- read.csv("e141biomassAboveBelowLitterRoot.csv")

#shorten variable names
biomass <- biomass %>% rename("above" ="AbovegroundTotal.Biomass..g.m.2.",
                              "below"="Total.root.biomass.0.20..g.m.2.",
                              "total"="Total.Biomass",
                              "litter"="Miscellaneous.Litter..g.m.2.")

#restrict to ambient CO2 treatment
ambCO2 <- biomass %>% filter(CO2.Treatment == "Camb")

#restrict to pairs of years with a spring burn followed by no burn 
# (these turn out to be even followed by odd years)
#ASSUMPTIONS:
## the spring burn resets litter levels to near zero
## August measurements of aboveground biomass give an upper estimate of litter production that year
## June measurements of litter the following year give a lower estimate of litter production

ambCO2burn <- ambCO2 %>% filter((year == "2003" & Season == "August") | 
                                  (year == "2004" & Season == "June") | 
                                  (year == "2005" & Season == "August") | 
                                  (year == "2006" & Season == "June") |
                                  (year == "2007" & Season == "August") | 
                                  (year == "2008" & Season == "June") |
                                  (year == "2009" & Season == "August") |
                                  (year == "2010" & Season == "June"))
#############################
### Exotic Species Group ####
#############################

#generate "exotics" sub-dataset from monospecies plots of Elymus (Agropyron) repens and Poa pratensis
exotics <- ambCO2burn %>% filter(monospecies == "Agropyron repens" | monospecies == "Poa pratensis")

#join plot data from adjacent years (Aug 2003 with June 2004, etc.)
exotics_Aug <- exotics %>% filter(Season == "August") %>% select(year, Season, Plot, monospecies, above, below, total) %>% mutate(augYear=year)
exotics_Jun <- exotics %>% filter(Season == "June") %>% select(year, Season, Plot, litter) %>% mutate(augYear=year-1)
exotics <- full_join(exotics_Aug, exotics_Jun, by=c("augYear", "Plot"))

#calculate upper and lower estimates of the parameter m
exotics <- exotics %>% mutate (mUpper = above / total) #August aboveground biomass per total biomass
exotics <- exotics %>% mutate (mLower = litter / total) #following June litter per previous year total biomass
#note: root biomass was estimated from samples 20cm deep, making "total" a likely underestimate and inflating mUpper, mLower.

#Based on the following boxplots, we select the value m=0.2 yr^(-1) 
#as reasonable for the exotic litter production parameter
ggplot(data = exotics, aes(y = mUpper, x = monospecies)) + geom_boxplot() 
ggplot(data = exotics, aes(y = mLower, x = monospecies)) + geom_boxplot()

#############################
### Native Species Group ####
#############################

#generate "natives" sub-dataset from monospecies plots of Schizachyrium scoparium and 16-species plots
#(the 16-species plots may contain some non-native species, but feature biodiversity, a linked attribute in our model)
natives <- ambCO2burn %>% filter(monospecies == "Schizachyrium scoparium" | CountOfSpecies == "16")

#join plot data from adjacent years (Aug 2003 with June 2004, etc.)
natives_Aug <- natives %>% filter(Season == "August") %>% select(year, Season, Plot, monospecies, above, below, total) %>% mutate(augYear=year)
natives_Jun <- natives %>% filter(Season == "June") %>% select(year, Season, Plot, litter) %>% mutate(augYear=year-1)
natives <- full_join(natives_Aug, natives_Jun, by=c("augYear", "Plot"))

#calculate upper and lower estimates of the parameter m
natives <- natives %>% mutate (mUpper = above / total) #August aboveground biomass per total biomass
natives <- natives %>% mutate (mLower = litter / total) #following June litter per previous year total biomass
#note: root biomass was estimated from samples 20cm deep, making "total" a likely underestimate and inflating mUpper, mLower.

#Based on the following boxplots, we select the value m=0.2 yr^(-1) 
#as a slightly low but still plausible value for the native litter production parameter
ggplot(data = natives, aes(y = mUpper, x = monospecies)) + geom_boxplot() 
ggplot(data = natives, aes(y = mLower, x = monospecies)) + geom_boxplot()
