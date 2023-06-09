## A mid-20th century estimate of global vegetation carbon stocks based on inventory statistics
## Authors: Bhan et al., Journal of Land Use Science
## Sensitivity Analysis

# We use the best-guess and range of values developed in Publication_UncertaintyAnalysis_Figures.R
# First to estimate sensitivity at the global level and then at the subcontinental for the best-guess case
# All carbon density estimates used here are described in SI Table 6

setwd("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles")

library(dplyr)
library(magrittr)
library(tidyr)
library(readxl)
library(ggplot2)
library(cowplot)
library(scales)
library(ggpubr)

options(scipen = 999)
options(digits = 5)


## Input files----

SubContinentalClassification <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_CountriesListSubContinentalClassification.csv")

LandUseEstimates <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/5_CountriesListLandUseBaselineWithKastnerExtraForest.csv") %>%
  mutate(accessibleforest_Mha = if_else(country == 'Indonesia', 64.4, accessibleforest_Mha),
         inaccessibleforest_Mha = if_else(country == 'Indonesia', 46.4, inaccessibleforest_Mha),
         totalforest_Mha = if_else(country == 'Indonesia', 
                                   accessibleforest_Mha + inaccessibleforest_Mha + extraforest_Mha, 
                                   totalforest_Mha),
         othervegetatedland_Mha = if_else(country == 'Indonesia', 45.49, othervegetatedland_Mha), # till here, British North Borneo correction
         
         accessibleforest_Mha = if_else(country == 'Malaya', 17.58, accessibleforest_Mha),
         inaccessibleforest_Mha = if_else(country == 'Malaya', 7.82, inaccessibleforest_Mha),
         extraforest_Mha = if_else(country == 'Malaya', 4.4, extraforest_Mha),
         totalforest_Mha = if_else(country == 'Malaya', 
                                   accessibleforest_Mha + inaccessibleforest_Mha + extraforest_Mha, 
                                   totalforest_Mha),
         othervegetatedland_Mha = if_else(country == 'Malaya', 0, othervegetatedland_Mha)) # till here, British North Borneo correction


SubcontinentalNPP <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_SubcontinentalListNPPValues.csv") %>%
  dplyr::select(subcontinent, NPPact_tChayr)

FRAForestCarbonDensities <- 
  read.csv('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_CountriesListAllForestCarbonDensities.csv')

PotentialForestCarbonStocks <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/3_PotentialForestCarbonStocks.csv")

#PotentialForestCarbonStocks <- read.csv('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/3_CountriesListCountryMeanPotentialCarbonStocks.csv')

PotentialCarbonStocks <- 
  read.csv('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/3_CountriesListCountryMeanPotentialCarbonStocksForOtherland.csv')

Erb2018CarbonDensities <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_CountriesListErb2018CarbonDensities.csv") %>%
  replace(is.na(.), 0) %>%
  mutate_all(funs(ifelse(.== '#DIV/0!', 0, .))) %>%
  mutate_all(funs(ifelse(.== '', 0, .)))

Xia2014RangelandCS <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/Xia2014GrasslandSubcontinentCS.csv") # measured at peak AGB

LandUseEstimatesWithShiftingCultivation <- 
  read.csv('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/5_CountriesListLandUseWithShiftingCultivationWithKastnerExtraForest.csv') %>%
  mutate(accessibleforest_Mha = if_else(country == 'Indonesia', 64.4, accessibleforest_Mha),
         inaccessibleforest_Mha = if_else(country == 'Indonesia', 46.4, inaccessibleforest_Mha),
         totalforest_Mha = if_else(country == 'Indonesia', 
                                   accessibleforest_Mha + inaccessibleforest_Mha + extraforest_Mha, 
                                   totalforest_Mha),
         othervegetatedland_Mha = if_else(country == 'Indonesia', 45.49, othervegetatedland_Mha), # till here, British North Borneo correction
         othervegetatedland_Mha = if_else(country == 'Indonesia', othervegetatedland_Mha - 2.008 - 2.938, othervegetatedland_Mha), # SC correction
         
         accessibleforest_Mha = if_else(country == 'Malaya', 17.58, accessibleforest_Mha),
         inaccessibleforest_Mha = if_else(country == 'Malaya', 7.82, inaccessibleforest_Mha),
         totalforest_Mha = if_else(country == 'Malaya', 29.8, totalforest_Mha),
         othervegetatedland_Mha = if_else(country == 'Malaya', 0, othervegetatedland_Mha), # till here, British North Borneo correction
         extraforest_Mha = if_else(country == 'Malaya', 3.735, extraforest_Mha),
         totalforest_Mha = if_else(country == 'Malaya', 
                                   accessibleforest_Mha + inaccessibleforest_Mha + extraforest_Mha, 
                                   totalforest_Mha)) # SC correction

SCCarbonDensities <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_CountriesListShiftingCultivationAreasCarbonDensitiesSilva2011.csv") %>%
  mutate(R.Sratio = 0.28,
         BGC_shortfallow_tCha = AGC_shortfallow_tCha * R.Sratio,
         totalC_shortfallow_tCha = AGC_shortfallow_tCha + BGC_shortfallow_tCha,
         BGC_longfallow_tCha = AGC_longfallow_tCha * R.Sratio,
         totalC_longfallow_tCha = AGC_longfallow_tCha + BGC_longfallow_tCha) %>%
  dplyr::select(country, totalC_shortfallow_tCha, totalC_longfallow_tCha) %>%
  filter(!(country == 'British North Borneo' |
             country == 'New Guinea (Neth.)' |
             country == 'Sarawak'))

CarbonDensitiesWhereActGreaterThanPot <- 
  read.csv('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_ForestAreasPotentialActualOfCountriesWithForestCarbonDensitiesActGreaterThanPot.csv') %>%
  dplyr::select(country, nongrowingstockcarbondensity_tCha)

PotentialForestArea <- 
  read_xlsx('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/3_CountriesListPotentialForestAreas_3kgC.xlsx') %>%
  mutate(potentialforest_Mha = as.numeric(potentialforest_Mha)) %>%
  dplyr::select(country, potentialforest_Mha)

NPPGrazing1950 <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_NPPactGrazingLand1950_Kastner2021.csv") %>%
  dplyr::select(subcontinent, NPPactgrazing1950_tCha)

NPPCropland1950 <-
  read_xlsx("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/NPP1950Krausmann2000Haberl.xlsx") %>%
  dplyr::select(subcontinent, totalNPPact1950_tCha)

## GLOBAL----

## Accessible forest with SC----

AccessibleForestWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Select accessible forest
  dplyr::select(country, accessibleforest_Mha, subcontinent, climatezone) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha, forestcarbondensityallFAO_tCha, forestcarbondensitytype), 
            by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  
  
  # Forest CS
  mutate(accessibleforestcarbonstocksactual_MtC = (accessibleforest_Mha * forestcarbondensity_tCha),
         accessibleforestcarbonstocksactual_PgC = accessibleforestcarbonstocksactual_MtC/1000,
         
         # all FAO
         accessibleforestcarbonstocksallFAO_MtC = (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
         accessibleforestcarbonstocksallFAO_PgC = accessibleforestcarbonstocksallFAO_MtC/1000,
         
         # Upper range
         accessibleforestcarbonstocksallFAOupperrange_MtC = 
           case_when(climatezone == 'temperate' & forestcarbondensitytype == 'FAO-derived' ~ 
                       1.25 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'FAO-derived' ~ 
                       1.25 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'Continental average' ~ 
                       1.4 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'Independent studies' ~ 
                       1.25 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'temperate' & forestcarbondensitytype == 'Independent studies' ~ 
                       1.25 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'temperate' & forestcarbondensitytype == 'Continental average' ~ 
                       1.4 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'FAO-derived' ~ 
                       1.4 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'Independent studies' ~ 
                       1.4 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'Continental average' ~ 
                       1.5 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha)),
         
         accessibleforestcarbonstocksallFAOupperrange_PgC = accessibleforestcarbonstocksallFAOupperrange_MtC/1000,
         
         # Lower range
         accessibleforestcarbonstocksallFAOlowerrange_MtC = 
           case_when(climatezone == 'temperate' & forestcarbondensitytype == 'FAO-derived' ~ 
                       0.75 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'FAO-derived' ~ 
                       0.75 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'Continental average' ~ 
                       0.6 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'Independent studies' ~ 
                       0.75 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'temperate' & forestcarbondensitytype == 'Independent studies' ~ 
                       0.75 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'temperate' & forestcarbondensitytype == 'Continental average' ~ 
                       0.6 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'FAO-derived' ~ 
                       0.6 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'Independent studies' ~ 
                       0.6 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'Continental average' ~ 
                       0.5 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha)),
         
         accessibleforestcarbonstocksallFAOlowerrange_PgC = accessibleforestcarbonstocksallFAOlowerrange_MtC/1000,
         
         #accessibleforestcarbonstockspotential_MtC = (accessibleforest_Mha * potentialforestcarbon_tCha),
         #accessibleforestcarbonstockspotential_PgC = accessibleforestcarbonstockspotential_MtC/1000
         
  ) %>%
  
  group_by(subcontinent) %>%
  summarise(
    carbonstocksactual_PgC = sum(accessibleforestcarbonstocksactual_PgC),
    carbonstocksallFAO_PgC = sum(accessibleforestcarbonstocksallFAO_PgC),
    carbonstocksallFAOupperrange_PgC = sum(accessibleforestcarbonstocksallFAOupperrange_PgC),
    carbonstocksallFAOlowerrange_PgC = sum(accessibleforestcarbonstocksallFAOlowerrange_PgC),
    # carbonstockspotential_PgC = sum(accessibleforestcarbonstockspotential_PgC)
  ) %>%
  
  summarise(globalcarbonstocksactual_PgC = sum(carbonstocksactual_PgC),
            globalcarbonstocksallFAO_PgC = sum(carbonstocksallFAO_PgC),
            globalcarbonstocksallFAOupperrange_PgC = sum(carbonstocksallFAOupperrange_PgC),
            globalcarbonstocksallFAOlowerrange_PgC = sum(carbonstocksallFAOlowerrange_PgC)) %>%
  mutate(region = 'Global') %>%
  mutate(sensitivityactual_PgC = globalcarbonstocksactual_PgC - globalcarbonstocksactual_PgC,
         sensitivityallFAO_PgC = globalcarbonstocksallFAO_PgC - globalcarbonstocksactual_PgC,
         sensitivityupperrange_PgC = globalcarbonstocksallFAOupperrange_PgC - globalcarbonstocksactual_PgC,
         sensitivitylowerrange_PgC = globalcarbonstocksallFAOlowerrange_PgC - globalcarbonstocksactual_PgC)



## ACCESSIBLE FOREST SENSITIVITY----

AccessibleForestSensitivity <- 
  data.frame(region = c('Global', 'Global', 'Global', 'Global'),
             focalcase = c('Best-guess', 'Best-guess', 'Best-guess', 'Best-guess'),
             accessibleforest = c(450.21 - AccessibleForestWithSC$sensitivityactual_PgC,
                                  450.21 - AccessibleForestWithSC$sensitivityallFAO_PgC,
                                  450.21 - AccessibleForestWithSC$sensitivityupperrange_PgC,
                                  450.21 - AccessibleForestWithSC$sensitivitylowerrange_PgC)) %>%
  mutate(sensitivitycategory = 'Accessible Forest') %>%
  mutate(sensitivity = 450.21 - accessibleforest) %>%
  dplyr::select(region, sensitivitycategory, sensitivity)


## Inaccessible forest with SC----

InaccessibleForestWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Select inaccessible forest
  dplyr::select(country, inaccessibleforest_Mha, subcontinent) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha), by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  left_join(Erb2018CarbonDensities, by = c('country')) %>%
  left_join(CarbonDensitiesWhereActGreaterThanPot, by = c('country')) %>%
  
  # Erb 2018
  mutate(naturalotherlandtcgreater10_tCha = as.numeric(naturalotherlandtcgreater10_tCha),
         wildforests_tCha = as.numeric(wildforests_tCha),
         wild_owl_tCha = as.numeric(wild_owl_tCha),
         naturalowltc5_10_tCha = as.numeric(naturalowltc5_10_tCha)) %>%
  
  # Inaccessible and extra forest carbon density
  mutate(nongrowingstockcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   potentialforestcarbon_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Inaccessible forest CS
  mutate(inaccessibleforestcarbonstocks_MtC = (inaccessibleforest_Mha * nongrowingstockcarbondensity_tCha)) %>%
  
  # Changing inaccessible forest CS for some countries
  mutate(inaccessibleforestcarbonstocks_MtC = if_else(country == 'Algeria' |
                                                        country == 'Australia' |
                                                        country == 'Austria' |
                                                        country == 'Canada' |
                                                        #country == 'China' |
                                                        country == 'Finland',
                                                      0.75 * (inaccessibleforest_Mha * forestcarbondensity_tCha), # can also assume it only reaches 75% as a measure of unproductive forests 
                                                      inaccessibleforestcarbonstocks_MtC)) %>%
  
  # total inaccessible forest CS
  mutate(inaccessibleforestcarbonstockspotential_PgC = inaccessibleforestcarbonstocks_MtC/1000,
         
         inaccessibleforestcarbonstocksactual_MtC = inaccessibleforest_Mha * forestcarbondensity_tCha,
         inaccessibleforestcarbonstocksactual_PgC = inaccessibleforestcarbonstocksactual_MtC/1000,
         
         inaccessibleforestcarbonstockserbwildforest_MtC = inaccessibleforest_Mha * wildforests_tCha,
         inaccessibleforestcarbonstockserbwildforest_PgC = inaccessibleforestcarbonstockserbwildforest_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstocksactual_PgC = sum(inaccessibleforestcarbonstocksactual_PgC),
            carbonstockspotential_PgC = sum(inaccessibleforestcarbonstockspotential_PgC),
            carbonstockserbwildforest_PgC = sum(inaccessibleforestcarbonstockserbwildforest_PgC)) %>%
  summarise(globalcarbonstocksactual_PgC = sum(carbonstocksactual_PgC),
            globalcarbonstockspotential_PgC = sum(carbonstockspotential_PgC),
            globalcarbonstockserbwildforest_PgC = sum(carbonstockserbwildforest_PgC)) %>%
  mutate(region = 'Global') %>%
  mutate(sensitivityactual_PgC = globalcarbonstocksactual_PgC - globalcarbonstockspotential_PgC,
         sensitivitypotential_PgC = globalcarbonstockspotential_PgC - globalcarbonstockspotential_PgC,
         sensitivityerbwildforest_PgC = globalcarbonstockserbwildforest_PgC - globalcarbonstockspotential_PgC)


## INACCESSIBLE FOREST SENSITIVITY----

InaccessibleForestSensitivity <- 
  data.frame(region = c('Global', 'Global', 'Global'),
             focalcase = c('Best-guess', 'Best-guess', 'Best-guess'),
             inaccessibleforest = c(450.21 - InaccessibleForestWithSC$sensitivityactual_PgC,
                                    450.21 - InaccessibleForestWithSC$sensitivitypotential_PgC,
                                    450.21 - InaccessibleForestWithSC$sensitivityerbwildforest_PgC)) %>%
  mutate(sensitivitycategory = 'Inaccessible Forest') %>%
  mutate(sensitivity = 450.21 - inaccessibleforest) %>%
  dplyr::select(region, sensitivitycategory, sensitivity) 

## Extra Forest with SC----

ExtraForestWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha), by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  left_join(Erb2018CarbonDensities, by = c('country')) %>%
  left_join(CarbonDensitiesWhereActGreaterThanPot, by = c('country')) %>%
  
  # Erb 2018
  mutate(naturalotherlandtcgreater10_tCha = as.numeric(naturalotherlandtcgreater10_tCha),
         wildforests_tCha = as.numeric(wildforests_tCha),
         wild_owl_tCha = as.numeric(wild_owl_tCha),
         naturalowltc5_10_tCha = as.numeric(naturalowltc5_10_tCha)) %>%
  
  # Extra forest carbon density
  mutate(nongrowingstockcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   potentialforestcarbon_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  # Extra Forest CS
  mutate(extraforestcarbonstockspotential_MtC = (extraforest_Mha * nongrowingstockcarbondensity_tCha),
         extraforestcarbonstockspotential_PgC = extraforestcarbonstockspotential_MtC/1000,
         
         extraforestcarbonstocksactual_MtC = extraforest_Mha * forestcarbondensity_tCha,
         extraforestcarbonstocksactual_PgC = extraforestcarbonstocksactual_MtC/1000,
         
         extraforestcarbonstockserbwildforest_MtC = extraforest_Mha * wildforests_tCha,
         extraforestcarbonstockserbwildforest_PgC = extraforestcarbonstockserbwildforest_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstocksactual_PgC = sum(extraforestcarbonstocksactual_PgC),
            carbonstockspotential_PgC = sum(extraforestcarbonstockspotential_PgC),
            carbonstockserbwildforest_PgC = sum(extraforestcarbonstockserbwildforest_PgC),) %>%
  summarise(globalcarbonstocksactual_PgC = sum(carbonstocksactual_PgC),
            globalcarbonstockspotential_PgC = sum(carbonstockspotential_PgC),
            globalcarbonstockserbwildforest_PgC = sum(carbonstockserbwildforest_PgC)) %>%
  mutate(region = 'Global') %>%
  mutate(sensitivityactual_PgC = globalcarbonstocksactual_PgC - globalcarbonstockspotential_PgC,
         sensitivitypotential_PgC = globalcarbonstockspotential_PgC - globalcarbonstockspotential_PgC,
         sensitivityerbwildforest_PgC = globalcarbonstockserbwildforest_PgC - globalcarbonstockspotential_PgC)


## EXTRA FOREST SENSITIVITY----

ExtraForestSensitivity <- 
  data.frame(region = c('Global', 'Global', 'Global'),
             focalcase = c('Best-guess', 'Best-guess', 'Best-guess'),
             extraforest = c(450.21 - ExtraForestWithSC$sensitivityactual_PgC,
                             450.21 - ExtraForestWithSC$sensitivitypotential_PgC,
                             450.21 - ExtraForestWithSC$sensitivityerbwildforest_PgC)) %>%
  mutate(sensitivitycategory = 'Extra Forest') %>%
  mutate(sensitivity = 450.21 - extraforest) %>%
  dplyr::select(region, sensitivitycategory, sensitivity) 


## OVL with SC----

OVLWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha), by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  left_join(PotentialCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  left_join(SubcontinentalNPP, by = c('subcontinent')) %>%
  left_join(Erb2018CarbonDensities, by = c('country')) %>%
  left_join(CarbonDensitiesWhereActGreaterThanPot, by = c('country')) %>%
  
  # OVL carbon density
  mutate(OVLcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   forestcarbondensity_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Inaccessible and extra forest carbon
  mutate(nongrowingstockcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   potentialforestcarbon_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Erb numeric
  mutate(naturalotherlandtcgreater10_tCha = as.numeric(naturalotherlandtcgreater10_tCha),
         wildforests_tCha = as.numeric(wildforests_tCha),
         wild_owl_tCha = as.numeric(wild_owl_tCha),
         naturalowltc5_10_tCha = as.numeric(naturalowltc5_10_tCha)) %>%
  
  # OVL CS
  mutate(OVLcarbonstockshalfgrowingstocks_MtC = othervegetatedland_Mha * (OVLcarbondensity_tCha/2),
         OVLcarbonstockshalfgrowingstocks_PgC = OVLcarbonstockshalfgrowingstocks_MtC/1000,
         
         OVLcarbonstocksfullgrowingstocks_MtC = othervegetatedland_Mha * (0.75 * OVLcarbondensity_tCha),
         OVLcarbonstocksfullgrowingstocks_PgC = OVLcarbonstocksfullgrowingstocks_MtC/1000,
         
         OVLcarbonstockshalfpotential_MtC = othervegetatedland_Mha * (nongrowingstockcarbondensity_tCha/2),
         OVLcarbonstockshalfpotential_PgC = OVLcarbonstockshalfpotential_MtC/1000,
         
         OVLcarbonstocksnaturalotherlandtcgreater10_MtC = othervegetatedland_Mha * (naturalotherlandtcgreater10_tCha),
         OVLcarbonstocksnaturalotherlandtcgreater10_PgC = OVLcarbonstocksnaturalotherlandtcgreater10_MtC/1000,
         
         #OVLcarbonstockswildowl_MtC = othervegetatedland_Mha * (wild_owl_tCha),
         #OVLcarbonstockswildowl_PgC = OVLcarbonstockswildowl_MtC/1000,
         
         OVLcarbonstocksnaturalowltc5_10_MtC = othervegetatedland_Mha * (naturalowltc5_10_tCha),
         OVLcarbonstocksnaturalowltc5_10_PgC = OVLcarbonstocksnaturalowltc5_10_MtC/1000,
         
         OVLcarbonstockswildforest_MtC = othervegetatedland_Mha * (0.5 * wildforests_tCha),
         OVLcarbonstockswildforest_PgC = OVLcarbonstockswildforest_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstockshalfgrowingstocks_PgC = sum(OVLcarbonstockshalfgrowingstocks_PgC),
            carbonstocksfullgrowingstocks_PgC = sum(OVLcarbonstocksfullgrowingstocks_PgC),
            carbonstockshalfpotential_PgC = sum(OVLcarbonstockshalfpotential_PgC),
            carbonstocksnaturalotherlandtcgreater10_PgC = sum(OVLcarbonstocksnaturalotherlandtcgreater10_PgC),
            #carbonstockswildowl_PgC = sum(OVLcarbonstockswildowl_PgC),
            carbonstocksnaturalowltc5_10_PgC = sum(OVLcarbonstocksnaturalowltc5_10_PgC),
            carbonstockswildforest_PgC = sum(OVLcarbonstockswildforest_PgC)) %>%
  
  summarise(globalcarbonstockshalfgrowingstocks_PgC = sum(carbonstockshalfgrowingstocks_PgC),
            globalcarbonstocksfullgrowingstocks_PgC = sum(carbonstocksfullgrowingstocks_PgC),
            globalcarbonstockshalfpotential_PgC = sum(carbonstockshalfpotential_PgC),
            globalcarbonstocksnaturalotherlandtcgreater10_PgC = sum(carbonstocksnaturalotherlandtcgreater10_PgC),
            globalcarbonstocksnaturalowltc5_10_PgC = sum(carbonstocksnaturalowltc5_10_PgC),
            globalcarbonstockswildforest_PgC = sum(carbonstockswildforest_PgC)) %>%
  mutate(region = 'Global') %>%
  mutate(sensitivityhalfgrowingstocks_PgC = globalcarbonstockshalfgrowingstocks_PgC - globalcarbonstockshalfgrowingstocks_PgC,
         sensitivityfullgrowingstocks_PgC = globalcarbonstocksfullgrowingstocks_PgC - globalcarbonstockshalfgrowingstocks_PgC,
         sensitivityhalfpotential_PgC =  globalcarbonstockshalfpotential_PgC - globalcarbonstockshalfgrowingstocks_PgC,
         sensitivitynaturalotherlandtcgreater10_PgC = globalcarbonstocksnaturalotherlandtcgreater10_PgC - globalcarbonstockshalfgrowingstocks_PgC,
         sensitivitynaturalowltc5_10_PgC = globalcarbonstocksnaturalowltc5_10_PgC - globalcarbonstockshalfgrowingstocks_PgC,
         sensitivitywildforest_PgC = globalcarbonstockswildforest_PgC - globalcarbonstockshalfgrowingstocks_PgC)


## OVL SENSITIVITY----

OVLSensitivity <- 
  data.frame(region = c('Global', 'Global', 'Global', 'Global', 'Global', 'Global'),
             focalcase = c('Best-guess', 'Best-guess', 'Best-guess', 'Best-guess', 'Best-guess', 'Best-guess'),
             OVL = c(450.21 - OVLWithSC$sensitivityhalfgrowingstocks_PgC,
                     450.21 - OVLWithSC$sensitivityfullgrowingstocks_PgC,
                     450.21 - OVLWithSC$sensitivityhalfpotential_PgC,
                     450.21 - OVLWithSC$sensitivitynaturalotherlandtcgreater10_PgC,
                     450.21 - OVLWithSC$sensitivitynaturalowltc5_10_PgC,
                     450.21 - OVLWithSC$sensitivitywildforest_PgC)) %>%
  mutate(sensitivitycategory = 'Other Vegetated Land') %>%
  mutate(sensitivity = 450.21 - OVL) %>%
  dplyr::select(region, sensitivitycategory, sensitivity) 


## Rangeland with SC----

RangelandWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha), by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  left_join(PotentialCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  left_join(SubcontinentalNPP, by = c('subcontinent')) %>%
  left_join(Erb2018CarbonDensities, by = c('country')) %>%
  left_join(Xia2014RangelandCS, by = c('subcontinent')) %>%
  left_join(NPPGrazing1950, by = c('subcontinent')) %>%
  left_join(NPPCropland1950, by = c('subcontinent')) %>%
  
  # Erb mutate
  mutate(naturalotherlandtcless5_tCha = as.numeric(naturalotherlandtcless5_tCha),
         wildtcless5_tCha = as.numeric(wildtcless5_tCha)) %>%
  
  # Rangeland CS
  mutate(
    #rangelandscarbonstocksNPP_MtC = rangelands_Mha * (NPPact_tChayr),
    #    rangelandscarbonstocksNPP_PgC = rangelandscarbonstocksNPP_MtC/1000,
    
    rangelandscarbonstocks2NPP_MtC = rangelands_Mha * (NPPactgrazing1950_tCha * 2),
    rangelandscarbonstocks2NPP_PgC = rangelandscarbonstocks2NPP_MtC/1000,
    
    rangelandscarbonstocksXia_MtC = rangelands_Mha * (rangelandCS_tCha),
    rangelandscarbonstocksXia_PgC = rangelandscarbonstocksXia_MtC/1000,
    
    rangelandscarbonstocksnaturalotherland_MtC = rangelands_Mha * (naturalotherlandtcless5_tCha),
    rangelandscarbonstocksnaturalotherland_PgC = rangelandscarbonstocksnaturalotherland_MtC/1000,
    
    rangelandscarbonstockswildotherland_MtC = rangelands_Mha * (wildtcless5_tCha),
    rangelandscarbonstockswildotherland_PgC = rangelandscarbonstockswildotherland_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstocks2NPP_PgC = sum(rangelandscarbonstocks2NPP_PgC),
            carbonstocksXia_PgC = sum(rangelandscarbonstocksXia_PgC),
            carbonstocksnaturalotherland_PgC = sum(rangelandscarbonstocksnaturalotherland_PgC),
            carbonstockswildotherland_PgC = sum(rangelandscarbonstockswildotherland_PgC)) %>%
  
  summarise(globalcarbonstocks2NPP_PgC = sum(carbonstocks2NPP_PgC),
            globalcarbonstocksXia_PgC = sum(carbonstocksXia_PgC),
            globalcarbonstocksnaturalotherland_PgC = sum(carbonstocksnaturalotherland_PgC),
            globalcarbonstockswildotherland_PgC = sum(carbonstockswildotherland_PgC)) %>%
  mutate(region = 'Global') %>%
  mutate(sensitivity2NPP_PgC = globalcarbonstocks2NPP_PgC - globalcarbonstocks2NPP_PgC,
         sensitivityXia_PgC = globalcarbonstocksXia_PgC - globalcarbonstocks2NPP_PgC,
         sensitivitynaturalotherland_PgC = globalcarbonstocksnaturalotherland_PgC - globalcarbonstocks2NPP_PgC,
         sensitivitywildotherland_PgC = globalcarbonstockswildotherland_PgC - globalcarbonstocks2NPP_PgC)


## RANGELAND SENSITIVITY----

RangelandSensitivity <- 
  data.frame(region = c('Global', 'Global', 'Global', 'Global'),
             focalcase = c('Best-guess', 'Best-guess', 'Best-guess', 'Best-guess'),
             rangeland = c(450.21 - RangelandWithSC$sensitivity2NPP_PgC,
                     450.21 - RangelandWithSC$sensitivityXia_PgC,
                     450.21 - RangelandWithSC$sensitivitynaturalotherland_PgC,
                     450.21 - RangelandWithSC$sensitivitywildotherland_PgC)) %>%
  mutate(sensitivitycategory = 'Rangeland') %>%
  mutate(sensitivity = 450.21 - rangeland) %>%
  dplyr::select(region, sensitivitycategory, sensitivity) 


## BIND ALL ROWS + PLOT----

(GlobalSensitivity <- rbind(AccessibleForestSensitivity, InaccessibleForestSensitivity,
                            ExtraForestSensitivity, OVLSensitivity, RangelandSensitivity) %>%
   filter(!(sensitivity = 0)) %>%
   mutate(sensitivitycategory = factor(sensitivitycategory,
                                       levels = c('Accessible Forest',
                                                  'Inaccessible Forest',
                                                  'Extra Forest',
                                                  'Other Vegetated Land',
                                                  'Rangeland'),
                                       labels = c('Ac\\nFor',
                                                  'Inac\\nFor',
                                                  'Ext\\nFor',
                                                  'OVL',
                                                  'Ran'))) %>%
   ggplot(aes(x = sensitivitycategory, y = sensitivity)) +
   geom_point(size = 2) +
   theme_bw() +
   scale_y_continuous(name = 'Carbon stocks (PgC)', limits = c(-50, 80), breaks = c(-50, 0, 50)) +
   theme(panel.grid = element_blank(),
         legend.title = element_blank(),
         axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
         axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
         axis.title.x = element_blank(),
         axis.title.y = element_text(size = 10, face = 'bold'),
         legend.text = element_text(size = 10, face = 'bold', color = 'black'), legend.position = 'none',
         legend.background = element_rect(fill = alpha(0.4))))

#ggsave(filename = 'SI Figure 5.jpg', GlobalSensitivity, dpi = 300, width = 5.5, height = 4, units = 'in')

## SUBCONTINENT----

## Accessible forest with SC----

AccessibleForestWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Select accessible forest
  dplyr::select(country, accessibleforest_Mha, subcontinent, climatezone) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha, forestcarbondensityallFAO_tCha, forestcarbondensitytype), 
            by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  
  # Forest CS
  mutate(accessibleforestcarbonstocksactual_MtC = (accessibleforest_Mha * forestcarbondensity_tCha),
         accessibleforestcarbonstocksactual_PgC = accessibleforestcarbonstocksactual_MtC/1000,
         
         # all FAO
         accessibleforestcarbonstocksallFAO_MtC = (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
         accessibleforestcarbonstocksallFAO_PgC = accessibleforestcarbonstocksallFAO_MtC/1000,
         
         # Upper range
         accessibleforestcarbonstocksallFAOupperrange_MtC = 
           case_when(climatezone == 'temperate' & forestcarbondensitytype == 'FAO-derived' ~ 
                       1.25 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'FAO-derived' ~ 
                       1.25 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'Continental average' ~ 
                       1.4 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'Independent studies' ~ 
                       1.25 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'temperate' & forestcarbondensitytype == 'Independent studies' ~ 
                       1.25 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'temperate' & forestcarbondensitytype == 'Continental average' ~ 
                       1.4 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'FAO-derived' ~ 
                       1.4 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'Independent studies' ~ 
                       1.4 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'Continental average' ~ 
                       1.5 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha)),
         
         accessibleforestcarbonstocksallFAOupperrange_PgC = accessibleforestcarbonstocksallFAOupperrange_MtC/1000,
         
         # Lower range
         accessibleforestcarbonstocksallFAOlowerrange_MtC = 
           case_when(climatezone == 'temperate' & forestcarbondensitytype == 'FAO-derived' ~ 
                       0.75 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'FAO-derived' ~ 
                       0.75 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'Continental average' ~ 
                       0.6 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'boreal' & forestcarbondensitytype == 'Independent studies' ~ 
                       0.75 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'temperate' & forestcarbondensitytype == 'Independent studies' ~ 
                       0.75 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'temperate' & forestcarbondensitytype == 'Continental average' ~ 
                       0.6 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'FAO-derived' ~ 
                       0.6 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'Independent studies' ~ 
                       0.6 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha),
                     
                     climatezone == 'tropical' & forestcarbondensitytype == 'Continental average' ~ 
                       0.5 * (accessibleforest_Mha * forestcarbondensityallFAO_tCha)),
         
         accessibleforestcarbonstocksallFAOlowerrange_PgC = accessibleforestcarbonstocksallFAOlowerrange_MtC/1000,
         
         #accessibleforestcarbonstockspotential_MtC = (accessibleforest_Mha * potentialforestcarbon_tCha),
         #accessibleforestcarbonstockspotential_PgC = accessibleforestcarbonstockspotential_MtC/1000
         
  ) %>%
  
  group_by(subcontinent) %>%
  summarise(
    carbonstocksactual_PgC = sum(accessibleforestcarbonstocksactual_PgC),
    carbonstocksallFAO_PgC = sum(accessibleforestcarbonstocksallFAO_PgC),
    carbonstocksallFAOupperrange_PgC = sum(accessibleforestcarbonstocksallFAOupperrange_PgC),
    carbonstocksallFAOlowerrange_PgC = sum(accessibleforestcarbonstocksallFAOlowerrange_PgC),
    # carbonstockspotential_PgC = sum(accessibleforestcarbonstockspotential_PgC)
  ) %>%
  
  mutate(landcategory = 'Accessible Forest') %>%
  mutate(sensitivityactual_PgC = carbonstocksactual_PgC - carbonstocksactual_PgC,
         sensitivityallFAO_PgC = carbonstocksallFAO_PgC - carbonstocksactual_PgC,
         sensitivityupperrange_PgC = carbonstocksallFAOupperrange_PgC - carbonstocksactual_PgC,
         sensitivitylowerrange_PgC = carbonstocksallFAOlowerrange_PgC - carbonstocksactual_PgC) %>%
  gather(sensitivitytype, estimate, sensitivityactual_PgC:sensitivitylowerrange_PgC) %>%
  dplyr::select(subcontinent, landcategory, sensitivitytype, estimate)


## Inaccessible forest with SC----

InaccessibleForestWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Select inaccessible forest
  dplyr::select(country, inaccessibleforest_Mha, subcontinent) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha), by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  left_join(Erb2018CarbonDensities, by = c('country')) %>%
  left_join(CarbonDensitiesWhereActGreaterThanPot, by = c('country')) %>%
  
  # Erb 2018
  mutate(naturalotherlandtcgreater10_tCha = as.numeric(naturalotherlandtcgreater10_tCha),
         wildforests_tCha = as.numeric(wildforests_tCha),
         wild_owl_tCha = as.numeric(wild_owl_tCha),
         naturalowltc5_10_tCha = as.numeric(naturalowltc5_10_tCha)) %>%
  
  # Inaccessible and extra forest carbon density
  mutate(nongrowingstockcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   potentialforestcarbon_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Inaccessible forest CS
  mutate(inaccessibleforestcarbonstocks_MtC = (inaccessibleforest_Mha * nongrowingstockcarbondensity_tCha)) %>%
  
  # Changing inaccessible forest CS for some countries
  mutate(inaccessibleforestcarbonstocks_MtC = if_else(country == 'Algeria' |
                                                        country == 'Australia' |
                                                        country == 'Austria' |
                                                        country == 'Canada' |
                                                        #country == 'China' |
                                                        country == 'Finland',
                                                      0.75 * (inaccessibleforest_Mha * forestcarbondensity_tCha), # can also assume it only reaches 75% as a measure of unproductive forests 
                                                      inaccessibleforestcarbonstocks_MtC)) %>%
  
  # total inaccessible forest CS
  mutate(inaccessibleforestcarbonstockspotential_PgC = inaccessibleforestcarbonstocks_MtC/1000,
         
         inaccessibleforestcarbonstocksactual_MtC = inaccessibleforest_Mha * forestcarbondensity_tCha,
         inaccessibleforestcarbonstocksactual_PgC = inaccessibleforestcarbonstocksactual_MtC/1000,
         
         inaccessibleforestcarbonstockserbwildforest_MtC = inaccessibleforest_Mha * wildforests_tCha,
         inaccessibleforestcarbonstockserbwildforest_PgC = inaccessibleforestcarbonstockserbwildforest_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstocksactual_PgC = sum(inaccessibleforestcarbonstocksactual_PgC),
            carbonstockspotential_PgC = sum(inaccessibleforestcarbonstockspotential_PgC),
            carbonstockserbwildforest_PgC = sum(inaccessibleforestcarbonstockserbwildforest_PgC)) %>%
  
  mutate(landcategory = 'Inaccessible Forest') %>%
  mutate(sensitivityactual_PgC = carbonstocksactual_PgC - carbonstockspotential_PgC,
         sensitivitypotential_PgC = carbonstockspotential_PgC - carbonstockspotential_PgC,
         sensitivityerbwildforest_PgC = carbonstockserbwildforest_PgC - carbonstockspotential_PgC) %>%
  gather(sensitivitytype, estimate, sensitivityactual_PgC:sensitivityerbwildforest_PgC) %>%
  dplyr::select(subcontinent, landcategory, sensitivitytype, estimate)




## Extra Forest with SC----

ExtraForestWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha), by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  left_join(Erb2018CarbonDensities, by = c('country')) %>%
  left_join(CarbonDensitiesWhereActGreaterThanPot, by = c('country')) %>%
  
  # Erb 2018
  mutate(naturalotherlandtcgreater10_tCha = as.numeric(naturalotherlandtcgreater10_tCha),
         wildforests_tCha = as.numeric(wildforests_tCha),
         wild_owl_tCha = as.numeric(wild_owl_tCha),
         naturalowltc5_10_tCha = as.numeric(naturalowltc5_10_tCha)) %>%
  
  # Extra forest carbon density
  mutate(nongrowingstockcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   potentialforestcarbon_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  # Extra Forest CS
  mutate(extraforestcarbonstockspotential_MtC = (extraforest_Mha * nongrowingstockcarbondensity_tCha),
         extraforestcarbonstockspotential_PgC = extraforestcarbonstockspotential_MtC/1000,
         
         extraforestcarbonstocksactual_MtC = extraforest_Mha * forestcarbondensity_tCha,
         extraforestcarbonstocksactual_PgC = extraforestcarbonstocksactual_MtC/1000,
         
         extraforestcarbonstockserbwildforest_MtC = extraforest_Mha * wildforests_tCha,
         extraforestcarbonstockserbwildforest_PgC = extraforestcarbonstockserbwildforest_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstocksactual_PgC = sum(extraforestcarbonstocksactual_PgC),
            carbonstockspotential_PgC = sum(extraforestcarbonstockspotential_PgC),
            carbonstockserbwildforest_PgC = sum(extraforestcarbonstockserbwildforest_PgC)) %>%
  
  mutate(landcategory = 'Extra Forest') %>%
  mutate(sensitivityactual_PgC = carbonstocksactual_PgC - carbonstockspotential_PgC,
         sensitivitypotential_PgC = carbonstockspotential_PgC - carbonstockspotential_PgC,
         sensitivityerbwildforest_PgC = carbonstockserbwildforest_PgC - carbonstockspotential_PgC) %>%
  gather(sensitivitytype, estimate, sensitivityactual_PgC:sensitivityerbwildforest_PgC) %>%
  dplyr::select(subcontinent, landcategory, sensitivitytype, estimate)



## OVL with SC----

OVLWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha), by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  left_join(PotentialCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  left_join(SubcontinentalNPP, by = c('subcontinent')) %>%
  left_join(Erb2018CarbonDensities, by = c('country')) %>%
  left_join(CarbonDensitiesWhereActGreaterThanPot, by = c('country')) %>%
  
  # OVL carbon density
  mutate(OVLcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   forestcarbondensity_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Inaccessible and extra forest carbon
  mutate(nongrowingstockcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   potentialforestcarbon_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Erb numeric
  mutate(naturalotherlandtcgreater10_tCha = as.numeric(naturalotherlandtcgreater10_tCha),
         wildforests_tCha = as.numeric(wildforests_tCha),
         wild_owl_tCha = as.numeric(wild_owl_tCha),
         naturalowltc5_10_tCha = as.numeric(naturalowltc5_10_tCha)) %>%
  
  # OVL CS
  mutate(OVLcarbonstockshalfgrowingstocks_MtC = othervegetatedland_Mha * (OVLcarbondensity_tCha/2),
         OVLcarbonstockshalfgrowingstocks_PgC = OVLcarbonstockshalfgrowingstocks_MtC/1000,
         
         OVLcarbonstocksfullgrowingstocks_MtC = othervegetatedland_Mha * (0.75 * OVLcarbondensity_tCha),
         OVLcarbonstocksfullgrowingstocks_PgC = OVLcarbonstocksfullgrowingstocks_MtC/1000,
         
         OVLcarbonstockshalfpotential_MtC = othervegetatedland_Mha * (nongrowingstockcarbondensity_tCha/2),
         OVLcarbonstockshalfpotential_PgC = OVLcarbonstockshalfpotential_MtC/1000,
         
         OVLcarbonstocksnaturalotherlandtcgreater10_MtC = othervegetatedland_Mha * (naturalotherlandtcgreater10_tCha),
         OVLcarbonstocksnaturalotherlandtcgreater10_PgC = OVLcarbonstocksnaturalotherlandtcgreater10_MtC/1000,
         
         #OVLcarbonstockswildowl_MtC = othervegetatedland_Mha * (wild_owl_tCha),
         #OVLcarbonstockswildowl_PgC = OVLcarbonstockswildowl_MtC/1000,
         
         OVLcarbonstocksnaturalowltc5_10_MtC = othervegetatedland_Mha * (naturalowltc5_10_tCha),
         OVLcarbonstocksnaturalowltc5_10_PgC = OVLcarbonstocksnaturalowltc5_10_MtC/1000,
         
         OVLcarbonstockswildforest_MtC = othervegetatedland_Mha * (0.5 * wildforests_tCha),
         OVLcarbonstockswildforest_PgC = OVLcarbonstockswildforest_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstockshalfgrowingstocks_PgC = sum(OVLcarbonstockshalfgrowingstocks_PgC),
            carbonstocksfullgrowingstocks_PgC = sum(OVLcarbonstocksfullgrowingstocks_PgC),
            carbonstockshalfpotential_PgC = sum(OVLcarbonstockshalfpotential_PgC),
            carbonstocksnaturalotherlandtcgreater10_PgC = sum(OVLcarbonstocksnaturalotherlandtcgreater10_PgC),
            #carbonstockswildowl_PgC = sum(OVLcarbonstockswildowl_PgC),
            carbonstocksnaturalowltc5_10_PgC = sum(OVLcarbonstocksnaturalowltc5_10_PgC),
            carbonstockswildforest_PgC = sum(OVLcarbonstockswildforest_PgC)) %>%
  
  mutate(landcategory = 'Other Vegetated Land') %>%
  
  mutate(sensitivityhalfgrowingstocks_PgC = carbonstockshalfgrowingstocks_PgC - carbonstockshalfgrowingstocks_PgC,
         sensitivityfullgrowingstocks_PgC = carbonstocksfullgrowingstocks_PgC - carbonstockshalfgrowingstocks_PgC,
         sensitivityhalfpotential_PgC =  carbonstockshalfpotential_PgC - carbonstockshalfgrowingstocks_PgC,
         sensitivitynaturalotherlandtcgreater10_PgC = carbonstocksnaturalotherlandtcgreater10_PgC - carbonstockshalfgrowingstocks_PgC,
         sensitivitynaturalowltc5_10_PgC = carbonstocksnaturalowltc5_10_PgC - carbonstockshalfgrowingstocks_PgC,
         sensitivitywildforest_PgC = carbonstockswildforest_PgC - carbonstockshalfgrowingstocks_PgC) %>%
  gather(sensitivitytype, estimate, sensitivityhalfgrowingstocks_PgC:sensitivitywildforest_PgC) %>%
  dplyr::select(subcontinent, landcategory, sensitivitytype, estimate)


## Rangeland with SC----

RangelandWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
  # Add continental classification
  left_join(SubContinentalClassification, by = c('country')) %>%
  
  # Add all CS densities
  left_join(FRAForestCarbonDensities %>%
              dplyr::select(country, forestcarbondensity_tCha), by = c('country')) %>%
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  rename(potentialforestcarbon_tCha = potential_tCha) %>%
  left_join(PotentialCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  left_join(SubcontinentalNPP, by = c('subcontinent')) %>%
  left_join(Erb2018CarbonDensities, by = c('country')) %>%
  left_join(Xia2014RangelandCS, by = c('subcontinent')) %>%
  left_join(NPPGrazing1950, by = c('subcontinent')) %>%
  left_join(NPPCropland1950, by = c('subcontinent')) %>%
  
  # Erb mutate
  mutate(naturalotherlandtcless5_tCha = as.numeric(naturalotherlandtcless5_tCha),
         wildtcless5_tCha = as.numeric(wildtcless5_tCha)) %>%
  
  # Rangeland CS
  mutate(
    #rangelandscarbonstocksNPP_MtC = rangelands_Mha * (NPPact_tChayr),
    #    rangelandscarbonstocksNPP_PgC = rangelandscarbonstocksNPP_MtC/1000,
    
    rangelandscarbonstocks2NPP_MtC = rangelands_Mha * (NPPactgrazing1950_tCha * 2),
    rangelandscarbonstocks2NPP_PgC = rangelandscarbonstocks2NPP_MtC/1000,
    
    rangelandscarbonstocksXia_MtC = rangelands_Mha * (rangelandCS_tCha),
    rangelandscarbonstocksXia_PgC = rangelandscarbonstocksXia_MtC/1000,
    
    rangelandscarbonstocksnaturalotherland_MtC = rangelands_Mha * (naturalotherlandtcless5_tCha),
    rangelandscarbonstocksnaturalotherland_PgC = rangelandscarbonstocksnaturalotherland_MtC/1000,
    
    rangelandscarbonstockswildotherland_MtC = rangelands_Mha * (wildtcless5_tCha),
    rangelandscarbonstockswildotherland_PgC = rangelandscarbonstockswildotherland_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstocks2NPP_PgC = sum(rangelandscarbonstocks2NPP_PgC),
            carbonstocksXia_PgC = sum(rangelandscarbonstocksXia_PgC),
            carbonstocksnaturalotherland_PgC = sum(rangelandscarbonstocksnaturalotherland_PgC),
            carbonstockswildotherland_PgC = sum(rangelandscarbonstockswildotherland_PgC)) %>%
  
  mutate(landcategory = 'Rangeland') %>%
  mutate(sensitivity2NPP_PgC = carbonstocks2NPP_PgC - carbonstocks2NPP_PgC,
         sensitivityXia_PgC = carbonstocksXia_PgC - carbonstocks2NPP_PgC,
         sensitivitynaturalotherland_PgC = carbonstocksnaturalotherland_PgC - carbonstocks2NPP_PgC,
         sensitivitywildotherland_PgC = carbonstockswildotherland_PgC - carbonstocks2NPP_PgC) %>%
  gather(sensitivitytype, estimate, sensitivity2NPP_PgC:sensitivitywildotherland_PgC) %>%
  dplyr::select(subcontinent, landcategory, sensitivitytype, estimate)

## BIND ALL ROWS + PLOT----

(SubcontinentSensitivity <- rbind(AccessibleForestWithSC,
                                  InaccessibleForestWithSC,
                                  ExtraForestWithSC,
                                  OVLWithSC,
                                  RangelandWithSC) %>%
   filter(!(landcategory == 'Accessible Forest' & estimate == 0 & sensitivitytype == 'sensitivityactual_PgC')) %>%
   filter(!(landcategory == 'Inaccessible Forest' & estimate == 0 & sensitivitytype == 'sensitivitypotential_PgC')) %>%
   filter(!(landcategory == 'Extra Forest' & estimate == 0 & sensitivitytype == 'sensitivitypotential_PgC')) %>%
   filter(!(landcategory == 'Other Vegetated Land' & estimate == 0 & sensitivitytype == 'sensitivityhalfgrowingstocks_PgC')) %>%
   filter(!(landcategory == 'Rangeland' & estimate == 0 & sensitivitytype == 'sensitivity2NPP_PgC')) %>%
   mutate(landcategory = factor(landcategory,
                                levels = c('Accessible Forest',
                                           'Inaccessible Forest',
                                           'Extra Forest',
                                           'Other Vegetated Land',
                                           'Rangeland'),
                                labels = c('Ac\\nFor',
                                           'Inac\\nFor',
                                           'Ext\\nFor',
                                           'OVL',
                                           'Ran'))) %>%
   ggplot(aes(x = landcategory, y = estimate)) +
   geom_point(size = 2) +
   theme_bw() +
   scale_y_continuous(name = 'Carbon stocks (PgC)') +
   facet_wrap(~ subcontinent, ncol = 5, scales = "free_y", labeller = label_wrap_gen(width = 20)) +
   theme(panel.grid = element_blank(),
         legend.title = element_blank(),
         axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
         axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
         axis.title.x = element_blank(),
         axis.title.y = element_text(size = 10, face = 'bold'),
         legend.text = element_text(size = 10, face = 'bold', color = 'black'), legend.position = 'none',
         legend.background = element_rect(fill = alpha(0.4))))

#ggsave(filename = 'SI Figure 6.jpg', SubcontinentSensitivity, dpi = 300, width = 11, height = 8.5, units = 'in')
