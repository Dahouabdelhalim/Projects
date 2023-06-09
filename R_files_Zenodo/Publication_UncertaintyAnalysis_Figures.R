## A mid-20th century estimate of global vegetation carbon stocks based on inventory statistics
## Authors: Bhan et al., Journal of Land Use Science
## Uncertainty Analysis and Figures

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

## LAND ACCOUNT WITHOUT SC----

## Accessible forest----

AccessibleForest <- LandUseEstimates %>%
  
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
         
         # accessibleforestcarbonstockspotential_MtC = (accessibleforest_Mha * potentialforestcarbon_tCha),
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
  gather(actpot, accessibleforestcarbonstocks_PgC, carbonstocksactual_PgC:carbonstocksallFAOlowerrange_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(accessibleforestcarbonstocks_PgC))
  dplyr::select(subcontinent, accessibleforestcarbonstocks_PgC)


#AccessibleForest <- AccessibleForest %>% mutate(landmanagement_PgC = carbonstockspotential_PgC - carbonstocksactual_PgC, landmanagementperc_PgC = landmanagement_PgC/carbonstockspotential_PgC)

## Inaccessible forest----

InaccessibleForest <- LandUseEstimates %>%
  
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
            carbonstockserbwildforest_PgC = sum(inaccessibleforestcarbonstockserbwildforest_PgC),) %>%
  gather(actpot, inaccessibleforestcarbonstocks_PgC, carbonstocksactual_PgC:carbonstockserbwildforest_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(inaccessibleforestcarbonstocks_PgC))
  dplyr::select(subcontinent, inaccessibleforestcarbonstocks_PgC)

## ACCESSIBLE + INACCESSIBLE FOREST----

AccessibleInaccessibleCarbonStocks <- 
  full_join(AccessibleForest, 
            InaccessibleForest, by = c('subcontinent'))

## Extra forest----

ExtraForest <- LandUseEstimates %>%
  
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
  gather(actpot, extraforestcarbonstocks_PgC, carbonstocksactual_PgC:carbonstockserbwildforest_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(extraforestcarbonstocks_PgC))
  dplyr::select(subcontinent, extraforestcarbonstocks_PgC)

## TOTAL FOREST----

TotalForestCarbonStocks <- 
  full_join(AccessibleInaccessibleCarbonStocks, 
            ExtraForest, by = c('subcontinent'))

# Other vegetated land----

OVL <- LandUseEstimates %>%
  
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
         
         OVLcarbonstocks75growingstocks_MtC = othervegetatedland_Mha * (0.75 * OVLcarbondensity_tCha),
         OVLcarbonstocks75growingstocks_PgC = OVLcarbonstocks75growingstocks_MtC/1000,
         
         OVLcarbonstockshalfpotential_MtC = othervegetatedland_Mha * (nongrowingstockcarbondensity_tCha/2),
         OVLcarbonstockshalfpotential_PgC = OVLcarbonstockshalfpotential_MtC/1000,
         
         OVLcarbonstocksnaturalotherlandtcgreater10_MtC = othervegetatedland_Mha * (naturalotherlandtcgreater10_tCha),
         OVLcarbonstocksnaturalotherlandtcgreater10_PgC = OVLcarbonstocksnaturalotherlandtcgreater10_MtC/1000,
         
         #OVLcarbonstockswildowl_MtC = othervegetatedland_Mha * (wild_owl_tCha),
         #OVLcarbonstockswildowl_PgC = OVLcarbonstockswildowl_MtC/1000,
         
         OVLcarbonstocksnaturalowltc5_10_MtC = othervegetatedland_Mha * (naturalowltc5_10_tCha),
         OVLcarbonstocksnaturalowltc5_10_PgC = OVLcarbonstocksnaturalowltc5_10_MtC/1000,
         
         OVLcarbonstocks50wildforest_MtC = othervegetatedland_Mha * (0.5 * wildforests_tCha),
         OVLcarbonstocks50wildforest_PgC = OVLcarbonstocks50wildforest_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstockshalfgrowingstocks_PgC = sum(OVLcarbonstockshalfgrowingstocks_PgC),
            carbonstocks75growingstocks_PgC = sum(OVLcarbonstocks75growingstocks_PgC),
            carbonstockshalfpotential_PgC = sum(OVLcarbonstockshalfpotential_PgC),
            carbonstocksnaturalotherlandtcgreater10_PgC = sum(OVLcarbonstocksnaturalotherlandtcgreater10_PgC),
            #carbonstockswildowl_PgC = sum(OVLcarbonstockswildowl_PgC),
            carbonstocksnaturalowltc5_10_PgC = sum(OVLcarbonstocksnaturalowltc5_10_PgC),
            carbonstocks50wildforest_PgC = sum(OVLcarbonstocks50wildforest_PgC)) %>%
  
  gather(actpot, OVLcarbonstocks_PgC, carbonstockshalfgrowingstocks_PgC:carbonstocks50wildforest_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(OVLcarbonstocks_PgC))
  dplyr::select(subcontinent, OVLcarbonstocks_PgC)

## ALL TREE-BEARING----

AllTreeBearingCarbonStocks <- 
  full_join(TotalForestCarbonStocks,
            OVL, by = c('subcontinent'))

## Rangelands----

Rangeland <- LandUseEstimates %>%
  
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
    #rangelandscarbonstocksNPP_PgC = rangelandscarbonstocksNPP_MtC/1000,
    
    rangelandscarbonstocks2NPP_MtC = rangelands_Mha * (NPPactgrazing1950_tCha * 2),
    rangelandscarbonstocks2NPP_PgC = rangelandscarbonstocks2NPP_MtC/1000,
    
    rangelandscarbonstocksXia_MtC = rangelands_Mha * (rangelandCS_tCha),
    rangelandscarbonstocksXia_PgC = rangelandscarbonstocksXia_MtC/1000,
    
    rangelandscarbonstocksnaturalotherland_MtC = rangelands_Mha * (naturalotherlandtcless5_tCha),
    rangelandscarbonstocksnaturalotherland_PgC = rangelandscarbonstocksnaturalotherland_MtC/1000,
    
    rangelandscarbonstockswildotherland_MtC = rangelands_Mha * (wildtcless5_tCha),
    rangelandscarbonstockswildotherland_PgC = rangelandscarbonstockswildotherland_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(#carbonstocksNPP_PgC = sum(rangelandscarbonstocksNPP_PgC),
    carbonstocks2NPP_PgC = sum(rangelandscarbonstocks2NPP_PgC),
    carbonstocksXia_PgC = sum(rangelandscarbonstocksXia_PgC),
    carbonstocksnaturalotherland_PgC = sum(rangelandscarbonstocksnaturalotherland_PgC),
    carbonstockswildotherland_PgC = sum(rangelandscarbonstockswildotherland_PgC)) %>%
  
  gather(actpot, rangelandcarbonstocks_PgC, carbonstocks2NPP_PgC:carbonstockswildotherland_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(rangelandcarbonstocks_PgC))
  dplyr::select(subcontinent, rangelandcarbonstocks_PgC)

## TREE-BEARING + RANGELAND----

TreeBearingRangeland <- 
  full_join(AllTreeBearingCarbonStocks,
            Rangeland, by = c('subcontinent'))

## Infra + Crops + Pasture----

CropPastureInfrastructure <- LandUseEstimates %>%
  
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
  left_join(NPPGrazing1950, by = c('subcontinent')) %>%
  left_join(NPPCropland1950, by = c('subcontinent')) %>%
  
  # Cropland CS
  mutate(annualcropscarbonstocks_MtC = (AnnualCrops_Mha * totalNPPact1950_tCha),
         permanentcropscarbonstocks_MtC = (0.5 * PermanentCrops_Mha) * 15 + (0.5 * PermanentCrops_Mha) * 30) %>%
  
  # pasture CS
  mutate(pasturescarbonstocks_MtC = pastures_Mha * NPPactgrazing1950_tCha) %>%
  
  # infrastructure CS
  mutate(infrastructurecarbonstocks_MtC = infrastructure_Mha * (potential_tCha/6)) %>%
  
  # total ecosystem CS
  mutate(CropPastureInfrastructure_MtC = 
           annualcropscarbonstocks_MtC + 
           permanentcropscarbonstocks_MtC +
           pasturescarbonstocks_MtC +
           infrastructurecarbonstocks_MtC,
         CropPastureInfrastructure_PgC = CropPastureInfrastructure_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(croppastureinfrastructure_PgC = sum(CropPastureInfrastructure_PgC))

#sum(CropPastureInfrastructure$croppastureinfrastructure_PgC) # 9.91

## Ecosystem----

EcosystemCarbonStocksWithoutSC <- 
  full_join(TreeBearingRangeland,
            CropPastureInfrastructure,
            by = c('subcontinent')) %>%
  mutate(ecosystemcarbonstocks_PgC = accessibleforestcarbonstocks_PgC +
           inaccessibleforestcarbonstocks_PgC +
           extraforestcarbonstocks_PgC +
           OVLcarbonstocks_PgC +
           rangelandcarbonstocks_PgC +
           croppastureinfrastructure_PgC)


## LAND ACCOUNT WITH SC----

# Accessible forest----

AccessibleForest <- LandUseEstimatesWithShiftingCultivation %>%
  
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
    carbonstocksallFAOlowerrange_PgC = sum(accessibleforestcarbonstocksallFAOlowerrange_PgC)
    #        carbonstockspotential_PgC = sum(accessibleforestcarbonstockspotential_PgC)
  ) %>%
  gather(actpot, accessibleforestcarbonstocks_PgC, carbonstocksactual_PgC:carbonstocksallFAOlowerrange_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(accessibleforestcarbonstocks_PgC))
  dplyr::select(subcontinent, accessibleforestcarbonstocks_PgC)


## Inaccessible forest----

InaccessibleForest <- LandUseEstimatesWithShiftingCultivation %>%
  
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
            carbonstockserbwildforest_PgC = sum(inaccessibleforestcarbonstockserbwildforest_PgC),) %>%
  gather(actpot, inaccessibleforestcarbonstocks_PgC, carbonstocksactual_PgC:carbonstockserbwildforest_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(inaccessibleforestcarbonstocks_PgC))
  dplyr::select(subcontinent, inaccessibleforestcarbonstocks_PgC)

## ACCESSIBLE + INACCESSIBLE FOREST----

AccessibleInaccessibleCarbonStocks <- 
  full_join(AccessibleForest, 
            InaccessibleForest, by = c("subcontinent"))

## Extra forest----

ExtraForest <- LandUseEstimatesWithShiftingCultivation %>%
  
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
  gather(actpot, extraforestcarbonstocks_PgC, carbonstocksactual_PgC:carbonstockserbwildforest_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(extraforestcarbonstocks_PgC))
  dplyr::select(subcontinent, extraforestcarbonstocks_PgC)

## TOTAL FOREST----

TotalForestCarbonStocks <- 
  full_join(AccessibleInaccessibleCarbonStocks, 
            ExtraForest, by = c('subcontinent'))

## SC land----

ShiftingCultivation <- LandUseEstimatesWithShiftingCultivation %>%
  
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
  left_join(SCCarbonDensities, by = c('country')) %>%
  
  # Shifting cultivation
  mutate(totalC_shortfallow_tCha = if_else(is.na(totalC_shortfallow_tCha), 0, totalC_shortfallow_tCha),
         totalC_longfallow_tCha = if_else(is.na(totalC_longfallow_tCha), 0, totalC_longfallow_tCha)) %>%
  
  mutate(SCshortfallowcarbonstocks_MtC = (shortfallow_Mha * totalC_shortfallow_tCha),
         SCshortfallowcarbonstocks_PgC = SCshortfallowcarbonstocks_MtC/1000,
         SClongfallowcarbonstocks_MtC = (longfallow_Mha * totalC_longfallow_tCha),
         SClongfallowcarbonstocks_PgC = SClongfallowcarbonstocks_MtC/1000,
         SCcarbonstocks_PgC = SCshortfallowcarbonstocks_PgC + SClongfallowcarbonstocks_PgC) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstocksSC_PgC = sum(SCcarbonstocks_PgC))

sum(ShiftingCultivation$carbonstocksSC_PgC) # 10.22

## SC + FOREST----

TotalForestSCCarbonStocks <- 
  full_join(TotalForestCarbonStocks, 
            ShiftingCultivation, by = c('subcontinent'))

# Other vegetated land----

OVL <- LandUseEstimatesWithShiftingCultivation %>%
  
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
         
         OVLcarbonstocks75growingstocks_MtC = othervegetatedland_Mha * (0.75 * OVLcarbondensity_tCha),
         OVLcarbonstocks75growingstocks_PgC = OVLcarbonstocks75growingstocks_MtC/1000,
         
         OVLcarbonstockshalfpotential_MtC = othervegetatedland_Mha * (nongrowingstockcarbondensity_tCha/2),
         OVLcarbonstockshalfpotential_PgC = OVLcarbonstockshalfpotential_MtC/1000,
         
         OVLcarbonstocksnaturalotherlandtcgreater10_MtC = othervegetatedland_Mha * (naturalotherlandtcgreater10_tCha),
         OVLcarbonstocksnaturalotherlandtcgreater10_PgC = OVLcarbonstocksnaturalotherlandtcgreater10_MtC/1000,
         
         #OVLcarbonstockswildowl_MtC = othervegetatedland_Mha * (wild_owl_tCha),
         #OVLcarbonstockswildowl_PgC = OVLcarbonstockswildowl_MtC/1000,
         
         OVLcarbonstocksnaturalowltc5_10_MtC = othervegetatedland_Mha * (naturalowltc5_10_tCha),
         OVLcarbonstocksnaturalowltc5_10_PgC = OVLcarbonstocksnaturalowltc5_10_MtC/1000,
         
         OVLcarbonstocks50wildforest_MtC = othervegetatedland_Mha * (0.5 * wildforests_tCha),
         OVLcarbonstocks50wildforest_PgC = OVLcarbonstocks50wildforest_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(carbonstockshalfgrowingstocks_PgC = sum(OVLcarbonstockshalfgrowingstocks_PgC),
            carbonstocks75growingstocks_PgC = sum(OVLcarbonstocks75growingstocks_PgC),
            carbonstockshalfpotential_PgC = sum(OVLcarbonstockshalfpotential_PgC),
            carbonstocksnaturalotherlandtcgreater10_PgC = sum(OVLcarbonstocksnaturalotherlandtcgreater10_PgC),
            #carbonstockswildowl_PgC = sum(OVLcarbonstockswildowl_PgC),
            carbonstocksnaturalowltc5_10_PgC = sum(OVLcarbonstocksnaturalowltc5_10_PgC),
            carbonstocks50wildforest_PgC = sum(OVLcarbonstocks50wildforest_PgC)) %>%
  
  gather(actpot, OVLcarbonstocks_PgC, carbonstockshalfgrowingstocks_PgC:carbonstocks50wildforest_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(OVLcarbonstocks_PgC))
  dplyr::select(subcontinent, OVLcarbonstocks_PgC)


## ALL TREE-BEARING----

AllTreeBearingCarbonStocks <- 
  full_join(TotalForestSCCarbonStocks,
            OVL, by = c('subcontinent'))

## Rangelands----

Rangeland <- LandUseEstimatesWithShiftingCultivation %>%
  
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
  summarise(
    #carbonstocksNPP_PgC = sum(rangelandscarbonstocksNPP_PgC),
    carbonstocks2NPP_PgC = sum(rangelandscarbonstocks2NPP_PgC),
    carbonstocksXia_PgC = sum(rangelandscarbonstocksXia_PgC),
    carbonstocksnaturalotherland_PgC = sum(rangelandscarbonstocksnaturalotherland_PgC),
    carbonstockswildotherland_PgC = sum(rangelandscarbonstockswildotherland_PgC)) %>%
  
  gather(actpot, rangelandcarbonstocks_PgC, carbonstocks2NPP_PgC:carbonstockswildotherland_PgC) %>%
  #group_by(actpot) %>%
  #summarise(carbonstocks_PgC = sum(rangelandcarbonstocks_PgC))
  dplyr::select(subcontinent, rangelandcarbonstocks_PgC)

## TREE-BEARING + RANGELAND----

TreeBearingRangeland <- 
  full_join(AllTreeBearingCarbonStocks,
            Rangeland, by = c('subcontinent'))

## Infra + Crops + Pasture----

CropPastureInfrastructure <- LandUseEstimatesWithShiftingCultivation %>%
  
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
  left_join(NPPGrazing1950, by = c('subcontinent')) %>%
  left_join(NPPCropland1950, by = c('subcontinent')) %>%
  
  # Cropland CS
  mutate(annualcropscarbonstocks_MtC = (AnnualCrops_Mha * totalNPPact1950_tCha),
         permanentcropscarbonstocks_MtC = (0.5 * PermanentCrops_Mha) * 15 + (0.5 * PermanentCrops_Mha) * 30) %>%
  
  # pasture CS
  mutate(pasturescarbonstocks_MtC = pastures_Mha * NPPactgrazing1950_tCha) %>%
  
  # infrastructure CS
  mutate(infrastructurecarbonstocks_MtC = infrastructure_Mha * (potential_tCha/6)) %>%
  
  # total ecosystem CS
  mutate(CropPastureInfrastructure_MtC = 
           annualcropscarbonstocks_MtC + 
           permanentcropscarbonstocks_MtC +
           pasturescarbonstocks_MtC +
           infrastructurecarbonstocks_MtC,
         CropPastureInfrastructure_PgC = CropPastureInfrastructure_MtC/1000) %>%
  
  group_by(subcontinent) %>%
  summarise(croppastureinfrastructure_PgC = sum(CropPastureInfrastructure_PgC))

## Ecosystem----

EcosystemCarbonStocksWithSC <- 
  full_join(TreeBearingRangeland,
            CropPastureInfrastructure,
            by = c('subcontinent')) %>%
  mutate(ecosystemcarbonstocks_PgC = accessibleforestcarbonstocks_PgC +
           inaccessibleforestcarbonstocks_PgC +
           extraforestcarbonstocks_PgC +
           OVLcarbonstocks_PgC +
           carbonstocksSC_PgC +
           rangelandcarbonstocks_PgC +
           croppastureinfrastructure_PgC)

## WITH + WITHOUT SC----

EcosystemCarbonStocksWithoutSC <- EcosystemCarbonStocksWithoutSC %>%
  mutate(carbonstocksSC_PgC = 0)

EcosystemCarbonStocks <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  group_by(subcontinent)

EcosystemCarbonStocksForest <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  group_by(subcontinent) %>%
  mutate(forestcarbonstocks_PgC = accessibleforestcarbonstocks_PgC + inaccessibleforestcarbonstocks_PgC + extraforestcarbonstocks_PgC,
         nonforestcarbonstocks_PgC = OVLcarbonstocks_PgC + rangelandcarbonstocks_PgC + croppastureinfrastructure_PgC + carbonstocksSC_PgC) %>%
  summarize(medianforestcarbonstocks_PgC = median(forestcarbonstocks_PgC, na.rm = TRUE),
            lowerforestcarbonstocks_PgC = quantile(forestcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            upperforestcarbonstocks_PgC = quantile(forestcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            mediannonforestcarbonstocks_PgC = median(nonforestcarbonstocks_PgC, na.rm = TRUE),
            lowernonforestcarbonstocks_PgC = quantile(nonforestcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            uppernonforestcarbonstocks_PgC = quantile(nonforestcarbonstocks_PgC, p = 0.75, na.rm = TRUE))

sum(EcosystemCarbonStocksForest$medianforestcarbonstocks_PgC) # 376.46
sum(EcosystemCarbonStocksForest$lowerforestcarbonstocks_PgC) # 318.23
sum(EcosystemCarbonStocksForest$upperforestcarbonstocks_PgC) # 441.85

sum(EcosystemCarbonStocksForest$mediannonforestcarbonstocks_PgC) # 131.03
sum(EcosystemCarbonStocksForest$lowernonforestcarbonstocks_PgC) # 109.21
sum(EcosystemCarbonStocksForest$uppernonforestcarbonstocks_PgC) # 166.26

EcosystemCarbonStocksMedian <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  group_by(subcontinent) %>%
  summarise(subcontinentmedian = median(ecosystemcarbonstocks_PgC, na.rm = TRUE)) %>%
  mutate(subcontinentabbreviation = 
           case_when(subcontinent == 'Central America & the Caribbean' ~ 'CAmC',
                     subcontinent == 'Eastern Africa' ~ 'EAf',
                     subcontinent == 'Eastern Asia' ~ 'EAs',
                     subcontinent == 'Eastern Europe' ~ 'EEu',
                     subcontinent == 'Northern Africa & Western Asia' ~ 'NAWA',
                     subcontinent == 'Northern America' ~ 'NAm',
                     subcontinent == 'Oceania' ~ 'Oc',
                     subcontinent == 'Southeastern Asia' ~ 'SEAs',
                     subcontinent == 'Southern Africa' ~ 'SAf',
                     subcontinent == 'Southern America' ~ 'SAm',
                     subcontinent == 'Southern Asia' ~ 'SAs',
                     subcontinent == 'Soviet Union' ~ 'SU',
                     subcontinent == 'Western Africa' ~ 'WAf',
                     subcontinent == 'Western Europe' ~ 'WEu'))

## GLOBAL SUMMARY STAT----

EcosystemCarbonStocksSummaryStat <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  group_by(subcontinent) %>%
  summarise(subcontinentmedian = median(ecosystemcarbonstocks_PgC, na.rm = TRUE),
            subcontinentmax = max(ecosystemcarbonstocks_PgC, na.rm = TRUE),
            subcontinentmin = min(ecosystemcarbonstocks_PgC, na.rm = TRUE),
            subcontinentlowerquantile = quantile(ecosystemcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentupperquantile = quantile(ecosystemcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            subcontinentmean = mean(ecosystemcarbonstocks_PgC, na.rm = TRUE))


sum(EcosystemCarbonStocksSummaryStat$subcontinentmedian) # 517.84 PgC
sum(EcosystemCarbonStocksSummaryStat$subcontinentmax) # 762.03 PgC
sum(EcosystemCarbonStocksSummaryStat$subcontinentmin) #  327.92 PgC
sum(EcosystemCarbonStocksSummaryStat$subcontinentupperquantile) # 584.02 PgC
sum(EcosystemCarbonStocksSummaryStat$subcontinentlowerquantile) # 443.69 PgC
sum(EcosystemCarbonStocksSummaryStat$subcontinentmean) # 518.05 PgC

#write.csv(EcosystemCarbonStocksSummaryStat, 'PublicationCarbonStocksSubcontinentSummary.csv')

## SUMMARY STAT BY LAND CATEGORY----

EcosystemCarbonStocksLandCategorySummaryStat <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  mutate(totalforestcarbonstocks_PgC = 
           accessibleforestcarbonstocks_PgC +
           inaccessibleforestcarbonstocks_PgC + 
           extraforestcarbonstocks_PgC,
         
         totalnonforestcarbonstocks_PgC = 
           OVLcarbonstocks_PgC + 
           rangelandcarbonstocks_PgC + 
           croppastureinfrastructure_PgC +
           carbonstocksSC_PgC) %>%
  
  group_by(subcontinent) %>%
  summarise(subcontinenttotalforestmedian = median(totalforestcarbonstocks_PgC, na.rm = TRUE),
            subcontinentforestlowerquantile = quantile(totalforestcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentforestupperquantile = quantile(totalforestcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            
            subcontinentextraforestmedian = median(extraforestcarbonstocks_PgC, na.rm = TRUE),
            subcontinentextraforestlowerquantile = quantile(extraforestcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentextraforestupperquantile = quantile(extraforestcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            
            subcontinenttotalnonforestmedian = median(totalnonforestcarbonstocks_PgC, na.rm = TRUE),
            subcontinentnonforestlowerquantile = quantile(totalnonforestcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentnonforestupperquantile = quantile(totalnonforestcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            
            subcontinentOVLmedian = median(OVLcarbonstocks_PgC, na.rm = TRUE),
            subcontinentOVLlowerquantile = quantile(OVLcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentOVLupperquantile = quantile(OVLcarbonstocks_PgC, p = 0.75, na.rm = TRUE))

sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinenttotalforestmedian) # 376.46 PgC
sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentforestlowerquantile) # 318.23
sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentforestupperquantile) # 441.85

sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinenttotalnonforestmedian) # 131.03
sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentnonforestlowerquantile) # 109.21
sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentnonforestupperquantile) # 166.26

sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentOVLmedian) # 91.8
sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentOVLlowerquantile) # 77.49
sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentOVLupperquantile) # 116.39

sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentextraforestmedian) # 33.43 PgC
sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentextraforestlowerquantile) # 28.64
sum(EcosystemCarbonStocksLandCategorySummaryStat$subcontinentextraforestupperquantile) # 65.69

#write.csv(EcosystemCarbonStocksLandCategorySummaryStat, 'PublicationCarbonStocksLandCategorySummary.csv')

## FIGURE 2: Total CS----                                                

setwd("N:/H73700/data/!pers/Manan/Paper3/JLUSRevision")

GlobalBestGuess <- data.frame(subcontinent = 'Global',
                              subcontinentmedian = sum(EcosystemCarbonWithShiftingCultivation$totalecosystemcarbonstocks_PgC),
                              value = 'Subcontinent')

EcosystemCarbonBestGuess <- EcosystemCarbonSubcontinentWithShiftingCultivation %>% # 5_EcosystemCarbonStocks.R
  mutate(value = 'Subcontinent') %>%
  rename(subcontinentmedian = subcontinentalcarbonstocks_PgC) %>%
  dplyr::select(subcontinent, subcontinentmedian, value) %>%
  rbind(GlobalBestGuess) %>%
  mutate(subcontinent =
           factor(subcontinent, 
                  levels = c('Global',
                             'Central America & the Caribbean',
                             'Eastern Africa',
                             'Eastern Asia',
                             'Eastern Europe',
                             'Northern Africa & Western Asia',
                             'Northern America',
                             'Oceania',
                             'Southeastern Asia',
                             'Southern Africa',
                             'Southern America',
                             'Southern Asia',
                             'Soviet Union',
                             'Western Africa',
                             'Western Europe'),
                  labels = c('Global',
                             'CAmC',
                             'EAf',
                             'EAs',
                             'EEu',
                             'NAWA',
                             'NAm',
                             'Oc',
                             'SEAs',
                             'SAf',
                             'SAm',
                             'SAs',
                             'SU',
                             'WAf',
                             'WEu'))) %>%
  filter(subcontinent != 'Global')

EcosystemCarbonStocksSummaryStatFigure <- EcosystemCarbonStocksSummaryStat %>%
  dplyr::select(subcontinent, subcontinentmedian, subcontinentlowerquantile, subcontinentupperquantile)

(GlobalCarbonStocks <- 
    data.frame(subcontinent = 'Global',
               subcontinentmedian = sum(EcosystemCarbonStocksSummaryStat$subcontinentmedian),
               subcontinentlowerquantile = sum(EcosystemCarbonStocksSummaryStat$subcontinentlowerquantile),
               subcontinentupperquantile = sum(EcosystemCarbonStocksSummaryStat$subcontinentupperquantile)) %>%
    rbind(EcosystemCarbonStocksSummaryStatFigure) %>%
    mutate(value = 'Subcontinent') %>%
    mutate(subcontinent =
             factor(subcontinent, 
                    levels = c('Global',
                               'Central America & the Caribbean',
                               'Eastern Africa',
                               'Eastern Asia',
                               'Eastern Europe',
                               'Northern Africa & Western Asia',
                               'Northern America',
                               'Oceania',
                               'Southeastern Asia',
                               'Southern Africa',
                               'Southern America',
                               'Southern Asia',
                               'Soviet Union',
                               'Western Africa',
                               'Western Europe'),
                    labels = c('Global',
                               'CAmC',
                               'EAf',
                               'EAs',
                               'EEu',
                               'NAWA',
                               'NAm',
                               'Oc',
                               'SEAs',
                               'SAf',
                               'SAm',
                               'SAs',
                               'SU',
                               'WAf',
                               'WEu'))) %>%
    filter(subcontinent != 'Global') %>%
    ggplot(aes(x = subcontinent, y = subcontinentmedian)) +
    #geom_bar(stat = 'identity', size = 0.5, fill = 'grey90', colour = 'black') +
    geom_bar(stat = 'identity', size = 0.5, aes(fill = subcontinent)) +
    geom_errorbar(aes(ymin = subcontinentlowerquantile, ymax = subcontinentupperquantile), size = 1, width = 0.1) +
    theme_bw() +
    geom_point(data = EcosystemCarbonBestGuess, aes(x = subcontinent, y = subcontinentmedian), colour = "#0072B2", size = 3) + # blue
    #facet_wrap(~ subcontinent, ncol = 5, scales = "free_y", labeller = label_wrap_gen(width = 20)) +
    scale_y_continuous(name = "Carbon stocks (PgC)", limits = c(0, 200)) + # , limits = c(0,NA)
    #scale_fill_manual(values = c("#0072B2", "#D55E00")) +
    scale_fill_manual(values = c(#"grey50", 
      "#004949","#009292","#ff6db6","#ffb6db",
      "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
      "#920000","#924900","#db6d00","#24ff24","#ffff6d")) +
    guides(fill = 'none') +
    theme(axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = 'bold', size = 10, color = "black"),
          strip.text.x = element_text(size = 10, face = 'bold'),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))


#ggsave(filename = 'Figure 2 Bar Graph.jpg', GlobalCarbonStocks, dpi = 300, width = 8.5, height = 5, units = 'in')


# Pie chart contribution----

(EcosystemCarbonPieChart <- EcosystemCarbonSubcontinentWithShiftingCultivation %>% # 5_EcosystemCarbonStocks.R
   mutate(value = 'Subcontinent') %>%
   rename(subcontinentmedian = subcontinentalcarbonstocks_PgC) %>%
   dplyr::select(subcontinent, subcontinentmedian, value) %>%
   rbind(GlobalBestGuess) %>%
   mutate(subcontinent =
            factor(subcontinent, 
                   levels = c('Global',
                              'Central America & the Caribbean',
                              'Eastern Africa',
                              'Eastern Asia',
                              'Eastern Europe',
                              'Northern Africa & Western Asia',
                              'Northern America',
                              'Oceania',
                              'Southeastern Asia',
                              'Southern Africa',
                              'Southern America',
                              'Southern Asia',
                              'Soviet Union',
                              'Western Africa',
                              'Western Europe'),
                   labels = c('Global',
                              'CAmC',
                              'EAf',
                              'EAs',
                              'EEu',
                              'NAWA',
                              'NAm',
                              'Oc',
                              'SEAs',
                              'SAf',
                              'SAm',
                              'SAs',
                              'SU',
                              'WAf',
                              'WEu'))) %>%
   filter(subcontinent !='Global') %>%
   mutate(subcontinentmedian = round(subcontinentmedian, 3)) %>%
   ggplot(aes(x = "", y = subcontinentmedian, fill = subcontinent)) +
   geom_bar(width = 1, stat = "identity") +
   scale_fill_manual(values = c("#004949","#009292","#ff6db6","#ffb6db",
                                "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                                "#920000","#924900","#db6d00","#24ff24","#ffff6d")) +
   coord_polar("y", start = 0) + 
   theme_bw() +
   #geom_text(aes(label = ifelse(subcontinentmedian > 40, (subcontinentmedian*100)/sum(subcontinentmedian), '')), hjust = 0.3, size = 3, position = position_stack(vjust = 0.5)) +
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(size = 10, color = "black"),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         axis.ticks = element_blank(),
         panel.border = element_blank(),
         legend.title = element_blank(),
         strip.text.x = element_text(size = 10, face = 'bold'),
         legend.text = element_text(size = 10, face = 'bold', color = "black"),
         legend.position = 'none',
         legend.background = element_rect(fill = "grey90")))

#ggsave(filename = 'Figure 2 Pie Chart.jpg', EcosystemCarbonPieChart, dpi = 300, width = 5, height = 5, units = 'in')
#rm(GlobalCarbonStocks, EcosystemCarbonStocksSummaryStatFigure, GlobalBestGuess, EcosystemCarbonBestGuess)


## FIGURE 3: Pie-chart----

GlobalDistribution <-
  data.frame(subcontinent = 'Global',
             landcategory = c('Forest', 
                              'Other vegetated land', 
                              'Rangeland', 
                              'Crop, pasture and infrastructure',
                              'Shifting cultivation'),
             estimate = c(sum(EcosystemCarbonLandCategoryWithShiftingCultivation$totalforestcarbonstocks_PgC)/
                            sum(EcosystemCarbonWithShiftingCultivation$totalecosystemcarbonstocks_PgC),
                          
                          sum(EcosystemCarbonLandCategoryWithShiftingCultivation$OVLcarbonstocks_PgC)/
                            sum(EcosystemCarbonWithShiftingCultivation$totalecosystemcarbonstocks_PgC),
                          
                          sum(EcosystemCarbonLandCategoryWithShiftingCultivation$rangelandscarbonstocks_PgC)/
                            sum(EcosystemCarbonWithShiftingCultivation$totalecosystemcarbonstocks_PgC),
                          
                          sum(EcosystemCarbonLandCategoryWithShiftingCultivation$croppastureinfrastructurecarbonstocks_PgC)/
                            sum(EcosystemCarbonWithShiftingCultivation$totalecosystemcarbonstocks_PgC),
                          
                          sum(EcosystemCarbonLandCategoryWithShiftingCultivation$SCcarbonstocks_PgC/
                                sum(EcosystemCarbonWithShiftingCultivation$totalecosystemcarbonstocks_PgC))))

(SubcontinentContributionCarbonStocks <- 
    EcosystemCarbonWithShiftingCultivation %>% # Best-guess
    mutate(totalforestcarbonstocks_MtC = 
             accessibleforestcarbonstocks_MtC + inaccessibleforestcarbonstocks_MtC + extraforestcarbonstocks_MtC) %>%
    mutate(croppastureinfrastructure_MtC = 
             infrastructurecarbonstocks_MtC + permanentcropscarbonstocks_MtC + annualcropscarbonstocks_MtC) %>%
    dplyr::select(subcontinent, totalforestcarbonstocks_MtC, OVLcarbonstocks_MtC, 
                  rangelandscarbonstocks_MtC, croppastureinfrastructure_MtC, SCcarbonstocks_MtC) %>%
    
    group_by(subcontinent) %>%
    summarise(totalforestcarbonstocks_PgC = sum(totalforestcarbonstocks_MtC, na.rm = TRUE)/1000,
              OVLcarbonstocks_PgC = sum(OVLcarbonstocks_MtC, na.rm = TRUE)/1000,
              rangelandscarbonstocks_PgC = sum(rangelandscarbonstocks_MtC, na.rm = TRUE)/1000,
              SCcarbonstocks_PgC = sum(SCcarbonstocks_MtC, na.rm = TRUE)/1000,
              croppastureinfrastructure_PgC = sum(croppastureinfrastructure_MtC, na.rm = TRUE)/1000) %>%
    mutate(forestfraction = totalforestcarbonstocks_PgC/
             (totalforestcarbonstocks_PgC + OVLcarbonstocks_PgC + rangelandscarbonstocks_PgC + SCcarbonstocks_PgC + croppastureinfrastructure_PgC),
           OVLfraction = OVLcarbonstocks_PgC/
             (totalforestcarbonstocks_PgC + OVLcarbonstocks_PgC + rangelandscarbonstocks_PgC + SCcarbonstocks_PgC + croppastureinfrastructure_PgC),
           rangelandfraction = rangelandscarbonstocks_PgC/
             (totalforestcarbonstocks_PgC + OVLcarbonstocks_PgC + rangelandscarbonstocks_PgC + SCcarbonstocks_PgC + croppastureinfrastructure_PgC),
           SCfraction = SCcarbonstocks_PgC/
             (totalforestcarbonstocks_PgC + OVLcarbonstocks_PgC + rangelandscarbonstocks_PgC + SCcarbonstocks_PgC + croppastureinfrastructure_PgC),
           croppastureinfrastructurefraction = croppastureinfrastructure_PgC/
             (totalforestcarbonstocks_PgC + OVLcarbonstocks_PgC + rangelandscarbonstocks_PgC + SCcarbonstocks_PgC + croppastureinfrastructure_PgC)) %>%
    
    dplyr::select(subcontinent, forestfraction, OVLfraction, rangelandfraction, SCfraction, croppastureinfrastructurefraction) %>%
    
    mutate(forestfraction = round(forestfraction, 3),
           OVLfraction = round(OVLfraction, 3), 
           rangelandfraction = round(rangelandfraction, 3), 
           SCfraction = round(SCfraction, 3),
           croppastureinfrastructurefraction = round(croppastureinfrastructurefraction, 3)) %>%
    gather(landcategory, estimate, forestfraction:croppastureinfrastructurefraction) %>%
    
    # Name all land categories for legend
    mutate(landcategory = factor(landcategory, 
                                 levels = c('forestfraction', 
                                            'OVLfraction',
                                            'rangelandfraction', 
                                            'croppastureinfrastructurefraction',
                                            'SCfraction'),
                                 labels = c('Forest', 
                                            'Other vegetated land', 
                                            'Rangeland', 
                                            'Crop, pasture and infrastructure',
                                            'Shifting cultivation'))) %>%
    
    rbind(GlobalDistribution) %>%
    
    mutate(estimate = round(estimate, 3)) %>%
    
    mutate(subcontinent =
             factor(subcontinent, 
                    levels = c('Global',
                               'Central America & the Caribbean',
                               'Eastern Africa',
                               'Eastern Asia',
                               'Eastern Europe',
                               'Northern Africa & Western Asia',
                               'Northern America',
                               'Oceania',
                               'Southeastern Asia',
                               'Southern Africa',
                               'Southern America',
                               'Southern Asia',
                               'Soviet Union',
                               'Western Africa',
                               'Western Europe'),
                    labels = c('Global',
                               'Central America & the Caribbean',
                               'Eastern Africa',
                               'Eastern Asia',
                               'Eastern Europe',
                               'Northern Africa & Western Asia',
                               'Northern America',
                               'Oceania',
                               'Southeastern Asia',
                               'Southern Africa',
                               'Southern America',
                               'Southern Asia',
                               'Soviet Union',
                               'Western Africa',
                               'Western Europe'))) %>%
    
    # Plot
    ggplot(aes(x = "", y = estimate, fill = landcategory)) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = c("#009E73", "#E69F00", "#56B4E9", "#F0E442", "#0072B2")) +
    coord_polar("y", start = 0) + 
    facet_wrap(~ subcontinent, ncol = 5, labeller = label_wrap_gen(width = 20)) +
    theme_bw() +
    geom_text(aes(label = ifelse(estimate > 0.1, percent(estimate), '')), hjust = 0.3, size = 3, position = position_stack(vjust = 0.5)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 10, face = 'bold'),
          legend.text = element_text(size = 10, face = 'bold', color = "black"),
          legend.position = 'bottom',
          legend.background = element_rect(fill = "grey90")))

#ggsave(filename = 'Figure 3.jpg', SubcontinentContributionCarbonStocks, dpi = 300, width = 11, height = 8.5, units = 'in')

rm(GlobalDistribution, SubcontinentContributionCarbonStocks)


## FIGURE 4: Total land distribution----

SubcontinentCarbonStocksDistribution <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  dplyr::select(-ecosystemcarbonstocks_PgC) %>%
  group_by(subcontinent) %>%
  summarise(accessibleforestmedian = median(accessibleforestcarbonstocks_PgC, na.rm = TRUE),
            inaccessibleforestmedian = median(inaccessibleforestcarbonstocks_PgC, na.rm = TRUE),
            extraforestmedian = median(extraforestcarbonstocks_PgC, na.rm = TRUE),
            OVLmedian = median(OVLcarbonstocks_PgC, na.rm = TRUE),
            rangelandmedian = median(rangelandcarbonstocks_PgC, na.rm = TRUE),
            shiftingcultivationmedianmedian = median(carbonstocksSC_PgC, na.rm = TRUE),
            croppastureinframedian = median(croppastureinfrastructure_PgC, na.rm = TRUE)) %>%
  gather(carbonstockstype, estimate, accessibleforestmedian:croppastureinframedian) %>%
  mutate(carbonstockstype =
           factor(carbonstockstype, 
                  levels = c('accessibleforestmedian',
                             'inaccessibleforestmedian', 
                             'extraforestmedian',
                             'OVLmedian',
                             'rangelandmedian',
                             'croppastureinframedian',
                             'shiftingcultivationmedianmedian'),
                  labels = c('Ac\\nFor',
                             'Inac\\nFor',
                             'Ext\\nFor',
                             'OVL',
                             'Ran',
                             'CrP\\nInfra',
                             'SC')))

BestGuessLandCategorySubcontinentCarbonStocks <- 
  EcosystemCarbonLandCategoryWithShiftingCultivation %>% # 5_EcosystemCarbonStocks.R
  dplyr::select(c(1:8)) %>%
  mutate(shiftingcultivationcarbonstocks_PgC = SCcarbonstocks_PgC) %>%
  dplyr::select(-7) %>%
  mutate(focalcase = 'Best-guess') %>%
  gather(carbonstockstype, estimate, accessibleforestcarbonstocks_PgC:shiftingcultivationcarbonstocks_PgC) %>%
  mutate(carbonstockstype = as.character(carbonstockstype),
         estimate = as.numeric(estimate)) %>%
  mutate(carbonstockstype =
           factor(carbonstockstype, 
                  levels = c('accessibleforestcarbonstocks_PgC',
                             'inaccessibleforestcarbonstocks_PgC', 
                             'extraforestcarbonstocks_PgC',
                             'OVLcarbonstocks_PgC',
                             'rangelandscarbonstocks_PgC',
                             'croppastureinfrastructurecarbonstocks_PgC',
                             'shiftingcultivationcarbonstocks_PgC'),
                  labels = c('Ac\\nFor',
                             'Inac\\nFor',
                             'Ext\\nFor',
                             'OVL',
                             'Ran',
                             'CrP\\nInfra',
                             'SC')))

(LandCategorySubcontinentCarbonStocks <- 
    rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
    group_by(subcontinent) %>%
    #mutate(totalforestcarbonstocks_PgC = accessibleforestcarbonstocks_PgC + inaccessibleforestcarbonstocks_PgC + extraforestcarbonstocks_PgC) %>%
    dplyr::select(-c(ecosystemcarbonstocks_PgC)) %>%
    rename(shiftingcultivationcarbonstocks_PgC = carbonstocksSC_PgC) %>%
    gather(carbonstockstype, estimate, accessibleforestcarbonstocks_PgC:shiftingcultivationcarbonstocks_PgC) %>%
    mutate(carbonstockstype =
             factor(carbonstockstype, 
                    levels = c('accessibleforestcarbonstocks_PgC',
                               'inaccessibleforestcarbonstocks_PgC', 
                               'extraforestcarbonstocks_PgC',
                               'OVLcarbonstocks_PgC',
                               'rangelandcarbonstocks_PgC',
                               'croppastureinfrastructure_PgC',
                               'shiftingcultivationcarbonstocks_PgC'),
                    labels = c('Ac\\nFor',
                               'Inac\\nFor',
                               'Ext\\nFor',
                               'OVL',
                               'Ran',
                               'CrP\\nInfra',
                               'SC'))) %>%
    #group_by(subcontinent, carbonstockstype) %>%
    #summarize(subcontinentcarbonstockstypemax = max(estimate), subcontinentcarbonstockstypemin = min(estimate)) %>%
    #gather(carbonstockstype, estimate, subcontinentcarbonstockstypemax: subcontinentcarbonstockstypemin) %>%
    
    ggplot(aes(x = carbonstockstype, y = estimate)) +
    #geom_errorbar(aes(ymax = subcontinentcarbonstockstypemax, ymin = subcontinentcarbonstockstypemin)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    scale_y_continuous(name = "Carbon stocks (PgC)") +
    
    geom_point(data = BestGuessLandCategorySubcontinentCarbonStocks, aes(x = carbonstockstype, y = estimate), colour = "#0072B2", size = 3) +
    facet_wrap(~ subcontinent, ncol = 5, scales = "free_y", labeller = label_wrap_gen(width = 20)) +# blue
    #geom_point(data = SubcontinentCarbonStocksDistribution, aes(x = carbonstockstype, y = estimate), colour = "#D55E00", size = 3) + # mustard
    theme(axis.text.x = element_text(face = 'bold', size = 8, color = "black"),
          axis.text.y = element_text(face = 'bold', size = 8, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, face = 'bold'),
          legend.title = element_blank(), legend.position = 'bottom',
          strip.text.x = element_text(size = 10, face = 'bold'),
          plot.title = element_text(face = 'bold', size = 10, family = '', hjust = 0.5, vjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

# Blue is without SC
# Mustard is with SC
#ggsave(filename = 'Figure 4 Subcontinent Stocks.jpg', LandCategorySubcontinentCarbonStocks, dpi = 300, width = 11, height = 8.5, units = 'in')

## FIGURE 4 Inset: Global distribution PgC----

GlobalCarbonStocksDistributionMedian <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  dplyr::select(-ecosystemcarbonstocks_PgC) %>%
  group_by(subcontinent) %>%
  summarise(subcontinentaccessibleforestmedian = median(accessibleforestcarbonstocks_PgC, na.rm = TRUE),
            subcontinentinaccessibleforestmedian = median(inaccessibleforestcarbonstocks_PgC, na.rm = TRUE),
            subcontinentextraforestmedian = median(extraforestcarbonstocks_PgC, na.rm = TRUE),
            subcontinentOVLmedian = median(OVLcarbonstocks_PgC, na.rm = TRUE),
            subcontinentrangelandmedian = median(rangelandcarbonstocks_PgC, na.rm = TRUE),
            subcontinentcroppastureinframedian = median(croppastureinfrastructure_PgC, na.rm = TRUE),
            subcontinentSCmedian = median(carbonstocksSC_PgC, na.rm = TRUE)) %>%
  summarise(AcFor = sum(subcontinentaccessibleforestmedian),
            InacFor = sum(subcontinentinaccessibleforestmedian),
            ExtFor = sum(subcontinentextraforestmedian),
            OVL = sum(subcontinentOVLmedian),
            Ran = sum(subcontinentrangelandmedian),
            CrPInfra = sum(subcontinentcroppastureinframedian),
            SC = sum(subcontinentSCmedian)) %>%
  mutate(subcontinent = 'Global') %>%
  gather(landcategory, median, AcFor:SC)

GlobalCarbonStocksDistributionLowerQuantile <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  dplyr::select(-ecosystemcarbonstocks_PgC) %>%
  group_by(subcontinent) %>%
  summarise(subcontinentaccessibleforestmedian = quantile(accessibleforestcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentinaccessibleforestmedian = quantile(inaccessibleforestcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentextraforestmedian = quantile(extraforestcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentOVLmedian = quantile(OVLcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentrangelandmedian = quantile(rangelandcarbonstocks_PgC, p = 0.25, na.rm = TRUE),
            subcontinentcroppastureinframedian = quantile(croppastureinfrastructure_PgC, p = 0.25, na.rm = TRUE),
            subcontinentSCmedian = quantile(carbonstocksSC_PgC, p = 0.25, na.rm = TRUE)) %>%
  summarise(AcFor = sum(subcontinentaccessibleforestmedian),
            InacFor = sum(subcontinentinaccessibleforestmedian),
            ExtFor = sum(subcontinentextraforestmedian),
            OVL = sum(subcontinentOVLmedian),
            Ran = sum(subcontinentrangelandmedian),
            CrPInfra = sum(subcontinentcroppastureinframedian),
            SC = sum(subcontinentSCmedian)) %>%
  mutate(subcontinent = 'Global') %>%
  gather(landcategory, lowerquantile, AcFor:SC)

GlobalCarbonStocksDistributionUpperQuantile <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  dplyr::select(-ecosystemcarbonstocks_PgC) %>%
  group_by(subcontinent) %>%
  summarise(subcontinentaccessibleforestmedian = quantile(accessibleforestcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            subcontinentinaccessibleforestmedian = quantile(inaccessibleforestcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            subcontinentextraforestmedian = quantile(extraforestcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            subcontinentOVLmedian = quantile(OVLcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            subcontinentrangelandmedian = quantile(rangelandcarbonstocks_PgC, p = 0.75, na.rm = TRUE),
            subcontinentcroppastureinframedian = quantile(croppastureinfrastructure_PgC, p = 0.75, na.rm = TRUE),
            subcontinentSCmedian = quantile(carbonstocksSC_PgC, p = 0.75, na.rm = TRUE)) %>%
  summarise(AcFor = sum(subcontinentaccessibleforestmedian),
            InacFor = sum(subcontinentinaccessibleforestmedian),
            ExtFor = sum(subcontinentextraforestmedian),
            OVL = sum(subcontinentOVLmedian),
            Ran = sum(subcontinentrangelandmedian),
            CrPInfra = sum(subcontinentcroppastureinframedian),
            SC = sum(subcontinentSCmedian)) %>%
  mutate(subcontinent = 'Global') %>%
  gather(landcategory, upperquantile, AcFor:SC)


BestGuessDistribution <- 
  data.frame(landcategory = c('AcFor', 'InacFor', 'ExtFor', 'OVL', 'Ran', 'CrPInfra', 'SC'),
             estimate = c(sum(EcosystemCarbonWithShiftingCultivation$accessibleforestcarbonstocks_MtC)/1000,
                          sum(EcosystemCarbonWithShiftingCultivation$inaccessibleforestcarbonstocks_MtC)/1000,
                          sum(EcosystemCarbonWithShiftingCultivation$extraforestcarbonstocks_MtC)/1000,
                          sum(EcosystemCarbonWithShiftingCultivation$OVLcarbonstocks_MtC)/1000,
                          sum(EcosystemCarbonWithShiftingCultivation$rangelandscarbonstocks_MtC)/1000,
                          (sum(EcosystemCarbonWithShiftingCultivation$annualcropscarbonstocks_MtC) +
                             sum(EcosystemCarbonWithShiftingCultivation$permanentcropscarbonstocks_MtC) +
                             sum(EcosystemCarbonWithShiftingCultivation$pasturescarbonstocks_MtC) +
                             sum(EcosystemCarbonWithShiftingCultivation$infrastructurecarbonstocks_MtC))/1000,
                          sum(EcosystemCarbonWithShiftingCultivation$SCcarbonstocks_MtC)/1000)) %>%
  mutate(landcategory = factor(landcategory,
                               levels = c('AcFor', 'InacFor', 'ExtFor', 'OVL', 'Ran', 'CrPInfra', 'SC'),
                               labels = c('Ac\\nFor', 'Inac\\nFor', 'Ext\\nFor', 'OVL', 'Ran', 'CrP\\nInfra', 'SC')))

(GlobalDistribution <- GlobalCarbonStocksDistributionMedian %>%
    left_join(GlobalCarbonStocksDistributionLowerQuantile, by = c('subcontinent', 'landcategory')) %>%
    left_join(GlobalCarbonStocksDistributionUpperQuantile, by = c('subcontinent', 'landcategory')) %>%
    mutate(landcategory = factor(landcategory,
                                 levels = c('AcFor', 'InacFor', 'ExtFor', 'OVL', 'Ran', 'CrPInfra', 'SC'),
                                 labels = c('Ac\\nFor', 'Inac\\nFor', 'Ext\\nFor', 'OVL', 'Ran', 'CrP\\nInfra', 'SC'))) %>%
    ggplot(aes(x = landcategory, y = median)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lowerquantile, ymax = upperquantile), size = 0.5, width = 0.1) +
    scale_y_continuous(name = "Carbon stocks (PgC)", limits = c(0, 300), breaks = c(0, 50, 100, 150, 200, 250, 300)) +
    scale_fill_manual(values = c("#0072B2", "#D55E00")) +
    theme_bw() +
    geom_point(data = BestGuessDistribution, aes(x = landcategory, y = estimate), colour = "#0072B2", size = 3) + 
    theme(axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
          axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = 'bold', size = 10, color = "black"),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(), legend.position = c(0.8, 0.75),
          legend.text = element_text(face = 'bold', size = 12, color = "black"),
          legend.background = element_rect(fill = "grey90"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 

#ggsave(filename = 'Figure 4 Global Stocks.jpg', GlobalDistribution, dpi = 300, width = 3, height = 3, units = 'in')


## SI FIGURE 4: Subcontinent tC/ha----

LandDistribution <- LandUseEstimatesWithShiftingCultivation %>%
  left_join(SubContinentalClassification, by = c('country')) %>%
  mutate(croppastureinfra_Mha = cropland_Mha + pastures_Mha + infrastructure_Mha,
         SC_Mha = shortfallow_Mha + longfallow_Mha) %>%
  group_by(subcontinent) %>%
  summarize(subcontinentaccessibleforest_Mha = sum(accessibleforest_Mha, na.rm = TRUE),
            subcontinentinaccessibleforest_Mha = sum(inaccessibleforest_Mha, na.rm = TRUE),
            subcontinentextraforest_Mha = sum(extraforest_Mha, na.rm = TRUE),
            subcontinentOVL_Mha = sum(othervegetatedland_Mha, na.rm = TRUE),
            subcontinentrangeland_Mha = sum(rangelands_Mha, na.rm = TRUE),
            subcontinentshiftingcultivation_Mha = sum(SC_Mha, na.rm = TRUE),
            subcontinentcroppastureinfra_Mha = sum(croppastureinfra_Mha, na.rm = TRUE))


BestGuessLandCategorySubcontinentCarbonDensity <- 
  EcosystemCarbonLandCategoryWithShiftingCultivation %>% # 5_EcosystemCarbonStocks.R
  dplyr::select(c(1:8)) %>%
  mutate(shiftingcultivationcarbonstocks_PgC = SCcarbonstocks_PgC) %>%
  dplyr::select(-7) %>%
  mutate(focalcase = 'Best-guess') %>%
  left_join(LandDistribution, by = c('subcontinent')) %>%
  mutate(accessibleforest_tCha = (accessibleforestcarbonstocks_PgC * 1000)/subcontinentaccessibleforest_Mha,
         inaccessibleforest_tCha = (inaccessibleforestcarbonstocks_PgC * 1000)/subcontinentinaccessibleforest_Mha,
         extraforest_tCha = (extraforestcarbonstocks_PgC * 1000)/subcontinentextraforest_Mha,
         OVL_tCha = (OVLcarbonstocks_PgC * 1000)/subcontinentOVL_Mha,
         rangeland_tCha = (rangelandscarbonstocks_PgC * 1000)/subcontinentrangeland_Mha,
         SC_tCha = (shiftingcultivationcarbonstocks_PgC * 1000)/subcontinentshiftingcultivation_Mha,
         croppastureinfra_tCha = (croppastureinfrastructurecarbonstocks_PgC * 1000)/subcontinentcroppastureinfra_Mha) %>%
  dplyr::select(subcontinent, accessibleforest_tCha, inaccessibleforest_tCha,
                extraforest_tCha, OVL_tCha, rangeland_tCha, SC_tCha, croppastureinfra_tCha) %>%
  gather(carbondensitytype, estimate, accessibleforest_tCha:croppastureinfra_tCha) %>%
  mutate(carbondensitytype = as.character(carbondensitytype),
         estimate = as.numeric(estimate)) %>%
  mutate(estimate = ifelse(is.nan(estimate), 0, estimate)) %>%
  mutate(carbondensitytype =
           factor(carbondensitytype, 
                  levels = c('accessibleforest_tCha',
                             'inaccessibleforest_tCha', 
                             'extraforest_tCha',
                             'OVL_tCha',
                             'rangeland_tCha',
                             'croppastureinfra_tCha',
                             'SC_tCha'),
                  labels = c('Ac\\nFor',
                             'Inac\\nFor',
                             'Ext\\nFor',
                             'OVL',
                             'Ran',
                             'CrP\\nInfra',
                             'SC')))

(LandCategorySubcontinentCarbonDensity <- 
    rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
    #group_by(subcontinent) %>%
    #mutate(totalforestcarbonstocks_PgC = accessibleforestcarbonstocks_PgC + inaccessibleforestcarbonstocks_PgC + extraforestcarbonstocks_PgC) %>%
    dplyr::select(-c(ecosystemcarbonstocks_PgC)) %>%
    left_join(LandDistribution, by = c('subcontinent')) %>%
    mutate(accessibleforest_tCha = (accessibleforestcarbonstocks_PgC * 1000)/subcontinentaccessibleforest_Mha,
           inaccessibleforest_tCha = (inaccessibleforestcarbonstocks_PgC * 1000)/subcontinentinaccessibleforest_Mha,
           extraforest_tCha = (extraforestcarbonstocks_PgC * 1000)/subcontinentextraforest_Mha,
           OVL_tCha = (OVLcarbonstocks_PgC * 1000)/subcontinentOVL_Mha,
           rangeland_tCha = (rangelandcarbonstocks_PgC * 1000)/subcontinentrangeland_Mha,
           SC_tCha = (carbonstocksSC_PgC * 1000)/subcontinentshiftingcultivation_Mha,
           croppastureinfra_tCha = (croppastureinfrastructure_PgC * 1000)/subcontinentcroppastureinfra_Mha) %>%
    dplyr::select(subcontinent, accessibleforest_tCha, inaccessibleforest_tCha,
                  extraforest_tCha, OVL_tCha, rangeland_tCha, SC_tCha, croppastureinfra_tCha) %>%
    gather(carbondensitytype, estimate, accessibleforest_tCha:croppastureinfra_tCha) %>%
    mutate(carbondensitytype = as.character(carbondensitytype),
           estimate = as.numeric(estimate)) %>%
    mutate(estimate = ifelse(is.nan(estimate), 0, estimate)) %>%
    mutate(carbondensitytype =
             factor(carbondensitytype, 
                    levels = c('accessibleforest_tCha',
                               'inaccessibleforest_tCha', 
                               'extraforest_tCha',
                               'OVL_tCha',
                               'rangeland_tCha',
                               'croppastureinfra_tCha',
                               'SC_tCha'),
                    labels = c('Ac\\nFor',
                               'Inac\\nFor',
                               'Ext\\nFor',
                               'OVL',
                               'Ran',
                               'CrP\\nInfra',
                               'SC'))) %>%
    
    #group_by(subcontinent, carbonstockstype) %>%
    #summarize(subcontinentcarbonstockstypemax = max(estimate), subcontinentcarbonstockstypemin = min(estimate)) %>%
    #gather(carbonstockstype, estimate, subcontinentcarbonstockstypemax: subcontinentcarbonstockstypemin) %>%
    
    ggplot(aes(x = carbondensitytype, y = estimate)) +
    #geom_errorbar(aes(ymax = subcontinentcarbonstockstypemax, ymin = subcontinentcarbonstockstypemin)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    scale_y_continuous(name = "Carbon density (tC/ha)") +
    
    geom_point(data = BestGuessLandCategorySubcontinentCarbonDensity, aes(x = carbondensitytype, y = estimate), colour = "#0072B2", size = 3) +
    facet_wrap(~ subcontinent, ncol = 5, 
               scales = "free_y", 
               labeller = label_wrap_gen(width = 20)) +# blue
    #geom_point(data = FC2LandCategorySubcontinentCarbonStocks, aes(x = carbonstockstype, y = estimate), colour = "#D55E00", size = 3) + # mustard
    theme(axis.text.x = element_text(face = 'bold', size = 8, color = "black"),
          axis.text.y = element_text(face = 'bold', size = 8, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, face = 'bold'),
          legend.title = element_blank(), legend.position = 'bottom',
          strip.text.x = element_text(size = 10, face = 'bold'),
          plot.title = element_text(face = 'bold', size = 10, family = '', hjust = 0.5, vjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

# Blue is without SC
# Mustard is with SC
#ggsave(filename = 'SI Figure 4 Subcontinent.jpg', LandCategorySubcontinentCarbonDensity, dpi = 300, width = 11, height = 8.5, units = 'in')


## SI FIGURE 4: Global tC/ha----

GlobalCarbonStocksDistributionMedian <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  dplyr::select(-ecosystemcarbonstocks_PgC) %>%
  left_join(LandDistribution, by = c('subcontinent')) %>%
  mutate(accessibleforest_tCha = (accessibleforestcarbonstocks_PgC * 1000)/subcontinentaccessibleforest_Mha,
         inaccessibleforest_tCha = (inaccessibleforestcarbonstocks_PgC * 1000)/subcontinentinaccessibleforest_Mha,
         extraforest_tCha = (extraforestcarbonstocks_PgC * 1000)/subcontinentextraforest_Mha,
         OVL_tCha = (OVLcarbonstocks_PgC * 1000)/subcontinentOVL_Mha,
         rangeland_tCha = (rangelandcarbonstocks_PgC * 1000)/subcontinentrangeland_Mha,
         SC_tCha = (carbonstocksSC_PgC * 1000)/subcontinentshiftingcultivation_Mha,
         croppastureinfra_tCha = (croppastureinfrastructure_PgC * 1000)/subcontinentcroppastureinfra_Mha) %>%
  dplyr::select(subcontinent, accessibleforest_tCha, inaccessibleforest_tCha,
                extraforest_tCha, OVL_tCha, rangeland_tCha, SC_tCha, croppastureinfra_tCha) %>%
  group_by(subcontinent) %>%
  summarise(subcontinentmedianaccessibleforest_tCha = median(accessibleforest_tCha, na.rm = TRUE),
            subcontinentmedianinaccessibleforest_tCha = median(inaccessibleforest_tCha, na.rm = TRUE),
            subcontinentmedianextraforest_tCha = median(extraforest_tCha, na.rm = TRUE),
            subcontinentmedianOVL_tCha = median(OVL_tCha, na.rm = TRUE),
            subcontinentmedianrangeland_tCha = median(rangeland_tCha, na.rm = TRUE),
            subcontinentmediancroppastureinfra_tCha = median(croppastureinfra_tCha, na.rm = TRUE),
            subcontinentmedianSC_tCha = median(SC_tCha, na.rm = TRUE)) %>%
  #replace(is.na(.), 0)
  summarise(AcFor = mean(subcontinentmedianaccessibleforest_tCha, na.rm = TRUE),
            InacFor = mean(subcontinentmedianinaccessibleforest_tCha, na.rm = TRUE),
            ExtFor = mean(subcontinentmedianextraforest_tCha, na.rm = TRUE),
            OVL = mean(subcontinentmedianOVL_tCha, na.rm = TRUE),
            Ran = mean(subcontinentmedianrangeland_tCha, na.rm = TRUE),
            CrPInfra = mean(subcontinentmediancroppastureinfra_tCha, na.rm = TRUE),
            SC = mean(subcontinentmedianSC_tCha, na.rm = TRUE)) %>%
  mutate(subcontinent = 'Global') %>%
  gather(landcategory, median, AcFor:SC)

GlobalCarbonStocksDistributionLowerQuantile <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  dplyr::select(-ecosystemcarbonstocks_PgC) %>%
  left_join(LandDistribution, by = c('subcontinent')) %>%
  mutate(accessibleforest_tCha = (accessibleforestcarbonstocks_PgC * 1000)/subcontinentaccessibleforest_Mha,
         inaccessibleforest_tCha = (inaccessibleforestcarbonstocks_PgC * 1000)/subcontinentinaccessibleforest_Mha,
         extraforest_tCha = (extraforestcarbonstocks_PgC * 1000)/subcontinentextraforest_Mha,
         OVL_tCha = (OVLcarbonstocks_PgC * 1000)/subcontinentOVL_Mha,
         rangeland_tCha = (rangelandcarbonstocks_PgC * 1000)/subcontinentrangeland_Mha,
         SC_tCha = (carbonstocksSC_PgC * 1000)/subcontinentshiftingcultivation_Mha,
         croppastureinfra_tCha = (croppastureinfrastructure_PgC * 1000)/subcontinentcroppastureinfra_Mha) %>%
  dplyr::select(subcontinent, accessibleforest_tCha, inaccessibleforest_tCha,
                extraforest_tCha, OVL_tCha, rangeland_tCha, SC_tCha, croppastureinfra_tCha) %>%
  group_by(subcontinent) %>%
  summarise(subcontinentquantileaccessibleforest_tCha = quantile(accessibleforest_tCha, p = 0.25, na.rm = TRUE),
            subcontinentquantileinaccessibleforest_tCha = quantile(inaccessibleforest_tCha, p = 0.25, na.rm = TRUE),
            subcontinentquantileextraforest_tCha = quantile(extraforest_tCha, p = 0.25, na.rm = TRUE),
            subcontinentquantileOVL_tCha = quantile(OVL_tCha, p = 0.25, na.rm = TRUE),
            subcontinentquantilerangeland_tCha = quantile(rangeland_tCha, p = 0.25, na.rm = TRUE),
            subcontinentquantilecroppastureinfra_tCha = quantile(croppastureinfra_tCha, p = 0.25, na.rm = TRUE),
            subcontinentquantileSC_tCha = quantile(SC_tCha, p = 0.25, na.rm = TRUE)) %>%
  #replace(is.na(.), 0)
  summarise(AcFor = mean(subcontinentquantileaccessibleforest_tCha, na.rm = TRUE),
            InacFor = mean(subcontinentquantileinaccessibleforest_tCha, na.rm = TRUE),
            ExtFor = mean(subcontinentquantileextraforest_tCha, na.rm = TRUE),
            OVL = mean(subcontinentquantileOVL_tCha, na.rm = TRUE),
            Ran = mean(subcontinentquantilerangeland_tCha, na.rm = TRUE),
            CrPInfra = mean(subcontinentquantilecroppastureinfra_tCha, na.rm = TRUE),
            SC = mean(subcontinentquantileSC_tCha, na.rm = TRUE)) %>%
  mutate(subcontinent = 'Global') %>%
  gather(landcategory, lowerquantile, AcFor:SC)


GlobalCarbonStocksDistributionUpperQuantile <- 
  rbind(EcosystemCarbonStocksWithoutSC, EcosystemCarbonStocksWithSC) %>%
  dplyr::select(-ecosystemcarbonstocks_PgC) %>%
  left_join(LandDistribution, by = c('subcontinent')) %>%
  mutate(accessibleforest_tCha = (accessibleforestcarbonstocks_PgC * 1000)/subcontinentaccessibleforest_Mha,
         inaccessibleforest_tCha = (inaccessibleforestcarbonstocks_PgC * 1000)/subcontinentinaccessibleforest_Mha,
         extraforest_tCha = (extraforestcarbonstocks_PgC * 1000)/subcontinentextraforest_Mha,
         OVL_tCha = (OVLcarbonstocks_PgC * 1000)/subcontinentOVL_Mha,
         rangeland_tCha = (rangelandcarbonstocks_PgC * 1000)/subcontinentrangeland_Mha,
         SC_tCha = (carbonstocksSC_PgC * 1000)/subcontinentshiftingcultivation_Mha,
         croppastureinfra_tCha = (croppastureinfrastructure_PgC * 1000)/subcontinentcroppastureinfra_Mha) %>%
  dplyr::select(subcontinent, accessibleforest_tCha, inaccessibleforest_tCha,
                extraforest_tCha, OVL_tCha, rangeland_tCha, SC_tCha, croppastureinfra_tCha) %>%
  group_by(subcontinent) %>%
  summarise(subcontinentquantileaccessibleforest_tCha = quantile(accessibleforest_tCha, p = 0.75, na.rm = TRUE),
            subcontinentquantileinaccessibleforest_tCha = quantile(inaccessibleforest_tCha, p = 0.75, na.rm = TRUE),
            subcontinentquantileextraforest_tCha = quantile(extraforest_tCha, p = 0.75, na.rm = TRUE),
            subcontinentquantileOVL_tCha = quantile(OVL_tCha, p = 0.75, na.rm = TRUE),
            subcontinentquantilerangeland_tCha = quantile(rangeland_tCha, p = 0.75, na.rm = TRUE),
            subcontinentquantilecroppastureinfra_tCha = quantile(croppastureinfra_tCha, p = 0.75, na.rm = TRUE),
            subcontinentquantileSC_tCha = quantile(SC_tCha, p = 0.75, na.rm = TRUE)) %>%
  #replace(is.na(.), 0)
  summarise(AcFor = mean(subcontinentquantileaccessibleforest_tCha, na.rm = TRUE),
            InacFor = mean(subcontinentquantileinaccessibleforest_tCha, na.rm = TRUE),
            ExtFor = mean(subcontinentquantileextraforest_tCha, na.rm = TRUE),
            OVL = mean(subcontinentquantileOVL_tCha, na.rm = TRUE),
            Ran = mean(subcontinentquantilerangeland_tCha, na.rm = TRUE),
            CrPInfra = mean(subcontinentquantilecroppastureinfra_tCha, na.rm = TRUE),
            SC = mean(subcontinentquantileSC_tCha, na.rm = TRUE)) %>%
  mutate(subcontinent = 'Global') %>%
  gather(landcategory, upperquantile, AcFor:SC)



BestGuessDistribution <- 
  data.frame(landcategory = c('AcFor', 'InacFor', 'ExtFor', 'OVL', 'Ran', 'CrPInfra', 'SC'),
             estimate = c(mean(sum(EcosystemCarbonWithShiftingCultivation$accessibleforestcarbonstocks_MtC)/sum(EcosystemCarbonWithShiftingCultivation$accessibleforest_Mha)),
                          mean(sum(EcosystemCarbonWithShiftingCultivation$inaccessibleforestcarbonstocks_MtC)/sum(EcosystemCarbonWithShiftingCultivation$inaccessibleforest_Mha)),
                          mean(sum(EcosystemCarbonWithShiftingCultivation$extraforestcarbonstocks_MtC)/sum(EcosystemCarbonWithShiftingCultivation$extraforest_Mha)),
                          mean(sum(EcosystemCarbonWithShiftingCultivation$OVLcarbonstocks_MtC)/sum(EcosystemCarbonWithShiftingCultivation$othervegetatedland_Mha)),
                          mean(sum(EcosystemCarbonWithShiftingCultivation$rangelandscarbonstocks_MtC)/sum(EcosystemCarbonWithShiftingCultivation$rangelands_Mha)),
                          mean((sum(EcosystemCarbonWithShiftingCultivation$annualcropscarbonstocks_MtC) + 
                                  sum(EcosystemCarbonWithShiftingCultivation$permanentcropscarbonstocks_MtC) +
                                  sum(EcosystemCarbonWithShiftingCultivation$pasturescarbonstocks_MtC) +
                                  sum(EcosystemCarbonWithShiftingCultivation$infrastructurecarbonstocks_MtC))/
                                 (sum(EcosystemCarbonWithShiftingCultivation$AnnualCrops_Mha) + 
                                    sum(EcosystemCarbonWithShiftingCultivation$PermanentCrops_Mha) +
                                    sum(EcosystemCarbonWithShiftingCultivation$pastures_Mha) +
                                    sum(EcosystemCarbonWithShiftingCultivation$infrastructure_Mha))),
                          mean(sum(EcosystemCarbonWithShiftingCultivation$SCcarbonstocks_MtC)/
                                 (sum(EcosystemCarbonWithShiftingCultivation$shortfallow_Mha) + sum(EcosystemCarbonWithShiftingCultivation$longfallow_Mha))))) %>%
  mutate(landcategory = factor(landcategory,
                               levels = c('AcFor', 'InacFor', 'ExtFor', 'OVL', 'Ran', 'CrPInfra', 'SC'),
                               labels = c('Ac\\nFor', 'Inac\\nFor', 'Ext\\nFor', 'OVL', 'Ran', 'CrP\\nInfra', 'SC')))

(GlobalDistribution <- GlobalCarbonStocksDistributionMedian %>%
    left_join(GlobalCarbonStocksDistributionLowerQuantile, by = c('subcontinent', 'landcategory')) %>%
    left_join(GlobalCarbonStocksDistributionUpperQuantile, by = c('subcontinent', 'landcategory')) %>%
    mutate(landcategory = factor(landcategory,
                                 levels = c('AcFor', 'InacFor', 'ExtFor', 'OVL', 'Ran', 'CrPInfra', 'SC'),
                                 labels = c('Ac\\nFor', 'Inac\\nFor', 'Ext\\nFor', 'OVL', 'Ran', 'CrP\\nInfra', 'SC'))) %>%
    ggplot(aes(x = landcategory, y = median)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lowerquantile, ymax = upperquantile), size = 0.5, width = 0.1) +
    scale_y_continuous(name = "Carbon density (tC/ha)", limits = c(0, 250), breaks = c(0, 50, 100, 150, 200, 250)) +
    scale_fill_manual(values = c("#0072B2", "#D55E00")) +
    theme_bw() +
    geom_point(data = BestGuessDistribution, aes(x = landcategory, y = estimate), colour = "#0072B2", size = 3) + 
    theme(axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
          axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = 'bold', size = 10, color = "black"),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(), legend.position = c(0.8, 0.75),
          legend.text = element_text(face = 'bold', size = 12, color = "black"),
          legend.background = element_rect(fill = "grey90"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())) 


#ggsave(filename = 'SI Figure 4 Global Inset.jpg', GlobalDistribution, dpi = 300, width = 3, height = 3, units = 'in')

## FIGURE 5: Forest area vs C density----

ForestAreaPerc <- LandUseEstimatesWithShiftingCultivation %>%
  left_join(SubContinentalClassification, by = c('country')) %>%
  mutate(totalforestarea_Mha = accessibleforest_Mha + inaccessibleforest_Mha + extraforest_Mha) %>%
  dplyr::select(country, subcontinent, landarea_Mha, totalforestarea_Mha) %>%
  group_by(subcontinent) %>%
  summarise(landarea_Mha = sum(landarea_Mha),
            forestarea_Mha = sum(totalforestarea_Mha)) %>%
  mutate(forestareashare = forestarea_Mha/landarea_Mha) %>%
  dplyr::select(subcontinent, landarea_Mha, forestareashare)

ForestAreaCarbonFraction <- EcosystemCarbonWithShiftingCultivation %>% ## Focal Case 2
  mutate(forestcarbonstocks_PgC = (accessibleforestcarbonstocks_MtC + inaccessibleforestcarbonstocks_MtC + extraforestcarbonstocks_MtC)/1000) %>%
  group_by(subcontinent) %>%
  summarise(totalforestarea_Mha = sum(totalforest_Mha, na.rm = TRUE),
            totalforestcarbonstocks_PgC = sum(forestcarbonstocks_PgC, na.rm = TRUE)) %>%
  mutate(forestcarbondensity_tCha = (totalforestcarbonstocks_PgC/totalforestarea_Mha)*1000) %>%
  dplyr::select(subcontinent, totalforestcarbonstocks_PgC, forestcarbondensity_tCha)

PublicationForestAreaCarbonDensityShare <- ForestAreaPerc %>%
  left_join(ForestAreaCarbonFraction, by = c('subcontinent')) %>%
  mutate(subcontinentabbreviation = 
           case_when(subcontinent == 'Central America & the Caribbean' ~ 'CAmC',
                     subcontinent == 'Eastern Africa' ~ 'EAf',
                     subcontinent == 'Eastern Asia' ~ 'EAs',
                     subcontinent == 'Eastern Europe' ~ 'EEu',
                     subcontinent == 'Northern Africa & Western Asia' ~ 'NAWA',
                     subcontinent == 'Northern America' ~ 'NAm',
                     subcontinent == 'Oceania' ~ 'Oc',
                     subcontinent == 'Southeastern Asia' ~ 'SEAs',
                     subcontinent == 'Southern Africa' ~ 'SAf',
                     subcontinent == 'Southern America' ~ 'SAm',
                     subcontinent == 'Southern Asia' ~ 'SAs',
                     subcontinent == 'Soviet Union' ~ 'SU',
                     subcontinent == 'Western Africa' ~ 'WAf',
                     subcontinent == 'Western Europe' ~ 'WEu'))

#write.csv(PublicationForestAreaCarbonDensityShare, 'PublicationForestAreaForestCarbonDensity.csv')
       

(Figure6 <- ForestAreaPerc %>%
    left_join(ForestAreaCarbonFraction, by = c('subcontinent')) %>%
    mutate(subcontinentabbreviation = 
             case_when(subcontinent == 'Central America & the Caribbean' ~ 'CAmC',
                       subcontinent == 'Eastern Africa' ~ 'EAf',
                       subcontinent == 'Eastern Asia' ~ 'EAs',
                       subcontinent == 'Eastern Europe' ~ 'EEu',
                       subcontinent == 'Northern Africa & Western Asia' ~ 'NAWA',
                       subcontinent == 'Northern America' ~ 'NAm',
                       subcontinent == 'Oceania' ~ 'Oc',
                       subcontinent == 'Southeastern Asia' ~ 'SEAs',
                       subcontinent == 'Southern Africa' ~ 'SAf',
                       subcontinent == 'Southern America' ~ 'SAm',
                       subcontinent == 'Southern Asia' ~ 'SAs',
                       subcontinent == 'Soviet Union' ~ 'SU',
                       subcontinent == 'Western Africa' ~ 'WAf',
                       subcontinent == 'Western Europe' ~ 'WEu')) %>%
    ggplot(aes(x = forestcarbondensity_tCha, y = forestareashare, label = subcontinentabbreviation)) +
    geom_point(aes(size = totalforestcarbonstocks_PgC, colour = subcontinent, stroke = 0.5), alpha = 1) +
    scale_colour_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", 
                                   "#999933", "#882255", "#661100", "#6699CC", "#888888", "#F0E442", "#CC79A7")) +
    guides(size = 'none', colour = 'none') +
    theme_bw() +
    #geom_text(check_overlap = TRUE, hjust = 1.3, vjust = 0.5) +
    scale_size(range = c(2, 10)) +
    scale_x_continuous(name = 'Forest carbon density (tC/ha)', limits = c(0, 170), breaks = c(0, 50, 100, 150)) +
    scale_y_continuous(name = 'Fraction of land area as forest', limits = c(0, 0.7), breaks = c(0, 0.35, 0.7)) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
          axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
          axis.title.x = element_text(size = 10, face = 'bold'),
          axis.title.y = element_text(size = 10, face = 'bold')))

# With label
#ggsave(filename = 'Figure 5.jpg', Figure6, dpi = 300, width = 5, height = 5.5, units = 'in')

# Without label
#ggsave(filename = 'Figure 5 Without Label.jpg', Figure6, dpi = 300, width = 5, height = 5.5, units = 'in')

## FIGURE 6: COUNTRY-LEVEL----

## FIGURE 6A: Top 3 per subcontinent----

Top3Subcontinent <-
  EcosystemCarbonWithShiftingCultivation %>%
  dplyr::select(country, subcontinent, totalecosystemcarbonstocks_PgC) %>% # Focal Case 2
  group_by(subcontinent) %>%
  arrange(desc(totalecosystemcarbonstocks_PgC)) %>%
  mutate(Index = 1:n(), country = ifelse(Index > 3, "Other", country)) %>%
  ungroup() %>%
  group_by(country, subcontinent) %>%
  summarise(totalecosystemcarbonstocks_PgC = sum(totalecosystemcarbonstocks_PgC)) %>%
  group_by(subcontinent) %>%
  mutate(fraction = totalecosystemcarbonstocks_PgC/sum(totalecosystemcarbonstocks_PgC))

#write.csv(Top3Subcontinent, 'PublicationTop3CountriesPerSubcontinent.csv')

unique(Top3Subcontinent$country)
unique(Top3Subcontinent$subcontinent)
table(Top3Subcontinent$country)


(Top3PerSubcontinent <- Top3Subcontinent %>%
    mutate(country = factor(country, 
                            levels = c("Honduras", "Mexico", "Nicaragua", 
                                       "Anglo-Egyptian Sudan", "Ethiopia", "Kenya",
                                       "China", 'Japan', "Mongol. Peopl. Rep.",
                                       'Poland', 'Rumania', 'Yugoslavia',
                                       "French Morocco", 'Iran', 'Turkey',
                                       'Canada', 'United States',
                                       'Australia', "New Guinea (Austr.)", 'New Zealand',
                                       'Indonesia', 'Malaya', 'Burma',
                                       'Angola', 'Madagascar', 'Northern Rhodesia',
                                       'Brazil', 'Colombia', 'Peru',
                                       'India', 'Pakistan', 'Ceylon',
                                       "U.S.S.R.",
                                       "Belgian Congo", "French Equat. Africa", "French West Africa",
                                       'Finland', 'France', 'Sweden', 
                                       'Other'))) %>%
    mutate(subcontinent = factor(subcontinent,
                                 levels = c("Western Europe",
                                            "Western Africa",
                                            "Soviet Union",
                                            "Southern Asia",
                                            "Southern America",
                                            "Southern Africa",
                                            "Southeastern Asia",
                                            "Oceania",
                                            "Northern America",
                                            "Northern Africa & Western Asia",
                                            "Eastern Europe",
                                            "Eastern Asia",
                                            "Eastern Africa",
                                            "Central America & the Caribbean"),
                                 labels = c("Western Europe (1.2%)",
                                            "Western Africa (16.6%)",
                                            "Soviet Union (12.7%)",
                                            "Southern Asia (2.5%)",
                                            "Southern America (27.7%)",
                                            "Southern Africa (4.1%)",
                                            "Southeastern Asia (9.7%)",
                                            "Oceania (3.8%)",
                                            "Northern America (10.5%)",
                                            "Northern Africa & Western Asia (0.9%)",
                                            "Eastern Europe (0.5%)",
                                            "Eastern Asia (4.0%)",
                                            "Eastern Africa (2.7%)",
                                            "Central America & the Caribbean (3.0%)"))) %>%
    
    mutate(color = if_else(country == 'Other', 'red', "white")) %>%
    #ungroup() %>%
    ggplot(aes(x = subcontinent, y = fraction, group = country, colour = 'black', fill = color)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = c('red' = "#999933", 'white' = 'white')) +
    scale_color_manual(values = c('black')) +
    labs(y = "Carbon stocks fraction") + 
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    #geom_text(aes(label = country), colour = "black", check_overlap = TRUE) +
    theme_bw() +
    theme(axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
          axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_text(face = 'bold', size = 10, color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_blank(), 
          legend.position = 'none',
          legend.background = element_rect(fill = alpha(0.4))) +
    coord_flip())

#ggsave(filename = 'Figure 6A.jpg', Top3PerSubcontinent, dpi = 300, width = 8.5, height = 5, units = 'in')

## FIGURE 6B: Top 20 countries----

PublicationTop20Countries <- 
  EcosystemCarbonWithShiftingCultivation %>% # Best-guess
  
  # Arrange in descending order
  arrange(desc(totalecosystemcarbonstocks_PgC)) %>%
  
  # Pick top 20 countries
  mutate(Index = 1:n(), country = ifelse(Index > 20, "Other", country)) %>%
  filter(!(country == 'Other')) %>%
  
  # Make forest and non-forest carbon stocks
  mutate(totalforestcarbonstocks_MtC = accessibleforestcarbonstocks_MtC +
           inaccessibleforestcarbonstocks_MtC +
           extraforestcarbonstocks_MtC,
         totalforestcarbonstocks_PgC = totalforestcarbonstocks_MtC/1000) %>%
  #mutate(othervegetatedlandcarbonstocks_PgC = OVLcarbonstocks_MtC/1000) %>%
  mutate(totalnonforestcarbonstocks_MtC = OVLcarbonstocks_MtC + 
           annualcropscarbonstocks_MtC +
           permanentcropscarbonstocks_MtC +
           rangelandscarbonstocks_MtC +
           pasturescarbonstocks_MtC +
           infrastructurecarbonstocks_MtC,
         totalnonforestcarbonstocks_PgC = totalnonforestcarbonstocks_MtC/1000,
         totalcarbonstocks_PgC = totalnonforestcarbonstocks_PgC + totalforestcarbonstocks_PgC,
         perccontribution = totalcarbonstocks_PgC/sum(totalcarbonstocks_PgC)) %>%
  dplyr::select(country, 
                totalforestcarbonstocks_PgC,
                #othervegetatedlandcarbonstocks_PgC,
                totalnonforestcarbonstocks_PgC,
                totalcarbonstocks_PgC,
                perccontribution) %>%
  gather(carbonstockstype, 
         estimate, 
         totalforestcarbonstocks_PgC, 
         #othervegetatedlandcarbonstocks_PgC, 
         totalnonforestcarbonstocks_PgC) %>%
  mutate(carbonstockstype = factor(carbonstockstype, 
                                   levels = c('totalforestcarbonstocks_PgC', 
                                              #'othervegetatedlandcarbonstocks_PgC',
                                              'totalnonforestcarbonstocks_PgC'),
                                   labels = c('Forest carbon stocks', 
                                              #'OVL carbon stocks', 
                                              'Non-forest carbon stocks'))) %>%
  mutate(country = if_else(country == 'New Guinea (Austr.)', 'Papua New Guinea',
                           if_else(country == 'Anglo-Egyptian Sudan', 'Sudan',
                                   if_else(country == 'Malaya', 'Malaysia', country)))) %>%
  mutate(country = factor(country,
                          levels = c('Brazil', 'U.S.S.R.', 'Belgian Congo', 'United States', 'Canada', 'Indonesia', 
                                     'French Equat. Africa', 'China', 'French West Africa', 'Peru', 'India', 'Colombia', 'Venezuela', 'Bolivia', 
                                     'Papua New Guinea', 'Australia', 'Mexico', 'Sudan', 'Angola', 'Malaysia')))

#write.csv(PublicationTop20Countries, 'PublicationTop20Countries.csv')

(EcosystemCarbonStocksTop20Countries <- 
    EcosystemCarbonWithShiftingCultivation %>% # Best-guess
    
    # Arrange in descending order
    arrange(desc(totalecosystemcarbonstocks_PgC)) %>%
    
    # Pick top 20 countries
    mutate(Index = 1:n(), country = ifelse(Index > 20, "Other", country)) %>%
    filter(!(country == 'Other')) %>%
    
    # Make forest and non-forest carbon stocks
    mutate(totalforestcarbonstocks_MtC = accessibleforestcarbonstocks_MtC +
             inaccessibleforestcarbonstocks_MtC +
             extraforestcarbonstocks_MtC,
           totalforestcarbonstocks_PgC = totalforestcarbonstocks_MtC/1000) %>%
    #mutate(othervegetatedlandcarbonstocks_PgC = OVLcarbonstocks_MtC/1000) %>%
    mutate(totalnonforestcarbonstocks_MtC = OVLcarbonstocks_MtC + 
             annualcropscarbonstocks_MtC +
             permanentcropscarbonstocks_MtC +
             rangelandscarbonstocks_MtC +
             pasturescarbonstocks_MtC +
             infrastructurecarbonstocks_MtC,
           totalnonforestcarbonstocks_PgC = totalnonforestcarbonstocks_MtC/1000,
           totalcarbonstocks_PgC = totalnonforestcarbonstocks_PgC + totalforestcarbonstocks_PgC,
           perccontribution = totalcarbonstocks_PgC/sum(totalcarbonstocks_PgC)) %>%
    dplyr::select(country, 
                  totalforestcarbonstocks_PgC,
                  #othervegetatedlandcarbonstocks_PgC,
                  totalnonforestcarbonstocks_PgC,
                  totalcarbonstocks_PgC,
                  perccontribution) %>%
    gather(carbonstockstype, 
           estimate, 
           totalforestcarbonstocks_PgC, 
           #othervegetatedlandcarbonstocks_PgC, 
           totalnonforestcarbonstocks_PgC) %>%
    mutate(carbonstockstype = factor(carbonstockstype, 
                                     levels = c('totalforestcarbonstocks_PgC', 
                                                #'othervegetatedlandcarbonstocks_PgC',
                                                'totalnonforestcarbonstocks_PgC'),
                                     labels = c('Forest carbon stocks', 
                                                #'OVL carbon stocks', 
                                                'Non-forest carbon stocks'))) %>%
    mutate(country = if_else(country == 'New Guinea (Austr.)', 'Papua New Guinea',
                             if_else(country == 'Anglo-Egyptian Sudan', 'Sudan',
                                     if_else(country == 'Malaya', 'Malaysia', country)))) %>%
    mutate(country = factor(country,
                            levels = c('Brazil', 'U.S.S.R.', 'Belgian Congo', 'United States', 'Canada', 'Indonesia', 
                                       'French Equat. Africa', 'China', 'French West Africa', 'Peru', 'India', 'Colombia', 'Venezuela', 'Bolivia', 
                                       'Papua New Guinea', 'Australia', 'Mexico', 'Sudan', 'Angola', 'Malaysia'))) %>%
    
    # Make plot
    ggplot(aes(x = country, y = estimate, 
               fill = carbonstockstype)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("#009E73", "#D55E00")) +
    theme_bw() +
    #geom_text(aes(label = totalcarbonstocks_PgC), vjust = -0.2) +
    #ggtitle('Top 20 countries with highest carbon stocks in 1950') +
    scale_y_continuous(name = "Carbon stocks (PgC)", limits = c(0, 80), breaks = c(0, 40, 80)) +
    theme(axis.text.x = element_text(face = 'bold', size = 10, color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10, face = 'bold'),
          legend.title = element_blank(), legend.position = c(0.6, 0.8),
          legend.text = element_text(face = 'bold', size = 10, color = "black"),
          legend.background = element_rect(fill = "grey90"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()))

#ggsave(filename = 'Figure 6B.jpg', EcosystemCarbonStocksTop20Countries, dpi = 300, width = 8.5, height = 5, units = 'in')

## 6A + 6B----

library(ggpubr)

(Figure6 <- ggarrange(Top3PerSubcontinent, EcosystemCarbonStocksTop20Countries, nrow = 2))

#ggsave(filename = 'Figure 6 without labels.jpg', Figure6, dpi = 300, width = 11, height = 8.5, units = 'in')

## SI FIGURE 7: Forest area/carbon fractions----

(SIForFigure6 <- EcosystemCarbonWithShiftingCultivation %>% ## Focal Case 2
   mutate(forestcarbonstocks_PgC = (accessibleforestcarbonstocks_MtC + inaccessibleforestcarbonstocks_MtC + extraforestcarbonstocks_MtC)/1000) %>%
   group_by(subcontinent) %>%
   summarise(totallandarea_Mha = sum(landarea_Mha, na.rm = TRUE),
             totalcarbonstocks_PgC = sum(totalecosystemcarbonstocks_PgC, na.rm = TRUE),
             totalforestarea_Mha = sum(totalforest_Mha, na.rm = TRUE),
             totalforestcarbonstocks_PgC = sum(forestcarbonstocks_PgC, na.rm = TRUE)) %>%
   mutate(forestcarbondensity_tCha = (totalforestcarbonstocks_PgC/totalforestarea_Mha)*1000,
          carbonstocksfraction = totalcarbonstocks_PgC/sum(totalcarbonstocks_PgC),
          forestareafraction = totalforestarea_Mha/sum(totalforestarea_Mha)) %>%
   mutate(subcontinentabbreviation = 
            case_when(subcontinent == 'Central America & the Caribbean' ~ 'CAmC',
                      subcontinent == 'Eastern Africa' ~ 'EAf',
                      subcontinent == 'Eastern Asia' ~ 'EAs',
                      subcontinent == 'Eastern Europe' ~ 'EEu',
                      subcontinent == 'Northern Africa & Western Asia' ~ 'NAWA',
                      subcontinent == 'Northern America' ~ 'NAm',
                      subcontinent == 'Oceania' ~ 'Oc',
                      subcontinent == 'Southeastern Asia' ~ 'SEAs',
                      subcontinent == 'Southern Africa' ~ 'SAf',
                      subcontinent == 'Southern America' ~ 'SAm',
                      subcontinent == 'Southern Asia' ~ 'SAs',
                      subcontinent == 'Soviet Union' ~ 'SU',
                      subcontinent == 'Western Africa' ~ 'WAf',
                      subcontinent == 'Western Europe' ~ 'WEu')) %>%
   ggplot(aes(x = forestcarbondensity_tCha, y = carbonstocksfraction, label = subcontinentabbreviation)) +
   geom_point(aes(size = totallandarea_Mha, colour = subcontinent, stroke = 0.5), alpha = 1) +
   scale_colour_manual(values = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", 
                                  "#999933", "#882255", "#661100", "#6699CC", "#888888", "#F0E442", "#CC79A7")) +
   guides(size = 'none', colour = 'none') +
   theme_bw() +
   geom_text(check_overlap = TRUE, hjust = - 0.5, vjust = 0.5) +
   scale_size(range = c(4, 8)) +
   scale_x_continuous(name = 'Forest carbon stocks (tC/ha)', limits = c(0, 170), breaks = c(0, 50, 100, 150)) +
   scale_y_continuous(name = 'Carbon stocks fraction', limits = c(0, 0.3), breaks = c(0, 0.1, 0.2, 0.3)) +
   theme(panel.grid = element_blank(),
         axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
         axis.text.y = element_text(face = 'bold', size = 10, color = "black"),
         axis.title.x = element_text(size = 10, face = 'bold'),
         axis.title.y = element_text(size = 10, face = 'bold')))

#ggsave(filename = 'SI Figure 7.jpg', SIForFigure6, dpi = 300, width = 5, height = 5.5, units = 'in')


## FIGURE 1: Conceptual----

## Estimate----

EcosystemCarbonWithoutSC <- LandUseEstimates %>%
  
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
  left_join(CarbonDensitiesWhereActGreaterThanPot, by = c('country')) %>%
  
  # OVL carbon density
  mutate(OVLcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   forestcarbondensity_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Inaccessible and extra forest carbon density
  mutate(nongrowingstockcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   potentialforestcarbon_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Forest CS
  mutate(accessibleforestcarbonstocks_MtC = (accessibleforest_Mha * forestcarbondensity_tCha),
         inaccessibleforestcarbonstocks_MtC = (inaccessibleforest_Mha * nongrowingstockcarbondensity_tCha),
         extraforestcarbonstocks_MtC = (extraforest_Mha * nongrowingstockcarbondensity_tCha)) %>%
  
  # Changing inaccessible forest CS for some countries
  mutate(inaccessibleforestcarbonstocks_MtC = if_else(country == 'Algeria' |
                                                        country == 'Australia' |
                                                        country == 'Austria' |
                                                        country == 'Canada' |
                                                        #country == 'China' |
                                                        country == 'Finland',
                                                      0.75 * (inaccessibleforest_Mha * forestcarbondensity_tCha), # can also assume it only reaches 75% as a measure of unproductive forests 
                                                      inaccessibleforestcarbonstocks_MtC)) %>%
  
  
  # Cropland CS
  mutate(annualcropscarbonstocks_MtC = (AnnualCrops_Mha * totalNPPact1950_tCha),
         permanentcropscarbonstocks_MtC = (0.5 * PermanentCrops_Mha) * 15 + (0.5 * PermanentCrops_Mha) * 30) %>%
  
  # OVL CS
  mutate(OVLcarbonstocks_MtC = othervegetatedland_Mha * (OVLcarbondensity_tCha/2)) %>%
  
  # rangeland CS
  mutate(rangelandscarbonstocks_MtC = rangelands_Mha * (NPPactgrazing1950_tCha * 2)) %>%
  
  # pasture CS
  mutate(pasturescarbonstocks_MtC = pastures_Mha * NPPactgrazing1950_tCha) %>%
  
  # infrastructure CS
  mutate(infrastructurecarbonstocks_MtC = infrastructure_Mha * (potential_tCha/6)) %>%
  
  # total ecosystem CS
  mutate(totalecosystemcarbonstocks_MtC = 
           accessibleforestcarbonstocks_MtC + inaccessibleforestcarbonstocks_MtC + extraforestcarbonstocks_MtC +
           annualcropscarbonstocks_MtC + permanentcropscarbonstocks_MtC +
           OVLcarbonstocks_MtC + 
           rangelandscarbonstocks_MtC +
           pasturescarbonstocks_MtC +
           infrastructurecarbonstocks_MtC,
         totalecosystemcarbonstocks_PgC = totalecosystemcarbonstocks_MtC/1000)

sum(EcosystemCarbonWithoutSC$totalecosystemcarbonstocks_PgC)

EcosystemCarbonWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
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
  left_join(SCCarbonDensities, by = c('country')) %>%
  left_join(CarbonDensitiesWhereActGreaterThanPot, by = c('country')) %>%
  left_join(NPPGrazing1950, by = c('subcontinent')) %>%
  left_join(NPPCropland1950, by = c('subcontinent')) %>%
  
  # OVL carbon density
  mutate(OVLcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   forestcarbondensity_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Inaccessible and extra forest carbon density
  mutate(nongrowingstockcarbondensity_tCha = 
           if_else(is.na(nongrowingstockcarbondensity_tCha),
                   potentialforestcarbon_tCha,
                   nongrowingstockcarbondensity_tCha)) %>%
  
  # Forest CS
  mutate(accessibleforestcarbonstocks_MtC = (accessibleforest_Mha * forestcarbondensity_tCha),
         inaccessibleforestcarbonstocks_MtC = (inaccessibleforest_Mha * nongrowingstockcarbondensity_tCha),
         extraforestcarbonstocks_MtC = (extraforest_Mha * nongrowingstockcarbondensity_tCha)) %>%
  
  # Changing inaccessible forest CS for some countries
  mutate(inaccessibleforestcarbonstocks_MtC = if_else(country == 'Algeria' |
                                                        country == 'Australia' |
                                                        country == 'Austria' |
                                                        country == 'Canada' |
                                                        country == 'China' |
                                                        country == 'Finland',
                                                      0.75 * (inaccessibleforest_Mha * forestcarbondensity_tCha), # can also assume it only reaches 75% as a measure of unproductive forests 
                                                      inaccessibleforestcarbonstocks_MtC)) %>%
  
  
  # Shifting cultivation
  mutate(totalC_shortfallow_tCha = if_else(is.na(totalC_shortfallow_tCha), 0, totalC_shortfallow_tCha),
         totalC_longfallow_tCha = if_else(is.na(totalC_longfallow_tCha), 0, totalC_longfallow_tCha)) %>%
  
  mutate(SCcarbonstocks_MtC = (shortfallow_Mha * totalC_shortfallow_tCha) +
           (longfallow_Mha * totalC_longfallow_tCha)) %>%
  
  # Cropland CS
  mutate(annualcropscarbonstocks_MtC = (AnnualCrops_Mha * totalNPPact1950_tCha),
         permanentcropscarbonstocks_MtC = (0.5 * PermanentCrops_Mha) * 15 + (0.5 * PermanentCrops_Mha) * 30) %>%
  
  # OVL CS
  mutate(OVLcarbonstocks_MtC = othervegetatedland_Mha * (OVLcarbondensity_tCha/2)) %>%
  
  # rangeland CS
  mutate(rangelandscarbonstocks_MtC = rangelands_Mha * (NPPactgrazing1950_tCha * 2)) %>%
  
  # pasture CS
  mutate(pasturescarbonstocks_MtC = pastures_Mha * NPPactgrazing1950_tCha) %>%
  
  # infrastructure CS
  mutate(infrastructurecarbonstocks_MtC = infrastructure_Mha * (potential_tCha/6)) %>%
  
  # total ecosystem CS
  mutate(totalecosystemcarbonstocks_MtC = 
           accessibleforestcarbonstocks_MtC + inaccessibleforestcarbonstocks_MtC + extraforestcarbonstocks_MtC +
           SCcarbonstocks_MtC +
           annualcropscarbonstocks_MtC + permanentcropscarbonstocks_MtC +
           OVLcarbonstocks_MtC +
           rangelandscarbonstocks_MtC +
           pasturescarbonstocks_MtC + 
           infrastructurecarbonstocks_MtC,
         totalecosystemcarbonstocks_PgC = totalecosystemcarbonstocks_MtC/1000)

sum(EcosystemCarbonWithSC$totalecosystemcarbonstocks_PgC) # 454 PgC

sum(EcosystemCarbonWithSC$AnnualCrops_Mha) # 1154 Mha
sum(EcosystemCarbonWithSC$PermanentCrops_Mha) # 67.1
sum(EcosystemCarbonWithoutSC$accessibleforest_Mha) - sum(EcosystemCarbonWithSC$accessibleforest_Mha) # 17.3
sum(EcosystemCarbonWithoutSC$inaccessibleforest_Mha) - sum(EcosystemCarbonWithSC$inaccessibleforest_Mha) # 56.5
sum(EcosystemCarbonWithoutSC$extraforest_Mha) - sum(EcosystemCarbonWithSC$extraforest_Mha) # 39.3
sum(EcosystemCarbonWithoutSC$othervegetatedland_Mha) - sum(EcosystemCarbonWithSC$othervegetatedland_Mha) # 140

## Plot----

(Figure1Horizontal <- 
   read_xlsx("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/Figure1Conceptual.xlsx") %>%
   mutate(percentage100 = percentage * 100) %>%
   mutate(plot = 'Plot') %>%
   mutate(label = factor(label,
                         levels = c('Accessible Forest',
                                    'SC in Accessible Forest',
                                    'Inaccessible Forest',
                                    'SC in Inaccessible Forest',
                                    'Extra Forest',
                                    'SC in Extra Forest',
                                    'Annual Cropland',
                                    'Permanent Cropland',
                                    'Pasture',
                                    'Rangeland',
                                    'Infrastructure',
                                    'Non-Productive Land',
                                    'Other Vegetated Land',
                                    'SC in Other Vegetated Land'))) %>%
   ggplot(aes(x = plot, y = percentage100)) +
   geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.5, aes(fill = label)) +
   coord_flip() +
   theme_bw() +
   scale_y_continuous(name = "Terrestrial land area (%)", limits = c(0, 100), breaks = seq(0, 100, 20)) + 
   scale_fill_manual(values = c("#117733", "#DDCC77", "#999933", "#DDCC77", "#44AA99", '#DDCC77',
                                '#F0E442', '#E69F00', '#56B4E9', '#D55E00', '#661100', "#999999", 
                                '#CC6677', "#DDCC77")) +
   theme(axis.text.x = element_text(face = 'bold', size = 10, color = "black"),
         axis.text.y = element_blank(),
         axis.title.x = element_text(face = 'bold', size = 10, color = "black"),
         axis.title.y = element_blank(),
         legend.title = element_blank(), 
         legend.position = 'none',
         #legend.text = element_text(face = 'bold', size = 6, color = "black"),
         axis.ticks = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank()))

#ggsave(filename = 'Figure 1.jpg', Figure1Horizontal, dpi = 300, width = 5, height = 2.5, units = 'in') 

## FOREST MANAGEMENT IMPACT----

AccessibleForestWithSC <- LandUseEstimatesWithShiftingCultivation %>%
  
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
  
  # Forest CS
  mutate(accessibleforestcarbonstocksactual_MtC = (accessibleforest_Mha * forestcarbondensity_tCha),
         accessibleforestcarbonstocksactual_PgC = accessibleforestcarbonstocksactual_MtC/1000,
         accessibleforestcarbonstockspotential_MtC = (accessibleforest_Mha * potentialforestcarbon_tCha),
         accessibleforestcarbonstockspotential_PgC = accessibleforestcarbonstockspotential_MtC/1000,) %>%
  
  group_by(subcontinent) %>%
  summarise(accessibleforestsubcontinentactual = sum(accessibleforestcarbonstocksactual_PgC, na.rm = TRUE),
            accessibleforestsubcontinentpotential = sum(accessibleforestcarbonstockspotential_PgC, na.rm = TRUE)) %>%
  mutate(forestmanagementimpact = (accessibleforestsubcontinentpotential- accessibleforestsubcontinentactual),
         forestmanagementimpactproportion = (accessibleforestsubcontinentpotential- accessibleforestsubcontinentactual)/accessibleforestsubcontinentpotential)


sum(AccessibleForestWithSC$accessibleforestsubcontinentactual) # 129.28 Mha
sum(AccessibleForestWithSC$accessibleforestsubcontinentpotential) # 175.58 PgC
