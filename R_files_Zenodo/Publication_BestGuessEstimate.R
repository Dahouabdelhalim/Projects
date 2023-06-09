## A mid-20th century estimate of global vegetation carbon stocks based on inventory statistics
## Authors: Bhan et al., Journal of Land Use Science
## Best-guess biomass carbon stocks dataset

# Land use areas: "N:\\H73700\\data\\!pers\\Manan\\Paper3\\CarbonStockTimeseries\\DataFiles\\1_CountriesListLandUseCarbonStocks_FAO1953HANPP.R"
# Forest carbon stock densities: "N:\\H73700\\data\\!pers\\Manan\\Paper3\\CarbonStockTimeseries\\DataFiles\\1_FRAForestCarbonDensitiesCountriesAndAveraged.R"

setwd("N:/H73700/data/!pers/Manan/Paper3/JLUSRevision")

library(dplyr)
library(magrittr)
library(tidyr)
library(readxl)

options(scipen = 999)
options(digits = 5)

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

sum(LandUseEstimates$accessibleforest_Mha) # 1779.6
sum(LandUseEstimates$inaccessibleforest_Mha) # 2021.9
sum(LandUseEstimates$extraforest_Mha) # 282.98
#4091

#write.csv(LandUseEstimates, 'CountriesListLandUseEstimates.csv')

ExtraForest <- LandUseEstimates %>%
  dplyr::select(country, extraforest_Mha) %>%
  filter(extraforest_Mha > 0)

#write.csv(ExtraForest, 'CountriesListLandUseEstimatesExtraForest.csv')

sum(LandUseEstimates$infrastructure_Mha) # 78.6
sum(LandUseEstimates$nps_Mha)
sum(LandUseEstimates$cropland_Mha)
sum(LandUseEstimates$rangelands_Mha)
sum(LandUseEstimates$pastures_Mha)
sum(LandUseEstimates$othervegetatedland_Mha) # 3327

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

sum(LandUseEstimatesWithShiftingCultivation$accessibleforest_Mha) # 1762.3
sum(LandUseEstimatesWithShiftingCultivation$inaccessibleforest_Mha) # 1965.4
sum(LandUseEstimatesWithShiftingCultivation$extraforest_Mha) # 243.73
# 3971.43

ExtraForest <- LandUseEstimatesWithShiftingCultivation %>%
  dplyr::select(country, extraforest_Mha)

sum(LandUseEstimatesWithShiftingCultivation$shortfallow_Mha) # 120.18 Mha
sum(LandUseEstimatesWithShiftingCultivation$longfallow_Mha) # 138.7 Mha
# SC = 258.88 Mha
sum(LandUseEstimatesWithShiftingCultivation$infrastructure_Mha) # 78.58 Mha
sum(LandUseEstimatesWithShiftingCultivation$AnnualCrops_Mha) # 1153.9 Mha
sum(LandUseEstimatesWithShiftingCultivation$PermanentCrops_Mha) # 67.1 Mha
sum(LandUseEstimatesWithShiftingCultivation$pastures_Mha) # 528.8 Mha
sum(LandUseEstimatesWithShiftingCultivation$rangelands_Mha) # 1762.1 Mha

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

NPPGrazing1950 <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_NPPactGrazingLand1950_Kastner2021.csv") %>%
  dplyr::select(subcontinent, NPPactgrazing1950_tCha)

NPPCropland1950 <-
  read_xlsx("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/NPP1950Krausmann2000Haberl.xlsx") %>%
  dplyr::select(subcontinent, totalNPPact1950_tCha)


# BEST-GUESS----

EcosystemCarbonWithShiftingCultivation <- LandUseEstimatesWithShiftingCultivation %>%
  
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

sum(EcosystemCarbonWithShiftingCultivation$totalecosystemcarbonstocks_PgC) # 450.21 PgC

## Data for publication----

PublicationBestGuess <- EcosystemCarbonWithShiftingCultivation %>%
 dplyr::select(country, subcontinent, accessibleforestcarbonstocks_MtC, inaccessibleforestcarbonstocks_MtC, extraforestcarbonstocks_MtC,
                SCcarbonstocks_MtC, annualcropscarbonstocks_MtC, permanentcropscarbonstocks_MtC,
               OVLcarbonstocks_MtC, rangelandscarbonstocks_MtC, pasturescarbonstocks_MtC, infrastructurecarbonstocks_MtC,
              totalecosystemcarbonstocks_MtC)

#write.csv(PublicationBestGuess, 'PublicationBestGuessEstimate.csv')

## PER SUBCONTINENT----

EcosystemCarbonSubcontinentWithShiftingCultivation <- EcosystemCarbonWithShiftingCultivation %>%
  group_by(subcontinent) %>%
  summarize(subcontinentalcarbonstocks_PgC = sum(totalecosystemcarbonstocks_PgC)) %>%
  mutate(subcontinentalcarbonstocks_perc = subcontinentalcarbonstocks_PgC * 100/sum(subcontinentalcarbonstocks_PgC))

## PER LAND CATEGORY----

EcosystemCarbonLandCategoryWithShiftingCultivation <- EcosystemCarbonWithShiftingCultivation %>%
  group_by(subcontinent) %>%
  summarize(accessibleforestcarbonstocks_PgC = sum(accessibleforestcarbonstocks_MtC)/1000,
            inaccessibleforestcarbonstocks_PgC = sum(inaccessibleforestcarbonstocks_MtC)/1000,
            extraforestcarbonstocks_PgC = sum(extraforestcarbonstocks_MtC)/1000,
            OVLcarbonstocks_PgC = sum(OVLcarbonstocks_MtC)/1000,
            rangelandscarbonstocks_PgC = sum(rangelandscarbonstocks_MtC)/1000,
            SCcarbonstocks_PgC = sum(SCcarbonstocks_MtC)/1000,
            croppastureinfrastructurecarbonstocks_PgC = 
              (sum(annualcropscarbonstocks_MtC) + sum(permanentcropscarbonstocks_MtC) + 
                 sum(pasturescarbonstocks_MtC) + sum(infrastructurecarbonstocks_MtC))/1000) %>%
  mutate(totalforestcarbonstocks_PgC = 
           accessibleforestcarbonstocks_PgC + inaccessibleforestcarbonstocks_PgC + extraforestcarbonstocks_PgC,
         totalnonforestcarbonstocks_PgC = 
           OVLcarbonstocks_PgC + rangelandscarbonstocks_PgC + SCcarbonstocks_PgC + croppastureinfrastructurecarbonstocks_PgC, 
         totalforestcarbonstocksperc = totalforestcarbonstocks_PgC  * 100/sum(totalforestcarbonstocks_PgC),
         totalnonforestcarbonstocksperc = totalnonforestcarbonstocks_PgC  * 100/sum(totalnonforestcarbonstocks_PgC),
         OVLcarbonstocksperc = OVLcarbonstocks_PgC * 100/(totalforestcarbonstocks_PgC + totalnonforestcarbonstocks_PgC),
         forestcontribution = totalforestcarbonstocks_PgC/(totalforestcarbonstocks_PgC + totalnonforestcarbonstocks_PgC),
         rangelandcontribution = rangelandscarbonstocks_PgC/(totalforestcarbonstocks_PgC + totalnonforestcarbonstocks_PgC),
         croppastureinfracontribution = croppastureinfrastructurecarbonstocks_PgC/(totalforestcarbonstocks_PgC + totalnonforestcarbonstocks_PgC),
         SCcontribution = SCcarbonstocks_PgC/(totalforestcarbonstocks_PgC + totalnonforestcarbonstocks_PgC))

sum(EcosystemCarbonLandCategoryWithShiftingCultivation$totalforestcarbonstocks_PgC) # 340.06 PgC
sum(EcosystemCarbonLandCategoryWithShiftingCultivation$totalnonforestcarbonstocks_PgC) # 110.15 PgC
sum(EcosystemCarbonLandCategoryWithShiftingCultivation$OVLcarbonstocks_PgC) # 76.58 PgC
sum(EcosystemCarbonLandCategoryWithShiftingCultivation$SCcarbonstocks_PgC) # 10.22 PgC

