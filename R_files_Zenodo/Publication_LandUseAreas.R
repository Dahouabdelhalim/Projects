## A mid-20th century estimate of global vegetation carbon stocks based on inventory statistics
## Authors: Bhan et al., Journal of Land Use Science
## Land use dataset 1950

# Here, I start with the file '1_CountriesListAllLandUse_FAO1953HANPPTimeseries.csv'
# Available at: "N:\\H73700\\data\\!pers\\Manan\\Paper3\\CarbonStockTimeseries\\DataFiles\\1_CountriesListAllLandUse_FAO1953HANPPTimeseriesData.R"

# Then one by one I make distinctions between different land uses, add OVL as a category
# ESTIMATE 1: WITHOUT SHIFTING CULTIVATION
# ESTIMATE 2: WITH SHIFTING CULTIVATION from '0_ShiftingCultivationAdjustment.R'

setwd("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles")

library(dplyr)
library(magrittr)
library(tidyr)

options(scipen = 999)
options(digits = 3)

CountriesListLandUse <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_CountriesListAllLandUseFAO1953HANPPTimeseries.csv")

SubContinentalClassification <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_CountriesListSubContinentalClassification.csv")

## OVL ADDITION----

# DUAL APPROACH: Land use correction by reducing rangeland/inaccessible forest areas, depending on country----

CountriesListWithOVL <- CountriesListLandUse %>%
  
  # Correcting for anomalous forest data (more than country land area)
  # Logic is that Bechuanaland already has no pasture, a country cannot have no pasture and no rangeland
  mutate(accessibleforest_Mha = if_else(country == 'Bechuanaland', # anomalous figure of forests, more than total land area itself
                                        landarea_Mha - (infrastructure_Mha + cropland_Mha + nps_Mha + rangelands_Mha + pastures_Mha),
                                        accessibleforest_Mha),
         totalforest_Mha = if_else(country == 'Bechuanaland',
                                   accessibleforest_Mha, totalforest_Mha)) %>%
  
  # Adding OVL
  mutate(othervegetatedland_Mha = round(landarea_Mha - 
                                          (infrastructure_Mha + cropland_Mha + nps_Mha + rangelands_Mha + pastures_Mha + totalforest_Mha), 3)) %>%
  
  # Correcting for negative OVL areas depending on country situation
  # First with rangelands
  mutate(rangelands_Mha = if_else(country == 'French Equat. Africa' |
                                    country == 'Anglo-Egyptian Sudan' |
                                    country == 'Eritrea' |
                                    country == 'Iran' |
                                    country == 'Mongol. Peopl. Rep.' |
                                    country == 'French Somaliland', rangelands_Mha + othervegetatedland_Mha, rangelands_Mha),
         othervegetatedland_Mha = if_else(country == 'French Equat. Africa' |
                                            country == 'Anglo-Egyptian Sudan' |
                                            country == 'Eritrea' |
                                            country == 'Iran' |
                                            country == 'Mongol. Peopl. Rep.' |
                                            country == 'French Somaliland', 0, othervegetatedland_Mha)) %>%
  
  # Then with inaccessible forests
  mutate(inaccessibleforest_Mha = if_else(country == 'Dominican Republic' |
                                            country == 'French Guiana' |
                                            country == 'Liberia', 
                                          inaccessibleforest_Mha + othervegetatedland_Mha, inaccessibleforest_Mha),
         othervegetatedland_Mha = if_else(country == 'Dominican Republic' |
                                            country == 'French Guiana' |
                                            country == 'Liberia', 
                                          0, othervegetatedland_Mha),
         totalforest_Mha = if_else(country == 'Dominican Republic' |
                                     country == 'French Guiana' |
                                     country == 'Liberia',
                                   accessibleforest_Mha + inaccessibleforest_Mha, totalforest_Mha)) %>%
  
  # With pastures
  mutate(pastures_Mha = if_else(country == 'Germany, Western', 
                                pastures_Mha + othervegetatedland_Mha, pastures_Mha),
         othervegetatedland_Mha = if_else(country == 'Germany, Western', 
                                          0, othervegetatedland_Mha)) %>%
  
  # Total land area
  mutate(landarea_Mha = infrastructure_Mha + cropland_Mha + nps_Mha + rangelands_Mha + pastures_Mha + totalforest_Mha + othervegetatedland_Mha)

rm(CountriesListLandUse)

## PERMANENT CROP FRACTION----

PermanentCropFraction <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_CountriesListFAOSTAT19611965PermanentCropFraction.csv")

CountriesListWithPermanentCropFraction  <- CountriesListWithOVL  %>%
  
  # Attach permanent crop fraction
  left_join(PermanentCropFraction, by = c('country')) %>%
  mutate(MeanPermanentFraction = as.numeric(MeanPermanentFraction)) %>%
  
  # Permanent and annual crops
  mutate(PermanentCrops_Mha = cropland_Mha * MeanPermanentFraction,
         AnnualCrops_Mha = cropland_Mha - PermanentCrops_Mha)

rm(CountriesListWithOVL, PermanentCropFraction)

## EXTRA FOREST AREA SUPPLEMENT----

# NEW, 22.11.2021: We use Kastner 2021 to define Extra Forests for 26 countries, 257 Mha
# See sheet: "N:\\H73700\\data\\!pers\\Manan\\Paper3\\CarbonStockTimeseries\\DataFiles\\1_CountriesListExtraForestCorrectionKastnerForestryWilderness.R"

ForestAreaSupplementCountriesKastner <-
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_CountriesListExtraForestCorrectionKastnerForestryWilderness.csv") %>%
  dplyr::select(country, extraforest_Mha)

# We identified 26 tropical countries with higher forest areas in 2000 than in 1950, but without evidence of any substantial and continued forest area gain
# We compared FAO estimates with (forestry + forested wilderness) from Kastner 2021
# Extra forest = Kastner 2021 - FAO 1953 Forest 1953
# This additional forest area added in 26 countries has to be subtracted from some other land category
# 26 countries, 297 Mha

CountriesListWithForestAreaAdjustment <- CountriesListWithPermanentCropFraction %>%
  left_join(ForestAreaSupplementCountriesKastner, by = c('country')) %>%
  mutate(extraforest_Mha = if_else(is.na(extraforest_Mha), 0, extraforest_Mha))

rm(CountriesListWithPermanentCropFraction)

# We subtract OVL and allocate them to extra forests

CountriesListWithForestAreaAdjustmentInOVL <- CountriesListWithForestAreaAdjustment %>%
  mutate(othervegetatedland_Mha = if_else(othervegetatedland_Mha > 0, 
                                          othervegetatedland_Mha - extraforest_Mha,
                                          othervegetatedland_Mha)) %>%
  mutate(extraforest_Mha = if_else(country == 'Malaya', 10.6, extraforest_Mha),
         othervegetatedland_Mha = if_else(country == 'Malaya', 0, othervegetatedland_Mha)) %>%
  mutate(totalforest_Mha = accessibleforest_Mha + inaccessibleforest_Mha + extraforest_Mha)

#write.csv(CountriesListWithForestAreaAdjustmentInOVL, '5_CountriesListLandUseBaselineWithKastnerExtraForest.csv')

rm(ForestAreaSupplementCountriesKastner, CountriesListWithForestAreaAdjustment)


# SHIFTING CULTIVATION ADJUSTMENT----

ShiftingCultivationAdjustment <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_CountriesListShiftingCultivationAdjustmentFearnsideAdjustment.csv") %>%
  dplyr::select(country, shortfallow_Mha, longfallow_Mha)

# SC database has 44 countries, by disaggregating Malaysia and Indonesia into its constituent units. Need to join that back
MalaysiaIndonesiaAdjustment <- ShiftingCultivationAdjustment %>%
  left_join(CountriesListWithForestAreaAdjustmentInOVL, by = c('country'))

ShiftingCultivationAdjustmentFinal <- ShiftingCultivationAdjustment %>%
  mutate(shortfallow_Mha = if_else(country == 'Malaya', 0.130 + 0.121, shortfallow_Mha),
         shortfallow_Mha = if_else(country == 'Indonesia', 1.326 + 0.067 + 0.615, shortfallow_Mha),
         longfallow_Mha = if_else(country == 'Malaya', 0.218 + 0.201, longfallow_Mha),
         longfallow_Mha = if_else(country == 'Indonesia', 2.21 + 0.113 + 0.615, longfallow_Mha))

rm(ShiftingCultivationAdjustment, MalaysiaIndonesiaAdjustment)


CountriesListWithShiftingCultivation <- CountriesListWithForestAreaAdjustmentInOVL %>%
  
  # Attach shifting cultivation adjustment
  left_join(ShiftingCultivationAdjustmentFinal, by = c('country')) %>%
  mutate(shortfallow_Mha = if_else(is.na(shortfallow_Mha), 0, shortfallow_Mha),
         longfallow_Mha = if_else(is.na(longfallow_Mha), 0, longfallow_Mha))

SCNot <- CountriesListWithShiftingCultivation %>%
  filter(shortfallow_Mha == 0 | country == 'Sierra Leone') # SL has irregaular statistics with SC and forest areas

# Adjustment in extra forest

SCInExtraForest <- CountriesListWithShiftingCultivation %>% 
  filter(shortfallow_Mha != 0 & 
           extraforest_Mha > (shortfallow_Mha + longfallow_Mha)) %>%
  mutate(extraforest_Mha = extraforest_Mha - (shortfallow_Mha + longfallow_Mha))

# Adjustment in OVL

SCDistributedInOVL <- CountriesListWithShiftingCultivation %>%
  filter(shortfallow_Mha != 0 & 
           extraforest_Mha < (shortfallow_Mha + longfallow_Mha) &
           othervegetatedland_Mha > (shortfallow_Mha + longfallow_Mha)) %>%
  mutate(othervegetatedland_Mha = othervegetatedland_Mha - (shortfallow_Mha + longfallow_Mha))

# Adjustment in inaccessible forest

SCDistributedInInaccessible <- CountriesListWithShiftingCultivation %>% 
  filter(shortfallow_Mha != 0 & 
           extraforest_Mha < (shortfallow_Mha + longfallow_Mha) &
           othervegetatedland_Mha < (shortfallow_Mha + longfallow_Mha) &
           inaccessibleforest_Mha > (shortfallow_Mha + longfallow_Mha)) %>%
  mutate(inaccessibleforest_Mha = inaccessibleforest_Mha - (shortfallow_Mha + longfallow_Mha))

# Adjustment in accessible forests

SCDistributedInAccessible <- CountriesListWithShiftingCultivation %>% 
  filter(shortfallow_Mha != 0 & 
           extraforest_Mha < (shortfallow_Mha + longfallow_Mha) &
           othervegetatedland_Mha < (shortfallow_Mha + longfallow_Mha) &
           inaccessibleforest_Mha < (shortfallow_Mha + longfallow_Mha)) %>%
  filter(!(country == 'Sierra Leone')) %>%
  mutate(accessibleforest_Mha = accessibleforest_Mha - (shortfallow_Mha + longfallow_Mha))


CountriesListWithShiftingCultivation <- 
  rbind(SCNot, SCInExtraForest, SCDistributedInOVL, SCDistributedInInaccessible, SCDistributedInAccessible) %>%
  mutate(totalforest_Mha = accessibleforest_Mha + inaccessibleforest_Mha + extraforest_Mha)

#write.csv(CountriesListWithShiftingCultivation, '5_CountriesListLandUseWithShiftingCultivationWithKastnerExtraForest.csv')

rm(CountriesListWithShiftingCultivation, CountriesListWithForestAreaAdjustmentInOVL, ShiftingCultivationAdjustment,
   ShiftingCultivationAdjustmentFinal, SubContinentalClassification,
   SCNot, SCInExtraForest, SCDistributedInOVL, SCDistributedInInaccessible, SCDistributedInAccessible)

## Data for publication----

PublicationLandUse <-
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/5_CountriesListLandUseBaselineWithKastnerExtraForest.csv") %>%
  dplyr::select(3:17)

#write.csv(PublicationLandUse, 'PublicationLandUseDataset.csv')
