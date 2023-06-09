## A mid-20th century estimate of global vegetation carbon stocks based on inventory statistics
## Authors: Bhan et al., Journal of Land Use Science
## Forest carbon densities

## Here we obtain CS estimates per country by combining:
# 1. Growing stocks-derived CS in FRA 1953, 1958 and 1963
# 2. Independent-derived estimates
# 3. Making act = pot if act > pot
# 4. Obtaining continental averages

## FRA 1953 estimates: "N:\\H73700\\data\\!pers\\Manan\\Paper3\\CarbonStockTimeseries\\DataFiles\\1_GrowingStocktoCarbonStocks.R"
## FRA 1958 and 1963 estimates: "N:\\H73700\\data\\!pers\\Manan\\Paper3\\CarbonStockTimeseries\\DataFiles\\0_CountriesListNonGrowingStocksCountries.R"

setwd("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles")

library(dplyr)
library(magrittr)
library(tidyr)

options(scipen = 999)
options(digits = 3)

SubContinentalClassification <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_CountriesListSubContinentalClassification.csv")

## FRA 1953----

ForestsInUse <- readxl::read_xlsx("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/FRA1953.xlsx", 
                                  sheet = 'accessible_use') %>%
  mutate(accessible = as.numeric(accessible),
         inuse = as.numeric(inuse),
         unexploited = as.numeric(unexploited)) %>%
  mutate(inuse_Mha = inuse/1000) %>%
  dplyr::select(country, inuse_Mha)

FRA1953CarbonStockDensities <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_CountryCarbonStockDensityOfGrowingStockCountries.csv") %>%
  rename(forestcarbondensity_tCha = carbondensity) %>%
  dplyr::select(country, forestcarbondensity_tCha) %>%
  left_join(ForestsInUse, by = c('country')) 
#%>%
#filter(!(country == 'Panama' | country == 'Belgian Congo' | country == 'Bechuanaland'))

FRA1953CarbonStockDensitiesCompleteCountries <- FRA1953CarbonStockDensities %>%
  filter(!(country == 'United States' | country == 'Hawaii' | country == 'Alaska' |
             country == 'Indonesia' | country == 'British North Borneo' | country == 'New Guinea (Neth.)' |
             #country == 'India' | country == 'French India' | country == 'Portuguese India' |
             country == 'Germany, Eastern' | country == 'Germany, Western' | country == 'Germany, Berlin' | country == 'Saar' |
             country == 'French Morocco' | country == 'Spanish Morocco' | country == 'Ifni' |
             country == 'Japan' | country == 'Ryukyu Islands' |
             #country == 'Svalbard' | country == 'Jan Mayen' |
             #country == 'Italy' | country == 'Trieste' |
             country == 'Malaya' | country == 'Sarawak' |
             #country == 'Panama' | country == 'Panama Canal Zone' |
             country == 'Somalia' | country == 'British Somaliland' |
             #country == 'Tanganyika' | country == 'Zanzibar' |
             #country == 'Yemen' | country == 'Aden Protectorate' |
             country == 'Great Britain' | country == 'Northern Ireland'
  ))



## United States, Alaska and Hawaii----

UnitedStatesAlaskaHawaii <- FRA1953CarbonStockDensities %>%
  filter(country == 'United States' | country == 'Alaska' | country == 'Hawaii') %>%
  mutate(country = 'United States') %>%
  group_by(country) %>%
  mutate(product = forestcarbondensity_tCha * inuse_Mha) %>%
  mutate(productsum = sum(product), sum = sum(inuse_Mha)) %>%
  mutate(forestcarbondensity_tCha = productsum/sum, inuse_Mha = sum(inuse_Mha)) %>%
  dplyr::select(-c(product, productsum, sum)) %>%
  distinct(country, .keep_all = TRUE)

## Indonesia, British North Borneo & New Guinea----

IndonesiaBritishNorthBorneoNewGuineaNeth <- FRA1953CarbonStockDensities %>%
  filter(country == 'Indonesia' | country == 'British North Borneo' | country == 'New Guinea (Neth.)') %>%
  mutate(country = 'Indonesia') %>%
  group_by(country) %>%
  mutate(product = forestcarbondensity_tCha * inuse_Mha) %>%
  mutate(productsum = sum(product), sum = sum(inuse_Mha)) %>%
  mutate(forestcarbondensity_tCha = productsum/sum, inuse_Mha = sum(inuse_Mha)) %>%
  dplyr::select(-c(product, productsum, sum)) %>%
  distinct(country, .keep_all = TRUE)

## NA: India, French and Portuguese----




## Germany and other states----

GermanyWesternEasternBerlinSaar <- FRA1953CarbonStockDensities  %>%
  filter(country == 'Germany, Eastern' | country == 'Germany, Western' | country == 'Germany, Berlin' | country == 'Saar') %>%
  mutate(country = 'Germany, Western') %>%
  group_by(country) %>%
  mutate(product = forestcarbondensity_tCha * inuse_Mha) %>%
  mutate(productsum = sum(product), sum = sum(inuse_Mha)) %>%
  mutate(forestcarbondensity_tCha = productsum/sum, inuse_Mha = sum(inuse_Mha)) %>%
  dplyr::select(-c(product, productsum, sum)) %>%
  distinct(country, .keep_all = TRUE)

## Morocco and other states----

MoroccoFrenchSpanishIfni <- FRA1953CarbonStockDensities %>%
  filter(country == 'French Morocco' | country == 'Spanish Morocco' | country == 'Ifni') %>%
  mutate(country = 'French Morocco') %>%
  group_by(country) %>%
  mutate(product = forestcarbondensity_tCha * inuse_Mha) %>%
  mutate(productsum = sum(product), sum = sum(inuse_Mha)) %>%
  mutate(forestcarbondensity_tCha = productsum/sum, inuse_Mha = sum(inuse_Mha)) %>%
  dplyr::select(-c(product, productsum, sum)) %>%
  distinct(country, .keep_all = TRUE)

## Japan and Ryukyu Islands----

JapanRyukyuIslands <- FRA1953CarbonStockDensities %>%
  filter(country == 'Japan' | country == 'Ryukyu Islands') %>%
  mutate(country = 'Japan') %>%
  group_by(country) %>%
  mutate(product = forestcarbondensity_tCha * inuse_Mha) %>%
  mutate(productsum = sum(product), sum = sum(inuse_Mha)) %>%
  mutate(forestcarbondensity_tCha = productsum/sum, inuse_Mha = sum(inuse_Mha)) %>%
  dplyr::select(-c(product, productsum, sum)) %>%
  distinct(country, .keep_all = TRUE)

## NA: Italy and Trieste----


## Malaysia and Sarawak----

MalaysiaSarawak <- FRA1953CarbonStockDensities %>%
  filter(country == 'Malaya' | country == 'Sarawak') %>%
  mutate(country = 'Malaya') %>%
  group_by(country) %>%
  mutate(product = forestcarbondensity_tCha * inuse_Mha) %>%
  mutate(productsum = sum(product), sum = sum(inuse_Mha)) %>%
  mutate(forestcarbondensity_tCha = productsum/sum, inuse_Mha = sum(inuse_Mha)) %>%
  dplyr::select(-c(product, productsum, sum)) %>%
  distinct(country, .keep_all = TRUE)


## NA: Panama and Panama Canal Zone----

## Somalia and British Somaliland----

SomaliaBritishSomaliland <- FRA1953CarbonStockDensities %>%
  filter(country == 'Somalia' | country == 'British Somaliland') %>%
  mutate(country = 'Somalia') %>%
  group_by(country) %>%
  mutate(product = forestcarbondensity_tCha * inuse_Mha) %>%
  mutate(productsum = sum(product), sum = sum(inuse_Mha)) %>%
  mutate(forestcarbondensity_tCha = productsum/sum, inuse_Mha = sum(inuse_Mha)) %>%
  dplyr::select(-c(product, productsum, sum)) %>%
  distinct(country, .keep_all = TRUE)

## United Kingdom----

UnitedKingdom <- FRA1953CarbonStockDensities %>%
  filter(country == 'Great Britain' | country == 'Northern Ireland') %>%
  mutate(country = 'Great Britain') %>%
  group_by(country) %>%
  mutate(product = forestcarbondensity_tCha * inuse_Mha) %>%
  mutate(productsum = sum(product), sum = sum(inuse_Mha)) %>%
  mutate(forestcarbondensity_tCha = productsum/sum, inuse_Mha = sum(inuse_Mha)) %>%
  dplyr::select(-c(product, productsum, sum)) %>%
  distinct(country, .keep_all = TRUE)

## NA: Tanganyika and Zanzibar----

## NA: Yemen and Aden----

## NA: Jan Mayen and Svalbard----

# FRA 1953 densities + forests in use---- 

FRA1953CarbonStockDensitiesInUseForests <- rbind(FRA1953CarbonStockDensitiesCompleteCountries, 
                                                 UnitedStatesAlaskaHawaii, IndonesiaBritishNorthBorneoNewGuineaNeth, 
                                                 GermanyWesternEasternBerlinSaar, MoroccoFrenchSpanishIfni, JapanRyukyuIslands,
                                                 MalaysiaSarawak, SomaliaBritishSomaliland, UnitedKingdom) %>%
  mutate(source = 'FAO 1953') %>%
  mutate_if(is.numeric, round, digits = 3)


# FRA 1953 densities----

FRA1953CarbonStockDensities <- rbind(FRA1953CarbonStockDensitiesCompleteCountries, 
                                     UnitedStatesAlaskaHawaii, IndonesiaBritishNorthBorneoNewGuineaNeth, 
                                     GermanyWesternEasternBerlinSaar, MoroccoFrenchSpanishIfni, JapanRyukyuIslands,
                                     MalaysiaSarawak, SomaliaBritishSomaliland, UnitedKingdom) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  dplyr::select(country, forestcarbondensity_tCha)

# n = 78, earlier 75 but right now not excluded 'Panama' | country == 'Belgian Congo' | country == 'Bechuanaland'

rm(FRA1953CarbonStockDensitiesCompleteCountries, ForestsInUse,
   UnitedStatesAlaskaHawaii, IndonesiaBritishNorthBorneoNewGuineaNeth, 
   GermanyWesternEasternBerlinSaar, MoroccoFrenchSpanishIfni, JapanRyukyuIslands,
   MalaysiaSarawak, SomaliaBritishSomaliland, UnitedKingdom)


## FRA 1958 and 1963 densities + forests in use----

FRA19581963CarbonStockDensitiesInUseForests <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_CountriesListForestCarbonStockDensitiesFromLaterAssessments.csv") %>%
  rename(forestcarbondensity_tCha = carbondensity) %>%
  filter(!(country == 'Zanzibar')) %>%
  rename(inuse_Mha = forestinuse_Mha) %>%
  dplyr::select(country, forestcarbondensity_tCha, inuse_Mha, source)


## FRA 1958 and 1963 densities----

# n = 22
# file: "N:\\H73700\\data\\!pers\\Manan\\Paper3\\CarbonStockTimeseries\\DataFiles\\0_CountriesListNonGrowingStocksCountries.R"
# made a mistake, need to put in values again from FRA 1958 and 1963 if required, but short list of countries now known
# remove Zanzibar because that is covered under Tanzania

FRA19581963CarbonStockDensities <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_CountriesListForestCarbonStockDensitiesFromLaterAssessments.csv") %>%
  rename(forestcarbondensity_tCha = carbondensity) %>%
  dplyr::select(country, forestcarbondensity_tCha) %>%
  filter(!(country == 'Zanzibar'))

## 1953 + 1958 + 1963 densities + forests in use----

FRAForestCarbonDensitiesInUseForests <- 
  rbind(FRA1953CarbonStockDensitiesInUseForests, FRA19581963CarbonStockDensitiesInUseForests)

#write.csv(FRAForestCarbonDensitiesInUseForests, '0_CountriesListFRA195319581963ForestCarbonDensitiesInUseForests.csv')

rm(FRA1953CarbonStockDensitiesInUseForests, FRA19581963CarbonStockDensitiesInUseForests, FRAForestCarbonDensitiesInUseForests)

## 1953 + 1958 + 1963 densities----

FRAForestCarbonDensities <- rbind(FRA1953CarbonStockDensities, FRA19581963CarbonStockDensities) %>%
  left_join(SubContinentalClassification, by = c('country')) # n = 100

## CS DENSITIES FOR ALL KNOWN COUNTRIES----

LandUseEstimate <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/5_CountriesListLandUseBaselineWithKastnerExtraForest.csv") %>%
  dplyr::select(country, landarea_Mha)

PotentialForestCarbonStocks <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/3_PotentialForestCarbonStocks.csv")

IndependentStudyCarbonDensities <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/4_CountryCarbonDensitiesComparison_UsedIndependentPotential.csv") %>%
  dplyr::select(country, countryderived_tCha) %>%
  filter(!(country == 'Canada')) # the Canada estimate is from Kurz 1989 for the 1980s

# Independent vs FRA----

ForestCarbonDensitiesIndependentFRA <- LandUseEstimate %>%
  
  # Join FRA 1953, 1958 and 1963 forest carbon densities
  left_join(FRAForestCarbonDensities, by = c('country')) %>%
  
  # Join potential CS
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  
  # Join independent CS studies
  left_join(IndependentStudyCarbonDensities, by = c('country')) %>%
  filter(!(is.na(countryderived_tCha))) %>%
  dplyr::select(country, forestcarbondensity_tCha, potential_tCha, countryderived_tCha)

#write.csv(ForestCarbonDensitiesIndependentFRA, '1_ForestCarbonDensitiesForIndependentComparedToFRA.csv')

rm(ForestCarbonDensitiesIndependentFRA)

# Final CS densities----

ForestCarbonDensitiesForKnownCountries <- LandUseEstimate %>%
  
  # Join FRA 1953, 1958 and 1963 forest carbon densities
  left_join(FRAForestCarbonDensities, by = c('country')) %>%
  
  # Join potential CS
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  
  # Join independent CS studies
  left_join(IndependentStudyCarbonDensities, by = c('country')) %>%
  filter(!(is.na(forestcarbondensity_tCha))) %>%
  
  mutate(forestcarbondensitytype = 'Independent studies') %>%
  
  # Considering only FRA-based values
  mutate(forestcarbondensityallFAO_tCha = forestcarbondensity_tCha) %>%
  
  # Integrate FAO-based and independent studies
  mutate(forestcarbondensity_tCha = if_else(is.na(countryderived_tCha), 
                                            forestcarbondensity_tCha, 
                                            countryderived_tCha)) %>%
  
  # Calibrate observed values with potential CS
  mutate(forestcarbondensitypotentialadjusted_tCha = if_else(forestcarbondensity_tCha > potential_tCha, 
                                                             potential_tCha,
                                                             forestcarbondensity_tCha)) %>%
  
  # Label these countries
  mutate(forestcarbondensitytype = if_else(is.na(countryderived_tCha), 
                                           'FAO-derived', 'Independent studies')) %>%
  
  dplyr::select(country, landarea_Mha, continent, climatezone, subcontinent, 
                potential_tCha, forestcarbondensity_tCha, forestcarbondensityallFAO_tCha,
                forestcarbondensitypotentialadjusted_tCha, forestcarbondensitytype)

#write.csv(ForestCarbonDensitiesForKnownCountries, '1_ForestCarbonDensitiesForKnownCountries.csv')

## Countries with act > pot----

CountriesWithActGreaterThanPot <- LandUseEstimate %>%
  
  # Join FRA 1953, 1958 and 1963 forest carbon densities
  left_join(FRAForestCarbonDensities, by = c('country')) %>%
  
  # Join potential CS
  left_join(PotentialForestCarbonStocks %>%
              dplyr::select(country, potential_tCha), by = c('country')) %>%
  
  # Join independent CS studies
  left_join(IndependentStudyCarbonDensities, by = c('country')) %>%
  filter(!(is.na(forestcarbondensity_tCha))) %>%
  
  # Integrate FAO-based and independent studies
  mutate(forestcarbondensity_tCha = if_else(is.na(countryderived_tCha), 
                                            forestcarbondensity_tCha, 
                                            countryderived_tCha)) %>%
  
  # Calibrate observed values with potential CS
  mutate(forestcarbondensityfinal_tCha = if_else(is.na(countryderived_tCha),
                                                 forestcarbondensity_tCha,
                                                 countryderived_tCha)) %>%
  
  # Filter for countries with act > pot
  filter(forestcarbondensityfinal_tCha > potential_tCha) %>% # n = 15 
  dplyr::select(country, potential_tCha, forestcarbondensityfinal_tCha)

#write.csv(CountriesWithActGreaterThanPot, '1_ForestCarbonDensitiesActGreaterThanPot.csv')

rm(CountriesWithActGreaterThanPot)

## Subcontinent average----

SubcontinentClimateZoneAveragedCarbonStocks <- ForestCarbonDensitiesForKnownCountries %>%
  group_by(subcontinent, climatezone) %>%
  summarize(avgcarbondensity_tCha = mean(forestcarbondensity_tCha, na.rm = TRUE), n = n()) %>%
  rename(samplecountries = n)

#write.csv(SubcontinentClimateZoneAveragedCarbonStocks, '1_SubcontinentClimateZoneAveragedCarbonStocks.csv')


## CS DENSITIES FOR ALL COUNTRIES----

CountriesListForestCarbonDensities <- LandUseEstimate %>%
  
  # Attach continent classification, known forest CS densities and average CS densities
  left_join(SubContinentalClassification, by = c('country')) %>%
  left_join(ForestCarbonDensitiesForKnownCountries %>%
              dplyr::select(country, 
                            forestcarbondensity_tCha, 
                            forestcarbondensityallFAO_tCha,
                            forestcarbondensitytype), 
            by = c('country')) %>%
  left_join(SubcontinentClimateZoneAveragedCarbonStocks, by = c('subcontinent', 'climatezone')) %>%
  
  # Yemen and Oman given same CS densities as other NAWA countries
  mutate(avgcarbondensity_tCha = if_else(climatezone == 'tropical' & subcontinent == 'Northern Africa & Western Asia',
                                         24.4, avgcarbondensity_tCha)) %>%
  
  # Attaching forest CS densities for all leftover countries
  mutate(forestcarbondensity_tCha = 
           if_else(is.na(forestcarbondensity_tCha), 
                   avgcarbondensity_tCha, 
                   forestcarbondensity_tCha)) %>%
  
  mutate(forestcarbondensityallFAO_tCha = 
           if_else(is.na(forestcarbondensityallFAO_tCha), 
                   avgcarbondensity_tCha, 
                   forestcarbondensityallFAO_tCha)) %>%
  
  mutate(forestcarbondensitytype = if_else(is.na(forestcarbondensitytype), 'Continental average', forestcarbondensitytype)) %>%
  dplyr::select(country, continent, climatezone, subcontinent, 
                forestcarbondensity_tCha, forestcarbondensityallFAO_tCha,
                forestcarbondensitytype)

#write.csv(CountriesListForestCarbonDensities, '1_CountriesListAllForestCarbonDensities.csv') 

rm(CountriesListForestCarbonDensities, SubcontinentClimateZoneAveragedCarbonStocks, 
   ContinentalClassification, ForestCarbonDensitiesForKnownCountries, FRA1953CarbonStockDensities, 
   FRA19581963CarbonStockDensities, FRAForestCarbonDensities, IndependentStudyCarbonDensities,
   LandUseEstimate, PotentialForestCarbonStocks, SubContinentalClassification, ForestCarbonDensitiesIndependentFRA)

## Publication forest carbon stocks----

PublicationCarbonDensities <- 
  read.csv('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_CountriesListAllForestCarbonDensities.csv') %>%
  dplyr::select(country, subcontinent, forestcarbondensity_tCha, forestcarbondensityallFAO_tCha, forestcarbondensitytype)

#write.csv(PublicationCarbonDensities, 'PublicationForestCarbonDensities.csv')
