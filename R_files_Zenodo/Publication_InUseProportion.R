## A mid-20th century estimate of global vegetation carbon stocks based on inventory statistics
## Authors: Bhan et al., Journal of Land Use Science
## SI Figure 3: In-use proportion

## Here we draw a graph where:
# x-axis: Growing stock-derived CS densities: "N:\\H73700\\data\\!pers\\Manan\\Paper3\\CarbonStockTimeseries\\DataFiles\\1_FRAForestCarbonDensitiesCountriesAndAveraged.R"
# y-axis: In-use forests as % of total

setwd("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles")

library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)

options(scipen = 999)
options(digits = 3)

SubContinentalClassification <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/0_CountriesListSubContinentalClassification.csv")

LandUseEstimatesWithShiftingCultivation <- 
  read.csv('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/5_CountriesListLandUseWithShiftingCultivationWithKastnerExtraForest.csv') %>%
  dplyr::select(country, totalforest_Mha)

sum(LandUseEstimatesWithShiftingCultivation$totalforest_Mha) # 3978 Mha

ForestCarbonDensities <- 
  read.csv('N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_ForestCarbonDensitiesForKnownCountries.csv') %>%
  dplyr::select(country, continent, climatezone, subcontinent, forestcarbondensity_tCha, forestcarbondensityallFAO_tCha, forestcarbondensitytype)
#%>%left_join(LandUseEstimatesWithShiftingCultivation, by = c('country'))

sum(ForestCarbonDensities$totalforest_Mha) # 3567 Mha
# forest % of countries with forest density data = 3567/3978 = 89.6%

LandUseEstimates <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/5_CountriesListLandUseBaselineWithKastnerExtraForest.csv") %>%
  dplyr::select(country, totalforest_Mha)

sum(LandUseEstimates$totalforest_Mha) # 4091 Mha
# forest % of countries with forest density data = 3567/4091 = 87.2% 

ForestsInUse <- 
  readxl::read_xlsx("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/FRA1953.xlsx", 
                                  sheet = 'accessible_use') %>%
  mutate(accessible = as.numeric(accessible),
         inuse = as.numeric(inuse),
         unexploited = as.numeric(unexploited)) %>%
  mutate(inuse_Mha = inuse/1000) %>%
  dplyr::select(country, inuse_Mha)

CountriesActGreaterThanPot <- 
  read.csv("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/DataFiles/1_ForestAreasPotentialActualOfCountriesWithForestCarbonDensitiesActGreaterThanPot.csv") %>%
  dplyr::select(country, forestinuse_Mha, source, inuseproportion)

InUseProportion <-
  ForestCarbonDensities %>%
  left_join(LandUseEstimatesWithShiftingCultivation, by = c('country')) %>%
  left_join(CountriesActGreaterThanPot, by = c('country')) %>%
  left_join(ForestsInUse, by = c('country')) %>%
  mutate(forestinuse_Mha = if_else(country == 'Pakistan', 2.469,
                                   if_else(country == 'Bahamas', 0.202,
                                           if_else(country == 'New Guinea (Austr.)', 0.25, 
                                                   if_else(country == 'Afghanistan', 0.35, forestinuse_Mha))))) %>%
  mutate(finalforestinuse_Mha = forestinuse_Mha,
         finalforestinuse_Mha = if_else(is.na(finalforestinuse_Mha), inuse_Mha, finalforestinuse_Mha),
         finalforestinuse_Mha = if_else(country == 'Bechuanaland', totalforest_Mha, finalforestinuse_Mha)) %>%
  mutate(finalinuseproportion = finalforestinuse_Mha/totalforest_Mha,
         finalinuseproportion = round(finalinuseproportion, 3)) %>%
  dplyr::select(country, continent, climatezone, subcontinent, forestcarbondensityallFAO_tCha, totalforest_Mha, finalforestinuse_Mha, finalinuseproportion) %>%
  #filter(!(country == 'Belgian Congo')) %>%
  mutate(climatezone = factor(climatezone, 
                              levels = c('boreal', 'temperate', 'tropical'),
                              labels = c('Boreal', 'Temperate', 'Tropical')))

(PlotInUseProportion <-
    ggplot(data = InUseProportion, aes(x = forestcarbondensityallFAO_tCha, y = finalinuseproportion)) +
    geom_point(aes(size = finalforestinuse_Mha, colour = climatezone, stroke = 1.5), alpha = 1) +
    scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
    guides(size = 'none') +
    theme_bw() +
    scale_size(range = c(2, 10)) +
    xlab('Forest carbon density (tC/ha)') +
    scale_x_continuous(limits = c(0, 600), breaks = seq(0, 600, 200)) +
    scale_y_continuous(name = 'Proportion of forests in use', limits = c(0, 1), breaks = c(0.5, 1)) +
    theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(face = 'bold', size = 14, color = "black"),
        axis.text.y = element_text(face = 'bold', size = 14, color = "black"),
        axis.title.x = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size = 14, face = 'bold', color = 'black'), legend.position = c(0.8, 0.9),
        legend.background = element_rect(fill = alpha(0.4))))

#ggsave(filename = '6_PlotInUseProportion.jpg', PlotInUseProportion, dpi = 300, width = 5.5, height = 5, units = 'in')

rm(PlotInUseProportion, InUseProportion, CountriesActGreaterThanPot, ForestsInUse, LandUseEstimates, ForestCarbonDensities, 
   LandUseEstimatesWithShiftingCultivation, SubContinentalClassification)
