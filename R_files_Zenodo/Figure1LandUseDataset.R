
setwd("N:/H73700/data/!pers/Manan/Paper3/CarbonStockTimeseries/Paper3Graphics")

library(dplyr)
library(magrittr)
library(tidyr)
library(readxl)
library(ggplot2)


## Best-guess----

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

## Without SC----

EcosystemCarbonWithoutShiftingCultivation <- LandUseEstimates %>%
  
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
  #mutate(totalC_shortfallow_tCha = if_else(is.na(totalC_shortfallow_tCha), 0, totalC_shortfallow_tCha),
   #      totalC_longfallow_tCha = if_else(is.na(totalC_longfallow_tCha), 0, totalC_longfallow_tCha)) %>%
  
  #mutate(SCcarbonstocks_MtC = (shortfallow_Mha * totalC_shortfallow_tCha) +
   #        (longfallow_Mha * totalC_longfallow_tCha)) %>%
  
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
           #SCcarbonstocks_MtC +
           annualcropscarbonstocks_MtC + permanentcropscarbonstocks_MtC +
           OVLcarbonstocks_MtC +
           rangelandscarbonstocks_MtC +
           pasturescarbonstocks_MtC + 
           infrastructurecarbonstocks_MtC,
         totalecosystemcarbonstocks_PgC = totalecosystemcarbonstocks_MtC/1000)

sum(EcosystemCarbonWithoutShiftingCultivation$totalecosystemcarbonstocks_PgC) # 460.67 PgC

## Figure 1: Land use dataset----

# NPS
sum(EcosystemCarbonWithShiftingCultivation$nps_Mha) # 1971.3 Mha
sum(EcosystemCarbonWithShiftingCultivation$nps_Mha) # 1971.3 Mha

# infra
sum(EcosystemCarbonWithShiftingCultivation$infrastructure_Mha) # 78.6 Mha
sum(EcosystemCarbonWithShiftingCultivation$infrastructure_Mha) # 78.6 Mha


# Cropland area
sum(EcosystemCarbonWithShiftingCultivation$AnnualCrops_Mha) # 1153.9 Mha
sum(EcosystemCarbonWithShiftingCultivation$PermanentCrops_Mha) # 67.1 Mha

# Grazing area
sum(EcosystemCarbonWithShiftingCultivation$pastures_Mha) # 528.8 Mha
sum(EcosystemCarbonWithShiftingCultivation$rangelands_Mha) # 1762.1 Mha

# Forest
sum(EcosystemCarbonWithoutShiftingCultivation$accessibleforest_Mha) # 1779.6 Mha
sum(EcosystemCarbonWithShiftingCultivation$accessibleforest_Mha) # 1762.3 Mha
# SC in accessible forest = 1779.6 - 1762.3 = 17.3 Mha

sum(EcosystemCarbonWithoutShiftingCultivation$inaccessibleforest_Mha) # 2021.9 Mha
sum(EcosystemCarbonWithShiftingCultivation$inaccessibleforest_Mha) # 1965.4 Mha
# SC in inaccessible forest = 2021.9 - 1965.4 = 56.5 Mha

sum(EcosystemCarbonWithoutShiftingCultivation$extraforest_Mha) # 282.9 Mha
sum(EcosystemCarbonWithShiftingCultivation$extraforest_Mha) # 243.7 Mha
# SC in extra forest = 282.9 - 243.7 = 39.2 Mha

sum(EcosystemCarbonWithoutShiftingCultivation$othervegetatedland_Mha) # 3333.5 Mha
sum(EcosystemCarbonWithShiftingCultivation$othervegetatedland_Mha) # 3193.2 Mha
# SC in inaccessible forest = 3333.5 - 3193.5 = 140.3 Mha

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
