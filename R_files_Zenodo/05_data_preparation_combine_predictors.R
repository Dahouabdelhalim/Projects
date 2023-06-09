library(tidyverse)
library(fasterize)
library(raster)
library(phylowood)
library(sf)
library(writexl)

# load data
load(file = "output/GIFT_species_numbers_per_island.rda")
load(file = "output/phylowood_species_numbers_per_island.rda")
load(file = "output/species_list_with_dw_species_identified.rda")
load(file = "output/herbivore_numbers_islands.rda")
load(file = "output/GIFT_island_characteristica.rda")
load(file = "output/GIFT_traits.rda")
load(file = "output/GIFT_island_climate.rda")
load(file = "output/GIFT_island_basic_info.rda")
load(file = "output/GIFT_entity_names.rda")
load(file = "output/GIFT_archipelago_ID.rda")
load(file = "output/GIFT_geology.rda")
load(file = "output/GIFT_entity_richness.rda")
polys <- st_read(dsn = "input", layer = "geoentities_simple")

# archipelago assignment for all islands
archipelago_assignment <- geoentities_env_misc%>% 
  left_join(arch_gift %>% 
              filter(Arch_level == 1) %>%
              dplyr::select(Arch_ID, Arch_name),
            by = c("arch_lvl_1" = "Arch_ID")) %>% 
  dplyr::select(entity_ID, Arch_name) %>% 
  left_join(entity_names %>% dplyr::select(entity_ID, geo_entity), by = "entity_ID") %>%
  distinct()

write_xlsx(archipelago_assignment, "output/archipelago_assignment")

# more archipelagos
geoentities_env_misc <- geoentities_env_misc%>% 
  left_join(arch_gift %>% filter(Arch_level == 1) %>% dplyr::select(Arch_ID, Arch_name), by = c("arch_lvl_1" = "Arch_ID"))

# select only mean values for climate
geoentities_env_raster <- geoentities_env_raster %>% 
  dplyr::select(entity_ID, contains("mean"))

# geology, standardize to broad categories
geol <- 
  geol %>% 
  mutate(geology_simple = dplyr::recode(geology_state,
                                        volcanic = "oceanic",
                                        atoll = "oceanic",
                                        floor = "oceanic", 
                                        `floor/volcanic` = "oceanic",
                                        `atoll/floor/volcanic` = "oceanic",
                                        `atoll/volcanic` = "oceanic",
                                        shelf = "continental",
                                        fragment = "continental",
                                        `shelf/volcanic` = "mixed",
                                        `atoll/shelf/volcanic` = "mixed",
                                        `fragment/volcanic` = "mixed",
                                        `atoll/floor/fragment/shelf/volcanic` = "mixed",
                                        `atoll/floor/fragment` = "mixed",
                                        `fragment/shelf/volcanic` = "mixed"))%>% 
  dplyr::select(entity_ID, age_Ma, geology_state, geology_simple)

#number of species per growth form
sps_gf <- sp_list %>% 
  # Rename NAs in growth form
  mutate(gf = ifelse(is.na(gf) | gf == "NA", "GIFT_list_gf_spnum_unknown_gf", gf)) %>% 
  # get the number of species per growth form
  group_by(entity_ID, gf) %>% 
  dplyr::summarize(n_species = n()) %>% 
  mutate(gf = dplyr::recode(gf,  
                            derived_woody = "GIFT_list_gf_spnum_dw", 
                            herb = "GIFT_list_gf_spnum_herb",
                            monocot = "GIFT_list_gf_spnum_monocot", 
                            shrub = "GIFT_list_gf_spnum_shrub", 
                            tree = "GIFT_list_gf_spnum_tree"))%>% 
  tidyr::spread(key = gf, value = n_species) %>% 
  replace_na(list(GIFT_list_gf_spnum_dw = 0,
                  GIFT_list_gf_spnum_herb = 0,
                  GIFT_list_gf_spnum_monocot = 0, 
                  GIFT_list_gf_spnum_shrub = 0,
                  GIFT_list_gf_spnum_tree = 0,
                  spnum_unknown_gf = 0)) %>% 
  mutate(GIFT_list_gf_spnum_ancestral_wood = GIFT_list_gf_spnum_shrub + GIFT_list_gf_spnum_tree)

herb_num <- herb_num %>% 
  mutate(entity_ID = as.numeric(entity_ID))


# join data sets
raw_data <- 
  # GIFT species numbers
  GIFT_species_numbers %>%
  # phylowood species numbers
  full_join(phylowood_species_numbers,  by = c("entity_ID" = "GIFT_entity_ID")) %>% 
  # species numbers by GIF growth form
  left_join(sps_gf) %>% 
  # Island names
  left_join(entity_names %>% dplyr::select(entity_ID, geo_entity, entity_class), by = "entity_ID") %>% # add entites' names
  # geographic island characteristics: area, archipelago, latitude, longitude, distance to continent, LGM area
  left_join(geoentities_env_misc, by = "entity_ID") %>% 
  # climatic variables
  left_join(geoentities_env_raster, by = "entity_ID") %>% 
  # geology and age (add and simplify)
  left_join(geol, by = "entity_ID") %>% 
  # add the number of large mammal herbivores
  left_join(herb_num, by = "entity_ID") %>% 
  #tropic/temperate
  mutate(tropical = ifelse(latitude > -23.44 & latitude < 23.44, "tropical", "temperate"))

# write to disk
save(raw_data, file = "output/GIFT_dataset_for_model_and_plots_not_checked.rda")
