# A script tp obtain mammal herbivore distribution from phylacine and a fossil list compiled from literature

# libraries
#remotes::install_github("crowcanyon/dssatr", force = TRUE)

library(tidyverse)
library(sf)
library(raster)
library(dssatr)
library(readxl)

# general data
load(file = "output/GIFT_archipelago_ID.rda")
load("output/GIFT_island_basic_info.rda")
load("output/GIFT_island_characteristica.rda")

arch <- geoentities %>% 
  dplyr::select(entity_ID, geo_entity) %>% 
  left_join(geoentities_env_misc %>%  dplyr::select(entity_ID, arch_lvl_1)) %>% 
  left_join(arch_gift %>% dplyr::select(Arch_name, Arch_ID), by = c("arch_lvl_1" = "Arch_ID"))

# PHYLACINE species list
# load data
load("output/GIFT_island_checklist.rda") 

geoent <- st_read(dsn = "input", layer = "geoentities_simple") %>% 
  filter(entt_ID %in% unique(checklists$entity_ID)) %>% 
  filter(area > 10)

save(geoent, file = "output/relevant_islands_phyllacine.rda")

load("output/relevant_islands_phyllacine.rda")

traits <- read_csv("input/PHYLACINE_1.2-master/Data/Traits/Trait_data.csv") %>% 
  dplyr::select(Family.1.2, Genus.1.2, Binomial.1.2, Terrestrial, Island.Endemicity, Diet.Plant, Mass.g, IUCN.Status.1.2) %>% 
  filter(Terrestrial == 1) %>% 
  filter(Diet.Plant >= 20) %>% 
  filter(Mass.g > 1000) %>% 
  filter(IUCN.Status.1.2 != "EP")

# get ranges phylacine
lis <-  list.files("input/PHYLACINE_1.2-master/Data/Ranges/Present_natural")
lis <- lis[gsub(".tif", "", lis) %in% traits$Binomial.1.2]

## load the data
ras <- lapply(paste("input/PHYLACINE_1.2-master/Data/Ranges/Present_natural/", lis, sep = ""), "raster")
ras <-  stack(ras)


for(i in 1:dim(ras)[3]){
  
  print(paste(i, dim(ras)[3], sep = " - "))
  
  # select relevant raster
  sub <- ras[[i]]
  
  if(sub@data@max > 0){
    # trim to only values with presence to save computing time
    sub[sub < 1] <-  NA
    sub <- trim(sub)
    
    # disaggregate by factor 20, c. 5x5 km
    sub <- raster::disaggregate(sub, fact = 20)
    
    #get points
    out <-  rasterToPoints(sub) %>% 
      as_tibble
    names(out)[3] <-  "layer"
    out <-  filter(out, layer > 0)
    out <- SpatialPoints(out[,c(1,2)], proj4string = CRS(proj4string(sub))) %>% 
      st_as_sf()
    
    # reproject to lat lon
    out <- out %>% st_transform(crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    
    # PIP
    out <- st_join(out, geoent, join = st_intersects) %>% 
      filter(!is.na(entt_ID))
    
    # Create output
    st_geometry(out) <- NULL
    
    out <- out %>% 
      dplyr::select(entt_ID, ge_ntty) %>% 
      distinct() %>% 
      mutate(species = names(ras[[i]]))
    
    
    if(i == 1){
      write_csv(out, path = "output/mammals_per_island_geoent.csv")
    }else{
      write_csv(out, path = "output/mammals_per_island_geoent.csv", append = TRUE)
    }
  }
}

# create final species list
phyl <- read_csv("output/mammals_per_island_geoent.csv") %>% 
  # merge in traits
  left_join(traits, by = c("species" = "Binomial.1.2")) %>% 
  # remove species that have been filtered out because of the low fraction of plant material 
  # in their diet or because they went extinct in prehistoric times
  filter(!is.na(Diet.Plant)) %>% 
  
  #merge in archipelago
  left_join(arch %>%  dplyr::select(entity_ID, Arch_name), by = c("entt_ID" = "entity_ID")) %>% 
  #add data source
  mutate(type = "present_natural", source = "phylacine1.2") %>% 
  dplyr::rename(entity_ID= entt_ID,
                tax_genus = Genus.1.2,
                tax_family = Family.1.2,
                tax_species = species,
                geo_entity = ge_ntty) %>% 
  dplyr::select(-Terrestrial, -Island.Endemicity)

arch_no_phyl <- phyl %>% 
  dplyr::select(Arch_name, tax_species) %>% 
  distinct() %>% 
  group_by(Arch_name) %>% 
  summarize(phyl_sp = n())

# FOSSIL SPECIES LIST ARCHIPELAGO
arch_no_fos <- read_excel("input/archipelagos_roberto_RR_herbivores_omnivores.xlsx") %>% 
  dplyr::select(Arch_name = geo_archipelago, 
                tax_family = Family, 
                tax_species = `Extinct_herbivores_omnivores_BM>1kg`) %>% 
  dplyr::select(Arch_name, tax_species) %>% 
  distinct() %>% 
  filter(!is.na(tax_species)) %>% 
  group_by(Arch_name) %>% 
  summarize(fossil_sp = n())


arch_no_comb <- read_excel("input/archipelagos_roberto_RR_herbivores_omnivores.xlsx") %>% 
  dplyr::select(Arch_name = geo_archipelago, tax_family = Family, tax_species = `Extinct_herbivores_omnivores_BM>1kg`) %>% 
  bind_rows(phyl %>%  dplyr::select(Arch_name, tax_species, tax_family)) %>% 
  dplyr::select(Arch_name, tax_species) %>% 
  distinct() %>% 
  filter(!is.na(tax_species)) %>% 
  group_by(Arch_name) %>% 
  summarize(combined_sp = n())
  
arch_herb <- arch_no_phyl %>% 
  left_join(arch_no_fos) %>% 
  left_join(arch_no_comb) %>% 
  replace_na(list(phyl_sp = 0, fossil_sp = 0))

save(arch_herb, file = "output/herbivore_numbers_archipelagos.rda")

# species list per island
fossils <- read_excel("input/list_of_islands_roberto_RR_herbivores_omnivores.xlsx")

fossils <- fossils %>% 
  dplyr::select(entity_ID, geo_entity, 
                Arch_name, 
                tax_family = Family,
                tax_species = `Extinct_herbivores_omnivores_BM>1_kg`,
                source = Source
  ) %>% 
  # remove islands without entries
  filter(!is.na(tax_species)) %>% 
  # add genus column
  mutate(tax_genus = str_split_fixed(tax_species, pattern = "_", n = 2)[,1]) %>% 
  mutate(type = "fossil")

# numbers are different between the files because not all islands were in the island file
test <- fossils %>% 
  distinct() %>% 
  group_by(Arch_name) %>% 
  summarize(spnum = n()) %>% 
  left_join(arch_herb) %>% 
  dplyr::select(Arch_name, island_list = spnum, archipelago_list = fossil_sp)

save(fossils, file = "output/fossil_herbivores_islands.rda")

# combine lists and write to disk
herbs <- bind_rows(phyl, fossils)

herbs <-  herbs %>% 
  mutate(ID = 1:nrow(.))

dupls <- herbs %>% 
  distinct(entity_ID, tax_species, .keep_all = TRUE)

herbs <- herbs %>% 
  mutate(marked_duplicate = ifelse(ID %in% dupls$ID, FALSE, TRUE)) %>% 
  arrange(geo_entity, tax_species)

write_csv(herbs,  file = "output/herbivore_list_islands.csv")

# This list will be manually checked by Roberto Rozzi for duplicates and missing entries
