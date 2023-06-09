# A script to prepare the data.frame for further analyses
library(tidyverse)
library(phylowood)

# load data
load("output/GIFT_island_checklist.rda") # load checklist data
data("phylowood_sp")
data("phylowood_syn")
load("output/GIFT_island_characteristica.rda")
load("output/GIFT_traits.rda")
load(file = "output/GIFT_entity_richness.rda")
load(file = "output/GIFT_archipelago_ID.rda")

# A vector to identify monocot species
monocots <-  c("Arecales", 
               "Acorales", 
               "Alismatale", 
               "Asparagales", 
               "Commelinales", 
               "Dasypogonales", 
               "Dioscoreales",
               "Liliales", 
               "Pandanales", 
               "Petrosaviales", 
               "Poales",
               "Zingiberales")

# Archipelago names
aid <- geoentities_env_misc%>% 
  left_join(arch_gift %>% 
              filter(Arch_level == 1) %>%
              dplyr::select(Arch_ID, Arch_name), 
            by = c("arch_lvl_1" = "Arch_ID")) %>% 
  dplyr::select(entity_ID, 
                Arch_name, 
                arch_lvl_1)

# join checklists and phylowood
synonyms <- phylowood_syn %>% 
  mutate(binomial_accepted = paste(str_split_fixed(accepted_name, n = 3, pattern = " ")[,1],
                                   str_split_fixed(accepted_name, n = 3, pattern = " ")[,2], sep = " ")) %>% 
  left_join(phylowood_sp %>%  dplyr::select(tax_name, geo_islandendemic, geo_insularwoody),
            by = c("binomial_accepted" = "tax_name"))

gift <- checklists %>% 
  dplyr::select(entity_ID, species, family, order, native)

dat <-  gift %>% 
  left_join(synonyms %>% 
              mutate(derived_wood = TRUE) %>%
              dplyr::select(name_binomial, 
                            derived_wood,
                            islandendemic = geo_islandendemic,
                            insularwoody = geo_insularwoody), 
            by = c("species" = "name_binomial")) %>% 
  mutate(derived_wood = ifelse(is.na(derived_wood), FALSE, derived_wood)) %>% 
  left_join(traits, by =  "species")%>% 
  mutate(gf = ifelse(derived_wood, "derived_woody", Growth_form_1)) %>% 
  mutate(gf = ifelse(is.na(gf) & !is.na(Growth_form_2), Growth_form_2, gf)) %>% 
  mutate(gf = ifelse(order %in% monocots, "monocot", gf)) %>% 
  mutate(gf = dplyr::recode(gf, forb = "herb", subshrub = "shrub", other = "NA")) %>% 
  distinct(entity_ID, species, .keep_all = TRUE) %>% 
  filter(native == 1) %>% 
  left_join(aid)

sum(dat$derived_wood)
dat %>% filter(derived_wood) %>% dplyr::select(species) %>% distinct() %>% nrow() #only 1600 matches species out of 7000, of which c 1900 should be on islands!

sp_list <- dat

# write out the merged list
save(sp_list, file = "output/species_list_with_dw_species_identified.rda")

## Total species
spnum <- dat %>%
  group_by(entity_ID)%>% 
  left_join(geoentities_env_misc, by = "entity_ID") %>% 
  dplyr::summarize(GIFT_list_spnum = n())

## non-monocot species
spnum <- dat %>% 
  filter(!order %in% monocots) %>% 
  group_by(entity_ID)%>% 
  dplyr::summarize(GIFT_list_spnum_no_monocots = n()) %>% 
  right_join(spnum, by = "entity_ID")

## herbaceous non-monocot species (as a proxy on the number of species that could evolve woodiness)
spnum <- dat %>% 
  filter(!order %in% monocots) %>% 
  filter(Growth_form_1 == "herb") %>% 
  group_by(entity_ID)%>% 
  dplyr::summarize(GIFT_list_no_monocots_herbacous = n()) %>% 
  right_join(spnum, by = "entity_ID")

## derived woody species 
spnum <- dat %>% 
  group_by(entity_ID)%>%
  filter(derived_wood) %>% 
  dplyr::summarize(GIFT_list_dwspnum = n()) %>% 
  right_join(spnum, by = "entity_ID")

## insular woody species
spnum <- dat %>% 
  group_by(entity_ID)%>%
  filter(insularwoody) %>% 
  dplyr::summarize(GIFT_list_iwspnum = n()) %>% 
  right_join(spnum, by = "entity_ID")

##  number of insular woody endemics per island
spnum <- dat %>% 
  group_by(entity_ID)%>%
  filter(derived_wood) %>% 
  filter(islandendemic) %>% 
  dplyr::summarize(GIFT_list_dwislandendemicspnum = n()) %>% 
  right_join(spnum, by = "entity_ID")

### total endemics
# number of total endemics
rich <- richness %>% dplyr::select(entity_ID, 
                                   GIFT_rich_spnum = total_raw,
                                   GIFT_rich_nativespnum = native_raw,
                                   GIFT_rich_naturalizedspnum = naturalized_raw,
                                   GIFT_rich_endemicspnum = endemic_raw_max) %>% 
  mutate(GIFT_rich_endemcispnum = ifelse(is.na(GIFT_rich_endemicspnum), 0, GIFT_rich_endemicspnum))

GIFT_species_numbers <- spnum %>% 
  left_join(rich)

# write species numbers to disk
save(GIFT_species_numbers, file = "output/GIFT_species_numbers_per_island.rda")
