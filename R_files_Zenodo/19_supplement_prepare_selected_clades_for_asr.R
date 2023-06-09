library(tidyverse)
library(phylowood)
library(readxl)
library(ape)

# load list by Frederic Lens, selected based on the relevance for IW and the quality of the phylogenetic tree
inp <- read_xlsx("input/selected_insular_woody_clades_for_asr.xlsx") %>% 
  separate_rows(genera)

# load tree
tr <- read.tree("input/GBMB_v01_20201027.tre")

# load derived woody trait
dw <- phylowood_sp %>% select(tax_name, geo_onisland, geo_archipelago) %>% 
  mutate(tax_name = gsub(" ", "_", tax_name))

# load gift island checklists
load(file = "output/GIFT_island_checklist_including_incomplete.rda")

gift <- checklists_broad %>% 
  filter(native == 1) %>% 
  select(species) %>% 
  mutate(species = gsub(" " , "_", species)) %>% 
  distinct()

# get species in phylogeny
sp_list <- 
  tr$tip.label %>% 
  as_tibble() %>% 
  mutate(genus = str_split_fixed(value, pattern = "_", n = 2)[,1]) %>% 
  filter(genus %in% inp$genera) %>% 
  left_join(inp %>%  select(genera, family, clade_nr), by = c("genus" = "genera")) %>% 
  left_join(dw, by = c("value" = "tax_name")) %>% 
  mutate(dwood = ifelse(is.na(geo_onisland), "a", "b")) %>% # merge derived woody information
  mutate(onisland = ifelse(geo_onisland, TRUE, FALSE)) %>% 
  mutate(onisland_gift = ifelse(value %in% gift$species, TRUE, FALSE)) %>% 
  mutate(onisland = ifelse(is.na(onisland), onisland_gift, onisland)) %>% 
  dplyr::select(species = value, genus, family, clade_nr, dwood, onisland)

# remove clades without derived woody species in the phylogeny
sel <- sp_list %>% 
  filter(dwood == "b") %>% 
  select(clade_nr) %>% 
  unlist() %>% 
  unique()

sp_list <- sp_list %>% 
  filter(clade_nr %in% sel)

# write to disk
save(sp_list, file = "output/species_list_islands_insularwoody.rda")


