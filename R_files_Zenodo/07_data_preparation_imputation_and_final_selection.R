library(tidyverse)
library(missForest)

# set seed for reproducibility
set.seed(2505)

# load data
load(file = "output/spnum_dataset_checked.rda")

#select relevant predictors for modelling
iwdat <- env_dat %>% 
  dplyr::select(
    geo_entity,
    entity_ID,
    entity_class,
    Arch_name,
    latitude,
    longitude,
    tropical,
    geology_simple,
    phylowood_iwspnum,
    GIFT_list_spnum,
    GIFT_list_iwspnum,
    GIFT_list_spnum_no_monocots,
    area,
    dist,
    mam_total,
    contains("mean"),
    phylowood_dwspnum
  )

# manually impute the missing distance from google earth
iwdat %>% filter(is.na(dist)) %>% View()
iwdat$dist[is.na(iwdat$dist)] <- 1670
iwdat %>% filter(entity_ID == 11587) %>% View()

#geology
iwdat %>% filter(is.na(geology_simple)) %>%  View()

## write the names and entity IDs without 
write_csv(iwdat %>% filter(is.na(geology_simple)) %>%  dplyr::select(geo_entity, entity_ID, Arch_name),
          file = "output/entities_with_missing_geology.csv")

## read manually code file and merge in
add <- read_csv("input/entities_with_missing_geology_added.csv") %>% 
  dplyr::select(entity_ID, geology_simple_add)

iwdat <- iwdat %>%
  left_join(add, by = "entity_ID") %>% 
  mutate(geology_simple = ifelse(is.na(geology_simple), geology_simple_add, geology_simple)) %>% 
  dplyr::select(-geology_simple_add)

# total species numbers
iwdat %>% filter(is.na(GIFT_list_spnum)) %>%  View()
## write the names and entity IDs without 
write_csv(iwdat %>% filter(is.na(GIFT_list_spnum)) %>%  dplyr::select(geo_entity, entity_ID, Arch_name),
          file = "output/entities_with_missing_total_species_richness.csv")

## load manually coded numbers and merge in
add <- read_csv("input/entities_with_missing_total_species_richness_additions.csv") %>% 
  filter(!is.na(estimate)) %>% 
  dplyr::select(entity_ID, estimate)

iw_dat <- iwdat %>%
  left_join(add, by = "entity_ID") %>% 
  mutate(GIFT_list_spnum = ifelse(is.na(GIFT_list_spnum), estimate, GIFT_list_spnum)) %>% 
  dplyr::select(-estimate)

# dataset for plotting
iw_plot <- iw_dat  %>% 
  #remove Australia
  filter(entity_ID != 10542)%>% 
  mutate(geology_simple = recode(geology_simple,
                                 continental = "continental/mixed", 
                                 mixed = "continental/mixed"))

save(iw_plot, file = "output/final_island_dataset_plots.rda")

# dataset for modelling
iw_model <- iw_dat %>%
  #remove Australia
  filter(entity_ID != 10542)%>%
  # rename south island class to island
  mutate(entity_class = ifelse(entity_ID == 12057, "Island", entity_class)) %>% 
  filter(entity_class == "Island") %>% 
  filter(entity_ID != 1058) # remove the philippines

# imputation
imp <- iw_model %>% 
  dplyr::select(-c(1:4), -contains("GIFT"), -contains("phylowood"), -tropical, -geology_simple)

colSums(apply(imp, 2, is.na))

imp_out <- missForest(xmis = as.matrix(imp))

iw_model <- bind_cols(iw_model %>%  dplyr::select(1:4,
                                      contains("GIFT"), 
                                      contains("phylowood"),
                                      tropical,
                                      geology_simple), 
                   as_tibble(imp_out$ximp)) %>% 
 # remove islands without total species number
  filter(!is.na(GIFT_list_spnum)) %>% 
  dplyr::select(-phylowood_dwspnum)

# write to disk to check again
write_csv(iw_model, "output/final_island_dataset_for_checks.csv")
save(iw_model, file = "output/final_island_dataset_models.rda")
