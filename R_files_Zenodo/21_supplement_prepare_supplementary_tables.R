library(tidyverse)
library(phylowood)
library(writexl)

# the number of species per genus/ family
data(phylowood_gen)
data(phylowood_sp)

##IW species per genus
sp <- phylowood_sp %>%
  filter(geo_insularwoody) %>% 
  left_join(phylowood_gen %>%  
              select(tax_genus, tax_apgIV), by = "tax_genus") %>% 
  select(tax_name, tax_genus,tax_family, tax_apgIV) %>% 
  group_by(tax_genus) %>%
  summarize(IW_species = n())

## minimum shifts
shifts <- phylowood_isl %>% 
  left_join(phylowood_gen %>%  select(tax_genus, tax_family, tax_apgIV), by = "tax_genus") %>% 
  select(tax_genus, shifts = trt_minshift_nr) %>% 
  group_by(tax_genus) %>% 
  summarize(minimum_transitions = sum(shifts))

refs <- phylowood_gen %>% 
  select(tax_genus, ref_id_free)

out <- left_join(sp, shifts) %>% 
  arrange(desc(minimum_transitions), desc(IW_species)) %>% 
  left_join(phylowood_gen %>% select(tax_genus, tax_family)) %>% 
  select(tax_genus, tax_family, IW_species, minimum_transitions) %>% 
  left_join(refs)

write_xlsx(out, "output/S6_number of insular_woody_species_per_genus.xlsx")

# Number of IW species and shifts per archipelago
## IW species
arch <- phylowood_sp %>% 
  filter(geo_insularwoody) %>% 
  separate_rows(geo_archipelago, sep = ",") %>% 
  mutate(geo_archipelago = trimws(geo_archipelago)) %>% 
  group_by(geo_archipelago) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  rename(number_of_insular_woody_species = n)

## number of shifts
shif <- phylowood_isl %>% 
  group_by(geo_archipelago) %>% 
  summarize(number_of_evolutionary_transitions = sum(trt_minshift_nr)) %>% 
  arrange(desc(number_of_evolutionary_transitions))

out <- full_join(arch, shif)
write_xlsx(out, "output/S7_number_insular_woody_species__and_transitions_per_archipelago.xlsx")  
  
# the final plot/model dataset
load(file = "output/final_island_dataset_plots.rda")
load(file = "output/final_island_dataset_models.rda")

out <- iw_model %>% 
  bind_rows(iw_plot %>% filter(!geo_entity %in% iw_model$geo_entity)) %>% 
  mutate(in_model = entity_ID %in% iw_model$entity_ID) %>% 
  select(-mean_wc2.0_bio_30s_01, -mean_wc2.0_bio_30s_12, -GIFT_list_iwspnum, -mean_wc2.0_bio_30s_17, -phylowood_dwspnum) %>% 
  arrange(desc(phylowood_iwspnum), geo_entity)

out[,10:22] <- round(out[,10:22],3) 

write_xlsx(out, "output/08_Dataset_S6_environment_and_iw_species_per_island.xlsx")
