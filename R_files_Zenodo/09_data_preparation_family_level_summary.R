library(phylowood)
library(tidyverse)

data(phylowood_sp)

dat <- phylowood_sp %>% filter(geo_insularwoody)

# total number of IW species
sp_tot <- dat %>%
  group_by(tax_family) %>% 
  count() %>% 
  rename(total_species = n)
  
#top 3 archipelago, fraction of species
arch_num <- dat %>% 
  dplyr::select(tax_family,geo_archipelago) %>% 
  separate_rows(contains("geo"), sep = ",") %>% 
  group_by(tax_family,geo_archipelago) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  group_by(tax_family) %>% 
  slice(1:3) %>% 
  left_join(sp_tot, by = "tax_family") %>% 
  mutate(frac = n / total_species * 100)

save(arch_num, file = "output/iw_perarchipelago_family.rda")

# top 3 genera number of species
gen_num <- dat %>% 
  dplyr::select(tax_family,tax_genus) %>% 
  group_by(tax_family,tax_genus) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  group_by(tax_family) %>% 
  slice(1:3) %>% 
  left_join(sp_tot, by = "tax_family") %>% 
  mutate(frac = n / total_species * 100)

save(gen_num, file = "output/iw_pergenus_family.rda")

# proportion of open/forest habitat
for_op <- dat %>% 
  dplyr::select(tax_family, geo_habitat) %>% 
  mutate(geo_habitat_broad = ifelse(geo_habitat == "forest", "forest", "open")) %>% 
  mutate(geo_habitat_broad = ifelse(grepl("forest,", geo_habitat) |
                                      grepl(",forest", geo_habitat), 
                                    "widespread", geo_habitat_broad)) %>% 
  mutate(geo_habitat_broad = ifelse(geo_habitat == "NA", NA, geo_habitat_broad)) %>% 
  group_by(tax_family,geo_habitat_broad) %>% 
  count() %>% 
  group_by(tax_family) %>% 
  arrange(tax_family, desc(n)) %>% 
  left_join(sp_tot, by = "tax_family") %>%  
  mutate(frac = n / total_species * 100)

save(for_op, file = "output/iw_perhabitat_family.rda")
  
# proportion arid/humid/mesic/na
moist <- dat %>% 
  dplyr::select(tax_family, geo_moisture_preference) %>%
  mutate(geo_moisture_preference = ifelse(grepl(",", geo_moisture_preference), 
         "widespread", 
         geo_moisture_preference)) %>% 
  group_by(tax_family, geo_moisture_preference) %>% 
  count()
  
save(moist, file = "output/iw_permoisturepreference_family.rda")

# shifts per family
data(phylowood_isl)

n_shift <- phylowood_isl %>% 
  dplyr::select(tax_genus, trt_minshift_nr) %>% 
  left_join(phylowood_sp %>% dplyr::select(tax_genus, tax_family) %>% distinct(), by = "tax_genus") %>% 
  group_by(tax_family) %>% 
  summarize(shifts = sum (trt_minshift_nr))

save(n_shift, file = "output/iw_perfamily_Shifts.rda")
