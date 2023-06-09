library(tidyverse)
library(readxl)

# load input data
## This is the manually verified list of mammal herbivores from Roberto Rozzi
herbs <- read_excel("input/herbivore_list_islands_Roberto_revision_included.xlsx")

# write the number of herbivore species to disk for the model
total_n <- herbs %>% 
  dplyr::select(entity_ID, tax_species) %>% 
  distinct() %>% 
  group_by(entity_ID) %>% 
  count()

herb_num <- herbs%>% 
  group_by(entity_ID, type) %>%
  count() %>% 
  pivot_wider(names_from = type, values_from = n) %>% 
  replace_na(list(fossil = 0,
                  present_natural = 0)) %>% 
  left_join(total_n) %>% 
  dplyr::select(entity_ID, mam_fossil = fossil, mam_present_natural = present_natural, mam_total = n)

# write to disk
save(herb_num, file = "output/herbivore_numbers_islands.rda")
save(herbs, file = "output/herbivore_list_islands.rda")
#Write out supplementary file
sup_out <- 
  herbs %>% 
  select(-IUCN.Status.1.2, -Arch_name, -ID, -marked_duplicate) %>% 
  rename(
    GIFT_entity_ID = entity_ID,
    GIFT_entity_name = geo_entity,
    Binomial = tax_species,
    Family = tax_family,
    Genus = tax_genus,
    "Proportion of plant mass in diet [%]" = Diet.Plant,
    "Body mass [g]" = Mass.g,
    "Type of occurrence record" = type,
    "source [see Supplementary material S1 for literature list]" = source
  )


write_xlsx(sup_out, path = "output/09_Dataset_S7_time_integrated_list_of_mammal_herbivores.xlsx")

# list of source for fossil data to get the references
write_csv(herbs %>% select(source) %>% filter(source!= "phylacine1.2") %>%  distinct(), 
          file = "output/references_fossil_occurrences.csv")
