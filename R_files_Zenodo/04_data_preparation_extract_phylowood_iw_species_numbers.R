library(phylowood)
library(tidyverse)
library(writexl)
library(readxl)
library(unheadr)

# load data
data(phylowood_sp)
load(file = "output/GIFT_entity_names.rda")
load(file = "output/GIFT_island_characteristica.rda")
load(file = "output/GIFT_archipelago_ID.rda")

# get list of insular woody specie s in the right format
phy <- phylowood_sp %>% 
  dplyr::select(tax_name,geo_archipelago, geo_onisland, geo_insularwoody, ft_distribution, ft_remarks) %>% 
  filter(geo_onisland) %>% 
  dplyr::select(geo_archipelago, 
                ft_distribution,
                tax_name,
                ft_remarks) %>% 
  arrange(geo_archipelago)

# get GIFT island characteristics
aid <- geoentities_env_misc%>% 
  left_join(arch_gift %>% 
              filter(Arch_level == 1) %>%
              dplyr::select(Arch_ID, Arch_name), 
            by = c("arch_lvl_1" = "Arch_ID")) %>% 
  dplyr::select(entity_ID, 
                Arch_name, 
                arch_lvl_1)

gift <- entity_names %>%
  dplyr::select(entity_ID, geo_entity) %>% 
  left_join(aid, by = "entity_ID") %>% 
  dplyr::select(Arch_name, geo_entity, entity_ID) %>% 
  arrange(Arch_name, geo_entity)
  
# print out for manual checks and corrections
write_csv(phy, "output/phylowood_insular_species_for_GIFT_coding.csv")
write_csv(gift, "output/gift_referenceids.csv")

# reload manually entered data
dat <- read_excel("input/phylowood_insular_species_GIFT_coded.xlsx") %>% 
  dplyr::select(tax_name, GIFT_entity_ID, GIFT_entity_precision)

# replace the IDs on archipelago level to only add those for islands that have at least one recorded dw species wih island precision
isl <- dat %>% 
  filter(GIFT_entity_precision == "island")%>% 
  separate_rows(GIFT_entity_ID, sep = ",")

arch <- dat %>% 
  filter(GIFT_entity_precision == "archipelago") %>% 
  separate_rows(GIFT_entity_ID, sep = ",") %>% 
  filter(GIFT_entity_ID %in% unique(isl$GIFT_entity_ID))

# merge in information if derived woody or insular woody
dat <- bind_rows(arch,isl) %>% 
  left_join(phylowood_sp %>% dplyr::select(tax_name, geo_insularwoody, geo_onisland))


# generate the list of IW species for the supplementary file
sup_out <- phylowood_sp %>% 
    filter(geo_insularwoody) %>% 
  select(tax_name, tax_aut, tax_genus, tax_family) %>% 
  left_join(dat %>% filter(GIFT_entity_precision == "island") %>% 
              select(tax_name, GIFT_entity_ID) %>% 
              unwrap_cols(groupingVar = tax_name, separator = ",")) %>% 
  rename(Binomial = tax_name,
         Authority = tax_aut,
         Genus = tax_genus,
         Family = tax_family)

write_xlsx(sup_out, "output/XX_Dataset_SX_list_of_insular_woody_species.xlsx")

# summarize numbers on island ID level
phylowood_species_numbers <- dat %>% 
  distinct(tax_name, GIFT_entity_ID, .keep_all = TRUE) %>% 
  group_by(GIFT_entity_ID) %>% 
  summarize(phylowood_iwspnum = sum(geo_insularwoody),
            phylowood_dwspnum = sum(geo_onisland)) %>% 
  mutate(GIFT_entity_ID = as.numeric(GIFT_entity_ID))

# write to disk
save(phylowood_species_numbers, file = "output/phylowood_species_numbers_per_island.rda")

# compare numbers between GIFT and phylowood
load(file = "output/GIFT_species_numbers_per_island.rda")

test <- GIFT_species_numbers %>% 
  dplyr::select(entity_ID, GIFT_list_spnum, GIFT_list_iwspnum,GIFT_list_dwspnum) %>% 
  full_join(phylowood_species_numbers, by = c("entity_ID" = "GIFT_entity_ID")) %>% 
  replace_na(list(GIFT_list_iwspnum = 0,
                  GIFT_list_dwspnum = 0,
                  phylowood_iwspnum = 0,
                  phylowood_dwspnum = 0)) %>% 
  mutate(diff_iw = GIFT_list_iwspnum - phylowood_iwspnum) %>% 
  mutate(diff_dw = GIFT_list_dwspnum - phylowood_dwspnum) %>% 
  left_join(entity_names %>% dplyr::select(entity_ID, geo_entity)) %>% 
  dplyr::select(entity_ID, 
         geo_entity, 
         GIFT_total = GIFT_list_spnum,
         GIFT_IW = GIFT_list_iwspnum, 
         phylowood_IW = phylowood_iwspnum,
         difference = diff_iw) %>% 
  arrange(difference)

# plot agaisnt each other
ggplot()+
  geom_point(data = test, aes(GIFT_IW, y = phylowood_IW))

# write table with differences
write_csv(test, file = "output/comparison_GIFT_phylowood.csv")
