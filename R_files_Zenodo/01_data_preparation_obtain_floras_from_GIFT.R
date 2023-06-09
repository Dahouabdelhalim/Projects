# A script to download checklists for islands from GIFT and obtain species trait, island characteristics and climate
library(tidyverse)
source("input/GIFT_export_tools_new.R")

#source("database_export_tools.R")

### Path to data on our network drive
wd_path <- ifelse(.Platform$OS.type == "windows", "M:/joint_projects/Checklist DB/","/media/Macroecology/joint_projects/Checklist DB/")

### Define user credentials and connection details. !!! Run this line seperately and not as part of a block of code !!!
# user: wf_idiv
# PW: Dshx4uz1IzIQRDcc
# IP: 134.76.19.22
# Port: 3306
get_credentials()

### Connect to GIFT-DB
conn <- DB_connect()

# Get species lists islands for spermatophyta
# checklists <- DB_get_checklists_conditional(entity_class = c("Island","Island Group","Island Part"), 
#                                             native_indicated = TRUE, natural_indicated = F, end_ref = F, end_list = F, 
#                                             type_ref = 1:11, ref_included = c(1,2,3,4), tax_group = 5, suit_geo = T, 
#                                             exclude_restricted = T, include_names_unique = F, return_query_only = F, 
#                                             complete_taxonomy = TRUE)

checklists <- DB_get_checklists_conditional(entity_class = c("Island"), 
                                            native_indicated = TRUE,
                                            natural_indicated = F, 
                                            end_ref = F, 
                                            end_list = F, 
                                            type_ref = 1:11, 
                                            ref_included = c(1,2,3,4),
                                            tax_group = 5,
                                            suit_geo = T, 
                                            exclude_restricted = T,
                                            include_names_unique = F,
                                            return_query_only = F, 
                                            complete_taxonomy = TRUE)

checklists_broad <- DB_get_checklists_conditional(entity_class = c("Island"), 
                                            native_indicated = TRUE,
                                            natural_indicated = F,
                                            end_ref = F, 
                                            end_list = F, 
                                            type_ref = 1:11,
                                            ref_included = c(1,2,3,4,5,7,8),
                                            tax_group = 5,
                                            suit_geo = FALSE, 
                                            exclude_restricted = T, 
                                            include_names_unique = FALSE, 
                                            return_query_only = F, 
                                            complete_taxonomy = FALSE)

dbGetQuery(conn, "SELECT * FROM reference_included")

# A list of all species for the phylogeny
checklist_raw <- DB_get_checklists_conditional(entity_class = c("Island","Mainland"), 
                                               native_indicated = F,
                                               natural_indicated = F, 
                                               end_ref = F, 
                                               end_list = F, 
                                               type_ref = 1:11, 
                                               ref_included = c(1,2,3,4,5,7,8), 
                                               tax_group = 5, 
                                               suit_geo = F, 
                                               exclude_restricted = T, 
                                               include_names_unique = F,
                                               return_query_only = F,
                                               complete_taxonomy = FALSE)
  
checklists_broad$family <- assign_higher_taxa(work_IDs=checklists_broad$work_ID,
                                              taxon_lvl="family",
                                              higher_lvl=FALSE,
                                              return_ID=FALSE)


# include_names_unique = TRUE adds the species names before our taxonomic standardization to the data set. Maybe you find some of your names in those original names
# suit_geo = FALSE adds lists that are obiously incomplete for the taxon they cover but they can still include valuable information
# complete_taxonomy = FALSE adds lists of islands that we only have data for a taxonomic subset of your group of interest (vascular plants) for. We might for example have a checklist for a cartain family
#                       this also adds data from KEW's World Checklist of Selected Plant families
# I added to ref_included lists that only include endemic, endangered species or other subsets, for example only trees (This adds data from BGCI tree search)


# some examples:
nrow(checklists)
nrow(checklists_broad)
length(which(checklists$ref_ID==10193)) # World Checklist data only for islands that are also covered by lists for all vascular plants
length(which(checklists_broad$ref_ID==10193)) # World Checklist data for all islands in World Checklist
length(which(checklists$ref_ID==10321)) # No data from BGCI tree search
length(which(checklists_broad$ref_ID==10321)) # BGCI tree search data for all islands in BGCI tree search
303 %in% checklists$entity_ID # Cuba not included because we only have a list for angiosperms but nothing for all vascular plants
1066 %in% checklists$entity_ID # Sumatra missing because we only have data from the World Checklist (i.e. only certain families)
303 %in% checklists_broad$entity_ID # Cuba included
1066 %in% checklists_broad$entity_ID # Sumatra included


## add family information to 
checklists$family <- assign_higher_taxa(work_IDs = checklists$work_ID, taxon_lvl="family", 
                                        higher_lvl=FALSE,
                                        return_ID=FALSE)
checklists$order <- assign_higher_taxa(work_IDs=checklists$work_ID, 
                                       taxon_lvl="order", 
                                       higher_lvl=FALSE, 
                                       return_ID=FALSE)

checklists_broad$family <- assign_higher_taxa(work_IDs=checklists_broad$work_ID,
                                              taxon_lvl="family",
                                              higher_lvl=FALSE, 
                                              return_ID=FALSE)
checklists_broad$order <- assign_higher_taxa(work_IDs=checklists_broad$work_ID,
                                             taxon_lvl="order",
                                             higher_lvl=FALSE, 
                                             return_ID=FALSE)


# write species list to disk
save(checklists, file = "output/GIFT_island_checklist.rda")
save(checklists_broad, file = "output/GIFT_island_checklist_including_incomplete.rda")

# get total GIFT species list
checklists2 <- DB_get_checklists_conditional(entity_class = c("Island","Mainland"), 
                                             native_indicated = TRUE, 
                                             natural_indicated = F, 
                                             end_ref = F,
                                             end_list = F, 
                                             type_ref = 1:11,
                                             ref_included = c(1,2,3,4), 
                                             tax_group = 5,
                                             suit_geo = T, 
                                             exclude_restricted = T,
                                             include_names_unique = F,
                                             return_query_only = F, 
                                             complete_taxonomy = TRUE)

checklists2$family <- assign_higher_taxa(work_IDs=checklists2$work_ID,
                                         taxon_lvl="family",
                                         higher_lvl=FALSE, 
                                         return_ID=FALSE)
checklists2$order <- assign_higher_taxa(work_IDs=checklists2$work_ID, 
                                        taxon_lvl="order", 
                                        higher_lvl=FALSE, 
                                        return_ID=FALSE)

gift_tot <- checklists2

save(gift_tot, file = "output/GIFT_total_checklist.rda")

# Species richness
dbGetQuery(conn, "DESCRIBE geoentities_specs")

richness <- dbGetQuery(conn, "SELECT * FROM geoentities_specs WHERE taxon_ID = 5")
save(richness, file = "output/GIFT_entity_richness.rda")

#get island IDS 
entity_names <- dbGetQuery(conn, "SELECT * FROM geoentities")
save(entity_names, file = "output/GIFT_entity_names.rda")

# Island characteristica
env_misc <- dbGetQuery(conn, "SELECT * FROM env_misc")

geoentities_env_misc <- dbGetQuery(conn, "SELECT entity_ID, area, biome, dist, latitude, longitude, 
                                   LGM_area, gmmc, arch_lvl_1, arch_lvl_2, arch_lvl_3, perimeter FROM geoentities_env_misc")

save(geoentities_env_misc, file = "output/GIFT_island_characteristica.rda")


# get archipelago ID, island geology and age and link to continent during ice age
arch_gift <- dbGetQuery(conn, "SELECT * FROM lookup_archipelagos")
save(arch_gift, file = "output/GIFT_archipelago_ID.rda")

# get archipelago ID, island geology and age and link to continent during ice age
dbGetQuery(conn, "DESCRIBE geoentities_geology")
dbGetQuery(conn, "DESCRIBE geology")

geol<- dbGetQuery(conn, "SELECT * FROM geoentities_geology")
geol_lookup <- dbGetQuery(conn, "SELECT * FROM geology")
geol <- geol %>% 
  dplyr::select(entity_ID,
                geology,
                age_Ma, 
                age_original) %>% 
  left_join(geol_lookup,
            by = c("geology" = "ID")) %>% 
  dplyr::select(entity_ID, 
                age_Ma, 
                age_original,
                geology_state = geology.y,
                geology_description = description)

save(geol, file = "output/GIFT_geology.rda")

#get traits per species
traits_meta <- dbGetQuery(conn, "SELECT * FROM traits_meta")
save(traits_meta, file = "output/GIFT_traits_meta.rda")

## get information on growth form
intrest <- c("1.2.1", "1.2.2", "1.6.2", "2.3.1", "2.4.1")

traits <- DB_get_traits(trait_IDs = intrest)

## Get species names for traits
working_names <- dbGetQuery(conn, "SELECT work_ID, species FROM names_work_unique")
traits <- join(traits,working_names, by="work_ID", type="left")

# write to disk
names(traits) <- c("work_id", traits_meta[traits_meta$Lvl3 %in% intrest,]$Trait2, "species")

save(traits, file = "output/GIFT_traits.rda")

# flora completeness
geoentities <- dbGetQuery(conn, "SELECT * FROM  geoentities")
save(geoentities, file = "output/GIFT_island_basic_info.rda")

# Island climate
dbGetQuery(conn, "DESCRIBE env_raster")

env_raster <- dbGetQuery(conn, "SELECT * FROM env_raster")

# wirte a table for the supplementary material, describing th predictors
meta_out <- env_raster %>% 
  filter(layer_name %in% c("wc2.0_bio_30s_01", 
                           "wc2.0_bio_30s_04", 
                           "wc2.0_bio_30s_12",
                           "wc2.0_bio_30s_1",
                           "wc2.0_bio_30s_15",
                           "wc2.0_bio_30s_17",
                           "wc2.0_bio_30s_18",
                           "elev", 
                           "ai_yr",
                           "b1_v0",
                           "b12_v0",
                           "CHELSA_NFD"))

write_csv(meta_out, "output/meta_data_table_climate_gift_raw.csv")


geoentities_env_raster <- DB_get_env_raster(entity_IDs = NULL, 
                                            layers = c("wc2.0_bio_30s_01", 
                                                       "wc2.0_bio_30s_04", 
                                                       "wc2.0_bio_30s_12",
                                                       "wc2.0_bio_30s_1",
                                                       "wc2.0_bio_30s_15",
                                                       "wc2.0_bio_30s_17",
                                                       "wc2.0_bio_30s_18",
                                                       "elev", 
                                                       "ai_yr",
                                                       "b1_v0",
                                                       "b12_v0",
                                                       "CHELSA_NFD"), 
                                            metrics = c("mean","med","sd","n"))

save(geoentities_env_raster, file = "output/GIFT_island_climate.rda")
