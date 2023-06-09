library(tidyverse)
library(raster)
library(elevatr)
library(fasterize)
library(sf)
library(cowplot)

# load data
load(file = "output/GIFT_dataset_for_model_and_plots_not_checked.rda")
polys <- st_read(dsn = "input", layer = "geoentities_simple")
load("output/GIFT_entity_names.rda")
load("output/GIFT_island_basic_info.rda")

# adjust to the minimum size thresholds for the analyses
env_dat <- raw_data %>% 
  filter(area >= 10) %>%
  filter(GIFT_list_spnum >= 20 | phylowood_iwspnum > 0)

# check the list of checklists manually and remove overlapping entities
env_dat %>%
  dplyr::select(entity_ID, geo_entity, Arch_name, GIFT_list_spnum) %>% 
  arrange(Arch_name, geo_entity) %>% 
  View()

rem <- data.frame(ID = c(732,
                         1036,
                         130,
                         1044,
                         215, 
                         1500,
                         291,
                         1033,
                         546,
                         371,
                         984,
                         726, 
                         10544,
                         10640, 	
                         12065, 
                         1315),
                  name = c("Channel Islands","Canary Islands", "Heard and McDonald Islands", "Japan",
                           "New Zealand", "Department of San AndrÃ©s, Providencia y Santa Catalina,", 
                           "Lesser Antilles","Baleares", "Singapore", "Palau", "New Caledonia", "Fiji", 
                           "Southwestern Pacific", "Comoros incl Mayotte", "Ibiza incl. Formentera", "North Aegean Islands"),
                  reason = c("spans multiple entities","spans multiple entities","spans multiple entities","spans multiple entities",
                             "spans multiple entities", "unclear political unit","spans multiple entities","spans multiple entities",
                             "human footprint","spans multiple entities","spans multiple entities","spans multiple entities",
                             "spans multiple entities",
                             "spans multiple entities", "duplicated","spans multiple entities"))
env_dat <- env_dat %>% 
  filter(!entity_ID %in% rem$ID)

# no name
missing_name <- env_dat%>% 
  filter(is.na(geo_entity))

polys %>% filter(entt_ID %in% missing_name$entity_ID)
entity_names %>% filter(entity_ID %in% missing_name$entity_ID)
geoentities %>% filter(entity_ID %in% missing_name$entity_ID)

env_dat <- env_dat %>% 
  mutate(geo_entity = ifelse(entity_ID == 12057, polys %>% filter(entt_ID == 12057) %>% pull(ge_ntty), geo_entity)) %>% 
  mutate(geo_entity = ifelse(entity_ID == 12058, polys %>% filter(entt_ID == 12058) %>% pull(ge_ntty), geo_entity)) %>% 
  mutate(geo_entity = ifelse(entity_ID == 12065, polys %>% filter(entt_ID == 12065) %>% pull(ge_ntty), geo_entity)) %>% 
  mutate(geo_entity = ifelse(entity_ID == 12063, polys %>% filter(entt_ID == 12063) %>% pull(ge_ntty), geo_entity))

env_dat%>% 
  filter(is.na(geo_entity))

# remove duplicated entities
env_dat[env_dat$geo_entity %in% env_dat[duplicated(env_dat$geo_entity),]$geo_entity,] %>% View()

ggplot()+
  geom_sf(data = polys %>% filter(entt_ID == 10824))+
  ggtitle(942)
ggplot()+
  geom_sf(data = polys %>% filter(entt_ID == 366))

# write island list to check
write_xlsx(env_dat, "output/model_data_for_manual_check.xlsx") #looks good for the number of species, number of DW species and entity names

# plot
isls <- polys %>% filter(entt_ID %in% env_dat$entity_ID)

env_dat %>% filter(!entity_ID %in% polys$entt_ID)

ggplot()+
  geom_sf(data = polys, aes(geometry = geometry))+
  geom_sf(data = isls, aes(geometry = geometry), color = "red", fill = "red", alpha = 0.3)

# Australia means minor islands surrounding Australia

# number of NAs per variable
colSums(apply(env_dat, 2, is.na))

# fix NA for number of species
#these values have been assinged NA if there were no species
env_dat <- env_dat %>% 
  replace_na(list(GIFT_list_dwislandendemicspnum  = 0,
                  GIFT_list_iwspnum = 0,
                  GIFT_list_dwspnum = 0,
                  phylowood_iwspnum = 0,
                  phylowood_dwspnum = 0,
                  mam_fossil = 0,
                  mam_present_natural = 0,
                  mam_total = 0))

# remove unnecessary columns with NA
env_dat <- env_dat %>% 
  dplyr::select(-arch_lvl_2, -arch_lvl_3)

# number of NAs per variable
colSums(apply(env_dat, 2, is.na))


#Nas for environmental
# climate
## identify enties with missing climate
sel_clim <- env_dat %>%
  dplyr::select(entity_ID, geo_entity, contains("mean")) %>%
  dplyr::select(-contains("_ai_"), -contains("_elev"))

sel_clim <- sel_clim[colSums(!apply(sel_clim, 1, is.na)) != ncol(sel_clim),]

sel_clim_poly <- polys %>%
  filter(entt_ID %in%sel_clim$entity_ID)

### get new climate varibales from world clim
clim_rep <- lapply(split(sel_clim_poly, f = sel_clim_poly$entt_ID), function(k){getData('worldclim', var='bio', res=0.5, lon=k$point_x, lat = k$point_y)})
clim_rep <- lapply(clim_rep, function(k){out <- k[[c(1,4,12,15,17,18)]]; return(out)})

clim_replace <- list()

for(i in 1:length(sel_clim_poly$entt_ID)){
  print(i)
  sub <- sel_clim_poly$entt_ID[i]

  clim_sub <- clim_rep[names(clim_rep) == sub]
  poly_sub <- polys %>%
    filter(entt_ID == sub) %>%
    fasterize(clim_sub[[1]][[1]])

  clim_out <- stack(mask(clim_sub[[1]], poly_sub))

  clim_replace[[i]] <- data.frame(entity_ID = sub,
                        mean_wc2.0_bio_30s_01 = mean(getValues(clim_out[[1]]), na.rm = TRUE)/ 10,
                        mean_wc2.0_bio_30s_04 = mean(getValues(clim_out[[2]]), na.rm = TRUE)/ 1000,
                        mean_wc2.0_bio_30s_12 = mean(getValues(clim_out[[3]]), na.rm = TRUE),
                        mean_wc2.0_bio_30s_15 = mean(getValues(clim_out[[4]]), na.rm = TRUE),
                        mean_wc2.0_bio_30s_17 = mean(getValues(clim_out[[5]]), na.rm = TRUE),
                        mean_wc2.0_bio_30s_18 = mean(getValues(clim_out[[6]]), na.rm = TRUE))
}

clim_replace <- bind_rows(clim_replace)
clim_replace[is.na(clim_replace)] <- NA

### replace the NAs in the original data.frame
for(i in 1:length(clim_replace$entity_ID)){
  sub <- clim_replace %>% filter(entity_ID == clim_replace$entity_ID[i])

  env_dat <- env_dat %>%
    mutate(mean_wc2.0_bio_30s_01 = ifelse(entity_ID == sub$entity_ID, sub$mean_wc2.0_bio_30s_01, mean_wc2.0_bio_30s_01)) %>%
    mutate(mean_wc2.0_bio_30s_04 = ifelse(entity_ID == sub$entity_ID, sub$mean_wc2.0_bio_30s_04, mean_wc2.0_bio_30s_04)) %>%
    mutate(mean_wc2.0_bio_30s_12 = ifelse(entity_ID == sub$entity_ID, sub$mean_wc2.0_bio_30s_12, mean_wc2.0_bio_30s_12)) %>%
    mutate(mean_wc2.0_bio_30s_15 = ifelse(entity_ID == sub$entity_ID, sub$mean_wc2.0_bio_30s_15, mean_wc2.0_bio_30s_15)) %>%
    mutate(mean_wc2.0_bio_30s_17 = ifelse(entity_ID == sub$entity_ID, sub$mean_wc2.0_bio_30s_17, mean_wc2.0_bio_30s_17)) %>%
    mutate(mean_wc2.0_bio_30s_18 = ifelse(entity_ID == sub$entity_ID, sub$mean_wc2.0_bio_30s_18, mean_wc2.0_bio_30s_18))
}

## elevation
# get the 250m resolution SRTM worldwide and extract coordinates

sel_ele <- env_dat %>% 
  dplyr::select(entity_ID, geo_entity, contains("mean")) %>% 
  dplyr::select(-contains("_ai_"), -contains("_bio"))

sel_ele <- sel_ele[colSums(!apply(sel_ele, 1, is.na)) != ncol(sel_ele),]

sel_ele_poly <- polys %>% 
  filter(entt_ID %in%sel_ele$entity_ID)

# use the elevatr package to get the elevation for those islands that are missing the information from gift. THis will take time 
ele_replace <- list()

for(i in 1:length(sel_ele_poly$entt_ID)){
  print(i)

  sub_ele <- try(get_elev_raster(locations = sel_ele_poly[i,], z = 8, override_size_check = TRUE))
  
  if(class(sub_ele) != "try-error"){
    sub_ele <- sub_ele %>% 
      mask(sel_ele_poly[i, 1])
    
    ele_replace[[i]] <- data.frame(entity_ID = sel_ele_poly[i,]$entt_ID,
                                   mean_elev = mean(getValues(sub_ele), na.rm = TRUE))
  }else{
    test <- st_cast(sel_ele_poly[i,],"POLYGON")
    test_ele <- list()
    
    for(j in 1:nrow(test)){
      print(j)
      sub <- try(get_elev_raster(locations = test[j,], z = 8))
      
      if(class(sub) != "try-error"){ 
          test_ele[[j]] <-  sub %>% 
            mask(test[j, 1])
      }
    }
    # summarize all test ele
    test_ele <- test_ele %>% discard(is.null)
    test_ele <- lapply(test_ele, "getValues")
    ele_replace[[i]] <- data.frame(entity_ID = sel_ele_poly[i,]$entt_ID,
                                   mean_elev = mean(unlist(test_ele), na.rm = TRUE))
  }
}

ele_replace <- bind_rows(ele_replace)
save(ele_replace, file = "output/elevation_replacement.rda")

### replace the NAs in the original data.frame
load(file = "output/elevation_replacement.rda")
for(i in 1:length(ele_replace$entity_ID)){
  sub <- ele_replace %>% filter(entity_ID == ele_replace$entity_ID[i])
  
  env_dat <- env_dat %>% 
    mutate( mean_elev = ifelse(entity_ID == sub$entity_ID, sub$mean_elev, mean_elev))
}

# number of NAs per variable
colSums(apply(env_dat, 2, is.na))

#plot relevant variables to find outliers
## IW GIFT vs total spnum
gif <- ggplot()+
  geom_point(data =  env_dat,
             aes(x = GIFT_list_spnum, y = GIFT_list_iwspnum, color = geology_simple ))+
  geom_label(data =  env_dat %>%  filter(GIFT_list_iwspnum > 15 | GIFT_list_spnum > 1500),
             aes(x = GIFT_list_spnum, y = GIFT_list_iwspnum, color = geology_simple, label = geo_entity ),
             size = 2)+
  theme_bw()+
  theme(legend.position = "bottom")

# IW phylowood against total spnum
phy <- ggplot()+
  geom_point(data =  env_dat,
             aes(x = GIFT_list_spnum, y = phylowood_iwspnum, color = geology_simple ))+
  geom_label(data =  env_dat %>%  filter(phylowood_iwspnum > 15 | GIFT_list_spnum > 1500),
             aes(x = GIFT_list_spnum, y = phylowood_iwspnum, color = geology_simple, label = geo_entity ),
             size = 2)+
  theme_bw()+
  theme(legend.position = "none")

plot_grid(phy,gif, ncol = 1)

# DW phylowood against sp num
hist(env_dat$phylowood_iwspnum)
hist(log(env_dat$phylowood_iwspnum))
table(env_dat$phylowood_iwspnum)

# species number vs size
# remove australia for the modelling!!!
plot(GIFT_list_spnum ~ area, data = env_dat, xlim = c(0, 100000))
text(GIFT_list_spnum ~ area,
     data = env_dat[env_dat$area > 30000 & env_dat$GIFT_list_spnum < 2500,],
     labels = env_dat[env_dat$area > 30000 & env_dat$GIFT_list_spnum < 2500,]$geo_entity)

plot(phylowood_iwspnum ~ area, data = env_dat)
text(phylowood_iwspnum ~ area,
     data = env_dat[env_dat$area > 30000 & env_dat$phylowood_iwspnum < 2500,],
     labels = env_dat[env_dat$area > 30000 & env_dat$phylowood_iwspnum < 2500,]$geo_entity)

# map vs mat
plot(mean_wc2.0_bio_30s_01 ~ mean_wc2.0_bio_30s_12, data = env_dat)
text(mean_wc2.0_bio_30s_01 ~ mean_wc2.0_bio_30s_12, 
     data = env_dat %>% filter(mean_wc2.0_bio_30s_01 > 100),
     labels = env_dat %>% filter(mean_wc2.0_bio_30s_01 > 100) %>%  pull(geo_entity))

plot(mean_wc2.0_bio_30s_04 ~ mean_wc2.0_bio_30s_12, data = env_dat)
text(mean_wc2.0_bio_30s_04 ~ mean_wc2.0_bio_30s_12, 
     data = env_dat %>% filter(mean_wc2.0_bio_30s_04 > 1000),
     labels = env_dat %>% filter(mean_wc2.0_bio_30s_04 > 1000) %>%  pull(geo_entity))

#elevation
hist(env_dat$mean_elev)
env_dat %>% filter(mean_elev <0) %>%  View()
env_dat %>% filter(mean_elev > 1000) %>%  View()

# for now reverse, wait for patricks reply if other fixes are necessary
env_dat <- env_dat %>% 
  mutate(mean_elev = ifelse(entity_ID == env_dat %>% filter(mean_elev <0) %>% pull(entity_ID), mean_elev * -1, mean_elev))

# other climate variables
hist(env_dat$mean_wc2.0_bio_30s_04)
hist(env_dat$mean_wc2.0_bio_30s_15)
hist(env_dat$mean_wc2.0_bio_30s_17)
hist(env_dat$mean_wc2.0_bio_30s_18)
hist(env_dat$mean_CHELSA_NFD)
hist(env_dat$mean_ai_yr)
hist(env_dat$latitude)
hist(env_dat$dist)
hist(env_dat$perimeter)
hist(env_dat$mean_b12_v0)
hist(env_dat$mean_b1_v0)

env_dat %>% filter(mean_wc2.0_bio_30s_04 > 10) %>%  View()

#herbivore diversity vs island size and dist
plot(mam_fossil ~ mam_total, data = env_dat)
hist(env_dat$mam_total)

plot(mam_total ~ area, data = env_dat, xlim = c(0, 100000))
text(mam_total ~ area, 
     data = env_dat[env_dat$area > 30000 | env_dat$mam_total > 30,],
     labels = env_dat[env_dat$area > 30000 | env_dat$mam_total > 30,]$geo_entity)

plot(mam_total ~ dist, data = env_dat)
text(mam_total ~ dist, 
     data = env_dat[env_dat$dist > 1000 | env_dat$mam_total > 30,],
     labels = env_dat[env_dat$dist > 1000 | env_dat$mam_total > 30,]$geo_entity)

# others
hist(env_dat$dist)
env_dat %>% filter(dist > 3000) %>% dplyr::select(entity_ID, geo_entity)

hist(env_dat$age_Ma)
env_dat %>% filter(age_Ma > 90) %>% dplyr::select(entity_ID, geo_entity)

# write to disk
save(env_dat, file = "output/spnum_dataset_checked.rda")

