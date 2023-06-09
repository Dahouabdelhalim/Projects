# a script to prepare the data.frame for the archipelago regression model
library(tidyverse)
library(phylowood)
library(missForest)

# set seed for reproducibility
set.seed(2505)

# load data
load("output/GIFT_island_characteristica.rda")
load("output/GIFT_island_basic_info.rda")
load(file = "output/GIFT_archipelago_ID.rda")
load(file = "output/GIFT_geology.rda")

# calculate the number of shifts for the island file in phylowood
shifts <- phylowood_isl %>% 
  group_by(geo_archipelago) %>% 
  dplyr::summarize(shifts = sum(trt_minshift_nr)) %>% 
  #remove those archipelagos without recorded shifts
  filter(shifts > 0) %>% 
  arrange(desc(shifts)) %>% 
  # change names to match with Gift data
  mutate(geo_archipelago = recode(geo_archipelago, 
                                  `Galapagos Islands` = "GalÃ¡pagos Islands",
                                  `Mozambique Channel Islands` = "Mozambique Channel Islands: Bassas da India, Europa, Glorioso Islands, Juan de Nova")) %>% #also in the archipelago file
  # remove shifts not clearly linked to a specific archipelago
  filter(geo_archipelago != "Ambiguous") 

# merge input and select relevant data
dat <- geoentities %>% 
  #island names and ID
  dplyr::select(entity_ID, geo_entity) %>% 
  # merge in area, distance to continent, latitudem and archipelago name
  left_join(geoentities_env_misc %>%  dplyr::select(entity_ID, area, dist, latitude, arch_lvl_1), by = "entity_ID") %>% 
  # merge in island age
  left_join(geol %>% dplyr::select(entity_ID, age_Ma, geology_state), by = "entity_ID") %>% 
  # merge in archipelago names
  left_join(arch_gift %>% dplyr::select(Arch_ID, Arch_name), by = c( "arch_lvl_1" = "Arch_ID"))%>% #also in the archipelago file
  #recode some of the names to match the shifts
  mutate(Arch_name = recode(Arch_name, 
                            `Galapagos Islands` = "GalÃ¡pagos Islands",
                            `Mozambique Channel Islands` = "Mozambique Channel Islands: Bassas da India, Europa, Glorioso Islands, Juan de Nova")) %>%
  filter(!Arch_name %in% c("Kerguelen Islands",
                           "Crozet",
                           "Prince Edward Islands", 
                           "Heard and MacDonald")) 
# write to disk for lineage age figure
isl_age <- dat
save(isl_age, file = "output/island_per_archipelago_age.rda")

dat <- dat%>% 
  group_by(Arch_name, arch_lvl_1) %>% 
  #summarize relevant variables per archipelago
  dplyr::summarize(age = max(age_Ma, na.rm = TRUE),
            dist = min(dist, na.rm = TRUE),
            lat = mean(latitude, na.rm = TRUE),
            area = sum(area, na.rm = TRUE)) %>%
  # code archipelagos without known age to NA
  mutate(age = ifelse(is.infinite(age), NA, age)) %>% 
  filter(!is.na(Arch_name))

#Add info on Desaventuradas from wikipedia
dat <- bind_rows(dat, data.frame(Arch_name = "Desventuradas Islands",
                                 arch_lvl_1 = NA,
                                 age = NA,
                                 dist = 850,
                                 lat = mean(-26.343611, -26.31, -26.291667, -26.273611),
                                 area = 5.36))
# add shifts
dat <- dat %>% 
  left_join(shifts, by = c("Arch_name" = "geo_archipelago"))

for_imp <- 
  dat %>% 
  # only include archipelagos with all predictors or at least one evolutionary shift
  filter(!(is.na(age) & is.na(shifts)))

# impute age
imp <- for_imp %>% 
  ungroup() %>% 
  dplyr::select(-1, -2, -7)

colSums(apply(imp, 2, is.na))

imp_out <- missForest(xmis = as.matrix(imp))

imp_out <- bind_cols(for_imp %>%  dplyr::select(1:2, 7), 
                 as_tibble(imp_out$ximp))

dat <- imp_out %>% 
  bind_rows(dat %>% filter(!Arch_name %in% imp_out$Arch_name))%>% 
  mutate(shifts = ifelse(is.na(shifts), 0, shifts))
  

# get oceanic/continental
geo <- geol %>% 
  left_join(geoentities_env_misc %>%  dplyr::select(entity_ID, area, dist, latitude, arch_lvl_1), by = "entity_ID") %>% 
  dplyr::select(arch_lvl_1, geology_state) 

li <- unique(geo$arch_lvl_1)
out <- list()

for(i in 1:length(li)){
  print(i)
  sub <- geo %>%  filter(arch_lvl_1 == li[i])
  
  type <- ifelse(any(sub$geology_state == "fragment")|
             any(sub$geology_state == "shelf")|
             any(sub$geology_state == "atoll/floor/fragment/shelf/volcanic")|
             any(sub$geology_state == "fragment/volcanic")|
             any(sub$geology_state == "atoll/shelf/volcanic")|
             any(sub$geology_state == "shelf/volcanic")|
             any(sub$geology_state == "atoll/floor/fragment")|
             any(sub$geology_state == "fragment/shelf/volcanic"),
             "continental",
             "oceanic")
  
  out[[i]] <- tibble(arch_lvl_1 = unique(sub$arch_lvl_1),
                     type = type
  )
}

arch_geo <- bind_rows(out)

dat <- left_join(dat, arch_geo)

# transform predictors
arch_pred <- dat %>% 
  mutate(log_area = log(1 + area)) %>% 
  mutate(log_dist = log(1 + dist)) %>%
  mutate(log_age = log(1+age)) %>% 
  mutate(abs_lat = abs(lat)) 

# Add zeros
arch_pred <- arch_pred %>% 
  replace_na(list(shifts = 0))

# plot predictors to each other
plo <- arch_pred %>% 
  dplyr::select(Arch_name, shifts, type, lat, contains("log")) %>% 
  pivot_longer(cols = c(lat, contains("log"))) %>% 
  mutate(name = recode(name, 
                       lat = "Absolute latitude",
                       log_age = "Minimum age [log Ma]",
                       log_area = "Total current area [log sqkm]",
                       log_dist = "Minimum distance to closest mainland [log km]",
                       log_herbivores = "Mammal herbivore richness [log]"))

ggplot(data = plo, aes(x= value, y = shifts, colour = type))+
  geom_point()+
  #geom_smooth(method = "glm")+
  facet_wrap(name~., ncol =  2, scales = "free_x")+
  ylab("")+
  scale_color_discrete(name = "Archipelago type")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

ggsave("output/supplementary_figure_archipelago_level_data.jpeg", height = 10, width = 8)

# write to disk
save(arch_pred, file = "output/data_for_archipelago_shift_model.rda")
