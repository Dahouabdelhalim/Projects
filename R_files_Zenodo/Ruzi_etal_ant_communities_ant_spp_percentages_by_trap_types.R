## to get the percent presence of each of the ant species in different sample types

## -- load libraries ####
library(here)
library(tidyverse)

## -- set file paths ####
raw_data_path <- here::here("data/raw_data/csv")

## -- load pitfall trap data ####
Ruzi_etal_freq_data_pitfalls_dry_wet_season <- read_csv(paste(raw_data_path, "Ruzi_etal_freq_data_pitfalls_dry_wet_season.csv", sep = "/"),
                                                        col_types = cols(
                                                          .default = col_double(),
                                                          ID = col_character(),
                                                          Plot_abbre = col_character(),
                                                          Season = col_character(),
                                                          Species = col_character(),
                                                          Experiment = col_character()
                                                        ))
Ruzi_etal_freq_data_pitfalls_dry_wet_season

## -- most common ants in pitfall traps ####
## To get the most common ant species present in pitfall traps (above-ground),
## based on averaging the frequency of presences at all sites (25, A, D, P, Z) 
## with both seasons pooled (Dry and Wet)

names(Ruzi_etal_freq_data_pitfalls_dry_wet_season) # the species start 2:60
unique(Ruzi_etal_freq_data_pitfalls_dry_wet_season$Sample_num) # all Sample_num = 16

ant_spp <- names(Ruzi_etal_freq_data_pitfalls_dry_wet_season[2:60])
ant_spp

pitfalls_all_samples_pooled_percentages <- Ruzi_etal_freq_data_pitfalls_dry_wet_season %>%
  pivot_longer(ant_spp, names_to = "ant_spp_names",
               values_to = "freq") %>%
  # grouping only by the ant species, that way it will pool all
  # sites and both seasons
  group_by(ID,ant_spp_names) %>%
  summarise(tally = freq*Sample_num) %>%
  mutate(Sample_num = 16) %>%
  ungroup(ID, ant_spp_names) %>%
  group_by(ant_spp_names) %>%
  summarise(n = sum(tally),
            new_Sample_Num = sum(Sample_num)) %>%
  mutate(freq = n/new_Sample_Num,
         percent = freq*100) %>%
  arrange(desc(percent))
pitfalls_all_samples_pooled_percentages  

## -- load in subterranean data file ####
Ruzi_etal_data_below_frequency_condensed <- read_csv(paste(raw_data_path, "Ruzi_etal_data_below_frequency_condensed.csv", sep = "/"),
                                                     col_types = cols(
                                                       .default = col_double(),
                                                       ID = col_character(),
                                                       Plot_abbre = col_character(),
                                                       Season = col_character(),
                                                       Species = col_character(),
                                                       Experiment = col_character()))

Ruzi_etal_data_below_frequency_condensed

## most common ant species in subterranean traps - all sites and all cache types ####

names(Ruzi_etal_data_below_frequency_condensed)  # the species start 2:60
unique(Ruzi_etal_data_below_frequency_condensed$Sample_num) # all Sample_num = 4
## when pooling everything together, there should be 160 samples total
# 5 sites/season x 2 trials/sites x 1 season x 8 cache_types/trial x 2 traps/cache_type = 160

below_ant_spp_all <- names(Ruzi_etal_data_below_frequency_condensed[2:21])
below_ant_spp_all



below_all_samples_pooled_percentages <- Ruzi_etal_data_below_frequency_condensed %>%
  pivot_longer(below_ant_spp_all, names_to = "ant_spp_names",
               values_to = "freq") %>%
  # grouping only by the ant species, that way it will pool all
  # sites and both seasons
  group_by(ID,ant_spp_names) %>%
  summarise(tally = freq*Sample_num) %>%
  ungroup(ID, ant_spp_names) %>%
  group_by(ant_spp_names) %>%
  summarise(n = sum(tally),
            new_Sample_Num = 160) %>%
  mutate(freq = n/new_Sample_Num,
         percent = freq*100) %>%
  arrange(desc(percent))
below_all_samples_pooled_percentages  


## most common ant species in subterranean traps - all seed caches ####

below_seeds <- Ruzi_etal_data_below_frequency_condensed %>%
  filter(Species != "G" & Species != "P")
below_seeds

# the ones that have 0's here then are the ones that need to be removed
freq_sums_b_seeds <- colSums(below_seeds[below_ant_spp_all])

# gives a logical, but it isn't a vector
logicalvect_b_seeds <- freq_sums_b_seeds == "0"

logical_vect_names_b_seeds <- names(logicalvect_b_seeds)
logical_vect_b_seeds <- as.vector(logicalvect_b_seeds)
logical_vect2_b_seeds <- logical_vect_b_seeds == FALSE # to get the ant spp that are present

# gives the column names that I need to get rid of for belowground wet season samples
colsums_w_0_b_seeds <- logical_vect_names_b_seeds[logical_vect_b_seeds]
colsums_w_samp_b_seeds <- logical_vect_names_b_seeds[logical_vect2_b_seeds]

below_seeds_condensed <- below_seeds %>%
  select(-colsums_w_0_b_seeds)
below_seeds_condensed

names(below_seeds_condensed) # 2:19
ant_spp_b <- names(below_seeds_condensed[2:19])

below_seeds_percent <- below_seeds_condensed %>%
  pivot_longer(ant_spp_b, names_to = "ant_spp_names",
               values_to = "freq") %>%
  # grouping only by the ant species, that way it will pool all
  # sites and both seasons
  group_by(ID,ant_spp_names) %>%
  summarise(tally = freq*Sample_num) %>%
  ungroup(ID, ant_spp_names) %>%
  group_by(ant_spp_names) %>%
  summarise(n = sum(tally),
            new_Sample_Num = 120) %>%
  mutate(freq = n/new_Sample_Num,
         percent = freq*100) %>%
  arrange(desc(percent))
below_seeds_percent



## most common ant species in subterranean traps - all control caches ####

below_controls <- Ruzi_etal_data_below_frequency_condensed %>%
  filter(Species == "G" | Species == "P")
below_controls

# the ones that have 0's here then are the ones that need to be removed
freq_sums_b_cont <- colSums(below_controls[below_ant_spp_all])

# gives a logical, but it isn't a vector
logicalvect_b_cont <- freq_sums_b_cont == "0"

logical_vect_names_b_cont <- names(logicalvect_b_cont)
logical_vect_b_cont <- as.vector(logicalvect_b_cont)
logical_vect2_b_cont <- logical_vect_b_cont == FALSE # to get the ant spp that are present

# gives the column names that I need to get rid of for belowground wet season samples
colsums_w_0_b_cont <- logical_vect_names_b_cont[logical_vect_b_cont]

# gives the column names that have samples
colsums_w_samp_b_cont <- logical_vect_names_b_cont[logical_vect2_b_cont]

below_controls_condensed <- below_controls %>%
  select(-colsums_w_0_b_cont)
below_controls_condensed

names(below_controls_condensed) # 2:11
ant_spp_b_cont <- names(below_controls_condensed[2:11])


below_cont_percent <- below_controls_condensed %>%
  pivot_longer(ant_spp_b_cont, names_to = "ant_spp_names",
               values_to = "freq") %>%
  # grouping only by the ant species, that way it will pool all
  # sites and both seasons
  group_by(ID,ant_spp_names) %>%
  summarise(tally = freq*Sample_num) %>%
  ungroup(ID, ant_spp_names) %>%
  group_by(ant_spp_names) %>%
  summarise(n = sum(tally),
            new_Sample_Num = 40) %>%
  mutate(freq = n/new_Sample_Num,
         percent = freq*100) %>%
  arrange(desc(percent))
below_cont_percent
