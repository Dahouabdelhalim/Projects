## to figure out which ant species captured in different trap types and which are in common
## hand samples = ants were observed removing seed(s) from seed caches
## pitfall samples
## subterranean traps - those associated with seeds and control traps (associated
## with beads or passive samples)

## -- load libraries ####
library(here)
library(tidyverse)

## -- set file paths ####
raw_data_path <- here::here("data/raw_data/csv")

## -- load hand sample data ####
# data on ants removing seeds

Ruzi_etal_seed_removal_data <- read_csv(paste(raw_data_path, "Ruzi_etal_seed_removal_data.csv",
                                              sep = "/"), col_types = cols(
                                                .default = col_character(),
                                                diff.between.this.timeperiod.and.next = col_double()))
                        
Ruzi_etal_seed_removal_data

## -- get frequencies of ants removing seeds and number of ant species collected removing seeds ####

hand_samples_raw2 <- Ruzi_etal_seed_removal_data %>%
  # make a new column for plot abbreviation
  mutate(plot_abbre = if_else(Plot == "AVA", "A", Plot)) %>%
  mutate(plot_abbre = if_else(Plot == "25Ha", "25", plot_abbre)) %>%
  mutate(plot_abbre = if_else(Plot == "Pearson", "P", plot_abbre)) %>%
  mutate(plot_abbre = if_else(Plot == "Drayton", "D", plot_abbre)) %>%
  mutate(plot_abbre = if_else(Plot == "Zetek", "Z", plot_abbre)) %>%
  # make a new collumn for the species shorter acronyms
  mutate(spp = if_else(Species == "ApMe", "Ape", Species)) %>%
  mutate(spp = if_else(Species == "CeLo", "Cec", spp),
         spp = if_else(Species == "CoVi", "Coc", spp),
         spp = if_else(Species == "FiIn", "Fic", spp),
         spp = if_else(Species == "GuUl", "Gua", spp),
         spp = if_else(Species == "HiAl", "Hie", spp),
         spp = if_else(Species == "JaCo", "Jac", spp),
         spp = if_else(Species == "LuSe", "Lue", spp),
         spp = if_else(Species == "OcPy", "Och", spp),
         spp = if_else(Species == "ZaEk", "Zan", spp)) %>%
  #  select(Species, spp)
  # need to unite into unique identifiers
  mutate(Season = season) %>%
  unite(plot_season, plot_abbre, season, sep = "") %>%
  unite(ID, spp, plot_season, sep = "_")
hand_samples_raw2

hand_samples_filtered <- hand_samples_raw2 %>%
  # daytime time periods where there were seeds remaining
  filter(Day != "3") %>%
  filter(zero.seeds.remaining == "no") %>%
  filter(Time.Point != "16") %>%
  # had a change in hourly seed count
  filter(diff.between.this.timeperiod.and.next != "0") %>%
  # observed ant present
  filter(Observed.ant.present != "0") %>%
  # ant observed removing seed
  filter(observed.ant.removing != "0")
hand_samples_filtered

names(hand_samples_filtered) # ant species are columns 11 to 29

ant_names <- names(hand_samples_filtered[11:29])
ant_names

total_possible_sampling_timeperiods <- hand_samples_filtered %>%
  count(ID) %>%
  rename(total_possible_samples = n)
total_possible_sampling_timeperiods

## get a tally of number of time periods ants were collected in by ID

hand_samples_removal_tally_all <- hand_samples_filtered %>%
  select(ID, Plot, Species, Season, ant_names) %>%
  pivot_longer(ant_names, names_to = "ant_spp_names", values_to = "present") %>%
  mutate(present = as.double(present)) %>%
  # now to group by plot, season, ant_spp
  group_by(ID, ant_spp_names) %>%
  # to get the number of traps in which the different ant species are present
  summarise(n = sum(present)) %>%
  # now to spread it back out so that the ant sppecies are spread out again
  pivot_wider(names_from = ant_spp_names, values_from = n) %>%
  # add the number of samples from each of the plots
  left_join(total_possible_sampling_timeperiods) %>%
  rename(Sample_num = total_possible_samples)
hand_samples_removal_tally_all 


## need to change this to frequencies now

hand_samples_removal_freq_all <- hand_samples_removal_tally_all %>%
  # to make it so that all the ants are in one column, easier to manipulate this way
  pivot_longer(ant_names, names_to = "ant_spp_names", values_to = "tally") %>%
  # so that it group by the plot_season unique identifier and the ants within those
  # categories
  group_by(ID, ant_spp_names) %>%
  # to get the frequency, take the tally and then divide it by the sample number
  summarise(freq = tally/Sample_num) %>%
  # to spread the ant columns back out which is what I need for the different
  # analyses
  pivot_wider(names_from = ant_spp_names, values_from = freq) %>%
  # adding back in the Sample_num column
  left_join(total_possible_sampling_timeperiods) %>%
  rename(Sample_num = total_possible_samples) %>%
  # making some other columns that are shared with the pitfall data
  mutate(Species = sapply(strsplit(ID, split = "\\\\_"), head, 1),
         Plot_season = sapply(strsplit(ID, split = "\\\\_"), tail, 1),
         Plot_abbre = sapply(strsplit(Plot_season, split = ""), head, 1),
         Season = str_extract(Plot_season, regex("
                                                 (w|d) # w or d
                                                 (e|r) # e or r
                                                 (t|y) # t or y
                                                 ", comments = TRUE))) #%>%
#select(ID, Sample_num, Species, Plot_season, Plot_abbre, Season)
hand_samples_removal_freq_all

## need to remove the columns of columns were ants were not observed removing
## seeds

names(hand_samples_removal_freq_all) # ants are 2:20

hand_ant_spp <- names(hand_samples_removal_freq_all[2:20])

# the ones that have 0's here then are the ones that need to be removed
freq_sums <- colSums(hand_samples_removal_freq_all[hand_ant_spp])

# gives a logical, but it isn't a vector
logicalvect <- freq_sums == "0"

logical_vect_names <- names(logicalvect)
logical_vect <- as.vector(logicalvect)
logical_vect2 <- logical_vect == FALSE # to get the ant spp that are present

# gives the column names that I need to get rid of for belowground wet season samples
colsums_w_0 <- logical_vect_names[logical_vect]

# gives the column names that I want to keep for all belowground wet season samples
colsums_greaterThan_0 <- logical_vect_names[logical_vect2]

## now to remove the columns with those names that have freq's that add to zero
hand_samples_removal_freq_condensed <- hand_samples_removal_freq_all %>%
  select(-colsums_w_0)
hand_samples_removal_freq_condensed

## not to get a condensed list of ant species observed for each tree species
names(hand_samples_removal_freq_condensed)

condensed_hand_ant_spp <- names(hand_samples_removal_freq_condensed[2:14])

treeSpp_antSpp_hand <- hand_samples_removal_freq_condensed %>%
  pivot_longer(condensed_hand_ant_spp, names_to = "ant_spp_names", values_to = "freq") %>%
  ungroup() %>%
  group_by(Species, ant_spp_names) %>%
  summarise(total_freq = sum(freq)) %>%
  filter(total_freq != "0")
treeSpp_antSpp_hand

hand_samples <- treeSpp_antSpp_hand %>%
  mutate(experiment3 = "hand_samples") %>%
  ungroup() %>%
  group_by(experiment3, ant_spp_names) %>%
  summarise(total_freq2 = sum(total_freq)) %>%
  rename(total_freq = total_freq2) %>%
  select(-total_freq)
hand_samples

treeSpp_antSpp_hand %>%
  ungroup() %>%
  count(ant_spp_names)# 13 ant species collected and observed removing seeds



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

## -- to get which ant species are in common between hand samples and pitfall traps ####

names(Ruzi_etal_freq_data_pitfalls_dry_wet_season) # the ant species are 2:60

pitfall_ant_spp <- names(Ruzi_etal_freq_data_pitfalls_dry_wet_season[2:60])

spp_pitfalls <- Ruzi_etal_freq_data_pitfalls_dry_wet_season %>%
  pivot_longer(pitfall_ant_spp, names_to = "ant_spp_names", values_to = "freq") %>%
  mutate(experiment2 = "pitfall") %>%
  group_by(experiment2, ant_spp_names) %>%
  summarise(total_freq = sum(freq)) %>%
  select(-total_freq)
spp_pitfalls

## to get a list of ant species that are in common between pitfalls and hand samples
# by using an inner join
pitfall_hand_in_common <- spp_pitfalls %>%
  inner_join(hand_samples)
pitfall_hand_in_common # 9 ants in common

## using anti_join to get what is unique in x tibble
pitfall_not_hand <- spp_pitfalls %>%
  anti_join(hand_samples)
pitfall_not_hand ## 50 ant species unique to pitfall traps 

hand_not_pitfall <- hand_samples %>%
  anti_join(spp_pitfalls)
hand_not_pitfall # 4 ants are unique to hand samples

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

## -- list of ant species associated with seed caches vs controls ####
names(Ruzi_etal_data_below_frequency_condensed) # ant species are 2:21

below_ant_spp_all <- names(Ruzi_etal_data_below_frequency_condensed[2:21])

## -- separate out the controls from the seed caches
## controls first
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


below_controls_list <- below_controls_condensed %>%
  pivot_longer(colsums_w_samp_b_cont, names_to = "ant_spp_names", values_to = "freq")%>%
  group_by(ant_spp_names) %>%
  summarise(total_freq = sum(freq)) %>%
  mutate(experiment3 = "below_controls") %>%
  ungroup() %>%
  group_by(experiment3, ant_spp_names) %>%
  summarise(total_freq2 = sum(total_freq)) %>%
  rename(total_freq = total_freq2) %>%
  select(-total_freq)
below_controls_list

## seeds

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


below_seeds_list <- below_seeds_condensed %>%
  pivot_longer(colsums_w_samp_b_seeds, names_to = "ant_spp_names", values_to = "freq")%>%
  group_by(ant_spp_names) %>%
  summarise(total_freq = sum(freq)) %>%
  mutate(experiment2 = "below_seeds") %>%
  ungroup() %>%
  group_by(experiment2, ant_spp_names) %>%
  summarise(total_freq2 = sum(total_freq)) %>%
  rename(total_freq = total_freq2) %>%
  select(-total_freq)
below_seeds_list

# determine which ant species are in common and which are unique between seed associated traps and control traps

b_seeds_control_in_common <- below_seeds_list %>%
  inner_join(below_controls_list)
b_seeds_control_in_common # 8 ants in common

b_control_not_seeds <- below_controls_list %>%
  anti_join(below_seeds_list)
b_control_not_seeds ## 2 unique ants in common

b_seeds_not_cont <- below_seeds_list %>%
  anti_join(below_controls_list)
b_seeds_not_cont ## 10 ant species unique to seed associated traps

# determine how many different caches ant species are found in
treeSpp_antSpp_below <- below_seeds_condensed %>%
  pivot_longer(colsums_w_samp_b_seeds, names_to = "ant_spp_names", values_to = "freq") %>%
  ungroup() %>%
  group_by(Species, ant_spp_names) %>%
  summarise(total_freq = sum(freq)) %>%
  filter(total_freq != "0")
treeSpp_antSpp_below

treeSpp_antSpp_below %>%
  ungroup() %>%
  count(ant_spp_names) # to get the number of different tree species caches associated with

treeSpp_antSpp_below %>%
  ungroup() %>%
  count(ant_spp_names) %>%
  filter(n > 1) # 8 ant species associated with more than one tree species cache
