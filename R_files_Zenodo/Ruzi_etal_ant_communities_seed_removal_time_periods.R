## to determine the number of timeperiods with a change in hourly seed count
## and which ant collections are associated with those timeperiods

## -- load libraries ####
library(here)
library(tidyverse)

## -- set paths ####
raw_data_path <- here::here("data/raw_data/csv")

## -- read in datafile ####

Ruzi_etal_seed_removal_data <- read_csv(paste(raw_data_path, "Ruzi_etal_seed_removal_data.csv",
                                              sep = "/"))

Ruzi_etal_seed_removal_data

names(Ruzi_etal_seed_removal_data)
unique(Ruzi_etal_seed_removal_data$diff.between.this.timeperiod.and.next) ## not sure if this column is needed

raw_data2 <- Ruzi_etal_seed_removal_data %>%
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
  unite(plot_season, plot_abbre, season, sep = "") %>%
  unite(ID, spp, plot_season, sep = "_") %>%
  select(ID, Day, Time.Point, Observed.ant.present, Observed.ant.present,
         sample.of.ant.observed.removing,observed.ant.removing, ant.sample, ant.identification,
         diff.between.this.timeperiod.and.next, zero.seeds.remaining, notes)
raw_data2

## -- i) to get all possible daytime removal time periods ####
# ie. no day 3, no hour 16 to the next morning, no time periods where there were zero seeds remaining

timeperiods <- raw_data2 %>%
  filter(Day != "3") %>%
  filter(zero.seeds.remaining == "no") %>%
  filter(Time.Point != "16")
timeperiods # 2684 rows -- 2684 possible daytime removal time periods

num_daytime_seed_count_change <- timeperiods %>%
  filter(diff.between.this.timeperiod.and.next != "0")
num_daytime_seed_count_change # 158 times -- 158 times there was an hourly change in seed count

158/2684 *100 # percent of time periods where there was an hourly change in seed count
# out of all possible daytime intervals with seeds remaining

## -- ii) to get the number of time periods where ants were present ####

num_daytime_seed_count_change_ants_present <- num_daytime_seed_count_change %>%
  filter(Observed.ant.present != "0")
num_daytime_seed_count_change_ants_present # 75 times ants were present

75/158 *100 # percentage of time periods where there was an hourly change in seed count and ants were present

## -- iii) to get the number of time periods where ants were observed removing seeds ####

num_daytime_seed_count_change_ants_observed_removing <- num_daytime_seed_count_change %>%
  filter(observed.ant.removing != "0")
num_daytime_seed_count_change_ants_observed_removing # 61 timeperiods where ants were observed removing seeds

61/158 * 100 # percentage of time periods where ants were observed removing seeds

## to get the most common ant species observed removing seeds



## -- iv) to get the total possible time periods per plant species ants could have been observed removing seeds ####
# defined as: the total number of samples of ants that could have been collected is based on the total number of time periods
# (both trials pooled) over which there was a change in seed count

total_possible_sampling_timeperiods <- num_daytime_seed_count_change %>%
  count(ID) %>%
  rename(total_possible_samples = n)
total_possible_sampling_timeperiods 
# these values were used as the Sample_number for the datafile used to determine if there were different ant
# communities observed removing seeds of Guazuma ulmifolia, Ochroma pyramidale, Zanthoxlum ekmanii

## -- v) the plant species with the most samples of ants observed removing them ####
amt_w_ant_samples <- num_daytime_seed_count_change_ants_observed_removing %>%
  count(ID) %>%
  rename(tp_ants_removed_seeds = n)
amt_w_ant_samples


# putting the two numbers together so that I can see which ones have samples, and then
## can pull just those with samples out

consolidated_list <- total_possible_sampling_timeperiods %>%
  # this keeps only the values that occur in the the amt_w_ant_samples file
  right_join(amt_w_ant_samples, key = ID) %>%
  # filter by at least 3 ant samples at a given site in a season
  filter(tp_ants_removed_seeds >= "3")
consolidated_list


## need to pull out the ant identifications then
Gua_25_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Gua_25dry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Gua_25_dry_ai

Gua_A_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Gua_Adry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Gua_A_dry_ai

Gua_Z_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Gua_Zdry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Gua_Z_dry_ai

Gua_Z_wet_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Gua_Zwet") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Gua_Z_wet_ai

Och_A_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Och_Adry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Och_A_dry_ai

Och_Z_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Och_Zdry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Och_Z_dry_ai

Och_Z_wet_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Och_Zwet") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Och_Z_wet_ai

Zan_25_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Zan_25dry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Zan_25_dry_ai

Zan_A_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Zan_Adry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Zan_A_dry_ai

Zan_D_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Zan_Ddry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Zan_D_dry_ai

Zan_Z_dry_ai <- num_daytime_seed_count_change_ants_observed_removing %>%
  filter(ID == "Zan_Zdry") %>%
  select(ID, ant.identification, observed.ant.removing, ant.sample, notes)
Zan_Z_dry_ai

## -- vi) list of different ant species that had foragers observed removing seeds ####
# used for Table 3 -- both the number of time period ants were observed removing seeds and the species to put in the hand samples part of the table

## number of timeperiods where ants were observed removing seeds of each tree species
num_daytime_seed_count_change_ants_observed_removing %>%
  separate(ID, into = c("species", "site_season"), sep = "_") %>%
  count(species) # these numbers match when making the list in the following code

# to get the lists
raw_data3 <- Ruzi_etal_seed_removal_data %>%
  filter(Day != "3") %>%
  filter(zero.seeds.remaining == "no") %>%
  filter(Time.Point != "16") %>%
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
  #unite(plot_season, plot_abbre, season, sep = "") %>%
  #unite(ID, spp, plot_season, sep = "_") %>%
  filter(diff.between.this.timeperiod.and.next != "0") %>%
  filter(observed.ant.removing != "0") %>%
  select(spp, plot_abbre, season, Day, Time.Point, Observed.ant.present, Observed.ant.present,
         sample.of.ant.observed.removing,observed.ant.removing, ant.sample, ant.identification,
         diff.between.this.timeperiod.and.next, zero.seeds.remaining, day.time, notes)
raw_data3

unique(raw_data3$spp)

Ape <- raw_data3 %>%
  filter(spp == "Ape") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Ape


Cec <- raw_data3 %>%
  filter(spp == "Cec") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Cec

Gua <- raw_data3 %>%
  filter(spp == "Gua") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Gua

Och <- raw_data3 %>%
  filter(spp == "Och") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Och

TrBr <- raw_data3 %>%
  filter(spp == "TrBr") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
TrBr

Zan <- raw_data3 %>%
  filter(spp == "Zan") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Zan

TrBl <- raw_data3 %>%
  filter(spp == "TrBl") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
TrBl

Fic <- raw_data3 %>%
  filter(spp == "Fic") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Fic


Coc <- raw_data3 %>%
  filter(spp == "Coc") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Coc

Hie <- raw_data3 %>%
  filter(spp == "Hie") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Hie

Jac <- raw_data3 %>%
  filter(spp == "Jac") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Jac

Lue <- raw_data3 %>%
  filter(spp == "Lue") %>%
  select(spp, plot_abbre, season, ant.identification, observed.ant.removing, day.time, notes)
Lue

## -- (vii) to find which ant species had foragers collected in the most number of time periods ####

unique(raw_data3$ant.identification)
raw_data3 %>%
  count(ant.identification) %>%
  arrange(desc(n))

ant_names <- c("PheSus", "EctRui", "SerAma", "PheMul", "WasAur", "CypRim", "ParCor", "MycIst", "MycZet", 
               "PheSim", "PheRug", "LabPra", "Phe025", "PhePug", "AphAra", "Phe091927", "CypCos", "ParBug", "Phe015")

ant_list <- Ruzi_etal_seed_removal_data %>%
  filter(Day != "3") %>%
  filter(zero.seeds.remaining == "no") %>%
  filter(Time.Point != "16") %>%
  filter(diff.between.this.timeperiod.and.next != "0") %>%
  filter(observed.ant.removing != "0") %>%
  # to make all the ant species columns into one column
  pivot_longer(cols = ant_names, names_to = "ant_spp") %>%
  # remove the times when those ant species were not collected
  filter(value != "0") %>%
  # the countes will be the number of time periods in which that foragers of that ant species
  # was collected when there was observed seed removal
  count(ant_spp) %>%
  # arrange in descending order
  arrange(desc(n)) 
ant_list


  


