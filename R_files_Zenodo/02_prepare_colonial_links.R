library(countrycode)
library(tidyverse)
library(readxl)


# load data
inp <- read_xls("input/dataverse_files/icowcol.xls")

# wrangel

out <- inp %>% 
  ## convert to iso codes
  mutate(iso2 = countrycode(State, origin = "cown", destination = "iso2c")) %>% 
  mutate(iso3 = countrycode(State, origin = "cown", destination = "iso3c")) %>% 
  ## year of independence
  mutate(indep_year = as.numeric(substr(Indep, start = 1, stop = 4))) %>% 
  ## time since independence
  mutate(timesince = 2021 - indep_year) %>% 
  ## convert code imperialistic power
  mutate(independence_from = countrycode(From, origin = "cown", destination = "iso2c")) %>% 
  # continent
  mutate(cont_ent = countrycode(iso3, origin = "iso3c", destination = "continent")) %>% 
  mutate(cont_suppr = countrycode(independence_from, origin = "iso2c", destination = "continent")) %>% 
  # select relevant columns
  select(iso2, iso3, cow_state = State, Name, indep_year, timesince, independence_from, cont_ent, cont_suppr) %>% 
  ## remove non-country entities
  filter(!is.na(iso2)) %>% 
  ## remove connections not described in the correlates of war dataset
  filter(!is.na(independence_from)) %>% 
  #remove relations on same continent
  filter(cont_ent != cont_suppr) %>% 
  #remove colonies within EUrope
  filter(cont_ent != "Europe") %>% 
  # only since 1800
  filter(indep_year>= 1800)

# check duplicated coutnries and remove
out[out$iso2 %in% out[duplicated(out$iso2),]$iso2,]

out <- out %>% 
  filter(!Name %in% c("Egypt (original)", "Syria (post-UAR)"))

# write to disk
write_csv(out, "output/colonial_ties.csv")
