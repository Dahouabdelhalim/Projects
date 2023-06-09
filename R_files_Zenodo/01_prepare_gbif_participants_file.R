library(tidyverse)
library(readxl)
library(countrycode)

inp <- read_xlsx("input/gbif_members_20210430_wwwgbiforg.xlsx")

out <- inp %>% 
  # remove non-country participants
  filter(status %in% c("Voting participant", "Associate country participant")) %>% 
  #get iso 2 and iso3 codes from country names
  mutate(iso2 = countrycode(name, origin = "country.name", destination = "iso2c")) %>% 
  mutate(iso3 = countrycode(name, origin = "country.name", destination = "iso3c")) %>% 
  select(iso3, iso2, name, status, member_since)

# write to disk for manual inspection
write_csv(out, "output/gbif_participating_countries.csv")
         
