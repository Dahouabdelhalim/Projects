library(tidyverse)
library(countrycode)
library(readxl)
library(states)

# load input data
old <- read_csv("input/vdem_variables.csv")
conf <- read_excel("input/ucdp-prio-acd-201.xlsx")
gbif <- read_excel("input/gbif_members_20210430_wwwgbiforg.xlsx")


# conflict
conf <- conf %>%  
  select(location, gwno_loc, year, intensity_level) %>% 
  separate_rows(gwno_loc, sep = ", ") %>% 
  mutate(country = countrycode(gwno_loc, origin = "gwn", destination = "iso3c")) %>% 
  mutate(conflict_hi = ifelse(intensity_level == 2, 1, 0)) %>% 
  mutate(conflict_li = ifelse(intensity_level == 1, 1, 0)) %>% 
  filter(year >= 1960) %>% 
  #add fe that were not matched manually
  mutate(country = ifelse(gwno_loc == 678, "YMN", country)) %>% 
  mutate(country = ifelse(gwno_loc == 680, "YMN", country)) %>% 
  mutate(country = ifelse(gwno_loc == 816, "VNM", country)) %>% 
  mutate(country = ifelse(gwno_loc == 55, "GRN", country)) %>% 
  mutate(country = ifelse(gwno_loc == 345, "SRB", country)) %>% 
  select(country, year, conflict_hi, conflict_li) %>% 
  distinct()%>% 
  mutate(year = as.character(year))

#remove minor conflicts if there is a major conflict already
hi <- conf %>% select(country,year, conflict_hi) %>% 
  filter(conflict_hi == 1)
li <- conf %>% select(country,year, conflict_li) %>% 
  filter(conflict_li == 1)

conf <- full_join(hi,li, by = c("country","year")) %>% 
  replace_na(list(conflict_hi = 0,
                  conflict_li = 0)) %>% 
  mutate(conflict_li = ifelse(conflict_hi == 1, 0, conflict_li)) %>% 
  mutate(year = as.numeric(year))

         # ceck if there are other missing ones, 751 is OK
conf %>% 
  select(country, gwno_loc) %>% 
  filter(is.na(country)) %>% 
  distinct()

# GBIF membership
gbif <- gbif %>% 
  filter(status != "Other associate participant") %>% 
  select(name, member_since) %>% 
  mutate(country = countrycode(name, origin = "country.name", destination = "iso3c")) %>% 
  select(-name, year = member_since) 


li <- unique(gbif$country)
gbif_out <- list()

for(i in 1:length(li)){
  sub <- gbif %>%  filter(country == li[i])
  gbif_out[[i]] <- tibble(year = sub$year:2019,
                     country = li[i], 
                     gbif_member = 1)
  
  
}
gbif_out <- bind_rows(gbif_out)

# merge with other vdem variables

out <- old %>% 
  select(-confl) %>% 
  left_join(conf, by = c("year", "country")) %>%
  mutate(conflict_hi = ifelse(is.na(conflict_hi) & year <= 2019, 0, conflict_hi)) %>% 
  mutate(conflict_li = ifelse(is.na(conflict_li) & year <= 2019, 0, conflict_li))# %>% 
 # left_join(gbif_out, by = c("year", "country")) %>% 
#  mutate(gbif_member = ifelse(is.na(gbif_member), 0, gbif_member))

# write to disk
write_csv(out, "output/vdem_variables_conflict.csv")


