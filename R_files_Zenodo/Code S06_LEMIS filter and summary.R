
library(dplyr)

# LEMIS data downloaded directly from the Zenodo 10.5281/zenodo.3565869 repo as
# a csv. dplyr v.1.1 did not work (missing function) with fstplyr a dependent
# for the lemis package.

lemisData <- read.csv("./Data/TradeDatabases/lemis_2000_2014_cleaned.csv")

names(lemisData)

lemisAraData <- lemisData[which(lemisData$class == "Arachnida"),]

write.csv(x = lemisAraData,
          file = "./Data/TradeDatabases/LEMIS_arachnid_data.csv",
          row.names = FALSE)

lemisAraData <- read.csv(file = "./Data/TradeDatabases/LEMIS_arachnid_data.csv",
                      stringsAsFactors = FALSE)

# there are this many genera listed only to genus level, check that all are
# covered by specific species
percentageGenusListings <- lemisAraData %>% 
  filter(!genus == "Non-CITES entry" & !is.na(genus)) %>% 
  mutate(genLevel = str_detect(species, "sp\\\\.$")) %>% 
  group_by(genus) %>%
  count(genLevel) %>% 
  summarise(perGenLevel = n[genLevel] / (n[!genLevel] + n[genLevel])*100 )

percentageGenusListings %>% 
  filter(perGenLevel > 99)
# This is suggesting that there are no genera that are only represented by sp.
# level listings.
percentageGenusListings %>% 
  ungroup() %>% 
  summarise(mean = mean(perGenLevel),
            min = min(perGenLevel))

# Raw species count from LEMIS, excluding genus only listings, but not tackling
# synonyms
lemisAraData %>% 
  mutate(lemisName = paste(genus, species)) %>% 
  filter(!str_detect(lemisName, "sp\\\\.$")) %>% # filter out genus only listings for species count
  pull(lemisName) %>% 
  unique() %>% 
  length()

# raw genus count
lemisAraData %>% 
  pull(genus) %>% 
  unique() %>% 
  length()
