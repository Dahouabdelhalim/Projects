#### Supplemental Code for Madsen, Vander Meiden & Shizuka: Social partners and temperature jointly affect morning foraging activity of small birds in winter#####

### PART 2: Run Gaussian Mixture Model#####

#runs Gaussian Mixture Model from all visitation data, then saves the output

#required libraries
library(tidyverse)
library(asnipe)
library(lubridate)


all_visits <- read.csv('all_visits.csv')
all_visits <- all_visits %>%
  ungroup() %>%
  mutate(index = as.numeric(Date- min(Date)+1)) %>%
  mutate(Timestamp = as.numeric(Datetime))

dat.dowo <- all_visits %>%
  ungroup() %>%
  filter(Species == "DOWO") %>%
  mutate(LoggerDate = paste0(Logger, Date))
dat.wbnu <- all_visits %>%
  ungroup() %>%
  filter(Species == "WBNU") %>%
  mutate(LoggerDate = paste0(Logger, Date))

#Run Gaussian Mixture Model. Note that this is a stochastic process so end result can be slightly different each time.
dowoflocks=gmmevents(dat.dowo$Timestamp, dat.dowo$RFID, dat.dowo$LoggerDate)
wbnuflocks=gmmevents(dat.wbnu$Timestamp, dat.wbnu$RFID, dat.wbnu$LoggerDate)

saveRDS(dowoflocks, "conspecificDOWOflocks_final.rds")
saveRDS(wbnuflocks, "conspecificWBNUflocks_final.rds")