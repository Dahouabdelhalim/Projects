### Code for processing raw data files from the cubic aerosol measurement platform (CAMP)
#
# Author: Christian Pilz, contact: pilz@tropos.de
#
# Valid for CAMP configuration on MOSAiC 
# Technical information about CAMP are provided by Pilz et al. 2022 (https://doi.org/10.5194/amt-15-6889-2022)
# Information about balloon measurements with CAMP on MOSAiC are provided by Pilz et al. 2023 (in preparation for Nature Scientific Data)
# 
# Single data files from CPCs, POPS, and STAP are combined into one file
# RecalcSTAPFiletoFile_Function_v2.R required for first level processing of STAP data
# Corrections for counting efficiency and sampling losses are applied
# STP correction to p = 1013.25 hPa and T = 273.15 K are applied for the CPCs and the POPS
# Sample air temperature used for STP correction occasionally missing in raw files, 25 C should be used instead
# Barometric altitude is calculated with ground temperatures at start of balloon flight taken from MET-City observations
#
# Code requires manual plausibility check of data!!! 
###


# required libraries to run the code
library(lubridate)
library(tidyverse)

# Select raw data file of the CPCs, e.g. "CAMP2_2020-06-29_07-38-31.txt"
cpc.raw <- read.table(choose.files(), dec = ".",header = T,sep = ",", skip = 7)

# Select raw data file of the POPS, e.g. "HK_20200629x002.csv"
pops.raw <- read.table(choose.files(),sep = ",",dec = ".", header = TRUE)

#STAP recalc files e.g. "200629A0_blk_stp_azumi_tau_K2_1s-avg_proc.dat"
stap.raw <- read_delim(choose.files(), delim="\\t")

# CPC data processing
file.end.time <- as.POSIXct(paste(cpc.raw$Date[nrow(cpc.raw)], cpc.raw$GNSS.time[nrow(cpc.raw)]), tz = "UTC", tryFormats = c("%d/%m/%Y %H:%M:%OS"))
cpc.raw$date <- file.end.time - c((nrow(cpc.raw)-1):0)

cpc.proc <- cpc.raw %>% 
  select(date, p_baro, T_out, RH_out, T_in, RH_in, Conc2, Conc1) %>%
  rename("N_12" = "Conc1", "N_8" = "Conc2") %>% 
  mutate(N_8 = ifelse(date < ymd_hms("2020-07-16 00:00:00"), round(N_8 / 0.76, digits = 0),
                      ifelse(between(date, ymd_hms("2020-07-16 00:00:00"), ymd_hms("2020-07-22 00:00:00")), 
                             round(N_8 / 0.6, digits = 0), 
                             round(N_8 / 0.65, digits = 0))),
         N_12 = round(N_12/0.86, digits = 0)) # corrections accounting for counting efficiency and system losses 

# POPS data processing
dt.POPS <- - 65    # delay of internal POPS time 65s
pops.raw$date <- ymd_hms(as.POSIXct(pops.raw$DateTime + dt.POPS, tz = "UTC", origin = "1970-01-01"))

pops.proc <- pops.raw %>% 
  select(c(45, 30:43)) %>%
  rowwise() %>%
  mutate(N_150 = sum(c_across(2:15), na.rm = T)) %>% # total counts per minute
  relocate(N_150, .after = date) %>%
  mutate_at(.vars = 2:16, .funs = list(~ round(./1.76, digits = 2))) # concentration per bin calculated with constant volumetric flow including counting efficiency and system losses
  
names(pops.proc)[3:16] <- c("151", "166", "182", "200", "219", "249", "296", "391", "530", "774", "1143", "1455", "1875", "2480")

POPS.bins <- data.frame(c(144, 158, 174, 191, 209, 229, 270, 324, 473, 594, 1009, 1294, 1637, 2148, 2864)) # bin limits
dlogDp <- log10(POPS.bins[2:15, 1] / POPS.bins[1:14, 1])

pops.proc[,3:16] <- round(t(t(pops.proc[,3:16])/dlogDp), digits = 1) #dN/dlogDp per bin

# STAP data processing
dt.stap <- 44  # internal STAP time 44s earlier
stap.raw$date <- ymd_hms(stap.raw$date) + dt.stap
stap.proc <- stap.raw %>% 
  select(date, blue, green, red) %>%
  mutate_at(.vars = 2:4, .funs = list(~ round(./0.75, digits = 2))) # sampling loss correction

names(stap.proc)[2:4] <- c("abs_450", "abs_525", "abs_624")

#create one combined file
CAMP <- merge(merge(cpc.proc, pops.proc, by = "date"), stap.proc, by = "date")

# calculate barometric height with ground temperatures from MET-City
T.tbl <- data.frame("date" = c("29.06.2020 08:35", "29.06.2020 12:41", "10.07.2020 10:11", "13.07.2020 11:51", "14.07.2020 11:39", "15.07.2020 07:47", "15.07.2020 13:15", "19.07.2020 11:44", "20.07.2020 13:11", "21.07.2020 07:35", "21.07.2020 13:30", "22.07.2020 05:10", "23.07.2020 09:01", "23.07.2020 12:36", "24.07.2020 07:38", "25.07.2020 14:00", "26.07.2020 13:43", "27.07.2020 07:39"), "T_0" = c(-1.70, -1.20, 0.30, 0.50, -1.00, -0.80, -0.40, 0.25, 0.20, 0.00, 0.25, 0.75, 0.20, -0.25, 0.25, 0.30, 1.10, 0.80))

p_GND <- max(CAMP$p_baro, na.rm = T)
T_GND <- T.tbl %>% 
  mutate(date = dmy_hm(date)) %>% 
  filter(between(date, CAMP$date[1], CAMP$date[nrow(CAMP)]))

CAMP.proc <- CAMP %>% 
  mutate(alt_baro = round(287.058*(T_GND$T_0 + 273.15) / 9.81 * log(p_GND / p_baro), digits = 2)) %>%
  relocate(alt_baro, .after = date) %>%
  mutate_at(.vars = 8:24, .funs = list(~ round(. * (p_baro / 1013.25 * 273.15 / (T_in + 273.15)), digits = 1))) # STP correction for CPCs and POPS (STAP corrected within recalc function)


# save file to directory
date <- as.character(CAMP.proc$date[1])

fname <- str_replace_all(date,"[:]", "-")

dname <- paste(fname,"_CAMP.txt", sep = "")

write.table(CAMP.proc, dname, sep = "\\t", dec = ".", row.names = FALSE)
