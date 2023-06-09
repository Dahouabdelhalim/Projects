## Calculation of RRT and RRN

# This script provides calculation of the warming and N-addition response ratios.

library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(ggpmisc) # needed to add fromula to plot (stat_poly_eq_())
library(ggpubr) # can be used to combine and arrange multiple plots
library(cowplot)
#--------------------------


#### RRT calculation ####
# Instead of  Q10 calculation according to Van't Hoff equation we calculate a simple response ratio.
# Calculate RRT of multiple columns at once
multipleRRT <- sample_data_unordered %>% 
  filter(X0N_.N==0) %>% group_by(site, landuse, replicate) %>% 
  summarise_at(c("cumul_resp_incub_ngC_g_total_time",  
                 "qCO2_ugC_mgCmic_h",
                 "cumul_resp_Corg",
                 "Cresp_CUE_ngC_gDW_h", 
                 "total_DNA_extracted_ugDNA_gDW",
                 "DNAprod_ng_g_h",
                 "Cgrowth_ngC_gDW_h",
                 "CUE",
                 "mass_specific_growth_rate_1perd",
                 "turnover_d",
                 "Cmic_ugC_gDW",
                 "nF_ugC_gDW",
                 "Nmic_ugN_gDW",
                 "nF_ugN_gDW",
                 "CN_mic"), 
               funs((.[treatment=="T20"]/.[treatment=="T10"]))) %>%
  ungroup()
multipleRRT$RR <- c("T")

#### RRN calculation ####
# Calculate RRN of multiple columns at once
multipleRRN <- sample_data_unordered %>% 
  filter(Temp==20) %>% group_by(site, landuse, replicate) %>% 
  summarise_at(c("cumul_resp_incub_ngC_g_total_time",  
                 "qCO2_ugC_mgCmic_h",
                 "cumul_resp_Corg",
                 "Cresp_CUE_ngC_gDW_h", 
                 "total_DNA_extracted_ugDNA_gDW",
                 "DNAprod_ng_g_h",
                 "Cgrowth_ngC_gDW_h",
                 "CUE",
                 "mass_specific_growth_rate_1perd",
                 "turnover_d",
                 "Cmic_ugC_gDW",
                 "nF_ugC_gDW",
                 "Nmic_ugN_gDW",
                 "nF_ugN_gDW",
                 "CN_mic"), 
               funs((.[treatment=="T20+N"]/.[treatment=="T20"]))) %>%
  ungroup()
multipleRRN$RR <- c("N")

# combine RRT and RRN table
multipleRRTN <- rbind(multipleRRT, multipleRRN)



#### summarise results ####
# descriptive statistics on RRT and RRN
multipleRRTN.summary <- multipleRRTN %>%  # have to rethink how to remove na's
  group_by(RR, landuse) %>%
  summarise_if(is.numeric, list(
    min = ~ min(.),
    max = ~ max(.),
    mean = ~ mean(.),
    sd = ~ sd(.)))

## fast overview on RRT and RRN data for results part
# CUE
# Cresp_CUE_ngC_gDW_h
# Cgrowth_ngC_gDW_h
# cumul_resp_incub_ngC_g_total_time
# Cmic_ugC_gDW   == Cmic.end
# CN_mic
# turnover_d
multipleRRTN %>%  
  group_by(RR, landuse, site) %>% 
  summarise(min = min(Cgrowth_ngC_gDW_h, na.rm=T),
            max = max(Cgrowth_ngC_gDW_h, na.rm=T),
            median = median(Cgrowth_ngC_gDW_h, na.rm=T),
            mean = mean(Cgrowth_ngC_gDW_h, na.rm=T),
            sd = sd(Cgrowth_ngC_gDW_h, na.rm=T),
            se = sd(Cgrowth_ngC_gDW_h, na.rm=T)/sqrt(length(Cgrowth_ngC_gDW_h)))

# Test if RRT and RRN are significantly different to 1
multipleRRTN.t.test <- multipleRRTN[,c(1:17,19)] %>%  
  group_by(RR, landuse) %>%
  summarise_if(is.numeric, list(P = ~ t.test(na.omit(.), mu = 1)$p.value))

# Create table indicating significances
multipleRRTN.t.test.sig <- multipleRRTN.t.test %>%  
  group_by(RR, landuse) %>%
  summarise_all(list(sig = ~ ifelse(. < 0.05, "*", "ns")))
multipleRRTN.t.test.sig
