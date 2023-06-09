## linking ant communities to seed removal - partial least squares regression analyses
## subterranean traps and below ground seed removal

## -- load in the libraries ####
library(here)
library(tidyverse)
library(plyr)
library(pls)
library(reshape)

## -- set folder paths ####
raw_data_path <- here::here("data/raw_data/csv")
figure_path <- here::here("figs")


## -- load in the data files ####

Ruzi_etal_data_seed_removal_below_totals <- read_csv(paste(raw_data_path, "Ruzi_etal_data_seed_removal_below_totals.csv", sep = "/"), 
                               col_types = cols(
                                 Species = col_character(),
                                 Season = col_character(),
                                 Dispersion = col_character(),
                                 Num.weeks = col_double(),
                                 Plot = col_character(),
                                 Direction = col_character(),
                                 Replicate = col_double(),
                                 Total.Recovered = col_double(),
                                 Total.Removed = col_double(),
                                 Partial.Recovered = col_double(),
                                 Vials.Present = col_character(),
                                 Dormancy = col_character()))
Ruzi_etal_data_seed_removal_below_totals


Ruzi_etal_data_below_frequency_condensed <- read_csv(paste(raw_data_path, "Ruzi_etal_data_below_frequency_condensed.csv", sep = "/"),
                                                     col_types = cols(
                                                       .default = col_double(),
                                                       ID = col_character(),
                                                       Plot_abbre = col_character(),
                                                       Season = col_character(),
                                                       Species = col_character(),
                                                       Experiment = col_character()))

Ruzi_etal_data_below_frequency_condensed


## -- process the seed removal data to get seed removal per species and per glass ####

sr_below <- Ruzi_etal_data_seed_removal_below_totals %>%
  # make a proportion and a percent column
  mutate(proportion = Total.Removed/(Total.Recovered + Total.Removed),
         percent = proportion * 100) %>%
  mutate(Plot_abbre = if_else(Plot == "25Ha", "25", Plot)) %>%
  mutate(Plot_abbre = if_else(Plot == "AVA", "A", Plot_abbre)) %>%
  mutate(Plot_abbre = if_else(Plot_abbre == "Drayton", "D", Plot_abbre)) %>%
  mutate(Plot_abbre = if_else(Plot_abbre == "Pearson", "P", Plot_abbre)) %>%
  mutate(Plot_abbre = if_else(Plot_abbre == "Zetek", "Z", Plot_abbre)) %>%
  mutate(Season = if_else(Season == "Wet", "wet", Season))
sr_below

## -- getting the seed removal percents by plot for each cache type

#- Ape

sr_Apeb <- sr_below %>%
  filter(Species == "ApMe")

sr_Apeb_plot <- ddply(sr_Apeb, "Plot_abbre", summarise,
                      n = sum(!is.na(percent)),
                      Mean = mean(percent, na.rm = TRUE))

sr_Apeb_plot$ID <- c(paste("Ape", sr_Apeb_plot$Plot_abbre, sep = "_"))
sr_Apeb_plot  


#- Cec

sr_Cecb <- sr_below %>%
  filter(Species == "CeLo")

sr_Cecb_plot <- ddply(sr_Cecb, "Plot_abbre", summarise,
                      n = sum(!is.na(percent)),
                      Mean = mean(percent, na.rm = TRUE))

sr_Cecb_plot$ID <- c(paste("Cec", sr_Cecb_plot$Plot_abbre, sep = "_"))
sr_Cecb_plot  

#- Jac

sr_Jacb <- sr_below %>%
  filter(Species == "JaCo")

sr_Jacb_plot <- ddply(sr_Jacb, "Plot_abbre", summarise,
                      n = sum(!is.na(percent)),
                      Mean = mean(percent, na.rm = TRUE))

sr_Jacb_plot$ID <- c(paste("Jac", sr_Jacb_plot$Plot_abbre, sep = "_"))
sr_Jacb_plot  

#- Och

sr_Ochb <- sr_below %>%
  filter(Species == "OcPy")

sr_Ochb_plot <- ddply(sr_Ochb, "Plot_abbre", summarise,
                      n = sum(!is.na(percent)),
                      Mean = mean(percent, na.rm = TRUE))

sr_Ochb_plot$ID <- c(paste("Och", sr_Ochb_plot$Plot_abbre, sep = "_"))
sr_Ochb_plot  

#- TrBl

sr_TrBlb <- sr_below %>%
  filter(Species == "TrBl")

sr_TrBlb_plot <- ddply(sr_TrBlb, "Plot_abbre", summarise,
                       n = sum(!is.na(percent)),
                       Mean = mean(percent, na.rm = TRUE))

sr_TrBlb_plot$ID <- c(paste("TrBl", sr_TrBlb_plot$Plot_abbre, sep = "_"))
sr_TrBlb_plot  

#- Zan

sr_Zanb <- sr_below %>%
  filter(Species == "ZaEk")

sr_Zanb_plot <- ddply(sr_Zanb, "Plot_abbre", summarise,
                      n = sum(!is.na(percent)),
                      Mean = mean(percent, na.rm = TRUE))

sr_Zanb_plot$ID <- c(paste("Zan", sr_Zanb_plot$Plot_abbre, sep = "_"))
sr_Zanb_plot  

#- Glass

sr_Glassb <- sr_below %>%
  filter(Species == "Glass")

sr_Glassb_plot <- ddply(sr_Glassb, "Plot_abbre", summarise,
                        n = sum(!is.na(percent)),
                        Mean = mean(percent, na.rm = TRUE))

sr_Glassb_plot$ID <- c(paste("Glass", sr_Apeb_plot$Plot_abbre, sep = "_"))
sr_Glassb_plot  

## -- process the ant community data ####

below_glass <- Ruzi_etal_data_below_frequency_condensed %>%
  # remove the passive samples (but leaving in the glass cache samples)
  filter(Species != "P") %>%
  # remove the current ID column which does not match the ID column
  # in the species by plot averages which calculated previously
  select(-ID) %>%
  # duplicate the Species and plot columns
  mutate(species = Species,
         plot = Plot_abbre) %>%
  # make a new ID column
  unite(ID, species, plot, sep = "_")
below_glass


## -- join the two datasets together ####

# first join all the of the species by plot information 
seed_removal_below <- sr_Apeb_plot %>%
  full_join(sr_Cecb_plot) %>%
  full_join(sr_Glassb_plot) %>%
  full_join(sr_Jacb_plot) %>%
  full_join(sr_Ochb_plot) %>%
  full_join(sr_TrBlb_plot) %>%
  full_join(sr_Zanb_plot)
seed_removal_below

# now need to do a left join so that only the cache types that have community
# information also get the seed removal information into the dataframe
below_glass_sr <- below_glass %>%
  left_join(seed_removal_below)
below_glass_sr

names(below_glass_sr) # ant community columns are 1:20


below_glass_mat <- as.matrix(below_glass_sr[,1:20])
dimnames(below_glass_mat)
dimnames(below_glass_mat)[[1]] <- below_glass_sr$ID

head(below_glass_mat)
ncol(below_glass_mat)

## putting the response and community together
below_glass_community_data <- data.frame(removal = below_glass_sr$Mean,
                                         community = I(below_glass_mat))
str(below_glass_community_data)

## -- run the PLS-R analyses ####
below_glass_plsr1 <- plsr(removal ~ community, data=below_glass_community_data,
                          validation = "LOO")
below_glass_plsr1
summary(below_glass_plsr1) # 5 comps based on criteria mentioned in manuscript

below_glass_plsr2 <- plsr(removal ~ community, ncomp=5, data=below_glass_community_data,
                          validation = "LOO") 
below_glass_plsr2
summary(below_glass_plsr2)
#Data: 	X dimension: 27 20 
#Y dimension: 27 1
#Fit method: kernelpls
#Number of components considered: 5

#VALIDATION: RMSEP
#Cross-validated using 27 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps  4 comps  5 comps
#CV           22.92    22.29    24.07    24.37    25.23    27.31
#adjCV        22.92    22.25    23.90    24.25    25.01    26.76

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps  4 comps  5 comps
#X          25.06    34.62    53.76    62.82    66.36
#removal    26.43    47.08    53.69    63.04    74.50


RMSEP(below_glass_plsr2)
#       (Intercept)  1 comps  2 comps  3 comps  4 comps  5 comps
#CV           22.92    22.29    24.07    24.37    25.23    27.31
#adjCV        22.92    22.25    23.90    24.25    25.01    26.76






