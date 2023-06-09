## linking ant communities to seed removal - partial least squares regression analyses
## pitfall traps and above ground seed removal

## -- load in the libraries ####
library(here)
library(tidyverse)
library(plyr)
library(pls)
library(reshape)

## -- set folder paths ####
raw_data_path <- here::here("data/raw_data/csv")
figure_path <- here::here("figs")

## -- load in the data file ####

Ruzi_etal_data_seed_removal_above_totals <- read_csv(paste(raw_data_path, "Ruzi_etal_data_seed_removal_above_totals.csv", sep = "/"), 
                                                     col_types =cols(
                                                       Species = col_character(),
                                                       Dispersion = col_character(),
                                                       Season = col_character(),
                                                       Plot = col_character(),
                                                       Replicate = col_double(),
                                                       Total.Seeds.Remaining = col_double(),
                                                       Total.Seeds.Removed = col_double(),
                                                       Position = col_double(),
                                                       Start.Date = col_character(),
                                                       Dormancy = col_character()))
Ruzi_etal_data_seed_removal_above_totals

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

## -- process the seed removal data to get all seed removal per site per season -- all tree species pooled ####

sr_all <- Ruzi_etal_data_seed_removal_above_totals %>%
  # to remove the silica bead controls from this dataset
  filter(Dispersion != "Control") %>%
  # make a proportion column and a percentage column
  mutate(proportion = Total.Seeds.Removed/(Total.Seeds.Remaining + Total.Seeds.Removed),
         percent = proportion * 100)
sr_all


# subsetting the wet and dry season out separately to that I can summarize seed removal by plot
sr_wet <- sr_all %>%
  filter(Season == "Wet")
sr_wet  

sr_dry <- sr_all %>%
  filter(Season == "Dry")
sr_dry

sr_wet_plot <- ddply(sr_wet, "Plot", summarise,
                     n=sum(!is.na(percent)),
                     Mean=mean(percent, na.rm=TRUE))
sr_wet_plot$Season <- c(rep("Wet",5))
sr_wet_plot

sr_dry_plot <- ddply(sr_dry, "Plot", summarise,
                     n=sum(!is.na(percent)),
                     Mean=mean(percent, na.rm=TRUE))
sr_dry_plot$Season <- c(rep("Dry",5))
sr_dry_plot

seed_removal_db <- sr_dry_plot %>%
  # join the two together
  full_join(sr_wet_plot) %>%
  mutate(Plot_abbre = if_else(Plot == "25Ha", "25", Plot)) %>%
  mutate(Plot_abbre = if_else(Plot == "AVA", "A", Plot_abbre)) %>%
  mutate(Plot_abbre = if_else(Plot_abbre == "Drayton", "D", Plot_abbre)) %>%
  mutate(Plot_abbre = if_else(Plot_abbre == "Pearson", "P", Plot_abbre)) %>%
  mutate(Plot_abbre = if_else(Plot_abbre == "Zetek", "Z", Plot_abbre)) %>%
  mutate(Season = if_else(Season == "Wet", "wet", Season)) %>%
  mutate(Season = if_else(Season == "Dry", "dry", Season)) %>%
  # make a unique identifier and add an order column and rearrange based on it
  mutate(Plot_season = paste(Plot_abbre, Season, sep = "_"),
         order = c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10)) %>%
  arrange(order)
seed_removal_db

## -- process the ant community data ####

Plot_season <- as.vector(Ruzi_etal_freq_data_pitfalls_dry_wet_season$ID)
Plot_season

names(Ruzi_etal_freq_data_pitfalls_dry_wet_season) #  the species start 2:60

pitfalls_mat <- as.matrix(Ruzi_etal_freq_data_pitfalls_dry_wet_season[,2:60])
dimnames(pitfalls_mat)
pitfallsnames_mat <- as.vector(Plot_season)
dimnames(pitfalls_mat)[[1]] <- c(pitfallsnames_mat)
summary(pitfalls_mat) 

head(pitfalls_mat)

## -- join the two datasets together ####

removal_community_data <- data.frame(removal=seed_removal_db$Mean)
removal_community_data
removal_community_data$community <- pitfalls_mat
str(removal_community_data)

## -- run the PLS-R analyses ####

sr1 <- plsr(removal ~ community, ncomp=8, data=removal_community_data,
            validation = "LOO")
sr1
summary(sr1) # using 3 comps 

sr2 <- plsr(removal ~ community, ncomp=3, data=removal_community_data,
            validation = "LOO")
summary(sr2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 3
#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps
#CV           15.57    15.01    18.88    17.79
#adjCV        15.57    14.90    17.94    16.93
#TRAINING: % variance explained
#         1 comps  2 comps  3 comps
#X          43.39    49.40    65.48
#removal    35.47    90.41    96.31

RMSEP(sr2)
#       (Intercept)  1 comps  2 comps  3 comps
#CV           15.57    15.01    18.88    17.79
#adjCV        15.57    14.90    17.94    16.93


## -- running plsr analyses per tree species  -- Ape ####
sr_Ape_dry <- sr_dry %>% filter(Species == "ApMe")
sr_Ape_wet <- sr_wet %>% filter(Species == "ApMe")

sr_Ape_dry_plot <- ddply(sr_Ape_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Ape_dry_plot$Season <- c(rep("dry",5))
sr_Ape_dry_plot

sr_Ape_wet_plot <- ddply(sr_Ape_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Ape_wet_plot$Season <- c(rep("wet",5))
sr_Ape_wet_plot

Ape_sr_dataframe <- rbind(sr_Ape_dry_plot, sr_Ape_wet_plot)
Ape_sr_dataframe$ID <- c(paste("Ape",Ape_sr_dataframe$Plot, sep="_",
                               Ape_sr_dataframe$Season))

Ape_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Ape_sr_dataframe <- arrange(Ape_sr_dataframe, order)
Ape_sr_dataframe

Ape_community_data <- data.frame(removal=Ape_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Ape_community_data)


Ape1 <- plsr(removal ~ community, ncomp=8, data=Ape_community_data,
             validation = "LOO")
Ape1
summary(Ape1)

### - going to use 4 comps for Ape

Ape2 <- plsr(removal ~ community, ncomp=4, data=Ape_community_data,
             validation = "LOO")
summary(Ape2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 4

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps  4 comps
#CV           21.15    19.48    18.64    24.51    23.57
#adjCV        21.15    19.31    18.55    23.46    22.39

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps  4 comps
#X          33.99    59.98    70.71    81.46
#removal    43.35    65.06    88.92    97.92

RMSEP(Ape2)
#       (Intercept)  1 comps  2 comps  3 comps  4 comps
#CV           21.15    19.48    18.64    24.51    23.57
#adjCV        21.15    19.31    18.55    23.46    22.39


## -- running plsr analyses per tree species  -- Cec ####
sr_Cec_dry <- subset(sr_dry, Species=="CeLo")
sr_Cec_wet <- subset(sr_wet, Species=="CeLo")

sr_Cec_dry_plot <- ddply(sr_Cec_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Cec_dry_plot$Season <- c(rep("dry",5))
sr_Cec_dry_plot

sr_Cec_wet_plot <- ddply(sr_Cec_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Cec_wet_plot$Season <- c(rep("wet",5))
sr_Cec_wet_plot

Cec_sr_dataframe <- rbind(sr_Cec_dry_plot, sr_Cec_wet_plot)
Cec_sr_dataframe$ID <- c(paste("Cec",Cec_sr_dataframe$Plot, sep="_",
                               Cec_sr_dataframe$Season))

Cec_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Cec_sr_dataframe <- arrange(Cec_sr_dataframe, order)
Cec_sr_dataframe

Cec_community_data <- data.frame(removal=Cec_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Cec_community_data)


Cec1 <- plsr(removal ~ community, ncomp=8, data=Cec_community_data,
             validation = "LOO")
Cec1
summary(Cec1)

## choosing 3 comps for this
Cec2 <- plsr(removal ~ community, ncomp=3, data=Cec_community_data,
             validation = "LOO")

summary(Cec2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 3

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps
#CV           36.53    42.06    48.46    45.36
#adjCV        36.53    42.19    46.55    43.25

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps
#X          36.43    52.39    66.03
#removal    35.83    81.27    93.95

RMSEP(Cec2)
#       (Intercept)  1 comps  2 comps  3 comps
#CV           36.53    42.06    48.46    45.36
#adjCV        36.53    42.19    46.55    43.25


## -- running plsr analyses per tree species  -- Coc ####
sr_Coc_dry <- subset(sr_dry, Species=="CoVi")
sr_Coc_wet <- subset(sr_wet, Species=="CoVi")

sr_Coc_dry_plot <- ddply(sr_Coc_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Coc_dry_plot$Season <- c(rep("dry",5))
sr_Coc_dry_plot

sr_Coc_wet_plot <- ddply(sr_Coc_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Coc_wet_plot$Season <- c(rep("wet",5))
sr_Coc_wet_plot

Coc_sr_dataframe <- rbind(sr_Coc_dry_plot, sr_Coc_wet_plot)
Coc_sr_dataframe$ID <- c(paste("Coc",Coc_sr_dataframe$Plot, sep="_",
                               Coc_sr_dataframe$Season))

Coc_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Coc_sr_dataframe <- arrange(Coc_sr_dataframe, order)
Coc_sr_dataframe

Coc_community_data <- data.frame(removal=Coc_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Coc_community_data)


Coc1 <- plsr(removal ~ community, ncomp=8, data=Coc_community_data,
             validation = "LOO")
Coc1
summary(Coc1)

## using 3 comps for this one

Coc2 <- plsr(removal ~ community, ncomp=3, data=Coc_community_data,
             validation = "LOO")
summary(Coc2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 3

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps
#CV           32.96    34.05    42.85    42.34
#adjCV        32.96    33.76    41.17    40.19

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps
#X          41.90    54.02    67.28
#removal    32.85    77.25    95.56


RMSEP(Coc2)
#       (Intercept)  1 comps  2 comps  3 comps
#CV           32.96    34.05    42.85    42.34
#adjCV        32.96    33.76    41.17    40.19


## -- running plsr analyses per tree species  -- Fic ####
sr_Fic_dry <- subset(sr_dry, Species=="FiIn")
sr_Fic_wet <- subset(sr_wet, Species=="FiIn")

sr_Fic_dry_plot <- ddply(sr_Fic_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Fic_dry_plot$Season <- c(rep("dry",5))
sr_Fic_dry_plot

sr_Fic_wet_plot <- ddply(sr_Fic_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Fic_wet_plot$Season <- c(rep("wet",5))
sr_Fic_wet_plot

Fic_sr_dataframe <- rbind(sr_Fic_dry_plot, sr_Fic_wet_plot)
Fic_sr_dataframe$ID <- c(paste("Fic",Fic_sr_dataframe$Plot, sep="_",
                               Fic_sr_dataframe$Season))

Fic_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Fic_sr_dataframe <- arrange(Fic_sr_dataframe, order)
Fic_sr_dataframe

Fic_community_data <- data.frame(removal=Fic_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Fic_community_data)


Fic1 <- plsr(removal ~ community, ncomp=8, data=Fic_community_data,
             validation = "LOO")
Fic1
summary(Fic1)
## using 3 comps for this one

Fic2 <- plsr(removal ~ community, ncomp=3, data=Fic_community_data,
             validation = "LOO")
summary(Fic2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 3

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps
#CV           30.85    44.40    48.05    47.75
#adjCV        30.85    43.01    45.87    45.30

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps
#X          28.09    49.46    67.16
#removal    40.87    85.85    97.09


RMSEP(Fic2)
#       (Intercept)  1 comps  2 comps  3 comps
#CV           30.85    44.40    48.05    47.75
#adjCV        30.85    43.01    45.87    45.30


## -- running plsr analyses per tree species  -- Gua ####
sr_Gua_dry <- subset(sr_dry, Species=="GuUl")
sr_Gua_wet <- subset(sr_wet, Species=="GuUl")

sr_Gua_dry_plot <- ddply(sr_Gua_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Gua_dry_plot$Season <- c(rep("dry",5))
sr_Gua_dry_plot

sr_Gua_wet_plot <- ddply(sr_Gua_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Gua_wet_plot$Season <- c(rep("wet",5))
sr_Gua_wet_plot

Gua_sr_dataframe <- rbind(sr_Gua_dry_plot, sr_Gua_wet_plot)
Gua_sr_dataframe$ID <- c(paste("Gua",Gua_sr_dataframe$Plot, sep="_",
                               Gua_sr_dataframe$Season))

Gua_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Gua_sr_dataframe <- arrange(Gua_sr_dataframe, order)
Gua_sr_dataframe

Gua_community_data <- data.frame(removal=Gua_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Gua_community_data)


Gua1 <- plsr(removal ~ community, ncomp=8, data=Gua_community_data,
             validation = "LOO")
Gua1
summary(Gua1)

## using 4 comps for this

Gua2 <- plsr(removal ~ community, ncomp=4, data=Gua_community_data,
             validation = "LOO")
summary(Gua2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 4

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps  4 comps
#CV           32.66    31.96    37.80    40.17    40.05
#adjCV        32.66    31.46    36.48    38.31    37.99

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps  4 comps
#X          38.39    62.35    71.27    81.52
#removal    57.31    76.47    94.12    99.48

RMSEP(Gua2)
#       (Intercept)  1 comps  2 comps  3 comps  4 comps
#CV           32.66    31.96    37.80    40.17    40.05
#adjCV        32.66    31.46    36.48    38.31    37.99

## -- running plsr analyses per tree species  -- Hie ####
sr_Hie_dry <- subset(sr_dry, Species=="HiAl")
sr_Hie_wet <- subset(sr_wet, Species=="HiAl")

sr_Hie_dry_plot <- ddply(sr_Hie_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Hie_dry_plot$Season <- c(rep("dry",5))
sr_Hie_dry_plot

sr_Hie_wet_plot <- ddply(sr_Hie_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Hie_wet_plot$Season <- c(rep("wet",5))
sr_Hie_wet_plot

Hie_sr_dataframe <- rbind(sr_Hie_dry_plot, sr_Hie_wet_plot)
Hie_sr_dataframe$ID <- c(paste("Hie",Hie_sr_dataframe$Plot, sep="_",
                               Hie_sr_dataframe$Season))

Hie_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Hie_sr_dataframe <- arrange(Hie_sr_dataframe, order)
Hie_sr_dataframe

Hie_community_data <- data.frame(removal=Hie_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Hie_community_data)


Hie1 <- plsr(removal ~ community, ncomp=8, data=Hie_community_data,
             validation = "LOO")
Hie1
summary(Hie1)

## using 3 comps for this

Hie2 <- plsr(removal ~ community, ncomp=3, data=Hie_community_data,
             validation = "LOO")
summary(Hie2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 3

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps
#CV           32.51    28.03    33.33    29.64
#adjCV        32.51    27.80    31.71    28.17

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps
#X          43.53    50.92    62.24
#removal    47.69    90.05    98.08

RMSEP(Hie2)
#       (Intercept)  1 comps  2 comps  3 comps
#CV           32.51    28.03    33.33    29.64
#adjCV        32.51    27.80    31.71    28.17


## -- running plsr analyses per tree species  -- Jac ####
sr_Jac_dry <- subset(sr_dry, Species=="JaCo")
sr_Jac_wet <- subset(sr_wet, Species=="JaCo")

sr_Jac_dry_plot <- ddply(sr_Jac_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Jac_dry_plot$Season <- c(rep("dry",5))
sr_Jac_dry_plot

sr_Jac_wet_plot <- ddply(sr_Jac_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Jac_wet_plot$Season <- c(rep("wet",5))
sr_Jac_wet_plot

Jac_sr_dataframe <- rbind(sr_Jac_dry_plot, sr_Jac_wet_plot)
Jac_sr_dataframe$ID <- c(paste("Jac",Jac_sr_dataframe$Plot, sep="_",
                               Jac_sr_dataframe$Season))

Jac_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Jac_sr_dataframe <- arrange(Jac_sr_dataframe, order)
Jac_sr_dataframe

Jac_community_data <- data.frame(removal=Jac_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Jac_community_data)


Jac1 <- plsr(removal ~ community, ncomp=8, data=Jac_community_data,
             validation = "LOO")
Jac1
summary(Jac1)

## using 4 comps for this one

Jac2 <- plsr(removal ~ community, ncomp=4, data=Jac_community_data,
             validation = "LOO")
summary(Jac2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 4

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps  4 comps
#CV            32.2    35.86    37.55    42.36    40.39
#adjCV         32.2    35.44    36.81    40.47    38.37

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps  4 comps
#X          33.36    57.34    68.44    81.61
#removal    37.74    64.64    90.39    97.99

RMSEP(Jac2)
#       (Intercept)  1 comps  2 comps  3 comps  4 comps
#CV            32.2    35.86    37.55    42.36    40.39
#adjCV         32.2    35.44    36.81    40.47    38.37



## -- running plsr analyses per tree species  -- Lue ####
sr_Lue_dry <- subset(sr_dry, Species=="LuSe")
sr_Lue_wet <- subset(sr_wet, Species=="LuSe")

sr_Lue_dry_plot <- ddply(sr_Lue_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Lue_dry_plot$Season <- c(rep("dry",5))
sr_Lue_dry_plot

sr_Lue_wet_plot <- ddply(sr_Lue_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Lue_wet_plot$Season <- c(rep("wet",5))
sr_Lue_wet_plot

Lue_sr_dataframe <- rbind(sr_Lue_dry_plot, sr_Lue_wet_plot)
Lue_sr_dataframe$ID <- c(paste("Lue",Lue_sr_dataframe$Plot, sep="_",
                               Lue_sr_dataframe$Season))

Lue_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Lue_sr_dataframe <- arrange(Lue_sr_dataframe, order)
Lue_sr_dataframe

Lue_community_data <- data.frame(removal=Lue_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Lue_community_data)


Lue1 <- plsr(removal ~ community, ncomp=8, data=Lue_community_data,
             validation = "LOO")
Lue1
summary(Lue1)

## using 3 comps for this one

Lue2 <- plsr(removal ~ community, ncomp=3, data=Lue_community_data,
             validation = "LOO")
summary(Lue2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 3

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps
#CV           18.29    19.11    25.37    23.74
#adjCV        18.29    18.81    24.14    22.54

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps
#X          43.73    50.21    60.90
#removal    45.18    88.19    96.69

RMSEP(Lue2)
#       (Intercept)  1 comps  2 comps  3 comps
#CV           18.29    19.11    25.37    23.74
#adjCV        18.29    18.81    24.14    22.54

## -- running plsr analyses per tree species  -- Och ####

sr_Och_dry <- subset(sr_dry, Species=="OcPy")
sr_Och_wet <- subset(sr_wet, Species=="OcPy")

sr_Och_dry_plot <- ddply(sr_Och_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Och_dry_plot$Season <- c(rep("dry",5))
sr_Och_dry_plot

sr_Och_wet_plot <- ddply(sr_Och_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Och_wet_plot$Season <- c(rep("wet",5))
sr_Och_wet_plot

Och_sr_dataframe <- rbind(sr_Och_dry_plot, sr_Och_wet_plot)
Och_sr_dataframe$ID <- c(paste("Och",Och_sr_dataframe$Plot, sep="_",
                               Och_sr_dataframe$Season))

Och_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Och_sr_dataframe <- arrange(Och_sr_dataframe, order)
Och_sr_dataframe

Och_community_data <- data.frame(removal=Och_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Och_community_data)


Och1 <- plsr(removal ~ community, ncomp=8, data=Och_community_data,
             validation = "LOO")
Och1
summary(Och1)

## using 3 comps for this one

Och2 <- plsr(removal ~ community, ncomp=3, data=Och_community_data,
             validation = "LOO")
summary(Och2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 3

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps
#CV            21.4    19.94    25.75    24.35
#adjCV         21.4    19.76    24.60    23.07

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps
#X          43.59    49.63    57.14
#removal    41.83    87.51    97.44

RMSEP(Och2)
#       (Intercept)  1 comps  2 comps  3 comps
#CV            21.4    19.94    25.75    24.35
#adjCV         21.4    19.76    24.60    23.07



## -- running plsr analyses per tree species  -- TrBl ####

sr_TrBl_dry <- subset(sr_dry, Species=="TrBl")
sr_TrBl_wet <- subset(sr_wet, Species=="TrBl")

sr_TrBl_dry_plot <- ddply(sr_TrBl_dry, "Plot", summarise,
                          n=sum(!is.na(percent)),
                          Mean=mean(percent, na.rm=TRUE))
sr_TrBl_dry_plot$Season <- c(rep("dry",5))
sr_TrBl_dry_plot

sr_TrBl_wet_plot <- ddply(sr_TrBl_wet, "Plot", summarise,
                          n=sum(!is.na(percent)),
                          Mean=mean(percent, na.rm=TRUE))
sr_TrBl_wet_plot$Season <- c(rep("wet",5))
sr_TrBl_wet_plot

TrBl_sr_dataframe <- rbind(sr_TrBl_dry_plot, sr_TrBl_wet_plot)
TrBl_sr_dataframe$ID <- c(paste("TrBl",TrBl_sr_dataframe$Plot, sep="_",
                                TrBl_sr_dataframe$Season))

TrBl_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

TrBl_sr_dataframe <- arrange(TrBl_sr_dataframe, order)
TrBl_sr_dataframe

TrBl_community_data <- data.frame(removal=TrBl_sr_dataframe$Mean, 
                                  community=I(pitfalls_mat))
str(TrBl_community_data)


TrBl1 <- plsr(removal ~ community, ncomp=8, data=TrBl_community_data,
              validation = "LOO")
TrBl1
summary(TrBl1)

## using 3 comps for this one

TrBl2 <- plsr(removal ~ community, ncomp=3, data=TrBl_community_data,
              validation = "LOO")
summary(TrBl2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 3

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps
#CV           18.93    18.90    15.85    17.70
#adjCV        18.93    18.18    15.52    16.84

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps
#X          23.72     63.3    69.15
#removal    73.60     80.0    96.72


RMSEP(TrBl2)
#       (Intercept)  1 comps  2 comps  3 comps
#CV           18.93    18.90    15.85    17.70
#adjCV        18.93    18.18    15.52    16.84



## -- running plsr analyses per tree species  -- TrBr ####
sr_TrBr_dry <- subset(sr_dry, Species=="TrBr")
sr_TrBr_wet <- subset(sr_wet, Species=="TrBr")

sr_TrBr_dry_plot <- ddply(sr_TrBr_dry, "Plot", summarise,
                          n=sum(!is.na(percent)),
                          Mean=mean(percent, na.rm=TRUE))
sr_TrBr_dry_plot$Season <- c(rep("dry",5))
sr_TrBr_dry_plot

sr_TrBr_wet_plot <- ddply(sr_TrBr_wet, "Plot", summarise,
                          n=sum(!is.na(percent)),
                          Mean=mean(percent, na.rm=TRUE))
sr_TrBr_wet_plot$Season <- c(rep("wet",5))
sr_TrBr_wet_plot

TrBr_sr_dataframe <- rbind(sr_TrBr_dry_plot, sr_TrBr_wet_plot)
TrBr_sr_dataframe$ID <- c(paste("TrBr",TrBr_sr_dataframe$Plot, sep="_",
                                TrBr_sr_dataframe$Season))

TrBr_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

TrBr_sr_dataframe <- arrange(TrBr_sr_dataframe, order)
TrBr_sr_dataframe

TrBr_community_data <- data.frame(removal=TrBr_sr_dataframe$Mean, 
                                  community=I(pitfalls_mat))
str(TrBr_community_data)


TrBr1 <- plsr(removal ~ community, ncomp=8, data=TrBr_community_data,
              validation = "LOO")
TrBr1
summary(TrBr1)

## using 5 components for this one

TrBr2 <- plsr(removal ~ community, ncomp=5, data=TrBr_community_data,
              validation = "LOO")
summary(TrBr2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 5

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps  4 comps  5 comps
#CV           22.33    31.52    29.81    32.23    33.35    39.75
#adjCV        22.33    30.45    28.57    30.99    31.79    37.78

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps  4 comps  5 comps
#X          18.00    36.35    70.27    75.22    82.36
#removal    46.77    72.93    80.46    92.53    97.55

RMSEP(TrBr2)
#       (Intercept)  1 comps  2 comps  3 comps  4 comps  5 comps
#CV           22.33    31.52    29.81    32.23    33.35    39.75
#adjCV        22.33    30.45    28.57    30.99    31.79    37.78


## -- running plsr analyses per tree species  -- Zan ####
sr_Zan_dry <- subset(sr_dry, Species=="ZaEk")
sr_Zan_wet <- subset(sr_wet, Species=="ZaEk")

sr_Zan_dry_plot <- ddply(sr_Zan_dry, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Zan_dry_plot$Season <- c(rep("dry",5))
sr_Zan_dry_plot

sr_Zan_wet_plot <- ddply(sr_Zan_wet, "Plot", summarise,
                         n=sum(!is.na(percent)),
                         Mean=mean(percent, na.rm=TRUE))
sr_Zan_wet_plot$Season <- c(rep("wet",5))
sr_Zan_wet_plot

Zan_sr_dataframe <- rbind(sr_Zan_dry_plot, sr_Zan_wet_plot)
Zan_sr_dataframe$ID <- c(paste("Zan",Zan_sr_dataframe$Plot, sep="_",
                               Zan_sr_dataframe$Season))

Zan_sr_dataframe$order <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10) 

Zan_sr_dataframe <- arrange(Zan_sr_dataframe, order)
Zan_sr_dataframe

Zan_community_data <- data.frame(removal=Zan_sr_dataframe$Mean, 
                                 community=I(pitfalls_mat))
str(Zan_community_data)


Zan1 <- plsr(removal ~ community, ncomp=8, data=Zan_community_data,
             validation = "LOO")
Zan1
summary(Zan1)

## using 4 comps for this one

Zan2 <- plsr(removal ~ community, ncomp=4, data=Zan_community_data,
             validation = "LOO")
summary(Zan2)
#Data: 	X dimension: 10 59 
#Y dimension: 10 1
#Fit method: kernelpls
#Number of components considered: 4

#VALIDATION: RMSEP
#Cross-validated using 10 leave-one-out segments.
#       (Intercept)  1 comps  2 comps  3 comps  4 comps
#CV           22.97    29.01    28.57    29.24    28.22
#adjCV        22.97    28.33    27.73    28.00    26.86

#TRAINING: % variance explained
#         1 comps  2 comps  3 comps  4 comps
#X          40.43    57.64    70.16    77.88
#removal    39.26    70.51    89.57    97.29

RMSEP(Zan2)
#       (Intercept)  1 comps  2 comps  3 comps  4 comps
#CV           22.97    29.01    28.57    29.24    28.22
#adjCV        22.97    28.33    27.73    28.00    26.86

