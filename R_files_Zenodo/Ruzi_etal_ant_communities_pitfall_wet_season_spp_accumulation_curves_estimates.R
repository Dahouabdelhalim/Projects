## to process raw data for pitfall - aboveground species accumulation curves for the wet season,
## sampling coverage, and species estimates

## libraries installed from cran except for SpadeR ####
# install_github('AnneChao/SpadeR')

## -- load libraries ####
library(tidyverse)
library(here)
library(iNEXT)
library(SpadeR)

## -- set paths ####
# set in a way that if download entire dataset package folder, should be able to use relative pathways to run these analyses
Rsource_path <- here::here("R")
raw_data_path <- here::here("data/raw_data/csv")
rdata_path <- here::here("data/rdata")
figure_path <- here::here("figs")

## -- source scripts ####
# to load in the function that will help convert from raw to processed data
source(paste(Rsource_path, "Ruzi_etal_ant_communities_data_processing_function.R", sep = "/"))

## -- read in the data

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

names(Ruzi_etal_freq_data_pitfalls_dry_wet_season)


## to manipulate the data into a format that can be used with iNEXT ####
## need to transform the data to count of frequencies instead of proportions
## so will need to multiple the ant spp columns to the sample num

# pull out the ant species names
spp_names <- names(Ruzi_etal_freq_data_pitfalls_dry_wet_season[2:60])
spp_names

temp <- Ruzi_etal_freq_data_pitfalls_dry_wet_season %>%
  select(ID, spp_names) %>%
  column_to_rownames(var = "ID")
temp

# to go from the proportion frequencies to frequencies
temp <- Ruzi_etal_freq_data_pitfalls_dry_wet_season$Sample_num*temp  
temp 

## need to bring row names back to ID though 
temp <- temp %>%
  rownames_to_column("ID")

## now to separate out the dry season samples
t2 <- temp %>%
  separate(ID, into = c("plot","season"), sep = "_")
head(t2)

t2_wet <- t2 %>%
  filter(season == "wet") %>%
  select(-season)
head(t2_wet)

## now to transpose the dry season tibble
class(t2_wet) 
t2_wet <- as_tibble(t2_wet)

t3_wet <- t2_wet %>%
  gather(spp_names, key = "spp", value = "presence", factor_key = TRUE) %>%
  spread(key = plot, value = presence)
t3_wet

### now need to separate each one of the sites into it's own list,
### and put the trap numer, which was 16 for all of them at the begining of the vector
## and can remove the 0's
num_traps <- 16

## wet season -- 25Ha
w25Ha <- t3_wet %>%
  filter(.$"25" > 0) %>%
  select("25")
w25Ha


w_25Ha <- c(num_traps, w25Ha$`25`) # add in the sample number
w_25Ha

## wet season -- AVA
wA <- t3_wet %>%
  filter(A > 0) %>%
  select(A)
wA


w_A <- c(num_traps, wA$A) # add in the sample number
w_A

## wet season -- Drayton
wD <- t3_wet %>%
  filter(D > 0) %>%
  select(D)
wD


w_D <- c(num_traps, wD$D) # add in the sample number
w_D

## wet season -- Pearson
wP <- t3_wet %>%
  filter(P > 0) %>%
  select(P)
wP


w_P <- c(num_traps, wP$P) # add in the sample number
w_P

## wet season -- Zetek
wZ <- t3_wet %>%
  filter(Z > 0) %>%
  select(Z)
wZ


w_Z <- c(num_traps, wZ$Z) # add in the sample number
w_Z

## wet season list 
# changing from site names to site numbers
# 25HA = 1
# Ava = 2
# Drayton = 3
# Pearson = 4
# Zetek = 5

wet_list <- list("1" = w_25Ha, "2" = w_A, "3" = w_D,
                 "4" = w_P, "5" = w_Z)
wet_list

## wet season rarefaction and extrapolation ####
# choosing to extrapolate out to double the sample size
t_pit <- seq(1,32) 
t_pit

wet_out <- iNEXT(wet_list, q = c(0), datatype = "incidence_freq",
                 size = t_pit) 
wet_out

## graph the iNEXT results - Figure S3e,f ####
pit_wet_iNEXT_graph2 <- ggiNEXT(wet_out, type = 1) +
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,55, by=5)) +
  scale_x_continuous(expand= c(0,0), breaks=seq(0,35, by=5)) + # expand = c(0,0)
  #scale_color_manual(values = c("red","blue","purple","green")) +
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        #panel.border=element_blank(), #gets rid of square going around the entire graph
        #axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="none",
        legend.key=element_blank())
pit_wet_iNEXT_graph2

pit_wet_CV_iNEXT_graph <- ggiNEXT(wet_out, type = 2) +
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,1, by=.1)) +
  scale_x_continuous(expand= c(0,0), breaks=seq(0,35, by=5)) + # expand = c(0,0)
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        #panel.border=element_blank(), #gets rid of square going around the entire graph
        #axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="none",
        #legend.title=element_blank(),
        #legend.text=element_text(face="bold"),
        #legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank())
pit_wet_CV_iNEXT_graph

## save the species accumulation curves and sampling coverage outputs ####
# commented out to not overwrite

#ggsave("Ruzi_etal_FigureS3e.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = pit_wet_iNEXT_graph2,
#       path = figure_path)
#ggsave("Ruzi_etal_FigureS3f.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = pit_wet_CV_iNEXT_graph,
#       path = figure_path)

## to get species estimates using SpadeR ####
wet_list

## choose k to be the cut off where "rare" is when found in less than 25% of the traps
## and common is then more than that

H25_wet_est <- SpadeR::ChaoSpecies(wet_list$'1', datatype = "incidence_freq", k = 4, conf = 0.95)

A_wet_est <- SpadeR::ChaoSpecies(wet_list$'2', datatype = "incidence_freq", k = 4, conf = 0.95)

D_wet_est <- SpadeR::ChaoSpecies(wet_list$'3', datatype = "incidence_freq", k = 4, conf = 0.95)

P_wet_est <- SpadeR::ChaoSpecies(wet_list$'4', datatype = "incidence_freq", k = 4, conf = 0.95)

Z_wet_est <- SpadeR::ChaoSpecies(wet_list$'5', datatype = "incidence_freq" , k = 4, conf = 0.95)


H25_wet_graphing <- as.data.frame(H25_wet_est$Species_table)
H25_wet_graphing$Experiment <- "wet_25"
A_wet_graphing <- as.data.frame(A_wet_est$Species_table)
A_wet_graphing$Experiment <- "wet_A"
D_wet_graphing <- as.data.frame(D_wet_est$Species_table)
D_wet_graphing$Experiment <- "wet_D"
P_wet_graphing <- as.data.frame(P_wet_est$Species_table)
P_wet_graphing$Experiment <- "wet_P"
Z_wet_graphing <- as.data.frame(Z_wet_est$Species_table)
Z_wet_graphing$Experiment <- "wet_Z"

# clean up estimator names
estimators <- rownames(H25_wet_graphing)
estimators2 <- c("Homogeneous Model", "Chao2", "Chao2-bc",
                 "iChao2", "ICE", "ICE-1", "1st order jackknife",
                 "2nd order jackknife")

# clean up column names
names(H25_wet_graphing)[3] <- "Lower"
names(A_wet_graphing)[3] <- "Lower"
names(D_wet_graphing)[3] <- "Lower"
names(P_wet_graphing)[3] <- "Lower"
names(Z_wet_graphing)[3] <- "Lower"

names(H25_wet_graphing)[2] <- "SE"
names(A_wet_graphing)[2] <- "SE"
names(D_wet_graphing)[2] <- "SE"
names(P_wet_graphing)[2] <- "SE"
names(Z_wet_graphing)[2] <- "SE"

names(H25_wet_graphing)[4] <- "Upper"
names(A_wet_graphing)[4] <- "Upper"
names(D_wet_graphing)[4] <- "Upper"
names(P_wet_graphing)[4] <- "Upper"
names(Z_wet_graphing)[4] <- "Upper"


wet_season_estimators <- H25_wet_graphing %>%
  full_join(A_wet_graphing) %>%
  full_join(D_wet_graphing) %>%
  full_join(P_wet_graphing) %>%
  full_join(Z_wet_graphing) %>%
  mutate(ID = rep(estimators2, 5)) %>%
  separate(Experiment, into = c("Season", "Site"), sep = "_") %>%
  mutate(site_num = if_else(Site == "25", "1", Site)) %>%
  mutate(site_num = if_else(Site == "A", "2", site_num)) %>%
  mutate(site_num = if_else(Site == "D", "3", site_num)) %>%
  mutate(site_num = if_else(Site == "P", "4", site_num)) %>%
  mutate(site_num = if_else(Site == "Z", "5", site_num))
wet_season_estimators

## graph Figure S4c and save it ####
wet_est_graph <- ggplot(wet_season_estimators, aes(x = ID, y = Estimate, color = site_num)) +
  geom_point(stat = "identity", position = position_dodge(width = .75)) + 
  geom_errorbar(aes(ymin = wet_season_estimators$Lower, 
                    ymax = wet_season_estimators$Upper,
                    color = wet_season_estimators$site_num), position = position_dodge(width = .75)) +
  coord_flip()+
  xlab("Estimator") +
  ylab("Number of Species Estimated") +
  labs(color = "Site") +
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        #panel.border=element_blank(), #gets rid of square going around the entire graph
        #axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top")
wet_est_graph

# commenting out to prevent overwriting
#ggsave("Ruzi_etal_FigureS4c.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = wet_est_graph,
#       path = figure_path)


## H25_wet_est output ####
#(1) BASIC DATA INFORMATION:

#  Variable Value
#Number of observed species                  D    15
#Number of sampling units                    T    16
#Total number of incidences                  U    30
#Coverage estimate for entire dataset        C 0.717
#CV for entire dataset                      CV 1.257

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq    20
#Number of observed species for infrequent group   D_infreq    14
#Estimated sample coverage for infrequent group    C_infreq 0.575
#Estimated CV for infrequent group in ICE         CV_infreq     0
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq     0
#Number of observed species for frequent group       D_freq     1

#Q1 Q2 Q3 Q4
#Incidence freq. counts  9  4  1  0


#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#                           Estimate  s.e. 95%Lower 95%Upper
#Homogeneous Model           25.340 7.432   17.919   51.630
#Chao2 (Chao, 1987)          24.492 8.489   17.113   57.648
#Chao2-bc                    21.750 6.075   16.492   45.544
#iChao2 (Chiu et al. 2014)   25.968 5.145   19.576   41.287
#ICE (Lee & Chao, 1994)      25.340 7.432   17.919   51.630
#ICE-1 (Lee & Chao, 1994)    25.340 7.432   17.919   51.630
#1st order jackknife         23.438 4.043   18.461   35.569
#2nd order jackknife         28.046 6.729   20.036   48.799


## A_wet_est output ####
#(1) BASIC DATA INFORMATION:

#  Variable Value
#Number of observed species                  D    18
#Number of sampling units                    T    16
#Total number of incidences                  U    67
#Coverage estimate for entire dataset        C 0.888
#CV for entire dataset                      CV 1.179

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq    19
#Number of observed species for infrequent group   D_infreq    13
#Estimated sample coverage for infrequent group    C_infreq 0.605
#Estimated CV for infrequent group in ICE         CV_infreq     0
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq     0
#Number of observed species for frequent group       D_freq     5

#Q1 Q2 Q3 Q4
#Incidence freq. counts  8  4  1  0


#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#                           Estimate  s.e. 95%Lower 95%Upper
#Homogeneous Model           26.478 6.361   20.288   49.413
#Chao2 (Chao, 1987)          25.500 7.049   19.577   53.671
#Chao2-bc                    23.250 5.048   19.075   43.637
#iChao2 (Chiu et al. 2014)   26.773 4.411   21.460   40.244
#ICE (Lee & Chao, 1994)      26.478 6.361   20.288   49.413
#ICE-1 (Lee & Chao, 1994)    26.478 6.361   20.288   49.413
#1st order jackknife         25.500 3.812   20.931   37.191
#2nd order jackknife         29.233 6.339   22.008   49.481



## D_wet_est output ####
#(1) BASIC DATA INFORMATION:

#  Variable Value
#Number of observed species                  D    16
#Number of sampling units                    T    16
#Total number of incidences                  U    35
#Coverage estimate for entire dataset        C 0.789
#CV for entire dataset                      CV 0.895

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq    21
#Number of observed species for infrequent group   D_infreq    14
#Estimated sample coverage for infrequent group    C_infreq 0.648
#Estimated CV for infrequent group in ICE         CV_infreq     0
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq     0
#Number of observed species for frequent group       D_freq     2

#Q1 Q2 Q3 Q4
#Incidence freq. counts  8  5  1  0


#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#  Estimate  s.e. 95%Lower 95%Upper
#Homogeneous Model           23.593 5.586   18.092   43.560
#Chao2 (Chao, 1987)          22.000 5.586   17.275   44.225
#Chao2-bc                    20.375 4.245   16.886   37.592
#iChao2 (Chiu et al. 2014)   23.185 2.996   19.278   31.750
#ICE (Lee & Chao, 1994)      23.593 5.586   18.092   43.560
#ICE-1 (Lee & Chao, 1994)    23.593 5.586   18.092   43.560
#1st order jackknife         23.500 3.812   18.931   35.191
#2nd order jackknife         26.417 6.327   19.473   47.242



## P_wet_est output ####
#(1) BASIC DATA INFORMATION:

#  Variable Value
#Number of observed species                  D    14
#Number of sampling units                    T    16
#Total number of incidences                  U    43
#Coverage estimate for entire dataset        C 0.897
#CV for entire dataset                      CV 1.188

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq    22
#Number of observed species for infrequent group   D_infreq    12
#Estimated sample coverage for infrequent group    C_infreq 0.799
#Estimated CV for infrequent group in ICE         CV_infreq     0
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq     0
#Number of observed species for frequent group       D_freq     2

#Q1 Q2 Q3 Q4
#Incidence freq. counts  5  5  1  1


#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#                          Estimate  s.e. 95%Lower 95%Upper
#Homogeneous Model           17.010 2.770   14.648   27.973
#Chao2 (Chao, 1987)          16.344 2.799   14.370   28.838
#Chao2-bc                    15.562 2.091   14.214   25.383
#iChao2 (Chiu et al. 2014)   16.919 1.417   15.185   21.191
#ICE (Lee & Chao, 1994)      17.010 2.770   14.648   27.973
#ICE-1 (Lee & Chao, 1994)    17.010 2.770   14.648   27.973
#1st order jackknife         18.688 3.014   15.480   28.846
#2nd order jackknife         18.979 4.974   14.975   39.427



## Z_wet_est output ####
# (1) BASIC DATA INFORMATION:
#   
#   Variable Value
# Number of observed species                  D    21
# Number of sampling units                    T    16
# Total number of incidences                  U    48
# Coverage estimate for entire dataset        C 0.693
# CV for entire dataset                      CV 1.575
# 
# Variable Value
# Cut-off point                                            k     4
# Total number of incidences in infrequent group    U_infreq    23
# Number of observed species for infrequent group   D_infreq    18
# Estimated sample coverage for infrequent group    C_infreq 0.359
# Estimated CV for infrequent group in ICE         CV_infreq 0.831
# Estimated CV1 for infrequent group in ICE-1     CV1_infreq 1.289
# Number of observed species for frequent group       D_freq     3

# Q1 Q2 Q3 Q4
# Incidence freq. counts 15  2  0  1


# (2) SPECIES RICHNESS ESTIMATORS TABLE:

#                            Estimate   s.e. 95%Lower 95%Upper
# Homogeneous Model           53.109 23.419   29.927  136.491
# Chao2 (Chao, 1987)          73.734 46.742   32.857  255.543
# Chao2-bc                    53.812 24.084   30.062  139.806
# iChao2 (Chiu et al. 2014)   73.734 42.339   34.235  231.125
# ICE (Lee & Chao, 1994)      81.926 47.514   36.779  256.243
# ICE-1 (Lee & Chao, 1994)   122.512 97.506   41.812  516.135
# 1st order jackknife         35.062  5.220   27.954   49.437
# 2nd order jackknife         46.554  8.727   34.327   69.999


