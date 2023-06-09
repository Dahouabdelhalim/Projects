## to process raw data for pitfall - aboveground species accumulation curves for the dry season,
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

t2_dry <- t2 %>%
  filter(season == "dry") %>%
  select(-season)
head(t2_dry)

## now to transpose the dry season tibble
class(t2_dry) 
t2_dry <- as_tibble(t2_dry)

t3_dry <- t2_dry %>%
  gather(spp_names, key = "spp", value = "presence", factor_key = TRUE) %>%
  spread(key = plot, value = presence)
t3_dry

### now need to separate each one of the sites into it's own list,
### and put the trap numer, which was 16 for all of them at the begining of the vector
## and can remove the 0's
num_traps <- 16

## dry season -- 25Ha
d25Ha <- 
  t3_dry %>%
  filter(.$"25" > 0) %>%
  select("25")
d25Ha


d_25Ha <- c(num_traps, d25Ha$`25`) # add in the sample number
d_25Ha

## dry season -- AVA
dA <- 
  t3_dry %>%
  filter(A > 0) %>%
  select(A)
dA


d_A <- c(num_traps, dA$A) # add in the sample number
d_A

## dry season -- Drayton
dD <- 
  t3_dry %>%
  filter(D > 0) %>%
  select(D)
dD


d_D <- c(num_traps, dD$D) # add in the sample number
d_D

## dry season -- Pearson
dP <- 
  t3_dry %>%
  filter(P > 0) %>%
  select(P)
dP


d_P <- c(num_traps, dP$P) # add in the sample number
d_P

## dry season -- Zetek
dZ <- 
  t3_dry %>%
  filter(Z > 0) %>%
  select(Z)
dZ


d_Z <- c(num_traps, dZ$Z) # add in the sample number
d_Z

## dry season list 
# changing from site names to site numbers
# 25HA = 1
# Ava = 2
# Drayton = 3
# Pearson = 4
# Zetek = 5

dry_list <- list("1" = d_25Ha, "2" = d_A, "3" = d_D,
                 "4" = d_P, "5" = d_Z)
dry_list

## dry season rarefaction and extrapolation ####
# choosing to extrapolate out to double the sample size
t_pit <- seq(1,32) 
t_pit

dry_out <- iNEXT(dry_list, q = c(0), datatype = "incidence_freq",
                 size = t_pit) 
dry_out

## graph the iNEXT results - Figure S3c,d ####
pit_dry_iNEXT_graph2 <- ggiNEXT(dry_out, type = 1) +
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
        legend.position="top",
        legend.key=element_blank())
pit_dry_iNEXT_graph2 # to get the legend

pit_dry_iNEXT_graph <- ggiNEXT(dry_out, type = 1) +
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
pit_dry_iNEXT_graph # without legend

pit_dry_CV_iNEXT_graph <- ggiNEXT(dry_out, type = 2) +
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,1, by=.1)) +
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
        #legend.title=element_blank(),
        #legend.text=element_text(face="bold"),
        #legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank())
pit_dry_CV_iNEXT_graph

## save the species accumulation curves and sampling coverage outputs ####
# commented out to not overwrite

#ggsave("Ruzi_etal_FigureS3c_legend.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = pit_dry_iNEXT_graph2,
#       path = figure_path)
#ggsave("Ruzi_etal_FigureS3c_no_legend.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = pit_dry_iNEXT_graph,
#       path = figure_path)
#ggsave("Ruzi_etal_FigureS3d.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = pit_dry_CV_iNEXT_graph,
#       path = figure_path)

## to get species estimates using SpadeR ####
dry_list

## choose k to be the cut off where "rare" is when found in less than 25% of the traps
## and common is then more than that

H25_dry_est <- SpadeR::ChaoSpecies(dry_list$'1', datatype = "incidence_freq", k = 4, conf = 0.95)

A_dry_est <- SpadeR::ChaoSpecies(dry_list$'2', datatype = "incidence_freq", k = 4, conf = 0.95)

D_dry_est <- SpadeR::ChaoSpecies(dry_list$'3', datatype = "incidence_freq", k = 4, conf = 0.95)

P_dry_est <- SpadeR::ChaoSpecies(dry_list$'4', datatype = "incidence_freq", k = 4, conf = 0.95)

Z_dry_est <- SpadeR::ChaoSpecies(dry_list$'5', datatype = "incidence_freq" , k = 4, conf = 0.95)


H25_dry_graphing <- as.data.frame(H25_dry_est$Species_table)
H25_dry_graphing$Experiment <- "dry_25"
A_dry_graphing <- as.data.frame(A_dry_est$Species_table)
A_dry_graphing$Experiment <- "dry_A"
D_dry_graphing <- as.data.frame(D_dry_est$Species_table)
D_dry_graphing$Experiment <- "dry_D"
P_dry_graphing <- as.data.frame(P_dry_est$Species_table)
P_dry_graphing$Experiment <- "dry_P"
Z_dry_graphing <- as.data.frame(Z_dry_est$Species_table)
Z_dry_graphing$Experiment <- "dry_Z"

# clean up estimator names
estimators <- rownames(H25_dry_graphing)
estimators2 <- c("Homogeneous Model", "Chao2", "Chao2-bc",
                 "iChao2", "ICE", "ICE-1", "1st order jackknife",
                 "2nd order jackknife")

# clean up column names
names(H25_dry_graphing)[3] <- "Lower"
names(A_dry_graphing)[3] <- "Lower"
names(D_dry_graphing)[3] <- "Lower"
names(P_dry_graphing)[3] <- "Lower"
names(Z_dry_graphing)[3] <- "Lower"

names(H25_dry_graphing)[2] <- "SE"
names(A_dry_graphing)[2] <- "SE"
names(D_dry_graphing)[2] <- "SE"
names(P_dry_graphing)[2] <- "SE"
names(Z_dry_graphing)[2] <- "SE"

names(H25_dry_graphing)[4] <- "Upper"
names(A_dry_graphing)[4] <- "Upper"
names(D_dry_graphing)[4] <- "Upper"
names(P_dry_graphing)[4] <- "Upper"
names(Z_dry_graphing)[4] <- "Upper"


dry_season_estimators <- H25_dry_graphing %>%
  full_join(A_dry_graphing) %>%
  full_join(D_dry_graphing) %>%
  full_join(P_dry_graphing) %>%
  full_join(Z_dry_graphing) %>%
  mutate(ID = rep(estimators2, 5)) %>%
  separate(Experiment, into = c("Season", "Site"), sep = "_") %>%
  mutate(site_num = if_else(Site == "25", "1", Site)) %>%
  mutate(site_num = if_else(Site == "A", "2", site_num)) %>%
  mutate(site_num = if_else(Site == "D", "3", site_num)) %>%
  mutate(site_num = if_else(Site == "P", "4", site_num)) %>%
  mutate(site_num = if_else(Site == "Z", "5", site_num))
dry_season_estimators

## graph Figure S4b and save it ####
dry_est_graph <- ggplot(dry_season_estimators, aes(x = ID, y = Estimate, color = site_num)) +
  geom_point(stat = "identity", position = position_dodge(width = .75)) + 
  geom_errorbar(aes(ymin = dry_season_estimators$Lower, 
                    ymax = dry_season_estimators$Upper,
                    color = dry_season_estimators$site_num), position = position_dodge(width = .75)) +
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
dry_est_graph

# commenting out to prevent overwriting
#ggsave("Ruzi_etal_FigureS4b.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = dry_est_graph,
#       path = figure_path)


## H25_dry_est output ####

#(1) BASIC DATA INFORMATION:
#  
#  Variable Value
#Number of observed species                  D    18
#Number of sampling units                    T    16
#Total number of incidences                  U    29
#Coverage estimate for entire dataset        C 0.569
#CV for entire dataset                      CV 1.251

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq    21
#Number of observed species for infrequent group   D_infreq    17
#Estimated sample coverage for infrequent group    C_infreq 0.405
#Estimated CV for infrequent group in ICE         CV_infreq     0
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq     0
#Number of observed species for frequent group       D_freq     1

#Q1 Q2 Q3 Q4
#Incidence freq. counts 13  4  0  0


#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#                           Estimate   s.e. 95%Lower 95%Upper
#Homogeneous Model           42.939 16.739   25.547  100.411
#Chao2 (Chao, 1987)          37.805 15.445   23.129   94.469
#Chao2-bc                    32.625 10.966   21.950   72.154
#iChao2 (Chiu et al. 2014)   37.805 11.197   25.055   73.596
#ICE (Lee & Chao, 1994)      42.939 16.739   25.547  100.411
#ICE-1 (Lee & Chao, 1994)    42.939 16.739   25.547  100.411
#1st order jackknife         30.188  4.859   23.741   43.874
#2nd order jackknife         38.296  8.104   27.551   61.131


## A_dry_est output ####
#(1) BASIC DATA INFORMATION:
#  
#  Variable Value
#Number of observed species                  D    20
#Number of sampling units                    T    16
#Total number of incidences                  U    57
#Coverage estimate for entire dataset        C 0.796
#CV for entire dataset                      CV 1.299

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq    21
#Number of observed species for infrequent group   D_infreq    16
#Estimated sample coverage for infrequent group    C_infreq 0.447
#Estimated CV for infrequent group in ICE         CV_infreq 0.301
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq 0.406
#Number of observed species for frequent group       D_freq     4
#Q1 Q2 Q3 Q4
#Incidence freq. counts 12  3  1  0
#(2) SPECIES RICHNESS ESTIMATORS TABLE:
#  Estimate   s.e. 95%Lower 95%Upper
#Homogeneous Model           39.794 13.604   25.849   86.985
#Chao2 (Chao, 1987)          42.500 18.974   25.350  114.630
#Chao2-bc                    35.469 12.145   23.977   80.164
#iChao2 (Chiu et al. 2014)   44.673 14.850   28.298   93.366
#ICE (Lee & Chao, 1994)      42.233 17.855   25.578  108.615
#ICE-1 (Lee & Chao, 1994)    44.229 22.660   25.123  134.581
#1st order jackknife         31.250  4.669   25.150   44.573
#2nd order jackknife         39.300  7.792   29.011   61.337
## D_dry_est output ####
#(1) BASIC DATA INFORMATION:
#  
#  Variable Value
#Number of observed species                  D    18
#Number of sampling units                    T    16
#Total number of incidences                  U    43
#Coverage estimate for entire dataset        C 0.724
#CV for entire dataset                      CV 1.059

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq    18
#Number of observed species for infrequent group   D_infreq    14
#Estimated sample coverage for infrequent group    C_infreq 0.341
#Estimated CV for infrequent group in ICE         CV_infreq 0.846
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq  1.33
#Number of observed species for frequent group       D_freq     4

#Q1 Q2 Q3 Q4
#Incidence freq. counts 12  0  2  0
#(2) SPECIES RICHNESS ESTIMATORS TABLE:
#                           Estimate   s.e. 95%Lower 95%Upper
#Homogeneous Model           45.018 22.253   24.594  128.706
#Chao2 (Chao, 1987)          79.875 37.412   38.719  202.780
#Chao2-bc                    79.875 38.167   38.332  206.303
#iChao2 (Chiu et al. 2014)   84.750 35.853   42.883  197.061
#ICE (Lee & Chao, 1994)      70.182 46.660   29.616  252.411
#ICE-1 (Lee & Chao, 1994)   107.191 91.495   34.927  487.974
#1st order jackknife         29.250  4.669   23.150   42.573
#2nd order jackknife         39.750  7.821   28.981   61.082
## P_dry_est output ####
#(1) BASIC DATA INFORMATION:
#  
#  Variable Value
#Number of observed species                  D    10
#Number of sampling units                    T    16
#Total number of incidences                  U    38
#Coverage estimate for entire dataset        C 0.875
#CV for entire dataset                      CV 1.304

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq     9
#Number of observed species for infrequent group   D_infreq     7
#Estimated sample coverage for infrequent group    C_infreq 0.473
#Estimated CV for infrequent group in ICE         CV_infreq     0
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq     0
#Number of observed species for frequent group       D_freq     3

#Q1 Q2 Q3 Q4
#Incidence freq. counts  5  2  0  0
#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#  Estimate  s.e. 95%Lower 95%Upper
#Homogeneous Model           17.812 7.964   11.494   50.860
#Chao2 (Chao, 1987)          15.859 7.106   10.908   47.791
#Chao2-bc                    13.125 3.897   10.468   30.855
#iChao2 (Chiu et al. 2014)   15.859 5.015   11.370   35.054
#ICE (Lee & Chao, 1994)      17.812 7.964   11.494   50.860
#ICE-1 (Lee & Chao, 1994)    17.812 7.964   11.494   50.860
#1st order jackknife         14.688 3.014   11.480   24.846
#2nd order jackknife         17.429 5.019   12.234   34.709
## Z_dry_est output ####
#(1) BASIC DATA INFORMATION:

#  Variable Value
#Number of observed species                  D    14
#Number of sampling units                    T    16
#Total number of incidences                  U    50
#Coverage estimate for entire dataset        C 0.907
#CV for entire dataset                      CV 0.851

#Variable Value
#Cut-off point                                            k     4
#Total number of incidences in infrequent group    U_infreq    15
#Number of observed species for infrequent group   D_infreq     9
#Estimated sample coverage for infrequent group    C_infreq 0.691
#Estimated CV for infrequent group in ICE         CV_infreq 0.436
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq 0.562
#Number of observed species for frequent group       D_freq     5

#Q1 Q2 Q3 Q4
#Incidence freq. counts  5  3  0  1


#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#                           Estimate  s.e. 95%Lower 95%Upper
#Homogeneous Model           18.018 3.801   14.838   33.256
#Chao2 (Chao, 1987)          17.906 4.604   14.627   38.332
#Chao2-bc                    16.344 2.989   14.342   30.084
#iChao2 (Chiu et al. 2014)   17.906 3.102   14.992   29.382
#ICE (Lee & Chao, 1994)      19.393 5.812   14.965   44.141
#ICE-1 (Lee & Chao, 1994)    20.303 7.783   14.956   55.563
#1st order jackknife         18.688 3.014   15.480   28.846
#2nd order jackknife         20.613 5.004   15.769   38.724