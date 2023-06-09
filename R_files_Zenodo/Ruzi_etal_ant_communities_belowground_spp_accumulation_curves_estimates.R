## to process raw data for belowground species accumulation curves, sampling coverage,
## and species estimates

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

## -- load datafile as a tibble ####

Ruzi_etal_raw_data_belowground <- read_csv(paste(raw_data_path, "Ruzi_etal_raw_data_belowground.csv", sep = "/"),
                    col_types = cols(
                      .default = col_double(),
                      Sample = col_character(),
                      Plot = col_character(),
                      Experiment = col_character(),
                      Replicate = col_character(),
                      Season = col_character(),
                      Position = col_character(),
                      Species = col_character(),
                      Data_type = col_character()
                    ))
Ruzi_etal_raw_data_belowground

names(Ruzi_etal_raw_data_belowground)

# these columns contain informatin about the traps and do not pertain to ant species
ID_cols <- c("Sample", "Plot", "Experiment", "Replicate", "Year", "Season", "Position",
             "Species", "Sample_num", "Data_type")

# these are the columns that contain the ant species
spp_cols <- names(Ruzi_etal_raw_data_belowground[11:113])

# filter out the control traps from seed associated traps ####
bel_controls <- Ruzi_etal_raw_data_belowground %>%
  filter(Species == "glass" | Species == "passive")
bel_controls

bel_seeds <- Ruzi_etal_raw_data_belowground %>%
  filter(Species != "glass" & Species != "passive")
bel_seeds


## getting the incidence frequences of the belowground samples ####

# doing controls separately from seed caches
# setting rem_spp to TRUE to remove ant species not collected at all
# settign freq to TRUE to output frequencies 

b_cont_inc_freq <- SBS_transform(bel_controls, trap_ID = ID_cols,
                                 spp_ID = spp_cols, rem_spp = TRUE, 
                                 freq = TRUE, warn = FALSE)
b_cont_inc_freq # the first number is the number of traps that are associated with controls (silica beads and passive/unbaited traps)


b_seeds_inc_freq <- SBS_transform(bel_seeds, trap_ID = ID_cols,
                                  spp_ID = spp_cols, rem_spp = TRUE, 
                                  freq = TRUE, warn = FALSE)
b_seeds_inc_freq # the first number is the number of traps that are traps associated with seed caches

## species accumulation curves ####

# join into a list
below_list <- list(cont = b_cont_inc_freq, seeds = b_seeds_inc_freq)

# set how long want to extrapolate out (going to 120 which is the number of traps associated with seed caches)
t <- seq(1,120)

# run iNEXT analysis
below_out <- iNEXT(below_list, q = c(0), datatype = "incidence_freq", size = t)

# graph the species accumulation curve -- Figure S3a
bel_iNEXT_graph <- ggiNEXT(below_out, type = 1) +
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,35, by=5)) +
  annotate("text", label="Associated with\\nSeed Traps:", x=20, y=25,
           size=5, colour="black") + # family = "Times"
  annotate("text", label="Yes", x=20, y=23, size=4, colour="black") + # family = "Times"
  annotate("segment", x=25, xend=30, y=23, yend=23, colour="black") + # family = "Times"
  annotate("text", label="No", x=20, y=21, size=4, colour="red") + # family = "Times"
  annotate("segment", x=25, xend=30, y=21, yend=21, colour="red") + # family = "Times"
  scale_x_continuous(expand= c(0,0), breaks=seq(0,130, by=10)) + # expand = c(0,0)
  scale_color_manual(values = c("red","black")) +
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
        legend.key=element_blank(), #remove the border around each item
        #axis.title.x = element_text(family="Arial", face="bold", size=14, color="black"), #size of x-axis title
        #axis.title.y = element_text(family="Arial", face="bold", size=14, color="black"), #size of y-axis title
        #axis.text.x = element_text(family="Arial", size=12, color="black"), #size of x-axis text
        #axis.text.y = element_text(family="Arial", size=12, color="black")
  )
bel_iNEXT_graph  

# graph sample coverage -- Figure S3b
bel_iNEXT_SC_graph <- ggiNEXT(below_out, type = 2) +
  annotate("text", label="Associated with\\nSeed Traps:", x=90, y=.4,
           size=5, colour="black") + # family = "Times"
  annotate("text", label="Yes", x=90, y=.3, size=4, colour="black") + # family = "Times"
  annotate("segment", x=95, xend=100, y=.3, yend=.3, colour="black") + # family = "Times"
  annotate("text", label="No", x=90, y=.25, size=4, colour="red") + # family = "Times"
  annotate("segment", x=95, xend=100, y=.25, yend=.25, colour="red") + # family = "Times"
  scale_y_continuous(expand = c(0, 0), breaks=seq(0,1, by=.1)) +
  scale_x_continuous(expand= c(0,0), breaks=seq(0,130, by=10)) + # expand = c(0,0)
  scale_color_manual(values = c("red","black")) +
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
        legend.key=element_blank(), #remove the border around each item
        #axis.title.x = element_text(family="Arial", face="bold", size=14, color="black"), #size of x-axis title
        #axis.title.y = element_text(family="Arial", face="bold", size=14, color="black"), #size of y-axis title
        #axis.text.x = element_text(family="Arial", size=12, color="black"), #size of x-axis text
        #axis.text.y = element_text(family="Arial", size=12, color="black")
  )
bel_iNEXT_SC_graph

## save Figures S3a,b to the fig folder ####
# commented out so as to not overwrite

#ggsave("Ruzi_etal_FigureS3a.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = bel_iNEXT_graph,
#       path = figure_path)
#ggsave("Ruzi_etal_FigureS3b.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = bel_iNEXT_SC_graph,
#       path = figure_path)


## to get species estimates using SpadeR ####
below_list

## using the default k for below-ground, keeping that k the same so can compare the two
## below_control amount (40 traps total so cut off is k = 10)
## below_seeds amount (120 traps total so cut off is k = 10)

bel_cont_est <- SpadeR::ChaoSpecies(below_list$cont, datatype = "incidence_freq", k = 10, conf = 0.95)

bel_seeds_est <- SpadeR::ChaoSpecies(below_list$seeds, datatype = "incidence_freq", k = 10, conf = 0.95)

## to graph the information from SpadeR ####
bel_cont_graphing <- as.data.frame(bel_cont_est$Species_table)
bel_cont_graphing$Experiment <- "bel_cont"
bel_seeds_graphing <- as.data.frame(bel_seeds_est$Species_table)
bel_seeds_graphing$Experiment <- "bel_seeds"

# to get the estimates
estimators <- rownames(bel_cont_graphing)
# clean up the estimates in the same order as estimators was
estimators2 <- c("Homogeneous Model", "Chao2", "Chao2-bc",
                 "iChao2", "ICE", "ICE-1", "1st order jackknife",
                 "2nd order jackknife")

# clean up column names
names(bel_cont_graphing)[3] <- "Lower"
names(bel_cont_graphing)[2] <- "SE"
names(bel_cont_graphing)[4] <- "Upper"

names(bel_seeds_graphing)[3] <- "Lower"
names(bel_seeds_graphing)[2] <- "SE"
names(bel_seeds_graphing)[4] <- "Upper"

below_graphing_estimates <- bel_cont_graphing %>%
  # join the two datasets together
  full_join(bel_seeds_graphing) %>%
  mutate(ID = rep(estimators2,2)) %>%
  # set up color
  mutate(line_col = if_else(Experiment == "bel_cont", "red", Experiment)) %>%
  mutate(line_col = if_else(line_col != "red", "black", line_col))
below_graphing_estimates


below_est_graph <- ggplot(below_graphing_estimates, aes(x = ID, y = Estimate, color = Experiment)) +
  geom_point(stat = "identity", position = position_dodge(width = .75)) + 
  geom_errorbar(aes(ymin = below_graphing_estimates$Lower, 
                    ymax = below_graphing_estimates$Upper,
                    color = below_graphing_estimates$Experiment), position = position_dodge(width = .75)) +
  coord_flip()+
  xlab("Estimator") +
  ylab("Number of Species Estimated") +
  scale_color_manual(values = c("red","black"),
                     name = "Associated with\\nSeed Traps",
                     labels = c("no", "yes")) +
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        #panel.border=element_blank(), #gets rid of square going around the entire graph
        #axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top")
below_est_graph


## save Figures S4a to the fig folder ####
# commented out so as to not overwrite

#ggsave("Ruzi_etal_FigureS4a.png", width = 8, height = 6,
#       units = "in", dpi = 300, plot = below_est_graph,
#       path = figure_path)

# bel_cont_est output ####
#(1) BASIC DATA INFORMATION:

#  Variable Value
#Number of observed species                  D    10
#Number of sampling units                    T    40
#Total number of incidences                  U    35
#Coverage estimate for entire dataset        C  0.83
#CV for entire dataset                      CV 0.957

#Variable Value
#Cut-off point                                            k    10
#Total number of incidences in infrequent group    U_infreq    35
#Number of observed species for infrequent group   D_infreq    10
#Estimated sample coverage for infrequent group    C_infreq  0.83
#Estimated CV for infrequent group in ICE         CV_infreq 0.986
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq 1.457
#Number of observed species for frequent group       D_freq     0

#Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 Q10
#Incidence freq. counts  6  0  0  0  1  0  1  1  1   0


#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#  Estimate   s.e. 95%Lower 95%Upper
#Homogeneous Model           12.044  2.067   10.394   20.592
#Chao2 (Chao, 1987)          24.625 13.216   13.217   76.490
#Chao2-bc                    24.625 13.681   13.092   79.180
#iChao2 (Chiu et al. 2014)   24.625 11.424   13.781   66.566
#ICE (Lee & Chao, 1994)      19.069  8.958   11.800   55.703
#ICE-1 (Lee & Chao, 1994)    27.383 21.409   12.645  124.263
#1st order jackknife         15.850  3.399   12.033   26.836
#2nd order jackknife         21.550  5.812   14.552   39.309


#(3) DESCRIPTION OF ESTIMATORS/MODELS:

#  Homogeneous Model: This model assumes that all species have the same incidence or detection probabilities. See Eq. (3.2) of Lee and Chao (1994) or Eq. (12a) in Chao and Chiu (2016b).

#Chao2 (Chao, 1987): This approach uses the frequencies of uniques and duplicates to estimate the number of undetected species; see Chao (1987) or Eq. (11a) in Chao and Chiu (2016b).

#Chao2-bc: A bias-corrected form for the Chao2 estimator; see Chao (2005).

#iChao2: An improved Chao2 estimator; see Chiu et al. (2014).

#ICE (Incidence-based Coverage Estimator): A non-parametric estimator originally proposed by Lee and Chao (1994) in the context of capture-recapture data analysis. The observed species are separated as frequent and infrequent species groups; only data in the infrequent group are used to estimate the number of undetected species. The estimated CV for species in the infrequent group characterizes the degree of heterogeneity among species incidence probabilities. See Eq. (12b) of Chao and Chiu (2016b), which is an improved version of Eq. (3.18) in Lee and Chao (1994). This model is also called Model(h) in capture-recapture literature where h denotes "heterogeneity".

#ICE-1: A modified ICE for highly-heterogeneous cases.

#1st order jackknife: It uses the frequency of uniques to estimate the number of undetected species; see Burnham and Overton (1978).

#2nd order jackknife: It uses the frequencies of uniques and duplicates to estimate the number of undetected species; see Burnham and Overton (1978).

#95% Confidence interval: A log-transformation is used for all estimators so that the lower bound of the resulting interval is at least the number of observed species. See Chao (1987).



## to get Figure S4a

# bel_seeds_est output ####
#(1) BASIC DATA INFORMATION:

#  Variable Value
#Number of observed species                  D    18
#Number of sampling units                    T   120
#Total number of incidences                  U   114
#Coverage estimate for entire dataset        C  0.93
#CV for entire dataset                      CV 1.402

#Variable Value
#Cut-off point                                            k    10
#Total number of incidences in infrequent group    U_infreq    27
#Number of observed species for infrequent group   D_infreq    14
#Estimated sample coverage for infrequent group    C_infreq 0.706
#Estimated CV for infrequent group in ICE         CV_infreq 0.694
#Estimated CV1 for infrequent group in ICE-1     CV1_infreq 0.944
#Number of observed species for frequent group       D_freq     4

#Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 Q10
#Incidence freq. counts  8  3  1  0  2  0  0  0  0   0


#(2) SPECIES RICHNESS ESTIMATORS TABLE:

#  Estimate   s.e. 95%Lower 95%Upper
#Homogeneous Model           23.842  4.398   19.571   39.724
#Chao2 (Chao, 1987)          28.578 10.189   20.161   69.765
#Chao2-bc                    24.942  6.606   19.439   51.493
#iChao2 (Chiu et al. 2014)   30.168  7.217   22.149   53.690
#ICE (Lee & Chao, 1994)      29.309  9.072   20.841   63.018
#ICE-1 (Lee & Chao, 1994)    33.940 14.915   21.368   93.427
#1st order jackknife         25.933  3.975   21.138   38.059
#2nd order jackknife         30.875  6.851   22.839   52.257
