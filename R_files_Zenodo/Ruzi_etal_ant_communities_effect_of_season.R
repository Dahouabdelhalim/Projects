## Answer the question: Does the foraging ant community differ between wet
## and dry seasons?

## -- load libraries ####
library(here)
library(tidyverse)
library(data.table)
library(vegan)
library(directlabels)

## -- set file paths ####
raw_data_path <- here::here("data/raw_data/csv")
rdata_path <- here::here("data/rdata")
figure_path <- here::here("figs")

## -- read in the data file ####
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

## -- convert data to a matrix ####
names(Ruzi_etal_freq_data_pitfalls_dry_wet_season)
# the ant species columns are 2 through 60

pitfalls_mat <- as.matrix(Ruzi_etal_freq_data_pitfalls_dry_wet_season[,2:60])
dimnames(pitfalls_mat)
pitfallsnames_mat <- as.vector(Ruzi_etal_freq_data_pitfalls_dry_wet_season$Experiment)
dimnames(pitfalls_mat)[[1]] <- c(pitfallsnames_mat)
head(pitfalls_mat)

## to check if there are any na's anywhere in the dataset
csum_mat <- colSums(pitfalls_mat)
any(is.na(csum_mat)) # FALSE

## -- make the initial NMDS ####
## will comment out to prevent overwriting in future

#pitfalls_NMDS_2dim <- metaMDS(pitfalls_mat, distance = "bray", k = 2,
#                              trymax = 50, autotransform = TRUE, pc = TRUE)

#pitfalls_NMDS_3dim <- metaMDS(pitfalls_mat, distance = "bray", k = 3, 
#                              trymax = 50, autotransform = TRUE, pc = TRUE)

#pitfalls_NMDS_3dim$stress #0.0396511 
#pitfalls_NMDS_2dim$stress #0.07502392
## adding another dimension did not reduce stress by more than 0.05 
## keeping the 2 dim solution

## -- save the metaMDS output ####
#save(pitfalls_NMDS_2dim, 
#     file = paste(rdata_path, "drypitfalls_wetpitfalls_NMDS_2dim.rda", sep = "/"))

## -- load the saved Rdata metaMDS output ####
load(paste(rdata_path,"drypitfalls_wetpitfalls_NMDS_2dim.rda", sep = "/"))

## -- conduct MRPP ####
stressplot(pitfalls_NMDS_2dim)
pitfalls.dist <- vegdist(pitfalls_mat, method = "bray")

pitfalls.mrpp <- mrpp(dat = pitfalls.dist, grouping = Ruzi_etal_freq_data_pitfalls_dry_wet_season$Experiment,
                      permutations = 999, distance = "bray", 
                      strata = Ruzi_etal_freq_data_pitfalls_dry_wet_season$Plot_abbre)
pitfalls.mrpp
pitfalls.mrpp$Pvalue # 0.0625
pitfalls.mrpp$distance # "bray"
pitfalls.mrpp$A # 0.003412522

## -- make the NMDS figure ####
pitfalls_NMDS1 <- pitfalls_NMDS_2dim$points[,1]
pitfalls_NMDS2 <- pitfalls_NMDS_2dim$points[,2]
pitfalls_species1 <- pitfalls_NMDS_2dim$species[,1]
pitfalls_species2 <- pitfalls_NMDS_2dim$species[,2]
pitfalls_graphing_file <- cbind(Ruzi_etal_freq_data_pitfalls_dry_wet_season, pitfalls_NMDS1, pitfalls_NMDS2)

pitfalls_species_graphing <- cbind(pitfalls_species1,pitfalls_species2)
pitfalls_plot_level <- as.data.table(dimnames(pitfalls_species_graphing)[[1]])

pitfalls_graph <- ggplot(pitfalls_graphing_file, 
                         aes(x=pitfalls_NMDS1,y=pitfalls_NMDS2, 
                             fill=Experiment, colour=Experiment))+ 
  geom_point()+ #position=position_jitter(.1), shape=3
  #geom_text(data=pitfalls_graphing_file, 
  #          aes(x=pitfalls_NMDS1,y=pitfalls_NMDS2, 
  #              label=Plot_abbre))+
  #geom_text(data=plot_level_species_graphing,
  #          aes(x=plot_level_species1,y=plot_level_species2,
  #              label=species_names_plot_level))+
  xlab("NMDS 1")+
  ylab("NMDS 2")+
  annotate("text", x=min(pitfalls_NMDS1)-.05, y=min(pitfalls_NMDS2)-.5,
           label=paste('Stress =',round(pitfalls_NMDS_2dim$stress,3)))+
  stat_ellipse(type='t',size=1, geom="polygon", alpha=.2)+ #t assumes multivariet t distribution
  scale_fill_manual(values=c("red","blue"),
                    breaks=c("pitfall_wet","pitfall_dry"),
                    labels=c("Pitfalls\\nWet Season","Pitfalls\\nDry Season"))+
  scale_colour_manual(values=c("red","blue"),
                      breaks=c("pitfall_wet","pitfall_dry"),
                      labels=c("Pitfalls\\nWet Season","Pitfalls\\nDry Season"))+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        legend.position="top",
        legend.title=element_blank(),
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"))
pitfalls_graph 

## -- save figure output ####

# commenting out so as to now overwrite
#ggsave("Ruzi_etal_FigureS5.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = pitfalls_graph,
#       path = figure_path) 


