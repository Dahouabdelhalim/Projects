## Do ant communities differ by microhabitat?

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

## -- read in the pitfall data file ####
Ruzi_etal_freq_data_pitfalls_wet_season_below <- read_csv(paste(raw_data_path, "Ruzi_etal_freq_data_pitfalls_wet_season_below.csv", sep = "/"),
                                                        col_types = cols(
                                                          .default = col_double(),
                                                          ID = col_character(),
                                                          Plot_abbre = col_character(),
                                                          Season = col_character(),
                                                          Species = col_character(),
                                                          Experiment = col_character()
                                                        ))
Ruzi_etal_freq_data_pitfalls_wet_season_below

## -- convert data to a matrix ####
names(Ruzi_etal_freq_data_pitfalls_wet_season_below) # than species are colums 2:49

raw_data_mat <- as.matrix(Ruzi_etal_freq_data_pitfalls_wet_season_below[,2:49])
dimnames(raw_data_mat)
rawnames_mat <- as.vector(Ruzi_etal_freq_data_pitfalls_wet_season_below$ID)
dimnames(raw_data_mat)[[1]] <- c(rawnames_mat)
head(raw_data_mat)

## to check if there are any na's anywhere in the dataset
csum_mat <- colSums(raw_data_mat)
any(is.na(csum_mat)) # FALSE

## -- make the initial NMDS ####
## will comment out to prevent overwriting in future

#raw_data_NMDS_2dim <- metaMDS(raw_data_mat, distance="bray", k=2, trymax=50, autotransform=TRUE, pc=TRUE)
#stressplot(raw_data_NMDS_2dim)

#raw_data_NMDS_3dim <- metaMDS(raw_data_mat, distance="bray", k=3, trymax=100, 
#                              autotransform=TRUE)

raw_data_NMDS_4dim <- metaMDS(raw_data_mat, distance="bray", k=4, trymax=100, 
                              autotransform=TRUE)
stressplot(raw_data_NMDS_4dim)


#raw_data_NMDS_2dim$stress # 0.1942893
#raw_data_NMDS_3dim$stress # 0.1221308
#raw_data_NMDS_4dim$stress # 0.07262162

# since 3dim reduces stress by over 0.05 will save the 3dim one, the 4dim doesn't quite reduce stress enough

## -- save the metaMDS output ####
#save(raw_data_NMDS_3dim, 
#     file = paste(rdata_path, "wet_pitfalls_vs_below_NMDS_3dim.rda", sep = "/"))

## -- load the saved Rdata metaMDS output ####
load(paste(rdata_path,"wet_pitfalls_vs_below_NMDS_3dim.rda", sep = "/"))
stressplot(raw_data_NMDS_3dim)
raw_data_NMDS_3dim$stress # 0.1221308

## -- conduct MRPP ####
raw.dist <- vegdist(raw_data_mat, method = "bray")

wetpitfall_below.mrpp <- mrpp(dat=raw.dist, grouping=Ruzi_etal_freq_data_pitfalls_wet_season_below$Experiment, permutations=999,
                              distance="bray", strata=Ruzi_etal_freq_data_pitfalls_wet_season_below$Plot_abbre)
wetpitfall_below.mrpp 
wetpitfall_below.mrpp$Pvalue # 0.001
wetpitfall_below.mrpp$distance #"bray"
wetpitfall_below.mrpp$A # 0.08343449

## -- make the NMDS figures ####

raw_data_NMDS1 <- raw_data_NMDS_3dim$points[,1]
raw_data_NMDS2 <- raw_data_NMDS_3dim$points[,2]
raw_data_NMDS3 <- raw_data_NMDS_3dim$points[,3]
raw_data_species1 <- raw_data_NMDS_3dim$species[,1]
raw_data_species2 <- raw_data_NMDS_3dim$species[,2]
raw_data_species3 <- raw_data_NMDS_3dim$species[,3]
raw_data_graphing_file <- cbind(Ruzi_etal_freq_data_pitfalls_wet_season_below, raw_data_NMDS1, raw_data_NMDS2, raw_data_NMDS3)

NMDS1_NMDS2 <- ggplot(raw_data_graphing_file, 
                      aes(x=raw_data_NMDS1,y=raw_data_NMDS2, fill=Experiment,
                          colour=Experiment))+ # if go by plot because have two trials per season then get ellipses
  geom_point()+ #position=position_jitter(.1),shape=3
  #geom_text(data=raw_data_graphing_file, 
  #          aes(x=raw_data_NMDS1,y=raw_data_NMDS2, 
  #              label=Plot_abbre))+
  xlab("NMDS 1")+
  ylab("NMDS 2")+
  annotate("text", x=min(raw_data_NMDS1)+0.25, y=min(raw_data_NMDS2)-.5,
           label=paste('Stress =',round(raw_data_NMDS_3dim$stress,3)))+
  stat_ellipse(type='t',size=1, geom="polygon", alpha=.2)+ #t assumes multivariet t distribution
  scale_fill_manual(values=c("brown","blue"),
                    breaks=c("below","pitfall_wet"),
                    labels=c("Below-ground\\nWet Season","Pitfalls\\nWet Season"))+
  scale_colour_manual(values=c("brown","blue"),
                      breaks=c("below","pitfall_wet"),
                      labels=c("Below-ground\\nWet Season","Pitfalls\\nWet Season"))+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        #panel.border=element_blank(), #gets rid of square going around the entire graph
        #axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top",
        legend.title=element_blank(),
        #legend.text=element_text(face="bold"),
        #legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title # family = "Arial"
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title # family = "Arial"
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text # family = "Arial"
        axis.text.y = element_text(size=12, color="black")) # family = "Arial"
NMDS1_NMDS2

NMDS1_NMDS3 <- ggplot(raw_data_graphing_file, 
                      aes(x=raw_data_NMDS1,y=raw_data_NMDS3, 
                          fill=Experiment, colour=Experiment))+ # if go by plot because have two trials per season then get ellipses
  geom_point()+ #position=position_jitter(.1),shape=3
  scale_fill_manual(values=c("brown","blue"),
                    breaks=c("below","pitfall_wet"),
                    labels=c("Below-ground\\nWet Season","Pitfalls\\nWet Season"))+
  #geom_text(data=raw_data_graphing_file, 
  #          aes(x=raw_data_NMDS1,y=raw_data_NMDS3, 
  #              label=Plot_abbre))+
  scale_colour_manual(values=c("brown","blue"),
                      breaks=c("below","pitfall_wet"),
                      labels=c("Below-ground\\nWet Season","Pitfalls\\nWet Season"))+
  #geom_text(data=plot_level_species_graphing,
  #          aes(x=plot_level_species1,y=plot_level_species2,
  #              label=species_names_plot_level))+
  xlab("NMDS 1")+
  ylab("NMDS 3")+
  annotate("text", x=min(raw_data_NMDS1)+.7, y=min(raw_data_NMDS3)-.5,
           label=paste('Stress =',round(raw_data_NMDS_3dim$stress,3)))+
  stat_ellipse(type='t',size=1, geom="polygon", alpha=.2)+ #t assumes multivariet t distribution
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        #panel.border=element_blank(), #gets rid of square going around the entire graph
        #axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top",
        legend.title=element_blank(),
        #legend.text=element_text(face="bold"),
        #legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title # family = "Arial"
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title # family = "Arial"
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text # family = "Arial"
        axis.text.y = element_text(size=12, color="black")) # family = "Arial"
NMDS1_NMDS3

NMDS3_NMDS2 <- ggplot(raw_data_graphing_file, 
                      aes(x=raw_data_NMDS3,y=raw_data_NMDS2, 
                          fill=Experiment, colour=Experiment))+ # if go by plot because have two trials per season then get ellipses
  geom_point()+ #position=position_jitter(.1),shape=3
  scale_fill_manual(values=c("brown","blue"),
                    breaks=c("below","pitfall_wet"),
                    labels=c("Below-ground\\nWet Season","Pitfalls\\nWet Season"))+
  #geom_text(data=raw_data_graphing_file, 
  #          aes(x=raw_data_NMDS3,y=raw_data_NMDS2, 
  #              label=Plot_abbre))+
  #geom_text(data=plot_level_species_graphing,
  #          aes(x=plot_level_species1,y=plot_level_species2,
  #              label=species_names_plot_level))+
  xlab("NMDS 3")+
  ylab("NMDS 2")+
  scale_colour_manual(values=c("brown","blue"),
                      breaks=c("below","pitfall_wet"),
                      labels=c("Below-ground\\nWet Season","Pitfalls\\nWet Season"))+
  annotate("text", x=min(raw_data_NMDS3)-.5, y=min(raw_data_NMDS2)-.5,
           label=paste('Stress =',round(raw_data_NMDS_3dim$stress,3)))+
  stat_ellipse(type='t',size=1, geom="polygon", alpha=.2)+ #t assumes multivariet t distribution
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_rect(colour="black", fill=NA),
        #panel.border=element_blank(), #gets rid of square going around the entire graph
        #axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top",
        legend.title=element_blank(),
        #legend.text=element_text(face="bold"),
        #legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title # family = "Arial"
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title # family = "Arial"
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text # family = "Arial"
        axis.text.y = element_text(size=12, color="black")) # family = "Arial"

NMDS1_NMDS2
NMDS1_NMDS3
NMDS3_NMDS2

## -- save figure output ####

# commenting out so as to now overwrite

NMDS1_NMDS2
#ggsave("Ruzi_etal_Figure1a.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = NMDS1_NMDS2,
#       path = figure_path) 

NMDS1_NMDS3
#ggsave("Ruzi_etal_Figure1b.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = NMDS1_NMDS3,
#       path = figure_path)

NMDS3_NMDS2
#ggsave("Ruzi_etal_Figure1c.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = NMDS3_NMDS2,
#       path = figure_path)




