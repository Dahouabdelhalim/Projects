## ant communities observed removing seeds of different tree species
## belowground (using subterranean traps associated with and without seed caches)

## -- load libraries ####
library(here)
library(vegan)
library(data.table)
library(tidyverse)
library(directlabels)

## -- set paths ####
raw_data_path <- here::here("data/raw_data/csv")
rdata_path <- here::here("data/rdata")
figure_path <- here::here("figs")

## -- read in data file ####

Ruzi_etal_data_below_frequency_condensed <- read_csv(paste(raw_data_path, "Ruzi_etal_data_below_frequency_condensed.csv", sep = "/"),
                                                     col_types = cols(
                                                       .default = col_double(),
                                                       ID = col_character(),
                                                       Plot_abbre = col_character(),
                                                       Season = col_character(),
                                                       Species = col_character(),
                                                       Experiment = col_character()))

Ruzi_etal_data_below_frequency_condensed


## -- convert data to a matrix ####

names(Ruzi_etal_data_below_frequency_condensed) # ants 2:21


below_mat <- as.matrix(Ruzi_etal_data_below_frequency_condensed[,2:21])
dimnames(below_mat)

unique(Ruzi_etal_data_below_frequency_condensed$ID)

belownames_mat <- as.vector(Ruzi_etal_data_below_frequency_condensed$ID)
dimnames(below_mat)[[1]] <- c(belownames_mat)
head(below_mat)

#to check if there are any na's anywhere in the dataset
csum_mat <- colSums(below_mat)
any(is.na(csum_mat)) # FALSE

## -- make the initial NMDS ####
## will comment out to prevent overwriting in future

#below_NMDS_2dim <- metaMDS(below_mat, distance="bray", k=2, trymax=50, autotransform=TRUE, pc=TRUE)
#stressplot(below_NMDS_2dim)
#below_NMDS_2dim$stress # 0.1717277

#below_NMDS_3dim <- metaMDS(below_mat, distance="bray", k=3, trymax=50, autotransform=TRUE, pc=TRUE)
#stressplot(below_NMDS_3dim)
#below_NMDS_3dim$stress #0.08737429

# saving the 3dim one

## -- save the metaMDS output ####

#save(below_NMDS_3dim, 
#        file = paste(rdata_path, "below_NMDS_3dim.rda", sep = "/"))

## -- load the saved Rdata metaMDS output ####

load(paste(rdata_path, "below_NMDS_3dim.rda", sep = "/"))

below_NMDS_3dim
below_NMDS_3dim$stress #0.08737429

## -- conduct MRPP ####

below.dist <- vegdist(below_mat, method='bray')

below.mrpp <- mrpp(dat=below.dist, grouping=Ruzi_etal_data_below_frequency_condensed$Species, permutations=999,
                   distance="bray", strata=Ruzi_etal_data_below_frequency_condensed$Plot_abbre) 
below.mrpp
below.mrpp$Pvalue #0.479 not significant
below.mrpp$distance
below.mrpp$A #-0.03932171

## -- prep to make the figures ####

below_NMDS1 <- below_NMDS_3dim$points[,1]
below_NMDS2 <- below_NMDS_3dim$points[,2]
below_NMDS3 <- below_NMDS_3dim$points[,3]
below_species1 <- below_NMDS_3dim$species[,1]
below_species2 <- below_NMDS_3dim$species[,2]
below_species3 <- below_NMDS_3dim$species[,3]
below_graphing_file <- cbind(Ruzi_etal_data_below_frequency_condensed, below_NMDS1, below_NMDS2, below_NMDS3)

below_species_graphing <- cbind(below_species1,below_species2, below_species3)
below_plot_level <- as.data.table(dimnames(below_species_graphing)[[1]])

names(below_graphing_file)
unique(below_graphing_file$Species)

below_graphing_file <- below_graphing_file %>%
  mutate(Species = if_else(Species == "G", "Beads", Species),
         Species = if_else(Species == "P", "Unbaited", Species))
below_graphing_file

## -- make the NMDS figure ####

below_NMDS1_NMDS2 <- ggplot(below_graphing_file, 
                            aes(x=below_NMDS1,y=below_NMDS2, 
                                fill=Species, colour=Species))+ # if go by plot because have two trials per season then get ellipses
  geom_text(data=below_graphing_file, 
            aes(x=below_NMDS1,y=below_NMDS2, 
                label=Species))+
  xlab("NMDS 1")+
  ylab("NMDS 2")+
  annotate("text", x=min(below_NMDS1)+.5, y=min(below_NMDS2)-.5,
           label=paste('Stress =',round(below_NMDS_3dim$stress,3)))+
  geom_polygon(alpha=0.2)+
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
below_NMDS1_NMDS2

below_NMDS1_NMDS3 <- ggplot(below_graphing_file, 
                            aes(x=below_NMDS1,y=below_NMDS3, 
                                fill=Species, colour=Species))+ # if go by plot because have two trials per season then get ellipses
  geom_text(data=below_graphing_file, 
            aes(x=below_NMDS1,y=below_NMDS3, 
                label=Species))+
  xlab("NMDS 1")+
  ylab("NMDS 3")+
  annotate("text", x=min(below_NMDS1)+.5, y=min(below_NMDS3)-.5,
           label=paste('Stress =',round(below_NMDS_3dim$stress,3)))+
  geom_polygon(alpha=0.2)+
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
below_NMDS1_NMDS3

below_NMDS3_NMDS2 <- ggplot(below_graphing_file, 
                            aes(x=below_NMDS3,y=below_NMDS2, 
                                fill=Species, colour=Species))+ # if go by plot because have two trials per season then get ellipses
  geom_text(data=below_graphing_file, 
            aes(x=below_NMDS3,y=below_NMDS2, 
                label=Species))+
  xlab("NMDS 3")+
  ylab("NMDS 2")+
  annotate("text", x=min(below_NMDS3)+.5, y=min(below_NMDS2)-.5,
           label=paste('Stress =',round(below_NMDS_3dim$stress,3)))+
  geom_polygon(alpha=0.2)+
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
below_NMDS3_NMDS2

## -- save figure output ####

# commenting out so as to now overwrite

#ggsave("Ruzi_etal_FigureS6a.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = below_NMDS1_NMDS2,
#       path = figure_path) 
#ggsave("Ruzi_etal_FigureS6b.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = below_NMDS1_NMDS3,
#       path = figure_path) 
#ggsave("Ruzi_etal_FigureS6c.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = below_NMDS3_NMDS2,
#       path = figure_path) 
