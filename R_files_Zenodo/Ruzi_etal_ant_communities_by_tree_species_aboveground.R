## ant communities observed removing seeds of different tree species
## aboveground (so data from observing seed removal caches and from pitfall traps)

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
Ruzi_etal_data_removal_pitfalls_met_criteria <- read_csv(paste(raw_data_path, "Ruzi_etal_data_removal_pitfalls_met_criteria.csv", sep = "/"),
                         col_types = cols(
                           .default = col_double(),
                           ID = col_character(),
                           Plot_abbre = col_character(),
                           Season = col_character(),
                           Species = col_character(),
                           Experiment = col_character(),
                           Plot_season = col_character()
                         ))

Ruzi_etal_data_removal_pitfalls_met_criteria

## -- convert data to a matrix ####
names(Ruzi_etal_data_removal_pitfalls_met_criteria) # ant species 2:57

removal_crit_mat <- as.matrix(Ruzi_etal_data_removal_pitfalls_met_criteria[,2:57])
dimnames(removal_crit_mat)[[1]]
removalnames_mat_crit <- as.vector(Ruzi_etal_data_removal_pitfalls_met_criteria$ID)
dimnames(removal_crit_mat)[[1]] <- c(removalnames_mat_crit)
head(removal_crit_mat)

#to check if there are any na's anywhere in the dataset
csum_mat <- colSums(removal_crit_mat)
any(is.na(csum_mat)) # FALSE

## -- make the initial NMDS ####
## will comment out to prevent overwriting in future

#removal_NMDS_2dim_crit <- metaMDS(removal_crit_mat, distance="bray", k=2, trymax=50, autotransform=TRUE, pc=TRUE)
#stressplot(removal_NMDS_2dim_crit)
#removal_NMDS_2dim_crit$stress # 0.07988232

#removal_NMDS_3dim_crit <- metaMDS(removal_crit_mat, distance="bray", k=2, trymax=50, autotransform=TRUE, pc=TRUE)
#stressplot(removal_NMDS_3dim_crit)
#removal_NMDS_3dim_crit$stress #  0.07988352 -- using the 2dim solution

# will keep the 2dim solution

## -- save the metaMDS output ####

#save(removal_NMDS_2dim_crit, 
#     file = paste(rdata_path, "removal_pitfalls_NMDS_2dim_crit.rda", sep = "/"))


## -- load the saved Rdata metaMDS output ####
load(paste(rdata_path,"removal_pitfalls_NMDS_2dim_crit.rda", sep = "/"))

removal_NMDS_2dim_crit
removal_NMDS_2dim_crit$stress # 0.07988232



## -- conduct MRPP ####
## -- tree species pooled versus pitfall traps

removal2_crit <- Ruzi_etal_data_removal_pitfalls_met_criteria %>%
  mutate(new_cat = if_else(Species == "P", "pitfall", Species)) %>%
  mutate(new_cat = if_else(Species != "P", "above", new_cat))
removal2_crit

removal_crit.dist <- vegdist(removal_crit_mat, method='bray')

removal_crit.mrpp <- mrpp(dat=removal_crit.dist, grouping=removal2_crit$new_cat, permutations=999,
                          distance="bray", strata=removal2_crit$Plot_abbre) 
removal_crit.mrpp
removal_crit.mrpp$Pvalue # significantly different 0.005
removal_crit.mrpp$distance # "bray"
removal_crit.mrpp$A # 0.09763321

## -- tree species versus each other and versus pitfall traps

removal.species_crit.mrpp <- mrpp(dat=removal_crit.dist, grouping=removal2_crit$Species, permutations=999,
                                  distance="Bray", strata=removal2_crit$Plot_abbre)
removal.species_crit.mrpp
removal.species_crit.mrpp$Pvalue #significant 0.006
removal.species_crit.mrpp$A # 0.126609

## -- prep to make the figures and fit environmental conditions ####

removal_NMDS1_crit <- removal_NMDS_2dim_crit$points[,1]
removal_NMDS2_crit <- removal_NMDS_2dim_crit$points[,2]
removal_species1_crit <- removal_NMDS_2dim_crit$species[,1]
removal_species2_crit <- removal_NMDS_2dim_crit$species[,2]
removal_graphing_file_crit <- cbind(removal2_crit, removal_NMDS1_crit, removal_NMDS2_crit)

removal_species_graphing_crit <- cbind(removal_species1_crit,removal_species2_crit)
removal_plot_level_crit <- as.data.table(dimnames(removal_species_graphing_crit)[[1]])

fit_plot_removal_crit <- envfit(removal_NMDS_2dim_crit, removal_crit_mat)
arrow <- data.frame(fit_plot_removal_crit$vectors$arrows,R=fit_plot_removal_crit$vectors$r,
                    P=fit_plot_removal_crit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p<-subset(arrow, P <= 0.05)
arrow.p # so 4 species contribute to the differences
#           NMDS1      NMDS2         R     P     FG
#AphAra  0.98582618 -0.1677699 0.7705608 0.007 AphAra
#EctRui -0.37401955 -0.9274208 0.8075026 0.001 EctRui
#ParBug -0.09282372  0.9956826 0.4784402 0.014 ParBug
#ParCor -0.02609240  0.9996595 0.6696161 0.007 ParCor

## -- make the NMDS figure ####
removal_graph_crit <- ggplot(removal_graphing_file_crit, 
                             aes(x=removal_NMDS1_crit,y=removal_NMDS2_crit, 
                                 fill=new_cat, colour=new_cat))+ # if go by plot because have two trials per season then get ellipses
  geom_text(data=removal_graphing_file_crit, 
            aes(x=removal_NMDS1_crit,y=removal_NMDS2_crit, 
                label=Species))+
  xlab("NMDS 1")+
  ylab("NMDS 2")+
  #xlim(c(-2,4))+
  annotate("text", x=max(removal_NMDS1_crit)-.8, y=min(removal_NMDS2_crit)-0.5,
           label=paste('Stress =',round(removal_NMDS_2dim_crit$stress,3)))+
  stat_ellipse(type='t',size=1, geom="polygon", alpha=.2)+ #t assumes multivariet t distribution
  scale_fill_manual(values=c("grey60","purple"),
                    breaks=c("pitfall","above"),
                    labels=c("Pitfall\\nTraps","Observed\\nRemoval"))+
  scale_colour_manual(values=c("grey60","purple"),
                      breaks=c("pitfall","above"),
                      labels=c("Pitfall\\nTraps","Observed\\nRemoval"))+
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
removal_graph_crit ## Figure 3A


removal_crit_ant_species <- ggplot(removal_graphing_file_crit, 
                                   aes(x=removal_NMDS1_crit,y=removal_NMDS2_crit, 
                                       fill=Species, colour=Species))+ # if go by plot because have two trials per season then get ellipses
  #geom_point(position=position_jitter(.1),shape=3)+ 
  geom_text(data=removal_graphing_file_crit, 
            aes(x=removal_NMDS1_crit,y=removal_NMDS2_crit, 
                label=Species))+ #label=Species
  geom_polygon(alpha=0.2)+
  #geom_text(data=plot_level_species_graphing,
  #          aes(x=plot_level_species1,y=plot_level_species2,
  #              label=species_names_plot_level))+
  xlab("NMDS 1")+
  ylab("NMDS 2")+
  xlim(c(-2,4))+
  guides(fill = FALSE, color = FALSE)+
  geom_segment(data=arrow.p, aes(x=0, y=0, xend=NMDS1, yend=NMDS2, lty=FG), #lty=FG, colour=FG
               arrow=arrow(length=unit(.2,"cm")*arrow.p$R), inherit.aes=FALSE)+
  geom_dl(data=arrow.p, aes(x=NMDS1-.05, y=NMDS2+0.1, label = FG), #colour=FG
          method = list(dl.trans(x = x + 0.2), "last.bumpup", cex = 0.60),
          inherit.aes=FALSE) +
  annotate("text", x=min(removal_NMDS1_crit)+.7, y=min(removal_NMDS2_crit)-.5,
           label=paste('Stress =',round(removal_NMDS_2dim_crit$stress,3)))+
  #annotate("text", x=min(pitfalls_NMDS1)+.5, y=min(pitfalls_NMDS2)-.65,
  #         label="perMANOVA: F=4.89, df=2, P=0.001")+
  #stat_ellipse(type='t',size=1, geom="polygon", alpha=.2)+ #t assumes multivariet t distribution
  #scale_fill_manual(values=c("red","blue"),
  #                  breaks=c("pitfall","above"),
  #                  labels=c("Pitfalls","Removal"))+
  #scale_colour_manual(values=c("red","blue"),
  #                    breaks=c("pitfall","above"),
  #                    labels=c("Pitfalls","Removal"))+
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
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"))
#axis.line.x	=element_line(colour="black"),
#axis.line.y =element_line(colour="black")) #size of y-axis text
removal_crit_ant_species ## Figure 3b


## -- save figure output ####

# commenting out so as to now overwrite
#ggsave("Ruzi_etal_Figure3a.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = removal_graph_crit,
#       path = figure_path)

#ggsave("Ruzi_etal_Figure3b.png", width = 6, height = 6,
#       units = "in", dpi = 300, plot = removal_crit_ant_species,
#       path = figure_path)

