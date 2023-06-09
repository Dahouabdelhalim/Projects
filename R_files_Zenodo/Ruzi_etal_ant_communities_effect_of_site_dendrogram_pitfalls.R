## ant communities by site
## dendrograms and hierarchical clustering for the pitfall samples

## -- load libraries ####
library(here)
library(vegan)
library(data.table)
library(tidyverse)
library(cluster)
library(ggdendro)
library(dendextend)

## -- set file paths ####
raw_data_path <- here::here("data/raw_data/csv")
figure_path <- here::here("figs")

## -- flexbeta function ####
flexbeta <- function (dis,beta) 
{
  alpha <- (1-beta)/2
  out <- agnes(dis,meth='flex',par.method=alpha)
  out
}

## -- load in the data file ####
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

## -- hierarchical clustering ####
names(Ruzi_etal_freq_data_pitfalls_dry_wet_season) # ant speices from 2:60

cluster_above_mat <- as.matrix(Ruzi_etal_freq_data_pitfalls_dry_wet_season[,2:60])
dimnames(cluster_above_mat)
cluster_above_mat_names <- as.vector(Ruzi_etal_freq_data_pitfalls_dry_wet_season$ID)
dimnames(cluster_above_mat)[[1]] <- c(cluster_above_mat_names)
head(cluster_above_mat)

cluster_above_mat.dist <- vegdist(cluster_above_mat, method='bray')
cluster_above_mat.dist

# base cluster dendrogram plot
# will be redrawn to be saved as a figure later in script
above_flex1 <- as.hclust(flexbeta(cluster_above_mat.dist,-0.25))
plot(above_flex1)
coef(above_flex1) # 0.5712407

## -- MRPP on groups based on the base cluster dendogram ####
Ruzi_etal_freq_data_pitfalls_dry_wet_season$ID
#[1] "25_dry" "25_wet" "A_dry"  "A_wet"  "D_dry"  "D_wet"  "P_dry"  "P_wet"  "Z_dry" 
#[10] "Z_wet"

# so based on that and the plot then the grouping are (following the order of ID)
Ruzi_etal_freq_data_pitfalls_dry_wet_season$group3 <- c(1,1,2,2,3,3,1,2,2,2)

#to see if groups (3) are different from each other
cluster_above.mrpp.plot_3groups <- mrpp(dat=cluster_above_mat.dist, grouping=Ruzi_etal_freq_data_pitfalls_dry_wet_season$group3,
                                        permutations=999, distance="bray", strata=Ruzi_etal_freq_data_pitfalls_dry_wet_season$Season)
cluster_above.mrpp.plot_3groups$Pvalue # 0.004
cluster_above.mrpp.plot_3groups$A # 0.1886003
# the groups are significant at the 3 level

## -- plotting the clusters in ggplot in a dendrogram ####
#extract the data of the dendrogram to then allow for nicer graphing
dend_above <- as.dendrogram(above_flex1)
dend_above_data <- dendro_data(dend_above, type="rectangle")
names(dend_above_data)

dend_above_data$segments
dend_above_data$labels
# labels2 follows the order of labels
dend_above_data$labels2 <- c("25 Dry","P Dry", "25 Wet","A Dry","A Wet",
                             "P Wet","Z Dry",
                             "Z Wet", "D Dry","D Wet")
dend_above_data$site_abbre <- c("25","P", "25","A","A",
                                "P","Z",
                                "Z", "D","D")
dend_above_data$season <- c("Dry","Dry", "Wet","Dry","Wet",
                            "Wet","Dry",
                            "Wet", "Dry","Wet")
# color blind palette source from https://jfly.uni-koeln.de/color/ 
# but pulled from https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/
dend_above_data$site_color <- c("#0072B2", "#009E73", 
                                "#0072B2", "#E69F00",
                                "#E69F00", "#009E73",
                                "#D55E00", "#D55E00",
                                "#56B4E9", "#56B4E9")
dend_above_data

customized_dendrogram_above_site_colors <- ggplot(dend_above_data$segments) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=1)+
  geom_text(data=dend_above_data$labels, aes(x,y,label=dend_above_data$season,
                                             color = dend_above_data$site_abbre),
            hjust=1, angle=0,size=5)+
  ylim(-.1,1.15)+
  ylab("Dissimilarity")+
  xlab("Communities")+
  #scale_x_discrete(labels=c("1" = "25 Dry",
  #                          "2" = "P Dry", 
  #                          "3" = "25 Wet",
  #                          "4" = "A Dry",
  #                          "5" = "A Wet",
  #                          "6" = "P Wet",
  #                          "7" = "Z Dry",
  #                          "8" = "Z Wet", 
  #                          "9" = "D Dry",
  #                          "10" = "D Wet"))+
  scale_x_discrete(labels=c("1" = "Dry",
                            "2" = "Dry", 
                            "3" = "Wet",
                            "4" = "Dry",
                            "5" = "Wet",
                            "6" = "Wet",
                            "7" = "Dry",
                            "8" = "Wet", 
                            "9" = "Dry",
                            "10" = "Wet"))+
  scale_color_manual(values = c("#0072B2", # 25Ha
                                "#E69F00", # Ava
                                "#56B4E9",  # Drayton
                                "#009E73", #Pearson
                                "#D55E00"), # Zetek
                     breaks = dend_above_data$site_abbre,
                     labels = c("1", "4", "1", "2", "2", "4", "5", "5", "3", "3")) +
  labs(color = "Sites") +
  #scale_y_discrete(name="Height", limits)+
  geom_hline(yintercept=.7, linetype="dashed", color = "grey54", size=1)+
  coord_flip()+
  theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        panel.border=element_blank(), #gets rid of square going around the entire graph
        axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks.y=element_line(colour = 'black', size = 0.5), #sets the tick lines
        axis.ticks.x=element_blank(),
        #legend.position=c(0.10,0.25), # moves the location of the legend
        legend.position="top",
        #legend.background=element_blank(), # removes the overall border
        #legend.background=element_rect(fill="white", colour="black"), #puts a black box around the legend
        #legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"),
        #axis.line.x	=element_blank(),
        axis.line.y	=element_blank())
customized_dendrogram_above_site_colors

## -- save the dendrogram ####
# commenting out to not overwrite it
# boxes over the figure added in powerpoint


#ggsave("Ruzi_etal_Figure2a.png", width = 8, height = 3,
#       units = "in", dpi = 300, plot = customized_dendrogram_above_site_colors,
#       path = figure_path) 
