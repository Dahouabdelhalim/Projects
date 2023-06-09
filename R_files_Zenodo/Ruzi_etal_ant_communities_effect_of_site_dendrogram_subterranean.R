## ant communities by site and cache type
## dendrograms and hierarchical clustering for the subterranean samples

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
Ruzi_etal_data_below_frequency_condensed <- read_csv(paste(raw_data_path, "Ruzi_etal_data_below_frequency_condensed.csv", sep = "/"),
                                                     col_types = cols(
                                                       .default = col_double(),
                                                       ID = col_character(),
                                                       Plot_abbre = col_character(),
                                                       Season = col_character(),
                                                       Species = col_character(),
                                                       Experiment = col_character()))

Ruzi_etal_data_below_frequency_condensed

## -- hierarchical clustering ####

names(Ruzi_etal_data_below_frequency_condensed)  # the species start 2:21

cluster_below_mat <- as.matrix(Ruzi_etal_data_below_frequency_condensed[,2:21])
dimnames(cluster_below_mat)
cluster_below_mat_names <- as.vector(paste(Ruzi_etal_data_below_frequency_condensed$Plot_abbre, Ruzi_etal_data_below_frequency_condensed$Species,sep="_"))
dimnames(cluster_below_mat)[[1]] <- c(cluster_below_mat_names)
head(cluster_below_mat)

cluster_below_mat.dist <- vegdist(cluster_below_mat, method='bray')

# base cluster dendrogram plot
# will be redrawn to be saved as a figure later in script
below_flex1 <- as.hclust(flexbeta(cluster_below_mat.dist,-0.25))
plot(below_flex1)
coef(below_flex1) #0.8504618

## -- MRPP on groups based on the base cluster dendogram ####

Ruzi_etal_data_below_frequency_condensed$ID

below1 <- Ruzi_etal_data_below_frequency_condensed %>%
  mutate(new_group = if_else(ID == "25_Cec" |
                               ID == "25_Och" |
                               ID == "Z_P" |
                               ID == "A_Ape" |
                               ID == "A_Och" |
                               ID == "A_Jac" |
                               ID == "D_Cec" |
                               ID == "A_Cec" |
                               ID == "D_Ape" |
                               ID == "Z_Zan", "1", ID)) %>%
  mutate(new_group = if_else(new_group == "25_P" |
                               new_group == "25_TrBl" |
                               new_group == "25_Zan" |
                               new_group == "A_Zan" |
                               new_group == "A_G" |
                               new_group == "A_TrBl" |
                               new_group == "D_Jac" |
                               new_group == "D_G" |
                               new_group == "D_TrBl" |
                               new_group == "D_Och" |
                               new_group == "D_P", "2", new_group)) %>%
  mutate(new_group = if_else(new_group != "1" &
                               new_group != "2", "3", new_group))
below1
below1$new_group <- as_factor(below1$new_group)

# to see if groups (3) are different from each other -- they are significantly different from each other
# but there isn't a clear pattern as to why
cluster_below_cache.mrpp.3groups <- mrpp(dat=cluster_below_mat.dist, group=below1$new_group,
                                         permutations=999, distance="bray", strata=below1$Plot_abbre)# without the strata
#the values below are the same
cluster_below_cache.mrpp.3groups$Pvalue # [1] 0.001
cluster_below_cache.mrpp.3groups$A # [1] 0.2103564


## -- plotting the clusters in ggplot in a dendrogram ####
#extract the data of the dendrogram to then allow for nicer graphing
dend_below <- as.dendrogram(below_flex1)
dend_below_data <- dendro_data(dend_below, type="rectangle")
names(dend_below_data)

dend_below_data$segments
dend_below_data$labels

dend_data_tibble <- as_tibble(dend_below_data$labels)

dend_data_tibble2 <- dend_data_tibble %>%
  separate(label, into = c("plot", "spp"), sep = "_") %>%
  unite(labels2, plot, spp, sep = " ") %>%
  select(labels2)
dend_data_tibble2

dend_data_tibble$label <- as.character(dend_data_tibble$label)

dend_data_tibble3 <- dend_data_tibble %>%
  mutate(new_group = label) %>%
  mutate(new_group = if_else(new_group == "25_Cec" |
                               new_group == "25_Och" |
                               new_group == "Z_P" |
                               new_group == "A_Ape" |
                               new_group == "A_Och" |
                               new_group == "A_Jac" |
                               new_group == "D_Cec" |
                               new_group == "A_Cec" |
                               new_group == "D_Ape" |
                               new_group == "Z_Zan", "black", new_group)) %>%
  mutate(new_group = if_else(new_group == "25_P" |
                               new_group == "25_TrBl" |
                               new_group == "25_Zan" |
                               new_group == "A_Zan" |
                               new_group == "A_G" |
                               new_group == "A_TrBl" |
                               new_group == "D_Jac" |
                               new_group == "D_G" |
                               new_group == "D_TrBl" |
                               new_group == "D_Och" |
                               new_group == "D_P", "grey40", new_group)) %>%
  mutate(new_group = if_else(new_group != "black" &
                               new_group != "grey40", "grey80", new_group)) %>%
  mutate(group3 = new_group) %>%
  select(group3)
dend_data_tibble3

dend_below_data$labels2 <- as_vector(dend_data_tibble2$labels2)
dend_below_data$group3 <- as_vector(dend_data_tibble3$group3)

dend_below_data$group3

dend_below_data$group3 <- as.factor(dend_below_data$group3)
dend_below_data$group3

# color blind palette source from https://jfly.uni-koeln.de/color/ 
# but pulled from https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/

dend_data_tibble4 <- dend_data_tibble %>%
  mutate(site = sapply(strsplit(label, split = "_"), head, 1),
         site_color = if_else(site == "25", "#0072B2", site)) %>%
  mutate(site_color = if_else(site == "A", "#E69F00", site_color)) %>%
  mutate(site_color = if_else(site == "D", "#56B4E9", site_color)) %>%
  mutate(site_color = if_else(site == "P", "#009E73", site_color)) %>%
  mutate(site_color = if_else(site == "Z", "#D55E00" , site_color))

dend_data_tibble5 <- dend_data_tibble %>%
  mutate(site_abbre = sapply(strsplit(label, split = "_"), head, 1))  %>%
  mutate(cache_type = sapply(strsplit(label, split = "_"), tail, 1)) %>%
  mutate(site_num = if_else(site_abbre == "25", "1", site_abbre)) %>%
  mutate(site_num = if_else(site_abbre == "A", "2", site_num)) %>%
  mutate(site_num = if_else(site_abbre == "D", "3", site_num)) %>%
  mutate(site_num = if_else(site_abbre == "P", "4", site_num)) %>%
  mutate(site_num = if_else(site_abbre == "Z", "5", site_num)) %>%
  select(site_abbre, site_num, cache_type)

dend_below_data$site_color <- as_vector(dend_data_tibble4$site_color)
dend_below_data$site_abbre <- as_vector(dend_data_tibble5$site_abbre)
dend_below_data$site_num <- as_vector(dend_data_tibble5$site_num)
dend_below_data$cache_type <- as_vector(dend_data_tibble5$cache_type)

levels(dend_below_data$site_num) <- c(1:5)
levels(dend_below_data$site_num)

customized_dendrogram_below_site_color <- ggplot(dend_below_data$segments) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=1)+
  geom_text(data=dend_below_data$labels, aes(x,y,label=dend_below_data$cache_type,
                                             colour=dend_below_data$site_abbre),
            hjust=1, angle=0,size=5)+
  ylim(-.1,2)+
  ylab("Dissimilarity")+
  xlab("Communities")+
  geom_hline(yintercept=1.25, linetype="dashed", color = "grey54", size=1)+
  #scale_x_discrete(labels=c("1"="25 Cec", "2"="25 Och", "3"="Z Unbait", "4"="A Ape",
  #                          "5"="A Och", "6"="A Cec", "7"="A Jac", "8"="D Cec",
  #                          "9"="D Ape", "10"="Z Zan", "11"="25 Unbait", "12"="25 TrBl",
  #                          "13"="25 Zan", "14"="A Zan", "15"="A Beads", "16"="A TrBl",
  #                          "17"="D Jac", "18"="D Beads", "19"="D TrBl", "20"="D Och",
  #                          "21"="D Unbait", "22"="25 Beads", "23"="Z Och", "24"="P Unbait",
  #                          "25"="25 Jac", "26"="P TrBl", "27"="Z TrBl", "28"="Z Beads",
  #                          "29"="Z Unbait", "30"="Z Cec", "31"="P Jac", "32"="P Och",
  #                          "33"="Z Jac", "34"="D Zan", "35"="P Cec", "36"="P Beads"))+
  scale_x_discrete(labels=c("1"="Cec", "2"="Och", "3"="Unbait", "4"="Ape",
                            "5"="Och", "6"="Cec", "7"="Jac", "8"="Cec",
                            "9"="Ape", "10"="Zan", "11"="Unbait", "12"="TrBl",
                            "13"="Zan", "14"="Zan", "15"="Beads", "16"="TrBl",
                            "17"="Jac", "18"="Beads", "19"="TrBl", "20"="Och",
                            "21"="Unbait", "22"="Beads", "23"="Och", "24"="Unbait",
                            "25"="Jac", "26"="TrBl", "27"="TrBl", "28"="Beads",
                            "29"="Unbait", "30"="Cec", "31"="Jac", "32"="Och",
                            "33"="Jac", "34"="Zan", "35"="Cec", "36"="Beads"))+
  scale_color_manual(values = c("#0072B2", # 25Ha
                                "#E69F00", # Ava
                                "#56B4E9",  # Drayton
                                "#009E73", #Pearson
                                "#D55E00"), # Zetek
                     breaks = dend_below_data$site_abbre,
                     labels = dend_below_data$site_num) +
  labs(color = "Sites") +
  #scale_y_discrete(name="Height", limits)+
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
        legend.key=element_blank(), #remove the border around each item
        axis.title.x = element_text(size=14, color="black"), #size of x-axis title
        axis.title.y = element_text(size=14, color="black"), #size of y-axis title
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.text.x = element_text(size=12, color="black"), #size of x-axis text
        axis.text.y = element_text(size=12, color="black"),
        #axis.line.x	=element_blank(),
        axis.line.y	=element_blank())
customized_dendrogram_below_site_color

## -- save the dendrogram ####
# commenting out to not overwrite it
# boxes over the figure added in powerpoint

#ggsave("Ruzi_etal_Figure2b.png", width = 8, height = 8,
#       units = "in", dpi = 300, plot = customized_dendrogram_below_site_color,
#       path = figure_path)
