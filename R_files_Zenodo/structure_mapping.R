## Making sctructure pie charts and plotting them on a map

library(ggplot2)
library(scatterpie)
library(maptools)
library(maps)
library(rgdal)
library(raster)
library(plyr)
library(dplyr)

setwd("~/Desktop/")
states <- map_data("state")


## Make the map

setwd("~/Desktop/")
pop_loc <- read.csv("27pops_latlon.csv", header = T, stringsAsFactors = F) 
# load population cluster assignment information
pop_df <- read.csv("~/STRUCTURE.outfile.csv", header = F)
pop_df$long <- pop_loc$lon
pop_df$lat <- pop_loc$lat
colnames(pop_df) <- c("Popnum", "cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "indivcount", "long", "lat")

ggplot(data = subset(states,long > -100 & long < -75)) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_polygon(aes(x = long, y = lat, group = group), fill=NA,color = "black",size=0.25) + 
  scale_fill_manual(values = c("dodgerblue1", "darkorange2", "firebrick2", "green2", "gold1")) +
  geom_scatterpie(aes(x = long, y = lat, r = 0.5), data = pop_df, cols = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")) + 
  coord_fixed(1.3) +
  xlab("longitude") +
  ylab("latitude") +
  theme(legend.position="none") +
  NULL

## Delta K plotting
evanno <- read.csv("~/evanno.csv", header = T, stringsAsFactors = F)

ggplot(data = evanno, aes(x = K, y = Delta.K))+
  geom_point()+
  labs(x = "K", y = expression(Delta*"K"))+
  theme_bw()+
  NULL
  
