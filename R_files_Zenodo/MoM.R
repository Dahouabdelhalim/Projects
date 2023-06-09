# code to estimate crown-group diversification

setwd("~/Documents/HYDIV")
rm(list = ls())

library(tidyverse)
library(ape)
library(BAMMtools)
library(coda)

# read in the family-level phylogeny

phy <- read.tree("data/full_tree.tre")

# tip names
tips <- phy$tip.label

nodes <- sapply(tips, function(x,y) which(y == x), y = phy$tip.label)

# named numeric of family names and edge lengths
edge.lengths <- setNames(phy$edge.length[sapply(nodes,
                                               function(x,y) which(y==x),y=phy$edge[,2])],names(nodes))
edge.lengths

# read in family size data
richness <- read.csv("data/FamSizes_use.csv")

# create new data frame with family names and branch lengths
MoM <- data.frame(Family = names(edge.lengths), Age = edge.lengths)
rownames(MoM) <- NULL 

# add richness info
MoM <- merge(MoM, richness, by.x = "Family", by.y = "taxon")

# calculate MoM diversification using epsilon = 0
MoM <- MoM %>% 
  rename(FamSize = n.taxa) %>% 
  mutate(divers0 = log(FamSize) / Age)
# check to make sure it looks okay
str(MoM)
head(MoM)

# now calculate MoM diversification using epsilon = 0.5 and epsilon = 0.9
MoM <- MoM %>% 
  mutate(divers0.5 = log(FamSize*(1-0.5) + 0.5)/(Age)) %>% 
  mutate(divers0.9 = log(FamSize*(1-0.9) + 0.9)/(Age)) 



