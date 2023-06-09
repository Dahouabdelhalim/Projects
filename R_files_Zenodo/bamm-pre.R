# get BAMM priors 

setwd("~/Documents/HYDIV")
rm(list = ls())

library(tidyverse)
library(ape)
library(BAMMtools)
library(coda)

# read in the family-level phylogeny

# get sampleProbs 1/number of species in the family

phy <- read.tree("data/full_tree.tre")
phy$tip.label

sampleProbs <- data.frame(speciesName = phy$tip.label, cladeName = phy$tip.label)

fam_use <- read.csv("data/FamSizes_use.csv", header = TRUE)

sampleProbs <- merge(sampleProbs, fam_use, by.x = "cladeName", by.y = "taxon")

library(dplyr)

# make samplePros file
sampleProbs <- sampleProbs %>% mutate(samplingFraction = 1/n.taxa)

write.table(sampleProbs, "sampleProbs.txt", sep = "\\t", row.names = FALSE)


## set BAMM priors
library(BAMMtools)

setBAMMpriors(phy, total.taxa = 315824)
