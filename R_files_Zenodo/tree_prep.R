## prepare tree data 

setwd("~/Documents/HYDIV")
rm(list = ls())

library(tidyverse)
library(ape)
library(phytools)

## qian tree
qian_tree <- read.tree("data/qian.tre")

to_drop <- c("Guamatelaceae", "Lophopyxidaceae", "Pennantiaceae",
             "Phellinaceae", "Thomandersiaceae", "Trimeniaceae", 
             "Nyssaceae", "Alzateaceae", "Amphorogynaceae", "Aptandraceae", 
             "Cervantesiaceae", "Comandraceae", "Coulaceae", "Crypteroniaceae",
             "Erythropalaceae", "Griseliniaceae", "Mazaceae", "Nanodeaceae",
             "Octoknemaceae", "Strombosiaceae", "Tetracarpaeaceae", 
             "Thesiaceae", "Viscaceae", "Ximeniaceae")

qian_trim <- drop.tip(qian_tree, to_drop)

# testo and sundue tree
ferns <- read.tree("data/testo_sundue.tre")
is.ultrametric(ferns)
is.binary(ferns)

# fern taxa: random species selected from each family
tips <- read.csv("data/fern_taxa.csv")

tips <- tips$Species
tips %in% ferns$tip.label

phylo_trim <- drop.tip(ferns, ferns$tip.label[-match(tips, ferns$tip.label)])

# bind the trees together
new_tree <- bind.tree(phylo_trim, qian_trim, where = 0, position = 0)
new_tree <- drop.tip(new_tree, "Welwitschia_mirabilis")

# write this tree
write.tree(new_tree, "data/prelim_tree.tre")

# this tree re-rooted and tip labels changed to families in Mesquite to produce the final tree
