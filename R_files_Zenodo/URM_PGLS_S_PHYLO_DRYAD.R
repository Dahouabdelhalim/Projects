+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Phylogenetic Generalized Least Square in ethnobiology
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Coe and Gaoue 2021 PGLS script
# Michael A. Coe | coem@hawaii.edu
# February 12, 2021
  
## Code : Phylogenetic Generalized Least Square in ethnobiology
  
## Data 1: First you need data provided to build the backbone of your phylogeny that you will later prune for your list of species. This is from Qian, H. & Y. Jin. (2016) An updated megaphylogeny of plants, a tool for generating plant phylogenies and an analysis of phylogenetic community structure. Journal of Plant Ecology 9(2): 233â€“239
  
  ## Qian_PhytoPhylo.tre
  ## Qian_nodes.csv
  ## Qian_S.phyloMaker.R

## Data 2: Second you need your own phylogeny data, here the list of species and family, and then your ethnobiology data
## data_2_phylo.csv
## data_3_shipibo_phylo.nex
## data_4_shipibo_redundancy.csv


# This script adopted and modified from an original R script file used for data analyses in the paper entitled: "Phylogeny explains why less therapeutically redundant plant species are not necessarily facing greater use pressure" in People and Nature Journal. In Press.

# The original analyses were conducted in 2018 using R version 3.4.3 "Kite-Eating Tree" (R Development Core Team, 2017) along with the APE (Analyses of Phylogenetics and Evolution) package version 5.0 (Paradis & Schliep, 2019).

rm(list=ls())

# Install S.phyloMaker function In the Functions File
source("") # set source

# install packages
library(ape)
library(caper)
library(geiger)
library(nlme)
library(phytools)
library(picante)

## data import
setwd("") #set working directory
splist<-read.csv("phylomaker.csv", header=T)
tree=read.tree("PhytoPhylo.tre")
nodes <- read.csv("nodes.csv")

#S3 adjusts branch lenght using BLADJ as in phylocom
phylo_result=S.PhyloMaker(splist, tree, nodes, output.splist = T, scenarios = c("S1"))

phylo_result$Species.list # check taxa in tree

write.tree(phylo_result$Scenario.1,"shipibo_urm_s_phylo.nex", digits= 2 ) # writing tree

shipibo = read.csv("URM_S_PHYLO.csv",row.names = 1) # Read in data

str(shipibo)

shipibo.tree = read.tree("shipibo_urm_s_phylo.nex") # read in tree to build phylogeny

class(shipibo.tree) # check class to ensure tree in phylo

plot(shipibo.tree, type="fan",cex=0.38) # plot tree

name.check(shipibo.tree,shipibo) # all names of taxa need to be consistent between tree and data file.

### Phylogenetic signal test. First, make sure data for redundancy, preference and use_pressure correspond to the order of the tip labels

shipibo.tree$tip.label # exact spelling, case, etc. is important in that all names of taxa shoud be consistent in both tree and data in CSV.

shipibo$tip_taxa # check tip taxa

# PGLS regressions

# First, testing for a correlation among predictors and removing them. In this case, we were looking at the correlation between uv and redundancy and between preference and redunancy

cor.test(shipibo$uv_total, shipibo$redundancy)

cor.test(preference, redundancy)


# PGLS model candidates using a Brownian correlation structure. Starting with a full saturated model then removing one predictor or interaction with each subsequent model.

urm1<- gls(use_pressure ~ preference+ redundancy+ preference*redundancy,data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

urm2<- gls(use_pressure ~ preference+ redundancy, data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

urm3<- gls(use_pressure ~ preference, data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

urm4<- gls(use_pressure ~ redundancy, data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

urm5<- gls(use_pressure ~ 1, data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

# install bbmle for AIC tab
install.packages("bbmle")
library(bbmle)

AICtab(urm1, urm2, urm3, urm4, urm5) # check model candidates for an AIC<2

# Viewing model summary for the best PGLS model candidate
summary(urm1)


# Without phylogenetic control generalized lineral model candidates. Starting with a full saturated model then removing one predictor or interaction with each subsequent model.

urmnc1<- glm(use_pressure ~ preference+ redundancy +preference*redundancy,  data=shipibo)

urmnc2<- glm(use_pressure ~ preference+ redundancy, data=shipibo)

urmnc3<- glm(use_pressure ~ preference,  data=shipibo)

urmnc4<- glm(use_pressure ~redundancy, data=shipibo)

urmnc5<- glm(use_pressure ~ 1, data=shipibo)


# Estimating AIC and selecting model candidates with a delta AIC<2

AICtab(urmnc1, urmnc2, urmnc3, urmnc4, urmnc5)

# Viewing model summary for the best model candidate without phylogenetic control
summary(urmnc3)

# Examining the predictive power or differenve in Akaike information criterion (AIC) between the best model's for PGLS and without phylogenetic control

 AIC(urm1,urmnc3)



