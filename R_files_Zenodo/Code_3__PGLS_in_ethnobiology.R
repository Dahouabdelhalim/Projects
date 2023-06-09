## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
## CODE R3:                                                                  ###
## Phylogenetic Generalized Least Square in ethnobiology                     ###
## Gaoue et al. (2021) Biological Reviews, in press                          ###
## ogaoue@utk.edu | February 16, 2021                                        ###
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ###

## This script to reproduce the analyses from:
## Gaoue, O. G.,J. K. Moutouama , M A. Coe, M. O. Bond, E. Green, N. B. Sero, B. S. Bezeng, and K. Yessoufou 2021. “Methodological advances for hypothesis-driven ethnobotany” Biological Review, in press.

## Code R3: Phylogenetic Generalized Least Square in ethnobiology

## Data 1: First you need data provided to build the backbone of your phylogeny that you will later prune for your list of species. This is from Qian, H. & Y. Jin. (2016) An updated megaphylogeny of plants, a tool for generating plant phylogenies and an analysis of phylogenetic community structure. Journal of Plant Ecology 9(2): 233–239

## Qian_PhytoPhylo.tre
## Qian_nodes.csv
## Qian_S.phyloMaker.R

## Data 2: Second you need your own phylogeny data, here the list of species and family, and then your ethnobiology data
## data_2_phylo.csv
## data_3_shipibo_phylo.nex
## data_4_shipibo_redundancy.csv

## This script is modified from an original R script file used for data analyses in the paper entitled: "Phylogeny explains why therapeutically redundant plant species are not necessarily facing greater use pressure" in People and Nature.

## The original analyses were conducted in 2018 using R version 3.4.3 "Kite-Eating Tree" (R Development Core Team, 2017) along with the APE (Analyses of Phylogenetics and Evolution) package version 5.0 (Paradis & Schliep, 2019).

rm(list=ls())

## load packages
library(ape)
library(caper)
library(geiger)
library(nlme)
library(phytools)
library(picante)

## Set your working directory
## setwd(/Users/Dropbox/Biological Review/Data/)

## 1. Build the phylogeny 
## Install the S.phyloMaker function that is used to build the phylogeny. This is coming from Qian, H. & Y. Jin. (2016) An updated megaphylogeny of plants, a tool for generating plant phylogenies and an analysis of phylogenetic community structure. Journal of Plant Ecology 9(2): 233–239.

source("Qian_S.phyloMaker.R")

## Data import: data files from the Phylomaker orginal source code
splist<-read.csv("data_2_phylo.csv", header=T)
tree=read.tree("Qian_PhytoPhylo.tre")
nodes <- read.csv("Qian_nodes.csv")

## S3 adjusts branch lenght using BLADJ as in phylocom
phylo_result=S.PhyloMaker(splist, tree, nodes, output.splist=T, scenarios=c("S1"))

## check taxa in tree
phylo_result$Species.list 

## Writing the phylogenetic tree out
write.tree(phylo_result$Scenario.1,"data_3_shipibo_phylo.nex", digits= 2 ) 

## reading the ethnobiology data
shipibo = read.csv("data_4_shipibo_redundancy.csv",row.names = 1) # Read in data
str(shipibo)

## read in the phylogeny tree
shipibo.tree = read.tree("data_3_shipibo_phylo.nex") 

class(shipibo.tree) ## check class to ensure tree in phylo

dev.off()
plot(shipibo.tree, type="fan", cex=.9) ## plot tree

name.check(shipibo.tree,shipibo) ## all names of taxa need to be consistent between tree and data file.

### PGLS test. First, make sure data for redundancy, preference and use_pressure correspond to the order of the tip labels

shipibo.tree$tip.label ## exact spelling, case, etc. is important in that all names of taxa shoud be consistent in both tree and data in CSV.

shipibo$tip_taxa ## check tip taxa

## Pylogenetics Generalized Least Squares (PGLS) regressions

## First, testing for a correlation among predictors and removing them. In this case, we were looking at the correlation between uv and redundancy and between preference and redunancy

cor.test(shipibo$uv_total, shipibo$redundancy)
cor.test(shipibo$preference, shipibo$redundancy)

## PGLS model candidates using a Brownian correlation structure. Starting with a full saturated model then removing one predictor or interaction with each subsequent model.

urm1<-gls(use_pressure ~ preference+ redundancy+ preference*redundancy,data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

urm2<-gls(use_pressure ~ preference+ redundancy, data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

urm3<-gls(use_pressure ~ preference, data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

urm4<-gls(use_pressure ~ redundancy, data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

urm5<-gls(use_pressure ~ 1, data=shipibo, correlation = corBrownian(1,form=~tip_taxa, phy=shipibo.tree))

## install bbmle for AIC tab
install.packages("bbmle")
library(bbmle)

AICtab(urm1, urm2, urm3, urm4, urm5) ## check model candidates for an AIC<2

## Viewing model summary for the best PGLS model candidate
summary(urm1)

## Without phylogenetic control generalized lineral model candidates. Starting with a full saturated model then removing one predictor or interaction with each subsequent model.

urmnc1<- glm(use_pressure ~ preference+ redundancy +preference*redundancy,  data=shipibo)

urmnc2<- glm(use_pressure ~ preference+ redundancy, data=shipibo)
urmnc3<- glm(use_pressure ~ preference,  data=shipibo)
urmnc4<- glm(use_pressure ~redundancy, data=shipibo)
urmnc5<- glm(use_pressure ~ 1, data=shipibo)

## Estimating AIC and selecting model candidates with a delta AIC<2

AICtab(urmnc1, urmnc2, urmnc3, urmnc4, urmnc5)

## Viewing model summary for the best model candidate without phylogenetic control
summary(urmnc3)

## Examining the predictive power or differenve in Akaike information criterion (AIC) between the best model's for PGLS and without phylogenetic control

AIC(urm1,urmnc3)
