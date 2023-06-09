## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
## CODE R4:                                                                  ###
## Phylogenetic signal analysis in Ethnobiology                              ###
## Gaoue et al. (2021) Biological Reviews, in press                          ###
## ogaoue@utk.edu | February 16, 2021                                        ###
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###

## This script to reproduce the analyses from:
## Gaoue, O. G.,J. K. Moutouama , M A. Coe, M. O. Bond, E. Green, N. B. Sero, B. S. Bezeng, and K. Yessoufou 2021. “Methodological advances for hypothesis-driven ethnobotany” Biological Review, in press.

## Code R4: tesing for phylogenetic signal in Ethnobiology  
## Dataset 3:
## data_5_shipibo_tree.txt
## data_4_shipibo_redundancy.csv
## data_6_ethno_shipibo.txt

rm(list=ls()) ## clearing the memory of the analysis

## Loading packages
library(caper)
library(picante)
library(phytools)

## Set your working directory
## setwd(/Users/Dropbox/Biological Reviews/Data/)

## Data inport
tree_rev <- read.tree("data_5_shipibo_tree.txt") # to read the phylogeny

tree_rev$node.label<-NULL ## to solve the error message "Labels duplicated between tips and nodes in phylogeny" (this doesn’t occur often)

data_rev <- read.csv("data_4_shipibo_redundancy.csv",header=TRUE) ## to read the dataframe where the ethnobotanical variables are found

## attach(data_rev)
names(data_rev) ## to check if the data_rev is actually the one I read in the line immediately above.

## Phylosignal test using Blomberg K test for continuous variables
## 1- Signal in UV?
phylosig(tree_rev, data_rev$uv_total, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list()) ## In this line, method can also be lambda if one wants to use lambda instead of K.

## Results: Phylogenetic signal K : 0.209566 
## P-value (based on 1000 randomizations): 0.469 # observed K is not different from rando==>meaning there is no signal in UV.

## 2-Signal in redundancy?
phylosig(tree_rev, data_rev$redundancy, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
## results: K=0.11, P=0.793 (1,000 randomisation) # observed K is not different from rando==>meaning there is no signal in UV.

## 3-signal in use pressure?
phylosig(tree_rev, data_rev$use_pressure, method="K", test=TRUE, nsim=1000, se=NULL, start=NULL, control=list())
## results: K=0.283, P=0.33 # observed K is not different from rando==>meaning there is no signal in UV.

## 4. Calculate NRI or NTI
sample_rev <- read.table("data_6_ethno_shipibo.txt", header=TRUE)

comm <- as.matrix(sample_rev)
prunedphy <- prune.sample(comm, tree_rev) ## making sure tips of the phylogeny match species names in the data
phydist <- cophenetic(tree_rev) ## phylogenetic distance between species

## Calculate NRI
ses.mpd.result <- ses.mpd(comm, phydist, null.model = "taxa.labels",
      abundance.weighted = FALSE, runs = 1000) ## calculate NRI (= -Z value)
ses.mpd.result

## Calculate NTI
ses.mntd.result <- ses.mntd(comm, phydist, null.model = "taxa.labels",
          abundance.weighted = FALSE, runs = 1000)
ses.mntd.result

## Figure: plot the tree showing in red preferred species by community in our case study

for (i in row.names(comm)) {
  plot(prunedphy, show.tip.label = FALSE, main = i,type = "fan")
  tiplabels(tip = which(prunedphy$tip.label %in% names(which(comm[i,] > 0
  ))), pch = 19, col="red",cex = 1)
}

