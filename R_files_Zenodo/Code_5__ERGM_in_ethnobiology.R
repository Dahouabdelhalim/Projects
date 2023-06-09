## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
## CODE R5:                                                                  ###
## Exponential random growth model (ERGM) in ethnobiology                    ###
## Gaoue et al. (2021) Biological Reviews, in press                          ###
## ogaoue@utk.edu | February 16, 2021                                        ###
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###

## This script to reproduce the analyses from:
## Gaoue, O. G.,J. K. Moutouama , M A. Coe, M. O. Bond, E. Green, N. B. Sero, B. S. Bezeng, and K. Yessoufou 2021. “Methodological advances for hypothesis-driven ethnobotany” Biological Review, in press.

## Code R5: Exponential random growth model (ERGM) in ethnobiology   
## Dataset 4: Data used in this example is available at: 
## Bond, Matthew; Gaoue, Orou (2020), Adjacency matrices and nodal attributes for prestige and homophily predict network structure for social learning of medicinal plant knowledge, Dryad, Dataset, https://datadryad.org/stash/dataset/doi:10.5061/dryad.cfxpnvx3q

## data_7_adj_knowledge.csv
## data_8_adj_married.csv
## data_9_node_attributes.csv

rm(list=ls())

## Loading packages
library(network)
library(intergraph)
library(igraph)
library(ergm)

## Set your working directory
## setwd(/Users/Dropbox/Biological Reviews/Data/)

## 1. Import data: load adjacency matrices as data frames
knowledge_sharing = read.csv("data_7_adj_knowledge.csv", row.names=1) 
knowledge_net = as.network.matrix(knowledge_sharing)  # convert to network
spouses = read.csv("data_8_adj_married.csv", row.names=1)   
spouse_net = as.network.matrix(spouses)   # convert to network

nodes_all = read.csv("data_9_node_attributes.csv", row.names=1)   # load data frame
nodes = nodes_all[!is.na(nodes_all$uses),]   # remove rows containing NA

## 2. Set node attributes
knowledge_net%v%"village" = nodes$village   
knowledge_net%v%"age" = nodes$age   
knowledge_net%v%"use_#" = nodes$uses
knowledge_net%v%"use_uniqueness" = nodes$ab

## 3. Convert network to igraph
k_net <- asIgraph(knowledge_net) 

## 4. Calculate node-level connectivity & save to the “nodes” data frame
nodes$betweenness = betweenness(k_net) 
nodes$closeness = closeness(k_net)
out_d_cen = centr_degree(k_net, mode="out")  #could also use mode="in"
nodes$out_degree = out_d_cen$res

## 5. Calculate network-level connectivity
transitivity(k_net, type="global")
edge_density(k_net, loops = FALSE)
diameter(k_net)

## these metrics can be used in future analyses or compared between groups (eg mean betweenness of each village)

## 6. run exponential random graph models (ERGMs), compare AIC
set.seed(99) # to make sure results are reproducible
## M1: Model with only configuration variables
m1 = ergm(knowledge_net ~ edges + isolates + gwodegree(0.5) + mutual) 
summary(m1) 

# M2: Model with configuration variables and node-level covariates
m2 = ergm(knowledge_net ~ edges + isolates + gwodegree(0.5) + mutual + nodeocov('age') + nodeocov('use_#') + nodeocov('use_uniqueness') +  edgecov(spouse_net) + absdiff('age') + absdiff('use_#'))
summary(m2) 

## We see that adding node-level covariates greatly improves the model fit (AIC reduced by 700). We see that use_uniqueness has very high MCMC % (99- normal values are <5), which means the model isn't fitting it well- let's drop it and try again

## M3: Model with the node 'use_uniqueness' droped
m3 = ergm(knowledge_net ~ edges + isolates + gwodegree(0.5) + mutual + nodeocov('age') + nodeocov('use_#') + edgecov(spouse_net) + absdiff('age') + absdiff('use_#'))
summary(m3) 

## This change drops AIC slightly and model statistics (below) are greatly improved

## Significant results provide support for several hypotheses (Figure 7, Table 1): edges- prestige/rarity, gwodegree- prestige, mutual- homophily, nodeocov.age - prestige, nodeocov.use_# - prestige, edgecov.spouse_net -homophily, absdiff.age - prestige, absdiff.use_#  - control- yes people who say they learned from someone have more similar knowledge to that person. 

## 7. Evaluate model fit
mcmc.diagnostics(m3, center=T) ## density plots should be symmetrical and centered on 0. m3 has a better fit than m2, removing the variable really helped the model

## Figure
par(mfrow = c(2, 3))
plot(ergm::gof(m3)) # thick black line represents observed network, grey histograms represent simulated networks. The better the model, the more the black line and histograms will align exactly. 

## convert log odds model coefficients into probability using the logistic function 
plogis(5.69) ## this means that there is a 99.7% chance of people sharing knowledge if they are spouses, supports a homophily effect (See Table 1)
plogis(0.02) ## this means that there is a 50.5% chance of people sharing knowledge if they are close in age- this is significantly greater than random (50%), but not much greater. Still, this supports a homophily effect (See Table 1)

