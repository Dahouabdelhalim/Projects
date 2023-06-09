# In the Name of God
# Comparing Methods R script
# Question: Why use 'LA' 
# in calculation of first k eigenvalues
# 'LA' mean largest algebraically
# is there a reason? //???????

#### Description ####
# Changed all variables to local. No global is left. 
# Also added comments to make the code more readable

# Version Dates:
# V1: 13980224
# V2: 13980225

#### Preparations ####
# setwd('C:/Users/Fava/Desktop/Dars/spm-mux/')

# Clearing the workspace 
remove(list = ls())

# Start timing
start_time = Sys.time()
# Setting working directory to access functions and data sets
setwd("~/GitHub/utlp-spm-mux")

# Adding needed libraries 
library(igraph)
library(LinkPrediction)
library(lattice)

# Adding needed source codes 
# Creates adjacency matrix from edge list 
source("form_net.R")
# assesment function that is capable of handling other layers as input
source("prune_recover_MultipleLayer.R")

source("lp_spm_mux_coeffs.R")
source('structural_similarity.R')
source('lp_rnd.R')
source('lp_mux_spm.R')
source('lp_spm_custom.R')

source('LinkPrediction_Algorithm0.R')


#### Preparing Data Sets ####
network_layer1_name <- 'aarhus_l1'
network_layer2_name <- 'aarhus_l2'
network_layer3_name <- 'aarhus_l3'
network_layer4_name <- 'aarhus_l4'
network_layer5_name <- 'aarhus_l5'

net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')
net_l3_fpath <- paste0('mux-data/', network_layer3_name, '.txt')
net_l4_fpath <- paste0('mux-data/', network_layer4_name, '.txt')
net_l5_fpath <- paste0('mux-data/', network_layer5_name, '.txt')

# The main layer to be predicted is always layer A. 
# All other layers act as auxiliary 
gA <- graph_from_adjacency_matrix(form_net(net_l1_fpath), mode = "undirected", diag = FALSE)
gB <- graph_from_adjacency_matrix(form_net(net_l2_fpath), mode = "undirected", diag = FALSE)
gC <- graph_from_adjacency_matrix(form_net(net_l3_fpath), mode = "undirected", diag = FALSE)
gD <- graph_from_adjacency_matrix(form_net(net_l4_fpath), mode = "undirected", diag = FALSE)
gE <- graph_from_adjacency_matrix(form_net(net_l5_fpath), mode = "undirected", diag = FALSE)

# Making sure that all graphes are the same size 
maximumVcount = 
  max(vcount(gA),vcount(gB),vcount(gC),vcount(gD),vcount(gE))
if(vcount(gA) < maximumVcount) {
  gA = add_vertices(gA, maximumVcount - vcount(gA))
}
if(vcount(gB) < maximumVcount) {
  gB = add_vertices(gB, maximumVcount - vcount(gB))
}
if(vcount(gC) < maximumVcount) {
  gC = add_vertices(gC, maximumVcount - vcount(gC))
}
if(vcount(gD) < maximumVcount) {
  gD = add_vertices(gD, maximumVcount - vcount(gD))
}
if(vcount(gE) < maximumVcount) {
  gE = add_vertices(gE, maximumVcount - vcount(gE))
}

gRA <<-erdos.renyi.game(nrow(as_adjacency_matrix(gA, type = "both")), graph.density(gA), type = c("gnp"), directed = FALSE, loops = FALSE)
gRB <<-erdos.renyi.game(nrow(as_adjacency_matrix(gB, type = "both")), graph.density(gB), type = c("gnp"), directed = FALSE, loops = FALSE)


#### Analysis ####
# prune_recover_MultipleLayer gets the complete graph,
# randomly splits it into train and test subgraphs 
# according to probes 
auxiliaryLayers = list(gD)
assessment <- 
  prune_recover_MultipleLayer(
    gA,
    "LinkPrediction_Algorithm0",
    probes = seq(0.1, 0.9, 0.1), 
    epochs = 30, 
    preserve_conn = FALSE, 
    use_weights = FALSE, 
    auxiliaryLayers = auxiliaryLayers # NA means SPM, auxiliaryLayers uses these layers
    )
# Description of different methods:
# lp_spm_custom: Basic SPM method. In other words, lp

# metric: recall_at_k, aupr, auroc, avg_prec
plot_lp_performance(assessment, metric = "auroc", err = "sd")


# Save an object to a file
saveRDS(assessment, file = "assessment_aarhus_A_D.rds")
# saveRDS(assessment, file = "assessment_aarhus_A_RB.rds")
# Restore the object
assessment_aarhus_A =
  readRDS(file = "assessment_aarhus_A.rds")
assessment_aarhus_A$method = "SPM"

assessment_aarhus_A_D =
  readRDS(file = "assessment_aarhus_A_D.rds")
assessment_aarhus_A_D$method = "Leisure"

assessment_aarhus_A_E =
  readRDS(file = "assessment_aarhus_A_E.rds")
assessment_aarhus_A_E$method = "Work"

assessment_aarhus_A_BCDE = 
  readRDS(file = "assessment_aarhus_A_BCDE.rds")
assessment_aarhus_A_BCDE$method = "F C L W"

assessment_aarhus_A_B = 
  readRDS(file = "assessment_aarhus_A_B.rds")
assessment_aarhus_A_B$method = "FB"

assessment_aarhus_A_C = 
  readRDS(file = "assessment_aarhus_A_C.rds")
assessment_aarhus_A_C$method = "Coauthor"

assessment_aarhus_A_BC = 
  readRDS(file = "assessment_aarhus_A_BC.rds")
assessment_aarhus_A_BC$method = "BC"

assessment_aarhus_A_RB = 
  readRDS(file = "assessment_aarhus_A_RB.rds")
assessment_aarhus_A_RB$method = "RB"

AggregatedAssesment = 
  bind_rows(assessment_aarhus_A,
            assessment_aarhus_A_BCDE,
            assessment_aarhus_A_B,
            assessment_aarhus_A_C,
            assessment_aarhus_A_D,
            assessment_aarhus_A_E)

plot_lp_performance(AggregatedAssesment, metric = "auroc", err = "sd")


end_time = Sys.time()
end_time - start_time