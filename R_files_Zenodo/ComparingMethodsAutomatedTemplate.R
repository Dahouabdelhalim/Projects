# In the Name of God

# This scripts show the case in which the ComparingMethodsCallerFunction is used

#### Preparations ####

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

source('LRM_Algorithm.R')
source('Masked_LRM_Algorithm.R')

source('ComparingMethodsCallerFunction.R')


#### Preparing Data Sets ####
network_layer1_name <- 'at_l1'
network_layer2_name <- 'at_l2'
network_layer3_name <- 'LondonTransport_l3'

net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')
net_l3_fpath <- paste0('mux-data/', network_layer3_name, '.txt')


# network_layer1_name <- 'aarhus_l1'
# network_layer2_name <- 'aarhus_l2'
# network_layer3_name <- 'aarhus_l3'
# network_layer4_name <- 'aarhus_l4'
# network_layer5_name <- 'aarhus_l5'
# 
# net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
# net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')
# net_l3_fpath <- paste0('mux-data/', network_layer3_name, '.txt')
# net_l4_fpath <- paste0('mux-data/', network_layer4_name, '.txt')
# net_l5_fpath <- paste0('mux-data/', network_layer5_name, '.txt')

# The main layer to be predicted is always layer A. 
# All other layers act as auxiliary 
gA <- graph_from_adjacency_matrix(form_net(net_l1_fpath), mode = "undirected", diag = FALSE)
gB <- graph_from_adjacency_matrix(form_net(net_l2_fpath), mode = "undirected", diag = FALSE)
gC <- graph_from_adjacency_matrix(form_net(net_l3_fpath), mode = "undirected", diag = FALSE)
# gD <- graph_from_adjacency_matrix(form_net(net_l4_fpath), mode = "undirected", diag = FALSE)
# gE <- graph_from_adjacency_matrix(form_net(net_l5_fpath), mode = "undirected", diag = FALSE)
# Randomized graphes
gRA <<-erdos.renyi.game(nrow(as_adjacency_matrix(gA, type = "both")), graph.density(gA), type = c("gnp"), directed = FALSE, loops = FALSE)
gRB <<-erdos.renyi.game(nrow(as_adjacency_matrix(gB, type = "both")), graph.density(gB), type = c("gnp"), directed = FALSE, loops = FALSE)
gRC <<-erdos.renyi.game(nrow(as_adjacency_matrix(gC, type = "both")), graph.density(gC), type = c("gnp"), directed = FALSE, loops = FALSE)

#### Analysis ####

Target_Layer_Graph = gA
permutatedgA = permute(gA, 69:1)
ReconstructionList_Graph = 
  list( 
    # list(gB), 
    # list(gC), 
    # list(gRB), 
    # list(gRC),
    list(gA)
    )
LegendList =
  list(
    # "LRM(Overground)",
    # "LRM(DLR)",
    # "leisure",
    # "LRM-rand(Overground)",
    "LRM(Self+Backup)"
  )

AggregatedAssesment = 
  ComparingMethodsCallerFunction(
    Target_Layer_Graph, 
    method = "LRM_Algorithm", 
    ReconstructionList_Graph = ReconstructionList_Graph, 
    LegendList,
    probes = seq(0.1, 0.9, 0.1), 
    epochs = 30, # Number of Train and Test Sets are calculated
    PowerSelector = NA, # Percent of power of eigenvalues to be used. Used by algegraic largest method. If powerselector has value, NumberOfLeadingEigenValues will be ignored. If number of eigenvalues is going to be specifically given, give NA as input for this arguement
    NumberOfLeadingEigenValues = NA, # if NA, all eigenvalues are used
    SelectorOfEigenValues = "LA" # common options include LA for algebraic largest, and LM for magnitude largest
    ) 

# Save an object to a file
saveRDS(AggregatedAssesment, file = "AggregatedAssesment_Air-Train-Air-k=NA-mirror-paper.rds")

plot_lp_performance(AggregatedAssesment, metric = "auroc", err = "sd")


end_time = Sys.time()
end_time - start_time