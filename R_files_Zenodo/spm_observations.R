# Preparations
# setwd('C:/Users/Fava/Desktop/Dars/spm-mux/')

remove(list = ls())
setwd("~/GitHub/utlp-spm-mux")
library(igraph)
library(LinkPrediction)
library(lattice)
source("form_net.R")
source("prune_recover_custom.R")
source("lp_spm_mux_coeffs.R")
source('structural_similarity.R')
# physicians_l2
# physicians_l1
# 'celegans_l1'
# 'celegans_l3'
network_layer1_name <- 'LondonTransport_l1'
network_layer2_name <- 'aarhus_l3'

net_l1_fpath <- paste0('mux-data/', network_layer1_name, '.txt')
net_l2_fpath <- paste0('mux-data/', network_layer2_name, '.txt')

# graphT_table <- read.table(paste0('mux-data/', network_layer1_name, '.txt'))
# gT <- graph.data.frame(graphT_table[, c(1, 2)], directed = FALSE)
# gT <- simplify(gT, remove.multiple = TRUE, remove.loops = TRUE)
gT <- graph_from_adjacency_matrix(form_net(net_l1_fpath), mode = "undirected", diag = FALSE)

# graphA_table <- read.table(paste0('mux-data/', network_layer2_name, '.txt'))
# gA <- graph.data.frame(graphA_table[, c(1, 2)], directed = FALSE)
# global_gA <<- simplify(gA, remove.multiple = TRUE, remove.loops = TRUE)
global_gA <<- graph_from_adjacency_matrix(form_net(net_l2_fpath), mode = "undirected", diag = FALSE)


# global_gA <<-erdos.renyi.game(69, 0.96, type = c("gnp"), directed = FALSE,
#                  loops = FALSE)

# erdos.renyi.game(n, p.or.m, type = c("gnp", "gnm"), directed = FALSE,
#                  loops = FALSE, ...)
# g <- jazz_collab

# ADJ <- as_adjacency_matrix(jazz_collab, type = "both")
# XLX <- eigen(ADJ, symmetric = T)

# Calc structural consistency
# sigma_c <- structural_consistency(gT)
# print(paste('Structural Consistency of', network_layer1_name, ':', mean(sigma_c)))

# source('spm_correlation.R')
# # spm_correlation(gT)
# mean(replicate(100, spm_correlation(gT)))
# 
# # 
global_mux_pH <<- 0.1
# source('spm_correlation_mux.R')
# mean(replicate(100, spm_correlation_mux(gT)))

global_k <<- 10


# NCGT = nrow(as_adjacency_matrix(gT, type = "both"))
# 
# print(NCGT)
# 
# DENSITY_GT = sum(as_adjacency_matrix(gT, type = "both"))/(NCGT*NCGT-NCGT)
# 
# print(DENSITY_GT)
# 
# NCGA = nrow(as_adjacency_matrix(global_gA, type = "both"))
# 
# print(NCGA)
# 
# DENSITY_GA = sum(as_adjacency_matrix(global_gA, type = "both"))/(NCGA*NCGA-NCGA)
# 
# print(DENSITY_GA)

# gT <<-erdos.renyi.game(nrow(as_adjacency_matrix(gT, type = "both")), DENSITY_GT, type = c("gnp"), directed = FALSE, loops = FALSE)
# global_gA <<-erdos.renyi.game(nrow(as_adjacency_matrix(gT, type = "both")), DENSITY_GA, type = c("gnp"), directed = FALSE, loops = FALSE)

# SS <- structural_similarity(gT)


# rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "Lab")


# levelplot(t(SS))


# h = levelplot(t(SS),col.regions=rgb.palette(120), scales=list(x=list(cex=.3), y=list(cex=.3)))

# CE_l1_l2_m=0.3455233
# CE_l1_l3_m=0.3302769
# CE_l2_l3_m=0.4726207
# CE_2R_l3_m=0.2471570
# CE_2R_3R_m=0.2419014

# CE_l1_l2_fixed_m=0.3455233
# CE_l1_l3_fixed_m=0.3302769
# CE_l2_l3_fixed_m=0.4726207
# CE_2R_l3_fixed_m=0.2471570
# CE_2R_3R_fixed_m=0.2419014

# structural_similarity(gT)
# Structural_Similarity_arxiv_l1_l2_k=10_measure=0.6599148
source('lp_rnd.R')
# source('lp_mux_spm.R')
source('lp_spm_custom.R')
# assessment <- prune_recover_custom(gT, "lp_spm_custom", "lp_spm_mux_coeffs")
# assessment <- prune_recover(gT, "lp_mux_spm", "lp_spm_custom")
# source("prune_recover_custom.R")
assessment <- prune_recover_custom(gT, "lp_spm_custom")
# metric: recall_at_k, aupr, auroc, avg_prec
plot_lp_performance(assessment, metric = "auroc", err = "sd")

