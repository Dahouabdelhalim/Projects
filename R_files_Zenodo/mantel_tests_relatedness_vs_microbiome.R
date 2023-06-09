#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# Does relatedness predict microbial distance between samples?

rm(list=ls())
graphics.off()

library(igraph)
library(vegan)
library(phyloseq)

focal_dir <- "~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/"
setwd(focal_dir)
getwd()

load('data/sifaka_trimmed_phyloseq_markedonly_normalized_openref.RData')
load('data/sifaka_braycurtis_marked_openref.RData') #bdist_marked_df
load('data/sifaka_unifrac_unweighted_marked_openref.RData') #udist_marked_mat
load('data/sifaka_unifrac_weighted_marked_openref.RData') #uwdsit_marked_mat

#######################################################################
# RELATEDNESS & BACTERIA COMMUNITY DISTANCE
#######################################################################
# microbiome matrices
bdist_marked_mat <- as.matrix(bdist_marked_df)
udist_marked_mat <- as.matrix(udist_marked_df)
uwdist_marked_mat <- as.matrix(uwdist_marked_df)

# social group matrix
obs_social_group <- read.csv("data/social_group_matrix_2012.csv")
obs_social_group <- obs_social_group[,-1]
rownames(obs_social_group) <- colnames(obs_social_group)
obs_social_group <- as.matrix(obs_social_group)

#relatedness matrix
relatedness <- read.csv("data/sifaka_relatedness.csv")
relatedness <- relatedness[,-1]
rownames(relatedness) <- colnames(relatedness)
relatedness <- as.matrix(relatedness)

#######################################################################
### Relatedness and Microbiome Distance
#######################################################################

mantel(relatedness, bdist_marked_mat, method="kendall")

mantel(relatedness, udist_marked_mat, method="kendall")

mantel(relatedness, uwdist_marked_mat, method="kendall")

mantel.partial(obs_social_group,bdist_marked_mat,relatedness,method="kendall") #control for relatedness

mantel.partial(relatedness,bdist_marked_mat,obs_social_group,method="kendall") #control for social group affiliation

mantel(obs_social_group, relatedness, method="kendall") #individuals in the same group tend to be related

