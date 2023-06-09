#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# Individual differences in diet: Shannon's evenness in plant parts and plant species consumed
# Between-group differences in plant parts and plant species consumed
# Does dietary distance predict microbial dissimilarity? 
rm(list=ls())
graphics.off()

home_dir <- "~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/"
setwd(home_dir)
getwd()

library("coin")
library("RColorBrewer")
library("lattice")
library("latticeExtra")
library("vegan")
library("ggplot2")
library("dplyr")
library("scales")
library("tidyr")

load("data/diet_diversity_df.RData")
diet_diversity_df

## relative abundance table of plant species consumed
load("data/sifaka_plant_spp_rel_abund.RData") #spp_rel
spp_rel

## relative abundance table of plant parts consumed
load("data/sifaka_plant_parts_rel_abund.RData") #food_rel
food_rel

## diet distance (plant spp)
bdist_diet_spp <- vegdist(spp_rel[,-14], method="bray")

## diet distance (plant parts)
food_part_bdist <-vegdist(food_rel[,-8], method="bray")

######################################################
#### Individual differences in diet
######################################################

load("data/full_table_for_poisson_glmm.RData") #poisson_glmm
colnames(poisson_glmm)
head(poisson_glmm)

diet_richness_glmm <- merge(poisson_glmm, diet_diversity_df, by.y="Individual")
diet_richness_glmm <- droplevels(diet_richness_glmm)
colnames(diet_richness_glmm)

save(diet_richness_glmm, file="data/diet_diversity_table_for_poisson_glmm.RData")

######################################################
#### Between-group differences in plant parts consumed
######################################################

######################################################
# Fruit
######################################################

kw_fruit <- kruskal_test(Fruit ~ Group, data=food_rel, conf.int=TRUE, distribution= approximate(B = 10000))
kw_fruit #NS
kw_fruit2 <- independence_test(Fruit~Group, data=food_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_fruit2, method="single-step"))
pvalue(kw_fruit2)

######################################################
# Mature Leaves
######################################################

kw_mat2 <- independence_test(Mature~Group, data=food_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_mat2, method="single-step"))
pvalue(kw_mat2)

######################################################
# Young Leaves
######################################################

kw_young2 <- independence_test(Young~Group, data=food_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_young2, method="single-step"))
pvalue(kw_young2)

######################################################
# Seeds
######################################################

kw_seeds2 <- independence_test(Seeds~Group, data=food_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_seeds2, method="single-step"))
pvalue(kw_seeds2)

######################################################
# Bark
######################################################

kw_bark2 <- independence_test(Bark~Group, data=food_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_bark2, method="single-step"))
pvalue(kw_bark2)

######################################################
# Flowers
######################################################

kw_flowers2 <- independence_test(Flowers~Group, data=food_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_flowers2, method="single-step"))
pvalue(kw_flowers2)

######################################################
# Stems
######################################################

kw_stems2 <- independence_test(Stems~Group, data=food_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_stems2, method="single-step"))
pvalue(kw_stems2)

######################################################
#### Between-group differences in common plant species consumed
######################################################
kw_met2 <- independence_test(Metampototsy~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_met2, method="single-step"))
pvalue(kw_met2)

kw_alim2 <- independence_test(AlimboroMahalao~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_alim2, method="single-step"))
pvalue(kw_alim2)

kw_tan2 <- independence_test(Tandredretsy~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_tan2, method="single-step"))
pvalue(kw_tan2)

kw_haz2 <- independence_test(Hazomboenga~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_haz2, method="single-step"))
pvalue(kw_haz2)

kw_alimb2 <- independence_test(Alimboro~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_alimb2, method="single-step"))
pvalue(kw_alimb2)

kw_aman2 <- independence_test(Amaninomby~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_aman2, method="single-step"))
pvalue(kw_alimb2)

kw_remo2 <- independence_test(Remotiny~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_remo2, method="single-step"))
pvalue(kw_remo2)

kw_MB2 <- independence_test(Magnary_Baomby~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_MB2, method="single-step"))
pvalue(kw_MB2)

kw_bag2 <- independence_test(Bagnaky~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_bag2, method="single-step"))
pvalue(kw_bag2)

kw_man2 <- independence_test(Manjakabetany ~ Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_man2, method="single-step"))
pvalue(kw_man2)

kw_sako2 <- independence_test(Sakoambanditsy~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_sako2, method="single-step"))
pvalue(kw_sako2)

kw_lova2 <- independence_test(Lovainjafy~Group, data=spp_rel,distribution = approximate(B = 10000), ytrafo = function(data)trafo(data, numeric_trafo = rank_trafo),xtrafo = mcp_trafo(Group = "Tukey"))
print(pvalue(kw_lova2, method="single-step"))
pvalue(kw_lova2)

#######################################################################
## Dietary Distance Analysis
#######################################################################

library(igraph)
library(vegan)
library(phyloseq)

source("social_data.R")
source('network_duration_function_dyad.R')

load('data/sifaka_trimmed_phyloseq_markedonly_normalized_openref.RData')
load('data/sifaka_braycurtis_marked_openref.RData')
load('data/sifaka_unifrac_weighted_marked_openref.RData')

#######################################################################
# create edgelist of bray curtis distances between dyads
bdist_marked_df  <- bdist_marked_df[order(names(bdist_marked_df)),order(names(bdist_marked_df))]
bdist_marked_mat <- as.matrix(bdist_marked_df)
colnames(bdist_marked_mat)
bdist_trimmed <- bdist_marked_mat[-c(1:4,33:35),-c(1:4,33:35)]

#######################################################################
# weighted unifrac
uwdist_marked_df <- uwdist_marked_df[order(names(uwdist_marked_df)),order(names(uwdist_marked_df))]
uwdist_marked_mat <- as.matrix(uwdist_marked_df)
uwdist_trimmed <- uwdist_marked_mat[-c(1:4,33:35),-c(1:4,33:35)]

#######################################################################
### Proximity network
#######################################################################

sifaka_prox_network <- net.foc.population(date.ref=as.Date("7/28/2012", "%m/%d/%Y"), duration=180, after=F, beh=proximity_only, edgewidth=50, direction="undirected", arrowsize=.5, labelcex=1, labeldist=.5, labels=T, margin=2, lay=99)
prox_paths_matrix <- as.matrix(sifaka_prox_network$paths)
prox_paths_matrix_trimmed <- prox_paths_matrix[-7,-7] 
prox_paths_matrix_trimmed <- prox_paths_matrix_trimmed[order(rownames(prox_paths_matrix_trimmed)),order(colnames(prox_paths_matrix_trimmed))]
#######################################################################
### Grooming network
#######################################################################

sifaka_groom_network <- net.foc.population(date.ref=as.Date("7/28/2012", "%m/%d/%Y"), duration=180, after=F, beh=groom, edgewidth=1000, direction="undirected", arrowsize=.5, labelcex=1, labeldist=.5, labels=T, margin=2, lay=99)
groom_paths_matrix <- as.matrix(sifaka_groom_network$paths)
groom_paths_matrix_trimmed <- groom_paths_matrix[-7,-7]
groom_paths_matrix_trimmed <- groom_paths_matrix_trimmed[order(rownames(groom_paths_matrix_trimmed)),order(colnames(groom_paths_matrix_trimmed))]

#######################################################################
### Diet matrix
#######################################################################
diet_part_dist_mat <- as.matrix(food_part_bdist)
diet_part_dist_mat <- diet_part_dist_mat[order(rownames(diet_part_dist_mat)),order(colnames(diet_part_dist_mat))]

diet_spp_dist_mat <- as.matrix(bdist_diet_spp)
diet_spp_dist_mat <- diet_spp_dist_mat[order(rownames(diet_spp_dist_mat)),order(colnames(diet_spp_dist_mat))]
#######################################################################
### Trimming other matrices
#######################################################################

groom_paths_matrix_trimmed_diet <- groom_paths_matrix_trimmed[-c(2,7,8,10,13,15,16,17,22,23,25,26),-c(2,7,8,10,13,15,16,17,22,23,25,26)]
prox_paths_matrix_trimmed_diet <- prox_paths_matrix_trimmed[-c(2,7,8,10,13,15,16,17,22,23,25,26),-c(2,7,8,10,13,15,16,17,22,23,25,26)]
bdist_trimmed_diet <- bdist_trimmed[-c(2,7,8,10,13,15,16,17,22,23,25,26),-c(2,7,8,10,13,15,16,17,22,23,25,26)]
uwdist_trimmed_diet <- uwdist_trimmed[-c(2,7,8,10,13,15,16,17,22,23,25,26),-c(2,7,8,10,13,15,16,17,22,23,25,26)]

#######################################################################
### Does diet (tree species) predict microbiome dissimilarity? NO
#######################################################################
mantel(diet_spp_dist_mat, bdist_trimmed_diet, method="kendall")

mantel(diet_spp_dist_mat, uwdist_trimmed_diet, method="kendall")

#######################################################################
### Does diet (plant parts) predict microbiome dissimilarity? NO
#######################################################################
mantel(diet_part_dist_mat, bdist_trimmed_diet, method="kendall")

mantel(diet_part_dist_mat, uwdist_trimmed_diet, method="kendall")

#######################################################################
### Does proximity predict diet? YES (for plant species)
#######################################################################

### proximity vs. diet spp distance 
mantel(prox_paths_matrix_trimmed_diet, diet_spp_dist_mat, method="kendall")
#sig
### proximity vs. diet part distance
mantel(prox_paths_matrix_trimmed_diet, diet_part_dist_mat, method="kendall")
#NS
#######################################################################
### Does proximity predict microbiome dissimilarity, while controlling for diet (tree species?) YES
#######################################################################
mantel.partial(prox_paths_matrix_trimmed_diet, bdist_trimmed_diet, diet_spp_dist_mat, method="kendall")
# significant 

mantel.partial(prox_paths_matrix_trimmed_diet, uwdist_trimmed_diet, diet_spp_dist_mat, method="kendall")
# significant

#######################################################################
### Does proximity predict microbiome dissimilarity, while controlling for diet (plant parts)? YES
#######################################################################
mantel.partial(prox_paths_matrix_trimmed_diet, bdist_trimmed_diet, diet_part_dist_mat, method="kendall")
#sig

mantel.partial(prox_paths_matrix_trimmed_diet, uwdist_trimmed_diet, diet_spp_dist_mat, method="kendall")
#sig
#######################################################################
### Does grooming predict diet? Yes, for tree species but no for plant parts
#######################################################################
mantel(groom_paths_matrix_trimmed_diet, diet_spp_dist_mat, method="kendall")
#Sig
mantel(groom_paths_matrix_trimmed_diet, diet_part_dist_mat, method="kendall")
#NS
#######################################################################
### Does grooming predict microbiome dissimilarity, while controlling for diet (tree species?) YES
#######################################################################
mantel.partial(groom_paths_matrix_trimmed_diet, bdist_trimmed_diet, diet_spp_dist_mat, method="kendall")
#sig

mantel.partial(groom_paths_matrix_trimmed_diet, uwdist_trimmed_diet, diet_spp_dist_mat, method="kendall")
#sig

#######################################################################
### Does grooming predict microbiome dissimilarity, while controlling for diet (plant parts)? YES
#######################################################################
mantel.partial(groom_paths_matrix_trimmed_diet, bdist_trimmed_diet, diet_part_dist_mat, method="kendall")
#sig

mantel.partial(groom_paths_matrix_trimmed_diet, uwdist_trimmed_diet, diet_part_dist_mat, method="kendall")
#sig