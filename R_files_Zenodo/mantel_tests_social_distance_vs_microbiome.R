#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# Does social distance (grooming path length) predict microbial dissimilarity? 

rm(list=ls())
graphics.off()

library(igraph)
library(vegan)
library(agricolae)
library(ggplot2)

focal_dir <- "~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/"
setwd(focal_dir)
getwd()
source("social_data.R")
source("network_duration_function_dyad.R")

load('data/sifaka_braycurtis_marked_openref.RData')
load('data/sifaka_unifrac_weighted_marked_openref.RData')

#######################################################################
### Grooming network 6 months
#######################################################################

sifaka_groom_network6 <- net.foc.population(date.ref=as.Date("7/28/2012", "%m/%d/%Y"), duration=180, after=F, beh=groom, edgewidth=1000, direction="undirected", arrowsize=.5, labelcex=1, labeldist=.5, labels=T, margin=2, lay=99)

## stats 
sum(rowSums(sifaka_groom_network6$mat3)) # 3192 number of grooming interactions
sum(diag(sifaka_groom_network6$weight_matrix)/20) #2621 number of scans
length(sifaka_groom_network6$ind) #29
mean(sifaka_groom_network6$hours[sifaka_groom_network6$hours>0]/60) #22.99123
sd(sifaka_groom_network6$hours[sifaka_groom_network6$hours>0]/60) #6.914485
sum(sifaka_groom_network6$hours[sifaka_groom_network6$hours>0]/60) #436.8333

groom_paths_matrix6 <- as.matrix(sifaka_groom_network6$paths) #used inverse edge weights to calculate paths

#######################################################################
### Grooming network 12 months
#######################################################################

sifaka_groom_network12 <- net.foc.population(date.ref=as.Date("7/28/2012", "%m/%d/%Y"), duration=365, after=F, beh=groom, edgewidth=1000, direction="undirected", arrowsize=.5, labelcex=1, labeldist=.5, labels=F, margin=2, lay=99)

# stats
sum(rowSums(sifaka_groom_network12$mat3)) # 6972 number of grooming interactions
sum(diag(sifaka_groom_network12$weight_matrix)/20) #5126 number of scans
length(diag(sifaka_groom_network12$weight_matrix)[diag(sifaka_groom_network12$weight_matrix)>0])
sifaka_groom_network12$ind #33
mean(sifaka_groom_network12$hours[sifaka_groom_network12$hours>0]/60) #38.83333
sd(sifaka_groom_network12$hours[sifaka_groom_network12$hours>0]/60) #15.11455
sum(sifaka_groom_network12$hours[sifaka_groom_network12$hours>0]/60) #854.3333

groom_paths_matrix12 <- as.matrix(sifaka_groom_network12$paths) #used inverse edge weights to calculate paths

#######################################################################
# Proximity 6 months
#######################################################################
sifaka_prox_6 <- net.foc.population(date.ref=as.Date("7/28/2012", "%m/%d/%Y"), duration=180, after=F, beh=proximity_only, edgewidth=50, direction="undirected", arrowsize=.5, labelcex=1, labeldist=.5, labels=T, margin=2, lay=99)

# stats
sum(rowSums(sifaka_prox_6$mat3)) # 7482 number of proximity interactions
length(sifaka_prox_6$ind) #29

prox_paths_matrix6 <- as.matrix(sifaka_prox_6$paths) #used inverse edge weights to calculate paths

#######################################################################
# Proximity 12 months
#######################################################################
sifaka_prox <- net.foc.population(date.ref=as.Date("7/28/2012", "%m/%d/%Y"), duration=365, after=F, beh=proximity_only, edgewidth=50, direction="undirected", arrowsize=.5, labelcex=1, labeldist=.5, labels=F, margin=2, lay=99)

#stats
sum(rowSums(sifaka_prox$mat3)) # 13340 number of proximity interactions
sifaka_prox$ind #34
sum(diag(sifaka_prox$weight_matrix)/20) #5126 number of scans

prox_paths_matrix12 <- as.matrix(sifaka_prox$paths) #used inverse edge weights to calculate paths

#######################################################################
## Bray Curtis 
#######################################################################

bdist_marked_df  <- bdist_marked_df[order(names(bdist_marked_df)),order(names(bdist_marked_df))]
bdist_marked_mat <- as.matrix(bdist_marked_df)

#######################################################################
## Weighted Unifrac
#######################################################################

uwdist_marked_df <- uwdist_marked_df[order(names(uwdist_marked_df)),order(names(uwdist_marked_df))]
uwdist_marked_mat <- as.matrix(uwdist_marked_df)

#######################################################################
# Partial Mantel tests: 6 months Grooming vs microbiome dissimilarity while controlling for proximity
#######################################################################

#######################################################################
#### Bray Curtis
#######################################################################
# make sure matrices have same individuals
prox_paths_matrix <- prox_paths_matrix6[-7,-7]
groom_paths_matrix <- groom_paths_matrix6[-7,-7]
bdist_mat <- bdist_marked_mat[-c(1,2,3,4,33,34,35), -c(1,2,3,4,33,34,35)]

mantel.partial(groom_paths_matrix, bdist_mat, prox_paths_matrix, method="kendall", permutations=1000)
#sig
#######################################################################
#### Weighted Unifrac
#######################################################################

# make sure matrices have same individuals
prox_paths_matrix <- prox_paths_matrix6[-7,-7]
groom_paths_matrix <- groom_paths_matrix6[-7,-7]
uwdist_mat <- uwdist_marked_mat[-c(1,2,3,4,33,34,35), -c(1,2,3,4,33,34,35)]

mantel.partial(groom_paths_matrix, uwdist_mat, prox_paths_matrix, method="kendall", permutations=1000)
#sig
#######################################################################
# Partial Mantel tests: 12 months Grooming vs BC dissimilarity while controlling for proximity
#######################################################################

#######################################################################
#### Bray Curtis
#######################################################################
prox_paths_matrix <- prox_paths_matrix12[-c(3,7,10,25,34), -c(3,7,10,25,34)]
groom_paths_matrix <- groom_paths_matrix12[-c(3,7,10,25), -c(3,7,10,25)]
bdist_mat <- bdist_marked_mat[-c(1,2,4,33,34,35), -c(1,2,4,33,34,35)]

mantel.partial(groom_paths_matrix, bdist_mat, prox_paths_matrix, method="kendall", permutations=1000)
#NS
#######################################################################
#### Weighted Unifrac
#######################################################################
prox_paths_matrix <- prox_paths_matrix12[-c(3,7,10,25,34), -c(3,7,10,25,34)]
groom_paths_matrix <- groom_paths_matrix12[-c(3,7,10,25), -c(3,7,10,25)]
uwdist_mat <- uwdist_marked_mat[-c(1,2,4,33,34,35), -c(1,2,4,33,34,35)]

mantel.partial(groom_paths_matrix, uwdist_mat, prox_paths_matrix, method="kendall", permutations=1000)
#NS