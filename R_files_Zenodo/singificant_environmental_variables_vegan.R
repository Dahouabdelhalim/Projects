library(tidyverse)
library(qiime2R)
library(phyloseq)
library(plyr)
library(ggplot2)
library(metacoder)
library(taxa)
library(dplyr)

#load data
#File location may need to be changed here

metadata<-read_tsv("MappingFile_mangue.csv")
feature_table<-read_qza("table-dada2.qza")
info_data <-feature_table$data
temptaxa <- read_qza('taxonomy.qza')
temptax<-temptaxa$data
rooted_tree<- read_qza("rooted_tree.qza")
root_tree <- rooted_tree$data
taxa_counts_transposed <- read.delim("taxa_counts_transposed.tab")
env_variables <- read.delim("env_variables.csv")

#Code for Figure 5
#Fig. 5. Significant correlation of salinity and organic matter with prokaryotic communities.
library(vegan)
my_varechem = env_variables
my_varespec = taxa_counts_transposed
rownames(my_varespec) <- my_varespec$OTU_ID
my_varespec <- my_varespec[-c(1)]
#Using Mantel
Salinity_DM <- vegdist(env_variables[,3])
Temperature_DM <- vegdist(env_variables[,4])
OrganicMatter_DM <- vegdist(env_variables[,5])
taxa_dist <- vegdist(my_varespec) # Creates distance matrix of the OTUs
taxa_salinity_Mantel <- mantel(taxa_dist, Salinity_DM, permutations = 999)
taxa_salinity_Mantel
taxa_temp_Mantel <- mantel(taxa_dist, Temperature_DM, permutations = 999)
taxa_temp_Mantel
taxa_OM_Mantel <- mantel(taxa_dist, OrganicMatter_DM, permutations = 999)
taxa_OM_Mantel

#Using NMDS, envfit
pdf('NMDS_taxa.pdf')
ord <- metaMDS((my_varespec), try=1000, k = 3)
(fit <- envfit(ord, my_varechem, perm = 999))
priSite <- diversity(my_varespec, index = "invsimpson", MARGIN = 1)
plot(ord)
orditorp(ord, display = "sites", priority = priSite, scaling = 3,
         col = "blue", cex = 1, pch = 19)
scores(fit, "vectors")
plot(fit)
plot(fit, p.max = 0.9, col = "blue")
plot(fit, p.max = 0.05, col = "red")
dev.off()
print(fit)

