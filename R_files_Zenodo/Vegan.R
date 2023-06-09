#We are now moving to R to handle the files generated above to generate additional figures.

#```{r library load}
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(plyr)
library(ggplot2)
library(metacoder)
library(taxa)
library(dplyr)
#File location may need to be changed here
setwd('C:/Gresham/Project_Gravimondo/BGS_final/aldex_results/')

metadata<-read_tsv("C:/Gresham/Project_Gravimondo/BGS_final/metadata/MappingFile_mangue.csv")
feature_table<-read_qza("C:/Gresham/Project_Gravimondo/BGS_final/qiime_results/table-dada2.qza")
info_data <-feature_table$data
temptaxa <- read_qza('C:/Gresham/Project_Gravimondo/BGS_final/qiime_results/taxonomy.qza')
temptax<-temptaxa$data
rooted_tree<- read_qza("C:/Gresham/Project_Gravimondo/BGS_final/qiime_results/rooted-tree.qza")
root_tree <- rooted_tree$data
#unzip tabe-dada2.gza
#biom convert -i feature-table.biom -o taxa_table.tsv --to-tsv
#taxa_counts <- read.delim("C:/Gresham/Project_Gravimondo/BGS_final/qiime_results/feature_taxonomy.tab")
taxa_counts <- read.delim("C:/Gresham/Project_Gravimondo/BGS_final/qiime_results/taxa_table.tsv", sep='\\t')
taxa_counts_transposed <- as.data.frame(t(taxa_counts))
env_variables <- read.delim("C:/Gresham/Project_Gravimondo/BGS_final/metadata/env_variables.csv")
#```

#Fig. 5. Significant correlation of salinity and organic matter with prokaryotic communities.

#```{r figure 5}
#Fig. 5. Significant correlation of salinity and organic matter with prokaryotic communities.
library(vegan)
my_varechem = env_variables
rownames(my_varechem)<-my_varechem$Sample
my_varechem<-my_varechem[-c(1)]

my_varespec = taxa_counts
rownames(my_varespec) <- my_varespec$taxa
my_varespec <- my_varespec[-c(1)]

my_varespec <- as.data.frame(t(my_varespec))

#Using NMDS, envfit
pdf('NMDS_taxa.pdf')
ord <- metaMDS(t(taxa_counts), engine = c("isoMDS"), try=1000, k = 3)
fit <- envfit(ord, my_varechem, perm = 999)
priSite <- diversity(taxa_counts, index = "invsimpson", MARGIN = 1)
plot(ord)
orditorp(ord, display = "sites", priority = priSite, scaling = 3,
         col = "blue", cex = 1, pch = 19)
scores(fit, "vectors")
plot(fit)
plot(fit, p.max = 0.4, col = "blue")
plot(fit, p.max = 0.05, col = "red")
dev.off()
print(fit)

###
Otu_veg <- t(info_data) # Samples become rows and OTUs are columns
Envdt <- env_variables[,2:5]

vare.cap <- capscale(Otu_veg ~ (Salinity + Organic_Matter + Temperature + Water_content) , Envdt,
                     dist="bray", sqrt.dist= TRUE,  metaMDS = TRUE)
vare.cap
pdf('vegan_capscale.pdf')
plot(vare.cap)
dev.off()
anova(vare.cap)
