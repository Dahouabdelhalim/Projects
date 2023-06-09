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

#Code for Figure 2
#Fig 2. Taxonomic abundances from all sample sites.

```{r figure 2}
#Fig 2. Taxonomic abundances from all sample sites.
#### Create phyloseq object
MG_data<-phyloseq(
  otu_table(info_data, taxa_are_rows = T),
  phy_tree(root_tree),
  tax_table(as.data.frame(taxtable_p1) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), # transforming the taxonomy for compatibility with phyloseq
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample-id")))
### Subsample the most abundant phyla
phylum.sum = tapply(taxa_sums(MG_data), tax_table(MG_data)[, "Phylum"], sum, na.rm=TRUE) #objeto sรณ com 6 filos mais abundantes
top6phyla = names(sort(phylum.sum, TRUE))[1:6]
MG_6p_new = prune_taxa((tax_table(MG_data)[, "Phylum"] %in% top6phyla), MG_data)
### Exclude OTUs that account for less than 1 percent of total (for more taxonomy levels)
Abund_otus1 <- filter_taxa(MG_data, function(x) sum(x > total*0.01) > 0, TRUE)
### Taxonomy barplots
plot_bar(Abund_otus1, "Phylum", fill="Class") + geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
plot_bar(Abund_otus1, "Class", fill="Order") + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")
plot_bar(Abund_otus1, "Order", fill="Family") + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")+ theme(legend.position="bottom")


