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

#Code for Figure 8
#Fig 8. Phylogenetic tree 


#Fig 8. Phylogenetic tree of life at the family level.
feature_taxonomy <- read.delim("feature_taxonomy.tab")
obj <- parse_tax_data(feature_taxonomy,
                      class_cols = "taxa", # the column that contains taxonomic information
                      class_sep = ";", # The character used to separate taxa in the classification
                      class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                      class_key = c(tax_rank = "info", # A key describing each regex capture group
                                    tax_name = "taxon_name"))
obj$data$tax_data <- zero_low_counts(obj, data = "tax_data", min_count = 100)
obj <- filter_obs(obj, target = "tax_data", ! no_reads, drop_taxa = TRUE)
set.seed(1) # This makes the plot appear the same each time it is run 
heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          output_file = "heat_tree.pdf") # Saves the plot as a pdf file)
