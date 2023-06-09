########## packges load ###################################################################################
library(vegan)
library(car)
library (DESeq2); packageVersion("DESeq2")
library(Rfast) #to take row maximums
library (ggthemes)
library(data.table); packageVersion("data.table")
library(adaptiveGPCA)
library(ggrepel)
library(phyloseq)

########## Data ###########################################################################################
data.in <- readRDS(file="Crane_CleanDataNewAproach.rds")
#data.in1 <- readRDS(file="Crane_CleanDataNewAproachFinal.rds")
# Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(data.in))

# Histogram of sample read counts, from http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
#Get the read numbers for each sample
read.number <- sample_sums(data.in)
########## Rearefy to 8000 #################################################################################
#Keep only samples with at least 8000 reads (before rerefaction)
filtered.data.pared.down.fall.2017.8000reads = prune_samples(sample_sums(data.in)>=8000, data.in)
filtered.data.pared.down.fall.2017.8000reads = filter_taxa(filtered.data.pared.down.fall.2017.8000reads, function(x) sum(x) > 0, TRUE)
# there are 14 with less then 8000 reads (7.5% of samples lost), total of 4666 ASVs

## Keeps 184 samples ## Kepps 5326 taxa
titles_reads_8000=c("sum","mean","st.dev", "median", "min","max")
values_reads_8000=print(c(sum(sample_sums(filtered.data.pared.down.fall.2017.8000reads)),
                          mean(sample_sums(filtered.data.pared.down.fall.2017.8000reads)),
                          sd(sample_sums(filtered.data.pared.down.fall.2017.8000reads)), 
                          median(sample_sums(filtered.data.pared.down.fall.2017.8000reads)),
                          min(sample_sums(filtered.data.pared.down.fall.2017.8000reads)),
                          max(sample_sums(filtered.data.pared.down.fall.2017.8000reads))))
reads_data_8000=data.frame(titles_reads_8000,values_reads_8000)
reads_data_8000
#titles_reads_8000 values_reads_8000
#1               sum       2717202.000
#2              mean         14767.402
#3            st.dev          3956.448
#4            median         14226.500
#5               min          8158.000
#6               max         25670.000
#rarefy_even_depth downsamples/normalizes all samples to the same depth and prunes OTUs that disappear from all samples as a result.
filtered.data.pared.down.fall.2017.8000reads.rarefied <- rarefy_even_depth(filtered.data.pared.down.fall.2017.8000reads, rngseed = 999)
filtered.data.pared.down.fall.2017.8000reads.rarefied #1038 OTUs removed
#otu_table()   OTU Table:         [4288 taxa and 184 samples]

# SAVE The clean data
saveRDS (filtered.data.pared.down.fall.2017.8000reads.rarefied, file="Crane_CleanData_rarefied_8000.rds")