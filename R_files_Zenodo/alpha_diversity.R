#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# Rarefy dataset 100 times to create summary table of rarefied richness values
# Determine differences in richness across social groups

rm(list=ls())
graphics.off()

library("phyloseq")
library("grid")
library("vegan")
library("DESeq2")
library("ggplot2")
library("agricolae")
setwd("~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/")
getwd()

load("data/sifaka_allotus_phyloseq_openref.RData") #kmnp_data (raw counts dataset for alpha diversity; singletons removed)
sample_data(kmnp_data)
#######################################################################
## Rarefaction and Alpha Diversity
#######################################################################

# Initialize matrices to store richness and evenness estimates
richness = matrix(nrow=47,ncol=100)
row.names(richness) <- sample_names(kmnp_data)
evenness = matrix(nrow=47,ncol=100)
row.names(evenness) <- sample_names(kmnp_data)
chao1 = matrix(nrow=47,ncol=100)
row.names(chao1) <- sample_names(kmnp_data)

set.seed(42)
# For 100 replications, rarefy the OTU table to 39533 reads and store the richness and evenness estimates. The default for the rarefy_even_depth command is to pick with replacement so I set it to false. Picking without replacement is more computationally intensive 
for (i in 1:100) {
  r=rarefy_even_depth(kmnp_data,sample.size= 39533,verbose=FALSE,replace = FALSE)
  rich= as.numeric(as.matrix(estimate_richness(r,measures="Observed")))
  richness[,i]=rich
  even=as.numeric(as.matrix(estimate_richness(r,measures="Shannon")))
  evenness[,i]=even
  chao=estimate_richness(r,measures="Chao1")
  chao_1=as.numeric(as.matrix(chao[,1]))
  chao1[,i]=chao_1
  print(i)
}

# Create a new matrix to hold the means and standard deviations of all the richness estimates
rich.stats = matrix(nrow=47,ncol=2)
rich.stats[,1] = apply(richness,1,mean)
rich.stats[,2] = apply(richness,1,sd)
rich.stats = data.frame(row.names(richness),rich.stats)
colnames(rich.stats) = c("samples","mean","sd")

# Create a new matrix to hold the means and standard deviations of the evenness estimates
even.stats = matrix(nrow=47,ncol=2)
even.stats[,1] = apply(evenness,1,mean)
even.stats[,2] = apply(evenness,1,sd)
even.stats = data.frame(row.names(evenness),even.stats)
colnames(even.stats) = c("samples","mean","sd")

chao.stats = matrix(nrow=47,ncol=2)
chao.stats[,1] = apply(chao1,1,mean)
chao.stats[,2] = apply(chao1,1,sd)
chao.stats = data.frame(row.names(richness),chao.stats)
colnames(chao.stats) = c("samples","mean","sd")

Gp = data.frame(SampleID=sample_data(kmnp_data)$X.SampleID,Name=sample_data(kmnp_data)$Name,Group=sample_data(kmnp_data)$Group, Natal=sample_data(kmnp_data)$Natal, Sex=sample_data(kmnp_data)$Sex, Immigrant=sample_data(kmnp_data)$Immigrant, Age=sample_data(kmnp_data)$Age, Chest=sample_data(kmnp_data)$Chest, Age_Rank=sample_data(kmnp_data)$Age_Rank)
head(Gp)

#re-order group levels
Gp$Group <- factor(Gp$Group, levels = c("I","II","III","IV","V","VI","Camp"))

# Rename the headers
colnames(rich.stats)[1] <- "SampleID"     
colnames(even.stats)[1] <- "SampleID"     
colnames(chao.stats)[1] <- "SampleID"     
rich.stats2 = merge(rich.stats, Gp,by="SampleID")
even.stats2 = merge(even.stats, Gp,by="SampleID")
chao.stats2 = merge(chao.stats, Gp,by="SampleID")

rich.stats2$Age <- as.factor(rich.stats2$Age)
even.stats2$Age <- as.factor(even.stats2$Age)
chao.stats2$Age <- as.factor(chao.stats2$Age)
rich.stats2$Sex <- as.factor(rich.stats2$Sex)
even.stats2$Sex <- as.factor(even.stats2$Sex)
chao.stats2$Sex <- as.factor(chao.stats2$Sex)
rich.stats2$Chest <- as.factor(rich.stats2$Chest)
even.stats2$Chest <- as.factor(even.stats2$Chest)
chao.stats2$Chest <- as.factor(chao.stats2$Chest)
rich.stats2$Dominance <- as.factor(rich.stats2$Age_Rank)
even.stats2$Dominance <- as.factor(even.stats2$Age_Rank)
chao.stats2$Dominance <- as.factor(chao.stats2$Age_Rank)

mean(rich.stats2$mean) 
min(rich.stats2$mean) 
max(rich.stats2$mean)
sd(rich.stats2$mean)
mean(even.stats2$mean) 
sd(even.stats2$mean) 
mean(chao.stats2$mean) 
sd(chao.stats2$mean) 

save(rich.stats2, file="data/sifaka_alpha_diversity_richstats2_table_openref.RData")
save(even.stats2, file="data/sifaka_alpha_diversity_evenstats2_table_openref.RData")
save(chao.stats2, file="data/sifaka_alpha_diversity_chaostats2_table_openref.RData")

load(file="data/sifaka_alpha_diversity_richstats2_table_openref.RData")
load(file="data/sifaka_alpha_diversity_evenstats2_table_openref.RData")
load(file="data/sifaka_alpha_diversity_chaostats2_table_openref.RData")

rich <- kruskal(rich.stats2$mean, rich.stats2$Group, p.adj="BH")
rich

even <- kruskal(even.stats2$mean, even.stats2$Group, p.adj="BH")
even

chao <- kruskal(chao.stats2$mean, chao.stats2$Group, p.adj="BH")
chao

#######################################################################
# Unique OTUs per sample 
########################################################################
# Normalize raw data set so that sequencing depth is accounted for
kmnp.deseq <- phyloseq_to_deseq2(kmnp_data, ~Group) #convert phyloseq object to deseq2 object

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(kmnp.deseq), 1, gm_mean)
kmnp.deseq= estimateSizeFactors(kmnp.deseq, geoMeans = geoMeans)
kmnp.deseq = estimateDispersions(kmnp.deseq)
diagvst = getVarianceStabilizedData(kmnp.deseq)
dim(diagvst)
#[1] 33433    47
sifaka2 <- kmnp_data #rename kmnp_processed to sifaka2 so that we still have the original OTU counts for raw reads
normalizedCounts <- t( t(counts(kmnp.deseq)) / sizeFactors(kmnp.deseq) )
otu_table(kmnp_data) <- otu_table(normalizedCounts, taxa_are_rows = TRUE)
head(otu_table(sifaka2))
head(otu_table(kmnp_data))
kmnp_data <- prune_taxa(taxa_sums(kmnp_data)>0,kmnp_data)
otu_table(kmnp_data) <- round(otu_table(kmnp_data), digits=0)
head(otu_table(kmnp_data))
kmnp_processed <- kmnp_data
save(kmnp_processed, file="data/sifaka_processed_phyloseq_normalized_openref.RData") #normalized, all OTUs (excluding singletons)
save(sifaka2, file="data/sifaka_processed_phyloseq_raw_openref.RData") #raw reads (not normalized), all OTUs (excluding singletons); same as kmnp_data
get_taxa_unique(sifaka2, "Phylum")

# unique otus for normalized data (all OTUs)
otu_normalized <- estimate_richness(kmnp_processed,split=T,measures="Observed")
otu_normalized <- data.frame(SampleID=rownames(otu_normalized), Observed=otu_normalized)
otu_normalized  <- otu_normalized [order(otu_normalized $SampleID),]
otu_normalized <- otu_normalized[,-1]

#unique otus for raw data (all OTUs)
otu_raw <- estimate_richness(sifaka2,split=T,measures="Observed")
otu_raw <- data.frame(SampleID=rownames(otu_raw), Observed=otu_raw)
otu_raw <- otu_raw[order(otu_raw$SampleID),]
otu_raw <- otu_raw[,-1]
otu_rarefied <- rich.stats2$mean
otu_chao <- chao.stats2$mean

#######################################################################
# Unique Genera per sample

# normalized dataset (all OTUs)
# create genus factor (but this discards unassigned OTUs)
# Create a factor corresponding to the Genera
genfac = factor(tax_table(kmnp_processed)[, "Genus"])
# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(kmnp_processed), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(gentab)[, 1:10]

# To get number of non-zero genera per sample, sum the values that are above
# your threshold, in your case, 1.
observationThreshold = 1
normalized_genera_table <- data.frame(apply(gentab > observationThreshold, 2, sum))
colnames(normalized_genera_table) <-  "Genera"
normalized_genera_table <- data.frame(SampleID=rownames(normalized_genera_table), Genera=normalized_genera_table)
normalized_genera_table <- normalized_genera_table[order(normalized_genera_table$SampleID),]
normalized_genera_table <- normalized_genera_table[,-1]


#raw dataset (all OTUs)
genfac = factor(tax_table(sifaka2)[, "Genus"])
# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(sifaka2), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(gentab)[, 1:10]

# To get number of non-zero genera per sample, sum the values that are above
# your threshold, in your case, 1.
observationThreshold = 1
raw_genera_table <- data.frame(apply(gentab > observationThreshold, 2, sum))
colnames(raw_genera_table) <-  "Genera"
raw_genera_table <- data.frame(SampleID=rownames(raw_genera_table), Genera=raw_genera_table)
raw_genera_table <- raw_genera_table[order(raw_genera_table$SampleID),]
raw_genera_table <- raw_genera_table[,-1]

alpha_diversity_df <- data.frame(SampleID=rich.stats2$SampleID, Rarefied_OTUs = otu_rarefied, Norm_OTUs=otu_normalized, Raw_OTUs=otu_raw, Chao_OTU = otu_chao, Norm_Genera=normalized_genera_table, Raw_Genera=raw_genera_table)

save(alpha_diversity_df, file="data/alpha_diversity_table.RData")
