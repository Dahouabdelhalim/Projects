#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# PERMANOVA analyses: which traits are associated with variation in microbial composition? 
# Group affiliation
# Age #individuals of known identity ('marked')
# Sex #individuals of known identity ('marked')
# Sex/dominance rank (adults in multi-male and/or multi-female groups only)

rm(list=ls())
graphics.off()

library("phyloseq")
library("vegan")
setwd("~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/") #set working directory to correct folder
getwd()

load("data/sifaka_trimmed_phyloseq_normalized_openref.RData") #kmnp_trimmed (normalized dataset for beta diversity analyses)
kmnp_trimmed
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 9037 taxa and 47 samples ]
# sample_data() Sample Data:       [ 47 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 9037 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 9037 tips and 9035 internal nodes ]
# refseq()      DNAStringSet:      [ 9037 reference sequences ]

sample_data(kmnp_trimmed)

load("data/sifaka_trimmed_phyloseq_markedonly_normalized_openref.RData") #kmnp_marked (normalized, but just including 'marked' individuals)
kmnp_marked
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8849 taxa and 35 samples ]
# sample_data() Sample Data:       [ 35 samples by 22 sample variables ]
# tax_table()   Taxonomy Table:    [ 8849 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8849 tips and 8847 internal nodes ]
# refseq()      DNAStringSet:      [ 8849 reference sequences ]

sample_data(kmnp_marked)
#######################################################################
## Adonis for Normalized OTU counts
#######################################################################

bdist <- phyloseq::distance(kmnp_trimmed,"bray")
uwdist <- phyloseq::distance(kmnp_trimmed,"unifrac", weighted=TRUE)
udist <- phyloseq::distance(kmnp_trimmed,"unifrac", weighted=FALSE)

adonis.group <- adonis(bdist~Group, as(sample_data(kmnp_trimmed),"data.frame"), permutations=999)
adonis.group

#######################################################################
# Weighted Unifrac

adonis.group <- adonis(udist~Group, as(sample_data(kmnp_trimmed),"data.frame"), permutations=999)
adonis.group

#######################################################################
# Unifrac

adonis.group <- adonis(udist~Group, as(sample_data(kmnp_trimmed),"data.frame"), permutations=999)
adonis.group

#######################################################################
## Adonis (marked individuals only)
#######################################################################

bdist_marked <- phyloseq::distance(kmnp_marked,"bray")
bdist_marked_mat <- as.matrix(bdist_marked)
colnames(bdist_marked_mat) <- sample_data(kmnp_marked)$Name
rownames(bdist_marked_mat) <- sample_data(kmnp_marked)$Name

uwdist_marked <- phyloseq::distance(kmnp_marked,"unifrac", weighted=TRUE)
uwdist_marked_mat <- as.matrix(uwdist_marked)
colnames(uwdist_marked_mat) <- sample_data(kmnp_marked)$Name
rownames(uwdist_marked_mat) <- sample_data(kmnp_marked)$Name

udist_marked <- phyloseq::distance(kmnp_marked,"unifrac", weighted=FALSE)
udist_marked_mat <- as.matrix(udist_marked)
colnames(udist_marked_mat) <- sample_data(kmnp_marked)$Name
rownames(udist_marked_mat) <- sample_data(kmnp_marked)$Name

#######################################################################
# Bray Curtis
adonis.sex <- adonis(bdist_marked~Sex, as(sample_data(kmnp_marked),"data.frame"), permutations=999)
adonis.sex

#######################################################################
# Weighted Unifrac

adonis.sex.uw <- adonis(uwdist_marked~Sex, as(sample_data(kmnp_marked),"data.frame"), permutations=999)
adonis.sex.uw 

#######################################################################
# Unweighted Unifrac

adonis.sex.u <- adonis(udist_marked~Sex, as(sample_data(kmnp_marked),"data.frame"), permutations=999)
adonis.sex.u 

#######################################################################
## Adonis by Age
#######################################################################

#######################################################################
# Bray Curtis

adonis.age <- adonis(bdist_marked~Age, as(sample_data(kmnp_marked),"data.frame"), permutations=999)
adonis.age

#######################################################################
# Weighted Unifrac

adonis.age.uw <- adonis(uwdist_marked~Age, as(sample_data(kmnp_marked),"data.frame"), permutations=999)
adonis.age.uw 

#######################################################################
# Unweighted Unifrac

adonis.age.u <- adonis(udist_marked~Age, as(sample_data(kmnp_marked),"data.frame"), permutations=999)
adonis.age.u 

#######################################################################
## Adonis by  Dominance/Sex
#######################################################################
# Bray Curtis
# only compare adults in multi-male/multi-female groups so take out adults in single-adult groups
sex_rank <- prune_samples(sample_data(kmnp_marked)$Age=="Adult", kmnp_marked) #limit analysis to adults
sex_rank <- prune_samples(sample_data(sex_rank)$Chest != "Female", sex_rank) #limit analysis to individuals of known identity
sex_rank <- prune_samples(sample_data(sex_rank)$Name != "V-F1", sex_rank) 
sex_rank <- prune_samples(sample_data(sex_rank)$Name != "V-M2", sex_rank)
sex_rank <- prune_samples(sample_data(sex_rank)$Name != "IV-M2", sex_rank)

bdist_marked2 <- phyloseq::distance(sex_rank,"bray")
uwdist_marked2 <- phyloseq::distance(sex_rank,"wunifrac")

adonis.rank <- adonis(bdist_marked2~Chest, as(sample_data(sex_rank),"data.frame"), permutations=999)
adonis.rank

adonis.rank <- adonis(uwdist_marked2~Chest, as(sample_data(sex_rank),"data.frame"), permutations=999)
adonis.rank
