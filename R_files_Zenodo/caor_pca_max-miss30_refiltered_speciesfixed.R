# Introduction ----
# caor_tassel_max-miss30_refiltered_speciesfixed_pca.R
# Started by Libby Natola on 7 December 2021.
# Based on R code written by Darren Irwin for Greenish Warbler analysis, and then the NA warbler analyses, and then the "TOWA_GBS_withMaddie_R_analysis_script.R", and finally "TOWA_BTNW_HEWA_8plates_R_analysis_script.R". Made using GBS reads prepped and fileterd using the Irwin paper processing steps.

# Initial setup ----

# set directory where files stored
setwd("~/Documents/UBC/Bioinformatics/rnrb/pcas")

#I think he used this because there are so many data points and we want to see them all
options(max.print = 10000)

source("./genomics_R_functions.R")
# PCA whole-genome ----
# Load vcf file containing only variable sites throughout genome;
# construct a PCA based on all sites passing an Fst threshold between the "groups" below;
# and all individuals in "groups.to.plot.PCA" according to colors in "group.colors.PCA"
# groups <- c("RNSA", "RBSA", "HYSA")  # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
# groups.to.plot.PCA <- c("RNSA", "RBSA", "HYSA") 
# group.colors.PCA <- c("blue","red","grey")

groups_and_colors <- rbind(c("RNSA", "blue"),
                           c("RBSA", "red"),
                           c("HYSA", "grey")
)

groups.to.plot.PCA <- groups_and_colors[,1]
group.colors.PCA <- groups_and_colors[,2]
groups <- c("RNSA", "RBSA", "HYSA")  # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)

base.file.name <- "../012NA/caor_max-miss30_filtered.snps_unrelated_filtered.birds.tab"
pos <- read.table(paste0(base.file.name, ".012.pos"), col.names = c("chrom", "position"))
column_names <- c("null", paste("c", pos$chrom, pos$position, sep="."))
geno <- read.table(paste0(base.file.name, ".012NA"), colClasses = "integer", col.names = column_names)
SNPnum <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
ind <- read.table(paste0(base.file.name, ".012.indv"))
locations <- read.table("caor_Fst_max-miss30_refiltered_fixedspp.group.txt", header=TRUE)
num_loc_cols <- ncol(locations)
ind_with_locations <- cbind(ind,locations)
combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])

# determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs, but I've already done this to 70 bash and vcf tools
X <- 75  # this is the percentage threshold 
#start 108 individuals
threshold_NA <- SNPnum * X/100
numNAs <- rowSums(is.na(combo[(num_loc_cols+1):ncol(combo)]))
numNAs_by_ID <- data.frame(combo$ID, numNAs)  # useful to see numNAs per sample: numNAs_by_ID
selection <- (numNAs < threshold_NA)
if(any(is.na(selection))) cat("selection contains NA values\\n")  # this is a check for noticing errors / bugs
combo.NApass.all <- combo[selection,]
combo$ID[which(selection==F)]
birds.pass.spp <- locations$group[which(selection==T)]

# finish 108 individuals, good didn't lose any

###starting SNPnum:29236
# filter out SNPs with too many missing genotypes, have also already done this to max 20% missing:
SNP_NAs <- colSums(is.na(combo.NApass.all[,(num_loc_cols+1):ncol(combo.NApass.all)]))
X <- 75  # this is the percentage threshold
threshold_SNP_NAs <- length(combo.NApass.all[,1]) * X/100
selection <- (SNP_NAs <= threshold_SNP_NAs)
if(any(is.na(selection))) cat("selection contains NA values\\n")  # this is a check for noticing errors / bugs
combo.NApass.subset <- combo.NApass.all[, c(rep(TRUE, times=num_loc_cols),selection)]
pos.subset <- pos[selection,]
### 29236, cool ...

# option to filter out all but selected chromosome (or set of them):
# can't actually use this everything is on chromosome 1 bc it is de novo
choose.chrom <- F
if (choose.chrom == TRUE) {
  chrom <- "PseudoWWNC01000009.1_Melanerpes_aurifrons_Melanerpes_aurifrons_OMNH24340_contig_9,_whole_genome_shotgun_sequence"
  # selection <- (pos.subset$chrom == chrom)
  selection <- (pos.subset$chrom == chrom)
  if(any(is.na(selection))) cat("selection contains NA values\\n")  # this is a check for noticing errors / bugs
  #pos.subset.one.chr <- pos.subset[selection,]
  #loci.selection <- c(rep(TRUE, times=num_loc_cols), selection)  # add placeholders for info columns
  # which(loci.selection == T)    # to check which entries are TRUE
  combo.NApass <- combo.NApass.subset[, c(rep(TRUE, times=num_loc_cols), selection)]
  pos.NApass <- pos.subset[selection,]
  region.text <- paste0("chr", chrom)	
}	else {
  region.text <- "whole_genome"
  combo.NApass <- combo.NApass.subset
  pos.NApass <- pos.subset
}

# Calculate allele freqs and sample sizes (use column Fst_group)
temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
freqs <- temp.list$freqs
sample_size <- temp.list$sample_size
rm(temp.list)
# calculate WC84_Fst 
temp.list <- getWC84Fst(freqs, sample_size, groups, among=TRUE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
WC84_Fst <- temp.list$WC84_Fst
rm(temp.list)

# make the figure:
Fst.filter <- F   # option to filter to high-Fst markers only, using cutoff below
Fst.cutoff <- 0.25  # has no effect if Fst.filter is FALSE
# choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Fst_among"
axes <- 3


PCA_results <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
                       groups.to.plot.PCA, group.colors.PCA, axes, flip1=F, flip2=T)

PCA_table <- data.frame(PCA_results$scores, row.names = PCA_results$data$ID)
CAOR_PCA_table <- PCA_table

PCA_results$var_explained
# [1] 0.04232219 0.06539801 0.08514751
PCA_results$var_explained[2]-PCA_results$var_explained[1]
# [1] 0.0204911

### Make it purty

library(ggplot2)
library(tidyr)
library(dplyr)

CAOR_PCA_table$spp <- birds.pass.spp

### save PCA_table so I don't have to rerun the whole script every time I want to tweak the figure
write.csv(CAOR_PCA_table, "caor_pca_table.csv")
# CAOR_PCA_table <- read.csv("caor_pca_table.csv")

### save it
pdf("./pca_figures/caor_max-miss30_refiltered_fixedspp_PCA.pdf")
theme_set(theme_bw())
caor_pca<-theme_update(text = element_text(size=20), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(color=rgb(33/255, 145/255, 140/255), size = 2))
ggplot(CAOR_PCA_table, aes(x=PC1, y=PC2, color = spp)) + geom_point(size=2.5) + xlab("PC1 4.2%") + ylab("PC2 2.0%") + ggtitle("CA/OR Genomic PCA") + labs(colour = "group") + scale_color_manual(values = c("#DE73FF", "#FF4343", "#4343FF"), labels = c("Hybrid", "Red-breasted", "Red-naped")) + coord_fixed(0.5)
dev.off()

####################################################################################################

