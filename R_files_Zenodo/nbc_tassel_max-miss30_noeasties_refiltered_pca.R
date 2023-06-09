# Introduction ----
# nbc_tassel_max-miss30_noeasties_refiltered_pca.R
# Started by Libby Natola on 7 December 2021.
# Based on R code written by Darren Irwin for Greenish Warbler analysis, and then the NA warbler analyses, and then the "TOWA_GBS_withMaddie_R_analysis_script.R", and finally "TOWA_BTNW_HEWA_8plates_R_analysis_script.R". Made using GBS reads prepped and fileterd using the Irwin paper processing steps.

# Initial setup ----

# set directory where files stored
setwd("~/Documents/UBC/Bioinformatics/rnrb/pcas")

#I think he used this because there are so many data points and we want to see them all
options(max.print = 10000)

source("genomics_R_functions.R")
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

base.file.name <- "../012NA/nbc_max-miss30_filtered.snps_unrelated_filtered.birds_noeasties_refiltered.tab"
pos <- read.table(paste0(base.file.name, ".012.pos"), col.names = c("chrom", "position"))
column_names <- c("null", paste("c", pos$chrom, pos$position, sep="."))
geno <- read.table(paste0(base.file.name, ".012NA"), colClasses = "integer", col.names = column_names)
SNPnum <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
ind <- read.table(paste0(base.file.name, ".012.indv"))
locations <- read.table("nbc_Fst_noeasties_max-miss30_refiltered.group.txt", header=TRUE)
num_loc_cols <- ncol(locations)
ind_with_locations <- cbind(ind,locations)
combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])

# determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs
X <- 60  # this is the percentage threshold 
#start 76 individuals
threshold_NA <- SNPnum * X/100
numNAs <- rowSums(is.na(combo[(num_loc_cols+1):ncol(combo)]))
numNAs_by_ID <- data.frame(combo$ID, numNAs)  # useful to see numNAs per sample: numNAs_by_ID
selection <- (numNAs < threshold_NA)
if(any(is.na(selection))) cat("selection contains NA values\\n")  # this is a check for noticing errors / bugs
combo.NApass.all <- combo[selection,]
combo$ID[which(selection==F)]
birds.pass.spp <- locations$group[which(selection==T)]
# finish 76 individuals

###starting SNPnum:66247
# filter out SNPs with too many missing genotypes:
SNP_NAs <- colSums(is.na(combo.NApass.all[,(num_loc_cols+1):ncol(combo.NApass.all)]))
X <- 60  # this is the percentage threshold
threshold_SNP_NAs <- length(combo.NApass.all[,1]) * X/100
selection <- (SNP_NAs <= threshold_SNP_NAs)
if(any(is.na(selection))) cat("selection contains NA values\\n")  # this is a check for noticing errors / bugs
combo.NApass.subset <- combo.NApass.all[, c(rep(TRUE, times=num_loc_cols),selection)]
pos.subset <- pos[selection,]
### 66247 perfect

# option to filter out all but selected chromosome (or set of them):
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
                       groups.to.plot.PCA, group.colors.PCA, axes, flip1=T, flip2=T)

PCA_table <- data.frame(PCA_results$scores, row.names = PCA_results$data$ID)
NBC_PCA_table <- PCA_table

PCA_results$var_explained
# [1] 0.05576705 0.07868128 0.09913703
PCA_results$var_explained[2]-PCA_results$var_explained[1]
# [1] 0.02291424

### Make it purty

library(ggplot2)
library(tidyr)
library(dplyr)

NBC_PCA_table$spp <- birds.pass.spp

### save PCA_table so I don't have to rerun the whole script every time I want to tweak the figure
write.csv(NBC_PCA_table, "nbc_pca_table.csv")
#NBC_PCA_table <- read.csv("nbc_pca_table.csv")

### save it
pdf("./pca_figures/nbc_max-miss30_noeasties_refiltered_PCA.pdf")
theme_set(theme_bw())
nbc_pca<-theme_update(text = element_text(size=20), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(color=rgb(253/255, 231/255, 37/255), size=2))
ggplot(NBC_PCA_table, aes(x=PC1, y=PC2, color = spp)) + geom_point(size=2.5) + xlab("PC1 5.6%") + ylab("PC2 2.3%") + ggtitle("NBC Genomic PCA") + labs(colour = "group") + scale_color_manual(values = c("#DE73FF", "#FF4343", "#4343FF"), labels = c("Hybrid", "Red-breasted", "Red-naped")) + coord_fixed(0.5)
dev.off()

####################################################################################################


# ##### redo it all with the admixture assignments
# 
# groups_and_colors <- rbind(c("RNSA", "blue"),
#                            c("RBSA", "red"),
#                            c("OTSA", "black"),
#                            c("YBSA", "gold")
# )
# 
# groups.to.plot.PCA <- groups_and_colors[,1]
# group.colors.PCA <- groups_and_colors[,2]
# groups <- c("RNSA", "RBSA", "OTSA", "YBSA")  # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
# 
# base.file.name <- "012NA/pg_variants_whole_genome_MQ20_Het06_NoOut_70Miss_NoDups_biallelic_GQ10_noIndel_maf005_MaxMiss70.tab"
# pos <- read.table(paste0(base.file.name, ".012.pos"), col.names = c("chrom", "position"))
# column_names <- c("null", paste("c", pos$chrom, pos$position, sep="."))
# geno <- read.table(paste0(base.file.name, ".012NA"), colClasses = "integer", col.names = column_names)
# SNPnum <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
# ind <- read.table(paste0(base.file.name, ".012.indv"))
# locations <- read.table("pg_genotype_groups_otherz.Fst_groups.txt", header=TRUE)
# num_loc_cols <- ncol(locations)
# ind_with_locations <- cbind(ind,locations)
# combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
# 
# # determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs
# X <- 60  # this is the percentage threshold  (# tried this at 60 but seemed to be sizeable plate effects, with Haley's plates missing more data and being in strange position in PCA)
# #start 104 individuals
# threshold_NA <- SNPnum * X/100
# numNAs <- rowSums(is.na(combo[(num_loc_cols+1):ncol(combo)]))
# numNAs_by_ID <- data.frame(combo$ID, numNAs)  # useful to see numNAs per sample: numNAs_by_ID
# selection <- (numNAs < threshold_NA)
# if(any(is.na(selection))) cat("selection contains NA values\\n")  # this is a check for noticing errors / bugs
# combo.NApass.all <- combo[selection,]
# combo$ID[which(selection==F)]
# 
# # [1] "PG_2019_TE07L01" "PG_2019_TE25L02" "PG_2019_TE25L03" "PG_2019_TE27L01" "PG_2019_TE28L03"
# # [6] "PG_2019_TE28L04" "PG_2019_TF13L01" "PG_2019_TF13L02" "PG_2019_TF14L02" "PG_2019_TF15L04"
# # [11] "PG_2019_TF23L01" "PG_2019_TF23L03" "PG_2020_UE15L01" "PG_2020_UE20L01"
# # # finish 91 individuals
# 
# ###starting SNPnum:122026
# # filter out SNPs with too many missing genotypes:
# SNP_NAs <- colSums(is.na(combo.NApass.all[,(num_loc_cols+1):ncol(combo.NApass.all)]))
# X <- 60  # this is the percentage threshold
# threshold_SNP_NAs <- length(combo.NApass.all[,1]) * X/100
# selection <- (SNP_NAs <= threshold_SNP_NAs)
# if(any(is.na(selection))) cat("selection contains NA values\\n")  # this is a check for noticing errors / bugs
# combo.NApass.subset <- combo.NApass.all[, c(rep(TRUE, times=num_loc_cols),selection)]
# pos.subset <- pos[selection,]
# ### 122020
# 
# # option to filter out all but selected chromosome (or set of them):
# choose.chrom <- F
# if (choose.chrom == TRUE) {
#   chrom <- "Z"
#   # selection <- (pos.subset$chrom == chrom)
#   selection <- (pos.subset$chrom == chrom)
#   if(any(is.na(selection))) cat("selection contains NA values\\n")  # this is a check for noticing errors / bugs
#   #pos.subset.one.chr <- pos.subset[selection,]
#   #loci.selection <- c(rep(TRUE, times=num_loc_cols), selection)  # add placeholders for info columns
#   # which(loci.selection == T)    # to check which entries are TRUE
#   combo.NApass <- combo.NApass.subset[, c(rep(TRUE, times=num_loc_cols), selection)]
#   pos.NApass <- pos.subset[selection,]
#   region.text <- paste0("chr", chrom)	
# }	else {
#   region.text <- "whole_genome"
#   combo.NApass <- combo.NApass.subset
#   pos.NApass <- pos.subset
# }
# 
# # Calculate allele freqs and sample sizes (use column Fst_group)
# temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
# freqs <- temp.list$freqs
# sample_size <- temp.list$sample_size
# rm(temp.list)
# # calculate WC84_Fst 
# temp.list <- getWC84Fst(freqs, sample_size, groups, among=TRUE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
# WC84_Fst <- temp.list$WC84_Fst
# rm(temp.list)
# 
# # jmake histograms of Fst
# rnrb <- hist(WC84_Fst[1,])
# rnot <- hist(WC84_Fst[2,])
# rnyb <- hist(WC84_Fst[3,])
# rbot <- hist(WC84_Fst[4,])
# rbyb <- hist(WC84_Fst[5,])
# otyb <- hist(WC84_Fst[6,])
# among <- hist(WC84_Fst[7,])
# 
# ### see heterozygosity
# ggplot(combo.NApass, )
# 
# library(tidyr)
# library(dplyr)
# boxplot(WC84_Fst[1:7,])
# 
# ### transpose data columns and rows
# WC84_Fst_better <- t(WC84_Fst)
# 
# ### make it a data frame
# WC84_Fst_better <- as.data.frame(WC84_Fst_better)
# 
# ### plot histograms with ggplot
# ggplot(WC84_Fst_better, aes(x=RNSA_RBSA)) + geom_histogram()
# 
# ### try violin plot
# ggplot(WC84_Fst_better, aes(x='count', y=RNSA_YBSA)) + geom_violin()
# 
# ### rats I deleted the instructions to get the locus column, something like this...
# # WC84_Fst_betterer %>%
# #   rownames_to_column('locus') %>%
# #   gather(comparison, Fst, RNSA_RBSA:Fst_among) %>%
# #   column_to_rownames('locus')
# 
# ### gather data so you can make one plot with all the comparisons violin plots
# WC84_Fst_betterer <- gather(WC84_Fst_better, comparison, Fst, RNSA_RBSA:Fst_among)
# 
# ### plot as violin plot
# ggplot(WC84_Fst_betterer, aes(x=comparison, y=Fst, color=comparison)) + geom_violin() + scale_color_manual(values=c("black", "goldenrod3", "dark red", "orange", "dark blue", "purple", "green"))
# 
# ### plot boxplot
# ggplot(WC84_Fst_betterer, aes(x=comparison, y=Fst, color=comparison)) + geom_boxplot() + scale_color_manual(values=c("black", "goldenrod3", "dark red", "orange", "dark blue", "purple", "green"))
# 
# ### remove the weird categories
# WC84_Fst_bettererer <- WC84_Fst_betterer %>% group_by(comparison) %>% filter(comparison != "HYSA_YBSA") %>% filter(comparison != "RBSA_HYSA") %>% filter(comparison != "Fst_among") %>% filter(comparison != "RNSA_HYSA")
# 
# ### plot
# ggplot(WC84_Fst_bettererer, aes(x=comparison, y=Fst, color=comparison)) + geom_violin() + scale_color_manual(values=c("orange", "purple", "green")) + stat_summary(fun.y=mean, geom="point", size=0.5, color="black")
# ggplot(WC84_Fst_bettererer, aes(x=comparison, y=Fst, color=comparison)) + geom_boxplot() + scale_color_manual(values=c("orange", "purple", "green"))
# 
# ### out of curiosity, how many fsts do we have that are 1?
# Fst_1 <- WC84_Fst_bettererer %>% filter(Fst == 1)
# Fst_0.95 <- WC84_Fst_bettererer %>% filter(Fst > 0.95)
# Fst_0ish <- WC84_Fst_bettererer %>% filter(Fst < 0.05)
# 
# ### plot Fst of 1 by spp comparison
# ggplot(Fst_0.95, aes(x=comparison, fill=comparison)) + geom_bar() +scale_fill_manual(values = c("orange","green"))
# ggplot(Fst_0ish, aes(x=comparison, fill=comparison)) + geom_bar() + scale_fill_manual(values = c("orange", "purple", "green"))
# 
# ### I wanna gander at these regions
# write.csv(Fst_1,"Fst_1.csv")
# 
# # make the figure:
# Fst.filter <- F   # option to filter to high-Fst markers only, using cutoff below
# Fst.cutoff <- 0.25  # has no effect if Fst.filter is FALSE
# # choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
# groups.to.compare <- "Fst_among"
# axes <- 3
# combo.NApass_noTF09L03 <- combo.NApass[-c(48), ] 
# 
# 
# PCA_results <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass_noTF09L03, num_loc_cols, region.text,
#                        groups.to.plot.PCA, group.colors.PCA, axes, flip1=F, flip2=T)
# 
# PCA_table_noTF09L03 <- data.frame(PCA_results$scores, row.names = PCA_results$data$ID)
# 
# PCA_results$var_explained
# #[1] 0.09839706 0.12255781 0.14347847
# 0.12255781 -0.09839706
# # [1] 0.0241607
# 
# ### Make it purty
# 
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# 
# pca_table_plot <- read.csv("../pca/pg_pca_table.csv")
# 
# a <- ggplot(pca_table_plot, aes(x=PC1, y=PC2)) + geom_point()
# a
# 
# theme_set(theme_bw())
# theme_update(plot.title = element_text(hjust = 0.5))
# ggplot(pca_table_plot, aes(x=PC1, y=PC2, color = group)) + geom_point() + xlab("PC1 9.8%") + ylab("PC2 2.4%") + ggtitle("Genomic PCA") + labs(colour = "group") + scale_color_manual(values=c("dark grey", "red", "blue", "yellow"))
# 
# ### PG_2019_TF09L01, PG_2019_TF09L02, PG_2019_TF09L03 are all pulling out on PC3, they're 2 parents and an offspring? Let's just pull PG_2019_TF09L03 (row 48) out and see what it looks like without it
# 
# no_fam <- read.csv("../pca/pg_pca_table_noTF09L03.csv")
# ggsave(file = "../pca/pca_professional.pdf", dpi =500, dev=pdf)
# ggplot(no_fam, aes(x=PC1, y=PC2, color = group)) + geom_point() + xlab("PC1 9.8%") + ylab("PC2 2.4%") + ggtitle("Genomic PCA") + labs(colour = "Species") + scale_color_manual(values=c("dark grey", "red", "blue", "yellow"))
# dev.off()
# 
