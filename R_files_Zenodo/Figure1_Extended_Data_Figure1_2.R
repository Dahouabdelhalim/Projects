###################################################################################################################
## This script contains the code to produce Figure 1 and Extended Data Figures 1 and 2
## Weigert, Hetzel et al. Dynamic antagonism between key repressive pathways maintains the placental epigenome 2023
## Author: Sara Hetzel
###################################################################################################################

library(ggplot2)
library(vioplot)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(WebGestaltR)
library(circlize)
library(data.table)
library(rtracklayer)
library(viridis)
library(gridExtra)

###################################################################################################################
## Figure 1a
###################################################################################################################

## Figure 1a was generated using the processed methylation rates and entropy tracks of TSC1 WT as well as
## Epiblast and ExE (GSE137337, replicate 1) visualized in IGV.

###################################################################################################################
## Figure 1b and Extended Data Figure 1e
###################################################################################################################

avg_features <- read.table(file = "Figure1b_Extended_Data_Figure1e_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_features_df <- melt(avg_features)
avg_features_df$feature <- factor(avg_features_df$feature, levels = c("HMD", "PMD", "exe_hyper_cgi"))
avg_features_df$variable <- factor(avg_features_df$variable, levels = rev(c("TSC4_WT", "TSC3_WT", "TSC2_WT", "TSC1_WT", "ExE_WT", "Epi_WT")))

pdf("Figure1b_Extended_Data_Figure1e.pdf", width = 7, height = 6)
par(mar=c(14, 4, 4, 2))
vioplot(value ~ variable, data = subset(avg_features_df, feature == "HMD"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), col = brewer.pal(3, "Paired")[1], main = "HMD", las = 2)
vioplot(value ~ variable, data = subset(avg_features_df, feature == "PMD"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), col = brewer.pal(3, "Paired")[2], main = "PMD", las = 2)
vioplot(value ~ variable, data = subset(avg_features_df, feature == "exe_hyper_cgi"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), col = brewer.pal(3, "Paired")[3], main = "ExE hyper CGI", las = 2)
dev.off()

###################################################################################################################
## Figure 1c and Extended Data Figure 1f
###################################################################################################################

avg_exe_hyper_cgi_entropy <- read.table(file = "Figure1c_Extended_Data_Figure1f_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

## Combine methylation and entropy
avg_exe_hyper_cgi_entropy_entropy <- avg_exe_hyper_cgi_entropy[,grepl("entropy", colnames(avg_exe_hyper_cgi_entropy))]
colnames(avg_exe_hyper_cgi_entropy_entropy) <- gsub("entropy_", "", colnames(avg_exe_hyper_cgi_entropy_entropy))
avg_exe_hyper_cgi_entropy_entropy_df <- melt(avg_exe_hyper_cgi_entropy_entropy)
colnames(avg_exe_hyper_cgi_entropy_entropy_df) <- c("variable", "entropy")

avg_exe_hyper_cgi_entropy_methylation <- avg_exe_hyper_cgi_entropy[,grepl("methylation", colnames(avg_exe_hyper_cgi_entropy))]
colnames(avg_exe_hyper_cgi_entropy_methylation) <- gsub("methylation_", "", colnames(avg_exe_hyper_cgi_entropy_methylation))

avg_exe_hyper_cgi_entropy_methylation_df <- melt(avg_exe_hyper_cgi_entropy_methylation)
colnames(avg_exe_hyper_cgi_entropy_methylation_df) <- c("variable", "methylation")

avg_exe_hyper_cgi_entropy_entropy_df$methylation <- avg_exe_hyper_cgi_entropy_methylation_df$methylation
avg_exe_hyper_cgi_entropy_entropy_df$variable <- factor(avg_exe_hyper_cgi_entropy_entropy_df$variable, levels = c("Epi_WT_Rep1", "ExE_WT_Rep1", "TSC1_WT", "TSC2_WT", "TSC3_WT", "TSC4_WT"))

pdf("Figure1c_Extended_Data_Figure1f.pdf", height = 5, width = 30)
ggplot(avg_exe_hyper_cgi_entropy_entropy_df, aes(x=methylation, y=entropy)) + geom_point(size = 0.3) + theme_classic() + xlab("Mean methylation") + ylab("Mean methylation entropy") + geom_density_2d(size = 0.5, color = "grey") + facet_grid(~variable)
dev.off()

###################################################################################################################
## Figure 1d and Extended Data Figure 2b
###################################################################################################################

entropy_methylation_clones <- read.table(file = "Figure1d_Extended_Data_Figure2b_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

annotation_samples <- data.frame(row.names = colnames(entropy_methylation_clones)[-c(1,55)], cell_type = factor(sapply(colnames(entropy_methylation_clones)[-c(1,55)], function(x) strsplit(x, "_WT")[[1]][1]), levels = c("mESC", "TSC1", "TSC2")))
annotation_samples$type <- ifelse(grepl("single_cell",colnames(entropy_methylation_clones)[-c(1,55)]), "single_cell_clone", "bulk")

entropy_methylation_clones_df <- melt(entropy_methylation_clones)
entropy_methylation_clones_df$cell_type <- factor(annotation_samples[as.character(entropy_methylation_clones_df$variable),"cell_type"], levels = c("mESC", "TSC1", "TSC2"))
entropy_methylation_clones_df$variable <- factor(entropy_methylation_clones_df$variable, levels = rownames(annotation_samples[order(annotation_samples$cell_type, annotation_samples$type),]))

colors_cell_type <- c("plum3", "orchid1", "lightsteelblue3")
names(colors_cell_type) <- c("TSC1", "TSC2", "mESC")

pdf("Figure1d_Extended_Data_Figure2b.pdf", width = 12, height = 8)
ggplot(data = subset(entropy_methylation_clones_df, measure == "entropy"), aes(x = variable, y = value, fill = cell_type)) + geom_boxplot(outlier.shape = NA) + theme_classic() + scale_fill_manual(values = colors_cell_type) + xlab("") + ylab("Entropy per 4mer in hyper CGIs") + theme(axis.text.x=element_text(size=12, hjust = 1, angle = 90), axis.text.y=element_text(size=12), axis.title=element_text(size=14)) + geom_hline(yintercept = 0.25, lty = 2, color = "darkgrey", size = 1) + ggtitle("Entropy")

ggplot(data = subset(entropy_methylation_clones_df, measure == "methylation"), aes(x = variable, y = value, fill = cell_type)) + geom_boxplot(outlier.shape = NA) + theme_classic() + scale_fill_manual(values = colors_cell_type) + xlab("") + ylab("Methylation per 4mer in hyper CGIs") + theme(axis.text.x=element_text(size=12, hjust = 1, angle = 90), axis.text.y=element_text(size=12), axis.title=element_text(size=14)) + geom_hline(yintercept = 0.25, lty = 2, color = "darkgrey", size = 1) + ggtitle("Methylation")
dev.off()

###################################################################################################################
## Figure 1e
###################################################################################################################

## The upper part of Figure 1e was generated using the processed methylation rates of TSC1 WT profiled
## with WGBS and Nanopore sequencing visualized in IGV.

###################################################################################################################
## Figure 1f
###################################################################################################################

exe_hyper_cgi_pair_ratios_with_random <- read.table(file = "Figure1f_Extended_Data_Figure2d_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

set.seed(42)
pdf("Figure1f.pdf", width = 5, height = 7)
boxplot(list(phased = exe_hyper_cgi_pair_ratios_with_random$phased, unphased = exe_hyper_cgi_pair_ratios_with_random$unphased), ylim = c(0,1), ylab = "Fraction concordant reads", outline = FALSE)

points(jitter(rep(1, nrow(exe_hyper_cgi_pair_ratios_with_random)), factor = 7), exe_hyper_cgi_pair_ratios_with_random$phased, pch = 19, cex = 0.6, col = alpha(c("HOXA" = "firebrick3", "Other" = "black")[exe_hyper_cgi_pair_ratios_with_random$pair_in_hoxa], 0.5))

points(jitter(rep(2, nrow(exe_hyper_cgi_pair_ratios_with_random)), factor = 7), exe_hyper_cgi_pair_ratios_with_random$unphased, pch = 19, cex = 0.6, col = alpha(c("HOXA" = "firebrick3", "Other" = "black")[exe_hyper_cgi_pair_ratios_with_random$pair_in_hoxa], 0.5))
dev.off()

###################################################################################################################
## Extended Data Figure 1a
###################################################################################################################

## Extended Data Figure 1a was generated using the processed methylation rates of TSC2, TSC3 and TSC4 WT
## visualized in IGV.

###################################################################################################################
## Extended Data Figure 1b
###################################################################################################################

exe_hyper_promoter_genes <- read.table(file = "Extended_Data_Figure1b_source_data.tsv", header = TRUE, stringsAsFactors = FALSE)

ora <- WebGestaltR(interestGene = exe_hyper_promoter_genes[,1], enrichMethod = "ORA", enrichDatabase = "geneontology_Biological_Process", organism = "mmusculus", referenceSet = "genome", isOutput = FALSE, interestGeneType="genesymbol", referenceGeneType="genesymbol", minNum = 10, maxNum = 500, sigMethod = "top", topThr = 20)

dotplot_enrichment <- function(data, top, selection = NULL, analysis = "ora")
{
    if (length(selection) > 0)
    {
        subset_data <- subset(data, geneSet %in% selection)
    }
    else
    {
        subset_data <- data
    }

    subset_data <- subset_data[order(subset_data$FDR),]
    subset_data <- head(subset_data, top)

    subset_data <- subset_data[order(subset_data$enrichmentRatio, decreasing = TRUE),]
    subset_data$description <- factor(subset_data$description, levels = rev(subset_data$description))
    subset_data$gene_ratio <- subset_data$overlap / subset_data$size

    p <- ggplot(subset_data, aes(x = enrichmentRatio, y = description, size = gene_ratio, color = FDR)) + geom_point() + scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.05)) + theme_light() + ylab("") + xlab("Enrichment ratio") + guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) + lims(size = c(0, 1))
    return(p)
}

pdf("Extended_Data_Figure1b.pdf", width = 8)
dotplot_enrichment(ora, top = 20)
dev.off()

###################################################################################################################
## Extended Data Figure 1c
###################################################################################################################

cgis_distance_to_tss_any_hyper <- read.table(file = "Extended_Data_Figure1c_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

pdf("Extended_Data_Figure1c.pdf", width = 10, height = 8)
ggplot(cgis_distance_to_tss_any_hyper, aes(x=distance, color = type, fill = type, group = type)) + geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]), ..count..[..group..==2]/sum(..count..[..group..==2]))), position = "identity", binwidth = 1, alpha = 0.2) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + xlab("Distance to TSS") + ylab("Fraction of CGIs") + coord_cartesian(xlim = c(-20, 20))
dev.off()

###################################################################################################################
## Extended Data Figure 2a
###################################################################################################################

avg_tile_1kb_single_sorted_clones <- read.table(file = "Extended_Data_Figure2a_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_tile_1kb_single_sorted_clones_df <- melt(avg_tile_1kb_single_sorted_clones)

pdf("Extended_Data_Figure2a.pdf", height = 7, width = 34)
par(mar=c(20, 4, 4, 2))
vioplot(value ~ variable, data = subset(avg_tile_1kb_single_sorted_clones_df, feature == "HMD"), xlab = "", ylab = "Mean methylation", main = "HMD", las = 2, col = rep(brewer.pal(3, "Paired")[1], 50), side = "left", plotCentre = "line")
vioplot(value ~ variable, data = subset(avg_tile_1kb_single_sorted_clones_df, feature == "PMD"), xlab = "", ylab = "Mean methylation", main = "PMD", col = rep(brewer.pal(3, "Paired")[2], 50), side = "right", add = TRUE, plotCentre = "line", las = 2)
dev.off()

###################################################################################################################
## Extended Data Figure 2c
###################################################################################################################

methyl_data_wgbs_nanopore <- read.table(file = "Extended_Data_Figure2c_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

pdf("Extended_Data_Figure2c.pdf", height = 5, width = 4.5)
smoothScatter(methyl_data_wgbs_nanopore[,"TSC1_WT_WGBS"], methyl_data_wgbs_nanopore[,"TSC1_WT_Nanopore"], xlab = "WGBS", ylab = "Nanopore", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0)
abline(0, 1)
abline(0.1, 1, lty = 2, col = "darkgrey")
abline(-0.1, 1, lty = 2, col = "darkgrey")
dev.off()

###################################################################################################################
## Extended Data Figure 2d
###################################################################################################################

pdf("Extended_Data_Figure2d.pdf")
plot(exe_hyper_cgi_pair_ratios_with_random$distance_pair, exe_hyper_cgi_pair_ratios_with_random$phased, pch = 19, xlab = "Distance between ExE hyper CGI pairs", ylab = "Fraction concordant reads", ylim = c(0,1), col = c("HOXA" = "firebrick3", "Other" = "black")[exe_hyper_cgi_pair_ratios_with_random$pair_in_hoxa])
dev.off()

###################################################################################################################
## Extended Data Figure 2e
###################################################################################################################

avg_features_late_placenta <- read.table(file = "Extended_Data_Figure2e_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_features_late_placenta_df <- melt(avg_features_late_placenta)
avg_features_late_placenta_df$variable <- factor(avg_features_late_placenta_df$variable, levels = c("ExE_WT", "E15_placental_junctional_zone_merged", "E15_placental_labyrinthine_zone_merged", "E18_placental_junctional_zone_merged", "E18_placental_labyrinthine_zone_merged"))

pdf("Extended_Data_Figure2e.pdf", height = 12, width = 12)
par(mar=c(18, 4, 4, 2))
layout(matrix(1:4, ncol = 2, byrow = TRUE))
vioplot(value ~ variable, data = subset(avg_features_late_placenta_df, feature == "HMD"), xlab = "", ylab = "Mean methylation", main = "HMD", las = 2, ylim = c(0,1))
vioplot(value ~ variable, data = subset(avg_features_late_placenta_df, feature == "PMD"), xlab = "", ylab = "Mean methylation", main = "PMD", las = 2, ylim = c(0,1))
vioplot(value ~ variable, data = subset(avg_features_late_placenta_df, feature == "exe_hyper_cgi"), xlab = "", ylab = "Mean methylation", main = "ExE hyper CGI", las = 2, ylim = c(0,1))
dev.off()
