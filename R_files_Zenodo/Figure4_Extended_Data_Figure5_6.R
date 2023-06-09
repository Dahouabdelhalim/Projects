###################################################################################################################
## This script contains the code to produce Figure 4 and Extended Data Figures 5 and 6
## Weigert, Hetzel et al. Dynamic antagonism between key repressive pathways maintains the placental epigenome 2023
## Author: Sara Hetzel
###################################################################################################################

library(ggplot2)
library(vioplot)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)
library(data.table)
library(ggridges)
library(ggrepel)
library(rtracklayer)
library(Vennerable)
library(gridExtra)
library(viridis)
library(xlsx)
library(ggrastr)
theme_set(theme_ridges())

###################################################################################################################
## Figure 4b
###################################################################################################################

avg_features <- read.table(file = "Figure4b_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_features_df <- melt(avg_features)
avg_features_df$variable <- factor(avg_features_df$variable, levels = c("TSC1_WT", "TSC2_3BKO", "TSC1_TET3KO", "TSC1_KDM2BKO", "TSC1_RNF2KO", "TSC3_EEDKO"))

pdf("Figure4b.pdf.pdf", width = 10, height = 7)
par(mar=c(12, 4, 4, 2))
vioplot(value ~ variable, data = subset(avg_features_df, feature == "HMD"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), main = "HMD/PMD", las = 2, plotCentre = "line", side = "left", col = brewer.pal(4, "Paired")[1])
vioplot(value ~ variable, data = subset(avg_features_df, feature == "PMD"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), add = TRUE, plotCentre = "line", side = "right", col = brewer.pal(4, "Paired")[2])

vioplot(value ~ variable, data = subset(avg_features_df, feature == "CGI"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), main = "CGI/ExE hyper CGI", las = 2, plotCentre = "line", side = "left", col = brewer.pal(4, "Paired")[4])
vioplot(value ~ variable, data = subset(avg_features_df, feature == "exe_hyper_cgi"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), add = TRUE, plotCentre = "line", side = "right", col = brewer.pal(4, "Paired")[3])
dev.off()

###################################################################################################################
## Figure 4c
###################################################################################################################

## Figure 4c was generated using the processed methylation rates and smoothed MINUTE-ChIP tracks of TSC1 WT
## as well as the processed methylation rates of TSC2 DNMT3B KO, TSC1 TET3 KO, TSC1 KDM2B KO, TSC1 RNF2 KO
## and TSC3 EED KO visualized in IGV.

###################################################################################################################
## Figure 4d and Extended Data Figure 5c,e
###################################################################################################################

methyl_data <- data.frame(fread("Figure4d_Extended_Data_Figure5c_5e_source_data.tsv"), stringsAsFactors = FALSE)

pdf("Figure4d_Extended_Data_Figure5c_5e_smooth_scatter.pdf", width = 24, height = 5)
layout(matrix(1:5, ncol = 5, byrow = TRUE))
smoothScatter(methyl_data[,"TSC1_WT"], methyl_data[,"TSC2_3BKO"], xlab = "TSC WT", ylab = "TSC DNMT3B KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
smoothScatter(methyl_data[,"TSC1_WT"], methyl_data[,"TSC1_TET3KO"], xlab = "TSC WT", ylab = "TSC TET3 KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
smoothScatter(methyl_data[,"TSC1_WT"], methyl_data[,"TSC1_KDM2BKO"], xlab = "TSC WT", ylab = "TSC KDM2B KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
smoothScatter(methyl_data[,"TSC1_WT"], methyl_data[,"TSC1_RNF2KO"], xlab = "TSC WT", ylab = "TSC RNF2 KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
smoothScatter(methyl_data[,"TSC1_WT"], methyl_data[,"TSC3_EEDKO"], xlab = "TSC WT", ylab = "TSC EED KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")

smoothScatter(methyl_data[,"TSC3_EEDKO"], methyl_data[,"TSC1_RNF2KO"], xlab = "TSC EED KO", ylab = "TSC RNF2 KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
smoothScatter(methyl_data[,"TSC3_EEDKO"], methyl_data[,"TSC1_KDM2BKO"], xlab = "TSC EED KO", ylab = "TSC KDM2B KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
smoothScatter(methyl_data[,"TSC1_RNF2KO"], methyl_data[,"TSC1_KDM2BKO"], xlab = "TSC RNF2 KO", ylab = "TSC KDM2B KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
dev.off()

## Barplots
diff_per_cpg <- data.frame(TSC2_3BKO = methyl_data[,"TSC2_3BKO"] - methyl_data[,"TSC1_WT"],
                           TSC1_TET3KO = methyl_data[,"TSC1_TET3KO"] - methyl_data[,"TSC1_WT"],
                           TSC1_KDM2BKO = methyl_data[,"TSC1_KDM2BKO"] - methyl_data[,"TSC1_WT"],
                           TSC1_RNF2KO = methyl_data[,"TSC1_RNF2KO"] - methyl_data[,"TSC1_WT"],
                           TSC3_EEDKO = methyl_data[,"TSC3_EEDKO"] - methyl_data[,"TSC1_WT"])

count_per_cpg_class <- matrix(NA, nrow = 3, ncol = 5)
colnames(count_per_cpg_class) <- c("TSC2_3BKO", "TSC1_TET3KO", "TSC1_KDM2BKO", "TSC1_RNF2KO", "TSC3_EEDKO")
rownames(count_per_cpg_class) <- c("hyper", "stable", "hypo")
count_per_cpg_class["hyper",] <- apply(diff_per_cpg, 2, function(x) length(which(x > 0.2)))
count_per_cpg_class["hypo",] <- apply(diff_per_cpg, 2, function(x) length(which(x < -0.2)))
count_per_cpg_class["stable",] <- apply(diff_per_cpg, 2, function(x) length(which(x >= -0.2 & x <= 0.2)))
count_per_cpg_class <- data.frame(count_per_cpg_class)
count_per_cpg_class$feature <- rownames(count_per_cpg_class)
count_per_cpg_class_df <- melt(count_per_cpg_class)
count_per_cpg_class_df$feature <- factor(count_per_cpg_class_df$feature, levels = c("hypo", "stable", "hyper"))

pdf("Figure4d_Extended_Data_Figure5c_barplot.pdf")
ggplot(data = count_per_cpg_class_df, aes(x = variable, fill = feature, y = value)) + geom_bar(position = "fill", stat = "identity") + theme_classic() + xlab("") + ylab("Fraction CpGs") + theme(axis.text.x=element_text(size=10, hjust = 1, angle = 45), axis.text.y=element_text(size=10), axis.title=element_text(size=12)) + scale_fill_manual(values = c("darkblue", "grey", "darkred"))
dev.off()

###################################################################################################################
## Figure 4e and Extended Data Figure 5g
###################################################################################################################

avg_tile_1kb_chip_wgbs <- read.table(file = "Figure4e_Extended_Data_Figure5g_5h_6a_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

pdf("Figure4e_Extended_Data_Figure5g.pdf", width = 20, height = 5)
layout(matrix(1:4, ncol = 4))
smoothScatter(avg_tile_1kb_chip_wgbs$diff_3BKO, avg_tile_1kb_chip_wgbs$log2_ratio_3BKO_H3K27me3, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H3K27me3 (KO vs WT)", nrpoints = 0, main = "3BKO - H3K27me3", ylim = c(-4.5,2.5), xlim = c(-1, 0.5))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_3BKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_3BKO_H3K27me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_KDM2BKO, avg_tile_1kb_chip_wgbs$log2_ratio_KDM2BKO_H3K27me3, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H3K27me3 (KO vs WT)", nrpoints = 0, main = "KDM2BKO - H3K27me3", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_KDM2BKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_KDM2BKO_H3K27me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_RNF2KO, avg_tile_1kb_chip_wgbs$log2_ratio_RNF2KO_H3K27me3, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H3K27me3 (KO vs WT)", nrpoints = 0, main = "RNF2KO - H3K27me3", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_RNF2KO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_RNF2KO_H3K27me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_EEDKO, avg_tile_1kb_chip_wgbs$log2_ratio_EEDKO_H3K27me3, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H3K27me3 (KO vs WT)", nrpoints = 0, main = "EEDKO - H3K27me3", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_EEDKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_EEDKO_H3K27me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)

smoothScatter(avg_tile_1kb_chip_wgbs$diff_3BKO, avg_tile_1kb_chip_wgbs$log2_ratio_3BKO_H2AK119ub, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H2AK119ub (KO vs WT)", nrpoints = 0, main = "3BKO - H2AK119ub", ylim = c(-4.5,2.5), xlim = c(-1, 0.5))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_3BKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_3BKO_H2AK119ub"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_KDM2BKO, avg_tile_1kb_chip_wgbs$log2_ratio_KDM2BKO_H2AK119ub, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H2AK119ub (KO vs WT)", nrpoints = 0, main = "KDM2BKO - H2AK119ub", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_KDM2BKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_KDM2BKO_H2AK119ub"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_RNF2KO, avg_tile_1kb_chip_wgbs$log2_ratio_RNF2KO_H2AK119ub, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H2AK119ub (KO vs WT)", nrpoints = 0, main = "RNF2KO - H2AK119ub", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_RNF2KO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_RNF2KO_H2AK119ub"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_EEDKO, avg_tile_1kb_chip_wgbs$log2_ratio_EEDKO_H2AK119ub, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H2AK119ub (KO vs WT)", nrpoints = 0, main = "EEDKO - H2AK119ub", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_EEDKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_EEDKO_H2AK119ub"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)

smoothScatter(avg_tile_1kb_chip_wgbs$diff_3BKO, avg_tile_1kb_chip_wgbs$log2_ratio_3BKO_H3K4me3, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H3K4me3 (KO vs WT)", nrpoints = 0, main = "3BKO - H3K4me3", ylim = c(-4.5,2.5), xlim = c(-1, 0.5))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_3BKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_3BKO_H3K4me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_KDM2BKO, avg_tile_1kb_chip_wgbs$log2_ratio_KDM2BKO_H3K4me3, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H3K4me3 (KO vs WT)", nrpoints = 0, main = "KDM2BKO - H3K4me3", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_KDM2BKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_KDM2BKO_H3K4me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_RNF2KO, avg_tile_1kb_chip_wgbs$log2_ratio_RNF2KO_H3K4me3, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H3K4me3 (KO vs WT)", nrpoints = 0, main = "RNF2KO - H3K4me3", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_RNF2KO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_RNF2KO_H3K4me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_EEDKO, avg_tile_1kb_chip_wgbs$log2_ratio_EEDKO_H3K4me3, xlab = "Delta DNAme (KO vs WT)", ylab = "Log2 ratio H3K4me3 (KO vs WT)", nrpoints = 0, main = "EEDKO - H3K4me3", ylim = c(-4.5,2.5), xlim = c(-0.5, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_EEDKO"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_EEDKO_H3K4me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
dev.off()

###################################################################################################################
## Figure 4f and Extended Data Figure 6d
###################################################################################################################

avg_cgi_chip_wgbs <- read.table(file = "Figure4f_Extended_Data_Figure6d_6e_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

pdf("Figure4f_Extended_Data_Figure6d.pdf", width = 11, height = 12)
p1 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = diff_EEDKO, y = diff_KDM2BKO, color = log2_cut_TSC1_KDM2BKO_H3K4me3)) + rasterise(geom_point(size = 0.5), dpi = 500) + theme_classic() + scale_color_gradientn(colors = viridis(10)) + coord_cartesian(xlim = c(-0.5, 1), ylim = c(-0.5, 1)) + xlab("Delta EEDKO") + ylab("Delta KDM2BKO") + theme(legend.position="bottom") + geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2)
p2 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = diff_EEDKO, y = diff_KDM2BKO, color = log2_cut_TSC3_EEDKO_H3K4me3)) + rasterise(geom_point(size = 0.5), dpi = 500) + theme_classic() + scale_color_gradientn(colors = viridis(10)) + coord_cartesian(xlim = c(-0.5, 1), ylim = c(-0.5, 1)) + xlab("Delta EEDKO") + ylab("Delta KDM2BKO") + theme(legend.position="bottom") + geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2)
p3 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = diff_EEDKO, y = diff_RNF2KO, color = log2_cut_TSC1_RNF2KO_H3K4me3)) + rasterise(geom_point(size = 0.5), dpi = 500) + theme_classic() + scale_color_gradientn(colors = viridis(10)) + coord_cartesian(xlim = c(-0.5, 1), ylim = c(-0.5, 1)) + xlab("Delta EEDKO") + ylab("Delta RNF2KO") + theme(legend.position="bottom") + theme(legend.position="bottom") + geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2)
p4 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = diff_EEDKO, y = diff_RNF2KO, color = log2_cut_TSC3_EEDKO_H3K4me3)) + rasterise(geom_point(size = 0.5), dpi = 500) + theme_classic() + scale_color_gradientn(colors = viridis(10)) + coord_cartesian(xlim = c(-0.5, 1), ylim = c(-0.5, 1)) + xlab("Delta EEDKO") + ylab("Delta RNF2KO") + theme(legend.position="bottom") + theme(legend.position="bottom") + geom_hline(yintercept = 0, lty = 2) + geom_vline(xintercept = 0, lty = 2)
grid.arrange(p1,p2,p3,p4, ncol = 2)
dev.off()

###################################################################################################################
## Figure 4g and Extended Data Figure 6c
###################################################################################################################

fig4g_extended_data_fig6c_source_data <- read.table(file = "Figure4g_Extended_Data_Figure6c_source_data.tsv", sep = "\\t", header = TRUE, stringsAsFactors = FALSE)

samples <- unique(fig4g_extended_data_fig6c_source_data$sample)

norm_mat_list <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(fig4g_extended_data_fig6c_source_data, sample == x)[,!colnames(fig4g_extended_data_fig6c_source_data) %in% c("sample", "gene_id")]), k_upstream=100, k_downstream=100, k_target=67, extend = c(5000, 5000), signal_name = x, target_name = "prc_hyper_cgi"))
names(norm_mat_list) <- samples

samples_for_heatmap <- c("TSC1_WT", "TSC2_3BKO", "TSC1_WT_pool2", "TSC1_KDM2BKO", "TSC1_RNF2KO", "TSC3_EEDKO")
samples_pool1 <- c("TSC1_WT_H2AK119ub", "TSC1_WT_H3K27me3", "TSC1_WT_H3K4me3", "TSC2_3BKO_H2AK119ub", "TSC2_3BKO_H3K27me3", "TSC2_3BKO_H3K4me3")

ht_list1 <- NULL
for(sample in paste(samples_for_heatmap, "H3K4me3", sep = "_"))
{
    if (sample %in% samples_pool1)
    {
        ht_list1 <- ht_list1 + EnrichedHeatmap(norm_mat_list[[sample]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 200), c("white", "seagreen3")), column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 50))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))

    }
    else
    {
        ht_list1 <- ht_list1 + EnrichedHeatmap(norm_mat_list[[sample]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 100), c("white", "seagreen3")), column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 50))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
    }
}
ht_list1 <- ht_list1 + EnrichedHeatmap(norm_mat_list[["CpG_density"]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 0.2), c("white", "black")), column_title = "CpG_density", name = "CpG_density", show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 0.2))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))


ht_list2 <- NULL
for(sample in paste(samples_for_heatmap, "H3K27me3", sep = "_"))
{
    if (sample %in% samples_pool1)
    {
        ht_list2 <- ht_list2 + EnrichedHeatmap(norm_mat_list[[sample]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 20), c("white", "violetred2")), column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 10))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))

    }
    else
    {
        ht_list2 <- ht_list2 + EnrichedHeatmap(norm_mat_list[[sample]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 10), c("white", "violetred2")), column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 5))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
    }
}
ht_list2 <- ht_list2 + EnrichedHeatmap(norm_mat_list[["CpG_density"]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 0.2), c("white", "black")), column_title = "CpG_density", name = "CpG_density", show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 0.2))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))


ht_list3 <- NULL
for(sample in paste(samples_for_heatmap, "H2AK119ub", sep = "_"))
{
    if (sample %in% samples_pool1)
    {
        ht_list3 <- ht_list3 + EnrichedHeatmap(norm_mat_list[[sample]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 20), c("white", "deepskyblue3")), column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 10))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
    }
    else
    {
        ht_list3 <- ht_list3 + EnrichedHeatmap(norm_mat_list[[sample]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 10), c("white", "deepskyblue3")), column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 5))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
    }
}
ht_list3 <- ht_list3 + EnrichedHeatmap(norm_mat_list[["CpG_density"]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 0.2), c("white", "black")), column_title = "CpG_density", name = "CpG_density", show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 0.2))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))

samples_for_heatmap_wbgs_double_wt <- c("WGBS_TSC1_WT", "WGBS_TSC2_3BKO", "WGBS_TSC1_WT_pool2", "WGBS_TSC1_KDM2BKO", "WGBS_TSC1_RNF2KO", "WGBS_TSC3_EEDKO")

actual_samples_for_heatmap_wbgs_double_wt <- c("WGBS_TSC1_WT", "WGBS_TSC2_3BKO", "WGBS_TSC1_WT", "WGBS_TSC1_KDM2BKO", "WGBS_TSC1_RNF2KO", "WGBS_TSC3_EEDKO")
names(actual_samples_for_heatmap_wbgs_double_wt) <- samples_for_heatmap_wbgs_double_wt

ht_list4 <- NULL
for(sample in samples_for_heatmap_wbgs_double_wt)
{
    ht_list4 <- ht_list4 + EnrichedHeatmap(norm_mat_list[[actual_samples_for_heatmap_wbgs_double_wt[sample]]], row_order = 1:nrow(norm_mat_list[[actual_samples_for_heatmap_wbgs_double_wt[sample]]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 1))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
}
ht_list4 <- ht_list4 + EnrichedHeatmap(norm_mat_list[["CpG_density"]], row_order = 1:nrow(norm_mat_list[[sample]]), use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = colorRamp2(c(0, 0.2), c("white", "black")), column_title = "CpG_density", name = "CpG_density", show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, 0.2))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))

pdf("Figure4g_Extended_Data_Figure6c.pdf", width = 20)
draw(ht_list1)
draw(ht_list2)
draw(ht_list3)
draw(ht_list4)
dev.off()

###################################################################################################################
## Extended Data Figure 5a
###################################################################################################################

## Extended Data Figure 5a was generated using the WGBS BAM files of TSC1 WT, TSC2 DNMT3B KO, TSC1 TET3 KO,
## TSC1 KDM2B KO, TSC1 RNF2 KO and TSC3 EED KO visualized in IGV.

###################################################################################################################
## Extended Data Figure 5b
###################################################################################################################

## With same WT
diff_features <- data.frame(
    row.names = rownames(avg_features),
    TSC2_3BKO_vs_TSC1_WT = avg_features$TSC2_3BKO - avg_features$TSC1_WT,
    TSC1_TET3KO_vs_TSC1_WT = avg_features$TSC1_TET3KO - avg_features$TSC1_WT,
    TSC1_KDM2BKO_vs_TSC1_WT = avg_features$TSC1_KDM2BKO - avg_features$TSC1_WT,
    TSC1_RNF2KO_vs_TSC1_WT = avg_features$TSC1_RNF2KO - avg_features$TSC1_WT,
    TSC3_EEDKO_vs_TSC1_WT = avg_features$TSC3_EEDKO - avg_features$TSC1_WT,
    feature = avg_features$feature
)

diff_feature_df <- melt(diff_features)
diff_feature_df$variable <- factor(diff_feature_df$variable, levels = rev(c("TSC2_3BKO_vs_TSC1_WT", "TSC1_TET3KO_vs_TSC1_WT", "TSC1_KDM2BKO_vs_TSC1_WT", "TSC1_RNF2KO_vs_TSC1_WT", "TSC3_EEDKO_vs_TSC1_WT")))

## With matching WT
avg_features_matching_wt <- read.table(file = "Extended_Data_Figure5b_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

diff_features_matching_wt <- data.frame(
    row.names = rownames(avg_features_matching_wt),
    TSC2_3BKO_vs_TSC2_WT = avg_features_matching_wt$TSC2_3BKO - avg_features_matching_wt$TSC2_WT,
    TSC1_TET3KO_vs_TSC1_WT = avg_features_matching_wt$TSC1_TET3KO - avg_features_matching_wt$TSC1_WT,
    TSC1_KDM2BKO_vs_TSC1_WT = avg_features_matching_wt$TSC1_KDM2BKO - avg_features_matching_wt$TSC1_WT,
    TSC1_RNF2KO_vs_TSC1_WT = avg_features_matching_wt$TSC1_RNF2KO - avg_features_matching_wt$TSC1_WT,
    TSC3_EEDKO_vs_TSC3_WT = avg_features_matching_wt$TSC3_EEDKO - avg_features_matching_wt$TSC3_WT,
    feature = avg_features_matching_wt$feature
)

diff_features_matching_wt_df <- melt(diff_features_matching_wt)
diff_features_matching_wt_df$variable <- factor(diff_features_matching_wt_df$variable, levels = rev(c("TSC2_3BKO_vs_TSC2_WT", "TSC1_TET3KO_vs_TSC1_WT", "TSC1_KDM2BKO_vs_TSC1_WT", "TSC1_RNF2KO_vs_TSC1_WT", "TSC3_EEDKO_vs_TSC3_WT")))

pdf("Extended_Data_Figure5b.pdf", height = 10, width = 12)
ggplot(subset(diff_feature_df, feature != "CGI"), aes(y=variable, x=value, fill = feature)) + geom_density_ridges(aes(fill = feature), rel_min_height=.01, alpha = 0.5, panel_scaling = FALSE, scale = 1) + coord_cartesian(xlim = c(-1, 1)) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + ylab("") + xlab("Delta mean methylation") + scale_fill_manual(values = brewer.pal(3, "Paired")) + geom_vline(xintercept = 0, lty = 2, col = "firebrick") + geom_vline(xintercept = c(-0.5, -0.25, 0.25, 0.5), lty = 2, col = "grey")

ggplot(subset(diff_features_matching_wt_df, feature != "CGI"), aes(y=variable, x=value, fill = feature)) + geom_density_ridges(aes(fill = feature), rel_min_height=.01, alpha = 0.5, panel_scaling = FALSE, scale = 1) + coord_cartesian(xlim = c(-1, 1)) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + ylab("") + xlab("Delta mean methylation") + scale_fill_manual(values = brewer.pal(3, "Paired")) + geom_vline(xintercept = 0, lty = 2, col = "firebrick") + geom_vline(xintercept = c(-0.5, -0.25, 0.25, 0.5), lty = 2, col = "grey")
dev.off()

###################################################################################################################
## Extended Data Figure 5d
###################################################################################################################

dna_modifications <- read.table(file = "Extended_Data_Figure5d_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1, colClasses = c("character", "character", rep("numeric", 24)))

dna_modifications$TSC_DNMT3B_KO_1 <- apply(dna_modifications[,c("TSC_DNMT3B_KO_1_1", "TSC_DNMT3B_KO_1_2", "TSC_DNMT3B_KO_1_3")], 1, mean)
dna_modifications$TSC_DNMT3B_KO_2 <- apply(dna_modifications[,c("TSC_DNMT3B_KO_2_1", "TSC_DNMT3B_KO_2_2", "TSC_DNMT3B_KO_2_3")], 1, mean)

dna_modifications$TSC_TET3_KO_1 <- apply(dna_modifications[,c("TSC_TET3_KO_1_1", "TSC_TET3_KO_1_2", "TSC_TET3_KO_1_3")], 1, mean)
dna_modifications$TSC_TET3_KO_2 <- apply(dna_modifications[,c("TSC_TET3_KO_2_1", "TSC_TET3_KO_2_2", "TSC_TET3_KO_2_3")], 1, mean)

dna_modifications$TSC_WT_1 <- apply(dna_modifications[,c("TSC_WT_1_1", "TSC_WT_1_2", "TSC_WT_1_3")], 1, mean)
dna_modifications$TSC_WT_2 <- apply(dna_modifications[,c("TSC_WT_2_1", "TSC_WT_2_2", "TSC_WT_2_3")], 1, mean)

dna_modifications$ESC_WT_1 <- apply(dna_modifications[,c("ESC_WT_1_1", "ESC_WT_1_2", "ESC_WT_1_3")], 1, mean)
dna_modifications$ESC_WT_2 <- apply(dna_modifications[,c("ESC_WT_2_1", "ESC_WT_2_2", "ESC_WT_2_3")], 1, mean)

rownames(dna_modifications) <- dna_modifications$Compound
dna_modifications <- dna_modifications[,-c(1:25)]

## Normalize by thymidine
dna_modifications <- data.frame(t(dna_modifications), check.names = FALSE)
dna_modifications <- dna_modifications / dna_modifications$Thymidine
dna_modifications$type <- gsub("_1|_2|_3|", "", rownames(dna_modifications))

dna_modifications_5mc_5hmc_df <- melt(dna_modifications[,c("5-hmdC", "5-mdC", "type")])
dna_modifications_5mc_5hmc_df$type <- factor(dna_modifications_5mc_5hmc_df$type, levels = c("ESC_WT", "TSC_WT", "TSC_TET3_KO", "TSC_DNMT3B_KO"))

pdf("Extended_Data_Figure5d.pdf", width = 5, height = 6)
ggplot(data = dna_modifications_5mc_5hmc_df, aes(x = type, y = value)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = type), position = position_jitterdodge(seed = 42), size = 2) + theme_classic() + xlab("") + ylab("Intensity normalized to thymidine") + facet_grid(variable~., scales = "free") + expand_limits(y = 0) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1))
dev.off()

###################################################################################################################
## Extended Data Figure 5h
###################################################################################################################

tiles_log2 <- list(diff_3BKO = melt(avg_tile_1kb_chip_wgbs[,c("log2_ratio_3BKO_H3K27me3", "log2_ratio_3BKO_H2AK119ub", "log2_ratio_3BKO_H3K4me3")]),
                   diff_KDM2BKO = melt(avg_tile_1kb_chip_wgbs[,c("log2_ratio_KDM2BKO_H3K27me3", "log2_ratio_KDM2BKO_H2AK119ub", "log2_ratio_KDM2BKO_H3K4me3")]),
                   diff_RNF2KO = melt(avg_tile_1kb_chip_wgbs[,c("log2_ratio_RNF2KO_H3K27me3", "log2_ratio_RNF2KO_H2AK119ub", "log2_ratio_RNF2KO_H3K4me3")]),
                   diff_EEDKO = melt(avg_tile_1kb_chip_wgbs[,c("log2_ratio_EEDKO_H3K27me3", "log2_ratio_EEDKO_H2AK119ub", "log2_ratio_EEDKO_H3K4me3")]))

tiles_log2$diff_3BKO$variable <- factor(gsub("log2_ratio_3BKO_", "", as.character(tiles_log2$diff_3BKO$variable)), levels = c("H3K27me3", "H2AK119ub", "H3K4me3"))
tiles_log2$diff_KDM2BKO$variable <- factor(gsub("log2_ratio_KDM2BKO_", "", as.character(tiles_log2$diff_KDM2BKO$variable)), levels = c("H3K27me3", "H2AK119ub", "H3K4me3"))
tiles_log2$diff_RNF2KO$variable <- factor(gsub("log2_ratio_RNF2KO_", "", as.character(tiles_log2$diff_RNF2KO$variable)), levels = c("H3K27me3", "H2AK119ub", "H3K4me3"))
tiles_log2$diff_EEDKO$variable <- factor(gsub("log2_ratio_EEDKO_", "", as.character(tiles_log2$diff_EEDKO$variable)), levels = c("H3K27me3", "H2AK119ub", "H3K4me3"))

pdf("Extended_Data_Figure5h.pdf", width = 14, height = 8)
layout(matrix(1:2, ncol = 2, byrow = TRUE))
vioplot(value ~ variable, data = tiles_log2$diff_3BKO, xlab = "", ylab = "Log2 ratio KO vs WT", main = "3BKO", ylim = c(-6, 6))
abline(h = 0, lty = 2)
abline(h = -1, lty = 2)
abline(h = 1, lty = 2)
vioplot(value ~ variable, data = tiles_log2$diff_KDM2BKO, xlab = "", ylab = "Log2 ratio KO vs WT", main = "KDM2BKO", ylim = c(-6, 6))
abline(h = 0, lty = 2)
abline(h = -1, lty = 2)
abline(h = 1, lty = 2)
vioplot(value ~ variable, data = tiles_log2$diff_RNF2KO, xlab = "", ylab = "Log2 ratio KO vs WT", main = "RNF2KO", ylim = c(-6, 6))
abline(h = 0, lty = 2)
abline(h = -1, lty = 2)
abline(h = 1, lty = 2)
vioplot(value ~ variable, data = tiles_log2$diff_EEDKO, xlab = "", ylab = "Log2 ratio KO vs WT", main = "EEDKO", ylim = c(-6, 6))
abline(h = 0, lty = 2)
abline(h = -1, lty = 2)
abline(h = 1, lty = 2)
dev.off()

###################################################################################################################
## Extended Data Figure 6a
###################################################################################################################

pdf("Extended_Data_Figure6a.pdf", width = 15, height = 20)
layout(matrix(1:12, ncol = 3, byrow = TRUE))
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_H3K27me3_Input, avg_tile_1kb_chip_wgbs$log2_ratio_3BKO_H3K27me3_Input, ylab = "3BKO", xlab = "WT", nrpoints = 0, main = "H3K27me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_H2AK119ub_Input, avg_tile_1kb_chip_wgbs$log2_ratio_3BKO_H2AK119ub_Input, ylab = "3BKO", xlab = "WT", nrpoints = 0, main = "H2AK119ub", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_H3K4me3_Input, avg_tile_1kb_chip_wgbs$log2_ratio_3BKO_H3K4me3_Input, ylab = "3BKO", xlab = "WT", nrpoints = 0, main = "H3K4me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)

smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H3K27me3_Input, avg_tile_1kb_chip_wgbs$log2_ratio_KDM2BKO_H3K27me3_Input, ylab = "KDM2BKO", xlab = "WT", nrpoints = 0, main = "H3K27me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H2AK119ub_Input, avg_tile_1kb_chip_wgbs$log2_ratio_KDM2BKO_H2AK119ub_Input, ylab = "KDM2BKO", xlab = "WT", nrpoints = 0, main = "H2AK119ub", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H3K4me3_Input, avg_tile_1kb_chip_wgbs$log2_ratio_KDM2BKO_H3K4me3_Input, ylab = "KDM2BKO", xlab = "WT", nrpoints = 0, main = "H3K4me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)

smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H3K27me3_Input, avg_tile_1kb_chip_wgbs$log2_ratio_RNF2KO_H3K27me3_Input, ylab = "RNF2KO", xlab = "WT", nrpoints = 0, main = "H3K27me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H2AK119ub_Input, avg_tile_1kb_chip_wgbs$log2_ratio_RNF2KO_H2AK119ub_Input, ylab = "RNF2KO", xlab = "WT", nrpoints = 0, main = "H2AK119ub", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H3K4me3_Input, avg_tile_1kb_chip_wgbs$log2_ratio_RNF2KO_H3K4me3_Input, ylab = "RNF2KO", xlab = "WT", nrpoints = 0, main = "H3K4me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)

smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H3K27me3_Input, avg_tile_1kb_chip_wgbs$log2_ratio_EEDKO_H3K27me3_Input, ylab = "EEDKO", xlab = "WT", nrpoints = 0, main = "H3K27me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H2AK119ub_Input, avg_tile_1kb_chip_wgbs$log2_ratio_EEDKO_H2AK119ub_Input, ylab = "EEDKO", xlab = "WT", nrpoints = 0, main = "H2AK119ub", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_WT_pool2_H3K4me3_Input, avg_tile_1kb_chip_wgbs$log2_ratio_EEDKO_H3K4me3_Input, ylab = "EEDKO", xlab = "WT", nrpoints = 0, main = "H3K4me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
dev.off()

###################################################################################################################
## Extended Data Figure 6b
###################################################################################################################

avg_cgi <- read.table(file = "Extended_Data_Figure6b_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

diff_avg_cgi <- data.frame(row.names = rownames(avg_cgi),
                           TSC1_KDM2BKO = avg_cgi$TSC1_KDM2BKO - avg_cgi$TSC1_WT,
                           TSC1_RNF2KO = avg_cgi$TSC1_RNF2KO - avg_cgi$TSC1_WT,
                           TSC3_EEDKO = avg_cgi$TSC3_EEDKO - avg_cgi$TSC1_WT)

hyper_cgis_prc_0.2 <- apply(diff_avg_cgi[,c("TSC1_KDM2BKO", "TSC1_RNF2KO", "TSC3_EEDKO")], 2, function(x) rownames(diff_avg_cgi)[which(x > 0.2)])
union_hyper_cgis_prc <- unique(do.call(c, hyper_cgis_prc_0.2))
union_hyper_cgis_prc_exe_hyper_cgi <- list(PRC = union_hyper_cgis_prc, ExE_hyper = rownames(subset(avg_cgi, feature == "exe_hyper_cgi")))

avg_cgi_prc_hyper_df <- melt(avg_cgi[union_hyper_cgis_prc,c("TSC1_WT", "TSC1_KDM2BKO", "TSC1_RNF2KO", "TSC3_EEDKO")])
avg_cgi_prc_hyper_df$variable <- factor(avg_cgi_prc_hyper_df$variable, levels = c("TSC1_WT", "TSC1_KDM2BKO", "TSC1_RNF2KO", "TSC3_EEDKO"))

pdf("Extended_Data_Figure6b.pdf")
plot(Venn(hyper_cgis_prc_0.2), doWeights = TRUE)
plot(Venn(union_hyper_cgis_prc_exe_hyper_cgi), doWeights = TRUE)

par(mar=c(16, 4, 4, 2))
vioplot(value ~ variable, data = avg_cgi_prc_hyper_df, xlab = "", ylab = "Mean methylation", ylim = c(0,1), main = "PRC hyper CGIs", las = 2)
dev.off()

###################################################################################################################
## Extended Data Figure 6e
###################################################################################################################

pdf("Extended_Data_Figure6e.pdf", width = 17, height = 5)
p1 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = TSC1_WT, y = log2_TSC1_WT_pool2_H3K4me3, color = log2_cut_TSC1_WT_pool2_H3K27me3)) + rasterise(geom_point(size = 0.5), dpi = 500) + theme_classic() + scale_color_gradientn(colors = plasma(10)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 10)) + xlab("Average DNAme (WT)") + ylab("Average RPGC H3K4me3 (WT)") + ggtitle("TSC WT") + theme(legend.position="bottom")
p2 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = TSC1_KDM2BKO, y = log2_TSC1_KDM2BKO_H3K4me3, color = log2_cut_TSC1_KDM2BKO_H3K27me3)) + rasterise(geom_point(size = 0.5), dpi = 500) + theme_classic() + scale_color_gradientn(colors = plasma(10)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 10)) + xlab("Average DNAme (KDM2BKO)") + ylab("Average RPGC H3K4me3 (KDM2BKO)") + ggtitle("TSC KDM2BKO") + theme(legend.position="bottom")
p3 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = TSC1_RNF2KO, y = log2_TSC1_RNF2KO_H3K4me3, color = log2_cut_TSC1_RNF2KO_H3K27me3)) + rasterise(geom_point(size = 0.5), dpi = 500) + theme_classic() + scale_color_gradientn(colors = plasma(10)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 10)) + xlab("Average DNAme (RNF2KO)") + ylab("Average RPGC H3K4me3 (RNF2KO)") + ggtitle("TSC RNF2KO") + theme(legend.position="bottom")
p4 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = TSC3_EEDKO, y = log2_TSC3_EEDKO_H3K4me3, color = log2_cut_TSC3_EEDKO_H3K27me3)) + rasterise(geom_point(size = 0.5), dpi = 500) + theme_classic() + scale_color_gradientn(colors = plasma(10)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 10)) + xlab("Average DNAme (EEDKO)") + ylab("Average RPGC H3K4me3 (EEDKO)") + ggtitle("TSC EEDKO") + theme(legend.position="bottom")
grid.arrange(p1,p2,p3,p4, ncol = 4)

p1 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = log2_TSC1_WT_pool2_H3K4me3)) + theme_classic() + xlab("Average H3K4me3 (WT)") + ylab("Fraction") + ggtitle("TSC WT") + theme(legend.position="bottom") + stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.5) + coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.3))
p2 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = log2_TSC1_KDM2BKO_H3K4me3)) + theme_classic() + xlab("Average H3K4me3 (KDM2BKO)") + ylab("Fraction") + ggtitle("TSC KDM2BKO") + theme(legend.position="bottom")+ stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.5) + coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.3))
p3 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = log2_TSC1_RNF2KO_H3K4me3)) + theme_classic() + xlab("Average H3K4me3 (RNF2KO)") + ylab("Fraction") + ggtitle("TSC RNF2KO") + theme(legend.position="bottom")+ stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.5) + coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.3))
p4 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = log2_TSC3_EEDKO_H3K4me3)) + theme_classic() + xlab("Average H3K4me3 (EEDKO)") + ylab("Fraction") + ggtitle("TSC EEDKO") + theme(legend.position="bottom")+ stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.5) + coord_cartesian(xlim = c(0, 10), ylim = c(0, 0.3))
grid.arrange(p1,p2,p3,p4, ncol = 4)

p1 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = TSC1_WT)) + theme_classic() + xlab("Average DNAme (WT)") + ylab("Fraction") + ggtitle("TSC WT") + theme(legend.position="bottom") + stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.05) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.4))
p2 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = TSC1_KDM2BKO)) + theme_classic() + xlab("Average DNAme (KDM2BKO)") + ylab("Fraction") + ggtitle("TSC KDM2BKO") + theme(legend.position="bottom")+ stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.05) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.4))
p3 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = TSC1_RNF2KO)) + theme_classic() + xlab("Average DNAme (RNF2KO)") + ylab("Fraction") + ggtitle("TSC RNF2KO") + theme(legend.position="bottom")+ stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.05) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.4))
p4 <- ggplot(data = avg_cgi_chip_wgbs, aes(x = TSC3_EEDKO)) + theme_classic() + xlab("Average DNAme (EEDKO)") + ylab("Fraction") + ggtitle("TSC EEDKO") + theme(legend.position="bottom")+ stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.05) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.4))
grid.arrange(p1,p2,p3,p4, ncol = 4)
dev.off()
