###################################################################################################################
## This script contains the code to produce Figure 2 and 3 and Extended Data Figures 3 and 4
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
library(rtracklayer)
library(viridis)
library(gridExtra)
library(karyoploteR)
library(xlsx)
library(plyr)

###################################################################################################################
## Figure 2a
###################################################################################################################

# hicPlotMatrix -m \\
# Figure2a_source_data_matrix.h5 \\
# --vMin -4 --vMax 4 \\
# --region chr1 \\
# --dpi 300 \\
# --colorMap "RdBu_r" \\
# --bigwig Figure2a_source_data_compartments.bw \\
# --bigwigAdditionalVerticalAxis \\
# --vMinBigwig -1 \\
# --vMaxBigwig 1 \\
# --increaseFigureHeight 3 \\
# -o Figure2a.pdf

###################################################################################################################
## Figure 2b
###################################################################################################################

df_ab_interaction_ratio <- read.table(file = "Figure2b.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

pdf("Figure2b.pdf", width = 5)
ggplot(df_ab_interaction_ratio, aes(x=compartment, y=ab_interaction_ratio, fill = cell_type)) + geom_boxplot() + theme_classic() + xlab("Compartment") + ylab("A/B interaction ratio") + theme(legend.position="bottom") + coord_cartesian(ylim = c(-6, 6)) + geom_hline(yintercept = 0, lty = 2, size = 0.5, col = "darkgrey")
dev.off()

wilcox.test(subset(df_ab_interaction_ratio, cell_type == "ESC" & compartment == "A")$ab_interaction_ratio, subset(df_ab_interaction_ratio, cell_type == "TSC" & compartment == "A")$ab_interaction_ratio) ## 0.0177
wilcox.test(subset(df_ab_interaction_ratio, cell_type == "ESC" & compartment == "B")$ab_interaction_ratio, subset(df_ab_interaction_ratio, cell_type == "TSC" & compartment == "B")$ab_interaction_ratio) ## < 2.2e-16

###################################################################################################################
## Figure 2c
###################################################################################################################

merged_compartments <- read.table(file = "Figure2c_2d_Extended_Data_Figure3a_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

pdf("Figure2c.pdf")
smoothScatter(merged_compartments$mESC_WT, merged_compartments$TSC1_WT, nrpoints = 0, xlab = "mESC", ylab = "TSC", xlim = c(-0.06, 0.06), ylim = c(-0.06, 0.06))
points(subset(merged_compartments, hyper_cgi == "yes")$mESC_WT, subset(merged_compartments, hyper_cgi == "yes")$TSC1_WT, pch = 19, col = brewer.pal(3, "Paired")[3], cex = 0.3)
abline(0, 1, lty = 2)
dev.off()

###################################################################################################################
## Figure 2d
###################################################################################################################

merged_compartments_hyper_cgis_df <- melt(subset(merged_compartments, hyper_cgi == "yes")[,c("TSC1_WT", "mESC_WT")])
merged_compartments_hyper_cgis_df$variable <- factor(merged_compartments_hyper_cgis_df$variable, levels = c("mESC_WT", "TSC1_WT"))

set.seed(42)
pdf("Figure2d.pdf", width = 6)
ggplot(merged_compartments_hyper_cgis_df, aes(x=variable, y=value)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2)) + theme_classic() + xlab("") + ylab("PC1") + geom_hline(yintercept = 0, lty = 2, col = "grey")
dev.off()

wilcox.test(subset(merged_compartments, hyper_cgi == "yes")[,"mESC_WT"], subset(merged_compartments, hyper_cgi == "yes")[,"TSC1_WT"])

###################################################################################################################
## Figure 2e
###################################################################################################################

figure2e_source_data <- read.table(file = "Figure2e_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

df_distribution_exe_hyper_cgis_compartments_esc <- data.frame(table(figure2e_source_data[,c("compartment_esc_wt")]))
df_distribution_exe_hyper_cgis_compartments_tsc <- data.frame(table(figure2e_source_data[,c("compartment_tsc_wt")]))

df_distribution_exe_hyper_cgis_compartments_esc$cell_type <- "ESC"
df_distribution_exe_hyper_cgis_compartments_tsc$cell_type <- "TSC"

df_distribution_exe_hyper_cgis_compartments <- rbind(df_distribution_exe_hyper_cgis_compartments_esc, df_distribution_exe_hyper_cgis_compartments_tsc)

pdf("Figure2e.pdf", width = 5)
ggplot(data = df_distribution_exe_hyper_cgis_compartments, aes(x = cell_type, y = Freq, fill = Var1)) + geom_bar(position = "fill", stat = "identity", alpha = 0.6) + theme_classic() + xlab("") + ylab("Fraction hyper CGIs")
dev.off()

contingency_table <- matrix(df_distribution_exe_hyper_cgis_compartments$Freq, ncol = 2, nrow = 2, byrow = TRUE)

chisq.test(contingency_table) ## 0.2063

###################################################################################################################
## Figure 2f
###################################################################################################################

## Figure 2f was generated using the processed methylation rates and smoothed MINUTE-ChIP tracks of TSC1
## and ESC WT visualized in IGV.

###################################################################################################################
## Figure 2g
###################################################################################################################

avg_tile_1kb_chip_wgbs <- read.table(file = "Figure2g_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

pdf("Figure2g.pdf", width = 15, height = 5)
layout(matrix(1:3, ncol = 3))
smoothScatter(avg_tile_1kb_chip_wgbs$diff_TSC_ESC, avg_tile_1kb_chip_wgbs$log2_ratio_TSC_ESC_H3K27me3, xlab = "Delta DNAme (TSC vs ESC)", ylab = "Log2 ratio H3K27me3 (TSC vs ESC)", nrpoints = 0, main = "H3K27me3", ylim = c(-5, 5), xlim = c(-1, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_TSC_ESC"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_TSC_ESC_H3K27me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_TSC_ESC, avg_tile_1kb_chip_wgbs$log2_ratio_TSC_ESC_H2AK119ub, xlab = "Delta DNAme (TSC vs ESC)", ylab = "Log2 ratio H2AK119ub (TSC vs ESC)", nrpoints = 0, main = "H2AK119ub", ylim = c(-5, 5), xlim = c(-1, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_TSC_ESC"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_TSC_ESC_H2AK119ub"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$diff_TSC_ESC, avg_tile_1kb_chip_wgbs$log2_ratio_TSC_ESC_H3K4me3, xlab = "Delta DNAme (TSC vs ESC)", ylab = "Log2 ratio H3K4me3 (TSC vs ESC)", nrpoints = 0, main = "H3K4me3", ylim = c(-5, 5), xlim = c(-1, 1))
points(subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"diff_TSC_ESC"], subset(avg_tile_1kb_chip_wgbs, overlap_exe_hyper_cgi == "yes")[,"log2_ratio_TSC_ESC_H3K4me3"], pch = 19, col = "#B2DF8A", cex = 0.3)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
dev.off()

###################################################################################################################
## Figure 3a
###################################################################################################################

figure_3a_source_data <- read.table(file = "Figure3a_source_data.tsv", sep = "\\t", header = TRUE, stringsAsFactors = FALSE)

samples <- unique(figure_3a_source_data$sample)

norm_mat_list <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(figure_3a_source_data, sample == x)[,!colnames(figure_3a_source_data) %in% "sample"]), k_upstream=100, k_downstream=100, k_target=67, extend = c(5000, 5000), signal_name = x, target_name = "exe_hyper_cgis"))
names(norm_mat_list) <- samples

col_fun <- list(
    "mESC_WT_H3K27me3_Thermo_Mnase_BS_merged" = colorRamp2(c(0, 40), c("white", "violetred3")),
    "TSC1_WT_H3K27me3_Thermo_Mnase_BS_merged" = colorRamp2(c(0, 20), c("white", "violetred3")),
    "mESC_WT_EED_ab240650_merged" = colorRamp2(c(0, 200), c("white", "mediumpurple1")),
    "TSC1_WT_EED_ab240650_merged" = colorRamp2(c(0, 40), c("white", "mediumpurple1")),
    "ChIP_BS_DNAme_mESC_WT_H3K27me3_Thermo_Mnase_BS_merged" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    "ChIP_BS_DNAme_TSC1_WT_H3K27me3_Thermo_Mnase_BS_merged" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    )

meta_height <- c(
    "mESC_WT_H3K27me3_Thermo_Mnase_BS_merged" = 30,
    "TSC1_WT_H3K27me3_Thermo_Mnase_BS_merged" = 5,
    "mESC_WT_EED_ab240650_merged" = 60,
    "TSC1_WT_EED_ab240650_merged" = 10,
    "ChIP_BS_DNAme_mESC_WT_H3K27me3_Thermo_Mnase_BS_merged" = 1,
    "ChIP_BS_DNAme_TSC1_WT_H3K27me3_Thermo_Mnase_BS_merged" = 1
    )

ht_list <- NULL
for(sample in samples)
{
    ht_list <- ht_list + EnrichedHeatmap(norm_mat_list[[sample]], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 4))
}

pdf("Figure3a.pdf", width = 24)
draw(ht_list)
dev.off()

###################################################################################################################
## Figure 3b
###################################################################################################################

## Figure 3b was generated using the processed methylation rates, H3K27me3 ChIP-BS-seq and EED ChIPseq
## tracks of TSC1 and ESC WT visualized in IGV. Single reads (output from RLM) were colored by the average
## methylation rates and also visualized using IGV.

###################################################################################################################
## Figure 3c
###################################################################################################################

avg_exe_hyper_cgis_wgbs_chip_bs <- read.table(file = "Figure3c_Extended_Data_Figure4b_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

pdf("Figure3c.pdf", width = 8.5, height = 10)
p1 <- ggplot(data = avg_exe_hyper_cgis_wgbs_chip_bs, aes(x = WGBS_mESC_WT, y = ChIP_BS_mESC_WT_merged, color = log2_cut_mESC_WT_H3K27me3_merged)) + geom_point(size = 0.5) + theme_classic() + scale_color_gradientn(colors = plasma(10)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + xlab("WGBS ESC") + ylab("ChIP-BS ESC") + geom_abline(intercept = 0, slope = 1, lty = 2)
p2 <- ggplot(data = avg_exe_hyper_cgis_wgbs_chip_bs, aes(x = WGBS_TSC1_WT, y = ChIP_BS_TSC1_WT_merged, color = log2_cut_TSC1_WT_H3K27me3_merged)) + geom_point(size = 0.5) + theme_classic() + scale_color_gradientn(colors = plasma(10)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + xlab("WGBS TSC") + ylab("ChIP-BS TSC") + geom_abline(intercept = 0, slope = 1, lty = 2)
grid.arrange(p1,p2, ncol = 1)

p1 <- ggplot(data = avg_exe_hyper_cgis_wgbs_chip_bs, aes(x = WGBS_mESC_WT)) + theme_classic() + xlab("WGBS ESC") + ylab("Fraction") + theme(legend.position="bottom") + stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.05) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.4))
p2 <- ggplot(data = avg_exe_hyper_cgis_wgbs_chip_bs, aes(x = ChIP_BS_mESC_WT_merged)) + theme_classic() + xlab("ChIP-BS ESC") + ylab("Fraction") + theme(legend.position="bottom")+ stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.05) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.4))
grid.arrange(p1,p2, ncol = 1)

p1 <- ggplot(data = avg_exe_hyper_cgis_wgbs_chip_bs, aes(x = WGBS_TSC1_WT)) + theme_classic() + xlab("WGBS TSC") + ylab("Fraction") + theme(legend.position="bottom") + stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.05) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.4))
p2 <- ggplot(data = avg_exe_hyper_cgis_wgbs_chip_bs, aes(x = ChIP_BS_TSC1_WT_merged)) + theme_classic() + xlab("ChIP-BS TSC") + ylab("Fraction") + theme(legend.position="bottom")+ stat_bin(aes(y=..count../sum(..count..)), binwidth = 0.05) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.4))
grid.arrange(p1,p2, ncol = 1)
dev.off()

###################################################################################################################
## Extended Data Figure 3a
###################################################################################################################

compartment_types <- rep(NA, nrow(merged_compartments))
compartment_types[merged_compartments$compartment_esc_wt == "A" & merged_compartments$compartment_tsc_wt == "A"] <- "esc_A_tsc_A"
compartment_types[merged_compartments$compartment_esc_wt == "B" & merged_compartments$compartment_tsc_wt == "B"] <- "esc_B_tsc_B"
compartment_types[merged_compartments$compartment_esc_wt == "A" & merged_compartments$compartment_tsc_wt == "B"] <- "esc_A_tsc_B"
compartment_types[merged_compartments$compartment_esc_wt == "B" & merged_compartments$compartment_tsc_wt == "A"] <- "esc_B_tsc_A"

compartment_types_range <- makeGRangesFromDataFrame(merged_compartments[,1:3])
compartment_types_range$type <- compartment_types

colors_compartments <- brewer.pal(8, "Paired")[c(1,2,6,8)]
names(colors_compartments) <- c("esc_A_tsc_A", "esc_B_tsc_B", "esc_A_tsc_B", "esc_B_tsc_A")

pdf("Extended_Data_Figure3a.pdf", height = 16, width = 12)
pp <- getDefaultPlotParams(plot.type=2)
pp$topmargin <- 300
pp$data1height <- 180
kp <- plotKaryotype(genome="mm10", chromosomes=c("autosomal"), plot.type=2, plot.params = pp)
kpPlotRegions(kp, compartment_types_range, col = colors_compartments[compartment_types_range$type], r0 = 0, r1 = 0.5)
legend(x = "right", legend = names(colors_compartments), fill = colors_compartments)
kpPlotBigWig(kp, data="Figure2a_source_data_compartments.bw", ymax="visible.region", r0=1.3, r1=1.9)
kpPlotBigWig(kp, data="Extended_Data_Figure3b_source_data_compartments.bw", ymax="visible.region", r0=0.6, r1=1.2)
kpAddLabels(kp, labels = "mESC WT", r0=1.3, r1=1.9, cex=0.8, label.margin = 0.035)
kpAddLabels(kp, labels = "TSC WT", r0=0.6, r1=1.2, cex=0.8, label.margin = 0.035)
dev.off()

###################################################################################################################
## Extended Data Figure 3b
###################################################################################################################

# hicPlotMatrix -m \\
# Extended_Data_Figure3b_3c_source_data_matrix_ESC.h5 \\
# --vMin 1 --vMax 100000 \\
# --region chr1 \\
# --log1p \\
# --dpi 300 \\
# --colorMap "Reds" \\
# --bigwig Figure2a_source_data_compartments.bw \\
# --bigwigAdditionalVerticalAxis \\
# --vMinBigwig -1 \\
# --vMaxBigwig 1 \\
# --increaseFigureHeight 3 \\
# -o Extended_Data_Figure3b_ESC.pdf
#
# hicPlotMatrix -m \\
# Extended_Data_Figure3b_3c_source_data_matrix_TSC.h5 \\
# --vMin 1 --vMax 100000 \\
# --region chr1 \\
# --log1p \\
# --dpi 300 \\
# --colorMap "Reds" \\
# --bigwig Extended_Data_Figure3b_source_data_compartments.bw \\
# --bigwigAdditionalVerticalAxis \\
# --vMinBigwig -1 \\
# --vMaxBigwig 1 \\
# --increaseFigureHeight 3 \\
# -o Extended_Data_Figure3b_TSC.pdf

###################################################################################################################
## Extended Data Figure 3c
###################################################################################################################

# hicPlotDistVsCounts -m \\
# Extended_Data_Figure3b_3c_source_data_matrix_ESC.h5 \\
# Extended_Data_Figure3b_3c_source_data_matrix_TSC.h5 \\
# -o Extended_Data_Figure3c.pdf \\
# --labels 'mESC_WT' 'TSC1_WT' \\
# --maxdepth 80000000 \\
# --chromosomeExclude chrX chrY \\
# --plotsize 7 4.2

###################################################################################################################
## Extended Data Figure 3d
###################################################################################################################

avg_tile_1kb_chip_wgbs <- read.table(file = "Extended_Data_Figure3d_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

pdf("Extended_Data_Figure3d.pdf", width = 20, height = 5)
layout(matrix(1:4, ncol = 4))
smoothScatter(avg_tile_1kb_chip_wgbs$mESC_WT, avg_tile_1kb_chip_wgbs$TSC1_WT, xlab = "ESC", ylab = "TSC", nrpoints = 0, main = "DNAme")
abline(0, 1, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_ESC_Input_H3K27me3, avg_tile_1kb_chip_wgbs$log2_ratio_TSC_Input_H3K27me3, xlab = "ESC", ylab = "TSC", nrpoints = 0, main = "H3K27me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_ESC_Input_H2AK119ub, avg_tile_1kb_chip_wgbs$log2_ratio_TSC_Input_H2AK119ub, xlab = "ESC", ylab = "TSC", nrpoints = 0, main = "H2AK119ub", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip_wgbs$log2_ratio_ESC_Input_H3K4me3, avg_tile_1kb_chip_wgbs$log2_ratio_TSC_Input_H3K4me3, xlab = "ESC", ylab = "TSC", nrpoints = 0, main = "H3K4me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
dev.off()

###################################################################################################################
## Extended Data Figure 3f
###################################################################################################################

histone_tails <- read.table(file = "Extended_Data_Figure3f_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)

histone_tails$ESC_WT_1 <- apply(histone_tails[,c("ESC_WT_1_1","ESC_WT_1_2")], 1, mean, na.rm = TRUE)
histone_tails$ESC_WT_2 <- apply(histone_tails[,c("ESC_WT_2_1","ESC_WT_2_2")], 1, mean, na.rm = TRUE)
histone_tails$TSC1_WT_PB3_1 <- apply(histone_tails[,c("TSC1_WT_PB3_1_1","TSC1_WT_PB3_1_2")], 1, mean, na.rm = TRUE)
histone_tails$TSC1_WT_PB3_2 <- apply(histone_tails[,c("TSC1_WT_PB3_2_1","TSC1_WT_PB3_2_2")], 1, mean, na.rm = TRUE)
histone_tails$TSC1_WT_IP <- apply(histone_tails[,c("TSC1_WT_IP_1","TSC1_WT_IP_2")], 1, mean, na.rm = TRUE)
histone_tails$ESC_WT <- apply(histone_tails[,c("ESC_WT_1","ESC_WT_2")], 1, mean, na.rm = TRUE)

histone_tails$lfc_TSC1_WT_PB3_1 <- log2(histone_tails$TSC1_WT_PB3_1 / histone_tails$ESC_WT)
histone_tails$lfc_TSC1_WT_PB3_2 <- log2(histone_tails$TSC1_WT_PB3_2 / histone_tails$ESC_WT)
histone_tails$lfc_TSC1_WT_IP <- log2(histone_tails$TSC1_WT_IP / histone_tails$ESC_WT)

histone_tails_df <- melt(histone_tails[,c("Compound", "lfc_TSC1_WT_PB3_1", "lfc_TSC1_WT_PB3_2", "lfc_TSC1_WT_IP")])

histone_tails_df_avg_sd <- ddply(histone_tails_df[order(histone_tails_df$Compound),], ~Compound, summarise, mean_lfc = mean(value), sd_lfc = sd(value))

pdf("Extended_Data_Figure3f.pdf", width = 14, height = 5)
ggplot(histone_tails_df_avg_sd) + geom_bar(aes(x=Compound, y=mean_lfc), stat="identity") + geom_errorbar(aes(x=Compound, ymin=mean_lfc-sd_lfc, ymax=mean_lfc+sd_lfc), width=0.4) + theme_light() + theme(axis.text.x=element_text(angle=90, hjust=1, size=12), axis.text.y=element_text(size=12), axis.title=element_text(size=14), strip.text.x=element_text(size=12)) + ylab("Log2 fold change vs ESC") + coord_cartesian(ylim = c(-2.5, 2.5)) + geom_point(data = histone_tails_df, aes(x = Compound, y = value, fill = "black"), position = position_jitterdodge(seed = 42), size = 1)
dev.off()

###################################################################################################################
## Extended Data Figure 3g
###################################################################################################################

extended_data_3g_source_data <- read.table(file = "Extended_Data_Figure3g_source_data_heatmap.tsv", sep = "\\t", header = TRUE, stringsAsFactors = FALSE)

samples <- unique(extended_data_3g_source_data$sample)

norm_mat_list_cgis <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(extended_data_3g_source_data, sample == x & feature == "CGI")[,!colnames(extended_data_3g_source_data) %in% c("sample", "feature")]), k_upstream=100, k_downstream=100, k_target=67, extend = c(5000, 5000), signal_name = x, target_name = "CGI"))
names(norm_mat_list_cgis) <- samples

norm_mat_list_exe_cgis <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(extended_data_3g_source_data, sample == x & feature == "exe_hyper_cgi")[,!colnames(extended_data_3g_source_data) %in% c("sample", "feature")]), k_upstream=100, k_downstream=100, k_target=67, extend = c(5000, 5000), signal_name = x, target_name = "exe_hyper_cgi"))
names(norm_mat_list_exe_cgis) <- samples

norm_mat_list_promoters <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(extended_data_3g_source_data, sample == x & feature == "promoter")[,!colnames(extended_data_3g_source_data) %in% c("sample", "feature")]), k_upstream=100, k_downstream=100, k_target=67, extend = c(5000, 5000), signal_name = x, target_name = "promoter"))
names(norm_mat_list_promoters) <- samples

norm_mat_list_iaps <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(extended_data_3g_source_data, sample == x & feature == "full_length_iap")[,!colnames(extended_data_3g_source_data) %in% c("sample", "feature")]), k_upstream=100, k_downstream=100, k_target=67, extend = c(5000, 5000), signal_name = x, target_name = "full_length_iap"))
names(norm_mat_list_iaps) <- samples

## Heatmaps
samples_esc_tsc <- c("mESC_WT_H3K27me3", "TSC1_WT_H3K27me3", "mESC_WT_H2AK119ub", "TSC1_WT_H2AK119ub", "mESC_WT_H3K4me3", "TSC1_WT_H3K4me3", "WGBS_mESC_WT", "WGBS_TSC1_WT", "CpG_density")

col_fun <- list(
    "mESC_WT_H3K27me3" = colorRamp2(c(0, 20), c("white", "violetred2")),
    "TSC1_WT_H3K27me3" = colorRamp2(c(0, 20), c("white", "violetred2")),
    "mESC_WT_H2AK119ub" = colorRamp2(c(0, 20), c("white", "deepskyblue3")),
    "TSC1_WT_H2AK119ub" = colorRamp2(c(0, 20), c("white", "deepskyblue3")),
    "mESC_WT_H3K4me3" = colorRamp2(c(0, 200), c("white", "seagreen3")),
    "TSC1_WT_H3K4me3" = colorRamp2(c(0, 200), c("white", "seagreen3")),
    "WGBS_mESC_WT" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    "WGBS_TSC1_WT" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    "CpG_density" = colorRamp2(c(0, 0.2), c("white", "black"))
    )

meta_height <- c(
    "mESC_WT_H3K27me3" = 20,
    "TSC1_WT_H3K27me3" = 20,
    "mESC_WT_H2AK119ub" = 20,
    "TSC1_WT_H2AK119ub" = 20,
    "mESC_WT_H3K4me3" = 100,
    "TSC1_WT_H3K4me3" = 100,
    "WGBS_mESC_WT" = 1,
    "WGBS_TSC1_WT" = 1,
    "CpG_density" = 0.2
    )

ht_list1 <- NULL
for(sample in samples_esc_tsc)
{
    ht_list1 <- ht_list1 + EnrichedHeatmap(norm_mat_list_exe_cgis[[sample]], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
}

ht_list2 <- NULL
for(sample in samples_esc_tsc)
{
    ht_list2 <- ht_list2 + EnrichedHeatmap(norm_mat_list_cgis[[sample]], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
}

samples_esc_tsc <- c("mESC_WT_H3K4me3", "TSC1_WT_H3K4me3", "mESC_WT_H3K27me3", "TSC1_WT_H3K27me3", "mESC_WT_H2AK119ub", "TSC1_WT_H2AK119ub", "WGBS_mESC_WT", "WGBS_TSC1_WT", "CpG_density")

ht_list3 <- NULL
for(sample in samples_esc_tsc)
{
    ht_list3 <- ht_list3 + EnrichedHeatmap(norm_mat_list_promoters[[sample]], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
}

col_fun <- list(
    "mESC_WT_H3K27me3" = colorRamp2(c(0, 20), c("white", "violetred2")),
    "TSC1_WT_H3K27me3" = colorRamp2(c(0, 20), c("white", "violetred2")),
    "mESC_WT_H2AK119ub" = colorRamp2(c(0, 20), c("white", "deepskyblue3")),
    "TSC1_WT_H2AK119ub" = colorRamp2(c(0, 20), c("white", "deepskyblue3")),
    "mESC_WT_H3K4me3" = colorRamp2(c(0, 20), c("white", "seagreen3")),
    "TSC1_WT_H3K4me3" = colorRamp2(c(0, 20), c("white", "seagreen3")),
    "WGBS_mESC_WT" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    "WGBS_TSC1_WT" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    "CpG_density" = colorRamp2(c(0, 0.2), c("white", "black"))
    )
meta_height <- c(
    "mESC_WT_H3K27me3" = 20,
    "TSC1_WT_H3K27me3" = 20,
    "mESC_WT_H2AK119ub" = 20,
    "TSC1_WT_H2AK119ub" = 20,
    "mESC_WT_H3K4me3" = 20,
    "TSC1_WT_H3K4me3" = 20,
    "WGBS_mESC_WT" = 1,
    "WGBS_TSC1_WT" = 1,
    "CpG_density" = 0.2
    )

samples_esc_tsc <- c("mESC_WT_H3K27me3", "TSC1_WT_H3K27me3", "mESC_WT_H2AK119ub", "TSC1_WT_H2AK119ub", "mESC_WT_H3K4me3", "TSC1_WT_H3K4me3", "WGBS_mESC_WT", "WGBS_TSC1_WT", "CpG_density")

ht_list4 <- NULL
for(sample in samples_esc_tsc)
{
    ht_list4 <- ht_list4 + EnrichedHeatmap(norm_mat_list_iaps[[sample]], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
}

pdf("Extended_Data_Figure3g_heatmap.pdf", width = 22, height = 5)
draw(ht_list2)
draw(ht_list1)
draw(ht_list3)
draw(ht_list4)
dev.off()

df_repeat_iap_int_esc_tsc <- read.table(file = "Extended_Data_Figure3g_source_data_repeat_expression.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)
df_repeat_iap_int_esc_tsc$sample <- factor(df_repeat_iap_int_esc_tsc$sample, levels = c("ESC_WT", "TSC_WT"))

pdf("Extended_Data_Figure3g_repeat_expression.pdf", width = 5, height = 6)
ggplot(df_repeat_iap_int_esc_tsc, aes(color=repeats, y=value, x=sample)) + geom_point(shape = 4, size = 3) + ylab("Reads overlapping repeats (normalized)") + xlab("") + theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1, size=12), axis.text.y=element_text(size=12), axis.title=element_text(size=14), strip.text.x=element_text(size=12)) + scale_color_manual(values = c(brewer.pal(9, "RdYlBu"), brewer.pal(8, "PRGn")))
dev.off()

###################################################################################################################
## Extended Data Figure 4a
###################################################################################################################

methyl_data <- data.frame(fread("Extended_Data_Figure4a_source_data.tsv"), stringsAsFactors = FALSE)

pdf("Extended_Data_Figure4a.pdf", width = 9, height = 5)
layout(matrix(1:2, ncol = 2, byrow = TRUE))
smoothScatter(methyl_data[,"WGBS_mESC_WT"], methyl_data[,"ChIP_BS_DNAme_mESC_WT_merged"], xlab = "ESC WGBS", ylab = "ESC ChIP-BS", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
smoothScatter(methyl_data[,"WGBS_TSC1_WT"], methyl_data[,"ChIP_BS_DNAme_TSC1_WT_merged"], xlab = "TSC WGBS", ylab = "TSC ChIP-BS", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
dev.off()

###################################################################################################################
## Extended Data Figure 4b
###################################################################################################################

## Vioplots ExE hyper CGIs
avg_exe_hyper_cgis_wgbs_chip_bs_df <- melt(avg_exe_hyper_cgis_wgbs_chip_bs[,c("WGBS_mESC_WT", "ChIP_BS_mESC_WT_merged", "WGBS_TSC1_WT", "ChIP_BS_TSC1_WT_merged")])
avg_exe_hyper_cgis_wgbs_chip_bs_df$variable <- factor(avg_exe_hyper_cgis_wgbs_chip_bs_df$variable, levels = c("WGBS_mESC_WT", "ChIP_BS_mESC_WT_merged", "WGBS_TSC1_WT", "ChIP_BS_TSC1_WT_merged"))

pdf("Extended_Data_Figure4b.pdf", width = 7, height = 8)
par(mar=c(16, 4, 4, 2))
vioplot(value ~ variable, data = avg_exe_hyper_cgis_wgbs_chip_bs_df, xlab = "", ylab = "Mean methylation", ylim = c(0,1), las = 2, cex.lab = 0.8)
dev.off()

###################################################################################################################
## Extended Data Figure 4c
###################################################################################################################

## Extended Data Figure 4c was generated using the processed methylation rates and smoothed ChIP-BS tracks
## of TSC1 and ESC WT visualized in IGV. Single reads (output from RLM) were colored by the average
## methylation rates and also visualized using IGV.
