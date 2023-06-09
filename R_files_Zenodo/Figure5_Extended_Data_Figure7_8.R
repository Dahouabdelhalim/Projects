###################################################################################################################
## This script contains the code to produce Figure 5 and Extended Data Figures 7 and 8
## Weigert, Hetzel et al. Dynamic antagonism between key repressive pathways maintains the placental epigenome 2023
## Author: Sara Hetzel
###################################################################################################################

library(ggplot2)
library(vioplot)
library(reshape2)
library(RColorBrewer)
library(WebGestaltR)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)
library(data.table)
library(ggridges)
library(ggrepel)
library(rtracklayer)
library(Vennerable)
theme_set(theme_ridges())

###################################################################################################################
## Figure 5a and Extended Data Figure 7c
###################################################################################################################

avg_features_dnmt1i_tsc1 <- read.table(file = "Figure5a_Extended_Data_Figure7c_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

samples_time_course_tsc1 <- c("TSC1_WT", "TSC1_WT_DNMT1i_1d", "TSC1_WT_DNMT1i_2d", "TSC1_WT_DNMT1i_3d", "TSC1_WT_DNMT1i_4d", "TSC1_WT_DNMT1i_5d", "TSC1_WT_DNMT1i_6d", "TSC1_WT_DMSO_7d", "TSC1_WT_DNMT1i_7d", "TSC1_WT_DNMT1i_7d_recovery_1d", "TSC1_WT_DNMT1i_7d_recovery_2d", "TSC1_WT_DNMT1i_7d_recovery_3d", "TSC1_WT_DNMT1i_7d_recovery_4d", "TSC1_WT_DNMT1i_7d_recovery_5d", "TSC1_WT_DNMT1i_7d_recovery_6d", "TSC1_WT_DMSO_7d_recovery_7d", "TSC1_WT_DNMT1i_7d_recovery_7d", "TSC1_WT_DMSO_7d_recovery_14d",  "TSC1_WT_DNMT1i_7d_recovery_14d", "TSC1_WT_DMSO_7d_recovery_28d", "TSC1_WT_DNMT1i_7d_recovery_28d")

avg_features_dnmt1i_tsc1_df <- melt(avg_features_dnmt1i_tsc1)
avg_features_dnmt1i_tsc1_df$variable <- factor(avg_features_dnmt1i_tsc1_df$variable, levels =samples_time_course_tsc1)
avg_features_dnmt1i_tsc1_df$treatment <- sapply(as.character(avg_features_dnmt1i_tsc1_df$variable), function(x) strsplit(x, "TSC1_WT_|_[1-9]d")[[1]][2])
avg_features_dnmt1i_tsc1_df$treatment[is.na(avg_features_dnmt1i_tsc1_df$treatment)] <- "untreated"
avg_features_dnmt1i_tsc1_df$treatment <- factor(avg_features_dnmt1i_tsc1_df$treatment, levels = c("untreated", "DMSO", "DNMT1i"))
avg_features_dnmt1i_tsc1_df$feature <- factor(avg_features_dnmt1i_tsc1_df$feature, levels = c("HMD", "PMD", "exe_hyper_cgi"))

pdf("Figure5a_Extended_Data_Figure7c.pdf", height = 7, width = 18)
ggplot(avg_features_dnmt1i_tsc1_df, aes(x=variable, y=value, fill = feature)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 1)) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + xlab("") + ylab("Mean methylation") + scale_fill_manual(values = c(brewer.pal(3, "Paired")[c(1,2,3)])) + facet_grid(~treatment, scales = "free_x", space = "free_x") + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1))
dev.off()

###################################################################################################################
## Figure 5b and Extended Data Figure 8a
###################################################################################################################

avg_features_ezh2i_long_time_course <- read.table(file = "Figure5b_Extended_Data_Figure8a_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_features_ezh2i_long_time_course_df <- melt(avg_features_ezh2i_long_time_course)

avg_features_ezh2i_long_time_course_df$experiment <- ifelse(grepl("E4", as.character(avg_features_ezh2i_long_time_course_df$variable)), "EZH2i_DNMT1i", "EZH2i")
avg_features_ezh2i_long_time_course_df$type <- ifelse(grepl("DMSO", as.character(avg_features_ezh2i_long_time_course_df$variable)), "control", "treatment")

pdf("Figure5b_Extended_Data_Figure8a.pdf", height = 8, width = 12)
ggplot(avg_features_ezh2i_long_time_course_df, aes(x=variable, y=value, fill = feature)) + geom_boxplot(outlier.shape = NA) + theme_classic() + theme(axis.text=element_text(size=10), axis.text.x=element_text(hjust = 1, angle = 45), axis.title=element_text(size=12)) + ylab("") + xlab("Mean methylation") + scale_fill_manual(values = brewer.pal(3, "Paired")) + facet_grid(type~experiment, scales = "free", space = "free_x")
dev.off()

###################################################################################################################
## Figure 5c
###################################################################################################################

avg_tile_1kb_chip <- read.table(file = "Figure5c_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_tile_1kb_chip$log2_ratio_TSC1_WT_H2AK119ub_Input <- log2((avg_tile_1kb_chip$TSC1_WT_H2AK119ub + 1) / (avg_tile_1kb_chip$TSC1_WT_Input + 1))
avg_tile_1kb_chip$log2_ratio_TSC1_WT_H3K27me3_Input <- log2((avg_tile_1kb_chip$TSC1_WT_H3K27me3 + 1) / (avg_tile_1kb_chip$TSC1_WT_Input + 1))
avg_tile_1kb_chip$log2_ratio_TSC1_WT_H3K4me3_Input <- log2((avg_tile_1kb_chip$TSC1_WT_H3K4me3 + 1) / (avg_tile_1kb_chip$TSC1_WT_Input + 1))
avg_tile_1kb_chip$log2_ratio_TSC1_WT_EZH2i_5w_recovery_4w_H2AK119ub_Input <- log2((avg_tile_1kb_chip$TSC1_WT_EZH2i_5w_recovery_4w_H2AK119ub + 1) / (avg_tile_1kb_chip$TSC1_WT_EZH2i_5w_recovery_4w_Input + 1))
avg_tile_1kb_chip$log2_ratio_TSC1_WT_EZH2i_5w_recovery_4w_H3K27me3_Input <- log2((avg_tile_1kb_chip$TSC1_WT_EZH2i_5w_recovery_4w_H3K27me3 + 1) / (avg_tile_1kb_chip$TSC1_WT_EZH2i_5w_recovery_4w_Input + 1))
avg_tile_1kb_chip$log2_ratio_TSC1_WT_EZH2i_5w_recovery_4w_H3K4me3_Input <- log2((avg_tile_1kb_chip$TSC1_WT_EZH2i_5w_recovery_4w_H3K4me3 + 1) / (avg_tile_1kb_chip$TSC1_WT_EZH2i_5w_recovery_4w_Input + 1))

pdf("Figure5c.pdf", width = 15, height = 5)
layout(matrix(1:3, ncol = 3, byrow = TRUE))
smoothScatter(avg_tile_1kb_chip$log2_ratio_TSC1_WT_H3K27me3_Input, avg_tile_1kb_chip$log2_ratio_TSC1_WT_EZH2i_5w_recovery_4w_H3K27me3_Input, ylab = "EZH2i 5w recovery 4w", xlab = "WT", nrpoints = 0, main = "H3K27me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip$log2_ratio_TSC1_WT_H2AK119ub_Input, avg_tile_1kb_chip$log2_ratio_TSC1_WT_EZH2i_5w_recovery_4w_H2AK119ub_Input, ylab = "EZH2i 5w recovery 4w", xlab = "WT", nrpoints = 0, main = "H2AK119ub1", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
smoothScatter(avg_tile_1kb_chip$log2_ratio_TSC1_WT_H3K4me3_Input, avg_tile_1kb_chip$log2_ratio_TSC1_WT_EZH2i_5w_recovery_4w_H3K4me3_Input, ylab = "EZH2i 5w recovery 4w", xlab = "WT", nrpoints = 0, main = "H3K4me3", xlim = c(-2,8), ylim = c(-2,8))
abline(0, 1, lty = 2)
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)
dev.off()

###################################################################################################################
## Figure 5e
###################################################################################################################

tsc_mass_spec_results <- read.table(file = "Figure5e_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

tsc_mass_spec_results$label <- tsc_mass_spec_results$gene
tsc_mass_spec_results[!tsc_mass_spec_results$label %in% c("Ezh2", "Eed", "Rnf2", "Dnmt3b", "Suz12", "Dnmt3a", "Kdm2b", "Tet1", "Tet2", "Tet3", "Dnmt3a;Dnmt3b", "Ezh1;Ezh2", "Aebp2", "Jarid2"),"label"] <- NA

group <- ifelse(tsc_mass_spec_results$significance_corrected == "+" & abs(tsc_mass_spec_results$lfc) > 1, "significant", "not_significant")
group[tsc_mass_spec_results$label %in% c("Ezh2", "Eed", "Rnf2", "Dnmt3b", "Suz12", "Dnmt3a", "Kdm2b", "Tet1", "Tet2", "Tet3", "Dnmt3a;Dnmt3b", "Ezh1;Ezh2", "Aebp2", "Jarid2")] <- "selected"
tsc_mass_spec_results$group <- group

tsc_mass_spec_results <- tsc_mass_spec_results[rev(order(tsc_mass_spec_results$label)),]

pdf("Figure5e.pdf", height = 7, width = 7)
ggplot(data = tsc_mass_spec_results, aes(x = lfc, y = neg_log_pval, label = label, color = group)) + geom_point(size = 1) + theme_classic() + ylab("-Log(P-value)") + xlab("Log2 fold change") + geom_text_repel(size = 3, segment.size = 0.3, box.padding = 0.5) + scale_color_manual(values = c("grey", "firebrick", "dodgerblue", "darkblue")) + geom_vline(xintercept = c(-1, 1), lty = 2) + geom_hline(yintercept = min(subset(tsc_mass_spec_results, group == "significant")$neg_log_pval), lty = 2)
dev.off()

###################################################################################################################
## Extended Data Figure 7a
###################################################################################################################

avg_features <- read.table(file = "Extended_Data_Figure7a_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_features_df <- melt(avg_features)
avg_features_df$variable <- factor(avg_features_df$variable, levels = c("Epi_WT", "Epi_KDM2BKO", "Epi_RNF2KO", "Epi_EEDKO", "ExE_WT", "ExE_KDM2BKO", "ExE_RNF2KO", "ExE_EEDKO",                 "TSC1_WT", "TSC1_KDM2BKO", "TSC1_RNF2KO", "TSC3_EEDKO"))

## Overview plots
pdf("Extended_Data_Figure7a.pdf", height = 7, width = 12)
par(mar=c(18, 4, 4, 2))
vioplot(value ~ variable, data = subset(avg_features_df, feature == "HMD"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), main = "HMD/PMD", las = 2, plotCentre = "line", side = "left", col = brewer.pal(4, "Paired")[1])
vioplot(value ~ variable, data = subset(avg_features_df, feature == "PMD"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), add = TRUE, plotCentre = "line", side = "right", col = brewer.pal(4, "Paired")[2])
vioplot(value ~ variable, data = subset(avg_features_df, feature == "CGI"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), main = "CGI/ExE hyper CGI", las = 2, plotCentre = "line", side = "left", col = brewer.pal(4, "Paired")[4])
vioplot(value ~ variable, data = subset(avg_features_df, feature == "exe_hyper_cgi"), xlab = "", ylab = "Mean methylation", ylim = c(0,1), add = TRUE, plotCentre = "line", side = "right", col = brewer.pal(4, "Paired")[3])
dev.off()

###################################################################################################################
## Extended Data Figure 7b
###################################################################################################################

## Extended Data Figure 7b was generated using the processed methylation rates of Epiblast and ExE WT,
## KDM2B KO, RNF2 KO and EED KO as well as TSC1 WT, TSC1 KDM2B KO, TSC1 RNF2 KO and TSC3 EED KO
## visualized in IGV.

###################################################################################################################
## Extended Data Figure 7d
###################################################################################################################

methyl_data <- data.frame(fread("Extended_Data_Figure7d_source_data.tsv"), stringsAsFactors = FALSE)

pdf("Extended_Data_Figure7d.pdf", width = 10, height = 5.5)
layout(matrix(1:2, ncol = 2, byrow = TRUE))
smoothScatter(methyl_data[,"TSC1_WT"], methyl_data[,"TSC1_WT_DNMT1i_7d"], xlab = "TSC WT", ylab = "TSC DNMT1i 7d", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
smoothScatter(methyl_data[,"TSC1_WT"], methyl_data[,"TSC1_WT_DNMT1i_7d_recovery_14d"], xlab = "TSC WT", ylab = "TSC DNMT1i 7d recovery 14d", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0, xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lty = 2)
abline(-0.2, 1, lty = 2, col = "grey30")
abline(0.2, 1, lty = 2, col = "grey30")
dev.off()

###################################################################################################################
## Extended Data Figure 7e
###################################################################################################################

avg_features_dnmt1i_wgbs <- read.table(file = "Extended_Data_Figure7e_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_features_dnmt1i_wgbs_df <- melt(avg_features_dnmt1i_wgbs)
avg_features_dnmt1i_wgbs_df$variable <- factor(avg_features_dnmt1i_wgbs_df$variable, levels = rev(levels(avg_features_dnmt1i_wgbs_df$variable)))

avg_features_dnmt1i_wgbs_methylated <- subset(avg_features_dnmt1i_wgbs, TSC1_WT > 0.2)
avg_features_dnmt1i_wgbs_methylated_df <- melt(avg_features_dnmt1i_wgbs_methylated)
avg_features_dnmt1i_wgbs_methylated_df$variable <- factor(avg_features_dnmt1i_wgbs_methylated_df$variable, levels = rev(levels(avg_features_dnmt1i_wgbs_methylated_df$variable)))

ratio_many_features_dnmt1i_wgbs_methylated <- data.frame(
TSC1_WT_DNMT1i_7d = avg_features_dnmt1i_wgbs_methylated$TSC1_WT_DNMT1i_7d / avg_features_dnmt1i_wgbs_methylated$TSC1_WT,
TSC1_WT_DNMT1i_7d_recovery_14d = avg_features_dnmt1i_wgbs_methylated$TSC1_WT_DNMT1i_7d_recovery_14d / avg_features_dnmt1i_wgbs_methylated$TSC1_WT,
feature = avg_features_dnmt1i_wgbs_methylated$feature)

ratio_many_features_dnmt1i_wgbs_methylated_df <- melt(ratio_many_features_dnmt1i_wgbs_methylated)
ratio_many_features_dnmt1i_wgbs_methylated_df$variable <- factor(ratio_many_features_dnmt1i_wgbs_methylated_df$variable, levels = c("TSC1_WT_DNMT1i_7d_recovery_14d", "TSC1_WT_DNMT1i_7d"))
ratio_many_features_dnmt1i_wgbs_methylated_df$feature <- factor(ratio_many_features_dnmt1i_wgbs_methylated_df$feature, levels = c("CGIshelf", "CGIshore", "exe_hyper_cgi", "CGI", "genebody", "promoter", "SINE", "LTR", "LINE", "PMD", "HMD"))

pdf("Extended_Data_Figure7e.pdf", height = 7, width = 10)
ggplot(ratio_many_features_dnmt1i_wgbs_methylated_df, aes(x=variable, y = value, fill = feature)) + geom_boxplot(outlier.shape = NA) + coord_flip(ylim = c(0, 1.5)) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + xlab("") + ylab("Ratio methylation DNMT1i vs DMSO") + scale_fill_manual(values = c("turquoise3", "turquoise1", brewer.pal(4, "Paired")[c(4,3)], "purple", "magenta", rev(brewer.pal(3, "Oranges")), brewer.pal(4, "Paired")[c(2,1)])) + geom_hline(yintercept = 1, lty = 2, col = "grey") + geom_hline(yintercept = 0.5, lty = 2, col = "grey") + geom_hline(yintercept = 0.25, lty = 2, col = "grey")
dev.off()

###################################################################################################################
## Extended Data Figure 7f
###################################################################################################################

avg_features_dnmt1i_key_time_points <- read.table(file = "Extended_Data_Figure7f_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

avg_features_dnmt1i_key_time_points_df <- melt(avg_features_dnmt1i_key_time_points)
avg_features_dnmt1i_key_time_points_df$sample <- factor(sapply(as.character(avg_features_dnmt1i_key_time_points_df$variable), function(x) strsplit(x, "_DNMT1i|_DMSO")[[1]][1]), levels = c("TSC1_WT", "TSC1_WT_early_passage", "TSC1_WT_late_passage", "TSC2_WT"))
avg_features_dnmt1i_key_time_points_df$treatment <- factor(gsub("TSC1_WT|TSC2_WT|early_passage|late_passage", "untreated", sapply(as.character(avg_features_dnmt1i_key_time_points_df$variable), function(x) rev(strsplit(x, "passage_|WT_")[[1]])[1])), levels = rev(c("untreated", "DMSO_7d", "DNMT1i_7d", "DNMT1i_7d_recovery_14d")))

pdf("Extended_Data_Figure7f.pdf", height = 8, width = 14)
ggplot(avg_features_dnmt1i_key_time_points_df, aes(y=treatment, x=value, fill = feature)) + geom_density_ridges(aes(fill = feature), rel_min_height=.01, scale=1, alpha = 0.5, panel_scaling = FALSE) + coord_cartesian(xlim = c(0, 1)) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + ylab("") + xlab("Mean methylation") + scale_fill_manual(values = brewer.pal(3, "Paired")) + facet_grid(~sample)
dev.off()

###################################################################################################################
## Extended Data Figure 7g
###################################################################################################################

ratio_features_dnmt1i_methylated_key_time_points <- read.table(file = "Extended_Data_Figure7g_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

ratio_features_dnmt1i_methylated_key_time_points_df <- melt(ratio_features_dnmt1i_methylated_key_time_points)
ratio_features_dnmt1i_methylated_key_time_points_df$sample <- factor(ratio_features_dnmt1i_methylated_key_time_points_df$sample, levels = rev(c("TSC1_WT", "TSC1_WT_early_passage", "TSC1_WT_late_passage", "TSC2_WT")))
ratio_features_dnmt1i_methylated_key_time_points_df$variable <- factor(ratio_features_dnmt1i_methylated_key_time_points_df$variable, levels = rev(c("DNMT1i_7d", "DNMT1i_7d_recovery_14d")))

pdf("Extended_Data_Figure7g.pdf", height = 5, width = 12)
ggplot(ratio_features_dnmt1i_methylated_key_time_points_df, aes(x=sample, y = value, fill = variable)) + geom_boxplot(outlier.shape = NA) + coord_flip(ylim = c(0, 1.5)) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + xlab("") + ylab("Ratio methylation DNMT1i vs DMSO") + scale_fill_manual(values = c(brewer.pal(10, "Paired")[c(9,6)])) + facet_grid(~feature) + geom_hline(yintercept = c(0.5, 1), lty = 2, col = "grey")
dev.off()

###################################################################################################################
## Extended Data Figure 8b
###################################################################################################################

figure_8b_source_data <- read.table(file = "Extended_Data_Figure8b_source_data.tsv", sep = "\\t", header = TRUE, stringsAsFactors = FALSE)

samples <- unique(figure_8b_source_data$sample)

norm_mat_list <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(figure_8b_source_data, sample == x)[,!colnames(figure_8b_source_data) %in% "sample"]), k_upstream=100, k_downstream=100, k_target=67, extend = c(5000, 5000), signal_name = x, target_name = "prc_hyper_cgi"))
names(norm_mat_list) <- samples

col_fun <- list(
    "H3K27me3_TSC1_WT_EZH2i_5w_recovery_4w" = colorRamp2(c(0, 10), c("white", "violetred2")),
    "H3K27me3_TSC1_WT" = colorRamp2(c(0, 10), c("white", "violetred2")),
    "H2AK119ub1_TSC1_WT_EZH2i_5w_recovery_4w" = colorRamp2(c(0, 10), c("white", "deepskyblue3")),
    "H2AK119ub1_TSC1_WT" = colorRamp2(c(0, 10), c("white", "deepskyblue3")),
    "H3K4me3_TSC1_WT_EZH2i_5w_recovery_4w" = colorRamp2(c(0, 100), c("white", "seagreen3")),
    "H3K4me3_TSC1_WT" = colorRamp2(c(0, 100), c("white", "seagreen3"))
    )

meta_height <- c(
    "H3K27me3_TSC1_WT_EZH2i_5w_recovery_4w" = 5,
    "H3K27me3_TSC1_WT" = 5,
    "H2AK119ub1_TSC1_WT_EZH2i_5w_recovery_4w" = 5,
    "H2AK119ub1_TSC1_WT" = 5,
    "H3K4me3_TSC1_WT_EZH2i_5w_recovery_4w" = 50,
    "H3K4me3_TSC1_WT" = 50
    )

ht_list <- NULL
for(sample in samples)
{
    ht_list <- ht_list + EnrichedHeatmap(norm_mat_list[[sample]], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col="black", lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6))
}

pdf("Extended_Data_Figure8b.pdf", width = 20)
draw(ht_list)
dev.off()

###################################################################################################################
## Extended Data Figure 8d
###################################################################################################################

## Extended Data Figure 8d was generated using the smoothed EED ChIPseq tracks of TSC1 WT untreated and
## treated with EZH2i visualized in IGV.

###################################################################################################################
## Extended Data Figure 8e
###################################################################################################################

figure_8e_source_data <- read.table(file = "Extended_Data_Figure8e_source_data.tsv", sep = "\\t", header = TRUE, stringsAsFactors = FALSE)

samples <- unique(figure_8e_source_data$sample)

norm_mat_list <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(figure_8e_source_data, sample == x)[,!colnames(figure_8e_source_data) %in% c("sample", "row_split", "cgi")]), k_upstream=100, k_downstream=100, extend = c(5000, 5000), signal_name = x, target_name = "eed_peaks"))
names(norm_mat_list) <- samples

col_fun <- list(
    "TSC1_WT_EED_ab240650_merged" = colorRamp2(c(0, 15), c("white", "mediumpurple1")),
    "TSC1_WT_EZH2i_5w4d_EED_ab240650_merged" = colorRamp2(c(0, 15), c("white", "mediumpurple1")),
    "WGBS_TSC1_WT" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")),
    "WGBS_TSC3_EEDKO" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    )

meta_height <- c(
    "TSC1_WT_EED_ab240650_merged" = 15,
    "TSC1_WT_EZH2i_5w4d_EED_ab240650_merged" = 15,
    "WGBS_TSC1_WT" = 1,
    "WGBS_TSC3_EEDKO" = 1
    )

ht_list <- NULL
for(sample in samples)
{
    ht_list <- ht_list + EnrichedHeatmap(norm_mat_list[[sample]], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, row_split = subset(figure_8e_source_data, sample == samples[1])$row_split, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col=c("black", "grey"), lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 4))
}
ht_list <- ht_list + Heatmap(as.matrix(subset(figure_8e_source_data, sample == samples[1])[,"cgi",drop=FALSE]), name = "CGI", show_row_names = FALSE, width = unit(3, "mm"), col = c("exe_hyper_cgi" = "palegreen3", "other" = "grey95"))

pdf("Extended_Data_Figure8e.pdf", width = 20)
draw(ht_list)
dev.off()

###################################################################################################################
## Extended Data Figure 8f
###################################################################################################################

avg_eed_peak_tsc_wt <- read.table(file = "Extended_Data_Figure8f_source_data.tsv", sep = "\\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

avg_eed_peak_tsc_wt_df <- melt(avg_eed_peak_tsc_wt)
avg_eed_peak_tsc_wt_df$type <- ifelse(grepl("DMSO", as.character(avg_eed_peak_tsc_wt_df$variable)), "control", "treatment")
avg_eed_peak_tsc_wt_df$peak_group <- factor(avg_eed_peak_tsc_wt_df$peak_group, levels = c("overlap_peak_ezh2i", "no_overlap_peak_ezh2i"))

pdf("Extended_Data_Figure8f.pdf", height = 8, width = 20)
par(mar=c(18, 4, 4, 2))
vioplot(value ~ variable, data = subset(avg_eed_peak_tsc_wt_df, peak_group == "overlap_peak_ezh2i"), xlab = "", ylab = "Mean methylation", main = "TSC WT peak overlap EZH1i peak / no overlap EZH2i peak", las = 2, side = "left", plotCentre = "line", col = "grey40")
vioplot(value ~ variable, data = subset(avg_eed_peak_tsc_wt_df, peak_group == "no_overlap_peak_ezh2i"), xlab = "", ylab = "Mean methylation", main = "overlap_peak_ezh2i", las = 2, add = TRUE, side = "right", plotCentre = "line", col = "lightgrey")
dev.off()

###################################################################################################################
## Extended Data Figure 8h
###################################################################################################################

esc_mass_spec_results <- read.table(file = "Extended_Data_Figure8h_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)

significant_interactions_esc_tsc <- list(ESC = as.character(unlist(sapply(as.character(subset(esc_mass_spec_results, significance_corrected == "+" & lfc > 1)$gene), function(x) strsplit(x, ";")[[1]]))), TSC = as.character(unlist(sapply(as.character(subset(tsc_mass_spec_results, significance_corrected == "+" & lfc > 1)$gene), function(x) strsplit(x, ";")[[1]]))))

pdf("Extended_Data_Figure8h.pdf")
plot(Venn(significant_interactions_esc_tsc), doWeights = TRUE)
dev.off()

###################################################################################################################
## Extended Data Figure 8i
###################################################################################################################

## GO term enrichment TSC genes
gene_sets <- list(ESC_only = sort(setdiff(significant_interactions_esc_tsc$ESC, significant_interactions_esc_tsc$TSC)), TSC_only = sort(setdiff(significant_interactions_esc_tsc$TSC, significant_interactions_esc_tsc$ESC)), TSC_ESC_both = sort(intersect(significant_interactions_esc_tsc$TSC, significant_interactions_esc_tsc$ESC)))

ora <- lapply(gene_sets, function(x) WebGestaltR(interestGene = x, enrichMethod = "ORA", enrichDatabase = "geneontology_Molecular_Function", organism = "mmusculus", referenceSet = "genome", isOutput = FALSE, interestGeneType="genesymbol", referenceGeneType="genesymbol", minNum = 10, maxNum = 500, sigMethod = "top", topThr = 20))
ora_detail <- lapply(gene_sets, function(x) WebGestaltR(interestGene = x, enrichMethod = "ORA", enrichDatabase = "geneontology_Molecular_Function", organism = "mmusculus", referenceSet = "genome", isOutput = FALSE, interestGeneType="genesymbol", referenceGeneType="genesymbol", minNum = 10, maxNum = 500, sigMethod = "fdr"))

top_gene_sets <- lapply(ora, function(x) subset(x, FDR < 0.05)$geneSet)
names(top_gene_sets) <- names(ora_detail)

top_gene_sets_combined <- unique(do.call(c, top_gene_sets))

ora_detail[["TSC_only"]]$sample <- "TSC_only"
ora_detail[["ESC_only"]]$sample <- "ESC_only"
ora_detail[["TSC_ESC_both"]]$sample <- "TSC_ESC_both"

ora_detail_combined <- do.call(rbind, lapply(ora_detail, function(x) x[,c("geneSet", "description", "FDR", "overlap", "size", "sample")]))
ora_detail_combined <- subset(ora_detail_combined, geneSet %in% top_gene_sets_combined)

top_gene_sets_combined_description <- data.frame(unique(ora_detail_combined[,c("geneSet", "description")]), row.names = 1)
top_gene_sets_combined_description <- top_gene_sets_combined_description[top_gene_sets_combined,,drop=FALSE]

ora_detail_combined$sample <- factor(ora_detail_combined$sample, levels = c("ESC_only", "TSC_only", "TSC_ESC_both"))
ora_detail_combined$description <- factor(ora_detail_combined$description, unique(top_gene_sets_combined_description$description))
ora_detail_combined$gene_ratio <- ora_detail_combined$overlap / ora_detail_combined$size

ora_detail_combined_mat <- data.frame(dcast(description ~ sample, data = ora_detail_combined[,c("sample", "description", "FDR")]), row.names = 1)

breaksList <- seq(0, 0.05, 0.001)

pdf("Extended_Data_Figure8i.pdf", width = 10, height = 7)
pheatmap::pheatmap(t(ora_detail_combined_mat), cluster_rows = FALSE, cluster_cols = FALSE, breaks = breaksList, color = rev(colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length(breaksList))), na_col = "white")
dev.off()
