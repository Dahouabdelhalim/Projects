###################################################################################################################
## This script contains the code to produce Figure 6 and Extended Data Figures 9 and 10
## Weigert, Hetzel et al. Dynamic antagonism between key repressive pathways maintains the placental epigenome 2023
## Author: Sara Hetzel
###################################################################################################################

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(WebGestaltR)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(DESeq2)
library(circlize)
library(Vennerable)
library(plyr)

###################################################################################################################
## Data preparation
###################################################################################################################

## Read source data tables
count_data <- read.table(file = "Figure6_Extended_Data_Figure9_source_data_counts.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)
tpm <- read.table(file = "Figure6_Extended_Data_Figure9_source_data_tpm.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)
tpm_with_noncoding <- read.table(file = "Extended_Data_Figure10_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, row.names = 1, header = TRUE)

annotation_genes <- tpm[,1:5]
annotation_genes_with_noncoding <- tpm_with_noncoding[,1:5]

tpm <- tpm[,-c(1:5)]
tpm_with_noncoding <- tpm_with_noncoding[,-c(1:5)]

## Remove counts and TPMs to autosomes for protein coding genes
annotation_genes <- subset(annotation_genes, !chr %in% c("chrX", "chrY", "chrM"))
count_data <- count_data[rownames(annotation_genes),]
tpm <- tpm[rownames(annotation_genes),]

## Metadata per sample
col_data <- data.frame(row.names = colnames(count_data), type = sapply(colnames(count_data), function(x) strsplit(x, "TSC[1-9]_|_Rep")[[1]][2]), replicate = sapply(colnames(count_data), function(x) rev(strsplit(x, "_")[[1]])[1]))
col_data$type <- factor(col_data$type, levels = c("WT", "3BKO", "RNF2KO", "KDM2BKO", "EEDKO", "WT_DMSO_7d_acute_response", "WT_DNMT1i_7d_acute_response", "WT_EZH2i_7d_acute_response", "WT_DNMT1i_EZH2i_7d_acute_response", "WT_DMSO_7d", "WT_DNMT1i_2d", "WT_DNMT1i_7d", "WT_DNMT1i_7d_recovery_14d", "WT_DMSO_5w", "WT_EZH2i_5w", "WT_EZH2i_5w_recovery_1w", "WT_EZH2i_5w_recovery_4w"))

## Re-order columns data, generate log2-transformed TPM
count_data <- count_data[,rownames(col_data)]
tpm <- tpm[,rownames(col_data)]
tpm_with_noncoding <- tpm_with_noncoding[,rownames(col_data)]

tpm_log2 <- log2(tpm + 1)
tpm_with_noncoding_log2 <- log2(tpm_with_noncoding + 1)

## Average replicates, generate z-score TPMs
cal_z_score <- function(x)
{
  (x - mean(x)) / sd(x)
}

tpm_avg <- sapply(as.character(unique(col_data$type)), function(x) apply(tpm[,rownames(subset(col_data, type == x))], 1, mean))
tpm_avg_log2 <- log2(tpm_avg + 1)
tpm_avg_log2_zscore <- t(apply(tpm_avg_log2, 1, cal_z_score))

tpm_avg_with_noncoding <- sapply(as.character(unique(col_data$type)), function(x) apply(tpm_with_noncoding[,rownames(subset(col_data, type == x))], 1, mean))
tpm_avg_with_noncoding_log2 <- log2(tpm_avg_with_noncoding + 1)

## Colors per type
colors_type <- c("#CCAA60", "#A94098", "#748DAA", "#6082BF", "grey40", rep("darkgrey", 3), brewer.pal(4, "Purples")[4:2], brewer.pal(4, "Purples")[3], brewer.pal(5, "Blues")[5:2], brewer.pal(4, "Greens")[4])
names(colors_type) <- c("3BKO", "EEDKO", "RNF2KO", "KDM2BKO", "WT", "WT_DMSO_7d", "WT_DMSO_7d_acute_response", "WT_DMSO_5w", "WT_DNMT1i_2d", "WT_DNMT1i_7d", "WT_DNMT1i_7d_recovery_14d","WT_DNMT1i_7d_acute_response", "WT_EZH2i_7d_acute_response", "WT_EZH2i_5w", "WT_EZH2i_5w_recovery_1w", "WT_EZH2i_5w_recovery_4w", "WT_DNMT1i_EZH2i_7d_acute_response")

###################################################################################################################
## Differentially expressed genes
###################################################################################################################

col_data_deg <- col_data
col_data_deg$type[grepl("DMSO", rownames(col_data_deg))] <- "WT"

inhibitor_acute_samples <- rownames(subset(col_data_deg, grepl("acute_response", type) | type == "WT"))

dds_inhibitor_acute <- DESeqDataSetFromMatrix(countData = count_data[,inhibitor_acute_samples], colData = droplevels(col_data_deg[inhibitor_acute_samples,]), design = ~type)
keep <- rowSums(counts(dds_inhibitor_acute)) >= 10
print(table(keep))
dds_inhibitor_acute <- dds_inhibitor_acute[keep,]
dds_inhibitor_acute <- DESeq(dds_inhibitor_acute)

active_genes_inhibitor_acute <- names(which(apply(tpm[,inhibitor_acute_samples], 1, function(x) length(which(x > 2))) >= 2))

results_combined_inhibitor_acute <- list(results(dds_inhibitor_acute, contrast = c("type", "WT_DNMT1i_7d_acute_response", "WT")),
                                         results(dds_inhibitor_acute, contrast = c("type", "WT_EZH2i_7d_acute_response", "WT")),
                                         results(dds_inhibitor_acute, contrast = c("type", "WT_DNMT1i_EZH2i_7d_acute_response", "WT")))

results_combined_inhibitor_acute_significant <- lapply(results_combined_inhibitor_acute, function(x) subset(x, (padj < 0.05) & (abs(log2FoldChange) > 2)))
results_combined_inhibitor_acute_significant_filtered <- lapply(results_combined_inhibitor_acute_significant, function(x) x[intersect(rownames(x), active_genes_inhibitor_acute),])
names(results_combined_inhibitor_acute_significant_filtered) <- c("WT_DNMT1i_7d_acute_response", "WT_EZH2i_7d_acute_response", "WT_DNMT1i_EZH2i_7d_acute_response")

unlist(lapply(results_combined_inhibitor_acute_significant_filtered, nrow))
# WT_DNMT1i_7d_acute_response        WT_EZH2i_7d_acute_response
#                         238                               377
# WT_DNMT1i_EZH2i_7d_acute_response
#                         851

results_combined_inhibitor_acute_significant_filtered_up <- lapply(results_combined_inhibitor_acute_significant_filtered, function(x) subset(x, log2FoldChange > 0))
results_combined_inhibitor_acute_significant_filtered_down <- lapply(results_combined_inhibitor_acute_significant_filtered, function(x) subset(x, log2FoldChange < 0))

###################################################################################################################
## Overrepresentation analysis
###################################################################################################################

gene_names_up <- lapply(results_combined_inhibitor_acute_significant_filtered_up, function(x) sort(annotation_genes[rownames(x),"gene_name"]))
gene_names_down <- lapply(results_combined_inhibitor_acute_significant_filtered_down, function(x) sort(annotation_genes[rownames(x),"gene_name"]))

names(gene_names_up) <- paste(names(gene_names_up), "up", sep = "_")
names(gene_names_down) <- paste(names(gene_names_down), "down", sep = "_")

gene_sets <- c(gene_names_up, gene_names_down)

ora <- lapply(gene_sets, function(x) WebGestaltR(interestGene = x, enrichMethod = "ORA", enrichDatabase = "geneontology_Biological_Process", organism = "mmusculus", referenceSet = "genome", isOutput = FALSE, interestGeneType="genesymbol", referenceGeneType="genesymbol", minNum = 10, maxNum = 500, sigMethod = "top", topThr = 10))

ora_detail <- lapply(gene_sets, function(x) WebGestaltR(interestGene = x, enrichMethod = "ORA", enrichDatabase = "geneontology_Biological_Process", organism = "mmusculus", referenceSet = "genome", isOutput = FALSE, interestGeneType="genesymbol", referenceGeneType="genesymbol", minNum = 10, maxNum = 500, sigMethod = "fdr"))

top_gene_sets <- lapply(ora, function(x) x$geneSet)
names(top_gene_sets) <- names(ora_detail)

top_10_up <- unique(do.call(c, top_gene_sets[c("WT_DNMT1i_7d_acute_response_up", "WT_EZH2i_7d_acute_response_up", "WT_DNMT1i_EZH2i_7d_acute_response_up")]))
top_10_down <- unique(do.call(c, top_gene_sets[c("WT_DNMT1i_7d_acute_response_down", "WT_EZH2i_7d_acute_response_down", "WT_DNMT1i_EZH2i_7d_acute_response_down")]))

ora_detail[["WT_DNMT1i_7d_acute_response_up"]]$sample <- "WT_DNMT1i_7d_acute_response_up"
ora_detail[["WT_DNMT1i_7d_acute_response_down"]]$sample <- "WT_DNMT1i_7d_acute_response_down"
ora_detail[["WT_EZH2i_7d_acute_response_up"]]$sample <- "WT_EZH2i_7d_acute_response_up"
ora_detail[["WT_EZH2i_7d_acute_response_down"]]$sample <- "WT_EZH2i_7d_acute_response_down"
ora_detail[["WT_DNMT1i_EZH2i_7d_acute_response_up"]]$sample <- "WT_DNMT1i_EZH2i_7d_acute_response_up"
ora_detail[["WT_DNMT1i_EZH2i_7d_acute_response_down"]]$sample <- "WT_DNMT1i_EZH2i_7d_acute_response_down"

###################################################################################################################
## Figure 6a
###################################################################################################################

source_data_Fig6a <- read.table(file = "Figure6a_source_data.tsv", sep = "\\t", stringsAsFactors = FALSE, header = TRUE)
source_data_Fig6a$day <- as.character(source_data_Fig6a$day)
source_data_Fig6a_df <- melt(source_data_Fig6a)
source_data_Fig6a_summary_df <- ddply(source_data_Fig6a_df, .(day, variable), summarize, avg = mean(value), sd = sd(value))
source_data_Fig6a_df$day <- as.numeric(source_data_Fig6a_df$day)
source_data_Fig6a_summary_df$day <- as.numeric(source_data_Fig6a_summary_df$day)

set.seed(42)
pdf("Figure6a.pdf", height = 5, width = 8)
ggplot(source_data_Fig6a_summary_df, aes(x=day, y=avg, color = variable, group = variable)) + geom_point(shape = 17, size = 2) + geom_line() + theme_classic() + ylab("Cell count") + geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2) + geom_point(data = source_data_Fig6a_df, aes(y = value), position = position_jitter(width = 0.2))
dev.off()

###################################################################################################################
## Figure 6b
###################################################################################################################

pdf("Figure6b.pdf")
plot(Venn(lapply(results_combined_inhibitor_acute_significant_filtered_up, rownames)), doWeights = TRUE) ## Up-regulated
plot(Venn(lapply(results_combined_inhibitor_acute_significant_filtered_down, rownames)), doWeights = TRUE) ## Down-regulated
dev.off()

###################################################################################################################
## Figure 6c
###################################################################################################################

deg_down_top_go_terms_double_treatment <- unique(as.character(sapply(subset(ora[["WT_DNMT1i_EZH2i_7d_acute_response_down"]], description == "regulation of mitotic cell cycle")$userId, function(x) strsplit(x, ";")[[1]])))

deg_down_top_go_terms_double_treatment_id <- sapply(deg_down_top_go_terms_double_treatment, function(x) rownames(subset(annotation_genes_with_noncoding, gene_name == x)))
order_deg_down_top_go_terms_double_treatment <- order(tpm_avg_with_noncoding_log2[deg_down_top_go_terms_double_treatment_id,"WT_DMSO_7d_acute_response"], decreasing = TRUE)
deg_down_top_go_terms_double_treatment <- deg_down_top_go_terms_double_treatment[order_deg_down_top_go_terms_double_treatment]
deg_down_top_go_terms_double_treatment_id <- deg_down_top_go_terms_double_treatment_id[order_deg_down_top_go_terms_double_treatment]

breaksList <- seq(0, 10, 0.1)
pdf("Figure6c.pdf", width = 5, height = 12)
pheatmap::pheatmap(tpm_avg_with_noncoding_log2[deg_down_top_go_terms_double_treatment_id,c("WT_DMSO_7d_acute_response", "WT_DNMT1i_7d_acute_response", "WT_EZH2i_7d_acute_response", "WT_DNMT1i_EZH2i_7d_acute_response")], labels_row = deg_down_top_go_terms_double_treatment, cluster_cols = FALSE, cluster_rows = FALSE, breaks = breaksList, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)))
dev.off()

###################################################################################################################
## Figure 6d
###################################################################################################################

## Generate z-score matrices
tpm_avg_log2_zscore_wt_dnmt1i <- t(apply(tpm_avg_log2[,c("WT_DMSO_7d", "WT_DNMT1i_2d", "WT_DNMT1i_7d", "WT_DNMT1i_7d_recovery_14d")], 1, cal_z_score))
tpm_avg_log2_zscore_wt_ezh2i <- t(apply(tpm_avg_log2[,c("WT_DMSO_5w", "WT_EZH2i_5w", "WT_EZH2i_5w_recovery_1w", "WT_EZH2i_5w_recovery_4w")], 1, cal_z_score))

col_fun_rna <- colorRamp2(c(-2, 0, 2), c("#1B7837", "white", "#762A83"))

## Take DEG from acute response
m1 <- as.matrix(tpm_avg_log2_zscore_wt_dnmt1i[rownames(results_combined_inhibitor_acute_significant_filtered_up[["WT_DNMT1i_7d_acute_response"]]),])
m2 <- as.matrix(tpm_avg_log2_zscore_wt_dnmt1i[rownames(results_combined_inhibitor_acute_significant_filtered_down[["WT_DNMT1i_7d_acute_response"]]),])

m3 <- as.matrix(tpm_avg_log2_zscore_wt_ezh2i[rownames(results_combined_inhibitor_acute_significant_filtered_up[["WT_EZH2i_7d_acute_response"]]),])
m4 <- as.matrix(tpm_avg_log2_zscore_wt_ezh2i[rownames(results_combined_inhibitor_acute_significant_filtered_down[["WT_EZH2i_7d_acute_response"]]),])

pdf("Figure6d.pdf", width = 4.5, height = 8)
Heatmap(as.matrix(tpm_avg_log2_zscore_wt_dnmt1i[rownames(results_combined_inhibitor_acute_significant_filtered[["WT_DNMT1i_7d_acute_response"]]),]), show_row_names = FALSE, row_split = factor(ifelse(results_combined_inhibitor_acute_significant_filtered[["WT_DNMT1i_7d_acute_response"]]$log2FoldChange > 0, "up", "down"), levels = c("up", "down")), col = col_fun_rna, cluster_columns = FALSE, top_annotation = HeatmapAnnotation(up = anno_boxplot(m1, gp = gpar(fill = "red"), box_width = 0.4, ylim = c(-3, 3), outline = FALSE), down = anno_boxplot(m2, gp = gpar(fill = "blue"), box_width = 0.4, ylim = c(-3, 3), outline = FALSE)))

Heatmap(as.matrix(tpm_avg_log2_zscore_wt_ezh2i[rownames(results_combined_inhibitor_acute_significant_filtered[["WT_EZH2i_7d_acute_response"]]),]), show_row_names = FALSE, row_split = factor(ifelse(results_combined_inhibitor_acute_significant_filtered[["WT_EZH2i_7d_acute_response"]]$log2FoldChange > 0, "up", "down"), levels = c("up", "down")), col = col_fun_rna, cluster_columns = FALSE, top_annotation = HeatmapAnnotation(up = anno_boxplot(m3, gp = gpar(fill = "red"), box_width = 0.4, ylim = c(-3, 3), outline = FALSE), down = anno_boxplot(m4, gp = gpar(fill = "blue"), box_width = 0.4, ylim = c(-3, 3), outline = FALSE)))
dev.off()

###################################################################################################################
## Extended Data Figure 9a
###################################################################################################################

## Combine enriched terms for heatmap
deg_combined_up <- do.call(rbind, lapply(ora_detail[c("WT_DNMT1i_7d_acute_response_up", "WT_EZH2i_7d_acute_response_up", "WT_DNMT1i_EZH2i_7d_acute_response_up")], function(x) x[,c("geneSet", "description", "FDR", "overlap", "size", "sample")]))
deg_combined_down <- do.call(rbind, lapply(ora_detail[c("WT_DNMT1i_7d_acute_response_down", "WT_EZH2i_7d_acute_response_down", "WT_DNMT1i_EZH2i_7d_acute_response_down")], function(x) x[,c("geneSet", "description", "FDR", "overlap", "size", "sample")]))

deg_combined_up <- subset(deg_combined_up, geneSet %in% top_10_up)
deg_combined_down <- subset(deg_combined_down, geneSet %in% top_10_down)

top_10_up_description <- data.frame(unique(deg_combined_up[,c("geneSet", "description")]), row.names = 1)
top_10_up_description <- top_10_up_description[top_10_up,,drop=FALSE]

top_10_down_description <- data.frame(unique(deg_combined_down[,c("geneSet", "description")]), row.names = 1)
top_10_down_description <- top_10_down_description[top_10_down,,drop=FALSE]

deg_combined_up$sample <- factor(deg_combined_up$sample, levels = c("WT_DNMT1i_7d_acute_response_up", "WT_EZH2i_7d_acute_response_up", "WT_DNMT1i_EZH2i_7d_acute_response_up"))
deg_combined_up$gene_ratio <- deg_combined_up$overlap / deg_combined_up$size
deg_combined_up$description <- factor(deg_combined_up$description, unique(top_10_up_description$description))

deg_combined_down$sample <- factor(deg_combined_down$sample, levels = c("WT_DNMT1i_7d_acute_response_down", "WT_EZH2i_7d_acute_response_down", "WT_DNMT1i_EZH2i_7d_acute_response_down"))
deg_combined_down$gene_ratio <- deg_combined_down$overlap / deg_combined_down$size
deg_combined_down$description <- factor(deg_combined_down$description, unique(top_10_down_description$description))

deg_combined_up_mat <- data.frame(dcast(description ~ sample, data = deg_combined_up[,c("sample", "description", "FDR")]), row.names = 1)
deg_combined_down_mat <- data.frame(dcast(description ~ sample, data = deg_combined_down[,c("sample", "description", "FDR")]), row.names = 1)

breaksList <- seq(0, 0.05, 0.001)

pdf("Extended_Data_Figure9a.pdf", width = 7, height = 10)
pheatmap::pheatmap(deg_combined_up_mat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = breaksList, color = rev(colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length(breaksList))), na_col = "white", main = "Up-regulated genes")
pheatmap::pheatmap(deg_combined_down_mat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = breaksList, color = rev(colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length(breaksList))), na_col = "white", main = "Down-regulated genes")
dev.off()

###################################################################################################################
## Extended Data Figure 9b
###################################################################################################################

active_genes_all_avg <- names(which(apply(tpm_avg, 1, function(x) length(which(x > 2))) >= 1))
breaksList <- seq(0.8, 1, by = 0.01)

pdf("Extended_Data_Figure9b.pdf", width = 12, height = 8)
pheatmap::pheatmap(cor(log2(tpm_avg[active_genes_all_avg,] + 1)), show_colnames = FALSE, breaks = breaksList, color = rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(length(breaksList))), clustering_distance_rows = as.dist(1 - cor(log2(tpm_avg[active_genes_all_avg,] + 1))), clustering_distance_cols = as.dist(1 - cor(log2(tpm_avg[active_genes_all_avg,] + 1))))
dev.off()

###################################################################################################################
## Extended Data Figure 9c
###################################################################################################################

order_samples <- c("WT_DMSO_7d_acute_response", "WT_DNMT1i_7d_acute_response", "WT_EZH2i_7d_acute_response", "WT_DNMT1i_EZH2i_7d_acute_response", "WT_DMSO_7d", "WT_DNMT1i_2d", "WT_DNMT1i_7d", "WT_DNMT1i_7d_recovery_14d", "WT_DMSO_5w", "WT_EZH2i_5w", "WT_EZH2i_5w_recovery_1w", "WT_EZH2i_5w_recovery_4w", "WT", "3BKO", "KDM2BKO", "RNF2KO", "EEDKO")

tpm_avg_log2_deg_groups_up_df <- melt(lapply(results_combined_inhibitor_acute_significant_filtered_up, function(x)
{
    data <- data.frame(tpm_avg_log2[rownames(x),order_samples], check.names = FALSE)
    data$deg <- "up"
    return(data)
}))

tpm_avg_log2_deg_groups_down_df <- melt(lapply(results_combined_inhibitor_acute_significant_filtered_down, function(x)
{
    data <- data.frame(tpm_avg_log2[rownames(x),order_samples], check.names = FALSE)
    data$deg <- "down"
    return(data)
}))

tpm_avg_log2_deg_groups_df <- rbind(tpm_avg_log2_deg_groups_up_df, tpm_avg_log2_deg_groups_down_df)
tpm_avg_log2_deg_groups_df$variable <- factor(tpm_avg_log2_deg_groups_df$variable, levels = order_samples)
tpm_avg_log2_deg_groups_df$deg <- factor(tpm_avg_log2_deg_groups_df$deg, levels = c("up", "down"))
tpm_avg_log2_deg_groups_df$L1 <- factor(tpm_avg_log2_deg_groups_df$L1, levels = c("WT_DNMT1i_7d_acute_response", "WT_EZH2i_7d_acute_response", "WT_DNMT1i_EZH2i_7d_acute_response"))

## Genes with ExE hyper CGIs overlapping its promoter
exe_hyper_promoter_genes <- c(
"4930447C04Rik", "6430628N08Rik", "A830005F24Rik", "Abcg8",
"Acan",          "Adra1d",        "Agbl2",         "Alx3",
"Alx4",          "Arhgef25",      "Arpp21",        "Arsi",
"Arsj",          "Asic2",         "Asic4",         "Astn1",
"Atp1a3",        "Avpr1a",        "BC035947",      "Barhl1",
"Bdnf",          "Brinp1",        "C030006K11Rik", "C1ql2",
"C1ql3",         "C1qtnf12",      "CK137956",      "Cacng6",
"Cadm3",         "Cartpt",        "Casp8",         "Catsperz",
"Ccdc114",       "Ccdc194",       "Cd200",         "Cd302",
"Cdh13",         "Cdh20",         "Cdh7",          "Cdh8",
"Cdo1",          "Cerkl",         "Chrm2",         "Cidea",
"Clcn4",         "Cldn11",        "Clvs2",         "Col12a1",
"Col26a1",       "Col2a1",        "Colca2",        "Corin",
"Cplx1",         "Cpxm2",         "Cpz",           "Creb3l1",
"Crhr1",         "Cryba2",        "Cspg4",         "Cul9",
"Cybrd1",        "Cyp24a1",       "Cyp7b1",        "Dab1",
"Dbx1",          "Dchs2",         "Dclk1",         "Dmrt1",
"Dmtn",          "Dnah10",        "Dpp6",          "Dpys",
"Dpysl3",        "Drd1",          "Ebf2",          "Ece1",
"Eef1a2",        "Efnb3",         "Egfem1",        "Egr4",
"Elavl2",        "Eln",           "Enpp2",         "Epdr1",
"Ephb1",         "Esrrg",         "Evc",           "Evc2",
"Evpl",          "Evx1",          "Evx2",          "Far2",
"Fbln5",         "Fbp1",          "Fev",           "Fgf17",
"Fndc1",         "Foxa2",         "Foxb1",         "Foxi2",
"Foxl1",         "Frem3",         "Fxyd7",         "Gad2",
"Galr1",         "Gdf1",          "Gfi1",          "Gfra1",
"Gfra2",         "Gjb6",          "Gldn",          "Gm13420",
"Gm19410",       "Gm5878",        "Gm6034",        "Gpc6",
"Gpr12",         "Gpr139",        "Gpr158",        "Gpr26",
"Gpr6",          "Gpr83",         "Gramd1b",       "Grid2ip",
"Grik1",         "Grm7",          "Gsc",           "Gsx2",
"H2-Q10",        "H2-Q6",         "Hck",           "Hcrtr2",
"Helt",          "Hhip",          "Hmcn2",         "Hoxa11",
"Hoxb1",         "Hoxb5",         "Hoxb7",         "Hoxc11",
"Hoxc12",        "Hoxc5",         "Hoxc8",         "Hoxd11",
"Hoxd12",        "Hrh2",          "Hspb9",         "Htr1f",
"Igfbp6",        "Igsf21",        "Il12rb2",       "Irx6",
"Itga4",         "Kat2a",         "Kcna6",         "Kcnc2",
"Kcnh5",         "Kcnh8",         "Kcnip2",        "Kcnj6",
"Kcnk9",         "Kcnv1",         "Kif1a",         "Klhdc9",
"Kremen2",       "Krt87",         "L3mbtl4",       "Lbx1",
"Leng9",         "Lepr",          "Lfng",          "Lgr5",
"Lhfpl4",        "Lhx5",          "Lin7a",         "Lrfn2",
"Lrmda",         "Lrrtm1",        "Lrtm2",         "Lvrn",
"Lyl1",          "Mab21l1",       "Magi2",         "Marc1",
"March10",       "March4",        "Mcoln2",        "Me3",
"Megf11",        "Mfap2",         "Mmp2",          "Myo16",
"Myocd",         "Myod1",         "Naprt",         "Nell1",
"Neurod1",       "Neurod2",       "Neurog3",       "Ngb",
"Ngf",           "Nkx1-1",        "Nkx1-2",        "Nkx2-1",
"Nkx2-5",        "Nkx2-9",        "Nlrc5",         "Nov",
"Npb",           "Npl",           "Npy5r",         "Nr2e1",
"Nrn1",          "Ntng2",         "Ntsr1",         "Nxph1",
"Odf3l1",        "Olig3",         "Ostm1",         "Otop3",
"Otp",           "Otx1",          "Otx2",          "Oxct2a",
"Oxct2b",        "Oxt",           "Pak7",          "Pax3",
"Pax5",          "Pax6",          "Pax7",          "Pax9",
"Pcdh7",         "Pcdh8",         "Pcdhac1",       "Pcsk9",
"Pde4b",         "Pdgfra",        "Pdzrn4",        "Pfdn4",
"Pgm5",          "Pitx1",         "Pitx2",         "Pkp1",
"Pld5",          "Plscr4",        "Pou3f3",        "Pou4f2",
"Pou4f3",        "Prdm12",        "Prmt8",         "Prok2",
"Prokr2",        "Prrxl1",        "Prune2",        "Ptf1a",
"Ptgis",         "Ptgs1",         "Ptpro",         "Rasl10b",
"Rax",           "Rbfox1",        "Reep2",         "Rfx6",
"Rgl3",          "Rhbdl1",        "Rhbg",          "Rhcg",
"Rhebl1",        "Ric3",          "Rspo2",         "Scara3",
"Scn3b",         "Sema6d",        "Sfrp2",         "Sfrp4",
"Sgk1",          "Shh",           "Six3",          "Six6",
"Skor2",         "Slc13a5",       "Slc16a7",       "Slc18a2",
"Slc18a3",       "Slc22a3",       "Slc26a10",      "Slc2a12",
"Slc32a1",       "Slc5a1",        "Slc5a7",        "Slc7a10",
"Slc7a14",       "Slit1",         "Slitrk3",       "Sntg1",
"Sorcs3",        "Sox14",         "Sox17",         "Sox9",
"Spag16",        "Spata18",       "Speg",          "Spon1",
"Srd5a2",        "Srsf12",        "Sv2b",          "Syt10",
"Syt6",          "Tac1",          "Tbx1",          "Tbx4",
"Tcf21",         "Tdh",           "Tdrd6",         "Tha1",
"Them7",         "Thsd7b",        "Thy1",          "Tmem132d",
"Tmem178b",      "Tmem235",       "Tnfrsf13c",     "Tnfrsf19",
"Trhde",         "Tril",          "Trim58",        "Trp73",
"Trpc4",         "Trpc6",         "Trpm3",         "Tspan11",
"Unc79",         "Vsx2",          "Vwc2",          "Wnt1",
"Wnt2",          "Wnt3",          "Wrap73",        "Zdhhc22",
"Zfp169",        "Zfp286",        "Zfp365",        "Zfp488",
"Zfp536",        "Zfp663",        "Zfp853",        "Zfp872",
"Zic4")

exe_hyper_genes_id <- sapply(exe_hyper_promoter_genes, function(x) rownames(subset(annotation_genes, gene_name == x)))

tpm_avg_log2_exe_hyper_df <- melt(tpm_avg_log2[exe_hyper_genes_id,order_samples])
tpm_avg_log2_exe_hyper_df$Var2 <- factor(tpm_avg_log2_exe_hyper_df$Var2, levels = order_samples)

## Combine gene sets and ExE hyper CGI genes
tpm_avg_log2_deg_groups_df_with_exe_hyper <- tpm_avg_log2_deg_groups_df
tpm_avg_log2_deg_groups_df_with_exe_hyper$gene_set <- paste(tpm_avg_log2_deg_groups_df_with_exe_hyper$L1, tpm_avg_log2_deg_groups_df_with_exe_hyper$deg, sep = "_")

tpm_avg_log2_exe_hyper_df$gene_set <- "ExE_hyper"
colnames(tpm_avg_log2_exe_hyper_df) <- c("gene", "variable", "value", "gene_set")

tpm_avg_log2_deg_groups_df_with_exe_hyper <- rbind(tpm_avg_log2_deg_groups_df_with_exe_hyper[,c("variable", "value", "gene_set")], tpm_avg_log2_exe_hyper_df[,c("variable", "value", "gene_set")])

tpm_avg_log2_deg_groups_df_with_exe_hyper$series <- factor(ifelse(as.character(tpm_avg_log2_deg_groups_df_with_exe_hyper$variable) %in% c("WT", "3BKO", "KDM2BKO", "RNF2KO", "EEDKO"), "KO", "inhibitor"), levels = c("inhibitor", "KO"))
tpm_avg_log2_deg_groups_df_with_exe_hyper$gene_set <- factor(tpm_avg_log2_deg_groups_df_with_exe_hyper$gene_set, levels = c("WT_DNMT1i_7d_acute_response_up", "WT_DNMT1i_7d_acute_response_down", "WT_EZH2i_7d_acute_response_up", "WT_EZH2i_7d_acute_response_down", "WT_DNMT1i_EZH2i_7d_acute_response_up", "WT_DNMT1i_EZH2i_7d_acute_response_down", "ExE_hyper"))

pdf("Extended_Data_Figure9c.pdf", width = 6, height = 20)
ggplot(data = tpm_avg_log2_deg_groups_df_with_exe_hyper, aes(x = variable, y = value, fill = variable)) + geom_boxplot(show.legend = FALSE) + theme_classic() + xlab("") + ylab("Log2 TPM") + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) + scale_fill_manual(values = colors_type) + geom_hline(yintercept = 1, lty = 2) + facet_grid(gene_set~series, scales = "free_x", space = "free_x")
dev.off()

###################################################################################################################
## Extended Data Figure 9d
###################################################################################################################

extended_data_figure_9d_source_data <- read.table(file = "Extended_Data_Figure9d_source_data.tsv", sep = "\\t", header = TRUE, stringsAsFactors = FALSE)

samples <- unique(extended_data_figure_9d_source_data$sample)

norm_mat_list <- lapply(samples, function(x) as.normalizedMatrix(as.matrix(subset(extended_data_figure_9d_source_data, sample == x)[,!colnames(extended_data_figure_9d_source_data) %in% c("sample", "gene_id")]), k_upstream=100, k_downstream=100, k_target=67, extend = c(5000, 5000), signal_name = x, target_name = "promoter"))
names(norm_mat_list) <- samples

col_fun <- list(
    "H3K27me3_TSC1_WT" = colorRamp2(c(0, 20), c("white", "violetred2")),
    "H2AK119ub1_TSC1_WT" = colorRamp2(c(0, 20), c("white", "deepskyblue3")),
    "H3K4me3_TSC1_WT" = colorRamp2(c(0, 50), c("white", "seagreen3")),
    "WGBS_TSC1_WT" = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    )

meta_height <- c(
    "H3K27me3_TSC1_WT" = 10,
    "H2AK119ub1_TSC1_WT" = 10,
    "H3K4me3_TSC1_WT" = 50,
    "WGBS_TSC1_WT" = 1
    )

genes_covered_minute_chip <- subset(extended_data_figure_9d_source_data, sample == samples[1])$gene_id
results_combined_inhibitor_acute_significant_filtered_heatmap <- lapply(results_combined_inhibitor_acute_significant_filtered, function(x) x[intersect(rownames(x), genes_covered_minute_chip),])
index_deg_dnmt1i <- sapply(rownames(results_combined_inhibitor_acute_significant_filtered_heatmap$WT_DNMT1i_7d_acute_response), function(x) which(genes_covered_minute_chip == x))
index_deg_ezh2i <- sapply(rownames(results_combined_inhibitor_acute_significant_filtered_heatmap$WT_EZH2i_7d_acute_response), function(x) which(genes_covered_minute_chip == x))
index_deg_dnmt1i_ezh2i <- sapply(rownames(results_combined_inhibitor_acute_significant_filtered_heatmap$WT_DNMT1i_EZH2i_7d_acute_response), function(x) which(genes_covered_minute_chip == x))

ht_list_deg_dnmt1i <- NULL
for(sample in samples)
{
    ht_list_deg_dnmt1i <- ht_list_deg_dnmt1i + EnrichedHeatmap(norm_mat_list[[sample]][index_deg_dnmt1i,], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col=c("dodgerblue", "firebrick"), lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6), row_split = ifelse(results_combined_inhibitor_acute_significant_filtered_heatmap$WT_DNMT1i_7d_acute_response$log2FoldChange > 0, "up", "down"))
}

ht_list_deg_ezh2i <- NULL
for(sample in samples)
{
    ht_list_deg_ezh2i <- ht_list_deg_ezh2i + EnrichedHeatmap(norm_mat_list[[sample]][index_deg_ezh2i,], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col=c("dodgerblue", "firebrick"), lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6), row_split = ifelse(results_combined_inhibitor_acute_significant_filtered_heatmap$WT_EZH2i_7d_acute_response$log2FoldChange > 0, "up", "down"))
}

ht_list_deg_dnmt1i_ezh2i <- NULL
for(sample in samples)
{
    ht_list_deg_dnmt1i_ezh2i <- ht_list_deg_dnmt1i_ezh2i + EnrichedHeatmap(norm_mat_list[[sample]][index_deg_dnmt1i_ezh2i,], use_raster = TRUE, raster_quality = 3, axis_name_rot= 90, col = col_fun[[sample]], column_title = sample, name = sample, show_heatmap_legend=T, top_annotation = HeatmapAnnotation(enrich=anno_enriched(show_heatmap_legend=F, gp=gpar(col=c("dodgerblue", "firebrick"), lwd=3), ylim=c(0, meta_height[sample]))), width=unit(4,"cm"), column_title_gp = gpar(fontsize = 6), row_split = ifelse(results_combined_inhibitor_acute_significant_filtered_heatmap$WT_DNMT1i_EZH2i_7d_acute_response$log2FoldChange > 0, "up", "down"))
}

pdf("Extended_Data_Figure9d.pdf", width = 20)
draw(ht_list_deg_dnmt1i)
draw(ht_list_deg_ezh2i)
draw(ht_list_deg_dnmt1i_ezh2i)
dev.off()

###################################################################################################################
## Extended Data Figure 10a
###################################################################################################################

## Placental marker genes
placental_marker_genes <- list(
ExE = c("Mapk1", "Cdx2", "Eomes", "Fgf4", "Fgfr2", "Esrrb", "Ascl2", "Rhox4b"),
TSC = c("Cdx2", "Esrrb", "Fgfr2", "Eomes", "Elf5"),
Ectoplacental_cone = c("Tfap2c", "Hand1", "Sdc1", "Brd3", "Tpbpa", "Flt1", "Ascl2", "Tcf12", "Tcf4", "Id1", "Id2", "Mapk1", "Ets2"),
ExE_chorion = c("Gcm1", "Esx1", "Dlx3", "Tead3", "Esrrb"),
Chorioallantoic_fusion = c("Bmp5", "Bmp7", "Itga4", "Lhx1", "Dnajb6", "T", "Vcam1"),
Chorioallantoic_branching = c("Dlx3", "Fgfr2", "Fosl1", "Fzd5", "Gab1", "Gcm1", "Grb2", "Gjc1", "Hgf", "Hsp90ab1", "Itgav", "Junb", "Lifr", "Map2k1", "Map2k2", "Map3k3", "Met", "Pdgfra", "Pdgfb", "Pparg", "Rxra", "Rxrb", "Sos1", "Vhl", "Wnt2"),
Syncytiotrophoblast = c("Gcm1", "Esx1", "Dlx3", "Tead3", "Esrrb", "Slc2a3", "Igfbp2"),
SynTI_precursor = c("Epha4", "Tgfa", "Tfrc", "Slc16a1", "Glis1", "Stra6"),
SynTI = c("Tfrc", "Slc16a1", "Glis1", "Stra6"),
SynTII_precursor = c("Ror2", "Lgr5"),
SynTII = c("Gcm1", "Synb", "Gcgr", "Vegfa", "Cdh1", "Igf1r"),
Vascularization_labyrinthe = c("Esx1", "Arnt", "Tfeb"),
Labyrinth_trophoblast_progenitor1 = c("Tcf7l1", "Ror2", "Lgr5", "Pvt1"),
Labyrinth_trophoblast_progenitor2 = c("Tcf7l1", "Pvt1", "Egfr"),
Junctional_zone_precursor1 = c("Cdh4", "Prune2", "Ncam1"),
Junctional_zone_precursor2 = c("Prune2", "Ncam1"),
TGC_precursor = c("Nos1ap", "Podxl", "Lepr", "Ctsq"),
TGC = c("Hand1", "Mdfi", "Cenpx", "Mdfi", "Prl3d1", "Prl3d3", "Ctsj", "Prl3b1", "Prl2c2", "Prl3d1", "Ctsq", "Lepr", "Nos1ap", "Tpbpa", "Pfpl", "Cts7", "Sult1e1", "Psg29", "Itga2", "Mmp9", "Plaur", "Psg18", "Mitf", "Flt1", "Prl8a9", "Plac8"),
Spongiotrophoblast_precursor = c("Tpbpa", "Flt1", "Ncoa6", "Hsf1", "Plac8", "Prl6a1", "Ascl2", "Phlda2", "Egfr", "Hsf1", "Mitf", "Flt1", "Prl8a9", "Slco2a1"),
Spongiotrophoblast = c("Adgre1", "Itgam", "Cd84", "Ctsb", "Hmox1", "Ccl5", "Cd74", "Cebpb", "Nr4a1"),
Glycogen_cell = c("Prune2", "Ncam1", "Igfbp7", "Pla2g4d", "Plac8", "Mitf", "Flt1"),
Prolactin_genes = c("Prl", "Prl2a1", "Prl2b1", "Prl2c2", "Prl2c3", "Prl2c5", "Prl3a1", "Prl3b1", "Prl3c1", "Prl3d1", "Prl3d2", "Prl3d3", "Prl4a1", "Prl5a1", "Prl6a1", "Prl7a1", "Prl7a2", "Prl7b1", "Prl7c1", "Prl7d1", "Prl8a1", "Prl8a2", "Prl8a6", "Prl8a8", "Prl8a9"),
DMR = c("A1cf", "Adcy7", "Ash2l", "Atf6b", "Cyb5r3", "D430020J02Rik", "Fpgs", "Gabpb1", "Gfod2", "Gm15760", "Gm16617", "Gsr", "Il22ra1", "Klra3", "Krt39", "Mir3069", "Mir3082", "Mir6391", "Mkln1os", "Piezo1", "Pstpip1", "Rgs7", "Rpl22", "Sap25", "Serpinb9c", "Srsf1", "Supt5", "Taf1c", "Tagap", "Tex14", "Tex19.1", "Tmem63c", "Vps72", "Zfp42")
)

placental_marker_genes_super_group <- list(
Early_progenitor = unique(c(placental_marker_genes$ExE, placental_marker_genes$TSC, placental_marker_genes$Ectoplacental_cone, placental_marker_genes$ExE_chorion)),
Syncytiotrophoblast = unique(c(placental_marker_genes$Syncytiotrophoblast, placental_marker_genes$SynTI_precursor, placental_marker_genes$SynTI, placental_marker_genes$SynTII_precursor, placental_marker_genes$SynTII)),
Labyrinth_trophoblast = unique(c(placental_marker_genes$Labyrinth_trophoblast_progenitor1, placental_marker_genes$Labyrinth_trophoblast_progenitor2)),
Junctional_zone = unique(c(placental_marker_genes$Junctional_zone_precursor1, placental_marker_genes$Junctional_zone_precursor2, placental_marker_genes$Glycogen_cell)),
TGC = unique(c(placental_marker_genes$TGC_precursor, placental_marker_genes$TGC)),
Prolactin_genes = placental_marker_genes$Prolactin_genes,
Spongiotrophoblast = unique(c(placental_marker_genes$Spongiotrophoblast_precursor, placental_marker_genes$Spongiotrophoblast)),
DMR = placental_marker_genes$DMR)

placental_marker_genes_super_group <- lapply(placental_marker_genes_super_group, function(x)
{
    curr_order_wt_expression <- order(tpm_avg_with_noncoding_log2[sapply(x, function(y) rownames(subset(annotation_genes_with_noncoding, gene_name == y))),"WT"], decreasing = TRUE)
    return(x[curr_order_wt_expression])
})

placental_marker_genes_super_group_id <- lapply(placental_marker_genes_super_group, function(x) sapply(x, function(y) rownames(subset(annotation_genes_with_noncoding, gene_name == y))))

breaksList <- seq(0, 12.2, by = 0.1)

pdf("Extended_Data_Figure10a.pdf", height = 8)
for (curr_set in names(placental_marker_genes_super_group_id))
{
    pheatmap::pheatmap(tpm_avg_with_noncoding_log2[placental_marker_genes_super_group_id[[curr_set]],order_samples], labels_row = placental_marker_genes_super_group[[curr_set]], cluster_rows = FALSE, cluster_cols = FALSE, breaks = breaksList, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)), main = curr_set, gaps_col = c(4,12))
}
dev.off()

###################################################################################################################
## Extended Data Figure 10b
###################################################################################################################

tpm_avg_log2_placental_marker_genes_super_group_df <- melt(do.call(rbind, lapply(names(placental_marker_genes_super_group_id), function(x)
{
    data <- data.frame(tpm_avg_with_noncoding_log2[placental_marker_genes_super_group_id[[x]],order_samples], check.names = FALSE)
    data$gene_set <- x
    return(data)
})))

tpm_avg_log2_placental_marker_genes_super_group_df$series <- factor(ifelse(as.character(tpm_avg_log2_placental_marker_genes_super_group_df$variable) %in% c("WT", "3BKO", "KDM2BKO", "RNF2KO", "EEDKO"), "KO", "inhibitor"), levels = c("inhibitor", "KO"))
tpm_avg_log2_placental_marker_genes_super_group_df$gene_set <- factor(tpm_avg_log2_placental_marker_genes_super_group_df$gene_set, levels = names(placental_marker_genes_super_group))

pdf("Extended_Data_Figure10b.pdf", width = 6, height = 10)
ggplot(data = subset(tpm_avg_log2_placental_marker_genes_super_group_df, gene_set %in% c("Early_progenitor", "TGC", "Prolactin_genes")), aes(x = variable, y = value, fill = variable)) + geom_boxplot(show.legend = FALSE) + theme_classic() + xlab("") + ylab("Log2 TPM") + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=90, hjust=1)) + scale_fill_manual(values = colors_type) + geom_hline(yintercept = 1, lty = 2) + facet_grid(gene_set~series, scales = "free_x", space = "free_x")
dev.off()
