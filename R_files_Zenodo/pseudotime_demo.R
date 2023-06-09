library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(pheatmap)
library(destiny)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(scran)
library(scater)
library(scuttle)
setwd('D:/GitHub/2021_cellreports_ontogeny')
# save.image('replot.rdata')
load('replot.rdata')
###
Fig.S7 <- readRDS('D:/Penn/Projects/ontogeny/Keren_pseudotime_20201228/dimplot.rds') # this is your downloaded dimplot.rds
DimPlot(Fig.S7, reduction = 'tsne')

###################################################################################################
# diffusion map
###################################################################################################
library(destiny)
top2000 <- head(VariableFeatures(Fig.S7), 2000)
top1000 <- head(VariableFeatures(Fig.S7), 1000)
# use SCTransform
Fig.S7 <- Seurat::SCTransform(Fig.S7)
# use seurat normalized counts.
Fig.S7 <- NormalizeData(Fig.S7, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(Fig.S7)
Fig.S7 <- ScaleData(Fig.S7, features = all.genes)

#################################################
## creat singlecellexperiment object.
Fig.S7_sce <- as.SingleCellExperiment(Fig.S7)
Fig.S7_sce <- computeSumFactors(Fig.S7_sce)
Fig.S7_sce <- logNormCounts(Fig.S7_sce)
# this is log normalized matrix, but subset with high variable genes.
dec <- modelGeneVar(Fig.S7_sce)
# Get the top 2000 genes.
top.hvgs2k <- getTopHVGs(dec, n=2000)
top.hvgs1k <- getTopHVGs(dec, n=1000)

Fig.S7.mtx <-t(as.matrix(GetAssayData(Fig.S7, assay = 'SCT', slot = 'scale.data')[which(rownames(Fig.S7) %in% top2000),]))
class(Fig.S7.mtx)
# this takes time.
sigmas <- find_sigmas(Fig.S7.mtx, verbose= TRUE)
optimal_sigma(sigmas)
# sigma small, more concentrated, sigma large, more disperse.
mtx_dm <-DiffusionMap(Fig.S7.mtx, n_pcs = 50, sigma = 20, density_norm = TRUE, verbose= TRUE) 
plot(mtx_dm)
##
plot(mtx_dm,1:2)
###################################################################################################
## calculate dpt, pseudotime.
dpt <- DPT(mtx_dm)
plot(dpt, col_by='branch', divide=3, dcs= c(-1,-3,2), pch=20)
###################################################################################################
library(pheatmap)
# create data.frame for easy plotting of DPT
tsne_df <- Embeddings(Fig.S7, reduction = 'tsne')
##
info_table <- cbind(tsne_df, mtx_dm@eigenvectors) %>% as.data.frame()

info_table$dpt8 <- dpt$DPT8


## so here we choose the 8th pseudotime which meets our expection.
ggplot2::ggplot(info_table, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(fill = dpt$DPT8), shape = 21,  size = 2) +
  scale_fill_viridis_b() +
  theme_classic()
###################################################################################################
### get the high variable genes with FDRã€‚
# Get all genes with FDR below 5%.
top.hvgs4 <- getTopHVGs(dec, fdr.threshold=0.05)
## order the barcode by dpt8
info_table <- info_table %>% rownames_to_column(var = 'barcode') %>% dplyr::arrange(dpt8)
summary(info_table$dpt8)

mtx_hm <- as.matrix(GetAssayData(Fig.S7, assay = 'SCT', slot = 'scale.data')[top.hvgs4, info_table$barcode])
head(colnames(mtx_hm))
## quantile the data to remove outliers
quantile(mtx_hm, probs = seq(0, 1, 0.05))
quantile(mtx_hm, 0.95)
mtx_hm[mtx_hm > quantile(mtx_hm, 0.95)] <- quantile(mtx_hm, 0.95)
mtx_hm[mtx_hm < quantile(mtx_hm, 0.05)] <- quantile(mtx_hm, 0.05)

## check the ordered plot.
hm_plot <- pheatmap(mtx_hm,
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    treeheight_row = 5,
                    show_rownames = FALSE, 
                    show_colnames = FALSE
                    )

pdf('heatmap_pseudotime.pdf', height = 6, width = 8)
hm_plot
dev.off()

###################################################################################################
## gene expression along the pseudotime

mtx_exp <- t(as.matrix(GetAssayData(Fig.S7, assay = 'SCT', slot = 'data')[, info_table$barcode]))
mtx_exp[1:6, 1:5]

mtx_exp_plot <- cbind(mtx_exp, info_table[, c('barcode', 'dpt8')])
colnames(mtx_exp_plot)

plot_genes <- c('EYA1', 'GATA4', 'GATA6', 'HOXB9', 'LHX9', 'OSR1')
mtx_exp_plot_gather <- tidyr::gather(mtx_exp_plot[,c(plot_genes, 'dpt8')], key = 'gene', value = 'LogExp', -dpt8)
head(mtx_exp_plot_gather)

ggplot(mtx_exp_plot_gather, aes(x=dpt8, y=LogExp, group=gene, color=gene)) +
  geom_smooth(se = FALSE, span = 0.1)+
  scale_color_brewer(palette = 'Set1') +
  ggtitle("Gene expression with pseudotime") +
  ylim(0,2)+
  theme_classic() +
  ylab("Log normalized expression")
  











