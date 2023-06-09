library(Seurat)
library(ggplot2)
library(DESeq2)
library(dplyr)
library(rtracklayer)
library(pheatmap)
library(amap)
library(BiocParallel)

bulk <- read.table("merged_htseq_count.tsv", sep = "\\t", header = TRUE)

head(bulk)
dim(bulk)

load("IDSymbol.Rdata")
IDsymbol[nrow(IDsymbol)+1,] <- c("eGFP", "eGFP")

bulk <- as.matrix(bulk)
bulk <- bulk[match(IDsymbol$ID, bulk[,1]),]
bulk <- bulk[,-1]
bulk <- apply(bulk,2,as.numeric)
rownames(bulk) <- IDsymbol$Symbol

bulk <- round(limma::avereps(bulk))
bulk <- as.data.frame(bulk)

sham <- bulk[,1:6]
tac7 <- bulk[,7:12]

sample_num <- 3

colnames(sham) <- c("sample1_neg", "sample2_neg", "sample3_neg",
                    "sample1_pos", "sample2_pos", "sample3_pos")

colnames(tac7) <- c("sample5_neg", "sample6_neg", "sample7_neg",
                    "sample5_pos", "sample6_pos", "sample7_pos")

endofb <- cbind(sham,tac7)

geneFilter <- apply(sham, 1, function(x){
  sum(x >= 10) > sample_num
})
sham <- sham[geneFilter,]

geneFilter <- apply(tac7, 1, function(x){
  sum(x >= 10) > sample_num
})
tac7 <- tac7[geneFilter,]

geneFilter <- apply(endofb, 1, function(x){
  sum(x >= 10) > sample_num
})
endofb <- endofb[geneFilter,]

condition_sham <- factor(c(rep("sham_neg", sample_num), rep("sham_pos", sample_num)))
condition_tac7 <- factor(c(rep("tac7_neg", sample_num), rep("tac7_pos", sample_num)))
condition_endofb <- factor(c(rep("sham_neg", sample_num), rep("sham_pos", sample_num),
                             rep("tac7_neg", sample_num), rep("tac7_pos", sample_num)))

sample_sham <- factor(unlist(strsplit(colnames(sham), "_"))[c(seq(1,4*sample_num,2))])
sample_tac7 <- factor(unlist(strsplit(colnames(tac7), "_"))[c(seq(1,4*sample_num,2))])
sample_endofb <- factor(unlist(strsplit(colnames(endofb), "_"))[c(seq(1,8*sample_num,2))])

coldata_sham <- data.frame(row.names=colnames(sham), condition_sham, sample_sham)
coldata_tac7 <- data.frame(row.names=colnames(tac7), condition_tac7, sample_tac7)
coldata_endofb <- data.frame(row.names=colnames(endofb), condition_endofb, sample_endofb)

dds_sham <- DESeqDataSetFromMatrix(countData=sham, colData=coldata_sham, 
                                   design= ~ sample_sham+condition_sham)

dds_tac7 <- DESeqDataSetFromMatrix(countData=tac7, colData=coldata_tac7, 
                                   design= ~ sample_tac7+condition_tac7)

dds_endofb <- DESeqDataSetFromMatrix(countData=endofb, colData=coldata_endofb, 
                                     design= ~ condition_endofb)

library(sva)
adjusted_sham <- ComBat_seq(counts = counts(dds_sham), batch = dds_sham$sample_sham)
dds_adjusted_sham <- DESeqDataSetFromMatrix(countData=adjusted_sham, colData=coldata_sham, 
                                            design= ~ sample_sham+condition_sham)
dds_adjusted_sham <- DESeq(dds_adjusted_sham)

adjusted_tac7 <- ComBat_seq(counts = counts(dds_tac7), batch = dds_tac7$sample_tac7)
dds_adjusted_tac7 <- DESeqDataSetFromMatrix(countData=adjusted_tac7, colData=coldata_tac7, 
                                            design= ~ sample_tac7+condition_tac7)
dds_adjusted_tac7 <- DESeq(dds_adjusted_tac7)

adjusted_endofb <- ComBat_seq(counts = counts(dds_endofb), batch = dds_endofb$sample_endofb)
dds_adjusted_endofb <- DESeqDataSetFromMatrix(countData=adjusted_endofb, colData=coldata_endofb, 
                                              design= ~ condition_endofb)
dds_adjusted_endofb <- DESeq(dds_adjusted_endofb)

# Figure 6B
rld_endofb <- rlog(dds_adjusted_endofb, blind = FALSE)

pdf(file = "endofb_qc-pca_modified.pdf", width = 20, height = 10, family = "ArialMT")
plotPCA(rld_endofb, intgroup=c("condition_endofb"), ntop=3000)+
  theme_bw()+
  theme(
    legend.background=element_blank(), legend.key=element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )
dev.off()

res <- results(dds_adjusted_sham)
sham_posvsneg <- merge(as.data.frame(res), counts(dds_adjusted_sham, normalized=TRUE), by="row.names", sort=FALSE)
names(sham_posvsneg)[1] <- "Gene"
sham_posvsneg <- arrange(sham_posvsneg, desc(log2FoldChange))
sham_posvsneg <- na.omit(sham_posvsneg)
sham_posvsneg_Sig <- subset(sham_posvsneg, sham_posvsneg$pvalue < 0.05)
sham_posvsneg_Sig_up <- subset(sham_posvsneg_Sig, sham_posvsneg_Sig$log2FoldChange >= 0.25)
sham_posvsneg_Sig_dw <- subset(sham_posvsneg_Sig, sham_posvsneg_Sig$log2FoldChange <= (-1)*0.25)


res <- results(dds_adjusted_tac7)
tac7_posvsneg <- merge(as.data.frame(res), counts(dds_adjusted_tac7, normalized=TRUE), by="row.names", sort=FALSE)
names(tac7_posvsneg)[1] <- "Gene"
tac7_posvsneg <- arrange(tac7_posvsneg, desc(log2FoldChange))
tac7_posvsneg <- na.omit(tac7_posvsneg)
tac7_posvsneg_Sig <- subset(tac7_posvsneg, tac7_posvsneg$pvalue < 0.05)
tac7_posvsneg_Sig_up <- subset(tac7_posvsneg_Sig, tac7_posvsneg_Sig$log2FoldChange >= 0.25)
tac7_posvsneg_Sig_dw <- subset(tac7_posvsneg_Sig, tac7_posvsneg_Sig$log2FoldChange <= (-1)*0.25)

sham_posvsneg_Sig_up_filtered <- sham_posvsneg_Sig_up[1:500,]
tac7_posvsneg_Sig_up_filtered <- tac7_posvsneg_Sig_up[1:500,]

# Figure 6C
sham_posvsneg$sig <- "nosig"
sham_posvsneg[which(sham_posvsneg$log2FoldChange > 0.25 & sham_posvsneg$pvalue < 0.05), ]$sig <- "up"
sham_posvsneg[which(sham_posvsneg$log2FoldChange < -0.25 & sham_posvsneg$pvalue < 0.05), ]$sig <- "down"

label_sham <- c("Thbs4", "Fmod", "Comp", "Postn",
                "Cilp", "Gsn", "Serping1", "Gpx3", "Entpd2", "Adamtsl2",
                "Ctnnb1", "eGFP")

sham_posvsneg$label <- ""

sham_posvsneg[which(sham_posvsneg$Gene %in% label_sham), ]$label <- sham_posvsneg[which(sham_posvsneg$Gene %in% label_sham), ]$Gene

pdf(file = "sham_posvsneg_volcano.pdf", family = "ArialMT")
ggplot(data = sham_posvsneg, aes(x = log2FoldChange, y = MinMax(-log10(pvalue), 0, 30), label=label)) + 
  geom_point(aes(colour = sig), alpha=0.9) +
  scale_color_manual(values=c("up" = "red","down" = "blue", "nosig" = "black"))+
  scale_x_continuous("log2FoldChange") +
  scale_y_continuous("-log10(pvalue)")+
  theme_bw()+
  theme(
    legend.background=element_blank(), legend.key=element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )+NoLegend()+ggrepel::geom_text_repel(max.overlaps = Inf, point.padding = 0.15, box.padding = 3)
dev.off()

tac7_posvsneg$sig <- "nosig"
tac7_posvsneg[which(tac7_posvsneg$log2FoldChange > 0.25 & tac7_posvsneg$pvalue < 0.05), ]$sig <- "up"
tac7_posvsneg[which(tac7_posvsneg$log2FoldChange < -0.25 & tac7_posvsneg$pvalue < 0.05), ]$sig <- "down"

label_tac7 <- c("Thbs4", "Fmod", "Cthrc1", "Col12a1", "Comp", "Postn", "Ddah1", 
                "Cilp", "Ifitm3",
                "Ctnnb1", "eGFP", "Wif1",
                "Cyp26a1", "Crabp2", "Gcnt1", "Col22a1", "Gm41541", "Fgf9", "Chad") 

tac7_posvsneg$label <- ""

tac7_posvsneg[which(tac7_posvsneg$Gene %in% label_tac7), ]$label <- tac7_posvsneg[which(tac7_posvsneg$Gene %in% label_tac7), ]$Gene

pdf(file = "tac7_posvsneg_volcano.pdf", family = "ArialMT")
ggplot(data = tac7_posvsneg, aes(x = log2FoldChange, y = MinMax(-log10(pvalue), 0, 30), label=label)) + 
  geom_point(aes(colour = sig), alpha=0.9) +
  scale_color_manual(values=c("up" = "red","down" = "blue", "nosig" = "black"))+
  scale_x_continuous("log2FoldChange") +
  scale_y_continuous("-log10(pvalue)")+
  theme_bw()+
  theme(
    legend.background=element_blank(), legend.key=element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )+NoLegend()+ggrepel::geom_text_repel(max.overlaps = Inf, point.padding = 0.15, box.padding = 3)
dev.off()

# Figure 6D
sham_GSEA <- dplyr::select(sham_posvsneg, Gene, log2FoldChange)[-c(1:2),]
tac7_GSEA <- dplyr::select(tac7_posvsneg, Gene, log2FoldChange)[-c(8,17),]

write.table(sham_GSEA, file = "./prerank_sham.rnk", sep = "\\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(tac7_GSEA, file = "./prerank_tac7.rnk", sep = "\\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

regional <- read.table("merged_htseq_count.tsv", sep = "\\t", header = TRUE)

head(regional)
dim(regional)

load("IDSymbol.Rdata")
IDsymbol[nrow(IDsymbol)+1,] <- c("eGFP", "eGFP")

regional <- as.matrix(regional)
regional <- regional[match(IDsymbol$ID, regional[,1]),]
regional <- regional[,-1]
regional <- apply(regional,2,as.numeric)
rownames(regional) <- IDsymbol$Symbol

regional <- round(limma::avereps(regional))
regional <- as.data.frame(regional)

lower <- regional[,1:6]
upper <- regional[,7:12]
neg <- cbind(lower[,1:3], upper[,1:3])
pos <- cbind(lower[,4:6], upper[,4:6])
sample_num <- 3
endofb <- cbind(lower,upper)

geneFilter <- apply(lower, 1, function(x){
  sum(x >= 10) > sample_num
})
lower <- lower[geneFilter,]

geneFilter <- apply(upper, 1, function(x){
  sum(x >= 10) > sample_num
})
upper <- upper[geneFilter,]

geneFilter <- apply(neg, 1, function(x){
  sum(x >= 10) > sample_num
})
neg <- neg[geneFilter,]

geneFilter <- apply(pos, 1, function(x){
  sum(x >= 10) > sample_num
})
pos <- pos[geneFilter,]

geneFilter <- apply(endofb, 1, function(x){
  sum(x >= 10) > sample_num
})
endofb <- endofb[geneFilter,]

condition_lower <- factor(c(rep("lower_neg", sample_num), rep("lower_pos", sample_num)))
condition_upper <- factor(c(rep("upper_neg", sample_num), rep("upper_pos", sample_num)))
condition_neg <- factor(c(rep("neg_lower", sample_num), rep("neg_upper", sample_num)))
condition_pos <- factor(c(rep("pos_lower", sample_num), rep("pos_upper", sample_num)))
condition_endofb <- factor(c(rep("lower_neg", sample_num), rep("lower_pos", sample_num),
                             rep("upper_neg", sample_num), rep("upper_pos", sample_num)))

sample_lower <- factor(rep(c("672232", "672510", "672514"), 2))
sample_upper <- factor(rep(c("672232", "672510", "672514"), 2))
sample_neg <- factor(rep(c("672232", "672510", "672514"), 2))
sample_pos <- factor(rep(c("672232", "672510", "672514"), 2))
sample_endofb <- factor(rep(c("672232", "672510", "672514"), 4))

coldata_lower <- data.frame(row.names=colnames(lower), condition_lower, sample_lower)
coldata_upper <- data.frame(row.names=colnames(upper), condition_upper, sample_upper)
coldata_neg <- data.frame(row.names=colnames(neg), condition_neg, sample_neg)
coldata_pos <- data.frame(row.names=colnames(pos), condition_pos, sample_pos)
coldata_endofb <- data.frame(row.names=colnames(endofb), condition_endofb, sample_endofb)

dds_lower <- DESeqDataSetFromMatrix(countData=lower, colData=coldata_lower, 
                                    design= ~ sample_lower + condition_lower)

dds_upper <- DESeqDataSetFromMatrix(countData=upper, colData=coldata_upper, 
                                    design= ~ sample_upper + condition_upper)

dds_neg <- DESeqDataSetFromMatrix(countData=neg, colData=coldata_neg, 
                                  design= ~ sample_neg + condition_neg)

dds_pos <- DESeqDataSetFromMatrix(countData=pos, colData=coldata_pos, 
                                  design= ~ sample_pos + condition_pos)

dds_endofb <- DESeqDataSetFromMatrix(countData=endofb, colData=coldata_endofb, 
                                     design= ~ sample_endofb + condition_endofb)

library(sva)
adjusted_lower <- ComBat_seq(counts = counts(dds_lower), batch = dds_lower$sample_lower)
dds_adjusted_lower <- DESeqDataSetFromMatrix(countData=adjusted_lower, colData=coldata_lower, 
                                             design= ~ sample_lower + condition_lower)
dds_adjusted_lower <- DESeq(dds_adjusted_lower)

adjusted_upper <- ComBat_seq(counts = counts(dds_upper), batch = dds_upper$sample_upper)
dds_adjusted_upper <- DESeqDataSetFromMatrix(countData=adjusted_upper, colData=coldata_upper, 
                                             design= ~ sample_upper + condition_upper)
dds_adjusted_upper <- DESeq(dds_adjusted_upper)

adjusted_neg <- ComBat_seq(counts = counts(dds_neg), batch = dds_neg$sample_neg)
dds_adjusted_neg <- DESeqDataSetFromMatrix(countData=adjusted_neg, colData=coldata_neg, 
                                           design= ~ sample_neg + condition_neg)
dds_adjusted_neg <- DESeq(dds_adjusted_neg)

adjusted_pos <- ComBat_seq(counts = counts(dds_pos), batch = dds_pos$sample_pos)
dds_adjusted_pos <- DESeqDataSetFromMatrix(countData=adjusted_pos, colData=coldata_pos, 
                                           design= ~ sample_pos + condition_pos)
dds_adjusted_pos <- DESeq(dds_adjusted_pos)

adjusted_endofb <- ComBat_seq(counts = counts(dds_endofb), batch = dds_endofb$sample_endofb)
dds_adjusted_endofb <- DESeqDataSetFromMatrix(countData=adjusted_endofb, colData=coldata_endofb, 
                                              design= ~ condition_endofb)
dds_adjusted_endofb <- DESeq(dds_adjusted_endofb)

# Figure 6G
rld_endofb <- rlog(dds_adjusted_endofb, blind = FALSE)

pdf(file = "endofb_qc-pca_modified.pdf", width = 20, height = 10, family = "ArialMT")
plotPCA(rld_endofb, intgroup=c("condition_endofb"), ntop=3000)+
  theme_bw()+
  theme(
    legend.background=element_blank(), legend.key=element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )
dev.off()

res <- results(dds_adjusted_lower)
lower_posvsneg <- merge(as.data.frame(res), counts(dds_adjusted_lower, normalized=TRUE), by="row.names", sort=FALSE)
names(lower_posvsneg)[1] <- "Gene"
lower_posvsneg <- arrange(lower_posvsneg, desc(log2FoldChange))
lower_posvsneg <- na.omit(lower_posvsneg)
lower_posvsneg_Sig <- subset(lower_posvsneg, lower_posvsneg$pvalue < 0.05)
lower_posvsneg_Sig_up <- subset(lower_posvsneg_Sig, lower_posvsneg_Sig$log2FoldChange >= 0.25)
lower_posvsneg_Sig_dw <- subset(lower_posvsneg_Sig, lower_posvsneg_Sig$log2FoldChange <= (-1)*0.25)


res <- results(dds_adjusted_upper)
upper_posvsneg <- merge(as.data.frame(res), counts(dds_adjusted_upper, normalized=TRUE), by="row.names", sort=FALSE)
names(upper_posvsneg)[1] <- "Gene"
upper_posvsneg <- arrange(upper_posvsneg, desc(log2FoldChange))
upper_posvsneg <- na.omit(upper_posvsneg)
upper_posvsneg_Sig <- subset(upper_posvsneg, upper_posvsneg$pvalue < 0.05)
upper_posvsneg_Sig_up <- subset(upper_posvsneg_Sig, upper_posvsneg_Sig$log2FoldChange >= 0.25)
upper_posvsneg_Sig_dw <- subset(upper_posvsneg_Sig, upper_posvsneg_Sig$log2FoldChange <= (-1)*0.25)


res <- results(dds_adjusted_neg)
neg_uppervslower <- merge(as.data.frame(res), counts(dds_adjusted_neg, normalized=TRUE), by="row.names", sort=FALSE)
names(neg_uppervslower)[1] <- "Gene"
neg_uppervslower <- arrange(neg_uppervslower, desc(log2FoldChange))
neg_uppervslower <- na.omit(neg_uppervslower)
neg_uppervslower_Sig <- subset(neg_uppervslower, neg_uppervslower$pvalue < 0.05)
neg_uppervslower_Sig_up <- subset(neg_uppervslower_Sig, neg_uppervslower_Sig$log2FoldChange >= 0.25)
neg_uppervslower_Sig_dw <- subset(neg_uppervslower_Sig, neg_uppervslower_Sig$log2FoldChange <= (-1)*0.25)


res <- results(dds_adjusted_pos)
pos_uppervslower <- merge(as.data.frame(res), counts(dds_adjusted_pos, normalized=TRUE), by="row.names", sort=FALSE)
names(pos_uppervslower)[1] <- "Gene"
pos_uppervslower <- arrange(pos_uppervslower, desc(log2FoldChange))
pos_uppervslower <- na.omit(pos_uppervslower)
pos_uppervslower_Sig <- subset(pos_uppervslower, pos_uppervslower$pvalue < 0.05)
pos_uppervslower_Sig_up <- subset(pos_uppervslower_Sig, pos_uppervslower_Sig$log2FoldChange >= 0.25)
pos_uppervslower_Sig_dw <- subset(pos_uppervslower_Sig, pos_uppervslower_Sig$log2FoldChange <= (-1)*0.25)

lower_posvsneg_Sig_up_filtered <- lower_posvsneg_Sig_up[1:500,]
upper_posvsneg_Sig_up_filtered <- upper_posvsneg_Sig_up[1:500,]

# Figure6 H,I
upper_posvsneg$sig <- "nosig"
upper_posvsneg[which(upper_posvsneg$log2FoldChange > 0.25 & upper_posvsneg$pvalue < 0.05), ]$sig <- "up"
upper_posvsneg[which(upper_posvsneg$log2FoldChange < -0.25 & upper_posvsneg$pvalue < 0.05), ]$sig <- "down"

label_upper <- c("Thbs4", "Fmod", "Cthrc1", "Col12a1", "Abi3bp", "Comp", "Postn", "Ddah1",
                 "Cilp", 
                 "Ctnnb1", "eGFP",
                 "Epyc", "Ism1", "Lonrf2", "Mobp", "Npy4r", "Tent5b", "Lrrc1", "Scin", "Gcnt1", "Tbx5", "Wnt2")

upper_posvsneg$label <- ""

upper_posvsneg[which(upper_posvsneg$Gene %in% label_upper), ]$label <- upper_posvsneg[which(upper_posvsneg$Gene %in% label_upper), ]$Gene

pdf(file = "upper_posvsneg_volcano.pdf", family = "ArialMT")
ggplot(data = upper_posvsneg, aes(x = log2FoldChange, y = MinMax(-log10(pvalue), 0, 30), label=label)) + 
  geom_point(aes(colour = sig), alpha=0.9) +
  scale_color_manual(values=c("up" = "red","down" = "blue", "nosig" = "black"))+
  scale_x_continuous("log2FoldChange") +
  scale_y_continuous("-log10(pvalue)")+
  theme_bw()+
  theme(
    legend.background=element_blank(), legend.key=element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )+NoLegend()+ggrepel::geom_text_repel(max.overlaps = Inf, point.padding = 0.15, box.padding = 3)
dev.off()

pos_uppervslower$sig <- "nosig"
pos_uppervslower[which(pos_uppervslower$log2FoldChange > 0.25 & pos_uppervslower$pvalue < 0.05), ]$sig <- "up"
pos_uppervslower[which(pos_uppervslower$log2FoldChange < -0.25 & pos_uppervslower$pvalue < 0.05), ]$sig <- "down"

label_upper <- c("Tnni3", "Hey1", "Bmp4", "Pecam1", "Egfl7", "Ptn", "Cdh5", "Tie1", "Cdh13", "Kdr", "Aplnr")

pos_uppervslower$label <- ""

pos_uppervslower[which(pos_uppervslower$Gene %in% label_upper), ]$label <- pos_uppervslower[which(pos_uppervslower$Gene %in% label_upper), ]$Gene

pdf(file = "pos_uppervslower_volcano.pdf", family = "ArialMT")
ggplot(data = pos_uppervslower, aes(x = log2FoldChange, y = MinMax(-log10(pvalue), 0, 30), label=label)) + 
  geom_point(aes(colour = sig), alpha=0.9) +
  scale_color_manual(values=c("up" = "red","down" = "blue", "nosig" = "black"))+
  scale_x_continuous("log2FoldChange") +
  scale_y_continuous("-log10(pvalue)")+
  theme_bw()+
  theme(
    legend.background=element_blank(), legend.key=element_blank(),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank()
  )+NoLegend()+ggrepel::geom_text_repel(max.overlaps = Inf, point.padding = 0.15, box.padding = 3)
dev.off()

# Supplementary B
library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(org.Mm.eg.db)

egoall_sham <- clusterProfiler::enrichGO(sham_posvsneg_Sig_up_filtered$Gene, 
                                         OrgDb = org.Mm.eg.db, ont = "BP", 
                                         pAdjustMethod = 'fdr', pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.05, 
                                         keyType = 'SYMBOL')
egosimp_sham <- clusterProfiler::simplify(egoall_sham)

egoall_tac7 <- clusterProfiler::enrichGO(tac7_posvsneg_Sig_up_filtered$Gene, 
                                         OrgDb = org.Mm.eg.db, ont = "BP", 
                                         pAdjustMethod = 'fdr', pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.05, 
                                         keyType = 'SYMBOL')
egosimp_tac7 <- clusterProfiler::simplify(egoall_tac7)

egoall_upper <- clusterProfiler::enrichGO(upper_posvsneg_Sig_up_filtered$Gene, OrgDb = org.Mm.eg.db, ont = "BP", 
                                          pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                                          keyType = 'SYMBOL')
egosimp_upper <- clusterProfiler::simplify(egoall_upper)

egoall_lower <- clusterProfiler::enrichGO(lower_posvsneg_Sig_up_filtered$Gene, OrgDb = org.Mm.eg.db, ont = "BP", 
                                          pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                                          keyType = 'SYMBOL')
egosimp_lower <- clusterProfiler::simplify(egoall_lower)

egoall_neg <- clusterProfiler::enrichGO(neg_uppervslower_Sig_up$Gene, OrgDb = org.Mm.eg.db, ont = "BP", 
                                        pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                                        keyType = 'SYMBOL')
egosimp_neg <- clusterProfiler::simplify(egoall_neg)

egoall_pos <- clusterProfiler::enrichGO(pos_uppervslower_Sig_up$Gene, OrgDb = org.Mm.eg.db, ont = "BP", 
                                        pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                                        keyType = 'SYMBOL')
egosimp_pos <- clusterProfiler::simplify(egoall_pos)

egoall <- list(sham=egosimp_sham, tac7=egosimp_tac7, lower=egosimp_lower, 
               upper=egosimp_upper, neg=egosimp_neg, pos=egosimp_pos)
egoall_merge <- merge_result(egoall)
egotree <- pairwise_termsim(egoall_merge)
save(egotree, file = "egotree.Rdata")

pdf(file = "enrichGO_treeplot.pdf", height=60, width=40, family="ArialMT")
treeplot(egotree, showCategory = 50)
dev.off()

# Supplementary C
length(tac7_posvsneg$Gene)
length(upper_posvsneg$Gene)
length(lower_posvsneg$Gene)
length(neg_uppervslower$Gene)
length(pos_uppervslower$Gene)

length(intersect(upper_posvsneg_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))
length(upper_posvsneg_Sig_up$Gene)-length(intersect(upper_posvsneg_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))
length(tac7_posvsneg_Sig_up$Gene)-length(intersect(upper_posvsneg_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))

length(intersect(lower_posvsneg_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))
length(lower_posvsneg_Sig_up$Gene)-length(intersect(lower_posvsneg_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))
length(tac7_posvsneg_Sig_up$Gene)-length(intersect(lower_posvsneg_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))

length(intersect(pos_uppervslower_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))
length(pos_uppervslower_Sig_up$Gene)-length(intersect(pos_uppervslower_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))
length(tac7_posvsneg_Sig_up$Gene)-length(intersect(pos_uppervslower_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))


length(intersect(neg_uppervslower_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))
length(neg_uppervslower_Sig_up$Gene)-length(intersect(neg_uppervslower_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))
length(tac7_posvsneg_Sig_up$Gene)-length(intersect(neg_uppervslower_Sig_up$Gene, tac7_posvsneg_Sig_up$Gene))

# Supplementary D; Figure 6 E,J
getcommongene <- function(pathway, group=2){
  discard <- c()
  keep <- c()
  if (group==6) {
    for (gene in pathway) {
      if (gene %in% sham_posvsneg$Gene) {
        if (gene %in% tac7_posvsneg$Gene) {
          if (gene %in% upper_posvsneg$Gene) {
            if (gene %in% lower_posvsneg$Gene) {
              keep <- c(gene, keep)
            }else{
              discard <- c(gene, discard)
            }
          }else{
            discard <- c(gene, discard)
          }
        }else{
          discard <- c(gene, discard)
        }
      }else{
        discard <- c(gene, discard)
      }
    }
    pathway <- keep
    return(pathway)
  }else{
    if (group==2) {
      for (gene in pathway) {
        if (gene %in% sham_posvsneg$Gene) {
          if (gene %in% tac7_posvsneg$Gene) {
            keep <- c(gene, keep)
          }else{
            discard <- c(gene, discard)
          }
        }else{
          discard <- c(gene, discard)
        }
      }
      pathway <- keep
      return(pathway)
    }else{
      if (group==4) {
        for (gene in pathway) {
          if (gene %in% upper_posvsneg$Gene) {
            if (gene %in% lower_posvsneg$Gene) {
              keep <- c(gene, keep)
            }else{
              discard <- c(gene, discard)
            }
          }else{
            discard <- c(gene, discard)
          }
        }
        pathway <- keep
        return(pathway)
      }
    }
  }
}

Wound_healing <- egotree@compareClusterResult[c(15,169,208,245,502,604),]$geneID
Wound_healing <- unique(unlist(strsplit(Wound_healing, "/")))
Wound_healing <- getcommongene(Wound_healing, group = 6)

Ossification <- egotree@compareClusterResult[c(8,152,213,255,507,553),]$geneID
Ossification <- unique(unlist(strsplit(Ossification, "/")))
Ossification <- getcommongene(Ossification, group = 6)

ECM_organization <- egotree@compareClusterResult[c(1,159,196,280),]$geneID
ECM_organization <- unique(unlist(strsplit(ECM_organization, "/")))
ECM_organization <- getcommongene(ECM_organization, group = 6)

Wnt_signaling <- egotree@compareClusterResult[c(46,174,240),]$geneID
Wnt_signaling <- unique(unlist(strsplit(Wnt_signaling, "/")))
Wnt_signaling <- getcommongene(Wnt_signaling, group = 6)

Cell_adhension <- egotree@compareClusterResult[c(11,387),]$geneID
Cell_adhension <- unique(unlist(strsplit(Cell_adhension, "/")))
Cell_adhension <- getcommongene(Cell_adhension, group = 6)

Cell_migration <- egotree@compareClusterResult[c(94,191,181,385,263),]$geneID
Cell_migration <- unique(unlist(strsplit(Cell_migration, "/")))
Cell_migration <- getcommongene(Cell_migration, group = 6)

MItochondrial <- egotree@compareClusterResult[c(147,187,163,189),]$geneID
MItochondrial <- unique(unlist(strsplit(MItochondrial, "/")))
MItochondrial <- getcommongene(MItochondrial, group = 2)

Cell_cycle <- egotree@compareClusterResult[c(172),]$geneID
Cell_cycle <- unique(unlist(strsplit(Cell_cycle, "/")))
Cell_cycle <- getcommongene(Cell_cycle, group = 2)

Ribosomal <- egotree@compareClusterResult[c(149,157,148),]$geneID
Ribosomal <- unique(unlist(strsplit(Ribosomal, "/")))
Ribosomal <- getcommongene(Ribosomal, group = 2)

Endo_cushion <- egotree@compareClusterResult[c(242),]$geneID
Endo_cushion <- unique(unlist(strsplit(Endo_cushion, "/")))
Endo_cushion <- getcommongene(Endo_cushion, group = 4)

Muscle_adaptation <- egotree@compareClusterResult[c(262,251),]$geneID
Muscle_adaptation <- unique(unlist(strsplit(Muscle_adaptation, "/")))
Muscle_adaptation <- getcommongene(Muscle_adaptation, group = 4)

Interleukin <- egotree@compareClusterResult[c(437,562),]$geneID
Interleukin <- unique(unlist(strsplit(Interleukin, "/")))
Interleukin <- getcommongene(Interleukin, group = 4)

Cardiac_development <- egotree@compareClusterResult[c(442,451,434, 430),]$geneID
Cardiac_development <- unique(unlist(strsplit(Cardiac_development, "/")))
Cardiac_development <- getcommongene(Cardiac_development, group = 4)

Angiogenesis <- egotree@compareClusterResult[c(425,545,427,558,426,422,424,465,547),]$geneID
Angiogenesis <- unique(unlist(strsplit(Angiogenesis, "/")))
Angiogenesis <- getcommongene(Angiogenesis, group = 4)

pathway <- list(Wound_healing=Wound_healing, Ossification=Ossification, 
                ECM_organization=ECM_organization,
                Wnt_signaling=Wnt_signaling, Cell_adhension=Cell_adhension, 
                Cell_migration=Cell_migration,
                MItochondrial=MItochondrial, Cell_cycle=Cell_cycle, Ribosomal=Ribosomal, 
                Endo_cushion=Endo_cushion,
                Muscle_adaptation=Muscle_adaptation, Interleukin=Interleukin, 
                Cardiac_development=Cardiac_development,
                Angiogenesis=Angiogenesis)

library(pheatmap)
plotheatmap <- function(pathway, group=2){
  if (group==2) {
    plotheat <- list()
    for (path in c(names(pathway)[1:9])) {
      gene <- pathway[[path]]
      sham <- as.data.frame(t(scale(t(sham_posvsneg[match(gene, sham_posvsneg$Gene), c(8:13)]))))
      sham$Gene <- sham_posvsneg[match(gene, sham_posvsneg$Gene), c(1)]
      sham$log2FoldChange <- sham_posvsneg[match(gene, sham_posvsneg$Gene), c(3)]
      sham <- sham[sham$log2FoldChange > 0,]
      rownames(sham) <- sham[,7]
      tac7 <- as.data.frame(t(scale(t(tac7_posvsneg[match(gene, tac7_posvsneg$Gene), c(8:13)]))))
      tac7$Gene <- tac7_posvsneg[match(gene, tac7_posvsneg$Gene), c(1)]
      tac7$log2FoldChange <- tac7_posvsneg[match(gene, tac7_posvsneg$Gene), c(3)]
      tac7 <- tac7[tac7$log2FoldChange > 0,]
      rownames(tac7) <- tac7[,7]
      keep <- intersect(rownames(sham), rownames(tac7))
      sham <- sham[keep,1:6]
      rownames(sham) <- keep
      tac7 <- tac7[keep,1:6]
      rownames(tac7) <- keep
      assign(paste0("heat", path), cbind(sham,tac7))
      plotheat[[path]] <- as.matrix(get(paste0("heat", path)))
    }
    return(plotheat)
    
  }else{
    if (group==4) {
      plotheat <- list()
      for (path in c(names(pathway)[c(1:6,10:11)])) {
        gene <- pathway[[path]]
        lower <- as.data.frame(t(scale(t(lower_posvsneg[match(gene, lower_posvsneg$Gene), c(8:13)]))))
        lower$Gene <- lower_posvsneg[match(gene, lower_posvsneg$Gene), c(1)]
        lower$log2FoldChange <- lower_posvsneg[match(gene, lower_posvsneg$Gene), c(3)]
        lower <- lower[lower$log2FoldChange > 0,]
        rownames(lower) <- lower$Gene
        upper <- as.data.frame(t(scale(t(upper_posvsneg[match(gene, upper_posvsneg$Gene), c(8:13)]))))
        upper$Gene <- upper_posvsneg[match(gene, upper_posvsneg$Gene), c(1)]
        upper$log2FoldChange <- upper_posvsneg[match(gene, upper_posvsneg$Gene), c(3)]
        upper <- upper[upper$log2FoldChange > 0,]
        rownames(upper) <- upper$Gene
        keep <- intersect(rownames(lower), rownames(upper))
        lower <- lower[keep,1:6]
        rownames(lower) <- keep
        upper <- upper[keep,1:6]
        rownames(upper) <- keep
        assign(paste0("heat", path), cbind(lower,upper))
        plotheat[[path]] <- as.matrix(get(paste0("heat", path)))
      }
      
      for (path in c(names(pathway)[c(12:14)])) {
        gene <- pathway[[path]]
        gene <- intersect(gene, neg_uppervslower$Gene)
        gene <- intersect(gene, pos_uppervslower$Gene)
        neg <- as.data.frame(t(scale(t(neg_uppervslower[match(gene, neg_uppervslower$Gene), c(8:13)]))))
        neg$Gene <- neg_uppervslower[match(gene, neg_uppervslower$Gene), c(1)]
        neg$log2FoldChange <- neg_uppervslower[match(gene, neg_uppervslower$Gene), c(3)]
        neg <- neg[neg$log2FoldChange > 0,]
        rownames(neg) <- neg$Gene
        pos <- as.data.frame(t(scale(t(pos_uppervslower[match(gene, pos_uppervslower$Gene), c(8:13)]))))
        pos$Gene <- pos_uppervslower[match(gene, pos_uppervslower$Gene), c(1)]
        pos$log2FoldChange <- pos_uppervslower[match(gene, pos_uppervslower$Gene), c(3)]
        pos <- pos[pos$log2FoldChange > 0,]
        rownames(pos) <- pos$Gene
        keep <- intersect(rownames(neg), rownames(pos))
        neg <- neg[keep,1:6]
        rownames(neg) <- keep
        pos <- pos[keep,1:6]
        rownames(pos) <- keep
        assign(paste0("heat", path), cbind(neg,pos))
        assign(paste0("heat", path), get(paste0("heat", path))[,c(1:3,7:9,4:6,10:12)])
        plotheat[[path]] <- as.matrix(get(paste0("heat", path)))
      }
      return(plotheat)
    }
  }
}
plot2 <- plotheatmap(pathway = pathway, group = 2)
plot4 <- plotheatmap(pathway = pathway, group = 4)

pdf(file = "heatplot2_woundhealing.pdf", height = 40, width = 10)
pheatmap(plot2$Wound_healing, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot2_Ossification.pdf", height = 40, width = 10)
pheatmap(plot2$Ossification, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot2_ECM.pdf", height = 40, width = 10)
pheatmap(plot2$ECM_organization, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot2_Wnt.pdf", height = 40, width = 10)
pheatmap(plot2$Wnt_signaling, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot2_celladhension.pdf", height = 40, width = 10)
pheatmap(plot2$Cell_adhension, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot2_cellmigration.pdf", height = 40, width = 10)
pheatmap(plot2$Cell_migration, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot2_Mito.pdf", height = 40, width = 10)
pheatmap(plot2$MItochondrial, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot2_CC.pdf", height = 40, width = 10)
pheatmap(plot2$Cell_cycle, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot2_Ribo.pdf", height = 40, width = 10)
pheatmap(plot2$Ribosomal, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))

dev.off()

pdf(file = "heatplot4_Wound.pdf", height = 40, width = 10)
pheatmap(plot4$Wound_healing, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_Ossification.pdf", height = 40, width = 10)
pheatmap(plot4$Ossification, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_ECM.pdf", height = 40, width = 10)
pheatmap(plot4$ECM_organization, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_Wnt.pdf", height = 40, width = 10)
pheatmap(plot4$Wnt_signaling, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_celladhension.pdf", height = 40, width = 10)
pheatmap(plot4$Cell_adhension, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_cellmigration.pdf", height = 40, width = 10)
pheatmap(plot4$Cell_migration, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_Endocushion.pdf", height = 40, width = 10)
pheatmap(plot4$Endo_cushion, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_Muscleadaptation.pdf", height = 40, width = 10)
pheatmap(plot4$Muscle_adaptation, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_Interleukin.pdf", height = 40, width = 10)
pheatmap(plot4$Interleukin, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_cardiacdevelop.pdf", height = 40, width = 10)
pheatmap(plot4$Cardiac_development, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()

pdf(file = "heatplot4_Angiogenesis.pdf", height = 40, width = 10)
pheatmap(plot4$Angiogenesis, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 35, cellheight = 5, annotation_names_row = TRUE, annotation_names_col = FALSE,
         border_color = FALSE, show_colnames = FALSE,
         color=colorRampPalette(c("#5000f1", "white", "#eb0048"))(1000))
dev.off()
