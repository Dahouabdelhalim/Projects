library(tximport)
library(DESeq2)

tx2gene <- read.table("hg38/hg38.v103.protein_coding.tx2gene.txt")
names(tx2gene) <- c("TXNAME", "GENEID")

files<-c("data/salmon/FB_RNAseq_rep1/quant.sf",
         "data/salmon/FB_RNAseq_rep2/quant.sf",
         "data/salmon/CMC_RNAseq_rep1/quant.sf",
         "data/salmon/CMC_RNAseq_rep2/quant.sf")
names(files) <- c("FB_rep1", "FB_rep2", "CMC_rep1", "CMC_rep2")

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

for (cell in c("FB_rep1", "FB_rep2", "CMC_rep1", "CMC_rep2")){
    tpm <- data.frame(target_id=rownames(txi$abundanc),tpm=txi$abundanc[,cell])
    write.table(tpm, paste0("data/", cell, "_TPM.txt"), sep="\\t", quote=FALSE, row.names=FALSE)
}

table <- data.frame(name=c("FB_rep1", "FB_rep2", "CMC_rep1", "CMC_rep2"),
                    cell=c("FB", "FB", "CMC", "CMC"),
                    condition=c("rep1", "rep2", "rep1", "rep2"))
dds <- DESeqDataSetFromTximport(txi, colData=table, design=~cell)
dds <- DESeq(dds)
res <- results(dds, contrast=c("cell", "CMC", "FB"))
resOrdered <- res[order(res$padj),]
resOrdered@listData$resid <- rownames(resOrdered)
fdf <- data.frame(resOrdered)[,c("resid", "log2FoldChange", "padj")]
fdf$log2FoldChange[is.na(fdf$log2FoldChange)] <- 0
fdf$padj[is.na(fdf$padj)] <- 1
write.table(fdf, "data/FB2CMC_degenes.csv", sep="\\t", row.names=FALSE, quote=FALSE)
