#Geneset enrichment analysis
library(data.table)
library(clusterProfiler)
term2gene <- fread("/data1/ggulati/mberger/validationSets/FullTables_validationVersion2.0/GenesetAnalysis/msigdb_Tirosh_Table.txt", header = F)
term2gene <- term2gene[,-1]

geneList <- res2$log2FoldChange
names(geneList) <- rownames(res2)
geneList <- sort(geneList, decreasing = T)

term2gene <- read.table("geneset_analysis_forNeo1.txt", header = T)
term2gene <- term2gene[,-1]
score <- GSEA(geneList, TERM2GENE = term2gene,pvalueCutoff = 1)


noeffect1 <- gseaScores(geneList, as.character(unlist(df$gene))[which(df$group == "c")], 
                        exponent = 1, mode ="graph")$runningScore
#increased1 <- gseaScores(geneList, as.character(unlist(df$gene))[which(df$group == "b")], 
#           exponent = 1, mode ="graph")$runningScore
decreased1 <- gseaScores(geneList, as.character(unlist(df$gene))[which(df$group == "b")], 
                         exponent = 1, mode ="graph")$runningScore


t1 <- unique(term2gene$V2[grep("HALLMARK", term2gene$V2)])

runHallmark <- function(geneList){
t1 <- unique(term2gene$V2[grep("HALLMARK", term2gene$V2)])
for(i in 1:length(t1)){
padj <- GSEA(geneList, verbose = F, TERM2GENE = term2gene[which(term2gene$V2 %in% t1[i]),],pvalueCutoff = 1)@result$p.adjust
NES <- GSEA(geneList, verbose = F, TERM2GENE = term2gene[which(term2gene$V2 %in% t1[i]),],pvalueCutoff = 1)@result$NES
f1 <- cbind.data.frame(padj, NES)
ifelse(i == 1, f2 <- f1, f2 <- rbind.data.frame(f2, f1))
}
rownames(f2) <- t1
return(f2)
}

runHallmark(sort(diffgenes[[1]], decreasing = T)) -> a1; runHallmark(sort(diffgenes[[2]], decreasing = T)) -> a2; runHallmark(sort(diffgenes[[3]], decreasing = T)) -> a3
a1$padj <- p.adjust(a1$padj); a2$padj <- p.adjust(a2$padj); a3$padj <- p.adjust(a3$padj)
a1 <- a1[order(a1$padj),]; a2 <- a2[order(a2$padj),] ; a3 <- a3[order(a3$padj),]




