
a<-read.table("yourbetaLRT_EVEtransposedmatrix_output.txt", sep='\\t', header=T)
LRT<-a$LRT
genes<-a$gene_ortholog
df = 1
P_diverge <- 1 - pchisq(LRT, df = df)
P_diverse <- pchisq(LRT, df = df)
FDR_diverge <- p.adjust(P_diverge, "fdr")

FDR_diverse <- p.adjust(P_diverse, "fdr")
EVE_diverge <- as.data.frame(cbind(genes, LRT, P_diverge, FDR_diverge))
EVE_diverse <- as.data.frame(cbind(genes, LRT, P_diverse, FDR_diverse))
colnames(EVE_diverge) <- c("Gene", "LRT", "P", "FDR")
colnames(EVE_diverse) <- c("Gene", "LRT", "P", "FDR")
rate<-0.05
EVE_diverge.sub <- subset(EVE_diverge, EVE_diverge$P < rate)
EVE_diverge.sub <- as.data.frame(EVE_diverge.sub)
colnames(EVE_diverge.sub)<- c("Gene", "LRT", "P", "FDR")
EVE_diverse.sub <- subset(EVE_diverse, EVE_diverse$P < rate)
EVE_diverse.sub <- as.data.frame(EVE_diverse.sub)
colnames(EVE_diverse.sub)<- c("Gene", "LRT", "P", "FDR")
FDR_diverge.len <- length(EVE_diverge.sub[,1])
EVE_diverge.len <- length(EVE_diverge[,1])
FDR_diverse.len <- length(EVE_diverse.sub[,1])
EVE_diverse.len <- length(EVE_diverse[,1])
print(paste0(FDR_diverge.len, " genes have higher variance among than within lineages at ", rate, "FDR."))
print(paste0(FDR_diverse.len, " genes have higher variance within than among lineages at ", rate, "FDR."))
write.table(EVE_diverge.sub, file ="EVEdiverge_p0.05_HOSTALL.txt", sep="\\t")
write.table(EVE_diverse.sub, file ="EVEdiverse_p0.05_HOSTALL.txt", sep="\\t")
