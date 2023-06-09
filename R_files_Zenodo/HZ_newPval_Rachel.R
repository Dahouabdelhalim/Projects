minCountsPerSample4Pval <- 10

noSamples <- 15253
outTable_zscores <- matrix(,nrow=0,ncol=(noSamples + 4))
outTable_ase <- matrix(,nrow=0,ncol=(noSamples + 4))
outTable_counts1 <- matrix(,nrow=0,ncol=(noSamples + 4))
outTable_counts2 <- matrix(,nrow=0,ncol=(noSamples + 4))

#compute final p-values
ase <- read.table("HZ_newPval_ase.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
zscores <- read.table("HZ_newPval_zscores.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")

#compute a z-score per gene based on the hybrid counts
counts1 <- read.table("HZ_newPval_counts1.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
counts2 <- read.table("HZ_newPval_counts2.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")

#remove samples with < X reads
cols <- c(5:ncol(ase))
for (col in cols) {
  idx <- which(counts1[,col] + counts2[,col] < minCountsPerSample4Pval)
  ase[idx,col] <- NA
}

#compute hybrid z-scores
counts_table_D50_H <- read.table("Rachel_humreffed_counts_filtered_D50_human_counts_clean.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
counts_table_D50_C <- read.table("Rachel_humreffed_counts_filtered_D50_chimp_counts_clean.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
counts_table_D100_H <- read.table("Rachel_humreffed_counts_filtered_D100_human_counts_clean.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
counts_table_D100_C <- read.table("Rachel_humreffed_counts_filtered_D100_chimp_counts_clean.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
counts_table_D150_H <- read.table("Rachel_humreffed_counts_filtered_D150_human_counts_clean.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
counts_table_D150_C <- read.table("Rachel_humreffed_counts_filtered_D150_chimp_counts_clean.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
zscores_hybrid_D50 <- matrix(,nrow=nrow(counts_table_D50_H),ncol=ncol(counts_table_D50_H))
zscores_hybrid_D100 <- matrix(,nrow=nrow(counts_table_D100_H),ncol=ncol(counts_table_D100_H))
zscores_hybrid_D150 <- matrix(,nrow=nrow(counts_table_D150_H),ncol=ncol(counts_table_D150_H))


for (gg in 1:nrow(counts_table_D50_H)) {
  for (i in 1:(ncol(counts_table_D50_H)-1)) {
    c1 <- as.numeric(counts_table_D50_H[gg,(i + 1)])
    totalCounts <- c1 + as.numeric(counts_table_D50_C[gg,(i + 1)])
    zscores_hybrid_D50[gg,(i + 1)] <- as.numeric((c1 - 0.5*totalCounts)/sqrt(totalCounts*0.5*0.5))
  }
  for (j in 1:(ncol(counts_table_D100_H)-1)) {
    c1 <- as.numeric(counts_table_D100_H[gg,(j + 1)])
    totalCounts <- c1 + as.numeric(counts_table_D100_C[gg,(j + 1)])
    zscores_hybrid_D100[gg,(j + 1)] <- as.numeric((c1 - 0.5*totalCounts)/sqrt(totalCounts*0.5*0.5))
  }
  for (k in 1:(ncol(counts_table_D150_H)-1)) {
    c1 <- as.numeric(counts_table_D150_H[gg,(k + 1)])
    totalCounts <- c1 + as.numeric(counts_table_D150_C[gg,(k + 1)])
    zscores_hybrid_D150[gg,(k + 1)] <- as.numeric((c1 - 0.5*totalCounts)/sqrt(totalCounts*0.5*0.5))
  }
}

#go over each gene, compute a p-value by comparing its ASE and z-score to the GTEx ASE and z-score
colnames <- c("gene","gene ID","pval","FDR")
ase_p_D50 <- matrix(NA,nrow=nrow(counts_table_D50_H),ncol=length(colnames))
ase_p_D100 <- matrix(NA,nrow=nrow(counts_table_D100_H),ncol=length(colnames))
ase_p_D150 <- matrix(NA,nrow=nrow(counts_table_D150_H),ncol=length(colnames))
zscores_p_D50 <- matrix(NA,nrow=nrow(counts_table_D50_H),ncol=length(colnames))
zscores_p_D100 <- matrix(NA,nrow=nrow(counts_table_D100_H),ncol=length(colnames))
zscores_p_D150 <- matrix(NA,nrow=nrow(counts_table_D150_H),ncol=length(colnames))
colnames(ase_p_D50) <- colnames
colnames(ase_p_D100) <- colnames
colnames(ase_p_D150) <- colnames
colnames(zscores_p_D50) <- colnames
colnames(zscores_p_D100) <- colnames
colnames(zscores_p_D150) <- colnames

#Need to download remaining files from David's folder and see what the format of ase table is.
#Also need to include ENSEMBL id information (sucks)
genes <- cbind(ase[,c(1:4)],matrix(,nrow=nrow(ase),ncol=1))

ASE_D50 <- unlist(read.table("Rachel_ASE_Values_Clean_D50.csv", header=T, sep=',', stringsAsFactors=F, quote="")$"D50_AvgLFC")
ASE_D100 <- unlist(read.table("Rachel_ASE_Values_Clean_D100.csv", header=T, sep=',', stringsAsFactors=F, quote="")$"D100_AvgLFC")
ASE_D150 <- unlist(read.table("Rachel_ASE_Values_Clean_D150.csv", header=T, sep=',', stringsAsFactors=F, quote="")$"D150_AvgLFC")

#clean Ensembl ID
for (gg in 1:nrow(genes)) {
  idxSep <- gregexpr(".",genes[gg,2],fixed=T)[[1]]
  genes[gg,5] <- substr(genes[gg,2],1,(idxSep-1))
}

ens <- read.table("gene_to_ens_ordered_rachel.txt", header=T, sep='\\t', stringsAsFactors=F, quote="")
for (gg in 1:nrow(counts_table_D50_H)) {
  print(paste0("Going over gene #",gg," out of ",nrow(counts_table_D50_H)," rows"))
  gene <- unlist(ens[gg, "ens"])
  idx_gene <- which(unlist(genes[,5]) == gene)
  if (length(idx_gene) > 0) {idx_gene <- idx_gene[1]}
  
  #compare z-scores using a U-test
  zscores_tempVec <- as.numeric(unlist(zscores[idx_gene,cols]))
  zscores_tempVec <- zscores_tempVec[which(!is.na(zscores_tempVec))]
  zscores_p_D50[gg,"gene"] <- unlist(counts_table_D50_H[gg,"gene"])
  zscores_p_D50[gg,"gene ID"] <- gene
  zscores_p_D100[gg,"gene"] <- unlist(counts_table_D100_H[gg,"gene"])
  zscores_p_D100[gg,"gene ID"] <- gene
  zscores_p_D150[gg,"gene"] <- unlist(counts_table_D150_H[gg,"gene"])
  zscores_p_D150[gg,"gene ID"] <- gene
  if (length(zscores_hybrid_D50) > 0 & length(zscores_tempVec) > 1) {
    if (sum(zscores_hybrid_D50[gg,c(2:8)],na.rm=T) != 0) {
      zscores_p_D50[gg,"pval"] <- wilcox.test(zscores_hybrid_D50[gg,c(2:8)],zscores_tempVec,alternative="two.sided")$"p.value"
    }
  }
  if (length(zscores_hybrid_D100) > 0 & length(zscores_tempVec) > 1) {
    if (sum(zscores_hybrid_D100[gg,c(2:10)],na.rm=T) != 0) {
      zscores_p_D100[gg,"pval"] <- wilcox.test(zscores_hybrid_D100[gg,c(2:10)],zscores_tempVec,alternative="two.sided")$"p.value"
    }
  }
  if (length(zscores_hybrid_D150) > 0 & length(zscores_tempVec) > 1) {
    if (sum(zscores_hybrid_D150[gg,c(2:8)],na.rm=T) != 0) {
      zscores_p_D150[gg,"pval"] <- wilcox.test(zscores_hybrid_D150[gg,c(2:8)],zscores_tempVec,alternative="two.sided")$"p.value"
    }
  }
  
  #compare ASE
  ase_tempVec <- as.numeric(ase[idx_gene,cols])
  idx <- which(!is.na(ase_tempVec))
  ase_tempVec <- ase_tempVec[idx]
  ase_p_D50[gg,"gene"] <- unlist(counts_table_D50_H[gg,"gene"])
  ase_p_D50[gg,"gene ID"] <- gene
  ase_p_D100[gg,"gene"] <- unlist(counts_table_D100_H[gg,"gene"])
  ase_p_D100[gg,"gene ID"] <- gene
  ase_p_D150[gg,"gene"] <- unlist(counts_table_D150_H[gg,"gene"])
  ase_p_D150[gg,"gene ID"] <- gene
  if (length(ASE_D50[gg]) > 0 & length(ase_tempVec) > 1) {
    if (!is.na(ASE_D50[gg])) {
      if (ASE_D50[gg] > 1) {
        ase_p_D50[gg,"pval"] <- (sum(ase_tempVec > ASE_D50[gg]) + sum((ase_tempVec) < 1/ASE_D50[gg]))/length(ase_tempVec)
      } else {
        ase_p_D50[gg,"pval"] <- (sum(ase_tempVec < ASE_D50[gg]) + sum((ase_tempVec) > 1/ASE_D50[gg]))/length(ase_tempVec)
      }
    }
  }
  if (length(ASE_D100[gg]) > 0 & length(ase_tempVec) > 1) {
    if (!is.na(ASE_D100[gg])) {
      if (ASE_D100[gg] > 1) {
        ase_p_D100[gg,"pval"] <- (sum(ase_tempVec > ASE_D100[gg]) + sum((ase_tempVec) < 1/ASE_D100[gg]))/length(ase_tempVec)
      } else {
        ase_p_D100[gg,"pval"] <- (sum(ase_tempVec < ASE_D100[gg]) + sum((ase_tempVec) > 1/ASE_D100[gg]))/length(ase_tempVec)
      }
    }
  }
  if (length(ASE_D150[gg]) > 0 & length(ase_tempVec) > 1) {
    if (!is.na(ASE_D150[gg])) {
      if (ASE_D150[gg] > 1) {
        ase_p_D150[gg,"pval"] <- (sum(ase_tempVec > ASE_D150[gg]) + sum((ase_tempVec) < 1/ASE_D150[gg]))/length(ase_tempVec)
      } else {
        ase_p_D150[gg,"pval"] <- (sum(ase_tempVec < ASE_D150[gg]) + sum((ase_tempVec) > 1/ASE_D150[gg]))/length(ase_tempVec)
      }
    }
  }
}

ase_p_D50[which(ase_p_D50[,3] == 0),3] <- 1/(ncol(ase)-4)
ase_p_D100[which(ase_p_D100[,3] == 0),3] <- 1/(ncol(ase)-4)
ase_p_D150[which(ase_p_D100[,3] == 0),3] <- 1/(ncol(ase)-4)


#FDR
idx <- which(!is.na(ase_p_D50[,"pval"]))
fdrs <- p.adjust(ase_p_D50[idx,"pval"], method = "BH")
ase_p_D50[idx,"FDR"] = fdrs

idx <- which(!is.na(ase_p_D100[,"pval"]))
fdrs <- p.adjust(ase_p_D100[idx,"pval"], method = "BH")
ase_p_D100[idx,"FDR"] = fdrs

idx <- which(!is.na(ase_p_D150[,"pval"]))
fdrs <- p.adjust(ase_p_D150[idx,"pval"], method = "BH")
ase_p_D150[idx,"FDR"] = fdrs

idx <- which(!is.na(zscores_p_D50[,"pval"]))
fdrs <- p.adjust(zscores_p_D50[idx,"pval"], method = "BH")
zscores_p_D50[idx,"FDR"] = fdrs

idx <- which(!is.na(zscores_p_D100[,"pval"]))
fdrs <- p.adjust(zscores_p_D100[idx,"pval"], method = "BH")
zscores_p_D100[idx,"FDR"] = fdrs

idx <- which(!is.na(zscores_p_D150[,"pval"]))
fdrs <- p.adjust(zscores_p_D150[idx,"pval"], method = "BH")
zscores_p_D150[idx,"FDR"] = fdrs

#WRITE
write.table(ase_p_D50, file="HZ_newPval_ase_Pvals_D50.txt", 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)

write.table(ase_p_D100, file="HZ_newPval_ase_Pvals_D100.txt", 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)

write.table(ase_p_D150, file="HZ_newPval_ase_Pvals_D150.txt", 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)

write.table(zscores_p_D50, file="HZ_newPval_zscores_Pvals_D50.txt", 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)

write.table(zscores_p_D100, file="HZ_newPval_zscores_Pvals_D100.txt", 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)

write.table(zscores_p_D150, file="HZ_newPval_zscores_Pvals_D150.txt", 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)
