outpath <- "/Users/dgokhman/Dropbox (Weizmann Institute)/POSTDOC/HUMANZEE/RESULTS/new pvals - Hunter/"
filepath <- "/Users/davidgo/Dropbox (Weizmann Institute)/POSTDOC/HUMANZEE/RESULTS/new pvals - Hunter/GTEx_v8_WASP/phASER_GTEx_v8_matrix_WASP.txt"
ASEfile_path <- "/Users/dgokhman/Dropbox (Weizmann Institute)/POSTDOC/HUMANZEE/ASE.xlsx"

noSamples <- 15253 #number of samples in the file
minCountsPerSample4Pval <- 10 #minimum counts per sample for the p-val computation

#prepare empty matrices
outTable_zscores <- matrix(,nrow=0,ncol=(noSamples + 4))
outTable_ase <- matrix(,nrow=0,ncol=(noSamples + 4))
outTable_counts1 <- matrix(,nrow=0,ncol=(noSamples + 4))
outTable_counts2 <- matrix(,nrow=0,ncol=(noSamples + 4))


#READ FILE, PARSE AND SAVE COUNTS
con = file(filepath, "r")
headers <- read.delim(con,nrows=1,header=F,sep="\\t")
ll <- 1
while (TRUE) {
  print(paste0("Going over line: ",ll))
  line <- read.delim(con,nrows=1,header=F,sep="\\t")
  if (length(line) == 0) {
    break
  }
  tempGeneData <- line[1:4]
  #take allele1 and 2 counts
  separator <- "|"
  counts1 <- matrix(,nrow=1, ncol=noSamples)
  counts2 <- matrix(,nrow=1, ncol=noSamples)
  for (sample in 1:noSamples) {
    sampCounts <- unlist(line[sample+4])
    idxSep <- gregexpr(separator,sampCounts,fixed=T)[[1]]
    c1 <- as.numeric(substr(sampCounts,1,(idxSep-1)))
    c2 <- as.numeric(substr(sampCounts,idxSep+1,nchar(as.character((sampCounts)))))
    totalCounts <- sum(c1,c2,na.rm=T)
    if (totalCounts > 0) {
      #pvals[1,sample] <- binom.test(counts1,totalCounts,p=0.5)$p.value
      #Binomial mean is n*p, i.e., totalCounts*0.5
      #Binomial variation is n*p*q, i.e., totalCounts*0.5*(1-0.5)
      #zscores[1,sample] <- (c1 - 0.5*totalCounts)/sqrt(totalCounts*0.5*0.5)
      #ase[1,sample] <- c1/c2
      counts1[1,sample] <- c1
      counts2[1,sample] <- c2
    }
  }
  #tempLine_zscores <- cbind(line[1:4],zscores)
  #outTable_zscores <- rbind(outTable_zscores,tempLine_zscores)
  tempLine_counts1 <- cbind(tempGeneData,counts1)
  outTable_counts1 <- rbind(outTable_counts1,tempLine_counts1)
  tempLine_counts2 <- cbind(tempGeneData,counts2)
  outTable_counts2 <- rbind(outTable_counts2,tempLine_counts2)
  ll <- ll + 1
}
close(con)

#outTable_counts1 <- read.table(paste0(outpath,"HZ_newPval_counts1.txt"), header=T, sep='\\t', stringsAsFactors=F, quote="")
#outTable_counts2 <- read.table(paste0(outpath,"HZ_newPval_counts2.txt"), header=T, sep='\\t', stringsAsFactors=F, quote="")


#deleted row 41600 because the pasting of files didn't match. Due to this, there's one less gene
#in the GTEx data
write.table(outTable_counts1, file=paste0(outpath,"HZ_newPval_counts1.txt"), 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)
write.table(outTable_counts2, file=paste0(outpath,"HZ_newPval_counts2.txt"), 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)



#COMPUTE Z-SCORES AND ASE
#pvals <- matrix(,nrow=1, ncol=noSamples)
zscores <- cbind(outTable_counts1[,c(1:4)],matrix(,nrow=nrow(outTable_counts1), ncol=noSamples))
ase <- cbind(outTable_counts1[,c(1:4)],matrix(,nrow=nrow(outTable_counts1), ncol=noSamples))
for (ll in 1:nrow(outTable_counts1)) {
  print(paste0("Going over line: ",ll))
  for (sample in 1:noSamples) {
    c1 <- as.numeric(outTable_counts1[ll,(sample+4)])
    c2 <- as.numeric(outTable_counts2[ll,(sample+4)])
    totalCounts <- sum(c1,c2,na.rm=T)
    if (totalCounts > 0) {
      #pvals[ll,(sample+4)] <- binom.test(counts1,totalCounts,p=0.5)$p.value
      #Binomial mean is n*p, i.e., totalCounts*0.5
      #Binomial variation is n*p*q, i.e., totalCounts*0.5*(1-0.5)
      zscores[ll,(sample+4)] <- (c1 - 0.5*totalCounts)/sqrt(totalCounts*0.5*0.5)
      ase[ll,(sample+4)] <- c1/c2
    }
  }
}

write.table(zscores, file=paste0(outpath,"HZ_newPval_zscores.txt"), 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)
write.table(ase, file=outpath(,"HZ_newPval_ase.txt"), 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)



#compute final p-values
ase <- read.table(paste0(outpath,"HZ_newPval_ase.txt"), header=T, sep='\\t', stringsAsFactors=F, quote="")
zscores <- read.table(paste0(outpath,"HZ_newPval_zscores.txt"), header=T, sep='\\t', stringsAsFactors=F, quote="")
library("readxl")
ASE_table <- read_excel(ASEfile_path,skip=3,col_names=T)

#compute a z-score per gene based on the hybrid counts
counts1 <- read.table(past0(outpath,"HZ_newPval_counts1.txt"), header=T, sep='\\t', stringsAsFactors=F, quote="")
counts2 <- read.table(paste0(outpath,"HZ_newPval_counts2.txt"), header=T, sep='\\t', stringsAsFactors=F, quote="")

#remove samples with < X reads
cols <- c(5:ncol(ase))
for (col in cols) {
  idx <- which(counts1[,col] + counts2[,col] < minCountsPerSample4Pval)
  ase[idx,col] <- NA
}

#compute hybrid z-scores
counts_table <- read_excel(ASEfile_path,skip=3,col_names=T,sheet="Counts")
colnames_iPSC <- c("gene","HL1-25_1","HL1-25_2","HL1-29_1","HL1-29_2","HL1-30_1","HL1-30_2","HL2-9_1","HL2-9_2","HL2-16_1","HL2-16_2")
colnames_CNCC <- c("gene","HL1-30_1","HL1-30_2","HL1-30_3")
zscores_hybrid_iPSC <- matrix(,nrow=nrow(counts_table),ncol=length(colnames_iPSC))
zscores_hybrid_CNCC <- matrix(,nrow=nrow(counts_table),ncol=length(colnames_CNCC))
counts_table_iPSC_H <- counts_table[,c("gene","ref_counts...19","ref_counts...23","ref_counts...27","ref_counts...31","ref_counts...35","ref_counts...39","ref_counts...43","ref_counts...47","ref_counts...51","ref_counts...55")]
counts_table_iPSC_C <- counts_table[,c("gene","alt_counts...20","alt_counts...24","alt_counts...28","alt_counts...32","alt_counts...36","alt_counts...40","alt_counts...44","alt_counts...48","alt_counts...52","alt_counts...56")]
counts_table_CNCC_H <- counts_table[,c("gene","ref_counts...107","ref_counts...111","ref_counts...115")]
counts_table_CNCC_C <- counts_table[,c("gene","alt_counts...108","alt_counts...112","alt_counts...116")]
for (gg in 1:nrow(counts_table)) {
  for (iPSC_samp in 1:10) {
    c1 <- as.numeric(counts_table_iPSC_H[gg,(iPSC_samp + 1)])
    totalCounts <- c1 + as.numeric(counts_table_iPSC_C[gg,(iPSC_samp + 1)])
    zscores_hybrid_iPSC[gg,(iPSC_samp + 1)] <- as.numeric((c1 - 0.5*totalCounts)/sqrt(totalCounts*0.5*0.5))
  }
  for (CNCC_samp in 1:3) {
    c1 <- as.numeric(counts_table_CNCC_H[gg,(CNCC_samp + 1)])
    totalCounts <- c1 + as.numeric(counts_table_CNCC_C[gg,(CNCC_samp + 1)])
    zscores_hybrid_CNCC[gg,(CNCC_samp + 1)] <- (c1 - 0.5*totalCounts)/sqrt(totalCounts*0.5*0.5)
  }
}

#go over each gene, compute a p-value by comparing its ASE and z-score to the GTEx ASE and z-score
colnames <- c("gene","gene ID","pval","FDR")
ase_p_iPSC <- matrix(NA,nrow=nrow(ASE_table),ncol=length(colnames))
ase_p_CNCC <- matrix(NA,nrow=nrow(ASE_table),ncol=length(colnames))
zscores_p_iPSC <- matrix(NA,nrow=nrow(ASE_table),ncol=length(colnames))
zscores_p_CNCC <- matrix(NA,nrow=nrow(ASE_table),ncol=length(colnames))
colnames(ase_p_iPSC) <- colnames
colnames(ase_p_CNCC) <- colnames
colnames(zscores_p_iPSC) <- colnames
colnames(zscores_p_CNCC) <- colnames
genes <- cbind(ase[,c(1:4)],matrix(,nrow=nrow(ase),ncol=1))

ASE_iPSC <- unlist(2^ASE_table[,"log2FoldChange_ASE_iPS"])
ASE_CNCC <- unlist(2^ASE_table[,"log2FoldChange_ASE_CNCC"])

#clean Ensembl ID
for (gg in 1:nrow(genes)) {
  idxSep <- gregexpr(".",genes[gg,2],fixed=T)[[1]]
  genes[gg,5] <- substr(genes[gg,2],1,(idxSep-1))
}


for (gg in 1:nrow(ASE_table)) {
  print(paste0("Going over gene #",gg," out of ",nrow(ASE_table)," rows"))
  gene <- unlist(ASE_table[gg,"Ensembl (combined - Ensembl + Gene ORGANizer)"])
  idx_gene <- which(unlist(genes[,5]) == gene)
  if (length(idx_gene) > 0) {idx_gene <- idx_gene[1]}
  
  #compare z-scores using a U-test
  zscores_tempVec <- as.numeric(unlist(zscores[idx_gene,cols]))
  zscores_tempVec <- zscores_tempVec[which(!is.na(zscores_tempVec))]
  zscores_p_iPSC[gg,"gene"] <- unlist(ASE_table[gg,"gene"])
  zscores_p_iPSC[gg,"gene ID"] <- gene
  zscores_p_CNCC[gg,"gene"] <- unlist(ASE_table[gg,"gene"])
  zscores_p_CNCC[gg,"gene ID"] <- gene
  if (length(zscores_hybrid_iPSC) > 0 & length(zscores_tempVec) > 1) {
    if (sum(zscores_hybrid_iPSC[gg,c(2:11)],na.rm=T) != 0) {
      zscores_p_iPSC[gg,"pval"] <- wilcox.test(zscores_hybrid_iPSC[gg,c(2:11)],zscores_tempVec,alternative="two.sided")$"p.value"
    }
  }
  if (length(zscores_hybrid_CNCC) > 0 & length(zscores_tempVec) > 1) {
    if (sum(zscores_hybrid_CNCC[gg,c(2:4)],na.rm=T) != 0) {
      zscores_p_CNCC[gg,"pval"] <- wilcox.test(zscores_hybrid_CNCC[gg,c(2:4)],zscores_tempVec,alternative="two.sided")$"p.value"
    }
  }
  
  #compare ASE
  ase_tempVec <- as.numeric(ase[idx_gene,cols])
  idx <- which(!is.na(ase_tempVec))
  ase_tempVec <- ase_tempVec[idx]
  ase_p_iPSC[gg,"gene"] <- unlist(ASE_table[gg,"gene"])
  ase_p_iPSC[gg,"gene ID"] <- gene
  ase_p_CNCC[gg,"gene"] <- unlist(ASE_table[gg,"gene"])
  ase_p_CNCC[gg,"gene ID"] <- gene
  if (length(ASE_iPSC[gg]) > 0 & length(ase_tempVec) > 1) {
    if (!is.na(ASE_iPSC[gg])) {
      if (ASE_iPSC[gg] > 1) {
        ase_p_iPSC[gg,"pval"] <- (sum(ase_tempVec > ASE_iPSC[gg]) + sum((ase_tempVec) < 1/ASE_iPSC[gg]))/length(ase_tempVec)
      } else {
        ase_p_iPSC[gg,"pval"] <- (sum(ase_tempVec < ASE_iPSC[gg]) + sum((ase_tempVec) > 1/ASE_iPSC[gg]))/length(ase_tempVec)
      }
    }
  }
  if (length(ASE_CNCC[gg]) > 0 & length(ase_tempVec) > 1) {
    if (!is.na(ASE_CNCC[gg])) {
      if (ASE_CNCC[gg] > 1) {
        ase_p_CNCC[gg,"pval"] <- (sum(ase_tempVec > ASE_CNCC[gg]) + sum(ase_tempVec < 1/ASE_CNCC[gg]))/length(ase_tempVec)
      } else {
        ase_p_CNCC[gg,"pval"] <- (sum(ase_tempVec < ASE_CNCC[gg]) + sum(ase_tempVec > 1/ASE_CNCC[gg]))/length(ase_tempVec)
      }
    }
  }
}

ase_p_iPSC[which(ase_p_iPSC[,3] == 0),3] <- 1/(ncol(ase)-4)
ase_p_CNCC[which(ase_p_CNCC[,3] == 0),3] <- 1/(ncol(ase)-4)


#FDR
idx <- which(!is.na(ase_p_iPSC[,"pval"]))
fdrs <- p.adjust(ase_p_iPSC[idx,"pval"], method = "BH")
ase_p_iPSC[idx,"FDR"] = fdrs

idx <- which(!is.na(ase_p_CNCC[,"pval"]))
fdrs <- p.adjust(ase_p_CNCC[idx,"pval"], method = "BH")
ase_p_CNCC[idx,"FDR"] = fdrs

idx <- which(!is.na(zscores_p_iPSC[,"pval"]))
fdrs <- p.adjust(zscores_p_iPSC[idx,"pval"], method = "BH")
zscores_p_iPSC[idx,"FDR"] = fdrs

idx <- which(!is.na(zscores_p_CNCC[,"pval"]))
fdrs <- p.adjust(zscores_p_CNCC[idx,"pval"], method = "BH")
zscores_p_CNCC[idx,"FDR"] = fdrs

#WRITE FILES
write.table(ase_p_iPSC, file=paste0(outpath,"HZ_newPval_ase_Pvals_iPSC.txt"), 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)

write.table(ase_p_CNCC, file=paste0(outpath,"HZ_newPval_ase_Pvals_CNCC.txt"), 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)

write.table(zscores_p_iPSC, file=paste0(outpath,"HZ_newPval_zscores_Pvals_iPSC.txt"), 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)

write.table(zscores_p_CNCC, file=paste0(outpath,"HZ_newPval_zscores_Pvals_CNCC.txt"), 
            sep="\\t", quote = F,na = "", row.names = F, col.names = T)




#gg <- 6842 #EVC2

#gene <- unlist(ASE_table[gg,"Ensembl (combined - Ensembl + Gene ORGANizer)"])
#idx_gene <- which(unlist(genes[,5]) == gene)
#if (length(idx_gene) > 0) {idx_gene <- idx_gene[1]}
#zscores_tempVec <- as.numeric(unlist(zscores[idx_gene,cols]))
#zscores_tempVec <- zscores_tempVec[which(!is.na(zscores_tempVec))]
#ase_tempVec <- as.numeric(ase[idx_gene,cols])
#idx <- which(!is.na(ase_tempVec))
#ase_tempVec <- ase_tempVec[idx]
#ase_tempVec[ase_tempVec > 10] <- 10

#hist(ase_tempVec,breaks=100)
#hist(zscores_tempVec,breaks=100)

#write.table(counts1_samp, file=paste0(outpath,"counts1_samp.txt"), 
#            sep="\\t", quote = F,na = "", row.names = F, col.names = T)



#original matlab script
#%s=std(log2((e+1)./repmat(sum(e+1),length(e(:,1)),1))')';
       
#       %log_ase=load('ips_ase.txt'); %log2 ratios
#       %r=correlate2(abs(log_ase),s)
#       %scatter(log_ase,s)
#       par_counts=load('human_chimp_parental_counts_in_order.txt'); %same order as other files
#       samples=length(par_counts(1,:))/4;
#       for i=1:samples
#       par_counts_all(:,i)=nansum(par_counts(:,(i-1)*4+1:i*4),2);
#       end
#       par=par_counts_all./repmat(nansum(par_counts_all),length(par_counts_all(:,1)),1); 
#       hipsci=(e./repmat(sum(e),length(e(:,1)),1));
#       hipsci_mean=mean(hipsci')';
#                        
#                        r=correlate2(hipsci_mean,par); 
#                        z=(par-repmat(hipsci_mean,1,length(par(1,:))))./repmat(std(hipsci')',1,length(par(1,:)));
#                                                                               [h,p]=ttest2(z(:,1:6)',z(:,7:end)');
#                                                                               
#                                                                               figure,scatter(mean(z(:,1:6)')',mean(z(:,7:end)')')
#                                                                                                   
#                                                                                                  par_diff=log2(nanmedian(par_counts_all(:,1:6),2)./nanmedian(par_counts_all(:,7:12),2));
#                                                                                                   
#                                                                                                   f=find(p<1e-9);
#                                                                                                   for i=1:length(f)
#                                                                                                   figure,histogram(hipsci(f(i),:))
#                                                                                                   hold on
#                                                                                                   histogram(par(f(i),1:6),10)
#                                                                                                   histogram(par(f(i),7:12:end),10)
#                                                                                                   end