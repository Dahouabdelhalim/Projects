GTm <- t(Phenig$G)
GTm[GTm %in% 9]<- NA
colnames(GTm)<- Phenig$SNP
rownames(GTm)<- Phenig$id
GTm_0 <- (GTm==0)+0 # transform GT==0 to 1, rest to 0
GTm_1 <- (GTm==1)+0 # transform GT==1 to 1, rest to 0
GTm_2 <- (GTm==2)+0 # transform GT==2 to 1, rest to 0

colnames(GTm_0) <- paste0(colnames(GTm_0),'_0') # add identifier to each marker name
colnames(GTm_1) <- paste0(colnames(GTm_1),'_1') # add identifier to each marker name
colnames(GTm_2) <- paste0(colnames(GTm_2),'_2') # add identifier to each marker name

# create sambada genotype input

sGTm <- cbind('Ind'=rownames(GTm),GTm_0,GTm_1,GTm_2)

write.table(sGTm, 'c:/Users/jsons/Dropbox/R_projects/Phenig_reanalysis/mol-data-Phenig2.txt', row.names=F, quote=F)

ncol(sGTm)
nrow(sGTm)
