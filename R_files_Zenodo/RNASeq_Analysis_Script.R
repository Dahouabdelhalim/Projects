setwd("~/RNASeq Experiment") #Set working directory

# load packages----
library(plyr)
library(dplyr) 
library(data.table)
library(magrittr)
library(countToFPKM) #version 1.0 https://cran.r-project.org/web/packages/countToFPKM/index.html
library(qdapRegex)
library(IsoformSwitchAnalyzeR)
library(DESeq2)
library(ballgown)
library(pcaExplorer)
library(factoextra)
library(VennDiagram)
library(RRHO)
library(RRHO2)
library(topGO)
library(GOplot)
library(genefilter)
library(devtools)
library(smatr)


########Final Transcript Filtering and Quantification###########################

#Read in count data 
countData<-as.matrix(read.csv("transcript_count_matrix.csv", row.names="transcript_id"))

#Following Bloch et al. 2018, Nature Ecology and Evolution https://www.nature.com/articles/s41559-018-0682-4?WT.feed_name=subjects_genomics....
#Filter out transcripts with <2 FPKM in more than half of samples for each cohort and treatment
#We need to calculate FPKMs per transcript per sample using transcript counts, transcript (feature) lengths, and mean fragment lengths per sample
featureLength<-read.csv("FeatureLengths.csv")
featureLength<-featureLength %>% dplyr::filter(featureLength$isoform_id %in% rownames(countData))
featureLength<- featureLength[order(match(featureLength$isoform_id, rownames(countData))),]
featureLength<-featureLength$length
#Get mean fragment length (insert size) per sample from Picard output
#reset working directory to folder containing insert sizes
setwd("~/RNASeq Experiment/Picard Insert Sizes")
data_all=list.files(pattern='.*.txt')
read<-function(x, mean_insert= "MEAN_INSERT_SIZE") { fread(x, select=mean_insert)}
meanFragmentLength<-lapply(data_all,read)
setwd("~/RNASeq Experiment")
meanFragmentLength<-as.numeric(unlist(meanFragmentLength))
meanFragmentLength<-as.data.frame(cbind(data_all,meanFragmentLength))
colnames(meanFragmentLength)<-c("Sample_ID","Mean_Frag_Length")
#Extract simple sample ID from filenames
meanFragmentLength$Sample_ID<-rm_between(meanFragmentLength$Sample_ID, 'JDSR_RNA_Seq_Libby_RNA_Pool_Schmitz-Libby_', '_I1032_L4_insert_size_metrics.txt', extract=TRUE)
meanFragmentLength$Sample_ID<-gsub("-",".",meanFragmentLength$Sample_ID)
meanFragmentLength$Sample_ID<-gsub("_Head","",meanFragmentLength$Sample_ID)
meanFragmentLength$Sample_ID<-gsub("_Larvae","",meanFragmentLength$Sample_ID)
#Take mean fragment length for each sample ID
meanFragmentLength<-meanFragmentLength %>% dplyr::filter(meanFragmentLength$Sample_ID %in% colnames(countData))
meanFragmentLength<- meanFragmentLength[order(match(meanFragmentLength$Sample_ID, colnames(countData))),]
meanFragmentLength<-meanFragmentLength$Mean_Frag_Length
meanFragmentLength<-as.numeric(as.character(meanFragmentLength))

#To run this code, need to load OLD (v. 1.0) not new version of package, countToFPKM
fpkm_matrix<-fpkm(countData,featureLength,meanFragmentLength) #Save all this as a matrix

#Do the filtering for each treatment within each cohort. The output should be a list of transcripts to 'keep'
FCStringCounts<-as.data.frame(table(c(rownames(fpkm_matrix[which(fpkm_matrix[,1]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,4]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,7]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,10]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,13]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,31]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,34]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,37]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,40]>2),]))))


FEStringCounts<-as.data.frame(table(c(rownames(fpkm_matrix[which(fpkm_matrix[,16]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,19]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,22]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,25]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,28]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,43]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,46]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,49]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,52]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,55]>2),]))))


#Keep only transcripts that are expressed in over half of samples
FCStringCounts<-FCStringCounts[which(FCStringCounts$Freq>4),]
FEStringCounts<-FEStringCounts[which(FEStringCounts$Freq>5),]
#merge the lists together
FStringCounts<-unique(c(as.character(FCStringCounts$Var1),as.character(FEStringCounts$Var1)))


MCStringCounts<-as.data.frame(table(c(rownames(fpkm_matrix[which(fpkm_matrix[,3]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,5]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,9]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,12]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,15]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,33]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,36]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,39]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,42]>2),]))))

MEStringCounts<-as.data.frame(table(c(rownames(fpkm_matrix[which(fpkm_matrix[,18]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,21]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,23]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,27]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,30]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,45]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,48]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,51]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,54]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,57]>2),]))))

MCStringCounts<-MCStringCounts[which(MCStringCounts$Freq>4),]
MEStringCounts<-MEStringCounts[which(MEStringCounts$Freq>5),]
MStringCounts<-unique(c(as.character(MCStringCounts$Var1),as.character(MEStringCounts$Var1)))


LCStringCounts<-as.data.frame(table(c(rownames(fpkm_matrix[which(fpkm_matrix[,2]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,4]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,8]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,11]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,14]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,32]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,35]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,38]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,41]>2),]))))

LEStringCounts<-as.data.frame(table(c(rownames(fpkm_matrix[which(fpkm_matrix[,17]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,20]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,22]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,26]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,29]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,44]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,47]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,50]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,53]>2),]),
                                      rownames(fpkm_matrix[which(fpkm_matrix[,56]>2),]))))

LCStringCounts<-LCStringCounts[which(LCStringCounts$Freq>4),]
LEStringCounts<-LEStringCounts[which(LEStringCounts$Freq>5),]
LStringCounts<-unique(c(as.character(LCStringCounts$Var1),as.character(LEStringCounts$Var1)))


all<-c(as.character(FCStringCounts$Var1),as.character(FEStringCounts$Var1),as.character(MCStringCounts$Var1),as.character(MEStringCounts$Var1), as.character(LCStringCounts$Var1),as.character(LEStringCounts$Var1))
all<-unique(all)
length(all) #16,229 total unique transcripts


#Read in gene annotations
geneannotations<-read.csv("AnnotationDF2.csv")
#Save list of genes expressed in each cohort
geneannotationsF<-geneannotations[geneannotations$isoform_id %in% c(as.character(FCStringCounts$Var1),as.character(FEStringCounts$Var1)),]
geneannotationsM<-geneannotations[geneannotations$isoform_id %in% c(as.character(MCStringCounts$Var1),as.character(MEStringCounts$Var1)),]
geneannotationsL<-geneannotations[geneannotations$isoform_id %in% c(as.character(LCStringCounts$Var1),as.character(LEStringCounts$Var1)),]
#And also take the consensus across cohorts
geneannotations<-geneannotations[geneannotations$isoform_id %in% all,]

#Now take just one column for isoform id and gene/protein product
ALLannotations<-geneannotations[,c("isoform_id","gene_id")]

#Return to count data table and retain only transcripts that passed the filtering criteria
countData<-countData[which(rownames(countData) %in% ALLannotations$isoform_id),]

#Convert transcripts to genes/merge isoforms using IsoformSwitchAnalyzer
countData.g<-isoformToGeneExp(countData,isoformGeneAnnotation=ALLannotations)

# read in study design 
design <- read.csv("studydesign.csv",header=TRUE,row.names="ID")
design$Member<-as.factor(design$Member)
design$Temp<-as.factor(design$Temp)
#Check all sample IDs in design are also in countData.g and match their orders
all(rownames(design) %in% colnames(countData.g))
countData.g <- countData.g[, rownames(design)]
all(rownames(design) == colnames(countData.g)) #Yes



##############Differential Expression Analysis##########################

#Use DESeq to normalize counts and perform differential expression
dds.full <- DESeqDataSetFromMatrix(countData = countData.g, colData = design, design = ~ Member+Temp+Social+Member*Temp*Social)
#Run PCA with prcomp
vst.full<-vst(dds.full) #apply variance stabilizing transformation
rv <- rowVars(assay(vst.full)) #select row variables
select <- order(rv, decreasing = TRUE)[seq_len(min(12406, length(rv)))] #all transcripts
mat <- t( assay(vst.full)[select, ] ) #generate input matrix
pca.full<-prcomp(mat) #perform the pca
summary(pca.full)
pca.full$sdev^2 #eigenvalues
loadings<-as.data.frame(pca.full$x[,1:2])
loadings$Sample<-rownames(design)
loadings<-cbind(loadings,design[,4:6])
anova(lm(PC1~Member*Temp*Social,data=loadings))
anova(lm(PC2~Member*Temp*Social,data=loadings))
Tab1<-as.data.frame(anova(lm(PC1~Member*Temp*Social,data=loadings)))
Tab2<-as.data.frame(anova(lm(PC2~Member*Temp*Social,data=loadings)))

#Perform DE separately for each cohort (starting with females) and execute contrasts of interest
dds.F <- DESeqDataSetFromMatrix(countData = countData.g[,c(1:19)], colData = design[1:19,], design = ~ Temp+Social)
#Run likelihood-ratio test to identify DEGs
dds.Ftemp <- DESeq(dds.F,test='LRT',full=~Temp+Social,reduced=~Social)
dds.Fsocial <- DESeq(dds.F,test='LRT',full=~Temp+Social,reduced=~Temp)
#Extract results
resFTemp<-results(dds.Ftemp,name="Temp_Hot_vs_Cold")
resFSocial<-results(dds.Fsocial,name="Social_Exp_vs_Control")
#Save significant DEGs (adjusted p<0.05) as a data frame
DEGs.tempF<-as.data.frame(cbind(rownames(subset(resFTemp,padj<0.05)),subset(resFTemp,padj<0.05)$log2FoldChange,subset(resFTemp,padj<0.05)$padj)) 
colnames(DEGs.tempF)<-c("gene_id","log2FoldChange","p")
#Merge in annotations for DEGs
DEGs.tempF<-merge(DEGs.tempF,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.tempF<-DEGs.tempF[!duplicated(DEGs.tempF$gene_id),]
DEGs.socialF<-as.data.frame(cbind(rownames(subset(resFSocial,padj<0.05)),subset(resFSocial,padj<0.05)$log2FoldChange,subset(resFSocial,padj<0.05)$padj)) 
colnames(DEGs.socialF)<-c("gene_id","log2FoldChange","p")
DEGs.socialF<-merge(DEGs.socialF,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.socialF<-DEGs.socialF[!duplicated(DEGs.socialF$gene_id),]
#Save these data frames for supplement
write.csv(DEGs.tempF,file="DEGS-tempF.csv")
write.csv(DEGs.socialF,file="DEGS-socialF.csv")

#Repeat the above for males and larvae
dds.M <- DESeqDataSetFromMatrix(countData = countData.g[,c(20:38)], colData = design[20:38,], design = ~ Temp+Social)
dds.Mtemp <- DESeq(dds.M,test='LRT',full=~Temp+Social,reduced=~Social)
dds.Msocial <- DESeq(dds.M,test='LRT',full=~Temp+Social,reduced=~Temp)
resMTemp<-results(dds.Mtemp,name="Temp_Hot_vs_Cold")
resMSocial<-results(dds.Msocial,name="Social_Exp_vs_Control")
DEGs.tempM<-as.data.frame(cbind(rownames(subset(resMTemp,padj<0.05)),subset(resMTemp,padj<0.05)$log2FoldChange,subset(resMTemp,padj<0.05)$padj)) 
colnames(DEGs.tempM)<-c("gene_id","log2FoldChange","p")
DEGs.tempM<-merge(DEGs.tempM,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.tempM<-DEGs.tempM[!duplicated(DEGs.tempM$gene_id),]
DEGs.socialM<-as.data.frame(cbind(rownames(subset(resMSocial,padj<0.05)),subset(resMSocial,padj<0.05)$log2FoldChange,subset(resMSocial,padj<0.05)$padj)) 
colnames(DEGs.socialM)<-c("gene_id","log2FoldChange","p")
DEGs.socialM<-merge(DEGs.socialM,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.socialM<-DEGs.socialM[!duplicated(DEGs.socialM$gene_id),]
write.csv(DEGs.tempM,file="DEGS-tempM.csv")
write.csv(DEGs.socialM,file="DEGS-socialM.csv")

dds.L <- DESeqDataSetFromMatrix(countData = countData.g[,c(39:57)], colData = design[39:57,], design = ~ Temp+Social)
dds.Ltemp <- DESeq(dds.L,test='LRT',full=~Temp+Social,reduced=~Social)
dds.Lsocial <- DESeq(dds.L,test='LRT',full=~Temp+Social,reduced=~Temp)
resLTemp<-results(dds.Ltemp,name="Temp_Hot_vs_Cold")
resLSocial<-results(dds.Lsocial,name="Social_Exp_vs_Control")
DEGs.tempL<-as.data.frame(cbind(rownames(subset(resLTemp,padj<0.05)),subset(resLTemp,padj<0.05)$log2FoldChange,subset(resLTemp,padj<0.05)$padj)) 
colnames(DEGs.tempL)<-c("gene_id","log2FoldChange","p")
DEGs.tempL<-merge(DEGs.tempL,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.tempL<-DEGs.tempL[!duplicated(DEGs.tempL$gene_id),]
DEGs.socialL<-as.data.frame(cbind(rownames(subset(resLSocial,padj<0.05)),subset(resLSocial,padj<0.05)$log2FoldChange,subset(resLSocial,padj<0.05)$padj)) 
colnames(DEGs.socialL)<-c("gene_id","log2FoldChange","p")
DEGs.socialL<-merge(DEGs.socialL,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.socialL<-DEGs.socialL[!duplicated(DEGs.socialL$gene_id),]
write.csv(DEGs.tempL,file="DEGS-tempL.csv")
write.csv(DEGs.socialL,file="DEGS-socialL.csv")




#######Overlap Between Family Members For Main Effects########################

#Look at overlap between family members
Socialvenn<-venn.diagram(list(DEGs.socialL$gene_id,DEGs.socialF$gene_id, DEGs.socialM$gene_id ), filename=NULL,fill=c("goldenrod1","red","dodgerblue1"),category.names=c("","",""),cat.just=list(c(1,1.4) , c(0.6,0.6) ,  c(0.6,0.6)),cat.cex=1.3,cex=1.3 ,fontfamily ="Arial" ,cat.fontfamily="Arial",margin=0.1)
grid.draw(Socialvenn)
dev.off()
Tempvenn<-venn.diagram(list(DEGs.tempL$gene_id,DEGs.tempF$gene_id, DEGs.tempM$gene_id ), filename=NULL,fill=c("goldenrod1","red","dodgerblue1"),category.names=c("","",""),cat.just=list(c(1,1.4) , c(0.6,0.6) ,  c(0.6,0.6)),cat.cex=1.3,cex=1.3 ,fontfamily ="Arial" ,cat.fontfamily="Arial",margin=0.1)
grid.draw(Tempvenn)
dev.off()


#Test significance of overlap with permutation
#Build a bootstrap function to do this iteratively.
perm.fun<-function(dat1,dat2,n1,n2,iterations) {
  output <- vector("double", 1000)  
  for (i in 1:iterations) {                     
  samp1<-sample(dat1$gene_id,n1)
  samp2<-sample(dat2$gene_id,n2)
  output[[i]] <- length(intersect(samp1,samp2))    
  }
  return(output)
}

null.socialMF<-perm.fun(geneannotationsM,geneannotationsF,nrow(DEGs.socialM),nrow(DEGs.socialF),10000)
#calculate proportion of null distribution that exceeds observed
length(which(null.socialMF>3)+1)/10001
null.socialFL<-perm.fun(geneannotationsF,geneannotationsL,nrow(DEGs.socialF),nrow(DEGs.socialL),10000)
length(which(null.socialFL>115)+1)/10001
null.socialML<-perm.fun(geneannotationsM,geneannotationsL,nrow(DEGs.socialM),nrow(DEGs.socialL),10000)
length(which(null.socialML>28)+1)/10001
null.tempMF<-perm.fun(geneannotationsM,geneannotationsF,nrow(DEGs.tempM),nrow(DEGs.tempF),10000)
length(which(null.tempMF>60)+1)/10001
null.tempFL<-perm.fun(geneannotationsF,geneannotationsL,nrow(DEGs.tempF),nrow(DEGs.tempL),10000)
length(which(null.tempFL>31)+1)/10001
null.tempML<-perm.fun(geneannotationsM,geneannotationsL,nrow(DEGs.tempM),nrow(DEGs.tempL),10000)
length(which(null.tempML>31)+1)/10001
#Overlap across all groups greater than expected by chance


#Examine overlap using rank-rank hypogeometric overlap (RRHO) Analysis
#Generate a list of gene names with a ranking based on -log10(adj_pval)*sign(logfoldchange)
temp.mums<-as.data.frame(cbind(rownames(resFTemp),-log10(resFTemp$pvalue)*sign(resFTemp$log2FoldChange)) )
temp.dads<-as.data.frame(cbind(rownames(resMTemp),-log10(resMTemp$pvalue)*sign(resMTemp$log2FoldChange)) )
temp.larv<-as.data.frame(cbind(rownames(resLTemp),-log10(resLTemp$pvalue)*sign(resLTemp$log2FoldChange)) )
temp.mums<-temp.mums[order(temp.mums$V1),]
temp.dads<-temp.dads[order(temp.dads$V1),]
temp.larv<-temp.larv[order(temp.larv$V1),]
temp.mums[which(is.na(temp.mums$V2)),]$V2=0 #DESeq sets p=NA for genes with zero expression. Change NAs to 0
temp.dads[which(is.na(temp.dads$V2)),]$V2=0
temp.larv[which(is.na(temp.larv$V2)),]$V2=0
temp.dads$V2<-as.numeric(as.character(temp.dads$V2))
temp.mums$V2<-as.numeric(as.character(temp.mums$V2))
temp.larv$V2<-as.numeric(as.character(temp.larv$V2))
#Run analysis
rrho_mftemp<-RRHO2_initialize(temp.mums,temp.dads)
par(mar=c(1,1,2,3))
RRHO2_heatmap((rrho_mftemp))
dev.off()
rrho_fltemp<-RRHO2_initialize(temp.mums,temp.larv)
RRHO2_heatmap((rrho_fltemp))
dev.off()
rrho_mltemp<-RRHO2_initialize(temp.dads,temp.larv)
RRHO2_heatmap((rrho_mltemp))
#Repeat for social
social.mums<-as.data.frame(cbind(rownames(resFSocial),-log10(resFSocial$pvalue)*sign(resFSocial$log2FoldChange)) )
social.dads<-as.data.frame(cbind(rownames(resMSocial),-log10(resMSocial$pvalue)*sign(resMSocial$log2FoldChange)) )
social.larv<-as.data.frame(cbind(rownames(resLSocial),-log10(resLSocial$pvalue)*sign(resLSocial$log2FoldChange)) )
social.mums<-social.mums[order(social.mums$V1),]
social.dads<-social.dads[order(social.dads$V1),]
social.larv<-social.larv[order(social.larv$V1),]
social.mums[which(is.na(social.mums$V2)),]$V2=0 #DESeq sets p=NA for genes with zero expression. Change NAs to 0
social.dads[which(is.na(social.dads$V2)),]$V2=0
social.larv[which(is.na(social.larv$V2)),]$V2=0
social.dads$V2<-as.numeric(as.character(social.dads$V2))
social.mums$V2<-as.numeric(as.character(social.mums$V2))
social.larv$V2<-as.numeric(as.character(social.larv$V2))

rrho_mfsoc<-RRHO2_initialize(social.mums,social.dads)
RRHO2_heatmap((rrho_mfsoc))
dev.off()
rrho_flsoc<-RRHO2_initialize(social.mums,social.larv)
RRHO2_heatmap((rrho_flsoc))
dev.off()
rrho_mlsoc<-RRHO2_initialize(social.dads,social.larv)
RRHO2_heatmap((rrho_mlsoc))
dev.off()



##########Interaction (Temp x Parenting) Effects#############################

#Examine differential expression in response to the combination of thermal stress and parenting
#First use the LRT approach including a statistical interation term
dds.Ffull <- DESeqDataSetFromMatrix(countData = countData.g[,c(1:19)], colData = design[1:19,], design = ~ Temp+Social+Temp:Social)
dds.Finteraction<-DESeq(dds.Ffull,test='LRT',full=~Temp+Social+Temp:Social,reduced=~Temp+Social)
resFinteraction<-results(dds.Finteraction,name="TempHot.SocialExp")
DEGs.tempxsocialF<-as.data.frame(cbind(rownames(subset(resFinteraction,padj<0.05)),subset(resFinteraction,padj<0.05)$log2FoldChange,subset(resFinteraction,padj<0.05)$padj)) 
colnames(DEGs.tempxsocialF)<-c("gene_id","log2FoldChange","p")
DEGs.tempxsocialF<-merge(DEGs.tempxsocialF,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.tempxsocialF<-DEGs.tempxsocialF[!duplicated(DEGs.tempxsocialF$gene_id),]
write.csv(DEGs.tempxsocialF,file="DEGS-tempsocialF.csv")

dds.Mfull <- DESeqDataSetFromMatrix(countData = countData.g[,c(20:38)], colData = design[20:38,], design = ~ Temp+Social+Temp:Social)
dds.Minteraction<-DESeq(dds.Mfull,test='LRT',full=~Temp+Social+Temp:Social,reduced=~Temp+Social)
resMinteraction<-results(dds.Minteraction,name="TempHot.SocialExp")
DEGs.tempxsocialM<-as.data.frame(cbind(rownames(subset(resMinteraction,padj<0.05)),subset(resMinteraction,padj<0.05)$log2FoldChange,subset(resMinteraction,padj<0.05)$padj)) 
colnames(DEGs.tempxsocialM)<-c("gene_id","log2FoldChange","p")
DEGs.tempxsocialM<-merge(DEGs.tempxsocialM,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.tempxsocialM<-DEGs.tempxsocialM[!duplicated(DEGs.tempxsocialM$gene_id),]
write.csv(DEGs.tempxsocialM,file="DEGS-tempsocialM.csv")

dds.Lfull <- DESeqDataSetFromMatrix(countData = countData.g[,c(39:57)], colData = design[39:57,], design = ~ Temp+Social+Temp:Social)
dds.Linteraction<-DESeq(dds.Lfull,test='LRT',full=~Temp+Social+Temp:Social,reduced=~Temp+Social)
resLinteraction<-results(dds.Linteraction,name="TempHot.SocialExp")
DEGs.tempxsocialL<-as.data.frame(cbind(rownames(subset(resLinteraction,padj<0.05)),subset(resLinteraction,padj<0.05)$log2FoldChange,subset(resLinteraction,padj<0.05)$padj)) 
colnames(DEGs.tempxsocialL)<-c("gene_id","log2FoldChange","p")
DEGs.tempxsocialL<-merge(DEGs.tempxsocialL,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
DEGs.tempxsocialL<-DEGs.tempxsocialL[!duplicated(DEGs.tempxsocialL$gene_id),]
write.csv(DEGs.tempxsocialL,file="DEGS-tempsocialL.csv")

#Only females show a large number of 'buffering' genes (n=79). Let's examine the patterns of differential expression of these genes across conditions
#Create a dataframe of normalized expression of the 189 parenting genes
vst.F<-vst(dds.F)
PlasticityFemales<-as.data.frame(t(assay(vst.F[rownames(assay(vst.F)) %in% DEGs.tempxsocialF$gene_id,])))
PlasticityFemales<-cbind(design[c(1:19),c(3,5,6)],PlasticityFemales)
#Summarize normalized expression for each gene for each group
PlasticityFemales$SocialTemp<-apply(PlasticityFemales[,3:2], 1 , paste , collapse = "-")
PlasticityFemales_summary<- as.data.frame(t(PlasticityFemales[,c(4:83)] %>%
                                              group_by(SocialTemp) %>%
                                              summarise_all("mean")))[2:80,]
colnames(PlasticityFemales_summary)<-c("Control-Cold","Control-Hot","Exp-Cold","Exp-Hot")

#Build forloop to extract the run a posthoc test to compare levels of expression between each pair of sample sets for each gene in the plasticity set
posthoc <- lapply(4:82, function(x) TukeyHSD(aov(PlasticityFemales[,x] ~ PlasticityFemales$SocialTemp)))
posthoc_ControlHotvControlCold<-t(as.data.frame(lapply(posthoc, function(x) x[[1]][1,c(1,4)] )))
rownames(posthoc_ControlHotvControlCold)<-colnames(PlasticityFemales[,c(4:82)])
colnames(posthoc_ControlHotvControlCold)<-c("ControlHotvControlCold_diff","ControlHotvControlCold_pval")
posthoc_ExpColdvControlCold<-t(as.data.frame(lapply(posthoc, function(x) x[[1]][2,c(1,4)] )))
rownames(posthoc_ExpColdvControlCold)<-colnames(PlasticityFemales[,c(4:82)])
colnames(posthoc_ExpColdvControlCold)<-c("ExpColdvControlCold_diff","ExpColdvControlCold_pval")
posthoc_ExpHotvControlCold<-t(as.data.frame(lapply(posthoc, function(x) x[[1]][3,c(1,4)] )))
rownames(posthoc_ExpHotvControlCold)<-colnames(PlasticityFemales[,c(4:82)])
colnames(posthoc_ExpHotvControlCold)<-c("ExpHotvControlCold_diff","ExpHotvControlCold_pval")
posthoc_ExpColdvControlHot<-t(as.data.frame(lapply(posthoc, function(x) x[[1]][4,c(1,4)] )))
rownames(posthoc_ExpColdvControlHot)<-colnames(PlasticityFemales[,c(4:82)])
colnames(posthoc_ExpColdvControlHot)<-c("ExpColdvControlHot_diff","ExpColdvControlHot_pval")
posthoc_ExpHotvControlHot<-t(as.data.frame(lapply(posthoc, function(x) x[[1]][5,c(1,4)] )))
rownames(posthoc_ExpHotvControlHot)<-colnames(PlasticityFemales[,c(4:82)])
colnames(posthoc_ExpHotvControlHot)<-c("ExpHotvControlHot_diff","ExpHotvControlHot_pval")
posthoc_ExpHotvExpCold<-t(as.data.frame(lapply(posthoc, function(x) x[[1]][6,c(1,4)] )))
rownames(posthoc_ExpHotvExpCold)<-colnames(PlasticityFemales[,c(4:82)])
colnames(posthoc_ExpHotvExpCold)<-c("ExpHotvExpCold_diff","ExpHotvExpCold_pval")
posthoc_df<-cbind(posthoc_ControlHotvControlCold,posthoc_ExpColdvControlCold,posthoc_ExpHotvControlCold,posthoc_ExpColdvControlHot,posthoc_ExpHotvControlHot,posthoc_ExpHotvExpCold)

posthoc_df<-cbind(PlasticityFemales_summary,posthoc_df)
posthoc_df$TrendMagnitudeBefore<-NA
posthoc_df[which(posthoc_df$ControlHotvControlCold_pval>0.05),]$TrendMagnitudeBefore<-"Same"
posthoc_df[which(posthoc_df$ControlHotvControlCold_pval<0.05&posthoc_df$ControlHotvControlCold_diff>0),]$TrendMagnitudeBefore<-"Higher at 24 C"
posthoc_df[which(posthoc_df$ControlHotvControlCold_pval<0.05&posthoc_df$ControlHotvControlCold_diff<0),]$TrendMagnitudeBefore<-"Lower at 24 C"
posthoc_df$TrendMagnitudeAfter<-NA
posthoc_df[which(posthoc_df$ExpHotvExpCold_pval>0.05),]$TrendMagnitudeAfter<-"Same"
posthoc_df[which(posthoc_df$ExpHotvExpCold_pval<0.05&posthoc_df$ExpHotvExpCold_diff>0),]$TrendMagnitudeAfter<-"Higher at 24 C"
posthoc_df[which(posthoc_df$ExpHotvExpCold_pval<0.05&posthoc_df$ExpHotvExpCold_diff<0),]$TrendMagnitudeAfter<-"Lower at 24 C"
posthoc_df$TrendSign20<-NA
posthoc_df[which(posthoc_df$ExpColdvControlCold_pval>0.05),]$TrendSign20<-"No Change"
posthoc_df[which(posthoc_df$ExpColdvControlCold_pval<0.05&posthoc_df$ExpColdvControlCold_diff>0),]$TrendSign20<-"Increases After Parenting"
posthoc_df[which(posthoc_df$ExpColdvControlCold_pval<0.05&posthoc_df$ExpColdvControlCold_diff<0),]$TrendSign20<-"Decreases After Parenting"
posthoc_df$TrendSign24<-NA
posthoc_df[which(posthoc_df$ExpHotvControlHot_pval>0.05),]$TrendSign24<-"No Change"
posthoc_df[which(posthoc_df$ExpHotvControlHot_pval<0.05&posthoc_df$ExpHotvControlHot_diff>0),]$TrendSign24<-"Increases After Parenting"
posthoc_df[which(posthoc_df$ExpHotvControlHot_pval<0.05&posthoc_df$ExpHotvControlHot_diff<0),]$TrendSign24<-"Decreases After Parenting"

write.csv(posthoc_df,file="Supp Table 11.csv")


#As a complementary approach, compare differential expression of parenting genes between the two thermal environments (i.e., Log2fold changes at 20*C vs 24*C)
#Need to redo the DE Analysis on subsets of females
dds.socialF_cold <- DESeqDataSetFromMatrix(countData = countData.g[,c(1:10)], 
                                           colData = design[c(1:10),], design = ~ Social)
dds.socialF_cold <- DESeq(dds.socialF_cold)
resFsocial_cold<-lfcShrink(dds.socialF_cold, coef="Social_Exp_vs_Control", type='ashr') 
Logfold_socialF_cold<-as.data.frame(cbind(rownames(resFsocial_cold),resFsocial_cold$log2FoldChange,resFsocial_cold$padj)) 
colnames(Logfold_socialF_cold)<-c("gene_id","log2FoldChange_cold","p")
Logfold_socialF_cold<-merge(Logfold_socialF_cold,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
Logfold_socialF_cold<-Logfold_socialF_cold[!duplicated(Logfold_socialF_cold$gene_id),]
dds.socialF_hot <- DESeqDataSetFromMatrix(countData = countData.g[,c(11:19)], 
                                          colData = design[c(11:19),], design = ~ Social)
dds.socialF_hot <- DESeq(dds.socialF_hot)
resFsocial_hot<-lfcShrink(dds.socialF_hot, coef="Social_Exp_vs_Control", type='ashr') 
Logfold_socialF_hot<-as.data.frame(cbind(rownames(resFsocial_hot),resFsocial_hot$log2FoldChange,resFsocial_hot$padj)) 
colnames(Logfold_socialF_hot)<-c("gene_id","log2FoldChange_hot","p")
Logfold_socialF_hot<-merge(Logfold_socialF_hot,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
Logfold_socialF_hot<-Logfold_socialF_hot[!duplicated(Logfold_socialF_hot$gene_id),]

Logfold_socialF<-merge(Logfold_socialF_cold[,c(1,2)],Logfold_socialF_hot[,c(1,2)],by="gene_id")
Logfold_socialF<-Logfold_socialF[which(Logfold_socialF$gene_id %in% DEGs.socialF$gene_id),]
Logfold_socialF$log2FoldChange_hot<-as.numeric(as.character(Logfold_socialF$log2FoldChange_hot))
Logfold_socialF$log2FoldChange_cold<-as.numeric(as.character(Logfold_socialF$log2FoldChange_cold))

#Perform major axis regression to estimate slope
sma(log2FoldChange_hot~log2FoldChange_cold, data=Logfold_socialF,method='MA',slope.test=1)


#Males
dds.socialM_cold <- DESeqDataSetFromMatrix(countData = countData.g[,c(20:29)], 
                                           colData = design[c(20:29),], design = ~ Social)
dds.socialM_cold <- DESeq(dds.socialM_cold)
resMsocial_cold<-lfcShrink(dds.socialM_cold, coef="Social_Exp_vs_Control", type='ashr') 
Logfold_socialM_cold<-as.data.frame(cbind(rownames(resMsocial_cold),resMsocial_cold$log2FoldChange,resMsocial_cold$padj)) 
colnames(Logfold_socialM_cold)<-c("gene_id","log2FoldChange_cold","p")
Logfold_socialM_cold<-merge(Logfold_socialM_cold,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
Logfold_socialM_cold<-Logfold_socialM_cold[!duplicated(Logfold_socialM_cold$gene_id),]
dds.socialM_hot <- DESeqDataSetFromMatrix(countData = countData.g[,c(30:38)], 
                                          colData = design[c(30:38),], design = ~ Social)
dds.socialM_hot <- DESeq(dds.socialM_hot)
resMsocial_hot<-lfcShrink(dds.socialM_hot, coef="Social_Exp_vs_Control", type='ashr') 
Logfold_socialM_hot<-as.data.frame(cbind(rownames(resMsocial_hot),resMsocial_hot$log2FoldChange,resMsocial_hot$padj)) 
colnames(Logfold_socialM_hot)<-c("gene_id","log2FoldChange_hot","p")
Logfold_socialM_hot<-merge(Logfold_socialM_hot,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
Logfold_socialM_hot<-Logfold_socialM_hot[!duplicated(Logfold_socialM_hot$gene_id),]

Logfold_socialM<-merge(Logfold_socialM_cold[,c(1,2)],Logfold_socialM_hot[,c(1,2)],by="gene_id")
Logfold_socialM<-Logfold_socialM[which(Logfold_socialM$gene_id %in% DEGs.socialM$gene_id),]
Logfold_socialM$log2FoldChange_hot<-as.numeric(as.character(Logfold_socialM$log2FoldChange_hot))
Logfold_socialM$log2FoldChange_cold<-as.numeric(as.character(Logfold_socialM$log2FoldChange_cold))

sma(log2FoldChange_hot~log2FoldChange_cold, data=Logfold_socialM,method='MA',slope.test=1)


#Larvae
dds.socialL_cold <- DESeqDataSetFromMatrix(countData = countData.g[,c(39:48)], 
                                           colData = design[c(39:48),], design = ~ Social)
dds.socialL_cold <- DESeq(dds.socialL_cold)
resLsocial_cold<-lfcShrink(dds.socialL_cold, coef="Social_Exp_vs_Control", type='ashr') 
Logfold_socialL_cold<-as.data.frame(cbind(rownames(resLsocial_cold),resLsocial_cold$log2FoldChange,resLsocial_cold$padj)) 
colnames(Logfold_socialL_cold)<-c("gene_id","log2FoldChange_cold","p")
Logfold_socialL_cold<-merge(Logfold_socialL_cold,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
Logfold_socialL_cold<-Logfold_socialL_cold[!duplicated(Logfold_socialL_cold$gene_id),]
dds.socialL_hot <- DESeqDataSetFromMatrix(countData = countData.g[,c(49:57)], 
                                          colData = design[c(49:57),], design = ~ Social)
dds.socialL_hot <- DESeq(dds.socialL_hot)
resLsocial_hot<-lfcShrink(dds.socialL_hot, coef="Social_Exp_vs_Control", type='ashr') 
Logfold_socialL_hot<-as.data.frame(cbind(rownames(resLsocial_hot),resLsocial_hot$log2FoldChange,resLsocial_hot$padj)) 
colnames(Logfold_socialL_hot)<-c("gene_id","log2FoldChange_hot","p")
Logfold_socialL_hot<-merge(Logfold_socialL_hot,geneannotations[,c(1,3:8)],by="gene_id",all.x=T)
Logfold_socialL_hot<-Logfold_socialL_hot[!duplicated(Logfold_socialL_hot$gene_id),]

Logfold_socialL<-merge(Logfold_socialL_cold[,c(1,2)],Logfold_socialL_hot[,c(1,2)],by="gene_id")
Logfold_socialL<-Logfold_socialL[which(Logfold_socialL$gene_id %in% DEGs.socialL$gene_id),]
Logfold_socialL$log2FoldChange_hot<-as.numeric(as.character(Logfold_socialL$log2FoldChange_hot))
Logfold_socialL$log2FoldChange_cold<-as.numeric(as.character(Logfold_socialL$log2FoldChange_cold))

sma(log2FoldChange_hot~log2FoldChange_cold, data=Logfold_socialL,method='MA',slope.test=1)




#############Functional Enrichment Analysis################################

#Use TopGO to set up input list of genes and their GO terms
GOList<-geneannotations[,c(1,28)]
GOList<-GOList[-which(is.na(GOList$GO.terms)),]
GOList<-GOList[!duplicated(GOList$gene_id),]
GOList$GO.terms<-gsub(";",",",GOList$GO.terms)
write.table(GOList, file = "GOList.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
geneID2GO<-readMappings(file="GOList.txt")
geneNames <- names(geneID2GO)

#Perform functional enrichment analysis of all genes differentially expressed 
#in response to temp (across all groups)
options(scipen = 999)
temp.dad<-DEGs.tempM$gene_id 
temp.mum<-DEGs.tempF$gene_id
temp.larv<-DEGs.tempL$gene_id
GOList.tempD<-GOList[GOList$gene_id %in% temp.dad,]
GOList.tempM<-GOList[GOList$gene_id %in% temp.mum,]
GOList.tempL<-GOList[GOList$gene_id %in% temp.larv,]
write.table(GOList.tempD, file = "GOList.tempD.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
write.table(GOList.tempM, file = "GOList.tempM.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
write.table(GOList.tempL, file = "GOList.tempL.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
tempD.GOs<-readMappings(file="GOList.tempD.txt")
tempM.GOs<-readMappings(file="GOList.tempM.txt")
tempL.GOs<-readMappings(file="GOList.tempL.txt")
tempDNames<-names(tempD.GOs)
tempMNames<-names(tempM.GOs)
tempLNames<-names(tempL.GOs)
geneList_tempD <- factor(as.integer(geneNames %in% tempDNames))
geneList_tempM <- factor(as.integer(geneNames %in% tempMNames))
geneList_tempL <- factor(as.integer(geneNames %in% tempLNames))
names(geneList_tempD) <- geneNames
names(geneList_tempM) <- geneNames
names(geneList_tempL) <- geneNames
GOdata_MF_tempD <- new("topGOdata", ontology = "MF", allGenes = geneList_tempD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF_tempM <- new("topGOdata", ontology = "MF", allGenes = geneList_tempM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF_tempL <- new("topGOdata", ontology = "MF", allGenes = geneList_tempL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_tempD <- new("topGOdata", ontology = "BP", allGenes = geneList_tempD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_tempM <- new("topGOdata", ontology = "BP", allGenes = geneList_tempM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_tempL <- new("topGOdata", ontology = "BP", allGenes = geneList_tempL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_tempD <- new("topGOdata", ontology = "CC", allGenes = geneList_tempD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_tempM <- new("topGOdata", ontology = "CC", allGenes = geneList_tempM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_tempL <- new("topGOdata", ontology = "CC", allGenes = geneList_tempL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
BPres_tempD<-runTest(GOdata_BP_tempD, algorithm = "classic", statistic = "fisher")
BPres_tempM<-runTest(GOdata_BP_tempM, algorithm = "classic", statistic = "fisher")
BPres_tempL<-runTest(GOdata_BP_tempL, algorithm = "classic", statistic = "fisher")
#The below will take the number of enriched GO terms that are significant (p<0.01) from each category. To see how many there are,return the 'BPres_tempD' object created above.
#Note: Weight = P VALUE!!! Smaller = greater weight!!
BP_obj_Test_res_tempD <- GenTable( GOdata_BP_tempD,weight   = BPres_tempD,orderBy  = "weight",topNodes=50,numChar=1000)
BP_obj_Test_res_tempM <- GenTable( GOdata_BP_tempM,weight   = BPres_tempM,orderBy  = "weight",topNodes=48,numChar=1000)
BP_obj_Test_res_tempL <- GenTable( GOdata_BP_tempL,weight   = BPres_tempL,orderBy  = "weight",topNodes=54,numChar=1000)
MFres_tempD<-runTest(GOdata_MF_tempD, algorithm = "classic", statistic = "fisher")
MFres_tempM<-runTest(GOdata_MF_tempM, algorithm = "classic", statistic = "fisher")
MFres_tempL<-runTest(GOdata_MF_tempL, algorithm = "classic", statistic = "fisher")
MF_obj_Test_res_tempD <- GenTable( GOdata_MF_tempD, weight   = MFres_tempD,orderBy  = "weight",topNodes=22,numChar=1000)
MF_obj_Test_res_tempM <- GenTable( GOdata_MF_tempM, weight   = MFres_tempM,orderBy  = "weight",topNodes=28,numChar=1000)
MF_obj_Test_res_tempL <- GenTable( GOdata_MF_tempL, weight   = MFres_tempL,orderBy  = "weight",topNodes=20,numChar=1000)
CCres_tempD<-runTest(GOdata_CC_tempD, algorithm = "classic", statistic = "fisher")
CCres_tempM<-runTest(GOdata_CC_tempM, algorithm = "classic", statistic = "fisher")
CCres_tempL<-runTest(GOdata_CC_tempL, algorithm = "classic", statistic = "fisher")
CC_obj_Test_res_tempD <- GenTable( GOdata_CC_tempD,weight   = CCres_tempD,orderBy  = "weight",topNodes=14,numChar=1000)
CC_obj_Test_res_tempM <- GenTable( GOdata_CC_tempM,weight   = CCres_tempM,orderBy  = "weight",topNodes=8,numChar=1000)
CC_obj_Test_res_tempL <- GenTable( GOdata_CC_tempL,weight   = CCres_tempL,orderBy  = "weight",topNodes=19,numChar=1000)
BP_obj_Test_res_tempD$category<-"BP"
MF_obj_Test_res_tempD$category<-"MF"
CC_obj_Test_res_tempD$category<-"CC"
BP_obj_Test_res_tempM$category<-"BP"
MF_obj_Test_res_tempM$category<-"MF"
CC_obj_Test_res_tempM$category<-"CC"
BP_obj_Test_res_tempL$category<-"BP"
MF_obj_Test_res_tempL$category<-"MF"
CC_obj_Test_res_tempL$category<-"CC"
allGORes_tempD <- data.frame(rbind(BP_obj_Test_res_tempD, CC_obj_Test_res_tempD, MF_obj_Test_res_tempD))
allGORes_tempM <- data.frame(rbind(BP_obj_Test_res_tempM, CC_obj_Test_res_tempM, MF_obj_Test_res_tempM))
allGORes_tempL <- data.frame(rbind(BP_obj_Test_res_tempL, CC_obj_Test_res_tempL, MF_obj_Test_res_tempL))
allGORes_tempD$weight <- as.numeric(allGORes_tempD$weight)
allGORes_tempM$weight <- as.numeric(allGORes_tempM$weight)
allGORes_tempL$weight <- as.numeric(allGORes_tempL$weight)
allGORes_tempD <- allGORes_tempD[order(allGORes_tempD$weight), ]
allGORes_tempM <- allGORes_tempM[order(allGORes_tempM$weight), ]
allGORes_tempL <- allGORes_tempL[order(allGORes_tempL$weight), ]
write.csv(allGORes_tempD, file = "allGORes_tempM.csv")
write.csv(allGORes_tempM, file = "allGORes_tempF.csv")
write.csv(allGORes_tempL, file = "allGORes_tempL.csv")


#Functional Enrichment of DEGs in response to parenting/being parented
options(scipen = 999)
social.mum<-DEGs.socialF$gene_id
social.dad<-DEGs.socialM$gene_id
social.larv<-DEGs.socialL$gene_id
GOList.socialM<-GOList[GOList$gene_id %in% social.mum,]
GOList.socialD<-GOList[GOList$gene_id %in% social.dad,]
GOList.socialL<-GOList[GOList$gene_id %in% social.larv,]
write.table(GOList.socialM, file = "GOList.socialM.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
write.table(GOList.socialD, file = "GOList.socialD.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
write.table(GOList.socialL, file = "GOList.socialL.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
socialM.GOs<-readMappings(file="GOList.socialM.txt")
socialD.GOs<-readMappings(file="GOList.socialD.txt")
socialL.GOs<-readMappings(file="GOList.socialL.txt")
socialMNames<-names(socialM.GOs)
socialDNames<-names(socialD.GOs)
socialLNames<-names(socialL.GOs)
geneList_socialM <- factor(as.integer(geneNames %in% socialMNames))
geneList_socialD <- factor(as.integer(geneNames %in% socialDNames))
geneList_socialL <- factor(as.integer(geneNames %in% socialLNames))
names(geneList_socialM) <- geneNames
names(geneList_socialD) <- geneNames
names(geneList_socialL) <- geneNames
GOdata_MF_socialM <- new("topGOdata", ontology = "MF", allGenes = geneList_socialM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF_socialD <- new("topGOdata", ontology = "MF", allGenes = geneList_socialD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF_socialL <- new("topGOdata", ontology = "MF", allGenes = geneList_socialL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_socialM <- new("topGOdata", ontology = "BP", allGenes = geneList_socialM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_socialD <- new("topGOdata", ontology = "BP", allGenes = geneList_socialD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_socialL <- new("topGOdata", ontology = "BP", allGenes = geneList_socialL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_socialM <- new("topGOdata", ontology = "CC", allGenes = geneList_socialM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_socialD <- new("topGOdata", ontology = "CC", allGenes = geneList_socialD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_socialL <- new("topGOdata", ontology = "CC", allGenes = geneList_socialL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
BPres_socialM<-runTest(GOdata_BP_socialM, algorithm = "classic", statistic = "fisher")
BPres_socialD<-runTest(GOdata_BP_socialD, algorithm = "classic", statistic = "fisher")
BPres_socialL<-runTest(GOdata_BP_socialL, algorithm = "classic", statistic = "fisher")
BP_obj_Test_res_socialM <- GenTable( GOdata_BP_socialM,weight   = BPres_socialM,orderBy  = "weight",topNodes=67,numChar=1000)
BP_obj_Test_res_socialD <- GenTable( GOdata_BP_socialD,weight   = BPres_socialD,orderBy  = "weight",topNodes=54,numChar=1000)
BP_obj_Test_res_socialL <- GenTable( GOdata_BP_socialL,weight   = BPres_socialL,orderBy  = "weight",topNodes=144,numChar=1000)
MFres_socialM<-runTest(GOdata_MF_socialM, algorithm = "classic", statistic = "fisher")
MFres_socialD<-runTest(GOdata_MF_socialD, algorithm = "classic", statistic = "fisher")
MFres_socialL<-runTest(GOdata_MF_socialL, algorithm = "classic", statistic = "fisher")
MF_obj_Test_res_socialM <- GenTable( GOdata_MF_socialM, weight   = MFres_socialM,orderBy  = "weight",topNodes=26,numChar=1000)
MF_obj_Test_res_socialD <- GenTable( GOdata_MF_socialD, weight   = MFres_socialD,orderBy  = "weight",topNodes=14,numChar=1000)
MF_obj_Test_res_socialL <- GenTable( GOdata_MF_socialL, weight   = MFres_socialL,orderBy  = "weight",topNodes=41,numChar=1000)
CCres_socialM<-runTest(GOdata_CC_socialM, algorithm = "classic", statistic = "fisher")
CCres_socialD<-runTest(GOdata_CC_socialD, algorithm = "classic", statistic = "fisher")
CCres_socialL<-runTest(GOdata_CC_socialL, algorithm = "classic", statistic = "fisher")
CC_obj_Test_res_socialM <- GenTable( GOdata_CC_socialM,weight   = CCres_socialM,orderBy  = "weight",topNodes=2,numChar=1000)
CC_obj_Test_res_socialD <- GenTable( GOdata_CC_socialD,weight   = CCres_socialD,orderBy  = "weight",topNodes=2,numChar=1000)
CC_obj_Test_res_socialL <- GenTable( GOdata_CC_socialL,weight   = CCres_socialL,orderBy  = "weight",topNodes=36,numChar=1000)
BP_obj_Test_res_socialM$category<-"BP"
MF_obj_Test_res_socialM$category<-"MF"
CC_obj_Test_res_socialM$category<-"CC"
BP_obj_Test_res_socialD$category<-"BP"
MF_obj_Test_res_socialD$category<-"MF"
CC_obj_Test_res_socialD$category<-"CC"
BP_obj_Test_res_socialL$category<-"BP"
MF_obj_Test_res_socialL$category<-"MF"
CC_obj_Test_res_socialL$category<-"CC"
allGORes_socialM <- data.frame(rbind(BP_obj_Test_res_socialM, CC_obj_Test_res_socialM, MF_obj_Test_res_socialM))
allGORes_socialD <- data.frame(rbind(BP_obj_Test_res_socialD,  CC_obj_Test_res_socialD,MF_obj_Test_res_socialD))
allGORes_socialL <- data.frame(rbind(BP_obj_Test_res_socialL, CC_obj_Test_res_socialL, MF_obj_Test_res_socialL))
allGORes_socialM$weight <- as.numeric(allGORes_socialM$weight)
allGORes_socialD$weight <- as.numeric(allGORes_socialD$weight)
allGORes_socialL$weight <- as.numeric(allGORes_socialL$weight)
allGORes_socialM <- allGORes_socialM[order(allGORes_socialM$weight), ]
allGORes_socialD <- allGORes_socialD[order(allGORes_socialD$weight), ]
allGORes_socialL <- allGORes_socialL[order(allGORes_socialL$weight), ]
write.csv(allGORes_socialM, file = "allGORes_socialF.csv")
write.csv(allGORes_socialL, file = "allGORes_socialL.csv")
write.csv(allGORes_socialD, file = "allGORes_socialD.csv")


#The interaction of temp x parenting
options(scipen = 999)
tempxsocial.dad<-DEGs.tempxsocialM$gene_id #we will combine these lists because there is high concordance
tempxsocial.mum<-DEGs.tempxsocialF$gene_id
tempxsocial.larv<-DEGs.tempxsocialL$gene_id
GOList.tempxsocialD<-GOList[GOList$gene_id %in% tempxsocial.dad,]
GOList.tempxsocialM<-GOList[GOList$gene_id %in% tempxsocial.mum,]
GOList.tempxsocialL<-GOList[GOList$gene_id %in% tempxsocial.larv,]
write.table(GOList.tempxsocialD, file = "GOList.tempxsocialD.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
write.table(GOList.tempxsocialM, file = "GOList.tempxsocialM.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
write.table(GOList.tempxsocialL, file = "GOList.tempxsocialL.txt", sep = "\\t",row.names = FALSE, col.names = FALSE)
tempxsocialD.GOs<-readMappings(file="GOList.tempxsocialD.txt")
tempxsocialM.GOs<-readMappings(file="GOList.tempxsocialM.txt")
tempxsocialL.GOs<-readMappings(file="GOList.tempxsocialL.txt")
tempxsocialDNames<-names(tempxsocialD.GOs)
tempxsocialMNames<-names(tempxsocialM.GOs)
tempxsocialLNames<-names(tempxsocialL.GOs)
geneList_tempxsocialD <- factor(as.integer(geneNames %in% tempxsocialDNames))
geneList_tempxsocialM <- factor(as.integer(geneNames %in% tempxsocialMNames))
geneList_tempxsocialL <- factor(as.integer(geneNames %in% tempxsocialLNames))
names(geneList_tempxsocialD) <- geneNames
names(geneList_tempxsocialM) <- geneNames
names(geneList_tempxsocialL) <- geneNames
GOdata_MF_tempxsocialD <- new("topGOdata", ontology = "MF", allGenes = geneList_tempxsocialD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF_tempxsocialM <- new("topGOdata", ontology = "MF", allGenes = geneList_tempxsocialM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_MF_tempxsocialL <- new("topGOdata", ontology = "MF", allGenes = geneList_tempxsocialL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_tempxsocialD <- new("topGOdata", ontology = "BP", allGenes = geneList_tempxsocialD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_tempxsocialM <- new("topGOdata", ontology = "BP", allGenes = geneList_tempxsocialM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_BP_tempxsocialL <- new("topGOdata", ontology = "BP", allGenes = geneList_tempxsocialL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_tempxsocialD <- new("topGOdata", ontology = "CC", allGenes = geneList_tempxsocialD,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_tempxsocialM <- new("topGOdata", ontology = "CC", allGenes = geneList_tempxsocialM,annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata_CC_tempxsocialL <- new("topGOdata", ontology = "CC", allGenes = geneList_tempxsocialL,annot = annFUN.gene2GO, gene2GO = geneID2GO)
BPres_tempxsocialD<-runTest(GOdata_BP_tempxsocialD, algorithm = "classic", statistic = "fisher")
BPres_tempxsocialM<-runTest(GOdata_BP_tempxsocialM, algorithm = "classic", statistic = "fisher")
BPres_tempxsocialL<-runTest(GOdata_BP_tempxsocialL, algorithm = "classic", statistic = "fisher")
BP_obj_Test_res_tempxsocialD <- GenTable( GOdata_BP_tempxsocialD,weight   = BPres_tempxsocialD,orderBy  = "weight",topNodes=8,numChar=1000)
BP_obj_Test_res_tempxsocialM <- GenTable( GOdata_BP_tempxsocialM,weight   = BPres_tempxsocialM,orderBy  = "weight",topNodes=43,numChar=1000)
BP_obj_Test_res_tempxsocialL <- GenTable( GOdata_BP_tempxsocialL,weight   = BPres_tempxsocialL,orderBy  = "weight",topNodes=8,numChar=1000)
MFres_tempxsocialD<-runTest(GOdata_MF_tempxsocialD, algorithm = "classic", statistic = "fisher")
MFres_tempxsocialM<-runTest(GOdata_MF_tempxsocialM, algorithm = "classic", statistic = "fisher")
MFres_tempxsocialL<-runTest(GOdata_MF_tempxsocialL, algorithm = "classic", statistic = "fisher",p.adjust=TRUE)
MF_obj_Test_res_tempxsocialD <- GenTable( GOdata_MF_tempxsocialD, weight   = MFres_tempxsocialD,orderBy  = "weight",topNodes=4,numChar=1000)
MF_obj_Test_res_tempxsocialM <- GenTable( GOdata_MF_tempxsocialM, weight   = MFres_tempxsocialM,orderBy  = "weight",topNodes=21,numChar=1000)
MF_obj_Test_res_tempxsocialL <- GenTable( GOdata_MF_tempxsocialL, weight   = MFres_tempxsocialL,orderBy  = "weight",topNodes=5,numChar=1000)
CCres_tempxsocialD<-runTest(GOdata_CC_tempxsocialD, algorithm = "classic", statistic = "fisher")
CCres_tempxsocialM<-runTest(GOdata_CC_tempxsocialM, algorithm = "classic", statistic = "fisher")
CCres_tempxsocialL<-runTest(GOdata_CC_tempxsocialL, algorithm = "classic", statistic = "fisher")
CC_obj_Test_res_tempxsocialD <- GenTable( GOdata_CC_tempxsocialD,weight   = CCres_tempxsocialD,orderBy  = "weight",topNodes=17,numChar=1000)
CC_obj_Test_res_tempxsocialM <- GenTable( GOdata_CC_tempxsocialM,weight   = CCres_tempxsocialM,orderBy  = "weight",topNodes=6,numChar=1000)
CC_obj_Test_res_tempxsocialL <- GenTable( GOdata_CC_tempxsocialL,weight   = CCres_tempxsocialL,orderBy  = "weight",topNodes=2,numChar=1000)
BP_obj_Test_res_tempxsocialD$category<-"BP"
MF_obj_Test_res_tempxsocialD$category<-"MF"
CC_obj_Test_res_tempxsocialD$category<-"CC"
BP_obj_Test_res_tempxsocialM$category<-"BP"
MF_obj_Test_res_tempxsocialM$category<-"MF"
CC_obj_Test_res_tempxsocialM$category<-"CC"
BP_obj_Test_res_tempxsocialL$category<-"BP"
MF_obj_Test_res_tempxsocialL$category<-"MF"
CC_obj_Test_res_tempxsocialL$category<-"CC"
allGORes_tempxsocialD <- data.frame(rbind(BP_obj_Test_res_tempxsocialD, CC_obj_Test_res_tempxsocialD, MF_obj_Test_res_tempxsocialD))
allGORes_tempxsocialM <- data.frame(rbind(BP_obj_Test_res_tempxsocialM, CC_obj_Test_res_tempxsocialM, MF_obj_Test_res_tempxsocialM))
allGORes_tempxsocialL <- data.frame(rbind(BP_obj_Test_res_tempxsocialL, CC_obj_Test_res_tempxsocialL, MF_obj_Test_res_tempxsocialL))
allGORes_tempxsocialD$weight <- as.numeric(allGORes_tempxsocialD$weight)
allGORes_tempxsocialM$weight <- as.numeric(allGORes_tempxsocialM$weight)
allGORes_tempxsocialL$weight <- as.numeric(allGORes_tempxsocialL$weight)
allGORes_tempxsocialD <- allGORes_tempxsocialD[order(allGORes_tempxsocialD$weight), ]
allGORes_tempxsocialM <- allGORes_tempxsocialM[order(allGORes_tempxsocialM$weight), ]
allGORes_tempxsocialL <- allGORes_tempxsocialL[order(allGORes_tempxsocialL$weight), ]
write.csv(allGORes_tempxsocialD, file = "allGORes_tempxsocialM.csv")
write.csv(allGORes_tempxsocialM, file = "allGORes_tempxsocialF.csv")
write.csv(allGORes_tempxsocialL, file = "allGORes_tempxsocialL.csv")


#Examine the directionality of change of enriched pathways. 
#This requires an input with a list of genes associated with each GO term to calculate the mean direction of the enrichment.
getgenes <-function(x) {
  genes<-data.frame()
  for (i in 1:length(names(x))) {
    name<-names(x)[i]
    item<-cbind(name, toString(gsub("\\"","", x[name]),sep="," ))
    genes<-rbind.data.frame(genes, item)
  }
  return(genes)
}

#Start with female temperature DEGs
allGO.tempM_MF = genesInTerm(GOdata_MF_tempM)
allGO.tempM_BP = genesInTerm(GOdata_BP_tempM)
allGO.tempM_CC = genesInTerm(GOdata_CC_tempM)
df_allGO.tempM_MF<-getgenes(allGO.tempM_MF)
df_allGO.tempM_BP<-getgenes(allGO.tempM_BP)
df_allGO.tempM_CC<-getgenes(allGO.tempM_CC)
df_allGO.tempM<-rbind(df_allGO.tempM_MF,df_allGO.tempM_BP,df_allGO.tempM_CC)
colnames(df_allGO.tempM)<-c("ID","genes")
colnames(allGORes_tempM)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_tempM<-merge(allGORes_tempM,df_allGO.tempM,by="ID",all.x=T)
allGORes_tempM$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_tempM$genes))
allGORes_tempM$genes<-gsub("\\\\\\\\", "", allGORes_tempM$genes)
allGORes_tempM$genes<-gsub("\\n", "", allGORes_tempM$genes)

#read in gene info
DEGs.tempF$gene_id<-gsub(" ","",DEGs.tempF$gene_id)
colnames(DEGs.tempF)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_tempF<-direcle_dat(allGORes_tempM,DEGs.tempF)
direc_tempF<-direc_tempF[!duplicated(direc_tempF),]
direc_tempF<-direc_tempF[which(!is.na(direc_tempF$logFC)),]
#Manually calculate z score as (up-down)/squrt(count) following Walter, SÃ¡nchez-Cabo, & Ricote (2015)
counts_tempF<-allGORes_tempM[,c(1,4)] #correct count column. It didn't calculate it correctly.
colnames(counts_tempF)<-c("ID","count")
direc_tempF<-merge(direc_tempF[,-4],counts_tempF,all.x=T)
direc_tempF_summarized<-direc_tempF %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_tempF_summarized$zscore<-(direc_tempF_summarized$up-direc_tempF_summarized$down)/direc_tempF_summarized$count_sqrt
direc_tempF2<-merge(direc_tempF_summarized[,c(1,5)],direc_tempF[,c(1:3,6,8)],all.x=T)
direc_tempF2<-direc_tempF2[!duplicated(direc_tempF2),]

#Repeat for other DEG lists
#Temp-Males
allGO.tempD_MF = genesInTerm(GOdata_MF_tempD)
allGO.tempD_BP = genesInTerm(GOdata_BP_tempD)
allGO.tempD_CC = genesInTerm(GOdata_CC_tempD)
df_allGO.tempD_MF<-getgenes(allGO.tempD_MF)
df_allGO.tempD_BP<-getgenes(allGO.tempD_BP)
df_allGO.tempD_CC<-getgenes(allGO.tempD_CC)
df_allGO.tempD<-rbind(df_allGO.tempD_MF,df_allGO.tempD_BP,df_allGO.tempD_CC)
colnames(df_allGO.tempD)<-c("ID","genes")
colnames(allGORes_tempD)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_tempD<-merge(allGORes_tempD,df_allGO.tempD,by="ID",all.x=T)
allGORes_tempD$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_tempD$genes))
allGORes_tempD$genes<-gsub("\\\\\\\\", "", allGORes_tempD$genes)
allGORes_tempD$genes<-gsub("\\n", "", allGORes_tempD$genes)
DEGs.tempM$gene_id<-gsub(" ","",DEGs.tempM$gene_id)
colnames(DEGs.tempM)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_tempM<-direcle_dat(allGORes_tempD,DEGs.tempM)
direc_tempM<-direc_tempM[!duplicated(direc_tempM),]
direc_tempM<-direc_tempM[which(!is.na(direc_tempM$logFC)),]
counts_tempM<-allGORes_tempD[,c(1,4)] 
colnames(counts_tempM)<-c("ID","count")
direc_tempM<-merge(direc_tempM[,-4],counts_tempM,all.x=T)
direc_tempM_summarized<-direc_tempM %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_tempM_summarized$zscore<-(direc_tempM_summarized$up-direc_tempM_summarized$down)/direc_tempM_summarized$count_sqrt
direc_tempM2<-merge(direc_tempM_summarized[,c(1,5)],direc_tempM[,c(1:3,6,8)],all.x=T)
direc_tempM2<-direc_tempM2[!duplicated(direc_tempM2),]

#Temp-Larvae
allGO.tempL_MF = genesInTerm(GOdata_MF_tempL)
allGO.tempL_BP = genesInTerm(GOdata_BP_tempL)
allGO.tempL_CC = genesInTerm(GOdata_CC_tempL)
df_allGO.tempL_MF<-getgenes(allGO.tempL_MF)
df_allGO.tempL_BP<-getgenes(allGO.tempL_BP)
df_allGO.tempL_CC<-getgenes(allGO.tempL_CC)
df_allGO.tempL<-rbind(df_allGO.tempL_MF,df_allGO.tempL_BP,df_allGO.tempL_CC)
colnames(df_allGO.tempL)<-c("ID","genes")
colnames(allGORes_tempL)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_tempL<-merge(allGORes_tempL,df_allGO.tempL,by="ID",all.x=T)
allGORes_tempL$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_tempL$genes))
allGORes_tempL$genes<-gsub("\\\\\\\\", "", allGORes_tempL$genes)
allGORes_tempL$genes<-gsub("\\n", "", allGORes_tempL$genes)
DEGs.tempL$gene_id<-gsub(" ","",DEGs.tempL$gene_id)
colnames(DEGs.tempL)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_tempL<-direcle_dat(allGORes_tempL,DEGs.tempL)
direc_tempL<-direc_tempL[!duplicated(direc_tempL),]
direc_tempL<-direc_tempL[which(!is.na(direc_tempL$logFC)),]
counts_tempL<-allGORes_tempL[,c(1,4)] 
colnames(counts_tempL)<-c("ID","count")
direc_tempL<-merge(direc_tempL[,-4],counts_tempL,all.x=T)
direc_tempL_summarized<-direc_tempL %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_tempL_summarized$zscore<-(direc_tempL_summarized$up-direc_tempL_summarized$down)/direc_tempL_summarized$count_sqrt
direc_tempL2<-merge(direc_tempL_summarized[,c(1,5)],direc_tempL[,c(1:3,6,8)],all.x=T)
direc_tempL2<-direc_tempL2[!duplicated(direc_tempL2),]

#Parenting-Females
allGO.socialM_MF = genesInTerm(GOdata_MF_socialM)
allGO.socialM_BP = genesInTerm(GOdata_BP_socialM)
allGO.socialM_CC = genesInTerm(GOdata_CC_socialM)
df_allGO.socialM_MF<-getgenes(allGO.socialM_MF)
df_allGO.socialM_BP<-getgenes(allGO.socialM_BP)
df_allGO.socialM_CC<-getgenes(allGO.socialM_CC)
df_allGO.socialM<-rbind(df_allGO.socialM_MF,df_allGO.socialM_BP,df_allGO.socialM_CC)
colnames(df_allGO.socialM)<-c("ID","genes")
colnames(allGORes_socialM)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_socialM<-merge(allGORes_socialM,df_allGO.socialM,by="ID",all.x=T)
allGORes_socialM$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_socialM$genes))
allGORes_socialM$genes<-gsub("\\\\\\\\", "", allGORes_socialM$genes)
allGORes_socialM$genes<-gsub("\\n", "", allGORes_socialM$genes)

#read in gene info
DEGs.socialF$gene_id<-gsub(" ","",DEGs.socialF$gene_id)
colnames(DEGs.socialF)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_socialF<-direcle_dat(allGORes_socialM,DEGs.socialF)
direc_socialF<-direc_socialF[!duplicated(direc_socialF),]
direc_socialF<-direc_socialF[which(!is.na(direc_socialF$logFC)),]
#Manually calculate z score as (up-down)/squrt(count)
counts_socialF<-allGORes_socialM[,c(1,4)] #correct count column. It didn't calculate it correctly.
colnames(counts_socialF)<-c("ID","count")
direc_socialF<-merge(direc_socialF[,-4],counts_socialF,all.x=T)
direc_socialF_summarized<-direc_socialF %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_socialF_summarized$zscore<-(direc_socialF_summarized$up-direc_socialF_summarized$down)/direc_socialF_summarized$count_sqrt
direc_socialF2<-merge(direc_socialF_summarized[,c(1,5)],direc_socialF[,c(1:3,6,8)],all.x=T)
direc_socialF2<-direc_socialF2[!duplicated(direc_socialF2),]

#Parenting-Males
allGO.socialD_MF = genesInTerm(GOdata_MF_socialD)
allGO.socialD_BP = genesInTerm(GOdata_BP_socialD)
allGO.socialD_CC = genesInTerm(GOdata_CC_socialD)
df_allGO.socialD_MF<-getgenes(allGO.socialD_MF)
df_allGO.socialD_BP<-getgenes(allGO.socialD_BP)
df_allGO.socialD_CC<-getgenes(allGO.socialD_CC)
df_allGO.socialD<-rbind(df_allGO.socialD_MF,df_allGO.socialD_BP,df_allGO.socialD_CC)
colnames(df_allGO.socialD)<-c("ID","genes")
colnames(allGORes_socialD)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_socialD<-merge(allGORes_socialD,df_allGO.socialD,by="ID",all.x=T)
allGORes_socialD$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_socialD$genes))
allGORes_socialD$genes<-gsub("\\\\\\\\", "", allGORes_socialD$genes)
allGORes_socialD$genes<-gsub("\\n", "", allGORes_socialD$genes)
DEGs.socialM$gene_id<-gsub(" ","",DEGs.socialM$gene_id)
colnames(DEGs.socialM)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_socialM<-direcle_dat(allGORes_socialD,DEGs.socialM)
direc_socialM<-direc_socialM[!duplicated(direc_socialM),]
direc_socialM<-direc_socialM[which(!is.na(direc_socialM$logFC)),]
counts_socialM<-allGORes_socialD[,c(1,4)] 
colnames(counts_socialM)<-c("ID","count")
direc_socialM<-merge(direc_socialM[,-4],counts_socialM,all.x=T)
direc_socialM_summarized<-direc_socialM %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_socialM_summarized$zscore<-(direc_socialM_summarized$up-direc_socialM_summarized$down)/direc_socialM_summarized$count_sqrt
direc_socialM2<-merge(direc_socialM_summarized[,c(1,5)],direc_socialM[,c(1:3,6,8)],all.x=T)
direc_socialM2<-direc_socialM2[!duplicated(direc_socialM2),]

#Parented-Larvae
allGO.socialL_MF = genesInTerm(GOdata_MF_socialL)
allGO.socialL_BP = genesInTerm(GOdata_BP_socialL)
allGO.socialL_CC = genesInTerm(GOdata_CC_socialL)
df_allGO.socialL_MF<-getgenes(allGO.socialL_MF)
df_allGO.socialL_BP<-getgenes(allGO.socialL_BP)
df_allGO.socialL_CC<-getgenes(allGO.socialL_CC)
df_allGO.socialL<-rbind(df_allGO.socialL_MF,df_allGO.socialL_BP,df_allGO.socialL_CC)
colnames(df_allGO.socialL)<-c("ID","genes")
colnames(allGORes_socialL)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_socialL<-merge(allGORes_socialL,df_allGO.socialL,by="ID",all.x=T)
allGORes_socialL$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_socialL$genes))
allGORes_socialL$genes<-gsub("\\\\\\\\", "", allGORes_socialL$genes)
allGORes_socialL$genes<-gsub("\\n", "", allGORes_socialL$genes)
DEGs.socialL$gene_id<-gsub(" ","",DEGs.socialL$gene_id)
colnames(DEGs.socialL)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_socialL<-direcle_dat(allGORes_socialL,DEGs.socialL)
direc_socialL<-direc_socialL[!duplicated(direc_socialL),]
direc_socialL<-direc_socialL[which(!is.na(direc_socialL$logFC)),]
counts_socialL<-allGORes_socialL[,c(1,4)] 
colnames(counts_socialL)<-c("ID","count")
direc_socialL<-merge(direc_socialL[,-4],counts_socialL,all.x=T)
direc_socialL_summarized<-direc_socialL %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_socialL_summarized$zscore<-(direc_socialL_summarized$up-direc_socialL_summarized$down)/direc_socialL_summarized$count_sqrt
direc_socialL2<-merge(direc_socialL_summarized[,c(1,5)],direc_socialL[,c(1:3,6,8)],all.x=T)
direc_socialL2<-direc_socialL2[!duplicated(direc_socialL2),]

#Females-Temp x Parenting
allGO.tempxsocialM_MF = genesInTerm(GOdata_MF_tempxsocialM)
allGO.tempxsocialM_BP = genesInTerm(GOdata_BP_tempxsocialM)
allGO.tempxsocialM_CC = genesInTerm(GOdata_CC_tempxsocialM)
df_allGO.tempxsocialM_MF<-getgenes(allGO.tempxsocialM_MF)
df_allGO.tempxsocialM_BP<-getgenes(allGO.tempxsocialM_BP)
df_allGO.tempxsocialM_CC<-getgenes(allGO.tempxsocialM_CC)
df_allGO.tempxsocialM<-rbind(df_allGO.tempxsocialM_MF,df_allGO.tempxsocialM_BP,df_allGO.tempxsocialM_CC)
colnames(df_allGO.tempxsocialM)<-c("ID","genes")
colnames(allGORes_tempxsocialM)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_tempxsocialM<-merge(allGORes_tempxsocialM,df_allGO.tempxsocialM,by="ID",all.x=T)
allGORes_tempxsocialM$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_tempxsocialM$genes))
allGORes_tempxsocialM$genes<-gsub("\\\\\\\\", "", allGORes_tempxsocialM$genes)
allGORes_tempxsocialM$genes<-gsub("\\n", "", allGORes_tempxsocialM$genes)

#read in gene info
DEGs.tempxsocialF$gene_id<-gsub(" ","",DEGs.tempxsocialF$gene_id)
colnames(DEGs.tempxsocialF)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_tempxsocialF<-direcle_dat(allGORes_tempxsocialM,DEGs.tempxsocialF)
direc_tempxsocialF<-direc_tempxsocialF[!duplicated(direc_tempxsocialF),]
direc_tempxsocialF<-direc_tempxsocialF[which(!is.na(direc_tempxsocialF$logFC)),]
#Manually calculate z score as (up-down)/squrt(count)
counts_tempxsocialF<-allGORes_tempxsocialM[,c(1,4)] #correct count column. It didn't calculate it correctly.
colnames(counts_tempxsocialF)<-c("ID","count")
direc_tempxsocialF<-merge(direc_tempxsocialF[,-4],counts_tempxsocialF,all.x=T)
direc_tempxsocialF_summarized<-direc_tempxsocialF %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_tempxsocialF_summarized$zscore<-(direc_tempxsocialF_summarized$up-direc_tempxsocialF_summarized$down)/direc_tempxsocialF_summarized$count_sqrt
direc_tempxsocialF2<-merge(direc_tempxsocialF_summarized[,c(1,5)],direc_tempxsocialF[,c(1:3,6,8)],all.x=T)
direc_tempxsocialF2<-direc_tempxsocialF2[!duplicated(direc_tempxsocialF2),]

#tempxparenting-Males
allGO.tempxsocialD_MF = genesInTerm(GOdata_MF_tempxsocialD)
allGO.tempxsocialD_BP = genesInTerm(GOdata_BP_tempxsocialD)
allGO.tempxsocialD_CC = genesInTerm(GOdata_CC_tempxsocialD)
df_allGO.tempxsocialD_MF<-getgenes(allGO.tempxsocialD_MF)
df_allGO.tempxsocialD_BP<-getgenes(allGO.tempxsocialD_BP)
df_allGO.tempxsocialD_CC<-getgenes(allGO.tempxsocialD_CC)
df_allGO.tempxsocialD<-rbind(df_allGO.tempxsocialD_MF,df_allGO.tempxsocialD_BP,df_allGO.tempxsocialD_CC)
colnames(df_allGO.tempxsocialD)<-c("ID","genes")
colnames(allGORes_tempxsocialD)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_tempxsocialD<-merge(allGORes_tempxsocialD,df_allGO.tempxsocialD,by="ID",all.x=T)
allGORes_tempxsocialD$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_tempxsocialD$genes))
allGORes_tempxsocialD$genes<-gsub("\\\\\\\\", "", allGORes_tempxsocialD$genes)
allGORes_tempxsocialD$genes<-gsub("\\n", "", allGORes_tempxsocialD$genes)
DEGs.tempxsocialM$gene_id<-gsub(" ","",DEGs.tempxsocialM$gene_id)
colnames(DEGs.tempxsocialM)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_tempxsocialM<-direcle_dat(allGORes_tempxsocialD,DEGs.tempxsocialM)
direc_tempxsocialM<-direc_tempxsocialM[!duplicated(direc_tempxsocialM),]
direc_tempxsocialM<-direc_tempxsocialM[which(!is.na(direc_tempxsocialM$logFC)),]
counts_tempxsocialM<-allGORes_tempxsocialD[,c(1,4)] 
colnames(counts_tempxsocialM)<-c("ID","count")
direc_tempxsocialM<-merge(direc_tempxsocialM[,-4],counts_tempxsocialM,all.x=T)
direc_tempxsocialM_summarized<-direc_tempxsocialM %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_tempxsocialM_summarized$zscore<-(direc_tempxsocialM_summarized$up-direc_tempxsocialM_summarized$down)/direc_tempxsocialM_summarized$count_sqrt
direc_tempxsocialM2<-merge(direc_tempxsocialM_summarized[,c(1,5)],direc_tempxsocialM[,c(1:3,6,8)],all.x=T)
direc_tempxsocialM2<-direc_tempxsocialM2[!duplicated(direc_tempxsocialM2),]

#tempxparented-Larvae
allGO.tempxsocialL_MF = genesInTerm(GOdata_MF_tempxsocialL)
allGO.tempxsocialL_BP = genesInTerm(GOdata_BP_tempxsocialL)
allGO.tempxsocialL_CC = genesInTerm(GOdata_CC_tempxsocialL)
df_allGO.tempxsocialL_MF<-getgenes(allGO.tempxsocialL_MF)
df_allGO.tempxsocialL_BP<-getgenes(allGO.tempxsocialL_BP)
df_allGO.tempxsocialL_CC<-getgenes(allGO.tempxsocialL_CC)
df_allGO.tempxsocialL<-rbind(df_allGO.tempxsocialL_MF,df_allGO.tempxsocialL_BP,df_allGO.tempxsocialL_CC)
colnames(df_allGO.tempxsocialL)<-c("ID","genes")
colnames(allGORes_tempxsocialL)<-c("ID","Term","Annotated","Significant","Expected","adj_pval","category")
allGORes_tempxsocialL<-merge(allGORes_tempxsocialL,df_allGO.tempxsocialL,by="ID",all.x=T)
allGORes_tempxsocialL$genes<-as.character(gsub("c\\\\(|\\\\)",'',allGORes_tempxsocialL$genes))
allGORes_tempxsocialL$genes<-gsub("\\\\\\\\", "", allGORes_tempxsocialL$genes)
allGORes_tempxsocialL$genes<-gsub("\\n", "", allGORes_tempxsocialL$genes)
DEGs.tempxsocialL$gene_id<-gsub(" ","",DEGs.tempxsocialL$gene_id)
colnames(DEGs.tempxsocialL)<-c("ID","logFC", "p","Parent_gene","Protein_ID","Accession", "GeneID","Locus","Protein_Name")
direc_tempxsocialL<-direcle_dat(allGORes_tempxsocialL,DEGs.tempxsocialL)
direc_tempxsocialL<-direc_tempxsocialL[!duplicated(direc_tempxsocialL),]
direc_tempxsocialL<-direc_tempxsocialL[which(!is.na(direc_tempxsocialL$logFC)),]
counts_tempxsocialL<-allGORes_tempxsocialL[,c(1,4)] 
colnames(counts_tempxsocialL)<-c("ID","count")
direc_tempxsocialL<-merge(direc_tempxsocialL[,-4],counts_tempxsocialL,all.x=T)
direc_tempxsocialL_summarized<-direc_tempxsocialL %>%
  group_by(ID) %>%
  summarise(up = sum(logFC>0),down = sum(logFC<0),count_sqrt=sqrt(mean(count)))
direc_tempxsocialL_summarized$zscore<-(direc_tempxsocialL_summarized$up-direc_tempxsocialL_summarized$down)/direc_tempxsocialL_summarized$count_sqrt
direc_tempxsocialL2<-merge(direc_tempxsocialL_summarized[,c(1,5)],direc_tempxsocialL[,c(1:3,6,8)],all.x=T)
direc_tempxsocialL2<-direc_tempxsocialL2[!duplicated(direc_tempxsocialL2),]

