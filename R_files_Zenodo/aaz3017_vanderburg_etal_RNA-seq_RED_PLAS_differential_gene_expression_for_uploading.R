
#Van der Burg et al.
#Genomic architecture of a genetically assimilated seasonal color pattern
#Supplemental script
#RNAseq analysis

library("DESeq2")
library(plyr)
library(ggplot2)
setwd("~/Box Sync/A.Cornell/Research/RNAseq")



counts <- as.data.frame(read.csv("RNAseq_RED_PLAS_HW_readcount.txt",sep=" "), header=TRUE)

ncol(counts)
names(counts)
sampleinfo<-NULL
sampleinfo$names<-names(counts)
sampleinfo$stage<-c('5th','5th','5th','pp','pp','pp','72h','72h','72h','d6','d6','d6','5th','5th','5th','pp','pp','pp','72h','72h','72h','d6','d6','d6')
sampleinfo$pop<-c("PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","RED","RED","RED","RED","RED","RED","RED","RED","RED","RED","RED","RED")
sampleinfo<-as.data.frame(sampleinfo)



dds_RNA_RED_PLAS <- DESeqDataSetFromMatrix(countData = counts,
                                colData = sampleinfo,
                                design = ~ stage + pop )
dds_RNA_RED_PLAS <- dds_RNA_RED_PLAS[ rowSums(counts(dds_RNA_RED_PLAS)) > 1, ]

dds_RP<-dds_RNA_RED_PLAS
dds_RP$group <- factor(paste0(dds_RP$stage, dds_RP$pop))
design(dds_RP) <- ~ group
dds_RP <- DESeq(dds_RP)




dds_RPres<-results(dds_RP)


dds_RPresord<-subset(dds_RPres[order(dds_RPres$padj),], padj < 0.05)

dds_RPres_5th<-results(dds_RP, contrast=c("group", "5thPLAS", "5thRED"))
dds_RPres_5th_ord<-subset(dds_RPres_5th[order(dds_RPres_5th$padj),], padj < 0.05)
length(dds_RPres_5th_ord$baseMean) ##300


dds_RPres_pp<-results(dds_RP, contrast=c("group", "ppPLAS", "ppRED"))
dds_RPres_pp_ord<-subset(dds_RPres_pp[order(dds_RPres_pp$padj),], padj < 0.05)
length(dds_RPres_pp_ord$baseMean) ##620

dds_RPres_72h<-results(dds_RP, contrast=c("group", "72hPLAS", "72hRED"))
dds_RPres_72h_ord<-subset(dds_RPres_72h[order(dds_RPres_72h$padj),], padj < 0.05)
length(dds_RPres_72h_ord$baseMean) ##664

dds_RPres_d6<-results(dds_RP, contrast=c("group", "d6PLAS", "d6RED"))
#save(dds_RPres_d6,file="dds_RNA_RED_PLAS_ommo-diff.rData")
dds_RPres_d6_ord<-subset(dds_RPres_d6[order(dds_RPres_d6$padj),], padj < 0.05)
length(dds_RPres_d6_ord$baseMean) ##2166



#plot all genes with p-values
genes_padj<-NULL
genes_padj$genes<-c(row.names(dds_RPres_5th),row.names(dds_RPres_pp),row.names(dds_RPres_72h),row.names(dds_RPres_d6))
genes_padj$padj<-c(dds_RPres_5th$padj,dds_RPres_pp$padj,dds_RPres_72h$padj,dds_RPres_d6$padj)
genes_padj<-as.data.frame(genes_padj)

plot(-log10(sort(genes_padj$padj,decreasing = TRUE)),
     col=ifelse(-log10(sort(genes_padj$padj,decreasing = TRUE))>-log10(0.05), "#ED7A34","black"))


genes_padj_05<-subset(genes_padj,genes_padj$padj<0.05)


rld_RP <- rlog(dds_RP, blind = FALSE)


##Plot different transformations
library("dplyr")
library("hexbin")
library("pheatmap")
library("RColorBrewer")


##calculate sample distances
sampleDists <- dist(t(assay(rld_RP)))
sampleDists




sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( sampleinfo$stage, sampleinfo$pop, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)

jpeg('RED_PLAS_heatmap.jpg', res = 300 , width = (3000*.8), height = (2250*.8))
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()


##DE gene expression plot

##Use summarySE function:

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}






##Get differential gene expression of genes close to SNPS (SNPS_closest_genes.txt file created in SNP analysis)

GOI<-read.table("SNPS_closest_genes.txt",header=F)
GOI<-GOI$V1[-which(GOI$V1=="JC_0011877-RA")]
GOI<-GOI[-which(GOI=="JC_0013043-RA")]
GOI<-GOI[-which(GOI=="JC_0013058-RA")]
GOI<-c(as.character(GOI),"JC_0019234-RA")


SNPS_GOI<-data.frame(matrix(0, ncol = 25, nrow = length(GOI)))
SNPS_GOI_cols<-c('gene','RNA_PLAS_5th_HW_S1', 'RNA_PLAS_5th_HW_S2', 'RNA_PLAS_5th_HW_S3', 
                     'RNA_PLAS_pp_HW_S1',  'RNA_PLAS_pp_HW_S2',  'RNA_PLAS_pp_HW_S3',
                     'RNA_PLAS_72h_HW_S1',  'RNA_PLAS_72h_HW_S2',  'RNA_PLAS_72h_HW_S3',  
                     'RNA_PLAS_d6_HW_S1',  'RNA_PLAS_d6_HW_S2',  'RNA_PLAS_d6_HW_S3',
                     'RNA_RED_5th_HW_S1', 'RNA_RED_5th_HW_S2', 'RNA_RED_5th_HW_S3', 
                     'RNA_RED_pp_HW_S1',  'RNA_RED_pp_HW_S2',  'RNA_RED_pp_HW_S3',
                     'RNA_RED_72h_HW_S1',  'RNA_RED_72h_HW_S2',  'RNA_RED_72h_HW_S3',  
                     'RNA_RED_d6_HW_S1',  'RNA_RED_d6_HW_S2',  'RNA_RED_d6_HW_S3')

SNPS_GOI<-data.frame(matrix(0, ncol = length(SNPS_GOI_cols), nrow = length(GOI)))
colnames(SNPS_GOI)<-SNPS_GOI_cols



for(i in 1:length(GOI)){
  SNPS_GOI$gene[i]<-as.character(GOI[i])
  y<-as.character(GOI[i])
  x<-plotCounts(dds_RP, gene=as.character(GOI[i]), intgroup=c('group'), returnData=TRUE)
  SNPS_GOI[i,2:25]<-x$count
}

dds_RP_GOI<-dds_RP[as.character(GOI),]

dds_RP_GOI_5th<-results(dds_RP_GOI, contrast=c("group", "5thPLAS", "5thRED"))
dds_RP_GOI_5th_subs<-dds_RP_GOI_5th[which(dds_RP_GOI_5th$padj<0.05),]

length(dds_RP_GOI_5th_subs$baseMean) ##1: JC_0019234-RA
dds_RP_GOI_5th_subs$log2FoldChange # -4.004187

dds_RP_GOI_pp<-results(dds_RP_GOI, contrast=c("group", "ppPLAS", "ppRED"))
dds_RP_GOI_pp_subs<-dds_RP_GOI_pp[which(dds_RP_GOI_pp$padj<0.05),]
length(dds_RP_GOI_pp_subs$baseMean) ##2: JC_0013047-RA, JC_0013081-RA
dds_RP_GOI_pp_subs$log2FoldChange  #-0.9286182 -1.3790571

dds_RP_GOI_72h<-results(dds_RP_GOI, contrast=c("group", "72hPLAS", "72hRED"))
dds_RP_GOI_72h_subs<-dds_RP_GOI_72h[which(dds_RP_GOI_72h$padj<0.05),]
length(dds_RP_GOI_72h_subs$baseMean) ##0

dds_RP_GOI_d6<-results(dds_RP_GOI, contrast=c("group", "d6PLAS", "d6RED"))
dds_RP_GOI_d6_subs<-dds_RP_GOI_d6[which(dds_RP_GOI_d6$padj<0.05),]
dds_RP_GOI_d6_subs2<-dds_RP_GOI_d6_subs[which(dds_RP_GOI_d6_subs$log2FoldChange>2|dds_RP_GOI_d6_subs$log2FoldChange<(-2)),]





