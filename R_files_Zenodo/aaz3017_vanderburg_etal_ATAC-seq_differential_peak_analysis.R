##RED_PLAS ATACseq analysis

library("DESeq2")

#Files created by "aaz3017_vanderburg_etal_ATACseq_RED_PLAS_peak_calling.txt"
P_count <- as.data.frame(read.csv("PLAS_RED_PLAS_HW_peaks.readcount",sep="\\t"), header=F)
R_count <- as.data.frame(read.csv("RED_RED_PLAS_HW.sorted.merged.readcount",sep="\\t"), header=FALSE)

colnames(R_count)<-c('contig','start','end','RED_5th_HW_S1_2','RED_5th_HW_S2_3','RED_5th_HW_S3_3','RED_pp_HW_S1',
'RED_pp_HW_S2','RED_pp_HW_S3_2','RED_72h_HW_S1_2','RED_72h_HW_S2','RED_72h_HW_S3_2',
'RED_d6_HW_S1_2','RED_d6_HW_S2_3','RED_d6_HW_S3_2')
colnames(P_count)<-c('contig','start','end','PLAS_5th_HW_S1_2','PLAS_5th_HW_S2_2','PLAS_5th_HW_S3_2',
'PLAS_72h_HW_S1_2','PLAS_72h_HW_S2_2','PLAS_72h_HW_S3_2','PLAS_d6_HW_S1','PLAS_d6_HW_S2','PLAS_d6_HW_S3',
'PLAS_pp_HW_S1','PLAS_pp_HW_S2','PLAS_pp_HW_S3')

colnames(P_count)
P_count2<-P_count[,c(1,2,3,4,5,6,13,14,15,7,8,9,10,11,12)]
P_count<-P_count2

P_count$peak <- paste(P_count$contig,":",P_count$start,"-",P_count$end,sep="")
R_count$peak <- paste(P_count$contig,":",P_count$start,"-",P_count$end,sep="")

P_count$contig<-NULL
P_count$end<-NULL
P_count$start<-NULL

R_count$contig<-NULL
R_count$end<-NULL
R_count$start<-NULL

RP_count <- merge(R_count, P_count, by = 'peak')
rownames(RP_count) <- RP_count$peak
RP_count$peak<-NULL
RP_count[1:5,]

RP_col<-matrix(c("RED","RED","RED","RED","RED","RED","RED","RED","RED","RED","RED","RED",
                     "PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS",
                     "PLAS","PLAS", '5th','5th','5th','pp','pp','pp','72h','72h','72h','d6',
                     'd6','d6','5th','5th','5th','pp','pp','pp','72h','72h','72h','d6','d6',
                     'd6'), nrow=24,ncol=2)
          
                   

rownames(RP_col)<- c("RED_5th_HW_S1_2","RED_5th_HW_S2_3","RED_5th_HW_S3_3", "RED_pp_HW_S1", "RED_pp_HW_S2",  
                     "RED_pp_HW_S3_2","RED_72h_HW_S1_2", "RED_72h_HW_S2", "RED_72h_HW_S3_2" , "RED_d6_HW_S1_2",
                     "RED_d6_HW_S2_3", "RED_d6_HW_S3_2","PLAS_5th_HW_S1_2", "PLAS_5th_HW_S2_2","PLAS_5th_HW_S3_2",
                     "PLAS_pp_HW_S1", "PLAS_pp_HW_S2", "PLAS_pp_HW_S3", "PLAS_72h_HW_S1_2", "PLAS_72h_HW_S2_2",
                     "PLAS_72h_HW_S3_2", "PLAS_d6_HW_S1", "PLAS_d6_HW_S2", "PLAS_d6_HW_S3")
colnames(RP_col)<-c("population","stage")


ddsRP <- DESeqDataSetFromMatrix(countData = RP_count,
                                colData = RP_col,
                                design = ~ stage + population)
ddsRP <- ddsRP[ rowSums(counts(ddsRP)) > 1, ]



dds_ATAC_RPres_5th<-results(ddsRPtest, contrast=c("group", "5thPLAS", "5thRED"))
dds_ATAC_RPres_5th_ord<-subset(dds_ATAC_RPres_5th[order(dds_ATAC_RPres_5th$padj),], padj < 0.05)
length(dds_ATAC_RPres_5th_ord$baseMean) ##520

row.names(dds_ATAC_RPres_5th[grep("000191F",row.names(dds_ATAC_RPres_5th)),])[700:800]
dds_ATAC_RPres_5th["000191F:515919-516171",]

dds_ATAC_RPres_pp<-results(ddsRPtest, contrast=c("group", "ppPLAS", "ppRED"))
dds_ATAC_RPres_pp_ord<-subset(dds_ATAC_RPres_pp[order(dds_ATAC_RPres_pp$padj),], padj < 0.05)
length(dds_ATAC_RPres_pp_ord$baseMean) ##2119

dds_ATAC_RPres_72h<-results(ddsRPtest, contrast=c("group", "72hPLAS", "72hRED"))
dds_ATAC_RPres_72h_ord<-subset(dds_ATAC_RPres_72h[order(dds_ATAC_RPres_72h$padj),], padj < 0.05)
length(dds_ATAC_RPres_72h_ord$baseMean) ##3712

dds_ATAC_RPres_d6<-results(ddsRPtest, contrast=c("group", "d6PLAS", "d6RED"))
dds_ATAC_RPres_d6_ord<-subset(dds_ATAC_RPres_d6[order(dds_ATAC_RPres_d6$padj),], padj < 0.05)
length(dds_ATAC_RPres_d6_ord$baseMean) ##45544


##write significant peaks to bedfile:
##5th instar
dds_ATAC_RPres_5th_ord_df<-dds_ATAC_RPres_5th_ord
dds_ATAC_RPres_5th_ord_df$peaks<-row.names(dds_ATAC_RPres_5th_ord)
dds_ATAC_RPres_5th_ord_df.vec<-as.vector(c(dds_ATAC_RPres_5th_ord_df$peaks))

dds_ATAC_RPres_5th_ord_df.vec.ord<-sort(dds_ATAC_RPres_5th_ord_df.vec)

dds_ATAC_RPres_5th.bed<-read.table(text=dds_ATAC_RPres_5th_ord_df.vec.ord,sep=":",colClasses = "character")
dds_ATAC_RPres_5th.bed2<-read.table(text=dds_ATAC_RPres_5th.bed$V2,sep="-",colClasses = "character")
dds_ATAC_RPres_5th.bed$V2<-dds_ATAC_RPres_5th.bed2$V1
dds_ATAC_RPres_5th.bed$V3<-dds_ATAC_RPres_5th.bed2$V2
write.table(dds_ATAC_RPres_5th.bed,file="dds_ATAC_RPres_5th.bed", quote=FALSE,col.names = FALSE,row.names=FALSE,sep = "\\t")


##pp
dds_ATAC_RPres_pp_ord_df<-dds_ATAC_RPres_pp_ord
dds_ATAC_RPres_pp_ord_df$peaks<-row.names(dds_ATAC_RPres_pp_ord)
dds_ATAC_RPres_pp_ord_df.vec<-as.vector(c(dds_ATAC_RPres_pp_ord_df$peaks))

dds_ATAC_RPres_pp_ord_df.vec.ord<-sort(dds_ATAC_RPres_pp_ord_df.vec)

dds_ATAC_RPres_pp.bed<-read.table(text=dds_ATAC_RPres_pp_ord_df.vec.ord,sep=":",colClasses = "character")
dds_ATAC_RPres_pp.bed2<-read.table(text=dds_ATAC_RPres_pp.bed$V2,sep="-",colClasses = "character")
dds_ATAC_RPres_pp.bed$V2<-dds_ATAC_RPres_pp.bed2$V1
dds_ATAC_RPres_pp.bed$V3<-dds_ATAC_RPres_pp.bed2$V2
write.table(dds_ATAC_RPres_pp.bed,file="dds_ATAC_RPres_pp.bed", quote=FALSE,col.names = FALSE,row.names=FALSE,sep = "\\t")


##72h
dds_ATAC_RPres_72h_ord_df<-dds_ATAC_RPres_72h_ord
dds_ATAC_RPres_72h_ord_df$peaks<-row.names(dds_ATAC_RPres_72h_ord)
dds_ATAC_RPres_72h_ord_df.vec<-as.vector(c(dds_ATAC_RPres_72h_ord_df$peaks))

dds_ATAC_RPres_72h_ord_df.vec.ord<-sort(dds_ATAC_RPres_72h_ord_df.vec)

dds_ATAC_RPres_72h.bed<-read.table(text=dds_ATAC_RPres_72h_ord_df.vec.ord,sep=":",colClasses = "character")
dds_ATAC_RPres_72h.bed2<-read.table(text=dds_ATAC_RPres_72h.bed$V2,sep="-",colClasses = "character")
dds_ATAC_RPres_72h.bed$V2<-dds_ATAC_RPres_72h.bed2$V1
dds_ATAC_RPres_72h.bed$V3<-dds_ATAC_RPres_72h.bed2$V2
write.table(dds_ATAC_RPres_72h.bed,file="dds_ATAC_RPres_72h.bed", quote=FALSE,col.names = FALSE,row.names=FALSE,sep = "\\t")

##d6
dds_ATAC_RPres_d6_ord_df<-dds_ATAC_RPres_d6_ord
dds_ATAC_RPres_d6_ord_df$peaks<-row.names(dds_ATAC_RPres_d6_ord)
dds_ATAC_RPres_d6_ord_df.vec<-as.vector(c(dds_ATAC_RPres_d6_ord_df$peaks))

dds_ATAC_RPres_d6_ord_df.vec.ord<-sort(dds_ATAC_RPres_d6_ord_df.vec)

dds_ATAC_RPres_d6.bed<-read.table(text=dds_ATAC_RPres_d6_ord_df.vec.ord,sep=":",colClasses = "character")
dds_ATAC_RPres_d6.bed2<-read.table(text=dds_ATAC_RPres_d6.bed$V2,sep="-",colClasses = "character")
dds_ATAC_RPres_d6.bed$V2<-dds_ATAC_RPres_d6.bed2$V1
dds_ATAC_RPres_d6.bed$V3<-dds_ATAC_RPres_d6.bed2$V2
write.table(dds_ATAC_RPres_d6.bed,file="dds_ATAC_RPres_d6.bed", quote=FALSE,col.names = FALSE,row.names=FALSE,sep = "\\t")



### get peaks in regions of interest
##Get 5th instar peaks
dds_ATAC_RPres_5th_ord_df<-dds_ATAC_RPres_5th_ord
dds_ATAC_RPres_5th_ord_df$peaks<-row.names(dds_ATAC_RPres_5th_ord)
dds_ATAC_RPres_5th_ord_df.vec<-as.vector(c(dds_ATAC_RPres_5th_ord_df$peaks))

dds_ATAC_RPres_5th_ord_df.vec.ord<-sort(dds_ATAC_RPres_5th_ord_df.vec)


dds_ATAC_RPres_5th_ord_df.vec.GOI<-
  c(dds_ATAC_RPres_5th_ord_df.vec.ord[grep("000044F",dds_ATAC_RPres_5th_ord_df.vec.ord)],dds_ATAC_RPres_5th_ord_df.vec.ord[grep("000090F",dds_ATAC_RPres_5th_ord_df.vec.ord)],
    dds_ATAC_RPres_5th_ord_df.vec.ord[grep("000152F",dds_ATAC_RPres_5th_ord_df.vec.ord)],dds_ATAC_RPres_5th_ord_df.vec.ord[grep("000191F",dds_ATAC_RPres_5th_ord_df.vec.ord)])

dds_ATAC_RPres_5th.GOI.bed<-read.table(text=dds_ATAC_RPres_5th_ord_df.vec.GOI,sep=":",colClasses = "character")
dds_ATAC_RPres_5th.GOI.bed2<-read.table(text=dds_ATAC_RPres_5th.GOI.bed$V2,sep="-",colClasses = "character")
dds_ATAC_RPres_5th.GOI.bed$V2<-dds_ATAC_RPres_5th.GOI.bed2$V1
dds_ATAC_RPres_5th.GOI.bed$V3<-dds_ATAC_RPres_5th.GOI.bed2$V2
write.table(dds_ATAC_RPres_5th.GOI.bed,file="dds_ATAC_RPres_5th.GOI.bed", quote=FALSE,col.names = FALSE,row.names=FALSE,sep = "\\t")

##pp peaks
dds_ATAC_RPres_pp_ord_df<-dds_ATAC_RPres_pp_ord
dds_ATAC_RPres_pp_ord_df$peaks<-row.names(dds_ATAC_RPres_pp_ord)
dds_ATAC_RPres_pp_ord_df.vec<-as.vector(c(dds_ATAC_RPres_pp_ord_df$peaks))

dds_ATAC_RPres_pp_ord_df.vec.ord<-sort(dds_ATAC_RPres_pp_ord_df.vec)


dds_ATAC_RPres_pp_ord_df.vec.GOI<-
  c(dds_ATAC_RPres_pp_ord_df.vec.ord[grep("000044F",dds_ATAC_RPres_pp_ord_df.vec.ord)],dds_ATAC_RPres_pp_ord_df.vec.ord[grep("000090F",dds_ATAC_RPres_pp_ord_df.vec.ord)],
    dds_ATAC_RPres_pp_ord_df.vec.ord[grep("000152F",dds_ATAC_RPres_pp_ord_df.vec.ord)],dds_ATAC_RPres_pp_ord_df.vec.ord[grep("000191F",dds_ATAC_RPres_pp_ord_df.vec.ord)])

dds_ATAC_RPres_pp.GOI.bed<-read.table(text=dds_ATAC_RPres_pp_ord_df.vec.GOI,sep=":",colClasses = "character")
dds_ATAC_RPres_pp.GOI.bed2<-read.table(text=dds_ATAC_RPres_pp.GOI.bed$V2,sep="-",colClasses = "character")
dds_ATAC_RPres_pp.GOI.bed$V2<-dds_ATAC_RPres_pp.GOI.bed2$V1
dds_ATAC_RPres_pp.GOI.bed$V3<-dds_ATAC_RPres_pp.GOI.bed2$V2
write.table(dds_ATAC_RPres_pp.GOI.bed,file="dds_ATAC_RPres_pp.GOI.bed", quote=FALSE,col.names = FALSE,row.names=FALSE,sep = "\\t")

##72h peaks

dds_ATAC_RPres_72h_ord_df<-dds_ATAC_RPres_72h_ord
dds_ATAC_RPres_72h_ord_df$peaks<-row.names(dds_ATAC_RPres_72h_ord)
dds_ATAC_RPres_72h_ord_df.vec<-as.vector(c(dds_ATAC_RPres_72h_ord_df$peaks))

dds_ATAC_RPres_72h_ord_df.vec.ord<-sort(dds_ATAC_RPres_72h_ord_df.vec)


dds_ATAC_RPres_72h_ord_df.vec.GOI<-
  c(dds_ATAC_RPres_72h_ord_df.vec.ord[grep("000044F",dds_ATAC_RPres_72h_ord_df.vec.ord)],dds_ATAC_RPres_72h_ord_df.vec.ord[grep("000090F",dds_ATAC_RPres_72h_ord_df.vec.ord)],
    dds_ATAC_RPres_72h_ord_df.vec.ord[grep("000152F",dds_ATAC_RPres_72h_ord_df.vec.ord)],dds_ATAC_RPres_72h_ord_df.vec.ord[grep("000191F",dds_ATAC_RPres_72h_ord_df.vec.ord)])

dds_ATAC_RPres_72h.GOI.bed<-read.table(text=dds_ATAC_RPres_72h_ord_df.vec.GOI,sep=":",colClasses = "character")
dds_ATAC_RPres_72h.GOI.bed2<-read.table(text=dds_ATAC_RPres_72h.GOI.bed$V2,sep="-",colClasses = "character")
dds_ATAC_RPres_72h.GOI.bed$V2<-dds_ATAC_RPres_72h.GOI.bed2$V1
dds_ATAC_RPres_72h.GOI.bed$V3<-dds_ATAC_RPres_72h.GOI.bed2$V2
write.table(dds_ATAC_RPres_72h.GOI.bed,file="dds_ATAC_RP_72h.GOI.bed", quote=FALSE,col.names = FALSE,row.names=FALSE,sep = "\\t")



##Get d6 GOI peaks
dds_ATAC_RPres_d6_ord_df<-dds_ATAC_RPres_d6_ord
dds_ATAC_RPres_d6_ord_df$peaks<-row.names(dds_ATAC_RPres_d6_ord)
dds_ATAC_RPres_d6_ord_df.vec<-as.vector(c(dds_ATAC_RPres_d6_ord_df$peaks))

dds_ATAC_RPres_d6_ord_df.vec.ord<-sort(dds_ATAC_RPres_d6_ord_df.vec)


dds_ATAC_RPres_d6_ord_df.vec.GOI<-
  c(dds_ATAC_RPres_d6_ord_df.vec.ord[grep("000044F",dds_ATAC_RPres_d6_ord_df.vec.ord)],dds_ATAC_RPres_d6_ord_df.vec.ord[grep("000090F",dds_ATAC_RPres_d6_ord_df.vec.ord)],
    dds_ATAC_RPres_d6_ord_df.vec.ord[grep("000152F",dds_ATAC_RPres_d6_ord_df.vec.ord)],dds_ATAC_RPres_d6_ord_df.vec.ord[grep("000191F",dds_ATAC_RPres_d6_ord_df.vec.ord)])

dds_ATAC_RPres_d6.GOI.bed<-read.table(text=dds_ATAC_RPres_d6_ord_df.vec.GOI,sep=":",colClasses = "character")
dds_ATAC_RPres_d6.GOI.bed2<-read.table(text=dds_ATAC_RPres_d6.GOI.bed$V2,sep="-",colClasses = "character")
dds_ATAC_RPres_d6.GOI.bed$V2<-dds_ATAC_RPres_d6.GOI.bed2$V1
dds_ATAC_RPres_d6.GOI.bed$V3<-dds_ATAC_RPres_d6.GOI.bed2$V2
write.table(dds_ATAC_RPres_d6.GOI.bed,file="dds_ATAC_RP_d6.GOI.bed", quote=FALSE,col.names = FALSE,row.names=FALSE,sep = "\\t")

##RED PLAS heatmap

names(RP_count)
sampleinfo<-NULL
sampleinfo$names<-names(RP_count)
sampleinfo$stage<-c('5th','5th','5th','pp','pp','pp','72h','72h','72h','d6','d6','d6','5th','5th','5th','pp','pp','pp','72h','72h','72h','d6','d6','d6')
sampleinfo$pop<-c("RED","RED","RED","RED","RED","RED","RED","RED","RED","RED","RED","RED","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS","PLAS")
sampleinfo<-as.data.frame(sampleinfo)


ATAC_vst_RP <- vst(ddsRPtest, blind = FALSE)


sampleDistsATAC <- dist(t(assay(ATAC_vst_RP)))
sampleDistsATAC

library("pheatmap")
library("RColorBrewer")

sampleDistMatrixATAC <- as.matrix( sampleDistsRNA )
rownames(sampleDistMatrixATAC) <- paste( sampleinfo$stage, sampleinfo$pop, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrixATAC,
         clustering_distance_rows = sampleDistsATAC,
         clustering_distance_cols = sampleDistsATAC,
         col = colors)

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddsRPtest)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( sampleinfo$stage, sampleinfo$pop, sep = " - " )
colnames(samplePoisDistMatrix) <- NULL

jpeg('RED_PLAS_ATAC_heatmap.jpg', res = 300 , width = (3000*.8), height = (2250*.8))
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
dev.off()


