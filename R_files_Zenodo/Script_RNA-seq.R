#30/07/15_Katharina.Wyschetzki@ur.de
#twofactorial:injury(injury, control) + lane


#load counts (see Htseq-count table)
qes20<-read.table('A02_qes20.txt',F,row.names=1)
qes24<-read.table('A03_qes24.txt',F,row.names=1)
qes25<-read.table('A04_qes25.txt',F,row.names=1)
qes26<-read.table('A05_qes26.txt',F,row.names=1)
qes27<-read.table('A06_qes27.txt',F,row.names=1)
qes36<-read.table('A08_qes36.txt',F,row.names=1)
qes41<-read.table('A09_qes41.txt',F,row.names=1)

colnames(qes20)<-"qes20"
colnames(qes24)<-"qes24"
colnames(qes25)<-"qes25"
colnames(qes26)<-"qes26"
colnames(qes27)<-"qes27"
colnames(qes36)<-"qes36"
colnames(qes41)<-"qes41"

qe5<-read.table('A10_qe5.txt',F,row.names=1)
qe6<-read.table('A11_qe6.txt',F,row.names=1)
qe9<-read.table('A12_qe9.txt',F,row.names=1)
qe10<-read.table('A13_qe10.txt',F,row.names=1)
qe13<-read.table('A14_qe13.txt',F,row.names=1)
qe18<-read.table('A16_qe18.txt',F,row.names=1)
qe22<-read.table('A17_qe22.txt',F,row.names=1)

colnames(qe5)<-"qe5"
colnames(qe6)<-"qe6"
colnames(qe9)<-"qe9"
colnames(qe10)<-"qe10"
colnames(qe13)<-"qe13"
colnames(qe18)<-"qe18"
colnames(qe22)<-"qe22"

qws4<-read.table('A19_qws4.txt',F,row.names=1)
qws11<-read.table('A20_qws11.txt',F,row.names=1)
qws27<-read.table('A22_qws27.txt',F,row.names=1)
qws33<-read.table('A23_qws33.txt',F,row.names=1)
qws38<-read.table('A24_qws38.txt',F,row.names=1)
qws43<-read.table('A25_qws43.txt',F,row.names=1)
qws15<-read.table('A26_qws15.txt',F,row.names=1)

colnames(qws4)<-"qws4"
colnames(qws11)<-"qws11"
colnames(qws27)<-"qws27"
colnames(qws33)<-"qws33"
colnames(qws38)<-"qws38"
colnames(qws43)<-"qws43"
colnames(qws15)<-"qws15"

qw3<-read.table('A27_qw3.txt',F,row.names=1)
qw5<-read.table('A28_qw5.txt',F,row.names=1)
qw14<-read.table('A29_qw14.txt',F,row.names=1)
qw22<-read.table('A30_qw22.txt',F,row.names=1)
qw29<-read.table('A31_qw29.txt',F,row.names=1)
qw34<-read.table('A32_qw34.txt',F,row.names=1)
qw42<-read.table('A33_qw42.txt',F,row.names=1)

colnames(qw3)<-"qw3"
colnames(qw5)<-"qw5"
colnames(qw14)<-"qw14"
colnames(qw22)<-"qw22"
colnames(qw29)<-"qw29"
colnames(qw34)<-"qw34"
colnames(qw42)<-"qw42"

#merge counts
countTable_injury<-cbind(qe5,qe6,qe9,qe10,qe13,qe18,qe22,
                         qw3,qw5,qw14,qw22,qw29,qw34,qw42,
                         qes20,qes24,qes25,qes26,qes27,qes36,qes41,
                         qws4,qws11,qws27,qws33,qws38,qws43,qws15)

Design=data.frame(row.names=colnames(countTable_injury),
                  treatment=c("control","control","control","control","control","control","control",
                              "control","control","control","control","control","control","control",
                              "injury","injury","injury","injury","injury","injury","injury",
                              "injury","injury","injury","injury","injury","injury","injury"),
                  lane=c("1","2","3","4","6","7","8",
                         "1","2","4","5","6","7","8",
                         "1","2","3","4","5","7","8",
                         "1","2","3","5","6","7","8"))

library("DESeq")
cds = newCountDataSet(countTable_injury, Design) 
detach(package:DESeq) 

#create DESeqDataSet for DESeq2
library("DESeq2")
counts <- counts(cds)
columns <- pData(cds)[c("treatment","lane")]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = columns, design = ~ lane + treatment)
dds <- DESeq(dds)


#extract result
res <- results(dds, contrast = c ("treatment","injury","control"))
res <- res[order(res$padj),]


---------------------------------------------------------------------------------------------------------------------------

#correct counts for lane effect with ComBat
library("sva")
modcombat = model.matrix(~ treatment, colData(dds))

vsd <- varianceStabilizingTransformation(dds)
vsd.counts <- assay(vsd)

vsd.combat = ComBat(dat=vsd.counts, batch=dds$lane, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


#PCA plot
vsd.combat -> assay(vsd)
plotPCA(vsd, intgroup="treatment")
data <- plotPCA(vsd, intgroup=c("treatment","lane"), returnData=TRUE)
treatment <- data$treatment
pc<-cbind(data$PC1, data$PC2)
par(pin=c(4.2,4.2))
plot(data$PC1, data$PC2, col=c("#56B4E9","#E69F00")[treatment], 
     xlab="PC1: 35% variance", ylab="PC2: 10% variance", pch=16)
library(vegan)
ordispider(pc, factor(data$treatment), label=TRUE)
ordihull(pc, factor(data$treatment), lty = "dotted")



#Gene clustering with 20 most variable genes
topVarGenes <- head(order(rowVars(assay(vsd)),decreasing=TRUE),20)
mat <- assay(vsd)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[c("treatment")])

row.names(mat) <- c("Hymenoptaecin (Cobs_04663)","Cobs_10979","yellow-g","Chorion protein b at 7F","rdgBbeta","Tep2 (Cobs_11839)",
"Dopamine transporter","CG9518 (Cobs_03171)","Cobs_12443","pale","yellow-g2","Cobs_00992","Es2","Cobs_01490",
"Pepck","CG4409","Ugt86De (Cobs_17854)","CG5618","CG6142 (Cobs_10405)","Karl (Cobs_01807)")

#colors
Treatment = c("#56B4E9","#E69F00")
names(Treatment) = c("control", "injury")
ann_colors = list(treatment = Treatment)

library(pheatmap)
pheatmap(mat, color = colorRampPalette(rev(c("#D55E00", "white", "#0072B2")))(100), 
         annotation_col=df, border_color=NA, show_colnames = F, annotation_colors=ann_colors, fontsize=10)
