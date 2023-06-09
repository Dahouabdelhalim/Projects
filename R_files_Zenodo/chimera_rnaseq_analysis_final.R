### R/bioconductor differential expression (DESeq2) workflow using counted reads (htseq-count)
## edited for publication: Oct 11, 2019

libraries <- c("DESeq2", "ggplot2", "gplots", "RColorBrewer", "gridExtra", "GenomicFeatures", "Gviz", "gage", "GOstats", "GSEABase", "ShortRead", "reshape2", "mgcv", "AICcmodavg", "GGally", "lsmeans")
lapply(libraries, require, character.only=T) 

# used for formatting figures, adjust parameters as necessary
manuscript_theme = theme_classic() + theme(axis.line.x = element_line(color="black", size = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=20))
manuscript_theme_grid = theme_classic() + theme(axis.line.x = element_line(color="black", size = .7),axis.line.y = element_line(color="black", size = .7),text=element_text(size=16))


#### QC of fastq and aligned bams ####
## summary stats on sequence quality in R using ShortRead
# quick quality assessment of sequenced reads, by default using sample of 1 mil reads per file
fls <- dir("./", "*clip.fastq$", full=TRUE)
qaSummary <- qa(fls, type="fastq")
report(qaSummary,type="html",dest="/personal/snoh/chimera2_redo/bamQAreport")
# Comments: ShortRead â€“ base call frequencies suggests a5,ab5,b5 should not be used

## summary stats on alignment from picard CollectAlignmentSummaryMetrics
bam.stat <- read.table("all.alignment.stats.txt", h=T, row.names=1, sep="\\t", na.strings ="")
summary(bam.stat)


#### Read in count data and generate metadata table ####
## full dataset including samples to be removed
counts.chi <- read.table("all.gmap.count.txt", sep="\\t", h=T, row.names=1)

sampleTable <- data.frame(shortname=c("a.3","a.4","a.5","ab.3","ab.4","ab.5","b.3","b.4","b.5","c.1","c.2","c.3","cd.1","cd.2","cd.3","d.1","d.2","d.3","e.1","e.2","e.3","ef.1","ef.2","ef.3","f.1","f.2","f.3","g.1","g.2","g.3","gh.1","gh.2","gh.3","h.1","h.2","h.3"),
                          condition=factor(c(rep(c("clonal","clonal","clonal","chimeric","chimeric","chimeric","clonal","clonal","clonal"),4)), levels=c("clonal","chimeric")),
                          rep=factor(c(rep("ab",9),rep("cd",9),rep("ef",9),rep("gh",9))),
                          batch=factor(c(rep(c("lane4","lane5","lane6"),12))),
                          clone=factor(c("a","a","a","ab","ab","ab","b","b","b","c","c","c","cd","cd","cd","d","d","d","e","e","e","ef","ef","ef","f","f","f","g","g","g","gh","gh","gh","h","h","h")))

## doublecheck outlier status of ab5 replicate with a basic DESeq2 analysis without removing these samples
dds.all <- DESeqDataSetFromMatrix(countData=counts.chi, colData=sampleTable, tidy=F, design=~batch+rep+condition)

cor.ab <- as.data.frame(cor(assay(dds.all)[,c(1:9)]))
cor.cd <- as.data.frame(cor(assay(dds.all)[,c(10:18)]))
cor.ef <- as.data.frame(cor(assay(dds.all)[,c(19:27)]))
cor.gh <- as.data.frame(cor(assay(dds.all)[,c(28:36)]))

# check how sample read counts relate to each other
sample.quality <- data.frame(in.ab = c(cor.ab["a.3","a.4"],cor.ab["b.3","b.4"]),
                             in.cd = c(mean(cor.cd["c.1","c.2"],cor.cd["c.1","c.3"],cor.cd["c.2","c.3"]),mean(cor.cd["d.1","d.2"],cor.cd["d.1","d.3"],cor.cd["d.2","d.3"])),
                             in.ef = c(mean(cor.ef["e.1","e.2"],cor.ef["e.1","e.3"],cor.ef["e.2","e.3"]),mean(cor.ef["f.1","f.2"],cor.ef["f.1","f.3"],cor.ef["f.2","f.3"])),
                             in.gh = c(mean(cor.gh["g.1","g.2"],cor.gh["g.1","g.3"],cor.gh["g.2","g.3"]),mean(cor.gh["h.1","h.2"],cor.gh["h.1","h.3"],cor.gh["h.2","h.3"])),
                             out.ab = c(mean(cor.ab["a.3","a.5"],cor.ab["a.4","a.5"]),mean(cor.ab["b.3","b.5"],cor.ab["b.4","b.5"]))
)
sample.quality

mean(unlist(sample.quality[,1:4])) # mean pairwise correlation within each pair across all included samples 0.93
mean(unlist(sample.quality[,5])) # mean pairwise correlation within ab with excluded samples 0.63
rm(dds.all, cor.ab, cor.cd, cor.ef, cor.gh, sample.quality)


#### Final differential expression analysis after ab5 removal ####
temp1 <- counts.chi[,-c(3,6,9)]
counts.chi <- temp1
temp2 <- sampleTable[-c(3,6,9),]
sampleTable <- temp2

dds.chi <- DESeqDataSetFromMatrix(countData=counts.chi, colData=sampleTable, tidy=F, design=~batch+rep+condition)
dds.chi <- estimateSizeFactors(dds.chi)
dds.chi <- DESeq(dds.chi)
res.chi <- results(dds.chi)
res.chi <- res.chi[order(res.chi$padj),]
summary(res.chi)
topGene <- rownames(subset(res.chi,padj<0.1))
topGene.05 <- rownames(subset(res.chi,padj<0.05))

# IMPORTANT doublecheck levels of chimeric vs clonal to interpret foldchange properly
mcols(res.chi)

# DE genes are slightly different if the version of DESeq2 is different. To recreate the exact analysis we have provided the original list of up- and down- regulated genes 
# export results tables
df.res.chi <- as.data.frame(res.chi)
summary(df.res.chi) #padj==NA are genes with no counts or those that don't pass the minimum baseMean level (4 reads)
write.table(df.res.chi, file=paste("chimera_deseq2_results",format(Sys.time(),"%Y%m%d"),"txt",sep="."), quote=F, sep="\\t")
write.table(rownames(subset(df.res.chi, padj<0.1 & log2FoldChange>0)), file=paste("chimera_deseq2_upreg_gene",format(Sys.time(),"%Y%m%d"),"txt",sep="."), quote=F, row.names=F, col.names=F) 
write.table(rownames(subset(df.res.chi, padj<0.1 & log2FoldChange<0)), file=paste("chimera_deseq2_downreg_gene",format(Sys.time(),"%Y%m%d"),"txt",sep="."), quote=F, row.names=F, col.names=F) 
write.table(rownames(subset(df.res.chi, padj!="NA")), file=paste("chimera_deseq2_universe_gene",format(Sys.time(),"%Y%m%d"),"txt",sep="."), quote=F, row.names=F, col.names=F) 


#### visualize results ####
# first set color palette for consistency
mypalette <- colorRampPalette(c("#e7003e","#ff8900","#5de100","#086fa1")) # http://paletton.com/#uid=75r0Q0kw0w0jyC+oRxVy4oIDfjr
hist(discoveries, col=mypalette(4))
mypalette.bw <- colorRampPalette(c("#e7003e","white"))
hist(discoveries, col=mypalette.bw(4))

# MA plot
DESeq2::plotMA(res.chi, ylim=c(-1,1))

# volcano plot
hist(discoveries, col=mypalette(3))
df.res.chi$sig <- df.res.chi$padj
df.res.chi$sig <- ifelse(df.res.chi$padj<=0.05, "FDR < 0.05","not sig.")
df.res.chi$sig <- factor(ifelse(df.res.chi$padj>0.05 & df.res.chi$padj<=0.1, "FDR < 0.1", df.res.chi$sig))

ggplot(df.res.chi, aes(x=log2FoldChange, y=-log10(pvalue), color=sig)) + geom_point(shape=1, size=1) + scale_color_manual(values=c(mypalette(3)[1],mypalette(3)[2],"grey"),name="") + manuscript_theme + geom_hline(yintercept = 3, col=mypalette(3)[3],lty=2) + ylab("-log10(P-value)")
ggplot(df.res.chi, aes(x=log2FoldChange, y=-log10(pvalue), color=sig)) + geom_point(shape=1, size=1) + scale_color_manual(values=c(mypalette(3)[1],mypalette(3)[2],"grey"),name="") + manuscript_theme + geom_hline(yintercept = 3, col=mypalette(3)[3],lty=2) + ylab("-log10(P-value)") + ylim(-1,13) + xlim(-6,6)

setEPS() # outlier is at log2FC=-2.19, -log10pvalue=70.7
postscript(file=paste("volcano_plot",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=8, height=6)
ggplot(df.res.chi, aes(x=log2FoldChange, y=-log10(pvalue), color=sig)) + geom_point(shape=1, size=1) + scale_color_manual(values=c(mypalette(3)[1],mypalette(3)[2],"grey"),name="") + manuscript_theme + geom_hline(yintercept = 3, col=mypalette(3)[3],lty=2) + ylab("-log10(P-value)") + ylim(-1,13) + xlim(-6,6)
dev.off()

# count comparison of individual DE genes mentioned in paper
normalized.data <- data.frame(plotCounts(dds.chi, gene=topGene[1], intgroup=c("condition","rep"), returnData=TRUE)[,1])
for (i in 2:length(topGene)) {
  normalized.data[,i] <- plotCounts(dds.chi, gene=topGene[i], intgroup=c("condition","rep"), returnData=TRUE)[,1]
}
colnames(normalized.data) <- topGene
temp <- cbind(sampleTable, normalized.data)
normalized.data <- temp

ggplot(go.fc.fig.liter, aes(x=gene_set, y=log2FC, fill=gene_set)) + theme_classic() + geom_dotplot(binaxis='y', stackdir='center') + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="blue") + scale_fill_manual(values=c(sm.palette(20)), guide=F) 

gene.a <- ggplot(normalized.data, aes(x=rep, y=DDB_G0289121, fill=condition)) + geom_boxplot() + scale_fill_manual(values=c(rev(mypalette.bw(2)))) + manuscript_theme + ggtitle("(a) pde4") + xlab("") + ylab("") 
gene.b <- ggplot(normalized.data, aes(x=rep, y=DDB_G0276893, fill=condition)) + geom_boxplot() + scale_fill_manual(values=c(rev(mypalette.bw(2)))) + manuscript_theme + ggtitle("(b) ctxB") + xlab("") + ylab("Normalized counts") 
gene.c <- ggplot(normalized.data, aes(x=rep, y=DDB_G0288179, fill=condition)) + geom_boxplot() + scale_fill_manual(values=c(rev(mypalette.bw(2)))) + manuscript_theme + ggtitle("(c) carB") + xlab("Replicate pair") + ylab("") 

grid.arrange(gene.a, gene.b, gene.c)

setEPS()
postscript(file=paste("genes_counts_by_condition",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=8, height=8)
grid.arrange(gene.a, gene.b, gene.c)
dev.off()

# transform data (used in PCA plots and heatmap)
rld.chi <- rlogTransformation(dds.chi)
plotPCA(rld.chi, intgroup = "condition", ntop = 500, returnData = FALSE)

# pca plot with ellipses
plotPCA.ellipse <-  function(object, intgroup="condition", ellipsegroup=NA, ntop=500, returnData=FALSE, pcs = c(1,2))
{
  # calculate the variance for each gene
  rv <- matrixStats::rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  percentVar <- round(percentVar, 3) * 100
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=intgroup)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  # specify color palette
  mypalette <- colorRampPalette(c("#e7003e","#ff8900","#5de100","#086fa1"))
  # find centroid of ellipse
  cent <- data.frame(do.call(rbind, lapply(split(d[,1:2], d[,3]), colMeans)))
  cent$grp <- factor(rownames(cent))
  # graph to output
  ggplot(data=d, aes(x=PC1, y=PC2, group=group, color=group)) + geom_point(shape=16, size=2, show.legend=F) + scale_color_manual(values=c(mypalette(12)),guide=F) + geom_point(data=cent, aes(group=grp, color=grp), shape=3, size=5, show.legend=F) + geom_text(data=cent, aes(label=grp, group=grp, color=grp), hjust = 0, nudge_x = 1, nudge_y = 3, size=8, fontface="italic", show.legend=F) + manuscript_theme + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + stat_ellipse(aes(group=ellipsegroup, lty=ellipsegroup))

}

plotPCA.ellipse(rld.chi, intgroup=sampleTable$clone, ellipsegroup=sampleTable$condition)
#plotPCA.ellipse(rld.chi, intgroup=sampleTable$clone, ellipsegroup=sampleTable$rep)

setEPS()
postscript(file=paste("pcplot_clone_labeled",format(Sys.time(),"%Y%m%d"),"eps",sep="."), width=8, height=6)
plotPCA.ellipse(rld.chi, intgroup=sampleTable$clone, ellipsegroup=sampleTable$condition)
dev.off()

# specify divergent palette for heatmaps
sm.palette <- colorRampPalette(c("black","white","#e7003e"))
sidecols <- c("white","#e7003e")[rld.chi$condition]
hist(discoveries, col=sm.palette(4))

# heatmap comparing samples, with top 500 variable genes
top500 <- rownames(res.chi[1:500,]) 
sampleDists500 <- dist(t(assay(rld.chi[top500,])))
dismatrix500 <- as.matrix(sampleDists500)
colnames(dismatrix500) <- colData(rld.chi)$shortname
rownames(dismatrix500) <- colData(rld.chi)$shortname
hc500 <- hclust(sampleDists500)

dev.off()
heatmap.2(dismatrix500, Rowv=as.dendrogram(hc500), col=sm.palette, dendrogram="row", symm=F, RowSideColors=sidecols, trace="none", cexRow=1, cexCol=1, keysize=1, key.xlab="Distance", key.ylab="Count", key.title=NA, sepwidth=c(0.001,0.001), sepcolor="gray80", colsep=1:ncol(dismatrix500), rowsep=1:nrow(dismatrix500))
par(lend = 1)  
legend("topleft", legend = c("clonal", "chimera"), col = c("#f5f5f5","#e7003e"), lty= 1, lwd = 10, bty = "n", cex=.8, horiz=TRUE, title="Sample condition")            

# heatmaps with top DE genes, overall and centered(d)
mat.gene <- assay(rld.chi[topGene,])
mat.gened <- mat.gene - rowMeans(mat.gene)
colnames(mat.gened) <- colData(rld.chi)$shortname

# dendrograms ordered by branch mean (default is branch sum), with same hclust as samples above
dev.off()    
heatmap.2(mat.gened, Colv=as.dendrogram(hc500), col=sm.palette, trace="none", dendrogram="row", ColSideColors=sidecols, cexRow=1, cexCol=1, keysize=1, key.xlab="Gene expression", key.ylab="Count", key.title=NA, reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), sepwidth=c(0.001,0.001), sepcolor="gray90", colsep=1:ncol(mat.gened), rowsep=1:nrow(mat.gened))
par(lend = 1)  
legend("topleft", legend = c("clonal", "chimera"), col = c("#f5f5f5","#e7003e"), lty= 1, lwd = 10, bty = "n", cex=.8, horiz=TRUE, title="Sample condition")            

# dendrograms ordered by branch mean (default is branch sum), with samples sorted by clonal vs. chimeric above
so <- c("a.3","a.4","b.3","b.4","c.1","c.2","c.3","d.1","d.2","d.3","e.1","e.2","e.3","f.1","f.2","f.3","g.1","g.2","g.3","h.1","h.2","h.3","ab.3","ab.4","cd.1","cd.2","cd.3","ef.1","ef.2","ef.3","gh.1","gh.2","gh.3")
dev.off()    
heatmap.2(mat.gened[,so], Colv=F, col=sm.palette, trace="none", dendrogram="row", ColSideColors=c(rep("white",22),rep("#e7003e",11)), cexRow=1, cexCol=1, margins=c(4,4), keysize=1, key.xlab="Gene expression", key.ylab="Count", key.title=NA, reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), sepwidth=c(0.001,0.001), sepcolor="gray90", colsep=1:ncol(mat.gened), rowsep=1:nrow(mat.gened))
par(lend = 1)  
legend("topleft", legend = c("clonal", "chimera"), col = c("#f5f5f5", "#e7003e"), lty= 1, lwd = 10, bty = "n", cex=.8, horiz=TRUE, title="Sample condition")            
dev.off()    

# export eps file of heatmaps
setEPS(width=9, height=9) # default is 7x7
postscript(file=paste("chimera_heatmap_top_genes",format(Sys.time(),"%Y%m%d"),"eps",sep="."))
heatmap.2(mat.gened[,so], Colv=F, trace="none", col=sm.palette, dendrogram="row", ColSideColors=c(rep("#f5f5f5",22),rep("#e7003e",11)), cexRow=0.6, cexCol=1, margins=c(4,6), keysize=1, key.xlab="Gene expression", key.ylab="Count", key.title=NA, reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), sepwidth=c(0.001,0.001), sepcolor="gray90", colsep=1:ncol(mat.gened), rowsep=1:nrow(mat.gened))
par(lend = 1)  
legend("topleft", legend = c("clonal", "chimera"), col = c("#f5f5f5", "#e7003e"), lty= 1, lwd = 10, bty = "n", cex=.8, horiz=TRUE, title="Sample condition")            
dev.off()


#### go term analysis with GOstats and KEGG, and GO visualization ####
## GOstats analysis
# read in go annotations
gocurateData <- read.table("gene_association.dictyBase.20180701.reformat.filter.txt",h=F)
goCurate=GOFrame(gocurateData,organism="Dictyostelium discoideum")
goCurFrame=GOAllFrame(goCurate)
gsc.c <- GeneSetCollection(goCurFrame,setType=GOCollection())

universe <- rownames(subset(res.chi, baseMean>=4))
universe.c <- universe[universe %in% gocurateData[,3]]

downgene <- rownames(subset(res.chi, padj<0.1 & log2FoldChange<0))
upgene <- rownames(subset(res.chi, padj<0.1 & log2FoldChange>0)) # too few to do enrichment analyses?
downgene.c <- downgene[downgene %in% gocurateData[,3]]
upgene.c <- upgene[upgene %in% gocurateData[,3]]

# up and down genes separately
bp.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc.c, geneIds=downgene.c, universeGeneIds=universe.c, ontology="BP",pvalueCutoff=0.05,conditional=T,testDirection="over")
bp.down.c <- hyperGTest(bp.params)
rm(bp.params)
mf.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc.c, geneIds=downgene.c, universeGeneIds=universe.c, ontology="MF",pvalueCutoff=0.05,conditional=T,testDirection="over")
mf.down.c <- hyperGTest(mf.params)
rm(mf.params)
bp.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc.c, geneIds=upgene.c, universeGeneIds=universe.c, ontology="BP",pvalueCutoff=0.05,conditional=T,testDirection="over")
bp.up.c <- hyperGTest(bp.params)
rm(bp.params)
mf.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc.c, geneIds=upgene.c, universeGeneIds=universe.c, ontology="MF",pvalueCutoff=0.05,conditional=T,testDirection="over")
mf.up.c <- hyperGTest(mf.params)
rm(mf.params)

summary(bp.down.c)
summary(mf.down.c)

summary(bp.up.c)
summary(mf.up.c)

# format GOSTAT dataframes for tables
temp <- subset(as.data.frame(summary(bp.down.c)), Size>8 & Size<400)[,c(7,1,4,5,6,2)]
df.bp.down.cur <- cbind(temp[,1:2], round(temp[,3:6], digits=3))
temp <- subset(as.data.frame(summary(mf.down.c)), Size>8 & Size<400)[,c(7,1,4,5,6,2)]
df.mf.down.cur <- cbind(temp[,1:2], round(temp[,3:6], digits=3))
temp <- subset(as.data.frame(summary(bp.up.c)), Size>8 & Size<400)[,c(7,1,4,5,6,2)]
df.bp.up.cur <- cbind(temp[,1:2], round(temp[,3:6], digits=3))
temp <- subset(as.data.frame(summary(mf.up.c)), Size>8 & Size<400)[,c(7,1,4,5,6,2)]
df.mf.up.cur <- cbind(temp[,1:2], round(temp[,3:6], digits=3))
#rm(temp,bp.over.c, mf.over.c, bp.down.c, mf.down.c, bp.up.c, mf.up.c)

# these data can be accessed as below
methods(class="GOHyperGResult")
geneIdsByCategory(bp.over.c)["GO:0006260"]
#$`GO:0006260`
#[1] "DDB_G0272760" "DDB_G0292958"

# get genes in each category: x is the data frame, y is the hypergtest object
get.members2 <- function(x,y) {
  for(i in as.numeric(rownames(x))) {
    print(geneIdsByCategory(y)[i])
  }
}

sink('DE_down_genes_GOBP_curated.txt')
df.bp.down.cur
get.members2(df.bp.down.cur, bp.down.c)
sink()

sink('DE_down_genes_GOMF_curated.txt')
df.mf.down.cur
get.members2(df.mf.down.cur, mf.down.c) 
sink()

sink('DE_up_genes_GOBP_curated.txt')
df.bp.up.cur
get.members2(df.bp.up.cur, bp.up.c)
sink()

sink('DE_up_genes_GOMF_curated.txt')
df.mf.up.cur
get.members2(df.mf.up.cur, mf.up.c)
sink()

## compare to KEGG
kg.ddi <- gage::kegg.gsets(species="ddi")

# manual hypergeometric test
# order: (gene list in target list, universe genes in target list, universe genes not in target list, gene list size) 
hg.test <- function(a,b,c,d) { min(1-cumsum(dhyper(0:(a-1),b,c,d))) }

#kg.ddi$kg.sets 
mylist <- downgene[downgene %in% unlist(kg.ddi)]
checklist <- universe[universe %in% unlist(kg.ddi)]
hg.1 <- do.call(rbind, lapply(kg.ddi$kg.sets, function(x) sum(x %in% mylist)))
hg.2 <- do.call(rbind, lapply(kg.ddi$kg.sets, function(x) sum(x %in% checklist)))
hg.3 <- length(checklist) - hg.2
hg.kg <- as.data.frame(cbind(hg.1, hg.2, hg.3))
hg.kg.nozero <- subset(hg.kg, V1!=0)
hg.kg.nozero$hg.4 <- rep(length(mylist), nrow(hg.kg.nozero))
hg.kg.nozero$raw.p <- apply(hg.kg.nozero, 1, function(x) hg.test(x[1],x[2],x[3],x[4]) )
hg.kg.nozero$adj.p <- p.adjust(hg.kg.nozero$raw.p, method="BH")
hg.kg.nozero$exp <- hg.kg.nozero[,2]*length(mylist)/length(checklist)
temp <- subset(hg.kg.nozero, adj.p<0.05)[,c(7,1,2,6)]
df.kegg.down <- round(temp, digits=3)
colnames(df.kegg.down) <- c("ExpCount","ObsCount","SizeTerm","P-value")
df.kegg.down

mylist <- upgene[upgene %in% unlist(kg.ddi)]
checklist <- universe[universe %in% unlist(kg.ddi)]
hg.1 <- do.call(rbind, lapply(kg.ddi$kg.sets, function(x) sum(x %in% mylist)))
hg.2 <- do.call(rbind, lapply(kg.ddi$kg.sets, function(x) sum(x %in% checklist)))
hg.3 <- length(checklist) - hg.2
hg.kg <- as.data.frame(cbind(hg.1, hg.2, hg.3))
hg.kg.nozero <- subset(hg.kg, V1!=0)
hg.kg.nozero$hg.4 <- rep(length(mylist), nrow(hg.kg.nozero))
hg.kg.nozero$raw.p <- apply(hg.kg.nozero, 1, function(x) hg.test(x[1],x[2],x[3],x[4]) )
hg.kg.nozero$adj.p <- p.adjust(hg.kg.nozero$raw.p, method="BH")
hg.kg.nozero$exp <- hg.kg.nozero[,2]*length(mylist)/length(checklist)
temp <- subset(hg.kg.nozero, adj.p<0.05)[,c(7,1,2,6)]
df.kegg.up <- round(temp, digits=3)
colnames(df.kegg.up) <- c("ExpCount","ObsCount","SizeTerm","P-value")
df.kegg.up

sink('DE_downgenes_KEGG.txt')
df.kegg.down
downgene[downgene %in% kg.ddi$kg.sets$"ddi03030 DNA replication"]
sink()


## GO visualization 
source("gobp_terms/gobp_terms_and_members_for_R.R")
vars <- read.csv("gobp_terms/gobp_variable_list.txt")
names(vars)

# show GO terms with descriptions, not numbers
dict <- read.table("gobp_terms/gobp_names.txt", h=F)
dict$term <- seq(1,length(vars))
names(dict) <- c("gene_set","description","term")

# get the log2FC for each of the genes in these terms
go.fc.list = list()
for(i in 1:length(vars)) {
  index <- rownames(df.res.chi) %in% get(names(vars[i]))
  go.fc.list[[i]] <- df.res.chi[index,][,2]
}

# make long form dataframe that combines observed log2FC and go term
go.fc.len <- do.call(rbind, lapply(go.fc.list, length))
go.fc.rep <- character()
for(i in 1:length(vars)) {
  go.fc.rep <- c(go.fc.rep, rep(i, go.fc.len[i]))
}
go.fc.for <- ldply(go.fc.list, data.frame)
go.fc.fig <- cbind(go.fc.for, go.fc.rep)
names(go.fc.fig) <- c("log2FC", "term") 
rm(go.fc.len, go.fc.rep, go.fc.for)
go.fc.fig$gene_set <- dict[match(go.fc.fig$term, dict$term),2]

# order all go terms (BP) and their FC means, but this is impossible to look at
go.fc.agg <- aggregate(go.fc.fig[1], go.fc.fig[3], mean, na.rm=T) 
go.fc.fig$gene_set <- factor(go.fc.fig$gene_set, levels=as.matrix(go.fc.agg[with(go.fc.agg, order(log2FC)),][1]))
#ggplot(go.fc.fig, aes(x=gene_set, y=log2FC)) + geom_boxplot() + theme(axis.text.x=element_text(angle=45,hjust=1))  +xlab("Gene set") 

# show with top and bottom 10 DE GO terms
go.fc.fig.liter <- subset(go.fc.fig,term %in% c(levels(go.fc.fig$term)[1:10], levels(go.fc.fig$term)[186:195]))
temp <- droplevels(go.fc.fig.liter)
go.fc.fig.liter <- temp
levels(go.fc.fig.liter$gene_set)

# version with dotplots across
ggplot(go.fc.fig.liter, aes(x=gene_set, y=log2FC, fill=gene_set)) + theme_classic() + geom_dotplot(binaxis='y', stackdir='center') + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="blue") + scale_fill_manual(values=c(sm.palette(20)), guide=F) + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position = "none", plot.margin=margin(t = .2, r = .2, b = .2, l = 1.2, unit = "cm")) + xlab("GO Biological Process") + geom_vline(xintercept=c(10.5), color="darkgrey")  

setEPS()
postscript(file=paste("chimera_gene_set_with_GO_dotplot",format(Sys.time(),"%Y%m%d"),"eps",sep="."), onefile=F, width=8, height=5)
ggplot(go.fc.fig.liter, aes(x=gene_set, y=log2FC, fill=gene_set)) + theme_classic() + geom_dotplot(binaxis='y', stackdir='center') + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="blue") + scale_fill_manual(values=c(sm.palette(20)), guide=F) + theme(axis.text.x=element_text(angle=45,hjust=1), legend.position = "none", plot.margin=margin(t = .2, r = .2, b = .2, l = 1.2, unit = "cm")) + xlab("GO Biological Process") + geom_vline(xintercept=c(10.5), color="darkgrey")
dev.off()

# version with boxplots stacked
ggplot(go.fc.fig.liter, aes(x=gene_set, y=log2FC, fill=gene_set)) + theme_classic() + geom_boxplot() + scale_fill_manual(values=c(sm.palette(20)), guide=F) + xlab("GO Biological Process") + coord_flip() + geom_vline(xintercept=c(10.5), color="darkgrey") 

setEPS()
postscript(file=paste("chimera_gene_set_with_GO_boxplot",format(Sys.time(),"%Y%m%d"),"eps",sep="."), onefile=F, width=8, height=5)
ggplot(go.fc.fig.liter, aes(x=gene_set, y=log2FC, fill=gene_set)) + theme_classic() + geom_boxplot() + scale_fill_manual(values=c(sm.palette(20)), guide=F) + xlab("GO Biological Process") + coord_flip() + geom_vline(xintercept=c(10.5), color="darkgrey") 
dev.off()


#### comparison with other social genes ####
## Are Santorelli/Ostrowski REMI cheater genes different in chimeric expression?
cheater <- scan("santorelli_intragenic_ostrowski.txt", what="", sep=",")
cheat.index <- rownames(df.res.chi) %in% cheater
cheat.exp <- df.res.chi[cheat.index,][,2]
# how many of these genes are in our dataset
table(!is.na(cheat.exp)) #84

nperm <- 10000
GE.cheat = rep(0, nperm)
for(i in 1:nperm) {
  x <- sample(df.res.chi[,2], length(cheat.exp), FALSE)
  GE.cheat[i] <- mean(x, na.rm=T)
}
GE.cheat <- as.data.frame(GE.cheat)
colnames(GE.cheat)[1] <- "log2FC"
# probability that a random set of genes have lower expression than REMI cheaters 
sum(unlist(GE.cheat)<=mean(cheat.exp, na.rm=T))/10000  # 0.26

ggplot(GE.cheat, aes(x=log2FC)) + geom_histogram() + manuscript_theme + geom_vline(xintercept=mean(cheat.exp, na.rm=T), colour="red", linetype="longdash")

##Are Hirose genes different in chimeric expression?
hirose.alld <- scan("hirose_down.txt", what="", sep="\\n")
hirose.allu <- scan("hirose_up.txt", what="", sep="\\n")

hirose.d <- scan("hirose_12down.txt", what="", sep="\\n")
hirose.u <- scan("hirose_12up.txt", what="", sep="\\n")

hup.index <- rownames(df.res.chi) %in% hirose.u
hup.exp <- df.res.chi[hup.index,][,2]
hdown.index <- rownames(df.res.chi) %in% hirose.d
hdown.exp <- df.res.chi[hdown.index,][,2]

table(!is.na(hup.exp)) #13
table(!is.na(hdown.exp)) #71

GE.hup = rep(0, nperm)
for(i in 1:nperm) {
  x <- sample(df.res.chi[,2], length(hup.exp), FALSE)
  GE.hup[i] <- mean(x, na.rm=T)
}
GE.hup <- as.data.frame(GE.hup)
colnames(GE.hup)[1] <- "log2FC"
# probability that a random set of genes have higher expression than hirose-up
sum(unlist(GE.hup)>=mean(hup.exp, na.rm=T))/10000 # 0.03

ggplot(GE.hup, aes(x=log2FC)) + geom_histogram() + manuscript_theme + geom_vline(xintercept=mean(hup.exp, na.rm=T), colour="red", linetype="longdash")


GE.hdown = rep(0, nperm)
for(i in 1:nperm) {
  x <- sample(df.res.chi[,2], length(hdown.exp), FALSE)
  GE.hdown[i] <- mean(x, na.rm=T)
}
GE.hdown <- as.data.frame(GE.hdown)
colnames(GE.hdown)[1] <- "log2FC"
# probability that a random set of genes have lower expression than hirose-down
sum(unlist(GE.hdown)<=mean(hdown.exp, na.rm=T))/10000 # < 0.0001

ggplot(GE.hdown, aes(x=log2FC)) + geom_histogram() + manuscript_theme + geom_vline(xintercept=mean(hdown.exp, na.rm=T), colour="red", linetype="longdash")


# observed
temp1 <- data.frame(log2FC=subset(df.res.chi, rownames(df.res.chi) %in% upgene)[,2], group=rep("chimera_up", length(upgene)))
temp2 <- data.frame(log2FC=subset(df.res.chi, rownames(df.res.chi) %in% downgene)[,2], group=rep("chimera_down", length(downgene)))
temp3 <- data.frame(log2FC=cheat.exp, group=rep("cheater", length(cheat.exp)))
temp4 <- data.frame(log2FC=hup.exp, group=rep("hirose_up_12hr", length(hup.exp)))
temp5 <- data.frame(log2FC=hdown.exp, group=rep("hirose_down_12hr", length(hdown.exp)))
fc.all <- rbind(temp1, temp2, temp3, temp4, temp5)
hist(discoveries, col=sm.palette(5))

ggplot(fc.all, aes(x=group, y=log2FC, fill=group)) + geom_violin(trim = F) + manuscript_theme + geom_boxplot(width=.1, fill="white") + scale_fill_manual(values=c(sm.palette(5)[5],sm.palette(5)[1],sm.palette(5)[3],sm.palette(5)[4],sm.palette(5)[2]),guide=F) + theme(axis.text.x=element_text(angle=45,hjust=1)) + geom_hline(yintercept=0, lty=2, color="grey") + geom_vline(xintercept=c(2.5,3.5), color="darkgrey") + xlab("Gene set") + scale_x_discrete("",labels=c("chimera_up"="chimera-biased (up)","chimera_down"="chimera-biased (down)","cheater"="cheater","hirose_up_12hr"="12-hr (up)","hirose_down_12hr"="12-hr (down)")) #+ ylab("Relative gene expression (log2FC)") 

setEPS()
postscript(file=paste("candidate gene_set_mean_expr",format(Sys.time(),"%Y%m%d"),"eps",sep="."), onefile=F, width=8, height=5)
ggplot(fc.all, aes(x=group, y=log2FC, fill=group)) + geom_violin(trim = F) + manuscript_theme + geom_boxplot(width=.1, fill="white") + scale_fill_manual(values=c(sm.palette(5)[5],sm.palette(5)[1],sm.palette(5)[3],sm.palette(5)[4],sm.palette(5)[2]),guide=F) + theme(axis.text.x=element_text(angle=45,hjust=1)) + geom_hline(yintercept=0, lty=2, color="grey") + geom_vline(xintercept=c(2.5,3.5), color="darkgrey") + xlab("Gene set") + scale_x_discrete("",labels=c("chimera_up"="chimera-biased (up)","chimera_down"="chimera-biased (down)","cheater"="cheater","hirose_up_12hr"="12-hr (up)","hirose_down_12hr"="12-hr (down)")) #+ ylab("Relative gene expression (log2FC)") 
dev.off()


#### what can we tell about the timing of chimeric development ####
counts.chi <- read.table("all.gmap.count.txt", sep="\\t", h=T, row.names=1)

sampleTable.time <- data.frame(shortname=c("rep1.00","rep1.01","rep1.02","rep1.03","rep1.04","rep1.05","rep1.06","rep1.07","rep1.08","rep1.09","rep1.10","rep1.11","rep1.12","rep1.14","rep1.16","rep1.18","rep1.20","rep1.22","rep1.24","rep2.00","rep2.01","rep2.02","rep2.03","rep2.04","rep2.05","rep2.06","rep2.07","rep2.08","rep2.09","rep2.10","rep2.11","rep2.12","rep2.14","rep2.16","rep2.18","rep2.20","rep2.22","rep2.24"),
                               condition=c("00hr","01hr","02hr","03hr","04hr","05hr","06hr","07hr","08hr","09hr","10hr","11hr","12hr","14hr","16hr","18hr","20hr","22hr","24hr","00hr","01hr","02hr","03hr","04hr","05hr","06hr","07hr","08hr","09hr","10hr","11hr","12hr","14hr","16hr","18hr","20hr","22hr","24hr"),
                               rep=factor(c(rep("rep1",19),rep("rep2",19))),
                               batch=factor(c(rep("lane1",19),rep("lane2",19))),
                               clone=factor(c(rep("ax4",38))))

sampleTable.long <- rbind(sampleTable, sampleTable.time)
sampleTable.long <- droplevels(sampleTable.long)
sampleTable.long$plot <- c(rep("wild",33),rep(c(rep("ax4",7),"07hr","08hr","09hr","10hr","11hr","12hr","14hr","16hr",rep("ax4",4)),2))
sampleTable.long$sample <- factor(c(rep("wild",33),rep("ax4",38)), levels=c("ax4","wild"))
row.names(sampleTable.long) <- seq(1:71)

counts.time <- read.table("rosengarten.gmap.count.txt", sep="\\t", h=T, row.names=1)
counts.ref <- cbind(counts.chi, counts.time)
  
# can only fit one factor, and batch seems to be the most nuisance 
dds.ref <- DESeqDataSetFromMatrix(countData=counts.ref, colData=sampleTable.long, tidy=F, design=~batch)
dds.ref <- estimateSizeFactors(dds.ref)
dds.ref <- DESeq(dds.ref)

res.ref <- results(dds.ref)
res.ref <- res.ref[order(res.ref$padj),]

rld.ref <- rlogTransformation(dds.ref, blind=F)
colnames(rld.ref) <- colnames(counts.ref)
plotPCA(rld.ref, intgroup="plot")
rld.ref.data <- plotPCA(rld.ref, intgroup=c("condition","rep","plot"), returnData=TRUE)
rld.ref.data$plot2 <- factor(c(rep("wild",33),rep(c(rep("earlier",7),"07 hr","08 hr","09 hr","10 hr","11 hr","12 hr","14 hr","16 hr",rep("later",4)),2)), levels=c("earlier","07 hr","08 hr","09 hr","10 hr","11 hr","12 hr","14 hr","16 hr","later","wild"))


setEPS()
postscript(file=paste("pca_wild_ax4_all_genes",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=6)
ggplot(rld.ref.data, aes(x=PC1, y=PC2, shape=plot2, fill=plot2)) + geom_point(size=4) + manuscript_theme + xlab("PC1: 54% variance") + ylab("PC2: 23% variance") + scale_fill_manual(values=c("#F8F8FF","#D9D9DF","#BABABF","#9B9B9F","#7C7C7F","#3399FF","#5D5D5F","#3E3E3F","#1F1F1F","#000000","#DD1153")) + scale_shape_manual(values=c(21, 22, 24, 23, 25, 22, 24, 23, 25, 21, 21)) +theme(legend.title=element_blank(), legend.text = element_text(size=20), axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.x = element_text(size=20),axis.text.y = element_text(size=20))
dev.off()


sampleDists.ref <- dist(t(assay(rld.ref)))
m <- data.frame(t(combn(rownames(t(assay(rld.ref))),2)), as.numeric(sampleDists.ref))
names(m) <- c("c1", "c2", "distance")

get_timeseries_distance_annotated <- function(m) {
  dist.a <- subset(m, grepl("^a\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.a$hr <- substr(dist.a$c2,6,7)
  
  dist.b <- subset(m, grepl("^b\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.b$hr <- substr(dist.b$c2,6,7)
  
  dist.ab <- subset(m, grepl("^ab\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.ab$hr <- substr(dist.ab$c2,6,7)
  
  dist.c <- subset(m, grepl("^c\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.c$hr <- substr(dist.c$c2,6,7)
  
  dist.d <- subset(m, grepl("^d\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.d$hr <- substr(dist.d$c2,6,7)
  
  dist.cd <- subset(m, grepl("^cd\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.cd$hr <- substr(dist.cd$c2,6,7)
  
  dist.e <- subset(m, grepl("^e\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.e$hr <- substr(dist.e$c2,6,7)
  
  dist.f <- subset(m, grepl("^f\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.f$hr <- substr(dist.f$c2,6,7)
  
  dist.ef <- subset(m, grepl("^ef\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.ef$hr <- substr(dist.ef$c2,6,7)
  
  dist.g <- subset(m, grepl("^g\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.g$hr <- substr(dist.g$c2,6,7)
  
  dist.h <- subset(m, grepl("^h\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.h$hr <- substr(dist.h$c2,6,7)
  
  dist.gh <- subset(m, grepl("^gh\\\\.", m$c1) & grepl("^rep", m$c2))
  dist.gh$hr <- substr(dist.gh$c2,6,7)
  
  dist.all <- rbind(dist.a, dist.b, dist.ab, dist.c, dist.d, dist.cd, dist.e, dist.f, dist.ef, dist.g, dist.h, dist.gh)
  dist.all$condition <- factor( c(rep("clonal",nrow(dist.a)+nrow(dist.b)), rep("chimeric",nrow(dist.ab)), rep("clonal",nrow(dist.c)+nrow(dist.d)), rep("chimeric",nrow(dist.cd)), rep("clonal",nrow(dist.e)+nrow(dist.f)), rep("chimeric",nrow(dist.ef)), rep("clonal",nrow(dist.g)+nrow(dist.h)), rep("chimeric",nrow(dist.gh)) ), levels=c("clonal","chimeric"))
  dist.all$rep <- factor(c(rep("ab",nrow(dist.a)+nrow(dist.b)+nrow(dist.ab)), rep("cd",nrow(dist.c)+nrow(dist.d)+nrow(dist.cd)), rep("ef",nrow(dist.e)+nrow(dist.f)+nrow(dist.ef)), rep("gh",nrow(dist.g)+nrow(dist.h)+nrow(dist.gh))))
  dist.all$hr <- as.factor(dist.all$hr)
  
  return(dist.all)
}

## distance matrix using all genes
dist.all <- get_timeseries_distance_annotated(m)

## alternate distance matrix using only genes related to development 
# rose.up and rose.down gene lists were previously read in
rose.up <- gage::readList("rosengarten_filter_up_by_hour.gmt")
rose.down <- gage::readList("rosengarten_filter_down_by_hour.gmt")
timing.genes <- unique(c(unlist(rose.up), unlist(rose.down)))

sub.assay <- assay(rld.ref)[rownames(assay(rld.ref)) %in% timing.genes,]
sampleDists.sub <- dist(t(sub.assay), method="euclidean")

m2 <- data.frame(t(combn(rownames(t(sub.assay)),2)), as.numeric(sampleDists.sub))
names(m2) <- c("c1", "c2", "distance")
dist.sub <- get_timeseries_distance_annotated(m2)

# quick visualization to compare these two datasets
dist.all.fig <- aggregate(distance~rep + hr + condition, data=dist.all, FUN=mean)
dist.sub.fig <- aggregate(distance~rep + hr + condition, data=dist.sub, FUN=mean)

temp.all <- ggplot(dist.all.fig, aes(x=hr, y=distance, group=condition, fill=condition)) + geom_point(shape=21) + geom_smooth(method="loess", se=F, aes(linetype=condition)) + manuscript_theme + scale_fill_manual(values = c("black","grey")) + xlab("gene expression at each hour from starvation") + ggtitle("(a) All genes") + theme(plot.title=element_text(hjust=0))
temp.sub <- ggplot(dist.sub.fig, aes(x=hr, y=distance, group=condition, fill=condition)) + geom_point(shape=21) + geom_smooth(method="loess", se=F, aes(linetype=condition)) + manuscript_theme + scale_fill_manual(values = c("black","grey")) + xlab("gene expression at each hour from starvation") + ggtitle("(b) Timing genes") + theme(plot.title=element_text(hjust=0))
grid.arrange(temp.all, temp.sub)
# difference between these two types of distances is not great but present


## use gam and derivatives to characterize interaction between condition and distance by hour
dum <- as.vector(rep(1,1254)) # dummy variable to extract random effects later
gam.dist.all.full <- gam(data=dist.all, distance ~ condition + s(as.numeric(hr), by=condition) + s(rep, bs="re", by=dum), method="REML")
gam.dist.all.null <- gam(data=dist.all, distance ~ condition + s(as.numeric(hr)) + s(rep, bs="re", by=dum), method="REML")
AICcmodavg::AICc(gam.dist.all.null)
AICcmodavg::AICc(gam.dist.all.full)
AICcmodavg::AICc(gam.dist.all.full) - AICcmodavg::AICc(gam.dist.all.null) 
anova(gam.dist.all.null, gam.dist.all.full, test="Chisq") # little evidence to suggest condition-specific curves
rm(gam.dist.all.full)

# Approximate the derivatives using finite differences. This is always possible but less accurate than the other methods
# predict for all genes 
pdat.all <- cbind(dist.all[,c(4:6)],dum)
pdat.all$hr <- as.numeric(pdat.all$hr)

step <- 0.1 # step size for finite difference http://stackoverflow.com/questions/14207250/determining-derivatives-from-gam-smooth-object 
pdat.all2 <- pdat.all
pdat.all2$hr <- pdat.all$hr + step
pdat.all3 <- pdat.all
pdat.all3$hr <- pdat.all$hr - step

x0 <- predict(gam.dist.all.null, newdata=pdat.all, type='lpmatrix')
x1 <- predict(gam.dist.all.null, newdata=pdat.all2, type='lpmatrix')
xp <- (x0 - x1)/step
pdat.all$fd1 <- as.vector(xp %*% coef(gam.dist.all.null)) # 1st deriv - slope of tangent line. When 0, the curve is potentially turning (local minimum is what I want so 2nd derivative must be positive, or zero)
x2 <- predict(gam.dist.all.null, newdata=pdat.all3, type='lpmatrix')
xpp <- (x2+x1-2*x0)/step^2
pdat.all$fd2 <- as.vector(xpp %*% coef(gam.dist.all.null)) # 2nd deriv - rate of change of tangent line.  
rm(x0,x1,xp,x2,xpp)

p.ci <- predict(gam.dist.all.null, newdata=pdat.all, type="link", se.fit=TRUE)
pdat.all$upr <- as.vector(gam.dist.all.null$family$linkinv(p.ci$fit + (2*p.ci$se.fit)))
pdat.all$lwr <- as.vector(gam.dist.all.null$family$linkinv(p.ci$fit - (2*p.ci$se.fit)))
#pdat.all$fitted <- as.vector(p1.null)
pdat.all$fitted <- fitted(gam.dist.all.null)
pdat.all$resid <- residuals(gam.dist.all.null, type = "response")
pdat.all$observed.y <- napredict(gam.dist.all.null$na.action, gam.dist.all.null$y)

# inspect the gam fit
pdat.all.fig <- aggregate(pdat.all$fitted, by=list(pdat.all$hr, pdat.all$condition), FUN=mean)
names(pdat.all.fig) <- c("hr","condition","fitted")
ggplot() + geom_boxplot(data=dist.all.fig, aes(x=as.numeric(hr), y=distance, group=condition:as.factor(hr), fill=condition))  + geom_line(data=pdat.all.fig, aes(x=hr, y=fitted, group=condition, linetype=condition)) + manuscript_theme + scale_fill_manual(values = c("black","grey")) + xlab("Time series (hr)") + ylab("Distance from AX4 expression") + theme(strip.background = element_rect(linetype = "blank")) + theme(strip.text.x = element_text(size=14))


# predict for timing genes 
dum <- as.vector(rep(1,1254)) # dummy variable to extract random effects later
gam.dist.sub.full <- gam(data=dist.sub, distance ~ condition + s(as.numeric(hr), by=condition) + s(rep, bs="re", by=dum), method="REML")
gam.dist.sub.null <- gam(data=dist.sub, distance ~ condition + s(as.numeric(hr)) + s(rep, bs="re", by=dum), method="REML")
AICcmodavg::AICc(gam.dist.sub.null)
AICcmodavg::AICc(gam.dist.sub.full) 
AICcmodavg::AICc(gam.dist.sub.full) - AICcmodavg::AICc(gam.dist.sub.null) 
anova(gam.dist.sub.null, gam.dist.sub.full, test="Chisq") # significant evidence to suggest condition-specific curves
summary(gam.dist.sub.full)

pdat.sub <- cbind(dist.sub[,c(4:6)],dum)
pdat.sub$hr <- as.numeric(pdat.sub$hr)

step <- 0.1 # step size for finite difference http://stackoverflow.com/questions/14207250/determining-derivatives-from-gam-smooth-object 
pdat.sub2 <- pdat.sub
pdat.sub2$hr <- pdat.sub$hr + step
pdat.sub3 <- pdat.sub
pdat.sub3$hr <- pdat.sub$hr - step

x0 <- predict(gam.dist.sub.full, newdata=pdat.sub, type='lpmatrix')
x1 <- predict(gam.dist.sub.full, newdata=pdat.sub2, type='lpmatrix')
xp <- (x0 - x1)/step
pdat.sub$fd1 <- as.vector(xp %*% coef(gam.dist.sub.full)) # 1st deriv - slope of tangent line. When 0, the curve is potentially turning (local minimum is what I want so 2nd derivative must be positive, or zero)
x2 <- predict(gam.dist.sub.full, newdata=pdat.sub3, type='lpmatrix')
xpp <- (x2+x1-2*x0)/step^2
pdat.sub$fd2 <- as.vector(xpp %*% coef(gam.dist.sub.full)) # 2nd deriv - rate of change of tangent line.  
rm(x0,x1,xp,x2,xpp)

p.ci <- predict(gam.dist.sub.full, newdata=pdat.sub, type="link", se.fit=TRUE)
pdat.sub$upr <- as.vector(gam.dist.sub.full$family$linkinv(p.ci$fit + (2*p.ci$se.fit)))
pdat.sub$lwr <- as.vector(gam.dist.sub.full$family$linkinv(p.ci$fit - (2*p.ci$se.fit)))
pdat.sub$fitted <- fitted(gam.dist.sub.full)
pdat.sub$resid <- residuals(gam.dist.sub.full, type = "response")
pdat.sub$observed.y <- napredict(gam.dist.sub.full$na.action, gam.dist.sub.full$y)

# inspect the gam fit
pdat.sub.fig <- aggregate(pdat.sub$fitted, by=list(pdat.sub$hr, pdat.sub$condition), FUN=mean)
names(pdat.sub.fig) <- c("hr","condition","fitted")

ggplot() + geom_line(data=pdat.sub.fig, aes(x=hr, y=fitted, group=condition, color=condition)) + geom_boxplot(data=dist.sub.fig, aes(x=as.numeric(hr), y=distance, group=condition:as.factor(hr), fill=condition)) + manuscript_theme + scale_fill_manual(values = c("white","#e7003e")) + scale_color_manual(values = c("black","#e7003e")) + xlab("Time series (hr)") + ylab("Distance from AX4 expression") + theme(strip.background = element_rect(linetype = "blank")) + theme(strip.text.x = element_text(size=14)) + scale_x_continuous(breaks=c(1, 6, 12, 19), labels=c("00","05","11","24"))
ggplot() + geom_line(data=pdat.sub.fig, aes(x=hr, y=fitted, group=condition, color=condition)) + manuscript_theme + scale_fill_manual(values = c("white","#e7003e")) + scale_color_manual(values = c("black","#e7003e")) + xlab("Time series (hr)") + ylab("Distance from AX4 expression") + theme(strip.background = element_rect(linetype = "blank")) + theme(strip.text.x = element_text(size=14)) + scale_x_continuous(breaks=c(1, 6, 12, 19), labels=c("00","05","11","24"))

# inspect the residuals
gam.check1 <- ggplot(pdat.sub, aes(sample=resid)) + stat_qq(shape=21) + manuscript_theme + xlab("theoretical quantiles") + ylab("residuals") + scale_fill_manual(values = c("black","grey")) 
gam.check2 <- ggplot(pdat.sub, aes(x=fitted, y=resid)) + geom_point(shape=21) + manuscript_theme + xlab("fitted values") + ylab("residuals")+ scale_fill_manual(values = c("black","grey")) 
gam.check3 <- ggplot(pdat.sub, aes(x=resid)) + geom_histogram() + manuscript_theme + xlab("residuals") + ylab("frequency") #+ scale_fill_manual(values = c("black","grey")) + facet_grid(condition~.)
gam.check4 <- ggplot(pdat.sub, aes(x=fitted, y=observed.y)) + geom_point(shape=21) + manuscript_theme + xlab("fitted values") + ylab("observed values")+ scale_fill_manual(values = c("black","grey"))
grid.arrange(gam.check1, gam.check2, gam.check3, gam.check4)

setEPS()
postscript(file=paste("gam_fit_timing_genes",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=8, height=6)
ggplot() + geom_line(data=pdat.sub.fig, aes(x=hr, y=fitted, group=condition, color=condition)) + geom_boxplot(data=dist.sub.fig, aes(x=as.numeric(hr), y=distance, group=condition:as.factor(hr), fill=condition)) + manuscript_theme + scale_fill_manual(values = c("white","#e7003e")) + scale_color_manual(values=rev(mypalette(2))) + xlab("Time series (hr)") + ylab("Distance from AX4 expression") + theme(strip.background = element_rect(linetype = "blank")) + theme(strip.text.x = element_text(size=14)) + scale_x_continuous(breaks=c(1, 6, 12, 19), labels=c("00","05","11","24"))
dev.off()


#### shape change during development ####
# read compiled imageJ data
chi.all <- read.table("development_photos_parsed_with_splits.txt", h=T)
str(chi.all)
summ_by_strain <- data.frame(cbind(aggregate(Area ~ clone * time, data = chi.all, FUN = mean), 
                                   perim = aggregate(Perim. ~ clone * time, data = chi.all, FUN = mean)[,3], 
                                   circ = aggregate(Circ. ~ clone * time, data = chi.all, FUN = mean)[,3], 
                                   AR = aggregate(AR ~ clone * time, data = chi.all, FUN = mean)[,3], 
                                   round = aggregate(Round ~ clone * time, data = chi.all, FUN = mean)[,3],
                                   solidity = aggregate(Solidity ~ clone * time, data = chi.all, FUN = mean)[,3],
                                   split = c(aggregate(split ~ clone * time, data = chi.all, FUN = mean)[,3], rep(NA,12)),
                                   num = c(aggregate(Area ~ clone * time, data = chi.all, FUN = length)[,3]),
                                   chi = factor(rep(c(rep("TRUE",4),rep("FALSE",8)),3))))
str(summ_by_strain)
GGally::ggpairs(summ_by_strain[,2:8], mapping=ggplot2::aes(colour = time))
GGally::ggpairs(summ_by_strain[,c(2:6,8)], mapping=ggplot2::aes(colour = time)) #Roundness eliminated
GGally::ggpairs(summ_by_strain[,c(2,3,5,6,8)], mapping=ggplot2::aes(colour = time)) #Roundness and perimeter eliminated

Y <- cbind(chi.all$Area, chi.all$Circ., chi.all$AR, chi.all$Solidity)
chi.pca <- prcomp(Y, retx=TRUE, center=TRUE, scale=TRUE)
chi.pca$rotation
summary(chi.pca)
chi.pca.scores <- as.data.frame(chi.pca$x)
chi.pca.all <- cbind(chi.all, chi.pca.scores[,1:2])
chi.pca.all$clone <- factor(chi.pca.all$clone, levels=c("qs6","chi6","qs160","qs4","chi4","qs174","qs18","chi18","qs154","qs17","chi17","qs157"))

# main stats results 1 - chi*time is significant
# could either analyze pc1-3 or correct each variable by its mean: results are same
shape.chi <- nlme::lme(PC1 ~ chi*time, random=~1|clone, data=chi.pca.all)
anova(shape.chi) # chi, time and chi:time significant
summary(shape.chi)

split.chi <- nlme::lme(split ~ chi:time + chi + time + PC1, random=~1|clone, data=subset(chi.pca.all, time!="t3" & split!="NA"))
anova(split.chi) # chi, time, chi:time, and PC1 significant
summary(split.chi)

# posthoc test just to look at area
library(lsmeans)
summ_by_strain$PC1 <- c(aggregate(PC1 ~ clone * time, data=chi.pca.all, FUN=mean)[,3])
area.chi <- nlme::lme(Area ~ chi+time, random=~1|clone, data=summ_by_strain)
contrast_area <- lsmeans(area.chi, ~chi + time)
contrast(contrast_area, alpha=0.05, method="pairwise",adjust=NULL) 

# visualize this difference in slope
split.all <- droplevels(subset(summ_by_strain, time!="t3"))
split1 <- ggplot(split.all, aes(x=time, y=split)) + geom_smooth(se=F, method="lm", size=.5, aes(group=chi, color=chi), position=position_dodge(0.75)) + geom_boxplot(aes(group=chi:time), color="grey") + geom_dotplot(binaxis='y', stackdir='center', dotsize=2, aes(fill=chi), position=position_dodge(0.75)) + scale_fill_manual(name="Treatment", values=rev(mypalette(2)),labels=c(`FALSE`="clonal",`TRUE`="chimeric")) + scale_color_manual(values=rev(mypalette(2))) + manuscript_theme + xlab("Time interval") + ylab("Proportion split") + scale_x_discrete(labels=c("t1" = "loose to tight agg.", "t2" = "tight to tipped agg.")) + guides(color=FALSE, lty=FALSE) + ggtitle("(a)")

split2 <- ggplot(summ_by_strain, aes(x=time, y=num)) + geom_smooth(se=F, method="lm", size=.5, aes(group=chi, color=chi), position=position_dodge(0.75)) + geom_boxplot(aes(group=chi:time), color="grey") + geom_dotplot(binaxis='y', stackdir='center', dotsize=2, aes(fill=chi), position=position_dodge(0.75)) + scale_fill_manual(name="Treatment", values=rev(mypalette(2)),labels=c(`FALSE`="clonal",`TRUE`="chimeric")) + scale_color_manual(values=rev(mypalette(2))) + manuscript_theme + xlab("Developmental stage") + ylab("Number") + scale_x_discrete(labels=c("t1" = "loose", "t2" = "tight", "t3"="tipped agg.")) + guides(color=FALSE, lty=FALSE) + ggtitle("(b)")

split3 <- ggplot(summ_by_strain, aes(x=time, y=Area)) + geom_smooth(se=F, method="lm", size=.5, aes(group=chi, color=chi), position=position_dodge(0.75)) + geom_boxplot(aes(group=chi:time), color="grey") + geom_dotplot(binaxis='y', stackdir='center', dotsize=2, aes(fill=chi), position=position_dodge(0.75)) + scale_fill_manual(name="Treatment", values=rev(mypalette(2)),labels=c(`FALSE`="clonal",`TRUE`="chimeric")) + scale_color_manual(values=rev(mypalette(2))) + manuscript_theme + xlab("Development stage") + ylab("Area") + scale_x_discrete(labels=c("t1" = "loose", "t2" = "tight", "t3"="tipped agg.")) + guides(color=FALSE, lty=FALSE) + ggtitle("(c)")

grid.arrange(split1, split2, split3)

setEPS()
postscript(file=paste("development_agg_splits",format(Sys.time(),"%Y%m%d"),"eps",sep="."), onefile=F, width=8, height=8)
grid.arrange(split1, split2, split3)
dev.off()
