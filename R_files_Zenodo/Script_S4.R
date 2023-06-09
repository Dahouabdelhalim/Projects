# This is code to replicate the analysis from Alleman et al. "Tandem-Running and Scouting Behavior is Characterized by Up-Regulation of Learning and Memory formation genes within the Ant Brain"
# Weighted Gene Co-expression Analysis (WGCNA) of T. longispinosus
# Transcriptome constructed using Trinity and filtered for those contigs that contain an open reading frame
# Expected read counts created using RSEM and filtered for contigs with at least 10 count in at least 4 samples and maximum Cook's distance below 38
# Script modified from Tutorials provided on https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/ by M. Stoldt

####################################################################################################################
## Loading required packages
####################################################################################################################

library(DESeq2)
library(ggpubr)
library(WGCNA)
options(stringsAsFactors = FALSE)

####################################################################################################################
## Establishment of directories
####################################################################################################################

# Specify working directory for project.
setwd("/My/Working/Directory/")

####################################################################################################################
## Preparation of data
####################################################################################################################

# Read in filtered countsmatrix
load("matrix_cook_Tlongi.RData")

# Transform matrix using variance stabilizing transformation provided in DESeq2
vsd <- varianceStabilizingTransformation(matrix_cook, blind=FALSE)

####################################################################################################################
## Cleaning the data
####################################################################################################################

# Transpose transformed matrix
datExpr0 = as.data.frame(t(vsd))

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# Create dataframe containing the behavioral group of each sample in the countsmatrix
caste=matrix(c("Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout"), ncol=1)
rownames(caste)=rownames(datExpr0)
colnames(caste)=c("Behavior")

# Cluster with information on behavioral groups
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = labels2colors(caste, colorSeq=c( "paleturquoise3","palevioletred2", "steelblue")) 
plotDendroAndColors(sampleTree2,traitColors,groupLabels = colnames(caste),main ="Sample dendrogram and trait heatmap", autoColorHeight=FALSE,colorHeight = 0.15)

####################################################################################################################
## Blockwise network contruction
####################################################################################################################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# threshold
abline(h=0.87,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# Pick the lowest value that touches the threshold as power
softPower = 4

# Run blockwise network construction
bwnet = blockwiseModules(datExpr0, maxBlockSize = 20000,
                         power = 4, TOMType = "unsigned", minModuleSize = 50,
                         reassignThreshold = 0, mergeCutHeight = 0.3,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Tlongi_unfilt_TOM-blockwise",
                         verbose = 3)

####################################################################################################################
# Relating modules to behavioral group
####################################################################################################################

moduleColors<-bwnet$colors
# Calculate module eigengenes 
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

group<-as.factor(c("Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout","Follower", "Leader", "Scout"))

# Create "empty" dataframe
results_all_tl<-data.frame(module=c(rep(3,length(MEs))), pvalue=c(rep(3,length(MEs))), p_fo_le=c(rep(3,length(MEs))), p_fo_sc=c(rep(3,length(MEs))), p_le_sc=c(rep(3,length(MEs))))

# Run over all modules and test whether their module eigengenes can be explained by behavioral group using Kruskal-Wallis Test
# followed by pairwise Wilcoxon test adjusting for multiple testing using Benjamini-Hochberg procedure
for (i in 1:length(MEs)){
  results_all_tam$module[i]<-as.character(names(MEs[i]))
  x<-cbind(group,MEs[i])
  colnames(x)<-c("group", "moduleeigengene")
  y<-kruskal.test(moduleeigengene ~group, data=x)
  results_all_tl$pvalue[i]<-y$p.value
  
  Wilc = pairwise.wilcox.test(x$moduleeigengene,x$group,p.adjust.method="BH")
  results_all_tl$p_fo_le[i]<-Wilc$p.value[1]
  results_all_tl$p_fo_sc[i]<-Wilc$p.value[2]
  results_all_tl$p_le_sc[i]<-Wilc$p.value[4]
  
  # Create according boxplot of module eigengene
  z<-ggboxplot(x, x = "group", y = "moduleeigengene", 
               color = "group", palette = c( "paleturquoise3","palevioletred2", "steelblue"),
               order = c("Follower", "Leader", "Scout"),
               ylab = paste("Module eigengene ",names(MEs[i]), sep = ""), xlab = "Caste", legend.title="Caste")+stat_compare_means()
  
  ggexport(z,filename = paste("Boxplot_wgcna_",names(MEs[i]),".png", sep=""))
  
  # Write the contigs belonging to each module in a separate file
  name_ME<-strsplit(names(MEs[i]), "E")[[1]][2]
  w<-names(datExpr0)[moduleColors==name_ME]
  write.table(w, paste("ModuleContigList_",names(MEs[i]),"_Tlongi.txt", sep=""), sep="\\t",  row.names = FALSE, quote = FALSE)
  
}

# Write the test results in a separate file
write.table(results_all_tl, file="Results_KruskalWallis_WGCNA_Tlongi.txt", sep="\\t", row.names = FALSE, quote = FALSE)






