library("DESeq2")
library("WGCNA")
##normalise data and set up input data for WGCNA
counts <- read.table("../counts_all.csv",sep=",",row.names=1,header=T)
#convert to integers:
for (i in 1:27){
 counts[,i] <- as.integer(counts[,i])
}
design <- read.table("../coldata.csv",sep=",",row.names=1,header=T)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design= ~nutrition)
dds <- DESeq(dds)

###for DE######
XXX <- as.data.frame(subset(results(dds),padj < 0.05))
###############

datExpr <- t(counts(dds, normalized=TRUE))
##remove genes with expression less than 10 in > 13 samples (50%)
datExpr <- datExpr[ ,colSums(datExpr >= 10) >= 11 ]

#dim(datExpr)
#[1]    27 10320
#################
worker_bin <- c(rep(1,3),rep(0,23),1)   
T0_bin <- c(rep(0,3),rep(1,4),rep(0,20))
T1_bin <- c(rep(0,7),rep(1,4),rep(0,16)) 
T2_bin <- c(rep(0,11),rep(1,2),rep(0,14)) 
T4_bin <- c(rep(0,13),rep(1,4),rep(0,10)) 
king_bin <- c(rep(0,17),rep(1,3),rep(0,7))  
T0ov_bin <- c(rep(0,20),rep(1,2),rep(0,5))  
T4ov_bin <- c(rep(0,22),rep(1,4),0)  

design <- cbind(worker_bin, T0_bin,T1_bin,T2_bin,T4_bin,king_bin,T0ov_bin, T4ov_bin)
row.names(design) <- row.names(datExpr)

#####################


##start WGCNA
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

##check for outliers
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

mad.datExpr = apply(datExpr,2,mad)
datExpr = datExpr[ , colnames(datExpr) %in% names(mad.datExpr[mad.datExpr > 0])]
#dim(datExpr)
#[1]    27 10318
powers = 1:20
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05) )
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##choose threshold of 14
softPower = 14
adjacency = adjacency(datExpr, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05));
TOM = TOMsimilarity(adjacency);
colnames(TOM) = colnames(adjacency)
rownames(TOM) = rownames(adjacency)
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);

minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
#table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

##merging of modules that are similar in expression profiles

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
MEDissThres = 0.5
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;



###relate modules to traits with aov
for (j in 1:length(colnames(design))){
trait = j
trait.aov = data.frame() #repeat 

for (i in 1:length(names(mergedMEs))  ){           
    trait.aov = rbind(trait.aov, unlist(summary(aov(mergedMEs[ , i] ~ design[ , trait]))))
}

names(trait.aov) = c("DF1", "DF2", "SS1", "SS2", "MS1", "MS2", "F1", "F2", "PR1", "PR2")
row.names(trait.aov) = names(mergedMEs)

assign(paste(colnames(design)[trait],".aov",sep=""),trait.aov)
}

trait_pvalues <- cbind(nutrition.aov$PR1, age.aov$PR1,caste.aov$PR1,gender.aov$PR1,tissue.aov$PR1,colony.aov$PR1)
colnames(trait_pvalues) <- colnames(design)
row.names(trait_pvalues) <- names(mergedMEs)

#get fdr
trait_fdr <- cbind(
p.adjust(trait_pvalues[,1],method = "fdr"),
p.adjust(trait_pvalues[,2],method = "fdr"),
p.adjust(trait_pvalues[,3],method = "fdr"),
p.adjust(trait_pvalues[,4],method = "fdr"),
p.adjust(trait_pvalues[,5],method = "fdr"),
p.adjust(trait_pvalues[,6],method = "fdr")
)
colnames(trait_fdr) <- colnames(trait_pvalues) 

trait_fdr_bin <- trait_fdr
for(i in 1:length(row.names(trait_fdr) )){
    for (j in 1:length(colnames(trait_fdr))){
        if (trait_fdr[i,j] >= 0.05){
            trait_fdr_bin[i,j] <- NA}
        else{
        trait_fdr_bin[i,j] <- trait_fdr[i,j]}
    }
}
##create heatmap
library(RColorBrewer)

pdf("signed_module_traits.pdf")
red_palette <- colorRampPalette(c("red","yellow"))
heatmap.2(trait_fdr_bin[,c(1,2,4)],Rowv = FALSE, Colv=FALSE, dendrogram="none",trace="none", col=red_palette(50))
dev.off()

####post-hoc test
#run tukeyhsd on each aov to get information on which levels of each trait are significant, e.g.:

TukeyHSD(aov(mergedMEs[,"MEyellow"]~ design[,"tissue"]))

###relate modules to traits
moduleColors = mergedColors

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs <- MEs[,c(1,4,3,2,8,5,6,7)]
moduleTraitCor = cor(MEs, design_fat[,c(7,14,15,9:13)], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf("mnat_all_module_traits_new.pdf")
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3), cex=0.5);
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = c("W","T0a","T0b","T1","T2","T3","T4","K"),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships Mnat"))

dev.off()

###fat only
pdf("mnat_all_module_traits_fat_only.pdf")
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor[ 1:8 , 1:6 ], 2), "\\n(",
signif(moduleTraitPvalue[ 1:8 , 1:6 ], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[ 1:8 , 1:6 ])
par(mar = c(6, 8.5, 3, 3), cex=0.5);
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor[ 1:8 , 1:6 ],
xLabels = c("W","T0","T1","T2","T4","K"),
yLabels = names(MEs[1:8]),
ySymbols = names(MEs[1:8]),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1)
)

dev.off()


###summarise gene info in a table
geneTraitSignificance = as.data.frame(cor(datExpr, worker_bin, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

modNames = substring(names(mergedMEs), 3)
geneModuleMembership = as.data.frame(signedKME(datExpr,  mergedMEs, corFnc = "bicor", corOptions = list(maxPOutliers = 0.05)));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneInfo0 = data.frame(genes = colnames(datExpr),
                    moduleColor = mergedColors,
                    geneTraitSignificance,
                    GSPvalue
                    )

modOrder = order(-abs(cor(mergedMEs,worker_bin, use = "p")))


for (mod in 1:ncol(geneModuleMembership))
{
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
MMPvalue[, modOrder[mod]]);
names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$V1));
geneInfo = geneInfo0[geneOrder, ]

##calculate connectivity
mnat.connectivity <- softConnectivity( datExpr, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(maxPOutliers = 0.05)  )
geneInfo <- cbind(geneInfo, connectivity = mnat.connectivity)


##intramodular connectivity
kIM <- intramodularConnectivity(adjacency, mergedColors, scaleByMax = TRUE) 
geneInfo <- cbind(geneInfo, IntraConnectivity = kIM[,2])

write.table(geneInfo, "geneInfo_mnat_all.tab", sep="\\t",quote=F, row.names=F)



####export to cytoscape

#annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("blue");
# Select module probes
probes = colnames(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = TFs$TF[match(modProbes, TFs$Gene)]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
weighted = TRUE,
#threshold = 0.1243,
nodeNames = modProbes,
#altNodeNames = modGenes,
nodeAttr = moduleColors[inModule]);

##for reduced data set - only genes with at least one TOM of 0.1, all modules (except blue)
TOM_red <- TOM[rowSums(TOM >= 0.2) > 1, colSums(TOM >= 0.2) > 1]
modules = c("skyblue","darkred","red","darkturquoise","cyan","black","orange","pink")
probes = as.character(geneInfo[ geneInfo$genes %in% colnames(TOM_red) & geneInfo$moduleColor %in% modules,1])
moduleColor_red <- as.character(geneInfo[ geneInfo$genes %in% probes,"moduleColor"])
#inModule = is.finite(match(moduleColor_red, modules));
#modProbes = probes[inModule];
modTOM = TOM_red[probes, probes];
#dimnames(modTOM) = list(modProbes, modProbes)
node_attr = annotations[ annotations$gene_id %in% probes, c(7,8,13:17,22:27)]
node_attr$modules = moduleColor_red
cyt = exportNetworkToCytoscape(modTOM,edgeFile = "Cytoscape_edges_test.txt", nodeFile = "Cytoscape_nodes_test.txt",weighted = TRUE, nodeNames = probes, nodeAttr = node_attr, threshold=0.15)
                                                                                                                                                                                                                                                                    


###get genes most connected to TFs
gene_list <- TFs[ TFs$TF == "C2H2" & TFs$Gene %in% geneInfo[ geneInfo$moduleColor == "blue",1],1]
tf_connected_genes <- c()                                        
for(i in 1:134){                                                 
 gene <- gene_list[i]                                             
 TOM_blue_sort <- TOM_blue[ , colnames(TOM_blue) == gene]         
 TOM_blue_sort <- TOM_blue_sort[order(- TOM_blue_sort)]           
 tf_connected_genes <- c(tf_connected_genes, names(TOM_blue_sort)[2:11])
}
##do go term enrichment on these!


