# WGCNA for ZEFI gonad cold vs control

# Set working directory
setwd("~/Dropbox/Rosvall_Postdoc/Heat stress/WGCNA/Heat Stress Gonad/Gonad cold vs control")

# Load the WGCNA package
library(WGCNA);
library("genefilter")
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

# Data formatting  (Part 1)
#Read in ZEFI gonad gene counts - can use raw or normalized counts
gonad = read.csv("ZF_gonad_normalizedCounts_PE_noheat.csv"); # This dataset is of normalized counts from PE reads

# Take a quick look at what is in the data set:
dim(gonad); # 17051    17
names(gonad); # ENSEMBL ID is the first column

# Reformat - get rid of unneccessary columns to only include count data
gonad_normcounts <- data.frame(gonad[,-1], row.names=gonad[,1])
dim(gonad_normcounts) # 17051    15
names(gonad_normcounts)
#We suggest removing features whose counts are consistently low (for example, removing all 
#features that have a count of less than say 10 in more than 90% of the samples) 
#because such low-expressed features tend to reflect noise and correlations based on 
#counts that are mostly zero aren't really meaningful. 
use = rowSums(gonad_normcounts > 10) >=13 # 90% must have at least a normalized count of 10
countMatrixFiltered = gonad_normcounts[ use, ]
dim(countMatrixFiltered) # 11812    30

#We then recommend a variance-stabilizing transformation. Start with normalized counts 
#(or RPKM/FPKM data) and log-transform them using log2(x+1).
#Whether one uses RPKM, FPKM, or simply normalized counts doesn't make a whole lot of difference for WGCNA analysis as long as all samples were processed the same way.
log2countMatrixFiltered = log2(countMatrixFiltered+1)

#Remove the genes with the lowest variance (lower 25 percentile).
log2countMatrixFiltered = as.matrix(log2countMatrixFiltered)
rv = rowVars(log2countMatrixFiltered)
q25 = quantile(rowVars(log2countMatrixFiltered), .25, na.rm=TRUE)
Filtered = log2countMatrixFiltered[rv > q25, ]
dim(Filtered) # Left with 8859 genes out of 17051
write.csv(Filtered,file="ZF_gonad_coldvscontrol_normCounts_logreduced25.csv")

# Read gene annotation file - so that we can match gene names with Ensembl IDs
#annot = read.csv(file = "ZF_gonad_annGeneNames.csv");
annot = read.csv("ZF_gonad_annGeneNames.csv")


########################################### Read in trait data ##########################################
traitData = read.csv("ZF_heat_traits_WGCNA_gonad_noheat.csv");
dim(traitData) # 15 31
names(traitData)

# remove columns that hold information we do not need.
library(tidyverse)
datTrait = subset(traitData, select=-c(Chamber,Mass.Before,Mass.After, Mass.Change,Tarsus.length, Euthanized,Caught,Wing.length, Bill.depth,Bill.width, Bill.length,
                                       Sac.Round,Brain, Gonad, Neutral.Posture, Groom, Wingspread, Out.of.View,Fluff,Move.Down, 
                                       Move.Up, Move.Total,Move,Stand.Tall, Time.Down.s, ratio, Mass.Tarsus.Resid,Thermal.Time,Drink,Eat,Time.aloft.s))
datTrait
dim(datTrait) # 15 11
names(datTrait)

# Form a data frame analogous to expression data that will hold the traits.
datTraits <- datTrait[,-1] # Get rid of first column... did this do anything?
rownames(datTraits) <- datTrait[,1]
datExpr<-read.csv("ZF_gonad_coldvscontrol_normCounts_logreduced25.csv")
dim(datExpr) # 8859   16

# Note that each row corresponds to a gene and column to a sample or auxiliary information.
# We now remove the auxiliary data and transpose the expression data for further analysis.
#datExpr0 = t(datExpr[, -1]) # this gets rid of Ensembl names...
datExpr0 = t(datExpr[,-1])
names(datExpr0) = datExpr$ENSEMBL.ID                    #Use the first column name after '$'
rownames(datExpr0) = names(datExpr)[-1]

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # TRUE

sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

############## Soft connectivity #######

# Choose a set of soft-thresholding powers
# R^2>0.80
# mean connectivity should be high so that the network contains enough information (e.g. for module detection)
#Scale-free topology fit needs to reach values above 0.8 for reasonable powers 
#(less than 15 for unsigned or signed hybrid networks, and less than 30 for signed) 
#and the mean connectivity remains relatively high (in the hundreds or above).
#If this doesn't happen, chances are that the data exhibit a strong driver that 
#makes a subset of the samples globally different from the rest. 
#The difference causes high correlation among large groups of genes which 
#invalidates the assumption of the scale-free topology approximation.
#If you can't meet these criteria,
# default choice of B: for an unsigned/signedhybrid; B=9, for a signed network B=12. for <20 samples (conservative)
#unsigned = positive or negative correlations equal a connection
#signed = negatively correlated genes are considered unconnected
#signed hybrid puts negative correlations at 0 for cosmetic reasons - not that different from signed
#the resulting matrix is always non-negative #s so best to use signed hybrid to avoid confusion
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, networkType="signed hybrid",corOptions = list(maxPOutliers =0.1),
                        powerVector = powers, verbose = 5,corFnc= "bicor")
#### Why doesn't this work with datExpr0, but does work with datExpr?


# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2),mar=c(4.5,4.5,0.5,0.5),oma=c(0,0,0,0));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#We now calculate the adjacencies:
#power of 6 vs 11 for gonad cold vs control
adjacencymatrix = adjacency(datExpr0, power=10, corFnc="bicor", 
                            corOptions = list(maxPOutliers =0.1),type="signed hybrid")

#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix,
#and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacencymatrix, TOMType="signed")
dissTOM = 1-TOM

#We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes.
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2,cutHeight = 0.99, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
#such modules since their genes are highly co-expressed.

### How do you know if your modules need to be merged?

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25 # Corresponds to similarity of.75 to merge
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Yes, looks like we should merge a few modules, e.g. cyan + green yellow, blue + turquoise, black + midnight blue

# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;

MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "ZF_gonad_coldvscontrol_normcounts_modules.RData")

#identify modules that are significantly associated with the measured clinical traits.
#Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external
#traits and look for the most significant associations:
# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = bicor(MEs, datTraits, use = "p",robustY = FALSE,maxPOutliers =0.1);

#moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# color changes for heat map
library(RColorBrewer)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colorRampPalette(brewer.pal(8, "Blues"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

table(moduleColors)
# Black module looks promising! Neg cor with pant and pos with neutral posture
#write.csv(MEs, file = "ZF_gonad_hotvscold_normcounts_nocold_eigengenes_modules.csv")

#focuse on one trait
# Define variable weight containing the weight column of datTrait
Piloerect = as.data.frame(datTraits$Piloerect);
names(Piloerect) = "Piloerect"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(bicor(datExpr0, Piloerect, use = "p",robustY = FALSE,maxPOutliers =0.1));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Piloerect), sep="");
names(GSPvalue) = paste("p.GS.", names(Piloerect), sep="");

#link trait with module
module = "lightgreen"
column = match(module, modNames);
moduleGenes = moduleColors==module;
#sizeGrWindow(7, 7);
#par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Piloerection",
                   main = paste("Module membership vs. gene significance\\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#GS and MM are highly correlated, illustrating that genes highly significantly
#associated with a trait are often also the most important (central) elements of modules associated with the trait

#### What genes are in each module??
#black.ID <- colnames(datExpr)[moduleColors=="black"] ## Couldn't get this to work
#write.csv(black.ID, file = "ZF_gonad_normcounts_black_Ensembl.ID.csv")

# Summary output of network analysis results 
annot = read.csv("ZF_gonad_annGeneNames.csv")
head(annot)
dim(annot) # 17051     2
#annot<-annot[,c(1,18,19)]
names(annot)

# Need to add column names back to dataframe
names(datExpr0) = datExpr$X
probes = datExpr$X
probes2annot = match(probes, annot$ENSEMBL.ID)

# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0

#We now create a data frame holding the following information for all probes: gene ID, gene symbol, Locus Link ID
#(Entrez code), module color, gene significance for weight, and module membership and p-values in all modules. The
#modules will be ordered by their significance for weight, with the most significant ones to the left.


# Create the starting data frame
geneInfo0 = data.frame(geneID = probes,
                       genename = annot$Gene.Name[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for trait of interest
#modOrder = order(-abs(cor(MEs, Pant, use = "p")))
#modOrder = order(-abs(cor(MEs, Neutral.Posture, use = "p")))
modOrder = order(-abs(cor(MEs, Piloerect, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
#geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Pant))
#geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Neutral.Posture))
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Piloerect))

geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "ZF_gonad_genename_piloerect.csv")

# Select modules
modules = c("lightgreen");
#modules = c("cyan");
#modules=c("pink");
# Select module probes
probes = annot$Gene.Name[probes2annot] # Was this the right way to replace ENSEMBL IDs with gene names?

inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$Symbol[match(modProbes, annot$Gene.Name)];
# Select the corresponding Topological Overlap
TOM = TOMsimilarityFromExpr(datExpr0, power = 16);

modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the netLwork into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

chooseTopHubInEachModule(datExpr0,colorh = "black") # can't get this to work... out put is 8168 for both
chooseTopHubInEachModule(datExpr0,colorh = "cyan")


# On Cytoscape 3.8.0, you can import WGCNA network by 
# File -> Import -> Network from File and selecting the module edge file. 
# On the import dialogue box, you will typically select the fromNode column as the Source Node 
# and the toNode column as the Target Node. The weight column should be left as an Edge Attribute. 
# The direction column should be changed to interaction type.
