###### Euglossa dilemma experimental disruption analysis #####################

#Variables for each behavioral group are represented by a single letter throughout all analyses
# Group ID: D = dominant, G = Natural Guard, R = Isolated Subordinate, S = Subordinate
# depending on the data set/analysis these groups may be represented by upper or lowercase letters, which will be specified before each analysis


### Brood size (number of brood cells) comparison between isolated subordinates and natural guards ######

#required packaged for analysis
library(car)
library(Steel.Dwass.test)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(vegan)
library(ecodist)
library(gplots)
library(MASS)
library(RVAideMemoire)
library(corrplot)
library(randomForest)

BroodData = read.csv("BroodSize.csv")
Behavior = BroodData$Behavior
Cells = BroodData$Cells


#Brood Cells vs Behavior
fit <- aov(Cells ~ Behavior , data=BroodData)
summary(fit)
leveneTest(Cells~Behavior, data =BroodData)
shapiro.test(Cells)

#data pass homogeneity of variances and normality so no need for kruskal test
kruskal.test(Cells~Behavior, data = BroodData)


#Boxplots

gp = ggplot(BroodData, aes(x = Behavior, y = Cells))+
        geom_boxplot(width=0.25)+
        #geom_point(aes())+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        geom_point(aes(fill = Behavior,shape = Behavior), size = 2, position=position_jitter(0.1)) +
        scale_shape_manual(values = c(22,23,23,24,25))+
        scale_fill_manual(breaks = c("G", "R"),
                          values=c("purple", "palegreen4"))+
        ylim(0,15)

gp

### Ovary size comparison among Reproductives, Natural guards, and Isolated subordinates ###

#Note: In this data set dominants and subordinates are lumped together and labeled "A" instead of the "D" and "S" throughout the other files

OvData = read.csv("OvaryData.csv")

Behavior = OvData$Reproductive_Behavior
Ovary_Index = OvData$TotalIndex

#Ovary Index vs Behavior
fit <- aov(Ovary_Index ~ Behavior , data=OvData)
summary(fit)
TukeyHSD(fit)
leveneTest(Ovary_Index~Behavior, data =OvData)
shapiro.test(Ovary_Index)

# homogeneity of variances and normality confirmed, kruskal test not needed
kruskal.test(Ovary_Index~Behavior, data = OvData)

##Boxplots ###

#behavior with dominants and subordinates separate for plotting
Behavior_fine = OvData$Behavior_finescale

gp = ggplot(OvData, aes(x = Behavior, y = Ovary_Index))+
        geom_boxplot(width=0.25)+
        geom_point(aes(fill = Behavior_fine,shape = Behavior_fine), size = 3, position=position_jitter(0.1)) +
        scale_shape_manual(values = c(21,22,23,24,25))+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gp + scale_fill_manual(breaks = c("D", "G", "R", "S"),
                       values=c("red", "purple", "palegreen4", "blue"))

## Body size test ##

BodySize = OvData$Intertegular


fit <- aov(BodySize ~ Behavior , data=OvData)
summary(fit)
TukeyHSD(fit)
leveneTest(BodySize~Behavior, data =OvData)
shapiro.test(BodySize)

#fails normality so kruskal instead
kruskal.test(BodySize~Behavior, data = OvData)
Steel.Dwass(BodySize, OvData$Reproductive_Behavior_numbered)

##body size boxplot
gp = ggplot(OvData, aes(x = Behavior, y = BodySize))+
        geom_boxplot(width=0.25)+
        geom_point(aes(fill = Behavior_fine,shape = Behavior_fine), size = 2, position=position_jitter(0.1)) +
        scale_shape_manual(values = c(21,22,23,24,25))+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gp + scale_fill_manual(breaks = c("D", "G", "R", "S"),
                       values=c("red", "purple", "palegreen4", "blue"))


##### Ovary size comparison with Saleh and Ramirez, 2019 ##########

## "OG" in the data file refers to the samples from Saleh and Ramirez, 2019 and "G" refers to samples from this study

OvData = read.csv("Ovaries_Cross_Study_Comparison.csv")
Behavior = OvData$Behavior
Behavior_fine = OvData$Behavior_fine
OvarySize = OvData$Ovary_Index

#Peaks vs Behavior
fit <- aov(OvarySize ~ Behavior , data=OvData)
summary(fit)

## Boxplot ###

gp = ggplot(OvData, aes(x = Behavior, y = OvarySize))+
        geom_boxplot(width=0.25)+
        geom_point(aes(fill = Behavior_fine,shape = Behavior_fine), size = 2, position=position_jitter(0.1)) +
        scale_shape_manual(values = c(21,22,23,24,25))+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gp

########################## CHC analysis ##########################
### NMDS first round is done with rare polymorphism samples present. Then the samples will be removed for Random Forest analysis as in the main text###########

CHC = read.csv("CHCs_all_samples.csv")
dim(CHC)

##find relative abundances (i.e. percents)
rsums <- rowSums (CHC[,c(4:20)])
CHC.norm <- CHC[,c(4:20)]/rsums
rowSums(CHC.norm)
head(CHC.norm)

##add back in descriptive columns
CHC <- cbind(CHC[,1:3], CHC.norm)
head(CHC)
#write.matrix(CHC, file = 'Relative_abundance.txt')
#########################
##heatmap that shows clustering of the four samples with the rarer polymorphism
heatmap.2(as.matrix(CHC[,4:20]), Colv=NA, col=terrain.colors(250), labRow=paste (CHC$Sample,sep="_"), main="Female CHCs", cexRow=2, cexCol=.5, trace="none")

##########################
##nMDS plots

##create triangular distance matrix
bc.dist <- bcdist(CHC[,4:20]) 
bc.dist

### nMDS
CHC
CHC.nmds <- nmds (bc.dist, mindim=2, maxdim=2, nits=10)
CHC.nmin <- nmds.min (CHC.nmds)

CHC.nmin <- cbind(CHC.nmin, CHC$Sample, CHC$Behavior, CHC$Behavior_numerical)
CHC.nmin
head(CHC.nmin)
colnames (CHC.nmin)[3]<- "Sample"
colnames (CHC.nmin)[4]<- "Behavior"
colnames (CHC.nmin)[5]<- "Behavior_numerical"
head(CHC.nmin)

##plot nMDS
#this is an exploratory plot, so symbols and numbers don't match the main text. The ggplot version from the main text is below
colors <- c("blue","red","purple","green", "orange", "olivedrab","cyan","firebrick","gray","magenta","lightgreen", "yellow", "tan", "thistle", "deepskyblue4","limegreen","wheat", "sienna","orchid4","black")

plot(CHC.nmin$X1, CHC.nmin$X2, cex=1, main="Female CHCs",col=colors[factor(CHC.nmin$Behavior)],pch=c(8,17,16,15,18,14,19,5,6,7)[factor(CHC.nmin$Behavior)], xlab="", ylab="")

text(CHC.nmin$X1, CHC.nmin$X2, pos=1, labels=CHC.nmin$Sample, cex=0.2)

##Nicer ggplot version

CHC.nmin$Behavior_numerical = as.factor(CHC.nmin$Behavior_numerical)
gp =ggplot(CHC.nmin, aes(x=CHC.nmin$X1, y=CHC.nmin$X2)) +
        geom_point(aes(fill = CHC.nmin$Behavior_numerical ,shape = CHC.nmin$Behavior_numerical), size = 4) +
        scale_shape_manual(values = c(21,24,22,23))+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
gp

gp + scale_fill_manual(breaks = c("1", "2", "3", "4","5"),
                       values=c("red", "blue", "purple", "palegreen4"))


##################### Random Forest analysis with polymorphism samples removed as in main text ############
#Code based on Bruckner and Heethoff, 2017 (see references in main text)

##Important note: in this dataset, groups have numerical labels. "1" = dominant, "2" = subordinate, "3" = natural guard, and "4" = isolated subordinate

my.data = read.csv("Relative_abundance_randomForest.csv")

my.data$group = as.factor(my.data$group)


#note: rule for mtry: mtry=sqrt(n_compounds)
my.data.rf<-randomForest(group~.,my.data, mtry=17,importance=TRUE, do.trace=1000, ntree=10000, proximity=TRUE, keep.forest=TRUE) 
plot(my.data.rf) # trace the groups errors related to ntree


##confusion matrix of the randomForests
my.data.rf 


#plot a MDS plot based on the classification of randomForests
#Important note: default colors of the points were edited and do not match the figure in the main text (shapes are the same). Be aware when comparing this plot to the main text. 
MDS = MDSplot(my.data.rf, my.data$group, cex=1.6, pch=c(16,17,15,18)[my.data$group]) 



# Variable importance plot
varImpPlot(my.data.rf, sort=TRUE, n.var=min(30, nrow(my.data.rf$importance)),
           type=NULL, class=NULL, scale=TRUE,
           main=deparse(substitute(my.data.rf)))

# gives the variable importance 
varimp<-importance(my.data.rf,scale=TRUE)[,6]
sort(varimp)

################ EdgeR and WGCNA transcriptome analysis for the brains #############
#required packages for DE and WGCNA analysis
library(edgeR)
library(MASS)
library(gplots)
library(RColorBrewer)
library(WGCNA)
library(igraph)
library(limma)
library(remotes)

## We will start with the brain samples. 

##Note: This run through has "o_SRNS33.txt" removed due to its outlier status. The analysis can be rerun by inserting this sample back into the blank spot in the column below if desired.
targets = read.delim("targets_brains.txt", stringsAsFactors = FALSE)
targets
d = readDGE(targets)
colnames(d) = c("a_SRNS1.txt",
                "b_SRNS3.txt",
                "c_SRNS5.txt",
                "d_SRNS7.txt",
                "e_SRNS9.txt",
                "f_SRNS11.txt",
                "g_SRNS13.txt",
                "h_SRNS15.txt",
                "i_SRNS17.txt",
                "j_SRNS21.txt",
                "k_SRNS23.txt",
                "l_SRNS25.txt",
                "m_SRNS29.txt",
                "n_SRNS31.txt",
               
                "p_SRNS35.txt",
                "q_SRNS37.txt",
                "r_SRNS39.txt",
                "s_SRNS41.txt",
                "t_SRNS43.txt",
                "u_SRNS45.txt",
                "v_SRNS47.txt",
                "w_SRNS49.txt",
                "x_SRNS53.txt",
                "y_SRNS55.txt",
                "z_SRNS57.txt"
                
                
                
)
d$samples
dim(d)
#filter out genes of low expression
keep <- rowSums(cpm(d)>1) >= 5
d <- d[keep, , keep.lib.sizes=FALSE]
d = calcNormFactors(d)

#Show sample info with normalization factors
d$samples

# MDS for batch check, outlier check, etc.
mds = plotMDS(d,method= "bcv")
plot(mds, col=c('red', 'blue', 'darkgreen', 'purple','orange')[d$samples$group], pch = 20, cex = 2)

# Run EdgeR-robust model. For simplicity variables for each behavioral group are represented by a single letter
# Group ID: d = dominant, g = Natural Guard, r = Isolated Subordinate, S = Subordinate

d = estimateGLMRobustDisp(y = d,design = model.matrix(~ d$sample$group))

#The contrast coefficients refer to the groups in alphabetical order. Change to run desired contrast.
# Reminder: variables for each behavioral group are represented by a single letter throughout all analyses
# Group ID: d = dominant, d = Natural Guard, r = Isolated Subordinate, s = Subordinate
fit = glmFit(d, design = model.matrix(~ 0+d$sample$group))
et = glmLRT(fit, contrast=c(0,0,1,-1))
                          # d g r s 

#Shows how many genes are up and down regulated in the chosen contrast.

summary(de <- decideTestsDGE(et, p=0.05))

#Shows the top 100 genes by Pvalue
topTags(et, n=100)


#Plot to visually investigate DE results for chosen contrast.
detags <- rownames(d)[as.logical(de)]
detags
plotSmear(et, de.tags=detags, cex = 0.7)
abline(h = c(-2, 2), col = "blue")

#Print and save a table showing the info (including FDR Pvalues) for all genes in the chosen contrast
toptable = topTags(et, n=Inf)
toptable
#write.table(toptable, file = 'Isolated_subordinates_vs_Subordinates_Brain.txt')


###### Make a heatmap of significant genes for chosen contrast

#create log counts per million using standard settings
logCPM <- cpm(d, prior.count=2, log=TRUE)

#apply group names to columns
colnames(logCPM) = paste(d$samples$group)

#order by 
o = order(et$table$PValue)

#only use the first 94 genes, which are the significant DEGS from the contrast of interest. Change value to reflect # of DEGS in chosen contrast
logCPM = logCPM[o[1:94],]

### scale values for heatmap
logCPM = t(scale(t(logCPM)))

# choose methods for distance and clustering 
dist.method <- function(logCPM) dist(logCPM, method="euclidean")
hclust.method <- function(logCPM) hclust(logCPM, method="ward.D2")

#change colors and visual settings for actual heatmap
col.cell <- c("red", "purple","palegreen4","blue")[d$samples$group]
col.pan <- colorpanel(1000, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", hclustfun = hclust.method, 
          key = 'TRUE', ColSideColors = col.cell )


### Using the approach above I found all DEGs among the 6 possible pairwise comparisons. The resulting heatmap shown in the main text is from those 132 DEGs.

logCPM_2 = read.csv("All_132_Brain_Behavior_DEGs.csv")
logCPM_3 = logCPM_2[2:26]
colnames(logCPM_3) = paste(d$samples$group)
logCPM_3

logCPM = as.matrix(logCPM_3)

dist.method <- function(logCPM) dist(logCPM, method="euclidean")
hclust.method <- function(logCPM) hclust(logCPM, method="ward.D2")
hclust.method


col.cell <- c("red", "purple","palegreen4","blue")[d$samples$group]
col.pan <- colorpanel(1000, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", hclustfun = hclust.method, 
          key = 'TRUE', ColSideColors = col.cell )


##################### WGCNA for Brains on Same Data Set ##################################
### WGCNA ####

# Load the WGCNA package

#this code is based on the WGCNA tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# This data is the filtered/normalized CPM data such as used in the heatmap above, except with the whole gene set and not just a selection as in the heatmap.
BrainData = read.csv("Brain_WGCNA.csv")
dim(BrainData)
names(BrainData)

datExprBRAIN = as.data.frame(t(BrainData[,-c(1:1)]));
names(datExprBRAIN) = BrainData$X;
rownames(datExprBRAIN) = names(BrainData)[-c(1:1)];
gsg = goodSamplesGenes(datExprBRAIN, verbose = 3);
gsg$allOK


###### brain data network construction #####

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprBRAIN, powerVector = powers, verbose = 5)
# Plot the results:
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
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

### We chose a threshold power 5 for proceeding

###### now brain data module detection #########

net = blockwiseModules(datExprBRAIN, power = 5, maxBlockSize = 11100,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "BrainTOM",
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "BRAIN-02-networkConstruction-auto.RData")

#save a table with the eigengene values for each module
write.table(MEs, file = "eigengenesBRAIN.txt")

#save a table containing the gene IDs for each module
head(net$colors)
write.table(net$colors, file = "geneIDSmodulesBRAIN.txt")


net$colors
#### kmes ###
datKME = signedKME(datExprBRAIN, MEs, outputColumnName = "kME")
head(datKME)

#save a table of kme values, indicating hub genes, for the different modules
write.csv(datKME, file = "KME_values_Brain_modules.csv")





########## Same analyses but on ovary transcriptomes ###########


targets = read.delim("targets_ovaries.txt", stringsAsFactors = FALSE)
targets
d = readDGE(targets)
colnames(d) = c("a_SRNS2.txt",
                "b_SRNS4.txt",
                "c_SRNS6.txt",
                "d_SRNS8.txt",
                "e_SRNS10.txt",
                "f_SRNS12.txt",
                "g_SRNS14.txt",
                "h_SRNS16.txt",
                "i_SRNS18.txt",
                "j_SRNS20.txt",
                "k_SRNS22.txt",
                "l_SRNS24.txt",
                "m_SRNS26.txt",
                "n_SRNS28.txt",
                "o_SRNS30.txt",
                "p_SRNS32.txt",
                "q_SRNS34.txt",
                "r_SRNS36.txt",
                "s_SRNS38.txt",
                "t_SRNS40.txt",
                "u_SRNS42.txt",
                "v_SRNS44.txt",
                "w_SRNS46.txt",
                "x_SRNS48.txt",
                "y_SRNS50.txt",
                "z_SRNS52.txt",
                "za_SRNS54.txt",
                "zb_SRNS56.txt",
                "zc_SRNS58.txt"
                
                
)
d$samples
dim(d)
#filter out genes of low expression
keep <- rowSums(cpm(d)>1) >= 5
d <- d[keep, , keep.lib.sizes=FALSE]
d = calcNormFactors(d)

#Show sample info with normalization factors
d$samples

# MDS for batch check, outlier check, etc.
mds = plotMDS(d,method= "bcv")
plot(mds, col=c('red', 'blue', 'darkgreen', 'purple','orange')[d$samples$group], pch = 20, cex = 2)

# Run EdgeR-robust model. For simplicity variables for each behavioral group are represented by a single letter
# Group ID: d = dominant, g = Natural Guard, r = Isolated Subordinate, S = Subordinate

d = estimateGLMRobustDisp(y = d,design = model.matrix(~ d$sample$group))

#Fit model. The contrast coefficients refer to the groups in alphabetical order. Change to run desired contrast.
# Reminder: variables for each behavioral group are represented by a single letter throughout all analyses
# Group ID: d = dominant, g = Natural Guard, r = Isolated Subordinate, s = Subordinate
fit = glmFit(d, design = model.matrix(~ 0+d$sample$group))
et = glmLRT(fit, contrast=c(0,0,1,-1))
                          # d g r s 

#Shows how many genes are up and down regulated in the chosen contrast.

summary(de <- decideTestsDGE(et, p=0.05))

#Shows the top 100 genes by Pvalue
topTags(et, n=100)


#Plot to visually investigate DE results for chosen contrast.
detags <- rownames(d)[as.logical(de)]
detags
plotSmear(et, de.tags=detags, cex = 0.7)
abline(h = c(-2, 2), col = "blue")

#Print and save a table showing the info (including FDR Pvalues) for all genes in the chosen contrast
toptable = topTags(et, n=Inf)
toptable
#write.table(toptable, file = 'Isolated_subordinates_vs_Subordinates_Ovaries.txt')


###### Make a heatmap of significant genes for chosen contrast

#create log counts per million using standard settings
logCPM <- cpm(d, prior.count=2, log=TRUE)

#apply group names to columns
colnames(logCPM) = paste(d$samples$group)

#order by 
o = order(et$table$PValue)

#only use the first 73 genes, which are the significant DEGS from the contrast of interest. Change value to reflect # of DEGS in chosen contrast
logCPM = logCPM[o[1:73],]

### scale values for heatmap
logCPM = t(scale(t(logCPM)))

# choose methods for distance and clustering 
dist.method <- function(logCPM) dist(logCPM, method="euclidean")
hclust.method <- function(logCPM) hclust(logCPM, method="ward.D2")

#change colors and visual settings for actual heatmap
col.cell <- c("red", "purple","palegreen4","blue")[d$samples$group]
col.pan <- colorpanel(1000, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", hclustfun = hclust.method, 
          key = 'TRUE', ColSideColors = col.cell )

### Using the approach above I found all DEGs among the 6 possible pairwise comparisons. The resulting heatmap shown in the main text is from those 412 DEGs.


logCPM_2 = read.csv("All_412_Ovary_Behavior_DEGs.csv")
logCPM_3 = logCPM_2[2:30]
colnames(logCPM_3) = paste(d$samples$group)
logCPM_3

logCPM = as.matrix(logCPM_3)

dist.method <- function(logCPM) dist(logCPM, method="euclidean")
hclust.method <- function(logCPM) hclust(logCPM, method="ward.D2")
hclust.method


col.cell <- c("red", "purple","palegreen4","blue")[d$samples$group]
col.pan <- colorpanel(1000, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",
          trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", hclustfun = hclust.method, 
          key = 'TRUE', ColSideColors = col.cell )

##################### WGCNA for Ovaries on Same Data Set ##################################
### WGCNA ####

# Load the WGCNA package

#this code is based on the WGCNA tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# This data is the filtered/normalized CPM data such as used in the heatmap above, except with the whole gene set and not just a selection as in the heatmap.
OvaryData = read.csv("Ovaries_WGCNA.csv")
dim(OvaryData)
names(OvaryData)

datExprOVA = as.data.frame(t(OvaryData[,-c(1:1)]));
names(datExprOVA) = OvaryData$X;
rownames(datExprOVA) = names(OvaryData)[-c(1:1)];
gsg = goodSamplesGenes(datExprOVA, verbose = 3);
gsg$allOK

###### Ovary data network construction #####

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprOVA, powerVector = powers, verbose = 5)
# Plot the results:
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
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

### We chose a threshold power 4 for proceeding

###### now ovary data module detection #########

net = blockwiseModules(datExprOVA, power = 4, maxBlockSize = 10132,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "OvaryTOM",
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "OVARY-02-networkConstruction-auto.RData")

#save a table with the eigengene values for each module
write.table(MEs, file = "eigengenesOVARIES.txt")

#save a table containing the gene IDs for each module
head(net$colors)
write.table(net$colors, file = "geneIDSmodulesOVARIES.txt")


net$colors
#### kmes ###
datKME = signedKME(datExprOVA, MEs, outputColumnName = "kME")
head(datKME)

#save a table of kme values, indicating hub genes, for the different modules
write.csv(datKME, file = "KME_values_Ovaries_modules.csv")

############ Now use the results of WGCNA analysis to look for correlations with ovary size index ############

####### ovary module eigengenes#######
eigen = read.csv("Ovary_module_eigengenes.csv")


module = eigen$ME6
Behavior = eigen$group
OvarySize = eigen$TotalOvaryIndex


#plot ovary size vs module eigengenes

gp =ggplot(eigen, aes(x=OvarySize, y=eigen$ME10)) +
        geom_point(aes(fill = Behavior ,shape = Behavior), size = 4) +
        scale_shape_manual(values = c(21,22,23,24))+
        geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
        stat_cor(method = "spearman", color = 1)+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("r", "g", "s", "d"),
                       values=c("palegreen4", "purple", "blue", "red"))


#perform correlation test without plot to generate p values for each module
cor.test(x=OvarySize, y=eigen$ME1, method = "spearman")

#collected all Pvalues for FDR adjusment. Will adjust for all 26 comparisons to ovary size after brain analysis as well
p =  c(0.3136,
       0.01207,
       0.02846,
       0.3168,
       0.00513,
       0.00000000008053,
       0.003737,
       0.9121,
       0.000004734,
       0.00000000006353,
       0.01146,
       0.009439,
       0.1949)

adjusted = p.adjust(p, method = "fdr", length(p))
adjusted


##### Brain module eigengenes ####


eigen = read.csv("Brain_module_eigengenes.csv")

module = eigen$ME3
Behavior = eigen$group
OvarySize = eigen$TotalOvaryIndex


#Plot ovary size vs module eigengenes
gp =ggplot(eigen, aes(x = OvarySize, y=eigen$ME3)) +
        geom_point(aes(fill = eigen$group ,shape = eigen$group), size = 4) +
        scale_shape_manual(values = c(21,22,23,24))+
        geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
        stat_cor(method = "spearman", color = 1)+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("r", "g", "s", "d"),
                       values=c("palegreen4", "purple", "blue", "red"))

#perform correlation test without plot to generate p values for each module

cor.test(x= OvarySize, y=eigen$ME1, method = "spearman")

#collected all Pvalues for FDR adjusment. Will adjust for all 26 comparisons to ovary size with brain and ovary data sets 
p =  c(0.0495,
       0.1822,
       0.006439,
       0.8152,
       0.1847,
       0.01922,
       0.9663,
       0.314,
       0.2707,
       0.6097,
       0.8815,
       0.07713,
       0.4299)

adjusted = p.adjust(p, method = "fdr", length(p))
adjusted


#### combined pvalues of 26 comparisons to ovary size for fdr adjustment ####
p =  c(0.0495,
       0.1822,
       0.006439,
       0.8152,
       0.1847,
       0.01922,
       0.9663,
       0.314,
       0.2707,
       0.6097,
       0.8815,
       0.07713,
       0.4299,
       0.3136,
       0.01207,
       0.02846,
       0.3168,
       0.00513,
       0.00000000008053,
       0.003737,
       0.9121,
       0.000004734,
       0.00000000006353,
       0.01146,
       0.009439,
       0.1949)

adjusted = p.adjust(p, method = "fdr", length(p))

#final pvalues for all comparisons to ovary size 
adjusted



################ Combined ovary and brain module eigengenes  ##########

#note: in the CSV brain modules are preceded by "B" ie brain module 3 is column "B3", ovary modules are preceded by "O"
eigen = read.csv("Brain_and_Ovaries_Combined_Eigengenes_and_Ovary_Size.csv")
Behavior = eigen$group


# plot of module eigengenes vs other eigengenes

gp =ggplot(eigen, aes(x=eigen$B3, y=eigen$O6)) +
        geom_point(aes(fill = Behavior ,shape = Behavior), size = 4) +
        scale_shape_manual(values = c(21,22,23,24))+
        geom_smooth(color = 1, method=lm, se=FALSE, fullrange=TRUE)+
        stat_cor(method = "spearman", color = 1)+
        theme_bw()+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
gp

gp + scale_fill_manual(breaks = c("r", "g", "s", "d"),
                       values=c("palegreen4", "purple", "blue", "red"))


##### module correlations with each other and with ovary size index ##

#read in and rename data
All = read.csv("Brain_ovary_eigengenes_combined.csv")

res1 <- cor.mtest(All)

m = cor(All, method = "pearson")

#make correlation plot with columns arranged by alphabet
corrplot(m, method = 'circle', order = "alphabet", addrect = 0, hclust.method = "ward.D2")

# make correlation plot with columns arranged by hierarchical clustering and pvalue = 0.05 (not fdr corrected). an "X" means pvalue is above p=0.05

m = cor(All, method = "spearman")
corrplot(m, method = 'circle', order = "hclust", addrect = 0, hclust.method = "ward.D2",p.mat = res1$p, sig.level = .05)

# make correlation plot without pvalue X's but with columns arranged by hierarchical clustering

m = cor(All, method = "spearman")
corrplot(m, method = 'circle', order = "hclust", addrect = 0, hclust.method = "ward.D2")

########## hypergeometric tests for gene list overlaps ############################

#Note: comparisons evaluated with a hypergeometric test are listed in the order that they are mentioned in the main text

#isolated subordinates vs subordinates overlap with Saleh and Ramirez, 2019 natural guards vs subordinates brains
phyper(49,3403,7095,88, lower.tail = FALSE)

#isolated subordinates vs subordinates overlap with Saleh and Ramirez, 2019 natural guards vs subordinates ovaries
phyper(41,2648,7166,67, lower.tail = FALSE)

# reproductive state contrast overlap with all DEGs detected in the ovaries among pairwise comparisons
phyper(318,514,9618,412, lower.tail = FALSE) 

# reproductive state contrast overlap with all DEGs detected in the brains among pairwise comparisons
phyper(64,132,10909,109, lower.tail = FALSE) 

# overlap between genes in ovary module 6 and brain module 10
phyper(71,1049,8661,350, lower.tail = FALSE) 

#overlap between genes upregulated in nonreproductive E. dilemma individuals and M. genalis workers
phyper(64,1647,4206,126, lower.tail = FALSE) 

#overlap between genes upregulated in nonreproductive E. dilemma individuals and nonlaying A. mellifera workers
phyper(35,481,5982,196, lower.tail = FALSE) 

#overlap between genes upregulated in reproductive E. dilemma individuals and M. genalis Queens
phyper(65,1984,3869,163, lower.tail = FALSE) 

#overlap between genes upregulated in reproductive E. dilemma individuals and laying A. mellifera workers
phyper(36,935,5528,227, lower.tail = FALSE) 

#overlap between ovary module 6 genes and DEGs between M. genalis queens and workers
phyper(159,3631,2222,227, lower.tail = FALSE)

#overlap between ovary module 6 genes and DEGs between A. mellifera nonlaying vs laying workers
phyper(89,1416,5047,266, lower.tail = FALSE)

#overlap between ovary module 10 genes and DEGs between M. genalis queens and workers
phyper(49,3631,2222,80, lower.tail = FALSE)

#overlap between ovary module 10 genes and DEGs between A. mellifera nonlaying vs laying workers
phyper(26,1416,5047,92, lower.tail = FALSE)

