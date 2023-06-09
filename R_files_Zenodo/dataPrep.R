# Data Prep Subsetting Data
# Author: Sathvik Palakurty
# Afkhami Lab

library("WGCNA")
library(devtools)
library(Biobase)
library(preprocessCore)
library(DESeq2)

enableWGCNAThreads()


# Read the raw count data
exprData1 <- read.csv("D:/Google Drive/Coexpression Network Analysis/Plants Raw Count Data/exprData1.txt", row.names=1, sep="")
coldata <- read.csv("D:/Google Drive/Coexpression Network Analysis/Plants Raw Count Data/coldata.csv", row.names=1)

#Subset the data into the individual treatment groups
exprBoth <- exprData1[c("M.R.21","M.R.40","M.R.42", "M.R.69","M.R.82","M.R.95","M.R.101","M.R.110")]
exprM <- exprData1[c("M.R.23","M.R.30","M.R.39","M.R.66","M.R.90","M.R.106","M.R.118")]
exprR <- exprData1[c("M.R.5","M.R.28","M.R.44","M.R.52","M.R.73","M.R.88","M.R.104","M.R.109")]
exprNone <- exprData1[c("M.R.2","M.R.43","M.R.47","M.R.51","M.R.77","M.R.91","M.R.113","M.R.116")]

#Find the genes that satisfy the filtering criterion for each treatment group
#assign to individual booleans
boolean <- c()
for(i in 1:nrow(exprBoth))
{
	if(length(which(exprBoth[i,] > 1.0)) >= 1.0)
	{
		boolean <- c(boolean, TRUE)
	}
	else
	{
		boolean <- c(boolean, FALSE)
	}
}
booleanBoth <- boolean

boolean <- c()
for(i in 1:nrow(exprM))
{
	if(length(which(exprM[i,] > 1.0)) >= 1.0)
	{
		boolean <- c(boolean, TRUE)
	}
	else
	{
		boolean <- c(boolean, FALSE)
	}
}
booleanM <- boolean

boolean <- c()
for(i in 1:nrow(exprR))
{
	if(length(which(exprR[i,] > 1.0)) >= 1.0)
	{
		boolean <- c(boolean, TRUE)
	}
	else
	{
		boolean <- c(boolean, FALSE)
	}
}
booleanR <- boolean

boolean <- c()
for(i in 1:nrow(exprNone))
{
	if(length(which(exprNone[i,] > 1.0)) >= 1.0)
	{
		boolean <- c(boolean, TRUE)
	}
	else
	{
		boolean <- c(boolean, FALSE)
	}
}
booleanNone <- boolean

# Combine all the treatment groups into a data frame to choose the genes that will
# go into the final analysis
booleanDF <- data.frame(booleanBoth, booleanM, booleanR, booleanNone)
boolean <- c()
for(i in 1:nrow(booleanDF))
{
    if(length(which(booleanDF[i,] == TRUE)) >= 1)
    {
        boolean <- c(boolean, TRUE)
    }
    else
    {
        boolean <- c(boolean, FALSE)
    }
 }

# Subset the original data and reassign the filtered data to the treatment group datasets
exprData1Subset <- exprData1[boolean,]
exprData1Subset <- exprData1Subset * 10
DESeqData <- DESeqDataSetFromMatrix(countData = exprData1Subset, colData = coldata , design = ~Myco + Rhizo , tidy = FALSE)
transformedData <- rlogTransformation(DESeqData)
exprData1SubsetTransformed <- assay(transformedData)
exprData1SubsetTransformedNormed <- normalize.quantiles(exprData1SubsetTransformed)

row.names(exprData1SubsetTransformedNormed) <- row.names(exprData1Subset)
colnames(exprData1SubsetTransformedNormed) <- colnames(exprData1Subset)
exprData1SubsetTransformedNormed <- as.data.frame(exprData1SubsetTransformedNormed)

exprBothTSubset <- t(exprData1SubsetTransformedNormed[c("M.R.21","M.R.40","M.R.42", "M.R.69","M.R.82","M.R.95","M.R.101","M.R.110")])
exprMTSubset <- t(exprData1SubsetTransformedNormed[c("M.R.23","M.R.30","M.R.39","M.R.66","M.R.90","M.R.106","M.R.118")])
exprRTSubset <- t(exprData1SubsetTransformedNormed[c("M.R.5","M.R.28","M.R.44","M.R.52","M.R.73","M.R.88","M.R.104","M.R.109")])
exprNoneTSubset <- t(exprData1SubsetTransformedNormed[c("M.R.2","M.R.43","M.R.47","M.R.51","M.R.77","M.R.91","M.R.113","M.R.116")])
exprData1SubsetTransformedNormedTransposed <- t(exprData1SubsetTransformedNormed)

save(exprData1SubsetTransformedNormedTransposed, file = "exprAll.RData")
save(exprBothTSubset, file = "exprBothTSubset.RData")
save(exprNoneTSubset, file = "exprNoneTSubset.RData")
save(exprMTSubset, file = "exprMTSubset.RData")
save(exprRTSubset, file = "exprRTSubset.RData")