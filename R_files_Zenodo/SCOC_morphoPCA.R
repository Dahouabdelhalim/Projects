## Sceloporus Morphological Analysis
## Written by Aryeh H. Miller, edited by Hayden R. Davis on 5 October 2021

#Load relevant packages
library(viridis)
library(RColorBrewer)
library(ggsci)
library(tidyverse)
library(viridis)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(car)
library(ggplot2)

# Load raw data of all mature individuals
setwd("~/<directory>")
scelops <- na.exclude(read.csv(file = "SCOC_morphological_data.csv"))
head(scelops)
nrow(scelops)
ncol(scelops)

#Mann-Whitney U Test on sexual dimorphism in raw data
wilcox.test(svl ~ sex, data = scelops)

#normalize data when dealing with meristic data
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

norm_scelop <- as.data.frame(lapply(scelops[11:16], min_max_norm))
scelops_new <- cbind(scelops, norm_scelop)
scelops_new <- scelops_new[,c(1:10, 17:22)]
head(scelops_new)

#Log transform data
TlogSVL <- log(scelops_new$svl)
TlogHeadL <- log(scelops_new$hl)
TlogHeadW <- log(scelops_new$hw) 
TlogLFing <- log(scelops_new$lfl) 
TlogRFing <- log(scelops_new$rfl) 
TlogLLongToe <- log(scelops_new$lrtl)
TlogRLongToe <- log(scelops_new$lltl)

#Correct for SVL and generate residuals. If you want to separate Male/Female, do this same 
#correction for both groups
TlogHeadLc <- residuals(lm(TlogHeadL ~ TlogSVL))
TlogHeadWc <- residuals(lm(TlogHeadW ~ TlogSVL))
TlogLFingc <- residuals(lm(TlogLFing ~ TlogSVL))
TlogRFingc <- residuals(lm(TlogRFing ~ TlogSVL))
TlogLLongTc <- residuals(lm(TlogLLongToe ~ TlogSVL))
TlogRLongTc <- residuals(lm(TlogRLongToe ~ TlogSVL))

#subset data that will not be size corrected
ms <- scelops_new[,15]
ds <- scelops_new[,16]
lfp <- scelops_new[,11]
rfp <- scelops_new[,12]
ltl <- scelops_new[,13]
rtl <- scelops_new[,14]

head(scelops_new)

#Combine into new dataframe
datameasurements <-  droplevels(subset(scelops_new, 
      select = c(id, sex, population,  svl,   hl,   hw,  lfl,  rfl, lrtl, lltl, ltl, rtl, rfp,
                 lfp, ms, ds)))

TTable <-  cbind(datameasurements,TlogSVL, TlogHeadL,TlogHeadW,
      TlogLFing,TlogRFing, TlogLLongToe, TlogRLongToe, 
      ltl, rtl, ms, ds, TlogHeadLc ,TlogHeadWc,TlogLFingc,TlogRFingc, 
      TlogLLongTc, TlogRLongTc)

head(TTable)

#Add species labels
PopLabels <- c("PUGn", "PUGs", "KITs", "KITw", "OLY")
head(TTable)

#PCA
pca_data  <-  TTable
pops <- TTable[,3] #Isolate species
data3 <- cbind(pops, TTable[24:33]) #Dataframe with OTU residual data and species identities
pca <- prcomp(data3[ ,2:ncol(data3)], scale = T) #Perform PCA
scores <- data.frame(pops, pca$x[ ,1:2]) #Dataframe for PC1 and PC2
summary(pca)

#Assemble single table
loadings <- pca$rotation
temp <- summary(pca) 
sum <- temp$importance 
eigen <- pca$sdev^2 
sumM <- rbind(sum, eigen, loadings) #Combine into single table
View(sumM) #View table
write.table(sumM, "pcaSummary_meristic_corr.csv", sep = "\\t")

#Extract PC1 and PC2 scores for plot
PC1 = scores$PC1
PC2 = scores$PC2

pmainM <- ggplot(scores, aes(x = PC1, y = PC2, group = pops)) +
    theme_bw() + theme(legend.key = element_blank())  + 
    scale_color_manual(values = c("indianred4", "darkgreen", "steelblue3", "mediumorchid4", "orangered2"))  + 
    geom_point((aes(color = pops)), size = 6) + 
    xlab("PC1") + ylab("PC2") + theme(legend.title = element_blank(), 
    axis.text = element_text(size = 12), axis.title = element_text(size = 12),
    legend.text = element_text(size = 12)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + labs(x = "PC1 (35.8%)", y = "PC2 (20.7%)")
pmainM

