library("ctv")
install.views("Phylogenetics")
library(ape)
library(caper)
library(geiger)
library(nlme)
library(phangorn)
library(phytools)
library(phylolm)
library(corHMM)
library(dplyr)
library(tidyverse)
library(readxl)
library(tidyr)
library(AICcmodavg)
library(Hmisc)
library(jtools)
library(binom)

#Reading in the species data#
Rallids<-read.csv("Data/RallidPhylogenyData.csv",header=TRUE,row.names=1)
Rallids
str(Rallids)

Rallid_CSV<-read.csv("Data/RallidLogData.csv",header=TRUE,row.names=1)
Rallid_CSV
str(Rallid_CSV)

head(Rallids)
summary(Rallids)

#Reading in the rooted tree#
RallidTree<-read.tree("Data/Rallid_Tree.tre")
RallidTree
plotTree(RallidTree,ftype="i",fsize=0.6,lwd=1,par(family="serif"))
Ntip(RallidTree)

#Checking for polytomies#
is.binary.tree(RallidTree)

#Checking for ultramectricity
is.ultrametric(RallidTree)
UltraTree<-force.ultrametric(RallidTree)
is.ultrametric(UltraTree)
print(UltraTree$edge.length)
plotTree(UltraTree,ftype="i",fsize=0.6,lwd=1,par(family="serif"))

#Comparing species in tree to data
NameCheck<-name.check(RallidTree,Rallids)
NameCheck

#Phylogenetic principal component analysis

#Drop categorical variables
RallidsCont <- Rallid_CSV %>% select(-c(Habitat:LogMass))
RallidsCont
str(RallidsCont)
ContData <- treedata(RallidTree, RallidsCont, sort = TRUE)
physignal(ContData$data, ContData$phy, 
          iter = 10000, print.progress = FALSE)

#Run the PCA without phylogenetic correction
Rallids_pca<- prcomp(RallidsCont, scale = TRUE)
Rallids_pca
summary(Rallids_pca)
Rallids_pca$rotation <- -1*Rallids_pca$rotation
summary(Rallids_pca$rotation)
Rallids_pca$x <- -1*Rallids_pca$x
head(Rallids_pca$x)
biplot(Rallids_pca)
round(Rallids_pca$sdev^2,2)

#View the structure and calculate percent total variance of each component
Rallids_pca$sdev^2 / sum(Rallids_pca$sdev^2/100)

#Run the PPCA after choosing the appropriate model
Rallids_non_rotated_ppca<- phyl.pca(RallidTree, RallidsCont, method="lambda", mode="cov")
Rallids_non_rotated_ppca
summary(Rallids_non_rotated_ppca)

#View the structure and calculate percent total variance of each component
str(Rallids_non_rotated_ppca)
percent.var.ppca <-(diag(Rallids_non_rotated_ppca$Eval)/sum(diag(Rallids_non_rotated_ppca$Eval)) * 100)
percent.var.ppca

# A scree plot is useful for understanding how variance is distributed among
#the principal components, and it should be the first step in analyzing a PCA.
#The scree plot is particularly critical for determining how many principal 
#components should be interpreted. Although this could be done by calling plot(pca),
#a better-annotated plot that plots percent of total variance for each principal
#component can be made as follows.

plot(Rallids_pca, type="l")
plot(Rallids_non_rotated_ppca, type="l")

barplot(percent.var.ppca, xlab='PC', ylab='Percent Variance',
names.arg=1:length(percent.var.ppca), las=1, col='gray')
#One guideline for the number of principal components to use is to
#accept all principal components that explain more than one variable's worth of data.
#If all the variables contributed the same variance, this cutoff would be 1/p, where
#p is the number of variables.
abline(h=1/ncol(RallidsCont)*100, col="red")

#A table of loadings should be examined next, as it shows which variables have
#high loadings (positive or negative) on each principal component, that is,
#which variables contribute most strongly to each PC. Examining this table can give
#you a good sense of what each principal component represents, in terms of the
#original data. A positive loading means that a variable correlates positively with
#the principal component; a negative loading indicates a negative correlation.
#Rounding loadings to 2-3 decimal places often makes these tables easier to parse.

# Since we will be considering only the first 2-3 PCs, we will display only those.
loadings0<-round(Rallids_pca$rotation, 2)[ , 1:3]
loadings0
loadings<-round(Rallids_non_rotated_ppca$L, 2)[ , 1:3]
loadings
 
### Note that there is a strong negative loading for frequency variables on PC1.
#Higher PC values correspond to lower values of the variables. 
#On PC2, Call length has a strong positive loading. 
#On PC3 the bandwidth variables have strong positive loadings.

#Get the loadings and Eigenvalues 
write.csv(Rallids_non_rotated_ppca$L, "Rallids Phylogenetic PCA Loadings.csv")
write.csv(Rallids_non_rotated_ppca$Eval, "Rallids Phylogenetic PCA Eigenvalues.csv")

#Get the pca scores for each taxon
Rallids_pca_scores <- as.data.frame(Rallids_pca$x)
Rallids_pca_scores$names <- row.names(Rallids_pca_scores)
head(Rallids_pca_scores)

Rallids_ppca_scores <- as.data.frame(Rallids_non_rotated_ppca$S)
Rallids_ppca_scores$names <- row.names(Rallids_ppca_scores)
head(Rallids_ppca_scores)

Rallids_pca_scores<-cbind(Rallid_CSV,Rallids_pca_scores)
Rallids_ppca_scores<-cbind(Rallid_CSV,Rallids_ppca_scores)

#Plot the first two components
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
plot(Rallids_ppca_scores$PC1,
     Rallids_ppca_scores$PC2,
     main="Scatterplot of Rallid PPCA Scores",
     xlab="PC1",
     ylab="PC2",
     pch=19)

#Add loadings
   scale <- 2.2
   arrows(0, 0, loadings[, 1]*scale, loadings[, 2]*scale,
   length=0.1, angle=20, col='red')
   
   labelScale <- 1.1
   text(loadings[, 1]*scale*labelScale, loadings[, 2]*scale*
   labelScale, rownames(loadings), col='red', cex=0.7)
   
#Look at some grouping, e.g. habitat
text(Rallids_pca_scores$PC1, Rallids_pca_scores$PC2+0.01, 
       labels=Rallids_pca_scores$Habitat, col= "green", cex= 0.7, pos= 4)

text(Rallids_pca_scores$PC1, Rallids_pca_scores$PC2+0.01, 
     labels=Rallids_pca_scores$Duets, col= "red", cex= 0.7, pos= 4)

### Run Phylogenetic MANOVA in geiger
dat_pca=Rallids_pca_scores[,19:27]
datshort_pca=Rallids_pca_scores[,19:21]
dat_pca
datshort_pca

dat=Rallids_ppca_scores[,19:27]
datshort=Rallids_ppca_scores[,19:21]
dat
datshort
#If you are doing this for the PCs, you specify the relevant columns too
d1=dat[,1]
d2=dat[,2]
d3=dat[,3]
#To specify one particular vector for ANOVA
names(d1)=rownames(dat)
names(d2)=rownames(dat)
names(d3)=rownames(dat)

#Habitat Effects on First Three PCAs Only
hgrp_0<-as.factor(Rallids_pca_scores$Habitat)
names(hgrp_0)=rownames(datshort_pca)
hab0=aov.phylo(datshort_pca~hgrp_0, RallidTree, nsim=1000, test="Wilks")
### MANOVA of habitat on PC scores 1-3.
summary(hab0)
summary.aov(hab0)

#Duet Effects on First Three PCAs Only
dgrp_0<-as.factor(Rallids_pca_scores$Duets)
names(dgrp_0)=rownames(datshort_pca)
duets0=aov.phylo(datshort_pca~dgrp_0, RallidTree, nsim=1000, test="Wilks")
### MANOVA of duets on PC scores 1-3.
summary(duets0)
summary.aov(duets0)

#Habitat Effects on First Three PPCAs Only
hgrp<-as.factor(Rallids_ppca_scores$Habitat)
names(hgrp)=rownames(datshort)
hab=aov.phylo(datshort~hgrp, RallidTree, nsim=1000, test="Wilks")
### MANOVA of habitat on PC scores 1-3.
summary(hab)
summary.aov(hab)

#Duet Effects on First Three PPCAs Only
dgrp<-as.factor(Rallids_ppca_scores$Duets)
names(dgrp)=rownames(datshort)
duets=aov.phylo(datshort~dgrp, RallidTree, nsim=1000, test="Wilks")
### MANOVA of duets on PC scores 1-3.
summary(duets)
summary.aov(duets)

### Comparison of ANOVAs of habitat on PC scores 1-3.
d1anov=aov.phylo(d1~hgrp, RallidTree, nsim=50, test="Wilks")  ### ANOVA on PC1
d2anov=aov.phylo(d2~hgrp, RallidTree, nsim=50, test="Wilks")  ### ANOVA on PC2
d3anov=aov.phylo(d3~hgrp, RallidTree, nsim=50, test="Wilks")  ### ANOVA on PC3

### ANOVA of habitat on PC scores 1-3, with Holm correction in post hoc tests.
phylANOVA(RallidTree,hgrp,d1,posthoc=TRUE,p.adj="holm")
phylANOVA(RallidTree,hgrp,d2,posthoc=TRUE,p.adj="holm")
phylANOVA(RallidTree,hgrp,d3,posthoc=TRUE,p.adj="holm")
### ANOVA of duets on PC scores 1-3, with Holm correction in post hoc tests.
phylANOVA(RallidTree,dgrp,d1,posthoc=TRUE,p.adj="holm")
phylANOVA(RallidTree,dgrp,d2,posthoc=TRUE,p.adj="holm")
phylANOVA(RallidTree,dgrp,d3,posthoc=TRUE,p.adj="holm")