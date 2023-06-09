# This R script tests whether pollination syndromes predict 
# floral trait factor analysis dimensionss 1 - 10 using phylogenetic 
# Anova.  It also makes a figure projecting phylogeny onto species 
# values for floral trait dimensions 1 and 2 (Figure 5). It uses the 
# phylogenetic tree file "chronogram_v8_2.tre" and the floral trait factor 
# analysis output data file "famd_coords.csv"

library(geiger) # match tree and data
library(phytools)
library(ggplot2)

setwd("~/Dropbox/Costus grant/Dena and Kathleen/FloralEvolution/Dryad_Deposit")

#chronogram
tree <- read.nexus("chronogram_v8_2.tre")
plot(tree, cex=0.4)

#floral dim characters
data <- read.csv("famd_coords.csv",row.names=1)

# match tree and data
(name.check(tree, data) -> phyOverlap)
drop.tip(tree, phyOverlap$tree_not_data) -> phyComparativeTree
data <- data[!(data[,1]  %in% phyOverlap$data_not_tree), ]
data<- data[,1:11]
names(data) <- c("syndrome", "Dim.1", "Dim.2","Dim.3", "Dim.4","Dim.5", "Dim.6","Dim.7", "Dim.8","Dim.9", "Dim.10")

# use random forest predicted syndromes
data["Costus_dirzoi_98079","syndrome"]<-"hummingbird"
data["Costus_sp_nov_19168","syndrome"]<-"hummingbird"
data["Costus_varzearum_19252","syndrome"]<-"hummingbird"

# Phy ANOVAs Dim1 by syndrome
dat = data$Dim.1
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant
# Analysis of Variance Table
#
# Response: dat
# Df Sum-Sq Mean-Sq F-value     Pr(>F) Pr(>F) given phy    
# group      1 350.98  350.98  126.11 2.8149e-15         0.000999 ***
#  Residuals 50 139.16    2.78       

# Phy ANOVAs Dim2 by syndrome
dat = data$Dim.2
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant
# Analysis of Variance Table
#
# Response: dat
# Df Sum-Sq Mean-Sq F-value  Pr(>F) Pr(>F) given phy
# group      1   7.90  7.9003  2.2991 0.13575           0.3626
# Residuals 50 171.82  3.4363               

# Phy ANOVAs Dim3 by syndrome
dat = data$Dim.3
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant

# Phy ANOVAs Dim4 by syndrome
dat = data$Dim.4
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant

# Phy ANOVAs Dim5 by syndrome
dat = data$Dim.5
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant

# Phy ANOVAs Dim6 by syndrome
dat = data$Dim.6
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant

# Phy ANOVAs Dim7 by syndrome
dat = data$Dim.7
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant

# Phy ANOVAs Dim8 by syndrome
dat = data$Dim.8
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant

# Phy ANOVAs Dim9 by syndrome
dat = data$Dim.9
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant

# Phy ANOVAs Dim10 by syndrome
dat = data$Dim.10
dat = as.data.frame(dat)
rownames(dat) = rownames(data)
grp = as.factor(data$syndrome)
names(grp)=rownames(data)
t=aov.phylo(dat ~ grp, phyComparativeTree) #not significant

##############
#FIGURE 5
##############

#get colors to paint tips
tip.col <- data[,1]
tip.col[tip.col=="bee"] <- "#1F639B"
tip.col[tip.col=="hummingbird"] <- "#ED553B"
names(tip.col)<- row.names(data)
cols<-c(tip.col[phyComparativeTree$tip.label],rep("black",phyComparativeTree$Nnode))
names(cols)<-1:(length(phyComparativeTree$tip)+phyComparativeTree$Nnode)


pdf(file="Figure5.pdf",width=5,height=5.5,paper='special') 
phylomorphospace(phyComparativeTree, data[,c(2,3)], ftype="off",node.by.map=TRUE,
                 xlab="Dim1 (33.5% PVE)", ylab="Dim2 (12.4 PVE)", bty="l",control=list(col.node=cols),
                 node.size=c(0,1) )
legend(x=-5, y=4, legend=c("Hummingbird","Bee"), 
       title="Pollination syndrome",col=c("black"), pch=21, pt.bg=c("#ED553B","#1F639B"), cex=0.7)

#add X symbol to taxa where RF syndrome was different than apriori
data.sub <- data[c("Costus_dirzoi_98079","Costus_sp_nov_19168","Costus_varzearum_19252"),]
points(data.sub$Dim.1,data.sub$Dim.2, pch=4, cex=0.45)

dev.off()


