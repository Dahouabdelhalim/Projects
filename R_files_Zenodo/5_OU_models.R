# This R script detects location of evolutionary shifts in the 
# optima of floral trait factor analysis dimension 1 and 2, and 
# tests for convergent regimes, i.e., convergent optima. It uses 
# the phylogenetic tree file "chronogram_v8_2.tre" and the floral 
# trait factor analysis output data file "famd_coords.csv"

library(geiger) # match tree and data
#library(devtools)
#install_github("glmgen/genlasso")
#install_github("khabbazian/l1ou")
library(l1ou)

setwd("~/Dropbox/Costus grant/Dena and Kathleen/FloralEvolution/Dryad_Deposit")

# chronogram
tree <- read.nexus("chronogram_v8_2.tre")
plot(tree, cex=0.4)

# floral Dim characters
data <- read.csv("famd_coords.csv",row.names=1)

# match tree and data
(name.check(tree, data) -> phyOverlap)
drop.tip(tree, phyOverlap$tree_not_data) -> phyComparativeTree
data <- data[!(row.names(data)  %in% phyOverlap$data_not_tree), ]

# Prep tree and data to meet requirements for estimate_shift_configuration
#dim1
dat <- as.data.frame(data[,2]) #dim 1 -5 (each explains >5% of variation in floral space)
rownames(dat) <- rownames(data)
myadj.dim1 <- adjust_data(phyComparativeTree, dat,normalize = TRUE, quietly = FALSE)
#dim2
dat <- as.data.frame(data[,3]) #dim 1 -5 (each explains >5% of variation in floral space)
rownames(dat) <- rownames(data)
myadj.dim2 <- adjust_data(phyComparativeTree, dat,normalize = TRUE, quietly = FALSE)

#detect evolutionary shifts under OU model and get convergent regimes
#Dim1
(eModel <- estimate_shift_configuration(myadj.dim1$tree, myadj.dim1$Y, 
                criterion="AICc", alpha.upper = 100))
(fit_conv <- estimate_convergent_regimes(eModel,criterion="pBIC"))
pdf("FigOUConvergenceDim1.pdf", width=10, height=10)
plot(fit_conv)
dev.off()

#Dim2
(eModel <- estimate_shift_configuration(myadj.dim2$tree, myadj.dim2$Y, criterion="AICc", alpha.upper=500))
(fit_conv <- estimate_convergent_regimes(eModel,criterion="AICc"))
pdf("FigOUConvergenceDim2.pdf", width=10, height=10)
plot(fit_conv)
dev.off()
