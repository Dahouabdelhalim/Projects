
##### install and load required packages

#install.packages("remotes")
#remotes::install_github("fmichonneau/phyloch")
library(phytools)
library(phyloch)
library(strap)

# # # # #     other possible option for plotting geo scale could be with packages ggtree and deeptime
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("ggtree")
#library(devtools)
#install_github("willgearty/deeptime")


# set working directory
setwd("/path/to/BEAST/output/")

# load the tree with phyloch package
phy<-phyloch::read.beast("out_mean_20p_MCC.tre")

# set the root age with function nodeHeights from phytools package
phy$root.time <- max(nodeHeights(phy))

# when exporting tree as pdf
pdffn = "BEASTtreeWithGeologicalScale.pdf"
pdf(pdffn, width=100, height=200)

# plot the tree with geological time scale (from package strap)
geoscalePhylo(tree = ladderize(phy, right = TRUE), units = c("Period"), boxes = "Period", tick.scale=50, font =1, cex.tip = 1.5, cex.age = 4, cex.ts = 8, width = 5, x.lim = c(-30,350), quat.rm=TRUE)

# plot node bars on the tree
HPDbars(phy, label = "height_95%_HPD", tab = NULL, col = "blue", lwd=8, broken = FALSE)

# required for PDFs
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)


