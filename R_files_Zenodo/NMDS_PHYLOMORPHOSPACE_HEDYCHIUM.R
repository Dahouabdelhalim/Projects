##load the packages
library(vegan)
library(phytools)

##load the tree file and the data file
morpho_hedy <- read.csv("H:/000_Since_October_2020/000_Floral_Evolution/0_Floral_Evolution_Hedychium/Complete_21.csv", header= T, row.names = 1)
h_tree <- read.tree("H:/000_Since_October_2020/000_Floral_Evolution/0_Floral_Evolution_Hedychium/C_tree.nwk")
plotTree(h_tree)

##perform NMDS analysis
morpho_nmds <- metaMDS(morpho_hedy, k=3)

##perform phylomorphospace analysis
phylomorpho <- phylomorphospace(h_tree, morpho_nmds$points)
phylomorpho3d <- phylomorphospace3d(h_tree, morpho_nmds$points)

##Plotting the phylomorpho
# 1. Open pdf file
pdf("phylomorpho.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
phylomorpho <- phylomorphospace(h_tree, morpho_nmds$points)
# 3. Close the file
dev.off()

##Plotting the phylomorpho3d
# 1. Open pdf file
pdf("phylomorpho3d.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
phylomorpho3d <- phylomorphospace(h_tree, morpho_nmds$points)
# 3. Close the file
dev.off()
