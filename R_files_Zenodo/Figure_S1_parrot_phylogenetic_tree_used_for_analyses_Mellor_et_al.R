#script for Figure S1 (parrot tree used in analyses)

library(ggplot2)
library(ggtree)
library(ape)
library(extrafont)
library(RColorBrewer)
#read in tree
parrottree <- read.nexus("consensus_parrot_tree.nex")
#now read data in
parrot_data<-read.csv("parrot_comp_data_tree.csv", header=TRUE)
attach(parrot_data)

#prune tree so it's just species we have outcome data for
pruned.tree<-drop.tip(parrottree,parrottree$tip.label[-match(Species_name, parrottree$tip.label)])
#have a look
plot(pruned.tree)
#make into a ggtree object
p<-ggtree(pruned.tree, size=1.5, aes(color=Outcome_data), continuous=FALSE)
p <- p%<+% parrot_data

#start adding aesthetics
tree<-p+theme_tree(bgcolor="white")+xlim(0, 1.2)
#make an open fan shape to better see the species names. This also sorting what the tip
#labels look like
tree1<-open_tree(tree, 0)+geom_tiplab2(aes(label=Species_lab), fontface="bold.italic", 
                                       family="Calibri", size=3, color="black", offset=0.01)
#now add colour to the branch lengths to show which species have data for SB only (green),
#hatch rate data only (blue) and both types (orange)
tree1<-tree1+scale_color_brewer(type="qual", palette = "Dark2", na.value = "gray34")+
                                      theme(legend.position = "none")

#save as a file
png(filename="parrot_tree1.png", type="cairo", units="in",
    width=15, height=8, pointsize = 20, res = 500)
plot(tree1)
dev.off()




  

