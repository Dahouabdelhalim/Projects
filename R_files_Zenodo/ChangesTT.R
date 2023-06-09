library(phytools)
library(geiger)
packageVersion("phytools")
setwd("/PATH/TO/FOLDER")
#Import data
tree <- read.nexus("rescale.tree")
habitat <- read.csv("Habitat.csv")
habitat.simmap <- habitat[,2]
names(habitat.simmap) <- habitat[,1]
tree <- ladderize(treedata(tree,habitat.simmap)$phy)

#match tips and states; define state colors
x <- habitat.simmap[tree$tip.label]


#Simulate stochastic character maps on a phylogenetic tree
trees<-make.simmap(tree,x,model="SYM",nsim=100)

#Generate (or simulates) a 'changes through time' plot from a set of stochastic map character histories

object<-ctt(trees,segments = 13)
object

#use the function sim.multiCtt to simulate various rather than a single CTT

Q<-trees[[1]]$Q
Q

nulo<-sim.multiCtt(tree,Q,nsim=100)

#Plot results
pdf("changes_Null_number.pdf")
plot(nulo,type="number")
plot(object,add=TRUE,type="number")
plotTree(ladderize(tree),add=TRUE,ftype="off",lwd=1,color=make.transparent("blue",0.1),
         mar=par()$mar)
dev.off()

pdf("changes_Null_rate.pdf")
plot(nulo,alpha=0.2,ylim=c(0,0.075))
plot(object,add=TRUE, type="l")
plotTree(ladderize(tree),add=TRUE,ftype="off",lwd=1,color=make.transparent("blue",0.1),
         mar=par()$mar)
dev.off()
