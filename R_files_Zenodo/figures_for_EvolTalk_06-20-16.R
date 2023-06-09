
library(ape)
library(phangorn)
library(phytools)

# run E (morph-only MrBayes vs Cohen & Bitner

tree1<-prepTree(trees["runE_rynchComb_05-05-15_MrB_MajRule"][[1]])
tree1$edge<-tree1$edge[c(1:10,14:19,11:13,20:26),]
tree2<-treeMolOrig
#
resObj<-cophylo(tree1,prepTree(tree2))
plot.cophylo(resObj,pts=FALSE)
treeContradiction(tree1,tree2)


# run D (mol only) versus run F (combined)

tree1<-ladderize(trees["runD_rynchComb_05-05-15_MrB_MajRule"][[1]])
tree2<-ladderize(trees["runF_rynchComb_05-05-15_MrB_MajRule"][[1]])
#
resObj<-cophylo(prepTree(tree1),prepTree(tree2))
plot.cophylo(resObj,pts=FALSE)
treeContradiction(tree1,tree2)



# run E versus run F 

tree1<-prepTree(trees["runE_rynchComb_05-05-15_MrB_MajRule"][[1]])
tree1$edge<-tree1$edge[c(1:10,14:19,11:13,20:26),]
tree2<-ladderize(trees["runF_rynchComb_05-05-15_MrB_MajRule"][[1]])
#
resObj<-cophylo(tree1,prepTree(tree2))
plot.cophylo(resObj,pts=FALSE)
treeContradiction(tree1,tree2)



# missing data

if(!identical(names(propMissingMol),names(propMissingMorph))){
  stop("propMissingMol and propMissingMorph have different taxon names")
}
oldPar<-par(no.readonly = T)
layout(matrix(1:2,1,2),
      width=c(2,1))
par(mar=c(4,11,1,1))
barplot(rev(propMissingMol*100),xlab="% Missing Data",
        horiz=TRUE,las=2,main="Molecular")
par(mar=c(4,0,1,1))
barplot(rev(propMissingMorph*100),xlab="% Missing Data",
        horiz=TRUE,las=2,names="",main="Morph")
layout(1)
#mtext(side=1,line=2.5,"    % Missing Data",cex=1.5)
par(oldPar)