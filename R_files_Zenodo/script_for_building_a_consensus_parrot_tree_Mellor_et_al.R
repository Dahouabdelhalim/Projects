#script for building consensus tree used for all but final hypothesis-testing analyses
library(ape)
library(caper)
library(phytools)

#read tree in 
parrottrees <- read.nexus("parrot_trees.nex")

#Note: most of this code is taken from the phytools blog post mentioned below, which is ran by
#Liam Revell. 

# MRC tree (from phytools blog, see here: http://blog.phytools.org/2016/03/method-to-compute-consensus-edge.html)
#p=0.5 specifies that the tree must be "majority rules consensus (MRC)"
consensus<-consensus.edges(parrottrees, consensus.tree=consensus(parrottrees,p=0.5)) 
plotTree(consensus, fsize=0.6)
#how does it look? too many zeros at the end? If so, you will need to run the code below....
print(consensus$edge.length) 
n <- length(consensus$edge.length)
#need to reduce the number of elements by 1... perhaps there is a bug in the code? We also add a very small amount
#to make sure no branch length is zero.
consensus$edge.length <- consensus$edge.length[1:n-1] + .0000001 
#how does it look now? One shorter, if you did it right...
print(consensus$edge.length) 
plotTree(consensus, fsize=0.6)

# First we take the consensus tree and make it ultrametric. Then we make that ultrametric tree into a 
# dichotomous ultrametric tree. Without doing this, the PhylANOVA and PIC analysis fails. I think this
# order makes the most sense. 
# Lambda is the rate-smoothing parameter. I think I want lambda to be small, so there is very little smoothing
#See here (citation at bottom): http://www.justinbagley.org/1226/update-new-functions-for-generating-starting-trees-for-beast-or-starbeast-in-r
is.ultrametric(consensus_ultra) #Check again, now it's ultrametric.
consensus_ultra=chronos(consensus, lambda=0) 
tree.unmatched <- multi2di(consensus_ultra, random=TRUE)

#Let's plot our new tree.
plotTree(tree.unmatched,fsize=0.6) 

tree.unmatched
#sorting out the polytomies in a random way
polysolved_tree<- multi2di(tree.unmatched, random = TRUE)

#now checking to see if tree now binary (comes back TRUE if so)
is.binary.tree(polysolved_tree)

#write into nexus file and export to folder
writeNexus(polysolved_tree, file = "consensus_parrot_tree.nex")
