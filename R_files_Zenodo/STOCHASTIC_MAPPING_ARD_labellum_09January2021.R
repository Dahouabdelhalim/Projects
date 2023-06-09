## first load packages
require(phytools)
## read data
X<-read.csv("labellum.csv",header=TRUE)
X
labellum <- X$labellum
names(labellum) <- X$Species
class(labellum)

labellum <- as.factor(labellum)
class(labellum)
## read tree
tree <- read.tree("C_tree.nwk")
tree

colr <- c("azure3", "gold", "chocolate1", "red")
names(colr) <- c(0,1,2,3)

## estimate ancestral states under an ARD model
fitARD<-ace(labellum,tree,model="ARD",type="discrete")
fitARD

## The element lik.anc gives us the marginal ancestral states, also known as the 'empirical Bayesian posterior probabilities.
fitARD$lik.anc

## It is fairly straightforwARD to overlay these posterior probabilities on the tree.
# 1. Open pdf file
pdf("labellum_ARD_1.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plotTree(tree,type="fan",fsize=0.7,ftype="i",lwd=1)
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitARD$lik.anc,piecol=colr,cex=0.4)
tiplabels(pie=to.matrix(labellum[tree$tip.label],
                        levels(labellum)),piecol=colr,cex=0.3)
add.simmap.legend(colors=colr,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)
# 3. Close the file
dev.off()

## An alternative tactic to the one outline above is to use an MCMC approach to sample character histories from their posterior probability labellum. 
## This is called stochastic character mapping (Huelsenbeck et al. 2003). 
## The model is the same but in this case we get a sample of unambiguous histories for our discrete character's evolution on the tree
## rather than a probability labellum for the character at nodes.
## simulate single stochastic character map using empirical Bayes method
mtreeARD<-make.simmap(tree,labellum,model="ARD")

# 1. Open pdf file
pdf("labellum_ARD_2.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(mtreeARD,colr,type="fan",fsize=0.7,ftype="i")
add.simmap.legend(colors=colr,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)
# 3. Close the file
dev.off()

## A single stochastic character map does not mean a whole lot in isolation - 
## we need to look at the whole labellum from a sample of stochastic maps. 
## This can be a bit overwhelming.

mtreesARD<-make.simmap(tree,labellum,model="ARD",nsim=1000)

# 1. Open pdf file
pdf("labellum_ARD_3.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
## For instance, the following code generates 100 stochastic character maps from our dataset and plots them in a grid
par(mfrow=c(10,10))
null<-sapply(mtreesARD,plotSimmap,colors=colr,lwd=1,ftype="off")
# 3. Close the file
dev.off()

## It is possible to summarize a set of stochastic maps in a much more meaningful way. 
## For instance, we can estimate the number of changes of each type, the proportion of time spent in each state, 
## and the posterior probabilities that each internal node is in each state, under our model. 
pdARD<-summary(mtreesARD)
pdARD

# 1. Open pdf file
pdf("labellum_ARD_4.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mfrow=c(1,1))
plot(pdARD,fsize=0.6,ftype="i",colors=colr,ylim=c(-2,Ntip(tree)))
add.simmap.legend(colors=colr,prompt=FALSE,x=0,y=-4,vertical=FALSE)
# 3. Close the file
dev.off()

## now let's plot a random map, and overlay the posterior probabilities
## Saving files in JPEG format
# 1. Open pdf file
pdf("labellum_ARD_5.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
plot(sample(mtreesARD,1)[[1]],colr,fsize=0.6,ftype="i",
     ylim=c(-2,Ntip(tree)))
nodelabels(pie=pdARD$ace,piecol=colr,cex=0.5)
add.simmap.legend(colors=colr,prompt=FALSE,x=0,y=-4,
                  vertical=FALSE)
# 3. Close the file
dev.off()


## For binary discrete characters, we can also use a method in phytools called densityMap 
## to visualize the posterior probability of being in each state across all the edges and nodes of the tree as follows:
# 1. Open pdf file
pdf("labellum_ARD_6.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
objARD<-densityMap(mtreesARD,states=levels(labellum)[2:1],plot=FALSE)
plot(objARD,fsize=c(0.6,1))
# 3. Close the file
dev.off()

## can generate a visualizing showing the posterior labellum of changes on the tree (rather than states). This is accomplished as follows:
# 1. Open pdf file
pdf("labellum_ARD_7.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
dotTree(tree,labellum,colors=colr,fsize=0.7,ftype="i",
        legend=FALSE)
add.simmap.legend(x=0,y=-4,colors=colr,prompt=FALSE,
                  vertical=FALSE,shape="circle")
nulo<-sapply(mtreesARD,markChanges,sapply(colr,
                                       make.transparent,0.1))
add.simmap.legend(colors=sapply(setNames(colr[2:1],
                                         c("1->2","2->1")),
                                make.transparent,0.1),prompt=FALSE,x=50,y=-4,
                  vertical=FALSE)
# 3. Close the file
dev.off()

## since we obtained these inferences under exactly the same model, let's compare the posterior probabilities from stochastic mapping with our marginal ancestral states.
## A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
## Saving files in PDF format
# 1. Open pdf file
pdf("labellum_ARD_8.pdf", width = 8.27, height = 11.69)
# 2. Create the plot
par(mar=c(4,4,3,3))
plot(fitARD$lik.anc,pdARD$ace,xlab="marginal ancestral states",
     ylab="posterior probabilities from stochastic mapping",
     pch=21,cex=1.4,bg="grey")
lines(c(0,1),c(0,1),lty="dashed",xlab="marginal ancestral states", ylab="posterior probabilities from stochastic mapping",col="red",lwd=2)
# 3. Close the file
dev.off()