# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #
#																					  					 #
#	Supplementary R code for:														  					 #
#																					  					 #
#				GROUP FORMATION AND THE EVOLUTIONARY PATHWAY TO SOCIAL COMPLEXITY IN BIRDS				 #
#																					  					 #
#	Philip A. Downing, Ashleigh S. Griffin and Charlie K. Cornwallis				  						 #
#																					  					 #
#	contact: philip.downing@biol.lu.se												  					 #
#																					  					 #
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

# packages

library(ape)
library(coda)
library(MCMCglmm)
library(QuantPsyc)
library(doBy)
library(nnet)
library(phytools)
library(metafor)
library(reshape2)


# functions

AncMultTree <- function(responselevels=NULL, model=NULL, trees=NULL, starttree=NULL, endtree=NULL){
Res <- data.frame(tree=vector(), ancestors=character(), descendents=character(), edge.length=numeric(), probAncNCoop=numeric(), probAncFam=numeric(), probAncNonFam=numeric(), probDesNCoop=numeric(), probDesFam=numeric(), probDesNonFam=numeric())
IJ <- (1/3) * (diag(2) + matrix(1, 2, 2))
Res2 <- data.frame(matrix(nrow=0,ncol=4))
colnames(Res2) <- c("node","X1","X2","X3")
nonodes <- (length(model$Sol[1,])-(responselevels-1))/2
ob1 <- model$Sol[,responselevels:length(model$Sol[1,])]
trees <- trees[starttree:endtree]
for(i in 1:nonodes){
  ob2 <- mcmc(cbind(ob1[,i],ob1[,nonodes+i]))
  Delta <- cbind(c(-1, 1, 0), c(-1, 0, 1))
  c2 <- (16 * sqrt(3)/(15 * pi))^2
  D <- ginv(Delta %*% t(Delta)) %*% Delta
  Int <- t(apply(ob2, 1, function(x) {
    D %*% (x/sqrt(1 + c2 * diag(IJ)))
  }))
  probs <- data.frame(node=colnames(ob1)[i], exp(Int)/rowSums(exp(Int)))
  probs$node <- gsub("traitgroupFormation.nonCoop.animal.", "", probs$node)
  Res2 <- data.frame(rbind(Res2, probs))  
}
Res2$tree <- rep(1:length(model$Sol[,1]), nonodes) 
Res2$node <- as.factor(Res2$node)
Res2$code <- paste(Res2$tree, Res2$node,sep=".")
for(i in 1:length(model$Sol[,1])) {
  tree <- trees[[i]]
  TreeDat <- as.data.frame(tree$edge)
  TreeDat$edge.length <- tree$edge.length
  TreeDat$V1name <- paste("animal.Node",(TreeDat$V1-length(tree$tip.label)),sep="")
  TreeDat$V2name <- paste("animal.Node",(TreeDat$V2-length(tree$tip.label)),sep="")
  treesp <- data.frame(tip.label=paste("animal.",tree$tip.label,sep=""),no=1:length(tree$tip.label))
  TreeDat <- data.frame(TreeDat, ancestors=treesp$tip.label[match(TreeDat$V1,treesp$no)], descendents=treesp$tip.label[match(TreeDat$V2,treesp$no)])
  TreeDat$ancestors <- as.character(TreeDat$ancestors)
  TreeDat$descendents <- as.character(TreeDat$descendents)
  TreeDat$ancestors <- as.character(ifelse(!is.na(TreeDat$ancestors), TreeDat$ancestors, TreeDat$V1name))
  TreeDat$descendents <- as.character(ifelse(!is.na(TreeDat$descendents),TreeDat$descendents,TreeDat$V2name))
  TreeDat$ancestors <- as.vector(gsub("animal.","",TreeDat$ancestors))
  TreeDat$descendents <- as.vector(gsub("animal.","",TreeDat$descendents))   
  TreeDat <- data.frame(tree=i, ancestors=TreeDat$ancestors,descendents=TreeDat$descendents,edge.length=TreeDat$edge.length)
  Res <- data.frame(rbind(Res,TreeDat)) 
}
Res$code <- paste(Res$tree, Res$ancestors,sep=".")
Res$code2 <- paste(Res$tree, Res$descendents,sep=".")
Res <- data.frame(Res, probAnc1=Res2$X1[match(Res$code,Res2$code)], probAnc2=Res2$X2[match(Res$code,Res2$code)], probAnc3=Res2$X3[match(Res$code,Res2$code)], probDes1=Res2$X1[match(Res$code2,Res2$code)], probDes2=Res2$X2[match(Res$code2,Res2$code)], probDes3=Res2$X3[match(Res$code2,Res2$code)])
results <- list(Res, Res2)
return(results)
}


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### DATA MANIPULATION ###

## polyandry ##
poly <- read.csv("data/polyandry.csv")
# collapse to one entry species
polyData <- summaryBy(perEPPbr + EPPbrN + perEGPbr + EGPbrN + groupFormation + Riehl ~ animal, data=poly, FUN=mean, keep.names=TRUE, na.rm=TRUE)
polyData$newG <- polyData$groupFormation
polyData$groupFormation <- "nonCoop"
polyData$groupFormation[which(polyData$newG == 1)] <- "Family"
polyData$groupFormation[which(polyData$newG == 3)] <- "nonFamily"
polyData$groupFormation <- factor(polyData$groupFormation)
# make EPPbr and EGPbr multinomial
polyData$nBrWithEPP <- round((polyData$perEPPbr/100) * polyData$EPPbrN)
polyData$nBrWithoutEPP <- round(polyData$EPPbrN - polyData$nBrWithEPP)
polyData$nBrWithEGP <- round((polyData$perEGPbr/100) * polyData$EGPbrN)
polyData$nBrWithoutEGP <- round(polyData$EGPbrN - polyData$nBrWithEGP)
# create family and nonfamily columns for recons
polyData$FamilyGroup <- ifelse(polyData$groupFormation == "Family", "Family", "nonCoop")
polyData$nonFamilyGroup <- ifelse(polyData$groupFormation == "nonFamily", "nonFamily", "nonCoop")

## reproductive division of labour ##
skew <- read.csv("data/skew.csv")
skew$gF2 <- ifelse(skew$groupFormation == "nonFamily", 1, 3)		# do this for plotting

## group size ##
groupSize <- read.csv("data/groupSize.csv")
groupSize$gF2 <- ifelse(groupSize$groupFormation == "nonFamily", 1, 3)		# do this for plotting

## specialization effect sizes ##
effects <- read.csv("data/effectSizes.csv")
# collapse to one effect size per species
uniqueEffects <- summaryBy(zRme + RmeN + varRme + zRmi + RmiN + varRmi + zRsx + RsxN + varRsx + groupFormation ~ animal, data=effects, keep.names=TRUE, na.rm=TRUE)
uniqueEffects$groupFormation <- ifelse(uniqueEffects$groupFormation == 2, "nonFamily", "Family")
uniqueEffects$groupFormation <- as.factor(uniqueEffects$groupFormation)
uniqueEffects$gF2 <- ifelse(uniqueEffects$groupFormation == "nonFamily", 1, 3)		# do this for plotting

## phylogenetic trees ##
trees <- read.nexus("data/BirdZilla1300.nex")
is.ultrametric(trees[sample(1:1300, 10)])


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #
#												MAIN ANALYSES											 #
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 1.1: ANCESTRAL STATE ESTIMATION ###

## (a). Estimating Breeding Systems using MCMCglmm ##

# these models were run on uppMax - took c. 30 days
IJ <- 10*(1/3)*(diag(2) + matrix(1, 2, 2))
pr.recon <- list(B = list(mu=rep(0,2), V=diag(2)*(1+pi^2/3)), R = list(V=IJ, fix=1), G = list(G1=list(V=IJ, nu=2, alpha.mu=rep(0,2), alpha.V=25^2*diag(2))))

INtree <- inverseA(trees[[1]], nodes="ALL")
reconModel.start <- MCMCglmm(groupFormation ~ trait-1, random=~us(trait):animal, rcov=~us(trait):units, data = polyData, ginverse=list(animal=INtree$Ainv), family ="categorical", nodes="ALL",  prior=pr.recon, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=0, verbose=FALSE)

reconModel <- reconModel.start
for(i in 1:1300){
  INtree <- inverseA(trees[[i]], nodes="ALL")
  start <- list(Liab=reconModel$Liab[1,], R=list(R1=matrix(ncol=2,nrow=2,reconModel$VCV[1,5:8])), G=list(G1=matrix(ncol=2,nrow=2,reconModel$VCV[1,1:4])))
  reconModel <- MCMCglmm(groupFormation ~ trait-1, random=~us(trait):animal, rcov=~us(trait):units, data = polyData, ginverse=list(animal=INtree$Ainv), family ="categorical", nodes="ALL", prior=pr.recon, pr=TRUE, pl=TRUE, slice=TRUE, nitt=100000, thin=1, burnin=99999, start=start, verbose=FALSE)
  if(i > 300){
    reconModel.start$VCV[i-300,] <- reconModel$VCV[1,]
    reconModel.start$Sol[i-300,] <- reconModel$Sol[1,]
    reconModel.start$Liab[i-300,] <- reconModel$Liab[1,]
  }
  print(i)
}
reconModel_A <- reconModel.start	
save(reconModel_A, file="reconModel_A")
load("results/reconModel_A")

# chain convergence
hist(reconModel_A$Liab)
plot(reconModel_A$VCV)     				# look well mixed
plot(reconModel_A$Sol[,1:3])   			# intercept estimates for nonCoop and nonFamily well mixed
autocorr(reconModel_A$VCV)   			# correlation between successive samples < 0.1 for all components
autocorr(reconModel_A$Sol[,1:5])  		# correlation between successive samples < 0.1 for all components
reconModelSols <- mcmc.list(list(reconModel_A$Sol, reconModel_B$Sol, reconModel_C$Sol))
plot(reconModelSols)
gelman.diag(reconModelSols)    			# upper CI = 1.00 for intercept suggesting convergence
heidel.diag(reconModel_A$VCV)   			# esimated variances passed halfwidth
heidel.diag(reconModel_A$Sol[,1:3])   	# intercepts passed halfwidth

# model parameters
summary(reconModel_A)
posterior.mode(reconModel_A$Sol[,1:2])		# nonCoop = 3.58, nonFamily = -2.44
HPDinterval(reconModel_A$Sol[,1:2])			# nonCoop = -0.61 to 7.47, nonFamily = -4.73 to 2.23

# find ancestors #
multM <- AncMultTree(responselevels=3, model=reconModel_A, trees=trees, starttree=300, endtree = 1300)
save(multM, file="results/multM")
load("results/multM")

# create Transition dataframes #
# split multM[[1]] up into a list: 1000 dataframes = trees 301 to 1300
transAnals <- split(multM[[1]], f = multM[[1]]$tree)
head(transAnals[[1]])

# assign ancestral transition categories (note - this was also done with 0.33 and 0.9 values)
transAnals2 <- lapply(transAnals, function(x) {x$ancms2 <- ifelse(x$probAnc1 > 0.67, "Family", "NA"); return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$ancms2 <- ifelse(x$probAnc2 > 0.67, "nonCoop", x$ancms2); return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$ancms2 <- ifelse(x$probAnc3 > 0.67, "nonFamily", x$ancms2); return(x)})

# assign descendent transition categories
transAnals2 <- lapply(transAnals2, function(x) {x$desms2 <- ifelse(x$probDes1 > 0.67, "Family", "NA"); return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$desms2 <- ifelse(x$probDes2 > 0.67, "nonCoop", x$desms2); return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$desms2 <- ifelse(x$probDes3 > 0.67, "nonFamily", x$desms2); return(x)})

# assign descendents that are species, their actual breeding systems (i.e. not the predicted ones)
transAnals2 <- lapply(transAnals2, function(x) {x$realValues <- polyData$groupFormation[match(x$descendents, polyData$animal)]; return(x)})

# work out the descendents of each node for desms2
transAnals2 <- lapply(transAnals2, function(x) {
obs <- data.frame(table(x$desms2, x$ancestors))
nodes <- unique(obs$Var2)
nodePatterns <- data.frame(node = unique(obs$Var2), patternA = NA, patternB = NA, patternC = NA, DES = NA)
for(i in 1:length(nodes)){nodePatterns$patternA[i] <- obs$Freq[which(obs$Var2 == nodes[i])][1]
	nodePatterns$patternB[i] <- obs$Freq[which(obs$Var2 == nodes[i])][3]
	nodePatterns$patternC[i] <- obs$Freq[which(obs$Var2 == nodes[i])][4]}	# have to skip [2] here which is "NA", use 1, 2, 3 for 0.33 cut-off
# the above can give 0,0,0 in the node pattern
nodePatterns$DES[nodePatterns$patternA == 2 & nodePatterns$patternB == 0 & nodePatterns$patternC == 0]	<- "onlyFamily"
nodePatterns$DES[nodePatterns$patternA == 0 & nodePatterns$patternB == 2 & nodePatterns$patternC == 0]	<- "onlynonCoop"
nodePatterns$DES[nodePatterns$patternA == 0 & nodePatterns$patternB == 0 & nodePatterns$patternC == 2]	<- "onlynonFamily"
nodePatterns$DES[nodePatterns$patternA == 1 & nodePatterns$patternB == 1 & nodePatterns$patternC == 0]	<- "nonCoop.Family"
nodePatterns$DES[nodePatterns$patternA == 0 & nodePatterns$patternB == 1 & nodePatterns$patternC == 1]	<- "nonCoop.nonFamily"
nodePatterns$DES[nodePatterns$patternA == 1 & nodePatterns$patternB == 0 & nodePatterns$patternC == 1]	<- "nonFamily.Family"
length(which(is.na(nodePatterns$DES)))		# 841 NAs (i.e. descendants unknown for these!)
# add these to transAnals2
x$DES2 <- nodePatterns$DES[match(x$ancestors, nodePatterns$node)]
x$CAT5 <- paste(x$ancms2, x$DES2, sep=".") ; return(x)})

# create CAT6: a factor with 6 levels in each object in transAnals2
# note - some of the trees will be missing certain types of transitions
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6 <- "other"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6[which(x$CAT5 == "nonCoop.onlynonCoop")] <- "ncTonc"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6[which(x$CAT5 == "nonCoop.nonCoop.Family")] <- "famOrigin"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6[which(x$CAT5 == "nonCoop.onlyFamily")] <- "famOrigin"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6[which(x$CAT5 == "Family.onlyFamily")] <- "famTofam"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6[which(x$CAT5 == "nonCoop.nonCoop.nonFamily")] <- "nonFamOrigin"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6[which(x$CAT5 == "nonCoop.onlynonFamily")] <- "nonFamOrigin"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6[which(x$CAT5 == "nonFamily.onlynonFamily")] <- "nonFamTononFam"; return(x)})

# turn each CAT6 into a factor
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6 <- factor(x$CAT6); return(x)})

# assign all six levels to each CAT6 (some trees may be missing certain transitions)
CAT6lengths <- lapply(transAnals2, function(x) {length(names(table(x$CAT6)))})
CAT6lengths2 <- unlist(CAT6lengths)
CAT6levels <- c("other", "ncTonc", "famOrigin", "famTofam", "nonFamOrigin", "nonFamTononFam")
transAnals2[[1]]$CAT6 <- factor(transAnals2[[1]]$CAT6, levels = CAT6levels)
transAnals2 <- lapply(transAnals2, function(x) {x$CAT6 <- factor(x$CAT6, levels = CAT6levels); return(x)})

# work out the average number of transitions across trees
# want four possible transition types (use CAT5 for this)
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5 <- "other"; return(x)})
# want "nonCoop.onlyFamily" and "nonCoop.nonCoop.Family" to be "nonCoopToFam"
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "nonCoop.onlyFamily")] <- "nonCoopToFam"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "nonCoop.nonCoop.Family")] <- "nonCoopToFam"; return(x)})
# want "nonCoop.onlynonFamily" and "nonCoop.nonCoop.nonFamily" to be "nonCoopToNonFam"
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "nonCoop.onlynonFamily")] <- "nonCoopToNonFam"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "nonCoop.nonCoop.nonFamily")] <- "nonCoopToNonFam"; return(x)})
# want "nonFamily.nonCoop.Family", "nonFamily.onlyFamily" and "nonFamily.nonFamily.Family" to be "nonFamToFam"
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "nonFamily.nonCoop.Family")] <- "nonFamToFam"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "nonFamily.onlyFamily")] <- "nonFamToFam"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "nonFamily.nonFamily.Family")] <- "nonFamToFam"; return(x)})
# want "Family.nonCoop.nonFamily", "Family.onlynonFamily" and "Family.nonFamily.Family" to be "famToNonFam"
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "Family.nonCoop.nonFamily")] <- "famToNonFam"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "Family.onlynonFamily")] <- "famToNonFam"; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$figCAT5[which(x$CAT5 == "Family.nonFamily.Family")] <- "famToNonFam"; return(x)})

# want family and nonFamily origins from each other and from nonCoop
results3 <- lapply(transAnals2, function(x) {table(x$figCAT5)})
transitions3 <- data.frame(nonCoopToFam = 1:1000, nonFamToFam = 1:1000, nonCoopToNonFam = 1:1000, famToNonFam = 1:1000)
for(i in 1:1000) {transitions3[i,] <- results3[[i]][match(names(transitions3), names(results3[[i]]))]}
# non cooperative to family
posterior.mode(mcmc(transitions3$nonCoopToFam), na.rm=TRUE)/2		# 132
HPDinterval(mcmc(transitions3$nonCoopToFam))/2						# 58 to 254
# non-family to family
posterior.mode(mcmc(transitions3$nonFamToFam), na.rm=TRUE)/2			# 1
HPDinterval(mcmc(transitions3$nonFamToFam))/2						# 1 to 16
# non cooperative to non-family
posterior.mode(mcmc(transitions3$nonCoopToNonFam), na.rm=TRUE)/2		# 14
HPDinterval(mcmc(transitions3$nonCoopToNonFam))/2					# 1 to 228
# family to non-family
posterior.mode(mcmc(transitions3$famToNonFam), na.rm=TRUE)/2			# 1
HPDinterval(mcmc(transitions3$famToNonFam))/2						# 1 to 12


## transition rate models using SCM with q = empirical ##
x <- polyData$groupFormation
names(x) <- polyData$animal
scmTrees <- trees[c(301:1300)]
# equal rates
scms <- make.simmap(scmTrees, x, model="ER", nsim=10)
describe.simmap(scms, plot=FALSE)		# 92 family group origins and 29 nonFamily origins
save(scms, file="results/SCMs")
load("results/SCMs")
# all rates different
scmARD <- make.simmap(scmTrees, x, model="ARD", nsim=10)
describe.simmap(scmARD, plot=FALSE)		# 71 family group origins and 29 nonFamily origins
save(scmARD, file="results/scmARD")



## (b). Estimating Polyandry Rates using MCMCglmm ##

# make the ancestors animal
transAnals2 <- lapply(transAnals2, function(x) {x$animal <- x$ancestors; return(x)})

# add perEPPbr
transAnals2 <- lapply(transAnals2, function(x) {x$nBrWithEPP <- polyData$nBrWithEPP[match(x$descendents, polyData$animal)]; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$nBrWithoutEPP <- polyData$nBrWithoutEPP[match(x$descendents, polyData$animal)]; return(x)})

# add perEGPbr
transAnals2 <- lapply(transAnals2, function(x) {x$nBrWithEGP <- polyData$nBrWithEGP[match(x$descendents, polyData$animal)]; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {x$nBrWithoutEGP <- polyData$nBrWithoutEGP[match(x$descendents, polyData$animal)]; return(x)})

# delete node 1 from transitions data frames
# note that this will mess stuff up if already deleted!
transAnals2 <- lapply(transAnals2, function(x) {x <- x[-which(x$animal == "Node1"),]; return(x)})
transAnals2 <- lapply(transAnals2, function(x) {rownames(x) <- NULL; return(x)})

save(transAnals2, file="results/transAnals2")
load("results/transAnals2")

# the first 300 trees were discarded in reconModel, therefore the first object in transAnals2 is tree 301 and the last is tree 1300
trees2 <- trees[301:1300]

# build the perEPPbr model with CAT6
pr.poly <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
INtree <- inverseA(trees2[[1]], nodes="ALL")

# this needs to be 700 iterations
polyModel.start <- MCMCglmm(cbind(nBrWithEPP, nBrWithoutEPP) ~ CAT6-1, random=~animal, data = transAnals2[[1]], ginverse=list(animal=INtree$Ainv), family ="multinomial2", nodes="ALL", prior=pr.poly, pr=TRUE, pl=TRUE, slice=TRUE, nitt=700, thin=1, burnin=0, verbos=FALSE, singular.ok=TRUE)

polyModel <- polyModel.start
for(i in 1:1000){
  INtree <- inverseA(trees2[[i]], nodes="ALL")
  start <- list(Liab=polyModel$Liab[1,], R=polyModel$VCV[1,2], G=list(G1=polyModel$VCV[1,1]))
  polyModel <- MCMCglmm(cbind(nBrWithEPP, nBrWithoutEPP) ~ CAT6-1, random=~animal, data = transAnals2[[i]], ginverse=list(animal=INtree$Ainv), family ="multinomial2", nodes="ALL", prior=pr.poly, pr=TRUE, pl=TRUE, slice=TRUE, nitt=50000, thin=1, burnin=49999, start=start, verbose=FALSE, singular.ok=TRUE)
  if(i > 300){
    polyModel.start$VCV[i-300,] <- polyModel$VCV[1,]
    polyModel.start$Sol[i-300,] <- polyModel$Sol[1,]
    polyModel.start$Liab[i-300,] <- polyModel$Liab[1,]
  }
  print(i)
}
polyModel3A <- polyModel.start
polyModel3B <- polyModel.start
polyModel3C <- polyModel.start
save(polyModel3A, file="results/polyModel3A")
save(polyModel3B, file="results/polyModel3B")
save(polyModel3C, file="results/polyModel3C")
load("results/polyModel3A")

# chain convergence
hist(polyModel3A$Liab)
plot(polyModel3A$VCV)     			# animal and units close to 0
plot(polyModel3A$Sol[,1:6])   		# intercept estimate well mixed
autocorr(polyModel3A$VCV)   			# correlation between successive samples < 0.1 for all components
autocorr(polyModel3A$Sol[,1:6])  	# correlation between successive samples < 0.1 for all components
polyModel3Sols <- mcmc.list(list(polyModel3A$Sol, polyModel3B$Sol, polyModel3C$Sol))
plot(polyModel3Sols)
gelman.diag(polyModel3Sols)     		# upper CI = 1.00 for intercept suggesting convergence
heidel.diag(polyModel3A$VCV)   		# animal and units passed halfwidth
heidel.diag(polyModel3A$Sol[,1:6])

summary(polyModel3A)
plot(polyModel3A$Sol[,1:6])		# 4/6 chains well mixed
posterior.mode(inv.logit(polyModel3A$Sol[,1:6]))
HPDinterval(inv.logit(polyModel3A$Sol[,1:6]))

# test for differences in ancestral polyandry rates
# origin of nonCoop (0.20) vs. family (0.17)
table(polyModel3A$Sol[, 2] > polyModel3A$Sol[, 3]) / length(polyModel3A$Sol[, 1])    # FALSE = 0.40, TRUE = 0.60
# origin of nonCoop (0.20) vs nonFamily (0.41)
table(polyModel3A$Sol[, 2] > polyModel3A$Sol[, 5]) / length(polyModel3A$Sol[, 1])   # FALSE = 0.95, TRUE = 0.05
# origin of family vs. origin of nonfamily
table(polyModel3A$Sol[,3] > polyModel3A$Sol[,5]) / length(polyModel3A$Sol[,1])    # FALSE = 0.95, TRUE = 0.05


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 1.2: GROUP FORMATION AND REPRODUCTIVE DIVISION OF LABOUR ###

# trim the trees
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, skew$animal)))]
skewTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
skewTrees <- lapply(skewTrees, makeNodeLabel, method = "number")
skew$animal[which((skew$animal %in% skewTrees[[1]]$tip.label) == FALSE)]
skewTrees[[1]]$tip.label[which((skewTrees[[1]]$tip.label %in% skew$animal) == FALSE)]
# reduced to 48 species

# multi-response model on multiple trees
pr.skew <- list(R = list(V=diag(2), n = 1.002), G = list(G1=list(V=diag(2), n=1.002)))	# check if need the 1 here any more!
INtree <- inverseA(skewTrees[[1]], nodes="TIPS")
skewModel.start <- MCMCglmm(cbind(cbind(nGrMixedF, nGrMonogF), cbind(nGrMixedM, nGrMonogM)) ~ trait:groupFormation-1, random = ~idh(trait):animal, rcov = ~idh(trait):units, family = c("multinomial2", "multinomial2"), data=skew, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, prior=pr.skew, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
skewModel <- skewModel.start
for(i in 1:1300){
  INtree <- inverseA(skewTrees[[i]], nodes="TIPS")
  start <- list(Liab=skewModel$Liab[1,], R=diag(skewModel$VCV[1,3:4]), G=list(G1=diag(skewModel$VCV[1,1:2])))
  skewModel <- MCMCglmm(cbind(cbind(nGrMixedF, nGrMonogF), cbind(nGrMixedM, nGrMonogM)) ~ trait:groupFormation-1, random = ~idh(trait):animal, rcov = ~idh(trait):units, family = c("multinomial2", "multinomial2"), data=skew, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, prior=pr.skew, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    skewModel.start$VCV[i-300,] <- skewModel$VCV[1,]
    skewModel.start$Sol[i-300,] <- skewModel$Sol[1,]
    skewModel.start$Liab[i-300,] <- skewModel$Liab[1,]
  }
}
skewModelA <- skewModel.start
skewModelB <- skewModel.start
skewModelC <- skewModel.start
save(skewModelA, file="results/skewModelA")
save(skewModelB, file="results/skewModelB")
save(skewModelC, file="results/skewModelC")
load("results/skewModelA")

# chain convergence
hist(polyModel3A$Liab)
plot(polyModel3A$VCV)     			# animal and units close to 0
plot(polyModel3A$Sol[,1:6])   		# intercept estimate well mixed
autocorr(polyModel3A$VCV)   			# correlation between successive samples < 0.1 for all components
autocorr(polyModel3A$Sol[,1:6])  	# correlation between successive samples < 0.1 for all components
polyModel3Sols <- mcmc.list(list(polyModel3A$Sol, polyModel3B$Sol, polyModel3C$Sol))
plot(polyModel3Sols)
gelman.diag(polyModel3Sols)     		# upper CI = 1.00 for intercept suggesting convergence
heidel.diag(polyModel3A$VCV)   		# animal and units passed halfwidth
heidel.diag(polyModel3A$Sol[,1:6])

summary(skewModelA)
posterior.mode(inv.logit(skewModelA$Sol[,1:4]))
HPDinterval(inv.logit(skewModelA$Sol[,1:4]))
# Family Fmixed = 0.00 < 0.06 < 0.24
# Family Mmixed = 0.01 < 0.07 < 0.17
# nonFamily Fmixed = 0.73 < 0.98 < 0.99
# nonFamily Mmixed = 0.57 < 0.83 < 0.95

# parameters are means, not differences, therefore use table(a > b)/length(a) to test
table(skewModelA$Sol[,1] > skewModelA$Sol[,3]) / length(skewModelA$Sol[,3])   # FALSE = 1, TRUE = 0
table(skewModelA$Sol[,2] > skewModelA$Sol[,4]) / length(skewModelA$Sol[,4])   # FALSE = 1, TRUE = 0


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 1.3: GROUP FORMATION AND GROUP SIZE ###

# trim the trees
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, groupSize$animal)))]
gsTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
gsTrees <- lapply(gsTrees, makeNodeLabel, method = "number")
groupSize$animal[which((groupSize$animal %in% gsTrees[[1]]$tip.label) == FALSE)]
gsTrees[[1]]$tip.label[which((gsTrees[[1]]$tip.label %in% groupSize$animal) == FALSE)]
# reduced to 127 species

pr.gs <- list(R = list(V=1, n = 0.002), G = list(G1=list(V=1, n=0.002)))

# mean model on multiple trees
INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
gsModel.start <- MCMCglmm(log(mean) ~ groupFormation-1, random = ~animal, family="gaussian", data=groupSize, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, prior=pr.gs, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
gsModel <- gsModel.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=gsModel$Liab[1,], R=gsModel$VCV[1,2], G=list(G1=gsModel$VCV[1,1]))
  gsModel <- MCMCglmm(log(mean) ~ groupFormation-1, random = ~animal, family="gaussian", data=groupSize, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, prior=pr.gs, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    gsModel.start$VCV[i-300,] <- gsModel$VCV[1,]
    gsModel.start$Sol[i-300,] <- gsModel$Sol[1,]
    gsModel.start$Liab[i-300,] <- gsModel$Liab[1,]
  }
}
gsMeanModelA <- gsModel.start
gsMeanModelB <- gsModel.start
gsMeanModelC <- gsModel.start
save(gsMeanModelA, file="results/gsMeanModelA")
save(gsMeanModelB, file="results/gsMeanModelB")
save(gsMeanModelC, file="results/gsMeanModelC")
load("results/gsMeanModelA")
summary(gsMeanModelA)
posterior.mode(exp(gsMeanModelA$Sol))
HPDinterval(exp(gsMeanModelA$Sol))
# Family mean = 2.1 < 3.4 < 5.0
# nonFamily mean = 1.8 < 2.5 < 4.2
# parameters are means, not differences, therefore use table(a > b)/length(a) to test
table(gsMeanModelA$Sol[,1] > gsMeanModelA$Sol[,2]) / length(gsMeanModelA$Sol[,1])   # FALSE = 0.04, TRUE = 0.96


# max model on multiple trees
INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
gsModel.start <- MCMCglmm(max ~ groupFormation-1, random = ~animal, family="poisson", data=groupSize, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, prior=pr.gs, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
gsModel <- gsModel.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=gsModel$Liab[1,], R=gsModel$VCV[1,2], G=list(G1=gsModel$VCV[1,1]))
  gsModel <- MCMCglmm(max ~ groupFormation-1, random = ~animal, family="poisson", data=groupSize, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, prior=pr.gs, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    gsModel.start$VCV[i-300,] <- gsModel$VCV[1,]
    gsModel.start$Sol[i-300,] <- gsModel$Sol[1,]
    gsModel.start$Liab[i-300,] <- gsModel$Liab[1,]
  }
}
gsMaxModelA <- gsModel.start
gsMaxModelB <- gsModel.start
gsMaxModelC <- gsModel.start
save(gsMaxModelA, file="results/gsMaxModelA")
save(gsMaxModelB, file="results/gsMaxModelB")
save(gsMaxModelC, file="results/gsMaxModelC")
load("results/gsMaxModelA")
summary(gsMaxModelA)
posterior.mode(exp(gsMaxModelA$Sol))
HPDinterval(exp(gsMaxModelA$Sol))
# Family mean = 4.4 < 7.0 < 11.4
# nonFamily mean = 2.7 < 4.6 < 7.2
# parameters are means, not differences, therefore use table(a > b)/length(a) to test
table(gsMaxModelA$Sol[,1] > gsMaxModelA$Sol[,2]) / length(gsMaxModelA$Sol[,1])   # FALSE = 0.00, TRUE = 1.00


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 1.4: GROUP FORMATION AND REPRODUCTIVE SPECIALIZATION ###

# trim the trees
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, uniqueEffects$animal)))]
dolTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
dolTrees <- lapply(dolTrees, makeNodeLabel, method = "number")
uniqueEffects$animal[which((uniqueEffects$animal %in% dolTrees[[1]]$tip.label) == FALSE)]
dolTrees[[1]]$tip.label[which((dolTrees[[1]]$tip.label %in% uniqueEffects$animal) == FALSE)]
# reduced to 60 species


## (a). Publication Bias in Effect Sizes ##

# publication bias tests (done for each Zr seperately)

# zRfecundity (in metafor)
zRmiModela <- rma(zRmi, vi=varRmi, data=uniqueEffects)
summary(zRmiModela)    			# zRmi = 0.06 < 0.18 < 0.29; I^2 = 87%; between-study var = 0.10
funnel(trimfill(zRmiModela))		# 8 studies missing
zRmiModelb <- rma(zRmi ~ groupFormation, vi=varRmi, data=uniqueEffects)
summary(zRmiModelb)    			# I^2 = 73%; between-study var = 0.05
uniqueEffects$invSErmi <- 1/sqrt(uniqueEffects$varRmi)		# create precision measure
zRmiEggMod <- lm(scale(zRmi) ~ invSErmi, data=uniqueEffects[-which(is.na(uniqueEffects$zRmi)),])		# Egger's test
summary(zRmiEggMod)				# intercept = 0.31 (p = 0.2), slope = -0.05 (p = 0.1)
zRmiEggMod2 <- lm(scale(zRmi) ~ invSErmi + groupFormation, data=uniqueEffects[-which(is.na(uniqueEffects$zRmi)),])
summary(zRmiEggMod2)			# invSE = -0.03 (p = 0.3)
# including phylogeny
uniqueEffects$ID <- uniqueEffects$animal
birdCor <- vcv.phylo(dolTrees[[1]], cor=TRUE)
zRmiModel2 <- rma.mv(zRmi ~ 1, V=varRmi, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects)
summary(zRmiModel2)    			# zRmi = -0.35 < -0.01 < 0.34; between-study var = 0.003, phylo var = 0.15
zRmiModel3 <- rma.mv(zRmi ~ groupFormation-1, V=varRmi, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects)
summary(zRmiModel3)				# family: 0.21 < 0.32 < 0.43; nonFamily = -0.34 < -0.18 < -0.02

# zRcare (in metafor)
zRmeModela <- rma(zRme, vi=varRme, data=uniqueEffects)
summary(zRmeModela)    			# -0.44 < -0.36 < -0.28; I^2 = 34%; between-study var = 0.02
funnel(trimfill(zRmeModela))		# 17 studies missing
zRmeModelb <- rma(zRme ~ groupFormation, vi=varRme, data=uniqueEffects)
summary(zRmeModelb)    			# I^2 = 34%; between-study var = 0.02
uniqueEffects$invSErme <- 1/sqrt(uniqueEffects$varRme)		# create precision measures
zRmeEggMod <- lm(scale(zRme) ~ invSErme, data=uniqueEffects[-which(is.na(uniqueEffects$zRme)),])		# Egger's test
summary(zRmeEggMod)				# intercept = -0.7 (p = 0.01), slope = 0.16 (p = 0.01)
zRmeEggMod2 <- lm(scale(zRme) ~ invSErme + groupFormation, data=uniqueEffects[-which(is.na(uniqueEffects$zRme)),])
summary(zRmeEggMod2)			# invSE = 0.16 (p = 0.01)
# including phylogeny
zRmeModel2 <- rma.mv(zRme ~ 1, V=varRme, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects)
summary(zRmeModel2)   			# zRmi = -0.46 < -0.33 < -0.20; between-study var = 0.02, phylo var = 0.01
zRmeModel3 <- rma.mv(zRme ~ groupFormation-1, V=varRme, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects)
summary(zRmeModel3)   			# family: -0.48 < -0.34 < -0.20; nonFamily = -0.51 < -0.31 < -0.10

# zRsurvival (in metafor)
zRsxModela <- rma(zRsx, vi=varRsx, data=uniqueEffects)
summary(zRsxModela)    			# 0.02 < 0.08 < 0.14; I^2 = 65%; between-study var = 0.01
funnel(trimfill(zRsxModela))		# 0 studies missing
zRsxModelb <- rma(zRsx ~ groupFormation, vi=varRsx, data=uniqueEffects)
summary(zRsxModelb)   			# I^2 = 59%; between-study var = 0.01
uniqueEffects$invSErsx <- 1/sqrt(uniqueEffects$varRsx)		# create precision measures
zRsxEggMod <- lm(scale(zRsx) ~ invSErsx, data=uniqueEffects[-which(is.na(uniqueEffects$zRsx)),])		# Egger's test
summary(zRsxEggMod)				# intercept = 0.29 (p = 0.4), slope = -0.02 (p = 0.3)
zRsxEggMod2 <- lm(scale(zRsx) ~ invSErsx + groupFormation, data=uniqueEffects[-which(is.na(uniqueEffects$zRsx)),])
summary(zRsxEggMod2)			# invSE = -0.03 (p = 0.2)
# including phylogeny
zRsxModel2 <- rma.mv(zRsx ~ 1, V=varRsx, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects)
summary(zRsxModel2)    			# zRmi = 0.001 < 0.07 < 0.14; between-study var = 0.01, phylo var = 0.001
zRsxModel3 <- rma.mv(zRsx ~ groupFormation-1, V=varRsx, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects)
summary(zRsxModel3) 			# family: -0.03 < 0.07 < 0.16; nonFamily = -0.08 < 0.07 < 0.23


## (b). Testing for Changes in Fecundity, Care and Survival with Group Size ##

# reformat the data to be able to use a single model
meltedzR <- melt(uniqueEffects[,c("animal", "groupFormation", "zRme", "zRmi", "zRsx")], id.vars=c("animal", "groupFormation"))
meltedVar <- melt(uniqueEffects[,c("animal", "groupFormation", "varRme", "varRmi", "varRsx")], id.vars=c("animal", "groupFormation"))
# note that animals appear in the same order in each melted data frame
meltedEffects <- data.frame(animal = meltedzR$animal, groupFormation = meltedzR$groupFormation, effectType = meltedzR$variable, effectSize = meltedzR$value, var = meltedVar$value)
# get rid of NAs
meltedEffects <- meltedEffects[-which(is.na(meltedEffects$effectSize)),]
tapply(meltedEffects$effectSize, list(meltedEffects$effectType, meltedEffects$groupFormation), FUN=mean, na.rm=T)

# random meta-analytic model on multiple trees
pr.zR <- list(R = list(V=diag(3), n = 2.002), G = list(G1=list(V=diag(3), n=2.002)))		# don't think you need the 1 here any more!
MEV <- meltedEffects$var
INtree <- inverseA(dolTrees[[1]], nodes="TIPS")
zRModel.start <- MCMCglmm(effectSize ~ effectType:groupFormation-1, random = ~idh(effectType):animal, rcov = ~idh(effectType):units, data=meltedEffects, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=pr.zR, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
zRModel <- zRModel.start
for(i in 1:1300){
  INtree <- inverseA(dolTrees[[i]], nodes="TIPS")
    start <- list(Liab=zRModel$Liab[1,], R=diag(zRModel$VCV[1,5:7]), G=list(G1=diag(zRModel$VCV[1,1:3])))
  zRModel <- MCMCglmm(effectSize ~ effectType:groupFormation-1, random = ~idh(effectType):animal, rcov = ~idh(effectType):units, data=meltedEffects, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=pr.zR, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    zRModel.start$VCV[i-300,] <- zRModel$VCV[1,]
    zRModel.start$Sol[i-300,] <- zRModel$Sol[1,]
    zRModel.start$Liab[i-300,] <- zRModel$Liab[1,]
  }
}
zRModelA <- zRModel.start
zRModelB <- zRModel.start
zRModelC <- zRModel.start
save(zRModelA, file="results/zRModelA")
save(zRModelB, file="results/zRModelB")
save(zRModelC, file="results/zRModelC")
load("results/zRModelC")

summary(zRModelC)	# gives 6 parameters
posterior.mode(zRModelC$Sol)
HPDinterval(zRModelC$Sol)
# parameters are means, not differences, therefore use table(a > b)/length(a) to test
# zRmi family (0.06 < 0.33 < 0.50) vs. nonFamily (-0.41 < -0.14 < 0.04)
table(zRModelC$Sol[,2] > zRModelC$Sol[,5]) / length(zRModelC$Sol[,2])   # FALSE = 0.0, TRUE = 1.0
# zRme family (-0.53 < -0.34 < -0.15) vs. nonFamily (-0.57 < -0.32 < -0.11)
table(zRModelC$Sol[,1] > zRModelC$Sol[,4]) / length(zRModelC$Sol[,2])   # FALSE = 0.60, TRUE = 40
# zRsx family (-0.09 < 0.06 < 0.21) vs. nonFamily (-0.12 < 0.10 < 0.29)
table(zRModelC$Sol[,3] > zRModelC$Sol[,6]) / length(zRModelC$Sol[,2])   # FALSE = 57, TRUE = 0.43

# note that these parameters agree with phylogenetic metafor models
summary(rma.mv(zRme ~ groupFormation, V=varRme, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects))			# fam = -0.34, diff = -0.17 < 0.04 < 0.24
summary(rma.mv(zRmi ~ groupFormation, V=varRmi, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects))			# fam = 0.31, diff = -0.69 < -0.50 < -0.31
summary(rma.mv(zRsx ~ groupFormation, V=varRsx, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueEffects))			# fam = 0.07, diff = -0.17 < 0.01 < 0.17



## (c). The Relationship Between Fecundity and Care in Family and Non-Family Groups ##

# trim data (have to do this - can't handle missing data in MEV term)
zRmezRmiData <- uniqueEffects[-c(which(is.na(uniqueEffects$zRme)), which(is.na(uniqueEffects$zRmi))),]	# 28 species: 21 fam, 7 nonfam
# run the model on multiple trees
pr.zRmezRmi <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
MEV <- zRmezRmiData$varRmi
INtree <- inverseA(dolTrees[[1]], nodes="TIPS")
zRmezRmiModel.start <- MCMCglmm(zRmi ~ zRme*groupFormation + log(RmeN), random= ~animal, family="gaussian", data=zRmezRmiData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=pr.zRmezRmi, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
zRmezRmiModel <- zRmezRmiModel.start
for(i in 1:1300){
  INtree <- inverseA(dolTrees[[i]], nodes="TIPS")
  start <- list(Liab=zRmezRmiModel$Liab[1,], R=zRmezRmiModel$VCV[1,3], G=list(G1=zRmezRmiModel$VCV[1,1]))
  zRmezRmiModel <- MCMCglmm(zRmi ~ zRme*groupFormation + log(RmeN), random= ~animal, family="gaussian", data=zRmezRmiData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=pr.zRmezRmi, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    zRmezRmiModel.start$VCV[i-300,] <- zRmezRmiModel$VCV[1,]
    zRmezRmiModel.start$Sol[i-300,] <- zRmezRmiModel$Sol[1,]
    zRmezRmiModel.start$Liab[i-300,] <- zRmezRmiModel$Liab[1,]
  }
}
zRmezRmiModelA <- zRmezRmiModel.start
zRmezRmiModelB <- zRmezRmiModel.start
zRmezRmiModelC <- zRmezRmiModel.start
save(zRmezRmiModelA, file="results/zRmezRmiModelA")
save(zRmezRmiModelB, file="results/zRmezRmiModelB")
save(zRmezRmiModelC, file="results/zRmezRmiModelC")
load("results/zRmezRmiModelA")

summary(zRmezRmiModelA)
posterior.mode(zRmezRmiModelA$Sol[,1:2])						# family: intercept = -0.16 slope = -1.02
HPDinterval(zRmezRmiModelA$Sol[,1])								# family intercept = -0.64 to 0.21
HPDinterval(zRmezRmiModelA$Sol[,2])								# family slope = -1.42 to -0.62
posterior.mode(zRmezRmiModelA$Sol[,1]+zRmezRmiModelA$Sol[,3])	# nonFamily: intercept = -0.25
posterior.mode(zRmezRmiModelA$Sol[,2]+zRmezRmiModelA$Sol[,5])	# nonFamily: slope = 0.11
HPDinterval(zRmezRmiModelA$Sol[,1:5])							# slope difference = 0.58 to 1.53
HPDinterval(zRmezRmiModelA$Sol[,1]+zRmezRmiModelA$Sol[,3])		# nonFamily intercept = -0.36 to 0.41
HPDinterval(zRmezRmiModelA$Sol[,2]+zRmezRmiModelA$Sol[,5])		# nonFamily slope = -0.36 to 0.41
table(zRmezRmiModelA$Sol[,2] > (zRmezRmiModelA$Sol[,2] + zRmezRmiModelA$Sol[,5])) / length(zRmezRmiModelA$Sol[,2])		# FALSE = 1
posterior.mode(zRmezRmiModelA$Sol[,4])							# sample size = 0.04
HPDinterval(zRmezRmiModelA$Sol[,4])								# sample size = -0.05 to 0.13


## (d). The Relationship Between Fecundity and Survival in Family Groups ##

# trim data (have to do this - can't handle missing data in MEV term)
zRmizRsxData <- uniqueEffects[-c(which(is.na(uniqueEffects$zRmi)), which(is.na(uniqueEffects$zRsx))),]		# 14 species: 11 fam, 3 nonfam
# delete non-family cooperative breeders
zRmizRsxData <- zRmizRsxData[-which(zRmizRsxData$groupFormation =="nonFamily"),]
# run the model on multiple trees
pr.zRmizRsx <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
MEV <- zRmizRsxData$varRsx
INtree <- inverseA(dolTrees[[1]], nodes="TIPS")
zRmizRsxModel.start <- MCMCglmm(zRmi ~  zRsx + log(RsxN), random= ~animal, family ="gaussian", data=zRmizRsxData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=pr.zRmizRsx, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
zRmizRsxModel <- zRmizRsxModel.start
for(i in 1:1300){
  INtree <- inverseA(dolTrees[[i]], nodes="TIPS")
  start <- list(Liab=zRmizRsxModel$Liab[1,], R=zRmizRsxModel$VCV[1,3], G=list(G1=zRmizRsxModel$VCV[1,1]))
  zRmizRsxModel <- MCMCglmm(zRmi ~ zRsx + log(RsxN), random= ~animal, family ="gaussian", data=zRmizRsxData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=pr.zRmizRsx, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    zRmizRsxModel.start$VCV[i-300,] <- zRmizRsxModel$VCV[1,]
    zRmizRsxModel.start$Sol[i-300,] <- zRmizRsxModel$Sol[1,]
    zRmizRsxModel.start$Liab[i-300,] <- zRmizRsxModel$Liab[1,]
  }
}
zRmizRsxModelA <- zRmizRsxModel.start
zRmizRsxModelB <- zRmizRsxModel.start
zRmizRsxModelC <- zRmizRsxModel.start
save(zRmizRsxModelA, file="results/zRmizRsxModelA")
save(zRmizRsxModelB, file="results/zRmizRsxModelB")
save(zRmizRsxModelC, file="results/zRmizRsxModelC")
load("results/zRmizRsxModelA")

summary(zRmizRsxModelA)
posterior.mode(zRmizRsxModelA$Sol)		# intercept = 0.71 slope = 1.15, sample size = -0.10
HPDinterval(zRmizRsxModelA$Sol)			# intercept = -0.21 to 1.20, slope = -0.11 to 2.17, sample size = -0.23 to 0.02



## (e). The Relationship Between Survival and Care in Family Groups ##

# trim data (have to do this - can't handle missing data in MEV term)
zRmezRsxData <- uniqueEffects[-c(which(is.na(uniqueEffects$zRme)), which(is.na(uniqueEffects$zRsx))),]		# 19 species: 17 fam, 2 nonfam
# delete non-family cooperative breeders
zRmezRsxData <- zRmezRsxData[-which(zRmezRsxData$groupFormation =="nonFamily"),]
# run the model on multiple trees
pr.zRmezRsx <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
MEV <- zRmezRsxData$varRsx
INtree <- inverseA(dolTrees[[1]], nodes="TIPS")
zRmezRsxModel.start <- MCMCglmm(zRsx ~ zRme + log(RmeN), random= ~animal, family="gaussian", data=zRmezRsxData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=pr.zRmezRsx, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
zRmezRsxModel <- zRmezRsxModel.start
for(i in 1:1300){
  INtree <- inverseA(dolTrees[[i]], nodes="TIPS")
  start <- list(Liab=zRmezRsxModel$Liab[1,], R=zRmezRsxModel$VCV[1,3], G=list(G1=zRmezRsxModel$VCV[1,1]))
  zRmezRsxModel <- MCMCglmm(zRsx ~ zRme + log(RmeN), random= ~animal, family="gaussian", data=zRmezRsxData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=pr.zRmezRsx, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    zRmezRsxModel.start$VCV[i-300,] <- zRmezRsxModel$VCV[1,]
    zRmezRsxModel.start$Sol[i-300,] <- zRmezRsxModel$Sol[1,]
    zRmezRsxModel.start$Liab[i-300,] <- zRmezRsxModel$Liab[1,]
  }
}
zRmezRsxModelA <- zRmezRsxModel.start
zRmezRsxModelB <- zRmezRsxModel.start
zRmezRsxModelC <- zRmezRsxModel.start
save(zRmezRsxModelA, file="results/zRmezRsxModelA")
save(zRmezRsxModelB, file="results/zRmezRsxModelB")
save(zRmezRsxModelC, file="results/zRmezRsxModelC")
load("results/zRmezRsxModelA")

summary(zRmezRsxModelA)
posterior.mode(zRmezRsxModelA$Sol)		# intercept = -0.16 slope = -0.37
HPDinterval(zRmezRsxModelA$Sol)			# slope = -0.72 to -0.05


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #
#											SENSITIVITY ANALYSES											 #
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 2.1: BREEDING SYSTEM CLASSIFICATION ###

# Is Riehl's category 3 nested in any of the other states and how are these related to each other?
# using SCMs
x <- polyData$Riehl
x[which(is.na(x))] <- 5
names(x) <- polyData$animal
scmTrees <- trees[301:1300]
RiehlscmARD <- make.simmap(scmTrees, x, model="ARD", nsim=1)
RiehlscmER <- make.simmap(scmTrees, x, model="ER", nsim=1)
save(RiehlscmARD, file="results/RiehlscmARD")
save(RiehlscmER, file="results/RiehlscmER")
describe.simmap(RiehlscmARD, plot=FALSE)	# 71 family group origins and 29 nonFamily origins
describe.simmap(RiehlscmER, plot=FALSE)		# 89 family group origins and 27 nonFamily origins
# details on stochastic character mapping - see Huelsenbeck + (2003), Bollback (2006), Yang (2006)


# How do Riehl's categories differ in terms of reprodcutive division of labour, group size and specialization?

## Reproductive division of labour ##
# add Riehl classification to skew data
skew$Riehl <- as.factor(polyData$Riehl[match(skew$animal, polyData$animal)])
table(skew$Riehl)		# cat 1 = 17, cat 2 = 9, cat 3 = 14, cat 4 = 7

# multi-response model on multiple trees
pr.skew <- list(R = list(V=diag(2), n = 1.002), G = list(G1=list(V=diag(2), n=1.002)))
INtree <- inverseA(skewTrees[[1]], nodes="TIPS")
skewCheckMod.start <- MCMCglmm(cbind(cbind(nGrMixedF, nGrMonogF), cbind(nGrMixedM, nGrMonogM)) ~ trait:Riehl-1, random = ~idh(trait):animal, rcov = ~idh(trait):units, family = c("multinomial2", "multinomial2"), data=skew, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, prior=pr.skew, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
skewCheckMod <- skewCheckMod.start
for(i in 1:1300){
  INtree <- inverseA(skewTrees[[i]], nodes="TIPS")
  start <- list(Liab=skewCheckMod$Liab[1,], R=diag(skewCheckMod$VCV[1,3:4]), G=list(G1=diag(skewCheckMod$VCV[1,1:2])))
  skewCheckMod <- MCMCglmm(cbind(cbind(nGrMixedF, nGrMonogF), cbind(nGrMixedM, nGrMonogM)) ~ trait:Riehl-1, random = ~idh(trait):animal, rcov = ~idh(trait):units, family = c("multinomial2", "multinomial2"), data=skew, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, prior=pr.skew, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    skewCheckMod.start$VCV[i-300,] <- skewCheckMod$VCV[1,]
    skewCheckMod.start$Sol[i-300,] <- skewCheckMod$Sol[1,]
    skewCheckMod.start$Liab[i-300,] <- skewCheckMod$Liab[1,]
  }
}
skewCheckModA <- skewCheckMod.start
skewCheckModB <- skewCheckMod.start
skewCheckModC <- skewCheckMod.start
save(skewCheckModA, file="results/skewCheckModA")
save(skewCheckModB, file="results/skewCheckModB")
save(skewCheckModC, file="results/skewCheckModC")
load("results/skewCheckModA")

summary(skewCheckModA)
posterior.mode(inv.logit(skewCheckModA$Sol[,1:8]))
HPDinterval(inv.logit(skewCheckModA$Sol[,1:8]))

# females: cats 1, 2, and 4 vs. cat 3
table(skewCheckModA$Sol[,1] > skewCheckModA$Sol[,5]) / length(skewCheckModA$Sol[,1])   # FALSE = 1, TRUE = 0
table(skewCheckModA$Sol[,3] > skewCheckModA$Sol[,5]) / length(skewCheckModA$Sol[,1])   # FALSE = 99, TRUE = 0.01
table(skewCheckModA$Sol[,7] > skewCheckModA$Sol[,5]) / length(skewCheckModA$Sol[,1])   # FALSE = 0.99, TRUE = 0.01
# males: cats 1, 2, and 4 vs. cat 3
table(skewCheckModA$Sol[,2] > skewCheckModA$Sol[,6]) / length(skewCheckModA$Sol[,1])   # FALSE = 1, TRUE = 0
table(skewCheckModA$Sol[,4] > skewCheckModA$Sol[,6]) / length(skewCheckModA$Sol[,1])   # FALSE = 1, TRUE = 0
table(skewCheckModA$Sol[,8] > skewCheckModA$Sol[,6]) / length(skewCheckModA$Sol[,1])   # FALSE = 0.99, TRUE = 0.1


## Group size ##
# add Riehl classification to group size data
groupSize$Riehl <- as.factor(polyData$Riehl[match(groupSize$animal, polyData$animal)])
table(groupSize$Riehl)		# cat 1 = 55, cat 2 = 31, cat 3 = 21, cat 4 = 17

# multi-response model on multiple trees
pr.gs <- list(R = list(V=diag(2), n = 1.002), G = list(G1=list(V=diag(2), n=1.002)))
INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
gsCheckMod.start <- MCMCglmm(cbind(log(mean), max) ~ trait:Riehl-1, random = ~idh(trait):animal, rcov = ~idh(trait):units, family=c("gaussian", "poisson"), data=groupSize, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, prior=pr.gs, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
gsCheckMod <- gsCheckMod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=gsCheckMod$Liab[1,], R=diag(gsCheckMod$VCV[1,3:4]), G=list(G1=diag(gsCheckMod$VCV[1,1:2])))
  gsCheckMod <- MCMCglmm(cbind(log(mean), max) ~ trait:Riehl-1, random = ~idh(trait):animal, rcov = ~idh(trait):units, family=c("gaussian", "poisson"), data=groupSize, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, prior=pr.gs, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    gsCheckMod.start$VCV[i-300,] <- gsCheckMod$VCV[1,]
    gsCheckMod.start$Sol[i-300,] <- gsCheckMod$Sol[1,]
    gsCheckMod.start$Liab[i-300,] <- gsCheckMod$Liab[1,]
  }
}
gsCheckModA <- gsCheckMod.start
gsCheckModB <- gsModel.start
gsCheckModC <- gsCheckMod.start
save(gsCheckModA, file="results/gsCheckModA")
save(gsCheckModB, file="results/gsCheckModB")
save(gsCheckModC, file="results/gsCheckModC")
load("results/gsCheckModA")
summary(gsCheckModA)
posterior.mode(exp(gsCheckModA$Sol[,c(1:8)]))
HPDinterval(exp(gsCheckModA$Sol[,c(1:8)]))

# parameters are means, not differences, therefore use table(a > b)/length(a) to test
# mean: cats 1, 2, and 4 vs. cat 3
table(gsCheckModA$Sol[,1] > gsCheckModA$Sol[,5]) / length(gsCheckModA$Sol[,1])   # FALSE = 0.21, TRUE = 0.79
table(gsCheckModA$Sol[,3] > gsCheckModA$Sol[,5]) / length(gsCheckModA$Sol[,1])   # FALSE = 0.09, TRUE = 0.91
table(gsCheckModA$Sol[,7] > gsCheckModA$Sol[,5]) / length(gsCheckModA$Sol[,1])   # FALSE = 0.01, TRUE = 0.99
# max: cats 1, 2, and 4 vs. cat 3
table(gsCheckModA$Sol[,2] > gsCheckModA$Sol[,6]) / length(gsCheckModA$Sol[,1])   # FALSE = 0.03, TRUE = 0.97
table(gsCheckModA$Sol[,4] > gsCheckModA$Sol[,6]) / length(gsCheckModA$Sol[,1])   # FALSE = 0.01, TRUE = 0.99
table(gsCheckModA$Sol[,8] > gsCheckModA$Sol[,6]) / length(gsCheckModA$Sol[,1])   # FALSE = 0.01, TRUE = 0.99


## Reproductive specialisation ##
# add Riehl classification to specialization data
uniqueEffects$Riehl <- as.factor(polyData$Riehl[match(uniqueEffects$animal, polyData$animal)])

# trim data
zRmiData <- uniqueEffects[-which(is.na(uniqueEffects$zRmi)),]	# 41 species
table(zRmiData$Riehl)		# cat 1 = 13, cat 2 = 9, cat 3 = 11, cat 4 = 7 (and 1 NA)

# run the model on multiple trees
pr.zRmi <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
MEV <- zRmiData$varRmi
INtree <- inverseA(dolTrees[[1]], nodes="TIPS")
zRmiCheckModel.start <- MCMCglmm(zRmi ~ Riehl-1, random= ~animal, family="gaussian", data=zRmiData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=pr.zRmi, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
zRmiCheckModel <- zRmiCheckModel.start
for(i in 1:1300){
  INtree <- inverseA(dolTrees[[i]], nodes="TIPS")
  start <- list(Liab=zRmiCheckModel$Liab[1,], R=zRmiCheckModel$VCV[1,3], G=list(G1=zRmiCheckModel$VCV[1,1]))
  zRmiCheckModel <- MCMCglmm(zRmi ~ Riehl-1, random= ~animal, family="gaussian", data=zRmiData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=pr.zRmi, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    zRmiCheckModel.start$VCV[i-300,] <- zRmiCheckModel$VCV[1,]
    zRmiCheckModel.start$Sol[i-300,] <- zRmiCheckModel$Sol[1,]
    zRmiCheckModel.start$Liab[i-300,] <- zRmiCheckModel$Liab[1,]
  }
}
zRmiCheckModelA <- zRmiCheckModel.start
zRmiCheckModelB <- zRmiCheckModel.start
zRmiCheckModelC <- zRmiCheckModel.start
save(zRmiCheckModelA, file="results/zRmiCheckModelA")
save(zRmiCheckModelB, file="results/zRmiCheckModelB")
save(zRmiCheckModelC, file="results/zRmiCheckModelC")
load("results/zRmiCheckModelA")

summary(zRmiCheckModelA)
posterior.mode(zRmiCheckModelA$Sol[,1:4])
HPDinterval(zRmiCheckModelA$Sol[,1:4])

# females: cats 1, 2, and 4 vs. cat 3
table(zRmiCheckModelA$Sol[,1] > zRmiCheckModelA$Sol[,3]) / length(zRmiCheckModelA$Sol[,1])   # FALSE = 0.01, TRUE = 0.99
table(zRmiCheckModelA$Sol[,2] > zRmiCheckModelA$Sol[,3]) / length(zRmiCheckModelA$Sol[,1])   # FALSE = 0.01, TRUE = 0.99
table(zRmiCheckModelA$Sol[,4] > zRmiCheckModelA$Sol[,3]) / length(zRmiCheckModelA$Sol[,1])   # FALSE = 0.02, TRUE = 0.98


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 2.2: DIFFERENT MEASURES OF FECUNDITY AND MATERNAL CARE ###

## (a). Fecundity Measures ##
effects$ID <- effects$animal
birdCor <- vcv.phylo(dolTrees[[1]], cor=TRUE)
fecundModel <- rma.mv(zRmi ~ zRmiMeasure*groupFormation, V=varRmi, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=effects)
summary(fecundModel)    # zRmi = -0.35 < -0.01 < 0.34; between-study var = 0.003, phylo var = 0.15

# correlations when different fecundity measures made on the same species
Feffects <- effects[,c("animal", "zRmi", "zRmiMeasure")]
Feffects <- Feffects[-which(is.na(Feffects$zRmi)),]
meltedClutch <- reshape(Feffects, idvar="animal", timevar="zRmiMeasure", direction="wide")
meltedClutch$gf <- as.factor(polyData$groupFormation[match(meltedClutch$animal, polyData$animal)])
names(meltedClutch) <- c("animal", "clutch", "reNest", "egg", "gf")
cor.test(meltedClutch$clutch[which(meltedClutch$gf == "Family")], meltedClutch$reNest[which(meltedClutch$gf == "Family")], method ="spearman")			# 0.60 (note - 0.24 using Spearman's rho)
# data too few for non-family groups: 2 species with clutch ~ reNest and 1 species with clutch ~ egg

# clutch size ~ group formation
clutchEffects <- effects[which(effects$zRmiMeasure == "clutch size"),]
clutchEffects$ID <- clutchEffects$animal		# to account for repeated measures
pr.clutch <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002), G2=list(V=1,nu=0.002)))
MEV <- clutchEffects$varRmi
INtree <- inverseA(dolTrees[[1]], nodes="TIPS")
clutchModelX <- MCMCglmm(zRmi ~ groupFormation-1, random= ~animal+ID, family="gaussian", data=clutchEffects, pl=TRUE, slice=TRUE, nitt=1100000, thin=1000, burnin=100000, mev=MEV, prior=pr.clutch, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
summary(clutchModelX)
posterior.mode(clutchModelX$Sol)	# family = 0.21, non-family = -0.20
HPDinterval(clutchModelX$Sol)		# family = 0.05 to 0.36, non-family = -0.40 to -0.06
table(clutchModelX$Sol[,1] > clutchModelX$Sol[,2]) / length(clutchModelX$Sol[,1])		# FALSE = 0, TRUE = 1

# clutch size ~ zRme in family and non-family groups
plot(zRmi ~ zRme, clutchEffects, pch=21, bg=ifelse(effects$groupFormation == "Family", "dodgerblue", "red"))
uniqueClutch <- summaryBy(zRme + RmeN + varRme + zRmi + RmiN + varRmi + zRsx + RsxN + varRsx + groupFormation ~ animal, data=clutchEffects, keep.names=TRUE, na.rm=TRUE)
uniqueClutch$ID <- uniqueClutch$animal
summary(lm(zRmi ~ zRme * groupFormation, uniqueClutch[-which(is.na(uniqueClutch$zRme)),]))
clutchModel <- rma.mv(zRmi ~ zRme*groupFormation+log(RmeN), V=varRmi, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=uniqueClutch)
summary(clutchModel)    # slope difference = 0.5 < 0.9 < 1.3; between-study var = 0.00, phylo var = 0.01
clutchEffects2 <- clutchEffects[-which(is.na(clutchEffects$zRme)),]
MEV <- clutchEffects2$varRmi
INtree <- inverseA(dolTrees[[1]], nodes="TIPS")
clutchModel2 <- MCMCglmm(zRmi ~ zRme*groupFormation + log(RmeN), random= ~animal+ID, family="gaussian", data=clutchEffects2, pl=TRUE, slice=TRUE, nitt=1100000, thin=1000, burnin=100000, mev=MEV, prior=pr.clutch, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
summary(clutchModel2)
posterior.mode(clutchModel2$Sol[,5])	# difference in slope = 0.66
HPDinterval(clutchModel2$Sol[,5])		# 0.21 < 0.66 < 1.01


## (b). Measures of Maternal Care ##

# correlations when different care measures made on the same species
Ceffects <- effects[,c("animal", "zRme", "zRmeMeasure")]
Ceffects <- Ceffects[-which(is.na(Ceffects$zRme)),]
meltedFeeding <- reshape(Ceffects, idvar="animal", timevar="zRmeMeasure", direction="wide")
meltedFeeding$gf <- as.factor(polyData$groupFormation[match(meltedFeeding$animal, polyData$animal)])
names(meltedFeeding) <- c("animal", "brooding", "feeding", "incubating", "gf")
cor.test(meltedFeeding$feeding, meltedFeeding$incubating, method ="spearman")			# 0.80 (note - 0.61 using Spearman's rho)

# feeding~ group formation
feedingEffects <- effects[which(effects$zRmeMeasure == "feeding"),]
feedingEffects$ID <- feedingEffects$animal		# to account for repeated measures
pr.feeding <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002), G2=list(V=1,nu=0.002)))
MEV <- feedingEffects$varRme
INtree <- inverseA(dolTrees[[1]], nodes="TIPS")
feedingModelX <- MCMCglmm(zRme ~ groupFormation-1, random= ~animal+ID, family="gaussian", data=feedingEffects, pl=TRUE, slice=TRUE, nitt=1100000, thin=1000, burnin=100000, mev=MEV, prior=pr.feeding, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
summary(feedingModelX)
posterior.mode(feedingModelX$Sol)	# family = -0.32, non-family = -0.10
HPDinterval(feedingModelX$Sol)		# family = -0.53 to -0.15, non-family = -0.51 to 0.04
table(feedingModelX$Sol[,1] > feedingModelX$Sol[,2]) / length(feedingModelX$Sol[,1])		# FALSE = 0.82, TRUE = 18


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 2.3: UNCERTAINTY IN ANCESTRAL STATE ESTIMATION AND DIFFERENT POLYANDRY MEASURES ###

## the following results were obtained by repeating parts 1.1(a) and 1.1(b) of the main analyses ##
## we varied the cut-off (lines 182 to 190) and changed our measure of polyandry (lines 297 to 337) ##

## EGP at 0.67 cut-off ##
load("results/polyModel3AEGP")
# chain convergence
hist(polyModel3AEGP$Liab)
plot(polyModel3AEGP$VCV)     		# animal and units close to 0
plot(polyModel3AEGP$Sol[,1:6])   	# 4/6 chains well mixed
autocorr(polyModel3AEGP$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(polyModel3AEGP$Sol)   		# correlation between successive samples < 0.1 for all components
polyModel3Sols <- mcmc.list(list(polyModel3AEGP$Sol, polyModel3BEGP$Sol, polyModel3CEGP$Sol))
plot(polyModel3EGPSols)
gelman.diag(polyModel3EGPSols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(polyModel3EGPA$VCV)    # animal failed halfwidth
heidel.diag(polyModel3EGPA$Sol)    # intercept failed halfwidth
summary(polyModel3AEGP)		
posterior.mode(inv.logit(polyModel3AEGP$Sol[,1:6]))
HPDinterval(inv.logit(polyModel3AEGP$Sol[,1:6]))		# wide ~ 0.0 to 1.0 for most parameters
# origin of nonCoop (0.19) vs. family (0.07)
table(polyModel3AEGP$Sol[, 2] > polyModel3AEGP$Sol[, 3]) / length(polyModel3AEGP$Sol[, 1])    # FALSE = 0.01, TRUE = 0.99
# origin of nonCoop (0.19) vs nonFamily (0.10)
table(polyModel3AEGP$Sol[, 2] > polyModel3AEGP$Sol[, 5]) / length(polyModel3AEGP$Sol[, 1])   # FALSE = 0.35, TRUE = 0.65
# origin of family vs. origin of nonfamily
table(polyModel3AEGP$Sol[,3] > polyModel3AEGP$Sol[,5]) / length(polyModel3AEGP$Sol[,1])    # FALSE = 0.75, TRUE = 0.25

## EPP at 0.33 cut-off ##
load("results/polyModelA")
hist(polyModelA$Liab)
plot(polyModelA$VCV)     			# animal > units
plot(polyModelA$Sol[,1:6])   		# 3/6 chains well mixed
autocorr(polyModelA$VCV)   			# correlation between successive samples < 0.1 for all components
autocorr(polyModelA$Sol[,1:6])  	# correlation between successive samples < 0.1 for all components
polyModelSols <- mcmc.list(list(polyModelA$Sol, polyModelB$Sol, polyModelC$Sol))
plot(polyModelSols)
gelman.diag(polyModelSols)     		# upper CI = 1.00 for intercept suggesting convergence
heidel.diag(polyModelA$VCV)    		# animal and units passed halfwidth
heidel.diag(polyModelA$Sol[,1:6])   # intercept failed halfwidth
summary(polyModelA)
posterior.mode(inv.logit(polyModelA$Sol[,1:6]))
HPDinterval(inv.logit(polyModelA$Sol[,1:6]))
# origin of nonCoop (0.18) vs. family (0.12)
table(polyModelA$Sol[, 2] > polyModelA$Sol[, 3]) / length(polyModelA$Sol[, 1])    # FALSE = 0.30, TRUE = 0.70
# origin of nonCoop (0.18) vs nonFamily (0.32)
table(polyModelA$Sol[, 2] > polyModelA$Sol[, 5]) / length(polyModelA$Sol[, 1])   # FALSE = 0.96, TRUE = 0.04
# origin of family vs. origin of nonfamily
table(polyModelA$Sol[,3] > polyModelA$Sol[,5]) / length(polyModelA$Sol[,1])    # FALSE = 0.95, TRUE = 0.05

## EGP at 0.33 cut-off ##
load("results/polyModelAEGP")
hist(polyModelAEGP$Liab)
plot(polyModelAEGP$VCV)     			# animal > units
plot(polyModelAEGP$Sol[,1:6])   		# 5/6 chains well mixed
autocorr(polyModelAEGP$VCV)   			# correlation between successive samples < 0.1 for all components
autocorr(polyModelAEGP$Sol[,1:6])  		# correlation between successive samples < 0.1 for all components
polyModelSols <- mcmc.list(list(polyModelAEGP$Sol, polyModelBEGP$Sol, polyModelCEGP$Sol))
plot(polyModelSols)
gelman.diag(polyModelSols)     			# upper CI = 1.00 for intercept suggesting convergence
heidel.diag(polyModelAEGP$VCV)    		# animal and units passed halfwidth
heidel.diag(polyModelAEGP$Sol[,1:6])   	# intercept passed halfwidth
summary(polyModelAEGP)
posterior.mode(inv.logit(polyModelAEGP$Sol[,1:6]))
HPDinterval(inv.logit(polyModelAEGP$Sol[,1:6]))
# origin of nonCoop (0.19) vs. family (0.9)
table(polyModelAEGP$Sol[, 2] > polyModelAEGP$Sol[, 3]) / length(polyModelAEGP$Sol[, 1])    # FALSE = 0.01, TRUE = 0.99
# origin of nonCoop (0.19) vs nonFamily (0.1)
table(polyModelAEGP$Sol[, 2] > polyModelAEGP$Sol[, 5]) / length(polyModelAEGP$Sol[, 1])   # FALSE = 0.37, TRUE = 0.63
# origin of family vs. origin of nonfamily
table(polyModelAEGP$Sol[,3] > polyModelAEGP$Sol[,5]) / length(polyModelA$Sol[,1])    # FALSE = 0.85, TRUE = 0.15

## EPP at 0.9 cut-off ##
load("results/polyModel2A")
# chain convergence
hist(polyModel2A$Liab)
plot(polyModel2A$VCV)     				# animal > units
plot(polyModel2A$Sol[,1:6])   			# 2/6 chains well mixed
autocorr(polyModel2A$VCV)   			# correlation between successive samples < 0.1 for all components
autocorr(polyModel2A$Sol[,1:6])  		# correlation between successive samples < 0.1 for all components
polyModel2Sols <- mcmc.list(list(polyModel2A$Sol, polyModel2B$Sol, polyModel2C$Sol))
plot(polyModel2Sols)
gelman.diag(polyModel2Sols)    	 		# upper CI = 1.00 for intercept suggesting convergence
heidel.diag(polyModel2A$VCV)   			# animal and units passed halfwidth
heidel.diag(polyModel2A$Sol[,1:6])   	# intercept passed halfwidth
summary(polyModel2A)
posterior.mode(inv.logit(polyModel2A$Sol[,1:6]))
HPDinterval(inv.logit(polyModel2A$Sol[,1:6]))		# wide ~ 0.0 to 1.0 for most parameters
# origin of nonCoop (0.20) vs. family (0.14)
table(polyModel2A$Sol[, 2] > polyModel2A$Sol[, 3]) / length(polyModel2A$Sol[, 1])    # FALSE = 0.41, TRUE = 0.59
# origin of nonCoop (0.20) vs nonFamily (0.99)
table(polyModel2A$Sol[, 2] > polyModel2A$Sol[, 5]) / length(polyModel2A$Sol[, 1])   # FALSE = 0.90, TRUE = 0.10
# origin of family vs. origin of nonfamily
table(polyModel2A$Sol[,3] > polyModel2A$Sol[,5]) / length(polyModel2A$Sol[,1])    # FALSE = 0.92, TRUE = 0.08

## EGP at 0.9 cut-off ##
load("results/polyModel2AEGP")
# chain convergence
hist(polyModel2AEGP$Liab)
plot(polyModel2AEGP$VCV)     				# animal > units
plot(polyModel2AEGP$Sol[,1:6])   			# 2/6 chains well mixed
autocorr(polyModel2AEGP$VCV)   				# correlation between successive samples < 0.1 for all components
autocorr(polyModel2AEGP$Sol[,1:6])  		# correlation between successive samples < 0.1 for all components
polyModel2Sols <- mcmc.list(list(polyModel2AEGP$Sol, polyModel2BEGP$Sol, polyModel2CEGP$Sol))
plot(polyModel2Sols)
gelman.diag(polyModel2Sols)    	 			# upper CI = 1.00 for intercept suggesting convergence
heidel.diag(polyModel2AEGP$VCV)   			# animal and units passed halfwidth
heidel.diag(polyModel2AEGP$Sol[,1:6])   	# intercept passed halfwidth
summary(polyModel2AEGP)
posterior.mode(inv.logit(polyModel2AEGP$Sol[,1:6]))
HPDinterval(inv.logit(polyModel2AEGP$Sol[,1:6]))		# wide ~ 0.0 to 1.0 for most parameters
# origin of nonCoop (0.22) vs. family (0.7)
table(polyModel2AEGP$Sol[, 2] > polyModel2AEGP$Sol[, 3]) / length(polyModel2AEGP$Sol[, 1])  # FALSE = 0.3, TRUE = 0.97
# origin of nonCoop (0.22) vs nonFamily (0.5)
table(polyModel2AEGP$Sol[, 2] > polyModel2AEGP$Sol[, 5]) / length(polyModel2AEGP$Sol[, 1])  # FALSE = 0.34, TRUE = 0.66
# origin of family vs. origin of nonfamily
table(polyModel2AEGP$Sol[,3] > polyModel2AEGP$Sol[,5]) / length(polyModel2AEGP$Sol[,1])		# FALSE = 0.70, TRUE = 0.30

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #
#                                                                                     					 #
#                            		END - thanks for reading! p. 										 #
#                                                                                     					 #  
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #