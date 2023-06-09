#MVPORPH Multivariate Model Fit#
library("phytools")
library("ape")
library("maps")
library("corpcor")
library("subplex")
library("spam")
library("grid")
library("geiger")
library("mvMORPH")
library("picante")
library("phangorn")

#Importing tree and data, and prunning tree

setwd("PATH/TO/FOLDER")
#Import data
mytree <- read.nexus("myTree.tre") #Import tree
mydata <- read.csv("pPCScores.csv", row.names = 1) #Import pPC scores

shiftTime<- 20 #expected time for a shift in mode of evolution

#Prune tree and data
pruned.data <- treedata(mytree, mydata)$data
pruned.tree <- treedata(mytree, mydata)$phy

#Generate a imput tree for mode shift analyses
treeShift <- make.era.map(pruned.tree, c(0,shiftTime))
#Select pPC axis to be analysed
MainPCs <- pruned.data[,1:12]

#Calculate fit of alternative models of evolution to a multivariate dataset of continuous traits 
fitBM1 <- mvBM(tree = pruned.tree, data =MainPCs, model = "BM1" )
fitOU1 <- mvOU(tree = pruned.tree, data =MainPCs, model = "OU1" )
fitEB <- mvEB(tree = pruned.tree, data =MainPCs )
fitEBOUtreeShift <- mvSHIFT(treeShift, data = MainPCs, model = "EBOUi")
fitBMOUtreeShift <- mvSHIFT(treeShift, data = MainPCs, model = "BMOUi")
fitEBBMtreeShift <- mvSHIFT(treeShift, data = MainPCs, model = "EBBMi")

#return the Akaike weights for a set of fitted models
AICwpPCA <- aicw(c(fitBM1$AIC, fitOU1$AIC, fitEB$AIC, fitEBOUtreeShift$AIC, fitBMOUtreeShift$AIC, fitEBBMtreeShift$AIC))

