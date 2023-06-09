#!/bin/R

args=commandArgs(TRUE)

seedR=as.numeric(args[1])
phenocolR=as.character(args[2])
outfileR=args[3]

# load library
library("qtl")

# read in data (previously copied it to scratch)
PNPxNOMdat <- read.cross("csv", "", "/cluster/scratch/annafell/QTLstuff/GenoPheno_PNPxNOM_7more.csv", map.function="kosambi")
PNPxNOMdat<-jittermap(PNPxNOMdat)

# prepare for mapping (calc.genoprob with stepsize 3 to make it faster)
PNPxNOMdat <- calc.genoprob(PNPxNOMdat, step=3, error.prob=0.05, map.function="kosambi")

# run two QTL model with 250 permutations 
set.seed(seedR)
out <- scantwo(PNPxNOMdat, pheno.col=phenocolR, model="normal", method="hk", n.perm=250, addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family))
# save output
saveRDS(out, file=paste(outfileR))
