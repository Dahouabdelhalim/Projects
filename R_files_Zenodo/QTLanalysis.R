# QTL mapping with R/qtl.
# For more details on R/qtl, see: Broman KW, Sen S (2009) A Guide to QTL Mapping with R/qtl Springer, New York.

# The codes for the QTL analysis is given for two phenotypes:
# - estsuc_2012 -> as an example for all phenotypes with a normal distribution (before or after quantile-normalization)
# - adultsurv_2012_bin -> the only phenotype with a binary distribution

#### Load the data ####
rm(list=ls())
library(qtl)
mydata <-read.cross("csv", dir="", "RIL_phenotypes_genotypes_QTLanalysis.csv", na.strings=NA, genotypes = c("a", "b"), alleles=c("It","Sw"))
rildata<-convert2riself(mydata)

#### Quantile-normalize subset-data ####
quantnorm<-function(x) {
  n=sum(!is.na(x),na.rm=T)
  x=rank(x)/(n+1)
  x=qnorm(x)
  x[is.infinite(x)]=NA
  x }

names(rildata$pheno)
rildata$pheno$qn.estsuc_2012<-quantnorm(rildata$pheno$estsuc_2012)


#### Run the permutations ####
# OBS! This code is written to run the permutations on 8 clusters
# The binary permutations take a very long time and were therefore divided in 10x1000 permutations

rildata<-calc.genoprob(rildata, step=2, error.prob=.0001, map.function="kosambi",stepwidth="max", off.end=0)


## Permutations on quantile-normalized subsets of data ##
perm2_qn.estsuc_2012 <-scantwo(rildata, pheno.col="qn.estsuc", model="normal",  method="hk",  n.perm=10000, 
                               addcovar=NULL, intcovar=NULL, weights=NULL,clean.output=FALSE, clean.nmar=FALSE, 
                               clean.distance=FALSE, incl.markers=TRUE,  assumeCondIndep=FALSE,  verbose=FALSE)
save(perm2_qn.estsuc_2012, file="qn.estsuc_2012_10k.RData")


## Binary permutations ##
# These permutation take a long computation time!
# These permutations were performed 10 times with 1000 permutations each time, and seed set was changed for each time
# After that, the permutations were put together using this code:
#perm2_adultsurv_2012_bin<-c(perm2_adultsurv_2012_bin_1, perm2_adultsurv_2012_bin_2, perm2_adultsurv_2012_bin_3, perm2_adultsurv_2012_bin_4, perm2_adultsurv_2012_bin_5,perm2_adultsurv_2012_bin_6, perm2_adultsurv_2012_bin_7, perm2_adultsurv_2012_bin_8, perm2_adultsurv_2012_bin_9, perm2_adultsurv_2012_bin_10)
#save(perm2_adultsurv_2012_bin, file="perm2_adultsurv_2012_bin_10k.RData")

set.seed(071510)
perm2_adultsurv_2012_bin<-scantwo(rildata, pheno.col="adultsurv_2012_bin", model="binary",  method="hk",  n.perm=1000, n.cluster=8,
                                  addcovar=NULL, intcovar=NULL, weights=NULL, clean.output=FALSE, clean.nmar=FALSE, 
                                  clean.distance=FALSE, incl.markers=TRUE,  assumeCondIndep=FALSE,  verbose=TRUE)
save(perm2_adultsurv_2012_bin, file="perm2_adultsurv_2012_bin_1.RData")

#### Calculate LOD threshold penalties for alpha=0.05 ####
pen_qn.estsuc_2012 <-calc.penalties(perm2_qn.estsuc_2012, alpha=.05)
pen_adultsurv_2012_bin<-calc.penalties(perm2_adultsurv_2012_bin, alpha=.05)

#### Run stepwiseqtl for multiple QTL model selection ####
# Model is binary for binary phenotypes
# Quantile-normalized phenotypes are used for phenotypes which did not have a normal distribution before quantile-normalization.
model_estsuc_2012 <- stepwiseqtl(rildata, pheno.col="qn.estsuc_2012", penalties=pen_qn.estsuc_2012, 
                               method="hk", model="normal", max.qtl=14, covar=NULL, scan.pairs=F, additive.only=F, 
                               keeplodprofile=T, keeptrace=T, refine.locations=T, verbose=F)
summary(model_estsuc_2012) 

model_adultsurv_2012_bin  <- stepwiseqtl(rildata, pheno.col="adultsurv_2012_bin", penalties=pen_adultsurv_2012_bin,
                                          method="hk", model="binary", max.qtl=14, covar=NULL, scan.pairs=F, additive.only=F, 
                                          keeplodprofile=T, keeptrace=T, refine.locations=T, verbose=F) 
summary(model_adultsurv_2012_bin) # No QTL detected


#### Fit the model ####
# using the untransformed phenotypes
fit_model_estsuc_2012 <- fitqtl(rildata, pheno.col="estsuc_2012", qtl=model_estsuc_2012,  
                              formula=formula(model_estsuc_2012), method="hk", 
                              model="normal", covar=NULL, get.ests=T, dropone=T)
summary(fit_model_estsuc_2012) 


#### Get the bayesian credible intervals ####
for (i in 1:length(model_estsuc_2012)) {
    print(bayesint(model_estsuc_2012, prob=0.95, qtl.index=i))}

