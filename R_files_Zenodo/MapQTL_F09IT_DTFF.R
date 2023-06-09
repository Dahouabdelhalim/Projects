

rm(list=ls())
library(qtl)

sweden_italy_cross<-read.cross("csvs",dir=".",
  genotypes=c("a","b"), 
  genfile="geno.csv", 
  phefile= "pheno.csv", 
  na.strings=c("NA","-"))

# Set the cross type to RIL
class(sweden_italy_cross)[1] <-"riself"

# Calculate genotype probabilities
sweden_italy.calc<-calc.genoprob(sweden_italy_cross,
step=2, error.prob=.0001, map.function="kosambi",
stepwidth="max", off.end=0)

# Do a one way scan
scan1_F09_IT_DTFF_noedge_Mean <- scanone(sweden_italy.calc, method="hk", pheno.col="F09_IT_DTFF_noedge_Mean")
perm_scan1_F09_IT_DTFF_noedge_Mean <- scanone(sweden_italy.calc, method="hk",n.perm=10000, pheno.col="F09_IT_DTFF_noedge_Mean")


# Summarize one way scan with permutations summary
summary(scan1_F09_IT_DTFF_noedge_Mean, perms=perm_scan1_F09_IT_DTFF_noedge_Mean, alpha=0.05)
np_F09_IT_DTFF_noedge_Mean <- summary(perm_scan1_F09_IT_DTFF_noedge_Mean)
np_F09_IT_DTFF_noedge_Mean

# Plot one way scan with LOD threshold from permutations
pdf(file="scanone_F09_IT_DTFF_noedge_Mean.pdf")
plot(scan1_F09_IT_DTFF_noedge_Mean, main="scan1_F09_IT_DTFF_noedge_Mean")
abline(h=np_F09_IT_DTFF_noedge_Mean[1], col="orange")
dev.off()

# Do a two way scan
scan2_QN_F09_IT_DTFF_noedge_Mean <- scantwo(sweden_italy.calc,
pheno.col="QN_F09_IT_DTFF_noedge_Mean", model="normal",  method="hk", addcovar=NULL, intcovar=NULL, weights=NULL,
clean.output=FALSE, clean.nmar=FALSE, clean.distance=FALSE,
incl.markers=TRUE,  assumeCondIndep=FALSE,  verbose=FALSE)

# Two way scan summary
summary(scan2_QN_F09_IT_DTFF_noedge_Mean)

# Plot the two way scan
pdf(file="scantwo_QN_F09_IT_DTFF_noedge_Mean.pdf")
plot(scan2_QN_F09_IT_DTFF_noedge_Mean, main="scan2_QN_F09_IT_DTFF_noedge_Mean")
dev.off()

# Save an R image contaning the objects thus far
save.image("QN_F09_IT_DTFF_noedge_Mean.RData")

# Perform permutations for two way scan
QN_F09_IT_DTFF_noedge_Mean_perm<-scantwo(sweden_italy.calc,
pheno.col="QN_F09_IT_DTFF_noedge_Mean", model="normal",  method="hk", n.perm=10000, n.cluster=20,
addcovar=NULL, intcovar=NULL, weights=NULL,
clean.output=FALSE, clean.nmar=FALSE, clean.distance=FALSE,
incl.markers=TRUE,  assumeCondIndep=FALSE,  verbose=FALSE)

# Summarize permutations
summary(QN_F09_IT_DTFF_noedge_Mean_perm)

# Calculate LOD threshold penalties at alpha=0.05 for additive QTL
QN_F09_IT_DTFF_noedge_Mean_pen1.05<-calc.penalties(QN_F09_IT_DTFF_noedge_Mean_perm, alpha=.05)

# Save an R image containing the objects thus far
save.image("QN_F09_IT_DTFF_noedge_Mean.RData")


# Perform the stepwise forward-backward regression procedure on quantile normalized data
#  at alpha=0.05 for additive and epistatic QTL, max.qtl is the maximum number of QTL allowed in the model,
#  should be set to be > actual number of QTL
EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean<- stepwiseqtl(sweden_italy.calc, method="hk", model="normal",
pheno.col="QN_F09_IT_DTFF_noedge_Mean", penalties=QN_F09_IT_DTFF_noedge_Mean_pen1.05[1:2], max.qtl=15, covar=NULL,
scan.pairs=FALSE, additive.only=F,keeplodprofile=TRUE, keeptrace=TRUE,
refine.locations=TRUE, verbose=TRUE)

# Summarize the stepwise model and plot the lod profiles of significant QTL
summary(EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean)
pdf(file="EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean.pdf")
plotLodProfile(EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean, main="EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean_LodProfile")
dev.off()

# Fit and refine the stepwise model on the quantile normalized data
EPI_HEAVY_fitqtl_QN_F09_IT_DTFF_noedge_Mean<-fitqtl(sweden_italy.calc, pheno.col="QN_F09_IT_DTFF_noedge_Mean",
method="hk", model="normal", covar=NULL,
qtl=EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean, formula=formula(EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean),
get.ests=TRUE, dropone=TRUE)

# Re-fit this model using the raw data
EPI_HEAVY_fitqtl_F09_IT_DTFF_noedge_Mean<-fitqtl(sweden_italy.calc, pheno.col="F09_IT_DTFF_noedge_Mean",
method="hk", model="normal", covar=NULL,
qtl=EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean, formula=formula(EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean),
get.ests=TRUE, dropone=TRUE)

# Save an R image containing the objects thus far
save.image("QN_F09_IT_DTFF_noedge_Mean.RData")

###############################################################
# Get ANOVA tables from fitqtl models
# Use the LOD scores and PVE from the quantile-normalized data
# Use the allelic effect sizes from the raw data (multiply these by 2 to get genotypic effects (RILs))
# Layout changes in R/qtl for n QTL >1 vs. n QTL = 1
# if statements extract the bits you want regardless of number of qtl found

if(length(EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean$name)>1) {summary(EPI_HEAVY_fitqtl_QN_F09_IT_DTFF_noedge_Mean)[1:2]} else {summary(EPI_HEAVY_fitqtl_QN_F09_IT_DTFF_noedge_Mean)[1]}

if(length(EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean$name)>1) {summary(EPI_HEAVY_fitqtl_F09_IT_DTFF_noedge_Mean)[3]} else {summary(EPI_HEAVY_fitqtl_F09_IT_DTFF_noedge_Mean)[2]}


# For loop to print Bayesian 95% credible intervals for each QTL in model

for (i in 1:length(EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean$name)) {
    
    print(bayesint(EPI_HEAVY_stepwise_QN_F09_IT_DTFF_noedge_Mean, prob=0.95, qtl.index=i, expandtomarkers=F))
    
}
