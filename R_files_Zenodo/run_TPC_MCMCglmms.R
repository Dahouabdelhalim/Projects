#!/usr/bin/Rscript

#######################################################################
# This R script fits all the MCMCglmms described in this study to our #
# entire TPC dataset. 2 chains (a and b) are run per model.           #
#                                                                     #
# Requirements for running:                                           #
#                                                                     #
#   ->   the TPC dataset in the working directory                     #
#                                                                     #
#   ->   the calibrated phylogeny in the working directory            #
#######################################################################

library(ape)
library(MCMCglmm)

# Read the TPC dataset.
dataset <- read.csv('dataset_for_TPC_MCMCglmms.csv', stringsAsFactors = FALSE)

# Read the calibrated tree.
tree <- read.nexus('final_calibrated_tree.phy')

# Drop species for which we do not have TPC parameter estimates.
tree_for_GLMs <- drop.tip(
	tree, 
	tree$tip.label[!(tree$tip.label %in% dataset$Species_stand)])
inv.phylo.1 <- inverseA(tree_for_GLMs, nodes="ALL", scale=TRUE)

# Collect all estimates of measurement error variance.
MEVs <- c(dataset$B_T_ref_1_4_sev, dataset$ln_E_sev, 
	dataset$T_pk_squared_sev, dataset$ln_B_pk_sev, 
	dataset$ln_E_D_sev, dataset$ln_W_op_sev)

# When TPC parameter estimates are missing (NA), set their measurement 
# error variance to 0. Otherwise, the code breaks.
MEVs[is.na(MEVs)] <- 0
MEVs[is.nan(MEVs)] <- 0

#######################
# Fit all the models. #
#######################

set.seed(1337)
model.1a_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(1604)
model.1b_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(1904)
model.2a_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(Latitude) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(0408)
model.2b_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(Latitude) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(1990)
model.3a_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:habitat + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(2017)
model.3b_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:habitat + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(1871)
model.4a_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:habitat + trait:abs(Latitude) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(1917)
model.4b_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:habitat + trait:abs(Latitude) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(26062017)
model.5a_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(Latitude, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(19072017)
model.5b_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(Latitude, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(19902017)
model.6a_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:habitat + trait:poly(Latitude, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

set.seed(19661990)
model.6b_MEV <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:habitat + trait:poly(Latitude, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo.1$Ainv),
	mev = MEVs,
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	rcov=~us(trait):units,
	nitt=100000000,
	burnin=10000000,
	thin=1000,
	verbose = TRUE)

# Save model fits to a file for analysis.
save.image(file = "TPC_MCMCglmms.Rda")
