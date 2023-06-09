#!/usr/bin/Rscript

####################################################################
# This R script fits the B_0 or B_pk vs volume regressions of this #
# study. 2 chains (a and b) are run per model.                     #
#                                                                  #
# Requirements for running:                                        #
#                                                                  #
#   ->   the TPC dataset in the working directory                  #
#                                                                  #
#   ->   the calibrated phylogeny in the working directory         #
####################################################################

library(ape)
library(MCMCglmm)

# Read the B_0, B_pk, and cell volume dataset.
dataset <- read.csv('dataset_of_B0_B_pk_and_cell_volume.csv', stringsAsFactors = FALSE)

# Read the calibrated tree.
tree <- read.nexus('final_calibrated_tree.phy')

# Drop species for which we do not have TPC parameter estimates.
tree_for_GLMs <- drop.tip(
	tree,
	tree$tip.label[!(tree$tip.label %in% dataset$Species_stand)]
)

# Calculate the inverse of the phylogenetic covariance matrix.
inv.phylo <- inverseA(tree_for_GLMs, nodes="ALL", scale=TRUE)

# Collect all estimates of measurement error variance.
MEVs_B_0 <- c(dataset$B_T_ref_1_4_sev)
MEVs_B_pk <- c(dataset$ln_B_pk_sev)

# When TPC parameter estimates are missing (NA), set their measurement 
# error variance to 0. Otherwise, the code breaks.
MEVs_B_0[is.na(MEVs_B_0)] <- 0
MEVs_B_0[is.nan(MEVs_B_0)] <- 0
MEVs_B_pk[is.na(MEVs_B_pk)] <- 0
MEVs_B_pk[is.nan(MEVs_B_pk)] <- 0

#######################
# Fit all the models. #
#######################

fit_B_0_1a <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_1b <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_1a_non_phy <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~Species_stand,
	family="gaussian",
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_1b_non_phy <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~Species_stand,
	family="gaussian",
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_2a <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~us(ln_volume):Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_2b <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~us(ln_volume):Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_2a_non_phy <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~us(ln_volume):Species_stand,
	family="gaussian",
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_2b_non_phy <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~us(ln_volume):Species_stand,
	family="gaussian",
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_3a <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~Species_stand + us(ln_volume):Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002), G2=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_3b <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~Species_stand + us(ln_volume):Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002), G2=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_3a_non_phy <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~Species_stand + us(ln_volume):Species_stand,
	family="gaussian",
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002), G2=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_0_3b_non_phy <- MCMCglmm(
	B_T_ref_1_4 ~ ln_volume,
	random=~Species_stand + us(ln_volume):Species_stand,
	family="gaussian",
	mev = MEVs_B_0,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002), G2=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_1a <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_1b <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_1a_non_phy <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~Species_stand,
	family="gaussian",
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_1b_non_phy <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~Species_stand,
	family="gaussian",
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_2a <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~us(ln_volume):Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_2b <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~us(ln_volume):Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_2a_non_phy <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~us(ln_volume):Species_stand,
	family="gaussian",
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_2b_non_phy <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~us(ln_volume):Species_stand,
	family="gaussian",
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_3a <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~Species_stand + us(ln_volume):Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002), G2=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_3b <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~Species_stand + us(ln_volume):Species_stand,
	family="gaussian",
	ginverse=list(Species_stand = inv.phylo$Ainv),
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002), G2=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_3a_non_phy <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~Species_stand + us(ln_volume):Species_stand,
	family="gaussian",
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002), G2=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

fit_B_pk_3b_non_phy <- MCMCglmm(
	ln_B_pk ~ ln_volume,
	random=~Species_stand + us(ln_volume):Species_stand,
	family="gaussian",
	mev = MEVs_B_pk,
	prior=list(G=list(G1=list(V=diag(1),nu=1.002), G2=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
	data=dataset,
	rcov=~units,
	nitt=3000000,
	burnin=300000,
	thin=1000,
	verbose = TRUE)

# Save model fits to a file for analysis.
save.image(file = "B_0_B_pk_vs_cell_volume_MCMCglmms.Rda")
