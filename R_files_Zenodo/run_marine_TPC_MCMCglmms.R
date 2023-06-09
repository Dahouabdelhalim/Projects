#!/usr/bin/Rscript

#######################################################################
# This R script fits all the MCMCglmms described in this study to the #
# TPC dataset of marine phytoplankton. 2 chains (a and b) are run per #
# model.                                                              #
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
dataset <- read.csv('dataset_for_marine_TPC_MCMCglmms.csv', stringsAsFactors = FALSE)

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

model.1a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.1b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.2a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(Latitude) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.2b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(Latitude) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.3a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(Latitude, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.3b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(Latitude, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.4a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:IQR_temp_Eul + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.4b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:IQR_temp_Eul + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.5a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:Median_temp_Eul + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.5b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:Median_temp_Eul + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.6a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:Median_temp_Eul + trait:IQR_temp_Eul + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.6b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:Median_temp_Eul + trait:IQR_temp_Eul + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.7a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.7b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.8a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.8b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.9a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_temps + trait:d_2.5_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.9b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_temps + trait:d_2.5_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.10a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.10b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.11a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.11b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.12a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_lats + trait:abs(d_2.5_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.12b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_lats + trait:abs(d_2.5_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.13a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.13b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.14a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_lats + trait:poly(d_2.5_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.14b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_50_w_med_iqr_lats + trait:poly(d_2.5_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.15a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.15b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.16a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.16b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.17a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_temps + trait:d_2.5_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.17b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_temps + trait:d_2.5_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.18a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.18b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.19a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.19b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.20a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_lats + trait:abs(d_2.5_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.20b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_lats + trait:abs(d_2.5_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.21a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.21b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.22a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_lats + trait:poly(d_2.5_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.22b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_150_w_med_iqr_lats + trait:poly(d_2.5_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.23a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.23b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.24a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.24b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.25a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_temps + trait:d_2.5_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.25b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_temps + trait:d_2.5_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.26a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.26b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.27a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.27b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.28a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_lats + trait:abs(d_2.5_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.28b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_lats + trait:abs(d_2.5_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.29a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.29b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.30a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_lats + trait:poly(d_2.5_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.30b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_250_w_med_iqr_lats + trait:poly(d_2.5_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.31a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.31b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.32a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.32b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.33a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_temps + trait:d_2.5_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.33b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_temps + trait:d_2.5_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.34a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.34b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.35a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.35b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.36a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_lats + trait:abs(d_2.5_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.36b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_lats + trait:abs(d_2.5_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.37a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.37b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.38a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_lats + trait:poly(d_2.5_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.38b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_350_w_med_iqr_lats + trait:poly(d_2.5_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.39a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.39b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.40a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.40b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.41a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_temps + trait:d_2.5_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.41b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_temps + trait:d_2.5_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.42a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.42b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.43a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.43b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_2.5_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.44a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_lats + trait:abs(d_2.5_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.44b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_lats + trait:abs(d_2.5_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.45a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.45b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.46a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_lats + trait:poly(d_2.5_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.46b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_2.5_500_w_med_iqr_lats + trait:poly(d_2.5_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.47a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.47b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.48a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.48b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.49a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_temps + trait:d_50_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.49b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_temps + trait:d_50_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.50a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.50b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.51a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.51b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.52a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_lats + trait:abs(d_50_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.52b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_lats + trait:abs(d_50_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.53a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.53b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.54a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_lats + trait:poly(d_50_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.54b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_50_w_med_iqr_lats + trait:poly(d_50_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.55a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.55b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.56a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.56b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.57a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_temps + trait:d_50_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.57b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_temps + trait:d_50_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.58a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.58b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.59a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.59b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.60a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_lats + trait:abs(d_50_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.60b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_lats + trait:abs(d_50_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.61a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.61b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.62a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_lats + trait:poly(d_50_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.62b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_150_w_med_iqr_lats + trait:poly(d_50_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.63a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.63b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.64a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.64b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.65a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_temps + trait:d_50_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.65b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_temps + trait:d_50_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.66a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.66b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.67a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.67b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.68a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_lats + trait:abs(d_50_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.68b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_lats + trait:abs(d_50_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.69a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.69b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.70a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_lats + trait:poly(d_50_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.70b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_250_w_med_iqr_lats + trait:poly(d_50_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.71a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.71b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.72a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.72b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.73a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_temps + trait:d_50_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.73b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_temps + trait:d_50_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.74a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.74b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.75a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.75b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.76a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_lats + trait:abs(d_50_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.76b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_lats + trait:abs(d_50_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.77a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.77b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.78a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_lats + trait:poly(d_50_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.78b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_350_w_med_iqr_lats + trait:poly(d_50_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.79a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.79b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.80a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.80b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.81a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_temps + trait:d_50_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.81b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_temps + trait:d_50_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.82a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.82b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.83a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.83b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_50_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.84a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_lats + trait:abs(d_50_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.84b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_lats + trait:abs(d_50_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.85a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.85b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.86a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_lats + trait:poly(d_50_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.86b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_50_500_w_med_iqr_lats + trait:poly(d_50_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.87a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.87b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.88a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.88b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.89a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_temps + trait:d_100_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.89b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_temps + trait:d_100_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.90a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.90b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.91a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.91b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.92a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_lats + trait:abs(d_100_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.92b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_lats + trait:abs(d_100_50_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.93a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.93b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.94a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_lats + trait:poly(d_100_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.94b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_50_w_med_iqr_lats + trait:poly(d_100_50_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.95a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.95b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.96a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.96b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.97a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_temps + trait:d_100_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.97b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_temps + trait:d_100_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.98a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.98b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.99a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.99b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.100a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_lats + trait:abs(d_100_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.100b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_lats + trait:abs(d_100_150_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.101a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.101b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.102a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_lats + trait:poly(d_100_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.102b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_150_w_med_iqr_lats + trait:poly(d_100_150_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.103a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.103b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.104a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.104b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.105a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_temps + trait:d_100_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.105b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_temps + trait:d_100_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.106a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.106b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.107a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.107b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.108a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_lats + trait:abs(d_100_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.108b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_lats + trait:abs(d_100_250_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.109a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.109b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.110a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_lats + trait:poly(d_100_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.110b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_250_w_med_iqr_lats + trait:poly(d_100_250_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.111a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.111b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.112a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.112b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.113a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_temps + trait:d_100_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.113b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_temps + trait:d_100_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.114a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.114b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.115a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.115b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.116a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_lats + trait:abs(d_100_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.116b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_lats + trait:abs(d_100_350_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.117a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.117b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.118a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_lats + trait:poly(d_100_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.118b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_350_w_med_iqr_lats + trait:poly(d_100_350_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.119a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.119b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.120a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.120b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.121a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_temps + trait:d_100_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.121b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_temps + trait:d_100_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.122a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.122b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_lats + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.123a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.123b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:abs(d_100_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.124a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_lats + trait:abs(d_100_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.124b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_lats + trait:abs(d_100_500_w_med_med_lats) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.125a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.125b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.126a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_lats + trait:poly(d_100_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.126b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:d_100_500_w_med_iqr_lats + trait:poly(d_100_500_w_med_med_lats, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

######################

model.127a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_50_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.127b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_50_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.128a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_50_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.128b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_50_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.129a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_150_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.129b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_150_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.130a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_150_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.130b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_150_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.131a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_250_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.131b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_250_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.132a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_250_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.132b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_250_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.133a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_350_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.133b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_350_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.134a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_350_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.134b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_350_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.135a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_500_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.135b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_500_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.136a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_500_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.136b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_2.5_500_w_med_iqr_temps, 2, raw = TRUE) + trait:d_2.5_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.137a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_50_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.137b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_50_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.138a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_50_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.138b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_50_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.139a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_150_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.139b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_150_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.140a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_150_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.140b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_150_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.141a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_250_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.141b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_250_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.142a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_250_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.142b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_250_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.143a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_350_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.143b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_350_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.144a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_350_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.144b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_350_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.145a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_500_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.145b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_500_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.146a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_500_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.146b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_50_500_w_med_iqr_temps, 2, raw = TRUE) + trait:d_50_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.147a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_50_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.147b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_50_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.148a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_50_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.148b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_50_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_50_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.149a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_150_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.149b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_150_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.150a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_150_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.150b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_150_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_150_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.151a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_250_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.151b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_250_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.152a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_250_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.152b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_250_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_250_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.153a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_350_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.153b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_350_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.154a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_350_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.154b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_350_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_350_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.155a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_500_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.155b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_500_w_med_iqr_temps, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.156a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_500_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.156b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(d_100_500_w_med_iqr_temps, 2, raw = TRUE) + trait:d_100_500_w_med_med_temps + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.157a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(IQR_temp_Eul, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.157b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:poly(IQR_temp_Eul, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.158a <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:Median_temp_Eul + trait:poly(IQR_temp_Eul, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

model.158b <- MCMCglmm(
	cbind(B_T_ref_1_4, ln_E, T_pk_squared, ln_B_pk, ln_E_D, ln_W_op) ~ trait:Median_temp_Eul + trait:poly(IQR_temp_Eul, 2, raw = TRUE) + trait - 1,
	random=~us(trait):Species_stand,
	family=rep("gaussian", 6),
	ginverse=list(Species_stand=inv.phylo$Ainv),
	prior=list(G=list(G1=list(V=diag(6),nu=1.002)), R=list(V=diag(6),nu=1.002)),
	data=dataset,
	mev = MEVs,
	rcov=~us(trait):units,
	nitt=60000000,
	burnin=6000000,
	thin=1000,
	verbose = TRUE)

# Save model fits to a file for analysis.
save.image(file = "marine_TPC_MCMCglmms.Rda")
