########################################################################################
# This script includes code used to fit phylogenetic generalized linear mixed models
# (PGLMMs) to test for relationships among trait divergence, sympatry, and time since
# common ancestry. 
########################################################################################

#load necessary libraries
library(ape)
library(MCMCglmm)
library(phytools)

#read in phylogeny and check that it is ultrametric
noto.tree<-read.tree("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/notothenioid_timetree.tre")
is.ultrametric(noto.tree)

# read in data files for each clade - each datafile includes all possible pairwise comparisons of species within each 
# clade along with sympatry/allopatry as binary variables, age of split, and pairwise Euclidean distances for each trait
chan_data <- read.csv("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/chan_PairwiseTraits.csv")
trem_data <- read.csv("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/trem_PairwiseTraits.csv")

############################################
# First complete analyses for icefishes
############################################

# square root (sqrt) transform response variables
c.sqrt_dietPC1<-sqrt(chan_data$dietPC1)
c.sqrt_dietPC2<-sqrt(chan_data$dietPC2)
c.sqrt_dietPC3<-sqrt(chan_data$dietPC3)
c.sqrt_shapePC1<-sqrt(chan_data$shapePC1)
c.sqrt_shapePC2<-sqrt(chan_data$shapePC2)
c.sqrt_shapePC3<-sqrt(chan_data$shapePC3)
c.sqrt_size<-sqrt(chan_data$size)

# sqrt- and z-transform covariate (age)
c.sqrt_age<-sqrt(chan_data$age)
c.Zage<-((c.sqrt_age - mean(c.sqrt_age))/(sd(c.sqrt_age)))

# create species identity random effects 
sp <- unique(chan_data$lineage1)
sp_id <- data.frame(1:length(sp), sp, stringsAsFactors=F)
names(sp_id) <- c("id", "sp")

chan_lineage_1 <- apply(chan_data, 1, function(x) sp_id$id[grep(x[1], sp_id$sp)])
chan_lineage_2 <- apply(chan_data, 1, function(x) sp_id$id[grep(x[2], sp_id$sp)])

# prune tree so that tips match trait data and retain one outgroup taxon to root phylogeny
chan_tree <- keep.tip(noto.tree, c(as.vector(sp), "20895_Gymnodraco"))
chan_tree <- root(chan_tree, "20895_Gymnodraco")
is.binary.tree(chan_tree)
is.ultrametric(chan_tree)
is.rooted(chan_tree)

# create node random effect from pruned tree
chan_node <- apply(chan_data, 1, function(x) getMRCA(chan_tree, c(x[1], x[2]))) 
chan_node <- chan_node-length(chan_tree$tip.label)
chan_node <- paste("Node", chan_node, sep="")

#create object with data on sympatry and allopatry
symp<-chan_data$sympatry

#add transformed variables, node IDs, and lineage IDs to single dataframes for each trait
chan_shape1 <- data.frame(c.sqrt_shapePC1, symp, c.Zage, chan_node, chan_lineage_1, chan_lineage_2, stringsAsFactors=F)
chan_shape1 <- na.omit(chan_shape1)
chan_shape2 <- data.frame(c.sqrt_shapePC2, symp, c.Zage, chan_node, chan_lineage_1, chan_lineage_2, stringsAsFactors=F)
chan_shape2 <- na.omit(chan_shape2)
chan_shape3 <- data.frame(c.sqrt_shapePC3, symp, c.Zage, chan_node, chan_lineage_1, chan_lineage_2, stringsAsFactors=F)
chan_shape3 <- na.omit(chan_shape3)
chan_diet1 <- data.frame(c.sqrt_dietPC1, symp, c.Zage, chan_node, chan_lineage_1, chan_lineage_2, stringsAsFactors=F)
chan_diet1 <- na.omit(chan_diet1)
chan_diet2 <- data.frame(c.sqrt_dietPC2, symp, c.Zage, chan_node, chan_lineage_1, chan_lineage_2, stringsAsFactors=F)
chan_diet2 <- na.omit(chan_diet2)
chan_diet3 <- data.frame(c.sqrt_dietPC3, symp, c.Zage, chan_node, chan_lineage_1, chan_lineage_2, stringsAsFactors=F)
chan_diet3 <- na.omit(chan_diet3)
chan_size <- data.frame(c.sqrt_size, symp, c.Zage, chan_node, chan_lineage_1, chan_lineage_2, stringsAsFactors=F)
chan_size <- na.omit(chan_size)

# running the model
# set inverse-gamma prior for random effects
prior.fixed <- list(G=list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002)), R=list(V=1, nu=0.002))

# create phylogenetic covariance matrix
phylovcv <- inverseA(chan_tree, nodes="ALL")

# set up and fit the models
chann_shape1_model <- MCMCglmm(c.sqrt_shapePC1 ~ symp + c.Zage + symp*c.Zage, random = ~chan_node + chan_lineage_1 + chan_lineage_2, data=chan_shape1, ginverse = list(chan_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)
chann_shape2_model <- MCMCglmm(c.sqrt_shapePC2 ~ symp + c.Zage + symp*c.Zage, random = ~chan_node + chan_lineage_1 + chan_lineage_2, data=chan_shape2, ginverse = list(chan_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)
chann_shape3_model <- MCMCglmm(c.sqrt_shapePC3 ~ symp + c.Zage + symp*c.Zage, random = ~chan_node + chan_lineage_1 + chan_lineage_2, data=chan_shape3, ginverse = list(chan_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)

chann_diet1_model <- MCMCglmm(c.sqrt_dietPC1 ~ symp + c.Zage + symp*c.Zage, random = ~chan_node + chan_lineage_1 + chan_lineage_2, data=chan_diet1, ginverse = list(chan_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)
chann_diet2_model <- MCMCglmm(c.sqrt_dietPC2 ~ symp + c.Zage + symp*c.Zage, random = ~chan_node + chan_lineage_1 + chan_lineage_2, data=chan_diet2, ginverse = list(chan_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)
chann_diet3_model <- MCMCglmm(c.sqrt_dietPC3 ~ symp + c.Zage + symp*c.Zage, random = ~chan_node + chan_lineage_1 + chan_lineage_2, data=chan_diet3, ginverse = list(chan_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)

chann_size_model <- MCMCglmm(c.sqrt_size ~ symp + c.Zage + symp*c.Zage, random = ~chan_node + chan_lineage_1 + chan_lineage_2, data=chan_size, ginverse = list(chan_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)

# check for convergence of the model
# first, plot the model diagnostics
par(mar=c(1,1,1,1))
plot(chann_shape1_model$Sol)
plot(chann_shape1_model$VCV)

# next calculate effective samples sizes (ESS)
effectiveSize(chann_shape1_model$Sol[, 1:chann_shape1_model$Fixed$nfl, drop = FALSE])

# finally check levels of autocorrelation
autocorr(chann_shape1_model$VCV)

# now check the outputs of the model
summary(chann_shape1_model)

# now for all the other models

plot(chann_shape2_model$Sol)
plot(chann_shape2_model$VCV)
effectiveSize(chann_shape2_model$Sol[, 1:chann_shape2_model$Fixed$nfl, drop = FALSE])
autocorr(chann_shape2_model$VCV)
summary(chann_shape2_model)

plot(chann_shape3_model$Sol)
plot(chann_shape3_model$VCV)
effectiveSize(chann_shape3_model$Sol[, 1:chann_shape3_model$Fixed$nfl, drop = FALSE])
autocorr(chann_shape3_model$VCV)
summary(chann_shape3_model)

plot(chann_diet1_model$Sol)
plot(chann_diet1_model$VCV)
effectiveSize(chann_diet1_model$Sol[, 1:chann_diet1_model$Fixed$nfl, drop = FALSE])
autocorr(chann_diet1_model$VCV)
summary(chann_diet1_model)

plot(chann_diet2_model$Sol)
plot(chann_diet2_model$VCV)
effectiveSize(chann_diet2_model$Sol[, 1:chann_diet2_model$Fixed$nfl, drop = FALSE])
autocorr(chann_diet2_model$VCV)
summary(chann_diet2_model)

plot(chann_diet3_model$Sol)
plot(chann_diet3_model$VCV)
effectiveSize(chann_diet3_model$Sol[, 1:chann_diet3_model$Fixed$nfl, drop = FALSE])
autocorr(chann_diet3_model$VCV)
summary(chann_diet3_model)

plot(chann_size_model$Sol)
plot(chann_size_model$VCV)
effectiveSize(chann_size_model$Sol[, 1:chann_size_model$Fixed$nfl, drop = FALSE])
autocorr(chann_size_model$VCV)
summary(chann_size_model)

############################################
# Next complete analyses for notoperches
############################################

#square root (sqrt) transform response variables
t.sqrt_buo<-sqrt(trem_data$buo)
t.sqrt_depth<-sqrt(trem_data$depth)
t.sqrt_dietPC1<-sqrt(trem_data$dietPC1)
t.sqrt_dietPC2<-sqrt(trem_data$dietPC2)
t.sqrt_dietPC3<-sqrt(trem_data$dietPC3)
t.sqrt_shapePC1<-sqrt(trem_data$shapePC1)
t.sqrt_shapePC2<-sqrt(trem_data$shapePC2)
t.sqrt_shapePC3<-sqrt(trem_data$shapePC3)
t.sqrt_size<-sqrt(trem_data$size)

#sqrt- and z-transform covariate (age)
t.sqrt_age<-sqrt(trem_data$age)
t.Zage<-((t.sqrt_age - mean(t.sqrt_age))/(sd(t.sqrt_age)))

#create species identity random effects 
trem_sp <- unique(trem_data$lineage1)
trem_sp_id <- data.frame(1:length(trem_sp), trem_sp, stringsAsFactors=F)
names(trem_sp_id) <- c("id", "sp")

trem_lineage_1 <- apply(trem_data, 1, function(x) trem_sp_id$id[grep(x[1], trem_sp_id$sp)])
trem_lineage_2 <- apply(trem_data, 1, function(x) trem_sp_id$id[grep(x[2], trem_sp_id$sp)])

#prune tree so that tips match trait data and retain one outgroup taxon to root phylogeny
trem_tree <- keep.tip(noto.tree, c(as.vector(trem_sp), "7738_Dissostichus_m"))
trem_tree <- root(trem_tree, "7738_Dissostichus_m")
is.binary.tree(trem_tree)
is.ultrametric(trem_tree)
is.rooted(trem_tree)

#create node random effect from pruned tree
trem_node <- apply(trem_data, 1, function(x) getMRCA(trem_tree, c(x[1], x[2]))) 
trem_node <- trem_node-length(trem_tree$tip.label)
trem_node <- paste("Node", trem_node, sep="")

#create object with data on sympatry and allopatry
symp<-trem_data$sympatry

#add transformed variables, node IDs, and lineage IDs to single dataframes for each trait
trem_buo <- data.frame(t.sqrt_buo, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_buo <- na.omit(trem_buo)
trem_depth <- data.frame(t.sqrt_depth, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_depth <- na.omit(trem_depth)
trem_shape1 <- data.frame(t.sqrt_shapePC1, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_shape1 <- na.omit(trem_shape1)
trem_shape2 <- data.frame(t.sqrt_shapePC2, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_shape2 <- na.omit(trem_shape2)
trem_shape3 <- data.frame(t.sqrt_shapePC3, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_shape3 <- na.omit(trem_shape3)
trem_diet1 <- data.frame(t.sqrt_dietPC1, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_diet1 <- na.omit(trem_diet1)
trem_diet2 <- data.frame(t.sqrt_dietPC2, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_diet2 <- na.omit(trem_diet2)
trem_diet3 <- data.frame(t.sqrt_dietPC3, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_diet3 <- na.omit(trem_diet3)
trem_size <- data.frame(t.sqrt_size, symp, t.Zage, trem_node, trem_lineage_1, trem_lineage_2, stringsAsFactors=F)
trem_size <- na.omit(trem_size)

#running the model
#set inverse-gamma prior for random effects
prior.fixed <- list(G=list(G1=list(V=1, nu=0.002), G2=list(V=1, nu=0.002), G3=list(V=1, nu=0.002)), R=list(V=1, nu=0.002))

#create phylogenetic covariance matrix
phylovcv <- inverseA(trem_tree, nodes="ALL")

#set up and fit the model

trem_buo_model <- MCMCglmm(t.sqrt_buo ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_buo, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)

trem_depth_model <- MCMCglmm(t.sqrt_depth ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_depth, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)

trem_shape1_model <- MCMCglmm(t.sqrt_shapePC1 ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_shape1, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)
trem_shape2_model <- MCMCglmm(t.sqrt_shapePC2 ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_shape2, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)
trem_shape3_model <- MCMCglmm(t.sqrt_shapePC3 ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_shape3, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)

trem_diet1_model <- MCMCglmm(t.sqrt_dietPC1 ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_diet1, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)
trem_diet2_model <- MCMCglmm(t.sqrt_dietPC2 ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_diet2, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)
trem_diet3_model <- MCMCglmm(t.sqrt_dietPC3 ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_diet3, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)

trem_size_model <- MCMCglmm(t.sqrt_size ~ symp + t.Zage + symp*t.Zage, random = ~trem_node + trem_lineage_1 + trem_lineage_2, data=trem_size, ginverse = list(trem_node = phylovcv$Ainv), prior=prior.fixed, nitt=20000000, burnin=2000000, thin=10000, pl=TRUE, pr=TRUE)

#check for convergence and check model outputs
par(mar=c(1,1,1,1))
plot(trem_buo_model$Sol)
plot(trem_buo_model$VCV)
effectiveSize(trem_buo_model$Sol[, 1:trem_buo_model$Fixed$nfl, drop = FALSE])
autocorr(trem_buo_model$VCV)
summary(trem_buo_model)

plot(trem_depth_model$Sol)
plot(trem_depth_model$VCV)
effectiveSize(trem_depth_model$Sol[, 1:trem_depth_model$Fixed$nfl, drop = FALSE])
autocorr(trem_depth_model$VCV)
summary(trem_depth_model)

plot(trem_shape1_model$Sol)
plot(trem_shape1_model$VCV)
effectiveSize(trem_shape1_model$Sol[, 1:trem_shape1_model$Fixed$nfl, drop = FALSE])
autocorr(trem_shape1_model$VCV)
summary(trem_shape1_model)

plot(trem_shape2_model$Sol)
plot(trem_shape2_model$VCV)
effectiveSize(trem_shape2_model$Sol[, 1:trem_shape2_model$Fixed$nfl, drop = FALSE])
autocorr(trem_shape2_model$VCV)
summary(trem_shape2_model)

plot(trem_shape3_model$Sol)
plot(trem_shape3_model$VCV)
effectiveSize(trem_shape3_model$Sol[, 1:trem_shape3_model$Fixed$nfl, drop = FALSE])
autocorr(trem_shape3_model$VCV)
summary(trem_shape3_model)

plot(trem_diet1_model$Sol)
plot(trem_diet1_model$VCV)
effectiveSize(trem_diet1_model$Sol[, 1:trem_diet1_model$Fixed$nfl, drop = FALSE])
autocorr(trem_diet1_model$VCV)
summary(trem_diet1_model)

plot(trem_diet2_model$Sol)
plot(trem_diet2_model$VCV)
effectiveSize(trem_diet2_model$Sol[, 1:trem_diet2_model$Fixed$nfl, drop = FALSE])
autocorr(trem_diet2_model$VCV)
summary(trem_diet2_model)

plot(trem_diet3_model$Sol)
plot(trem_diet3_model$VCV)
effectiveSize(trem_diet3_model$Sol[, 1:trem_diet3_model$Fixed$nfl, drop = FALSE])
autocorr(trem_diet3_model$VCV)
summary(trem_diet3_model)

plot(trem_size_model$Sol)
plot(trem_size_model$VCV)
effectiveSize(trem_size_model$Sol[, 1:trem_size_model$Fixed$nfl, drop = FALSE])
autocorr(trem_size_model$VCV)
summary(trem_size_model)
