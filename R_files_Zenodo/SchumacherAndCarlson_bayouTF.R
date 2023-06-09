#############################################################################################
#	Brain Mass Analysis for Schumacher and Carlson 2022										#
#																							#
#                Combines new otophysan and osteoglossiform Sukhum et al 2018 data          #
#				 with published actinopterygian data from Tsuboi et al 2018					#
#				 and Tsuboi 2021															#
#																							#
#				 TF Model																	#
#																							#
#				 Requires: SchumacherAndCarlson_SupplementaryFile2.csv						#
#             			   actinopt_12k_treePL.tre											#
#				 Revised May 2022															#
#																							#
#############################################################################################

## R script to run multiple bayou models on the cluster

setwd("/storage1/fs1/carlson.bruce/Active/ELS/bayou");

library(plyr)
library(tidyverse)
library(car)
library(phytools)
library(ape)
library(nlme)
library(bayou)
library(ggplot2)


############################
#    phylogeny and data    #
############################

## actinopterygii phylogeny from Rabosky et al 2018
actino_rab_tree<-read.tree("actinopt_12k_treePL.tre")

## adjustments to phylo tip labels
## correct genus for B niger,
idx<-grep("Brienomyrus_niger",actino_rab_tree$tip.label)
actino_rab_tree$tip.label[idx]<-"Brevimyrus_niger"
## use the krypto with the shortest branch length from the genus node as k vitreolus
idx<-grep("Kryptopterus_limpok",actino_rab_tree$tip.label)
actino_rab_tree$tip.label[idx]<-"Kryptopterus_vitreolus"

## reading in data and subsetting tree/data to only incl spp in both
tel_data <- read.csv("SchumacherAndCarlson_SupplementaryFile2.csv")
tel_data[,11:12] <- log10(tel_data[,11:12])
tel_data<-subset(tel_data, Genus_Species != "Carassius_carassius")					## removing C carassius following Tsuboi 21

## used rng to choose which spp of spp pairs containng very short terminal branches to drop
tel_data<-subset(tel_data, Genus_Species != "Aphareus_furca")					
tel_data<-subset(tel_data, Genus_Species != "Bodianus_perditio")
tel_data<-subset(tel_data, Genus_Species != "Tanganicodus_irsacae")
tel_data<-subset(tel_data, Genus_Species != "Symphodus_melanocercus")
tel_data<-subset(tel_data, Genus_Species != "Chaetodon_mertensii")
tel_data<-subset(tel_data, Genus_Species != "Anoplogaster_cornuta")
tel_data<-subset(tel_data, Genus_Species != "Sufflamen_albicaudatum")
tel_data<-subset(tel_data, Genus_Species != "Gaidropsarus_ensis")					
tel_data<-subset(tel_data, Genus_Species != "Alosa_alosa")
tel_data<-subset(tel_data, Genus_Species != "Alosa_immaculata")


tel_tree_rab<-drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, tel_data$Genus_Species))  #tree with all study species

phylo_tel_data <- subset(tel_data, Genus_Species %in% tel_tree_rab$tip.label)		


###################################
#   bayesian allometry modeling   #
###################################

## phylo_tel_data has all the allometric data
## tel_tree_rab is the phylogeny

## reordering the tree to be in postorder which makes it easier to match with bayou output
tel_tree_rab <- reorder(tel_tree_rab, "postorder")

## calculating mean and standard error across species
## if N=1, SE=mean(SE)
tel_brain <- ddply(phylo_tel_data, "Genus_Species", summarise,
			   N    = length(Brain.weight..g.),
               mean = mean(Brain.weight..g.),
               sd   = sd(Brain.weight..g.),
               se   = sd / sqrt(N))
if (any(is.na(tel_brain$se))) {
            tel_brain$se[which(is.na(tel_brain$se))] <- mean(tel_brain$se, na.rm = TRUE)}
tel_body <- ddply(phylo_tel_data, "Genus_Species", summarise,
			   N    = length(Body.weight..g.),
               mean = mean(Body.weight..g.),
               sd   = sd(Body.weight..g.),
               se   = sd / sqrt(N))
if (any(is.na(tel_body$se))) {
            tel_body$se[which(is.na(tel_body$se))] <- mean(tel_body$se, na.rm = TRUE)}
tel_avg_data<-data.frame(setNames(tel_body$mean,tel_body$Genus_Species),setNames(tel_brain$mean,tel_brain$Genus_Species))
colnames(tel_avg_data)<-c("BodyMass","BrainMass")
#rownames(tel_avg_data)<-tel_body$Genus_Species
tel_se_data<-data.frame(setNames(tel_body$se,tel_body$Genus_Species),setNames(tel_brain$se,tel_brain$Genus_Species))
colnames(tel_se_data)<-c("BodyMass","BrainMass")
#rownames(tel_se_data)<-tel_body$Genus_Species

tel_brain_se <- tel_brain$se
names(tel_brain_se) <- tel_brain$Genus_Species

## Ordering data to match tip labels.
tel_avg_data <- tel_avg_data[match(tel_tree_rab$tip.label,rownames(tel_avg_data)),]
tel_se_data <- tel_se_data[match(tel_tree_rab$tip.label,rownames(tel_se_data)),]

BodyMass <- setNames(tel_avg_data$BodyMass, rownames(tel_avg_data))
BrainMass <- setNames(tel_avg_data$BrainMass, rownames(tel_avg_data))
BrainSE <- setNames(tel_se_data$BrainMass, rownames(tel_se_data))

## 37 orders total w data, 13 with >= 9 spp
## calculate pgls for orders (N=8) with spp >= 20
#anguillid_names <- unique(subset(phylo_tel_data, Order == 'Anguilliformes')$Genus_Species)
#berycid_names <- unique(subset(phylo_tel_data, Order == 'Beryciformes')$Genus_Species)
#gadid_names <- unique(subset(phylo_tel_data, Order == 'Gadiformes')$Genus_Species)
#osteoglossid_names <- unique(subset(phylo_tel_data, Order == 'Osteoglossiformes')$Genus_Species)
#percid_names <- unique(subset(phylo_tel_data, Order == 'Perciformes')$Genus_Species)
#scorpaenid_names <- unique(subset(phylo_tel_data, Order == 'Scorpaeniformes')$Genus_Species)
#syngnathid_names <- unique(subset(phylo_tel_data, Order == 'Syngnathiformes')$Genus_Species)
#tetraodontid_names <- unique(subset(phylo_tel_data, Order == 'Tetraodontiformes')$Genus_Species)
## orders (N=4) with spp >= 10
#belonid_names <- unique(subset(phylo_tel_data, Order == 'Beloniformes')$Genus_Species)
#gymnotid_names <- unique(subset(phylo_tel_data, Order == 'Gymnotiformes')$Genus_Species)
#silurid_names <- unique(subset(phylo_tel_data, Order == 'Siluriformes')$Genus_Species)
#stomiid_names <- unique(subset(phylo_tel_data, Order == 'Stomiiformes')$Genus_Species)
## spp = 9
#clupeid_names <- unique(subset(phylo_tel_data, Order == 'Clupeiformes')$Genus_Species)

#anguillid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% anguillid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, anguillid_names))))
#berycid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% berycid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, berycid_names))))
#gadid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% gadid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, gadid_names))))
#osteoglossid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% osteoglossid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, osteoglossid_names))))
#percid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% percid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, percid_names))))
#scorpaenid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% scorpaenid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, scorpaenid_names))))
#syngnathid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% syngnathid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, syngnathid_names))))
#tetraodontid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% tetraodontid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, tetraodontid_names))))
#alltel_pgls<-gls(BrainMass~BodyMass, data=tel_avg_data, correlation=corBrownian(phy=tel_tree_rab))
#belonid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% belonid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, belonid_names))))
#gymnotid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% gymnotid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, gymnotid_names))))
#silurid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% silurid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, silurid_names))))
#stomiid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% stomiid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, stomiid_names))))
#clupeid_pgls<-gls(BrainMass~BodyMass, data=subset(tel_avg_data, rownames(tel_avg_data) %in% clupeid_names), correlation=corBrownian(phy=drop.tip(actino_rab_tree, setdiff(actino_rab_tree$tip.label, clupeid_names))))

#ord_ints <- c(anguillid_pgls$coefficients[1],berycid_pgls$coefficients[1],gadid_pgls$coefficients[1],osteoglossid_pgls$coefficients[1],percid_pgls$coefficients[1],scorpaenid_pgls$coefficients[1],syngnathid_pgls$coefficients[1],tetraodontid_pgls$coefficients[1],alltel_pgls$coefficients[1],belonid_pgls$coefficients[1],gymnotid_pgls$coefficients[1],silurid_pgls$coefficients[1],stomiid_pgls$coefficients[1],clupeid_pgls$coefficients[1])
#ord_ints_mean<-mean(ord_ints)		## = -1.653
#ord_ints_sd<-sd(ord_ints)			## = .306 <- raise to .4
## Converting Tsuboi intercepts from ln to log10 (log10(e^intercept#)) -> osetoglossids = -1.771, anguillids = -2.118, beryciforms = -1.494, gadiformes = -1.664, perciformes = -1.750, scorpaniformes = -1.711, sygnathiformes = -2.087, tetraodontiforms = -1.647
## Mean tsuboi intercepts = -1.78025, sd = 0.2161929;

#ord_slps <- c(anguillid_pgls$coefficients[2],berycid_pgls$coefficients[2],gadid_pgls$coefficients[2],osteoglossid_pgls$coefficients[2],percid_pgls$coefficients[2],scorpaenid_pgls$coefficients[2],syngnathid_pgls$coefficients[2],tetraodontid_pgls$coefficients[2],alltel_pgls$coefficients[2],belonid_pgls$coefficients[2],gymnotid_pgls$coefficients[2],silurid_pgls$coefficients[2],stomiid_pgls$coefficients[2],clupeid_pgls$coefficients[2])
#ord_slps_mean<-mean(ord_slps)		## = .505
#ord_slps_sd<-sd(ord_slps)			## = .114


## "This is a good reminder to *always try to use measurement error in your analyses*. OU models especially are affected by measurement error. 
## This is because OU models have the effect of "erasing" evolutionary history with increasing *alpha*. If you don't account for measurement error, 
## then that measurement error will be transferred to the evolutionary process model. You can make a Brownian Motion model look very OU like if there is a lot of measurement error."

## Defining the priors: "prior function is going to take our parameters and output the *prior probability* of our parameter values. It represents our initial degree of belief in what values the parameters will take."
## "*make.prior* tries to be reasonable and not make you type everything out, but **do not be comfortable with defaults**. One trick to make sure your prior functions are reasonable is to simulate a bunch of values 
## and take the quantiles of the distribution. We are using a half-Cauchy distribution for *alpha* and *sig2*, which is a good weakly informative prior for scale parameters."

## theta = mean trait value estimate
## alpha = rate at which theta changes
## sigma = variance in trait value over time
## total # shifts = conditional poisson prior; mean = 1% total branches, max = 5% total branches

max_shft <- nrow(tel_avg_data)*.05
avg_shft <- nrow(tel_avg_data)*.01			## tried for runs 300s as Smaers had lambda = 10 and max = 100

## Run models w 1. global intercept and slopes (GG), 2. separate intercept and global slope (SG), and 3. separate intercepts and slopes (SS)
## And use model comparison to see which has a better fit (likely the third model which allows both slope and intercept to vary)
## Defined priors following Uyeda et al 2017

## Using .599 as slope prior mean and .2 as sd; .599 = mean of 8 reported tsuboi slopes, sd of tsuboi slope = .18 -> raised to .2
## Using -1.653 as theta prior mean and .4 as sd; mean and sd of 13 orders (spp>=9) and fit across all tel

########################
#   parameter tuning   #
########################

## Set tuning parameters:
## "There is a bit of art to tuning the parameters, which may require making multiple runs and trying to get the acceptance ratios in the right 
## region (0.2-0.4)...If the acceptance ratio for a certain parameter is too high, increase the tuning parameter for that variable. If the 
## acceptance ratio is too low, decrease it. The scale of the regression coefficient, for example, should give you some idea of what these parameters should be."

DSS.47 = list(alpha=7, sig2=.7, beta_BodyMass=0.5, k=1, theta=3, slide=1)


##########################
#   fixed shift chains   #
##########################

## First, make cache object to pull info about brnch nums for bayou sampling function
tree <- tel_tree_rab
cache <- bayou:::.prepare.ou.univariate(tree, setNames(tel_avg_data$BrainMass, tree$tip.label), SE=BrainSE)
pred <- as.matrix(BodyMass)
colnames(pred) <- "BodyMass"

## SS Model, Hypothesis Tuberous
## k = 2 shifts
## mormyroid sb = 97
## gymnotes sb = 164

k <- 2
sb <- c(97, 164)

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(351)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(351)
model.fixed.SS.351 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.351$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.351 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.351$model, prior=prior.fixed.SS, startpar=model.fixed.SS.351$startpar, file.dir=getwd(), outname="TUBfixedSS_r351", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.351$run(2000000)
chain.SS.351 <- fixed.mcmc.SS.351$load()


set.seed(352)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(352)
model.fixed.SS.352 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.352$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.352 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.352$model, prior=prior.fixed.SS, startpar=model.fixed.SS.352$startpar, file.dir=getwd(), outname="TUBfixedSS_r352", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.352$run(2000000)
chain.SS.352 <- fixed.mcmc.SS.352$load()


set.seed(353)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(353)
model.fixed.SS.353 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.353$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.353 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.353$model, prior=prior.fixed.SS, startpar=model.fixed.SS.353$startpar, file.dir=getwd(), outname="TUBfixedSS_r353", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.353$run(2000000)
chain.SS.353 <- fixed.mcmc.SS.353$load()


set.seed(354)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(354)
model.fixed.SS.354 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.354$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.354 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.354$model, prior=prior.fixed.SS, startpar=model.fixed.SS.354$startpar, file.dir=getwd(), outname="TUBfixedSS_r354", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.354$run(2000000)
chain.SS.354 <- fixed.mcmc.SS.354$load()


set.seed(355)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(355)
model.fixed.SS.355 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.355$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.355 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.355$model, prior=prior.fixed.SS, startpar=model.fixed.SS.355$startpar, file.dir=getwd(), outname="TUBfixedSS_r355", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.355$run(2000000)
chain.SS.355 <- fixed.mcmc.SS.355$load()


set.seed(356)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(356)
model.fixed.SS.356 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.356$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.356 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.356$model, prior=prior.fixed.SS, startpar=model.fixed.SS.356$startpar, file.dir=getwd(), outname="TUBfixedSS_r356", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.356$run(2000000)
chain.SS.356 <- fixed.mcmc.SS.356$load()


set.seed(357)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(357)
model.fixed.SS.357 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.357$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.357 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.357$model, prior=prior.fixed.SS, startpar=model.fixed.SS.357$startpar, file.dir=getwd(), outname="TUBfixedSS_r357", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.357$run(2000000)
chain.SS.357 <- fixed.mcmc.SS.357$load()


set.seed(358)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(358)
model.fixed.SS.358 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.358$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.358 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.358$model, prior=prior.fixed.SS, startpar=model.fixed.SS.358$startpar, file.dir=getwd(), outname="TUBfixedSS_r358", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.358$run(2000000)
chain.SS.358 <- fixed.mcmc.SS.358$load()


set.seed(359)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(359)
model.fixed.SS.359 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.359$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.359 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.359$model, prior=prior.fixed.SS, startpar=model.fixed.SS.359$startpar, file.dir=getwd(), outname="TUBfixedSS_r359", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.359$run(2000000)
chain.SS.359 <- fixed.mcmc.SS.359$load()


set.seed(360)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.fixed.SS <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="fixed", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb="fixed",
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, sb=startpar$sb, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(360)
model.fixed.SS.360 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.47)
model.fixed.SS.360$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.360 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.360$model, prior=prior.fixed.SS, startpar=model.fixed.SS.360$startpar, file.dir=getwd(), outname="TUBfixedSS_r360", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.360$run(2000000)
chain.SS.360 <- fixed.mcmc.SS.360$load()


#####################
#   chain results   #
#####################

#chain.SS.351 <- set.burnin(fixed.mcmc.SS.351$load(), 0.3)
#out.SS.351 <- summary(chain.SS.351)
#chain.SS.352 <- set.burnin(fixed.mcmc.SS.352$load(), 0.3)
#out.SS.352 <- summary(chain.SS.352)
#chain.SS.353 <- set.burnin(fixed.mcmc.SS.353$load(), 0.3)
#out.SS.353 <- summary(chain.SS.353)
#chain.SS.354 <- set.burnin(fixed.mcmc.SS.354$load(), 0.3)
#out.SS.354 <- summary(chain.SS.354)
#chain.SS.355 <- set.burnin(fixed.mcmc.SS.355$load(), 0.3)
#out.SS.355 <- summary(chain.SS.355)
#chain.SS.356 <- set.burnin(fixed.mcmc.SS.356$load(), 0.3)
#out.SS.356 <- summary(chain.SS.356)
#chain.SS.357 <- set.burnin(fixed.mcmc.SS.357$load(), 0.3)
#out.SS.357 <- summary(chain.SS.357)
#chain.SS.358 <- set.burnin(fixed.mcmc.SS.358$load(), 0.3)
#out.SS.358 <- summary(chain.SS.358)
#chain.SS.359 <- set.burnin(fixed.mcmc.SS.359$load(), 0.3)
#out.SS.359 <- summary(chain.SS.359)
#chain.SS.360 <- set.burnin(fixed.mcmc.SS.360$load(), 0.3)
#out.SS.360 <- summary(chain.SS.360)

#pdf("TUBfixedSS.350s_trace_plots.pdf") 
# Creating a plot
#plot(chain.SS.351)
#plot(chain.SS.352)
#plot(chain.SS.353)
#plot(chain.SS.354)
#plot(chain.SS.355)
#plot(chain.SS.356)
#plot(chain.SS.357)
#plot(chain.SS.358)
#plot(chain.SS.359)
#plot(chain.SS.360)
# Closing the graphical device
#dev.off() 


chain.SS.351 <- fixed.mcmc.SS.351$load()
chain.SS.352 <- fixed.mcmc.SS.352$load()
chain.SS.353 <- fixed.mcmc.SS.353$load()
chain.SS.354 <- fixed.mcmc.SS.354$load()
chain.SS.355 <- fixed.mcmc.SS.355$load()
chain.SS.356 <- fixed.mcmc.SS.356$load()
chain.SS.357 <- fixed.mcmc.SS.357$load()
chain.SS.358 <- fixed.mcmc.SS.358$load()
chain.SS.359 <- fixed.mcmc.SS.359$load()
chain.SS.360 <- fixed.mcmc.SS.360$load()


## All 10 chains combined
SS.chains <- list(chain.SS.351,chain.SS.352,chain.SS.353,chain.SS.354,chain.SS.355,chain.SS.356,chain.SS.357,chain.SS.358,chain.SS.359,chain.SS.360)
SS.350s.comb <- combine.chains(SS.chains, burnin.prop = 0.3)

print("TUBfixedSS.350s.comb")
#summary(SS.350s.comb)
#pdf("TUBfixedSS.350s_comb_trace_plot.pdf")
#plot(SS.350s.comb)
#dev.off()

## summary output of comb chain:
#bayou MCMC chain: 2000001 generations
#140059 samples, first 0 samples discarded as burnin
#
#Summary statistics for parameters:
#                            Mean           SD     Naive SE Time-series SE
#lnL                 4.640454e+02 1.649963e+00 4.408781e-03   5.901977e-03
#prior               4.023470e+00 1.069918e+00 2.858874e-03   2.878100e-03
#alpha               4.993556e-04 4.293778e-04 1.147319e-06   1.404520e-06
#sig2                4.397999e-04 2.877593e-05 7.689070e-08   1.084567e-07
#k                   2.000000e+00 0.000000e+00 0.000000e+00   0.000000e+00
#ntheta              3.000000e+00 0.000000e+00 0.000000e+00   0.000000e+00
#root.theta         -1.821096e+00 1.634007e-01 4.366149e-04   4.609388e-04
#root.beta_BodyMass  5.096440e-01 7.348635e-03 1.963592e-05   3.154702e-05
#all theta          -1.708506e+00 3.474963e-01           NA             NA
#all beta_BodyMass   6.241526e-01 1.219132e-01           NA             NA
#                   Effective Size    HPD95Lower    HPD95Upper
#lnL                      78154.34  4.608077e+02  4.666510e+02
#prior                   138193.98  1.891352e+00  5.523159e+00
#alpha                    93459.51  6.110418e-09  1.349714e-03
#sig2                     70395.71  3.847510e-04  4.966082e-04
#k                            0.00  2.000000e+00  2.000000e+00
#ntheta                       0.00  3.000000e+00  3.000000e+00
#root.theta              125667.10 -2.145804e+00 -1.501791e+00
#root.beta_BodyMass       54262.11  4.952916e-01  5.240525e-01
#all theta                      NA            NA            NA
#all beta_BodyMass              NA            NA            NA
#
#Branches with posterior probabilities higher than 0.1:
#    pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
#97   1           -1.659071        0.001065443    0.5016795
#164  1           -1.645352        0.001060314    0.5000418


## Adding stepping stone to estimate marginal likelihoods and compare bayes factors for model comparison
## Reference code from Uyeda et al
require(foreach)
require(doParallel)
registerDoParallel(cores=10)

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

currentTime <- Sys.time()
print(currentTime)

ss.SS.350s<-fixed.mcmc.SS.351$steppingstone(500000, SS.350s.comb, Bk = Bk, burnin=0.3, plot=FALSE)
ss.SS.350s$lnr
## lnr == 446.9542
pdf("TUBfixedSS.350s.comb_ss_plot.pdf")
plot(ss.SS.350s)
dev.off()
print("TUBfixedss.SS.350s")
ss.SS.350s$lnr
currentTime <- Sys.time()
print(currentTime)

