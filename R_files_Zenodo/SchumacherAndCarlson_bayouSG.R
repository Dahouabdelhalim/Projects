#############################################################################################
#	Brain Mass Analysis for Schumacher and Carlson 2022										#
#																							#
#                Combines new otophysan and osteoglossiform Sukhum et al 2018 data          #
#				 with published actinopterygian data from Tsuboi et al 2018					#
#				 and Tsuboi 2021															#
#																							#
#				 SG Model																	#
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

DSG.47 = list(alpha=7.5, sig2=.7, beta_BodyMass=0.08, k=1, theta=4, slide=1)


#################
#   SG models   #
#################

## First initialize the MCMC chain with 10k gens of no birth death proposals
## i.e. run a separate MCMC chain first with fixed shift locations (randomly drawn from priors)
## then use the output of that chain as the starting values for the rjMCMC <- allows other params to optimize briefly before adding rj to improve fit
## Primer section based off analysis code from Smaers et al 21 -- Thanks for sharing JS!

## First, make cache object to pull info about brnch nums for bayou sampling function
tree <- tel_tree_rab
cache <- bayou:::.prepare.ou.univariate(tree, setNames(tel_avg_data$BrainMass, tree$tip.label), SE=BrainSE)
pred <- as.matrix(BodyMass)
colnames(pred) <- "BodyMass"

prior.SG <- make.prior(tel_tree_rab, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                       			  dbeta_BodyMass=list(mean=0.599, sd=0.2),
                                  dsb=list(bmax=1, prob=1),
                                  dk=list(lambda=avg_shft, kmax=max_shft),
                                  dtheta=list(mean=-1.653, sd=.4))
)


#######
# 471 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(471)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(471)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(471)
model.primer.SG.471 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.471$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.471 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.471$model, prior=prior.primer.SG, startpar=model.primer.SG.471$startpar, file.dir=getwd(), outname="primerSG_r471", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.471$run(1000000)
primer.chain.SG <- primer.mcmc.SG.471$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.471$model)


## Building the rjMCMC models
set.seed(471)
model.SG.471 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.471 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.471$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r471", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(471)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.471$run(2000000)			


#######
# 472 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(472)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(472)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(472)
model.primer.SG.472 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.472$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.472 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.472$model, prior=prior.primer.SG, startpar=model.primer.SG.472$startpar, file.dir=getwd(), outname="primerSG_r472", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.472$run(1000000)
primer.chain.SG <- primer.mcmc.SG.472$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.472$model)



## Building the rjMCMC models
set.seed(472)
model.SG.472 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.472 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.472$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r472", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(472)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.472$run(2000000)			


#######
# 473 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(473)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(473)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(473)
model.primer.SG.473 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.473$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.473 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.473$model, prior=prior.primer.SG, startpar=model.primer.SG.473$startpar, file.dir=getwd(), outname="primerSG_r473", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.473$run(1000000)
primer.chain.SG <- primer.mcmc.SG.473$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.473$model)



## Building the rjMCMC models
set.seed(473)
model.SG.473 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.473 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.473$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r473", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(473)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.473$run(2000000)			


#######
# 474 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(474)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(474)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(474)
model.primer.SG.474 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.474$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.474 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.474$model, prior=prior.primer.SG, startpar=model.primer.SG.474$startpar, file.dir=getwd(), outname="primerSG_r474", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.474$run(1000000)
primer.chain.SG <- primer.mcmc.SG.474$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.474$model)


## Building the rjMCMC models
set.seed(474)
model.SG.474 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.474 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.474$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r474", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(474)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.474$run(2000000)			


#######
# 475 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(475)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(475)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(475)
model.primer.SG.475 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.475$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.475 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.475$model, prior=prior.primer.SG, startpar=model.primer.SG.475$startpar, file.dir=getwd(), outname="primerSG_r475", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.475$run(1000000)
primer.chain.SG <- primer.mcmc.SG.475$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.475$model)


## Building the rjMCMC models
set.seed(475)
model.SG.475 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.475 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.475$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r475", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(475)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.475$run(2000000)			


#######
# 476 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(476)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(476)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(476)
model.primer.SG.476 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.476$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.476 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.476$model, prior=prior.primer.SG, startpar=model.primer.SG.476$startpar, file.dir=getwd(), outname="primerSG_r476", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.476$run(1000000)
primer.chain.SG <- primer.mcmc.SG.476$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.476$model)


## Building the rjMCMC models
set.seed(476)
model.SG.476 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.476 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.476$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r476", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(476)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.476$run(2000000)			


#######
# 477 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(477)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(477)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(477)
model.primer.SG.477 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.477$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.477 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.477$model, prior=prior.primer.SG, startpar=model.primer.SG.477$startpar, file.dir=getwd(), outname="primerSG_r477", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.477$run(1000000)
primer.chain.SG <- primer.mcmc.SG.477$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.477$model)


## Building the rjMCMC models
set.seed(477)
model.SG.477 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.477 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.477$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r477", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(477)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.477$run(2000000)			


#######
# 478 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(478)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(478)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(478)
model.primer.SG.478 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.478$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.478 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.478$model, prior=prior.primer.SG, startpar=model.primer.SG.478$startpar, file.dir=getwd(), outname="primerSG_r478", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.478$run(1000000)
primer.chain.SG <- primer.mcmc.SG.478$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.478$model)


## Building the rjMCMC models
set.seed(478)
model.SG.478 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.478 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.478$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r478", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(478)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.478$run(2000000)			


#######
# 479 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(479)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(479)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(479)
model.primer.SG.479 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.479$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.479 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.479$model, prior=prior.primer.SG, startpar=model.primer.SG.479$startpar, file.dir=getwd(), outname="primerSG_r479", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.479$run(1000000)
primer.chain.SG <- primer.mcmc.SG.479$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.479$model)


## Building the rjMCMC models
set.seed(479)
model.SG.479 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.479 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.479$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r479", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(479)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.479$run(2000000)			


#######
# 480 #
#######

## SG Model Primer
## Initializing with a smaller shft num than 2.5% total brnchs to help improve fit
## Randomly draw a number of shifts from the prior distribution
set.seed(480)
k <-rpois(1, lambda=5)

## Next, "sample K # of branches randomly drawn from the phylogeny with a probability proportional to the size of the clade descended from that branch"
sb <- NULL
sb <- sample(((1:length(cache$bdesc))), k, replace=FALSE, prob = sapply(cache$bdesc, length))

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(480)
startpar <- list(alpha=exp(rnorm(1, -3, 0.5)), 
                 sig2=rnorm(1,0.75, 0.1), 
                 k=k, 
                 ntheta=k+1, 
                 beta_BodyMass=rnorm(1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SG <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dk="fixed", dsb="dsb", dtheta="dnorm", dloc="dloc"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=.599, sd=.2),
                                  dk="fixed", 
                                  dsb=list(bmax=1, prob=1),
                                  dtheta=list(mean=-1.653, sd=.4)),
                       fixed=list(ntheta=startpar$ntheta, k=startpar$k, t2=startpar$t2)
                       )
                       
## Building the primer model
set.seed(480)
model.primer.SG.480 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SG, startpar=startpar, D=DSG.47)
model.primer.SG.480$model$moves$k<-".splitmergebd"
prior.primer.SG(startpar)

chainN <- "Primer"

primer.mcmc.SG.480 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SG.480$model, prior=prior.primer.SG, startpar=model.primer.SG.480$startpar, file.dir=getwd(), outname="primerSG_r480", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SG.480$run(1000000)
primer.chain.SG <- primer.mcmc.SG.480$load()

startpar.SG <- pull.pars(11, primer.chain.SG, model=model.primer.SG.480$model)


## Building the rjMCMC models
set.seed(480)
model.SG.480 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SG, startpar=startpar.SG, D=DSG.47)

## Making the MCMC objects
mcmc.SG.480 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SG.480$model, prior=prior.SG, startpar=startpar.SG, file.dir=getwd(), outname="modelSG_r480", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(480)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.SG.480$run(2000000)			


#####################
#   Chain results   #
#####################

## Chain results
#chain.SG.471 <- set.burnin(mcmc.SG.471$load(), 0.3)
#chain.SG.472 <- set.burnin(mcmc.SG.472$load(), 0.3)
#chain.SG.473 <- set.burnin(mcmc.SG.473$load(), 0.3)
#chain.SG.474 <- set.burnin(mcmc.SG.474$load(), 0.3)
#chain.SG.475 <- set.burnin(mcmc.SG.475$load(), 0.3)
#chain.SG.476 <- set.burnin(mcmc.SG.476$load(), 0.3)
#chain.SG.477 <- set.burnin(mcmc.SG.477$load(), 0.3)
#chain.SG.478 <- set.burnin(mcmc.SG.478$load(), 0.3)
#chain.SG.479 <- set.burnin(mcmc.SG.479$load(), 0.3)
#chain.SG.480 <- set.burnin(mcmc.SG.480$load(), 0.3)

chain.SG.471 <- mcmc.SG.471$load()
chain.SG.472 <- mcmc.SG.472$load()
chain.SG.473 <- mcmc.SG.473$load()
chain.SG.474 <- mcmc.SG.474$load()
chain.SG.475 <- mcmc.SG.475$load()
chain.SG.476 <- mcmc.SG.476$load()
chain.SG.477 <- mcmc.SG.477$load()
chain.SG.478 <- mcmc.SG.478$load()
chain.SG.479 <- mcmc.SG.479$load()
chain.SG.480 <- mcmc.SG.480$load()


#out.SG.471 <- summary(chain.SG.471)
#out.SG.472 <- summary(chain.SG.472)
#out.SG.473 <- summary(chain.SG.473)
#out.SG.474 <- summary(chain.SG.474)
#out.SG.475 <- summary(chain.SG.475)
#out.SG.476 <- summary(chain.SG.476)
#out.SG.477 <- summary(chain.SG.477)
#out.SG.478 <- summary(chain.SG.478)
#out.SG.479 <- summary(chain.SG.479)
#out.SG.480 <- summary(chain.SG.480)

#pdf("SG.470s_trace_plots.pdf") 
# Creating a plot
#plot(chain.SG.471)
#plot(chain.SG.472)
#plot(chain.SG.473)
#plot(chain.SG.474)
#plot(chain.SG.475)
#plot(chain.SG.476)
#plot(chain.SG.477)
#plot(chain.SG.478)
#plot(chain.SG.479)
#plot(chain.SG.480)
# Closing the graphical device
#dev.off() 


## All 10 chains combined
SG.chains <- list(chain.SG.471,chain.SG.472,chain.SG.473,chain.SG.474,chain.SG.475,chain.SG.476,chain.SG.477,chain.SG.478,chain.SG.479,chain.SG.480)
SG.470s.comb <- combine.chains(SG.chains, burnin.prop = 0.3)

print("SG.470s.comb")
#summary(SG.470s.comb)
#pdf("SG.470s_comb_trace_plot.pdf")
#plot(SG.470s.comb)
#dev.off()

## summary output of comb chain:
#bayou MCMC chain: 2000001 generations
#140047 samples, first 0 samples discarded as burnin
#
#
#Summary statistics for parameters:
#                       Mean           SD     Naive SE Time-series SE
#lnL            4.579521e+02 1.496372e+00 3.998553e-03   1.176613e-02
#prior         -2.887068e+01 1.347535e+01 3.600834e-02   4.121131e-01
#alpha          5.146245e-04 4.466227e-04 1.193449e-06   2.139235e-06
#sig2           4.447372e-04 2.941369e-05 7.859821e-08   1.012359e-07
#beta_BodyMass  5.122124e-01 7.349651e-03 1.963947e-05   2.012490e-05
#k              4.411176e+00 2.060699e+00 5.506527e-03   6.352324e-02
#ntheta         5.411176e+00 2.060699e+00 5.506527e-03   6.352324e-02
#root.theta    -1.820398e+00 1.642240e-01 4.388335e-04   5.900117e-04
#all theta     -1.683557e+00 3.743247e-01           NA             NA
#              Effective Size    HPD95Lower    HPD95Upper
#lnL                16173.806  4.551346e+02  4.600113e+02
#prior               1069.171 -5.495832e+01 -6.014970e+00
#alpha              43587.766  2.967294e-09  1.391360e-03
#sig2               84417.072  3.881844e-04  5.027796e-04
#beta_BodyMass     133372.386  4.979692e-01  5.268124e-01
#k                   1052.359  1.000000e+00  8.000000e+00
#ntheta              1052.359  2.000000e+00  9.000000e+00
#root.theta         77473.326 -2.144484e+00 -1.497093e+00
#all theta                 NA            NA            NA
#
#
#Branches with posterior probabilities higher than 0.1:
#[1] pp                  magnitude.of.theta2 naive.SE.of.theta2 
#[4] rel.location       
#<0 rows> (or 0-length row.names)


## Adding stepping stone to estimate marginal likelihoods and compare bayes factors for model comparison
## Reference code from Uyeda et al
require(foreach)
require(doParallel)
registerDoParallel(cores=10)

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

currentTime <- Sys.time()
print(currentTime)

ss.SG.470s<-mcmc.SG.471$steppingstone(500000, SG.470s.comb, Bk = Bk, burnin=0.3, plot=FALSE)
ss.SG.470s$lnr
## lnr == 442.5791
pdf("SG.470s.comb_ss_plot.pdf")
plot(ss.SG.470s)
dev.off()
print("ss.SG.470s")
ss.SG.470s$lnr
currentTime <- Sys.time()
print(currentTime)

