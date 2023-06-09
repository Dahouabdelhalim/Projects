#############################################################################################
#	Brain Mass Analysis for Schumacher and Carlson 2022										#
#																							#
#                Combines new otophysan and osteoglossiform Sukhum et al 2018 data          #
#				 with published actinopterygian data from Tsuboi et al 2018					#
#				 and Tsuboi 2021															#
#																							#
#				 SS Model																	#
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


#################
#   SS models   #
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

prior.SS <- make.prior(tel_tree_rab, plot.prior = FALSE, 
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

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.471 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.471$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.471 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.471$model, prior=prior.primer.SS, startpar=model.primer.SS.471$startpar, file.dir=getwd(), outname="primerSS_r471", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.471$run(1000000)
primer.chain.SS <- primer.mcmc.SS.471$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.471$model)

## Building the rjMCMC models
set.seed(471)
model.SS.471 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.471 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.471$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r471", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(471)
mcmc.SS.471$run(2000000)		


#######
# 472 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.472 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.472$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.472 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.472$model, prior=prior.primer.SS, startpar=model.primer.SS.472$startpar, file.dir=getwd(), outname="primerSS_r472", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.472$run(1000000)
primer.chain.SS <- primer.mcmc.SS.472$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.472$model)


## Building the rjMCMC models
set.seed(472)
model.SS.472 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.472 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.472$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r472", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(472)
mcmc.SS.472$run(2000000)		


#######
# 473 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.473 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.473$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.473 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.473$model, prior=prior.primer.SS, startpar=model.primer.SS.473$startpar, file.dir=getwd(), outname="primerSS_r473", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.473$run(1000000)
primer.chain.SS <- primer.mcmc.SS.473$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.473$model)


## Building the rjMCMC models
set.seed(473)
model.SS.473 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.473 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.473$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r473", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(473)
mcmc.SS.473$run(2000000)		


#######
# 474 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.474 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.474$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.474 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.474$model, prior=prior.primer.SS, startpar=model.primer.SS.474$startpar, file.dir=getwd(), outname="primerSS_r474", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.474$run(1000000)
primer.chain.SS <- primer.mcmc.SS.474$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.474$model)


## Building the rjMCMC models
set.seed(474)
model.SS.474 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.474 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.474$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r474", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(474)
mcmc.SS.474$run(2000000)		


#######
# 475 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.475 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.475$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.475 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.475$model, prior=prior.primer.SS, startpar=model.primer.SS.475$startpar, file.dir=getwd(), outname="primerSS_r475", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.475$run(1000000)
primer.chain.SS <- primer.mcmc.SS.475$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.475$model)


## Building the rjMCMC models
set.seed(475)
model.SS.475 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.475 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.475$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r475", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(475)
mcmc.SS.475$run(2000000)		


#######
# 476 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.476 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.476$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.476 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.476$model, prior=prior.primer.SS, startpar=model.primer.SS.476$startpar, file.dir=getwd(), outname="primerSS_r476", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.476$run(1000000)
primer.chain.SS <- primer.mcmc.SS.476$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.476$model)


## Building the rjMCMC models
set.seed(476)
model.SS.476 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.476 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.476$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r476", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(476)
mcmc.SS.476$run(2000000)		


#######
# 477 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.477 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.477$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.477 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.477$model, prior=prior.primer.SS, startpar=model.primer.SS.477$startpar, file.dir=getwd(), outname="primerSS_r477", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.477$run(1000000)
primer.chain.SS <- primer.mcmc.SS.477$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.477$model)


## Building the rjMCMC models
set.seed(477)
model.SS.477 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.477 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.477$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r477", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(477)
mcmc.SS.477$run(2000000)		


#######
# 478 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.478 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.478$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.478 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.478$model, prior=prior.primer.SS, startpar=model.primer.SS.478$startpar, file.dir=getwd(), outname="primerSS_r478", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.478$run(1000000)
primer.chain.SS <- primer.mcmc.SS.478$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.478$model)


## Building the rjMCMC models
set.seed(478)
model.SS.478 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.478 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.478$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r478", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(478)
mcmc.SS.478$run(2000000)		


#######
# 479 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.479 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.479$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.479 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.479$model, prior=prior.primer.SS, startpar=model.primer.SS.479$startpar, file.dir=getwd(), outname="primerSS_r479", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.479$run(1000000)
primer.chain.SS <- primer.mcmc.SS.479$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.479$model)


## Building the rjMCMC models
set.seed(479)
model.SS.479 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.479 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.479$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r479", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(479)
mcmc.SS.479$run(2000000)		


#######
# 480 #
#######

## SS Model Primer
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
                 beta_BodyMass=rnorm(k+1, 0.599, 0.2), 
                 theta=rnorm(k+1, -1.653, .4), 
                 sb=sb, 
                 loc=rep(0, k), 
                 t2=(2:(k+1)))

bmax <- prob <- rep(1, nrow(cache$edge))
bmax[which(cache$edge.length == .Machine$double.eps)] <- 0

## Define primer priors (identical to priors below but with fixed shft locations randomly decided above)
prior.primer.SS <- make.prior(tree, plot.prior = FALSE, 
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
model.primer.SS.480 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.primer.SS, startpar=startpar, D=DSS.47)
model.primer.SS.480$model$moves$k<-".splitmergebd"
prior.primer.SS(startpar)

chainN <- "Primer"

primer.mcmc.SS.480 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.primer.SS.480$model, prior=prior.primer.SS, startpar=model.primer.SS.480$startpar, file.dir=getwd(), outname="primerSS_r480", plot.freq=NULL, ticker.freq=1000, samp = 100)

primer.mcmc.SS.480$run(1000000)
primer.chain.SS <- primer.mcmc.SS.480$load()

startpar.SS <- pull.pars(11, primer.chain.SS, model=model.primer.SS.480$model)


## Building the rjMCMC models
set.seed(480)
model.SS.480 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c("theta", "BodyMass"),  
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.SS, startpar=startpar.SS, D=DSS.47)

## Making the MCMC objects
mcmc.SS.480 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.SS.480$model, prior=prior.SS, startpar=startpar.SS, file.dir=getwd(), outname="modelSS_r480", plot.freq=NULL, ticker.freq=1000, samp = 100)

## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; 
set.seed(480)
mcmc.SS.480$run(2000000)		


#####################
#   Chain results   #
#####################

## Chain results
#chain.SS.471 <- set.burnin(mcmc.SS.471$load(), 0.3)
#chain.SS.472 <- set.burnin(mcmc.SS.472$load(), 0.3)
#chain.SS.473 <- set.burnin(mcmc.SS.473$load(), 0.3)
#chain.SS.474 <- set.burnin(mcmc.SS.474$load(), 0.3)
#chain.SS.475 <- set.burnin(mcmc.SS.475$load(), 0.3)
#chain.SS.476 <- set.burnin(mcmc.SS.476$load(), 0.3)
#chain.SS.477 <- set.burnin(mcmc.SS.477$load(), 0.3)
#chain.SS.478 <- set.burnin(mcmc.SS.478$load(), 0.3)
#chain.SS.479 <- set.burnin(mcmc.SS.479$load(), 0.3)
#chain.SS.480 <- set.burnin(mcmc.SS.480$load(), 0.3)

chain.SS.471 <- mcmc.SS.471$load()
chain.SS.472 <- mcmc.SS.472$load()
chain.SS.473 <- mcmc.SS.473$load()
chain.SS.474 <- mcmc.SS.474$load()
chain.SS.475 <- mcmc.SS.475$load()
chain.SS.476 <- mcmc.SS.476$load()
chain.SS.477 <- mcmc.SS.477$load()
chain.SS.478 <- mcmc.SS.478$load()
chain.SS.479 <- mcmc.SS.479$load()
chain.SS.480 <- mcmc.SS.480$load()


#out.SS.471 <- summary(chain.SS.471)
#out.SS.472 <- summary(chain.SS.472)
#out.SS.473 <- summary(chain.SS.473)
#out.SS.474 <- summary(chain.SS.474)
#out.SS.475 <- summary(chain.SS.475)
#out.SS.476 <- summary(chain.SS.476)
#out.SS.477 <- summary(chain.SS.477)
#out.SS.478 <- summary(chain.SS.478)
#out.SS.479 <- summary(chain.SS.479)
#out.SS.480 <- summary(chain.SS.480)

#pdf("SS.470s_trace_plots.pdf") 
# Creating a plot
#plot(chain.SS.471)
#plot(chain.SS.472)
#plot(chain.SS.473)
#plot(chain.SS.474)
#plot(chain.SS.475)
#plot(chain.SS.476)
#plot(chain.SS.477)
#plot(chain.SS.478)
#plot(chain.SS.479)
#plot(chain.SS.480)
# Closing the graphical device
#dev.off() 


## All 10 chains combined
SS.chains <- list(chain.SS.471,chain.SS.472,chain.SS.473,chain.SS.474,chain.SS.475,chain.SS.476,chain.SS.477,chain.SS.478,chain.SS.479,chain.SS.480)
SS.470s.comb <- combine.chains(SS.chains, burnin.prop = 0.3)

print("SS.470s.comb")
#summary(SS.470s.comb)
#pdf("SS.470s_comb_trace_plot.pdf")
#plot(SS.470s.comb)
#dev.off()

## summary output of comb chain:
#bayou MCMC chain: 2000001 generations
#140047 samples, first 0 samples discarded as burnin
#
#
#Summary statistics for parameters:
#                            Mean           SD     Naive SE Time-series SE
#lnL                 4.980759e+02 7.571068e+00 2.023114e-02   2.620453e-01
#prior              -4.234599e+01 1.195275e+01 3.193971e-02   4.393920e-01
#alpha               5.468204e-04 4.917752e-04 1.314104e-06   3.917682e-06
#sig2                4.016560e-04 2.823501e-05 7.544857e-08   2.771445e-07
#k                   6.805330e+00 1.916978e+00 5.122479e-03   7.214336e-02
#ntheta              7.805330e+00 1.916978e+00 5.122479e-03   7.214336e-02
#root.theta         -1.879346e+00 1.605488e-01 4.290128e-04   1.192547e-03
#root.beta_BodyMass  5.268315e-01 1.402416e-02 3.747486e-05   5.813447e-04
#all theta          -1.660805e+00 3.919596e-01           NA             NA
#all beta_BodyMass   6.046931e-01 1.620217e-01           NA             NA
#                   Effective Size    HPD95Lower    HPD95Upper
#lnL                      834.7598  4.835484e+02  5.129733e+02
#prior                    739.9996 -6.819868e+01 -2.248792e+01
#alpha                  15757.0482  7.606009e-09  1.536815e-03
#sig2                   10379.1848  3.494511e-04  4.595564e-04
#k                        706.0598  3.000000e+00  1.000000e+01
#ntheta                   706.0598  4.000000e+00  1.100000e+01
#root.theta             18124.3925 -2.195011e+00 -1.567698e+00
#root.beta_BodyMass       581.9510  5.018367e-01  5.510874e-01
#all theta                      NA            NA            NA
#all beta_BodyMass              NA            NA            NA
#
#
#Branches with posterior probabilities higher than 0.1:
#            pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
#14   0.9779503           -1.648736        0.001084798    0.4989339
#1730 0.4784822           -1.544820        0.001534337    0.4897422
#144  0.4359251           -1.652659        0.001615017    0.5067931
#1059 0.4299200           -1.457929        0.001644126    0.5000082
#1035 0.2605411           -1.680531        0.002091286    0.5002670
#86   0.2292802           -1.656584        0.002228004    0.4935275
#537  0.2246674           -1.652785        0.002266366    0.5029421
#145  0.2223039           -1.654455        0.002248587    0.4930654
#1061 0.1759267           -1.426565        0.002603732    0.5034896

## Adding stepping stone to estimate marginal likelihoods and compare bayes factors for model comparison
## Reference code from Uyeda et al
require(foreach)
require(doParallel)
registerDoParallel(cores=10)

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

currentTime <- Sys.time()
print(currentTime)

ss.SS.470s<-mcmc.SS.471$steppingstone(500000, SS.470s.comb, Bk = Bk, burnin=0.3, plot=FALSE)
ss.SS.470s$lnr
## lnr == 463.6921
pdf("SS.470s.comb_ss_plot.pdf")
plot(ss.SS.470s)
dev.off()
print("ss.SS.470s")
ss.SS.470s$lnr
currentTime <- Sys.time()
print(currentTime)


#########################
#   identified shifts   #
#########################

## using pp > .2 following Uyeda et al
shiftsum.SS.470s <- shiftSummaries(SS.470s.comb, mcmc.SS.471, pp.cutoff=0.2)

## identified shifts
## 14 == "Gymnothorax_meleagris" [turkey moray]
## 86 == "Isichthys_henryi"         "Brienomyrus_brachyistius"
## 144 == "Synodontis_multipunctatus"
## 145 == "Synodontis_petricola"
## 537 == "Arothron_nigropunctatus" [black spotted puffer]
## 1035 == Gasterosteiformes [sticklebacks]; Perciforms: Stichaeidae [pricklebacks], Pholidae [gunnels]; Scorpaeniformes: Anarhichadidae [wolffish], Agonidae [poachers], Cottidae [cottids], Cyclopteridae [lumpsuckers], Liparidae [snailfish]
## 1059 == Lophiiformes [anglerfish]; Tetraodontiformes [puffers, boxfish, filefish, triggerfish etc]; Acanthuriformes [surgeonfish, moorish angels]; Perciformes: Cepolidae [bandfishes], Ephippidae [spadefishes], Siganidae [rabbitfishes], Caproidae [boarfishes], Scatophagidae [scats], Priacanthidae [bigeyes], Ammodytidae [sandlance], Caesionidae [fusiliers], Centrarchidae [sunfish], Chaetodontidae [butterflyfish], Channichthyidae [crocodile icefish], Cirrhitidae [hawkfish], Epigonidae [deepwater cardinal fish], Gerreidae [mojarras], Haemulidae [grunts], Kuhliidae [flagtails], Nemipteridae [threadfin breams], Nototheniidae [cod icefish], Percidae [perch + relatives], Pinguipedidae [sandperch], Sparidae [porgies], Terapontidae [grunters], Trachinidae [weevers]; Scorpeaniformes: Congiopodidae [pigfish], Hoplichthyidae [ghostflatheads], Platycephalidae [flatheads], Scorpaenidae [scorpionfish], Sebastidae [rockfishes], Synanceiidae [stonefish], Triglidae [sea robins], 
## 1730 == osteoglossiformes

#bayou fitted regressions
#         theta beta_BodyMass
#root -1.880559     0.5278382
#14   -1.648572     0.8131417
#86   -1.653568     0.4611389
#144  -1.650950     0.7365777
#145  -1.656488     0.7847541
#537  -1.653491     0.5974406
#1035 -1.681612     0.3852765
#1059 -1.450751     0.4657941
#1730 -1.540541     0.7215232

subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

subroot<-rbind(subroot,sub145)

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub1035,sub1059,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))


## Pagel pgls [fitted lambda = 0.9859765], grade N > 2
sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species))
lam<-as.numeric(sub_pgls$modelStruct)

pgls0<-gls(BrainMass~BodyMass, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species))
pgls1035<-gls(BrainMass~BodyMass, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species))
pgls1059<-gls(BrainMass~BodyMass, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species))
pgls1730<-gls(BrainMass~BodyMass, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species))

anova(sub_pgls)

## posthoc intercept
pairs(emmeans(sub_pgls, "code")) ## does a tukey pairwise comparison!

## posthoc slope 
pairs(emtrends(sub_pgls,"code",var="BodyMass",mode = "df.error"))

## output of this:
#> anova(sub_pgls)
#Denom. DF: 856 
#              numDF  F-value p-value
#(Intercept)       1   15.814  0.0001
#BodyMass          1 6079.447  <.0001
#code              3    4.581  0.0034
#BodyMass:code     3   18.496  <.0001
#> pairs(emmeans(sub_pgls, "code"))
#NOTE: Results may be misleading due to involvement in interactions
# contrast    estimate     SE  df t.ratio p.value
# 0 - 1035       0.216 0.1451 410   1.491  0.4442
# 0 - 1059      -0.084 0.0888 402  -0.946  0.7802
# 0 - 1730      -0.687 0.2091 400  -3.287  0.0060
# 1035 - 1059   -0.300 0.1165 419  -2.577  0.0503
# 1035 - 1730   -0.903 0.2540 403  -3.556  0.0024
# 1059 - 1730   -0.603 0.2267 400  -2.661  0.0403
#
#Degrees-of-freedom method: satterthwaite 
#P value adjustment: tukey method for comparing a family of 4 estimates 
#> pairs(emtrends(sub_pgls,"code",var="BodyMass",mode = "df.error"))
# contrast    estimate     SE  df t.ratio p.value
# 0 - 1035      0.1500 0.0316 854   4.752  <.0001
# 0 - 1059      0.0649 0.0145 854   4.483  <.0001
# 0 - 1730     -0.1831 0.0494 854  -3.710  0.0013
# 1035 - 1059  -0.0851 0.0327 854  -2.604  0.0462
# 1035 - 1730  -0.3331 0.0574 854  -5.802  <.0001
# 1059 - 1730  -0.2480 0.0501 854  -4.951  <.0001
#
#Degrees-of-freedom method: df.error 
#P value adjustment: tukey method for comparing a family of 4 estimates 


####################
#     figure 3     #
####################

## 14 == "Gymnothorax_meleagris" [turkey moray]
## 86 == "Isichthys_henryi"         "Brienomyrus_brachyistius"
## 144 == "Synodontis_multipunctatus"
## 145 == "Synodontis_petricola"
## 537 == "Arothron_nigropunctatus" [black spotted puffer]
## 1035 == Gasterosteiformes [sticklebacks]; Perciforms: Stichaeidae [pricklebacks], Pholidae [gunnels]; Scorpaeniformes: Anarhichadidae [wolffish], Agonidae [poachers], Cottidae [cottids], Cyclopteridae [lumpsuckers], Liparidae [snailfish]
## 1059 == Lophiiformes [anglerfish]; Tetraodontiformes [puffers, boxfish, filefish, triggerfish etc]; Acanthuriformes [surgeonfish, moorish angels]; Perciformes: Cepolidae [bandfishes], Ephippidae [spadefishes], Siganidae [rabbitfishes], Caproidae [boarfishes], Scatophagidae [scats], Priacanthidae [bigeyes], Ammodytidae [sandlance], Caesionidae [fusiliers], Centrarchidae [sunfish], Chaetodontidae [butterflyfish], Channichthyidae [crocodile icefish], Cirrhitidae [hawkfish], Epigonidae [deepwater cardinal fish], Gerreidae [mojarras], Haemulidae [grunts], Kuhliidae [flagtails], Nemipteridae [threadfin breams], Nototheniidae [cod icefish], Percidae [perch + relatives], Pinguipedidae [sandperch], Sparidae [porgies], Terapontidae [grunters], Trachinidae [weevers]; Scorpeaniformes: Congiopodidae [pigfish], Hoplichthyidae [ghostflatheads], Platycephalidae [flatheads], Scorpaenidae [scorpionfish], Sebastidae [rockfishes], Synanceiidae [stonefish], Triglidae [sea robins], 
## 1730 == osteoglossiformes

subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

subroot<-rbind(subroot,sub145)

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub1059,sub1730,sub1035,sub14,sub86,sub144,sub537)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

Species<-rownames(subroot)
subroot<-cbind(subroot,Species)
Species<-rownames(sub1035)
sub1035<-cbind(sub1035,Species)
Species<-rownames(sub1059)
sub1059<-cbind(sub1059,Species)
Species<-rownames(sub1730)
sub1730<-cbind(sub1730,Species)


## Pagel pgls [fitted lambda = 0.9859765]
sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corPagel(0.9859765,phy=subtree,form=~Species,fixed=TRUE), control = list(singular.ok = TRUE))
lam<-as.numeric(sub_pgls$modelStruct)

pgls0<-gls(BrainMass~BodyMass, data=subroot, correlation=corPagel(lam,phy=drop.tip(subtree, setdiff(subtree$tip.label, rownames(subroot))),form=~Species,fixed=TRUE))
pgls1035<-gls(BrainMass~BodyMass, data=sub1035, correlation=corPagel(lam,phy=drop.tip(subtree, setdiff(subtree$tip.label, rownames(sub1035))),form=~Species,fixed=TRUE))
pgls1059<-gls(BrainMass~BodyMass, data=sub1059, correlation=corPagel(lam,phy=drop.tip(subtree, setdiff(subtree$tip.label, rownames(sub1059))),form=~Species,fixed=TRUE))
pgls1730<-gls(BrainMass~BodyMass, data=sub1730, correlation=corPagel(lam,phy=drop.tip(subtree, setdiff(subtree$tip.label, rownames(sub1730))),form=~Species,fixed=TRUE))


code_colors<-c("0" = "grey", "14" = 'lightgreen',
			   "86" = '#CCBB44',
			   "144" = '#FFAABB',
			   "537" = "#D55E00",
		   	   "1035" = '#332288',
			   "1059" = '#88CCEE',
			   "1730" = '#AA4499')
code_shapes<-c("0" = 21, "14" = 22,
			   "86" = 22,
			   "144" = 22,
			   "537" = 22,
		   	   "1035" = 21,
			   "1059" = 21,
			   "1730" = 21)


actino_brain_size_plot<-ggplot(subpgls, aes(x=BodyMass, y=BrainMass, fill=code, colour=code, shape=code),size=.1) + geom_point() + scale_shape_manual(values=code_shapes) + theme_classic() + theme(legend.position='None') + labs(x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_fill_manual(values = code_colors) + scale_colour_manual(values = c("black","black","black","black","black","black","black","black")) +
	geom_abline(slope = pgls0$coefficients[2], intercept = pgls0$coefficients[1],color="grey") + 
	geom_abline(slope = pgls1035$coefficients[2], intercept = pgls1035$coefficients[1],color='#332288') +
	geom_abline(slope = pgls1059$coefficients[2], intercept = pgls1059$coefficients[1],color='#88CCEE') +
	geom_abline(slope = pgls1730$coefficients[2], intercept = pgls1730$coefficients[1],color='#AA4499')

pdf("bayou_grades_pgls_brainsize_plot_updatedApr25.pdf", width = 8, height = 5)
actino_brain_size_plot
dev.off()

#pdf("bayou_grades_pgls_brainsize_plot_updated.pdf", width = 5, height = 5)
#actino_brain_size_plot
#dev.off()

#actino_brain_size_plot2<-ggplot(subpgls, aes(x=BodyMass, y=BrainMass, colour=code),size=.1) + geom_point() + theme_classic() + theme(legend.position='None') + labs(x=bquote('Log Body Mass' ~(g)), y = bquote('Log Brain Mass' ~(g))) + scale_colour_manual(values = code_colors) +
#	geom_abline(slope = pgls0$coefficients[2], intercept = pgls0$coefficients[1],color="grey") + 
#	geom_abline(slope = pgls1035$coefficients[2], intercept = pgls1035$coefficients[1],color='#332288') +
#	geom_abline(slope = pgls1059$coefficients[2], intercept = pgls1059$coefficients[1],color='#88CCEE') +
#	geom_abline(slope = pgls1730$coefficients[2], intercept = pgls1730$coefficients[1],color='#AA4499')

#pdf("bayou_grades_pgls_brainsize_plot_noout.pdf", width = 4, height = 4)
#actino_brain_size_plot2
#dev.off()


## cladogram of ensen phenos
#clado_type<-c("tubegen","tubegen","tubegen","tubegen","tubegen","out","oamp","out","out","tubegen","tubegen","tubegen","tubegen","tubegen","out","tubegen","tubegen","oamp","oamp","tubegen","out","tubegen","out","tubegen","tubegen","tubegen","tubegen","ampegen","ampegen","ampegen","ampegen","oamp","out","out","out","out")
#names(clado_type)<-full_clado_names

cls <- list(c1=rownames(subroot),
            c2=rownames(sub14),
            c3=rownames(sub86),
            c4=rownames(sub144),
            c5=rownames(sub537),
            c6=rownames(sub1035),
            c7=rownames(sub1059),
            c8=rownames(sub1730))
subtree <- groupOTU(subtree, cls)
clado_group_vals<-c(c1="grey", c2='lightgreen', c3='#CCBB44', c4='#FFAABB', c5="#D55E00", c6='#332288', c7='#88CCEE', c8='#AA4499')
#clado_group_outs<-c(c1="magenta", c2="magenta", c3="black", c4="darkolivegreen3", c5="darkolivegreen3", c6="black", c7="black", c8="black")

#library(ggnewscale)

## labels internal nodes
#paper_clado<-ggtree(subtree, size=.45, layout="circular") +
#	geom_tree(layout="circular", aes(color=group), size=.35) + 
#	geom_tiplab(aes(angle=angle), size=.4,offset=.2) + scale_color_manual(values=c(clado_group_vals)) +
#	geom_text(aes(label=node),size=.5)

paper_clado<-ggtree(subtree, size=.45, layout="circular") +
	geom_tree(layout="circular", aes(color=group), size=.35) + 
	geom_tiplab(aes(angle=angle), size=.4,offset=.2) + scale_color_manual(values=c(clado_group_vals)) +
	geom_text(aes(label=node),size=.5) +
	geom_cladelabel(node=1690, label="Osteoglossiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=870, label="Polypteriformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1739, label="Acipenseriformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=866, label="Lepisosteiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=867, label="Amiiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=865, label="Notacanthiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1718, label="Anguilliformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1689, label="Alepocephaliformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1683, label="Clupeiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1654, label="Cypriniformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1674, label="Gymnotiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1671, label="Characiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1660, label="Siluriformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1650, label="Argentiniformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1649, label="Esociformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1644, label="Salmoniformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=765, label="Osmeriformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1633, label="Stomiiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1628, label="Aulopiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1623, label="Myctophiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=748, label="Zeiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1590, label="Gadiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=690, label="Stephanoberyciformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1573, label="Beryciformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1575, label="Holocentriformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1570, label="Ophidiiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1207, label="Gempylidae+Centrolophidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1194, label="Scombridae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1189, label="Fistulariidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=304, label="Pegasidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1190, label="Dactylopteridae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1185, label="Mullidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1188, label="Callionymidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=295, label="Aulostomidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=294, label="Centriscidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1165, label="Syngnathidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1154, label="Kurtiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1074, label="Gobiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=187, label="Anabantiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1041, label="Sphyraenidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=186, label="Polynemidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=185, label="Menidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1071, label="Istiophoridae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=184, label="Xiphiidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1047, label="Carangidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=180, label="Coryphaenidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=181, label="Echeneidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1038, label="Mugiliformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=101, label="Grammatidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=100, label="Opistognathidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=99, label="Gobiesociformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=984, label="Blenniiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=989, label="Pomacentriidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=902, label="Atheriniformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=892, label="Cyprinodontiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=893, label="Beloniformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=905, label="Cichliformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1498, label="Labriformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=613, label="Pempheridae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=614, label="Glaucosomatidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1495, label="Lophiiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=609, label="Cepolidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1462, label="Tetraodontiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=560, label="Ephippidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1449, label="Siganidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=570, label="Caproidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1459, label="Scatophagidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1458, label="Priacanthidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1418, label="Acanthuriformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1238, label="Triglidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=364, label="Platycephalidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=365, label="Hoplichthyidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=362, label="Synanceiidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=368, label="Congiopodidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1246, label="Sebastidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1251, label="Scorpaenidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1236, label="Gasterosteiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=337, label="Stichaeidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=338, label="Pholidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1234, label="Anarhichidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1230, label="Agonidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1227, label="Cottidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=333, label="Cyclopteridae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=332, label="Liparidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1295, label="Channichthyidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1293, label="Nototheniidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1260, label="Serranidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=420, label="Ammodytidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1338, label="Caesionidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1304, label="Centrarchidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1356, label="Chaetodontidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1305, label="Cirrhitidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=418, label="Epigonidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1312, label="Gerreidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1316, label="Haemulidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=416, label="Kuhliidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=417, label="Kyphosidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1408, label="Lethrinidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1323, label="Lutjanidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1340, label="Lutjanidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=424, label="Moronidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1405, label="Nemipteridae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1299, label="Percidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=419, label="Pinguipedidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1343, label="Pomacanthidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1394, label="Sparidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=415, label="Terapontidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2) +
	geom_cladelabel(node=1300, label="Trachinidae", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5, size=2)

#rotate(1067)
#rotate(1245)
#rotate(1294)
#rotate(1322)
#	+
#	new_scale_color() +
#	xlim(0, 30) + 
#	geom_cladelabel(node=63, label="Osteoglossiformes", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5) +
#	geom_cladelabel(node=66, label="Mormyroids", angle=270,  offset = 12, offset.text=.7, barsize=1.5, hjust=.5) +
#	geom_cladelabel(node=42, label="Otophysans", angle=270,  offset = 14, offset.text=.7, barsize=1.5, hjust=.5) +
#	geom_cladelabel(node=51, label="Gymnotiformes", angle=270,  offset = 12, offset.text=.7, barsize=1.5, hjust=.5) +
#	geom_cladelabel(node=45, label="Siluriformes", angle=270,  offset = 12, offset.text=.7, barsize=1.5, hjust=.5)
#paper_clado2<-paper_clado %>% rotate(37) %>% rotate(38) %>% rotate(61) %>% rotate(50) %>% rotate(55) 


pdf("bayou_grades_phylo_test.pdf", width = 8.5, height = 11)
paper_clado
dev.off()

paper_clado2 <- ggtree(subtree, size=.45, aes(color=group), layout="circular") +
#	geom_tiplab(aes(angle=angle), size=.4,offset=.2) + 
	scale_color_manual(values=c(clado_group_vals)) +
	geom_cladelabel(node=1690, label="Osteoglossiformes", angle=287,  offset = .1, offset.text=1.8, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=870, label="Polypteriformes", angle=0,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1739, label="Acipenseriformes", angle=0,  offset = .1, offset.text=25, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=866, label="Lepisosteiformes", angle=0,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=867, label="Amiiformes", angle=0,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=865, label="Notacanthiformes", angle=1,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1718, label="Anguilliformes", angle=283,  offset = .1, offset.text=1.8, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1689, label="Alepocephaliformes", angle=20,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1683, label="Clupeiformes", angle=23,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1654, label="Cypriniformes", angle=25,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1674, label="Gymnotiformes", angle=30,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1671, label="Characiformes", angle=32,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1660, label="Siluriformes", angle=36,  offset = .1, offset.text=16, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1650, label="Argentiniformes", angle=40,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1649, label="Esociformes", angle=40,  offset = .1, offset.text=16, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1644, label="Salmoniformes", angle=42,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=765, label="Osmeriformes", angle=45,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1633, label="Stomiiformes", angle=50,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1628, label="Aulopiformes", angle=52,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1623, label="Myctophiformes", angle=55,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=743, label="Zeiformes", angle=58,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1590, label="Gadiformes", angle=330,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=690, label="Stephanoberyciformes", angle=69,  offset = .1, offset.text=25, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1573, label="Beryciformes", angle=72,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1575, label="Holocentriformes", angle=76,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1570, label="Ophidiiformes", angle=78,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1207, label="Gempylidae+Centrolophidae", angle=80,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1194, label="Scombridae", angle=83,  offset = .1, offset.text=16, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1189, label="Fistulariidae", angle=89,  offset = .1, offset.text=16, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=304, label="Pegasidae", angle=88,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1190, label="Dactylopteridae", angle=86,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1185, label="Mullidae", angle=90,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1188, label="Callionymidae", angle=89,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=295, label="Aulostomidae", angle=90,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=294, label="Centriscidae", angle=90,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1165, label="Syngnathidae", angle=92,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1154, label="Kurtiformes", angle=100,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1074, label="Gobiformes", angle=31,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=187, label="Anabantiformes", angle=-45,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1041, label="Sphyraenidae", angle=-45,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=186, label="Polynemidae", angle=-42,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=185, label="Menidae", angle=-42,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1071, label="Istiophoridae", angle=-40,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=184, label="Xiphiidae", angle=-42,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1047, label="Carangidae", angle=54,  offset = .1, offset.text=1.8, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=180, label="Coryphaenidae", angle=-40,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=181, label="Echeneidae", angle=-40,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1038, label="Mugiliformes", angle=-32,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=101, label="Grammatidae", angle=-32,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=100, label="Opistognathidae", angle=-32,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=99, label="Gobiesociformes", angle=-32,  offset = .1, offset.text=23, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=984, label="Blenniiformes", angle=-30,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=989, label="Pomacentriidae", angle=72,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=902, label="Atheriniformes", angle=-8,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=892, label="Cyprinodontiformes", angle=-4,  offset = .1, offset.text=23, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=893, label="Beloniformes", angle=0,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=905, label="Cichliformes", angle=108,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1498, label="Labriformes", angle=-45,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=613, label="Pempheridae", angle=57,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=614, label="Glaucosomatidae", angle=59,  offset = .1, offset.text=22, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1495, label="Lophiiformes", angle=61,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=609, label="Cepolidae", angle=62,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1462, label="Tetraodontiformes", angle=-20,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=560, label="Ephippidae", angle=86,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1449, label="Siganidae", angle=84,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=570, label="Caproidae", angle=82,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1459, label="Scatophagidae", angle=80,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1458, label="Priacanthidae", angle=80,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1418, label="Acanthuriformes", angle=0,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1238, label="Triglidae", angle=-26,  offset = .1, offset.text=13, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=364, label="Platycephalidae", angle=-30,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=365, label="Hoplichthyidae", angle=-30,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=362, label="Synanceiidae", angle=-30,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=363, label="Congiopodidae", angle=-30,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1246, label="Sebastidae", angle=-28,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1251, label="Scorpaenidae", angle=-28,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1249, label="Scorpaenidae", angle=-27,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1236, label="Gasterosteiformes", angle=-23,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=337, label="Stichaeidae", angle=-23,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=338, label="Pholidae", angle=-23,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1234, label="Anarhichidae", angle=-21,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1230, label="Agonidae", angle=-18,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1227, label="Cottidae", angle=-18,  offset = .1, offset.text=13, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=333, label="Cyclopteridae", angle=-20,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=332, label="Liparidae", angle=-20,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1295, label="Channichthyidae", angle=-15,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1296, label="Nototheniidae", angle=-15,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1260, label="Serranidae", angle=263,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=420, label="Ammodytidae", angle=-38,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1338, label="Caesionidae", angle=-40,  offset = .1, offset.text=18, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1304, label="Centrarchidae", angle=-34,  offset = .1, offset.text=19, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1356, label="Chaetodontidae", angle=24,  offset = .1, offset.text=1.5, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1305, label="Cirrhitidae", angle=-33,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=418, label="Epigonidae", angle=-36,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1312, label="Gerreidae", angle=-36,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1316, label="Haemulidae", angle=-49,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=416, label="Kuhliidae", angle=-36,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=417, label="Kyphosidae", angle=-36,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1408, label="Lethrinidae", angle=-80,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1323, label="Lutjanidae", angle=-45,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1340, label="Lutjanidae", angle=-45,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=424, label="Moronidae", angle=-38,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1405, label="Nemipteridae", angle=-78,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1299, label="Percidae", angle=-31,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=419, label="Pinguipedidae", angle=-38,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1343, label="Pomacanthidae", angle=-52,  offset = .1, offset.text=19, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1394, label="Sparidae", angle=-72,  offset = .1, offset.text=13, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=415, label="Terapontidae", angle=-36,  offset = .1, offset.text=20, barsize=1, hjust=.5, fontsize=1) +
	geom_cladelabel(node=1300, label="Trachinidae", angle=-30,  offset = .1, offset.text=15, barsize=1, hjust=.5, fontsize=1)

paper_clado3 <- paper_clado2 %>% rotate(1067) %>% rotate(1245) %>% rotate(1294) %>% rotate(1322) %>% rotate(1337)


#pdf("bayou_grades_phylo_noout.pdf", width = 8.5, height = 11)
#paper_clado2
#dev.off()

pdf("bayou_grades_phylo_fig3a.pdf", width = 8, height = 8)
paper_clado3
dev.off()


########################################################################################
# confirming with pgls by systematically excluding grades and calc model fit using AIC #
########################################################################################

## just with grades N >= 3, pagel's lambda
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub1035,sub1059,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species), na.action = na.omit, method = "ML")
AIC(sub_pgls)
## lambda = 0.9856621
## AIC = -1057.634
## loglik = 538.8168


## collapsing gr 1035 to ancestral 1059
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(1059,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub1035,sub1059,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species), na.action = na.omit, method = "ML")
AIC(sub_pgls)
## lambda = 0.986091
## AIC = -1048.929
## loglik = 532.4643


## collapsing gr 1059 to ancestral root
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(0,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub1035,sub1059,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species), na.action = na.omit, method = "ML")
AIC(sub_pgls)
## lambda = 0.9847622
## AIC = -1041.134
## loglik = 528.5672


## collapsing gr 1730 to ancestral root
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(0,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub1035,sub1059,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species), na.action = na.omit, method = "ML")
AIC(sub_pgls)
## lambda = 0.9849083
## AIC = -1039.097
## loglik = 527.5486


## collapsing all grades to ancestral root, i.e. no shifts
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(0,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(0,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(0,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub1035,sub1059,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass, data=subpgls, correlation=corPagel(.8,phy=subtree,form=~Species), na.action = na.omit, method = "ML")
AIC(sub_pgls)
## lambda = 0.9840213 
## AIC = -1002.688
## loglik = 505.3438



## All Grades

## fit incl all grades
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -1019.356
## loglik = 528.6781

#X <- model.matrix(BrainMass~BodyMass+code+BodyMass:code, data=subpgls) 
#summary(X)


## collapse gr 14 to ancestral root
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(0,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -995.9932
## loglik = 514.9966


## collapsing gr 86 to ancestral 1730
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(1730,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -1010.933
## loglik = 522.4665


## collapsing gr 144 to ancestral root
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(0,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -1012.728
## loglik = 523.364


## collapsing gr 145 to ancestral root
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(0,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -1023.133
## loglik = 528.5667


## collapsing gr 537 to ancestral 1059
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(1059,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -985.5852
## loglik = 509.7926


## collpasing gr 1035 to ancestral 1059
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1059,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -1011.758
## loglik = 522.8789


## collapsing gr 1059 to ancestral root
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(0,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -999.7714
## loglik = 516.8857


## collapsing gr 1730 to ancestral root
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(0,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass+code+BodyMass:code, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML", control = list(singular.ok = TRUE))
AIC(sub_pgls)
## AIC = -996.3034
## loglik = 515.1517


## collapsing all grades to the ancestral root, i.e. no shifts
subroot<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[1]])
sub14<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[2]])
sub86<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[3]])
sub144<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[4]])
sub145<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[5]])
sub537<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[6]])
sub1035<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[7]])
sub1059<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[8]])
sub1730<-subset(tel_avg_data, rownames(tel_avg_data) %in% shiftsum.SS.470s$descendents[[9]])

code<-as.factor(rep(0,length(subroot$BodyMass)))
subroot<-cbind(subroot,code)
code<-as.factor(rep(14,length(sub14$BodyMass)))
sub14<-cbind(sub14,code)
code<-as.factor(rep(86,length(sub86$BodyMass)))
sub86<-cbind(sub86,code)
code<-as.factor(rep(144,length(sub144$BodyMass)))
sub144<-cbind(sub144,code)
code<-as.factor(rep(145,length(sub145$BodyMass)))
sub145<-cbind(sub145,code)
code<-as.factor(rep(537,length(sub537$BodyMass)))
sub537<-cbind(sub537,code)
code<-as.factor(rep(1035,length(sub1035$BodyMass)))
sub1035<-cbind(sub1035,code)
code<-as.factor(rep(1059,length(sub1059$BodyMass)))
sub1059<-cbind(sub1059,code)
code<-as.factor(rep(1730,length(sub1730$BodyMass)))
sub1730<-cbind(sub1730,code)
subpgls<-rbind(subroot,sub14,sub86,sub144,sub145,sub537,sub1035,sub1059,sub1730)
#subpgls<-rbind(subroot,sub86,sub1035,sub1059,sub1061,sub1730)
Species<-rownames(subpgls)
subpgls<-cbind(subpgls,Species)
subtree<-drop.tip(tel_tree_rab, setdiff(tel_tree_rab$tip.label, rownames(subpgls)))

sub_pgls<-gls(BrainMass~BodyMass, data=subpgls, correlation=corBrownian(phy=subtree,form=~Species), na.action = na.omit, method = "ML")
AIC(sub_pgls)
## AIC = -886.7848
## loglik = 446.3924
