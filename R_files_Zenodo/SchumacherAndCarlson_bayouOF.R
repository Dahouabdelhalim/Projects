#############################################################################################
#	Brain Mass Analysis for Schumacher and Carlson 2022										#
#																							#
#                Combines new otophysan and osteoglossiform Sukhum et al 2018 data          #
#				 with published actinopterygian data from Tsuboi et al 2018					#
#				 and Tsuboi 2021															#
#																							#
#				 OF Model																	#
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

DSS.48 = list(alpha=7, sig2=.7, beta_BodyMass=0.4, k=1, theta=3, slide=1)


##########################
#   fixed shift chains   #
##########################

## First, make cache object to pull info about brnch nums for bayou sampling function
tree <- tel_tree_rab
cache <- bayou:::.prepare.ou.univariate(tree, setNames(tel_avg_data$BrainMass, tree$tip.label), SE=BrainSE)
pred <- as.matrix(BodyMass)
colnames(pred) <- "BodyMass"

## SS Model, Hypothesis Osteoglossiform
## k = 1 shift at osteoglossiforms
## osteoglossiform sb = 1730

k <- 1
sb <- 1730

## draw primer chain starting params from prior dists using randmly gened shft num and locs from above
## Updated alpha and sig2 priors to reflect those in Ksepka et al 20 and Smaers et al 21, which modeled brain allometries specifically (scale = 0.1 from 1 used in Uyeda et al)
set.seed(161)
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
set.seed(161)
model.fixed.SS.161 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.161$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.161 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.161$model, prior=prior.fixed.SS, startpar=model.fixed.SS.161$startpar, file.dir=getwd(), outname="OSTfixedSS_r161", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.161$run(2000000)
chain.SS.161 <- fixed.mcmc.SS.161$load()


set.seed(162)
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
set.seed(162)
model.fixed.SS.162 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.162$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.162 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.162$model, prior=prior.fixed.SS, startpar=model.fixed.SS.162$startpar, file.dir=getwd(), outname="OSTfixedSS_r162", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.162$run(2000000)
chain.SS.162 <- fixed.mcmc.SS.162$load()


set.seed(163)
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
set.seed(163)
model.fixed.SS.163 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.163$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.163 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.163$model, prior=prior.fixed.SS, startpar=model.fixed.SS.163$startpar, file.dir=getwd(), outname="OSTfixedSS_r163", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.163$run(2000000)
chain.SS.163 <- fixed.mcmc.SS.163$load()


set.seed(164)
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
set.seed(164)
model.fixed.SS.164 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.164$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.164 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.164$model, prior=prior.fixed.SS, startpar=model.fixed.SS.164$startpar, file.dir=getwd(), outname="OSTfixedSS_r164", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.164$run(2000000)
chain.SS.164 <- fixed.mcmc.SS.164$load()


set.seed(165)
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
set.seed(165)
model.fixed.SS.165 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.165$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.165 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.165$model, prior=prior.fixed.SS, startpar=model.fixed.SS.165$startpar, file.dir=getwd(), outname="OSTfixedSS_r165", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.165$run(2000000)
chain.SS.165 <- fixed.mcmc.SS.165$load()


set.seed(166)
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
set.seed(166)
model.fixed.SS.166 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.166$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.166 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.166$model, prior=prior.fixed.SS, startpar=model.fixed.SS.166$startpar, file.dir=getwd(), outname="OSTfixedSS_r166", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.166$run(2000000)
chain.SS.166 <- fixed.mcmc.SS.166$load()


set.seed(167)
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
set.seed(167)
model.fixed.SS.167 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.167$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.167 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.167$model, prior=prior.fixed.SS, startpar=model.fixed.SS.167$startpar, file.dir=getwd(), outname="OSTfixedSS_r167", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.167$run(2000000)
chain.SS.167 <- fixed.mcmc.SS.167$load()


set.seed(168)
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
set.seed(168)
model.fixed.SS.168 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.168$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.168 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.168$model, prior=prior.fixed.SS, startpar=model.fixed.SS.168$startpar, file.dir=getwd(), outname="OSTfixedSS_r168", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.168$run(2000000)
chain.SS.168 <- fixed.mcmc.SS.168$load()


set.seed(169)
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
set.seed(169)
model.fixed.SS.169 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.169$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.169 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.169$model, prior=prior.fixed.SS, startpar=model.fixed.SS.169$startpar, file.dir=getwd(), outname="OSTfixedSS_r169", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.169$run(2000000)
chain.SS.169 <- fixed.mcmc.SS.169$load()


set.seed(170)
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
set.seed(170)
model.fixed.SS.170 <- makeBayouModel(BrainMass ~ BodyMass, rjpars=c("theta","BodyMass"), tree=tree, dat=BrainMass, pred=pred, SE=BrainSE, prior=prior.fixed.SS, startpar=startpar, D=DSS.48)
model.fixed.SS.170$model$moves$k<-".splitmergebd"
prior.fixed.SS(startpar)

fixed.mcmc.SS.170 <- bayou.makeMCMC(tree, dat=BrainMass, pred=as.matrix(BodyMass), SE=BrainSE, model=model.fixed.SS.170$model, prior=prior.fixed.SS, startpar=model.fixed.SS.170$startpar, file.dir=getwd(), outname="OSTfixedSS_r170", plot.freq=NULL, ticker.freq=1000, samp = 100)

fixed.mcmc.SS.170$run(2000000)
chain.SS.170 <- fixed.mcmc.SS.170$load()



#####################
#   chain results   #
#####################

#chain.SS.161 <- set.burnin(fixed.mcmc.SS.161$load(), 0.3)
#out.SS.161 <- summary(chain.SS.161)
#chain.SS.162 <- set.burnin(fixed.mcmc.SS.162$load(), 0.3)
#out.SS.162 <- summary(chain.SS.162)
#chain.SS.163 <- set.burnin(fixed.mcmc.SS.163$load(), 0.3)
#out.SS.163 <- summary(chain.SS.163)
#chain.SS.164 <- set.burnin(fixed.mcmc.SS.164$load(), 0.3)
#out.SS.164 <- summary(chain.SS.164)
#chain.SS.165 <- set.burnin(fixed.mcmc.SS.165$load(), 0.3)
#out.SS.165 <- summary(chain.SS.165)
#chain.SS.166 <- set.burnin(fixed.mcmc.SS.166$load(), 0.3)
#out.SS.166 <- summary(chain.SS.166)
#chain.SS.167 <- set.burnin(fixed.mcmc.SS.167$load(), 0.3)
#out.SS.167 <- summary(chain.SS.167)
#chain.SS.168 <- set.burnin(fixed.mcmc.SS.168$load(), 0.3)
#out.SS.168 <- summary(chain.SS.168)
#chain.SS.169 <- set.burnin(fixed.mcmc.SS.169$load(), 0.3)
#out.SS.169 <- summary(chain.SS.169)
#chain.SS.170 <- set.burnin(fixed.mcmc.SS.170$load(), 0.3)
#out.SS.170 <- summary(chain.SS.170)


#pdf("OSTfixedSS.160s_trace_plots.pdf") 
# Creating a plot
#plot(chain.SS.161)
#plot(chain.SS.162)
#plot(chain.SS.163)
#plot(chain.SS.164)
#plot(chain.SS.165)
#plot(chain.SS.166)
#plot(chain.SS.167)
#plot(chain.SS.168)
#plot(chain.SS.169)
#plot(chain.SS.170)
# Closing the graphical device
#dev.off() 


chain.SS.161 <- fixed.mcmc.SS.161$load()
chain.SS.162 <- fixed.mcmc.SS.162$load()
chain.SS.163 <- fixed.mcmc.SS.163$load()
chain.SS.164 <- fixed.mcmc.SS.164$load()
chain.SS.165 <- fixed.mcmc.SS.165$load()
chain.SS.166 <- fixed.mcmc.SS.166$load()
chain.SS.167 <- fixed.mcmc.SS.167$load()
chain.SS.168 <- fixed.mcmc.SS.168$load()
chain.SS.169 <- fixed.mcmc.SS.169$load()
chain.SS.170 <- fixed.mcmc.SS.170$load()


## All 10 chains combined
SS.chains <- list(chain.SS.161,chain.SS.162,chain.SS.163,chain.SS.164,chain.SS.165,chain.SS.166,chain.SS.167,chain.SS.168,chain.SS.169,chain.SS.170)
SS.160s.comb <- combine.chains(SS.chains, burnin.prop = 0.3)

print("OSTfixedSS.160s.comb")
#summary(SS.160s.comb)
#pdf("OSTfixedSS.160s_comb_trace_plot.pdf")
#plot(SS.160s.comb)
#dev.off()

## summary output of comb chain:
#bayou MCMC chain: 2000001 generations
#140065 samples, first 0 samples discarded as burnin
#
#
#Summary statistics for parameters:
#                            Mean           SD     Naive SE Time-series SE
#lnL                 4.657582e+02 1.539310e+00 4.113023e-03   5.221548e-03
#prior               4.046540e+00 7.711508e-01 2.060509e-03   2.088345e-03
#alpha               6.853875e-04 5.616382e-04 1.500693e-06   2.058816e-06
#sig2                4.416138e-04 3.034089e-05 8.107061e-08   1.221601e-07
#k                   1.000000e+00 0.000000e+00 0.000000e+00   0.000000e+00
#ntheta              2.000000e+00 0.000000e+00 0.000000e+00   0.000000e+00
#root.theta         -1.846249e+00 1.590350e-01 4.249404e-04   4.358038e-04
#root.beta_BodyMass  5.091518e-01 7.384422e-03 1.973112e-05   2.470480e-05
#all theta          -1.701225e+00 3.348171e-01           NA             NA
#all beta_BodyMass   6.141180e-01 1.113882e-01           NA             NA
#                   Effective Size    HPD95Lower    HPD95Upper
#lnL                      86906.74  4.627134e+02  4.680600e+02
#prior                   136355.99  2.522520e+00  4.959985e+00
#alpha                    74418.00  2.792187e-08  1.785116e-03
#sig2                     61687.62  3.845335e-04  5.022915e-04
#k                            0.00  1.000000e+00  1.000000e+00
#ntheta                       0.00  2.000000e+00  2.000000e+00
#root.theta              133169.15 -2.158376e+00 -1.534531e+00
#root.beta_BodyMass       89345.00  4.947281e-01  5.236921e-01
#all theta                      NA            NA            NA
#all beta_BodyMass              NA            NA            NA
#
#
#Branches with posterior probabilities higher than 0.1:
#     pp magnitude.of.theta2 naive.SE.of.theta2 rel.location
#1730  1             -1.5562         0.00105822    0.4956748


## Adding stepping stone to estimate marginal likelihoods and compare bayes factors for model comparison
## Reference code from Uyeda et al
require(foreach)
require(doParallel)
registerDoParallel(cores=10)

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

currentTime <- Sys.time()
print(currentTime)

ss.SS.160s<-fixed.mcmc.SS.161$steppingstone(500000, SS.160s.comb, Bk = Bk, burnin=0.3, plot=FALSE)
print("OSTfixedss.SS.160s")
ss.SS.160s$lnr
## lnr == 449.2375
pdf("OSTfixedSS.160s.comb_ss_plot.pdf")
plot(ss.SS.160s)
dev.off()
print("OSTfixedss.SS.160s")
ss.SS.160s$lnr
currentTime <- Sys.time()
print(currentTime)
