#############################################################################################
#	Brain Mass Analysis for Schumacher and Carlson 2022										#
#																							#
#                Combines new otophysan and osteoglossiform Sukhum et al 2018 data          #
#				 with published actinopterygian data from Tsuboi et al 2018					#
#				 and Tsuboi 2021															#
#																							#
#				 GG Model																	#
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

BrainMass <- setNames(tel_avg_data$BrainMass, rownames(tel_avg_data))
BrainSE <- setNames(tel_se_data$BrainMass, rownames(tel_se_data))

## "This is a good reminder to *always try to use measurement error in your analyses*. OU models especially are affected by measurement error. 
## This is because OU models have the effect of "erasing" evolutionary history with increasing *alpha*. If you don't account for measurement error, 
## then that measurement error will be transferred to the evolutionary process model. You can make a Brownian Motion model look very OU like if there is a lot of measurement error."

## Defining the priors: "prior function is going to take our parameters and output the *prior probability* of our parameter values. It represents our initial degree of belief in what values the parameters will take."
## "*make.prior* tries to be reasonable and not make you type everything out, but **do not be comfortable with defaults**. One trick to make sure your prior functions are reasonable is to simulate a bunch of values 
## and take the quantiles of the distribution. We are using a half-Cauchy distribution for *alpha* and *sig2*, which is a good weakly informative prior for scale parameters."

## theta = mean trait value estimate
## alpha = rate at which theta changes
## sigma = variance in trait value over time

## Run models w 1. global intercept and slopes (GG), 2. separate intercept and global slope (SG), and 3. separate intercepts and slopes (SS)
## And use model comparison to see which has a better fit (likely the third model which allows both slope and intercept to vary)
## Defined priors following Uyeda et al 2017

## Using .599 as slope prior mean and .2 as sd; .599 = mean of 8 reported tsuboi slopes, sd of tsuboi slope = .18 -> raised to .2

## For runs 51-70, I updated the intercept (theta) prior to have a mean of -2 instead of 0 since a positive intercept value would in no way be expected here. -2 was arbitrarily chosen based on the plot of brain v body; sd remains at a wide value of 1.
prior.GG <- make.prior(tel_tree_rab, plot.prior = FALSE, 
                       dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dbeta_BodyMass="dnorm",
                                  dsb="fixed", dk="fixed", dtheta="dnorm"), 
                       param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                  dbeta_BodyMass=list(mean=0.599, sd=0.2),
                                  dtheta=list(mean=-1.682, sd=.4)),
                       fixed=list(k=0, sb=numeric(0))
)

## Set tuning parameters:
## "There is a bit of art to tuning the parameters, which may require making multiple runs and trying to get the acceptance ratios in the right 
## region (0.2-0.4). But these should work well for these models and data. If the acceptance ratio for a certain parameter is too high, increase 
## the tuning parameter for that variable. If the acceptance ratio is too low, decrease it. The scale of the regression coefficient, for example, 
## should give you some idea of what these parameters should be."

DGG.67 = list(alpha=7, sig2=.5, beta_BodyMass=0.1, k=1, theta=2, slide=1)


## Building the MCMC models
set.seed(671)
model.GG.671 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.671 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.671$model, prior=prior.GG, startpar=model.GG.671$startpar, file.dir=getwd(), outname="modelGG_r671", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(671)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.671$run(2000000)		

## Building the MCMC models
set.seed(672)
model.GG.672 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.672 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.672$model, prior=prior.GG, startpar=model.GG.672$startpar, file.dir=getwd(), outname="modelGG_r672", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(672)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.672$run(2000000)		

## Building the MCMC models
set.seed(673)
model.GG.673 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.673 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.673$model, prior=prior.GG, startpar=model.GG.673$startpar, file.dir=getwd(), outname="modelGG_r673", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(673)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.673$run(2000000)		


## Building the MCMC models
set.seed(674)
model.GG.674 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.674 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.674$model, prior=prior.GG, startpar=model.GG.674$startpar, file.dir=getwd(), outname="modelGG_r674", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(674)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.674$run(2000000)		


## Building the MCMC models
set.seed(675)
model.GG.675 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.675 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.675$model, prior=prior.GG, startpar=model.GG.675$startpar, file.dir=getwd(), outname="modelGG_r675", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(675)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.675$run(2000000)		


## Building the MCMC models
set.seed(676)
model.GG.676 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.676 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.676$model, prior=prior.GG, startpar=model.GG.676$startpar, file.dir=getwd(), outname="modelGG_r676", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(676)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.676$run(2000000)		

## Building the MCMC models
set.seed(677)
model.GG.677 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.677 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.677$model, prior=prior.GG, startpar=model.GG.677$startpar, file.dir=getwd(), outname="modelGG_r677", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(677)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.677$run(2000000)		

## Building the MCMC models
set.seed(678)
model.GG.678 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.678 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.678$model, prior=prior.GG, startpar=model.GG.678$startpar, file.dir=getwd(), outname="modelGG_r678", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(678)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.678$run(2000000)		

## Building the MCMC models
set.seed(679)
model.GG.679 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.679 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.679$model, prior=prior.GG, startpar=model.GG.679$startpar, file.dir=getwd(), outname="modelGG_r679", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(679)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.679$run(2000000)		

## Building the MCMC models
set.seed(680)
model.GG.680 <- makeBayouModel(BrainMass ~ BodyMass, rjpars = c(), 
                           tree=tel_tree_rab, dat=BrainMass, pred=tel_avg_data, SE=BrainSE, prior=prior.GG, D=DGG.67)
## Making the MCMC objects
mcmc.GG.680 <- bayou.makeMCMC(tel_tree_rab, BrainMass, pred=tel_avg_data, SE=BrainSE, model=model.GG.680$model, prior=prior.GG, startpar=model.GG.680$startpar, file.dir=getwd(), outname="modelGG_r680", plot.freq=NULL)
## Running MCMC models:
## Ksepka ran 10 parallel chains @ 2 mil iterations, burn-in proportion = .3; sampling freq from chain??
set.seed(680)		## sets seed for rng, change btwn runs, but keep reproducable
mcmc.GG.680$run(2000000)		


## Chain results
#chain.GG.671 <- set.burnin(mcmc.GG.671$load(), 0.3)
#chain.GG.672 <- set.burnin(mcmc.GG.672$load(), 0.3)
#chain.GG.673 <- set.burnin(mcmc.GG.673$load(), 0.3)
#chain.GG.674 <- set.burnin(mcmc.GG.674$load(), 0.3)
#chain.GG.675 <- set.burnin(mcmc.GG.675$load(), 0.3)
#chain.GG.676 <- set.burnin(mcmc.GG.676$load(), 0.3)
#chain.GG.677 <- set.burnin(mcmc.GG.677$load(), 0.3)
#chain.GG.678 <- set.burnin(mcmc.GG.678$load(), 0.3)
#chain.GG.679 <- set.burnin(mcmc.GG.679$load(), 0.3)
#chain.GG.680 <- set.burnin(mcmc.GG.680$load(), 0.3)

chain.GG.671 <- mcmc.GG.671$load()
chain.GG.672 <- mcmc.GG.672$load()
chain.GG.673 <- mcmc.GG.673$load()
chain.GG.674 <- mcmc.GG.674$load()
chain.GG.675 <- mcmc.GG.675$load()
chain.GG.676 <- mcmc.GG.676$load()
chain.GG.677 <- mcmc.GG.677$load()
chain.GG.678 <- mcmc.GG.678$load()
chain.GG.679 <- mcmc.GG.679$load()
chain.GG.680 <- mcmc.GG.680$load()


#out.GG.671 <- summary(chain.GG.671)
#out.GG.672 <- summary(chain.GG.672)
#out.GG.673 <- summary(chain.GG.673)
#out.GG.674 <- summary(chain.GG.674)
#out.GG.675 <- summary(chain.GG.675)
#out.GG.676 <- summary(chain.GG.676)
#out.GG.677 <- summary(chain.GG.677)
#out.GG.678 <- summary(chain.GG.678)
#out.GG.679 <- summary(chain.GG.679)
#out.GG.680 <- summary(chain.GG.680)

#pdf("GG.670s_trace_plots_Mar8.pdf") 
# Creating a plot
#plot(chain.GG.671)
#plot(chain.GG.672)
#plot(chain.GG.673)
#plot(chain.GG.674)
#plot(chain.GG.675)
#plot(chain.GG.676)
#plot(chain.GG.677)
#plot(chain.GG.678)
#plot(chain.GG.679)
#plot(chain.GG.680)
# Closing the graphical device
#dev.off() 


## All 10 chains combined
GG.chains <- list(chain.GG.671,chain.GG.672,chain.GG.673,chain.GG.674,chain.GG.675,chain.GG.676,chain.GG.677,chain.GG.678,chain.GG.679,chain.GG.680)
GG.670s.comb <- combine.chains(GG.chains, burnin.prop = 0.3)

#print("GG.670s.comb")
#summary(GG.670s.comb)
#pdf("GG.670s_comb_trace_plot.pdf")
#plot(GG.670s.comb)
#dev.off()

## output of comb chain:
#bayou MCMC chain: 2000001 generations
#1400180 samples, first 0 samples discarded as burnin
#
#
#Summary statistics for parameters:
#                       Mean           SD     Naive SE Time-series SE
#lnL            4.579242e+02 1.425718e+00 1.204874e-03   2.628556e-03
#prior          4.147576e+00 1.881440e-01 1.590005e-04   1.901848e-04
#alpha          5.069351e-04 4.355569e-04 3.680891e-07   8.164631e-07
#sig2           4.444677e-04 2.937021e-05 2.482076e-08   6.480755e-08
#beta_BodyMass  5.122141e-01 7.352377e-03 6.213493e-06   1.215749e-05
#theta         -1.824832e+00 1.635987e-01 1.382572e-04   1.622981e-04
#k              0.000000e+00 0.000000e+00 0.000000e+00   0.000000e+00
#ntheta         1.000000e+00 0.000000e+00 0.000000e+00   0.000000e+00
#              Effective Size    HPD95Lower    HPD95Upper
#lnL                 294193.2  4.551225e+02  4.599798e+02
#prior               978653.3  3.761217e+00  4.327832e+00
#alpha               284588.0  5.893848e-10  1.373826e-03
#sig2                205382.2  3.879973e-04  5.024976e-04
#beta_BodyMass       365736.2  4.980938e-01  5.268844e-01
#theta              1016091.5 -2.145273e+00 -1.502743e+00
#k                        0.0  0.000000e+00  0.000000e+00
#ntheta                   0.0  1.000000e+00  1.000000e+00
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
registerDoParallel(cores=5)

Bk <- qbeta(seq(0,1, length.out=50), 0.3,1)

currentTime <- Sys.time()
print(currentTime)

ss.GG.670s<-mcmc.GG.671$steppingstone(2000000, GG.670s.comb, Bk = Bk, burnin=0.3, plot=FALSE)
pdf("GG.670s.comb_ss_plot.pdf")
plot(ss.GG.670s)
dev.off()
print("ss.GG.670s")
ss.GG.670s$lnr
## ss.GG.670s$lnr = 442.5746
currentTime <- Sys.time()
print(currentTime)
