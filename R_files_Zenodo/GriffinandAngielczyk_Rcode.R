#install.packages("phytools")
library(phytools)
#install.packages("geiger")
library(geiger)
#install.packages("phylolm")
library(phylolm) 
#install.packages("nlme")
library(nlme)
#install.packages("caper")
library(caper)
#install.packages("bayou")
library(bayou)
getwd()
 
conc.tree<-read.tree("Olroydetal_2018_timecalibrated.tre") #read in consensus tree from Olroyd et al. (2018)
#plotTree(conc.tree, color="black") #plot tree

##ANCESTRAL STATE RECONSTRUCTION FOR SACRUM COUNT##

sacrumtip<-c("Galepus", "Colobodectes", "Rhachiocephalus", "Kitchinganomodon", 
             "Keyseria", "L_murrayi","L_curvatus", "Daptocephalus", "Dinanomodon", "Brachyprosopus") #Terminal taxa without sacrum data

sactree<-drop.tip(conc.tree, sacrumtip) #drop taxa with missing data
plotTree(sactree, color="black")#plot tree


sac<-read.csv("dicynodont_sacrum_states.csv", row.names = 1) #read in sacrum data
sac2<-as.matrix(sac)[,1] #transform into matrix


##Fit models for sacrum evolution

meristicARD.matrix<-matrix(c(0, 1, 0, 0, 2, 0, 3, 
                             0, 0, 4, 0, 5, 0, 0, 6, 0),nrow=4) #make index rates matrix for meristic ARD model

meristicARD.fit<-fitDiscrete(sactree, sac2, model=meristicARD.matrix); meristicARD.fit #meristic RAD model
meristic.fit<-fitDiscrete(sactree, sac2, model="meristic"); meristic.fit #standard meristic model
eqrates.fit<-fitDiscrete(sactree, sac2, model="ER"); eqrates.fit #equal-rates model
symmet.fit<-fitDiscrete(sactree, sac2, model="SYM"); symmet.fit #symmetrical-rates model
diffrat.fit<-fitDiscrete(sactree, sac2, model="ARD"); diffrat.fit #All-rates-different model

##AICc and other model support metrics
aicdisc<-setNames(
  c(meristicARD.fit$opt$aicc, meristic.fit$opt$aicc, eqrates.fit$opt$aicc, symmet.fit$opt$aicc, diffrat.fit$opt$aicc),
  c("meristic ARD","meristic","ER", "SYM", "ARD"))
aicw(aicdisc) 

##Plot rate parameters for transitions for all models 
plot(meristicARD.fit, show.zeros=F, cex.rates=0.7, 
     main= "Meristic all-rates-different model")
plot(meristic.fit, show.zeros=F, cex.rates=0.7,
     main= "Meristic model")
plot(eqrates.fit, show.zeros=F, cex.rates=0.7, 
     main= "Equal-rates model")
plot(symmet.fit, show.zeros=F, cex.rates=0.7, 
     main= "Symmetrical rates model")
plot(diffrat.fit, show.zeros=F, cex.rates=0.7, 
     main= "All-rates-different model")

##Ancestral state reconstruction

mer.matrix<-matrix(c(0, 1, 0, 0, 1, 0, 2, 
                     0, 0, 2, 0, 3, 0, 0, 3, 0),nrow=4) #make index rates matrix for a standard meristic model

fitmer.ARD<-ace(sac2,sactree,model=meristicARD.matrix, type="discrete"); fitmer.ARD #meristic ARD model
fitmer<-ace(sac2,sactree,model=mer.matrix, type="discrete"); fitmer #standard meristic model
fitER<-ace(sac2,sactree,model="ER", type="discrete"); fitER #equal-rates model
fitSYM<-ace(sac2, sactree, model="SYM", type="discrete"); fitSYM #symmetrical-rates model
fitARD<-ace(sac2,sactree,model="ARD", type="discrete"); fitARD #all-rates-different model

#Plots of ancestral states
cols<-setNames(palette(c("lightblue1","blue",
                         "darkorange1","firebrick2"))[1:4],sort(unique(sac2))); cols

#Meristic model
plotTree(sactree, color="black", fsize=1,ftype="i")
nodelabels(node=1:sactree$Nnode+Ntip(sactree),pie=fitmer$lik.anc,piecol=cols,cex=0.75)
add.simmap.legend(leg = c("3 Sacrals",  "4 Sacrals", "5 Sacrals", "6+ Sacrals"), 
                  colors=cols, prompt = F, x = 40, y = 30, shape = "circle")
title("Meristic Model", line = -2)

#Meristic ARD model
plotTree(sactree, color="black", fsize=1,ftype="i")
nodelabels(node=1:sactree$Nnode+Ntip(sactree),pie=fitmer.ARD$lik.anc,piecol=cols,cex=0.75)
add.simmap.legend(leg = c("3 Sacrals",  "4 Sacrals", "5 Sacrals", "6+ Sacrals"), 
                  colors=cols, prompt = F, x = 40, y = 30, shape = "circle")
title("Meristic All-Rates-Different Model", line = -2)

#ER model
plotTree(sactree, color="black", fsize=1,ftype="i")
nodelabels(node=1:sactree$Nnode+Ntip(sactree),pie=fitER$lik.anc,piecol=cols,cex=0.75)
add.simmap.legend(leg = c("3 Sacrals",  "4 Sacrals", "5 Sacrals", "6+ Sacrals"), 
                  colors=cols, prompt = F, x = 40, y = 30, shape = "circle")
title("Equal-Rates Model", line = -2)

#SYM model
plotTree(sactree, color="black", fsize=1,ftype="i")
nodelabels(node=1:sactree$Nnode+Ntip(sactree),pie=fitSYM$lik.anc,piecol=cols,cex=0.75)
add.simmap.legend(leg = c("3 Sacrals",  "4 Sacrals", "5 Sacrals", "6+ Sacrals"), 
                  colors=cols, prompt = F, x = 40, y = 30, shape = "circle")
title("Symmetrical-Rates Model", line = -2)

#ARD model
plotTree(sactree, color="black", fsize=1,ftype="i")
nodelabels(node=1:sactree$Nnode+Ntip(sactree),pie=fitARD$lik.anc,piecol=cols,cex=0.75)
add.simmap.legend(leg = c("3 Sacrals",  "4 Sacrals", "5 Sacrals", "6+ Sacrals"), 
                  colors=cols, prompt = F, x = 40, y = 30, shape = "circle")
title("All-Rates-Different Model", line = -2)


##ANCESTRAL STATE RECONSTRUCTION FOR BODY SIZE##

bsl<-read.csv("dicynodont_sizes.csv", row.names = 1) #read in basal skull length data

bsl2<-as.matrix(bsl)[,1]; bsl2 #turn into matrix

#Model support for BM vs OU models
BMfit<-fitContinuous(conc.tree, bsl2, model ="BM"); BMfit
OUfit<-fitContinuous(conc.tree, bsl2, model ="OU"); OUfit

aiccont<-setNames(
  c(BMfit$opt$aicc, OUfit$opt$aicc),
  c("BM","OU"))
aicw(aiccont) 

#Ancestral state reconstruction with BM model
bsfit<-anc.ML(conc.tree,bsl2); bsfit #run ancestral state reconstruction under Brownian motion

?anc.ML
plot(conc.tree); nodelabels(frame="circle", cex=.70)

#Plot ancestral states for basal skull length
obj<-contMap(conc.tree, bsl2, res=1000, lwd=6, legend=0.5*max(nodeHeights(conc.tree)), 
             ftype=c("bi","reg"),fsize=c(0.9,0.9), outline=T, plot=F) #start plot

plot(obj,legend=0.7*max(nodeHeights(conc.tree)),fsize=c(0.9,0.9), lwd=5, col = "black") #plot ancestral states


##CORRELATION BETWEEN BODY SIZE AND SACRAL COUNT##

bothvars<-c("Galepus", "Colobodectes", "Rhachiocephalus", "Kitchinganomodon", 
            "Keyseria", "L_murrayi","L_curvatus", "Daptocephalus", "Dinanomodon","Brachyprosopus") #Terminal taxa lacking either basal skull length data or sacrum count data

bothtree<-drop.tip(conc.tree, bothvars) #drop taxa with missing data
plotTree(bothtree, color="black") #plot reduced tree

sizes<-read.csv("dicynodont_sizes.csv") #read in size data
size<-sizes[-37:-46,] #drop Brachyprosopus, which doesn't have sacrum data
  
sac<-read.csv("dicynodont_sacrum_states.csv") #read in sacrum states

sac$binary<-as.numeric(sac$state>1) #Add column of binary sacrum data ("0" is < 5, "1" is => 4)

bothdata<-data.frame(size, sac$state, sac$binary, row.names = 1) #make data frame sacrum states

##Phylogenetic signal D-statistic for binary trait
dicynodont.cd <- comparative.data(phy = bothtree, 
                                  data = sac,  names.col = "species.sac", vcv = TRUE, 
                                  na.omit = FALSE, warn.dropped = TRUE) #convert data into comparative data object

sacrum.D <- phylo.d(data = dicynodont.cd, 
                    phy = bothtree, binvar = binary, permut = 1000); sacrum.D #Find D-statistic for phylogenetic signal

#D-statistic for just Bidentalia
bident<-c("Suminia", "Patranomodon", "Galepus", "Eodicynodon", "Colobodectes", "Diictodon", 
          "Robertia", "Eosimops", "Pristerodon", "Brachyprosopus", "Endothiodon", "Niassodon", 
          "Dicynodontoides", "Emydops", "Myosaurus", "Cistecephalus", "Rhachiocephalus", "Kitchinganomodon", "Keyseria", 
          "L_murrayi","L_curvatus", "Daptocephalus", "Dinanomodon") #Non-bidentalians, or bidentalians lacking either basal skull length data or sacrum count data

bident.tree <- drop.tip(conc.tree, bident) #Make Bidentalia-only tree

bident.cd <- comparative.data(phy = bident.tree, 
                              data = sac,  names.col = "species.sac", vcv = TRUE, 
                              na.omit = FALSE, warn.dropped = TRUE) #convert data into comparative data object
bident.D <- sacrum.D <- phylo.d(data = bident.cd, 
                                phy = bident.tree, binvar = state, permut = 1000); bident.D #Find D-statistic for phylogenetic signal

##Regression analyses

#Phylogenetic logistic regression
plog<-phyloglm(bothdata$sac.binary ~ bothdata$bsl,
             data=bothdata, phy=bothtree, btol=500); plog #run phylogenetic logistic regression on size vs. binary sacrum state
summary(plog)
plog$aic

#Standard logistic regression
reglog<-glm(bothdata$sac.binary ~ bothdata$bsl,
            data=bothdata, family = binomial(link = logit)); reglog
summary(reglog)


#Plot phylogenetic logistic regression
par(mar=c(4,4,4,4))
plot(size$bsl,jitter(sac$binary,factor=0,amount=0.02), axes=F, 
     xlab="Basal Skull Length (mm)",ylab="response",cex=1.75, pch=21,
     bg = ifelse(bothdata$sac.state <= 0, "lightblue1", 
                 ifelse(bothdata$sac.state <= 1, "blue", 
                        ifelse(bothdata$sac.state <= 2, "darkorange1", "firebrick2"))))
axis(2)
axis(1,at=seq(0, 700, by=50))
box()
coeff.plog <- coef(plog) # make log curve a plottable curve
curve(plogis(coeff.plog[1]+coeff.plog[2]*x),add=TRUE, lwd=2) #add phylogenetic logistic regression curve to plot

coeff.reglog <- coef(reglog) # make log curve a plottable curve
curve(plogis(coeff.reglog[1]+coeff.reglog[2]*x),add=TRUE, lwd=2, lty = "dashed") #add standard logistic regression curve to plot



#Other regressions
phylo.pois<-phyloglm(bothdata$sac.state ~ bothdata$bsl, data= bothdata, 
                     phy=bothtree, method = "poisson_GEE"); phylo.pois #run phylogenetic poisson regression on size vs. sacrum state
summary(phylo.pois)

pois<-glm(bothdata$sac.state ~ bothdata$bsl, 
          data= bothdata, family = "poisson"); pois #run poisson regression on size vs. sacrum state
summary(pois)


genlm<-gls(sac.state ~ bsl, data=bothdata, 
           correlation=corBrownian(phy=bothtree)); genlm #run phylogenetic generalized least squares regresson on size vs. sacrum state
summary(genlm)


#Plot other regressions
par(mar=c(4,4,4,4))
plot(bothdata$bsl, jitter(bothdata$sac.state,factor=0,amount=0.02), axes=F, 
     xlab="Basal Skull Length (mm)",ylab="Number of Sacral Vertebrae",cex=1.75, pch=21,
     bg = ifelse(bothdata$sac.state <= 0, "lightblue1", 
                 ifelse(bothdata$sac.state <= 1, "blue", 
                        ifelse(bothdata$sac.state <= 2, "darkorange1", "firebrick2"))))
axis(2)
axis(1,at=seq(0, 700, by=50))
box()

coeff <- coef(pois)
xsort <- sort(bothdata$bsl)
logpois=coeff[1]+coeff[2]*xsort
lines(xsort, exp(logpois), lwd = 2,lty = "dashed") #Add normal Poisson regression line

coeff2 <- coef(phylo.pois)
log.phypois=coeff2[1]+coeff2[2]*xsort
lines(xsort, exp(log.phypois), lwd = 2) #Add phylogenetic Poisson regression line

coeff3<- coef(genlm)
phygls=coeff3[1]+coeff3[2]*xsort
lines(xsort, phygls, lwd = 1.5, lty = "dotted") #Add PGLS line

### Combined plot
par(mfrow=c(2,1))

par(mar=c(4,4,4,4))
plot(size$bsl,jitter(sac$binary,factor=0,amount=0.02), axes=F, 
     xlab="Basal Skull Length (mm)",ylab="response",cex=1.75, pch=21,
     bg = ifelse(bothdata$sac.state <= 0, "lightblue1", 
                 ifelse(bothdata$sac.state <= 1, "blue", 
                        ifelse(bothdata$sac.state <= 2, "darkorange1", "firebrick2"))))
axis(2)
axis(1,at=seq(0, 700, by=50))
box()
coeff.plog <- coef(plog) # make log curve a plottable curve
curve(plogis(coeff.plog[1]+coeff.plog[2]*x),add=TRUE, lwd=2) #add phylogenetic logistic regression curve to plot

coeff.reglog <- coef(reglog) # make log curve a plottable curve
curve(plogis(coeff.reglog[1]+coeff.reglog[2]*x),add=TRUE, lwd=2, lty = "dashed") #add phylogenetic logistic regression curve to plot

par(mar=c(4,4,4,4))
plot(bothdata$bsl, jitter(bothdata$sac.state,factor=0,amount=0.02), axes=F, 
     xlab="Basal Skull Length (mm)",ylab="Number of Sacral Vertebrae",cex=1.75, pch=21,
     bg = ifelse(bothdata$sac.state <= 0, "lightblue1", 
                 ifelse(bothdata$sac.state <= 1, "blue", 
                        ifelse(bothdata$sac.state <= 2, "darkorange1", "firebrick2"))))
axis(2)
axis(1,at=seq(0, 700, by=50))
box()

coeff <- coef(pois)
xsort <- sort(bothdata$bsl)
logpois=coeff[1]+coeff[2]*xsort
lines(xsort, exp(logpois), lwd = 2, lty = "dashed") #Add normal Poisson regression line

coeff2 <- coef(phylo.pois)
log.phypois=coeff2[1]+coeff2[2]*xsort
lines(xsort, exp(log.phypois), lwd = 2) #Add phylogenetic Poisson regression line

coeff3<- coef(genlm)
phygls=coeff3[1]+coeff3[2]*xsort
lines(xsort, phygls, lwd = 1.5, lty = "dotted") #Add PGLS line

##FITTING A MULTI-OPTIMA ORNSTEIN-UHLENBECK MODEL TO DATA TO DETECT SHIFTS IN BODY SIZE EVOLUTION##

MEvar <- 0.01 #Make measurement error factor for a standard error prior

priorOU <- make.prior(conc.tree, 
                      dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", 
                                 dk="cdpois", dtheta="dnorm"),
                      param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                 dk=list(lambda=3, kmax=10), dsb=list(bmax=1, prob=1), 
                                 dtheta=list(mean=mean(bsl2), sd=1.5*sd(bsl2)))
) #Set the prior probability for body size as 1.5 times the standard deviation of the measured body sizes

par(mfrow=c(1,1))
startpars <- priorSim(priorOU, conc.tree, plot=TRUE)$pars[[1]] #run this a few different times to see different possible starting points

priorOU(startpars) #run MCMC chain with multiple starting values taken from our distribution of priors 


set.seed(1)

## Replace "new.dir=TRUE" with the directory you want files to be put if you would like to keep these files. If new.dir=TRUE, then it's put in a temporary directory that will be deleted by your computer eventually
mcmcOU <- bayou.makeMCMC(conc.tree, bsl2, SE=MEvar, prior=priorOU, 
                         new.dir=F, outname="modelOU_r001", plot.freq=NULL) # Set up the MCMC
mcmcOU$run(1000000) # Run the MCMC. We ran for 1,000,000 generations to ensure that we get convergence

chainOU <- mcmcOU$load(saveRDS = TRUE) #Load back MCMC result files into R, save as a *.rds file in working directory


#chainOU<-readRDS(file="modelOU_r001.chain.rds") #if necessary, read *.rds file back into R (e.g., if session was interrupted)

chainOU <- set.burnin(chainOU, 0.03) #set a "burn-in" parameter that tells the package coda to discard the first 0.3 entries of the chain
summary(chainOU) #View summary of results

par(mfrow=c(2,2))
plot(chainOU, auto.layout=F) #plot results, the traces look good


par(mfrow=c(1,1))
plotSimmap.mcmc(chainOU, burnin = 0.03, show.tip.label=T, lwd = 2, pp.cutoff = 0.3, edge.type = "regimes") #Diameters of circles at nodes are proportional to the posterior probability of shifts of phenotypic optima at those locations. 
phenogram.density(conc.tree, bsl2, burnin = 0.03, chainOU, pp.cutoff = 0.3) #phenogram of basal skull length vs. time recovered by bayou, with ddensity of phenotypic optima


#Now we will run another MCMC chain from an independent starting position and compare with our previous chain to make sure we have convergence in our parameters 

startpars2 <- priorSim(priorOU, conc.tree, plot=TRUE)$pars[[1]] #run this a few different times to see different starting values

priorOU(startpars2) #run MCMC chain with multiple starting values taken from our distribution of priors 
set.seed(1)
mcmcOU2 <- bayou.makeMCMC(conc.tree, bsl2, SE=MEvar, prior=priorOU, 
                          new.dir=F, outname="modelOU2_r002", plot.freq=NULL) # Set up the new MCMC
mcmcOU2$run(1000000) # Run the new MCMC

chainOU2 <- mcmcOU2$load(saveRDS = TRUE) #Load back new MCMC result files into R

#chainOU2<-readRDS(file="modelOU2_r002.chain.rds") #if necessary, read *.rds file back into R (e.g., if session was interrupted)
chainOU2 <- set.burnin(chainOU2, 0.03) #set burn-in

#Compare branch posterior probabilities of the two chains. If converged, they should fall along the y=x line
L1 <- Lposterior(chainOU,conc.tree, burnin=0.03)
L2 <- Lposterior(chainOU2,conc.tree, burnin=0.03)
plot(L1$pp,L2$pp, xlab="Chain 1", ylab="Chain 2")
curve(1*x, add=TRUE, lty=2)

