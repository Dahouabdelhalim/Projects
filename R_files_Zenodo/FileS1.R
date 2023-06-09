# Code for Phylogenetic evidence from freshwater crayfishes that cave adaptation is not an evolutionary dead-end
# Simplified to include only reported results 
# Much of the code is adapted from :
# Smith, Stacey D., and Emma E. Goldberg. "Tempo and mode of flower color evolution." American Journal of Botany 102.7 (2015): 1014-1025.
# AND
# Hamm, Christopher A., and James A. Fordyce. "Greater host breadth still not associated with increased diversification rate in the Nymphalidaeâ€”A response to Janz et al." Evolution 70.5 (2016): 1156-1160.

library(diversitree)
library(parallel) 
library(phytools)

##########
#  Read in trees and data
##########
# setwd('<Data_folder>')
originaltrees <- read.tree("DatasetS4.tre")
taxa <- originaltrees[[1]]$tip.label

# Read in habitat data derived from TableS2.csv
temp <- read.csv("FileS3.habitats_cave.geosse.csv", header=FALSE, as.is=TRUE)
states <- structure(temp[[2]], names=temp[[1]])
to.drop <- setdiff(taxa, names(states))
trees <- lapply(originaltrees, drop.tip, to.drop)
# new version of ape may not think these trees are perfectly ultrametric
# if so run this function from phytools: 
#trees <- lapply(trees, force.ultrametric)

#set up GeoSSE sampling for Cambaridae subtree
# 0=both, 1=cave, 2=surface
sampling <- table(states)
sampling[1]<- 2/3
sampling[2]<- 30/43
sampling[3]<- 277/366

##########
#  Run MCMC over bootstrap trees
##########
outdir <- "mcmc_results/"
if (!file.exists(outdir))
    dir.create(outdir)
zfill <- function(x, n)
{
    nc <- nchar(x)
    zeros <- paste(rep(0, n), collapse = "")
    paste(substring(zeros, nchar(x) + 1, n), substring(x, 1, nchar(x)), sep = "")
}

run.mcmc <- function(treenum)
{
    tree <- trees[[treenum]]

	lik.full <- make.geosse(tree, states, sampling.f=sampling)
    lik.names <- argnames(lik.full)
    #get good starting point with ML
    p <- starting.point.geosse(tree)
    ml1 <- find.mle(lik.full, p)
    p <- coef(ml1)

    #small run to get step sizes
    prior <- make.prior.exponential(1/2)
    set.seed(1)
    tmp <- mcmc(lik.full, p, nsteps=250, prior=prior, w=1, print.every=10)
    w <- diff(sapply(tmp[2:8], quantile, c(0.025, 0.975)))

    message(paste("Done with prelim for tree", treenum))

    # real chain
    ans <- mcmc(lik.full, p, nsteps=5000, prior=prior, w=w, print.every=100)

    # write the results for this tree
    ans <- data.frame(treeid=treenum, ans)
    outfile <- paste(outdir, "geosse_tree", zfill(treenum, 4), ".csv", sep="")
    write.csv(ans, file=outfile, row.names=F)

    message(paste("Done fitting tree", treenum))
}

# Adjust subset of bootstrap trees to analyze this run
treesubset <- 1:1000
# Adjust number of cores available
junk <- mclapply(treesubset, run.mcmc, mc.cores=16)


############
# Summarize GeoSSE MCMC output
############
library(coda)
library(RSkittleBrewer)

# Function to discard first 10% of each run
burnin <- function(csv) {
    burnstart <- floor(0.1 * nrow(csv))
    postburn <- csv[burnstart:nrow(csv), ]
}
# Function to load data and burn-in each run
load_data <- function(path) { 
  files <- dir(path, pattern = '.csv', full.names = TRUE)
  tables <- lapply(files, read.csv)
  post <- lapply(tables, burnin)
  do.call(rbind, post)
}

dat <- load_data("mcmc_results")

#check for convergence
plot(dat$p ~ dat$i)
effectiveSize(dat$p)
HPDinterval(mcmc(dat))
colMeans(dat)
marginalL <- 1/mean(1/dat$p)
marginalL

### Test for evolutionary dead end (do cave lineages go extinct before speciating or leaving?)
# Is extinction greater than speciation?
dif <- with(dat, data.frame(dif=xA-sA))
pp.dif <- length(which(dif>0))/nrow(dif)
pp.dif
rat <- with(dat, data.frame(rat=xA/sA))
mean(rat$rat)
# Is extinction greater than dispersal?
dif <- with(dat, data.frame(dif=xA-dA))
pp.dif <- length(which(dif>0))/nrow(dif)
pp.dif
rat <- with(dat, data.frame(rat=xA/dA))
mean(rat$rat)
# Is extinction greater than speciation+dispersal?
dif <- with(dat, data.frame(dif=xA-(dA+sA)))
pp.dif <- length(which(dif>0))/nrow(dif)
pp.dif
rat <- with(dat, data.frame(rat=xA/(dA+sA)))
mean(rat$rat)

### Test for evolutionary dead end (are cave lineages stuck in state?)
# Is dispersal lower than speciation?
dif <- with(dat, data.frame(dif=sA-dA))
pp.dif <- length(which(dif>0))/nrow(dif)
pp.dif
rat <- with(dat, data.frame(rat=sA/dA))
mean(rat$rat)

# Test for unequal dispersal rates between cave and surface lineages
dif.disp <- with(dat, data.frame(dif.disp=dB-dA))
HPDinterval(mcmc(dif.disp$dif.disp))

# Test for unequal diversification rates
# Net diversification
net.div <- with(dat, data.frame(div.A=sA-xA, div.B=sB-xB)) 
write.csv(net.div, 'net_div.csv')
net.div.dif <- with(net.div, data.frame(dif=div.B-div.A))
pp.div <- length(which(net.div.dif>0))/nrow(net.div.dif)
pp.div
HPDinterval(mcmc(net.div.dif$dif))
ratio.div <- with(net.div, data.frame(rat.div=div.A/div.B))
write.csv(ratio.div, 'div_ratio.csv')
mean(ratio.div$rat.div)
# Speciation
dif.sp <- with(dat, data.frame(dif.sp=sB-sA))
HPDinterval(mcmc(dif.sp$dif.sp))
ratio.sp <- with(dat, data.frame(ratio.sp=sB/sA))
mean(ratio.sp$ratio.sp)
# Extinction
dif.ext <- with(dat, data.frame(dif.ext=xB-xA))
HPDinterval(mcmc(dif.ext$dif.ext))
ratio.ext <- with(dat, data.frame(ratio.ext=sB/sA))
mean(ratio.ext$ratio.ext)

########
## Marginal parameter distribution plots
########
##Rates in Cave Lineages##
cols = c("dodgerblue","darkorange","brown")
profiles.plot(c(dat[3],dat[6],dat[8]), col.line=cols, xlab="", ylab="", xlim=c(0,0.15), n.br=500)
legend("topright", c('Speciation','Extinction', 'Dispersal'), col=cols, lty=1)
title(xlab="Rate in Cave Lineages", ylab="posterior probability density", outer=T, line=-1)

##Difference in Net Diversification (Surface- Cave) and Dispersal (Surface-Cave)#
cols = c("darkslateblue","brown")
profiles.plot(c(net.div.dif, dif.disp), col.line=cols, xlab="", ylab="", n.br=500)
legend("topright", c('Net Diversification','Dispersal'), col=cols, lty=1)
title(xlab="Surface-Cave", ylab="posterior probability density", outer=T, line=-1)

##Difference in Rates(Surface-Cave)#
cols = c("dodgerblue","darkorange","brown")
profiles.plot(c(dif.sp,dif.ext,dif.disp), col.line=cols, xlab="", ylab="", xlim=c(-0.1,0.05), n.br=500)
legend("topright", c('Speciation','Extinction', 'Dispersal'), col=cols, lty=1)
title(xlab="Difference in Rates (Surface-Cave)", ylab="posterior probability density", outer=T, line=-1)

### Habitat diversification comparison plot, GeoSSE, assuming the above analysis has been run for each habitat
#cols = c("dodgerblue3","limegreen","yellow","firebrick")
#lentic <- read.csv('div_ratio.lentic.csv', row.names=1)
#lentic <- subset(lentic, rat.div>-20 & rat.div<20) 
#lotic <- read.csv('div_ratio.lotic.csv', row.names=1)
#cave <- read.csv('div_ratio.cave.csv', row.names=1)
#burrow <- read.csv('div_ratio.burrow.csv', row.names=1)
#profiles.plot(c(lentic,lotic,cave,burrow), col.line=cols, xlab="", ylab="", xlim=c(-0.5,2.5),n.br=500)
#legend("topright", c('Lentic','Lotic', 'Cave','Burrow'), col=cols, lty=1)
#title(xlab="Relative Net Diversification Rates", ylab="posterior probability density", outer=T, line=-1)

#############
#   Model Adequacy for GeoSSE
#############
library(phylometrics)
library(dplyr)
# Simulate phylogenies using the posterior samples with diversitree
phys <- c()
repeat{
	pars <- as.numeric(sample_n(dat, 1)[3:9])
	#pars <- c(pars[1:5],0,pars[6])
	t <- trees(pars, type="geosse", max.taxa=412, n=1,x0=2, include.extinct=FALSE)
	if (length(which(t[[1]]$tip.state==1))>0){
		phys <- c(phys, t)
		print(length(phys))
		}
	if (length(phys)==1000){
		break
		}
	}
# Count number of tips surviving to present in each state in simulated trees
A <- c()
B <- c()
AB <- c()
for (tree in phys) {
    tmp <- tree$tip.state
    A <- c(A, length(which(tmp==1)))
}
for (tree in phys) {
    tmp <- tree$tip.state
    B <- c(B, length(which(tmp==2)))
}
for (tree in phys) {
    tmp <- tree$tip.state
    AB <- c(AB, length(which(tmp==0)))
}
chi_dat <- data.frame(c(mean(A), mean(B), mean(AB)), c(43, 366,3))
chisq.test(chi_dat)

# Calculate SSCD for bootstrap trees considering obligate cave-dwellers as state 0 and all other state 1
stlist <- names(which(states==1))
sscd.emp <- c()
	for (tree in trees) {
	n <- treestat(tree, stlist=stlist,func=sscd)
	sscd.emp <- c(sscd.emp, n)
	}

# Backward simulation using posterior samples to keep state prevalence fixed
pars <- as.numeric(sample_n(dat, 1)[3:9])
pars <- c(pars[1:2],pars[4:7])
#pars <- c(pars[1:2],pars[4:5],0,pars[6])
phys_fixed <- treesim(pars=pars,N0=30,N1=279,sampling.f=c(30/43,279/369),max.t=500)
repeat{
	try({
	pars <- as.numeric(sample_n(dat, 1)[3:9])
	pars <- c(pars[1:2],pars[4:7])
	#pars <- c(pars[1:2],pars[4:5],0,pars[6])
	t <- treesim(pars=pars,N0=30,N1=279,sampling.f=c(30/43,279/369),max.t=500)
	phys_fixed <- c(phys_fixed, t)
	print(length(phys_fixed))
	if (length(phys_fixed)==1000){
		break
	}
	})
}
# Calculate SSCD for simulated phylogenies
sscd.sim <- c()
for (t in phys_fixed) {
    n <- treestat(t, func=sscd)
    sscd.sim <- c(sscd.sim, n)
}

#Proportion of empirical values in in 50% quantile of simulated ones
sum(sscd.emp > quantile(sscd.sim, c(0.25, 0.75))[1] & sscd.emp  < quantile(sscd.sim, c(0.25, 0.75))[2] )/1000
#Proportion of empirical values in in 95% quantile of simulated ones
sum(sscd.emp > quantile(sscd.sim, c(0.025, 0.975))[1] & sscd.emp  < quantile(sscd.sim, c(0.025, 0.975))[2] )/1000
# Student's t-test
t.test(sscd.emp, sscd.sim)

##########
#   HiSSE Model Comparison
##########
library(hisse)
outdir <- "ML_results/"
if (!file.exists(outdir))
    dir.create(outdir)

outfile1 <- paste(outdir, "Full.csv", sep="")
outfile2 <- paste(outdir, "BiSSE.csv", sep="")
outfile3 <- paste(outdir, "HiddenCave.csv", sep="")
outfile4 <- paste(outdir, "HiddenSurface.csv", sep="")
outfile5 <- paste(outdir, "CID_4.csv", sep="")
outfile6 <- paste(outdir, "CID_2.csv", sep="")
outfile7 <- paste(outdir, "GeoSSE.csv", sep="")


# Read in state data 
dat <- read.delim("FileS4.habitats_Cave.hisse.tsv")

run.MLE <- function(treenum)
{
    tree <- trees[[treenum]]

    # Set up the transition rate matrices for state-dependent models
    rate.matrix <- TransMatMaker(hidden.states = TRUE) #Full HiSSE
    Bis.rates.matrix <- TransMatMaker(hidden.states=FALSE) #Full BiSSE
    His.matrix.surface <- ParDrop(rate.matrix, c(2, 5, 7, 8, 9, 12)) #HiSSE with hidden state for surface lineages
    His.matrix.cave <- ParDrop(rate.matrix, c(3, 6, 9, 10, 11, 12)) #HiSSE with hidden state for cave lineages
    
    # Fit HiSSE models
    bisse.fit <-  hisse(tree, dat, f = c(0.71111, 0.7568), hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=Bis.rates.matrix, output.type="turnover", bounded.search=TRUE)
    His.cave <- hisse(tree, dat, f = c(0.71111, 0.7568), turnover.anc = c(1, 2, 3, 0), eps.anc = c(1, 2, 3, 0), trans.rate = His.matrix.cave, output.type = "turnover", hidden.states = TRUE, bounded.search=TRUE)
    His.surface <- hisse(tree, dat, f = c(0.71111, 0.7568), turnover.anc = c(1, 2, 0, 3), eps.anc = c(1, 2, 0, 3), trans.rate = His.matrix.surface, output.type = "turnover", hidden.states = TRUE, bounded.search=TRUE)
    His.full <- hisse(tree, dat, f = c(0.71111, 0.7568), turnover.anc = c(1, 2, 3, 4), eps.anc = c(1, 2, 3, 4), trans.rate = rate.matrix, output.type = "turnover", hidden.states = TRUE, bounded.search=TRUE)
	
	# Set up CID2
    trans.rates.hisse <- TransMatMaker(hidden.states = TRUE) #Full HiSSE
    trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))
    trans.rates.hisse[!is.na(trans.rates.hisse) & !trans.rates.hisse == 0] = 1
    
    # Fit state idependent models
    cid4.hisse <- hisse.null4(tree, dat, f= c(0.71111,0.7568), turnover.anc=rep(c(1,2,3,4),2), eps.anc=rep(c(1,2,3,4),2), trans.type="equal", bounded.search=TRUE)
    cid2.hisse <- hisse(tree, dat, f= c(0.71111,0.7568), hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.hisse, bounded.search=TRUE)
               
    # geosse
	lik.full.geosse <- make.geosse(tree, states, sampling.f=sampling)
	#get good starting point
    p <- starting.point.geosse(tree)
    geosse.fit <- find.mle(lik.full.geosse, p)
       
    aics <- data.frame(c("Full.HiSSE", "BiSSE", "Hidden.Cave", "Hidden.Surface", "CID4","CID2","GeoSSE"), c(His.full$AIC, bisse.fit$AIC, His.cave$AIC, His.surface$AIC, cid4.hisse$AIC, cid2.hisse$AIC, AIC(geosse.fit)), row.names = NULL)
	colnames(aics) <- c("model", "AIC")
	rownames(aics) <- c("Full.HiSSE", "BiSSE", "Hidden.Cave", "Hidden.Surface", "CID4","CID2", "GeoSSE")
	aics <- aics[order(aics$AIC), ]

	#Calculate AICw for each model for this tree
	for(i in 1:dim(aics)[1]){ 
	aics$delta[i] <- aics$AIC[i] - aics$AIC[1]
	} 
	aics$W <- (exp(-0.5 * aics$delta) / sum(exp(-0.5 * aics$delta)))

    # write the results for this tree
    full <- data.frame(treeid=treenum, His.full$loglik, aics["Full.HiSSE","delta"], aics["Full.HiSSE","W"])
    bisse <- data.frame(treeid=treenum, bisse.fit$loglik, aics["BiSSE","delta"], aics["BiSSE","W"])
    hidden.cave <- data.frame(treeid=treenum, His.cave$loglik, aics["Hidden.Cave","delta"], aics["Hidden.Cave","W"])
    hidden.surface <- data.frame(treeid=treenum, His.surface$loglik, aics["Hidden.Surface","delta"], aics["Hidden.Surface","W"])
    cid4 <- data.frame(treeid=treenum, cid4.hisse$loglik, aics["CID4","delta"], aics["CID4","W"])
    cid2 <- data.frame(treeid=treenum, cid2.hisse$loglik, aics["CID2","delta"], aics["CID2","W"])
    geosse <- data.frame(treeid=treenum, geosse.fit$lnLik, aics["GeoSSE","delta"], aics["GeoSSE","W"])

    
    write.table(full, file=outfile1, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
    write.table(bisse, file=outfile2, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
    write.table(hidden.cave, file=outfile3, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
    write.table(hidden.surface, file=outfile4, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
    write.table(cid4, file=outfile5, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
    write.table(cid2, file=outfile6, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
    write.table(geosse, file=outfile7, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)


    message(paste("Done fitting tree", treenum))
}

#----------
# Run the analysis
#----------

treesubset <- 1:1000

# 
junk <- mclapply(treesubset, run.MLE, mc.cores=16)