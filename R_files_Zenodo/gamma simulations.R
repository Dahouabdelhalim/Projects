### under Yule model, parameters estimated from original tree

norig <- 115  ## number of tips in original tree (115 in overall, 67 in pruned)
pm <- 20  ## percent missing taxa
nterm <- ceiling(norig/(1-pm/100)) ## total terminals including deleted taxa
nreps <- 10000
b <- 0.534 ## birth rate; 0.534 for Yule model for original tree, 0.295 for bd model for original tree
d <- 0 ## death rate; equals zero under Yule model, 0.116 for bd model for original tree (=extinction fraction (a) * speciation (lambda))

glist <- vector("numeric", nreps)

for (i in 1:nreps) {
	to_drop <- sample(1:nterm,size=nterm-norig)
	sample_tree <- drop.tip(rphylo(nterm,b,d),to_drop)
	glist[i] <- gammaStat(sample_tree)
	}
glist <- sort(round(glist,4))
n_cutoff <- ceiling(nreps*0.05+1)
cat("For",nreps,"simulated trees of size",nterm,"with",nterm-norig,"missing tips (",pm,"%):\\n\\tgamma (mean)= ",mean(glist),"\\n\\tgamma(95%)= ",glist[n_cutoff],"\\n\\t( b=",b,"d=",d,")\\n")
