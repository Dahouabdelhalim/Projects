### under Yule model, parameters estimated from original tree

norig <- 115  ## number of tips in original tree (115 in overall, 67 in pruned)
pm <- 20  ## percent missing taxa
nterm <- ceiling(norig/(1-pm/100)) ## total terminals including deleted taxa
nreps <- 1000
b <- 0.534 ## birth rate; 0.534 for Yule model for original tree, 0.295 for bd model for original tree
d <- 0 ## death rate; equals zero under Yule model, 0.116 for bd model for original tree (=extinction fraction (a) * speciation (lambda))

dlist <- vector("numeric", nreps)  ## list with deltaAIC values (min rate constant model-min DD model)

for (i in 1:nreps) {
	to_drop <- sample(1:nterm,size=nterm-norig)
	sample_tree <- drop.tip(rphylo(nterm,b,d),to_drop)
	fit.res <- fitdAICrc(branching.times(sample_tree),model=c("pureBirth","bd","DDL","DDX"))$dAIC
	dlist[i] <- min(fit.res[1:2]) - min(fit.res[3:4]) ## difference in min rate constant model minus min DD model
	}
dlist <- sort(round(dlist,4),decreasing=TRUE)
n_cutoff <- ceiling(nreps*0.05+1)
cat("For",nreps,"simulated trees of size",nterm,"with",nterm-norig,"missing tips (",pm,"%):\\n\\tdAIC (mean)= ",mean(dlist),"\\n\\tdAIC(95%)= ",dlist[n_cutoff],"\\n\\t( b=",b,"d=",d,")\\n\\t(positive values of dAIC > the 95% value favor density denpent models)")
