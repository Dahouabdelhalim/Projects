# Computes the G matrices from the mcmc chains
# The script does not run the model on the diallel data, 
# this is not the topic from the current paper. 

datadir <- "./data"

# the way to collect the data is not entirely clean. The chain is 
# stored in one or several Rdata files for each population.
# Chains are in a "mod" variable in independent R sessions. 

get.VCV <- function(pop, excluding.self=FALSE) {
	library(coda) # for as.mcmc

	ans <- NULL
	MCMC.files <- list.files(path=datadir, pattern=paste0(pop, "_MCMC", if (excluding.self) "_excludingself" else "", "\\\\.RData$"), full=TRUE)
	for (chainfile in MCMC.files) {
		load(chainfile)
		ans <- rbind(ans, mod$VCV)
		rm("mod")
	}
	return(as.mcmc(ans))
}


get.matrices <- function(pop, scale=c("data", "mcmc")[1], excluding.self=FALSE) {
	VCV <- get.VCV(pop, excluding.self=excluding.self)
	
	genet.VCV <- grep(colnames(VCV), pattern="animal")
	resid.VCV <- grep(colnames(VCV), pattern="units")
	indiv.VCV <- grep(colnames(VCV), pattern="ind")
	
	G.VCV <- VCV[,genet.VCV]
	E.VCV <- VCV[,resid.VCV]+VCV[,indiv.VCV]
	P.VCV <- VCV[,genet.VCV]+VCV[,resid.VCV]+VCV[,indiv.VCV]
	
	# In the MCMC procedure, phenotypes are multiplied by 100, and variances by 100*100
	scale.factor <- 1
	if (scale == "data") scale.factor <- 1/10000
	
	ans <- list(
		G.VCV = G.VCV*scale.factor,
		E.VCV = E.VCV*scale.factor,
		P.VCV = P.VCV*scale.factor,
		G = matrix(apply(G.VCV*scale.factor, 2, median), ncol=2),
		E = matrix(apply(E.VCV*scale.factor, 2, median), ncol=2),
		P = matrix(apply(P.VCV*scale.factor, 2, median), ncol=2),
		Gv = matrix(apply(G.VCV*scale.factor, 2, var), ncol=2),
		Ev = matrix(apply(E.VCV*scale.factor, 2, var), ncol=2))
}
