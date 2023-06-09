### Functions used in analyses to assist with iterating analyses over datasets. 
# Also includes functions published by Adams 2014 Syst. Biol. (phylogenetic transformation) 
# and Moen et al. 2013 (two-block partial least-squares analysis). 

###	From Evolution submission: "Testing adaptive radiation: a new test with Madagascar frogs"

### Phylogenetic transformation of data, following Adams 2014 (Syst. Biol.) and Garland and Ives (2000 AmNat)
	# INPUT:
		# (1) 'tree' is a phylogeny, assumed to have exactly the same taxa as the data
		# (2) 'X' are the continuous data, which are assumed to be a matrix (i.e. more than one character)
	# OUTPUT:
		# (1) 'U' is the data matrix rotated to show no Brownian covariance fitting the tree
	phylo.transform <- function(tree, X){
		# Phylogenetic VCV
			C <- vcv.phylo(tree)
		# Formatting
			X <- as.matrix(X)
		# Need to make a matrix of the phylogenetic means
			invC <- solve(C)
			aX <- colSums(invC %*% X)/sum(invC) # Ancestral values, from 'phytools' and based 
				# on Rohlf 2001 Evolution (p. 2147). Just reversing order of matrix computation, since numerator
				# is t(ones) %*% invC %*% X, and denominator is t(ones) %*% invC %*% ones
			aX <- rep(1,nrow(X)) %*% t(aX) # Turn into n x p matrix for subtraction from data
			difs <- X[tree$tip.label,] - aX
		# SVD of C to get D (from Blankers et al. 2012 J. Evol. Biol. Appendix 1)
			sing <- svd(C) # Break down phylogenetic VCV
			D <- solve(sing$u %*% diag(sqrt(sing$d)) %*% t(sing$u)) # Determine transformation matrix
			# head(round(D %*% phy.mat %*% t(D), 5)) # Dummy check - should be n x n identity matrix
		# Transformation
			U <- D %*% difs # Adams 2014 Eqn 4, reversed to be consistent with Garland and Ives 2000 (and to make it work!)
			U <- as.data.frame(U)
			rownames(U) <- rownames(X[tree$tip.label,]) # Label following taxon order of analysis (i.e. tree)
			U <- U[rownames(X),] # Reorder following original data matrix
			return(U)
	}
	

###  Two-block partial least-squares analysis following Rohlf and Corti 2000 Syst. Biol. 49:740-753
	# INPUT:
		# (1) "dat" is a dataframe with only continuous data (i.e., that which is to be analyzed and nothing more)
		# (2) "nvar1" is the number of variables in the first set (function assumes that all of "data" will be analysed)
		# (3) "correlation" is a boolean option that asks whether you want to do the analysis on the correlation (T) or covariance (F) matrix
		# (4) "give.scores" is a boolean option for whether you want to output scores on the 2B-PLS axes
	# OUTPUT:
		# (1) a correlation (or covariance, depending on the option you choose) matrix of the original variables
		# (2) the maximum possible covariance between all variables
		# (3) the proportion of the maximum possible covariance that is described by the analysis (this is often quite low, 
		# 	even with strong relationships among variables)
		# (4) the eigenvalues (as "d")
		# (5) the eigenvectors that describe how the original variables load onto the 2B-PLS variables 
		#	("u" for variable set 1, and "v" for variable set 2)
		# (6) the correlation between the 2B-PLS variables in set one and two, give for each axis
	tbpls <- function(dat, nvar1, correlation = T, give.scores = T){
		if (correlation){
			cov.data <- cor(dat)
		} else {
			cov.data <- cov(dat)
		}
		max.covar <- sum(sqrt(diag(cov.data[1:nvar1,1:nvar1]))) * sum(sqrt(diag(cov.data[(nvar1+1):ncol(dat),(nvar1+1):ncol(dat)]))) # Maximum possible covariance among sets of variables
		tbpls.data <- svd(cov.data[1:nvar1,(nvar1+1):ncol(dat)]) 
		prop.covar <- sum(tbpls.data$d)/max.covar #gives amount of potential covariance that is described by 2B-PLS
		# tbpls.data$u <- tbpls.data$u*(-1) #for some reason R turns out negative values (or at least opposite to Rohlf and Corti)
		# tbpls.data$v <- tbpls.data$v*(-1) #for some reason R turns out negative values (or at least opposite to Rohlf and Corti)
		rownames(tbpls.data$u) <- colnames(dat)[1:nvar1] 
		rownames(tbpls.data$v) <- colnames(dat)[(nvar1+1):ncol(dat)] 
		ndim <- min(c(nvar1,(ncol(dat) - nvar1))) #number of 2B-PLS dimensions
		scores <- matrix(ncol=(2*ndim),nrow=nrow(dat)) 
		rownames(scores) <- rownames(dat) 
		colnames(scores) <- c(paste0('set1.',seq(ndim)), paste0('set2.',seq(ndim)))
		for (j in seq(ndim)){
			for (k in 1:nrow(dat)){
				scores[k,j]<-sum(t(tbpls.data$u[,j]) * dat[k,1:nvar1]) # New latent variable from set 1
				scores[k,(j+ndim)]<-sum(t(tbpls.data$v[,j]) * dat[k,(nvar1+1):ncol(dat)]) # New latent variable from set 2
			}
		}
		correlation <- numeric(ndim) # Calculate the correlation for each axis of the 2B-PLS latent variables
		for (j in seq(ndim)){
			correlation[j] <- cor(scores[,j], scores[,(j+ndim)]) 
		}
		if(give.scores){
			out <- list(cov.matrix=cov.data,maximum.covariance=max.covar,prop.covar.explained=prop.covar,two.block=tbpls.data,axes.correlation=correlation,scores=scores)
		} else {
			out <- list(cov.matrix=cov.data,maximum.covariance=max.covar,prop.covar.explained=prop.covar,two.block=tbpls.data,axes.correlation=correlation)
		}
		return(out) 
	}


### Test of significance of individual 2B-PLS axes, as in Rohlf and Corti 2000
	# INPUT:
	# 	(1) Most input is the same as the above function
	#	(2) "nperms" is the number of permutations you choose to do for calculating the P-value
	#	(3) "sumtest" is a boolean that allows one to choose one of two tests: the default
	#		("sumtest=FALSE") compares the axes independently - in other words, it does not 
	#		consider how many variation has already been explained in lower PC dimensions.
	#		This seems reasonable but can give strange results, as in PC2 not being significant
	#		even when higher axes are significant. The "sumtest=TRUE" option sums lower dimensions
	#		and asks whether each PC dimension and all lower dimensions together explain 
	#		more variation than random.  This can be overly liberal because PC1 often explains
	#		so much variation that when lumped with higher dimensions, everything becomes significant.
	# OUTPUT:
	#	(1) Singular values (on which the test of significance is based); this should be the same as obtained above
	#	(2) P-values, given for each axis
	choose.axes.2BPLS <- function(data, nvar1, correlation = T, nperms = 1000, sumtest = FALSE){
		nvar2 <- ncol(data) - nvar1
		if (correlation){
			data <- scale(data)
			}
		ntax <- nrow(data) 
		data.1 <- as.data.frame(data[,1:nvar1]) 
		data.2 <- as.data.frame(data[,(nvar1+1):(nvar1+nvar2)]) 
		sing.values <- matrix(ncol = min(c(nvar1,nvar2)),nrow = nperms) 
		sing.values[1,] <- svd(cov(data)[1:nvar1, (1+nvar1):(nvar1+nvar2)])$d
		for (j in seq(nperms - 1)){
			rand.list <- sample(1:ntax, ntax) 
			permute.data <- as.matrix(cbind(data.1, data.2[rand.list,])) 
			big.R <- cov(permute.data) 
			R12 <- big.R[1:nvar1, (1+nvar1):(nvar1+nvar2)] 
			sing.values[j+1,] <- svd(R12)$d 
		}
		P <- vector(len = ncol(sing.values))
		if (!sumtest){ 
			for (i in seq(ncol(sing.values))){
				P[i] <- (length(sing.values[sing.values[,i] > sing.values[1,i],i]) + 1) / nperms 
			}
		} else { 
			mod.sing <- matrix(ncol = min(c(nvar1, nvar2)), nrow = nperms) 
			### First get sums of eigenvalues
			for (i in seq(nperms)){ 
				for (j in seq(ncol(sing.values))){
					mod.sing[i,j] <- sum(sing.values[i, seq(j)])
				}
			}
			### Now get P-values
			for (i in seq(ncol(mod.sing))){
				P[i] <- (length(mod.sing[mod.sing[,i] > mod.sing[1,i],i]) + 1) / nperms 
			}
		}	
		list(sample.sing.values = sing.values[1,], P.values = P) 
	}
	
### Simple function to loop over many diversification rates when conducting phylogenetic ANOVAs
	# INPUT:
	#	(1) 'dat' is a data frame with a grouping factor 'gp' for the ANOVA, and one or more diversification-rate variables
	#	(2) 'phy' is a phylogeny that matches 'dat' by rownames
	#	(3) 'var.name' is a single variable name that will in analyzed in an ANOVA test (i.e. one of these diversification-rate variables in 'dat')
	#	(4) 'sims' is the number of simulation replicates you want to use in the phylANOVA; default is same as 'phylANOVA()'
	# OUTPUT:
	#	(1) A single P-value from the phylogenetic ANOVA
	loop.ANOVA.pval <- function(dat, phy, var.name, sims = 1000){
		subdat <- dat[dat[, var.name] > 0, c('gp', var.name)] # Remove any rates == 0 or even negative (it happens with crown estimates of single-species families)
		subdat <- subdat[complete.cases(subdat),] # Remove families missing data
		phy <- keep.tip(phy, rownames(subdat))
		x <- subdat[phy$tip.label, 'gp'] 
		y <- subdat[phy$tip.label, var.name]
		names(x) <- names(y) <- phy$tip.label # Just to be absolutely sure there is no indexing error in 'phylANOVA'
		out <- phylANOVA(phy, x, y, nsim = sims)$'Pf' 
		return(out)
	}

### Function for rate simulations: subsampling a large tree
	# INPUT:
	#	(1) 'phy' is the phylogeny you want to use
	#	(2) 'node' is where you want to stretch the tree (Mantellidae here)
	#	(3) 'subsamp' is a two-element vector of the numbers of taxa you want to subsample from the focal clade (that indicated by 'node') and the outgroup, respectively
	#	(3) 'p' is the number of traits you want to simulate
	#	(4) 'R' is the trait correlation you want in your simulation (same correlation assigned to all pairwise combinations of traits)
	# 	(5) 'rates' is a vector of the rate ratios you want to consider
	#	(6) 'iter' is the number of simulation reps you want 
	#	(7) 'para.proc' is whether you want to do parallel processing with package 'parallel'
	#	(8) 'cores' is the number of parallel cores you want to use when analyzing your simulation replicates
	# OUTPUT:
	#	(1) 'summary' is a table of means and standard deviations for parameters across simulation replicates. It 
	#		also includes Type I error rate (for rate ratio of 1.0) or power for each rate tested (all rate ratios != 1.0)
	#	(2) 'results_objs' is a list of matrices, one for each rate ratio, that have parameter estimates from model fits for each simulation replicate 
	sigmamult_subsample_sims <- function(phy, node = 3466, subsamp, p, R, rates, iter = 100, para.proc = FALSE, cores = 2, pv = power.verbose){
		# Packages
			library(geomorph) # Version 3.3.1
			library(geiger) # Version 2.0.7
			library(motmot) # Version 2.1.3
			library(parallel) # Version 4.0.2
		# Taxa to sample
			focal <- extract.clade(phy, node)$tip.label
			other <- drop.tip(phy, focal)$tip.label
		# Initial setup of output matrix
			results <- matrix(ncol = 9, nrow = length(rates)) # Output matrix for results
			params <- c('focal','other','common','ratio')
			rownames(results) <- paste0('rr_',rates)
			colnames(results) <- c(paste0(params, '_mean'), paste0(params, '_sd'), 'power')
		# Run analysis
			all_fits <- list()
			for (j in seq(rates)){
				stretch.phy <- transformPhylo(phy, model = 'clade', rateType = 'clade', nodeIDs = node, cladeRates = rates[j])
				r <- matrix(R, ncol = p, nrow = p)
				diag(r) <- 1 # Force all rates to be 1.0, since simulating on stretched tree will
					# produce rate difference on stretched branches (when rates are analyze on the unstretched tree)
				x <- sim.char(stretch.phy,r,iter)
			# Subsample taxa for analysis
				# Choose taxa
					focal_set <- replicate(iter, sample(focal, subsamp[1])) # Each column is a rep of 'iter'
					other_set <- replicate(iter, sample(other, subsamp[2])) # Each column is a rep of 'iter'
				# Make vectors for grouping ('focal' and 'other'), then get data subset and tree for these taxa
					tmp <- list()
					for (i in seq(iter)){
						states <- c(rep('focal', subsamp[1]), rep('other', subsamp[2]))
						names(states) <- c(focal_set[,i], other_set[,i])
						dat <- x[names(states),,i]
						tree <- keep.tip(phy, names(states))
						tmp[[i]] <- list(states = states, dat = dat, tree = tree)
					}
			# Run
				if(para.proc){
					sx <- mclapply(tmp, function(y) compare.evol.rates(y$dat, y$tree, y$states, method='simulation', print.progress=F), mc.cores = cores)
				} else {
					sx <- lapply(tmp, function(y) compare.evol.rates(y$dat, y$tree, y$states, method='simulation', print.progress=F))
				}
			# Extract results from output
				out <- matrix(NA, ncol = 5, nrow = iter)
				colnames(out) <- c('focal.rate','other.rate','common.rate','rate.ratio','P-val')
				for (i in seq(iter)) out[i,] <- with(sx[[i]], c(sigma.d.gp, sigma.d.all, sigma.d.ratio, P.value))
				all_fits[[j]] <- out
			# Get mean ratios, power
				results[j,1:4] <- colMeans(out[,1:4])
				results[j,5:8] <- apply(out[,1:4], 2, sd)
				results[j,'power'] <- mean(out[,'P-val'] < 0.05) 
				if(pv) print(paste('Done with rate.ratio =',rates[j],'at',date()))
			}
		# Kick out results
			results <- as.data.frame(results)
			names(all_fits) <- paste0('rate_ratio_',rates)
			return(list(summary = results, results_objs = all_fits))
	}