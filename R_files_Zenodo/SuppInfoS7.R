### Function for testing adaptive radiation, via phenotypic rate comparisons
###	From Evolution submission: "Testing adaptive radiation: a new test with Madagascar frogs"

# Note: requires prior installation of packages 'geomorph', 'motmot', and 'phytools'. 
#	For parallel processing of power simulations, uses package 'parallel'
	# INPUT:
	#	(a) 'tree': Phylogeny used for analyses, as class 'phylo'
	#	(b) 'ar.taxa': Vector of taxa considered the adaptive radiation
	#	(c) 'dat' : A data frame of continuous phenotypic data, with row names that match the phylogeny's tip labels
	#	(d) 'rr' : Vector of rate ratios to test for power
	#	(e) 'sims' : the number of simulation replicates you want to use in power simulations
	#	(f) 'power.verbose' : whether you want the function to return a statement after it evaluates each rate ratio in power simulations
	#	(g) 'parallel.proc': If using parallel processing with package 'parallel', this is the number of parallel nodes you want to use 
	# OUTPUT:
	#	(a) 'results.summary' : a vector of predicted group differences in rate, observed differences, common rate for the whole phylogeny, and P-value of Adams 2014 test
	#	(b) 'power.curve' : a matrix of hypothetical rate ratios and their statistical power on your phylogeny, under the designated candidate adaptive radiation
	#	(c) 'fitted.sigma.mult' : a fitted object from the geomorph function 'compare.evol.rates'
candidate.adap.rad <- function(tree, ar.taxa, dat, rr = seq(1.5, 4, by = 0.5), sims = 1000, power.verbose = T, parallel.proc = NULL){
	# Load any necessary libraries
		library(geomorph) # Version 3.3.1
		library(phytools) # Version 0.7-70
		library(motmot) # Version 2.1.3
		if(!is.null(parallel.proc)) library(parallel) # Version 4.0.2
	# Designate auxiliary functions
		### Calculating mean trait correlation (from 'evolvcv.lite' in phytools)
			# Input
			#	(1) 'tr' : phylogeny
			#	(2) 'd' : data
			# Output: mean trait correlation, calculated with Brownian motion on phylogeny (Revell and Collar 2009 Evolution)
			cor.ml <- function(tr, d){
				# Prep calculations
					X <- as.matrix(d)
					rownames(X) <- rownames(d)
					colnames(X) <- colnames(X)
					X <- X[tr$tip.label, ]
					n <- nrow(X)
					m <- ncol(X)
					ones <- diag(diag(n))
					D <- matrix(0, n * m, m)
					for (i in 1:(n * m)) for (j in 1:m) if ((j - 1) * n < i && i <= j * n) D[i, j] = 1
					y <- as.matrix(as.vector(X))
					C <- vcv.phylo(tr)
				# Get PIC means as estimates of ancestral-state values
					sv <- vector(length = m)
					for (i in seq(p)){
						sv[i] <- mean(pic(X[, i], tr)^2)
					}
				# Calculate correlations
					sig <- (t(X - ones %*% t(sv)) %*% solve(C) %*% (X - ones %*% t(sv)))/n
					cor.sig <- matrix(ncol=m, nrow = m)
					for (i in seq(m)){
						for (j in seq(m)){
							cor.sig[i,j] <- sig[i,j] / sqrt(sig[i,i]*sig[j,j])
						}
					}
					out <- mean(cor.sig[lower.tri(cor.sig)])
					return(out) # mean correlation across traits
			}			
		### Function for power simulations
			# INPUT:
			#	(1) 'phy' is the phylogeny you want to use
			#	(2) 'ar.taxa' is a taxon list of your candidate adaptive radiation
			#	(3) 'p' is the number of traits you want to simulate
			#	(4) 'R' is the average trait correlation you want in your simulation
			# 	(5) 'rates' is a vector of the rate ratios you want to consider
			#	(6) 'iter' is the number of simulation reps you want 
			#	(7) 'para.proc' is whether you want to do parallel processing with package 'parallel'
			#	(8) 'nodes' is the number of parallel nodes you want to use when analyzing your simulation replicates
			sigma.mult.powersims <- function(phy, ar.taxa, p, R, rates, iter=1000, para.proc = FALSE, nodes=2, pv = power.verbose){
				# Initial setup of output matrix
					results <- matrix(ncol = 5, nrow = length(rates)) # Output matrix for results
					rownames(results) <- paste0('rr.',rates)
					colnames(results) <- c('focal.rate','other.rate','common.rate','rate.ratio','power')
				# Tree stuff (adaptive radiation on this tree is represented by 'ar.taxa')
					ar.node <- findMRCA(phy, ar.taxa)
					states <- rep('other',length(phy$tip.label))
					names(states) <- phy$tip.label
					states[ar.taxa] <- 'AR' 
				# Run analysis
					for (j in seq(rates)){
						stretch.phy <- transformPhylo(phy, model = 'clade', rateType = 'clade', nodeIDs = ar.node, cladeRates = rates[j])
						r <- matrix(R, ncol = p, nrow = p)
						diag(r) <- 1 # Force all rates to be 1.0, since simulating on stretched tree will
							# produce rate difference on stretched branches (when rates are analyze on the unstretched tree)
						x <- sim.char(stretch.phy,r,iter)
					# Prep list for cluster/lapply
						tmp <- list()
						for (k in seq(iter)) tmp[[k]] <- x[,,k]
					# Run
						if(para.proc){
							sx <- mclapply(tmp, function(y) compare.evol.rates(y,phy,states,method='simulation',print.progress=F), mc.cores = nodes)
						} else {
							sx <- lapply(tmp, function(y) compare.evol.rates(y,phy,states,method='simulation',print.progress=F))
						}
					# Extract results from output
						out <- matrix(NA, ncol = 5, nrow = iter)
						colnames(out) <- c('focal.rate','other.rate','common.rate','rate.ratio','P-val')
						for (k in seq(iter)) out[k,] <- with(sx[[k]], c(sigma.d.gp, sigma.d.all, sigma.d.ratio, P.value))
					# Get mean ratios, power
						results[j,1:4] <- colMeans(out[,1:4])
						results[j,'power'] <- mean(out[,'P-val'] < 0.05) 
						if(pv) print(paste('Done with rate.ratio =',rates[j],'at',date()))
					}
				# Kick out results
					return(as.data.frame(results))
			}
	# Number of traits
		p <- ncol(dat)		
	# Mean trait correlation
		R <- cor.ml(tree, dat) # See code below for function 'cor.ml'
	# Define vector for taxon assignment
		states <- rep('other',length(tree$tip.label))
		names(states) <- tree$tip.label
		states[ar.taxa] <- 'AR' 			
	# Calculate base rate of sigma.mult for whole tree
		sigma.mult.out <- compare.evol.rates(dat, tree, states, iter = 1000, print.progress = F)
		sigma.mult.all <- sigma.mult.out$sigma.d.all
	# Then power analyses:
		# Set up parallel-processing options
			pp <- ifelse(is.null(parallel.proc), FALSE, TRUE)
			nd <- parallel.proc
		# Do power analyses
			pow <- sigma.mult.powersims(phy = tree, ar.taxa = ar.taxa, p = p, R = R, rates = rr, iter=sims, para.proc=pp, nodes=nd)[,c('common.rate','power')]
		# Check whether your range of rate ratios was sufficient to get power of 0.80 (for predictions)
			if(sum(pow[,'power'] >= 0.8) == 0){
				print('Your rate ratios did not produce power sufficient power to make a prediction')
				power.predict <- FALSE
			} else power.predict <- TRUE 
	# Empirical analyses output
		if(power.predict){
			emp.results <- with(sigma.mult.out, c(sigma.d.gp, sigma.d.ratio, sigma.d.all, P.value))
			i <- which.max(pow[,'power'] >= 0.8)
			r.pred <- as.numeric(unlist(strsplit(rownames(pow)[i],'rr.'))[2])
			pred.rates <- sigma.mult.all * c(r.pred,1)/pow[i,'common.rate'] # Multiple empirical common rate by the proportionality of group rates vs. common rate in power simulations
			pred.rates <- c(pred.rates, r.pred) # Add rate ratio
			emp.results <- c(pred.rates,emp.results)						
		} else {
			emp.results <- with(sigma.mult.out, c(sigma.d.gp, sigma.d.ratio, sigma.d.all, P.value))
		}
		index <- ifelse(power.predict,1,4):8
		names(emp.results) <- c('AR.pred','out.pred','rr.pred','AR.obs','out.obs','rr.obs','common.rate','P-value')[index] # Only use names for the data
	# Output
		power.out <- as.data.frame(cbind(rr, pow[,'power']))
		colnames(power.out) <- c('rate.ratio','power')
		out <- list(results.summary = emp.results, power.curve = power.out, fitted.sigma.mult = sigma.mult.out)
		return(out)
}	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	