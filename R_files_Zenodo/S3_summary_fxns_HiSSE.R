### Functions for HiSSE analyses for Evolution submission, "Improving inference and avoiding 
#	over-interpretation of hidden-state diversification models: specialized plant breeding 
#	has no effect on diversification in frogs"

### Function to summarize models
	# INPUT:
	#	(1) 'out_mods' is a list of fitted HiSSE models, with names that correspond to the models
	#	(2) 'type' only affects parameter numbers: 'std' applies to most models (including BiSSE), 'CID4' applies that model only
	#	(3) 'params' is a vector of parameter names you would like to extract from 'out_mods'. All parameters MUST be in all models. 
	#		If you leave this NULL you will only get lnL, AICc, AICc weights, and number of free parameters
	#	(4) 'spec' specifies whether you want speciation and extinction rates or not
	#	(5) 'net.div' specifies whether you want net diversification rates or not
	# OUTPUT:
	#	(1) A matrix of log-likelihoods, AICc values, AICc weights, parameter numbers, and parameter estimates
	mod_summary <- function(out_mods, type = c('std','CID4'), params = NULL, spec = FALSE, net.div = FALSE){
		if(length(type) == 2) stop('You need to specify which type of model was estimated, between the options "std" and "CID4"')
		b <- !is.null(params) # Makes subsequent code much simpler
		if(b){
			tu.nm <- grep('turn', params, value = T) # Names of turnover parameters
			ex.nm <- grep('eps', params, value = T) # Names of extinction fraction parameters	
			st.nm <- gsub('eps', '', ex.nm) # State combinations, for later naming speciation, extinction, and net diversification rates
		}
		npars <- 5 # Baseline output: model name, lnL, AICc, weight, k
		cntmp <- c('model', 'lnL', 'AICc', 'AICcwgt','npars') # First set of column names
		if(b){
			npars <- npars + length(params)
			cntmp <- c(cntmp, params)
		}	
		if(spec){
			sp.nm <- paste0('lam',st.nm)
			mu.nm <- paste0('mu',st.nm)
			if(type == 'std'){
				npars <- npars + 8 # Adding 4 speciation rates and 4 extinction rates
			} else {
				npars <- npars + 16 # Adding 8 speciation rates and 8 extinction rates
			}
			cntmp <- c(cntmp, sp.nm, mu.nm)
		}
		if(net.div){
			nd.nm <- paste0('r',st.nm)
			if(type == 'std'){
				npars <- npars + 4 # Four net diversification rates for 0A, 1A, 0B, 1B
			} else {
				npars <- npars + 8 # Eight rates for 0A, 1A, 0B, etc.
			}
			cntmp <- c(cntmp, nd.nm)
		}
		if(class(out_mods) == "hisse.fit") out_mods <- list(out_mods) # Allows us to avoid duplicating code below ('lapply()' functions)
		x <- as.data.frame(matrix(ncol = npars, nrow = length(out_mods)))
		if(length(out_mods) == 1){
			x <- x[,-1] # Don't need 'model' column
			colnames(x) <- cntmp[-1] # Remove 'model' term
		} else {	
			colnames(x) <- cntmp
			x[['model']] <- names(out_mods)
		}
		x[['lnL']] <- unlist(lapply(out_mods, function(y) y$loglik))
		x[['AICc']] <- unlist(lapply(out_mods, function(y) y$AICc))
		x[['AICcwgt']] <- aic.weight(x[['AICc']])
		if(type == 'CID4') x[['npars']] <- unlist(lapply(out_mods, function(y) max(y$index.par)))
		if(type == 'std') x[['npars']] <- unlist(lapply(out_mods, function(y) length(y$starting.vals))) # Each parameter has a starting value
		# Extract parameter values, assuming your output argument in HiSSE was 'turnover'
		if(b) x[params] <- t(sapply(out_mods, function(y) y$solution[params]))	
		# Give matrix for turnover and extinction rates (used in next two steps)
			if(b & spec | net.div){
				tu <- x[tu.nm] # Turnover rates
				ex <- x[ex.nm] # Extinction fractions
				lam <- tu / (ex + 1) # Speciation rates
				mu <- (ex * tu) / (ex + 1) # Extinction rates
			}
		# Assign speciation and extinction rates
			if(b & spec){
				x[sp.nm] <- lam
				x[mu.nm] <- mu
			}
		# Calculate net diversification rates
			if(b & net.div) x[nd.nm] <- lam - mu
		# Kick it out
			return(x)
	}

### Function to calculate AIC weights (following Burnham and Anderson 2002)
# Input is a vector of AIC values for various models (from the same dataset)
# Output is a vector of weights in the same order, plus AIC differences and ranks of models if you prefer
aic.weight <- function(x, AICdif = FALSE, rank.mods = FALSE, rd.digits = 3){
	l <- min(x) # Calculate the lowest AIC score
	dif <- x - l # Differences
	d <- exp(-0.5*(dif)) # Exponentiate
	if(is.null(rd.digits)){
		wg <- d / sum(d)
	} else {
		wg <- round(d / sum(d), digits = rd.digits)
		dif <- round(dif, digits = rd.digits)		
	}
	ranks <- rank(dif)
	if (AICdif & rank.mods) return(as.data.frame(list(AICdif = dif, AICweight = wg, AICrank = ranks)))
	else if(AICdif) return(as.data.frame(list(AICdif = dif, AICweight = wg)))
	else if(rank.mods) return(as.data.frame(list(AICweight = wg, AICrank = ranks)))
	else return(wg)
}


### Summary function for 'SupportRegionfHiSSE()'
#	Calculates speciation, extinction, and net diversification rates. It also summarizes across a list of such analyses
#	if looped under the same model (e.g. to allow parallelization)
	# INPUT:
	#	(1) 'ob' is either a single fitted object from 'SupportRegionfHiSSE()' or a list of such objects conducted under the same model
	#	(2) 'loop_list' is the default FALSE if it is a single output from 'SupportRegionfHiSSE' or TRUE if it is a list of such outputs
	#	(3) 'spec' specifies whether you want speciation and extinction rates or not
	#	(4) 'net.div' specifies whether you want net diversification rates or not
	# OUTPUT:
	#	(1) 'support_region' is a matrix of the support region with bounds and MLE
	#	(2) 'all_values' is a matrix of all diversification parameters pulled from the output(s)
	sum_div_supreg <- function(ob, loop_list = FALSE, net.div = T, spec = T){
		# Choose only diversification variables
			library(dplyr)
			library(tibble)
			if(loop_list){
				x <- as_tibble(ob[[1]]$points.within.region) %>%
					select(lnLik, contains(c('0A','1A','0B','1B'))) %>%
					select(lnLik, starts_with('turnover'), starts_with('eps'))			
				for(i in seq(ob)[-1]){
					tmp <- as_tibble(ob[[i]]$points.within.region) %>%
						select(lnLik, contains(c('0A','1A','0B','1B'))) %>%
						select(lnLik, starts_with('turnover'), starts_with('eps'))
					x <- as_tibble(rbind(x, tmp[-1,])) # First row is always MLE, which is repeated across searches
				}
				x <- arrange(x, lnLik) # Put in order of -lnLik
			} else {
				x <- as_tibble(ob$points.within.region) %>%
					select(lnLik, contains(c('0A','1A','0B','1B'))) %>%
					select(lnLik, starts_with('turnover'), starts_with('eps'))			
			}
		# Give matrix for turnover and extinction rates (used in next two steps)
			if(spec | net.div){
				tu.nm <- str_subset(names(x), 'turn') # Names of turnover parameters
				ex.nm <- str_subset(names(x), 'eps') # Names of extinction fraction parameters	
				st.nm <- gsub('eps', '', ex.nm) # State combinations, for later naming speciation, extinction, and net diversification rates
				tu <- x[tu.nm] # Turnover rates
				ex <- x[ex.nm] # Extinction fractions
				lam <- tu / (ex + 1) # Speciation rates
				mu <- (ex * tu) / (ex + 1) # Extinction rates
				colnames(mu) <- paste0('mu',st.nm)
				colnames(lam) <- paste0('lam',st.nm)
			}
		# Calculate net diversification rates
			if(net.div){
				nd <- lam - mu
				colnames(nd) <- paste0('net.div',st.nm)
				x <- as_tibble(cbind(x, nd))
			}
		# Assign speciation and extinction rates
			if(spec) x <- as_tibble(cbind(x, lam, mu))
		# Calculate bounds for all parameters
			sr <- matrix(nrow = 3, ncol = ncol(x))
			sr[1,] <- apply(x, 2, min) # Minimum value within set
			sr[2,] <- unlist(x[1,]) # Values at MLE
			sr[3,] <- apply(x, 2, max) # Maximum value within set
			colnames(sr) <- names(x)				
			sr <- as_tibble(sr) %>%
				add_column(what = c('lower','MLE','upper')) %>%
				mutate(lnLik = -1 * lnLik) %>%
				select(what, everything())
		# Finish
			x <- mutate(x, lnLik = -1 * lnLik) # At end, since it simplifies confidence-interval calculations 
			out <- list(support_region = sr, all_values = x)
			return(out)
	}	