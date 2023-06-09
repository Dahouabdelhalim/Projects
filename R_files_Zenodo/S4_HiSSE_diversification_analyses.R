### HiSSE diversification analyses for Evolution submission, "Improving inference and avoiding 
#	over-interpretation of hidden-state diversification models: specialized plant breeding 
#	has no effect on diversification in frogs"


############################################
### Load relevant packages and functions ###
############################################

	library(parallel)
	library(hisse) # Used version 1.9.8; note that newer versions use somewhat different 
				   # function names and have different default likelihood search options
	source('SI3_summary_fxns_HiSSE.R') # Loads functions 'mod_summary' used below

#################
### Load data ###
#################

	tree <- read.tree('Tonini.tree.1579sp.tre')
	dat <- read.csv('Tonini.phyto.dat.1579sp.csv')
	frac <- c(0.4994893, 0.6666667) # Sampling fraction of the two states, based on TEA's full data matrix of 3,105 species with state information

#####################################
### Set up diversification models ###
#####################################

### (1) Transition matrices
	# Base matrices
		trans_rates <- TransMatMaker(hidden.states = TRUE)
		trans_rates_nodual <- ParDrop(trans_rates, c(3,5,8,10))
		trans_rates_bisse <- TransMatMaker(hidden.states = FALSE)
	# CID4, which is a bit more complicated
		trans_rates_CID4 <- TransMatMakerfHiSSE(hidden.traits = 3) # 3 hidden (binary) traits gives us 4 hidden states
		# Set all transitions from 0 -> 1 to be governed by a single rate:
			trans_8tr <- ParEqual(trans_rates_CID4, c(2,4,2,6,2,8))
		# Set all transitions from 1 -> 0 to be governed by a single rate:
			trans_8tr <- ParEqual(trans_8tr, c(1,3,1,4,1,5))
		# Set all transitions from A -> B and B -> A to be governed by a single rate:
		#	to_change <- cbind(c(1,3,2,4), c(3,1,4,2))
		#	trans_8tr[to_change] <- 3 # Already this way
		# Set all transitions from A -> C and C -> A to be governed by a single rate:
			to_change <- cbind(c(1,5,2,6), c(5,1,6,2))
			trans_8tr[to_change] <- 4
		# Set all transitions from A -> D and D -> A to be governed by a single rate:
			to_change <- cbind(c(1,7,2,8), c(7,1,8,2))
			trans_8tr[to_change] <- 5
		# Set all transitions from B -> C and C -> B to be governed by a single rate:
			to_change <- cbind(c(3,5,4,6), c(5,3,6,4))
			trans_8tr[to_change] <- 6
		# Set all transitions from B -> D and D -> B to be governed by a single rate:
			to_change <- cbind(c(3,7,4,8), c(7,3,8,4))
			trans_8tr[to_change] <- 7
		# Set all transitions from C -> D and D -> C to be governed by a single rate:
			to_change <- cbind(c(5,7,6,8), c(7,5,8,6))
			trans_8tr[to_change] <- 8
	# Put in list for easier organization
		tr <- list(
			bisse_dif = trans_rates_bisse,
			bisse_equal = ParEqual(trans_rates_bisse, c(1,2)),
			bisse_irrev = ParDrop(trans_rates_bisse, 1),
			nodual = trans_rates_nodual,
			trans_8tr = trans_8tr
		)
	# Clean up workspace
		rm(trans_rates, trans_rates_nodual, trans_rates_bisse, trans_rates_CID4, trans_8tr, to_change)

### (2) Diversification-rate vectors (turnover and extinction fraction have same variants, so just do once)
	dr <- list(
		Bnull = c(1,2,0,0), # Hidden state has no effect (same as BiSSE, unequal rates; note 
		#	that this will NOT work with 'hisse.new()', which I avoid because it throws errors 
		#	when using simulated annealing for BiSSE)
		Aeq_Bnull = c(1,1,0,0), # Equal rates, no hidden-rate variation (same as BiSSE, equal rates; 
		#	same implementation note as previous BiSSE model in lines 63â€“65)
		alldif = c(1,2,3,4), # All turnover rates different
		Adif_Beq = c(1,2,3,3), # Only the B values are the same
		cid2 = c(1,1,2,2), # Only hidden state affects turnover rates
		cid4 = rep(1:4, each = 2) # All hidden states different, focal states the same
	)

### (3) Set up models (easier to parallelize as part of a list)
	models <- list()
	models$full_bisse <- list(dr$Bnull, dr$Bnull, tr$bisse_dif)
	models$baseline_bisse <- list(dr$Aeq_Bnull, dr$Aeq_Bnull, tr$bisse_equal)
	models$deadend_bisse <- list(dr$Aeq_Bnull, dr$Bnull, tr$bisse_equal)
	models$suicide_bisse <- list(dr$Aeq_Bnull, dr$Bnull, tr$bisse_dif)
	models$lonely_bisse <- list(dr$Bnull, dr$Aeq_Bnull, tr$bisse_equal)
	models$irrev_bisse <- list(dr$Aeq_Bnull, dr$Aeq_Bnull, tr$bisse_irrev)
	models$mk2_bisse <- list(dr$Aeq_Bnull, dr$Aeq_Bnull, tr$bisse_dif)
	models$full2 <- list(dr$alldif, dr$alldif, tr$nodual)
	models$full3 <- list(dr$Adif_Beq, dr$Adif_Beq, tr$nodual)
	models$suicide4 <- list(dr$cid2, dr$alldif, tr$nodual)
	models$suicide6 <- list(dr$cid2, dr$Adif_Beq, tr$nodual)
	models$CID2 <- list(dr$cid2, dr$cid2, tr$nodual)
	models$CID4 <- list(dr$cid4, dr$cid4, tr$trans_8tr)

### (4) Name parts of models		
	for (i in seq(models)) names(models[[i]]) <- c('turnover','ext.frac','transitions')

### (5) Calculate number of parameters (to ensure I specified them correctly)
	param_num <- as.data.frame(matrix(ncol = 5, nrow = length(models)))
	param_num[,1] <- names(models)
	for(i in seq(models)){
		x <- models[[i]]
		param_num[i, 3:5] <- c(max(x$turnover), max(x$ext.frac), max(x$transitions, na.rm = T))
	}
	param_num[,2] <- rowSums(param_num[,3:5])
	colnames(param_num) <- c('model', 'total', 'turnover', 'ext.frac', 'transitions')

### (6) Clean up workspace
	rm(i, x, dr, tr, param_num)

#################################################################################################################################
### Estimate BiSSE models: default search options, bounding parameter space, simulated annealing, and the latter two together ###
#################################################################################################################################

### General idea: Parallelize searches over 7 BiSSE models using 'mclapply()', then
#	use a 'for' loop to run analyses across 4 likelihood search options

### (1) Specify models and number of 'for' loop runs
	test_mods <- models[grep('_bisse', names(models), value = T) ]
	cores <- length(test_mods) # 1 core / model; this will be an argument to 'mclapply()'; here so one can manually reduce if necessary
	hidden <- FALSE
	
### (2) Runs across different likelihood search options
	sc <- list(
			default = list(FALSE, NULL),
			parbounds = list(FALSE, c(10, 10, 100)),
			sann = list(TRUE, NULL),
			sann_pb = list(TRUE, c(10, 10, 100))
		)
	all_out <- list()	
	for (i in seq(sc)){
			y <- sc[[i]]		
		### Specify what kind of analysis you want to do
			mod_type <- ifelse(hidden, 'hisse', 'bisse')
			sim_ann <- y[[1]]
			par_bounds <- y[[2]] # Either 'NULL' or a vector of length 3 with bounds in this order: turnover, extinction fraction, and transition rate
			start_vals <- NULL # Either 'NULL' or a vector of length 3 with starting vals in this order: turnover, extinction fraction, and transition rate
			print(paste0('Running a ',
					mod_type, ' analysis ',
					'with sann = ',
					sim_ann,', ',
					ifelse(is.null(start_vals),"default starting params","user's starting params"),', and ',
					ifelse(is.null(par_bounds),"default parameter bounds","user's bounds for params")
				))
		### Implement search conditions
			# Using older version of 'hisse()', as simulated annealing doesn't seem to work for BiSSE with 'hisse.new()'
			if(is.null(par_bounds)){
				wrapper <- function(x){
					hisse(tree, as.data.frame(dat), f = frac, hidden.states = hidden, 
						turnover.anc = x$turnover, eps.anc = x$ext.frac, 
						trans.rate = x$transitions, sann = sim_ann, starting.vals = start_vals)
				}		
			} else {
				wrapper <- function(x){
					hisse(tree, as.data.frame(dat), f = frac, hidden.states = hidden, 
						turnover.anc = x$turnover, eps.anc = x$ext.frac, 
						trans.rate = x$transitions, sann = sim_ann, turnover.upper = par_bounds[[1]],
						eps.upper = par_bounds[[2]], trans.upper = par_bounds[[3]], starting.vals = start_vals)
				}
			}
		### Estimate models
			out_hisse <- mclapply(test_mods, wrapper, mc.cores = cores)
			names(out_hisse) <- names(test_mods)
			# Save intermediate results in case cluster times out
			save(out_hisse, file = paste0(
					mod_type,
					'_sann'[sim_ann],
					'_parbounds'[!is.null(par_bounds)],
					'_',cores,'cores',
					'_', Sys.Date(), '.Rdata')
			)	
		all_out[[i]] <- out_hisse # Add to total results
	}
	
### (3) Save concatenated results
	names(all_out) <- names(sc) # A two-layered list, with the first layer the 4 search type and the second the fits for the 7 BiSSE models within each search types
	save(all_out, file = paste0(
		mod_type,
		'_',length(test_mods),'mods',
		'_',length(sc),'searchtypes',
		'_',Sys.Date(),
		'.Rdata')
	)
	
### (4) Summarize results (Table 1)
	bisse_summary <- lapply(all_out, mod_summary, type = 'std')
	names(bisse_summary) <- names(sc)	
	table1 <- data.frame(	bisse_summary$default$model, 
							bisse_summary$default$npars, 
							TEA = c(-6719.070,-6738.642,-6733.113,-6728.748,-6733.571,-6772.158,-6681.715), # Manually adding values from Tonini et al. 2020
							bisse_summary$default$lnL,
							bisse_summary$parbounds$lnL,
							bisse_summary$sann$lnL,
							bisse_summary$sann$lnL - bisse_summary$default$lnL,
							bisse_summary$sann$AICc,
							bisse_summary$sann$AICcwgt
						)
	table1[,4:9] <- round(table1[,4:9], 3)
	colnames(table1) <- c('Model', 'k', 'TEA', 'Default search', 'Bounded parameters', 'Simulated annealing', 'Difference', 'AICc', 'wi')
	table1
	rm(bisse_summary, all_out, table1) # Clean workspace, as future objects will take name 'all_out'


###################################################################################################################
### Estimate top HiSSE models: default search, bounding parameters, sim. annealing, and the latter two together ###
###################################################################################################################

### General idea: Parallelize searches over 7 HiSSE models using 'mclapply()', then
#	use a 'for' loop to run analyses across 4 likelihood search options

### (1) Specify models and number of 'for' loop runs
	test_mods <- models[c('full2','full3','suicide4','suicide6','CID2','CID4')]
	cores <- length(test_mods) # 1 core / model; this will be an argument to 'mclapply()'; here so one can manually reduce if necessary
	hidden <- TRUE
	
### (2) Runs across different likelihood search options
	sc <- list(
			default = list(FALSE, NULL),
			parbounds = list(FALSE, c(10, 10, 100)),
			sann = list(TRUE, NULL),
			sann_pb = list(TRUE, c(10, 10, 100))
		)
	all_out <- list()	
	for (i in seq(sc)){
			y <- sc[[i]]		
		### Specify what kind of analysis you want to do
			mod_type <- ifelse(hidden, 'hisse', 'bisse')
			sim_ann <- y[[1]]
			par_bounds <- y[[2]] # Either 'NULL' or a vector of length 3 with bounds in this order: turnover, extinction fraction, and transition rate
			start_vals <- NULL # Either 'NULL' or a vector of length 3 with starting vals in this order: turnover, extinction fraction, and transition rate
			print(paste0('Running a ',
					mod_type, ' analysis ',
					'with sann = ',
					sim_ann,', ',
					ifelse(is.null(start_vals),"default starting params","user's starting params"),', and ',
					ifelse(is.null(par_bounds),"default parameter bounds","user's bounds for params")
				))
		### Implement search conditions
			# Using newer version of HiSSE for all other analyses
			if(is.null(par_bounds)){
				wrapper <- function(x){
					hisse.new(tree, as.data.frame(dat), f = frac, hidden.states = TRUE, 
						turnover = x$turnover, eps = x$ext.frac, trans.rate = x$transitions,
						sann = sim_ann, starting.vals = start_vals)
				}
			} else {
				wrapper <- function(x){
					hisse.new(tree, as.data.frame(dat), f = frac, hidden.states = TRUE, 
						turnover = x$turnover, eps = x$ext.frac, trans.rate = x$transitions,
						sann = sim_ann, turnover.upper = par_bounds[[1]], eps.upper = par_bounds[[2]],
						trans.upper = par_bounds[[3]], starting.vals = start_vals)
				}
			}
		### Estimate models
			out_hisse <- mclapply(test_mods, wrapper, mc.cores = cores)
			names(out_hisse) <- names(test_mods)
			# Save intermediate results in case cluster times out
			save(out_hisse, file = paste0(
					mod_type,
					'_sann'[sim_ann],
					'_parbounds'[!is.null(par_bounds)],
					'_',cores,'cores',
					'_', Sys.Date(), '.Rdata')
			)	
		all_out[[i]] <- out_hisse # Add to total results
	}
	
### (3) Save concatenated results
	names(all_out) <- names(sc) # A two-layered list, with the first layer the 4 search type and the second the fits for the 7 BiSSE models within each search types
	save(all_out, file = paste0(
		mod_type,
		'_',length(test_mods),'mods',
		'_',length(sc),'searchtypes',
		'_',Sys.Date(),
		'.Rdata')
	)
	
### (4) Summarize results (Table 2)
	hs <- lapply(all_out, mod_summary, type = 'std')
	table2 <- data.frame(	hs$default$model, 
							hs$default$npars, 
							hs$default$lnL,
							hs$parbounds$lnL,
							hs$sann$lnL,
							hs$sann_pb$lnL 
						)
	table2[,3:6] <- round(table2[,3:6], 3)
	colnames(table2) <- c('Model', 'k', 'Default search', 'Bounded parameters', 'Simulated annealing', 'Bounds and sim. annealing')
	table2
	rm(hs, all_out, table2) # Clean workspace, as future objects will take name 'all_out'

#################################################################
### Estimate top HiSSE models under variables starting values ###
#################################################################

### General idea: Use a 'for' loop to search across 6 different models, and within each loop
#	parallelize searches over different starting values with 'mclapply()'
#	Note: I only provide code for one of two analyses (bounded parameters); the searches 
#	with simulated annealing can be easily specified in section (1) below

### (1) Specify models and number of 'for' loop runs
	test_mods <- models[c('full2','full3','suicide4','suicide6','CID2','CID4')]
	nruns <- length(test_mods)
	cores <- 24 # User inputted based on computational resources
	hidden <- TRUE
	sim_ann <- FALSE # Change to 'TRUE' for second set of analyses
	par_bounds <- c(10,10,100)
	
### (2) Load starting parameter values, then put into a list, which is easier to parallelize
	sv <- read.csv('starting_vals_first24_exp_last24_unif.csv')
	start_vals <- list()
	for (j in seq(nrow(sv))) start_vals[[j]] <- sv[j,]
	rm(sv)
	
### (3) Runs across different models
	all_out <- list()	
	for (i in seq(test_mods)){
		### Specify model
			y <- test_mods[[i]]		
		### Specify what kind of analysis you want to do
			mod_type <- ifelse(hidden, 'hisse', 'bisse')
			print(paste0('Running a ',
					mod_type, ' analysis ',
					'with sann = ',
					sim_ann,', ',
					"user's starting params, and ",
					ifelse(is.null(par_bounds),"default parameter bounds","user's bounds for params")
				))
		### Implement search conditions (assumes variable starting parameters are input as a list)
			if(is.null(par_bounds)){
				wrapper <- function(x){
					hisse.new(tree, as.data.frame(dat), f = frac, hidden.states = TRUE, 
						turnover = y$turnover, eps = y$ext.frac, 
						trans.rate = y$transitions, sann = sim_ann, starting.vals = unlist(x))
				}
			} else {
				wrapper <- function(x){
					hisse.new(tree, as.data.frame(dat), f = frac, hidden.states = TRUE, 
						turnover = y$turnover, eps = y$ext.frac, 
						trans.rate = y$transitions, sann = sim_ann, turnover.upper = par_bounds[[1]],
						eps.upper = par_bounds[[2]], trans.upper = par_bounds[[3]], starting.vals = unlist(x))
				}
			}		
		### Estimate models
			out_hisse <- mclapply(start_vals, wrapper, mc.cores = cores)
			names(out_hisse) <- paste0('run',seq(start_vals))
			# Save intermediate results in case cluster times out
			save(out_hisse, file = paste0(
					mod_type,
					'_',names(test_mods)[i],
					'_sann'[sim_ann],
					'_parbounds'[!is.null(par_bounds)],
					'_var.start.vals',
					'_', Sys.Date(), '.Rdata')
			)	
		all_out[[i]] <- out_hisse # Add to total results
	}

### (4) Save all results
	names(all_out) <- names(test_mods) # A two-layered list, with the first layer the 7 models, and the deeper layer the 48 model fits from different starting values
	save(all_out, file = paste0(
		'hisse',
		'_',length(test_mods),'mods',
		'_sann'[sim_ann],
		'_parbounds'[!is.null(par_bounds)],
		'_var.start.vals',
		'_',Sys.Date(),
		'.Rdata')
	)
	
### (5) Summary done in a separate file (SI 13: figure code)

#########################################################################################################
### Calculate model-average parameter estimates for two models with nearly all AICc support (Table 4) ###
#########################################################################################################

### (1) Extract just the two best models for summary
	mods <- global_opt_models[c('full2','suicide4')]

### (2) Set parameters that I want to extract
	# Note that net diversification rates, speciation rates, and extinction rates are not specified here; they come out of 'mod_summary'
		states <- paste0(rep(c(0, 1), 2), rep(c('A', 'B'), each = 2))
		params_to_extract <- c(paste0('turnover', states), paste0('eps', states))
	
### (3) Extract all diversification parameters, then calculate model-averaged estimates
	out <- mod_summary(mods, type = 'std', params = params_to_extract, spec = T, net.div = T)
	avg <- colSums(out$AICcwgt * out[,6:ncol(out)])
	out[3,1] <- 'mod.avg'
	out[3,6:ncol(out)] <- avg

### (4) Make Table 4
	table4 <- data.frame(
						Rate = c('Turnover','Ext.frac','Speciation','Extinction','Net.div'), 
						matrix(round(avg, 3), ncol = 4, nrow = 5, byrow = TRUE)
						)
	colnames(table4)[2:5] <- c('NonA','PhytoA','NonB','PhytoB')
	table4 <- table4[c(3:5,1:2),] # Reorder to fit MS
	table4












