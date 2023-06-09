### R code for replicating figures for Evolution submission, 
#	"Improving inference and avoiding over-interpretation of hidden-state diversification 
#	models: specialized plant breeding has no effect on diversification in frogs"

############################################
### Load relevant packages and functions ###
############################################

	library(tidyverse)
	filter <- dplyr::filter # Force default to be 'dplyr' version
	library(hisse) # Used version 1.9.8; note that newer versions use somewhat different 
				   # function names and have different default likelihood search options
	source('SI3_summary_fxns_HiSSE.R') # Includes 'aic.weight' for summarizing model support

##################################################################################################################
### Figure 1: support for models across different likelihood search options, default starting parameter values ###
##################################################################################################################

### (1) Load data: a single 4-element list (for 4 types of searches), with each element having 6 model fits (the six models tested)
	load('SI7_fitted_HiSSEmods_defaultstartvals.Rdata') # Loads object 'all_out'
	
### (2) Calculate AICc weights for each model under each condition, including each model's MLE across the four searches
	hs <- lapply(all_out, mod_summary, type = 'std')
	aw <- data.frame(	hs$default$model,
						hs$default$AICc,
						hs$parbounds$AICc,
						hs$sann$AICc,
						hs$sann_pb$AICc
					)
	aw[,6] <- apply(aw[,2:5], 1, min)
	colnames(aw) <- c('model',paste0('AICc_',names(all_out)), 'AICc_min')
	aw <- as_tibble(aw) %>%
			select(model, starts_with('AICc_')) %>%
			mutate_if(.predicate = is.numeric, ~aic.weight(.x)) %>%
			rename_with(~gsub('AICc', 'AICcwgt', .x))

### (3) Make the bar graph
		names(aw) <- gsub('AICcwgt_', '', names(aw))
		aw <- pivot_longer(aw, -1, names_to = 'search', values_to = 'weight')
		aw <- aw %>% 
				mutate(search = factor(search), model = factor(model)) %>%
				mutate(search = fct_relevel(search, c('default','parbounds','sann','sann_pb','min')),
						model = fct_relevel(model, c('full2','full3','suicide4','suicide6','CID2','CID4'))
					) %>%
				mutate(model = fct_recode(model, 
								'Full 2' = 'full2',
								'Full 3' = 'full3',
								'Suicide 4' = 'suicide4',
								'Suicide 6' = 'suicide6'
								)
						)
		names(aw) <- gsub('(^\\\\w)', '\\\\U\\\\1', names(aw), perl = T)
		ggplot(aw, aes(Search, Weight, fill = Model)) +
			geom_col(width = 0.95) + # Nearly touching
			theme_classic() +
			theme(
				axis.text = element_text(color = 'black'),
				axis.ticks = element_line(color = 'black')
			) +
			scale_x_discrete(labels = c('Default', 'Bounded\\n parameters', 'Simulated\\n annealing','Sim. ann.\\n + bounds','Minimum\\n AICc')) +
			scale_fill_brewer(palette = 'Set1')
		# Looks a bit better in an external plotting device:
			ggsave(filename = paste0('four_searches_weights_',Sys.Date(),'.pdf'), width = 8, height = 6, units = 'in')

### (4) Clean up
	rm(hs, aw, all_out)

#####################################################################
### Figure 2: support for models across starting parameter values ###
#####################################################################

### (1) Load data from variable starting value analyses; each is a 6-element list (for 6 models), with each element having 48 model fits 
	load('SI8_list_modfits_pbounds_varstart.Rdata') # Loads 'bp_list'
	load('SI9_list_modfits_simann_pbounds_varstart.Rdata') # Loads 'sa_bp_list'
	
### (2) Summarize support for each model at each starting value
	# (a) Parameter bounds only
		aicc_pb <- t(sapply(bp_list, function(x) unlist(lapply(x, function(y) y$AICc))))
	# (b) Simulated annealing + parameter bounds
		aicc_s <- t(sapply(sa_bp_list, function(x) unlist(lapply(x, function(y) y$AICc))))
	
### (3) Choose min AICc / max lnL for each combination of model and starting values, across the two types of searches
	sum(aicc_pb > aicc_s) # 192, the number of times across 48 * 6 = 288 searches that sann=T improved the search (about 2/3 of the time)
	best_aicc <- aicc_s * (aicc_pb >= aicc_s) + aicc_pb * (aicc_pb < aicc_s) # Pulling those that were best, then merging them into a single matrix
	best_wgt <- apply(best_aicc, 2, aic.weight, rd.digits = NULL) # Don't round weights, as this causes some to sum > 1 (screws up graphics)
	colnames(best_wgt) <- paste0('run', seq(ncol(best_wgt)))

### (4) Graph weights
	# (a) Make into tibble
		best_wgt <- data.frame(best_wgt)
		best_wgt$model <- rownames(best_wgt)
		best_wgt <- as_tibble(best_wgt) %>%
			select(model, everything())
	# (b) Stretch for plotting
		bw <- pivot_longer(best_wgt, -1, names_to = 'run', values_to = 'weight')
		bw <- bw %>%
			mutate(
				run = factor(run), 
				model = factor(model)
			) %>%
			mutate(
				model = fct_relevel(model, c('full2','full3','suicide4','suicide6','CID2','CID4'))
			) %>%
			mutate(model = fct_recode(model, 
				'Full2' = 'full2',
				'Full3' = 'full3',
				'Suicide4' = 'suicide4',
				'Suicide6' = 'suicide6'
				)
			)
		names(bw) <- gsub('(^\\\\w)', '\\\\U\\\\1', names(bw), perl = T)
	# Put in order of the largest peak weight, regardless of model:
		pw <- bw %>%
			group_by(Run) %>%
			mutate(peakwgt = max(Weight))
	# Cannot figure out how to order bars better than this, so plot and later fix by hand:
	# Round values back to 3 digits, so that working with the figure manually is less of a pain.
		ggplot(pw, aes(fct_reorder(Run, peakwgt, .desc = T), round(Weight,3), fill = Model)) +
			geom_col(width = 0.95) + # Touching
			scale_fill_brewer(palette = 'Set1') +
			scale_y_continuous(name = 'AICc weight', expand = c(0,0)) +
			xlab('Starting values') +
			theme_classic() +
			theme(
				axis.text = element_text(color = 'black'),
				axis.text.x = element_blank(),
				axis.ticks = element_line(color = 'black'),
				axis.ticks.x = element_blank()
			)
		# Note: looks rough in Rstudio. Much better is stretched horizontally in an external plotting device:
			ggsave(filename = paste0('startingparams_weights_',Sys.Date(),'.pdf'), width = 5, height = 2, units = 'in')

### (5) Clean up
	rm(aicc_pb, aicc_s, best_aicc, best_wgt, bp_list, bw, pw, sa_bp_list)

###########################################################################
### Figure 3: effect of focal character states on diversification rates ###
###########################################################################

### (1) Load data
	# (a) Load ancestral-state estimates, obtained by running 'MarginReconfHiSSE()' with 'pars'. 
	# 	argument from MLEs of each model, hidden.states = 2, and all other default arguments
		load('SI12_ancstates_full2suicide4.Rdata') # Loads 'final_recon'
		final_recon # Ancestral-state estimates (and other information) for 'full2' and 'suicide4'
	# (b) Phenotypic data
		dat <- read.csv('Tonini.phyto.dat.1579sp.csv')

### (2) Plotting model-averaged ancestral-state estimates on tree, after new model fitting in this paper
	# plot.hisse.states(final_recon, rate.param = "net.div", show.tip.label = FALSE,
	#		rate.colors = c('black','white'), state.colors = c('blue', 'gray', 'red'),
	#		edge.width = 0.5, width.factor = 1.5) 
	# Nearly identical to Figure 2 in Tonini et al. 2020, so not worth putting in manuscript
			
### (3) Figure 3: model-averaged values, across Full2 and Suicide4, the models accounting for >98% of model weight
	# (a) Calculate tip estimates for diversification rates
		rates <- lapply(final_recon, function(pp) pp$tip.mat[,-1] %*% t(pp$rates.mat)) # Assuming multiple models with non-zero weights
		aicc_wgt <- aic.weight(unlist(map(final_recon, ~.$aic))) # Extract from the model fit
		rates <- rates[[1]] * aicc_wgt[[1]] + rates[[2]] * aicc_wgt[[2]]
		div_dat <- as_tibble(cbind(dat, rates)) %>%
					mutate(states = factor(states)) %>%
					mutate(dstate = fct_recode(states, 'non-phytotelm' = '0', 'phytotelm' = '1'))
		plot_dat <- div_dat %>%
			pivot_longer(cols = c(net.div, speciation, extinction), names_to = 'type', values_to = 'value') %>%
			mutate(type = factor(type, c('net.div', 'speciation', 'extinction'))) # Just looking at these three rates, as the typical diversification rates
	# (b) Plot histograms of net diversification rate, speciation rate, and extinction rate
		ggplot(plot_dat, aes(value, fill = dstate, color = dstate)) +
			geom_histogram(binwidth = 0.002) + 
			scale_x_continuous(
				name = 'Rate (lineages / Myr)',
				breaks = c(0, 0.05, 0.1)
			) +
			scale_y_continuous(name = 'Number of species', expand = c(0,0)) +
			scale_fill_manual(
				name = 'Focal state', 
				values = c('non-phytotelm' = 'gray', 'phytotelm' = 'black'),
				labels = c('Non-phytotelm', 'Phytotelm')
			) +
			scale_color_manual(
				name = 'Focal state', 
				values = c('non-phytotelm' = 'gray', 'phytotelm' = 'black'),
				labels = c('Non-phytotelm', 'Phytotelm')
			) +
			theme_classic() +
			theme(
				axis.title = element_text(size = 8),
				axis.text = element_text(size = 6, color = 'black'),
				axis.ticks = element_line(color = 'black'),
				legend.title = element_text(size = 8),
				legend.text = element_text(size = 6),
				strip.text = element_text(size = 6)
			) +
			facet_wrap(~type, nrow = 1, scales = 'free_y')
	# (c) Much nicer as an exported object
		ggsave(paste0('modavg_rateshistogram_grey_',Sys.Date(),'.pdf'), width = 5, height = 2)
	
### (4) Clean up workspace
	rm(rates, div_dat, plot_dat, aicc_wgt, dat, final_recon)


