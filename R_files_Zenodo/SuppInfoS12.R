### Full R code for subsampling rate simulations
###	From Evolution submission: "Testing adaptive radiation: a new test with Madagascar frogs"

### All analyses done on R version 4.0.2

### (0) Load relevant libraries
	library(tidyverse) # Version 1.3.0
	library(devtools) # Version 2.3.0
	# devtools::install_github('geomorphR/geomorph', ref = 'Develop') # Body-size analyses will only work with this version of 'geomorph' (Fall 2020)
	library(geomorph) # Version 3.3.1
	library(phytools) # Version 0.7-70
	library(geiger) # Version 2.0.7
	library(motmot) # Version 2.1.3
	library(parallel) # Version 4.0.2

### (1) Load tree
	tree <- read.tree('SuppInfoS13.3449-species.phylogeny.tre')	
	# Make tree unit length
		tree$edge.length <- tree$edge.length / max(nodeHeights(tree))
		max(nodeHeights(tree)) # Looks good

### (2) Function for rate simulations
	source('SuppInfoS11.helper.functions.R') # Loads function 'sigmamult_subsample_sims', which does most of the work here

### (3) Set up looped simulations
	# Set simulation parameters
		p <- c(1, 3, 9) # Number of different traits
		corp <- c(0, 0.5) # Consider no correlation and one of 0.5
		rates <- seq(1,2.5,0.25) # Seven rate ratios
		tree_size <- list(small = c(25, 55), large = c(36, 181))

	# Idea: loop over number of traits, correlation, and tree size. 5 total loops, since p == 1 has no correlation.
		conditions <- list()
		for (i in seq(tree_size)){
			for (j in seq(corp)){
				for (k in seq(p)){
					ind <- 6*i + 3*j + k - 9
					tmp <- list()
					tmp[[1]] <- tree_size[[i]]
					tmp[[2]] <- corp[j]
					tmp[[3]] <- p[k]
					conditions[[ind]] <- tmp
					names(conditions)[ind] <- paste0(
												names(tree_size)[i],'tree_',
												paste0(corp[j],'corr_'),
												paste0(p[k],'traits')
												)
				}
			}
		}
		for (i in seq(conditions)) names(conditions[[i]]) <- c('tree_size', 'correlation', 'ntraits')
		conditions[grep('0.5corr_1trait', names(conditions))] <- NULL # Remove single-trait 'correlation' setup
		names(conditions) # Looks good
		
	# Now full analyses with 1000 iterations. Note that 'mclapply()' from package 'parallel' is called within the function 'sigmamult_subsample_sims()'
		cores <- 32 # Adjust for your computer; I ran this on a 32-core cluster node
		nsims <- 1000 # Number of simulation replicates run for paper; adjust down to test code
		all_out <- lapply(conditions, 
					function(y) sigmamult_subsample_sims(tree, node = 3466, subsamp = y$tree_size, 
						p = y$ntraits, R = y$correlation, rates, iter = nsims, para.proc = TRUE, 
						cores = cores, pv = TRUE)
					)
		names(all_out) <- names(conditions)
		save(all_out, file = paste0('subsampling_rate_sims_',Sys.Date(),'.Rdata'))
	# Clean up workspace for processing output below
		rm(conditions, tmp, tree, tree_size, corp, i, j, k, p, rates, ind)		
	# Load 'all_out' to skip rerunning the simulations (e.g. to test code)
		load('SuppInfoS14.subsamplingsims.output.Rdata') # Loads 'all_out', from line 58 above

### (4) Plot power curves
	# (a) First extract summary tables, put into a list, and extract only power and analysis conditions
		pow_out <- map(all_out, ~.$summary)
		names(pow_out) <- str_replace(names(pow_out), 'tree', '') %>%
							str_replace('corr', '') %>%
							str_replace('traits', '')
		for (i in seq(pow_out)){
			pow_out[[i]]$rate <- str_replace(rownames(pow_out[[i]]), 'rr_', '')
			pow_out[[i]]$analysis <- names(pow_out)[i]
			pow_out[[i]] <- pow_out[[i]] %>% as_tibble() %>% select(analysis, rate, power)
		}
	# (b) Write function to loop 'full_join' on a list greater than length 2
		multi_join <- function(x){
			out <- full_join(x[[1]], x[[2]])
			for (i in 3:length(x)) out <- full_join(out, x[[i]])
			return(out)	
		}
	# (c) Now join them all 
		all_pow <- multi_join(pow_out)		
		all_pow # 70 rows (7 rates, 10 other conditions: 3 trait numbers, two tree sizes, and two correlations, minus for the single trait)				
	# (d) Format 'analysis' and 'rate' variables
		all_pow <- all_pow %>%
					separate(analysis, into = c('treesize','corr','ntrait'), sep = '_') %>%
					mutate(
						treesize = factor(treesize), 
						corr = factor(corr), 
						ntrait = factor(ntrait),
						rate = as.double(rate)
					) 
	# (e) Function to reduce subsequent code
        power_curve_theme <- function(){
			list(
				theme_classic(),
				scale_x_continuous(
					name = 'Rate ratio',
					breaks = seq(1, 2.5, 0.25)
				),
				scale_y_continuous(
					name = 'Power',
					limits = c(0, 1),
					breaks = seq(0, 1, 0.2),
					expand = c(0,0)		
				),
				theme(
					axis.text.y = element_text(angle = 90, hjust = 0.5),
					axis.text = element_text(color = 'black'),
					axis.ticks = element_line(color = 'black')
				)			
			)
        }
	# (f) Concatenate treesize and corr to facilitate plotting
		pow_plot <- all_pow %>% unite('size_cor', treesize, corr, remove = F)
	# (g) Plot power curve 
		ggplot(pow_plot, aes(rate, power)) +
			geom_hline(yintercept = 0.05, linetype = 'dashed') +
			geom_line(aes(group = size_cor, linetype = corr)) + 
			geom_point(aes(fill = treesize), shape = 21, size = 2) +
			scale_fill_manual(values = c('black','white')) + 
			scale_linetype_manual(values = c('solid', 'dotted')) +
			power_curve_theme() +
			facet_wrap(~ntrait, nrow = 1, scales = 'free_y') # Y on same scale, but this puts a scale next to each facet
	# (h) A bit nicer as a separate object		
		ggsave(filename = paste0('subsample_power_curves_',Sys.Date(),'.pdf'), width = 10, height = 3, units = 'in')
	# (d) Clean up space
		rm(all_pow, pow_out, pow_plot, i, power_curve_theme)
		
### (5) Make matrix for plotting rate estimates
	# (a) First extract summary tables, put into a list, and extract only power and analysis conditions
		par_out <- map(all_out, ~.$summary)
		names(par_out) <- str_replace(names(par_out), 'tree', '') 
		for (i in seq(par_out)){
			par_out[[i]]$rate <- str_replace(rownames(par_out[[i]]), 'rr_', '')
			par_out[[i]]$analysis <- names(par_out)[i]
			par_out[[i]] <- par_out[[i]] %>% as_tibble() %>% select(analysis, rate, starts_with('focal'), starts_with('other'))
		}
	# (b) Now join them all 
		all_par <- multi_join(par_out) # Function 'multi_join()' defined above in lines 82â€“86
		all_par # 70 rows (7 rates, 10 other conditions: 3 trait numbers, two tree sizes, and two correlations, minus for the single trait)				
	# (c) Format 'analysis' and 'rate' variables
		all_par <- all_par %>%
					separate(analysis, into = c('treesize','corr','ntrait'), sep = '_') %>%
					mutate(
						treesize = factor(treesize, c('small','large')), # Ensures they plot in this order
						corr = factor(corr, c('0corr','0.5corr')), # Ensures they plot in this order
						ntrait = factor(ntrait),
						rate = as.double(rate)
					) 					
	# (d) Pivot longer to allow plotting of focal and other on same figure
		long_par <- all_par %>%
					pivot_longer(ends_with('mean'), names_to = 'group', values_to = 'mean') %>%
					pivot_longer(ends_with('sd'), names_to = 'group2', values_to = 'sd') %>%
					separate(group, into = c('group','dummy')) %>%
					separate(group2, into = c('group2', 'dummy2')) %>%
					select(-starts_with('dummy')) %>%
					filter(group == group2) %>%
					select(-group2)
	# (e) Function to reduce subsequent code
		param_theme <- function(){
			list(
				theme_classic(),
				scale_x_continuous(
					name = 'Simulated rate ratio',
					breaks = seq(1, 2.5, 0.5)
				),
				scale_y_continuous(
					name = 'Estimated rate',
					breaks = seq(0.5, 3.0, 0.5)
				),
				theme(
					axis.text.y = element_text(angle = 90, hjust = 0.5),
					axis.text = element_text(color = 'black'),
					axis.ticks = element_line(color = 'black')
				)			
			)
		}
	# (f) Plot simulated rate by estimated rate, with error bars via 'sd'. 
		ggplot(long_par, aes(rate, mean)) +
			geom_hline(yintercept = 1.0, linetype = 'solid') + # Line at 1.0 for 'other' rate
			geom_abline(slope = 1.0, intercept = 0.0, linetype = 'solid') + # 1:1 line
			geom_linerange(aes(ymin = mean - sd, ymax = mean + sd)) + # Standard deviation bars
			geom_point(aes(fill = group), shape = 21, size = 2) + # Means of rate estimates
			scale_fill_manual(values = c('black','white')) + 
			param_theme() +
			coord_fixed() + 
			facet_grid(treesize ~ corr + ntrait) # Do both number of traits and trait correlation in the same dimension
	# (g) A bit nicer as a separate object		
		ggsave(filename = paste0('subsample_param_curves_',Sys.Date(),'.pdf'), width = 12, height = 4, units = 'in')