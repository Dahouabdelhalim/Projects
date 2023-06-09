### R code for data analysis, written by Daniel Moen for the paper: "Phylogenetic analysis 
#	of adaptation in comparative physiology and biomechanics: overview and a case study of 
#	thermal physiology in treefrogs" by Moen, Cabrera-Guzmán, Caviedes-Solis, González-Bernal, and Hanna



##########################################################################################
### (0) Load libraries and custom functions
##########################################################################################

library(OUwie) # Our principal package for OU analysis. Version 2.6. 
library(ouch) # An alternative package for OU analysis. Version 2.17.
library(phytools) # Used here for stochastic mapping. Version 0.7-70.
library(tidyverse) # Suite of packages for data manipulation. Key included packages used 
# are dplyr, tidyr, ggplot2, and purrr. Version 1.3.0.
library(lemon) # Used for plotting TPCs. Version 0.4.5
filter <- dplyr::filter # Force the function, as the same name is used by another package
source('AppendixS5_custom_functions.R') # All custom-written functions for analyses. Includes:
# (1) 'OU_singlebranch_sim()' for simulations that produce Fig. 1
# (2) 'aicc_wgt_calc()' for calculating AICc weights
# (3) 'OUwie_hylids()' for looping OU analyses
# (4) 'paraboot_OUwie()' for conducting parametric bootstrap model tests with package OUwie
# These functions are all described (e.g. arguments, output) within the source file. 



##########################################################################################  
### (1) Run Fig. 1 simulations ----
##########################################################################################

# (a) Alpha calculator. Inputs are time at which regimes switch ('st') and number of half-lives 
# 	traversed in final time period ('num_hl'). Returns the exact alpha needed to produce the desired expected value.
	calc_alpha <- function(st, num_hl) return(num_hl * log(2) / (1 - st))

# (b) Run a simulation where a species switches late (sw = 0.8) but has a strong alpha (calculate).
	# Conditions
		nsims <- 25
		switch_sim <- 0.8
		sigma_sim <- c(2, 2)
		time_sim <- 1
		opt_alpha <- calc_alpha(switch_sim, 2) # Necessary alpha to move about halfway to new optimum over last 20% of simulation
	# Simulate
		set.seed(1978)
		strong_alpha <- OU_singlebranch_sim(alpha = rep(opt_alpha,2), sigma_sq = sigma_sim, 
							theta = c(0,1), num_steps = 100, switch_cond = switch_sim, 
							nsim = nsims, time = time_sim)
	# Graph  
		theta_out <- unique(strong_alpha$thetas)
		time_thresholds <- c(0, max(strong_alpha[strong_alpha$thetas==theta_out[1], 'times']), 
							min(strong_alpha[strong_alpha$thetas==theta_out[2], 'times']), 
							max(strong_alpha[strong_alpha$thetas==theta_out[2], 'times']))
		theta_draw_strong <- tibble(regime = c('1st', '2nd'), 
							 x = time_thresholds[c(1,3)], 
							 xend = time_thresholds[c(2,4)],
							 y = unique(strong_alpha$thetas),
							 yend = unique(strong_alpha$thetas))
		strong_alpha <- as_tibble(strong_alpha) %>%
		  pivot_longer(-1:-3)
		strong_alpha %>%
		  ggplot(aes(times, value, group = name)) +
			geom_line() + 
			geom_segment(data = theta_draw_strong, aes(x = x, xend = xend, y = y, 
				yend = yend, group = 'regime'), color = 'red')
		
# (c) Run a simulation where a species switches early (sw = 0.2) but has a weak alpha (calculate)
	# Set parameters
		switch_sim <- 0.2
		opt_alpha <- calc_alpha(switch_sim, 2) # Necessary alpha to move about halfway to new optimum over last 1/4 of simulation
	# Run simulation
		set.seed(1978)
		weak_alpha <- OU_singlebranch_sim(alpha = rep(opt_alpha,2), sigma_sq = sigma_sim,
						theta = c(0,1), num_steps = 100, switch_cond = switch_sim, nsim = nsims, time = time_sim)
		theta_out <- unique(weak_alpha$thetas)
		time_thresholds <- c(0, max(weak_alpha[weak_alpha$thetas==theta_out[1], 'times']),
							min(weak_alpha[weak_alpha$thetas==theta_out[2], 'times']),
							max(weak_alpha[weak_alpha$thetas==theta_out[2], 'times']))
		theta_draw_weak <- tibble(regime = c('1st', '2nd'), 
						   x = time_thresholds[c(1,3)], 
						   xend = time_thresholds[c(2,4)],
						   y = unique(weak_alpha$thetas),
						   yend = unique(weak_alpha$thetas))
		weak_alpha <- as_tibble(weak_alpha) %>%
		  pivot_longer(-1:-3)
		weak_alpha %>%
		  ggplot(aes(times, value, group = name)) +
		  geom_line() + 
		  geom_segment(data = theta_draw_weak, aes(x = x, xend = xend, y = y, yend = yend, group = 'regime'), color = 'red')
	  
# (d) Combine output matrices for figure
	weak_alpha <- mutate(weak_alpha, sim_run = '(A) Earlier shift, weaker alpha')
	strong_alpha <- mutate(strong_alpha, sim_run = '(B) Later shift, stronger alpha')
	theta_draw_weak <- mutate(theta_draw_weak, sim_run = '(A) Earlier shift, weaker alpha')
	theta_draw_strong <- mutate(theta_draw_strong, sim_run = '(B) Later shift, stronger alpha')
	all_alpha <- full_join(weak_alpha, strong_alpha)
	all_theta <- full_join(theta_draw_weak, theta_draw_strong)
	all_alpha %>%
	  ggplot(aes(times, value, group = name)) +
	  geom_line(color = 'black', alpha = .2, size = 0.5) + 
	  geom_segment(data = all_theta, aes(x = x, xend = xend, y = y, yend = yend, group = 'regime'), 
				   color = 'black',
				   linetype = 'longdash',
				   size = 0.75) +
	  geom_segment(aes(x = 1.07, xend = 1.02, y = 0.75, yend = 0.75), 
				   arrow = arrow(length = unit(0.25, 'cm'), type = 'open'),
				   lineend = 'round',
				   linejoin = 'round',
				   size = 0.5) +
	  scale_x_continuous(name = 'Time', breaks = seq(0,1,0.2), limits = c(0, 1.1), expand = c(0,0)) +
	  scale_y_continuous(name = 'Trait value', breaks = seq(0,1,0.25), limits = c(-0.2, 1.2), expand = c(0,0)) +
	  theme_classic() +
	  theme(
		axis.title = element_text(size = 20),
		axis.text = element_text(size = 16, color = 'black'),
		strip.text = element_text(size = 20, hjust = 0),
		strip.background = element_blank()
	  ) +
	  facet_wrap(~sim_run, nrow = 2, scales = 'free')

# (e) Save (looks much nicer this way, as I chose plotting parameters for PDF, not for displaying within RStudio
	dims <- c(8, 8) # Width by height, in inches
	dir.create('Figures') # Make new output directory for figures
	fig_wd <- './Figures/'
	ggsave(paste0(fig_wd,'Fig1_OUsims_',Sys.Date(),'.pdf'), width = dims[1], height = dims[2], units = 'in')
    
# (f) Clean up workspace
	rm(all_alpha, all_theta, strong_alpha, weak_alpha, theta_draw_strong, theta_draw_weak,
	   dims, nsims, opt_alpha, sigma_sim, switch_sim, theta_out, time_sim, calc_alpha,
	   time_thresholds, OU_singlebranch_sim)



##########################################################################################
### (2) Fit models for calculating L80 (and related variables)
##########################################################################################

# (a) Load intraspecific data
	all_dat <- tibble(read.csv('AppendixS1_hylid_jumping_data.csv'))
	# Put species names in the order we would like for plotting thermal performance curves
	all_dat <- all_dat %>% mutate(species = factor(species, c(
		'Acris blanchardi', 'Pseudacris crucifer', 'Pseudacris fouquettei',
		'Hyla arenicolor', 'Hyla avivoca', 'Hyla cinerea',
		'Charadrahyla nephila', 'Exerodonta abdivita', 'Ptychohyla zophodes',
		'Smilisca baudinii', 'Smilisca cyanosticta', 'Tlalocohyla smithii'
	))) 
		
# (b) Loop across all species. In each comparison, extract the AICc and the AICc weight
	# Function for assessing model support (we also use 'aicc_wgt_calc()', previously loaded)
		aicc_calc <- function(n, k, aic) return(aic + 2*k*(k+1) / (n - k - 1))
	# Set up outputs
		mod_sum <- data.frame(matrix(nrow = nlevels(all_dat$species), ncol = 5))
		colnames(mod_sum) <- c('species','poly.aicc','gaus.aicc','poly.wgt','gaus.wgt')
		mod_coefs <- data.frame(matrix(nrow = nlevels(all_dat$species), ncol = 6)) 
		colnames(mod_coefs) <- c('species', paste0('L',c(70,80,90)), 'peak.temp', 'breadth')
		mod_sum[,1] <- mod_coefs[,1] <- levels(all_dat$species)
	# Run loop across species	
		dir.create('TPCs') # Make a directory for keeping model fits for each species
		for (i in levels(all_dat$species)){
		# Set data
			ind <- which(mod_sum$species == i)
			tmp <- all_dat %>% filter(species == i)
		# Estimate models
			tpo <- with(tmp, lm(sv ~ poly(Tb, 2)))
			tga <- with(tmp, nls(sv ~ a*exp(-0.5 * (abs(Tb - b) / c)^2), start=list(a=1,
					b=32,c=20), control = list(warnOnly = T)))
			mods <- list(polynomial2 = tpo, Gaussian = tga)
		# Save species-specific model fits for later predicting values in figure code.
			save(mods, file = paste0('./TPCs/model.fits.',i,'.',Sys.Date(),'.Rdata'))
		# Extract results
			aic <- unlist(lapply(mods, AIC))		
			aicc <- aicc_calc(nrow(tmp), c(4,4), aic) # AICc values
			mod_sum[ind,2:3] <- aicc
			mod_sum[ind,4:5] <- aicc_wgt_calc(aicc) # AICc weights
		# Extract optimal model, then find L70, L80, L90, peak.temp, and breadth (temp. difference between peak 
		#	predicted temperature at peak velocity and L80)
			# Predict the whole curve across the full range of temps. 
				x <- mod_coefs[ind,] # Easier bookeeping
				pmod <- mods[[which.min(aicc)]]
				pred <- data.frame(Tb = seq(0.1,35,0.01)) # Derive curve from 0.1–35 Celsius
				pred$perf <- predict(pmod, pred)
			# Clip values above the peak (since we don't want the H80, for example)
				pred <- pred[1:which.max(pred$perf),]
			# Now extract results	
				l <- c(0.7, 0.8, 0.9)
				for (j in l){
					x[, 1 + which(l %in% j)] <- pred$Tb[which.min(abs(j - pred$perf))]
				}
				x[,'peak.temp'] <- pred$Tb[which.max(pred$perf)]
				x[,'breadth'] <- x$peak.temp - x$L80
				x[,-1] <- round(x[,-1], 2)
				mod_coefs[ind,] <- x
		}
	# Note that some of our Gaussian fits did not fit particularly well, producing error messages. 	
		mod_sum # Statistical support for models across all species
		mod_coefs # Estimated values for species
	# Clean workspace
		rm(ind, tmp, tpo, tga, mods, aic, aicc, x, pmod, pred, l, i, j)	

# (c) Calculate a set of predicted values for plotting all types of curves 
	mod_sum <- tibble(mod_sum) # Format model-fitting results
	# Make a loop to save only values within range of observations per species
	sp <- unique(mod_sum$species)
	for(i in sp){
		ind <- which.min(mod_sum[mod_sum$species == i,2:3])
		load(paste0('./TPCs/',grep(i, dir('./TPCs/'), value = T))) # Loads 'mods'
		pmod <- mods[[ind]]
		pred <- data.frame(Tb = seq(0.1,40,0.1))
		pred$perf <- predict(pmod, pred)
		rg <- all_dat %>% filter(species == i) %>% select(Tb) %>% range
		pred <- tibble(pred) %>% filter(Tb >= rg[1] & Tb <= rg[2])
	# Reshape 'pred' for output
		zone <- all_dat %>% filter(species == i) %>% pull(region) %>% unique()
		pred <- pred %>% 
			mutate(species = i, region = zone) %>%
			rename(sv = perf) %>%
			select(species, region, Tb, sv)
		if(i == sp[1]){
			pred_out <- pred
		} else {
			pred_out <- bind_rows(pred_out, pred)
		}
	}
	# Briefly examine our coordinates for plotting lines
		pred_out 
	# Clean workspace
		rm(sp, ind, pmod, pred, rg, zone, mods, i)



##########################################################################################
### (3) Plot TPCs ------------
##########################################################################################

# (a) Subset data for calculating species means:
	#	We removed some data (19 of 432 datapoints) whose observed 
	#	temperature was closer to the next experimental temperature than the intended temperature. 
	#	These data were still used in regression models, as we used observed temperatures for
	#	fitting curves. However, we excluded them for plotting species means at different 
	#	temperatures. These means were NOT used as data; they are only plotted for 
	#	visualization of species means. Note that 'eTb' for many of these reflect an adjusted
	#	'eTb', meaning that in essence these individuals have two 'eTb', whereas the observed
	#	temps are different (see full dataset). These datapoints were thus removed so as to not
	#	pseudoreplicate individuals when calculating species means for plotting. 
	mdat <- all_dat %>% filter(plot_temps) # Variable 'plot_dat' is logical, meaning we keep on 'TRUE'

# (b) Temperate curves: format data
	dat <- all_dat %>% 
			filter(region == 'temperate')
	reg <- as.character(unique(dat$region))
	sp_means <- mdat %>% 
		filter(species %in% unique(dat$species)) %>%
		mutate(species = factor(species)) %>%
		group_by(species, eTb) %>%
		summarise_at(vars(Tb:sv), mean, na.rm = T)
	spec_labs <- tibble(mod_coefs) %>%
		filter(species %in% sp_means$species) %>%
		mutate(species = factor(species)) # This last one to remove the extra tropical species levels
	pred_out_cr <- pred_out %>% 
		filter(region == reg) %>% 
		mutate(species = factor(species))
			
# (c) Temperate curves: plot
	ggplot(dat, aes(Tb, sv)) +
		# Plot curves based on predicted values for optimal models
		geom_line(data = pred_out_cr, size = 1.5, color = 'grey70') +
		# Intraspecific points
		geom_point(size = 2) +
		# Species means
		geom_point(data = sp_means, size = 3, shape = 21, fill = 'white') +
		# Add species names (instead of labeling the facets; see 'theme' adjustments below
		geom_text(data = spec_labs, aes(x = 35, y = 0.40, label = species), hjust = 'inward', 
			size = 5, fontface = 'italic') +
		scale_x_continuous('Body temperature (°C)', breaks = c(8,14,20,26,32), limits = c(5,36.5), expand = c(0,0)) +
		scale_y_continuous('Standardized jumping velocity', breaks = seq(0.4,1.0,0.2), limits = c(0.35, 1.05), expand = c(0,0)) +
		theme_classic() +
		theme(
			axis.title = element_text(size = 20),
			axis.text = element_text(size = 14, color = 'black'),
			axis.ticks = element_line(color = 'black', size = 0.75, length),
			axis.ticks.length = unit(0.0625, 'in'),
			strip.text = element_blank(),
			strip.background = element_blank()
		) +
		facet_rep_wrap(~species, nrow = 2, repeat.tick.labels = 'none')
	# Now save plot (plotting parameters chosen for PDF, not RStudio plotting device)
		ggsave(paste0(fig_wd, 'Fig3_temperate_TPCs_',Sys.Date(),'.pdf'), height = 8, width = 12)
			
# (d) Tropical curves: format data
	dat <- all_dat %>% 
			filter(region == 'tropical') %>% 
			mutate(species = factor(species)) # Alphabetical order is fine
	reg <- as.character(unique(dat$region))
	sp_means <- mdat %>% 
		filter(species %in% unique(dat$species)) %>%
		mutate(species = factor(species)) %>%
		group_by(species, eTb) %>%
		summarise_at(vars(Tb:sv), mean, na.rm = T)
	spec_labs <- tibble(mod_coefs) %>%
		filter(species %in% sp_means$species) %>%
		mutate(species = factor(species)) # This last one to remove the extra tropical species levels
	pred_out_cr <- pred_out %>% 
		filter(region == reg) %>% 
		mutate(species = factor(species))
		
# (e) Tropical curves: plot
	ggplot(dat, aes(Tb, sv)) +
		# Plot curves based on predicted values for optimal models
		geom_line(data = pred_out_cr, size = 1.5, color = 'grey70') +
		# Intraspecific points
		geom_point(size = 2) +
		# Species means
		geom_point(data = sp_means, size = 3, shape = 21, fill = 'white') +
		# Add species names (instead of labeling the facets; see 'theme' adjustments below
		geom_text(data = spec_labs, aes(x = 35, y = 0.40, label = species), hjust = 'inward', 
			size = 5, fontface = 'italic') +
		scale_x_continuous('Body temperature (°C)', breaks = c(8,14,20,26,32), limits = c(5,36.5), expand = c(0,0)) +
		scale_y_continuous('Standardized jumping velocity', breaks = seq(0.4,1.0,0.2), limits = c(0.35, 1.05), expand = c(0,0)) +
		theme_classic() +
		theme(
			axis.title = element_text(size = 20),
			axis.text = element_text(size = 14, color = 'black'),
			axis.ticks = element_line(color = 'black', size = 0.75, length),
			axis.ticks.length = unit(0.0625, 'in'),
			strip.text = element_blank(),
			strip.background = element_blank()
		) +
		facet_rep_wrap(~species, nrow = 2, repeat.tick.labels = 'none')
	# Now save plot (plotting parameters chosen for PDF, not RStudio plotting device)
		ggsave(paste0(fig_wd, 'Fig4_tropical_TPCs_',Sys.Date(),'.pdf'), height = 8, width = 12)
  # Clean up workspace
		rm(reg, mod_coefs, mod_sum, all_dat, mdat, pred_out, pred_out_cr, sp_means, 
		   spec_labs, dat)


##########################################################################################
### (4) Load and format data for phylogenetic comparative analyses ----
##########################################################################################

# (a) Load data
	tree <- read.tree('AppendixS3_12taxon.mcc.tre') 
	dat <- read_csv('AppendixS2_species_means.csv')
	
# (b) Manually assign ancestral states to new consensus trees (temp-trop, temp-highland-lowland). States follow Moen et al. (2009). 
# Set the tropical-temperate internal states for OU2. 
	temptrop <- c('tropical',rep('tropical',4), rep('temperate',2), rep('tropical',2), rep('temperate',2))
	plot(tree, cex = 0.75, no.margin = T)
	nodelabels(temptrop) # 'temptrop' is in the correct node order
	tree_ou2 <- tree
	tree_ou2$node.label <- temptrop	# OUwie needs internal node labels for ancestral regimes
# Now set highland-lowland for OU3
	highlow <- c('lowland',rep('highland',3),'lowland',rep('temperate',2),rep('lowland',2),rep('temperate',2))
	plot(tree, cex = 0.75)
	nodelabels(highlow) #  'highlow' is in the correct node order
	tree_ou3 <- tree
	tree_ou3$node.label <- highlow
# Compile both to run in function
	trees_ou <- list(OU2 = tree_ou2,
					OU3 = tree_ou3)
	class(trees_ou) <- 'multiPhylo'
# (c) Clean up
	rm(tree_ou2, tree_ou3, temptrop, highlow)



##########################################################################################
### (5) OUwie analyses ---------------
##########################################################################################

# (a) Run models: all data, with and w/o SEs for CTmin (first result is with SEs) ----
	traits <- names(dat)[c(5,5,7:11)]	
	SEs <- c(TRUE, rep(FALSE, length(traits) - 1))
	dir.create('Results') # Make new output directory for figures
	res_wd <- './Results/'
	out <- map2(traits, SEs, ~OUwie_hylids(dat, trees_ou, resp_var = .x, SE = .y, 
				wd_out = res_wd, OUwie_args = list(scaleHeight = T, algorithm = 'invert',
				quiet = T, warn = F, root.station = T)))
	names(out) <- paste0(traits, '_', ifelse(SEs,'wSE','noSE'))

# (b) Run L80 analyses (again) without Tlalocohyla smithii ----
	dat_11 <- dat %>%
		filter(species != 'Tlalocohyla smithii')
	tlac_wd <- './Results/noTlalocohyla_'
	# Save this first one
	out_11 <- OUwie_hylids(dat_11, trees_ou, resp_var = 'L80', SE = F, wd_out = tlac_wd, 
				get.root.theta = FALSE, OUwie_args = list(scaleHeight = T, 
				algorithm = 'invert', quiet = T, warn = F, root.station = T))
	# Do not save this second one, since I've already done it. 
	# Just running it to more easily compare data when making tables below.
	out_12 <- OUwie_hylids(dat, trees_ou, resp_var = 'L80', SE = F, wd_out = NULL, 
				get.root.theta = FALSE, OUwie_args = list(scaleHeight = T, 
				algorithm = 'invert', quiet = T, warn = F, root.station = T))
				
# (c) Plot results: species means and adaptive optima for OU2 (CTmin only)
	# Upload optima for OU2, CTmin
		theta <- read_csv(paste0('./Results/',grep('CTmin_withSE_model_parameters_',dir('./Results/'),value=T)))
		theta <- theta %>% 
			filter(model == 'temptrop') %>% 
			select(theta.temp_MLE, theta.high_MLE) %>%
			as.numeric()
		theta <- tibble(region =  c('temperate', 'tropical'), value = theta, variable = 'CTmin')	

	# Adjust data for plotting ----
		plot_dat <- dat %>%
			select(species, region, elevation, CTmin, L80) %>%
			pivot_longer(cols = c('CTmin', 'L80'), names_to = 'variable') %>%
			mutate(region = fct_recode(region, 'Temperate' = 'temperate', 'Tropical' = 'tropical'))
		theta_plot <- mutate(theta, 
			region = fct_recode(region, 'Temperate' = 'temperate', 'Tropical' = 'tropical'))  
		my_labeller <- as_labeller(c(CTmin = '(A) CTmin', L80 = '(B) L80'))

	# Plot ----
		ggplot(plot_dat, aes(region, value)) +
		  geom_boxplot(aes(fill = region), coef = 100, show.legend = F, na.rm = T) + # Extreme 'coef' to ensure all data are included in the whiskers
		  geom_point(data = theta_plot, aes(color = region), show.legend = F, shape = 8, size = 4) +
		  scale_fill_manual(values = c(Temperate = 'grey70', Tropical = 'black')) +
		  scale_color_manual(values = c(Temperate = 'black', Tropical = 'gray60')) +
		  scale_y_continuous(limits = c(-5,30), expand = c(0,0), breaks = seq(-5,30,5)) +
		  ylab(quote(Temperature~~(degree*C))) +
		  theme_classic() +
		  theme(
			axis.title.y = element_text(size = 16),
			axis.title.x = element_blank(),
			axis.text.y = element_text(size = 12, color = 'black'),
			axis.text.x = element_text(size = 16, color = 'black', vjust = -0.5),
			axis.ticks.y = element_line(color = 'black'),
			axis.ticks.length.y = unit(5,'points'),
			axis.ticks.length.x = unit(0,'cm'),
			strip.text = element_text(size = 20, hjust = 0.5, vjust = 2, face = 'bold'),
			strip.background = element_blank(),
			panel.spacing = unit(13, 'points')
		  ) + 
		  facet_wrap(~variable, label = my_labeller)

	# Save PDF (plotting parameters optimized for PDF, not RStudio plotting device)
		wd <- c(8, 6)
		ggsave(paste0(fig_wd,'Fig5_bothvars_boxplot.',Sys.Date(),'.pdf'), width = wd[1], height = wd[2], unit = 'in')



##########################################################################################
### (6) Compile results into tables ----------------------
##########################################################################################

# (a) Short function to put together results for tables. Based on your input 'var', it
#	reads the results file based on the names in 'tabs', chooses the table variables,
#	then returns it in a format that allows joining with other similar tables.
	read_join_morph <- function(var, tabs){
		out <- read_csv(paste0('./Results/',grep(var, tabs, value = T))) %>% 
		  select(model:AICc, wgt) %>%
		  mutate(type = all_of(var)) %>%
		  unite(mod_var, model, type)
		return(out)
	}
# (b) Table 2: main MS results for CTmin and L80
	tabs <- grep('model_comparison',dir('./Results/'),value=T)
	tabs <- tabs[c(3,5)] # Keep those for Table 1
	ct <- read_join_morph('CTmin', tabs)
	table2 <- read_join_morph('L80', tabs) %>%
			full_join(ct)
	table2
			
# (c) Table S3: comparing top models for all data extracted from TPCs
	tabs <- grep('CTmin', dir('./Results/'), value = T, invert = T)
	tabs <- grep('model_comparison', tabs, value = T)
	tabs <- grep('noTlalocohyla', tabs, value = T, invert = T)
	names(tabs) <- seq(tabs) # Dummy names to make the next function work
	vars <- t(map_dfr(list(tabs), ~strsplit(.x, split = '_')))[,1]
	vars <- vars[c(1,5:2)]
	tableS3 <- read_join_morph(vars[1], tabs)
	for(i in vars[-1]) tableS3 <- read_join_morph(i, tabs) %>% 
	  full_join(tableS3)
	tableS3 <- tableS3 %>%
	          separate(mod_var, c('model','variable'), sep = '_')
	tableS3
	
# (d) Table S4: comparing CTmin with and w/o SEs
	tabs <- grep('CTmin', dir('./Results/'), value = T)
	tabs <- grep('model_comparison', tabs, value = T)
	ct <- read_join_morph('without', tabs)
	tableS4 <- read_join_morph('withSE', tabs) %>%
		full_join(ct)
	tableS4
	
# (e) Table S5: L80 results with and without Tlalocohyla smithii
	tableS5 <- out_11$mod_summary %>%
		select(model:AICc, wgt) %>%
		mutate(type = 'noTlalo') %>%
		unite(mod_var, model, type)
	tableS5 <- out_12$mod_summary %>%
		select(model:AICc, wgt) %>%
		mutate(type = 'withTlalo') %>%
		unite(mod_var, model, type) %>% 
		full_join(tableS5)
	tableS5

# (f) Clean workspace
	rm(list = ls(pattern = 'table'))
	rm(ct, dat_11, out, out_11, out_12, plot_dat, theta, theta_plot, i, my_labeller,
	   SEs, tabs, tlac_wd, traits, vars, wd, read_join_morph)

##########################################################################################
### (7) Parametric bootstrapping for supplementary materials and methods ----------------
##########################################################################################

### (a) Set up OUwie data frames ----
	dat_OU2_CTmin <- dat %>%
		select(species_phylo, region, CTmin, CTmin_SE) %>%
		dplyr::filter(!is.na(CTmin)) %>%
		data.frame()
	dat_OU3_CTmin <- dat %>%
		select(species_phylo, elevation, CTmin, CTmin_SE) %>%
		dplyr::filter(!is.na(CTmin)) %>%
		data.frame()
	# Only analyzing OU3 for CTmin, so simply make OU2 for L80
	dat_OU2_L80 <- dat %>%
		select(species_phylo, region, L80) %>%
		dplyr::filter(!is.na(L80)) %>%
		data.frame()

### (b) Set up phylogenies
	trees_ou_CTmin <- lapply(trees_ou, drop.tip, tip = 'Tlalocohyla_smithii')
	trees_ou_L80 <- trees_ou
	
  
### (c) Compare highest-ranked models for CTmin (BM1 vs. OU2 vs. OU3) ----
	nsims <- 1000
	OUwie_args <- list(scaleHeight = T, root.station = TRUE, algorithm = 'invert', quiet = T, warn = F)
	BMvOU2_CTmin <- paraboot_OUwie(
						trees = list(trees_ou_CTmin[['OU2']], trees_ou_CTmin[['OU2']]), 
						dat_list = list(dat_OU2_CTmin, dat_OU2_CTmin),
						mods = c('BM1', 'OUM'), 
						nsim = nsims, 
						seed = 1978,
						OUwie_args = OUwie_args, 
						ncores = NULL)
	ou2v3_CTmin <- paraboot_OUwie(
						trees = list(trees_ou_CTmin[['OU2']], trees_ou_CTmin[['OU3']]), 
						dat_list = list(dat_OU2_CTmin, dat_OU3_CTmin),
						mods = c('OUM', 'OUM'), 
						nsim = nsims, 
						seed = 1978,
						OUwie_args = OUwie_args, 
						ncores = NULL)
	BMvOU2_CTmin$stat_props # Type I error = 0.065; power = 0.992
	ou2v3_CTmin$stat_props # Type I error = 0.268; power = 0.598
	results_CTmin <- bind_rows(
		data.frame(comparison = "BMvOU2_CTmin", 
			null = BMvOU2_CTmin$LR_sim[,1], 
			test = BMvOU2_CTmin$LR_sim[,2], 
			lr = BMvOU2_CTmin$LR_emp, 
			sig = BMvOU2_CTmin$LR_sig
			),
		data.frame(comparison = "ou2v3_CTmin", 
			null = ou2v3_CTmin$LR_sim[,1], 
			test = ou2v3_CTmin$LR_sim[,2], 
			lr = ou2v3_CTmin$LR_emp, 
			sig = ou2v3_CTmin$LR_sig
		)) %>%
		pivot_longer(c(null, test), names_to = 'test_type', values_to = 'LR_ratio')
  
### (c) Compare highest-ranked models for L80 (BM1 vs. OU2 vs. OU3) ----
	nsims <- 1000
	OUwie_args <- list(scaleHeight = T, root.station = TRUE, algorithm = 'invert', quiet = T, warn = F)
	BMvOU1_L80 <- paraboot_OUwie(
						trees = list(trees_ou_L80[['OU2']], trees_ou_L80[['OU2']]), 
						dat_list = list(dat_OU2_L80, dat_OU2_L80),
						mods = c('BM1', 'OU1'), 
						nsim = nsims, 
						seed = 1978,
						OUwie_args = OUwie_args, 
						ncores = NULL)
	BMvOU2_L80 <- paraboot_OUwie(
						trees = list(trees_ou_L80[['OU2']], trees_ou_L80[['OU2']]), 
						dat_list = list(dat_OU2_L80, dat_OU2_L80),
						mods = c('BM1', 'OUM'), 
						nsim = nsims, 
						seed = 1977,
						OUwie_args = OUwie_args, 
						ncores = NULL)
	BMvOU1_L80$stat_props # Type I error = 0.030; power = 0.129
	BMvOU2_L80$stat_props # Type I error = 0.079; power = 0.250
	results_L80 <- bind_rows(
		data.frame(comparison = "BMvOU1_L80", 
					null = BMvOU1_L80$LR_sim[,1], 
					test = BMvOU1_L80$LR_sim[,2], 
					lr = BMvOU1_L80$LR_emp, 
					sig = qchisq(0.95, BMvOU1_L80$dof)
					),
		data.frame(comparison = "BMvOU2_L80", 
					null = BMvOU2_L80$LR_sim[,1], 
					test = BMvOU2_L80$LR_sim[,2], 
					lr = BMvOU2_L80$LR_emp, 
					sig = qchisq(0.95, BMvOU2_L80$dof))
		) %>%
		pivot_longer(c(null, test), names_to = 'test_type', values_to = 'LR_ratio') 

### (d) Link the two sets of results in order to plot together
	all_results <- rbind(results_CTmin, results_L80)
	all_results %>%
		separate(comparison, into = c('mods', 'var'), sep = '_') %>%
		unite(col = 'comparison', var, mods, sep = '_') %>% # Reversing the order
		mutate(comparison = factor(comparison)) %>%
		mutate(comparison = fct_recode(comparison,
							'(A) CTmin: BM vs. OU2' = 'CTmin_BMvOU2',
							'(B) CTmin: OU2 vs. OU3' = 'CTmin_ou2v3',
							'(C) L80: BM vs. OU1' = 'L80_BMvOU1',
							'(D) L80: BM vs. OU2' = 'L80_BMvOU2'
		)) %>%
		ggplot(aes(LR_ratio, fill = test_type)) + # Plot
			geom_density(alpha = 0.75) +
			geom_vline(aes(xintercept = lr), color = 'black', linetype = 'solid', size = 1) +
			geom_vline(aes(xintercept = sig), color = 'black', linetype = '31', size = 1) +
			scale_fill_manual(values = c('grey60', 'black')) + 
			scale_y_continuous(name = 'Density', expand = c(0,0)) +
			scale_x_continuous(name = 'Likelihood ratio', expand = c(0,0)) +
			theme_classic() +
			theme(
			  axis.title = element_text(size = 12),
			  axis.text = element_text(size = 9, color = 'black'),
			  strip.text = element_text(size = 12, hjust = 0.5, vjust = 2, face = 'bold'),
			  strip.background = element_blank(),
			  panel.spacing = unit(13, 'points')
			) +
			facet_wrap(~comparison, nrow = 2, scales = 'free')

# (e) Save (plotting parameters optimized for PDF, not RStudio plotting device)
	dims <- c(7.5, 6) # Width by height, in inches
	ggsave(paste0(fig_wd,'FigS1_paraboot_',Sys.Date(),'.pdf'), width = dims[1], height = dims[2], units = 'in')