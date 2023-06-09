
## ----setup, include=FALSE ------
### These lines are only necessary for recreating the R Markdown tutorial. 
### Do not run them if simply implementing the code. 
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(cache = TRUE) # Change to 'FALSE' when debugging
knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 75), tidy = TRUE) 


## ----load_data, message = F ------
library(OUwie) # Our principal package for OU analysis. Version 2.6. 
library(ouch) # An alternative package for OU analysis. Version 2.17.
library(phytools) # Used here for stochastic mapping. Version 0.7-70.
library(tidyverse) # Suite of packages for data manipulation. Key included packages used 
# in this tutorial are dplyr, tidyr, ggplot2, and purrr. Version 1.3.0.
source('AppendixS6_tutorial_functions.R') # Custom functions we wrote for this tutorial.
load('AppendixS7_OU_tutorial_data.RData') # R workspace with objects 'dat' and 'tree' 


## ----print_data ------
dat # Print data frame


## ----ancstates, dependson = 'load_data' ------
# Assign the values tropical-temperate to internal nodes
temptrop <- c(rep('tropical',5), rep('temperate',2), rep('tropical',2),rep('temperate',2))
tree_ou2 <- tree
tree_ou2$node.label <- temptrop	

# Assign the values highland-lowland-temperate to internal nodes
highlow <- c('lowland',rep('highland',3),'lowland',rep('temperate',2),
             rep('lowland',2),rep('temperate',2))
tree_ou3 <- tree
tree_ou3$node.label <- highlow


## ----ancstates_plot, echo = F, dependson = 'ancstates' ------
# OU2: temperate vs. tropical
par(mfcol=c(1,2), oma=c(1,1,1,0.5))
plot(tree, cex = 0.5, no.margin = T, label.offset = 1)
nodelabels(tree_ou2$node.label, cex = 0.5, adj = c(0, 0.5))

# OU3: temperate vs. tropical highland vs. tropical lowland
plot(tree, cex = 0.5, label.offset = 1)
nodelabels(tree_ou3$node.label, cex = 0.5, adj = c(0, 0.5))


## ----prune_trees, dependson = 'ancstates' ------
tree_ou2 <- drop.tip(tree_ou2, 'Tlalocohyla_smithii')
tree_ou3 <- drop.tip(tree_ou3, 'Tlalocohyla_smithii')


## ----data_setup, dependson = 'load_data' ------
dat_OU2 <- dat %>%
  dplyr::filter(!is.na(CTmin)) %>%
  select(-elevation) %>%
  data.frame()
dat_OU3 <- dat %>%
  dplyr::filter(!is.na(CTmin)) %>%
  select(-region) %>%
  data.frame()


## ----OUwie_analysis, message = F, warning = F, results = 'hide', dependson = c('data_setup', 'prune_trees', 'ancstates') ------
# (a) Brownian motion
  brown_OUwie <- OUwie(phy = tree_ou2, data = dat_OU2, model = 'BM1', mserr = 'known',
                       scaleHeight = 1, root.station = F, algorithm = 'invert', diagn = T)
# (b) Single-optimum OU
  OU1_OUwie <- OUwie(phy = tree_ou2, data = dat_OU2, model = 'OU1', mserr = 'known',
                     scaleHeight = 1, root.station = T, algorithm = 'invert', diagn = T)
# (c) Temperate-tropical OU2 
  OU2_OUwie <- OUwie(phy = tree_ou2, data = dat_OU2, model = 'OUM', mserr = 'known',
                     scaleHeight = 1, root.station = T, algorithm = 'invert', diagn = T)
# (d) Highland tropical, lowland tropical, and temperate OU3    
  OU3_OUwie <- OUwie(phy = tree_ou3, data = dat_OU3, model = 'OUM', mserr = 'known',
                     scaleHeight = 1, root.station = T, algorithm = 'invert', diagn = T)


## ----OUwie_stats, dependson = 'OUwie_analysis' ------
# Put output objects into a list to more easily summarize results
  all_out <- list(brown_OUwie, OU1_OUwie, OU2_OUwie, OU3_OUwie)
  names(all_out) <- c('BM', 'OU1', 'OU2', 'OU3')
  
# Summarize model fits
  CTmin_fits <- tibble(model = names(all_out),
                       k = map_dbl(all_out, ~.$param.count),
                       loglik = map_dbl(all_out, ~.$loglik),
                       AICc = map_dbl(all_out, ~.$AICc),
                       wgt = aicc_wgt_calc(AICc)
                )
  CTmin_fits # Print results table


## ----OUwie_params, dependson = 'OUwie_analysis' ------
# Summary of most model parameters, including MLEs and SEs of optima and MLEs of alpha and sigma squared
OU2_OUwie


# Standard errors of alpha and sigma squared, as calculated from the curvature of the likelihood surface
OU2_OUwie$solution.se 


## ----OUwie_eigen, dependson = 'OUwie_analysis' ------
OU2_OUwie$eigval


## ----ouch_setup, dependson = c('data_setup', 'prune_trees', 'ancstates') ------
# OU2 (temperate vs. tropical)
  dat_ouch_OU2 <- data.frame(dat_OU2, row.names = 'species_phylo')  
  dat_ouch_OU2 <- dat_ouch_OU2[,-3] # ouch does not use SEs
# Now we add the phylogeny and internal-state estimates
  dat_ouch_OU2 <- divers_to_ouch(tree_ou2, tree_ou2$node.label, dat_ouch_OU2)

# OU3 (temperate vs. tropical lowland vs. tropical highland)
  dat_ouch_OU3 <- data.frame(dat_OU3, row.names = 'species_phylo')  
  dat_ouch_OU3 <- dat_ouch_OU3[,-3]
  dat_ouch_OU3 <- divers_to_ouch(tree_ou3, tree_ou3$node.label, dat_ouch_OU3)


## ----ouch_plot_1, echo = F, dependson = 'ouch_setup' ------
# OU2: temperate vs. tropical
#par(mfcol=c(1,2), oma=c(1,1,1,0.5))
#par(mfcol = c(1,1))
ouch_tree <- with(dat_ouch_OU2, ouchtree(nodes, ancestors, times, labels)) 
plot(ouch_tree, regimes = dat_ouch_OU2['model'], cex = 0.75)


## ----ouch_plot_2, echo = F, dependson = 'ouch_setup' ------
# OU3: temperate vs. tropical highland vs. tropical lowland
ouch_tree <- with(dat_ouch_OU3, ouchtree(nodes, ancestors, times, labels))
plot(ouch_tree, regimes = dat_ouch_OU3['model'], cex = 0.75)


## ----ouch_analysis, dependson = c('ouch_setup', 'ouch_plot') ------
# Set ouchtree for all analyses
  ouch_tree <- with(dat_ouch_OU2, ouchtree(nodes, ancestors, times, labels)) 
# (a) Brownian motion
  brown_ouch <- brown(dat_ouch_OU2['CTmin'], ouch_tree)
# (b) Single-optimum OU 
  dat_ouch_OU2$global <- as.factor('global') # Set single regime for OU1
  OU1_ouch <- hansen(dat_ouch_OU2['CTmin'], ouch_tree, regimes = dat_ouch_OU2['global'],
                     sqrt.alpha = 1, sigma = 1)
# (c) Temperate-tropical OU2  
  OU2_ouch <- hansen(dat_ouch_OU2['CTmin'], ouch_tree, regimes = dat_ouch_OU2['model'],
                     sqrt.alpha = 1, sigma = 1)
# (d) Highland tropical, lowland tropical, and temperate OU3    
  OU3_ouch <- hansen(dat_ouch_OU3['CTmin'], ouch_tree, regimes = dat_ouch_OU3['model'],
                     sqrt.alpha = 1, sigma = 1)


## ----ouch_stats, dependson = 'ouch_analysis' ------
# Put output objects into a list to more easily summarize results
  all_out <- list(brown_ouch, OU1_ouch, OU2_ouch, OU3_ouch)
  names(all_out) <- c('BM', 'OU1', 'OU2', 'OU3')
# Summarize model fits
  ouch_sum <- map(all_out, summary)
  CTmin_fits <- tibble(model = names(all_out),
                       k = map_dbl(ouch_sum, ~.$dof),
                       loglik = map_dbl(ouch_sum, ~.$loglik),
                       AICc = map_dbl(ouch_sum, ~.$aic.c),
                       wgt = aicc_wgt_calc(AICc)
  )
  CTmin_fits # Print summary


## ----ouch_params, dependson = 'ouch_analysis' ------
coef(OU2_ouch)


## ----OU2_simmap_setup, message = F, warning = F, dependson = c('data_setup', 'prune_trees', 'ancstates') ------
# OU2: set up data
  # Make a vector of tip states and name them with their species' names
  states_ou2 <- setNames(dat_OU2[,2], dat_OU2[,1])
  # Put the data in the same order as the phylogeny
  states_ou2 <- states_ou2[tree_ou2$tip.label]
# OU2: Compare models
  AIC(ace(states_ou2, tree_ou2, type = 'discrete', model = 'ER'))
  AIC(ace(states_ou2, tree_ou2, type = 'discrete', model = 'ARD'))


## ----OU3_simmap_setup, dependson = c('data_setup', 'prune_trees', 'ancstates') ------    
# OU3: set up the data
  states_ou3 <- setNames(dat_OU3[,2], dat_OU3[,1]) 
  states_ou3 <- states_ou3[tree_ou3$tip.label] 
# OU3: compare models    
  ou3_er <- ace(states_ou3, tree_ou3, type = 'discrete', model = 'ER')
  ou3_sym <- ace(states_ou3, tree_ou3, type = 'discrete', model = 'SYM')
  ou3_ard <- ace(states_ou3, tree_ou3, type = 'discrete', model = 'ARD')
  map_dbl(list(ou3_er, ou3_sym, ou3_ard), AIC)


## ----simmap_analysis, message = F, results = 'hide', dependson = c('OU2_simmap_setup', 'OU3_simmap_setup') ------    
set.seed(1978) # To ensure your results look like ours
nsims <- 100 # Number of SIMMAP replicates
simmap_ou2 <- make.simmap(tree_ou2, states_ou2, model = 'ER', nsim = nsims)
simmap_ou3 <- make.simmap(tree_ou3, states_ou3, model = 'ER', nsim = nsims)


## ----simmap_plot, echo = FALSE, dependson = 'simmap_analysis' ------
# Plot results to examine a typical result
par(mfcol=c(1,2), mar=c(4,4,0.5,0.5), oma=c(1,2,1,1))
plotSimmap(simmap_ou2[[1]],fsize=0.75,ftype='i',colors=c(temperate='gray',
                                                         tropical='black'))
plotSimmap(simmap_ou3[[1]],fsize=0.75,ftype='i',colors=c(temperate='gray',
                                                         highland='black',
                                                         lowland='dodgerblue'))  


## ----OUwie_simmap_analysis, message = F, dependson = 'simmap_analysis' ------  
# Temperate-tropical OU2 (loop over SIMMAP trees)
OU2_OUwie_sm <- lapply(simmap_ou2, OUwie, data = dat_OU2, model = 'OUM', mserr = 'known',
                       simmap.tree = T, scaleHeight = T, root.station = T,
                       warn = F, quiet = T, algorithm = 'invert') 

# Temperate-highland-lowland OU3 (loop over SIMMAP trees)
OU3_OUwie_sm <- lapply(simmap_ou3, OUwie, data = dat_OU3, model = 'OUM', mserr = 'known',
                       simmap.tree = T, scaleHeight = T, root.station = T,
                       warn = F, quiet = T, algorithm = 'invert') 


## ----OUwie_simmap_summary, dependson = 'OUwie_simmap_analysis' ------  
OU2_fits <- tibble(run = seq(nsims), 
                   lnL = map_dbl(OU2_OUwie_sm, ~.$loglik), 
                   AICc = map_dbl(OU2_OUwie_sm, ~.$AICc))
OU3_fits <- tibble(run = seq(nsims),
                   lnL = map_dbl(OU3_OUwie_sm, ~.$loglik),
                   AICc = map_dbl(OU3_OUwie_sm, ~.$AICc))
out_table <- data.frame(matrix(nrow = 4, ncol = 4))  
colnames(out_table) <- c('model', 'AICc_mean', 'AICc_2.5', 'AIC_97.5')
out_table[,1] <- c('BM', 'OU1', 'OU2', 'OU3')
out_table[1:2, 'AICc_mean'] <- map_dbl(list(brown_OUwie, OU1_OUwie), ~.$AICc)
out_table[3, 'AICc_mean'] <- mean(OU2_fits$AICc)
out_table[3, c('AICc_2.5', 'AIC_97.5')] <- quantile(OU2_fits$AICc, probs = c(0.025,0.975))
out_table[4, 'AICc_mean'] <- mean(OU3_fits$AICc)
out_table[4, c('AICc_2.5', 'AIC_97.5')] <- quantile(OU3_fits$AICc, probs = c(0.025,0.975))
out_table # Print to see results


## ----OUwie_simmap_params, dependson = 'OUwie_simmap_analysis' ------
OU2_params <- data.frame(cbind(t(sapply(OU2_OUwie_sm, function(y) y$solution[1:2,1])),
                               t(sapply(OU2_OUwie_sm, function(y) y$theta[,1]))))
param_table <- data.frame(matrix(nrow = ncol(OU2_params), ncol = 4))
colnames(param_table) <- c('param', 'mean', 'low95CI', 'high95CI')
param_table[,1] <- c('alpha', 'sigma_sq', 'temp_theta', 'trop_theta')
param_table[,2] <- apply(OU2_params, 2, mean)
param_table[,3:4] <- t(apply(OU2_params, 2, quantile, probs = c(0.025, 0.975)))
param_table # Print to see results