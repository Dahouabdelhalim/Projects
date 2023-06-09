### R custom functions written by Daniel Moen for the paper: "Phylogenetic analysis 
#	of adaptation in comparative physiology and biomechanics: overview and a case study of 
#	thermal physiology in treefrogs" by Moen, Cabrera-Guzmán, Caviedes-Solis, González-Bernal, and Hanna


### OU simulation function for Fig. 1 ----
# (1) Input:
#   (a) 'alpha' = Alpha value before and after change (two-element numeric vector)
#   (b) 'sigma_sq' = Sigma-squared value before and after change (two-element numeric vector)
#   (c) 'theta' = Theta value before and after change (two-element numeric vector; default = c(0,0))
#   (d) 'nsim' = Number of simulation replicates (default = 1)
#   (e) 'time' = Total length of simulation (default = 1)
#   (f) 'switch_cond' = Proportion of total time at which simulation conditions change (default = 0.5)
#   (g) 'num_steps' = Number of time steps (break down with seq(), total time, and switch point). Default = 100.
# (2) Output:
#   (a) A matrix with time step in the first column and 'nsim' additional columns with output

### Function
  OU_singlebranch_sim <- function(alpha = NULL, sigma_sq = NULL, theta = c(0,0), nsim = 1, time = 1, switch_cond = 0.5, num_steps = 100){
    switch_cond <- switch_cond * time # Put the switching point on an absolute scale
    # Start with a single simulation run  
      onesim <- function(alpha, theta, sigma_sq, time, switch_cond, num_steps){
        times <- seq(0, time, length.out = num_steps)
        BM_changes <- c(rnorm(sum(times < switch_cond), 0, sigma_sq[1]), rnorm(sum(times >= switch_cond), 0, sigma_sq[2]))
        alphas <- ifelse(times < switch_cond, alpha[1], alpha[2])
        thetas <- ifelse(times < switch_cond, theta[1], theta[2])
        out <- numeric(length(times))
        for(i in 2:length(times)){
          out[i] <- out[i-1] + (BM_changes[i] - alphas[i] * (out[i-1] - thetas[i])) * (time / num_steps) # Scale down by time-step size, otherwise simulation will be biased by number of time steps
        }
        return(out)  
      }
   # Now run all simulations       
    all_out <- sapply(1:nsim, function(y) onesim(alpha, theta, sigma_sq, time, switch_cond, num_steps))
    if(nsim > 1){
      times <- seq(0, time, length.out = num_steps)
      alphas <- ifelse(times < switch_cond, alpha[1], alpha[2])
      thetas <- ifelse(times < switch_cond, theta[1], theta[2])
      all_out <- as.data.frame(cbind(times, alphas, thetas, all_out))
      colnames(all_out) <- c('times', 'alphas', 'thetas', paste0('sim',seq(nsim)))
    } else {
      all_out <- as.data.frame(cbind(times, all_out))
      colnames(all_out) <- c('times', paste0('sim',seq(nsim)))
    }
    return(all_out)
  }


### Function for calculating AICc weight
    aicc_wgt_calc <- function(y){ # Input 'y' is a list of AICc scores
      dif <- y - min(y) # Calculate differences between models and best model (i.e. lowest AIC score)
      d <- exp(-0.5 * (dif)) # Exponentiate to produce weights
      wg <- d / sum(d) # Normalize weights
      return(wg) # Return weights
    }

### Function for conducting and summarizing OUwie analyses, with root.station = F for BM analyses when specified TRUE for OU analyses
# This function runs OUwie across all input response variables, filtering by inclusion or not of Tlalocohyla smithii, then pruning the tree
# INPUT:
#	(1) 'dat' is a tibble with species name in the column 'species.phylo', columns 'region', 'elevation', and one fitting the variable in 'resp_var', and possibly an SE column
#	(2) 'phy' is a two-element 'multiPhylo' object with the first tree being a mapping of region, and the second a mapping of elevation
#	(3) 'resp_var' is the response trait name
#	(4) 'SE' is a Boolean for whether you are analyzing with standard errors
# 	(5) 'wd_out' is a directory address (full) if you want to specify where the models and summary matrices are output. Needs a terminal '/' to work properly.
# 	(6) 'get.root.theta' is a Boolean for whether you want to estimate a separate root optimum state
# 	(7) 'OUwie_args' is a list of optional arguments for OUwie
# OUTPUT:
#	(1) A list of four models fits (BM, OU1, OU2, OU3)
### Function
  OUwie_hylids <- function(dat, phy, resp_var, SE = FALSE, wd_out = NULL, get.root.theta = FALSE, OUwie_args = NULL){
    library(OUwie)
    library(tidyverse)
    # Prep data
    dat1 <- dat %>%
      dplyr::filter(!is.na(get(resp_var))) %>%
      select(species.phylo, region, all_of(resp_var), CTmin_SE) %>% 
      data.frame()
    dat2 <- dat %>%
      filter(!is.na(get(resp_var))) %>%
      select(species.phylo, elevation, all_of(resp_var), CTmin_SE) %>% 
      data.frame()
    # Remove CTmin SE if not using it
    ifelse(SE, sb <- c(1:3,4), sb <- 1:3)
    dat1 <- dat1[,sb]	
    dat2 <- dat2[,sb]
    if(nrow(dat1) != length(phy[[1]]$tip.label)){
      phy[[1]] <- keep.tip(phy[[1]], dat1$species.phylo)
      phy[[2]] <- keep.tip(phy[[2]], dat1$species.phylo)
    }
    mse <- ifelse(SE, 'known', 'none')
    data_list <- list(	BM = list(phy[[1]], dat1, 'BM1'),
                       OU1 = list(phy[[1]],dat1, 'OU1'),
                       OU2 = list(phy[[1]], dat1, 'OUM'),
                       OU3 = list(phy[[2]], dat2, 'OUM')
    )
    # Model fits here, starting by pruning the trees if necessary. Note that some arguments are default, but based 
    # on changes in OUwie since May 2020 (when we started these analyses), I want to make them explicit
    #	mod_fit <- map(data_list, ~OUwie(.[[1]], .[[2]], .[[3]], scaleHeight=T, mserr = mse, root.station = T,
    #	                                 shift.point = 0.5, diagn = T, warn = F, quiet = T, algorithm = 'invert', get.root.theta = get.root.theta))
    mod_fit <- list()
    BMargs <- OUwie_args
    BMargs$root.station <- FALSE
    mod_fit[[1]] <- do.call(OUwie, c(list(data_list[[1]][[1]], data_list[[1]][[2]], model = 'BM1', mserr = mse, get.root.theta = get.root.theta, diagn = T), BMargs))
    mod_fit[2:4] <- map(data_list[2:4], ~do.call(OUwie, c(list(phy = .[[1]], data = .[[2]], model = .[[3]], mserr = mse, get.root.theta = get.root.theta, diagn = T), OUwie_args)))
    # Summarize model comparison
    aicc_wgt_calc <- function(y){ # Input 'y' is a list of AIC scores
      dif <- y - min(y) # Calculate differences between models and best model (i.e. lowest AIC score)
      d <- exp(-0.5 * (dif)) # Exponentiate to produce weights
      wg <- d / sum(d) # Normalize weights
      out <- matrix(c(dif, wg), ncol = 2)
      return(out) # Return weights
    }
    mod_sum <- data.frame(matrix(ncol = 6, nrow = length(mod_fit)))
    names(mod_sum) <- c('model','k','lnL','AICc','dif','wgt')
    mod_sum[,1] <- names(data_list)
    mod_sum[,2] <- map_dbl(mod_fit, ~.$param.count)
    mod_sum[,3] <- round(map_dbl(mod_fit, ~.$loglik), 3)
    mod_sum[,4] <- round(map_dbl(mod_fit, ~.$AICc), 3)
    mod_sum[,5:6] <- round(aicc_wgt_calc(mod_sum[,4]), 3)
    # Extract parameters and standard errors from model fits
    mod_pars <- data.frame(matrix(ncol = 13, nrow = length(mod_fit)))
    params <- c('sigma','alpha','theta.temp','theta.high','theta.low')
    names(mod_pars) <- c('model', paste0(rep(params, each = 2), c('_MLE','_SE')), 'eigval1', 'eigval2')
    mod_pars[,1] <- c('BM','OU1','temptrop','elevation')
    mod_pars[,2] <- map_dbl(mod_fit, ~.$solution[2,1])
    mod_pars[,3] <- map_dbl(mod_fit, ~.$solution.se[2,1])
    mod_pars[,4] <- map_dbl(mod_fit, ~.$solution[1,1])
    mod_pars[,5] <- map_dbl(mod_fit, ~.$solution.se[1,1])			
    if(get.root.theta){
      mod_pars[1,6:7] <- mod_fit[[1]]$theta[1,]
      mod_pars[2,6:7] <- mod_fit[[2]]$theta[2,]
      #mod_pars[1:2,6] <- map_dbl(mod_fit[1:2], ~.$theta[1,])
      #mod_pars[1:2,7] <- map_dbl(mod_fit[1:2], ~.$theta[1,2])
      mod_pars[3,c(6,8)] <- mod_fit[[3]]$theta[-1,1]
      mod_pars[3,c(7,9)] <- mod_fit[[3]]$theta[-1,2]
      mod_pars[4,c(6,8,10)] <- mod_fit[[4]]$theta[c(4,2,3),1]
      mod_pars[4,c(7,9,11)] <- mod_fit[[4]]$theta[c(4,2,3),2]
    } else {
      mod_pars[1:2,6] <- map_dbl(mod_fit[1:2], ~.$theta[1,1])
      mod_pars[1:2,7] <- map_dbl(mod_fit[1:2], ~.$theta[1,2])
      mod_pars[3,c(6,8)] <- mod_fit[[3]]$theta[,1]
      mod_pars[3,c(7,9)] <- mod_fit[[3]]$theta[,2]
      mod_pars[4,c(6,8,10)] <- mod_fit[[4]]$theta[c(3,1,2),1]
      mod_pars[4,c(7,9,11)] <- mod_fit[[4]]$theta[c(3,1,2),2]
    }
    names(mod_fit) <- seq(mod_fit) # Necessary for next function to work
    mod_pars[,12:13] <- t(map_df(mod_fit, ~.$eigval))
    mod_pars[,-1] <- round(mod_pars[,-1], 3)			
    # Save
    if(!is.null(wd_out)){
      save(mod_fit, file = paste0(wd_out, resp_var,'_',ifelse(SE,'withSE_','withoutSE_'),'model_fits_',Sys.Date(),'.Rdata'))
      write_csv(mod_sum, path = paste0(wd_out, resp_var,'_',ifelse(SE,'withSE_','withoutSE_'),'model_comparison_',Sys.Date(),'.csv'))
      write_csv(mod_pars, path = paste0(wd_out, resp_var,'_',ifelse(SE,'withSE_','withoutSE_'),'model_parameters_',Sys.Date(),'.csv'))
    }
    # Return a list of model fits for later summarization
    return(list(fitted_obs = mod_fit, mod_summary = mod_sum, mod_parameters = mod_pars))
}


### Parametric bootstrapping function for OUwie, for model testing (i.e. OUwie.boot() is for calculating 95% CIs of parameters only)
# INPUT:
# (1) 'trees' is a list of two trees that you will use for the simulations (tree[[1]] has node labels for the simpler model, tree[[2]] for the more complex model)
# (2) 'dat_list' is a list of two data matrices, with columns 1 and 3 identical (same taxa and continuous data) and column 2 different, corresponding to the two regime models. 
#   You cannot have a fourth column for SEs, since OUwie.sim() does not support them. Also, the two data frames can be identical if the simpler model is either Brownian motion or OU1.
# (3) 'mods' is a two-element vector of the models you want to simulate and test. These should be identical to those in OUwie. 
# (4) 'criterion' is what you want for output: likelihood-ratio test, AICc comparisons, or both (default)
# (5) 'nsim' is the number of simulation replicates you want.
# (6) 'OUwie_args' is a list of any additional arguments to pass on to 'OUwie()'. Cannot include 'mserr', because OUwie simulations do not incorporate standard errors. 
# (7) 'ncores' is the number of cores you want to parallelize the analysis over. The default NULL will use the package 'parallel' to detect the maximum number on your computer. 
# (8) 'seed' if you want to set the random-number generator for producing repeatable results
# OUTPUT:
# A list with degrees of freedom of the LR-test, the empirical LR-test, a 'null' matrix of LR tests and AICc comparisons when the data were simulated
# under the first model (mods[[1]]), and a 'test' matrix of LR tests and AICc comparisons when the data were simulated under the second model (mods[[2]])
paraboot_OUwie <- function(trees, dat_list, mods, criterion = c('LR-test', 'AICc'), nsim = 100, OUwie_args = NULL, ncores = NULL, seed = NULL){
  require(OUwie)
  if(is.null(ncores)){
  	require(parallel)
  	ncores <- parallel::detectCores()
  }
  if(ncores > 1) require(parallel)
  if(!is.null(seed)) set.seed(seed)
  if(!any(mods %in% c('BM1','OU1','OUM'))) stop('The only OUwie mods tested for this function are BM1, OU1, and OUM')
  # if(length(criterion) == 2) stop('You need to choose a single criterion')
  # Run MLE fits for data
  mod_fit <- list()
  if(mods[[1]] == 'BM1'){
    BM_args <- OUwie_args
    BM_args$root.station <- FALSE # Returns incorrect results for BM1 if root.station == TRUE
    mod_fit[[1]] <- do.call(OUwie, c(list(phy = trees[[1]], data = dat_list[[1]], model = mods[[1]], mserr = 'none'), BM_args))
  } else {
    mod_fit[[1]] <- do.call(OUwie, c(list(phy = trees[[1]], data = dat_list[[1]], model = mods[[1]], mserr = 'none'), OUwie_args))  
  }
  mod_fit[[2]] <- do.call(OUwie, c(list(phy = trees[[2]], data = dat_list[[2]], model = mods[[2]], mserr = 'none'), OUwie_args))
  # Simulation
  sim_dat <- list()
  for (i in seq(mods)){
    tmp <- mod_fit[[i]]
    if(mods[[i]] %in% c('BM1','OU1')){
    # Run BM simulations
      if(mods[[i]] == 'BM1'){
        # BM simulation
        sim_dat[[i]] <- lapply(1:nsim, function(y) OUwie.sim(phy = trees[[1]], data = dat_list[[1]], simmap.tree = tmp$simmap.tree, 
                                                             root.age = NULL, scaleHeight = tmp$scaleHeight, shift.point = tmp$shift.point,
                                                             alpha = c(0,0), sigma.sq = rep(tmp$solution[2,1], 2), 
                                                             theta0 = tmp$theta[1,1], theta = c(0,0)))
      } else {
        # OU1 simulation
        sim_dat[[i]] <- lapply(1:nsim, function(y) OUwie.sim(phy = trees[[1]], data = dat_list[[1]], simmap.tree = tmp$simmap.tree, 
                                                             root.age = NULL, scaleHeight = tmp$scaleHeight, shift.point = tmp$shift.point,
                                                             alpha = rep(tmp$solution[1,1], 2), sigma.sq = rep(tmp$solution[2,1], 2), 
                                                             theta0 = tmp$theta[1,1], theta = rep(tmp$theta[1,1], 2)))
      }
    } else {
      # Any other type of simulation
      sim_dat[[i]] <- lapply(1:nsim, function(y) OUwie.sim(fitted.object = tmp, simmap.tree = tmp$simmap.tree, scaleHeight = tmp$scaleHeight, 
                                                           shift.point = tmp$shift.point))
    }
    # Since 'OUwie.sim()' puts a numeric in place of the named regimes in the original data, we need to extract the first two columns of the original data
    sim_dat[[i]] <- lapply(sim_dat[[i]], function(y) data.frame(dat_list[[1]][,1:2], y$X))
  }
  simA <- sim_dat[[1]] # Makes code below a bit simpler
  simB <- sim_dat[[2]] # Makes code below a bit simpler
  # Analyze simulation replicates
  if(ncores == 1){
    # Non-parallel analyses  
    if(mods[[1]] == 'BM1'){
      AA <- lapply(simA, function(y) do.call(OUwie, c(list(phy = trees[[1]], data = y, model = mods[[1]], mserr = 'none'), BM_args)))  
      BA <- lapply(simB, function(y) do.call(OUwie, c(list(phy = trees[[1]], data = y, model = mods[[1]], mserr = 'none'), BM_args)))
    } else {
      AA <- lapply(simA, function(y) do.call(OUwie, c(list(phy = trees[[1]], data = y, model = mods[[1]], mserr = 'none'), OUwie_args)))
      BA <- lapply(simB, function(y) do.call(OUwie, c(list(phy = trees[[1]], data = y, model = mods[[1]], mserr = 'none'), OUwie_args)))
    }
    AB <- lapply(simA, function(y) do.call(OUwie, c(list(phy = trees[[2]], data = y, model = mods[[2]], mserr = 'none'), OUwie_args)))
    BB <- lapply(simB, function(y) do.call(OUwie, c(list(phy = trees[[2]], data = y, model = mods[[2]], mserr = 'none'), OUwie_args)))
    
    # Four searches, each time extracting AICc and lnL: AA, AB, BA, BB
    # Then get DOF difference between modA and modB for LR-test
  } else {
    # Parallel analyses
    if(mods[[1]] == 'BM1'){
      AA <- mclapply(simA, function(y) do.call(OUwie, c(list(phy = trees[[1]], data = y, model = mods[[1]], mserr = 'none'), BM_args)), mc.cores = ncores)  
      BA <- mclapply(simB, function(y) do.call(OUwie, c(list(phy = trees[[1]], data = y, model = mods[[1]], mserr = 'none'), BM_args)), mc.cores = ncores)
    } else {
      AA <- mclapply(simA, function(y) do.call(OUwie, c(list(phy = trees[[1]], data = y, model = mods[[1]], mserr = 'none'), OUwie_args)), mc.cores = ncores)
      BA <- mclapply(simB, function(y) do.call(OUwie, c(list(phy = trees[[1]], data = y, model = mods[[1]], mserr = 'none'), OUwie_args)), mc.cores = ncores)
    }
    AB <- mclapply(simA, function(y) do.call(OUwie, c(list(phy = trees[[2]], data = y, model = mods[[2]], mserr = 'none'), OUwie_args)), mc.cores = ncores)
    BB <- mclapply(simB, function(y) do.call(OUwie, c(list(phy = trees[[2]], data = y, model = mods[[2]], mserr = 'none'), OUwie_args)), mc.cores = ncores)
  }
  # Summarize results  
  dof <- mod_fit[[2]]$param.count - mod_fit[[1]]$param.count
  emp_lr <- 2 * (mod_fit[[2]]$loglik - mod_fit[[1]]$loglik)
  sig_lr <- qchisq(0.95, dof)
  if('LR-test' %in% criterion){
    ll_AA <- unlist(lapply(AA, function(y) y$loglik))
    ll_AB <- unlist(lapply(AB, function(y) y$loglik))
    ll_BA <- unlist(lapply(BA, function(y) y$loglik))
    ll_BB <- unlist(lapply(BB, function(y) y$loglik))
    lr_test_sim <- data.frame(null = 2 * (ll_AB - ll_AA), test = 2 * (ll_BB - ll_BA))
    stats <- setNames(c(mean(lr_test_sim$null > sig_lr), mean(lr_test_sim$test > sig_lr)), c('TypeI_error', 'Power'))
  } else {
    lr_test_sim <- NULL
    stats <- NULL
  }
  if('AICc' %in% criterion){
    aicc_AA <- unlist(lapply(AA, function(y) y$AICc))
    aicc_AB <- unlist(lapply(AB, function(y) y$AICc))
    aicc_BA <- unlist(lapply(BA, function(y) y$AICc))
    aicc_BB <- unlist(lapply(BB, function(y) y$AICc))
    aicc_sim <- data.frame(null = aicc_AB - aicc_AA, test = aicc_BA - aicc_BB)
  } else {
    aicc_sim <- NULL
  }
  # Kick out:
  return(list(dof = dof, LR_emp = emp_lr, LR_sig = sig_lr, stat_props = stats, LR_sim = lr_test_sim, AICc_sim = aicc_sim,
              mod1_fit = mod_fit[[1]], mod2_fit = mod_fit[[2]]))
}