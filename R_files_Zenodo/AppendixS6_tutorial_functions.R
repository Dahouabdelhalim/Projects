### R code associated with the R Markdown tutorial. This file has the custom functions that are 
#	used in the tutorial.

### Calculate AICc weight
    aicc_wgt_calc <- function(y){ # Input 'y' is a list of AICc scores
      dif <- y - min(y) # Calculate differences between models and best model (i.e. lowest AIC score)
      d <- exp(-0.5 * (dif)) # Exponentiate to produce weights
      wg <- d / sum(d) # Normalize weights
      return(wg) # Return weights
    }

### Convert ape/diversitree format of ancestral states to OUCH formatted data matrix, which specifies regime states of internal nodes
### INPUT
	# 'tree' should be a standard phylogeny of class 'phylo'
	# 'node_states' should be a vector of internal states as given by the R package 'diversitree' (e.g. the 'joint' estimate). The follows the phylogeny format of R package 'ape'. 
	# 'dat' is a data frame of 1 discrete character (1st column) and as many continuous characters as you want, with row names corresponding to tip names on the phylogeny. 
    # Assumes all species in the phylogeny are in the data, and vice versa. 
	# NOTE: if you want to run a single-optimum OU model, just insert vectors of the appropriate length where all elements are 'global'
### OUTPUT
	# A 'ouch' class data matrix

### Function (will need to be called for each individual model)
divers_to_ouch <- function(tree, node_states, dat){
    require(ouch)
    require(geiger)
    if (geiger:::name.check(tree,dat) != "OK") stop('Something is weird with the tree or tip.states')
    dat <- dat[tree$tip.label,]
    colnames(dat)[1] <- 'model'
    dat[,1] <- as.factor(dat[,1]) # Otherwise next chunk of code on 'node_states' won't work
    # Check to see if the input is a data frame or vector; needs to be a vector
    if(!is.null(dim(node_states))) node_states <- node_states[,1]
    # Ensure that 'node_states' are numeric
    if(!is.numeric(node_states)){ 
      node_states <- as.factor(node_states)
      # The following code ensures that we do not skip any states due to them being missing in the internal nodes; simply 
      # calling 'as.numeric' would number nodes only as high as the number of levels within the internal nodes
      node_states <- unlist(lapply(levels(node_states),function(x) which(levels(dat[,1])==x)))[as.numeric(node_states)]
      ### Old code below to show how it works (not much clearer...)
      #state <- unlist(lapply(levels(node_states),function(x) which(levels(tip.states[,1])==x)))
      #tmp <- as.numeric(node_states)
      #node_states <- state[tmp]
    }
    tree$node.label <- node_states # Here is where the states need to be numeric, so as to not screw up later steps in building the ouch file
    ouch.tree <- ape2ouch(tree)
    otd <- as(ouch.tree,"data.frame")
    ### in these data, it so happens that the rownames correspond to node names
    ### we will exploit this correspondence in the 'merge' operation:
    dat$labels <- rownames(dat)
    otd <- merge(otd,dat,by="labels",all=TRUE)
    rownames(otd) <- otd$nodes
    ouch.tree <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
    ### Now put together regimes
    ### Multi-state model, using node labels as those states with the highest marginal likelihood
    nd <- tree$Nnode # Number of nodes
    int.lab <- as.factor(otd[1:nd,]$labels) # Seemingly numeric internal node labels, these are actually factors
    mi <- max(as.numeric(int.lab)) # The number of different states at internal nodes
    states <- levels(dat[,1])[as.numeric(levels(int.lab))[1:mi]] # Gives only the states that apply to internal nodes
    otd$model[1:nd] <- states[as.numeric(int.lab)]
    return(otd) # Returns the 'ready to go' matrix for OUCH analyses
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