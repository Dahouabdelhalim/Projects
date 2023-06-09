library(tidyverse)
library(magrittr)
library(vegan)
library(adespatial)
library(reshape2)
library(abind)
library(patchwork)

options(dplyr.summarise.inform = FALSE) # repress dplyr summarise() info

# Color Scheme

# Color scheme inspired by Neon Genesis Evangelion, Ep8 "Asuka Arrives in Japan", Dir. Anno, Hideaki

color_scheme = c("#5A72D1", "#F2CEA2", "#72B3CF", "#D97437", "#C4353A")
color_gradient = c("#5A72D1", "#838edb", "#a7ace6", "#c9caef")

# Caution: some (most) of the code is not "elegant" in any sense, but it works.
# What that means: https://www.reddit.com/r/ProgrammerHumor/comments/f4dbmj/saw_this_on_dank_memes/


########## MODEL SETUP ##########


# Initialize Metacommunity
# Input: size: the size of metacommunity (integer)
#        nspp: number of species (integer)
#        init_pop: initial population (vector)
#        random: whether to randomize population (boolean, default to F)
#        seed: optional random seed (integer, default to NULL)
# Output: a metacommunity with species as rows and patch as columns (matrix)
# This function generates initial metacommunities of Scenario 2 and 3

init_metacom = function(size, nspp, init_pop, each = F, random = F, seed = NULL) {
  
  # Set seed and initialize the metacommunity
  if (!is.null(seed)) {set.seed(seed)}
  metacom_matrix = matrix(0, nrow = nspp, ncol = size)
  
  # If the initial population is specified for each species
  if (length(init_pop) > 1) {metacom_matrix[, ] = init_pop}
  
  # If the initial population is the same for all species
  else {metacom_matrix[, ] = init_pop}
  
  # If each patch contains only one species, set all others to 0
  if (each) {
    p = sample(rep(1:nspp, size/nspp))
    for (i in 1:length(p)) {
      metacom_matrix[-p[i], i] = 0
    }
  } 

  # Randomization
  if (random) {
    return(apply(metacom_matrix, 1:2, function(x) ifelse(x > 0, rnorm(1, x, 15), 0)))
  }
  
  return(metacom_matrix)
}


# Calculate Competition Coefficients
# Input: B: the maximum of the sigmoidal scaling function (numeric)
#        scal: scaling constant of the sigmoidal scaling function (numeric)
#        deltp: relative emergence time (numeric)
#        type: type of priority effects (1 or 0)
# Output: a competition coefficient (numeric)
# If type == 1, priority effects are trait-dependent, and the competition coefficient is determined by deltp; 
# if type == 0, priority effects are frequency-dependent, and the competition coefficient is constant and equals to B/2

comp_coeff = function(B, scal, deltp, type) {
  
  # Reciprocal scaling function assuming early-arrival advantage
  if (type == 1) {return(B/(1+exp((deltp)/scal)))}
  
  # Constant scaling function
  if (type == 0) {return(B/2)}
}


# Get Mean Times of Emergence
# Input: nspp: number of species (integer)
#        max_time: the maximum (latest) emergence time (numeric)
# Output: a list of mean times of emergence (vector)

get_mean_times = function(nspp, max_time) {
  return(seq(0, max_time, max_time/(nspp-1)))
}


# Initialize Matrix of Delta p
# Input: nspp: number of species (integer), 
#        max_time: the maximum (latest) emergence time (numeric)
#        var: variation around mean times of emergence (numeric)
#        seed: optional random seed (integer, default to NULL)
# Output: a matrix of pairwise differences in times of emergence (delta p) (matrix)

init_deltp_matrix = function(nspp, max_time, var, seed = NULL) {
  
  # Set seed
  if (!is.null(seed)) {set.seed(seed)}
  
  # Initialize the matrix
  deltp_matrix = matrix(0, nrow = nspp, ncol = nspp)
  
  # If no difference in emergence times, just return all 0
  if (max_time == 0) {return(deltp_matrix)}
  
  # Get mean emergence times
  mean_times = get_mean_times(nspp, max_time)
  
  # Randomly draw actual emergence times from a normal distribution
  mean_times = sapply(mean_times, function(x) rnorm(1, x, sd = var))
  
  # Calculate differences in emergence times and assign to the matrix
  for (i in 1:length(mean_times)) {
    for (j in 1:length(mean_times)) {
      deltp_matrix[i, j] = mean_times[j]-mean_times[i]
    }
  }
  
  return(deltp_matrix)
}


# Initialize Interaction Matrix (Alpha Matrix)
# Input: nspp: number of species (integer)
#        deltp_matrix: matrix of delta p (matrix)
#        intra: intraspecific competition coefficients (vector)
#        B: the maximum of the sigmoidal scaling function (numeric)
#        scal: scaling constant of the sigmoidal scaling function (numeric)
#        type: type of priority effects (1 or 0)
# Output: a matrix of pairwise competition coefficients (matrix)

init_amat = function(nspp, deltp_matrix, intra, B, scal, type) {
  
  # Compute interspecific competition according to Eqn. 2
  amat = apply(deltp_matrix, 1:2, function(x) comp_coeff(B, scal, x, type))
  
  # Set intraspecific competition as the diagonal
  diag(amat) = intra
  
  return(amat)
}


# Density-dependent Population Growth
# Input: N: population of all species (vector)
#        amat: alpha matrix (matrix)
#        lambda: intrinsic growth rates of all species (vector)
# Output: new population of all species (vector)

competition = function(N, amat, lambda) {
  
  # Calculate population dynamics according to Eqn. 1
  new_N = lambda*N/(1+(amat %*% (N*lambda)))
  
  return(new_N)
}


# Dispersal
# Input: metacom_matrix: metacommunity matrix (matrix)
#        disp_rates: dispersal rates of all species (vector)
#        immi_patches: number of patches with dispersal (integer)
# Output: a list of two items: new metacommunity matrix after dispersal (matrix) and the proportion of dispersed population
#         of each species

dispersal = function(metacom_matrix, disp_rates, immi_patches) {
  
  # Extract size and number of species from the input metacommunity
  size = ncol(metacom_matrix)
  nspp = nrow(metacom_matrix)
  
  # If patches with dispersal is 0, just return the original metacommunity (no dispersal)
  if (immi_patches == 0) {return(list(metacom_matrix, rep(0, nspp)))}
  
  # Set up emigration matrix and immigration matrix
  E = matrix(0, nrow = nspp, ncol = immi_patches)
  I = matrix(0, nrow = nspp, ncol = size)

  # Emigration: draw emigrants of each species and stage from randomly selected patch
  source = sample(1:size, immi_patches)
  
  for (i in 1:length(source)) {
    N = metacom_matrix[, source[i]]
    
    # Draw emigrant from a binomial distribution; note that
    # this method does not allow dispersal if source population is less than 1
    e = rbinom(n = nspp, size = floor(N), prob = t(disp_rates))
    E[, i] = e
    
    metacom_matrix[, source[i]] = N-e
  }
  
  # Pool all emigrants
  pool = rowSums(E)
  
  # Warning message if the input patches of dispersal is larger than the size of the metacommunity
  if (immi_patches > size) stop ("Number of patches receiving immigrants exceeds metacommunity size")
  
  # Randomly select patches that accept immigrants
  destination = sample(1:size, immi_patches)
  
  # Immigration: randomly distribute emigrants to patches
  for (spp in 1:nspp) {
    I[spp, destination] = diff(c(0, sort(runif(immi_patches-1, 0, pool[spp])), pool[spp]))
  }
  
  return(list("metacom" = metacom_matrix+I, "dispersal" = rowSums(I)/rowSums(metacom_matrix+I)))
}


########## TEMPORAL DYNAMICS ##########


# Simulating Population Over Seasons
# Input: metacom_matrix: metacommunity matrix (matrix)
#        generation: number of generations in a season (integer)
#        season_length: length of each season in generations (integer)
#        n_seasons: number of seasons (integer), 
#        disp_rates: dispersal rates of all species (vector) 
#        immi_patches: number of patches with dispersal (integer)
#        disp_freq: dispersal frequency within a season (integer, minimum 1, maximum season_length)
#        disp_freq_seasons: dispersal frequency between seasons (integer, default to 1)
#        lambda: intrinsic growth rates of all species (vector)
#        intra: intraspecific competition coefficients (vector)
#        B: the maximum of the sigmoidal scaling function (numeric)
#        scal: scaling constant of the sigmoidal scaling function (numeric)
#        deltp: relative emergence time (numeric)
#        type: type of priority effects (1 or 0)
#        max_time: the maximum (latest) time of emergence (numeric)
#        var: variation around mean times of emergence (numeric)
#        burnin: a burn-in time with no dispersal for Scenario 1 (integer, default to 0)
#        seed: optional random seed (integer, default to NULL)
# Output: a list of three items: population dynamics over time (all_out; array with dimensions: number of species, time, number of patches),
#         all pairwise delta p over time (out_deltp; array with dimensions: number of species, number of species, time) and
#         all proportions of dispersed population of each species (out_disp; matrix)

seasons = function(metacom_matrix, season_length, n_seasons, # space and time
                   disp_rates, immi_patches, disp_freq, disp_freq_seasons = 1, # dispersal
                   lambda, intra, B, scal, type, # vital rates
                   max_time, var, n_disturb = 0, # delta p and disturbance
                   burnin = 0, seed = NULL) {
  
  # Warning message if the frequency of dispersal is smaller than 1
  if (disp_freq < 1) {stop("disp_freq should be an integer larger than 0")}
  
  # Extract size and number of species from the input metacommunity
  size = ncol(metacom_matrix)
  nspp = nrow(metacom_matrix)
  
  # Calculate total time of the simulation
  time = season_length*(n_seasons)+1

  # Set up time points for end of seasons
  seasons_tp = seq(2, time, season_length)

  # Set up time points for dispersal
  disp_tp = c()
  
  # If dispersal happens once over several seasons
  if (disp_freq_seasons >= 1) {
    index = seq(1, n_seasons+1, disp_freq_seasons)
    disp_tp = seasons_tp[index]
  }
  
  # If dispersal happens several times within one season
  else if (disp_freq > 1) {
    disp_tp = c()
    for (i in 1:(length(seasons_tp)-1)) {
      within_seasons = c(seasons_tp[i], sort(sample((seasons_tp[i]+1):(seasons_tp[i+1]-1), disp_freq-1)), seasons_tp[i+1])
      disp_tp = c(disp_tp, within_seasons)
    }
    disp_tp = unique(disp_tp)
  }
  
  # Recording population dynamics, delta p and dispersal during simulation
  out_metacom = array(NA, dim = c(nspp, time, size), dimnames = list(1:nspp, 1:time, 1:size))
  out_metacom[, 1, ] = metacom_matrix
  out_deltp = array(NA, dim = c(nspp, nspp, size, length(seasons_tp)), 
                    dimnames = list(1:nspp, 1:nspp, 1:size, 1:length(seasons_tp)))
  out_disp = matrix(NA, nrow = nspp, ncol = length(disp_tp))
  
  # Run the model
  for (t in 2:time) {
    for (patch in 1:size) {
      
      # Initialize competition matrix; if not at the beginning of the season, delta p matrix = 0
      deltp_matrix = init_deltp_matrix(nspp, max_time, var, seed)
      if (t %in% seasons_tp) {out_deltp[, , patch, match(t, seasons_tp)] = deltp_matrix}
      else {deltp_matrix = deltp_matrix*0}
      a_matrix = init_amat(nspp, deltp_matrix, intra, B, scal, type)

      # Run the Beverton-Holt model
      new_N = metacom_matrix[, patch]
      new_N = competition(new_N, a_matrix, lambda)
      
      # Update new population and record results
      metacom_matrix[, patch] = new_N
    }
    
    # Disturbance, if needed
    if (n_disturb != 0) {metacom_matrix[, sample(1:size, n_disturb)] = 0}
    
    # Dispersal, if needed
    if (t %in% disp_tp) {
      
      # For Scenario 1, turn off dispersal for the first 10 seasons
      # Change the following number to 0 for Scenario 2
      if (t > burnin) {
        disp = dispersal(metacom_matrix, disp_rates, immi_patches)
        metacom_matrix = disp[[1]]
        out_disp[, match(t, disp_tp)] = disp[[2]]
      } else {
        out_disp[, match(t, disp_tp)] = 0
      }
    }
    
    # If at the end of the season, multiply by survival rates
    if (t %in% seasons_tp) {
      end_season = rnorm(nspp, 0.5, 0.05)
      metacom_matrix = metacom_matrix*end_season
    }
    
    # Record the metacommunity
    out_metacom[, t, ] = metacom_matrix
  }
  
  return(list(out_metacom, out_deltp, out_disp))
}


########## SIMULATIONS ##########


# Run Multiple Replications of the Model
# Input: metacom_matrix: metacommunity matrix (matrix)
#        generation: number of generations in a season (integer)
#        season_length: length of each season in generations (integer)
#        n_seasons: number of seasons (integer), 
#        disp_rates: dispersal rates of all species (vector) 
#        immi_patches: number of patches with dispersal (integer)
#        disp_freq: dispersal frequency within a season (integer, minimum 1, maximum season_length)
#        disp_freq_seasons: dispersal frequency between seasons (integer, default to 1)
#        lambda: intrinsic growth rates of all species (vector)
#        intra: intraspecific competition coefficients (vector)
#        B: the maximum of the sigmoidal scaling function (numeric)
#        scal: scaling constant of the sigmoidal scaling function (numeric)
#        deltp: relative emergence time (numeric)
#        type: type of priority effects (1 or 0)
#        max_time: the maximum (latest) time of emergence (numeric)
#        var: variation around mean times of emergence (numeric)
#        rep: number of repetitions (integer)
#        burnin: a burn-in time with no dispersal for Scenario 1 (integer, default to 0)
#        seed: optional random seed (integer, default to NULL)
# Output: a list of three items: population dynamics over time (all_out; array with dimensions: number of species, time, number of patches
#         repetition), all pairwise delta p over time (out_deltp; array with dimensions: number of species, number of species, time, 
#         repetition) and all proportions of dispersed population of each species (out_disp; array with dimensions: number of species, time,
#         repetition)

multiple_reps = function(metacom_matrix, season_length, n_seasons, # space and time
                         disp_rates, immi_patches, disp_freq, disp_freq_seasons, # dispersal
                         lambda, intra, B, scal, type, # vital rates
                         max_time, var, rep, n_disturb = 0, # delta p and disturbance
                         burnin = 0, seed = NULL) {
  
  # Extract size and number of species from the input metacommunity
  nspp = nrow(metacom_matrix)
  size = ncol(metacom_matrix)
  
  # Calculate total time of the simulation
  time = season_length*n_seasons+1
  
  # Set up arrays for recording population dynamics, emergence times and dispersal
  all_dynamics = array(dim = c(nspp, season_length*n_seasons+1, size, rep))
  all_deltp = array(dim = c(nspp, nspp, size, n_seasons, rep))
  all_disp = array(dim = c(nspp, n_seasons*disp_freq+1, rep))

  for (r in 1:rep) {
    
    # Randomly draw a intraspecific competition coefficient from a normal distribution
    new_intra = rnorm(nspp, mean(intra), 0.01)

    # Run simulations
    dynamics = seasons(metacom_matrix, season_length, n_seasons, 
                       disp_rates, immi_patches, disp_freq, disp_freq_seasons, 
                       lambda, new_intra, B, scal, type, 
                       max_time, var, n_disturb, burnin, seed)
    
    # Record results
    all_dynamics[, , , r] = dynamics[[1]]
    all_deltp[, , , , r] = dynamics[[2]]
    all_disp[, , r] = dynamics[[3]]
  }
  
  return(list(all_dynamics, all_deltp, all_disp))
}


# Parse Arrays of Population Dynamics into Data Frame
# Input: dynamics: an array of population dynamics from multiple_rep (array)
# Output: a data frame with population dynamics of each species over time and all repetitions

met_df_multirep = function(dynamics) {
  
  out = plyr::adply(dynamics, c(1, 2, 3, 4))
  colnames(out) = c("species", "time", "patch", "rep", "population")
  out$time = as.numeric(out$time)-1
  
  return(out)
}


########## CALCULATING DIVERSITY METRICS ##########


# Calculating Beta Diversity
# Input: all_dynamics: array of population dynamics from multiple_reps (array)
# Output: a data frame with dissimilarity (Bray-Curtis) values of all repetitions

dissimilarity = function(all_dynamics) {
  
  # Get the end time of the simulation
  end_time = dim(all_dynamics)[2]
  
  # Set up the data frame for recording dissimilarities
  all_dissim = tibble()
  
  for (r in 1:dim(all_dynamics)[4]) {
    metacom_matrix = t(all_dynamics[, end_time, , r])
    metacom_matrix = metacom_matrix[rowSums(metacom_matrix)>0, ]
    
    # Calculate beta diversity as the mean pairwise dissimilarity between all patches
    all_dissim = rbind(all_dissim, tibble(dissim = mean(vegdist(metacom_matrix, "bray")), rep = r))
  }
  return(all_dissim)
}


# Calculating Species Richness (Simpson's Index)
# Input: all_dynamics: array of population dynamics from multiple_reps (array) 
#        regional: the range of calculation (default to F, calculating local species richness; 
#                                            if T, calculate regional species richness)
# Output: a data frame of local (alpha diversity) or regional (gamma diversity) species richness of all repetitions
# Note: species richness is calculated by Simpson's index (function diversity() in the package vegan)

alphadiv = function(all_dynamics, regional = F) {
  
  # Get the end metacommunity of all simulations
  end_time = dim(all_dynamics)[2]
  end_array = all_dynamics[, end_time, , ]
  
  # Set up the vector for recording richness
  rch = numeric()
  
  # Go through all repetitions
  for (i in 1:dim(end_array)[3]) {
    
    # Calculate regional richness (gamma diversity)
    if (regional) {rch = c(rch, diversity(colSums(t(end_array[, , i]))))} 
    
    # Calculate local richness (alpha diversity), then average across all patches
    else {rch = c(rch, mean(diversity(end_array[, , i], MARGIN = 2)))}
  }
  
  richness = data.frame(richness = rch, rep = 1:dim(end_array)[3])
  
  return(richness)
}


# Calculate Temporal Beta Diversity
# Input: all_dynamics: array of population dynamics from multiple_reps (array) 
#        step: length between each time steps from which values will be calculated (integer)
# Output: a data frame with species gain, loss and dissimilarity over time calculated by TBI() in the package adespatial
#         from dynamics of all repetitions

temporal_beta = function(all_dynamics, step) {
  
  # Set up the data frame for recording metrics
  bcd = tibble()
  
  # Time points for calculating metrics
  tp = seq(1, dim(all_dynamics)[2], step)
  
  # Go through all repetitions
  for (r in 1:(dim(all_dynamics)[4])) {
    
    # Go through all time points
    for (i in 1:(length(tp)-1)) {
      
      # Actual calculation
      b = TBI(t(all_dynamics[, tp[i], , r]), t(all_dynamics[, tp[i+1], , r]))
      avg_B = mean(b$BCD.mat[, 1])
      avg_C = mean(b$BCD.mat[, 2])
      avg_D = mean(b$BCD.mat[, 3])
      bcd = rbind(bcd, tibble(B = avg_B, C = avg_C, D = avg_D, time = i*10+1, rep = r))
    }
  }
  
  return(bcd)
}


# Calculate Diversity Metrics with Multiple Dispersal Parameters
# Input: metacom_matrix: metacommunity matrix (matrix)
#        generation: number of generations in a season (integer)
#        season_length: length of each season in generations (integer)
#        n_seasons: number of seasons (integer), 
#        disp_rates_lst: dispersal rates (vector) 
#        immi_patches: number of patches with dispersal (integer)
#        disp_freq: dispersal frequency within a season (integer, minimum 1, maximum season_length)
#        disp_freq_seasons: dispersal frequency between seasons (integer, default to 1)
#        lambda: intrinsic growth rates of all species (vector)
#        intra: intraspecific competition coefficients (vector)
#        B: the maximum of the sigmoidal scaling function (numeric)
#        scal: scaling constant of the sigmoidal scaling function (numeric)
#        deltp: relative emergence time (numeric)
#        type: type of priority effects (1 or 0)
#        max_time: the maximum (latest) time of emergence (numeric)
#        var: variation around mean times of emergence (numeric)
#        rep: number of repetitions (integer)
#        burnin: a burn-in time with no dispersal for Scenario 1 (integer, default to 0)
#        seed: optional random seed (integer, default to NULL)
# Output: a data frame with diversity metrics over different dispersal parameters and all repetitions

multiple_disp = function(metacom_matrix, season_length, n_seasons, # space and time
                         disp_rates_lst, immi_patches_lst, disp_freq, disp_freq_seasons = 1, # dispersal
                         lambda, intra, B, scal, type, # vital rates
                         max_time, var = 0, rep, n_disturb = 0, 
                         burnin = 0, seed = NULL) {
  
  # Set up the data frame for recording values
  all_out = tibble()
  
  # Go through all dispersal rates
  for (d in disp_rates_lst) {
    
    # Go through all patches with dispersal
    for (i in immi_patches_lst) {
      
      # Run simulations
      res = multiple_reps(metacom_matrix, season_length, n_seasons, 
                          disp_rates = rep(d, nspp), immi_patches = i, disp_freq, disp_freq_seasons,
                          lambda, intra, B, scal, type, 
                          max_time, var = var, rep = rep, n_disturb = n_disturb, burnin = burnin, seed = seed)[[1]]
      
      local = alphadiv(res)
      dis = dissimilarity(res)
      regional = alphadiv(res, regional = T)
      
      # Calculate biodiversity and record values
      all_out = rbind(all_out, tibble(rep = local$rep, 
                                      alpha = local$richness, 
                                      beta = dis$dissim, 
                                      gamma = regional$richness, 
                                      disp_rate = d, immi_patches = i))
    }
  }
  
  return(all_out)
}

