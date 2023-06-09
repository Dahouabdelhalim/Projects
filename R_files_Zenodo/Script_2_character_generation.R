################################################################################ 

# SCRIPT 2: GENERATING DISCREET DATASETS - STANDARD MARKOV, HIDDEN RATES,
# AND AMPLIFIED HIDDEN RATES.

# AUTHOR: Jeremy D. Wilson

# ARTICLE: Wilson et al. 2022 - Chronogram or phylogram for ancestral state
# estimation? Model-fit statistics indicate the branch lengths underlying a
# binary characterâ€™s evolution


################################################################################

# LOAD PACKAGES & IMPORT DATA---------------------------------------------------


require(tidyverse)
require(castor) 
load("trees.Rdata")


# SIMULATING BINARY CHARACTERS WITH MARKOV 'ARD' MODEL--------------------------


data$character_branch_lengths <- sample(c("chronogram", "phylogram"),
                                        nrow(data),
                                        replace = T
)


# Generate transition matrices (Q).

data$Q <- replicate(nrow(data),
                    get_random_mk_transition_matrix(Nstates = 2,
                                                    rate_model = "ARD"),
                    simplify = F
)

data$Q_rate_1 <- map_dbl(data$Q, function(x){x[1,2]})
data$Q_rate_2 <- map_dbl(data$Q, function(x){x[2,1]})


# Function to generate character.

sim.binary.markov <- function(Chrono, Phylo, BranchLengths, Q){
  states_present <- 1
  min_frequency <- 0
  if(BranchLengths == "chronogram"){
    while(states_present < 2 | min_frequency < 0.05){
      character <- simulate_mk_model(tree = Chrono, 
                                     Q = Q,
                                     root_probabilities = "stationary",
                                     Nsimulations = 1,
                                     drop_dims = F
      )
      states_present <- length(unique(character$tip_states))
      min_frequency <- min(table(character$tip_states)/
                             length(character$tip_states))
    }
  } else {
    while(states_present < 2 | min_frequency < 0.05){
      character <- simulate_mk_model(tree = Phylo, 
                                     Q = Q,
                                     root_probabilities = "stationary",
                                     Nsimulations = 1,
                                     drop_dims = F
      )
      states_present <- length(unique(character$tip_states))
      min_frequency <- min(table(character$tip_states)/
                             length(character$tip_states))
    }
  }
  names(character$tip_states) <- Chrono$tip.label
  return(character)
}

data$character <- purrr::pmap(list(data$chronogram,
                                   data$phylogram,
                                   data$character_branch_lengths,
                                   data$Q),
                              sim.binary.markov)


# Extract node states and tip states to columns in dataframe.

data <- unnest_wider(data, character)

# SIMULATE BINARY CHARACTER WITH HIDDEN RATE CATAGORIES-------------------------


hidden.rates.matrix <- function(NStates, RateModel){
  HRM <- get_random_mk_transition_matrix(Nstates = NStates,
                                         rate_model = RateModel)
  HRM[c(4, 7, 10, 13)] <- 0
  HRM[c(1, 6, 11, 16)] <-
    0 - HRM[c(5, 2, 3, 4)] - HRM[c(9, 10, 7, 8)] - HRM[c(13, 14, 15, 12)]
  return(HRM)
}

data$HRM_Q <- replicate(nrow(data),
                        hidden.rates.matrix(4, "ARD"),
                        simplify = F
)


# Generate binary character.

sim.binary.markov.HRM <- function(Chrono, Phylo, BranchLengths, HRMQ){
  states_present <- 1
  min_frequency <- 0
  if(BranchLengths == "chronogram"){
    while(states_present < 2 | min_frequency < 0.05){
      character <- simulate_mk_model(tree = Chrono, 
                                     Q = HRMQ,
                                     root_probabilities = "stationary",
                                     Nsimulations = 1,
                                     drop_dims = F
      )
      character$tip_states[character$tip_states == 2] <- 1
      character$tip_states[character$tip_states > 2] <- 2
      character$node_states[character$node_states ==2] <- 1
      character$node_states[character$node_states > 2] <- 2
      states_present <- length(unique(character$tip_states))
      min_frequency <- min(table(character$tip_states)/
                             length(character$tip_states))
    }
  } else {
    while(states_present < 2 | min_frequency < 0.05){
      character <- simulate_mk_model(tree = Phylo, 
                                     Q = HRMQ,
                                     root_probabilities = "stationary",
                                     Nsimulations = 1,
                                     drop_dims = F
      )
      character$tip_states[character$tip_states == 2] <- 1
      character$tip_states[character$tip_states > 2] <- 2
      character$node_states[character$node_states ==2] <- 1
      character$node_states[character$node_states > 2] <- 2
      states_present <- length(unique(character$tip_states))
      min_frequency <- min(table(character$tip_states)/
                             length(character$tip_states))
    }
  }
  names(character$tip_states) <- Chrono$tip.label
  return(character)
}


# Generate HRM characters.

data$HRM_character <- purrr::pmap(list(data$chronogram,
                                       data$phylogram,
                                       data$character_branch_lengths,
                                       data$HRM_Q),
                                  sim.binary.markov.HRM)


# Extract node states and tip states to columns in dataframe.

data <- hoist(data,
              HRM_character,
              HRM_tip_states = "tip_states",
              HRM_node_states = "node_states")


# SIMULATE BINARY CHARACTER WITH SCALED HIDDEN RATE CATAGORIES------------------


# Create transition matrix for with slow and fast (x100) hidden rate catagories

hidden.rates.matrix.scaled <- function(NStates, RateModel){
  HRM <- get_random_mk_transition_matrix(Nstates = NStates,
                                         rate_model = RateModel)
  HRM[c(4, 7, 10, 13)] <- 0
  HRM[c(8, 14)] <- HRM[c(8, 14)]*100
  HRM[c(1, 6, 11, 16)] <-
    0 - HRM[c(5, 2, 3, 4)] - HRM[c(9, 10, 7, 8)] - HRM[c(13, 14, 15, 12)]
  return(HRM)
}

data$HRMs_Q <- replicate(nrow(data),
                         hidden.rates.matrix.scaled(4, "ARD"),
                         simplify = F
)


# Generate Amplified Hidden Rates characters.

data$HRMs_character <- purrr::pmap(list(data$chronogram,
                                        data$phylogram,
                                        data$character_branch_lengths,
                                        data$HRMs_Q),
                                   sim.binary.markov.HRM)


# Extract node states and tip states to columns in dataframe.

data <- hoist(data,
              HRMs_character,
              HRMs_tip_states = "tip_states",
              HRMs_node_states = "node_states")

save(data, file = "trees&characters.RData")