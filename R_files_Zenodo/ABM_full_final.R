library(ggplot2)
library(reshape2)
library(Rmisc)
library(ggthemes)
library(grid)
library(gridExtra)


#### Global simulation set-up for all regimes ####
set.seed(5)

# Number of Realizations
num_real <- 100
# Number of Time Steps
num_steps <- 150
# Number of initial genotypes in Slow Environments
slow_num_geno <- 5000
# Number of initial genotypes in Fast Environments
fast_num_geno <- 5000
# Numerical state at which an individual becomes reproductive
sexually_mature <- 1.0
# Numberical state at which an individual dies of age
death_state <-  3.0

# Disturbance-associated stage-specific mortality rates
dist_juvmort <- 0.3
dist_adultmort <- 0.00 # adults unaffected by disturbance in current version

# Parent-offspring trait error rate (called 'mutation' just for shorthand)
mu_mutation_rate <- 0.01
f_mutation_rate <- 0.01

# Disturbance Periodicity
freq_slow <- 14 # deterministic slow period
freq_fast <- 3 # deterministic fast period

# Function to create pre-determined set of stochastic slow disturbance sequences
# NB: total # of disturbance equals that of the SLOW regime
slowstoch <- matrix(NA, nrow = num_real, floor(num_steps/freq_slow))
for (i in 1: nrow(slowstoch)){
  slowstoch[i,] <- c(sort(sample(1:num_steps, floor(num_steps/freq_slow))))
}

# Function to created pre-determined set of stochastic disturbance sequences
# NB: total # of disturbance equals that of the FAST regime
faststoch <- matrix(NA, nrow = num_real, floor(num_steps/freq_fast))
for (i in 1: nrow(faststoch)){
  faststoch[i,] <- c(sort(sample(1:num_steps, floor(num_steps/freq_fast))))
}



#### DETERMINISTIC SLOW FLUCTUATION SIMULATION ####

# empty lists to track means and variances of mu and f through time, for all realizations
slowfluc_rlzs_muMEAN <- list()
slowfluc_rlzs_muVAR <- list()
slowfluc_rlzs_fMEAN <- list()
slowfluc_rlzs_fVAR <- list()

# FOR EACH REALIZATION
for (r in 1:num_real){
  
  # sample from fixed range of phenotypic values to set initial distribution of traits
  slow_mu_init <- sample(seq(0.02, 0.04, length.out = slow_num_geno*10), slow_num_geno, replace = T)
  slow_f_init <- c(1/(slow_mu_init*4))
  # sample states such that there is an initial 25:1 ratio of juveniles:adults
  slow_states_init <-  rep(c(runif(25, 0, sexually_mature), runif(1, sexually_mature, death_state)), length.out = slow_num_geno)
  
  sim_slowfluc.ls <- list()
  
  # first time step distribution of traits and states (initial individuals)
  sim_slowfluc.ls[[1]] <- cbind(slow_mu_init, slow_f_init, slow_states_init)
  
  # EVENTS AT EACH TIME STEP FROM HERE ONWARD
  for (i in 1:num_steps){
    mu_distr <- c(sim_slowfluc.ls[[i]][,1]) # copy distribution of mu genotypes from previous time step
    f_distr <- c(sim_slowfluc.ls[[i]][,2]) # copy distribution of f genotypes from previous time step
    states_distr <- c(sim_slowfluc.ls[[i]][,3] + sim_slowfluc.ls[[i]][,1]) # update state of each genotype by growth rate mu
    
    current <- cbind(mu_distr, f_distr, states_distr) #collection of surviving genotypes and their states BEFORE death or repro events
   
    # Juvenile mu & Adult f trading off with respective survival probabilities IN THE ABSENCE OF DISTURBANCE
    juv_indices <- which(current[,3] < sexually_mature)
    juv_mus <- current[juv_indices, 1, drop = FALSE]
    # draw from poisson distribution with lambda = d (juvenile background mortality, which increases as mu increases)
    juv_fates <- cbind(juv_indices, juv_mus, rpois(length(juv_mus), juv_mus*5)) # mu X 5 is roughly the range of juv intrinsic mortality in original analytical model
    juv_fates[,3][which(juv_fates[,3] > 0)] <- 1 # any number of poisson 'events' (non-zero value) counts as death, so just call it '1'
    surviving_juvs <- juv_fates[which(juv_fates[,3] == 0), 1] # indices of surviving juveniles after these intrinsic mortality events
    
    adult_indices <- which(current[,3] >= sexually_mature)
    adult_fs <- current[adult_indices, 2, drop = FALSE]
    # draw from poisson distribution with lambda = d (adult background mortality, which increases as f increases)
    adult_fates <- cbind(adult_indices, adult_fs, rpois(length(adult_fs), adult_fs*0.001)) # f X 0.001 is roughly the range of adult intrinsic mortality in original analytical model
    adult_fates[,3][which(adult_fates[,3] > 0)] <- 1
    surviving_adults <- adult_fates[which(adult_fates[,3] == 0), 1] # indices of surviving adults after these intrinsic mortality events
    
    current <- current[c(surviving_juvs, surviving_adults), , drop = FALSE] # update population with surviving individuals after intrinsic mortality 
    
    
    # Density-dependent mortality
    density <- nrow(current) # current density
    density_mort <- density * (0.001 * (1 + (2 * (density/500000)))) # base mortality * (1 + (density-dependence * (density / scaling factor)))
    density_mort_fates <- sample(density, density_mort) # indices of individuals that die from density-dep mortality
    # update current population with individuals that died from density-dep mortality, if there were any
    if (length(density_mort_fates) > 0){
      current <- current[-density_mort_fates, , drop = FALSE]
    }
    
    
    # Disturbance-caused death events
    if (i%%freq_slow == 0){ # at every time step corresponding to disturbance periodicity
      # juvenile deaths
      too_young <- which(current[,3] < sexually_mature) # indices of individuals that are JUVENILES at current 
      juv_deaths <- sample(too_young, dist_juvmort*length(too_young)) # random dist_juvmort% sample of juveniles to die
      if (length(juv_deaths) != 0){ # as long as there is a non-zero number of juveniles to die,
        current <- current[-juv_deaths, , drop = FALSE] # excise that random % of juveniles
      }
      
      # adult deaths
      old_enough <- which(current[,3] >= sexually_mature) # indices of individuals that are ADULTS at current 
      adult_deaths <- sample(old_enough, dist_adultmort*length(old_enough)) # random dist_adultmort% sample of adults to die
      if (length(adult_deaths) != 0){ # as long as there is a non-zero number of adults to die,
        current <- current[-adult_deaths, , drop = FALSE] # excise that random % of adults
      }
    }
 
    # End of life
    end_of_life <- which(current[,3] >= death_state) # indices of individuals that have reached death state
    if (length(end_of_life) != 0){ # as long as there is a non-zero number of individuals that have reached max life span 
      current <- current[-end_of_life, , drop = FALSE] # excise those individuals
    }
    
    # Birth events
    old_enough <- which(current[,3] >= sexually_mature) # indices of individuals that have reached age of first reproduction & survived
    if (length(old_enough) != 0){ # as long as there are any such old_enough individuals
      fecund_adults <- current[old_enough, , drop = FALSE] # get traits and states of those old_enough individuals
      # NB: due to inheritance error, there are occassionally individuals with 'negative f', which causes errors in Poisson sampling in the next line; therefore make them f = 0
      fecund_adults[,2][which(fecund_adults[,2] < 0)] <- 0 
      fecund_adults <- cbind(fecund_adults, rpois(nrow(fecund_adults), fecund_adults[,2])) # reproduction = make N copies of old_enough individuals, where N = Poisson draw with rate f
      new_births <- fecund_adults[rep(seq(nrow(fecund_adults)), fecund_adults[,4]),] 
      new_births <- new_births[,c(1:3), drop = FALSE]
      
      # Add error to mu and f values in new offspring
      # for each entry in new_births, add normally distributed error with mean 0 and sd that scales to each entry
      new_births[,2] <- 1/(new_births[,1]*4) 
      new_births[,1] <- new_births[,1] + rnorm(length(new_births[,1]), mean = 0, sd = mu_mutation_rate*mean(slow_mu_init)) # add noise to mu of offspring
      new_births[,2] <- new_births[,2] + rnorm(length(new_births[,2]), mean = 0, sd = f_mutation_rate*mean(slow_f_init)) # add noise to f of offspring
      
      new_births[,3] <- 0 # all offspring start at stage = 0
      
      current <- rbind(current, new_births) # append newly birthed offspring to population vector
    }
    
    # after all above events, below is the final current state
    sim_slowfluc.ls[[i+1]] <- current
  }
  
  # store the tracked means and variances for each realization
  slowfluc_rlzs_muMEAN[[r]]<- sapply(sim_slowfluc.ls, function(x) mean(x[,1]))
  slowfluc_rlzs_muVAR[[r]]<- sapply(sim_slowfluc.ls, function(x) var(x[,1]))
  slowfluc_rlzs_fMEAN[[r]] <- sapply(sim_slowfluc.ls, function(x) mean(x[,2]))
  slowfluc_rlzs_fVAR[[r]] <- sapply(sim_slowfluc.ls, function(x) var(x[,2]))
  
}

####------------------------------------------------------------------------------------------------####
# All below simulations follow above example, therefore refer to detailed comments within above simulation
# All below simulations only differ from above in fluctuation regime
#------------------------------------------------------------------------------------------------------#

#### DETERMINISTIC FAST FLUCTUATION SIMULATION #### 
fastfluc_rlzs_muMEAN <- list()
fastfluc_rlzs_muVAR <- list()
fastfluc_rlzs_fMEAN <- list()
fastfluc_rlzs_fVAR <- list()

for (r in 1:num_real){
  
  fast_mu_init <- sample(seq(0.02, 0.04, length.out = fast_num_geno*10), fast_num_geno, replace = T)
  fast_f_init <- c(1/(fast_mu_init*4))
  fast_states_init <-  rep(c(runif(25, 0, sexually_mature), runif(1, sexually_mature, death_state)), length.out = fast_num_geno)
  
  sim_fastfluc.ls <- list()
  sim_fastfluc.ls[[1]] <- cbind(fast_mu_init, fast_f_init, fast_states_init)
  
  for (i in 1:num_steps){
    mu_distr <- c(sim_fastfluc.ls[[i]][,1]) 
    f_distr <- c(sim_fastfluc.ls[[i]][,2]) 
    states_distr <- c(sim_fastfluc.ls[[i]][,3] + sim_fastfluc.ls[[i]][,1]) 
    
    current <- cbind(mu_distr, f_distr, states_distr) 
    
    # Juvenile mu & Adult f trading off with respective survival probabilities IN THE ABSENCE OF DISTURBANCE
    juv_indices <- which(current[,3] < sexually_mature)
    juv_mus <- current[juv_indices, 1, drop = FALSE]
    juv_fates <- cbind(juv_indices, juv_mus, rpois(length(juv_mus), juv_mus*5)) 
    juv_fates[,3][which(juv_fates[,3] > 0)] <- 1 
    surviving_juvs <- juv_fates[which(juv_fates[,3] == 0), 1]
    
    adult_indices <- which(current[,3] >= sexually_mature)
    adult_fs <- current[adult_indices, 2, drop = FALSE]
    adult_fates <- cbind(adult_indices, adult_fs, rpois(length(adult_fs), adult_fs*0.001))
    adult_fates[,3][which(adult_fates[,3] > 0)] <- 1
    surviving_adults <- adult_fates[which(adult_fates[,3] == 0), 1]
    
    current <- current[c(surviving_juvs, surviving_adults), , drop = FALSE]
    
    # Density-dependent mortality
    density <- nrow(current)
    density_mort <- density * (0.001 * (1 + (2 * (density/500000)))) 
    density_mort_fates <- sample(density, density_mort)
    
    if (length(density_mort_fates) > 0){
      current <- current[-density_mort_fates, , drop = FALSE]
    }
    
    # Disturbance-caused death events
    if (i%%freq_fast == 0){ 
      # juvenile deaths
      too_young <- which(current[,3] < sexually_mature) 
      juv_deaths <- sample(too_young, dist_juvmort*length(too_young)) 
      if (length(juv_deaths) != 0){ 
        current <- current[-juv_deaths, , drop = FALSE] 
      }
      
      # adult deaths
      old_enough <- which(current[,3] >= sexually_mature) 
      adult_deaths <- sample(old_enough, dist_adultmort*length(old_enough)) 
      if (length(adult_deaths) != 0){ 
        current <- current[-adult_deaths, , drop = FALSE] 
      }
    }
    
    # End of life
    end_of_life <- which(current[,3] >= death_state) 
    if (length(end_of_life) != 0){
      current <- current[-end_of_life,, drop = FALSE] 
    }
    
    # Birth events
    old_enough <- which(current[,3] >= sexually_mature) 
    if (length(old_enough) != 0){ 
      fecund_adults <- current[old_enough, , drop = FALSE] 
      fecund_adults[,2][which(fecund_adults[,2] < 0)] <- 0 
      fecund_adults <- cbind(fecund_adults, rpois(nrow(fecund_adults), fecund_adults[,2]))
      new_births <- fecund_adults[rep(seq(nrow(fecund_adults)), fecund_adults[,4]), ] 
      new_births <- new_births[,1:3, drop = FALSE]
      
      # add error to mu and f values in new offspring
      new_births[,2] <- 1/(new_births[,1]*4) 
      new_births[,1] <- new_births[,1] + rnorm(length(new_births[,1]), mean = 0, sd = mu_mutation_rate*mean(fast_mu_init)) 
      new_births[,2] <- new_births[,2] + rnorm(length(new_births[,2]), mean = 0, sd = f_mutation_rate*mean(fast_f_init))
      
      new_births[,3] <- 0 
      
      current <- rbind(current, new_births)
    }
    
    sim_fastfluc.ls[[i+1]] <- current
  }
  
  fastfluc_rlzs_muMEAN[[r]] <- sapply(sim_fastfluc.ls, function(x) mean(x[,1]))
  fastfluc_rlzs_muVAR[[r]] <- sapply(sim_fastfluc.ls, function(x) var(x[,1]))
  fastfluc_rlzs_fMEAN[[r]] <- sapply(sim_fastfluc.ls, function(x) mean(x[,2]))
  fastfluc_rlzs_fVAR[[r]] <- sapply(sim_fastfluc.ls, function(x) var(x[,2]))
  
}

#------------------------------------------------------------------------------------------------------#

#### STOCHASTIC SLOW FLUCTUATION SIMULATION ####
slowstoch_rlzs_muMEAN <- list()
slowstoch_rlzs_muVAR <- list()
slowstoch_rlzs_fMEAN <- list()
slowstoch_rlzs_fVAR <- list()

for (r in 1:num_real){
  
  slowstoch_mu_init <- sample(seq(0.02, 0.04, length.out = slow_num_geno*10), slow_num_geno, replace = T)
  slowstoch_f_init <- c(1/(slow_mu_init*4))
  slowstoch_states_init <-  rep(c(runif(25, 0, sexually_mature), runif(1, sexually_mature, death_state)), length.out = slow_num_geno)
  
  sim_slowstoch.ls <- list()
  sim_slowstoch.ls[[1]] <- cbind(slowstoch_mu_init, slowstoch_f_init, slowstoch_states_init)
  
  for (i in 1:num_steps){
    mu_distr <- c(sim_slowstoch.ls[[i]][,1]) 
    f_distr <- c(sim_slowstoch.ls[[i]][,2]) 
    states_distr <- c(sim_slowstoch.ls[[i]][,3] + sim_slowstoch.ls[[i]][,1]) 
    
    current <- cbind(mu_distr, f_distr, states_distr) 
    
    # Juvenile mu & Adult f trading off with respective survival probabilities IN THE ABSENCE OF DISTURBANCE
    juv_indices <- which(current[,3] < sexually_mature)
    juv_mus <- current[juv_indices, 1, drop = FALSE]
    juv_fates <- cbind(juv_indices, juv_mus, rpois(length(juv_mus), juv_mus*5)) 
    juv_fates[,3][which(juv_fates[,3] > 0)] <- 1 
    surviving_juvs <- juv_fates[which(juv_fates[,3] == 0), 1]
    
    adult_indices <- which(current[,3] >= sexually_mature)
    adult_fs <- current[adult_indices, 2, drop = FALSE]
    adult_fates <- cbind(adult_indices, adult_fs, rpois(length(adult_fs), adult_fs*0.001))
    adult_fates[,3][which(adult_fates[,3] > 0)] <- 1
    surviving_adults <- adult_fates[which(adult_fates[,3] == 0), 1]
    
    current <- current[c(surviving_juvs, surviving_adults), , drop = FALSE]
    
    # Density-dependent mortality
    density <- nrow(current)
    density_mort <- density * (0.001 * (1 + (2 * (density/500000)))) 
    density_mort_fates <- sample(density, density_mort)
    
    if (length(density_mort_fates) > 0){
      current <- current[-density_mort_fates, , drop = FALSE]
    }
    
    # Disturbance-caused death events
    if (i %in% slowstoch[r,]){ 
      #juvenile deaths
      too_young <- which(current[,3] < sexually_mature) 
      juv_deaths <- sample(too_young, dist_juvmort*length(too_young)) 
      if (length(juv_deaths) != 0){ 
        current <- current[-juv_deaths, , drop = FALSE] 
      }
      
      #adult deaths
      old_enough <- which(current[,3] >= sexually_mature) 
      adult_deaths <- sample(old_enough, dist_adultmort*length(old_enough))
      if (length(adult_deaths) != 0){ 
        current <- current[-adult_deaths, , drop = FALSE] 
      }
    }
    
    # End of life
    end_of_life <- which(current[,3] >= death_state) 
    if (length(end_of_life) != 0){ 
      current <- current[-end_of_life, , drop = FALSE] 
    }
    
    # Birth events
    old_enough <- which(current[,3] >= sexually_mature) 
    if (length(old_enough) != 0){ 
      fecund_adults <- current[old_enough, , drop = FALSE] 
      fecund_adults[,2][which(fecund_adults[,2] < 0)] <- 0 
      fecund_adults <- cbind(fecund_adults, rpois(nrow(fecund_adults), fecund_adults[,2]))
      new_births <- fecund_adults[rep(seq(nrow(fecund_adults)), fecund_adults[,4]),] 
      if(nrow(new_births) > 0){
        new_births <- new_births[,c(1:3), drop = FALSE]
      }
      
      # add error to mu and f values in new offspring
      new_births[,2] <- 1/(new_births[,1]*4) 
      new_births[,1] <- new_births[,1] + rnorm(length(new_births[,1]), mean = 0, sd = mu_mutation_rate*mean(slow_mu_init)) 
      new_births[,2] <- new_births[,2] + rnorm(length(new_births[,2]), mean = 0, sd = f_mutation_rate*mean(slow_f_init))
      
      new_births[,3] <- 0
      
      current <- rbind(current, new_births)
    }
    
    sim_slowstoch.ls[[i+1]] <- current
  }
  
  slowstoch_rlzs_muMEAN[[r]]<- sapply(sim_slowstoch.ls, function(x) mean(x[,1]))
  slowstoch_rlzs_muVAR[[r]]<- sapply(sim_slowstoch.ls, function(x) var(x[,1]))
  slowstoch_rlzs_fMEAN[[r]] <- sapply(sim_slowstoch.ls, function(x) mean(x[,2]))
  slowstoch_rlzs_fVAR[[r]] <- sapply(sim_slowstoch.ls, function(x) var(x[,2]))
  
}


#------------------------------------------------------------------------------------------------------#

#### STOCHASTIC FAST FLUCTUATION SIMULATION ####
faststoch_rlzs_muMEAN <- list()
faststoch_rlzs_muVAR <- list()
faststoch_rlzs_fMEAN <- list()
faststoch_rlzs_fVAR <- list()

for (r in 1:num_real){
  
  faststoch_mu_init <- sample(seq(0.02, 0.04, length.out = fast_num_geno*10), fast_num_geno, replace = T)
  faststoch_f_init <- c(1/(fast_mu_init*4))
  faststoch_states_init <-  rep(c(runif(25, 0, sexually_mature), runif(1, sexually_mature, death_state)), length.out = fast_num_geno)
  
  sim_faststoch.ls <- list()
  sim_faststoch.ls[[1]] <- cbind(faststoch_mu_init, faststoch_f_init, faststoch_states_init)
  
  for (i in 1:num_steps){
    mu_distr <- c(sim_faststoch.ls[[i]][,1]) 
    f_distr <- c(sim_faststoch.ls[[i]][,2]) 
    states_distr <- c(sim_faststoch.ls[[i]][,3] + sim_faststoch.ls[[i]][,1]) 
    
    current <- cbind(mu_distr, f_distr, states_distr) 
    
    # Juvenile mu & Adult f trading off with respective survival probabilities IN THE ABSENCE OF DISTURBANCE
    juv_indices <- which(current[,3] < sexually_mature)
    juv_mus <- current[juv_indices, 1, drop = FALSE]
    juv_fates <- cbind(juv_indices, juv_mus, rpois(length(juv_mus), juv_mus*5)) 
    juv_fates[,3][which(juv_fates[,3] > 0)] <- 1 
    surviving_juvs <- juv_fates[which(juv_fates[,3] == 0), 1]
    
    adult_indices <- which(current[,3] >= sexually_mature)
    adult_fs <- current[adult_indices, 2, drop = FALSE]
    adult_fates <- cbind(adult_indices, adult_fs, rpois(length(adult_fs), adult_fs*0.001))
    adult_fates[,3][which(adult_fates[,3] > 0)] <- 1
    surviving_adults <- adult_fates[which(adult_fates[,3] == 0), 1]
    
    current <- current[c(surviving_juvs, surviving_adults), , drop = FALSE]
    
    # Density-dependent mortality
    density <- nrow(current)
    density_mort <- density * (0.001 * (1 + (2 * (density/500000)))) 
    density_mort_fates <- sample(density, density_mort)
    
    if (length(density_mort_fates) > 0){
      current <- current[-density_mort_fates, , drop = FALSE]
    }
    
    # Disturbance-caused death events
    if (i %in% faststoch[r,]){ 
      #juvenile deaths
      too_young <- which(current[,3] < sexually_mature) 
      juv_deaths <- sample(too_young, dist_juvmort*length(too_young)) 
      if (length(juv_deaths) != 0){ 
        current <- current[-juv_deaths, , drop = FALSE] 
      }
      
      #adult deaths
      old_enough <- which(current[,3] >= sexually_mature) 
      adult_deaths <- sample(old_enough, dist_adultmort*length(old_enough))
      if (length(adult_deaths) != 0){ 
        current <- current[-adult_deaths, , drop = FALSE] 
      }
    }
    
    # End of life
    end_of_life <- which(current[,3] >= death_state) 
    if (length(end_of_life) != 0){ 
      current <- current[-end_of_life, , drop = FALSE]
    }
    
    # Birth events
    old_enough <- which(current[,3] >= sexually_mature) 
    if (length(old_enough) != 0){ 
      fecund_adults <- current[old_enough, , drop = FALSE] 
      fecund_adults[,2][which(fecund_adults[,2] < 0)] <- 0
      fecund_adults <- cbind(fecund_adults, rpois(nrow(fecund_adults), fecund_adults[,2]))
      new_births <- fecund_adults[rep(seq(nrow(fecund_adults)), fecund_adults[,4]),] 
      if(nrow(new_births) > 0){
        new_births <- new_births[,c(1:3), drop = FALSE]
      }
      
      # add error to mu and f values in new offspring
      new_births[,2] <- 1/(new_births[,1]*4) 
      new_births[,1] <- new_births[,1] + rnorm(length(new_births[,1]), mean = 0, sd = mu_mutation_rate*mean(slow_mu_init)) 
      new_births[,2] <- new_births[,2] + rnorm(length(new_births[,2]), mean = 0, sd = f_mutation_rate*mean(slow_f_init))
      
      new_births[,3] <- 0 
      
      current <- rbind(current, new_births)
    }
   
    sim_faststoch.ls[[i+1]] <- current
  }
  
  faststoch_rlzs_muMEAN[[r]]<- sapply(sim_faststoch.ls, function(x) mean(x[,1]))
  faststoch_rlzs_muVAR[[r]]<- sapply(sim_faststoch.ls, function(x) var(x[,1]))
  faststoch_rlzs_fMEAN[[r]] <- sapply(sim_faststoch.ls, function(x) mean(x[,2]))
  faststoch_rlzs_fVAR[[r]] <- sapply(sim_faststoch.ls, function(x) var(x[,2]))
  
}

#------------------------------------------------------------------------------------------------------#

#### Simulation data cleaning ####

# Maturation (mu) sim data #
# Deterministic Slow Sim - MU means and variances through time across realizations
slowfluc_rlzs_muMEAN.df <- data.frame(matrix(unlist(slowfluc_rlzs_muMEAN), nrow = num_real, byrow = T))
colnames(slowfluc_rlzs_muMEAN.df) <- c(1:(num_steps+1))
slowfluc_rlzs_muMEAN.df <- cbind("realizations" = c(1:num_real), slowfluc_rlzs_muMEAN.df)
slowfluc_rlzs_muMEAN_melt <- melt(slowfluc_rlzs_muMEAN.df, id.vars = "realizations")

slowfluc_rlzs_muVAR.df <- data.frame(matrix(unlist(slowfluc_rlzs_muVAR), nrow = num_real, byrow = T))
colnames(slowfluc_rlzs_muVAR.df) <- c(1:(num_steps+1))
slowfluc_rlzs_muVAR.df <- cbind("realizations" = c(1:num_real), slowfluc_rlzs_muVAR.df)
slowfluc_rlzs_muVAR_melt <- melt(slowfluc_rlzs_muVAR.df, id.vars = "realizations")

# Deterministic Fast Sim - MU means and variances through time across realizations
fastfluc_rlzs_muMEAN.df <- data.frame(matrix(unlist(fastfluc_rlzs_muMEAN), nrow = num_real, byrow = T))
colnames(fastfluc_rlzs_muMEAN.df) <- c(1:(num_steps+1))
fastfluc_rlzs_muMEAN.df <- cbind("realizations" = c(1:num_real), fastfluc_rlzs_muMEAN.df)
fastfluc_rlzs_muMEAN_melt <- melt(fastfluc_rlzs_muMEAN.df, id.vars = "realizations")

fastfluc_rlzs_muVAR.df <- data.frame(matrix(unlist(fastfluc_rlzs_muVAR), nrow = num_real, byrow = T))
colnames(fastfluc_rlzs_muVAR.df) <- c(1:(num_steps+1))
fastfluc_rlzs_muVAR.df <- cbind("realizations" = c(1:num_real), fastfluc_rlzs_muVAR.df)
fastfluc_rlzs_muVAR_melt <- melt(fastfluc_rlzs_muVAR.df, id.vars = "realizations")

# Stochastic Slow Sim - MU means and variances through time across realizations
slowstoch_rlzs_muMEAN.df <- data.frame(matrix(unlist(slowstoch_rlzs_muMEAN), nrow = num_real, byrow = T))
colnames(slowstoch_rlzs_muMEAN.df) <- c(1:(num_steps+1))
slowstoch_rlzs_muMEAN.df <- cbind("realizations" = c(1:num_real), slowstoch_rlzs_muMEAN.df)
slowstoch_rlzs_muMEAN_melt <- melt(slowstoch_rlzs_muMEAN.df, id.vars = "realizations")

slowstoch_rlzs_muVAR.df <- data.frame(matrix(unlist(slowstoch_rlzs_muVAR), nrow = num_real, byrow = T))
colnames(slowstoch_rlzs_muVAR.df) <- c(1:(num_steps+1))
slowstoch_rlzs_muVAR.df <- cbind("realizations" = c(1:num_real), slowstoch_rlzs_muVAR.df)
slowstoch_rlzs_muVAR_melt <- melt(slowstoch_rlzs_muVAR.df, id.vars = "realizations")


# Stochastic Fast Sim - MU means and variances through time across realizations
faststoch_rlzs_muMEAN.df <- data.frame(matrix(unlist(faststoch_rlzs_muMEAN), nrow = num_real, byrow = T))
colnames(faststoch_rlzs_muMEAN.df) <- c(1:(num_steps+1))
faststoch_rlzs_muMEAN.df <- cbind("realizations" = c(1:num_real), faststoch_rlzs_muMEAN.df)
faststoch_rlzs_muMEAN_melt <- melt(faststoch_rlzs_muMEAN.df, id.vars = "realizations")

faststoch_rlzs_muVAR.df <- data.frame(matrix(unlist(faststoch_rlzs_muVAR), nrow = num_real, byrow = T))
colnames(faststoch_rlzs_muVAR.df) <- c(1:(num_steps+1))
faststoch_rlzs_muVAR.df <- cbind("realizations" = c(1:num_real), faststoch_rlzs_muVAR.df)
faststoch_rlzs_muVAR_melt <- melt(faststoch_rlzs_muVAR.df, id.vars = "realizations")


# Fecundity (f) sim data #
# Deterministic Slow Sim - F means and variances through time across realizations
slowfluc_rlzs_fMEAN.df <- data.frame(matrix(unlist(slowfluc_rlzs_fMEAN), nrow = num_real, byrow = T))
colnames(slowfluc_rlzs_fMEAN.df) <- c(1:(num_steps+1))
slowfluc_rlzs_fMEAN.df <- cbind("realizations" = c(1:num_real), slowfluc_rlzs_fMEAN.df)
slowfluc_rlzs_fMEAN_melt <- melt(slowfluc_rlzs_fMEAN.df, id.vars = "realizations")

slowfluc_rlzs_fVAR.df <- data.frame(matrix(unlist(slowfluc_rlzs_fVAR), nrow = num_real, byrow = T))
colnames(slowfluc_rlzs_fVAR.df) <- c(1:(num_steps+1))
slowfluc_rlzs_fVAR.df <- cbind("realizations" = c(1:num_real), slowfluc_rlzs_fVAR.df)
slowfluc_rlzs_fVAR_melt <- melt(slowfluc_rlzs_fVAR.df, id.vars = "realizations")

# Deterministic Fast Sim - F means and variances through time across realizations
fastfluc_rlzs_fMEAN.df <- data.frame(matrix(unlist(fastfluc_rlzs_fMEAN), nrow = num_real, byrow = T))
colnames(fastfluc_rlzs_fMEAN.df) <- c(1:(num_steps+1))
fastfluc_rlzs_fMEAN.df <- cbind("realizations" = c(1:num_real), fastfluc_rlzs_fMEAN.df)
fastfluc_rlzs_fMEAN_melt <- melt(fastfluc_rlzs_fMEAN.df, id.vars = "realizations")

fastfluc_rlzs_fVAR.df <- data.frame(matrix(unlist(fastfluc_rlzs_fVAR), nrow = num_real, byrow = T))
colnames(fastfluc_rlzs_fVAR.df) <- c(1:(num_steps+1))
fastfluc_rlzs_fVAR.df <- cbind("realizations" = c(1:num_real), fastfluc_rlzs_fVAR.df)
fastfluc_rlzs_fVAR_melt <- melt(fastfluc_rlzs_fVAR.df, id.vars = "realizations")

# Stochastic Slow Sim - F means and variances through time across realizations
slowstoch_rlzs_fMEAN.df <- data.frame(matrix(unlist(slowstoch_rlzs_fMEAN), nrow = num_real, byrow = T))
colnames(slowstoch_rlzs_fMEAN.df) <- c(1:(num_steps+1))
slowstoch_rlzs_fMEAN.df <- cbind("realizations" = c(1:num_real), slowstoch_rlzs_fMEAN.df)
slowstoch_rlzs_fMEAN_melt <- melt(slowstoch_rlzs_fMEAN.df, id.vars = "realizations")

slowstoch_rlzs_fVAR.df <- data.frame(matrix(unlist(slowstoch_rlzs_fVAR), nrow = num_real, byrow = T))
colnames(slowstoch_rlzs_fVAR.df) <- c(1:(num_steps+1))
slowstoch_rlzs_fVAR.df <- cbind("realizations" = c(1:num_real), slowstoch_rlzs_fVAR.df)
slowstoch_rlzs_fVAR_melt <- melt(slowstoch_rlzs_fVAR.df, id.vars = "realizations")

# Stochastic Fast Sim - F means and variances through time across realizations
faststoch_rlzs_fMEAN.df <- data.frame(matrix(unlist(faststoch_rlzs_fMEAN), nrow = num_real, byrow = T))
colnames(faststoch_rlzs_fMEAN.df) <- c(1:(num_steps+1))
faststoch_rlzs_fMEAN.df <- cbind("realizations" = c(1:num_real), faststoch_rlzs_fMEAN.df)
faststoch_rlzs_fMEAN_melt <- melt(faststoch_rlzs_fMEAN.df, id.vars = "realizations")

faststoch_rlzs_fVAR.df <- data.frame(matrix(unlist(faststoch_rlzs_fVAR), nrow = num_real, byrow = T))
colnames(faststoch_rlzs_fVAR.df) <- c(1:(num_steps+1))
faststoch_rlzs_fVAR.df <- cbind("realizations" = c(1:num_real), faststoch_rlzs_fVAR.df)
faststoch_rlzs_fVAR_melt <- melt(faststoch_rlzs_fVAR.df, id.vars = "realizations")

#### Simulation visualization plots ####

# Tracking mean of MU through time across all realizations in DETERMINISTIC Slow vs. Fast regimes
rlzs_muMEAN_deterministic.p <- ggplot() + 
  geom_line(data = slowfluc_rlzs_muMEAN_melt, aes(x = variable, y = value, group = realizations, color = "Deterministic Slow"), size = 0.3, alpha = 0.6) +
  geom_line(data = fastfluc_rlzs_muMEAN_melt, aes(x = variable, y = value, group = realizations, color = "Deterministic Fast"), size = 0.3, alpha = 0.6) +
  scale_x_discrete(breaks = seq(0, length(slowfluc_rlzs_muMEAN_melt$variable), 10)) +
  scale_color_manual(name = "Fluctuation Regimes", values = c("#EB95FF", "#5AFFC3")) +
  guides(color = guide_legend(title.position = 'top', override.aes=list(size=3))) +
  ggtitle(expression(paste("MEAN of ", italic("\\u03BC"), " within populations"))) +
  ylab(expression(paste("Mean (", italic("\\u03BC"), " )"))) +
  xlab("Time") +
  ylim(min(min(fastfluc_rlzs_muMEAN_melt$value, na.rm = T), min(faststoch_rlzs_muMEAN_melt$value, na.rm = T), min(slowfluc_rlzs_muMEAN_melt$value, na.rm = T), min(slowstoch_rlzs_muMEAN_melt$value, na.rm = T)),
       max(max(fastfluc_rlzs_muMEAN_melt$value, na.rm = T), max(faststoch_rlzs_muMEAN_melt$value, na.rm = T), max(slowfluc_rlzs_muMEAN_melt$value, na.rm = T), max(slowstoch_rlzs_muMEAN_melt$value, na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white", family = "Avenir"),
        legend.position = c(0.25,0.80),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, color = "white", family = "Avenir", face = "bold"),
        legend.text = element_text(size = 13, colour = "white", family = "Avenir", face = "italic"),
        axis.title = element_text(size = 16, color = "white", family = "Avenir"),
        axis.text.x = element_text(size = 11, colour = "gray70", face = "italic"),
        axis.text.y = element_text(size = 14, colour = "gray70", face = "italic"),
        axis.ticks = element_line(colour = "gray30"),
        axis.line = element_line(colour = "gray30"),
        panel.background = element_rect("black"),
        plot.background = element_rect("black"),
        plot.margin = margin(10,15,10,10),
        panel.grid = element_blank())

# Tracking mean of MU through time across all realizations in STOCHASTIC Slow vs. Fast regimes
rlzs_muMEAN_stochastic.p <- ggplot() + 
  geom_line(data = slowstoch_rlzs_muMEAN_melt, aes(x = variable, y = value, group = realizations, color = "Stochastic Slow"), size = 0.3, alpha = 0.6) +
  geom_line(data = faststoch_rlzs_muMEAN_melt, aes(x = variable, y = value, group = realizations, color = "Stochastic Fast"), size = 0.3, alpha = 0.6) +
  scale_x_discrete(breaks = seq(0, length(slowstoch_rlzs_muMEAN_melt$variable), 10)) +
  scale_color_manual(name = "Fluctuation Regimes", values = c("#BA78F3", "#36F2F9")) +
  guides(color = guide_legend(title.position = 'top', override.aes=list(size=3))) +
  ggtitle(expression(paste("MEAN of ", italic("\\u03BC"), " within populations"))) +
  ylab(expression(paste("Mean (", italic("\\u03BC"), " )"))) +
  xlab("Time") +
  ylim(min(min(fastfluc_rlzs_muMEAN_melt$value, na.rm = T), min(faststoch_rlzs_muMEAN_melt$value, na.rm = T), min(slowfluc_rlzs_muMEAN_melt$value, na.rm = T), min(slowstoch_rlzs_muMEAN_melt$value, na.rm = T)),
       max(max(fastfluc_rlzs_muMEAN_melt$value, na.rm = T), max(faststoch_rlzs_muMEAN_melt$value, na.rm = T), max(slowfluc_rlzs_muMEAN_melt$value, na.rm = T), max(slowstoch_rlzs_muMEAN_melt$value, na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white", family = "Avenir"),
        legend.position = c(0.25,0.80),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, color = "white", family = "Avenir", face = "bold"),
        legend.text = element_text(size = 13, colour = "white", family = "Avenir", face = "italic"),
        axis.title = element_text(size = 16, color = "white", family = "Avenir"),
        axis.text.x = element_text(size = 11, colour = "gray70", face = "italic"),
        axis.text.y = element_text(size = 14, colour = "gray70", face = "italic"),
        axis.ticks = element_line(colour = "gray30"),
        axis.line = element_line(colour = "gray30"),
        panel.background = element_rect("black"),
        plot.background = element_rect("black"),
        plot.margin = margin(10,15,10,10),
        panel.grid = element_blank())

# Tracking variance of MU through time across all realizations in DETERMINISTIC Slow vs. Fast regimes
rlzs_muVAR_deterministic.p <- ggplot() + 
  geom_line(data = slowfluc_rlzs_muVAR_melt, aes(x = variable, y = value, group = realizations, color = "Deterministic Slow"), size = 0.1, alpha = 0.6) +
  geom_line(data = fastfluc_rlzs_muVAR_melt, aes(x = variable, y = value, group = realizations, color = "Deterministic Fast"), size = 0.1, alpha = 0.6) +
  scale_x_discrete(breaks = seq(0, length(slowfluc_rlzs_muVAR_melt$variable), 10)) +
  scale_color_manual(name = "Fluctuation Regimes", values = c("#EB95FF", "#5AFFC3")) +
  guides(color = guide_legend(title.position = 'top', override.aes=list(size=3))) +
  ggtitle(expression(paste("VARIANCE of ", italic("\\u03BC"), " within populations"))) +
  ylab(expression(paste("Var (", mu, ")")))+
  xlab("Time") +
  ylim(min(min(fastfluc_rlzs_muVAR_melt$value, na.rm = T), min(faststoch_rlzs_muVAR_melt$value, na.rm = T), min(slowfluc_rlzs_muVAR_melt$value, na.rm = T), min(slowstoch_rlzs_muVAR_melt$value, na.rm = T)),
       max(max(fastfluc_rlzs_muVAR_melt$value, na.rm = T), max(faststoch_rlzs_muVAR_melt$value, na.rm = T), max(slowfluc_rlzs_muVAR_melt$value, na.rm = T), max(slowstoch_rlzs_muVAR_melt$value, na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white", family = "Avenir"),
        legend.position = c(0.30,0.80),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, color = "white", family = "Avenir", face = "bold"),
        legend.text = element_text(size = 13, colour = "white", family = "Avenir", face = "italic"),
        axis.title = element_text(size = 16, color = "white", family = "Avenir"),
        axis.text.x = element_text(size = 12, colour = "gray70", face = "italic"),
        axis.text.y = element_text(size = 14, colour = "gray70", face = "italic"),
        axis.ticks = element_line(colour = "gray30"),
        axis.line = element_line(colour = "gray30"),
        panel.background = element_rect("black"),
        plot.background = element_rect("black"),
        plot.margin = margin(10,15,10,10),
        panel.grid = element_blank())

# Tracking variance of MU through time across all realizations in STOCHASTIC Slow vs. Fast regimes
rlzs_muVAR_stochastic.p <- ggplot() + 
  geom_line(data = slowstoch_rlzs_muVAR_melt, aes(x = variable, y = value, group = realizations, color = "Stochastic Slow"), size = 0.3, alpha = 0.6) +
  geom_line(data = faststoch_rlzs_muVAR_melt, aes(x = variable, y = value, group = realizations, color = "Stochastic Fast"), size = 0.3, alpha = 0.6) +
  scale_x_discrete(breaks = seq(0, length(slowstoch_rlzs_muVAR_melt$variable), 10)) +
  scale_color_manual(name = "Fluctuation Regimes", values = c("#BA78F3", "#36F2F9")) +
  guides(color = guide_legend(title.position = 'top', override.aes=list(size=3))) +
  ggtitle(expression(paste("VARIANCE of ", italic("\\u03BC"), " within populations"))) +
  ylab(expression(paste("Var (", mu, ")")))+
  xlab("Time") +
  ylim(min(min(fastfluc_rlzs_muVAR_melt$value, na.rm = T), min(faststoch_rlzs_muVAR_melt$value, na.rm = T), min(slowfluc_rlzs_muVAR_melt$value, na.rm = T), min(slowstoch_rlzs_muVAR_melt$value, na.rm = T)),
       max(max(fastfluc_rlzs_muVAR_melt$value, na.rm = T), max(faststoch_rlzs_muVAR_melt$value, na.rm = T), max(slowfluc_rlzs_muVAR_melt$value, na.rm = T), max(slowstoch_rlzs_muVAR_melt$value, na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white", family = "Avenir"),
        legend.position = c(0.30,0.80),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, color = "white", family = "Avenir", face = "bold"),
        legend.text = element_text(size = 13, colour = "white", family = "Avenir", face = "italic"),
        axis.title = element_text(size = 16, color = "white", family = "Avenir"),
        axis.text.x = element_text(size = 12, colour = "gray70", face = "italic"),
        axis.text.y = element_text(size = 14, colour = "gray70", face = "italic"),
        axis.ticks = element_line(colour = "gray30"),
        axis.line = element_line(colour = "gray30"),
        panel.background = element_rect("black"),
        plot.background = element_rect("black"),
        plot.margin = margin(10,15,10,10),
        panel.grid = element_blank())




# Tracking mean of F through time across all realizations in DETERMINISTIC Slow vs. Fast regimes
rlzs_fMEAN_deterministic.p <- ggplot() + 
  geom_line(data = slowfluc_rlzs_fMEAN_melt, aes(x = variable, y = value, group = realizations, color = "Deterministic Slow"), size = 0.3, alpha = 0.6) +
  geom_line(data = fastfluc_rlzs_fMEAN_melt, aes(x = variable, y = value, group = realizations, color = "Deterministic Fast"), size = 0.3, alpha = 0.6) +
  scale_x_discrete(breaks = seq(0, length(slowfluc_rlzs_fMEAN_melt$variable), 10)) +
  scale_color_manual(name = "Fluctuation Regimes", values = c("#EB95FF", "#5AFFC3")) +
  guides(color = guide_legend(title.position = 'top', override.aes=list(size=3))) +
  ggtitle(expression(paste("MEAN of ", italic(f), " within populations"))) +
  ylab(expression(paste("Mean (", italic("f"), " )"))) +
  xlab("Time") +
  ylim(min(min(fastfluc_rlzs_fMEAN_melt$value, na.rm = T), min(faststoch_rlzs_fMEAN_melt$value, na.rm = T), min(slowfluc_rlzs_fMEAN_melt$value, na.rm = T), min(slowstoch_rlzs_fMEAN_melt$value, na.rm = T)),
       max(max(fastfluc_rlzs_fMEAN_melt$value, na.rm = T), max(faststoch_rlzs_fMEAN_melt$value, na.rm = T), max(slowfluc_rlzs_fMEAN_melt$value, na.rm = T), max(slowstoch_rlzs_fMEAN_melt$value, na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white", family = "Avenir"),
        legend.position = c(0.25,0.80),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, color = "white", family = "Avenir", face = "bold"),
        legend.text = element_text(size = 13, colour = "white", family = "Avenir", face = "italic"),
        axis.title = element_text(size = 16, color = "white", family = "Avenir"),
        axis.text.x = element_text(size = 11, colour = "gray70", face = "italic"),
        axis.text.y = element_text(size = 14, colour = "gray70", face = "italic"),
        axis.ticks = element_line(colour = "gray30"),
        axis.line = element_line(colour = "gray30"),
        panel.background = element_rect("black"),
        plot.background = element_rect("black"),
        plot.margin = margin(10,15,10,10),
        panel.grid = element_blank())

# Tracking mean of F through time across all realizations in STOCHASTIC Slow vs. Fast regimes
rlzs_fMEAN_stochastic.p <- ggplot() + 
  geom_line(data = slowstoch_rlzs_fMEAN_melt, aes(x = variable, y = value, group = realizations, color = "Stochastic Slow"), size = 0.3, alpha = 0.6) +
  geom_line(data = faststoch_rlzs_fMEAN_melt, aes(x = variable, y = value, group = realizations, color = "Stochastic Fast"), size = 0.3, alpha = 0.6) +
  scale_x_discrete(breaks = seq(0, length(slowstoch_rlzs_fMEAN_melt$variable), 10)) +
  scale_color_manual(name = "Fluctuation Regimes", values = c("#BA78F3", "#36F2F9")) +
  guides(color = guide_legend(title.position = 'top', override.aes=list(size=3))) +
  ggtitle(expression(paste("MEAN of ", italic(f), " within populations"))) +
  ylab(expression(paste("Mean (", italic("f"), " )"))) +
  xlab("Time") +
  ylim(min(min(fastfluc_rlzs_fMEAN_melt$value, na.rm = T), min(faststoch_rlzs_fMEAN_melt$value, na.rm = T), min(slowfluc_rlzs_fMEAN_melt$value, na.rm = T), min(slowstoch_rlzs_fMEAN_melt$value, na.rm = T)),
       max(max(fastfluc_rlzs_fMEAN_melt$value, na.rm = T), max(faststoch_rlzs_fMEAN_melt$value, na.rm = T), max(slowfluc_rlzs_fMEAN_melt$value, na.rm = T), max(slowstoch_rlzs_fMEAN_melt$value, na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white", family = "Avenir"),
        legend.position = c(0.25,0.80),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, color = "white", family = "Avenir", face = "bold"),
        legend.text = element_text(size = 13, colour = "white", family = "Avenir", face = "italic"),
        axis.title = element_text(size = 16, color = "white", family = "Avenir"),
        axis.text.x = element_text(size = 11, colour = "gray70", face = "italic"),
        axis.text.y = element_text(size = 14, colour = "gray70", face = "italic"),
        axis.ticks = element_line(colour = "gray30"),
        axis.line = element_line(colour = "gray30"),
        panel.background = element_rect("black"),
        plot.background = element_rect("black"),
        plot.margin = margin(10,15,10,10),
        panel.grid = element_blank())

# Tracking variance of F through time across all realizations in DETERMINISTIC Slow vs. Fast regimes
rlzs_fVAR_deterministic.p <- ggplot() + 
  geom_line(data = slowfluc_rlzs_fVAR_melt, aes(x = variable, y = value, group = realizations, color = "Deterministic Slow"), size = 0.3, alpha = 0.6) +
  geom_line(data = fastfluc_rlzs_fVAR_melt, aes(x = variable, y = value, group = realizations, color = "Deterministic Fast"), size = 0.3, alpha = 0.6) +
  scale_x_discrete(breaks = seq(0, length(slowfluc_rlzs_fVAR_melt$variable), 10)) +
  scale_color_manual(name = "Fluctuation Regimes", values = c("#EB95FF", "#5AFFC3")) +
  guides(color = guide_legend(title.position = 'top', override.aes=list(size=3))) +
  ggtitle(expression(paste("VARIANCE of ", italic(f), " within populations"))) +
  ylab("Var (f)")+
  xlab("Time") +
  ylim(min(min(fastfluc_rlzs_fVAR_melt$value, na.rm = T), min(faststoch_rlzs_fVAR_melt$value, na.rm = T), min(slowfluc_rlzs_fVAR_melt$value, na.rm = T), min(slowstoch_rlzs_fVAR_melt$value, na.rm = T)),
       max(max(fastfluc_rlzs_fVAR_melt$value, na.rm = T), max(faststoch_rlzs_fVAR_melt$value, na.rm = T), max(slowfluc_rlzs_fVAR_melt$value, na.rm = T), max(slowstoch_rlzs_fVAR_melt$value, na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white", family = "Avenir"),
        legend.position = c(0.30,0.80),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, color = "white", family = "Avenir", face = "bold"),
        legend.text = element_text(size = 13, colour = "white", family = "Avenir", face = "italic"),
        axis.title = element_text(size = 16, color = "white", family = "Avenir"),
        axis.text.x = element_text(size = 12, colour = "gray70", face = "italic"),
        axis.text.y = element_text(size = 14, colour = "gray70", face = "italic"),
        axis.ticks = element_line(colour = "gray30"),
        axis.line = element_line(colour = "gray30"),
        panel.background = element_rect("black"),
        plot.background = element_rect("black"),
        plot.margin = margin(10,15,10,10),
        panel.grid = element_blank())

# Tracking variance of F through time across all realizations in STOCHASTIC Slow vs. Fast regimes
rlzs_fVAR_stochastic.p <- ggplot() + 
  geom_line(data = slowstoch_rlzs_fVAR_melt, aes(x = variable, y = value, group = realizations, color = "Stochastic Slow"), size = 0.3, alpha = 0.6) +
  geom_line(data = faststoch_rlzs_fVAR_melt, aes(x = variable, y = value, group = realizations, color = "Stochastic Fast"), size = 0.3, alpha = 0.6) +
  scale_x_discrete(breaks = seq(0, length(slowfluc_rlzs_fMEAN_melt$variable), 10)) +
  scale_color_manual(name = "Fluctuation Regimes", values = c("#BA78F3", "#36F2F9")) +
  guides(color = guide_legend(title.position = 'top', override.aes=list(size=3))) +
  ggtitle(expression(paste("VARIANCE of ", italic(f), " within populations"))) +
  ylab("Var (f)")+
  xlab("Time") +
  ylim(min(min(fastfluc_rlzs_fVAR_melt$value, na.rm = T), min(faststoch_rlzs_fVAR_melt$value, na.rm = T), min(slowfluc_rlzs_fVAR_melt$value, na.rm = T), min(slowstoch_rlzs_fVAR_melt$value, na.rm = T)),
       max(max(fastfluc_rlzs_fVAR_melt$value, na.rm = T), max(faststoch_rlzs_fVAR_melt$value, na.rm = T), max(slowfluc_rlzs_fVAR_melt$value, na.rm = T), max(slowstoch_rlzs_fVAR_melt$value, na.rm = T)))+
  theme(plot.title = element_text(hjust = 0.5, size = 17, color = "white", family = "Avenir"),
        legend.position = c(0.30,0.80),
        legend.key = element_rect(color = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 15, color = "white", family = "Avenir", face = "bold"),
        legend.text = element_text(size = 13, colour = "white", family = "Avenir", face = "italic"),
        axis.title = element_text(size = 16, color = "white", family = "Avenir"),
        axis.text.x = element_text(size = 12, colour = "gray70", face = "italic"),
        axis.text.y = element_text(size = 14, colour = "gray70", face = "italic"),
        axis.ticks = element_line(colour = "gray30"),
        axis.line = element_line(colour = "gray30"),
        panel.background = element_rect("black"),
        plot.background = element_rect("black"),
        plot.margin = margin(10,15,10,10),
        panel.grid = element_blank())


# Showing all MEAN plots together
grid.arrange(rlzs_muMEAN_deterministic.p, rlzs_fMEAN_deterministic.p,
             rlzs_muMEAN_stochastic.p, rlzs_fMEAN_stochastic.p,
             ncol = 2)




#### ANALYSIS of ENDPOINTS of simulations ####

# Get endpoint means and variances across simulation realizations for all four regimes
mumean_fastfluc_end <- sapply(fastfluc_rlzs_muMEAN, function(x) tail(x,1))
mumean_slowfluc_end <- sapply(slowfluc_rlzs_muMEAN, function(x) tail(x,1))
mumean_slowstoch_end <- sapply(slowstoch_rlzs_muMEAN, function(x) tail(x,1))
mumean_faststoch_end <- sapply(faststoch_rlzs_muMEAN, function(x) tail(x,1))

fmean_fastfluc_end <- sapply(fastfluc_rlzs_fMEAN, function(x) tail(x,1))
fmean_slowfluc_end <- sapply(slowfluc_rlzs_fMEAN, function(x) tail(x,1))
fmean_slowstoch_end <- sapply(slowstoch_rlzs_fMEAN, function(x) tail(x,1))
fmean_faststoch_end <- sapply(faststoch_rlzs_fMEAN, function(x) tail(x,1))

muvar_fastfluc_end <- sapply(fastfluc_rlzs_muVAR, function(x) tail(x,1))
muvar_slowfluc_end <- sapply(slowfluc_rlzs_muVAR, function(x) tail(x,1))
muvar_slowstoch_end <- sapply(slowstoch_rlzs_muVAR, function(x) tail(x,1))
muvar_faststoch_end <- sapply(faststoch_rlzs_muVAR, function(x) tail(x,1))

fvar_fastfluc_end <- sapply(fastfluc_rlzs_fVAR, function(x) tail(x,1))
fvar_slowfluc_end <- sapply(slowfluc_rlzs_fVAR, function(x) tail(x,1))
fvar_slowstoch_end <- sapply(slowstoch_rlzs_fVAR, function(x) tail(x,1))
fvar_faststoch_end <- sapply(faststoch_rlzs_fVAR, function(x) tail(x,1))

# Collect endpoint MEANS of mu across all realizations for the four regimes
mumean_end <- melt(as.data.frame(cbind(mumean_fastfluc_end,
                                       mumean_slowfluc_end,
                                       mumean_faststoch_end,
                                       mumean_slowstoch_end
)), na.rm = T)

# Summarize endpoint MEANS of mu across all realizations for the four regimes
mumean_end_summ <- summarySE(mumean_end, measurevar = "value", groupvars = "variable")
mumean_end_summ <- mumean_end_summ[,c(1,3,4)]
colnames(mumean_end_summ) <- c("Env", "mean", "sd")
mumean_end_summ <- melt(mumean_end_summ, id.vars = c("Env", "mean"))

# Plot endpoint MEANS of mu across all realizations for the four regimes
mumean_end.p <- ggplot(mumean_end_summ) +
  geom_point(aes(Env, mean, col = Env, size = 4)) +
  geom_errorbar(aes(x = Env, ymin = mean - value, ymax = mean + value, col = Env), width = 0.1, size = 1) +
  scale_color_manual(values = c("#EB95FF", "#5AFFC3", "#BA78F3", "#36F2F9")) +
  scale_fill_manual(values = c("#EB95FF", "#5AFFC3", "#BA78F3", "#36F2F9")) +
  scale_x_discrete(labels = c("Fast\\nDeterministic", "Slow\\nDeterministic", "Fast\\nStochastic", "Slow\\nStochastic")) +
  ylim(0.02,0.04) +
  ggtitle(expression(paste("Population ", bold("means"), " of ", italic("\\u03BC"), " at ends of SIMULATIONS"))) +
  xlab("Fluctuation Regime") +
  ylab(expression(paste("Mean(", italic("\\u03BC"), ")"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Avenir",  size = 14, face = "bold"),
        axis.title.y = element_text( vjust = 0.5, family = "Avenir"),
        legend.position = "none",
        panel.grid = element_line(color = "gray90"),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(family = "Avenir", size = 13, face = "bold"),
        axis.text = element_text(family = "Avenir", size = 12)) 


# Collect endpoint MEANS of f across all realizations for the four regimes
fmean_end <- melt(as.data.frame(cbind(fmean_fastfluc_end,
                                      fmean_slowfluc_end,
                                      fmean_faststoch_end,
                                      fmean_slowstoch_end
)), na.rm = T)

# Summarize endpoint MEANS of f across all realizations for the four regimes
fmean_end_summ <- summarySE(fmean_end, measurevar = "value", groupvars = "variable")
fmean_end_summ <- fmean_end_summ[,c(1,3,4)]
colnames(fmean_end_summ) <- c("Env", "mean", "sd")
fmean_end_summ <- melt(fmean_end_summ, id.vars = c("Env", "mean"))

# Plot endpoint MEANS of f across all realizations for the four regimes
fmean_end.p <- ggplot(fmean_end_summ) +
  geom_point(aes(Env, mean, col = Env), size = 4) +
  geom_errorbar(aes(x = Env, ymin = mean - value, ymax = mean + value, col = Env), width = 0.2, size = 1) +
  scale_color_manual(values = c("#EB95FF", "#5AFFC3", "#BA78F3", "#36F2F9")) +
  scale_fill_manual(values = c("#EB95FF", "#5AFFC3", "#BA78F3", "#36F2F9")) +
  scale_x_discrete(labels = c("Fast\\nDeterministic", "Slow\\nDeterministic", "Fast\\nStochastic", "Slow\\nStochastic")) +
  ylim(8,13) +
  ggtitle(expression(paste("Population ", bold("means"), " of ", italic(f), " at ends of SIMULATIONS"))) +
  xlab("Fluctuation Regime") +
  ylab("Mean(f)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Avenir",  size = 14, face = "bold"),
        axis.title.y = element_text( vjust = 0.5, family = "Avenir"),
        legend.position = "none",
        panel.grid = element_line(color = "gray90"),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(family = "Avenir", size = 13, face = "bold"),
        axis.text = element_text(family = "Avenir", size = 12)) 


# Collect endpoint VARIANCES of mu across all realizations for the four regimes
muvar_end <- melt(as.data.frame(cbind(muvar_fastfluc_end,
                                      muvar_slowfluc_end,
                                      muvar_faststoch_end,
                                      muvar_slowstoch_end
)), na.rm = T)

# Summarize endpoint VARIANCES of mu across all realizations for the four regimes
muvar_end_summ <- summarySE(muvar_end, measurevar = "value", groupvars = "variable")
muvar_end_summ <- muvar_end_summ[,c(1,3,4)]
colnames(muvar_end_summ) <- c("Env", "mean", "sd")
muvar_end_summ <- melt(muvar_end_summ, id.vars = c("Env", "mean"))

# Plot endpoint VARIANCES of mu across all realizations for the four regimes
muvar_end.p <- ggplot(muvar_end_summ) +
  geom_errorbar(aes(x = Env, ymin = mean - value, ymax = mean + value, col = Env), width = 0.1, size = 1.2) +
  geom_point(aes(Env, mean, col = Env, size = 4), shape = 22, size = 4, stroke = 2, fill = "black") +
  scale_color_manual(values = c("#EB95FF", "#5AFFC3", "#BA78F3", "#36F2F9")) +
  scale_fill_manual(values = c("#EB95FF", "#5AFFC3", "#BA78F3", "#36F2F9")) +
  scale_x_discrete(labels = c("Fast\\nDeterministic", "Slow\\nDeterministic", "Fast\\nStochastic", "Slow\\nStochastic")) +
  #ylim(1e-05, 5e-05) +
  ggtitle(expression(atop("SIMULATION", paste("Population ", bold("variances"), " of ", italic("\\u03BC"))))) +
  xlab("Fluctuation Regime") +
  ylab(expression(paste("Var(", mu, ")"))) +
  theme_minimal() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Avenir",  size = 14, face = "bold"),
        axis.title = element_text(family = "Avenir", size = 13, face = "bold"),
        axis.text = element_text(family = "Avenir", size = 12),
        #axis.text.x = element_text(angle = 45, vjust = 0.8),
        legend.position = "none",
        panel.background = element_rect(color = "black", size = 1, fill = "gray40"),
        panel.grid = element_line(color = "gray70"),
        panel.grid.major.x = element_blank()) 


# Collect endpoint VARIANCES of f across all realizations for the four regimes
fvar_end <- melt(as.data.frame(cbind(fvar_fastfluc_end,
                                     fvar_slowfluc_end,
                                     fvar_faststoch_end,
                                     fvar_slowstoch_end
)), na.rm = T)

# Summarize endpoint VARIANCES of f across all realizations for the four regimes
fvar_end_summ <- summarySE(fvar_end, measurevar = "value", groupvars = "variable")
fvar_end_summ <- fvar_end_summ[,c(1,3,4)]
colnames(fvar_end_summ) <- c("Env", "mean", "sd")
fvar_end_summ <- melt(fvar_end_summ, id.vars = c("Env", "mean"))

# Plot endpoint VARIANCES of f across all realizations for the four regimes
fvar_end.p <- ggplot(fvar_end_summ) +
  geom_errorbar(aes(x = Env, ymin = mean - value, ymax = mean + value, col = Env), width = 0.1, size = 1.2) +
  geom_point(aes(Env, mean, col = Env, size = 4), shape = 22, size = 4, stroke = 2, fill = "black") +
  scale_color_manual(values = c("#EB95FF", "#5AFFC3", "#BA78F3", "#36F2F9")) +
  scale_fill_manual(values = c("#EB95FF", "#5AFFC3", "#BA78F3", "#36F2F9")) +
  scale_x_discrete(labels = c("Fast\\nDeterministic", "Slow\\nDeterministic", "Fast\\nStochastic", "Slow\\nStochastic")) +
  #ylim(0,3.5)+
  ggtitle(expression(atop("SIMULATION", paste("Population ", bold("variances"), " of ", italic("f"))))) +
  xlab("Fluctuation Regime") +
  ylab("Var(f)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, family = "Avenir",  size = 14, face = "bold"),
        axis.title = element_text(family = "Avenir", size = 13, face = "bold"),
        axis.text = element_text(family = "Avenir", size = 12),
        #axis.text.x = element_text(angle = 45, vjust = 0.8),
        legend.position = "none",
        panel.background = element_rect(color = "black", size = 1, fill = "gray40"),
        panel.grid = element_line(color = "gray70"),
        panel.grid.major.x = element_blank()) 


# Pairwise combinations of Two-way t-tests between endpoint means or variances

t.test(log(fvar_slowfluc_end), log(fvar_fastfluc_end), var.equal = F, na.action = na.omit)

t.test(log(fvar_slowstoch_end), log(fvar_faststoch_end), var.equal = F, na.action = na.omit)

t.test(log(fvar_slowfluc_end), log(fvar_slowstoch_end), var.equal = F, na.action = na.omit)

t.test(log(fvar_fastfluc_end), log(fvar_faststoch_end), var.equal = F, na.action = na.omit)


t.test(log(muvar_slowfluc_end), log(muvar_fastfluc_end), var.equal = F, na.action = na.omit)

t.test(log(muvar_slowstoch_end), log(muvar_faststoch_end), var.equal = F, na.action = na.omit)

t.test(log(muvar_slowfluc_end), log(muvar_slowstoch_end), var.equal = F, na.action = na.omit)

t.test(log(muvar_fastfluc_end), log(muvar_faststoch_end), var.equal = F, na.action = na.omit)



### Visualizing demographic dynamics ####

# Juvenile, Adult, and Total Population size plots of the four regimes from one example simulation realization
slow_juv <- sapply(sim_slowfluc.ls, function(x) length(which(x[,3] < sexually_mature)))
slow_adult <- sapply(sim_slowfluc.ls, function(x) length(which(x[,3] > sexually_mature)))
slow_totalpopsize <- sapply(sim_slowfluc.ls, function(x) nrow(x))

slowstoch_juv <- sapply(sim_slowstoch.ls, function(x) length(which(x[,3] < sexually_mature)))
slowstoch_adult <- sapply(sim_slowstoch.ls, function(x) length(which(x[,3] > sexually_mature)))
slowstoch_totalpopsize <- sapply(sim_slowstoch.ls, function(x) nrow(x))

fast_juv <- sapply(sim_fastfluc.ls, function(x) length(which(x[,3] < 1.0)))
fast_adult <- sapply(sim_fastfluc.ls, function(x) length(which(x[,3] > 1.0)))
fast_totalpopsize <- sapply(sim_fastfluc.ls, function(x) nrow(x))

faststoch_juv <- sapply(sim_faststoch.ls, function(x) length(which(x[,3] < sexually_mature)))
faststoch_adult <- sapply(sim_faststoch.ls, function(x) length(which(x[,3] > sexually_mature)))
faststoch_totalpopsize <- sapply(sim_faststoch.ls, function(x) nrow(x))

  
# Deterministic Slow regime Juvenile vs. Adult number trajectory
ggplot() + geom_line(aes(x = 0:num_steps, y = log(slow_juv), col = "Juvenile")) +
  geom_line(aes(x = 0:num_steps, y = log(slow_adult), col = "Adult")) +
  scale_color_manual(name = "Stage",
                     values = c("blue", "red")) +
  scale_x_continuous(limits = c(0,150)) +
  ggtitle("Stage-Structure Dynamics in Slow Environment") +
  xlab("Time Steps") +
  ylab("Log (count)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 
  
# Stochastic Slow regime Juvenile vs. Adult number trajectory
ggplot() + geom_line(aes(x = 0:num_steps, y = log(slowstoch_juv), col = "Juvnile")) +
  geom_line(aes(x = 0:num_steps, y = log(slowstoch_adult), col = "Adult")) +
  scale_color_manual(name = "Stage",
                     values = c("blue", "red")) +
  scale_x_continuous(limits = c(0,150)) +
  ggtitle("Stage-Structure Dynamics in Slow Stochastic Environment") +
  xlab("Time Steps") +
  ylab("Log (count)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Deterministic Fast regime Juvenile vs. Adult number trajectory
ggplot() + geom_line(aes(x = 0:num_steps, y = log(fast_juv), col = "Juvenile")) +
  geom_line(aes(x = 0:num_steps, y = log(fast_adult), col = "Adult")) +
  scale_color_manual(name = "Stage",
                     values = c("blue", "red")) +
  scale_x_continuous(limits = c(0,150)) +
  ggtitle("Stage-Structure Dynamics in Fast Environment") +
  xlab("Time Steps") +
  ylab("Log (count)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Stochastic Fast regime Juvenile vs. Adult number trajectory
ggplot() + geom_line(aes(x = 0:num_steps, y = log(faststoch_juv), col = "Juvenile")) +
  geom_line(aes(x = 0:num_steps, y = log(faststoch_adult), col = "Adult")) +
  scale_color_manual(name = "Stage",
                     values = c("blue", "red")) +
  scale_x_continuous(limits = c(0,150)) +
  ggtitle("Stage-Structure Dynamics in Fast Stochastic Environment") +
  xlab("Time Steps") +
  ylab("Log (count)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Deterministic Slow and Fast regimes total population trajectory
ggplot() + geom_line(aes(x = 0:num_steps, y = log(slow_totalpopsize), col = "Slow"))+
  geom_line(aes(x = 0:num_steps, y = log(fast_totalpopsize), col = "Fast")) +
  scale_color_manual(name = "Environment",
                     values = c("black", "gray")) +
  scale_x_continuous(limits = c(0,150)) +
  ggtitle("Total Population Size") +
  xlab("Time Steps") +
  ylab("Log (count)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Stochastic Slow and Fast regimes total population trajectory
ggplot() + geom_line(aes(x = 0:num_steps, y = log(slowstoch_totalpopsize), col = "Slow Stochastic"))+
  geom_line(aes(x = 0:num_steps, y = log(faststoch_totalpopsize), col = "Fast Stochastic")) +
  scale_color_manual(name = "Environment",
                     values = c("black", "gray")) +
  scale_x_continuous(limits = c(0,150)) +
  ggtitle("Total Population Size") +
  xlab("Time Steps") +
  ylab("Log (count)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 
