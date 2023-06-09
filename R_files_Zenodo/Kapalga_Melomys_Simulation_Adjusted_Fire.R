# Load packages
library(steps)
library(raster)
library(viridis)
library(future)
library(foreach)

# Source functions:
source("custom_functions.R")

# Load in spatial layers:
load(file = "spatial_rasters")

# Create initial populations:

melomys_juvenile_density <- 0.2365
melomys_adult_density <- 0.1935

total_juveniles <- round(melomys_juvenile_density * sum(HS[!is.na(HS)]), 0)
total_adults <- round(melomys_adult_density * sum(HS[!is.na(HS)]), 0)

Melomys_init_densities <- stack(HS, HS)
names(Melomys_init_densities) <- c("Juvenile", "Adult")

Melomys_init_densities[!is.na(Melomys_init_densities)] <- 0

Melomys_init_densities[[1]][sampleRandom(Melomys_init_densities[[1]], total_juveniles, cells = TRUE)[,1]] <- 1
Melomys_init_densities[[2]][sampleRandom(Melomys_init_densities[[2]], total_adults, cells = TRUE)[,1]] <- 1

#########################
### Transition matrix ###
#########################

# Female-only transition matrix adapted from Tony's thesis:

melomys_trans_mat <- matrix(c(0.36,0.78,
                              0.20,0.77),
                            nrow = 2, ncol = 2, byrow = TRUE)
colnames(melomys_trans_mat) <- rownames(melomys_trans_mat) <- c('Juvenile','Adult')


# Stochasticity matrix

melomys_stoch_mat <- matrix(c(0.05,0.23,
                              0.04,0.05),
                            nrow = 2, ncol = 2, byrow = TRUE)
colnames(melomys_stoch_mat) <- rownames(melomys_stoch_mat) <- c('Juvenile','Adult')


# Initial growth rate is slightly positive and lower than Tony's thesis

(Rmax_pop <- abs(eigen(melomys_trans_mat)$values[1])) #1.01

# Carrying capacity

max_ind <- 1.08

# x <- seq(0,1,.1)
# 
# plot(x,
#      max_ind / (1 + exp(-(x - 0.5) / 0.05)),
#      type = 'l',
#      xlab = "Habitat Quality",
#      ylab = "Max Individuals")

###################
### Steps setup ###
###################

ca_dispersal <- cellular_automata_dispersal(max_cells = c(10,10),
                                            dispersal_proportion = density_dependence_dispersing(maximum_proportions = c(1,0.1)))

# Specify pop dynamics and landscape objects:

melomys_landscape <- landscape(population = Melomys_init_densities,
                               suitability = HS,
                               carrying_capacity = k_function,
                               "Early_fire" = Early_fire,
                               "Late_fire" = Late_fire)


cap_population <- ceiling_density(stages = c(1, 2))

pd <- population_dynamics(change = growth(melomys_trans_mat,
                                          transition_function = modified_transition_custom(fire_stack_early = "Early_fire",
                                                                                           fire_stack_late = "Late_fire",
                                                                                           early_surv_mult = 0.79,
                                                                                           late_surv_mult = 0.825,
                                                                                           early_tran_mult = 0.79,
                                                                                           late_tran_mult = 0.825,
                                                                                           early_recr_mult = 1.45,
                                                                                           late_recr_mult = 1.13),
                                          global_stochasticity = melomys_stoch_mat),
                          dispersal = ca_dispersal,
                          density_dependence = cap_population)

# Run simulation

plan(multisession, workers = 40)

system.time(
  melomys_AF_pop_data <- foreach (i = 1:25, .combine = rbind) %do% {
    set.seed(333)
    Melomys2 <- simulation(landscape = melomys_landscape,
                           population_dynamics = pd,
                           habitat_dynamics = list(habitat_suitability_mod(fire_stack_early = "Early_fire",
                                                                           fire_stack_late = "Late_fire",
                                                                           reduction_early = 0.7,
                                                                           reduction_late = 0.4,
                                                                           noise = 0.1)),
                           demo_stochasticity = "none",
                           timesteps = 126,
                           replicates = 40,
                           future.globals = list(modified_transition_custom = modified_transition_custom,
                                                 melomys_trans_mat = melomys_trans_mat,
                                                 melomys_stoch_mat = melomys_stoch_mat,
                                                 max_ind = max_ind))
    
    pop_data <- return_pop_data(Melomys2)
    rm(Melomys2)
    gc()
    pop_data
  }
)

save(melomys_AF_pop_data, file = "melomys_AF_pop_data.RData")

# Plot possum population trajectory of each life-stage
plot(Melomys)

# Plot total possum population trajectory
plot(Melomys2, stage = 0, ylim = c(0,15000))

# Plot possum abundance by timestep
plot(Melomys[1], type = "raster", stage = 0, timesteps = 1:126, animate = TRUE)
