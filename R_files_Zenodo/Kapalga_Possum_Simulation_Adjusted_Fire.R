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

possum_juvenile_density <- 0.3744
possum_adult_density <- 0.4056

total_juveniles <- round(possum_juvenile_density * sum(HS[!is.na(HS)]), 0)
total_adults <- round(possum_adult_density * sum(HS[!is.na(HS)]), 0)

Possum_init_densities <- stack(HS, HS)
names(Possum_init_densities) <- c("Juvenile", "Adult")

Possum_init_densities[!is.na(Possum_init_densities)] <- 0

Possum_init_densities[[1]][sampleRandom(Possum_init_densities[[1]], total_juveniles, cells = TRUE)[,1]] <- 1
Possum_init_densities[[2]][sampleRandom(Possum_init_densities[[2]], total_adults, cells = TRUE)[,1]] <- 1

#########################
### Transition matrix ###
#########################

# Female-only transition matrix adapted from Tony's thesis

possum_trans_mat <- matrix(c(0.55,0.46,
                             0.25,0.79),
                           nrow = 2, ncol = 2, byrow = TRUE)
colnames(possum_trans_mat) <- rownames(possum_trans_mat) <- c('Juvenile','Adult')


# Stochasticity matrix

possum_stoch_mat <- matrix(c(0.09,0.26,
                             0.11,0.08),
                           nrow = 2, ncol = 2, byrow = TRUE)
colnames(possum_stoch_mat) <- rownames(possum_stoch_mat) <- c('Juvenile','Adult')


# Initial growth rate is slightly positive and lower than Tony's thesis

(Rmax_pop <- abs(eigen(possum_trans_mat)$values[1])) #1.029722

# Carrying capacity

max_ind <- 2.46

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

ca_dispersal <- cellular_automata_dispersal(max_cells = c(3,3),
                                            dispersal_proportion = density_dependence_dispersing(maximum_proportions = c(1,0.1)))

# Specify pop dynamics and landscape objects:

poss_landscape <- landscape(population = Possum_init_densities,
                            suitability = HS,
                            carrying_capacity = k_function,
                            "Early_fire" = Early_fire,
                            "Late_fire" = Late_fire)

cap_population <- ceiling_density(stages = c(1, 2))

pd <- population_dynamics(change = growth(possum_trans_mat,
                                          transition_function = modified_transition_custom(fire_stack_early = "Early_fire",
                                                                                           fire_stack_late = "Late_fire",
                                                                                           early_surv_mult = 0.915,
                                                                                           late_surv_mult = 0.85,
                                                                                           early_tran_mult = 0.915,
                                                                                           late_tran_mult = 0.85,
                                                                                           early_recr_mult = 0.73,
                                                                                           late_recr_mult = 0.51),
                                          global_stochasticity = possum_stoch_mat),
                          dispersal = ca_dispersal,
                          density_dependence = cap_population)

# Run simulation

plan(multisession, workers = 40)

system.time(
  possum_AF_pop_data <- foreach (i = 1:25, .combine = rbind) %do% {
    set.seed(333)
    Poss2 <- simulation(landscape = poss_landscape,
                        population_dynamics = pd,
                        habitat_dynamics = list(habitat_suitability_mod(fire_stack_early = "Early_fire",
                                                                        fire_stack_late = "Late_fire",
                                                                        reduction_early = 0.9,
                                                                        reduction_late = 0.75,
                                                                        noise = 0.1)),
                        demo_stochasticity = "none",
                        timesteps = 126,
                        replicates = 40,
                        future.globals = list(modified_transition_custom = modified_transition_custom,
                                              possum_trans_mat = possum_trans_mat,
                                              possum_stoch_mat = possum_stoch_mat,
                                              max_ind = max_ind))
    
    pop_data <- return_pop_data(Poss2)
    rm(Poss2)
    gc()
    pop_data
  }
)

save(possum_AF_pop_data, file = "possum_AF_pop_data.RData")

# Plot possum population trajectory of each life-stage
plot(Poss)

# Plot total possum population trajectory
plot(Poss2, stage = 0, ylim = c(0,50000))

# Plot possum abundance by timestep
plot(Poss[1], type = "raster", stage = 0, timesteps = 1:126, animate = TRUE)
