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

# Change late season fires to early season: note this overwrites the layer stacks
for (i in 1:nlayers(Early_fire)) {
  idx <- which(Late_fire[[i]][] == 1)
  if (length(idx) != 0) {
    Early_fire[[i-2]][idx] <- 1
    Late_fire[[i]][idx] <- 0
  }
}

# Create initial populations:

bandicoot_juvenile_density <- 0.3713
bandicoot_adult_density <- 0.4187

total_juveniles <- round(bandicoot_juvenile_density * sum(HS[!is.na(HS)]), 0)
total_adults <- round(bandicoot_adult_density * sum(HS[!is.na(HS)]), 0)

Bandicoot_init_densities <- stack(HS, HS)
names(Bandicoot_init_densities) <- c("Juvenile", "Adult")

Bandicoot_init_densities[!is.na(Bandicoot_init_densities)] <- 0

Bandicoot_init_densities[[1]][sampleRandom(Bandicoot_init_densities[[1]], total_juveniles, cells = TRUE)[,1]] <- 1
Bandicoot_init_densities[[2]][sampleRandom(Bandicoot_init_densities[[2]], total_adults, cells = TRUE)[,1]] <- 1

#########################
### Transition matrix ###
#########################

# Female-only transition matrix adapted from Tony's thesis

bandicoot_trans_mat <- matrix(c(0.55,0.45,
                                0.36,0.75),
                              nrow = 2, ncol = 2, byrow = TRUE)
colnames(bandicoot_trans_mat) <- rownames(bandicoot_trans_mat) <- c('Juvenile','Adult')


# Stochasticity matrix

bandicoot_stoch_mat <- matrix(c(0.07,0.12,
                                0.10,0.11),
                              nrow = 2, ncol = 2, byrow = TRUE)
colnames(bandicoot_stoch_mat) <- rownames(bandicoot_stoch_mat) <- c('Juvenile','Adult')


# Initial growth rate is slightly positive and lower than Tony's thesis

(Rmax_pop <- abs(eigen(bandicoot_trans_mat)$values[1])) #1.064

# Carrying capacity

max_ind <- 2.0

# Plot relationship of carrying capacity to habitat suitability,

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

ca_dispersal <- cellular_automata_dispersal(max_cells = c(50,50),
                                            dispersal_proportion = density_dependence_dispersing(maximum_proportions = c(1,0.1)))

# Specify pop dynamics and landscape objects:

coot_landscape <- landscape(population = Bandicoot_init_densities,
                            suitability = HS,
                            carrying_capacity = k_function,
                            "Early_fire" = Early_fire,
                            "Late_fire" = Late_fire)

cap_population <- ceiling_density(stages = c(1, 2))

pd <- population_dynamics(change = growth(bandicoot_trans_mat,
                                          transition_function = modified_transition_custom(fire_stack_early = "Early_fire",
                                                                                           fire_stack_late = "Late_fire",
                                                                                           early_surv_mult = 0.93,
                                                                                           late_surv_mult = 0.8,
                                                                                           early_tran_mult = 0.93,
                                                                                           late_tran_mult = 0.8,
                                                                                           early_recr_mult = 1.0,
                                                                                           late_recr_mult = 0.605),
                                          global_stochasticity = bandicoot_stoch_mat),
                          dispersal = ca_dispersal,
                          density_dependence = cap_population)

# Run simulation

plan(multisession, workers = 40)

system.time(
  bandicoot_AF_AEF_pop_data <- foreach (i = 1:25, .combine = rbind) %do% {
    set.seed(333)
    Coot4 <- simulation(landscape = coot_landscape,
                        population_dynamics = pd,
                        habitat_dynamics = list(habitat_suitability_mod(fire_stack_early = "Early_fire",
                                                                        fire_stack_late = "Late_fire",
                                                                        reduction_early = 0.5,
                                                                        reduction_late = 0.3,
                                                                        noise = 0.1)),
                        demo_stochasticity = "none",
                        timesteps = 126,
                        replicates = 40,
                        future.globals = list(modified_transition_custom = modified_transition_custom,
                                              bandicoot_trans_mat = bandicoot_trans_mat,
                                              bandicoot_stoch_mat = bandicoot_stoch_mat,
                                              max_ind = max_ind))
    
    pop_data <- return_pop_data(Coot4)
    rm(Coot4)
    gc()
    pop_data
  }
)

save(bandicoot_AF_AEF_pop_data, file = "bandicoot_AF_AEF_pop_data.RData")

# Plot bandicoot population trajectory of each life-stage
plot(Coot)

# Plot total bandicoot population trajectory
plot(Coot4, stage = 0, ylim = c(0,35000))

# Plot bandicoot abundance by timestep
plot(Coot[1], type = "raster", stage = 0, timesteps = 1:126, animate = TRUE)
