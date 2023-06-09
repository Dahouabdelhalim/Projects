library(deSolve) # For solving differential equations
library(tidyverse) # For efficient data manipulation and plotting
library(rcartocolor) # Colorblind-friendly color palette
theme_set(theme_bw()) # Set basic plotting theme
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Evaluate the right-hand side of the model equations
# Input:
# - time: Moment to evaluate right-hand side at
# - p: Vector of state variables (here: site occupancies)
# - pars: A list of model parameters, with the following entries:
#     $c: Vector of colonization rates
#     $e: Vector of mortality rates (or single mortality, if equal for all species)
#     $D: The extent of the disturbance (i.e., fraction of density removed)
#     $period: Disturbance period
#     $H: Competitive matrix; H[i,j] = prob. that individual of species i displaces j
# Output:
# - The vector of right-hand sides, coerced into a list
modelEqs <- function(time, p, pars) {
  recruitment <- pars$c * p * (1 - sum(p)) # c_i * p_i * (1 - sum_j p_j)
  # Mortality: e_i * p_i plus forcing term, -log(1 - D) / period, which
  # ensures that p_i drops to (1 - D) * p_i over one period
  mortality <- pars$e * p - (p * log(1 - pars$D) / pars$period)
  # Displacement term:
  displacement <- as.numeric(diag(pars$c * p) %*% (pars$H %*% p) -
                               diag(p) %*% (t(pars$H) %*% (pars$c * p)))
  return(list(recruitment - mortality + displacement))
}


# Integrate the model equations and return solution in a tidy table
# Input:
# - pars: A list with the following entries:
#     $n: Number of species
#     $D: Extent of disturbance (fraction of density removed)
#     $period: Disturbance period
#     $c: Vector of colonization rates
#     $e: Vector of mortality rates (or single mortality, if equal for all species)
#     $H: Competitive matrix; H[i,j] = prob. that individual of species i displaces j
#     $tmax: Integration time
# Output:
# - A tidy data frame with three columns: time, species, and p (site occupancy)
runModel <- function(pars) {
  ic <- rep(1 / pars$n, pars$n) # Initial conditions
  tseq <- seq(0, pars$tmax, by = pars$tstep) # Sampling points in time
  ode(func = modelEqs, y = ic, times = tseq, parms = pars, method = "bdf") %>%
    as.data.frame() %>% # Convert solution to data frame (needed for next step)
    as_tibble() %>% # Convert solution to tibble
    pivot_longer(cols = 2:ncol(.), names_to = "species", values_to = "p") # Tidy up data
}


# Obtain the matrix H for n species
# Input:
# - n: The number of species (so H will be an n-by-n matrix)
# - rseed: Random seed (for reproducibility of randomly generated matrix entries)
# - u: Limit for upper-triangular entries (which go between 1-u and 1)
# - l: Limit for lower-triangular entries (which go between 0 and l)
# Output:
# - An n-by-n matrix whose (i,j)th entry is the probability that an individual of
#   species i displaces an individual of species j in competition for a site
Hmat <- function(n, rseed, u, l) {
  set.seed(rseed) # Set random seed
  H <- matrix(0, n, n) # Initialize H matrix to 0
  H[upper.tri(H)] <- runif(choose(n, 2), 1-u, 1) # Upper triangular entries from [1-u, 1]
  H[lower.tri(H)] <- runif(choose(n, 2), 0, l) # Lower triangular entries from [0, l]
  return(H)
}


# Uniformly drawn colonization rates, sorted in increasing order
# Input:
# - rseed: Random seed (for reproducibility)
# - n: The number of species (so the function returns a vector with n entries)
# - clower: Lower limit for randomly drawn colonization rates
# - cupper: Upper limit for randomly drawn colonization rates
# Output:
# - A vector with n entries, sorted in increasing order but otherwise randomly and
#   independently drawn from the interval [clower, cupper].
randomColon <- function(rseed, n, clower, cupper) {
  set.seed(rseed) # Set random seed
  list(sort(runif(n, clower, cupper))) # Return sorted random colonization rates
}


# Create table of a factorial numerical experiment; each row contains the parameters
# for one model run
# Input:
# - n: Number of species
# - D: Extent of disturbance (fraction of density removed)
# - period: Disturbance period
# - tmax: Number of time units to integrate equations for
# - tstep: Time step for ODE integration output
# - clower: Lower limit for colonization rates
# - cupper: Upper limit for colonization rates
# - e: Vector of mortality rates (or single mortality, if equal for all species)
# - rseed: Random seed (for reproducibility of randomly generated values)
# - randc: TRUE if colonization rates are randomized; FALSE if they are evenly spaced
# - u: Limit for upper-triangular entries of H (which will go between 1-u and 1)
# - l: Limit for lower-triangular entries of H (which will go between 0 and l)
# Output:
# - A tibble with the parameters of a single numerical scenario in each row
paramTable <- function(n = 3:6, D = seq(0, 0.5, by = 0.005), period = 1, tmax = 20000,
                       tstep = 20, clower = 0.45, cupper = 0.8, e = 0.2, rseed = 57891,
                       randc = FALSE, u = 0, l = 0) {
  expand_grid(n = n, D = D, period = period, tmax = tmax, u = u, l = l) %>%
    mutate(rseed = rseed + 100 * (n + l) + 10 * u) %>% # Random seeds
    rowwise() %>% # For each row of this table:
    mutate(c = if_else(randc, # If colonization rates are randomized:
                       randomColon(rseed, n, clower, cupper), # Random rates
                       list(seq(clower, cupper, l = n))), # Else: evenly spaced rates
           e = list(rep(e, n)), # Extinction rates
           H = list(Hmat(n, rseed = rseed, u = u, l = l))) %>% # H matrix
    ungroup() %>% # Drop the grouping information introduced by rowwise() above
    mutate(pars = mapply(list, n = n, D = D, period = period, tmax = tmax, tstep = tstep,
                         c = c, e = e, H = H,
                         SIMPLIFY = FALSE)) # Merge all parameters into single list
}



# Figure 1
paramTable(rseed = 0) %>% # Assemble table of parameter combinations
  bind_rows(paramTable(clower = 0.25, cupper = 1, rseed = 1e4)) %>% # Add to that table
  mutate(result = map(pars, runModel)) %>% # Integrate model for each row of data
  unnest(cols = c(result)) %>% # Unpack solution into the main table
  filter(time >= tmax * 3 / 4) %>% # Drop transients
  group_by(n, D, period, tmax, rseed, species) %>%
  summarise(meanp = mean(p), .groups = "drop_last") %>% # Mean steady-state occupancies
  filter(meanp > 1e-15) %>% # Drop species that are effectively at zero density
  mutate(meanp = meanp / sum(meanp)) %>% # Express occupancies as fraction, to calculate:
  summarise(`Shannon index` = -sum(meanp * log(meanp)), # Shannon diversity index
            `Inverse Simpson index` = 1 / sum(meanp^2), # Inverse Simpson index
            `Species richness` = n(), .groups="drop") %>% # Raw species richness
  na_if(Inf) %>%
  pivot_longer(cols = c("Shannon index", "Inverse Simpson index", "Species richness"),
               names_to = "index", values_to = "diversity") %>% # Tidy up table
  mutate(n = as.character(n), # Manipulate data to show better on plot
         index = fct_relevel(index, "Shannon index", "Inverse Simpson index",
                             "Species richness"),
         `Spread of colonization rates` = if_else(rseed < 1000, "Low", "High"),
         `Spread of colonization rates` = fct_relevel(`Spread of colonization rates`,
                                               "Low", "High")) %>%
  ggplot() + # Create plot
  aes(x = D, y = diversity, colour = `Spread of colonization rates`) +
  geom_line() +
  facet_grid(index ~ n, labeller = label_bquote(cols = italic(n) == .(n)),
             scales = "free_y") +
  scale_x_continuous(name = expression(paste("Disturbance extent (",italic(D),")"))) +
  scale_y_continuous(name = "Species diversity") +
  scale_colour_manual(values = c("steelblue", "goldenrod"))


# Figure 2
paramTable(rseed = 0) %>%
  bind_rows(paramTable(clower = 0.25, cupper = 1, rseed = 1e4)) %>%
  mutate(result = map(pars, runModel)) %>%
  unnest(cols = c(result)) %>%
  filter(time >= tmax * 3 / 4, rseed < 1000) %>%
  group_by(n, D, period, tmax, species) %>%
  summarise(meanp = mean(p), .groups = "drop_last") %>%
  filter(meanp > 1e-15) %>%
  mutate(meanp = meanp / sum(meanp)) %>%
  ungroup() %>%
  mutate(Species = fct_reorder(species, as.numeric(species))) %>%
  ggplot() +
  aes(x = D, y = meanp, colour = Species) +
  geom_line() +
  facet_wrap(~ n, labeller = label_bquote(cols = italic(n) == .(n))) +
  scale_x_continuous(name = expression(paste("Disturbance extent (",italic(D),")"))) +
  scale_y_continuous(name = "Relative abundances at equilibrium") +
  scale_colour_manual(values = carto_pal(name = "Safe"))


# Figure 3
paramTable(rseed = 0, u = 0.25, l = 0.25) %>%
  bind_rows(paramTable(randc = TRUE, rseed = 1e5)) %>%
  mutate(result = map(pars, runModel)) %>%
  unnest(cols = c(result)) %>%
  filter(time >= tmax * 3 / 4) %>%
  group_by(n, D, period, tmax, u, species) %>%
  summarise(meanp = mean(p), .groups = "drop_last") %>%
  filter(meanp > 1e-15) %>%
  mutate(meanp = meanp / sum(meanp)) %>%
  summarise(`Shannon index` = -sum(meanp * log(meanp)),
            `Inverse Simpson index` = 1 / sum(meanp^2),
            `Species richness` = n(), .groups="drop") %>%
  na_if(Inf) %>%
  pivot_longer(cols = c("Shannon index", "Inverse Simpson index", "Species richness"),
               names_to = "index", values_to = "diversity") %>%
  mutate(n = as.character(n),
         index = fct_relevel(index, "Shannon index", "Inverse Simpson index",
                             "Species richness"),
         Model = if_else(u == 0, "Perturbed colonization rates", "Perturbed hierarchy")) %>%
  ggplot() +
  aes(x = D, y = diversity, colour = Model) +
  geom_line() +
  facet_grid(index ~ n, labeller = label_bquote(cols = italic(n) == .(n)),
             scales = "free_y") +
  scale_x_continuous(name = expression(paste("Disturbance extent (",italic(D),")"))) +
  scale_y_continuous(name = "Species diversity") +
  scale_colour_manual(values = c("steelblue", "goldenrod"))


# Figure 4
paramTable(n = 25, tstep = 100, l = 0.1, u = 0.1, rseed = 98765, randc = TRUE) %>%
  bind_rows(paramTable(n = 25, clower = 0.25, tstep = 50, cupper = 1,
                       rseed = 200, randc = TRUE, l = 0.1, u = 0.1)) %>%
  mutate(result = map(pars, runModel)) %>%
  unnest(cols = c(result)) %>%
  filter(time >= tmax * 3 / 4) %>%
  group_by(n, D, period, tmax, rseed, species) %>%
  summarise(meanp = mean(p), .groups = "drop_last") %>%
  filter(meanp > 1e-15) %>%
  mutate(meanp = meanp / sum(meanp)) %>%
  summarise(`Shannon index` = -sum(meanp * log(meanp)),
            `Inverse Simpson index` = 1 / sum(meanp^2),
            `Species richness` = n(), .groups="drop") %>%
  na_if(Inf) %>%
  pivot_longer(cols = c("Shannon index", "Inverse Simpson index", "Species richness"),
               names_to = "index", values_to = "diversity") %>%
  mutate(n = as.character(n),
         index = fct_relevel(index, "Shannon index", "Inverse Simpson index",
                             "Species richness"),
         `Spread of colonization rates` = if_else(rseed < 1e4, "High", "Low"),
         `Spread of colonization rates` = fct_relevel(`Spread of colonization rates`,
                                               "Low", "High")) %>%
  ggplot() +
  aes(x = D, y = diversity, colour = `Spread of colonization rates`) +
  geom_line() +
  facet_grid(index ~ n, labeller = label_bquote(cols = italic(n) == .(n)),
             scales = "free_y") +
  scale_x_continuous(name = expression(paste("Disturbance extent (",italic(D),")"))) +
  scale_y_continuous(name = "Species diversity") +
  scale_colour_manual(values = c("steelblue", "goldenrod"))
