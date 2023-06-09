# Code to reproduce figures in Miller and Allesina (2022)
# Habitat heterogeneity, environmental feedbacks, and species coexistence across timescales
# zachmiller@uchicago.edu

library("tidyverse")
library("reshape2")
library("ggpubr")
source("functions.R")


##### FIG 2 #####

set.seed(34)

### Set / draw parameters

n <- 3 # number of species (must be 3 for visualization with ternary plot)

# Sample P matrix with y > 0 (required for feasibility) and stable at small m and/or d values (for endo dynamics)
stab <- FALSE # flag for stability
while(!stab){
  
  #Sample P with y > 0 and find feasible m
  P <- sample_P(n, check.feas = TRUE, check.stab = FALSE, # NOTE: we will check stability "by hand" after adding a constant to P (below)
                m = "proportion", eps = 0.1)$P # eps is an arbitrary small number here
  
  rsPinv <- solve(P, rep(1, n))
  m <- 0.1 * (1 / sum(rsPinv)) # set m to be a small fraction (e.g. 10%) of max value
  
  P <- P + m # ensure that every species has net positive growth in every patch type
  
  # Check stability condition (from SI)
  M <- diag(rsPinv) %*% t(P)
  lambdas <- eigen(M)$val
  stab <- all((Re(lambdas) / Mod(lambdas)^2)[-1] <= 1)
}

# Pre-compute matrices for evaluating x* (eq. 3)
A <- solve(M)
b <- -m * rowSums(solve(t(P)))

### Compute data

# Loop over all w values 
n_points <- 1001
out_full <- matrix(nrow = n_points^2, ncol = 5) # initialize matrix to hold w1, w2, landscape entropy, species diversity, and feasibility
counter <- 1
for(w1 in seq(0.001, 0.999, length.out = n_points)){
  for(w2 in seq(0.001, 1 - w1, length.out = n_points)){
    
    w3 <- 1 - w1 - w2
    w <- c(w1, w2, w3)
    
    x <- A  %*% w + b # compute x* given w
    
    out_full[counter, ] <- c(w1, w2, 
                             -sum(w * log(w)), # landscape entropy
                             if(all(x > 0)) -sum(x * log(x)) else NA, # species diversity (entropy) if all coexisting 
                             all(x > 0)) # feasiblity
    counter <- counter + 1
  }
}
out_full <- out_full %>% as_tibble() %>% setNames(c("w1", "w2", "entropy", "diversity", "feasible")) # transform to tibble and name cols

# Compute w* (eq. S21)
perron <- Re(eigen(M)$vectors[, 1])
perron <- perron / sum(perron)
w_star <- m * rsPinv + (1 - m * sum(rsPinv)) * perron

# Compute x and diversity along a transect from low to high landscape entropy (chosen to include w*)
n_points <- 1000
out_transect <- matrix(nrow = n_points, ncol = 9) # initialize matrix to hold p, w's, entropy, diversity, feasiblity, and x's
counter <- 1
for(p in seq(0, 5, length.out = n_points)){ # range of exponents
    
    w <- w_star^p / sum(w_star^p) # generate w's of varying entropy
    
    x <- A %*% w + b # compute x* given w
    
    out_transect[counter, ] <- c(p, # exponent used
                                 w[1], w[2], 
                                 -sum(w * log(w)), # landscape entropy
                                 if(all(x > 0)) -sum(x * log(x)) else NA, # species diversity (entropy) if all coexisting 
                                 all(x > 0), # feasibility
                                 x[1], x[2], x[3]) # abundance of each species
    counter <- counter + 1
}
out_transect <- out_transect %>% as_tibble() %>% 
  setNames(c("p", "w1", "w2", "entropy", "diversity", "feasible", "x1", "x2", "x3")) # transform to tibble and name cols

# Compute range of p over which x is feasible
range_p <- out_transect %>% filter(feasible == 1) %>% pull(p) %>% range()
extend_range <- c(0.8 * range_p[1], 1.3 * range_p[2]) # extend range slightly (for plotting)

# Compute trajectories for endogenous feedback model
max_time <- 1000
times <- seq(0, max_time, by = 1) # vector of times for numerical integration

dynamics_out <- tibble(d = numeric(), i = numeric(), t = numeric(), w1 = numeric(), w2 = numeric()) # initialize tibble to hold time series
for(i in 1:n){ # loop over different initial conditions (here, one IC corresponding to each species)
  
  # Set up initial conditions
  # Chosen so that w_i = 0.9, w_j = w_k = 0.05 for each i
  X_init <- matrix(0.01, n, n)
  X_init[, i] <- 0.18
  
  x_init <- rowSums(X_init)
  y_init <- 2 * X_init[1, ]
  w_init <- y_init + colSums(X_init)
  
  # Integrate dynamics
  time_series <- ode(y = c(x_init, y_init, w_init), times = times, 
                     parms = list(n = n, P = P, m = m, d = d), 
                     func = general_feedback_model, method = "ode45")
  
  dynamics_out <- dynamics_out %>% add_row(d = d,
                                           i = i,
                                           t = time_series[, 1],
                                           w1 = time_series[, 2 * n + 2],
                                           w2 = time_series[, 2 * n + 3])
}


### Make plots

# Define coordinate transforms for ternary plotting
x_p <- function(x, y) 0.5 * (1 - x + y)
y_p <- function(x, y) sqrt(3)/2 * (1 - x - y)

# Create common plot elements
p <- out_full %>% filter(feasible == 1) %>% # only plot feasible points
  mutate(x_pos = x_p(w1, w2), # transform to ternary coordinates
         y_pos = y_p(w1, w2)) %>%
  ggplot() + 
  aes(x = x_pos, y = y_pos) + 
  annotate("polygon", x = c(0, 0.5, 1, 0), y = c(0, sqrt(3)/2, 0, 0), fill = "gray90", color = "black") + # create ternary plot background
  annotate("text", label = "w[1]", 
           x = x_p(1, 0) - 0.02, y = y_p(1, 0) - 0.04, # adjust label coordinates slightly for nice looking figure
           size = 5, parse = TRUE) + 
  annotate("text", label = "w[2]", 
           x = x_p(0, 1) + 0.03, y = y_p(0, 1) - 0.04, # adjust label coordinates slightly for nice looking figure
           size = 5 ,parse = TRUE) + 
  annotate("text", label = "w[3]", 
           x = x_p(0, 0), y = y_p(0, 0) + 0.03, # adjust label coordinates slightly for nice looking figure
           size = 5, parse = TRUE) + 
  geom_point(aes(color = diversity), size = 0.5) + # color points by species diversity
  scale_color_viridis_c(name = "Species diversity") + 
  theme_void()

# Plot feasiblity region and species diversity as a function of w's (for exogeneous dynamics)
p_1 <- p + 
  geom_path(data = out_transect %>% filter(p > extend_range[1], p < extend_range[2]) %>% # add transect line
              mutate(x_pos = x_p(w1, w2),
                     y_pos = y_p(w1, w2)), linetype = "dashed") + 
  theme(legend.position = "none")

# Pre-compute ternary coordinates for dynamical trajectories (with endogenous heterogeneity)
dynamics_out <- dynamics_out %>% mutate(x_pos = x_p(w1, w2),
                                        y_pos = y_p(w1, w2))

# Plot feasiblity region, species diversity, and w trajectories (for endogenous dynamics)
p_2 <- p + # plot trajectories 3 times each (up to t = 3, 25 or full) to get intermediate arrows
  geom_path(data = dynamics_out %>% filter(i == 1, t < 3), color = "red", 
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
  geom_path(data = dynamics_out %>% filter(i == 1, t < 25), color = "red",
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
  geom_path(data = dynamics_out %>% filter(i == 1), color = "red") + 
  geom_path(data = dynamics_out %>% filter(i == 2, t < 3), color = "purple",
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
  geom_path(data = dynamics_out %>% filter(i == 2, t < 25), color = "purple",
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
  geom_path(data = dynamics_out %>% filter(i == 2), color = "purple") + 
  geom_path(data = dynamics_out %>% filter(i == 3, t < 3), color = "darkorange",
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
  geom_path(data = dynamics_out %>% filter(i == 3, t < 25), color = "darkorange",
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
  geom_path(data = dynamics_out %>% filter(i == 3), color = "darkorange") + 
  annotate("point", x = x_p(w_star[1], w_star[2]), y = y_p(w_star[1], w_star[2]), size = 2) + # add w*
  theme(legend.position = c(0.85, 0.78),
        legend.key.size = unit(0.5, 'cm'))

# Inset plot showing landscape vs. species diversity
p_inset <- out_transect %>% filter(feasible == 1) %>%
  ggplot() + 
  aes(x = entropy, y = diversity) + 
  geom_vline(xintercept = out_transect %>% filter(p %in% range_p) %>% pull(entropy), # add upper and lower limits
             linetype = "dashed", color = "gray") +
  geom_line(size = 1) + 
  xlab("Landscape entropy") + 
  ylab("Species diversity") + 
  theme_classic() + 
  theme(plot.background = element_rect(fill = "transparent", # make plot background transparent to avoid plot overlap
                            color = NA_character_),
        panel.background = element_rect(color = "black", size = 1)) # add black blox around plot

# combine plots
ggarrange(plotlist = list(p_1, NULL, p_2), nrow = 1, widths = c(1, 0.05, 1)) + # add small space in between ternary plots
  annotation_custom(ggplotGrob(p_inset), # add inset plot
  xmin = 0.33, xmax = 0.63, ymin = 0.48, ymax =0.95) + 
  annotate("segment", x = 0.25, xend = 0.39, y = 0.23, yend = 0.52, # add an arrow pointing from transect line to inset plot
           arrow = arrow(length = unit(0.25, "cm"))) + 
  annotate("text", label = "A", x = 0.15, y = 0.95, size = 6, parse = TRUE) + # add panel labels
  annotate("text", label = "B", x = 0.31, y = 0.95, size = 6, parse = TRUE) + 
  annotate("text", label = "C", x = 0.68, y = 0.95, size = 6, parse = TRUE)


##### FIG 3 #####

set.seed(140)

n <- 5 # number of species

# Define a function to sample diagonal elements from a distribution with larger values (to ensure species are habitat specialists for example figure)
diag.dist <- function(n) {
  return(runif(n, 1, 3))
}

P <- sample_P(n, check.feas = TRUE, check.stab = FALSE, m = "proportion", eps = 0.5,
              distribution = runif, diag.distribution = diag.dist)$P
m <- 0.5 * (1 / sum(solve(P))) # set m to half of max value

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
w <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose w according to Eq. S28 to ensure feasibility

# Image of P matrix
tidy_P <- melt(P)
colnames(tidy_P) <- c("row", "col", "val")
p1 <- tidy_P %>% ggplot() + 
  aes(x = row, y = col, fill = val) + 
  geom_tile(color = "white", size = 1) + 
  scale_y_reverse() + 
  scale_fill_viridis_c() + 
  guides(fill = "none") +
  theme_void()

show(p1)

# Image of x %o% y matrix
x <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% w - m * rep(1, n))
y <- m * rowSums(solve(P))
xy <- as.vector(x) %*% t(as.vector(y))
tidy_xy <- melt(xy)
colnames(tidy_xy) <- c("row", "col", "val")
p2 <- tidy_xy %>% ggplot() + 
  aes(x = row, y = col, fill = val) + 
  geom_tile(color = "white", size = 1) +  
  scale_y_reverse() + 
  scale_fill_viridis_c() + 
  guides(fill = "none") +
  theme_void()

show(p2)

# Image of X matrix (element-wise product of P and x %o% y)
X <- P * xy
tidy_X <- melt(X)
colnames(tidy_X) <- c("row", "col", "val")
p3 <- tidy_X %>% ggplot() + 
  aes(x = row, y = col, fill = val) + 
  geom_tile(color = "white", size = 1) + 
  scale_y_reverse() + 
  scale_fill_viridis_c() + 
  guides(fill = "none") +
  theme_void()

show(p3)


##### FIG 4 #####

set.seed(10)

n <- 3 # number of species
m <- 0.1

# Sample P to be feasible and stable (so species coexist in example dynamics)
P <- sample_P(n, check.feas = TRUE, check.stab = TRUE,
              symmetric = TRUE, m = 0.5)$P

d <- 0.00001 # choose d small

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
v <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose v according to Eq. S28 to ensure feasibility

# Sample initial w near v to ensure feasiblity of the x(w(0))
x_of_w <- -1 
while(any(x_of_w < 0)) { # resample until feasible parameters are found
  w_init <- v * runif(n, 0.5, 1.5) # perturb equilibrium setting
  w_init <- (w_init / sum(w_init)) # normalize
  x_of_w <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% w_init - m * rep(1, n)) # compute x(w(0))
}

# Choose random initial x and y
y_init <- 0.5 * w_init
x_init <- runif(n)
x_init <- (x_init / sum(x_init)) * (1 - sum(y_init)) # normalize

# Set up numerical integration
times <- seq(0, 200000, 5) # vector of times
pars <- list(n = n, P = P, m = m, d = d)
init <- c(x_init, y_init, w_init)

# Integrate the general model (Eq. 5)
out <- ode(y = init, times = times, parms = pars, 
           func = general_feedback_model, method = "ode45")

# Define slow system for w
A <- solve(diag(rowSums(solve(P))) %*% t(P)) - diag(n)
b <- -m * rowSums(solve(t(P)) - solve(P))

# Integrate slow system (Eq. 7)
out_slow <- ode(y = w_init, times = times, 
                parms = list(n = n, P = P, m = m, d = d, A = A, b = b),
                func = slow_system, method = "ode45")

# Collect w values from full model and slow system
w_predicted <- cbind(1:nrow(out), out_slow[, 2:(n + 1)])
colnames(w_predicted) <- c("time", 1:n)
tidy_w_predicted <- pivot_longer(as_tibble(w_predicted), -time)
colnames(tidy_w_predicted) <- c("time", "species", "predicted")

w_actual <- cbind(1:nrow(out), out[, (2 * n + 2):(3 * n + 1)])
colnames(w_actual) <- c("time", 1:n)
tidy_w_actual <- pivot_longer(as_tibble(w_actual), -time)
colnames(tidy_w_actual) <- c("time", "species", "actual")

dat_w <- inner_join(tidy_w_actual, tidy_w_predicted) %>% mutate(type = "Patch states (w)") # combine for plotting

# Collect x values from full model and predicted from slow system
w <- w_predicted[, -1]
x_predicted <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% t(w) - m)

x_predicted <- cbind(1:nrow(out), t(x_predicted))
colnames(x_predicted) <- c("time", 1:n)
tidy_x_predicted <- pivot_longer(as_tibble(x_predicted), -time)
colnames(tidy_x_predicted) <- c("time", "species", "predicted")

x_actual <- cbind(1:nrow(out), out[, 2:(n + 1)])
colnames(x_actual) <- c("time", 1:n)
tidy_x_actual <- pivot_longer(as_tibble(x_actual), -time)
colnames(tidy_x_actual) <- c("time", "species", "actual")

dat_x <- inner_join(tidy_x_actual, tidy_x_predicted) %>% mutate(type = "Species (x)") # combine for plotting

dat <- dat_w %>% add_row(dat_x) # combine w and x for plotting
dat$type = factor(dat$type, levels=c('Species (x)', 'Patch states (w)'))

# Plot all time series for comparison with limiting case model
scale_round <- function(x) sprintf("%.2f", x)
dat %>%
  ggplot() +
  aes(x = time) + 
  geom_line(aes(y = actual, group = species, color = species), size = 1, alpha = 0.9) +
  geom_line(aes(y = predicted, group = species), lty = "dashed") +
  xlab("Time") + 
  ylab("Frequency") + 
  scale_x_log10() + 
  scale_y_continuous(labels = scale_round) + 
  guides(color = "none") + 
  facet_grid(type~., scales = "free") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white"))


##### FIG 5 #####

### Panel A

set.seed(18)

n <- 3 # number of species
m <- 0.5

# Sample P to be feasible and stable (so species coexist in example dynamics)
P <- sample_P(n, check.feas = TRUE, check.stab = TRUE, symmetric = TRUE, m = 0.5)$P

d <- 0.00001 # choose d small

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
v <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose v according to Eq. S28 to ensure feasibility

# Choose w(0) as a random perturbation of v (equilibrium setting) without checking feasibility
w_init <- v * runif(n, 0.5, 1.5)
w_init <- (w_init / sum(w_init))

# Choose random initial x and y
y_init <- 0.8 * w_init
x_init <- runif(n)
x_init <- (x_init / sum(x_init)) * (1 - sum(y_init))

# Set up numerical integration
times <- seq(0, 1000000, 20)
pars <- list(n = n, P = P, m = m, d = d)
init <- c(x_init, y_init, w_init)

# Integrate the general model (Eq. 5)
out <- ode(y = init, times = times, parms = pars, 
           func = general_feedback_model, method = "ode45")

dat <- cbind(1:nrow(out), out[, 2:(n + 1)])
colnames(dat) <- c("time", 1:n)
dat <- pivot_longer(as_tibble(dat), -time)
colnames(dat) <- c("time", "species", "value")

# Plot x time series
scale_round <- function(x) sprintf("%.2f", x)
p_1 <- dat %>% filter(time > 30) %>%
  ggplot() +
  aes(x = time) + 
  geom_line(aes(y = value, group = species, color = species), size = 1, alpha = 0.9) +
  xlab("Time") + 
  ylab("Frequency (x)") + 
  geom_hline(yintercept = 0, color = "gray", size = 1.2, alpha = 0.7) + 
  scale_y_continuous(labels = scale_round) + 
  guides(color = "none") + 
  theme_bw() + 
  theme(panel.background = element_rect(color = "black", size = 1)) # add black blox around plot

### Panel B

set.seed(30)

n <- 3 # number of species
m <- 0.5

# Sample P without regard for stability
P <- sample_P(n, check.feas = TRUE, check.stab = FALSE, symmetric = TRUE, m = 0.5)$P
P <- P + diag(n) # increase diagonal to make coexistence unstable

d <- 0.00001 # choose d small

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
v <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose v according to Eq. S28 to ensure feasibility

# Choose w near v (equilibrium setting) to avoid switching dynamics
w_init <- v * runif(n, 0.98, 1.02)
w_init <- (w_init / sum(w_init))

# Start dynamics on slow manifold to avoid distracting transients
x_init <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% w_init - m * rep(1, n))
y_init <- m %*% rowSums(solve(P))
y_init <- (y_init / sum(y_init)) * (1 - sum(x_init)) # normalize

# Integrate for some time before perturbation
times <- seq(0, 20000, 5)
pars <- list(n = n, P = P, m = m, d = d)
init <- c(x_init, y_init, w_init)

# Integrate the general model (Eq. 5)
out1 <- ode(y = init, times = times, parms = pars, 
            func = general_feedback_model, method = "ode45")

# Perturb current species abundances
perturb_idx <- 2 # perturb the species that will go extinct (NOTE: chosen by hand for this seed)
current <- out1[nrow(out1), -1] # take the current state of the system
y <- current[(n+1):(2*n)]
x <- current[perturb_idx] 
perturb_dist <- y * P[perturb_idx, ] # distribution of patches occupied by the perturbed species
current[perturb_idx] <- x * 0.05 # reduce frequency by 95%
current[(n+1):(2*n)] <- y + x * (1 - 0.05) * (perturb_dist / sum(perturb_dist)) # add newly emptied patches to y 

# Integrate for longer following perturbation
out2 <- ode(y = current, times = seq(0, 50000, 5), parms = pars, 
            func = general_feedback_model, method = "ode45")

out <- rbind(out1[, -1],
             matrix(current, nrow = 50, ncol = 3 * n, byrow = TRUE),
             out2[, -1]) # combine time series

# Reformat for plotting
dat <- cbind(1:nrow(out), out[, 1:n])
colnames(dat) <- c("time", 1:n)
dat <- pivot_longer(as_tibble(dat), -time)
colnames(dat) <- c("time", "species", "value")

# Plot combined time series to show pre- and post-perturbation
scale_round <- function(x) sprintf("%.2f", x)
p_2 <- dat %>% filter(time > 1) %>%
  ggplot() +
  aes(x = time) + 
  geom_line(aes(y = value, group = species, color = species), size = 1, alpha = 0.9) +
  geom_vline(xintercept = length(times), lty = "dashed") + 
  xlab("Time") + 
  ylab("Frequency (x)") + 
  scale_y_continuous(labels = scale_round) + 
  geom_hline(yintercept = 0, color = "gray", size = 1.2, alpha = 0.7) +
  guides(color = "none") + 
  theme_bw() + 
  theme(panel.background = element_rect(color = "black", size = 1)) # add black blox around plot

# Combine panels
ggarrange(plotlist = list(p_1, p_2), labels = "AUTO", font.label = list(face = "plain"))

  
##### FIG S1 #####

plot_circle <- function(center, radius){
  
  angles <- seq(0, 2 * pi, 0.0001)
  
  out <- tibble(x = radius * cos(angles) + center[1],
                y = radius * sin(angles) + center[2])
  return(out)
}

outer_circ <- plot_circle(c(0,0), 1)
outer_circ %>% ggplot() + aes(x = x, y = y) + 
  geom_vline(xintercept = 0, alpha = 0.5) + 
  geom_hline(yintercept = 0, alpha = 0.5) + 
  geom_polygon(alpha = 0.2, size = 1, fill = "blue", color = "black") + 
  geom_polygon(data = plot_circle(c(1/2, 0), 1/2), 
               alpha = 0.4, size = 1, fill = "red", color = "black") + 
  geom_point(x = 1, y = 0, pch = 4, size = 3) + 
  xlab("Re") + 
  ylab("Im") +
  coord_equal() + 
  theme_bw()
