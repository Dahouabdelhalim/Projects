################################################################
### 01_fit_model_eDNA.R
###     Fit the tracer model to the eDNA survey dataset.
################################################################

### Set working directory
setwd("/specify/your/directory")


### Load data
load("tracer_matrix.Rdata")                        # Tracer model matrix
cell_attributes <- read.csv("cell_attributes.csv") # Cell attributes
eDNA <- read.csv("eDNA_measurements.csv")          # eDNA data


### Process data, and define some constants and variables
eDNA  <- eDNA[eDNA$conc > 0, ]          # Omit negative samples
N     <- nrow(eDNA)                     # Number of eDNA samples
I     <- nrow(A)                        # Number of grid cells
log_y <- log(eDNA$conc)                 # Logarithms of eDNA concentration
B     <- A[eDNA$cell, ]                 # The design matrix B -- see the article
v     <- cell_attributes$vol_m.3 * 1000 # Water volume of each grid cell (unit: L)


### Fit the model using RStan
library(rstan)

## Data 
data <- list(N     = N,
             I     = I,
             log_y = log_y,
             B     = B,
             v     = v)

## Initial values
init <- function() {
    list(mu    = 2E-9  + rnorm(1, sd = 1E-10),
         log_x = rnorm(I, -20.7, 1),
         tau   = 7.4 + rnorm(1, sd = 1E-1),
         sigma = exp(0.2) + rnorm(1, sd = 1E-2))
}

## Parameters to save
pars <- c("mu", "tau", "sigma", "log_x", "abundance")

## MCMC Settings
n.chains <- 3
n.thin   <- 10
n.adapt  <- 1000
n.iter   <- n.adapt + 1000 * n.thin
max_treedepth <- 30

## Model fitting
rstan_options(auto_write = TRUE)
options(mc.cores = n.chains)
start   <- Sys.time()
fit     <- stan(file = "model_eDNA.stan", data = data, init = init, pars = pars,
                chains = n.chains, iter = n.iter, warmup = n.adapt, thin = n.thin,
                control = list(max_treedepth = max_treedepth))
end     <- Sys.time()
runtime <- end - start
errmsg  <- warnings()

## Save results
save(fit, runtime, errmsg, file = "res_model_eDNA.Rdata")


### Posterior convergence
## R-hat statistics
sort(summary(fit)$summary[, "Rhat"])

## Number of divergent transitions
library(bayesplot)
np_fit <- nuts_params(fit)
table(subset(np_fit, np_fit$Parameter == "divergent__")$Value)


### Goodness-of-fit assessment
post <- apply(extract(fit, c("sigma", "log_x"), permuted = FALSE), 3, c)
chi2_rep <- chi2_obs <- vector(length = nrow(post))
chi2 <- function(y, mu) sum((y - mu)^2 / mu)

for (i in seq_len(nrow(post))) {
    # Simulate replicated data
    mu_y      <- log(c(B %*% exp(post[i, -1])))
    log_y_rep <- rnorm(length(log_y), mean = mu_y, sd = post[i, "sigma"])

    # Calculate chi-squared discrepancy measure
    chi2_rep[i] <- chi2(log_y_rep, mu_y) # for simulated data set
    chi2_obs[i] <- chi2(log_y, mu_y)     # for observed data set
}

## Bayesian p-value
sum(chi2_obs < chi2_rep) / nrow(post)

