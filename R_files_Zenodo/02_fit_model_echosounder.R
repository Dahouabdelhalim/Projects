################################################################
### 02_fit_model_echosounder.R
###     Fit a linear mixed model to the acoustic survey dataset.
################################################################

### Set working directory
setwd("/specify/your/directory")


### Load data
d <- read.csv("acoustic_measurements.csv")         # Acoustic data
cell_attributes <- read.csv("cell_attributes.csv") # Cell attributes


### Process data, and define some constants and variables
## Logarithms of jack mackerel density (individuals / m^3)
logD <- (d$Sv / 10) * log(10) - (-59.6 / 10) * log(10)

## Adjust the scale of density to individuals / L
logD <- logD - log(1000)

## Constants
I   <- nrow(cell_attributes)          # Number of grid cells
N   <- length(logD[is.finite(logD)])  # Number of positive acoustic measurements 
n2i <- d$cell[is.finite(logD)]        # Cell ID for each acoustic measurements
R   <- length(unique(n2i))            # Number of grid cells with acoustic measurements
S   <- I - length(unique(n2i))        # Number of grid cells without acoustic measurements
v   <- cell_attributes$vol_m.3 * 1000 # Water volume of each grid cell (unit: L)

## Arrange cells with acoustic measurements in the first half and the rest in the second half
# Rearrange water volume vector
v        <- c(v[unique(n2i)], v[-unique(n2i)])
names(v) <- c(seq_len(I)[unique(n2i)], seq_len(I)[-unique(n2i)])

# Define a new cell indicator for each observation
n2j  <- vector(length = N)            # Cell ID for each acoustic measurements (starts with 1)
last <- 0
tab  <- table(n2i)
for (r in seq_len(R)) {
    n2j[(last + 1):(last + tab[r])] <- r
    last <- last + tab[r]
}


### Fit the model using RStan
library(rstan)

## Data 
data <- list(N     = N,
             I     = I,
             R     = R,
             n2j   = n2j,
             log_d = logD[is.finite(logD)],
             v     = v)

## Initial values
init <- function() {
    list(mu    = exp(mean(logD[is.finite(logD)])) + rnorm(1, sd = 1E-10),
         log_x = rnorm(R, mean(logD[is.finite(logD)]), 1),
         tau   = 7.4 + rnorm(1, sd = 1E-1),
         sigma = 2 + rnorm(1, sd = 1E-2))
}

## Parameters to save
pars <- c("mu", "tau", "sigma", "log_x_new", "abundance")

## MCMC settings
n.chains <- 3
n.thin   <- 1
n.adapt  <- 1000
n.iter   <- n.adapt + 1000 * n.thin

## Model fitting
rstan_options(auto_write = TRUE)
options(mc.cores = n.chains)
start   <- Sys.time()
fit     <- stan(file = "model_echosounder.stan", data = data, init = init, pars = pars,
                chains = n.chains, iter = n.iter, warmup = n.adapt, thin = n.thin)
end     <- Sys.time()
runtime <- end - start
errmsg  <- warnings()

save(fit, runtime, errmsg, file = "res_model_echosounder.Rdata")


### Posterior convergence
## R-hat statistics
sort(summary(fit)$summary[, "Rhat"])

## Number of divergent transitions
library(bayesplot)
np_fit <- nuts_params(fit)
table(subset(np_fit, np_fit$Parameter == "divergent__")$Value)


### Goodness-of-fit assessment
post <- apply(extract(fit, c("sigma", "log_x_new"), permuted = FALSE), 3, c)
chi2_rep <- chi2_obs <- vector(length = nrow(post))
chi2 <- function(y, mu) sum((y - mu)^2 / mu)

for (i in seq_len(nrow(post))) {
    # Simulate replicated data
    mu_logD  <- post[i, n2j + 1]
    logD_rep <- rnorm(length(logD[is.finite(logD)]), mean = mu_logD, sd = post[i, "sigma"])

    # Calculate chi-squared discrepancy measure
    chi2_rep[i] <- chi2(logD_rep, mu_logD)
    chi2_obs[i] <- chi2(logD[is.finite(logD)], mu_logD)
}

## Bayesian p-value
sum(chi2_obs < chi2_rep) / nrow(post)

