################################################################
### 04_fit_coarser_models.R
###     Fit the tracer model to the eDNA survey dataset in which
###     fish density is estimated in coarser grid blocks.
################################################################

### Set working directory
setwd("/specify/your/directory")


### Set the spatial resolution of the inference
resolution <- 200  # 200 m grid cells
#resolution <- 300 # 300 m grid cells
#resolution <- 400 # 400 m grid cells


### Load data
load("tracer_matrix.Rdata")                        # Tracer model matrix
cell_attributes <- read.csv("cell_attributes.csv") # Cell attributes
eDNA <- read.csv("eDNA_measurements.csv")          # eDNA data
coarser_cells <- read.csv("coarser_cells.csv")     # Coarse cell ID


### Process data, and define some constants and variables
eDNA  <- eDNA[eDNA$conc > 0, ]          # Omit negative samples
N     <- nrow(eDNA)                     # Number of eDNA samples
I     <- nrow(A)                        # Number of grid cells
log_y <- log(eDNA$conc)                 # Logarithms of eDNA concentration
B     <- A[eDNA$cell, ]                 # The design matrix B -- see the article
v     <- cell_attributes$vol_m.3 * 1000 # Water volume of each grid cell (unit: L)

# Extract an integer vector mapping each original cell to coarser cells
if (resolution == 200) {
    map <- coarser_cells[[1]]
} else if (resolution == 300) {
    map <- coarser_cells[[2]]
} else if (resolution == 400) {
    map <- coarser_cells[[3]]
}


### Fit the model using RStan
library(rstan)

## Data 
data <- list(N     = N,
             I     = I,
             J     = max(map),
             j_i   = map,
             log_y = log_y,
             B     = B,
             v     = v)

## Initial values
init <- function() {
    list(mu      = 2E-9  + rnorm(1, sd = 1E-10),
         epsilon = rnorm(max(map)),
         tau     = 7.4 + rnorm(1, sd = 1E-1),
         sigma   = exp(0.2) + rnorm(1, sd = 1E-2))
}

## Parameters to save
pars <- c("mu", "tau", "sigma", "log_x", "abundance")

## MCMC settings
n.chains <- 3
n.thin   <- 5
n.adapt  <- 1000
n.iter   <- n.adapt + 1000 * n.thin
max_treedepth <- 10

## Model fitting
rstan_options(auto_write = TRUE)
options(mc.cores = n.chains)
start   <- Sys.time()
fit     <- stan(file = "model_eDNA_coarser.stan", data = data, init = init, pars = pars,
                chains = n.chains, iter = n.iter, warmup = n.adapt, thin = n.thin,
                control = list(max_treedepth = max_treedepth))
end     <- Sys.time()
runtime <- end - start
errmsg  <- warnings()

## Save results
save(fit, runtime, errmsg, file = sprintf("res_model_eDNA_%s.Rdata", resolution))


