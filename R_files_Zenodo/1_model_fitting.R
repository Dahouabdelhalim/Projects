################################################################################
### 1_model_fitting.R
###     Fit the multispecies site occupancy model for the Kasumigaura fish eDNA
###     metabarcoding dataset.
###
###     The following libraries are used:
###         * JAGS: https://mcmc-jags.sourceforge.io/
###         * jagsUI: https://CRAN.R-project.org/package=jagsUI
###         * coda: https://CRAN.R-project.org/package=coda
###         * ggplot2: https://CRAN.R-project.org/package=ggplot2
################################################################################

setwd("<specify your directory>")


### Load data
load("data.Rdata")


### Set constants
I <- dim(y)[1]              # Number of species
J <- dim(y)[2]              # Number of sites
K <- dim(y)[3]              # Number of replicates
N <- apply(y, c(2, 3), sum) # Sequence depth
M <- 4                      # Order of species effects

# Set y to NA and N to 1 when sequence depth is zero (i.e., missing replicates)
for (j in 1:J) {
    for (k in 1:K) {
        if (N[j, k] == 0) {
            y[, j, k] <- NA
            N[j, k]   <- 1
        }
    }
}


### Set data list
data <- list(I        = I,
             J        = J,
             K        = K,
             y        = y,
             N        = N,
             M        = M,
             cov_site = riverbank,
             cov_sp   = mismatch,
             itau     = 1E-3,
             U        = 1E3)


### Set initial values
inits <- function() {
    list(z = matrix(1, I, J),
         u = array(1, dim = c(I, J, K)),
         x = array(rnorm(I * J * K, mean = 1, sd = 0.1), dim = c(I, J, K)),
         species_effects = matrix(rnorm(M * I, sd = 0.1), M, I),
         Mu     = rnorm(M, sd = 0.1),
         alpha1 = rnorm(1, sd = 0.1),
         sigma  = rep(1, M),
         rho    = matrix(0, M, M))
}


### Parameters monitored
params <- c("Mu", "sigma", "rho", "alpha1", "gamma", 
            "psi", "theta", "phi", "z", "pi")


### Run MCMC in JAGS (approx. 12 hours to complete)
library(jagsUI)
res <- jags(data, inits, params, "model.jags",
            n.chains = 6,
            n.adapt  = 1000,
            n.burnin = 30000,
            n.iter   = 280000,
            n.thin   = 500,
            parallel = TRUE)


### Save results
save(res, file = "result.Rdata")


### ------------------------------------------------------------------------ ###

#load("data.Rdata")
#I <- dim(y)[1]              # Number of species
#J <- dim(y)[2]              # Number of sites
#K <- dim(y)[3]              # Number of replicates
#N <- apply(y, c(2, 3), sum) # Sequence depth
#M <- 4                      # Order of species effects
#for (j in 1:J) {
#    for (k in 1:K) {
#        if (N[j, k] == 0) {
#            y[, j, k] <- NA
#            N[j, k]   <- 1
#        }
#    }
#}
#load("result.Rdata")
library(coda)
library(ggplot2)
source("functions.R")


### Posterior predictive check using the Freeman-Tukey statistics
gof(res, y)


### Posterior estimates of some key quantities
# mu_gamma1: effect of the absence of vegetation on site occupancy
median(res$sims.list$Mu[, 4])
HPDinterval(as.mcmc(res$sims.list$Mu[, 4]))

# alpha1: effect of primer-template mismatches on sequence dominance
median(res$sims.list$alpha1)
HPDinterval(as.mcmc(res$sims.list$alpha1))

# species ratio of sequence relative dominance phi, contrasting between the highest and lowest species
median(res$sims.list$phi[, dimnames(y)[[1]] == "Rhinogobius spp."] / res$sims.list$phi[, dimnames(y)[[1]] == "Salangichthys microdon"])
HPDinterval(as.mcmc(res$sims.list$phi[, dimnames(y)[[1]] == "Rhinogobius spp."] / res$sims.list$phi[, dimnames(y)[[1]] == "Salangichthys microdon"]))


### Correlation matrix of the four species random effects (Table 1)
cormat <- array(dim = c(4, 4, 3))
dimnames(cormat) <- list(c("alpha_0i", "logit theta_i", "gamma_0i", "gamma_1i"),
                          c("alpha_0i", "logit theta_i", "gamma_0i", "gamma_1i"),
                          c("post median", "95% HPDI lower limit", "95% HPDI upper limit"))
for (i in 2:4) {
    for (j in 1:(i - 1)) {
        cormat[i, j, 1] <- median(res$sims.list$rho[, j, i])
        cormat[i, j, 2:3] <- HPDinterval(as.mcmc(res$sims.list$rho[, j, i]))[1:2]
    }
}
cormat


## Species-level parameters (Figure 2)
plot_species_params(res, pal_rb_or)


## Posterior probabilities of species site occupancy (Figure 3)
plot_post_z(res, rev(pal_rb_or))


