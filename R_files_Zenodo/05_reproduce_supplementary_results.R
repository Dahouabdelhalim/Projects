################################################################
### 05_reproduce_supplementary_results.R
###     Post-process the model fit objects to obtain the results
###     shown in the supplementary file.
################################################################

### Set working directory
setwd("/specify/your/directory")


### Load data
cell_attributes <- read.csv("cell_attributes.csv")       # Cell attributes
coarser_cells   <- read.csv("coarser_cells.csv")         # Coarse cell ID
d               <- read.csv("acoustic_measurements.csv") # Acoustic data


### Process data, and define some constants and variables
I   <- nrow(cell_attributes)           # Number of grid cells
v   <- cell_attributes$vol_m.3 * 1000  # Water volume of each grid cell (unit: L)
vm3 <- cell_attributes$vol_m.3         # Water volume of each grid cell (m^3)
r   <- unique(d$cell[is.finite(d$Sv)]) # Unique set of cells with acoustic measurements


### Load packages and functions
library(rstan)
library(coda)
library(RColorBrewer)
library(tagcloud)
library(celestial)
source("functions.R")


### Load model fit objects
load("res_model_eDNA.Rdata")
fit_100 <- fit

load("res_model_eDNA_200.Rdata")
fit_200 <- fit

load("res_model_eDNA_300.Rdata")
fit_300 <- fit

load("res_model_eDNA_400.Rdata")
fit_400 <- fit

load("res_model_echosounder.Rdata")
fit_echo <- fit


### Bay-scale abundance estimates and correlation (Table S1)

## Cell-level density
logx_100 <- apply(as.array(fit_100)[, , sprintf("log_x[%s]", 1:I)], 3, c)
logx_200 <- apply(as.array(fit_200)[, , sprintf("log_x[%s]", 1:I)], 3, c)
logx_300 <- apply(as.array(fit_300)[, , sprintf("log_x[%s]", 1:I)], 3, c)
logx_400 <- apply(as.array(fit_400)[, , sprintf("log_x[%s]", 1:I)], 3, c)

logx_echo <- apply(as.array(fit_echo)[, , sprintf("log_x_new[%s]", 1:I)], 3, c)
# Arrange cells in the original order
tmp       <- matrix(nrow = nrow(logx_echo), ncol = ncol(logx_echo))
tmp[, r]  <- logx_echo[, seq_along(r)]
tmp[, -r] <- logx_echo[, (length(r) + 1):I]
logx_echo <- tmp

## Cell-level abundance
abundance_100  <- t(exp(t(logx_100)) * v)
abundance_200  <- t(exp(t(logx_200)) * v)
abundance_300  <- t(exp(t(logx_300)) * v)
abundance_400  <- t(exp(t(logx_400)) * v)
abundance_echo <- t(exp(t(logx_echo)) * v)

## Bay-level abundance (not corrected)
est_100 <- cbind(median(apply(abundance_100, 1, sum)), HPDinterval(as.mcmc(apply(abundance_100, 1, sum))))
est_200 <- cbind(median(apply(abundance_200, 1, sum)), HPDinterval(as.mcmc(apply(abundance_200, 1, sum))))
est_300 <- cbind(median(apply(abundance_300, 1, sum)), HPDinterval(as.mcmc(apply(abundance_300, 1, sum))))
est_400 <- cbind(median(apply(abundance_400, 1, sum)), HPDinterval(as.mcmc(apply(abundance_400, 1, sum))))

rbind(est_100, est_200, est_300, est_400)   # Table S1 column 2

## Bay-level abundance (corrected)
# Vertically aggregate abundance estimates
aggregated_abundance_100 <- aggregate_abundance(abundance_100, cell_attributes)
aggregated_abundance_200 <- aggregate_abundance(abundance_200, cell_attributes)
aggregated_abundance_300 <- aggregate_abundance(abundance_300, cell_attributes)
aggregated_abundance_400 <- aggregate_abundance(abundance_400, cell_attributes)

# Identify fish market cells
market_100 <- 61
market_200 <- which(matrix(coarser_cells[[1]], ncol = 9) == coarser_cells[[1]][61])
market_300 <- which(matrix(coarser_cells[[2]], ncol = 9) == coarser_cells[[2]][61])
market_400 <- which(matrix(coarser_cells[[3]], ncol = 9) == coarser_cells[[3]][61])

est_100c <- cbind(median(apply(aggregated_abundance_100[, -market_100], 1, sum)),
                  HPDinterval(as.mcmc(apply(aggregated_abundance_100[, -market_100], 1, sum))))
est_200c <- cbind(median(apply(aggregated_abundance_200[, -market_200], 1, sum)),
                  HPDinterval(as.mcmc(apply(aggregated_abundance_200[, -market_200], 1, sum))))
est_300c <- cbind(median(apply(aggregated_abundance_300[, -market_300], 1, sum)),
                  HPDinterval(as.mcmc(apply(aggregated_abundance_300[, -market_300], 1, sum))))
est_400c <- cbind(median(apply(aggregated_abundance_400[, -market_400], 1, sum)),
                  HPDinterval(as.mcmc(apply(aggregated_abundance_400[, -market_400], 1, sum))))

rbind(est_100c, est_200c, est_300c, est_400c)   # Table S1 column 3

## Correlation between acoustic estimates
# Obtain vertically aggregated density
aggregated_density_100  <- aggregate_density(abundance_100, cell_attributes, vm3)
aggregated_density_200  <- aggregate_density(abundance_200, cell_attributes, vm3)
aggregated_density_300  <- aggregate_density(abundance_300, cell_attributes, vm3)
aggregated_density_400  <- aggregate_density(abundance_400, cell_attributes, vm3)
tmp                     <- aggregate_density(abundance_echo, cell_attributes, vm3, r, "acoustic")
aggregated_density_echo <- tmp[[1]]
avail_echo              <- tmp[[2]]

# Table S1 column 4
rbind(cor(x = aggregated_density_100[-market_100][avail_echo[-market_100]],
          y = aggregated_density_echo[-market_100][avail_echo[-market_100]]),
      cor(x = aggregated_density_200[-market_200][avail_echo[-market_200]],
          y = aggregated_density_echo[-market_200][avail_echo[-market_200]]),
      cor(x = aggregated_density_300[-market_300][avail_echo[-market_300]],
          y = aggregated_density_echo[-market_300][avail_echo[-market_300]]),
      cor(x = aggregated_density_400[-market_400][avail_echo[-market_400]],
          y = aggregated_density_echo[-market_400][avail_echo[-market_400]]))


### Estimated spatial distribution (Figure S1)
par(mfrow = c(2, 2))
draw_map(fit_200, cell_attributes, market_200, "coarser", 200)
draw_map(fit_300, cell_attributes, market_300, "coarser", 300)
draw_map(fit_400, cell_attributes, market_400, "coarser", 400)

