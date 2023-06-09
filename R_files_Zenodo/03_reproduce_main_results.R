################################################################
### 03_reproduce_main_results.R
###     Post-process the model fit objects to obtain the results
###     shown in the main text.
################################################################

### Set working directory
setwd("/specify/your/directory")


### Load data
load("tracer_matrix.Rdata")                        # Tracer model matrix
cell_attributes <- read.csv("cell_attributes.csv") # Cell attributes
eDNA <- read.csv("eDNA_measurements.csv")          # eDNA data
d    <- read.csv("acoustic_measurements.csv")      # Acoustic data


### Process data, and define some constants and variables
eDNA  <- eDNA[eDNA$conc > 0, ]          # Omit negative samples
N     <- nrow(eDNA)                     # Number of eDNA samples
I     <- nrow(A)                        # Number of grid cells
log_y <- log(eDNA$conc)                 # Logarithms of eDNA concentration
B     <- A[eDNA$cell, ]                 # The design matrix B -- see the article
v     <- cell_attributes$vol_m.3 * 1000 # Water volume of each grid cell (unit: L)
n2i   <- d$cell[is.finite(d$Sv)]        # Cell ID for each acoustic measurements
r     <- unique(n2i)                    # Unique set of cells with acoustic measurements


### Load packages and functions
library(rstan)
library(coda)
library(RColorBrewer)
library(tagcloud)
library(celestial)
source("functions.R")


### Load model fit objects
load("res_model_eDNA.Rdata")
fit_edna <- fit

load("res_model_echosounder.Rdata")
fit_echo <- fit


### Vertically aggregate abundance estimates
logx <- apply(as.array(fit_edna)[, , sprintf("log_x[%s]", 1:I)], 3, c)
abundance <- t(exp(t(logx)) * v)
aggregated_abundance <- aggregate_abundance(abundance, cell_attributes)


### Posterior estimates of bay-level abundance (Table 1)
## eDNA estimate (not corrected)
est_eDNA1 <- cbind(summary(fit_edna)$summary["abundance", "50%"],
                   HPDinterval(as.mcmc(c(as.array(fit_edna)[, , "abundance"]))))

## eDNA estimate, the fish market cells removed
est_eDNA2 <- cbind(median(apply(aggregated_abundance[, -61], 1, sum)),
                   HPDinterval(as.mcmc(apply(aggregated_abundance[, -61], 1, sum))))

## acoustic estimate
est_acoustic <- cbind(summary(fit_echo)$summary["abundance", "50%"],
                      HPDinterval(as.mcmc(c(as.array(fit_echo)[, , "abundance"]))))

## Table 1
rbind(est_eDNA1, est_eDNA2, est_acoustic)


### Posterior estimates of aggregated abundance (Figure 2b)
med  <- apply(log10(aggregated_abundance), 2, median)
lci  <- apply(log10(aggregated_abundance), 2, function(x) HPDinterval(as.mcmc(x))[1, 1])
uci  <- apply(log10(aggregated_abundance), 2, function(x) HPDinterval(as.mcmc(x))[1, 2])

plot(sort(med, decreasing = TRUE), pch = 16, type = "l",
     xlab = "Abundance rank", ylab = "log10 fish abundance", ylim = range(c(lci, uci)))
arrows(seq_along(med), y0 = lci[order(med, decreasing = TRUE)],
       y1 = uci[order(med, decreasing = TRUE)], length = 0, col = "gray")
points(seq_along(med), sort(med, decreasing = TRUE), pch = 21, col = "gray40",
       bg = c("red", rep("gray80", length(med))))


### Predicted vs observed eDNA concentration (Figure 4)
## Calculate predicted eDNA concentration
logy_pred <- matrix(nrow = dim(fit_edna)[1] * dim(fit_edna)[2], ncol = nrow(B))
for (i in seq_len(nrow(logx))) {
    logy_pred[i, ] <- c(log(B %*% exp(logx[i, ])))
}

## Draw figure
col <- brewer.pal(6, "Pastel1")[2]
xlim <- range(apply(log10(exp(logy_pred)), 2, quantile, 0.025),
              apply(log10(exp(logy_pred)), 2, quantile, 0.975))
plot(log10(exp(log_y)) ~ log10(exp(apply(logy_pred, 2, median))),
     xlim = xlim, type = "n", 
     xlab = "Predicted eDNA concentration [copies/L]",
     ylab = "Observed eDNA concentration [copies/L]")

il <- apply(log10(exp(logy_pred)), 2, function(x) HPDinterval(as.mcmc(x))[1])
iu <- apply(log10(exp(logy_pred)), 2, function(x) HPDinterval(as.mcmc(x))[2])
arrows(x0 = il, x1 = iu, y0 = log10(exp(log_y)), length = 0, col = "gray")

par(new = TRUE)
plot(log10(exp(log_y)) ~ log10(exp(apply(logy_pred, 2, median))),
     xlim = xlim, bg = col, pch = 21, ann = FALSE, axes = FALSE)

abline(a = 0, b = 1)

## Posterior estimate of the correlation
post_rho <- vector(length = dim(fit_edna)[1] * dim(fit_edna)[2])
for (i in seq_along(post_rho)) {
    post_rho[i] <- cor(x = log10(exp(log_y)), y = log10(exp(logy_pred[i, ])))
}

# Posterior median and 95% credible interval
c(median(post_rho), HPDinterval(as.mcmc(post_rho)))


## Estimated spatial distribution (Figure 5)
par(mfrow = c(1, 2))
draw_map(fit_edna, cell_attributes, 61, "eDNA")
draw_map(fit_echo, cell_attributes, NA, "acoustic")

