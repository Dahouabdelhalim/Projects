library(ubms)

# Load data
grouse <- read.csv("grouse_data.csv")

# Stack site covariates
site_covs <- rbind(grouse[,1:2], grouse[,1:2], grouse[,1:2])
site_covs$site <- factor(site_covs$site)
site_covs$year <- factor(rep(1:3, each=nrow(grouse)))

# Stack count data
y1 <- grouse[,3:7]
y2 <- grouse[,8:12]
y3 <- grouse[,13:17]
colnames(y1) <- colnames(y2) <- colnames(y3) <- paste0("count",1:5)
y <- rbind(y1, y2, y3)

# Stack date data
d1 <- grouse[,18:22]
d2 <- grouse[,23:27]
d3 <- grouse[,28:32]
colnames(d1) <- colnames(d2) <- colnames(d3) <- paste0("date",1:5)
date <- rbind(d1, d2, d3)

# Build unmarkedFrame
umf <- unmarkedFramePCount(y=y, siteCovs=site_covs, obsCovs=list(date=date))

# Keep only sites without missing covariate values
has_NA <- is.na(site_covs$pctaspen)
umf <- umf[which(!has_NA),]

# Fit models
# All models include date as a covariate on detection
# All covariates standardized with scale()

# Model with no covariates or random effect on abundance
mod_null <- stan_pcount(~scale(date) ~1, umf, K=20,
                        chains=3, iter=2000, seed=123)

# Model with random intercepts by site
mod_rand <- stan_pcount(~scale(date) ~(1|site), umf, K=20,
                        chains=3, iter=2000, seed=123)

# Model with random intercepts and percent aspen as a covariate on abundance
mod_aspen <- stan_pcount(~scale(date) ~scale(pctaspen) + (1|site), umf, K=20,
                         chains=3, iter=2000, seed=123)

# Save model objects
save(mod_null, mod_rand, mod_aspen, file="grouse_models.Rdata")
load("grouse_models.Rdata")

# Diagnostics and model fit

# Look at effective sample size and Rhat for mod_aspen
mod_aspen

# Look at traceplots for mod_aspen
traceplot(mod_aspen)

# Residual plots
plot(mod_aspen)

# Posteriors
plot_posteriors(mod_aspen)

# Goodness-of-fit
set.seed(123)
(fit <- gof(mod_aspen, nsim=1000))
plot(fit)

# Inference

# LOO and WAIC for mod_aspen
# See ?loo::pareto_k_values for more info on pareto-k warnings
loo(mod_aspen)
waic(mod_aspen)

# Rank candidate models
round(modSel(fitList(mod_null, mod_rand, mod_aspen)),3)

# Parameter estimates and 95% credible intervals
mod_aspen

# Extract a posterior for the effect of aspen
# Look at parameter names
names(mod_aspen)
beta_aspen <- ubms::extract(mod_aspen, "beta_state[scale(pctaspen)]")[[1]]

# Plot posterior distribution
hist(beta_aspen)

# Plot marginal effects and 95% credible intervals of covariates
plot_marginal(mod_aspen, "state")
plot_marginal(mod_aspen, "det")

# Get predicted mean abundance (lambda) estimates for each site-year
abund <- predict(mod_aspen, "state")
head(abund)

# Generate a map of abundance based on percent aspen across a landscape
library(raster)
asp <- raster('pctaspen.tif')
abun_raster <- predict(mod_aspen, "state", newdata=asp, re.form=NA)
plot(abun_raster)

# Get posterior distribution of latent abundance z (niter x nsite-year)
z <- posterior_predict(mod_aspen, "z")
hist(rowMeans(z))

# Look at the Stan code used to fit the model
# Note that because a common Stan file is used to fit several types of models
# (to improve compilation times), there are many parts of the output that
# will not be used in a given model
# Thus this code is much more complicated than what you would need to
# fit your own version of an N-mixture model in Stan
cat(get_stancode(mod_aspen))

# You can use other Stan-related tools directly on the output from rstan,
# which is contained in the @stanfit slot
# For example, model diagnostics with package shinystan:
library(shinystan)
launch_shinystan(mod_aspen@stanfit)

# Illustrate models with multiple random effects are possible
# Random effects of both site and year
# Since we only have 3 unique years this is probably not a good idea, but
# serves as an illustration
mod_year <- stan_pcount(~scale(date) ~(1|year) + (1|site), umf, K=20,
                         chains=3, iter=500, seed=123)
mod_year
