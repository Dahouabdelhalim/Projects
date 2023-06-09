library(ubms)
data(crossbill) # dataset comes from unmarked package
?crossbill      # info on dataset

# We'll fit a single-season model, so we need to select a subset of
# the multi-year dataset
y <- crossbill[,c("det021","det022","det023")] # presence/absence for 1 year
sc <- crossbill[,c("ele","forest")]            # site covariates

# Create unmarkedFrame
umf <- unmarkedFrameOccu(y=y, siteCovs=sc)
head(umf)

# Fit series of models

# Null model
mod_null <- stan_occu(~1 ~1, umf, chains=3, iter=1000)

# Model with date effect on detection and elevation on occupancy
# NOTE: should always scale continuous covariates
# ubms will warn you if it thinks you haven't
mod_ele <- stan_occu(~1 ~scale(ele), umf, chains=3, iter=1000)

# Date and forest
mod_for <- stan_occu(~1 ~scale(forest), umf, chains=3, iter=1000)

# Date, forest, elevation
mod_both <- stan_occu(~1 ~scale(ele)+scale(forest), umf, chains=3, iter=1000)

# Combine models into fitList
fl <- fitList(mod_null, mod_ele, mod_for, mod_both)

# Model selection
# Focus on elpd_diff and se_diff columns
# models where abs(elpd_diff) >> se_diff have less predictive power than top model
# In this case it looks like mod_both and mod_for are basically equivalent
# but both are better than mod_ele and mod_null
round(modSel(fl),3)

# Compare to unmarked
unm_null <- occu(~1~1, umf)
unm_ele <- occu(~1~scale(ele), umf)
unm_for <- occu(~1~scale(forest), umf)
unm_both <- occu(~1~scale(ele)+scale(forest), umf)

# Same results
modSel(fitList(null=unm_null, ele=unm_ele, forest=unm_for, both=unm_both))

# Look at traceplots
traceplot(mod_both) # look good

# Look at posteriors
plot_posteriors(mod_both)
plot_posteriors(mod_both, density=TRUE)

# Check state residuals of top model
# (no detection covariates so can't make residuals for that submodel)
plot_residuals(mod_both, "state")

# Check goodness of fit, posterior predictive check with MacKenzie-Bailey chi-square
(g <- gof(mod_both))
plot(g) # Doesn't look like the fit is good as P should be ~ 0.5

# Look at output of top model
# Both covariates have positive effect on occupancy of crossbill
mod_both

# Compare to frequentist estimates from unmarked; similar
(occu(~1~scale(ele)+scale(forest), umf))

# Look at marginal effects plot
plot_marginal(mod_both, "state")

# Predict occupancy for a new site
newdata <- data.frame(ele=400, forest=50)
predict(mod_both, 'state', newdata=newdata)

# Predict occupancy for a raster map
# Adapted from code written by Richard Chandler

# Spatial landscape data corresponding with crossbill dataset
data(Switzerland)

# Convert to RasterStack and project
library(raster)
ele <- rasterFromXYZ(Switzerland[,c("x","y","elevation")],
    crs="+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333
    +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel
    +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")
names(ele) <- "ele"

forest <- rasterFromXYZ(Switzerland[,c("x","y","forest")],
    crs="+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333
    +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel
    +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")

# Create stack; double check that covariate names match names in ubms model
(spatial_data <- stack(ele, forest))
plot(spatial_data)

# Generate predicted occupancy map
crossbill_map <- predict(mod_both, "state", newdata=spatial_data)
plot(crossbill_map)
