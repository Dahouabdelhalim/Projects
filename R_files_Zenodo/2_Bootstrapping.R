#2: Generating bootstrapped confidence intervals for beta regression model
#Written by Ariane Castagner and Alexander Brenning, finalized in June 2022

# Load required packages:
library("boot")
library("xts")
library("betareg")
library("dplyr")
library("plyr")
library("purrr")
library("rcompanion") # for exploratory stats: Cramer's V
library("pROC")       # for exploratory stats: AUROC
boot_method <- "bca" 

#Set working directory (please change to your path to the project folder)

setwd("...//ExcessIceBetaRegression")

# Import estimate_coefficients function

source("src//estimate_coefficients.R")
source("src//boot_extract.R")
source("src//depthmean.R")
source("src//estimate_depthmean.R")
source("src//estimate_depthmean2.R")
source("src//mywrapper.R")
source("src//estimate_depth_cumul.R")
source("src//mywrapper_cumul.R")

# Import data
data <- read.csv("data//SentinelBHdata.csv")

# Convert variables to factors and set reference levels

data$Material.simplified <- relevel(factor(data$Material.simplified), ref = "Coarse")
data$Deposit.type <- relevel(factor(data$Deposit.type), ref = "Moraine")
data$Deposit.type.conf <- relevel(factor(data$Deposit.type.conf), ref = "2")

# Set up formula for all subsequent analyses:
fo <- EITrans ~ Depth + Material.simplified + Deposit.type + Deposit.type.conf + Mid_VI

# Fit a beta regression on the entire data set:
fit <- betareg(fo, data = data)

# The coefficient names are needed for the estimate_coefficients function:
COEF_NAMES <- names(coef(fit))

###########################################################################
### Test the estimation function using random subsamples of the data
### resampled at the observation or site level:
###########################################################################

# only a subset - you will get an NA somewhere every now and then:
estimate_coefficients(formula = fo, data = data, 
                      x = sample(1:nrow(data), size = floor(nrow(data)/3), replace = TRUE),
                      verbose = 2)
# try sampling at the Site.ID level:
estimate_coefficients(formula = fo, data = unique(data$Site.ID), 
                      actual_data = data, group_by = "Site.ID",
                      x = sample(1:length(unique(data$Site.ID)), replace = TRUE),
                      verbose = 2)



###########################################################################
### Bootstrap estimation of model coefficients
### using ordinary and grouped resampling:
###########################################################################

# Ordinary bootstrap:
set.seed(789)
results_ordinary <- boot(statistic = estimate_coefficients, formula = fo, 
                         data = data, R = 10000)

# Site.ID-level bootstrap (this method was used for the publication):
set.seed(789)
results_grouped <- boot(statistic = estimate_coefficients, formula = fo, 
                        group_by = "Site.ID",
                        data = unique(data$Site.ID), 
                        actual_data = data, R = 10000)

saveRDS(results_grouped, file = "data//results_grouped.rds")

# Instead of running the previous (time-consuming) function for this data, you may import the results:

results_grouped <- readRDS("src//results_grouped.rds")
results <- results_grouped

###########################################################################
### Extract bootstrap confidence intervals at the 95% level
###########################################################################

# collect bootstrap confidence intervals;
# results table also includes coefficient estimates:

coef_estimate <- estimate_coefficients(formula = fo, data = data, x = NULL, verbose = 2)
coef_estimate <- coef_estimate[c(1:9)] # Remove Phi

coef_intv <- data.frame(boot.lwr = rep(NA, length(coef_estimate)),
                        est = coef_estimate,
                        boot.est = colMeans(results$t, na.rm = TRUE),
                        boot.upr = rep(NA, length(coef_estimate)))

for (i in 1:length(coef_estimate))
  coef_intv[i,c("boot.lwr","boot.upr")] <- boot.ci(boot.out = results, type = "basic", 
                                                   index=i)[["basic"]][4:5]
round(coef_intv, 3)

###########################################################################
### Predict EI profiles on large cryostratigraphic dataset
### and create prediction intervals
###########################################################################

# Applying the coefficients to large cryostratigraphic dataset to generate average EI profile with CIs

# Import cryostratigraphic data with simplified materials classification
cryodata <- read.csv("data//CryostratigraphyBHdata.csv")

# There's not enough Deposit.type.conf == 1 in the training set,
# therefore we decided to merge them into Deposit.type.conf == 2:
cryodata$Deposit.type.conf[ cryodata$Deposit.type.conf == 1 ] <- 2
cryodata$Deposit.type.conf <- as.factor(cryodata$Deposit.type.conf)

# Create a clean data.frame for model application:
pred_data <- data.frame(
  Site.ID = cryodata$Borehole.ID,
  From..m. = cryodata$Top.depth,
  To..m. = cryodata$Bottom.depth,
  Depth = cryodata$Depth,
  Material.simplified = cryodata$Material.simplified,
  Deposit.type = cryodata$Deposit.type,
  Deposit.type.conf = cryodata$Deposit.type.conf,
  Mid_VI = cryodata$Mid_VI
)

# backtransform z = (y(n-1)+0.5)/n:
# y = (z*n - 0.5) / (n-1)
backtransform <- function(z, n = nrow(data)) {
  (z*n - 0.5) / (n - 1)
}

pred <- backtransform(predict(fit, newdata = pred_data), n = nrow(data))

# Try out mean prediction at Depth = 1:
depthmean(1, pred_data, fit = fit)

# Try out new estimation function:
estimate_depthmean(formula = fo, data = unique(data$Site.ID), 
                   actual_data = data, group_by = "Site.ID",
                   x = sample(1:length(unique(data$Site.ID)), replace = TRUE),
                   verbose = 2, pred_data = pred_data, Depth = 1)

# Bootstrap the mean excess ice content at Depth = 1
set.seed(789)
pred_grouped <- boot(statistic = estimate_depthmean, formula = fo, 
                     group_by = "Site.ID",
                     data = unique(data$Site.ID), 
                     actual_data = data, 
                     pred_data = pred_data, Depth = 1,
                     R = 1000)

all_pred <- list()
depths <- seq(0.5, 10, by = 0.25)
Nrep <- 1000
for (i in 1:length(depths)) {
  # purrr::quietly() hides all error messages
  all_pred[[i]] <- boot(statistic = estimate_depthmean, formula = fo, 
                        group_by = "Site.ID",
                        data = unique(data$Site.ID), 
                        actual_data = data, 
                        pred_data = pred_data, Depth = depths[i],
                        R = Nrep, verbose = 0)
  # I occasionally (rarely) get an 'optim' error, which I assume is a non-convergence
  # in the betareg fitting process due to some extravagant training sample
  # -> nothing to worry about, given the low frequency of these exceptions
  # (like one in a thousand models...).
}
names(all_pred) <- as.character(depths)


## Parallelizing this time-consuming loop:
## allpred <- lapply(depths, mywrapper)

Nrep = 100
depths <- seq(.5, 10, by = .5)

library("parallel")
cl <- parallel::makeCluster(getOption("cl.cores", 2))
clusterExport(cl, c("estimate_coefficients", "mywrapper", "estimate_depthmean",
                    "estimate_depthmean2", "backtransform",
                    "fo", "data", "pred_data", "Nrep", "COEF_NAMES",
                    "depths"))
#parallel::clusterCall(cl = cl, fun = mywrapper, depth = 1)
all_pred <- parallel::parLapply(cl, X = depths, mywrapper)
stopCluster(cl)

save(all_pred, Nrep, depths, file = "output//data//depth_curve.Rdata")

# Or it can be imported:
# (load("src//depth_curve.Rdata"))

# Compile bootstrap results and model predictions into a data.frame:
pred_boot <- as.data.frame(bind_rows(lapply(all_pred, boot_extract)))

pred_boot$estimate <- sapply(depths, depthmean, data = pred_data, fit = fit)
pred_boot$Depth <- depths

# Correct negative lower CI values to 0 if needed:
pred_boot$boot.corr.lwr <- pred_boot$boot.lwr %>% purrr::map_dbl(~max(.x,0))

# Export pred_boot

write.csv(pred_boot, "output//data//pred_boot.csv")

# Simple plot of results and CIs
par(mfrow = c(1,1))
plot(Depth ~ estimate, data = pred_boot, type = "l", col = "blue", lwd = 2, 
     xlim = c(0, 0.35), ylim = c(10,0),
     xlab = "Volumetric excess ice content")
lines(Depth ~ boot.mean, data = pred_boot, col = "blue", lty = "dotted") # slight bias may be related to asymmetric error distribution
lines(Depth ~ boot.corr.lwr, data = pred_boot, col = "blue", lty = "dashed")
lines(Depth ~ boot.upr, data = pred_boot, col = "blue", lty = "dashed")


###########################################################################
### Predict CUMULATIVE EI profiles on large cryostratigraphic dataset
### and create prediction intervals
###########################################################################

Nrep = 1000
maxdepths <- seq(1, 10, by = 1)

library("parallel")
cl <- parallel::makeCluster(getOption("cl.cores", 2))
clusterExport(cl, c("estimate_coefficients", "mywrapper_cumul", 
                    "estimate_depthmean2",
                    "backtransform", "estimate_depth_cumul",
                    "fo", "data", "pred_data", "Nrep", "COEF_NAMES",
                    "maxdepths"))
all_pred <- parallel::parLapply(cl, X = maxdepths, mywrapper_cumul)
stopCluster(cl)

save(all_pred, Nrep, maxdepths, file = "output//data//depth_cumul.Rdata")

# Or it can be imported:
# (load("src//depth_cumul.Rdata"))


# Compile bootstrap results and model predictions into a data.frame:
pred_cum <- as.data.frame(bind_rows(lapply(all_pred, boot_extract)))
pred_cum$estimate <- sapply(maxdepths, estimate_depth_cumul, data = data, x = NULL, from_Depth = 0.5, step_Depth = 0.1,
                             pred_data = pred_data, formula = fo)
pred_cum$Depth <- maxdepths

pred_cum

# Change negative values to 0 if needed:
pred_cum$boot.corr.lwr <- pred_cum$boot.lwr %>% purrr::map_dbl(~max(.x,0))

# Export cumulative pred_boot (pred_cum):
write.csv(pred_cum, "output//data//pred_cum.csv")

# Simple plot of cumulative results and CIs
par(mfrow = c(1,1))
plot(Depth ~ estimate, data = pred_boot, type = "l", col = "blue", lwd = 2, 
     xlim = c(0, 3), ylim = c(10,0),
     xlab = "Cumulative volumetric excess ice content")
lines(Depth ~ boot.mean, data = pred_boot, col = "blue", lty = "dotted") # slight bias may be related to asymmetric error distribution
lines(Depth ~ boot.corr.lwr, data = pred_boot, col = "blue", lty = "dashed")
lines(Depth ~ boot.upr, data = pred_boot, col = "blue", lty = "dashed")


