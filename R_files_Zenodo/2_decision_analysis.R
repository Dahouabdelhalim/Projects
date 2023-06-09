################################################################################
### 2_decision_analysis.R
###     Obtain profiles of the species detection effectiveness by using the
###     results of the model fitting.
###
###     The following R packages are used:
###         * foreach: https://CRAN.R-project.org/package=foreach
###         * doSNOW: https://CRAN.R-project.org/package=doSNOW
###         * DirichletReg: https://CRAN.R-project.org/package=DirichletReg
###         * ggplot2: https://CRAN.R-project.org/package=ggplot2
###         * RColorBrewer: https://CRAN.R-project.org/package=RColorBrewer
################################################################################

setwd("<specify your directory>")


### Load functions and posterior samples
source("functions.R")
load("result.Rdata")


### Profiles for local species diversity assessment
# Set constants:
#   B: amount of budget
#   lambda1: cost per sequence read for HTS
#   lambda2: cost per replicate for library preparation
B       <- c(437.5E3, 875E3, 1750E3, 2625E3, 3500E3)
lambda1 <- 0.01
lambda2 <- 5000

# Extract posterior samples
post <- list(z = res$sims.list$z, theta = res$sims.list$theta, phi = res$sims.list$phi)

# Calculate species detection effectiveness; results are output to file
# "result_da1_B_x.Rdata" where x is the value of the budget level.
# Approx. 15 mins to complete.
for (i in seq_along(B)) {
    profile_eutil1(B[i], lambda1, lambda2, post)
}


### Profiles for regional species diversity assessment (Appendix S1)
# Set constants:
#   B: amount of budget
#   lambda1: the cost per sequence read for HTS
#   lambda2: the cost per replicate for library preparation
#   lambda3: the cost per site for visiting
B       <- 1125000
lambda1 <- 0.01
lambda2 <- 5000
lambda3 <- 5000

# Extract posterior samples
post <- list(z = res$sims.list$z, psi = plogis(res$sims.list$gamma[, 1, ]),
             theta = res$sims.list$theta, phi = res$sims.list$phi)

# Calculate species detection effectiveness; results are output to file
# "result_da2_B_x.Rdata" where x is the value of the budget level.
# Approx. 30 mins to complete; increase rep to obtain smoother profiles.
profile_eutil2(B, lambda1, lambda2, lambda3, post, rep = 4)


### Profiles for local species diversity assessment, varying sequence depth (Appendix S2)
# Set constants:
#   B: amount of budget
#   lambda1: the cost per sequence read for HTS
#   lambda2: the cost per replicate for library preparation
#   CV: coefficient of variation of sequence depth
B       <- c(437.5E3, 875E3, 1750E3, 2625E3, 3500E3)
lambda1 <- 0.01
lambda2 <- 5000
CV      <- c(0.5, 1, 2)

# Extract posterior samples
post <- list(z = res$sims.list$z, theta = res$sims.list$theta, phi = res$sims.list$phi)

# Calculate species detection effectiveness; results are output to file
# "result_da1x_CV_x.Rdata" where x is the value of the CV of sequence depth.
# Approx. 25 mins to complete for each CV value; increase repN to obtain smoother profiles.
for (i in seq_along(CV)) {
    profile_eutil1x(B, lambda1, lambda2, post, CV = CV[i], repN = 4)
}


### Profiles for regional species diversity assessment, varying sequence depth (Appendix S2)
# Set constants:
#   B: amount of budget
#   lambda1: the cost per sequence read for HTS
#   lambda2: the cost per replicate for library preparation
#   lambda3: the cost per site for visiting
#   CV: coefficient of variation of sequence depth
B       <- 1125000
lambda1 <- 0.01
lambda2 <- 5000
lambda3 <- 5000
CV <- c(0.5, 1, 2)

# Extract posterior samples
post <- list(z = res$sims.list$z, psi = plogis(res$sims.list$gamma[, 1, ]),
             theta = res$sims.list$theta, phi = res$sims.list$phi)

# Calculate species detection effectiveness; results are output to file
# "result_da2x_CV_x.Rdata" where x is the value of the CV of sequence depth.
# Approx. 30 mins to complete for CV = 0.5 and 1, but longer for CV = 2.
# Increase repN to obtain smoother profiles.
for (i in seq_along(CV)) {
    profile_eutil2x(B, lambda1, lambda2, lambda3, post, CV = CV[i], repN = 4)
}


### ------------------------------------------------------------------------ ###

#source("functions.R")
library(ggplot2)
library(RColorBrewer)


### Draw the profiles for local species diversity assessment
B <- c(437.5E3, 875E3, 1750E3, 2625E3, 3500E3)
load(sprintf("result_da1_B_%s.Rdata", B[1]))

res <- result
for (i in seq_along(B)[-1]) {
    load(sprintf("result_da1_B_%s.Rdata", B[i]))
    res <- rbind(res, result)
}

d   <- as.data.frame(res)
d$B <- as.factor(c(437.5, rep(875, 3), rep(1750, 6), rep(2625, 10), rep(3500, 13)))
colnames(d) <- c("Budget", "lambda1", "lambda2", "replicate", "N", "detection")
g  <- ggplot(d, aes(x = replicate, y = detection, colour = Budget))
g + geom_line() + geom_point() + theme_light() + 
    scale_colour_manual(values = c(pal_rb_or[10], "black", pal_rb_or[3], pal_rb_or[1], "red"))


### Draw the profiles for regional species diversity assessment (Appendix S1)
load("result_da2.Rdata")

res <- data.frame(budget = result[, "B"],
                  site = result[, "J"],
                  replicate = result[, "K"],
                  detection = apply(result[, -c(1:7)], 1, mean, na.rm = TRUE))

d <- as.data.frame(res)
d$replicate <- as.factor(d$replicate)
g <- ggplot(d, aes(x = site, y = detection, colour = replicate))
g + geom_line(size = 0.9) + theme_light() + scale_colour_manual(values = c(pal_rb_or[10:8], pal_rb_or[3:1]))


### Draw profiles for local species diversity assessment, varying sequence depth (Appendix S2)
B  <- c(437.5E3, 875E3, 1750E3, 2625E3, 3500E3)
CV <- c(0.5, 1, 2)

load(sprintf("result_da1_B_%s.Rdata", B[1]))
tmp <- result
for (i in seq_along(B)[-1]) {
    load(sprintf("result_da1_B_%s.Rdata", B[i]))
    tmp <- rbind(tmp, result)
}

res <- data.frame(budget = tmp[, "B"],
                  replicate = tmp[, "K"],
                  CV = 0,
                  detection = tmp[, "U"])

for (i in seq_along(CV)) {
    load(sprintf("result_da1x_CV_%s.Rdata", CV[i]))
    tmp <- data.frame(budget = result[, "B"],
                      replicate = result[, "K"],
                      CV = CV[i],
                      detection = apply(result[, -c(1:6)], 1, mean, na.rm = TRUE))
    res <- merge(res, tmp, all = TRUE)
}

d <- as.data.frame(res)
d$budget <- as.factor(rep(c(437.5, rep(875, 3), rep(1750, 6), rep(2625, 10), rep(3500, 13)), each = 4))
d$CV <- factor(d$CV, levels = c(2, 1, 0.5, 0))
g <- ggplot(d, aes(x = replicate, y = detection, colour = CV, shape = CV))
g + geom_line(size = 0.5) + geom_point(size = 1.2) + theme_light() + facet_grid(~ budget) +
    scale_colour_manual(values = brewer.pal(6, "Dark2")[c(4, 6, 3, 1)]) +
    scale_shape_manual(values = c(16, 15, 17, 18))


### Draw profiles for regional species diversity assessment, varying sequence depth (Appendix S2)
CV <- c(0.5, 1, 2)

load("result_da2.Rdata")
res <- data.frame(budget = result[, "B"],
                  site = result[, "J"],
                  replicate = result[, "K"],
                  CV = 0,
                  detection = apply(result[, -c(1:7)], 1, mean, na.rm = TRUE))

for (i in seq_along(CV)) {
    load(sprintf("result_da2x_CV_%s.Rdata", CV[i]))
    tmp <- data.frame(budget = result[, "B"],
                      site = result[, "J"],
                      replicate = result[, "K"],
                      CV = CV[i],
                      detection = apply(result[, -c(1:8)], 1, mean, na.rm = TRUE))
    res <- merge(res, tmp, all = TRUE)
}

d <- as.data.frame(res)
d$replicate <- as.factor(d$replicate)
d$CV <- factor(d$CV, levels = c(2, 1, 0.5, 0))
g <- ggplot(d, aes(x = site, y = detection, colour = CV))
g + geom_line(size = 0.5) + theme_light() + facet_grid(~ replicate) + scale_colour_manual(values = brewer.pal(6, "Dark2")[c(4, 6, 3, 1)])

