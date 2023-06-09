## Marta Shocket, Indiana University & UCLA
## Code for Can hot temperatures limit disease transmission? A test of mechanisms in a zooplankton-fungus system in Functional Ecology.
##
## Purpose: Analyze foraging rate at 30C and transmission rate & spore infectivity from beta + u measurement assay where
##          (1) hosts were exposed to spores and (2) infections were allowed to establish in constant 20C or cycling 32C conditions
##          (16:8 hour cycle of 32 degree C/20 degree C). Four treatments (exposure/infection establishment): 20/20, 20/32, 32/20, & 32/32C. 
##
## Table of Contents:
##      1. Set working directory, load libraries and data
##      2. Calculate transmission rate (beta) for each temperature treatment
##      3. Bootstrap transmission rate (beta) for each temperature treatment
##      4. Calculating likelihood profiles for beta in 3232 and 3220 treatments (& replacing 100% infected bootstraps)
##      5. Calculate and bootstrap foraging rate (f) at 30C
##      6. Calculate and boostrap spore infectivity (u) for each temperature treatment
##      7. Randomization tests comparing beta and u from different temperature treatments
##      8. Figure 2


######
##### 1. Set working directory, load libraries and data
######

# Set wd
setwd("~/Dropbox/Research Hall Lab/Temperature/Upper Thermal Constraints/Final Code and Data/For Dryad")

# Load libraries
library(bbmle)

# Get Beta + u measurement assay, 30C foraging rate assay data, and bootstrapped foraging rate 
# from Shocket et al. 2018 in Ecology ("Parasite Rearing and Infection Temperatures Jointly Influence...")
beta.data <- read.csv("UTC_BetaUAssay_Prev.csv")
f.data <- read.csv("UTC_HighTempForagingAssay.csv")
bs.f.Shocket2018 <- read.csv("MLE_f_points.csv")

## Subset transmission rate data by temperature treatment
beta.data.2020 <- subset(beta.data, Treatment == "20/20")
beta.data.2032 <- subset(beta.data, Treatment == "20/32")
beta.data.3220 <- subset(beta.data, Treatment == "32/20")
beta.data.3232 <- subset(beta.data, Treatment == "32/32")


######
##### 2. Calculate transmission rate for each temperature treatment
######

# Hosts were exposed to 45 sp / ml in 15 ml for 1 day
# Initially fit in beaker units, but needs to be transformed to L later

beta.2020 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*45*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = beta.data.2020, method="Brent", lower = 0, upper = .1)

beta.2032 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*45*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = beta.data.2032, method="Brent", lower = 0, upper = .1)

beta.3220 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*45*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = beta.data.3220, method="Brent", lower = 0, upper = .1)

beta.3232 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*45*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = beta.data.3232, method="Brent", lower = 0, upper = .1)


beta.u.out <- data.frame(exp.temp = c(20, 20, 32, 32), dev.temp = c(20, 32, 20, 32), 
                       beta.calc = numeric(4), beta.lower = numeric(4), beta.upper = numeric(4),
                       u.calc = numeric(4), u.lower = numeric(4), u.upper = numeric(4))
row.names(beta.u.out) <- c("20/20","20/32", "32/20", "32/32")

# Store estimates in output dataframe - transform into correct units by dividing by host density = 1/0.015 L 
beta.u.out$beta.calc[1] <- coef(beta.2020)/(1/0.015)
beta.u.out$beta.calc[2] <- coef(beta.2032)/(1/0.015)
beta.u.out$beta.calc[3] <- coef(beta.3220)/(1/0.015)
beta.u.out$beta.calc[4] <- coef(beta.3232)/(1/0.015)


######
##### 3. Bootstrap transmission rate for each temperature treatment
######

# Create data sets for sampling for bootstrapping
bsdata.2020.U <- rep(0, 16)
bsdata.2020.I <- rep(1, 30)
bsdata.2020 <- c(bsdata.2020.U, bsdata.2020.I)

bsdata.2032.U <- rep(0, 25)
bsdata.2032.I <- rep(1, 24)
bsdata.2032 <- c(bsdata.2032.U, bsdata.2032.I)

bsdata.3220.U <- rep(0, 3)
bsdata.3220.I <- rep(1, 42)
bsdata.3220 <- c(bsdata.3220.U, bsdata.3220.I)

bsdata.3232.U <- rep(0, 6)
bsdata.3232.I <- rep(1, 41)
bsdata.3232 <- c(bsdata.3232.U, bsdata.3232.I)

################## Function - BootstrapBeta - needs Vol, Spore.Dose, Time, and Density loaded separately
# Arguments: 
  # data.set          vector: 0/1 infection data
  # nBoots            integer: number of bootstraping iterations
# Returns: 
  # function.coefs    data frame: bootstrapped fit coefficients

BootstrapBeta = function(data.set, nBoots) {
  
  # Create data frame for output, simulated dataset
  output = data.frame(beta = numeric(length(nBoots)), prev = numeric(length(nBoots)))
  
  # Loop through iterations
  for (i in 1:nBoots) {
    
    # Sample data
    sample.list <- sample(data.set, length(data.set), replace = TRUE)
    
    # Calculate number of uninfected and total from sample
    simulated.dataset = data.frame(uninf = numeric(1), total = numeric(1))
    simulated.dataset[1, 1] <- length(sample.list) - sum(sample.list)
    simulated.dataset[1, 2] <- length(sample.list)
    
    # Calculate beta of simulated data set
    beta <- mle2(uninf ~ dbinom(size=total, prob=exp(-Beta*Vol*Spore.Dose*Time)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = simulated.dataset, method="Brent", lower = 0, upper = .1)
    
    # Store beta parameter and prev: divide by host density (1/L) to get correct units (L instead of tube)
    output[i, 1] = coef(beta) / Density # beta
    output[i, 2] = (simulated.dataset[1,2] - simulated.dataset[1,1]) / simulated.dataset[1,2] # prevalence
    
  }
  
  # Return output
  output
}

# Load constants for this data set
Vol <- 15 # ml
Spore.Dose <- 45 # sp/ml
Time <- 1 # day
Density <- 1/0.015 # hosts/L

# Bootstrap betas
bs.beta.2020 <- BootstrapBeta(bsdata.2020, 10000)
bs.beta.2032 <- BootstrapBeta(bsdata.2032, 10000)
bs.beta.3220 <- BootstrapBeta(bsdata.3220, 10000)
bs.beta.3232 <- BootstrapBeta(bsdata.3232, 10000)

# Calculate and store quantiles in output dataframe
beta.u.out$beta.lower[1] <- quantile(bs.beta.2020$beta, probs = c(0.025))
beta.u.out$beta.upper[1] <- quantile(bs.beta.2020$beta, probs = c(0.975))
beta.u.out$beta.lower[2] <- quantile(bs.beta.2032$beta, probs = c(0.025))
beta.u.out$beta.upper[2] <- quantile(bs.beta.2032$beta, probs = c(0.975))
beta.u.out$beta.lower[3] <- quantile(bs.beta.3220$beta, probs = c(0.025))
beta.u.out$beta.upper[3] <- quantile(bs.beta.3220$beta, probs = c(0.975))
beta.u.out$beta.lower[4] <- quantile(bs.beta.3232$beta, probs = c(0.025))
beta.u.out$beta.upper[4] <- quantile(bs.beta.3232$beta, probs = c(0.975))


#####
#### 4. Calculating likelihood profiles for beta in 3232 and 3220 treatments (& replacing 100% infected bootstraps)
#####

# The bootstrapping method in section 3 works for intermediate levels of infection, but doesn't accurately capture
# the upper confidence interval for beta when prevalence is very high and you sometimes sample 100% prevalence
# giving an inflated value for beta = 0.0015 (as occurs in the 32/32 and 32/20 treatments where only 6/47 and 
# 3/45 animals remained uninfected, respectively).

# To better capture the upper interval, we correct the beta value for these 100% prevalence samples using the 
# likelihood profile method. This approach calculates the likelihood profile and finds the beta value where 
# the loglikelihood is 1.92 units lower than our estimated value, the highest reasonable value according to
# the MLE. Then we substitute it in our bootstrapped samples that by chance had 100% prevalence. (1.92 is half 
# the 3.84 chi-square value with 1 degree of freedom.) See Cole et al. 2013 American Journal of Epidemiology. 
# "Maximum Likelihood, Profile Likelihood, and Penalized Likelihood: A Primer."

# Step 1: Calculate the best parameter for the probability of remaining uninfected (p) for the binomial distribution.
# The untransformed beta for 3232 = 0.003049464, for 3220 = 0.004011926 (from the mle2() estimates above).
p.3232 = exp(-45*15*0.003049464) # this value is also the raw prevalence = 6/47
p.3220 = exp(-45*15*0.004011926) # this value is also the raw prevalence = 3/45

# Calculate log-likelihoods according to the binomial distribution: log(p^i * q^(n-i)),
# where q = 1-p, i is the # of successes (non-infections), and n is the # of trials
log(p.3232^6*(1-p.3232)^41) # the maximum loglikelihood
log(p.3232^6*(1-p.3232)^41)-1.92 # the threshold value 1.92 units away

log(p.3220^3*(1-p.3220)^42) # the maximum loglikelihood
log(p.3220^3*(1-p.3220)^42)-1.92 # the threshold value 1.92 units away

########### Iteratively solve for the p > the best estimate that gives you the threshold value - 3220
# (Keep adding decimals points to the list until you zero in on the threshold value)
list3220 <- seq(0.0060351, 0.0060352, 0.00000001)
nll.3220 <- log(exp(-45*15*list3220)^3*(1-exp(-45*15*list3220))^42)
nll.3220
# Answer for 3220 is untransformed beta = 0.00603511 (beaker units)
# 0.00603511*0.015 = 9.052665e-05 (transformed into L units)

########### Iteratively solve for the p > the best estimate that gives you the threshold value - 3232
# (Keep adding decimals points until you zero in on the threshold value)
list3232 <- seq(0.00435676, 0.0043568, 0.00000001)
nll.3232 <- log(exp(-45*15*list3232)^6*(1-exp(-45*15*list3232))^41)
nll.3232
# Answer for 3232 is untransformed beta = 0.00435677 (beaker units)
# 0.00435677*0.015 = 6.535155e-05 (transformed into L units)

########### Replacing values in bootstrapped samples
# To be consistent, we'll replace 100% prevalence bootstraps in both treatments with the same value (9.052665e-05).
# It's the MLE estimated upper CI for the 3220 treatment (the higher prevalence of the two treatments), and it's 
# greater than all of the calculated beta values (the highest is prev = 0.9777778 -> beta = 8.459246e-05). If we
# used the 6.535155e-05 value for the 3232 treatment, it would actually be lower than the highest reasonable 
# calculated beta value (prev = 0.9787234 -> beta = 8.555884e-05).

# historgrams show small peak at 0.0015, and all other values below 10^-4
hist(bs.beta.3220$beta)
hist(bs.beta.3232$beta)

# Replace values > 0.001 with the likelihood-based maximum beta (for some reason R isn't recognizing the exact value)
bs.3220.rep <- replace(bs.beta.3220$beta, bs.beta.3220$beta > 0.001, 9.052665e-05)
bs.3232.rep <- replace(bs.beta.3232$beta, bs.beta.3232$beta > 0.001, 9.052665e-05)

# now histograms look reasonable
hist(bs.3220.rep) 
hist(bs.3232.rep)

# Replace the columns in the original data frame
bs.beta.3220$beta <- bs.3220.rep
bs.beta.3232$beta <- bs.3232.rep

# Replace the upper CI estimates for the 3220 treatment (3232 was already okay b/c 100% prev samples < 2.5%,
# but we still had to replace the values in 3232 because we'll use them to calculate transmission potential later.).
beta.u.out$beta.upper[3] <- quantile(bs.beta.3220$beta, probs = c(0.975))

bs.beta <- data.frame(beta.2020 = bs.beta.2020$beta, beta.2032 = bs.beta.2032$beta, 
                      beta.3220 = bs.beta.3220$beta, beta.3232 = bs.beta.3232$beta)

save(beta.u.out, file = "beta.u.out.Rsave")
save(bs.beta, file = "bs.beta.Rsave")
load("bs.beta.Rsave")

######
##### 5. Calculate and bootstrap foraging rate (f) at 30C
######

# Best estimate for foraging rate at 30 C
f.30 <- mle2(f ~ dnorm(mean = f.hat * (Size ^ n), sd = sd),
              start = list(f.hat = 0.005, n = 2, sd = 0.01), data = f.data)
summary(f.30)

################## Function - BootstrapF 
# Arguments: 
  # data.set     data frame: foraging rate data set
  # nBoots       integer: number of bootstraping iterations
# Returns: 
  # f.out        data frame: bootstrapped fit coefficients

BootstrapF = function(data.set, nBoots){
  
  f.out <- data.frame(f.hat = numeric(nBoots), n = numeric(nBoots), sd = numeric(nBoots), adult.f = numeric(nBoots))
  
  for (i in 1:nBoots){
    
    # Create a sample of random row numbers
    boot.rows = sample(nrow(data.set), size=nrow(data.set), replace=TRUE)
    
    # Create a data frame and fill it with the first animal from random sample of row numbers
    simulated.data.set = data.set[boot.rows[1], ]
    
    #Finish filling in the data frame with the rest of the animals from random sample of row numbers
    for (iRow in 2:length(boot.rows)){
      simulated.data.set = rbind(simulated.data.set, data.set[boot.rows[iRow], ])
    }      
    
    f.fit <- mle2(f ~ dnorm(mean = f.hat * (Size ^ n), sd = sd),
                 start = list(f.hat = 0.005, n = 2, sd = 0.01), data = simulated.data.set)
    
    f.out[i, 1] <- coef(f.fit)[1]
    f.out[i, 2] <- coef(f.fit)[2]
    f.out[i, 3] <- coef(f.fit)[3]
    f.out[i, 4] <- coef(f.fit)[1]*1.5^coef(f.fit)[2]
    
  }

  f.out # return output
  
}

# Bootstrap f at 30C
bs.f.30 <- BootstrapF(f.data, 10000)

quantile(bs.f.30$n, prob = c(0.025, 0.975))
par(mfrow = c(1,1))
hist(bs.f.30$n)

# Create dataframe of bootstrapped foraging rate for all three temps
bs.f <- data.frame(f.20 = bs.f.Shocket2018$fr20, f.25 = bs.f.Shocket2018$fr25, f.30 = bs.f.30$adult.f)

# Calculate p-values by comparing best estimates to cumulative probability distributions of bootstrapped values of different temperatures
#    Note: For correct p-value use the actual cumulative probability density if the distribution mean is > the estimate
#          but use "1 -" the cumulative probability density if the estimate is > the distribution mean
1 - ecdf(bs.f$f.20)(f.out$f.calc[2]) # 20 vs. 25
ecdf(bs.f$f.25)(f.out$f.calc[1]) # 20 vs. 25

1 - ecdf(bs.f$f.20)(f.out$f.calc[3]) # 20 vs. 30
ecdf(bs.f$f.30)(f.out$f.calc[1]) # 20 vs. 30

ecdf(bs.f$f.25)(f.out$f.calc[3]) # 25 vs. 30
1 - ecdf(bs.f$f.30)(f.out$f.calc[2]) # 30 vs. 25


# Create output dataframe for f
f.out <- data.frame(temp = c(20, 25, 30), f.calc = numeric(3), f.upper = numeric(3), f.lower = numeric(3))

# Store the best estimate and quantiles for foraging rate for 1.5 mm adults at 30C
f.out$f.calc[3] <- coef(f.30)[1]*1.5^coef(f.30)[2]
f.out$f.lower[3] <- quantile(bs.f.30$adult.f, 0.025)
f.out$f.upper[3] <- quantile(bs.f.30$adult.f, 0.975)

# Do the same for foraging rate for 1.5 mm adults at 20 and 25 C from Shocket et al. 2018 Ecology
f.out$f.calc[1] <- mean(bs.f$f.20)
f.out$f.lower[1] <- quantile(bs.f$f.20, 0.025)
f.out$f.upper[1] <- quantile(bs.f$f.20, 0.975)

f.out$f.calc[2] <- mean(bs.f$f.25)
f.out$f.lower[2] <- quantile(bs.f$f.25, 0.025)
f.out$f.upper[2] <- quantile(bs.f$f.25, 0.975)

save(bs.f, file = "bs.f.Rsave")
save(f.out, file = "f.out.Rsave")
load("bs.f.Rsave")
load("f.out.Rsave")


######
##### 6. Calculate and boostrap spore infectivity (u) for each temperature treatment
######

bs.u <- data.frame("u.2020" = numeric(10000), "u.2032" = numeric(10000),"u.3220" = numeric(10000),"u.3232" = numeric(10000))

# Calculate bootstrapped u from bootstrapped beta and f
bs.u$u.2020 <- bs.beta$beta.2020 / bs.f$f.20
bs.u$u.2032 <- bs.beta$beta.2032 / bs.f$f.20
bs.u$u.3220 <- bs.beta$beta.3220 / bs.f$f.30
bs.u$u.3232 <- bs.beta$beta.3232 / bs.f$f.30

# Store best estimate in outputdataframe
beta.u.out$u.calc[1] <- beta.u.out$beta.calc[1] / f.out$f.calc[1]
beta.u.out$u.calc[2] <- beta.u.out$beta.calc[2] / f.out$f.calc[1]
beta.u.out$u.calc[3] <- beta.u.out$beta.calc[3] / f.out$f.calc[3]
beta.u.out$u.calc[4] <- beta.u.out$beta.calc[4] / f.out$f.calc[3]

# Store upper and lower quantiles in output dataframe
beta.u.out$u.lower[1] <- quantile(bs.u$u.2020, 0.025)
beta.u.out$u.lower[2] <- quantile(bs.u$u.2032, 0.025)
beta.u.out$u.lower[3] <- quantile(bs.u$u.3220, 0.025)
beta.u.out$u.lower[4] <- quantile(bs.u$u.3232, 0.025)

beta.u.out$u.upper[1] <- quantile(bs.u$u.2020, 0.975)
beta.u.out$u.upper[2] <- quantile(bs.u$u.2032, 0.975)
beta.u.out$u.upper[3] <- quantile(bs.u$u.3220, 0.975)
beta.u.out$u.upper[4] <- quantile(bs.u$u.3232, 0.975)

save(beta.u.out, file = "beta.u.out.Rsave")
save(bs.u, file = "bs.u.Rsave")


######
##### 7. Randomization tests comparing beta and u from different temperature treatments
######

##### Make data frames for all the randomization comparisons (Need to load 0/1 infection datasets in section 3)
# 20C exposures (2020 vs. 2032)
exps20.inf <- c(bsdata.2020, bsdata.2032)
exps20.treat <- c(rep(2020, length(bsdata.2020)), rep(2032, length(bsdata.2032)))
exps20.f <- c(rep(f.out$f.calc[1], length(bsdata.2020)), rep(f.out$f.calc[1], length(bsdata.2032)))
exps20.df <- data.frame(infect = exps20.inf, Treatment = exps20.treat, f = exps20.f)

# 32C exposures (3220 vs. 3232)
exps32.inf <- c(bsdata.3220, bsdata.3232)
exps32.treat <- c(rep(3220, length(bsdata.3220)), rep(3232, length(bsdata.3232)))
exps32.f <- c(rep(f.out$f.calc[3], length(bsdata.3220)), rep(f.out$f.calc[3], length(bsdata.3232)))
exps32.df <- data.frame(infect = exps32.inf, Treatment = exps32.treat, f = exps32.f)

# 20C infection developments (2020 vs. 3220)
devs20.inf <- c(bsdata.2020, bsdata.3220)
devs20.treat <- c(rep(2020, length(bsdata.2020)), rep(3220, length(bsdata.3220)))
devs20.f <- c(rep(f.out$f.calc[1], length(bsdata.2020)), rep(f.out$f.calc[3], length(bsdata.3220)))
devs20.df <- data.frame(infect = devs20.inf, Treatment = devs20.treat, f = devs20.f)

# 32C infection developments (2032 vs. 3232)
devs32.inf <- c(bsdata.2032, bsdata.3232)
devs32.treat <- c(rep(2032, length(bsdata.2032)), rep(3232, length(bsdata.3232)))
devs32.f <- c(rep(f.out$f.calc[1], length(bsdata.2032)), rep(f.out$f.calc[3], length(bsdata.3232)))
devs32.df <- data.frame(infect = devs32.inf, Treatment = devs32.treat, f = devs32.f)

# Diagonal (exposure =  infection development) (2020 vs. 3232)
diag.inf <- c(bsdata.2020, bsdata.3232)
diag.treat <- c(rep(2020, length(bsdata.2020)), rep(3232, length(bsdata.3232)))
diag.f <- c(rep(f.out$f.calc[1], length(bsdata.2020)), rep(f.out$f.calc[3], length(bsdata.3232)))
diag.df <- data.frame(infect = diag.inf, Treatment = diag.treat, fr = diag.f)


################## Function - BetaUDiffRand
# Arguments: 
#   data.set          dataframe: 0/1 infection data
#   nBoots            integer: number of randomizations
# Returns: 
#   betaudiff.out      dataframe:  bootstrapped fit coefficients

################## Function - NOTE: This version replaces 100% infected samples with highest value from MLE profiling
BetaUDiffRand = function(data.set, treat.rand.1, treat.rand.2, nBoots) {
  
  # Create data frame for output
  betaudiff.out = data.frame(beta.1 = numeric(nBoots), beta.2 = numeric(nBoots), beta.diff = numeric(nBoots),
                            u.1 = numeric(nBoots),  u.2 = numeric(nBoots),  u.diff = numeric(nBoots))
  
  # Loop through iterations
  for (i in 1:nBoots) {
    
    # Add randomized dates to data frame
    simulated.dataset <- data.set
    simulated.dataset$Treat.rand <- sample(simulated.dataset$Treatment, replace = FALSE)    
    
    # Subset data frames
    simulated.dataset.treat.1 <- subset(simulated.dataset, Treat.rand == treat.rand.1)
    simulated.dataset.treat.2 <- subset(simulated.dataset, Treat.rand == treat.rand.2)
    
    # Flatten each date for mle calculation
    sim.dataset.treat.1.flat <- data.frame(Uninf = numeric(1), Total = numeric(1))
    sim.dataset.treat.1.flat[1, 1] <- nrow(simulated.dataset.treat.1) - sum(simulated.dataset.treat.1$infect)
    sim.dataset.treat.1.flat[1, 2] <- nrow(simulated.dataset.treat.1)
    
    sim.dataset.treat.2.flat <- data.frame(Uninf = numeric(1), Total = numeric(1))
    sim.dataset.treat.2.flat[1, 1] <- nrow(simulated.dataset.treat.2) - sum(simulated.dataset.treat.2$infect)
    sim.dataset.treat.2.flat[1, 2] <- nrow(simulated.dataset.treat.2)
    
    # Calculate betas and mean f of randomized data sets
    beta.treat.1 <- mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*Vol*Spore.Dose*1)), start=list(Beta = 0.00004), 
                         control=list(parscale = c(Beta = 0.00001)), data = sim.dataset.treat.1.flat, method="Brent", lower = 0, upper = .1)
    beta.treat.2 <- mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*Vol*Spore.Dose*1)), start=list(Beta = 0.00004), 
                         control=list(parscale = c(Beta = 0.00001)), data = sim.dataset.treat.2.flat, method="Brent", lower = 0, upper = .1)
    f.treat.1 <- mean(simulated.dataset.treat.1$f)
    f.treat.2 <- mean(simulated.dataset.treat.2$f)
    
    # Calculate beta parameters - divide by host density (1/L) to get correct units (L instead of tube)
    cf.beta.treat.1 <- coef(beta.treat.1)/Density
    cf.beta.treat.2 <- coef(beta.treat.2)/Density
    # Replace 100% infected samples with highest beta value from likelihood profile
    if(cf.beta.treat.1 > 1e-04) cf.beta.treat.1 = 9.054065e-05
    if(cf.beta.treat.2 > 1e-04) cf.beta.treat.2 = 9.054065e-05
    
    # Store betas/us and calculate the absolute value of their differences (so positive & negative extreme values both count)
    betaudiff.out$beta.1[i] = cf.beta.treat.1
    betaudiff.out$beta.2[i] = cf.beta.treat.2
    betaudiff.out$beta.diff[i] = abs(betaudiff.out$beta.1[i] - betaudiff.out$beta.2[i])
    betaudiff.out$u.1[i] = cf.beta.treat.1 / f.treat.1
    betaudiff.out$u.2[i] = cf.beta.treat.2 / f.treat.2
    betaudiff.out$u.diff[i] = abs(betaudiff.out$u.1[i] - betaudiff.out$u.2[i])
    
  }
  
  betaudiff.out # Return output
  
}


# Run the randomizations
beta.u.diff.20exp <- BetaUDiffRand(exps20.df, 2020, 2032, 10000)
beta.u.diff.32exp <- BetaUDiffRand(exps32.df, 3220, 3232, 10000)
beta.u.diff.20dev <- BetaUDiffRand(devs20.df, 2020, 3220, 10000)
beta.u.diff.32dev <- BetaUDiffRand(devs32.df, 2032, 3232, 10000)
beta.u.diff.diag <- BetaUDiffRand(diag.df, 2020, 3232, 10000)

# Check the output
hist(beta.u.diff.20exp$beta.diff)
hist(beta.u.diff.20dev$beta.diff)
hist(beta.u.diff.32exp$beta.diff)
hist(beta.u.diff.32dev$beta.diff)
hist(beta.u.diff.diag$beta.diff)

hist(beta.u.diff.20exp$u.diff)
hist(beta.u.diff.20dev$u.diff)
hist(beta.u.diff.32exp$u.diff)
hist(beta.u.diff.32dev$u.diff)
hist(beta.u.diff.diag$u.diff)

# Collect the differences in one data frame
beta.u.diff.df <- data.frame(beta.exp20 = beta.u.diff.20exp$beta.diff, beta.exp32 = beta.u.diff.32exp$beta.diff, 
                           beta.dev20 = beta.u.diff.20dev$beta.diff, beta.dev32 = beta.u.diff.32dev$beta.diff,
                           beta.diag = beta.u.diff.diag$beta.diff,
                           u.exp20 = beta.u.diff.20exp$u.diff, u.exp32 = beta.u.diff.32exp$u.diff, 
                           u.dev20 = beta.u.diff.20dev$u.diff, u.dev32 = beta.u.diff.32dev$u.diff,
                           u.diag = beta.u.diff.diag$u.diff)

# Calculate differences between best estimates of beta
bd.exp20 <- beta.u.out$beta.calc[1] - beta.u.out$beta.calc[2]
bd.exp32 <- beta.u.out$beta.calc[3] - beta.u.out$beta.calc[4]
bd.dev20 <- beta.u.out$beta.calc[3] - beta.u.out$beta.calc[1]
bd.dev32 <- beta.u.out$beta.calc[4] - beta.u.out$beta.calc[2]
bd.diag <- beta.u.out$beta.calc[4] - beta.u.out$beta.calc[1]

# Calculate differences between best estimates of u
ud.exp20 <- beta.u.out$u.calc[1] - beta.u.out$u.calc[2]
ud.exp32 <- beta.u.out$u.calc[3] - beta.u.out$u.calc[4]
ud.dev20 <- beta.u.out$u.calc[3] - beta.u.out$u.calc[1]
ud.dev32 <- beta.u.out$u.calc[4] - beta.u.out$u.calc[2]
ud.diag <- beta.u.out$u.calc[4] - beta.u.out$u.calc[1]

# Calculate p-values for beta randomization tests from the cumulative probability distribution density
1 - ecdf(beta.u.diff.df$beta.exp20)(bd.exp20)
1 - ecdf(beta.u.diff.df$beta.exp32)(bd.exp32)
1 - ecdf(beta.u.diff.df$beta.dev20)(bd.dev20)
1 - ecdf(beta.u.diff.df$beta.dev32)(bd.dev32)
1 - ecdf(beta.u.diff.df$beta.diag)(bd.diag)

# Calculate p-values for u randomization tests from the cumulative probability distribution density
1 - ecdf(beta.u.diff.df$u.exp20)(ud.exp20)
1 - ecdf(beta.u.diff.df$u.exp32)(ud.exp32)
1 - ecdf(beta.u.diff.df$u.dev20)(ud.dev20)
1 - ecdf(beta.u.diff.df$u.dev32)(ud.dev32)
1 - ecdf(beta.u.diff.df$u.diag)(ud.diag)

save(beta.u.diff.df, file = "beta.u.diff.df.Rsave")
load("beta.u.diff.df.Rsave")


######
##### 8. Figure 2
######

# Load saved results for plotting
load("f.out.Rsave")
load("beta.u.out.Rsave")

exp.temp.offset <- c(-.25, .25, -.25, .25)
beta.u.out$exp.temp <- beta.u.out$exp.temp + exp.temp.offset
beta.u.out.byexp <- rbind(beta.u.out[1, ], beta.u.out[3, ], beta.u.out[2, ], beta.u.out[4, ])
beta.u.out <- beta.u.out.byexp 

## plot settings
par(mfrow = c(3,1), las = 1, oma = c(0, 0, 0, 0), mar = c(5, 5, 1, 1))

##### Panel A: beta
plot(beta.calc*10000 ~ exp.temp, xlab = expression(paste("Exposure max. temp. (",degree,"C)")),  cex.axis = 1,
     ylab = expression(paste("Trans. rate - ",beta," (L spore"^-1,"day"^-1,"10"^-4,")")), ylim = c(0, 1),
     data = beta.u.out, type = "p", pch = 21, col = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4, cex.axis = 1.2)
axis(side = 1, at = c(20, 32), cex.axis = 1.2)
arrows(beta.u.out$exp.temp, beta.u.out$beta.lower*10000, beta.u.out$exp.temp, beta.u.out$beta.upper*10000, angle = 90, length = 0, code = 3, lwd = 1.4)
points(beta.calc*10000 ~ exp.temp, data = beta.u.out[1:2,], type = "o", pch = 21, bg = "white", cex = 1.75)
points(beta.calc*10000 ~ exp.temp, data = beta.u.out[3:4,], type = "o", pch = 21, bg = "grey40", lty = 3, cex = 1.75)

a <- expression(paste("20",degree,"C"))
b <- expression(paste("32",degree,"C"))
legend(24, 0.95, bty = "n", pch = 21, pt.bg = c("white", "grey40"), legend = c(a, b), pt.cex = 1.75, lty = c(1, 3), cex = 1.1)
text(x = 26, y = 0.95, labels = "Development max. temp.:", cex = 1.1)
legend(x = 17, y = 1.05, legend = "A", cex = 1.25, bty = "n", adj = 1)
text(x = 33.4, y = .61, labels = "a", cex = 1.1)
text(x = 18.6, y = .24, labels = "b", cex = 1.1)
text(x = 33.4, y = .46, labels = "a", cex = 1.1)
text(x = 18.6, y = .15, labels = "b", cex = 1.1)
text(x = 20.2, y = 0.02, labels = "+/- 95% CIs", cex = 1.1)

##### Panel B: f
plot(f.calc*100 ~ temp, xlab = expression(paste("Temperature (",degree,"C)")),  cex.axis = 1.2,
     ylab = expression(paste("Foraging rate - ",italic(f)," (L day"^-1,"10"^-2,")")), ylim = c(0.9, 2.6),
     data = f.out, type = "p", pch = 21, bg = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4, cex.axis = 1.2)
axis(side = 1, at = c(20, 25, 30), cex.axis = 1.2)

arrows(f.out$temp, f.out$f.lower*100, f.out$temp, f.out$f.upper*100, angle = 90, length = 0, code = 3, lwd = 1.4)
points(f.calc*100 ~ temp, data = f.out, type = "p", pch = 21, bg = c("white", "grey80", "grey40"), cex = 1.75)
legend(x = 17, y = 2.72, legend = "B", cex = 1.25, bty = "n", adj = 1)
text(x = 20.2, y = 0.95, labels = "+/- 95% CIs", cex = 1.1)
text(x = 28.6, y = 2.0, labels = "a", cex = 1.1)
text(x = 23.6, y = 2.12, labels = "a", cex = 1.1)
text(x = 18.6, y = 1.34, labels = "b", cex = 1.1)

##### Panel C: u
plot(u.calc*1000 ~ exp.temp, xlab = expression(paste("Exposure max. temp. (",degree,"C)")),  cex.axis = 1,
     ylab = expression(paste("Spore infectivity - ",italic(u)," (spore"^-1,"10"^-3,")")), ylim = c(0, 5),
     data = beta.u.out, type = "p", pch = 21, col = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4, cex.axis = 1.2)
axis(side = 1, at = c(20, 32), cex.axis = 1.2)
arrows(beta.u.out$exp.temp, beta.u.out$u.lower*1000, beta.u.out$exp.temp, beta.u.out$u.upper*1000, angle = 90, length = 0, code = 3, lwd = 1.4)
points(u.calc*1000 ~ exp.temp, data = beta.u.out[1:2,], type = "o", pch = 21, bg = "white", cex = 1.75)
points(u.calc*1000 ~ exp.temp, data = beta.u.out[3:4,], type = "o", pch = 21, bg = "grey40", lty = 3, cex = 1.75)

legend(24, 4.7, bty = "n", pch = 21, pt.bg = c("white", "grey40"), legend = c(a, b), pt.cex = 1.75, lty = c(2, 3), cex = 1.1)
text(x = 26, y = 4.75, labels = "Development max. temp.:", cex = 1.1)
legend(x = 17, y = 5.25, legend = "C", cex = 1.25, bty = "n", adj = 1)
text(x = 33.4, y = 3, labels = "a", cex = 1.1)
text(x = 18.6, y = 1.85, labels = "b,c", cex = 1.1)
text(x = 33.4, y = 2.2, labels = "a,b", cex = 1.1)
text(x = 18.6, y = 1.125, labels = "c", cex = 1.1)
text(x = 20.2, y = .1, labels = "+/- 95% CIs", cex = 1.1)