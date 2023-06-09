## Marta Shocket, Indiana University & UCLA
## Code for Can hot temperatures limit disease transmission? A test of mechanisms in a zooplankton-fungus system in Functional Ecology.
##
## Purpose: Analyze rearing effects and effects on free-living spores from high temperatures (26C and 32C).
##          Data come from common-garden experiments where hosts at 20C were exposed to different spore treatments.
##          Spores analyzed for rearing effect came from two experiments: (1) beta + u measurement assay
##          (20/20 and 32/32 treatments only) and (2) within-host parasite growth assay (20, 26, and 32 treatments).
##          Free-living spores were incubated in 20, 25 and 32C for 1 and 7 days. 
##
## Table of Contents:
##      1. Set working directory, load libraries and data
##      2. Rearing effect: calculate transmission rate
##      3. Rearing effect: bootstrap transmission rate
##      4. Rearing effect: transmission rate randomizations
##      5. Free-living spore effect: calculate transmission rate
##      6. Free-living spore effect: bootstrap transmission rate
##      7. Free-living spore effect: transmission rate randomizations
##      8. Transform rho into a relative measure
##      9. Transform phi into a relative, time-weighted measure
##      10. Transform phi into a relative, time-weighted measure - sensitivity analyses for Appendix
##      11. Figure 4


######
##### 1. Set working directory, load libraries and data
######

# Set wd
setwd("~/Dropbox/Research Hall Lab/Temperature/Upper Thermal Constraints/Final Code and Data")

# Load library
library(bbmle)

## Get data
re.data <- read.csv("UTC_SporeRearingEffectAssay.csv")
fe.data <- read.csv("UTC_FreeLivingSporeAssay.csv")

######
##### 2. Rearing effect: calculate transmission rate
######

## Subset data
beta.u.data <- subset(re.data, SporeSource == "beta.u.assay")
WHG.data <- subset(re.data, SporeSource == "WHG.assay")

WHG.data.20 <- subset(WHG.data, RearTemp == 20)
WHG.data.26 <- subset(WHG.data, RearTemp == 26)
WHG.data.32 <- subset(WHG.data, RearTemp == 32)

beta.u.data.20 <- subset(beta.u.data, RearTemp == 20)
beta.u.data.32 <- subset(beta.u.data, RearTemp == 32)

# Calculate beta values for each treatment - exposed to 30 sp / ml in 15 ml for 1 day
beta.bu.20 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = beta.u.data.20, method="Brent", lower = 0, upper = .1)

beta.bu.32 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = beta.u.data.32, method="Brent", lower = 0, upper = .1)

beta.WHG.20 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = WHG.data.20, method="Brent", lower = 0, upper = .1)

beta.WHG.26 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = WHG.data.26, method="Brent", lower = 0, upper = .1)

beta.WHG.32 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*1)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = WHG.data.32, method="Brent", lower = 0, upper = .1)

# Create output dataframe
re.out <- data.frame(Exp = c("beta.u", "beta.u", "WHG", "WHG", "WHG"), RearTemp = c(20, 32, 20, 26, 32), beta.calc = numeric(5), lower = numeric(5), upper = numeric(5))

# Store calculation - transform from beaker units to volume units (divide by host density = 1 / 15 ml)
re.out$beta.calc[1] <- coef(beta.bu.20)/(1/0.015)
re.out$beta.calc[2] <- coef(beta.bu.32)/(1/0.015)
re.out$beta.calc[3] <- coef(beta.WHG.20)/(1/0.015)
re.out$beta.calc[4] <- coef(beta.WHG.26)/(1/0.015)
re.out$beta.calc[5] <- coef(beta.WHG.32)/(1/0.015)

# Difference calculations
# 20 vs 26
(re.out$beta.calc[4] - re.out$beta.calc[1])/re.out$beta.calc[1]
(re.out$beta.calc[4] - re.out$beta.calc[3])/re.out$beta.calc[3]

# 20 vs 32
(re.out$beta.calc[1] - re.out$beta.calc[2])/re.out$beta.calc[1] 
(re.out$beta.calc[3] - re.out$beta.calc[5])/re.out$beta.calc[3] 

# 26 vs 32
(re.out$beta.calc[4] - re.out$beta.calc[2])/re.out$beta.calc[4] 
(re.out$beta.calc[4] - re.out$beta.calc[5])/re.out$beta.calc[4]


######
##### 3. Rearing effect: boostrap transmission rate
######

# Set up datasets for sampling
bsdata.bu.20.U <- rep(0, 29)
bsdata.bu.20.I <- rep(1, 16)
bsdata.bu.20 <- c(bsdata.bu.20.U, bsdata.bu.20.I)

bsdata.bu.32.U <- rep(0, 36)
bsdata.bu.32.I <- rep(1, 9)
bsdata.bu.32 <- c(bsdata.bu.32.U, bsdata.bu.32.I)

bsdata.WHG.20.U <- rep(0, 24)
bsdata.WHG.20.I <- rep(1, 19)
bsdata.WHG.20 <- c(bsdata.WHG.20.U, bsdata.WHG.20.I)

bsdata.WHG.26.U <- rep(0, 17)
bsdata.WHG.26.I <- rep(1, 28)
bsdata.WHG.26 <- c(bsdata.WHG.26.U, bsdata.WHG.26.I)

bsdata.WHG.32.U <- rep(0, 35)
bsdata.WHG.32.I <- rep(1, 10)
bsdata.WHG.32 <- c(bsdata.WHG.32.U, bsdata.WHG.32.I)

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

# Load constants
Vol <- 15 # ml
Spore.Dose <- 30 # sp/ml
Time <- 1 # day
Density <- 1/0.015 # hosts/L

# Bootstrap beta
bs.beta.bu.20 <- BootstrapBeta(bsdata.bu.20, 10000)
bs.beta.bu.32 <- BootstrapBeta(bsdata.bu.32, 10000)
bs.beta.WHG.20 <- BootstrapBeta(bsdata.WHG.20, 10000)
bs.beta.WHG.26 <- BootstrapBeta(bsdata.WHG.26, 10000)
bs.beta.WHG.32 <- BootstrapBeta(bsdata.WHG.32, 10000)

# Calculate quantiles and store output
re.out$upper[1] <- quantile(bs.beta.bu.20$beta, 0.025)
re.out$upper[2] <- quantile(bs.beta.bu.32$beta, 0.025)
re.out$upper[3] <- quantile(bs.beta.WHG.20$beta, 0.025)
re.out$upper[4] <- quantile(bs.beta.WHG.26$beta, 0.025)
re.out$upper[5] <- quantile(bs.beta.WHG.32$beta, 0.025)

re.out$lower[1] <- quantile(bs.beta.bu.20$beta, 0.975)
re.out$lower[2] <- quantile(bs.beta.bu.32$beta, 0.975)
re.out$lower[3] <- quantile(bs.beta.WHG.20$beta, 0.975)
re.out$lower[4] <- quantile(bs.beta.WHG.26$beta, 0.975)
re.out$lower[5] <- quantile(bs.beta.WHG.32$beta, 0.975)

#Save output
save(re.out, file = "re.out.Rsave")

# Gather and save bootstrapped values
bs.re <- data.frame(bu.20 = bs.beta.bu.20$beta, bu.32 = bs.beta.bu.32$beta, 
                    WHG.20 = bs.beta.WHG.20$beta, WHG.26 = bs.beta.WHG.26$beta, WHG.32 = bs.beta.WHG.32$beta)
save(bs.re, file = "bs.re.Rsave")


######
##### 4. Rearing effect: transmission rate randomizations
######

# Set up datasets for label shuffling
inf.2632 <- c(bsdata.WHG.26, bsdata.WHG.32)
treat.2632 <- c(rep(26, length(bsdata.WHG.26)), rep(32, length(bsdata.WHG.32)))
df.2632 <- data.frame(infect = inf2632, Treatment = treat2632)

inf.2026 <- c(bsdata.WHG.20, bsdata.WHG.26)
treat.2026 <- c(rep(20, length(bsdata.WHG.20)), rep(26, length(bsdata.WHG.26)))
df.2026 <- data.frame(infect = inf2026, Treatment = treat2026)

inf.2032 <- c(bsdata.WHG.20, bsdata.WHG.32)
treat.2032 <- c(rep(20, length(bsdata.WHG.20)), rep(32, length(bsdata.WHG.32)))
df.2032 <- data.frame(infect = inf2032, Treatment = treat2032)

inf.2032.1 <- c(bsdata.bu.20, bsdata.bu.32)
treat.2032.1 <- c(rep(20, length(bsdata.bu.20)), rep(32, length(bsdata.bu.32)))
df.2032.1 <- data.frame(infect = inf2032.1, Treatment = treat2032.1)

inf2026.1 <- c(bsdata.bu.20, bsdata.WHG.26)
treat2026.1 <- c(rep(20, length(bsdata.bu.20)), rep(26, length(bsdata.WHG.26)))
df2026.1 <- data.frame(infect = inf2026.1, Treatment = treat2026.1)

inf2632.1 <- c(bsdata.WHG.26, bsdata.bu.32)
treat2632.1 <- c(rep(26, length(bsdata.WHG.26)), rep(32, length(bsdata.bu.32)))
df2632.1 <- data.frame(infect = inf2632.1, Treatment = treat2632.1)

inf20s <- c(bsdata.bu.20, bsdata.WHG.20)
treat20s <- c(rep(1.20, length(bsdata.bu.20)), rep(3.20, length(bsdata.WHG.20)))
df20s <- data.frame(infect = inf20s, Treatment = treat20s)

inf32s <- c(bsdata.bu.32, bsdata.WHG.32)
treat32s <- c(rep(1.32, length(bsdata.bu.32)), rep(3.32, length(bsdata.WHG.32)))
df32s <- data.frame(infect = inf32s, Treatment = treat32s)

################## Function - BetaDiffRand
# Arguments: 
#   data.set          dataframe: 0/1 infection data
#   nBoots            integer: number of randomizations
# Returns: 
#   betadiff.out      dataframe:  bootstrapped fit coefficients

################## Function - NOTE: This version replaces 100% infected samples with highest value from MLE profiling
BetaDiffRand = function(data.set, treat.rand.1, treat.rand.2, nBoots) {
  
  # Create data frame for output
  betadiff.out = data.frame(beta.1 = numeric(nBoots), beta.2 = numeric(nBoots), beta.diff = numeric(nBoots))
  
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
    
    # Calculate beta parameters - divide by host density (1/L) to get correct units (L instead of tube)
    cf.beta.treat.1 <- coef(beta.treat.1)/Density
    cf.beta.treat.2 <- coef(beta.treat.2)/Density
    
    # Store betas/us and calculate the absolute value of their differences (so positive & negative extreme values both count)
    betadiff.out$beta.1[i] = cf.beta.treat.1
    betadiff.out$beta.2[i] = cf.beta.treat.2
    betadiff.out$beta.diff[i] = abs(betaudiff.out$beta.1[i] - betaudiff.out$beta.2[i])
    
  }
  
  betadiff.out # Return output
  
}


# Run the randomizations
betadiff.2632 <- BetaDiffRand(df2632, 26, 32, 10000)
betadiff.2026 <- BetaDiffRand(df2026, 20, 26, 10000)
betadiff.2032 <- BetaDiffRand(df2032, 20, 32, 10000)
betadiff.2032.1 <- BetaDiffRand(df2032.1, 20, 32, 10000)
betadiff.2026.1 <- BetaDiffRand(df2026.1, 20, 26, 10000)
betadiff.2632.1 <- BetaDiffRand(df2632.1, 26, 32, 10000)

betadiff.20s <- BetaDiffRand(df20s, 1.20, 3.20, 10000)
betadiff.32s <- BetaDiffRand(df32s, 1.32, 3.32, 10000)

hist(betadiff.2632$beta.diff)# <- abs(betadiff.2632$beta.diff)
hist(betadiff.2026$beta.diff)# <- abs(betadiff.2026$beta.diff)
hist(betadiff.2032$beta.diff)# <- abs(betadiff.2032$beta.diff)
hist(betadiff.2032.1$beta.diff)# <- abs(betadiff.2032.1$beta.diff)
hist(betadiff.2026.1$beta.diff)# <- abs(betadiff.2026.1$beta.diff)
hist(betadiff.2632.1$beta.diff)# <- abs(betadiff.2632.1$beta.diff)

hist(betadiff.20s$beta.diff)# <- abs(betadiff.20s$beta.diff)
hist(betadiff.32s$beta.diff)# <- abs(betadiff.32s$beta.diff)

coef(beta.WHG.26)*0.015 - coef(beta.WHG.32)*0.015
coef(beta.WHG.20)*0.015 - coef(beta.WHG.26)*0.015
coef(beta.WHG.20)*0.015 - coef(beta.WHG.32)*0.015
coef(beta.bu.20)*0.015 - coef(beta.bu.32)*0.015
coef(beta.bu.20)*0.015 - coef(beta.WHG.26)*0.015
coef(beta.WHG.26)*0.015 - coef(beta.bu.32)*0.015
coef(beta.bu.20)*0.015 - coef(beta.WHG.20)*0.015
coef(beta.bu.32)*0.015 - coef(beta.WHG.32)*0.015

1 - ecdf(betadiff.2632$beta.diff)(2.407113e-05)
1 - ecdf(betadiff.2026$beta.diff)(1.301009e-05)
1 - ecdf(betadiff.2032$beta.diff)(1.106104e-05)
1 - ecdf(betadiff.2032.1$beta.diff)(7.207452e-06)
1 - ecdf(betadiff.2026.1$beta.diff)(1.78027e-05)
1 - ecdf(betadiff.2632.1$beta.diff)(2.501015e-05)
1 - ecdf(betadiff.20s$beta.diff)(0.000004792609)
1 - ecdf(betadiff.32s$beta.diff)(0.0000009390217)

sporerear.rand.list <- list(betadiff.2026, betadiff.2632, betadiff.2032, betadiff.2026.1, betadiff.2632.1, betadiff.2032.1, betadiff.20s, betadiff.32s)
save(sporerear.rand.list, file = "sporerear.rand.list.RSave")
load("sporerear.rand.list.Rsave")


######
##### 5. Free-living effect: calculate transmission rate
######

## Subset data
fe.data.20.1 <- subset(fe.data, Duration == 1 & Temp == 20)
fe.data.25.1 <- subset(fe.data, Duration == 1 & Temp == 25)
fe.data.30.1 <- subset(fe.data, Duration == 1 & Temp == 30)
fe.data.20.7 <- subset(fe.data, Duration == 7 & Temp == 20)
fe.data.25.7 <- subset(fe.data, Duration == 7 & Temp == 25)
fe.data.30.7 <- subset(fe.data, Duration == 7 & Temp == 30)

# Calculate beta for each treatment - exposed to 30 sp / ml in 15 ml for 13 hours
beta.20.1 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*13/24)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = fe.data.20.1, method="Brent", lower = 0, upper = .1)

beta.20.7 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*13/24)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = fe.data.20.7, method="Brent", lower = 0, upper = .1)

beta.25.1 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*13/24)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = fe.data.25.1, method="Brent", lower = 0, upper = .1)

beta.25.7 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*13/24)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = fe.data.25.7, method="Brent", lower = 0, upper = .1)

beta.30.1 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*13/24)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = fe.data.30.1, method="Brent", lower = 0, upper = .1)

beta.30.7 = mle2(Uninf ~ dbinom(size=Total, prob=exp(-Beta*15*30*13/24)), start=list(Beta = 0.00004), 
                 control=list(parscale = c(Beta = 0.00001)), data = fe.data.30.7, method="Brent", lower = 0, upper = .1)

# Create output dataframe
fe.out <- data.frame(Duration = c(1, 1, 1, 7, 7, 7), IncTemp = c(20, 25, 30, 20, 25, 30), 
                          beta.calc = numeric(6), lower = numeric(6), upper = numeric(6))

# Store calculation - transform from beaker units to volume units (divide by host density = 1 / 15 ml)
fe.out$beta.calc[1] <- coef(beta.20.1)/(1/0.015)
fe.out$beta.calc[2] <- coef(beta.25.1)/(1/0.015)
fe.out$beta.calc[3] <- coef(beta.30.1)/(1/0.015)
fe.out$beta.calc[4] <- coef(beta.20.7)/(1/0.015)
fe.out$beta.calc[5] <- coef(beta.25.7)/(1/0.015)
fe.out$beta.calc[6] <- coef(beta.30.7)/(1/0.015)


######
##### 6. Free-living effect: bootstrap transmission rate
######

# Set up datasets for sampling
bsdata.20.1.U <- rep(0, 21)
bsdata.20.1.I <- rep(1, 19)
bsdata.20.1 <- c(bsdata.20.1.U, bsdata.20.1.I)

bsdata.25.1.U <- rep(0, 19)
bsdata.25.1.I <- rep(1, 19)
bsdata.25.1 <- c(bsdata.25.1.U, bsdata.25.1.I)

bsdata.30.1.U <- rep(0, 18)
bsdata.30.1.I <- rep(1, 20)
bsdata.30.1 <- c(bsdata.30.1.U, bsdata.30.1.I)

bsdata.20.7.U <- rep(0, 11)
bsdata.20.7.I <- rep(1, 28)
bsdata.20.7 <- c(bsdata.20.7.U, bsdata.20.7.I)

bsdata.25.7.U <- rep(0, 18)
bsdata.25.7.I <- rep(1, 20)
bsdata.25.7 <- c(bsdata.25.7.U, bsdata.25.7.I)

bsdata.30.7.U <- rep(0, 36)
bsdata.30.7.I <- rep(1, 4)
bsdata.30.7 <- c(bsdata.30.7.U, bsdata.30.7.I)

############# Load function BootstrapBeta() in section # 2 above

###### Experiment constants for bootstrapping function
Vol <- 15 # ml
Spore.Dose <- 30 # sp/ml
Time <- 13/24 # day
Density <- 1/0.015 # hosts/L

# Bootstrap beta
bs.beta.20.1 <- BootstrapBeta(bsdata.20.1, 10000)
bs.beta.25.1 <- BootstrapBeta(bsdata.25.1, 10000)
bs.beta.30.1 <- BootstrapBeta(bsdata.30.1, 10000)
bs.beta.20.7 <- BootstrapBeta(bsdata.20.7, 10000)
bs.beta.25.7 <- BootstrapBeta(bsdata.25.7, 10000)
bs.beta.30.7 <- BootstrapBeta(bsdata.30.7, 10000)

# Calculate quantiles and store output
fe.out$upper[1] <- quantile(bs.beta.20.1$beta, 0.025)
fe.out$upper[2] <- quantile(bs.beta.25.1$beta, 0.025)
fe.out$upper[3] <- quantile(bs.beta.30.1$beta, 0.025)
fe.out$upper[4] <- quantile(bs.beta.20.7$beta, 0.025)
fe.out$upper[5] <- quantile(bs.beta.25.7$beta, 0.025)
fe.out$upper[6] <- quantile(bs.beta.30.7$beta, 0.025)

fe.out$lower[1] <- quantile(bs.beta.20.1$beta, 0.975)
fe.out$lower[2] <- quantile(bs.beta.25.1$beta, 0.975)
fe.out$lower[3] <- quantile(bs.beta.30.1$beta, 0.975)
fe.out$lower[4] <- quantile(bs.beta.20.7$beta, 0.975)
fe.out$lower[5] <- quantile(bs.beta.25.7$beta, 0.975)
fe.out$lower[6] <- quantile(bs.beta.30.7$beta, 0.975)

#Save output
save(fe.out, file = "fe.out.Rsave")

# Gather and save bootstrapped values
bs.fe <- data.frame(fe.20.1 = bs.beta.20.1$beta, fe.25.1 = bs.beta.25.1$beta, fe.30.1 = bs.beta.30.1$beta, 
                    fe.20.7 = bs.beta.20.7$beta, fe.25.7 = bs.beta.25.7$beta, fe.30.7 = bs.beta.30.7$beta)
save(bs.fe, file = "bs.fe.Rsave")


######
##### 7. Free-living effect: transmission rate randomizations
######

# Set up datasets for label shuffling
inf2025.7 <- c(bsdata.20.7, bsdata.25.7)
treat2025.7 <- c(rep(20.7, length(bsdata.20.7)), rep(25.7, length(bsdata.25.7)))
df2025.7 <- data.frame(infect = inf2025.7, Treatment = treat2025.7)

inf2030.7 <- c(bsdata.20.7, bsdata.30.7)
treat2030.7 <- c(rep(20.7, length(bsdata.20.7)), rep(30.7, length(bsdata.30.7)))
df2030.7 <- data.frame(infect = inf2030.7, Treatment = treat2030.7)

inf2530.7 <- c(bsdata.25.7, bsdata.30.7)
treat2530.7 <- c(rep(25.7, length(bsdata.25.7)), rep(30.7, length(bsdata.30.7)))
df2530.7 <- data.frame(infect = inf2530.7, Treatment = treat2530.7)

inf2025.1 <- c(bsdata.20.1, bsdata.25.1)
treat2025.1 <- c(rep(20.1, length(bsdata.20.1)), rep(25.1, length(bsdata.25.1)))
df2025.1 <- data.frame(infect = inf2025.1, Treatment = treat2025.1)

inf2030.1 <- c(bsdata.20.1, bsdata.30.1)
treat2030.1 <- c(rep(20.1, length(bsdata.20.1)), rep(30.1, length(bsdata.30.1)))
df2030.1 <- data.frame(infect = inf2030.1, Treatment = treat2030.1)

inf2530.1 <- c(bsdata.25.1, bsdata.30.1)
treat2530.1 <- c(rep(25.1, length(bsdata.25.1)), rep(30.1, length(bsdata.30.1)))
df2530.1 <- data.frame(infect = inf2530.1, Treatment = treat2530.1)

inf20 <- c(bsdata.20.1, bsdata.20.7)
treat20 <- c(rep(20.1, length(bsdata.20.1)), rep(20.7, length(bsdata.20.7)))
df20 <- data.frame(infect = inf20, Treatment = treat20)

betadiff.2025.7 <- BetaDiffRand(df2025.7, 20.7, 25.7, 10000)
betadiff.2030.7 <- BetaDiffRand(df2030.7, 20.7, 30.7, 10000)
betadiff.2530.7 <- BetaDiffRand(df2530.7, 25.7, 30.7, 10000)

betadiff.2025.1 <- BetaDiffRand(df2025.1, 20.1, 25.1, 10000)
betadiff.2030.1 <- BetaDiffRand(df2030.1, 20.1, 30.1, 10000)
betadiff.2530.1 <- BetaDiffRand(df2530.1, 25.1, 30.1, 10000)

betadiff.20 <- BetaDiffRand(df20, 20.1, 20.7, 10000)

betadiff.2025.1$beta.diff <- abs(betadiff.2025.1$beta.diff)
betadiff.2025.7$beta.diff <- abs(betadiff.2025.7$beta.diff)
betadiff.2030.1$beta.diff <- abs(betadiff.2030.1$beta.diff)
betadiff.2030.7$beta.diff <- abs(betadiff.2030.7$beta.diff)
betadiff.2530.1$beta.diff <- abs(betadiff.2530.1$beta.diff)
betadiff.2530.7$beta.diff <- abs(betadiff.2530.7$beta.diff)

hist(betadiff.2025.1$beta.diff)
hist(betadiff.2025.7$beta.diff)
hist(betadiff.2030.1$beta.diff)
hist(betadiff.2030.7$beta.diff)
hist(betadiff.2530.1$beta.diff)
hist(betadiff.2530.7$beta.diff)

hist(betadiff.20$beta.diff)

spore.decay.rand.list <- list(betadiff.2025.7, betadiff.2030.7, betadiff.2530.7, betadiff.2025.1, betadiff.2030.1, betadiff.2530.1)
save(spore.decay.rand.list, file = "spore.decay.rand.list.Rsave")

# actual difference between estimated betas
coef(beta.20.7)*0.015 - coef(beta.25.7)*0.015
coef(beta.25.7)*0.015 - coef(beta.30.7)*0.015
coef(beta.20.7)*0.015 - coef(beta.30.7)*0.015

coef(beta.20.1)*0.015 - coef(beta.25.1)*0.015
coef(beta.25.1)*0.015 - coef(beta.30.1)*0.015
coef(beta.20.1)*0.015 - coef(beta.30.1)*0.015

coef(beta.20.1)*0.015 - coef(beta.20.7)*0.015

# Calculating p-values
1 - ecdf(betadiff.2025.7$beta.diff)(3.190479e-05)
1 - ecdf(betadiff.2530.7$beta.diff)(3.949869e-05)
1 - ecdf(betadiff.2030.7$beta.diff)(7.140347e-05)

1 - ecdf(betadiff.2025.1$beta.diff)(3.002467e-06)
1 - ecdf(betadiff.2530.1$beta.diff)(3.327206e-06)
1 - ecdf(betadiff.2030.1$beta.diff)(6.329672e-06)

1 - ecdf(betadiff.20$beta.diff)(3.823446e-05)


######
##### 8. Transform rho into a relative measure
######

# For the temperature treatments that have two spores sources (20C, 32C), make a combined dataset with 50% from each experiment
#    (for 26C, just use the data from the one experiment)
re.20 <- c(bs.re$bu.20[1:5000], bs.re$WHG.20[1:5000])
re.26 <- bs.re$WHG.26
re.32 <- c(bs.re$bu.32[1:5000], bs.re$WHG.32[1:5000])

# Shuffle the two datasets to mix them (and save a version of the 20 to use later)
re.20.shuffle <- sample(re.20, size = length(re.20), replace = FALSE)
re.20 <- sample(re.20, size = length(re.20), replace = FALSE)
re.32 <- sample(re.32, size = length(re.32), replace = FALSE)

# Calculate transmission rate relative to 20 degrees (using shuffled data for 20C treatment to add in measurement error)
bs.re$rel.20 <- re.20.shuffle / re.20
bs.re$rel.26 <- re.26 / re.20
bs.re$rel.32 <- re.32 / re.20

save(bs.re, file = "bs.re.Rsave")


######
##### 9. Transform phi into a relative, time-weighted measure
######

# Create data set to fit linear model coefficients for declining beta over time
#    (We assume that spores start at the 20C day 7 infectivity, since it was the highest value measured)
fe.tw <- data.frame(Temp = c(20, 20, 25, 25, 30, 30), Time = c(1, 7, 1, 7, 1, 7), 
                    Beta = c(fe.out$beta.calc[4], fe.out$beta.calc[4], 
                             fe.out$beta.calc[4], fe.out$beta.calc[5], 
                             fe.out$beta.calc[4], fe.out$beta.calc[6]))

# Assume the spores are removed over time (by host consumption) according to: y = exp(-.25*x)
#     So most spores are consumed quickly (70% within one-week), in line with Shocket et al. 2018 Ecology
curve(exp(-.25*x), from = 0, to = 7, ylim = c(0,1), ylab = "Proportion of spores remaining", xlab = "Time (days)")
legend("topleft", bty = "n", legend = "B")

# Meanwhile, transmission rate is declining over time from decreasing spore infectivity
par(mfrow = c(1,2), oma = c(0,0,0,0), mar = c(5, 5, 1.5, 1.5))
plot(Beta*10000 ~ Time, data = fe.tw[5:6,], type = "o", col = "red", ylim =c(0, 1), xlim = c(1, 7), lwd = 1.5, ylab = expression(paste("Transmission rate (L spore"^-1,"day"^-1,"10"^4,")")), xlab = "Time (days)")
points(Beta*10000 ~ Time, data = fe.tw[3:4,], type = "o", col = "purple", lwd = 1.5)
points(Beta*10000 ~ Time, data = fe.tw[1:2,], type = "o", col = "blue", lwd = 1.5)
a <- expression(paste("20",degree,"C"))
b <- expression(paste("25",degree,"C"))
c <- expression(paste("30",degree,"C"))
legend("bottomleft", legend = c(a, b, c), lwd = 1.5, col = c("blue", "purple", "red"), bty = "n")
legend("topleft", bty = "n", legend = "A")

# Fit linear models for slope coefficients
lm.20 <- lm(Beta ~ Time, data = fe.tw[1:2,])
lm.25 <- lm(Beta ~ Time, data = fe.tw[3:4,])
lm.30 <- lm(Beta ~ Time, data = fe.tw[5:6,])

# Create sequence of days / x-values
xs <- seq(1,7)

# Create sequence of spore proportions / y-values
ys <- c(1, exp(-.25*xs))
ys.prop <- numeric(7)
for(i in 1:6){ ys.prop[i] <- ys[i] - ys[i+1]} # calculate the difference in the proportion remaining on each day
ys.prop[7] <- ys[7] # for the last day, count all the remaining spores instead of the difference
sum(ys.prop) # check that it sums to 1 = total spores

# Calculate declining beta over time from the linear models
betas.30 <- coef(lm.30)[2]*xs + coef(lm.30)[1]
betas.25 <- coef(lm.25)[2]*xs + coef(lm.25)[1]
betas.20 <- coef(lm.20)[2]*xs + coef(lm.20)[1]

# Calculate the beta-weights based on the proportion of spores remaining each day
beta.30.weights <- betas.30*ys.prop
beta.25.weights <- betas.25*ys.prop
beta.20.weights <- betas.20*ys.prop

# Take the sum to calculate the time-weighted beta
w.beta.30 <- sum(beta.30.weights)
w.beta.25 <- sum(beta.25.weights)
w.beta.20 <- sum(beta.20.weights)

fe.out$beta.calc.tw <- c(NA, NA, NA, w.beta.20, w.beta.25, w.beta.30)
save(fe.out, file = "fe.out.Rsave")

# Create output for bootstrapped time-weighted free-living spore effect
bs.fe.tw <- data.frame(fe.tw.20 = numeric(10000), fe.tw.25 = numeric(10000), fe.tw.30 = numeric(10000))

###### Loop to bootstrap time-weighted free-living spore effect
for(i in 1:10000) {
  
  # Get beta values for each temp and put in data frame
  l.fe.tw <- data.frame(Temp = c(20, 20, 25, 25, 30, 30), Time = c(1, 7, 1, 7, 1, 7), 
                        Beta = c(bs.fe$fe.20.7[i], bs.fe$fe.20.7[i], 
                                 bs.fe$fe.20.7[i], bs.fe$fe.25.7[i], 
                                 bs.fe$fe.20.7[i], bs.fe$fe.30.7[i]))
  
  # Fit linear models for slope coefficients
  l.lm.20 <- lm(Beta ~ Time, data = l.fe.tw[1:2,])
  l.lm.25 <- lm(Beta ~ Time, data = l.fe.tw[3:4,])
  l.lm.30 <- lm(Beta ~ Time, data = l.fe.tw[5:6,])
  
  # Calculate declining beta over time from the linear models
  l.betas.30 <- coef(l.lm.30)[2]*xs + coef(l.lm.30)[1]
  l.betas.25 <- coef(l.lm.25)[2]*xs + coef(l.lm.25)[1]
  l.betas.20 <- coef(l.lm.20)[2]*xs + coef(l.lm.20)[1]
  
  # Calculate the beta-weights based on the proportion of spores remaining each day
  l.beta.30.weights <- l.betas.30*ys.prop
  l.beta.25.weights <- l.betas.25*ys.prop
  l.beta.20.weights <- l.betas.20*ys.prop
  
  # Take the sum to calculate the time-weighted beta
  bs.fe.tw$fe.tw.30[i] <- sum(l.beta.30.weights)
  bs.fe.tw$fe.tw.25[i] <- sum(l.beta.25.weights)
  bs.fe.tw$fe.tw.20[i] <- sum(l.beta.20.weights)
  
}

# Add the bootstrapped time-weighted measures to the dataframe
bs.fe$fe.20.tw <- bs.fe.tw$fe.tw.20
bs.fe$fe.25.tw <- bs.fe.tw$fe.tw.25
bs.fe$fe.30.tw <- bs.fe.tw$fe.tw.30

# Calculate and add the **relative** bootstrapped time-weighted measures to the dataframe
fe.tw.20.shuffle <- sample(bs.fe.tw$fe.tw.20, size = length(bs.fe.tw$fe.tw.20), replace = FALSE) 
bs.fe$fe.20.tw.rel <- fe.tw.20.shuffle / bs.fe.tw$fe.tw.20
bs.fe$fe.25.tw.rel <- bs.fe.tw$fe.tw.25 / bs.fe.tw$fe.tw.20
bs.fe$fe.30.tw.rel <- bs.fe.tw$fe.tw.30 / bs.fe.tw$fe.tw.20

save(bs.fe, file = "bs.fe.Rsave")
load("bs.fe.Rsave")
load("bs.re.Rsave")

# Calculate quantiles for relative parameters rho and phi
#   take the mean of bootstrapped samples (for rho, to combine both spore sources; for phi, to use time-weighted beta)
par.re.fe.out <- data.frame(Par = c("rho", "rho", "rho", "phi", "phi", "phi"), Temp = c(20, 26, 32, 20, 25, 30), 
                       mean = c(mean(bs.re$rel.20), mean(bs.re$rel.26), mean(bs.re$rel.32),
                                    mean(bs.fe$fe.20.tw.rel), mean(bs.fe$fe.25.tw.rel), mean(bs.fe$fe.30.tw.rel)),
                       lower = c(quantile(bs.re$rel.20, 0.025), quantile(bs.re$rel.26, 0.025), quantile(bs.re$rel.32, 0.025),
                                 quantile(bs.fe$fe.20.tw.rel, 0.025), quantile(bs.fe$fe.25.tw.rel, 0.025), quantile(bs.fe$fe.30.tw.rel, 0.025)), 
                       upper = c(quantile(bs.re$rel.20, 0.975), quantile(bs.re$rel.26, 0.975), quantile(bs.re$rel.32, 0.975),
                                 quantile(bs.fe$fe.20.tw.rel, 0.975), quantile(bs.fe$fe.25.tw.rel, 0.975), quantile(bs.fe$fe.30.tw.rel, 0.975)))
save("par.re.fe.out.Rsave")


######
##### 10. Transform phi into a relative, time-weighted measure - sensitivity analyses for Appendix
######

load("fe.out.Rsave")
load("bs.fe.Rsave")

# Function to create sequence of spore proportions / y-values for a given exponential decay model
spore.ys = function(x.seq, cf.exp){
  
  cf.exp <- cf.exp
  ys <- c(1, exp(cf.exp*xs))
  ys.prop <- numeric(7)
  for(i in 1:6){ ys.prop[i] <- ys[i] - ys[i+1]} # calculate the difference in the proportion remaining on each day
  ys.prop[7] <- ys[7] # for the last day, count all the remaining spores instead of the difference
  
  ys.prop # return
}

# Get spore proportions / y-values for each exponential decay model
xs <- seq(1,7)
spore.ys.exp10 <- spore.ys(xs, -0.1)
spore.ys.exp15 <- spore.ys(xs, -0.15)
spore.ys.exp25 <- spore.ys(xs, -0.25)
spore.ys.exp50 <- spore.ys(xs, -0.5)

# Function to calculate time-weighted beta for each exponential decay model
fe.bs.slopes = function(ys.prop){
  
  out <- data.frame(fe.tw.20 = numeric(10000), fe.tw.25 = numeric(10000), fe.tw.30 = numeric(10000))
  
  ###### Loop to bootstrap time-weighted free-living spore effect
  for(i in 1:10000) {
    
    # Get beta values for each temp and put in data frame
    l.fe.tw <- data.frame(Temp = c(20, 20, 25, 25, 30, 30), Time = c(1, 7, 1, 7, 1, 7), 
                          Beta = c(bs.fe$fe.20.7[i], bs.fe$fe.20.7[i], 
                                   bs.fe$fe.20.7[i], bs.fe$fe.25.7[i], 
                                   bs.fe$fe.20.7[i], bs.fe$fe.30.7[i]))
    
    # Fit linear models for slope coefficients
    l.lm.20 <- lm(Beta ~ Time, data = l.fe.tw[1:2,])
    l.lm.25 <- lm(Beta ~ Time, data = l.fe.tw[3:4,])
    l.lm.30 <- lm(Beta ~ Time, data = l.fe.tw[5:6,])
    
    # Calculate declining beta over time from the linear models
    l.betas.30 <- coef(l.lm.30)[2]*xs + coef(l.lm.30)[1]
    l.betas.25 <- coef(l.lm.25)[2]*xs + coef(l.lm.25)[1]
    l.betas.20 <- coef(l.lm.20)[2]*xs + coef(l.lm.20)[1]
    
    # Calculate the beta-weights based on the proportion of spores remaining each day
    l.beta.30.weights <- l.betas.30*ys.prop
    l.beta.25.weights <- l.betas.25*ys.prop
    l.beta.20.weights <- l.betas.20*ys.prop
    
    # Take the sum to calculate the time-weighted beta
    out$fe.tw.30[i] <- sum(l.beta.30.weights)
    out$fe.tw.25[i] <- sum(l.beta.25.weights)
    out$fe.tw.20[i] <- sum(l.beta.20.weights)
    
  }
  
  out #return
}

# Calculate time weighted beta for each exponential decay model
fe.tw.betas.exp10 <- fe.bs.slopes(spore.ys.exp10)
fe.tw.betas.exp15 <- fe.bs.slopes(spore.ys.exp15)
fe.tw.betas.exp25 <- fe.bs.slopes(spore.ys.exp25)
fe.tw.betas.exp50 <- fe.bs.slopes(spore.ys.exp50)
fe.tw.betas.max <- data.frame(fe.tw.20 = bs.fe$fe.20.7, fe.tw.25 = bs.fe$fe.25.7, fe.tw.30 = bs.fe$fe.30.7) # this is just the 7-day beta value

# Calculate **relative** time weighted beta for each exponential decay model
# exp = 0.1 model
fe.tw.20.exp10.shuffle <- sample(fe.tw.betas.exp10$fe.tw.20, size = length(fe.tw.betas.exp10$fe.tw.20), replace = FALSE) 
fe.tw.betas.exp10$fe.tw.20.rel <- fe.tw.20.exp10.shuffle / fe.tw.betas.exp10$fe.tw.20
fe.tw.betas.exp10$fe.tw.25.rel <- fe.tw.betas.exp10$fe.tw.25 / fe.tw.betas.exp10$fe.tw.20
fe.tw.betas.exp10$fe.tw.30.rel <- fe.tw.betas.exp10$fe.tw.30 / fe.tw.betas.exp10$fe.tw.20

# exp = 0.15 model
fe.tw.20.exp15.shuffle <- sample(fe.tw.betas.exp15$fe.tw.20, size = length(fe.tw.betas.exp15$fe.tw.20), replace = FALSE) 
fe.tw.betas.exp15$fe.tw.20.rel <- fe.tw.20.exp15.shuffle / fe.tw.betas.exp15$fe.tw.20
fe.tw.betas.exp15$fe.tw.25.rel <- fe.tw.betas.exp15$fe.tw.25 / fe.tw.betas.exp15$fe.tw.20
fe.tw.betas.exp15$fe.tw.30.rel <- fe.tw.betas.exp15$fe.tw.30 / fe.tw.betas.exp15$fe.tw.20

# exp = 0.25 model
fe.tw.20.exp25.shuffle <- sample(fe.tw.betas.exp25$fe.tw.20, size = length(fe.tw.betas.exp25$fe.tw.20), replace = FALSE) 
fe.tw.betas.exp25$fe.tw.20.rel <- fe.tw.20.exp25.shuffle / fe.tw.betas.exp25$fe.tw.20
fe.tw.betas.exp25$fe.tw.25.rel <- fe.tw.betas.exp25$fe.tw.25 / fe.tw.betas.exp25$fe.tw.20
fe.tw.betas.exp25$fe.tw.30.rel <- fe.tw.betas.exp25$fe.tw.30 / fe.tw.betas.exp25$fe.tw.20

# exp = 0.5 model
fe.tw.20.exp50.shuffle <- sample(fe.tw.betas.exp50$fe.tw.20, size = length(fe.tw.betas.exp50$fe.tw.20), replace = FALSE) 
fe.tw.betas.exp50$fe.tw.20.rel <- fe.tw.20.exp50.shuffle / fe.tw.betas.exp50$fe.tw.20
fe.tw.betas.exp50$fe.tw.25.rel <- fe.tw.betas.exp50$fe.tw.25 / fe.tw.betas.exp50$fe.tw.20
fe.tw.betas.exp50$fe.tw.30.rel <- fe.tw.betas.exp50$fe.tw.30 / fe.tw.betas.exp50$fe.tw.20

# max model (7 day assay values)
fe.tw.20.max.shuffle <- sample(fe.tw.betas.max$fe.tw.20, size = length(fe.tw.betas.max$fe.tw.20), replace = FALSE) 
fe.tw.betas.max$fe.tw.20.rel <- fe.tw.20.max.shuffle / fe.tw.betas.max$fe.tw.20
fe.tw.betas.max$fe.tw.25.rel <- fe.tw.betas.max$fe.tw.25 / fe.tw.betas.max$fe.tw.20
fe.tw.betas.max$fe.tw.30.rel <- fe.tw.betas.max$fe.tw.30 / fe.tw.betas.max$fe.tw.20

# Function to calculate quantiles and save output as a list
calcQuantFeTwBetas = function(input){
  
  output = data.frame(temp = c(20, 25, 30), mean = numeric(3), upperCI = numeric(3), lowerCI = numeric(3))
  
  output$mean[1] <- mean(input$fe.tw.20.rel)
  output$mean[2] <- mean(input$fe.tw.25.rel)
  output$mean[3] <- mean(input$fe.tw.30.rel)
  
  output$upperCI[1] <- quantile(input$fe.tw.20.rel, probs = 0.975)
  output$upperCI[2] <- quantile(input$fe.tw.25.rel, probs = 0.975)
  output$upperCI[3] <- quantile(input$fe.tw.30.rel, probs = 0.975)
  
  output$lowerCI[1] <- quantile(input$fe.tw.20.rel, probs = 0.025)
  output$lowerCI[2] <- quantile(input$fe.tw.25.rel, probs = 0.025)
  output$lowerCI[3] <- quantile(input$fe.tw.30.rel, probs = 0.025)
  
  output #return
}

# Calculate quantiles
fe.tw.betas.exp10.out <- calcQuantFeTwBetas(fe.tw.betas.exp10)
fe.tw.betas.exp15.out <- calcQuantFeTwBetas(fe.tw.betas.exp15)
fe.tw.betas.exp25.out <- calcQuantFeTwBetas(fe.tw.betas.exp25)
fe.tw.betas.exp50.out <- calcQuantFeTwBetas(fe.tw.betas.exp50)
fe.tw.betas.max.out <- calcQuantFeTwBetas(fe.tw.betas.max)

fe.tw.betas.sense.out.list <- list(fe.tw.betas.exp10.out, fe.tw.betas.exp15.out, fe.tw.betas.exp25.out, 
                                   fe.tw.betas.exp50.out, fe.tw.betas.max.out)

fe.tw.betas.sense.list <- list(fe.tw.betas.exp10, fe.tw.betas.exp15, fe.tw.betas.exp25, 
                                   fe.tw.betas.exp50, fe.tw.betas.max)

save(fe.tw.betas.sense.list, file = "fe.tw.betas.sense.list.Rsave")
save(fe.tw.betas.sense.out.list, file = "fe.tw.betas.sense.out.list.Rsave")

### Transmission potential calculations for sensitivity analysis in UTC_TransmissionPotential.R


######
##### 11. Figure 4
######

load("re.out.Rsave")
re.out.bu <- subset(re.out, Exp == "beta.u")
re.out.WHG <- subset(re.out, Exp == "WHG")
re.out.bu$RearTemp <- re.out.bu$Temp - 0.4
re.out.WHG$RearTemp <- re.out.WHG$Temp + 0.4

load("fe.out.Rsave")
fe.out.1 <- subset(fe.out, Duration == 1)
fe.out.7 <- subset(fe.out, Duration == 7)
fe.out.7$IncTemp <- fe.out.7$IncTemp + .4
fe.out.1$IncTemp <- fe.out.1$IncTemp - .4

load("par.re.fe.out.Rsave")
re.out.rel <- subset(par.re.fe.out, Par == "rho")
fe.out.rel <- subset(par.re.fe.out, Par == "phi")

# plot settings
par(mfrow = c(2,2), las = 1, oma = c(0, 0, 3, 0), mar = c(4.5, 4.5, 1, 1), mgp = c(3, 1, 0))

##### Panel A: rearing effect transmission rate
plot(beta.calc*10000 ~ RearTemp, xlab = expression(paste("Maximum rearing temperature (",degree,"C)")), 
     ylab = expression(paste("Transmission rate - ",beta," (L spore"^-1,"day"^-1,"10"^-4,")")), ylim = c(0, .5),
     data = re.out.WHG, type = "n", pch = 21, bg = "white", xaxt = "n", xlim = c(17.5,34.5), cex.lab = 1.15, cex.axis = 0.9)
axis(side = 1, at = c(20, 26, 32))
mtext(side = 3, text = expression(paste("Rearing effect - ",rho)), line = 0.25, cex = 1.25)
arrows(re.out.bu$RearTemp, re.out.bu$lower*10000, re.out.bu$RearTemp, re.out.bu$upper*10000, angle = 90, length = 0, code = 3, lwd = 1.4)
arrows(re.out.WHG$RearTemp, re.out.WHG$lower*10000, re.out.WHG$RearTemp, re.out.WHG$upper*10000, angle = 90, length = 0, code = 3, lwd = 1.4)
points(beta.calc*10000 ~ RearTemp, data = re.out.bu, type = "p", pch = 22, col = "black", bg = c("white", "grey40"), cex = 1.5)
points(beta.calc*10000 ~ RearTemp, data = re.out.WHG, type = "p", pch = 23, col = "black", bg = c("white", "grey80", "grey40"), lty = 2, cex = 1.5)

a <- expression(paste(italic(beta + u)," assay")) 
legend(27.7, 0.45, bty = "n", pch = c(22, 23), pt.bg = c("white"), legend = c(a, "WHPG assay"), pt.cex = 1.5, cex = 1)
text(x = 30.5, y = 0.455, labels = "Spore source:", cex = 1)
text(x = 20.8, y = 0.01, labels = "+/- 95% CIs", cex = 1)
text(x = 18.4, y = 0.195, labels = "a,b", cex = 1)
text(x = 18, y = 0.15, labels = "b,c", cex = 1)
text(x = 25.2, y = 0.325, labels = "a", cex = 1)
text(x = 34, y = 0.085, labels = "c", cex = 1)
text(x = 33.4, y = 0.07, labels = "c", cex = 1)
text(x = 18.4, y = 0.49, labels = "A", cex = 1.1)

##### Panel B: free living effect transmission rate
plot(beta.calc*10000 ~ IncTemp, xlab = expression(paste("Spore incubation temperature (",degree,"C)")), 
     ylab = expression(paste("Transmission rate - ",beta," (L spore"^-1,"day"^-1,"10"^-4,")")), ylim = c(0, 1.2),
     data = fe.out.1, type = "n", pch = 21, bg = "white", xaxt = "n", xlim = c(18.8,31.2), cex.lab = 1.15, cex.axis = 0.9)
axis(side = 1, at = c(20, 25, 30))
mtext(side = 3, text = expression(paste("Free-living spore effect - ",phi)), line = 0.25, cex = 1.25)
arrows(fe.out.1$IncTemp, fe.out.1$lower*10000, fe.out.1$IncTemp, fe.out.1$upper*10000, angle = 90, length = 0, code = 3, lwd = 1.4)
arrows(fe.out.7$IncTemp, fe.out.7$lower*10000, fe.out.7$IncTemp, fe.out.7$upper*10000, angle = 90, length = 0, code = 3, lwd = 1.4)
points(beta.calc*10000 ~ IncTemp, data = fe.out.1, type = "o", pch = 22, bg = c("white", "grey80", "grey40"), lty = 2, cex = 1.5)
points(beta.calc*10000 ~ IncTemp, data = fe.out.7, type = "o", pch = 23, bg = c("white", "grey80", "grey40"), lty = 1, cex = 1.5)

legend(26.2, 1.07, bty = "n", pch = c(22, 23), pt.bg = c("white"), legend = c("1 day", "7 days"), pt.cex = 1.5, lty = c(2, 1), cex = 1)
text(x = 27.5, y = 1.08, labels = "Spore incubation length:", cex = 1)
text(x = 21.35, y = 0.02, labels = "+/- 95% CIs", cex = 1)
text(x = 19.6, y = .78, labels = "a", cex = 1)
text(x = 18.9, y = .4, labels = "b", cex = 1)
text(x = 24, y = .36, labels = "b", cex = 1)
text(x = 26.1, y = .52, labels = "b", cex = 1)
text(x = 30.4, y = .47, labels = "b", cex = 1)
text(x = 31.2, y = .07, labels = "c", cex = 1)
text(x = 19.2, y = 1.18, labels = "B", cex = 1.1)

##### Panel C: rearing effect parameter
plot(mean ~ Temp, xlab = expression(paste("Maximum rearing temperature (",degree,"C)")), 
     ylab = expression(paste("Parameter value - ",rho)), ylim = c(0, 4.2),
     data = re.out.rel, type = "n", pch = 21, bg = "white", xaxt = "n", xlim = c(17.5,34.5), cex.lab = 1.15, cex.axis = 0.9)
axis(side = 1, at = c(20, 26, 32))
arrows(re.out.rel$Temp, re.out.rel$lower, re.out.rel$Temp, re.out.rel$upper, angle = 90, length = 0, code = 3, lwd = 1.4)
points(mean ~ Temp, data = re.out.rel, type = "p", pch = 21, col = "black", bg = c("white", "grey80", "grey40"), cex = 1.5)
text(x = 20.8, y = 0.08, labels = "+/- 95% CIs", cex = 1)
text(x = 18.4, y = 4.1, labels = "C", cex = 1.1)

##### Panel D: free living effect parameter
plot(mean ~ Temp, xlab = expression(paste("Spore incubation temperature (",degree,"C)")), 
     ylab = expression(paste("Parameter value - ",phi)), ylim = c(0.25, 2),
     data = fe.out.rel, type = "n", pch = 21, bg = "white", xaxt = "n", xlim = c(18.8,31.2), cex.lab = 1.15, cex.axis = 0.9)
axis(side = 1, at = c(20, 25, 30))
arrows(fe.out.rel$Temp, fe.out.rel$lower, fe.out.rel$Temp, fe.out.rel$upper, angle = 90, length = 0, code = 3, lwd = 1.4)
points(mean ~ Temp, data = fe.out.rel, type = "p", pch = 21, bg = c("white", "grey80", "grey40"), lty = 2, cex = 1.5)
text(x = 21.35, y = 0.3, labels = "+/- 95% CIs", cex = 1)
text(x = 19.35, y = 1.95, labels = "D", cex = 1.1)
