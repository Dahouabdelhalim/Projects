##### Marta Shocket, Stanford University
##### Updated April 2018

##### Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for mosquito and pathogen traits 
##             to be used as priors for RRV model
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit trait thermal reponses to data
##           5) Fit gamma distributions & save prior hyperparameters


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

# Load trait data for priors
data.prior <- read.csv("RRVPriorData.csv")

# Subset data
data.a.prior <- subset(data.prior, trait.name == "a")
data.EFD.prior <- subset(data.prior, trait.name == "EFD")
data.PDR.prior <- subset(data.prior, trait.name == "PDR")
data.pEA.prior <- subset(data.prior, trait.name == "pEA")
data.MDR.prior <- subset(data.prior, trait.name == "MDR")

# Lifespan data needs to be formatted
data.ls1.prior <- subset(data.prior, trait.name == "1/mu")
data.ls2.prior <- subset(data.prior, trait.name == "mu")
data.ls2.prior$trait <- 1/data.ls2.prior$trait # Convert Yang et al. data to lifespan not mu
data.ls.prior <- rbind(data.ls1.prior, data.ls2.prior) # combine the two data lifespan sets


##########
###### 2. JAGS Models
##########

############## Quadratic Model with uniform priors

sink("quad.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()

############## Briere Model with uniform priors

sink("briere.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 20)
    cf.Tm ~ dunif(28, 45)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### 3. Shared settings for all models
##########

##### inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Temp sequence for derived quantity calculations
Temp.xs <- seq(5, 45, 0.2)
N.Temp.xs <-length(Temp.xs)


##########
###### 4. Fits for priors
##########

############## 1/mu = lifespan for Ae. aegypti (to widen artifically narrow fit from overfitting) - quadratic

##### Get data
data <- data.ls.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
ls.prior.out<-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
ls.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(ls.prior.out)

save(ls.prior.out, file = "jagsout_ls_prior.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(0, 45), ylim = c(0,42), data = data.ls.prior, ylab = "Lifespan Prior", xlab = "Temperature")
lines(ls.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(ls.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(ls.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Calculate optimum for lifespan: 23.0 C
Temp.xs[which.max(as.vector(ls.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## EFD - Briere

##### Get data
data <- data.EFD.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EFD.prior.out <-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EFD.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(EFD.prior.out)

save(EFD.prior.out, file = "jagsout_EFD_prior.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,25), data = data.EFD.prior, ylab = "Fecundity Prior", xlab = "Temperature")
lines(EFD.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Calculate optimum for EFD: 28.8 C
Temp.xs[which.max(as.vector(EFD.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## a = biting rate for priors - Briere

##### Get data
data <- data.a.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
a.prior.out <-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(a.prior.out)

save(a.prior.out, file = "jagsout_a_prior.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.6), data = data.a.prior, ylab = "Biting rate Prior", xlab = "Temperature")
lines(a.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Calculate optimum for a: 35.0 C
Temp.xs[which.max(as.vector(a.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]

############## MDR  for 3 spp. - Briere

##### Get data
data <- data.MDR.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.prior.out <-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.prior.out)

save(MDR.prior.out, file = "jagsout_MDR_prior.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.3), data = data.MDR.prior, ylab = "MDR Prior", xlab = "Temperature")
lines(MDR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Calculate optimum for MDR: 35.0 C
Temp.xs[which.max(as.vector(MDR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## pEA - Briere

##### Get data
data <- data.pEA.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pEA.prior.out <-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pEA.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pEA.prior.out)

save(pEA.prior.out, file = "jagsout_pEA_prior.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pEA.prior, ylab = "pEA Prior", xlab = "Temperature")
lines(pEA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pEA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pEA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Calculate optimum for pEA: 35.0 C
Temp.xs[which.max(as.vector(pEA.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## PDR = 1/EIP for Ae. aegypti and albopictus - Briere

##### Get data
data <- data.PDR.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.prior.out)

save(PDR.prior.out, file = "jagsout_PDR_prior.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.prior, ylab = "PDR for Ae. aegypti + albopictus", xlab = "Temperature")
lines(PDR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Calculate optimum for PDR: 35.0 C
Temp.xs[which.max(as.vector(PDR.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


##########
###### 5. Fit gamma distributions & save prior hyperparameters
##########

# Load Prior Fits
load("jagsout_ls_prior.Rdata")
load("jagsout_EFD_prior.Rdata")
load("jagsout_a_prior.Rdata")
load("jagsout_MDR_prior.Rdata")
load("jagsout_pEA_prior.Rdata")
load("jagsout_PDR_prior.Rdata")

######## lifespan

# Get the posterior dists for all 4 parameters into a data frame
ls.prior.cf.dists <- data.frame(q = as.vector(ls.prior.out$BUGSoutput$sims.list$cf.q),
                               T0 = as.vector(ls.prior.out$BUGSoutput$sims.list$cf.T0), 
                               Tm = as.vector(ls.prior.out$BUGSoutput$sims.list$cf.Tm), 
                               sigma = as.vector(ls.prior.out$BUGSoutput$sims.list$cf.sigma))

# Fit gamma distributions for each parameter posterior dists
ls.prior.gamma.fits = apply(ls.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

######## EFD

# Get the posterior dists for all 4 parameters into a data frame
EFD.prior.dists <- data.frame(q = as.vector(EFD.prior.out$BUGSoutput$sims.list$cf.q),
                               T0 = as.vector(EFD.prior.out$BUGSoutput$sims.list$cf.T0), 
                               Tm = as.vector(EFD.prior.out$BUGSoutput$sims.list$cf.Tm), 
                               sigma = as.vector(EFD.prior.out$BUGSoutput$sims.list$cf.sigma))

# Fit gamma distributions for each parameter posterior dists
EFD.prior.gamma.fits = apply(EFD.prior.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

######## a

# Get the posterior dists for all 4 parameters into a data frame
a.prior.cf.dists <- data.frame(q = as.vector(a.prior.out$BUGSoutput$sims.list$cf.q),
                               T0 = as.vector(a.prior.out$BUGSoutput$sims.list$cf.T0),
                               Tm = as.vector(a.prior.out$BUGSoutput$sims.list$cf.Tm),
                               sigma = as.vector(a.prior.out$BUGSoutput$sims.list$cf.sigma))

# Fit gamma distributions for each parameter posterior dists - '2' indicates apply to columns not rows in df
a.prior.gamma.fits = apply(a.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

######## pEA

# Get the posterior dists for all 4 parameters into a data frame
pEA.prior.cf.dists <- data.frame(q = as.vector(pEA.prior.out$BUGSoutput$sims.list$cf.q),
                                 T0 = as.vector(pEA.prior.out$BUGSoutput$sims.list$cf.T0),
                                 Tm = as.vector(pEA.prior.out$BUGSoutput$sims.list$cf.Tm),
                                 sigma = as.vector(pEA.prior.out$BUGSoutput$sims.list$cf.sigma))

# Fit gamma distributions for each parameter posterior dists - '2' indicates apply to columns not rows in df
pEA.prior.gamma.fits = apply(pEA.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

######## MDR

# Get the posterior dists for all 4 parameters into a data frame
MDR.prior.cf.dists <- data.frame(q = as.vector(MDR.prior.out$BUGSoutput$sims.list$cf.q),
                                 T0 = as.vector(MDR.prior.out$BUGSoutput$sims.list$cf.T0),
                                 Tm = as.vector(MDR.prior.out$BUGSoutput$sims.list$cf.Tm),
                                 sigma = as.vector(MDR.prior.out$BUGSoutput$sims.list$cf.sigma))

# Fit gamma distributions for each parameter posterior dists - '2' indicates apply to columns not rows in df
MDR.prior.gamma.fits = apply(MDR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

######## PDR

# Get the posterior dists for all 4 parameters into a data frame
PDR.prior.cf.dists <- data.frame(q = as.vector(PDR.prior.out$BUGSoutput$sims.list$cf.q),
                               T0 = as.vector(PDR.prior.out$BUGSoutput$sims.list$cf.T0),
                               Tm = as.vector(PDR.prior.out$BUGSoutput$sims.list$cf.Tm),
                               sigma = as.vector(PDR.prior.out$BUGSoutput$sims.list$cf.sigma))

# Fit gamma distributions for each parameter posterior dists - '2' indicates apply to columns not rows in df
PDR.prior.gamma.fits = apply(PDR.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

prior.fits.list <- list(ls.prior.gamma.fits, EFD.prior.gamma.fits, a.prior.gamma.fits,
                        pEA.prior.gamma.fits, MDR.prior.gamma.fits, PDR.prior.gamma.fits)

save(prior.fits.list, file = "PriorFitsList.Rdata")