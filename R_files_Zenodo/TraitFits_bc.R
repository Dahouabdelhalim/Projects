## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for vector competence tarits (bc, b, c) for six arboviruses in many Culex and Aedes species
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit bc thermal responses with uniform priors
##           5) Fit bc thermal responses for priors
##           6) Fit gamma distributions to bc prior thermal responses
##           7) Fit bc thermal responses with data-informed priors
##           8) Calculate treatment averages for plotting


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("~/Fitting Traits")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')

# Load Data
data.bc <- read.csv("TraitData_bc.csv") # Data from database for most traits (except below)
unique(data.bc$joint.code)

# Subset Data
data.c.CpipWNV <- subset(data.bc, joint.code == "CpipWNV" & trait.name == "c")
data.bc.CpipWNV <- subset(data.bc, joint.code == "CpipWNV" & trait.name == "bc")

data.b.CtarWNV <- subset(data.bc, joint.code == "CtarWNV" & trait.name == "b")
data.bc.CuniWNV <- subset(data.bc, joint.code == "CuniWNV")

data.bc.CtarWEEV <- subset(data.bc, joint.code == "CtarWEEV" & trait.name == "bc")
data.c.CtarWEEV <- subset(data.bc, joint.code == "CtarWEEV" & trait.name == "c")
data.b.CtarWEEV <- subset(data.bc, joint.code == "CtarWEEV" & trait.name == "b")
data.c.CtarSLEV <- subset(data.bc, joint.code == "CtarSLEV" & trait.name == "c")
data.b.CtarSLEV <- subset(data.bc, joint.code == "CtarSLEV" & trait.name == "b")

data.c.AtaeSINV <- subset(data.bc, joint.code == "AtaeSINV")
data.c.CpipSINV <- subset(data.bc, joint.code == "CpipSINV")

data.bc.AtaeRVFV <- subset(data.bc, joint.code == "AtaeRVFV")
data.bc.AtriEEEV <- subset(data.bc, joint.code == "AtriEEEV")

data.bc.only <- subset(data.bc, trait.name == "bc")
data.b.only <- subset(data.bc, trait.name == "b")
data.c.only <- subset(data.bc, trait.name == "c")

par(mfrow = c(1,1))
plot(trait ~ T, xlim = c(5, 45), data = data.bc.only, ylab = "bc", xlab = "Temperature")
points(trait ~ T, data = data.bc.CpipWNV, col = "grey")
points(trait ~ T, data = data.bc.CuniWNV, col = "orange")
points(trait ~ T, data = data.bc.CtarWEEV, col = "blue")
points(trait ~ T, data = data.bc.AtriEEEV, col = "violet")
points(trait ~ T, data = data.bc.AtaeRVFV, col = "green")

plot(trait ~ T, xlim = c(5, 45), data = data.c.only, ylab = "b", xlab = "Temperature")
points(trait ~ T, data = data.c.CpipWNV, col = "grey")
points(trait ~ T, data = data.c.CtarWEEV, col = "dodgerblue")
points(trait ~ T, data = data.c.CtarSLEV, col = "navyblue")
points(trait ~ T, data = data.c.CpipSINV, col = "grey30")
points(trait ~ T, data = data.c.AtaeSINV, col = "darkgreen")

plot(trait ~ T, xlim = c(5, 45), data = data.b.only, ylab = "b", xlab = "Temperature")
points(trait ~ T, data = data.b.CtarWNV, col = "blue")
points(trait ~ T, data = data.b.CtarWEEV, col = "dodgerblue")
points(trait ~ T, data = data.b.CtarSLEV, col = "navyblue")


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
    cf.Tm ~ dunif(26, 50)
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


############## Quadratic Model with uniform priors - derived quantities always =< 1 (i.e., for probabilities)

sink("quadprob.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(26, 50)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
    }
    
    } # close model
    ",fill=T)
sink()


############## Quadratic Model with gamma priors (except sigma) - truncated for probabilities

sink("quadprob_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.q * (temp[i] - cf.T0) * (temp[i] - cf.Tm) * (cf.Tm > temp[i]) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i]) * (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) < 1) + (-1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) > 1)
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
# For actual fits
Temp.xs <- seq(1, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(5, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##########
###### 4. Fit bc thermal responses with uniform priors
##########

###################################### b for WNV Cx tarsalis - quadratic

##### Set data
data <- data.b.CtarWNV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
b.CtarWNV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarWNV.out$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarWNV.out)

save(b.CtarWNV.out, file = "jagsout_b_CtarWNV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarWNV, ylab = "b for WNV in Cx tarsalis", xlab = "Temperature")
lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### bc for WEEV Cx tarsalis - quadratic

##### Set data
data <- data.bc.CtarWEEV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.CtarWEEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.CtarWEEV.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.CtarWEEV.out)

save(bc.CtarWEEV.out, file = "jagsout_bc_CtarWEEV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CtarWEEV, ylab = "bc for WEEV in Cx tarsalis", xlab = "Temperature")
lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### c for WEEV Cx tarsalis - quadratic

##### Set data
data <- data.c.CtarWEEV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.CtarWEEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CtarWEEV.out$BUGSoutput$summary[1:5,]
mcmcplot(c.CtarWEEV.out)

save(c.CtarWEEV.out, file = "jagsout_c_CtarWEEV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CtarWEEV, ylab = "c for WEEV in Cx tarsalis", xlab = "Temperature")
lines(c.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### b for WEEV Cx tarsalis - quadratic

##### Set data
data <- data.b.CtarWEEV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
b.CtarWEEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarWEEV.out$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarWEEV.out)

save(b.CtarWEEV.out, file = "jagsout_b_CtarWEEV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarWEEV, ylab = "b for WEEV in Cx tarsalis", xlab = "Temperature")
lines(b.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### c for SLEV Cx tarsalis - quadratic

##### Set data
data <- data.c.CtarSLEV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.CtarSLEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CtarSLEV.out$BUGSoutput$summary[1:5,]
mcmcplot(c.CtarSLEV.out)

save(c.CtarSLEV.out, file = "jagsout_c_CtarSLEV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CtarSLEV, ylab = "c for SLEV in Cx tarsalis", xlab = "Temperature")
lines(c.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



###################################### b for SLEV Cx tarsalis - quadratic

##### Set data
data <- data.b.CtarSLEV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
b.CtarSLEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarSLEV.out$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarSLEV.out)

save(b.CtarSLEV.out, file = "jagsout_c_CtarSLEV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CtarSLEV, ylab = "b for SLEV in Cx tarsalis", xlab = "Temperature")
lines(b.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### bc for WNV in Cx. pipiens - quadratic

##### Set data
data <- data.bc.CpipWNV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.CpipWNV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.CpipWNV.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.CpipWNV.out)

save(bc.CpipWNV.out, file = "jagsout_bc_CpipWNV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CpipWNV, ylab = "bc for WNV in Cx. pipiens", xlab = "Temperature")
lines(bc.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### c for WNV in Cx. pipiens - quadratic

##### Set data
data <- data.c.CpipWNV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.CpipWNV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CpipWNV.out$BUGSoutput$summary[1:5,]
mcmcplot(c.CpipWNV.out)

save(c.CpipWNV.out, file = "jagsout_c_CpipWNV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CpipWNV, ylab = "c for WNV Cx pipiens", xlab = "Temperature")
lines(c.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### bc for WNV in Cx. univittatus - quadratic

##### Set data
data <- data.bc.CuniWNV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.CuniWNV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.CuniWNV.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.CuniWNV.out)

save(bc.CuniWNV.out, file = "jagsout_bc_CuniWNV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CuniWNV, ylab = "bc for WNV in Cx. univittatus", xlab = "Temperature")
lines(bc.CuniWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CuniWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CuniWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### bc for RVFV in Cx. pipiens - quadratic

##### Set data
data <- data.bc.AtaeRVFV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.AtaeRVFV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.AtaeRVFV.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.AtaeRVFV.out)

save(bc.AtaeRVFV.out, file = "jagsout_bc_CpipRVFV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.AtaeRVFV, ylab = "bc for RVFV Ae taeniorynchus", xlab = "Temperature")
lines(bc.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### bc for SINV in Ae. taeniorhynchus - quadratic

##### Set data
data <- data.c.AtaeSINV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.AtaeSINV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.AtaeSINV.out$BUGSoutput$summary[1:5,]
mcmcplot(c.AtaeSINV.out)

save(c.AtaeSINV.out, file = "jagsout_bc_AtaeSINV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.AtaeSINV, ylab = "c for SINV in Ae. taeniorhynchus", xlab = "Temperature")
lines(c.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### c for SINV in Cx. pipiens - quadratic

##### Set data
data <- data.c.CpipSINV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.CpipSINV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CpipSINV.out$BUGSoutput$summary[1:5,]
mcmcplot(c.CpipSINV.out)

save(c.CpipSINV.out, file = "jagsout_bc_CpipSINV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CpipSINV, ylab = "c for SINV in Cx. pipiens", xlab = "Temperature")
lines(c.CpipSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### bc for EEEV in Ae. triseriatus - quadratic

##### Set data
data <- data.bc.AtriEEEV

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.AtriEEEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.AtriEEEV.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.AtriEEEV.out)

save(bc.AtriEEEV.out, file = "jagsout_bc_AtriEEEV.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.AtriEEEV, ylab = "bc for EEEV in Ae. triseriatus", xlab = "Temperature")
lines(bc.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)




##########
###### 5. Fit thermal responses for bc priors using a leave-one-out approach
##########

# Subset Data
data.bc.only <- subset(data.bc, trait.name == "bc")
data.b.only <- subset(data.bc, trait.name == "b")
data.c.only <- subset(data.bc, trait.name == "c")
  
data.c.CpipWNV.prior <- subset(data.c.only, joint.code != "CpipWNV")
data.bc.CpipWNV.prior <- subset(data.bc.only, joint.code != "CpipWNV")

data.b.CtarWNV.prior <- subset(data.b.only, joint.code != "CtarWNV")
data.bc.CuniWNV.prior <- subset(data.bc.only, joint.code != "CuniWNV")

data.bc.CtarWEEV.prior <- subset(data.bc.only, joint.code != "CtarWEEV")
data.c.CtarWEEV.prior <- subset(data.c.only, joint.code != "CtarWEEV")
data.b.CtarWEEV.prior <- subset(data.b.only, joint.code != "CtarWEEV")
data.c.CtarSLEV.prior <- subset(data.c.only, joint.code != "CtarSLEV")
data.b.CtarSLEV.prior <- subset(data.b.only, joint.code != "CtarSLEV")

data.c.AtaeSINV.prior <- subset(data.c.only, joint.code != "AtaeSINV")
data.c.CpipSINV.prior <- subset(data.c.only, joint.code != "CpipSINV")

data.bc.AtaeRVFV.prior <- subset(data.bc.only, joint.code != "AtaeRVFV")
data.bc.AtriEEEV.prior <- subset(data.bc.only, joint.code != "AtriEEEV")


###################################### bc for WNV in Cpip prior

##### Set data
data <- data.bc.CpipWNV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.CpipWNV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.CpipWNV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.CpipWNV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CpipWNV.prior, ylab = "bc for WNV in Cpip prior", xlab = "Temperature")
lines(bc.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### c for WNV in Cpip prior

##### Set data
data <- data.c.CpipWNV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.CpipWNV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CpipWNV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(c.CpipWNV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CpipWNV.prior, ylab = "c for WNV in Cpip prior", xlab = "Temperature")
lines(c.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### b for WNV in Ctar prior

##### Set data
data <- data.b.CtarWNV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
b.CtarWNV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarWNV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarWNV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarWNV.prior, ylab = "b for WNV in Ctar prior", xlab = "Temperature")
lines(b.CtarWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### bc for WNV in Cuni prior

##### Set data
data <- data.bc.CuniWNV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.CuniWNV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.CuniWNV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.CuniWNV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CuniWNV.prior, ylab = "bc for WNV in Cuni prior", xlab = "Temperature")
lines(bc.CuniWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CuniWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CuniWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### bc for WEEV in Ctar prior

##### Set data
data <- data.bc.CtarWEEV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.CtarWEEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.CtarWEEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.CtarWEEV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CtarWEEV.prior, ylab = "bc for WEEV in Ctar prior", xlab = "Temperature")
lines(bc.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### c for WEEV in Ctar prior

##### Set data
data <- data.c.CtarWEEV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.CtarWEEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CtarWEEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(c.CtarWEEV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CtarWEEV.prior, ylab = "c for WEEV in Ctar prior", xlab = "Temperature")
lines(c.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### b for WEEV in Ctar prior

##### Set data
data <- data.b.CtarWEEV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
b.CtarWEEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarWEEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarWEEV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarWEEV.prior, ylab = "b for WEEV in Ctar prior", xlab = "Temperature")
lines(b.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### c for SLEV in Ctar prior

##### Set data
data <- data.c.CtarSLEV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.CtarSLEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CtarSLEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(c.CtarSLEV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CtarSLEV.prior, ylab = "c for SLEV in Ctar prior", xlab = "Temperature")
lines(c.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### b for SLEV in Ctar prior

##### Set data
data <- data.b.CtarSLEV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
b.CtarSLEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarSLEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarSLEV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarSLEV.prior, ylab = "b for SLEV in Ctar prior", xlab = "Temperature")
lines(b.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### c for SINV in Cpip prior

##### Set data
data <- data.c.CpipSINV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.CpipSINV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CpipSINV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(c.CpipSINV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CpipSINV.prior, ylab = "c for SINV in Cpip prior", xlab = "Temperature")
lines(c.CpipSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### c for SINV in Atae prior

##### Set data
data <- data.c.AtaeSINV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
c.AtaeSINV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.AtaeSINV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(c.AtaeSINV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.AtaeSINV.prior, ylab = "c for SINV in Atae prior", xlab = "Temperature")
lines(c.AtaeSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.AtaeSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.AtaeSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### bc for RVFV in Atae prior

##### Set data
data <- data.bc.AtaeRVFV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.AtaeRVFV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.AtaeRVFV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.AtaeRVFV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.AtaeRVFV.prior, ylab = "bc for RVFV in Atae prior", xlab = "Temperature")
lines(bc.AtaeRVFV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtaeRVFV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtaeRVFV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### bc for EEEV in Atri prior

##### Set data
data <- data.bc.AtriEEEV.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
bc.AtriEEEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.AtriEEEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(bc.AtriEEEV.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.AtriEEEV.prior, ylab = "bc for EEEV in Atri prior", xlab = "Temperature")
lines(bc.AtriEEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtriEEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtriEEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 6. Fit gamma distributions to bc prior thermal responses
##########

###################################### bc for WNV in Cpip prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
bc.CpipWNV.prior.cf.dists <- data.frame(q = as.vector(bc.CpipWNV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(bc.CpipWNV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(bc.CpipWNV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
bc.CpipWNV.prior.gamma.fits = apply(bc.CpipWNV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### c for WNV in Cpip prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
c.CpipWNV.prior.cf.dists <- data.frame(q = as.vector(c.CpipWNV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(c.CpipWNV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(c.CpipWNV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
c.CpipWNV.prior.gamma.fits = apply(c.CpipWNV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### b for WNV in Ctar prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
b.CtarWNV.prior.cf.dists <- data.frame(q = as.vector(b.CtarWNV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(b.CtarWNV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(b.CtarWNV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
b.CtarWNV.prior.gamma.fits = apply(b.CtarWNV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### bc for WNV in Cuni prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
bc.CuniWNV.prior.cf.dists <- data.frame(q = as.vector(bc.CuniWNV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(bc.CuniWNV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(bc.CuniWNV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
bc.CuniWNV.prior.gamma.fits = apply(bc.CuniWNV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### bc for WEEV in Ctar prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
bc.CtarWEEV.prior.cf.dists <- data.frame(q = as.vector(bc.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(bc.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(bc.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
bc.CtarWEEV.prior.gamma.fits = apply(bc.CtarWEEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### c for WEEV in Ctar prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
c.CtarWEEV.prior.cf.dists <- data.frame(q = as.vector(c.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(c.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(c.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
c.CtarWEEV.prior.gamma.fits = apply(c.CtarWEEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### b for WEEV in Ctar prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
b.CtarWEEV.prior.cf.dists <- data.frame(q = as.vector(b.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(b.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(b.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
b.CtarWEEV.prior.gamma.fits = apply(b.CtarWEEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### c for SLEV in Ctar prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
c.CtarSLEV.prior.cf.dists <- data.frame(q = as.vector(c.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(c.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(c.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
c.CtarSLEV.prior.gamma.fits = apply(c.CtarSLEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### b for SLEV in Ctar prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
b.CtarSLEV.prior.cf.dists <- data.frame(q = as.vector(b.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(b.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(b.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
b.CtarSLEV.prior.gamma.fits = apply(b.CtarSLEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### c for SINV in Cpip prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
c.CpipSINV.prior.cf.dists <- data.frame(q = as.vector(c.CpipSINV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(c.CpipSINV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(c.CpipSINV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
c.CpipSINV.prior.gamma.fits = apply(c.CpipSINV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### c for SINV in Atae prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
c.AtaeSINV.prior.cf.dists <- data.frame(q = as.vector(c.AtaeSINV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(c.AtaeSINV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(c.AtaeSINV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
c.AtaeSINV.prior.gamma.fits = apply(c.AtaeSINV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### bc for RVFV in Atae prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
bc.AtaeRVFV.prior.cf.dists <- data.frame(q = as.vector(bc.AtaeRVFV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(bc.AtaeRVFV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(bc.AtaeRVFV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
bc.AtaeRVFV.prior.gamma.fits = apply(bc.AtaeRVFV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### bc for EEEV in Atri prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
bc.AtriEEEV.prior.cf.dists <- data.frame(q = as.vector(bc.AtriEEEV.prior.out$BUGSoutput$sims.list$cf.q),
                                        T0 = as.vector(bc.AtriEEEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                        Tm = as.vector(bc.AtriEEEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
bc.AtriEEEV.prior.gamma.fits = apply(bc.AtriEEEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


bc.hypers <- list(c.CpipWNV.prior.gamma.fits, bc.CpipWNV.prior.gamma.fits, 
                  b.CtarWNV.prior.gamma.fits, bc.CuniWNV.prior.gamma.fits,
                  c.CtarWEEV.prior.gamma.fits, b.CtarWEEV.prior.gamma.fits, bc.CtarWEEV.prior.gamma.fits, 
                  c.CtarSLEV.prior.gamma.fits, b.CtarSLEV.prior.gamma.fits,
                  bc.AtriEEEV.prior.gamma.fits, bc.AtaeRVFV.prior.gamma.fits,
                  c.AtaeSINV.prior.gamma.fits, c.CpipSINV.prior.gamma.fits)
save(bc.hypers, file = "bchypers.Rsave")


##########
###### 7. Fit bc thermal responses with data-informed priors
##########

load("bchypers.Rsave")
c.CpipWNV.prior.gamma.fits <- bc.hypers[[1]]
bc.CpipWNV.prior.gamma.fits <- bc.hypers[[2]]
b.CtarWNV.prior.gamma.fits <- bc.hypers[[3]]
bc.CuniWNV.prior.gamma.fits <- bc.hypers[[4]]
c.CtarWEEV.prior.gamma.fits <- bc.hypers[[5]]
b.CtarWEEV.prior.gamma.fits <- bc.hypers[[6]]
bc.CtarWEEV.prior.gamma.fits <- bc.hypers[[7]]
c.CtarSLEV.prior.gamma.fits <- bc.hypers[[8]]
b.CtarSLEV.prior.gamma.fits <- bc.hypers[[9]]
bc.AtriEEEV.prior.gamma.fits <- bc.hypers[[10]]
bc.AtaeRVFV.prior.gamma.fits <- bc.hypers[[11]]
c.AtaeSINV.prior.gamma.fits <- bc.hypers[[12]]
c.CpipSINV.prior.gamma.fits <- bc.hypers[[13]]

###################################### b for WNV in Cx tarsalis - quadratic

##### Set data
data <- data.b.CtarWNV
hypers <- b.CtarWNV.prior.gamma.fits * 1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
b.CtarWNV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarWNV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarWNV.out.inf)

save(b.CtarWNV.out.inf, file = "jagsout_b_CtarWNV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarWNV, ylab = "b for WNV in Cx tarsalis", xlab = "Temperature")
lines(b.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

# Get optimum for b: 26.4 C
Temp.xs[which.max(as.vector(b.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


###################################### bc for WEEV Cx tarsalis - quadratic

##### Set data
data <- data.bc.CtarWEEV
hypers <- bc.CtarWEEV.prior.gamma.fits * 0.5

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
bc.CtarWEEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine infput
bc.CtarWEEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(bc.CtarWEEV.out.inf)

save(bc.CtarWEEV.out.inf, file = "jagsout_bc_CtarWEEV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CtarWEEV, ylab = "bc for WEEV in Cx tarsalis", xlab = "Temperature")
lines(bc.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)

# Get optimum for bc: 21.4 C
Temp.xs[which.max(as.vector(bc.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### c for WEEV Cx tarsalis - quadratic

##### Set data
data <- data.c.CtarWEEV
hypers <- c.CtarWEEV.prior.gamma.fits * .01

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
c.CtarWEEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CtarWEEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(c.CtarWEEV.out.inf)

save(c.CtarWEEV.out.inf, file = "jagsout_c_CtarWEEV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CtarWEEV, ylab = "c for WEEV in Cx tarsalis", xlab = "Temperature")
lines(c.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)

# Get optimum for c: 20.6 C
Temp.xs[which.max(as.vector(c.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))] # mean captures it better than median because it collapses at 1

###################################### b for WEEV Cx tarsalis - quadratic

##### Set data
data <- data.b.CtarWEEV
hypers <- b.CtarWEEV.prior.gamma.fits * 0.1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
b.CtarWEEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarWEEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarWEEV.out.inf)

save(b.CtarWEEV.out.inf, file = "jagsout_b_CtarWEEV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarWEEV, ylab = "b for WEEV in Cx tarsalis", xlab = "Temperature")
lines(b.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

# Get optimum for b: 21.0 C
Temp.xs[which.max(as.vector(b.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### c for SLEV Cx tarsalis - quadratic

##### Set data
data <- data.c.CtarSLEV
hypers <- c.CtarSLEV.prior.gamma.fits * 0.01
hypers[,3] <- c.CtarSLEV.prior.gamma.fits[,3] * 0.1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
c.CtarSLEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CtarSLEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(c.CtarSLEV.out.inf)

save(c.CtarSLEV.out.inf, file = "jagsout_c_CtarSLEV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CtarSLEV, ylab = "c for SLEV in Cx tarsalis", xlab = "Temperature")
lines(c.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.2)

# Get optimum for c: 26.2 C
Temp.xs[which.max(as.vector(c.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### b for SLEV Cx tarsalis - quadratic

##### Set data
data <- data.b.CtarSLEV
hypers <- b.CtarSLEV.prior.gamma.fits * 0.5

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
b.CtarSLEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
b.CtarSLEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(b.CtarSLEV.out.inf)

save(b.CtarSLEV.out.inf, file = "jagsout_b_CtarSLEV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarSLEV, ylab = "b for SLEV in Cx tarsalis", xlab = "Temperature")
lines(b.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(b.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for b: 26.2 C
Temp.xs[which.max(as.vector(b.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### bc for WNV in Cx. pipiens - quadratic

##### Set data
data <- data.bc.CpipWNV
hypers <- bc.CpipWNV.prior.gamma.fits * .5

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
bc.CpipWNV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.CpipWNV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(bc.CpipWNV.out.inf)

save(bc.CpipWNV.out.inf, file = "jagsout_bc_CpipWNV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CpipWNV, ylab = "bc for WNV in Cx. pipiens", xlab = "Temperature")
lines(bc.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

# Get optimum for bc: 27.9 C
Temp.xs[which.max(as.vector(bc.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### c for WNV in Cx. pipiens - quadratic

##### Set data
data <- data.c.CpipWNV
hypers <- c.CpipWNV.prior.gamma.fits * 1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
c.CpipWNV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CpipWNV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(c.CpipWNV.out.inf)

save(c.CpipWNV.out.inf, file = "jagsout_c_CpipWNV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CpipWNV, ylab = "c for WNV Cx pipiens", xlab = "Temperature")
lines(c.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

# Get optimum for c: 33.6 C
Temp.xs[which.max(as.vector(c.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### bc for WNV in Cx. univittatus - quadratic

##### Set data
data <- data.bc.CuniWNV
hypers <- bc.CuniWNV.prior.gamma.fits * .01
hypers[,3] <- c.CtarSLEV.prior.gamma.fits[,3] * 0.1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
bc.CuniWNV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.CuniWNV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(bc.CuniWNV.out.inf)

save(bc.CuniWNV.out.inf, file = "jagsout_bc_CuniWNV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CuniWNV, ylab = "bc for WNV in Cx. univittatus", xlab = "Temperature")
lines(bc.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)

# Get optimum for bc: 24.6 C
Temp.xs[which.max(as.vector(bc.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### bc for RVFV in Ae. taeniorhynchus - quadratic

##### Set data
data <- data.bc.AtaeRVFV
hypers <- bc.AtaeRVFV.prior.gamma.fits * 2

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
bc.AtaeRVFV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.AtaeRVFV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(bc.AtaeRVFV.out.inf)

save(bc.AtaeRVFV.out.inf, file = "jagsout_bc_AtaeRVFV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.AtaeRVFV, ylab = "bc for RVFV Ae taeniorynchus", xlab = "Temperature")
lines(bc.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "G", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for bc: 24.7 C
Temp.xs[which.max(as.vector(bc.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### c for SINV in Ae. taeniorhynchus - quadratic

##### Set data
data <- data.c.AtaeSINV
hypers <- c.AtaeSINV.prior.gamma.fits * .1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
c.AtaeSINV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.AtaeSINV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(c.AtaeSINV.out.inf)

save(c.AtaeSINV.out.inf, file = "jagsout_c_AtaeSINV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.AtaeSINV, ylab = "c for SINV in Ae. taeniorhynchus", xlab = "Temperature")
lines(c.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.2)

# Get optimum for c: 25.2 C
Temp.xs[which.max(as.vector(c.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### c for SINV in Cx. pipiens - quadratic

##### Set data
data <- data.c.CpipSINV
hypers <- c.CpipSINV.prior.gamma.fits * 0.01

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
c.CpipSINV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
c.CpipSINV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(c.CpipSINV.out.inf)

save(c.CpipSINV.out.inf, file = "jagsout_c_CpipSINV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CpipSINV, ylab = "c for SINV in Cx. pipiens", xlab = "Temperature")
lines(c.CpipSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(c.CpipSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#lines(c.CpipSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"] ~ Temp.xs, col = "red")
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.2)

# Get optimum for c: 17.2 C
Temp.xs[which.max(as.vector(c.CpipSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

###################################### bc for EEEV in Ae. triseriatus - quadratic

##### Set data
data <- data.bc.AtriEEEV
hypers <- bc.AtriEEEV.prior.gamma.fits * 3
hypers[,3] <- bc.AtriEEEV.prior.gamma.fits[,3] * .01

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
bc.AtriEEEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
bc.AtriEEEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(bc.AtriEEEV.out.inf)

save(bc.AtriEEEV.out.inf, file = "jagsout_bc_AtriEEEV_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.AtriEEEV, ylab = "bc for EEEV in Ae. triseriatus", xlab = "Temperature")
lines(bc.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)

# Get optimum for bc: 28.8 C
Temp.xs[which.max(as.vector(bc.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


##########
###### 8. Calculate Treatment averages for plotting
##########

# Function to trait averages for plotting only (JAGS models fit to raw data)
CalcTraitAvg = function(data.set){
  data.set <- data.set
  temp.list <- unique(data.set$T)
  out <- data.frame(Temp = numeric(length(temp.list)), mean = numeric(length(temp.list)), sd = numeric(length(temp.list)), SE = numeric(length(temp.list)), 
                    n = numeric(length(temp.list)), upper = numeric(length(temp.list)), lower = numeric(length(temp.list)))
  for(i in 1:length(temp.list)){
    data.sub <- subset(data.set, T == temp.list[i])
    out$Temp[i] <- temp.list[i]
    out$mean[i] <- mean(data.sub$trait)
    out$sd[i] <- sd(data.sub$trait)
    out$n[i] <- nrow(data.sub)
    out$SE[i] <- sd(data.sub$trait) / sqrt(nrow(data.sub))
    out$lower[i] <- mean(data.sub$trait) - sd(data.sub$trait) / sqrt(nrow(data.sub))
    out$upper[i] <- mean(data.sub$trait) + sd(data.sub$trait) / sqrt(nrow(data.sub))
  }
  out
}


CpipWNV.c.avg <- CalcTraitAvg(data.c.CpipWNV)
CpipWNV.bc.avg <- CalcTraitAvg(data.bc.CpipWNV)

CtarWNV.b.avg <- CalcTraitAvg(data.b.CtarWNV)
CuniWNV.bc.avg <- CalcTraitAvg(data.bc.CuniWNV)

CtarWEEV.c.avg <- CalcTraitAvg(data.c.CtarWEEV)
CtarWEEV.b.avg <- CalcTraitAvg(data.b.CtarWEEV)
CtarWEEV.bc.avg <- CalcTraitAvg(data.bc.CtarWEEV)
CtarSLEV.c.avg <- CalcTraitAvg(data.c.CtarSLEV)
CtarSLEV.b.avg <- CalcTraitAvg(data.b.CtarSLEV)

CpipSINV.c.avg <- CalcTraitAvg(data.c.CpipSINV)
AtaeSINV.c.avg <- CalcTraitAvg(data.c.AtaeSINV)
AtaeRVFV.bc.avg <- CalcTraitAvg(data.bc.AtaeRVFV)
AtriEEEV.bc.avg <- CalcTraitAvg(data.bc.AtriEEEV)


plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.c.avg, ylab = "c for WNV in Cpip", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipWNV.bc.avg, ylab = "bc for WNv in Cpip", xlab = "Temperature")

plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CtarWNV.c.avg, ylab = "c for WNV in Cpip", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CtarWNV.b.avg, ylab = "b for WNV in Cpip", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CuniWNV.bc.avg, ylab = "bc for WNV in Cuni", xlab = "Temperature")

plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CtarWEEV.c.avg, ylab = "c for WEEV in Ctar", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CtarWEEV.b.avg, ylab = "b for WEEV in Ctar", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CtarWEEV.bc.avg, ylab = "bc for WEEV in Ctar", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CtarSLEV.c.avg, ylab = "c for SLEV in Ctar", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CtarSLEV.b.avg, ylab = "b for SLEV in Ctar", xlab = "Temperature")

plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = CpipSINV.c.avg, ylab = "c for SINV in Cpip", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = AtaeSINV.c.avg, ylab = "c for SINV in Atae", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = AtaeRVFV.bc.avg, ylab = "bc for RVFV in Atae", xlab = "Temperature")
plot(mean ~ Temp, xlim = c(5, 45), ylim = c(0,1), data = AtriEEEV.bc.avg, ylab = "bc for EEEV in Atri", xlab = "Temperature")


##########
###### 9. Plot Preliminary Figures
##########

plot(trait ~ T, xlim = c(5, 45), data = data.bc, ylab = "bc for Cx. mosquitoes", xlab = "Temperature")
points(trait ~ T, data = data.bc.Cpipmol, col = "blue")
points(trait ~ T, data = data.bc.Cpippip, col = "dodgerblue")
points(trait ~ T, data = data.bc.Cqui, col = "red")
points(trait ~ T, data = data.bc.Ctar, col = "darkgreen")

lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "dodgerblue")
lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "dodgerblue")
lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "dodgerblue")

lines(c.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(c.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(c.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "blue")

lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(bc.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(bc.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(bc.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

lines(bc.CquiWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "darkgreen")
lines(bc.CquiWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "darkgreen")
lines(bc.CquiWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "darkgreen")


################### Figures with original data points

par(mfrow = c(3,3), mar = c(3, 4.5, 2, 1), oma = c(2, 0, 0, 0))

##### c for WNV in Cx. pipiens
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CpipWNV, xaxt = "n", pch = 19,
     ylab = "Infection Probability (c)", xlab = "", main = expression(paste("WNV in ",italic(Cx.)," ",italic(pipiens))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(c.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(c.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(c.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### b for WNV in Cx. pipiens
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CpipWNV, xaxt = "n", pch = 19,
     ylab = "Transmission Probability (b)", xlab = "", main = expression(paste("WNV in ",italic(Cx.)," ",italic(pipiens))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(b.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(b.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(b.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### c for WNV in Cx. tarsalis
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.c.CtarWNV, xaxt = "n", pch = 19,
     ylab = "Infection Probability (c)", xlab = "", main = expression(paste("WNV in ",italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(c.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(c.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(c.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### b for WNV in Cx. tarsalis
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.b.CtarWNV, xaxt = "n", pch = 19,
     ylab = "Transmission Probability (b)", xlab = "", main = expression(paste("WNV in ",italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(b.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### bc for WEEV in Cx. tarsalis
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CtarWEEV, xaxt = "n", pch = 19,
     ylab = "Transmission Probability (bc)", xlab = "", main = expression(paste("WEEV in ",italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### bc for SLEV in Cx. tarsalis
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CtarSLEV, xaxt = "n", pch = 19,
     ylab = "Transmission Probability (bc)", xlab = "", main = expression(paste("SLEV in ",italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(bc.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### bc for RVFV in Ae. taeniorhynchus
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.AtaeRVFV, xaxt = "n", pch = 19,
     ylab = "Transmission Probability (bc)", xlab = "", main = expression(paste("RVFV in ",italic(Cx.)," ",italic(pipiens))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(bc.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### bc for EEEV in Ae. triseriatus
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.AtriEEEV, xaxt = "n", pch = 19,
     ylab = "Transmission Probability (bc)", xlab = "", main = expression(paste("EEEV in ",italic(Ae.)," ",italic(triseriatus))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(bc.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### bc for SINV in Cx. pipiens
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.bc.CpipSINV, xaxt = "n", pch = 19,
     ylab = "Transmission Probability (bc)", xlab = "", main = expression(paste("SINV in ",italic(Cx.)," ",italic(pipiens))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 0.9)
