## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for fecundity related traits for many Culex and Aedes species
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit EFD thermal responses with uniform priors
##           5) Fit EFD thermal responses for priors
##           6) Fit gamma distributions to EFD prior thermal responses
##           7) Fit EFD thermal responses with data-informed priors
##           8) Calculate treatment averages for plotting


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("~/Fitting Traits")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')
library('MASS')

# Load Data
data.E <- read.csv("TraitData_EFD.csv")

# Subset Data by trait
data.EFOC.Cpip <- subset(data.E, Trait.Name == "EFOC")
data.EPR <- subset(data.E, Trait.Name == "EPR")
data.pO <- subset(data.E, Trait.Name == "pO")

# Subset Data
data.EPR.Cqui <- subset(data.EPR, host.code == "Cqui")
data.EPR.Cmol <- subset(data.EPR, host.code == "Cmol")
data.EPR.Cpal <- subset(data.EPR, host.code == "Cpal")
data.EPR.Ador <- subset(data.EPR, host.code == "Ador")

# Subset Data
data.pO.Cpip <- subset(data.pO, host.code == "Cpip")
data.pO.Cqui <- subset(data.pO, host.code == "Cqui")
data.pO.Cmol <- subset(data.pO, host.code == "Cmol")
data.pO.Cpal <- subset(data.pO, host.code == "Cpal")
data.pO.Ador <- subset(data.pO, host.code == "Ador")
data.pO.Cmel <- subset(data.pO, host.code == "Cmel") # This data was added later, not included in priors for other species

par(mfrow = c(1,1))
plot(trait ~ T, xlim = c(5, 45), data = data.EFOC, ylab = "EFOC for Cx. pipiens", xlab = "Temperature")

plot(trait ~ T, xlim = c(5, 45), data = data.EPR, ylab = "EPR for Cx. mosquitoes", xlab = "Temperature")
points(trait ~ T, data = data.EPR.Cqui, col = "red")
points(trait ~ T, data = data.EPR.Cmol, col = "violet")
points(trait ~ T, data = data.EPR.Cpal, col = "darkorange")
points(trait ~ T, data = data.EPR.Ador, col = "darkgreen")

plot(trait ~ T, xlim = c(5, 45), data = data.pO, ylab = "pO for Cx. mosquitoes", xlab = "Temperature")
points(trait ~ T, data = data.pO.Cpip, col = "grey")
points(trait ~ T, data = data.pO.Cqui, col = "red")
points(trait ~ T, data = data.pO.Cmol, col = "violet")
points(trait ~ T, data = data.pO.Ador, col = "darkgreen")


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
    cf.Tm ~ dunif(25, 50)
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

############## Quadratic Model with gamma priors

sink("quad_inf.txt")
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
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


############## Briere Model with gamma priors (except sigma) and predicted values <= 1 (for probabilities )

sink("briereprob.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 50)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i]) * (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) < 1) + (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) > 1)
    }
    
    } # close model
    ",fill=T)
sink()

############## Briere Model with gamma priors (except sigma) and predicted values <= 1 (for probabilities )

sink("briereprob_inf.txt")
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
    trait.mu[i] <- cf.q * temp[i] * (temp[i] - cf.T0) * sqrt((cf.Tm - temp[i]) * (cf.Tm > temp[i])) * (cf.T0 < temp[i])
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i]) * (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) < 1) + (cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) > 1)
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
###### 4. Fit EFD thermal responses with uniform priors
##########

###################################### EFD for Cx pipiens - quadratic

##### Set data
data <- data.EFOC.Cpip

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EFOC.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EFOC.Cpip.out$BUGSoutput$summary[1:5,]
mcmcplot(EFOC.Cpip.out)

save(EFOC.Cpip.out, file = "jagsout_EFOC_Cpip.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,250), data = data.EFOC.Cpip, ylab = "EFOC for Cx pipiens", xlab = "Temperature")
lines(EFOC.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFOC.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFOC.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### EPR for Cx quinquefasciatus - quadratic

##### Set data
data <- data.EPR.Cqui

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EPR.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EPR.Cqui.out$BUGSoutput$summary[1:5,]
mcmcplot(EPR.Cqui.out)

save(EPR.Cqui.out, file = "jagsout_EPR_Cqui.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,300), data = data.EPR.Cqui, ylab = "EPR for Cx quinquefasciatus", xlab = "Temperature")
lines(EPR.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EPR.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EPR.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pO for Cx quinquefasciatus - briere

##### Set data
data <- data.pO.Cqui

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pO.Cqui.out.b <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briereprob.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
pO.Cqui.out.q <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
pO.Cqui.out.b$BUGSoutput$DIC # briere is better
pO.Cqui.out.q$BUGSoutput$DIC

##### Examine Output
pO.Cqui.out$BUGSoutput$summary[1:5,]
mcmcplot(pO.Cqui.out)

save(pO.Cqui.out, file = "jagsout_pO_Cqui.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cqui, ylab = "pO for Cx quinquefasciatus", xlab = "Temperature")
lines(pO.Cqui.out.q$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cqui.out.q$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cqui.out.q$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pO.Cqui.out.b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(pO.Cqui.out.b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(pO.Cqui.out.b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "blue")


###################################### pO for Cx pipiens - quadratic

##### Set data
data <- data.pO.Cpip

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pO.Cpip.out.b <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
pO.Cpip.out.q <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                      n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
pO.Cpip.out.b$BUGSoutput$DIC
pO.Cpip.out.q$BUGSoutput$DIC  # quadratic is slightly better

##### Examine Output
pO.Cpip.out$BUGSoutput$summary[1:5,]
mcmcplot(pO.Cpip.out)

save(pO.Cpip.out, file = "jagsout_pO_Cpip.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cpip, ylab = "pO for Cx pipiens", xlab = "Temperature")
lines(pO.Cpip.out.q$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cpip.out.q$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cpip.out.q$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
lines(pO.Cpip.out.b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(pO.Cpip.out.b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(pO.Cpip.out.b$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "blue")


##########
###### 5. Fit thermal responses for EFD priors
##########

# Subset Data - leave-one-out
data.EFOC.Cpip.prior <- data.EPR
data.EPR.Cqui.prior <- rbind(subset(data.EPR, host.code != "Cqui"), data.EFOC.Cpip)
data.pO.Cpip.prior <- subset(data.pO, host.code != "Cpip")
data.pO.Cqui.prior <- subset(data.pO, host.code != "Cqui")
data.pO.Cqui.prior <- subset(data.pO, host.code != "Cqui")
data.pO.Cmel.prior <- subset(data.pO, host.code != "Cmel")

###################################### EFOC for Cpip prior - quadratic

##### Set data
data <- data.EFOC.Cpip.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EFOC.Cpip.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EFOC.Cpip.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(EFOC.Cpip.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,300), data = data.EFOC.Cpip.prior, ylab = "EFOC for Cx. pipiens prior", xlab = "Temperature")
lines(EFOC.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFOC.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFOC.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### EPR for Cqui prior - quadratic

##### Set data
data <- data.EPR.Cqui.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EPR.Cqui.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EPR.Cqui.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(EPR.Cqui.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,300), data = data.EPR.Cqui.prior, ylab = "EPR for Cx quinquefasciatus prior", xlab = "Temperature")
lines(EPR.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EPR.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EPR.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pO for Cqui prior - quadratic

##### Set data
data <- data.pO.Cqui.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pO.Cqui.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briereprob.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pO.Cqui.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pO.Cqui.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cqui.prior, ylab = "pO for Cx quinquefasciatus prior", xlab = "Temperature")
lines(pO.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pO for Cpip prior - quadratic

##### Set data
data <- data.pO.Cpip.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pO.Cpip.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pO.Cpip.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pO.Cpip.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cpip.prior, ylab = "pO for Cx pipiens prior", xlab = "Temperature")
lines(pO.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pO for Cmel prior - quadratic

##### Set data
data <- data.pO.Cmel.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pO.Cmel.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pO.Cmel.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pO.Cmel.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cmel.prior, ylab = "pO for Cx pipiens prior", xlab = "Temperature")
lines(pO.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



##########
###### 6. Fit gamma distributions to EFD prior thermal responses
##########

###################################### EFOC for Cpip prior

# Get the posterior dists for all 4 parameters into a data frame
EFOC.Cpip.prior.cf.dists <- data.frame(q = as.vector(EFOC.Cpip.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(EFOC.Cpip.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(EFOC.Cpip.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EFOC.Cpip.prior.gamma.fits = apply(EFOC.Cpip.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### EPR for Cqui prior

# Get the posterior dists for all 4 parameters into a data frame
EPR.Cqui.prior.cf.dists <- data.frame(q = as.vector(EPR.Cqui.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(EPR.Cqui.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(EPR.Cqui.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EPR.Cqui.prior.gamma.fits = apply(EPR.Cqui.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pO for Cpip prior

# Get the posterior dists for all 4 parameters into a data frame
pO.Cpip.prior.cf.dists <- data.frame(q = as.vector(pO.Cpip.prior.out$BUGSoutput$sims.list$cf.q),
                                       T0 = as.vector(pO.Cpip.prior.out$BUGSoutput$sims.list$cf.T0), 
                                       Tm = as.vector(pO.Cpip.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pO.Cpip.prior.gamma.fits = apply(pO.Cpip.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pO for Cqui prior

# Get the posterior dists for all 4 parameters into a data frame
pO.Cqui.prior.cf.dists <- data.frame(q = as.vector(pO.Cqui.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pO.Cqui.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pO.Cqui.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pO.Cqui.prior.gamma.fits = apply(pO.Cqui.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pO for Cmel prior

# Get the posterior dists for all 4 parameters into a data frame
pO.Cmel.prior.cf.dists <- data.frame(q = as.vector(pO.Cmel.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(pO.Cmel.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(pO.Cmel.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pO.Cmel.prior.gamma.fits = apply(pO.Cmel.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


E.hypers <- list(EFOC.Cpip.prior.gamma.fits, EPR.Cqui.prior.gamma.fits, 
                 pO.Cpip.prior.gamma.fits, pO.Cqui.prior.gamma.fits, pO.Cmel.prior.gamma.fits)
save(E.hypers, file = "Ehypers.Rsave")


##########
###### 7. Fit EFD thermal responses with data-informed priors
##########

load("Ehypers.Rsave")
EFOC.Cpip.prior.gamma.fits <- E.hypers[[1]]
EPR.Cqui.prior.gamma.fits <- E.hypers[[2]]
pO.Cpip.prior.gamma.fits <- E.hypers[[3]]
pO.Cqui.prior.gamma.fits <- E.hypers[[4]]
pO.Cmel.prior.gamma.fits <- E.hypers[[5]]

############## EFOC for Cx. pipiens - quadratic

##### Set data 
data <- data.EFOC.Cpip
hypers <- EFOC.Cpip.prior.gamma.fits * 3

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
EFOC.Cpip.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EFOC.Cpip.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(EFOC.Cpip.out.inf)

save(EFOC.Cpip.out.inf, file = "jagsout_EFOC_Cpip_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,320), data = data.EFOC.Cpip, ylab = "EFGC for Cx. pipiens ", xlab = "Temperature", pch = 1)
lines(EFOC.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFOC.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFOC.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

# Get optimum for EFOC: 22.2 C
Temp.xs[which.max(as.vector(EFOC.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## EPR for Cx. quinquefasciatus - quadratic

##### Set data 
data <- data.EPR.Cqui
hypers <- EPR.Cqui.prior.gamma.fits * 1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
EPR.Cqui.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EPR.Cqui.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(EPR.Cqui.out.inf)

save(EPR.Cqui.out.inf, file = "jagsout_EPR_Cqui_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,320), data = data.EPR.Cqui, ylab = "ER for Cx. quinquefasciatus", xlab = "Temperature", pch = 1)
lines(EPR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EPR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EPR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

# Get optimum for EPR: 21.5 C
Temp.xs[which.max(as.vector(EPR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## pO for Cx. pipiens - quadratic

##### Set data 
data <- data.pO.Cpip
hypers <- pO.Cpip.prior.gamma.fits * 0.5

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pO.Cpip.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pO.Cpip.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pO.Cpip.out.inf)

save(pO.Cpip.out.inf, file = "jagsout_pO_Cpip_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cpip, ylab = "pO for Cx. pipiens ", xlab = "Temperature", pch = 1)
lines(pO.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)

# Get optimum for pO: 20.8 C
Temp.xs[which.max(as.vector(pO.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## pO for Cx. quinquefasciatus - quadratic

##### Set data 
data <- data.pO.Cqui
hypers <- pO.Cqui.prior.gamma.fits * 1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pO.Cqui.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briereprob_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pO.Cqui.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pO.Cqui.out.inf)

save(pO.Cqui.out.inf, file = "jagsout_pO_Cqui_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cqui, ylab = "pO for Cx. quinquefasciatus", xlab = "Temperature", pch = 1)
lines(pO.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)

# Get optimum for pO: 24.9 C
Temp.xs[which.max(as.vector(pO.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## pO for Cs. melanura - quadratic

##### Set data 
data <- data.pO.Cmel
hypers <- pO.Cmel.prior.gamma.fits * 1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pO.Cmel.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pO.Cmel.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pO.Cmel.out.inf)

save(pO.Cmel.out.inf, file = "jagsout_pO_Cmel_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pO.Cmel, ylab = "pO for Cs. melanura ", xlab = "Temperature", pch = 1)
lines(pO.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pO.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.2)

# Get optimum for pO: 21.1 C
Temp.xs[which.max(as.vector(pO.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


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

Cqui.EPR.avg <- CalcTraitAvg(data.EPR.Cqui)
Cpip.pO.avg <- CalcTraitAvg(data.pO.Cpip)
Cqui.pO.avg <- CalcTraitAvg(data.pO.Cqui)


#########
##### Figures
########

## See a (biting rate) R script for joint preliminary figures