## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for juvenile survival (larval-to-adult, pLA) for many Culex and Aedes species
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit pLA thermal responses with uniform priors
##           5) Fit pLA thermal responses for priors
##           6) Fit gamma distributions to pLA prior thermal responses
##           7) Fit pLA thermal responses with data-informed priors
##           8) Calculate Treatment averages for plotting
##           9) Plot preliminary figures

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
data.pLA <- read.csv("TraitData_pLA.csv") # Data from database for most traits (except below)
unique(data.pLA$host.code)

# Subset Data
data.pLA.Cpip <- subset(data.pLA, host.code == "Cpip")
data.pLA.Cqui <- subset(data.pLA, host.code == "Cqui")
data.pLA.Ctar <- subset(data.pLA, host.code == "Ctar")
data.pLA.Cmol <- subset(data.pLA, host.code == "Cmol")
data.pLA.Cpal <- subset(data.pLA, host.code == "Cpal")
data.pLA.Cres <- subset(data.pLA, host.code == "Cres")
data.pLA.Csal <- subset(data.pLA, host.code == "Csal")
data.pLA.Atri <- subset(data.pLA, host.code == "Atri")
data.pLA.Avex <- subset(data.pLA, host.code == "Avex")
data.pLA.Asol <- subset(data.pLA, host.code == "Asol")
data.pLA.Anig <- subset(data.pLA, host.code == "Anig")
data.pLA.Cmel <- subset(data.pLA, host.code == "Cmel") # This data was added later, not included in priors for other species

par(mfrow = c(1,1))
plot(trait ~ T, xlim = c(5, 45), data = data.pLA, ylab = "pLA for Cx. mosquitoes", xlab = "Temperature")
points(trait ~ T, data = data.pLA.Cpip, col = "grey")
points(trait ~ T, data = data.pLA.Cqui, col = "red")
points(trait ~ T, data = data.pLA.Ctar, col = "blue")
points(trait ~ T, data = data.pLA.Cmol, col = "violet")
points(trait ~ T, data = data.pLA.Cpal, col = "darkorange")
points(trait ~ T, data = data.pLA.Cres, col = "violetred")
points(trait ~ T, data = data.pLA.Csal, col = "goldenrod2")
points(trait ~ T, data = data.pLA.Atri, col = "darkgreen")
points(trait ~ T, data = data.pLA.Avex, col = "cyan3")
points(trait ~ T, data = data.pLA.Asol, col = "green")
points(trait ~ T, data = data.pLA.Anig, col = "seagreen")


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


############## Quadratic Model with gamma priors (except sigma)

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
###### 4. Fit pLA thermal responses with uniform priors
##########

###################################### pLA for Cx pipiens - quadratic

##### Set data
data <- data.pLA.Cpip

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cpip.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cpip.out)

save(pLA.Cpip.out, file = "jagsout_pLA_Cpip.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cpip, ylab = "pLA for Cx pipiens", xlab = "Temperature")
lines(pLA.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Cx pipiens molestus - quadratic

##### Set data
data <- data.pLA.Cpipmol

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cpipmol.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cpipmol.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cpipmol.out)

save(pLA.Cpipmol.out, file = "jagsout_pLA_Cpipmol.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cpipmol, ylab = "pLA for Cx pipiens molestus", xlab = "Temperature")
lines(pLA.Cpipmol.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpipmol.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpipmol.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Cx quinquefasciatus - quadratic

##### Set data
data <- data.pLA.Cqui

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cqui.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cqui.out)

save(pLA.Cqui.out, file = "jagsout_pLA_Cqui.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cqui, ylab = "pLA for Cx quinquefasciatus", xlab = "Temperature")
lines(pLA.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Cx tarsalis - quadratic

##### Set data
data <- data.pLA.Ctar

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Ctar.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Ctar.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Ctar.out)

save(pLA.Ctar.out, file = "jagsout_pLA_Ctar.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Ctar, ylab = "pLA for Cx tarsalis", xlab = "Temperature")
lines(pLA.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pLA for Ae. triseriatus - quadratic

##### Set data
data <- data.pLA.Atri

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Atri.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Atri.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Atri.out)

save(pLA.Atri.out, file = "jagsout_pLA_Atri.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Atri, ylab = "pLA for Ae triseriatus", xlab = "Temperature")
lines(pLA.Atri.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Atri.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Atri.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pLA for Ae. vexans - quadratic

##### Set data
data <- data.pLA.Avex

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Avex.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Avex.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Avex.out)

save(pLA.Avex.out, file = "jagsout_pLA_Avex.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Avex, ylab = "pLA for Ae vexans", xlab = "Temperature")
lines(pLA.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 5. Fit thermal responses for pLA priors using a leave-one-out approach
##########

# Subset Data - leave-one-out
data.pLA.Cpip.prior <- subset(data.pLA, host.code != "Cpip")
data.pLA.Cqui.prior <- subset(data.pLA, host.code != "Cqui")
data.pLA.Ctar.prior <- subset(data.pLA, host.code != "Ctar")
data.pLA.Cmol.prior <- subset(data.pLA, host.code != "Cmol")
data.pLA.Cpal.prior <- subset(data.pLA, host.code != "Cpal")
data.pLA.Cres.prior <- subset(data.pLA, host.code != "Cres")
data.pLA.Csal.prior <- subset(data.pLA, host.code != "Csal")
data.pLA.Atri.prior <- subset(data.pLA, host.code != "Atri")
data.pLA.Avex.prior <- subset(data.pLA, host.code != "Avex")
data.pLA.Asol.prior <- subset(data.pLA, host.code != "Asol")
data.pLA.Anig.prior <- subset(data.pLA, host.code != "Anig")
data.pLA.Cmel.prior <- subset(data.pLA, host.code != "Cmel")

###################################### pLA for Cpip prior - quadratic

##### Set data
data <- data.pLA.Cpip.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cpip.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cpip.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cpip.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cpip.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Cqui prior - quadratic

##### Set data
data <- data.pLA.Cqui.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cqui.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cqui.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cqui.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cqui.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Ctar prior - quadratic

##### Set data
data <- data.pLA.Ctar.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Ctar.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Ctar.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Ctar.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Ctar.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Cmol prior - quadratic

##### Set data
data <- data.pLA.Cmol.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cmol.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cmol.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cmol.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cmol.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Cmol.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cmol.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cmol.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Cpal prior - quadratic

##### Set data
data <- data.pLA.Cpal.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cpal.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cpal.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cpal.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cpal.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Cpal.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpal.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpal.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pLA for Cres prior - quadratic

##### Set data
data <- data.pLA.Cres.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cres.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cres.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cres.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cres.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Cres.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cres.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cres.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pLA for Csal prior - quadratic

##### Set data
data <- data.pLA.Csal.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Csal.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Csal.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Csal.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Csal.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Csal.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Csal.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Csal.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Avex prior - quadratic

##### Set data
data <- data.pLA.Avex.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Avex.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Avex.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Avex.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Avex.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Atri prior - quadratic

##### Set data
data <- data.pLA.Atri.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Atri.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Atri.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Atri.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Atri.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Atri.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Atri.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Atri.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pLA for Asol prior - quadratic

##### Set data
data <- data.pLA.Asol.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Asol.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Asol.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Asol.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Asol.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Asol.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Asol.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Asol.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pLA for Anig prior - quadratic

##### Set data
data <- data.pLA.Anig.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Anig.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Anig.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Anig.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Anig.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Anig.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Anig.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Anig.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### pLA for Anig prior - quadratic

##### Set data
data <- data.pLA.Anig.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Anig.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Anig.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Anig.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Anig.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Anig.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Anig.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Anig.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### pLA for Cmel prior - quadratic

##### Set data
data <- data.pLA.Cmel.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
pLA.Cmel.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cmel.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cmel.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cmel.prior, ylab = "pLA for all species", xlab = "Temperature")
lines(pLA.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 6. Fit gamma distributions to pLA prior thermal responses
##########

###################################### pLA for Cpip prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Cpip.prior.cf.dists <- data.frame(q = as.vector(pLA.Cpip.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(pLA.Cpip.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(pLA.Cpip.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Cpip.prior.gamma.fits = apply(pLA.Cpip.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Cqui prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Cqui.prior.cf.dists <- data.frame(q = as.vector(pLA.Cqui.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Cqui.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Cqui.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Cqui.prior.gamma.fits = apply(pLA.Cqui.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Ctar prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Ctar.prior.cf.dists <- data.frame(q = as.vector(pLA.Ctar.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Ctar.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Ctar.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Ctar.prior.gamma.fits = apply(pLA.Ctar.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Cmol prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Cmol.prior.cf.dists <- data.frame(q = as.vector(pLA.Cmol.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Cmol.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Cmol.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Cmol.prior.gamma.fits = apply(pLA.Cmol.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Cpal prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Cpal.prior.cf.dists <- data.frame(q = as.vector(pLA.Cpal.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Cpal.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Cpal.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Cpal.prior.gamma.fits = apply(pLA.Cpal.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Cres prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Cres.prior.cf.dists <- data.frame(q = as.vector(pLA.Cres.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Cres.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Cres.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Cres.prior.gamma.fits = apply(pLA.Cres.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Csal prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Csal.prior.cf.dists <- data.frame(q = as.vector(pLA.Csal.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Csal.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Csal.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Csal.prior.gamma.fits = apply(pLA.Csal.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


###################################### pLA for Avex prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Avex.prior.cf.dists <- data.frame(q = as.vector(pLA.Avex.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Avex.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Avex.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Avex.prior.gamma.fits = apply(pLA.Avex.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Atri prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Atri.prior.cf.dists <- data.frame(q = as.vector(pLA.Atri.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Atri.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Atri.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Atri.prior.gamma.fits = apply(pLA.Atri.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Asol prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Asol.prior.cf.dists <- data.frame(q = as.vector(pLA.Asol.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Asol.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Asol.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Asol.prior.gamma.fits = apply(pLA.Asol.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Anig prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Anig.prior.cf.dists <- data.frame(q = as.vector(pLA.Anig.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Anig.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Anig.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Anig.prior.gamma.fits = apply(pLA.Anig.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### pLA for Cmel prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
pLA.Cmel.prior.cf.dists <- data.frame(q = as.vector(pLA.Cmel.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(pLA.Cmel.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(pLA.Cmel.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
pLA.Cmel.prior.gamma.fits = apply(pLA.Cmel.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


pLA.hypers <- list(pLA.Cpip.prior.gamma.fits, pLA.Cqui.prior.gamma.fits, pLA.Ctar.prior.gamma.fits, pLA.Cmol.prior.gamma.fits, 
                   pLA.Cpal.prior.gamma.fits, pLA.Cres.prior.gamma.fits, pLA.Csal.prior.gamma.fits, pLA.Avex.prior.gamma.fits, 
                   pLA.Atri.prior.gamma.fits, pLA.Asol.prior.gamma.fits, pLA.Anig.prior.gamma.fits, pLA.Cmel.prior.gamma.fits)
save(pLA.hypers, file = "pLAhypers.Rsave")


##########
###### 7. Fit pLA thermal responses with data-informed priors
##########

load("pLAhypers.Rsave")
pLA.Cpip.prior.gamma.fits <- pLA.hypers[[1]]
pLA.Cqui.prior.gamma.fits <- pLA.hypers[[2]]
pLA.Ctar.prior.gamma.fits <- pLA.hypers[[3]]
pLA.Cmol.prior.gamma.fits <- pLA.hypers[[4]]
pLA.Cpal.prior.gamma.fits <- pLA.hypers[[5]]
pLA.Cres.prior.gamma.fits <- pLA.hypers[[6]]
pLA.Csal.prior.gamma.fits <- pLA.hypers[[7]]
pLA.Avex.prior.gamma.fits <- pLA.hypers[[8]]
pLA.Atri.prior.gamma.fits <- pLA.hypers[[9]]
pLA.Asol.prior.gamma.fits <- pLA.hypers[[10]]
pLA.Anig.prior.gamma.fits <- pLA.hypers[[11]]
pLA.Cmel.prior.gamma.fits <- pLA.hypers[[12]]

############## pLA for Cx. pipiens - quadratic

##### Set data 
data <- data.pLA.Cpip
hypers <- pLA.Cpip.prior.gamma.fits * .1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pLA.Cpip.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cpip.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cpip.out.inf)

save(pLA.Cpip.out.inf, file = "jagsout_pLA_Cpip_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cpip, ylab = "pLA for Cx. pipiens", xlab = "Temperature", pch = 1)
lines(pLA.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

# Get optimum for pLA: 23.1 C
Temp.xs[which.max(as.vector(pLA.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## pLA for Cx. quinquefasciatus - quadratic

##### Set data 
data <- data.pLA.Cqui
hypers <- pLA.Cqui.prior.gamma.fits * .1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pLA.Cqui.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cqui.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cqui.out.inf)

save(pLA.Cqui.out.inf, file = "jagsout_pLA_Cqui_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cqui, ylab = "pLA for Cx. quinquefasciatus", xlab = "Temperature", pch = 1)
lines(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

# Get optimum for pLA: 23.3 C
Temp.xs[which.max(as.vector(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## pLA for Cx. tarsalis - quadratic

##### Set data 
data <- data.pLA.Ctar
hypers <- pLA.Ctar.prior.gamma.fits * .025

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pLA.Ctar.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Ctar.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Ctar.out.inf)

save(pLA.Ctar.out.inf, file = "jagsout_pLA_Ctar_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Ctar, ylab = "pLA for Cx. tarsalis", xlab = "Temperature", pch = 1)
lines(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)

# Get optimum for pLA: 24.6 C
Temp.xs[which.max(as.vector(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## pLA for Ae. vexans - quadratic

##### Set data 
data <- data.pLA.Avex
hypers <- pLA.Avex.prior.gamma.fits * .05

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pLA.Avex.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Avex.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Avex.out.inf)

save(pLA.Avex.out.inf, file = "jagsout_pLA_Avex_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Avex, ylab = "pLA for Ae. vexans", xlab = "Temperature", pch = 1)
lines(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)

# Get optimum for pLA: 25.0 C
Temp.xs[which.max(as.vector(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## pLA for Ae. triseriatus - quadratic

##### Set data 
data <- data.pLA.Atri
hypers <- pLA.Atri.prior.gamma.fits * .05

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pLA.Atri.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Atri.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Atri.out.inf)

save(pLA.Atri.out.inf, file = "jagsout_pLA_Atri_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Atri, ylab = "pLA for Ae. triseriatus", xlab = "Temperature", pch = 1)
lines(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for pLA: 21.9 C
Temp.xs[which.max(as.vector(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## pLA for Cs. melanura - quadratic

##### Set data 
data <- data.pLA.Cmel
hypers <- pLA.Cmel.prior.gamma.fits * .05

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pLA.Cmel.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cmel.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cmel.out.inf)

save(pLA.Cmel.out.inf, file = "jagsout_pLA_Cmel_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cmel, ylab = "pLA for Cs. melanura", xlab = "Temperature", pch = 1)
lines(pLA.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for pLA: 23.5 C
Temp.xs[which.max(as.vector(pLA.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


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

Cpip.pLA.avg <- CalcTraitAvg(data.pLA.Cpip)
Cqui.pLA.avg <- CalcTraitAvg(data.pLA.Cqui)
Ctar.pLA.avg <- CalcTraitAvg(data.pLA.Ctar)
Atri.pLA.avg <- CalcTraitAvg(data.pLA.Atri)


##########
###### 9. Plot preliminary figures
##########

##################### Figure - all 3 Cx. spp. together

# Make polygon coordinates
Cpip.cord.y <- c(pLA.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Cqui.cord.y <- c(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
Ctar.cord.y <- c(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"], rev(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]))
cord.x <- c(Temp.xs, rev(Temp.xs))

# Plot pLA for the 3 Culex spp. 
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1.05), data = data.pLA, xaxt = "n", pch = 19,
     ylab = "Survival Probability", xlab = "", cex.lab = 1.15, col = "white")
polygon(cord.x, Cqui.cord.y, col = rgb(0.9, 0.25, 0.25, 0.2), border = NA)
polygon(cord.x, Ctar.cord.y, col = rgb(0.25, 0.25, 0.9, 0.2), border = NA)
lines(pLA.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 1, col = "darkgrey", lwd = 1.5)
lines(pLA.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 1, col = "darkgrey", lwd = 1.5)
arrows(Cpip.pLA.avg$Temp, Cpip.pLA.avg$upper, Cpip.pLA.avg$Temp, Cpip.pLA.avg$lower, length = 0, lwd = 1.5, col = "grey40")
arrows(Cqui.pLA.avg$Temp, Cqui.pLA.avg$upper, Cqui.pLA.avg$Temp, Cqui.pLA.avg$lower, length = 0, lwd = 1.5, col = "red")
arrows(Ctar.pLA.avg$Temp, Ctar.pLA.avg$upper, Ctar.pLA.avg$Temp, Ctar.pLA.avg$lower, length = 0, lwd = 1.5, col = "blue")
points(Cpip.pLA.avg$mean ~ Cpip.pLA.avg$Temp, pch = 19, col = "grey40")
points(Cqui.pLA.avg$mean ~ Cqui.pLA.avg$Temp, pch = 19, col = "red")
points(Ctar.pLA.avg$mean ~ Ctar.pLA.avg$Temp, pch = 19, col = "blue")
legend(x = 40, y = 1, bty = "n", legend = c("pip", "qui", "tar"), pch = 19, col = c("grey40", "red", "blue"))
legend(x = 38, y = 1, bty = "n", legend = c("", "", ""), fill = c("white", rgb(0.9, 0.25, 0.25, 0.3), rgb(0.25, 0.25, 0.9, 0.3)), border = c("darkgrey", "white", "white"))


##################### Figure - each spp. separate

par(mfrow = c(5,2), mar = c(1, 4.5, 1.75, 1), oma = c(3.5, 3, 1.75, 0))
par(mfrow = c(3,2), mar = c(1, 4.5, 1.75, 1), oma = c(3.5, 3, 1.75, 0))

##### pLA for Cx. pipiens pipiens
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cpip, xaxt = "n", pch = 19,
     ylab = "Survival Probability", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(pLA.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Cx. pipiens", side = 2, line = 5, font = 4)
mtext(text = "Larval-to-adult survival (pLA)", side = 3, line = 1.25, font = 2, cex = 0.95)

##### MDR for Cx. pipiens pipiens
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Cpippip, xaxt = "n", pch = 19,
     ylab = expression(paste("Development rate (day"^-1,")")), xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(MDR.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Cpippip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "B", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Mosquito development rate (MDR)", side = 3, line = 1.25, font = 2, cex = 0.95)

##### pLA for Cx. quinquefasciatus
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cqui, xaxt = "n", pch = 19,
     ylab = "Survival Probability", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "C", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Cx. quinquefasciatus", side = 2, line = 5, font = 4)

##### MDR for Cx. quinquefasciatus
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Cqui, xaxt = "n", pch = 19,
     ylab = expression(paste("Development rate (day"^-1,")")), xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "D", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### pLA for Cx. tarsalis
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Ctar, xaxt = "n", pch = 19,
     ylab = "Survival Probability", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "E", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Cx. tarsalis", side = 2, line = 5, font = 4)
mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 1)

##### MDR for Cx. tarsalis
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Ctar, xaxt = "n", pch = 19,
     ylab = expression(paste("Development rate (day"^-1,")")), xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "F", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 1)






##### pLA for Ae. triseriatus
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Atri, xaxt = "n", pch = 19,
     ylab = "Survival Probability", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "G", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Ae. triseriatus", side = 2, line = 5, font = 4)

##### MDR for Ae. triseriatus
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Atri, xaxt = "n", pch = 19,
     ylab = expression(paste("Development rate (day"^-1,")")), xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "H", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### pLA for Ae. vexans
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Avex, xaxt = "n", pch = 19,
     ylab = "Survival Probability", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "I", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 1)
mtext(text = "Ae. vexans", side = 2, line = 5, font = 4)

##### MDR for Ae. vexans
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Avex, xaxt = "n", pch = 19,
     ylab = expression(paste("Development rate (day"^-1,")")), xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "J", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 1)