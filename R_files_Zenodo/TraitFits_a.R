## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for biting rate (a) for many Culex and Aedes species
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit a thermal responses with uniform priors
##           5) Fit a thermal responses for priors
##           6) Fit gamma distributions to a prior thermal responses
##           7) Fit a thermal responses with data-informed priors
##           8) Calculate treatment averages for plotting
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

# Load Data - NOTE: these data are in GCD (gonotrophic cycle duration) format
data.a <- read.csv("TraitData_a.csv")

# Subset Data
data.a.Cpip <- subset(data.a, host.code == "Cpip") 
data.a.Cqui <- subset(data.a, host.code == "Cqui")
data.a.Ctar <- subset(data.a, host.code == "Ctar")
data.a.Cpal <- subset(data.a, host.code == "Cpal")
data.a.Cmel <- subset(data.a, host.code == "Cmel") # This data was added later, not included in priors for other species

# Plot Data
par(mfrow = c(1,1))
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a, ylab = "a for Cx. mosquitoes", xlab = "Temperature")
points(1/trait ~ T, data = data.a.Cpip, col = "grey")
points(1/trait ~ T, data = data.a.Cqui, col = "red")
points(1/trait ~ T, data = data.a.Ctar, col = "blue")
points(1/trait ~ T, data = data.a.Cpal, col = "yellow")


##########
###### 2. JAGS Models
##########

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

############## Briere Model with gamma priors (except sigma)

sink("briere_inf.txt")
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
    z.trait.mu.pred[i] <- cf.q * Temp.xs[i] * (Temp.xs[i] - cf.T0) * sqrt((cf.Tm - Temp.xs[i]) * (cf.Tm > Temp.xs[i])) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


############## Briere Model with gamma priors (except sigma and Tmax) - for Cx. tarsalis fit only

sink("briere_comp.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dunif(25, 45)
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
# For actual fits
Temp.xs <- seq(1, 45, 0.1)
N.Temp.xs <-length(Temp.xs)
# For priors - fewer temps for derived calculations makes it go faster
Temp.xs <- seq(5, 45, 0.5)
N.Temp.xs <-length(Temp.xs)


##########
###### 4. Fit a thermal responses with uniform priors
##########

###################################### a for Cx pipiens - Briere

##### Set data
data <- data.a.Cpip

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
a.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Cpip.out$BUGSoutput$summary[1:5,]
mcmcplot(a.Cpip.out)

save(a.Cpip.out, file = "jagsout_a_Cpip.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Cpip, ylab = "a for Cx pipiens", xlab = "Temperature")
lines(a.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### a for Cx quinquefasciatus - Briere

##### Set data
data <- data.a.Cqui

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
a.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Cqui.out$BUGSoutput$summary[1:5,]
mcmcplot(a.Cqui.out)

save(a.Cqui.out, file = "jagsout_a_Cqui.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.a.Cqui, ylab = "a for Cx quinquefasciatus", xlab = "Temperature")
lines(a.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### a for Cx tarsalis - Briere

##### Set data
data <- data.a.Ctar

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
a.Ctar.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Ctar.out$BUGSoutput$summary[1:5,]
mcmcplot(a.Ctar.out)

save(a.Ctar.out, file = "jagsout_a_Ctar.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Ctar, ylab = "a for Cx tarsalis", xlab = "Temperature")
lines(a.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 5. Fit a thermal responses for a priors
##########

# Create leave-one-out subsetted data
data.a.Cpip.prior <- subset(data.a, host.code != "Cpip")
data.a.Cqui.prior <- subset(data.a, host.code != "Cqui")
data.a.Ctar.prior <- subset(data.a, host.code != "Ctar")
data.a.Cmel.prior <- subset(data.a, host.code != "Cmel")

###################################### a for Cx. pipiens prior - Briere

##### Set data (and swith Avex data from a to 1/a to match the rest)
data <- data.a.Cpip.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
a.Cpip.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
a.Cpip.prior.out.T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Cpip.prior.out$BUGSoutput$summary[1:5,]
a.Cpip.prior.out.T$BUGSoutput$summary[1:5,]
mcmcplot(a.Cpip.prior.out)

# Plot data + fit
par(mfrow = c(2,1))
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.3), data = data.a.Cpip.prior, ylab = "a for Cx. pipiens prior", xlab = "Temperature")
lines(a.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.3), data = data.a.Cpip.prior, ylab = "a for Cx. pipiens prior", xlab = "Temperature")
lines(a.Cpip.prior.out.T$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cpip.prior.out.T$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cpip.prior.out.T$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### a for Cx. quinquefasciatus prior - Briere

##### Set data (and swith Avex data from a to 1/a to match the rest)
data <- data.a.Cqui.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
a.Cqui.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
a.Cqui.prior.out.T <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_T.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())


##### Examine Output
a.Cqui.prior.out$BUGSoutput$summary[1:5,]
a.Cqui.prior.out.T$BUGSoutput$summary[1:5,]
mcmcplot(a.Cqui.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.4), data = data.a.Cqui.prior, ylab = "a for Cx. quinquefasciatus prior", xlab = "Temperature")
lines(a.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.4), data = data.a.Cqui.prior, ylab = "a for Cx. quinquefasciatus prior", xlab = "Temperature")
lines(a.Cqui.prior.out.T$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cqui.prior.out.T$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cqui.prior.out.T$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### a for Cx. tarsalis prior - Briere

##### Set data (and swith Avex data from a to 1/a to match the rest)
data <- data.a.Ctar.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
a.Ctar.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Ctar.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(a.Ctar.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 50), ylim = c(0,0.5), data = data.a.Ctar.prior, ylab = "a for Cx. tarsalis prior", xlab = "Temperature")
lines(a.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### a for Cs. melanura prior - Briere

##### Set data (and swith Avex data from a to 1/a to match the rest)
data <- data.a.Cmel.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
a.Cmel.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Cmel.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(a.Cmel.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 50), ylim = c(0,0.5), data = data.a.Cmel.prior, ylab = "a for Cx. tarsalis prior", xlab = "Temperature")
lines(a.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 6. Fit gamma distributions to a prior thermal responses
##########

###################################### a for Cx pipiens prior

# Get the posterior dists for all 4 parameters into a data frame
a.Cpip.prior.cf.dists <- data.frame(q = as.vector(a.Cpip.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(a.Cpip.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(a.Cpip.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
a.Cpip.prior.gamma.fits = apply(a.Cpip.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### a for Cx quinquefasciatus prior

# Get the posterior dists for all 4 parameters into a data frame
a.Cqui.prior.cf.dists <- data.frame(q = as.vector(a.Cqui.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(a.Cqui.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(a.Cqui.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
a.Cqui.prior.gamma.fits = apply(a.Cqui.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### a for Cx tarsalis prior

# Get the posterior dists for all 4 parameters into a data frame
a.Ctar.prior.cf.dists <- data.frame(q = as.vector(a.Ctar.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(a.Ctar.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(a.Ctar.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
a.Ctar.prior.gamma.fits = apply(a.Ctar.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### a for Cs melanura prior

# Get the posterior dists for all 4 parameters into a data frame
a.Cmel.prior.cf.dists <- data.frame(q = as.vector(a.Cmel.prior.out$BUGSoutput$sims.list$cf.q),
                                    T0 = as.vector(a.Cmel.prior.out$BUGSoutput$sims.list$cf.T0), 
                                    Tm = as.vector(a.Cmel.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
a.Cmel.prior.gamma.fits = apply(a.Cmel.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

a.hypers <- list(a.Cpip.prior.gamma.fits, a.Cqui.prior.gamma.fits, a.Ctar.prior.gamma.fits, a.Cmel.prior.gamma.fits)
save(a.hypers, file = "ahypers.Rsave")


##########
###### 7. Fit a thermal responses with data-informed priors
##########

load("ahypers.Rsave")
a.Cpip.prior.gamma.fits <- a.hypers[[1]]
a.Cqui.prior.gamma.fits <- a.hypers[[2]]
a.Ctar.prior.gamma.fits <- a.hypers[[3]]
a.Cmel.prior.gamma.fits <- a.hypers[[4]]

############## a for Cx. pipiens - Briere

##### Set data 
data <- data.a.Cpip
hypers <- a.Cpip.prior.gamma.fits * 0.5

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
a.Cpip.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Cpip.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(a.Cpip.out.inf)

save(a.Cpip.out.inf, file = "jagsout_a_Cpip_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Cpip, ylab = "a for Cx. pipiens", xlab = "Temperature", pch = 1)
lines(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

# Get optimum for a: 32.7  C
Temp.xs[which.max(as.vector(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## a for Cx. quinquefasciatus - Briere

##### Set data 
data <- data.a.Cqui
hypers <- a.Cqui.prior.gamma.fits * 0.1
hypers[,3] <- a.Cqui.prior.gamma.fits[,3]

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
a.Cqui.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Cqui.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(a.Cqui.out.inf)

save(a.Cqui.out.inf, file = "jagsout_a_Cqui_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Cqui, ylab = "a for Cx. quinquefasciatus", xlab = "Temperature", pch = 1)
lines(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

# Get optimum for a: 31.9  C
Temp.xs[which.max(as.vector(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## a for Cx. tarsalis - Briere

##### Set data 
data <- data.a.Ctar
hypers <- a.Ctar.prior.gamma.fits * 0.05
hypers[,3] <- a.Cqui.prior.gamma.fits[,3]

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
a.Ctar.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_comp.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Ctar.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(a.Ctar.out.inf)

save(a.Ctar.out.inf, file = "jagsout_a_Ctar_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.51), data = data.a.Ctar, ylab = "a for Cx. tarsalis", xlab = "Temperature", pch = 1)
lines(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for a: 25.9 C
Temp.xs[which.max(as.vector(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## a for Cs. melanura - Briere

##### Set data 
data <- data.a.Cmel
hypers <- a.Ctar.prior.gamma.fits * .75
hypers[,3] <- a.Ctar.prior.gamma.fits[,3] * .1

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
a.Cmel.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                       n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
a.Cmel.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(a.Cmel.out.inf)

save(a.Cmel.out.inf, file = "jagsout_a_Cmel_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Cmel, ylab = "a for Cs. melanura", xlab = "Temperature", pch = 1)
lines(a.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(a.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(a.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for a: 26.3  C
Temp.xs[which.max(as.vector(a.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


##########
###### 8. Calculate Treatment averages for plotting
##########

# Function to calc trait averages for plotting only (JAGS models fit to raw data)
CalcTraitAvgRecip = function(data.set){
  data.set <- data.set
  temp.list <- unique(data.set$T)
  out <- data.frame(Temp = numeric(length(temp.list)), mean = numeric(length(temp.list)), sd = numeric(length(temp.list)), SE = numeric(length(temp.list)), 
                    n = numeric(length(temp.list)), upper = numeric(length(temp.list)), lower = numeric(length(temp.list)))
  for(i in 1:length(temp.list)){
    data.sub <- subset(data.set, T == temp.list[i])
    out$Temp[i] <- temp.list[i]
    out$mean[i] <- mean(1/data.sub$trait)
    out$sd[i] <- sd(1/data.sub$trait)
    out$n[i] <- nrow(data.sub)
    out$SE[i] <- sd(1/data.sub$trait) / sqrt(nrow(data.sub))
    out$lower[i] <- mean(1/data.sub$trait) - sd(1/data.sub$trait) / sqrt(nrow(data.sub))
    out$upper[i] <- mean(1/data.sub$trait) + sd(1/data.sub$trait) / sqrt(nrow(data.sub))
  }
  out
}

Cpip.a.avg <- CalcTraitAvgRecip(data.a.Cpip)
Cqui.a.avg <- CalcTraitAvgRecip(data.a.Cqui)


##########
###### 9. Plot Preliminary Figures
##########

par(mfrow = c(3,3), mar = c(1, 4.5, 1.75, 1), oma = c(3.5, 3, 1.75, 0))

##### a for Cx. pipiens
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Cpip, xaxt = "n", pch = 19,
     ylab = "Biting rate (1/day)", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(a.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Cx. pipiens", side = 2, line = 5, font = 4)
mtext(text = "Biting Rate (a)", side = 3, line = 1.25, font = 2, cex = 0.95)

##### EFD for Cx. pipiens
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,360), data = data.EFD.Cpip, xaxt = "n", pch = 19,
     ylab = "Eggs per female per day", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(EFD.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(EFD.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(EFD.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "B", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Fecundity (EFD)", side = 3, line = 1.25, font = 2, cex = 0.95)

##### lf for Cx. pipiens
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Cpip, xaxt = "n", pch = 19,
     ylab = "Lifespan (days)", xlab = "", cex.lab = 1.15, col = "white")
legend("topleft", legend = "C", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Lifespan (lf)", side = 3, line = 1.25, font = 2, cex = 0.95)

##### a for Cx. quinquefasciatus
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Cqui, xaxt = "n", pch = 19,
     ylab = "Biting rate (1/day)", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(a.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "D", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Cx. pipiens", side = 2, line = 5, font = 4)

##### EFD for Cx. quinquefasciatus
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,360), data = data.EFD.Cqui, xaxt = "n", pch = 19,
     ylab = "Eggs per female per day", xlab = "", cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(EFD.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(EFD.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(EFD.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "E", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### lf for Cx. quinquefasciatus
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Cpip, xaxt = "n", pch = 19,
     ylab = "Lifespan (days)", xlab = "", cex.lab = 1.15, col = "white")
legend("topleft", legend = "F", bty= "n", cex = 1.6, adj = c(1.5, 0))
mtext(text = "Lifespan (lf)", side = 3, line = 1.25, font = 2, cex = 0.95)

##### a for Cx. tarsalis
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.a.Ctar, xaxt = "n", pch = 19,
     ylab = "Biting rate (1/day)", xlab = "", main = expression(paste(italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(a.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "G", bty= "n", cex = 1.6, adj = c(1.5, 0))