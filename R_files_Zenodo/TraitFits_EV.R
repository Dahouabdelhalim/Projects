## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for egg viability for many Culex and Aedes species
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit EV thermal responses with uniform priors
##           5) Fit EV thermal responses for priors
##           6) Fit gamma distributions to EV prior thermal responses
##           7) Fit EV thermal responses with data-informed priors
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
data.EV <- read.csv("TraitData_EV.csv") # Data from database for most traits (except below)
unique(data.EV$host.code)

# Subset Data
data.EV.Cpip <- subset(data.EV, host.code == "Cpip")
data.EV.Cqui <- subset(data.EV, host.code == "Cqui")
data.EV.Cthe <- subset(data.EV, host.code == "Cthe")
data.EV.Cmol <- subset(data.EV, host.code == "Cmol")
data.EV.Avex <- subset(data.EV, host.code == "Avex")
data.EV.Ador <- subset(data.EV, host.code == "Ador")
data.EV.Anig <- subset(data.EV, host.code == "Anig")

par(mfrow = c(1,1), mar = c(4.5, 4.5, 1, 1))
plot(trait ~ T, xlim = c(5, 45), data = data.EV, ylab = "EV for Cx. mosquitoes", xlab = "Temperature")
points(trait ~ T, data = data.EV.Cpip, col = "grey")
points(trait ~ T, data = data.EV.Cqui, col = "red")
points(trait ~ T, data = data.EV.Cmol, col = "violet")
points(trait ~ T, data = data.EV.Cthe, col = "purple")
points(trait ~ T, data = data.EV.Avex, col = "cyan3")
points(trait ~ T, data = data.EV.Ador, col = "green")
points(trait ~ T, data = data.EV.Anig, col = "seagreen")


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
###### 4. Fit EV thermal responses with uniform priors
##########

###################################### EV for Cx pipiens - quadratic

##### Set data
data <- data.EV.Cpip

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EV.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Cpip.out$BUGSoutput$summary[1:5,]
mcmcplot(EV.Cpip.out)

save(EV.Cpip.out, file = "jagsout_EV_Cpip.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cpip, ylab = "EV for Cx pipiens", xlab = "Temperature")
lines(EV.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### EV for Cx quinquefasciatus - quadratic

##### Set data
data <- data.EV.Cqui

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EV.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Cqui.out$BUGSoutput$summary[1:5,]
mcmcplot(EV.Cqui.out)

save(EV.Cqui.out, file = "jagsout_EV_Cqui.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cqui, ylab = "EV for Cx pipiens", xlab = "Temperature")
lines(EV.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### EV for Ae. vexans - quadratic

##### Set data
data <- data.EV.Avex

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EV.Avex.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Avex.out$BUGSoutput$summary[1:5,]
mcmcplot(EV.Avex.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Avex, ylab = "EV for Cx pipiens", xlab = "Temperature")
lines(EV.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 5. Fit thermal responses for EV priors
##########

# Subset Data - leave-one-out
data.EV.Cpip.prior <- subset(data.EV, host.code != "Cpip")
data.EV.Cqui.prior <- subset(data.EV, host.code != "Cqui")
data.EV.Cthe.prior <- subset(data.EV, host.code != "Cthe")
data.EV.Avex.prior <- subset(data.EV, host.code != "Avex")

###################################### EV for Cx pipiens prior - quadratic

##### Set data
data <- data.EV.Cpip.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EV.Cpip.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Cpip.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(EV.Cpip.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cpip.prior, ylab = "EV for Cx pipiens prior", xlab = "Temperature")
lines(EV.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### EV for Cx quinquefasciatus prior - briere

##### Set data
data <- data.EV.Cqui.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EV.Cqui.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Cqui.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(EV.Cqui.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cqui.prior, ylab = "EV for Cx quinquefasciatus prior", xlab = "Temperature")
lines(EV.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### EV for Cx theileri prior - quadratic

##### Set data
data <- data.EV.Cthe.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EV.Cthe.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Cthe.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(EV.Cthe.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cthe.prior, ylab = "EV for Cx theileri prior", xlab = "Temperature")
lines(EV.Cthe.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cthe.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cthe.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### EV for Ae. vexans prior - quadratic

##### Set data
data <- data.EV.Avex.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
EV.Avex.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Avex.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(EV.Avex.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Avex.prior, ylab = "EV for Ae vexans prior", xlab = "Temperature")
lines(EV.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 6. Fit gamma distributions to EV prior thermal responses
##########

###################################### EV for Cpip prior

# Get the posterior dists for all 4 parameters into a data frame
EV.Cpip.prior.cf.dists <- data.frame(q = as.vector(EV.Cpip.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(EV.Cpip.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(EV.Cpip.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EV.Cpip.prior.gamma.fits = apply(EV.Cpip.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### EV for Cqui prior

# Get the posterior dists for all 4 parameters into a data frame
EV.Cqui.prior.cf.dists <- data.frame(q = as.vector(EV.Cqui.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EV.Cqui.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EV.Cqui.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EV.Cqui.prior.gamma.fits = apply(EV.Cqui.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### EV for Cthe prior

# Get the posterior dists for all 4 parameters into a data frame
EV.Cthe.prior.cf.dists <- data.frame(q = as.vector(EV.Cthe.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EV.Cthe.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EV.Cthe.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EV.Cthe.prior.gamma.fits = apply(EV.Cthe.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


###################################### EV for Avex prior

# Get the posterior dists for all 4 parameters into a data frame
EV.Avex.prior.cf.dists <- data.frame(q = as.vector(EV.Avex.prior.out$BUGSoutput$sims.list$cf.q),
                                     T0 = as.vector(EV.Avex.prior.out$BUGSoutput$sims.list$cf.T0), 
                                     Tm = as.vector(EV.Avex.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
EV.Avex.prior.gamma.fits = apply(EV.Avex.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

EV.hypers <- list(EV.Cpip.prior.gamma.fits, EV.Cqui.prior.gamma.fits, EV.Cthe.prior.gamma.fits, EV.Avex.prior.gamma.fits)
save(EV.hypers, file = "EVhypers.Rsave")


##########
###### 7. Fit EV thermal responses with data-informed priors
##########

load("EVhypers.Rsave")
EV.Cpip.prior.gamma.fits <- EV.hypers[[1]]
EV.Cqui.prior.gamma.fits <- EV.hypers[[2]]
EV.Cthe.prior.gamma.fits <- EV.hypers[[3]]
EV.Avex.prior.gamma.fits <- EV.hypers[[4]]

############## EV for Cx. pipiens - quadratic

##### Set data 
data <- data.EV.Cpip
hypers <- EV.Cpip.prior.gamma.fits * 0.2

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
EV.Cpip.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Cpip.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(EV.Cpip.out.inf)

save(EV.Cpip.out.inf, file = "jagsout_EV_Cpip_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cpip, ylab = "EV for Cx. pipiens ", xlab = "Temperature", pch = 1)
lines(EV.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.2)

# Get optimum for EV: 23.2 C
Temp.xs[which.max(as.vector(EV.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## EV for Cx. quinquefasciatus - briere

##### Set data 
data <- data.EV.Cqui
data <- data.EV.Cqui[4:22,] # Remove Oda data - it shows no temp-dep, and briere fit doesn't work with it (T0 = 0)
hypers <- EV.Cqui.prior.gamma.fits * 0.1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
EV.Cqui.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Cqui.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(EV.Cqui.out.inf)

save(EV.Cqui.out.inf, file = "jagsout_EV_Cqui_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cqui[4:22,], ylab = "EV for Cx. quinquefasciatus ", xlab = "Temperature", pch = 1)
lines(EV.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "G", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for EV: 32 C
Temp.xs[which.max(as.vector(EV.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


############## EV for Cx. theileri - quadratic

##### Set data 
data <- data.EV.Cthe
hypers <- EV.Cthe.prior.gamma.fits * 1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
EV.Cthe.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Cthe.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(EV.Cthe.out.inf)

save(EV.Cthe.out.inf, file = "jagsout_EV_Cthe_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Cthe, ylab = "EV for Cx. theileri ", xlab = "Temperature", pch = 1)
lines(EV.Cthe.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cthe.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Cthe.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "H", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for EV: 23.7 C
Temp.xs[which.max(as.vector(EV.Cthe.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## EV for Ae. vexans - quadratic

##### Set data 
data <- data.EV.Avex
hypers <- EV.Avex.prior.gamma.fits * .01

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
EV.Avex.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EV.Avex.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(EV.Avex.out.inf)

save(EV.Avex.out.inf, file = "jagsout_EV_Avex_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.EV.Avex, ylab = "EV for Ae. vexans ", xlab = "Temperature", pch = 1)
lines(EV.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EV.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EV.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "I", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for EV: 28 C
Temp.xs[which.max(as.vector(EV.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


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

Cqui.EV.avg <- CalcTraitAvg(data.EV.Cqui)
Cthe.EV.avg <- CalcTraitAvg(data.EV.Cthe)
Avex.EV.avg <- CalcTraitAvg(data.EV.Avex)