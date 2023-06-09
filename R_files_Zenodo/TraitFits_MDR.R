## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for mosquito development rate (MDR) for 
##          many Culex and Aedes species - 1) with uniform priors and 2) fitting priors using a leave-one-out approach
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit MDR thermal responses with uniform priors
##           5) Fit MDR thermal responses for priors
##           6) Fit gamma distributions to MDR prior thermal responses
##           7) Fit MDR thermal responses with data-informed priors
##           8) Calculate Treatment averages for plotting


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
data.MDR <- read.csv("TraitData_MDR.csv") # Data from database for most traits (except below)
unique(data.MDR$host.code)

# Subset Data
data.MDR.Cpip <- subset(data.MDR, host.code == "Cpip")
data.MDR.Cqui <- subset(data.MDR, host.code == "Cqui")
data.MDR.Ctar <- subset(data.MDR, host.code == "Ctar")
data.MDR.Cmol <- subset(data.MDR, host.code == "Cmol")
data.MDR.Cres <- subset(data.MDR, host.code == "Cres")
data.MDR.Csal <- subset(data.MDR, host.code == "Csal")
data.MDR.Cpal <- subset(data.MDR, host.code == "Cpal")
data.MDR.Atri <- subset(data.MDR, host.code == "Atri")
data.MDR.Avex <- subset(data.MDR, host.code == "Avex")
data.MDR.Asol <- subset(data.MDR, host.code == "Asol")
data.MDR.Anig <- subset(data.MDR, host.code == "Anig")
data.MDR.Cmel <- subset(data.MDR, host.code == "Cmel") # This data was added later, not included in priors for other species

par(mfrow = c(1,1))
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR, ylab = "MDR for Cx. mosquitoes", xlab = "Temperature")
points(1/trait ~ T, data = data.MDR.Cpip, col = "grey")
points(1/trait ~ T, data = data.MDR.Cqui, col = "red")
points(1/trait ~ T, data = data.MDR.Ctar, col = "blue")
points(1/trait ~ T, data = data.MDR.Cmol, col = "violet")
points(1/trait ~ T, data = data.MDR.Cpal, col = "darkorange")
points(1/trait ~ T, data = data.MDR.Cres, col = "violetred")
points(1/trait ~ T, data = data.MDR.Csal, col = "goldenrod2")
points(1/trait ~ T, data = data.MDR.Atri, col = "darkgreen")
points(1/trait ~ T, data = data.MDR.Avex, col = "cyan3")
points(1/trait ~ T, data = data.MDR.Asol, col = "green")
points(1/trait ~ T, data = data.MDR.Anig, col = "seagreen")


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
###### 4. Fit MDR thermal responses with uniform priors
##########

###################################### MDR for Cx pipiens - Briere

##### Set data
data <- data.MDR.Cpip

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cpip.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cpip.out)

save(MDR.Cpip.out, file = "jagsout_MDR_Cpip.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Cpip, ylab = "MDR for Cx pipiens pipiens", xlab = "Temperature")
lines(MDR.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cpip.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### MDR for Cx quinquefasciatus - Briere

##### Set data
data <- data.MDR.Cqui

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cqui.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cqui.out)

save(MDR.Cqui.out, file = "jagsout_MDR_Cqui.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Cqui, ylab = "MDR for Cx quinquefasciatus", xlab = "Temperature")
lines(MDR.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cqui.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### MDR for Cx tarsalis - Briere

##### Set data
data <- data.MDR.Ctar

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Ctar.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Ctar.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Ctar.out)

save(MDR.Ctar.out, file = "jagsout_MDR_Ctar.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Ctar, ylab = "MDR for Cx tarsalis", xlab = "Temperature")
lines(MDR.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Ctar.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### MDR for Ae. triseriatus - Briere

##### Set data
data <- data.MDR.Atri

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Atri.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Atri.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Atri.out)

save(MDR.Atri.out, file = "jagsout_MDR_Atri.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Atri, ylab = "MDR for Ae triseriatus", xlab = "Temperature")
lines(MDR.Atri.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Atri.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Atri.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### MDR for Ae. vexans - Briere

##### Set data
data <- data.MDR.Avex

##### Organize Data for JAGS
trait <- data$trait # Avex is the only data in regular MDR format (so not 1/data$trait)
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Avex.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Avex.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Avex.out)

save(MDR.Avex.out, file = "jagsout_MDR_Avex.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Avex, ylab = "MDR for Ae vexans", xlab = "Temperature")
lines(MDR.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

Temp.xs[which.max(as.vector(MDR.Avex.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


##########
###### 5. Fit thermal responses for MDR priors using a leave-one-out approach
##########

# Subset Data - leave-one-out
data.MDR.Cpip.prior <- subset(data.MDR, host.code != "Cpip")
data.MDR.Cqui.prior <- subset(data.MDR, host.code != "Cqui")
data.MDR.Ctar.prior <- subset(data.MDR, host.code != "Ctar")
data.MDR.Atri.prior <- subset(data.MDR, host.code != "Atri")
data.MDR.Avex.prior <- subset(data.MDR, host.code != "Avex")
data.MDR.Cmel.prior <- subset(data.MDR, host.code != "Cmel")

###################################### MDR for Cpip prior - Briere

##### Set data
data <- data.MDR.Cpip.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Cpip.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cpip.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cpip.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Cpip.prior, ylab = "MDR for Cx pipiens prior", xlab = "Temperature")
lines(MDR.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cpip.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### MDR for Cqui prior - Briere

##### Set data
data <- data.MDR.Cqui.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Cqui.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cqui.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cqui.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Cqui.prior, ylab = "MDR for Cx. quinquefasciatus prior", xlab = "Temperature")
lines(MDR.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cqui.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### MDR for Ctar prior - Briere

##### Set data
data <- data.MDR.Ctar.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Ctar.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Ctar.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Ctar.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Ctar.prior, ylab = "MDR for Cx tarsalis prior", xlab = "Temperature")
lines(MDR.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Ctar.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### MDR for Avex prior - Briere

##### Set data
data <- data.MDR.Avex.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Avex.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Avex.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Avex.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Avex.prior, ylab = "MDR for Ae vexans prior", xlab = "Temperature")
lines(MDR.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Avex.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### MDR for Atri prior - Briere

##### Set data
data <- data.MDR.Atri.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Atri.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Atri.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Atri.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Atri.prior, ylab = "MDR for Ae. triseriatus prior", xlab = "Temperature")
lines(MDR.Atri.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Atri.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Atri.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### MDR for Cmel prior - Briere

##### Set data
data <- data.MDR.Cmel.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
MDR.Cmel.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                           n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cmel.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cmel.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.MDR.Cmel.prior, ylab = "MDR for Ae. triseriatus prior", xlab = "Temperature")
lines(MDR.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cmel.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 6. Fit gamma distributions to MDR prior thermal responses
##########

###################################### MDR for Cpip prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.Cpip.prior.cf.dists <- data.frame(q = as.vector(MDR.Cpip.prior.out$BUGSoutput$sims.list$cf.q),
                                T0 = as.vector(MDR.Cpip.prior.out$BUGSoutput$sims.list$cf.T0), 
                                Tm = as.vector(MDR.Cpip.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.Cpip.prior.gamma.fits = apply(MDR.Cpip.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### MDR for Cqui prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.Cqui.prior.cf.dists <- data.frame(q = as.vector(MDR.Cqui.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(MDR.Cqui.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(MDR.Cqui.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.Cqui.prior.gamma.fits = apply(MDR.Cqui.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### MDR for Ctar prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.Ctar.prior.cf.dists <- data.frame(q = as.vector(MDR.Ctar.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(MDR.Ctar.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(MDR.Ctar.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.Ctar.prior.gamma.fits = apply(MDR.Ctar.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### MDR for Avex prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.Avex.prior.cf.dists <- data.frame(q = as.vector(MDR.Avex.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(MDR.Avex.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(MDR.Avex.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.Avex.prior.gamma.fits = apply(MDR.Avex.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### MDR for Avex prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.Avex.prior.cf.dists <- data.frame(q = as.vector(MDR.Avex.prior.out$BUGSoutput$sims.list$cf.q),
                                      T0 = as.vector(MDR.Avex.prior.out$BUGSoutput$sims.list$cf.T0), 
                                      Tm = as.vector(MDR.Avex.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.Avex.prior.gamma.fits = apply(MDR.Avex.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)


###################################### MDR for Cmel prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
MDR.Cmel.prior.cf.dists <- data.frame(q = as.vector(MDR.Cmel.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(MDR.Cmel.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(MDR.Cmel.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
MDR.Cmel.prior.gamma.fits = apply(MDR.Cmel.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

MDR.hypers <- list(MDR.Cpip.prior.gamma.fits, MDR.Cqui.prior.gamma.fits, MDR.Ctar.prior.gamma.fits,
                   MDR.Avex.prior.gamma.fits, MDR.Atri.prior.gamma.fits, MDR.Cmel.prior.gamma.fits)
save(MDR.hypers, file = "MDRhypers.Rsave")


##########
###### 7. Fit MDR thermal responses with data-informed priors
##########

load("MDRhypers.Rsave")
MDR.Cpip.prior.gamma.fits <- MDR.hypers[[1]]
MDR.Cqui.prior.gamma.fits <- MDR.hypers[[2]]
MDR.Ctar.prior.gamma.fits <- MDR.hypers[[3]]
MDR.Avex.prior.gamma.fits <- MDR.hypers[[4]]
MDR.Atri.prior.gamma.fits <- MDR.hypers[[5]]
MDR.Cmel.prior.gamma.fits <- MDR.hypers[[6]]

############## MDR for Cx. pipiens - Briere

##### Set data 
data <- data.MDR.Cpip
hypers <- MDR.Cpip.prior.gamma.fits * 0.1

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
MDR.Cpip.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cpip.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cpip.out.inf)

save(MDR.Cpip.out.inf, file = "jagsout_MDR_Cpip_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Cpip, ylab = "MDR for Cx. pipiens", xlab = "Temperature", pch = 1)
lines(MDR.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

# Get optimum for MDR: 30.9 C
Temp.xs[which.max(as.vector(MDR.Cpip.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## MDR for Cx. quinquefasciatus - Briere

##### Set data 
data <- data.MDR.Cqui
hypers <- MDR.Cqui.prior.gamma.fits * 0.1

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
MDR.Cqui.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cqui.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cqui.out.inf)

save(MDR.Cqui.out.inf, file = "jagsout_MDR_Cqui_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Cqui, ylab = "MDR for Cx. quinquefasciatus", xlab = "Temperature", pch = 1)
lines(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

# Get optimum for MDR: 31 C
Temp.xs[which.max(as.vector(MDR.Cqui.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## MDR for Cx. tarsalis - Briere

##### Set data 
data <- data.MDR.Ctar
hypers <- MDR.Ctar.prior.gamma.fits * 0.1

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
MDR.Ctar.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Ctar.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Ctar.out.inf)

save(MDR.Ctar.out.inf, file = "jagsout_MDR_Ctar_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Ctar, ylab = "MDR for Cx. tarsalis", xlab = "Temperature", pch = 1)
lines(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)

# Get optimum for MDR: 32.3 C
Temp.xs[which.max(as.vector(MDR.Ctar.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## MDR for Ae. vexans - Briere

##### Set data 
data <- data.MDR.Avex
hypers <- MDR.Avex.prior.gamma.fits * 0.5

##### Organize Data for JAGS
trait <- 1/data$trait # Avex is the only data in regular MDR format (so not 1/data$trait)
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
MDR.Avex.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Avex.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Avex.out.inf)

save(MDR.Avex.out.inf, file = "jagsout_MDR_Avex_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Avex, ylab = "MDR for Ae. vexans", xlab = "Temperature", pch = 1)
lines(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)

# Get optimum for MDR: 30.8 C
Temp.xs[which.max(as.vector(MDR.Avex.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## MDR for Ae. triseriatus - Briere

##### Set data 
data <- data.MDR.Atri
hypers <- MDR.Atri.prior.gamma.fits * 0.2

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
MDR.Atri.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Atri.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Atri.out.inf)

save(MDR.Atri.out.inf, file = "jagsout_MDR_Atri_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Atri, ylab = "MDR for Ae. triseriatus", xlab = "Temperature", pch = 1)
lines(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for MDR: 29.0 C
Temp.xs[which.max(as.vector(MDR.Atri.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## MDR for Cs. melanura - Briere

##### Set data 
data <- data.MDR.Cmel
hypers <- MDR.Cmel.prior.gamma.fits * .1

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
MDR.Cmel.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cmel.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cmel.out.inf)

save(MDR.Cmel.out.inf, file = "jagsout_MDR_Cmel_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Cmel, ylab = "MDR for Cs. melanura", xlab = "Temperature", pch = 1)
lines(MDR.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for MDR: 30.8 C
Temp.xs[which.max(as.vector(MDR.Cmel.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


##########
###### 8. Calculate Treatment averages for plotting
##########

# Function to trait averages for plotting only (JAGS models fit to raw data)
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

Cpip.MDR.avg <- CalcTraitAvgRecip(data.MDR.Cpip)
Cqui.MDR.avg <- CalcTraitAvgRecip(data.MDR.Cqui)
Ctar.MDR.avg <- CalcTraitAvgRecip(data.MDR.Ctar)

#########
##### Figures
########

## See pLA (larval survival) R script for joint preliminary figures