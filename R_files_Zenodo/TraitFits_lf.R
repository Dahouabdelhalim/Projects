## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for lifespan for Culex species with linear thermal responses
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit lf thermal responses with uniform priors
##           5) Fit lf thermal responses for priors
##           6) Fit gamma distributions to lf prior thermal responses
##           7) Fit lf thermal responses with data-informed priors
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

# Load Data
data.lf <- read.csv("TraitData_lf.csv") # Data from database for most traits (except below)
unique(data.lf$host.code)

# Subset Data
data.lf.Cpip <- subset(data.lf, host.code == "Cpip")
data.lf.Cqui <- subset(data.lf, host.code == "Cqui")
data.lf.Ctar <- subset(data.lf, host.code == "Ctar")
data.lf.Cmol <- subset(data.lf, host.code == "Cmol")
data.lf.Cpal <- subset(data.lf, host.code == "Cpal")
data.lf.Cres <- subset(data.lf, host.code == "Cres")
data.lf.Atae <- subset(data.lf, host.code == "Atae")

par(mfrow = c(1,1))
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,150), data = data.lf, ylab = "lf for Cx. mosquitoes", xlab = "Temperature")
points(trait ~ T, data = data.lf.Cpip, col = "grey")
points(trait ~ T, data = data.lf.Cqui, col = "red")
points(trait ~ T, data = data.lf.Ctar, col = "blue")
points(trait ~ T, data = data.lf.Cmol, col = "purple")
points(trait ~ T, data = data.lf.Cpal, col = "magenta")
points(trait ~ T, data = data.lf.Cres, col = "green")
points(trait ~ T, data = data.lf.Atae, col = "violet")

summary(data.lf.Cpip$T)
summary(data.lf.Cqui$T)
summary(data.lf.Ctar$T)
summary(data.lf.Atae$T)


##########
###### 2. JAGS Models
##########

############## Linear Model with uniform priors

# NOTE: The derived quantity has an extra term to prevent it from going negative
# NOTE 2: The slope is specified as negative, because we use a gamma distribution to fit the priors, and gamma can only be positive

sink("linear.txt")
cat("
    model{
    
    ## Priors
    cf.m ~ dunif(-10, 10)
    cf.b ~ dunif(0, 250)
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.m * temp[i] + cf.b
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (-1 * cf.m * Temp.xs[i] + cf.b) * (Temp.xs[i] < cf.b / cf.m )
    }
    
    } # close model
    ",fill=T)
sink()


############## Linear Model with gamma priors - slope has an added negative sign because the slope is negative but gamma prior can only be positive

# NOTE: The derived quantity has an extra term to prevent it from going negative
# NOTE 2: The slope is specified as negative, because we use a gamma distribution to fit the priors, and gamma can only be positive

sink("linear_inf.txt")
cat("
    model{
    
    ## Priors
    cf.m ~ dgamma(hypers[1,1], hypers[2,1])
    cf.b ~ dgamma(hypers[1,2], hypers[2,2])
    cf.sigma ~ dunif(0, 1000)
    cf.tau <- 1 / (cf.sigma * cf.sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- -1 * cf.m * temp[i] + cf.b
    trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }
    
    ## Derived Quantities and Predictions
    for(i in 1:N.Temp.xs){
    z.trait.mu.pred[i] <- (-1 * cf.m * Temp.xs[i] + cf.b) * (Temp.xs[i] < cf.b / cf.m)
    }
    
    } # close model
    ",fill=T)
sink()


##########
###### 3. Shared settings for all models
##########

##### inits Function for linear
inits<-function(){list(
  cf.m = 2,
  cf.b = 35,
  cf.sigma = rlnorm(1))}

##### Parameters to Estimate for linear
parameters <- c("cf.m", "cf.b", "cf.sigma", "z.trait.mu.pred")

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
###### 4. Fit traits with uniform priors
##########

###################################### lifespan for Cx pipiens

##### Set data
data <- data.lf.Cpip

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.Cpip.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Cpip.out$BUGSoutput$summary[1:5,]
mcmcplot(lf.Cpip.out)

save(lf.Cpip.out, file = "jagsout_lf_Cpip.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,150), data = data.lf.Cpip, ylab = "lifespan for Cx pipiens", xlab = "Temperature")
lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### lifespan for Cx quinquefasciatus

##### Set data
data <- data.lf.Cqui

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.Cqui.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Cqui.out$BUGSoutput$summary[1:5,]
mcmcplot(lf.Cqui.out)

save(lf.Cqui.out, file = "jagsout_lf_Cqui.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,150), data = data.lf.Cqui, ylab = "lifespan for Cx quinquefacsiatus", xlab = "Temperature")
lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### lifespan for Cx tarsalis

##### Set data
data <- data.lf.Ctar

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.Ctar.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Ctar.out$BUGSoutput$summary[1:5,]
mcmcplot(lf.Ctar.out)

save(lf.Ctar.out, file = "jagsout_lf_Ctar.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,90), data = data.lf.Ctar, ylab = "lifespan for Cx tarsalis", xlab = "Temperature")
lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### lifespan for Ae. taenorhynchus

##### Set data
data <- data.lf.Atae

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.Atae.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Atae.out$BUGSoutput$summary[1:5,]
mcmcplot(lf.Atae.out)

save(lf.Atae.out, file = "jagsout_lf_Atae.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,80), data = data.lf.Atae, ylab = "lifespan for Cx tarsalis", xlab = "Temperature")
lines(lf.Atae.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Atae.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Atae.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 5. Fit thermal responses for lifespan priors using a leave-one-out approach
##########

# Subset Data - leave-one-out
data.lf.Cpip.prior <- subset(data.lf, host.code != "Cpip" | host.code == "Cpip" & Trait == "Male")
data.lf.Cqui.prior <- subset(data.lf, host.code != "Cqui" | host.code == "Cqui" & Trait == "Male")
data.lf.Ctar.prior <- subset(data.lf, host.code != "Ctar" | host.code == "Ctar" & Trait == "Male")
data.lf.Atae.prior <- subset(data.lf, host.code != "Atae" | host.code == "Atae" & Trait == "Male")

###################################### lifespan for Cx pipiens prior

##### Set data
data <- data.lf.Cpip.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.Cpip.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Cpip.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(lf.Cpip.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,120), data = data.lf.Cpip.prior, ylab = "lifespan for Cx pipiens prior", xlab = "Temperature")
lines(lf.Cpip.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cpip.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cpip.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### lifespan for Cx quinquefasciatus prior

##### Set data
data <- data.lf.Cqui.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.Cqui.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Cqui.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(lf.Cqui.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,120), data = data.lf.Cqui.prior, ylab = "lifespan for Cx quinqefasciatus prior", xlab = "Temperature")
lines(lf.Cqui.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cqui.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cqui.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### lifespan for Cx tarsalis prior

##### Set data
data <- data.lf.Ctar.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.Ctar.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Ctar.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(lf.Ctar.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,120), data = data.lf.Ctar.prior, ylab = "lifespan for Cx tarsalis prior", xlab = "Temperature")
lines(lf.Ctar.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Ctar.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Ctar.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


###################################### lifespan for Ae taeniorhynchus prior

##### Set data
data <- data.lf.Atae.prior

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
lf.Atae.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear.txt",
                          n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Atae.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(lf.Atae.prior.out)

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,120), data = data.lf.Atae.prior, ylab = "lifespan for Ae taeniorhynchus prior", xlab = "Temperature")
lines(lf.Atae.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Atae.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Atae.prior.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 6. Fit gamma distributions to lf prior thermal responses
##########
# Note that slope is getting multiplied by negative 1 because gamma dist must be positive

###################################### lf for Cpip prior

# Get the posterior dists for all 4 parameters into a data frame
lf.Cpip.prior.cf.dists <- data.frame(m = as.vector(lf.Cpip.prior.out$BUGSoutput$sims.list$cf.m),
                                     b = as.vector(lf.Cpip.prior.out$BUGSoutput$sims.list$cf.b))

# Fit gamma distributions for each parameter posterior dists
lf.Cpip.prior.gamma.fits = apply(lf.Cpip.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### lf for Cqui prior

# Get the posterior dists for all 4 parameters into a data frame
lf.Cqui.prior.cf.dists <- data.frame(m = as.vector(lf.Cqui.prior.out$BUGSoutput$sims.list$cf.m),
                                     b = as.vector(lf.Cqui.prior.out$BUGSoutput$sims.list$cf.b))

# Fit gamma distributions for each parameter posterior dists
lf.Cqui.prior.gamma.fits = apply(lf.Cqui.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### lf for Cqui prior

# Get the posterior dists for all 4 parameters into a data frame
lf.Ctar.prior.cf.dists <- data.frame(m = as.vector(lf.Ctar.prior.out$BUGSoutput$sims.list$cf.m),
                                     b = as.vector(lf.Ctar.prior.out$BUGSoutput$sims.list$cf.b))

# Fit gamma distributions for each parameter posterior dists
lf.Ctar.prior.gamma.fits = apply(lf.Ctar.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### lf for Atae prior

# Get the posterior dists for all 4 parameters into a data frame
lf.Atae.prior.cf.dists <- data.frame(m = as.vector(lf.Atae.prior.out$BUGSoutput$sims.list$cf.m),
                                     b = as.vector(lf.Atae.prior.out$BUGSoutput$sims.list$cf.b))

# Fit gamma distributions for each parameter posterior dists
lf.Atae.prior.gamma.fits = apply(lf.Atae.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

lf.hypers <- list(lf.Cpip.prior.gamma.fits, lf.Cqui.prior.gamma.fits, lf.Ctar.prior.gamma.fits, lf.Atae.prior.gamma.fits)
save(lf.hypers, file = "lfhypers.Rsave")


##########
###### 7. Fit lf thermal responses with data-informed priors
##########

load("lfhypers.Rsave")
lf.Cpip.prior.gamma.fits <- lf.hypers[[1]]
lf.Cqui.prior.gamma.fits <- lf.hypers[[2]]
lf.Ctar.prior.gamma.fits <- lf.hypers[[3]]
lf.Atae.prior.gamma.fits <- lf.hypers[[4]]

# Exclude males
data.lf.Cqui.fem <- subset(data.lf, host.code == "Cqui" & Trait != "Male")
data.lf.Ctar.fem <- subset(data.lf, host.code == "Ctar" & Trait != "Male")
data.lf.Atae.fem <- subset(data.lf, host.code == "Atae" & Trait != "Male")
data.lf.Cpip.fem <- subset(data.lf, host.code == "Cpip" & Trait != "Male") # For some reason this is also excluding NAs, so work around below
which(data.lf.Cpip$Trait == "Male")
nrow(data.lf.Cpip)
data.lf.Cpip.fem <- rbind(data.lf.Cpip[1:15,], data.lf.Cpip[21:37,])

nrow(subset(data.lf, host.code == "Cpip" & Trait == "Male"))
nrow(subset(data.lf, host.code == "Cqui" & Trait == "Male"))
nrow(subset(data.lf, host.code == "Ctar" & Trait == "Male"))
nrow(subset(data.lf, host.code == "Atae" & Trait == "Male"))

###################################### lifespan for Cx pipiens

##### Set data
data <- data.lf.Cpip.fem
hypers <- lf.Cpip.prior.gamma.fits * 0.01

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
lf.Cpip.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear_inf.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Cpip.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(lf.Cpip.out.inf)

save(lf.Cpip.out.inf, file = "jagsout_lf_Cpip_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,150), data = data.lf.Cpip.fem, ylab = "lifespan for Cx pipiens", xlab = "Temperature")
lines(lf.Cpip.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cpip.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cpip.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)


###################################### lifespan for Cx quinquefasciatus

##### Set data
data <- data.lf.Cqui.fem
hypers <- lf.Cqui.prior.gamma.fits * 0.01

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
lf.Cqui.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Cqui.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(lf.Cqui.out.inf)

save(lf.Cqui.out.inf, file = "jagsout_lf_Cqui_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,150), data = data.lf.Cqui.fem, ylab = "lifespan for Cx quinquefacsiatus", xlab = "Temperature")
lines(lf.Cqui.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cqui.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Cqui.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)


###################################### lifespan for Cx tarsalis

##### Set data
data <- data.lf.Ctar.fem
hypers <- lf.Ctar.prior.gamma.fits * 0.01

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
lf.Ctar.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Ctar.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(lf.Ctar.out.inf)

save(lf.Ctar.out.inf, file = "jagsout_lf_Ctar_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,150), data = data.lf.Ctar.fem, ylab = "lifespan for Cx tarsalis", xlab = "Temperature")
lines(lf.Ctar.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Ctar.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Ctar.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)


###################################### lifespan for Ae taeniorhynchus

##### Set data
data <- data.lf.Atae.fem
hypers <- lf.Atae.prior.gamma.fits * 0.1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
lf.Atae.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="linear_inf.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
lf.Atae.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(lf.Atae.out.inf)

save(lf.Atae.out.inf, file = "jagsout_lf_Atae_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,150), data = data.lf.Atae.fem, ylab = "lifespan for Ae taeniorhynchus", xlab = "Temperature")
lines(lf.Atae.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Atae.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Atae.out.inf$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)


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

Cpip.lf.avg <- CalcTraitAvg(data.lf.Cpip.fem)
Cqui.lf.avg <- CalcTraitAvg(data.lf.Cqui.fem)
Ctar.lf.avg <- CalcTraitAvg(data.lf.Ctar.fem)
Atae.lf.avg <- CalcTraitAvg(data.lf.Atae.fem)


##########
###### 9. Plot Figures
##########

plot(trait ~ T, xlim = c(5, 45), ylim = c(0,80), data = data.lf.Ctar, ylab = "lifespan for Cx spp", xlab = "Temperature")
points(trait ~ T, data = data.lf.Cpip, col = "red")
points(trait ~ T, data = data.lf.Cqui, col = "blue")

lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "blue")

par(mfrow = c(2,3), mar = c(3, 4.5, 2, 1), oma = c(2, 0, 0, 0))

##### lf for Cx. pipiens pipiens
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,82), data = data.lf.Cpippip, xaxt = "n", pch = 19,
     ylab = "lifespan (days)", xlab = "", main = expression(paste(italic(Cx.)," ",italic(pipiens)," ",italic(pipiens))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(lf.Cpippip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Cpippip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Cpippip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### lf for Cx. pipiens molestus
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,150), data = data.lf.Cpipmol, xaxt = "n", pch = 19,
     ylab = "lifespan (days)", xlab = "", main = expression(paste(italic(Cx.)," ",italic(pipiens)," ",italic(molestus))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(lf.Cpipmol.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Cpipmol.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Cpipmol.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### lf for Cx. pipiens (both)
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,82), data = data.lf.Cpip, xaxt = "n", pch = 19,
     ylab = "lifespan (days)", xlab = "", main = expression(paste(italic(Cx.)," ",italic(pipiens)," (both)")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Cpip.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### lf for Cx. quinquefasciatus 
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,82), data = data.lf.Cqui, xaxt = "n", pch = 19,
     ylab = "lifespan (days)", xlab = "", main = expression(paste(italic(Cx.)," ",italic(quinquefasciatus))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Cqui.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### lf for Cx. tarsalis 
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,82), data = data.lf.Ctar, xaxt = "n", pch = 19,
     ylab = "lifespan (days)", xlab = "", main = expression(paste(italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Ctar.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### lf for Ae. taenorhynchus 
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,82), data = data.lf.Atae, xaxt = "n", pch = 19,
     ylab = "lifespan (days)", xlab = "", main = expression(paste(italic(Ae.)," ",italic(taenorhynchus))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(lf.Atae.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Atae.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(lf.Atae.out$BUGSoutput$summary[5:(5 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))