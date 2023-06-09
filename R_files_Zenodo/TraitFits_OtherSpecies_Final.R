## Marta Shocket, Stanford University
## Updated April 2018
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for mosquito and pathogen traits for RRV model
##            from alternate vector and virus species
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit trait responses for alternate vector and virus species with uniform and data-informed priors
##           5) Plot Appendix Figures S2, S3, and S4


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("~/Dropbox/Research Mordecai Lab/VBD Temp Project/RRV/Fitting Traits")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')

# Load Data for fitting traits
data.RRV <- read.csv("RRVTraitData.csv") # Data from database for all mosquito traits, RRV, and MVE

##### Subset data
data.pLA.Anot <- subset(data.RRV, trait.name == "pEA" & host.code == "Anot")
data.MDR.Anot <- subset(data.RRV, trait.name == "1/MDR" & host.code == "Anot")
data.pLA.Ocam <- subset(data.RRV, trait.name == "pEA" & host.code == "Ocam")
data.MDR.Ocam <- subset(data.RRV, trait.name == "1/MDR" & host.code == "Ocam")
data.bc.MVE <- subset(data.RRV, trait.name == "bc" & paras.code == "MVE")
data.PDR.MVE <- subset(data.RRV, trait.name == "EIP" & paras.code == "MVE")

# Load saved prior fits (gamma dists fit to parameter posteriors)
load("PriorFitsList.Rdata")
pEA.prior.gamma.fits  <- as.data.frame(prior.fits.list[4])
MDR.prior.gamma.fits  <- as.data.frame(prior.fits.list[5])
PDR.prior.gamma.fits  <- as.data.frame(prior.fits.list[6])


##########
###### 2. JAGS Models
##########

# NOTE: only new models unique to this file are listed here: basic models are in the Uniform Priors and DataInformedPriors R code files

############## Quadratic Model with gamma priors and truncated so derived quantities always =< 1

sink("quadprob_inf.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dgamma(hypers[1,1], hypers[2,1])
    cf.T0 ~ dgamma(hypers[1,2], hypers[2,2])
    cf.Tm ~ dgamma(hypers[1,3], hypers[2,3])
    cf.sigma ~ dgamma(hypers[1,4], hypers[2,4])
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

############## Quadratic Model with uniform priors and higher Tm upper limit (for MVEV bc)

sink("quadmax.txt")
cat("
    model{
    
    ## Priors
    cf.q ~ dunif(0, 1)
    cf.T0 ~ dunif(0, 24)
    cf.Tm ~ dunif(25, 55)
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
Temp.xs <- seq(5, 45, 0.2)
N.Temp.xs <-length(Temp.xs)


##########
###### 4. Fit Aedes notoscriptus, Ochlerotatus camptorhynchus, and MVEV traits with uniform and data-informed priors
##########

############## pEA (pLA - larval-adult survival) for Ae. notoscriptus - quadratic, uninformative

##### Set data
data <- data.pLA.Anot
hypers <- pEA.prior.gamma.fits * 0.1

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data and Run JAGS for uniform priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

pLA.Anot.out<-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Bundle Data and Run JAGS for data-informed priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

pLA.Anot.out.inf<-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob_inf.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Anot.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Anot.out)

pLA.Anot.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Anot.out.inf)

save(pLA.Anot.out, file = "jagsout_pLA_Anot.Rdata")
save(pLA.Anot.out.inf, file = "jagsout_pLA_Anot_inf.Rdata")

# Plot data + fits
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLA.Anot, ylab = "pLA for Ae. notoscriptus", xlab = "Temperature")
lines(pLA.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(pLA.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLA.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLA.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

# Calculate optimum for pLA: 22.8 C
Temp.xs[which.max(as.vector(pLA.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## MDR for Ae. notoscriptus - Briere

##### Set data
data <- data.MDR.Anot
hypers <- MDR.prior.gamma.fits

##### Organize Data for JAGS
trait <- 1/data$trait # take the inverse since it's given in development time and we want rate
N.obs <- length(trait)
temp <- data$T

##### Bundle Data and Run JAGS for uniform priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

MDR.Anot.out<-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Bundle Data and Run JAGS for data-informed priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

MDR.Anot.out.inf<-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Anot.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Anot.out)

MDR.Anot.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Anot.out.inf)

save(MDR.Anot.out, file = "jagsout_MDR_Anot.Rdata")
save(MDR.Anot.out.inf, file = "jagsout_MDR_Anot_inf.Rdata")

# Plot data + fits
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,.25), data = data.MDR.Anot, ylab = "MDR for Ae. notoscriptus", xlab = "Temperature")
lines(MDR.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(MDR.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(MDR.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(MDR.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

# Calculate optimum for MDR: 32.0 C
Temp.xs[which.max(as.vector(MDR.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## pEA (pLA - larval-adult survival) for O. camptorhynchus - quadratic

##### Set data
data <- data.pLA.Ocam
hypers <- pEA.prior.gamma.fits * 0.01

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data and Run JAGS for uniform priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

pLA.Ocam.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quadprob.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Bundle Data and Run JAGS for data-informed priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

pLA.Ocam.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Ocam.out$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Ocam.out)

pLA.Ocam.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Ocam.out.inf)

save(pLA.Ocam.out, file = "jagsout_pLA_Ocam.Rdata")
save(pLA.Ocam.out.inf, file = "jagsout_pLA_Ocam_inf.Rdata")

# Plot trait data + fits
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1.2), data = data.pLA.Ocam, ylab = "pLA for O. camptorhynchus", xlab = "Temperature")
lines(pLA.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(pLA.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLA.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(pLA.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

# Calculate optimum for pLA: 22.8 C
Temp.xs[which.max(as.vector(pLA.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## MDR for Ae. camptorhynchus - briere

##### Set data
data <- data.MDR.Ocam
hypers <- MDR.prior.gamma.fits * 0.1
hypers[,3] <- MDR.prior.gamma.fits[,3]

##### Organize Data for JAGS
trait <- 1/data$trait # take the inverse since it's given in development time and we want rate
N.obs <- length(trait)
temp <- data$T

##### Bundle Data and Run JAGS for uniform priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

MDR.Ocam.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Bundle Data and Run JAGS for data-informed priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

MDR.Ocam.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Ocam.out$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Ocam.out)

MDR.Ocam.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Ocam.out.inf)

save(MDR.Ocam.out, file = "jagsout_MDR_Ocam.Rdata")
save(MDR.Ocam.out.inf, file = "jagsout_MDR_Ocam_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.2), data = data.MDR.Ocam, ylab = "MDR for O. camptorhynchus", xlab = "Temperature")
lines(MDR.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(MDR.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(MDR.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(MDR.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

# Calculate optimum for MDR: 32.2 C
Temp.xs[which.max(as.vector(MDR.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## PDR for MVE in C. annulirostris - Briere

##### Set data
data <- data.PDR.MVE
hypers <- PDR.prior.gamma.fits

##### Organize Data for JAGS
trait <- 1/data$trait # take the inverse since it's given in development time and we want rate
N.obs <- length(trait)
temp <- data$T

##### Bundle Data and Run JAGS for uniform priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

PDR.MVE.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Bundle Data and Run JAGS for data-informed priors
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

PDR.MVE.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.MVE.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.MVE.out)

PDR.MVE.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.MVE.out.inf)

save(PDR.MVE.out, file = "jagsout_PDR_MVE_Cann.Rdata")
save(PDR.MVE.out.inf, file = "jagsout_PDR_MVE_Cann_inf.Rdata")

# plot trait data
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,.4), data = data.PDR.MVE, ylab = "PDR for MVE in C. annulirostrus", xlab = "Temperature")
lines(PDR.MVE.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MVE.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.MVE.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

lines(PDR.MVE.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.MVE.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.MVE.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

# Calculate optimum for PDR: 34.8 C
Temp.xs[which.max(as.vector(PDR.MVE.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


##########
###### 5. Plot Figures
##########

data.PDR.MVE <- subset(data.RRV, trait.name == "EIP" & paras.code == "MVE")
data.bc.MVE <- subset(data.RRV, trait.name == "bc" & paras.code == "MVE")
load("jagsout_PDR_MVE_Cann.Rdata")
load("jagsout_PDR_MVE_Cann_inf.Rdata")
load("jagsout_bc_MVE_Cann.Rdata")
load("jagsout_MDR_Ocam.Rdata")
load("jagsout_MDR_Ocam_inf.Rdata")
load("jagsout_pLA_Ocam.Rdata")
load("jagsout_pLA_Ocam_inf.Rdata")
load("jagsout_MDR_Anot.Rdata")
load("jagsout_MDR_Anot_inf.Rdata")
load("jagsout_pLA_Anot.Rdata")
load("jagsout_pLA_Anot_inf.Rdata")

############### Figure S2: Ae. camptorhynchus fits

par(mfrow = c(2,2), mar = c(3, 4.5, 2, 1), oma = c(1.5, 0, 3, 0))

# Ae. camptorhynchus MDR data-informed priors
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.125), data = data.MDR.Ocam, pch = 19, cex.lab = 1.15, xaxt = "n",
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Mosquito Development Rate (",italic(MDR),")")))
axis(1, at = seq(5, 45, 5))
lines(MDR.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 3, text = "Data Informed Priors", line = 2, cex = 1.2)
mtext(side = 3, text = expression(paste(italic(Ae.)," ",italic(camptorhynchus))), line = 3.3, cex = 1.4, adj = 3.8)
legend("topleft", legend = "A", bty= "n", cex = 1.2, adj = c(1.5, 0))

# Ae. camptorhynchus MDR uniform priors
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.125), data = data.MDR.Ocam, pch = 19, cex.lab = 1.15, xaxt = "n",
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Mosquito Development Rate (",italic(MDR),")")))
axis(1, at = seq(5, 45, 5))
lines(MDR.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 3, text = "Uniform Priors", line = 1.8, cex = 1.2)
legend("topleft", legend = "B", bty= "n", cex = 1.2, adj = c(1.5, 0))

# Ae. camptorhynchus pLA data-informed priors
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1.1), data = data.pLA.Ocam, pch = 19, cex.lab = 1.15, xaxt = "n",
     ylab = "Probability", xlab = "", main = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")))
axis(1, at = seq(5, 45, 5))
lines(pLA.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Ocam.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 1, text = expression(paste("Temperature (",degree,"C)")), line = 3)
legend("topleft", legend = "C", bty= "n", cex = 1.2, adj = c(1.5, 0))

# Ae. camptorhynchus pLA uniform priors
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1.1), data = data.pLA.Ocam, pch = 19, cex.lab = 1.15, xaxt = "n",
     ylab = "Probability", xlab = "", main = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")))
axis(1, at = seq(5, 45, 5))
lines(pLA.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Ocam.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 1, text = expression(paste("Temperature (",degree,"C)")), line = 3)
legend("topleft", legend = "D", bty= "n", cex = 1.2, adj = c(1.5, 0))

############### Figures S3: Ae. notoscriptus fits

par(mfrow = c(2,2), mar = c(3, 4.5, 2, 1), oma = c(1.5, 0, 3, 0))

# Ae. notoscriptus MDR data-informed priors
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0, .2), data = data.MDR.Anot, pch = 19, cex.lab = 1.15, xaxt = "n",
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Mosquito Development Rate (",italic(MDR),")")))
axis(1, at = seq(5, 45, 5))
lines(MDR.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 3, text = "Data Informed Priors", line = 2, cex = 1.2)
mtext(side = 3, text = expression(paste(italic(Ae.)," ",italic(notoscriptus))), line = 3.3, cex = 1.4, adj = 2.3)
legend("topleft", legend = "A", bty= "n", cex = 1.2, adj = c(1.5, 0))

# Ae. notoscriptus MDR uniform priors
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0, .2), data = data.MDR.Anot, pch = 19, cex.lab = 1.15, xaxt = "n",
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Mosquito Development Rate (",italic(MDR),")")))
axis(1, at = seq(5, 45, 5))
lines(MDR.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 3, text = "Uniform Priors", line = 2, cex = 1.2)
legend("topleft", legend = "B", bty= "n", cex = 1.2, adj = c(1.5, 0))

# Ae. notoscriptus pLA data-informed priors
plot(trait ~ T, xlim = c(5, 45), ylim = c(0, 1.1), data = data.pLA.Anot, pch = 19, cex.lab = 1.15, xaxt = "n", 
     ylab = "Probability", xlab = "", main = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")))
axis(1, at = seq(5, 45, 5))
lines(pLA.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Anot.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 1, text = expression(paste("Temperature (",degree,"C)")), line = 3)
legend("topleft", legend = "C", bty= "n", cex = 1.2, adj = c(1.5, 0))

# Ae. notoscriptus pLA uniform priors
plot(trait ~ T, xlim = c(5, 45), ylim = c(0, 1.1), data = data.pLA.Anot, pch = 19, cex.lab = 1.15, xaxt = "n", 
     ylab = "Probability", xlab = "", main = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")))
axis(1, at = seq(5, 45, 5))
lines(pLA.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Anot.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 1, text = expression(paste("Temperature (",degree,"C)")), line = 3)
legend("topleft", legend = "D", bty= "n", cex = 1.2, adj = c(1.5, 0))

############### Figure S4: Murray Valley Encephalitis virus fits

data1 <- subset(data.bc.MVE, T == 20)
data2 <- subset(data.bc.MVE, T == 27)
data3 <- subset(data.bc.MVE, T == 33.5)

bc.means <- data.frame(T = c(20, 27, 33.5), trait = c(mean(data1$trait), mean(data2$trait), mean(data3$trait)))
points(trait ~ T, data = bc.means, col = "red")

par(mfrow = c(2,2), mar = c(3, 4.5, 2, 1), oma = c(1.5, 0, 3, 0))

# PDR data-informed priors
plot(1/trait ~ T, xlim = c(7, 42), ylim = c(0,0.3), data = data.PDR.MVE, cex.lab = 1.15, xaxt = "n",  pch = 19,
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Parasite Development Rate (",italic(PDR),")")))
axis(1, at = seq(5, 45, 5))
lines(PDR.MVE.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.MVE.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.MVE.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 3, text = "Data Informed Priors", line = 2, cex = 1.2)
mtext(side = 3, text = "Murray Valley Encephalitis virus", line = 3.3, cex = 1.4, adj = -9.6)
legend("topleft", legend = "A", bty= "n", cex = 1.2, adj = c(1.5, 0))

# PDR uniform priors
plot(1/trait ~ T, xlim = c(7, 42), ylim = c(0,0.3), data = data.PDR.MVE, cex.lab = 1.15, xaxt = "n",  pch = 19,
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Parasite Development Rate (",italic(PDR),")")))
axis(1, at = seq(5, 45, 5))
lines(PDR.MVE.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.MVE.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.MVE.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
mtext(side = 3, text = "Uniform Priors", line = 2, cex = 1.2)
mtext(side = 1, text = expression(paste("Temperature (",degree,"C)")), line = 3)
legend("topleft", legend = "B", bty= "n", cex = 1.2, adj = c(1.5, 0))

# bc no fit
plot(trait ~ T, xlim = c(7, 42), ylim = c(0,1.2), data = data.bc.MVE, xaxt = "n", pch = 19,
     ylab = "Probability", xlab = "", main = expression(paste("Transmission Probability (",italic(bc),")")), cex.lab = 1.15)
points(trait ~ T, data = bc.means, pch = 21, bg = "gray")
axis(1, at = seq(5, 45, 5))
mtext(side = 1, text = expression(paste("Temperature (",degree,"C)")), line = 3)
legend("topleft", legend = "C", bty= "n", cex = 1.2, adj = c(1.5, 0))