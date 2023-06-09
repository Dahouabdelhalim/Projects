## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for parasite development rate (PDR) for six temperature arboviruses
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit PDR thermal responses with uniform priors
##           5) Fit PDR thermal responses for priors
##           6) Fit gamma distributions to PDR prior thermal responses
##           7) Fit PDR thermal responses with data-informed priors
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

# Load Data
data.PDR <- read.csv("TraitData_PDR.csv") # Data from database for most traits (except below)
unique(data.PDR$joint.code)

# Subset Data
data.PDR.CpipWNV <- subset(data.PDR, joint.code == "CpipWNV")
data.PDR.CtarWNV <- subset(data.PDR, joint.code == "CtarWNV")
data.PDR.CquiWNV <- subset(data.PDR, joint.code == "CquiWNV")
data.PDR.CuniWNV <- subset(data.PDR, joint.code == "CuniWNV")

data.PDR.CtarWEEV <- subset(data.PDR, joint.code == "CtarWEEV")
data.PDR.CtarSLEV <- subset(data.PDR, joint.code == "CtarSLEV")

data.PDR.AtaeSINV <- subset(data.PDR, joint.code == "AtaeSINV")
data.PDR.AtaeRVFV <- subset(data.PDR, joint.code == "AtaeRVFV")
data.PDR.AtriEEEV <- subset(data.PDR, joint.code == "AtriEEEV")

par(mfrow = c(1,1))
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.4), data = data.PDR, ylab = "PDR", xlab = "Temperature")
points(1/trait ~ T, data = data.PDR.CpipWNV, col = "grey")
points(1/trait ~ T, data = data.PDR.CtarWNV, col = "blue")
points(1/trait ~ T, data = data.PDR.CquiWNV, col = "red")
points(1/trait ~ T, data = data.PDR.CuniWNV, col = "orange")
points(1/trait ~ T, data = data.PDR.CtarWEEV, col = "dodgerblue")
points(1/trait ~ T, data = data.PDR.CtarSLEV, col = "navyblue")
points(1/trait ~ T, data = data.PDR.AtaeSINV, col = "green")
points(1/trait ~ T, data = data.PDR.AtaeRVFV, col = "darkgreen")
points(1/trait ~ T, data = data.PDR.AtriEEEV, col = "purple")


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
###### 4. Fit PDR thermal responses with uniform priors
##########

###################################### PDR for WNV in Cx pipiens - Briere

##### Set data
data <- data.PDR.CpipWNV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CpipWNV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CpipWNV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CpipWNV.out)

save(PDR.CpipWNV.out, file = "jagsout_PDR_CpipWNV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CpipWNV, ylab = "PDR for WNV in Cx pipiens", xlab = "Temperature")
lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for WNV in Cx quinquefasciatus - Briere

##### Set data
data <- data.PDR.CquiWNV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CquiWNV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CquiWNV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CquiWNV.out)

save(PDR.CquiWNV.out, file = "jagsout_PDR_CquiWNV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CquiWNV, ylab = "PDR for WNV in Cx pipiens", xlab = "Temperature")
lines(PDR.CquiWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CquiWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CquiWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for WNV in Cx tarsalis - Briere

##### Set data
data <- data.PDR.CtarWNV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CtarWNV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarWNV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarWNV.out)

save(PDR.CtarWNV.out, file = "jagsout_PDR_CtarWNV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CtarWNV, ylab = "PDR for WNV in Cx tarsalis", xlab = "Temperature")
lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for WEEV Cx tarsalis - Briere

##### Set data
data <- data.PDR.CtarWEEV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CtarWEEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarWEEV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarWEEV.out)

save(PDR.CtarWEEV.out, file = "jagsout_PDR_CtarWEEV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CtarWEEV, ylab = "PDR for WEEV in Cx tarsalis", xlab = "Temperature")
lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for SLEV in Cx tarsalis - Briere

##### Set data
data <- data.PDR.CtarSLEV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CtarSLEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarSLEV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarSLEV.out)

save(PDR.CtarSLEV.out, file = "jagsout_PDR_CtarSLEV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CtarSLEV, ylab = "PDR for SLEV in Cx tarsalis", xlab = "Temperature")
lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for WNV in Cx univittatus - Briere

##### Set data
data <- data.PDR.CuniWNV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CuniWNV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                        n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CuniWNV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CuniWNV.out)

save(PDR.CuniWNV.out, file = "jagsout_PDR_CuniWNV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CuniWNV, ylab = "PDR for WNV in Cx univittatus", xlab = "Temperature")
lines(PDR.CuniWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CuniWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CuniWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for SINV in Ae taenorhynchus - Briere - fit is not good! Do not use this fit

##### Set data
data <- data.PDR.AtaeSINV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.AtaeSINV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())


##### Examine Output
PDR.AtaeSINV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.AtaeSINV.out)

save(PDR.AtaeSINV.out, file = "jagsout_PDR_AtaeSINV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.AtaeSINV, ylab = "PDR for SINV in Ae. taenorhynchus", xlab = "Temperature")
lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for RVFV in Ae taenorhynchus - Briere

##### Set data
data <- data.PDR.AtaeRVFV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.AtaeRVFV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.AtaeRVFV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.AtaeRVFV.out)

save(PDR.AtaeRVFV.out, file = "jagsout_PDR_AtaeRVFV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.PDR.AtaeRVFV, ylab = "PDR for RVFV in Ae. taenorhynchus", xlab = "Temperature")
lines(PDR.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeRVFV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for EEEV in Ae. triseriatus - Briere

##### Set data
data <- data.PDR.AtriEEEV

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.AtriEEEV.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                     n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.AtriEEEV.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.AtriEEEV.out)

save(PDR.AtriEEEV.out, file = "jagsout_PDR_AtriEEEV.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.AtriEEEV, ylab = "PDR for EEEV in Ae triseriatus", xlab = "Temperature")
lines(PDR.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtriEEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)



##########
###### 5. Fit thermal responses for PDR priors using a leave-one-out approach
##########

# Subset Data
data.PDR.CpipWNV.prior <- subset(data.PDR, joint.code != "CpipWNV")
data.PDR.CquiWNV.prior <- subset(data.PDR, joint.code != "CquiWNV")
data.PDR.CtarWNV.prior <- subset(data.PDR, joint.code != "CtarWNV")
data.PDR.CuniWNV.prior <- subset(data.PDR, joint.code != "CuniWNV")

data.PDR.CtarWEEV.prior <- subset(data.PDR, joint.code != "CtarWEEV")
data.PDR.CtarSLEV.prior <- subset(data.PDR, joint.code != "CtarSLEV")

data.PDR.AtaeSINV.prior <- subset(data.PDR, joint.code != "AtaeSINV")
data.PDR.AtaeRVFV.prior <- subset(data.PDR, joint.code != "AtaeRVFV")
data.PDR.AtriEEEV.prior <- subset(data.PDR, joint.code != "AtriEEEV")

###################################### PDR for WNV in Cpip prior - Briere

##### Set data
data <- data.PDR.CpipWNV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CpipWNV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CpipWNV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CpipWNV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.CpipWNV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CpipWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for WNV in Cqui prior - Briere

##### Set data
data <- data.PDR.CquiWNV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CquiWNV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CquiWNV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CquiWNV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.CquiWNV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.CquiWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CquiWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CquiWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for WNV in Ctar prior - Briere

##### Set data
data <- data.PDR.CtarWNV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CtarWNV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarWNV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarWNV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.CtarWNV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.CtarWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for WEEV in Ctar prior - Briere

##### Set data
data <- data.PDR.CtarWEEV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CtarWEEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarWEEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarWEEV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.CtarWEEV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for SLEV in Ctar prior - Briere

##### Set data
data <- data.PDR.CtarSLEV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CtarSLEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                               n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarSLEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarSLEV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.CtarSLEV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarSLEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for WNV in Cuni prior - Briere

##### Set data
data <- data.PDR.CuniWNV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.CuniWNV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CuniWNV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CuniWNV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.CuniWNV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.CuniWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CuniWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CuniWNV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for EEEV in Atri prior - Briere

##### Set data
data <- data.PDR.AtriEEEV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.AtriEEEV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.AtriEEEV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.AtriEEEV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.AtriEEEV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.AtriEEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtriEEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtriEEEV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for SINV in Atae prior - Briere

##### Set data
data <- data.PDR.AtaeSINV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.AtaeSINV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                               n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.AtaeSINV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.AtaeSINV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.AtaeSINV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.AtaeSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeSINV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

###################################### PDR for RVFV in Atae prior - Briere

##### Set data
data <- data.PDR.AtaeRVFV.prior

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS
PDR.AtaeRVFV.prior.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                               n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.AtaeRVFV.prior.out$BUGSoutput$summary[1:5,]
mcmcplot(PDR.AtaeRVFV.prior.out)

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR.AtaeRVFV.prior, ylab = "PDR for all species", xlab = "Temperature")
lines(PDR.AtaeRVFV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeRVFV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeRVFV.prior.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


##########
###### 6. Fit gamma distributions to PDR prior thermal responses
##########

###################################### PDR for CpipWNV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.CpipWNV.prior.cf.dists <- data.frame(q = as.vector(PDR.CpipWNV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.CpipWNV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.CpipWNV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.CpipWNV.prior.gamma.fits = apply(PDR.CpipWNV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### PDR for CquiWNV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.CquiWNV.prior.cf.dists <- data.frame(q = as.vector(PDR.CquiWNV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.CquiWNV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.CquiWNV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.CquiWNV.prior.gamma.fits = apply(PDR.CquiWNV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### PDR for CtarWNV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.CtarWNV.prior.cf.dists <- data.frame(q = as.vector(PDR.CtarWNV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.CtarWNV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.CtarWNV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.CtarWNV.prior.gamma.fits = apply(PDR.CtarWNV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### PDR for CuniWNV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.CuniWNV.prior.cf.dists <- data.frame(q = as.vector(PDR.CuniWNV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.CuniWNV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.CuniWNV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.CuniWNV.prior.gamma.fits = apply(PDR.CuniWNV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### PDR for CtarSLEV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.CtarSLEV.prior.cf.dists <- data.frame(q = as.vector(PDR.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.CtarSLEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.CtarSLEV.prior.gamma.fits = apply(PDR.CtarSLEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### PDR for CtarWEEV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.CtarWEEV.prior.cf.dists <- data.frame(q = as.vector(PDR.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.CtarWEEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.CtarWEEV.prior.gamma.fits = apply(PDR.CtarWEEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### PDR for AtriEEEV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.AtriEEEV.prior.cf.dists <- data.frame(q = as.vector(PDR.AtriEEEV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.AtriEEEV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.AtriEEEV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.AtriEEEV.prior.gamma.fits = apply(PDR.AtriEEEV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### PDR for AtaeSINV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.AtaeSINV.prior.cf.dists <- data.frame(q = as.vector(PDR.AtaeSINV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.AtaeSINV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.AtaeSINV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.AtaeSINV.prior.gamma.fits = apply(PDR.AtaeSINV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

###################################### PDR for AtaeRVFV prior

# Get the posterior dists for 3 main parameters (not sigma) into a data frame
PDR.AtaeRVFV.prior.cf.dists <- data.frame(q = as.vector(PDR.AtaeRVFV.prior.out$BUGSoutput$sims.list$cf.q),
                                         T0 = as.vector(PDR.AtaeRVFV.prior.out$BUGSoutput$sims.list$cf.T0), 
                                         Tm = as.vector(PDR.AtaeRVFV.prior.out$BUGSoutput$sims.list$cf.Tm))

# Fit gamma distributions for each parameter posterior dists
PDR.AtaeRVFV.prior.gamma.fits = apply(PDR.AtaeRVFV.prior.cf.dists, 2, function(df) fitdistr(df, "gamma")$estimate)

PDR.hypers <- list(PDR.CpipWNV.prior.gamma.fits, PDR.CquiWNV.prior.gamma.fits, PDR.CtarWNV.prior.gamma.fits, PDR.CuniWNV.prior.gamma.fits,
                   PDR.CtarWEEV.prior.gamma.fits, PDR.CtarSLEV.prior.gamma.fits, 
                   PDR.AtriEEEV.prior.gamma.fits, PDR.AtaeSINV.prior.gamma.fits, PDR.AtaeRVFV.prior.gamma.fits)
save(PDR.hypers, file = "PDRhypers.Rsave")


##########
###### 7. Fit PDR thermal responses with data-informed priors
##########

load("PDRhypers.Rsave")
PDR.CpipWNV.prior.gamma.fits <- PDR.hypers[[1]]
PDR.CquiWNV.prior.gamma.fits <- PDR.hypers[[2]]
PDR.CtarWNV.prior.gamma.fits <- PDR.hypers[[3]]
PDR.CuniWNV.prior.gamma.fits <- PDR.hypers[[4]]
PDR.CtarWEEV.prior.gamma.fits <- PDR.hypers[[5]]
PDR.CtarSLEV.prior.gamma.fits <- PDR.hypers[[6]]
PDR.AtriEEEV.prior.gamma.fits <- PDR.hypers[[7]]
PDR.AtaeSINV.prior.gamma.fits <- PDR.hypers[[8]]
PDR.AtaeRVFV.prior.gamma.fits <- PDR.hypers[[9]]


############## PDR for WNV in Cx. pipiens - Briere

##### Set data 
data <- data.PDR.CpipWNV
hypers <- PDR.CpipWNV.prior.gamma.fits * 2
hypers[,2] <- PDR.CpipWNV.prior.gamma.fits[,2] * 5

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.CpipWNV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CpipWNV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CpipWNV.out.inf)

save(PDR.CpipWNV.out.inf, file = "jagsout_PDR_CpipWNV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CpipWNV, ylab = "PDR for WNV in Cx. pipiens", xlab = "Temperature", pch = 1)
lines(PDR.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.2)

# Get optimum for PDR: 37.6 C
Temp.xs[which.max(as.vector(PDR.CpipWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## PDR for WNV in Cx. quinquefasciatus - Briere

##### Set data 
data <- data.PDR.CquiWNV
hypers <- PDR.CquiWNV.prior.gamma.fits * 1

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.CquiWNV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CquiWNV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CquiWNV.out.inf)

save(PDR.CquiWNV.out.inf, file = "jagsout_PDR_CquiWNV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CquiWNV, ylab = "PDR for WNV in Cx. quinquefasciatus", xlab = "Temperature", pch = 1)
lines(PDR.CquiWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CquiWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CquiWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.2)

# Get optimum for PDR: 37.7 C
Temp.xs[which.max(as.vector(PDR.CquiWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## PDR for WNV in Cx. tarsalis - Briere

##### Set data 
data <- data.PDR.CtarWNV
hypers <- PDR.CtarWNV.prior.gamma.fits * 2

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.CtarWNV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarWNV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarWNV.out.inf)

save(PDR.CtarWNV.out.inf, file = "jagsout_PDR_CtarWNV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CtarWNV, ylab = "PDR for WNV in Cx. tarsalis", xlab = "Temperature", pch = 1)
lines(PDR.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.2)

# Get optimum for PDR: 36.9 C
Temp.xs[which.max(as.vector(PDR.CtarWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## PDR for WNV in Cx. univittatus - Briere

##### Set data 
data <- data.PDR.CuniWNV
hypers <- PDR.CuniWNV.prior.gamma.fits * 1
hypers[,2] <- PDR.CuniWNV.prior.gamma.fits[,2] * 3
hypers[,3] <- PDR.CuniWNV.prior.gamma.fits[,3] * .2

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.CuniWNV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CuniWNV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CuniWNV.out.inf)

save(PDR.CuniWNV.out.inf, file = "jagsout_PDR_CuniWNV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CuniWNV, ylab = "PDR for WNV in Cx. univittatus", xlab = "Temperature", pch = 1)
lines(PDR.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#lines(PDR.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"] ~ Temp.xs, col = "red")
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.2)

# Get optimum for PDR: 27.6 C
Temp.xs[which.max(as.vector(PDR.CuniWNV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## PDR for SLEV in Cx. tarsalis - Briere

##### Set data 
data <- data.PDR.CtarSLEV
hypers <- PDR.CtarSLEV.prior.gamma.fits * 1

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.CtarSLEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                            n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarSLEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarSLEV.out.inf)

save(PDR.CtarSLEV.out.inf, file = "jagsout_PDR_CtarSLEV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CtarSLEV, ylab = "PDR for SLEV in Cx. tarsalis", xlab = "Temperature", pch = 1)
lines(PDR.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.2)

# Get optimum for PDR: 37.5 C
Temp.xs[which.max(as.vector(PDR.CtarSLEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]

############## PDR for WEEV in Cx. tarsalis - Briere

##### Set data 
data <- data.PDR.CtarWEEV
hypers <- PDR.CtarWEEV.prior.gamma.fits * 1
hypers[,2] <- PDR.CtarWEEV.prior.gamma.fits[,2] * 0.05

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.CtarWEEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.CtarWEEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.CtarWEEV.out.inf)

save(PDR.CtarWEEV.out.inf, file = "jagsout_PDR_CtarWEEV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.CtarWEEV, ylab = "PDR for WEEV in Cx. tarsalis", xlab = "Temperature", pch = 1)
lines(PDR.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#lines(PDR.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"] ~ Temp.xs, col = "red")
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for PDR: 35.3 C
Temp.xs[which.max(as.vector(PDR.CtarWEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]
hist(PDR.CtarWEEV.out.inf$BUGSoutput$sims.list$cf.T0)

############## PDR for EEEV in Ae. triseriatus - Briere

##### Set data 
data <- data.PDR.AtriEEEV
hypers <- PDR.AtriEEEV.prior.gamma.fits * 2

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.AtriEEEV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.AtriEEEV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.AtriEEEV.out.inf)

save(PDR.AtriEEEV.out.inf, file = "jagsout_PDR_AtriEEEV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.AtriEEEV, ylab = "PDR for EEEV in Ae. triseriatus", xlab = "Temperature", pch = 1)
lines(PDR.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
#lines(PDR.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"] ~ Temp.xs, col = "red")
legend("topleft", legend = "G", bty = "n", adj = 1, cex = 1.2)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)

# Get optimum for PDR: 37.4 C
Temp.xs[which.max(as.vector(PDR.AtriEEEV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


# ############## PDR for SINV in Ae. taeniorhynchus - NO THERMAL RESPONSE FITS WELL
# 
# ##### Set data 
# data <- data.PDR.AtaeSINV
# hypers <- PDR.AtaeSINV.prior.gamma.fits * 1
# 
# ##### Organize Data for JAGS
# trait <- 1/data$trait
# N.obs <- length(trait)
# temp <- data$T
# 
# ##### Bundle Data
# jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)
# 
# ##### Run JAGS
# PDR.AtaeSINV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
#                              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
# 
# ##### Examine Output
# PDR.AtaeSINV.out.inf$BUGSoutput$summary[1:5,]
# mcmcplot(PDR.AtaeSINV.out.inf)
# 
# # Plot data + fit
# plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.PDR.AtaeSINV, ylab = "PDR for SINV in Ae. taeniorhynchus", xlab = "Temperature", pch = 19)
# lines(PDR.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
# lines(PDR.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
# lines(PDR.AtaeSINV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)


############## PDR for RVFV in Ae. taeniorhynchus - Briere

##### Set data 
data <- data.PDR.AtaeRVFV
hypers <- PDR.AtaeRVFV.prior.gamma.fits * 1
hypers[,1] <- PDR.AtaeRVFV.prior.gamma.fits[,1] * .2
hypers[,2] <- PDR.AtaeRVFV.prior.gamma.fits[,2] * 2
hypers[,3] <- PDR.AtaeRVFV.prior.gamma.fits[,3] * 2

##### Organize Data for JAGS
trait <- 1/data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.AtaeRVFV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                             n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.AtaeRVFV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.AtaeRVFV.out.inf)

save(PDR.AtaeRVFV.out.inf, file = "jagsout_PDR_AtaeRVFV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.6), data = data.PDR.AtaeRVFV, ylab = "PDR for RVFV in Ae. taeniorhynchus", xlab = "Temperature", pch = 1)
lines(PDR.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Get optimum for PDR: 37.5 C
Temp.xs[which.max(as.vector(PDR.AtaeRVFV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "50%"]))]


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

CpipWNV.PDR.avg <- CalcTraitAvgRecip(data.PDR.CpipWNV)
CtarWEEV.PDR.avg <- CalcTraitAvgRecip(data.PDR.CtarWEEV)


##########
###### 6. Plot preliminary figures
##########

plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.22), data = data.PDR, ylab = "PDR", xlab = "Temperature")
points(1/trait ~ T, data = data.PDR.CpipWNV, col = "red")
points(1/trait ~ T, data = data.PDR.CtarWNV, col = "dodgerblue")
points(1/trait ~ T, data = data.PDR.CtarWEEV, col = "blue")
points(1/trait ~ T, data = data.PDR.CtarSLEV, col = "darkgreen")
points(1/trait ~ T, data = data.PDR.AtaeSINV, col = "purple")


lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red")
lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "red")

lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "dodgerblue")
lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "dodgerblue")
lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "dodgerblue")

lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "blue")
lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "blue")

lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "darkgreen")
lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "darkgreen")
lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "darkgreen")

lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "purple")
lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "purple")
lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, col = "purple")


par(mfrow = c(2,3), mar = c(3, 4.5, 2, 1), oma = c(2, 0, 0, 0))

##### PDR for WNV in Cx. pipiens 
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.32), data = data.PDR.CpipWNV, xaxt = "n", pch = 19,
     ylab = "rate (1/day)", xlab = "", main = expression(paste("WNV in ",italic(Cx.)," ",italic(pipiens)," ")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.CpipWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### PDR for WNV Cx. tarsalis
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.32), data = data.PDR.CtarWNV, xaxt = "n", pch = 19,
     ylab = "rate (1/day)", xlab = "", main = expression(paste("WNV in ",italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.CtarWNV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### PDR for WEEV in Cx. tarsalis
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.32), data = data.PDR.CtarWEEV, xaxt = "n", pch = 19,
     ylab = "rate (1/day)", xlab = "", main = expression(paste("WEEV in ",italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.CtarWEEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### PDR for SLEV Cx. tarsalis
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.32), data = data.PDR.CtarSLEV, xaxt = "n", pch = 19,
     ylab = "rate (1/day)", xlab = "", main = expression(paste("SLEV in ",italic(Cx.)," ",italic(tarsalis))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.CtarSLEV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

##### PDR for SINV Ae. taenorhynchus
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.32), data = data.PDR.AtaeSINV, xaxt = "n", pch = 19,
     ylab = "rate (1/day)", xlab = "", main = expression(paste("SINV in ",italic(Ae.)," ",italic(taenorhynchus))), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.AtaeSINV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
legend("topleft", legend = "A", bty= "n", cex = 1.6, adj = c(1.5, 0))

mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 0.9)