## Marta Shocket, Stanford University
## Updated April 2018
##
## Purpose: Use Bayesian Inference (JAGS) to fit temperature-dependent functions for mosquito and pathogen traits for RRV model
##            with informative priors fit from other mosquito data
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) JAGS models
##           3) Shared settings for all models
##           4) Fit Cx. annulirostris trait thermal responses with data-informed priors
##           5) Plot Manuscript Figure 2


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("")

# Load libraties for fitting traits
library('R2jags')
library('mcmcplots')

# Load Data
data.RRV <- read.csv("RRVTraitData.csv") # Data from database for most traits (except below)
data.surv.proc <- read.csv("McDonaldSurvivalDataExpanded_ForJAGS.csv") # Processed survival data from McDonald 1980
data.EFD.proc <- read.csv("McDonaldEFDDataExpanded_ForJAGS.csv") # Processed fecundity data from McDonald 1980

#### EFD mean calcs for plotting
data.EFD.20 <- subset(data.EFD.proc, T == 20)
data.EFD.25 <- subset(data.EFD.proc, T == 25)
data.EFD.30 <- subset(data.EFD.proc, T == 30)
EFD.means <- data.frame(T = c(20, 25, 30), mean = c(mean(data.EFD.20$trait), mean(data.EFD.25$trait), mean(data.EFD.30$trait)), 
                        SE = c(sd(data.EFD.20$trait)/sqrt(nrow(data.EFD.20)), sd(data.EFD.25$trait)/sqrt(nrow(data.EFD.25)), sd(data.EFD.30$trait)/sqrt(nrow(data.EFD.30))))

#### Lifespan mean calcs for plotting
data.surv.20 <- subset(data.surv.proc, T == 20)
data.surv.25 <- subset(data.surv.proc, T == 25)
data.surv.30 <- subset(data.surv.proc, T == 30)
surv.means <- data.frame(T = c(20, 25, 30), mean = c(mean(data.surv.20$trait), mean(data.surv.25$trait), mean(data.surv.30$trait)), 
                         SE = c(sd(data.surv.20$trait)/sqrt(nrow(data.surv.20)), sd(data.surv.25$trait)/sqrt(nrow(data.surv.25)), sd(data.surv.30$trait)/sqrt(nrow(data.surv.30))))

# Load saved prior fits (gamma dists fit to parameter posteriors)
load("PriorFitsList.Rdata")
ls.prior.gamma.fits <- as.data.frame(prior.fits.list[[1]])
EFD.prior.gamma.fits  <- as.data.frame(prior.fits.list[[2]])
a.prior.gamma.fits  <- as.data.frame(prior.fits.list[3])
pEA.prior.gamma.fits  <- as.data.frame(prior.fits.list[4])
MDR.prior.gamma.fits  <- as.data.frame(prior.fits.list[5])
PDR.prior.gamma.fits  <- as.data.frame(prior.fits.list[6])


##########
###### 2. JAGS Models
##########



sink("quad_inf.txt")
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
    z.trait.mu.pred[i] <- -1 * cf.q * (Temp.xs[i] - cf.T0) * (Temp.xs[i] - cf.Tm) * (cf.Tm > Temp.xs[i]) * (cf.T0 < Temp.xs[i])
    }
    
    } # close model
    ",fill=T)
sink()


############## Briere Model with gamma priors

sink("briere_inf.txt")
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
###### 4. Fit Cx. Annulirostris / RRV Traits with data informed priors
##########

############## Lifespan for Cx. annulirostris - quadratic

##### Get data (reconstructed individual data from McDonald et al.) and priors 
data.ls.Cann <- data.surv.proc
data <- data.ls.Cann
hypers <- ls.prior.gamma.fits * 0.7

##### Organize data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
ls.Cann.out.inf<-jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                  n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
ls.Cann.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(ls.Cann.out.inf)

save(ls.Cann.out.inf, file = "jagsout_ls_Cann_inf.Rdata")

# Plot data + fit
plot(mean ~ T, xlim = c(7, 42), ylim = c(0,29), data = surv.means, pch = 19,
     ylab = "Lifespan (days)", xlab = expression(paste("Temperature (",degree,"C)")), main = "Adult Lifespan", cex.lab = 1.15)
arrows(x0 = surv.means$T, y0 = surv.means$mean + surv.means$SE, x1 = surv.means$T, y1 = surv.means$mean - surv.means$SE, length = 0, angle = 0)
lines(ls.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(ls.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(ls.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Caculate optimum for lifespan: 23.4 C
Temp.xs[which.max(as.vector(ls.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## EFD for Cx. annulirostris - Briere

##### Get data (processed/binned data from McDonald et al.) and priors
data <- data.EFD.proc
hypers <- EFD.prior.gamma.fits * 7.5

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
EFD.Cann.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
EFD.Cann.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(EFD.Cann.out.inf)

save(EFD.Cann.out.inf, file = "jagsout_EFD_Cann_inf.Rdata")

# Plot data + fit
plot(mean ~ T, xlim = c(5, 45), ylim = c(0,6), data = EFD.means, ylab = "EFD for Cx. annulirostris", xlab = "Temperature", pch = 19)
arrows(x0 = EFD.means$T, y0 = EFD.means$mean + EFD.means$SE, x1 = EFD.means$T, y1 = EFD.means$mean - EFD.means$SE, length = 0, angle = 0)
lines(EFD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(EFD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(EFD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Get optimum for EFD: 27.0 C
Temp.xs[which.max(as.vector(EFD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## a (= 1/GCD) for Cx. annulirostris - Briere

##### Get data and priors
data.GCD.Cann <- subset(data.RRV, trait.name == "GCD")
data <- data.GCD.Cann
hypers <- a.prior.gamma.fits*0.01
hypers[,3] <- a.prior.gamma.fits[,3]

##### Organize Data for JAGS
trait <- 1/data$trait # take the inverse since it's given in GCD and we want a
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
GCD.Cann.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
GCD.Cann.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(GCD.Cann.out.inf)

save(GCD.Cann.out.inf, file = "jagsout_GCD_Cann_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.GCD.Cann, ylab = "1/GCD = a for Cx. annulirostrus", xlab = "Temperature")
lines(GCD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(GCD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(GCD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Get optimum for GCD: 31.8 C
Temp.xs[which.max(as.vector(GCD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## pLA for Cx. annulirostris - quadratic

##### Get data and priors
data.pLA.Cann <- subset(data.RRV, trait.name == "pEA" & host.code == "Cann" & notes == "pLA")
data <- data.pLA.Cann
hypers <- pEA.prior.gamma.fits * 0.3

##### Organize Data for JAGS
trait <- data$trait
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
pLA.Cann.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
pLA.Cann.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(pLA.Cann.out.inf)

save(pLA.Cann.out.inf, file = "jagsout_pLA_Cann_inf.Rdata")

# Plot data + fit
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1), data = data.pLA.Cann, ylab = "pLA for Cx. annulirostrus", xlab = "Temperature")
lines(pLA.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(pLA.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Get optimum for pLA: 27.0 C
Temp.xs[which.max(as.vector(pLA.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## MDR for Cx. annulirostris - Briere

##### Get data and priors
data.MDR.Cann <- subset(data.RRV, trait.name == "1/MDR" & host.code == "Cann")
data <- data.MDR.Cann
hypers <- MDR.prior.gamma.fits

##### Organize Data for JAGS
trait <- 1/data$trait  # take the inverse since it's given in development time and we want rate
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
MDR.Cann.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
MDR.Cann.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(MDR.Cann.out.inf)

save(MDR.Cann.out.inf, file = "jagsout_MDR_Cann_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.2), data = data.MDR.Cann, ylab = "MDR for Cx. annulirostrus", xlab = "Temperature")
lines(MDR.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(MDR.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Get optimum for MDR: 32.6 C
Temp.xs[which.max(as.vector(MDR.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


############## PDR for RRV in O. vigilax - Briere

##### Get data and priors
data.PDR.RRV <- subset(data.RRV, trait.name == "EIP" & paras.code == "RRV")
data <- data.PDR.RRV
hypers <- PDR.prior.gamma.fits * 0.3

##### Organize Data for JAGS
trait <- 1/data$trait  # take the inverse since it's given in development time and we want rate
N.obs <- length(trait)
temp <- data$T

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, hypers = hypers)

##### Run JAGS
PDR.RRV.out.inf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_inf.txt",
                         n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

##### Examine Output
PDR.RRV.out.inf$BUGSoutput$summary[1:5,]
mcmcplot(PDR.RRV.out.inf)

save(PDR.RRV.out.inf, file = "jagsout_PDR_RRV_inf.Rdata")

# Plot data + fit
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.5), data = data.PDR.RRV, ylab = "PDR = 1/EIP for RRV in O. vigilax", xlab = "Temperature")
lines(PDR.RRV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(PDR.RRV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(PDR.RRV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Get optimum for PDR: 33 C
Temp.xs[which.max(as.vector(PDR.RRV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]))]


##########
###### 5. Plot Manuscript Figure 2: Traits with data-informed priors when available and uniform priors otherwise
##########

# Load/subset data, saved fits, and parameters for plotting

data.pLA.Cann <- subset(data.RRV, trait.name == "pEA" & host.code == "Cann" & notes == "pLA")
data.pRH.Cann <- subset(data.RRV, trait.name == "pEA" & host.code == "Cann" & notes == "pRH")
data.nLR.Cann <- subset(data.RRV, trait.name == "pEA" & host.code == "Cann" & notes == "nLR")
data.MDR.Cann <- subset(data.RRV, trait.name == "1/MDR" & host.code == "Cann")
data.EFD.Cann <- subset(data.RRV, trait.name == "EFD" & host.code == "Cann")
data.ls.Cann <- subset(data.RRV, trait.name == "1/mu" & host.code == "Cann")
data.GCD.Cann <- subset(data.RRV, trait.name == "GCD" & host.code == "Cann")
data.PDR.RRV <- subset(data.RRV, trait.name == "EIP" & paras.code == "RRV")
data.bc.RRV <- subset(data.RRV, trait.name == "bc" & paras.code == "RRV")

load("jagsout_bc_RRV_Ovig.Rdata")
load("jagsout_EFD_Cann_inf.Rdata")
load("jagsout_GCD_Cann_inf.Rdata")
load("jagsout_ls_Cann_inf.Rdata")
load("jagsout_MDR_Cann_inf.Rdata")
load("jagsout_nLR_Cann.Rdata")
load("jagsout_PDR_RRV_inf.Rdata")
load("jagsout_pLA_Cann_inf.Rdata")
load("jagsout_pRH_Cann.Rdata")

Temp.xs <- seq(5, 45, 0.2)
N.Temp.xs <-length(Temp.xs)

#### EFD mean calcs for plotting
data.EFD.20 <- subset(data.EFD.proc, T == 20)
data.EFD.25 <- subset(data.EFD.proc, T == 25)
data.EFD.30 <- subset(data.EFD.proc, T == 30)
EFD.means <- data.frame(T = c(20, 25, 30), mean = c(mean(data.EFD.20$trait), mean(data.EFD.25$trait), mean(data.EFD.30$trait)), 
                        SE = c(sd(data.EFD.20$trait)/sqrt(nrow(data.EFD.20)), sd(data.EFD.25$trait)/sqrt(nrow(data.EFD.25)), sd(data.EFD.30$trait)/sqrt(nrow(data.EFD.30))))

#### Lifespan mean calcs for plotting
data.surv.20 <- subset(data.surv.proc, T == 20)
data.surv.25 <- subset(data.surv.proc, T == 25)
data.surv.30 <- subset(data.surv.proc, T == 30)
surv.means <- data.frame(T = c(20, 25, 30), mean = c(mean(data.surv.20$trait), mean(data.surv.25$trait), mean(data.surv.30$trait)), 
                         SE = c(sd(data.surv.20$trait)/sqrt(nrow(data.surv.20)), sd(data.surv.25$trait)/sqrt(nrow(data.surv.25)), sd(data.surv.30$trait)/sqrt(nrow(data.surv.30))))

################## The Figure

par(mfrow = c(3,3), mar = c(3, 4.5, 2, 1), oma = c(2, 0, 0, 0))

##### biting rate
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.35), data = data.GCD.Cann, xaxt = "n", pch = 19,
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Biting Rate (",italic(a),")")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(GCD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(GCD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(GCD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
points(1/trait ~ T, data = data.GCD.Cann, pch = 19)
legend("topleft", legend = "A", bty= "n", cex = 1.3, adj = c(1, 0.5))

##### Transmission rate for RRV
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1.2), data = data.bc.RRV, xaxt = "n", pch = 19,
     ylab = "Transmission probability", xlab = "", main = expression(paste("Vector Competence (",italic(bc),")")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(bc.RRV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.RRV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(bc.RRV.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
points(trait ~ T, data = data.bc.RRV, pch = 19)
legend("topleft", legend = "B", bty= "n", cex = 1.3, adj = c(1, 0.5))

##### Lifespan / Survival
plot(mean ~ T, xlim = c(5, 45), ylim = c(0,29), data = surv.means, xaxt = "n", pch = 19,
     ylab = "Lifespan (days)", xlab = "", main = expression(paste("Adult Lifespan (",italic(lf)," = ",italic(mu)^-1,")")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(ls.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(ls.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(ls.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
arrows(x0 = surv.means$T, y0 = surv.means$mean + surv.means$SE, x1 = surv.means$T, y1 = surv.means$mean - surv.means$SE, length = 0, angle = 0)
points(mean ~ T, data = surv.means, pch = 19)
legend("topleft", legend = "C", bty= "n", cex = 1.3, adj = c(1, 0.5))

##### EIP for RRV
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.45), data = data.PDR.RRV, xaxt = "n", pch = 19,
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Parasite Development Rate (",italic(PDR),")")), cex.lab = 1.15, lwd = 1.5)
axis(1, at = seq(5, 45, 5))
lines(PDR.RRV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.RRV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(PDR.RRV.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)
points(1/trait ~ T, data = data.PDR.RRV, pch = 19)
legend("topleft", legend = "D", bty= "n", cex = 1.3, adj = c(1, 0.5))

##### Eggs per female per day
plot(mean ~ T, xlim = c(5, 45), ylim = c(0,5), data = EFD.means,  xaxt = "n", pch = 19,
     ylab = "Eggs per female per day", xlab = "", main = expression(paste("Fecundity (",italic(EFD),")")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(EFD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(EFD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(EFD.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
arrows(x0 = EFD.means$T, y0 = EFD.means$mean + EFD.means$SE, x1 = EFD.means$T, y1 = EFD.means$mean - EFD.means$SE, length = 0, angle = 0)
points(mean ~ T, data = EFD.means, pch = 19)
legend("topleft", legend = "E", bty= "n", cex = 1.3, adj = c(1, 0.5))

##### % of rafts that hatch
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1.1), data = data.pRH.Cann,  xaxt = "n", pch = 19,
     ylab = "Hatching probability", xlab = "", main = expression(paste("Egg Raft Viability (",italic(pRH),")")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(pRH.Cann.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pRH.Cann.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pRH.Cann.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
points(trait ~ T, data = data.pRH.Cann, pch = 19)
legend("topleft", legend = "F", bty= "n", cex = 1.3, adj = c(1, 0.5))

##### # of larvae per raft
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,300), data = data.nLR.Cann, xaxt = "n", pch = 19,
     ylab = "Emerged larvae", xlab = "", main = expression(paste("No. Viable Eggs per Raft (",italic(nLR),")")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(nLR.Cann.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"]*1000 ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(nLR.Cann.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"]*1000 ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(nLR.Cann.out$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"]*1000 ~ Temp.xs, lwd = 1.5)
points(trait ~ T, data = data.nLR.Cann, pch = 19)
legend("topleft", legend = "G", bty= "n", cex = 1.3, adj = c(1, 0.5))
mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 0.9)

##### Larval to adult survival
plot(trait ~ T, xlim = c(5, 45), ylim = c(0,1.1), data = data.pLA.Cann, xaxt = "n", pch = 19,
     ylab = "Survival probability", xlab = "", main = expression(paste("Larval-to-Adult Survival (",italic(pLA),")")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(pLA.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(pLA.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
points(trait ~ T, data = data.pLA.Cann, pch = 19)
legend("topleft", legend = "H", bty= "n", cex = 1.3, adj = c(1, 0.5))
mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 0.9)

##### Mosquito dev. rate
plot(1/trait ~ T, xlim = c(5, 45), ylim = c(0,0.18), data = data.MDR.Cann, xaxt = "n", pch = 19,
     ylab = "Rate (1/day)", xlab = "", main = expression(paste("Mosquito Development Rate (",italic(MDR),")")), cex.lab = 1.15)
axis(1, at = seq(5, 45, 5))
lines(MDR.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2, col = "red", lwd = 1.5)
lines(MDR.Cann.out.inf$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs, lwd = 1.5)
points(1/trait ~ T, data = data.MDR.Cann, pch = 19)
legend("topleft", legend = "I", bty= "n", cex = 1.3, adj = c(1, 0.5))
mtext(expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, las = 1, cex = 0.9)
