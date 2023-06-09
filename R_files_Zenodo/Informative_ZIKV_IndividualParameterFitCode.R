### This file will process the temperature data and create the individual parameter 
### samples for all the traits that are included in the DENV R0 model. 

wd = "/Users/erinmordecai/Documents/projects/Stanford projects/temperature and vector-borne disease/R0 for Ae aegypti and Ae albopictus/Final Code Package"

## Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')
library('MASS')

# This file contains tools for analysis and visualization.
source(paste(wd, "mcmc_utils_all.R", sep = "/")) 

# This file contains the thermal response functions
source(paste(wd, "temp_functions_all.R", sep = "/")) 

# Load posterior distributions of traits previously fit for Aedes aegypti and DENV
load(paste(wd, "LifespanFits_2016-03-30-informative.Rsave", sep = "/"))

load(paste(wd, "Informative_Aegypti_DENV_ParameterFits_2016-03-30.Rsave", sep = "/"))

# Fit gamma distributions of priors for each previous trait from which we want to derive priors

gamma.fits.PDR = apply(PDR.samps, 2, function(df) fitdistr(df, "gamma")$estimate)

lf.samps.pos = lf.DENV.samps[,c(1:4,6)]
lf.samps.pos$qd <- -lf.samps.pos$qd
gamma.fits.lf = apply(lf.samps.pos, 2, function(df) fitdistr(df, "gamma")$estimate)

# the previous work fit b and c separately
# Since c has a broader distribution, use it for a prior
gamma.fits.c = apply(c.samps, 2, function(df) fitdistr(df, "gamma")$estimate)

# Creating a small constant to keep denominators from being zero.
ec<-0.000001 

## Loading the data
data.all <- read.csv("zikv_traits.csv", header=TRUE)
head(data.all)
str(data.all)

## Now the code will choose all the temperature sensitive traits and fit either a
## Briere or Quadratic model to it using MCMC sampling.
# Specifing the parameters that control the MCMC (these will be used throughout the code). 

n.chains <- 5
n.adapt <- 5000
n.samps <- 5000

###########################################
## The next trait we'll be choosing is the vector competence, b*c.
## Due to the data that was collected it will be necessary
## to decompose it into its two parts, b and c.

# Choose b, the probability a human will be bitten and infected by
# an infectious mosquito (ie. transmission).

data <- data.all[ which(data.all$trait.name=="bc"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

# Since the data are somewhat ambiguous, and previous work fit Briere functions for this trait
# Fit an informative Briere here

hypers = gamma.fits.c

jags <- jags.model(paste(wd, 'jags-briere-informative.bug', sep = "/"),
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T), 'hypers' = hypers),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
b.samps <- samps

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
priors1<-list()
priors1$names<-c( "T0", "Tm", "c","tau")
priors1$fun<-c( "gamma", "gamma", "gamma","gamma")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(hypers[1,1], hypers[1,2], NA)
priors1$hyper[,2]<-c(hypers[1,2], hypers[2,2], NA)
priors1$hyper[,3]<-c(hypers[1,3], hypers[2,3], NA)
priors1$hyper[,4]<-c(hypers[1,5], hypers[2,5], NA)
plot.hists(b.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", b.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(15, 40), ylim = c(0,0.4),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

# The Briere function really doesn't fit well for these data
# Stick with the uninformative Quadratic fit

#########################################
## Fit lifespan
## Lifespan estimated by calculating Kaplan-Meyer survival rates from censored data
## Then fitting Gompertz curves to the survival rates
## Then estimating area under the curve
## For 38 C survival was very low and Gompertz curves didn't fit well, 
## so area under curve estimated directly from data

data = subset(data.all, trait.name=="lf")

plot(trait ~ T, data = data, xlim = c(10, 40))
points(trait ~ T, data = subset(data, rep %in% c("1-inf", "2-inf")), pch = 16, col = 2)

hypers = gamma.fits.lf

# Note that estimated lifespans are very long for infected mosquitoes at near-optimal temperatures

# Given the data we've chosen to use the Quadratic function. 

jags <- jags.model(paste(wd, 'jags-quad-neg-informative.bug', sep = "/"),
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T), 'hypers' = hypers),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, n.qd=0.005),
                   n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('Tm', 'T0', 'qd', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.quad.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
lf.samps <- samps

priors1<-list()
priors1$names<-c( "T0", "Tm", "qd","tau")
priors1$fun<-c( "gamma", "gamma", "gamma","gamma")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(hypers[1,1], hypers[2,1], NA)
priors1$hyper[,2]<-c(hypers[1,2], hypers[2,2], NA)
priors1$hyper[,3]<-c(hypers[1,3], hypers[2,3], NA)
priors1$hyper[,4]<-c(hypers[1,4], hypers[2,4], NA)

## the plot.hists command can be used to plot histograms of posterior
## samples of parameters with overlying priors
plot.hists(cbind(lf.samps[,1:2], -lf.samps[,3], lf.samps[,4]), my.par=c(2,2), n.hists=4, priors=priors1)


# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="quad", lf.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(0, 40), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)


#########################################
## The final parameter to fit is PDR (EIR), the inverse of the extrinsic incubation period
## Calculated as the time it takes for the population to get to half of its maximum proportion infected

data <- subset(data.all, trait.name=="EIR")

# Plot the data to see which function, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

# Tm is very uncertain because of the strongly bimodal data
# Add a prior that is stronger for Tm than for other parameters
hypers = gamma.fits.PDR*0.3
hypers[,2] = hypers[,2]*3

# Given the data the Briere fuction is chosen. Jag-briere.bug contains the specifics of 
# the Briere model with the default priors.

jags <- jags.model(paste(wd, 'jags-briere-informative.bug', sep = "/"),
                   data = list('Y' = data$trait, 'T' = data$T, 'N' = length(data$T), 'hypers' = hypers),
                   n.chains = n.chains, inits = list(Tm = 38, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt) 

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c','Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps <- make.briere.samps(coda.samps, nchains=n.chains, samp.lims=c(1, n.samps))
samps$tau <- 1/samps$sigma
PDR.samps <-  samps

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
priors1<-list()
priors1$names<-c( "T0", "Tm", "c","tau")
priors1$fun<-c( "gamma", "gamma", "gamma","gamma")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(hypers[1,1], hypers[1,2], NA)
priors1$hyper[,2]<-c(hypers[1,2], hypers[2,2], NA)
priors1$hyper[,3]<-c(hypers[1,3], hypers[2,3], NA)
priors1$hyper[,4]<-c(hypers[1,5], hypers[2,5], NA)

plot.hists(PDR.samps[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)

# Plot the fits against the data
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="briere", PDR.samps, Temps, thinned=seq(1,n.samps, length=1000))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim = c(5, 45), 
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab=data$trait.name[1],
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)


## This code is just save the MCMC samples for further analysis in the R0 model.
save(b.samps, lf.samps, PDR.samps,
     file = "zikv_trait_fits_informative.Rsave")

