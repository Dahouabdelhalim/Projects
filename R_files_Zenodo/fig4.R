rm(list = ls())
source("R/mechanistic_mods.R")
source("R/likelihoods.R")
source("R/fit_model.R")
source("R/helpers.R")
library(knitr)

png(file="figures/fig4.png", width=700, height=500, units="px", pointsize=12)
par(mfrow=c(2, 2))

source('R/karvonen2.R')
xvals <- seq(min(d$Cerc), max(d$Cerc), length.out=100)
plot(jitter(d$Cerc, .2), d$Metacerc, 
     xlab="Number of cercariae", ylab="Number of metacercariae")
lines(xvals, 
      y=powerC(model_fits[['Power C']]$par, 
               Np=xvals, 
               H=unique(d$host), t = unique(d$time)),
      lwd=2)
legend('topleft', col=1, lty=1, lwd=2, legend=c("Power C"), bty='n')

source('R/karvonen1.R')
xvals <- seq(min(d$Cerc), max(d$Cerc), length.out=200)
plot(jitter(d$Cerc, .2), 
     d$Meta, 
     xlab="Number of cercariae", ylab="Number of metacercariae")
best_models
lines(xvals, 
      y = NB1(model_fits[['Negative binomial 1']]$par, 
              Np = xvals, t = unique(d$time)), 
      col=1, lty=1, lwd=2)
lines(xvals, 
      y = NB2(model_fits[['Negative binomial 2']]$par, 
              Np=xvals, H=unique(d$host), t=unique(d$time)), 
      col=2, lty=2, lwd=2)
lines(xvals, 
      y=powerC(model_fits[['Power C']]$par, 
               Np=xvals, 
               H=unique(d$host), t = unique(d$time)),
      col=3, lty=3, lwd=2)
legend('topleft', col=1:3, 
       lty=1:3, lwd=2, 
       legend=c("Neg. Binom. 1", "Neg. Binom. 2", "Power C"), 
       bty='n')

source('R/paller2.R')
xvals <- seq(min(d$Cerc), max(d$Cerc), length.out=200)
plot(jitter(d$Cerc, .2), d$Meta, 
     xlab="Number of cercariae", ylab="Number of metacercariae")
best_models
lines(xvals, 
      y = NB1(model_fits[['Negative binomial 1']]$par, 
              Np = xvals, t = unique(d$time)), 
      col=1, lty=1, lwd=2)
lines(xvals, 
      y = NB2(model_fits[['Negative binomial 2']]$par, 
              Np=xvals, H=unique(d$host), t=unique(d$time)), 
      col=2, lty=2, lwd=2)
legend("topleft", lty=c(1, 2), col=1:2, lwd=c(2, 2), 
       legend=c("Neg. Binom. 1", "Neg. Binom. 2"), bty='n')

source('R/paller1.R')
xvals <- seq(min(d$Cerc), max(d$Cerc), length.out=200)
plot(jitter(d$Cerc, .2), d$Meta, 
     xlab="Number of cercariae", ylab="Number of metacercariae")
best_models
lines(xvals, 
      y = NB1(model_fits[['Negative binomial 1']]$par, 
              Np = xvals, t = unique(d$time)), 
      col=1, lty=1, lwd=2)
lines(xvals, 
      y = NB2(model_fits[['Negative binomial 2']]$par, 
              Np=xvals, H=unique(d$host), t=unique(d$time)), 
      col=2, lty=2, lwd=2)
legend("topleft", lty=c(1, 2), col=1:2, lwd=c(2, 2), 
       legend=c("Neg. Binom. 1", "Neg. Binom. 2"), bty='n')

dev.off()