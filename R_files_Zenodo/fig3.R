rm(list = ls())
source("R/mechanistic_mods.R")
source("R/likelihoods.R")
source("R/fit_model.R")
source("R/helpers.R")
library(knitr)
source('R/behavior_volume2.R')

png(file = "figures/fig3.png", width = 450, height = 300, units = "px", 
    pointsize = 12)
xvals <- seq(min(d$Number), max(d$Number), length.out=200)
plot(jitter(d$Number, .2), d$Meta, 
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
