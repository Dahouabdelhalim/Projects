rm(list = ls())
source("R/mechanistic_mods.R")
source("R/likelihoods.R")
source("R/fit_model.R")
source("R/helpers.R")
library(knitr)

png(file="figures/fig2.png", width=700, height=500, units="px", pointsize=12)
par(mfrow=c(2, 2))
source('R/volume.R')
plot(jitter(d$Number), d$Meta, 
     xlab="Number of Cercariae", ylab="Number of metacercariae")
xvals <- seq(min(d$Number), max(d$Number), length.out=100)
lines(xvals, 
      y = NB1(model_fits[['Negative binomial 1']]$par, 
            Np = xvals, t = unique(d$time)), 
      col = 1, lty = 1, lwd = 2)
lines(xvals, 
      y = NB2(model_fits[['Negative binomial 2']]$par, 
              Np=xvals, H=unique(d$host), t=unique(d$time)), 
      col=2, lty=2, lwd=2)
legend("topleft", legend=c("Neg. Binom. 1", "Neg. Binom. 2"), 
       col=1:2, 
       lty=1:2,
       lwd=2,
       bty="n")
mtext("A", at=0)

source('R/host_number.R')
xvals <- seq(min(d$host), max(d$host), length.out=200)
plot(jitter(d$host, .2), d$Meta, 
     xlab="Number of Hosts", ylab="Number of metacercariae",
     xaxt="n")
axis(1, at=1:max(d$host))
lines(xvals, 
      y=powerC(model_fits[['Power C']]$par, 
               Np=unique(d$Number), 
               H=xvals, t = unique(d$time)), 
      col=1, lty=1, lwd=2)
lines(xvals, 
      y=powerH(model_fits[['Power H']]$par, 
                Np=unique(d$Number), 
                H=xvals, t = unique(d$time)), 
      lty=2, lwd=2, col=2)
legend("topleft", legend=c("Power C","Power H"), col=1:2, 
       lty=1:2, lwd=rep(2, 2), bty="n")
mtext("B", at=0)

source('R/time.R')
xvals <- seq(min(d$time), max(d$time), length.out=200)
plot(jitter(d$time), d$Meta, xlab="Exposure duration (min)", 
     ylab="Number of metacercariae", 
     ylim=c(0, max(d$Number) + 3))
lines(xvals, 
      y=powerC(model_fits[['Power C']]$par, 
               Np=unique(d$Number), 
               H=unique(d$host), t = xvals),
      lwd=2)
legend("topleft", legend=c("Power C"), 
       col=1, lty=1, lwd=2, bty="n")
mtext("C", at=0)

source('R/density.R')
plot_best(model_fits, best_models, d)
mtext("D", at=0)
dev.off()
