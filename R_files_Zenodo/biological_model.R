###########################################################################
# model
###########################################################################
require(jagsUI)
sink(file="ms_rd.jags")
cat("
    model {
    
    ######################################################################################
    # survival and breeding probability
    ######################################################################################
    alpha.phi ~ dnorm(2, 0.01)
    beta.phi ~ dnorm(0, 0.01)
    
    mean.gamma ~ dbeta(72,28)
    mu.gamma <- log(mean.gamma/(1 - mean.gamma))
    
    beta.gamma ~ dnorm(0, 0.01)
    beta.tar.phi ~ dnorm(0, 0.01)
    
    sigma.phi ~ dunif(0,5)
    sigma.gamma ~ dunif(0,5)
    tau.phi <- pow(sigma.phi, -2)
    tau.gamma <- pow(sigma.gamma, -2)
    
    
    for (t in 1:(n.years-1)){
    mu.phi[t] <- alpha.phi + beta.phi * t
    eps.phi[t] ~ dnorm(mu.phi[t], tau.phi)
    eps.gamma[t] ~ dnorm(0, tau.gamma)
    }
    
    for (i in 1:n.ind){
    for (t in 1:(n.years-1)){
    logit(phi[i,t]) <- eps.phi[t] + (beta.tar.phi * ltar[i])
    logit(gamma[i,t]) <- mu.gamma + (beta.gamma * ltar[i]) + eps.gamma[t]
    }
    }
    
    
    ######################################################################################
    # Secondary occasions p's
    ######################################################################################
    for (j in 1:3){
    mean.p[j] ~ dnorm(0, 0.01)
    sig.p[j] ~ dunif(0,5)
    tau.p[j] <- pow(sig.p[j], -2)
    }
    
    
    for (t in 1:n.years){
    for (j in 1:n.sec[t]){
    
    eps.p[t,j] ~ dnorm(0, tau.p[j])
    
    logit(p[t,j]) <- mean.p[j] + eps.p[t,j]
    
    yes[t,j] ~ dbin(p[t,j], total[t,j])
    
    }
    }   
    
    
    ######################################################################################
    # Primary occasions p's
    ######################################################################################    
    for (t in 1:n.years){
    pstar[t] <- 1 - prod(1 - p[t,])
    }
    
    
    ###########################################################################
    # likelihood
    ###########################################################################
    for (i in 1:n.ind){
    
    z[i,first[i]] <- ch[i,first[i]]
    
    
    for (t in (first[i]+1):n.years){
    
    mu1[i,t] <- z[i,t-1] * phi[i,t-1]
    mu2[i,t] <- z[i,t] * gamma[i,t-1] * pstar[t]
    
    z[i,t] ~ dbern(mu1[i,t])
    ch[i,t] ~ dbern(mu2[i,t])
    
    } 
    
    }
    
    ####################################################################################
    # derived parameters
    ####################################################################################
    for (t in 1:(n.years - 1)){
    logit(PHI[t]) <- eps.phi[t]
    logit(GAMMA[t]) <- mu.gamma + eps.gamma[t]
    }
    
    for (k in 1:100){
    logit(GAMMA.SIZE[k]) <- mu.gamma + beta.gamma * size[k]
    logit(PHI.SIZE[k]) <- alpha.phi + (beta.phi * 13) + beta.tar.phi * size[k]
    }
    
    ####################################################################################
    # end model
    ####################################################################################
    
    }
    
    ", fill=TRUE)
sink()


######################################################################################
# add data
######################################################################################
n.sec.occ <- 3
n.years <- 27

min(zltar)
max(zltar)

size <- seq(min(zltar), max(zltar), length.out = 100)

######################################################################################
# bundle data
######################################################################################
dat <- list(first = first, ch = tut, ltar = zltar,
            n.sec = rep(n.sec.occ,n.years), n.years = ncol(tut), n.ind = nrow(tut),
            yes = yes, total = total, size = size)







######################################################################################
# Initial Values
######################################################################################
inits <- function(){list(z = z.init)}  


######################################################################################
# parameters to monitor
######################################################################################
pars <- c('pstar','mu.phi','mu.gamma','eps.phi','eps.gamma',
          'beta.phi','beta.gamma','p','sig.p','sigma.phi','sigma.gamma',
          'PHI','GAMMA','GAMMA.SIZE', 'alpha.phi','beta.tar.phi','PHI.SIZE')

str(msrd$mean)
######################################################################################
# MCMC settings
######################################################################################
n.chains <- 1
n.thin <- 5
n.adapt <- 100
n.iter <- 5000
n.burnin <- 2500


######################################################################################
# compile model! 
######################################################################################
msrd <- jags(dat, inits, pars, "ms_rd.jags", 
             n.chains = n.chains, n.adapt = n.adapt, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)

summary(msrd)
boxplot(msrd$sims.list$beta.tar.phi)
str(msrd)
msrd$f$beta.tar.phi
msrd$mean$beta.tar.phi
msrd$mean$beta.phi

msrd <- update(msrd, n.iter = 5000)
summary(msrd$mean$beta.tar.phi)
summary(msrd$f$beta.tar.phi)

summary(msrd$mean$beta.phi)
summary(msrd$f$beta.phi)

summary(msrd$mean$beta.tar.phi)
summary(msrd$f$beta.tar.phi)
######################################################################################
# PLOTS!
######################################################################################
par(family = 'serif', mfrow = c(1,1))

# hist(ltar)

######################################################################################
# detection
######################################################################################

# par(family = 'serif', mfrow = c(1,1))

Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.23/bin/gswin64c.exe")
bitmap("W:\\\\Ph_D\\\\Manuscripts\\\\robust_design\\\\MEE_submission\\\\Resubmission\\\\in_press_version\\\\Figure4.tiff",
       width = 12, height = 12, units = 'in', res = 300, type = 'tifflzw')
par(mfrow = c(2,2), mar = c(5.1,5.1,2.1,2.1))
plot(msrd$mean$pstar ~ seq(1988,2014,1), ylab = expression('p*'), xlab = 'Year',
     cex.lab = 2, pch = 19, ylim = c(0.2,0.8))
arrows(seq(1988,2014,1), msrd$q2.5$pstar, seq(1988,2014,1), msrd$q97.5$pstar, code = 3, length = 0, angle = 90, lty = 2)
plot(msrd$mean$p[,1] ~ seq(1988,2014,1), ylab = expression(p['nest']), xlab = 'Year',
     cex.lab = 2, pch = 19, ylim = c(0.2,0.8), las = 1)
arrows(seq(1988,2014,1), msrd$q2.5$p[,1], seq(1988,2014,1), msrd$q97.5$p[,1], code = 3, length = 0, angle = 90, lty = 2)
plot(msrd$mean$p[,2] ~ seq(1988,2014,1), ylab = expression(p['tower']), xlab = 'Year',
     cex.lab = 2, pch = 19, ylim = c(0,0.5), las = 1)
arrows(seq(1988,2014,1), msrd$q2.5$p[,2], seq(1988,2014,1), msrd$q97.5$p[,2], code = 3, length = 0, angle = 90, lty = 2)
plot(msrd$mean$p[,3] ~ seq(1988,2014,1), ylab = expression(p['banding']), xlab = 'Year',
     cex.lab = 2, pch = 19, ylim = c(0,0.5), las = 1)
arrows(seq(1988,2014,1), msrd$q2.5$p[,3], seq(1988,2014,1), msrd$q97.5$p[,3], code = 3, length = 0, angle = 90, lty = 2)
dev.off()

# par(family = 'serif', mfrow = c(1,1))
# plot(msrd$mean$pstar ~ seq(1988,2014,1), ylab = expression('p*'), xlab = 'Year',
#      cex.lab = 2, pch = 19, ylim = c(0.2,0.8))
# arrows(seq(1988,2014,1), msrd$q2.5$pstar, seq(1988,2014,1), msrd$q97.5$pstar, code = 3, length = 0, angle = 90, lty = 2)



######################################################################################
# phenotype simple
######################################################################################
par(family = 'sans', mfrow = c(1,1))
######################################################################################
str(msrd)
x <- size
x <- (x * 37.0001 + 655.1193)/10
print(x)

msrd$mean$GAMMA.SIZE
Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.23/bin/gswin64c.exe")
bitmap("W:\\\\Ph_D\\\\Manuscripts\\\\robust_design\\\\MEE_submission\\\\Resubmission\\\\in_press_version\\\\Figure6.tiff",
       width = 12, height = 12, units = 'in', res = 300, type = 'tifflzw')
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(msrd$mean$GAMMA.SIZE ~ x, ylim = c(0.4,1), xlim = c(50, 79.7), las = 1,
     xlab = 'Tarsus at ~30 days (mm)', ylab = expression(gamma[italic(i)]), cex.lab = 2, pch = 19, type = 'n')
lines(msrd$mean$GAMMA.SIZE ~ x, lwd = 2)
lines(msrd$q2.5$GAMMA.SIZE ~ x, lwd = 2, lty = 2)
lines(msrd$q97.5$GAMMA.SIZE ~ x, lwd = 2, lty = 2)
dev.off()







g <- c(msrd$mean$GAMMA.SMALL, msrd$mean$GAMMA.TYPICAL, msrd$mean$GAMMA.LARGE)
gl <- c(msrd$q2.5$GAMMA.SMALL, msrd$q2.5$GAMMA.TYPICAL, msrd$q2.5$GAMMA.LARGE)
gu <- c(msrd$q97.5$GAMMA.SMALL, msrd$q97.5$GAMMA.TYPICAL, msrd$q97.5$GAMMA.LARGE)
x <- c(1.25, 2, 2.75)

hist(ltar/10, main = NULL, xlab = 'Tarsus as a Gosling (mm)', col = 'red', cex.lab = 2)

plot(g ~ x, ylim = c(0.6,0.9), xlim = c(1, 3), xaxt = 'n', las = 1,
     xlab = 'Size at ~30 days', ylab = expression(gamma), cex.lab = 2, pch = 19)
arrows(x, gl, x, gu, code = 3, length = 0.02, angle = 90)
axis(1, at = c(1.25,2,2.75), labels = c('Small','Typical','Large'), cex.axis = 1.25)



######################################################################################
# temporal variation
######################################################################################
par(family = 'sans', mfrow = c(1,1))
str(msrd)

bitmap("W:\\\\Ph_D\\\\Manuscripts\\\\robust_design\\\\MEE_submission\\\\Resubmission\\\\in_press_version\\\\Figure7.tiff",
       width = 12, height = 12, units = 'in', res = 300, type = 'tifflzw')
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(msrd$mean$PHI ~ seq(1988,2013,1), ylab = expression(phi), xlab = 'Year',
     cex.lab = 2, pch = 19, ylim = c(0.65,0.925))
arrows(seq(1988,2013,1), msrd$q2.5$PHI, seq(1988,2013,1), msrd$q97.5$PHI, code = 3, length = 0, angle = 90, lty = 2)
lines(plogis(msrd$mean$mu.phi) ~ seq(1988,2013,1), lwd = 2)
lines(plogis(msrd$q2.5$mu.phi) ~ seq(1988,2013,1), lty = 2, lwd = 2)
lines(plogis(msrd$q97.5$mu.phi) ~ seq(1988,2013,1), lty = 2, lwd = 2)
dev.off()




bitmap("W:\\\\Ph_D\\\\Manuscripts\\\\robust_design\\\\MEE_submission\\\\Resubmission\\\\in_press_version\\\\Figure5.tiff",
       width = 12, height = 12, units = 'in', res = 300, type = 'tifflzw')
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(msrd$mean$GAMMA ~ seq(1989,2014,1), ylab = expression(gamma), xlab = 'Year',
     cex.lab = 2, pch = 19, ylim = c(0.2,1))
arrows(seq(1989,2014,1), msrd$q2.5$GAMMA, seq(1989,2014,1), msrd$q97.5$GAMMA, code = 3, length = 0, angle = 90, lty = 2)
dev.off()