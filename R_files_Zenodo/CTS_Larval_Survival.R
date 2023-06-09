###California Tiger Salamander (Ambystoma californiense) larval survival given egg density
###Arianne Messerman
###February 2021

#setwd("")

#Data already natural log-transformed
lar<-read.csv("Larval-Survival-Density.csv")
lar<-lar[order(lar$ln.dens, decreasing = FALSE),]
lar$dens.cent<-lar$ln.dens-mean(lar$ln.dens)
str(lar)

plot(lar$ln.dens, lar$ln.surv)

library(jagsUI)
library(mcmcplots)
set.seed(423467)

##########################################################
## Linear regression of larval survival given egg density
## DIC=28.55
##########################################################

sink("lar.surv.LR.jags")
cat("
 model {
 
    #Priors
    alpha ~ dnorm(0, 0.001) #Intercept
    beta ~ dnorm(0, 0.001) #linear  slope of survival
    sigma ~ dunif(0, 100)
    
    #Likelihood
    for (i in 1:n){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta*x[i]
    }
    
    #Derived quantities
    tau <- 1/(sigma*sigma)
    
    #Assess model fit using sum-of-squares-type discrepancy
    for (i in 1:n){
    residual[i] <- y[i]-mu[i]     #Residuals for observed data
    predicted[i] <- mu[i]         #Predicted values
    sq[i] <- pow(residual[i],2)   #Squared residuals for observed data
    
    #Generate replicate data and compute fit stats for them
    y.new[i] ~ dnorm(mu[i], tau)  #One new data set at each MCMC iteration
    sq.new[i] <- pow(y.new[i]-predicted[i],2) #Squared residuals for new data
    }
    
    fit <- sum(sq[])              #Sum of squared residuals for actual data set
    fit.new <- sum(sq.new[])      #Sum of squared residuals for new data set
    test <- step(fit.new - fit)   #Test whether nes data set more extreme
    bpvalue <- mean(test)         #Bayesian p-value
 }
", fill=TRUE)
sink()

#Bundle data
jags.data <- list(y = c(lar$ln.surv), n = length(lar[,1]),  x = c(lar$dens.cent))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(1), sigma = runif(1,0,2))}

#Parameters to estimate
parameters <- c("alpha", "beta", "sigma", "fit", "fit.new", "mu",
                "bpvalue", "residual", "predicted")

#MCMC settings
ni <- 3000
nt <- 2
nb <- 1000
nc <- 3

# Call JAGS from R
lar.surv.LR <- jags(jags.data, inits, parallel=TRUE, parameters, "lar.surv.LR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(lar.surv.LR, digits = 3)#DIC=28.36

#Residual plot
plot(lar.surv.LR$mean$predicted, lar.surv.LR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(lar.surv.LR$sims.list$fit, lar.surv.LR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(lar.surv.LR$sims.list$fit.new>lar.surv.LR$sims.list$fit)#Bayesian p-value=0.57

#Examine posteriors
plot(lar.surv.LR)

##Larval density-dependent survival equation
pred.lar <- lar.surv.LR$mean$alpha + lar.surv.LR$mean$beta[1]*lar$dens.cent

library(dplyr)
larv <- dplyr::mutate(lar, mu = lar.surv.LR$mean$mu, #mu are the predicted values
                       q2.5_mu = lar.surv.LR$q2.5$mu, 
                       q97.5_mu = lar.surv.LR$q97.5$mu)

#Check that mu and function closely match
plot(lar$ln.dens, lar$ln.surv)
points(lar$ln.dens, pred.lar, type="l", lwd=2)
points(lar$ln.dens, larv$mu, type="l", col=2, lwd=2)

library(ggplot2)
library(cowplot)
#Figure showing 95% Credible intervals around mean estimate
ggplot(larv) + 
  geom_ribbon(aes(x = dens.cent, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
  geom_path(aes(x = dens.cent, y = mu), color = "red") +
  geom_point(aes(x = dens.cent, y = ln.surv)) +
  scale_y_continuous("Ln(larval survival)")+
  scale_x_continuous("Ln(egg density (eggs/m^3))")+
  theme_cowplot()

#Allow survival curves to sample from residual distribution around the function mean
plot(lar$ln.dens, lar.surv.LR$mean$residual)
resid.lar<-c(lar.surv.LR$mean$residual)
plot(density(resid.lar))

mean(lar$ln.dens)#3.38

##########################################################
##Quadratic regression of larval survival given egg density
#DIC=29.42
##########################################################
sink("lar.surv.QR.jags")
cat("
 model {
 
    #Priors
    alpha ~ dnorm(0, 0.001) #Intercept
    for (i in 1:2){
    beta[i] ~ dnorm(0, 0.001) #linear and quadratic slopes of survival
    }
    sigma ~ dunif(0, 100)
    
    #Likelihood
    for (i in 1:n){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta[1]*x[i] + beta[2]*pow(x[i],2)
    }
    
    #Derived quantities
    tau <- 1/(sigma*sigma)
    
    #Assess model fit using sum-of-squares-type discrepancy
    for (i in 1:n){
    residual[i] <- y[i]-mu[i]     #Residuals for observed data
    predicted[i] <- mu[i]         #Predicted values
    sq[i] <- pow(residual[i],2)   #Squared residuals for observed data
    
    #Generate replicate data and compute fit stats for them
    y.new[i] ~ dnorm(mu[i], tau)  #One new data set at each MCMC iteration
    sq.new[i] <- pow(y.new[i]-predicted[i],2) #Squared residuals for new data
    }
    
    fit <- sum(sq[])              #Sum of squared residuals for actual data set
    fit.new <- sum(sq.new[])      #Sum of squared residuals for new data set
    test <- step(fit.new - fit)   #Test whether nes data set more extreme
    bpvalue <- mean(test)         #Bayesian p-value
 }
", fill=TRUE)
sink()

#Bundle data
jags.data <- list(y = c(lar$ln.surv), n = length(lar[,1]),  x = c(lar$dens.cent))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(2, 0, 1), sigma = runif(1,0,2))}

#Parameters to estimate
parameters <- c("alpha", "beta", "sigma", "fit", "fit.new", "mu",
                "bpvalue", "residual", "predicted")

#MCMC settings
ni <- 3000
nt <- 2
nb <- 1000
nc <- 3

# Call JAGS from R
lar.surv.QR <- jags(jags.data, inits, parallel=TRUE, parameters, "lar.surv.QR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(lar.surv.QR, digits = 3)#DIC=28.96

#Residual plot
plot(lar.surv.QR$mean$predicted, lar.surv.QR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(lar.surv.QR$sims.list$fit, lar.surv.QR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(lar.surv.QR$sims.list$fit.new>lar.surv.QR$sims.list$fit)#Bayesian p-value=0.58

#Examine posteriors
plot(lar.surv.QR)

##lar density-dependent survival equation
pred.lar2 <- lar.surv.QR$mean$alpha + lar.surv.QR$mean$beta[1]*lar$ln.dens + lar.surv.QR$mean$beta[2]*(lar$ln.dens*lar$ln.dens)

library(dplyr)
lars2 <- dplyr::mutate(lar, mu = lar.surv.QR$mean$mu, #mu are the predicted values
                        q2.5_mu = lar.surv.QR$q2.5$mu, 
                        q97.5_mu = lar.surv.QR$q97.5$mu)


#Check that mu and function closely match
plot(lar$ln.dens, lar$ln.surv)
points(lar$ln.dens, pred.lar2, type="l", lwd=2)
points(lar$ln.dens, lars2$mu, type="l", col=2, lwd=2)

#Figure showing 95% Credible intervals around mean estimate
ggplot(lars2) + 
  geom_ribbon(aes(x = ln.dens, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
  geom_path(aes(x = ln.dens, y = mu), color = "red") +
  geom_point(aes(x = ln.dens, y = ln.surv)) +
  scale_y_continuous("Ln(larval survival)")+
  scale_x_continuous("Ln(egg density (eggs/m^3))")+
  theme_cowplot()

#Allow survival curves to sample from residual distribution around the function mean
plot(lar$ln.dens, lar.surv.QR$mean$residual)
resid.lar.qr<-c(lar.surv.QR$mean$residual)
plot(density(resid.lar))

#Pull 500 sets of iterations from posteriors of coefficients to .CSV
mcmc.sample<-lar.surv.QR$mcmc.info$n.samples#Same across life stage models
sub.set <- sort(sample(1:mcmc.sample, size=500))

lar.inter1 <- lar.slope1 <- lar.sigma1  <- c(1:length(sub.set))
lar.inter <- lar.l.slope <- lar.q.slope <- lar.sigma <-c(1:length(sub.set))

for(i in 1:length(sub.set)){
  lar.inter1[i] <- lar.surv.LR$sims.list$alpha[sub.set[i]]
  lar.slope1[i] <- lar.surv.LR$sims.list$beta[sub.set[i]]
  lar.sigma1[i] <- lar.surv.LR$sims.list$sigma[sub.set[i]]
  lar.inter[i] <- lar.surv.QR$sims.list$alpha[sub.set[i]]
  lar.l.slope[i] <- lar.surv.QR$sims.list$beta[sub.set[i],1]
  lar.q.slope[i] <- lar.surv.QR$sims.list$beta[sub.set[i],2]
  lar.sigma[i] <- lar.surv.QR$sims.list$sigma[sub.set[i]]
}

mean(lar.inter1)
mean(lar.slope1)
ests<-data.frame(cbind(lar.inter1, lar.slope1, lar.sigma1,
                       lar.inter, lar.l.slope, lar.q.slope, lar.sigma))
str(ests)
write.csv(ests,"larval-survival-posterior-samples.csv")

min(exp(-2.38))#0.093

plot(density(lar.surv.LR$sims.list$alpha))
str(lar.surv.LR$sims.list$alpha)
