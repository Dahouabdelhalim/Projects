###California Tiger Salamander (Ambystoma californiense) breeding females given December-January precipitation
###Arianne Messerman
###February 2021

#setwd("")

fem <-read.delim("Breeding_Females.txt")
mean(fem$Precip)#221.09
fem$precip.cent<-fem$Precip-mean(fem$Precip)
fem$log.females<-log(fem$Females)

hist(fem$Females)
hist(fem$log.females)

plot(fem$precip.cent, fem$log.females)

library(jagsUI)
library(mcmcplots)
set.seed(423467)

##########################################################
## Linear regression of females breeding given December-January precipitation
## DIC=23.22
##########################################################

sink("fem.LR.jags")
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
jags.data <- list(y = c(fem$log.females), n = length(fem[,1]),  x = c(fem$precip.cent))

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
fem.LR <- jags(jags.data, inits, parallel=TRUE, parameters, "fem.LR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(fem.LR, digits = 3)#DIC=23.23

#Residual plot
plot(fem.LR$mean$predicted, fem.LR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(fem.LR$sims.list$fit, fem.LR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(fem.LR$sims.list$fit.new>fem.LR$sims.list$fit)#Bayesian p-value=0.59

#Examine posteriors
plot(fem.LR)

##Proportion females breeding equation
pred.fem <- fem.LR$mean$alpha + fem.LR$mean$beta[1]*fem$precip.cent

library(dplyr)
fems <- dplyr::mutate(fem, mu = fem.LR$mean$mu, #mu are the predicted values
                       q2.5_mu = fem.LR$q2.5$mu, 
                       q97.5_mu = fem.LR$q97.5$mu)

#Check that mu and function closely match
plot(fem$precip.cent, fem$log.females)
points(fem$precip.cent, pred.fem, type="l", lwd=2)
points(fem$precip.cent, fems$mu, type="l", col=2, lwd=2)

library(ggplot2)
library(cowplot)
#Figure showing 95% Credible intervals around mean estimate
ggplot(fems) + 
  geom_ribbon(aes(x = Precip, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
  geom_path(aes(x = Precip, y = mu), color = "red") +
  geom_point(aes(x = Precip, y = log.females)) +
  scale_y_continuous("Log(proportion of females breeding)")+
  scale_x_continuous("December-January precipitation (mm)")+
  theme_cowplot()

##########################################################
##Quadratic regression of females breeding given precipitation
#DIC=DIC=29.84
##########################################################
sink("fem.QR.jags")
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
jags.data <- list(y = c(fem$log.females), n = length(fem[,1]),  x = c(fem$precip.cent))

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
fem.QR <- jags(jags.data, inits, parallel=TRUE, parameters, "fem.QR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(fem.QR, digits = 3)#DIC=29.84

#Residual plot
plot(fem.QR$mean$predicted, fem.QR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(fem.QR$sims.list$fit, fem.QR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(fem.QR$sims.list$fit.new>fem.QR$sims.list$fit)#Bayesian p-value=0.59

#Examine posteriors
plot(fem.QR)

#Pull 500 sets of iterations from posteriors of coefficients to .CSV
mcmc.sample<-fem.LR$mcmc.info$n.samples#Same across life stage models
sub.set <- sort(sample(1:mcmc.sample, size=500))

fem.inter1 <- fem.slope1 <- fem.sigma1  <- c(1:length(sub.set))

for(i in 1:length(sub.set)){
  fem.inter1[i] <- fem.LR$sims.list$alpha[sub.set[i]]
  fem.slope1[i] <- fem.LR$sims.list$beta[sub.set[i]]
  fem.sigma1[i] <- fem.LR$sims.list$sigma[sub.set[i]]
}

fem.LR$mean$alpha
fem.LR$mean$beta
mean(fem.inter1)
mean(fem.slope1)
ests<-data.frame(cbind(fem.inter1, fem.slope1, fem.sigma1))
str(ests)
write.csv(ests,"females-precip-posterior-samples-CENTERED.csv")
