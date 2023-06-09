###California Tiger Salamander (Ambystoma californiense) Mass Growth Function
###Arianne Messerman
###October 2020

#setwd("")

#Data already log transformed
meta <-read.table("Metamorph_Growth.txt",header=T)
adult <-read.table("Adult_Growth.txt",header=T)
meta$size.cent<-meta$Size-mean(meta$Size, na.rm=TRUE)
adult$size.cent<-adult$Size-mean(adult$Size, na.rm=TRUE)
str(meta)
str(adult)

plot(adult$size.cent, adult$SizeNext, ylim=c(0,4))
points(meta$size.cent, meta$SizeNext, col=2)

library(ggplot2)
library(cowplot)

ggplot(meta, aes(x=size.cent, y=SizeNext))+
    geom_smooth(method = "lm", alpha=0.15)+
    theme_cowplot()+
    labs(y="Ln(Mass in year t+1 (g))", x="Ln(Mass in year t (g))")

ggplot(adult, aes(x=size.cent, y=SizeNext))+
    geom_smooth(method = "lm", alpha=0.15)+
    theme_cowplot()+
    labs(y="Ln(Mass in year t+1 (g))", x="Ln(Mass in year t (g))")

hist(meta$size.cent)
hist(meta$SizeNext)
hist(adult$size.cent)
hist(adult$SizeNext)


library(jagsUI)
library(mcmcplots)
set.seed(423467)

##########################################################
##Linear regression of metamorph size at t+1 by size at t
#DIC=185.81
##########################################################

sink("meta.growth.LR.jags")
cat("
 model {
 
    #Priors
    alpha ~ dnorm(0, 0.001) #Intercept
    beta ~ dnorm(0, 0.001) #linear  slope of growth
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
jags.data <- list(y = c(meta$SizeNext), n = length(meta[,1]),  x = c(meta$size.cent))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(1), sigma = runif(1,0,2))}

#Parameters to estimate
parameters <- c("alpha", "beta", "sigma", "fit", "fit.new", "mu",
                "bpvalue", "residual", "predicted")

#MCMC settings
ni <- 2000
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R
meta.growth.LR <- jags(jags.data, inits, parallel=TRUE, parameters, "meta.growth.LR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(meta.growth.LR, digits = 3)#DIC=185.81

#Residual plot
plot(meta.growth.LR$mean$predicted, meta.growth.LR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(meta.growth.LR$sims.list$fit, meta.growth.LR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(meta.growth.LR$sims.list$fit.new>meta.growth.LR$sims.list$fit)#Bayesian p-value=0.52

#Examine posteriors
plot(meta.growth.LR)

##Metamorph growth equation
pred.meta <- meta.growth.LR$mean$alpha + meta.growth.LR$mean$beta[1]*meta$size.cent

library(dplyr)
metas <- dplyr::mutate(meta, mu = meta.growth.LR$mean$mu, #mu are the predicted values
                       q2.5_mu = meta.growth.LR$q2.5$mu, 
                       q97.5_mu = meta.growth.LR$q97.5$mu)

#Check that mu and function closely match
plot(meta$size.cent, meta$SizeNext)
points(meta$size.cent, pred.meta, type="l", lwd=2)
points(meta$size.cent, metas$mu, type="l", col=2, lwd=2)

#Figure showing 95% Credible intervals around mean estimate
ggplot(metas) + 
    geom_ribbon(aes(x = size.cent, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
    geom_path(aes(x = size.cent, y = mu), color = "red") +
    geom_point(aes(x = size.cent, y = SizeNext)) +
    scale_y_continuous("Ln(First-year mass (g))")+
    scale_x_continuous("Ln(metamorph mass (g))")+
    theme_cowplot()

mean(meta$Size)#2.39
##########################################################
##Quadratic regression of metamorph size at t+1 by size at t
#DIC=182.56
##########################################################

sink("meta.growth.QR.jags")
cat("
 model {
 
    #Priors
    alpha ~ dnorm(0, 0.001) #Intercept
    for (i in 1:2){
    beta[i] ~ dnorm(0, 0.001) #linear and quadratic slopes of growth
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
jags.data <- list(y = c(meta$SizeNext), n = length(meta[,1]),  x = c(meta$size.cent))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(2, 0, 1), sigma = runif(1,0,2))}

#Parameters to estimate
parameters <- c("alpha", "beta", "sigma", "fit", "fit.new", "mu",
                "bpvalue", "residual", "predicted")

#MCMC settings
ni <- 2000
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R
meta.growth.QR <- jags(jags.data, inits, parallel=TRUE, parameters, "meta.growth.QR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(meta.growth.QR, digits = 3)#DIC=182.56

#Residual plot
plot(meta.growth.QR$mean$predicted, meta.growth.QR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(meta.growth.QR$sims.list$fit, meta.growth.QR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(meta.growth.QR$sims.list$fit.new>meta.growth.QR$sims.list$fit)#Bayesian p-value=0.53

#Examine posteriors
plot(meta.growth.QR)

##Metamorph growth equation
pred.meta <- meta.growth.QR$mean$alpha + meta.growth.QR$mean$beta[1]*meta$size.cent + meta.growth.QR$mean$beta[2]*(meta$size.cent*meta$size.cent)

library(dplyr)
metas <- dplyr::mutate(meta, mu = meta.growth.QR$mean$mu, #mu are the predicted values
                         q2.5_mu = meta.growth.QR$q2.5$mu, 
                         q97.5_mu = meta.growth.QR$q97.5$mu)

#Check that mu and function match
plot(meta$size.cent, meta$SizeNext)
points(meta$size.cent, pred.meta, type="l", lwd=2)
points(meta$size.cent, metas$mu, type="l", col=2, lwd=2)

#Figure showing 95% Credible intervals around mean estimate
ggplot(metas) + 
    geom_ribbon(aes(x = size.cent, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
    geom_path(aes(x = size.cent, y = mu), color = "red") +
    geom_point(aes(x = size.cent, y = SizeNext)) +
    scale_y_continuous("Ln(First-year mass (g))")+
    scale_x_continuous("Ln(metamorph mass (g))")+
    theme_cowplot()

##########################################################
##Linear regression of adult size at t+1 by size at t
#DIC=20.75
##########################################################

sink("adult.growth.LR.jags")
cat("
 model {
 
    #Priors
    alpha ~ dnorm(0, 0.001) #Intercept
    beta ~ dnorm(0, 0.001) #linear slope of growth
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
jags.data <- list(y = c(adult$SizeNext), n = length(adult[,1]),  x = c(adult$size.cent))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(1), sigma = runif(1,0,2))}

#Parameters to estimate
parameters <- c("alpha", "beta", "sigma", "fit", "fit.new", "mu",
                "bpvalue", "residual", "predicted")

#MCMC settings
ni <- 2000
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R
adult.growth.LR <- jags(jags.data, inits, parallel=TRUE, parameters, "adult.growth.LR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(adult.growth.LR, digits = 3)#DIC=20.75

#Residual plot
plot(adult.growth.LR$mean$predicted, adult.growth.LR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(adult.growth.LR$sims.list$fit, adult.growth.LR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(adult.growth.LR$sims.list$fit.new>adult.growth.LR$sims.list$fit)#Bayesian p-value=0.51

#Examine posteriors
plot(adult.growth.LR)

##adultmorph growth equation
pred.adult <- adult.growth.LR$mean$alpha + adult.growth.LR$mean$beta[1]*adult$size.cent

library(dplyr)
adults <- dplyr::mutate(adult, mu = adult.growth.LR$mean$mu, #mu are the predicted values
                        q2.5_mu = adult.growth.LR$q2.5$mu, 
                        q97.5_mu = adult.growth.LR$q97.5$mu)

#Check that mu and function closely match
plot(adult$size.cent, adult$SizeNext)
points(adult$size.cent, pred.adult, type="l", lwd=2)
points(adult$size.cent, adults$mu, type="l", col=2, lwd=2)

#Figure showing 95% Credible intervals around mean estimate
ggplot(adults) + 
    geom_ribbon(aes(x = size.cent, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
    geom_path(aes(x = size.cent, y = mu), color = "red") +
    geom_point(aes(x = size.cent, y = SizeNext)) +
    scale_y_continuous("Ln(mass in year t+1 (g))")+
    scale_x_continuous("Ln(mass in year t (g))")+
    theme_cowplot()

mean(adult$Size)#2.95
##########################################################
##Quadratic regression of adult size at t+1 by size at t
#DIC=22.85
##########################################################

sink("adult.growth.QR.jags")
cat("
 model {
 
    #Priors
    alpha ~ dnorm(0, 0.001) #Intercept
    for (i in 1:2){
    beta[i] ~ dnorm(0, 0.001) #linear and quadratic slopes of growth
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
jags.data <- list(y = c(adult$SizeNext), n = length(adult[,1]),  x = c(adult$size.cent))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(2, 0, 1), sigma = runif(1,0,2))}

#Parameters to estimate
parameters <- c("alpha", "beta", "sigma", "fit", "fit.new", "mu",
                "bpvalue", "residual", "predicted")

#MCMC settings
ni <- 2000
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R
adult.growth.QR <- jags(jags.data, inits, parallel=TRUE, parameters, "adult.growth.QR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(adult.growth.QR, digits = 3)#DIC=22.85

#Residual plot
plot(adult.growth.QR$mean$predicted, adult.growth.QR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(adult.growth.QR$sims.list$fit, adult.growth.QR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(adult.growth.QR$sims.list$fit.new>adult.growth.QR$sims.list$fit)#Bayesian p-value=0.51

#Examine posteriors
plot(adult.growth.QR)

##adult growth equation
pred.adult <- adult.growth.QR$mean$alpha + adult.growth.QR$mean$beta[1]*adult$size.cent + adult.growth.QR$mean$beta[2]*(adult$size.cent*adult$size.cent)

library(dplyr)
adults <- dplyr::mutate(adult, mu = adult.growth.QR$mean$mu, #mu are the predicted values
                       q2.5_mu = adult.growth.QR$q2.5$mu, 
                       q97.5_mu = adult.growth.QR$q97.5$mu)

#Check that mu and function closely match
plot(adult$size.cent, adult$SizeNext)
points(adult$size.cent, pred.adult, type="l", lwd=2)
points(adult$size.cent, adults$mu, type="l", col=2, lwd=2)

#Figure showing 95% Credible intervals around mean estimate
ggplot(adults) + 
    geom_ribbon(aes(x = size.cent, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
    geom_path(aes(x = size.cent, y = mu), color = "red") +
    geom_point(aes(x = size.cent, y = SizeNext)) +
    scale_y_continuous("Ln(mass in year t+1 (g))")+
    scale_x_continuous("Ln(mass in year t (g))")+
    theme_cowplot()



#Pull 500 sets of iterations from posteriors of coefficients to .CSV
mcmc.sample<-adult.growth.QR$mcmc.info$n.samples#Same across life stage models
sub.set <- sort(sample(1:mcmc.sample, size=500))

meta.inter1 <- meta.slope1 <- meta.sigma1  <- c(1:length(sub.set))
meta.inter <- meta.l.slope <- meta.q.slope <- meta.sigma  <- c(1:length(sub.set))
ad.inter1 <- ad.slope1 <- ad.sigma1  <- c(1:length(sub.set))
ad.inter <- ad.l.slope <- ad.q.slope <- ad.sigma <-c(1:length(sub.set))

for(i in 1:length(sub.set)){
    meta.inter1[i] <- meta.growth.LR$sims.list$alpha[sub.set[i]]
    meta.slope1[i] <- meta.growth.LR$sims.list$beta[sub.set[i]]
    meta.sigma1[i] <- meta.growth.LR$sims.list$sigma[sub.set[i]]
    ad.inter1[i] <- adult.growth.LR$sims.list$alpha[sub.set[i]]
    ad.slope1[i] <- adult.growth.LR$sims.list$beta[sub.set[i]]
    ad.sigma1[i] <- adult.growth.LR$sims.list$sigma[sub.set[i]]
    meta.inter[i] <- meta.growth.QR$sims.list$alpha[sub.set[i]]
    meta.l.slope[i] <- meta.growth.QR$sims.list$beta[sub.set[i],1]
    meta.q.slope[i] <- meta.growth.QR$sims.list$beta[sub.set[i],2]
    meta.sigma[i] <- meta.growth.QR$sims.list$sigma[sub.set[i]]
    ad.inter[i] <- adult.growth.QR$sims.list$alpha[sub.set[i]]
    ad.l.slope[i] <- adult.growth.QR$sims.list$beta[sub.set[i],1]
    ad.q.slope[i] <- adult.growth.QR$sims.list$beta[sub.set[i],2]
    ad.sigma[i] <- adult.growth.QR$sims.list$sigma[sub.set[i]]
}

ests<-data.frame(cbind(meta.inter1, meta.slope1, meta.sigma1,
                       meta.inter, meta.l.slope, meta.q.slope, meta.sigma,
                       ad.inter1, ad.slope1, ad.sigma1,
                       ad.inter, ad.l.slope, ad.q.slope, ad.sigma))
str(ests)
write.csv(ests,"growth-posterior-samples-CENTERED.csv")
