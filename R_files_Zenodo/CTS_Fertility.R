###California Tiger Salamander (Ambystoma californiense) Fertility Function
###Arianne Messerman
###February 2021

#setwd("")

#Female body mass data when entering breeding pond is "Size" and is already log-transformed
fert <-read.table("Fertility_Function_Hastings.txt",header=T)
fert$size.cent<-fert$Size-mean(fert$Size, na.rm=TRUE)
str(fert)

#Estimated average clutch size in Blomquist population is 814 (Trenham and Shaffer 2005)
x.CM<-mean(fert$ClutchMass) #Average clutch mass = 6.80 g
x.EM<-x.CM/814 #Average egg mass = 0.0084 g

plot(fert$size.cent, fert$ClutchMass)
hist(fert$size.cent)
hist(fert$ClutchMass)

#Convert response variable of the fertility function to clutch size rather than clutch mass by dividing by the average mass of a CTS egg
fert$ClutchSize <- fert$ClutchMass/x.EM

plot(fert$size.cent, fert$ClutchSize)
abline(summary(lm(ClutchSize~size.cent, fert)))
hist(fert$ClutchSize)
hist(fert$Eggs)

library(ggplot2)
library(cowplot)
library(jagsUI)
library(mcmcplots)

set.seed(423467)

##########################################################
##Linear regression of clutch size by female size
#DIC=1755
##########################################################

sink("fert-LR.jags")
cat("
 model {
 
    #Priors
    alpha ~ dnorm(0, 0.0000001) #Intercept
    beta ~ dnorm(0, 0.0000001) #linear slope of female size on clutch size
    sigma ~ dunif(0, 1000)
    
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
jags.data <- list(y = c(fert$ClutchSize), n = length(fert[,1]),  x = c(fert$size.cent))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(1), sigma = runif(1,0,10))}

#Parameters to estimate
parameters <- c("alpha", "beta", "sigma", "fit", "fit.new",
                "bpvalue", "residual", "predicted",  "mu")

#MCMC settings
ni <- 4000
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R
fert.LR <- jags(jags.data, inits, parallel=TRUE, parameters, "fert-LR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(fert.LR, digits = 3)#DIC=1755.82

plot(fert.LR)

#Residual plot
plot(fert.LR$mean$predicted, fert.LR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(fert.LR$sims.list$fit, fert.LR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(fert.LR$sims.list$fit.new>fert.LR$sims.list$fit)#Bayesian p-value=0.52

#Examine posteriors
plot(fert.LR)

##linear fertility equation
pred <- fert.LR$mean$alpha + fert.LR$mean$beta*fert$size.cent

library(dplyr)
fert1 <- dplyr::mutate(fert, mu = fert.LR$mean$mu, #mu are the predicted values
                       q2.5_mu = fert.LR$q2.5$mu, 
                       q97.5_mu = fert.LR$q97.5$mu)
fert1 <- fert1[order(fert1[, "mu"]),]

#Check that mu and function match
plot(fert$size.cent, fert$ClutchSize)
abline(a=-3109.07, b=1126.77)#From lm()
points(fert$size.cent, pred, type="l", lwd=2)
points(fert1$size.cent, fert1$mu, type="l", col=2, lwd=2)


#Figure showing 95% Credible intervals around mean estimate
ggplot(fert1) + 
    geom_ribbon(aes(x = size.cent, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
    geom_path(aes(x = size.cent, y = mu), color = "red") +
    geom_point(aes(x = size.cent, y = ClutchSize)) +
    scale_y_continuous("Clutch size")+
    scale_x_continuous("Ln(female mass (g))")+
    theme_cowplot()

mean(log(fert$ClutchSize))#3.48

##########################################################
##Quadratic regression of clutch size by female size
#DIC=1753.814
##########################################################

sink("fert-QR.jags")
cat("
 model {
 
    #Priors
    alpha ~ dnorm(0, 0.0000001) #Intercept
    for (i in 1:2){
    beta[i] ~ dnorm(0, 0.0000001) #linear and quadratic slopes of female size on clutch size
    }
    sigma ~ dunif(0, 1000)
    
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
jags.data <- list(y = c(fert$ClutchSize), n = length(fert[,1]),  x = c(fert$size.cent))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(2, 0, 1), sigma = runif(1,0,2))}

#Parameters to estimate
parameters <- c("alpha", "beta", "sigma", "fit", "fit.new", "mu",
                "bpvalue", "residual", "predicted")

#MCMC settings
ni <- 4000
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R
fert.QR <- jags(jags.data, inits, parallel=TRUE, parameters, "fert-QR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(fert.QR, digits = 3)#DIC=1753.814

#Residual plot
plot(fert.QR$mean$predicted, fert.QR$mean$residual, las=1, xlab="Predicted values", ylab="Residuals")
abline(h=0)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(fert.QR$sims.list$fit, fert.QR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(fert.QR$sims.list$fit.new>fert.QR$sims.list$fit)#Bayesian p-value=0.52

#Examine posteriors
plot(fert.QR)

##fertility quadratic equation
pred <- fert.QR$mean$alpha + fert.QR$mean$beta[1]*fert$size.cent + fert.QR$mean$beta[2]*(fert$size.cent*fert$size.cent)

library(dplyr)
fert2 <- dplyr::mutate(fert, mu = fert.QR$mean$mu, #mu are the predicted values
                         q2.5_mu = fert.QR$q2.5$mu, 
                         q97.5_mu = fert.QR$q97.5$mu)
fert2 <- fert2[order(fert2[, "mu"]),]

#Check that mu and function match
plot(fert$size.cent, fert$ClutchSize)
lines(x=fert$size.cent, y=pred, lwd=2)
points(fert2$Size, fert2$mu, type="l", col=2, lwd=2)

#Figure showing 95% Credible intervals around mean estimate
ggplot(fert2) + 
    geom_ribbon(aes(x = size.cent, ymin = q2.5_mu, ymax = q97.5_mu), fill = "grey90") +
    geom_path(aes(x = size.cent, y = mu), color = "red") +
    geom_point(aes(x = size.cent, y = ClutchSize)) +
    scale_y_continuous("Clutch size")+
    scale_x_continuous("Ln(female mass (g))")+
    theme_cowplot()

#Allow individual growth trajectories to sample from residual distribution around the function mean
plot(fert$size.cent, fert.QR$mean$residual)
resid.fert<-c(fert.QR$mean$residual)
plot(density(resid.fert))

#Pull 500 sets of iterations from posteriors of coefficients to .CSV
mcmc.sample<-fert.QR$mcmc.info$n.samples
sub.set <- sort(sample(1:mcmc.sample, size=500))

inter <- l.slope <- q.slope <- sigma <- inter1 <- slope1 <- sigma1 <- c(1:length(sub.set))

for(i in 1:length(sub.set)){
    inter1[i] <- fert.LR$sims.list$alpha[sub.set[i]]
    slope1[i] <- fert.LR$sims.list$beta[sub.set[i]]
    sigma1[i] <- fert.LR$sims.list$sigma[sub.set[i]]
    inter[i] <- fert.QR$sims.list$alpha[sub.set[i]]
    l.slope[i] <- fert.QR$sims.list$beta[sub.set[i],1]
    q.slope[i] <- fert.QR$sims.list$beta[sub.set[i],2]
    sigma[i] <- fert.QR$sims.list$sigma[sub.set[i]]
}

ests<-data.frame(cbind(inter1, slope1, sigma1, inter, l.slope, q.slope, sigma))
str(ests)
write.csv(ests,"fertility-posterior-samples-CENTERED.csv")
