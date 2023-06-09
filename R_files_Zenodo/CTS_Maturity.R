###California Tiger Salamander (Ambystoma californiense) Maturity Function
###Arianne Messerman
###February 2021

#setwd("")

#"size" and is already log-transformed
mat <-read.table("Maturity Function.txt",header=T)
mat$size.cent<-mat$size-mean(mat$size, na.rm=TRUE)
str(mat)

library(ggplot2)
library(cowplot)
library(jagsUI)
library(mcmcplots)

set.seed(423467)

plot(mat$size.cent, mat$maturity)
m1<-glm(mat$maturity~mat$size.cent, family="binomial")
summary(m1)

m2<-glm(maturity~size+I(size^2), data=mat, family="binomial")
summary(m2)

##########################################################
##Bernoulli model of maturity by body size -linear structure
#DIC=2952.24
##########################################################
sink("mat-LR.jags")
cat("
    model {
    
    #Priors
    alpha ~ dnorm(0, 0.001) #Intercept
    beta ~ dnorm(0, 0.001) #linear slope of body size on maturity
    
    #Likelihood
    for (i in 1:n){
    C[i] ~ dbern(p[i]) #Bernoulli noise around expected maturation
    logit(p[i]) <- alpha + beta*x[i]
    }
    
    #Assess model fit using sum-of-squares-type discrepancy
    for (i in 1:n){
    residual[i] <- C[i]-p[i]     #Residuals for observed data
    predicted[i] <- p[i]         #Predicted values
    sq[i] <- pow(residual[i],2)   #Squared residuals for observed data
    
    #Generate replicate data and compute fit stats for them
    C.new[i] ~ dbern(p[i])  #One new data set at each MCMC iteration
    sq.new[i] <- pow(C.new[i]-predicted[i],2) #Squared residuals for new data
    }
    
    fit <- sum(sq[])              #Sum of squared residuals for actual data set
    fit.new <- sum(sq.new[])      #Sum of squared residuals for new data set
    test <- step(fit.new - fit)   #Test whether nes data set more extreme
    bpvalue <- mean(test)         #Bayesian p-value
    }
    ", fill=TRUE)
sink()

#Bundle data
jags.data <- list(n = length(mat[,1]), x = c(mat$size.cent), C=c(mat$maturity))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(1))}

#Parameters to estimate
parameters <- c("alpha", "beta", "fit", "fit.new","bpvalue", "p")

#MCMC settings
ni <- 15000
nt <- 5
nb <- 7000
nc <- 3

# Call JAGS from R
mat.LR <- jags(jags.data, inits, parallel=TRUE, parameters, "mat-LR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(mat.LR, digits = 3)#DIC=2952.24

plot(mat.LR)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(mat.LR$sims.list$fit, mat.LR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(mat.LR$sims.list$fit.new>mat.LR$sims.list$fit)#Bayesian p-value=0.50


##linear matility equation
pred <- plogis(mat.LR$mean$alpha + mat.LR$mean$beta*mat$size.cent)

library(dplyr)
mat1 <- dplyr::mutate(mat, p = mat.LR$mean$p, #p are the predicted values
                       q2.5_p = mat.LR$q2.5$p, 
                       q97.5_p = mat.LR$q97.5$p)
mat1 <- mat1[order(mat1[, "p"]),]

#Check that p and function match
plot(mat$size.cent, mat$maturity)
points(mat$size.cent, pred, type="l", lwd=2)
points(mat1$size.cent, mat1$p, type="l", col=2, lwd=2)


#Figure showing 95% Credible intervals around mean estimate
ggplot(mat1) + 
  geom_ribbon(aes(x = size, ymin = q2.5_p, ymax = q97.5_p), fill = "grey90") +
  geom_path(aes(x = size, y = p), color = "red") +
  geom_point(aes(x = size, y = maturity)) +
  scale_y_continuous("Maturity probability")+
  scale_x_continuous("Ln(mass (g))")+
  theme_cowplot()

##########################################################
##Bernoulli model of maturity by body size -quadratic structure
#DIC=3009.68
##########################################################
sink("mat-QR.jags")
cat("
    model {
    
    #Priors
    alpha ~ dnorm(0,1) #Intercept
    for (i in 1:2){
    beta[i] ~ dnorm(0,1) #slopes of body size on maturity
    }
    
    #Likelihood
    for (i in 1:n){
    C[i] ~ dbern(p[i]) #Bernoulli noise around expected maturation
    logit(p[i]) <- alpha + beta[1]*x[i] + beta[2]*pow(x[i],2)
    }
    }
    ", fill=TRUE)
sink()

#Bundle data
jags.data <- list(n = length(mat[,1]), x = c(mat$size.cent), C=c(mat$maturity))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(2, 0, 1))}

#Parameters to estimate
parameters <- c("alpha", "beta", "p")

#MCMC settings
ni <- 18000
nt <- 5
nb <- 8000
nc <- 3

# Call JAGS from R
mat.QR <- jags(jags.data, inits, parallel=TRUE, parameters, "mat-QR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(mat.QR, digits = 3)#DIC=3009.68

plot(mat.QR)

##linear matility equation
pred <- plogis(mat.QR$mean$alpha + mat.QR$mean$beta[1]*mat$size.cent + mat.QR$mean$beta[2]*(mat$size.cent^2))

library(dplyr)
mat2 <- dplyr::mutate(mat, p = mat.QR$mean$p, #p are the predicted values
                      q2.5_p = mat.QR$q2.5$p, 
                      q97.5_p = mat.QR$q97.5$p)
mat2 <- mat2[order(mat2[, "p"]),]

#Check that p and function closely match
plot(mat$size.cent, mat$maturity)
points(mat$size.cent, pred, type="l", lwd=2)
points(mat2$size.cent, mat2$p, type="l", col=2, lwd=2)


#Figure showing 95% Credible intervals around mean estimate
ggplot(mat2) + 
  geom_ribbon(aes(x = size.cent, ymin = q2.5_p, ymax = q97.5_p), fill = "grey90") +
  geom_path(aes(x = size.cent, y = p), color = "red") +
  geom_point(aes(x = size.cent, y = maturity)) +
  scale_y_continuous("Maturity probability")+
  scale_x_continuous("Ln(mass (g))")+
  theme_cowplot()

#Pull 500 sets of iterations from posteriors of coefficients to .CSV
mcmc.sample<-mat.LR$mcmc.info$n.samples
sub.set <- sort(sample(1:mcmc.sample, size=500))

inter <- slope <- c(1:length(sub.set))

for(i in 1:length(sub.set)){
  inter[i] <- mat.LR$sims.list$alpha[sub.set[i]]
  slope[i] <- mat.LR$sims.list$beta[sub.set[i]]
}

ests<-data.frame(cbind(inter, slope))
str(ests)
write.csv(ests,"maturity-posterior-samples-CENTERED.csv")
