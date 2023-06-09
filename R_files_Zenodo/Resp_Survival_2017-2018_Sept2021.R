######################################################################################################################################
## Script by Arianne F. Messerman
## 06 September, 2021
## This script calculates mean residual respiratory surface area water loss (RSAWL), log-body mass, and
##   standard metabolic rate (SMR) values for individual juvenile spotted and marbled salamanders (Ambystoma maculatum and A. opacum).
##   These data are then used in Bayesian analyses to identify the relationships between individual phenotypes and 
##   known survival at three time points that may be indicative of selective pressure. 
## Please see the associated manuscript for full study details: 
##   Messerman, A.F. and M. Leal. The contributions of individual traits to survival among terrestrial juvenile 
##   pond-breeding salamanders. Functional Ecology.
######################################################################################################################################

#setwd()

##Load RSAWL data
data<-read.csv("WVP_RSAWL_2017-2018_Survival.csv",header=TRUE)
set.seed(354)

##Format data
data$round<-as.factor(data$round)
data$chamber<-as.factor(data$chamber)
data$channel<-as.factor(data$channel)
data$tag<-as.factor(data$tag)
data$pen<-as.factor(data$pen)
data$block<-as.factor(data$block)
data$log.svl<-log10(data$svl.test)
data$log.tl<-log10(data$tl.test)
data$log.wide<-log10(data$wide.test)
data$log.mass<-log10(data$premass)

##Add column categorizing site latitude
data$pop [100:150]
data$lat<-ifelse(data$pop=="FLW", "2mid", "1north")
data$lat<-ifelse(data$pop=="Mingo", "3south", data$lat)
data$lat<-as.factor(data$lat)
data$lat [100:150]

##Column of spp*lat interaction
data$sppvlat<-interaction(data$spp,data$lat)

##Column of continuous latitudinal value
data$lat.vec<-ifelse(data$pop=="FLW", 37.7415, data$pop)
data$lat.vec<-ifelse(data$pop=="Mingo", 36.9615, data$lat.vec)
data$lat.vec<-ifelse(data$pop=="Forum", 38.9222, data$lat.vec)
data$lat.vec<-ifelse(data$pop=="Danville", 38.8717, data$lat.vec)
data$lat.vec<-ifelse(data$pop=="DBCA", 38.7861, data$lat.vec)
data$lat.vec<-as.numeric(data$lat.vec)
data$pop [200:300]
data$lat.vec [200:300]

#########################################################
## PCA for body size morphological data with rounds 3-5
#########################################################
# Remove round 1 and 2 in new data frame as flush
data$round<-as.integer(data$round)
unique(data$round)
data3<-subset(data, data$round>=3)
data3$round<-as.factor(data3$round)

pca3 <- prcomp(~ tl.test + svl.test + log.mass + wide.test, data = data3, center = TRUE, scale = TRUE) #The "center" and "scale" arguments are a method of preprossessing, or autoscaling, the data. Often variables in PCA are not reported in the same units or on the same scale. Centering the data substracts the variable mean from each individual observation within a variable, and scaling divides each centered data point by it's respective standard deviation.
pca3 #Scores

#extract eigenvalues and examine Kaiser's Rule
eigen <- pca3$sdev^2

plot(eigen)
abline(h = 1, col = "blue")

summary(pca3)#0.82 proportion variance explained on PC1

#Format data for PCA
pc1and2 <- as.data.frame(pca3$x[,c(1,2)])

pc1.3<-c(pca3$x[,1])#coordinates per individual on PC1
data3$pc1<-pc1.3

#Visualize loadings of individual morphometrics 
biplot(pca3, choices = c(1,2), xlabs=rep(".", nrow(data3)), cex=2)

#PCA Plots
library(ggplot2)
library(dplyr)
pc1plot <- ggplot(pc1and2, aes(PC1, PC2, color = data3$spp)) + geom_point() + labs(x = "PC1", y = "PC2")
pc1plot

unique(data3$spp)
plot(data3$mean.rwl~data3$pc1, col=as.factor(data3$spp), ylab="RSA Water Loss (mg cm^-2 h^-1)", xlab="Coordinates on PC1")
abline(6.96914,0.57622)
legend(3,6,legend=c("AMMA", "AMOP"), pch=1, col=c(1,2))
mod.pc.3<-lm(data3$mean.rwl~data3$pc1)
summary(mod.pc.3)
resid.pc.3<-resid(mod.pc.3)
data3$resid.pc<-resid.pc.3

###############################################
##Build dataframe of averaged RSAWL rounds 3-5
###############################################
names(data3)
mean.df<-aggregate(data3, list(data3$ID), mean)
str(mean.df)
cleaned<-unique.data.frame(data3[,c(1:13,53:54)], incomparables = FALSE)
cleaned<-cleaned[order(cleaned$ID),]
str(cleaned)
rwl<-cbind(cleaned, mean.df$mean.rwl)
rwl<-cbind(rwl, mean.df$days)
rwl<-cbind(rwl, mean.df$premass)
rwl<-cbind(rwl, mean.df$log.mass)
rwl<-cbind(rwl, mean.df$svl.test)
rwl<-cbind(rwl, mean.df$tl.test)
rwl<-cbind(rwl, mean.df$wide.test)
rwl<-cbind(rwl, mean.df$log.svl)
rwl<-cbind(rwl, mean.df$log.tl)
rwl<-cbind(rwl, mean.df$log.wide)
rwl<-cbind(rwl, mean.df$lat.vec)
names(rwl)
str(rwl)
names(rwl)[16]<-paste("mean.rwl")
names(rwl)[17]<-paste("days")
names(rwl)[18]<-paste("premass")
names(rwl)[19]<-paste("log.mass")
names(rwl)[20]<-paste("svl")
names(rwl)[21]<-paste("tl")
names(rwl)[22]<-paste("widest")
names(rwl)[23]<-paste("log.svl")
names(rwl)[24]<-paste("log.tl")
names(rwl)[25]<-paste("log.wide")
names(rwl)[26]<-paste("lat.vec")
str(rwl)
names(rwl)

#PCA Function and Summary
pca <- prcomp(~ tl + svl + log.mass + widest, data = rwl, center = TRUE, scale = TRUE) #The "center" and "scale" arguments are a method of preprossessing, or autoscaling, the data. Often variables in PCA are not reported in the same units or on the same scale. Centering the data substracts the variable mean from each individual observation within a variable, and scaling divides each centered data point by it's respective standard deviation.
pca #Scores

#extract eigenvalues and examine Kaiser's Rule
eigen <- pca$sdev^2

plot(eigen)
abline(h = 1, col = "blue")

summary(pca) #0.82 proportion variance explained on PC1

#Format data for PCA
pc1and2 <- as.data.frame(pca$x[,c(1,2)])

pc1<-c(pca$x[,1])#coordinates per individual on PC1
rwl$pc1<-pc1
str(rwl)

#Visualize loadings of individual morphometrics 
biplot(pca, choices = c(1,2))

#PCA Plots
pc1plot <- ggplot(pc1and2, aes(PC1, PC2, color = rwl$spp)) + geom_point() + labs(x = "PC1", y = "PC2")
pc1plot

plot(rwl$mean.rwl~rwl$pc1, col=as.factor(rwl$spp), ylab="RSA Water Loss (mg cm^-2 h^-1)", xlab="Coordinates on PC1")
abline(6.96955,0.57692)
legend(3,6,legend=c("AMMA", "AMOP"), pch=1, col=c(1,2))
mod.pc<-lm(rwl$mean.rwl~rwl$pc1)
summary(mod.pc)
resid.pc<-resid(mod.pc)
rwl$resid.pc<-resid.pc

species<-c("Spotted Salamander", "Marbled Salamander")
spp<-rwl$spp
unique(spp)
tab<-table(rwl$spp)

plot(rwl$pc1, rwl$mean.rwl, pch=c(as.factor(spp)), col=as.factor(spp),ylab="Respiratory Surface Area Water Loss (mg cm^-2 h^-1)", xlab="Coordinates on PC1", cex=1.2)
legend(2, 5.5,legend=species, col=c(1,2), pch=c(1:5), cex=1.1)
abline(lm(rwl$mean.rwl~rwl$pc1), lwd=2)

plot(as.factor(rwl$spp), rwl$resid.pc)
rwl.amma<-subset(rwl, rwl$spp=="AMMA")
rwl.amop<-subset(rwl, rwl$spp=="AMOP")
mean(rwl.amma$resid.pc)
sd(rwl.amma$resid.pc)
mean(rwl.amop$resid.pc)
sd(rwl.amop$resid.pc)

###############################################
#Add Metabolic Rate data
###############################################
data.mr<-read.csv("Released_Rounds3-5_RSAWL_MINVCO2_2017.csv", header=TRUE)
data.mr$days<-as.factor(data.mr$days)
data.mr$chamber<-as.factor(data.mr$chamber)

###Find mass-specific MR FOR ALL 2017 ANIMALS (Messerman and Leal. 2020. Oecologia)
plot(data.mr$log.mass,data.mr$log.vco2, col=as.factor(data.mr$spp), xlab="Log10 (Body Mass [g])", ylab="Log10 (VCO2 [ml/h])")
abline(-0.91418, 0.86824, col="black")
legend(0.405, -1.4, legend=unique(data.mr$spp), col=as.factor(unique(data.mr$spp)), pch=1, cex=.7)
int.A<-10^-0.91418 #Mass scaling of all 2017 animals coefficient back transformed from log 10

data.mr$mr<-data.mr$min.vco2/(data.mr$premass^0.86824) #Calculate mass-specific MR as VCO2/(W^b)
hist(data.mr$mr)

#Add mr to rwl
rwl<-rwl[order(rwl$ID),]
data.mr<-data.mr[order(data.mr$ID),]

rwl.mr<-cbind(rwl, data.mr[,16:28])
str(rwl.mr)
rwl.mr<-rwl.mr[order(rwl.mr$ID),]
head(rwl.mr)
rwl.mr$ID[1:50]
rwl$ID[1:50]
data.mr$ID[1:50]

#Inspect data
hist(rwl.mr$resid.pc) #Approximately normal
hist(rwl.mr$log.vco2)
hist(rwl.mr$lat.vec)
hist(rwl.mr$premass)

rwl.mr$mr.trans<-(rwl.mr$mr)^(1/2) #Transform mr to normalize predictor
hist(rwl.mr$mr.trans)
rwl.mr$log.mass<-log10(rwl.mr$premass) #Transform body mass to normalize predictor
hist(rwl.mr$log.mass)

#Regression of SMR on PCA coordinates
plot(rwl.mr$pc1, rwl.mr$log.vco2, pch=c(as.factor(spp)), col=as.factor(spp),ylab="SMR", xlab="Coordinates on PC1", cex=1.2)
legend(3, -0.3,legend=species, col=c(1,2), pch=c(1:5), cex=1.1)
abline(lm(rwl.mr$log.vco2~rwl.mr$pc1), lwd=2)

mod.pc.mr<-lm(rwl.mr$log.vco2~rwl.mr$pc1)
summary(mod.pc.mr)
resid.mr<-resid(mod.pc.mr)
rwl.mr$resid.mr<-resid.mr

hist(rwl.mr$resid.mr)

###Center and scale covariates
for (i in 1:length(rwl.mr$lat.vec)) {
  rwl.mr$lat.vec.std[i] <- (rwl.mr$lat.vec[i]-mean(rwl.mr$lat.vec[]))/sd(rwl.mr$lat.vec[])
  rwl.mr$resid.std[i] <- (rwl.mr$resid.pc[i]-mean(rwl.mr$resid.pc[]))/sd(rwl.mr$resid.pc[])
  rwl.mr$mr.std[i] <- (rwl.mr$mr.trans[i]-mean(rwl.mr$mr.trans[]))/sd(rwl.mr$mr.trans[])
  rwl.mr$resid.mr.std[i] <- (rwl.mr$resid.mr[i]-mean(rwl.mr$resid.mr[]))/sd(rwl.mr$resid.mr[])
  rwl.mr$mass.std[i] <- (rwl.mr$log.mass[i]-mean(rwl.mr$log.mass[]))/sd(rwl.mr$log.mass[])
}

####################################################################
#Build Bayesian Bernoulli Distributed Mixed Effects Models 
####################################################################
library(jagsUI)

surv.apr<-as.numeric(rwl.mr$surv.may)
surv.nov<-as.numeric(rwl.mr$surv.nov)
surv.oct<-as.numeric(rwl.mr$surv.oct)
spp<-as.numeric(as.factor(rwl.mr$spp))
pop<-as.numeric(as.factor(rwl.mr$pop))
block<-as.numeric(rwl.mr$block)
pen<-as.numeric(rwl.mr$pen)
n<-length(rwl.mr$ID) #Number of individuals
k.apr<-sum(rwl.mr$surv.may=="1") #Number of individuals that survived to end of study
k.nov<-sum(rwl.mr$surv.nov=="1") #Number of individuals that survived to winter
k.oct<-sum(rwl.mr$surv.oct=="1") #Number of individuals that survived initial two weeks
mass<-as.numeric(rwl.mr$mass.std)#Standardized log10(initial mass (g))
rwl<-as.numeric(rwl.mr$resid.std)#Standardized residual RSAWL on PC1 coordinates
mr<-as.numeric(rwl.mr$resid.mr.std)#Standardized residual SMR on PC1 coordinates

################
#APRIL
################
############################################
#Quadratic GLMM April: Species+Mass+RSAWL+MR Fixed effects
#No population effect 
#With individual random effect
#Block & Pen random effects
############################################
sink("glmm-apr-full2.jags")
cat("
    model{
    
    # Priors
    mu ~ dnorm(0,0.001)I(-15,15)
    mean <- 1/(1+exp(-mu))           # back-transformed grand mean/intercept
    
    beta1[1] <- 0          # Corner constraint on AMMA to avoid overparameterization
    beta1[2] ~  dnorm(0,0.0001)I(-10,10)    # Prior for difference in AMMA and AMOP survival
    
    for(i in 1:6){
    beta[i]~dnorm(0,0.0001)I(-10,10) #Priors for mass at release, residual RSAWL, and residual SMR
    b.est[i]<-1/(1+exp(-beta[i]))
    }
    
    for(i in 1:nblock){
    alpha[i]~dnorm(mu.alpha, tau.alpha) #Prior for random block effect
    a.est[i]<-1/(1+exp(-alpha[i]))
    }
    mu.alpha~dnorm(0,1)
    tau.alpha<-1/(sigma.alpha*sigma.alpha)
    sigma.alpha~dunif(0,10)
    sigma2.alpha <- pow(sigma.alpha, 2)
    
    for(i in 1:npen){
    gamma[i]~dnorm(mu.gamma, tau.gamma) #Prior for random pen effect
    g.est[i]<-1/(1+exp(-gamma[i]))
    }
    mu.gamma~dnorm(0,1)
    tau.gamma<-1/(sigma.gamma*sigma.gamma)
    sigma.gamma~dunif(0,10)
    sigma2.gamma <- pow(sigma.gamma, 2)
    
    for(i in 1:n){
    epsilon[i]~dnorm(mu.epsilon, tau.epsilon) #Prior for residual variation of individual
    }
    mu.epsilon~dnorm(0,1)
    tau.epsilon<-1/(sigma.epsilon*sigma.epsilon)
    sigma.epsilon~dunif(0,10)
    sigma2.epsilon <- pow(sigma.epsilon, 2)
    
    
    #Likelihood
    for(i in 1:n){
    C[i] ~ dbern(p[i]) #Bernoulli noise around expected survival
    p[i] <- (1/(1 + exp(-(mu + beta1[spp[i]] + beta[1]*mass[i] + beta[2]*pow(mass[i],2) + 
        beta[3]*rwl[i] + beta[4]*pow(rwl[i],2) + beta[5]*mr[i] +beta[6]*pow(mr[i],2) + alpha[block[i]] + 
        gamma[pen[i]] + epsilon[i]))))
    
    #Computation of fit statistic (Bayesian p-value)
    Presi[i] <- abs(C[i]-p[i])  #Absolute Pearson residuals
    C.new[i] ~ dbern(p[i])
    Presi.new[i] <- abs(C.new[i]-p[i])
    }
    
    fit <-sum(Presi[])  #Discrepancy for actual data set
    fit.new <-sum(Presi.new[])  #Discrepancy for replicate data set
    }
    ", fill=TRUE)

sink()

#Bundle data
jags.data <- list(C=surv.apr, n=n, spp=spp, mass=mass, rwl=rwl, mr=mr, 
                  nblock = length(unique(block)), block = as.numeric(block), 
                  npen = length(unique(pen)), pen = as.numeric(pen))

#Inits function
inits <- function(){list(mu = runif(1, 0, 1), beta1 = c(NA, rnorm(1,0,1)), 
                         beta2 = c(NA, rnorm(1,0,1), rnorm(1,0,1), rnorm(1,0,1)), 
                         beta = runif(6, -3, 3))}

#Parameters to estimate
parameters <- c("mu", "mean", "beta1","beta", "b.est", 
                "alpha", "gamma", "sigma2.epsilon", "Presi", "fit", "fit.new", "p")

# MCMC settings
ni <- 150000
nt <- 5
nb <- 80000
nc <- 3

# Call JAGS from R
glmm.apr.full2 <- jags(jags.data, parallel=TRUE, inits, parameters, "glmm-apr-full2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

plot(glmm.apr.full2$mean$Presi, ylab="Absolute Residual")
plot(glmm.apr.full2$sims.list$fit, glmm.apr.full2$sims.list$fit.new, main="", xlab="Discrepancy for actual data set",
     ylab="Discrepancy for perfect data set", las=1)
abline(0,1,lwd=2)

mean(glmm.apr.full2$sims.list$fit.new>glmm.apr.full2$sims.list$fit)#Bayesian p-value=0.47

print(glmm.apr.full2, digits = 3)
#DIC=80.80

glmm.apr.full2$mean$deviance
plot(glmm.apr.full2)


min(glmm.apr.full2$mean$alpha)
max(glmm.apr.full2$mean$alpha)
mean(glmm.apr.full2$sims.list$alpha)
sd(glmm.apr.full2$sims.list$alpha)

min(glmm.apr.full2$mean$gamma)
max(glmm.apr.full2$mean$gamma)
mean(glmm.apr.full2$sims.list$gamma)
sd(glmm.apr.full2$sims.list$gamma)

################
#NOVEMBER
################
############################################
#Quadratic GLMM November: Species+Mass+RSAWL+MR Fixed effects
#No population effect 
#With individual random effect
#Block & Pen random effects
#Add truncation on mu
############################################
sink("glmm-nov-full2.jags")
cat("
    model{
    
    # Priors
    mu ~ dnorm(0,0.001)I(-15,15)
    mean <- 1/(1+exp(-mu))           # back-transformed grand mean/intercept
    
    beta1[1] <- 0          # Corner constraint on AMMA to avoid overparameterization
    beta1[2] ~  dnorm(0,0.000001)I(-10,10)    # Prior for difference in AMMA and AMOP survival
    
    for(i in 1:6){
    beta[i]~dnorm(0,0.00001)I(-10,10) #Priors for mass at release, residual RSAWL, and residual SMR
    b.est[i]<-1/(1+exp(-beta[i]))
    }
    
    for(i in 1:nblock){
    alpha[i]~dnorm(mu.alpha, tau.alpha) #Prior for random block effect
    a.est[i]<-1/(1+exp(-alpha[i]))
    }
    mu.alpha~dnorm(0,1)
    tau.alpha<-1/(sigma.alpha*sigma.alpha)
    sigma.alpha~dunif(0,10)
    sigma2.alpha <- pow(sigma.alpha, 2)
    
    for(i in 1:npen){
    gamma[i]~dnorm(mu.gamma, tau.gamma) #Prior for random pen effect
    g.est[i]<-1/(1+exp(-gamma[i]))
    }
    mu.gamma~dnorm(0,1)
    tau.gamma<-1/(sigma.gamma*sigma.gamma)
    sigma.gamma~dunif(0,10)
    sigma2.gamma <- pow(sigma.gamma, 2)
    
    for(i in 1:n){
    epsilon[i]~dnorm(mu.epsilon, tau.epsilon) #Prior for residual variation of individual
    }
    mu.epsilon~dnorm(0,1)
    tau.epsilon<-1/(sigma.epsilon*sigma.epsilon)
    sigma.epsilon~dunif(0,10)
    sigma2.epsilon <- pow(sigma.epsilon, 2)
    
    
    #Likelihood
    for(i in 1:n){
    C[i] ~ dbern(p[i]) #Bernoulli noise around expected survival
    p[i] <- (1/(1 + exp(-(mu + beta1[spp[i]] + beta[1]*mass[i] + beta[2]*pow(mass[i],2) + 
        beta[3]*rwl[i] + beta[4]*pow(rwl[i],2) + beta[5]*mr[i] +beta[6]*pow(mr[i],2) + alpha[block[i]] + 
        gamma[pen[i]] + epsilon[i]))))
    
    #Computation of fit statistic (Bayesian p-value)
    Presi[i] <- abs(C[i]-p[i])  #Absolute Pearson residuals
    C.new[i] ~ dbern(p[i])
    Presi.new[i] <- abs(C.new[i]-p[i])
    }
    
    fit <-sum(Presi[])  #Discrepancy for actual data set
    fit.new <-sum(Presi.new[])  #Discrepancy for replicate data set
    }
    ", fill=TRUE)

sink()

#Bundle data
jags.data <- list(C=surv.nov, n=n, spp=spp, mass=mass, rwl=rwl, mr=mr, 
                  nblock = length(unique(block)), block = as.numeric(block), 
                  npen = length(unique(pen)), pen = as.numeric(pen))

#Inits function
inits <- function(){list(mu = runif(1, 0, 1), beta1 = c(NA, rnorm(1,0,1)), 
                         beta2 = c(NA, rnorm(1,0,1), rnorm(1,0,1), rnorm(1,0,1)), 
                         beta = runif(6, -3, 3))}

#Parameters to estimate
parameters <- c("mu", "mean", "beta1","beta", "b.est", 
                "alpha", "gamma", "sigma2.epsilon", "Presi", "fit", "fit.new", "p")

# MCMC settings
ni <- 150000
nt <- 5
nb <- 80000
nc <- 3

# Call JAGS from R
glmm.nov.full2 <- jags(jags.data, parallel=TRUE, inits, parameters, "glmm-nov-full2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

plot(glmm.nov.full2$mean$Presi, ylab="Absolute Residual")
plot(glmm.nov.full2$sims.list$fit, glmm.nov.full2$sims.list$fit.new, main="", xlab="Discrepancy for actual data set",
     ylab="Discrepancy for perfect data set", las=1)
abline(0,1,lwd=2)

mean(glmm.nov.full2$sims.list$fit.new>glmm.nov.full2$sims.list$fit)#Bayesian p-value=0.48

print(glmm.nov.full2, digits = 3)
#DIC=97.75

glmm.nov.full2$mean$deviance
plot(glmm.nov.full2)

min(glmm.nov.full2$mean$alpha)
max(glmm.nov.full2$mean$alpha)
mean(glmm.nov.full2$sims.list$alpha)
sd(glmm.nov.full2$sims.list$alpha)

min(glmm.nov.full2$mean$gamma)
max(glmm.nov.full2$mean$gamma)
mean(glmm.nov.full2$sims.list$gamma)
sd(glmm.nov.full2$sims.list$gamma)


################
#OCTOBER
################
############################################
#Quadratic GLMM October: Species+Mass+RSAWL+MR Fixed effects
#No population effect 
#With individual random effect
#Block & Pen random effects
#Add truncation on mu
############################################
sink("glmm-oct-full2.jags")
cat("
    model{
    
    # Priors
    mu ~ dnorm(0,0.001)I(-10,10)
    mean <- 1/(1+exp(-mu))           # back-transformed grand mean/intercept
    
    beta1[1] <- 0          # Corner constraint on AMMA to avoid overparameterization
    beta1[2] ~  dnorm(0,0.000001)I(-10,10)    # Prior for difference in AMMA and AMOP survival
    
    for(i in 1:6){
    beta[i]~dnorm(0,0.00001)I(-10,10) #Priors for mass at release, residual RSAWL, and residual SMR
    b.est[i]<-1/(1+exp(-beta[i]))
    }
    
    for(i in 1:nblock){
    alpha[i]~dnorm(mu.alpha, tau.alpha) #Prior for random block effect
    a.est[i]<-1/(1+exp(-alpha[i]))
    }
    mu.alpha~dnorm(0,1)
    tau.alpha<-1/(sigma.alpha*sigma.alpha)
    sigma.alpha~dunif(0,10)
    sigma2.alpha <- pow(sigma.alpha, 2)
    
    for(i in 1:npen){
    gamma[i]~dnorm(mu.gamma, tau.gamma) #Prior for random pen effect
    g.est[i]<-1/(1+exp(-gamma[i]))
    }
    mu.gamma~dnorm(0,1)
    tau.gamma<-1/(sigma.gamma*sigma.gamma)
    sigma.gamma~dunif(0,10)
    sigma2.gamma <- pow(sigma.gamma, 2)
    
    for(i in 1:n){
    epsilon[i]~dnorm(mu.epsilon, tau.epsilon) #Prior for residual variation of individual
    }
    mu.epsilon~dnorm(0,1)
    tau.epsilon<-1/(sigma.epsilon*sigma.epsilon)
    sigma.epsilon~dunif(0,10)
    sigma2.epsilon <- pow(sigma.epsilon, 2)
    
    
    #Likelihood
    for(i in 1:n){
    C[i] ~ dbern(p[i]) #Bernoulli noise around expected survival
    p[i] <- (1/(1 + exp(-(mu + beta1[spp[i]] + beta[1]*mass[i] + beta[2]*pow(mass[i],2) + 
        beta[3]*rwl[i] + beta[4]*pow(rwl[i],2) + beta[5]*mr[i] +beta[6]*pow(mr[i],2) + alpha[block[i]] + 
        gamma[pen[i]] + epsilon[i]))))
    
    #Computation of fit statistic (Bayesian p-value)
    Presi[i] <- abs(C[i]-p[i])  #Absolute Pearson residuals
    C.new[i] ~ dbern(p[i])
    Presi.new[i] <- abs(C.new[i]-p[i])
    }
    
    fit <-sum(Presi[])  #Discrepancy for actual data set
    fit.new <-sum(Presi.new[])  #Discrepancy for replicate data set
    }
    ", fill=TRUE)

sink()

#Bundle data
jags.data <- list(C=surv.oct, n=n, spp=spp, mass=mass, rwl=rwl, mr=mr, 
                  nblock = length(unique(block)), block = as.numeric(block), 
                  npen = length(unique(pen)), pen = as.numeric(pen))

#Inits function
inits <- function(){list(mu = runif(1, 0, 1), beta1 = c(NA, rnorm(1,0,1)), 
                         beta2 = c(NA, rnorm(1,0,1), rnorm(1,0,1), rnorm(1,0,1)), 
                         beta = runif(6, -3, 3))}

#Parameters to estimate
parameters <- c("mu", "mean", "beta1","beta", "b.est", 
                "alpha", "gamma", "sigma2.epsilon", "Presi", "fit", "fit.new", "p")

# MCMC settings
ni <- 150000
nt <- 5
nb <- 80000
nc <- 3

# Call JAGS from R
glmm.oct.full2 <- jags(jags.data, parallel=TRUE, inits, parameters, "glmm-oct-full2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

plot(glmm.oct.full2$mean$Presi, ylab="Absolute Residual")
plot(glmm.oct.full2$sims.list$fit, glmm.oct.full2$sims.list$fit.new, main="", xlab="Discrepancy for actual data set",
     ylab="Discrepancy for perfect data set", las=1)
abline(0,1,lwd=2)

mean(glmm.oct.full2$sims.list$fit.new>glmm.oct.full2$sims.list$fit)#Bayesian p-value=0.48

print(glmm.oct.full2, digits = 3)
#DIC=100.83

glmm.oct.full2$mean$deviance
plot(glmm.oct.full2)

min(glmm.oct.full2$mean$alpha)
max(glmm.oct.full2$mean$alpha)
mean(glmm.oct.full2$sims.list$alpha)
sd(glmm.oct.full2$sims.list$alpha)

min(glmm.oct.full2$mean$gamma)
max(glmm.oct.full2$mean$gamma)
mean(glmm.oct.full2$sims.list$gamma)
sd(glmm.oct.full2$sims.list$gamma)

###########################################################
##Predictive Logistic Regression across all groups
###########################################################
mean.mass<-mean(rwl.mr$log.mass)
mean.rwl<-mean(rwl.mr$resid.pc)
mean.mr<-mean(rwl.mr$resid.mr)
sd.mass<-sd(rwl.mr$log.mass)
sd.rwl<-sd(rwl.mr$resid.pc)
sd.mr<-sd(rwl.mr$resid.mr)

mcmc.sample<-glmm.apr.full2$mcmc.info$n.samples
original.mass.pred<-seq(min(rwl.mr$log.mass), max(rwl.mr$log.mass), length.out=length(rwl.mr$log.mass))
original.rwl.pred<-seq(min(rwl.mr$resid.pc), max(rwl.mr$resid.pc), length.out=length(rwl.mr$resid.pc))
original.mr.pred<-seq(min(rwl.mr$resid.mr), max(rwl.mr$resid.mr), length.out=length(rwl.mr$resid.mr))
mass.pred<-(original.mass.pred-mean.mass)/sd.mass
rwl.pred<-(original.rwl.pred-mean.rwl)/sd.rwl
mr.pred<-(original.mr.pred-mean.mr)/sd.mr
sub.set <- sort(sample(1:mcmc.sample, size=200))

par(mfrow=c(3,3))

#For October using full quadratic model
array.p.pred.mass2<-array.p.pred.rwl2<-array.p.pred.mr2<-array(NA, dim=c(length(mass.pred), mcmc.sample))
for(i in 1:mcmc.sample){
  array.p.pred.mass2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i] + glmm.oct.full2$sims.list$beta[i,1]*mass.pred + glmm.oct.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i] + glmm.oct.full2$sims.list$beta[i,3]*rwl.pred + glmm.oct.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i] + glmm.oct.full2$sims.list$beta[i,5]*mr.pred + glmm.oct.full2$sims.list$beta[i,6]*mr.pred^2)
}
p.pred.mass2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta[1]*mass.pred+glmm.oct.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta[3]*rwl.pred+glmm.oct.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta[5]*mr.pred+glmm.oct.full2$mean$beta[6]*mr.pred^2)

plot(original.mass.pred, p.pred.mass2, ylab="Survival Probability", xlab="Mass", type="l",lwd=3, ylim=c(0,1), col=4)
for( i in sub.set){
  lines(original.mass.pred, array.p.pred.mass2[,i], type="l", lwd=1, col="gray")
}
lines(original.mass.pred, p.pred.mass2, ylab="Survival Probability", xlab="Mass", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$mass.std, rwl.mr$surv.oct)
plot(original.rwl.pred, p.pred.rwl2, ylab="Survival Probability", xlab="RSAWL", type="l",lwd=3, ylim=c(0,1), col=4)
for( i in sub.set){
  lines(original.rwl.pred, array.p.pred.rwl2[,i], type="l", lwd=1, col="gray")
}
lines(original.rwl.pred, p.pred.rwl2, ylab="Survival Probability", xlab="RSAWL", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$resid.std, rwl.mr$surv.oct)
plot(original.mr.pred, p.pred.mr2, ylab="Survival Probability", xlab="SMR", type="l",lwd=3, ylim=c(0,1), col=4)
for( i in sub.set){
  lines(original.mr.pred, array.p.pred.mr2[,i], type="l", lwd=1, col="gray")
}
lines(original.mr.pred, p.pred.mr2, ylab="Survival Probability", xlab="SMR", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$resid.mr.std, rwl.mr$surv.oct)

#For November using full quadratic model
array.p.pred.mass1<-array.p.pred.rwl1<-array.p.pred.mr1<-array(NA, dim=c(length(mass.pred), mcmc.sample))
for(i in 1:mcmc.sample){
  array.p.pred.mass1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i] + glmm.nov.full2$sims.list$beta[i,1]*mass.pred + glmm.nov.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i] + glmm.nov.full2$sims.list$beta[i,3]*rwl.pred + glmm.nov.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i] + glmm.nov.full2$sims.list$beta[i,5]*mr.pred + glmm.nov.full2$sims.list$beta[i,6]*mr.pred^2)
}
p.pred.mass1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta[1]*mass.pred+glmm.nov.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta[3]*rwl.pred+glmm.nov.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta[5]*mr.pred+glmm.nov.full2$mean$beta[6]*mr.pred^2)

plot(original.mass.pred, p.pred.mass1, ylab="Survival Probability", xlab="Mass", type="l",lwd=3, ylim=c(0,1), col=4)
for( i in sub.set){
  lines(original.mass.pred, array.p.pred.mass1[,i], type="l", lwd=1, col="gray")
}
lines(original.mass.pred, p.pred.mass1, ylab="Survival Probability", xlab="Mass", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$mass.std, rwl.mr$surv.nov)
plot(original.rwl.pred, p.pred.rwl1, ylab="Survival Probability", xlab="RSAWL", type="l",lwd=3, ylim=c(0,1), col=4)
for( i in sub.set){
  lines(original.rwl.pred, array.p.pred.rwl1[,i], type="l", lwd=1, col="gray")
}
lines(original.rwl.pred, p.pred.rwl1, ylab="Survival Probability", xlab="RSAWL", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$resid.std, rwl.mr$surv.nov)
plot(original.mr.pred, p.pred.mr1, ylab="Survival Probability", xlab="SMR", type="l",lwd=3, ylim=c(0,1), col=4)
for( i in sub.set){
  lines(original.mr.pred, array.p.pred.mr1[,i], type="l", lwd=1, col="gray")
}
lines(original.mr.pred, p.pred.mr1, ylab="Survival Probability", xlab="SMR", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$resid.mr.std, rwl.mr$surv.nov)

#For April using full quadratic model
array.p.pred.mass<-array.p.pred.rwl<-array.p.pred.mr<-array(NA, dim=c(length(mass.pred), mcmc.sample))
for(i in 1:mcmc.sample){
  array.p.pred.mass[,i] <- plogis(glmm.apr.full2$sims.list$mu[i] + glmm.apr.full2$sims.list$beta[i,1]*mass.pred + glmm.apr.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl[,i] <- plogis(glmm.apr.full2$sims.list$mu[i] + glmm.apr.full2$sims.list$beta[i,3]*rwl.pred + glmm.apr.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr[,i] <- plogis(glmm.apr.full2$sims.list$mu[i] + glmm.apr.full2$sims.list$beta[i,5]*mr.pred + glmm.apr.full2$sims.list$beta[i,6]*mr.pred^2)
}

p.pred.mass<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta[1]*mass.pred+glmm.apr.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta[3]*rwl.pred+glmm.apr.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta[5]*mr.pred+glmm.apr.full2$mean$beta[6]*mr.pred^2)

plot(original.mass.pred, p.pred.mass, ylab="Survival Probability", xlab="Mass", type="l",lwd=3, ylim=c(0,1), frame.plot = FALSE)
for( i in sub.set){
  lines(original.mass.pred, array.p.pred.mass[,i], type="l", lwd=1, col="gray")
}
lines(original.mass.pred, p.pred.mass, ylab="Survival Probability", xlab="Mass", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$mass.std, rwl.mr$surv.may)
plot(original.rwl.pred, p.pred.rwl, ylab="Survival Probability", xlab="RSAWL", type="l",lwd=3, ylim=c(0,1), col=4)
for( i in sub.set){
  lines(original.rwl.pred, array.p.pred.rwl[,i], type="l", lwd=1, col="gray")
}
lines(original.rwl.pred, p.pred.rwl, ylab="Survival Probability", xlab="RSAWL", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$resid.std, rwl.mr$surv.may)
plot(original.mr.pred, p.pred.mr, ylab="Survival Probability", xlab="SMR", type="l",lwd=3, ylim=c(0,1), col=4)
for( i in sub.set){
  lines(original.mr.pred, array.p.pred.mr[,i], type="l", lwd=1, col="gray")
}
lines(original.mr.pred, p.pred.mr, ylab="Survival Probability", xlab="SMR", type="l",lwd=3, ylim=c(0,1), col=4)
points(rwl.mr$resid.mr.std, rwl.mr$surv.may)


###############################################################################
##Plot relationship between covariates and survival probability by species
###############################################################################
amma<-subset(rwl.mr, rwl.mr$spp=="AMMA")
amop<-subset(rwl.mr, rwl.mr$spp=="AMOP")

#For October using full quadratic model
array.p.pred.mass.amma2<-array.p.pred.rwl.amma2<-array.p.pred.mr.amma2<-array(NA, dim=c(length(mass.pred), mcmc.sample))
array.p.pred.mass.amop2<-array.p.pred.rwl.amop2<-array.p.pred.mr.amop2<-array(NA, dim=c(length(mass.pred), mcmc.sample))
for(i in 1:mcmc.sample){
  array.p.pred.mass.amma2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i] + glmm.oct.full2$sims.list$beta1[i,1] + glmm.oct.full2$sims.list$beta[i,1]*mass.pred + glmm.oct.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl.amma2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i]+ glmm.oct.full2$sims.list$beta1[i,1] + glmm.oct.full2$sims.list$beta[i,3]*rwl.pred + glmm.oct.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr.amma2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i]+ glmm.oct.full2$sims.list$beta1[i,1] + glmm.oct.full2$sims.list$beta[i,5]*mr.pred + glmm.oct.full2$sims.list$beta[i,6]*mr.pred^2)
  array.p.pred.mass.amop2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i] + glmm.oct.full2$sims.list$beta1[i,2] + glmm.oct.full2$sims.list$beta[i,1]*mass.pred + glmm.oct.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl.amop2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i]+ glmm.oct.full2$sims.list$beta1[i,2] + glmm.oct.full2$sims.list$beta[i,3]*rwl.pred + glmm.oct.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr.amop2[,i] <- plogis(glmm.oct.full2$sims.list$mu[i]+ glmm.oct.full2$sims.list$beta1[i,2] + glmm.oct.full2$sims.list$beta[i,5]*mr.pred + glmm.oct.full2$sims.list$beta[i,6]*mr.pred^2)
}

p.pred.mass.amma2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta1[1]+glmm.oct.full2$mean$beta[1]*mass.pred+glmm.oct.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl.amma2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta1[1]+glmm.oct.full2$mean$beta[3]*rwl.pred+glmm.oct.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr.amma2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta1[1]+glmm.oct.full2$mean$beta[5]*mr.pred+glmm.oct.full2$mean$beta[6]*mr.pred^2)
p.pred.mass.amop2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta1[2]+glmm.oct.full2$mean$beta[1]*mass.pred+glmm.oct.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl.amop2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta1[2]+glmm.oct.full2$mean$beta[3]*rwl.pred+glmm.oct.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr.amop2<-plogis(glmm.oct.full2$mean$mu+glmm.oct.full2$mean$beta1[2]+glmm.oct.full2$mean$beta[5]*mr.pred+glmm.oct.full2$mean$beta[6]*mr.pred^2)

max(amma$mass.std)#2.285
min(amma$mass.std)#-2.587
max(amop$mass.std)#1.525=[142]
min(amop$mass.std)#-1.741=[30]

#Mass all spp
set.seed(34)
pdf("Fig2.pdf", width = 12, height = 8)
#tiff("Fig2-resid.tiff", width = 6400, height =5600, units = 'px', res=800, compression = 'lzw')
par(mfrow=c(3,3))
par(mar = c(3.7, 3.7, 0.5, 0.5) + 1)
plot(original.mass.pred, p.pred.mass, ylab="Survival probability", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white",
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext(side=3,at = -.2,"a", cex=1.5)
for( i in sub.set){
  lines(original.mass.pred, array.p.pred.mass.amma2[,i], type="l", lwd=1, col="lightblue1")
  lines(original.mass.pred, array.p.pred.mass.amop2[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.mass.pred, p.pred.mass.amma2, ylab="Survival probability", xlab="Mass", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
lines(original.mass.pred[30:142], p.pred.mass.amop2[30:142], ylab="Survival probability", xlab="Mass", type="l",lwd=4, ylim=c(0,1), col="salmon1")
points(jitter(amma$log.mass, amount=.01), jitter(amma$surv.oct, amount=.1), col="deepskyblue3")
points(jitter(amop$log.mass, amount=.01), jitter(amop$surv.oct, amount=.1), col="salmon1")


#RSAWL all spp
max(amma$resid.std)#1.715
min(amma$resid.std)#-2.042
max(amop$resid.std)#3.042
min(amop$resid.std)#-2.388

plot(original.rwl.pred, p.pred.rwl, ylab="", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white",
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext(side=3,at = -1.1,"b", cex=1.5)
for( i in sub.set){
  lines(original.rwl.pred, array.p.pred.rwl.amma2[,i], type="l", lwd=1, col="lightblue1")
  lines(original.rwl.pred, array.p.pred.rwl.amop2[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.rwl.pred, p.pred.rwl.amop2, ylab="Survival probability", xlab="rwl", type="l",lwd=4, ylim=c(0,1), col="salmon1")
lines(original.rwl.pred[12:127], p.pred.rwl.amma2[12:127], ylab="Survival probability", xlab="rwl", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
points(jitter(amma$resid.pc, amount=.1), jitter(amma$surv.oct, amount=.1), col="deepskyblue3")
points(jitter(amop$resid.pc, amount=.1), jitter(amop$surv.oct, amount=.1), col="salmon1")

#SMR all spp
max(amma$resid.mr.std)#1.899=[168]
min(amma$resid.mr.std)#-2.955=[38]
max(amop$resid.mr.std)#1.520=[157]
min(amop$resid.mr.std)#-4.301=[1]

plot(original.mr.pred, p.pred.mr, ylab="", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white",
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext(side=3,at = -1.1,"c", cex=1.5)
for( i in sub.set){
  lines(original.mr.pred, array.p.pred.mr.amma2[,i], type="l", lwd=1, col="lightblue1")
  lines(original.mr.pred, array.p.pred.mr.amop2[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.mr.pred[38:168], p.pred.mr.amma2[38:168], ylab="Survival probability", xlab="mr", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
lines(original.mr.pred[1:157], p.pred.mr.amop2[1:157], ylab="Survival probability", xlab="mr", type="l",lwd=4, ylim=c(0,1), col="salmon1")
points(jitter(amma$resid.mr, amount=.01), jitter(amma$surv.oct, amount=.1), col="deepskyblue3")
points(jitter(amop$resid.mr, amount=.01), jitter(amop$surv.oct, amount=.1), col="salmon1")


#For November using full quadratic model
array.p.pred.mass.amma1<-array.p.pred.rwl.amma1<-array.p.pred.mr.amma1<-array(NA, dim=c(length(mass.pred), mcmc.sample))
array.p.pred.mass.amop1<-array.p.pred.rwl.amop1<-array.p.pred.mr.amop1<-array(NA, dim=c(length(mass.pred), mcmc.sample))
for(i in 1:mcmc.sample){
  array.p.pred.mass.amma1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i] + glmm.nov.full2$sims.list$beta1[i,1] + glmm.nov.full2$sims.list$beta[i,1]*mass.pred + glmm.nov.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl.amma1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i]+ glmm.nov.full2$sims.list$beta1[i,1] + glmm.nov.full2$sims.list$beta[i,3]*rwl.pred + glmm.nov.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr.amma1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i]+ glmm.nov.full2$sims.list$beta1[i,1] + glmm.nov.full2$sims.list$beta[i,5]*mr.pred + glmm.nov.full2$sims.list$beta[i,6]*mr.pred^2)
  array.p.pred.mass.amop1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i] + glmm.nov.full2$sims.list$beta1[i,2] + glmm.nov.full2$sims.list$beta[i,1]*mass.pred + glmm.nov.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl.amop1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i]+ glmm.nov.full2$sims.list$beta1[i,2] + glmm.nov.full2$sims.list$beta[i,3]*rwl.pred + glmm.nov.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr.amop1[,i] <- plogis(glmm.nov.full2$sims.list$mu[i]+ glmm.nov.full2$sims.list$beta1[i,2] + glmm.nov.full2$sims.list$beta[i,5]*mr.pred + glmm.nov.full2$sims.list$beta[i,6]*mr.pred^2)
}

p.pred.mass.amma1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta1[1]+glmm.nov.full2$mean$beta[1]*mass.pred+glmm.nov.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl.amma1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta1[1]+glmm.nov.full2$mean$beta[3]*rwl.pred+glmm.nov.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr.amma1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta1[1]+glmm.nov.full2$mean$beta[5]*mr.pred+glmm.nov.full2$mean$beta[6]*mr.pred^2)
p.pred.mass.amop1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta1[2]+glmm.nov.full2$mean$beta[1]*mass.pred+glmm.nov.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl.amop1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta1[2]+glmm.nov.full2$mean$beta[3]*rwl.pred+glmm.nov.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr.amop1<-plogis(glmm.nov.full2$mean$mu+glmm.nov.full2$mean$beta1[2]+glmm.nov.full2$mean$beta[5]*mr.pred+glmm.nov.full2$mean$beta[6]*mr.pred^2)

#Mass all spp
plot(original.mass.pred, p.pred.mass, ylab="Survival probability", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white", 
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext(side=3,at = -.2,"d", cex=1.5)
for( i in sub.set){
  lines(original.mass.pred, array.p.pred.mass.amma1[,i], type="l", lwd=1, col="lightblue1")
  lines(original.mass.pred, array.p.pred.mass.amop1[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.mass.pred, p.pred.mass.amma1, ylab="Survival probability", xlab="Mass", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
lines(original.mass.pred[30:142], p.pred.mass.amop1[30:142], ylab="Survival probability", xlab="Mass", type="l",lwd=4, ylim=c(0,1), col="salmon1")
points(jitter(amma$log.mass, amount=.01), jitter(amma$surv.nov, amount=.1), col="deepskyblue3")
points(jitter(amop$log.mass, amount=.01), jitter(amop$surv.nov, amount=.1), col="salmon1")

#RSAWL all spp
plot(original.rwl.pred, p.pred.rwl, ylab="", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white", 
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext(side=3,at = -1.1,"e", cex=1.5)
for( i in sub.set){
  lines(original.rwl.pred, array.p.pred.rwl.amma1[,i], type="l", lwd=1, col="lightblue1")
  lines(original.rwl.pred, array.p.pred.rwl.amop1[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.rwl.pred, p.pred.rwl.amop1, ylab="Survival probability", xlab="rwl", type="l",lwd=4, ylim=c(0,1), col="salmon1")
lines(original.rwl.pred[12:127], p.pred.rwl.amma1[12:127], ylab="Survival probability", xlab="rwl", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
points(jitter(amma$resid.pc, amount=.1), jitter(amma$surv.nov, amount=.1), col="deepskyblue3")
points(jitter(amop$resid.pc, amount=.1), jitter(amop$surv.nov, amount=.1), col="salmon1")

#SMR all spp
plot(original.mr.pred, p.pred.mr, ylab="", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white", 
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext(side=3,at = -1.1,"f", cex=1.5)
for( i in sub.set){
  lines(original.mr.pred, array.p.pred.mr.amma1[,i], type="l", lwd=1, col="lightblue1")
  lines(original.mr.pred, array.p.pred.mr.amop1[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.mr.pred[38:168], p.pred.mr.amma1[38:168], ylab="Survival probability", xlab="mr", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
lines(original.mr.pred[1:157], p.pred.mr.amop1[1:157], ylab="Survival probability", xlab="mr", type="l",lwd=4, ylim=c(0,1), col="salmon1")
points(jitter(amma$resid.mr, amount=.01), jitter(amma$surv.nov, amount=.1), col="deepskyblue3")
points(jitter(amop$resid.mr, amount=.01), jitter(amop$surv.nov, amount=.1), col="salmon1")

#For April using full quadratic model
array.p.pred.mass.amma<-array.p.pred.rwl.amma<-array.p.pred.mr.amma<-array(NA, dim=c(length(mass.pred), mcmc.sample))
array.p.pred.mass.amop<-array.p.pred.rwl.amop<-array.p.pred.mr.amop<-array(NA, dim=c(length(mass.pred), mcmc.sample))
for(i in 1:mcmc.sample){
  array.p.pred.mass.amma[,i] <- plogis(glmm.apr.full2$sims.list$mu[i] + glmm.apr.full2$sims.list$beta1[i,1] + glmm.apr.full2$sims.list$beta[i,1]*mass.pred + glmm.apr.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl.amma[,i] <- plogis(glmm.apr.full2$sims.list$mu[i]+ glmm.apr.full2$sims.list$beta1[i,1] + glmm.apr.full2$sims.list$beta[i,3]*rwl.pred + glmm.apr.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr.amma[,i] <- plogis(glmm.apr.full2$sims.list$mu[i]+ glmm.apr.full2$sims.list$beta1[i,1] + glmm.apr.full2$sims.list$beta[i,5]*mr.pred + glmm.apr.full2$sims.list$beta[i,6]*mr.pred^2)
  array.p.pred.mass.amop[,i] <- plogis(glmm.apr.full2$sims.list$mu[i] + glmm.apr.full2$sims.list$beta1[i,2] + glmm.apr.full2$sims.list$beta[i,1]*mass.pred + glmm.apr.full2$sims.list$beta[i,2]*mass.pred^2)
  array.p.pred.rwl.amop[,i] <- plogis(glmm.apr.full2$sims.list$mu[i]+ glmm.apr.full2$sims.list$beta1[i,2] + glmm.apr.full2$sims.list$beta[i,3]*rwl.pred + glmm.apr.full2$sims.list$beta[i,4]*rwl.pred^2)
  array.p.pred.mr.amop[,i] <- plogis(glmm.apr.full2$sims.list$mu[i]+ glmm.apr.full2$sims.list$beta1[i,2] + glmm.apr.full2$sims.list$beta[i,5]*mr.pred + glmm.apr.full2$sims.list$beta[i,6]*mr.pred^2)
}

p.pred.mass.amma<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta1[1]+glmm.apr.full2$mean$beta[1]*mass.pred+glmm.apr.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl.amma<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta1[1]+glmm.apr.full2$mean$beta[3]*rwl.pred+glmm.apr.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr.amma<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta1[1]+glmm.apr.full2$mean$beta[5]*mr.pred+glmm.apr.full2$mean$beta[6]*mr.pred^2)
p.pred.mass.amop<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta1[2]+glmm.apr.full2$mean$beta[1]*mass.pred+glmm.apr.full2$mean$beta[2]*mass.pred^2)
p.pred.rwl.amop<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta1[2]+glmm.apr.full2$mean$beta[3]*rwl.pred+glmm.apr.full2$mean$beta[4]*rwl.pred^2)
p.pred.mr.amop<-plogis(glmm.apr.full2$mean$mu+glmm.apr.full2$mean$beta1[2]+glmm.apr.full2$mean$beta[5]*mr.pred+glmm.apr.full2$mean$beta[6]*mr.pred^2)

#Mass all spp
plot(original.mass.pred, p.pred.mass, ylab="Survival probability", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white", 
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext("Log(body mass (g))", side=1, line=3, cex=1.2)
mtext(side=3,at = -0.2,"g", cex=1.5)
for( i in sub.set){
  lines(original.mass.pred, array.p.pred.mass.amma[,i], type="l", lwd=1, col="lightblue1")
  lines(original.mass.pred, array.p.pred.mass.amop[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.mass.pred, p.pred.mass.amma, ylab="Survival probability", xlab="Mass", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
lines(original.mass.pred[30:142], p.pred.mass.amop[30:142], ylab="Survival probability", xlab="Mass", type="l",lwd=4, ylim=c(0,1), col="salmon1")
points(jitter(amma$log.mass, amount=.01), jitter(amma$surv.may, amount=.1), col="deepskyblue3")
points(jitter(amop$log.mass, amount=.01), jitter(amop$surv.may, amount=.1), col="salmon1")

#RSAWL all spp
plot(original.rwl.pred, p.pred.rwl, ylab="", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white", 
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext("Mean residual RSAWL", side=1, line=3, cex=1.2)#expression(RSAWL ~ (mg ~ cm^{-2} ~ h^{-1}))
mtext(side=3,at = -1.1, "h", cex=1.5)
for( i in sub.set){
  lines(original.rwl.pred, array.p.pred.rwl.amma[,i], type="l", lwd=1, col="lightblue1")
  lines(original.rwl.pred, array.p.pred.rwl.amop[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.rwl.pred, p.pred.rwl.amop, ylab="Survival probability", xlab="rwl", type="l",lwd=4, ylim=c(0,1), col="salmon1")
lines(original.rwl.pred[12:127], p.pred.rwl.amma[12:127], ylab="Survival probability", xlab="rwl", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
points(jitter(amma$resid.pc, amount=.01), jitter(amma$surv.may, amount=.1), col="deepskyblue3")
points(jitter(amop$resid.pc, amount=.01), jitter(amop$surv.may, amount=.1), col="salmon1")

#SMR all spp
plot(original.mr.pred, p.pred.mr, ylab="", xlab="", type="l",lwd=4, ylim=c(-0.1,1.1), col="white",
     cex=2, cex.lab=2,cex.axis=1.7, bty="l")
mtext("Residual SMR", side=1, line=3.5, cex=1.2)
mtext(side=3,at = -1.1,"i", cex=1.5)
for( i in sub.set){
  lines(original.mr.pred, array.p.pred.mr.amma[,i], type="l", lwd=1, col="lightblue1")
  lines(original.mr.pred, array.p.pred.mr.amop[,i], type="l", lwd=1, col="rosybrown1")
}
lines(original.mr.pred[38:168], p.pred.mr.amma[38:168], ylab="Survival probability", xlab="mr", type="l",lwd=4, ylim=c(0,1), col="deepskyblue3")
lines(original.mr.pred[1:157], p.pred.mr.amop[1:157], ylab="Survival probability", xlab="mr", type="l",lwd=4, ylim=c(0,1), col="salmon1")
points(jitter(amma$resid.mr,amount=.01), jitter(amma$surv.may, amount=.1), col="deepskyblue3")
points(jitter(amop$resid.mr, amount=.01), jitter(amop$surv.may, amount=.1), col="salmon1")
dev.off()

#############################################################################
#Check whether growth/days correlated with individual SMR
#############################################################################
mod.mass<-lm(log(mass.day+1)~resid.mr, data=rwl.mr, na.action=na.omit)
summary(mod.mass)
plot(rwl.mr$resid.mr, log(rwl.mr$mass.day+1))
abline(0.014339, -0.01315)

mod.svl<-lm(log(svl.day+1)~resid.mr, data=rwl.mr)
summary(mod.svl)
plot(rwl.mr$resid.mr, log(rwl.mr$svl.day+1))
abline(0.0246, 0.04058)

mod.tl<-lm(log(tl.day+1)~resid.mr, data=rwl.mr)
summary(mod.tl)
plot(rwl.mr$resid.mr, log(rwl.mr$tl.day+1))
abline(0.01779, 0.01012)

#############################################################################
##Descriptive plot of abiotic conditions
#############################################################################
abiotic<-read.csv("abiotic.csv", header=TRUE)
abiotic$date<-as.Date(abiotic$date, format="%m/%d/%Y")

for (i in 1:length(abiotic$date)) {
  abiotic$temp.std[i] <- (abiotic$x.temp[i]-mean(abiotic$x.temp[]))/sd(abiotic$x.temp[])
  abiotic$temp.sd.std[i] <- (abiotic$sd.temp[i]-mean(abiotic$sd.temp[]))/sd(abiotic$sd.temp[])
  abiotic$precip.std[i] <- (abiotic$accum.precip[i]-mean(abiotic$accum.precip[]))/sd(abiotic$accum.precip[])
}
str(abiotic)

pdf("Fig1.pdf", width = 10, height = 8)
#tiff("Abiotic-300dpi.tiff", width = 10, height = 8, units = 'in', res=300, compression = 'none')
par(mar = c(5, 4, 4, 4) + 1.5)
plot(abiotic$date, abiotic$x.temp, ylim= c(-10,25), lwd=2, cex=2, pch=1, col=2, bty="u",
     xlab="Date", ylab=expression(Mean ~ temperature ~ (degree~C)), cex.axis=2, cex.lab=2)
lines(abiotic$date, abiotic$x.temp,lty=1, lwd=2, col=2)
segments(abiotic$date, abiotic$x.temp+abiotic$sd.temp, abiotic$date, abiotic$x.temp-abiotic$sd.temp, col=2, lwd=2)
par(new = TRUE)
plot(abiotic$date, abiotic$accum.precip, pch=2, col=4, 
     axes = FALSE, bty = "n", xlab = "", ylab = "", lwd=2, cex=2, cex.axis=2, cex.lab=2)
lines(abiotic$date, abiotic$accum.precip, lty=2, col=4, lwd=2)
axis(side=4, at = pretty(range(abiotic$accum.precip)), cex.axis=2)
mtext("Interval-specific precipitation (mm)", las=3, side=4, line=3, cex=2, cex.axis=2, cex.lab=2) 
dev.off()

#########################################################
#Examine potential interactions between species and traits 
  #For October time period only
  #Simplified model structure
#########################################################
library(lme4)
library(car)

mod.inter<-glm(surv.oct ~ spp + mass.std + I(mass.std^2) + 
                resid.mr.std + I(resid.mr.std^2) + 
                spp:mass.std + spp:resid.mr.std, 
                family="binomial", data=rwl.mr)
summary(mod.inter)
Anova(mod.inter, type="II")

mod.inter1<-glm(surv.oct ~ spp + mass.std + resid.mr.std + 
                 spp:mass.std + spp:resid.mr.std, 
                 family="binomial", data=rwl.mr)
summary(mod.inter1)
Anova(mod.inter1, type="II")

mod.inter2<-glm(surv.oct ~ spp + mass.std + I(mass.std^2) + 
                 resid.mr.std + I(resid.mr.std^2) + 
                 spp:mass.std + spp:resid.mr.std +
                 spp:I(mass.std^2) + spp:I(resid.mr.std^2),
                 family="binomial", data=rwl.mr)
summary(mod.inter2)
Anova(mod.inter2, type="II")