###California Tiger Salamander (Ambystoma californiense) adult replacement and reproductive failure functions
###Arianne Messerman
###April 2021

#setwd("")

#################################################################################
#Load needed components from IPM script, if not already loaded
#################################################################################
source("CTS_Source.R")

#Set up IPM parameter vectors to inform functions
surv<-read.csv("survival-posterior-samples-CENTERED.csv", header = T)
grow<-read.csv("growth-posterior-samples-CENTERED.csv", header = T)
mat<-read.csv("maturity-posterior-samples-CENTERED.csv", header = T)
fert<-read.csv("fertility-posterior-samples-CENTERED.csv", header = T)
Ndraws<-500
cts<-matrix(NA,nrow=22,ncol=Ndraws)

## Params 01-03: Metamorph growth
cts[1,]<-grow$meta.inter1      	## growth intercept
cts[2,]<-grow$meta.slope1				## growth linear slope
cts[3,]<-grow$meta.sigma1				## SD of metamorph growth function

## Params 04-06: Adult growth
cts[4,]<-grow$ad.inter1      	  ## growth intercept
cts[5,]<-grow$ad.slope1				  ## growth linear slope
cts[6,]<-grow$ad.sigma1					## SD of adult growth function

## Params 07-11: Survival params
cts[7,]<-surv$mean.meta.eta.phi ## mean metamorph survival intercept for years 1,2,6, and 7
cts[8,]<-surv$mass.alpha				## survival on mass slope
cts[9,]<-surv$ad.eta.phi1			  ## adult survival intercept
cts[10,]<-surv$fidel.m	        ## metamorph site fidelity
cts[11,]<-surv$fidel.a	        ## adult site fidelity

## Params 12-13: Maturity params
cts[12,]<-mat$inter   		      ## maturity intercept
cts[13,]<-mat$slope			        ## maturity slope

## Params 14-15: Fertility params (clutch size)
cts[14,]<-fert$inter1    	      ## fertility intercept 
cts[15,]<-fert$slope1		        ## fertility linear slope

## Params 16-27: Misc params (bounds of continuous size domain in units of mean-centered ln(g))
cts[16,]<- -0.85                ## min size overall
cts[17,]<- 5.2                  ## max size overall
cts[18,]<- 2.24                 ## mean size of new metamorphs at low density
cts[19,]<- 0.345                ## SD size of new metamorphs at low density
cts[20,]<- 0.50                 ## Proportion of females in breeding population (Smith and Voss 2009)
cts[21,]<- 0.341                ## Proportion of mature females that return to breed in a given year (Trenham and Shaffer 2005)
cts[22,]<- 0.092                ## Maximum larval survival probability at low density (Trenham and Shaffer 2005)

matsize<-122
lower<-cts[16]
upper<-cts[17]

n<-matsize
L<-lower
U<-upper
h<-(U-L)/n          #Bin size =0.05
b<-L+c(0:n)*h       #Lower boundaries of bins 
y<-0.5*(b[1:n]+b[2:(n+1)])  #Bin midpoints

Mpmat<-array(0, dim=c(length(y), length(y), Ndraws))
for(i in 1:Ndraws){
  for(s in 1:n){
    #Metamorph to one-year-old growth and survival (without bin size correction)
    Mpmat[,,i]<-t(outer(y,y,pxy.m,params=cts[,i]))
  }
}

################################################################################################
#Size distribution of metamorphs in each cohort multiplied by the growth and survival functions 
#to determine the percentage of metamorphs expected to make it to maturity
################################################################################################
#Load size distribution of each cohort across the 122 body size bins (equivalent to 'rec' objects in IPM script)
#where bins in which no metamorphs were measured are in interpolated from the overall cohort size distribution
#Complete reproductive failure (zero meatmophs produced) in 2007 and 2008
m05<-as.matrix(read.table("Olcott2005.txt", header=T))
m06<-as.matrix(read.table("Olcott2006.txt", header=T))
m09<-as.matrix(read.table("Olcott2009.txt", header=T))
m10<-as.matrix(read.table("Olcott2010.txt", header=T))
m11<-as.matrix(read.table("Olcott2011.txt", header=T))
m12<-as.matrix(read.table("Olcott2012.txt", header=T))
m13<-as.matrix(read.table("Olcott2013.txt", header=T))

#Examine number of metamorphs recruited in each study year
sum(m05)
sum(m06)
sum(m09)
sum(m10)
sum(m11)
sum(m12)
sum(m13)

#Calculate recruitment probability for each cohort at size y
rec.prob05<-as.matrix(m05/sum(m05))
rec.prob06<-as.matrix(m06/sum(m06))
rec.prob09<-as.matrix(m09/sum(m09))
rec.prob10<-as.matrix(m10/sum(m10))
rec.prob11<-as.matrix(m11/sum(m11))
rec.prob12<-as.matrix(m12/sum(m12))
rec.prob13<-as.matrix(m13/sum(m13))

#Calculate proportion of metamorphs expected to make it to maturity
M05<-M06<-M07<-M08<-M09<-M10<-M11<-M12<-M13<-array(0, dim=c(length(y), Ndraws))
for(i in 1:Ndraws){
  for(s in 1:n){
    #Size distribution of metamorphs in each cohort multiplied by growth*survival
    M05[s,i]<-(sum(Mpmat[,s,i]))*h*rec.prob05[s,]
    M06[s,i]<-(sum(Mpmat[,s,i]))*h*rec.prob06[s,]
    M09[s,i]<-(sum(Mpmat[,s,i]))*h*rec.prob09[s,]
    M10[s,i]<-(sum(Mpmat[,s,i]))*h*rec.prob10[s,]
    M11[s,i]<-(sum(Mpmat[,s,i]))*h*rec.prob11[s,]
    M12[s,i]<-(sum(Mpmat[,s,i]))*h*rec.prob12[s,]
    M13[s,i]<-(sum(Mpmat[,s,i]))*h*rec.prob13[s,]
  }
}

#Calculate total cohort proportion reaching maturity across body sizes
perc05<-perc06<-perc07<-perc08<-perc09<-perc10<-perc11<-perc12<-perc13<-array(0, dim=c(1, Ndraws))
for (i in 1:Ndraws){
  perc05[,i]<-sum(M05[,i])
  perc06[,i]<-sum(M06[,i])
  perc07[,i]<-sum(M07[,i])
  perc08[,i]<-sum(M08[,i])
  perc09[,i]<-sum(M09[,i])
  perc10[,i]<-sum(M10[,i])
  perc11[,i]<-sum(M11[,i])
  perc12[,i]<-sum(M12[,i])
  perc13[,i]<-sum(M13[,i])
}

#Multiply this proportion by the total number of metamorphs in the cohort
n05<-n06<-n07<-n08<-n09<-n10<-n11<-n12<-n13<-array(0, dim=c(1, Ndraws))
for (i in 1:Ndraws){
  n05[,i]<-perc05[,i]*sum(m05)
  n06[,i]<-perc06[,i]*sum(m06)
  n07[,i]<-perc07[,i]*0
  n08[,i]<-perc08[,i]*0
  n09[,i]<-perc09[,i]*sum(m09)
  n10[,i]<-perc10[,i]*sum(m10)
  n11[,i]<-perc11[,i]*sum(m11)
  n12[,i]<-perc12[,i]*sum(m12)
  n13[,i]<-perc13[,i]*sum(m13)
}

#Mean surviving one-year-olds across 500 random samples from growth and survival distributions
mean(n05)
mean(n06)
mean(n07)
mean(n08)
mean(n09)
mean(n10)
mean(n11)
mean(n12)
mean(n13)

#The upper cutoff of the relationship between precipitation and reproductive success is the inflection point 
#from a logistic regression between October-June precipitation and the probability that metamorph recruitment is
#above the replacement rate (adult carrying capacity * annual adult mortality rate).
#i.e., that the number of adults lost is less than the number of metamorphs recruited.

#Median adults-only lost after one year starting from carrying capacity from the IPM script
#Corrected for drift fence area of ~16% of Olcott's shoreline
A.lost<-529.03

#Is metamorph recruitment >= equivalent adult mortality at K?
mean(n05)-A.lost
mean(n06)-A.lost
mean(n07)-A.lost 
mean(n08)-A.lost 
mean(n09)-A.lost
mean(n10)-A.lost
mean(n11)-A.lost
mean(n12)-A.lost
mean(n13)-A.lost

repl<-c(1,1,0,0,0,1,1,0,0) #Vector indicating replacement of dead adults, where 1=yes, 0=no per study year

############################################################################
##Logistic regression of reproductive success/failure given precipitation
############################################################################
##Load precipitaion data
precip<-read.csv("precip.csv")
mean(precip$oct.jun.precip)#58.7
precip$precip.cent<-precip$oct.jun.precip-mean(precip$oct.jun.precip)

precip$fail<-c(1,1,0,0,1,1,1,1,1)#Add column indicating reproductive failure (0) due to pond drying as observed in the study
precip$replace<-repl #Add column indicating replacement success (1) or failure (0) due to pond drying

plot(precip$oct.jun.precip, precip$fail)
plot(precip$precip.cent, precip$replace)

library(jagsUI)
library(mcmcplots)
set.seed(423467)

#############################################################################
##Bernoulli model of reproductive success by precipitation -linear structure
#DIC=12.16
#############################################################################
sink("pre-LR.jags")
cat("
    model {
    
    #Priors
    alpha ~ dnorm(0, 0.001) #Intercept
    beta ~ dnorm(0, 0.001) #linear slope of reproductive success on precipitation
    
    #Likelihood
    for (i in 1:n){
    C[i] ~ dbern(p[i]) #Bernoulli noise around expected reproduction
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
jags.data <- list(n = length(precip[,1]), x = c(precip$precip.cent), C=c(precip$fail))

#Inits function
inits <- function(){list(alpha=rnorm(1), beta=rnorm(1))}

#Parameters to estipree
parameters <- c("alpha", "beta", "fit", "fit.new","bpvalue", "p")

#MCMC settings
ni <- 15000
nt <- 5
nb <- 7000
nc <- 3

# Call JAGS from R
pre.LR <- jags(jags.data, inits, parallel=TRUE, parameters, "pre-LR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(pre.LR, digits = 3)

plot(pre.LR)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(pre.LR$sims.list$fit, pre.LR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(pre.LR$sims.list$fit.new>pre.LR$sims.list$fit)#Bayesian p-value=0.19

##linear reproductive failure equation
pred <- plogis(pre.LR$mean$alpha + (pre.LR$mean$beta*precip$precip.cent))
pred <- sort(pred)
precip1 <- precip[order(precip$precip.cent),]

library(dplyr)
pre1 <- dplyr::mutate(precip, p = pre.LR$mean$p, #p are the predicted values
                       q2.5_p = pre.LR$q2.5$p, 
                       q97.5_p = pre.LR$q97.5$p)
pre1 <- pre1[order(pre1[, "p"]),]
pre1$pred <- pred

#Examine predictions
plot(precip1$precip.cent, precip1$fail)
points(precip1$precip.cent, pred, type="l", lwd=2)
points(pre1$precip.cent, pre1$p, type="l", col=2, lwd=2)

#Find inflection point
library(inflection)
check_curve(pre1$oct.jun.precip, pre1$pred)
bede(pre1$oct.jun.precip, pre1$pred, 0)#40.45 cm

library(ggplot2)
library(cowplot)
#Figure showing 95% Credible intervals around mean probability of reproductive success
ggplot(pre1) + 
  geom_ribbon(aes(x = oct.jun.precip, ymin = q2.5_p, ymax = q97.5_p), fill = "grey90") +
  geom_path(aes(x = oct.jun.precip, y = p), color = "red") +
  geom_point(aes(x = oct.jun.precip, y = fail)) +
  scale_y_continuous("Reproductive success probability")+
  scale_x_continuous("October-June precipitation (cm)")+
  theme_cowplot()

##########################################################
##Bernoulli model of replacement success by precipitation -linear structure
#DIC=0.63
##########################################################
sink("rep-LR.jags")
cat("
    model {
    
    #Priors (**altered priors and initial values to achieve convergence of beta)
    alpha ~ dnorm(0, 0.01) #Intercept 
    beta ~ dnorm(0, 0.1)I(-10,10) #linear slope of replacement success on precipitation
    
    #Likelihood
    for (i in 1:n){
    C[i] ~ dbern(p[i]) #Bernoulli noise around expected replacement
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
jags.data <- list(n = length(precip[,1]), x = c(precip$precip.cent), C=c(precip$replace))

#Inits function
inits <- function(){list(alpha=rlnorm(1), beta=rlnorm(1))}

#Parameters to estipree
parameters <- c("alpha", "beta", "fit", "fit.new","bpvalue", "p")

#MCMC settings
ni <- 2550000
nt <- 5
nb <- 500000
nc <- 3

# Call JAGS from R
rep.LR <- jags(jags.data, inits, parallel=TRUE, parameters, "rep-LR.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(rep.LR, digits = 3)

plot(rep.LR)

#Bayesian p-value plot, where values closer to 0.5 = better model fit
plot(rep.LR$sims.list$fit, rep.LR$sims.list$fit.new, las=1, xlab="SSQ for actual data set",
     ylab="SSQ for ideal data set")
abline(0,1)
mean(rep.LR$sims.list$fit.new>rep.LR$sims.list$fit)#Bayesian p-value=0.05


##Mean linear reproductive failure equation
pred1 <- plogis(rep.LR$mean$alpha + rep.LR$mean$beta*precip$precip.cent)
pred1 <- sort(pred1)

library(dplyr)
pre2 <- dplyr::mutate(precip, p = rep.LR$mean$p, #p are the predicted values
                      q2.5_p = rep.LR$q2.5$p, 
                      q97.5_p = rep.LR$q97.5$p)
pre2 <- pre2[order(pre2[, "p"]),]
pre2$pred1 <- pred1


#Examine predictions
plot(precip1$precip.cent, precip1$replace)
points(precip1$precip.cent, pred1, type="l", lwd=2)
points(pre2$precip.cent, pre2$p, type="l", col=2, lwd=2)

#Find inflection point
library(inflection)
check_curve(pre2$oct.jun.precip, pre2$pred1)
bede(pre2$oct.jun.precip, pre2$pred1, 0)#55.65 cm

library(ggplot2)
library(cowplot)
#Figure showing 95% Credible intervals around mean estipree
ggplot(pre2) + 
  geom_ribbon(aes(x = oct.jun.precip, ymin = q2.5_p, ymax = q97.5_p), fill = "grey90") +
  geom_path(aes(x = oct.jun.precip, y = p), color = "red") +
  geom_point(aes(x = oct.jun.precip, y = replace)) +
  scale_y_continuous("Reproductive success probability")+
  scale_x_continuous("October-June precipitation (cm)")+
  theme_cowplot()


###################################################################################
#Pull 500 sets of iterations from posteriors of coefficients to .CSV
###################################################################################
##For reproduction
mcmc.sample<-pre.LR$mcmc.info$n.samples
sub.set <- sort(sample(1:mcmc.sample, size=500))

inter <- slope <- c(1:length(sub.set))

for(i in 1:length(sub.set)){
  inter[i] <- pre.LR$sims.list$alpha[sub.set[i]]
  slope[i] <- pre.LR$sims.list$beta[sub.set[i]]
}

ests<-data.frame(cbind(inter, slope))
str(ests)
write.csv(ests,"repro-success-posterior-samples-CENTERED.csv")

##For replacement
mcmc.sample1<-rep.LR$mcmc.info$n.samples
sub.set1 <- sort(sample(1:mcmc.sample1, size=500))

inter <- slope <- c(1:length(sub.set1))

for(i in 1:length(sub.set1)){
  inter[i] <- rep.LR$sims.list$alpha[sub.set1[i]]
  slope[i] <- rep.LR$sims.list$beta[sub.set1[i]]
}

ests<-data.frame(cbind(inter, slope))
str(ests)
write.csv(ests,"replacement-success-posterior-samples-CENTERED.csv")

#############################################################################
## Find inflection points for each MCMC iteration
#############################################################################
replace1<-read.csv("replacement-success-posterior-samples-CENTERED.csv", header=TRUE)
repro1<-read.csv("repro-success-posterior-samples-CENTERED.csv", header=TRUE)

p.replace<-array(NA, dim=c(1, length(precip1[,1]), rep.LR$mcmc.info$n.samples))
p.repro<-array(NA, dim=c(1, length(precip1[,1]), pre.LR$mcmc.info$n.samples))

for(i in 1:pre.LR$mcmc.info$n.sample){
  p.repro[,,i] <- plogis(pre.LR$sims.list$alpha[i] + (pre.LR$sims.list$beta[i]*precip$precip.cent))
  p.repro[,,i]<- sort(p.repro[,,i])
}

for(i in 1:rep.LR$mcmc.info$n.sample){
  p.replace[,,i] <- plogis(rep.LR$sims.list$alpha[i] + (rep.LR$sims.list$beta[i]*precip$precip.cent))
  p.replace[,,i]<- sort(p.replace[,,i])
}

par(mfrow=c(1,1))

plot(precip1$oct.jun.precip, precip1$fail)
for(i in 1:pre.LR$mcmc.info$n.samples){
  lines(precip1$oct.jun.precip, p.repro[,,i])
}

plot(precip1$oct.jun.precip, precip1$replace)
for(i in 1:rep.LR$mcmc.info$n.samples){
  lines(precip1$oct.jun.precip, p.replace[,,i])
}

library(inflection)

#Calculate inflection point for each iteration of each function
inf.repro<-array(NA, dim=c(1, pre.LR$mcmc.info$n.samples))
inf.replace<-array(NA, dim=c(1, rep.LR$mcmc.info$n.samples))
for (i in 1: pre.LR$mcmc.info$n.samples){
  inf1<-bede(pre2$oct.jun.precip, p.repro[,,i], 0)
  inf.repro[,i]<-inf1$iters$EDE[1]
}
for (i in 1: rep.LR$mcmc.info$n.samples){
  inf2<-bede(pre2$oct.jun.precip, p.replace[,,i], 0)
  inf.replace[,i]<-inf2$iters$EDE[1]
}

hist(inf.repro)
hist(inf.replace)
median(inf.repro, na.rm = TRUE)#40.45
median(inf.replace, na.rm = TRUE)#55.65

inf.repro1<-x<-inf.repro[!is.na(inf.repro)]
inf.replace1<-x<-inf.replace[!is.na(inf.replace)]

##Save inflection estimates for 500 random draws
##For reproduction
sub.set <- sort(sample(1:length(inf.repro1), size=500))

inflect <- c(1:length(sub.set))
for(i in 1:length(sub.set)){
  inflect[i] <- inf.repro1[sub.set[i]]
}

inflect

write.csv(inflect,"repro-success-infection-samples-CENTERED.csv")

##For replacement
sub.set1 <- sort(sample(1:length(inf.replace1), size=500))

inflect1 <- c(1:length(sub.set))
for(i in 1:length(sub.set)){
  inflect1[i] <- inf.replace1[sub.set[i]]
}

inflect1

write.csv(inflect1,"replace-success-infection-samples-CENTERED.csv")

cm<-c(40.45, 55.65)
probs<-c(0,1)
summary(lm(probs~cm))
plot(cm, probs)
abline(summary(lm(probs~cm)))
