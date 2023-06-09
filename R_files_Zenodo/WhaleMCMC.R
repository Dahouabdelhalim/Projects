##uses jags to estimate the functional response of humpback whales
##From Ferguson, J.M., Hopkins III, J.B, and Witteveen, B.H. (2018) "Integrating abundance and diet data to improve inferences of food web dynamics"
##Derivation of the model models used here are given in the main manusucript and in the Supplemental materials

library(runjags)
library(readxl)
library(dplyr)

Diets 		<- read_excel(path="AbundanceDiet_data.xlsx", sheet=2, col_names=TRUE, range="A3:G9")
Abundances 	<- read_excel(path="AbundanceDiet_data.xlsx", sheet=1, col_names=TRUE, range="A3:G9")
 

##Read in datasets
N.whale     <- pull(Abundances, var=2) #Whale abundances
SE.whale    <- pull(Abundances, var=3) #Whale abundances
N.fish    	<- pull(Abundances, var=4) #Fish density
SE.fish    	<- pull(Abundances, var=5) #Fish density
N.zoop    	<- pull(Abundances, var=6) #zooplankton density
SE.zoop    	<- pull(Abundances, var=7) #zooplankton density

DP    		<- pull(Diets, var=3)/100 #proportion of fish in the diet
N	  		<- pull(Diets, var=2) #sample size of diet data


#Abundance data model for humpback whales and their prey. Does not incorporate diet estimates
AbundanceModel <- "model {

	#generate unobserved true states
	logPtrue[1] 	~ dnorm(log(P[1]), P[1]^2/PSE[1]^2)
	logN1true[1] 	~ dnorm(log(N1obs[1]), N1obs[1]^2/N1SE[1]^2)
	logN2true[1] 	~ dnorm(log(N2obs[1]), N2obs[1]^2/N2SE[1]^2)
	
	Ptrue[1] 	<- exp(logPtrue[1])
	N1true[1] 	<- exp(logN1true[1])
	N2true[1] 	<- exp(logN2true[1])
	
	for(i in 2:N) {
		logPtrue[i] ~ dnorm(log(P[i]), P[i]^2/PSE[i]^2)
		Ptrue[i] 	<- exp(logPtrue[i])
	
		logN1true[i] ~ dnorm(log(N1true[i-1]) + r1*tstep[i] - exp(c1)*tstep[i]*Ptrue[i-1], procE1/(tstep[i]))
		N1true[i]  <- exp(logN1true[i])
		
		logN2true[i] ~ dnorm(log(N2true[i-1]) + r2*tstep[i] - exp(c2)*tstep[i]*Ptrue[i-1], procE2/(tstep[i]))
		N2true[i]  <- exp(logN2true[i])
	}
	
	#generate observed states
	for(i in 2:N) {
		N1obs[i] ~ dnorm(N1true[i], 1/N1SE[i]^2)
		N2obs[i] ~ dnorm(N2true[i], 1/N2SE[i]^2)
	}

	r1 ~ dexp(0.01) #dunif(0.01, 1)
	r2 ~ dexp(0.01) #dunif(0.01, 1)
	c1 ~ dunif(-40, 0) #dexp(1) #
	c2 ~ dunif(-40, 0) #dexp(1) #dunif(0,1)
	procE1 ~ dgamma(0.01, 0.01)
	procE2 ~ dgamma(0.01, 0.01)
}"


#Integrated data model for humpback whales and their prey. Incorporates both diet and abundance data.
IntegratedModel <- "model {

	#generate unobserved true states
	logPtrue[1] 	~ dnorm(log(P[1]), P[1]^2/PSE[1]^2)
	logN1true[1] 	~ dnorm(log(N1obs[1]), N1obs[1]^2/N1SE[1]^2)
	logN2true[1] 	~ dnorm(log(N2obs[1]), N2obs[1]^2/N2SE[1]^2)
	
	Ptrue[1] 	<- exp(logPtrue[1])
	N1true[1] 	<- exp(logN1true[1])
	N2true[1] 	<- exp(logN2true[1])
	
	for(i in 2:N) {
		logPtrue[i] ~ dnorm(log(P[i]), P[i]^2/PSE[i]^2)
		Ptrue[i] 	<- exp(logPtrue[i])
	
		logN1true[i] ~ dnorm(log(N1true[i-1]) + r1*tstep[i] - ec1*tstep[i]*Ptrue[i-1], procE1/(tstep[i]))
		N1true[i]  <- exp(logN1true[i])
		
		logN2true[i] ~ dnorm(log(N2true[i-1]) + r2*tstep[i] - ec2*tstep[i]*Ptrue[i-1], procE2/(tstep[i]))
		N2true[i]  <- exp(logN2true[i])
	}
	
	#generate observed states
	for(i in 2:N) {
		N1obs[i] ~ dnorm(N1true[i], 1/N1SE[i]^2)
		N2obs[i] ~ dnorm(N2true[i], 1/N2SE[i]^2)
	}

	r1 ~ dexp(0.01) #dgamma(0.1, 0.1) #dunif(0.01, 2)
	r2 ~ dexp(0.01) #dgamma(0.1, 0.1) #dunif(0.01, 2)
	c1 ~ dunif(-40, 0) #dexp(1) #
	c2 ~ dunif(-40, 0) #dexp(1) #dunif(0,1)
	ec1 <- exp(c1)
	ec2 <- exp(c2)
	
	procE1  ~ dgamma(0.1, 0.1)
	procE2  ~ dgamma(0.1, 0.1)
	tauDraw ~ dgamma(0.1, 0.1)
	
	phi[1] 		~ dnorm(0, tauDraw)
	mu[1] 		<- bm1*ec1*mean(Ptrue)*mean(N1true)/(bm1*ec1*mean(Ptrue)*mean(N1true) + bm2*ec2*mean(Ptrue)*mean(N2true))
	alpha[1] 	<- mu[1] * exp(phi[1])
  	beta[1]  	<- (1-mu[1]) * exp(phi[1])
  	DP[1] 		~ dbeta(alpha[1], beta[1])
  	
  	tMax <- 1

	for(i in 2:N) {

		phi[i] ~ dnorm(0, tauDraw)
		    	
    	x1[i] <- bm1*(ec1*exp(r1*(-1 + tstep[i]) - ec1*Ptrue[i-1]*tstep[i] + lambda - tstep[i]*lambda)*(-exp(r1) + exp(ec1*Ptrue[i] + lambda))*N1true[i]*Ptrue[i-1])/(ec1*Ptrue[i-1] - r1 + lambda)
    	x2[i] <- bm2*(ec2*exp(r2*(-1 + tstep[i]) - ec2*Ptrue[i-1]*tstep[i] + lambda - tstep[i]*lambda)*(-exp(r2) + exp(ec2*Ptrue[i] + lambda))*N2true[i]*Ptrue[i-1])/(ec2*Ptrue[i-1] - r2 + lambda)
    	
    	mu[i] <- x1[i]/(x1[i]+x2[i]) 
    	
  		alpha[i] 	<- mu[i] * exp(phi[i])
  		beta[i]  	<- (1-mu[i]) * exp(phi[i])
  		DP[i] 		~ dbeta(alpha[i], beta[i])
  		
    	
	}

}"

#define some parameters for the models
bm1 = 137*1*10^3 #biomass of a fish adjusted for digestability
bm2 = 15.5*0.93*10^-3
tstep <- c(NA, 1, 2, 5, 1, 1) #number of years between observations
turnover.rate <- log(2)/(7/365)

#data to pass into the models
datTable <- list(N1obs=N.fish, N2obs=N.zoop, P=N.whale, PSE=SE.whale, N1SE=SE.fish, N2SE=SE.zoop, DP=DP, DPlogit= log(DP/(1-DP)), N=length(tstep), tstep=tstep, lambda=turnover.rate, bm1=bm1, bm2=bm2)


#run the abundance state space model
par.est  <- c('r1', 'r2', 'c1', 'c2', 'procE1', 'procE2',' cadeviance') 
abund.ss  <- run.jags(AbundanceModel, burnin=5e4, sample=1e4, n.chains=4, thin=10, adapt=1e3, monitor=par.est, data=datTable, method="parallel")
print(abund.ss)

#run the integrated state space model
par.est  <- c('r1', 'r2', 'c1', 'c2', 'procE1', 'procE2', 'tauDraw', 'deviance') 
int.ss  <- run.jags(IntegratedModel, burnin=5e4, sample=1e4, n.chains=4, thin=10, adapt=1e3, monitor=par.est, data=datTable, method="parallel")

print(int.ss)


