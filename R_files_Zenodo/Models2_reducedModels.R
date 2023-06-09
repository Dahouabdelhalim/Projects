### load packages
library(MCMCglmm)  #also loads coda
library(parallel)
library(lattice)
library(tictoc)

### obtain the data
load("objects/cebus.RData")

cebus$ENSO_ <- 'Average/Neutral'
cebus$ENSO_[cebus$ENSO<=-0.5] <- 'Cool/La_Niña'
cebus$ENSO_[cebus$ENSO>=0.5] <- 'Warm/El_Niño'

### create sine and cosine of month to capture annual seasonality
require(plyr)
cebus$rad_12 <- mapvalues(
	cebus$Month, 
	from=c("01/Jan", "02/Feb", "03/Mar", "04/Apr", "05/May",
				 "06/Jun", "07/Jul", "08/Aug", "09/Sep", "10/Oct",
				 "11/Nov", "12/Dec"), 
	to=seq(12)*(2*pi/12))
cebus$rad_12 <- as.numeric(cebus$rad_12)

### format factor variables
cols <- c("id", "animal", "Mother", "Group", "Sex", "Year", "Month", "ENSO_")
cebus[cols] <- lapply(cebus[cols], factor)

### Pedigree
load("objects/cebusPed.RData")
cols <- c("animal", "Mother", "Father")
Ped[cols] <- lapply(Ped[cols], factor)

rm(cols)

###############################################################################
###############################################################################
###############################################################################

### multinomial2 models ###

### NOTE: to set the prior
# one G structure needs to be specified for each random effect
# one R structure (for residuals) needs to be specified; bivariate model has two, univariate has one

### Parameter expanded prior ###
### taken from Tutorial by Pierre de Villemereuil
### Estimation of a biological trait Heritability using the animal model 
#### and MCMCglmm (version 2), 2021-09-22

prior_m9a <- list(
	R=list(V=1, nu=0.002),  
	G=list(
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),  
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))

### Intercept only

set.seed(1980)  #The Empire Strikes Back
tic() 
m09a.0 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 1,  #intercept only
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 1000, burnin = 50000, nitt = 2050000,
		prior = prior_m9a)
}, mc.cores=4)
toc()  #119562.022 sec elapsed

### Age

set.seed(1983)  #Return of the Jedi
tic() 
m09a.1 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE),  
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000,
		prior = prior_m9a)
}, mc.cores=4)
toc()  #146878.393 sec elapsed

### Age*Sex

set.seed(2015)  #The Force Awakens
tic() 
m09a.2 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex, 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000,
		prior = prior_m9a)
}, mc.cores=4)
toc()  #56076.345 sec elapsed

### Age*Sex + group size

set.seed(2017)  #The Last Jedi
tic() 
m09a.3 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000,
		prior = prior_m9a)
}, mc.cores=4)
toc()  #61084.651 sec elapsed

### Age*Sex + group size + seasonality

set.seed(2019)  #The Rise of Skywalker
tic() 
m09a.4 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000,
		prior = prior_m9a)
}, mc.cores=4)
toc()  #65323.258 sec elapsed

###############################################################################
###############################################################################
###############################################################################

## weakly informative inverse-Wishart prior (V=1, nu=0.002)
## performs poorly if variance is close to zero
prior_m9b <- list(R=list(V=1, nu=0.002),  
									G=list(G1=list(V=1, nu=0.002),  
												 G2=list(V=1, nu=0.002),
												 G3=list(V=1, nu=0.002),
												 G4=list(V=1, nu=0.002),
												 G5=list(V=1, nu=0.002),
												 G6=list(V=1, nu=0.002),
												 G7=list(V=1, nu=0.002),
												 G8=list(V=1, nu=0.002),
												 G9=list(V=1, nu=0.002)))

### Intercept only

set.seed(1980)  #The Empire Strikes Back
tic() 
m09b.0 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 1,  #intercept only
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9b)
}, mc.cores=4)
toc()  #50706.757 sec elapsed

### Age

set.seed(1983)  #Return of the Jedi
tic() 
m09b.1 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE),  
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9b)
}, mc.cores=4)
toc()  #55679.073 sec elapsed

### Age*Sex

set.seed(2015)  #The Force Awakens
tic() 
m09b.2 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex, 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9b)
}, mc.cores=4)
toc()  #52126.754 sec elapsed

### Age*Sex + group size

set.seed(2017)  #The Last Jedi
tic() 
m09b.3 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9b)
}, mc.cores=4)
toc()  #60529.528 sec elapsed

### Age*Sex + group size + seasonality

set.seed(2019)  #The Rise of Skywalker
tic() 
m09b.4 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9b)
}, mc.cores=4)
toc()  #65323.258 sec elapsed

###############################################################################
###############################################################################
###############################################################################

### Inverse Wishart prior, but slightly stronger belief parameter
prior_m9c <- list(R=list(V=1, nu=0.02),  
									G=list(G1=list(V=1, nu=0.02),  
												 G2=list(V=1, nu=0.02),
												 G3=list(V=1, nu=0.02),
												 G4=list(V=1, nu=0.02),
												 G5=list(V=1, nu=0.02),
												 G6=list(V=1, nu=0.02),
												 G7=list(V=1, nu=0.02),
												 G8=list(V=1, nu=0.02),
												 G9=list(V=1, nu=0.02)))

### Intercept only

set.seed(1980)  #The Empire Strikes Back
tic() 
m09c.0 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 1,  #intercept only
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9c)
}, mc.cores=4)
toc()  #49074.051 sec elapsed

### Age

set.seed(1983)  #Return of the Jedi
tic() 
m09c.1 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE),  
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9c)
}, mc.cores=4)
toc()  #57637.06 sec elapsed

### Age*Sex

set.seed(2015)  #The Force Awakens
tic() 
m09c.2 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex, 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9c)
}, mc.cores=4)
toc()  #64051.073 sec elapsed

### Age*Sex + group size

set.seed(2017)  #The Last Jedi
tic() 
m09c.3 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9c)
}, mc.cores=4)
toc()  #66759.551 sec elapsed

### Age*Sex + group size + seasonality

set.seed(2019)  #The Rise of Skywalker
tic() 
m09c.4 <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social, n-social) ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "multinomial2", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 1050000,
		prior = prior_m9c)
}, mc.cores=4)
toc()  #60832.799 sec elapsed

###############################################################################
###############################################################################
###############################################################################

## Priors for poisson models ##
###############################

### Parameter expanded prior ###
### taken from Tutorial by Pierre de Villemereuil
### Estimation of a biological trait Heritability using the animal model 
#### and MCMCglmm (version 2), 2021-09-22

k <- 2 # intercept plus log(n)
prior_p9a.0 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),  
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
prior_p9a.0$B$mu[k] <- 1 
prior_p9a.0$B$V[k,k] <- 1e-7

k <- 5 # intercept, log(n), age
prior_p9a.1 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),  
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
prior_p9a.1$B$mu[k] <- 1 
prior_p9a.1$B$V[k,k] <- 1e-7

k <- 9 # intercept, log(n), age*sex
prior_p9a.2 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),  
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
prior_p9a.2$B$mu[k-3] <- 1 
prior_p9a.2$B$V[k-3,k-3] <- 1e-7

k <- 10 # intercept, log(n), age*sex, grp_size
prior_p9a.3 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),  
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
prior_p9a.3$B$mu[k-3] <- 1 
prior_p9a.3$B$V[k-3,k-3] <- 1e-7

k <- 12 # number of fixed effects plus intercept
prior_p9a.4 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),  
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G3=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G4=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G5=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G6=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G7=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G8=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G9=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
prior_p9a.4$B$mu[k-3] <- 1 
prior_p9a.4$B$V[k-3,k-3] <- 1e-7

### Intercept only

set.seed(1982)  #Blade Runner :) 
tic() 
p09a.0 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			log(n),  #intercept and log(n) only
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000, 
		prior = prior_p9a.0)
}, mc.cores=4)
toc()  #140260.496 sec elapsed

### Age

set.seed(1982)  #Blade Runner :) 
tic() 
p09a.1 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE) + 
			log(n),  
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000, 
		prior = prior_p9a.1)
}, mc.cores=4)
toc()  #153117.004 sec elapsed

### Age*Sex

set.seed(1982)  #Blade Runner :) 
tic() 
p09a.2 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000, 
		prior = prior_p9a.2)
}, mc.cores=4)
toc()  #161707.664 sec elapsed

### Age*Sex + group size

set.seed(1982)  #Blade Runner :) 
tic() 
p09a.3 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000, 
		prior = prior_p9a.3)
}, mc.cores=4)
toc()  #169555.753 sec elapsed

### Age*Sex + group size + seasonality

set.seed(1982)  #Blade Runner :) 
tic() 
p09a.4 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 500, burnin = 50000, nitt = 2050000, 
		prior = prior_p9a.4)
}, mc.cores=4)
toc()  #198792.738 sec elapsed

###############################################################################
###############################################################################
###############################################################################

## Inverse Wishart prior ##
###############################

k <- 2 # intercept plus log(n)
prior_p9b.0 <- list(B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
										R=list(V=1, nu=0.002),  
										G=list(
											G1=list(V=1, nu=0.002),  
											G2=list(V=1, nu=0.002),
											G3=list(V=1, nu=0.002),
											G4=list(V=1, nu=0.002),
											G5=list(V=1, nu=0.002),
											G6=list(V=1, nu=0.002),
											G7=list(V=1, nu=0.002),
											G8=list(V=1, nu=0.002),
											G9=list(V=1, nu=0.002)))
prior_p9b.0$B$mu[k-3] <- 1 
prior_p9b.0$B$V[k-3,k-3] <- 1e-7

k <- 5 # intercept, log(n), age
prior_p9b.1 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V=1, nu=0.002),  
		G2=list(V=1, nu=0.002),
		G3=list(V=1, nu=0.002),
		G4=list(V=1, nu=0.002),
		G5=list(V=1, nu=0.002),
		G6=list(V=1, nu=0.002),
		G7=list(V=1, nu=0.002),
		G8=list(V=1, nu=0.002),
		G9=list(V=1, nu=0.002)))
prior_p9b.1$B$mu[k] <- 1 
prior_p9b.1$B$V[k,k] <- 1e-7

k <- 9 # intercept, log(n), age*sex
prior_p9b.2 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V=1, nu=0.002),  
		G2=list(V=1, nu=0.002),
		G3=list(V=1, nu=0.002),
		G4=list(V=1, nu=0.002),
		G5=list(V=1, nu=0.002),
		G6=list(V=1, nu=0.002),
		G7=list(V=1, nu=0.002),
		G8=list(V=1, nu=0.002),
		G9=list(V=1, nu=0.002)))
prior_p9b.2$B$mu[k-3] <- 1 
prior_p9b.2$B$V[k-3,k-3] <- 1e-7

k <- 10 # intercept, log(n), age*sex, grp_size
prior_p9b.3 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V=1, nu=0.002),  
		G2=list(V=1, nu=0.002),
		G3=list(V=1, nu=0.002),
		G4=list(V=1, nu=0.002),
		G5=list(V=1, nu=0.002),
		G6=list(V=1, nu=0.002),
		G7=list(V=1, nu=0.002),
		G8=list(V=1, nu=0.002),
		G9=list(V=1, nu=0.002)))
prior_p9b.3$B$mu[k-3] <- 1 
prior_p9b.3$B$V[k-3,k-3] <- 1e-7

k <- 12 # number of fixed effects plus intercept
prior_p9b.4 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V=1, nu=0.002),  
		G2=list(V=1, nu=0.002),
		G3=list(V=1, nu=0.002),
		G4=list(V=1, nu=0.002),
		G5=list(V=1, nu=0.002),
		G6=list(V=1, nu=0.002),
		G7=list(V=1, nu=0.002),
		G8=list(V=1, nu=0.002),
		G9=list(V=1, nu=0.002)))
prior_p9b.4$B$mu[k-3] <- 1 
prior_p9b.4$B$V[k-3,k-3] <- 1e-7

### Intercept only

set.seed(1982)  #Blade Runner :) 
tic() 
p09b.0 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			log(n),  #intercept and log(n) only
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9b.0)
}, mc.cores=4)
toc()  #29773.753 sec elapsed

### Age

set.seed(1982)  #Blade Runner :) 
tic() 
p09b.1 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE) + 
			log(n),  
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9b.1)
}, mc.cores=4)
toc()  #34112.004 sec elapsed

### Age*Sex

set.seed(1982)  #Blade Runner :) 
tic() 
p09b.2 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9b.2)
}, mc.cores=4)
toc()  #27160.523 sec elapsed

### Age*Sex + group size

set.seed(1982)  #Blade Runner :) 
tic() 
p09b.3 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9b.3)
}, mc.cores=4)
toc()  #28120.869 sec elapsed

### Age*Sex + group size + seasonality

set.seed(1982)  #Blade Runner :) 
tic() 
p09b.4 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9b.4)
}, mc.cores=4)
toc()  #38543.441 sec elapsed

###############################################################################
###############################################################################
###############################################################################

## Inverse Wishart priors, but slightly stronger belief parameter ##

k <- 2 # intercept plus log(n)
prior_p9c.0 <- list(B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
										R=list(V=1, nu=1),  
										G=list(
											G1=list(V=1, nu=0.02),  
											G2=list(V=1, nu=0.02),
											G3=list(V=1, nu=0.02),
											G4=list(V=1, nu=0.02),
											G5=list(V=1, nu=0.02),
											G6=list(V=1, nu=0.02),
											G7=list(V=1, nu=0.02),
											G8=list(V=1, nu=0.02),
											G9=list(V=1, nu=0.02)))
prior_p9c.0$B$mu[k-3] <- 1 
prior_p9c.0$B$V[k-3,k-3] <- 1e-7

k <- 5 # intercept, log(n), age
prior_p9c.1 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V=1, nu=0.02),  
		G2=list(V=1, nu=0.02),
		G3=list(V=1, nu=0.02),
		G4=list(V=1, nu=0.02),
		G5=list(V=1, nu=0.02),
		G6=list(V=1, nu=0.02),
		G7=list(V=1, nu=0.02),
		G8=list(V=1, nu=0.02),
		G9=list(V=1, nu=0.02)))
prior_p9c.1$B$mu[k] <- 1 
prior_p9c.1$B$V[k,k] <- 1e-7

k <- 9 # intercept, log(n), age*sex
prior_p9c.2 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V=1, nu=0.02),  
		G2=list(V=1, nu=0.02),
		G3=list(V=1, nu=0.02),
		G4=list(V=1, nu=0.02),
		G5=list(V=1, nu=0.02),
		G6=list(V=1, nu=0.02),
		G7=list(V=1, nu=0.02),
		G8=list(V=1, nu=0.02),
		G9=list(V=1, nu=0.02)))
prior_p9c.2$B$mu[k-3] <- 1 
prior_p9c.2$B$V[k-3,k-3] <- 1e-7

k <- 10 # intercept, log(n), age*sex, grp_size
prior_p9c.3 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V=1, nu=0.02),  
		G2=list(V=1, nu=0.02),
		G3=list(V=1, nu=0.02),
		G4=list(V=1, nu=0.02),
		G5=list(V=1, nu=0.02),
		G6=list(V=1, nu=0.02),
		G7=list(V=1, nu=0.02),
		G8=list(V=1, nu=0.02),
		G9=list(V=1, nu=0.02)))
prior_p9c.3$B$mu[k-3] <- 1 
prior_p9c.3$B$V[k-3,k-3] <- 1e-7

k <- 12 # number of fixed effects plus intercept
prior_p9c.4 <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(
		G1=list(V=1, nu=0.02),  
		G2=list(V=1, nu=0.02),
		G3=list(V=1, nu=0.02),
		G4=list(V=1, nu=0.02),
		G5=list(V=1, nu=0.02),
		G6=list(V=1, nu=0.02),
		G7=list(V=1, nu=0.02),
		G8=list(V=1, nu=0.02),
		G9=list(V=1, nu=0.02)))
prior_p9c.4$B$mu[k-3] <- 1 
prior_p9c.4$B$V[k-3,k-3] <- 1e-7

### Intercept only

set.seed(1982)  #Blade Runner :) 
tic() 
p09c.0 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			log(n),  #intercept and log(n) only
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9c.0)
}, mc.cores=4)
toc()  #29522.443 sec elapsed

### Age

set.seed(1982)  #Blade Runner :) 
tic() 
p09c.1 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE) + 
			log(n),  
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9c.1)
}, mc.cores=4)
toc()  #33640.63 sec elapsed

### Age*Sex

set.seed(1982)  #Blade Runner :) 
tic() 
p09c.2 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9c.2)
}, mc.cores=4)
toc()  #36266.344 sec elapsed

### Age*Sex + group size

set.seed(1982)  #Blade Runner :) 
tic() 
p09c.3 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9c.3)
}, mc.cores=4)
toc()  #38421.41 sec elapsed

### Age*Sex + group size + seasonality

set.seed(1982)  #Blade Runner :) 
tic() 
p09c.4 <- mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(scale(age), degree = 3, raw = TRUE)*Sex + 
			scale(grp_size) + 
			sin(rad_12) + cos(rad_12) + 
			log(n), 
		random = ~ animal + id + id:Year + 
			Mother + Mother:Group:Year +  
			Group + Group:Year + 
			Year + Month:Year,
		family = "poisson", 
		data = cebus, pedigree = Ped, 
		thin = 250, burnin = 50000, nitt = 550000, 
		prior = prior_p9c.4)
}, mc.cores=4)
toc()  #41540.663 sec elapsed

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# save all models objects (reduced models)
save.image("Models2_reducedModels.RData")
