### load packages
library(MCMCglmm)
library(parallel)
library(ggplot2)
library(ggpubr)
library(tidyverse)

### can load previously run models
#load("Models_Preliminary.RData")

### load data
load("cebus.RData")

### create variable for El Ni単o/La Ni単a southern oscillation patterns
cebus$ENSO_ <- 'Average/Neutral'
cebus$ENSO_[cebus$ENSO<=-0.5] <- 'Cool/La_Ni単a'
cebus$ENSO_[cebus$ENSO>=0.5] <- 'Warm/El_Ni単o'

### create sine and cosine of month to capture annual seasonality
require(plyr)
cebus$rad_12 <- mapvalues(cebus$Month, 
												 from=c("01/Jan", "02/Feb", "03/Mar", "04/Apr", "05/May",
												 			 "06/Jun", "07/Jul", "08/Aug", "09/Sep", "10/Oct",
												 			 "11/Nov", "12/Dec"), 
												 to=seq(12)*(2*pi/12))
cebus$rad_12 <- as.numeric(cebus$rad_12)
cebus$sin_12 <- sin(cebus$rad_12)
cebus$cos_12 <- cos(cebus$rad_12)

### format factor variables
cols <- c("id", "Sex", "Month", "ENSO_")
cebus[cols] <- lapply(cebus[cols], factor)

rm(cols)

##############################################################
##############################################################
##############################################################

### NOTE: to set the prior
# one G structure needs to be specified for each random effect
# one R structure needs to be specified, bivariate model has two, univariate has one

### parameters
# default: thin=10, burnin=3000
# reduce autocorrelation (below 0.1), effective sample size of at least 1000

### multinomial2 priors, parameter-expanded

p1a = list(
	R=list(V=1, nu=0.02), 
	G=list(
		G1=list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

p2a = list(
	R=list(V=1, nu=0.02), 
	G=list(
		G1=list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1),
		G1=list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1)))

### poisson priors, parameter-expanded

k <- 23 # number of fixed effects plus intercept
p1b <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(  #prior for random effect
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
p1b$B$mu[k-3] <- 1 
p1b$B$V[k-3,k-3] <- 1e-7

k <- 14 # number of fixed effects plus intercept
p2b <- list(
	B=list(V=diag(k)*1e7, mu=rep(0,k)),  #priors for fixed effects
	R=list(V=1, nu=1),  #prior for response
	G=list(  #prior for random effect
		G1=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
		G2=list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
p2b$B$mu[k-3] <- 1 
p2b$B$V[k-3,k-3] <- 1e-7

################################################################################
################################################################################

## Random: id
set.seed(5)
season1a <-  mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social,n-social) ~ 
			Month + 
			poly(z_age, degree=3, raw=TRUE)*Sex + 
			z_grp_size + ENSO_,
		random = ~id,
		thin = 100, burnin = 10000, nitt = 210000, 
		family = "multinomial2", data = cebus,
		prior = p1a, verbose = F)
}, mc.cores=4)

mod <- season1a[[1]]
heidel.diag(mod$VCV) # diagnostic test of convergence
# Stationarity start     p-value
# test         iteration        
# id    passed       1         0.892  
# units passed       1         0.453  
# 
# Halfwidth Mean  Halfwidth
# test                     
# id    passed    0.133 0.000547 
# units passed    0.182 0.000168 
autocorr.diag(mod$VCV)
#                   id      units
# Lag 0    1.000000000 1.00000000
# Lag 100  0.044866578 0.02997358
# Lag 500  0.005303712 0.01362369
# Lag 1000 0.003575660 0.01152172
# Lag 5000 0.012739117 0.00481630
plot(mod$VCV)
summary(mod)

set.seed(5)
season1b <-  mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			Month + 
			poly(z_age, degree=3, raw=TRUE)*Sex + 
			z_grp_size + ENSO_ + log(n),
		random = ~id,
		thin = 20, burnin = 10000, nitt = 50000, 
		family = "poisson", data = cebus,
		prior = p1b, verbose = F)
}, mc.cores=4)

mod <- season1b[[1]]
heidel.diag(mod$VCV) # diagnostic test of convergence
# Stationarity start     p-value
# test         iteration        
# id    passed       1         0.144  
# units passed       1         0.132  
# 
# Halfwidth Mean   Halfwidth
# test                      
# id    passed    0.0716 2.58e-04 
# units passed    0.1249 7.89e-05 
autocorr.diag(mod$VCV)
#                   id        units
# Lag 0     1.00000000  1.000000000
# Lag 20    0.01047829  0.014244778
# Lag 100   0.01370727 -0.014637687
# Lag 200  -0.02547719 -0.039401877
# Lag 1000  0.01741752  0.007793999
plot(mod$VCV)
summary(mod)

################################################################################

## Random: id, Month (with sinewave)
set.seed(5)
season2a <-  mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social,n-social) ~ 
			sin_12 + cos_12 + 
			poly(z_age, degree=3, raw=TRUE)*Sex + z_grp_size + ENSO_,
		random = ~id + Month,
		family = "multinomial2", data = cebus,
		thin = 100, burnin = 10000, nitt = 210000, 
		prior = p2a, verbose = F)
}, mc.cores=4)

mod <- season2a[[1]]
heidel.diag(mod$VCV) # diagnostic test of convergence
# Stationarity start     p-value
# test         iteration        
# id    passed       1         0.963  
# Month passed       1         0.776  
# units passed       1         0.968  
# 
# Halfwidth Mean    Halfwidth
# test                       
# id    passed    0.13364 0.000486 
# Month passed    0.00364 0.000104 
# units passed    0.18192 0.000165 
autocorr.diag(mod$VCV)
#                     id        Month        units
# Lag 0     1.0000000000  1.000000000  1.000000000
# Lag 100   0.0263731153  0.005283400 -0.008434326
# Lag 500  -0.0147736857  0.008680901 -0.013897594
# Lag 1000  0.0003535233  0.008077985 -0.031072076
# Lag 5000 -0.0302676205 -0.015296557 -0.004293213
plot(mod$VCV)
summary(mod)

set.seed(5)
season2b <-  mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			sin_12 + cos_12 + 
			poly(z_age, degree=3, raw=TRUE)*Sex + z_grp_size + ENSO_ + log(n),
		random = ~id + Month,
		family = "poisson", data = cebus,
		thin = 100, burnin = 10000, nitt = 210000, 
		prior = p2b, verbose = F)
}, mc.cores=4)

mod <- season2b[[1]]
heidel.diag(mod$VCV) # diagnostic test of convergence
# Stationarity start     p-value
# test         iteration        
# id    passed       401       0.0749 
# Month passed         1       0.7824 
# units passed         1       0.1875 
# 
# Halfwidth Mean    Halfwidth
# test                       
# id    passed    0.07133 2.96e-04 
# Month passed    0.00298 1.16e-04 
# units passed    0.12501 8.15e-05 
autocorr.diag(mod$VCV)
# id      Month       units
# Lag 0     1.000000000 1.00000000 1.000000000
# Lag 100   0.020824608 0.09509899 0.029953906
# Lag 500   0.014151958 0.01586306 0.013150900
# Lag 1000 -0.018849478 0.01046619 0.004003934
# Lag 5000  0.002932262 0.01184480 0.029839231
plot(mod$VCV)
summary(mod)

################################################################################

## Random: id, Month (without sinewave, no Month as fixed effect)
set.seed(5)
season3a <- mclapply(1:4, function(i) {
	MCMCglmm(
		cbind(social,n-social) ~ 
			poly(z_age, degree=3, raw=TRUE)*Sex + 
			z_grp_size + ENSO_,
		random = ~id + Month,
		family = "multinomial2", data = cebus,
		thin = 100, burnin = 10000, nitt = 210000, 
		prior = p2a, verbose = F)
}, mc.cores=4)

mod <- season3a[[1]]
heidel.diag(mod$VCV) # diagnostic test of convergence
autocorr.diag(mod$VCV)
plot(mod$VCV)
summary(mod)

set.seed(5)
season3b <-  mclapply(1:4, function(i) {
	MCMCglmm(
		partners ~ 
			poly(z_age, degree=3, raw=TRUE)*Sex + 
			z_grp_size + ENSO_ + log(n),
		random = ~id + Month,
		family = "multinomial2", data = cebus,
		thin = 20, burnin = 10000, nitt = 50000, 
		prior = p2c, verbose = F)
}, mc.cores=4)

mod <- season3b[[1]]
heidel.diag(mod$VCV) # diagnostic test of convergence
autocorr.diag(mod$VCV)
plot(mod$VCV)
summary(mod)

################################################################################

## variance explained by Month (with sine wave)
vmonth_lat1a <- season2a[[1]][["VCV"]][,"Month"]/rowSums(season2a[[1]][["VCV"]])
posterior.mode(vmonth_lat1a)  #0.006898364 
HPDinterval(vmonth_lat1a)
#            lower      upper
# var1 0.002448188 0.02459136
# attr(,"Probability")
# [1] 0.95

## variance explained by Month (without sine wave)
vmonth_lat2a <- season3a[[1]][["VCV"]][,"Month"]/rowSums(season3a[[1]][["VCV"]])
posterior.mode(vmonth_lat2a)  #0.0733448 
HPDinterval(vmonth_lat2a)
#           lower     upper
# var1 0.03833803 0.1914481
# attr(,"Probability")
# [1] 0.95

## variance explained by Month (with sine wave)
vmonth_lat1b <- season2b[[1]][["VCV"]][,"Month"]/rowSums(season2b[[1]][["VCV"]])
posterior.mode(vmonth_lat1b)  # 0.009433306 
HPDinterval(vmonth_lat1b)
#            lower      upper
# var1 0.003659318 0.03436415
# attr(,"Probability")
# [1] 0.95

## variance explained by Month (without sine wave)
vmonth_lat2a <- season3b[[1]][["VCV"]][,"Month"]/rowSums(season3b[[1]][["VCV"]])
posterior.mode(vmonth_lat2b)  #0.0733448 
HPDinterval(vmonth_lat2b)

################################################################################

### SI Figure 1 ###
###################

plotFE <- function(mod, title){
	tmp <- as.data.frame(summary(mod)$solutions)
	tmp$effects <- factor(rownames(tmp),levels=rownames(tmp))
	colnames(tmp)[2:3] <- c("lowerCI","upperCI")
	p <- ggplot(tmp[-1,], aes(y=post.mean, x=effects)) +
		geom_point(color="black") +
		geom_linerange(aes(ymin=lowerCI, ymax=upperCI)) +
		geom_hline(yintercept = 0, lty=3) +
		scale_x_discrete(limits = rev(levels(tmp$effects)[-1]), 
										 labels = function(effects) str_wrap(effects, width = 24)) + 
		ylab("") + xlab("") + labs(title=title) + coord_flip() + theme_light()
	p
}
fixed.p1a <- plotFE(season1a[[1]], "a) Preliminary Model 1")
fixed.p2a <- plotFE(season2a[[1]], "b) Preliminary Model 2")
fixed.p1b <- plotFE(season1b[[1]], "a) Preliminary Model 1")
fixed.p2b <- plotFE(season2b[[1]], "b) Preliminary Model 2")

pdf("supp/Supp_Fig1_FE_seasonality1.pdf")
ggarrange(fixed.p1a, fixed.p2a, 
					ncol = 2, nrow = 1)
dev.off()

pdf("supp/Supp_Fig2_FE_seasonality2.pdf")
ggarrange(fixed.p1b, fixed.p2b, 
					ncol = 2, nrow = 1)
dev.off()

rm(fixed.p1a, fixed.p2a, fixed.p1b, fixed.p2b, mod, p1a, p2a, p1b, p2b)
rm(mod, k, plotFE, vmonth_lat1a, vmonth_lat1b, vmonth_lat2a)

save.image("Models_Preliminary.RData")

################################################################################
################################################################################
