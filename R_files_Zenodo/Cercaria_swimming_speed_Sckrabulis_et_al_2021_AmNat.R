# Cercaria_swimming_speed_Sckrabulis_et_al_2021_AmNat.R

# Cercaria swimming speed analysis 

#####
# Libraries
library(nls.multstart)
library(nlstools)
library(modelr)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(Hmisc)

# Constants
k <- 8.62*10^-5
K <- 273.15

# Use 'Cercaria_swimming_speed_Sckrabulis_et_al_2021_AmNat.csv' in this repository
swimData <- read.csv(file.choose(), header=TRUE)
swimData$PerfTempK <- swimData$Temperature+K
swimData$logavgSpeedMM <- log(swimData$avgSpeedMM+1)
plot(swimData$Temperature, swimData$avgSpeedMM)

#####
# Xiao et al. 2011 AIC calculation
# Reference: Xiao, X., E. P. White, M. B. Hooten, and S. L. Durham. 2011. On the use of log-transformation vs. nonlinear regression for analyzing biological power laws. Ecology 92: 1887-1894. doi:10.1890/11-0538.1
AIC.Xiao <- function(numParams, loglik, n){
	2*numParams-2*loglik
}

#####
# Testing for best fit error distribution

#Boltzmann-Arrhenius equation
swimBA <- function(x, iTo, Ei, Toi){
	Toi <- Toi+K
	y <- iTo*exp(-Ei/k*(1/(x+K)-1/Toi))
}

# Model with normal error distribution
BAnorm <- nls(avgSpeedMM ~ swimBA(x=Temperature, iTo, Ei, Toi=19),	start=c(iTo=1, Ei=.6), data=swimData, control=c(warnOnly=TRUE))

# Plotting function using estimates from 'BAnorm' model
normPlotBA <- function(x){
	Toi <- 19+K
	iTo <- 2.6159
	Ei <- 0.09557
	y <- iTo*exp(-Ei/k*(1/(x+K)-1/Toi))
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd((swimData$avgSpeedMM)-normPlotBA(swimData$Temperature))
ll_norm <- sum(log(dnorm(swimData$avgSpeedMM, normPlotBA(swimData$Temperature), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=2, loglik=ll_norm, n=length(swimData$Temperature))
normAIC

# Model with lognormal error error distribution
BAlnorm <- nls(log(avgSpeedMM+1) ~ log(swimBA(x=Temperature, iTo, Ei, Toi=19)), start=c(iTo=1, Ei=.6), data=swimData, control=c(warnOnly=TRUE))

# Plotting function using estimates from 'BAlnorm' model
lnormPlotBA <- function(x){
	Toi <- 19+K
	iTo <- 3.3505
	Ei <- .05907
	y <- log(iTo*exp(-Ei/k*(1/(x+K)-1/Toi)))
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd((swimData$logavgSpeedMM)-lnormPlotBA(swimData$Temperature))
ll_lnorm <- sum(log(dlnorm(swimData$avgSpeedMM+1, lnormPlotBA(swimData$Temperature), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=2, loglik=ll_lnorm, n=length(swimData$Temperature))
lnormAIC

#__
# Sharpe-Schoolfield equation
swimSS <- function(x, iTo, Ei, Ehi, Thi, Toi){
	Toi <- Toi+K
	y <- iTo*exp(-Ei/k*(1/(x+K)-1/Toi))*(1+exp(Ehi/k*(1/Thi-1/(x+K))))^-1
}

# Model with normal error distribution
SSnorm <- nls(avgSpeedMM ~ swimSS(x=Temperature, iTo, Ei, Ehi, Thi, Toi=19), start=c(iTo=2, Ei=.6, Ehi=3, Thi=303), data=swimData, control=c(warnOnly=TRUE))

# Plotting function using estimates from 'SSnorm' model
normPlotSS <- function(x){
	Toi <- 19+K
	iTo <- 2.63454
	Ei <- .52545
	Ehi <- 4.24226
	Thi <- 304.05194
	y <- iTo*exp(-Ei/k*(1/(x+K)-1/Toi))*(1+exp(Ehi/k*(1/Thi-1/(x+K))))^-1
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd((swimData$avgSpeedMM)-normPlotSS(swimData$Temperature))
ll_norm <- sum(log(dnorm(swimData$avgSpeedMM, normPlotSS(swimData$Temperature), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=4, loglik=ll_norm, n=length(SwimData$Temperature))
normAIC

#Model with lognormal error distribution
SSlnorm <- nls(log(avgSpeedMM+1) ~ log(swimSS(x=Temperature, iTo, Ei, Ehi, Thi, Toi=19)), start=c(iTo=2, Ei=.6, Ehi=3, Thi=303), data=swimData, control=c(warnOnly=TRUE))

# Plotting function using estimates from 'SSlnorm' model
lnormPlotSS <- function(x){
	Toi <- 19+K
	iTo <- 3.5399
	Ei <- .4073
	Ehi <- 3.7302
	Thi <- 304.7443
	y <- log(iTo*exp(-Ei/k*(1/(x+K)-1/Toi))*(1+exp(Ehi/k*(1/Thi-1/(x+K))))^-1)
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd((swimData$logavgSpeedMM)-lnormPlotSS(swimData$Temperature))
ll_lnorm <- sum(log(dlnorm(swimData$avgSpeedMM+1, lnormPlotSS(swimData$Temperature), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=4, loglik=ll_lnorm, n=length(swimData$Temperature))
lnormAIC

#####
# Using nls.multstart package to generate confidence bands (Padfield & Matheson 2018)
# Reference: Padfield, D., and G. Matheson. 2018. Nls.multstart: An R package for robust non-linear regression using AIC scores. v1.0.0
swimMult <- function(lniTo, Ei, Ehi, Thi, temp, Toi) {
	Toi <- Toi+K
	boltzmann.term <- lniTo + log(exp(Ei/k*(1/Toi - 1/temp)))
	inactivation.term <- log(1/(1 + exp(Ehi/k*(1/Thi - 1/temp))))
	return(boltzmann.term + inactivation.term)
}

# One fit
fit <- nls_multstart(avgSpeedMM ~ swimMult(lniTo, Ei, Ehi, Thi, temp = PerfTempK, Toi = 19),
	data = swimData,iter = 500,
	start_lower = c(lniTo = -10, Ei = 0.1, Ehi = 0.2, Thi = 285),
	start_upper = c(lniTo = 10, Ei = 2, Ehi = 5, Thi = 330),
	supp_errors = 'Y',na.action = na.omit,
	lower = c(lniTo = -10, Ei = 0, Ehi = 0, Thi = 0))

params <- tidy(fit)
preds <- augment(fit)

# Bootstrapped fits
fit_boots <- swimData %>% 
	modelr::bootstrap(n = 200, id = 'boot_num') %>%
	group_by(boot_num) %>%
	mutate(fit = map(strap, ~nls_multstart(log(avgSpeedMM+1) ~ swimMult(lniTo, Ei, Ehi, Thi, temp = PerfTempK, Toi = 19),
		data = data.frame(.),iter = 100,
		start_lower = c(lniTo = -10, Ei = 0.1, Ehi = 0.2, Thi = 285),
		start_upper = c(lniTo = 10, Ei = 2, Ehi = 5, Thi = 330),
		lower = c(lniTo=-10, Ei=0, Ehi=0, Thi=0),
		supp_errors = 'Y')
	))
fit_boots

# Original example from Padfield & Matheson 2018 uses 'unnest()', but this function does not function as inteded with recent update to tidyr package. 'unnest_legacy()' is included in current versions that operates as intended 
params_boot <- fit_boots %>%
	unnest_legacy(fit %>% map(tidy)) %>%
	ungroup()
preds_boot <- fit_boots %>%
	unnest_legacy(fit %>% map(augment)) %>%
	ungroup()
  
# Create new data frame of predictions
new_preds <- swimData %>%
	do(., data.frame(PerfTempK = seq(min(.$PerfTempK), max(.$PerfTempK), length.out = 250), stringsAsFactors = FALSE))

# Create smoother predictions for best fit model
preds <- augment(fit, newdata = new_preds)

# Create smoother predictions for bootstrapped replicates
preds <- fit_boots %>%
	unnest_legacy(fit %>% map(augment, newdata = new_preds)) %>%
	# Group by each value of PerfTempK and get quantiles
	group_by(., PerfTempK) %>%
		summarise(lwr_CI = quantile(.fitted, 0.025),
		upr_CI = quantile(.fitted, 0.975)) %>%
	ungroup() %>%
	merge(., preds, by = 'PerfTempK')

#####
# Figure S4
snail1 <- subset(swimData, Snail==1)
temps1 <- c(13, 16, 19, 22, 25, 28, 31)
s1_13 <- subset(snail1, Temperature==13)
s1_16 <- subset(snail1, Temperature==16)
s1_19 <- subset(snail1, Temperature==19)
s1_22 <- subset(snail1, Temperature==22)
s1_25 <- subset(snail1, Temperature==25)
s1_28 <- subset(snail1, Temperature==28)
s1_31 <- subset(snail1, Temperature==31)

mean_s1_13 <- mean(s1_13$logavgSpeedMM)
mean_s1_16 <- mean(s1_16$logavgSpeedMM)
mean_s1_19 <- mean(s1_19$logavgSpeedMM)
mean_s1_22 <- mean(s1_22$logavgSpeedMM)
mean_s1_25 <- mean(s1_25$logavgSpeedMM)
mean_s1_28 <- mean(s1_28$logavgSpeedMM)
mean_s1_31 <- mean(s1_31$logavgSpeedMM)
means_s1 <- c(mean_s1_13, mean_s1_16, mean_s1_19, mean_s1_22, mean_s1_25, mean_s1_28, mean_s1_31)

se_s1_13 <- sd(s1_13$logavgSpeedMM)/sqrt(length(s1_13$logavgSpeedMM))
se_s1_16 <- sd(s1_16$logavgSpeedMM)/sqrt(length(s1_16$logavgSpeedMM))
se_s1_19 <- sd(s1_19$logavgSpeedMM)/sqrt(length(s1_19$logavgSpeedMM))
se_s1_22 <- sd(s1_22$logavgSpeedMM)/sqrt(length(s1_22$logavgSpeedMM))
se_s1_25 <- sd(s1_25$logavgSpeedMM)/sqrt(length(s1_25$logavgSpeedMM))
se_s1_28 <- sd(s1_28$logavgSpeedMM)/sqrt(length(s1_28$logavgSpeedMM))
se_s1_31 <- sd(s1_31$logavgSpeedMM)/sqrt(length(s1_31$logavgSpeedMM))
ses_s1 <- c(se_s1_13, se_s1_16, se_s1_19, se_s1_22, se_s1_25, se_s1_28, se_s1_31)
lower_s1 <- means_s1-ses_s1
upper_s1 <- means_s1+ses_s1

snail2 <- subset(swimData, Snail==2)
temps2 <- c(13, 16, 19, 22, 25, 28, 34)
s2_13 <- subset(snail2, Temperature==13)
s2_16 <- subset(snail2, Temperature==16)
s2_19 <- subset(snail2, Temperature==19)
s2_22 <- subset(snail2, Temperature==22)
s2_25 <- subset(snail2, Temperature==25)
s2_28 <- subset(snail2, Temperature==28)
s2_34 <- subset(snail2, Temperature==34)

mean_s2_13 <- mean(s2_13$logavgSpeedMM)
mean_s2_16 <- mean(s2_16$logavgSpeedMM)
mean_s2_19 <- mean(s2_19$logavgSpeedMM)
mean_s2_22 <- mean(s2_22$logavgSpeedMM)
mean_s2_25 <- mean(s2_25$logavgSpeedMM)
mean_s2_28 <- mean(s2_28$logavgSpeedMM)
mean_s2_34 <- mean(s2_34$logavgSpeedMM)
means_s2 <- c(mean_s2_13, mean_s2_16, mean_s2_19, mean_s2_22, mean_s2_25, mean_s2_28, mean_s2_34)

se_s2_13 <- sd(s2_13$logavgSpeedMM)/sqrt(length(s2_13$logavgSpeedMM))
se_s2_16 <- sd(s2_16$logavgSpeedMM)/sqrt(length(s2_16$logavgSpeedMM))
se_s2_19 <- sd(s2_19$logavgSpeedMM)/sqrt(length(s2_19$logavgSpeedMM))
se_s2_22 <- sd(s2_22$logavgSpeedMM)/sqrt(length(s2_22$logavgSpeedMM))
se_s2_25 <- sd(s2_25$logavgSpeedMM)/sqrt(length(s2_25$logavgSpeedMM))
se_s2_28 <- sd(s2_28$logavgSpeedMM)/sqrt(length(s2_28$logavgSpeedMM))
se_s2_34 <- sd(s2_34$logavgSpeedMM)/sqrt(length(s2_34$logavgSpeedMM))
ses_s2 <- c(se_s2_13, se_s2_16, se_s2_19, se_s2_22, se_s2_25, se_s2_28, se_s2_34)
lower_s2 <- means_s2-ses_s2
upper_s2 <- means_s2+ses_s2

errbar(temps2, means_s2, upper_s2, lower_s2, xlim=c(12,35), ylim=c(0,2),
	xlab=expression(paste("Performance temperature (",degree,"C)")),
	ylab="Ln Cercaria swimming speed (mm/s)", pch=1)
errbar(temps1, means_s1, upper_s1, lower_s1, pch=16, add=TRUE)
curve(lnormPlotSS, from=13, to=34, lwd=2, add=TRUE)
shade <- rgb(105, 105, 105, alpha=100, maxColorValue=255)
polygon(x=c(preds$PerfTempK-K,rev(preds$PerfTempK-K)), y=c(preds$lwr_CI,rev(preds$upr_CI)), col=shade, border=NA)