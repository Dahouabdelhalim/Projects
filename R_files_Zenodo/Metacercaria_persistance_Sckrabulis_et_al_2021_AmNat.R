# 'Metacercaria_persistence_Sckrabulis_et_al_2021_AmNat.R'

# Metacercaria clearance analysis

# Use Altman et al. 2016 data here as a .csv
# Reference: Altman, K. A., S. H. Paull, P. T. J. Johnson, M. N. Golembieski, J. P. Stephens, B. E. LaFonte, and T. R. Raffel. 2016. Host and parasite thermal acclimation responses depend on the stage of infection. Journal of Animal Ecology 85(4): 1014-1024. doi:10.1111/1365-2656.12510
# Dryad repository doi: 10.5061/dryad.f3k8p
propClearData <- read.csv(file.choose(), header=TRUE)
propClearData$PerfTempK <- propClearData$PerfTemp+K
propClearData$AccTempK <- propClearData$AccTemp+K
propClearData$AccCenter <- propClearData$AccTemp-mean(propClearData$AccTemp)
propClearData <- subset(propClearData, PropCleared1!="NA")
propClearData$logPropClear <- log(propClearData$PropCleared1+1)

# Constants
k<-8.62*10^-5
K<-273.15

#_
# Initial test of normal and lognormal error distribution & model type

# Xiao et al. 2011 AIC calculation
# Reference: Xiao, X., E. P. White, M. B. Hooten, and S. L. Durham. 2011. On the use of log-transformation vs. nonlinear regression for analyzing biological power laws. Ecology 92: 1887-1894. doi:10.1890/11-0538.1
AIC.Xiao <- function(numParams, loglik, n){
	2*numParams-2*loglik
}

#_
# Boltzmann-Arrhenius

# BA equation 
ba <- function(x, Toc, cTo, Ec){
	Toc <- Toc+K
	y <- cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc))))
}

# Normal error distribution
baModel <- nls(PropCleared1 ~ ba(x=PerfTemp, Toc=19, cTo, Ec), start=c(cTo=1, Ec=.6), data=propClearData, control=c(warnOnly=TRUE))

# Plotting function using parameter estimates from 'baModel', also used for model predictions
baPlot <- function(x){
	Toc <- 19+K
	cTo <- .46322
	Ec <- .43064
	y <- cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc))))
}
plot(propClearData$PerfTemp, propClearData$PropCleared1)
curve(baPlot, from=10, to=40, add=TRUE)

# Standard deviation and log-likelihood calculations
sdnorm <- sd(propClearData$PropCleared1-baPlot(propClearData$PerfTemp))
ll_norm <- sum(log(dnorm(propClearData$PropCleared1, baPlot(propClearData$PerfTemp), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=2, loglik=ll_norm, n=length(propClearData$PerfTemp))
normAIC

# Lognormal error distribution
logbaModel <- nls(logPropClear ~ log(ba(x=PerfTemp, Toc=19, cTo, Ec)), start=c(cTo=1, Ec=.6), data=propClearData, control=c(warnOnly=TRUE))

# Plotting function using parameter estimates from 'logbaModel', also used for model predictions
logbaPlot <- function(x){
	Toc <- 19+K
	cTo <- 1.44366
	Ec <- .14767
	y <- log(cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc)))))
}
plot(propClearData$PerfTemp, propClearData$logPropClear)
curve(logbaPlot, from=10, to=40, add=TRUE)

# Standard deviation and log-likelihood calculations
sdlnorm <- sd(propClearData$logPropClear-logbaPlot(propClearData$PerfTemp))
ll_lnorm <- sum(log(dlnorm(propClearData$PropCleared1+1, logbaPlot(propClearData$PerfTemp), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=2, loglik=ll_lnorm, n=length(propClearData$PerfTemp))
lnormAIC

#_
# Sharpe-Schoolfield

# Sharpe-Schoolfield equation
ss <- function(x, Toc, cTo, Ec, Ehc, Thc){
	Toc <- Toc+K
	y <- cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc))))*(1+exp(Ehc/k*(1/Thc-1/(x+K))))^-1
}

# Normal error distribution
ssModel <- nls(PropCleared1 ~ ss(x=PerfTemp, Toc=19, cTo, Ec, Ehc=3.25, Thc), start=c(cTo=1, Ec=.6, Thc=300), data=propClearData, control=c(warnOnly=TRUE))

# Plotting function using estimates from 'ssModel', also used for model predictions
ssPlot <- function(x){
	Toc <- 19+K
	cTo <- .46815
	Ec <- .54951
	Ehc <- 3.25
	Thc <- 304.98628
	y<-cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc))))*(1+exp(Ehc/k*(1/Thc-1/(x+K))))^-1
}
plot(propClearData$PerfTemp, propClearData$PropCleared1)
curve(ssPlot, from=10, to=40, add=TRUE)

# Standard deviation and log-likelihood calculations
sdnorm <- sd(propClearData$PropCleared1-ssPlot(propClearData$PerfTemp))
ll_norm <- sum(log(dnorm(propClearData$PropCleared1, ssPlot(propClearData$PerfTemp), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=3, loglik=ll_norm, n=length(propClearData$PerfTemp))
normAIC

# Lognormal error distribution - convergence error
logssModel <- nls(logPropClear ~ log(ss(x=PerfTemp, Toc=19, cTo, Ec, Ehc=3.25, Thc)),	start=c(cTo=1, Ec=.1, Thc=300), data=propClearData, control=c(warnOnly=TRUE), trace=TRUE)

#_

# Model selection - BA

# BA model equation with acclimation effects
clearanceBA <- function(x, Tacc, Toc, bcTo, mcTo, qcTo, bEc, mEc, qEc){
	K <- 273.15
	k <- 8.62*10^-5
	Toc <- Toc+K
	cTo <- function(Tacc, bcTo, mcTo, qcTo){
		cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
	}
	Ec <- function(Tacc, bEr, mEr, qEr){
		Er <- bEr+mEr*(Tacc)+qEr*(Tacc)^2
	}
	y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc))))
}

# No acclimation effects
noAcc <- nls(PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=27.3, bcTo, mcTo=0, qcTo=0, bEc, mEc=0, qEc=0), data=propClearData, start=c(bcTo=1, bEc=.6), control=c(warnOnly=TRUE))
summary(noAcc)
AIC(noAcc)

# Linear (l) acclimation effects

# cTo
lcTo <- nls(PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=27.3, bcTo, mcTo, qcTo=0, bEc, mEc=0, qEc=0), data=propClearData,start=c(bcTo=.4,mcTo=0,bEc=.6),control=c(warnOnly=TRUE))
summary(lcTo)
AIC(lcTo)

# Ec
# Optimize To
# We comment our code here, but omit these notes in all further instances of optimizing To

# Sequence of possible To values over temperature range, we later divide by 10 such that 50/10 = 5C
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))

# For loop to iterate through every possible To value generated above
for(j in range){
	
	# Clearance model equation with acclimation effects
	clearance <- function(x, Tacc, Toc, bcTo, mcTo, qcTo, bEc, mEc, qEc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))
	}
	
	# Model for each To iterated over
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0),	start=c(bcTo=.4, bEc=.6, mEc=0), data=propClearData, control=c(warnOnly=TRUE))
	
	# Insert AIC of model into vector
	vectAIC[j] <- AIC(model)
}

# Create data frame of To and AIC values
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])

# Plot data frame to identify any local minima visually, then search for that area within the data frame for best fit model and corresponding To value
plot(df$temp, df$AIC)

#
# To= 29.9
#
# Run final model used for model selection for this parameter combination
lEc <- nls(PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=29.9, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0), control=c(warnOnly=TRUE))
summary(lEc)
AIC(lEc)

# cTo + Ec
# Optimize To
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, bcTo, mcTo, qcTo, bEc, mEc, qEc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc,bEc,mEc,qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc,bThc,mThc,qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), bcTo, mcTo, qcTo=0, bEc, mEc, qEc=0), start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To
#
lcTo_lEc <- nls(PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=27.3, bcTo, mcTo, qcTo=0, bEc, mEc, qEc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0), control=c(warnOnly=TRUE))
summary(lcTo_lEc)
AIC(lcTo_lEc)

# Quadratic (qd) acclimation effects

# cTo
qdcTo <- nls(PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=27.3, bcTo, mcTo, qcTo, bEc, mEc=0, qEc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6), control=c(warnOnly=TRUE))
summary(qdcTo)
AIC(qdcTo)

# Ec
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, bcTo, mcTo, qcTo, bEc, mEc, qEc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc), start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To= 27.3
#
qdEc <- nls(PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=27.3, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc), start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0), data=propClearData, control=c(warnOnly=TRUE))
summary(qdEc)
AIC(qdEc)

# cTo + Ec
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, bcTo, mcTo, qcTo, bEc, mEc, qEc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), bcTo, mcTo, qcTo, bEc, mEc, qEc), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, qEc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To= 25.9
#
qdcTo_qdEc <- nls(PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=(25.9), bcTo, mcTo, qcTo, bEc, mEc, qEc), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, qEc=0), data=propClearData, control=c(warnOnly=TRUE))
summary(qdcTo_qdEc)
AIC(qdcTo_qdEc)

# Linear (l) and quadratic (qd) acclimation effect combinations

# qdcTo + lEc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, bcTo, mcTo, qcTo, bEc, mEc, qEc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}

		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), bcTo, mcTo, qcTo, bEc, mEc, qEc=0), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To= between 22.8 & 22.9 = 22.85
#
qdcTo_lEc <- nls(formula=PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=(22.85), bcTo, mcTo, qcTo, bEc, mEc, qEc=0), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0),	data=propClearData, control=c(warnOnly=TRUE))
summary(qdcTo_lEc)
AIC(qdcTo_lEc)

# qdEc + lcTo
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, bcTo, mcTo, qcTo, bEc, mEc, qEc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), bcTo, mcTo, qcTo=0, bEc, mEc, qEc), start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, qEc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To= 27.1
#
qdEc_lcTo <- nls(formula=PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=(27.1), bcTo, mcTo, qcTo=0, bEc, mEc, qEc),	start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, qEc=0), data=propClearData, control=c(warnOnly=TRUE))
summary(qdEc_lcTo)
AIC(qdEc_lcTo)

# Test best fit model normal and lognormal error distribution

# Normal error distribution
qdEc <- nls(PropCleared1 ~ clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=27.3, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc), start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0), data=propClearData, control=c(warnOnly=TRUE))

# Plotting and predictions function
qdEcPlot <- function(x, Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	Toc <- 27.3
	cTo <- .743678
	Ec <- .146023+(-.015732*(Tacc))+(.013028*(Tacc)^2)
	y <- cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc+K))))
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd(propClearData$PropCleared1-qdEcPlot(x=propClearData$PerfTemp, Tacc=propClearData$AccCenter))
ll_norm <- sum(log(dnorm(propClearData$PropCleared1, qdEcPlot(x=propClearData$PerfTemp, Tacc=propClearData$AccCenter), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=4, loglik=ll_norm, n=length(propClearData$PerfTemp))
normAIC

# Lognormal error distribution
logqdEc <- nls(logPropClear ~ log(clearanceBA(x=PerfTemp, Tacc=AccCenter, Toc=27.3, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc)), start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0), data=propClearData, control=c(warnOnly=TRUE))
logqdEcPlot <- function(x, Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	Toc <- 27.3
	cTo <- 1.6968172
	Ec <- .0706991+(-.0040683*(Tacc))+(.0029153*(Tacc)^2)
	y <- log(cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc+K)))))
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd(propClearData$logPropClear-logqdEcPlot(x=propClearData$PerfTemp, Tacc=propClearData$AccCenter))
ll_lnorm <- sum(log(dlnorm(propClearData$PropCleared1+1, logqdEcPlot(x=propClearData$PerfTemp, Tacc=propClearData$AccCenter), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=4, loglik=ll_lnorm, n=length(propClearData$PerfTemp))
lnormAIC

#_
# Model selection - SS

# Acclimation model
clearanceSS <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
	K <- 273.15
	k <- 8.62*10^-5
	Toc <- Toc+K
	cTo <- function(Tacc, bcTo, mcTo, qcTo){
		cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
	}
	Ec <- function(Tacc, bEc, mEc, qEc){
		Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
	}
	Thc <- function(Tacc, bThc, mThc, qThc){
		Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
	}
	y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
}

# No acclimation effects
noAcc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc=0, qEc=0, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=1, bEc=.6, bThc=300), control=c(warnOnly=TRUE))
summary(noAcc)
AIC(noAcc)

# Linear (l) acclimation effects

# cTo
lcTo <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc=0, qEc=0, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, bThc=303), control=c(warnOnly=TRUE))
summary(lcTo)
AIC(lcTo)

# Ec
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc,bcTo,mcTo,qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0, bThc, mThc=0, qThc=0), start=c(bcTo=.4, bEc=.6, mEc=0, bThc=303), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To= between 29.3 and 29.4 = 29.35
#
lEc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=29.35, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, bThc=303), control=c(warnOnly=TRUE))
summary(lEc)
AIC(lEc)

# Thc
lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc=0, qEc=0, bThc, mThc, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(lThc)
AIC(lThc)

# cTo + Ec
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
	for(j in range){
		clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
			K <- 273.15
			k <- 8.62*10^-5
			cTo <- function(Tacc, bcTo, mcTo, qcTo){
				cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
			}
			Ec <- function(Tacc, bEc, mEc, qEc){
				Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
			}
			Thc <- function(Tacc, bThc, mThc, qThc){
				Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
			}
			y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
		}
		model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc=0, bThc, mThc=0, qThc=0), start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, bThc=303), data=propClearData, control=c(warnOnly=TRUE))
		vectAIC[j] <- AIC(model)
	}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To cannot be optimized
#
lcTo_lEc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=29.35, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc=0, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, bThc=303), control=c(warnOnly=TRUE))
summary(lcTo_lEc)
AIC(lcTo_lEc)

# cTo + Thc
lcTo_lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc=0, qEc=0, bThc, mThc, qThc=0),	data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(lcTo_lThc)
AIC(lcTo_lThc)

# Ec + Thc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0, bThc, mThc, qThc=0), start=c(bcTo=.4, bEc=.6, mEc=0, bThc=303, mThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 23
#
lEc_lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=23, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0, bThc, mThc, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(lEc_lThc)
AIC(lEc_lThc)

# cTo + Ec + Thc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc,bThc,mThc,qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc=0, bThc, mThc, qThc=0), start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, bThc=303, mThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To cannot be optimized, no local minimum
#
lcTo_lEc_lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=23, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc=0, bThc, mThc, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(lcTo_lEc_lThc)
AIC(lcTo_lEc_lThc)

# Quadratic (qd) acclimation effects
# cTo
qdcTo <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc=0, qEc=0, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, bThc=303), control=c(warnOnly=TRUE))
summary(qdcTo)
AIC(qdcTo)

# Ec
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc,bcTo,mcTo,qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc,bEc,mEc,qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc,bThc,mThc,qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 26.9
#
qdEc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303), control=c(warnOnly=TRUE))
summary(qdEc)
AIC(qdEc)

# Thc
qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc=0, qEc=0, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, bEc=.6, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(qdThc)
AIC(qdThc)

# cTo + Ec
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc=0, qThc=0), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 25.2
#
adcTo_qdEc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=25.2, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303), control=c(warnOnly=TRUE))
summary(adcTo_qdEc)
AIC(adcTo_qdEc)

# cTo + Thc
qdcTo_qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc=0, qEc=0, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(qdcTo_qdThc)
AIC(qdcTo_qdThc)

# Ec + Thc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc, qThc), start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0, qThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 24.5
#
qdEc_qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=24.5, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(qdEc_qdThc)
AIC(qdEc_qdThc)

# cTo + Ec + Thc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0, qThc=0) ,data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 25.3
#
qdcTo_qdEc_qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=25.3, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(qdcTo_qdEc_qdThc)
AIC(qdcTo_qdEc_qdThc)

# Linear (l) and quadratic (qd) acclimation effects

# qdcTo + lEc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc=0, bThc, mThc=0, qThc=0), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, bThc=303), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 22.8
#
qdcTo_lEc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=22.8, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc=0, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, bThc=303), control=c(warnOnly=TRUE))
summary(qdcTo_lEc)
AIC(qdcTo_lEc)

# qdcTo + lThc
qdcTo_lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc=0, qEc=0, bThc, mThc, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(qdcTo_lThc)
AIC(qdcTo_lThc)

# qdcTo + lEc + lThc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc=0, bThc, mThc, qThc=0), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, bThc=303, mThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 21.8
#
qdcTo_lEc_lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=21.8, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc=0, bThc, mThc, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(qdcTo_lEc_lThc)
AIC(qdcTo_lEc_lThc)

# qdcTo + qdEc + lThc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc=0), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 27.7
#
qdcTo_qdEc_lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=27.7, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(qdcTo_qdEc_lThc)
AIC(qdcTo_qdEc_lThc)

# qdcTo + lEc + qdThc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc=0, bThc, mThc, qThc), start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, bThc=303, mThc=0, qThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 15.7
#
qdcTo_lEc_qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=15.7, Ehc=3.25, bcTo, mcTo, qcTo, bEc, mEc, qEc=0, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, mcTo=0, qcTo=0, bEc=.6, mEc=0, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(qdcTo_lEc_qdThc)
AIC(qdcTo_lEc_qdThc)

# qdEc + lcTo
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 26.7
#
lcTo_qdEc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.7, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303), control=c(warnOnly=TRUE))
summary(lcTo_qdEc)
AIC(lcTo_qdEc)

# qdEc + lThc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc, qThc=0), start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model) 
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 26.6
#
qdEc_lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.6, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(qdEc_lThc)
AIC(qdEc_lThc)

# lcTo + qdEc + lThc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc, bThc, mThc, qThc=0), start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 25.1
#
lcTo_qdEc_lThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=25.1, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc, bThc, mThc, qThc=0), data=propClearData,  start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0), control=c(warnOnly=TRUE))
summary(lcTo_qdEc_lThc)
AIC(lcTo_qdEc_lThc)

# lcTo + qdEc + qdThc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc, bThc, mThc, qThc), start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0, qThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 24.6
#
lcTo_qdEc_qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=24.6, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, qEc=0, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(lcTo_qdEc_qdThc)
AIC(lcTo_qdEc_qdThc)

# lcTo + qdThc
lcTo_qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc=0, qEc=0, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(lcTo_qdThc)
AIC(lcTo_qdThc)

# lEc + qdThc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0, bThc, mThc, qThc), start=c(bcTo=.4, bEc=.6, mEc=0, bThc=303, mThc=0, qThc=0) ,data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = between 22.7 and 22.8 = 22.75
#
lEc_qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=22.75, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(lEc_qdThc)
AIC(lEc_qdThc)

# lcTo + lEc + qdThc
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	clearance <- function(x, Tacc, Toc, Ehc, bcTo, mcTo, qcTo, bEc, mEc, qEc, bThc, mThc, qThc){
		K <- 273.15
		k <- 8.62*10^-5
		cTo <- function(Tacc, bcTo, mcTo, qcTo){
			cTo <- bcTo+mcTo*(Tacc)+qcTo*(Tacc)^2
		}
		Ec <- function(Tacc, bEc, mEc, qEc){
			Ec <- bEc+mEc*(Tacc)+qEc*(Tacc)^2
		}
		Thc <- function(Tacc, bThc, mThc, qThc){
			Thc <- bThc+mThc*(Tacc)+qThc*(Tacc)^2
		}
		y <- cTo(Tacc,bcTo,mcTo,qcTo)*exp(-(Ec(Tacc,bEc,mEc,qEc)/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc(Tacc,bThc,mThc,qThc)-1/(x+K))))^-1
	}
	model <- nls(formula=PropCleared1 ~ clearance(x=PerfTemp, Tacc=AccCenter, Toc=(j/10), Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc=0, bThc, mThc, qThc), start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, bThc=303, mThc=0, qThc=0), data=propClearData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# cannot optimize To, To = 26.9
#
lcTo_lEc_qdThc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo, qcTo=0, bEc, mEc, qEc=0, bThc, mThc, qThc), data=propClearData, start=c(bcTo=.4, mcTo=0, bEc=.6, mEc=0, bThc=303, mThc=0, qThc=0), control=c(warnOnly=TRUE))
summary(lcTo_lEc_qdThc)
AIC(lcTo_lEc_qdThc)

# Best fit model test normal and lognormal error distributions

# Normal error distribution
qdEc <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303), control=c(warnOnly=TRUE))

# Plotting and predictions function
qdEcPlot <- function(x, Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	cTo <- .85786
	Ec <- .260621+(-.016666*(Tacc))+(.014539*(Tacc)^2)
	Thc <- 304.382479
	Ehc <- 3.25
	Toc <- 26.9
	y <- cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc-1/(x+K))))^-1
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd(propClearData$PropCleared1-qdEcPlot(x=propClearData$PerfTemp, Tacc=propClearData$AccCenter))
ll_norm <- sum(log(dnorm(propClearData$PropCleared1, qdEcPlot(x=propClearData$PerfTemp, Tacc=propClearData$AccCenter), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=5, loglik=ll_norm, n=length(propClearData$PerfTemp))
normAIC

# Lognormal error distribution - convergence error
logqdEc <- nls(logPropClear ~ log(clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0)), data=propClearData, start=c(bcTo=2, bEc=.6, mEc=0, qEc=0, bThc=303), control=c(warnOnly=TRUE), trace=TRUE)

# Final model ANOVA
quad <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303), control=c(warnOnly=TRUE))
lin <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc=0, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, bThc=303), control=c(warnOnly=TRUE))
reduced <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc=0, qEc=0, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, bThc=303), control=c(warnOnly=TRUE))
anova(reduced, lin)
anova(lin, quad)
summary(reduced)
summary(quad)

# Check that final model is the same as modeling the mismatch itself
qdEcMismatch <- nls(1-PropCleared1 ~ 1-clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303), control=c(warnOnly=TRUE))

#####
# Plots

#_
# Figure 2: Full model plot 6 panels by acclimation temperature

# Libraries
library(Hmisc)

# Use Altman et al. 2016 data here as a .csv
# Reference: Altman, K. A., S. H. Paull, P. T. J. Johnson, M. N. Golembieski, J. P. Stephens, B. E. LaFonte, and T. R. Raffel. 2016. Host and parasite thermal acclimation responses depend on the stage of infection. Journal of Animal Ecology 85(4): 1014-1024. doi:10.1111/1365-2656.12510
# Dryad repository doi: 10.5061/dryad.f3k8p
propClearData <- read.csv(file.choose(), header=TRUE)
propClearData$PerfTempK <- propClearData$PerfTemp+K
propClearData$AccTempK <- propClearData$AccTemp+K
propClearData$AccCenter <- propClearData$AccTemp-mean(propClearData$AccTemp)
propClearData <- subset(propClearData, PropCleared1!="NA")
propClearData$logPropClear <- log(propClearData$PropCleared1+1)
propClearData$propSurvived <- 1-propClearData$PropCleared1

# Constants
k <- 8.62*10^-5
K <- 273.15

# Final model
full <- nls(PropCleared1 ~ clearanceSS(x=PerfTemp, Tacc=AccCenter, Toc=26.9, Ehc=3.25, bcTo, mcTo=0, qcTo=0, bEc, mEc, qEc, bThc, mThc=0, qThc=0), data=propClearData, start=c(bcTo=.4, bEc=.6, mEc=0, qEc=0, bThc=303), control=c(warnOnly=TRUE))

# Calculate Ec for each acclimation temperature
AccTemps <- c(13,16,19,22,25,28)
AccTempsCenter <- AccTemps-mean(AccTemps)
Ec <- .260621+(-.016666*(AccTempsCenter))+(.014539*(AccTempsCenter)^2)

# Interaction model (infect-resist)
fullPlot <- function(x, Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	cTo <- .85786
	Ec <- .260621+(-.016666*(Tacc))+(.014539*(Tacc)^2)
	Thc <- 304.382479
	Ehc <- 3.25
	Toc <- 26.9
	inf <- 1
	res <- cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc-1/(x+K))))^-1
	y <- inf-res
	ifelse(y>1, 1, y)
}

# Infectivity curve
infectPlot <- function(x){
	y <- 1+0*x
}

# Resistance curve
resistPlot <- function(x,Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	cTo <- .85786
	Ec <- .260621+(-.016666*(Tacc))+(.014539*(Tacc)^2)
	Thc <- 304.382479
	Ehc <- 3.25
	Toc <- 26.9
	y <- cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc-1/(x+K))))^-1
	ifelse(y>1, 1, y)
}

# Vector setup for plotting means
temps <- c(13, 16, 19, 22, 25, 28)
center <- mean(propClearData$AccTemp)
accs <- c(13-center, 16-center, 19-center, 22-center, 25-center, 28-center)
means <- as.numeric(length(temps))
upps <- as.numeric(length(temps))
lows <- as.numeric(length(temps))

# Plot parameters
par(mfrow=c(3,2), omi=c(.4,.4,0,0))
	#A: Acclimation temperature = 13
	accData <- subset(propClearData, AccTemp==temps[1])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[1]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", xaxt="n")
	curve(infectPlot(x),from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[1]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	axis(side=1, labels=FALSE)
	errbar(temps, 1-means, 1-upps, 1-lows, add=TRUE, pch=1)
	text(x=34, y=1.1, expression(italic("i"[T])), cex=1.5)
	text(x=39, y=.35, expression(italic("c"[T])), cex=1.5)
	text(x=36, y=.8, expression(italic(lambda[T])), cex=1.5)
	mtext(expression(paste("A: ",italic(T[Acc])," = 13",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	#B: Acclimation temperature = 16
	accData <- subset(propClearData, AccTemp==temps[2])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))

	curve(fullPlot(x, Tacc=accs[2]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", xaxt="n", yaxt="n")
	curve(infectPlot(x),from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[2]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	axis(side=1, labels=FALSE)
	axis(side=2, labels=FALSE)
	errbar(temps, 1-means, 1-upps, 1-lows, add=TRUE, pch=1)
	mtext(expression(paste("B: ",italic(T[Acc])," = 16",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	#C: Acclimation temperature: 19
	accData <- subset(propClearData, AccTemp==temps[3])
	perf13 <-subset(accData, PerfTemp==13)
	perf16 <-subset(accData, PerfTemp==16)
	perf19 <-subset(accData, PerfTemp==19)
	perf22 <-subset(accData, PerfTemp==22)
	perf25 <-subset(accData, PerfTemp==25)
	perf28 <-subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[3]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", xaxt="n")
	curve(infectPlot(x),from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[3]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	axis(side=1, labels=FALSE)
	errbar(temps, 1-means, 1-upps, 1-lows, add=TRUE, pch=1)
	mtext(expression(paste("C: ",italic(T[Acc])," = 19",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	#D: Acclimation temperature = 22
	accData <- subset(propClearData, AccTemp==temps[4])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[4]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", xaxt="n", yaxt="n")
	curve(infectPlot(x),from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[4]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	axis(side=1, labels=FALSE)
	axis(side=2, labels=FALSE)
	errbar(temps, 1-means, 1-upps, 1-lows, add=TRUE, pch=1)
	mtext(expression(paste("D: ",italic(T[Acc])," = 22",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	#E: Acclimation temperature = 25
	accData <- subset(propClearData, AccTemp==temps[5])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[5]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="")
	curve(infectPlot(x),from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[5]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	errbar(temps, 1-means, 1-upps, 1-lows, add=TRUE, pch=1)
	mtext(expression(paste("E: ",italic(T[Acc])," = 25",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	#F: Acclimation temperature = 28
	accData <- subset(propClearData, AccTemp==temps[6])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[6]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", yaxt="n")
	curve(infectPlot(x),from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[6]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	axis(side=2, labels=FALSE)
	errbar(temps, 1-means, 1-upps, 1-lows, add=TRUE, pch=1)
	mtext(expression(paste("F: ",italic(T[Acc])," = 28",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	mtext(expression(paste("Performance temperature (",degree,"C)")), side=1, line=1.5, outer=TRUE)
	mtext("Proportion metacercariae survived", side=2, line=1.5, outer=TRUE)

#_
# Figure 3: Respiration and clearance comparison

# Libraries
library(Hmisc)

# Use Altman et al. 2016 data here as a .csv
# Reference: Altman, K. A., S. H. Paull, P. T. J. Johnson, M. N. Golembieski, J. P. Stephens, B. E. LaFonte, and T. R. Raffel. 2016. Host and parasite thermal acclimation responses depend on the stage of infection. Journal of Animal Ecology 85(4): 1014-1024. doi:10.1111/1365-2656.12510
# Dryad repository doi: 10.5061/dryad.f3k8p
propClearData <- read.csv(file.choose(), header=TRUE)
propClearData$PerfTempK <- propClearData$PerfTemp+K
propClearData$AccTempK <- propClearData$AccTemp+K
propClearData$AccCenter <- propClearData$AccTemp-mean(propClearData$AccTemp)
propClearData <- subset(propClearData, PropCleared1!="NA")
propClearData$logPropClear <- log(propClearData$PropCleared1+1)

# Use 'Uninfected_tadpole_respiration_Sckrabulis_et_al_2021_AmNat.csv' in this repository
respData <- read.csv(file.choose(), header=TRUE)
respData$PerfTempK <- respData$PerfTemp+K
respData$AccTempK <- respData$AccTemp+K
respData <- subset(respData, AccTemp!="NA"&corO2.Time.Mass!="NA")
respData$logcorO2.Time.Mass <- log(respData$corO2.Time.Mass+1)
respData$AccCenter <- respData$AccTemp-mean(respData$AccTemp)

# Use 'Activation_energy_bootstrap_Sckrabulis_et_al_2021_AmNat.csv' in this repository
eaData <- read.csv(file.choose(), header=TRUE)

# Constants
k <- 8.62*10^-5
K <- 273.15

# Plot parameters
par(mfrow=c(2,2))

#A: Respiration model

	# Plotting function using estimates from 'Uninfected_tadpole_respiration_Sckrabulis_et_al_2021_AmNat.R' final model
	fullPlot <- function(x, Tacc){
		k <- 8.62*10^-5
		K <- 273.15
		RTo <- .001287
		EhR <- 2.894
		ThR <-308.2
		ToR <- 14.55
		ER <- .607+.00001964*Tacc+.002124*Tacc^2
		y <- RTo*exp(-(ER/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR-1/(x+K))))^-1
	}

	# Vector setup for plotting means
	temps <- c(13, 16, 19, 22, 25, 28, 31, 34)
	center <- mean(respData$AccTemp)
	accs <- c(13-center, 22-center, 28-center)
	means <- as.numeric(length(temps))
	upps <- as.numeric(length(temps))
	lows <- as.numeric(length(temps))

	# Acclimation temperature = 13
	accData <- subset(respData, AccTemp==temps[1])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	perf31 <- subset(accData, PerfTemp==31)
	perf34 <- subset(accData, PerfTemp==34)

	mean13 <- mean(perf13$corO2.Time.Mass)
	mean16 <- mean(perf16$corO2.Time.Mass)
	mean19 <- mean(perf19$corO2.Time.Mass)
	mean22 <- mean(perf22$corO2.Time.Mass)
	mean25 <- mean(perf25$corO2.Time.Mass)
	mean28 <- mean(perf28$corO2.Time.Mass)
	mean31 <- mean(perf31$corO2.Time.Mass)
	mean34 <- mean(perf34$corO2.Time.Mass)

	se13 <- sd(perf13$corO2.Time.Mass)/sqrt(length(perf13$corO2.Time.Mass))
	se16 <- sd(perf16$corO2.Time.Mass)/sqrt(length(perf16$corO2.Time.Mass))
	se19 <- sd(perf19$corO2.Time.Mass)/sqrt(length(perf19$corO2.Time.Mass))
	se22 <- sd(perf22$corO2.Time.Mass)/sqrt(length(perf22$corO2.Time.Mass))
	se25 <- sd(perf25$corO2.Time.Mass)/sqrt(length(perf25$corO2.Time.Mass))
	se28 <- sd(perf28$corO2.Time.Mass)/sqrt(length(perf28$corO2.Time.Mass))
	se31 <- sd(perf31$corO2.Time.Mass)/sqrt(length(perf31$corO2.Time.Mass))
	se34 <- sd(perf34$corO2.Time.Mass)/sqrt(length(perf34$corO2.Time.Mass))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28, mean31, mean34)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28, mean31+se31, mean34+se34)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28, mean31-se31, mean34-se34)
	
	# Panel parameters
	par(mai=c(.8,.8,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[1]), from=10, to=40, lwd=2, ylim=c(0,0.007),
		xlab=expression(paste("Performance temperature (",degree,"C)")),
		ylab=expression(paste("Tadpole respiration rate (mg O"[2]*" min"^-1*" g"^-1*")")),
		col="blue", lty=1)
	errbar(jitter(temps), means, upps, lows, add=TRUE, pch=2, col="blue")
	mtext(LETTERS[1], side=3, adj=.05, line=-1.5, cex=1)

	# Acclimation temperature = 22
	accData <- subset(respData, AccTemp==temps[4])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	perf31 <- subset(accData, PerfTemp==31)
	perf34 <- subset(accData, PerfTemp==34)

	mean13 <- mean(perf13$corO2.Time.Mass)
	mean16 <- mean(perf16$corO2.Time.Mass)
	mean19 <- mean(perf19$corO2.Time.Mass)
	mean22 <- mean(perf22$corO2.Time.Mass)
	mean25 <- mean(perf25$corO2.Time.Mass)
	mean28 <- mean(perf28$corO2.Time.Mass)
	mean31 <- mean(perf31$corO2.Time.Mass)
	mean34 <- mean(perf34$corO2.Time.Mass)

	se13 <- sd(perf13$corO2.Time.Mass)/sqrt(length(perf13$corO2.Time.Mass))
	se16 <- sd(perf16$corO2.Time.Mass)/sqrt(length(perf16$corO2.Time.Mass))
	se19 <- sd(perf19$corO2.Time.Mass)/sqrt(length(perf19$corO2.Time.Mass))
	se22 <- sd(perf22$corO2.Time.Mass)/sqrt(length(perf22$corO2.Time.Mass))
	se25 <- sd(perf25$corO2.Time.Mass)/sqrt(length(perf25$corO2.Time.Mass))
	se28 <- sd(perf28$corO2.Time.Mass)/sqrt(length(perf28$corO2.Time.Mass))
	se31 <- sd(perf31$corO2.Time.Mass)/sqrt(length(perf31$corO2.Time.Mass))
	se34 <- sd(perf34$corO2.Time.Mass)/sqrt(length(perf34$corO2.Time.Mass))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28, mean31, mean34)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28, mean31+se31, mean34+se34)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28, mean31-se31, mean34-se34)

	curve(fullPlot(x, Tacc=accs[2]), from=10, to=40, lwd=2, col="purple", add=TRUE, lty=2)
	errbar(jitter(temps), means, upps, lows, add=TRUE, pch=0, col="purple")

	# Acclimation temperature = 28
	accData <- subset(respData, AccTemp==temps[6])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	perf31 <- subset(accData, PerfTemp==31)
	perf34 <- subset(accData, PerfTemp==34)

	mean13 <- mean(perf13$corO2.Time.Mass)
	mean16 <- mean(perf16$corO2.Time.Mass)
	mean19 <- mean(perf19$corO2.Time.Mass)
	mean22 <- mean(perf22$corO2.Time.Mass)
	mean25 <- mean(perf25$corO2.Time.Mass)
	mean28 <- mean(perf28$corO2.Time.Mass)
	mean31 <- mean(perf31$corO2.Time.Mass)
	mean34 <- mean(perf34$corO2.Time.Mass)

	se13 <- sd(perf13$corO2.Time.Mass)/sqrt(length(perf13$corO2.Time.Mass))
	se16 <- sd(perf16$corO2.Time.Mass)/sqrt(length(perf16$corO2.Time.Mass))
	se19 <- sd(perf19$corO2.Time.Mass)/sqrt(length(perf19$corO2.Time.Mass))
	se22 <- sd(perf22$corO2.Time.Mass)/sqrt(length(perf22$corO2.Time.Mass))
	se25 <- sd(perf25$corO2.Time.Mass)/sqrt(length(perf25$corO2.Time.Mass))
	se28 <- sd(perf28$corO2.Time.Mass)/sqrt(length(perf28$corO2.Time.Mass))
	se31 <- sd(perf31$corO2.Time.Mass)/sqrt(length(perf31$corO2.Time.Mass))
	se34 <- sd(perf34$corO2.Time.Mass)/sqrt(length(perf34$corO2.Time.Mass))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28, mean31, mean34)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28, mean31+se31, mean34+se34)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28, mean31-se31, mean34-se34)

	curve(fullPlot(x, Tacc=accs[3]), from=10, to=40, lwd=2, col="red", add=TRUE, lty=3)
	errbar(jitter(temps), means, upps, lows, add=TRUE, pch=1, col="red")

#B: Clearance model summary with grouped acclimation temperatures for datapoints

	# Plotting function using estimates from best fit model above
	clearPlot <- function(x, Tacc){
		k <- 8.62*10^-5
		K <- 273.15
		cTo <- .85786
		Ec <- .260621+(-.016666*(Tacc))+(.014539*(Tacc)^2)
		Thc <- 304.382479
		Ehc <- 3.25
		Toc <- 26.9
		y <- cTo*exp(-(Ec/k*(1/(x+K)-1/(Toc+K))))*(1+exp(Ehc/k*(1/Thc-1/(x+K))))^-1
		ifelse(y>1, 1, y)
	}
	
	# Vector setup for plotting means
	temps <- c(13, 16, 19, 22, 25, 28)
	center <- mean(propClearData$AccTemp)
	accs <- c(13-center, 16-center, 19-center, 22-center, 25-center, 28-center)
	means <- as.numeric(length(temps))
	upps <- as.numeric(length(temps))
	lows <- as.numeric(length(temps))
	
	#Acclimation temperature = 13+16)
	accData <- subset(propClearData, AccTemp==13|AccTemp==16)
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)
	
	# Jitter is saved to replicate points in Panel D
	jitterLow <- jitter(temps)

	# Panel parameters
	par(mai=c(.8,.8,.2,.2))
	
	# Only plot curve from acclimation temperature 13
	curve(clearPlot(x, Tacc=accs[1]), from=10, to=40, ylim=c(0,1), lwd=2,
		xlab=expression(paste("Performance temperature (",degree,"C)")),
		ylab="Proportion of metacercariae cleared", col="blue")
	errbar(jitterLow, means, upps, lows, add=TRUE, pch=2, col="blue")
	mtext(LETTERS[2], side=3, adj=.05, line=-1.5, cex=1)

	#Acclimation temperature 19+22
	accData <- subset(propClearData, AccTemp==19|AccTemp==22)
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)
	
	# Jitter is saved to replicate points in Panel D
	jitterMid <- jitter(temps)

	# Only plot curve from acclimation temperature = 19
	curve(clearPlot(x, Tacc=accs[3]), from=10, to=40, ylim=c(0,1), lty=2, lwd=2, col="purple", add=TRUE)
	errbar(jitterMid, means, upps, lows, add=TRUE, pch=0, col="purple")

	# Acclimation temperature = 25
	accData <- subset(propClearData, AccTemp==25|AccTemp==28)
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)
	
	# Jitter is saved to replicate points in Panel D
	jitterHigh <- jitter(temps)

	# Only plot curve from acclimation temperature = 28
	curve(clearPlot(x, Tacc=accs[6]), from=10, to=40, ylim=c(0,1), lty=3, lwd=2, col="red", add=TRUE)
	errbar(jitterHigh, means, upps, lows, add=TRUE, pch=1, col="red")

#C: Activation energy comparisons

	# Panel parameters
	par(mai=c(.8,.8,.2,.2))
	
	plot(eaData$AccTemp, eaData$resp, type="l",
		xlab=expression(paste("Acclimation temperature (",degree,"C)")),xlim=c(13,28),
		ylab="Activation energy (eV)", ylim=c(0,1.8), lty=2, lwd=2)
	lines(eaData$AccTemp, eaData$rlow,lty=2)
	lines(eaData$AccTemp, eaData$rhigh,lty=2)
		shade<-rgb(105, 105, 105, alpha=100, maxColorValue=255)
	polygon(x=c(eaData$AccTemp, rev(eaData$AccTemp)), y=c(eaData$rlow, rev(eaData$rhigh)), col=shade, border=NA)

	lines(eaData$AccTemp, eaData$clear, lwd=2)
	lines(eaData$AccTemp, eaData$clow)
	lines(eaData$AccTemp, eaData$chigh)
	polygon(x=c(eaData$AccTemp, rev(eaData$AccTemp)), y=c(eaData$clow, rev(eaData$chigh)), col=shade, border=NA)
	mtext(LETTERS[3], side=3, adj=.05, line=-1.5, cex=1)

#D: Clearance data described by respiration parameters

	# Plotting function with estimates for EhR, ThR, ER from final model in 'Uninfected tadpole respiration - Sckrabulis et al 2021 AmNat.R'
	# Estimates for A and Tor are from the best fit clearance model above, in order to only show the predictive effect of activation energy
	clearWrespPlot <- function(x, Tacc){
		k <- 8.62*10^-5
		K <- 273.15
		RTo <- .82675
		EhR <- 2.894
		ThR <- 308.2
		ToR <- 26.9
		ER <- .607+.00001964*Tacc+.002124*Tacc^2
		y <- RTo*exp(-(ER/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR-1/(x+K))))^-1
		ifelse(y>1, 1, y)
	}
	# Acclimation temperature = 13
	accData <- subset(propClearData, AccTemp==13|AccTemp==16)
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.8,.8,.2,.2))
	
	# Plot only curve from acclimation temperature = 13
	curve(clearWrespPlot(x, Tacc=accs[1]), from=10, to=40, ylim=c(0,1), lwd=2,
		xlab=expression(paste("Performance temperature (",degree,"C)")),
		ylab="Proportion of metacercariae cleared", col="blue")
	errbar(jitterLow, means, upps, lows, add=TRUE, pch=2, col="blue")
	mtext(LETTERS[4], side=3, adj=.05, line=-1.5, cex=1)

	# Acclimation temperature = 22
	accData <-subset(propClearData, AccTemp==19|AccTemp==22)
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Plot on curve from acclimation temperature = 22
	curve(clearWrespPlot(x, Tacc=accs[3]), from=10, to=40, ylim=c(0,1), lty=2, lwd=2, col="purple", add=TRUE)
	errbar(jitterMid, means, upps, lows, add=TRUE, pch=0, col="purple")

	# Acclimation temperature = 28
	accData <- subset(propClearData, AccTemp==25|AccTemp==28)
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)

	mean13 <- mean(perf13$PropCleared1)
	mean16 <- mean(perf16$PropCleared1)
	mean19 <- mean(perf19$PropCleared1)
	mean22 <- mean(perf22$PropCleared1)
	mean25 <- mean(perf25$PropCleared1)
	mean28 <- mean(perf28$PropCleared1)

	se13 <- sd(perf13$PropCleared1)/sqrt(length(perf13$PropCleared1))
	se16 <- sd(perf16$PropCleared1)/sqrt(length(perf16$PropCleared1))
	se19 <- sd(perf19$PropCleared1)/sqrt(length(perf19$PropCleared1))
	se22 <- sd(perf22$PropCleared1)/sqrt(length(perf22$PropCleared1))
	se25 <- sd(perf25$PropCleared1)/sqrt(length(perf25$PropCleared1))
	se28 <- sd(perf28$PropCleared1)/sqrt(length(perf28$PropCleared1))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	curve(clearWrespPlot(x, Tacc=accs[6]), from=10, to=40, ylim=c(0,1), lty=3, lwd=2, col="red", add=TRUE)
	errbar(jitterHigh, means, upps, lows, add=TRUE, pch=1, col="red")