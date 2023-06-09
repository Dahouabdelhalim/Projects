# 'Uninfected_tadpole_respiration_Sckrabulis_et_al_2021_AmNat.R'

# Uninfected tadpole respitation analysis

# Libraries
library(Hmisc)

# Use 'Uninfected_tadpole_respiration_Sckrabulis_et_al_2021_AmNat.csv' in this repository
respData <- read.csv(file.choose(), header=TRUE)
respData$PerfTempK <- respData$PerfTemp+K
respData$AccTempK <- respData$AccTemp+K
respData <- subset(respData, AccTemp!="NA"&corO2.Time.Mass!="NA")
respData$logcorO2.Time.Mass <- log(respData$corO2.Time.Mass+1)
respData$AccCenter <- respData$AccTemp-mean(respData$AccTemp)

# Constants
k <- 8.62*10^-5
K <- 273.15

#_ Initial test of normal and lognormal error distribution and model type

# Xiao et al. 2011 AIC calculation
# Reference: Xiao, X., E. P. White, M. B. Hooten, and S. L. Durham. 2011. On the use of log-transformation vs. nonlinear regression for analyzing biological power laws. Ecology 92: 1887-1894. doi:10.1890/11-0538.1
AIC.Xiao <- function(numParams, loglik, n){
	2*numParams-2*loglik
}

#_ Boltzmann-Arrhenius

# BA equation
ba <- function(x, ToR, RTo, ER){
	ToR <- ToR+K
	y <- RTo*exp(-(ER/k*(1/(x+K)-1/(ToR))))
}

# Normal error distribution
baNorm <- nls(corO2.Time.Mass ~ ba(x=PerfTemp, ToR=19, RTo, ER), data=respData, start=c(RTo=1, ER=.6), control=c(warnOnly=TRUE))

baNormPlot <- function(x){
	k <- 8.62*10^-5
	K <- 273.15
	RTo <- .0021595
	ER <- .4631972
	ToR <- 19+K
	y <- RTo*exp(-(ER/k*(1/(x+K)-1/(ToR))))
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd(respData$corO2.Time.Mass-baNormPlot(x=respData$PerfTemp))
ll_norm <- sum(log(dnorm(respData$corO2.Time.Mass, baNormPlot(x=respData$PerfTemp), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=2, loglik=ll_norm, n=length(respData$PerfTemp))
normAIC

# Lognormal error distribution
baLognorm <- nls(logcorO2.Time.Mass ~ log(ba(x=PerfTemp, ToR=19, RTo, ER)), data=respData, start=c(RTo=1, ER=.6), control=c(warnOnly=TRUE))

baLognormPlot <- function(x){
	k <- 8.62*10^-5
	K <- 273.15
	RTo <- 1.002
	ER <- .001478
	ToR <- 19+K
	y <- log(RTo*exp(-(ER/k*(1/(x+K)-1/(ToR)))))
}

# Standard deviation and log-likelihood calculation
sdlnorm <- sd(respData$logcorO2.Time.Mass-baLognormPlot(x=respData$PerfTemp))
ll_lnorm <- sum(log(dlnorm(respData$corO2.Time.Mass+1, baLognormPlot(x=respData$PerfTemp), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=2, loglik=ll_lnorm, n=length(respData$PerfTemp))
lnormAIC

#_ Sharpe-Schoolfield

# SS equation
ss <- function(x, ToR, RTo, ER, EhR, ThR){
	ToR <- ToR+K
	y <- RTo*exp(-(ER/k*(1/(x+K)-1/(ToR))))*(1+exp(EhR/k*(1/ThR-1/(x+K))))^-1
}

# Normal error distribution
ssNorm <- nls(corO2.Time.Mass ~ ss(x=PerfTemp, ToR=19, RTo, ER, EhR, ThR), start=c(RTo=.1, ER=.65, EhR=3, ThR=300), data=respData, control=c(warnOnly=TRUE))

ssNormPlot <- function(x){
	k <- 8.62*10^-5
	K <- 273.15
	RTo <- .001967
	ThR <- 308.2
	ER <- .6878
	ToR <- 19+K
	EhR <- 3.111
	y <- RTo*exp(-(ER/k*(1/(x+K)-1/ToR)))*(1+exp(EhR/k*(1/ThR-1/(x+K))))^-1
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd((respData$corO2.Time.Mass)-(ssNormPlot(respData$PerfTemp)))
ll_norm <- sum(log(dnorm(respData$corO2.Time.Mass, ssNormPlot(respData$PerfTemp), sdnorm)))

# AICc calculation
normAIC <- AIC.Xiao(numParams=4, loglik=ll_norm, n=length(respData$PerfTemp))
normAIC

# Lognormal error distribution - convergence warning (step factor)
ssLognorm <- nls(logcorO2.Time.Mass ~ log(ss(x=PerfTemp, ToR=19, RTo, ER, EhR, ThR)),	start=c(RTo=.1, ER=.65, EhR=4, ThR=40+K), data=respData, control=c(warnOnly=TRUE))

ssLognormPlot<-function(x){
	RTo <- 1.002
	ThR <- 309.3
	ER <- .001546
	ToR <- 19+K
	EhR <- 20.07
	y <- log(RTo*exp(-(ER/k*(1/(x+K)-1/ToR)))*(1+exp(EhR/k*(1/ThR-1/(x+K))))^-1)
}

# Standard deviation and log-likelihood calculation
sdlnorm <- sd(respData$logcorO2.Time.Mass-(ssLognormPlot(respData$PerfTemp)))
ll_lnorm <- sum(log(dlnorm(respData$corO2.Time.Mass+1, ssLognormPlot(respData$PerfTemp), sdlnorm)))

# AICc calculation
lnormAIC <- AIC.Xiao(numParams=4, loglik=ll_lnorm, n=length(respData$PerfTemp))
lnormAIC

#_ 

# Model selection - BA

#BA model equation with acclimation effects
respBA <- function(x, Tacc, ToR, bRTo, mRTo, qRTo, bER, mER, qER){
	K <- 273.15
	k <- 8.62*10^-5
	RTo <- function(Tacc, bRTo, mRTo, qRTo){
		RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
	}
	ER <- function(Tacc, bER, mER, qER){
		Ec <- bER+mER*(Tacc)+qER*(Tacc)^2
	}
	y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))
}

# No acclimation effects
noAcc <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo=0, qRTo=0, bER, mER=0, qER=0), data=respData, start=c(bRTo=.1, bER=.6))
summary(noAcc)
AIC(noAcc)

# Linear (l) acclimation effects
# RTo
lRTo <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo, qRTo=0, bER, mER=0, qER=0), data=respData, start=c(bRTo=.1, mRTo=0, bER=.6))
summary(lRTo)
AIC(lRTo)

# ER
# Optimize To
# We comment our code here but omit these notes in all further instances of optimizing To

# Sequence of possible To values over temperature range, we later divide by 10 such that 50/10 = 5C
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))

# For loop to iterate through every possible To value generated above
for(j in range){

	# Respiration model equation with acclimation effects
	resp <- function(x, Tacc, ToR, bRTo, mRTo, qRTo, bER, mER, qER){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))
	}
	
	# Model for each To iterated over
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), bRTo, mRTo=0, qRTo=0, bER, mER, qER=0), start=c(bRTo=.4, bER=.6, mER=0), data=respData, control=c(warnOnly=TRUE))
	
	# Insert AIC of model into vector
	vectAIC[j] <- AIC(model)
}

# Create data frame of To and AIC values
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])

# Plot data frame to identify any local minima, then search for that area within the data frame for best fit model and corresponding To value
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 12.95
#

# Run final model used for model selection for this parameter combination
lER <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo=0, qRTo=0, bER, mER, qER=0), data=respData, start=c(bRTo=.1, bER=.6, mER=0))
summary(lER)
AIC(lER)

# RTo + ER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, bRTo, mRTo, qRTo, bER, mER, qER){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), bRTo, mRTo, qRTo=0, bER, mER, qER=0), start=c(bRTo=.4, mRTo=0, bER=.6, mER=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp,df$AIC)
#
# Cannot optimize To, To = 12.95
#
lRTo_lER <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo, qRTo=0, bER, mER, qER=0), data=respData, start=c(bRTo=.1,mRTo=0,bER=.6,mER=0))
summary(lRTo_lER)
AIC(lRTo_lER)

# Quadratic (qd) acclimation effects
# RTo
qdRTo <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo, qRTo, bER, mER=0, qER=0), data=respData, start=c(bRTo=.1, mRTo=0, qRTo=0, bER=.6))
summary(qdRTo)
AIC(qdRTo)

#ER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, bRTo, mRTo, qRTo, bER, mER, qER){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), bRTo, mRTo=0, qRTo=0, bER, mER, qER), start=c(bRTo=.4, bER=.6, mER=0, qER=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
#To between 12.9 and 13.0, To= 12.95
#
qdER <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo=0, qRTo=0, bER, mER, qER),	data=respData, start=c(bRTo=.1, bER=.6, mER=0, qER=0))
summary(qdER)
AIC(qdER)

# RTo + ER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, bRTo, mRTo, qRTo, bER, mER, qER){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), bRTo, mRTo, qRTo, bER, mER, qER), start=c(bRTo=.4, mRTo=0, qRTo=0, bER=.6, mER=0, qER=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 12.95
#
qdRTo_qdER <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo, qRTo, bER, mER, qER), data=respData, start=c(bRTo=.1, mRTo=0, qRTo=0, bER=.6, mER=0, qER=0), control=c(warnOnly=TRUE))
summary(qdRTo_qdER)
AIC(qdRTo_qdER)

# Quadratic (qd) and linear (l) acclimation effects
# qdRTo + lER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, bRTo, mRTo, qRTo, bER, mER, qER){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), bRTo, mRTo, qRTo, bER, mER, qER=0), start=c(bRTo=.4, mRTo=0, qRTo=0, bER=.6, mER=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 12.95
#
qdRTo_lER <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo, qRTo, bER, mER, qER=0), data=respData, start=c(bRTo=.1, mRTo=0, qRTo=0, bER=.6, mER=0))
summary(qdRTo_lER)
AIC(qdRTo_lER)

# lRTo + qdER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, bRTo, mRTo, qRTo, bER, mER, qER){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), bRTo, mRTo, qRTo=0, bER, mER, qER), start=c(bRTo=.4, mRTo=0, bER=.6, mER=0, qER=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 17
#
lRTo_qdER <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=17, bRTo, mRTo, qRTo=0, bER, mER, qER), data=respData, start=c(bRTo=.4, mRTo=0, bER=.6, mER=0, qER=0), control=c(warnOnly=TRUE))
summary(lRTo_qdER)
AIC(lRTo_qdER)

# Test best fit model normal and lognormal error distribution

# Normal error distribution
qdERNormBA <- nls(corO2.Time.Mass ~ respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo=0, qRTo=0, bER, mER, qER), data=respData, start=c(bRTo=.1, bER=.6, mER=0, qER=0))

# Plotting and predictions function
qdERNormBAPlot <- function(x, Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	RTo <- .0014658
	ER <- .3896329+(-.0002459*Tacc)+(.0018773*Tacc^2)
	ToR <- 12.95+K
	y <- RTo*exp(-ER/k*(1/(x+K)-1/(ToR)))
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd(respData$corO2.Time.Mass-qdERNormBAPlot(x=respData$PerfTemp, Tacc=respData$AccCenter))
ll_norm <- sum(log(dnorm(respData$corO2.Time.Mass, qdERNormBAPlot(x=respData$PerfTemp, Tacc=respData$AccCenter), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=4, loglik=ll_norm, n=length(respData$PerfTemp))
normAIC

# Lognormal error distribution
qdERLognormBA <- nls(logcorO2.Time.Mass ~ log(respBA(x=PerfTemp, Tacc=AccCenter, ToR=12.95, bRTo, mRTo=0, qRTo=0, bER, mER, qER)), data=respData, start=c(bRTo=.1, bER=.6, mER=0, qER=0))
	
# Plotting and predictions function
qdERLognormBAPlot <- function(x, Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	RTo <- 1.001
	ER <- .001215+(-.000004093*Tacc)+(.000007139*Tacc^2)
	ToR <- 12.95
	y <- log(RTo*exp(-(ER/k*(1/(x+K)-1/(ToR+K)))))
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd(respData$logcorO2.Time.Mass-qdERLognormBAPlot(x=respData$PerfTemp, Tacc=respData$AccCenter))
ll_lnorm <- sum(log(dlnorm(respData$corO2.Time.Mass+1, qdERLognormBAPlot(x=respData$PerfTemp, Tacc=respData$AccCenter), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=4, loglik=ll_lnorm, n=length(respData$PerfTemp))
lnormAIC

#_ 
# Model selection - SS

# Acclimation model
respSS <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
	K <- 273.15
	k <- 8.62*10^-5
	RTo <- function(Tacc, bRTo, mRTo, qRTo){
		RTo <- bRTo+mRTo*(Tacc+K)+qRTo*(Tacc)^2
	}
	ER <- function(Tacc, bER, mER, qER){
		Ec <- bER+mER*(Tacc+K)+qER*(Tacc)^2
	}
	ThR <- function(Tacc, bThR, mThR, qThR){
		ThR <- bThR+mThR*(Tacc+K)+qThR*(Tacc)^2
	}
	y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
}

# No acclimation effects
noAcc <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER=0, qER=0, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, bThR=300))
summary(noAcc)
AIC(noAcc)

# Linear (l) acclimation effects
# RTo
lRTo <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo=0, bER, mER=0, qER=0, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, mRTo=0, bER=.6, bThR=300))
summary(lRTo)
AIC(lRTo)

#ER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER=0, bThR, mThR=0, qThR=0), start=c(EhR=3, bRTo=.4, bER=.6, mER=0, bThR=305), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# cannot optimize To, To = 14.55
#
lER <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER=0, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, mER=0, bThR=300))
summary(lER)
AIC(lER)

# ThR
lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER=0, qER=0, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, bThR=300, mThR=0))
summary(lThR)
AIC(lThR)

# RTo + ER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo=0, bER, mER, qER=0, bThR, mThR=0, qThR=0), start=c(EhR=3, bRTo=.4, mRTo=0, bER=.6, mER=0, bThR=305), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
lRTo_lER <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo=0, bER, mER, qER=0, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, mRTo=0, bER=.6, mER=0, bThR=300))
summary(lRTo_lER)
AIC(lRTo_lER)

# RTo + ThR
lRTo_lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo=0, bER, mER=0, qER=0, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, mRTo=0, bER=.6, bThR=310, mThR=0))
summary(lRTo_lThR)
AIC(lRTo_lThR)

# ER + ThR
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER=0, bThR, mThR, qThR=0), start=c(EhR=3, bRTo=.4, bER=.6, mER=0, bThR=305, mThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
lER_lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER=0, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.001, bER=.6, mER=0, bThR=310, mThR=0))
summary(lER_lThR)
AIC(lER_lThR)

# RTo + ER + ThR
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo=0, bER, mER, qER=0, bThR, mThR, qThR=0), start=c(EhR=3, bRTo=.4, mRTo=0, bER=.6, mER=0, bThR=305, mThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
lRTo_lER_lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo=0, bER, mER, qER=0, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.001, mRTo=0, bER=.6, mER=0, bThR=310, mThR=0), control=c(warnOnly=TRUE))
summary(lRTo_lER_lThR)
AIC(lRTo_lER_lThR)

# Quadratic (qd) acclimation effects
# RTo
qdRTo <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER=0, qER=0, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, mRTo=0, qRTo=0, bER=.6, bThR=300))
summary(qdRTo)
AIC(qdRTo)

# ER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR=0, qThR=0), start=c(EhR=3, bRTo=.4, bER=.6, mER=0, qER=0, bThR=305), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To between 14.5 and 14.6, To = 14.55
#
qdER <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, mER=0, qER=0, bThR=300))
summary(qdER)
AIC(qdER)

# ThR
qdThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER=0, qER=0, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, bThR=300, mThR=0, qThR=0))
summary(qdThR)
AIC(qdThR)

# RTo + ER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo ,qRTo, bER, mER, qER, bThR, mThR=0, qThR=0), start=c(EhR=3, bRTo=.4, mRTo=0, qRTo=0, bER=.6, mER=0, qER=0, bThR=305), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
qdRTo_qdER <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, mRTo=0, qRTo=0, bER=.6, mER=0, qER=0, bThR=300))
summary(qdRTo_qdER)
AIC(qdRTo_qdER)

# RTo + ThR
qdRTo_qdThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER=0, qER=0, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.1, mRTo=0, qRTo=0, bER=.6, bThR=300, mThR=0, qThR=0), control=c(warnOnly=TRUE))
summary(qdRTo_qdThR)
AIC(qdRTo_qdThR)

# ER + ThR
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR, qThR), start=c(EhR=3, bRTo=.4, bER=.6, mER=0, qER=0, bThR=305, mThR=0, qThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 11.6
#
qdER_qdThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=11.6, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.4, bER=.6, mER=0, qER=0, bThR=305, mThR=0, qThR=0))
summary(qdER_qdThR)
AIC(qdER_qdThR)

# RTo + ER + ThR
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR), start=c(EhR=3, bRTo=.4, mRTo=0, qRTo=0, bER=.6, mER=0, qER=0, bThR=305, mThR=0, qThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
qdRTo_qdER_qdThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, qRTo=0, bER=.6, mER=0, qER=0, bThR=308, mThR=0, qThR=0))
summary(qdRTo_qdER_qdThR)
AIC(qdRTo_qdER_qdThR)

# Linear (l) and quadratic (qd) acclimation effect combinations
# qdRTo + lER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo, bER, mER, qER=0, bThR, mThR=0, qThR=0), start=c(EhR=3, bRTo=.4, mRTo=0, qRTo=0, bER=.6, mER=0, bThR=305), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
qdRTo_lER <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER, qER=0, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, qRTo=0, bER=.6, mER=0, bThR=308))
summary(qdRTo_lER)
AIC(qdRTo_lER)

# qdRTo + lThR
qdRTo_lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER=0, qER=0, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, qRTo=0, bER=.6, bThR=308, mThR=0))
summary(qdRTo_lThR)
AIC(qdRTo_lThR)

# qdRTo + lER + lThR
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo, bER, mER, qER=0, bThR, mThR, qThR=0), start=c(EhR=3, bRTo=.4, mRTo=0, qRTo=0, bER=.6, mER=0, bThR=305, mThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
qdRTo_lER_lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER, qER=0, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, qRTo=0, bER=.6, mER=0, bThR=308, mThR=0))
summary(qdRTo_lER_lThR)
AIC(qdRTo_lER_lThR)

# qdRTo + qdER + mThR
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR=0), start=c(EhR=3, bRTo=.4, mRTo=0, qRTo=0, bER=.6, mER=0, qER=0, bThR=305, mThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10,AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
qdRTo_qdER_lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, qRTo=0, bER=.6, mER=0, qER=0, bThR=308, mThR=0))
summary(qdRTo_qdER_lThR)
AIC(qdRTo_qdER_lThR)

# qdRTo + qdThR + lER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo, bER, mER, qER=0, bThR, mThR, qThR), start=c(EhR=3, bRTo=.4, mRTo=0, qRTo=0, bER=.6, mER=0, bThR=305, mThR=0, qThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# Cannot optimize To, To = 14.55
#
qdRTo_qdThR_lER <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo, bER, mER, qER=0, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, qRTo=0, bER=.6, mER=0, bThR=308, mThR=0,qThR=0))
summary(qdRTo_qdThR_lER)
AIC(qdRTo_qdThR_lER)

# qdER + lRTo
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, Edr, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo=0, bER, mER, qER, bThR, mThR=0, qThR=0), start=c(EhR=3, bRTo=.4, mRTo=0, bER=.6, mER=0, qER=0, bThR=305), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To between 18.3 and 18.4, To = 18.35
#
qdER_lRTo <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=18.35, EhR, bRTo, mRTo, qRTo=0, bER, mER, qER, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, bER=.6, mER=0, qER=0, bThR=308))
summary(qdER_lRTo)
AIC(qdER_lRTo)

#qdER + lThR
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR, qThR=0), start=c(EhR=3, bRTo=.4, bER=.6, mER=0, qER=0, bThR=305, mThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To between 14.3 and 14.4 = 14.35
#
qdER_lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.35, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.01, bER=.6, mER=0, qER=0, bThR=308, mThR=0))
summary(qdER_lThR)
AIC(qdER_lThR)

# qdER + lRTo + lThR
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo=0, bER, mER, qER, bThR, mThR, qThR=0), start=c(EhR=3, bRTo=.4, mRTo=0, bER=.6, mER=0, qER=0, bThR=305, mThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 18.4
#
qdER_lRTo_lThR <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=18.4, EhR, bRTo, mRTo, qRTo=0, bER, mER, qER, bThR, mThR, qThR=0), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, bER=.6, mER=0, qER=0, bThR=308, mThR=0))
summary(qdER_lRTo_lThR)
AIC(qdER_lRTo_lThR)

#qdER + qThR + lRTo
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo=0, bER, mER, qER, bThR, mThR, qThR), start=c(EhR=3, bRTo=.4, mRTo=0, bER=.6, mER=0, qER=0, bThR=305, mThR=0, qThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To between 19.1 and 19.2, To = 19.15
#
qdER_qdThR_lRTo <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=19.15, EhR, bRTo, mRTo, qRTo=0, bER, mER, qER, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, bER=.6, mER=0, qER=0, bThR=308, mThR=0, qThR=0))
summary(qdER_qdThR_lRTo)
AIC(qdER_qdThR_lRTo)

#qdThR + lRTo
qdThR_lRTo <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo, qRTo=0, bER, mER=0, qER=0, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, bER=.6, bThR=308, mThR=0, qThR=0), control=c(warnOnly=TRUE))
summary(qdThR_lRTo)
AIC(qdThR_lRTo)

#qdThR + lER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER=0, bThR, mThR, qThR), start=c(EhR=3, bRTo=.4, bER=.6, mER=0, bThR=305, mThR=0, qThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 32
#
qdThR_lER <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=32, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER=0, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.01, bER=.6, mER=0, bThR=308, mThR=0, qThR=0), control=c(warnOnly=TRUE))
summary(qdThR_lER)
AIC(qdThR_lER)

# qdThR + lRTo + lER
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	resp <- function(x, Tacc, ToR, EhR, bRTo, mRTo, qRTo, bER, mER, qER, bThR, mThR, qThR){
		K <- 273.15
		k <- 8.62*10^-5
		RTo <- function(Tacc, bRTo, mRTo, qRTo){
			RTo <- bRTo+mRTo*(Tacc)+qRTo*(Tacc)^2
		}
		ER <- function(Tacc, bER, mER, qER){
			ER <- bER+mER*(Tacc)+qER*(Tacc)^2
		}
		ThR <- function(Tacc, bThR, mThR, qThR){
			ThR <- bThR+mThR*(Tacc)+qThR*(Tacc)^2
		}
		y <- RTo(Tacc,bRTo,mRTo,qRTo)*exp(-(ER(Tacc,bER,mER,qER)/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR(Tacc,bThR,mThR,qThR)-1/(x+K))))^-1
	}
	model <- nls(formula=corO2.Time.Mass ~ resp(x=PerfTemp, Tacc=AccCenter, ToR=(j/10), EhR, bRTo, mRTo, qRTo=0, bER, mER, qER=0, bThR, mThR, qThR), start=c(EhR=3, bRTo=.4, mRTo=0, bER=.6, mER=0, bThR=305, mThR=0, qThR=0), data=respData, control=c(warnOnly=TRUE))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To between 32.0 and 32.1, To = 32.05
#
qdThR_lRTo_lER <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=32.05, EhR, bRTo, mRTo, qRTo=0, bER, mER, qER=0, bThR, mThR, qThR), data=respData, start=c(EhR=3, bRTo=.01, mRTo=0, bER=.6, mER=0, bThR=308, mThR=0, qThR=0), control=c(warnOnly=TRUE))
summary(qdThR_lRTo_lER)
AIC(qdThR_lRTo_lER)

# Test best fit model normal and lognormal error distribution

# Normal error distribution
qdERNormSS <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, mER=0, qER=0, bThR=300))

# Plotting and predictions function
qdERNormSSPlot <- function(x, Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	RTo <- .001287
	EhR <- 2.894
	ThR <- 308.2
	ToR <- 14.55
	ER <- .607+.00001964*Tacc+.002124*Tacc^2
	y <- RTo*exp(-(ER/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR-1/(x+K))))^-1
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd(respData$corO2.Time.Mass-qdERNormSSPlot(x=respData$PerfTemp, Tacc=respData$AccCenter))
ll_norm <- sum(log(dnorm(respData$corO2.Time.Mass, qdERNormSSPlot(x=respData$PerfTemp, Tacc=respData$AccCenter), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=6, loglik=ll_norm, n=length(respData$PerfTemp))
normAIC

# Lognormal error distribution
qdERLognormSS <- nls(logcorO2.Time.Mass ~ log(respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR=0, qThR=0)), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, mER=0, qER=0, bThR=300), control=c(warnOnly=TRUE))
qdERLognormSSPlot <- function(x, Tacc){
	k <- 8.62*10^-5
	K <- 273.15
	RTo <- 1.001
	EhR <- 17.96
	ThR <- 309.9
	ToR <- 14.55
	ER <- .002112+(-.000003362)*Tacc+.000008261*Tacc^2
	y <- log(RTo*exp(-(ER/k*(1/(x+K)-1/(ToR+K))))*(1+exp(EhR/k*(1/ThR-1/(x+K))))^-1)
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd(respData$logcorO2.Time.Mass-qdERLognormSSPlot(x=respData$PerfTemp, Tacc=respData$AccCenter))
ll_lnorm <- sum(log(dlnorm(respData$corO2.Time.Mass+1, qdERLognormSSPlot(x=respData$PerfTemp, Tacc=respData$AccCenter), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=6, loglik=ll_lnorm, n=length(respData$PerfTemp))
lnormAIC

# Final model ANOVA
quad <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, mER=0, qER=0, bThR=300))
lin <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER, qER=0, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, mER=0, bThR=300))
reduced <- nls(corO2.Time.Mass ~ respSS(x=PerfTemp, Tacc=AccCenter, ToR=14.55, EhR, bRTo, mRTo=0, qRTo=0, bER, mER=0, qER=0, bThR, mThR=0, qThR=0), data=respData, start=c(EhR=3, bRTo=.1, bER=.6, bThR=300))
anova(reduced, lin)
anova(lin, quad)
summary(quad)
summary(reduced)