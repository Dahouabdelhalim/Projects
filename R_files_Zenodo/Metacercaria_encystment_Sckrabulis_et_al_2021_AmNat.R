# Metacercaria encystment analysis

# Use Altman et al. 2016 data here as a .csv
# Reference: Altman, K. A., S. H. Paull, P. T. J. Johnson, M. N. Golembieski, J. P. Stephens, B. E. LaFonte, and T. R. Raffel. 2016. Host and parasite thermal acclimation responses depend on the stage of infection. Journal of Animal Ecology 85(4): 1014-1024. doi:10.1111/1365-2656.12510
# Dryad repository doi: 10.5061/dryad.f3k8p
cystData <- read.csv(file.choose(), header=TRUE)
cystData <- cystData[order(cystData$PerfTemp),]
cystData$logAdjPropEncysted <- log(cystData$AdjPropEncysted)
cystData$AccCenter <- cystData$AccTemp-mean(cystData$AccTemp)

# Constants
k<-8.62*10^-5
K<-273.15

#_
# Test upper and lower thresholds on infectivity

# Sharpe-Schoolfield infectivity functions with and without thresholding, using estimates from best fit model in 'Cercaria swimming speed.txt'

# No thresholds
infectNoThresh <- function(x, iTo){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1
}

# Upper threshold (conditional programming) and lower threshold (through the use of translational constant C and conditional programming)
infectThresh <- function(x, iTo, C){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	ifelse(y>1, 1, ifelse(y<0,0,y))
}

# Boltzmann-Arrhenius resistance equation
ba <- function(x, Tor, rTo, Er){
	Tor <- Tor+K
	y <- rTo*exp(-Er/k*(1/(x+K)-1/Tor))
}

baModel <- nls(AdjPropEncysted ~ infectNoThresh(x=PerfTemp, iTo)-ba(x=PerfTemp, Tor=19, rTo, Er),start=c(iTo=.6, rTo=.1, Er=.6), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(0,0,0), upper=c(10,10,10))
baUpper <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo, C=0)-ba(x=PerfTemp, Tor=19, rTo, Er), start=c(iTo=.6, rTo=.1, Er=.6), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(0,0,0), upper=c(10,10,10))
baUpperC <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo, C)-ba(x=PerfTemp, Tor=19, rTo, Er), start=c(iTo=.6, C=0, rTo=.1, Er=.6), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(0,-5,0,0), upper=c(10,10,10,10))

AIC(baModel)
AIC(baUpper)
AIC(baUpperC)


# Plot all threshold models

# Plot parameters
par(mfrow=c(3,1))

# No thresholds

# Interaction plot function
baPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- .8450
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1
	Tor <- 19+K
	rTo <- (.1108)
	Er <- 1.2477
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
	y <- infect-resist
}

# Resistance plot function
resistbaPlot <- function(x){
	Tor <- 19+K
	rTo <- (.1108)
	Er <- 1.2477
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
}

# Infectivity plot function
infectbaPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- .8450
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1
}
plot(cystData$PerfTemp, cystData$AdjPropEncysted, xlim=c(10,40), ylim=c(0,1.5))
curve(baPlot, from=10, to=40, add=TRUE)
curve(resistbaPlot, from=10, to=40, add=TRUE, col="green")
curve(infectbaPlot, from=10, to=40, add=TRUE, col="red")

# Upper threshold

# Interaction plot function
baUpperPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 1.00649
	C <- 0
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.20979)
	Er <- .22307
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
	y <- thresh-resist
}

# Resistance plot function
resistbaUpperPlot <- function(x){
	Tor <- 19+K
	rTo <- (.20979)
	Er <- .22307
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
}

# Infectivity plot function
infectbaUpperPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 1.00649
	C <- 0
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
}
plot(cystData$PerfTemp, cystData$AdjPropEncysted, xlim=c(10,40), ylim=c(0,1.5))
curve(baUpperPlot, from=10, to=40, add=TRUE)
curve(resistbaUpperPlot, from=10, to=40, add=TRUE, col="green")
curve(infectbaUpperPlot, from=10, to=40, add=TRUE, col="red")

# Upper and lower thresholds

# Interaction plot function
baUpperCPlot<-function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.58397
	C <- (1.28218)
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.17965)
	Er <- .4351
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
	y <- thresh-resist
	ifelse(y<0, 0, y)
}

# Resistance plot function
resistbaUpperCPlot<-function(x){
	Tor<-19+K
	rTo<-(.17965)
	Er<-.4351
	resist<-rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
}

# Infectivity plot function
infectbaUpperCPlot<-function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.58397
	C <- (1.28218)
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, ifelse(infect<0, 0, infect))
}
plot(cystData$PerfTemp, cystData$AdjPropEncysted, xlim=c(10,40), ylim=c(0,1.5))
curve(baUpperCPlot, from=10, to=40, add=TRUE)
curve(resistbaUpperCPlot, from=10, to=40, add=TRUE, col="green")
curve(infectbaUpperCPlot, from=10, to=40, add=TRUE, col="red")

#_
# Initial test of normal and lognormal error distribution & model type

# Xiao et al. 2011 AIC calculation
# Reference: Xiao, X., E. P. White, M. B. Hooten, and S. L. Durham. 2011. On the use of log-transformation vs. nonlinear regression for analyzing biological power laws. Ecology 92: 1887-1894. doi:10.1890/11-0538.1
AIC.Xiao <- function(numParams, loglik, n){
	2*numParams-2*loglik
}

#_
# Boltzmann-Arrhenius resistance

# Boltzmann-Arrhenius equation
ba <- function(x, Tor, rTo, Er){
	Tor <- Tor+K
	y <- rTo*exp(-(Er/k*(1/(x+K)-1/(Tor))))
}

# Normal error distribution
baUpperC <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo, C)-ba(x=PerfTemp, Tor=19, rTo, Er), start=c(iTo=.6, C=0, rTo=.1, Er=.6), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(0,-5,0,0), upper=c(10,10,10,10))

# Plotting function using parameter estimates from 'baUpperC' and best fit model in 'Cercaria_swimming_speed_Sckrabulis_et_al_2021_AmNat.R', also used for model predictions
upperCPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.58397
	C <- 1.28218
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.17965)
	Er <- .4351
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
	y <- thresh-resist
	ifelse(y<0, 0, y)
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd((cystData$AdjPropEncysted)-upperCPlot(cystData$PerfTemp))
ll_norm <- sum(log(dnorm(cystData$AdjPropEncysted, upperCPlot(cystData$PerfTemp), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=4, loglik=ll_norm, n=length(cystData$PerfTemp))
normAIC

# Lognormal error distribution
logbaUpperC <- nls(logAdjPropEncysted ~ log(infectThresh(x=PerfTemp, iTo, C)-ba(x=PerfTemp, Tor=19, rTo, Er)), start=c(iTo=.6, C=0, rTo=.1, Er=.6), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(0,-5,0,0), upper=c(10,10,10,10))

# Plotting function using estimates from 'logbaUpperC' and best fit model in 'Cercaria_swimming_speed_Sckrabulis_et_al_2021_AmNat.R', also used for model predictions
logupperCPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.82035
	C <- 1.49456
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.19581)
	Er <- .42382
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
	y <- log(thresh-resist)
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd((cystData$logAdjPropEncysted)-logupperCPlot(cystData$PerfTemp))
ll_lnorm <- sum(log(dlnorm(cystData$AdjPropEncysted, logupperCPlot(cystData$PerfTemp), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=4, loglik=ll_lnorm, n=length(cystData$PerfTemp))
lnormAIC

#_
# Sharpe-Schoolfield

# Sharpe-Schoolfield equation
ss <- function(x, Tor, rTo, Er, Ehr, Thr){
	Tor <- Tor+K
	y <- rTo*exp(-(Er/k*(1/(x+K)-1/(Tor))))*(1+exp(Ehr/k*(1/Thr-1/(x+K))))^-1
}

# Normal error distribution
ssUpperC <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo, C)-ss(x=PerfTemp, Tor=19, rTo, Er, Ehr, Thr), start=c(iTo=1, C=0, rTo=.17, Er=.6, Ehr=3, Thr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-5,-5,-5,0,0,0), upper=c(10,10,10,10,10,400))

# Plotting function using parameter estimates from 'ssUpperC', also used for model predictions
ssUpperCPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.6682
	C <- 1.382
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.1832)
	Er <- .8957
	Ehr <- 1.8584
	Thr <- 301.9951
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr-1/(x+K)))))^-1
	y <- thresh-resist
	ifelse(y<0, 0, y)
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd((cystData$AdjPropEncysted)-ssUpperCPlot(cystData$PerfTemp))
ll_norm <- sum(log(dnorm(cystData$AdjPropEncysted, ssUpperCPlot(cystData$PerfTemp), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=6, loglik=ll_norm, n=length(cystData$PerfTemp))
normAIC

# Lognormal error distribution
logssUpperC <- nls(logAdjPropEncysted ~ log(infectThresh(x=PerfTemp, iTo, C)-ss(x=PerfTemp, Tor=19, rTo, Er, Ehr, Thr)), start=c(iTo=1, C=0, rTo=.17, Er=.6, Ehr=3, Thr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-5,-5,-5,0,0,0), upper=c(10,10,10,10,10,400))

logssUpperCPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.82
	C <- 1.495
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.1958)
	Er <- .423
	Ehr <- 1.791
	Thr <- 332.9
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr-1/(x+K)))))^-1
	y <- log(thresh-resist)
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd((cystData$logAdjPropEncysted)-logssUpperCPlot(cystData$PerfTemp))
ll_lnorm <- sum(log(dlnorm(cystData$AdjPropEncysted, logssUpperCPlot(cystData$PerfTemp), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=6, loglik=ll_lnorm, n=length(cystData$PerfTemp))
lnormAIC

#_
# Model selection - BA

# BA model with acclimation effects

# Infectivity
infectThresh <- function(x,iTo,C){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	ifelse(y>1, 1, y)
}

# Resistance
encystBA <- function(x, Tacc, Tor, brTo, mrTo, bEr, mEr){
	Tor <- Tor+K
	rTo <- function(Tacc, brTo, mrTo){
		rTo <- brTo+(Tacc)*mrTo
	}
	Er <- function(Tacc, bEr, mEr){
		Er <- bEr+(Tacc)*mEr
	}
	y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))
}

# No acclimation effects
NoAcc <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=19, brTo, mrTo=0, bEr, mEr=0), start=c(brTo=.1, bEr=.6), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(0,0), upper=c(10,10))
AIC(NoAcc)

# Linear (l) acclimation effects
# rTo
lrTo <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=19, brTo, mrTo, bEr, mEr=0), start=c(brTo=.1, mrTo=0, bEr=.6), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-5,0,0), upper=c(10,10,10))
AIC(lrTo)

# Er
# Optimize To
# We comment out code here, but omit these notes in all further instances of optimizing To

# Sequence of possible To values over temperature range, we later divide by 10 such that 50/10 = 5C
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))

# For loop to iterate through every possible To value generated above
for(j in range){

	# Encystment interaction model (infect-resist) with acclimation effects
	infect <- function(x, iTo, C){
		Ei <- .40734
		Ehi <- 3.73023
		Thi <- 304.74427
		Toi <- 19+K
		y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
		ifelse(y>=1, 1, y)
	}
	resistMbaAcc <- function(x, Tacc, Tor, brTo, mrTo, bEr, mEr){
		Tor <- Tor+K
		rTo <- function(Tacc, brTo, mrTo){
			rTo <- brTo+mrTo*(Tacc)
		}
		Er <- function(Tacc, bEr, mEr){
			Er<-bEr+mEr*(Tacc)
		}
		y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))
	}
	
	# Model for each To iterated over
	model <- nls(formula=AdjPropEncysted ~ infect(x=PerfTemp, iTo=2.58397, C=1.28218)-resistMbaAcc(x=PerfTemp, Tacc=AccCenter, Tor=(j/10), brTo, mrTo=0, bEr, mEr), start=c(brTo=.4, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-20,0), upper=c(10,10,10))
	vectAIC[j] <- AIC(model)
}
	
# Plot data frame to identify any local minima visually, then search for that area within the data frame for best fit model and corresponding To value
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 19, local minimum
#
lEr <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=19, brTo, mrTo=0, bEr, mEr), start=c(brTo=.1, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-20,0), upper=c(10,10,10))
summary(lEr)
AIC(lEr)

# rTo + Er
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	infect <- function(x, iTo, C){
		Ei <- .40734
		Ehi <- 3.73023
		Thi <- 304.74427
		Toi <- 19+K
		y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
		ifelse(y>=1, 1, y)
	}
	resistMbaAcc <- function(x, Tacc, Tor, brTo, mrTo, bEr, mEr){
		Tor <- Tor+K
		rTo <- function(Tacc, brTo, mrTo){
			rTo <- brTo+mrTo*(Tacc)
		}
		Er <- function(Tacc, bEr, mEr){
			Er <- bEr+mEr*(Tacc)
		}
		y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))
	}
	model <- nls(formula=AdjPropEncysted ~ infect(x=PerfTemp, iTo=2.58397, C=1.28218)-resistMbaAcc(x=PerfTemp, Tacc=AccCenter, Tor=(j/10), brTo, mrTo, bEr, mEr), start=c(brTo=.4, mrTo=0, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-10,-20,0), upper=c(10,10,10,10))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 7.75, local minimum
# To = 33.7
#
lrTo_lEr1 <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=7.75, brTo, mrTo, bEr, mEr), start=c(brTo=.1, mrTo=0, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-10,-20,0), upper=c(10,10,10,10))
summary(lrTo_lEr1)
AIC(lrTo_lEr1)

lrTo_lEr2 <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccTemp, Tor=33.7, brTo, mrTo, bEr, mEr), start=c(brTo=.1, mrTo=0, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-10,-20,0), upper=c(10,10,10,10), trace=TRUE)
summary(lrTo_lEr2)
AIC(lrTo_lEr2)

# Test best fit model normal and lognormal

# Normal error distribution
lEr <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=19, brTo, mrTo=0, bEr, mEr), start=c(brTo=.1, bEr=.6,  mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-20,0), upper=c(10,10,10))

# Plotting and predictions function
lErPlot <- function(x, Tacc){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.6682
	C <- 1.382
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.18001)
	Er <- (.41197)+(Tacc)*.03605
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
	y <- thresh-resist
	ifelse(y<0, 0, y)
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd((cystData$AdjPropEncysted)-lErPlot(x=cystData$PerfTemp,Tacc=cystData$AccCenter))
ll_norm <- sum(log(dnorm(cystData$AdjPropEncysted, lErPlot(x=cystData$PerfTemp, Tacc=cystData$AccCenter), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=3, loglik=ll_norm, n=length(cystData$PerfTemp))
normAIC

# Lognormal error distribution
loglEr <- nls(logAdjPropEncysted ~ log(infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=19, brTo, mrTo=0, bEr, mEr)), start=c(brTo=.1, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-20,0), upper=c(10,10,10))

# Plotting and predictions function
loglErPlot <- function(x,Tacc){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.6682
	C <- 1.382
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.21613)
	Er <- (.28122)+(Tacc)*.0351
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
	y <- log(thresh-resist)
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd((cystData$logAdjPropEncysted)-loglErPlot(x=cystData$PerfTemp, Tacc=cystData$AccCenter))
ll_lnorm <- sum(log(dlnorm(cystData$AdjPropEncysted,loglErPlot(x=cystData$PerfTemp, Tacc=cystData$AccCenter), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=3, loglik=ll_lnorm, n=length(cystData$PerfTemp))
lnormAIC

#_
# Model selection - SS

# Acclimation model

# Infectivity
infectThresh <- function(x, iTo, C){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	ifelse(y>=1, 1, y)
}

# Resistance
encystSS <- function(x, Tacc, Tor, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr){
	Tor <- Tor+K
	rTo <- function(Tacc, brTo, mrTo){
		rTo <- brTo+mrTo*(Tacc)
	}
	Er <- function(Tacc, bEr, mEr){
		Er <- bEr+mEr*(Tacc)
	}
	Thr <- function(Tacc, bThr, mThr){
		Thr <- bThr+mThr*(Tacc)
	}
	y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr(Tacc,bThr,mThr)-1/(x+K)))))^-1
}

# No acclimation effect
noAcc <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=19.2, brTo, mrTo=0, bEr, mEr=0, Ehr, bThr, mThr=0), start=c(brTo=.2, bEr=.6, Ehr=3, bThr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-30,-30,-30,-30), upper=c(20,20,20,400))
summary(noAcc)
AIC(noAcc)

# rTo
lrTo <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=19.2, brTo, mrTo, bEr, mEr=0, Ehr, bThr, mThr=0), start=c(brTo=.5, mrTo=0, bEr=.6, Ehr=3, bThr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-10,0,0,0), upper=c(10,10,10,10,400))
summary(lrTo)
AIC(lrTo)

# Er
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	infectM <- function(x, iTo, C){
		Ei <- .40734
		Ehi <- 3.73023
		Thi <- 304.74427
		Toi <- 19+K
		y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
		ifelse(y>=1, 1, y)
	}
	resistMssAcc <- function(x, Tacc, Tor, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr){
		Tor <- Tor+K
		rTo <- function(Tacc, brTo, mrTo){
			rTo <- brTo+mrTo*(Tacc)
		}
		Er <- function(Tacc, bEr, mEr){
			Er <- bEr+mEr*(Tacc)
		}
		Thr <- function(Tacc, bThr, mThr){
			Thr <- bThr+mThr*(Tacc)
		}
		y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr(Tacc,bThr,mThr)-1/(x+K)))))^-1
	}
	model <- nls(formula=AdjPropEncysted ~ infectM(x=PerfTemp, iTo=2.6682, C=1.382)-resistMssAcc(x=PerfTemp, Tacc=AccCenter, Tor=(j/10), brTo, mrTo=0, bEr, mEr, Ehr, bThr, mThr=0), start=c(brTo=.4, bEr=.6, mEr=0, Ehr=3, bThr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,0,0,0), upper=c(20,20,20,20,400))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To= 19.2
#
lEr <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=19.2, brTo, mrTo=0, bEr, mEr, Edr, bThr, mThr=0), start=c(brTo=.5, bEr=.6, mEr=0, Ehr=3, bThr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,0,0,0), upper=c(20,20,20,20,400))
summary(lEr)
AIC(lEr)

# Thr
lThr <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=19.2, brTo, mrTo=0, bEr, mEr=0, Ehr, bThr, mThr), start=c(brTo=.5, bEr=.6, Ehr=3, bThr=300, mThr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-10,0,0,0), upper=c(10,10,10,400,400))
summary(lThr)
AIC(lThr)

# rTo + Er
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	infectM <- function(x, iTo, C){
		Ei <- .40734
		Ehi <- 3.73023
		Thi <- 304.74427
		Toi <- 19+K
		y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
		ifelse(y>=1, 1, y)
	}
	resistMssAcc<-function(x, Tacc, Tor, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr){
		Tor <- Tor+K
		rTo <- function(Tacc, brTo, mrTo){
			rTo <- brTo+mrTo*(Tacc)
		}
		Er <- function(Tacc, bEr, mEr){
			Er <- bEr+mEr*(Tacc)
		}
		Thr <- function(Tacc, bThr, mThr){
			Thr <- bThr+mThr*(Tacc)
		}
		y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr(Tacc,bThr,mThr)-1/(x+K)))))^-1
	}
	model <- nls(formula=AdjPropEncysted ~ infectM(x=PerfTemp, iTo=2.6682, C=1.382)-resistMssAcc(x=PerfTemp, Tacc=AccCenter, Tor=(j/10), brTo, mrTo, bEr, mEr, Ehr, bThr, mThr=0), start=c(brTo=.4, mrTo=0, bEr=.6, mEr=0, Ehr=3, bThr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,-20,0,0,0), upper=c(20,20,20,20,20,400))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To = 9.3, local minimum
# To = 32.7, local minimum
#
lrTo_lEr_1 <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=9.3, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr=0), start=c(brTo=.5, mrTo=0, bEr=.6, mEr=0, Ehr=3, bThr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,-20,0,0,0), upper=c(20,20,20,20,20,400))
summary(lrTo_lEr_1)
AIC(lrTo_lEr_1)

lrTo_lEr_2 <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=32.7, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr=0), start=c(brTo=.5, mrTo=0, bEr=.6, mEr=0, Ehr=3, bThr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,-20,0,0,0), upper=c(20,20,20,20,20,400))
summary(lrTo_lEr_2)
AIC(lrTo_lEr_2)

# rTo + Thr
lrTo_lThr <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=19.2, brTo, mrTo, bEr, mEr=0, Ehr, bThr, mThr), start=c(brTo=.5, mrTo=0, bEr=.6, Ehr=3, bThr=300, mThr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-10,0,0,0,0), upper=c(10,10,10,10,400,400))
summary(lrTo_lThr)
AIC(lrTo_lThr)

# Er + Thr
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	infectM <- function(x, iTo, C){
		Ei <- .40734
		Ehi <- 3.73023
		Thi <- 304.74427
		Toi <- 19+K
		y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
		ifelse(y>=1, 1, y)
	}
	resistMssAcc <- function(x, Tacc, Tor, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr){
		Tor <- Tor+K
		rTo <- function(Tacc, brTo, mrTo){
			rTo <- brTo+mrTo*(Tacc)
		}
		Er <- function(Tacc, bEr, mEr){
			Er <- bEr+mEr*(Tacc)
		}
		Thr <- function(Tacc, bThr, mThr){
			Thr <- bThr+mThr*(Tacc)
		}
		y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr(Tacc,bThr,mThr)-1/(x+K)))))^-1
	}
	model <- nls(formula=AdjPropEncysted ~ infectM(x=PerfTemp, iTo=2.6682, C=1.382)-resistMssAcc(x=PerfTemp, Tacc=AccCenter, Tor=(j/10), brTo, mrTo=0, bEr, mEr, Ehr, bThr, mThr), start=c(brTo=.4, bEr=.6, mEr=0, Ehr=3, bThr=300, mThr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,-20,0,0,0), upper=c(20,20,20,20,400,400))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To= 19.2 local minimum
#
lEr_lThr <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=19.2, brTo, mrTo=0, bEr, mEr, Ehr, bThr, mThr), start=c(brTo=.5, bEr=.6, mEr=0, Ehr=3, bThr=300, mThr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port" ,lower=c(-20,-20,-20,0,0,0), upper=c(20,20,20,20,400,400))
summary(lEr_lThr)
AIC(lEr_lThr)

# rTo + Er + Thr
range <- seq(from=50, to=400, by=1)
vectAIC <- numeric(length(range))
for(j in range){
	infectM <- function(x, iTo, C){
		Ei <- .40734
		Ehi <- 3.73023
		Thi <- 304.74427
		Toi <- 19+K
		y <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
		ifelse(y>=1, 1, y)
	}
	resistMssAcc <- function(x, Tacc, Tor, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr){
		Tor <- Tor+K
		rTo <- function(Tacc, brTo, mrTo){
			rTo <- brTo+mrTo*(Tacc)
		}
		Er <- function(Tacc, bEr, mEr){
			Er <- bEr+mEr*(Tacc)
		}
		Thr <- function(Tacc, bThr, mThr){
			Thr <- bThr+mThr*(Tacc)
		}
		y <- rTo(Tacc,brTo,mrTo)*exp(-(Er(Tacc,bEr,mEr)/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr(Tacc,bThr,mThr)-1/(x+K)))))^-1
	}
	model <- nls(formula=AdjPropEncysted ~ infectM(x=PerfTemp, iTo=2.6682, C=1.382)-resistMssAcc(x=PerfTemp, Tacc=AccCenter, Tor=(j/10), brTo, mrTo, bEr, mEr, Ehr, bThr, mThr), start=c(brTo=.4, mrTo=0, bEr=.6, mEr=0, Ehr=3, bThr=300, mThr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,-20,-20,0,0,-20), upper=c(20,20,20,20,20,400,400))
	vectAIC[j] <- AIC(model)
}
df <- data.frame(temp=range/10, AIC=vectAIC[vectAIC!=0])
plot(df$temp, df$AIC)
#
# To= 14.7 & 26.7 LOCAL minimum
#
lrTo_lEr_lThr1 <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=14.7, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr), start=c(brTo=.5, mrTo=0, bEr=.6, mEr=0, Ehr=3, bThr=300, mThr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,-20,0,0,0,-20), upper=c(20,20,20,20,20,400,400))
summary(lrTo_lEr_lThr1)
AIC(lrTo_lEr_lThr1)

lrTo_lEr_lThr2 <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.668, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=26.7, brTo, mrTo, bEr, mEr, Ehr, bThr, mThr), start=c(brTo=.5, mrTo=0, bEr=.6, mEr=0, Ehr=3, bThr=300, mThr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,-20,0,0,0,-20), upper=c(20,20,20,20,20,400,400))
summary(lrTo_lEr_lThr2)
AIC(lrTo_lEr_lThr2)

# Test best fit model normal and lognormal error distribution

# Normal error distribution
lEr<-nls(AdjPropEncysted~infectThresh(x=PerfTemp,iTo=2.6683,C=1.382)-encystSS(x=PerfTemp,Tacc=AccCenter,Tor=19.2,brTo,mrTo=0,bEr,mEr,Edr,bThr,mThr=0),
	start=c(brTo=.5,bEr=.6,mEr=0,Edr=3,bThr=300),data=cystData,control=c(warnOnly=TRUE),
	algorithm="port",lower=c(-20,-20,0,0,0),upper=c(20,20,20,20,400))

# Plotting and predictions function
lErPlot <- function(x, Tacc){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.6682
	C <- 1.382
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.17701)
	Er <- (.81712)+(Tacc)*.03888
	Ehr <- 2.53463
	Thr <- 302.54768
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr-1/(x+K)))))^-1
	y <- thresh-resist
	ifelse(y<0, 0, y)
}

# Standard deviation and log-likelihood calculations
sdnorm <- sd((cystData$AdjPropEncysted)-lErPlot(x=cystData$PerfTemp, Tacc=cystData$AccCenter))
ll_norm <- sum(log(dnorm(cystData$AdjPropEncysted, lErPlot(x=cystData$PerfTemp, Tacc=cystData$AccCenter), sdnorm)))

# AIC calculation
normAIC <- AIC.Xiao(numParams=5, loglik=ll_norm, n=length(cystData$PerfTemp))
normAIC

# Lognormal error distribution
loglEr <- nls(logAdjPropEncysted ~ log(infectThresh(x=PerfTemp, iTo=2.6682, C=1.382)-encystSS(x=PerfTemp, Tacc=AccCenter, Tor=19.2, brTo, mrTo=0, bEr, mEr, Ehr, bThr, mThr=0)), start=c(brTo=.5, bEr=.6, mEr=0, Ehr=3, bThr=300), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-20,-20,0,0,0), upper=c(20,20,20,30,400))
loglErPlot <- function(x, Tacc){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.6682
	C <- 1.382
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.2043)
	Er <- (-10.33)+(Tacc+K)*.03674
	Ehr <- 22.21
	Thr <- 301.7
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))*(1+exp((Ehr/k*(1/Thr-1/(x+K)))))^-1
	y <- log(thresh-resist)
}

# Standard deviation and log-likelihood calculations
sdlnorm <- sd((cystData$logAdjPropEncysted)-loglErPlot(x=cystData$PerfTemp, Tacc=cystData$AccTemp))
ll_lnorm <- sum(log(dlnorm(cystData$AdjPropEncysted, loglErPlot(x=cystData$PerfTemp, Tacc=cystData$AccTemp), sdlnorm)))

# AIC calculation
lnormAIC <- AIC.Xiao(numParams=5, loglik=ll_lnorm, n=length(cystData$PerfTemp))
lnormAIC


# Final model ANOVA
full <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=19, brTo, mrTo=0, bEr, mEr), start=c(brTo=.1, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-20,0), upper=c(10,10,10))
reduced <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=19, brTo, mrTo=0, bEr, mEr=0), start=c(brTo=.1, bEr=.6), data=cystData, control=c(warnOnly=TRUE), algorithm="port",lower=c(-10,-20),upper=c(10,10))
anova(reduced, full)
summary(full)
summary(reduced)

#####
# Calculate temperatures at which thresholds are active

# Plotting function
infectCurve <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.58397
	C <- 1.28218
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	y <- ifelse(infect>=1, 1, infect)
	ifelse(y<0, 0, y)
}
curve(infectCurve, from=0, to=40)
speedTemps <- seq(from=0, to=40, by=.1)
infectCurve(x=speedTemps)

# Temperatures at which lower and upper thresholds are active
speedTemps[69]	#6.8
speedTemps[168]	#16.7
speedTemps[322]	#32.1
speedTemps[344]	#34.3

# Swimming speed equation of best fit model from 'Cercaria_swimming_speed_Sckrabulis_et_al_2021_AmNat.R'
speedCurve <- function(x){
	lniTo <- 1.2640962
	Ei <- 0.4073374
	Ehi <- 3.7302859
	Thi <- 304.7443314
	Toi <- 19+K
	y <- lniTo + log(exp(Ei/k*(1/Toi-1/(x+K)))) + log(1/(1+exp(Ehi/k*(1/Thi-1/(x+K)))))
	y
}

# Calculate swimming speed at low temp and high temp thresholds
lowSpeed <- exp(speedCurve(6.8))
highSpeed <- exp(speedCurve(32.1))

#####
# Plot

#_
# Figure 1: Full model plot 6 panels by acclimation temperature

# Libraries
library(Hmisc)

# Use Altman et al. 2016 data here as a .csv
# Reference: Altman, K. A., S. H. Paull, P. T. J. Johnson, M. N. Golembieski, J. P. Stephens, B. E. LaFonte, and T. R. Raffel. 2016. Host and parasite thermal acclimation responses depend on the stage of infection. Journal of Animal Ecology 85(4): 1014-1024. doi:10.1111/1365-2656.12510
# Dryad repository doi: 10.5061/dryad.f3k8p
cystData<-read.csv(file.choose(), header=TRUE)
cystData<-cystData[order(cystData$PerfTemp),]
cystData$logAdjPropEncysted<-log(cystData$AdjPropEncysted)
cystData$AccCenter<-cystData$AccTemp-mean(cystData$AccTemp)

# Final model
full <- nls(AdjPropEncysted ~ infectThresh(x=PerfTemp, iTo=2.58397, C=1.28218)-encystBA(x=PerfTemp, Tacc=AccCenter, Tor=19, brTo, mrTo=0, bEr, mEr), start=c(brTo=.1, bEr=.6, mEr=0), data=cystData, control=c(warnOnly=TRUE), algorithm="port", lower=c(-10,-20,0), upper=c(10,10,10))

# Calculate Er for each acclimation temperature
AccTemps <- c(13,16,19,22,25,28)
AccTempsCenter <- AccTemps-mean(AccTemps)
Er <- (.41197)+(AccTempsCenter)*.03605

# Plotting functions using estimates from 'full' model and best fit model in 'Cercaria_swimming_speed_Sckrabulis_et_al_2021_AmNat.R'

# Interaction model (infect-resist)
fullPlot <- function(x, Tacc){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.58397
	C <- 1.28218
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	thresh <- ifelse(infect>=1, 1, infect)
	Tor <- 19+K
	rTo <- (.18001)
	Er <- (.41197)+(Tacc)*.03605
	resist <- rTo*exp(-(Er/k*(1/(x+K)-1/Tor)))
	y <- thresh-resist
	ifelse(y<0, 0, y)
}

# Infectivity curve
infectPlot <- function(x){
	Ei <- .40734
	Ehi <- 3.73023
	Thi <- 304.74427
	Toi <- 19+K
	iTo <- 2.58397
	C <- 1.28218
	infect <- iTo*exp(-(Ei/k*(1/(x+K)-1/Toi)))*(1+exp((Ehi/k*(1/Thi-1/(x+K)))))^-1-C
	y <- ifelse(infect>=1, 1, infect)
	ifelse(y<0, 0, y)
}

# Resistance curve
resistPlot <- function(x,Tacc){
	Tor <- 19+K
	rTo <- (.18001)
	Er <- (.41197)+(Tacc)*.03605
	y <- rTo*exp(-Er/k*(1/(x+K)-1/Tor))
	ifelse(y>1, 1, y)
}

# Vector setup for plotting means
temps <- c(13,16,19,22,25,28)
accs <- temps-mean(temps)
means <- as.numeric(length(temps))
upps <- as.numeric(length(temps))
lows <- as.numeric(length(temps))

# Plot parameters
par(mfrow=c(3,2), omi=c(.4,.4,0,0))

	#A: Acclimation temperature = 13
	accData <- subset(cystData, AccTemp==temps[1])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	
	mean13 <- mean(perf13$AdjPropEncysted)
	mean16 <- mean(perf16$AdjPropEncysted)
	mean19 <- mean(perf19$AdjPropEncysted)
	mean22 <- mean(perf22$AdjPropEncysted)
	mean25 <- mean(perf25$AdjPropEncysted)
	mean28 <- mean(perf28$AdjPropEncysted)

	se13 <- sd(perf13$AdjPropEncysted)/sqrt(length(perf13$AdjPropEncysted))
	se16 <- sd(perf16$AdjPropEncysted)/sqrt(length(perf16$AdjPropEncysted))
	se19 <- sd(perf19$AdjPropEncysted)/sqrt(length(perf19$AdjPropEncysted))
	se22 <- sd(perf22$AdjPropEncysted)/sqrt(length(perf22$AdjPropEncysted))
	se25 <- sd(perf25$AdjPropEncysted)/sqrt(length(perf25$AdjPropEncysted))
	se28 <- sd(perf28$AdjPropEncysted)/sqrt(length(perf28$AdjPropEncysted))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[1]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", xaxt="n")
	axis(side=1, labels=FALSE)
	curve(infectPlot, from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x,Tacc=accs[1]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	errbar(temps, means, upps, lows, add=TRUE, pch=1)
	text(x=34, y=1, expression(italic("i"[T])), cex=1.5)
	text(x=37.5, y=.35, expression(italic("r"[T])), cex=1.5)
	text(x=23, y=.7, expression(italic(lambda[T])), cex=1.5)
	mtext(expression(paste("A: ",italic(T[Accl])," = 13",degree,"C")), side=3, adj=.05, line=-2, cex=1)
	
	#B: Acclimation temperature = 16
	accData <- subset(cystData, AccTemp==temps[2])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	
	mean13 <- mean(perf13$AdjPropEncysted)
	mean16 <- mean(perf16$AdjPropEncysted)
	mean19 <- mean(perf19$AdjPropEncysted)
	mean22 <- mean(perf22$AdjPropEncysted)
	mean25 <- mean(perf25$AdjPropEncysted)
	mean28 <- mean(perf28$AdjPropEncysted)

	se13 <- sd(perf13$AdjPropEncysted)/sqrt(length(perf13$AdjPropEncysted))
	se16 <- sd(perf16$AdjPropEncysted)/sqrt(length(perf16$AdjPropEncysted))
	se19 <- sd(perf19$AdjPropEncysted)/sqrt(length(perf19$AdjPropEncysted))
	se22 <- sd(perf22$AdjPropEncysted)/sqrt(length(perf22$AdjPropEncysted))
	se25 <- sd(perf25$AdjPropEncysted)/sqrt(length(perf25$AdjPropEncysted))
	se28 <- sd(perf28$AdjPropEncysted)/sqrt(length(perf28$AdjPropEncysted))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[2]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", xaxt="n", yaxt="n")
	axis(side=1, labels=FALSE)
	axis(side=2, labels=FALSE)
	curve(infectPlot, from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[2]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	errbar(temps, means, upps, lows, add=TRUE, pch=1)
	mtext(expression(paste("B: ",italic(T[Accl])," = 16",degree,"C")), side=3, adj=.05, line=-2, cex=1)
	
	#C: Acclimation temperature = 19
	accData <- subset(cystData, AccTemp==temps[3])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	
	mean13 <- mean(perf13$AdjPropEncysted)
	mean16 <- mean(perf16$AdjPropEncysted)
	mean19 <- mean(perf19$AdjPropEncysted)
	mean22 <- mean(perf22$AdjPropEncysted)
	mean25 <- mean(perf25$AdjPropEncysted)
	mean28 <- mean(perf28$AdjPropEncysted)

	se13 <- sd(perf13$AdjPropEncysted)/sqrt(length(perf13$AdjPropEncysted))
	se16 <- sd(perf16$AdjPropEncysted)/sqrt(length(perf16$AdjPropEncysted))
	se19 <- sd(perf19$AdjPropEncysted)/sqrt(length(perf19$AdjPropEncysted))
	se22 <- sd(perf22$AdjPropEncysted)/sqrt(length(perf22$AdjPropEncysted))
	se25 <- sd(perf25$AdjPropEncysted)/sqrt(length(perf25$AdjPropEncysted))
	se28 <- sd(perf28$AdjPropEncysted)/sqrt(length(perf28$AdjPropEncysted))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[3]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", xaxt="n")
	axis(side=1, labels=FALSE)
	curve(infectPlot, from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[3]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	errbar(temps, means, upps, lows, add=TRUE, pch=1)
	mtext(expression(paste("C: ",italic(T[Accl])," = 19",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	#D: Acclimation temperature = 22
	accData <- subset(cystData, AccTemp==temps[4])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	
	mean13 <- mean(perf13$AdjPropEncysted)
	mean16 <- mean(perf16$AdjPropEncysted)
	mean19 <- mean(perf19$AdjPropEncysted)
	mean22 <- mean(perf22$AdjPropEncysted)
	mean25 <- mean(perf25$AdjPropEncysted)
	mean28 <- mean(perf28$AdjPropEncysted)

	se13 <- sd(perf13$AdjPropEncysted)/sqrt(length(perf13$AdjPropEncysted))
	se16 <- sd(perf16$AdjPropEncysted)/sqrt(length(perf16$AdjPropEncysted))
	se19 <- sd(perf19$AdjPropEncysted)/sqrt(length(perf19$AdjPropEncysted))
	se22 <- sd(perf22$AdjPropEncysted)/sqrt(length(perf22$AdjPropEncysted))
	se25 <- sd(perf25$AdjPropEncysted)/sqrt(length(perf25$AdjPropEncysted))
	se28 <- sd(perf28$AdjPropEncysted)/sqrt(length(perf28$AdjPropEncysted))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[4]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", xaxt="n", yaxt="n")
	axis(side=1, labels=FALSE)
	axis(side=2, labels=FALSE)
	curve(infectPlot, from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[4]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	errbar(temps, means, upps, lows, add=TRUE, pch=1)
	mtext(expression(paste("D: ",italic(T[Accl])," = 22",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	#E: Acclimation temperature = 25
	accData <- subset(cystData, AccTemp==temps[5])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	
	mean13 <- mean(perf13$AdjPropEncysted)
	mean16 <- mean(perf16$AdjPropEncysted)
	mean19 <- mean(perf19$AdjPropEncysted)
	mean22 <- mean(perf22$AdjPropEncysted)
	mean25 <- mean(perf25$AdjPropEncysted)
	mean28 <- mean(perf28$AdjPropEncysted)

	se13 <- sd(perf13$AdjPropEncysted)/sqrt(length(perf13$AdjPropEncysted))
	se16 <- sd(perf16$AdjPropEncysted)/sqrt(length(perf16$AdjPropEncysted))
	se19 <- sd(perf19$AdjPropEncysted)/sqrt(length(perf19$AdjPropEncysted))
	se22 <- sd(perf22$AdjPropEncysted)/sqrt(length(perf22$AdjPropEncysted))
	se25 <- sd(perf25$AdjPropEncysted)/sqrt(length(perf25$AdjPropEncysted))
	se28 <- sd(perf28$AdjPropEncysted)/sqrt(length(perf28$AdjPropEncysted))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[5]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="")
	curve(infectPlot, from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[5]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	errbar(temps, means, upps, lows, add=TRUE, pch=1)
	mtext(expression(paste("E: ",italic(T[Accl])," = 25",degree,"C")), side=3, adj=.05, line=-2, cex=1)

	#F: Acclimation temperature = 28
	accData <- subset(cystData, AccTemp==temps[6])
	perf13 <- subset(accData, PerfTemp==13)
	perf16 <- subset(accData, PerfTemp==16)
	perf19 <- subset(accData, PerfTemp==19)
	perf22 <- subset(accData, PerfTemp==22)
	perf25 <- subset(accData, PerfTemp==25)
	perf28 <- subset(accData, PerfTemp==28)
	
	mean13 <- mean(perf13$AdjPropEncysted)
	mean16 <- mean(perf16$AdjPropEncysted)
	mean19 <- mean(perf19$AdjPropEncysted)
	mean22 <- mean(perf22$AdjPropEncysted)
	mean25 <- mean(perf25$AdjPropEncysted)
	mean28 <- mean(perf28$AdjPropEncysted)

	se13 <- sd(perf13$AdjPropEncysted)/sqrt(length(perf13$AdjPropEncysted))
	se16 <- sd(perf16$AdjPropEncysted)/sqrt(length(perf16$AdjPropEncysted))
	se19 <- sd(perf19$AdjPropEncysted)/sqrt(length(perf19$AdjPropEncysted))
	se22 <- sd(perf22$AdjPropEncysted)/sqrt(length(perf22$AdjPropEncysted))
	se25 <- sd(perf25$AdjPropEncysted)/sqrt(length(perf25$AdjPropEncysted))
	se28 <- sd(perf28$AdjPropEncysted)/sqrt(length(perf28$AdjPropEncysted))

	means <- c(mean13, mean16, mean19, mean22, mean25, mean28)
	upps <- c(mean13+se13, mean16+se16, mean19+se19, mean22+se22, mean25+se25, mean28+se28)
	lows <- c(mean13-se13, mean16-se16, mean19-se19, mean22-se22, mean25-se25, mean28-se28)

	# Panel parameters
	par(mai=c(.2,.2,.2,.2))
	
	curve(fullPlot(x, Tacc=accs[6]), from=10, to=40, ylim=c(0,1.2), lwd=2, xlab="", ylab="", yaxt="n")
	axis(side=2, labels=FALSE)
	curve(infectPlot, from=10, to=40, add=TRUE, lty=4, lwd=2, col="grey30")
	curve(resistPlot(x, Tacc=accs[6]), from=10, to=40, add=TRUE, lty=2, lwd=2)
	errbar(temps, means, upps, lows, add=TRUE, pch=1)
	mtext(expression(paste("F: ",italic(T[Accl])," = 28",degree,"C")), side=3, adj=.05, line=-2, cex=1)
	
	mtext(expression(paste("Performance temperature (",degree,"C)")), side=1, line=1.5, outer=TRUE)
	mtext("Proportion of encysted metacercariae", side=2, line=1.5, outer=TRUE)