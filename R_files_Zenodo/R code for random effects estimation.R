#{The code snippet proposed here presents a function programmed in R to estimate the EBLUP random effects values in order to calibrate the model to a new dataset. The function is developed for the stem profile model presented in this study, which has 5 fixed parameters and 3 random effects. It can be easily adapted to any model with a different number of fixed and random parameters. Similarly, the variance-covariance matrix of the random effects (R) was considered as a diagonal matrix but more complex forms can be specified. Be aware that this code supposes that the estimation of the random effects can be realized within 20 iterations. Literature says that 20 iterations is enough and this proved true in the examples covered by the study, but checking that the convergence was obtained within the 20 iterations is mandatory â€“and easy enough.}

#{The input dataset profil contains the relative height (hr), relative diameter (dr), tree top height (H) for each tree and plot (plot). The model parameters and the variance components necessary to compute EBLUPS are estimated based on this dataset. The dataframe aps is the calibration dataset containing the (new) data for which the EBLUPS are estimated. It has the same structure and units as profil.}

# 1. The function
# Random effects estimation function
reef <- function(mix, aps){

#--- Create the dataframe D containing variance-covariance terms from the fit
vcvm <- VarCorr(mix, rdig=9)
vcv <- matrix(nrow=3, ncol=3)
vcv[1,1] <- as.numeric(vcvm[1,1])
vcv[2,2] <- as.numeric(vcvm[2,1])
vcv[3,3] <- as.numeric(vcvm[3,1])
vcv[1,2] <- as.numeric(vcvm[2,3])
vcv[1,3] <- as.numeric(vcvm[3,3])
vcv[2,1] <- as.numeric(vcvm[2,3])
vcv[2,3] <- as.numeric(vcvm[3,4])
vcv[3,1] <- as.numeric(vcvm[3,3])
vcv[3,2] <- as.numeric(vcvm[3,4])
D <- vcv

#--- Get the fixed effects parameters
B <- fixef(mix)

#--- Compute the derivates
x <- aps$hr
H_calib <- aps$H
a1 <- B[1]
a2 <- B[2]
a3 <- B[3]
b <- B[4]
c <- B[5]
d <- B[6]
e <- B[7]

der1 <- (a1-(H_calib*a2)/(H_calib+a3))*(-x*(c*exp(-x*d)+1) + x**e)
der2 <- (a1-(H_calib*a2)/(H_calib+a3))*(1-b*x)*exp(-d*x)
der3 <- -(a1-(H_calib*a2)/(H_calib+a3))*c*x*(1-b*x)*exp(-d*x)

Z <- matrix(nrow=dim(aps)[1], ncol=3)
Z[,1] <- der1
Z[,2] <- der2
Z[,3] <- der3

#--- Create the variance-covariance matrix of the random errors
R <- diag(1, nrow=dim(aps)[1], ncol=dim(aps)[1])*(mix$sigma)^2

#--- Initialisation of the EBLUP, u[1,] = FO estimation (zero expectation)
u <- c(0,0,0)
err <- aps$dr - (a1-(H_calib*a2)/(H_calib+a3)) * ((1-(a4+u[1])*x)*(1+(a5+u[2])*exp(-(a6+u[3])*x))-(1-(a4+u[1]))*(x^a7))
u <- D %*% t(Z) %*% solve(Z %*% D %*% t(Z)+R) %*% err

#--- Iterative EBLUP estimation, u[>1,] = FOCE estimation (conditional expectation)
estranef <- data.frame(iter=vector(length=20), u1=vector(length=20), u2=vector(length=20), u3=vector(length=20))
	estranef$iter[1] <- 1
	estranef$u1[1] <- u[1,1]
	estranef$u2[1] <- u[2,1]
	estranef$u3[1] <- u[3,1]

for(k in 2:20){
	err <- aps$dr - (a1-(H_calib*a2)/(Hx_calib+a3)) * ((1-(a4+u[1])*x)*(1+(a5+u[2])*exp(-(a6+u[3])*x))-(1-(a4+u[1]))*(x^a7))
	u <- D %*% t(Z) %*% solve(Z %*% D %*% t(Z)+R) %*% (err + Z %*% u)
	estranef$iter[k] <- k
	estranef$u1[k] <- u[1,1]
	estranef$u2[k] <- u[2,1]
	estranef$u3[k] <- u[3,1]	
}

return(estranef)
}



#2. The use of the function

# Step 1: fit the mixed-effect model

mix <- nlme(dr ~ (a1-(H*a2)/(H+a3)) * ((1-b*hr)*(1+c*exp(-d*hr))-(1-b)*hr**(e)),
		fixed=a1+a2+a3+b+c+d+e~1,
		random= b+c+d~1|plot,
		start=c(a1=1.6, a2=0.7, a3=2.6, b=0.5, c=0.3, d=70, e=8),
		data=profil )

# Step 2: load the application set, run the function, use the random effects estimations

estranef <- reef(mix, aps)

	# Fixed predictions (no EBLUPs)
dpredFix<-(a1-(aps$H*a2)/(aps$H+a3)) * ((1-b*x)*(1+c*exp(-d*x))-(1-b)*(x^e))

	# FO predictions
u1 <- estranef[1,2]
u2 <- estranef[1,3]
u3 <- estranef[1,4]
dpredFO<-(a1-(aps$H*a2)/(aps$H+a3)) * ((1-(b+u1)*x)*(1+(c+u2)*exp(-(d+u3)*x))-(1-(b+u1))*(x^e))

	# FOCE predictions
u1 <- estranef[20,2]
u2 <- estranef[20,3]
u3 <- estranef[20,4]
dpredFOCE<-(a1-(aps$H*a2)/(aps$H+a3)) * ((1-(b+u1)*x)*(1+(c+u2)*exp(-(d+u3)*x))-(1-(b+u1))*(x^e))
