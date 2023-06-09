# Activation_energy_bootstrap_Sckrabulis_et_al_2021_AmNat.R

# Bootstrapped activation energy (Ea) acclimation effects
# Predictions and confidence intervals generated here were combined into a single .csv file (we make our generated dataset available as 'Activation energy bootstrap - Sckrabulis et al 2021 AmNat.csv', but due to random sampling that data cannot be exactly reproduced)

# Clearance Data (reduced to 4 columns for centered Acc effect: 1=PerfTemp, 2=Response, 3=AccTemp, 4=AccCenter) for use in bootstrap function
# Use Altman et al. 2016 data here as a .csv
# Reference: Altman, K. A., S. H. Paull, P. T. J. Johnson, M. N. Golembieski, J. P. Stephens, B. E. LaFonte, and T. R. Raffel. 2016. Host and parasite thermal acclimation responses depend on the stage of infection. Journal of Animal Ecology 85(4): 1014-1024. doi:10.1111/1365-2656.12510
# Dryad repository doi: 10.5061/dryad.f3k8p
propClearData2 <- read.csv(file.choose(), header=TRUE)
propClearData2 <- subset(propClearData2, PropCleared1!="NA")
propClearData2$AccCenter <- propClearData2$AccTemp-mean(propClearData2$AccTemp)
reducedClear <- cbind(perf=propClearData2$PerfTemp, response=propClearData2$PropCleared1, acc=propClearData2$AccTemp, center=propClearData2$AccCenter)

# Respiration (reduced to 4 columns for centered Acc effect: 1=PerfTemp, 2=Response, 3=AccTemp, 4=AccCenter) for use in bootstrap function
# Use 'Uninfected_tadpole_respiration_Sckrabulis_et_al_2021_AmNat.csv' in this repository
respData2 <- read.csv(file.choose(), header=TRUE)
respData2 <- subset(respData2, AccTemp!="NA")
respData2$AccCenter <- respData2$AccTemp-mean(respData2$AccTemp)
reducedResp <- cbind(perf=respData2$PerfTemp, response=respData2$corO2.Time.Mass, acc=respData2$AccTemp, center=respData2$AccCenter)

# Function to bootstrap SS model data for quadratic acclimation effects (centered) on Ea with the arguments:
# data = dataset to bootstrap
# iterations = number of iterations for bootstrap sampling
# To = value for T0 from optimized final model
# starting = vector of starting values for parameters

bootEaCentered <- function(data, iterations, To, starting){
	plot(data[,1], data[,2])
	
	# Vector of x values for the entire range of AccTemp
	ciX<-seq(from=min(data[,4]), to=max(data[,4]),by=.1)

	# Number of columns for predictions and for complete dataframe (3 parameters in quadratic model)
	colPredict <- length(ciX)
	colTotal <- colPredict+3
	n <- length(data[,1])
	
	# Empty dataframe for storing parameters and predictions
	params <- data.frame(matrix(ncol=colTotal, nrow=iterations))
	conf_interval <- data.frame(matrix(ncol=3, nrow=length(ciX)))
	names(params) <- c("b", "m", "q")

	# NLS model equation as a function
	# We replace 'k' from all other .R files with 'boltz' here to facilitate for() loop nomenclature
	quadSS <- function(x, Tacc, A, bE, mE, qE, Th){
		K <- 273.15
		boltz <- 8.62*10^-5
		T0 <- To+K
		Ed <- 3.25
		Ea <- function(Tacc,bE,mE,qE){
			Ea <- bE+mE*(Tacc)+qE*(Tacc)^2
		}
		y <- A*exp(-Ea(Tacc,bE,mE,qE)/boltz*(1/(x+K)-1/T0))*(1+exp(Ed/boltz*(1/Th-1/(x+K))))^-1
	}

	# For loop to iterate bootstrap and predictions until X iterations are complete (accounting for convergence errors)
	for(i in 1:iterations){
		skip <- FALSE
		
		# Empty dataframe and generate bootstrap data with replacement
		bootData <- data.frame(matrix(ncol=4, nrow=n))
		bootRows <- sample(1:length(data[,1]), size=length(data[,1]), replace=TRUE)

		# For loop to insert data into bootData for modeling in next step
		for(j in 1:length(bootRows)){
			bootData[j,] <- data[bootRows[j],]
		}
		names(bootData) <- c("perf", "response", "acc", "center")
		plot(bootData[,1], bootData[,2])

		# Wrap model fitting in error catching function to circumvent fatal convergence errors and warnings
		attempt <- tryCatch(
			{
				model <- nls(response ~ quadSS(x=perf, Tacc=center, A, bE, mE, qE, Th), start=starting, data=bootData, control=c(warnOnly=TRUE))
				summary(model)
			},
			error = function(e){
				skip <<- TRUE
			}
		)
			
		if(skip == TRUE){
			params[i,] <- "NA"
			next
		} else if(skip == FALSE){
			model <- nls(response ~ quadSS(x=perf, Tacc=center, A, bE, mE, qE, Th), start=starting, data=bootData, control=c(warnOnly=TRUE))
			params[i,1] <- summary(model)$coef[2]
			params[i,2] <- summary(model)$coef[3]
			params[i,3] <- summary(model)$coef[4]
			
			# Predict Ea curve over entire range of ciX 
			for(k in 1:length(ciX)){
				params[i,k+3] <- summary(model)$coef[2]+summary(model)$coef[3]*(ciX[k])+summary(model)$coef[4]*(ciX[k])^2
			}
		}
	}

	# Calculate confidence intervals for each new value in ciX
	for(p in 1:colPredict){
		newParams <- subset(params, params$b!="NA")
		values <- as.numeric(newParams[,p+3])
		lowCI <- sort(values)[trunc(0.025*length(newParams[,p+3]))]
		highCI <- sort(values)[length(newParams[,p+3])-trunc(0.025*length(newParams[,p+3]))]
		conf_interval[p,1] <- ciX[p]
		conf_interval[p,2] <- lowCI
		conf_interval[p,3] <- highCI
	}
	
	# Plot best fit model Ea equation and confidence bands
	raw <- as.data.frame(data)
	model2 <- nls(response ~ quadSS(x=perf, Tacc=center, A, bE, mE, qE, Th), start=starting, data=raw, control=c(warnOnly=TRUE))
	int2 <- summary(model2)$coef[2]
	line2 <- summary(model2)$coef[3]
	quad2 <- summary(model2)$coef[4]

	# Predictions for full Ea curve over all values of ciX
	fullEa <- numeric(length(ciX))
	for(q in 1:length(ciX)){
		fullEa[q] <- int2+line2*ciX[q]+quad2*ciX[q]^2
	}

	plot(ciX+mean(raw[,3]), fullEa,xlim=c(10,30), ylim=c(0,2), col="lightblue", type="l")
	lines(conf_interval[,1]+mean(raw[,3]), conf_interval[,2], col="blue", lty=2)
	lines(conf_interval[,1]+mean(raw[,3]), conf_interval[,3], col="blue", lty=2)
	
	df <- as.data.frame(cbind(AccTemp=ciX+mean(raw[,3]), Pred=fullEa, low=conf_interval[,2], high=conf_interval[,3]))

	# Save data for later combined plot (replace your_path_here with complete file path and name you desire)
	write.csv(df,"your_path_here", row.names=FALSE)
}

# Running the function may result in instances of failed convergence and will display a 'subscript out of bounds' error. This is a result of the tryCatch() function in line 73 and that failed model is excluded from the end result.
bootEaCentered(data=reducedClear, iterations=10000, To=26.9, starting=c(A=.4, bE=.6, mE=0, qE=0, Th=303))

#####
# Our final model for tadpole respiration contained Eh as a free parameter, but allowing Eh to vary while bootstrapping caused too many errors and unrealistic (i.e., negative) Eh estimates
# Therefore we fixed Eh = 2.894 (the estimate from the best fit model)
# The following function is identical to the above bootEaCentered() but reflects this change
# Code comments in this function were omitted here as they are identical to those above
bootEaCentered2 <- function(data, iterations, To, starting){
	plot(data[,1], data[,2])
	
	ciX <- seq(from=min(data[,4]), to=max(data[,4]), by=.1)

	colPredict <- length(ciX)
	colTotal <- colPredict+3
	n <- length(data[,1])
	
	params <- data.frame(matrix(ncol=colTotal, nrow=iterations))
	conf_interval <- data.frame(matrix(ncol=3, nrow=length(ciX)))
	names(params) <- c("b", "m", "q")

	quadSS <- function(x, Tacc, A, bE, mE, qE, Th){
		K <- 273.15
		boltz <- 8.62*10^-5
		T0 <- To+K
		Ed <- 2.894
		Ea <- function(Tacc,bE,mE,qE){
			Ea <- bE+mE*(Tacc)+qE*(Tacc)^2
		}
		y <- A*exp(-Ea(Tacc,bE,mE,qE)/boltz*(1/(x+K)-1/T0))*(1+exp(Ed/boltz*(1/Th-1/(x+K))))^-1
	}

	for(i in 1:iterations){
		skip <- FALSE

		bootData <- data.frame(matrix(ncol=4, nrow=n))
		bootRows <- sample(1:length(data[,1]), size=length(data[,1]), replace=TRUE)

		for(j in 1:length(bootRows)){
			bootData[j,] <- data[bootRows[j],]
		}
		names(bootData) <- c("perf", "response", "acc", "center")
		plot(bootData[,1], bootData[,2])

		attempt <- tryCatch(
			{
				model <- nls(response ~ quadSS(x=perf, Tacc=center, A, bE, mE, qE, Th), start=starting, data=bootData, control=c(warnOnly=TRUE))
				summary(model)
			},
			error = function(e){
				skip <<- TRUE
			}
		)
			
		if(skip == TRUE){
			params[i,] <- "NA"
			next
		} else if(skip == FALSE){
			model <- nls(response ~ quadSS(x=perf, Tacc=center, A, bE, mE, qE, Th), start=starting, data=bootData, control=c(warnOnly=TRUE))
			params[i,1] <- summary(model)$coef[2]
			params[i,2] <- summary(model)$coef[3]
			params[i,3] <- summary(model)$coef[4]
			
			for(k in 1:length(ciX)){
				params[i,k+3] <- summary(model)$coef[2]+summary(model)$coef[3]*(ciX[k])+summary(model)$coef[4]*(ciX[k])^2
			}
		}
	}

	for(p in 1:colPredict){
		newParams <- subset(params, params$b!="NA")
		values <- as.numeric(newParams[,p+3])
		lowCI <- sort(values)[trunc(0.025*length(newParams[,p+3]))]
		highCI <- sort(values)[length(newParams[,p+3])-trunc(0.025*length(newParams[,p+3]))]
		conf_interval[p,1] <- ciX[p]
		conf_interval[p,2] <- lowCI
		conf_interval[p,3] <- highCI
	}
	
	raw <- as.data.frame(data)
	model2 <- nls(response ~ quadSS(x=perf, Tacc=center, A, bE, mE, qE, Th), start=starting, data=raw, control=c(warnOnly=TRUE))
	int2 <- summary(model2)$coef[2]
	line2 <- summary(model2)$coef[3]
	quad2 <- summary(model2)$coef[4]

	fullEa <- numeric(length(ciX))
	for(q in 1:length(ciX)){
		fullEa[q] <- int2+line2*ciX[q]+quad2*ciX[q]^2
	}

	plot(ciX+mean(raw[,3]), fullEa, xlim=c(10,30), ylim=c(0,2), col="lightblue", type="l")
	lines(conf_interval[,1]+mean(raw[,3]), conf_interval[,2], col="blue", lty=2)
	lines(conf_interval[,1]+mean(raw[,3]), conf_interval[,3], col="blue", lty=2)
	
	df <- as.data.frame(cbind(AccTemp=ciX+mean(raw[,3]), Pred=fullEa, low=conf_interval[,2], high=conf_interval[,3]))

	write.csv(df,"your_path_here2",row.names=FALSE)
}

bootEaCentered2(data=reducedResp, iterations=10000, To=14.55, starting=c(A=.1,bE=.6,mE=0,qE=0,Th=300))