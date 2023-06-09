library(stats4)
library(sra)


# Most of these functions are similar to those provided in the Selection Response Analysis (sra) R package
# These functions have been modified in order to account for bivariate data
# The algorithm assumes that there is no selection on the second trait (y)


# Very quick introduction to the algorithm:
#
# A model consists in a series of nested functions:
# * A model function (two possibilities: sraCstvar.bivar or sraCstvar.bivar.asym), to be called by the user
# * A time series function (sraTimeSeries.bivar or staTimeSeries.bivar.asym), called by the model function 
#   to simulate a time series as a function of parameter values
# * A full model (all time series) likelyhood function: sraMinuslogL.bivar (common to all models)
# * A time series likelyhood function: sraAutoregMinuslogL.bivar (common to all models)


####################################################
#### AICc for SRA (bivar or univar)
AICc.SRA<-function(SRAobj){
  LL<-logLik(SRAobj)[1]
  Nob<-nrow(SRAobj$data)-2 # number of contrasts taken as N  (if two lines starting form the same P population)
  K<-length(SRAobj$coefficients)
  AIC <- -2*LL + 2*K
  AICc<- AIC + 2*K*(K+1)/(Nob-K-1) 
return(AICc)
}


sraStartingvalues <- function(parameter, sradata, rand=0)
	# Provides stupid, educated, or random guess about starting values. A bit messy, but not important as long as the model converges. 
	# This is a specific function derived from the sra package, but able to deal with bivariate data
{
    "sraStartingvalues.mu0" <- function(sradata, colnm="mean")
        { return(mean(sradata[sradata[,"gen"]==1,colnm]) + rand*rnorm(1,0,abs(mean(sradata[sradata[,"gen"]==1,colnm])))) }
    "sraStartingvalues.logvarA0" <- function(sradata, colnm="var")
        { return(log(0.2*mean(sradata[,colnm]*exp(rand*rnorm(1,0,0.2*mean(sradata[,colnm], na.rm=TRUE))), na.rm=TRUE))) }
    "sraStartingvalues.logvarE0" <- function(sradata, colnm="var")
        { return(log(0.8*mean(sradata[,colnm]*exp(rand*rnorm(1,0,0.8*mean(sradata[,colnm], na.rm=TRUE))), na.rm=TRUE))) }
    
    if (parameter == "mu0")
        { return(sraStartingvalues.mu0(sradata, "mean")) }
    if (parameter == "mu0.A")
		{ return(sraStartingvalues.mu0(sradata, "phen.X.mean")) }
    if (parameter == "mu0.B")
		{ return(sraStartingvalues.mu0(sradata, "phen.Y.mean")) }		
    if (parameter == "logvarA0")
        { return(sraStartingvalues.logvarA0(sradata)) } 
    if (parameter == "logvarA0.pos")
        { return(sraStartingvalues.logvarA0(sradata)) } 
    if (parameter == "logvarA0.neg")
        { return(sraStartingvalues.logvarA0(sradata)) }                 
    if (parameter == "logvarA0.A")
        { return(sraStartingvalues.logvarA0(sradata, "phen.X.var")) }
    if (parameter == "logvarA0.A.pos")
        { return(sraStartingvalues.logvarA0(sradata, "phen.X.var")) }
    if (parameter == "logvarA0.A.neg")
        { return(sraStartingvalues.logvarA0(sradata, "phen.X.var")) }        
    if (parameter == "logvarA0.B")
        { return(sraStartingvalues.logvarA0(sradata, "phen.Y.var")) }                 
    if (parameter == "logvarE0")
        { return(sraStartingvalues.logvarE0(sradata)) }
    if (parameter == "logvarE0.A")
        { return(sraStartingvalues.logvarE0(sradata, "phen.X.var")) }
    if (parameter == "logvarE0.B")
        { return(sraStartingvalues.logvarE0(sradata, "phen.Y.var")) }                
    if (parameter == "covarA0.AB") 
		{ return(0) }
    if (parameter == "covarA0.AB.pos") 
		{ return(0) }    
	if (parameter == "covarA0.AB.neg") 
		{ return(0) }		
    if (parameter == "logvarME")
        { return(log(exp(sraStartingvalues.logvarE0(sradata))/2)) } 
    if (parameter == "logIA0")
        { return(log(exp(sraStartingvalues.logvarA0(sradata)) / (sraStartingvalues.mu0(sradata)^2)))}
    if (parameter == "logIE0")
        { return(log(exp(sraStartingvalues.logvarE0(sradata)) / (sraStartingvalues.mu0(sradata)^2)))}       
    if (parameter == "logith20")
        { a <- exp(sraStartingvalues.logvarA0(sradata)); e <- exp(sraStartingvalues.logvarE0(sradata));
          h2 <- a/(a+e); return(log(h2/(1-h2))) }       
    if (parameter == "logvarP0")
        { return(log(exp(sraStartingvalues.logvarA0(sradata))+exp(sraStartingvalues.logvarE0(sradata)))) }
    if (parameter == "o")
        { return(sraStartingvalues.mu0(sradata)) }
    if (parameter == "s")
        { return(1 + rand*rnorm(1,0,1)) }
    if (parameter == "logepsilon")
        { return(-10 + rand*rnorm(1,0,10)) }
    if (parameter == "logminusepsilon")
        { return(-10 + rand*rnorm(1,0,10)) }
    if (parameter == "logvarepsilon")
        { return(0 + rand*rnorm(1,0,1)) }
    if (parameter == "kc")
        { return(0 + rand*rnorm(1,0,1)) }
    if (parameter == "kg")
        { return(0 + rand*rnorm(1,0,1)) }
    if (parameter == "logvarM")
        { return(-20 + rand*rnorm(1,0,10)) }
    if (parameter == "logNe")
        { return(log(100) + rand*rnorm(1,0,2)) }
    if (parameter == "relativekA0")
        { return(0 + rand*rnorm(1,0,0.5)) }
    if (parameter == "relativekE0")
        { return(0 + rand*rnorm(1,0,0.5)) }
    if (parameter == "kA1")
        { return(1 + rand*rnorm(1,0,0.5)) }
    if (parameter == "kE1")
        { return(1 + rand*rnorm(1,0,0.5)) }
    if (parameter == "kA2")
        { return(0 + rand*rnorm(1,0,0.5)) } 
    if (parameter == "kE2")
        { return(0 + rand*rnorm(1,0,0.5)) } 
    if (parameter == "kA3")
        { return(0 + rand*rnorm(1,0,0.5)) } 
    if (parameter == "kE3")
        { return(0 + rand*rnorm(1,0,0.5)) } 
    if (parameter == "logrelativekA0")
        { return(-1 + rand*rnorm(1,0,1)) }  
    if (parameter == "logrelativekE0")
        { return(-1 + rand*rnorm(1,0,1)) }  
    if (parameter == "logkA1")
        { return(0 + rand*rnorm(1,0,1)) }   
    if (parameter == "logkE1")
        { return(0 + rand*rnorm(1,0,1)) }   
    if (parameter == "logkA2")
        { return(-3 + rand*rnorm(1,0,2)) }
    if (parameter == "logkE2")
        { return(-3 + rand*rnorm(1,0,2)) }  
    if (parameter == "logkA3")
        { return(-3 + rand*rnorm(1,0,2)) }  
    if (parameter == "logkE3")
        { return(-3 + rand*rnorm(1,0,2)) }      
    stop("Unknown parameter ", parameter, ".")
}


sraMakeObject.bivar <- function(sradata, model, start, fixed, FUNtimeseries)
{
	# Produces a summary object containing the data and the model fit. This is enough for
	# plotting, analysis, etc. 
    ans <- list()
    ans$data <- sradata
    ans$model <- model
    ans$start <- start
    ans$vcov <- model@vcov
    ans$coefficients <- coef(model)[names(start)]
    # Confidence intervals: approximated from standard errors
    ans$confint <- cbind(summary(ans$model)@coef[,1]-2*summary(ans$model)@coef[,2], 
            summary(ans$model)@coef[,1]+2*summary(ans$model)@coef[,2])  
    ans$pred <- list()
    ans$residuals <- list()
    ans$vresiduals <- list()
    for(pop in split(ans$data, ans$data$rep))
    {
        range <- 1:(nrow(pop)-1)
        detP <- pop$phen.X.var[range] * pop$phen.Y.var[range] - (pop$phen.XY.covar[range])^2
        s.A <- pop$sel.X.mean[range] - pop$phen.X.mean[range]
        s.B <- pop$sel.Y.mean[range] - pop$phen.Y.mean[range]
        beta.A <- (s.A*pop$phen.Y.var[range] - s.B * pop$phen.XY.covar[range]) / detP
        beta.B <- (s.B*pop$phen.X.var[range] - s.A * pop$phen.XY.covar[range]) / detP
        delta.A <- (pop$sel.X.var[range] - pop$phen.X.var[range])/pop$phen.X.var[range]
        delta.B <- (pop$sel.Y.var[range] - pop$phen.Y.var[range])/pop$phen.Y.var[range]
               
        cc <- c(list(beta.A=beta.A, beta.B=beta.B, delta.A=delta.A, delta.B=delta.B))
        cc <- c(cc, ans$coefficients, unlist(fixed))
        ts <- do.call(what=FUNtimeseries, args=cc)
            
        r <- pop$rep[1] # All should be the same 
        ans$pred[[as.character(r)]] <- list()
        ans$pred[[as.character(r)]]$Gen <- names(ts$mean.A)
        ans$pred[[as.character(r)]]$phen.A <- ts$mean.A
        ans$pred[[as.character(r)]]$phen.B <- ts$mean.B
        ans$pred[[as.character(r)]]$varP.A <- ts$varP.A
        ans$pred[[as.character(r)]]$varP.B <- ts$varP.B
        ans$pred[[as.character(r)]]$covarP <- ts$covarP
        ans$pred[[as.character(r)]]$varA.A <- ts$varA.A
        ans$pred[[as.character(r)]]$varA.B <- ts$varA.B
        
#~         ans$residuals[[as.character(r)]] <- ans$data[ans$data[,"rep"]==as.character(r),"mean"] - ts$mean
#~         ans$vresiduals[[as.character(r)]] <- ans$data[ans$data[,"rep"]==as.character(r),"var"] - ts$varP
    }
    class(ans) <- c("srafit", "list")
    return(ans)
}


############### The real algorithm starts here ##########################

### Functions 
sraMinuslogL.bivar <- function(sradata, FUNtimeseries=sraAutoregTimeseries.bivar, Bulmer=TRUE, ...)
{
	# Computes the likelihood of a full data set (all generations, several repetitions, etc). 
	# This generic function does not need to know about the model parameters, they are passed through the
	# ... argument. 

    ss <- split(sradata, sradata$rep)
    minuslogL <- 0
    for(pop in ss)
    {
        range <- 1:(nrow(pop)-1)
        detP <- pop$phen.X.var[range] * pop$phen.Y.var[range] - (pop$phen.XY.covar[range])^2
        s.A <- pop$sel.X.mean[range] - pop$phen.X.mean[range]
        s.B <- pop$sel.Y.mean[range] - pop$phen.Y.mean[range]
        beta.A <- (s.A*pop$phen.Y.var[range] - s.B * pop$phen.XY.covar[range]) / detP
        beta.B <- (s.B*pop$phen.X.var[range] - s.A * pop$phen.XY.covar[range]) / detP
        delta.A <- (pop$sel.X.var[range] - pop$phen.X.var[range])/pop$phen.X.var[range]
        delta.B <- (pop$sel.Y.var[range] - pop$phen.Y.var[range])/pop$phen.Y.var[range]
                
        tsr <- do.call(what=FUNtimeseries, args=c(list(
            beta.A=beta.A, beta.B=beta.B, delta.A=delta.A, delta.B=delta.B,
            ...)))
        minuslogL <- minuslogL + sraAutoregMinuslogL.bivar(
			X.mean=pop$phen.X.mean, X.var=pop$phen.X.var, 
			Y.mean=pop$phen.Y.mean, Y.var=pop$phen.Y.var, 
			XY.cov=pop$phen.XY.covar, 
			data.N=pop$N, 
			theor.X.mean=tsr$mean.A, theor.Y.mean=tsr$mean.B, 
			theor.X.var=tsr$varP.A, theor.Y.var=tsr$varP.B, 
			theor.XY.cov=tsr$covarP)
    }
    return(minuslogL)       
}

sraAutoregMinuslogL.bivar <- function(X.mean, X.var, Y.mean, Y.var, XY.cov, data.N, 
	theor.X.mean, theor.X.var, theor.Y.mean, theor.Y.var, theor.XY.cov)
	#Â This is the core likelihood function. it computes the likelihood of the bivariate normal
	# distribution. Parameters are the observations (means, variances, covariance) and the
	# theoretical correspondinf values
{
  minuslogL <- 0
  
  for (t in 1:(length(theor.X.mean))) {
    if (is.infinite(theor.X.var[t]) || is.nan(theor.X.var[t]) || (theor.X.var[t] < 0) ||
		is.infinite(theor.Y.var[t]) || is.nan(theor.Y.var[t]) || (theor.Y.var[t] < 0) ||
		is.infinite(theor.XY.cov[t]) || is.nan(theor.XY.cov[t])) {
      return(Inf)
    }
    if(is.na(X.mean[t]) || is.na(X.var[t]) || is.na(Y.mean[t]) || is.na(Y.var[t]) || is.na(XY.cov[t]) || is.na(data.N[t]))
		next
    rho <- theor.XY.cov[t]/(sqrt(theor.X.var[t])*sqrt(theor.Y.var[t]))
    if ((rho < -1) || (rho > 1))
    { 
		return(Inf)
	}
    K <- -data.N[t]*log(2*pi*sqrt(theor.X.var[t]*theor.Y.var[t]*(1-rho^2)))
    minuslogL <- minuslogL - K + data.N[t]*((X.var[t]+(X.mean[t]-theor.X.mean[t])^2)/theor.X.var[t] + (Y.var[t]+(Y.mean[t]-theor.Y.mean[t])^2)/theor.Y.var[t] - 2*rho*(XY.cov[t]+X.mean[t]*Y.mean[t] - theor.X.mean[t]*Y.mean[t] - theor.Y.mean[t]*X.mean[t]+theor.X.mean[t]*theor.Y.mean[t])/(sqrt(theor.X.var[t]*theor.Y.var[t])))/(2*(1-rho^2))
  }
  return(minuslogL)          
}

############ Symmetric response

sraTimeseries.bivar <- function(beta.A, beta.B, delta.A = rep(0, length(beta.A)), delta.B = rep(0, length(beta.B)), 
	mu0.A=0, mu0.B=0, logvarA0.A=0, 
	logvarE0.A=0, logvarE0.B=0, covarA0.AB=0, logNe=log(100), 
	G0.boost=FALSE)
	# Produces a theoretical time series from a set of parameters, and selection strengths (selection on the means --beta--
	# and the variances --delta--)
{
	# The model considers that betaB=0, and does not need the additive genetic variance for trait B
    ans <- list()
        
    Ne <- exp(logNe)

    mu.A <- mu0.A
    mu.B <- mu0.B
    
    d.A <- 0
    vara.A <- exp(logvarA0.A)  # Genic variance
    varA.A <- vara.A + d.A     # Genetic variance (accounting for linkage disequilibrium)
    covarA <- covarA0.AB
        
    varE.A <- exp(logvarE0.A)
    varE.B <- exp(logvarE0.B) 

    for (t in 1:(length(beta.A)))
    {
        V.correct <- if (G0.boost && t == 1) 1.5 else 1
        # Lande's equations
        mu.A <- c(mu.A, mu.A[t] + V.correct * varA.A[t] * beta.A[t]) # + covarA[t] * beta.B[t] : assuming beta.B = 0
        mu.B <- c(mu.B, mu.B[t] + V.correct * covarA[t] * beta.A[t]) # + varA.B[t] * beta.B[t]
    
        vara.A.tp1 <- vara.A[t] * (1 - 1/(2*Ne))
        d.A.tp1 <- 0.5*(1-1/Ne)*(d.A[t]+delta.A[t]*(varA.A[t]**2)/(varA.A[t]+varE.A[t]))
        varA.A.tp1 <- vara.A.tp1 + d.A.tp1
        covarA.tp1 <- covarA[t] * (1 - 1/(2*Ne))
        
        varE.A.tp1 <- varE.A[t]
        varE.B.tp1 <- varE.B[t]

        # To make the exploration of the likelihood surface more robust: handle infinite, zero, or 1/0 variances
        if (is.nan(varA.A.tp1) || is.infinite(varA.A.tp1) || is.na(varA.A.tp1)) 
            {varA.A.tp1 <- 0}
        if (varA.A.tp1 < 0)
            {varA.A.tp1 <- 0}
        if (is.nan(varE.A.tp1) || is.infinite(varE.A.tp1))
            {varE.A.tp1 <- 0} 
        if (is.nan(varE.B.tp1) || is.infinite(varE.B.tp1))
            {varE.B.tp1 <- 0}             
        if (varE.A.tp1 < 0)
            {varE.A.tp1 <- 0 }    
        if (varE.B.tp1 < 0)
            {varE.B.tp1 <- 0 } 
        
        vara.A <- c(vara.A, vara.A.tp1)
        d.A <- c(d.A, d.A.tp1)
        varA.A <- c(varA.A, varA.A.tp1)
        covarA <- c(covarA, covarA.tp1)
        
        varE.A <- c(varE.A, varE.A.tp1)
        varE.B <- c(varE.B, varE.B.tp1)
    }
    
    ans$mean.A <- mu.A
    ans$mean.B <- mu.B
    ans$varA.A <- varA.A
    ans$varA.B <- NA
    ans$covarA <- covarA
    ans$varE.A <- varE.A
    ans$varE.B <- varE.B
    
    ans$varP.A <- ans$varA.A + ans$varE.A
    ans$varP.B <- ans$varE.B # ans$varA.B + ans$varE.B
    ans$covarP <- ans$covarA
    
    return(ans)     
}

sraCstvar.bivar <-function (sradata, start = NULL, fixed = NULL, 
          Bulmer = TRUE, G0.boost=FALSE,...) 
{
  if (!Bulmer) {
    sradata$sel.X.var <- sradata$phen.X.var
    sradata$sel.Y.var <- sradata$phen.Y.var
  }
  default.start <- list(mu0.A = NA, mu0.B=NA, logvarA0.A = NA, logvarE0.A = NA, logvarE0.B = NA, covarA0.AB = NA)
  default.fixed <- list(logNe = log(1e+10))

  default.start[names(fixed)] <- NULL
  default.start[names(start)] <- start
  default.fixed[names(start)] <- NULL
  default.fixed[names(fixed)] <- fixed
  start <- default.start
  fixed <- default.fixed
  start[is.na(start)] <- sapply(names(start[is.na(start)]), 
                                sraStartingvalues, sradata = sradata)
  mlewrapper <- function(mu0.A, mu0.B, logvarA0.A,  logvarE0.A, logvarE0.B, covarA0.AB, logNe) {
    sraMinuslogL.bivar(sradata = sradata, FUNtimeseries = sraTimeseries.bivar, 
                 mu0.A = mu0.A, mu0.B = mu0.B, logvarA0.A = logvarA0.A,
                 logvarE0.A = logvarE0.A, logvarE0.B = logvarE0.B, covarA0.AB = covarA0.AB,
                 logNe = logNe, G0.boost=G0.boost)
  }
  fit <- mle(minuslogl = mlewrapper, start = start, fixed = fixed, 
             ...)
  return(sraMakeObject.bivar(sradata = sradata, model = fit, start = start, 
                       fixed = fixed, FUNtimeseries = sraTimeseries.bivar))
}


############ Asymmetric response

sraTimeseries.bivar.asym <- function(beta.A, beta.B, delta.A = rep(0, length(beta.A)), delta.B = rep(0, length(beta.B)), 
	mu0.A=0, mu0.B=0, logvarA0.A.pos=0, logvarA0.A.neg=0,
	logvarE0.A=0, logvarE0.B=0, covarA0.AB.pos=0, covarA0.AB.neg=0, logNe=log(100),
	G0.boost=FALSE)
{   # Almost the same than for symmetric response, the difference lies in the initial parameter names
    ans <- list()
    
    Ne <- exp(logNe)

    mu.A <- mu0.A
    mu.B <- mu0.B
    d.A <- 0
    vara.A <- if (beta.A[1] > 0) exp(logvarA0.A.pos) else exp(logvarA0.A.neg) 
    varA.A <- vara.A + d.A
    varE.A <- exp(logvarE0.A)
    varE.B <- exp(logvarE0.B) 
    covarA <- if (beta.A[1] > 0) covarA0.AB.pos else covarA0.AB.neg
    
    for (t in 1:(length(beta.A)))
    {
        V.correct <- if (G0.boost && t == 1) 1.5 else 1
        # Lande's equations
        mu.A <- c(mu.A, mu.A[t] + V.correct * varA.A[t] * beta.A[t])
        mu.B <- c(mu.B, mu.B[t] + V.correct * covarA[t] * beta.A[t])
    
		vara.A.tp1 <- vara.A[t] * (1 - 1/(2*Ne))
        d.A.tp1 <- 0.5*(1-1/Ne)*(d.A[t]+delta.A[t]*(varA.A[t]**2)/(varA.A[t]+varE.A[t]))
        covarA.tp1 <- covarA[t] * (1 - 1/(2*Ne))
                
        varA.A.tp1 <- vara.A.tp1 + d.A.tp1
        varE.A.tp1 <- varE.A[t]
        varE.B.tp1 <- varE.B[t]

        
        if (is.nan(varA.A.tp1) || is.infinite(varA.A.tp1)) 
            {varA.A.tp1 <- 0}           
        if (varA.A.tp1 < 0)
            {varA.A.tp1 <- 0}          
        if (is.nan(varE.A.tp1) || is.infinite(varE.A.tp1))
            {varE.A.tp1 <- 0} 
        if (is.nan(varE.B.tp1) || is.infinite(varE.B.tp1))
            {varE.B.tp1 <- 0}             
        if (varE.A.tp1 < 0)
            {varE.A.tp1 <- 0 }    
        if (varE.B.tp1 < 0)
            {varE.B.tp1 <- 0 } 
        
        vara.A <- c(vara.A, vara.A.tp1)
        d.A <- c(d.A, d.A.tp1)
        varA.A <- c(varA.A, varA.A.tp1)
        covarA <- c(covarA, covarA.tp1)
        
        varE.A <- c(varE.A, varE.A.tp1)
        varE.B <- c(varE.B, varE.B.tp1)
    }
    
    ans$mean.A <- mu.A
    ans$mean.B <- mu.B
    ans$varA.A <- varA.A
    ans$varA.B <- NA
    ans$covarA <- covarA
    ans$varE.A <- varE.A
    ans$varE.B <- varE.B
    
    ans$varP.A <- ans$varA.A + ans$varE.A
    ans$varP.B <- ans$varE.B # ans$varA.B + ans$varE.B
    ans$covarP <- ans$covarA
    return(ans)     
}

sraCstvar.bivar.asym <-function (sradata, start = NULL, fixed = NULL, 
          Bulmer = TRUE, G0.boost=FALSE, ...) 
{
  if (!Bulmer) {
    sradata$sel.X.var <- sradata$phen.X.var
    sradata$sel.Y.var <- sradata$phen.Y.var
  }
  default.start <- list(mu0.A = NA, mu0.B=NA, logvarA0.A.pos = NA, logvarA0.A.neg=NA, logvarE0.A = NA, logvarE0.B = NA, covarA0.AB.pos = NA, covarA0.AB.neg = NA)
  default.fixed <- list(logNe = log(1e+10))

  default.start[names(fixed)] <- NULL
  default.start[names(start)] <- start
  default.fixed[names(start)] <- NULL
  default.fixed[names(fixed)] <- fixed
  start <- default.start
  fixed <- default.fixed
  start[is.na(start)] <- sapply(names(start[is.na(start)]), 
                                sraStartingvalues, sradata = sradata)
  mlewrapper <- function(mu0.A, mu0.B, logvarA0.A.pos, logvarA0.A.neg, logvarE0.A, logvarE0.B, covarA0.AB.pos, covarA0.AB.neg, logNe) {
    sraMinuslogL.bivar(sradata = sradata, FUNtimeseries = sraTimeseries.bivar.asym, 
                 mu0.A = mu0.A, mu0.B = mu0.B, 
                 logvarA0.A.pos = logvarA0.A.pos, logvarA0.A.neg = logvarA0.A.neg,
                 logvarE0.A = logvarE0.A, logvarE0.B = logvarE0.B, 
                 covarA0.AB.pos = covarA0.AB.pos, covarA0.AB.neg = covarA0.AB.neg, 
                 logNe = logNe, G0.boost=G0.boost)
  }
  fit <- mle(minuslogl = mlewrapper, start = start, fixed = fixed, 
             ...)
  return(sraMakeObject.bivar(sradata = sradata, model = fit, start = start, 
                       fixed = fixed, FUNtimeseries = sraTimeseries.bivar.asym))
}


