### shifted distributions:

##################################################################
##################################################################
### Gamma
rsgamma <- function(n, shape, rate = 1, x0=0, ...) { 
  x0 + rgamma(n, shape = shape, rate = rate, ...) 
} 

dsgamma <- function(x, shape, rate = 1, log = FALSE, x0=0, ...) { 
  dgamma(x-x0, shape = shape, rate = rate, log = log, ...) 
} 

psgamma <- function(x, shape, rate = 1, log = FALSE, x0=0, ...) { 
  pgamma(x-x0, shape = shape, rate = rate, log = log, ...) 
}


##################################################################
##################################################################
### Exponential

rsexp <- function(n, rate = 1, x0=0) { 
  x0 + rexp(n, rate=rate) 
} 

dsexp <- function(x, rate = 1, log = FALSE, x0=0) { 
  dexp(x-x0, rate=rate, log=log) 
} 

psexp <- function(x, rate = 1, log = FALSE, x0=0) { 
  pexp(x-x0, rate=rate, log=log) 
}


MLEsexp <-function(x, x0=0){

    return(list(lam.hat=1/mean(x-x0)))

}


##################################################################
##################################################################
### Pareto -- modified from the "PtProcess" library

rpareto<-function(n, lambda, a) 
{
    if (any(lambda <= 0)) 
        stop("All values of lambda must be positive")
    if (length(lambda) != 1 & length(lambda) != n) 
        stop("Argument lambda is wrong length")
    if (a <= 0) 
        stop("Argument a must be positive")
    return(a * exp(rexp(n, rate = lambda)))
}

ppareto<-function(q, lambda, a, lower.tail = TRUE, log.p = FALSE) 
{
    if (any(lambda <= 0)) 
        stop("All values of lambda must be positive")
    if (length(lambda) != 1 & length(lambda) != length(q)) 
        stop("Argument lambda is wrong length")
    if (a <= 0) 
        stop("Argument a must be positive")
    p <- 1 - (a/q)^lambda
    if (any(q < a)) p[which(q < a)] <- 0

    if (lower.tail == FALSE) 
        p <- 1 - p
    if (log.p == TRUE) 
        p <- log(p)
    return(p)
}

dpareto<-function(x, lambda, a, log = FALSE) 
{
    if (any(lambda <= 0)) 
        stop("All values of lambda must be positive")
    if (length(lambda) != 1 & length(lambda) != length(x)) 
        stop("Argument lambda is wrong length")
    if (a <= 0) 
        stop("Argument a must be positive")
    d <- lambda/a * (a/x)^(lambda + 1)
    if (any(x < a))  d[which(x < a)] <- 0
    if (log == TRUE) 
        d <- log(d)
    return(d)
}

MLEpareto<-function(x, a=1){

    if(is.null(a) | is.na(a)){
        a.hat<-min(x)
    }else a.hat<-a

    lam.hat<-length(x)/(sum(log(x)-log(a.hat)))

    return(par=list(lam.hat=lam.hat, a.hat=a.hat))
    
    
}

MLEpareto.counts<-function(dat, a=1){

    if(is.null(a) | is.na(a)){
        a.hat<-min(x)
    }else a.hat<-a

    x<-NULL
    js<-unique(dat$J)
    for(i in length(js)) x<-c(x, rep(js[i], length=dat$count[i]))

    lam.hat<-length(x)/(sum(log(x)-log(a.hat)))

    return(par=list(lam.hat=lam.hat, a.hat=a.hat))
    
    
}


nlog.lik.pareto<-function(obs, params, dist="pareto", mi=30){

    print(paste("params = ", round(params, digits=4), sep=" "))

    ddist<-eval(as.name(paste("d", dist, sep="")))

    formals(ddist)[2] <- params[1]
    formals(ddist)$a <- mi
    
    temp<-ddist(obs, log=TRUE)
    
    return(-sum(temp))
	
}





##################################################################
##################################################################
### Tsallis q-exponential, implemented based on wikipedia info

qE<-function(x, q){

    if(q!=1){
        (1+(1-q)*x)^(1/(1-q))
    }else exp(x)
    
}

qLn<-function(x, q){
    if(q==1){
        log(x)
    }else{
        (x^(1-q) - 1)/(1-q)
    }    

}

rsqexp<-function(n, shape = 1, rate=1, x0=0){
  
    q <- shape
    lam<-rate
    if(q<1 | q>=2) stop("shape parameter (q) must be 1 <= q < 2")
    if(lam<0) stop("rate parameter (lambda) must be positive")
    ss<-runif(n)
    qq<-1/(2-q)
    return( x0 + -qq*qLn(ss, qq)/lam )
  
}


psqexp<-function(x, shape = 1, rate=1, x0=0, lower.tail=TRUE, log.p=FALSE){
  
    q <- shape
    lam<-rate
    xx<-x-x0
    
    if(q<1 | q>=2) stop("shape parameter (q) must be 1 <= q < 2")
    if(lam<0) stop("rate parameter (lambda) must be positive")
    
    
    qq<-1/(2-q)
    
    p <- 1 - qE(-lam*xx/qq, qq)
    if (any(xx < 0)) p[which(xx < 0)] <- 0
    
    if (lower.tail == FALSE) 
      p <- 1 - p
    if (log.p == TRUE) 
      p <- log(p)
    return(p)
  
}


dsqexp<-function(x, shape = 1, rate=1, log=FALSE, x0=0){
  
    q <- shape
    lam<-rate
    xx<-x-x0
    if(q<1 | q>=2) stop("shape parameter (q) must be 1 <= q < 2")
    if(lam<0) stop("rate parameter (lambda) must be positive")
    
    
    d <-  (2-q)*lam * qE(-lam*xx, q)
    if (any(xx < 0))  d[which(xx < 0)] <- 0
    if (log == TRUE) 
      d <- log(d)
    return(d)
  
}


##################################################################
##################################################################
## likelihoods for the data w/o an observation process

nlog.lik<-function(obs, params, dist="gamma", unitconv=NULL, w.conv=1, ...){

    print(paste("params = ", round(params, digits=4), sep=" "))

    if(!is.null(unitconv)) params[w.conv]<-params[w.conv]/unitconv

    ddist<-eval(as.name(paste("d", dist, sep="")))

    formals(ddist)[2] <- params[1]
    if(length(params) == 2) formals(ddist)[3] <- params[2]
    
    temp<-ddist(obs, log=TRUE, ...)

    return(-sum(temp))
	
}

