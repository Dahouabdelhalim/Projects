rm(list=ls()) # clear all

source("common_functions.R") # loads common functions

# calculates gamma2 = E[w^(N-1)]/muN
calculate.gamma2 <- function(dist,n,pn,w,muN){
	if(dist == "constant"){
		gamma2 <- w^(muN-1)/muN
	} else {
		gamma2 <- sum(pn*w^n)/(w*muN)
	}
}

calculate.degree.of.cooperation <- function(w,gamma,muN,n,pn,dist,gamma1,gamma2){
	if( (gamma > gamma1) & (gamma > gamma2) ){ # dominant defection; only x0=0 is stable; the degree of cooperation is x0=0
		dc <- 0
	} else if ( (gamma2 < gamma) & (gamma < gamma1) ){ # co-existence; only xF is stable; the degree of cooperation is xF
		dc <- xF_pggds(w,gamma,muN,n,pn,dist)
	} else if ( (gamma < gamma1) & (gamma < gamma2) ){ # dominant cooperation; only x1=1 stable; the degree of cooperation is x1=1
		dc <- 1
	} else if ( (gamma1 < gamma) & (gamma < gamma2) ){ # bi-stability; both x0=0 and x1=1 stable; the degree of cooperation is 1-xF
		dc <- 1 - xF_pggds(w,gamma,muN,n,pn,dist)
	}
	dc
}

dists <- c("constant","poisson","geometric","waring") # distributions

ws <- seq(0.01,2,len=200) # values of w

gammamuNs <- 10^seq(-2,2,len=200) # values of gamma*muN

nmin <- 2 # min group size
nmax <- 100 # max group size
muN <- 5 # mean group size

results <- expand.grid(w=ws,gammamuN=gammamuNs)

results$constant <- rep(NA,dim(results)[[1]])
results$poisson <- rep(NA,dim(results)[[1]])
results$geometric <- rep(NA,dim(results)[[1]])
results$waring <- rep(NA,dim(results)[[1]])

npoints <- dim(results)[[1]]

for(j in 1:length(dists)){
	dist <- dists[j]
	ncol <- j+2
	gamma1 <- 1/muN
	if(dist == "constant"){
		n <- NA
		pn <- NA
	} else {
		param <- find.param(muN,nmin,nmax,dist)
		pn <- truncated.dist(param,nmin,nmax,dist)
		n <- nmin:nmax
	}
	x <- as.matrix(results)
	for(i in 1:npoints){
		w <- x[i,1]
		gammamuN <- x[i,2]
		gamma <- gammamuN/muN
		gamma2 <- calculate.gamma2(dist,n,pn,w,muN)
		x[i,ncol] <- calculate.degree.of.cooperation(w,gamma,muN,n,pn,dist,gamma1,gamma2)
	}
	results <- as.data.frame(x)	
}

save(results, file="results_pggds.RData")

