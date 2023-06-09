library(SuppDists) # (for function dghyper)

# calculate the pmf of group size N with support n in {nmin, nmin+1, ... nmax}
# param is the free parameter of the distribution dist, i.e.
# lambda for pois, p for geom and a for waring
truncated.dist <- function(param,nmin,nmax,dist){
	k <- 0:(nmax-nmin)
	if(dist == "poisson"){
		phik <- dpois(x=k,lambda=param)
	} else if(dist == "geometric"){
		phik <- dgeom(x=k,p=param)		
	} else if(dist == "waring"){
		phik <- dghyper(x=k,a=-1,k=-param,N=1)
	}
	pn <- phik/sum(phik)
	pn
}

# function to optimize (the zero of f.opt is the value of param such that sum(n*pn) = muN)
f.opt <- function(param,muN,nmin,nmax,dist){
	n <- nmin:nmax
	pn <- truncated.dist(param,nmin,nmax,dist)
	f.opt <- sum(pn*n) - muN
	f.opt
} 

# numerically finds param so that sum(n*pn) = muN
find.param <- function(muN,nmin,nmax,dist){
	# initialize interval (parammin,parammax)
	parammin <- .Machine$double.eps
	if(dist == "poisson" | dist == "waring"){
		parammax <- nmax
	} else if(dist == "geometric"){
		parammax <- 1-.Machine$double.eps		
	}
	zero <- uniroot(f.opt,interval=c(parammin,parammax),muN=muN,nmin=nmin,nmax=nmax,dist=dist,tol=.Machine$double.eps)
	param <- zero$root
	param
}

# function F for the PGGDS with normalization, i.e. F*muN/b, which leads to the same equilibrium points that F
F_pggds <- function(x,w,gamma,n,pn){
	F_pggds <- sum( pn * ( (1-x+w*x)^(n-1) - gamma*n ) )	
}

# function F for the PGGDS, without normalization (b and c are parameters)
F_pggds_b_c <- function(x,w,b,c,muN,n,pn){
	F_pggds_b_c <- (b/muN)* sum( pn * (1-x+w*x)^(n-1) ) - c	
}


# numerically finds xF for the PGGDS
xF_pggds <- function(w,gamma,muN,n,pn,dist){
	if(dist == "constant"){
		xF_pggds <- ( 1-(gamma*muN)^(1/(muN-1)) ) / ( 1-w )
	} else {
		zero <- 10*.Machine$double.eps
 		one <- 1-.Machine$double.eps
		sol <- uniroot(F_pggds,interval=c(zero,one),w=w,gamma=gamma,n=n,pn=pn,tol=.Machine$double.eps)
		xF_pggds <- sol$root
	}
}

# function F for the VD with normalization, i.e. F*muN/b, which leads to the same equilibrium points that F
F_vd <- function(x,gamma,n,pn){
	F_vd <- sum( pn * ( n*(1-x)^(n-1) - gamma*n ) )	
}

# function F for the PGGDS, without normalization (b and c are parameters)
F_vd_b_c <- function(x,b,c,muN,n,pn){
	F_vd_b_c <- (b/muN)* sum( pn * n*(1-x)^(n-1) ) - c	
}

# numerically finds xF for the VD
xF_vd <- function(gamma,muN,n,pn,dist){
	if(dist == "constant"){
		xF_vd <- 1-gamma^(1/(muN-1))
	} else {
		zero <- 10*.Machine$double.eps
#		zero <- .Machine$double.eps
 		one <- 1-.Machine$double.eps
		sol <- uniroot(F_vd,interval=c(zero,one),gamma=gamma,n=n,pn=pn,tol=.Machine$double.eps)
		xF_vd <- sol$root
	}
}

