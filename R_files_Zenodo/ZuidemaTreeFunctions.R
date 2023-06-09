## These functions include a 'ceiling' on growth to avoid eviction
## Search on 'Ceiling' to find where it is imposed. 
if(!exists("zmax")) zmax=10^6; 

###### Species names 
species.names<-c(
expression(italic("Annamocarya")),
expression(italic("Calocedrus")),
expression(italic("Dacrydium")),
expression(italic("Pinus")),
expression(italic("Manglietia")),
expression(italic("Parashorea"))) 

sppNames<-c("Annamocarya","Calocedrus","Dacrydium","Pinus","Manglietia","Parashorea");  

species.col<-c("black","dark grey","red","blue","green","purple")

######################################### 
# Survival function 
#########################################
sx <- function(dbh,pars) {
    dvals=c(pars[2],pars[3],pars[5],pars[6]); 
    svals=c(pars[1],pars[4],pars[4],pars[7]); 
    sfun=approxfun(dvals,svals,rule=2)
    return(sfun(dbh)) 
} 

#############################################################################
# Probability to reproduce, conditional on survival to t+1
#############################################################################
p.r <- function(dbh,pars) {
    u <- pars[9]+pars[10]*dbh;  
    pvals <- 1/(1+exp(-u)); 
    pvals[dbh<pars[8]] <- 0;  
    return(pvals); 
}    

###########################################################
# Expected number of seedlings, conditional on breeding 
###########################################################
bx <- function(dbh,pars) {
    return(p.r(dbh,pars)*pars[11])
}    

#########################################
# Size transition probability density 
#########################################
delta.dbh <- function(dbh,pars) {
    a <- pars[12]; b<- pars[13]; c<-pars[14];
    num <- b*c*dbh^(c-1);
    den <- (b+(dbh^c)/a)^2;
    return(num/den); 
}

# dbh=seq(0,130,length=1066); 
# plot(dbh,delta.dbh(dbh,param[[6]]),type="l"); 

g_x1x <- function(x1,x,pars) {
    xc <- pmin(x,zmax);     #################### Ceiling 
    x1bar <- xc + delta.dbh(xc,pars);
    return(dnorm(x1,mean=x1bar,sd=pars[15])); 
}

## CV of growth (size change) 
CVgro <- function(dbh,pars) {
    Edz = delta.dbh(dbh,pars) 
    return(pars[15]/Edz)
}    

#############################################
# Survival/growth kernel density function 
#############################################
p_x1x <- function(x1,x,pars) {
    return( sx(x,pars)*g_x1x(x1,x,pars) )  
}

#################################################
# Integrated kernel functions 
# Based on ipm_book/Rcode/c6/TreeKernel-intB2B.R 
################################################

### Integral of Growth 
intg_x1x <- function(x1,x,pars) {
    xc <- pmin(x,zmax);   #################### Ceiling 
    x1bar <- xc + delta.dbh(xc,pars);
    return(pnorm(x1,mean=x1bar,sd=pars[15])); 
}

### Integral of p(x1,x); 
intp_x1x <- function(x1,x,pars) {sx(x,pars)*intg_x1x(x1,x,pars)}

###########################################################################
#  Function to make IPM part of the P-matrix using integrated 
#  P(x',x) and Gaussian quadrature of specified order, over the initial bin.  
#  NOTE: P is modified so the largest size class can die but cannot
#        grow or shrink. This greatly reduces eviction, but is only OK 
#        if U is well above the maximum size that a tree reaches. 
##########################################################################

# Gauss-Legendre quadrature, nodes & weights for interval (L,U). 
# Uses gauss.quad from statmod package, and rescales as needed. 
gaussQuadInt <- function(L,U,order=5) {
	# nodes and weights on [-1,1]
	out <- gauss.quad(order); # G-L is the default 
    w <- out$weights; x <- out$nodes;  
    weights=0.5*(U-L)*w; 
    nodes=0.5*(U+L) + 0.5*(U-L)*x; 
	return(list(weights=weights,nodes=nodes)); 
}

## Note: this is a slow, loopy method. It's here just for debugging, to
## verify that the more vectorized version below gives the same results 
mk_P_int <- function(m,L,U,order,pars) {
		h <- (U - L)/m; meshpts <- L + ((1:m) - 1/2) * h
  		out=gaussQuadInt(-h/2,h/2,order); 
  		nodes=out$nodes; weights=out$weights; 
  		P <- matrix(0,m,m);
		for(i in 1:m){
			for(j in 1:m){
				  pvals1=intp_x1x(meshpts[i]-h/2,meshpts[j]+out$nodes,pars);	
			      pvals2=intp_x1x(meshpts[i]+h/2,meshpts[j]+out$nodes,pars);	
				  P[i,j]=sum(weights*(pvals2-pvals1))
			}
		}
		P<- P/h; # scale integral to get the average over the initial bin
 
        # modify P so that the largest size class can die or shrink, but not grow!  
        # NOTE: this is OK for these tree IPMs where growth is nearly deterministic
        # but a bad idea in general.  
        P[m,m] <- sx(meshpts[m],pars)-sum(P[1:(m-1),m]); 

		return(list(meshpts = meshpts, P = P))
}

## This function uses expand.grid to vectorize ALL P(z',z) evaluations for 
## integrated bin-to-bin construction of the P iteration matrix. 
## Impressed? Then thank Dylan Childs for this great idea.   
mk_P_int_noloop <- function(m, L, U, order,pars,fixEvict=TRUE) {
  h <- (U - L) / m 
  meshpts <- L + ((1:m) - 1/2) * h 
  out <- gaussQuadInt(-h/2, h/2, order) 
  quad <- expand.grid(x1=meshpts, x=meshpts, map=seq.int(order))
  quad <- transform(quad, nodes=out$nodes[map], weights=out$weights[map])
  P <- with(quad, {
    fvals1 <- intp_x1x(x1 - h/2, x + nodes,pars) 
    fvals2 <- intp_x1x(x1 + h/2, x + nodes,pars) 
    weights * (fvals2 - fvals1) 
  })
  dim(P) <- c(m, m, order)
  P <- apply(P, c(1,2), sum) / h
  
   # modify P so that the largest size class can die or shrink, but not grow!  
   # NOTE: this is OK for these tree IPMs where growth is nearly deterministic
   # but a bad idea in general.  
  P[m,m] <- sx(meshpts[m],pars)-sum(P[1:(m-1),m]); # to block eviction, if desired 
  
  # compensate for eviction out the bottom in the smallest sizes
    if(fixEvict){
     k=round(length(meshpts)/10,digits=0); k=pmax(k,1);  
     for(j in 1:k) { 
       sbar = 0.25*(sx(meshpts[j]-h/2,pars)+ 2*sx(meshpts[j],pars) + sx(meshpts[j]+h/2,pars))
       P[,j] = sbar*P[,j]/sum(P[,j]) 
     } 
  }	   
  
  return(list(meshpts = meshpts, P = P))
}

## For backwards compatibility, just in case. 
mk_P_int <- mk_P_int_noloop; 

######################################################################
#  Functions to make the full P-matrix, including seedlings (4 height  
#  classes and saplings/trees (IPM iteration matrix).  
######################################################################

### All new saplings go into the smallest size bin. 
### Loopy construction of P is used 
mk_full_P <- function(m,L,U,order,pars,A1) {
    out<- mk_P_int(m,L,U,order,pars)
    P2 <- rbind( matrix(0,4,m), out$P); 
    P1 <- rbind(A1[1:4,1:4], matrix(0,m,4));
    P <- cbind(P1,P2); 
    P[5,4] <- A1[5,4];
    return(list(P=P,meshpts=out$meshpts)); 
}  

### All new saplings go into the smallest size bin. 
### No-loop construction of P is used. 
mk_full_P_noloop <- function(m,L,U,order,pars,A1) {
    out<- mk_P_int_noloop(m,L,U,order,pars)
    P2 <- rbind( matrix(0,4,m), out$P); 
    P1 <- rbind(A1[1:4,1:4], matrix(0,m,4));
    P <- cbind(P1,P2); 
    P[5,4] <- A1[5,4];
    return(list(P=P,meshpts=out$meshpts)); 
} 

### The Zuidema et al. IPM as published: truncated Gaussian
### size distribution for new saplings. This causes them to "jump"
### above dbh=1, so that dbh is larger for a new sapling, than for
### a sapling that was a sapling last year with dbh=1.    
mk_full_P_noloop_saplingJump <- function(m,L,U,order,pars,A1) {
    out<- mk_P_int_noloop(m,L,U,order,pars)
    P2 <- rbind( matrix(0,4,m), out$P); 
    P1 <- rbind(A1[1:4,1:4], matrix(0,m,4));
    P <- cbind(P1,P2); 
    mu5 = pars[17]; sd5=pars[18]; 
    x=L+(0:m)*(U-L)/m;
    probs=diff(pnorm(x,mu5,sd5)); probs=probs/sum(probs);
    P[5:(m+4),4]=A1[5,4]*probs; 
    return(list(P=P,meshpts=out$meshpts,probs=probs)); 
} 
 
### Modified size distribution for new saplings. Instead of jumping,
### they are assigned the size distribution of a sapling of size x=L.
### As with other size classes, the distribution is truncated at L and
### rescaled to eliminate eviction out the bottom of the size range.     
mk_full_P_noloop_saplingSpread <- function(m,L,U,order,pars,A1) {
    out<- mk_P_int_noloop(m,L,U,order,pars)
    P2 <- rbind( matrix(0,4,m), out$P); 
    P1 <- rbind(A1[1:4,1:4], matrix(0,m,4));
    P <- cbind(P1,P2); 
    mu4 = L + delta.dbh(L,pars); sd4=pars[15]; 
    x=L+(0:m)*(U-L)/m;
    probs=diff(pnorm(x,mu4,sd4)); probs=probs/sum(probs);
    P[5:(m+4),4]=A1[5,4]*probs; 
    return(list(P=P,meshpts=out$meshpts,probs=probs)); 
} 

## set the default 
mk_full_P <- mk_full_P_noloop <- mk_full_P_noloop_saplingSpread  

#######################################################
# Fecundity distribution function conditional on survival. 
# Assumes Poisson distribution of seedlings, given the mean
# Inputs: parent size (dbh): ONE number, not a vector 
#         parameter vector (pars) 
#         range (M) 
#         normalize (default TRUE): scale so that sum=1? 
# Value:  vector of prob(#kids=0:M),length=M+1   
#######################################################
f_kx <- function(dbh,pars,M,normalize=TRUE){
    p1 <- p.r(dbh,pars); # prob of reproducing 
    lambda <- pars[11];  # mean #seedlings, if reproducing  
    pvals <- p1*dpois(0:M,lambda=pars[11]);
    pvals[1] <- pvals[1]+(1-p1); 
    if(normalize) pvals <- pvals/sum(pvals); 
    return(pvals)
}    

######################################################################
#  Function to make the B matrix, including seedlings (4 height  
#  classes) and saplings/trees (IPM iteration matrix).  
#  B[i,j] is the probability that a size-class-j individual has i-1 kids
#    conditional on survival (in this IPM, F includes S, so )
#  Inputs: P matrix, produced by mk_full_P
#		   meshpts, for the IPM part of P
# 		   pars = tree IPM parameters
#		   M = maximum possible #kids
#		   normalize (default=TRUE): scale so each column sums to 1?
#		   k = number of non-IPM height classes in P (default=4)
#  NOTE: for use with the functions to make the (size,kids) transition
#	 matrix, M should be set to at least the maximum lifetime total
#    reproductive output in the (size,kids) model. Even with M very
#    large, this is fast here (no nested loops and f_kx is loopless)	
######################################################################
mk_full_B <- function(P,meshpts,pars,M,normalize=TRUE,k=4) {
	m <- length(meshpts)
	B <- matrix(0,M+1,k+m)
	B[1,1:k]=1;
	for(j in 1:m) {
		probs <- f_kx(meshpts[j],pars,M,normalize=normalize)
		B[,k+j] <- probs;
	}
	return(B)
}

######################################################################
#  Make the B matrix, including seedlings and saplings,  
#  where T is the number of times that flowering occurs, NOT number of saplings.  
#  B[1,j] is the probability of not flowering, B[2,j] is the probability of flowering
#  Both of these are probabilities conditional on survival, which is what 
#  is needed for the p_xT function below that multiplies B*P to give
#  (size x kids) transitions
#
#  Inputs: P matrix, produced by mk_full_P
#		   meshpts, for the IPM part of P
# 		   pars = tree IPM parameters
#		   k = number of non-IPM height classes in P (default=4)
#  NOTE: for use with the functions to make the (size,kids) transition
#	 matrix, M should be set to at least the maximum lifetime total
#    reproductive output in the (size,kids) model. Even with M very
#    large, this is fast here: no nested loops, and f_kx is loopless
######################################################################
mk_small_B <- function(P,meshpts,pars,M,k=4) {
	m <- length(meshpts)
	B <- matrix(0,2,k+m)
	B[1,1:k]=1; B[1,(k+1):(k+m)] = 1-p.r(meshpts,pars); 
    B[2,]=1-B[1,]; 
    B <- rbind(B,matrix(0,M-2,ncol(B))); 
    return(B)
}
