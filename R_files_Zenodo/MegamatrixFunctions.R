## Convert state matrix to state vector
vec <- function(nmat) matrix(nmat,ncol=1) 

## Convert state vector to state matrix
unvec <- function(nvec,nrow=NULL,ncol=NULL){ 
	if(is.null(nrow)) return(matrix(nvec,ncol=ncol)); 
	if(is.null(ncol)) return(matrix(nvec,nrow=nrow)); 
}

## Convert iteration array to iteration matrix 
flatten <- function(A4) {
	dim(A4)<-rep(sqrt(length(A4)),2)
	return(A4)
}

## Convert iteration matrix to iteration array
unfold <- function(A2,dim) {
    dim(A2) <- dim; 
    return(A2)
} 

################################################
#  Utility functions for a size-quality model 
#  info is a list with parts mx,mq,hx,hq,yx,yq 
###############################################

####### Matrix multiplication without typing %*%  
mTimes <- function(a,b) {a%*%b} 
	
####### x (size) distribution function for whole population 
xDist=function(nt,info) {
    nmat=matrix(nt,info$mx,info$mq)
    return(apply(nmat,1,sum))
}

####### q (quality) distribution function for whole population 
qDist=function(nt,info) {
    nmat=matrix(nt,info$mx,info$mq)
    return(apply(nmat,2,sum))
}

###### Mean size as a function of q 
meanSize = function(nt,info) {
	nmat=matrix(nt,info$mx,info$mq)
	return( (matrix(info$yx,nrow=1)%*%nmat)/apply(nmat,2,sum) )
}

 

######################################################################
# Function to take B and P matrices, and compute the transition 
# probabilities from size-class i and j total kids, to all size classes
# and l total kids. This returns a vector of zeros if (l-j) is < 0
# or above the assumed maxnumber of kids per year,  ncol(B) - 1. 
######################################################################
p_xT <- function(l,i,j,B,P) {
	mx <- ncol(P); M <- nrow(B)-1;
	newKids <- (l-j); 
	if((newKids<0)|(newKids>M)) {
		return(rep(0,mx))
	}else{
		return(P[,i]*B[newKids+1,i])
	}
}

##############################################################################
# Function to make the 2D iteration matrix A for a size-kids model
# based on the P and B matrices summarizing a size-structured
# IPM. Apart from B and P the only input is mT, dimension for T 
# (so range of T is 0 to mT-1). The 4-D iteration array K is also returned. 
#
# Iteration matrix is modifited so individuals who get to the maximum
# values of T in the matrix stay there, but continue to grow/shrink/die
#############################################################################   
make_AxT <- function(B,P,mT) {
   mx=ncol(P); Kvals=array(0,c(mx,mT,mx,mT));  
   for(i in 1:mx){ # initial size
	  for(j in 1:mT){ # initial T 
	  	for(k in 1:mT){ # final T 
			Kvals[,k,i,j]=p_xT(k,i,j,B,P)
        }
    }
#	cat(i,"\\n"); 
  }
  # make kids-class mT absorbing: stay there with prob=1
  Kvals[1:mx,1:mT,1:mx,mT] <- 0; 
  Kvals[1:mx,mT,1:mx,mT] <- P; 
  A <- Kvals; dim(A) <- c(mx*mT,mx*mT); 
  return(list(A=A,K=Kvals)) 

}  


