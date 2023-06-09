## Functions

## calcurate SE
se  <-  function(x){
  y  <-  x[!is.na(x)]  #  remove  the  missing  values
  sqrt(var(as.vector(y))/length(y))
}

## Create a dataframe with mean, sd, se
data_summrize <- function(d, variable, cname){
  l <- length(cname)
  if(l>=2){
    df <- data.frame(as.vector(tapply(d[,cname][,1], d[,cname], unique)))
    for(i in 2:l){
      df <- data.frame(df,
                       as.vector(tapply(d[,cname][,i], d[,cname], unique))
      )
    }
  } else {
    df <- data.frame(as.vector(tapply(d[,cname], d[,cname], unique)))
  }
  df <- data.frame(df,
                   as.vector(tapply(d[,variable], d[,cname], mean)),
                   as.vector(tapply(d[,variable], d[,cname], sd)),
                   as.vector(tapply(d[,variable], d[,cname], se)))
  colnames(df) <- c(cname, paste0(variable, ".", c("mean", "sd", "se")))
  return(df)
}

## Phylogenetically corrected mean
phyMean<- function(x=trait.value, VV=cvc)	{
  vec1<- rep(1, length(x))
  invVV<- solve(VV)
  aa<- solve(t(vec1)%*%invVV %*% vec1) %*%t (vec1) %*% invVV %*% x
  aat<- t(aa)
  return(aat)
}

##### Log likelihood function #####
## Unbiased Brownian Motion
logL.uBM<- function(p, x, VV, meserr) # parameter set, data, Variance-Covariance matrix, sampling error
{
  # p[1] = log(variance), p[2] = ancestral state
  VVt<- exp(p[1])*VV + diag(meserr^2)  # add sampling error to VCV matrix
  logl<- dmnorm(x, mean=p[2], varcov=VVt, log=TRUE)
  return(logl)	
}

## OU
logL.OU <- function(p, x, meserr, nt, VV) {
  # p[1] = log(variance), p[2] = alpha, p[3] = ancestral state
  ou.mat <- ouMatrix(VV, exp(p[2]), nt)
  VVt <- exp(p[1]) * ou.mat
  diag(VVt) <- diag(VVt) + meserr^2
  logl <- dmnorm(x, p[3], VVt, log = T)
  return(logl)
}
ouMatrix <- function (vcvMatrix, alpha, nt) 
{
  vcvDiag <- diag(vcvMatrix)
  diagi <- matrix(vcvDiag, nrow = nt, ncol = nt)
  diagj <- matrix(vcvDiag, nrow = nt, ncol = nt, byrow = T)
  Tij = diagi + diagj - (2 * vcvMatrix) # The time separating the species i and j (i.e., phylogenetic distance)
  vcvRescaled = (1/(2 * alpha)) * exp(-alpha * Tij) * (1 -  exp(-2 * alpha * vcvMatrix))
  return(vcvRescaled)
}

## Biased Brownian Motion
logL.trend<- function(p, x, VV, meserr) # parameter set, data, Variance-Covariance matrix, sampling error
{
  # p[1] = log(variance), p[2] = mean, p[3] = ancestral state
  VVt<- exp(p[1])*VV + diag(meserr^2) 
  MM <- p[3] + p[2] * diag(VV)
  logl<- dmnorm(x, mean=MM, varcov=VVt, log=TRUE)
  return(logl)	
}

## White noise
lnl.noise <- function(p, x, VV, meserr) {
  # p[1] = variance, p[2] = ancestral state
  n <- length(x)
  VV <- diag(exp(p[1]), nrow = n)
  diag(VV) <- diag(VV) + meserr^2
  logl<- dmnorm(x, p[2], VV, log = TRUE)
  return(logl)
}

## Pagel's lambda
lnl.lambda <- function(p, x, VV, meserr) {
  # p[1] = variance, p[2] = lambda, p[3] = ancestral state
  n <- length(x)
  index = matrix(TRUE, n, n)
  diag(index) = FALSE
  VVT <- VV
  VVT[index] <- VVT[index] * (p[2])
  VVT <- VVT * exp(p[1])
  diag(VVT) <- diag(VVT) + meserr^2
  logl<- dmnorm(x, p[3], VVT, log = TRUE)
  return(logl)
}

## Pagel's kappa
lnl.kappa <- function(p, x, VV, meserr, tree) {
  # p[1] = variance, p[2] = lambda, p[3] = ancestral state
  t <- rescale(tree, "kappa", kappa  = exp(p[2]))
  VVt <- exp(p[1]) * vcv(t) + diag(meserr^2)
  logl <- dmnorm(x, mean=p[3], varcov=VVt, log = T)
  return(logl)
}

## Pagel's delta
lnl.delta <- function(p, x, VV, meserr, tree) {
  # p[1] = variance, p[2] = lambda, p[3] = ancestral state
  t <- rescale(tree, "delta", delta  = exp(p[2]))
  VVt <- exp(p[1]) * vcv(t) + diag(meserr^2)
  logl <- dmnorm(x, mean=p[3], varcov=VVt, log = T)
  return(logl)
}

## Shift from biased BM to unbiased BM (separate variance)
logL.trendShift.vs<- function(p, x, mmi, CC, meserr, bn)
{
  nreg<- length(bn)+1  # number of regimes = no. of shift points + 1
  ms<- c(p[1],0)
  vs<- exp(p[2:3])
  x0<- p[4]
  
  # get vector of means for tips from this, then logL
  MM<- x0 + colSums(ms*mmi)
  VVt<- vs[1]*CC[[1]]
  VVt<- VVt + vs[2]*CC[[2]]
  VVt<- VVt + diag(meserr^2)  # add sampling error to VCV matrix
  
  logl<- dmnorm(x, mean=MM, varcov=VVt, log=TRUE)
  return(logl)	
}

