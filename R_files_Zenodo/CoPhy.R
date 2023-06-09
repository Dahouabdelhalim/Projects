#Co-evolution of traits for two interacting lineages
  #DC Adams. December, 2017

#Function performs a phylogenetic regression for two traits (x & Y), 
#obtained from two co-evolving lineages (e.g., hosts and their 
#parasites). The user provides a formula describing the linear model 
#to be examined (e.g., y~x1+x2), where the dependent variables (Y) are 
#derived from one set of interacting species (e.g., parasites) while 
#the independent variables are derived from a second set interacting 
#species (e.g., hosts). Phylogenetic transformation of both X and Y is 
#performed and the model is statistically evaluated using multivariate 
#permutation PGLS procedures (see Adams 2014; Evolution; Adams and 
#Collyer 2015; Evolution).

#The method is capable of accommodating perfect 1:1 correspondence 
#between X and Y taxa, in which case each species in X is matched to  
#one (and only one) species in Y.  Additionally, patterns in data with 
#imperfect species correspondence may also be evaluated: where a 
#parasite inhabits more than one host; where more than one parasite is 
#found on a single host, or cases where some host or parasites are not 
#matched to a species in the other lineage (see below for details). 

#Parameters
# @param f1 A formula for the linear model (e.g., y ~ x1 + x2)

# @param dataX A list containing at least two elements: 1: the 
  #phylogeny (of class 'phylo') for the species represented by the 
  #independent (X) variables, and 2: one or more vectors containing 
  #the values for the independent (X) variables, whose element names 
  #correspond to the species names in the X-phylogeny.  The name of 
  #each vector must correspond to the name of one of the independent 
  #variables as listed in f1 (e.g., x1, x2, etc.). 

# @param dataY A list containing two elements: 1: the phylogeny (of 
  #class 'phylo') for the species represented by the dependent (Y) 
  #variables, and 2: a vector (or matrix) containing the dependent (Y) 
  #variables, whose names/rownames correspond to the species names in 
  #the Y-phylogeny. The name of this matrix must correspond to the 
  #name of the independent variable as listed in f1 (e.g., Y).

# @param matchlist A matrix describing which species for X correspond 
  #to species for Y, in X-Y order, i.e., the dependent (X) variable 
  #species are first, followed by the independent (Y) species. If more 
  #than one 'host' (X-species) contains more than one 'parasite' (Y-
  #species), that species name is repeated in the matchlist at the 
  #corresponding locations for its multiple associations.  For the 
  #converse scenario (where a parasites is found on more than one 
  #host, its name is repeated in the matchlist in the appropriate 
  #locations. Finally, for any species that is not matched with a   
  #species in the other lineage, 'NA' is used to designate its 'match' 
  #in the second lineage.

# @param iter The number of iterations for significance testing

##############
CoPhy.pgls<-function(f1,dataX=NULL,dataY=NULL, matchlist, iter=999){  
  library(ape)
  Pmat<-function(phy,x){
    C <- vcv.phylo(phy, anc.nodes = FALSE)
    C <- C[rownames(x), rownames(x)]
    eigC <- eigen(C)
    lambda <- zapsmall(eigC$values)
    if (any(lambda == 0)) {   lambda = lambda[lambda > 0]  }
    eigC.vect = eigC$vectors[, 1:(length(lambda))]
    Pmat<-eigC.vect %*% diag(sqrt(lambda)) %*% t(eigC.vect)
    Px<-if (det(Pmat) > 1e-08) qr.solve(Pmat) else fast.ginv(Pmat)
    rownames(Px) <- colnames(Px) <- colnames(C)
    Px
  }
  fast.ginv<-function (X, tol = sqrt(.Machine$double.eps)) {
    k <- ncol(X)
    Xsvd <- La.svd(X, k, k)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    rtu <- ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    v <- t(Xsvd$vt)[, Positive, drop = FALSE]
    v %*% rtu
  }
  pval<-function(s){
    p = length(s)
    r = rank(s)[1] - 1
    pv = 1 - r/p
    pv
  }
  if(is.null(dataX))stop("No dataX provided.")
  if(is.null(dataY))stop("No dataY provided.") 
  NA.vals<-length(which(is.na(matchlist)==TRUE))
  if(NA.vals>0){
    NA.rows<-unlist(apply(matchlist,2,function (x) which(is.na(x))))
    matchlist<-matchlist[-NA.rows,]  
  }
  for(i in 1:length(dataX)) assign(names(dataX)[i], dataX[[i]])  
  checkX <- sapply(dataX, class)
  if(is.na(match("phylo",checkX)))stop("No phylogeny of class Phylo in list 'dataX'.") else{
    phyX<-dataX[[match("phylo",checkX)]]    
  }  
  y.tmp<-rep(1,length(phyX$tip.label))
  tmp.f<-as.formula(paste("y.tmp", formula(f1)[3],sep="~"))  
  x.mod<-model.matrix(terms(formula(tmp.f)))
  Terms<-terms(tmp.f)
  term.labels <- attr(Terms, "term.labels")
  X.k <- attr(x.mod, "assign")
  QRx <- qr(x.mod)
  x.mod <- x.mod[, QRx$pivot, drop = FALSE]
  x.mod <- x.mod[, 1:QRx$rank, drop = FALSE]
  X.k <- X.k[QRx$pivot][1:QRx$rank]
  uk <- unique(c(0,X.k))
  k <- length(attr(Terms, "term.labels"))
  Xs <- lapply(1:length(uk), function(j)  Xj <- x.mod[, X.k %in% uk[1:j]])
  Xrs <- Xs[1:k]
  Xfs <- Xs[2:(k+1)]
  Px<-Pmat(phy = phyX,x = as.matrix(dataX[[which(checkX !="phylo")[1]]]))
  nX<-nrow(Px)
  Xr <- lapply(Xrs, function(x) crossprod(Px, as.matrix(x))) 
  Xf <- lapply(Xfs, function(x) crossprod(Px, as.matrix(x))) 
  Xr.m<-lapply(Xr, function(x) as.matrix(x[match(matchlist[,1],rownames(x)),])) 
  Xf.m<-lapply(Xf, function(x) as.matrix(x[match(matchlist[,1],rownames(x)),]))  
  for(i in 1:length(dataY)) assign(names(dataY)[i], dataY[[i]])  
  checkY <- sapply(dataY, class)
  if(is.na(match("phylo",checkY)))stop("No phylogeny of class Phylo in list 'dataY'.") else{
    phyY<-dataY[[match("phylo",checkY)]]    
  }   
  Py<-Pmat(x = as.matrix(dataY[[which(checkY !="phylo")[1]]]),phy = phyY)
  Y<-as.matrix(get(names(dataY)[which(checkY!="phylo")]))  
  ind <- c(list(1:nrow(Y)), (Map(function(x) sample.int(nrow(Y), nrow(Y)), 1:iter)))
  perms <- length(ind)
  Y.pr<-crossprod(Py,Y)  
  Y.pr.m<-as.matrix(Y.pr[match(matchlist[,2],rownames(Y.pr)),])     
  N.xy<-n<-nrow(Y.pr.m)
  p<-ncol(Y.pr.m[[1]])
  gls.init<-lapply(1:length(Xr.m), function(j) lm.fit(Xr.m[[j]],Y.pr.m) ) 
  fitted<-lapply(1:length(Xr.m), function(j) as.matrix(gls.init[[j]]$fitted.values))
  res<-lapply(1:length(Xr.m), function(j) as.matrix(gls.init[[j]]$residuals))
  SS.p <- lapply(1: perms, function(j){ 
    Yi<-Map(function(f, r) f + r[ind[[j]], ], fitted, res)
    c(Map(function(y, ur, uf) sum(.lm.fit(ur,y)$residuals^2  - .lm.fit(uf,y)$residuals^2 ),
          Yi, Xr.m, Xf.m),
      sum(.lm.fit(Xr.m[[k]],Yi[[k]])$residuals^2)- sum(.lm.fit(Xr.m[[k]],Yi[[k]])$residuals^2  - .lm.fit(Xf.m[[k]],Yi[[k]])$residuals^2 ),
      sum(.lm.fit(Xr.m[[1]],Yi[[1]])$residuals^2) )
  })
  SS.p <- matrix(unlist(SS.p), k+2, perms)
  SS <- SS.p[1:k,]
  SSE <- SS.p[k+1,] 
  SSY <- SS.p[k+2,]
  df<-sapply(1:k, function(j) qr(Xf.m[[j]])$rank - qr(Xr.m[[j]])$rank)
  dfE <- n- sum(df) -1
  df <- c(df, dfE, n-1)  
  MS <- SS/df[1:k]
  MSE <- SSE/df[k+1]
  SSE.mat <- matrix(SSE, k, length(SSE), byrow = TRUE)
  MSE.mat <- matrix(MSE, k, length(MSE), byrow = TRUE)
  if(is.matrix(SS)){
    Fs <- (SS/df[1:k])/MSE.mat
    P.val <- apply(Fs, 1, pval)
  } else {
    MSE <- SSE/df[2]
    Fs <- (SS/df[1])/MSE
    P.val <- pval(Fs)
  }
  SS.obs<-SS.p[,1]
  if(k==1){MS.obs <- c(MS[1],MSE[[1]],NA)} else MS.obs <- c(MS[,1],MSE[[1]],NA)
  R2.obs <- c((SS.obs/SSY[[1]])[1:k],NA,NA)
  if(k==1){Fs <- c(Fs[1],NA,NA)} else Fs <- c(Fs[,1],NA,NA)
  Pr<-c(P.val,NA,NA)
  anova.tab <- data.frame(df,SS=SS.obs,MS=MS.obs,Rsq=R2.obs,F=Fs, Pr=Pr)
  rownames(anova.tab) <- c(term.labels, "Residuals", "Total")
  b <- lm(Y.pr.m~Xf.m[[k]]-1)$coefficients
  out <- list(anova.table = anova.tab, reg.coef = b)
  out
}

