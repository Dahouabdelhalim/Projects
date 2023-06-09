
# Online Resource 1 - R code for calculating all measures of pairwise and multiple site 
                # beta diveristy presented in the manuscript.
				# Included are scripts to generate the scenarios in Fig. 1 and for the 
				# simulation used to generate Fig. 2. The code presented here builds on 
				# previously published work, including Baselga (2010, 2012, 2013), Baselga
				# et al. (2013) and Carvalho et al. (2013)


# Supplementary information to: Ensing, D.J.* and J. Pither. (XXXX). A novel multiple-site 
	    # extension to pairwise partitioned taxonomic beta diversity. Ecological Complexity.


# *corresponding author

# email: djensing@gmail.com
# Tel: +1.613.533.6000

# Present Address:
# Department of Biology
# QueenÕs University
# Biosciences Complex, 116 Barrie Street
# Kingston, ON, Canada
# K7L 3N6


# this work requires Baselga et al.'s (2013) 'betapart' package:

# load required packages
library(betapart); 


##########################################################################################
# A list of the tools used and developed in this work, as presented below
#
#'carv.beta' - a modification of the functions provided by Carvalho et al. (2013) to 
#				produce a list output of pairwise beta diversity and its components under 
#				their parition, similar to the output of the beta.pair function of Baselga 
#				et al. (2013).
###
#
#'beta.multi.carv' - a modification of the beta.multi function from the betapart package 
#				(Baselga et al., 2013) which includes "carvalho" as an index family, it is
#				 the default. "jaccard" or "sorensen" may still be specified
##########################################################################################

# first a modification of Carvalho et al.'s (2013) taxonomic beta partition to give a list
# output as Baselga et al.'s (2013) beta.pair function does

carv.beta <- function(site.spec){# where site.spec is the site by species matrix you want 
								 # to calculate beta diversity on
beta.cc<-function (x)
{
x <- ifelse(x > 0, 1, 0)
d <- tcrossprod(x)
a <- as.dist(d)
S <- diag(d)
N <- length(S)
b <- as.dist(matrix(rep(S, N), nrow = N)) - a
c <- as.dist(matrix(rep(S, each = N), nrow = N)) - a
out = (b+c)/(a+b+c)
out
}
beta.3<-function (x)
{
x <- ifelse(x > 0, 1, 0)
d <- tcrossprod(x)
a <- as.dist(d)
S <- diag(d)
N <- length(S)
b <- as.dist(matrix(rep(S, N), nrow = N)) - a
c <- as.dist(matrix(rep(S, each = N), nrow = N)) - a
out = 2*pmin(b,c)/(a+b+c)
out
}
beta.rich<-function (x)
{
x <- ifelse(x > 0, 1, 0)
d <- tcrossprod(x)
a <- as.dist(d)
S <- diag(d)
N <- length(S)
b <- as.dist(matrix(rep(S, N), nrow = N)) - a
c <- as.dist(matrix(rep(S, each = N), nrow = N)) - a
out = abs(b-c)/(a+b+c)
out
}

#generate a list of the distance matrices of each diversity measure
all.carv.beta <-list(beta.3=beta.3(site.spec),
					betarich=beta.rich(site.spec),
					beta.cc=beta.cc(site.spec))
return(all.carv.beta)#return the list of distance matrices
}

#################################################


# Load function for multiple site measures of beta diversity

## Note that this function is capable of calculatig all three partitions (i.e. it includes 
## Baselga's multiple-site partition), but these must be specified using the index.family 
## argument of Baselga's beta.multi (i.e. 'carvalho' is default). This is directly 
## modified from the betapart package of Baselga et al. (2013), see also: 
## http://127.0.0.1:20382/library/betapart/html/betapart-package.html

beta.multi.carv <- function (x, index.family = "carvalho") 
{library(betapart)
    index.family <- match.arg(index.family, c("jaccard", "sorensen","carvalho"))
    if (!inherits(x, "betapart")) {
        x <- betapart.core(x)
    }
    maxbibj <- sum(x$max.not.shared[lower.tri(x$max.not.shared)])
    minbibj <- sum(x$min.not.shared[lower.tri(x$min.not.shared)])
    switch(index.family, sorensen = {
        beta.sim <- minbibj/(minbibj + x$a)
        beta.sne <- (x$a/(minbibj + x$a)) * ((maxbibj - minbibj)/((2 * 
            x$a) + maxbibj + minbibj))
        beta.sor <- (minbibj + maxbibj)/(minbibj + maxbibj + 
            (2 * x$a))
        multi <- list(beta.SIM = beta.sim, beta.SNE = beta.sne, 
            beta.SOR = beta.sor)
    }, jaccard = {
        beta.jtu <- (2 * minbibj)/((2 * minbibj) + x$a)
        beta.jne <- (x$a/((2 * minbibj) + x$a)) * ((maxbibj - 
            minbibj)/((x$a) + maxbibj + minbibj))
        beta.jac <- (minbibj + maxbibj)/(minbibj + maxbibj + 
            x$a)
        multi <- list(beta.JTU = beta.jtu, beta.JNE = beta.jne, 
            beta.JAC = beta.jac)
    }, carvalho = {
        beta.3 <- (2 * minbibj)/(minbibj + maxbibj + x$a)
        beta.rich <- ((maxbibj - minbibj)/(minbibj + maxbibj + x$a))
        beta.cc <- (minbibj + maxbibj)/(minbibj + maxbibj + x$a)
        multi <- list(beta.3multi = beta.3, beta.RICH = beta.rich, 
            beta.CC = beta.cc)
    })
    return(multi)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Appendix to: 
# Legendre, P. (2014) Interpreting the replacement and richness difference   
# components of beta diversity. Global Ecology and Biogeography, 23, xxx–xxx.  
#                                  Appendix S5 
# 
# R function to compute LCBD indices from a dissimilarity matrix (D) or from  
# beta diversity component matrices (Repl, RichDiff/AbDiff or Nes).  
LCBD.comp <- function(x, sqrt.x=TRUE) 
  # 
  # Description -- 
  # 
  # Computes LCBD indices (Legendre and De Cáceres 2013) from a dissimilarity  
  # matrix (D) or beta div. component matrices (Repl, RichDiff/AbDiff or Nes). 
  # 
  # Arguments -- 
  # 
  ###IMPORTANT NOTE: INPUT MUST BE A DISTANCE MATRIX!!!
  # x : D or beta diversity component matrix, class=dist. 
  # sqrt.x : Take sqrt() of components before computing LCBD.comp. Use 
#     sqrt.x=TRUE for the replacement and richness/abundance difference indices  
#     computed by beta.div.comp(), as well as for the corresponding D matrices. 
# 
# Reference -- 
# 
# Legendre, P. & De Cáceres, M. (2013) Beta diversity as the variance of  
# community data: dissimilarity coefficients and partitioning. Ecology  
# Letters 16: 951–963.  
# 
# License: GPL-2  
# Author:: Pierre Legendre, August 2013 
{ 
  ### Internal function 
  centre <- function(D,n) 
    # Centre a square matrix D by matrix algebra 
    # mat.cen = (I - 11'/n) D (I - 11'/n) 
  {  One <- matrix(1,n,n) 
  mat <- diag(n) - One/n 
  mat.cen <- mat %*% D %*% mat 
  } 
  ### 
  n <- nrow(as.matrix(x)) 
  
  if(sqrt.x) { 
    # x = sqrt(x) 
    SStotal <- sum(x)/n        # eq. 8 
    BDtotal <- SStotal/(n-1)   # eq. 3 
    G <- centre(as.matrix(-0.5*x), n)     # Gower-centred matrix 
  } else { 
    SStotal <- sum(x^2)/n      # eq. 8 
    BDtotal <- SStotal/(n-1)   # eq. 3 
    G <- centre(as.matrix(-0.5*x^2), n)   # Gower-centred matrix 
  } 
  LCBD <- diag(G)/SStotal   # Legendre & De Caceres (2013), eq. 10b out <- list(SStotal_BDtotal=c(SStotal,BDtotal), LCBD=LCBD, D=x) 
}


# Appendix S3 from Legendre et al., 2014
#
# R function to compute the Podani- or Baselga-family decomposition of the 
# Jaccard or SÃ¸rensen groups into replacement and richness difference 
# (or nestedness) components, for species presence-absence or abundance data.
beta.div <- function(Y, method="hellinger", sqrt.D=FALSE, samp=TRUE, nperm=999, save.D=FALSE, clock=FALSE)
  #
  # Compute estimates of total beta diversity as the total variance in Y, 
  # for 20 dissimilarity coefficients or analysis of raw data (not recommended). 
  # LCBD indices are tested by permutation within columns of Y.
  # This version includes direct calculation of the Jaccard, Sorensen and Ochiai 
  # coefficients for presence-absence data.
  #
  # Arguments --
  # 
  # Y : community composition data matrix.
  # method : name of one of the 20 dissimilarity coefficients, or "none" for
#          direct calculation on Y (also the case with method="euclidean").
# sqrt.D : If sqrt.D=TRUE, the distances in matrix D are square-rooted before 
#          computation of SStotal, BDtotal and LCBD. 
# samp : If samp=TRUE, the abundance-based distances (ab.jaccard, ab.sorensen,
#        ab.ochiai, ab.simpson) are computed for sample data. If samp=FALSE, 
#        they are computed for true population data.
# nperm : Number of permutations for test of LCBD.
# save.D : If save.D=TRUE, the distance matrix will appear in the output list.
# clock : If clock=TRUE, the computation time is printed in the R console.
#
# Reference --
#
# Legendre, P. and M. De CÃ¡ceres. 2013. Beta diversity as the variance of 
# community data: dissimilarity coefficients and partitioning. 
# Ecology Letters 16: 951-963. 
#
# License: GPL-2 
# Author:: Pierre Legendre, December 2012, April-May 2013
{
  ### Internal functions
  centre <- function(D,n)
    # Centre a square matrix D by matrix algebra
    # mat.cen = (I - 11'/n) D (I - 11'/n)
  {	One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  mat.cen <- mat %*% D %*% mat
  }
  ###
  BD.group1 <- function(Y, method, save.D, per, n)
  {
    if(method=="profiles") Y = decostand(Y, "total")
    if(method=="hellinger") Y = decostand(Y, "hellinger")
    if(method=="chord") Y = decostand(Y, "norm")
    if(method=="chisquare") Y = decostand(Y, "chi.square")
    #
    s <- scale(Y, center=TRUE, scale=FALSE)^2   # eq. 1
    SStotal <- sum(s)          # eq. 2
    BDtotal <- SStotal/(n-1)   # eq. 3
    if(!per) { SCBD<-apply(s,2,sum)/SStotal }else{ SCBD<-NA }  # eqs. 4a and 4b
    LCBD <- apply(s, 1, sum)/SStotal  # eqs. 5a and 5b
    #
    D <- NA
    if(!per & save.D)   D <- dist(Y)
    #
    out <- list(SStotal_BDtotal=c(SStotal,BDtotal), SCBD=SCBD, LCBD=LCBD, 
                method=method, D=D)
  }
  ###
  BD.group2 <- function(Y, method, sqrt.D, n)
  {
    if(method == "divergence") {
      D = D11(Y)		
      
    } else if(any(method == 
                  c("jaccard","sorensen","ochiai"))) 
    {
      if(method=="jaccard") D = dist.binary(Y, method=1) # ade4 takes sqrt(D)
      if(method=="sorensen")  D = dist.binary(Y, method=5) #ade4 takes sqrt(D)
      if(method=="ochiai") D = dist.binary(Y, method=7) # ade4 takes sqrt(D)
      
    } else if(any(method == 
                  c("manhattan","canberra","whittaker","percentagedifference","wishart"))) 
    {
      if(method=="manhattan") D = vegdist(Y, "manhattan")
      if(method=="canberra")  D = vegdist(Y, "canberra")
      if(method=="whittaker") D = vegdist(decostand(Y,"total"),"manhattan")/2
      if(method=="percentagedifference") D = vegdist(Y, "bray")
      if(method=="wishart")   D = WishartD(Y)
    } else {
      if(method=="modmeanchardiff") D = D19(Y)
      if(method=="kulczynski")  D = vegdist(Y, "kulczynski")
      if(method=="ab.jaccard")  D = chao(Y, coeff="Jaccard", samp=samp)
      if(method=="ab.sorensen") D = chao(Y, coeff="Sorensen", samp=samp)
      if(method=="ab.ochiai")   D = chao(Y, coeff="Ochiai", samp=samp)
      if(method=="ab.simpson")  D = chao(Y, coeff="Simpson", samp=samp)
      #ST: additional method added for Horn's information overlap index 
      if(method=="horn"){  
        dfHORN = sim.table (Y, q = 1)
        colnames(dfHORN) = rownames(Y)
        rownames(dfHORN) = rownames(Y)			
        D =  (1-as.dist(dfHORN))}
    }
    #	
    if(sqrt.D) D = sqrt(D)
    SStotal <- sum(D^2)/n      # eq. 8
    BDtotal <- SStotal/(n-1)   # eq. 3
    delta1 <- centre(as.matrix(-0.5*D^2), n)   # eq. 9
    LCBD <- diag(delta1)/SStotal               # eq. 10b
    #
    out <- list(SStotal_BDtotal=c(SStotal,BDtotal), LCBD=LCBD, 
                method=method, D=D)
  }
  ###
  ###
  epsilon <- sqrt(.Machine$double.eps)
  method <- match.arg(method, c("euclidean", "manhattan", "modmeanchardiff", "profiles", "hellinger", "chord", "chisquare", "divergence", "canberra", "whittaker", "percentagedifference", "wishart", "kulczynski", "ab.jaccard", "ab.sorensen","ab.ochiai","ab.simpson","jaccard","sorensen","ochiai","none","horn"))
  #
  if(any(method == c("profiles", "hellinger", "chord", "chisquare", "manhattan", "modmeanchardiff", "divergence", "canberra", "whittaker", "percentagedifference", "kulczynski"))) require(vegan)
  if(any(method == c("jaccard","sorensen","ochiai"))) require(ade4)
  if(any(method == c("horn"))){ 
    require(vegetarian)  
    require(ade4)}
  #
  if(is.table(Y)) Y <- Y[1:nrow(Y),1:ncol(Y)]    # In case class(Y) is "table"
  n <- nrow(Y)
  if((n==2)&(dist(Y)[1]<epsilon)) stop("Y contains two identical rows, hence BDtotal = 0")
  #
  aa <- system.time({
    if(any(method == 
           c("euclidean", "profiles", "hellinger", "chord", "chisquare","none"))) {
      note <- "Info -- This coefficient is Euclidean"
      res <- BD.group1(Y, method, save.D, per=FALSE, n)
      #
      # Permutation test for LCBD indices, distances group 1
      if(nperm>0) {
        p <- ncol(Y)
        nGE.L = rep(1,n)
        for(iperm in 1:nperm) {
          Y.perm = apply(Y,2,sample)
          res.p <- BD.group1(Y.perm, method, save.D, per=TRUE, n)
          ge <- which(res.p$LCBD+epsilon >= res$LCBD)
          nGE.L[ge] <- nGE.L[ge] + 1
        }
        p.LCBD <- nGE.L/(nperm+1)
      } else { p.LCBD <- NA }
      #
      if(save.D) { D <- res$D } else { D <- NA }
      #
      out <- list(SStotal_BDtotal=res$SStotal_BDtotal, SCBD=res$SCBD, 
                  LCBD=res$LCBD, p.LCBD=p.LCBD, method=method, note=note, D=D)
      
    } else {
      #
      res <- BD.group2(Y, method, sqrt.D, n)
      #
      # Permutation test for LCBD indices, distances group 2
      if(nperm>0) {
        nGE.L = rep(1,n)
        for(iperm in 1:nperm) {
          Y.perm = apply(Y,2,sample)
          res.p <- BD.group2(Y.perm, method, sqrt.D, n)
          ge <- which(res.p$LCBD+epsilon >= res$LCBD)
          nGE.L[ge] <- nGE.L[ge] + 1
        }
        p.LCBD <- nGE.L/(nperm+1)
      } else { p.LCBD <- NA }
      
      #		
      if(method == "divergence") {
        note = "Info -- This coefficient is Euclidean"
      } else if(any(method == c("jaccard","sorensen","ochiai"))) {
        note = c("Info -- This coefficient is Euclidean because dist.binary ",
                 "of ade4 computes it as sqrt(D). Use beta.div with option sqrt.D=FALSE")
      } else if(any(method == 
                    c("manhattan","canberra","whittaker","percentagedifference","wishart"))) {
        if(sqrt.D) {
          note = "Info -- This coefficient, in the form sqrt(D), is Euclidean"
        } else {
          note = c("Info -- For this coefficient, sqrt(D) would be Euclidean", 
                   "Use is.euclid(D) of ade4 to check Euclideanarity of this D matrix")
        }
      } else if(any(method == c("horn"))) {
        if(is.euclid(sqrt(res$D))) {			
          note = "Info -- This coefficient, in the form sqrt(D), is Euclidean"
        } else if(is.euclid(res$D)){
          note = "Info -- This coefficient, in untransformed form, is Euclidean"
        }
        else {
          note = "Info -- This coefficient is not Euclidean"
        }
      } else {
        note = c("Info -- This coefficient is not Euclidean", 
                 "Use is.euclid(D) of ade4 to check Euclideanarity of this D matrix")
      }
      
      #
      if(sqrt.D) note.sqrt.D<-"sqrt.D=TRUE"  else  note.sqrt.D<-"sqrt.D=FALSE"
      if(save.D) { D <- res$D } else { D <- NA }
      #
      out <- list(SStotal_BDtotal=res$SStotal_BDtotal, LCBD=res$LCBD,  
                  p.LCBD=p.LCBD, method=c(method,note.sqrt.D), note=note, D=D)
    }
    #
  })
  aa[3] <- sprintf("%2f",aa[3])
  if(clock) cat("Time for computation =",aa[3]," sec\\n")
  #
  class(out) <- "beta.div"
  out
}

D11 <- function(Y, algo=1)
  #
  # Compute Clark's coefficient of divergence. 
  # Coefficient D11 in Legendre and Legendre (2012, eq. 7.51).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2011
{
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  # Prepare to divide by pp = (p-d) = no. species present at both sites
  Y.ap <- 1 - decostand(Y, "pa")
  d <- Y.ap %*% t(Y.ap)
  pp <- p-d   # n. species present at the two compared sites
  #
  if(algo==1) {   # Faster algorithm
    D <- matrix(0, n, n)
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        num <- (Y[i,]-Y[j,])
        den <- (Y[i,]+Y[j,])
        sel <- which(den > 0)
        D[i,j] = sqrt(sum((num[sel]/den[sel])^2)/pp[i,j])
      }
    }
    #
  } else {   # Slower algorithm 
    D <- matrix(0, n, n)
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        temp = 0
        for(p2 in 1:p) {
          den = Y[i,p2] + Y[j,p2]
          if(den > 0) {
            temp = temp + ((Y[i,p2] - Y[j,p2])/den)^2
          }
        }
        D[i,j] = sqrt(temp/pp[i,j])
      }
    }
    #
  }	
  DD <- as.dist(D)
}

D19 <- function(Y)
  #
  # Compute the Modified mean character difference.
  # Coefficient D19 in Legendre and Legendre (2012, eq. 7.46).
  # Division is by pp = number of species present at the two compared sites
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, April 2011
{
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  # Prepare to divide by pp = (p-d) = n. species present at both sites
  Y.ap <- 1 - decostand(Y, "pa")
  d <- Y.ap %*% t(Y.ap)
  pp <- p-d   # n. species present at the two compared sites
  #
  D <- vegdist(Y, "manhattan")
  DD <- as.dist(as.matrix(D)/pp)
}

WishartD <- function(Y)
  #
  # Compute dissimilarity - 1 - Wishart similarity ratio (Wishart 1969).
  #
  # License: GPL-2 
  # Author:: Pierre Legendre, August 2012
{
  CP = crossprod(t(Y))
  SS = apply(Y^2,1,sum)
  n = nrow(Y)
  mat.sq = matrix(0, n, n)
  for(i in 2:n) {
    for(j in 1:(n-1)) { mat.sq[i,j] = CP[i,j]/(SS[i] + SS[j] - CP[i,j]) }
  }
  mat = 1 - as.dist(mat.sq)
}

chao <- function(mat, coeff="Jaccard", samp=TRUE)
  #
  # Compute Chao et al. (2006) abundance-based indices.
  #
  # Arguments -
  # mat = data matrix, species abundances
  # coef = "Jaccard" : modified abundance-based Jaccard index
  #        "Sorensen": modified abundance-based SÃ¸rensen index
  #        "Ochiai"  : modified abundance-based Ochiai index
  #        "Simpson" : modified abundance-based Simpson index
  # samp=TRUE : Compute dissimilarities for sample data
  #     =FALSE: Compute dissimilarities for true population data
#
# Details -
# For coeff="Jaccard", the output values are identical to those
# produced by vegan's function vegdist(mat, "chao").
#
# Help received from A. Chao and T. C. Hsieh in July 2012 for the computation  
# of dissimilarities for true population data is gratefully acknowledged.
#
# Reference --
# Chao, A., R. L. Chazdon, R. K. Colwell and T. J. Shen. 2006. 
# Abundance-based similarity indices and their estimation when there 
# are unseen species in samples. Biometrics 62: 361â€“371.
#
# License: GPL-2 
# Author:: Pierre Legendre, July 2012
{
  require(vegan)
  nn = nrow(mat)
  res = matrix(0,nn,nn)
  if(samp) {   # First for sample data
    for(k in 2:nn) {
      for(j in 1:(k-1)) {
        #cat("k =",k,"  j =",j,"\\n")
        v1 = mat[j,]   # Vector 1
        v2 = mat[k,]   # Vector 2
        v1.pa = decostand(v1,"pa")   # Vector 1 in presence-absence form
        v2.pa = decostand(v2,"pa")   # Vector 2 in presence-absence form
        N.j = sum(v1)   # Sum of abundances in vector 1
        N.k = sum(v2)   # Sum of abundances in vector 2
        shared.sp = v1.pa * v2.pa   # Vector of shared species ("pa")
        if(sum(shared.sp) == 0) { 
          res[k,j] = 1
        } else {
          C.j = sum(shared.sp * v1)   # Sum of shared sp. abundances in v1
          C.k = sum(shared.sp * v2)   # Sum of shared sp. abundances in v2
          # a1.j = sum(shared.sp * v1.pa)
          # a1.k = sum(shared.sp * v2.pa)
          a1.j = length(which((shared.sp * v2) == 1)) # Singletons in v2
          a1.k = length(which((shared.sp * v1) == 1)) # Singletons in v1
          a2.j = length(which((shared.sp * v2) == 2)) # Doubletons in v2
          if(a2.j == 0) a2.j <- 1
          a2.k = length(which((shared.sp * v1) == 2)) # Doubletons in v1
          if(a2.k == 0) a2.k <- 1
          # S.j = sum(v1[which(v2 == 1)]) # Sum abund. in v1 for singletons in v2
          # S.k = sum(v2[which(v1 == 1)]) # Sum abund. in v2 for singletons in v1
          sel2 = which(v2 == 1)
          sel1 = which(v1 == 1)
          if(length(sel2)>0) S.j = sum(v1[sel2]) else S.j = 0
          if(length(sel1)>0) S.k = sum(v2[sel1]) else S.k = 0
          
          U.j = (C.j/N.j) + ((N.k-1)/N.k) * (a1.j/(2*a2.j)) * (S.j/N.j) # Eq. 11
          if(U.j > 1) U.j <- 1
          U.k = (C.k/N.k) + ((N.j-1)/N.j) * (a1.k/(2*a2.k)) * (S.k/N.k) # Eq. 12
          if(U.k > 1) U.k <- 1
          
          if(coeff == "Jaccard") {                     # "Jaccard"
            res[k,j] = 1 - (U.j*U.k/(U.j + U.k - U.j*U.k))
          } else if(coeff == "Sorensen") {         # "Sorensen"
            res[k,j] = 1 - (2*U.j*U.k/(U.j + U.k))
          } else if(coeff == "Ochiai") {           # "Ochiai"
            res[k,j] = 1 - (sqrt(U.j*U.k))
          } else if(coeff == "Simpson") { 
            # Simpson (1943), or Lennon et al. (2001) in Chao et al. (2006)
            res[k,j] = 1 -
              (U.j*U.k/(U.j*U.k+min((U.j-U.j*U.k),(U.k-U.j*U.k))))
          } else { # 
            stop("Incorrect coefficient name")
          }
        }
      }
    }
    
  } else {   # Now for complete population data
    
    for(k in 2:nn) {
      for(j in 1:(k-1)) {
        v1 = mat[j,]   # Vector 1
        v2 = mat[k,]   # Vector 2
        v1.pa = decostand(v1,"pa")   # Vector 1 in presence-absence form
        v2.pa = decostand(v2,"pa")   # Vector 2 in presence-absence form
        shared.sp = v1.pa * v2.pa    # Vector of shared species ("pa")
        if(sum(shared.sp) == 0) { 
          res[k,j] = 1
        } else {
          N1 = sum(v1)   # Sum of abundances in vector 1
          N2 = sum(v2)   # Sum of abundances in vector 2
          U = sum(shared.sp * v1)/N1   # Sum of shared sp. abundances in v1
          V = sum(shared.sp * v2)/N2   # Sum of shared sp. abundances in v2
          
          if(coeff == "Jaccard") {                     # "Jaccard"
            res[k,j] = 1 - (U*V/(U + V - U*V))
          } else if(coeff == "Sorensen") {         # "Sorensen"
            res[k,j] = 1 - (2*U*V/(U + V))
          } else if(coeff == "Ochiai") {           # "Ochiai"
            res[k,j] = 1 - (sqrt(U*V))
          } else if(coeff == "Simpson") { # "Simpson"
            res[k,j] = 1 - (U*V/(U*V+min((U-U*V),(V-U*V)))) # Eq. ?
          } else { # 
            stop("Incorrect coefficient name")
          }
        }
      }
    }
  }
  res <- as.dist(res)
}

######## End of beta.div function


beta.div.comp <- function(mat, coef="J", quant=FALSE, save.abc=FALSE)
  #
  # Description --
  # 
  # Podani-family and Baselga-family decompositions of the Jaccard and SÃ¸rensen 
  # dissimilarity coefficients into replacement and richness difference 
  # components, for species presence-absence or abundance data, as described  
  # in Legendre (2014).
  #
  # Usage --
  #
  # beta.div.comp(mat, coef="J", quant=FALSE, save.abc=FALSE)
#
# Arguments --
#
# mat : Data in matrix or data.frame form.
# coef : Family of coefficients to be computed --
#        "S" or "Sorensen": Podani family, SÃ¸rensen-based indices
#        "J" or "Jaccard" : Podani family, Jaccard-based indices
#        "BS" : Baselga family, SÃ¸rensen-based indices
#        "BJ" : Baselga family, Jaccard-based indices
#        "N" : Podani & Schmera (2011) relativized nestedness index.
#        The quantitative form in SÃ¸rensen family is the percentage difference.
#        The quantitative form in the Jaccard family is the Ruzicka index.
#
# quant=TRUE : Compute the quantitative form of replacement, nestedness and D.
#      =FALSE: Compute the presence-absence form of the coefficients.
# save.abc=TRUE : Save the matrices of parameters a, b and c used in the
#      presence-absence calculations.
#
# Details --
#
#    For species presence-absence data, the distance coefficients are 
# Jaccard=(b+c)/(a+b+c) and SÃ¸rensen=(b+c)/(2*a+b+c) with usual abc notation.
#
#    For species abundance data, the distance coefficients are 
# the Ruzicka index = (B+C)/(A+B+C) and Odum's percentage difference 
# (incorrectly called Bray-Curtis) = (B+C)/(2A+B+C), where  
# A = sum of the intersections (or minima) of species abundances at two sites,
# B = sum at site 1 minus A, 
# C = sum at site 2 minus A.
#
#    The binary (quant=FALSE) and quantitative (quant=TRUE) forms of the S and  
# J indices return the same values when computed for presence-absence data.
#
# Value --
#
# repl : Replacement matrix, class = 'dist'.
# rich : Richness/abundance difference or nestedness matrix, class = 'dist'.
#        With options "BJ", "BS" and "N", 'rich' contains nestedness indices.
#        With option "N", the 'repl' and 'rich' values do not add up to 'D'.
# D    : Dissimilarity matrix, class = 'dist'.
# part : Beta diversity partitioning -- 
#        1. Total beta div. = sum(D.ij)/(n*(n-1)) (Legendre & De CÃ¡ceres 2013)
#        2. Total replacement diversity 
#        3. Total richness difference diversity (or nestedness)
#        4. Total replacement div./Total beta div.
#        5. Total richness difference div. (or nestedness)/Total beta div.
# Note : Name of the dissimilarity Scoefficient.
#
# References --
#
# Baselga, A. (2010) Partitioning the turnover and nestedness components of beta 
# diversity. Global Ecology and Biogeography, 19, 134â€“143.
#
# Baselga, A. (2012) The relationship between species replacement, dissimilarity 
# derived from nestedness, and nestedness. Global Ecology and Biogeography, 21, 
# 1223â€“1232. 
#
# Baselga, A. (2013) Separating the two components of abundance-based 
# dissimilarity: balanced changes in abundance vs. abundance gradients. Methods 
# in Ecology and Evolution, 4, 552â€“557.
#
# Carvalho, J.C., Cardoso, P., Borges, P.A.V., Schmera, D. & Podani, J. (2013)
# Measuring fractions of beta diversity and their relationships to nestedness: 
# a theoretical and empirical comparison of novel approaches. Oikos, 122, 
# 825â€“834.
#
# Legendre, P. 2014. Interpreting the replacement and richness difference   
# components of beta diversity. Global Ecology and Biogeography, 23, 1324-1334.
#
# Podani, J., Ricotta, C. & Schmera, D. (2013) A general framework for analyzing 
# beta diversity, nestedness and related community-level phenomena based on 
# abundance data. Ecological Complexity, 15, 52-61.
#
# Podani, J. & Schmera, D. 2011. A new conceptual and methodological framework 
# for exploring and explaining pattern in presence-absence data. Oikos, 120, 
# 1625â€“1638.
#
# License: GPL-2 
# Author:: Pierre Legendre
{
  coef <- pmatch(coef, c("S", "J", "BS", "BJ", "N"))
  #if(coef==5 & quant) stop("coef='N' and quant=TRUE: combination not programmed") #2016-02-02: is now programmed
  mat <- as.matrix(mat)
  n <- nrow(mat)
  if(is.null(rownames(mat))) noms <- paste("Site",1:n,sep="")
  else noms <- rownames(mat)
  #
  if(!quant) {      # Binary data provided, or make the data binary
    if(coef==1) form="Podani family, Sorensen" 
    if(coef==2) form="Podani family, Jaccard"
    if(coef==3) form="Baselga family, Sorensen" 
    if(coef==4) form="Baselga family, Jaccard"
    if(coef==5) form="Podani & Schmera (2011) relativized nestedness"
    mat.b <- ifelse(mat>0, 1, 0)
    a <- mat.b %*% t(mat.b)
    b <- mat.b %*% (1 - t(mat.b))
    c <- (1 - mat.b) %*% t(mat.b)
    min.bc <- pmin(b,c)
    #
    if(coef==1 || coef==2) {
      repl <- 2*min.bc   # replacement, turnover, beta-3
      rich <- abs(b-c)   # nestedness, richness diff., beta-rich
      #
      # Add the denominators
      if(coef==1) {                # SÃ¸rensen-based components
        repl <- repl/(2*a+b+c)
        rich <- rich/(2*a+b+c)
        D <- (b+c)/(2*a+b+c)
      } else if(coef==2) {     # Jaccard-based components
        repl <- repl/(a+b+c)
        rich <- rich/(a+b+c)
        D <- (b+c)/(a+b+c)
      }
    } else if(coef==3) {     # Baselga 2010 components based on SÃ¸rensen
      D <- (b+c)/(2*a+b+c)             # SÃ¸rensen dissimilarity
      repl <- min.bc/(a+min.bc)        # replacement, turnover
      rich <- D-repl                   # richness difference
      
    } else if(coef==4) {      # Baselga 2012 components based on Jaccard
      D <- (b+c)/(a+b+c)               # Jaccard dissimilarity
      repl <- 2*min.bc/(a+2*min.bc)    # replacement, turnover
      rich <- D-repl                   # richness difference
    } else if(coef==5) {      # rich = Podani N = nestdness based on Jaccard
      repl <- 2*min.bc/(a+b+c)
      D <- (b+c)/(a+b+c)
      rich <- matrix(0,n,n)
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          aa = a[i,j]; bb = b[i,j]; cc = c[i,j]
          if(a[i,j] == 0)  rich[i,j] <- 0  
          else  rich[i,j] <- (aa + abs(bb-cc))/(aa+bb+cc) 
        }
      }
    }
    #added code Sven Teurlincx (2015 08 17)
    #this code snippet deals with problems created by entire rows with no data
    D[is.nan(D)]=1
    repl[is.nan(repl)]=1
    rich[is.nan(rich)]=1
    
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    D <- as.dist(D)
    repl <- as.dist(repl)
    rich <- as.dist(rich)
    total.div <- sum(D)/(n*(n-1))
    repl.div <- sum(repl)/(n*(n-1))
    rich.div <- sum(rich)/(n*(n-1))
    part <- c(total.div,repl.div,rich.div,repl.div/total.div,rich.div/total.div)
    #
    if(save.abc) {
      res <- list(repl=repl, rich=rich, D=D, part=part, Note=form, 
                  a=as.dist(a), b=as.dist(b), c=as.dist(c))
    } else { 
      res <- list(repl=repl, rich=rich, D=D, part=part, Note=form)
    }
    #
  } else {      # Quantitative data
    # Calculations based on individuals.within.species
    if(coef==1) form<-"Podani family, percentage difference" 
    if(coef==2) form<-"Podani family, Ruzicka"
    if(coef==3) form<-"Baselga family, percentage difference"
    if(coef==4) form<-"Baselga family, Ruzicka"
    if(coef==5) form<-"Podani & Schmera (2013) relativized nestedness"
    # Baselga (2013) notation:
    # A = W = sum of minima in among-site comparisons
    # B = site.1 sum - W = K.1 - W
    # C = site.2 sum - W = K.2 - W
    K <- vector("numeric", n)   # site (row) sums
    W <- matrix(0,n,n)
    repl <- matrix(0,n,n)
    rich <- matrix(0,n,n)
    D <- matrix(0,n,n)
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    K <- apply(mat,1,sum)         # Row sums
    for(i in 2:n) for(j in 1:(i-1)) W[i,j] <- sum(pmin(mat[i,], mat[j,]))
    #
    # Quantitative extensions of the S and J decompositions
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        repl[i,j] <- 2*(min(K[i],K[j])-W[i,j]) # 2*min(B,C)
        rich[i,j] <- abs(K[i]-K[j])            # abs(B-C)
      }
    }
    #
    # Add the denominators
    if(coef==1) {         # SÃ¸rensen-based (% difference) components
      for(i in 2:n) {
        for(j in 1:(i-1)) {	                        # Baselga 2013 notation:
          repl[i,j] <- repl[i,j]/(K[i]+K[j])          # 2min(B,C)/(2A+B+C)
          rich[i,j] <- rich[i,j]/(K[i]+K[j])          # abs(B-C)/(2A+B+C)
          # cat(K[i], K[j], W[i,j],"\\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j])  # (B+C)/(2A+B+C)
        }
      }
    } else if(coef==2) {    # Jaccard-based (Ruzicka) components
      for(i in 2:n) {
        for(j in 1:(i-1)) {                         # Baselga 2013 notation:
          repl[i,j] <- repl[i,j]/(K[i]+K[j]-W[i,j])   # 2min(B,C)/(A+B+C)
          rich[i,j] <- rich[i,j]/(K[i]+K[j]-W[i,j])   # abs(B-C)/(A+B+C)
          # cat(K[i], K[j], W[i,j],"\\n")
          D[i,j]<-(K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j]) # (B+C)/(A+B+C)
        }
      }
    }
    #
    # Baselga (2013): quantitative extensions of the Baselga (2010) indices
    if(coef==3) {   # Baselga (2013) indices decomposing percentage difference
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          repl[i,j] <- (min(K[i],K[j])-W[i,j])/min(K[i],K[j])
          rich[i,j] <- abs(K[i]-K[j])*W[i,j]/((K[i]+K[j])*min(K[i],K[j]))
          # cat(K[i], K[j], W[i,j],"\\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j])
        }
      }
    }	
    if(coef==4) {   # Decomposing Ruzicka in the spirit of Baselga 2013
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          repl[i,j] <- 
            2*(min(K[i],K[j])-W[i,j])/(2*min(K[i],K[j])-W[i,j])
          rich[i,j] <- abs(K[i]-K[j])*W[i,j]/
            ((K[i]+K[j]-W[i,j])*(2*min(K[i],K[j])-W[i,j]))
          # cat(K[i], K[j], W[i,j],"\\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j])
        }
      }
    }	
    #Podani relative nestedness 
    #added by Sven Teurlincx 2016-02-02
    if(coef==5){      # rich = Podani N = nestdness based on Marczewski-Steinhaus index (1-Ruzicka: quantitive Jaccard) 
      #formula may be found in Podani et al., 2013: A general framework for analyzing beta diversity, nestedness and related community-level phenomena based on abundance data
      for(i in 2:n) {
        for(j in 1:(i-1)) {                         
          repl[i,j] <- repl[i,j]/(K[i]+K[j]-W[i,j])   				
          if(W[i,j] == 0)  rich[i,j] <- 0  
          else  rich[i,j] <- (abs(K[i]-K[j])+W[i,j])/(K[i]+K[j]-W[i,j])  				
          D[i,j]<-(K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j]) 
        }
      }
    }
    #
    #added code Sven Teurlincx (2015 08 17)
    #this code snippet deals with problems created by entire rows with no data
    D[is.nan(D)]=1
    repl[is.nan(repl)]=1
    rich[is.nan(rich)]=1
    
    repl <- as.dist(repl)
    rich <- as.dist(rich)
    D <- as.dist(D)
    repl.div <- sum(repl)/(n*(n-1))
    rich.div <- sum(rich)/(n*(n-1))
    total.div <- sum(D)/(n*(n-1))
    part <- c(total.div,repl.div,rich.div,repl.div/total.div,rich.div/total.div)
    #
    res <- list(repl=repl, rich=rich, D=D, part=part, Note=form)
  }
  res
}

#################################################
# END loading tools
#


#---------------------------------
#---DEFINE DIVERSITY PARTITIONS---
#---------------------------------

ST_divpart = function(dfCOMM, vAREAS_FAC, lAREAS, nQ=0, lSPEC_CONSERV=NA, sINDEX='J'){ 
	require(vegan)
	source("P:/R/Scripts/beta-diversity/beta.div.comp.R")

	vALPHA1				=	c()
	vBETA1_MULT			=	c()
	vBETA1_ADD			=	c()
	vGAMMA1				=	c()
	vBETA1_D_MULT		=	c()
	vBETA1_REPL_MULT	=	c()
	vBETA1_RICH_MULT	=	c()
	vBETA1_REPL_REL_MULT=	c()
	vBETA1_RICH_REL_MULT=	c()
	vBETA1_D			=	c()
	vBETA1_REPL			=	c()
	vBETA1_RICH			=	c()
	vBETA1_REPL_REL		=	c()
	vBETA1_RICH_REL		=	c()
	vBETA1_REPL_SE		=	c()
	vBETA1_RICH_SE		=	c()
	vBETA1_REPL_REL_SE	=	c()
	vBETA1_RICH_REL_SE	=	c()
	vBETA1_TOTAL		=	c()
	vBETA1_CONS_T		=	c()
	vBETA1_CONS_REPL	=	c()
	vBETA1_CONS_RICH	=	c()
	vBETA1_MULT_CONS_D		=	c()
	vBETA1_MULT_CONS_REPL	=	c()
	vBETA1_MULT_CONS_RICH	=	c()
	vALPHA2				=	c()
	vBETA2_MULT			=	c()
	vBETA2_ADD			=	c()
	vGAMMA2				=	c()
	vBETA2_D_MULT		=	c()
	vBETA2_REPL_MULT	=	c()
	vBETA2_RICH_MULT	=	c()
	vBETA2_REPL_REL_MULT=	c()
	vBETA2_RICH_REL_MULT=	c()
	vBETA2_D			=	c()
	vBETA2_REPL			=	c()
	vBETA2_RICH			=	c()
	vBETA2_REPL_REL		=	c()
	vBETA2_RICH_REL		=	c()
	vBETA2_REPL_SE		=	c()
	vBETA2_RICH_SE		=	c()
	vBETA2_REPL_REL_SE	=	c()
	vBETA2_RICH_REL_SE	=	c()
	vBETA2_TOTAL		=	c()
	vBETA2_CONS_T		=	c()
	vBETA2_CONS_REPL	=	c()
	vBETA2_CONS_RICH	=	c()
	vBETA2_MULT_CONS_D		=	c()
	vBETA2_MULT_CONS_REPL	=	c()
	vBETA2_MULT_CONS_RICH	=	c()
	
	
	for(sAREA in lAREAS){
		dbSPEC_AREA = subset(dfCOMM, vAREAS_FAC %in% sAREA)
		if(length(which(colSums(dbSPEC_AREA) == 0))>0){dbSPEC_AREA = dbSPEC_AREA[,-(which(colSums(dbSPEC_AREA) == 0)),drop=F]} #remove all zero columns		
			#---compute alpha, beta and gamma diversity at within polder level---
			#binary data (presence/absence)
			if(nQ == 0){
				#create binary matrix for use as presence absence data
				dbSPEC_SEL = (dbSPEC_AREA >0) + 0 #presence/absence matrix
				# calculation of gamma diversity (=sum of species over all sites)		
				fGAMMA1 		= 	sum((apply(dbSPEC_SEL>0,2,sum))>0)

				# Calculation of species richness(alpha diversity)		
				fALPHA1 		= 	mean(apply(dbSPEC_SEL>0,1,sum))		
						
				#--Compute the partitions of beta diversity--
				if(nrow(dbSPEC_AREA)>1){
					#compute multipart indices (as per Baselga 2013 / Ensing 2015)
					#	only works for presence/absence data
					#	not implemented for sorensen based podani indeces
					if(sINDEX=='BJ'){					
						bdSPEC_S_MULT 		= beta.multi.carv(dbSPEC_SEL, index='jaccard')
						fBETA1_D_MULT 		= bdSPEC_S_MULT[[3]]
						fBETA1_REPL_MULT 	= bdSPEC_S_MULT[[1]]
						fBETA1_RICH_MULT 	= bdSPEC_S_MULT[[2]]
					}else if(sINDEX=='BS'){
						bdSPEC_S_MULT 		= beta.multi.carv(dbSPEC_SEL, index='sorensen')
						fBETA1_D_MULT 		= bdSPEC_S_MULT[[3]]
						fBETA1_REPL_MULT 	= bdSPEC_S_MULT[[1]]
						fBETA1_RICH_MULT 	= bdSPEC_S_MULT[[2]]
					}else if(sINDEX=='J'){
						bdSPEC_S_MULT 		= beta.multi.carv(dbSPEC_SEL, index='carvalho')
						fBETA1_D_MULT 		= bdSPEC_S_MULT[[3]]
						fBETA1_REPL_MULT 	= bdSPEC_S_MULT[[1]]
						fBETA1_RICH_MULT 	= bdSPEC_S_MULT[[2]]
					}else{
						fBETA1_D_MULT 		= NA
						fBETA1_REPL_MULT 	= NA
						fBETA1_RICH_MULT 	= NA
					}
					
					#Using podani family for partitions with a Sorensen-based index on presence/absence data
					bdSPEC_S = beta.div.comp(dbSPEC_SEL, coef=sINDEX, quant=FALSE, save.abc=TRUE) 
					fBETA1_TOTAL	=	mean(bdSPEC_S$b+bdSPEC_S$c)# #equivalent to multiplicative beta diversity (see Legendre paper)
					fBETA1_D= bdSPEC_S$part[1]#is the same as sum(bdSPEC_S$D)/(nrow(as.matrix(bdSPEC_S$D))*(nrow(as.matrix(bdSPEC_S$D))-1))
					fBETA1_D_SE = sd(bdSPEC_S$D)/(nrow(as.matrix(bdSPEC_S$D))*(nrow(as.matrix(bdSPEC_S$D))-1))
					fBETA1_REPL=bdSPEC_S$part[2]
					fBETA1_REPL_SE = sd(bdSPEC_S$repl)/sqrt(nrow(as.matrix(bdSPEC_S$repl))*(nrow(as.matrix(bdSPEC_S$repl))-1))
					fBETA1_REPL_REL_SE = (sd(bdSPEC_S$repl/bdSPEC_S$D)/sqrt(nrow(as.matrix(bdSPEC_S$D))*(nrow(as.matrix(bdSPEC_S$D))-1)))*100
					fBETA1_RICH=bdSPEC_S$part[3]
					fBETA1_RICH_SE = sd(bdSPEC_S$rich)/sqrt(nrow(as.matrix(bdSPEC_S$D))*(nrow(as.matrix(bdSPEC_S$D))-1))
					fBETA1_RICH_REL_SE = (sd(bdSPEC_S$rich/bdSPEC_S$D)/sqrt(nrow(as.matrix(bdSPEC_S$D))*(nrow(as.matrix(bdSPEC_S$D))-1)))*100
				}else{
					fBETA1_TOTAL		=	NA	
					fBETA1_D			=	NA
					fBETA1_REPL			=	NA
					fBETA1_REPL_SE 		= 	NA
					fBETA1_REPL_REL_SE 	=	NA
					fBETA1_RICH			=	NA
					fBETA1_RICH_SE 		= 	NA
					fBETA1_RICH_REL_SE 	=	NA
					fBETA1_D_MULT 		= NA
					fBETA1_REPL_MULT 	= NA
					fBETA1_RICH_MULT 	= NA
				}
				
				if(is.na(lSPEC_CONSERV)[1]==FALSE){
					#CODE for splitting Legendre partitions of beta diversity into conservation species pattern and non-conservation species pattern
					dbSPEC_SEL_CONSERV = dbSPEC_SEL[, which(colnames(dbSPEC_SEL) %in% lSPEC_CONSERV), drop=F]
					bdSPEC_S = beta.div.comp(dbSPEC_SEL, coef=sINDEX, quant=FALSE, save.abc=TRUE) 
					bdSPEC_S_CONS = beta.div.comp(dbSPEC_SEL_CONSERV, coef=sINDEX, quant=FALSE, save.abc=TRUE) 				
					#--numerator replacement difference: 2*min(b,c)
					mSPEC_S_REPLFORM1=2*pmin(as.matrix(bdSPEC_S$b), as.matrix(bdSPEC_S$c))
					distSPEC_S_REPLFORM1 = as.dist(mSPEC_S_REPLFORM1)
					mSPEC_S_CONS_REPLFORM1=2*pmin(as.matrix(bdSPEC_S_CONS$b), as.matrix(bdSPEC_S_CONS$c))
					distSPEC_S_CONS_REPLFORM1 = as.dist(mSPEC_S_CONS_REPLFORM1)
					#--operational part richness difference: |b-c|
					mSPEC_S_RICHFORM1 = abs(as.matrix(bdSPEC_S$b) - as.matrix(bdSPEC_S$c))
					distSPEC_S_RICHFORM1 = as.dist(mSPEC_S_RICHFORM1)
					mSPEC_S_CONS_RICHFORM1 = abs(as.matrix(bdSPEC_S_CONS$b) - as.matrix(bdSPEC_S_CONS$c))
					distSPEC_S_CONS_RICHFORM1 = as.dist(mSPEC_S_CONS_RICHFORM1)
					
				
					#calculate community comparison components a, b and c (Legendre) and components a', b' and c' based on conservation relevant community only				
					bprime=as.matrix(bdSPEC_S_CONS$b)
					cprime=as.matrix(bdSPEC_S_CONS$c)
					b=as.matrix(bdSPEC_S$b)
					c=as.matrix(bdSPEC_S$c)
					
					#--REPL: Conservation species fraction--
					#calculate chance of finding conservation species in the richest site (site with overlap of rich and repl)
				#	(pmax(bprime, cprime)/pmax(b,c))
					#multiply by the number of species in richest site to find the expected number of conservation species in the richest site
				#	(pmax(bprime, cprime)/pmax(b,c))*pmin(b,c)
					#add the identity specific species that are replaced from the poor site (1*min(b',c')) to get the total number of expected conservation species in the REPL partition
					dfNC_REPL1 = ((pmax(bprime, cprime)/pmax(b,c))*pmin(b,c))+1*pmin(bprime,cprime)#expected number of conservation species in REPL
					dfNC_REPL1_FRAC = (((pmax(bprime, cprime)/pmax(b,c))*pmin(b,c))+1*pmin(bprime,cprime))/(2*pmin(b,c)) #fraction of conservation species to total species in replacement partition
					dfNC_REPL1[is.na(dfNC_REPL1)]=0#1 because if there are no species of conservation interest and no species where different we can say that all of the species that could be of conservation interest are identical

					
					#--RICH: Conservation species fraction--
					#calculate chance of finding conservation species in the richest site (site with overlap of rich and repl)
				#	(pmax(bprime, cprime)/pmax(b,c))
					#Get richness difference of total community
				#	abs(b-c)
					#multiply conservation species chance with the total species involved in the richness fraction
					dfNC_RICH1 = (pmax(bprime, cprime)/pmax(b,c))*abs(b-c)
					dfNC_RICH1_FRAC = ((pmax(bprime, cprime)/pmax(b,c))*abs(b-c))/abs(b-c) #fraction of conservation species to total species in replacement partition
					dfNC_RICH1[is.na(dfNC_RICH1)]=0
					
					dfNC_BETA1 = bprime+cprime
					#calculate mean and se
					#dfNC_REPL_FRAC_SITE = apply(dfNC_REPL_FRAC,1,FUN=mean, na.rm=TRUE)
					#mean(dfNC_REPL_FRAC_SITE)
					#sd(dfNC_REPL_FRAC_SITE)/sqrt(length(dfNC_REPL_FRAC_SITE))
					#dfNC_RICH_FRAC_SITE = apply(dfNC_RICH_FRAC,1,FUN=mean, na.rm=TRUE)
					#mean(dfNC_RICH_FRAC_SITE)
					#sd(dfNC_RICH_FRAC_SITE)/sqrt(length(dfNC_RICH_FRAC_SITE))
					fBETA1_CONS_T		=	mean(as.dist(dfNC_BETA1))
					fBETA1_CONS_REPL 	= 	mean(as.dist(dfNC_REPL1))
					fBETA1_CONS_RICH 	= 	mean(as.dist(dfNC_RICH1))
					
					#calculate multi site metrics for conservation species
					if(sINDEX=='BJ'){					
						bdSPEC_S_CONS_MULT 		=	beta.multi.carv(dbSPEC_SEL_CONSERV, index='jaccard')
						fBETA1_MULT_CONS_D 		=	bdSPEC_S_CONS_MULT[[3]]
						fBETA1_MULT_CONS_REPL 	=	bdSPEC_S_CONS_MULT[[1]]
						fBETA1_MULT_CONS_RICH 	=	bdSPEC_S_CONS_MULT[[2]]
					}else if(sINDEX=='BS'){	
						bdSPEC_S_CONS_MULT 		=	beta.multi.carv(dbSPEC_SEL_CONSERV, index='sorensen')
						fBETA1_MULT_CONS_D 		=	bdSPEC_S_CONS_MULT[[3]]
						fBETA1_MULT_CONS_REPL 	=	bdSPEC_S_CONS_MULT[[1]]
						fBETA1_MULT_CONS_RICH 	=	bdSPEC_S_CONS_MULT[[2]]
					}else if(sINDEX=='J'){	
						bdSPEC_S_CONS_MULT 		=	beta.multi.carv(dbSPEC_SEL_CONSERV, index='carvalho')
						fBETA1_MULT_CONS_D 		=	bdSPEC_S_CONS_MULT[[3]]
						fBETA1_MULT_CONS_REPL 	=	bdSPEC_S_CONS_MULT[[1]]
						fBETA1_MULT_CONS_RICH 	=	bdSPEC_S_CONS_MULT[[2]]
					}else{
						fBETA1_MULT_CONS_D 		=	NA
						fBETA1_MULT_CONS_REPL 	=	NA
						fBETA1_MULT_CONS_RICH 	=	NA
					}
				
				}else{
					fBETA1_CONS_T			=	NA
					fBETA1_CONS_REPL 		= 	NA
					fBETA1_CONS_RICH 		= 	NA
					fBETA1_MULT_CONS_D 		=	NA
					fBETA1_MULT_CONS_REPL 	=	NA
					fBETA1_MULT_CONS_RICH 	=	NA
				}
				
			}else if(nQ == 1){	
				# conversion data to fractions
				dbSPEC_SEL=decostand(dbSPEC_AREA,'total')
				if(length(which(colSums(dbSPEC_SEL) == 0))>0){dbSPEC_SEL = dbSPEC_SEL[,-(which(colSums(dbSPEC_SEL) == 0))]}  #remove all zero columns
				if(length(which(rowSums(dbSPEC_SEL) == 0))>0){dbSPEC_SEL = dbSPEC_SEL[-(which(rowSums(dbSPEC_SEL) == 0)),]}  #remove all zero rows	
				nSITES=nrow(dbSPEC_SEL)

				# calculation of D_GAMMA (Jost 2007 equation 17b)
				fGAMMA1=exp(sum(-(colSums(dbSPEC_SEL)/nSITES)*log(colSums(dbSPEC_SEL)/nSITES)))

				# calculation of D_ALPHA (Jost 2007 equation 14 or 17a)
				fALPHA1=exp(1/nSITES*sum(vegan::diversity(dbSPEC_SEL)))
				
				#--Compute the partitions of beta diversity--
				#Using Baselga family for partitions with a Jaccard-based index
				#As we can't easily do this on a Horn index.
				#Hence we have to assume that the baselga method yields similar results in terms of fractions.
				#Nor should we, see Baselga 2012 for comments on Horn partitioning not having been researched and being problematic
				fBETA1_TOTAL	=	NA
				if(nrow(dbSPEC_AREA)>1){
					#compute multipart indices (as per Baselga 2013 / Ensing 2015)
					#	only works for presence/absence data
					#	not implemented for sorensen based podani indeces
					fBETA1_D_MULT 		= NA
					fBETA1_REPL_MULT 	= NA
					fBETA1_RICH_MULT 	= NA
				
					bdSPEC_D = beta.div.comp(dbSPEC_AREA, coef=sINDEX, quant=TRUE, save.abc=TRUE)
					fBETA1_D= bdSPEC_D$part[1]#is the same as sum(bdSPEC_D$D)/(nrow(as.matrix(bdSPEC_D$D))*(nrow(as.matrix(bdSPEC_D$D))-1))
					fBETA1_D_SE = sd(bdSPEC_D$D)/(nrow(as.matrix(bdSPEC_D$D))*(nrow(as.matrix(bdSPEC_D$D))-1))
					fBETA1_REPL=bdSPEC_D$part[2]
					fBETA1_REPL_SE = sd(bdSPEC_D$repl)/sqrt(nrow(as.matrix(bdSPEC_D$repl))*(nrow(as.matrix(bdSPEC_D$repl))-1))
					fBETA1_REPL_REL_SE = (sd(bdSPEC_D$repl/bdSPEC_D$D)/sqrt(nrow(as.matrix(bdSPEC_D$D))*(nrow(as.matrix(bdSPEC_D$D))-1)))*100
					fBETA1_RICH=bdSPEC_D$part[3]
					fBETA1_RICH_SE = sd(bdSPEC_D$rich)/sqrt(nrow(as.matrix(bdSPEC_D$D))*(nrow(as.matrix(bdSPEC_D$D))-1))
					fBETA1_RICH_REL_SE = (sd(bdSPEC_D$rich/bdSPEC_D$D)/sqrt(nrow(as.matrix(bdSPEC_D$D))*(nrow(as.matrix(bdSPEC_D$D))-1)))*100
				}else{
					fBETA1_D			=	NA
					fBETA1_REPL			=	NA
					fBETA1_REPL_SE 		= 	NA
					fBETA1_REPL_REL_SE 	=	NA
					fBETA1_RICH			=	NA
					fBETA1_RICH_SE 		= 	NA
					fBETA1_RICH_REL_SE 	=	NA
					fBETA1_D_MULT 		= 	NA
					fBETA1_REPL_MULT 	= 	NA
					fBETA1_RICH_MULT 	= 	NA
				}
				fBETA1_CONS_T		=	NA
				fBETA1_CONS_REPL 	= 	NA
				fBETA1_CONS_RICH 	= 	NA
				fBETA1_MULT_CONS_D 		=	NA
				fBETA1_MULT_CONS_REPL 	=	NA
				fBETA1_MULT_CONS_RICH 	=	NA
			}		
		
		fBETA1_MULT = fGAMMA1/fALPHA1
		fBETA1_ADD = fGAMMA1-fALPHA1
		fBETA1_REPL_REL = (fBETA1_REPL / fBETA1_D) *100
		fBETA1_RICH_REL = (fBETA1_RICH / fBETA1_D) *100	
		
		vALPHA1				=	c(vALPHA1			,	fALPHA1				)
		vBETA1_MULT			=	c(vBETA1_MULT		,	fBETA1_MULT			)
		vBETA1_ADD			=	c(vBETA1_ADD		,	fBETA1_ADD			)
		vGAMMA1				=	c(vGAMMA1			,	fGAMMA1				)
		vBETA1_D_MULT		=	c(vBETA1_D_MULT		,	fBETA1_D_MULT		)
		vBETA1_REPL_REL_MULT=	c(vBETA1_REPL_REL_MULT	,	(fBETA1_REPL_MULT/fBETA1_D_MULT)*100	)
		vBETA1_RICH_REL_MULT=	c(vBETA1_RICH_REL_MULT	,	(fBETA1_RICH_MULT/fBETA1_D_MULT)*100	)
		vBETA1_REPL_MULT	=	c(vBETA1_REPL_MULT	,	fBETA1_REPL_MULT	)
		vBETA1_RICH_MULT	=	c(vBETA1_RICH_MULT	,	fBETA1_RICH_MULT	)
		vBETA1_D			=	c(vBETA1_D			,	fBETA1_D			)
		vBETA1_REPL			=	c(vBETA1_REPL		,	fBETA1_REPL			)
		vBETA1_RICH			=	c(vBETA1_RICH		,	fBETA1_RICH			)
		vBETA1_REPL_REL		=	c(vBETA1_REPL_REL	,	fBETA1_REPL_REL		)
		vBETA1_RICH_REL		=	c(vBETA1_RICH_REL	,	fBETA1_RICH_REL		)
		vBETA1_REPL_SE		=	c(vBETA1_REPL_SE	,	fBETA1_REPL_SE		)
		vBETA1_RICH_SE		=	c(vBETA1_RICH_SE	,	fBETA1_RICH_SE		)
		vBETA1_REPL_REL_SE	=	c(vBETA1_REPL_REL_SE,	fBETA1_REPL_REL_SE	)
		vBETA1_RICH_REL_SE	=	c(vBETA1_RICH_REL_SE,	fBETA1_RICH_REL_SE	)
		vBETA1_TOTAL		=	c(vBETA1_TOTAL		,	fBETA1_TOTAL		)
		vBETA1_CONS_T		=	c(vBETA1_CONS_T		,	fBETA1_CONS_T		)
		vBETA1_CONS_REPL	=	c(vBETA1_CONS_REPL	,	fBETA1_CONS_REPL	)
		vBETA1_CONS_RICH	=	c(vBETA1_CONS_RICH	,	fBETA1_CONS_RICH	)
		vBETA1_MULT_CONS_D		=	c(vBETA1_MULT_CONS_D	,	fBETA1_MULT_CONS_D 		)
		vBETA1_MULT_CONS_REPL	=	c(vBETA1_MULT_CONS_REPL	,	fBETA1_MULT_CONS_REPL	)
		vBETA1_MULT_CONS_RICH	=	c(vBETA1_MULT_CONS_RICH	,	fBETA1_MULT_CONS_RICH	)
		
	}#end of for sAREA loop

#-----BETWEEN POLDER ALPHA2 BETA2 GAMMA2----

	#---calculate between polder area alpha beta and gamma---
		#make polder aggregated dfs
		dbSPEC_SEL_AREAS =as.data.frame(dfCOMM)
		dbSPEC_SEL_AREAS$AREA=vAREAS_FAC
		dfSPEC_POLDER = aggregate(dbSPEC_SEL_AREAS[,-ncol(dbSPEC_SEL_AREAS)], list(dbSPEC_SEL_AREAS$AREA), sum)
		rownames(dfSPEC_POLDER) = dfSPEC_POLDER[,1]
		dfSPEC_POLDER1 = dfSPEC_POLDER[,-c(1)]
		if(length(which(colSums(dfSPEC_POLDER1) == 0))>0){dfSPEC_POLDER1 = dfSPEC_POLDER1[,-(which(colSums(dfSPEC_POLDER1) == 0)),drop=F]} #remove all zero columns		
		
		#---compute alpha, beta and gamma diversity at BETWEEN polder level---
		#binary data (presence/absence)
		if(nQ == 0){
			#create binary matrix for use as presence absence data
			dbSPEC_SEL2 = (dfSPEC_POLDER1 >0) + 0 #presence/absence matrix
			# calculation of gamma diversity (=sum of species over all sites)		
			fGAMMA2 	= 	sum((apply(dbSPEC_SEL2>0,2,sum))>0)

			# Calculation of species richness(alpha diversity)		
			fALPHA2 	= 	mean(apply(dbSPEC_SEL2>0,1,sum))		
			#--Compute the partitions of beta diversity--
			
			#compute multipart indices (as per Baselga 2013 / Ensing 2015)
			#	only works for presence/absence data
			#	not implemented for sorensen based podani indeces
			if(sINDEX=='BJ'){					
				bdSPEC_S_MULT 		= beta.multi.carv(dbSPEC_SEL2, index='jaccard')
				fBETA2_D_MULT 		= bdSPEC_S_MULT[[3]]
				fBETA2_REPL_MULT 	= bdSPEC_S_MULT[[1]]
				fBETA2_RICH_MULT 	= bdSPEC_S_MULT[[2]]
			}else if(sINDEX=='BS'){
				bdSPEC_S_MULT 		= beta.multi.carv(dbSPEC_SEL2, index='sorensen')
				fBETA2_D_MULT 		= bdSPEC_S_MULT[[3]]
				fBETA2_REPL_MULT 	= bdSPEC_S_MULT[[1]]
				fBETA2_RICH_MULT 	= bdSPEC_S_MULT[[2]]
			}else if(sINDEX=='J'){
				bdSPEC_S_MULT 		= beta.multi.carv(dbSPEC_SEL2, index='carvalho')
				fBETA2_D_MULT 		= bdSPEC_S_MULT[[3]]
				fBETA2_REPL_MULT 	= bdSPEC_S_MULT[[1]]
				fBETA2_RICH_MULT 	= bdSPEC_S_MULT[[2]]
			}else{
				fBETA2_D_MULT 		= NA
				fBETA2_REPL_MULT 	= NA
				fBETA2_RICH_MULT 	= NA
			}

			#Using podani family for partitions with a Sorensen-based index on presence/absence data
			 bdSPEC_S2 			= 	beta.div.comp(dbSPEC_SEL2, coef=sINDEX, quant=FALSE, save.abc=TRUE) 
			 fBETA2_TOTAL		=	mean(bdSPEC_S2$b+bdSPEC_S2$c)#bdSPEC_S$part[1] #equivalent to multiplicative beta diversity (see Legendre paper)
			 fBETA2_D= bdSPEC_S2$part[1]#is the same as sum(bdSPEC_S2$D)/(nrow(as.matrix(bdSPEC_S2$D))*(nrow(as.matrix(bdSPEC_S2$D))-1))
			 fBETA2_D_SE = sd(bdSPEC_S2$D)/(nrow(as.matrix(bdSPEC_S2$D))*(nrow(as.matrix(bdSPEC_S2$D))-1))
			 fBETA2_REPL=bdSPEC_S2$part[2]
			 fBETA2_REPL_SE = sd(bdSPEC_S2$repl)/sqrt(nrow(as.matrix(bdSPEC_S2$repl))*(nrow(as.matrix(bdSPEC_S2$repl))-1))
			 fBETA2_REPL_REL_SE = (sd(bdSPEC_S2$repl/bdSPEC_S2$D)/sqrt(nrow(as.matrix(bdSPEC_S2$D))*(nrow(as.matrix(bdSPEC_S2$D))-1)))*100
			 fBETA2_RICH=bdSPEC_S2$part[3]
			 fBETA2_RICH_SE = sd(bdSPEC_S2$rich)/sqrt(nrow(as.matrix(bdSPEC_S2$D))*(nrow(as.matrix(bdSPEC_S2$D))-1))
			 fBETA2_RICH_REL_SE = (sd(bdSPEC_S2$rich/bdSPEC_S2$D)/sqrt(nrow(as.matrix(bdSPEC_S2$D))*(nrow(as.matrix(bdSPEC_S2$D))-1)))*100
		    
			if(is.na(lSPEC_CONSERV[1])==FALSE){			
				 #CODE for splitting Legendre partitions of beta diversity into conservation species pattern and non-conservation species pattern
				 dbSPEC_SEL2_CONSERV = dbSPEC_SEL2[, which(colnames(dbSPEC_SEL2) %in% lSPEC_CONSERV), drop=F]
				 bdSPEC_S2 = beta.div.comp(dbSPEC_SEL2, coef=sINDEX, quant=FALSE, save.abc=TRUE) 
				 bdSPEC_S2_CONS = beta.div.comp(dbSPEC_SEL2_CONSERV, coef=sINDEX, quant=FALSE, save.abc=TRUE) 				
				 #--numerator replacement difference: 2*min(b,c)
				 mSPEC_S2_REPLFORM1=2*pmin(as.matrix(bdSPEC_S2$b), as.matrix(bdSPEC_S2$c))
				 distSPEC_S2_REPLFORM1 = as.dist(mSPEC_S2_REPLFORM1)
				 mSPEC_S2_CONS_REPLFORM1=2*pmin(as.matrix(bdSPEC_S2_CONS$b), as.matrix(bdSPEC_S2_CONS$c))
				 distSPEC_S2_CONS_REPLFORM1 = as.dist(mSPEC_S2_CONS_REPLFORM1)
				 #--operational part richness difference: |b-c|
				 mSPEC_S2_RICHFORM1 = abs(as.matrix(bdSPEC_S2$b) - as.matrix(bdSPEC_S2$c))
				 distSPEC_S2_RICHFORM1 = as.dist(mSPEC_S2_RICHFORM1)
				 mSPEC_S2_CONS_RICHFORM1 = abs(as.matrix(bdSPEC_S2_CONS$b) - as.matrix(bdSPEC_S2_CONS$c))
				 distSPEC_S2_CONS_RICHFORM1 = as.dist(mSPEC_S2_CONS_RICHFORM1)
				 
				 #calculate community comparison components a, b and c (Legendre) and components a', b' and c' based on conservation relevant community only				
				 bprime=as.matrix(bdSPEC_S2_CONS$b)
				 cprime=as.matrix(bdSPEC_S2_CONS$c)
				 b=as.matrix(bdSPEC_S2$b)
				 c=as.matrix(bdSPEC_S2$c)
				 
				 #--REPL: Conservation species fraction--
				 #calculate chance of finding conservation species in the richest site (site with overlap of rich and repl)
				# (pmax(bprime, cprime)/pmax(b,c))
				 #multiply by the number of species in richest site to find the expected number of conservation species in the richest site
				# (pmax(bprime, cprime)/pmax(b,c))*pmin(b,c)
				 #add the identity specific species that are replaced from the poor site (1*min(b',c')) to get the total number of expected conservation species in the REPL partition
				 dfNC_REPL2 = ((pmax(bprime, cprime)/pmax(b,c))*pmin(b,c))+1*pmin(bprime,cprime)#expected number of conservation species in REPL
				 dfNC_REPL2_FRAC = (((pmax(bprime, cprime)/pmax(b,c))*pmin(b,c))+1*pmin(bprime,cprime))/(2*pmin(b,c)) #fraction of conservation species to total species in replacement partition
				 dfNC_REPL2[is.na(dfNC_REPL2)]=0#1 because if there are no species of conservation interest and no species where different we can say that all of the species that could be of conservation interest are identical
				 
				 #--RICH: Conservation species fraction--
				 #calculate chance of finding conservation species in the richest site (site with overlap of rich and repl)
				# (pmax(bprime, cprime)/pmax(b,c))
				 #Get richness difference of total community
				# abs(b-c)
				 #multiply conservation species chance with the total species involved in the richness fraction
				 dfNC_RICH2 = (pmax(bprime, cprime)/pmax(b,c))*abs(b-c)
				 dfNC_RICH2_FRAC = ((pmax(bprime, cprime)/pmax(b,c))*abs(b-c))/abs(b-c) #fraction of conservation species to total species in replacement partition
				 dfNC_RICH2[is.na(dfNC_RICH2)]=0
				 
				 dfNC_BETA2 = bprime+cprime	
				 
				 
				fBETA2_CONS_T		=	mean(as.dist(dfNC_BETA2))
				fBETA2_CONS_REPL 	= 	mean(as.dist(dfNC_REPL2))
				fBETA2_CONS_RICH 	= 	mean(as.dist(dfNC_RICH2))
				 
				#calculate multi site metrics for conservation species
				if(sINDEX=='BJ'){					
					bdSPEC_S2_CONS_MULT 		=	beta.multi.carv(dbSPEC_SEL2_CONSERV, index='jaccard')
					fBETA2_MULT_CONS_D 		=	bdSPEC_S2_CONS_MULT[[3]]
					fBETA2_MULT_CONS_REPL 	=	bdSPEC_S2_CONS_MULT[[1]]
					fBETA2_MULT_CONS_RICH 	=	bdSPEC_S2_CONS_MULT[[2]]
				}else if(sINDEX=='BS'){	
					bdSPEC_S2_CONS_MULT 		=	beta.multi.carv(dbSPEC_SEL2_CONSERV, index='sorensen')
					fBETA2_MULT_CONS_D 		=	bdSPEC_S2_CONS_MULT[[3]]
					fBETA2_MULT_CONS_REPL 	=	bdSPEC_S2_CONS_MULT[[1]]
					fBETA2_MULT_CONS_RICH 	=	bdSPEC_S2_CONS_MULT[[2]]
				}else if(sINDEX=='J'){	
					bdSPEC_S2_CONS_MULT 		=	beta.multi.carv(dbSPEC_SEL2_CONSERV, index='carvalho')
					fBETA2_MULT_CONS_D 		=	bdSPEC_S2_CONS_MULT[[3]]
					fBETA2_MULT_CONS_REPL 	=	bdSPEC_S2_CONS_MULT[[1]]
					fBETA2_MULT_CONS_RICH 	=	bdSPEC_S2_CONS_MULT[[2]]
				}else{
					fBETA2_MULT_CONS_D 		=	NA
					fBETA2_MULT_CONS_REPL 	=	NA
					fBETA2_MULT_CONS_RICH 	=	NA
				}
				 
			}else{
				fBETA2_CONS_T			=	NA
				fBETA2_CONS_REPL 		= 	NA
				fBETA2_CONS_RICH 		= 	NA
				fBETA2_MULT_CONS_D 		=	NA
				fBETA2_MULT_CONS_REPL 	=	NA
				fBETA2_MULT_CONS_RICH 	=	NA
			}
		
		}else if(nQ == 1){	
			# conversion data to fractions
			dbSPEC_SEL2=decostand(dfSPEC_POLDER1,'total')
			if(length(which(colSums(dbSPEC_SEL2) == 0))>0){dbSPEC_SEL2 = dbSPEC_SEL2[,-(which(colSums(dbSPEC_SEL2) == 0))]}  #remove all zero columns
			if(length(which(rowSums(dbSPEC_SEL2) == 0))>0){dbSPEC_SEL2 = dbSPEC_SEL2[-(which(rowSums(dbSPEC_SEL2) == 0)),]}  #remove all zero rows	
			nSITES2=nrow(dbSPEC_SEL2)

			# calculation of D_GAMMA (Jost 2007 equation 17b)
			fGAMMA2	=	exp(sum(-(colSums(dbSPEC_SEL2)/nSITES2)*log(colSums(dbSPEC_SEL2)/nSITES2)))

			# calculation of D_ALPHA (Jost 2007 equation 14 or 17a)
			fALPHA2		=	exp(1/nSITES2*sum(vegan::diversity(dbSPEC_SEL2)))
			
			#compute multipart indices (as per Baselga 2013 / Ensing 2015)
			#	only works for presence/absence data
			#	not implemented for sorensen based podani indeces
			fBETA2_D_MULT 		= NA
			fBETA2_REPL_MULT 	= NA
			fBETA2_RICH_MULT 	= NA
			
			#--Compute the partitions of beta diversity--
			#Using Baselga family for partitions with a Sorensen-based index
			#As we can't easily do this on a Horn index.
			#Hence we have to assume that the baselga method yields similar results in terms of fractions. !!!Should think more about ways to do this in the future.
			#We can however calculate the horn index and use it here for total beta diversity and then multiply the fractions by the horn index totalBD
			if(nrow(dbSPEC_SEL2)>1){
				bdSPEC_D2 			= 	beta.div.comp(dfSPEC_POLDER1, coef=sINDEX, quant=TRUE, save.abc=TRUE)
				fBETA2_TOTAL		=	NA
				fBETA2_D			= 	bdSPEC_D2$part[1]#is the same as sum(bdSPEC_D2$D)/(nrow(as.matrix(bdSPEC_D2$D))*(nrow(as.matrix(bdSPEC_D2$D))-1))
				fBETA2_D_SE 		= 	sd(bdSPEC_D2$D)/(nrow(as.matrix(bdSPEC_D2$D))*(nrow(as.matrix(bdSPEC_D2$D))-1))
				fBETA2_REPL			=	bdSPEC_D2$part[2]
				fBETA2_REPL_SE 		= 	sd(bdSPEC_D2$repl)/sqrt(nrow(as.matrix(bdSPEC_D2$repl))*(nrow(as.matrix(bdSPEC_D2$repl))-1))
				fBETA2_REPL_REL_SE 	= 	(sd(bdSPEC_D2$repl/bdSPEC_D2$D)/sqrt(nrow(as.matrix(bdSPEC_D2$D))*(nrow(as.matrix(bdSPEC_D2$D))-1)))*100
				fBETA2_RICH			=	bdSPEC_D2$part[3]
				fBETA2_RICH_SE 		= 	sd(bdSPEC_D2$rich)/sqrt(nrow(as.matrix(bdSPEC_D2$D))*(nrow(as.matrix(bdSPEC_D2$D))-1))
				fBETA2_RICH_REL_SE 	= 	(sd(bdSPEC_D2$rich/bdSPEC_D2$D)/sqrt(nrow(as.matrix(bdSPEC_D2$D))*(nrow(as.matrix(bdSPEC_D2$D))-1)))*100
			}else{
				bdSPEC_D2 			= NA
				fBETA2_TOTAL		= NA
				fBETA2_D			= NA
				fBETA2_D_SE 		= NA
				fBETA2_REPL			= NA
				fBETA2_REPL_SE 		= NA
				fBETA2_REPL_REL_SE 	= NA
				fBETA2_RICH			= NA
				fBETA2_RICH_SE 		= NA
				fBETA2_RICH_REL_SE 	= NA
			}
			fBETA2_CONS_T		=	NA
			fBETA2_CONS_REPL 	= 	NA
			fBETA2_CONS_RICH 	= 	NA
			fBETA2_MULT_CONS_D 		=	NA
			fBETA2_MULT_CONS_REPL 	=	NA
			fBETA2_MULT_CONS_RICH 	=	NA
		
		}	
		fBETA2_MULT = fGAMMA2/fALPHA2
		fBETA2_ADD = fGAMMA2-fALPHA2
		fBETA2_REPL_REL = (fBETA2_REPL / fBETA2_D) *100
		fBETA2_RICH_REL = (fBETA2_RICH / fBETA2_D) *100	
		
		
	vALPHA2				=	c(vALPHA2			,	fALPHA2				)
	vBETA2_MULT			=	c(vBETA2_MULT		,	fBETA2_MULT			)
	vBETA2_ADD			=	c(vBETA2_ADD		,	fBETA2_ADD			)
	vGAMMA2				=	c(vGAMMA2			,	fGAMMA2				)
	vBETA2_D_MULT		=	c(vBETA2_D_MULT		,	fBETA2_D_MULT		)
	vBETA2_REPL_MULT	=	c(vBETA2_REPL_MULT	,	fBETA2_REPL_MULT	)
	vBETA2_RICH_MULT	=	c(vBETA2_RICH_MULT	,	fBETA2_RICH_MULT	)
	vBETA2_REPL_REL_MULT	=	c(vBETA2_REPL_REL_MULT	,	(fBETA2_REPL_MULT/fBETA2_D_MULT)*100	)
	vBETA2_RICH_REL_MULT	=	c(vBETA2_RICH_REL_MULT	,	(fBETA2_RICH_MULT/fBETA2_D_MULT)*100	)
	
	vBETA2_D			=	c(vBETA2_D			,	fBETA2_D			)
	vBETA2_REPL			=	c(vBETA2_REPL		,	fBETA2_REPL			)
	vBETA2_RICH			=	c(vBETA2_RICH		,	fBETA2_RICH			)
	vBETA2_REPL_REL		=	c(vBETA2_REPL_REL	,	fBETA2_REPL_REL		)
	vBETA2_RICH_REL		=	c(vBETA2_RICH_REL	,	fBETA2_RICH_REL		)
	
	vBETA2_REPL_SE		=	c(vBETA2_REPL_SE	,	fBETA2_REPL_SE		)
	vBETA2_RICH_SE		=	c(vBETA2_RICH_SE	,	fBETA2_RICH_SE		)
	vBETA2_REPL_REL_SE	=	c(vBETA2_REPL_REL_SE,	fBETA2_REPL_REL_SE	)
	vBETA2_RICH_REL_SE	=	c(vBETA2_RICH_REL_SE,	fBETA2_RICH_REL_SE	)
	
	vBETA2_TOTAL		=	c(vBETA2_TOTAL		,	fBETA2_TOTAL		)
	vBETA2_CONS_T		=	c(vBETA2_CONS_T		,	fBETA2_CONS_T		)
	vBETA2_CONS_REPL	=	c(vBETA2_CONS_REPL	,	fBETA2_CONS_REPL	)
	vBETA2_CONS_RICH	=	c(vBETA2_CONS_RICH	,	fBETA2_CONS_RICH	)
	
	vBETA2_MULT_CONS_D		=	c(vBETA2_MULT_CONS_D	,	fBETA2_MULT_CONS_D 		)
	vBETA2_MULT_CONS_REPL	=	c(vBETA2_MULT_CONS_REPL	,	fBETA2_MULT_CONS_REPL	)
	vBETA2_MULT_CONS_RICH	=	c(vBETA2_MULT_CONS_RICH	,	fBETA2_MULT_CONS_RICH	)
	

	return(
		list(ALPHA1				=	vALPHA1				,		
	         BETA1_MULT			= 	vBETA1_MULT			,
			 BETA1_ADD			= 	vBETA1_ADD			,
             GAMMA1				= 	vGAMMA1				,
			 BETA1_D_MULT		= 	vBETA1_D_MULT		,
             BETA1_REPL_MULT	= 	vBETA1_REPL_MULT	,
             BETA1_RICH_MULT	= 	vBETA1_RICH_MULT	,
			 BETA1_REPL_REL_MULT	= 	vBETA1_REPL_REL_MULT	,
             BETA1_RICH_REL_MULT	= 	vBETA1_RICH_REL_MULT	,
             BETA1_D			= 	vBETA1_D			,
             BETA1_REPL			= 	vBETA1_REPL			,
             BETA1_RICH			= 	vBETA1_RICH			,
             BETA1_REPL_REL		= 	vBETA1_REPL_REL		,
             BETA1_RICH_REL		= 	vBETA1_RICH_REL		,
             BETA1_REPL_SE		=	vBETA1_REPL_SE	    ,
			 BETA1_RICH_SE		=	vBETA1_RICH_SE	    ,
			 BETA1_REPL_REL_SE	=	vBETA1_REPL_REL_SE  ,
			 BETA1_RICH_REL_SE	=	vBETA1_RICH_REL_SE	,
			 BETA1_TOTAL		= 	vBETA1_TOTAL		,
             BETA1_CONS_T		= 	vBETA1_CONS_T		,
             BETA1_CONS_REPL	= 	vBETA1_CONS_REPL	,
             BETA1_CONS_RICH	= 	vBETA1_CONS_RICH	,
			 BETA1_MULT_CONS_D		= 	vBETA1_MULT_CONS_D		,
             BETA1_MULT_CONS_REPL	= 	vBETA1_MULT_CONS_REPL	,
             BETA1_MULT_CONS_RICH	= 	vBETA1_MULT_CONS_RICH	,
			 ALPHA2				=	vALPHA2				,		
	         BETA2_MULT			= 	vBETA2_MULT			,
			 BETA2_ADD			= 	vBETA2_ADD			,
             GAMMA2				= 	vGAMMA2				,
			 BETA2_D_MULT		= 	vBETA2_D_MULT		,
             BETA2_REPL_MULT	= 	vBETA2_REPL_MULT	,
             BETA2_RICH_MULT	= 	vBETA2_RICH_MULT	,
			 BETA2_REPL_REL_MULT	= 	vBETA2_REPL_REL_MULT	,
             BETA2_RICH_REL_MULT	= 	vBETA2_RICH_REL_MULT	,
             BETA2_D			= 	vBETA2_D			,
             BETA2_REPL			= 	vBETA2_REPL			,
             BETA2_RICH			= 	vBETA2_RICH			,
             BETA2_REPL_REL		= 	vBETA2_REPL_REL		,
             BETA2_RICH_REL		= 	vBETA2_RICH_REL		,
             BETA2_REPL_SE		=	vBETA2_REPL_SE	    ,
			 BETA2_RICH_SE		=	vBETA2_RICH_SE	    ,
			 BETA2_REPL_REL_SE	=	vBETA2_REPL_REL_SE  ,
			 BETA2_RICH_REL_SE	=	vBETA2_RICH_REL_SE	,		 
			 BETA2_TOTAL		= 	vBETA2_TOTAL		,
             BETA2_CONS_T		= 	vBETA2_CONS_T		,
             BETA2_CONS_REPL	= 	vBETA2_CONS_REPL	,
             BETA2_CONS_RICH	= 	vBETA2_CONS_RICH	,
			 BETA2_MULT_CONS_D		= 	vBETA2_MULT_CONS_D		,
             BETA2_MULT_CONS_REPL	= 	vBETA2_MULT_CONS_REPL	,
             BETA2_MULT_CONS_RICH	= 	vBETA2_MULT_CONS_RICH	
		)
	)

}

#function for pairwise comparison (convenience function)
pairwise_meanse_pvalue = 
  function(mean_v,se_v, g, alternative='two.sided', method = 'fdr')
  {
    n = length(unique(g))
    N = n*(n-1)/2
    d = data.frame(mean_v, se_v,  g = g)
    dfDATA_OUT=data.frame(matrix(NA,n,n))	
    colnames(dfDATA_OUT)=1:n
    row.names(dfDATA_OUT)=1:n  
    Z = data.frame(Sample1=rep("A", N),
                   Sample2=rep("A", N),
                   W=rep(NA, N),
                   p.value=rep(NA, N),
                   p_adj=rep(NA, N),
                   stringsAsFactors=FALSE)
    
    k=0               
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        k=k+1
        Namea = as.character(unique(g)[i])
        Nameb = as.character(unique(g)[j])
        Datax = subset(d, g==unique(g)[i])
        Datay = subset(d, g==unique(g)[j])
        Dataz = rbind(Datax, Datay)
        Dataz$g2 = factor(Dataz$g)
        #t-ratio
        t_ratio=(Datax$mean_v-Datay$mean_v)  / (Datax$se_v+Datay$se_v)
        pval=2*pnorm(-abs(t_ratio))	
        P = signif(pval, digits=4)
        S = signif(t_ratio, digits=4)
        P.adjust = NA                       
        Z[k,] =c(Namea, Nameb, 
                 S, P, P.adjust)
      }
    } 
    Z$p_adj = signif(p.adjust(Z$p.value, method = method), digits=4) 
    
    #populate a distance matrix
    for(nROW in 1:nrow(Z)){
      dfDATA_OUT[as.numeric(Z[nROW,1]),as.numeric(Z[nROW,2])]=Z[nROW,5]
      dfDATA_OUT[as.numeric(Z[nROW,2]),as.numeric(Z[nROW,1])]=Z[nROW,5]
    }
    
    return(list(long=Z,p_adj_dist=dfDATA_OUT))
  }

