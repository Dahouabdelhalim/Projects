################################################################################
## Filename: bcispeciescalcs.r
## Authors: Helene C. Muller-Landau and Marco D. Visser
## Article title: How do lianas and vines influence competitive differences 
##                and niche differences among tree species?
##                Concepts and a case study in a tropical forest 
## Purpose of this code: Calculation of species-level liana-tree interaction 
##         statistics for BCI tree species 
################################################################################

rm(list=ls()) # clear workspace

################################################################################
# directory name and file names for input files 
infiledir <- "./objects/"
bcispecdatafn <- paste0(infiledir,"BCIspeciesinputdata.csv")
lambdabsfn <- paste0(infiledir,"BootstrappedLambdas.csv")

# directory name and file names for output files
outfiledir <- "./objects/"
statlianacsvfn <- paste0(outfiledir,"treespstatsrelianas.csv")
statbsfn <- paste0(outfiledir,"bootstraplianastats.rds")


################################################################################
# DEFINE CONSTANTS
ALLCOVERS <- c(0,1/8,3/8,5/8,7/8) # proportional liana covers for classes 0 to 4

ALPHA <- 0.05 # p-values for confidence intervals
PLO <- ALPHA/2
PHI <- 1-ALPHA/2


################################################################################
# function to get key statistics from bootstrapped parameter value datasets
# varbs is a matrix with the species in rows, and bootstrapped replicates in columns 
################################################################################
getstatsfrombs <- function(varbs,varname) {
  thisdf <- data.frame(SD=apply(varbs,1,sd,na.rm=TRUE),
                    LoCI=apply(varbs,1,quantile,probs=PLO,na.rm=TRUE),
                    HiCI=apply(varbs,1,quantile,probs=PHI,na.rm=TRUE))
  names(thisdf) <- paste0(varname,names(thisdf))
  return(thisdf)
} # end getstatsfrombs



################################################################################
# function to calculate liana load, lambdaTot, lambdaL, lambda4-lambdaTot, 
# burden (loglambdaTot-loglambda0), and tolerance  
################################################################################
# Nbycover is a matrix with rows for species and columns for the number of individuals by liana cover class
# lambdas is a matrix with rows for species and columns for lambdas for different liana cover classes
calclianastats <- function(Nbycover,lambdas) {
  loglambdas <- log(lambdas)
  nsp <- nrow(Nbycover)
  nlc <- ncol(Nbycover) # number of liana cover classes
  Ntot <- apply(Nbycover,1,sum,na.rm=TRUE)
  Ninf <- apply(Nbycover[,2:(nlc)],1,sum,na.rm=TRUE)
  ps <- Nbycover/Ntot
  outstats <- data.frame(prevalence=Ninf/Ntot)
  outstats$load <- apply(Nbycover[,2:nlc],1,function(X) weighted.mean(ALLCOVERS[2:nlc],w=X))
  outstats$tolerance <- apply(loglambdas,1,cov,y=ALLCOVERS)/var(ALLCOVERS) # formula for linear regression slope
  lambdaTot <- apply(ps*lambdas,1,sum)  
  outstats$loglambdaTot <- log(lambdaTot)
  outstats$burden <- outstats$loglambdaTot-loglambdas[,1]
  outstats$loglambdadifC4 <- outstats$loglambdaTot-loglambdas[,nlc]
  outstats$loglambda0 <- loglambdas[,1] 
  outstats$loglambda4 <- loglambdas[,nlc]
  return(outstats)
} # end calclianastats

################################################################################

## LOAD FILES 

# load species level data on liana cover, eigenvalues, and shadetolerance 
totDat <- read.csv(bcispecdatafn)

# sets of bootstrapped eigenvalues (note these are correlated across species and liana cover classes)
# from  Visser et al. 2018. Journal Of Ecology 106: 781-794.
lambdabs2d <- read.csv(lambdabsfn)
bootstrappedLambdas <- array(0,dim=c(33,5,100))
for (i in 1:5) {
  bootstrappedLambdas[,i,] <- t(lambdabs2d[,(33*i-32):(33*i)])  
}
dimnames(bootstrappedLambdas) <- list(substr(dimnames(lambdabs2d)[[2]][1:33],3,8),
                        paste0("lambda",seq(0,4)),seq(1,100))


################################################################################

# do calculations on main datasets
totDat$Ninf  <-  totDat$N1+totDat$N2+totDat$N3+totDat$N4
totDat$Ntot  <-  totDat$N0+totDat$Ninf
totDat$loglambda0 <- log(totDat$lambda0)
totDat$loglambda1 <- log(totDat$lambda1)
totDat$loglambda2 <- log(totDat$lambda2)
totDat$loglambda3 <- log(totDat$lambda3)
totDat$loglambda4 <- log(totDat$lambda4)

# use previous 100 bootstraps of the lambdas to obtain the CIs on lambda0, loglambda0, etc.
for (i in 1:5) {
  thislogstats <- getstatsfrombs(log(bootstrappedLambdas[,i,]),paste0("loglambda",i-1))
  totDat <- cbind(totDat,thislogstats)
}

lianastats <- calclianastats(Nbycover=as.matrix(totDat[,c("N0","N1","N2","N3","N4")]),
                          lambdas=as.matrix(totDat[,c("lambda0","lambda1","lambda2","lambda3","lambda4")]))
totDat <- cbind(totDat,lianastats[,is.na(match(names(lianastats),names(totDat)))])


################################################################################
## Now independently bootstrap liana cover distributions and the lambda sets
# to obtain confidence intervals on the load, tolerance, loglambdaTot, 
# burden <- loglambdaTot-loglambda0, and loglambdadif4C <- loglambda4-loglambdaTot
nstat <- ncol(lianastats)
Nboots  <-  1000
neigensets <- dim(bootstrappedLambdas)[3]
nsp <- nrow(totDat)
statbs <- array(NA, dim=c(nsp,nstat,Nboots),
             dimnames=list(sp=totDat$sp,stats=names(lianastats),
                           bootstraps=seq(1,Nboots)))
truewts <- t(totDat[,c("N0","N1","N2","N3","N4")]/totDat$Ntot)
for(i in 1:Nboots){
    # resample the distribution of liana covers for each species
    thisNbycover  <-  t(sapply(1:nsp,function(X) 
      rmultinom(1,totDat$Ntot[X],truewts[,X])))

    # sample one of the sets of eigenvalues in the bootstrapped set 
    thisLambdaSp <- bootstrappedLambdas[,,sample(neigensets,1)]
    
    # calculate statistics 
    thislianastats <- calclianastats(thisNbycover,thisLambdaSp)     
    
    statbs[,,i] <- as.matrix(thislianastats)
    cat("\\r",round(i/Nboots,3)*100,"% \\r") # this takes a while to run - let the user know something is happening
}
saveRDS(statbs,file=statbsfn)
cat("\\n")
newstatbs  <-  statbs[,is.na(match(paste0(dimnames(statbs)[[2]],"SD"),names(totDat))),]
for (i in 1:dim(newstatbs)[2]) {
  thesestats <- getstatsfrombs(newstatbs[,i,],dimnames(newstatbs)[[2]][i])
  totDat <- cbind(totDat,thesestats)
}
########################## end bootstrap calculations ######################################

write.csv(totDat,file=statlianacsvfn,row.names=FALSE)








