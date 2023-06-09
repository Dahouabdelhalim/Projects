####################################################################
##
## KittiwakePartitioningNewLuckPluck13: Partitions luck and pluck and
## breeding status according to age for the kittiwake model.
##
#######################################################################

rm(list=ls(all=TRUE)); graphics.off(); 
require(Matrix); require(statmod);

## Change path as needed.
path ="/home/robin/ipm/drivers/Manuscripts/Luckpartition/code/"
out=try(setwd(path),silent=TRUE);

if(class(out)=="try-error") { 
    path=ifelse(.Platform$OS.type=="windows","c:/repos/drivers/Manuscripts/luckPartition","~/repos/drivers/Manuscripts/LuckPartition"); 
    setwd(path); 
  }
getwd(); 

source ("KittiwakeMatrices.R");
## Save the original fecundity matrix as origF.
origF = F;
source("Standard Graphical Pars.R")
source("Utilities.R") 
source("MegamatrixFunctions.R");

mz = 25;   # size of the age/stage transition matrix
bigmz = 2*mz  ## size of the expanded state space for kids/no-kids

### maxAge is really maximum age + 1, because age starts at 0 
maxAge = 100 # good for plotting. But it takes maxAge=300 for (sum over ages) \\approx total.

## Prob. of breeding (is 1 for breeding classes and zero otherwise
## because the real probability of breeding is given by the transition
## rates into the breeding classes.)
pb = c(0, 0, 0, 1, 1) # prob. of breeding
pbMegavec = rep(pb, 5)
bigpbMegavec = rep(pbMegavec, 2)  ## expanded state space

## Variance in per capita # of recruits. 
## Assumes those with 2 or 3 recruits have equal prob of 2 or 3. 
sigbsq = c(0, 0, 0, 0, 0.25)
sigbsqMegavec = rep(c(0, 0, 0, 0, 0.25), 5)
bigsigbsqMegavec = rep(sigbsqMegavec, 2)

## Make cross-classified vectors
## expected per capita # of recruits
bMegavec = rep(b, 5)
bigbMegavec = rep(bMegavec, 2)

## Set up birth size distribution - all are immature, age 1.  
c0 = c(1, rep(0, mz-1));
bigc0 = c(1, rep(0, 2*mz-1))

## How finely grained is the trait variation?
numTraitVals = 51

## function to calculate unconditional P ######
makeP = function (newpars, breedingProbFactor) {
  
  ## Modify reproductive stage transition matrices
  newM3 = M3; newM4 = M4; newM5 = M5; 
  
  newM5[4:5,] = newM5[4:5,] * breedingProbFactor
  for (j in 1:5) {
    B = sum(newM5[4:5,j]); A = sum(newM5[1:3,j]); bfac = (1-B)/A; 
    newM5[1:3,j]=newM5[1:3,j]*bfac; 
  }
  
  s = c(1, exp(newpars)/(1 + exp(newpars)))
  P = matrix (0, 25, 25)
  P[5*1 + 1:5, 5*0 + 1:5] = M1 * matrix (rep(s1, 5), 5, 5, byrow=TRUE)
  P[5*2 + 1:5, 5*1 + 1:5] = M2 * matrix (rep(s2, 5), 5, 5, byrow=TRUE)
  P[5*3 + 1:5, 5*2 + 1:5] = newM3 * matrix (rep(s3, 5), 5, 5, byrow=TRUE)
  P[5*4 + 1:5, 5*3 + 1:5] = newM4 * matrix (rep(s, 5), 5, 5, byrow=TRUE)
  P[5*4 + 1:5, 5*4 + 1:5] = newM5 * matrix (rep(s, 5), 5, 5, byrow=TRUE)

  if (sum(P < 0)) stop
  (c("Negative entries in P, breedingProbFactor = ",
  breedingProbFactor, "\\n"))
  return (P)
}

## function to calculate P conditional on breeding ###################
## Note: this is the big P in eq. S15
makePCondBreed = function (newpars, breedingProbFactor) {

  P = makeP (newpars, breedingProbFactor)
  bigF = bigP = P0 = matrix (0, bigmz, bigmz)

  indicatorS = rep(c(1, 1, 1, 0, 0), 5)
  indicatorM = rep(c(0, 0, 0, 1, 1), 5)

  ## make F for the extended state space Z1 \\cup Z2.  Note that the
  ## only non-zero entries are from Z2 to Z1.

  bigF[1, (mz+1):bigmz] = origF[1,]

  ## Make P for the extended state space Z1 \\cup Z2
  ## Z1 to Z1
  for (j in 1:mz) 
    bigP[j, 1:mz] = P[j,]*indicatorS[j]
  ## Z1 to Z2
  for (j in 1:mz)
    bigP[mz+j, 1:mz] = P[j,]*indicatorM[j]
  ## Z2 to Z2
  for (j in 1:mz)
    bigP[mz+j, (mz+1):(2*mz)] = P[j,]

  ## Make P0, as in KittwakePartitioningNewLuckPluck9.R.
  for (j in 1:bigmz)
    P0[j,] = (1 - bigpbMegavec)*bigP[j,]

  N0 = solve (diag(bigmz) - P0)
  probEverBreed = bigpbMegavec %*% N0

  PCondBreed = matrix (NA, bigmz, bigmz)
  for (j in 1:bigmz) {
    PCondBreed[,j] = bigP[,j]*probEverBreed/probEverBreed[j]
  }

  out = list (PCondBreed=PCondBreed,
              probEverBreed=probEverBreed,
              bigF=bigF)
  return (out)
}

calcLuckPluck = function () {

  expZRhoZeroACondXZ = matrix(0,maxAge,numTraitVals);
  expZStateTrajecLuckAgeCondXZ = expZFecLuckAgeCondXZ =
    matrix (0, maxAge, numTraitVals)
  expRCondXZ = rep (NA, mz)

  ## We need Fbar and Pbar for pluck.  Average kernels over trait
  ## conditional on breeding.  For kittiwakes, F does not depend on x,
  ## so Fbar = F.
  Pbar = matrix (0, mz, mz)

  for (ii in 1:numTraitVals) {
    ## Set survival qualities, breeding probs.
    newpars = survQualityDist[ii,]
    breedingProbFactor = breedingProbFactorDist[ii]

    P = makeP (newpars, breedingProbFactor)
    
    Pbar = Pbar + survQualityProb[ii] * P
  }
  Fbar = F  ## no trait dependence in F
  Nbar = solve (diag(mz) - Pbar)

  for (ii in 1:numTraitVals) {

    ## Set survival qualities, breeding probs.
    newpars = survQualityDist[ii,]
    breedingProbFactor = breedingProbFactorDist[ii]

    ## Make P conditional on breeding and never breeding
    P = makeP (newpars, breedingProbFactor)

    ## Fundamental matrix (N-tilde in the ms)
    N = solve (diag(mz) - P)

    ## Now that we have P, compute Pplus
    Pplus = cbind (P, rep(0, mz))
    deathProb = 1 - apply (Pplus, 2, sum)
    Pplus = rbind (Pplus, deathProb)
    
    ## expected lifetime reproduction conditional on state (rho1)
    expRCondXZ = apply (F %*% N, 2, sum)
    ## rho1 includes RO for the state of being dead
    rho1 = c(expRCondXZ, 0)
    
    ## V (eq. 7)
    V = rho1^2 %*% Pplus - (rho1 %*% Pplus)^2
    ## All we really need is V-tilde.  Chop off the last entry.
    V = V[-(mz+1)]
      
    ## There is no birth luck since everyone is born the same size.
    ## Var_0 (\\Ex [R | x, z_0]) = 0.
      
    ## time to loop over age
    ## record luck and pluck terms as functions of age and x
    Pa = diag(1, mz)
    PaC0 = c0

    stateTrajecLuckAgeCondXZ = V %*% Pa
    expZStateTrajecLuckAgeCondXZ[1,ii] = stateTrajecLuckAgeCondXZ %*%
      c0

    fecLuckAgeCondXZ = sigbsqMegavec %*% Pa
    expZFecLuckAgeCondXZ[1,ii] = fecLuckAgeCondXZ %*% c0

    eT = matrix(1,1,mz); 
    Ax = eT %*% F %*% N;
    bT = eT %*% (Fbar %*% Nbar - F %*% N) %*% P;
    rhoZeroACondXZ = eT%*%(F + Fbar%*%Nbar%*%P);
    expZRhoZeroACondXZ[1,ii] = rhoZeroACondXZ %*% c0
    
    ## a is actually age + 1, since the first year of life is age 0
    for (a in 2:maxAge) {
      Pa = P %*% Pa
      PaC0 = Pa %*% c0
      
      stateTrajecLuckAgeCondXZ = V %*% Pa
      expZStateTrajecLuckAgeCondXZ[a,ii] =
        stateTrajecLuckAgeCondXZ %*% c0
      
      fecLuckAgeCondXZ = sigbsqMegavec %*% Pa
      expZFecLuckAgeCondXZ[a,ii] = fecLuckAgeCondXZ %*% c0

      rhoZeroACondXZ = Ax + bT %*% Pa;
      expZRhoZeroACondXZ[a,ii] = rhoZeroACondXZ %*% c0

    }  ## end loop over age

    cat ("ii = ", ii, "\\n")
    
  }  ## end loop over trait values

  ## Perform appropriate averages over trait. If variation in survival
  ## and breeding probability is perfectly correlated, then P(surv,
  ## breedingProb) = P(surv).
  
  ## NEW LUCK ##################################################

  traitAve = function (x) {
    return (sum(survQualityProb*x))
  }

  traitVar = function (x) {
    return (sum(survQualityProb*x^2) - sum(survQualityProb*x)^2)
  }

  ## Average over x: Ex_{x,B} (Var(R | x, B))
  stateTrajecLuckAge = expXZStateTrajecLuckAgeCondXZ =
    apply (expZStateTrajecLuckAgeCondXZ, 1, traitAve)
  fecLuckAge = expXZFecLuckCondAgeXZ =
    apply (expZFecLuckAgeCondXZ, 1, traitAve)

  ## Take variance over x 
  varXExpZRhoZeroA = apply (expZRhoZeroACondXZ, 1, traitVar)
  ## Find marginal increase over age
  pluckAge = diff (c(0, varXExpZRhoZeroA))

##############################################################

  out = list (
      stateTrajecLuckAge=stateTrajecLuckAge,
      fecLuckAge=fecLuckAge,
      pluckAge=pluckAge)

  return (out)
}

calcLuckPluckCondBreeding = function () {

  expZRhoZeroACondXBZ = matrix(0,maxAge,numTraitVals);
  expZBStateTrajecLuckAgeCondXBZ = expZBFecLuckAgeCondXBZ =
    matrix (0, maxAge, numTraitVals)
  expRCondXBZ = rep (NA, bigmz)
  pXCondB = rep(0, numTraitVals)

  ## We need Fbar and Pbar for pluck.  Average kernels over trait
  ## conditional on breeding.  For kittiwakes, F does not depend on x,
  ## so Fbar = F.
  Pbar = matrix (0, bigmz, bigmz)

  pBCondZ = 0
  for (ii in 1:numTraitVals) {
    ## Set survival qualities, breeding probs.
    newpars = survQualityDist[ii,]
    breedingProbFactor = breedingProbFactorDist[ii]

    out = makePCondBreed (newpars, breedingProbFactor)
    P = out$PCondBreed
    F = out$bigF
    pBCondXZ = out$probEverBreed
    
    Pbar = Pbar + survQualityProb[ii] * P
    ## P(M=1 | z) = \\sum_x P(M=1 | x, z) * P(x)
    pBCondZ = pBCondZ + pBCondXZ*survQualityProb[ii]
  }
  Fbar = F  ## no trait dependence in F
  Nbar = solve (diag(bigmz) - Pbar)
  ## P(M = 1) = \\sum_z P(M=1 | z) * P(z)
  pB = sum(pBCondZ * bigc0)

  for (ii in 1:numTraitVals) {

    ## Set survival qualities, breeding probs.
    newpars = survQualityDist[ii,]
    breedingProbFactor = breedingProbFactorDist[ii]

    ## Make P conditional on breeding and never breeding
    out = makePCondBreed (newpars, breedingProbFactor)
    P = out$PCondBreed
    F = out$bigF
    ## P(M=1 | x, z)
    pBCondXZ = out$probEverBreed
    ## P(x | M=1, z) = P(M=1 | x, z) * P(x) / P(M=1 | z)
    ## pXCondBZ = pBCondXZ * survQualityProb / pBCondZ
    ## P(M=1 | x) = \\sum_z P(M=1 | x, z) * P(z)
    pBCondX = sum(pBCondXZ * bigc0)
    ## P(x | M=1) = P(M=1 | x) * P(x) / P(M=1)
    pXCondB[ii] = pBCondX * survQualityProb[ii] / pB

    ## Find c0 conditional on breeding
    bigc0CondBreed = bigc0 * pBCondXZ / pBCondX

    ## Fundamental matrix (N-tilde in the ms)
    N = solve (diag(bigmz) - P)

    ## Now that we have P, compute Pplus
    Pplus = cbind (P, rep(0, bigmz))
    deathProb = 1 - apply (Pplus, 2, sum)
    Pplus = rbind (Pplus, deathProb)
    
    ## expected lifetime reproduction conditional on state (rho1)
    expRCondXBZ = apply (F %*% N, 2, sum)
    ## rho1 includes RO for the state of being dead
    rho1 = c(expRCondXBZ, 0)
    
    ## V (eq. 7)
    V = rho1^2 %*% Pplus - (rho1 %*% Pplus)^2
    ## All we really need is V-tilde.  Chop off the last entry.
    V = V[-(bigmz+1)]
      
    ## There is no birth luck since everyone is born the same size.
    ## Var_0 (\\Ex [R | x, z_0]) = 0.
      
    ## time to loop over age
    ## record luck and pluck terms as functions of age and x
    Pa = diag(1, bigmz)
    PaC0 = bigc0

    stateTrajecLuckAgeCondXBZ = V %*% Pa
    expBStateTrajecLuckAgeCondXBZ = pBCondXZ*stateTrajecLuckAgeCondXBZ
    expZBStateTrajecLuckAgeCondXBZ[1,ii] = expBStateTrajecLuckAgeCondXBZ %*%
      bigc0

    fecLuckAgeCondXBZ = bigsigbsqMegavec %*% Pa
    expBFecLuckAgeCondXBZ = pBCondXZ*fecLuckAgeCondXBZ
    expZBFecLuckAgeCondXBZ[1,ii] = expBFecLuckAgeCondXBZ %*% bigc0

    eT = matrix(1,1,bigmz); 
    Ax = eT %*% F %*% N;
    bT = eT %*% (Fbar %*% Nbar - F %*% N) %*% P;
    rhoZeroACondXBZ = eT%*%(F + Fbar%*%Nbar%*%P);
    expZRhoZeroACondXBZ[1,ii] = rhoZeroACondXBZ %*% t(bigc0CondBreed)
    
    ## a is actually age + 1, since the first year of life is age 0
    for (a in 2:maxAge) {
      Pa = P %*% Pa
      PaC0 = Pa %*% bigc0
      
      stateTrajecLuckAgeCondXBZ = V %*% Pa
      expBStateTrajecLuckAgeCondXBZ = pBCondXZ*stateTrajecLuckAgeCondXBZ
      expZBStateTrajecLuckAgeCondXBZ[a,ii] =
        expBStateTrajecLuckAgeCondXBZ %*% bigc0
      
      fecLuckAgeCondXBZ = bigsigbsqMegavec %*% Pa
      expBFecLuckAgeCondXBZ = pBCondXZ*fecLuckAgeCondXBZ
      expZBFecLuckAgeCondXBZ[a,ii] = expBFecLuckAgeCondXBZ %*% bigc0

      rhoZeroACondXBZ = Ax + bT %*% Pa;
      expZRhoZeroACondXBZ[a,ii] = rhoZeroACondXBZ %*% t(bigc0CondBreed)

    }  ## end loop over age

    cat ("ii = ", ii, "\\n")
    
  }  ## end loop over trait values

  ## Perform appropriate averages over trait. If variation in survival
  ## and breeding probability is perfectly correlated, then P(surv,
  ## breedingProb) = P(surv).
  
  ## NEW LUCK ##################################################

  traitAve = function (x) {
    return (sum(survQualityProb*x))
  }

  traitVar = function (x) {
    return (sum(survQualityProb*x^2) - sum(survQualityProb*x)^2)
  }

  traitVarCondBreed = function (x) {
    return (sum(pXCondB*x^2) - sum(pXCondB*x)^2)
  }

  ## Average over x: Ex_{x,B} (Var(R | x, B))
  stateTrajecLuckAgeWithinBreeding = expXZBStateTrajecLuckAgeCondXBZ =
    apply (expZBStateTrajecLuckAgeCondXBZ, 1, traitAve)
  fecLuckAgeWithinBreeding = expXZBFecLuckCondAgeXBZ =
    apply (expZBFecLuckAgeCondXBZ, 1, traitAve)

  ## Take variance over x conditional on breeding: Var_{x|B=1}
  varXExpZRhoZeroA = apply (expZRhoZeroACondXBZ, 1, traitVarCondBreed)
  ## Average over B
  expBVarXExpZRhoZeroA = pB * varXExpZRhoZeroA
  pluckAgeWithinBreeding = diff (c(0, expBVarXExpZRhoZeroA))

##############################################################

  out = list (
      stateTrajecLuckAgeWithinBreeding=stateTrajecLuckAgeWithinBreeding,
      fecLuckAgeWithinBreeding=fecLuckAgeWithinBreeding,
      pluckAgeWithinBreeding=pluckAgeWithinBreeding)

  return (out)
}

## estimate CV of adult survivals from Cam et al. 2002, Fig. 1
q75 = log(0.9/0.1);   #75th percentile on inverse-logit scale ("quality")
q25 = log(0.7/0.3)    #25th percentile on inverse-logit scale ("quality")
IQR=q75-q25;
sigma=IQR/1.35;
mu = (q75+q25)/2;
CVSurv = sigma/mu;

## What is the estimated real CV for transitions to breeding classes?
## estimate CV of breeding prob. from Cam et al. 2002
q25 = 0.94
q75 = 0.98
IQR=q75-q25;
sigma=IQR/1.35;
mu = (q75+q25)/2;
CVBreedingProb = sigma/mu;

## Make matrix of inverse logit survivals (survival "quality")
sdSurvQuality = CVSurv * pars
minSurvQuality = pars - 2.5*sdSurvQuality 
maxSurvQuality = pars + 2.5*sdSurvQuality 
survQualityDist = matrix (0, numTraitVals, 4)
for (i in 1:4) {
  survQualityDist[,i] = seq(minSurvQuality[i],
                            maxSurvQuality[i], length.out=numTraitVals)
}
## vector of associated probabilities.

## Why are we using pars[1] as the mean and sdSurvQuality[1] as the
## sd?  Those are inverse logit (survival) mean and sd for stage 2,
## ages 4 and up.  Why stage 2?  I guess we have to choose
## something...
survQualityProb = dnorm(survQualityDist[,1], mean=pars[1],
                        sd=sdSurvQuality[1])
survQualityProb = survQualityProb/sum(survQualityProb)

## Now work on breeding prob.

## This is the range of breedingProbFactor that we need to get a +/-
## 1.5 std dev. variation in breeding probability (the "real"
## std. dev.) when we're using the "real" CV of breeding prob.
minBreedingProbFactor = 1 - 25.92*CVBreedingProb   
maxBreedingProbFactor = 1 + 25.92*CVBreedingProb   

## Breeding prob. adjustment factors 
breedingProbFactorDist =
  seq(minBreedingProbFactor, maxBreedingProbFactor, length.out=numTraitVals)

## why 18*sigma?  We want breeding probability to vary over an
## interval of 3*sigma=0.0888, which corresponds to varying the
## breeding probability factor over an interval of 1.6.  If we want
## the breeding probability factor interval to correspond to 3 times
## the standard deviation of the breeding probability factor
## distribution, then it seems like the sd of the breeding probability
## factor should be 1.6/0.0888 * sigma.
breedingProbFactorProb = dnorm (breedingProbFactorDist, mean=1,
                                sd=18*sigma)
breedingProbFactorProb = breedingProbFactorProb/sum(breedingProbFactorProb) 

######################################################################
## Now that trait variation is set up, start calculating
## partitionings.
######################################################################

## Calculate age partitionings of luck and pluck without conditioning
## on breeding.
out = calcLuckPluck ()
NewPluckAge = out$pluckAge
fecLuckAge = out$fecLuckAge
stateTrajecLuckAge = out$stateTrajecLuckAge
NewLuckAge = stateTrajecLuckAge + fecLuckAge

## Calculate within- and between-breeding age partitionings of luck
## and pluck.
out = calcLuckPluckCondBreeding ()
stateTrajecLuckAgeWithinBreeding = out$stateTrajecLuckAgeWithinBreeding
fecLuckAgeWithinBreeding = out$fecLuckAgeWithinBreeding
pluckAgeWithinBreeding = out$pluckAgeWithinBreeding
stateTrajecLuckAgeBetweenBreeding = stateTrajecLuckAge -
  stateTrajecLuckAgeWithinBreeding
fecLuckAgeBetweenBreeding = fecLuckAge -
  fecLuckAgeWithinBreeding
pluckAgeBetweenBreeding = NewPluckAge -
  pluckAgeWithinBreeding

###############################################################
## Make plots.
###############################################################
require(viridis); cols=plasma(100);

graphics.off(); dev.new(width=10,height=10);
par(fig=c(0,1,0,1),yaxs="i",bty="l",xaxs="i",mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4,mar=c(4,4,3,2),mfrow=c(2,2)) 

add_panel_label <- function(ltype="a") {
    text <- paste(LETTERS[letters==ltype], ")", sep="")
    mtext(text=text, side=3, adj=0)
}

## Plot a
plot (0:(maxAge-1), NewLuckAge[1:maxAge], type="l", xlab="Age", ylab="Contribution to Var(LRO)",lwd=2)
lines (0:(maxAge-1), NewPluckAge[1:maxAge], col="red",lwd=2)
lines (0:(maxAge-1), stateTrajecLuckAge[1:maxAge], col="black",lwd=2, lty=2)
lines (0:(maxAge-1), fecLuckAge[1:maxAge], col="black",lwd=2, lty=3)
legend ("topright",inset=0,col=c("black", "red", "black", "black"),
        lty=c(1,1,2,3), bty="n", lwd=2,cex=1.5,
        legend=c("Luck", "Quality", "State trajec. luck", "Fec. luck"))
add_panel_label ("a")

## Plot b
## Luck and Pluck are both 0 at ages 0 and 1.  Define the ratio to be 0.
plot (0:(maxAge-1),
      c(0, 0, (NewPluckAge/(NewLuckAge + NewPluckAge))[3:maxAge]), type="l", xlab="Age", ylab="Fraction contribution from quality", lwd=2)
add_panel_label ("b")

## Plot c
## What are the mean, 10%, and 90% bins for quality?
cumSurvQualityProb = rep (0, numTraitVals)
cumSurvQualityProb[1] = survQualityProb[1]
for (ii in 2:numTraitVals)
  cumSurvQualityProb[ii] = cumSurvQualityProb[ii-1] + survQualityProb[ii]
## Which is the bin closest to the 10th %ile?, 90th?
index10 = which.min (abs(cumSurvQualityProb - 0.1))
index90 = which.min (abs(cumSurvQualityProb - 0.9))
index50 = which.min (abs(cumSurvQualityProb - 0.5))
survProb = matrix (0, maxAge, 3)
survProb[1,] = 1
lowP = makeP (newpars=survQualityDist[index10],
              breedingProbFactorDist[index10])
medianP = makeP (newpars=survQualityDist[index50],
                 breedingProbFactorDist[index50])
highP = makeP (newpars=survQualityDist[index90],
               breedingProbFactorDist[index90])
lowPa = medianPa = highPa = diag(mz)
for (a in 2:maxAge) {
  lowPa = lowPa %*% lowP
  medianPa = medianPa %*% medianP
  highPa = highPa %*% highP
  survProb[a,1] = apply (lowPa, 2, sum) %*% c0
  survProb[a,2] = apply (medianPa, 2, sum) %*% c0
  survProb[a,3] = apply (highPa, 2, sum) %*% c0
}
for (k in 1:nrow(survProb)) survProb[k,] =
                              survProb[k,]/sum(survProb[k,])

plot (0:(maxAge-1), rep(1, maxAge), xlab="Age", ylim=c(0,1),
      ylab="Fraction among survivors", type="n", lwd=2, col="grey")
    lines (0:(maxAge-1), survProb[,1], lwd=2, col="grey")
    lines (0:(maxAge-1), survProb[,1]+survProb[,2], lwd=2, col="grey")
    polygon(c(0:(maxAge-1), (maxAge-1):0),
            c(survProb[,1], rep(0, maxAge)), col="red", border=NA)
    polygon(c(0:(maxAge-1), (maxAge-1):0),
            c(survProb[,1]+survProb[,2], rev(survProb[,1])),
            col="orange", border=NA) 
    polygon(c(0:(maxAge-1), (maxAge-1):0),
            c(rep(1, maxAge), rev(survProb[,1]+survProb[,2])),
            col="yellow", border=NA) 
  text (x=7, y=0.1, labels="10%ile", col="black", cex=1.4)
  text (x=17, y=0.2, labels="50%ile", col="black", cex=1.4)
  text (x=50, y=0.6, labels="90%ile quality", col="black", cex=1.4)
add_panel_label ("c")

## plot d
plot (0:(maxAge-1), stateTrajecLuckAge[1:maxAge], lwd=2, type="l",
      xlab="Age", ylab="Contributions to Var(R)")
lines (0:(maxAge-1), stateTrajecLuckAgeBetweenBreeding[1:maxAge], lwd=2,
       lty=2)
lines (0:(maxAge-1), NewPluckAge[1:maxAge], lwd=2, col="red")
lines (0:(maxAge-1), pluckAgeBetweenBreeding[1:maxAge], lwd=2,
      lty=2, col="red")
legend ("topright", lty=rep(1:2, 2), col=c(rep("black", 2), rep("red",
      2)), lwd=2, bty="n",
        legend=c("Total state trajectory luck", "Breed yes/no",
      "Total quality", "Breed yes/no"), cex=1.5)
add_panel_label ("d")

##dev.copy2pdf(file="KittiwakePartitioningNewLuckPluck13.pdf")

