##############################################################################
##
## Sneetches are classified as Small, Medium, Large, and Xtra-large. 
## Small and Medium sneeches do not breed (p_b =0). 
## Large and Xtra-large sneetches have p_b>0, and Poisson offspring number
## conditional on attempting to breed. 
##
## Note that this code uses B = 1, 0 as an indicator variable for
## whether or not an individual has bred, where the paper uses M.  pB
## = Pr(M = 1).
## 
##############################################################################

rm(list=ls(all=TRUE)); graphics.off();

############ CHANGE AS NEEDED to set working directory 
root ="/home/robin/ipm/drivers"; out=try(setwd(root),silent=TRUE);
if(class(out)=="try-error") { 
    root=ifelse(.Platform$OS.type=="windows","c:/repos/drivers","~/ipm/drivers"); 
    setwd(root); 
  }
setwd("Manuscripts/Luckpartition/code");

## What life history do you want?
fastDevel = FALSE
slowDevel = TRUE

shortRepro = FALSE
longRepro = TRUE

highJuvMort = FALSE
lowJuvMort = TRUE

### Parameter values, and functions to make P and conditional-P matrices 
source("SneetchMatrices.R")
## Helpful graphics stuff
source("Standard Graphical Pars.R")

calcLuckPluckCondBreeding = function (bigmz, Plist, F, pBCondXZMat, bigpb, bigb,
                                      bigsigbsq, bigc0) {
  expZRhoZeroACondXBZ = matrix(0,maxAge,numTraitVals);
  stateTrajecLuckAge = fecLuckAge = rep (0, maxAge)
  expZBStateTrajecLuckAgeCondXBZ = expZBFecLuckAgeCondXBZ =
    matrix (0, maxAge, numTraitVals)
  pXCondB = rep(0, numTraitVals)
  
  expRCondXBZ = varRCondXBZ = rep(NA, bigmz)
  expZBVarRCondXBZ = expZVarBExpRCondXBZ = numeric(numTraitVals)

  ## We need Fbar and Pbar for pluck.  Average kernels over trait
  ## conditional on breeding.  
  Pbar = matrix (0, bigmz, bigmz)
  pBCondZ = 0
  for (ii in 1:numTraitVals) {
    P = Plist[[ii]]
    Pbar = Pbar + starDist[ii] * P
  }
  ## For sneetches, F does not depend on x, so Fbar = F.
  Fbar = F  
  Nbar = solve (diag(bigmz) - Pbar)

  ## P(M=1 | z) = \\sum_x P(M=1 | x, z) * P(x)
  pBCondZ = starDist %*% pBCondXZMat
  ## P(M = 1) = \\sum_z P(M=1 | z) * P(z)
  pB = sum(pBCondZ * bigc0)

  for (ii in 1:numTraitVals) {
    ## Transition matrix for this trait value
    P = Plist[[ii]]
    ## Fundamental matrix (N-tilde in the ms)
    N = solve (diag(bigmz) - P)

    pBCondXZ = pBCondXZMat[ii,]

    ## P(M=1 | x) = \\sum_z P(M=1 | x, z) * P(z)
    pBCondX = sum(pBCondXZ * bigc0)
    ## P(x | M=1) = P(M=1 | x) * P(x) / P(M=1)
    pXCondB[ii] = pBCondX * starDist[ii] / pB

    ## Find c0 conditional on breeding
    bigc0CondBreed = bigc0 * pBCondXZ / pBCondX

    ## Now that we have P, compute Pplus, defined in "Background and
    ## Assumptions." 
    Pplus = cbind (P, rep(0, bigmz))
    deathProb = 1 - apply (Pplus, 2, sum)
    Pplus = rbind (Pplus, deathProb)
    
    ## expected lifetime reproduction conditional on state (rho1)
    expRCondXBZ = apply (F %*% N, 2, sum)
    ## rho1 includes RO for the state of being dead
    rho1 = c(expRCondXBZ, 0)

    ## V = variance of rho1 (eq. 5)
    V = rho1^2 %*% Pplus - (rho1 %*% Pplus)^2
    ## All we really need is V-tilde.  Chop off the last entry.
    V = V[-(bigmz+1)]

    ## There is no birth luck since everyone is born the same.
    ## Var_0 (\\Ex [R | x, z_0]) = 0.
      
    ## time to loop over age
    ## record luck and pluck terms as functions of age and x

    ## Pa = P^a, PaC0 = P^a %*% c_0
    Pa = diag(1, bigmz)
    PaC0 = bigc0

    ## eq. 6
    stateTrajecLuckAgeCondXBZ = V %*% Pa
    ## Take expectation over breeding prob. cond. on trait and initial
    ## size (eq. 20)
    expBStateTrajecLuckAgeCondXBZ =
      pBCondXZ*stateTrajecLuckAgeCondXBZ
    ## Take expectation over initial size (in theory, conditional on
    ## trait, but here traits don't affect initial size)
    expZBStateTrajecLuckAgeCondXBZ[1,ii] = expBStateTrajecLuckAgeCondXBZ %*%
      bigc0

    ## eq. 7
    fecLuckAgeCondXBZ = bigsigbsq %*% Pa
    ## Take expectation over breeding prob. cond. on trait and initial
    ## size (eq. 20)
    expBFecLuckAgeCondXBZ = pBCondXZ*fecLuckAgeCondXBZ
    ## Take expectation over initial size (in theory, conditional on
    ## trait, but here traits don't affect initial size)
    expZBFecLuckAgeCondXBZ[1,ii] = expBFecLuckAgeCondXBZ %*% bigc0

    eT = matrix(1,1,bigmz);
    ## eq. 26
    Ax = eT %*% F %*% N;
    bT = eT %*% (Fbar %*% Nbar - F %*% N) %*% P;
    rhoZeroACondXBZ = eT%*%(F + Fbar%*%Nbar%*%P);
    ## Take expectation over initial state conditional on breeding
    ## (and trait, in theory)
    expZRhoZeroACondXBZ[1,ii] = rhoZeroACondXBZ %*% bigc0CondBreed

    ## a is actually age + 1, since the first year of life is age 0
    for (a in 2:maxAge) {
      ## update P^a and P^a %*% c_0
      Pa = P %*% Pa
      PaC0 = Pa %*% bigc0

      ## eq. 6
      stateTrajecLuckAgeCondXBZ = V %*% Pa
      ## Take expectation over breeding prob. cond. on trait and initial
      ## size (eq. 20)
      expBStateTrajecLuckAgeCondXBZ =
        pBCondXZ*stateTrajecLuckAgeCondXBZ
      ## Take expectation over initial size (in theory, conditional on
      ## trait, but here traits don't affect initial size)
      expZBStateTrajecLuckAgeCondXBZ[a,ii] =
        expBStateTrajecLuckAgeCondXBZ %*% bigc0

      ## eq. 7
      fecLuckAgeCondXBZ = bigsigbsq %*% Pa
      ## Take expectation over breeding prob. cond. on trait and initial
      ## size (eq. 20)
      expBFecLuckAgeCondXBZ = pBCondXZ*fecLuckAgeCondXBZ
      ## Take expectation over initial size (in theory, conditional on
      ## trait, but here traits don't affect initial size)
      expZBFecLuckAgeCondXBZ[a,ii] = expBFecLuckAgeCondXBZ %*% bigc0

      ## eq. 26
      rhoZeroACondXBZ = Ax + bT %*% Pa;
      ## Take expectation over initial state conditional on breeding
      ## (and trait, in theory)
      expZRhoZeroACondXBZ[a,ii] = rhoZeroACondXBZ %*% bigc0CondBreed

    }  ## end loop over age

  }  ## end loop over trait values

  traitAve = function (x) {
    return (sum(starDist*x))
    ## return (sum(probXCondB*x))
  }

  traitVarCondBreed = function (x) {
    return (sum(pXCondB*x^2) - sum(pXCondB*x)^2)
  }

  ## Average over x: Ex_{x,M} (Var(R | x, M))  (eq. 21)
  stateTrajecLuckAgeWithinBreeding = expXZBStateTrajecLuckAgeCondXBZ =
    apply (expZBStateTrajecLuckAgeCondXBZ, 1, traitAve)
  fecLuckAgeWithinBreeding = expXZBFecLuckCondAgeXBZ =
    apply (expZBFecLuckAgeCondXBZ, 1, traitAve)

  ## Take variance over x conditional on breeding: Var_{x|M=1}  (eq. 28)
  varXExpZRhoZeroA = apply (expZRhoZeroACondXBZ, 1, traitVarCondBreed)
  ## Average over M  (eq. 25)
  expBVarXExpZRhoZeroA = pB * varXExpZRhoZeroA
  ## Find marginal value of each age (eq. 16)
  pluckAgeWithinBreeding = diff (c(0, expBVarXExpZRhoZeroA))


  out = list (stateTrajecLuckAgeWithinBreeding=stateTrajecLuckAgeWithinBreeding,
              fecLuckAgeWithinBreeding=fecLuckAgeWithinBreeding,
              pluckAgeWithinBreeding=pluckAgeWithinBreeding)

  return (out)
}

calcLuckPluck = function (Plist, F, b, sigbsq) {
  expZRhoZeroACondXZ = matrix(0,maxAge,numTraitVals)
  stateTrajecLuckAge = fecLuckAge = rep (0, maxAge)
  expZStateTrajecLuckAgeCondXZ = expZFecLuckAgeCondXZ =
    matrix (0, maxAge, numTraitVals)
  
  expRCondXZ = varRCondXZ = rep(NA, mz)
  expZVarRCondXZ = expZExpRCondXZ = numeric(numTraitVals)

  ## conditional on breeding.  
  Pbar = matrix (0, mz, mz)
  for (ii in 1:numTraitVals) {
    P = Plist[[ii]]
    Pbar = Pbar + starDist[ii] * P
  }
  ## For sneetches, F does not depend on x, so Fbar = F.
  Fbar = F  
  Nbar = solve (diag(mz) - Pbar)

  for (ii in 1:numTraitVals) {
    ## Transition matrix for this trait value
    P = Plist[[ii]]
    ## Fundamental matrix (N-tilde in the ms)
    N = solve (diag(mz) - P)

    ## Now that we have P, compute Pplus, defined in "Background and
    ## Assumptions." 
    Pplus = cbind (P, rep(0, mz))
    deathProb = 1 - apply (Pplus, 2, sum)
    Pplus = rbind (Pplus, deathProb)
    
    ## expected lifetime reproduction conditional on state (rho1)
    expRCondXZ = apply (F %*% N, 2, sum)
    ## rho1 includes RO for the state of being dead
    rho1 = c(expRCondXZ, 0)

    ## V = variance of rho1 (eq. 5)
    V = rho1^2 %*% Pplus - (rho1 %*% Pplus)^2
    ## All we really need is V-tilde.  Chop off the last entry.
    V = V[-(mz+1)]

    ## There is no birth luck since everyone is born the same.
    ## Var_0 (\\Ex [R | x, z_0]) = 0.
      
    ## time to loop over age
    ## record luck and pluck terms as functions of age and x

    ## Pa = P^a, PaC0 = P^a %*% c_0
    Pa = diag(1, mz)
    PaC0 = c0

    ## eq. 6
    stateTrajecLuckAgeCondXZ = V %*% Pa
    ## Take expectation over initial size (in theory, conditional on
    ## trait, but here traits don't affect initial size)
    expZStateTrajecLuckAgeCondXZ[1,ii] = stateTrajecLuckAgeCondXZ %*%
    c0

    ## eq. 7
    fecLuckAgeCondXZ = sigbsq %*% Pa
    ## Take expectation over initial size (in theory, conditional on
    ## trait, but here traits don't affect initial size)
    expZFecLuckAgeCondXZ[1,ii] = fecLuckAgeCondXZ %*% c0

    eT = matrix(1,1,mz);
    ## eq. 19
    Ax = eT %*% F %*% N;
    bT = eT %*% (Fbar %*% Nbar - F %*% N) %*% P;
    rhoZeroACondXZ = eT%*%(F + Fbar%*%Nbar%*%P);
    ## Take expectation over initial state (in theory conditional on
    ## trait) 
    expZRhoZeroACondXZ[1,ii] = rhoZeroACondXZ %*% c0


    ## a is actually age + 1, since the first year of life is age 0
    for (a in 2:maxAge) {
      ## update P^a and P^a %*% c_0
      Pa = P %*% Pa
      PaC0 = Pa %*% c0

      ## eq. 6
      stateTrajecLuckAgeCondXZ = V %*% Pa
      ## Take expectation over breeding prob. cond. on trait and initial
      ## size (eq. 20)
      expZStateTrajecLuckAgeCondXZ[a,ii] =
        stateTrajecLuckAgeCondXZ %*% c0

      ## eq. 7
      fecLuckAgeCondXZ = sigbsq %*% Pa
      ## Take expectation over breeding prob. cond. on trait and initial
      ## size (eq. 20)
      expZFecLuckAgeCondXZ[a,ii] = fecLuckAgeCondXZ %*% c0

      ## eq. 19
      rhoZeroACondXZ = Ax + bT %*% Pa;
      ## Take expectation over breeding prob. cond. on trait and initial
      ## size (eq. 20)
      expZRhoZeroACondXZ[a,ii] = rhoZeroACondXZ %*% c0

    }  ## end loop over age

  }  ## end loop over trait values

  traitAve = function (x) {
    return (sum(starDist*x))
  }
  
  traitVar = function (x) {
    return (sum(starDist*x^2) - sum(starDist*x)^2)
  }
  
  ## Average over x: Ex_x (Var(R | x))
  stateTrajecLuckAge = expXZStateTrajecLuckAgeCondXZ =
    apply (expZStateTrajecLuckAgeCondXZ, 1, traitAve)
  fecLuckAge = expXZFecLuckCondAgeXZ =
    apply (expZFecLuckAgeCondXZ, 1, traitAve)

  ## Take variance over x (eq. 16)
  varXExpZRhoZeroA = apply (expZRhoZeroACondXZ, 1, traitVar)
  ## Find marginal increase over age (eq. 17)
  pluckAge = diff (c(0, varXExpZRhoZeroA))

  out = list (stateTrajecLuckAge=stateTrajecLuckAge,
              fecLuckAge=fecLuckAge,
              pluckAge=pluckAge)

  return (out)
}

## Allocate space for results
pluckWithin = stateTrajecLuckWithin = fecLuckWithin = matrix(0, 3, maxAge)

####### Plain luck, no breeding partition ###############
Plist = list(Ug1, Ug2)
sigbsq = c(0, 0,
           pb[3]*(1 - pb[3])*theta[3]^2,
           pb[4]*(1 - pb[4])*theta[4]^2)
b = c(0, 0, pb[3]*theta[3], pb[4]*theta[4])
out = calcLuckPluck (Plist, F, b, sigbsq)
stateTrajecLuckAge = out$stateTrajecLuckAge
fecLuckAge = out$fecLuckAge
pluckAge = out$pluckAge

########################################################
## Breeder def. 1
#######################################################
mzS = 2
bigmz = mzS + mz

## Allocate memory for transition matrix list.
Plist = list(matrix (0, bigmz, bigmz), matrix (0, bigmz, bigmz))

## Allocate memory for probability of breeding
pBCondXZ = matrix (0, 2, bigmz)

## expanded state space init. distribution
bigc0 = c(1, 0, rep(0, mz))

## mean immediate RO
bigb = c(0, 0, 0, 0, pb[3]*theta[3], pb[4]*theta[4])

## variance in immediate RO
bigsigbsq = c(0, 0, 0, 0,
                pb[3]*(1 - pb[3])*theta[3]^2,
                pb[4]*(1 - pb[4])*theta[4]^2)
bigpb = c(0, 0, 0, 0, pb[3], pb[4])

out = makePCondBreedDef1(Ug1)
Plist[[1]] = out$PCondBreed
bigF = out$bigF
pBCondXZ[1,] = out$probEverBreed

out = makePCondBreedDef1 (Ug2)
Plist[[2]] = out$PCondBreed
pBCondXZ[2,] = out$probEverBreed

out = calcLuckPluckCondBreeding (bigmz, Plist, bigF, pBCondXZ, bigpb, bigb,
                            bigsigbsq, bigc0)
stateTrajecLuckWithin[1,] = out$stateTrajecLuckAgeWithinBreeding
fecLuckWithin[1,] = out$fecLuckAgeWithinBreeding
pluckWithin[1,] = out$pluckAgeWithinBreeding

## Sanity check ######
if (sum(apply(Plist[[1]], 2, sum) > 1.00000001))
  cat ("One or more columns of P1 exceed 1!\\n")
if (sum(apply(Plist[[2]], 2, sum) > 1.00000001))
  cat ("One or more columns of P2 exceed 1!\\n")

########################################################
## Breeder def. 2
#######################################################

bigmz = 3*mz
## Allocate memory for transition matrix list.
Plist = list(matrix (0, bigmz, bigmz), matrix (0, bigmz, bigmz))
## Allocate memory for probability of breeding
pBCondXZ = matrix (0, 2, bigmz)

## expanded state space init. distribution
bigc0 = c(1, 0, 0, 0, rep(0, mz), rep(0, mz))
## mean immediate RO
bigb = c(0, 0, 0, 0, 0, 0, theta[3], theta[4], 0, 0, pb[3]*theta[3],
         pb[4]*theta[4])
## variance in immediate RO
bigsigbsq = c(0, 0, 0, 0, 0, 0, theta[3], theta[4], 0, 0,
                pb[3]*(1 - pb[3])*theta[3]^2,
                pb[4]*(1 - pb[4])*theta[4]^2)
bigpb = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, pb[3], pb[4])

out = makePCondBreedDef23(Ug1,2)
Plist[[1]] = out$PCondBreed
bigF = out$bigF
pBCondXZ[1,] = out$probEverBreed

out = makePCondBreedDef23(Ug2,2)
Plist[[2]] = out$PCondBreed
pBCondXZ[2,] = out$probEverBreed

out = calcLuckPluckCondBreeding (bigmz, Plist, bigF, pBCondXZ, bigpb, bigb,
                            bigsigbsq, bigc0)
stateTrajecLuckWithin[2,] = out$stateTrajecLuckAgeWithinBreeding
fecLuckWithin[2,] = out$fecLuckAgeWithinBreeding
pluckWithin[2,] = out$pluckAgeWithinBreeding

## Sanity check ######
if (sum(apply(Plist[[1]], 2, sum) > 1.00000001))
  cat ("One or more columns of P1 exceed 1!\\n")
if (sum(apply(Plist[[2]], 2, sum) > 1.00000001))
  cat ("One or more columns of P2 exceed 1!\\n")

########################################################
## Breeder def. 3
#######################################################

bigmz = 3*mz
## Allocate memory for transition matrix list.
Plist = list(matrix (0, bigmz, bigmz), matrix (0, bigmz, bigmz))
## Allocate memory for probability of breeding
pBCondXZ = matrix (0, 2, bigmz)

## expanded state space init. distribution
bigc0 = c(1, 0, 0, 0, rep(0, mz), rep(0, mz))
## mean immediate RO
truncb3 = theta[3]/(1 - exp(-theta[3]))
truncb4 = theta[4]/(1 - exp(-theta[4]))
bigb = c(0, 0, 0, 0, 0, 0, truncb3, truncb4, 0, 0, theta[3], theta[4])
## variance in immediate RO
bigsigbsq = c(0, 0, 0, 0, 0, 0,
              truncb3*(1 + theta[3] - truncb3),
              truncb4*(1 + theta[4] - truncb4),
              0, 0, theta[3], theta[4])
bigpb = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, pb[3], pb[4])

out = makePCondBreedDef23(Ug1,3)
Plist[[1]] = out$PCondBreed
bigF = out$bigF
pBCondXZ[1,] = out$probEverBreed

out = makePCondBreedDef23(Ug2,3)
Plist[[2]] = out$PCondBreed
pBCondXZ[2,] = out$probEverBreed

out = calcLuckPluckCondBreeding (bigmz, Plist, bigF, pBCondXZ, bigpb, bigb,
                            bigsigbsq, bigc0)
stateTrajecLuckWithin[3,] = out$stateTrajecLuckAgeWithinBreeding
fecLuckWithin[3,] = out$fecLuckAgeWithinBreeding
pluckWithin[3,] = out$pluckAgeWithinBreeding

## Sanity check ######
if (sum(apply(Plist[[1]], 2, sum) > 1.00000001))
  cat ("One or more columns of P1 exceed 1!\\n")
if (sum(apply(Plist[[2]], 2, sum) > 1.00000001))
  cat ("One or more columns of P2 exceed 1!\\n")

#####################################################
## make some plots
##################################################

dev.new (width=10, height=10)
par(fig=c(0,1,0,1),yaxs="i",bty="l",xaxs="i",mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4,
    cex=1.4, mar=c(4,4,2,2), mfrow=c(2,2))

## Plot A
matplot (0:(maxAge-1), t(stateTrajecLuckWithin), type="l", lwd=2,
         xlab="Age", ylab="'Within' state trajec. luck",
         col=c("black", "red", "blue"))
legend (x="topright", col=c("black", "red", "blue"), lty=1:3, lwd=2,
        legend=c("Breeder def. 1", "Breeder def. 2",
                 "Breeder def. 3"), bty="n", cex=1.4)
add_panel_label ("a")

## Plot B
matplot (0:(maxAge-1), t(fecLuckWithin), type="l", lwd=2,
         col=c("black", "red", "blue"),
         xlab="Age", ylab="'Within' fecundity luck")
add_panel_label ("b")

## Plot C
totLuck = stateTrajecLuckAge + fecLuckAge
totLuckWithin = stateTrajecLuckWithin[3,] + fecLuckWithin[3,]
ymax = max(totLuck)
plot (0:(maxAge-1), totLuck, type="l", lwd=2,
      xlab="Age", ylab="Total luck", ylim=c(0, ymax))
lines (0:(maxAge-1), totLuck - totLuckWithin,
       lty=2, lwd=2)
legend (x="topright", legend=c("Within breeder, def. 3",
                               "Breed or not, def. 3"),
        col="black", lty=1:2, lwd=c(2,2), bty="n", cex=1.4)
add_panel_label ("c")

## Plot D
plot (0:(maxAge-1), pluckAge, col="black", type="l",
      lwd=2, xlab="Age", ylab="Pluck", ylim=c(0, max(pluckAge)))
lines (0:(maxAge-1), pluckAge - pluckWithin[3,], col="black",
       lty=2, lwd=2)
legend (x="topright", legend=c("Within breeder, def. 3",
                               "Breed or not, def. 3"),
        col="black", lty=1:2, lwd=c(2,2), bty="n", cex=1.4)
add_panel_label ("d")

dev.copy2pdf (file="SneetchCalculations3.pdf")

## Print some stuff to help us think about the effects of life history.
## What is the mean lifespan of each type?
N1 = solve (diag(4) - Ug1); N2 = solve (diag(4) - Ug2)
meanLifespan = rep(0, 2)
meanLifespan[1] = apply (N1, 2, sum)[1]
meanLifespan[2] = apply (N2, 2, sum)[1]
## What is the population mean lifespan, given trait distribution?
popMeanLifespan = sum(meanLifespan * starDist)

foo = which.max(stateTrajecLuckAge) 
cat ("State trajectory luck peaks at age",
     foo, "=", foo/popMeanLifespan, "of mean lifespan.\\n")
foo = which.max(stateTrajecLuckAge - stateTrajecLuckWithin[3,]) 
cat ("Breed-or-not state trajectory luck peaks at age",
     foo, "=", foo/popMeanLifespan, "of mean lifespan.\\n")
foo = which.max(stateTrajecLuckWithin[3,]) 
cat ("Within breeder state trajectory luck peaks at age",
     foo, "=", foo/popMeanLifespan, "of mean lifespan.\\n")
cat ("Max.(state trajectory luck) =", max(stateTrajecLuckAge), "\\n")

