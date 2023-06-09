#################################################################
##
## roePartitioningNewLuckPluck4: Partitions luck and pluck and
## breeding status according to age for the roe deer model.
##
#################################################################

rm(list=ls(all=TRUE)); graphics.off(); 
require(Matrix); require(statmod);

path = "/home/robin/ipm/drivers/Manuscripts/Luckpartition/code/"

if(.Platform$OS.type=="windows") {
  path = "c:/repos/drivers/Manuscripts/luckPartitioning"
}

library(spam)
library(boot)

##maxAge = 200  ## Good for getting sums to converge but don't plot that high
maxAge = 50  ## Realistic max. age for roe deer

## Read in params, make P and F.
cat ("Making P.\\n")
source ("roeMakeP.R")

## Number of cross-classified states = (# size classes)*(#
## birthdates)*(# age classes), but we never use this.
## totStates = nb*nd*n.age

## Number of states per birthdate = yearlings with the right bday and
## adults with the right bday
mz = 2*nb

## Get stable dist.
cat ("Finding eigenstuff.\\n")
out <- get.eigen.stuff(P+F)
evec1 <- as.vector(out[[2]]) ###stable distribution
cat ("Done finding eigenstuff.\\n")

yearlings = 1:(nb*nd)
adults = (nb*nd + 1):(nb*nd*n.age)
## Get stable birthdate dist.
bdDist = apply (matrix(evec1[yearlings], nrow=nb, ncol=nd), 2, sum)
bdDist = bdDist / sum(bdDist)

############ Normal luck/pluck calculation ##############################
calcLuckPluckStochToD = function () {

  expRCondBDZ = varRCondBDZ = matrix (0, mz, nd)
  stateTrajecLuckAgeBD = fecLuckAgeBD = expRCondBDFnAge = matrix (0, maxAge, nd)

  stateTrajecLuckAge = fecLuckAge = PluckAge = rep (0, maxAge)
  expRCondBD = rep(0, nd)
  expBDProductAge = rep(0, maxAge)
  survAgeBD = matrix (0, maxAge, nd)

  birthLuck = rep(0, nd)

  ## Find average kernels.
  Pbar = matrix(0, mz, mz)
  Fbar = matrix (0, nd*nb, mz)
  for (bd in 1:nd) {
    ## Yearlings with the right birthdate
    yearlingsBD = ((bd-1)*nb)+(1:nb)
    ## Adults with the right birthdate
    adultsBD = nb*nd + ((bd-1)*nb)+(1:nb)

    ## Subset P appropriately --- moms and kids with the right bday
    bdP = P[c(yearlingsBD, adultsBD), c(yearlingsBD, adultsBD)]
    ## Subset F appropriately --- moms with the right bday, all kids
    bdF = F[yearlings, c(yearlingsBD, adultsBD)]  

    Fbar = Fbar + bdF; Pbar = Pbar + bdP
  }

  Fbar = Fbar / nd; Pbar = Pbar / nd
  Nbar = solve(diag(mz) - Pbar)

## Loop over birthdates, get a P for each birth date. 

  rhoZeroA = matrix(0,maxAge,nd);
  expRCondBDStochToD = varRCondBDStochToD = matrix (0, maxAge, nd)
  expZVarRCondBDZ = varZExpRCondBDZ = rep (0, nd)

  for (bd in 1:nd) {
    varRCondBDZStochToD = matrix (0, maxAge, mz)
    
    cat ("bd = ", bd, "\\n")
    ## Yearlings with the right birthdate
    yearlingsBD = ((bd-1)*nb)+(1:nb)
    ## Adults with the right birthdate
    adultsBD = nb*nd + ((bd-1)*nb)+(1:nb)

    ## Subset P appropriately --- moms and kids with the right bday
    bdP = P[c(yearlingsBD, adultsBD), c(yearlingsBD, adultsBD)]
    ## Subset F appropriately --- moms with the right bday, all kids
    bdF = F[yearlings, c(yearlingsBD, adultsBD)]

    ## initial condition is kids with the right bday, following the
    ## stationary birth mass distribution
    c0 = evec1[yearlingsBD]
    c0 = c0 / sum(c0)
    ## Add zeros for adult states
    c0 = c(c0, rep(0, nb))

    Pplus = cbind (bdP, rep(0, mz))
    deathProb = 1 - apply (Pplus, 2, sum)
    Pplus = rbind (Pplus, deathProb)
    ## Fundamental matrix (N-tilde in the ms)
    N = solve (diag(mz) - bdP)

    ## Plard et al. don't describe breeding prob.  All is subsumed into
    ## F, which they call R.  Let pb be 1 for all adults.
    pb = c(rep(0, nb), rep(1, nb))

    ## In this model, everyone is assumed to produce twins.  Remember to
    ## divide by 2 since we're tracking females only.
    b = apply (bdF, 2, sum)
    ## We assume each mom produces 1 daughter, so the column sum of F =
    ## the survival probability.
    survProb = b  
    ## Variance in number of surviving daughters is np(1 - p), where n =
    ## litter size = 2, p = (sex ratio)*(survival prob.)
    pp = 0.5*survProb
    sigbsq = 2*pp*(1 - pp)

    ## expected lifetime reproduction conditional on state (rho1)
    expRCondBDZ[,bd] = apply (bdF %*% N, 2, sum)
    ## rho1 includes RO for the state of being dead
    rho1 = c(expRCondBDZ[,bd], 0)

    ## V (eq. 7)
    V = rho1^2 %*% Pplus - (rho1 %*% Pplus)^2
    ## All we really need is V-tilde.  Chop off the last entry.
    V = V[-(mz+1)]

    ## Birth luck = Var_0 (\\Ex [R | x, z_0]) = 0.
    birthLuck[bd] = sum(expRCondBDZ[,bd]^2 * c0) -
      sum(expRCondBDZ[,bd] * c0)^2

    ## time to loop over age
    ## record luck and pluck terms as functions of age and x
    Pa = diag(1, mz)
    PaC0 = c0
    survAgeBD[1,bd] = sum (PaC0)
    stateTrajecLuckAgeBD[1,bd] = V %*% PaC0
    fecLuckAgeBD[1,bd] = sigbsq %*% PaC0
    expRCondBDFnAge[1,bd] = (pb*b) %*% PaC0
    eT = matrix(1,1,nd*nb); 
    Ax = eT %*% bdF %*% N %*% c0;  
    bT = eT %*% (Fbar %*% Nbar - bdF %*% N) %*% bdP; 
    rhoZeroA[1,bd] =  eT%*%(bdF + Fbar%*%Nbar%*%bdP)%*%c0;
    PaSeries = Pa = NM = diag(1, mz)
    FD = outer (c0, expRCondBDZ[,bd])
    bigrbarPib = c(apply (FD%*%bdP, 2, sum), rep(0, mz))

    ## rho1Mmz the first mz entries of rho1M
    rho1Mmz = apply(bdF, 2, sum) + apply (FD%*% bdP, 2, sum)
    expRCondBDZStochToD = rho1Mmz
    expRCondBDStochToD[1,bd] = expRCondBDZStochToD %*% c0
    rho1 = apply (bdF%*%N, 2, sum)

    ## NMmz = the first mz-wide column of NM
    NMmz = rbind (diag(1, mz), bdP)
    sigbsqM = c(sigbsq, rep(0, mz))
    pbM = c(pb, rep(1, mz))
    bM = c(b, rho1)
    ## r2Mmz = the first mz entries of r2M
    r2Mmz = (sigbsqM + (pbM*bM)^2 + 2*pbM*bM*bigrbarPib) %*% NMmz
    ##  varRCondBDStochToD[1,bd] = (r2Mmz - rho1Mmz^2) %*% c0
    varRCondBDZStochToD[1,] = r2Mmz - rho1Mmz^2
    ## Var(R | x) = Ex_z(Var(R | x, z)) + Var_z(Ex(R | x, z))
    varRCondBDStochToD[1,bd] = varRCondBDZStochToD[1,]%*%c0 +
      c0%*%(expRCondBDZStochToD^2) - (c0 %*% expRCondBDZStochToD)^2
    
    ## a is actually age + 1, since the first year of life is age 0
    for (a in 2:maxAge) {
      Pa = bdP %*% Pa
      PaC0 = Pa %*% c0
      survAgeBD[a,bd] = sum (PaC0)
      stateTrajecLuckAgeBD[a,bd] = V %*% PaC0
      fecLuckAgeBD[a,bd] = sigbsq %*% PaC0
      ## apply (expRCondBDFnAge, 2, sum) should = expRCondBD and it does.
      expRCondBDFnAge[a,bd] = (pb*b) %*% PaC0
      rhoZeroA[a,bd] = Ax + bT %*% PaC0;
      PaSeries = PaSeries + Pa  ## I + P + P^2 + ... + P^(a-1)
      NMmz = rbind (NMmz, Pa%*%bdP)
      ## The first mz entries of rho1
      rho1Mmz = apply (bdF%*%PaSeries, 2, sum) + apply (FD%*%Pa%*%bdP, 2, sum)
      expRCondBDZStochToD = rho1Mmz
      expRCondBDStochToD[a,bd] = expRCondBDZStochToD %*% c0

      pbM = c(rep(pb, a), rep(1, mz))
      sigbsqM = c(rep(sigbsq, a), rep(0, mz))
      bM = c(rep(b, a), rho1)
      newbit = apply (bdF%*%(PaSeries - diag(1, mz)), 2, sum) +
        apply (FD%*%Pa%*%bdP, 2, sum)
      bigrbarPib = c(newbit, bigrbarPib)
      ## The first mz entries of r2M
      r2Mmz = (sigbsqM +
               (pbM*bM)^2 + 2*pbM*bM*bigrbarPib) %*% NMmz
      varRCondBDZStochToD[a,] = r2Mmz - rho1Mmz^2
      ## Var(R | x) = Ex_z(Var(R | x, z)) + Var_z(Ex(R | x, z))
      varRCondBDStochToD[a,bd] = varRCondBDZStochToD[a,]%*%c0 +
        c0%*%(expRCondBDZStochToD^2) - (c0 %*% expRCondBDZStochToD)^2
    }
  }  ## end loop over birthdate

  ## Perform appropriate averages over birthdate.
  
  ## Age-dependent ######
  stateTrajecLuckAge = stateTrajecLuckAgeBD %*% bdDist
  fecLuckAge = fecLuckAgeBD %*% bdDist

  ## new pluck #############
  bdVar=function(x) {
    xbar=sum(bdDist*x); x2bar=sum(bdDist*x^2);
    return(x2bar-xbar^2)
  }    
  vaPluck = apply(rhoZeroA,1,bdVar); 
  NewPluckAge=diff(c(0,vaPluck));  

  ## new luck #############
  bdAve = function (x) {
    return (sum(bdDist*x))
  }
  exLuck = apply (varRCondBDStochToD, 1, bdAve)
  NewLuckAge = diff (c(0, exLuck))

  out = list (NewLuckAge=NewLuckAge,
              NewPluckAge=NewPluckAge,
              stateTrajecLuckAge=stateTrajecLuckAge,
              fecLuckAge=fecLuckAge)
  return (out)
}

calcLuckPluck = function () {

  expZRhoZeroACondXZ = matrix(0,maxAge,nd);
  expZStateTrajecLuckAgeCondXZ = expZFecLuckAgeCondXZ =
    matrix (0, maxAge, nd)
  birthLuck = numeric(nd)
  expRCondXZ = rep (NA, mz)

  ## We need Fbar and Pbar for pluck.  Average kernels over trait
  ## conditional on breeding.  
  Pbar = matrix (0, mz, mz)
  Fbar = matrix (0, nd*nb, mz)

  for (ii in 1:nd) {
    ## Yearlings with the right birthdate
    yearlingsBD = ((ii-1)*nb)+(1:nb)
    ## Adults with the right birthdate
    adultsBD = nb*nd + ((ii-1)*nb)+(1:nb)

    ## Subset P appropriately --- moms and kids with the right bday
    bdP = P[c(yearlingsBD, adultsBD), c(yearlingsBD, adultsBD)]
    ## Subset F appropriately --- moms with the right bday, all kids
    bdF = F[yearlings, c(yearlingsBD, adultsBD)]  

    Fbar = Fbar + bdF; Pbar = Pbar + bdP
  }

  Fbar = Fbar / nd; Pbar = Pbar / nd
  Nbar = solve (diag(mz) - Pbar)

  for (ii in 1:nd) {
    cat ("ii = ", ii, "\\n")
    ## Yearlings with the right birthdate
    yearlingsBD = ((ii-1)*nb)+(1:nb)
    ## Adults with the right birthdate
    adultsBD = nb*nd + ((ii-1)*nb)+(1:nb)

    ## Subset P appropriately --- moms and kids with the right bday
    bdP = P[c(yearlingsBD, adultsBD), c(yearlingsBD, adultsBD)]
    ## Subset F appropriately --- moms with the right bday, all kids
    bdF = F[yearlings, c(yearlingsBD, adultsBD)]

    ## Fundamental matrix (N-tilde in the ms)
    N = solve (diag(mz) - bdP)
    
    ## initial condition is kids with the right bday, following the
    ## stationary birth mass distribution
    c0 = evec1[yearlingsBD]
    c0 = c0 / sum(c0)
    ## Add zeros for adult states
    c0 = c(c0, rep(0, nb))

    Pplus = cbind (bdP, rep(0, mz))
    deathProb = 1 - apply (Pplus, 2, sum)
    Pplus = rbind (Pplus, deathProb)

    ## Plard et al. don't describe breeding prob.  All is subsumed into
    ## F, which they call R.  Let pb be 1 for all adults.
    pb = c(rep(0, nb), rep(1, nb))

    ## In this model, everyone is assumed to produce twins.  Remember to
    ## divide by 2 since we're tracking females only.
    b = apply (bdF, 2, sum)
    ## We assume each mom produces 1 daughter, so the column sum of F =
    ## the survival probability.
    survProb = b  
    ## Variance in number of surviving daughters is np(1 - p), where n =
    ## litter size = 2, p = (sex ratio)*(survival prob.)
    pp = 0.5*survProb
    sigbsq = 2*pp*(1 - pp)

    ## expected lifetime reproduction conditional on state (rho1)
    expRCondXZ = apply (bdF %*% N, 2, sum)
    ## rho1 includes RO for the state of being dead
    rho1 = c(expRCondXZ, 0)
    
    ## V (eq. 7)
    V = rho1^2 %*% Pplus - (rho1 %*% Pplus)^2
    ## All we really need is V-tilde.  Chop off the last entry.
    V = V[-(mz+1)]
      
    ## Birth luck = Var_0 (\\Ex [R | x, z_0])
    birthLuck[ii] = sum(expRCondXZ^2 * c0) -
      sum(expRCondXZ * c0)^2
      
    ## time to loop over age
    ## record luck and pluck terms as functions of age and x
    Pa = diag(1, mz)
    PaC0 = c0

    stateTrajecLuckAgeCondXZ = V %*% Pa
    expZStateTrajecLuckAgeCondXZ[1,ii] = stateTrajecLuckAgeCondXZ %*%
      c0

    fecLuckAgeCondXZ = sigbsq %*% Pa
    expZFecLuckAgeCondXZ[1,ii] = fecLuckAgeCondXZ %*% c0

    eT = matrix(1,1,nd*nb); 
    Ax = eT %*% bdF %*% N;
    bT = eT %*% (Fbar %*% Nbar - bdF %*% N) %*% bdP;
    rhoZeroACondXZ = eT%*%(bdF + Fbar%*%Nbar%*%bdP);
    expZRhoZeroACondXZ[1,ii] = rhoZeroACondXZ %*% c0
    
    ## a is actually age + 1, since the first year of life is age 0
    for (a in 2:maxAge) {
      Pa = bdP %*% Pa
      PaC0 = Pa %*% c0
      
      stateTrajecLuckAgeCondXZ = V %*% Pa
      expZStateTrajecLuckAgeCondXZ[a,ii] =
        stateTrajecLuckAgeCondXZ %*% c0
      
      fecLuckAgeCondXZ = sigbsq %*% Pa
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
    return (sum(bdDist*x))
  }

  traitVar = function (x) {
    return (sum(bdDist*x^2) - sum(bdDist*x)^2)
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
      birthLuck=birthLuck,
      pluckAge=pluckAge)

  return (out)
}

out = calcLuckPluck ()
fecLuckAge = out$fecLuckAge
stateTrajecLuckAge = out$stateTrajecLuckAge
birthLuck = out$birthLuck
NewLuckAge = fecLuckAge + stateTrajecLuckAge
NewPluckAge = out$pluckAge

############################################################
## Make plots.
############################################################

require(viridis); cols=plasma(100);

add_panel_label <- function(ltype="a") {
    text <- paste(LETTERS[letters==ltype], ")", sep="")
    mtext(text=text, side=3, adj=0)
}

graphics.off(); dev.new(width=10,height=5);
##graphics.off(); dev.new(width=10,height=10);
par(fig=c(0,1,0,1),yaxs="i",bty="l",xaxs="i",mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4,mar=c(4,4,2,2), mfrow=c(1,2))


## Plot a
plot (0:50, NewLuckAge[1:51], type="l", xlab="Age", ylab="Contribution to Var(LRO)",lwd=2)
lines (0:50, NewPluckAge[1:51], col="red",lwd=2)
lines (0:50, stateTrajecLuckAge[1:51], col="black",lwd=2, lty=2)
lines (0:50, fecLuckAge[1:51], col="black",lwd=2, lty=3)
legend ("topright",inset=0,col=c("black", "red", "black", "black"),
        lty=c(1,1,2,3), bty="n", lwd=2,cex=1.5,
        legend=c("Luck", "Birthdate", "State trajec. luck", "Fec. luck"))
add_panel_label ("a")

## Plot b
plot (0:50, (NewPluckAge/(NewPluckAge + NewLuckAge))[1:51], type="l",
      xlab="Age", ylab="Fraction contribution from birthdate", lwd=2)
add_panel_label ("b")

##dev.copy2pdf(file="roePartitioningNewLuckPluck4.pdf");
