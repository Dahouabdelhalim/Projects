#####################################################################
# 
# Sneetches are classified as Small, Medium, Large, and Xtra-large. 
# Small and Medium sneeches do not breed (p_b =0). 
# Large and Xtra-large sneetches have p_b>0, and Poisson offspring number
# conditional on attempting to breed. 
#
# To be somewhat parameter-sparse, we give each size a probability of survival 
# and probability of stasis (no change in size). If changing size, S grows to M,
# M and L shrink or grow with equal odds, and X shrinks to L. 
#
# Starred sneetches (group 2) have lower mortality than starless in all stages,
# with the same fractional reduction in all stages. Otherwise they are identical
# to starless. 
#
# We use the van Daalen & Caswell notation in which A = U + F
# is the projection matrix on living states, while P (which we call
# P+ in the paper) is the projection matrix including \\omega = dead.  
# Other notation follows the Fulmar scripts. 
#  
#####################################################################

maxAge = 100
numTraitVals = 2

## Which classes are S and which are M?
indicatorS = c(1, 1, 0, 0)
indicatorM = 1-indicatorS

## Make up some parameter values for right now.

surv = c(0.5, 0.6, 0.7, 0.8)
starFac = 0.8
stay = c(0.3, 0.3, 0.6, 0.6)
theta = c(0, 0, 1, 2)
pb = c(0, 0, 0.3, 0.5)
c0 = c(1, 0, 0, 0)

stay[3:4] = rep(0.9, 2)
if (slowDevel) {
  stay[1:2] = rep(0.9, 2)
} else if (fastDevel) {
  stay[1:2] = rep(0.1, 2)
}

if (highJuvMort) {
  surv[1:2] = rep(0.1, 2)
} else if (lowJuvMort) {
  surv[1:2] = rep(0.9, 2)
}

if (shortRepro) {
  surv[3:4] = rep(0.1, 2)
} else if (longRepro) {
  surv[3:4] = rep(0.9, 2)
}


## The proportions of c(plain-bellied, star-bellied) sneetches. 
starDist = c(0.6, 0.4)

##################################################################### 
# Make the projection matrices U and F for each group. 
# We use van Daalen & Caswell notation for matrix sizes 
#####################################################################
mz = 4; s = mz+1;

## survival for star-belly sneetches
surv2 = 1 - starFac*(1 - surv)

Ug1 = Ug2 = matrix(0,mz,mz);
diag(Ug1) = surv*stay;  diag(Ug2) = surv2*stay
go = 1-stay; 
Ug1[2,1] = surv[1]*go[1]; 
Ug1[1,2] = Ug1[3,2] = surv[2]*go[2]/2;
Ug1[2,3] = Ug1[4,3] = surv[3]*go[3]/2;
Ug1[3,4]= surv[4]*go[4];  

Ug2[2,1] = surv2[1]*go[1]; 
Ug2[1,2] = Ug2[3,2] = surv2[2]*go[2]/2;
Ug2[2,3] = Ug2[4,3] = surv2[3]*go[3]/2;
Ug2[3,4]= surv2[4]*go[4]; 

cat(all.equal(apply(Ug1,2,sum),surv),"\\n"); # sanity check is passed
cat(all.equal(apply(Ug2,2,sum),surv2),"\\n"); # sanity check is passed 

F = matrix(0,mz,mz); 
F[1,] = theta; 

######################################################### 
# Function to make big P, Definition 1 
#########################################################
makePDef1 = function (U) {
  mzS = 2; bigmz = mzS + mz
  bigP = matrix (0, bigmz, bigmz)
  
  ## Make P for the extended state space Z1 \\cup Z2
  ## Z1 to Z1
  bigP[1:mzS, 1:mzS] = U[1:mzS, 1:mzS]
  ## Z1 to Z2
  for (j in 1:mz) bigP[mzS+j, 1:mzS] = U[j,1:mzS]*indicatorM[j]
  ## Z2 to Z2
  bigP[(mzS+1):bigmz, (mzS+1):bigmz] = U
  return(bigP) 
} 

Ug1; makePDef1(Ug1);  # check, everything is as it should be 


######################################################### 
# Function to make big P, Definition 2  
#########################################################
makePDef2 = function (U) {
  bigmz = 3*mz;
  bigP = matrix (0, bigmz, bigmz);
  
  ## Make P for the extended state space Z1 \\cup Z2 \\cup Z3
  # Z1 to Z1
  for (j in 1:mz) bigP[j, 1:mz] = U[j,]*(1 - pb[j])
  # Z1 to Z2 
  for (j in 1:mz) bigP[mz+j, 1:mz] = U[j,]*pb[j]
  # Z2 to Z3
  bigP[(2*mz+1):bigmz, (mz+1):(2*mz)] = U
  # Z3 to Z3
  bigP[(2*mz+1):bigmz, (2*mz+1):bigmz] = U
  return(bigP); 
}
Ug2; bigP=makePDef2(Ug2); bigP; apply(bigP,2,sum)-rep(surv2,3); #everything checks out 


######################################################### 
# Function to make big P, Definition 3  
#########################################################
makePDef3 = function (U) {
  bigmz = 3*mz
  bigF = bigP = P0 = matrix (0, bigmz, bigmz)

  ## calculate pd, the probability of producing at least one offspring
  ## in the current year.
  pd = rep(0, mz)
  pd[3] = 1 - exp(-theta[3])
  pd[4] = 1 - exp(-theta[4])
  
  ## Make P for the extended state space Z1 \\cup Z2 \\cup Z3
  # Z1 to Z1 
  for (j in 1:mz)  bigP[j, 1:mz] = U[j,]*(1 - pd[j])
  # Z1 to Z2
  for (j in 1:mz)  bigP[mz+j, 1:mz] = U[j,]*pd[j]
  # Z2 to Z3
  bigP[(2*mz+1):bigmz, (mz+1):(2*mz)] = U
  # Z3 to Z3
  bigP[(2*mz+1):bigmz, (2*mz+1):bigmz] = U
  
  return(bigP)
}  
Ug2; bigP=makePDef3(Ug2); bigP; apply(bigP,2,sum)-rep(surv2,3); #checks out 


########################################################################## 
# Function to make big P conditional on breeding, Definition 1 
##########################################################################  
makePCondBreedDef1 = function (U) {
  mzS = 2; bigmz = mzS + mz;
  
  bigP = makePDef1(U); 
  bigF = matrix (0, bigmz, bigmz)
  
  ## expanded state space init. distribution
  bigc0 = c(1, 0, rep(0, mz))

  ## make F for the extended state space Z1 \\cup Z2.  Note that the
  ## only non-zero entries are from Z2 to Z1.
  bigF[1, (mzS+1):bigmz] = F[1,]
  
  ## Calculate aM, the probability of going from Z1 to Z2 in one step
  aM = apply (U[(mzS+1):mz, 1:mzS], 2, sum)

  ## Q1 = U restricted to Z1, N1 = (I - Q1)^{-1}
  Q1 = U[1:mzS, 1:mzS]
  N1 = solve(diag(mzS) - Q1)

  ## Calculate probEverBreed (called B(z-hat) in the ms), the
  ## probability of reaching Z2 before death.
  probEverBreed = c(aM %*% N1, rep(1, mz))
  pM = sum(bigc0*probEverBreed);  
  
  PCondBreed = matrix (NA, bigmz, bigmz)
  for (j in 1:bigmz) {
    PCondBreed[,j] = bigP[,j]*probEverBreed/probEverBreed[j]
  }
  
  out = list (PCondBreed=PCondBreed, probEverBreed=probEverBreed,
  				pM=pM, bigF=bigF)
  return (out)
}

out=makePCondBreedDef1(Ug1); 
apply(out$PCondBreed,2,sum); # what it should be; if still in Z1 you don't die. 


########################################################################## 
# Function to make big P conditional on breeding, Definitions 2 and 3
##########################################################################  
makePCondBreedDef23 = function (U,def) {
  bigmz = 3*mz
  bigF = matrix (0, bigmz, bigmz)
  
  bigP=NULL; 
  if(def==2) bigP=makePDef2(U)
  if(def==3) bigP=makePDef3(U)

  ## expanded state space init. distribution
  bigc0 = c(1, 0, 0, 0, rep(0, 2*mz))

  ## make F for the extended state space Z1 \\cup Z2 \\cup Z3.  Note that the
  ## only non-zero entries are from Z2 and Z3 to Z1.
  bigF[1, (mz+1):bigmz] = F[1,]
  
  
  ## Calculate aM, the probability of going from Z1 to Z2 in one step
  aM = apply (bigP[(mz+1):(2*mz), 1:mz], 2, sum)
  
  ## Q2 = P restricted to Z1, N2 = (I - Q2)^{-1}
  Q2 = bigP[1:mz,1:mz];

  N2 = solve(diag(mz) - Q2)

  ## Calculate probEverBreed (called B(z-hat) in the ms), the
  ## probability of reaching Z2 before death.
  probEverBreed = c(aM %*% N2, rep(1, 2*mz))
  
  PCondBreed = matrix (NA, bigmz, bigmz)
  for (j in 1:bigmz) {
    PCondBreed[,j] = bigP[,j]*probEverBreed/probEverBreed[j]
  }

  out = list (PCondBreed=PCondBreed, probEverBreed=probEverBreed, 
  				pM=sum(bigc0*probEverBreed), bigF=bigF)
  return (out)
}

out=makePCondBreedDef23(Ug1,2); 
apply(out$PCondBreed,2,sum); # what it should be; if still in Z1 you don't die. 

