#########################################################################
##
## DacrydiumPartitioning: Partitions Var(R) by age and size for
## Dacrydium elatum. Model from Zuidema et al. 2010, J Ecology 98: 345–355
## and previously used by Snyder and Ellner (2016). 
##
#########################################################################
rm(list=ls(all=TRUE)); 
require(Matrix); 
require(statmod); 

############ Edit as needed to set working directory 
path ="/home/robin/ipm/drivers/Manuscripts/luckPartition/"
out=try(setwd(path),silent=TRUE);
if(class(out)=="try-error") { 
    path=ifelse(.Platform$OS.type=="windows","c:/repos/drivers/Manuscripts/luckPartition","~/repos/drivers/Manuscripts/LuckPartition"); 
    setwd(path); 
  }

### These are copied from Snyder & Ellner (Am Nat, 2016). 
### For explanation of how they work, look here; for here you just 
### need to know that they build the iteration matrices for and IPM with
### some discrete "sapling" classes and trees classified by size
source("ZuidemaParams.R"); source ("ZuidemaTreeFunctions.R")

source("Standard Graphical Pars.R"); source("domEig.R")
source("Utilities.R")

###############################################################################
# Construct IPM matrices
###############################################################################

maxKids = 120   ## Set max. number of kids
maxAge = 600    ## Set max. age
species = 3     ## Choose which species we're focusing on

# IPM parameters
zmin=1; zmax=125; m=700

# Construct growth/survival matrix P for classes defined only by size
out <- mk_full_P_noloop(m=m,L=zmin,U=zmax,order=9,pars=param[[species]],A1=A1s[[species]]);
P <- out$P; meshpts <- out$meshpts; 
allMesh <- c(-3,-2,-1, 0, meshpts) # artificial sizes for saplings -- just used in plotting 
mz=ncol(P);
Pplus = cbind (P, rep(0, mz))
deathProb = 1 - apply (Pplus, 2, sum)
Pplus = rbind (Pplus, deathProb)
## Fundamental matrix (N-tilde in the ms)
N = solve (diag(mz) - P)

## Contruct F; 
F=matrix(0,mz,mz); pars=param[[species]]; 
F[1,] <- p.r(allMesh,pars)*pars[11]; 

## Find variance in number of kids.
## Distribution is Poisson(lambda) conditional on flowering, else 0 
pars = param[[species]]; lambda = pars[11]; 
sigbsq = pr*lambda + pr*(1-pr)*lambda^2; 

## Set up birth size distribution 
c0 = c(1, rep(0, mz-1));

stateTrajecLuckAge = fecLuckAge = PluckAge = Survival = rep(0, maxAge)
StageDist = matrix(0,mz,maxAge); 

###########################################################
# Compute mean and variance given state z 
###########################################################
expRCondZ = varRCondZ = rep (0, mz)

## expected lifetime reproduction conditional on state (rho1)
expRCondZ = apply (F %*% N, 2, sum)
## rho1 includes RO for the state of being dead
rho1 = c(expRCondZ, 0)

## V, eqn. {eqn:varbarr}
V = rho1^2 %*% Pplus - (rho1 %*% Pplus)^2
## All we really need is V-tilde.  Chop off the last entry.
V = V[-(mz+1)]

## There is no birth luck since everyone is born the same size.
## Var_0 (\\Ex [R | x, z_0]) = 0.

####################################################################
## Elasticity analysis and "manage for luck" 
####################################################################
A = P+F; ev=eigen(A); 
lambda=Re(ev$values[1]); w = ev$vectors[,1]; w=Re(w); w=w/sum(w); 
v=eigen(t(A))$vectors[,1];v=Re(v)/Re(v[1]); 
sens=outer(v,w)/sum(v*w);
elas=A*sens/lambda;
matrix.image(elas); 
which(elas==max(elas),arr.ind=TRUE); 
# [1,]   4   4  <--- survival of largest sapling stage 

which(V==max(V)); 

#############################################################################
# Loop over age to compute the age-specific terms in the 
# age decomposition, eqn. 7 in the paper. 
#############################################################################
Pa = diag(1, mz); PaC0 = c0
stateTrajecLuckAge[1] = V %*% PaC0
fecLuckAge[1] = sigbsq %*% PaC0
Survival[1]=1; 
StageDist[,1]=c0; 

## a is actually age + 1, since the first year of life is age 0
for (a in 2:maxAge) {
  Pa = P %*% Pa
  PaC0 = Pa %*% c0
  stateTrajecLuckAge[a] = V %*% PaC0
  fecLuckAge[a] = sigbsq %*% PaC0
  Survival[a]=sum(PaC0); 
  StageDist[,a]=PaC0/sum(PaC0); 
  if(a%%10==0) cat(a,"\\n"); 
}
h=(zmax-zmin)/m; 
StageDist[5:mz,] = StageDist[5:mz,]/h; 


## record luck and pluck terms as functions of z 
Nc0 = N %*% c0; expR = aveKids %*% Nc0 ## This is sum(F %*% N %*% c0)
stateTrajecLuckZ = V * Nc0
fecLuckZ = sigbsq * Nc0

## For sanity check: variance of total LRO from now on, conditional on state
rbarPib = expRCondZ %*% P     # \\bar{r} \\pi_b
r2 = (sigbsq + (aveKids)^2 + 2*aveKids*rbarPib) %*% N
varRCondZ = r2 - expRCondZ^2

## Sanity check 
cat(sum(stateTrajecLuckAge + fecLuckAge),sum(stateTrajecLuckZ + fecLuckZ), varRCondZ[1],"\\n") 

Zterms=stateTrajecLuckAge;
Fterms=fecLuckAge; 

##################################################################################
# Decompose again, this time for breeders only 
# We use definition 2, based on the published model's p.r function
# Notation follows the Sneetch example, see SneetchMatrices.R and SneetchCalculations.R
###################################################################################

# Make transition matrix again, to make sure we have the right starting point
pars = param[[species]]; 
out <- mk_full_P_noloop(m=m,L=zmin,U=zmax,order=9,pars=param[[species]],A1=A1s[[species]]);
U <- out$P;  meshpts <- out$meshpts; allMesh = c(-3,-2,-1, 0,meshpts) 
pb = p.r(allMesh,pars); mz = length(allMesh); 

# Make the unconditional \\mathbf{P} for extended state space Z1 \\cup Z2 \\cup Z3
# Follows the Sneetch code very closely! 
bigmz = 3*mz; bigP = matrix (0, bigmz, bigmz);
for (j in 1:mz) bigP[j, 1:mz] = U[j,]*(1 - pb[j]);   # Z1 to Z1
for (j in 1:mz) bigP[mz+j, 1:mz] = U[j,]*pb[j];      # Z1 to Z2 
bigP[(2*mz+1):bigmz, (mz+1):(2*mz)] = U;             # Z2 to Z3
bigP[(2*mz+1):bigmz, (2*mz+1):bigmz] = U;            # Z3 to Z3

## make F for the extended state space Z1 \\cup Z2 \\cup Z3.  Note that the
## only non-zero entries are from Z2 and Z3 to Z1.
bigaveKids = rep(NA,bigmz); 
bigaveKids[1:mz] = 0;                 # on Z1 there is no flowering 
bigaveKids[(mz+1):(2*mz)] = lambda;    # on Z2, flowering is certain 
bigaveKids[(2*mz+1):bigmz] = aveKids;  # on Z3, the standard rules apply 
bigF = matrix (0, bigmz, bigmz)
bigF[1,] = bigaveKids; 

## Also find variance in immediate number of kids on expanded state space 
bigsigbsq = rep(NA,bigmz); 
bigsigbsq[1:mz] = 0;                # on Z1 there is no reproduction 
bigsigbsq[(mz+1):(2*mz)] = lambda   # on Z2, flowering is certain, immediate RO is Poisson(lambda). 
bigsigbsq[(2*mz+1):bigmz] = sigbsq; # on Z3, the standard rules apply 

# make state-at-birth distribution on extended state space 
bigc0 = rep(0,bigmz); bigc0[1]=1; 

## Make conditional \\mathbf{P} on extended state space 
aM = apply (bigP[(mz+1):(2*mz), 1:mz], 2, sum)
Q = bigP[1:mz,1:mz]; 
N2 = solve(diag(mz)-Q); 

probEverBreed = c(aM %*% N2, rep(1, 2*mz))
pM = probEverBreed[1]; # all individuals are born into stage 1 

PCondBreed = matrix (NA, bigmz, bigmz)
  for (j in 1:bigmz) {
    PCondBreed[,j] = bigP[,j]*probEverBreed/probEverBreed[j]
  }

PplusB = cbind(PCondBreed, rep(0, bigmz))
deathProbB = 1 - apply (PplusB, 2, sum)
PplusB = rbind (PplusB, deathProbB)

## Fundamental matrix (N-tilde in the ms)
NB = solve (diag(bigmz) - PCondBreed)

###########################################################
# Compute mean and variance given state z 
###########################################################
expRCondZ = varRCondZ = rep (0, bigmz)

## expected lifetime reproduction conditional on state (rho1)
expRCondZ = apply (bigF %*% NB, 2, sum)
## rho1 includes RO for the state of being dead
rho1 = c(expRCondZ, 0)

## V (eq. 7)
VB = rho1^2 %*% PplusB - (rho1 %*% PplusB)^2

## All we really need is V-tilde.  Chop off the last entry.
VB = VB[-(bigmz+1)]

## There is no birth luck since everyone is born the same size.
## Var_0 (\\Ex [R | x, z_0]) = 0.

##################################################################################
# Loop over age to compute age-specific decomposition conditional on breeding
##################################################################################
stateTrajecLuckAgeB = fecLuckAgeB = rep(0, maxAge)
Pa = diag(1, bigmz); PaC0 = bigc0
stateTrajecLuckAge[1] = VB %*% PaC0
fecLuckAge[1] = bigsigbsq %*% PaC0

## a is actually age + 1, since the first year of life is age 0
for (a in 2:maxAge) {
  Pa = PCondBreed %*% Pa
  PaC0 = Pa %*% bigc0
  stateTrajecLuckAgeB[a] = VB %*% PaC0
  fecLuckAgeB[a] = bigsigbsq %*% PaC0
  if(a%%20==0) cat(a,"\\n"); 
}

ZtermsB=stateTrajecLuckAgeB; FtermsB=fecLuckAgeB; 
LuckByAge=Zterms+Fterms; 
LuckByAgeB = pM*(ZtermsB+FtermsB); # decomposition for breeders 
LuckByAgeNB = LuckByAge-LuckByAgeB; # decomposition for non-breeders 

Nc0 = NB %*% bigc0 
## This is also sum(F %*% N %*% c0), but this way is a little faster
expR = bigaveKids %*% Nc0

## record luck and pluck terms as functions of z 
stateTrajecLuckZB = VB * Nc0
fecLuckZB = bigsigbsq * Nc0

## For sanity check: variance of total LRO from now on, 
## conditional on state
rbarPib = expRCondZ %*% bigP     # \\bar{r} \\pi_b
r2 = (bigsigbsq + (bigaveKids)^2 + 2*bigaveKids*rbarPib) %*% NB
varRCondZ = r2 - expRCondZ^2

## Sanity check 
cat(sum(stateTrajecLuckAgeB + fecLuckAgeB),sum(stateTrajecLuckZB + fecLuckZB), varRCondZ[1],"\\n") 

################################################################# 
#  Plotting some results 
#################################################################
require(viridis); 
graphics.off(); dev.new(height=10.5,width=8); 

par(mfrow=c(3,2),bty="l",yaxs="i",xaxs="i",cex.axis=1.6,cex.lab=1.6,mgp=c(2.3,1,0),mar=c(4,4,2,1)); 

cols=plasma(100); 
matplot(cbind(LuckByAge[1:350],LuckByAgeB[1:350],LuckByAgeNB[1:350]),xlab="Age",ylab="Total and Conditional Luck", 
type="l",lwd=2,ylim=c(0,max(Zterms)),col=c("black","blue","red")); 
legend("topright",legend=c("Total","Within-Breeders","Breed or Not"),lwd=2,lty=c(1,2,3),col=c("black","blue","red"),bty="n",cex=1.2,inset=0.02); 

add_panel_label("a"); 

plot(Fterms[1:350],xlab="Age",ylab="Fecundity contribution", type="l",lwd=2,ylim=c(0,max(Fterms))); 
add_panel_label("b"); 

plot(allMesh,stateTrajecLuckZ,xlab="Size (dbh, cm)",ylab="State-trajectory contribution", type="l",lwd=2,ylim=c(0,max(stateTrajecLuckZ))); 
add_panel_label("c"); 
par(xpd=TRUE); 
points(allMesh[1:4],stateTrajecLuckZ[1:4],type="p",pch=16,col="red",cex=2); 
par(xpd=FALSE); 

plot(allMesh,fecLuckZ,xlab="Size (dbh, cm)",ylab="Fecundity contribution", type="l",lwd=2,ylim=c(0,max(fecLuckZ))); 
add_panel_label("d"); 

plot(log10(Survival[1:350]),xlab="Age",ylab="log10(Survival)",type="l",lwd=2,);  
add_panel_label("e"); 

image(1:350,allMesh,t(StageDist[,1:350])^0.25,xlab="Age",ylab="Size",col=plasma(100)); 
jmax=which((Zterms+Fterms)==max(Zterms+Fterms)); abline(v=jmax,col="red",lty=2); 
jmax=which((stateTrajecLuckZ + fecLuckZ)==max(stateTrajecLuckZ + fecLuckZ)); abline(h=allMesh[jmax],col="red",lty=2); 

add_panel_label("f"); 

dev.copy2pdf (file="figures/DacrydiumPartitioning.pdf")

