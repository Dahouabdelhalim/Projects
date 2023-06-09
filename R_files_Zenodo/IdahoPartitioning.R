##############################################################################
# Age partitioning of luck in lifespan for A. tripartita or P. spicata. 
# Uses the single-species IPM implemented by IdahoIPMFunctions.R 
#
# This file does the calculations for one of the two species. To replicate the
# plots in the paper, run this file and save the results as an .Rdata file, 
# either ARTRpartitioning.Rdata or PSSPpartitioning.Rdata. With both of those
# made, run IdahoPartitioning-Plots.R.
###############################################################################

rm(list=ls(all=TRUE)); 
graphics.off();

# Edit as needed to set working directory 
path ="/home/robin/ipm/drivers/Manuscripts/luckPartition/"
out=try(setwd(path),silent=TRUE);
if(class(out)=="try-error") { 
    path=ifelse(.Platform$OS.type=="windows","c:/repos/drivers/Manuscripts/luckPartition","~/repos/drivers/Manuscripts/LuckPartition"); 
    setwd(path); 
  }

source("Utilities.R"); 
source("IdahoIPMFunctions.R"); 
source("Standard Graphical Pars.R") 

sppList<-sort(c("ARTR","HECO","POSE","PSSP"));

species = 1 #ARTR 
pars = params[species, ]; names(pars)<-names(params); 
AGpars = AdultSurvGparams[species,]; SGpars=SeedSurvGparams[species,] 
Lx = -3; Ux = 9; mx = round(5*(Ux-Lx));  
LW = -7; UW = 5; mW = round(5*(UW - LW));
maxAge = 300; 


#species = 4; # PSSP
#pars = params[species, ]; names(pars)<-names(params); 
#AGpars = AdultSurvGparams[species,]; SGpars=SeedSurvGparams[species,] 
#Lx = -3; Ux = 7; mx = round(5*(Ux-Lx));  
#LW = -6; UW = 3; mW = round(5*(UW - LW));
#maxAge = 250; 

cat(b_x(Ux,pars)," should equal 1", "\\n"); # better equal 1 

############ make all the iteration matrices and fundamental operators 
bigPlist = Klist =  Nlist = list(6); 
for(Group in 1:6) {
    out = mk_K(mx=mx, mW=mW, Lx=Lx, Ux=Ux, LW = LW, UW=UW, pars=pars,Group)
    Klist[[Group]]=out; 
    P1 = out$P1; P2 = out$P2; xF2 = out$F2; 
    cat ("Finished making P1, P2, Fx", Group, "\\n")
    bigP = cbind(P1,P2); 
    bigP = rbind( matrix(0,ncol(P1),ncol(bigP)), bigP); 
    mz = nrow(bigP); 
    N = solve (diag(mz) - bigP)
    bigPlist[[Group]] = bigP; 
    Nlist[[Group]] = N;
}

# Discrete probability distributions for initial size as Non-seedling, and competition 
# Same for all groups 
yx = Klist[[1]]$yx; yW = Klist[[1]]$yW; 
pZ = c2_x1(yx,pars); pZ = pZ/sum(pZ);   # size 
pW = c1_w1(yW,pars); pW = pW/sum(pW);   # competition 

## F matrix. Fictitious: everyone has one offspring, of type 1. 
## By this device, formulas for LRO are instead giving result for lifespan. 
## This is the same of plants in all quadrat groups. 
pr = aveKids = rep(1,mz); 
F = matrix (0, mz, mz); F[1,]=aveKids; 

## Variance in number of kids. Fictitious, same for all Groups 
sigbsq = rep(0,mz); 

## Birth state distribution: seedlings with different W values. Same for all Groups 
c0 = rep(0, mz); c0[1:mW]=pW; 
c0tilde=c(c0,0);

stateTrajecLuckAge = PluckAge = meanRAge = Survival = matrix(NA, maxAge,6)
preNatal=numeric(6); 

### Eqn. (14) says that we do the age-partition separately for each value of
### the trait x; here the trait is quadrat group. 
for(Group in 1:6) {
    out=Klist[[Group]];
    yx = out$yx; yW = out$yW; hx=out$hx; hW = out$hW; 
    bigP = bigPlist[[Group]]; 
    N = Nlist[[Group]]; 
    mz = nrow(bigP); 

#################### Expanded matrix including \\omega = dead.  
    Pplus = cbind (bigP, rep(0, mz))
    deathProb = 1 - apply (Pplus, 2, sum)
    Pplus = rbind (Pplus, deathProb)

###########################################################
# Compute mean and variance given state z 
###########################################################
    expRCondZ = varRCondZ = rep (0, mz);

## expected lifetime reproduction conditional on state (rho1)
    expRCondZ = apply (F %*% N, 2, sum)
## rho1 includes RO for the state of being dead
    rho1 = c(expRCondZ, 0)

## V, eqn. \\label{eqn:varbarr}
    V = rho1^2 %*% Pplus - (rho1 %*% Pplus)^2
## All we really need is V-tilde.  Chop off the last entry.
    V = V[-(mz+1)]

## Pre-natal luck Var_z0 (\\Ex [R | x, z_0])  = <c0 , rho1^2> - <c0, rho1>^2
    preNatal[Group] = sum(c0tilde*rho1^2) - (sum(c0tilde*rho1))^2 

##########################################################################
# Loop over age to compute the age-specific Luck terms (eqn. 7 in paper)
# Note, there is no fecundity luck in this model!  
##########################################################################
    Pa = diag(1, mz); PaC0 = c0;
    stateTrajecLuckAge[1,Group] = sum(V * PaC0)
    Survival[1,Group]=1; 

## a is actually age + 1, since the first year of life is age 0
    for (a in 2:maxAge) {
        Pa = bigP %*% Pa
        PaC0 = Pa %*% c0
        stateTrajecLuckAge[a,Group] = sum(V * PaC0) 
        meanRAge[a,Group] = sum(expRCondZ*PaC0); 
        Survival[a,Group]=sum(PaC0); 
        if(a%%10==0) cat(a,"\\n"); 
    }

## For sanity check: variance of total LRO from now on, conditional on state
    rbarPib = expRCondZ %*% bigP     # \\bar{r} \\pi_b
    r2 = (sigbsq + (aveKids)^2 + 2*aveKids*rbarPib) %*% N
    varRCondZ = r2 - expRCondZ^2;

## Unconditional variance of total LRO from now on
    EvarRCondZ = pmean(varRCondZ,c0); 
    VarexpRCondZ = pvar(expRCondZ,c0) 

## Sanity check 
    cat(sum(stateTrajecLuckAge[,Group]) + preNatal[Group], EvarRCondZ + VarexpRCondZ,"\\n") 

} 

#################################################################### 
# Overall Luck and Pluck partition
####################################################################

# Luck = average of Var(R) across trait values (eqn. 14 in the paper) 
VarR = numeric(6); for (j in 1:6) VarR[j] = sum(stateTrajecLuckAge[,j]) + preNatal[j]; 
Luck = pmean(VarR,rep(1,6)/6); 

# Pluck = variance of E(R) across groups (eqn. 13); here R is longevity
ER = numeric(6); 
for(j in 1:6) {
    N = Nlist[[j]];
    Lifespan=apply(N,2,sum);
    ER[j]=pmean(Lifespan,c0); 
}
Pluck = pvar(ER,rep(1,6)/6);     

