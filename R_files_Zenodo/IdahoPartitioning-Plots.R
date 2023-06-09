##############################################################################
# Uses the single-species IPM implemented by IdahoIPMFunctions.R 
# Run IdahoPartitioning.R first -- and see the header for that script. 
###############################################################################

rm(list=ls(all=TRUE)); 
graphics.off();

######### Edit as needed to set the working directory 
path ="/home/robin/ipm/drivers/Manuscripts/luckPartition/code"
out=try(setwd(path),silent=TRUE);
if(class(out)=="try-error") { 
    path=ifelse(.Platform$OS.type=="windows","c:/repos/drivers/Manuscripts/luckPartition/code","~/repos/drivers/Manuscripts/LuckPartition/code"); 
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

#################################################################### 
# Overall Luck and Pluck partition
####################################################################

load("ARTRpartitioning.Rdata"); 

# Luck = average of Var(R) across trait values
VarR = numeric(6); for (j in 1:6) VarR[j] = sum(stateTrajecLuckAge[,j]) + preNatal[j]; 
Luck_A = pmean(VarR,rep(1,6)/6); 

# Pluck = variance of E(R) across groups; here R is longevity;
ER = numeric(6); 
for(j in 1:6) {
    N = Nlist[[j]];
    Lifespan=apply(N,2,sum);
    ER[j]=pmean(Lifespan,c0); 
}
Pluck_A = pvar(ER,rep(1,6)/6);

LuckByAgeGroup_A = rbind(preNatal,stateTrajecLuckAge); 
LuckByAge_A = apply(LuckByAgeGroup_A, 1, mean); 
 
load("PSSPpartitioning.Rdata");  

VarR = numeric(6); for (j in 1:6) VarR[j] = sum(stateTrajecLuckAge[,j]) + preNatal[j]; 
Luck_P = pmean(VarR,rep(1,6)/6); 

# Pluck = variance of E(R) across groups; here R is longevity;
ER = numeric(6); 
for(j in 1:6) {
    N = Nlist[[j]];
    Lifespan=apply(N,2,sum);
    ER[j]=pmean(Lifespan,c0); 
}
Pluck_P = pvar(ER,rep(1,6)/6);

LuckByAgeGroup_P = rbind(preNatal,stateTrajecLuckAge); 
LuckByAge_P = apply(LuckByAgeGroup_P, 1, mean); 

graphics.off(); dev.new(width=8,height=5.5); 
par(bty="l",yaxs="i",xaxs="i",cex.axis=1.5,cex.lab=1.5,mgp=c(2.3,1,0),mar=c(4,4,2,1)); 

matplot((-1):60,cbind(LuckByAge_A[1:62],LuckByAge_P[1:62]),ylim=c(0,1.02*max(LuckByAge_P)),
xlab="Age",ylab="Contribution to Total Luck", type="l",lwd=2,col=c("black","red"),lty=c(2,1)); 

legend("top",legend=c("Artemisia","Pseudoroegneria"),lwd=2,lty=c(2,1),col=c("black","red"),bty="n",cex=1.3,inset=0.05); 

dev.copy2pdf(file="../figures/IdahoPartitioning.pdf"); 
