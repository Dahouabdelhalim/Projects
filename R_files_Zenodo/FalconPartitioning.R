##################################################################################
# Luck by Age and Stage, using the Peregrine Falcon model from Altwegg et al. (2014)
# Luck by Age is sub-partitioned by Breeders versus non-Breeders 
##################################################################################

rm(list=ls(all=TRUE))
graphics.off(); 

############ CHANGE AS NEEDED to set working directory 
root ="/home/robin/ipm/drivers"; out=try(setwd(root),silent=TRUE);
if(class(out)=="try-error") { 
    root=ifelse(.Platform$OS.type=="windows","c:/repos/drivers","~/repos/drivers"); 
    setwd(root); 
  }
setwd("Manuscripts/luckPartition"); 

add_panel_label <- function(ltype="a"){
    text <- paste(LETTERS[letters==ltype], ")", sep="")
    mtext(text=text, side=3, adj=0)
}

######################################################################
# Original code to load Peregrine Falcon model from COMADRE
######################################################################
# load("COMADRE_v.3.0.0.RData")
# D2 <- subsetDB(comadre,MatrixComposite == "Individual"&MatrixTreatment == "Unmanipulated")
# Birds <- subsetDB(D2,OrganismType=="Aves")
# Falcon<-subsetDB(Birds,CommonName=="Peregrine falcon"); 
# nyrs=nrow(Falcon$metadata);    
# U = Falcon$mat[[1]]$matU
# F = Falcon$mat[[1]]$matF[1,]
# for(j in 2:nyrs) {
#     F=F+Falcon$mat[[j]]$matF[1,]
#     U=U+Falcon$mat[[j]]$matU
# }
# FalconF=F/nyrs; FalconU=U/nyrs; 
# phi2=apply(FalconU,2,sum)[5];

######################################################################
# Load Peregrine Falcon model, saved in .Rdata file
######################################################################
load ("FalconMatrix.Rdata")
U = FalconU; 

##################################################################### 
# Make the rest of the matrices 
# We use the van Daalen & Caswell notation in which A = U + F
# is the projection matrix on living states, while P (which we call
# P+ in the paper) is the projection matrix including \\omega = dead.  
#####################################################################

P=rbind(U,1-apply(U,2,sum)); s=nrow(P); P=cbind(P,rep(0,s)); P[s,s]=1; # this is what we call P+ 

# Convert fecundities to a pre-breeding census. Only stage 5 breeds. 
fT=matrix(0,1,s-1); fT[s-1]=FalconF[s-1]/phi2; fT=cbind(fT,0); 

vec1s = matrix(1,s,1); 
R1 = (vec1s%*%fT); 
R2 = R1 + R1*R1; # Poisson number of offspring 

# Do the calculations for total LRO variance using the formulas from
# van Daalen & Caswell 
tau = s-1; alpha = 1; 
Z = cbind(diag(tau),matrix(0,tau,alpha)); 
N = solve(diag(tau)-U); 

tildeRho1 = t(N) %*% Z %*% t(P*R1) %*% vec1s; #check, same as their Table 3
tildeR1 = Z %*% R1 %*% t(Z) 
tildeRho2 = t(N) %*% ( Z %*% t(P*R2) %*% vec1s + 2*t(U*tildeR1)%*%tildeRho1) #check, same as their Table 3 

#########################################################
# Decompose by age, our way (eqn. 7 in the paper) 
#########################################################

# first term is zero: Var_{z0}(E(R|z0) = 0 because everyone is born the same.    

# the next bunch of terms, involving state trajectory variance  
Rho1 = c(tildeRho1,0);          # expand E(R) to include \\omega 
V = (Rho1^2) %*% P - (Rho1 %*% P)^2; 
c0 = rep(0,s); c0[1]=1; 

Amax=100; 
Zterms=Fterms = Survival = MeanStage = numeric(Amax); 
StageDist = matrix(0,tau,Amax); 
x = c0/sum(c0); 
for(k in 1:Amax) {
    Zterms[k]= V%*%x; 
    Fterms[k]= fT%*%x; # variance=mean for Poisson 
    xlive=x[-length(x)]; 
    Survival[k]=sum(xlive); 
    plive = xlive/sum(xlive); 
    MeanStage[k]= sum(plive*c(1:length(plive))); 
    x = P %*% x; # next higher power of P^a in the sum 
    StageDist[,k]=plive; 
}  
# and now, does it work? Sum over age should equal the total 
sum(Zterms+Fterms); tildeRho2[1]-tildeRho1[1]^2; 

######################################################
# Manage for Luck, vs. eigenvalue elasticity 
######################################################
cbind(V[1:5], tildeRho1); 
#[1,] 0.9258453 0.5185959
#[2,] 2.1742114 2.3038882  <- highest V 
#[3,] 1.4968500 2.8956259
#[4,] 1.5590376 2.9958305
#[5,] 1.5590376 3.5162330

# stage 2 is non-breeder age 1, can go to 3 (non-breeder age 3),
# or 5 (breeder) or die. 

A = FalconU; A[1,1:5]=fT[1:5];  
ev=eigen(A); 
lambda=Re(ev$values[1]); w = ev$vectors[,1]; w=Re(w); w=w/sum(w); 
v=eigen(t(A))$vectors[,1];v=Re(v)/Re(v[1]); 
sens=outer(v,w)/sum(v*w);
elas=A*sens/lambda;
# Highest elasticity is 5 -> 5, by far, e[5,5]= 0.73, more than 2 x everything else together.   
# This reflects that w[5] = 57% of population, 10% in 2 through 4.  


######## Some plots of the results 
require(viridis); 
graphics.off(); dev.new(height=7,width=8); 
par(mfrow=c(2,2),bty="l",yaxs="i",xaxs="i",cex.axis=1.3,cex.lab=1.3,mgp=c(2.2,1,0),mar=c(4,4,2,2)); 
jmax=which((Zterms+Fterms)==max(Zterms+Fterms)); 

plot(Zterms[1:30],xlab="Age",ylab="State-trajectory contribution", type="l",lwd=2,ylim=c(0,max(Zterms))); 
abline(v=jmax,col="red",lty=2); add_panel_label("a"); 

plot(Fterms[1:30],xlab="Age",ylab="Fecundity contribution", type="l",lwd=2,ylim=c(0,max(Fterms))); 
abline(v=jmax,col="red",lty=2); add_panel_label("b"); 

plot(log(Survival[1:30]),xlab="Age",ylab="log(Survival)",type="l",lwd=2,);  
abline(v=jmax,col="red",lty=2); abline(v=jmax,col="red",lty=2); add_panel_label("c"); 

image(1:30,1:5,t(StageDist[,1:30]),xlab="Age",ylab="Stage",col=plasma(6)); abline(v=jmax,col="red",lty=2); 
add_panel_label("d"); 

#########################################################
# Decompose by stage 
#########################################################
Zfunc = V[1:5]*(N%*%c0[1:5]); 
Ffunc = fT[1:5]*(N%*%c0[1:5]);
dev.new(height=7,width=8); 
par(mfrow=c(2,2),bty="l",mar=c(4,5,2,0),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4); 
barplot(as.vector(Zfunc),xlab="",ylab="Contribution to Var(LRO)",names.arg=as.character(1:5),main="State trajectory-driven luck")
add_panel_label("a"); 
barplot(as.vector(Ffunc),xlab="",ylab="Contribution to Var(LRO)",names.arg=as.character(1:5),main="Fecundity-driven luck")
add_panel_label("b"); 
barplot(as.vector(V[1:5]),xlab="Stage",ylab="Contribution per visit",names.arg=as.character(1:5))
add_panel_label("c"); 
barplot(as.vector(fT[1:5]),xlab="Stage",ylab="Contribution per visit",names.arg=as.character(1:5))
add_panel_label("d"); 

################################################################################
# Decompose by age for breeders only 
################################################################################

# Make the matrices conditional on breeding; see Appendix S4 for the underlying math 
Q = P[1:4,1:4]; 
NS = solve(diag(4)-Q); 
aB = P[5,1:4]; 
B = aB%*%NS; B=c(B,rep(1,1)); # chance to become a breeder given state
pM = B[1]; # all individuals are born into stage 1 
UB = U; for(z in 1:5) {UB[,z]=B*U[,z]/B[z]} # this is the U conditional on breeding 
PB=rbind(UB,1-apply(UB,2,sum)); s=nrow(PB); PB=cbind(PB,rep(0,s)); PB[s,s]=1; # this is the conditional P+ 

# Do the calculations for total LRO variance ###########################
tau = s-1; alpha = 1; 
Z = cbind(diag(tau),matrix(0,tau,alpha)); 
NB = solve(diag(tau)-UB); 

tildeRho1 = t(NB) %*% Z %*% t(PB*R1) %*% vec1s; 
tildeR1 = Z %*% R1 %*% t(Z) 
tildeRho2 = t(NB) %*% ( Z %*% t(PB*R2) %*% vec1s + 2*t(UB*tildeR1)%*%tildeRho1) 

# first term is zero: Var_{z0}(E(R|z0) = 0 because everyone is born the same.    
# the next bunch of terms, involving state trajectory variance  
Rho1 = c(tildeRho1,0);              # expand E(R) to include \\omega 
V = (Rho1^2) %*% PB - (Rho1 %*% PB)^2; 
c0 = rep(0,s); c0[1]=1; 

Amax=100; 
ZtermsB = FtermsB = numeric(Amax); 
x = c0/sum(c0); 
for(k in 1:Amax) {
    ZtermsB[k]= V%*%x; 
    FtermsB[k]= fT%*%x; # variance=mean for Poisson distribution of offspring 
    x = PB %*% x;
}  
LuckByAgeB = pM*(ZtermsB+FtermsB); # decomposition for breeders 
LuckByAgeNB = (Zterms+Fterms)-(LuckByAgeB); # decomposition for non-breeders 

sum(LuckByAgeB)/sum(Zterms+Fterms);

######## Plots of results for the manuscript 
graphics.off(); dev.new(height=7,width=8); 
par(mfrow=c(2,2),bty="l",yaxs="i",xaxs="i",cex.axis=1.2,cex.lab=1.3, mgp=c(2.3,1,0),mar=c(4,4,2,2)); 

jmax=which((Zterms+Fterms)==max(Zterms+Fterms)); 

matplot(1:30,cbind(Zterms[1:30],Fterms[1:30]) ,xlab="Age",ylab="Contribution to Var(R)", type="l",lty=c(1,2), 
lwd=2,ylim=c(0,max(Zterms)),col="black"); 
add_panel_label("a"); 
legend("topright",legend=c("State-trajectory Luck","Fecundity Luck"),lwd=2,lty=c(1,2),col="black",bty="n",cex=1.2,inset=0.02); 

py = cbind(LuckByAgeB,LuckByAgeNB);
fracs=apply(py,2,sum); fracs=fracs/sum(fracs); fracs; 
 
matplot(1:30,py[1:30,],xlab="Age",ylab="Conditional Contribution",col="black",type="l",lwd=2,ylim=c(0,1.02*max(py)));  
legend("topright",legend=c("Among-Breeders (57.1%)","Breed or Not (42.9%)"),lwd=2,lty=c(1,2),col="black",bty="n",cex=1.2,inset=0.02); 
add_panel_label("b"); 

image(1:25,1:5,t(StageDist[,1:25]),xlab="Age",ylab="Stage",col=plasma(20)); 
add_panel_label("c"); 

barplot(as.vector(Zfunc),xlab="Stage",ylab="State-trajectory contribution",names.arg=as.character(1:5))
add_panel_label("d"); 

dev.copy2pdf(file="figures/FalconLuckByAgeStage.pdf");







