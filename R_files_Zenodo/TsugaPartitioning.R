##################################################################################
# Luck by Age and Stage, using the Tsuga model from van Daalen and Caswell (2017),
# originally from COMPADRE. 
#
# Luck by Age is sub-partitioned by Breeders versus non-Breeders 
##################################################################################

rm(list=ls(all=TRUE))

graphics.off(); require(viridis); 

############ CHANGE AS NEEDED to set working directory ################
root ="/home/robin/ipm/drivers"; out=try(setwd(root),silent=TRUE);
if(class(out)=="try-error") { 
    root=ifelse(.Platform$OS.type=="windows","c:/repos/drivers","~/ipm/drivers/"); 
    setwd(root); 
  }
setwd("Manuscripts/Luckpartition"); 
########################################################################

source("Utilities.R"); 

add_panel_label <- function(ltype="a"){
    text <- paste(LETTERS[letters==ltype], ")", sep="")
    mtext(text=text, side=3, adj=0)
}

######################################################################
# Hemlock model from van Daalen and Caswell (2017) 
######################################################################

# Load the matrix from COMPADRE ######################################
# Note: the matrix displayed by van Daalen & Caswell paper is rounded 

################# 
# This commented-out code shows how the projection matrix was extracted
# from COMPADRE but for convenience, we have stored it in an .Rdata file 
## require(Mage)
## load("c:/repos/renew/COMPADRE_v.5.0.0.RData")
## D2 <- subsetDB(compadre,MatrixComposite == "Individual"&MatrixTreatment == "Unmanipulated")
## index <- which(D2$metadata$SpeciesAuthor=="Tsuga_canadensis")

## tsuga= list(length(index)) 
## for(j in 1:length(index)) {show=index[j]; tsuga[[j]]=D2$mat[[show]]$matA;}
## A = 0.5*(tsuga[[1]]+tsuga[[2]]);

load ("TsugaMatrix.Rdata")

ev=eigen(A); 
lambda=Re(ev$values[1]); w = ev$vectors[,1]; w=Re(w); w=w/sum(w); 
v=eigen(t(A))$vectors[,1];v=Re(v)/Re(v[1]); 
sens=outer(v,w)/sum(v*w);
elas=A*sens/lambda;
matrix.image(elas); 

##################################################################### 
# Make the rest of the matrices need for the analysis. 
# We use here the van Daalen & Caswell notation in which A = U + F
# is the projection matrix on living states, while P (which we call
# P+ in the paper) is the projection matrix including \\omega = dead.  
#####################################################################
U=A; U[1,-1]=0; # remove fecundity entries 

P=rbind(U,1-apply(U,2,sum)); s=nrow(P); 
P=cbind(P,rep(0,s)); P[s,s]=1; # this is what we call P+ 

fT = matrix(A[1,],1,s-1); fT[1]=0; fT=cbind(fT,0); 

vec1s = matrix(1,s,1); 
R1 = (vec1s%*%fT); 
R2 = R1 + R1*R1; # formula for Poisson number of offspring 

################################################################
# Do the calculations for total LRO variance using the formulas
# in van Daalen & Caswell 
################################################################
tau = s-1; alpha = 1; 
Z = cbind(diag(tau),matrix(0,tau,alpha)); 
N = solve(diag(tau)-U); 

tildeRho1 = t(N) %*% Z %*% t(P*R1) %*% vec1s; #check, same as their Table 3
  
tildeR1 = Z %*% R1 %*% t(Z) 
tildeRho2 = t(N) %*% ( Z %*% t(P*R2) %*% vec1s + 2*t(U*tildeR1)%*%tildeRho1) #check, same as their Table 3 

#########################################################
# Decompose by age, our way: eqn. (7)  
#########################################################

# first term is zero: Var_{z0}(E(R|z0) = 0 because everyone is born the same.    

# the next bunch of terms, involving state trajectory variance  
Rho1 = c(tildeRho1,0);          # expand E(R) to include \\omega 
V = (Rho1^2) %*% P - (Rho1 %*% P)^2; 
c0 = rep(0,s); c0[1]=1; 

Amax=2000; 
Zterms=Fterms = Survival = MeanStage = numeric(Amax); 
StageDist = matrix(0,tau,Amax); 
x = c0/sum(c0); # first term in the sum 
for(k in 1:Amax) {
    Zterms[k]= V%*%x; 
    Fterms[k]= fT%*%x; # variance=mean for Poisson 
    xlive=x[-length(x)]; 
    Survival[k]=sum(xlive); 
    plive = xlive/sum(xlive); 
    MeanStage[k]= sum(plive*c(1:length(plive))); 
    x = P %*% x; # move on to the next power of age a
    StageDist[,k]=plive; 
}  
# and now, does it work? The sum over ages should equal the whole 
sum(Zterms+Fterms); tildeRho2[1]-tildeRho1[1]^2; 

######## Plot results by age only 
if(FALSE) {
graphics.off(); dev.new(height=7,width=8); 
par(mfrow=c(2,2),bty="l",yaxs="i",xaxs="i",cex.axis=1.3,cex.lab=1.3,mgp=c(2.2,1,0),mar=c(4,4,2,2)); 
jmax=which((Zterms+Fterms)==max(Zterms+Fterms)); 

plot(Zterms[1:800],xlab="Age",ylab="State-trajectory contribution", type="l",lwd=2,ylim=c(0,max(Zterms))); 
abline(v=jmax,col="red",lty=2); add_panel_label("a"); 

plot(Fterms[1:800],xlab="Age",ylab="Fecundity contribution", type="l",lwd=2,ylim=c(0,max(Fterms))); 
abline(v=jmax,col="red",lty=2); add_panel_label("b"); 

plot(log(Survival[1:800]),xlab="Age",ylab="log(Survival)",type="l",lwd=2,);  
abline(v=jmax,col="red",lty=2); abline(v=jmax,col="red",lty=2); add_panel_label("c"); 

image(1:800,1:6,t(StageDist[,1:800]),xlab="Age",ylab="Stage",col=plasma(6)); abline(v=jmax,col="red",lty=2); 
add_panel_label("d"); 
dev.copy2pdf(file="TsugaLuckByAge.pdf");
}

#########################################################
# Decompose by stage and plot  (eqn. 10) 
#########################################################
Zfunc = V[1:6]*(N%*%c0[1:6]); 
Ffunc = fT[1:6]*(N%*%c0[1:6]);
graphics.off(); dev.new(height=7,width=8); 
par(mfrow=c(2,2),bty="l",mar=c(4,5,2,0),mgp=c(2.2,1,0),cex.axis=1.4,cex.lab=1.4); 
barplot(as.vector(Zfunc),xlab="",ylab="Contribution to Var(LRO)",names.arg=as.character(1:6),main="State trajectory-driven luck")
add_panel_label("a"); 
barplot(as.vector(Ffunc),xlab="",ylab="Contribution to Var(LRO)",names.arg=as.character(1:6),main="Fecundity-driven luck")
add_panel_label("b"); 
barplot(as.vector(V[1:6]),xlab="Stage",ylab="Contribution per visit",names.arg=as.character(1:6))
add_panel_label("c"); 
barplot(as.vector(fT[1:6]),xlab="Stage",ylab="Contribution per visit",names.arg=as.character(1:6))
add_panel_label("d"); 
# dev.copy2pdf(file="TsugaLuckByStage.pdf");

################################################################################
# Decompose by age for breeders only 
################################################################################

# Make the matrices conditional on breeding #########################
# Here S={1,2}, M={3,4,5,6}. As there is no shrinkage from M to S
# we can take Z_1=S, Z_2=M. See Appendix S4 for the underlying math  

Q1 = U[1:2,1:2]; N1 = solve(diag(2)-Q1); 
aM = apply(P[3:6,1:2],2,sum); 
probEverBreed = aM %*% N1; probEverBreed=c(probEverBreed,rep(1,4)); # chance to become a breeder, given state

pM = probEverBreed[1]; # all individuals are born into stage 1 

# Make the U conditional on breeding 
PCondBreed = matrix(NA,6,6); 
for(j in 1:6) {PCondBreed[,j]=U[,j]*probEverBreed/probEverBreed[j]} # eqn. S22

# Make the conditional P+ 
PB=rbind(PCondBreed,1-apply(PCondBreed,2,sum)); s=nrow(PB); PB=cbind(PB,rep(0,s)); PB[s,s]=1; 

### These have been computed already, but just to be safe do it again 
fT = matrix(A[1,],1,s-1); fT[1]=0; fT=cbind(fT,0); 
vec1s = matrix(1,s,1); R1 = (vec1s%*%fT); 
R2 = R1 + R1*R1; # Poisson number of offspring 

# Do the calculations for total LRO variance ###########################
tau = s-1; alpha = 1; 
Z = cbind(diag(tau),matrix(0,tau,alpha)); 
NB = solve(diag(tau)-PCondBreed); 

tildeRho1 = t(NB) %*% Z %*% t(PB*R1) %*% vec1s; 
tildeR1 = Z %*% R1 %*% t(Z) 
tildeRho2 = t(NB) %*% ( Z %*% t(PB*R2) %*% vec1s + 2*t(PCondBreed*tildeR1)%*%tildeRho1) 

# first term is zero: Var_{z0}(E(R|z0) = 0 because everyone is born the same.    
# the next bunch of terms, involving state trajectory variance  
Rho1 = c(tildeRho1,0);              # expand E(R) to include \\omega 
V = (Rho1^2) %*% PB - (Rho1 %*% PB)^2; 
c0 = rep(0,s); c0[1]=1; 

Amax=2000; 
ZtermsB = FtermsB = numeric(Amax); 
x = c0/sum(c0); 
for(k in 1:Amax) {
    ZtermsB[k]= V%*%x; 
    FtermsB[k]= fT%*%x; # variance=mean for Poisson 
    x = PB %*% x;
}  
LuckByAgeB = pM*(ZtermsB+FtermsB); # decomposition for breeders 
LuckByAgeNB = (Zterms+Fterms)-(LuckByAgeB); # decomposition for non-breeders 

sum(LuckByAgeB)/sum(Zterms+Fterms);

######################################################
# Manage for Luck, vs. eigenvalue elasticity 
######################################################
cbind(V[1:6], tildeRho1); # highest V in 5 and 6; help 5 survive & grow, help 6 live. 

ev=eigen(A); 
lambda=Re(ev$values[1]); w = ev$vectors[,1]; w=Re(w); w=w/sum(w); 
v=eigen(t(A))$vectors[,1];v=Re(v)/Re(v[1]); 
sens=outer(v,w)/sum(v*w);
elas=A*sens/lambda;
# Highest elasticity is 6 -> 6, second (2/3 as large) is 4 -> 4. 
# This reflects that w[5] < w[4]/3, elasticity is based on affecting all individuals 

######################################################
#  Plots of results for the manuscript 
######################################################

dev.new(height=10.5,width=8); 
par(mfrow=c(3,2),bty="l",yaxs="i",xaxs="i",cex.axis=1.4,cex.lab=1.5, mgp=c(2.3,1,0),mar=c(4,4,2,2)); 

jmax=which((Zterms+Fterms)==max(Zterms+Fterms)); 

plot(Zterms[1:800],xlab="Age",ylab="State-trajectory contribution", type="l",lwd=2,ylim=c(0,max(Zterms))); 
abline(v=jmax,col="red",lty=2); add_panel_label("a"); 

plot(Fterms[1:800],xlab="Age",ylab="Fecundity contribution", type="l",lwd=2,ylim=c(0,max(Fterms))); 
abline(v=jmax,col="red",lty=2); add_panel_label("b"); 

py = cbind(LuckByAgeB[1:800],LuckByAgeNB[1:800]); 
matplot(py,xlab="Age",ylab="Conditional Contribution",col="black",type="l",lwd=2,ylim=c(0,1.02*max(py)));  
legend("topright",legend=c("Among Breeders (88.4%)","Breed or Not (11.6%)"),lwd=2,lty=c(1,2),col="black",bty="n",cex=1.2,inset=0.02); 
abline(v=jmax,col="red",lty=2); 
add_panel_label("c"); 

image(1:800,1:6,t(StageDist[,1:800]),xlab="Age",ylab="Stage",col=plasma(6)); abline(v=jmax,col="red",lty=2); 
add_panel_label("d"); 

barplot(as.vector(Zfunc),xlab="Stage",ylab="State-trajectory contribution",names.arg=as.character(1:6))
add_panel_label("e"); 
barplot(as.vector(Ffunc),xlab="Stage",ylab="Fecundity contribution",names.arg=as.character(1:6))
add_panel_label("f"); 
dev.copy2pdf(file="figures/TsugaLuckByAgeStage.pdf");


