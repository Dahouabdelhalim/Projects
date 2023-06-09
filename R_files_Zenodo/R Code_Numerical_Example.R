###################################################
### R code for Numerical example
###  2022-11-22
###################################################

#rm(list=ls())

#BtheB data from HSAUR2 package
library(HSAUR2)

#BtheB Dataset
data(BtheB)

#For BtheB group (treatment)
#Baseline
X_T=BtheB$bdi.pre[BtheB$treatment=='BtheB']
#Post-treatment measurements
Y1_T=BtheB$bdi.2m[BtheB$treatment=='BtheB']
Y2_T=BtheB$bdi.3m[BtheB$treatment=='BtheB']
Y3_T=BtheB$bdi.5m[BtheB$treatment=='BtheB']
Y4_T=BtheB$bdi.8m[BtheB$treatment=='BtheB']

#For TAU group (control)
#Baseline
X_C=BtheB$bdi.pre[BtheB$treatment=='TAU']
#Post-treatment measurements
Y1_C=BtheB$bdi.2m[BtheB$treatment=='TAU']
Y2_C=BtheB$bdi.3m[BtheB$treatment=='TAU']
Y3_C=BtheB$bdi.5m[BtheB$treatment=='TAU']
Y4_C=BtheB$bdi.8m[BtheB$treatment=='TAU']

#Assume treatment and control group have the same variance and correlation
sigmaX=sd(c(X_C,X_T))
sigmaY=sqrt((var(c(Y1_C,Y1_T),na.rm=TRUE)+var(c(Y2_C,Y2_T),na.rm=TRUE)+
var(c(Y3_C,Y3_T),na.rm=TRUE)+var(c(Y4_C,Y4_T),na.rm=TRUE))/4)

sigmaX^2 #117.5163
sigmaY^2 #116.7616

#correlation between post-post
rhoY=(cor(c(Y1_C,Y1_T),c(Y2_C,Y2_T),use = "pairwise.complete.obs")+
cor(c(Y1_C,Y1_T),c(Y3_C,Y3_T),use = "pairwise.complete.obs")+
cor(c(Y1_C,Y1_T),c(Y4_C,Y4_T),use = "pairwise.complete.obs")+
cor(c(Y2_C,Y2_T),c(Y3_C,Y3_T),use = "pairwise.complete.obs")+
cor(c(Y2_C,Y2_T),c(Y4_C,Y4_T),use = "pairwise.complete.obs")+
cor(c(Y3_C,Y3_T),c(Y4_C,Y4_T),use = "pairwise.complete.obs"))/6

rhoY #0.7714201

#correlation between pre-post
rhoXY=(cor(c(X_C,X_T),c(Y1_C,Y1_T),use = "pairwise.complete.obs")+
cor(c(X_C,X_T),c(Y2_C,Y2_T),use = "pairwise.complete.obs")+
cor(c(X_C,X_T),c(Y3_C,Y3_T),use = "pairwise.complete.obs")+
cor(c(X_C,X_T),c(Y4_C,Y4_T),use = "pairwise.complete.obs"))/4

rhoXY #0.5186458

#Assume correlation between pre-pre equals rhoY
rhoX=rhoY

#The treatment effect
delta=mean(c(Y1_C, Y2_C,Y3_C,Y4_C),na.rm=TRUE)-mean(c(Y1_T, Y2_T,Y3_T,Y4_T),na.rm=TRUE)
delta #5.373436

#For given alpha=0.05, 1-beta=0.8, the sample size per group
n=function(S,T)
{
   2*(qnorm(1-0.05/2)+qnorm(0.8))^2*sigmaY^2/delta^2*((1+rhoX*(S-1))*(1+rhoY*(T-1))-rhoXY^2*S*T)/(T*(1+rhoX*(S-1)))
}

round(n(1,1),digits=0) #46
round(n(2,1),digits=0) #44
round(n(1,2),digits=0) #39
round(n(1,4),digits=0) #36
round(n(4,1),digits=0) #43
round(n(2,3),digits=0) #35
round(n(3,2),digits=0) #36
round(n(2,4),digits=0) #33
round(n(4,2),digits=0) #36

1-round(n(2,1),digits=0)/round(n(1,1),digits=0) #required sample size decreases 4.3% comparing with $S=1, T=1$ design
1-round(n(1,2),digits=0)/round(n(1,1),digits=0) #required sample size decreases 15.2% comparing with $S=1, T=1$ design
1-round(n(1,4),digits=0)/round(n(1,1),digits=0) #required sample size decreases 21.7% comparing with $S=1, T=1$ design
1-round(n(4,1),digits=0)/round(n(1,1),digits=0) #required sample size decreases 6.5% comparing with $S=1, T=1$ design
1-round(n(2,3),digits=0)/round(n(1,1),digits=0) #required sample size decreases 23.9% comparing with $S=1, T=1$ design
1-round(n(3,2),digits=0)/round(n(1,1),digits=0) #required sample size decreases 21.7% comparing with $S=1, T=1$ design
1-round(n(2,4),digits=0)/round(n(1,1),digits=0) #required sample size decreases 28.3% comparing with $S=1, T=1$ design
1-round(n(4,2),digits=0)/round(n(1,1),digits=0) #required sample size decreases 21.7% comparing with $S=1, T=1$ design

###################################################
######### Optimization function of Sopt
###################################################

rhoX*rhoY-rhoXY^2  #>0
sqrt((1-rhoY)/((1-rhoX)*rhoXY^2) )+1  #2.928098<M

Sopt=function(M)
{
   (M*rhoXY^2*(rhoX-1)-rhoX*(1-rhoX)*(1-rhoY)+rhoXY*(M*rhoX-rhoX+1)*sqrt((1-rhoX)*(1-rhoY)))/
      (rhoX*(rhoXY^2-rhoX*rhoY)-rhoXY^2+rhoX^2)
}

Sopt(M=5)  #1.832968

round(n(1,4),digits=0) #36
round(n(2,3),digits=0) #35