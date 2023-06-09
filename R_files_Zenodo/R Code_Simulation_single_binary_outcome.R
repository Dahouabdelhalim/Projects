######################################################################################
### R Code Simulation Study of Single Binary Outcomes
### 2022-11-25
### R version 4.1.3
### dichotomized probability u=0.4
######################################################################################

#rm(list=ls())
setwd("~/shiyangm")
library(copula)
library(mvtnorm)

#Small sample n0=n1=50, only run models 1 and 2
power_RB=function(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0,u=0.4,n0=50,n1=50){
  
  #u is the success probability (P(Y=1)) for the control group
  mu_C=rep(u,3)
  mu_T=rep(NA,3)
  mu_T[1:2]=rep(u,2)
  mu_T[3]=exp(beta+log(u/(1-u)))/(1+exp(beta+log(u/(1-u))))  
  
  cor_Y_X1=rep(NA, nsim)
  cor_Y_X2=rep(NA, nsim)
  cor_X1_X2=rep(NA, nsim)
  
  pvalue1=rep(NA, nsim)
  pvalue2=rep(NA, nsim)
  
  set.seed(12345)
  
  for (i in 1: nsim){
    
    #Generate the data  
    myCop=normalCopula(param = c(rhoX, rhoXY, rhoXY),dim = 3, dispstr = "un")
    U=rCopula(n0+n1,myCop)
    
    #Generate the correlated binary outcomes from the Gaussian Coupla
    X1=c(U[1:n0,1]<=mu_C[1],U[(n0+1):(n0+n1),1]<=mu_T[1])
    X2=c(U[1:n0,2]<=mu_C[2],U[(n0+1):(n0+n1),2]<=mu_T[2])
    Y=c(U[1:n0,3]<=mu_C[3],U[(n0+1):(n0+n1),3]<=mu_T[3])
    
    cor_Y_X1[i]=cor(Y,X1)
    cor_Y_X2[i]=cor(Y,X2)
    cor_X1_X2[i]=cor(X1,X2)
    
    X=X1+X2
    X_log=log((X+1/2)/(2-X+1/2))
    
    Treat=c(rep('Control',n0),rep('Treatment',n1))
    
    fit1=glm(Y ~ Treat+X2,family=binomial)
    pvalue1[i]=coef(summary(fit1))[2,4]
    
    fit2=glm(Y ~ Treat+X_log,family=binomial())
    pvalue2[i]=coef(summary(fit2))[2,4]
    
  }
  
  power=c(sum(pvalue1<=0.05)/nsim,
          sum(pvalue2<=0.05)/nsim)
  
  return(power)  
}

#n0=n1=50
TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n50=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0,u=0.4,n0=50,n1=50)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n50=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=50,n1=50)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n50=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=50,n1=50)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n50=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=50,n1=50)

TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n50=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=50,n1=50)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n50=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=50,n1=50)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n50=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=50,n1=50)

TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n50=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=50,n1=50)
TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n50=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=50,n1=50)

TypeI_u_0.4_n50=rbind(TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n50,
                      TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n50,
                      TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n50,
                      TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n50,
                      TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n50,
                      TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n50,
                      TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n50,
                      TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n50,
                      TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n50)

TypeI_u_0.4_n50
saveRDS(TypeI_u_0.4_n50, file="TypeI_u_0.4_n50.Rds")

Power_u_0.4_rhoXY_0.5_rhoX_0.6_n50=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)
Power_u_0.4_rhoXY_0.5_rhoX_0.7_n50=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)
Power_u_0.4_rhoXY_0.5_rhoX_0.8_n50=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)
Power_u_0.4_rhoXY_0.5_rhoX_0.9_n50=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)

Power_u_0.4_rhoXY_0.6_rhoX_0.7_n50=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)
Power_u_0.4_rhoXY_0.6_rhoX_0.8_n50=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)
Power_u_0.4_rhoXY_0.6_rhoX_0.9_n50=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)

Power_u_0.4_rhoXY_0.7_rhoX_0.8_n50=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)
Power_u_0.4_rhoXY_0.7_rhoX_0.9_n50=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=50,n1=50)

Power_u_0.4_n50=rbind(Power_u_0.4_rhoXY_0.5_rhoX_0.6_n50,
                      Power_u_0.4_rhoXY_0.5_rhoX_0.7_n50,
                      Power_u_0.4_rhoXY_0.5_rhoX_0.8_n50,
                      Power_u_0.4_rhoXY_0.5_rhoX_0.9_n50,
                      Power_u_0.4_rhoXY_0.6_rhoX_0.7_n50,
                      Power_u_0.4_rhoXY_0.6_rhoX_0.8_n50,
                      Power_u_0.4_rhoXY_0.6_rhoX_0.9_n50,
                      Power_u_0.4_rhoXY_0.7_rhoX_0.8_n50,
                      Power_u_0.4_rhoXY_0.7_rhoX_0.9_n50)

Power_u_0.4_n50
saveRDS(Power_u_0.4_n50, file="Power_u_0.4_n50.Rds")
#!No warning messages for u=0.4. But for u=0.3, 0.2, there is warnings messages for perfectly separation.

#larger sample size at n0=n1=75,100,125,150, run models 1,2,3
power_RB=function(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0,u=0.4,n0=50,n1=50){
  
  mu_C=rep(u,3)
  mu_T=rep(NA,3)
  mu_T[1:2]=rep(u,2)
  mu_T[3]=exp(beta+log(u/(1-u)))/(1+exp(beta+log(u/(1-u))))  
  
  pvalue1=rep(NA, nsim)
  pvalue2=rep(NA, nsim)
  pvalue3=rep(NA, nsim)
  
  set.seed(12345)
  
  for (i in 1: nsim){
    
    #Generate the data  
    myCop=normalCopula(param = c(rhoX, rhoXY, rhoXY),dim = 3, dispstr = "un")
    U=rCopula(n0+n1,myCop)
    
    #Generate the correlated binary outcomes from the Gaussian Coupla
    X1=c(U[1:n0,1]<=mu_C[1],U[(n0+1):(n0+n1),1]<=mu_T[1])
    X2=c(U[1:n0,2]<=mu_C[2],U[(n0+1):(n0+n1),2]<=mu_T[2])
    Y=c(U[1:n0,3]<=mu_C[3],U[(n0+1):(n0+n1),3]<=mu_T[3])
    
    X=X1+X2
    X_log=log((X+1/2)/(2-X+1/2))
    X_C=as.factor(X)
    
    Treat=c(rep('Control',n0),rep('Treatment',n1))
    
    fit1=glm(Y ~ Treat+X2,family=binomial)
    pvalue1[i]=coef(summary(fit1))[2,4]
    
    fit2=glm(Y ~ Treat+X_log,family=binomial())
    pvalue2[i]=coef(summary(fit2))[2,4]
    
    fit3=glm(Y ~ Treat+X_C,family=binomial)
    pvalue3[i]=coef(summary(fit3))[2,4]
  }
  
  power=c(sum(pvalue1<=0.05)/nsim,
          sum(pvalue2<=0.05)/nsim,
          sum(pvalue3<=0.05)/nsim)
  
  return(power)  
}

#n0=n1=75
TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n75=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0,u=0.4,n0=75,n1=75)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n75=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=75,n1=75)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n75=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=75,n1=75)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n75=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=75,n1=75)

TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n75=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=75,n1=75)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n75=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=75,n1=75)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n75=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=75,n1=75)

TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n75=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=75,n1=75)
TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n75=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=75,n1=75)

TypeI_u_0.4_n75=rbind(TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n75,
                      TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n75,
                      TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n75,
                      TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n75,
                      TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n75,
                      TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n75,
                      TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n75,
                      TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n75,
                      TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n75)

TypeI_u_0.4_n75
saveRDS(TypeI_u_0.4_n75, file="TypeI_u_0.4_n75.Rds")

Power_u_0.4_rhoXY_0.5_rhoX_0.6_n75=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)
Power_u_0.4_rhoXY_0.5_rhoX_0.7_n75=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)
Power_u_0.4_rhoXY_0.5_rhoX_0.8_n75=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)
Power_u_0.4_rhoXY_0.5_rhoX_0.9_n75=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)

Power_u_0.4_rhoXY_0.6_rhoX_0.7_n75=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)
Power_u_0.4_rhoXY_0.6_rhoX_0.8_n75=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)
Power_u_0.4_rhoXY_0.6_rhoX_0.9_n75=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)

Power_u_0.4_rhoXY_0.7_rhoX_0.8_n75=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)
Power_u_0.4_rhoXY_0.7_rhoX_0.9_n75=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=75,n1=75)

Power_u_0.4_n75=rbind(Power_u_0.4_rhoXY_0.5_rhoX_0.6_n75,
                      Power_u_0.4_rhoXY_0.5_rhoX_0.7_n75,
                      Power_u_0.4_rhoXY_0.5_rhoX_0.8_n75,
                      Power_u_0.4_rhoXY_0.5_rhoX_0.9_n75,
                      Power_u_0.4_rhoXY_0.6_rhoX_0.7_n75,
                      Power_u_0.4_rhoXY_0.6_rhoX_0.8_n75,
                      Power_u_0.4_rhoXY_0.6_rhoX_0.9_n75,
                      Power_u_0.4_rhoXY_0.7_rhoX_0.8_n75,
                      Power_u_0.4_rhoXY_0.7_rhoX_0.9_n75)

Power_u_0.4_n75
saveRDS(Power_u_0.4_n75, file="Power_u_0.4_n75.Rds")

#n0=n1=100
TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n100=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0,u=0.4,n0=100,n1=100)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n100=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=100,n1=100)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n100=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=100,n1=100)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n100=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=100,n1=100)

TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n100=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=100,n1=100)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n100=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=100,n1=100)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n100=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=100,n1=100)

TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n100=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=100,n1=100)
TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n100=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=100,n1=100)

TypeI_u_0.4_n100=rbind(TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n100,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n100,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n100,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n100,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n100,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n100,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n100,
                       TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n100,
                       TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n100)

TypeI_u_0.4_n100
saveRDS(TypeI_u_0.4_n100, file="TypeI_u_0.4_n100.Rds")

Power_u_0.4_rhoXY_0.5_rhoX_0.6_n100=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)
Power_u_0.4_rhoXY_0.5_rhoX_0.7_n100=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)
Power_u_0.4_rhoXY_0.5_rhoX_0.8_n100=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)
Power_u_0.4_rhoXY_0.5_rhoX_0.9_n100=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)

Power_u_0.4_rhoXY_0.6_rhoX_0.7_n100=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)
Power_u_0.4_rhoXY_0.6_rhoX_0.8_n100=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)
Power_u_0.4_rhoXY_0.6_rhoX_0.9_n100=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)

Power_u_0.4_rhoXY_0.7_rhoX_0.8_n100=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)
Power_u_0.4_rhoXY_0.7_rhoX_0.9_n100=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=100,n1=100)

Power_u_0.4_n100=rbind(Power_u_0.4_rhoXY_0.5_rhoX_0.6_n100,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.7_n100,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.8_n100,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.9_n100,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.7_n100,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.8_n100,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.9_n100,
                       Power_u_0.4_rhoXY_0.7_rhoX_0.8_n100,
                       Power_u_0.4_rhoXY_0.7_rhoX_0.9_n100)

Power_u_0.4_n100
saveRDS(Power_u_0.4_n100, file="Power_u_0.4_n100.Rds")
#n0=n1=125
TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n125=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0,u=0.4,n0=125,n1=125)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n125=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=125,n1=125)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n125=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=125,n1=125)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n125=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=125,n1=125)

TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n125=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=125,n1=125)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n125=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=125,n1=125)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n125=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=125,n1=125)

TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n125=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=125,n1=125)
TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n125=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=125,n1=125)

TypeI_u_0.4_n125=rbind(TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n125,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n125,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n125,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n125,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n125,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n125,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n125,
                       TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n125,
                       TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n125)

TypeI_u_0.4_n125
saveRDS(TypeI_u_0.4_n125, file="TypeI_u_0.4_n125.Rds")

Power_u_0.4_rhoXY_0.5_rhoX_0.6_n125=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)
Power_u_0.4_rhoXY_0.5_rhoX_0.7_n125=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)
Power_u_0.4_rhoXY_0.5_rhoX_0.8_n125=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)
Power_u_0.4_rhoXY_0.5_rhoX_0.9_n125=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)

Power_u_0.4_rhoXY_0.6_rhoX_0.7_n125=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)
Power_u_0.4_rhoXY_0.6_rhoX_0.8_n125=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)
Power_u_0.4_rhoXY_0.6_rhoX_0.9_n125=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)

Power_u_0.4_rhoXY_0.7_rhoX_0.8_n125=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)
Power_u_0.4_rhoXY_0.7_rhoX_0.9_n125=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=125,n1=125)

Power_u_0.4_n125=rbind(Power_u_0.4_rhoXY_0.5_rhoX_0.6_n125,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.7_n125,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.8_n125,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.9_n125,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.7_n125,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.8_n125,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.9_n125,
                       Power_u_0.4_rhoXY_0.7_rhoX_0.8_n125,
                       Power_u_0.4_rhoXY_0.7_rhoX_0.9_n125)

Power_u_0.4_n125
saveRDS(Power_u_0.4_n125, file="Power_u_0.4_n125.Rds")
#n0=n1=150
TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n150=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0,u=0.4,n0=150,n1=150)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n150=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=150,n1=150)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n150=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=150,n1=150)
TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n150=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=150,n1=150)

TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n150=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0,u=0.4,n0=150,n1=150)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n150=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=150,n1=150)
TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n150=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=150,n1=150)

TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n150=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0,u=0.4,n0=150,n1=150)
TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n150=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0,u=0.4,n0=150,n1=150)

TypeI_u_0.4_n150=rbind(TypeI_u_0.4_rhoXY_0.5_rhoX_0.6_n150,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.7_n150,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.8_n150,
                       TypeI_u_0.4_rhoXY_0.5_rhoX_0.9_n150,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.7_n150,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.8_n150,
                       TypeI_u_0.4_rhoXY_0.6_rhoX_0.9_n150,
                       TypeI_u_0.4_rhoXY_0.7_rhoX_0.8_n150,
                       TypeI_u_0.4_rhoXY_0.7_rhoX_0.9_n150)

TypeI_u_0.4_n150
saveRDS(TypeI_u_0.4_n150, file="TypeI_u_0.4_n150.Rds")

Power_u_0.4_rhoXY_0.5_rhoX_0.6_n150=power_RB(rhoXY=0.5,rhoX=0.6,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)
Power_u_0.4_rhoXY_0.5_rhoX_0.7_n150=power_RB(rhoXY=0.5,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)
Power_u_0.4_rhoXY_0.5_rhoX_0.8_n150=power_RB(rhoXY=0.5,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)
Power_u_0.4_rhoXY_0.5_rhoX_0.9_n150=power_RB(rhoXY=0.5,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)

Power_u_0.4_rhoXY_0.6_rhoX_0.7_n150=power_RB(rhoXY=0.6,rhoX=0.7,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)
Power_u_0.4_rhoXY_0.6_rhoX_0.8_n150=power_RB(rhoXY=0.6,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)
Power_u_0.4_rhoXY_0.6_rhoX_0.9_n150=power_RB(rhoXY=0.6,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)

Power_u_0.4_rhoXY_0.7_rhoX_0.8_n150=power_RB(rhoXY=0.7,rhoX=0.8,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)
Power_u_0.4_rhoXY_0.7_rhoX_0.9_n150=power_RB(rhoXY=0.7,rhoX=0.9,nsim=20000,beta=0.8,u=0.4,n0=150,n1=150)

Power_u_0.4_n150=rbind(Power_u_0.4_rhoXY_0.5_rhoX_0.6_n150,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.7_n150,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.8_n150,
                       Power_u_0.4_rhoXY_0.5_rhoX_0.9_n150,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.7_n150,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.8_n150,
                       Power_u_0.4_rhoXY_0.6_rhoX_0.9_n150,
                       Power_u_0.4_rhoXY_0.7_rhoX_0.8_n150,
                       Power_u_0.4_rhoXY_0.7_rhoX_0.9_n150)

Power_u_0.4_n150
saveRDS(Power_u_0.4_n150, file="Power_u_0.4_n150.Rds")