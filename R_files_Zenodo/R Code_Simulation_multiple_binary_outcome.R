######################################################################################
### R Code Simulation Study of Multiple Binary Outcomes
### 2022-11-27
### R version 4.1.3
### dichotomized probability (sp) u=0.5
######################################################################################

#rm(list=ls())
setwd("~/shiyangm")
#.libPaths()

library(mvtnorm)
library(geepack)
library(copula)

#small samples n0=n1=50
Sopt_GEE_small=function(rhoXY=0.5,rhoX=0.6,n0=50,n1=50,nsim=20000,u=0.4){
  #rhoXY=0.5;rhoX=0.6;n0=50;n1=50;nsim=100;u=0.4
  
  #Small sample, only consider GEE model 2
  set.seed(12345)
  #u is the success probability (P(Y=1)) for the control group
  M=10
  rhoY=rhoX
  #S=1, 2, 3, ..., M-1
  TypeI_orginal=rep(NA, M-1)
  TypeI=rep(NA, M-1)
  Power_orginal=rep(NA, M-1)
  Power=rep(NA, M-1)
  
  for (s in 1:(M-1)){
    S=(1:(M-1))[s]
    
    A=diag(1,nrow=S)
    A[row(A)!=col(A)]=rhoX
    
    B=diag(1,nrow=M-S)
    B[row(B)!=col(B)]=rhoY
    
    C=matrix(rhoXY,nrow=S,ncol=M-S)
    R=rbind(cbind(A,C), cbind(t(C),B))
    
    param=c(R[1,2:10],
            R[2,3:10],
            R[3,4:10],
            R[4,5:10],
            R[5,6:10],
            R[6,7:10],
            R[7,8:10],
            R[8,9:10],
            R[9,10:10])
    
    #Type I error rate
    beta=0
    mu_C=rep(u, M) 
    mu_T=rep(NA,M)
    mu_T[1:S]=rep(u,S)
    mu_T[(S+1):M]=exp(beta+log(u/(1-u)))/(1+exp(beta+log(u/(1-u))))  
    
    ##################### For Type I error rate ##########################
    zvalue_TypeI=rep(NA, nsim)
    
    for (i in 1: nsim){
      
      ######Generate the data ########
      myCop=normalCopula(param, dim = M, dispstr = "un")
      U=rCopula(n0+n1,myCop)
      
      #Generate the correlated binary outcomes from the Gaussian Coupla
      #pre-treatment measurements
      if (S==1){
        X=c(U[1:n0,1:S]<=mu_C[1:S],U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=rbind(U[1:n0,(S+1):M]<=mu_C[(S+1):M],
                U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        YY=c(t(Y))
        X_sum=as.numeric(rep(X,rep(M-S, (n0+n1)))) 
        X_sum_C=as.factor(X_sum)
        X_log=log((as.numeric(X_sum)+1/2)/(S-as.numeric(X_sum)+1/2))
        
      }else if (S==9){
        X=rbind(U[1:n0,1:S]<=mu_C[1:S],
                U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=c(U[1:n0,(S+1):M]<=mu_C[(S+1):M],U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        X_sum=rowSums(X) 
        X_sum_C=as.factor(X_sum)
        X_log=log((X_sum+1/2)/(S-X_sum+1/2))
        
      } else{  
        X=rbind(U[1:n0,1:S]<=mu_C[1:S],
                U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=rbind(U[1:n0,(S+1):M]<=mu_C[(S+1):M],
                U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        YY=c(t(Y))
        X_sum=rep(rowSums(X),rep(M-S, (n0+n1))) 
        X_sum_C=as.factor(X_sum)
        X_log=log((X_sum+1/2)/(S-X_sum+1/2))
      }
      
      Treat=c(rep('Control',n0*(M-S)),rep('Treatment',n1*(M-S)))
      subject=rep(seq(1,(n0+n1)),rep((M-S),(n0+n1)))
      
      if (S==9){   
        fit=glm(Y~Treat+X_log,family=binomial)
        zvalue_TypeI[i]=coef(summary(fit))[2,3]
        
      }else{
        fit_GEE=geeglm(YY~ Treat+X_log, family=binomial,id =subject, corstr = 'exchangeable')
        zvalue_TypeI[i]=coef(summary(fit_GEE))[2,1]/coef(summary(fit_GEE))[2,2]
      } 
    }
    
    #2.5% and 97.5% quantiles
    lower_cutoff=quantile(zvalue_TypeI, probs =0.025)
    upper_cutoff=quantile(zvalue_TypeI, probs =0.975)
    
    TypeI_orginal[s]=(sum(abs(zvalue_TypeI)>qnorm(0.975)))/nsim
    #type I error calibration
    TypeI[s]=(sum(zvalue_TypeI<lower_cutoff)+sum(zvalue_TypeI>upper_cutoff))/nsim
    
    ###############################################################################
    #Power
    beta=0.5
    mu_C=rep(u, M) 
    mu_T=rep(NA,M)
    mu_T[1:S]=rep(u,S)
    mu_T[(S+1):M]=exp(beta+log(u/(1-u)))/(1+exp(beta+log(u/(1-u))))  
    
    ##################### For Power ##########################   
    zvalue_Power=rep(NA, nsim)
    
    for (i in 1: nsim){
      
      ######Generate the data ########
      myCop=normalCopula(param, dim = M, dispstr = "un")
      U=rCopula(n0+n1,myCop)
      
      #Generate the correlated binary outcomes from the Gaussian Coupla
      #pre-treatment measurements
      if (S==1){
        X=c(U[1:n0,1:S]<=mu_C[1:S],U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=rbind(U[1:n0,(S+1):M]<=mu_C[(S+1):M],
                U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        YY=c(t(Y))
        X_sum=as.numeric(rep(X,rep(M-S, (n0+n1)))) 
        X_sum_C=as.factor(X_sum)
        X_log=log((as.numeric(X_sum)+1/2)/(S-as.numeric(X_sum)+1/2))
        
      }else if (S==9){
        X=rbind(U[1:n0,1:S]<=mu_C[1:S],
                U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=c(U[1:n0,(S+1):M]<=mu_C[(S+1):M],U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        X_sum=rowSums(X) 
        X_sum_C=as.factor(X_sum)
        X_log=log((X_sum+1/2)/(S-X_sum+1/2))
        
      } else{  
        X=rbind(U[1:n0,1:S]<=mu_C[1:S],
                U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=rbind(U[1:n0,(S+1):M]<=mu_C[(S+1):M],
                U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        YY=c(t(Y))
        X_sum=rep(rowSums(X),rep(M-S, (n0+n1))) 
        X_sum_C=as.factor(X_sum)
        X_log=log((X_sum+1/2)/(S-X_sum+1/2))
      }
      
      Treat=c(rep('Control',n0*(M-S)),rep('Treatment',n1*(M-S)))
      subject=rep(seq(1,(n0+n1)),rep((M-S),(n0+n1)))
      
      if (S==9){   
        fit=glm(Y~Treat+X_log,family=binomial)
        zvalue_Power[i]=coef(summary(fit))[2,3]
        
      }else{
        fit_GEE=geeglm(YY~ Treat+X_log, family=binomial,id =subject, corstr = 'exchangeable')
        zvalue_Power[i]=coef(summary(fit_GEE))[2,1]/coef(summary(fit_GEE))[2,2]
      } 
    }
    
    Power_orginal[s]=(sum(abs(zvalue_Power)>qnorm(0.975)))/nsim
    Power[s]=(sum(zvalue_Power<lower_cutoff)+sum(zvalue_Power>upper_cutoff))/nsim
    
  } 
  
  return(rbind(TypeI_orginal,TypeI,Power_orginal,Power))
  
}

######################## n0=n1=50 ####################################
### 1 ###
start_time <- Sys.time()
GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.5,rhoX=0.6,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n50
GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.5,rhoX=0.7,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n50
GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.5,rhoX=0.8,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n50
GEE_sp_0.4_n50_1=rbind(GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n50,
                       GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n50,
                       GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n50)
saveRDS(GEE_sp_0.4_n50_1, file="GEE_sp_0.4_n50_1.Rds")
end_time <- Sys.time() 
end_time - start_time #4.991398 hours

### 2 ###
start_time <- Sys.time()
GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.5,rhoX=0.9,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n50
GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.6,rhoX=0.7,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n50
GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.6,rhoX=0.8,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n50
GEE_sp_0.4_n50_2=rbind(GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n50,
                       GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n50,
                       GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n50)
saveRDS(GEE_sp_0.4_n50_2, file="GEE_sp_0.4_n50_2.Rds")
end_time <- Sys.time() 
end_time - start_time #5.305044 hours

### 3 ###
start_time <- Sys.time()
GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.6,rhoX=0.9,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n50
GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.7,rhoX=0.8,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n50
GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n50=Sopt_GEE_small(rhoXY=0.7,rhoX=0.9,n0=50,n1=50,nsim=20000,u=0.4) 
GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n50
GEE_sp_0.4_n50_3=rbind(GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n50,
                       GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n50,
                       GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n50)
saveRDS(GEE_sp_0.4_n50_3, file="GEE_sp_0.4_n50_3.Rds")
end_time <- Sys.time() 
end_time - start_time #

##n0=n1=100,150
Sopt_GEE=function(rhoXY=0.5,rhoX=0.6,n0=100,n1=100,nsim=20000,u=0.4){
  
  set.seed(12345)
  #sp is the success probability (P(Y=1)) for the control group
  M=10
  rhoY=rhoX
  #S=1, 2, 3, ..., M-1
  TypeI1_orginal=rep(NA, M-1)
  TypeI2_orginal=rep(NA, M-1)
  
  TypeI1=rep(NA, M-1)
  TypeI2=rep(NA, M-1)
  
  Power1_orginal=rep(NA, M-1)
  Power2_orginal=rep(NA, M-1)
  
  Power1=rep(NA, M-1)
  Power2=rep(NA, M-1)
  
  for (s in 1:(M-1)){
    
    S=(1:(M-1))[s]
    
    A=diag(1,nrow=S)
    A[row(A)!=col(A)]=rhoX
    
    B=diag(1,nrow=M-S)
    B[row(B)!=col(B)]=rhoY
    
    C=matrix(rhoXY,nrow=S,ncol=M-S)
    R=rbind(cbind(A,C), cbind(t(C),B))
    
    param=c(R[1,2:10],
            R[2,3:10],
            R[3,4:10],
            R[4,5:10],
            R[5,6:10],
            R[6,7:10],
            R[7,8:10],
            R[8,9:10],
            R[9,10:10])
    
    #Type I error rate
    beta=0
    mu_C=rep(u, M) 
    mu_T=rep(NA,M)
    mu_T[1:S]=rep(u,S)
    mu_T[(S+1):M]=exp(beta+log(u/(1-u)))/(1+exp(beta+log(u/(1-u))))  
    
    #exp(beta)*u/(1-u+exp(beta)*u)  
    ##################### For Type I error rate ##########################
    zvalue1_TypeI=rep(NA, nsim)
    zvalue2_TypeI=rep(NA, nsim)
    
    for (i in 1: nsim){
      
      ######Generate the data ########
      myCop=normalCopula(param, dim = M, dispstr = "un")
      U=rCopula(n0+n1,myCop)
      
      #Generate the correlated binary outcomes from the Gaussian Coupla
      #pre-treatment measurements
      if (S==1){
        X=c(U[1:n0,1:S]<=mu_C[1:S],U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=rbind(U[1:n0,(S+1):M]<=mu_C[(S+1):M],
                U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        YY=c(t(Y))
        X_sum=as.numeric(rep(X,rep(M-S, (n0+n1)))) 
        X_log=log((as.numeric(X_sum)+1/2)/(S-as.numeric(X_sum)+1/2))
        
      }else if (S==9){
        X=rbind(U[1:n0,1:S]<=mu_C[1:S],
                U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=c(U[1:n0,(S+1):M]<=mu_C[(S+1):M],U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        X_sum=rowSums(X) 
        X_log=log((X_sum+1/2)/(S-X_sum+1/2))
        
      } else{  
        X=rbind(U[1:n0,1:S]<=mu_C[1:S],
                U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=rbind(U[1:n0,(S+1):M]<=mu_C[(S+1):M],
                U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        YY=c(t(Y))
        X_sum=rep(rowSums(X),rep(M-S, (n0+n1))) 
        X_log=log((X_sum+1/2)/(S-X_sum+1/2))
      }
      
      Treat=c(rep('Control',n0*(M-S)),rep('Treatment',n1*(M-S)))
      subject=rep(seq(1,(n0+n1)),rep((M-S),(n0+n1)))
      
      if (S==9){   
        fit1=glm(Y~Treat+X_sum,family=binomial)
        zvalue1_TypeI[i]=coef(summary(fit1))[2,3]
        fit2=glm(Y~Treat+X_log,family=binomial)
        zvalue2_TypeI[i]=coef(summary(fit2))[2,3]
        
      }else{
        fit1_GEE=geeglm(YY~ Treat+X_sum, family=binomial,id =subject, corstr = 'exchangeable')
        fit2_GEE=geeglm(YY~ Treat+X_log, family=binomial,id =subject, corstr = 'exchangeable')
        
        zvalue1_TypeI[i]=coef(summary(fit1_GEE))[2,1]/coef(summary(fit1_GEE))[2,2]
        zvalue2_TypeI[i]=coef(summary(fit2_GEE))[2,1]/coef(summary(fit2_GEE))[2,2]
      } 
    }
    
    lower_cutoff1=quantile(zvalue1_TypeI, probs =0.025)
    upper_cutoff1=quantile(zvalue1_TypeI, probs =0.975)
    lower_cutoff2=quantile(zvalue2_TypeI, probs =0.025)
    upper_cutoff2=quantile(zvalue2_TypeI, probs =0.975)
    
    TypeI1_orginal[s]=(sum(abs(zvalue1_TypeI)>qnorm(0.975)))/nsim
    TypeI2_orginal[s]=(sum(abs(zvalue2_TypeI)>qnorm(0.975)))/nsim
    
    TypeI1[s]=(sum(zvalue1_TypeI<lower_cutoff1)+sum(zvalue1_TypeI>upper_cutoff1))/nsim
    TypeI2[s]=(sum(zvalue2_TypeI<lower_cutoff2)+sum(zvalue2_TypeI>upper_cutoff2))/nsim
    
    ###############################################################################
    #Power
    beta=0.5
    mu_C=rep(u, M) 
    mu_T=rep(NA,M)
    mu_T[1:S]=rep(u,S)
    mu_T[(S+1):M]=exp(beta+log(u/(1-u)))/(1+exp(beta+log(u/(1-u))))  
    
    ##################### For Power ##########################   
    zvalue1_Power=rep(NA, nsim)
    zvalue2_Power=rep(NA, nsim)
    
    for (i in 1: nsim){
      
      ######Generate the data ########
      myCop=normalCopula(param, dim = M, dispstr = "un")
      U=rCopula(n0+n1,myCop)
      
      #Generate the correlated binary outcomes from the Gaussian Coupla
      #pre-treatment measurements
      if (S==1){
        X=c(U[1:n0,1:S]<=mu_C[1:S],U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=rbind(U[1:n0,(S+1):M]<=mu_C[(S+1):M],
                U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        YY=c(t(Y))
        X_sum=as.numeric(rep(X,rep(M-S, (n0+n1)))) 
        X_log=log((as.numeric(X_sum)+1/2)/(S-as.numeric(X_sum)+1/2))
        
      }else if (S==9){
        X=rbind(U[1:n0,1:S]<=mu_C[1:S],
                U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=c(U[1:n0,(S+1):M]<=mu_C[(S+1):M],U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        X_sum=rowSums(X) 
        X_log=log((X_sum+1/2)/(S-X_sum+1/2))
        
      } else{  
        X=rbind(U[1:n0,1:S]<=mu_C[1:S],
                U[(n0+1):(n0+n1),1:S]<=mu_T[1:S])
        Y=rbind(U[1:n0,(S+1):M]<=mu_C[(S+1):M],
                U[(n0+1):(n0+n1),(S+1):M]<=mu_T[(S+1):M])
        YY=c(t(Y))
        X_sum=rep(rowSums(X),rep(M-S, (n0+n1))) 
        X_log=log((X_sum+1/2)/(S-X_sum+1/2))
      }
      
      Treat=c(rep('Control',n0*(M-S)),rep('Treatment',n1*(M-S)))
      subject=rep(seq(1,(n0+n1)),rep((M-S),(n0+n1)))
      
      if (S==9){   
        fit1=glm(Y~Treat+X_sum,family=binomial)
        zvalue1_Power[i]=coef(summary(fit1))[2,3]
        fit2=glm(Y~Treat+X_log,family=binomial)
        zvalue2_Power[i]=coef(summary(fit2))[2,3]
        
      }else{
        fit1_GEE=geeglm(YY~ Treat+X_sum, family=binomial,id =subject, corstr = 'exchangeable')
        fit2_GEE=geeglm(YY~ Treat+X_log, family=binomial,id =subject, corstr = 'exchangeable')
        
        zvalue1_Power[i]=coef(summary(fit1_GEE))[2,1]/coef(summary(fit1_GEE))[2,2]
        zvalue2_Power[i]=coef(summary(fit2_GEE))[2,1]/coef(summary(fit2_GEE))[2,2]
      } 
    }
    
    Power1_orginal[s]=(sum(abs(zvalue1_Power)>qnorm(0.975)))/nsim
    Power2_orginal[s]=(sum(abs(zvalue2_Power)>qnorm(0.975)))/nsim
    
    Power1[s]=(sum(zvalue1_Power<lower_cutoff1)+sum(zvalue1_Power>upper_cutoff1))/nsim
    Power2[s]=(sum(zvalue2_Power<lower_cutoff2)+sum(zvalue2_Power>upper_cutoff2))/nsim
    
  } 
  
  return(rbind(TypeI1_orginal,TypeI2_orginal,TypeI1,TypeI2,
               Power1_orginal,Power2_orginal,Power1,Power2))
  
}


######################## n0=n1=100 ####################################
### 4 ###
start_time <- Sys.time()
GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n100=Sopt_GEE(rhoXY=0.5,rhoX=0.6,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n100
GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n100=Sopt_GEE(rhoXY=0.5,rhoX=0.7,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n100
GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n100=Sopt_GEE(rhoXY=0.5,rhoX=0.8,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n100

GEE_sp_0.4_n100_1=rbind(GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n100,
                        GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n100,
                        GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n100)
saveRDS(GEE_sp_0.4_n100_1, file="GEE_sp_0.4_n100_1.Rds")

end_time <- Sys.time() 
end_time - start_time #14.4965 hours  

### 5 ###
start_time <- Sys.time()
GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n100=Sopt_GEE(rhoXY=0.5,rhoX=0.9,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n100
GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n100=Sopt_GEE(rhoXY=0.6,rhoX=0.7,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n100
GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n100=Sopt_GEE(rhoXY=0.6,rhoX=0.8,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n100

GEE_sp_0.4_n100_2=rbind(GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n100,
                        GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n100,
                        GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n100)
saveRDS(GEE_sp_0.4_n100_2, file="GEE_sp_0.4_n100_2.Rds")

end_time <- Sys.time() 
end_time - start_time  #14.3823 hours

### 6 ###
start_time <- Sys.time()
GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n100=Sopt_GEE(rhoXY=0.6,rhoX=0.9,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n100
GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n100=Sopt_GEE(rhoXY=0.7,rhoX=0.8,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n100
GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n100=Sopt_GEE(rhoXY=0.7,rhoX=0.9,n0=100,n1=100,nsim=20000,u=0.4) 
GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n100

GEE_sp_0.4_n100_3=rbind(GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n100,
                        GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n100,
                        GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n100)
saveRDS(GEE_sp_0.4_n100_3, file="GEE_sp_0.4_n100_3.Rds")
end_time <- Sys.time() 
end_time - start_time  #14.63787 hours

######################## n0=n1=150 ####################################
### 7 ###
start_time <- Sys.time()
GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n150=Sopt_GEE(rhoXY=0.5,rhoX=0.6,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n150
GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n150=Sopt_GEE(rhoXY=0.5,rhoX=0.7,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n150
GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n150=Sopt_GEE(rhoXY=0.5,rhoX=0.8,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n150

GEE_sp_0.4_n150_1=rbind(GEE_rhoXY_0.5_rhoX_0.6_sp_0.4_n150,
                        GEE_rhoXY_0.5_rhoX_0.7_sp_0.4_n150,
                        GEE_rhoXY_0.5_rhoX_0.8_sp_0.4_n150)
saveRDS(GEE_sp_0.4_n150_1, file="GEE_sp_0.4_n150_1.Rds")

end_time <- Sys.time()
end_time - start_time #19.86245 hours

### 8 ###
start_time <- Sys.time()
GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n150=Sopt_GEE(rhoXY=0.5,rhoX=0.9,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n150
GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n150=Sopt_GEE(rhoXY=0.6,rhoX=0.7,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n150
GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n150=Sopt_GEE(rhoXY=0.6,rhoX=0.8,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n150

GEE_sp_0.4_n150_2=rbind(GEE_rhoXY_0.5_rhoX_0.9_sp_0.4_n150,
                        GEE_rhoXY_0.6_rhoX_0.7_sp_0.4_n150,
                        GEE_rhoXY_0.6_rhoX_0.8_sp_0.4_n150)
saveRDS(GEE_sp_0.4_n150_2, file="GEE_sp_0.4_n150_2.Rds")

end_time <- Sys.time()
end_time - start_time  #20.14013 hours

### 9 ###
start_time <- Sys.time()
GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n150=Sopt_GEE(rhoXY=0.6,rhoX=0.9,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n150
GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n150=Sopt_GEE(rhoXY=0.7,rhoX=0.8,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n150
GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n150=Sopt_GEE(rhoXY=0.7,rhoX=0.9,n0=150,n1=150,nsim=20000,u=0.4)
GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n150

GEE_sp_0.4_n150_3=rbind(GEE_rhoXY_0.6_rhoX_0.9_sp_0.4_n150,
                        GEE_rhoXY_0.7_rhoX_0.8_sp_0.4_n150,
                        GEE_rhoXY_0.7_rhoX_0.9_sp_0.4_n150)
saveRDS(GEE_sp_0.4_n150_3, file="GEE_sp_0.4_n150_3.Rds")
end_time <- Sys.time()
end_time - start_time  #19.28277 hours