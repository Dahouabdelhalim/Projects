######################################################################################
### R Code Simulation Study of Multiple Continuous Outcomes
### 2023-1-30
### R version 3.6.0
######################################################################################

#rm(list=ls())
#setwd("~/shiyangm")
#.libPaths()

library(mvtnorm)

##n0=n1=50,100,150
Sopt_ANCOVA=function(rhoXY=0.5,rhoX=0.6,n0=50,n1=50,nsim=20000){
  
  set.seed(12345)
  M=10
  rhoY=rhoX
  #S=1, 2, 3, ..., M-1
  TypeI=rep(NA, M-1)
  Power=rep(NA, M-1)
  
  for (s in 1:(M-1)){
    
    S=(1:(M-1))[s]
    
    A=diag(1,nrow=S)
    A[row(A)!=col(A)]=rhoX
    
    B=diag(1,nrow=M-S)
    B[row(B)!=col(B)]=rhoY
    
    C=matrix(rhoXY,nrow=S,ncol=M-S)
    R=rbind(cbind(A,C), cbind(t(C),B))
  
    ############  Type I error rate ##########################
    delta=0
    mu_C=rep(0, M) 
    mu_T=rep(NA,M)
    mu_T[1:S]=rep(0,S)
    mu_T[(S+1):M]=delta

    pvalue_TypeI=rep(NA, nsim)
    for (i in 1: nsim){
      
      ######Generate the data ########
      #placebo
      P=rmvnorm(n0, mu_C, sigma=R)
      #treatment
      T=rmvnorm(n1, mu_T, sigma=R)
      
      Treat=c(rep('Control',n0),rep('Treatment',n1))
      All=rbind(P,T)

      #pre-treatment measurements
      X=All[,1:S]
      Y=All[,(S+1):10]

      if (S==1){
        pre_mean=X
        post_mean=rowMeans(Y)
      }else if (S==9){
        pre_mean=rowMeans(X)
        post_mean=Y
      } else{  
        pre_mean=rowMeans(X)
        post_mean=rowMeans(Y)
      }
      
      #Fitting Centered ANCOVA Model
      cpre_mean=pre_mean-mean(pre_mean)
      model <- lm(post_mean ~Treat+cpre_mean)
      pvalue_TypeI[i]=summary(model)$coefficients[2,4]
    }  
    TypeI[s]=sum(pvalue_TypeI<=0.05)/nsim
    
    ############  Power ##########################
    delta=0.25
    mu_C=rep(0, M) 
    mu_T=rep(NA,M)
    mu_T[1:S]=rep(0,S)
    mu_T[(S+1):M]=delta
  
    pvalue_Power=rep(NA, nsim)
    for (i in 1: nsim){
      
      ######Generate the data ########
      #placebo
      P=rmvnorm(n0, mu_C, sigma=R)
      #treatment
      T=rmvnorm(n1, mu_T, sigma=R)
      
      Treat=c(rep('Control',n0),rep('Treatment',n1))
      All=rbind(P,T)
      
      #pre-treatment measurements
      X=All[,1:S]
      Y=All[,(S+1):10]
      
      if (S==1){
        pre_mean=X
        post_mean=rowMeans(Y)
      }else if (S==9){
        pre_mean=rowMeans(X)
        post_mean=Y
      } else{  
        pre_mean=rowMeans(X)
        post_mean=rowMeans(Y)
      }
      
      #Fitting Centered ANCOVA Model
      cpre_mean=pre_mean-mean(pre_mean)
      model <- lm(post_mean ~Treat+cpre_mean)
      pvalue_Power[i]=summary(model)$coefficients[2,4]
    }  
    Power[s]=sum(pvalue_Power<=0.05)/nsim
  } 
  
  return(rbind(TypeI,Power))
  
}

######################## n0=n1=50 ####################################
### 1 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.5_rhoX_0.6_n50=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.6,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.5_rhoX_0.6_n50
ANCOVA_rhoXY_0.5_rhoX_0.7_n50=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.7,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.5_rhoX_0.7_n50
ANCOVA_rhoXY_0.5_rhoX_0.8_n50=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.8,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.5_rhoX_0.8_n50

ANCOVA_n50_1=rbind(ANCOVA_rhoXY_0.5_rhoX_0.6_n50,
                   ANCOVA_rhoXY_0.5_rhoX_0.7_n50,
                   ANCOVA_rhoXY_0.5_rhoX_0.8_n50)
saveRDS(ANCOVA_n50_1, file="ANCOVA_n50_1.Rds")

end_time <- Sys.time()
end_time - start_time #43.95138 mins

### 2 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.5_rhoX_0.9_n50=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.9,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.5_rhoX_0.9_n50
ANCOVA_rhoXY_0.6_rhoX_0.7_n50=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.7,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.6_rhoX_0.7_n50
ANCOVA_rhoXY_0.6_rhoX_0.8_n50=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.8,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.6_rhoX_0.8_n50

ANCOVA_n50_2=rbind(ANCOVA_rhoXY_0.5_rhoX_0.9_n50,
                   ANCOVA_rhoXY_0.6_rhoX_0.7_n50,
                   ANCOVA_rhoXY_0.6_rhoX_0.8_n50)
saveRDS(ANCOVA_n50_2, file="ANCOVA_n50_2.Rds")

end_time <- Sys.time()
end_time - start_time  #46.48444 mins

### 3 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.6_rhoX_0.9_n50=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.9,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.6_rhoX_0.9_n50
ANCOVA_rhoXY_0.7_rhoX_0.8_n50=Sopt_ANCOVA(rhoXY=0.7,rhoX=0.8,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.7_rhoX_0.8_n50
ANCOVA_rhoXY_0.7_rhoX_0.9_n50=Sopt_ANCOVA(rhoXY=0.7,rhoX=0.9,n0=50,n1=50,nsim=20000)
ANCOVA_rhoXY_0.7_rhoX_0.9_n50

ANCOVA_n50_3=rbind(ANCOVA_rhoXY_0.6_rhoX_0.9_n50,
                   ANCOVA_rhoXY_0.7_rhoX_0.8_n50,
                   ANCOVA_rhoXY_0.7_rhoX_0.9_n50)
saveRDS(ANCOVA_n50_3, file="ANCOVA_n50_3.Rds")
end_time <- Sys.time()
end_time - start_time  #45.86458 mins

######################## n0=n1=100 ####################################
### 4 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.5_rhoX_0.6_n100=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.6,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.5_rhoX_0.6_n100
ANCOVA_rhoXY_0.5_rhoX_0.7_n100=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.7,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.5_rhoX_0.7_n100
ANCOVA_rhoXY_0.5_rhoX_0.8_n100=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.8,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.5_rhoX_0.8_n100

ANCOVA_n100_1=rbind(ANCOVA_rhoXY_0.5_rhoX_0.6_n100,
                        ANCOVA_rhoXY_0.5_rhoX_0.7_n100,
                        ANCOVA_rhoXY_0.5_rhoX_0.8_n100)
saveRDS(ANCOVA_n100_1, file="ANCOVA_n100_1.Rds")

end_time <- Sys.time() 
end_time - start_time #48.75023 mins

### 5 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.5_rhoX_0.9_n100=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.9,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.5_rhoX_0.9_n100
ANCOVA_rhoXY_0.6_rhoX_0.7_n100=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.7,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.6_rhoX_0.7_n100
ANCOVA_rhoXY_0.6_rhoX_0.8_n100=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.8,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.6_rhoX_0.8_n100

ANCOVA_n100_2=rbind(ANCOVA_rhoXY_0.5_rhoX_0.9_n100,
                        ANCOVA_rhoXY_0.6_rhoX_0.7_n100,
                        ANCOVA_rhoXY_0.6_rhoX_0.8_n100)
saveRDS(ANCOVA_n100_2, file="ANCOVA_n100_2.Rds")

end_time <- Sys.time() 
end_time - start_time  #50.47581 mins

### 6 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.6_rhoX_0.9_n100=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.9,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.6_rhoX_0.9_n100
ANCOVA_rhoXY_0.7_rhoX_0.8_n100=Sopt_ANCOVA(rhoXY=0.7,rhoX=0.8,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.7_rhoX_0.8_n100
ANCOVA_rhoXY_0.7_rhoX_0.9_n100=Sopt_ANCOVA(rhoXY=0.7,rhoX=0.9,n0=100,n1=100,nsim=20000) 
ANCOVA_rhoXY_0.7_rhoX_0.9_n100

ANCOVA_n100_3=rbind(ANCOVA_rhoXY_0.6_rhoX_0.9_n100,
                        ANCOVA_rhoXY_0.7_rhoX_0.8_n100,
                        ANCOVA_rhoXY_0.7_rhoX_0.9_n100)
saveRDS(ANCOVA_n100_3, file="ANCOVA_n100_3.Rds")
end_time <- Sys.time() 
end_time - start_time  #50.70022 mins

######################## n0=n1=150 ####################################
### 7 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.5_rhoX_0.6_n150=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.6,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.5_rhoX_0.6_n150
ANCOVA_rhoXY_0.5_rhoX_0.7_n150=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.7,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.5_rhoX_0.7_n150
ANCOVA_rhoXY_0.5_rhoX_0.8_n150=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.8,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.5_rhoX_0.8_n150

ANCOVA_n150_1=rbind(ANCOVA_rhoXY_0.5_rhoX_0.6_n150,
                        ANCOVA_rhoXY_0.5_rhoX_0.7_n150,
                        ANCOVA_rhoXY_0.5_rhoX_0.8_n150)
saveRDS(ANCOVA_n150_1, file="ANCOVA_n150_1.Rds")

end_time <- Sys.time()
end_time - start_time #50.64473 mins

### 8 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.5_rhoX_0.9_n150=Sopt_ANCOVA(rhoXY=0.5,rhoX=0.9,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.5_rhoX_0.9_n150
ANCOVA_rhoXY_0.6_rhoX_0.7_n150=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.7,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.6_rhoX_0.7_n150
ANCOVA_rhoXY_0.6_rhoX_0.8_n150=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.8,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.6_rhoX_0.8_n150

ANCOVA_n150_2=rbind(ANCOVA_rhoXY_0.5_rhoX_0.9_n150,
                        ANCOVA_rhoXY_0.6_rhoX_0.7_n150,
                        ANCOVA_rhoXY_0.6_rhoX_0.8_n150)
saveRDS(ANCOVA_n150_2, file="ANCOVA_n150_2.Rds")

end_time <- Sys.time()
end_time - start_time  #53.58718 mins

### 9 ###
start_time <- Sys.time()
ANCOVA_rhoXY_0.6_rhoX_0.9_n150=Sopt_ANCOVA(rhoXY=0.6,rhoX=0.9,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.6_rhoX_0.9_n150
ANCOVA_rhoXY_0.7_rhoX_0.8_n150=Sopt_ANCOVA(rhoXY=0.7,rhoX=0.8,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.7_rhoX_0.8_n150
ANCOVA_rhoXY_0.7_rhoX_0.9_n150=Sopt_ANCOVA(rhoXY=0.7,rhoX=0.9,n0=150,n1=150,nsim=20000)
ANCOVA_rhoXY_0.7_rhoX_0.9_n150

ANCOVA_n150_3=rbind(ANCOVA_rhoXY_0.6_rhoX_0.9_n150,
                        ANCOVA_rhoXY_0.7_rhoX_0.8_n150,
                        ANCOVA_rhoXY_0.7_rhoX_0.9_n150)
saveRDS(ANCOVA_n150_3, file="ANCOVA_n150_3.Rds")
end_time <- Sys.time()
end_time - start_time  #52.57698 mins