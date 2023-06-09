######################################################################################
### R Code Simulation Study of Single Continuous Outcomes
### 2023-1-30
### R version 3.6.0
######################################################################################

# qrsh
# qlogin -l mem=1G
# module load R/3.6.0
# R

#rm(list=ls())
#setwd("~/shiyangm")
#install.packages('mvtnorm')
library(mvtnorm)

#different sample sizes at n0=n1=50,75,100,125,150
power_RB_continuous=function(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0,n0=50,n1=50){
  
  Sigma <- matrix(c(1, rhoX, rhoXY,
                    rhoX, 1, rhoXY,
                    rhoXY, rhoXY,1),
                  ncol=3)
  
  pvalue1=rep(NA, nsim)
  pvalue2=rep(NA, nsim)
  
  set.seed(12345)
  
  for (i in 1: nsim){
    
    #Generate the data  
    #placebo
    P=rmvnorm(n0, mean=c(0,0,0), sigma=Sigma)
    #treatment
    T=rmvnorm(n1, mean=c(0,0,delta), sigma=Sigma)

    X1=c(P[,1],T[,1])
    X2=c(P[,2],T[,2])
    Y=c(P[,3],T[,3])
    Treat=c(rep('Control',n0),rep('Treatment',n1))
    
    #Fit 1, single baseline
    pre_mean1=X1
    post_mean=Y

    #Fitting Centered ANCOVA Model
    cpre_mean1=pre_mean1-mean(pre_mean1)
    model1 <- lm(post_mean ~Treat+cpre_mean1)
    pvalue1[i]=summary(model1)$coefficients[2,4]
    
    #Fit 2, screening and baseline
    pre_mean2=rowMeans(cbind(X1,X2))
    
    #Fitting Centered ANCOVA Model
    cpre_mean2=pre_mean2-mean(pre_mean2)
    model2 <- lm(post_mean ~Treat+cpre_mean2)
    pvalue2[i]=summary(model2)$coefficients[2,4]
  }
  
  power=c(sum(pvalue1<=0.05)/nsim,
          sum(pvalue2<=0.05)/nsim)
  
  return(power)  
}

#n0=n1=50
TypeI_rhoXY_0.5_rhoX_0.6_n50=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0,n0=50,n1=50)
TypeI_rhoXY_0.5_rhoX_0.7_n50=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0,n0=50,n1=50)
TypeI_rhoXY_0.5_rhoX_0.8_n50=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0,n0=50,n1=50)
TypeI_rhoXY_0.5_rhoX_0.9_n50=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0,n0=50,n1=50)

TypeI_rhoXY_0.6_rhoX_0.7_n50=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0,n0=50,n1=50)
TypeI_rhoXY_0.6_rhoX_0.8_n50=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0,n0=50,n1=50)
TypeI_rhoXY_0.6_rhoX_0.9_n50=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0,n0=50,n1=50)

TypeI_rhoXY_0.7_rhoX_0.8_n50=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0,n0=50,n1=50)
TypeI_rhoXY_0.7_rhoX_0.9_n50=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0,n0=50,n1=50)

TypeI_n50=rbind(TypeI_rhoXY_0.5_rhoX_0.6_n50,
                TypeI_rhoXY_0.5_rhoX_0.7_n50,
                TypeI_rhoXY_0.5_rhoX_0.8_n50,
                TypeI_rhoXY_0.5_rhoX_0.9_n50,
                TypeI_rhoXY_0.6_rhoX_0.7_n50,
                TypeI_rhoXY_0.6_rhoX_0.8_n50,
                TypeI_rhoXY_0.6_rhoX_0.9_n50,
                TypeI_rhoXY_0.7_rhoX_0.8_n50,
                TypeI_rhoXY_0.7_rhoX_0.9_n50)

TypeI_n50
saveRDS(TypeI_n50, file="TypeI_n50.Rds")

Power_rhoXY_0.5_rhoX_0.6_n50=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0.3,n0=50,n1=50)
Power_rhoXY_0.5_rhoX_0.7_n50=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0.3,n0=50,n1=50)
Power_rhoXY_0.5_rhoX_0.8_n50=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0.3,n0=50,n1=50)
Power_rhoXY_0.5_rhoX_0.9_n50=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0.3,n0=50,n1=50)

Power_rhoXY_0.6_rhoX_0.7_n50=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0.3,n0=50,n1=50)
Power_rhoXY_0.6_rhoX_0.8_n50=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0.3,n0=50,n1=50)
Power_rhoXY_0.6_rhoX_0.9_n50=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0.3,n0=50,n1=50)

Power_rhoXY_0.7_rhoX_0.8_n50=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0.3,n0=50,n1=50)
Power_rhoXY_0.7_rhoX_0.9_n50=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0.3,n0=50,n1=50)

Power_n50=rbind(Power_rhoXY_0.5_rhoX_0.6_n50,
                Power_rhoXY_0.5_rhoX_0.7_n50,
                Power_rhoXY_0.5_rhoX_0.8_n50,
                Power_rhoXY_0.5_rhoX_0.9_n50,
                Power_rhoXY_0.6_rhoX_0.7_n50,
                Power_rhoXY_0.6_rhoX_0.8_n50,
                Power_rhoXY_0.6_rhoX_0.9_n50,
                Power_rhoXY_0.7_rhoX_0.8_n50,
                Power_rhoXY_0.7_rhoX_0.9_n50)

Power_n50
saveRDS(Power_n50, file="Power_n50.Rds")

#n0=n1=75
TypeI_rhoXY_0.5_rhoX_0.6_n75=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0,n0=75,n1=75)
TypeI_rhoXY_0.5_rhoX_0.7_n75=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0,n0=75,n1=75)
TypeI_rhoXY_0.5_rhoX_0.8_n75=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0,n0=75,n1=75)
TypeI_rhoXY_0.5_rhoX_0.9_n75=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0,n0=75,n1=75)

TypeI_rhoXY_0.6_rhoX_0.7_n75=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0,n0=75,n1=75)
TypeI_rhoXY_0.6_rhoX_0.8_n75=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0,n0=75,n1=75)
TypeI_rhoXY_0.6_rhoX_0.9_n75=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0,n0=75,n1=75)

TypeI_rhoXY_0.7_rhoX_0.8_n75=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0,n0=75,n1=75)
TypeI_rhoXY_0.7_rhoX_0.9_n75=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0,n0=75,n1=75)

TypeI_n75=rbind(TypeI_rhoXY_0.5_rhoX_0.6_n75,
                      TypeI_rhoXY_0.5_rhoX_0.7_n75,
                      TypeI_rhoXY_0.5_rhoX_0.8_n75,
                      TypeI_rhoXY_0.5_rhoX_0.9_n75,
                      TypeI_rhoXY_0.6_rhoX_0.7_n75,
                      TypeI_rhoXY_0.6_rhoX_0.8_n75,
                      TypeI_rhoXY_0.6_rhoX_0.9_n75,
                      TypeI_rhoXY_0.7_rhoX_0.8_n75,
                      TypeI_rhoXY_0.7_rhoX_0.9_n75)

TypeI_n75
saveRDS(TypeI_n75, file="TypeI_n75.Rds")

Power_rhoXY_0.5_rhoX_0.6_n75=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0.3,n0=75,n1=75)
Power_rhoXY_0.5_rhoX_0.7_n75=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0.3,n0=75,n1=75)
Power_rhoXY_0.5_rhoX_0.8_n75=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0.3,n0=75,n1=75)
Power_rhoXY_0.5_rhoX_0.9_n75=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0.3,n0=75,n1=75)

Power_rhoXY_0.6_rhoX_0.7_n75=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0.3,n0=75,n1=75)
Power_rhoXY_0.6_rhoX_0.8_n75=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0.3,n0=75,n1=75)
Power_rhoXY_0.6_rhoX_0.9_n75=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0.3,n0=75,n1=75)

Power_rhoXY_0.7_rhoX_0.8_n75=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0.3,n0=75,n1=75)
Power_rhoXY_0.7_rhoX_0.9_n75=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0.3,n0=75,n1=75)

Power_n75=rbind(Power_rhoXY_0.5_rhoX_0.6_n75,
                      Power_rhoXY_0.5_rhoX_0.7_n75,
                      Power_rhoXY_0.5_rhoX_0.8_n75,
                      Power_rhoXY_0.5_rhoX_0.9_n75,
                      Power_rhoXY_0.6_rhoX_0.7_n75,
                      Power_rhoXY_0.6_rhoX_0.8_n75,
                      Power_rhoXY_0.6_rhoX_0.9_n75,
                      Power_rhoXY_0.7_rhoX_0.8_n75,
                      Power_rhoXY_0.7_rhoX_0.9_n75)

Power_n75
saveRDS(Power_n75, file="Power_n75.Rds")

#n0=n1=100
TypeI_rhoXY_0.5_rhoX_0.6_n100=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0,n0=100,n1=100)
TypeI_rhoXY_0.5_rhoX_0.7_n100=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0,n0=100,n1=100)
TypeI_rhoXY_0.5_rhoX_0.8_n100=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0,n0=100,n1=100)
TypeI_rhoXY_0.5_rhoX_0.9_n100=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0,n0=100,n1=100)

TypeI_rhoXY_0.6_rhoX_0.7_n100=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0,n0=100,n1=100)
TypeI_rhoXY_0.6_rhoX_0.8_n100=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0,n0=100,n1=100)
TypeI_rhoXY_0.6_rhoX_0.9_n100=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0,n0=100,n1=100)

TypeI_rhoXY_0.7_rhoX_0.8_n100=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0,n0=100,n1=100)
TypeI_rhoXY_0.7_rhoX_0.9_n100=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0,n0=100,n1=100)

TypeI_n100=rbind(TypeI_rhoXY_0.5_rhoX_0.6_n100,
                       TypeI_rhoXY_0.5_rhoX_0.7_n100,
                       TypeI_rhoXY_0.5_rhoX_0.8_n100,
                       TypeI_rhoXY_0.5_rhoX_0.9_n100,
                       TypeI_rhoXY_0.6_rhoX_0.7_n100,
                       TypeI_rhoXY_0.6_rhoX_0.8_n100,
                       TypeI_rhoXY_0.6_rhoX_0.9_n100,
                       TypeI_rhoXY_0.7_rhoX_0.8_n100,
                       TypeI_rhoXY_0.7_rhoX_0.9_n100)

TypeI_n100
saveRDS(TypeI_n100, file="TypeI_n100.Rds")

Power_rhoXY_0.5_rhoX_0.6_n100=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0.3,n0=100,n1=100)
Power_rhoXY_0.5_rhoX_0.7_n100=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0.3,n0=100,n1=100)
Power_rhoXY_0.5_rhoX_0.8_n100=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0.3,n0=100,n1=100)
Power_rhoXY_0.5_rhoX_0.9_n100=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0.3,n0=100,n1=100)

Power_rhoXY_0.6_rhoX_0.7_n100=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0.3,n0=100,n1=100)
Power_rhoXY_0.6_rhoX_0.8_n100=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0.3,n0=100,n1=100)
Power_rhoXY_0.6_rhoX_0.9_n100=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0.3,n0=100,n1=100)

Power_rhoXY_0.7_rhoX_0.8_n100=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0.3,n0=100,n1=100)
Power_rhoXY_0.7_rhoX_0.9_n100=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0.3,n0=100,n1=100)

Power_n100=rbind(Power_rhoXY_0.5_rhoX_0.6_n100,
                       Power_rhoXY_0.5_rhoX_0.7_n100,
                       Power_rhoXY_0.5_rhoX_0.8_n100,
                       Power_rhoXY_0.5_rhoX_0.9_n100,
                       Power_rhoXY_0.6_rhoX_0.7_n100,
                       Power_rhoXY_0.6_rhoX_0.8_n100,
                       Power_rhoXY_0.6_rhoX_0.9_n100,
                       Power_rhoXY_0.7_rhoX_0.8_n100,
                       Power_rhoXY_0.7_rhoX_0.9_n100)

Power_n100
saveRDS(Power_n100, file="Power_n100.Rds")
#n0=n1=125
TypeI_rhoXY_0.5_rhoX_0.6_n125=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0,n0=125,n1=125)
TypeI_rhoXY_0.5_rhoX_0.7_n125=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0,n0=125,n1=125)
TypeI_rhoXY_0.5_rhoX_0.8_n125=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0,n0=125,n1=125)
TypeI_rhoXY_0.5_rhoX_0.9_n125=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0,n0=125,n1=125)

TypeI_rhoXY_0.6_rhoX_0.7_n125=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0,n0=125,n1=125)
TypeI_rhoXY_0.6_rhoX_0.8_n125=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0,n0=125,n1=125)
TypeI_rhoXY_0.6_rhoX_0.9_n125=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0,n0=125,n1=125)

TypeI_rhoXY_0.7_rhoX_0.8_n125=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0,n0=125,n1=125)
TypeI_rhoXY_0.7_rhoX_0.9_n125=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0,n0=125,n1=125)

TypeI_n125=rbind(TypeI_rhoXY_0.5_rhoX_0.6_n125,
                       TypeI_rhoXY_0.5_rhoX_0.7_n125,
                       TypeI_rhoXY_0.5_rhoX_0.8_n125,
                       TypeI_rhoXY_0.5_rhoX_0.9_n125,
                       TypeI_rhoXY_0.6_rhoX_0.7_n125,
                       TypeI_rhoXY_0.6_rhoX_0.8_n125,
                       TypeI_rhoXY_0.6_rhoX_0.9_n125,
                       TypeI_rhoXY_0.7_rhoX_0.8_n125,
                       TypeI_rhoXY_0.7_rhoX_0.9_n125)

TypeI_n125
saveRDS(TypeI_n125, file="TypeI_n125.Rds")

Power_rhoXY_0.5_rhoX_0.6_n125=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0.3,n0=125,n1=125)
Power_rhoXY_0.5_rhoX_0.7_n125=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0.3,n0=125,n1=125)
Power_rhoXY_0.5_rhoX_0.8_n125=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0.3,n0=125,n1=125)
Power_rhoXY_0.5_rhoX_0.9_n125=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0.3,n0=125,n1=125)

Power_rhoXY_0.6_rhoX_0.7_n125=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0.3,n0=125,n1=125)
Power_rhoXY_0.6_rhoX_0.8_n125=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0.3,n0=125,n1=125)
Power_rhoXY_0.6_rhoX_0.9_n125=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0.3,n0=125,n1=125)

Power_rhoXY_0.7_rhoX_0.8_n125=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0.3,n0=125,n1=125)
Power_rhoXY_0.7_rhoX_0.9_n125=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0.3,n0=125,n1=125)

Power_n125=rbind(Power_rhoXY_0.5_rhoX_0.6_n125,
                       Power_rhoXY_0.5_rhoX_0.7_n125,
                       Power_rhoXY_0.5_rhoX_0.8_n125,
                       Power_rhoXY_0.5_rhoX_0.9_n125,
                       Power_rhoXY_0.6_rhoX_0.7_n125,
                       Power_rhoXY_0.6_rhoX_0.8_n125,
                       Power_rhoXY_0.6_rhoX_0.9_n125,
                       Power_rhoXY_0.7_rhoX_0.8_n125,
                       Power_rhoXY_0.7_rhoX_0.9_n125)

Power_n125
saveRDS(Power_n125, file="Power_n125.Rds")
#n0=n1=150
TypeI_rhoXY_0.5_rhoX_0.6_n150=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0,n0=150,n1=150)
TypeI_rhoXY_0.5_rhoX_0.7_n150=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0,n0=150,n1=150)
TypeI_rhoXY_0.5_rhoX_0.8_n150=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0,n0=150,n1=150)
TypeI_rhoXY_0.5_rhoX_0.9_n150=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0,n0=150,n1=150)

TypeI_rhoXY_0.6_rhoX_0.7_n150=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0,n0=150,n1=150)
TypeI_rhoXY_0.6_rhoX_0.8_n150=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0,n0=150,n1=150)
TypeI_rhoXY_0.6_rhoX_0.9_n150=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0,n0=150,n1=150)

TypeI_rhoXY_0.7_rhoX_0.8_n150=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0,n0=150,n1=150)
TypeI_rhoXY_0.7_rhoX_0.9_n150=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0,n0=150,n1=150)

TypeI_n150=rbind(TypeI_rhoXY_0.5_rhoX_0.6_n150,
                       TypeI_rhoXY_0.5_rhoX_0.7_n150,
                       TypeI_rhoXY_0.5_rhoX_0.8_n150,
                       TypeI_rhoXY_0.5_rhoX_0.9_n150,
                       TypeI_rhoXY_0.6_rhoX_0.7_n150,
                       TypeI_rhoXY_0.6_rhoX_0.8_n150,
                       TypeI_rhoXY_0.6_rhoX_0.9_n150,
                       TypeI_rhoXY_0.7_rhoX_0.8_n150,
                       TypeI_rhoXY_0.7_rhoX_0.9_n150)

TypeI_n150
saveRDS(TypeI_n150, file="TypeI_n150.Rds")

Power_rhoXY_0.5_rhoX_0.6_n150=power_RB_continuous(rhoXY=0.5,rhoX=0.6,nsim=20000,delta=0.3,n0=150,n1=150)
Power_rhoXY_0.5_rhoX_0.7_n150=power_RB_continuous(rhoXY=0.5,rhoX=0.7,nsim=20000,delta=0.3,n0=150,n1=150)
Power_rhoXY_0.5_rhoX_0.8_n150=power_RB_continuous(rhoXY=0.5,rhoX=0.8,nsim=20000,delta=0.3,n0=150,n1=150)
Power_rhoXY_0.5_rhoX_0.9_n150=power_RB_continuous(rhoXY=0.5,rhoX=0.9,nsim=20000,delta=0.3,n0=150,n1=150)

Power_rhoXY_0.6_rhoX_0.7_n150=power_RB_continuous(rhoXY=0.6,rhoX=0.7,nsim=20000,delta=0.3,n0=150,n1=150)
Power_rhoXY_0.6_rhoX_0.8_n150=power_RB_continuous(rhoXY=0.6,rhoX=0.8,nsim=20000,delta=0.3,n0=150,n1=150)
Power_rhoXY_0.6_rhoX_0.9_n150=power_RB_continuous(rhoXY=0.6,rhoX=0.9,nsim=20000,delta=0.3,n0=150,n1=150)

Power_rhoXY_0.7_rhoX_0.8_n150=power_RB_continuous(rhoXY=0.7,rhoX=0.8,nsim=20000,delta=0.3,n0=150,n1=150)
Power_rhoXY_0.7_rhoX_0.9_n150=power_RB_continuous(rhoXY=0.7,rhoX=0.9,nsim=20000,delta=0.3,n0=150,n1=150)

Power_n150=rbind(Power_rhoXY_0.5_rhoX_0.6_n150,
                       Power_rhoXY_0.5_rhoX_0.7_n150,
                       Power_rhoXY_0.5_rhoX_0.8_n150,
                       Power_rhoXY_0.5_rhoX_0.9_n150,
                       Power_rhoXY_0.6_rhoX_0.7_n150,
                       Power_rhoXY_0.6_rhoX_0.8_n150,
                       Power_rhoXY_0.6_rhoX_0.9_n150,
                       Power_rhoXY_0.7_rhoX_0.8_n150,
                       Power_rhoXY_0.7_rhoX_0.9_n150)

Power_n150
saveRDS(Power_n150, file="Power_n150.Rds")