#Clear the environment
rm(list=ls())

#Load necessary libraries
library(deSolve)
library(Deriv)
library(rootSolve)
library(pracma)
library(fitdistrplus)
library(DescTools)
library(mclust)

#Now we want something of the form Beta=a*f^b
#betaA=a*fA_dep^b
#betaB=a*fB_dep^b

fA_dep=0.01231144
fB_dep=0.009
betaA=2.475841e-6
betaB=7.216188e-7

b=log(betaB/betaA)/log(fB_dep/fA_dep)
a=betaA/fA_dep^b

#Use a simpler tradeoff function
Beta=function(fmax,fmin,fscale,f){
  fmax*f^fmin
}

#Set the working directory to wherever you downloaded the data files
#setwd("C:/Users/Jason/OneDrive - University of Pittsburgh/Documents/Thesis Chapter 2")

#Need to do this for the castration case
#Define the equations
{
  #For the one clone case
  I_deq=function(w,K,f,e,th,d,B,v,s,m,I,S,R,Z){
    B*S*Z-(d+v)*I
  }
  S_deq=function(w,K,f,e,th,d,B,v,s,m,I,S,R,Z){
    e*f*R*(S+th*I)-B*S*Z-d*S
  }
  R_deq=function(w,K,f,e,th,d,B,v,s,m,I,S,R,Z){
    w*R*(1-R/K)-f*R*(S+I)
  }
  Z_deq=function(w,K,f,e,th,d,B,v,s,m,I,S,R,Z){
    s*(d+v)*I-m*Z
  }
  
  #For the one clone, castration case
  I_deq_cast=function(w,K,f,e,d,B,v,s,m,I,S,R,Z){
    B*S*Z-(d+v)*I
  }
  S_deq_cast=function(w,K,f,e,d,B,v,s,m,I,S,R,Z){
    e*f*R*S-B*S*Z-d*S
  }
  R_deq_cast=function(w,K,f,e,d,B,v,s,m,I,S,R,Z){
    w*R*(1-R/K)-f*R*(S+I)
  }
  Z_deq_cast=function(w,K,f,e,d,B,v,s,m,I,S,R,Z){
    s*(d+v)*I-m*Z
  }
  
  #For the two clone case
  I_deq1=function(w,K,f,f2,e,th,d,B,B2,v,s,m,I,I2,S,S2,R,Z){
    B*S*Z-(d+v)*I
  }
  I_deq2=function(w,K,f,f2,e,th,d,B,B2,v,s,m,I,I2,S,S2,R,Z){
    B2*S2*Z-(d+v)*I2
  }
  S_deq1=function(w,K,f,f2,e,th,d,B,B2,v,s,m,I,I2,S,S2,R,Z){
    e*f*R*(S+th*I)-B*S*Z-d*S
  }
  S_deq2=function(w,K,f,f2,e,th,d,B,B2,v,s,m,I,I2,S,S2,R,Z){
    e*f2*R*(S2+th*I2)-B2*S2*Z-d*S2
  }
  R_deq1=function(w,K,f,f2,e,th,d,B,B2,v,s,m,I,I2,S,S2,R,Z){
    w*R*(1-R/K)-f*R*(S+I)-f2*R*(S2+I2)
  }
  Z_deq1=function(w,K,f,f2,e,th,d,B,B2,v,s,m,I,I2,S,S2,R,Z){
    s*(d+v)*(I+I2)-m*Z
  }
}

#Take symbolic derivatives and define jacobian juctions
{
  J11=Deriv(I_deq,"I")
  J12=Deriv(I_deq,"S")
  J13=Deriv(I_deq,"R")
  J14=Deriv(I_deq,"Z")
  J21=Deriv(S_deq,"I")
  J22=Deriv(S_deq,"S")
  J23=Deriv(S_deq,"R")
  J24=Deriv(S_deq,"Z")
  J31=Deriv(R_deq,"I")
  J32=Deriv(R_deq,"S")
  J33=Deriv(R_deq,"R")
  J34=Deriv(R_deq,"Z")
  J41=Deriv(Z_deq,"I")
  J42=Deriv(Z_deq,"S")
  J43=Deriv(Z_deq,"R")
  J44=Deriv(Z_deq,"Z")
  
  J11_cast=Deriv(I_deq_cast,"I")
  J12_cast=Deriv(I_deq_cast,"S")
  J13_cast=Deriv(I_deq_cast,"R")
  J14_cast=Deriv(I_deq_cast,"Z")
  J21_cast=Deriv(S_deq_cast,"I")
  J22_cast=Deriv(S_deq_cast,"S")
  J23_cast=Deriv(S_deq_cast,"R")
  J24_cast=Deriv(S_deq_cast,"Z")
  J31_cast=Deriv(R_deq_cast,"I")
  J32_cast=Deriv(R_deq_cast,"S")
  J33_cast=Deriv(R_deq_cast,"R")
  J34_cast=Deriv(R_deq_cast,"Z")
  J41_cast=Deriv(Z_deq_cast,"I")
  J42_cast=Deriv(Z_deq_cast,"S")
  J43_cast=Deriv(Z_deq_cast,"R")
  J44_cast=Deriv(Z_deq_cast,"Z")
  
  Jacobian_symbolic=list(c(J11,J12,J13,J14),c(J21,J22,J23,J24),c(J31,J32,J33,J34),c(J41,J42,J43,J44))
  Jacobian_symbolic_cast=list(c(J11_cast,J12_cast,J13_cast,J14_cast),c(J21_cast,J22_cast,J23_cast,J24_cast),c(J31_cast,J32_cast,J33_cast,J34_cast),c(J41_cast,J42_cast,J43_cast,J44_cast))
  
  Jacobian_numeric=function(w,K,f,e,th,d,B,v,s,m,I,S,R,Z){
    Jacobian_holder=array(NA,dim=c(4,4))
    for(i in 1:4){
      for(j in 1:4){
        Jacobian_holder[i,j]=Jacobian_symbolic[[i]][[j]](w,K,f,e,th,d,B,v,s,m,I,S,R,Z)
      }
    }
    return(Jacobian_holder)
  }
  
  Jacobian_numeric_cast=function(w,K,f,e,d,B,v,s,m,I,S,R,Z){
    Jacobian_holder=array(NA,dim=c(4,4))
    for(i in 1:4){
      for(j in 1:4){
        Jacobian_holder[i,j]=Jacobian_symbolic_cast[[i]][[j]](w,K,f,e,d,B,v,s,m,I,S,R,Z)
      }
    }
    return(Jacobian_holder)
  }
  
  J11_2=Deriv(I_deq1,"I")
  J12_2=Deriv(I_deq1,"I2")
  J13_2=Deriv(I_deq1,"S")
  J14_2=Deriv(I_deq1,"S2")
  J15_2=Deriv(I_deq1,"R")
  J16_2=Deriv(I_deq1,"Z")
  
  J21_2=Deriv(I_deq2,"I")
  J22_2=Deriv(I_deq2,"I2")
  J23_2=Deriv(I_deq2,"S")
  J24_2=Deriv(I_deq2,"S2")
  J25_2=Deriv(I_deq2,"R")
  J26_2=Deriv(I_deq2,"Z")
  
  J31_2=Deriv(S_deq1,"I")
  J32_2=Deriv(S_deq1,"I2")
  J33_2=Deriv(S_deq1,"S")
  J34_2=Deriv(S_deq1,"S2")
  J35_2=Deriv(S_deq1,"R")
  J36_2=Deriv(S_deq1,"Z")
  
  J41_2=Deriv(S_deq2,"I")
  J42_2=Deriv(S_deq2,"I2")
  J43_2=Deriv(S_deq2,"S")
  J44_2=Deriv(S_deq2,"S2")
  J45_2=Deriv(S_deq2,"R")
  J46_2=Deriv(S_deq2,"Z")
  
  J51_2=Deriv(R_deq1,"I")
  J52_2=Deriv(R_deq1,"I2")
  J53_2=Deriv(R_deq1,"S")
  J54_2=Deriv(R_deq1,"S2")
  J55_2=Deriv(R_deq1,"R")
  J56_2=Deriv(R_deq1,"Z")
  
  J61_2=Deriv(Z_deq1,"I")
  J62_2=Deriv(Z_deq1,"I2")
  J63_2=Deriv(Z_deq1,"S")
  J64_2=Deriv(Z_deq1,"S2")
  J65_2=Deriv(Z_deq1,"R")
  J66_2=Deriv(Z_deq1,"Z")
  
  Jacobian_symbolic_2=list(c(J11_2,J12_2,J13_2,J14_2,J15_2,J16_2),c(J21_2,J22_2,J23_2,J24_2,J25_2,J26_2),c(J31_2,J32_2,J33_2,J34_2,J35_2,J36_2),c(J41_2,J42_2,J43_2,J44_2,J45_2,J46_2),c(J51_2,J52_2,J53_2,J54_2,J55_2,J56_2),c(J61_2,J62_2,J63_2,J64_2,J65_2,J66_2))
  
  Jacobian_numeric_2=function(w,K,f,f2,e,th,d,B,B2,v,s,m,I,I2,S,S2,R,Z){
    if(!is.na(w*K*f*f2*e*th*d*B*B2*v*s*m*I*I2*S*S2*R*Z)){
    Jacobian_holder=array(NA,dim=c(6,6))
    for(i in 1:6){
      for(j in 1:6){
        Jacobian_holder[i,j]=Jacobian_symbolic_2[[i]][[j]](w,K,f,f2,e,th,d,B,B2,v,s,m,I,I2,S,S2,R,Z)
      }
    }
    }else{
      Jacobian_holder=rbind(c(1,0,0,0,0,0),rep(0,6),rep(0,6),rep(0,6),rep(0,6),rep(0,6))
    }

    return(Jacobian_holder)
  }
}

Equilibria=function(w,K,f,e,th,d,B,v,s,m){
  rbind(
  #Expressions for I,S,R,Z
  c(0,
    0,
    0,
    -(d*m)/(B*d*s + B*s*v),
    -(d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*d^2*s*w + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*K*e*f^2*th*(d*s + s*v)),
    -(B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*K*e*f^2*th*(d*s + s*v))),
  c(0,
    0,
    -(w*(d - K*e*f))/(K*e*f^2),
    m/(B*s),
    m/(B*s),
    m/(B*s)),
  c(0,
    K,
    d/(e*f),
  0,
  (d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*d^2*s*w + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*e*f*th*(d*s*w + s*v*w)) - (K*f*m - B*K*s*w)/(B*s*w),
  (B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*e*f*th*(d*s*w + s*v*w)) - (K*f*m - B*K*s*w)/(B*s*w)),
  c(0,
    0,
    0,
    -d/B,
    -(d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*d^2*s*w + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*K*e*f^2*m*th),
    -(B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*K*e*f^2*m*th))
  )
}

Equilibria_cast=function(w,K,f,e,d,B,v,s,m){
  rbind(
    #Expressions for I,S,R,Z
    c(0,
      -(m*(K*e*m*f^2 - B*K*e*s*w*f + B*d*s*w))/(s*(B^2*d*s*w + B^2*s*v*w + B*K*e*f^2*m)),
      0,
      0,
      -(d*m)/(B*d*s + B*s*v)),
    c(0,
      m/(B*s),
      0,
      -(w*(d - K*e*f))/(K*e*f^2),
      m/(B*s)),
    c(0,
      (B*K*d*s*w - K*f*m*v + B*K*s*v*w)/(K*e*m*f^2 + B*s*v*w + B*d*s*w),
      K,
      d/(e*f),
      0),
    c(0,
      -(B*s*w*d^2 + K*e*m*d*f^2 - B*K*e*s*w*d*f + B*s*v*w*d + K*e*m*v*f^2 - B*K*e*s*v*w*f)/(B^2*d*s*w + B^2*s*v*w + B*K*e*f^2*m),
      0,
      0,
      -d/B))
}
#Outputs as I, I2, S, S2, R, Z


#Check which, if any, equilibria are stable, including the two disease free boundary equilibria (not trivial)
Eigen_out=function(w,K,f,e,th,d,B,v,s,m){
  Equilibria_numeric=Equilibria(w,K,f,e,th,d,B,v,s,m)
  eigen_length=length(Equilibria_numeric)/4
  Max_eigen=array(NA,dim=eigen_length)
  for(i in 1:eigen_length){
    Max_eigen[i]=max(Re(eigen(Jacobian_numeric(w,K,f,e,th,d,B,v,s,m,Equilibria_numeric[1,i],Equilibria_numeric[2,i],Equilibria_numeric[3,i],Equilibria_numeric[4,i]))$values))
  }
  ifelse(length(which(Max_eigen<0))==1,return(Equilibria_numeric[,which(Max_eigen<0)]),return(length(which(Max_eigen<0))))
  #Outputs as I, S, R, Z
}

#Do the same for the castration case
Eigen_out_cast=function(w,K,f,e,d,B,v,s,m){
  Equilibria_numeric=Equilibria_cast(w,K,f,e,d,B,v,s,m)
  eigen_length=length(Equilibria_numeric)/4
  Max_eigen=array(NA,dim=eigen_length)
  for(i in 1:eigen_length){
    Max_eigen[i]=max(Re(eigen(Jacobian_numeric_cast(w,K,f,e,d,B,v,s,m,Equilibria_numeric[1,i],Equilibria_numeric[2,i],Equilibria_numeric[3,i],Equilibria_numeric[4,i]))$values))
  }
  ifelse(length(which(Max_eigen<0))==1,return(Equilibria_numeric[,which(Max_eigen<0)]),return(length(which(Max_eigen<0))))
  #Outputs as I, S, R, Z
}

#Equilibrium comes from the following diff eq
#dR=w*R*(1-R/K)-f*R*(S+I)
#dS=e*f*R*(S+th*I)-d*S-B*S*Z
#dI=B*S*Z-(d+v)*I
#dZ=s*(d+v)*I-m*Z
#where B (Beta) = f/fscale*log((fmax-fmin)/(fmax-f)) (eq. 2)
#Beta=function(fmax,fmin,fscale,f){
#  f/fscale*log((fmax-fmin)/(fmax-f))
#}


#The only stable, endemic equilibrium imported from results
#of Symbolic_prep_alt_tradeoff_AmNat
Rpartial<-function(w,K,f,e,th,d,B,v,s,m){
  (B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*e*f*th*(d*s*w + s*v*w)) - (K*f*m - B*K*s*w)/(B*s*w)
}

Spartial<-function(w,K,f,e,th,d,B,v,s,m){
  m/(B*s)
}

Ipartial<-function(w,K,f,e,th,d,B,v,s,m){
  -(B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*K*e*f^2*th*(d*s + s*v))
}

Zpartial<-function(w,K,f,e,th,d,B,v,s,m){
  -(B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*K*e*f^2*m*th)
}

#Write the castration equilibria separately to avoid division by 0 issues
#if thetas are imperfectly cancelled by the computer
Rcast<-function(w,K,f,e,d,B,v,s,m){
  (B*K*d*s*w - K*f*m*v + B*K*s*v*w)/(K*e*m*f^2 + B*s*v*w + B*d*s*w)
}

Scast<-function(w,K,f,e,d,B,v,s,m){
  m/(B*s)
}

Icast<-function(w,K,f,e,d,B,v,s,m){
  -(m*(K*e*m*f^2 - B*K*e*s*w*f + B*d*s*w))/(s*(B^2*d*s*w + B^2*s*v*w + B*K*e*f^2*m))
}

Zcast<-function(w,K,f,e,d,B,v,s,m){
  -(B*s*w*d^2 + K*e*m*d*f^2 - B*K*e*s*w*d*f + B*s*v*w*d + K*e*m*v*f^2 - B*K*e*s*v*w*f)/(B^2*d*s*w + B^2*s*v*w + B*K*e*f^2*m)
}

#Calculate a proxy for fitness according to NGM theory
#using fM (mutant trait) and fR (resident trait)
Fit_proxy=function(fM,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m){
  #Written simply, efR/(d+BZ)+efR*th/(d+v)*BZ/(d+BZ)
  e*fM*Rpartial(w,K,fR,e,th,d,Beta(fmax,fmin,fscale,fR),v,s,m)/(d+Beta(fmax,fmin,fscale,fM)*Zpartial(w,K,fR,e,th,d,Beta(fmax,fmin,fscale,fR),v,s,m))+
  e*fM*th*Rpartial(w,K,fR,e,th,d,Beta(fmax,fmin,fscale,fR),v,s,m)/(d+v)*Beta(fmax,fmin,fscale,fM)*Zpartial(w,K,fR,e,th,d,Beta(fmax,fmin,fscale,fR),v,s,m)/(d+Beta(fmax,fmin,fscale,fM)*Zpartial(w,K,fR,e,th,d,Beta(fmax,fmin,fscale,fR),v,s,m))
}

Fit_proxy(f_temp*1.02,f_temp,w_ex,K_temp,e_ex,th_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex)

Fit_proxy_cast=function(fM,fR,w,K,e,d,fmax,fmin,fscale,v,s,m){
  (e*fM*Rcast(w,K,fR,e,d,Beta(fmax,fmin,fscale,fR),v,s,m))/(d+Beta(fmax,fmin,fscale,fM)*Zcast(w,K,fR,e,d,Beta(fmax,fmin,fscale,fR),v,s,m))
}


#Sometimes useful to write equilibrium expressions this way.
Spartial2=function(w,K,f,e,th,d,B,v,s,m){
  m/(B*s)
}
Rpartial2=function(w,K,f,e,th,d,B,v,s,m){
  (B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*e*f*th*(d*s*w + s*v*w)) - (K*f*m - B*K*s*w)/(B*s*w)
}

Ipartial2=function(w,K,f,e,th,d,B,v,s,m){
  -(B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*K*e*f^2*th*(d*s + s*v))
}
Zpartial2=function(w,K,f,e,th,d,B,v,s,m){
  -(B*d^2*s*w - v*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) - d*(B^2*K^2*e^2*f^2*s^2*th^2*w^2 - 2*B^2*K*d*e*f*s^2*th*w^2 - 2*B^2*K*e*f*s^2*th*v*w^2 + B^2*d^2*s^2*w^2 + 2*B^2*d*s^2*v*w^2 + B^2*s^2*v^2*w^2 - 2*B*K^2*e^2*f^3*m*s*th^2*w + 2*B*K^2*e^2*f^3*m*s*th*w - 2*B*K*d*e*f^2*m*s*th*w + 2*B*K*d*e*f^2*m*s*w + 2*B*K*e*f^2*m*s*th*v*w + 2*B*K*e*f^2*m*s*v*w + K^2*e^2*f^4*m^2*th^2 - 2*K^2*e^2*f^4*m^2*th + K^2*e^2*f^4*m^2)^(1/2) + B*s*v^2*w + 2*B*d*s*v*w + K*d*e*f^2*m + K*e*f^2*m*v + K*d*e*f^2*m*th + K*e*f^2*m*th*v - B*K*d*e*f*s*th*w - B*K*e*f*s*th*v*w)/(2*B*K*e*f^2*m*th)
}

Hpartial2=function(w,K,f,e,th,d,B,v,s,m){
  Ipartial2(w,K,f,e,th,d,B,v,s,m)+Spartial2(w,K,f,e,th,d,B,v,s,m)
}

#Can use this simple one when we're in the stable endemic only region
prevpartial2=function(w,K,f,e,th,d,B,v,s,m){
  Ipartial2(w,K,f,e,th,d,B,v,s,m)/Hpartial2(w,K,f,e,th,d,B,v,s,m)
}

Fit_proxy2=function(fM,fR,BM,BR,w,K,e,th,d,v,s,m){
  #Written simply, (efR(d+v+th*B*Z))/((d+v)*(d+B*Z))
  (e*fM*Rpartial2(w,K,fR,e,th,d,BR,v,s,m)*(d+v)+
     e*fM*Rpartial2(w,K,fR,e,th,d,BR,v,s,m)*th*BM*Zpartial2(w,K,fR,e,th,d,BR,v,s,m))/((d+v)*(d+BM*Zpartial2(w,K,fR,e,th,d,BR,v,s,m)))
}


NGM_root=function(fM,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m){
  Fit_proxy(fM,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m)-1
}

NGM_root_cast=function(fM,fR,w,K,e,d,fmax,fmin,fscale,v,s,m){
  Fit_proxy_cast(fM,fR,w,K,e,d,fmax,fmin,fscale,v,s,m)-1
}

#Take key derivatives with respect to the invader's fitness proxy
NGM_difffM=Deriv(NGM_root,"fM")
NGM_diff2fM=Deriv(NGM_difffM,"fM")
NGM_difffR=Deriv(NGM_root,"fR")
NGM_diff2fR=Deriv(NGM_difffR,"fR")

NGM_difffM_cast=Deriv(NGM_root_cast,"fM")
NGM_diff2fM_cast=Deriv(NGM_difffM_cast,"fM")
NGM_difffR_cast=Deriv(NGM_root_cast,"fR")
NGM_diff2fR_cast=Deriv(NGM_difffR_cast,"fR")

#Write functions to quickly and easily evaluate these derivatives
#at a point when fR=fM
#Need to use the variable name m_parm instead of m here since
#m means something else for the uniroot.all function that will be applied
#to this function
Diff1fM_NGM=function(fM,w,K,e,th,d,fmax,fmin,fscale,v,s,m_parm){
  NGM_difffM(fM,fM,w,K,e,th,d,fmax,fmin,fscale,v,s,m_parm)
}

Diff2fM_NGM=function(fM,w,K,e,th,d,fmax,fmin,fscale,v,s,m){
  NGM_diff2fM(fM,fM,w,K,e,th,d,fmax,fmin,fscale,v,s,m)
}

Diff1fR_NGM=function(fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m){
  NGM_difffR(fR,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m)
}

Diff2fR_NGM=function(fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m){
  NGM_diff2fR(fR,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m)
}

Diff1fM_NGM_cast=function(fM,w,K,e,d,fmax,fmin,fscale,v,s,m_parm){
  NGM_difffM_cast(fM,fM,w,K,e,d,fmax,fmin,fscale,v,s,m_parm)
}

Diff2fM_NGM_cast=function(fM,w,K,e,d,fmax,fmin,fscale,v,s,m){
  NGM_diff2fM_cast(fM,fM,w,K,e,d,fmax,fmin,fscale,v,s,m)
}

Diff1fR_NGM_cast=function(fR,w,K,e,d,fmax,fmin,fscale,v,s,m){
  NGM_difffR_cast(fR,fR,w,K,e,d,fmax,fmin,fscale,v,s,m)
}

Diff2fR_NGM_cast=function(fR,w,K,e,d,fmax,fmin,fscale,v,s,m){
  NGM_diff2fR_cast(fR,fR,w,K,e,d,fmax,fmin,fscale,v,s,m)
}

#Return the f value of the CSS for a given parameter set
NGM_checkedroots=function(w_func,K_func,e_func,th_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func){
  #Find evolutionary singular points
  roots_storage=uniroot.all(Diff1fM_NGM,interval=c(d_func/(e_func*K_func),1e2),tol=1e-8,n=1e6,w=w_func,K=K_func,e=e_func,th=th_func,d=d_func,fmax=fmax_func,fmin=fmin_func,fscale=fscale_func,v=v_func,s=s_func,m_parm=m_func)
  fchecks=array(NA,dim=c(length(roots_storage),2))
  for(j in 1:length(roots_storage)){
    #Evaluate evolutionary stability
    fchecks[j,1]=Diff2fM_NGM(roots_storage[j],w_func,K_func,e_func,th_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)
    #Evaluate convergence stability
    fchecks[j,2]=Diff2fR_NGM(roots_storage[j],w_func,K_func,e_func,th_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)-Diff2fM_NGM(roots_storage[j],w_func,K_func,e_func,th_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)
  }
  #Only return CSSs
  return(roots_storage[which(fchecks[,1]<0&fchecks[,2]>0)])
}

#When I just want to return singular points
NGM_singulars=function(w_func,K_func,e_func,th_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func){
  #Find evolutionary singular points
  return(uniroot.all(Diff1fM_NGM,interval=c(d_func/(e_func*K_func),1e2),tol=1e-8,n=1e6,w=w_func,K=K_func,e=e_func,th=th_func,d=d_func,fmax=fmax_func,fmin=fmin_func,fscale=fscale_func,v=v_func,s=s_func,m_parm=m_func))
}

#Same for castration case
NGM_checkedroots_cast=function(w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func){
  roots_storage=uniroot.all(Diff1fM_NGM_cast,interval=c(d_func/(e_func*K_func),1e2),tol=1e-8,n=1e6,w=w_func,K=K_func,e=e_func,d=d_func,fmax=fmax_func,fmin=fmin_func,fscale=fscale_func,v=v_func,s=s_func,m_parm=m_func)
  fchecks=array(NA,dim=c(length(roots_storage),2))
  for(j in 1:length(roots_storage)){
    fchecks[j,1]=Diff2fM_NGM_cast(roots_storage[j],w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)
    fchecks[j,2]=Diff2fR_NGM_cast(roots_storage[j],w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)-Diff2fM_NGM_cast(roots_storage[j],w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)
  }
  return(roots_storage[which(fchecks[,1]<0&fchecks[,2]>0)])
}
NGM_singulars_cast=function(w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func){
  return(roots_storage=uniroot.all(Diff1fM_NGM_cast,interval=c(d_func/(e_func*K_func),1e2),tol=1e-8,n=1e6,w=w_func,K=K_func,e=e_func,d=d_func,fmax=fmax_func,fmin=fmin_func,fscale=fscale_func,v=v_func,s=s_func,m_parm=m_func))
}

#Set some standard parameters as well.
Standard_parms<-c(0.9,30,.18,.65,.05,0.02068322,0.0051666,2983.647,.05,100000,1.5)
w_ex=Standard_parms[1]
K_ex=Standard_parms[2]
e_ex=Standard_parms[3]
th_ex=Standard_parms[4]
d_ex=Standard_parms[5]

#Have to parameterize the tradeoff a little differently than normal bc it
#is a different function.
fmax=a
fmin=b
fscale=Standard_parms[8]
v_ex=Standard_parms[9]
s_ex=Standard_parms[10]
m_ex=Standard_parms[11]
parms_ex<-c(w_ex,K_ex,e_ex,th_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex)

#Normally 200
K_length=200
K_series=10^linspace(-3,log10(150),K_length)


#Can run these lengthy computations or read in their results below to save time,
#skiping this entire section.
{
#noncast_singulars=mapply(FUN=NGM_singulars,w_func=w_ex,K_func=K_series,e_func=e_ex,th_func=th_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)
#cast_singulars=mapply(FUN=NGM_singulars_cast,w_func=w_ex,K_func=K_series,e_func=e_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)

#Due to some numerical sensitivity issues, the solver can return some
#"singular points" with very similar values. So I've added a step
#that only counts something as a second singular point if it is at least 
#5% different from previous CSSes that were found.

#Accommodate up to 1000 singular points for a given value of K, but that
#of course won't happen.
#noncast_singular_points_temp<-array(NA,dim=c(length(K_series),1000))
#cast_singular_points<-array(NA,dim=c(length(K_series),1000))
#for(i in 1:length(K_series)){
#  temp_length=NA
#  temp_length=length(noncast_singulars[[i]])
#  for(j in 1:temp_length){
#    if(length(noncast_singulars[[i]])>0&&(j==1||(min(abs(noncast_singulars[[i]][j]-noncast_singular_points_temp[i,]),na.rm=T)>0.05*noncast_singular_points_temp[i,which.min(abs(noncast_singulars[[i]][j]-noncast_singular_points_temp[i,]))]))){
#      noncast_singular_points_temp[i,j]=sort(noncast_singulars[[i]])[j]       
#    }
#  }
#  temp_length_cast=NA
#  temp_length_cast=length(cast_singulars[[i]])
#  for(j in 1:temp_length_cast){
#    if(length(cast_singulars[[i]])>0&&(j==1||(min(abs(cast_singulars[[i]][j]-cast_singular_points[i,]),na.rm=T)>0.05*cast_singular_points[i,which.min(abs(cast_singulars[[i]][j]-cast_singular_points[i,]))]))){
#      cast_singular_points[i,j]=sort(cast_singulars[[i]])[j]       
#    }
#  }
#}

#noncast_singular_points_temp=noncast_singular_points_temp[,which(colSums(noncast_singular_points_temp,na.rm=T)>0)]
#cast_singular_points=cast_singular_points[,which(colSums(cast_singular_points,na.rm=T)>0)]

#Can see that the very low singular points are unfeasible.
#prevpartial2(w_ex,K_series[187:200],noncast_singular_points_temp[187:200,1],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,noncast_singular_points_temp[187:200,1]),v_ex,s_ex,m_ex)
#So eliminate them.
#noncast_singular_points_temp[187:200,1]<-NA

#noncast_singular_points=array(NA,dim=c(length(K_series),2))
#Reduce down to just two columns
#for(i in 1:length(K_series)){
#  if(length(which(!is.na(noncast_singular_points_temp[i,])))>0){
#    noncast_singular_points[i,1:length(which(!is.na(noncast_singular_points_temp[i,])))]=sort(noncast_singular_points_temp[i,which(!is.na(noncast_singular_points_temp[i,]))]) 
#  }
#}

#Save these to avoid redoing the computation, which can be lengthy
#saveRDS(noncast_singular_points,"noncast_singular_points_alt_tradeoff.RDS")
#saveRDS(cast_singular_points,"cast_singular_points_alt_tradeoff.RDS")
}
  
noncast_singular_points=readRDS("noncast_singular_points_alt_tradeoff.RDS")
cast_singular_points=readRDS("cast_singular_points_alt_tradeoff.RDS")

#Can read these back in next time the file is run, as long as you
#are using the same working directory you saved the RDS files in
ES_checker<-function(parms,f,K){
  w=parms[1]
  e=parms[3]
  th=parms[4]
  d=parms[5]
  fmax=parms[6]
  fmin=parms[7]
  fscale=parms[8]
  v=parms[9]
  s=parms[10]
  m=parms[11]
  #1 if ES, 0 if not
  return(ifelse(Diff2fM_NGM(f,w,K,e,th,d,fmax,fmin,fscale,v,s,m)<0,1,0))
}

CS_checker<-function(parms,f,K){
  w=parms[1]
  e=parms[3]
  th=parms[4]
  d=parms[5]
  fmax=parms[6]
  fmin=parms[7]
  fscale=parms[8]
  v=parms[9]
  s=parms[10]
  m=parms[11]
  #1 if ES, 0 if not
  return(ifelse(Diff2fR_NGM(f,w,K,e,th,d,fmax,fmin,fscale,v,s,m)-Diff2fM_NGM(f,w,K,e,th,d,fmax,fmin,fscale,v,s,m)>0,1,0))
}

CSS_checker<-function(parms,f,K){
  ES_checker(parms,f,K)*CS_checker(parms,f,K)
}


ES_checker_cast<-function(parms,f,K){
  w=parms[1]
  e=parms[3]
  d=parms[5]
  fmax=parms[6]
  fmin=parms[7]
  fscale=parms[8]
  v=parms[9]
  s=parms[10]
  m=parms[11]
  #1 if ES, 0 if not
  return(ifelse(Diff2fM_NGM_cast(f,w,K,e,d,fmax,fmin,fscale,v,s,m)<0,1,0))
}

CS_checker_cast<-function(parms,f,K){
  w=parms[1]
  e=parms[3]
  d=parms[5]
  fmax=parms[6]
  fmin=parms[7]
  fscale=parms[8]
  v=parms[9]
  s=parms[10]
  m=parms[11]
  #1 if ES, 0 if not
  return(ifelse(Diff2fR_NGM_cast(f,w,K,e,d,fmax,fmin,fscale,v,s,m)-Diff2fM_NGM_cast(f,w,K,e,d,fmax,fmin,fscale,v,s,m)>0,1,0))
}

CSS_checker_cast<-function(parms,f,K){
  ES_checker_cast(parms,f,K)*CS_checker_cast(parms,f,K)
}


Stability_status<-array(NA,dim=dim(noncast_singular_points))
Stability_status_cast<-array(NA,dim=length(cast_singular_points))
for(i in 1:dim(noncast_singular_points)[1]){
  for(j in 1:dim(noncast_singular_points)[2]){
    ES_checked<-NA
    CS_checked<-NA
    ES_checked<-ES_checker(parms_ex,noncast_singular_points[i,j],K_series[i])
    CS_checked<-CS_checker(parms_ex,noncast_singular_points[i,j],K_series[i])
    #Use 0 for neither, 1 for CS not ES, 2 for ES not CS, 3 for CSS
    if(!is.na(ES_checked*CS_checked)){
    if(ES_checked==1&CS_checked==1){
      Stability_status[i,j]<-3
    }
    if(ES_checked==1&CS_checked==0){
      Stability_status[i,j]<-2
    }
    if(ES_checked==0&CS_checked==1){
      Stability_status[i,j]<-1
    }
    if(ES_checked==0&CS_checked==0){
      Stability_status[i,j]<-0
    }
    }
  }
}
for(i in 1:length(cast_singular_points)){
    ES_checked_cast<-NA
    CS_checked_cast<-NA
    ES_checked_cast<-ES_checker_cast(parms_ex,cast_singular_points[i],K_series[i])
    CS_checked_cast<-CS_checker_cast(parms_ex,cast_singular_points[i],K_series[i])
    #Use 0 for neither, 1 for CS not ES, 2 for ES not CS, 3 for CSS
    if(!is.na(ES_checked_cast*CS_checked_cast)){
      if(ES_checked_cast==1&CS_checked_cast==1){
        Stability_status_cast[i]<-3
      }
      if(ES_checked_cast==1&CS_checked_cast==0){
        Stability_status_cast[i]<-2
      }
      if(ES_checked_cast==0&CS_checked_cast==1){
        Stability_status_cast[i]<-1
      }
      if(ES_checked_cast==0&CS_checked_cast==0){
        Stability_status_cast[i]<-0
      }
   }
}

Stability_status
Stability_status_cast

#H and prev functions
{
  Hpartial_final=function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
  temp=Eigen_out(w,K,f,e,th,d,Beta(fmax,fmin,fscale,f),v,s,m)
  return(sum(temp[1:2]))
  }
  Hcast_final=function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    temp=Eigen_out_cast(w,K,f,e,d,Beta(fmax,fmin,fscale,f),v,s,m)
    return(sum(temp[1:2]))
  }
  prevpartial_final=function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    temp=Eigen_out(w,K,f,e,th,d,Beta(fmax,fmin,fscale,f),v,s,m)
    return(temp[1]/sum(temp[1:2]))
  }
  prevcast_final=function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    temp=Eigen_out_cast(w,K,f,e,d,Beta(fmax,fmin,fscale,f),v,s,m)
    return(temp[1]/sum(temp[1:2]))
  }
  Ipartial_final=function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    Eigen_out(w,K,f,e,th,d,Beta(fmax,fmin,fscale,f),v,s,m)[1]
  }
  Icast_final=function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    Eigen_out_cast(w,K,f,e,d,Beta(fmax,fmin,fscale,f),v,s,m)[1]
  }
  Spartial_final=function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    Eigen_out(w,K,f,e,th,d,Beta(fmax,fmin,fscale,f),v,s,m)[2]
  }
  Scast_final=function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    Eigen_out_cast(w,K,f,e,d,Beta(fmax,fmin,fscale,f),v,s,m)[2]
  }
  Rpartial_final=function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    Eigen_out(w,K,f,e,th,d,Beta(fmax,fmin,fscale,f),v,s,m)[3]
  }
  Rcast_final=function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    Eigen_out_cast(w,K,f,e,d,Beta(fmax,fmin,fscale,f),v,s,m)[3]
  }
  Zpartial_final=function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    Eigen_out(w,K,f,e,th,d,Beta(fmax,fmin,fscale,f),v,s,m)[4]
  }
  Zcast_final=function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    Eigen_out_cast(w,K,f,e,d,Beta(fmax,fmin,fscale,f),v,s,m)[4]
  }
}

fnum1=0.02068322
fnum2=0.0051666
fden=2983.647

f_examples=linspace(fnum2*1.01,fnum1*.99,n=1000)
Beta_examples1=a*f_examples^b
Beta_examples2=f_examples/fden*log((fnum1-fnum2)/(fnum1-f_examples))

#Check the fitness gradient above the repellor, up to Beta = 10^12
f_Beta_1012=((10^12)/(fmax))^(1/fmin)
for(i in 1:length(K_series)){
  if(!is.na(noncast_singular_points[i,2])){
    temp_fs=linspace(noncast_singular_points[i,2]*1.01,f_Beta_1012,10000)
    temp_grads=mapply(Diff1fM_NGM,fM=temp_fs,w=w_ex,K=K_series[i],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
    print(c(i,min(temp_grads)))
  }
}
#So we can see that, above the repellor, we only find selection for
#increasing Beta, all the way up to 10^12.


#png("Ch2_app_alttradeoff.png",width=8,height=5*1.5,units="in",res=600)
#Alt tradeoff figure, Figure A2.
{
cex_smallest_text=.5
cex_minor_text=1
cex_major_text=1
lwd_minor=1
lwd_major=2.5
par(mfrow=c(2,2),mar=c(4,5,.1,.5),mgp=c(1.8,.4,0),cex.lab=cex_major_text,oma=c(0,0,0,0),cex.axis=cex_minor_text,bg="white",cex=1)
plot(f_examples,1e6*Beta_examples1,type="l",lty=2,xlab="",ylab=expression(atop("Transmission rate:",italic(beta)~" (L parasite"^"-1"~"day"^"-1"~"x 10"^"-6"~")")),lwd=lwd_major,col="darkgray")
mtext(side=1,line=3,expression(atop("Foraging rate:",italic(f)~"(L host"^"-1"~"day"^"-1"~")")),lwd=lwd_major)
points(f_examples,1e6*Beta_examples2,type="l",lwd=lwd_major)
points(c(fA_dep,fB_dep),1e6*c(betaA,betaB),pch=16,cex=2)
text(c(fA_dep+.003,fB_dep),c(0,3)+1e6*c(betaA,betaB),cex=cex_major_text,c(expression(atop("High","resistance")),expression(atop("Low","resistance"))))
text((par("usr")[2]-par("usr")[1])*.075+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_minor_text)
legend(0.005,15,lty=c(1,2),lwd=lwd_major,c("Main text tradeoff","Power law tradeoff"),col=c("black","darkgray"))

plot(K_series,log10(Beta(fmax,fmin,fscale,noncast_singular_points[,1])),type="l",ylim=range(log10(Beta(fmax,fmin,fscale,noncast_singular_points)),na.rm=T),lwd=lwd_major)
points(K_series,log10(Beta(fmax,fmin,fscale,noncast_singular_points[,2])),type="l",lty=2,lwd=lwd_major)
text((par("usr")[2]-par("usr")[1])*.075+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_minor_text)
mtext(side=1,line=3,expression(atop("Carrying capacity:",italic(K)~"("~mu~"g chl "~italic(a)~" L"^"-1"~")")))
legend(50,5,lty=c(1,2),lwd=lwd_major,c("CSS","Repellor"))

plot(K_series,prevpartial2(w_ex,K_series,noncast_singular_points[,1],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,noncast_singular_points[,1]),v_ex,s_ex,m_ex),type="l",xlab="",ylim=c(0.1,1.05),ylab=expression(atop("Prevalence:",italic(p)^"*"~"(unitless)")),lwd=lwd_major,lty=1)
points(K_series,prevpartial2(w_ex,K_series,noncast_singular_points[,2],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,noncast_singular_points[,2]),v_ex,s_ex,m_ex),type="l",xlab="",ylim=c(0.1,1.05),ylab=expression(atop("Prevalence:",italic(p)^"*"~"(unitless)")),lwd=lwd_major,lty=2)
text((par("usr")[2]-par("usr")[1])*.075+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_minor_text)
mtext(side=1,line=3,expression(atop("Carrying capacity:",italic(K)~"("~mu~"g chl "~italic(a)~" L"^"-1"~")")))

plot(K_series,Hpartial2(w_ex,K_series,noncast_singular_points[,1],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,noncast_singular_points[,1]),v_ex,s_ex,m_ex),type="l",xlab="",ylab=expression(atop("Host density:",italic(H)^"*"~"(Hosts L"^"-1"~")")),lwd=lwd_major,lty=1)
points(K_series,Hpartial2(w_ex,K_series,noncast_singular_points[,2],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,noncast_singular_points[,2]),v_ex,s_ex,m_ex),type="l",xlab="",ylab=expression(atop("Host density:",italic(H)^"*"~"(Hosts L"^"-1"~")")),lwd=lwd_major,lty=2)
text((par("usr")[2]-par("usr")[1])*.075+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"D",cex=cex_minor_text)
mtext(side=1,line=3,expression(atop("Carrying capacity:",italic(K)~"("~mu~"g chl "~italic(a)~" L"^"-1"~")")))
}
#dev.off()

