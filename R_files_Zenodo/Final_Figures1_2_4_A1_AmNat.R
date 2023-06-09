#Clear the environment
rm(list=ls())

#Note that Greek letters in the manuscript are notated differently here in the
#code thus:
#d in the code corresponds to delta in the manuscript.
#s in the code corresponds to sigma in the manuscript.
#B in the code corresponds to beta in the manuscript
#th in the code corresponds to theta in the manuscript

#Load necessary libraries
library(deSolve)
library(Deriv)
library(rootSolve)
library(pracma)
library(fitdistrplus)
library(DescTools)
library(mclust)

#Set the working directory to wherever you downloaded the data files
#setwd("C:/Users/Jason/OneDrive - University of Pittsburgh/Documents/Thesis Chapter 2")

#Define the differential equations
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

#Take symbolic derivatives and define jacobian functions
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

#Give expressions for the equilibria. These were found using
#Symbolic_prep_AmNat.m.
Equilibria=function(w,K,f,e,th,d,B,v,s,m){
  #Each column is a different equilibrium.
  #Different rows hold expressions for I,S,R,Z
  rbind(
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

#Same formatting but for the castration case.
Equilibria_cast=function(w,K,f,e,d,B,v,s,m){
  rbind(
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

#Calculate all equilibria in the two clone case

#This function is very large so we have already saved it as an R object
#and read it in here to save space. In the output, each of the 6 rows is a state 
#variable, ordered I,I2,S,S2,R,Z, and each of 12 columns corresponds to an
#equilibrium
Equilibria2=readRDS("Equilibria2.RDS")
#Outputs as I, I2, S, S2, R, Z


#Check which, if any, equilibria are stable, including the two disease free boundary equilibria
Which_stable=function(w,K,f,e,th,d,B,v,s,m){
  Equilibria_numeric=Equilibria(w,K,f,e,th,d,B,v,s,m)
  eigen_length=length(Equilibria_numeric)/4
  Max_eigen=array(NA,dim=eigen_length)
  for(i in 1:eigen_length){
    Max_eigen[i]=max(Re(eigen(Jacobian_numeric(w,K,f,e,th,d,B,v,s,m,Equilibria_numeric[1,i],Equilibria_numeric[2,i],Equilibria_numeric[3,i],Equilibria_numeric[4,i]))$values))
  }
  ifelse(length(which(Max_eigen<0))==1,return(which(Max_eigen<0)),return(paste("Error: number of stable equilibria is",length(which(Max_eigen<0)),sep="")))
}

#Return the stable equilibrium values
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
Which_stable_cast=function(w,K,f,e,d,B,v,s,m){
  Equilibria_numeric=Equilibria_cast(w,K,f,e,d,B,v,s,m)
  eigen_length=length(Equilibria_numeric)/4
  Max_eigen=array(NA,dim=eigen_length)
  for(i in 1:eigen_length){
    Max_eigen[i]=max(Re(eigen(Jacobian_numeric_cast(w,K,f,e,d,B,v,s,m,Equilibria_numeric[1,i],Equilibria_numeric[2,i],Equilibria_numeric[3,i],Equilibria_numeric[4,i]))$values))
  }
  ifelse(length(which(Max_eigen<0))==1,return(which(Max_eigen<0)),return(paste("Error: number of stable equilibria is",length(which(Max_eigen<0)),sep="")))
}

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


#Do the same for the two clone case
#Need to also check for feasibility as there can occasionally, in the two clone case,
#be an equilibrium that is stable but unfeasible
Eigen_out2=function(w,K,f,f2,e,th,d,B,B2,v,s,m){
  Equilibria_numeric=Equilibria2(w,K,f,f2,e,th,d,B,B2,v,s,m)
  eigen_length=length(Equilibria_numeric)/6
  Max_eigen=array(NA,dim=eigen_length)
  Feasible=array(NA,dim=eigen_length)
  for(i in 1:eigen_length){
    Max_eigen[i]=max(Re(eigen(Jacobian_numeric_2(w,K,f,f2,e,th,d,B,B2,v,s,m,Equilibria_numeric[1,i],Equilibria_numeric[2,i],Equilibria_numeric[3,i],Equilibria_numeric[4,i],Equilibria_numeric[5,i],Equilibria_numeric[6,i]))$values))
    Feasible[i]=min(Equilibria_numeric[,i])>=0
  }
  
  ifelse(length(which(Max_eigen<0&Feasible))==1,return(Equilibria_numeric[,which(Max_eigen<0&Feasible)]),return(length(which(Max_eigen<0&Feasible))))
  #Outputs as I, I2, S, S2, R, Z
}

#Get two clone trait values.
fA_dep=0.01231144
fB_dep=0.009
betaA=2.475841e-6
betaB=7.216188e-7

#Set some standard parameters as well.
Standard_parms<-c(0.9,30,.18,.65,.05,0.02068322,0.0051666,2983.647,.05,100000,1.5)
w_ex=Standard_parms[1]
K_ex=Standard_parms[2]
e_ex=Standard_parms[3]
th_ex=Standard_parms[4]
d_ex=Standard_parms[5]
fmax=Standard_parms[6]
fmin=Standard_parms[7]
fscale=Standard_parms[8]
v_ex=Standard_parms[9]
s_ex=Standard_parms[10]
m_ex=Standard_parms[11]

#Numerically check that the stable equilibrium is always 6 (endemic), 3 (only hosts), or 2 (only resources) in the parameter
#range we use in the case without castration

#Note that this takes about one hour to run so it's commented out unless you
#really wish to run it.
#check_length=1000
#K_extremes=linspace(14,150,check_length)
#f_extremes=linspace(0.0065,0.02065,check_length)
#check_results=array(NA,dim=c(check_length,check_length))
#tic()
#for(k in 1:check_length){
#  check_results[,k]=mapply(Which_stable,w=w_ex,K=K_extremes[k],f=f_extremes,e=e_ex,th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,f_extremes),v=v_ex,s=s_ex,m=m_ex)
#  print(k)
#}
#toc()
#Indeed, the only stable equilibria were 2, 3, and 6
#unique(c(check_results))

#Calculate two genotype outcomes along a range of carrying capacities.
K_example=linspace(25,150,200)
Two_clone_test=array(NA,dim=c(length(K_example),6))
for(i in 1:length(K_example)){
  Two_clone_test[i,]=Eigen_out2(w_ex,K_example[i],fA_dep,fB_dep,e_ex,th_ex,d_ex,betaA,betaB,v_ex,s_ex,m_ex)
  if((Two_clone_test[i,1]==2)&(Two_clone_test[i,2]==2)&(Two_clone_test[i,3]==2)&(Two_clone_test[i,4]==2)&(Two_clone_test[i,5]==2)&(Two_clone_test[i,6]==2)){
    Two_clone_test[i,5]=NA
    Two_clone_test[i,6]=NA
    Two_clone_test[i,c(1,3)]=Eigen_out(w_ex,K_example[i],fA_dep,e_ex,th_ex,d_ex,betaA,v_ex,s_ex,m_ex)[1:2]
    Two_clone_test[i,c(2,4)]=Eigen_out(w_ex,K_example[i],fB_dep,e_ex,th_ex,d_ex,betaB,v_ex,s_ex,m_ex)[1:2]
  }
}

#Also calculate outcomes for each genotype alone.
SoloA=array(NA,dim=c(length(K_example),4))
SoloB=array(NA,dim=c(length(K_example),4))
for(i in 1:length(K_example)){
  SoloA[i,]=Eigen_out(w_ex,K_example[i],fA_dep,e_ex,th_ex,d_ex,betaA,v_ex,s_ex,m_ex)
  SoloB[i,]=Eigen_out(w_ex,K_example[i],fB_dep,e_ex,th_ex,d_ex,betaB,v_ex,s_ex,m_ex)
}

#Calculte frequencies, prevalences, and total host densities from what
#was already calculated above.
freq_calculator=array(NA,dim=c(length(K_example),2))
prev_calculator=array(NA,dim=c(length(K_example),2))
H_calculator=array(NA,dim=c(length(K_example),2))
for(i in 1:length(K_example)){
  if(!is.na(Two_clone_test[i,6])){
    freq_calculator[i,1]=(Two_clone_test[i,1]+Two_clone_test[i,3])/(Two_clone_test[i,1]+Two_clone_test[i,2]+Two_clone_test[i,3]+Two_clone_test[i,4])
    freq_calculator[i,2]=(Two_clone_test[i,2]+Two_clone_test[i,4])/(Two_clone_test[i,1]+Two_clone_test[i,2]+Two_clone_test[i,3]+Two_clone_test[i,4]) 
  }else{
    freq_calculator[i,1]=1
    freq_calculator[i,2]=1
  }

}

#Get indices for which one genotype is dominant, they coexist,
#or there is alternative stable states
A_dominant_1_index=1:(min(which(Two_clone_test[,2]>0))-1)
coex_index=which(Two_clone_test[,1]*Two_clone_test[,2]&!is.na(Two_clone_test[,6]))
B_dominant_index=which(Two_clone_test[,1]==0)
ASS_index=which(is.na(Two_clone_test[,6]))
A_dominant_2_index=(max(ASS_index)+1):length(K_example)


#Read in files that can be reproduced in the Final_FiguresA3_A4_AmNat.R file
#from simulations. I have attached all .RDS files so you don't have to run
#Final_FiguresA3_A4_AmNat.R if you do not want to.
freq_av=readRDS("freq_av.RDS")
prev_av=readRDS("prev_av.RDS")
H_av=readRDS("H_av.RDS")

prev_av_solo=readRDS("prev_av_solo.RDS")
H_av_solo=readRDS("H_av_solo.RDS")

#Define some points for plotting polygon and points on Figure 2.
start_point=which.min(freq_av)
K_start=K_example[start_point]

H_line0=c(H_av_solo[start_point,2],H_av[start_point],H_av_solo[start_point,1])
prev_line0=c(prev_av_solo[start_point,2],prev_av[start_point],prev_av_solo[start_point,1])

H_line=c(H_av_solo[200,2],H_av[200],H_av_solo[200,1])
prev_line=c(prev_av_solo[200,2],prev_av[200],prev_av_solo[200,1])
beta_line=c(betaB,freq_av[200]*betaA+(1-freq_av[200])*betaB,betaA)

H_lm=lm(H_line~beta_line)

library(betareg)
prev_curve=betareg(prev_line~beta_line)

new_data=data.frame(beta_line=freq_av[start_point]*betaA+(1-freq_av[start_point])*betaB)
new_data2=data.frame(beta_line=linspace(min(beta_line)/2,max(beta_line)*1.5,1e4))
predict_H=predict(H_lm,new_data)
predict_prev=predict(prev_curve,new_data)
predict_prev_pretty=predict(prev_curve,new_data2)

summary(H_lm)
int=6.216e1
slope=-7.505e6

polygon_Ks=c(K_example[which.min(freq_av)],200,200,K_example[which.min(freq_av)])
polygon_freqs=c(min(freq_av),min(freq_av),max(freq_av),max(freq_av))
polygon_prevs=c(prev_av[which.min(freq_av)],prev_av[which.min(freq_av)],1,1)
polygon_Hs=c(H_av[which.min(freq_av)],H_av[which.min(freq_av)],max(c(H_av,rowSums(Two_clone_test[,1:4]))),max(c(H_av,rowSums(Two_clone_test[,1:4]))))

#Ch2_Mesocosm.png. Makes manuscript Figure 2.
#For all figures, we have commented out the portion that exports the figure as 
#a .png file. Instead, the figure will simply display in Rstudio.


#pdf("AmNat_Fig2.pdf",width=14,height=9,family="ArialMT",useDingbats=FALSE)
{
  cex_smallest_text=.5*1.2
  cex_minor_text=1*1.2
  cex_major_text=1.5*1.2
  cex_big_points=3
  lwd_minor=1.5*1.2
  lwd_major=2*1.2
  dens_choice=10
  par(mar=c(5,8.5,1,.75),cex.lab=cex_major_text,cex.axis=cex_major_text,cex=1,oma=c(0,0,0,0),font.lab=2,xaxs="i")
  
  m<-rbind(c(1,1,2,2,3,3),c(6,4,4,5,5,7))
  layout(m)
  
  aesthetic_stop=which(diff(freq_calculator[,1])==1)

  
  plot(K_example[1:aesthetic_stop],freq_calculator[1:aesthetic_stop,1],type="l",ylab="",lty=2,ylim=c(0,1.05),xlab="",lwd=lwd_major,xlim=range(K_example)*c(1,1.02))
  points(K_example[(aesthetic_stop+1):length(K_example)],freq_calculator[(aesthetic_stop+1):length(K_example),1],type="l",ylab="",lty=2,ylim=c(0,1.05),xlab="",lwd=lwd_major,xlim=range(K_example)*c(1,1.02))
  points(K_example[ASS_index],0*K_example[ASS_index],lty=2,lwd=lwd_major,type="l")
  points(K_example,freq_av,type="l",lwd=lwd_major)
  mtext(side=2,expression(atop("Freq. of low","resistance genotype")),cex=cex_major_text,line=2.6)
  polygon(polygon_Ks,polygon_freqs,col=rgb(0,0,1,0.25,maxColorValue=1),border=rgb(0,0,1,0.5,maxColorValue=1),fillOddEven="non-zero")
  arrows(x0=120,y0=0.52,x1=120,y1=0.65,col=rgb(0,0,1,0.5,maxColorValue=1),lwd=lwd_major*2)
  text(120,.45,expression(atop("Resistance","is futile")),cex=cex_major_text,col=rgb(0,0,1,0.5,maxColorValue=1),adj=c(0.5,0.5))  
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_major_text)
  par(xpd=NA)
  legend(-6,-0.1,c("Low resistance","Evolving","High resistance"),pch=c(15,21,24),pt.bg=c("black","darkgrey","white"),cex=cex_big_points*0.8,bty="n")
  par(xpd=F)
  
  poly_height=0.03
  polygon(x=c(0,K_example[max(A_dominant_1_index)],K_example[max(A_dominant_1_index)],0),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",border=NA)
  polygon(x=c(K_example[min(coex_index)],K_example[max(coex_index)],K_example[max(coex_index)],K_example[min(coex_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=0,col="black",border=NA)
  polygon(x=c(K_example[min(coex_index)],K_example[max(coex_index)],K_example[max(coex_index)],K_example[min(coex_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",border=NA)
  polygon(x=c(K_example[min(B_dominant_index)],K_example[max(B_dominant_index)],K_example[max(B_dominant_index)],K_example[min(B_dominant_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=0,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(ASS_index)],K_example[max(ASS_index)],K_example[max(ASS_index)],K_example[min(ASS_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=-45,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(ASS_index)],K_example[max(ASS_index)],K_example[max(ASS_index)],K_example[min(ASS_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=45,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(A_dominant_2_index)],200,200,K_example[min(A_dominant_2_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",lty=1,border=NA)
  
  legend(85,0.22,c("Simulation","Equilibrium"),lty=c(1,2),cex=cex_major_text,bty="n")

    
  plot(K_example,(Two_clone_test[,1]+Two_clone_test[,2])/(Two_clone_test[,1]+Two_clone_test[,2]+Two_clone_test[,3]+Two_clone_test[,4])*Two_clone_test[,5]/Two_clone_test[,5],type="l",ylim=c(0,1),ylab="",lwd=lwd_major,lty=2,xlab="",xlim=range(K_example)*c(1,1.02))
  mtext(side=1,expression("Carrying capacity:"~italic(K)~"("~mu~"g chl"~italic(a)~"L"^"-1"~")"),line=4.1,cex=cex_major_text)
  points(K_example,SoloA[,1]/(SoloA[,1]+SoloA[,2]),type="l",col="darkgray",lty=2,lwd=lwd_major)
  points(K_example,SoloB[,1]/(SoloB[,1]+SoloB[,2]),type="l",col="darkgray",lty=3,lwd=lwd_major)
  points(K_example[ASS_index],Two_clone_test[ASS_index,1]/(Two_clone_test[ASS_index,1]+Two_clone_test[ASS_index,3]),type="l",lwd=lwd_major,lty=2)
  points(K_example[ASS_index],Two_clone_test[ASS_index,2]/(Two_clone_test[ASS_index,2]+Two_clone_test[ASS_index,4]),type="l",lwd=lwd_major,lty=2)
  points(K_example,(Two_clone_test[,1]+Two_clone_test[,2])/(Two_clone_test[,1]+Two_clone_test[,2]+Two_clone_test[,3]+Two_clone_test[,4])*Two_clone_test[,5]/Two_clone_test[,5],type="l",ylim=c(0,1),ylab="",lwd=lwd_major,lty=2,xlab=expression("Carrying capacity:"~italic(K)~"("~mu~"g chl"~italic(a)~"L"^"-1"~")"))
  points(K_example,prev_av,type="l",lwd=lwd_major)
  mtext(side=2,expression(atop("Prevalence:",italic(p)~" (unitless)")),line=2.6,cex=cex_major_text)
  points(c(K_start,150),c(prev_av_solo[start_point,1],prev_av_solo[200,1]),cex=cex_big_points,pch=15,col="black")
  points(c(K_start,150),c(prev_av_solo[start_point,2],prev_av_solo[200,2]),pch=24,cex=cex_big_points,bg="white",lwd=lwd_major)
  points(c(K_start,150),c(prev_av[start_point],prev_av[200]),pch=21,bg="darkgrey",lwd=lwd_major,cex=cex_big_points)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_major_text)  
  text(80,0.1,expression(atop("High resistance","alone")),cex=cex_major_text,col="gray30")
  text(55,0.8,expression(atop("Low resistance","alone")),cex=cex_major_text,col="gray30")
  polygon(x=c(0,K_example[max(A_dominant_1_index)],K_example[max(A_dominant_1_index)],0),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",border=NA)
  polygon(x=c(K_example[min(coex_index)],K_example[max(coex_index)],K_example[max(coex_index)],K_example[min(coex_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=0,col="black",border=NA)
  polygon(x=c(K_example[min(coex_index)],K_example[max(coex_index)],K_example[max(coex_index)],K_example[min(coex_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",border=NA)
  polygon(x=c(K_example[min(B_dominant_index)],K_example[max(B_dominant_index)],K_example[max(B_dominant_index)],K_example[min(B_dominant_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=0,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(ASS_index)],K_example[max(ASS_index)],K_example[max(ASS_index)],K_example[min(ASS_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=-45,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(ASS_index)],K_example[max(ASS_index)],K_example[max(ASS_index)],K_example[min(ASS_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=45,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(A_dominant_2_index)],200,200,K_example[min(A_dominant_2_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",lty=1,border=NA)
    
  plot(K_example,(Two_clone_test[,1]+Two_clone_test[,2]+Two_clone_test[,3]+Two_clone_test[,4])*Two_clone_test[,5]/Two_clone_test[,5],type="l",ylim=range(c(SoloA[,1]+SoloA[,2],SoloB[,1]+SoloB[,2])),ylab="",lwd=2,lty=2,xlab="",xlim=range(K_example)*c(1,1.02))
  points(K_example,(SoloA[,1]+SoloA[,2]),type="l",col="darkgray",lty=2,lwd=2)
  points(K_example,(SoloB[,1]+SoloB[,2]),type="l",col="darkgray",lty=3,lwd=2)
  mtext(side=2,expression(atop("Host density:",italic(H)~"(Hosts L"^"-1"~")")),line=2.6,cex=cex_major_text)
  points(K_example[ASS_index],(Two_clone_test[ASS_index,1]+Two_clone_test[ASS_index,3]),type="l",lwd=2,lty=2)
  points(K_example[ASS_index],(Two_clone_test[ASS_index,2]+Two_clone_test[ASS_index,4]),type="l",lwd=2,lty=2)
  points(K_example,(Two_clone_test[,1]+Two_clone_test[,2]+Two_clone_test[,3]+Two_clone_test[,4])*Two_clone_test[,5]/Two_clone_test[,5],type="l",ylim=range(c(SoloA[,1]+SoloA[,2],SoloB[,1]+SoloB[,2])),ylab="",lwd=2,lty=2,xlab="")
  points(K_example,H_av,type="l",lwd=2)
  points(c(K_start,150),c(H_av_solo[start_point,1],H_av_solo[200,1]),cex=cex_big_points,pch=15,col="black")
  points(c(K_start,150),c(H_av_solo[start_point,2],H_av_solo[200,2]),pch=24,cex=cex_big_points,bg="white",lwd=lwd_major)
  points(c(K_start,150),c(H_av[start_point],H_av[200]),pch=21,bg="darkgrey",lwd=lwd_major,cex=cex_big_points)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_major_text)
  polygon(x=c(0,K_example[max(A_dominant_1_index)],K_example[max(A_dominant_1_index)],0),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",border=NA)
  polygon(x=c(K_example[min(coex_index)],K_example[max(coex_index)],K_example[max(coex_index)],K_example[min(coex_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=0,col="black",border=NA)
  polygon(x=c(K_example[min(coex_index)],K_example[max(coex_index)],K_example[max(coex_index)],K_example[min(coex_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",border=NA)
  polygon(x=c(K_example[min(B_dominant_index)],K_example[max(B_dominant_index)],K_example[max(B_dominant_index)],K_example[min(B_dominant_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=0,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(ASS_index)],K_example[max(ASS_index)],K_example[max(ASS_index)],K_example[min(ASS_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=-45,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(ASS_index)],K_example[max(ASS_index)],K_example[max(ASS_index)],K_example[min(ASS_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=45,col="black",lty=1,border=NA)
  polygon(x=c(K_example[min(A_dominant_2_index)],200,200,K_example[min(A_dominant_2_index)]),y=c(-10,-10,(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3],(par("usr")[4]-par("usr")[3])*poly_height+par("usr")[3]),density=20,angle=90,col="black",lty=1,border=NA)
    
  plot(1e6*new_data2$beta_line,as.numeric(predict_prev_pretty),xlab="",ylab="",type="l",lwd=lwd_major)
  mtext(side=2,expression(atop("Prevalence:",italic(p)~" (unitless)")),line=2.6,cex=cex_major_text)
  mtext(side=1,expression("Transmission rate:"~italic(beta)~"(L parasite"^"-1"~" day"^"-1"~" x 10"^"-6"~")"),line=4.1,adj=-.3,cex=cex_major_text)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.97+par("usr")[3],"D",cex=cex_major_text)
  points(beta_line[c(1,3)]*1e6,prev_line[c(1,3)],pch=c(24,15),bg=c("white","black"),cex=cex_big_points,col="black")
  points(beta_line[2]*1e6,prev_line[2],pch=21,bg="darkgrey",cex=cex_big_points,col="black")
  arrows(x0=new_data$beta_line*1e6,y0=predict_prev,x1=beta_line[2]*1e6,y1=prev_line[2],col=rgb(0,0,1,0.5,maxColorValue=1),lwd=lwd_major*2)
  points(new_data$beta_line*1e6,predict_prev,pch=21,bg=rgb(0,0,1,0.5,maxColorValue=1),cex=cex_big_points*1.2)
  text(2.5,0.775,"Consequences",cex=cex_major_text,col=rgb(0,0,1,0.5,maxColorValue=1))
  
  
  #Todo: Figure out why this host density plot is ALL wrong
  plot(1e6*new_data2$beta_line,new_data2$beta_line*slope+int,xlab="",ylab="",type="l",lwd=lwd_major)
  mtext(side=2,expression(atop("Host density:",italic(H)~" (hosts L"^"-1"~")")),line=2.6,cex=cex_major_text)
  text((par("usr")[2]-par("usr")[1])*.95+par("usr")[1],(par("usr")[4]-par("usr")[3])*.97+par("usr")[3],"E",cex=cex_major_text)
  points(beta_line[c(1,3)]*1e6,H_line[c(1,3)],pch=c(24,15),bg=c("white","black"),cex=cex_big_points,col="black")
  points(beta_line[2]*1e6,H_line[2],pch=21,bg="darkgrey",cex=cex_big_points,col="black")
  arrows(x0=as.numeric(new_data*1e6),y0=predict_H,x1=beta_line[2]*1e6,y1=H_line[2],col=rgb(0,0,1,0.5,maxColorValue=1),lwd=lwd_major*2)
  points(new_data*1e6,predict_H,pch=21,bg=rgb(0,0,1,0.5,maxColorValue=1),cex=cex_big_points*1.2)
  text(2.5,49.5,"Consequences",cex=cex_major_text,col=rgb(0,0,1,0.5,maxColorValue=1))
}
#dev.off()

  #Define the tradeoff function linking foraging rate, f, and 
  #transmission rate, beta.
  Beta=function(fmax,fmin,fscale,f){
    f/fscale*log((fmax-fmin)/(fmax-f))
  }
  
  #The only stable, endemic equilibrium imported from results
  #of Symbolic_prep_AmNat
  Rpartial<-function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    (d^2*s*w*log(-(fmax - fmin)/(f - fmax)) - v*(d^2*s^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + s^2*v^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + K^2*e^2*f^2*fscale^2*m^2 + 2*d*s^2*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K^2*e^2*f^2*fscale^2*m^2*th + K^2*e^2*f^2*fscale^2*m^2*th^2 + K^2*e^2*f^2*s^2*th^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*d*e*f*s^2*th*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*e*f*s^2*th*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 + 2*K^2*e^2*f^2*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) - 2*K^2*e^2*f^2*fscale*m*s*th^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*d*e*f*fscale*m*s*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*v*w*log(-(fmax - fmin)/(f - fmax)) - 2*K*d*e*f*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))^(1/2) - d*(d^2*s^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + s^2*v^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + K^2*e^2*f^2*fscale^2*m^2 + 2*d*s^2*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K^2*e^2*f^2*fscale^2*m^2*th + K^2*e^2*f^2*fscale^2*m^2*th^2 + K^2*e^2*f^2*s^2*th^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*d*e*f*s^2*th*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*e*f*s^2*th*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 + 2*K^2*e^2*f^2*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) - 2*K^2*e^2*f^2*fscale*m*s*th^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*d*e*f*fscale*m*s*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*v*w*log(-(fmax - fmin)/(f - fmax)) - 2*K*d*e*f*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))^(1/2) + s*v^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*d*s*v*w*log(-(fmax - fmin)/(f - fmax)) + K*d*e*f*fscale*m + K*e*f*fscale*m*v + K*d*e*f*fscale*m*th + K*e*f*fscale*m*th*v - K*d*e*f*s*th*w*log(-(fmax - fmin)/(f - fmax)) - K*e*f*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))/(2*e*f*th*log(-(fmax - fmin)/(f - fmax))*(d*s*w + s*v*w)) - (K*fscale*m - K*s*w*log(-(fmax - fmin)/(f - fmax)))/(s*w*log(-(fmax - fmin)/(f - fmax)))
  }
  
  Spartial<-function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    (fscale*m)/(f*s*log(-(fmax - fmin)/(f - fmax)))
  }
  
  Ipartial<-function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    -(d^2*s*w*log(-(fmax - fmin)/(f - fmax)) - v*(d^2*s^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + s^2*v^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + K^2*e^2*f^2*fscale^2*m^2 + 2*d*s^2*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K^2*e^2*f^2*fscale^2*m^2*th + K^2*e^2*f^2*fscale^2*m^2*th^2 + K^2*e^2*f^2*s^2*th^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*d*e*f*s^2*th*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*e*f*s^2*th*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 + 2*K^2*e^2*f^2*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) - 2*K^2*e^2*f^2*fscale*m*s*th^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*d*e*f*fscale*m*s*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*v*w*log(-(fmax - fmin)/(f - fmax)) - 2*K*d*e*f*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))^(1/2) - d*(d^2*s^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + s^2*v^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + K^2*e^2*f^2*fscale^2*m^2 + 2*d*s^2*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K^2*e^2*f^2*fscale^2*m^2*th + K^2*e^2*f^2*fscale^2*m^2*th^2 + K^2*e^2*f^2*s^2*th^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*d*e*f*s^2*th*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*e*f*s^2*th*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 + 2*K^2*e^2*f^2*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) - 2*K^2*e^2*f^2*fscale*m*s*th^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*d*e*f*fscale*m*s*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*v*w*log(-(fmax - fmin)/(f - fmax)) - 2*K*d*e*f*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))^(1/2) + s*v^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*d*s*v*w*log(-(fmax - fmin)/(f - fmax)) + K*d*e*f*fscale*m + K*e*f*fscale*m*v + K*d*e*f*fscale*m*th + K*e*f*fscale*m*th*v - K*d*e*f*s*th*w*log(-(fmax - fmin)/(f - fmax)) - K*e*f*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))/(2*K*e*f^2*th*log(-(fmax - fmin)/(f - fmax))*(d*s + s*v))
  }
  
  Zpartial<-function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    -(d^2*s*w*log(-(fmax - fmin)/(f - fmax)) - v*(d^2*s^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + s^2*v^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + K^2*e^2*f^2*fscale^2*m^2 + 2*d*s^2*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K^2*e^2*f^2*fscale^2*m^2*th + K^2*e^2*f^2*fscale^2*m^2*th^2 + K^2*e^2*f^2*s^2*th^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*d*e*f*s^2*th*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*e*f*s^2*th*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 + 2*K^2*e^2*f^2*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) - 2*K^2*e^2*f^2*fscale*m*s*th^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*d*e*f*fscale*m*s*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*v*w*log(-(fmax - fmin)/(f - fmax)) - 2*K*d*e*f*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))^(1/2) - d*(d^2*s^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + s^2*v^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 + K^2*e^2*f^2*fscale^2*m^2 + 2*d*s^2*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K^2*e^2*f^2*fscale^2*m^2*th + K^2*e^2*f^2*fscale^2*m^2*th^2 + K^2*e^2*f^2*s^2*th^2*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*d*e*f*s^2*th*w^2*log(-(fmax - fmin)/(f - fmax))^2 - 2*K*e*f*s^2*th*v*w^2*log(-(fmax - fmin)/(f - fmax))^2 + 2*K^2*e^2*f^2*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) - 2*K^2*e^2*f^2*fscale*m*s*th^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*d*e*f*fscale*m*s*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*v*w*log(-(fmax - fmin)/(f - fmax)) - 2*K*d*e*f*fscale*m*s*th*w*log(-(fmax - fmin)/(f - fmax)) + 2*K*e*f*fscale*m*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))^(1/2) + s*v^2*w*log(-(fmax - fmin)/(f - fmax)) + 2*d*s*v*w*log(-(fmax - fmin)/(f - fmax)) + K*d*e*f*fscale*m + K*e*f*fscale*m*v + K*d*e*f*fscale*m*th + K*e*f*fscale*m*th*v - K*d*e*f*s*th*w*log(-(fmax - fmin)/(f - fmax)) - K*e*f*s*th*v*w*log(-(fmax - fmin)/(f - fmax)))/(2*K*e*f^2*m*th*log(-(fmax - fmin)/(f - fmax)))
  }
  
  #Write the castration equilibria separately to avoid division by 0 issues
  #if thetas are imperfectly cancelled by the computer
  Rcast<-function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    (K*d*s*w*log(-(fmax - fmin)/(f - fmax)) - K*fscale*m*v + K*s*v*w*log(-(fmax - fmin)/(f - fmax)))/(s*v*w*log(-(fmax - fmin)/(f - fmax)) + d*s*w*log(-(fmax - fmin)/(f - fmax)) + K*e*f*fscale*m)
  }
  
  Scast<-function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    (fscale*m)/(f*s*log(-(fmax - fmin)/(f - fmax)))
  }
  
  Icast<-function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    -(m*(d*fscale*s*w*log(-(fmax - fmin)/(f - fmax)) + K*e*f*fscale^2*m - K*e*f*fscale*s*w*log(-(fmax - fmin)/(f - fmax))))/(s*(d*f*s*w*log(-(fmax - fmin)/(f - fmax))^2 + f*s*v*w*log(-(fmax - fmin)/(f - fmax))^2 + K*e*f^2*fscale*m*log(-(fmax - fmin)/(f - fmax))))
  }
  
  Zcast<-function(w,K,f,e,d,fmax,fmin,fscale,v,s,m){
    -(d^2*fscale*s*w*log(-(fmax - fmin)/(f - fmax)) + K*d*e*f*fscale^2*m + K*e*f*fscale^2*m*v + d*fscale*s*v*w*log(-(fmax - fmin)/(f - fmax)) - K*d*e*f*fscale*s*w*log(-(fmax - fmin)/(f - fmax)) - K*e*f*fscale*s*v*w*log(-(fmax - fmin)/(f - fmax)))/(d*f*s*w*log(-(fmax - fmin)/(f - fmax))^2 + f*s*v*w*log(-(fmax - fmin)/(f - fmax))^2 + K*e*f^2*fscale*m*log(-(fmax - fmin)/(f - fmax)))
  }
  
  #Calculate a proxy for fitness according to NGM theory
  #using fM (mutant trait) and fR (resident trait)
  Fit_proxy=function(fM,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m){
    #Written simply, (efR(d+v+th*B*Z))/((d+v)*(d+B*Z))
    (e*fM*Rpartial(w,K,fR,e,th,d,fmax,fmin,fscale,v,s,m)*(d+v)+
       e*fM*Rpartial(w,K,fR,e,th,d,fmax,fmin,fscale,v,s,m)*th*fM/fscale*log((fmax-fmin)/(fmax-fM))*Zpartial(w,K,fR,e,th,d,fmax,fmin,fscale,v,s,m))/((d+v)*(d+fM/fscale*log((fmax-fmin)/(fmax-fM))*Zpartial(w,K,fR,e,th,d,fmax,fmin,fscale,v,s,m)))
  }
  
  #Same but for castration case.
  Fit_proxy_cast=function(fM,fR,w,K,e,d,fmax,fmin,fscale,v,s,m){
    (e*fM*Rcast(w,K,fR,e,d,fmax,fmin,fscale,v,s,m))/(d+fM/fscale*log((fmax-fmin)/(fmax-fM))*Zcast(w,K,fR,e,d,fmax,fmin,fscale,v,s,m))
  }
  
  
#It can also be useful to have these endemic expressions with B instead of
#f and the tradeoff parameters.
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

Beta=function(fmax,fmin,fscale,f){
  f/fscale*log((fmax-fmin)/(fmax-f))  
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
  roots_storage=uniroot.all(Diff1fM_NGM,interval=c(fmin_func*1.00001,fmax_func*.99999),tol=1e-8,n=1e6,w=w_func,K=K_func,e=e_func,th=th_func,d=d_func,fmax=fmax_func,fmin=fmin_func,fscale=fscale_func,v=v_func,s=s_func,m_parm=m_func)
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
  return(uniroot.all(Diff1fM_NGM,interval=c(fmin_func*1.00001,fmax_func*.99999),tol=1e-8,n=1e6,w=w_func,K=K_func,e=e_func,th=th_func,d=d_func,fmax=fmax_func,fmin=fmin_func,fscale=fscale_func,v=v_func,s=s_func,m_parm=m_func))
}

#Same for castration case
NGM_checkedroots_cast=function(w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func){
  roots_storage=uniroot.all(Diff1fM_NGM_cast,interval=c(fmin_func*1.00001,fmax_func*.99999),tol=1e-8,n=1e6,w=w_func,K=K_func,e=e_func,d=d_func,fmax=fmax_func,fmin=fmin_func,fscale=fscale_func,v=v_func,s=s_func,m_parm=m_func)
  fchecks=array(NA,dim=c(length(roots_storage),2))
  for(j in 1:length(roots_storage)){
    fchecks[j,1]=Diff2fM_NGM_cast(roots_storage[j],w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)
    fchecks[j,2]=Diff2fR_NGM_cast(roots_storage[j],w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)-Diff2fM_NGM_cast(roots_storage[j],w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func)
  }
  return(roots_storage[which(fchecks[,1]<0&fchecks[,2]>0)])
}
NGM_singulars_cast=function(w_func,K_func,e_func,d_func,fmax_func,fmin_func,fscale_func,v_func,s_func,m_func){
  return(roots_storage=uniroot.all(Diff1fM_NGM_cast,interval=c(fmin_func*1.00001,fmax_func*.99999),tol=1e-8,n=1e6,w=w_func,K=K_func,e=e_func,d=d_func,fmax=fmax_func,fmin=fmin_func,fscale=fscale_func,v=v_func,s=s_func,m_parm=m_func))
}

#Define a K gradient along which to find the roots.
K_length=200
K_series=linspace(14,150,K_length)

#Can run these lengthy computations or read in their results below to save time
{
#noncast_singulars=mapply(FUN=NGM_singulars,w_func=w_ex,K_func=K_series,e_func=e_ex,th_func=th_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)
#cast_singulars=mapply(FUN=NGM_singulars_cast,w_func=w_ex,K_func=K_series,e_func=e_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)

#noncast_singular_points<-array(NA,dim=c(length(K_series),10))
#for(i in 1:length(K_series)){
#  temp_length=NA
#  temp_length=length(noncast_singulars[[i]])
#  for(j in 1:temp_length){
#    noncast_singular_points[i,j]=sort(noncast_singulars[[i]])[j] 
#  }
#}

#noncast_singular_points=t(noncast_singulars)
#cast_singular_points=t(cast_singulars)
}

#Save these to avoid redoing the computation, which can be lengthy
#saveRDS(noncast_singular_points,"noncast_singular_points.RDS")
#saveRDS(cast_singular_points,"cast_singular_points.RDS")

noncast_singular_points=readRDS("noncast_singular_points.RDS")
cast_singular_points=readRDS("cast_singular_points.RDS")

CSS=cbind(K_series,noncast_singular_points[,c(2,4)])
CSS_cast=cbind(K_series,cast_singular_points[,2])

#Can read these back in next time the file is run, as long as you
#are using the same working directory you saved the RDS files in

#Can check the precise stability properties of a given singular point
#with the following functions.
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


CS_checker_numeric<-function(parms,f,K){
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
  Diff2fR_NGM(f,w,K,e,th,d,fmax,fmin,fscale,v,s,m)-Diff2fM_NGM(f,w,K,e,th,d,fmax,fmin,fscale,v,s,m)
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

CS_checker_cast_numeric<-function(parms,f,K){
  w=parms[1]
  e=parms[3]
  d=parms[5]
  fmax=parms[6]
  fmin=parms[7]
  fscale=parms[8]
  v=parms[9]
  s=parms[10]
  m=parms[11]
  
  Diff2fR_NGM_cast(f,w,K,e,d,fmax,fmin,fscale,v,s,m)-Diff2fM_NGM_cast(f,w,K,e,d,fmax,fmin,fscale,v,s,m)
}

CSS_checker_cast<-function(parms,f,K){
  ES_checker_cast(parms,f,K)*CS_checker_cast(parms,f,K)
}

parms_ex<-c(w_ex,K_ex,e_ex,th_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex)

Stability_status<-array(NA,dim=dim(noncast_singular_points))
Stability_status_cast<-array(NA,dim=dim(cast_singular_points))
CS_status<-array(NA,dim=dim(noncast_singular_points))
CS_status_cast<-array(NA,dim=dim(noncast_singular_points))
for(i in 1:dim(noncast_singular_points)[1]){
  for(j in 1:dim(noncast_singular_points)[2]){
    ES_checked<-NA
    CS_checked<-NA
    ES_checked<-ES_checker(parms_ex,noncast_singular_points[i,j],K_series[i])
    CS_checked<-CS_checker(parms_ex,noncast_singular_points[i,j],K_series[i])
    
    CS_status[i,j]<-CS_checker_numeric(parms_ex,noncast_singular_points[i,j],K_series[i])
    
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
for(i in 1:dim(cast_singular_points)[1]){
  for(j in 1:dim(cast_singular_points)[2]){
    ES_checked_cast<-NA
    CS_checked_cast<-NA
    ES_checked_cast<-ES_checker_cast(parms_ex,cast_singular_points[i,j],K_series[i])
    CS_checked_cast<-CS_checker_cast(parms_ex,cast_singular_points[i,j],K_series[i])
    
    CS_status_cast[i,j]<-CS_checker_cast_numeric(parms_ex,cast_singular_points[i,j],K_series[i])
    #Use 0 for neither, 1 for CS not ES, 2 for ES not CS, 3 for CSS
    if(!is.na(ES_checked_cast*CS_checked_cast)){
      if(ES_checked_cast==1&CS_checked_cast==1){
        Stability_status_cast[i,j]<-3
      }
      if(ES_checked_cast==1&CS_checked_cast==0){
        Stability_status_cast[i,j]<-2
      }
      if(ES_checked_cast==0&CS_checked_cast==1){
        Stability_status_cast[i,j]<-1
      }
      if(ES_checked_cast==0&CS_checked_cast==0){
        Stability_status_cast[i,j]<-0
      }
   }
}
}

#We can see the stability properties here. There is a very low
#singular point that corresponds to biologically-unreasonably low
#values of transmission rate. Other than that, one can see a single
#CSS and then a CSS that splits from a repellor for this case
Stability_status

#And a single CSS for this case.
Stability_status_cast

#Can easily see that the bottom singular point is fixed and the same for noncast_singular_points
#and cast_singular_points, as well as unreasonably low.
Beta(fmax,fmin,fscale,unique(c(cast_singular_points[,1],noncast_singular_points[,1])))

#H and prev functions which will take into account stability and ensure
#that the equilibrium we are working with everywhere, which we think is the
#stable one, is in fact the stable one.
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



equilibria=array(NA,dim=c(length(K_series)))
equilibria_cast=array(NA,dim=c(length(K_series)))
for(i in 1:length(K_series)){
  equilibria[i]=prevpartial_final(w_ex,K_series[i],noncast_singular_points[i,1],e_ex,th_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex)  
  equilibria_cast[i]=prevcast_final(w_ex,K_series[i],cast_singular_points[i,1],e_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex)  
}

#Disease does not persist at these very low foraging rates so we will ignore them.

#Can find an interesting case of evolutionary branching:
parms_branch=c(0.9,102,0.25,1,0.05,0.0039,0.0024,1e4,0.037,100000,0.5)
test_singulars=NGM_singulars(w_func=0.9,K_func=102,e_func=0.25,th_func=1,d_func=0.05,fmax_func=0.0039,fmin_func=0.0024,fscale_func=1e4,v_func=0.037,s_func=100000,m_func=0.5)

CS_checker(parms_branch,test_singulars[1],102)
ES_checker(parms_branch,test_singulars[1],102)

CS_checker(parms_branch,test_singulars[2],102)
ES_checker(parms_branch,test_singulars[2],102)

CS_checker(parms_branch,test_singulars[3],102)
ES_checker(parms_branch,test_singulars[3],102)

CS_checker(parms_branch,test_singulars[4],102)
ES_checker(parms_branch,test_singulars[4],102)

#Shows that the very lowest one is not feasible.
prevpartial_final(w=0.9,K=102,f=test_singulars[1],e=0.25,th=1,d=0.05,fmax=0.0039,fmin=0.0024,fscale=1e4,v=0.037,s=100000,m=0.5)
prevpartial_final(w=0.9,K=102,f=test_singulars[2],e=0.25,th=1,d=0.05,fmax=0.0039,fmin=0.0024,fscale=1e4,v=0.037,s=100000,m=0.5)
prevpartial_final(w=0.9,K=102,f=test_singulars[3],e=0.25,th=1,d=0.05,fmax=0.0039,fmin=0.0024,fscale=1e4,v=0.037,s=100000,m=0.5)
prevpartial_final(w=0.9,K=102,f=test_singulars[4],e=0.25,th=1,d=0.05,fmax=0.0039,fmin=0.0024,fscale=1e4,v=0.037,s=100000,m=0.5)

#Sort the CSSes, for plotting purposes, in the case where alternate CSSes arise.
CSS_sorted=array(NA,dim=dim(CSS))
for(i in 1:length(CSS[,1])){
  CSS_sorted[i,1]=CSS[i,1]
  CSS_sorted[i,2]=-sort(-CSS[i,c(2,3)])[1]
  CSS_sorted[i,3]=-sort(-CSS[i,c(2,3)])[2]
}

CSS_bottom=CSS_sorted[which(!is.na(CSS_sorted[,3])),c(1,3)]

CSS_repellor=noncast_singular_points[,3]

example_bottom=min(cast_singular_points[,2],na.rm=T)
example_middle=min(CSS_sorted[,2],na.rm=T)
example_top=max(CSS_sorted[,2],na.rm=T)

Beta(fmax,fmin,fscale,example_middle)
Beta(fmax,fmin,fscale,example_top)

prev_cast=mapply(prevcast_final,w=w_ex,K=CSS_cast[,1],f=CSS_cast[,2],e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
H_cast=mapply(Hcast_final,w=w_ex,K=CSS_cast[,1],f=CSS_cast[,2],e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
S_cast=mapply(Scast_final,w=w_ex,K=CSS_cast[,1],f=CSS_cast[,2],e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)


bottom_prev_cast=mapply(prevcast_final,w=w_ex,K=CSS_cast[,1],f=example_bottom,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
middle_prev_cast0=mapply(prevcast_final,w=w_ex,K=CSS_cast[,1],f=example_middle,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
top_prev_cast0=mapply(prevcast_final,w=w_ex,K=CSS_cast[,1],f=example_top,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)

Icast(w_ex,K_series[40],example_top,e_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex)
Scast(w_ex,K_series[40],example_top,e_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex)

bottom_H_cast=mapply(Hcast_final,w=w_ex,K=CSS_cast[,1],f=example_bottom,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
middle_H_cast0=mapply(Hcast_final,w=w_ex,K=CSS_cast[,1],f=example_middle,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
top_H_cast0=mapply(Hcast_final,w=w_ex,K=CSS_cast[,1],f=example_top,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)

bottom_prev_cast[which(bottom_H_cast==0)]=0
middle_prev_cast0[which(middle_H_cast0==0)]=0
top_prev_cast0[which(top_H_cast0==0)]=0

#Check to ensure these are using the stable equilibrium we expect for the 
#castration case
#number 3 (no hosts), number 4 (hosts but no parasites), number 2 (endemic)
mapply(Which_stable_cast,w=w_ex,K=K_series,f=CSS_cast[,2],e=e_ex,d=d_ex,
              B=Beta(fmax,fmin,fscale,CSS_cast[,2]),v=v_ex,s=s_ex,m=m_ex)
mapply(Which_stable_cast,w=w_ex,K=K_series,f=example_bottom,e=e_ex,d=d_ex,
              B=Beta(fmax,fmin,fscale,example_bottom),v=v_ex,s=s_ex,m=m_ex)
mapply(Which_stable_cast,w=w_ex,K=K_series,f=example_middle,e=e_ex,d=d_ex,
              B=Beta(fmax,fmin,fscale,example_middle),v=v_ex,s=s_ex,m=m_ex)
mapply(Which_stable_cast,w=w_ex,K=K_series,f=example_top,e=e_ex,d=d_ex,
              B=Beta(fmax,fmin,fscale,example_top),v=v_ex,s=s_ex,m=m_ex)

#Save which K values need simulations for medium f and for high f
middle_Ks=which(is.na(middle_prev_cast0))
top_Ks=which(is.na(top_prev_cast0))

bottom_S_cast=mapply(Scast_final,w=w_ex,K=CSS_cast[,1],f=example_bottom,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
middle_S_cast0=mapply(Scast_final,w=w_ex,K=CSS_cast[,1],f=example_middle,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
top_S_cast0=mapply(Scast_final,w=w_ex,K=CSS_cast[,1],f=example_top,e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)

#Import simulation results for the two gray curves
#that are unstable in the castration case.
middle_cast=readRDS("Fig1_outputs_middle.RDS")
top_cast=readRDS("Fig1_outputs_top.RDS")

#Overwrite the curves that have stability issues with simulation results
middle_prev_cast=middle_prev_cast0
middle_H_cast=middle_H_cast0
middle_S_cast=middle_S_cast0

top_prev_cast=top_prev_cast0
top_H_cast=top_H_cast0
top_S_cast=top_S_cast0

middle_prev_cast[which(is.na(middle_prev_cast0))]=middle_cast[which(!is.na(middle_cast[,1])),1]
middle_H_cast[which(is.na(middle_prev_cast0))]=middle_cast[which(!is.na(middle_cast[,2])),2]
middle_S_cast[which(is.na(middle_prev_cast0))]=middle_cast[which(!is.na(middle_cast[,3])),3]
top_prev_cast[which(is.na(top_prev_cast0))]=top_cast[which(!is.na(top_cast[,1])),1]
top_H_cast[which(is.na(top_prev_cast0))]=top_cast[which(!is.na(top_cast[,1])),2]
top_S_cast[which(is.na(top_prev_cast0))]=top_cast[which(!is.na(top_cast[,1])),3]

prev_repellor=mapply(prevpartial_final,w=w_ex,K=CSS[which(!is.na(noncast_singular_points[,3])),1],f=CSS_repellor[which(!is.na(CSS_repellor))],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
H_repellor=mapply(Hpartial_final,w=w_ex,K=CSS[which(!is.na(noncast_singular_points[,3])),1],f=CSS_repellor[which(!is.na(CSS_repellor))],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
S_repellor=mapply(Spartial_final,w=w_ex,K=CSS[which(!is.na(noncast_singular_points[,3])),1],f=CSS_repellor[which(!is.na(CSS_repellor))],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)
K_repellor=CSS[which(!is.na(noncast_singular_points[,3])),1]


polygon_xs=c(CSS[which.min(CSS_sorted[,2]),1],300,300,CSS[which.min(CSS_sorted[,2]),1])
f_ys=c(min(CSS_sorted[,2]),min(CSS_sorted[,2]),max(CSS_sorted[,2]),max(CSS_sorted[,2]))
beta_ys=mapply(Beta,fmax=fmax,fmin=fmin,fscale=fscale,f=f_ys)
prev_ys=c(rep(mapply(prevpartial_final,w=w_ex,K=polygon_xs[1],f=f_ys[1],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),2),rep(mapply(prevpartial_final,w=w_ex,K=max(K_series),f=f_ys[3],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),2))
H_ys=c(rep(mapply(Hpartial_final,w=w_ex,K=polygon_xs[1],f=f_ys[1],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),2),rep(mapply(Hpartial_final,w=w_ex,K=max(K_series),f=f_ys[3],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),2))
S_ys=c(rep(mapply(Spartial_final,w=w_ex,K=polygon_xs[1],f=f_ys[1],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),2),rep(mapply(Spartial_final,w=w_ex,K=max(K_series),f=f_ys[3],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),2))

#Confirms we are working at the stable equilibrium
range(H_cast*prev_cast-mapply(Icast,w=w_ex,K=CSS_cast[,1],f=CSS_cast[,2],e=e_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex))


#Makes Figure 1 of the manuscript
pdf("AmNatFig1.pdf",width=5,height=7.5,family="ArialMT",useDingbats=FALSE)
{
  cex_smallest_text=.5
  cex_minor_text=.7
  cex_major_text=.8
  lwd_minor=1
  lwd_major=1.5
  dens_choice=10
  par(mfrow=c(4,2),mar=c(.25,.25,.25,.25),mgp=c(1.8,.4,0),cex.lab=cex_major_text,oma=c(3,3,1.5,0),cex.axis=cex_minor_text,bg="white",cex=1)
  plot(CSS_cast[,1],0*CSS_cast[,1],type="l",ylim=1e6*c(min(Beta(fmax,fmin,fscale,example_bottom*.9),na.rm=T),1.1*Beta(fmax,fmin,fscale,example_top)),xlim=range(K_series),col="white",xaxt="n",ylab="")
  abline(h=1e6*Beta(fmax,fmin,fscale,example_bottom),col="gray40",lwd=lwd_minor,lty=3)
  abline(h=1e6*Beta(fmax,fmin,fscale,example_middle),col="gray40",lwd=lwd_minor,lty=2)
  abline(h=1e6*Beta(fmax,fmin,fscale,example_top),col="gray40",lwd=lwd_minor)
  points(CSS_cast[,1],1e6*Beta(fmax,fmin,fscale,CSS_cast[,2]),type="l",lwd=lwd_major)
  mtext(side=3,expression("Castration:"~italic(theta)~"= 0"),padj=-.5,cex=cex_major_text)
  mtext(side=2,cex=cex_major_text,expression(atop("Transmission rate:",italic(beta)~" (L parasite"^"-1"~"day"^"-1"~"x 10"^"-6"~")")),padj=-.5)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_minor_text)
  text(115,38,expression("Less resistant"),cex=cex_minor_text,col="darkgray")
  text(115,4,expression("More resistant"),cex=cex_minor_text,col="darkgray")
  
  plot(CSS[,1],CSS[,1]*0,col="white",xaxt="n",yaxt="n",ylim=c(-1,1e6*1.1*Beta(fmax,fmin,fscale,example_top)),xlim=range(K_series),ylab="")
  polygon(polygon_xs,beta_ys*1e6,col=rgb(0,0,1,0.25,maxColorValue=1),fillOddEven="non-zero",lwd=lwd_minor,border=rgb(0,0,1,.5,maxColorValue=1))
  abline(h=1e6*Beta(fmax,fmin,fscale,example_bottom),col="gray40",lwd=lwd_minor,lty=3)
  abline(h=1e6*Beta(fmax,fmin,fscale,example_middle),col="gray40",lwd=lwd_minor,lty=2)
  abline(h=1e6*Beta(fmax,fmin,fscale,example_top),col="gray40",lwd=lwd_minor)
  mtext(side=3,expression("Fecundity reduction:"~italic(theta)~"= 0.65"),padj=-.5,cex=cex_major_text)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"E",cex=cex_minor_text)
  points(CSS[,1],1e6*Beta(fmax,fmin,fscale,CSS_sorted[,2]),type="l",lwd=lwd_major)
  points(CSS[,1],1e6*Beta(fmax,fmin,fscale,CSS_sorted[,3]),type="l",lwd=lwd_major)
  points(CSS[,1],1e6*Beta(fmax,fmin,fscale,CSS_repellor),type="l",lwd=lwd_major,lty=2)
  text(105,10,"Resistance is futile",cex=cex_minor_text)
  points(CSS[127,1],1e6*Beta(fmax,fmin,fscale,CSS[127,3]),pch=15,cex=cex_major_text)
  points(CSS[127,1],1e6*Beta(fmax,fmin,fscale,CSS[127,2]),pch=24,cex=cex_major_text,bg="white",lwd=lwd_major)
  
  plot(CSS_cast[,1],prev_cast,ylab="",type="l",ylim=c(0,1),lwd=lwd_major,xaxt="n")
  points(CSS_cast[,1],bottom_prev_cast,ylab="",type="l",col="gray40",lwd=lwd_minor,lty=3)
  points(CSS_cast[,1],middle_prev_cast,ylab="",type="l",col="gray40",lwd=lwd_minor,lty=2)
  points(CSS_cast[,1],top_prev_cast,ylab="",type="l",col="gray40",lwd=lwd_minor,lty=1)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_minor_text)
  mtext(side=2,expression(atop("Prevalence:",italic(p)^"*"~" (unitless)")),cex=cex_major_text,padj=-.5)
  text(75,.9,expression("Less resistant"),cex=cex_minor_text,col="darkgray")
  text(75,.06,expression("More resistant"),cex=cex_minor_text,col="darkgray")
  
  plot(CSS[,1],mapply(prevpartial_final,w=w_ex,K=CSS[,1],f=CSS_sorted[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",lwd=lwd_major,xaxt="n",yaxt="n",ylim=c(0,1))
  polygon(polygon_xs,c(mapply(prevpartial_final,w=w_ex,K=CSS[,1],f=CSS[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)[which.min(CSS_sorted[,2])],mapply(prevpartial_final,w=w_ex,K=CSS[,1],f=CSS[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex)[which.min(CSS_sorted[,2])],1,1),col=rgb(0,0,1,0.25,maxColorValue=1),fillOddEven="non-zero",lwd=lwd_minor,border=rgb(0,0,1,.5,maxColorValue=1))
  points(CSS[,1],mapply(prevpartial_final,w=w_ex,K=CSS[,1],f=example_bottom,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor,lty=3)
  points(CSS[,1],mapply(prevpartial_final,w=w_ex,K=CSS[,1],f=example_middle,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor,lty=2)
  points(CSS[,1],mapply(prevpartial_final,w=w_ex,K=CSS[,1],f=example_top,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"F",cex=cex_minor_text)
  points(CSS_bottom[,1],mapply(prevpartial_final,w=w_ex,K=CSS_bottom[,1],f=CSS_bottom[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",ylim=c(0,1),lwd=lwd_major,xaxt="n",yaxt="n")
  points(K_repellor,prev_repellor,type="l",lwd=lwd_major,lty=2)
  points(CSS[127,1],prevpartial_final(w=w_ex,K=CSS[127,1],f=CSS[127,3],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),pch=15,cex=cex_major_text)
  points(CSS[127,1],prevpartial_final(w=w_ex,K=CSS[127,1],f=CSS[127,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),pch=24,cex=cex_major_text,bg="white",lwd=lwd_major)
  
    
  plot(CSS_cast[,1],H_cast,ylab="",type="l",ylim=c(-2,max(H_cast)),lwd=lwd_major,xaxt="n")
  points(CSS_cast[,1],bottom_H_cast,ylab="",type="l",col="gray40",lwd=lwd_minor,lty=3)
  points(CSS_cast[,1],middle_H_cast,ylab="",type="l",col="gray40",lwd=lwd_minor,lty=2)
  points(CSS_cast[,1],top_H_cast,ylab="",type="l",col="gray40",lwd=lwd_minor,lty=1)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_minor_text)
  mtext(side=2,expression(atop("Host density:",italic(H)^"*"~"(Hosts L"^"-1"~")")),cex=cex_major_text,padj=-.5)
  text(75,85,expression("More resistant"),cex=cex_minor_text,col="darkgray")
  text(120,-.5,expression("Less resistant"),cex=cex_minor_text,col="darkgray")
  
  plot(CSS[,1],mapply(Hpartial_final,w=w_ex,K=CSS[,1],f=CSS_sorted[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",lwd=lwd_major,xaxt="n",yaxt="n",ylim=c(-2,max(H_cast)))
  polygon(polygon_xs,H_ys,col=rgb(0,0,1,0.25,maxColorValue=1),fillOddEven="non-zero",lwd=lwd_minor,border=rgb(0,0,1,.5,maxColorValue=1))
  points(CSS[,1],mapply(Hpartial_final,w=w_ex,K=CSS[,1],f=example_bottom,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor,lty=3)
  points(CSS[,1],mapply(Hpartial_final,w=w_ex,K=CSS[,1],f=example_middle,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor,lty=2)
  points(CSS[,1],mapply(Hpartial_final,w=w_ex,K=CSS[,1],f=example_top,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"G",cex=cex_minor_text)
  points(CSS[,1],mapply(Hpartial_final,w=w_ex,K=CSS[,1],f=CSS_sorted[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",ylim=c(0,1),lwd=lwd_major,xaxt="n",yaxt="n")
  points(CSS_bottom[,1],mapply(Hpartial_final,w=w_ex,K=CSS_bottom[,1],f=CSS_bottom[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",ylim=c(0,1),lwd=lwd_major,xaxt="n",yaxt="n")
  points(K_repellor,H_repellor,type="l",lwd=lwd_major,lty=2)
  points(CSS[127,1],Hpartial_final(w=w_ex,K=CSS[127,1],f=CSS[127,3],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),pch=15,cex=cex_major_text)
  points(CSS[127,1],Hpartial_final(w=w_ex,K=CSS[127,1],f=CSS[127,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),pch=24,cex=cex_major_text,bg="white",lwd=lwd_major)
    
  plot(CSS_cast[,1],S_cast,ylab="",type="l",lwd=lwd_major,ylim=c(-2,max(H_cast)))
  points(CSS_cast[,1],bottom_S_cast,ylab="",type="l",col="gray40",lwd=lwd_minor,lty=3)
  #Will need simulations to fill these out
  points(CSS_cast[,1],middle_S_cast,ylab="",type="l",col="gray40",lwd=lwd_minor,lty=2)
  points(CSS_cast[,1],top_S_cast,ylab="",type="l",col="gray40",lwd=lwd_minor)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"D",cex=cex_minor_text)
  mtext(side=2,expression(atop("Susceptible host density:",italic(S)^"*"~"(Hosts L"^"-1"~")")),cex=cex_major_text,padj=-.5)
  text(120,11,expression("Less resistant"),cex=cex_minor_text,col="darkgray")
  text(75,85,expression("More resistant"),cex=cex_minor_text,col="darkgray")
  mtext(side=1,cex=cex_major_text,expression("Carrying capacity: "~italic(K)~"("~mu~"g chl"~italic(a)~"L"^"-1"~")"),padj=1.5,adj=-10)
  
  
  plot(CSS[,1],mapply(Spartial_final,w=w_ex,K=CSS[,1],f=CSS_sorted[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",lwd=lwd_major,yaxt="n",ylim=c(-2,max(H_cast)))
  polygon(polygon_xs,S_ys,col=rgb(0,0,1,0.25,maxColorValue=1),fillOddEven="non-zero",lwd=lwd_minor,border=rgb(0,0,1,.5,maxColorValue=1))
  points(CSS[,1],mapply(Spartial_final,w=w_ex,K=CSS[,1],f=example_bottom,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor,lty=3)
  points(CSS[,1],mapply(Spartial_final,w=w_ex,K=CSS[,1],f=example_middle,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor,lty=2)
  points(CSS[,1],mapply(Spartial_final,w=w_ex,K=CSS[,1],f=example_top,e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",col="gray40",lwd=lwd_minor)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"H",cex=cex_minor_text)
  points(CSS[,1],mapply(Spartial_final,w=w_ex,K=CSS[,1],f=CSS_sorted[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",ylim=c(0,1),lwd=lwd_major,xaxt="n",yaxt="n")
  points(CSS_bottom[,1],mapply(Spartial_final,w=w_ex,K=CSS_bottom[,1],f=CSS_bottom[,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),ylab="",type="l",ylim=c(0,1),lwd=lwd_major,xaxt="n",yaxt="n")
  points(K_repellor,S_repellor,type="l",lwd=lwd_major,lty=2)
  points(CSS[127,1],Spartial_final(w=w_ex,K=CSS[127,1],f=CSS[127,3],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),pch=15,cex=cex_major_text)
  points(CSS[127,1],Spartial_final(w=w_ex,K=CSS[127,1],f=CSS[127,2],e=e_ex,th=th_ex,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_ex,s=s_ex,m=m_ex),pch=24,cex=cex_major_text,bg="white",lwd=lwd_major)
}
dev.off()

#Begin fitting model distributions to field distributions to produce Figure 4.
#Initializing loop
run_length=2e3

#This determine how many "genotypes" are present in each generation
#of the genetic algorithm
simul_trials=50
#This determines the number of generations for which the genetic
#algorithm runs
num_gen=100

#Read in field data. This requires running the 20221031_FigureA7_AmNat.R file first
#setwd("C:/Users/Jason/OneDrive - University of Pittsburgh/Documents/Thesis Chapter 2/Final code files")
field_prev_max_0=readRDS("field_prev_max_0.RDS")

#Can vary this prevalence threshold to fill out the rest of
#Table A1
#Use these different set seeds for the different scenarios
#For max prevalence and 0.1 threshold
set.seed(1234)
#For max prevalence and 0.01 threshold
#set.seed(2341)
#For mean prevalence and 0.1 threshold
#set.seed(3412)
#For mean prevalence and 0.01 threshold
#set.seed(4123)

#Can use this to set the threshold to 0.01 or 0.1.
prev_thresh=0.1
final_prev<-field_prev_max_0[field_prev_max_0>prev_thresh]

#Can also read in mean data (instead of max) to fill out the rest of Table A1
#These .RDS files can be recreated from raw data in FigureA5_AmNat.R

#Read in either of these.
#final_prev<-readRDS("field_prev_mean_0.1.RDS")
#final_prev<-readRDS("field_prev_mean_0.01.RDS")


#Find these descriptive statistics for the field data
#Use the Mclust package to find the bimodality of the field data
model_master_temp1=Mclust(final_prev,G=1)
model_master_temp2=Mclust(final_prev,G=2)
model_master_temp3=Mclust(final_prev,G=32)
if(length(model_master_temp2$parameters$variance$sigmas)==1){
  temp_denom_master=sqrt(model_master_temp2$parameters$variance$sigmasq+model_master_temp2$parameters$variance$sigmasq)
}
if(length(model_master_temp2$parameters$variance$sigmas)==2){
  temp_denom_master=sqrt(model_master_temp2$parameters$variance$sigmasq[1]+model_master_temp2$parameters$variance$sigmasq[2])
}
master_ashman=(2^(1/2)*abs(model_master_temp2$parameters$mean[1]-model_master_temp2$parameters$mean[2]))/temp_denom_master
master_Skew=Skew(final_prev)
master_IQR=IQR(final_prev)

#Can use this to show better BIC with bimodal for both cases
plot(c(summary(model_master_temp1)$bic,summary(model_master_temp2)$bic,summary(model_master_temp3)$bic))

hist(final_prev,breaks=20,xlim=c(0,1),freq=F)
temp1=density(rnorm(1e5,mean=model_master_temp2$parameters$mean[1],sd=sqrt(model_master_temp2$parameters$variance$sigmasq)))
temp2=density(rnorm(1e5,mean=model_master_temp2$parameters$mean[2],sd=sqrt(model_master_temp2$parameters$variance$sigmasq)))
points(temp1$x,temp1$y*model_master_temp2$parameters$pro[1],col="blue",type="l")
points(temp2$x,temp2$y*model_master_temp2$parameters$pro[2],col="red",type="l")

#Evaluate each distribution on its ability to fit master_ashman, master_Skew, and master_IQR of prevalence

#Store the progress of the genetic algorithm.
#Each row represents a generation. Column 1 is the mean K of best
#fit. Column 2 is the sd K of the best fit. Mean and sd are on 
#the log scale. Column 3 is the summed deviance of the best fit
#for that generation.

#In general, a 1 suffix will denote the constrained case and 2
#the 'RiF' case.
#GA_progress1 stores results for the constrained case.
GA_progress1=array(NA,dim=c(num_gen,3))
#GA_progress2 stores results for the 'RiF' case.
GA_progress2=array(NA,dim=c(num_gen,3))

#Construct the trait values for the constrained case
bottom_index=which.min(CSS[,2])

#Construct the constrained case.
constrained=CSS_cast[,2]

#Only allow epidemics in this range, which is the only feasible range anyway
min_thresh=prev_thresh
max_thresh=1

#Mutation size somewhat shapes how the genetic algorithm progresses
mutation_size=0.002

#Starting values for the K distributions.
mean1_seed=30
mean2_seed=30
sd1_seed=10
sd2_seed=10


#Begin loop
#generation is just a counter that doesn't show up anywhere as only final results are preserved
for(generation in 1:num_gen){

#Mutate mean1_seed (either start value or best value from previous generation)
  #depending on whether or not this is generation 1. Creates a vector of mean values
  #that each represent a "genotype" in the genetic algorithm for the constrained case
  #this generation. sd values are the other half of each strategy.
mean1=c(mean1_seed,mean1_seed*10^rnorm(simul_trials-1,mean=0,sd=mutation_size))
sd1=c(sd1_seed,sd1_seed*10^rnorm(simul_trials-1,mean=0,sd=mutation_size))
#Do for the 'RiF' case as well.
mean2=c(mean2_seed,mean2_seed*10^rnorm(simul_trials-1,mean=0,sd=mutation_size))
sd2=c(sd2_seed,sd2_seed*10^rnorm(simul_trials-1,mean=0,sd=mutation_size))

#Create a distribution of K values for each mean and sd
K_values1=array(NA,dim=c(simul_trials,run_length))
K_values2=array(NA,dim=c(simul_trials,run_length))
for(j in 1:simul_trials){
#K_values1[j,]=rlnorm(run_length,meanlog=mean1[j],sdlog=sd1[j])
#K_values2[j,]=rlnorm(run_length,meanlog=mean2[j],sdlog=sd2[j])
  K_values1[j,]=rnorm(run_length,mean=mean1[j],sd=sd1[j])
  K_values2[j,]=rnorm(run_length,mean=mean2[j],sd=sd2[j])
}


#Create some empty arrays to hold upcoming outputs. This also
#functions to clear these arrays from the last generation of the genetic
#algorithm.
K_indices1=array(NA,dim=c(simul_trials,run_length))
K_indices2=array(NA,dim=c(simul_trials,run_length))
prev_evo_constrained=array(NA,dim=c(simul_trials,run_length))
prev_evo_RiF=array(NA,dim=c(simul_trials,run_length,2))
for(j in 1:simul_trials){
for(i in 1:run_length){

  #Find the index of the already-calculated outcome of host evolution
  #at the K value closest to each K value generated in the distribution  
  K_indices1[j,i]=which.min(abs(CSS[,1]-K_values1[j,i]))
  K_indices2[j,i]=which.min(abs(CSS[,1]-K_values2[j,i]))

  #Return the eco-evolutionary prevalence value at each of these K_indices
  
  


  #Only write here if the K_value in question is legitimate and the prev value is feasible
  #Otherwise, write NA. Most likely issue here could be negative prevalence for disease
  #being unfeasible. That would mean K lower than min(K_series) and get turned into a 0 here.
    holder1=prevpartial2(w_ex,CSS[K_indices1[j,i],1],constrained[K_indices1[j,i]],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,constrained[K_indices1[j,i]]),v_ex,s_ex,m_ex)
    prev_evo_constrained[j,i]=ifelse(holder1>min_thresh&holder1<max_thresh,holder1,NA)

    if(!is.na(CSS[K_indices2[j,i],3])){
      holder2a=prevpartial2(w_ex,CSS[K_indices2[j,i],1],CSS[K_indices2[j,i],2],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,CSS[K_indices2[j,i],2]),v_ex,s_ex,m_ex)
      holder2b=prevpartial2(w_ex,CSS[K_indices2[j,i],1],CSS[K_indices2[j,i],3],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,CSS[K_indices2[j,i],3]),v_ex,s_ex,m_ex)
      prev_evo_RiF[j,i,1]=ifelse(holder2a>min_thresh&holder2a<max_thresh,holder2a,NA)
      prev_evo_RiF[j,i,2]=ifelse(holder2b>min_thresh&holder2a<max_thresh,holder2b,NA)      
    }else{
      holder2a=prevpartial2(w_ex,CSS[K_indices2[j,i],1],CSS[K_indices2[j,i],2],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,CSS[K_indices2[j,i],2]),v_ex,s_ex,m_ex)
      prev_evo_RiF[j,i,1]=ifelse(holder2a>min_thresh&holder2a<max_thresh,holder2a,NA)
      prev_evo_RiF[j,i,2]=ifelse(holder2a>min_thresh&holder2a<max_thresh,holder2a,NA)      
    }
  }
}

#List of simul_trials lists. In each list, the first element is the logLik of bimodal fit, the second element is the df of bimodal fit
#the third element is logLik of unimodal fit, the fourth element is df of unimodal fit, the fifth
#element is the probability this distribution would arise from a unimodal fit?
#the sixth element is abs(skew data - master skew)/master skew
#the seventh element is abs(CV data-CV master)/CV master 
eval_full1=list(NA)
eval_full2=list(NA)

for(j in 1:simul_trials){
eval_full1[[j]]=list(NA,NA,NA,NA)
eval_full2[[j]]=list(NA,NA,NA,NA)
}

Ash_D1=array(NA,dim=simul_trials)
Ash_D2=array(NA,dim=simul_trials)
Skew_1=array(NA,dim=simul_trials)
Skew_2=array(NA,dim=simul_trials)

min_num=run_length*.8
for(j in 1:simul_trials){
  if(length(which(!is.na(prev_evo_constrained[j,])))>min_num){
    model1_1=Mclust(prev_evo_constrained[j,which(!is.na(prev_evo_constrained[j,]))],G=1,verbose=F)
    model1_2=Mclust(prev_evo_constrained[j,which(!is.na(prev_evo_constrained[j,]))],G=2,verbose=F)
    eval_full1[[j]][1]=1-pchisq(logLik(model1_2)-logLik(model1_1),df=model1_2$df-model1_1$df)
    if(length(model1_2$parameters$variance$sigmas)==1){
      temp_denom1=sqrt(model1_2$parameters$variance$sigmasq+model1_2$parameters$variance$sigmasq)
    }
    if(length(model1_2$parameters$variance$sigmas)==2){
      temp_denom1=sqrt(model1_2$parameters$variance$sigmasq[1]+model1_2$parameters$variance$sigmasq[2])
    }
    Ash_D1[j]=(2^(1/2)*abs(model1_2$parameters$mean[1]-model1_2$parameters$mean[2]))/temp_denom1
    eval_full1[[j]][2]=abs(Ash_D1[j]-master_ashman)/master_ashman
    Skew_1[j]=Skew(prev_evo_constrained[j,which(!is.na(prev_evo_constrained[j,]))])
    eval_full1[[j]][3]=abs(Skew_1[j]-master_Skew)/master_Skew
    eval_full1[[j]][4]=abs(IQR(prev_evo_constrained[j,which(!is.na(prev_evo_constrained[j,]))])-master_IQR)/master_IQR
  }else{
    Ash_D1[j]=NA
    eval_full1[[j]][2]=NA
    Skew_1[j]=NA
    eval_full1[[j]][3]=NA
    eval_full1[[j]][4]=NA
  }

  #Include the lower CSS
  if(length(which(!is.na(prev_evo_RiF[j,,1])))>min_num){
    model2_1=Mclust(c(prev_evo_RiF[j,which(!is.na(prev_evo_RiF[j,,1])),]),G=1,verbose=F)
    model2_2=Mclust(c(prev_evo_RiF[j,which(!is.na(prev_evo_RiF[j,,1])),]),G=2,verbose=F)
    eval_full2[[j]][1]=1-pchisq(logLik(model2_2)-logLik(model2_1),df=model2_2$df-model2_1$df)
    if(length(model2_2$parameters$variance$sigmas)==1){
      temp_denom2=sqrt(model2_2$parameters$variance$sigmasq+model2_2$parameters$variance$sigmasq)
    }
    if(length(model2_2$parameters$variance$sigmas)==2){
      temp_denom2=sqrt(model2_2$parameters$variance$sigmasq[1]+model2_2$parameters$variance$sigmasq[2])
    }
    Ash_D2[j]=(2^(1/2)*abs(model2_2$parameters$mean[1]-model2_2$parameters$mean[2]))/temp_denom2
    eval_full2[[j]][2]=abs(Ash_D2[j]-master_ashman)/master_ashman
    Skew_2[j]=Skew(c(prev_evo_RiF[j,which(!is.na(prev_evo_RiF[j,,1])),]))
    eval_full2[[j]][3]=abs(Skew_2[j]-master_Skew)/master_Skew
    eval_full2[[j]][4]=abs(IQR(c(prev_evo_RiF[j,which(!is.na(prev_evo_RiF[j,,1])),]))-master_IQR)/master_IQR
  }else{
    Ash_D2[j]=NA
    eval_full2[[j]][2]=NA
    Skew_2[j]=NA
    eval_full2[[j]][3]=NA
    eval_full2[[j]][4]=NA
  }
}

#Score the fits, here bimod is the score, not the parameter, etc
bimod1=array(NA,dim=simul_trials)
skew1=array(NA,dim=simul_trials)
IQR1=array(NA,dim=simul_trials)
bimod2=array(NA,dim=simul_trials)
skew2=array(NA,dim=simul_trials)
IQR2=array(NA,dim=simul_trials)

for(j in 2:simul_trials){
  bimod1[j]=as.numeric(eval_full1[[j]][2])
  skew1[j]=as.numeric(eval_full1[[j]][3])
  IQR1[j]=as.numeric(eval_full1[[j]][4])
  
  bimod2[j]=as.numeric(eval_full2[[j]][2])
  skew2[j]=as.numeric(eval_full2[[j]][3])
  IQR2[j]=as.numeric(eval_full2[[j]][4])
}

total1=bimod1+skew1+IQR1
total2=bimod2+skew2+IQR2

#total1=bimod1+IQR1
#total2=bimod2+IQR2

best_index1=which.min(total1)
best_index2=which.min(total2)

#Store progress of the algorithm but only if it HAS made progress
if(generation==1||total1[best_index1]<GA_progress1[generation-1,3]){
  mean1_seed<-mean1[best_index1]
  sd1_seed<-sd1[best_index1]
  GA_progress1[generation,]<-c(mean1[best_index1],sd1[best_index1],total1[best_index1])  
}else{
  GA_progress1[generation,]<-GA_progress1[generation-1,]
}
if(generation==1||total2[best_index2]<GA_progress2[generation-1,3]){
  mean2_seed<-mean2[best_index2]
  sd2_seed<-sd2[best_index2]  
  GA_progress2[generation,]<-c(mean2[best_index2],sd2[best_index2],total2[best_index2])
}else{
  GA_progress2[generation,]<-GA_progress2[generation-1,]
}

#Plot current best model distributions.
par(mfrow=c(1,2),mar=c(4,4,.5,.5),cex.axis=.8,cex.lab=1)
if(length(best_index1)==1){
hist(prev_evo_constrained[best_index1,],freq=F,main="",xlab="Evo constrained prev",breaks=20)
text((par("usr")[2]-par("usr")[1])*.10+par("usr")[1],(par("usr")[4]-par("usr")[3])*.15+par("usr")[3],generation,cex=1)
text((par("usr")[2]-par("usr")[1])*.80+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],GA_progress1[generation,3],cex=1)
text((par("usr")[2]-par("usr")[1])*.80+par("usr")[1],(par("usr")[4]-par("usr")[3])*.90+par("usr")[3],GA_progress1[generation,1],cex=1)
text((par("usr")[2]-par("usr")[1])*.80+par("usr")[1],(par("usr")[4]-par("usr")[3])*.85+par("usr")[3],GA_progress1[generation,2],cex=1)
}

if(length(best_index2)==1){
#From using both CSSs
hist(c(prev_evo_RiF[best_index2,!is.na(prev_evo_RiF[best_index2,,1]),]),freq=F,main="",xlab="Evo RiF prev",breaks=20)
text((par("usr")[2]-par("usr")[1])*.10+par("usr")[1],(par("usr")[4]-par("usr")[3])*.15+par("usr")[3],generation,cex=1)
text((par("usr")[2]-par("usr")[1])*.80+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],GA_progress2[generation,3],cex=1)
text((par("usr")[2]-par("usr")[1])*.80+par("usr")[1],(par("usr")[4]-par("usr")[3])*.90+par("usr")[3],GA_progress2[generation,1],cex=1)
text((par("usr")[2]-par("usr")[1])*.80+par("usr")[1],(par("usr")[4]-par("usr")[3])*.85+par("usr")[3],GA_progress2[generation,2],cex=1)
}

print(generation)
}

#Progress of genetic algorithm
par(mfrow=c(1,1))
plot(seq(1,100,1),GA_progress1[,3],type="l",ylim=c(0,6),ylab="Difference from observed distribution",xlab="Generation of genetic algorithm")
points(seq(1,100,1),GA_progress2[,3],type="l",col="blue")
text(50,1.6,"No RiF")
text(50,.4,"RiF",col="blue")

#2 is ashman's D, 3 is Skew, 4 is IQR
IQR(prev_evo_RiF[best_index2,!is.na(prev_evo_RiF[best_index2,,1]),],na.rm=T)
IQR(prev_evo_constrained[best_index1,],na.rm=T)
Ash_D2[best_index2]
Ash_D1[best_index1]
Skew_2[best_index2]
Skew_1[best_index1]


master_IQR
master_ashman
master_Skew

abs(Ash_D2[best_index2]-master_ashman)/master_ashman+abs(Skew_2[best_index2]-master_Skew)/master_Skew+abs(IQR(prev_evo_RiF[best_index2,!is.na(prev_evo_RiF[best_index2,,1]),],na.rm=T)-master_IQR)/master_IQR
abs(Ash_D1[best_index1]-master_ashman)/master_ashman+abs(Skew_1[best_index1]-master_Skew)/master_Skew+abs(IQR(prev_evo_constrained[best_index1,],na.rm=T)-master_IQR)/master_IQR

#For prev_thresh = 0.1 and max field prev
#mean1_seed=28.25593
#sd1_seed=9.994347
#mean2_seed=30.06419
#sd2_seed=9.938546

#Putting K values and prevalence values in convenient objects for plotting

#Generate the K values according to the best fit distribution for each case
Ks_1=rnorm(run_length,mean=mean1_seed,sd=sd1_seed)
Ks_2=rnorm(run_length,mean=mean2_seed,sd=sd2_seed)


#Create arrays to store loop outputs
K_index1=array(NA,dim=run_length)
K_index2=array(NA,dim=run_length)
prev_1=array(NA,dim=run_length)
prev_2=array(NA,dim=c(run_length,2))
for(i in 1:run_length){
  #Get the index of the closest K value for which evolution was pre-calculated
  K_index1[i]=which.min(abs(CSS[,1]-Ks_1[i]))
  K_index2[i]=which.min(abs(CSS[,1]-Ks_2[i]))
  #Calculate prevalence
  holder1=prevpartial2(w_ex,CSS[K_index1[i],1],constrained[K_index1[i]],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,constrained[K_index1[i]]),v_ex,s_ex,m_ex)

  #Only write here if the K_value in question is legitimate

    prev_1[i]=ifelse(holder1>min_thresh&holder1<max_thresh,holder1,NA)    
  
    if(!is.na(CSS[K_index2[i],3])){
      holder2a=prevpartial2(w_ex,CSS[K_index2[i],1],CSS[K_index2[i],2],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,CSS[K_index2[i],2]),v_ex,s_ex,m_ex)
      holder2b=prevpartial2(w_ex,CSS[K_index2[i],1],CSS[K_index2[i],3],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,CSS[K_index2[i],3]),v_ex,s_ex,m_ex)
      prev_2[i,1]=ifelse(holder2a>min_thresh&holder2a<max_thresh,holder2a,NA)
      prev_2[i,2]=ifelse(holder2b>min_thresh&holder2a<max_thresh,holder2b,NA)      
    }else{
      holder2a=prevpartial2(w_ex,CSS[K_index2[i],1],CSS[K_index2[i],2],e_ex,th_ex,d_ex,Beta(fmax,fmin,fscale,CSS[K_index2[i],2]),v_ex,s_ex,m_ex)
      prev_2[i,1]=ifelse(holder2a>min_thresh&holder2a<max_thresh,holder2a,NA)
      prev_2[i,2]=ifelse(holder2a>min_thresh&holder2a<max_thresh,holder2a,NA)      
    }
}

order1=order(Ks_1)
order2=order(Ks_2)
Ks_1_ordered=Ks_1[order1]
prev_1_ordered=prev_1[order1]
Ks_2_ordered=Ks_2[order2]
prev_2_ordered=prev_2[order2]
plot(Ks_1_ordered,prev_1_ordered,type="l")
plot(Ks_2_ordered,prev_2_ordered,type="l")

hist(prev_1,xlim=c(0,1))
hist(prev_2,xlim=c(0,1))







#Make Figure 4.
#pdf("AmNat_Fig4.pdf",width=14,height=9,family="ArialMT",useDingbats=FALSE)
  {
  cex_smallest_text=.5
  cex_minor_text=1
  cex_major_text=1
  lwd_minor=1
  lwd_major=2.5
  par(mfrow=c(2,2),mar=c(4,5,.1,.5),mgp=c(1.8,.4,0),cex.lab=cex_major_text,oma=c(0,0,0,0),cex.axis=cex_minor_text,bg="white",cex=1)
  top=1
  num_breaks=20
  bottom=0
  low_y_space=.1
  mid_y_space=.05
  top_y_space=.1
  left_space=.1
  mid_x_space=.1
  right_space=.1
  left=0
  right=1
  x_wid=(right-left-left_space-mid_x_space-right_space)/2
  y_wid=(top-bottom-low_y_space-mid_y_space-top_y_space)/2

  hist(final_prev,freq=F,breaks=num_breaks,xlim=c(0,max(prev_2,na.rm=T)),main="",xlab=expression("Observed maximum"~italic(p)),ylab="Probability density")
  box(which="plot",lty=1)
  points(density(final_prev),type="l",lwd=lwd_major)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_minor_text)

  plot(CSS[,1],1e6*Beta(fmax,fmin,fscale,CSS[,2]),type="l",col=rgb(0,0,1,0.5,maxColorValue=1),xlim=c(min(CSS[,1]),max(c(Ks_1,Ks_2))),
       ylab=expression(atop("Transmission rate:",italic(beta)~" (L parasite"^"-1"~"day"^"-1"~"x 10"^"-6"~")")),
       xlab=expression(paste("Carrying capacity: ",italic(K)," (",mu,"g chl ",italic(a)," L"^"-1",")",sep="")),lwd=lwd_major,ylim=c(0,max(1e6*Beta(fmax,fmin,fscale,CSS[,2]))))
  bh2=30
  bh1=28
  points(CSS[,1],1e6*Beta(fmax,fmin,fscale,CSS_cast[,2]),type="l",col="red",lwd=lwd_major)
  points(mean2_seed,bh2,col=rgb(0,0,1,0.5,maxColorValue=1),pch=16)
  segments(x0=mean2_seed+sd2_seed,x1=mean2_seed-sd2_seed,y0=bh2,y1=bh2,col=rgb(0,0,1,0.5,maxColorValue=1))
  points(mean1_seed,bh1,col="red",pch=16)
  segments(x0=mean1_seed+sd1_seed,x1=mean1_seed-sd1_seed,y0=bh1,y1=bh1,col="red")
  text(54,16,cex=cex_minor_text,col=rgb(0,0,1,.5,maxColorValue=1),expression(atop("Resistance","is futile")))
  text(54,4.5,cex=cex_minor_text,col="red","Constrained")
  text(28,32.5,expression(paste("Mean","\\u00B1","sd",sep="")))
  text((par("usr")[2]-par("usr")[1])*.075+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_minor_text)

  
hist(prev_1,freq=F,main="",xlab=expression("Evolved prevalence:"~italic(p)^"*"~""["CSS"]),breaks=num_breaks,ylab="Probability density",xlim=c(0,max(prev_2,na.rm=T)))
box(which="plot",lty=1)
points(density(c(prev_1)[which(!is.na(c(prev_1)))]),type="l",lwd=lwd_major,col="red")
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_minor_text)
text((par("usr")[2]-par("usr")[1])*.8+par("usr")[1],(par("usr")[4]-par("usr")[3])*.9+par("usr")[3],cex=cex_minor_text,col="red","Constrained")
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)

hist(prev_2,freq=F,main="",xlab=expression("Evolved prevalence:"~italic(p)^"*"~""["CSS"]),breaks=num_breaks,xlim=c(0,max(prev_2,na.rm=T)),ylab="Probability density")
text((par("usr")[2]-par("usr")[1])*.8+par("usr")[1],(par("usr")[4]-par("usr")[3])*.9+par("usr")[3],cex=cex_minor_text,expression(atop("Resistance","is futile")),col=rgb(0,0,1,0.5,maxColorValue=1))
box(which="plot",lty=1)
points(density(c(prev_2)[which(!is.na(c(prev_2)))]),type="l",lwd=lwd_major,col=rgb(0,0,1,0.5,maxColorValue=1))
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"D",cex=cex_minor_text)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
}
#dev.off()


parm_length=50
theta_range=10^linspace(0,log10(.1),parm_length)
v_range=linspace(v_ex/5,10*v_ex,parm_length)

temp_K_length=180

#Store the K_ranges for each parameter set
theta_K_ranges=array(NA,dim=c(parm_length,temp_K_length))
v_K_ranges=array(NA,dim=c(parm_length,temp_K_length))

theta_K_RiF=array(NA,dim=parm_length)
v_K_RiF=array(NA,dim=parm_length)
theta_f_RiF=array(NA,dim=parm_length)
v_f_RiF=array(NA,dim=parm_length)

unlister=function(input,maxwidth){
  if(typeof(input)=="list"){
    output=array(NA,dim=c(length(input),maxwidth))
    for(i in 1:length(input)){
      tempwidth=NA
      tempwidth=length(input[[i]])
      for(j in 1:tempwidth){
        output[i,j]=sort(input[[i]])[j]
      }
    }
  }else{
    output=input
  }
  return(output)
}

K_RiF_finder=function(matrix){
  rowMaxlength=ifelse(!is.null(dim(matrix)),dim(matrix)[1],length(matrix))
  rowMax=array(NA,dim=rowMaxlength)
  if(length(matrix)>rowMaxlength){
  for(i in 1:rowMaxlength){
    tempmax=NA
    tempmax=max(matrix[i,],na.rm=T)
    if(tempmax<0){
      tempmax<-NA
    }
    rowMax[i]=tempmax
  }
  }else{
    rowMax=matrix
  }
  if(length(which(diff(rowMax)>0))>0){
    output=min(which(diff(rowMax)>0))+1
  }else{
    output=0
  }
  return(output)
}

K_basal=linspace(15,16,20)
basal_test=mapply(NGM_checkedroots,w_func=w_ex,K_func=K_basal,e_func=e_ex,th_func=theta_range[1],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[1],s_func=s_ex,m_func=m_ex)

K_start=K_basal[K_RiF_finder(unlister(basal_test,10))]
max(basal_test[[K_RiF_finder(unlister(basal_test,10))]])

theta_K_ranges[1,]=linspace(15,16,temp_K_length)
v_K_ranges[1,]=linspace(15,16,temp_K_length)

theta_K_RiF[1]=K_start
v_K_RiF[1]=K_start
theta_f_RiF[1]=max(basal_test[[K_RiF_finder(unlister(basal_test,10))]])
v_f_RiF[1]=max(basal_test[[K_RiF_finder(unlister(basal_test,10))]])


#Can run this calculation or read in the results below to save a great deal
#of time
for(i in 2:parm_length){
  tic()
  theta_K_ranges[i,]=10^linspace(log10(theta_K_RiF[i-1]),log10(theta_K_RiF[i-1]*10),temp_K_length)
  v_K_ranges[i,]=10^linspace(log10(v_K_RiF[i-1]),log10(v_K_RiF[i-1]*10),temp_K_length)
  
  theta_CSS_list<-NA
  v_CSS_list<-NA
  theta_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=theta_K_ranges[i,],e_func=e_ex,th_func=theta_range[i],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[1],s_func=s_ex,m_func=m_ex)  
  v_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=v_K_ranges[i,],e_func=e_ex,th_func=theta_range[1],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[i],s_func=s_ex,m_func=m_ex)  

  theta_CSS_matrix=unlister(theta_CSS_list,10)
  v_CSS_matrix=unlister(v_CSS_list,10)
  
  #If RiF was found, repeat in a narrower interval.
  #Else, repeat in a higher interval.
  if(K_RiF_finder(theta_CSS_matrix)>0){
    theta_K_ranges[i,]=linspace(theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)-1],theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)],temp_K_length)
    theta_CSS_list<-NA
    theta_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=theta_K_ranges[i,],e_func=e_ex,th_func=theta_range[i],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[1],s_func=s_ex,m_func=m_ex)  
    theta_CSS_matrix=unlister(theta_CSS_list,10)
    theta_K_RiF[i]=theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)]
    theta_f_RiF[i]=max(theta_CSS_matrix[[K_RiF_finder(theta_CSS_matrix)]])
  }else{
    theta_K_ranges[i,]=10^linspace(log10(theta_K_RiF[i-1])+1,log10(theta_K_RiF[i-1])+2,temp_K_length)
    theta_CSS_list<-NA
    theta_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=theta_K_ranges[i,],e_func=e_ex,th_func=theta_range[i],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[1],s_func=s_ex,m_func=m_ex)  
    theta_CSS_matrix=unlister(theta_CSS_list,10)
    
    #Now to deal with whether or not RiF was found in higher interval
    if(K_RiF_finder(theta_CSS_matrix)>0){
      theta_K_ranges[i,]=linspace(theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)-1],theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)],temp_K_length)
      theta_CSS_list<-NA
      theta_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=theta_K_ranges[i,],e_func=e_ex,th_func=theta_range[i],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[1],s_func=s_ex,m_func=m_ex)  
      theta_CSS_matrix=unlister(theta_CSS_list,10)
      theta_K_RiF[i]=theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)]
      theta_f_RiF[i]=max(theta_CSS_matrix[[K_RiF_finder(theta_CSS_matrix)]])
    }else{
      theta_K_ranges[i,]=10^linspace(log10(theta_K_RiF[i-1])+2,log10(theta_K_RiF[i-1])+3,temp_K_length)
      theta_CSS_list<-NA
      theta_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=theta_K_ranges[i,],e_func=e_ex,th_func=theta_range[i],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[1],s_func=s_ex,m_func=m_ex)  
      theta_CSS_matrix=unlister(theta_CSS_list,10)
      
      #Now to deal with whether or not RiF was found in the much higher interval
      if(K_RiF_finder(theta_CSS_matrix)>0){
        theta_K_ranges[i,]=linspace(theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)-1],theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)],temp_K_length)
        theta_CSS_list<-NA
        theta_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=theta_K_ranges[i,],e_func=e_ex,th_func=theta_range[i],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[1],s_func=s_ex,m_func=m_ex)  
        theta_CSS_matrix=unlister(theta_CSS_list,10)
        theta_K_RiF[i]=theta_K_ranges[i,K_RiF_finder(theta_CSS_matrix)]
        theta_f_RiF[i]=max(theta_CSS_matrix[[K_RiF_finder(theta_CSS_matrix)]])
      }
    }
  }
  
  #Same for v
  if(K_RiF_finder(v_CSS_matrix)>0){
    v_K_ranges[i,]=linspace(v_K_ranges[i,K_RiF_finder(v_CSS_matrix)-1],v_K_ranges[i,K_RiF_finder(v_CSS_matrix)],temp_K_length)
    v_CSS_list<-NA
    v_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=v_K_ranges[i,],e_func=e_ex,th_func=theta_range[1],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[i],s_func=s_ex,m_func=m_ex)  
    v_CSS_matrix=unlister(v_CSS_list,10)
    v_K_RiF[i]=v_K_ranges[i,K_RiF_finder(v_CSS_matrix)]
    v_f_RiF[i]=max(v_CSS_matrix[[K_RiF_finder(v_CSS_matrix)]])
  }else{
    v_K_ranges[i,]=10^linspace(log10(v_K_RiF[i-1])+1,log10(v_K_RiF[i-1])+2,temp_K_length)
    v_CSS_list<-NA
    v_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=v_K_ranges[i,],e_func=e_ex,th_func=theta_range[1],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[i],s_func=s_ex,m_func=m_ex)  
    v_CSS_matrix=unlister(v_CSS_list,10)
    
    #Now to deal with whether or not RiF was found in higher interval
    if(K_RiF_finder(v_CSS_matrix)>0){
      v_K_ranges[i,]=linspace(v_K_ranges[i,K_RiF_finder(v_CSS_matrix)-1],v_K_ranges[i,K_RiF_finder(v_CSS_matrix)],temp_K_length)
      v_CSS_list<-NA
      v_CSS_list=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=v_K_ranges[i,],e_func=e_ex,th_func=theta_range[1],d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_range[i],s_func=s_ex,m_func=m_ex)  
      v_CSS_matrix=unlister(v_CSS_list,10)
      v_K_RiF[i]=v_K_ranges[i,K_RiF_finder(v_CSS_matrix)]
      v_f_RiF[i]=max(v_CSS_matrix[[K_RiF_finder(v_CSS_matrix)]])
    }
  }
  
  par(mfrow=c(2,2),xaxs="r",yaxs="r")
  plot(theta_range[1:i],theta_K_RiF[1:i],ylab="log10(K) of RiF",xlab="Theta")
  plot(v_range[1:i],v_K_RiF[1:i],ylab="log10(K) of RiF",xlab="v")
  plot(theta_range[1:i],mapply(prevpartial_final,w=w_ex,K=theta_K_RiF[1:i],f=theta_f_RiF[1:i],e=e_ex,th=theta_range[1:i],d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_range[1],s=s_ex,m=m_ex),ylab="prev of RiF",xlab="Theta")
  plot(v_range[1:i],mapply(prevpartial_final,w=w_ex,K=v_K_RiF[1:i],f=v_f_RiF[1:i],e=e_ex,th=theta_range[1],d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_range[1:i],s=s_ex,m=m_ex),ylab="prev of RiF",xlab="v")
  print(i)
  toc()
}

#saveRDS(theta_K_RiF,"theta_K_RiF.RDS")
#saveRDS(v_K_RiF,"v_K_RiF.RDS")
#saveRDS(theta_f_RiF,"theta_f_RiF.RDS")
#saveRDS(v_f_RiF,"v_f_RiF.RDS")


theta_K_RiF=readRDS("theta_K_RiF.RDS")
v_K_RiF=readRDS("v_K_RiF.RDS")
theta_f_RiF=readRDS("theta_f_RiF.RDS")
v_f_RiF=readRDS("v_f_RiF.RDS")

#Might also want to plot something about what RiF means for host density
#and see if that differs between theta and v. High theta seems to lure
#hosts toward low density. Does low v do the same?

#Figure A1 of the manuscript
#png("Ch2_RiF_summaries.png",width=8,height=8,units="in",res=300)
#{
  cex_smallest_text=.5*1.2
  cex_minor_text=1*1.2
  cex_major_text=1.5
  cex_big_points=3
  lwd_minor=1.5*1.2
  lwd_major=2*1.2
  dens_choice=10
  yline=2
  xline=5
  par(mfrow=c(2,2),mar=c(1.5,1.5,1,1),cex.lab=cex_major_text,cex.axis=cex_minor_text,cex=1,oma=c(5,5,0,0),font.lab=2,xaxs="i")
plot(1-theta_range,log10(theta_K_RiF),ylab="",xlab="",type="l",lwd=lwd_major,xlim=c(0,0.92),xaxt="n",yaxt="n")
axis(side=2,at=c(2,3),labels=c(expression("10"^"2"),expression("10"^"3")))
mtext(side=2,expression(atop("Carrying capacity","of RiF:"~italic(K)[RiF]~"("~mu~"g chl"~italic(a)~"L"^"-1"~")")),cex=cex_major_text,line=yline)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_major_text)

plot(v_range,v_K_RiF,ylab="",xlab="",type="l",lwd=lwd_major,xaxt="n",yaxt="n")
axis(side=2,at=c(20,60,100,140),labels=c("20","60","100","140"))
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_major_text)

plot(1-theta_range,mapply(prevpartial_final,w=w_ex,K=theta_K_RiF,f=theta_f_RiF,e=e_ex,th=theta_range,d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_range[1],s=s_ex,m=m_ex),ylab="",type="l",lwd=lwd_major,xlim=c(0,0.92),ylim=c(0,0.65))
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_major_text)
mtext(side=2,expression(atop("Prevalence of RiF:",italic(p)[RiF]^"*"~"(unitless)")),cex=cex_major_text,line=yline)
mtext(side=1,expression(atop("Fecundity reduction:","1-"~italic(theta)~"(unitless)")),cex=cex_major_text,line=xline)

plot(v_range,mapply(prevpartial_final,w=w_ex,K=v_K_RiF,f=v_f_RiF,e=e_ex,th=theta_range[1],d=d_ex,fmax=fmax,fmin=fmin,fscale=fscale,v=v_range,s=s_ex,m=m_ex),ylab="",xlab="",type="l",lwd=lwd_major,ylim=c(0,0.65))
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"D",cex=cex_major_text)
mtext(side=1,expression(atop("Mortality virulence:",italic(v)~"(day"^"-1"~")")),cex=cex_major_text,line=xline)
#}
#dev.off()

#Just make sure we are dealing with equilibrium 6 everywhere.

#This might just be an issue of finding the CSS accurately enough
mapply(Which_stable,w=w_ex,K=theta_K_RiF,f=theta_f_RiF,e=e_ex,th=theta_range,d=d_ex,B=Beta(fmax,fmin,fscale,theta_f_RiF),v=v_ex,s=s_ex,m=m_ex)

#Dedicate a little more numerical precision to finding the CSS here since
#the numerical results are least certain there at the end.
extra_accuracy=NGM_checkedroots(w_func=w_ex,K_func=theta_K_RiF[50],e_func=e_ex,th_func=theta_range[50],d_func=d_ex,fmax_func =fmax,fmin_func = fmin,fscale_func = fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)

#Yes so we're working with the expected equilibrium everywhere.
mapply(Which_stable,w=w_ex,K=theta_K_RiF[50],f=extra_accuracy,e=e_ex,th=theta_range[50],d=d_ex,B=Beta(fmax,fmin,fscale,extra_accuracy),v=v_ex,s=s_ex,m=m_ex)

#This looks great
mapply(Which_stable,w=w_ex,K=v_K_RiF,f=v_f_RiF,e=e_ex,th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,v_f_RiF),v=v_range,s=s_ex,m=m_ex)
