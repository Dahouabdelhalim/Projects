#Clear the environment
rm(list=ls())

#Load necessary libraries
library(Deriv)
library(rootSolve)
library(pracma)

#Set the working directory to wherever you downloaded the data files
#setwd("C:/Users/Jason/OneDrive - University of Pittsburgh/Documents/Thesis Chapter 2")

#where B (Beta) = f/fscale*log((fmax-fmin)/(fmax-f)) (eq. 2)
Beta=function(fmax,fmin,fscale,f){
  f/fscale*log((fmax-fmin)/(fmax-f))
}
#The only stable, endemic equilibrium imported from results
#of Symbolic_prep_AmNat

#Take symbolic derivatives and define jacobian functions
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
  
  Jacobian_symbolic=list(c(J11,J12,J13,J14),c(J21,J22,J23,J24),c(J31,J32,J33,J34),c(J41,J42,J43,J44))
  
  Jacobian_numeric=function(w,K,f,e,th,d,B,v,s,m,I,S,R,Z){
    Jacobian_holder=array(NA,dim=c(4,4))
    for(i in 1:4){
      for(j in 1:4){
        Jacobian_holder[i,j]=Jacobian_symbolic[[i]][[j]](w,K,f,e,th,d,B,v,s,m,I,S,R,Z)
      }
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

#Check which, if any, equilibria are stable, including the two disease free boundary equilibria (not trivial)
Which_stable=function(w,K,f,e,th,d,B,v,s,m){
  if(!is.na(w*K*f*e*th*d*B*v*s*m)){
  Equilibria_numeric=Equilibria(w,K,f,e,th,d,B,v,s,m)
  eigen_length=length(Equilibria_numeric)/4
  Max_eigen=array(NA,dim=eigen_length)
  for(i in 1:eigen_length){
    Max_eigen[i]=max(Re(eigen(Jacobian_numeric(w,K,f,e,th,d,B,v,s,m,Equilibria_numeric[1,i],Equilibria_numeric[2,i],Equilibria_numeric[3,i],Equilibria_numeric[4,i]))$values))
  }
  ifelse(length(which(Max_eigen<0))==1,return(which(Max_eigen<0)),return(paste("Error: number of stable equilibria is",length(which(Max_eigen<0)),sep="")))
}else{
  return(NA)
}
}

#Return the stable equilibrium values
Eigen_out=function(w,K,f,e,th,d,B,v,s,m){
  if(!is.na(w*K*f*e*th*d*B*v*s*m)){
  Equilibria_numeric=Equilibria(w,K,f,e,th,d,B,v,s,m)
  eigen_length=length(Equilibria_numeric)/4
  Max_eigen=array(NA,dim=eigen_length)
  for(i in 1:eigen_length){
    Max_eigen[i]=max(Re(eigen(Jacobian_numeric(w,K,f,e,th,d,B,v,s,m,Equilibria_numeric[1,i],Equilibria_numeric[2,i],Equilibria_numeric[3,i],Equilibria_numeric[4,i]))$values))
  }
  ifelse(length(which(Max_eigen<0))==1,return(Equilibria_numeric[,which(Max_eigen<0)]),return(length(which(Max_eigen<0))))
  }else{
   return(c(NA,NA,NA,NA))
 }
  #Outputs as I, S, R, Z
}

#w,K,e,th,d,fmax,fmin,fscale,v,s,m
Standard_parms<-c(.9,100,.25,1,.05,.0039,.0024,9918.944,.037,100000,.5)

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

#Calculate a proxy for fitness according to NGM theory
#using fM (mutant trait) and fR (resident trait)
Fit_proxy=function(fM,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m){
  #Written simply, (efR(d+v+th*B*Z))/((d+v)*(d+B*Z))
  (e*fM*Rpartial(w,K,fR,e,th,d,fmax,fmin,fscale,v,s,m)*(d+v)+e*fM*Rpartial(w,K,fR,e,th,d,fmax,fmin,fscale,v,s,m)*th*fM/fscale*log((fmax-fmin)/(fmax-fM))*Zpartial(w,K,fR,e,th,d,fmax,fmin,fscale,v,s,m))/((d+v)*(d+fM/fscale*log((fmax-fmin)/(fmax-fM))*Zpartial(w,K,fR,e,th,d,fmax,fmin,fscale,v,s,m)))
}

#Subtract 1 such that this is positive or negative to match sign of invader fitness
NGM_root=function(fM,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m){
  Fit_proxy(fM,fR,w,K,e,th,d,fmax,fmin,fscale,v,s,m)-1
}

#Take key derivatives with respect to the invader's fitness proxy
NGM_difffM=Deriv(NGM_root,"fM")
NGM_diff2fM=Deriv(NGM_difffM,"fM")
NGM_difffR=Deriv(NGM_root,"fR")
NGM_diff2fR=Deriv(NGM_difffR,"fR")

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

#Make some large, logarithmically spaced range of non-focal parameter values that
#includes the standard values
run_length=200

fmax=0.02069553
fmin=0.0051666
fscale=2983.647
Standard_parms<-c(0.9,40,.18,.65,.05,fmax,fmin,fscale,.05,100000,1.5)
#For w_ex=0.9,K_ex=30,e_ex=0.18,th_ex=0.65,d_ex=0.05,fmax,fmin,fscale,v_ex=0.045,s_ex=1e5,m_ex=1.5
w_ex=Standard_parms[1]
K_ex=Standard_parms[2]
e_ex=Standard_parms[3]
th_ex=Standard_parms[4]
d_ex=Standard_parms[5]
v_ex=Standard_parms[9]
s_ex=Standard_parms[10]
m_ex=Standard_parms[11]


#Construct logarithmically spaced gradients in these variables.
for(i in 1:run_length){
  w_series=10^linspace(log10(w_ex)-1,log10(w_ex)+1,run_length)
  e_series=10^linspace(log10(e_ex)-1,log10(e_ex)+1,run_length)
  th_series=10^linspace(log10(th_ex)-1,log10(th_ex)+1,run_length)
  d_series=10^linspace(log10(d_ex)-1,log10(d_ex)+1,run_length)
  v_series=10^linspace(log10(v_ex)-1,log10(v_ex)+1,run_length)
  s_series=10^linspace(log10(s_ex)-1,log10(s_ex)+1,run_length)
  m_series=10^linspace(log10(m_ex)-1,log10(m_ex)+1,run_length)
}

#Calculate the CSSs for each parameter range. Right now, this is whole
#section is commented
#out as it takes a long time and you can just read in the results as .RDS files.
{
#evo_w_series=mapply(FUN=NGM_checkedroots,w_func=w_series,K_func=K_ex,e_func=e_ex,th_func=th_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)
#evo_e_series=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=K_ex,e_func=e_series,th_func=th_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)
#evo_d_series=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=K_ex,e_func=e_ex,th_func=th_ex,d_func=d_series,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)
#evo_v_series=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=K_ex,e_func=e_ex,th_func=th_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_series,s_func=s_ex,m_func=m_ex)
#evo_s_series=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=K_ex,e_func=e_ex,th_func=th_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_series,m_func=m_ex)
#evo_m_series=mapply(FUN=NGM_checkedroots,w_func=w_ex,K_func=K_ex,e_func=e_ex,th_func=th_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_series)

#Get nonCSS singular points for these two gradients.
#singular_e_series=mapply(FUN=NGM_singulars,w_func=w_ex,K_func=K_ex,e_func=e_series,th_func=th_ex,d_func=d_ex,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)
#singular_d_series=mapply(FUN=NGM_singulars,w_func=w_ex,K_func=K_ex,e_func=e_ex,th_func=th_ex,d_func=d_series,fmax_func=fmax,fmin_func=fmin,fscale_func=fscale,v_func=v_ex,s_func=s_ex,m_func=m_ex)

#Be prepared to store up to 10 CSSes, though there will certainly be
#fewer than that.
#evolved_w_series=array(NA,dim=c(run_length,10))
#evolved_e_series=array(NA,dim=c(run_length,10))
#evolved_d_series=array(NA,dim=c(run_length,10))
#singular_es=array(NA,dim=c(run_length,10))
#singular_ds=array(NA,dim=c(run_length,10))
#evolved_v_series=array(NA,dim=c(run_length,10))
#evolved_s_series=array(NA,dim=c(run_length,10))
#evolved_m_series=array(NA,dim=c(run_length,10))

#Store the highest CSS in the first column, second-highest in the second column, etc.
#This code turns the resulting lists into arrays that I find easier to work with
for(i in 1:run_length){
  for(j in 1:10){
#    evolved_w_series[i,j]=-sort(-evo_w_series[[i]])[j]
#    evolved_e_series[i,j]=-sort(-evo_e_series[[i]])[j]
#    evolved_d_series[i,j]=-sort(-evo_d_series[[i]])[j]
#    singular_es[i,j]=-sort(-singular_e_series[[i]])[j]
#    singular_ds[i,j]=-sort(-singular_d_series[[i]])[j]
#    evolved_v_series[i,j]=-sort(-evo_v_series[[i]])[j]
#    evolved_s_series[i,j]=-sort(-evo_s_series[[i]])[j]
#    evolved_m_series[i,j]=-sort(-evo_m_series[[i]])[j]
  }
}
}

ES_checker<-function(w,e,th,d,fmax,fmin,fscale,v,s,m,f,K){
  #1 if ES, 0 if not
  return(ifelse(Diff2fM_NGM(f,w,K,e,th,d,fmax,fmin,fscale,v,s,m)<0,1,0))
}
CS_checker<-function(w,e,th,d,fmax,fmin,fscale,v,s,m,f,K){
  #1 if ES, 0 if not
  return(ifelse(Diff2fR_NGM(f,w,K,e,th,d,fmax,fmin,fscale,v,s,m)-Diff2fM_NGM(f,w,K,e,th,d,fmax,fmin,fscale,v,s,m)>0,1,0))
}

#Confirmed that these are repellors with same not CS, not ES properties as before.
#Then check equiilibrium stabilities.
#CS_checker(w_ex,e_series[which(!is.na(singular_es[,3]))],th_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex,singular_es[which(!is.na(singular_es[,3])),2],K_ex)
#CS_checker(w_ex,e_ex,th_ex,d_series[which(!is.na(singular_ds[,3]))],fmax,fmin,fscale,v_ex,s_ex,m_ex,singular_ds[which(!is.na(singular_ds[,3])),2],K_ex)

#ES_checker(w_ex,e_series[which(!is.na(singular_es[,3]))],th_ex,d_ex,fmax,fmin,fscale,v_ex,s_ex,m_ex,singular_es[which(!is.na(singular_es[,3])),2],K_ex)
#ES_checker(w_ex,e_ex,th_ex,d_series[which(!is.na(singular_ds[,3]))],fmax,fmin,fscale,v_ex,s_ex,m_ex,singular_ds[which(!is.na(singular_ds[,3])),2],K_ex)

#Save the model data to avoid redoing the long calculation
#saveRDS(evolved_w_series,"evolved_w_series.RDS")
#saveRDS(evolved_e_series,"evolved_e_series.RDS")
#saveRDS(evolved_d_series,"evolved_d_series.RDS")
#saveRDS(singular_es,"singular_es.RDS")
#saveRDS(singular_ds,"singular_ds.RDS")
#saveRDS(evolved_v_series,"evolved_v_series.RDS")
#saveRDS(evolved_s_series,"evolved_s_series.RDS")
#saveRDS(evolved_m_series,"evolved_m_series.RDS")

#Can read this data in if you have already done the computation
evolved_w_series=readRDS("evolved_w_series.RDS")
evolved_e_series=readRDS("evolved_e_series.RDS")
evolved_d_series=readRDS("evolved_d_series.RDS")
evolved_v_series=readRDS("evolved_v_series.RDS")
evolved_s_series=readRDS("evolved_s_series.RDS")
evolved_m_series=readRDS("evolved_m_series.RDS")
singular_es=readRDS("singular_es.RDS")
singular_ds=readRDS("singular_ds.RDS")


#H and prev functions defined
{
  Hpartial=function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    temp=Eigen_out(w,K,f,e,th,d,B,v,s,m)
    return(sum(temp[1:2]))
  }
  prevpartial=function(w,K,f,e,th,d,fmax,fmin,fscale,v,s,m){
    temp=Eigen_out(w,K,f,e,th,d,B,v,s,m)
    return(temp[1]/sum(temp[1:2]))
  }
}

#Correct some incorrect, low, irrelevant values.
singular_es[singular_es[,2]<0.006,2]<-NA

#Fix the issues in here with old Standard_parms[6] as a reference to fmax, etc.
#Re-check with higher v resolution
#Makes Figure A5 of the manuscript
#png("Appendix_plot_increase_parms.png",width=8,height=8,units="in",res=300)
{
  cex_smallest_text=.5
  cex_minor_text=.75
  cex_major_text=.9
  lwd_minor=1
  lwd_major=1.5
  par(mfrow=c(3,3),mar=c(2,2,0,0),mgp=c(1.8,.4,0),cex.lab=cex_major_text,oma=c(2,2,0,0),cex.axis=cex_minor_text,bg="white",cex=1)

  yB_all=c(0,45)
  yH_all=c(0,125)
    
plot(log10(e_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],evolved_e_series[,1]),type="l",lty=1,xlab="",ylim=yB_all,ylab="",lwd=lwd_major)
points(log10(e_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],singular_es[,3]),type="l",lwd=lwd_major)
points(log10(e_series)[which(!is.na(singular_es[,3]))],1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],singular_es[which(!is.na(singular_es[,3])),2]),type="l",lwd=lwd_major,lty=2)
mtext(side=2,expression(atop("Evolved transmission rate: ",italic(beta)[CSS]~"(L parasite"^"-1"~"day"^"-1"~" x 10"^"-6"~")")),cex=cex_major_text,padj=-.75)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_minor_text)

plot(log10(w_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],evolved_w_series[,1]),type="l",lty=1,xlab="",ylim=yB_all,ylab="",lwd=lwd_major)
points(log10(w_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],evolved_w_series[,2]),type="l",lwd=lwd_major)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"D",cex=cex_minor_text)

plot(log10(s_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],evolved_s_series[,1]),type="l",lty=1,xlab="",ylim=yB_all,ylab="",lwd=lwd_major)
points(log10(s_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],evolved_s_series[,2]),type="l",lwd=lwd_major)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"G",cex=cex_minor_text)

plot(log10(e_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_e_series[,1],e=e_series,th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",ylim=c(0,1),xlab="",ylab="",lwd=lwd_major)
points(log10(e_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=singular_es[,3],e=e_series,th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",lwd=lwd_major)
points(log10(e_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=singular_es[,2],e=e_series,th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",lwd=lwd_major,lty=2)
mtext(side=2,expression(atop("Prevalence: ",italic(p)[CSS]^"*"~"(unitless)")),cex=cex_major_text,padj=-.75)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_minor_text)

plot(log10(w_series),mapply(prevpartial,w=w_series,K=Standard_parms[2],f=evolved_w_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",ylim=c(0,1),xlab="",ylab="",lwd=lwd_major)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"E",cex=cex_minor_text)

plot(log10(s_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_s_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=s_series,m=Standard_parms[11]),type="l",ylim=c(0,1),xlab="",ylab="",lwd=lwd_major)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"H",cex=cex_minor_text)

plot(log10(e_series),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_e_series[,1],e=e_series,th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",xlab="",ylab="",ylim=yH_all,lwd=lwd_major,yaxt="n")
axis(side=2,at=c(0,50,100,150,200,250,300),labels=c("0","","100","","200","","300"))
points(log10(e_series),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=singular_es[,3],e=e_series,th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",lwd=lwd_major)
points(log10(e_series),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=singular_es[,2],e=e_series,th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",lwd=lwd_major,lty=2)
mtext(side=2,expression(atop("Host density: ",italic(H)[CSS]^"*"~"(Hosts L"^"-1"~")")),cex=cex_major_text,padj=-.75)
mtext(side=1,expression(atop("Conversion efficiency:","log"[10]~"("~italic(e)~") (host"~mu~"g chl"~italic(a)^"-1"~")")),cex=cex_major_text,padj=1.5)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_minor_text)

plot(log10(w_series),mapply(Hpartial,w=w_series,K=Standard_parms[2],f=evolved_w_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",xlab="",ylab="",ylim=yH_all,lwd=lwd_major)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"F",cex=cex_minor_text)
mtext(side=1,expression(atop("Intrinsic rate of increase of resource:","log"[10]~"("~italic(w)~") (day"^"-1"~")")),cex=cex_major_text,padj=1.5)

plot(log10(s_series),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_s_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=s_series,m=Standard_parms[11]),type="l",xlab="",ylab="",lwd=lwd_major,ylim=yH_all)
text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"I",cex=cex_minor_text)
mtext(side=1,expression(atop("Parasites released:","log"[10]~"("~italic(sigma)~") (parasites host"^"-1"~")")),cex=cex_major_text,padj=1.5)
}
#dev.off()

#Continue here by implementing more mapply()
#Makes Figure A3 of the manuscript
#png("Appendix_plot_decrease_parms.png",width=8,height=8,units="in",res=300)
{
  cex_smallest_text=.5
  cex_minor_text=.75
  cex_major_text=.9
  lwd_minor=1
  lwd_major=1.5
  par(mfrow=c(3,3),mar=c(2,2,0,0),mgp=c(1.8,.4,0),cex.lab=cex_major_text,oma=c(2,2,0,0),cex.axis=cex_minor_text,bg="white",cex=1)
  
  plot(log10(d_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],evolved_d_series[,1]),type="l",lty=1,xlab="",ylim=yB_all,ylab="",lwd=lwd_major)
  points(log10(d_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],singular_ds[,3]),type="l",lwd=lwd_major)
  points(log10(d_series)[which(!is.na(singular_ds[,3]))],1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],singular_ds[which(!is.na(singular_ds[,3])),2]),type="l",lwd=lwd_major,lty=2)
  mtext(side=2,expression(atop("Evolved transmission rate: ",italic(beta)[CSS]~"(L parasite"^"-1"~"day"^"-1"~" x 10"^"-6"~")")),cex=cex_major_text,padj=-.75)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"A",cex=cex_minor_text)
  
  plot(log10(m_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],evolved_m_series[,1]),type="l",lty=1,xlab="",ylim=yB_all,ylab="",lwd=lwd_major)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"D",cex=cex_minor_text)
  
  plot(log10(v_series),1e6*Beta(Standard_parms[6],Standard_parms[7],Standard_parms[8],evolved_v_series[,1]),type="l",lty=1,xlab="",ylim=yB_all,ylab="",lwd=lwd_major)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"G",cex=cex_minor_text)
  
  plot(log10(d_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_d_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=d_series,fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",ylim=c(0,1.1),xlab="",ylab="",lwd=lwd_major)
  points(log10(d_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=singular_ds[,3],e=Standard_parms[3],th=Standard_parms[4],d=d_series,fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",lwd=lwd_major)
  points(log10(d_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=singular_ds[,2],e=Standard_parms[3],th=Standard_parms[4],d=d_series,fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",lwd=lwd_major,lty=2)
  mtext(side=2,expression(atop("Prevalence: ",italic(p)[CSS]^"*"~"(unitless)")),cex=cex_major_text,padj=-.75)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"B",cex=cex_minor_text)
  
  plot(log10(m_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_m_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=m_series),type="l",ylim=c(0,1.1),xlab="",ylab="",lwd=lwd_major)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"E",cex=cex_minor_text)
  
  plot(log10(v_series),mapply(prevpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_v_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=v_series,s=Standard_parms[10],m=Standard_parms[11]),type="l",ylim=c(0,1.1),xlab="",ylab="",lwd=lwd_major)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"H",cex=cex_minor_text)
  
  plot(log10(d_series),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_d_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=d_series,fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",xlab="",ylab="",ylim=yH_all,lwd=lwd_major)
  points(log10(d_series),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=singular_ds[,3],e=Standard_parms[3],th=Standard_parms[4],d=d_series,fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",lwd=lwd_major)
  points(log10(d_series[which(!is.na(singular_ds[,3]))]),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=singular_ds[which(!is.na(singular_ds[,3])),2],e=Standard_parms[3],th=Standard_parms[4],d=d_series[which(!is.na(singular_ds[,3]))],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=Standard_parms[11]),type="l",lwd=lwd_major,lty=2)
  mtext(side=2,expression(atop("Host density: ",italic(H)[CSS]^"*"~"(Hosts L"^"-1"~")")),cex=cex_major_text,padj=-.75)
  mtext(side=1,expression(atop("Host background mortality:","log"[10]~"("~italic(d)~") ( day"^"-1"~")")),cex=cex_major_text,padj=1.5)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"C",cex=cex_minor_text)
  
  plot(log10(m_series),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_m_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=Standard_parms[9],s=Standard_parms[10],m=m_series),type="l",xlab="",ylab="",lwd=lwd_major,ylim=yH_all)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"F",cex=cex_minor_text)
  mtext(side=1,expression(atop("Loss rate of parasites:","log"[10]~"("~italic(m)~") (day"^"-1"~")")),cex=cex_major_text,padj=1.5)
  
  plot(log10(v_series),mapply(Hpartial,w=Standard_parms[1],K=Standard_parms[2],f=evolved_v_series[,1],e=Standard_parms[3],th=Standard_parms[4],d=Standard_parms[5],fmax=Standard_parms[6],fmin=Standard_parms[7],fscale=Standard_parms[8],v=v_series,s=Standard_parms[10],m=Standard_parms[11]),type="l",xlab="",ylab="",lwd=lwd_major,ylim=yH_all)
  text((par("usr")[2]-par("usr")[1])*.05+par("usr")[1],(par("usr")[4]-par("usr")[3])*.95+par("usr")[3],"I",cex=cex_minor_text)
  mtext(side=1,expression(atop("Mortality virulence:","log"[10]~"("~italic(v)~") (day"^"-1"~")")),cex=cex_major_text,padj=1.5)
}
#dev.off()


#Check which equilibrium is stable. All of these return equilibrium 6, as expected.
#wedvsm
w_stable=unique(mapply(Which_stable,w=w_series,K=K_ex,f=evolved_w_series[,1],e=e_ex,th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,evolved_w_series[,1]),v=v_ex,s=s_ex,m=m_ex))
e_stable1=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=evolved_e_series[which(!is.na(evolved_e_series[,1])),1],e=e_series[which(!is.na(evolved_e_series[,1]))],th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,evolved_e_series[which(!is.na(evolved_e_series[,1])),1]),v=v_ex,s=s_ex,m=m_ex))
e_stable2=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=evolved_e_series[which(!is.na(evolved_e_series[,2])),2],e=e_series[which(!is.na(evolved_e_series[,2]))],th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,evolved_e_series[which(!is.na(evolved_e_series[,2])),2]),v=v_ex,s=s_ex,m=m_ex))

e_stable3=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=singular_es[which(singular_es[,2]>0.006),2],e=e_series[which(singular_es[,2]>0.006)],th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,singular_es[which(singular_es[,2]>0.006),2]),v=v_ex,s=s_ex,m=m_ex))

d_stable1=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=evolved_d_series[which(!is.na(evolved_d_series[,1])),1],e=e_ex,th=th_ex,d=d_series[which(!is.na(evolved_d_series[,1]))],B=Beta(fmax,fmin,fscale,evolved_d_series[which(!is.na(evolved_d_series[,1])),1]),v=v_ex,s=s_ex,m=m_ex))
d_stable2=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=evolved_d_series[which(!is.na(evolved_d_series[,2])),2],e=e_ex,th=th_ex,d=d_series[which(!is.na(evolved_d_series[,2]))],B=Beta(fmax,fmin,fscale,evolved_d_series[which(!is.na(evolved_d_series[,2])),2]),v=v_ex,s=s_ex,m=m_ex))

d_stable3=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=singular_ds[which(singular_ds[,2]>0.006),2],e=e_ex,th=th_ex,d=d_series[which(singular_ds[,2]>0.006)],B=Beta(fmax,fmin,fscale,singular_ds[which(singular_ds[,2]>0.006),2]),v=v_ex,s=s_ex,m=m_ex))

v_stable=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=evolved_v_series[,1],e=e_ex,th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,evolved_v_series[,1]),v=v_series,s=s_ex,m=m_ex))
s_stable=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=evolved_s_series[,1],e=e_ex,th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,evolved_s_series[,1]),v=v_ex,s=s_series,m=m_ex))
m_stable=unique(mapply(Which_stable,w=w_ex,K=K_ex,f=evolved_m_series[,1],e=e_ex,th=th_ex,d=d_ex,B=Beta(fmax,fmin,fscale,evolved_m_series[,1]),v=v_ex,s=s_ex,m=m_series))
