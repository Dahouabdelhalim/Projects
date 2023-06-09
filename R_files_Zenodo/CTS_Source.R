###Source script for California Tiger Salamander (Ambystoma californiense) IPM
###Arianne Messerman
###September 2021
###Code adapted from Eldered and Miller (2016)

#setwd("")

##To calculate inverse logit functions
invlogit<-function(x){exp(x)/(1+exp(x))}

###########################################################################
##Vital rate functions. Parameter indices are hard-coded infunctions 
## and must correspond to rows of MCMC matrix as defined in IPM script
###########################################################################
##Metamorph growth from size x to size y
gxy.m<-function(x,y,params){
  xb=pmin(pmax(x,params[16]),params[17]) #Transforms all values below/above limits in min/max size
  xg=params[1] + params[2]*(xb-2.39) #Creates one mean growth value per body size bin per iteration of posterior draws, centered
  return(dnorm(y,mean=xg, sd=params[3]))
}

##Adult growth from size x to size y
gxy.a<-function(x,y,params){
  xb=pmin(pmax(x,params[16]),params[17]) 
  xg=params[4] + params[5]*(xb-2.95) 
  return(dnorm(y,mean=xg,sd=params[6]))
}

##Metamorph survival at size x
sx.m<-function(x,params){
  xb=pmin(pmax(x,params[16]),params[17])
  smeta<-(invlogit(params[7] + params[8]*(xb-2.28))/params[10])
  smeta<-ifelse(smeta>1, 1, smeta)
  return(smeta)
}

##Adult/juvenile survival at size x
sx.a<-function(x,params){
  xb=pmin(pmax(x,params[16]),params[17])
  sadj<-(invlogit(params[9] + params[8]*(xb-2.49))/params[11])
  sadj<-ifelse(sadj>1, 1, sadj)
  return(sadj)
}

#Metamorph survival*growth
pxy.m<-function(x,y,params){
  xb=pmin(pmax(x,params[16]),params[17])
  sx.m(xb,params)*gxy.m(xb,y,params)
}

#Adult survival*growth
pxy.a<-function(x,y,params){
  xb=pmin(pmax(x,params[16]),params[17])
  sx.a(xb,params)*gxy.a(xb,y,params)
}

#Fecundity (fx) is the product of maturity (mx), fertility (fer), percent females breeding (34%), 
# percent population that is female (50%), and maximum possible embryonic/larval survival
# representing the number of juveniles produced by a female of a size class (y)
mx<-function(y,params){
  yb=pmin(pmax(y,params[16]),params[17])
  invlogit(params[12] + params[13]*(yb-2.45))
}

fer<-function(y,params){
  yb=pmin(pmax(y,params[16]),params[17])
  return(params[14] + params[15]*(yb-3.48)) 
}

fx<-function(y,params){
  yb=pmin(pmax(y,params[16]),params[17])
  nfem<-params[20]
  nfem.b<-params[21]
  sx.l<-params[22]
  return(mx(yb, params)*fer(yb, params)*nfem*nfem.b*sx.l)
}
  
#Size distribution of new metamorphs under low densities
recruits<-function(y,params){
  yb=pmin(pmax(y,params[16]),params[17])
  dnorm(x=yb,mean=params[18],sd=params[19])
}


###Fecundity and body size at metamorphosis under density dependence
#Fecundity without embryonic/larval survival accounted for
fx1<-function(y,params){
  yb=pmin(pmax(y,params[16]),params[17])
  nfem<-params[20]
  nfem.b<-params[21]
  return(mx(yb, params)*fer(yb, params)*nfem*nfem.b)
}

#Function for metamorph body size given egg density
met<-function(d,params){
  return(params[24] + params[23]*(d-1.23)) #Mean-centered on Searcy et al 2015 data
}

#Function for larval survival to metamorphosis given egg density
lar<-function(d,params){
  lar.0<-params[25] + params[26]*(d-3.38) #Mean-centered on observed data
  return(ifelse(lar.0>=0, 0, lar.0))
}

#Metamorph survival*growth applied to current population size distribution (m)
pxy.m1<-function(x,y,m,params){
  xb=pmin(pmax(x,params[16]),params[17])
  (sx.m(xb,params)*m)*gxy.m(xb,y,params)
}

#Adult survival*growth applied to current population size distribution (a)
pxy.a1<-function(x,y,a,params){
  xb=pmin(pmax(x,params[16]),params[17])
  (sx.a(xb,params)*a)*gxy.a(xb,y,params)
}
#Size distribution of new metamorphs given egg densities (d)
recruits2<-function(y,d,params){ 
  yb=pmin(pmax(y,params[16]),params[17])
  return(dnorm(x=yb, mean = (params[24] + params[23]*(d-1.23)), sd=0.227))#0.227 set variance, mean-centered from Searcy et al 2015 data
}

###Fecundity under environmental and density dependence
#Fecundity without embryonic/larval survival accounted for
#Where percent females breeding is now predicted by Dec-Jan rainfall
nfem.b1<-function(e,params){
  return(min(exp(params[27] + params[28]*(e-221.09)), 1)) 
  #Where 221.09 mm is mean Dec-Jan precip from the study, and output is back-transformed from log and constrained to a maximum probability of 1
}

fx2<-function(y,e,params){
  yb=pmin(pmax(y,params[16]),params[17])
  nfem<-params[20]
  return(mx(yb, params)*fer(yb, params)*nfem*nfem.b1(e, params))
}
