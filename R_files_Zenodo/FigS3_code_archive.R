####Fig S3 simulation code

####Fig S3A: Vmax(t) and uncertainty as a function of keff
n=10000 #number of draws for bootstrapping

mu.Km=1042.4590 #this is for R. pseudoacacia (mL)
#Km is 467.1998 for A. rubra, 862.4747 for G. sepium, 856.7773 for M. cerifera

sigma.Km=62.60416 #this is for R. pseudoacacia (mL)
#sigma Km is 83.93578 for A. rubra, 87.17464 for G. sepium, 48.88134 for M. cerifera

shape.Km=mu.Km^2/sigma.Km^2
rate.Km=mu.Km/sigma.Km^2
Km=rgamma(n,shape.Km,rate.Km)

mu.A=817 #(mL) 2% injection of acetylene
sigma.A=10 #uncertainty (mL)
A=rnorm(n,mu.A,sigma.A)

t=0 #for simulations using time of zero

Vmax=688.4068 #example Vmax(t) for Fig S3A,B

E<-200 #(ppb) typical concentration of ethylene at 2% lab-made acetylene

#Solve for CV of keff to keep it constant for simulations across keff values
mu.keff<-0.004460 #(1/hr)
sigma.keff<-0.002518 #(1/hr)
CV.keff<-sigma.keff/mu.keff

keff.seq<-seq(0.0001,1,0.0001) #keff sequence (1/hr)
sigma.keff<-keff.seq*CV.keff #sigma keff sequence (1/hr)

#make keff distributions as a function of mean keff estimate with constant CV
keff<-matrix(NA,n,length(keff.seq))
for(i in 1:length(keff.seq)){
  shape.keff=keff.seq[i]^2/sigma.keff[i]^2
  rate.keff=keff.seq[i]/sigma.keff[i]^2
  keff[,i]<-rgamma(n,shape.keff,rate.keff)
}

#solve for how dE/dt changes as a function of keff so Vmax(t) is constant
dEdt.b<-rep(NA,length(keff.seq))
for(i in 1:length(keff.seq)){
  dEdt.b[i]<-quantile(Vmax/(1+Km/(A*exp(-keff.seq[i]*t)))-keff.seq[i]*E,.5)
}

sigma.slope=0.5 #assume same uncertainty for all dE/dt values

Vmax.keff.mean<-rep(NA,length(keff.seq))
Vmax.keff.low<-rep(NA,length(keff.seq))
Vmax.keff.high<-rep(NA,length(keff.seq))

for (i in 1:length(keff.seq)){
  slope=rnorm(n,dEdt.b[i],sigma.slope)
  Vmax.b=(slope+E*keff[,i])*(Km/(A*exp(-t*keff[,i]))+1)
  Vmax.keff.mean[i]<-mean(Vmax.b)
  Vmax.keff.low[i]<-quantile(Vmax.b,.025)
  Vmax.keff.high[i]<-quantile(Vmax.b,.975)
}

#assemble mean and bounds of 95% CI for Vmax(t) as a function of keff
Vmax.keff<-cbind(Vmax.keff.mean,Vmax.keff.low,Vmax.keff.high)
colnames(Vmax.keff)<-c("Mean","Low","High")

#########################################################################################################

####Fig S3B: Vmax(t) and uncertainty as a function of sigma keff

sigma.keff<-seq(0.0001,0.01,0.0001) #sigma keff sequence (1/hr)

#make keff distributions as a function of sigma keff (1/hr)
keff<-matrix(NA,n,length(sigma.keff))
for(i in 1:length(sigma.keff)){
  shape.keff=mu.keff^2/sigma.keff[i]^2
  rate.keff=mu.keff/sigma.keff[i]^2
  keff[,i]<-rgamma(n,shape.keff,rate.keff)
}

#solve for how dE/dt (ppb/hr) changes as a function of sigma keff so Vmax(t) is constant
dEdt.b<-quantile(Vmax/(1+Km/(A*exp(-mu.keff*t)))-mu.keff*E,.5)
dEdt.b<-dEdt.b[[1]]


Vmax.keff_SD.mean<-rep(NA,length(sigma.keff))
Vmax.keff_SD.low<-rep(NA,length(sigma.keff))
Vmax.keff_SD.high<-rep(NA,length(sigma.keff))

for (i in 1:length(sigma.keff)){
  slope=rnorm(n,dEdt.b,sigma.slope)
  Vmax.b=(slope+E*keff[,i])*(Km/(A*exp(-t*keff[,i]))+1)
  Vmax.keff_SD.mean[i]<-mean(Vmax.b)
  Vmax.keff_SD.low[i]<-quantile(Vmax.b,.025)
  Vmax.keff_SD.high[i]<-quantile(Vmax.b,.975)
}


#assemble mean and bounds of 95% CI for Vmax(t) as a function of sigma keff
Vmax.keff_SD<-cbind(Vmax.keff_SD.mean,Vmax.keff_SD.low,Vmax.keff_SD.high)
colnames(Vmax.keff_SD)<-c("Mean","Low","High")

######################################################################################################

####Fig S3C: Detection threshold as a function of keff

slope.seq<-seq(-1000,10,1) #(ppb/hr) slope sequence that works well with keff sequence

keff.seq<-c(seq(0.0001,0.0009,0.0001),seq(0.001,1,0.001)) #keff sequence

sigma.keff<-keff.seq*CV.keff #sigma keff sequence

#make keff distributions as a function of keff
keff<-matrix(NA,n,length(keff.seq))
for(i in 1:length(keff.seq)){
  shape.keff=keff.seq[i]^2/sigma.keff[i]^2
  rate.keff=keff.seq[i]/sigma.keff[i]^2
  keff[,i]<-rgamma(n,shape.keff,rate.keff)
}

#Find Vmax(t) where lower bound of 95% CI crosses zero as a function of keff
Vmax.d.m<-matrix(NA,length(keff.seq),length(slope.seq)) #mean Vmax(t) estimate
Vmax.d.low<-matrix(NA,length(keff.seq),length(slope.seq)) #lower bound of 95% CI
for (j in 1:length(slope.seq)){
  for (i in 1:length(keff.seq)){
    slope.d=rnorm(n,slope.seq[j],sigma.slope)
    Vmax.d=(slope.d+E*keff[,i])*(Km/(A*exp(-t*keff[,i]))+1)
    Vmax.d.low[i,j]<-quantile(Vmax.d,.025)
    Vmax.d.m[i,j]<-mean(Vmax.d)
  }
}

zero_Vmax.d.low<-rep(NA,length(keff.seq))
for (i in 1:length(keff.seq)){
  zero_Vmax.d.low[i]<-min(which(Vmax.d.low[i,]>0))
}

zero_Vmax.d.m<-rep(NA,length(keff.seq))
for (i in 1:length(keff.seq)){
  zero_Vmax.d.m[i]<-Vmax.d.m[i,zero_Vmax.d.low[i]]
}

#assemble detection limit of Vmax(t) as a function of keff
detection.data<-cbind(keff.seq,zero_Vmax.d.m)
colnames(detection.data)<-c("keff","Vmax.Detection")

######################################################################################################

####Fig S3D: Detection threshold as a function of sigma keff

slope.seq<-seq(-100,10,0.1) #(ppb/hr) slope sequence that works well with sigma keff sequence

sigma.keff<-seq(0.0001,0.1,0.0001) #sigma keff sequence

#make keff distributions as a function of sigma keff
keff<-matrix(NA,n,length(sigma.keff))
for(i in 1:length(sigma.keff)){
  shape.keff=mu.keff^2/sigma.keff[i]^2
  rate.keff=mu.keff/sigma.keff[i]^2
  keff[,i]<-rgamma(n,shape.keff,rate.keff)
}

#Find Vmax(t) where lower bound of 95% CI crosses zero as a function of sigma keff
Vmax.d.m<-matrix(NA,length(sigma.keff),length(slope.seq)) #mean Vmax(t) estimate
Vmax.d.low<-matrix(NA,length(sigma.keff),length(slope.seq)) #lower bound of 95% CI
for (j in 1:length(slope.seq)){
  for (i in 1:length(sigma.keff)){
    slope.d=rnorm(n,slope.seq[j],sigma.slope)
    Vmax.d=(slope.d+E*keff[,i])*(Km/(A*exp(-t*keff[,i]))+1)
    Vmax.d.low[i,j]<-quantile(Vmax.d,.025)
    Vmax.d.m[i,j]<-mean(Vmax.d)
  }
}

zero_Vmax.d.low<-rep(NA,length(sigma.keff))
for (i in 1:length(sigma.keff)){
  zero_Vmax.d.low[i]<-min(which(Vmax.d.low[i,]>0))
}

zero_Vmax.d.m<-rep(NA,length(sigma.keff))
for (i in 1:length(sigma.keff)){
  zero_Vmax.d.m[i]<-Vmax.d.m[i,zero_Vmax.d.low[i]]
}

#assemble detection limit of Vmax(t) as a function of sigma keff
detection.data<-cbind(sigma.keff,zero_Vmax.d.m)
colnames(detection.data)<-c("SD.keff","Vmax.Detection")