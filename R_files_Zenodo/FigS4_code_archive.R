####Fig S4 simulation code

####Fig S4A: Vmax(t) and uncertainty as a function of Km
n=10000 #number of draws for bootstrapping

t=0 #for simulations using time of zero

#leak rate coefficient
shape.keff<-3.137647916
rate.keff<-703.5331092
keff=rgamma(n,shape.keff,rate.keff)


mu.A=817 #(mL) 2% injection of acetylene
sigma.A=10 #uncertainty (mL)
A=rnorm(n,mu.A,sigma.A)

Vmax=688.4068 #example Vmax(t) for Fig S4A,B
E=200 #(ppb) typical concentration of ethylene at 2% lab-made acetylene

Km.seq<-seq(40,4100,10) #(mL)
sigma.Km=62.60416 #this is for R. pseudoacacia (mL)
#sigma Km is 83.93578 for A. rubra, 87.17464 for G. sepium, 48.88134 for M. cerifera

#make Km distributions as a function of mean Km and species specific sigma Km
Km<-matrix(NA,n,length(Km.seq))
for(i in 1:length(Km.seq)){
  shape.Km=Km.seq[i]^2/sigma.Km^2
  rate.Km=Km.seq[i]/sigma.Km^2
  Km[,i]<-rgamma(n,shape.Km,rate.Km)
}

sigma.slope=0.5 #assume same uncertainty for all dE/dt values

#solve for how dE/dt changes as a function of Km so Vmax(t) is constant
dEdt<-rep(NA,length(Km.seq))
for(i in 1:length(Km.seq)){
  dEdt[i]<-quantile(Vmax/(1+Km.seq[i]/(A*exp(-keff*t)))-keff*E,.5)
}

Vmax.Km.mean<-rep(NA,length(Km.seq))
Vmax.Km.low<-rep(NA,length(Km.seq))
Vmax.Km.high<-rep(NA,length(Km.seq))

for (i in 1:length(Km.seq)){
  slope=rnorm(n,dEdt[i],sigma.slope)
  Vmax.b=(slope+E*keff)*(Km[,i]/(A*exp(-t*keff))+1)
  Vmax.Km.mean[i]<-mean(Vmax.b)
  Vmax.Km.low[i]<-quantile(Vmax.b,.025)
  Vmax.Km.high[i]<-quantile(Vmax.b,.975)
}

#assemble mean and bounds of 95% CI for Vmax(t) as a function of Km
Vmax.Km<-cbind(Vmax.Km.mean,Vmax.Km.low,Vmax.Km.high)
colnames(Vmax.Km)<-c("Mean","Low","High")

#################################################################################################

####Fig S4B: Vmax(t) and uncertainty as a function of sigma Km

sigma.Km<-seq(1,1000,1) #(mL) sigma Km sequence
mu.Km=1042.4590 #this is for R. pseudoacacia (mL)
#Km is 467.1998 for A. rubra, 862.4747 for G. sepium, 856.7773 for M. cerifera

#make Km distributions as a function of sigma Km and species specific Km
Km<-matrix(NA,n,length(sigma.Km))
for(i in 1:length(sigma.Km)){
  shape.Km=mu.Km^2/sigma.Km[i]^2
  rate.Km=mu.Km/sigma.Km[i]^2
  Km[,i]<-rgamma(n,shape.Km,rate.Km)
}

#solve for how dE/dt changes as a function of sigma Km so Vmax(t) is constant
dEdt.Km_SD<-quantile(Vmax/(1+mu.Km/(A*exp(-keff*t)))-keff*E,.5)
dEdt.Km_SD<-dEdt.Km_SD[[1]]

Vmax.Km_SD.mean<-rep(NA,length(sigma.Km))
Vmax.Km_SD.low<-rep(NA,length(sigma.Km))
Vmax.Km_SD.high<-rep(NA,length(sigma.Km))

for (i in 1:length(sigma.Km)){
  slope.Km_SD=rnorm(n,dEdt.Km_SD,sigma.slope)
  Vmax.Km_SD=(slope.Km_SD+E*keff)*(Km[,i]/(A*exp(-t*keff))+1)
  Vmax.Km_SD.mean[i]<-mean(Vmax.Km_SD)
  Vmax.Km_SD.low[i]<-quantile(Vmax.Km_SD,.025)
  Vmax.Km_SD.high[i]<-quantile(Vmax.Km_SD,.975)
}

#assemble mean and bounds of 95% CI for Vmax(t) as a function of sigma Km
Vmax.Km_SD<-cbind(Vmax.Km_SD.mean,Vmax.Km_SD.low,Vmax.Km_SD.high)
colnames(Vmax.Km_SD)<-c("Mean","Low","High")

####################################################################################################

####Fig S4C: Detection threshold as a function of Km

slope.seq<-seq(-10,10,0.1) #(ppb/hr) slope sequence that works well with Km sequence

sigma.Km=62.60416 #this is for R. pseudoacacia (mL)

#make Km distributions as a function of mean Km and species specific sigma Km
Km<-matrix(NA,n,length(Km.seq))
for(i in 1:length(Km.seq)){
  shape.Km=Km.seq[i]^2/sigma.Km^2
  rate.Km=Km.seq[i]/sigma.Km^2
  Km[,i]<-rgamma(n,shape.Km,rate.Km)
}

#Find Vmax(t) where lower bound of 95% CI crosses zero as a function of Km
Vmax.d.m<-matrix(NA,length(Km.seq),length(slope.seq))
Vmax.d.low<-matrix(NA,length(Km.seq),length(slope.seq))
for (j in 1:length(slope.seq)){
  for (i in 1:length(Km.seq)){
    slope=rnorm(n,slope.seq[j],sigma.slope)
    Vmax.d=(slope+E*keff)*(Km[,i]/(A*exp(-t*keff))+1)
    Vmax.d.low[i,j]<-quantile(Vmax.d,.025)
    Vmax.d.m[i,j]<-mean(Vmax.d)
  }
}

zero_Vmax.d.low<-rep(NA,length(Km.seq))
for (i in 1:length(Km.seq)){
  zero_Vmax.d.low[i]<-min(which(Vmax.d.low[i,]>0))
}

zero_Vmax.d.m<-rep(NA,length(Km.seq))
for (i in 1:length(Km.seq)){
  zero_Vmax.d.m[i]<-Vmax.d.m[i,zero_Vmax.d.low[i]]
}

#assemble detection threshold as a function of Km
detection.data<-cbind(Km.seq,zero_Vmax.d.m)
colnames(detection.data)<-c("Km","Vmax.Detection")

#########################################################################################################

####Fig S4D: Detection threshold as a function of sigma Km

#solve for how dE/dt changes as a function of sigma Km so Vmax(t) is constant
Km<-matrix(NA,n,length(sigma.Km))
for(i in 1:length(sigma.Km)){
  shape.Km=mu.Km^2/sigma.Km[i]^2
  rate.Km=mu.Km/sigma.Km[i]^2
  Km[,i]<-rgamma(n,shape.Km,rate.Km)
}

#Find Vmax(t) where lower bound of 95% CI crosses zero as a function of sigma Km
Vmax.d.m<-matrix(NA,length(sigma.Km),length(slope.seq))
Vmax.d.low<-matrix(NA,length(sigma.Km),length(slope.seq))
for (j in 1:length(slope.seq)){
  for (i in 1:length(sigma.Km)){
    slope=rnorm(n,slope.seq[j],sigma.slope)
    Vmax.d=(slope+E*keff)*(Km[,i]/(A*exp(-t*keff))+1)
    Vmax.d.low[i,j]<-quantile(Vmax.d,.025)
    Vmax.d.m[i,j]<-mean(Vmax.d)
  }
}

zero_Vmax.d.low<-rep(NA,length(sigma.Km))
for (i in 1:length(sigma.Km)){
  zero_Vmax.d.low[i]<-min(which(Vmax.d.low[i,]>0))
}

zero_Vmax.d.m<-rep(NA,length(sigma.Km))
for (i in 1:length(sigma.Km)){
  zero_Vmax.d.m[i]<-Vmax.d.m[i,zero_Vmax.d.low[i]]
}

#assemble detection threshold as a function of sigma Km
detection.data<-cbind(sigma.Km,zero_Vmax.d.m)
colnames(detection.data)<-c("SD.Km","Vmax.Detection")