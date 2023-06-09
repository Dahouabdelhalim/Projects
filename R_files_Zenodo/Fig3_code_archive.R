####Fig 3 simulation code

####Fig 3A: Vmax(t) and uncertainty as a function of ethylene
n=10000 #number of draws for bootstrapping

#leak rate coefficient
shape.keff<-3.137647916
rate.keff<-703.5331092
keff=rgamma(n,shape.keff,rate.keff) #(1/hr)

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

Vmax=688.4068 #example Vmax(t) for Fig 3A,B

E.seq<-seq(0,200000,100) #(ppb) ethylene sequence

#solve for how dE/dt changes as a function of ethylene so Vmax(t) is constant
dEdt.b<-rep(NA,length(E.seq))
for(i in 1:length(E.seq)){
  dEdt.b[i]<-quantile(Vmax/(1+Km/(A*exp(-keff*t)))-keff*E.seq[i],.5)
}

sigma.slope=0.5 #assume same uncertainty for all dE/dt values

Vmax.E.mean<-rep(NA,length(E.seq))
Vmax.E.low<-rep(NA,length(E.seq))
Vmax.E.high<-rep(NA,length(E.seq))

for (i in 1:length(E.seq)){
  slope=rnorm(n,dEdt.b[i],sigma.slope)
  Vmax.b=(slope+E.seq[i]*keff)*(Km/(A*exp(-t*keff))+1)
  Vmax.E.mean[i]<-mean(Vmax.b)
  Vmax.E.low[i]<-quantile(Vmax.b,.025)
  Vmax.E.high[i]<-quantile(Vmax.b,.975)
}

#assemble mean and bounds of 95% CI for Vmax(t) as a function of ethylene
Vmax.E<-cbind(Vmax.E.mean,Vmax.E.low,Vmax.E.high)
colnames(Vmax.E)<-c("Mean","Low","High")

##########################################################################################

####Fig 3B: Vmax(t) and uncertainty as a function of acetylene

E<-200 #(ppb) typical concentration of ethylene at 2% lab-made acetylene

A.seq<-seq(40,4100,10) #(mL) sequence going from 0.1 to just past 10% acetylene

#sigma.A: 10mL below 1500mL, 20mL 1500-3000mL, and 30mL for 3000-4100mL (based on 1500 mL syringe used for injections)
sigma.A<-rep(NA,length(A.seq))
for (i in 1:length(A.seq)){
  if (A.seq[i]<=1500){
    sigma.A[i]<-10
  }
  else if (A.seq[i]>1500 & A.seq[i]<=3000){
    sigma.A[i]<-20
  }
  else{
    sigma.A[i]<-30
  }
}

A<-matrix(NA,n,length(A.seq))

for(i in 1:length(A.seq)){
  A[,i]=rnorm(n,A.seq[i],sigma.A[i])
}

E<-A.seq*0.25 #800ml acetylene roughly equals 200 ppb ethylene

#solve for how dE/dt changes as a function of acetylene so Vmax(t) is constant
dEdt.b<-rep(NA,length(A.seq))
for(i in 1:length(A.seq)){
  dEdt.b[i]<-quantile(Vmax/(1+Km/(A.seq[i]*exp(-keff*t)))-keff*E[i],.5)
}

Vmax.A.mean<-rep(NA,length(A.seq))
Vmax.A.low<-rep(NA,length(A.seq))
Vmax.A.high<-rep(NA,length(A.seq))

for (i in 1:length(A.seq)){
  slope=rnorm(n,dEdt.b[i],sigma.slope)
  Vmax.b=(slope+E[i]*keff)*(Km/(A[,i]*exp(-t*keff))+1)
  Vmax.A.mean[i]<-mean(Vmax.b)
  Vmax.A.low[i]<-quantile(Vmax.b,.025)
  Vmax.A.high[i]<-quantile(Vmax.b,.975)
}

#assemble mean and bounds of 95% CI for Vmax(t) as a function of acetylene
Vmax.A<-cbind(Vmax.A.mean,Vmax.A.low,Vmax.A.high)
colnames(Vmax.A)<-c("Mean","Low","High")

#####################################################################################

###Fig 3C Detection threshold as a function of ethylene

slope.seq<-seq(-1000,10,1) #(ppb/hr) slope sequence that works well with ethylene sequence

sigma.A=10 #(mL) acetylene uncertainty at 2% acetylene

E.seq<-seq(10,388,1)^2 #(ppb) ethylene sequence

A=rnorm(n,mu.A,sigma.A)

Vmax.d.m<-matrix(NA,length(E.seq),length(slope.seq)) #mean Vmax(t) estimate
Vmax.d.low<-matrix(NA,length(E.seq),length(slope.seq)) #lower bound of 95% CI

#Find Vmax(t) where lower bound of 95% CI crosses zero as a function of ethylene
for (j in 1:length(slope.seq)){
  for (i in 1:length(E.seq)){
    slope=rnorm(n,slope.seq[j],sigma.slope)
    Vmax.d=(slope+E.seq[i]*keff)*(Km/(A*exp(-t*keff))+1)
    Vmax.d.low[i,j]<-quantile(Vmax.d,.025)
    Vmax.d.m[i,j]<-mean(Vmax.d)
  }
}

zero_Vmax.d.low<-rep(NA,length(E.seq))
for (i in 1:length(E.seq)){
  zero_Vmax.d.low[i]<-min(which(Vmax.d.low[i,]>0))
}

zero_Vmax.d.m<-rep(NA,length(E.seq))
for (i in 1:length(E.seq)){
  zero_Vmax.d.m[i]<-Vmax.d.m[i,zero_Vmax.d.low[i]]
}

#assemble detection limit of Vmax(t) as a function of ethylene
detection.data.E<-cbind(E.seq,zero_Vmax.d.m)
colnames(detection.data.E)<-c("Ethylene","Vmax.Detection")

#####################################################################################

###Fig 3D Detection threshold as a function of acetylene

slope.seq<-seq(-10,10,0.1) #slope sequence that works well with acetylene sequence

#sigma.A: 10mL below 1500mL, 20mL 1500-3000mL, and 30mL for 3000-4100mL (based on 1500 mL syringe used for injections)
sigma.A<-rep(NA,length(A.seq))
for (i in 1:length(A.seq)){
  if (A.seq[i]<=1500){
    sigma.A[i]<-10
  }
  else if (A.seq[i]>1500 & A.seq[i]<=3000){
    sigma.A[i]<-20
  }
  else{
    sigma.A[i]<-30
  }
}

A<-matrix(NA,n,length(A.seq))

for(i in 1:length(A.seq)){
  A[,i]=rnorm(n,A.seq[i],sigma.A[i])
}

Vmax.d.m<-matrix(NA,length(A.seq),length(slope.seq))
Vmax.d.low<-matrix(NA,length(A.seq),length(slope.seq))

for (j in 1:length(slope.seq)){
  for (i in 1:length(A.seq)){
    slope.d=rnorm(n,slope.seq[j],sigma.slope)
    Vmax.d=(slope.d+E[i]*keff)*(Km/(A[,i]*exp(-t*keff))+1)
    Vmax.d.low[i,j]<-quantile(Vmax.d,.025)
    Vmax.d.m[i,j]<-mean(Vmax.d)
  }
}

zero_Vmax.d.low<-rep(NA,length(A.seq))
for (i in 1:length(A.seq)){
  zero_Vmax.d.low[i]<-min(which(Vmax.d.low[i,]>0))
}

zero_Vmax.d.m<-rep(NA,length(A.seq))
for (i in 1:length(A.seq)){
  zero_Vmax.d.m[i]<-Vmax.d.m[i,zero_Vmax.d.low[i]]
}

#assemble detection limit of Vmax(t) as a function of acetylene
detection.data.A<-cbind(A.seq,zero_Vmax.d.m)
colnames(detection.data.A)<-c("Acetylene","Vmax.Detection")