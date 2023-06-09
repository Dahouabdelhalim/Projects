##Robert Dunn
## Three stage urchin baseline model


require(deSolve)
require(lattice)
rm(list=ls())
setwd("~PhD/Community Dyn. Model/Plots")

##Define system of ODEs, A=kelp, Us=small urchin, Ul=large urchin, L=spiny lobster
threestage=function(t, state, parms){
  with(as.list(c(state,parms)),{
    dA= (r*(1-A/Ka)-(deltaUs*Us+deltaUm*Um+deltaUl*Ul))*A        #kelp dynamics (logistic growth, type I func resp)
    
    dUs = (((aM*deltaUm*Um+aL*deltaUl*Ul)*A) +phi)*(1-sigma+sigma*(Ul/Kul))-(gammaS + (L*deltaLs/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+MUs)*Us #small urchin dynamics
    #^^reproduction from large urchins via kelp conversion, loss due to growth, type II func resp, natural mortality), e is external recruitment
    
    dUm = (gammaS*Us)- (gammaM + (L*deltaLm/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul)) + MUm)*Um 
    
    dUl = (gammaM*Um)-((L*deltaLl/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+FU+MUl)*Ul   #large urchin dynamics
    #^^growth from small urchins, loss via type II func resp, fishing & natural mortality)
    
    dL = ((b*(deltaLs*Us+deltaLm*Um+deltaLl*Ul)/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))-ML-FL)*L          #spiny lobster dynamics     
    #^^conversion of urchins to lobsters via type II func resp, loss via fishing & natural mortality)
    return(list(c(dA, dUs, dUm, dUl, dL))) # return dn/dt as a list with each state variable as a column
  })  
}


##################################################
###----------parameters----------------------###
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system- kelp forest
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
phi=10            #external recruitment constant- how much to include? See Ebert et al. 1994, Okamoto dissertation 2014
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2 #0.001 for identical sim  #rate of small urchin consumption by lobsters .2 corresponds to 12mm urchin
deltaLm=0.15 #0.001 for identical sim #rate of med urchin consumption by lobsters  .0071 correponds to 40mm urchin
deltaLl=0.05 #0.001 for identical sim #rate of large urchin consumption by lobsters .0009 correponds to 82.5mm urchin
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL, omega=omega) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run model as a check
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#basic plot check the timeseries
plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen",
     ylim=c(0,100),xlim=c(0,25), main="3 urchin stages")
lines(out$time, out$Us, type="l", col="purple")
lines(out$time,out$Um, type="l", col="blue")
lines(out$time,out$Ul, type="l", col="darkred")
lines(out$time,out$L, type="l", col="indianred2")

####
###Now start in Urchin Barren state
####
state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Alternate barrens state

out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#basic plot
plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen",
     ylim=c(0,100), xlim=c(0,25),main="3 urchin stages")
lines(out$time, out$Us, type="l", col="purple")
lines(out$time,out$Um, type="l", col="blue")
lines(out$time,out$Ul, type="l", col="darkred")
lines(out$time,out$L, type="l", col="indianred2")



###########################
##### Hysteresis ##########
######  FL x FU  ##########

##Forward Path for Lobster mortality
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity
phi=10            #external recruitment
sigma= 0.5      #recruitment facilitation-  [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.0     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff= 1e-3  # difference between final equilibria and new starting conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)

FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1)
tail(out,5)

#####################
###Reverse path for Lobster mortality######

FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matlines(FLrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")


#plot
matplot(FLseq,Nf[,,1], type="l", cex=0.5, lty=3, lwd=2, col=1:10)
matlines(FLrevseq,Nr[,,1], type="l", cex=0.5, lty=5, lwd=2, col=11:20)

Nf1<-Nf         #FU = 0.0
Nr1<-Nr         #FU=0.0

####FU=0.1
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity-
phi=10
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff= 1e-3  # difference between final equilibria and new starting conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 

##Forward Path for Lobster mortality
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1)
tail(out,5)

#####################
###Reverse path for Lobster mortality######

FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matlines(FLrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")

Nf2<-Nf           #FU=0.10
Nr2<-Nr          #FU=0.10


###FU=0.25
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity
phi=10
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.25     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff= 1e-3  # difference between final equilibria and new starting conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 

##Forward Path for Lobster mortality
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1)
tail(out,5)

#####################
###Reverse path for Lobster mortality######

FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matlines(FLrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")

Nf3<-Nf         #FU=0.25
Nr3<-Nr          #FU=0.25


#####################################
##PLOTS#################
##############################
###Zoomed In, updated colors/points
pdf("HysteresisFL_FUexternal.pdf")
#pdf("test.pdf")
par(mfrow=c(3,2))
#Kelp
matplot(FLseq,Nf1[,,1], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Kelp (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,1000))
matlines(FLseq,Nf1[,,1], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,1], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,1], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,1], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,1], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,1], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,1], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,1], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Small Urchin
matplot(FLseq,Nf1[,,2], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Small urchins (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,40))
matlines(FLseq,Nf1[,,2], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,2], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,2], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,2], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,2], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,2], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,2], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,2], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Med Urchin
matplot(FLseq,Nf1[,,3], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Medium urchins (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,20))
matlines(FLseq,Nf1[,,3], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,3], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,3], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,3], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,3], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,3], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,3], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,3], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Large Urchin
matplot(FLseq,Nf1[,,4], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Large urchins (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,20))
matlines(FLseq,Nf1[,,4], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,4], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,4], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,4], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,4], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,4], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,4], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,4], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Lobster
matplot(FLseq,Nf1[,,5], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Lobsters (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,20))
matlines(FLseq,Nf1[,,5], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,5], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,5], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,5], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,5], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,5], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,5], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,5], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")
text(0.66, 18, labels="FU=0.0, Blue", cex=0.75)
text(0.66, 15, labels="FU=0.1, Red", cex=0.75)
text(0.66, 12, labels="FU=0.25, Black", cex=0.75)
dev.off()


##############################
######Hysteresis ############
####### FU x FL  #############

state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity- 
phi=10
sigma= 0.5      #recruitment facilitation-[Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of FU to run model over, goes funky with FU starting at 0.0
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FUseq,Nf[,,1], type="p", cex=0.5, pch=1)
tail(out,5)

#####################
###Reverse path for Urchin mortality######

FUrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matlines(FUrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")

Nf1<-Nf         #FL = 0.25
Nr1<-Nr         #FL=0.25

###FL=0.0##########
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity-
phi=10
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.0     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FUseq,Nf[,,1], type="p", cex=0.5, pch=1)
tail(out,5)

#####################
###Reverse path for Urchin mortality######
FUrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matlines(FUrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")

Nf2<-Nf           #FL=0.0
Nr2<-Nr          #FL=0.0

###FL=0.65###
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity- 
phi=10
sigma= 0.5      #recruitment facilitation-  [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.65     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FUseq,Nf[,,1], type="p", cex=0.5, pch=1, ylim=c(0,2500))
tail(out,5)

#####################
###Reverse path for Urchin mortality######
FUrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matlines(FUrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")

Nf3<-Nf         #FL=0.65
Nr3<-Nr          #FL=0.65

#Plots
pdf("HysteresisFU_FL_external_Kelp.pdf")
#pdf("test.pdf")
par(mfrow=c(3,2))
#Kelp
matplot(FUseq,Nf2[,,1], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Kelp (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,2500))
matlines(FUseq,Nf2[,,1], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,1], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,1], type="l", cex=0.5, lty=1, lwd=.5, col="lightblue3")
matlines(FUseq,Nf1[,,1], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,1], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,1], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,1], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,1], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Small Urchin
matplot(FUseq,Nf2[,,2], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Small urchins (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,40))
matlines(FUseq,Nf2[,,2], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,2], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,2], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FUseq,Nf1[,,2], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,2], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,2], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,2], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,2], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Med Urchin
matplot(FUseq,Nf2[,,3], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Medium urchins (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,20))
matlines(FUseq,Nf2[,,3], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,3], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,3], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FUseq,Nf1[,,3], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,3], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,3], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,3], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,3], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Large Urchin
matplot(FUseq,Nf2[,,4], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Large urchins (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,20))
matlines(FUseq,Nf2[,,4], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,4], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,4], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FUseq,Nf1[,,4], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,4], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,4], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,4], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,4], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Lobster
matplot(FUseq,Nf2[,,5], type="p", pch=6, cex=.75, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Lobsters (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,20))
matlines(FUseq,Nf2[,,5], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,5], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,5], type="l", cex=0.5, lty=1, lwd=.5, col="lightblue3")
matlines(FUseq,Nf1[,,5], type="p",  pch=6,cex=.75, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,5], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,5], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,5], type="p",  pch=6, cex=.75, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,5], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")
text(0.74, 58, labels="FL=0.0, Blue", cex=0.75)
text(0.74, 51, labels="FL=0.3, Red", cex=0.75)
text(0.74, 45, labels="FL=0.4, Black", cex=0.75)
dev.off()

##################################
###Hysteresis FU x FL start from barrens state

###Only run 40 points instead of 100 so plot is less cluttered
state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity-
phi=10
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=40)   #vector of FU to run model over, goes funky with FU starting at 0.0
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FUseq,Nf[,,1], type="p", cex=0.5, pch=1)
tail(out,5)

#####################
###Reverse path for Urchin mortality######

FUrevseq<-seq(1.0,0.0,length.out=40)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matlines(FUrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")

Nf1<-Nf         #FL = 0.25
Nr1<-Nr         #FL=0.25

###FL=0.0##########
state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity-
phi=10
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.0     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=40)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FUseq,Nf[,,1], type="p", cex=0.5, pch=1)
tail(out,5)

#####################
###Reverse path for Urchin mortality######
FUrevseq<-seq(1.0,0.0,length.out=40)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matlines(FUrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")

Nf2<-Nf           #FL=0.0
Nr2<-Nr          #FL=0.0

###FL=0.65###
state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
phi=10
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2   #rate of small urchin consumption by lobsters
deltaLm=0.15  #rate of med urchin consumption by lobsters
deltaLl=0.05   #rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.65     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=40)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FUseq,Nf[,,1], type="p", cex=0.5, pch=1, ylim=c(0,2500))
tail(out,5)

#####################
###Reverse path for Urchin mortality######
FUrevseq<-seq(1.0,0.0,length.out=40)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matlines(FUrevseq,Nr[,,1],type="p",cex=0.5,pch=1, col="red")

Nf3<-Nf         #FL=0.65
Nr3<-Nr          #FL=0.65

#Plots
pdf("HysteresisFU_FL_external_Barren.pdf")
#pdf("test.pdf")
par(mfrow=c(3,2))
#Kelp
matplot(FUseq,Nf2[,,1], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Kelp (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,2500))
matlines(FUseq,Nf2[,,1], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,1], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,1], type="l", cex=0.5, lty=1, lwd=.5, col="lightblue3")
matlines(FUseq,Nf1[,,1], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,1], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,1], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,1], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,1], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Small Urchin
matplot(FUseq,Nf2[,,2], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Small urchins (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,40))
matlines(FUseq,Nf2[,,2], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,2], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,2], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FUseq,Nf1[,,2], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,2], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,2], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,2], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,2], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Med Urchin
matplot(FUseq,Nf2[,,3], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Medium urchins (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,20))
matlines(FUseq,Nf2[,,3], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,3], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,3], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FUseq,Nf1[,,3], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,3], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,3], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,3], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,3], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Large Urchin
matplot(FUseq,Nf2[,,4], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Large urchins (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,20))
matlines(FUseq,Nf2[,,4], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,4], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,4], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FUseq,Nf1[,,4], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,4], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,4], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,4], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,4], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Lobster
matplot(FUseq,Nf2[,,5], type="p", pch=6, cex=.75, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Lobsters (kg)", xlab="Urchin Fishing Mortality", xlim=c(0,1), ylim=c(0,20))
matlines(FUseq,Nf2[,,5], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FUrevseq,Nr2[,,5], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FUrevseq,Nr2[,,5], type="l", cex=0.5, lty=1, lwd=.5, col="lightblue3")
matlines(FUseq,Nf1[,,5], type="p",  pch=6,cex=.75, lty=1, lwd=.5, col="darkred")
matlines(FUseq,Nf1[,,5], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FUrevseq,Nr1[,,5], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUrevseq,Nr1[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FUseq,Nf3[,,5], type="p",  pch=6, cex=.75, lty=1, lwd=.5, col="black")
matlines(FUseq,Nf3[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FUrevseq,Nr3[,,5], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FUrevseq,Nr3[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")
text(0.74, 58, labels="FL=0.0, Blue", cex=0.75)
text(0.74, 51, labels="FL=0.3, Red", cex=0.75)
text(0.74, 45, labels="FL=0.4, Black", cex=0.75)
dev.off()


###########
##############################
######Hysteresis ############
####### FL x Sigma  #############

state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
#state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Alternate barrens state
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity- 
phi=10
sigma= 0.0      #recruitment facilitation-  [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2  #rate of small urchin consumption by lobsters .2 corresponds to 12mm urchin
deltaLm=0.15  ##for identical sim #rate of med urchin consumption by lobsters  .0071 correponds to 40mm urchin
deltaLl=0.05  # for identical sim #rate of large urchin consumption by lobsters .0009 correponds to 82.5mm urchin
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for FL with 3 values of Urchin recruitment facilitation
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1, ylim=c(0,1500))
tail(out,5)

#####################
###Reverse path for Lobster mortality######

FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matpoints(FLrevseq,Nr[,,1],col="red", type="p",cex=0.5, pch=1)
tail(out,5)
#plot
#matplot(sigmaseq,Nf[,,1], type="l", cex=0.5, lty=3, lwd=2, col=1:10)
#matlines(sigmarevseq,Nr[,,1], type="l", cex=0.5, lty=5, lwd=2, col=11:20)

Nf1<-Nf         #sigma = 0.0
Nr1<-Nr         #sigma=0.0

#####Sigma =0.5###########
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity-
phi=10
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2 #0.001 for identical sim  #rate of small urchin consumption by lobsters .2 corresponds to 12mm urchin
deltaLm=0.15 #0.001 for identical sim #rate of med urchin consumption by lobsters  .0071 correponds to 40mm urchin
deltaLl=0.05 #0.001 for identical sim #rate of large urchin consumption by lobsters .0009 correponds to 82.5mm urchin
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

##Forward Path for FL with 3 values of Urchin recruitment facilitation
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1, ylim=c(0,1500))
tail(out,5)

#####################
###Reverse path for Lobster mortality######

FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matpoints(FLrevseq,Nr[,,1],col="red", type="p",cex=0.5, pch=1)
tail(out,5)

Nf2<-Nf           #sigma=0.5
Nr2<-Nr          #sigma=0.5


#######Sigma = 0.95#######
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity- 
phi=10
sigma= 0.95      #recruitment facilitation-  [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2 #0.001 for identical sim  #rate of small urchin consumption by lobsters .2 corresponds to 12mm urchin
deltaLm=0.15 #0.001 for identical sim #rate of med urchin consumption by lobsters  .0071 correponds to 40mm urchin
deltaLl=0.05 #0.001 for identical sim #rate of large urchin consumption by lobsters .0009 correponds to 82.5mm urchin
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but that was too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but that was too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

##Forward Path for FL with 3 values of Urchin recruitment facilitation
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nf[i,,1]=out[191:200,2]   #kelp pop matrix, forward path
  Nf[i,,2]=out[191:200,3]   #Us pop matrix
  Nf[i,,3]=out[191:200,4]   #Um pop matrix
  Nf[i,,4]=out[191:200,5]   #Ul pop matrix
  Nf[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
}
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1, ylim=c(0,1500))
tail(out,5)

#####################
###Reverse path for Lobster mortality######
FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, phi=phi, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))
  Nr[i,,1]=out[191:200,2]   #kelp pop matrix, reverse path
  Nr[i,,2]=out[191:200,3]   #Us pop matrix
  Nr[i,,3]=out[191:200,4]   #Um pop matrix
  Nr[i,,4]=out[191:200,5]   #Ul pop matrix
  Nr[i,,5]=out[191:200,6]   #lobster pop matrix
  state<-diff + c(A=out[200,2], Us=out[200,3], Um=out[200,4], Ul=out[200,5], L=out[200,6])
  
}
matpoints(FLrevseq,Nr[,,1],col="red", type="p",cex=0.5, pch=1)
tail(out,5)

Nf3<-Nf         #sigma=0.95
Nr3<-Nr          #sigma=0.95

#Plots
pdf("Hysteresis_FL_Sigma_External.pdf")
#pdf("test.pdf")
par(mfrow=c(3,2))
#Kelp
matplot(FLseq,Nf1[,,1], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Kelp (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,1000))
matlines(FLseq,Nf1[,,1], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,1], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,1], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,1], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,1], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,1], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,1], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,1], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,1], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Small Urchin
matplot(FLseq,Nf1[,,2], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Small urchins (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,50))
matlines(FLseq,Nf1[,,2], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,2], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,2], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,2], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,2], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,2], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,2], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,2], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,2], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Med Urchin
matplot(FLseq,Nf1[,,3], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Medium urchins (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,20))
matlines(FLseq,Nf1[,,3], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,3], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,3], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,3], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,3], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,3], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,3], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,3], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,3], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Large Urchin
matplot(FLseq,Nf1[,,4], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Large urchins (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,20))
matlines(FLseq,Nf1[,,4], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,4], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,4], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,4], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,4], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,4], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,4], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,4], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,4], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")

#Lobster
matplot(FLseq,Nf1[,,5], type="p", pch=6, cex=1, lty=1, lwd=.5, col="midnightblue", bty="L", las=1,
        ylab="Lobsters (kg)", xlab="Lobster Fishing Mortality", xlim=c(0.5,0.75), ylim=c(0,20))
matlines(FLseq,Nf1[,,5], type="l", cex=0.75, lty=1, lwd=.5, col="midnightblue")
matlines(FLrevseq,Nr1[,,5], type="p", pch=2, cex=0.75, lty=1, lwd=.5, col="lightblue3")
matlines(FLrevseq,Nr1[,,5], type="l", cex=0.5, lty=1, lwd=1.1, col="lightblue3")
matlines(FLseq,Nf2[,,5], type="p",  pch=6,cex=1, lty=1, lwd=.5, col="darkred")
matlines(FLseq,Nf2[,,5], type="l", cex=0.75, lty=1, lwd=.5, col="darkred")
matlines(FLrevseq,Nr2[,,5], type="p",  pch=2, cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLrevseq,Nr2[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="indianred1")
matlines(FLseq,Nf3[,,5], type="p",  pch=6, cex=1, lty=1, lwd=.5, col="black")
matlines(FLseq,Nf3[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="black")
matlines(FLrevseq,Nr3[,,5], type="p",  pch=2,cex=0.75, lty=1, lwd=.5, col="darkgrey")
matlines(FLrevseq,Nr3[,,5], type="l",  cex=0.75, lty=1, lwd=.5, col="darkgrey")
text(0.67, 18, labels="sigma=0.0, Blue", cex=0.75)
text(0.67, 15, labels="sigma=0.5, Red", cex=0.75)
text(0.67, 12, labels="sigma=0.95, Black", cex=0.75)
dev.off()
