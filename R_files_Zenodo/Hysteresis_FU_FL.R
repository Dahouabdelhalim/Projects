##Robert Dunn
## Three stage urchin model: Test for hysteresis (same method as in Baskett & Salomon 2010)

require(deSolve)
require(lattice)

rm(list=ls())

##Define system of ODEs, A=kelp, Us=small urchin, Ul=large urchin, L=spiny lobster
threestage=function(t, state, parms){
  with(as.list(c(state,parms)),{
    dA= (r*(1-A/Ka)-(deltaUs*Us+deltaUm*Um+deltaUl*Ul))*A        #kelp dynamics (logistic growth, type I func resp)
    
    dUs = ((aM*deltaUm*Um+aL*deltaUl*Ul)*A)*(1-sigma+sigma*(Ul/Kul))-
      (gammaS + (L*deltaLs/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+MUs)*Us #small urchin dynamics
    #^^reproduction from large urchins via kelp conversion, loss due to growth, type II func resp, natural mortality)
    
    dUm = (gammaS*Us)- (gammaM + (L*deltaLm/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul)) + MUm)*Um 
    
    dUl = (gammaM*Um)-((L*deltaLl/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+FU+MUl)*Ul   #large urchin dynamics
    #^^growth from small urchins, loss via type II func resp, fishing & natural mortality)
    
    dL = ((b*(deltaLs*Us+deltaLm*Um+deltaLl*Ul)/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))-ML-FL)*L #spiny lobster dynamics     
    #^^conversion of lobsters to urchins via type II func resp, loss via fishing & natural mortality)
    return(list(c(dA, dUs, dUm, dUl, dL))) # return dn/dt as a list with each state variable as a column
  })  
}
# recruitment facilitation term:   *(1-sigma+sigma*(Ul/Kul))
#To estimate carrying capacity for large urchins, remove the recruitment facilitation term
# and run the model with no lobsters present. Equilibrium value of Ul is the new 
# value of Kul (large urchin carrying capacity).
######
##Starting in Kelp Forest State
######

###----------parameters----------------------###
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
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
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#basic plot
#plot(out$time,out$A, type="l", lwd=2, xlab="t", ylab="N", col="forestgreen", ylim=c(0,100),main="3 urchin stages")
#lines(out$time, out$Us, type="l", lwd=2, col="purple")
#lines(out$time,out$Um, type="l",lwd=2, col="blue")
#lines(out$time,out$Ul, type="l", lwd=2, col="darkred")
#lines(out$time,out$L, type="l", lwd=2, col="indianred2")

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of FU to run model over, goes funky with FU starting at 0.0
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
FL= 0.0     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
Kul=40        # large urchin carrying capacity- 
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

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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


##################################
###########PLOTS##################
##################################
pdf("HysteresisFU_FL_Kelp.pdf")
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


######
###Starting in Urchin Barren State
######
###----------parameters----------------------###
state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
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
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#basic plot
#plot(out$time,out$A, type="l", lwd=2, xlab="t", ylab="N", col="forestgreen", ylim=c(0,100),main="3 urchin stages")
#lines(out$time, out$Us, type="l", lwd=2, col="purple")
#lines(out$time,out$Um, type="l",lwd=2, col="blue")
#lines(out$time,out$Ul, type="l", lwd=2, col="darkred")
#lines(out$time,out$L, type="l", lwd=2, col="indianred2")

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of FU to run model over, goes funky with FU starting at 0.0
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
###Reverse path for urchin mortality######


FUrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
Kul=40        # large urchin carrying capacity- 
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
FL= 0.0     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
###Reverse path for Lobster mortality######
FUrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
r= 10          #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
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

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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

##################################
###########PLOTS##################
##################################
pdf("HysteresisFU_FL_Barrens.pdf")
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


################################################################
####Same thing, starting in urchin barren, but with fewer#######
####values of FU plotted so that the figure is better###########
################################################################
###Starting in Urchin Barren State##############################
################################################################
###----------parameters----------------------###
state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
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
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#basic plot
#plot(out$time,out$A, type="l", lwd=2, xlab="t", ylab="N", col="forestgreen", ylim=c(0,100),main="3 urchin stages")
#lines(out$time, out$Us, type="l", lwd=2, col="purple")
#lines(out$time,out$Um, type="l",lwd=2, col="blue")
#lines(out$time,out$Ul, type="l", lwd=2, col="darkred")
#lines(out$time,out$L, type="l", lwd=2, col="indianred2")

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=40)   #vector of FU to run model over, goes funky with FU starting at 0.0
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
###Reverse path for Lobster mortality######

FUrevseq<-seq(1.0,0.0,length.out=40)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
Kul=40        # large urchin carrying capacity- 
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
FL= 0.0     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=40)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
###Reverse path for urchin mortality######
FUrevseq<-seq(1.0,0.0,length.out=40)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
r= 10          #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
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

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda

###############
##Forward Path for Urchin mortality
FUseq<-seq(0.0,1.0,length=40)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FUseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUseq)){
  FU=FUseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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
###Reverse path for urchin mortality######
FUrevseq<-seq(1.0,0.0,length.out=40)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FUrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FUrevseq)){
  FU=FUrevseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
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


##################################
###########PLOTS##################
##################################
pdf("HysteresisFU_FL_Barrensv2.pdf")
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


############################################################
############################################################
##Old plots
##Forward path is darker line, reverse path is lighter shade###
#Kelp Plot
#pdf("HysteresisFU3xKelp.pdf")
matplot(FUseq,Nf2[,,1], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Kelp", xlab="Urchin Fishing Mortality", ylim=c(0,2500))
matlines(FUrevseq,Nr2[,,1], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FUseq,Nf1[,,1], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FUrevseq,Nr1[,,1], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FUseq,Nf3[,,1], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FUrevseq,Nr3[,,1], type="l", cex=0.5, lty=1, lwd=2, col="orange")
text(0.8, 700, labels="FL=0.0, Blue", cex=0.75)
text(0.8,550, labels="FL=0.25, Red", cex=0.75)
text(0.8,400, labels="FL=0.5, Orange", cex=0.75)
#dev.off()

###Hystersis plots 
pdf("HysteresisFU3x_0.65.pdf")
par(mfrow=c(3,2))
#Kelp
matplot(FUseq,Nf2[,,1], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Kelp", xlab="Urchin Fishing Mortality", ylim=c(0,2500))
matlines(FUrevseq,Nr2[,,1], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FUseq,Nf1[,,1], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FUrevseq,Nr1[,,1], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FUseq,Nf3[,,1], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FUrevseq,Nr3[,,1], type="l", cex=0.5, lty=1, lwd=2, col="orange")

#Small Urchin
matplot(FUseq,Nf2[,,2], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Small Urchins", xlab="Urchin Fishing Mortality", ylim=c(0,40))
matlines(FUrevseq,Nr2[,,2], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FUseq,Nf1[,,2], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FUrevseq,Nr1[,,2], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FUseq,Nf3[,,2], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FUrevseq,Nr3[,,2], type="l", cex=0.5, lty=1, lwd=2, col="orange")
#Med Urchin
matplot(FUseq,Nf2[,,3], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1, 
        ylab="Medium Urchins", xlab="Urchin Fishing Mortality", ylim=c(0,20))
matlines(FUrevseq,Nr2[,,3], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FUseq,Nf1[,,3], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FUrevseq,Nr1[,,3], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FUseq,Nf3[,,3], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FUrevseq,Nr3[,,3], type="l", cex=0.5, lty=1, lwd=2, col="orange")
#Large Urchin
matplot(FUseq,Nf2[,,4], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Large Urchins", xlab="Urchin Fishing Mortality", ylim=c(0,20))
matlines(FUrevseq,Nr2[,,4], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FUseq,Nf1[,,4], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FUrevseq,Nr1[,,4], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FUseq,Nf3[,,4], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FUrevseq,Nr3[,,4], type="l", cex=0.5, lty=1, lwd=2, col="orange")
#Lobster
matplot(FUseq,Nf2[,,5], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Lobsters", xlab="Urchin Fishing Mortality", ylim=c(0,20))
matlines(FUrevseq,Nr2[,,5], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FUseq,Nf1[,,5], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FUrevseq,Nr1[,,5], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FUseq,Nf3[,,5], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FUrevseq,Nr3[,,5], type="l", cex=0.5, lty=1, lwd=2, col="orange")
text(0.8, 15.5, labels="FL=0.0, Red", cex=0.75)
text(0.8, 14, labels="FL=0.25, Blue", cex=0.75)
text(0.8, 12.5, labels="FL=0.65, Orange", cex=0.75)
dev.off()

