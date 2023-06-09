##Robert Dunn
## Three stage urchin model: Test for hysteresis (same method as in Baskett&Salomon2010)

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
###----------parameters----------------------###
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
#state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Alternate barrens state
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
sigma= 0.0      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2  #rate of small urchin consumption by lobsters 
deltaLm=0.15  ## #rate of med urchin consumption by lobsters  
deltaLl=0.05  #  #rate of large urchin consumption by lobsters 
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
plot(out$time,out$A, type="l", lwd=2, xlab="t", ylab="N", col="forestgreen", ylim=c(0,100),main="3 urchin stages")
lines(out$time, out$Us, type="l", lwd=2, col="purple")
lines(out$time,out$Um, type="l",lwd=2, col="blue")
lines(out$time,out$Ul, type="l", lwd=2, col="darkred")
lines(out$time,out$L, type="l", lwd=2, col="indianred2")

###############
##Forward Path for FL with 3 values of Urchin recruitment facilitation
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
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
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1, ylim=c(0,1500))
tail(out,5)

#####################
###Reverse path for Lobster mortality######

FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
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
sigma= 0.5      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2 #  #rate of small urchin consumption by lobsters 
deltaLm=0.15 # #rate of med urchin consumption by lobsters  
deltaLl=0.05 # #rate of large urchin consumption by lobsters 
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

##Forward Path for FL with 3 values of Urchin recruitment facilitation
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
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
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1, ylim=c(0,1500))
tail(out,5)

#####################
###Reverse path for Lobster mortality######

FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
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
matpoints(FLrevseq,Nr[,,1],col="red", type="p",cex=0.5, pch=1)
tail(out,5)

Nf2<-Nf           #sigma=0.5
Nr2<-Nr          #sigma=0.5


#######Sigma = 0.95#######
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
sigma= 0.95      #recruitment facilitation- [Unstable for values much >0.95]
deltaUs= 0.2   #rate of kelp consumption by small urchins
deltaUm= 0.2   #rate  of kelp consumption by med. urchins
deltaUl=0.2    #rate of kelp consumption by large urchins
aM= 0.1      #conversion of kelp to urchins by med. urchins
aL= 0.1      #conversion of kelp to urchins by large urchins
b=0.1        #conversion of urchins to lobsters
gammaS = 0.15  #growth from small to med size class 
gammaM= 0.1  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
deltaLs=0.2 #  #rate of small urchin consumption by lobsters 
deltaLm=0.15 # #rate of med urchin consumption by lobsters  
deltaLl=0.05 # #rate of large urchin consumption by lobsters 
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

##Forward Path for FL with 3 values of Urchin recruitment facilitation
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
Nf <-array(NA,dim=c(length(FLseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLseq)){
  FL=FLseq[i]
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
matplot(FLseq,Nf[,,1], type="p", cex=0.5, pch=1, ylim=c(0,1500))
tail(out,5)

#####################
###Reverse path for Lobster mortality######
FLrevseq<-seq(1.0,0.0,length.out=100)  #vector of lobster fishing mortality rates to run model over
Nr <-array(NA,dim=c(length(FLrevseq),10,5))  #create vector with NA's to store equilibrium state var. population size
for(i in 1:length(FLrevseq)){
  FL=FLrevseq[i]
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
matpoints(FLrevseq,Nr[,,1],col="red", type="p",cex=0.5, pch=1)
tail(out,5)

Nf3<-Nf         #sigma=0.95
Nr3<-Nr          #sigma=0.95

########################################
#PLOTS
########################################
###Zoomed in, new colors, points on Hysteresis plots 
pdf("HysteresisFL3x_SigmaZoomv2.pdf")
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





##Plots: Forward path is darker line, reverse path is lighter shade###
#Kelp Plot
#pdf("HysteresisFU3xKelp.pdf")
matplot(FLseq,Nf1[,,1], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Kelp", xlab="Lobster fishing mortality", ylim=c(0,2500))
matlines(FLrevseq,Nr1[,,1], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FLseq,Nf2[,,1], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FLrevseq,Nr2[,,1], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FLseq,Nf3[,,1], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FLrevseq,Nr3[,,1], type="l", cex=0.5, lty=1, lwd=2, col="orange")
text(0.8, 700, labels="sigma=0.0, Blue", cex=0.75)
text(0.8,550, labels="sigma=0.5, Red", cex=0.75)
text(0.8,400, labels="sigma=0.95, Orange", cex=0.75)
#dev.off()

#Small Urchin Plot
pdf("HysteresisFLSig3xSmallUrch.pdf")
matplot(FLseq,Nf1[,,2], type="l", cex=0.5, lty=1, lwd=2, col=1, bty="L", las=1,
        ylab="Small Urchins", xlab="Lobster Fishing Mortality", ylim=c(0,40))
matlines(FLrevseq,Nr1[,,2], type="l", cex=0.5, lty=2, lwd=2, col=2)
matlines(FLseq,Nf2[,,2], type="l", cex=0.5, lty=1, lwd=2, col=3)
matlines(FLrevseq,Nr2[,,2], type="l", cex=0.5, lty=2, lwd=2, col=4)
matlines(FLseq,Nf3[,,2], type="l", cex=0.5, lty=1, lwd=2, col=5)
matlines(FLrevseq,Nr3[,,2], type="l", cex=0.5, lty=2, lwd=2, col=6)
text(0.1, 35, labels="FU=0.0, Green/Blue", cex=0.75)
text(0.1,32.5, labels="FU=0.1, Black/Red", cex=0.75)
text(0.1,30, labels="FU=0.25, Sky/Magenta", cex=0.75)
dev.off()

#Med. Urchin Plot
pdf("HysteresisFU3xMedUrchins.pdf")
matplot(sigmaseq,Nf1[,,3], type="l", cex=0.5, lty=1, lwd=2, col=1, bty="L", las=1, 
        ylab="Medium Urchins", xlab="Recruitment Facilitation", ylim=c(0,40))
matlines(sigmarevseq,Nr1[,,3], type="l", cex=0.5, lty=2, lwd=2, col=2)
matlines(sigmaseq,Nf2[,,3], type="l", cex=0.5, lty=1, lwd=2, col=3)
matlines(sigmarevseq,Nr2[,,3], type="l", cex=0.5, lty=2, lwd=2, col=4)
matlines(sigmaseq,Nf3[,,3], type="l", cex=0.5, lty=1, lwd=2, col=5)
matlines(sigmarevseq,Nr3[,,3], type="l", cex=0.5, lty=2, lwd=2, col=6)
text(0.9, 8, labels="FU=0.0, Green/Blue", cex=0.75)
text(0.9,7, labels="FU=0.1, Black/Red", cex=0.75)
text(0.9,6, labels="FU=0.25, Sky/Magenta", cex=0.75)
dev.off()

#Large Urchin Plot
pdf("HysteresisFU3xLargeUrchin.pdf")
matplot(sigmaseq,Nf1[,,4], type="l", cex=0.5, lty=1, lwd=2, col=1, bty="L", las=1,
        ylab="Large Urchins", xlab="Recruitment Facilitation", ylim=c(0,40))
matlines(sigmarevseq,Nr1[,,4], type="l", cex=0.5, lty=2, lwd=2, col=2)
matlines(sigmaseq,Nf2[,,4], type="l", cex=0.5, lty=1, lwd=2, col=3)
matlines(sigmarevseq,Nr2[,,4], type="l", cex=0.5, lty=2, lwd=2, col=4)
matlines(sigmaseq,Nf3[,,4], type="l", cex=0.5, lty=1, lwd=2, col=5)
matlines(sigmarevseq,Nr3[,,4], type="l", cex=0.5, lty=2, lwd=2, col=6)
text(0.1, 16, labels="FU=0.0, Green/Blue", cex=0.75)
text(0.1,15, labels="FU=0.1, Black/Red", cex=0.75)
text(0.1,14, labels="FU=0.25, Sky/Magenta", cex=0.75)
dev.off()

#Lobster Plot
pdf("HysteresisFU3xLobster.pdf")
matplot(sigmaseq,Nf2[,,5], type="l", cex=0.5, lty=1, lwd=2, col=1, bty="L", las=1,
        ylab="Lobsters", xlab="Recruitment Facilitation", ylim=c(0,16))
matlines(sigmarevseq,Nr2[,,5], type="l", cex=0.5, lty=2, lwd=2, col=2)
matlines(sigmaseq,Nf1[,,5], type="l", cex=0.5, lty=1, lwd=2, col=3)
matlines(sigmarevseq,Nr1[,,5], type="l", cex=0.5, lty=2, lwd=2, col=4)
matlines(sigmaseq,Nf3[,,5], type="l", cex=0.5, lty=1, lwd=2, col=5)
matlines(sigmarevseq,Nr3[,,5], type="l", cex=0.5, lty=2, lwd=2, col=6)
text(0.9, 10, labels="FU=0.0, Green/Blue", cex=0.75)
text(0.9,9, labels="FU=0.1, Black/Red", cex=0.75)
text(0.9,8, labels="FU=0.25, Sky/Magenta", cex=0.75)
dev.off()


###Hystersis plots 
pdf("HysteresisFL3x_SigmaKelp.pdf")
par(mfrow=c(3,2))
#Kelp
matplot(FLseq,Nf1[,,1], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Kelp (kg)", xlab="Lobster Fishing Mortality", ylim=c(0,2500))
matlines(FLrevseq,Nr1[,,1], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FLseq,Nf2[,,1], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FLrevseq,Nr2[,,1], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FLseq,Nf3[,,1], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FLrevseq,Nr3[,,1], type="l", cex=0.5, lty=1, lwd=2, col="orange")

#Small Urchin
matplot(FLseq,Nf1[,,2], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Small Urchins (kg)", xlab="Lobster Fishing Mortality", ylim=c(0,40))
matlines(FLrevseq,Nr1[,,2], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FLseq,Nf2[,,2], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FLrevseq,Nr2[,,2], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FLseq,Nf3[,,2], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FLrevseq,Nr3[,,2], type="l", cex=0.5, lty=1, lwd=2, col="orange")
#Med Urchin
matplot(FLseq,Nf1[,,3], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1, 
        ylab="Medium Urchins (kg)", xlab="Lobster Fishing Mortality", ylim=c(0,20))
matlines(FLrevseq,Nr1[,,3], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FLseq,Nf2[,,3], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FLrevseq,Nr2[,,3], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FLseq,Nf3[,,3], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FLrevseq,Nr3[,,3], type="l", cex=0.5, lty=1, lwd=2, col="orange")
#Large Urchin
matplot(FLseq,Nf1[,,4], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Large Urchins (kg)", xlab="Lobster Fishing Mortality", ylim=c(0,20))
matlines(FLrevseq,Nr1[,,4], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FLseq,Nf2[,,4], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FLrevseq,Nr2[,,4], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FLseq,Nf3[,,4], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FLrevseq,Nr3[,,4], type="l", cex=0.5, lty=1, lwd=2, col="orange")
#Lobster
matplot(FLseq,Nf1[,,5], type="l", cex=0.5, lty=1, lwd=2, col="midnightblue", bty="L", las=1,
        ylab="Lobsters (kg)", xlab="Lobster Fishing Mortality", ylim=c(0,20))
matlines(FLrevseq,Nr1[,,5], type="l", cex=0.5, lty=1, lwd=2, col="lightblue3")
matlines(FLseq,Nf2[,,5], type="l", cex=0.5, lty=1, lwd=2, col="darkred")
matlines(FLrevseq,Nr2[,,5], type="l", cex=0.5, lty=1, lwd=2, col="indianred1")
matlines(FLseq,Nf3[,,5], type="l", cex=0.5, lty=1, lwd=2, col="darkorange3")
matlines(FLrevseq,Nr3[,,5], type="l", cex=0.5, lty=1, lwd=2, col="orange")
text(0.8, 15.5, labels="sigma=0.0, Blue", cex=0.75)
text(0.8, 14, labels="sigma=0.5, Red", cex=0.75)
text(0.8, 12.5, labels="sigma=0.95, Orange", cex=0.75)
dev.off()


