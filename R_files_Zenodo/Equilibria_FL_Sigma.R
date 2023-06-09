##Robert Dunn
##Three stage urchin model,  Lobster F and Sigma

require(deSolve)
require(lattice)
rm(list=ls())
##Define system of ODEs, A=kelp, Us=small urchin, Ul=large urchin, L=spiny lobster
threestage=function(t, state, parms){
  with(as.list(c(state,parms)),{
    dA= (r*(1-A/Ka)-(deltaUs*Us+deltaUm*Um+deltaUl*Ul))*A        #kelp dynamics (logistic growth, type I func resp)
    
    dUs = ((aM*deltaUm*Um+aL*deltaUl*Ul)*A)*(1-sigma+sigma*(Ul/Kul))-(gammaS + (L*deltaLs/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+MUs)*Us #small urchin dynamics
    #^^reproduction from large urchins via kelp conversion, loss due to growth, type II func resp, natural mortality)
    
    dUm = (gammaS*Us)- (gammaM + (L*deltaLm/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul)) + MUm)*Um 
    
    dUl = (gammaM*Um)-((L*deltaLl/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+FU+MUl)*Ul   #large urchin dynamics
    #^^growth from small urchins, loss via type II func resp, fishing & natural mortality)
    
    dL = ((b*(deltaLs*Us+deltaLm*Um+deltaLl*Ul)/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))-ML-FL)*L          #spiny lobster dynamics     
    #^^conversion of lobsters to urchins via type II func resp, loss via fishing & natural mortality)
    return(list(c(dA, dUs, dUm, dUl, dL))) # return dn/dt as a list with each state variable as a column
  })  
}
# recruitment facilitation term:  *(1-sigma+sigma*(Ul/Kul)) #KUl estimated as ~=40 large urchins
#To estimate carrying capacity for large urchins, remove the recruitment facilitation term
# and run the model with no lobsters present. Equilibrium value of Ul is the new 
# value of Kul (large urchin carrying capacity).

###----------parameters----------------------###
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity- 
sigma= 0.5      #recruitment facilitation- 
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
MUl=0.1     #natural mortality of large urchins  (from Hilb.&Gut. stock assess.)
MUm= 0.1    #natural mortality of med. urchins  (from Hilb.&Gut. stock assess.)
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-07      #handling time of small urchins by lobsters
taum= 1e-08    #handling time of med urchins by lobsters
taul=1e-08    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run simulation to check that everything works
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))   
tail(out,10)  #check output

#basic timeseries plot again to check if everything works
plot(out$time,out$A, type="l", lwd=2, xlab="t", ylab="N", col="forestgreen", ylim=c(0,100),main="3 urchin stages")
lines(out$time, out$Us, type="l", lwd=2, col="purple")
lines(out$time,out$Um, type="l",lwd=2, col="blue")
lines(out$time,out$Ul, type="l", lwd=2, col="darkred")
lines(out$time,out$L, type="l", lwd=2, col="indianred2")

###############
## Lobster F at 6 values of recruitment facilitation, keep final equilibrium pop value 
###############

## Sigma=0.0 ##
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
N1 <-matrix(NA,length(FLseq),5)  #create vector with NA's to store equilibrium state var. population size
sigma=0.0                        #set level of recruitment facilitation
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N1[i,1]=out[200,2]   #kelp pop vector
  N1[i,2]=out[200,3]   #Us pop vector
  N1[i,3]=out[200,4]   #Um pop vector
  N1[i,4]=out[200,5]   #Ul pop vector
  N1[i,5]=out[200,6]   #lobster pop vector
}
#plot(FLseq,N1[,5], type="p", cex=0.5, pch=1)

## sigma=0.2 ##
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
N2 <-matrix(NA,length(FLseq),5)  #create vector with NA's to store equilibrium state var. population size
sigma=0.2
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N2[i,1]=out[200,2]   #kelp pop matrix
  N2[i,2]=out[200,3]   #Us pop matrix
  N2[i,3]=out[200,4]   #Um pop matrix
  N2[i,4]=out[200,5]   #Ul pop matrix
  N2[i,5]=out[200,6]   #lobster pop matrix
}

## sigma=0.4 ##
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
N3 <-matrix(NA,length(FLseq),5)  #create vector with NA's to store equilibrium state var. population size
sigma=0.4
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N3[i,1]=out[200,2]   #kelp pop matrix
  N3[i,2]=out[200,3]   #Us pop matrix
  N3[i,3]=out[200,4]   #Um pop matrix
  N3[i,4]=out[200,5]   #Ul pop matrix
  N3[i,5]=out[200,6]   #lobster pop matrix
}

## sigma=0.6 ##
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
N4 <-matrix(NA,length(FLseq),5)  #create vector with NA's to store equilibrium state var. population size
sigma=0.6
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N4[i,1]=out[200,2]   #kelp pop matrix
  N4[i,2]=out[200,3]   #Us pop matrix
  N4[i,3]=out[200,4]   #Um pop matrix
  N4[i,4]=out[200,5]   #Ul pop matrix
  N4[i,5]=out[200,6]   #lobster pop matrix
}

## sigma=0.8 ##
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
N5 <-matrix(NA,length(FLseq),5)  #create vector with NA's to store equilibrium state var. population size
sigma=0.8
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N5[i,1]=out[200,2]   #kelp pop matrix
  N5[i,2]=out[200,3]   #Us pop matrix
  N5[i,3]=out[200,4]   #Um pop matrix
  N5[i,4]=out[200,5]   #Ul pop matrix
  N5[i,5]=out[200,6]   #lobster pop matrix
}

## sigma=0.9 (setting sigma=1.0 crashes the function) ##
FLseq<-seq(0.0,1.0,length=100)   #vector of lobster fishing mortality rates to run model over
N6 <-matrix(NA,length(FLseq),5)  #create vector with NA's to store equilibrium state var. population size
sigma=0.9
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N6[i,1]=out[200,2]   #kelp pop matrix
  N6[i,2]=out[200,3]   #Us pop matrix
  N6[i,3]=out[200,4]   #Um pop matrix
  N6[i,4]=out[200,5]   #Ul pop matrix
  N6[i,5]=out[200,6]   #lobster pop matrix
}

############################
#plot all 6 lines (for each state variable separately)
############################

pdf("Equilibria_FL_sigma.pdf")
par(mfrow=c(3,2))
plot(FLseq,N1[,1], type="p", cex=0.5, pch="o", col="black", lwd=2, xlab="Lobster F",
     ylab="Kelp (kg)", bty="L", las=1)
legend(0.9,2450,xjust=0.5,cex=0.75,title="sigma",c("0.0","0.2","0.4","0.6","0.8","0.9"), bty="n", c("black", "blue4","purple", "red",
                                                                                              "orange", "yellow"))
lines(FLseq,N2[,1], type="p", cex=0.5, pch="o", col="blue4", lwd=2)
lines(FLseq,N3[,1], type="p", cex=0.5, pch="o", col="purple", lwd=2)
lines(FLseq,N4[,1], type="p", cex=0.5, pch="o", col="red", lwd=2)
lines(FLseq,N5[,1], type="p", cex=0.5, pch="o", col="orange", lwd=2)
lines(FLseq,N6[,1], type="p", cex=0.5, pch="o", col="yellow", lwd=2)

plot(FLseq,N1[,2], type="p", cex=0.5, pch="o", col="black", lwd=2, ylim=c(0,40),xlab="Lobster F",
     ylab="Small Urchins (kg)", bty="L", las=1)
lines(FLseq,N2[,2], type="p", cex=0.5, pch="o", col="blue4", lwd=2)
lines(FLseq,N3[,2], type="p", cex=0.5, pch="o", col="purple", lwd=2)
lines(FLseq,N4[,2], type="p", cex=0.5, pch="o", col="red", lwd=2)
lines(FLseq,N5[,2], type="p", cex=0.5, pch="o", col="orange", lwd=2)
lines(FLseq,N6[,2], type="p", cex=0.5, pch="o", col="yellow", lwd=2)

plot(FLseq,N1[,3], type="p", cex=0.5, pch="o", col="black", lwd=2, ylim=c(0,40),xlab="Lobster F",
     ylab="Med Urchins (kg)", bty="L", las=1)
lines(FLseq,N2[,3], type="p", cex=0.5, pch="o", col="blue4", lwd=2)
lines(FLseq,N3[,3], type="p", cex=0.5, pch="o", col="purple", lwd=2)
lines(FLseq,N4[,3], type="p", cex=0.5, pch="o", col="red", lwd=2)
lines(FLseq,N5[,3], type="p", cex=0.5, pch="o", col="orange", lwd=2)
lines(FLseq,N6[,3], type="p", cex=0.5, pch="o", col="yellow", lwd=2)

plot(FLseq,N1[,4], type="p", cex=0.5, pch="o", col="black", lwd=2, ylim=c(0,40),xlab="Lobster F",
     ylab="Large Urchins (kg)", bty="L", las=1)
lines(FLseq,N2[,4], type="p", cex=0.5, pch="o", col="blue4", lwd=2)
lines(FLseq,N3[,4], type="p", cex=0.5, pch="o", col="purple", lwd=2)
lines(FLseq,N4[,4], type="p", cex=0.5, pch="o", col="red", lwd=2)
lines(FLseq,N5[,4], type="p", cex=0.5, pch="o", col="orange", lwd=2)
lines(FLseq,N6[,4], type="p", cex=0.5, pch="o", col="yellow", lwd=2)

plot(FLseq,N1[,5], type="p", cex=0.5, pch="o", col="black", lwd=2, xlab="Lobster F", ylim=c(0,15),
     ylab="Lobsters (kg)", bty="L", las=1)
lines(FLseq,N2[,5], type="p", cex=0.5, pch="o", col="blue4", lwd=2)
lines(FLseq,N3[,5], type="p", cex=0.5, pch="o", col="purple", lwd=2)
lines(FLseq,N4[,5], type="p", cex=0.5, pch="o", col="red", lwd=2)
lines(FLseq,N5[,5], type="p", cex=0.5, pch="o", col="orange", lwd=2)
lines(FLseq,N6[,5], type="p", cex=0.5, pch="o", col="yellow", lwd=2)
dev.off()
####################
