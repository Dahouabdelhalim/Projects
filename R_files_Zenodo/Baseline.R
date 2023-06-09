##Robert Dunn
## Three stage urchin baseline model


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

#To estimate carrying capacity for large urchins, remove the recruitment facilitation term
# and run the model with no lobsters present. Equilibrium value of Ul is the new 
# value of Kul (large urchin carrying capacity).
## Estimated large urchin carrying capacity by setting recruitment facilitation (sigma)
## =0 and with no lobsters present. KUl is ~40 large urchins.
##################################################
###----------parameters----------------------###
state<-c(A=1000, Us=70, Um=70, Ul=70, L=20)   #Initial state of the system- kelp forest
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC? [Unstable for values much >0.95]
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

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run model as a check
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#basic plot check the timeseries
plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen",
     ylim=c(0,500),main="3 urchin stages")
lines(out$time, out$Us, type="l", col="purple")
lines(out$time,out$Um, type="l", col="blue")
lines(out$time,out$Ul, type="l", col="darkred")
lines(out$time,out$L, type="l", col="indianred2")

#plot the time series with 2 axes: left axis for kelp, right for other state variables
pdf("BaselineKelp.pdf")
plot(out$time,out$A, type="l", xlab="t (yr)", lwd=2, col="forestgreen", 
     main="", axes=F, bty="L", ylab="")
legend(x=155,y=2350,cex=0.75, c("Kelp","Small Urchin","Medium Urchin","Large Urchin","Lobster"), 
       bty="n", c("forestgreen","purple", "blue", "darkred","indianred2"))
axis(2, ylim=c(0,2050), col="black", las=1)
#mtext(2, line=2)
par(new=T)
plot(out$time,out$Us, type="l", col="purple", lwd=2, axes=F, ylab="", xlab="")
axis(4, ylim=c(0,100), col="black", las=1)
mtext(4, text="N", line=2)
lines(out$time,out$Ul, type="l", col="darkred",ylab="", xlab="",lwd=2)
lines(out$time, out$L, type="l", col="indianred2",ylab="",xlab="",lwd=2)
lines(out$time,out$Um, type="l", col="blue", ylab="", xlab="", lwd=2)
axis(1,xlim=c(0,60), las=1)
dev.off()


####
###Now start in Urchin Barren state
####

state<-c(A=1, Us=1e2, Um=1e2, Ul=1e2, L=5)  #Alternate barrens state
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC? [Unstable for values much >0.95]
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

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#basic plot
plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen",
     ylim=c(0,500),main="3 urchin stages")
lines(out$time, out$Us, type="l", col="purple")
lines(out$time,out$Um, type="l", col="blue")
lines(out$time,out$Ul, type="l", col="darkred")
lines(out$time,out$L, type="l", col="indianred2")

#plot the time series with 2 axes: left axis for kelp, right for other state variables
pdf("BaselineBarren.pdf")
plot(out$time,out$A, type="l", xlab="t (yr)", lwd=2, col="forestgreen", 
     main="", axes=F, bty="L", ylab="")
legend(x=155,y=2350,cex=0.75, c("Kelp","Small Urchin","Medium Urchin","Large Urchin","Lobster"), 
       bty="n", c("forestgreen","purple", "blue", "darkred","indianred2"))
axis(2, ylim=c(0,2500), col="black", las=1)
#mtext(2, line=2)
par(new=T)
plot(out$time,out$Us, type="l", col="purple", lwd=2, axes=F, ylab="", xlab="")
axis(4, ylim=c(0,70), col="black", las=1)
mtext(4, text="N", line=2)
lines(out$time,out$Ul, type="l", col="darkred",ylab="", xlab="",lwd=2)
lines(out$time, out$L, type="l", col="indianred2",ylab="",xlab="",lwd=2)
lines(out$time,out$Um, type="l", col="blue", ylab="", xlab="", lwd=2)
axis(1,xlim=c(0,60), las=1)
dev.off()
