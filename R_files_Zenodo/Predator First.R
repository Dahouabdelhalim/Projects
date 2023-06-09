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
state<-c(A=2459.864, Us=8.209526, Um=0.6845794, Ul=0.1081592, L=10.65874) #Unexploited equilibrium (no fishing), inits for trajectory sims
#state<-c(A=2722.933, Us=8.089573, Um=0.8283963, Ul=0.1588088, L=8.431772) #Unexploited equilibrium with NCEs
#State<-c(A=2617.935, Us=6.061423 , Um=.2810842, Ul=.0252874, L=20.2311)   #unexploited equilibrium with sigma=0
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
deltaLs=0.2  ##rate of small urchin consumption by lobsters 
deltaLm=0.15 ##rate of med urchin consumption by lobsters  
deltaLl=0.05 ##rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment   #0.75 for barrens inits
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run model to set up harvested equilibrium
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#Now, set up an MPA with 0 fishing mortality, use harvested equilibrium population values from above
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
FL=0.0    #close lobster fishery
#deltaUs= 0.1   #create NCEs
#deltaUm= 0.1   #when applicable

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestage, parms=parms))  
tail(out2,7)

#now close the urchin fishery
newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200],Ls = out2$Ls[200], L=out2$L[200])
FU=0.00         
parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
out3=as.data.frame(lsoda(y=newerstate, times=times, func=threestage, parms=parms))  
tail(out3,7)

y<-rbind(out2,out3)
y$time<-seq(1,400)
#time series of raw biomass
plot(y$time,y$A, ylim=c(0,3000), type="l", lwd=2, col="forestgreen")
plot(y$time, y$Us, ylim=c(0,30), type="l", lwd=2, col="purple")
lines(y$Um, type="l", lwd=2, col="blue")
lines(y$Ul, type="l", lwd=2, col="darkred")
lines(y$L, type="l", lwd=2, col="indianred2")


#Calculate community volatility##
#First, sum biomass of large urchins and large lobsters for the full time series
y$CommDens<-y$Us+y$Um+y$Ul+y$L
y$FishedDens<-y$Ul+y$L

maxfished<-max(y$FishedDens)  #13.157  With NCE: 11.190   #starting in barrens: 16.965
fishedvol<-(maxfished-y[400,8])/y[400,8] #0.222  With NCE: 0.302     #starting in barrens:0.575
print(fishedvol)
which.max(y$FishedDens) #4
PredNCE<-y$FishedDens
maxcomm<-max(y$CommDens) #   30.595
commvol<-(maxcomm - y[400,7])/y[400,7] #0.556, With NCE: 0.747
print(commvol)
#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria

#Fished community biomass plot
plot(FishedDens~time, data=y, type="l", bty="L", las=1, ylim=c(7,15),lwd=2, ylab="Fished community biomass (kg)", xlab="Time (y)")

###################################
####Urchin barren fishing levels###
###################################
#Bi-stable fishing rates: fu=0.15, fl=0.65
#Barren fishing rates: fu=0.25, fl=0.8
###----------parameters----------------------###
state<-c(A=2459.864, Us=8.209526, Um=0.6845794, Ul=0.1081592, L=10.65874) #Unexploited equilibrium (no fishing), inits for trajectory sims
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
deltaLs=0.2  ##rate of small urchin consumption by lobsters 
deltaLm=0.15 ##rate of med urchin consumption by lobsters  
deltaLl=0.05 ##rate of large urchin consumption by lobsters
MUl=0.1     #natural mortality of large urchins  (0.003 from Hilb.&Gut. stock assess., but too small in the GSA, adjusted it higher)
MUm= 0.1    #natural mortality of med. urchins  (0.05 from Hilb.&Gut. stock assess. but too small in the GSA, adjusted it higher))
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.25     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= .8     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run model to set up harvested equilibrium
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#Now, set up an MPA with 0 fishing mortality, use harvested equilibrium population values from above
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
FL=0.0    #close lobster fishery
#deltaUs= 0.1   #create NCEs
#deltaUm= 0.1   

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestage, parms=parms))  
tail(out2,7)

#now close the urchin fishery
newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200],Ls = out2$Ls[200], L=out2$L[200])
FU=0.00         
parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
out3=as.data.frame(lsoda(y=newerstate, times=times, func=threestage, parms=parms))  
tail(out3,7)

y<-rbind(out2,out3)
y$time<-seq(1,400)
#time series of raw biomass
plot(y$time,y$A, ylim=c(0,3000), type="l", lwd=2, col="forestgreen")
plot(y$time, y$Us, ylim=c(0,30), type="l", lwd=2, col="purple")
lines(y$Um, type="l", lwd=2, col="blue")
lines(y$Ul, type="l", lwd=2, col="darkred")
lines(y$L, type="l", lwd=2, col="indianred2")


#Calculate community volatility##
#First, sum biomass of large urchins and large lobsters for the full time series
y$CommDens<-y$Us+y$Um+y$Ul+y$L
y$FishedDens<-y$Ul+y$L
which.max(y$FishedDens)   #5 for bi-stable, 47 for barren
maxfished<-max(y$FishedDens)  #13.157  #in barrens: 15.604  #bistable: 15.34
fishedvol<-(maxfished-y[400,8])/y[400,8] #0.222  #in bistable: 0.425  #in barrens: 0.449
maxcomm<-max(y$CommDens)
commvol<-(maxcomm-y[400,7])/y[400,7]

#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria
