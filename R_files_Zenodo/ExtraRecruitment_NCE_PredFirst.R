##Robert Dunn
## Three stage urchin baseline model
##With an extra recruitment term for small urchins and lobsters (phiU or phiL)


require(deSolve)
require(lattice)
rm(list=ls())

##Define system of ODEs, A=kelp, Us=small urchin, Ul=large urchin, L=spiny lobster
threestage=function(t, state, parms){
  with(as.list(c(state,parms)),{
    dA= (r*(1-A/Ka)-(deltaUs*Us+deltaUm*Um+deltaUl*Ul))*A        #kelp dynamics (logistic growth, type I func resp)
    
    dUs = (((aM*deltaUm*Um+aL*deltaUl*Ul)*A) +phiU)*(1-sigma+sigma*(Ul/Kul))-(gammaS + (L*deltaLs/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+MUs)*Us #small urchin dynamics
    #^^reproduction from large urchins via kelp conversion, loss due to growth, type II func resp, natural mortality), e is external recruitment
    
    dUm = (gammaS*Us)- (gammaM + (L*deltaLm/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul)) + MUm)*Um 
    
    dUl = (gammaM*Um)-((L*deltaLl/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))+FU+MUl)*Ul   #large urchin dynamics
    #^^growth from small urchins, loss via type II func resp, fishing & natural mortality)
    
    dL = phiL + ((b*(deltaLs*Us+deltaLm*Um+deltaLl*Ul)/(1+taus*deltaLs*Us+taum*deltaLm*Um+taul*deltaLl*Ul))-ML-FL)*L          #spiny lobster dynamics     
    #^^conversion of urchins to lobsters via type II func resp, loss via fishing & natural mortality)
    return(list(c(dA, dUs, dUm, dUl, dL))) # return dn/dt as a list with each state variable as a column
  })  
}

##################################################
###----------parameters----------------------###
state<-c(A=2685.046 , Us=4.886611, Um=0.3219738, Ul=0.04064373, L=13.84371)   #Unexploited equilibrium
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
phiU=10            #Urchin external recruitment constant- how much to include? See Ebert et al. 1994, Okamoto dissertation 2014
phiL=1          #Lobster external recruitment constant
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

parms=c(r=r, Ka=Ka, Kul=Kul, phiU=phiU, phiL=phiL, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Set up Exploited equilibrium
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#Now, set up an MPA with 0 fishing mortality, use equilibrium population values from above
# and look at trajectories.
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200], L=out$L[200])
FL=0.0    #close lobster fishery
parms=c(r=r, Ka=Ka, Kul=Kul, phiU=phiU, phiL=phiL, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
tnce = 10
ncetimes<-1:tnce
out2=as.data.frame(lsoda(y=newstate, times=ncetimes, func=threestage, parms=parms))  

#Time-delayed NCEs 
ncestate<-c(A=out2$A[10],Us=out2$Us[10],Um=out2$Um[10],Ul=out2$Ul[10], L=out2$L[10])
deltaUs= 0.1   #create NCEs
deltaUm= 0.1 
parms= c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
         deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
         deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
         MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)

out3=as.data.frame(lsoda(y=ncestate, times=times, func=threestage, parms=parms))  

#now close the urchin fishery
newerstate<-c(A=out3$A[200],Us=out3$Us[200],Um=out3$Um[200],Ul=out3$Ul[200], L=out3$L[200])
FU=0.00         #close urchin fishery
parms=c(r=r, Ka=Ka, Kul=Kul, phiU=phiU, phiL=phiL, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
out4=as.data.frame(lsoda(y=newerstate, times=times, func=threestage, parms=parms))

y<-rbind(out2,out3, out4)
y$time<-seq(1,410)

#Calculate community volatility##
#First, sum biomass of urchins and lobsters for the full time series
y$CommDens<-y$Us+y$Um+y$Ul+y$L
y$FishedDens<-y$Ul+y$L

maxfished<-max(y$FishedDens)  # 15.425  #15.42 with delayed NCE
fishedvol<-(maxfished-y[410,8])/y[410,8] #0.111   #with delayed NCE: 0.309
maxcomm<-max(y$CommDens) #27.07
commvol<-(maxcomm - y[410,7])/y[410,7] #0.418   #with delayed NCE: 0.657

#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria

plot(FishedDens ~ time, data=y, lwd=2, type="l")
################################
###Overharvested initial state##
################################

state<-c(A=2685.046 , Us=4.886611, Um=0.3219738, Ul=0.04064373, L=13.84371)   #Unexploited equilibrium
r= 10                #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        # large urchin carrying capacity- 
phiU=10            #Urchin external recruitment constant- how much to include? See Ebert et al. 1994, Okamoto dissertation 2014
phiL=1          #Lobster external recruitment constant
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
FU=0.25     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.8     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment
diff=1e-3  #difference between previous equilibrium N and new initial conditions

parms=c(r=r, Ka=Ka, Kul=Kul, phiU=phiU, phiL=phiL, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

#Run model as a check
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#Now, set up an MPA with 0 fishing mortality, use equilibrium population values from above
# and look at trajectories.
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200], L=out$L[200])
FL=0.0    #close lobster fishery

parms=c(r=r, Ka=Ka, Kul=Kul, phiU=phiU, phiL=phiL, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestage, parms=parms))  #run again, no fishing, new inits
tail(out2,7)

#now stop fishing for second species
newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200], L=out2$L[200])
FU=0.00         #close urchin fishery
parms=c(r=r, Ka=Ka, Kul=Kul, phiU=phiU, phiL=phiL, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
out3=as.data.frame(lsoda(y=newerstate, times=times, func=threestage, parms=parms))
tail(out3,7)

y<-rbind(out2,out3)
y$time<-seq(1,400)

#Calculate community volatility##
#First, sum biomass of urchins and lobsters for the full time series
y$CommDens<-y$Us+y$Um+y$Ul+y$L
y$FishedDens<-y$Ul+y$L

maxfished<-max(y$FishedDens)  #     #overharvested: 17.389
fishedvol<-(maxfished-y[400,8])/y[400,8] #    #overharvested: 0.252
maxcomm<-max(y$CommDens) #
commvol<-(maxcomm - y[400,7])/y[400,7] #   #overharvested:1.383

plot(CommDens~time, data=y)

#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria
