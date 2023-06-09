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

###SYNCHRONOUS RECOVERY with EXOGENOUS RECRUITMENT#######
###And variable size structured predation by lobsters ##
##################################################
###----------parameters----------------------###
#state<-c(A=2496.937, Us=7.848038, Um=0.512968, Ul=0.023371, L=15.75106)   #Unexploited equilibrium (no size structure)
#state<-c(A=2753.96, Us=3.625995, Um=0.3046916, Ul=0.1699775, L=15.85082)   #Unexploited equilibrium (extra structure)
state<-c(A=2326.84, Us=10.56738, Um=0.5887021, Ul=0.06324361, L=16.61697)  #Unexploited equilibirum (hump shaped)
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
deltaLs=0.1 ##rate of small urchin consumption by lobsters .2 corresponds to 12mm urchin, 0.133=no-size struct, 0.1 hump
deltaLm=0.15 #rate of med urchin consumption by lobsters  .0071 correponds to 40mm urchin   0.133=no-size struct  0.15 hump
deltaLl=0.05 #rate of large urchin consumption by lobsters .0009 correponds to 82.5mm urchin  0.133=no-size struct  0.05 hump
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

#Set up exploited equilibrium
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#basic plot check the timeseries
plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen",
     ylim=c(0,40),xlim=c(0,50), main="3 urchin stages")
lines(out$time, out$Us, type="l", col="purple")
lines(out$time,out$Um, type="l", col="blue")
lines(out$time,out$Ul, type="l", col="darkred")
lines(out$time,out$L, type="l", col="indianred2")

#Now, set up an MPA with 0 fishing mortality, use equilibrium population values from above and look at trajectories.
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200],Ls = out$Ls[200], L=out$L[200])
FU=0.00         #close urchin fishery
FL=0.0    #close lobster fishery

parms=c(r=r, Ka=Ka, Kul=Kul, phiU=phiU, phiL=phiL, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestage, parms=parms))  #run again, no fishing, new inits
tail(out2,7)

y<-out2
y$time<-seq(1,200)

#Calculate community volatility##
#First, sum biomass of urchins and lobsters for the full time series
y$CommDens<-y$Us+y$Um+y$Ul+y$L
y$FishedDens<-y$Ul+y$L

maxfished<-max(y$FishedDens)  # 
fishedvol<-(maxfished-y[200,8])/y[200,8] #.2037   hump: 0.075
print(fishedvol)
maxcomm<-max(y$CommDens) #
commvol<-(maxcomm - y[200,7])/y[200,7] #0.273    hump: 0.5039
print(commvol)
#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria

####################################################################
#####Predator First Recovery with Exogenous Recruitment##########
################################################################
#state<-c(A=2496.937, Us=7.848038, Um=0.512968, Ul=0.023371, L=15.75106)   #Unexploited equilibrium (no size structure)
#state<-c(A=2753.96, Us=3.625995, Um=0.3046916, Ul=0.1699775, L=15.85082)   #Unexploited equilibrium (extra structure)
state<-c(A=2326.84, Us=10.56738, Um=0.5887021, Ul=0.06324361, L=16.61697)  #Unexploited equilibirum (hump shaped)
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
deltaLs=0.1 ##rate of small urchin consumption by lobsters .2 corresponds to 12mm urchin, 0.133=no-size struct, 0.1 hump
deltaLm=0.15 #rate of med urchin consumption by lobsters  .0071 correponds to 40mm urchin   0.133=no-size struct  0.15 hump
deltaLl=0.05 #rate of large urchin consumption by lobsters .0009 correponds to 82.5mm urchin  0.133=no-size struct  0.05 hump
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

#Set up exploited equilibrium
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

maxfished<-max(y$FishedDens)  # 
fishedvol<-(maxfished-y[400,8])/y[400,8] #.1517   hump: 0.0677
print(fishedvol)
maxcomm<-max(y$CommDens) #
commvol<-(maxcomm - y[400,7])/y[400,7] #0.2615     hump: 0.5039
print(commvol)
#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria


#######################################################
######Prey first recovery with exogenous recruitment###
########################################################
#state<-c(A=2496.937, Us=7.848038, Um=0.512968, Ul=0.023371, L=15.75106)   #Unexploited equilibrium (no size structure)
#state<-c(A=2753.96, Us=3.625995, Um=0.3046916, Ul=0.1699775, L=15.85082)   #Unexploited equilibrium (extra structure)
state<-c(A=2326.84, Us=10.56738, Um=0.5887021, Ul=0.06324361, L=16.61697)  #Unexploited equilibirum (hump shaped)
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
deltaLs=0.1 ##rate of small urchin consumption by lobsters .2 corresponds to 12mm urchin, 0.133=no-size struct, 0.1 hump
deltaLm=0.15 #rate of med urchin consumption by lobsters  .0071 correponds to 40mm urchin   0.133=no-size struct  0.15 hump
deltaLl=0.05 #rate of large urchin consumption by lobsters .0009 correponds to 82.5mm urchin  0.133=no-size struct  0.05 hump
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

#Set up exploited equilibrium
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
tail(out,7)

#Now, set up an MPA with 0 fishing mortality, use equilibrium population values from above
# and look at trajectories.
newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200], L=out$L[200])
FU=0.0    #close urchin fishery

parms=c(r=r, Ka=Ka, Kul=Kul, phiU=phiU, phiL=phiL, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
out2=as.data.frame(lsoda(y=newstate, times=times, func=threestage, parms=parms))  #run again, no fishing, new inits
tail(out2,7)

#now stop fishing for second species
newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200], L=out2$L[200])
FL=0.00         #close lobster fishery

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

maxfished<-max(y$FishedDens)  # 
fishedvol<-(maxfished-y[400,8])/y[400,8] #0.0902, extra: 0.282, hump:0.077
print(fishedvol)
maxcomm<-max(y$CommDens) #
commvol<-(maxcomm - y[400,7])/y[400,7] #.4849, extra: 0.337, hump: 0.513
print(commvol)
#RETURN TIMES###
#Calculate manually by looking at y dataframe, add up duration of first and 2nd transient periods until within 10% of equilibria
