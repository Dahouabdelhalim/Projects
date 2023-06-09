##Robert Dunn
##Calculate MSY and F_msy for three stage urchin baseline model

require(deSolve)
require(lattice)
rm(list=ls())
setwd("~/PhD/Community Dyn. Model/Plots") #my working directory for this paper

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
    
    dFu = FU*Ul                #urchin fishery harvest biomass
    
    dFl = FL*L                 #lobster fishery harvest biomass
    return(list(c(dA, dUs, dUm, dUl, dL, dFu, dFl))) # return dn/dt as a list with each state variable as a column
  })  
}
# recruitment facilitation term:   *(1-sigma+sigma*(Ul/Kul))
#To estimate carrying capacity for large urchins, remove the recruitment facilitation term
# and run the model with no lobsters present. Equilibrium value of Ul is the new 
# value of Kul (large urchin carrying capacity).
###----------parameters----------------------###

r= 10               #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #guess at large urchin carrying capacity- ask Marissa
sigma= 0.5      #recruitment facilitation- how much is a first guess? see paper from BC?
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
MUl=0.1      #0.003     #natural mortality of large urchins  (from Hilb.&Gut. stock assess.)
MUm=0.1     # 0.05    #natural mortality of med. urchins  (from Hilb.&Gut. stock assess.)
MUs=0.1    #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
FU=0.1     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
taus= 1e-08      #handling time of small urchins by lobsters
taum= 1e-07    #handling time of med urchins by lobsters
taul=1e-06    #handling time of large urchins by lobsters
ML= 0.175      # natural mortality of lobsters: 0.17 from stock assessment
FL= 0.25     # fishing mortality of lobsters: Fmsy=0.25 from stock assessment

parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
        deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
        deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
        MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) #aggregate parameters into 1 vector to pass to lsoda
tf= 200   #run time
times<-1:tf     #times vector

state<-c(A=1e3, Us=70, Um=70, Ul=70, L=20, Fu=0, Fl=0) #Initial state of the system, kelp forest
#Run baseline simulation as check
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))   #numerically integrate the logistic function
tail(out,10)  #check output

#Time Series Plot

plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen",
     ylim=c(0,50),main="3 urchin stages")
lines(out$time, out$Us, type="l", col="purple")
lines(out$time,out$Um, type="l", col="blue")
lines(out$time,out$Ul, type="l", col="darkred")
lines(out$time,out$L, type="l", col="indianred2")


#harvest plots
plot(out$time, out$Fu, type="l", xlab="t", ylab="Harvest", col="purple")
plot(out$time, out$Fl, type="l", col="indianred2")

#####Calculate MSY for Lobsters when FU=0.1######
options(digits=9)
tf=200
times=1:tf
FLseq<-seq(0,1.0, length=101)
N<-array(NA, dim=c(length(FLseq), 10, 7))
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  state<-c(A=1e3, Us=70, Um=70, Ul=70, L=20, Fu=0, Fl=0)
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N[i,,1]=out[191:200,2]   #kelp pop matrix
  N[i,,2]=out[191:200,3]   #Us pop matrix
  N[i,,3]=out[191:200,4]   #Um pop matrix
  N[i,,4]=out[191:200,5]   #Ul pop matrix
  N[i,,5]=out[191:200,6]   #lobster pop matrix
  N[i,,6]=out[191:200,7]
  N[i,,7]=out[191:200,8]
}
N[1:101,,7]
N[12,7,7]-N[12,6,7]   #what are these?

Catch_0.1<-rep(NA,101)
Catch_0.1<-N[1:101,10,7]-N[1:101,9,7] #calculates annual catch
plot(Catch_0.1 ~ FLseq)  #plot confirms catch is maximized at FL= 0.5 when FU=0.1
Catch_0.1
max(Catch_0.1)

plot(Catch_0.1 ~ FLseq, xlab="Lobster fishery harvest", ylab="Equilibrium biomass yield", 
     type="p", cex=1.0, pch="o", col="black", lwd=2, bty="L", las=1)

###############
#####Calculate MSY for Lobsters when FU=0.0######
FU=0.0
N<-array(NA, dim=c(length(FLseq), 10, 7))
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  state<-c(A=1e3, Us=70, Um=70, Ul=70, L=20, Fu=0, Fl=0)
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N[i,,1]=out[191:200,2]   #kelp pop matrix
  N[i,,2]=out[191:200,3]   #Us pop matrix
  N[i,,3]=out[191:200,4]   #Um pop matrix
  N[i,,4]=out[191:200,5]   #Ul pop matrix
  N[i,,5]=out[191:200,6]   #lobster pop matrix
  N[i,,6]=out[191:200,7]
  N[i,,7]=out[191:200,8]
}

Catch_0.0<-rep(NA,101)
Catch_0.0<-N[1:101,10,7]-N[1:101,9,7] #calculates annual catch
max(Catch_0.0)
Catch_0.0

plot(Catch_0.0 ~ FLseq)  #plot confirms catch is maximized at FL=0.5 when FU=0.0
points(Catch_0.1 ~ FLseq, col="red")

####################
#####Calculate MSY for Lobsters when FU=0.25######
FU=0.25
N<-array(NA, dim=c(length(FLseq), 10, 7))
for(i in 1:length(FLseq)){
  FL=FLseq[i]
  state<-c(A=1e3, Us=70, Um=70, Ul=70, L=20, Fu=0, Fl=0)
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N[i,,1]=out[191:200,2]   #kelp pop matrix
  N[i,,2]=out[191:200,3]   #Us pop matrix
  N[i,,3]=out[191:200,4]   #Um pop matrix
  N[i,,4]=out[191:200,5]   #Ul pop matrix
  N[i,,5]=out[191:200,6]   #lobster pop matrix
  N[i,,6]=out[191:200,7]
  N[i,,7]=out[191:200,8]
}

Catch_0.25<-rep(NA,101)
Catch_0.25<-N[1:101,10,7]-N[1:101,9,7] #calculates annual catch
max(Catch_0.25)
Catch_0.25
plot(Catch_0.25 ~ FLseq)  #plot confirms catch is maximized at FL=0.5 when FU=0.3
points(Catch_0.1 ~ FLseq, col="blue")
points(Catch_0.0 ~ FLseq, col="red")

###Big plot with all 3####
pdf("LobsterYield_All3.pdf")
plot(Catch_0.0 ~ FLseq, xlab="Lobster fishery harvest", ylab="Equilibrium lobster biomass yield (kg)", 
 type="l", cex=1.0, pch="o", col="blue", lwd=3, bty="L", las=1)
lines(Catch_0.1 ~ FLseq, col="red", lwd=3)
lines(Catch_0.25 ~ FLseq, col="black", lwd=3)
dev.off()