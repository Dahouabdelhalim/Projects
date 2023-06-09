##3 stage urchin community dynamics model
##Robert Dunn

##GSA using Community Volatility as the model outcome: Predator First RECOVERY

##Create range of values of parameters from which to sample (uniform dist. randomly sampled)
##Make 4k sample combinations, store them, then run simulations for each
##Store the community volatility for each simulation
##See the random Forest vignette for specific questions about coding, functions

rm(list=ls())
require(deSolve)      #ODE solver
require(randomForest) #GSA package
require(rpart)        #Classification/regression tree package
require(rpart.plot)   #makes nicer trees

##Define system of ODEs, A=kelp, Us=small urchin, Ul=large urchin, L=spiny lobster
threestage=function(t, state, parms){
  with(as.list(c(state,parms)),{
    #print(t)
    #print(state)
    #state = state*(state>=0) # no negative population sizes allowed
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

##Randomly sample from uniform distributions for a range of each parameter value
##store in the parameters matrix
sims=4e3
parameters <- data.frame(rbind(runif(sims,6,12),    #r
                               runif(sims,1e3,6e3), #Ka
                               runif(sims,10,70),   #Kul
                               runif(sims,0.0,0.75), # sigma #originally was 0.975 but that gave some odd volatility values.
                               runif(sims,0.15,0.25), #deltaUs
                               runif(sims,0.15,0.25), #deltaUm
                               runif(sims,0.15,0.25), #deltaUl
                               runif(sims,0.05,0.15),  #aM
                               runif(sims,0.05,0.15),  #aL
                               runif(sims,0.05,0.15),  #b
                               runif(sims,0.05,0.25), #gammaS
                               runif(sims,0.05,0.25), #gammaM
                               runif(sims,0.01,0.3),  #deltaLs
                               runif(sims,0.01,0.3),  #deltaLm
                               runif(sims,0.01,0.2),  #deltaLl
                               runif(sims,0.1,0.2),   #Mus
                               runif(sims,0.1,0.2),   #Mum
                               runif(sims,0.1,0.2),   #Mul
                               runif(sims,0.0,0.4),   #FU
                               runif(sims,1e-10,1e-6), #tauS
                               runif(sims,1e-10,1e-6), #tauM
                               runif(sims,1e-10,1e-6), #tauL
                               runif(sims,0.1,0.4),    #ML
                               runif(sims,0,0.7)))     #FL

##Problems?   ##try RK4 solver    ##try reducing the step size/tolerance

###----Run a simulation with each set of parameters chosen above-------###
tf= 200          #run time
times<-1:tf      #times vector

N <-matrix(NA,1,ncol(parameters))  #create vector with NA's to store equilibrium state var. population size
state<-c(A=1000, Us=2, Um=2, Ul=2, L=5)     #Kelp forest inits
#state<-c(A=25, Us=25, Um=19, Ul=5, L=0.5)  #Barren Inits

for(j in 1:ncol(parameters)){
  r= parameters[1,j]                #growth rate of kelp
  Ka=parameters[2,j]        #carrying capacity of kelp
  Kul=parameters[3,j]         # large urchin carrying capacity- 
  sigma= parameters[4,j]       #recruitment facilitation- 
  deltaUs= parameters[5,j]    #rate of kelp consumption by small urchins
  deltaUm= parameters[6,j]    #rate  of kelp consumption by med. urchins
  deltaUl=parameters[7,j]     #rate of kelp consumption by large urchins
  aM=parameters[8,j]       #conversion of kelp to urchins by med. urchins
  aL= parameters[9,j]      #conversion of kelp to urchins by large urchins
  b=parameters[10,j]         #conversion of urchins to lobsters
  gammaS =parameters[11,j]   #growth from small to med size class 
  gammaM= parameters[12,j]  #growth of med to big size class  (maturation rate from Hilb&Gut stock assess.)
  deltaLs=parameters[13,j]   #rate of small urchin consumption by lobsters
  deltaLm=parameters[14,j]   #rate of med urchin consumption by lobsters
  deltaLl=parameters[15,j]    #rate of large urchin consumption by lobsters
  MUl=parameters[16,j]      #natural mortality of large urchins  (from Hilb.&Gut. stock assess.)
  MUm= parameters[17,j]    #natural mortality of med. urchins  (from Hilb.&Gut. stock assess.)
  MUs=parameters[18,j]     #natural mortality of small urchins (from Hilb.&Gut. stock assess.)
  FU=parameters[19,j]     # fishing mortality of urchins   ( Hilb.&Gut. stock assess. set F~0.001, but Loo disagrees)
  taus= parameters[20,j]       #handling time of small urchins by lobsters
  taum= parameters[21,j]     #handling time of med urchins by lobsters
  taul=parameters[22,j]     #handling time of large urchins by lobsters
  ML= parameters[23,j]      # natural mortality of lobsters: 0.17 from stock assessment
  FL= parameters[24,j]       #fishing mortality of urchins
  
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, #aggregate parameters
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM,   #for input into lsoda
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)  
  
  #Run model to set up harvested equilibrium
  #out=lsoda(y=state, times=times, func=threestage, parms=parms) #, method="lsode",maxsteps=1e4, verbose=TRUE)
  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  

  #Now, set up an MPA with 0 fishing mortality, use harvested equilibrium population values from above
  newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200], L=out$L[200])
  FL=0.0    #close lobster fishery
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
  out2=as.data.frame(lsoda(y=newstate, times=times, func=threestage, parms=parms))  
  
  #now stop fishing for urchins
  newerstate<-c(A=out2$A[200],Us=out2$Us[200],Um=out2$Um[200],Ul=out2$Ul[200], L=out2$L[200])
  FU=0.00         
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
  out3=as.data.frame(lsoda(y=newerstate, times=times, func=threestage, parms=parms))
  
  y<-rbind(out2,out3)
  
  #Calculate community volatility##
  #sum biomass of large urchins and large lobsters for the full time series
  FishedDens<-y$Ul+y$L
  
  maxfished<-max(FishedDens)  #calculate max value of aggregate fished community biomass
  fishedvol<-(maxfished-FishedDens[400])/FishedDens[400]  #calculate fished community volatility
  
  N[1,j]=fishedvol   #Volatility value for a given simulation
}

##Set up data for GSA
predGSA<-matrix(NA,sims,25)   #dummy matrix to hold state variables and sampled parameter values
predGSA[,1]<-t(N)          #input state variable values from simulations
predGSA[,2:25]<-t(parameters) #input randomly chosen parameter values
predGSA<-as.data.frame(predGSA)
colnames(predGSA)<-c("Volatility", "r","Ka","Kul","sigma",    #define parameter labels
                     "deltaUs","deltaUm", "deltaUl",
                     "aM","aL", "b", "gammaS", "gammaM", "deltaLs", "deltaLm",
                     "deltaLl","MUs","MUm","MUl", "FU","taus","taum","taul",
                     "ML", "FL")

#Run Predator First GSA- random forest and then CART
##Random forest on volatility
predvol.rf<-randomForest(predGSA[,1] ~ .,
                     data       = predGSA[,2:25],
                     importance = TRUE)
print(predvol.rf)
round(importance(predvol.rf),2)
predvol.imp = pmax (importance (predvol.rf,scale = FALSE) [, 1], 0)
predvol.imp = predvol.imp/sum(predvol.imp)

#Bar plot of normalized importance values
#pdf("ImportancePred.pdf")
par(mai=c(1.02,1.02,0.62,0.42))
barplot (predvol.imp,
         horiz = T,
         main  = 'Synchronous volatility',
         xlim  = c (0, .5),  las=1, col="springgreen3")
#dev.off()

#hi resolution plot
#tiff("testPredFirst.tif", width=3, height=3, units="in", res=600, pointsize=4)
par(mar=c(5,5,4,4), xaxs="i", yaxs="i", cex.axis=1.3, cex.lab=1.3)
barplot (predvol.imp,
         horiz = T,
         main  = 'Synchronous volatility',
         xlim  = c (0, .5),  las=1, xlab="Importance value")
#dev.off()

#Sorted bar plot of normalized importance values in base graphics 
Pred_importance<-as.data.frame(predvol.imp)
predvol.imp<-sort(predvol.imp, decreasing=F)
barplot(predvol.imp, horiz=T, xlim=c(0,0.5), las=1, xlab="Importance value") #sorted plot without greek symbols
#updated plot with sorted labels
tiff("PredSorted_HiRes.tif", width=6, height=6, units="in", res=800, pointsize=6, compression="lzw")
par(mai=c(1.32,1.42,0.42,0.62), cex.axis=1.5, cex.lab=1.5)
barplot(predvol.imp,
        names.arg=c(expression(italic(K[UL])),
                    expression(italic(M[Um])), 
                    expression(paste(tau[italic(S)])),
                    expression(paste(tau[italic(L)])),
                    expression(paste(alpha[italic(L)])),
                    expression(italic(M[UL])),
                    expression(paste(gamma[italic(M)])),
                    expression(paste(delta[italic(UL)])), 
                    expression(paste(delta[italic(Um)])),
                    expression(paste(tau[italic(M)])),
                    expression(paste(delta[italic(Us)])),
                    expression(paste(delta[italic(LL)])),
                    expression(paste(alpha[italic(m)])),
                    expression(paste(gamma[italic(s)])),
                    expression(italic(M[Us])),
                    expression(paste(sigma)), 
                    expression(paste(delta[italic(Lm)])),
                    expression(italic(r)), 
                    expression(italic(K[A])),
                    expression(italic(M[L])),
                    expression(italic(F[U])),
                    expression(paste(beta)),
                    expression(italic(F[L])),
                    expression(paste(delta[italic(Ls)]))),
        horiz=T, xlim=c(0,0.4), las=1, xlab="Importance Value")
dev.off()

#Create histogram of volatility for all 4K parameter combinations
#pdf("HistogramPredFirst.pdf")
par(mai=c(1.02,1.22,0.862,0.62))
hist(predGSA[,1], breaks=seq(0,4,by=0.15), las=1, col="grey", main="", xlab="Volatility", labels=T)
#dev.off()

###########################################################
###############Create CART tree figure####################
###########################################################
#Classification tree for mean trophic level
cartpredvol<-rpart(Volatility~r+Ka+Kul+sigma+deltaUs+deltaUm+deltaUl+
                 aM+aL+b+gammaS+gammaM+deltaLs+deltaLm+
                 deltaLl+MUs+MUm+MUl+FU+taus+taum+taul+
                 ML+FL, data=predGSA, control=rpart.control(xval=10,minsplit=20))
#, method="pois",  #use Poisson dist.
print(cartpredvol)
plot(cartpredvol)
text(cartpredvol, use.n=T, cex=0.75)
printcp(cartpredvol)
plotcp(cartpredvol)#look at cost-complexity parameter plot to set threshold cp for pruned tree
prunedpredvol<-prune(cartpredvol,cp=0.012) #use cp value from cp plot to prune tree

#CART plot
#tiff("Pred_GSA_Cart.tif", width=5, height=5, units="in", res=800, pointsize=4.8, compression="lzw")
par(mai=c(1.32,1.42,0.42,0.62), cex.axis=1.5, cex.lab=1.5)
prp(prunedpredvol, extra=1, yesno=T, cex=1.5)
#dev.off()



#####################  Haven't updated the 2 panel plot code below  #################
#two-panel GSA plot
#tiff("Synch_GSA.tif", width=8, height=10, units="in", res=800, pointsize=6)
par(mfrow=c(2,1),cex.axis=1.5, cex.lab=1.5)
par(mai=c(0.5, 1.5, 0.5, 1.5))
barplot(mtlkelp.imp, 
        names.arg=c(expression(italic(K[UL])),
                    expression(paste(delta[italic(LL)])),
                    expression(italic(M[UL])),
                    expression(paste(tau[italic(M)])),
                    expression(paste(gamma[italic(M)])),
                    expression(paste(delta[italic(UL)])),
                    expression(italic(F[U])), 
                    expression(paste(tau[italic(L)])),
                    expression(paste(tau[italic(S)])), 
                    expression(italic(M[Um])),
                    expression(paste(alpha[italic(L)])),
                    expression(italic(M[Us])),
                    expression(paste(delta[italic(Um)])), 
                    expression(paste(delta[italic(Us)])),
                    expression(paste(alpha[italic(m)])),
                    expression(paste(delta[italic(Lm)])), 
                    expression(paste(gamma[italic(s)])),
                    expression(italic(r)), 
                    expression(italic(K[A])), 
                    expression(paste(sigma)), 
                    expression(italic(M[L])), 
                    expression(paste(beta)),
                    expression(paste(delta[italic(Ls)])),
                    expression(italic(F[L]))),
        horiz=T, xlim=c(0,0.5), las=1, xlab="Importance Value")
par(mai=c(1, 2, 1, 2))
prp(prunedmtl, extra=1, yesno=F, cex=1.5)
dev.off()

