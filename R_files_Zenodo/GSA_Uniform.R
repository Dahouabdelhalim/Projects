##3 stage urchin community dynamics model
##Robert Dunn

##GSA using Mean trophic level as the model outcome

##Create range of values of urchin parameters from which to sample (uniform dist. randomly sampled)

##Make 4k sample combinations, store them, then run 200y simulations for each

##Store the equilibrium pop. size of each state variable for each simulation

##See the random Forest vignette for specific questions about coding, functions

#Did GSA with Random Forest over two different sets of parameter value/simulation
#combinations with 4000 parameter combos, also did randomForest GSA with a different
#set of 2000 parameter combos. All analyses provided nearly exactly the same results.

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
                               runif(sims,0.0,0.975), # sigma
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

##Problems?
##try RK4 solver
##try reducing the step size/tolerance

###----Run a simulation with each set of parameters chosen above-------###
tf= 200          #run time
times<-1:tf      #times vector

N <-matrix(NA,5,ncol(parameters))  #create vector with NA's to store equilibrium state var. population size
state<-c(A=1000, Us=2, Um=2, Ul=.5, L=5)     #Kelp forest inits

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
  
  out=lsoda(y=state, times=times, func=threestage, parms=parms)#, method="lsode",maxsteps=1e4, verbose=TRUE)
  
  N[1,j]=out[tf,2]   #kelp pop matrix
  N[2,j]=out[tf,3]   #Us pop matrix
  N[3,j]=out[tf,4]   #Um pop matrix
  N[4,j]=out[tf,5]   #Ul pop matrix
  N[5,j]=out[tf,6]   #lobster pop matrix
}

##Set up data for kelp GSA
kelpGSA<-matrix(NA,sims,29)   #dummy matrix to hold state variables and sampled parameter values
kelpGSA[,1:5]<-t(N)          #input state variable values from simulations
kelpGSA[,6:29]<-t(parameters) #input randomly chosen parameter values
kelpGSA<-as.data.frame(kelpGSA)
colnames(kelpGSA)<-c("Kelp","SmallUrch","MedUrch","LargeUrch","Lobster",      #define parameter labels
                     "r","Ka","Kul","sigma","deltaUs","deltaUm", "deltaUl",
                     "aM","aL", "b", "gammaS", "gammaM", "deltaLs", "deltaLm",
                     "deltaLl","MUs","MUm","MUl", "FU","taus","taum","taul",
                     "ML", "FL")
##Calculate MTL for kelp
MTLKelp<- matrix(NA,1,nrow(kelpGSA))  #dummy df for mean trophic level output
MTLcalcdata<-as.data.frame(t(kelpGSA[,1:5]))
MTLcalcdata<-MTLcalcdata[-1,]  #drop kelp biomass from pop size df
MTLcalcdata<-rbind(MTLcalcdata,colSums(MTLcalcdata)) #add total biomass of inverts as new row
## now calculate mean trophic level
MTLKelp<- (2*MTLcalcdata[1,]/MTLcalcdata[5,]) +(2*MTLcalcdata[2,]/MTLcalcdata[5,])+
  (2*MTLcalcdata[3,]/MTLcalcdata[5,]) + (2.5*MTLcalcdata[4,]/MTLcalcdata[5,])
MTLKelp<-t(MTLKelp)  #need to transpose MTL to match the structure of the GSA parameter object

#Run Kelp GSA
MTLkelpGSA<-cbind(MTLKelp,kelpGSA) #combine MTL values with parameter values
colnames(MTLkelpGSA)[1] <- "MTL"

##Run random forest on mean trophic level data
mtlkelp.rf<-randomForest(MTLkelpGSA[,1] ~ .,
                     data       = MTLkelpGSA[,7:30],
                     importance = TRUE)
print(mtlkelp.rf)
round(importance(mtlkelp.rf),2)
mtlkelp.imp = pmax (importance (mtlkelp.rf,scale = FALSE) [, 1], 0)
mtlkelp.imp = mtlkelp.imp/sum(mtlkelp.imp)

#Bar plot of normalized importance values
pdf("ImportanceMTLKelp.pdf")
par(mai=c(1.02,1.02,0.62,0.42))
barplot (mtlkelp.imp,
         horiz = T,
         main  = 'Mean trophic level',
         xlim  = c (0, .5),  las=1)
dev.off()

#hi resolution plot
tiff("testGSAfig.tif", width=3, height=3, units="in", res=600, pointsize=4)
par(mar=c(5,5,4,4), xaxs="i", yaxs="i", cex.axis=1.3, cex.lab=1.3)
barplot (mtlkelp.imp,
         horiz = T,
         main  = 'Mean trophic level',
         xlim  = c (0, .5),  las=1, xlab="Importance value")
dev.off()

#Sorted bar plot of normalized importance values in base graphics 
importance<-as.data.frame(mtlkelp.imp)
mtlkelp.imp<-sort(mtlkelp.imp, decreasing=F)
barplot(mtlkelp.imp, horiz=T, xlim=c(0,0.5), las=1, xlab="Importance value") #sorted plot without greek symbols
#updated plot with sorted labels
tiff("ImportanceMTLKelp2.tif", width=6, height=6, units="in", res=600, pointsize=6)
par(mai=c(1.32,1.42,0.42,0.62), cex.axis=1.5, cex.lab=1.5)
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
dev.off()



#Create histogram of MTL for all 4K parameter combinations
pdf("HistogramMTLKelp.pdf")
hist(MTLkelpGSA[,1], breaks=20, las=1, col="grey", main="")#, labels=T)
dev.off()

###########################################################
###############Create CART tree figure####################
###########################################################
#Classification tree for mean trophic level
cartkelpmtl<-rpart(MTL~r+Ka+Kul+sigma+deltaUs+deltaUm+deltaUl+
                 aM+aL+b+gammaS+gammaM+deltaLs+deltaLm+
                 deltaLl+MUs+MUm+MUl+FU+taus+taum+taul+
                 ML+FL, data=MTLkelpGSA, control=rpart.control(xval=10,minsplit=20))
#, method="pois",  #use Poisson dist.
print(cartkelpmtl)
plot(cartkelpmtl)
text(cartkelpmtl, use.n=T, cex=0.75)
printcp(cartkelpmtl)
plotcp(cartkelpmtl)#look at cost-complexity parameter plot to set threshold cp for pruned tree
prunedmtl<-prune(cartkelpmtl,cp=0.012) #use cp value from cp plot to prune tree

#CART plot
tiff("CARTkelpmtl2.tif", width=4, height=4, units="in", res=600, pointsize=4)
par(mai=c(1.32,1.42,0.42,0.62), cex.axis=1.5, cex.lab=1.5)
prp(prunedmtl, extra=1, yesno=F, cex=1)
dev.off()



############################################
###Run simulations for barrens state inits##
############################################
N <-matrix(NA,5,ncol(parameters)) 
state<-c(A=1, Us=30, Um=10, Ul=5, L=1) #barrens state inits 
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
  
  out=lsoda(y=state, times=times, func=threestage, parms=parms)#, method="lsode",maxsteps=1e4, verbose=TRUE)
  
  N[1,j]=out[tf,2]   #kelp pop matrix
  N[2,j]=out[tf,3]   #Us pop matrix
  N[3,j]=out[tf,4]   #Um pop matrix
  N[4,j]=out[tf,5]   #Ul pop matrix
  N[5,j]=out[tf,6]   #lobster pop matrix
}

##Set up data for barrens GSA
barGSA<-matrix(NA,sims,29)   #dummy matrix to hold state variables and sampled parameter values
barGSA[,1:5]<-t(N)          #input state variable values from simulations
barGSA[,6:29]<-t(parameters) #input randomly chosen parameter values
barGSA<-as.data.frame(barGSA)
colnames(barGSA)<-c("Kelp","SmallUrch","MedUrch","LargeUrch","Lobster",      #define parameter labels
                     "r","Ka","Kul","sigma","deltaUs","deltaUm", "deltaUl",
                     "aM","aL", "b", "gammaS", "gammaM", "deltaLs", "deltaLm",
                     "deltaLl","MUs","MUm","MUl", "FU","taus","taum","taul",
                     "ML", "FL")
##Calculate MTL for barrens
MTLbar<- matrix(NA,1,nrow(barGSA))  #dummy df for mean trophic level output
MTLcalcdata<-as.data.frame(t(barGSA[,1:5]))
MTLcalcdata<-MTLcalcdata[-1,]  #drop kelp biomass from pop size df
MTLcalcdata<-rbind(MTLcalcdata,colSums(MTLcalcdata)) #add total biomass of inverts as new row
## now calculate mean trophic level
MTLbar<- (2*MTLcalcdata[1,]/MTLcalcdata[5,]) +(2*MTLcalcdata[2,]/MTLcalcdata[5,])+
  (2*MTLcalcdata[3,]/MTLcalcdata[5,]) + (2.5*MTLcalcdata[4,]/MTLcalcdata[5,])
MTLbar<-t(MTLbar)  #need to transpose MTL to match the structure of the GSA parameter object

MTLbarGSA<-cbind(MTLbar,barGSA) #combine MTL values with parameter values
colnames(MTLbarGSA)[1] <- "MTL"

#Run RF on barrens simulations
mtlbar.rf<-randomForest(MTLbarGSA[,1] ~ .,
                     data       = MTLbarGSA[,7:30],
                     importance = TRUE)
print(mtlbar.rf)
round(importance(mtlbar.rf),2)
mtlbar.imp = pmax (importance (mtlbar.rf,scale = FALSE) [, 1], 0)
mtlbar.imp = mtlbar.imp/sum(mtlbar.imp)


#Bar plot of normalized importance values
par(mai=c(1.02,1.02,0.62,0.42))
barplot (mtlbar.imp,
         horiz = T,
         main  = 'Mean trophic level',
         xlim  = c (0, .5),  las=1)
dev.off()

### MAY NEED TO RE_ORDER PARAMETER NAMES
#Sorted bar plot of normalized importance values in base graphics
importance<-as.data.frame(mtlbar.imp)
mtlbar.imp<-sort(mtlbar.imp, decreasing=F)
barplot(mtlbar.imp,horiz=T,xlim=c(0,0.5), las=1, xlab="Importance value")
pdf("ImportanceMTLBarren.pdf")
par(mai=c(1.32,1.42,0.42,0.62))
barplot(mtlbar.imp, 
        names.arg=c(expression(italic(K[UL])),
                    expression(paste(delta[italic(UL)])),
                    expression(paste(delta[italic(LL)])),
                    expression(italic(M[Um])),
                    expression(paste(tau[italic(M)])),
                    expression(paste(tau[italic(L)])),
                    expression(paste(tau[italic(S)])), 
                    expression(italic(M[Us])),
                    expression(paste(alpha[italic(L)])),
                    expression(italic(F[U])), 
                    expression(paste(gamma[italic(M)])),
                    expression(paste(delta[italic(Um)])), 
                    expression(italic(M[UL])),
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
dev.off()

#Create histogram of MTL for all 4K parameter combinations
pdf("HistogramMTLBar.pdf")
hist(MTLbarGSA[,1], breaks=20, las=1, xlab="Mean Trophic Level", col="grey", main="")#, labels=T)
dev.off()

###########################################################
###############Create CART tree figure####################
###########################################################
#Classification tree for mean trophic level
cartbarmtl<-rpart(MTL~r+Ka+Kul+sigma+deltaUs+deltaUm+deltaUl+
                     aM+aL+b+gammaS+gammaM+deltaLs+deltaLm+
                     deltaLl+MUs+MUm+MUl+FU+taus+taum+taul+
                     ML+FL, data=MTLbarGSA, control=rpart.control(xval=10,minsplit=20))
#, method="pois",  #use Poisson dist.
print(cartbarmtl)
plot(cartbarmtl)
text(cartbarmtl, use.n=T, cex=0.75)
printcp(cartbarmtl)
plotcp(cartbarmtl)#look at cost-complexity parameter plot to set threshold cp for pruned tree
prunedmtl<-prune(cartbarmtl,cp=0.012) #use cp value from cp plot to prune tree
pdf("CARTbarmtl.pdf")
prp(prunedmtl, main="Mean Trophic Level", extra=1, yesno=F, cex=0.74)
dev.off()

##Attempt at doing a GSA with an A.S.S. output metric, ie, for simulations when ASS are possible, what are 
##the parameter sensitivities?
##Problem: Too few runs (55 out of 4K) demonstrated alt. stable states, not enough reps to run random forest

States<-cbind(MTLkelpGSA[,1],MTLbarGSA[,1])
colnames(States)<-c("KelpMTL", "BarrensMTL")
Diffs<-States[,1]-States[,2]
States<-cbind(States,Diffs)
Sums<-States[,1]+States[,2]
States<-cbind(States,Sums)


Master<-cbind(States,MTLkelpGSA[,-1])
subMTLbar<-MTLbarGSA[,-c(1,7:30)]
Master<-cbind(Master,subMTLbar)     #Biomass from Kelp inits are in columns 5:9, biomass from barrens inits in columns 34:38
Bistable<-subset(Master,Diffs>0.001)
