##3 stage urchin community dynamics model
##Robert Dunn

##Predator-first recovery
##Management Case-study
##Testing community outcomes at 5, 10 and 15 years post-urchin fishery closure.

##Create range of values of parameters from which to sample (uniform dist. randomly sampled)
##Make 4k sample combinations, store them, then run simulations for each
##Simulations: 1: Run to fished equililbrium.
## 2. Close lobster fishery, run for 2 years.
## 3. Close urchin fishery, run for 15 years.
##Extract the fished community biomass at 5, 10 and 15 years post urchin fishery closure
##Plot histograms of fished community biomass to demonstrate the range of possible outcomes at each time horizon, given
##parameter uncertainty.
#setwd()
par(mar=c(5,4.5,2,2))  #set plot margins c(bottom, left, top, right)
rm(list=ls())
require(deSolve)      #ODE solver
require(nationalparkcolors) #nice Nat Park palettes
require(ggplot2)
voy<-park_palette("Voyageurs") #for kelp forest plots
arch<-park_palette("Arches")  #for urchin barren plots

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
#write.table(parameters,"parameters.txt", sep="\\t", row.names=F, col.names=TRUE)

###----Run a simulation with each set of parameters chosen above-------###
tf= 200          #run time for initial harvested equilibrium
times<-1:tf      #times vector
Ltimes<-1:2    #time scale to run lobster closure
Utimes<-1:15   #time scale to run full closure

N <-matrix(NA,17,ncol(parameters))  #create vector with NA's to store 17 years of population trajectories (2+15)
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
#  out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))  
  
  #Now, set up an MPA with 0 lobster fishing mortality, use harvested equilibrium population values from above
#  newstate<-c(A=out$A[200],Us=out$Us[200],Um=out$Um[200],Ul=out$Ul[200], L=out$L[200])
  FL=0.0    #close lobster fishery
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
  out2=as.data.frame(lsoda(y=state, times=Ltimes, func=threestage, parms=parms))  #run for 2 years  #changed from y=newstate
  
  #now stop fishing for urchins
  newerstate<-c(A=out2$A[2],Us=out2$Us[2],Um=out2$Um[2],Ul=out2$Ul[2], L=out2$L[2])
  FU=0.00    #close urchin fishery     
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL) 
  out3=as.data.frame(lsoda(y=newerstate, times=Utimes, func=threestage, parms=parms)) #now run for 15 years
  
  y<-rbind(out2,out3)  #not sure this is necessary.
  
  #Calculate community volatility##
  #sum biomass of large urchins and large lobsters for the full time series
  FishedDens<-y$Ul+y$L
 # maxfished<-max(FishedDens)  #calculate max value of aggregate fished community biomass
 # fishedvol<-(maxfished-FishedDens[400])/FishedDens[400]  #calculate fished community volatility
  
  N[,j]=FishedDens[1:17]   #recovery trajectory for a given simulation  #to keep the full 17 year time course for each sim, change to 1:17 and make N have 17 rows.
}
write.csv(N, file="PredFirst_RecoveryCaseStudy_BarrenInits.csv",row.names=F)
#can now just call the text files from working directory folder
K<-N   #now have a new matrix for the kelp initial condition results, whichever you ran initially
#U<-N   #now have a matrix for urchin barren init results, whichever you ran initially


#########################
#############PLOTTING#########
################################
par(mar=c(5,4.5,2,2))
mean_trajU<-rowMeans(U)  #calculate mean value for each time point
mean_trajK<-rowMeans(K)
matplot(y=U[,], type="l", lwd=1, lty=1, col=arch, bty="L", las=1, cex.lab=1.5, cex.axis=1.25,
        xlab="Time", ylab="Recovering biomass (kg)", ylim=c(0,120))
lines(mean_trajU, lwd=4, col="black")  #plot mean time series in dark blue
abline(v=2, col="grey", lty=2, lwd=2)

matplot(y=K[,], type="l", lwd=1, lty=1, col=voy, bty="L", las=1, cex.lab=1.5, cex.axis=1.25,
        xlab="Time", ylab="Recovering biomass (kg)", ylim=c(0,120))
lines(mean_trajK, lwd=4, col="black")  #plot mean time series in dark blue
abline(v=2, col="grey", lty=2, lwd=2)

#only plot subset of randomly selected recovery trajectories
tN<-t(U)  #transpose N matrix
randtrajs<-tN[sample(nrow(tN),50),] #select 50 trajectories to show
randtrajs<-t(randtrajs)
matplot(randtrajs[,],type="l", lty=1, lwd=2, col=arch, bty="L", las =1, ylim=c(0,50),
        cex.lab=1.5, cex.axis=1.25, xlab="Time (y)", ylab="Recovering biomass (kg)")
lines(mean_trajU, lwd=4, col="black")  #plot mean time series in dark blue
abline(v=2, col="grey", lty=2, lwd=2)

tK<-t(K)
randtrajK<-tK[sample(nrow(tK),50),]
randtrajK<-t(randtrajK)
matplot(randtrajK[,], type="l", lty=1, lwd=2, col=voy, bty="L", las=1, ylim=c(0,50),
        cex.lab=1.5, cex.axis=1.25, xlab="Time (y)", ylab="Recovering biomass (kg)")
lines(mean_trajK, lwd=4, col="black")  #plot mean time series in dark blue
abline(v=2, col="grey", lty=2, lwd=2)

### plot for change relative to initial conditions (unfished biomass)####
deltaFish<-K-7
matplot(y=deltaFish[,], type="l", lwd=1, lty=1, col=voy, bty="L", las=1, cex.lab=1.5, cex.axis=1.25,
        xlab="Time", ylab="Biomass change (kg)",ylim=c(-20,120))
abline(v=2, col="grey", lty=2, lwd=2)
deltaBarren<-U-5.5
matplot(y=deltaBarren[,], type="l", lwd=1, lty=1, col=arch, bty="L", las=1, cex.lab=1.5, cex.axis=1.25,
        xlab="Time", ylab="Biomass change (kg)", ylim=c(-20,120))
abline(v=2, col="grey", lty=2, lwd=2)
############################################
############################################
#Histograms for 5, 10, 15 year distributions
hist(K[7,],breaks=seq(0,200,1), xlim=c(0,80),ylim=c(0,600), main="",xlab="Recovering biomass (kg), 5 y after closure",
     cex.lab=1.5, cex.axis=1.2, las=1, col=scales::alpha('darkseagreen4', .75), border=F)
hist(U[7,],breaks=seq(0,200,1), add=T,col=scales::alpha('darkorchid3',.3), border=F )
legend(65, 500 ,legend=c("Kelp", "Barren"), box.lty=0, fill=c("darkseagreen3", "darkorchid1"), cex=1.25)

hist(K[12,],breaks=seq(0,200,1),xlim=c(0,80),ylim=c(0,600),main="",xlab="Recovering biomass (kg), 10 y after closure",
     cex.lab=1.5, cex.axis=1.2,las=1, col=scales::alpha('darkseagreen4', .75), border=F)
hist(U[12,], breaks=seq(0,200,1), add=T, col=scales::alpha('darkorchid3',.3), border=F)

hist(K[17,],breaks=seq(0,200,1), xlim=c(0,80),ylim=c(0,600),main="",xlab="Recovering biomass (kg), 15 y after closure",
     cex.lab=1.5, cex.axis=1.2,las=1, col=scales::alpha('darkseagreen4', .75), border=F)
hist(U[17,], breaks=seq(0,200,1), add=T, col=scales::alpha('darkorchid3', .3), border=F)

############################################
##density plots instead of histograms#
#################################
kelp5<-data.frame(K[7,])  ###First have to convert everything to long form to use in ggplot2#
urch5<-data.frame(U[7,])
kelp5$State<-"Kelp"
urch5$State<-"Barren"
names(kelp5)[1]<-"FishedDens"
names(urch5)[1]<-"FishedDens"
hist5<-rbind(kelp5,urch5)


d5<-ggplot(hist5, aes(FishedDens, col=State, fill=State)) + geom_density(alpha= 0.3, lwd=.75) #need to expand x axis, set colors
d5 +  ylab("Density") + xlab("Recovering biomass (kg) after 5 years") + theme_bw() + 
  guides(fill=guide_legend(title=NULL)) + guides(color=guide_legend(title=NULL))+
  theme(plot.margin=unit(c(.75,.75,.5,.5), "cm")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,80)) + scale_y_continuous(expand=c(0,0), limits=c(0,0.15))+
  scale_fill_manual(values=c("darkorchid3","darkseagreen4")) +
  scale_color_manual(values=c("darkorchid3","darkseagreen4")) +
  theme(legend.position=c(.85,.85))+ theme(legend.text=element_text(size=15))+
  theme(panel.border = element_blank(), axis.line.y=element_line(color="black", size=.5)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(axis.line=element_line(color="black", size=.05)) +
  theme(axis.text.x=element_text(color="black", size=15)) +
  theme(axis.text.y=element_text(color="black", size=15)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18))
  

###10 year data
kelp10<-data.frame(K[12,])
urch10<-data.frame(U[12,])
kelp10$State<-"Kelp"
urch10$State<-"Barren"
names(kelp10)[1]<-"FishedDens"
names(urch10)[1]<-"FishedDens"
hist10<-rbind(kelp10,urch10)

d10<-ggplot(hist10, aes(FishedDens, col=State, fill=State)) + geom_density(alpha= 0.3, lwd=.75) #need to expand x axis, set colors
d10 +  ylab("Density") + xlab("Recovering biomass (kg) after 10 years") + theme_bw() + 
  guides(fill=FALSE) + guides(color=FALSE) +
  theme(plot.margin=unit(c(.75,.75,.5,.5), "cm")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,80)) + scale_y_continuous(expand=c(0,0), limits=c(0,0.15))+
  scale_fill_manual(values=c("darkorchid3","darkseagreen4")) +
  scale_color_manual(values=c("darkorchid3","darkseagreen4")) +
  theme(panel.border = element_blank(), axis.line.y=element_line(color="black", size=.5)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(axis.line=element_line(color="black", size=.05)) +
  theme(axis.text.x=element_text(color="black", size=15)) +
  theme(axis.text.y=element_text(color="black", size=15)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18))

###15 year data
kelp15<-data.frame(K[17,])
urch15<-data.frame(U[17,])
kelp15$State<-"Kelp"
urch15$State<-"Barren"
names(kelp15)[1]<-"FishedDens"
names(urch15)[1]<-"FishedDens"
hist15<-rbind(kelp15,urch15)

d15<-ggplot(hist15, aes(FishedDens, col=State, fill=State)) + geom_density(alpha= 0.3, lwd=.75) #need to expand x axis, set colors
d15 +  ylab("Density") + xlab("Recovering biomass (kg) after 15 years") + theme_bw() + 
  guides(fill=FALSE) + guides(color=FALSE) +
  theme(plot.margin=unit(c(.75,.75,.5,.5), "cm")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,80)) + scale_y_continuous(expand=c(0,0), limits=c(0,0.15))+
  scale_fill_manual(values=c("darkorchid3","darkseagreen4")) +
  scale_color_manual(values=c("darkorchid3","darkseagreen4")) +
  theme(panel.border = element_blank(), axis.line.y=element_line(color="black", size=.5)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) + 
  theme(axis.line=element_line(color="black", size=.05)) +
  theme(axis.text.x=element_text(color="black", size=15)) +
  theme(axis.text.y=element_text(color="black", size=15)) +
  theme(axis.title.x=element_text(size=18)) +
  theme(axis.title.y=element_text(size=18))
