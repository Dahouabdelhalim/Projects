##Robert Dunn
##Three stage urchin model
##Show parameter space of FL x FU with stable/bi-stable regions

library(deSolve)


#####################################################################
################Parameter space calculations#########################
##################in R with deSolve##################################
#####################################################################
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

###----------parameters----------------------###
r= 10               #growth rate of kelp
Ka=3000        #carrying capacity of kelp
Kul=40        #large urchin carrying capacity
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

state<-c(A=1e3, Us=70, Um=70, Ul=70, L=20) #Initial state of the system
out=as.data.frame(lsoda(y=state, times=times, func=threestage, parms=parms))   #numerically integrate the logistic function
tail(out,10)  #check output

#Time Series Plot, check plot
plot(out$time,out$A, type="l", xlab="t", ylab="N", col="forestgreen",
     ylim=c(0,50),main="3 urchin stages")
lines(out$time, out$Us, type="l", col="purple")
lines(out$time,out$Um, type="l", col="blue")
lines(out$time,out$Ul, type="l", col="darkred")
lines(out$time,out$L, type="l", col="indianred2")

################################################################
###Define parameter space for FL x FU at sigma=0.0, 0.5 and 0.95
###Loop over FL from 0.0 to 1.0 and FU from 0.0 to 1.0 
###Run loops from kelp forest and barrens initial conditions
###Calculate difference in equilibrium abundance between kelp and barrens inits.
###Those with large diffs. have bi-stability
################################################################


##################Sigma=0.0#######################
#Kelp forest 
sigma=0.0
state<-c(A=1e3, Us=70, Um=70, Ul=70, L=20)
FLseq<-seq(0.0,1.0, length=41)
FUseq<-seq(0.0,1.0, length=41)
N<-array(NA, dim=c(length(FLseq), length(FUseq), 5))
for(i in 1:length(FLseq)){
  for(j in 1:length(FUseq)){
  FL=FLseq[i]
  FU=FUseq[j]
  parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
          deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
          deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
          MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
  out=lsoda(y=state, times=times, func=threestage, parms=parms)
  N[i,j,1]=out[200,2]   #kelp pop matrix
  N[i,j,2]=out[200,3]   #Us pop matrix
  N[i,j,3]=out[200,4]   #Um pop matrix
  N[i,j,4]=out[200,5]   #Ul pop matrix
  N[i,j,5]=out[200,6]   #lobster pop matrix
}}
sig00kelp<-N
#sig00kelp[,,1]

##############
#Urchin barren
sigma=0.0
state<-c(A=.1, Us=20, Um=20, Ul=20, L=.1)
N<-array(NA, dim=c(length(FLseq), length(FUseq), 5))
for(i in 1:length(FLseq)){
  for(j in 1:length(FUseq)){
    FL=FLseq[i]
    FU=FUseq[j]
    parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
            deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
            deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
            MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
    out=lsoda(y=state, times=times, func=threestage, parms=parms)
    N[i,j,1]=out[200,2]   #kelp pop matrix
    N[i,j,2]=out[200,3]   #Us pop matrix
    N[i,j,3]=out[200,4]   #Um pop matrix
    N[i,j,4]=out[200,5]   #Ul pop matrix
    N[i,j,5]=out[200,6]   #lobster pop matrix
  }}
sig00bar<-N
#sig00bar[,,1]

#################Sigma=0.5##############
#Kelp Forest
sigma=0.5
state<-c(A=1e3, Us=70, Um=70, Ul=70, L=20)
N<-array(NA, dim=c(length(FLseq), length(FUseq), 5))
for(i in 1:length(FLseq)){
  for(j in 1:length(FUseq)){
    FL=FLseq[i]
    FU=FUseq[j]
    parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
            deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
            deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
            MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
    out=lsoda(y=state, times=times, func=threestage, parms=parms)
    N[i,j,1]=out[200,2]   #kelp pop matrix
    N[i,j,2]=out[200,3]   #Us pop matrix
    N[i,j,3]=out[200,4]   #Um pop matrix
    N[i,j,4]=out[200,5]   #Ul pop matrix
    N[i,j,5]=out[200,6]   #lobster pop matrix
  }}
sig05kelp<-N

###############
#Urchin barren
state<-c(A=.1, Us=20, Um=20, Ul=20, L=.1)
N<-array(NA, dim=c(length(FLseq), length(FUseq), 5))
for(i in 1:length(FLseq)){
  for(j in 1:length(FUseq)){
    FL=FLseq[i]
    FU=FUseq[j]
    parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
            deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
            deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
            MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
    out=lsoda(y=state, times=times, func=threestage, parms=parms)
    N[i,j,1]=out[200,2]   #kelp pop matrix
    N[i,j,2]=out[200,3]   #Us pop matrix
    N[i,j,3]=out[200,4]   #Um pop matrix
    N[i,j,4]=out[200,5]   #Ul pop matrix
    N[i,j,5]=out[200,6]   #lobster pop matrix
  }}
sig05bar<-N
#sig05bar[,,1]

###############Sigma=0.95##############
#Kelp Forest
sigma=0.95
state<-c(A=1e3, Us=70, Um=70, Ul=70, L=20)
N<-array(NA, dim=c(length(FLseq), length(FUseq), 5))
for(i in 1:length(FLseq)){
  for(j in 1:length(FUseq)){
    FL=FLseq[i]
    FU=FUseq[j]
    parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
            deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
            deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
            MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
    out=lsoda(y=state, times=times, func=threestage, parms=parms)
    N[i,j,1]=out[200,2]   #kelp pop matrix
    N[i,j,2]=out[200,3]   #Us pop matrix
    N[i,j,3]=out[200,4]   #Um pop matrix
    N[i,j,4]=out[200,5]   #Ul pop matrix
    N[i,j,5]=out[200,6]   #lobster pop matrix
  }}
sig095kelp<-N

###########
#Urchin barren
state<-c(A=.1, Us=20, Um=20, Ul=20, L=.1)
N<-array(NA, dim=c(length(FLseq), length(FUseq), 5))
for(i in 1:length(FLseq)){
  for(j in 1:length(FUseq)){
    FL=FLseq[i]
    FU=FUseq[j]
    parms=c(r=r, Ka=Ka, Kul=Kul, sigma=sigma,deltaUs=deltaUs, deltaUm=deltaUm, 
            deltaUl=deltaUl,aM=aM, aL=aL, b=b, gammaS=gammaS, gammaM=gammaM, 
            deltaLs=deltaLs, deltaLm=deltaLm,deltaLl=deltaLl, MUl=MUl, MUm=MUm,
            MUs=MUs, FU=FU, taus=taus, taum=taum, taul=taul, ML=ML, FL=FL)
    out=lsoda(y=state, times=times, func=threestage, parms=parms)
    N[i,j,1]=out[200,2]   #kelp pop matrix
    N[i,j,2]=out[200,3]   #Us pop matrix
    N[i,j,3]=out[200,4]   #Um pop matrix
    N[i,j,4]=out[200,5]   #Ul pop matrix
    N[i,j,5]=out[200,6]   #lobster pop matrix
  }}
sig095bar<-N

###################################################################
##Calculate differences for kelp forest state for each value of sigma
#####################################################################

kelpdiff00<-abs(sig00kelp[,,1]-sig00bar[,,1])
kelpdiff00<-round(kelpdiff00,3)
kelpdiff05<-abs(sig05kelp[,,1]-sig05bar[,,1])
kelpdiff05<-round(kelpdiff05,3)
kelpdiff095<-abs(sig095kelp[,,1]-sig095bar[,,1])
kelpdiff095<-round(kelpdiff095,3)

############################
####subsetting for plots###
##########################
###Sigma=0.0
plot00bi<-subset(kelpdiff00>50)  #determine where each is bistable
plot00bi<-matrix(as.integer(plot00bi), dim(plot00bi))
plot00bi[plot00bi==1]<-2  #look at plot00bi for where bistability occurs

plot00kelp<-subset(sig00kelp[,,1]>200) #pull out only those in the kelp forest state
plot00kelp<-matrix(as.integer(plot00kelp), dim(plot00kelp)) #turn T/F into numbers

plot00kelp[22,1]<-2       #from plot00bi, manually input where bistability occurs
plot00kelp[23,1:2]<-2
plot00kelp[24,1:3]<-2
plot00kelp[25,1:5]<-2
plot00kelp[26,1:8]<-2
plot00kelp[27,1:15]<-2
plot00kelp[28,1:34]<-2
plot00kelp[29,4:41]<-2
plot00kelp[30,32:41]<-2
#plot
colors<-c("white", "darkgrey", "black")  #set up color vector for black=kelp, gray=bistable, white=barren
sig0plot<-image(plot00kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")

############
#Sigma=0.5
plot05bi<-subset(kelpdiff05>50) #Check where is bistability
plot05bi<-matrix(as.integer(plot05bi), dim(plot05bi))
plot05bi[plot05bi==1]<-2

plot05kelp<-subset(sig05kelp[,,1]>200)
plot05kelp<-matrix(as.integer(plot05kelp), dim(plot05kelp))

plot05kelp[22,1]<-2   #check plot05bi, input region of bistability
plot05kelp[23,1:2]<-2
plot05kelp[24,1:3]<-2
plot05kelp[25,1:5]<-2
plot05kelp[26,1:9]<-2
plot05kelp[27,1:17]<-2
plot05kelp[28,5:41]<-2
plot05kelp[29,21:41]<-2
#plot
sig05plot<-image(plot05kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")

##########
#######Sigma=0.95
plot095bi<-subset(kelpdiff095>50)
plot095bi<-matrix(as.integer(plot095bi), dim(plot095bi))
plot095bi[plot095bi==1]<-2

plot095kelp<-subset(sig095kelp[,,1]>200)
plot095kelp<-matrix(as.integer(plot095kelp), dim(plot095kelp))

plot095kelp[22,1]<-2   #check plot095bi, input bistable region
plot095kelp[23,1:2]<-2
plot095kelp[24,1:4]<-2
plot095kelp[25,1:7]<-2
plot095kelp[26,7:16]<-2


#plot
sig095plot<-image(plot095kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")

#hi resolution plot of all 3
tiff("ParameterSpace.tif", width=6, height=4, units="in", res=600, pointsize=5)
par(mfrow=c(3,2), mar=c(5,5,4,4), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5)
image(plot00kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")
box(which="plot", lty="solid")
image(plot05kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")
box(which="plot", lty="solid")
image(plot095kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")
box(which="plot", lty="solid")
dev.off()


#PDF of all 3
pdf("ParameterSpace.pdf")
par(mfrow=c(3,2))
image(plot00kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")
box(which="plot", lty="solid")
image(plot05kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")
box(which="plot", lty="solid")
image(plot095kelp, col=colors, las=1, xlab="Lobster fishing mortality", ylab="Urchin fishing mortality", bty="O")
box(which="plot", lty="solid")
dev.off()
