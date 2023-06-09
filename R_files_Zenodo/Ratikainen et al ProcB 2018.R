#############################################################################################################
### Code for:
### Differential allocation of parental investment and the trade-off between size and number of offspring
### Proceeding of the Royal Society B. 2018. doi: 10.1098/rspb.2018.1074
### Ratikainen et al. 
### Code for both stochastic dynamic optimisation figures and simulations
#############################################################################################################


#Preamble----
#rm(list=ls())
library(abind)

Horizon <- 100 	#Number of time steps
Xmax <- 100	# Maximum energetic state
Xcrit <- 2	# Minimum energetic state for survival
Xmin <- Xcrit+1
y <- 30 	# value of food
lambda <- 0.8 	# probability of finding food
alpha <- 0.8 	# probability of surviving a season
I <- seq(0,Xmax-Xmin) # Investment possibilities
MaxNo <- 15 #2 	# maximum number of offspring, o_max

Pm <- c(0.3,0.4,0.3) #Probabilities of meeting different males. Baseline: c(0.3,0.4,0.3), also try c(0.6,0.2,0.2) and c(0.1,0.5,0.4)
Nm <- length(Pm)

# Cost function parameters
Ca <- c(0,0,0) #Use (10,5,0) for 2a, 3A: (0,5,10), Baseline: (0,0,0) 
Cc <- c(0,0,0) #Use (5,0,-5) for 2c, Baseline: (0,0,0) otherwise
a <- c(1.3,1.3,1.3)
b <- c(0.2,0.2,0.2)  #Change for scenario 2B c(0.22,0.2,0.18). Baseline: c(0.2,0.2,0.2) #Scenario 3A: (0.17,0.2,0.23) 
Cpars <- rbind(a,b)

# Offspring fitness function parameters
Wa <- c(-6,-3,0) 	#Use (-6,-3,0) for Scenario 1A, (0,0,0) otherwise
K <- c(11,11,11) 	#Upper asymptote c(11,11,11)   scenario 1D: c(5,8, 12)
A <-  c(-11,-11,-11) 	#Lower asymptote c(-11,-11,-11)  scenario 1D: c(-5,-8,-12)
B <- c(0.2,0.2,0.2) 	#g Scenario 1B: c(0.05, 0.125, 0.2) Baseline: c(0.2,0.2,0.2) scenario 1D: c(0.15,0.10,0.05) 3A:  c(0.1, 0.15, 0.2)
M <- c(5,5,5)  		#q Inflection point c(5,5,5), Scenario 1C: c(15, 10, 5) Scenario 1D: c(6,6,6)
Wpars <- rbind(K,A,B,M) 

F.vectors <- array(0,dim=c(Xmax,Nm,3)) # Terminal fitness is zero.
W <- c(0,0,0) 		# Empty fitness vector 
Highest <- array(0,dim=c(Nm,2)) # Empty vector for the investment loop

RHS <- array(0,dim=c(length(I),MaxNo,Nm))

#Fitness function----
FITNESS <- function(X,No,Xcrit,I,Xmax,Xmin,y,lambda,alpha,F.vectors) {
  InvestO<-I/No 	#Resources pr offspring
  for(Mqual in 1:Nm) {  #Loop over all male qualities
    Cost <- Cpars[2,Mqual]*(I+Cc[Mqual])+Cpars[1,Mqual]^(Cpars[2,Mqual]*(I+Cc[Mqual]))+Ca[Mqual]  #+No*No     #Add last part if there is to be increasings costs with number of offspring
    Xi <- X-Cost 	#Resources left after investment
    Lower <- min(floor(Xi),Xmax)
    Upper <- min(Lower+1,Xmax)
    p2 <- Xi-Lower
    p1 <- 1-p2
    if(Lower>Xcrit) { 	#only if there is enough resources to survive
      Lfood <- min(Lower+y,Xmax)
      Ufood <- min(Upper+y,Xmax)
      #Future fitness expectations
      Fut <- p1*alpha*(lambda*(F.vectors[Lfood,2,1]*Pm[1] + F.vectors[Lfood,2,2]*Pm[2] + F.vectors[Lfood,2,3]*Pm[3])+(1-lambda)*(F.vectors[Lower,2,1]*Pm[1] + F.vectors[Lower,2,2]*Pm[2] + F.vectors[Lower,2,3]*Pm[3]))+
        p2*alpha*(lambda*(F.vectors[Ufood,2,1]*Pm[1] + F.vectors[Ufood,2,2]*Pm[2] + F.vectors[Ufood,2,3]*Pm[3])+(1-lambda)*(F.vectors[Upper,2,1]*Pm[1] + F.vectors[Upper,2,2]*Pm[2] + F.vectors[Upper,2,3]*Pm[3]))
      if(I==0) {   	#Need to make fitness 0 if there is 0 investment
        Woff<-0
      } else {
        Woff <- Wa[Mqual]+Wpars[2,Mqual]+(Wpars[1,Mqual]-Wpars[2,Mqual])/(1+exp(-Wpars[3,Mqual]*(InvestO-Wpars[4,Mqual]))) #Offspring fitness function 
      }
      W[Mqual] <- Woff*No+Fut 
    } else {
      W[Mqual] <- 0}   
  }
  return(W)
}

#Program----
OVER.INV <- function(X,Xcrit,I,Xmax,Xmin,y,lambda,alpha,F.vectors) {
  for(i in 1:length(I)) {
    for (No in 1:MaxNo){  
      RHS[i,No,] <- FITNESS(X,No,Xcrit,I[i],Xmax,Xmin,y,lambda,alpha,F.vectors)
    }
  }
  for(j in 1:Nm) {
    F.vectors[X,1,j]<- max(RHS[,,j]) #The highest W value
    Highest[j,] <- which(RHS[,,j]==max(RHS[,,j]), arr.ind = TRUE)[1,] #The Inv that gave the highest W value, just pick lowest investment if there are several that give max fitness
    
  }
  Temp <- abind(F.vectors[X,1,],Highest)
  return(Temp)
}

OVER.STATES <- function(F.vectors,Xcrit,Xmax,Xmin,y,lambda,alpha) {
  Store <- array(0,dim=c(Xmax,Nm,3))
  for(X in Xmin:Xmax) {
    Temp <- OVER.INV(X,Xcrit,I,Xmax,Xmin,y,lambda,alpha,F.vectors)
    #n <- nrow(Temp)-1
    F.vectors[X,1,] <- Temp[,1]
    Store[X,,] <- Temp 
  }
  Temp <- abind(F.vectors,Store)
  return(Temp)
}

FxtT <- array(0,dim=c(Horizon,Xmax,Nm))
Best.Inv <- array(0,dim=c(Horizon,Xmax,Nm,2))

Time <- Horizon
while(Time>1) {
  Time <- Time-1
  Temp <- OVER.STATES(F.vectors,Xcrit,Xmax,Xmin,y,lambda,alpha)
  TempF <- Temp[,,1:Nm]
  F.vectors[Xmin:Xmax,2,] <- TempF[Xmin:Xmax,1,]
  Best.Inv[Time,,,] <- Temp[,,5:6]
  FxtT[Time,,] <- Temp[,,4]
}

X <- seq(1,Xmax)
Best.Inv[Horizon,,,] <- X
FxtT[Horizon,,] <- X


##Let's produce the graphs!

#Plots optimal investments from stochastic dynamic model----
par(mfrow=c(1,3), mar=c(4,4,1,1),oma=c(0,0,0,0))
plot(Xmin:Xmax,I[Best.Inv[1,Xmin:Xmax,1,1]],type="l",col="gold",lwd=2,
     ylim=c(0,Xmax),xlab="State",ylab="Total investment",cex.lab=1.3,main="A",bty="n")
lines(Xmin:Xmax,I[Best.Inv[1,Xmin:Xmax,2,1]],col="orange",lwd=2)
lines(Xmin:Xmax,I[Best.Inv[1,Xmin:Xmax,3,1]],col="red",lwd=2)

test<-which(Best.Inv[1,Xmin:Xmax,,1]==1, arr.ind = TRUE)
for(w in 1:length(test[,1])){
  Best.Inv[1,test[w,1],test[w,2],2]<-0
}

plot(Xmin:Xmax,Best.Inv[1,Xmin:Xmax,1,2],type="l",col="gold",lwd=2,
     ylim=c(0,MaxNo),xlab="State",ylab="Number of offspring",cex.lab=1.3,main="B",bty="n")
lines(Xmin:Xmax,Best.Inv[1,Xmin:Xmax,2,2],col="orange",lwd=2)
lines(Xmin:Xmax,Best.Inv[1,Xmin:Xmax,3,2],col="red",lwd=2)

plot(Xmin:Xmax,I[Best.Inv[1,Xmin:Xmax,1,1]]/Best.Inv[1,Xmin:Xmax,1,2],type="l",col="gold",lwd=2,
     ylim=c(0,60),xlab="State",ylab="Investment pr offspring",cex.lab=1.3,main="C",bty="n")
lines(Xmin:Xmax,I[Best.Inv[1,Xmin:Xmax,2,1]]/Best.Inv[1,Xmin:Xmax,2,2],col="orange",lwd=2)
lines(Xmin:Xmax,I[Best.Inv[1,Xmin:Xmax,3,1]]/Best.Inv[1,Xmin:Xmax,3,2],col="red",lwd=2)



# Fitness and cost function plots ----
#First do calculations to get numbers to plots
Mqual<-1:3

WoffPlot<-array(0,dim=c(Nm,length(I)))
CostPlot<-array(0,dim=c(Nm,length(I)))
Invested<-c(0,0,0)
MaxFit<-c(0,0,0)
InvestedTot<-c(0,0,0)
OptCost<-c(0,0,0)
for (m in Mqual){
  WoffPlot[m,] <- Wa[m]+Wpars[2,m]+(Wpars[1,m]-Wpars[2,m])/(1+exp(-Wpars[3,m]*(I-Wpars[4,m]))) 
  CostPlot[m,] <- Cpars[2,m]*(I+Cc[m])+Cpars[1,m]^(Cpars[2,m]*(I+Cc[m]))+Ca[m]
  Invested[m]<-mean(I[Best.Inv[1,20:Xmax,m,1]])/mean(Best.Inv[1,20:Xmax,m,2])
  InvestedTot[m]<-mean(I[Best.Inv[1,20:Xmax,m,1]])
  MaxFit[m]<-Wa[m]+Wpars[2,m]+(Wpars[1,m]-Wpars[2,m])/(1+exp(-Wpars[3,m]*(Invested[m]-Wpars[4,m])))
  OptCost[m]<-Cpars[2,m]*(InvestedTot[m]+Cc[m])+Cpars[1,m]^(Cpars[2,m]*(InvestedTot[m]+Cc[m]))+Ca[m] 
}


#The actual plots
## First plot----
###NB will not produce correct graphs for scenarios with only 2 offspring. Also: colours adjustments needed to reproduce graphs in manuscript

par(mfrow=c(1,2))

plot(I,WoffPlot[3,], yaxs="i", xaxs="i",xlim=c(0,50),ylim=c(0,15),type="l",xlab="Investment pr offspring",ylab="Fitness",col="red",lwd=2,bty="n",main="A) Offspring fitness functions")
lines(I,WoffPlot[2,],col="orange",lwd=2)
lines(I,WoffPlot[1,],col="gold",lwd=2)
lines(c(0,1.5*Invested[1]),c(0,1.5*MaxFit[1]),lwd=2, col="grey",lty=2)
points(Invested[1],MaxFit[1],pch=8)
lines(c(0,1.5*Invested[2]),c(0,1.5*MaxFit[2]),lwd=2, col="grey",lty=2)
points(Invested[2],MaxFit[2],pch=8)
lines(c(0,1.5*Invested[3]),c(0,1.5*MaxFit[3]),lwd=2, col="grey",lty=2)
points(Invested[3],MaxFit[3],pch=8)#,col="orange",cex=2,lwd=2)

plot(I,CostPlot[1,], yaxs="i", xaxs="i",ylim=c(0,100),xlim=c(0,100),type="l",col="gold",xlab="total investment",ylab="Cost",lwd=2,bty="n",main="B) Parental cost function")
lines(I,CostPlot[2,],col="orange",lwd=2)
lines(I,CostPlot[3,],col="black",lwd=2)
points(InvestedTot[1],OptCost[1],pch=8,col="gold",cex=2, lwd=2)
points(InvestedTot[2],OptCost[2],pch=8,col="orange",cex=2, lwd=2)
points(InvestedTot[3],OptCost[3],pch=8,col="red",cex=2, lwd=2)


#################
#To plot offspring fitness function in the Scenario with increasing costs with number of offspring 
#the following numbers and Cost plot is needed ----
Mqual<-1:3
Nb<-1:MaxNo
WoffPlot<-array(0,dim=c(Nm,length(I)))
CostPlot<-array(0,dim=c(Nm,length(I),MaxNo))
for (m in Mqual){
  WoffPlot[m,] <- Wa[m]+Wpars[2,m]+(Wpars[1,m]-Wpars[2,m])/(1+exp(-Wpars[3,m]*(I-Wpars[4,m]))) #Offspring fitness function for use in Scenario BenPos, CostElev, CostSlope
  for (N in Nb){ 
    CostPlot[m,,N] <- Cpars[2,m]*(I+Cc[m])+Cpars[1,m]^(Cpars[2,m]*(I+Cc[m]))+Ca[m]+N*N 
    InvestedTot[m]<-mean(I[Best.Inv[1,20:Xmax,m,1]])
    OptCost[m]<-Cpars[2,m]*(InvestedTot[m]+Cc[m])+Cpars[1,m]^(Cpars[2,m]*(InvestedTot[m]+Cc[m]))+Ca[m]+(mean(Best.Inv[1,20:Xmax,m,2]))^2 
  }
}

plot(I,CostPlot[1,,1], yaxs="i", xaxs="i",ylim=c(0,100),xlim=c(0,100),type="l",col="grey",xlab="total investment",ylab="Cost",lwd=2,bty="n",main="B) Parental cost function")
lines(I,CostPlot[1,,3],col="black",lwd=2)
lines(I,CostPlot[1,,5],col="grey",lwd=2)
points(InvestedTot[1],OptCost[1],pch=8,col="gold",cex=2, lwd=2)
points(InvestedTot[2],OptCost[2],pch=8,col="orange",cex=2, lwd=2)
points(InvestedTot[3],OptCost[3],pch=8,col="red",cex=2, lwd=2)


##########################
#Lots of simulations
#########################

# Can only run after the stochastic dynamic optimisation model!

Nreplicates<-10000

Output <- matrix(0,(Horizon-1)*1000,6)

Time <- seq(1,Horizon-1)
reptime <- 1
for(Replicate in 1:Nreplicates) {
  X <- round(rnorm(1,mean=Xmax/1.5,sd=8)) #Starting state
  for(i in 1:(Horizon-1))
    if(X>Xcrit) {
      M <- sample(x=c(1:Nm), size=1, prob=Pm) #Choose male quality          
      Inv <- Best.Inv[i,X,M,1] #Find how much to invest in total  Best.Inv[1,Xmin:Xmax,2,1]
      No <- Best.Inv[i,X,M,2] #Find how many offspring you should invest in
      Xleft <- X-(Cpars[2,M]*(Inv+Cc[M])+Cpars[1,M]^(Cpars[2,M]*(Inv+Cc[M]))+Ca[M]) # +No*No) Add last part if costs are increasing with number of offspring
      Index <- 0
      if(runif(1)<lambda) { #check whether female finds food
        Index <- 1
      }
      Xleft <- Xleft+y*Index
      Xleft <- min(Xleft,Xmax)
      if(runif(1)>(Xleft-floor(Xleft))){
        Xleft<-floor(Xleft)
      } else{
        Xleft<-floor(Xleft)+1
      }
      if(Xleft<Xcrit){Xleft <- Xcrit}
      if(runif(1)>alpha){ #Check whether female gets eaten
        Xleft <- 0
      }
      Output[reptime,] <- c(Replicate,i,Inv,No,X,M)
      reptime <- reptime+1
      X<-Xleft
    }
}

#Simulation figure
#This is needed for the figure----
Output <- as.data.frame(Output)
names(Output) <-c ("Replicate","timestep","Inv","Number","X","M")
Survivors <- Output[Output$X>0,]
palette(c("Gold","Orange","Red"))

par(mfrow=c(2,2),mar=c(4,4,2,2))

##First panel----
#Means and sds for states
plot(Survivors$X,Survivors$Inv,xlab="Female state",ylab="Total investment",main="a)",cex.main=0.8,ylim=c(0,100),type="n")
#mtext(side=4,text="Average state",cex=0.8,col="gray44")
cols <- c("gold","orange","red")
avgs <- matrix(0,Nm,Xmax)
sds <- matrix(0,Nm,Xmax)
bars <- rep(0,Xmax)
for(i in 1:Xmax){
  points(i,length(Survivors$Inv[Survivors$X==i])/100,type="h",col="gray88",lwd=4)
  for(j in 1:Nm){
    avgs[j,i] <- mean(Survivors$Inv[Survivors$M==j&Survivors$X==i])
    sds[j,i] <- sd(Survivors$Inv[Survivors$M==j&Survivors$X==i])
    segments(i,avgs[j,i]-sds[j,i],i,avgs[j,i]+sds[j,i],lwd=1)
    segments(i-0.1,avgs[j,i]-sds[j,i],i+0.1,avgs[j,i]-sds[j,i],lwd=1)
    segments(i-0.1,avgs[j,i]+sds[j,i],i+0.1,avgs[j,i]+sds[j,i],lwd=1)
  }
}
for(i in 1:Nm){
  points(1:Xmax,avgs[i,])
  lines(1:Xmax,avgs[i,],col=cols[i],lwd=2)
}

## Second fig----
#Means and sds for time steps
plot(Survivors$timestep,Survivors$Inv,xlab="Time step",main="b)",xlim=c(0,50),cex.main=0.8,ylab="Total investment,",ylim=c(0,100),type="n")
#mtext(side=4,text="# Survivors",cex=0.8,col="gray44")
avgs <- matrix(0,Nm,Horizon)
sds <- matrix(0,Nm,Horizon)
line <- rep(0,Horizon)
for(i in 1:Horizon){
  points(i,length(Survivors$Inv[Survivors$timestep==i])/100,type="h",col="gray88",lwd=4)
  for(j in 1:Nm){
    avgs[j,i] <- mean(Survivors$Inv[Survivors$M==j&Survivors$timestep==i])
    sds[j,i] <- sd(Survivors$Inv[Survivors$M==j&Survivors$timestep==i])
    segments(i,avgs[j,i]-sds[j,i],i,avgs[j,i]+sds[j,i],lwd=1)
    segments(i-0.1,avgs[j,i]-sds[j,i],i+0.1,avgs[j,i]-sds[j,i],lwd=1)
    segments(i-0.1,avgs[j,i]+sds[j,i],i+0.1,avgs[j,i]+sds[j,i],lwd=1)
  }
}
for(i in 1:Nm){
  points(1:Horizon,avgs[i,])
  lines(1:Horizon,avgs[i,],col=cols[i],lwd=2)
}
lines(1:Horizon,line,col="gray44")


##Third fig----
#Means and sds for states
plot(Survivors$X,Survivors$Number,xlab="Female state",ylab="Offspring number",main="c)",cex.main=0.8,ylim=c(0,15),type="n")
#mtext(side=4,text="Average state",cex=0.8,col="gray44")
cols <- c("gold","orange","red")
avgs <- matrix(0,Nm,Xmax)
sds <- matrix(0,Nm,Xmax)
bars <- rep(0,Xmax)
for(i in 1:Xmax){
  points(i,length(Survivors$Number[Survivors$X==i])/(Nreplicates/20),type="h",col="gray88",lwd=4)
  for(j in 1:Nm){
    avgs[j,i] <- mean(Survivors$Number[Survivors$M==j&Survivors$X==i])
    sds[j,i] <- sd(Survivors$Number[Survivors$M==j&Survivors$X==i])
    segments(i,avgs[j,i]-sds[j,i],i,avgs[j,i]+sds[j,i],lwd=1)
    segments(i-0.1,avgs[j,i]-sds[j,i],i+0.1,avgs[j,i]-sds[j,i],lwd=1)
    segments(i-0.1,avgs[j,i]+sds[j,i],i+0.1,avgs[j,i]+sds[j,i],lwd=1)
  }
}
for(i in 1:Nm){
  points(1:Xmax,avgs[i,])
  lines(1:Xmax,avgs[i,],col=cols[i],lwd=2)
}

#Fourth panel----
#Means and sds for time steps
plot(Survivors$timestep,Survivors$Number,xlab="Time step",main="d)",xlim=c(0,50),cex.main=0.8,ylab="Offspring number",ylim=c(0,15),type="n")
avgs <- matrix(0,Nm,Horizon)
sds <- matrix(0,Nm,Horizon)
line <- rep(0,Horizon)
for(i in 1:Horizon){
  points(i,length(Survivors$Number[Survivors$timestep==i])/(Nreplicates/15),type="h",col="gray88",lwd=4)
  for(j in 1:Nm){
    avgs[j,i] <- mean(Survivors$Number[Survivors$M==j&Survivors$timestep==i])
    sds[j,i] <- sd(Survivors$Number[Survivors$M==j&Survivors$timestep==i])
    segments(i,avgs[j,i]-sds[j,i],i,avgs[j,i]+sds[j,i],lwd=1)
    segments(i-0.1,avgs[j,i]-sds[j,i],i+0.1,avgs[j,i]-sds[j,i],lwd=1)
    segments(i-0.1,avgs[j,i]+sds[j,i],i+0.1,avgs[j,i]+sds[j,i],lwd=1)
  }
}
for(i in 1:Nm){
  points(1:Horizon,avgs[i,])
  lines(1:Horizon,avgs[i,],col=cols[i],lwd=2)
}
lines(1:Horizon,line,col="gray44")
