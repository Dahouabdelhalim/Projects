##---------------------------------------------------------------------##
# Code supplement for “Evolutionary history mediates population response to 
# rapid environmental change through within-generational and transgenerational 
# plasticity”, accepted by the American Naturalist.
# Run using R version 4.1.0 (2021-05-18)
# Last Edited: 11-21-22
##---------------------------------------------------------------------##

##---------------------------------------------------------------------##
# Load packages
library(lhs)
# Version 1.1.1
library(MASS)
# Version 7.3.54
library("lattice")
# Version 0.20.44
library("RColorBrewer")
# Version 1.1.2
##---------------------------------------------------------------------##

##---------------------------------------------------------------------##
# Note on environmental parameters
##---------------------------------------------------------------------##
#
# Parms: A list of environmental parameters with the following names,
#        l- environmental autocorrelation
#        K- maximum population growth rate. Note that this notation is outmoded. 
#           We use r in the manuscript, where K = exp(r). In all cases r = 1.
#        ve- environmental variance
#        vj- individual-level juvenile cue error
#        vme- individual-level maternal environment cue error
#        vmp- individual-level maternal phenotype cue error
#        mutheta- environmental mean
#        yj- population-level juvenile cue error
#        yme- population-level maternal environment cue error
#        muj- juvenile cue bias
#        mume- maternal environment cue bias
#        mump- maternal phenotype cue bias (not used in the manuscript; always set to zero)
#
##---------------------------------------------------------------------##

##---------------------------------------------------------------------##
# Note on cue weights
##---------------------------------------------------------------------##
# 
# CueWeights: A vector of cue weights described in order of entry
#         wj- The juvenile cue weight, given as w_j in the manuscript
#         wme- The maternal environment cue weight, given as w_me in the manuscript
#         wmp- The maternal phenotype cue weight, given as w_mp in the manuscript
#         wr- The randomization weight, given as w_r in the manuscript
#
##---------------------------------------------------------------------##


##---------------------------------------------------------------------##
# Functions
##---------------------------------------------------------------------##

Variance<-function(CueWeights=c(0,0,0,0,0),Parms){
  
  ###______________________________________________________________________________________________________
  ### Function purpose: Solve for equilibrium phenotypic variance
  ###______________________________________________________________________________________________________
  
  ### Inputs
  ###______________________________________________________________________________________________________
  ### v:           The value at which the polynomial is being evaluated.
  ###
  ### Parms:       A named list of environmental parameters
  ###
  ### Cue Weights: A vector of cue weights.
  ###______________________________________________________________________________________________________
  
  
  with(as.list(Parms),{
    # Round cue weights to prevent floating point errors
    wj<-round(CueWeights[2],10)
    wme<-round(CueWeights[3],10)
    wmp<-round(CueWeights[4],10)
    wr<-round(CueWeights[5],10)
    # Calculate the polynomial of v
    n2<-wj^2*vj+wme^2*vme+wmp^2*vmp+wr^2
    out<-(n2+wmp^2-1+sqrt((n2+wmp^2+1)^2-4*wmp^2))/2
    return(out)
  })
  
  ### Outputs
  ###______________________________________________________________________________________________________
  ### out:         Phenotypic Variance
  ###______________________________________________________________________________________________________
  
}

Fitness<-function(CueWeights=c(mu=0,wj=0,wme=0,wmp=0,wr=0),Parms,Info=F){
  ###______________________________________________________________________________________________________
  ### Function purpose: Determine fitness and components of fitness at stationarity
  ###______________________________________________________________________________________________________
  
  ### Inputs
  ###______________________________________________________________________________________________________
  ### Parms:       A named list of environmental parameters
  ###
  ### Cue Weights: A vector of cue weights.
  ###
  ### Info:        A logical scalar that determines whether extra fitness information is returned.
  ###______________________________________________________________________________________________________
  
  with(c(as.list(Parms),as.list(CueWeights)),{
    
    #Calculate stationary phenotypic variance
    v<-sqrt(Variance(CueWeights=CueWeights,Parms=Parms)+1)
    
    # Calculate the geometric mean fitness according to the formula
    # given in Supplement 1.
    
    A<-l*(wj-1)+wme+wmp
    B<-wmp/v^2

    # Equation S11- Covariance between the environment and phenotype-environment 
    # mismatch, Cov[Theta,D0]
    covTD<-ve/(1-l^2)*(l*(wme+wmp)+wj-1)/(1-l*B)
    
    # Variance of mean phenotype-environment mismatch (without population-level cue error)
    # Taken from the unlabeled equation at the bottom of page 6 of the supplement
    varD<-((A^2*ve)/(1-l^2)+(1-wj)^2*ve+A/(1-l)
           *(2*B*(1-l))*covTD)/(1-B^2)
    
    # Equation S13- Variance of mean phenotype-environment mismatch
    # (due to population-level cue error). Given in the text as
    # Var[Y'] = Var[Y]/(1-B^2).
    varY<-(wj^2*yj+wme^2*yme)/(1-B^2)
    
    # Equation S9- Expected mean phenotype-environment mismatch, E[D]
    eD<-(mu+(wme+wj+wmp)*(mutheta)+muj*wj+mume*wme-mutheta)/(1-wmp*v^-2)
      
    # Equation S16- Per-capita growth rate
    Fit<-log(K)-(1/2)*log(v^2)-(varD+eD^2+varY)/(2*v^2)
      
    # Return fitness and all fitness components
    if(Info){return(list("Fit"=Fit,"VX"=v^2-1,"PVL"= (1/2)*log(v^2),"EFL"=varD/(2*v^2),
                           "PEL"=varY/(2*v^2),"ESL"=eD^2/(2*v^2)))
      }else{return(Fit)}
  })
  
  ### Outputs
  ###______________________________________________________________________________________________________
  ### Fit:  Stationary mean log geometric fitness of the population.
  ###
  ### VX:   Stationary phenotypic variance of the population.
  ###
  ### PVL:  The phenotypic variance load
  ###
  ### EFL:  The environmental fluctuation load
  ###
  ### PEL:  The population-level error load
  ###
  ### ESL:  The environmental shift load
  ###______________________________________________________________________________________________________
  
}

GenLHS <- function(Parameters,NumSamples=50){
  ###______________________________________________________________________________________________________
  ### Function purpose: Generate Latin Hypercube Sample
  ###______________________________________________________________________________________________________
  
  ### Inputs
  ###______________________________________________________________________________________________________
  ### Parameters: An n x 2 matrix, where n is the number of parameters being sampled,
  ###             giving the upper and lower bounds over which each parameter is sampled.
  ###             This matrix is used to re-scale the uniform (0,1) sample generated by the randomLHS function
  ### NumSamples: Number of points in parameter space to sample.
  ###______________________________________________________________________________________________________
  
  smpl <- randomLHS(n= NumSamples,k=dim(Parameters)[1])
  
  for(i in 1:dim(Parameters)[1]){
    smpl[,i]<-Parameters[i,1]+(Parameters[i,2]-Parameters[i,1])*smpl[,i]
  }
  
  ### Outputs
  ###______________________________________________________________________________________________________
  ### smpl: The points in parameter space generated by the Latin hypercube sample.
  ###______________________________________________________________________________________________________
  
  return(smpl)
}


##---------------------------------------------------------------------##
##
## This section of code construct a series of contour plots on 
## the axes of environmental variation (VE) and population-level juvenile cue error (YJ).
##
## The contour plots include Figures 3a and 3c, all panels of Figure 5
## and Figure S2a and S2b.
##
##---------------------------------------------------------------------##

# Set grid of environmental parameters
ParMat4<-matrix(NA,2,17^2)
ParMat4[1,]<-rep(seq(0,4,0.25),each=17)
ParMat4[2,]<-rep(seq(0,4,0.25),17)

# Initialize matrices for data storage
Fit<-matrix(NA,17^2,7)
Fit2<-matrix(NA,17^2,4)
CueOpt3<-matrix(NA,17^2,5)
colnames(CueOpt3)<-c("mu","wj","wme","wmp","wr")
TempPar<-matrix(NA,10,5)
TempVal<-matrix(NA,10,1)

set.seed(1)
# For each set of parameters in the parameter grid...
for(i in 1:length(ParMat4[1,])){
  # ... set the parameters.
  Parms<-list(l=0.7,K=exp(1),ve=ParMat4[2,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  # Generate 10 initial conditions for the optimizer, using Latin Hypercube Sampling
  Inits<-GenLHS(matrix(c(-10,-10,-10,-10,0,10,10,10,10,10),5,2),10)
  colnames(Inits)<-c("mu","wj","wme","wmp","wr")
  # For each set of initial conditions...
  for(j in 1:dim(Inits)[1]){
    # ...find the optimal cue weights for the given environmental parameters.
    Solve <- optim(par=Inits[j,],fn=Fitness,Parms=Parms,method="L-BFGS-B",
                   control=list(fnscale=-1),lower=c(-10,-10,-10,-10,0),
                   upper=c(10,10,10,10,10))
    TempPar[j,]<-Solve$par
    TempVal[j]<-Solve$value
    # Round the value of mu to eliminate to numerical error
    TempPar[j,1]<-round(TempPar[j,1],3)
    
  }
  # If different intial conditions produce different optimal cue weights...
  if(sum((t(TempPar)-colMeans(TempPar))^2)>0.0001){
    # ...print the optimal cue weights for each initial condition.
    print(TempPar)
  }
  # Use the set of cue weights that correspond to the highest fitness value.
  CueOpt3[i,] <- TempPar[which.max(TempVal),]
  
  # Record fitness in pre-change environment under optimal cue weights.
  Fit[i,1]<-max(TempVal)
  
  # Calculate fitness after increase in mean environment by 1
  Parms<-list(l=0.7,K=exp(1),ve=ParMat4[2,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=1,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  Temp1 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,5]<-Temp1$ESL
  
  # Calculate fitness after increase in environmental variance by 0.4
  Parms<-list(l=0.7,K=exp(1),ve=ParMat4[2,i]+0.4,vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  Temp2 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,2] <- (Fit[i,1]-Temp2$Fit)
  Fit2[i,2]<-Temp2$Fit
  
  # Calculate fitness after increase in autocorrlation by 0.1
  Parms<-list(l=0.7-0.1,K=exp(1),ve=ParMat4[2,i]*(1-(0.7-0.1)^2)/(1-0.7^2),vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  Temp3 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,3] <- (Fit[i,1]-Temp3$Fit)
  Fit2[i,3]<-Temp3$Fit
  
  # Calculate fitness after increase in population-level juvenile cue error by 0.4
  Parms<-list(l=0.7,K=exp(1),ve=ParMat4[2,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i]+0.4,yme=0.25,muj=0,mume=0,mump=0)
  Temp4 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,4] <- (Fit[i,1]-Temp4$Fit)
  Fit2[i,4]<-Temp4$Fit
  
  # Calculate fitness after increase in individual-level juvenile cue error by 0.4
  Parms<-list(l=0.7,K=exp(1),ve=ParMat4[2,i],vj=0.25+0.4,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  Temp2 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,6] <- (Fit[i,1]-Temp2$Fit)
  
  # Calculate fitness after increase in juvenile cue bias by 1
  Parms<-list(l=0.7,K=exp(1),ve=ParMat4[2,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=1,mume=0,mump=0)
  Temp2 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,7] <- (Fit[i,1]-Temp2$Fit)
  
  print(i)
}

par(mar=c(5.1, 5, 4.1, 2.1))


#FitChange<-matrix(Fit[,5],17,17)
#levelplot(x=FitChange[-1,],xlab="Environmental Variance",cex.lab=1.75,cex.main=1.5,
#                    ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to change in ",mu[theta])),color=topo.colors)
#grid.arrange(Figure5a,ncol=1)

# Generate Figure 5a
pdf(file="Figure5a.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,5],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to change in ",mu[theta])),color=topo.colors,family="sans")
legend(c(0.025*4,0.1*4),c(0.975*4,0.85*4),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75,cex=2,bty="o",bg="white")
dev.off()

# Generate Figure 5b
pdf(file="Figure5b.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,7],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to change in ",mu[J])),color=topo.colors,family="sans")
legend(c(0.025*4,0.1*4),c(0.975*4,0.85*4),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Generate Figure 5c
pdf(file="Figure5c.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,2],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to increase in ",sigma[theta]^{2})),color=topo.colors,family="sans")
legend(c(0.025*4,0.1*4),c(0.975*4,0.85*4),
       c("C"),col=c(1),y.intersp=0,x.intersp=-0.75,cex=2,bty="o",bg="white")
dev.off()

# Generate Figure 5d
pdf(file="Figure5d.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,3],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to decrease in ",lambda)),color=topo.colors,family="sans")
legend(c(0.025*4,0.1*4),c(0.975*4,0.85*4),
       c("D"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()



# Generate Figure 3c
pdf(file="Figure3c.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,1],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=1.75,cex.main=1.75,
               ylab="Population Juvenile Cue Error",main=expression(paste(main="Per-Capita Growth Rate")),color=topo.colors,family="sans")
legend(c(0.025*4,0.1*4),c(0.975*4,0.85*4),
       c("C"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Generate Figure S2b
pdf(file="FigureS2b.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,4],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to increase in ",(sigma[J]^{pop})^{2})),color=topo.colors,family="sans")
legend(c(0.025*4,0.1*4),c(0.975*4,0.85*4),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Generate Figure S2a
pdf(file="FigureS2a.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,6],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to increase in ",(sigma[J]^{ind})^{2})),color=topo.colors,family="sans")
legend(c(0.025*4,0.1*4),c(0.975*4,0.85*4),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Generate Figure 3a
pdf(file="Figure3a.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(round(CueOpt3[,2],2),17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=1.75,cex.main=1.75,
               ylab="Population Juvenile Cue Error",main=expression(paste("Juvenile Cue Weight")),color=heat.colors,family="sans")
points(c(-0.5,2.85),c(1,1),type="l")
points(c(-0.5,2.85),c(3,3),type="l")
legend(c(0.025*4,0.1*4),c(0.975*4,0.85*4),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##---------------------------------------------------------------------##
##
## The section of code plots the cue weights as a function of environmental variance
##
## The plots include Figure 3b and 3d.
##
##---------------------------------------------------------------------##

#################################################
# Code for Figure 3d
#################################################

# Set vector of environmental variances
VE<-seq(0.25,4,0.25)

# Generate matrix to store optimal cue weights
CueOpt<-matrix(NA,length(VE),5)
colnames(CueOpt)<-c("mu","wj","wme","wmp","wr")

# For each value of environmental variance...
for(i in 1:length(VE)){
# ... set environmental parameters.
  Parms<-list(l=0.7,K=exp(1),ve=VE[i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=1,yme=0.25,muj=0,mume=0,mump=0)
# Solve for optimal cue weights
  Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                 control=list(fnscale=-1),lower=c(-10,-10,-10,-10,0),
                 upper=c(10,10,10,10,10))
  CueOpt[i,] <- Solve$par
}

# The reference phenotype is zero, so replace it with total cue weight for the purpose
# of visualization.
CueOpt[,1]<-rowSums(CueOpt[,2:4])

#Generate Figure 3d
pdf(file="Figure3d.pdf",width = (400/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
matplot(VE,CueOpt,type="l",ylab="Cue Weights",xlab="Environment Variance",lwd=2,ylim=c(0,1),cex.lab=1.75,cex.main=1.75)
legend(c(0.43*max(VE),0.7*max(VE)),c(0.65,0.3),
       c("Total Plasticity","Juvenile Cue","Mat. Env. Cue","Mat. Pheno. Cue","Randomization"),
       col=1:5,lty=1:5,lwd=2,y.intersp=0.75,cex=1.25,bty="n")
title(expression(paste((sigma[J]^{pop})^{2}," = 1")),cex.main=1.75)
legend(c(0.05*max(VE),0.15*max(VE)),c(1,0.85),
       c("D"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

#################################################
# Code for Figure 3b
#################################################

# For each value of environmental variance...
for(i in 1:length(VE)){
# ... set environmental parameters. Note that yj is now 3 rather than 1.
  Parms<-list(l=0.7,K=exp(1),ve=VE[i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=3,yme=0.25,muj=0,mume=0,mump=0)
# Solve for optimal cue weights
  Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                 control=list(fnscale=-1),lower=c(-10,-10,-10,-10,0),
                 upper=c(10,10,10,10,10))
  CueOpt[i,] <- Solve$par
}

# Calculate total cue weight
CueOpt[,1]<-rowSums(CueOpt[,2:4])

#Generate Figure 3b
pdf(file="Figure3b.pdf",width = (400/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
matplot(VE,CueOpt,type="l",ylab="Cue Weights",xlab="Environment Variance",lwd=2,ylim=c(0,1),cex.lab=1.75,cex.main=1.75,family="sans")
legend(c(0.2*max(VE),0.4*max(VE)),c(1.0,0.5),
       c("Total Plasticity","Juvenile Cue","","Mat. Env. Cue","Mat. Phen. Cue","Randomization"),
       col=c(1,2,0,3,4,5),lty=c(1,2,0,3,4,5),lwd=2,y.intersp=0.75,cex=1.25,bty="n")
title(expression(paste((sigma[J]^{pop})^{2}," = 3")),cex.main=1.75)
legend(c(0.05*max(VE),0.15*max(VE)),c(1,0.85),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##---------------------------------------------------------------------##
##
## This section of code constructs a series of contour plots on 
## the axes of autocorrelation (L) and population-level juvenile cue error (YJ).
##
## The counter plots include Figure 4a and 4c, all panels of Figure 6, 
## and Figure S2c and S2d.
##
##---------------------------------------------------------------------##

# Set grid of environmental parameters. Note that environmental variant (VE) 
# changes with autocorrelation (Lam) in order to keep the stationary variance 
# of the environmental process constant.
Lam<-seq(0,0.95,0.05)
VE<-(1-Lam^2)/(1-0.7^2)
ParMat4<-matrix(NA,3,17*20)
ParMat4[1,]<-rep(seq(0,4,0.25),each=20)
ParMat4[2,]<-rep(Lam,17)
ParMat4[3,]<-rep(VE,17)

# Initialize matrices for fitness and optimal cue weights
Fit<-matrix(NA,17*20,7)
CueOpt3<-matrix(NA,17*20,5)
colnames(CueOpt3)<-c("mu","wj","wme","wmp","wr")
TempPar<-matrix(NA,10,5)
TempVal<-matrix(NA,10,1)

set.seed(2)
# For each set of parameters in the parameter grid...
for(i in 1:length(ParMat4[1,])){
  # ... set the parameters.
  Parms<-list(l=ParMat4[2,i],K=exp(1),ve=ParMat4[3,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  # Generate 10 initial conditions for the optimizer, using Latin Hypercube Sampling.
  Inits<-GenLHS(matrix(c(-10,-10,-10,-1,0,10,10,10,10,10),5,2),10)
  colnames(Inits)<-c("mu","wj","wme","wmp","wr")
  # For each set of initial conditions...
  for(j in 1:dim(Inits)[1]){
    # ...solve for the optimal cue weights given the environmental parameters
    Solve <- optim(par=Inits[j,],fn=Fitness,Parms=Parms,method="L-BFGS-B",
                   control=list(fnscale=-1),lower=c(-10,-10,-10,-1,0),
                   upper=c(10,10,10,10,10))
    TempPar[j,]<-Solve$par
    TempVal[j]<-Solve$value
  }
  # If the initial conditions produce different optimal cue weights, 
  # print optimal cue weights for all initial conditions.
  if(sum((t(TempPar)-colMeans(TempPar))^2)>0.0001){
    print(TempPar)
  }
  # Use the set of cue weights that generates the highest fitness value
  CueOpt3[i,] <- TempPar[which.max(TempVal),]
  
  # Correct for floating point errors in mu.
  CueOpt3[,1]<-round(CueOpt3[,1],7)
  
  # Record fitness in pre-change environment under optimal cue weights.
  Fit[i,1]<-max(TempVal)
  
  # Calculate fitness after increase in mean environment by 1
  Parms<-list(l=ParMat4[2,i],K=exp(1),ve=ParMat4[3,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=1,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  Temp1 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,5]<-Temp1$ESL
  
  # Calculate fitness after increase in environmental variance by 0.4
  Parms<-list(l=ParMat4[2,i],K=exp(1),ve=ParMat4[3,i]+0.4,vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  Temp2 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,2] <- (Fit[i,1]-Temp2$Fit)
  
  # Calculate fitness after decrease in autocorrelation by 0.1. When autocorrelation is 0, record fitness without change.
  if(ParMat4[2,i]==0){
    Fit[i,3] <- Fit[i,1]
  }  else {
    Parms<-list(l=ParMat4[2,i]-0.1,K=exp(1),ve=ParMat4[3,i]*(1-(ParMat4[2,i]-0.1)^2)/(1-ParMat4[2,i]^2),vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
    Temp3 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
    Fit[i,3] <- (Fit[i,1]-Temp3$Fit)
  }
  
  # Calculate fitness after increase in population-level juvenile cue error by 0.4
  Parms<-list(l=ParMat4[2,i],K=exp(1),ve=ParMat4[3,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i]+0.4,yme=0.25,muj=0,mume=0,mump=0)
  Temp4 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,4] <- (Fit[i,1]-Temp4$Fit)
  
  # Calculate fitness after increase in individual-level juvenile cue error by 0.4
  Parms<-list(l=ParMat4[2,i],K=exp(1),ve=ParMat4[3,i],vj=0.25+0.4,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  Temp5 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,6] <- (Fit[i,1]-Temp5$Fit)
  
  # Calculate fitness after increase in juvenile cue bias by 1
  Parms<-list(l=ParMat4[2,i],K=exp(1),ve=ParMat4[3,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=1,mume=0,mump=0)
  Temp6 <- Fitness(CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,7] <- (Fit[i,1]-Temp6$Fit)
  
  print(i)
}

# Create Figure 4c
pdf(file="Figure4c.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,1],20,17)
filled.contour(Lam,seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=1.75,cex.main=1.75,
               ylab="Population Juvenile Cue Error",main=expression(paste("Per-Capita Growth Rate")),color=topo.colors,family="sans")
legend(c(-0.05,0.025),c(0.975*4,0.85*4),
       c("C"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Create Figure 6a
pdf(file="Figure6a.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,5],20,17)
filled.contour(Lam,seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to change in ",mu[theta])),color=topo.colors,family="sans")
legend(c(-0.05,0.025),c(0.975*4,0.85*4),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Create Figure 6b
pdf(file="Figure6b.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,7],20,17)
filled.contour(Lam,seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to change in ",mu[J])),color=topo.colors,family="sans")
legend(c(-0.05,0.025),c(0.975*4,0.85*4),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Create Figure 6c
pdf(file="Figure6c.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,2],20,17)
filled.contour(seq(0,0.95,0.05),seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to increase in ",sigma[theta]^{2})),color=topo.colors,family="sans")
legend(c(-0.05,0.025),c(0.975*4,0.85*4),
       c("C"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Create Figure 6d
pdf(file="Figure6d.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,3],20,17)
FitChange<-FitChange[-c(1,2),]
filled.contour(seq(0.1,0.95,0.05),seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to decrease in ",lambda)),color=topo.colors,family="sans")
legend(c(0.05,0.12),c(0.975*4,0.85*4),
       c("D"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Create Figure S2d
pdf(file="FigureS2d.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,4],20,17)
filled.contour(Lam,seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to increase in ",(sigma[J]^{pop})^{2})),color=topo.colors,family="sans")
legend(c(-0.05,0.025),c(0.975*4,0.85*4),
       c("D"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Create Figure S2c
pdf(file="FigureS2c.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(Fit[,6],20,17)
Borked<-ifelse(FitChange<0,-1000,FitChange)
filled.contour(Lam,seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=1.75,cex.main=1.5,
               ylab="Population Juvenile Cue Error",main=expression(paste("Fitness loss due to increase in ",(sigma[J]^{ind})^{2})),color=topo.colors,family="sans")
legend(c(-0.05,0.025),c(0.975*4,0.85*4),
       c("C"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Create Figure 4a
pdf(file="Figure4a.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix(CueOpt3[,2],20,17)
filled.contour(Lam,seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=1.75,cex.main=1.75,
               ylab="Population Juvenile Cue Error",main=expression(paste("Juvenile Cue Weight")),color=heat.colors,family="sans")
points(c(-0.5,0.657),c(1,1),type="l")
points(c(-0.5,0.657),c(3,3),type="l")
legend(c(-0.05,0.025),c(0.975*4,0.85*4),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##---------------------------------------------------------------------##
##
## The section of code plots the cue weights as a function of autocorrelation
##
## The plots include Figure 4b and 4d.
##
##---------------------------------------------------------------------##

#################################################
# Code for Figure 4d
#################################################

# Set vectors of autocorrelation and environmental variances
Lam<-seq(0.975,0,-0.025)
VE<-(1-Lam^2)/(1-0.7^2)

# Initialize matrix for optimal cue weights
CueOpt<-matrix(NA,length(Lam),5)
colnames(CueOpt)<-c("mu","wj","wme","wmp","wr")

# For each value of autocorrelation...
for(i in 1:length(Lam)){
  # ... set environmental parameters.
  Parms<-list(l=Lam[i],K=exp(1),ve=VE[i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=1,yme=0.25,muj=0,mume=0,mump=0)
  # Solve for optimal cue weights 
  Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                 control=list(fnscale=-1),lower=c(-10,-10,-10,-1,0),
                 upper=c(10,10,10,10,10))
  CueOpt[i,] <- Solve$par
}

# Calculate total cue weight.
CueOpt[,1]<-rowSums(CueOpt[,2:4])

# Generate Figure 4d
pdf(file="Figure4d.pdf",width = (400/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
matplot(Lam,CueOpt,type="l",ylab="Cue Weights",xlab="Autocorrelation",lwd=2,ylim=c(0,1),cex.lab=1.75,family="sans")
legend(c(0.05*max(Lam),0.25*max(Lam)),c(0.55,0.35),
       c("Total Plasticity","Juvenile Cue","Mat. Env. Cue","Mat. Phen. Cue","Randomization"),
       col=1:5,lty=1:5,lwd=2,y.intersp=0.75,cex=1.25,bty="n")
title(expression(paste((sigma[J]^{pop})^{2}," = 1")),cex.main=1.75)
legend(c(0*max(Lam),0.1*max(Lam)),c(1,0.85),
       c("D"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

#################################################
# Code for Figure 4b
#################################################

# Set vectors of autocorrelation and environmental variances
Lam<-seq(0.975,0,-0.025)
VE<-(1-Lam^2)/(1-0.7^2)

# Initialize matrix for optimal cue weights
CueOpt<-matrix(NA,length(Lam),5)
colnames(CueOpt)<-c("mu","wj","wme","wmp","wr")

# For each value of autocorrelation...
for(i in 1:length(Lam)){
  # ... set environmental parameters.
  Parms<-list(l=Lam[i],K=exp(1),ve=VE[i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=3,yme=0.25,muj=0,mume=0,mump=0)
  # Solve for optimal cue weights 
  Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=1),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                 control=list(fnscale=-1),lower=c(0,-10,-10,-1,0),
                 upper=c(10,10,10,10,10))
  CueOpt[i,] <- Solve$par
}

# Display total cue weight.
CueOpt[,1]<-rowSums(CueOpt[,2:4])

# Generate Figure 4b
pdf(file="Figure4b.pdf",width = (400/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
matplot(Lam,CueOpt,type="l",ylab="Cue Weights",xlab="Autocorrelation",lwd=2,ylim=c(0,1),cex.lab=1.75,family="sans")
legend(c(0.1*max(Lam),0.45*max(Lam)),c(1.07,0.67),x.intersp=0.4,
       c("Total Plasticity","Juvenile Cue","Mat. Env. Cue","Mat. Phen. Cue","Randomization"),
       col=1:5,lty=1:5,lwd=2,y.intersp=0.8,cex=1.25,bty="n")
title(expression(paste((sigma[J]^{pop})^{2}," = 3")),cex.main=1.5)
legend(c(0*max(Lam),0.1*max(Lam)),c(1,0.85),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##---------------------------------------------------------------------##
##
## Code for Figure S1: Comparison of how the optimal cue weights
## change as a function of individual vs population cue error.
##
##---------------------------------------------------------------------##

#################################################
# Code for Figure S1b: Population-level cue error
#################################################

# Set vector of population-level juvenile cue error
YJ<-seq(0.25,4,0.25)

# Initialize matrix for optimal cue weights
CueOpt<-matrix(NA,length(YJ),5)
colnames(CueOpt)<-c("mu","wj","wme","wmp","wr")

# For each value of population-level juvenile cue error...
for(i in 1:length(YJ)){
# ... set environmental parameters.
  Parms<-list(l=0.7,K=exp(1),ve=2,vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=YJ[i],yme=0.25,muj=0,mume=0,mump=0)
# Solve for optimal cue weights.
  Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                 control=list(fnscale=-1),lower=c(-10,-10,-10,-1,0),
                 upper=c(10,10,10,10,10))
  CueOpt[i,] <- Solve$par
}

# Calculate total cue weight.
CueOpt[,1]<-rowSums(CueOpt[,2:4])

# Generate Figure S1b
pdf(file="FigureS1b.pdf",width = (400/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
matplot(YJ,CueOpt,type="l",ylab="Cue Weights",xlab="Population Juvenile Cue Error",lwd=2,ylim=c(0,1),cex.lab=1.75,cex.main=1.75,family="sans")
legend(c(0.3*max(YJ),0.5*max(YJ)),c(1.025,0.675),
       c("Sum of Cue Weights","Juvenile Cue","","Mat. Env. Cue","Mat. Phen. Cue","Randomization"),
       col=c(1,2,0,3,4,5),lty=c(1,2,0,3,4,5),lwd=2,y.intersp=0.75,cex=1.25,bty="n")
legend(c(0.05*max(YJ),0.15*max(YJ)),c(1,0.85),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

#################################################
# Code for Figure S1a: Individual-level cue error
#################################################

# Set vector of individual-level juvenile cue error
VJ<-seq(0.25,4,0.25)

# Initialize matrix for optimal cue weights
CueOpt<-matrix(NA,length(VJ),5)
colnames(CueOpt)<-c("mu","wj","wme","wmp","wr")

# For each value of individual-level juvenile cue error...
for(i in 1:length(VJ)){
# ... set environmental parameters.
  Parms<-list(l=0.7,K=exp(1),ve=2,vj=VJ[i],vme=0.25,vmp=0.25,mutheta=0,yj=0.25,yme=0.25,muj=0,mume=0,mump=0)
# Solve for optimal cue weights.
  Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                 control=list(fnscale=-1),lower=c(-10,-10,-10,-1,0),
                 upper=c(10,10,10,10,10))
  CueOpt[i,] <- Solve$par
}

# Calculate total cue weight.
CueOpt[,1]<-rowSums(CueOpt[,2:4])

# Generate Figure S1a
pdf(file="FigureS1a.pdf",width = (400/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
matplot(VJ,CueOpt,type="l",ylab="Cue Weights",xlab="Individual Juvenile Cue Error",lwd=2,ylim=c(0,1),cex.lab=1.75,cex.main=1.75,family="sans")
legend(c(0.3*max(YJ),0.5*max(YJ)),c(1.025,0.675),
       c("Sum of Cue Weights","Juvenile Cue","","Mat. Env. Cue","Mat. Phen. Cue","Randomization"),
       col=c(1,2,0,3,4,5),lty=c(1,2,0,3,4,5),lwd=2,y.intersp=0.75,cex=1.25,bty="n")
legend(c(0.05*max(YJ),0.15*max(YJ)),c(1,0.85),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##-------------------------------------------------------------------------##
##-------------------------------------------------------------------------##
##
## Code for sensitivity analysis and figures contained in Supplements 4 & 5 
##
##-------------------------------------------------------------------------##
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
##
## Sensitivity with environmental variance (VE) as the independent variable
##
## This code generates Figure S4a and S4b, and Figure S6c
##
##-------------------------------------------------------------------------##

# Latin Hypercube Sample environmental parameters, except for VE
set.seed(10)
l<-c(0,1); ve<-c(0,4); vj<-c(0,4) 
vme<-c(0,4); vmp<-c(0,4); yj<-c(0,4); yme<-c(0,4)
ParmLim<-rbind(l,ve,vj,vme,vmp,yj,yme)
ParmSamp<-GenLHS(ParmLim,100)

l<-ParmSamp[,1]; vj<-ParmSamp[,3] 
vme<-ParmSamp[,4]; vmp<-ParmSamp[,5]; yj<-ParmSamp[,6]; yme<-ParmSamp[,7]

# For each set of sample parameters, calculate optimal cue weights along a gradient of VE
ve<-seq(0.1,4,length=50)

#Store cue optimum and fitness values
CueOpt11<-matrix(NA,100*50,5)
colnames(CueOpt11)<-c("mu","wj","wme","wmp","wr")
Fit11<-matrix(NA,100*50,6)

# For each Latin Hyper-cube Sample of parameters...
for(i in 1:100){
  # and for each value of environmental variance...
  for(j in 1:length(ve)){
    # ... set environmental parameters.
    Parms<-list(l=l[i],K=exp(1),ve=ve[j],vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[i],yme=yme[i],muj=0,mume=0,mump=0)
    # Solve for optimal cue weights.
    Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                   control=list(fnscale=-1),lower=c(-10,-10,-10,-10,0),
                   upper=c(10,10,10,10,10))
    CueOpt11[(i-1)*50+j,] <- Solve$par
    
    #### Correct for floating point errors in mu.
    CueOpt11[(i-1)*50+j,1]<-round(CueOpt11[(i-1)*50+j,1],7)
    ##########
    
    # Record fitness value for environmental parameters
    Fit11[(i-1)*50+j,1] <- Fitness(CueOpt11[(i-1)*50+j,],Parms=Parms,Info=F)
    
    # Calculate fitness loss after an increase in environmental variance by 0.4
    Parms<-list(l=l[i],K=exp(1),ve=ve[j]+0.4,vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[i],yme=yme[i],muj=0,mume=0,mump=0)
    Temp2 <- Fitness(CueOpt11[(i-1)*50+j,],Parms=Parms,Info=T)
    Fit11[(i-1)*50+j,2] <- (Fit11[(i-1)*50+j,1]-Temp2$Fit)
    
    #print(j)
  }
  print(i)
}

# Format and save results
#ParmMat<-cbind(rep(l,each=50),rep(vj,each=50),rep(vme,each=50),rep(vmp,each=50),rep(yj,each=50),rep(yme,each=50))
#OutPut<-cbind(ParmMat,CueOpt11,Fit11)
#write.csv(OutPut,"Sensitivity_Output_VE_Revised.csv")

# Generate Figure S4b
pdf(file="FigureS4b.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(c(-1,-1),xlim=c(0,4),ylim=c(0,1),xlab="Environmental Variance",ylab="TGP Cue Weights",cex.lab=2,family="sans")
# Plot values with autocorrelation < 0.5 as red and > 0.5 as black
for(i in 1:100){
  if(l[i]<0.5){
    Red<-2
  }else{
    Red<-1
  }
  points(ve,CueOpt11[((i-1)*50+1):(i*50),3]+CueOpt11[((i-1)*50+1):(i*50),4],type="l",col=Red)
}
legend(c(1,1),c(0.6,1),c("Autocorrelation > 0.5","Autocorrelation < 0.5"),col=c(1,2),lwd=1,bty="n",cex=1.5,y.intersp=0.75,x.intersp=0.75)
legend(c(0,0.08*4),c(0.975*1,0.825*1),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()


# Generate Figure S4a
pdf(file="FigureS4a.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(c(-1,-1),xlim=c(0,4),ylim=c(0,1),xlab="Environmental Variance",ylab="Juvenile Cue Weight",cex.lab=2,family="sans")
# Plot values with population-level juvenile cue error > 2 as red and < 2 as black
for(i in 1:100){
  if(yj[i]>2){
    Red<-2
  }else{
    Red<-1
  }
  points(ve,CueOpt11[((i-1)*50+1):(i*50),2],type="l",col=Red)
}
legend(c(1,4),c(0.05,0.25),c("Juvenile Cue Error < 2","Juvenile Cue Error > 2"),col=c(1,2),lwd=1,bty="n",cex=1.5,y.intersp=0.75,x.intersp=0.75)
legend(c(0,0.08*4),c(0.975*1,0.825*1),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Generate Figure S6c
pdf(file="FigureS6c.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(c(-1,-1),xlim=c(0,1),ylim=c(0,0.7),xlab=expression(paste("Juvenile Cue Weight (",w[J],")")),ylab=expression(paste("Fitness loss from increase in ",sigma[theta]^{2})),cex.lab=1.75,cex.main=1.25,
     main=expression(paste("Fitness loss and ",w[J] ," as a function of pre-change ",sigma[theta]^{2})),family="sans")
for(i in 1:100){
  points(CueOpt11[((i-1)*50+1):(i*50),2],Fit11[((i-1)*50+1):(i*50),2],type="l",col=1)
}
legend(c(0,0.08),c(0.975*0.7,0.825*0.7),
       c("C"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##-------------------------------------------------------------------------##
##
## Sensitivity with autocorrelation (L) as the independent variable 
##
## This code generates Figure S6a and S6d
##
##-------------------------------------------------------------------------##

# Latin Hypercube Sample environmental parameters, except for L
set.seed(11)
l<-c(0,1);ve<-c(0,4); vj<-c(0,4) 
vme<-c(0,4); vmp<-c(0,4); yj<-c(0,4); yme<-c(0,4)
ParmLim<-rbind(l,ve,vj,vme,vmp,yj,yme)
ParmSamp<-GenLHS(ParmLim,100)

ve<-ParmSamp[,2];  vj<-ParmSamp[,3] 
vme<-ParmSamp[,4]; vmp<-ParmSamp[,5]; yj<-ParmSamp[,6]; yme<-ParmSamp[,7]

# For each set of sample parameters, calculate optimal cue weights along a gradient of L

l<-seq(0.1,0.98,length=50)

#Store cue optimum and fitness values
CueOpt12<-matrix(NA,100*50,5)
colnames(CueOpt12)<-c("mu","wj","wme","wmp","wr")
Fit12<-matrix(NA,100*50,6)

# For each Latin Hyper-cube Sample of parameters...
for(i in 1:100){
  # for each value of autocorrelation...
  for(j in 1:length(l)){
    # ... set environmental parameters.
    Parms<-list(l=l[j],K=exp(1),ve=ve[i],vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[i],yme=yme[i],muj=0,mume=0,mump=0)
    # Solve for optimal cue weights.
    Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                   control=list(fnscale=-1),lower=c(-10,-10,-10,-1,0),
                   upper=c(10,10,10,10,10))
    CueOpt12[(i-1)*50+j,] <- Solve$par
    
    #### Correct for floating point errors in mu.
    CueOpt12[(i-1)*50+j,1]<-round(CueOpt12[(i-1)*50+j,1],7)
    ##########
    
    # Record fitness value for environmental parameters
    Fit12[(i-1)*50+j,1] <- Fitness(CueOpt12[(i-1)*50+j,],Parms=Parms,Info=F)
    
    # Calculate fitness loss after a decrease in autocorrelation by 0.1
    Parms<-list(l=l[j]-0.1,K=exp(1),ve=ve[i]*(1-(l[j]-0.1)^2)/(1-l[j]^2),vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[i],yme=yme[i],muj=0,mume=0,mump=0)
    Temp2 <- Fitness(CueOpt12[(i-1)*50+j,],Parms=Parms,Info=T)
    Fit12[(i-1)*50+j,2] <- (Fit12[(i-1)*50+j,1]-Temp2$Fit)
    
    # Calculate fitness loss after an increase in environmental mean by 0.1
    Parms<-list(l=l[j],K=exp(1),ve=ve[i],vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=1,yj=yj[i],yme=yme[i],muj=0,mume=0,mump=0)
    Temp2 <- Fitness(CueOpt12[(i-1)*50+j,],Parms=Parms,Info=T)
    Fit12[(i-1)*50+j,3]<-Temp2$ESL
    
    #print(j)
  }
  print(i)
}

# Format and save results
#ParmMat<-cbind(rep(ve,each=50),rep(vj,each=50),rep(vme,each=50),rep(vmp,each=50),rep(yj,each=50),rep(yme,each=50))
#OutPut<-cbind(ParmMat,CueOpt12,Fit12)
#write.csv(OutPut,"Sensitivity_Output_L.csv")

# Generate Figure S6d
pdf(file="FigureS6d.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(c(-1,-1),xlim=c(0,1),ylim=c(0,1.5),xlab=expression(paste("TGP Cue Weights")),ylab=expression(paste("Fitness loss due to decrease in ",lambda)),cex.lab=1.75,cex.main=1.25,
     main=expression(paste("Fitness loss and TGP as a function of pre-change ",lambda)),family="sans")
# Plot values with environment variance < 1 as red and > 1 as black
for(i in 1:100){
  if(ve[i]<1){
    Red<-2
  }else{
    Red<-1
  }
  points(CueOpt12[((i-1)*50+1):(i*50),3]+CueOpt12[((i-1)*50+1):(i*50),4],Fit12[((i-1)*50+1):(i*50),2],type="l",col=Red)
}
legend(c(0.7,1),c(1.35,1.55),c(expression(paste(sigma^{2}," > 1")),expression(paste(sigma^{2}," < 1"))),col=c(1,2),lwd=1,bty="n",cex=1.5,y.intersp=1,x.intersp=0.5)
legend(c(0,0.08),c(0.975*1.5,0.825*1.5),
       c("D"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Generate Figure S6a
pdf(file="FigureS6a.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(c(-1,-1),xlim=c(0,1),ylim=c(0,0.5),xlab="Total Cue Weight",ylab=expression(paste("Fitness loss due to change in ",mu[theta])),cex.lab=1.75,cex.main=1.25,
     main=expression(paste("Fitness loss and plasticity as a function of pre-change ",lambda)),family="sans")
for(i in 1:100){
  points(CueOpt12[((i-1)*50+1):(i*50),2]+CueOpt12[((i-1)*50+1):(i*50),3]+CueOpt12[((i-1)*50+1):(i*50),4],Fit12[((i-1)*50+1):(i*50),3],type="l",col=1)
}
legend(c(0,0.08),c(0.975*0.5,0.825*0.5),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##-------------------------------------------------------------------------##
##
## Sensitivity with Population-level Juvenile Cue Error (YJ) 
## as the independent variable.
##
## This code generates Figure S6b
##
##-------------------------------------------------------------------------##

# Latin Hypercube Sample environmental parameters, except for YJ
set.seed(10)
l<-c(0,1); ve<-c(0,4); vj<-c(0,4) 
vme<-c(0,4); vmp<-c(0,4); yme<-c(0,4); yj<-c(0,4)
ParmLim<-rbind(l,ve,vj,vme,vmp,yj,yme)
ParmSamp<-GenLHS(ParmLim,100)

l<-ParmSamp[,1]; ve<-ParmSamp[,2]; vj<-ParmSamp[,3] 
vme<-ParmSamp[,4]; vmp<-ParmSamp[,5]; yme<-ParmSamp[,7]

# For each set of sample parameters, calculate optimal cue weights along a gradient of yj

yj<-seq(0,4,length=50)

#Store cue optimum and fitness values
CueOpt13<-matrix(NA,100*50,5)
colnames(CueOpt13)<-c("mu","wj","wme","wmp","wr")
Fit13<-matrix(NA,100*50,6)

# For each Latin Hyper-cube Sample of parameters...
for(i in 1:100){
  # for each value of population-level juvenile cue error...
  for(j in 1:length(yj)){
    #... set environmental parameters
    Parms<-list(l=l[i],K=exp(1),ve=ve[i],vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[j],yme=yme[i],muj=0,mume=0,mump=0)
    # Solve for optimal cue weights.
    Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                   control=list(fnscale=-1),lower=c(-10,-10,-10,-1,0),
                   upper=c(10,10,10,10,10))
    CueOpt13[(i-1)*50+j,] <- Solve$par
    
    #### Fix floating point errors in mu.
    CueOpt13[(i-1)*50+j,1]<-round(CueOpt13[(i-1)*50+j,1],7)
    ##########
    
    # Record fitness value for environmental parameters
    Fit13[(i-1)*50+j,1] <- Fitness(CueOpt13[(i-1)*50+j,],Parms=Parms,Info=F)
    
    # Calculate fitness loss after an increase in population level juevnile cue error by 0.4
    Parms<-list(l=l[i],K=exp(1),ve=ve[i],vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[j]+0.4,yme=yme[i],muj=0,mume=0,mump=0)
    Temp2 <- Fitness(CueOpt13[(i-1)*50+j,],Parms=Parms,Info=T)
    Fit13[(i-1)*50+j,2] <- (Fit13[(i-1)*50+j,1]-Temp2$Fit)
    
    # Calculate fitness loss after an increase in juvenile cue bias by 1
    Parms<-list(l=l[i],K=exp(1),ve=ve[i],vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[j],yme=yme[i],muj=1,mume=0,mump=0)
    Temp2 <- Fitness(CueOpt13[(i-1)*50+j,],Parms=Parms,Info=T)
    Fit13[(i-1)*50+j,3] <- (Fit13[(i-1)*50+j,1]-Temp2$Fit)
    
    #print(j)
  }
  print(i)
}

# Format and record data
#ParmMat<-cbind(rep(l,each=50),rep(ve,each=50),rep(vj,each=50),rep(vme,each=50),rep(vmp,each=50),rep(yme,each=50))
#OutPut<-cbind(ParmMat,CueOpt13,Fit13)
#write.csv(OutPut,"Sensitivity_Output_YJ.csv")

# Generate Figure S6b
pdf(file="FigureS6b.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(c(-1,-1),xlim=c(0,1),ylim=c(0,0.5),xlab=expression(paste("Juvenile Cue Weight (",w[J],")")),
     ylab=expression(paste("Fitness loss due to increase in ",mu[J])),cex.lab=1.75,cex.main=1.25,
     main=expression(paste("Fitness loss and ",w[J]," as a function of pre-change ",(sigma[J]^{pop})^{2})),family="sans")
for(i in 1:100){
  points(CueOpt13[((i-1)*50+1):(i*50),2],Fit13[((i-1)*50+1):(i*50),3],type="l",col=1)
}
legend(c(0,0.08),c(0.975*0.5,0.825*0.5),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##-------------------------------------------------------------------------##
##
## Sensitivity in presence of randomization to environmental parameters
##
## This code generates figure S4c
##
##-------------------------------------------------------------------------##

# Latin hypercube sample environmental parameters. Restrict parameter space to 
# low individual-level juvenile cue error.

set.seed(15)
l<-c(0,1); ve<-c(0,4); vj<-c(0,2) 
vme<-c(0,4); vmp<-c(0,4); yj<-c(0,4); yme<-c(0,4)
ParmLim<-rbind(l,ve,vj,vme,vmp,yj,yme)
ParmSamp<-GenLHS(ParmLim,1000)

l<-ParmSamp[,1]; ve<-ParmSamp[,2]; vj<-ParmSamp[,3] 
vme<-ParmSamp[,4]; vmp<-ParmSamp[,5]; yj<-ParmSamp[,6]; yme<-ParmSamp[,7]

#Store cue optimum
CueOpt30<-matrix(NA,1000,5)
colnames(CueOpt11)<-c("mu","wj","wme","wmp","wr")


for(i in 1:1000){
  Parms<-list(l=l[i],K=exp(1),ve=ve[i],vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[i],yme=yme[i],muj=0,mume=0,mump=0)
  Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),fn=Fitness,Parms=Parms,method="L-BFGS-B",
                 control=list(fnscale=-1),lower=c(-10,-10,-10,-10,0),
                 upper=c(10,10,10,10,10))
  CueOpt30[i,] <- Solve$par
  
  #### Fix floating point errors in mu.
  CueOpt30[,1]<-round(CueOpt30[,1],7)
  ##########
  
  print(i)
}

# Create Figure S3d
pdf(file="FigureS3d.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
Red<-ifelse(CueOpt30[,5]>0,0,1)
plot(ve[Red==1],yj[Red==1],col=rgb(0,0,0,0.3),xlab="Environmental Variance",ylab="Population Juvenile Cue Error",cex.lab=2,family="sans")
points(ve[Red==0],yj[Red==0],col=rgb(1,0,0))
legend(c(0,2.5),c(0,1),c("Randomization Present","Randomization Absent"),col=c(2,1),pch=1,bty="o",cex=1.5,y.intersp = 1)
legend(c(0,0.08*4),c(0.975*4,0.825*4),
       c("C"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##---------------------------------------------------------------------##
##
## Tradeoff between TGP and WGP across environments
##
## This code generates Figure S5
##
##---------------------------------------------------------------------##

# Latin hypercube sample environmental parameters.
set.seed(17)
l<-c(0,1); ve<-c(0,4); vj<-c(0,4) 
vme<-c(0,4); vmp<-c(0,4); yj<-c(0,4); yme<-c(0,4)
ParmLim<-rbind(l,ve,vj,vme,vmp,yj,yme)
ParmSamp<-GenLHS(ParmLim,1000)

l<-ParmSamp[,1]; ve<-ParmSamp[,2]; vj<-ParmSamp[,3] 
vme<-ParmSamp[,4]; vmp<-ParmSamp[,5]; yj<-ParmSamp[,6]; yme<-ParmSamp[,7]

#Store cue optimum
CueOpt7<-matrix(NA,1000,5)
colnames(CueOpt7)<-c("mu","wj","wme","wmp","wr")

for(i in 1:1000){
  Parms<-list(l=l[i],K=exp(1),ve=ve[i],vj=vj[i],vme=vme[i],vmp=vmp[i],mutheta=0,yj=yj[i],yme=yme[i],muj=0,mume=0,mump=0)
  Solve <- optim(par=c(mu=0,wj=0,wme=0,wmp=0,wr=0),Parms=Parms,fn=Fitness,method="L-BFGS-B",
                 control=list(fnscale=-1),lower=c(-10,-10,-10,-10,0),
                 upper=c(10,10,10,10,10))
  CueOpt7[i,] <- Solve$par
  print(i)
}

# Core Result: Tradeoff across environments between WGP and TGP cues so long as environmental variance 
# and autocorrelation are not too low. 

# Partition parameters space into low autocorrelation, high autocorrelation + low environmental 
# variance, and high autocorrelation + high environmental variance. Low autocorrelation is not plotted below.
CueOpt7Low<-CueOpt7[(ve<1)&(l>0.5),]
CueOpt7High<-CueOpt7[ve>=1&(l>0.5),]
CueOpt7Med<-CueOpt7[l<=0.5,]


# Generate Figure S5.
pdf(file="FigureS5.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
plot(CueOpt7High[,2],CueOpt7High[,3]+CueOpt7High[,4],xlab="Within-Generational Plasticity",ylab="Transgenerational Plasticity",xlim=c(0,1),ylim=c(0,1),cex.lab=2,family="sans")
points(CueOpt7Low[,2],CueOpt7Low[,3]+CueOpt7Low[,4],xlab="WGP",ylab="TGP",col=2)
legend(c(0.16,1),c(0.85,1.05),c("Environmental Variance > 1","Environmental Variance < 1"),col=c(1,2),pch=1,bty="n",cex=1.5,y.intersp=1)
dev.off()

##-------------------------------------------------------------------------##
##
## Exploration of the nonlinearity in the phenotypic variance and in the 
## environment fluctuation load caused by the maternal phenotype cue weight
## 
## This code generates Figure S3
##
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
# Code for Figures S3a and S3c
##-------------------------------------------------------------------------##


# Set grid of environmental parameters
ParMat4<-matrix(NA,2,17^2)
ParMat4[1,]<-rep(seq(0,4,0.25),each=17)
ParMat4[2,]<-rep(seq(0,4,0.25),17)

# Initialize matrices for data storage
Fit<-matrix(NA,17^2,7)
Fit2<-matrix(NA,17^2,4)
FitMod2<-matrix(NA,17^2,7)
CueOpt3<-matrix(NA,17^2,5)
colnames(CueOpt3)<-c("mu","wj","wme","wmp","wr")
TempPar<-matrix(NA,10,5)
TempVal<-matrix(NA,10,1)
CueOptMod2<-matrix(NA,17*20,5)
colnames(CueOptMod2)<-c("mu","wj","wme","wmp","wr")

set.seed(1)
# For each set of parameters in the parameter grid...
for(i in 1:length(ParMat4[1,])){
  # ... set the parameters.
  Parms<-list(l=0.7,K=exp(1),ve=ParMat4[2,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0.25,muj=0,mume=0,mump=0)
  # Generate 10 initial conditions for the optimizer, using Latin Hypercube Sampling
  Inits<-GenLHS(matrix(c(0,0,0,0,0,10,10,10,10,10),5,2),10)
  colnames(Inits)<-c("mu","wj","wme","wmp","wr")
  # For each set of initial conditions...
  for(j in 1:dim(Inits)[1]){
    # ...find the optimal cue weights for the given environmental parameters.
    Solve <- optim(par=Inits[j,],fn=Fitness,Parms=Parms,method="L-BFGS-B",
                   control=list(fnscale=-1),lower=c(-10,-10,-10,-10,0),
                   upper=c(10,10,10,10,10))
    TempPar[j,]<-Solve$par
    TempVal[j]<-Solve$value
  }
  # If different intial conditions produce different optimal cue weights...
  if(sum((t(TempPar)-colMeans(TempPar))^2)>0.0001){
    # ...print the optimal cue weights for each initial condition.
    print(TempPar)
  }
  # Use the set of cue weights that correspond to the highest fitness value.
  CueOpt3[i,] <- TempPar[which.max(TempVal),]
  
  # Correct for floating point errors in mu.
  CueOpt3[,1]<-round(CueOpt3[,1],7)
  
  # Calculate and record each component of fitness
  
  Temp1<-Fitness(CueWeight=CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,1]<-Temp1$VX
  Fit[i,2]<-Temp1$EFL
  Fit[i,3]<-Temp1$PEL
  Fit[i,4]<-Temp1$Fit

  print(i)
}


# Calculate the approximate phenotypic variance, environmental fluctuation
# load, and fitness assuming that the maternal environment and phenotype
# cues are equivalent.

ApproxVX<-CueOpt3[,2]^2*0.25+CueOpt3[,3]^2*0.25+CueOpt3[,4]^2*0.25+CueOpt3[,5]^2
ApproxEFL<-ParMat4[2,]/(1-0.7^2)*((1-CueOpt3[,2])^2+(CueOpt3[,3]+CueOpt3[,4])^2-2*0.7*(1-CueOpt3[,2])*(CueOpt3[,3]+CueOpt3[,4]))/(2*(1+ApproxVX))
ApproxFit<-1-1/2*log(1+ApproxVX)-ApproxEFL-CueOpt3[,2]^2*ParMat4[1,]/(2*(1+ApproxVX))

# Generate Figure S3a
pdf(file="FigureS3a.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix((ApproxVX-Fit[,1])/Fit[,1],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=2,cex.main=2,
               ylab="Population Juvenile Cue Error",main=expression(paste("Approximation Error in ",sigma[X]^{2})),color=topo.colors,family="sans")
legend(c(0.02*4,0.09*4),c(0.975*4,0.85*4),
       c("A"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Generate Figure S3c
pdf(file="FigureS3c.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix((ApproxEFL-Fit[,2])/Fit[,2],17,17)
filled.contour(seq(0.25,4,0.25),seq(0,4,0.25),FitChange[-1,],nlevels=20,xlab="Environmental Variance",cex.lab=2,cex.main=2,
               ylab="Population Juvenile Cue Error",main=expression(paste(main="Approximation Error in ",F)),color=topo.colors,family="sans")
legend(c(0.02*4,0.09*4),c(0.975*4,0.85*4),
       c("C"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

##-------------------------------------------------------------------------##
# Code for Figures S3b and S3d
##-------------------------------------------------------------------------##

# Set grid of environmental parameters. Note that environmental variant (VE) 
# changes with autocorrelation (Lam) in order to keep the stationary variance 
# of the environmental process constant.
Lam<-seq(0,0.95,0.05)
VE<-(1-Lam^2)/(1-0.7^2)
ParMat4<-matrix(NA,3,17*20)
ParMat4[1,]<-rep(seq(0,4,0.25),each=20)
ParMat4[2,]<-rep(Lam,17)
ParMat4[3,]<-rep(VE,17)

# Initialize matrices for fitness and optimal cue weights
Fit<-matrix(NA,17*20,7)
FitMod<-matrix(NA,17*20,7)
FitMod2<-matrix(NA,17*20,7)
CueOpt3<-matrix(NA,17*20,5)
colnames(CueOpt3)<-c("mu","wj","wme","wmp","wr")
TempPar<-matrix(NA,10,5)
TempVal<-matrix(NA,10,1)

CueOptMod<-matrix(NA,17*20,5)
colnames(CueOptMod)<-c("mu","wj","wme","wmp","wr")

CueOptMod2<-matrix(NA,17*20,5)
colnames(CueOptMod2)<-c("mu","wj","wme","wmp","wr")


set.seed(2)
# For each set of parameters in the parameter grid...
for(i in 1:length(ParMat4[1,])){
  # ... set the parameters.
  Parms<-list(l=ParMat4[2,i],K=exp(1),ve=ParMat4[3,i],vj=0.25,vme=0.25,vmp=0.25,mutheta=0,yj=ParMat4[1,i],yme=0,muj=0,mume=0,mump=0)
  # Generate 10 initial conditions for the optimizer, using Latin Hypercube Sampling.
  Inits<-GenLHS(matrix(c(0,0,0,0,0,10,10,10,10,10),5,2),10)
  colnames(Inits)<-c("mu","wj","wme","wmp","wr")
  # For each set of initial conditions...
  for(j in 1:dim(Inits)[1]){
    # ...solve for the optimal cue weights given the environmental parameters
    Solve <- optim(par=Inits[j,],fn=Fitness,Parms=Parms,method="L-BFGS-B",
                   control=list(fnscale=-1),lower=c(-10,-10,-10,-10,0),
                   upper=c(10,10,10,10,10))
    TempPar[j,]<-Solve$par
    TempVal[j]<-Solve$value
  }
  # If the initial conditions produce different optimal cue weights, 
  # print optimal cue weights for all initial conditions.
  if(sum((t(TempPar)-colMeans(TempPar))^2)>0.0001){
    print(TempPar)
  }
  # Use the set of cue weights that generates the highest fitness value
  CueOpt3[i,] <- TempPar[which.max(TempVal),]
  
  # Correct for floating point errors in mu.
  CueOpt3[,1]<-round(CueOpt3[,1],7)

  # Calculate and record each component of fitness
  Temp1<-Fitness(CueWeight=CueOpt3[i,],Parms=Parms,Info=T)
  Fit[i,1]<-Temp1$VX
  Fit[i,2]<-Temp1$EFL
  Fit[i,3]<-Temp1$PEL
  Fit[i,4]<-Temp1$Fit

  
  print(i)
}

par(mar=c(5.1, 5, 4.1, 2.1))

# Calculate the approximate phenotypic variance, environmental fluctuation
# load, and fitness assuming that the maternal environment and phenotype
# cues are equivalent.

ApproxVX<-CueOpt3[,2]^2*0.25+CueOpt3[,3]^2*0.25+CueOpt3[,4]^2*0.25+CueOpt3[,5]^2
ApproxEFL<-ParMat4[3,]/(1-ParMat4[2,]^2)*((1-CueOpt3[,2])^2+(CueOpt3[,3]+CueOpt3[,4])^2-2*ParMat4[2,]*(1-CueOpt3[,2])*(CueOpt3[,3]+CueOpt3[,4]))/(2*(1+ApproxVX))
ApproxFit<-1-1/2*log(1+ApproxVX)-ApproxEFL-CueOpt3[,2]^2*ParMat4[1,]/(2*(1+ApproxVX))

# Generate Figure S3b
pdf(file="FigureS3b.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix((ApproxVX-Fit[,1])/Fit[,1],20,17)
filled.contour(Lam,seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=2,cex.main=2,
               ylab="Population Juvenile Cue Error",main=expression(paste("Approximation Error in ",sigma[X]^{2})),color=topo.colors,family="sans")
legend(c(-0.03,0.045),c(0.975*4,0.85*4),
       c("B"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()

# Generate Figure S3d
pdf(file="FigureS3d.pdf",width = (475/75),height=(400/75),useDingbats=FALSE)
par(mar=c(5.1, 5, 4.1, 2.1))
FitChange<-matrix((ApproxEFL-Fit[,2])/Fit[,2],20,17)
filled.contour(Lam,seq(0,4,0.25),FitChange,nlevels=20,xlab="Autocorrelation",cex.lab=2,cex.main=2,
               ylab="Population Juvenile Cue Error",main=expression(paste(main="Approximation Error in F")),color=topo.colors,family="sans")
legend(c(-0.03,0.045),c(0.975*4,0.85*4),
       c("D"),col=c(1),y.intersp=0,x.intersp=-0.75, cex=2,bty="o",bg="white")
dev.off()


