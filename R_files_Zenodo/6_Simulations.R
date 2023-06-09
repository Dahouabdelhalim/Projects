#######################################################################################
###This script will run and analyze three sets of simulations in different scenarios###
#######################################################################################

require(phytools)
require(hypervolume)
require(geiger)
require(geometry)
require(TreeSim)
require(caper)

nwdata<-read.csv('WH_CladeData.csv')
time<-na.omit(nwdata$Crown)
nwdata$Species<-log(nwdata$Species) #ln transform species richness
nwdata$Hyp<-nwdata$Hyp^0.25 #Transform hypervolumes
nwdata$Chull<-nwdata$Chull^0.25 #Transform convex hulls too


############################
###Set up the simulations###
#Write function to simulate trees with a minimum of 5 tips#
get.bd.age.tree<-function(pars){
  success<-F
  while(!success){
    try<-sim.bd.age(pars["t"],1,pars["b"],pars["d"],frac=1,mrca=T,complete=FALSE)
    success<-Ntip(try[[1]])>=5
  }
  return(try[[1]])
}

sigsq<-c(0.05108108,0.04915916,0.05852853) #Optimized sigma squared rates from ABC analyses
var.props<-c(76.52,10.71,5.62,2.74)/sum(c(76.52,10.71,5.62,2.74)) #Scaling factor for sigma squared rates

sim.mat<-matrix(nrow=length(time),ncol=4)
colnames(sim.mat)<-c('Age','Tips','TreeLength','CHull')
sim.posterior<-NULL

sim.models<-NULL #Blank list for storing models
############################

############################################
###Run the uniform birth-death simulation###
for(n in 1:1000){
    
    #Generate trees and store tip number and tree length#
    trees<-list(NULL) #Create blank object to store trees
    for (i in 1:length(time)){
      trees[[i]]<-get.bd.age.tree(c(b=0.396805,d=0.2798954,t=time[i]))
      print(paste("Tree",i,"done"))
    }
    sim.mat[,1]<-time
    for (t in 1:length(time)){sim.mat[t,2]<-Ntip(trees[[t]])}
    for (t in 1:length(time)){sim.mat[t,3]<-sum(trees[[t]]$edge.length)}
    
    #Simulate trait data for each tree#
    chulls<-NULL
    for (i in 1:length(time)){
      tree.data<-matrix(nrow=Ntip(trees[[i]]),ncol=4) #Create blank matrix to store tip data for 4 "PC axes"
      for (j in 1:4){ #Simulate the trait data given a set of sigma squared rates
        tree.data[,j]<-fastBM(trees[[i]],a=0,sig2=sigsq[1]*var.props[j],internal=F)
      }
      print(paste("Data for tree",i,"done"))
      chulls[[i]]<-convhulln(tree.data,option='FA')$vol
    }
    
    sim.mat[,4]<-chulls
    
    sim.posterior[[n]]<-sim.mat
    
    #Fit the model and store it#
    sim.models[[n]]<-lm(sim.mat[,4]^0.25~sim.mat[,1]*log(sim.mat[,2]))
}
save.image('UniformBD_Sims.RData')
############################################


#############################################
###Run the declining speciation simulation###
#Rewrite tree simulation function to add K parameter#
get.kbd.age.tree<-function(pars){
  success<-F
  while(!success){
    try<-sim.bd.age(pars["t"],1,pars["b"],0,frac=1,mrca=T,complete=FALSE,K=pars["K"])
    success<-Ntip(try[[1]])>=5
  }
  return(try[[1]])
}

for(n in 1:1000){
    
    #Generate trees and store tip number and tree length#
    trees<-list(NULL) #Create blank object to store trees
    for (i in 1:length(time)){
      trees[[i]]<-get.kbd.age.tree(c(b=0.2655675,t=time[i],K=1205))
      print(paste("Tree",i,"done"))
    }
    sim.mat[,1]<-time
    for (t in 1:length(time)){sim.mat[t,2]<-Ntip(trees[[t]])}
    for (t in 1:length(time)){sim.mat[t,3]<-sum(trees[[t]]$edge.length)}
    
    #Simulate trait data for each tree#
    chulls<-NULL
    for (i in 1:length(time)){
      tree.data<-matrix(nrow=Ntip(trees[[i]]),ncol=4) #Create blank matrix to store tip data for 4 "PC axes"
      for (j in 1:4){ #Simulate the trait data given a set of sigma squared rates
        tree.data[,j]<-fastBM(trees[[i]],a=0,sig2=sigsq[2]*var.props[j],internal=F)
      }
      print(paste("Data for tree",i,"done"))
      chulls[[i]]<-convhulln(tree.data,option='FA')$vol
    }
    
    sim.mat[,4]<-chulls
    
    sim.posterior[[n]]<-sim.mat
    
    #Fit the model and store it#
    sim.models[[n]]<-lm(sim.mat[,4]^0.25~sim.mat[,1]*log(sim.mat[,2]))
}

save.image('K-BD_Sims.RData')
#############################################


################################################
###Run the outlier diversification simulation###
lambdas<-c(0.6894,rep(0.3968050,6),0.4997,rep(0.3968050,8),0.2165,rep(0.3968050,9))

for(n in 1:1000){
    
    #Generate trees and store tip number and tree length#
    trees<-list(NULL) #Create blank object to store trees
    #for(i in 1:1){
    for (i in 1:length(time)){
      #trees[[i]]<-pbtree(b=lambdas[i],d=0.2798954,t=time[i]) ###GENERATING THE TREES IS FAILING WITH EMBERIZOIDEA'S RATE - TREE IS UNREASONABLY LARGE AND NEVER FULLY CALCULATES
      trees[[i]]<-get.kbd.age.tree(c(b=lambdas[i],d=0.2798954,t=time[i],K=4000))
      print(paste("Tree",i,"done"))
    }
    sim.mat[,1]<-time
    for (t in 1:length(time)){sim.mat[t,2]<-Ntip(trees[[t]])}
    for (t in 1:length(time)){sim.mat[t,3]<-sum(trees[[t]]$edge.length)}
    
    #Simulate trait data for each tree#
    chulls<-NULL
    for (i in 1:length(time)){
      tree.data<-matrix(nrow=Ntip(trees[[i]]),ncol=4) #Create blank matrix to store tip data for 4 "PC axes"
      for (j in 1:4){ #Simulate the trait data given a set of sigma squared rates
        tree.data[,j]<-fastBM(trees[[i]],a=0,sig2=sigsq[3]*var.props[j],internal=F)
      }
      print(paste("Data for tree",i,"done"))
      chulls[[i]]<-convhulln(tree.data,option='FA')$vol
    }
    
    sim.mat[,4]<-chulls
    
    sim.posterior[[n]]<-sim.mat
    
    #Fit the model and store it#
    sim.models[[n]]<-lm(sim.mat[,4]^0.25~sim.mat[,1]*log(sim.mat[,2]))
    
}
  
save.image('OutlierBD_Sims.RData')
################################################

#####################################################################
###Compare the real models to the simulated distribution of models###
backbone<-read.tree('PGLS_Backbone.tre') #Read in backbone tree with 1 tip/clade
comp.data<-comparative.data(backbone,nwdata,'Lineage')

real.mod<-pgls(Hyp~Crown*Species,data=comp.data) #Generate best-fitting model
summary(real.mod)

real.vec<-vector(length=4,mode='numeric')
names(real.vec)<-c('Crown','Species','Interaction','RSq')
real.vec[1]<-real.mod$model$coef[2]
real.vec[2]<-real.mod$model$coef[3]
real.vec[3]<-real.mod$model$coef[4]
real.vec[4]<-summary(real.mod)$adj.r.squared
real.vec

###Uniform Birth-death###
model.mat<-matrix(nrow=1000,ncol=4)
colnames(model.mat)<-names(real.vec)

load('UniformBD_Sims.RData')

for(n in 1:1000){
  model.mat[n,1]<-sim.models[[n]]$coefficients[[2]]
  model.mat[n,2]<-sim.models[[n]]$coefficients[[3]]
  model.mat[n,3]<-sim.models[[n]]$coefficients[[4]]
  model.mat[n,4]<-summary(sim.models[[n]])$adj.r.squared
}

###Calculate test statistics for differences between real and simulated model parameters###
#Crown parameter#
length(model.mat[,1][model.mat[,1]<(median(model.mat[,1])-abs(median(model.mat[,1])-real.vec[1]))|
                       model.mat[,1]>(median(model.mat[,1])+abs(median(model.mat[,1])-real.vec[1]))])/1000
#Calculate proportion of obvs on either end of distribution more extreme than difference of median and observed parameter value
#p = 0.219
plot(density(model.mat[,1]),lwd=2,
     xlab='Crown Parameter',main=NA,cex.lab=1.5)
polygon(density(model.mat[,1]),bg='black',col='gray')
abline(v=median(model.mat[,1]),lwd=3,lty=2)
abline(v=real.vec[1],lwd=3,col='red')

#Species parameter#
length(model.mat[,2][model.mat[,2]<(median(model.mat[,2])-abs(median(model.mat[,2])-real.vec[2]))|
                       model.mat[,2]>(median(model.mat[,2])+abs(median(model.mat[,2])-real.vec[2]))])/1000
#p = 0.004
plot(density(model.mat[,2]),lwd=2,
     xlab='Species Parameter',main='Uniform Diversification',cex.main=2,cex.lab=1.5)
polygon(density(model.mat[,2]),bg='black',col='gray')
abline(v=median(model.mat[,2]),lwd=3,lty=2)
abline(v=real.vec[2],lwd=3,col='red')

#Interaction parameter#
length(model.mat[,3][model.mat[,3]<(median(model.mat[,3])-abs(median(model.mat[,3])-real.vec[3]))|
                       model.mat[,3]>(median(model.mat[,3])+abs(median(model.mat[,3])-real.vec[3]))])/1000
#p = 0.674
plot(density(model.mat[,3]),lwd=2,
     xlab='Interaction Parameter',main=NA,cex.lab=1.5)
polygon(density(model.mat[,3]),bg='black',col='gray')
abline(v=median(model.mat[,3]),lwd=3,lty=2)
abline(v=real.vec[3],lwd=3,col='red')

###Make 3D plots with regression surfaces###
bdsim<-do.call(rbind,sim.posterior) #Do call is magical
bdsim<-as.data.frame(bdsim)
bdsim$Tips<-log(bdsim$Tips)
bdsim$CHull<-bdsim$CHull^0.25
plot3d(bdsim[,1],bdsim[,2],bdsim[,4])

bdmod<-lm(CHull~Age*Tips,data=bdsim)
summary(bdmod)

fit.age<-seq(from=0,to=30,length.out=100)
fit.tip<-seq(from=0,to=7,length.out=100)
grd<-expand.grid(Age=fit.age,Tips=fit.tip)
grd$pred<-predict(bdmod,newdata=grd)
plot3d(bdsim[,1],bdsim[,2],bdsim[,4],xlab='Crown Age (Ma)',
       ylab='ln(Species Richness)',zlab='Transformed Hypervolume')
points3d(nwdata$Crown,log(nwdata$Species),nwdata$Chull^0.25,size=10,col='red') #Add the real clade disparities too
persp3d(x=unique(grd[[1]]),y=unique(grd[[2]]), 
        z=matrix(grd[[3]],100,100),add=TRUE,col='red',alpha=0.5)

###Compare real vs simulated convex hull disparities###
real.chull<-nwdata$Chull[c(1:19,21:25,27:28)]^0.25

#Write function to generate a summary statistic for differences b/w real & sims CHulls#
get.stat<-function(vector,matrix){
  stats<-vector(length=length(vector),mode='numeric')
  for(i in 1:length(vector)){
    stats[i]<-length(matrix[,i][matrix[,i]<(median(matrix[,i])-abs(median(matrix[,i])-vector[i]))|
                                  matrix[,i]>(median(matrix[,i])+abs(median(matrix[,i])-vector[i]))])/nrow(matrix)
  }
  stats
}

bd.chull<-matrix(nrow=1000,ncol=26)
colnames(bd.chull)<-nwdata$Lineage[c(1:19,21:25,27:28)]
for(i in 1:1000){
  bd.chull[i,]<-sim.posterior[[i]][,4]^0.25
}

#Assemble results matrix#
bd.results<-data.frame(real.chull,apply(bd.chull,2,median),get.stat(real.chull,bd.chull))
colnames(bd.results)<-c('RealCHull','SimMedian','Stat')
View(bd.results)
write.csv(bd.results,'UniformBD_CHulls.csv')
################################################################


###K Decelerating Speciation Model###
load('K-BD_Sims.RData')

for(n in 1:1000){
  model.mat[n,1]<-sim.models[[n]]$coefficients[[2]]
  model.mat[n,2]<-sim.models[[n]]$coefficients[[3]]
  model.mat[n,3]<-sim.models[[n]]$coefficients[[4]]
  model.mat[n,4]<-summary(sim.models[[n]])$adj.r.squared
}

#Crown parameter#
length(model.mat[,1][model.mat[,1]<(median(model.mat[,1])-abs(median(model.mat[,1])-real.vec[1]))|
                       model.mat[,1]>(median(model.mat[,1])+abs(median(model.mat[,1])-real.vec[1]))])/1000
#Calculate number of obvs on either end of distribution more extreme than difference of median and observed parameter value
#p = 0.206
plot(density(model.mat[,1]),lwd=2,
     xlab='Crown Parameter',main=NA,cex.lab=1.5)
polygon(density(model.mat[,1]),bg='black',col='gray')
abline(v=median(model.mat[,1]),lwd=3,lty=2)
abline(v=real.vec[1],lwd=3,col='red')

#Species parameter#
length(model.mat[,2][model.mat[,2]<(median(model.mat[,2])-abs(median(model.mat[,2])-real.vec[2]))|
                       model.mat[,2]>(median(model.mat[,2])+abs(median(model.mat[,2])-real.vec[2]))])/1000
#p = 0.002
plot(density(model.mat[,2]),lwd=2,
     xlab='Species Parameter',main='Decelerating Speciation',cex.main=2,cex.lab=1.5)
polygon(density(model.mat[,2]),bg='black',col='gray')
abline(v=median(model.mat[,2]),lwd=3,lty=2)
abline(v=real.vec[2],lwd=3,col='red')

#Interaction parameter#
length(model.mat[,3][model.mat[,3]<(median(model.mat[,3])-abs(median(model.mat[,3])-real.vec[3]))|
                       model.mat[,3]>(median(model.mat[,3])+abs(median(model.mat[,3])-real.vec[3]))])/1000
#p = 0.211
plot(density(model.mat[,3]),lwd=2,
     xlab='Interaction Parameter',main=NA,cex.lab=1.5)
polygon(density(model.mat[,3]),bg='black',col='gray')
abline(v=median(model.mat[,3]),lwd=3,lty=2)
abline(v=real.vec[3],lwd=3,col='red')


#Make 3D plots#
ksim<-do.call(rbind,sim.posterior)
ksim<-as.data.frame(ksim)
ksim$Tips<-log(ksim$Tips)
ksim$CHull<-ksim$CHull^0.25
plot3d(ksim[,1],ksim[,2],ksim[,4])

kmod<-lm(CHull~Age*Tips,data=ksim) #-1 removes the intercept term
summary(kmod)

fit.age<-seq(from=0,to=30,length.out=100)
fit.tip<-seq(from=0,to=8,length.out=100)
grd<-expand.grid(Age=fit.age,Tips=fit.tip)
grd$pred<-predict(kmod,newdata=grd)
plot3d(ksim[,1],ksim[,2],ksim[,4],xlab='Crown Age (Ma)',
       ylab='ln(Species Richness)',zlab='Transformed Hypervolume')
points3d(nwdata$Crown,log(nwdata$Species),nwdata$Chull^0.25,size=10,col='red')
persp3d(x=unique(grd[[1]]),y=unique(grd[[2]]), 
        z=matrix(grd[[3]],100,100),add=TRUE,col='red',alpha=0.5)


###Compare real vs simulated convex hull disparities###
k.chull<-matrix(nrow=1000,ncol=26)
colnames(k.chull)<-nwdata$Lineage[c(1:19,21:25,27:28)]
for(i in 1:1000){
  k.chull[i,]<-sim.posterior[[i]][,4]^0.25
}

#Assemble results matrix#
k.results<-data.frame(real.chull,apply(k.chull,2,median),get.stat(real.chull,k.chull))
colnames(k.results)<-c('RealCHull','SimMedian','Stat')
View(k.results)
write.csv(k.results,'KBD_CHulls.csv')
################################################################


###Outlier Birth-death###
load('OutlierBD_Sims.RData')

for(n in 1:1000){
  model.mat[n,1]<-sim.models[[n]]$coefficients[[2]]
  model.mat[n,2]<-sim.models[[n]]$coefficients[[3]]
  model.mat[n,3]<-sim.models[[n]]$coefficients[[4]]
  model.mat[n,4]<-summary(sim.models[[n]])$adj.r.squared
}

#Crown parameter#
length(model.mat[,1][model.mat[,1]<(median(model.mat[,1])-abs(median(model.mat[,1])-real.vec[1]))|
                       model.mat[,1]>(median(model.mat[,1])+abs(median(model.mat[,1])-real.vec[1]))])/1000
#Calculate number of obvs on either end of distribution more extreme than difference of median and observed parameter value
#p = 0.115
plot(density(model.mat[,1]),lwd=2,
     xlab='Crown Parameter',main=NA,cex.lab=1.5)
polygon(density(model.mat[,1]),bg='black',col='gray')
abline(v=median(model.mat[,1]),lwd=3,lty=2)
abline(v=real.vec[1],lwd=3,col='red')

#Species parameter#
length(model.mat[,2][model.mat[,2]<(median(model.mat[,2])-abs(median(model.mat[,2])-real.vec[2]))|
                       model.mat[,2]>(median(model.mat[,2])+abs(median(model.mat[,2])-real.vec[2]))])/1000
#p = 0.000
plot(density(model.mat[,2]),lwd=2,xlim=c(0.06,0.19),
     xlab='Species Parameter',main='Outlier Diversification',cex.main=2,cex.lab=1.5)
polygon(density(model.mat[,2]),bg='black',col='gray')
abline(v=median(model.mat[,2]),lwd=3,lty=2)
abline(v=real.vec[2],lwd=3,col='red')

#Interaction parameter#
length(model.mat[,3][model.mat[,3]<(median(model.mat[,3])-abs(median(model.mat[,3])-real.vec[3]))|
                       model.mat[,3]>(median(model.mat[,3])+abs(median(model.mat[,3])-real.vec[3]))])/1000
#p = 0
plot(density(model.mat[,3]),lwd=2,xlim=c(0.0018,0.008),
     xlab='Interaction Parameter',main=NA,cex.lab=1.5)
polygon(density(model.mat[,3]),bg='black',col='gray')
abline(v=median(model.mat[,3]),lwd=3,lty=2)
abline(v=real.vec[3],lwd=3,col='red')

#Make 3D plots of regression surface#
outsim<-do.call(rbind,sim.posterior)
outsim<-as.data.frame(outsim)
outsim$Tips<-log(outsim$Tips)
outsim$CHull<-outsim$CHull^0.25
plot3d(outsim[,1],outsim[,2],outsim[,4])

outmod<-lm(CHull~Age*Tips,data=outsim) #-1 removes the intercept term
summary(outmod)

fit.age<-seq(from=0,to=30,length.out=100)
fit.tip<-seq(from=0,to=8,length.out=100)
grd<-expand.grid(Age=fit.age,Tips=fit.tip)
grd$pred<-predict(outmod,newdata=grd)
plot3d(outsim[,1],outsim[,2],outsim[,4],xlab='Crown Age (Ma)',
       ylab='ln(Species Richness)',zlab='Transformed Hypervolume')
points3d(nwdata$Crown,log(nwdata$Species),nwdata$Chull^0.25,size=10,col='red')
persp3d(x=unique(grd[[1]]),y=unique(grd[[2]]), 
        z=matrix(grd[[3]],100,100),add=TRUE,col='red',alpha=0.5)

###Compare real vs simulated convex hull disparities###
out.chull<-matrix(nrow=1000,ncol=26)
colnames(out.chull)<-nwdata$Lineage[c(1:19,21:25,27:28)]
for(i in 1:1000){
  out.chull[i,]<-sim.posterior[[i]][,4]^0.25
}

#Assemble results matrix#
out.results<-data.frame(real.chull,apply(out.chull,2,median),get.stat(real.chull,out.chull))
colnames(out.results)<-c('RealCHull','SimMedian','Stat')
View(out.results)
write.csv(out.results,'OutlierBD_CHulls.csv')
################################################################
