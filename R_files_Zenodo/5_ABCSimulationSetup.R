##########################################################################################
###This script will perform the ABC simulations to obtain optimized sigma squared rates###
##########################################################################################

##########################################################
###Set up Uniform Birth-death ABC simulation and run it###
sim.fun.bd<-function(x){
require(phytools)
require(geometry)
require(geiger)
require(TreeSim)

  clade.data<-read.csv('WH_CladeData.csv')
  real.times<-clade.data$Crown
  real.tips<-log(clade.data$Species)
  real.chulls<-clade.data$Chull^0.25
  

  get.bd.age.tree<-function(pars){
    success<-F
    while(!success){
      try<-sim.bd.age(pars["t"],1,pars["b"],pars["d"],frac=1,mrca=T,complete=FALSE)
      success<-Ntip(try[[1]])>=5
   }
    return(try[[1]])
  }

  set.seed(x[1])
  #simulate trees
  #set times
  #Optimal parameters: r=0.1169096, e=0.7053727, lnL=-150.8963
  #Translated: lambda=0.3968050, mu=0.2798954
  trees<-list(NULL)
  for (i in 1:length(real.times)){
    trees[[i]]<-get.bd.age.tree(c(b=0.396805,d=0.2798954,t=real.times[i]))
    print(paste("Tree",i,"done"))
  }
  ntips<-log(sapply(trees,Ntip))

  #calculate species.diff statistic
  #species.diff<-sum(abs(log(ntips)-real.tips))
  #print(paste("Net difference in log species diversity:",species.diff))
  
  #simulate trait data
  chulls<-NULL
  var.props<-c(76.52,10.71,5.62,2.74)/sum(c(76.52,10.71,5.62,2.74))
  for (i in 1:length(real.times)){
    tree.data<-matrix(nrow=Ntip(trees[[i]]),ncol=4) #Create blank matrix to store tip data for 4 "PC axes"
    for (j in 1:4){ #Simulate the trait data given a set of sigma squared rates
      tree.data[,j]<-fastBM(trees[[i]],a=0,sig2=x[2]*var.props[j],internal=F)
    }
    print(paste("Data for tree",i,"done"))
    chulls[[i]]<-convhulln(tree.data,option='FA')$vol
  }
  
  #data to compare:
  real.data<-cbind(real.times,real.tips,real.chulls)
  sim<-cbind(real.times,ntips,chulls^0.25)
  dists<-NULL
  #below is the Euclidean distance between the real and simulated data
  total.dist<-sum(sqrt(diag((real.data-sim)%*%t(real.data-sim))))
  
  print(paste("Overall Euclidean distance:",total.dist))
  
  #return summary statistic (total Euclidean distance)
  total.dist
}

###Run it###
priors<-list(c("unif",0.01,0.25))

abc.rej.bdfix<-ABC_rejection(model=sim.fun.bd,prior=priors,nb_simul=10000,tol=0.2,summary_stat_target=0,verbose=T,use_seed=T,n_cluster=4)

save.image('ABC_SimResults_UniformBD.RData')

plot(abc.rej.bdfix$param,abc.rej.bdfix$stats,pch=19,xlab='SigmaSq Parameter',
     ylab='ABC Stat')

##########################################################


##############################################################
###Set up Decelerating Speciation ABC simulation and run it###
sim.fun.K<-function(x){
  require(phytools)
  require(geometry)
  require(geiger)
  require(caper)
  require(TreeSim)
  
  clade.data<-read.csv('WH_CladeData.csv')
  real.times<-clade.data$Crown
  real.tips<-log(clade.data$Species)
  real.chulls<-clade.data$Chull^0.25
  
  
  get.bd.age.tree<-function(pars){
    success<-F
    while(!success){
      try<-sim.bd.age(pars["t"],1,pars["b"],0,frac=1,mrca=T,complete=FALSE,K=pars["K"])
      success<-Ntip(try[[1]])>=5
    }
    return(try[[1]])
  }
  
  set.seed(x[1])
  #simulate trees
  trees<-list(NULL)
  #ML values for declining speciation: l0=0.2655675, k=0.0410325, lnL=-151.1335
  #starting speciation: 0.2655675, ending speciation: 0.08347893
  #max diversity is 826, so set simulations to have lambda of 0.0835 at 826 spp, which corresponds to K of 1205
  for (i in 1:length(real.times)){
    trees[[i]]<-get.bd.age.tree(c(b=0.2655675,t=real.times[i],K=1205))
    print(paste("Tree",i,"done"))
  }
  ntips<-log(sapply(trees,Ntip))
  
  #calculate species.diff statistic
  species.diff<-sum(abs(log(ntips)-real.tips))
  print(paste("Net difference in log species diversity:",species.diff))
  
  #simulate trait data
  chulls<-NULL
  var.props<-c(76.52,10.71,5.62,2.74)/sum(c(76.52,10.71,5.62,2.74))
  for (i in 1:length(real.times)){
    tree.data<-matrix(nrow=Ntip(trees[[i]]),ncol=4) #Create blank matrix to store tip data for 4 "PC axes"
    for (j in 1:4){ #Simulate the trait data given a set of sigma squared rates
      tree.data[,j]<-fastBM(trees[[i]],a=0,sig2=x[2]*var.props[j],internal=F)
    }
    print(paste("Data for tree",i,"done"))
    chulls[[i]]<-convhulln(tree.data,option='FA')$vol
  }
  
  #data to compare:
  real.data<-cbind(real.times,real.tips,real.chulls)
  sim<-cbind(real.times,ntips,chulls^0.25)
  dists<-NULL
  
  #below is the Euclidean distance between the real and simulated data
  total.dist<-sum(sqrt(diag((real.data-sim)%*%t(real.data-sim))))
  
  print(paste("Overall Euclidean distance:",total.dist))
  
  #return summary statistic (total Euclidean distance)
  total.dist
  
}

###Run it###
priors<-list(c("unif",0.01,0.25))

abc.rej.K<-ABC_rejection(model=sim.fun.K,prior=priors,nb_simul=10000,tol=0.2,summary_stat_target=0,verbose=T,use_seed=T,n_cluster=4)

save.image('ABC_SimResults_K.RData')

plot(abc.rej.K$param,abc.rej.K$stats,pch=19,xlab='SigmaSq Parameter',
     ylab='ABC Stat')

##############################################################


##########################################################
###Set up Outlier Birth-death ABC simulation and run it###
sim.fun.bd<-function(x){
  require(phytools)
  require(geometry)
  require(geiger)
  require(TreeSim)
  
  clade.data<-na.omit(read.csv('WH_CladeData.csv'))
  real.times<-clade.data$Crown
  real.tips<-log(clade.data$Species)
  real.chulls<-clade.data$Chull^0.25
  lambdas<-c(0.6894,rep(0.3968050,6),0.4997,rep(0.3968050,8),0.2165,rep(0.3968050,9)) #Assign unique speciation rates to diversification outliers

  get.bd.age.tree<-function(pars){
    success<-F
    while(!success){
      try<-sim.bd.age(pars["t"],1,pars["b"],pars["d"],frac=1,mrca=T,complete=FALSE)
      success<-Ntip(try[[1]])>=5
    }
    return(try[[1]])
  }
  
  set.seed(x[1])
  #simulate trees
  #set times
  #Optimal parameters: r=0.1169096, e=0.7053727, lnL=-150.8963
  #Translated: lambda=0.3968050, mu=0.2798954
  trees<-list(NULL)
  for (i in 1:length(real.times)){
    trees[[i]]<-get.bd.age.tree(c(b=lambdas[i],d=0.2798954,t=real.times[i]))
    print(paste("Tree",i,"done"))
  }
  ntips<-log(sapply(trees,Ntip))
  
  #calculate species.diff statistic
  #species.diff<-sum(abs(log(ntips)-real.tips))
  #print(paste("Net difference in log species diversity:",species.diff))
  
  #simulate trait data
  chulls<-NULL
  var.props<-c(76.52,10.71,5.62,2.74)/sum(c(76.52,10.71,5.62,2.74))
  for (i in 1:length(real.times)){
    tree.data<-matrix(nrow=Ntip(trees[[i]]),ncol=4) #Create blank matrix to store tip data for 4 "PC axes"
    for (j in 1:4){ #Simulate the trait data given a set of sigma squared rates
      tree.data[,j]<-fastBM(trees[[i]],a=0,sig2=x[2]*var.props[j],internal=F)
    }
    print(paste("Data for tree",i,"done"))
    chulls[[i]]<-convhulln(tree.data,option='FA')$vol
  }
  
  #data to compare:
  real.data<-cbind(real.times,real.tips,real.chulls)
  sim<-cbind(real.times,ntips,chulls^0.25)
  dists<-NULL
  #below is the Euclidean distance between the real and simulated data
  total.dist<-sum(sqrt(diag((real.data-sim)%*%t(real.data-sim))))
  
  print(paste("Overall Euclidean distance:",total.dist))
  
  #return summary statistic (total Euclidean distance)
  total.dist
}

###Run it###
priors<-list(c("unif",0.01,0.25))

abc.rej.bdfix<-ABC_rejection(model=sim.fun.bd,prior=priors,nb_simul=10000,tol=0.2,summary_stat_target=0,verbose=T,use_seed=T,n_cluster=4)

save.image('ABC_SimResults_OutlierBD.RData')

plot(abc.rej.bdfix$param,abc.rej.bdfix$stats,pch=19,xlab='SigmaSq Parameter',
     ylab='ABC Stat')

##########################################################
