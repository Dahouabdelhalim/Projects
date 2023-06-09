phy.imbalance<-function(phy,b=1,d=0,nslice,nsim,show.plot=NULL,DVI.range=NULL,sim.trees=NULL,opacity=0.2){
  # phy is the tree
  # b is the speciation rate, can be a function of time, see TESS
  # d is the extinction rate
  # nslice is the number of time points to evaluate variance in des lineages
  # nsim is the number of simulated trees used to evaluated DVI
  # DVI range is the time space over which to evaluate DVI, can be used to exclude for incomplete
  # sim.trees can be a list/multiphylo object of user generated trees, b and d are ignored but nsim is still necessary
  # show.plot=NULL produces no plot, "lines" shows the result of each sim, "poly" produces a shaded 95% interval/range
  
  #function outputs DVI and associated p-value, and values for recreating plots
  require(geiger)
  require(phytools)
  require(svMisc)
  area.between.curves<-geiger:::.area.between.curves
  dtt.polygon<-geiger:::.dtt.polygon
  if(is.null(DVI.range)){DVI.range<-c(0,max(branching.times(phy)))}
  
  
  PI<-vector(length=nslice)
  oldest<-max(branching.times(phy))
  youngest<-min(branching.times(phy))
  slices<-seq(0,(oldest-youngest+0.0001),length.out=(nslice+1))
  slices<-slices[-1]
  
  for(i in 1:nslice){
    subtrees<-treeSlice(phy,slices[i],trivial=TRUE)
    PI[i]<-var(log10(sapply(subtrees,Ntip)))
  }
  
  
  PI.sim<-matrix(nrow=nsim,ncol=nslice)
  
  if(!is.null(sim.trees)){
    for(i in 1:length(sim.trees)){
      progress(i,max.value=length(sim.trees))
      for(m in 1:nslice){
        sub.sims<-treeSlice(sim.trees[[i]],slices[m],trivial=TRUE)
        PI.sim[i,m]<-var(log10(sapply(sub.sims,Ntip)))
      }
    }
    
  } else {
    
    
    
    sim.trees<-tess.sim.taxa.age(n=nsim,nTaxa=length(phy$tip.label),age=max(branching.times(phy)),lambda=b,mu=d)
    for(i in 1:nsim){
      progress(i,max.value=nsim)
      for(m in 1:nslice){
        sub.sims<-treeSlice(sim.trees[[i]],slices[m],trivial=TRUE)
        PI.sim[i,m]<-var(log10(sapply(sub.sims,Ntip)))
      }
    }
    
  }
  median.PI<-apply(PI.sim,2,median)
  
  if(!(is.null(show.plot))){
  if(show.plot=="lines"){
    plot(PI.sim[1,]~slices,type="l",col=rgb(0.5,0.5,0.5,opacity),ylim=c(0,max(PI.sim)),ylab="Phylogenetic Imbalance",xlab="Time since root")
    for(i in 2:nsim){
      lines(PI.sim[i,]~slices,type="l",col=rgb(0.5,0.5,0.5,opacity))
    }
    lines(PI~slices,type="l",lwd=2)
    lines(median.PI~slices,type="l",lty=2,lwd=2)
  } 
  
  if(show.plot=="poly"){
    plot.new()
    poly<-dtt.polygon(t(PI.sim),slices,alpha=0.05)
    plot(PI~slices,type="l",lwd=2,ylim=c(0,max(poly[,2])))
    polygon(poly[,1], poly[,2], col = "grey", border = NA)
    lines(PI~slices,type="l",lwd=2)
    lines(median.PI~slices,type="l",lty=2,lwd=2)
  }
  }   
  
  
  DVI<-area.between.curves(slices,median.PI,PI,sort(DVI.range))
  
  
  DVI.sims<-vector(length=nsim)
  for(i in 1:nsim){
    DVI.sims[i]<-area.between.curves(slices,median.PI,PI.sim[i,],sort(DVI.range))
  }
  
  p.value<-sum(ifelse(DVI.sims>=DVI,1,0))/nsim
  
  results<-list(PI=PI,time.slices=slices,PI.sim=PI.sim,DVI=DVI,p.value=p.value)
  return(results)
  
}
