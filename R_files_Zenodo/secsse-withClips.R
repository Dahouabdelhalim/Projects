# SecSSE script
rm(list=ls())
setwd("/scratch/osalehzi/phy/Final/")

require(secsse)
require(phytools)
require(data.table)
require(geiger)

###########################################################################
# Setup data and phylogenies,
###########################################################################

# load phylogenies and data,
phyl<-read.tree(file="input/Aphidoidea-MCC.tre")
data<-fread("input/aphiddata.csv",header=T,drop="V1")
phy.names<-fread("input/phy-names.csv",header = T,drop = "V1") # this matches the order of the mcc tre 
data$species<-gsub(" ","_",data$species)

# reformat state data: change data to mono-polymorphisms
data$poly<-"NA"
data[wings=="dimorphic"]$poly<-1
data[wings=="alate"]$poly<-0
data[wings=="apterous"]$poly<-2

# get states of our mcc tree
phyl$tip.label<-phy.names$tip_labels

# subset out macrosiphini from other aphididae (node 589)
phyl<-drop.tip(phyl,getDescendants(phyl,589)[getDescendants(phyl,589)<=584])
# subset in only macrosiphini clade,
#phyl<-drop.tip(phyl,-getDescendants(hyl,589)[getDescendants(phyl,589)<=584])

st<-unlist(lapply(phyl$tip.label,function(x) data[species==x]$poly))
x<-setNames(as.integer(st),phyl$tip.label)

# non-macrosiphini have 2 dimorphic, 107 apterous, and 72 winged species, 53 unknowns (tknown=181, total=234)


# Setup parameterization for SecSSE

idparslist <- id_paramPos(x, num_concealed_states = 3)

# make dual transitions 0
idparslist[[3]][1,c(5,6,8,9)] <- 0
idparslist[[3]][2,c(4,6,7,9)] <- 0
idparslist[[3]][3,c(4,5,7,8)] <- 0
idparslist[[3]][4,c(2,3,8,9)] <- 0
idparslist[[3]][5,c(1,3,7,9)] <- 0
idparslist[[3]][6,c(1,2,7,8)] <- 0
idparslist[[3]][7,c(2,3,5,6)] <- 0
idparslist[[3]][8,c(1,3,4,6)] <- 0
idparslist[[3]][9,c(1,2,4,5)] <- 0

# make monomorphic transitions 0 (q02=q20=0)
idparslist[[3]][1,c(3,6,9)]<-0
idparslist[[3]][3,c(1,4,7)]<-0
idparslist[[3]][4,c(3,6,9)]<-0
idparslist[[3]][6,c(1,4,7)]<-0
idparslist[[3]][7,c(3,6,9)]<-0
idparslist[[3]][0,c(1,4,7)]<-0

# and per secsse paper, 
# match examined and concealed transitions
idparslist[[3]][7,1]<-0
idparslist[[3]][8,2]<-0
idparslist[[3]][9,3]<-0
idparslist[[3]][9,7]<-0
idparslist[[3]][1,7]<-0
idparslist[[3]][2,8]<-0
idparslist[[3]][3,9]<-0

idparslist[[2]][]<-10 # all models will have a fixed trait independent extinction rate
diag(idparslist[[3]]) <- NA

# setup parameters, 
c<-na.omit(idparslist$Q[idparslist$Q>0])
attributes(c)$na.action<-NULL
#initparsopt<-c(rep(1.2,9),rep(0.1,9),rep(0.25,length(c)))
# idparsopt<-c(1:18,c) # optimizing speciation and extinction rates, 
#sampling.f<-rep(584/4644,3)
sampling.f<-rep(length(phyl$tip.label)/4644,3)
library(DDD)

if(FALSE){
  # CR
  idparslist[[1]][]<-1 # CR model has a fixed speciation rate independent of either concealed or examined states
  
  idparsopt<-c(c) # optimizing transition rates
  idparsfix<-c(0,1,10)
  startingpoint <- bd_ML(brts = ape::branching.times(phyl))
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  initparsopt <- c(rep(0.25,length(c)))
  parsfix<-c(0,rep(intGuessLamba,1),rep((intGuessLamba/5),1))
  
  out<-secsse_ml(phyl,x, num_concealed_states=3, idparslist, idparsopt, initparsopt, idparsfix, parsfix, cond="proper_cond",root_state_weight = "proper_weights", sampling_fraction=sampling.f, tol = c(1e-04, 1e-05, 1e-07), maxiter = 1000 * round((1.25)^length(idparsopt)), use_fortran=TRUE,methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=FALSE)
  Sys.time()
  out
  saveRDS(out, file="musse/SecSSE/secsse-clippingMacrosiphini/CRclipped.RDS")
}

if(FALSE){
  # CTD
  idparslist[[1]][]<-c(rep(1,3),rep(2,3),rep(3,3)) # CTD has fixed rates for concealed states
  
  idparsopt<-c(c(1,2,3),c)
  idparsfix<-c(0,10)
  startingpoint <- bd_ML(brts = ape::branching.times(phyl))
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  initparsopt <- c(rep(intGuessLamba,3),rep(0.25,length(c)))
  parsfix<-c(0,rep((intGuessLamba/5),1)) 
  
  out<-secsse_ml(phyl,x, num_concealed_states=3, idparslist, idparsopt, initparsopt, idparsfix, parsfix, cond="proper_cond",root_state_weight = "proper_weights", sampling_fraction=sampling.f, tol = c(1e-04, 1e-05, 1e-07), maxiter = 1000 * round((1.25)^length(idparsopt)), use_fortran=TRUE,methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=FALSE)
  
  saveRDS(out, file="musse/SecSSE/secsse-clippingMacrosiphini/CTDclipped.RDS")
}

if(TRUE){
  # ETD
  idparslist[[1]][]<-c(rep(c(1,2,3),3)) # ETD has fixed rates for examined states
  
  idparsopt<-c(c(1,2,3),c)
  idparsfix<-c(0,10)
  startingpoint <- bd_ML(brts = ape::branching.times(phyl))
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  initparsopt <- c(rep(intGuessLamba,3),rep(0.25,length(c)))
  parsfix<-c(0,rep((intGuessLamba/5),1)) 
  
  out<-secsse_ml(phyl,x, num_concealed_states=3, idparslist, idparsopt, initparsopt, idparsfix, parsfix, cond="proper_cond",root_state_weight = "proper_weights", sampling_fraction=sampling.f, tol = c(1e-04, 1e-05, 1e-07), maxiter = 1000 * round((1.25)^length(idparsopt)), use_fortran=TRUE,methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=FALSE)
  
  saveRDS(out, file="musse/SecSSE/secsse-clippingMacrosiphini/ETDclipped2.RDS")
  save(out,file="musse/SecSSE/secsse-clippingMacrosiphini/etdClipped2.RData")
}

if(FALSE){
  # June 7, 2022: intermediate model where winged=poly != wingless to discern signal.
  idparslist[[1]][]<-c(rep(c(1,1,2),3))
  
  idparsopt<-c(c(1,2),c)
  idparsfix<-c(0,10)
  startingpoint <- bd_ML(brts = ape::branching.times(phyl))
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  initparsopt <- c(rep(intGuessLamba,2),rep(0.25,length(c)))
  parsfix<-c(0,rep((intGuessLamba/5),1)) 
  
  out<-secsse_ml(phyl,x, num_concealed_states=3, idparslist, idparsopt, initparsopt, idparsfix, parsfix, cond="proper_cond",root_state_weight = "proper_weights", sampling_fraction=sampling.f, tol = c(1e-04, 1e-05, 1e-07), maxiter = 1000 * round((1.25)^length(idparsopt)), use_fortran=TRUE,methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=FALSE)
  
  saveRDS(out, file="musse/SecSSE/inter-ETD2.RDS")
  save(out,file="musse/SecSSE/inter-etd2.RData")
  print(Sys.time())
}

if(FALSE){
  idparslist[[1]][]<-c(rep(c(1,2,3),3)) # ETD has fixed rates for examined states
  
  idparsopt<-c(c(1,2,3),c)
  idparsfix<-c(0,10)
  startingpoint <- bd_ML(brts = ape::branching.times(phyl))
  intGuessLamba <- startingpoint$lambda0
  intGuessMu <- startingpoint$mu0
  initparsopt <- c(rep(intGuessLamba,3),rep(0.25,length(c)))
  parsfix<-c(0,rep((intGuessLamba/5),1)) 
  x<-setNames(sample(as.integer(st)),phyl$tip.label) # randomize tips
  
  out<-secsse_ml(phyl,x, num_concealed_states=3, idparslist, idparsopt, initparsopt, idparsfix, parsfix, cond="proper_cond",root_state_weight = "proper_weights", sampling_fraction=sampling.f, tol = c(1e-04, 1e-05, 1e-07), maxiter = 1000 * round((1.25)^length(idparsopt)), use_fortran=TRUE,methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=FALSE)
  
  saveRDS(out, file="musse/SecSSE/ETD-random.RDS")
}

if(FALSE){
  cr<-readRDS("musse/SecSSE/secsse-clippingMacrosiphini/CRclipped.RDS")
  ctd<-readRDS("musse/SecSSE/secsse-clippingMacrosiphini/CTDclipped.RDS")
  etd<-readRDS("musse/SecSSE/secsse-clippingMacrosiphini/ETDclipped.RDS")
  etdr<-readRDS("musse/SecSSE/secsse-clippingMacrosiphini/ETDclipped2.RDS")
  
  mods<-list(cr,ctd,etd,etdr) 
  m<-c("cr","ctd","etd","etdr") # random min/max
  logs<-setNames(c(cr$ML,ctd$ML,etd$ML,etdr$ML),m)
  k<-unlist(lapply(mods,function(x) length(unique(unlist(x$MLpars)))-2)) # -2 for NA and 0s
  #install.packages("AICcmodavg")
  require(AICcmodavg)
  aic<-unlist(lapply(1:length(mods),function(x) AICcCustom(logs[[x]],k[[x]],nobs=584)))
  #aic<-unlist(lapply(1:length(mods),function(x) AIC(mods[[x]],k=k[x])))
  
  dt<-data.table(model=m,k=k,logs,aic,aic.w=signif(aic.w(aic),4))
  dt<-dt[order(dt$k,decreasing = T)]
  dt$dAIC<-dt$aic-dt[model=="ctd"]$aic # 'best' fit model is extinction
  dt$er<-dt$aic.w/dt[model=="ctd"]$aic.w
  dt
  
  # is that all? bayesian mcmc sim on best model parameter estimates to get distributions?
  # do we care about the transition rates?
  q<-list(etd$MLpars[[3]][1:3,1:3],etd$MLpars[[3]][4:6,4:6],etd$MLpars[[3]][7:9,7:9])
  q.etd<-apply(simplify2array(q), 1:2, mean) # means across examined state-dependent transition rates
  q.etd
  
  q<-list(ctd$MLpars[[3]][1:3,1:3],ctd$MLpars[[3]][4:6,4:6],ctd$MLpars[[3]][7:9,7:9])
  q.ctd<-apply(simplify2array(q), 1:2, mean) # means across examined state-dependent transition rates
  q.ctd
  
  
  pars<-c(etd$MLpars[[1]][1:3],etd$MLpars[[2]][1],q[which(!is.na(q) & q>0)])
  pars<-setNames(pars,c("lambda1","lambda2","lambda3","mu1","q21","q12","q32","q23")) # 1,2,3 is 0,1,2 resp
  #save(pars,file="musse/SecSSE/pars.RData")
}