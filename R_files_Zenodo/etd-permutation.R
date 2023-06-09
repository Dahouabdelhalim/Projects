# batch job submission
rm(list=ls())
setwd("/scratch/osalehzi/phy/Final/")

XXX<-commandArgs(trailingOnly = TRUE)

require(secsse)
library(DDD)
require(phytools)
require(data.table)

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
st<-unlist(lapply(phyl$tip.label,function(x) data[species==x]$poly))

# container list of out models,
out<-list()

# tree run,
if(TRUE){
  for(i in 1:100){
    tree<-phyl
    st<-sample(st)
    x<-setNames(as.integer(st),phyl$tip.label)
    
    # Load in SecSSE parameters,
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
    # match examined and concealed transitions
    idparslist[[3]][7,1]<-0
    idparslist[[3]][8,2]<-0
    idparslist[[3]][9,3]<-0
    idparslist[[3]][9,7]<-0
    idparslist[[3]][1,7]<-0
    idparslist[[3]][2,8]<-0
    idparslist[[3]][3,9]<-0
    idparslist[[2]][]<-10   # all models will have a fixed trait independent extinction rate
    diag(idparslist[[3]]) <- NA
    # final params,
    c<-na.omit(idparslist$Q[idparslist$Q>0])
    attributes(c)$na.action<-NULL
    sampling.f<-rep(584/4644,3)
    
    # Setting up ETD parameters,
    idparslist[[1]][]<-c(rep(c(1,2,3),3)) # ETD has fixed rates for examined states
    idparsopt<-c(c(1,2,3),c)
    idparsfix<-c(0,10) # fix 0 and extinction (10)
    startingpoint <- bd_ML(brts = ape::branching.times(tree))
    intGuessLamba <- startingpoint$lambda0
    intGuessMu <- startingpoint$mu0
    initparsopt <- c(rep(intGuessLamba,3),rep(0.25,length(c)))
    parsfix<-c(0,rep((intGuessLamba/5),1)) 
    
    t1<-paste("Starting permuation : ",i," at: ",Sys.time())
    t1
    out[[i]]<-secsse_ml(tree,x, num_concealed_states=3, 
                        idparslist, idparsopt, initparsopt, idparsfix, parsfix, 
                        cond="proper_cond",root_state_weight = "proper_weights", 
                        sampling_fraction=sampling.f, tol = c(1e-04, 1e-05, 1e-07), 
                        maxiter = 1000 * round((1.25)^length(idparsopt)), use_fortran=TRUE,
                        methode="ode45", optimmethod = "simplex", num_cycles = 1, run_parallel=TRUE)
    
    t2<-paste("Finsihed permutation: ",i," at: ",Sys.time())
    print(t2)
    file<-paste("musse/SecSSE/random100etd",XXX,".RDS",sep='')
    saveRDS(out, file=file)
  }
}

require(plyr)
require(ggplot2)

etd<-readRDS("musse/SecSSE/ETD.RDS") # empirical observation
b<-etd$ML

# sloppy way of getting it
files<-grep("100.*RDS",list.files("musse/SecSSE/Permutation files/"),value=T)
rand<-lapply(files,function(x) readRDS(paste("musse/SecSSE/Permutation files/",x,sep='',collapse='')))

# saveRDS(rand,file="musse/SecSSE/Permutation files/randomSecSSEs.RDS")
rand<-readRDS("musse/SecSSE/Permutation files/randomSecSSEs.RDS")

a<-ldply(do.call(c,rand),function(x) x$ML)
a

ggplot(a,aes(V1))+geom_density()+geom_vline(aes(xintercept=b),color="red",linetype="dashed")

