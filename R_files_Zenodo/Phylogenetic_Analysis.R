# This script reads in trees and wing data, reconstructs ancestral states, simulates character histories, and plots the data.

setwd("path/to/files")

# Load Packages
require(data.table)
require(phytools)
require(stringr)
require(plyr)

# Import datasets and phylogenies
phyl<-read.tree(file="input/Aphidoidea-MCC.tre")
phy.names<-fread("input/phy-names.csv",header = T,drop = "V1") # this matches the order of the mcc tre; and has species names corrections
data<-fread("input/aphiddata.csv")
trees<-read.tree(file="input/Aphidoidea-100HPP.tre")

##################
# Primping 
# 
# Using the Hardy et al., 2015 phylogeny data available publically, we're just reformating the naming schemes and some species names,
# and pruning trees of species that don't have male wing data.
#
##################

# Changing tip label format
for(i in 1:length(trees)){
  trees[[i]]$tip.label<-gsub(" ","_",word(gsub("_"," ",trees[[i]]$tip.label),2,3))
}

# Implementing a more consistent taxonomic naming scheme,
cor<-fread("input/phyl-name-corrections.csv",header=T)
cor<-cor[,c("name","correct")]
for(i in 1:length(trees)){
  wr<-setdiff(trees[[i]]$tip.label,phy.names$tip_labels) # use phy.names to find typo tip labels
  ri<-gsub(" ","_",cor[match(gsub("_"," ",wr),cor$name)]$correct) # order cor in the same way as the order of wr then use cor to subset tips with the correct name
  trees[[i]]$tip.label[which(trees[[i]]$tip.label %in% wr)]<-ri # subset the wrong tip labels on the tree and replace the right ones
}

species<-data[!is.na(wings)]$species
# remove species not in our tree,
for(i in 1:length(trees)){
  poopy<-setdiff(trees[[i]]$tip.label,gsub(" ","_",species))
  trees[[i]]<-drop.tip(trees[[i]],poopy)
}

# loading background functions
treeName<-function(tips){ # function to get a formatted vector of tips and their states
  st<-unlist(lapply(tips,function(x) data[species==gsub("_"," ",x)]$wings))
  sta<-setNames(as.factor(st),tips)
  return(sta)
}

treeName<-function(tips){ # function to get names of a tree
  st<-unlist(lapply(tips,function(x) data[species==gsub("_"," ",x)]$wings))
  st<-gsub("alate","A",gsub("apterous","B",gsub("dimorphic","A+B",st)))
  sta<-setNames(as.factor(st),tips)
  return(sta)
}

stateGetter<-function(t){ # t is the tree input
  state<-as.factor(unlist(lapply(gsub("_"," ",t$tip.label),function(x) data[wings!="" & species==gsub("_"," ",x)]$wings)))
  return(state)
}

# Rewrite trees (multi-phylo object) with mapped states included in the tree files,
signal<-FALSE
if(signal){
  mtrees<-rtree(20) # seed it with a random tree to be removed later
  for(i in 1:length(trees)){
    print(i)
    state<-stateGetter(trees[[i]])
    names(state)<-trees[[i]]$tip.label
    x<-state
    t<-make.simmap(trees[[i]],x)
    mtrees<-append(mtrees,t)
  }
  mtrees<-mtrees[-1]
  lapply(mtrees,function(x) write.simmap(x,file="input/hpptrees-mapped.tre",format='phylip',append = TRUE))
}
if(!signal){
  trees<-read.simmap('input/hpptrees-mapped.tre',format='phylip')}


###############################
## Ancestral Reconstructions ##
###############################

# Make ancestral reconstructions on the maximum credibility clade tree, 
phyl$tip.label<-phy.names$tip_labels
species<-data[!is.na(wings)]$species
tree<-drop.tip(phyl,phyl$tip.label[-match(species,phy.names$tips)])
x<-treeName(tree$tip.label)
x<-factor(x,levels=c("A","B","A+B")) # reorder to match matrix indices below

# specify matrices for different models,
# standard discrete trait models,
er<-matrix(c(0,1,1,1,0,1,1,1,0),3)
sym<-matrix(c(0,1,2,1,0,3,2,3,0),3)
ard<-matrix(c(0,1,2,3,0,4,5,6,0),3)
# restricting transitions between monomorphic states directly,
er.tr<-matrix(c(0,0,1,0,0,1,1,1,0),3)
sym.tr<-matrix(c(0,0,1,0,0,2,1,2,0),3)
ard.tr<-matrix(c(0,0,1,0,0,2,3,4,0),3)
# assumes the only difference is the loss/gain rates of polymorphisms
tr.tr<-matrix(c(0,0,1,0,0,1,2,2,0),3) 

# fit the various models
fitER<-ace(x,tree,model=er,type="discrete") 
fitSYM<-ace(x,tree,model=sym,type="discrete") 
fitARD<-ace(x,tree,model=ard,type="discrete")
fitER.tr<-ace(x,tree,model=er.tr,type="discrete") 
#fitSYM.tr<-ace(x,tree,model=sym.tr,type="discrete") # "system is computationally singular"
fitARD.tr<-ace(x,tree,model=ard.tr,type="discrete") 
fitTR.tr<-ace(x,tree,model=tr.tr,type="discrete")

# construct a table to assess model fits,
m<-c("ER","SYM","ARD","ER.tr","ARD.tr","TR.tr") # ,"SYM.tr"
logs<-setNames(c(fitER$loglik,fitSYM$loglik,fitARD$loglik,fitER.tr$loglik,fitARD.tr$loglik,fitTR.tr$loglik),m) #fitSYM.tr$loglik,
mods<-list(fitER,fitSYM,fitARD,fitER.tr,fitARD.tr,fitTR.tr) # fitSYM.tr,
k<-unlist(lapply(list(er,sym,ard,er.tr,ard.tr,tr.tr),function(x) length(unique(c(x)))-1))# -1 to account for zeros # sym.tr
aic<-unlist(lapply(1:length(mods),function(x) AIC(mods[[x]],k=k[x])))
dt<-data.table(model=m,logs,aic,aic.w=aic.w(aic),k=k) 
dt[order(dt$aic.w,decreasing = T)]

save.image(file = "phylo/pre-loading.RData")

############################
## Stochastic simulations ##
############################

# Here we simulate character histories using stochastic mapping across100 sampled trees from the posterior distribution.

# We use make.simmap with Q='empirical' with the standard models, and use the transition matrix from the model fits for the others,
paste(Sys.time(),"Starting empirical ER",sep=": ")
mtrees.er<-make.simmap(trees,x,model="ER",Q='empirical',nsim=1000,message=FALSE)
save(mtrees.er,file="phylo/mtrees-er.RData")
paste(Sys.time(),"Starting empirical SYM",sep=": ")
mtrees.sym<-make.simmap(trees,x,model="SYM",Q='empirical',nsim=1000,message=FALSE)
save(mtrees.sym,file="phylo/mtrees-sym.RData")
paste(Sys.time(),"Starting empirical ARD",sep=": ")
mtrees.ard<-make.simmap(trees,x,model="ARD",Q='empirical',nsim=1000,message=FALSE)
save(mtrees.ard,file="phylo/mtrees-ard.RData")

# transient models,
paste(Sys.time(),"Starting empirical ER.tr",sep=": ")
mtrees.ertr<-make.simmap(trees,x,model=er.tr,Q='empirical',nsim=1000,message=FALSE)
save(mtrees.ertr,file="phylo/mtrees-ertr.RData")
#paste(Sys.time(),"Starting empirical SYM.tr",sep=": ") # this model runs into a computational error, and considered undoable.
#mtrees.symtr<-make.simmap(trees,x,model=sym.tr,Q='empirical',nsim=1000,message=FALSE)
#save(mtrees.symtr,file="phylo/mtrees-symtr.RData")
paste(Sys.time(),"Starting empirical ARD.tr",sep=": ")
mtrees.ardtr<-make.simmap(trees,x,model=ard.tr,Q='empirical',nsim=1000,message=FALSE)
save(mtrees.ardtr,file="phylo/mtrees-ardtr.RData")
paste(Sys.time(),"Starting empirical TR.tr",sep=": ")
mtrees.trtr<-make.simmap(trees,x,model=tr.tr,Q='empirical',nsim=1000,message=FALSE)
save(mtrees.trtr,file="phylo/mtrees-trtr.RData")

# maximum parsimony reconstructions - jan 2 2022
mtrees.mp.er<-make.simmap(trees,x,model=er,Q='empirical',nsim=1000,message=FALSE,prior=list(alpha=1,beta=1e4))
save(mtrees.mp.er,file="phylo/mtrees-mper.RData")
mtrees.mp.sym<-make.simmap(trees,x,model=sym,Q='empirical',nsim=1000,message=FALSE,prior=list(alpha=1,beta=1e4))
save(mtrees.mp.sym,file="phylo/mtrees-mpsym.RData")
mtrees.mp.ard<-make.simmap(trees,x,model=ard,Q='empirical',nsim=1000,message=FALSE,prior=list(alpha=1,beta=1e4))
save(mtrees.mp.ard,file="phylo/mtrees-mpard.RData")
# for ertr, there are some tree topologies that cause compuation time to become problematic (won't finish after weeks compared to days),
#   and these represent odd cases we exclude. For er.tr, its the 19th tree, trees[[19]] which we excluded. 
#   Because these make.simmap jobs take days to run, in a script not included here, we block off job runs on a subset of trees to make the exclusion...
#   (not shown here, but something like 'trees<-trees[c(1:18,20:100)]' before running the make.simmap )
mtrees.ertr<-make.simmap(trees,x,model=er.tr,Q='empirical',nsim=1000,message=TRUE,prior=list(alpha=1,beta=1e4)) # leaving messages on to track progress,
save(mtrees.ertr,file="phylo/mtrees-mpertr.RData")
# for ardtr, there are also problematic tree topologies (feb 27 2022)
# trees<-trees[c(76:90,92,93,95:100)] # removing 91 and 94, the problematic trees,
mtrees.ardtr<-make.simmap(trees,x,model=ard.tr,Q='empirical',nsim=1000,message=TRUE,prior=list(alpha=1,beta=1e4)) 
save(mtrees.ardtr,file="phylo/mtrees-mpardtr.RData")

##########################
## Analysis of mappings ##
##########################

describer<-function(simtree){
  simtree<-simtree[!unlist(lapply(simtree,is.null))] # remove missing lists indices,
  treelength<-length(simtree)
  metrics<-list()
  for(i in 1:(treelength/1000)){
    ii<-seq(1,treelength,1000)
    j<-seq(1000,treelength,1000)
    print(i)
    d1<-summary(simtree[ii[i]:j[i]])# for every real tree (100) there are 1000 'trees' or simulation replicates,
    # we want two things from each tree: the summary metrics simulation variances for transition counts and times for each state
    # get transition counts
    d1.count<-d1$count[,2:ncol(d1$count)]
    # get times,
    d1.time<-d1$times[,1:3]
    d<-data.table(cbind(d1.count,d1.time))
    d$tree<-i
    metrics[[i]]<-d
  }
  return(metrics)
}

# Run describer() on each mtrees, because this is compuationally exhaustive, each of these blocks were run as their own jobs.
# er
load(file="phylo/mtrees-er.RData")
er.full<-describer(mtrees.er)
rm(mtrees.er)
paste(Sys.time(),": finished ER model")
save(er.full,file="phylo/er-desc.RData")

#sym
load(file="phylo/mtrees-sym")
sym.full<-describer(mtrees.sym)
rm(mtrees.sym)
paste(Sys.time(),": finished SYM model")
save(sym.full,file="phylo/sym-desc.RData")

# ard
load(file="phylo/mtrees-ard.RData")
ard.full<-describer(mtrees.ard)
rm(mtrees.ard)
paste(Sys.time(),": finished ARD model")
save(ard.full,file="phylo/ard-desc.RData")

# ertr
load(file="phylo/mtrees-er-tr.RData")
ertr.full<-describer(mtrees.ertr)
rm(mtrees.ertr)
paste(Sys.time(),": finished ER-TR model")
save(ertr.full,file="phylo/ertr-desc.RData")

# ardtr
load(file="phylo/mtrees-ardtr.RData")
ardtr.full<-describer(mtrees.ardtr)
rm(mtrees.ardtr)
paste(Sys.time(),": finished ARD-TR model")
save(ardtr.full,file="phylo/ardtr-desc.RData")

# trtr
load(file="phylo/mtrees-trtr.RData")
tr.tr.full<-describer(mtrees.trtr)
rm(mtrees.tr)
paste(Sys.time(),": finished tr.tr")
save(tr.tr.full,file="phylo/trtr-desc.RData")

# load descibed datasets,
load("phylo/er-desc.RData")
load("phylo/sym-desc.RData")
load("phylo/ard-desc.RData")
load("phylo/ertr-desc.RData")
load("phylo/ardtr-desc.RData")
load("phylo/trtr-desc.RData")

models<-list(er.full,sym.full,ard.full,er.tr.full,ardtr.full,tr.tr.full)
models<-lapply(models,function(x) rbindlist(x))
mname<-c("er","sym","ard","ertr","ardtr","trtr")
models<-lapply(1:length(models),function(x) cbind(models[[x]],rep(mname[x],nrow(models[[x]]))))
models<-rbindlist(models)
save(models,file="phylo/models-full.RData")


# Now compile maximum parsimony models,
load("phylo/MP-ER-desc.RData")
load("phylo/MP-SYM-desc.RData")
load("phylo/MP-ARD-desc.RData")
load("phylo/mp-ertr-desc.RData") 
load("phylo/MP-ARDTR-descs.RData")
load("phylo/MP-trtr-desc.RData")


models2<-list(mp.er.full,mp.sym.full,mp.ard.full,mp.ertr.full,mp.ardtr.full,ptr.tr.full)
models2<-lapply(models2,function(x) rbindlist(x))
mname2<-c("mper","mpsym","mpard","mpertr","mpardtr","mptrtr")
models2<-lapply(1:length(models2),function(x) cbind(models2[[x]],rep(mname2[x],nrow(models2[[x]]))))
models2<-rbindlist(models2)
save(models2,file="phylo/models-full-MP.RData") # a temp version for now

