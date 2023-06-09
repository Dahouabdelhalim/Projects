# This script estimates the number of reversions that occur on our tree.

setwd("/scratch/osalehzi/phy/Final/")

# Load Packages
require(data.table)
require(phytools)
require(stringr)
require(plyr)

load('phylo/pre-loading.RData')

# Reversion fuction,
reversion<-function(t,estimates,x){
  # Step 1: get tips with a non-WL state,
  ts<-which(names(x) %in% names(x[x %in% c("A","A+B")]))
  dest<-data.table(1:t$Nnode+Ntip(t),estimates)
  # Step 2: for each non-WL tip, check every ancestor until the root to look for a high support WL node
  revs<-data.table(tip='',node='')
  for(i in 1:length(ts)){
    j<-ts[i]
    while(j!=(Ntip(t)+1)){
      j<-getParent(t,j) # get parent
      l<-dest[V1==j][,2:4] # subset probabilities
      #m<-names(l)[which(l==max(l))] # get state of MAX probability
      m<-names(l)[which(l>0.6)] # get state of higher than alpha threshold
      if(length(m)>0){
        if(m=="B"){
          revs<-rbind(revs,data.table(tip=ts[i],node=j))
        }}
    }
}
revs<-revs[!1,]
return(revs)
}


case.mcc<-FALSE
if(case.mcc){
  # MCC tree example
  revs.er<-reversion(tree,fitER$lik.anc,x)
  revs.sym<-reversion(tree,fitSYM$lik.anc,x)
  revs.ard<-reversion(tree,fitARD$lik.anc,x)
  revs.er.tr<-reversion(tree,fitER.tr$lik.anc,x)
  revs.ard.tr<-reversion(tree,fitARD.tr$lik.anc,x)
  revs.tr.tr<-reversion(tree,fitTR.tr$lik.anc,x)
  
  fin.er<-data.table(names(x)[as.numeric(unique(revs.er$tip))],rep("ER",length(as.numeric(unique(revs.er$tip)))))
  fin.sym<-data.table(names(x)[as.numeric(unique(revs.sym$tip))],rep("SYM",length(as.numeric(unique(revs.sym$tip)))))
  fin.ard<-data.table(names(x)[as.numeric(unique(revs.ard$tip))],rep("ARD",length(as.numeric(unique(revs.ard$tip)))))
  fin.er.tr<-data.table(names(x)[as.numeric(unique(revs.er.tr$tip))],rep("ER.tr",length(as.numeric(unique(revs.er.tr$tip)))))
  fin.ard.tr<-data.table(names(x)[as.numeric(unique(revs.ard.tr$tip))],rep("ARD.tr",length(as.numeric(unique(revs.ard.tr$tip)))))
  fin.tr.tr<-data.table(names(x)[as.numeric(unique(revs.tr.tr$tip))],rep("TR.TR",length(as.numeric(unique(revs.tr.tr$tip)))))
  fin<-data.table(rbind(fin.er,fin.sym,fin.ard,fin.er.tr,fin.ard.tr,fin.tr.tr))
  colnames(fin)<-c("rev.taxa","model")
  
  # just checking stuff,
  lapply(unique(fin$model),function(x) length(fin[model==x]$rev.taxa))
  pinetree<-count(table(fin$rev.taxa))
  pinetree[order(pinetree$x.Freq,decreasing = T),]
  
  # Plot the mcc run
  estimates<-fitER$lik.anc
  plotTree(tree,fsize=0.5)
  #
  colors<-setNames(c("blue","red","green3"),colnames(estimates))
  tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=colors,cex=0.3)
  nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.3,piecol=colors)
  
  nulo<-sapply(fin.er$V1,add.arrow,tree=tree,col="red",hedl=1,
               arrl=10,lwd=1,lend=1,offset=1)
  

  # Let's zoom into our various groups to make sure the counts are accurate.
  # for now, I'll use getnode() to get our zoom ins, 3 I think will work. 
  
  # Run these each separate after plotting. We don't even need the reconstructions to tell!
  
  # Case 1: node 881
  tips<-names(x[getDescendants(tree,881)[getDescendants(tree,881)<Ntip(tree)]])
  zoom(tree,tips)
  x2<-x[which(tree$tip.label %in% tips)]
  tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
  
  # Case 2: node 699
  tips<-names(x[getDescendants(tree,699)[getDescendants(tree,699)<Ntip(tree)]])
  zoom(tree,tips)
  x2<-x[which(tree$tip.label %in% tips)]
  tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
  
  # Case 3: node 592
  tips<-names(x[getDescendants(tree,592)[getDescendants(tree,592)<Ntip(tree)]])
  zoom(tree,tips)
  x2<-x[which(tree$tip.label %in% tips)]
  tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
  
  # Case 4: node 467
  tips<-names(x[getDescendants(tree,467)[getDescendants(tree,467)<Ntip(tree)]])
  zoom(tree,tips)
  x2<-x[which(tree$tip.label %in% tips)]
  tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
  
  # Figuring out the aphis node better,
  tips<-names(x[getDescendants(tree,596)[getDescendants(tree,596)<Ntip(tree)]])
  zoom(tree,tips)
  x2<-x[which(tree$tip.label %in% tips)]
  tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
  
  # get the internal nodes from estimates 
  a.nodes<-findMRCA(tree,which(tree$tip.label %in% tips))-1
  a.chillin<-lapply(a.nodes,function(x) c(x,getDescendants(tree,x)[getDescendants(tree,x)>Ntip(tree)]))
  a.chillin<-lapply(a.chillin,function(x) sort(x))
  a.estimate<-data.table(estimates)
  a.est2<-lapply(a.chillin,function(x) estimates[x-Ntip(tree),])
  a.estz<-as.matrix(ldply(a.est2,cbind))
  nodelabels(pie=a.estz,cex=0.4,piecol=colors)
  
  # figuring out the cavariella node better (697)
  tips<-names(x[getDescendants(tree,697)[getDescendants(tree,697)<Ntip(tree)]])
  zoom(tree,tips)
  x2<-x[which(tree$tip.label %in% tips)]
  tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
  
}

# I guess we should try and run this on all our different trees still,

errorhand<-function(x){
  tryCatch(revertortron(x),
           warning = function(w) {print("warning "); revertortron(x)},
           error = function(e) {print("error "); NaN}) 
}

revertortron<-function(tr){
  x<-tr$tip.label
  
  print(paste("fitting models:",Sys.time(),sep=""))
  fitER<-ace(x,tr,type="discrete",model="ER")$lik.anc
  print(paste("fitting models:",Sys.time(),sep=""))
  fitSYM<-ace(x,tr,type="discrete",model="SYM")$lik.anc
  fitARD<-ace(x,tr,type="discrete",model="ARD")$lik.anc
  fitER.tr<-ace(x,tr,type="discrete",model=er.tr)$lik.anc
  fitARD.tr<-ace(x,tr,type="discrete",model=ard.tr)$lik.anc
  fitTR.tr<-ace(x,tr,type="discrete",model=tr.tr)$lik.anc
  
  print(paste("starting reversions:",Sys.time(),sep=""))
  revs.er<-reversion(tree,fitER$lik.anc,x)
  revs.sym<-reversion(tree,fitSYM$lik.anc,x)
  revs.ard<-reversion(tree,fitARD$lik.anc,x)
  revs.er.tr<-reversion(tree,fitER.tr$lik.anc,x)
  revs.ard.tr<-reversion(tree,fitARD.tr$lik.anc,x)
  revs.tr.tr<-reversion(tree,fitTR.tr$lik.anc,x)
  
  print(paste("done with reversions:",Sys.time(),sep=""))
  fin.er<-data.table(names(x)[as.numeric(unique(revs.er$tip))],rep("ER",length(as.numeric(unique(revs.er$tip)))))
  fin.sym<-data.table(names(x)[as.numeric(unique(revs.sym$tip))],rep("SYM",length(as.numeric(unique(revs.sym$tip)))))
  fin.ard<-data.table(names(x)[as.numeric(unique(revs.ard$tip))],rep("ARD",length(as.numeric(unique(revs.ard$tip)))))
  fin.er.tr<-data.table(names(x)[as.numeric(unique(revs.er.tr$tip))],rep("ER.tr",length(as.numeric(unique(revs.er.tr$tip)))))
  fin.ard.tr<-data.table(names(x)[as.numeric(unique(revs.ard.tr$tip))],rep("ARD.tr",length(as.numeric(unique(revs.ard.tr$tip)))))
  fin.tr.tr<-data.table(names(x)[as.numeric(unique(revs.tr.tr$tip))],rep("TR.TR",length(as.numeric(unique(revs.tr.tr$tip)))))
  fin<-data.table(rbind(fin.er,fin.sym,fin.ard,fin.er.tr,fin.ard.tr,fin.tr.tr))
  colnames(fin)<-c("rev.taxa","model")
  
  return(fin)
}

dodo<-FALSE
if(dodo){
  revses<-list(list())
  for(i in 1:length(trees)){
    print(paste(Sys.time(),"starting: ",i))
    revses[[i]]<-errorhand(trees[[i]])
  }
  save.image(file="revs/all-trees.RData")
}

load("revs/all-trees.RData")

# add tree id,
revs<-data.table(ldply(1:length(revses),function(x) cbind(revses[[x]],tree=rep(x,nrow(revses[[x]])))))

# Find consensus reversing taxa,
r<-data.table(ldply(unique(revs$tree),function(x) cbind(count(table(revs[tree==x]$rev.taxa)),tree=x)))
rr<-data.table(ldply(unique(r$tree),function(x) r[tree==x & x.Freq>5]))
rrr<-data.table(count(table(rr$x.Var1))) # levels of factor account for the 0 frequencies
r4<-rrr[x.Freq>0][order(rrr[x.Freq>0]$x.Freq,decreasing = T),][,1:2]
colnames(r4)<-c("species","freq.trees")
r4

# looking at only best fit model reversions,
q<-data.table(count(table(revs[model=="ARD.tr"]$rev.taxa)))
q[order(q$x.Freq,decreasing = T)]

save(r4,file="revs/consensus.RData")
