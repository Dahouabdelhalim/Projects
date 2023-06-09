# This script generates our various plots, in rough form for most, to be then editted in Illustrator as figures.
require(phytools)
setwd("/scratch/osalehzi/phy/Final/")
load(file = "phylo/pre-loading.RData")

# Generate MCC phylogenies with ancestral reconstructions for all models
# This series will be the raw, all nodes plots for the models ER,SYM,ARD,ER.tr,ARD.tr,TR.tr (6 plots)
if(FALSE){
  pdf("plots/ER.pdf")
  estimates<-fitER$lik.anc
  plotTree(tree,ftype="i",lwd=1,fsize=1.2)
  #
  colors<-setNames(c("blue","red","green3"),colnames(estimates))
  nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.2,piecol=colors)
  dev.off()
  
  pdf("plots/SYM.pdf")
  estimates<-fitSYM$lik.anc
  plotTree(tree,ftype="i",lwd=1,fsize=1.2)
  #
  colors<-setNames(c("blue","red","green3"),colnames(estimates))
  nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.2,piecol=colors)
  dev.off()
  
  pdf("plots/ARD.pdf")
  estimates<-fitARD$lik.anc
  plotTree(tree,ftype="i",lwd=1,fsize=1.2)
  #
  colors<-setNames(c("blue","red","green3"),colnames(estimates))
  nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.2,piecol=colors)
  dev.off()
  
  pdf("plots/fitARD-tr.pdf")
  estimates<-fitARD.tr$lik.anc
  plotTree(tree,ftype="i",lwd=1,fsize=1.2)
  #
  colors<-setNames(c("blue","red","green3"),colnames(estimates))
  nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.2,piecol=colors)
  dev.off()
  
  pdf("plots/fitER-tr.pdf")
  estimates<-fitER.tr$lik.anc
  plotTree(tree,ftype="i",lwd=1,fsize=1.2)
  #
  colors<-setNames(c("blue","red","green3"),colnames(estimates))
  nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.2,piecol=colors)
  dev.off()
  
  pdf("plots/fitTR-tr.pdf")
  estimates<-fitTR.tr$lik.anc
  plotTree(tree,ftype="i",lwd=1,fsize=1.2)
  #
  colors<-setNames(c("blue","red","green3"),colnames(estimates))
  nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.2,piecol=colors)
  dev.off()
}
# Plotting Figure 1, MCC tree with >90% nodes and ecological traits,

# This script is an advanced plotting of phylogenies with multiple tip states for different traits. This is built off of Liam's phytools page: http://blog.phytools.org/2018/02/annotating-plotted-phylogeny-with.html

require(data.table)
library(RColorBrewer)
require(phytools)
# construct trait matrix,

# modifying treeName() to customize traits obtained,
treeName.multi<-function(tips){ # function to get names of a tree
  st<-matrix(unlist(lapply(tips,function(x) c(data[species==gsub("_"," ",x)]$wings,data[species==gsub("_"," ",x)]$host.breadth,data[species==gsub("_"," ",x)]$ant,data[species==gsub("_"," ",x)]$woodiness))),ncol = 4,byrow = T)
  st<-data.frame(st)
  row.names(st)<-tips
  return(st)
}
treeName<-function(tips){ # function to get names of a tree
  st<-unlist(lapply(tips,function(x) data[species==gsub("_"," ",x)]$wings))
  sta<-setNames(as.factor(st),tips)
  return(sta)
}

# Curate data
phyl<-read.tree(file="input/Aphidoidea-MCC.tre")
phy.names<-fread("input/phy-names.csv",header = T,drop = "V1") # this matches the order of the mcc tre 
data<-fread("input/aphiddata.csv",header=T,drop="V1") # male data

# gotta make sure that we plot the unknowns, then we can specify the conversions in the methods
moda<-data
moda$wings<-as.numeric(gsub("dimorphic",2,gsub("apterous","0",gsub("alate","1",data$wings))))                  # wings: 0=WL, 1=W, 2=both, NA=unknown   {4}
moda$host<-as.numeric(gsub("hetero",1,gsub("mono",0,moda$host.alt)))                                           # host.alt: 0=mono,1=hetero,NA=unknown   {3}
moda$gyno<-as.numeric(gsub("gyno",0,gsub("sexu",1,moda$`gyn-sex`)))                                            # host.gyno: 0=gyno,1=sexu,NA=unknown  {3}
moda$woody<-as.numeric(gsub("variable",2,gsub("W",1,gsub("H",0,moda$host.type))))                              # wood: 0=herbaceous,1=woody,2=variable  {4}
moda$ant<-as.numeric(gsub("not",0,gsub("ant",1,moda$ant)))                                                     # ant: 0=not,1=tended,NA=unknown         {3}
moda$breadth<-as.numeric(gsub("polyphagous",2,gsub("specialist",0,gsub("restricted",1,moda$host.breadth))))    # breadth: 0=specialist,1=restricted,2=polyphagous,NA=unknown {4}
moda$fem<-as.numeric(gsub("ala",2,gsub("apt",1,gsub("nei",0,gsub("both",3,moda$polyph)))))                     # fem: 0=neither, 1=apterous, 2=alate, 3=both, NA=unknown   {5}

# make db discrete for plotting, seq(1,40,4) is 10 increments to bin into
require(OneR)
moda$db<-bin(moda$db,11,na.omit = F)

data<-moda[,c("species","wings","host","gyno","db","woody","ant","breadth","fem","family")]                    # family:                                {10}


phyl$tip.label<-phy.names$tip_labels
species<-data[wings!=""]$species
tree<-drop.tip(phyl,phyl$tip.label[-match(species,phy.names$tips)])
x<-treeName(tree$tip.label)

# rotate node, 762 and 807 to get the aphid families grouped properly
tree<-rotateNodes(tree,762)
tree<-rotateNodes(tree,807)

# use best fit model,
estimates<-fitARD.tr$lik.anc

d<-data
d$species<-gsub(" ","_",d$species)

plotter<-function(tree){
  par(family="serif") # set font
  
  tips<-tree$tip.label
  st<-data.frame(d[match(tips,d$species)])
  st$species<-NULL
  st$phyl<-NULL
  #st$family<-NULL
  row.names(st)<-tips
  st$wings<-as.factor(st$wings)
  st$host<-as.factor(st$host)
  st$gyno<-as.factor(st$gyno)
  st$db<-as.factor(st$db) # i'll manually make the gradient look continuous in Illustrator by scaling gradient scale numerically ugh
  st$woody<-as.factor(st$wood)
  st$ant<-as.factor(st$ant)
  st$breadth<-as.factor(st$breadth)
  st$fem<-as.factor(st$fem)
  st$family<-as.factor(st$family)
  
  X<-st
  
  #colnames(X)<-c("wings","family","alt","db","wood","ant")
  tree<-reorder(tree,"cladewise")
  X<-X[tree$tip.label,]
  plotTree(tree,plot=FALSE)
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  plotTree(tree,lwd=1,ylim=c(0,obj$y.lim[2]*1.05),xlim=c(0,obj$x.lim[2]*1.3),
           ftype="off")
  nodes<-(1:tree$Nnode+Ntip(tree))[apply(estimates,1,function(x,alpha) any(x>(1-alpha)),alpha=0.1)]
  colors<-setNames(c("blue","red","green3"),colnames(estimates))
  nodelabels(node=nodes,pie=estimates[nodes-Ntip(tree),],piecol=colors,cex=0.2)
  #nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.2,piecol=colors)
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  h<-max(obj$xx)
  fsize<-0.3
  for(i in 1:Ntip(tree)){
    lines(c(obj$xx[i],h),rep(obj$yy[i],2),lty="dotted")
    text(h,obj$yy[i],tree$tip.label[i],cex=fsize,pos=4,font=3,offset=0.1)
  }
  s<-max(fsize*strwidth(tree$tip.label))
  start.x<-1.05*h+s
  palettes<-c("OrRd","PuOr","RdYlBu","Paired","Spectral","OrRd",c("Paired","Set3"),"Set3") # add more colors to match our sets
  cols<-list()
  for(i in 1:ncol(X)){ 
    text(start.x,max(obj$yy)+1,paste(colnames(X)[i]),pos=4,srt=60,
         cex=0.4,offset=0)
    cols[[i]]<-setNames(sample(c(brewer.pal(max(3,length(levels(X[[i]]))),
                                            palettes[i]),"black"),length(levels(X[[i]]))),levels(X[[i]]))
    for(j in 1:nrow(X)){
      xy<-c(start.x,obj$yy[j])
      y<-c(xy[2]-0.5,xy[2]+0.5,xy[2]+0.5,xy[2]-0.5)
      asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
        par()$pin[2]/par()$pin[1]
      x<-c(xy[1]-0.5*asp,xy[1]-0.5*asp,xy[1]+0.5*asp,xy[1]+0.5*asp)
      polygon(x,y,col=cols[[i]][as.character(X[[i]][j])])
    }
    start.x<-start.x+2*asp
  }
  start.y<-max(obj$yy)
  for(i in 1:ncol(X)){
    text(start.x,start.y,paste("trait",colnames(X)[i]),pos=4,cex=0.5,
         offset=3)
    add.simmap.legend(colors=cols[[i]],shape="square",prompt=FALSE,
                      x=start.x,y=start.y-2*strheight("W")*0.9,fsize=0.5)
    start.y<-start.y-1.5*0.9*strheight("W")*(length(cols[[i]])-1)-6
  }
}
pdf("plots/Figure1.pdf")
plotter(tree)
dev.off()

# Plotting Figure 2, stochastic simulations (and supplementary simulation plots)
require(ggplot2)
require(data.table)
require(plyr)
require(phytools)
setwd("/scratch/osalehzi/phy/Final/")
load("phylo/models-full.RData")
source("misc/multiplot.R")
# here's a basic plotter for all the data,
plotter<-function(model,tree.ind,trait,poly=FALSE){
  d.count<-models[V2==model & tree %in% tree.ind][,!c("A","B","A+B","A,B","B,A","A+B,A","A+B,B")]
  if(poly){
    # poly ==TRUE includes the A+B,A and A+B,B polymorphism loss transitions, 
    # only useful in the ardtr model; used for figure s8
    d.count<-models[V2==model & tree %in% tree.ind][,!c("A","B","A+B","A,B","B,A")]
  }
  m.count<-melt(d.count,id.vars = c("tree","V2"))
  mean.count<-ddply(m.count,"variable",summarise,grp.mean=mean(value))
  cnt<-ggplot(m.count[tree %in% tree.ind],aes(value,color=variable))+geom_density(alpha=0.5)+
    geom_vline(data=mean.count,aes(xintercept=grp.mean,color=variable),linetype="dashed")+theme_classic()+facet_grid(.~tree)
  if(trait=="count"){
    return(cnt)}
  # plot times
  d.times<-models[V2==model & tree==tree.ind][,c("A","B","A+B","tree")]
  # get proportion of time, not total;
  d.times<-data.table(t(apply(d.times[,1:3],1,function(x) x/sum(x))),tree=d.times$tree)
  m.times<-melt(d.times,id.vars = c("tree"))
  mean.times<-ddply(m.times,"variable",summarise,grp.mean=mean(value))
  tim<-ggplot(m.times,aes(value,color=variable))+geom_density(alpha=0.5)+
    geom_vline(data=mean.times,aes(xintercept=grp.mean,color=variable),linetype="dashed")+theme_classic()
  if(trait=="times"){
    return(tim)}
}

# best fit model, 
p1<-plotter("trtr",1,"count")
q1<-plotter("trtr",1,"times")

pdf('plots/fig2.pdf')
multiplot(p1,q1)
dev.off()

# plotting maximum parsimony models,
models<-models2
p2<-plotter("mptrtr",1,"count")
q2<-plotter("mptrtr",1,"times")
gridExtra::grid.arrange(p2,q2)

# plotting fig s8, where the ardtr poly gain AND loss are shown for count means,
pdf('plots/figs8.pdf')
plotter("ardtr",1,"count",TRUE) # for tree 1 of 100,
dev.off()

# here's a neato way of getting our supplementary tree means from simulations to summarize figure 2 primary,
suppplotter<-function(model,trait){
  # this is for counts,
  #model="mptrtr"
  a<-models[V2==model][,!c("A","B","A+B","A,B","B,A","A+B,A","A+B,B")]
  x<-data.table()
  for(i in unique(a$tree)){
    b<-a[tree==i]
    c<-melt(b,id.vars = c("tree","V2"))
    d<-ddply(c,"variable",summarise,grp.mean=mean(value))
    e<-data.table(d)
    e$i<-i
    colnames(e)<-c("dir","mean","tree")
    x<-rbind(x,e)
  }
  
  # this is for times,
  a<-models[V2==model][,c("A","B","A+B","tree")]
  a<-data.table(t(apply(a[,1:3],1,function(x) x/sum(x))),tree=a$tree)
  a<-a[`A+B`<.2] # filters out outliers, run without this to ensure outliers are consistent
  y<-data.table()
  for(i in unique(a$tree)){
    b<-a[tree==i]
    c<-melt(b,id.vars = c("tree"))
    d<-ddply(c,"variable",summarise,grp.mean=mean(value))
    e<-data.table(d)
    e$i<-i
    colnames(e)<-c("dir","mean","tree") # i'm keeping dir for now instead of pro
    y<-rbind(y,e)
  }
  
  if(model %in% c("ertr","ardtr","mpertr","mpardtr")){
    bads<-a[`A+B`>0.3]$tree
    x<-x[!tree %in% bads]
    y<-y[!tree %in% bads]
  }
  
  p<-ggplot(x[mean<1000], aes(x=dir, y=mean, group=tree)) +
    geom_point(aes(colour=dir), size=3, position=position_dodge(width=0.1)) +
    geom_line(size=1, alpha=0.3, position=position_dodge(width=0.1))+xlab('Wing gain/loss') +
    ylab('Mean count value ') +
    scale_colour_manual(values=c("#009E73", "#D55E00"), guide=FALSE) + 
    theme_bw()
  q<-ggplot(y, aes(x=dir, y=mean, group=tree)) +
    geom_point(aes(colour=dir), size=2, position=position_dodge(width=0.1)) +
    geom_line(size=.5, alpha=0.1, position=position_dodge(width=0.1))+xlab('Wing gain/loss') +
    ylab('Mean count value ') +
    scale_colour_manual(values=c("#009E73", "#D55E00","#7D2E68"), guide=FALSE) + 
    theme_bw()
  
  if(trait=="count"){
    return(p)
  }
  if(trait=="times"){
    return(q)
  }
}

# supplementary figures, 
p1<-suppplotter("er","count")
q1<-suppplotter("er","times")
p2<-suppplotter("sym","count")
q2<-suppplotter("sym","times")
p3<-suppplotter("ard","count")
q3<-suppplotter("ard","times")
p4<-suppplotter("ertr","count")
q4<-suppplotter("ertr","times")
p5<-suppplotter("ardtr","count")
q5<-suppplotter("ardtr","times")
p6<-suppplotter("trtr","count")
q6<-suppplotter("trtr","times")

pdf('plots/figs2.pdf')
multiplot(p1,q1,p4,q4,p2,q2,p5,q5,p3,q3,p6,q6,cols = 3)
dev.off()

sp1<-suppplotter("mper","count")
sq1<-suppplotter("mper","times")
sp2<-suppplotter("mpsym","count")
sq2<-suppplotter("mpsym","times")
sp3<-suppplotter("mpard","count")
sq3<-suppplotter("mpard","times")
sp4<-suppplotter("mpertr","count")
sq4<-suppplotter("mpertr","times")
sp5<-suppplotter("mpardtr","count")
sq5<-suppplotter("mpardtr","times")
sp6<-suppplotter("mptrtr","count")
sq6<-suppplotter("mptrtr","times")

pdf('plots/figs3.pdf')
multiplot(sp1,sq1,sp4,sq4,sp2,sq2,sp5,sq5,sp3,sq3,sp6,sq6,cols=3)
dev.off()

p1
ggsave("plots/figS2-1.pdf")
p2
ggsave("plots/figS2-2.pdf")
p3
ggsave("plots/figS2-3.pdf")
p4
ggsave("plots/figS2-4.pdf")
p5
ggsave("plots/figS2-5.pdf")
p6
ggsave("plots/figS2-6.pdf")

# figure s8, including polymorphism loss with gain,
suppplotters8<-function(){
  # this is for counts,
  model="ardtr" # only model this works for
  a<-models[V2==model][,!c("A","B","A+B","A,B","B,A")] #,"A+B,A","A+B,B")]
  x<-data.table()
  for(i in unique(a$tree)){
    b<-a[tree==i]
    c<-melt(b,id.vars = c("tree","V2"))
    d<-ddply(c,"variable",summarise,grp.mean=mean(value))
    e<-data.table(d)
    e$i<-i
    colnames(e)<-c("dir","mean","tree")
    x<-rbind(x,e)
  }
  
  p<-ggplot(x[mean<1000], aes(x=dir, y=mean, group=tree)) +
    geom_point(aes(colour=dir), size=3, position=position_dodge(width=0.1)) +
    geom_line(size=1, alpha=0.3, position=position_dodge(width=0.1))+
    xlab('Polymorphism gain/loss') +
    ylab('Mean count value ') +
    #scale_colour_manual(values=c("#009E73", "#D55E00"), guide=FALSE) + 
    theme_bw()
  
    return(p)
}

pdf("plots/figs8b.pdf")
suppplotters8()
dev.off()

# Plotting SecSSE permutation tests,
require(plyr)
require(ggplot2)

etd<-readRDS("musse/SecSSE/ETD.RDS")
b<-etd$ML

# sloppy way of getting it
files<-grep("100.*RDS",list.files("musse/SecSSE/"),value=T)
rand<-lapply(files,function(x) readRDS(paste("musse/SecSSE/",x,sep='',collapse='')))
a<-ldply(do.call(c,rand),function(x) x$ML)
a

ran.plot<-ggplot(a,aes(V1))+geom_histogram(aes(y=..density..),binwidth=2,color="black",fill="white")+geom_density(alpha=0.2,fill="#FF6666")+geom_vline(aes(xintercept=b),color="blue",linetype="dashed")+theme_classic()+xlim(-2980,-2810)+
  ylab("Density")+xlab("Maximum Likelihoods")

pdf("plots/secsse-etd-permutes.pdf")
ran.plot
dev.off()

################################
# Plot Reversions ##############
################################
load("revs/consensus.RData")
con.sp<-r4[freq.trees>50]$species
estimates<-fitARD.tr$lik.anc

# for nodes 699, 467, 596, and 881

# aphis
pdf("plots/revs-aphis.pdf")
plotTree(tree,fsize=0.5)
#
colors<-setNames(c("blue","red","green3"),colnames(estimates))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=colors,cex=0.3)
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.3,piecol=colors)
tips<-names(x[getDescendants(tree,596)[getDescendants(tree,596)<Ntip(tree)]])
zoom(tree,tips)
x2<-x[which(tree$tip.label %in% tips)]
tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
dev.off()

# brachycaudus
pdf("plots/revs-brachy.pdf")
plotTree(tree,fsize=0.5)
#
colors<-setNames(c("blue","red","green3"),colnames(estimates))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=colors,cex=0.3)
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.3,piecol=colors)

tips<-names(x[getDescendants(tree,467)[getDescendants(tree,467)<Ntip(tree)]])
zoom(tree,tips)
x2<-x[which(tree$tip.label %in% tips)]
tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
dev.off()

# longistigma
pdf("plots/revs-longistigma.pdf")
plotTree(tree,fsize=0.5)
#
colors<-setNames(c("blue","red","green3"),colnames(estimates))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=colors,cex=0.3)
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.3,piecol=colors)

tips<-names(x[getDescendants(tree,881)[getDescendants(tree,881)<Ntip(tree)]])
zoom(tree,tips)
x2<-x[which(tree$tip.label %in% tips)]
tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
dev.off()

# cavariella
pdf("plots/revs-cavariella.pdf")
plotTree(tree,fsize=0.5)
#
colors<-setNames(c("blue","red","green3"),colnames(estimates))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=colors,cex=0.3)
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.3,piecol=colors)

tips<-names(x[getDescendants(tree,699)[getDescendants(tree,699)<Ntip(tree)]])
zoom(tree,tips)
x2<-x[which(tree$tip.label %in% tips)]
tiplabels(pie=to.matrix(x2,sort(unique(x2))),piecol=colors,cex=0.4)
dev.off()


if(FALSE){# get the internal nodes from estimates 
  a.nodes<-findMRCA(tree,which(tree$tip.label %in% tips))-1
  a.chillin<-lapply(a.nodes,function(x) c(x,getDescendants(tree,x)[getDescendants(tree,x)>Ntip(tree)]))
  a.chillin<-lapply(a.chillin,function(x) sort(x))
  a.estimate<-data.table(estimates)
  a.est2<-lapply(a.chillin,function(x) estimates[x-Ntip(tree),])
  a.estz<-as.matrix(ldply(a.est2,cbind))
  nodelabels(pie=a.estz,cex=0.4,piecol=colors)
}

# Figure S2: Reversions, brachycaudus tree,
load("revs/brach.RData")
tree<-brach[[1]]
x<-brach[[2]]
estimates<-brach[[3]]

pdf('plots/S2.pdf')
plotTree(tree,ftype="i",lwd=1,fsize=.6)
#
colors<-setNames(c("green3","blue","red","black"),colnames(estimates))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=colors,cex=.2)
nodelabels(node=1:tree$Nnode+Ntip(tree),pie=estimates,cex=0.2,piecol=colors)
dev.off()

# Figure S4: Co-distributions of the male female wing states,
load("models/pagels/pagels-len.RData")
load("phylo/pre-loading.RData") # to get the proper tree format
x<-as.factor(comp.dat$wings)
y<-as.factor(comp.dat$host)
names(x)<-names(y)<-comp.dat$animal

pdf("plots/S4-len.pdf")
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(x[tree$tip.label],c("1","0")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("1","0")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(y[tree$tip.label],c("1","0")),piecol=c("blue","red"),
          cex=0.3)
dev.off()

load("models/pagels/pagels-str.RData")
load("phylo/pre-loading.RData") # to get the proper tree format
x<-as.factor(comp.dat$wings)
y<-as.factor(comp.dat$host)
names(x)<-names(y)<-comp.dat$animal

pdf("plots/S4-str.pdf")
par(mfrow=c(1,2))
plot(tree,show.tip.label=FALSE,no.margin=TRUE)
par(fg="transparent")
tiplabels(pie=to.matrix(x[tree$tip.label],c("1","0")),piecol=c("blue","red"),
          cex=0.3)
par(fg="black")
add.simmap.legend(colors=setNames(c("blue","red"),c("1","0")),prompt=FALSE,
                  x=0,y=10,fsize=1)
par(fg="transparent")
plot(tree,show.tip.label=FALSE,no.margin=TRUE,direction="leftwards")
tiplabels(pie=to.matrix(y[tree$tip.label],c("1","0")),piecol=c("blue","red"),
          cex=0.3)
dev.off()
