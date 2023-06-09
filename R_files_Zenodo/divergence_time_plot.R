#this script generates the time-calibrated phylogeny in Figure 3
#import nexus file from figtree in R
library(ape)
library(phylotools)
library(phytools)
library(phylotate)
library(plotrix)
library(diagram)

#in FigTree, you MUST:
#date the root to desired time (1 Mya)
#reverse axis
#export nexus file with annotations

#read nexus file
setwd("~/Documents")
df <- read_annotated("May3tree.nex", format="nexus")
highlight <- df$node.comment[167+5]
df <- read_annotated("collapsed_tree.nex", format="nexus")

#function to relabel tips
foo<-function(df){
  fsize<-80*par()$pin[2]/par()$pin[1]/Ntip(df)
  plotTree(df,fsize=fsize,lwd=1)
}

#plot tree
plotTree(df)
foo(df)

#df$edge is dataframe of all branches connecting nodes
#numbers 1-62 are the terminal taxa, in order of df$tip.label
#numbers 63 and above are nodes in the tree (==node.label+62)

#inspect node labels
df <- makeNodeLabel(df, method = "number", prefix=NULL, nodeList = list())
plot(df, show.node.label = TRUE)

#generate CI table for mean, 95%HPD 
CI <- data.frame(lower=numeric(123), upper=numeric(123), stringsAsFactors = F)
for(i in 1:length(highlight))
{
  grab <- highlight[i]
  grab2 <- sub(".*height_95%_HPD=","",grab)
  grab3 <- unlist(strsplit(grab2, split="}"))[1]
  pile <- unlist(strsplit(grab3, split=""))
  pile <- pile[2:length(pile)]
  smash <- paste(pile, collapse="")
  bottom <- as.numeric(unlist(strsplit(smash, split=","))[1])
  top <- as.numeric(unlist(strsplit(smash, split=","))[2])
  CI[62+5,] <- c(bottom, top)
}

#insert 0 for NAs
#CI[is.na(CI)] <- 0
#CI <- sapply(CI, as.numeric)
#CI <- as.data.frame(CI)
#recall...just use the ancestors
#CI <- CI[168:nrow(CI),]

#identify scalar to convert all nodes according to 1 Mya root
scalar  <- 1000000/max(nodeHeights(df))

#only keep interesting divergence time(s) -- 5th node is the western ancestor
CI <- CI*scalar

## plot tree with error bars around divergence times at nodes
## written by Liam J. Revell 2017
#http://blog.phytools.org/2017/03/function-to-plot-tree-with-error-bars.html
tree <- df
tree <- ladderize(tree,right=T)
tree$Ntip <- 62
plotTree.errorbars<-function(tree,CI,...){
  args<-list(...)
  if(!is.null(args$gridlines)){ 
    gridlines<-args$gridlines
    args$gridlines<-NULL
  } else gridlines<-TRUE
  if(is.null(args$mar)) args$mar<-c(4.1,1.1,1.1,1.1)
  if(is.null(args$ftype)) args$ftype<-"i"
  fsize<-if(!is.null(args$fsize)) args$fsize else 1
  if(is.null(args$direction)) args$direction<-"leftwards"
  if(!is.null(args$bar.width)){
    bar.width<-args$bar.width
    args$bar.width<-NULL
  } else bar.width<-11
  if(!is.null(args$cex)){
    cex<-args$cex
    args$cex<-NULL
  } else cex<-1.2
  if(!is.null(args$bar.col)){
    bar.col<-args$bar.col
    args$bar.col<-NULL
  } else bar.col<-"blue"
  par(mar=args$mar)
  plot.new()      
  th<-max(nodeHeights(tree))
  h<-max(th,max(CI))
  if(is.null(args$xlim)){
    m<-min(min(nodeHeights(tree)),min(CI))
    d<-diff(c(m,h))
    pp<-par("pin")[1]
    sw<-fsize*(max(strwidth(tree$tip.label,units="inches")))+
      1.37*fsize*strwidth("W",units="inches")
    alp<-optimize(function(a,d,sw,pp) (a*1.04*d+sw-pp)^2,
                  d=d,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    args$xlim<-if(args$direction=="leftwards") c(h,m-sw/alp) else 
      c(m,h+sw/alp)
  }
  if(is.null(args$at)) at<-seq(0,h,by=h/5)
  else {
    at<-args$at
    args$at<-NULL
  }
  args$tree<-tree
  args$add<-TRUE
  do.call(plotTree,args=args)
  if(gridlines) abline(v=at,lty="dashed",
                       col=make.transparent("grey",0.5))
  #axis(1,at=at,labels=signif(at,3))
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    lines(x=c(CI[67,1],CI[67,2]),
          y=rep(obj$yy[67],2),lwd=4,lend=0,
          col=make.transparent(bar.col,0.4))
  points(obj$xx[67],
         obj$yy[67],pch=19,col=bar.col,cex=1)
}

#rescale tree edges
for(i in 1:length(tree$edge.length)){tree$edge.length[i] <- tree$edge.length[i]*scalar}

#plot the tree with error bars
stored <- tree$tip.label
tree$tip.label[c(1:length(tree$tip.label))] <- ""
plotTree.errorbars(tree,CI)
axis(side=1,at=c(1250000,1000000,750000,500000,250000,0),labels=c("1.25","1.0","0.75","0.50","0.25","0"))
mtext("Million years before present",1,at=625000,line=2.5)

###for bar-plotting structure results for each tree leaf###
#import .csv with admixture proportions, Q
#i made sure my names with "_" sort identically to Carly's IDs in excel
nom <- read.csv("names_k3.csv")

#write now the position in tree$tip.label for later sorting
tree$tip.label <- stored
nom$pos <- 0
for(i in 1:nrow(nom)){
nom[i,]$pos <- which(tree$tip.label==nom[i,]$plant)
}

#sort dataframe by nom$pos so that it comports with tree$tip.label
nom <- nom[order(nom$pos),]
q <- nom
rownames(q) <- NULL
rownames(q) <- q[,1]
q <- q[,-1]
q <- q[,-ncol(q)]

#phylo.barplot(phylo=tree,dat=q,SE=NULL,var.lab=c("cluster 1","cluster 2"),cex=0.7)
plotTree.barplot(tree,q,
args.barplot=list(col=c("firebrick2","darkorange2","gold1"), border=NA,space=0,xlim=c(0,1)),args.axis=list(at=seq(0,1,by=0.25)))
mtext("Q",1,at=0.5,line=2.5)

#from: http://blog.phytools.org/2017/01/plottreebarplot-with-more-user-options.html
plotTree.barplot(tree,q,tip.labels=F,args.barplot=list(col=c("blue","red"),space=0,xlim=c(0,1.3)),args.axis=list(at=seq(0,1,by=0.25)))
mtext("Q",1,at=0.5,line=2.5)