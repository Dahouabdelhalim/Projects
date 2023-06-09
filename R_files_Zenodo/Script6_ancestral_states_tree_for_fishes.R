###### Library ######

library(diversitree)
library(Matrix)
library(ggplot2)
library(colorRamps)
library(ggpubr)
library(extrafont)

###### Data preparation for Otophysi ######

load("phy_o_neot_YK2021.Robj")
load("phy_o_otop_YK2021.Robj")

###### Transformation matrix from state to chr and arm ######

chrmax<-35
snum<-(chrmax+1)*(chrmax+2)/2-1
stch.mat<-matrix(rep(0),nrow=snum,ncol=(chrmax*2))
star.mat<-matrix(rep(0),nrow=snum,ncol=(chrmax*2))
stch.vec<-c()
star.vec<-c()
rnum<-1
for(d in 1:chrmax){
  for(j in d:(2*d)){
    stch.vec<-c(stch.vec,d)
    star.vec<-c(star.vec,j)
    stch.mat[rnum,d]<-1
    star.mat[rnum,j]<-1
    rnum<-rnum+1
  }
}

##### Drawing Eurypterygii tree ######

## Font preparation
subset(fonttable(), FamilyName == "Arial") ##Check not more than one "Arial" in your system
loadfonts(device = "postscript")

## Figure drawing
postscript("Fig_tree_with_ancestral_states_Eurypterygii_XXXXXX.eps", height = 7.4, width = 7.4,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)

#parameters
ew<-0.2 # edge width
panel_label_par<-c(-150,140,-127,116,3) #parameters used for labeling A, B, C, D
par(oma=c(1.5,2,2,0)) #margin

#layout 4 panels
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T),
       widths=c(3.65,3.65),heights=c(3.65,3.65,0.2))

#color vector using matlab.like 
colvec<-matlab.like((chrmax*2-1)*10+1)

###### Eurypterygii tree ######

## M0

# Load the R object for ASR marginal of nodes in M0 model, st.nodes.neot
load("st_nodes_neot_M0_YK2021.Robj")

## Chromosome number

# Change to probability in chromosome number
ch.node.mat<-t(st.nodes.neot) %*% stch.mat

# For nodes
# Get mean chromosome number
nodestate<-as.vector(ch.node.mat %*% 1:70)
# Get the index of color
nodecolindex<-round(nodestate,1)*10-9
# Get the color name
nodecol<-colvec[nodecolindex]

# For tips
tipcol<-colvec[stch.vec[phy.o.neot$tip.state]*10-9]

# Plot
plot(phy.o.neot,"fan",no.margin=T,edge.width=0.4,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=0.6)
nodelabels(pch=16,col=nodecol,frame="n",cex=1)
text(panel_label_par[1],panel_label_par[2],labels="A",cex=panel_label_par[5])
mtext("Chromosome number evolution",side=3,cex=1)
mtext("Eurypterygii (M0)",side=2,cex=1)

## Arm number

ar.node.mat<-t(st.nodes.neot) %*% star.mat
nodestate<-as.vector(ar.node.mat %*% 1:70)
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[star.vec[phy.o.neot$tip.state]*10-9]

plot(phy.o.neot,"fan",no.margin=T,edge.width=0.4,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=0.6)
nodelabels(pch=16,col=nodecol,frame="n",cex=1)
text(panel_label_par[1],panel_label_par[2],labels="B",cex=panel_label_par[5])
mtext("Arm number evolution",side=3,cex=1)

## M2

load("st_nodes_neot_M2_YK2021.Robj")

## Chromosome number

ch.node.mat<-t(st.nodes.neot) %*% stch.mat
nodestate<-as.vector(ch.node.mat %*% 1:70)
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[stch.vec[phy.o.neot$tip.state]*10-9]

plot(phy.o.neot,"fan",no.margin=T,edge.width=0.4,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=0.6)
nodelabels(pch=16,col=nodecol,frame="n",cex=1)
text(panel_label_par[1],panel_label_par[2],labels="C",cex=panel_label_par[5])
mtext("Eurypterygii (M2)",side=2,cex=1)

## Arm number

ar.node.mat<-t(st.nodes.neot) %*% star.mat
nodestate<-as.vector(ar.node.mat %*% 1:70)
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[star.vec[phy.o.neot$tip.state]*10-9]

plot(phy.o.neot,"fan",no.margin=T,edge.width=0.4,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=0.6)
nodelabels(pch=16,col=nodecol,frame="n",cex=1)
text(panel_label_par[1],panel_label_par[2],labels="D",cex=panel_label_par[5])

## Indexes

indexvec<-matlab.like(chrmax*2)
barplot(rep(2,length(indexvec)),
        axes = F, 
        border=NA,
        space = 0,
        col = indexvec)
axis(1,lwd=1,cex.axis=1,mgp=c(0,0.5,0))

dev.off()

##### Otophysi tree ######

## Parameters

postscript("Fig_tree_with_ancestral_states_Otophysi_XXXXXX.eps", height = 7.4, width = 7.4,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)
ew<-0.2
panel_label_par<-c(-150,140,-127,116,3)

par(oma=c(1.5,2,2,0))
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T),
       widths=c(3.65,3.65),heights=c(3.65,3.65,0.2))
colvec<-matlab.like((chrmax*2-1)*10+1)

## M0

load("st_nodes_otop_M0_YK2021.Robj")

## Chromosome number

ch.node.mat<-t(st.nodes.otop) %*% stch.mat
nodestate<-as.vector(ch.node.mat %*% 1:70)
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[stch.vec[phy.o.otop$tip.state]*10-9]

plot(phy.o.otop,"fan",no.margin=T,edge.width=0.4,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=0.6)
nodelabels(pch=16,col=nodecol,frame="n",cex=1)
text(panel_label_par[3],panel_label_par[4],labels="A",cex=panel_label_par[5])
mtext("Chromosome number evolution",side=3,cex=1)
mtext("Otophysi (M0)",side=2,cex=1)

## Arm number

ar.node.mat<-t(st.nodes.otop) %*% star.mat
nodestate<-as.vector(ar.node.mat %*% 1:70)
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[star.vec[phy.o.otop$tip.state]*10-9]

plot(phy.o.otop,"fan",no.margin=T,edge.width=0.4,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=0.6)
nodelabels(pch=16,col=nodecol,frame="n",cex=1)
text(panel_label_par[3],panel_label_par[4],labels="B",cex=panel_label_par[5])
mtext("Arm number evolution",side=3,cex=1)

## M2

load("st_nodes_otop_M2_YK2021.Robj")

## Chromosome number

ch.node.mat<-t(st.nodes.otop) %*% stch.mat
nodestate<-as.vector(ch.node.mat %*% 1:70)
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[stch.vec[phy.o.otop$tip.state]*10-9]

plot(phy.o.otop,"fan",no.margin=T,edge.width=0.4,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=0.6)
nodelabels(pch=16,col=nodecol,frame="n",cex=1)
text(panel_label_par[3],panel_label_par[4],labels="C",cex=panel_label_par[5])
mtext("Otophysi (M2)",side=2,cex=1)

## Arm number

ar.node.mat<-t(st.nodes.otop) %*% star.mat
nodestate<-as.vector(ar.node.mat %*% 1:70)
nodecolindex<-round(nodestate,1)*10-9
nodecol<-colvec[nodecolindex]
tipcol<-colvec[star.vec[phy.o.otop$tip.state]*10-9]

plot(phy.o.otop,"fan",no.margin=T,edge.width=0.4,show.tip.label=F) #cex,label size on the tip; edge.width, branch width 
tiplabels(pch=16,col=tipcol,frame="n",cex=0.6)
nodelabels(pch=16,col=nodecol,frame="n",cex=1)
text(panel_label_par[3],panel_label_par[4],labels="D",cex=panel_label_par[5])

#Index

indexvec<-matlab.like(chrmax*2)
barplot(rep(2,length(indexvec)),
        axes = F, 
        border=NA,
        space = 0,
        col = indexvec)
axis(1,lwd=1,cex.axis=1,mgp=c(0,0.5,0))

dev.off()
