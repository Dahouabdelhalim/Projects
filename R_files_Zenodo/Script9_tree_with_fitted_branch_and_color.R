###### Library ######

library(diversitree)
library(extrafont)

###### Font preparation ######

subset(fonttable(), FamilyName == "Arial") ##Check not more than one "Arial" in your system
loadfonts(device = "postscript")

###### Data preparation ######

load("Dataset/phy_fitted_bl_neot_YK2021.Robj")
load("Dataset/edge_col_neot_YK2021.Robj")
phy_neot<-phy_fitted_bl
ecol_neot<-edge_col
tlen_neot<-max(node.depth.edgelength(phy_neot))
load("Dataset/phy_fitted_bl_otop_YK2021.Robj")
load("Dataset/edge_col_otop_YK2021.Robj")
phy_otop<-phy_fitted_bl
ecol_otop<-edge_col
tlen_otop<-max(node.depth.edgelength(phy_otop))

######  Phylogenetic trees with fitted branch length and coloring for branch classification (Fig 4) ######
postscript("Fig_trees_with_fbl_and_colors_XXXXXX.eps", height = 7.5, width = 7.5,
           family = "Arial", paper = "special", onefile = FALSE,
           horizontal = FALSE)

## Common parameters
nratio<-tlen_neot/(tlen_neot+tlen_otop)
layout(matrix(c(1,2),1,2,byrow=T),
      widths=c(7.5*nratio,7.5*(1-nratio)),heights=c(7.3))
par(oma=c(0,0,0.5,0))
tcol<-"white"

## Eurypterygii

top10fitbrind<-order(phy_neot$edge.length,decreasing=T)[1:10]

# Plot phylogeny of Eurypterygii
plot(phy_neot,cex=0.1,no.margin=T,edge.width=1,label.offset=1,edge.color=ecol_neot,show.tip.label=F)

# Labels of the top 10 edges
edgelabels(pch=16,edge=top10fitbrind[1],col=ecol_neot[top10fitbrind[1]],frame="n",cex=2,adj=c(0,4))
edgelabels(pch=16,edge=top10fitbrind[2],col=ecol_neot[top10fitbrind[2]],frame="n",cex=2,adj=c(1300,4))
edgelabels(pch=16,edge=top10fitbrind[3],col=ecol_neot[top10fitbrind[3]],frame="n",cex=2,adj=c(0,-4))
#phy_neot$edge[top10fitbrind[4],2]
segments(2500,569,2700,549,lwd=1,col=ecol_neot[top10fitbrind[4]])
edgelabels(pch=16,edge=top10fitbrind[4],col=ecol_neot[top10fitbrind[4]],frame="n",cex=2,adj=c(970,-20))
#phy_neot$edge[top10fitbrind[5],2]
segments(800,290,1000,310,lwd=1,col=ecol_neot[top10fitbrind[5]])
edgelabels(pch=16,edge=top10fitbrind[5],col=ecol_neot[top10fitbrind[5]],frame="n",cex=2,adj=c(140,20))
edgelabels(pch=16,edge=top10fitbrind[6],col=ecol_neot[top10fitbrind[6]],frame="n",cex=2,adj=c(500,0))
edgelabels(pch=16,edge=top10fitbrind[7],col=ecol_neot[top10fitbrind[7]],frame="n",cex=2)
edgelabels(pch=16,edge=top10fitbrind[8],col=ecol_neot[top10fitbrind[8]],frame="n",cex=2,adj=c(0,-4))
edgelabels(pch=16,edge=top10fitbrind[9],col=ecol_neot[top10fitbrind[9]],frame="n",cex=2)
edgelabels(text=1,edge=top10fitbrind[1],col=tcol,frame="n",cex=0.7,adj=c(0.49,-0.02))
edgelabels(text=2,edge=top10fitbrind[2],col=tcol,frame="n",cex=0.7,adj=c(-21.05,-0.02))
edgelabels(text=3,edge=top10fitbrind[3],col=tcol,frame="n",cex=0.7,adj=c(0.49,0.9))
edgelabels(text=4,edge=top10fitbrind[4],col=tcol,frame="n",cex=0.7,adj=c(-15.6,2.6))
edgelabels(text=5,edge=top10fitbrind[5],col=tcol,frame="n",cex=0.7,adj=c(-1.78,-1.72))
edgelabels(text=6,edge=top10fitbrind[6],col=tcol,frame="n",cex=0.7,adj=c(-7.83,0.35))
edgelabels(text=7,edge=top10fitbrind[7],col=tcol,frame="n",cex=0.7,adj=c(0.30,0.4))
edgelabels(text=8,edge=top10fitbrind[8],col=tcol,frame="n",cex=0.7,adj=c(0.49,0.9))
edgelabels(text=9,edge=top10fitbrind[9],col=tcol,frame="n",cex=0.7,adj=c(0.48,0.4))
mtext("A",side=3,cex=2,adj=0.02,padj=1)
mtext("Eurypterygii",side=3,cex=1.5,adj=0.2,padj=1.3)

## Otophysi

top10fitbrind<-order(phy_otop$edge.length,decreasing=T)[1:10]

# Plot phylogeny of Otophysi
plot(phy_otop,cex=0.1,no.margin=T,edge.width=1,label.offset=1,edge.color=ecol_otop,show.tip.label=F)

# Labels of the top 10 edges
edgelabels(pch=16,edge=top10fitbrind[1],col=ecol_otop[top10fitbrind[1]],frame="n",cex=2)
edgelabels(pch=16,edge=top10fitbrind[2],col=ecol_otop[top10fitbrind[2]],frame="n",cex=2)
phy_neot$edge[top10fitbrind[2],2]
segments(1400,110,1600,98,lwd=1,col=ecol_otop[top10fitbrind[3]])
edgelabels(pch=16,edge=top10fitbrind[3],col=ecol_otop[top10fitbrind[3]],frame="n",cex=2,adj=c(390,-12))
phy_neot$edge[top10fitbrind[4],2]
segments(450,116.5,650,128.5,lwd=1,col=ecol_otop[top10fitbrind[4]])
edgelabels(pch=16,edge=top10fitbrind[4],col=ecol_otop[top10fitbrind[4]],frame="n",cex=2,adj=c(0,12))
segments(1500,113.5,1700,125.5,lwd=1,col=ecol_otop[top10fitbrind[4]])
edgelabels(pch=16,edge=top10fitbrind[5],col=ecol_otop[top10fitbrind[5]],frame="n",cex=2,adj=c(240,12))
edgelabels(pch=16,edge=top10fitbrind[6],col=ecol_otop[top10fitbrind[6]],frame="n",cex=2)
edgelabels(pch=16,edge=top10fitbrind[7],col=ecol_otop[top10fitbrind[7]],frame="n",cex=2)
segments(700,107,900,95,lwd=1,col=ecol_otop[top10fitbrind[8]])
edgelabels(pch=16,edge=top10fitbrind[8],col=ecol_otop[top10fitbrind[8]],frame="n",cex=2,adj=c(320,-12))
segments(700,112,900,124,lwd=1,col=ecol_otop[top10fitbrind[9]])
edgelabels(pch=16,edge=top10fitbrind[9],col=ecol_otop[top10fitbrind[9]],frame="n",cex=2,adj=c(320,12))
edgelabels(pch=16,edge=top10fitbrind[10],col=ecol_otop[top10fitbrind[10]],frame="n",cex=2)
edgelabels(text=1,edge=top10fitbrind[1],col=tcol,frame="n",cex=0.7,adj=c(0.38,0.4))
edgelabels(text=2,edge=top10fitbrind[2],col=tcol,frame="n",cex=0.7,adj=c(0.38,0.4))
edgelabels(text=3,edge=top10fitbrind[3],col=tcol,frame="n",cex=0.7,adj=c(-6,2.6))
edgelabels(text=4,edge=top10fitbrind[4],col=tcol,frame="n",cex=0.7,adj=c(0.38,-1.6))
edgelabels(text=5,edge=top10fitbrind[5],col=tcol,frame="n",cex=0.7,adj=c(-3.6,-1.6))
edgelabels(text=6,edge=top10fitbrind[6],col=tcol,frame="n",cex=0.7,adj=c(0.38,0.4))
edgelabels(text=7,edge=top10fitbrind[7],col=tcol,frame="n",cex=0.7,adj=c(0.38,0.4))
edgelabels(text=8,edge=top10fitbrind[8],col=tcol,frame="n",cex=0.7,adj=c(-4.9,2.6))
edgelabels(text=9,edge=top10fitbrind[9],col=tcol,frame="n",cex=0.7,adj=c(-4.85,-1.6))
edgelabels(text=10,edge=top10fitbrind[10],col=tcol,frame="n",cex=0.7,adj=c(0.47,0.4))

# Scale bar
segments(1600,490,1800,490,lwd=3)
text(1700,496,labels="200 my",cex=0.7)

mtext("B",cex=2,adj=0.02,padj=1)
mtext("Otophysi",cex=1.55,adj=0.19,padj=1.3)

dev.off()
