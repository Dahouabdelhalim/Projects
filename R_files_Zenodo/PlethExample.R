#Empirical Example: Plethodon body proportions

devtools::install_github('geomorphR/geomorph',ref="Develop")  #new version of method until manuscript acceptance
library(ape)
library(geomorph)
library(geiger)

#Tree and data
tree.best<-read.nexus("Consensus of 1000 salamander trees.nex") #Maximum Credible Tree
plot(tree.best)
plethdata<-read.csv("meandata-CinGlutOnly.csv",row.names=1, header=TRUE)
plethtree<-drop.tip(tree.best,setdiff(tree.best$tip.label,rownames(plethdata)))
plot(plethtree, edge.width = 3)
axisPhylo(1)
 
#isometric size adjustment to get to relative proportions
group<-as.factor(plethdata[,1]); names(group)<-rownames(plethdata)
size<-as.matrix(plethdata[,2]); rownames(size)<-rownames(plethdata)
  body<-as.matrix(plethdata[,-(1:3)])
Y<-apply(body,2,function(x) x/size); rownames(Y)<-rownames(body)
  gdf<-geomorph.data.frame(Y=Y, BW = Y[,3],
                           size=size,group=group)

##Group Aggregation
C<-vcv.phylo(plethtree)
X<-model.matrix(~group); rownames(X) <- names(group)

two.b.pls(C, X,iter=9999)  #Near perfect correlation: PLS-r = 0.99, P = 0.0001
  
##Analyses  
procD.lm(Y~size, data=gdf)
procD.lm(Y~group, data=gdf)
procD.lm(Y~group*size, data=gdf, SS.type = "II")

# is it the correlation between X variables or the C matrix causing the different results?
# repeat these steps a few times

newsize <- sim.char(plethtree, par=1)[,,1]
gdf$size <- newsize
procD.pgls(Y~size*group, data=gdf, phy = plethtree, SS.type = "I")
procD.pgls(Y~group*size, data=gdf, phy = plethtree, SS.type = "I")


##Interpreting
res<-procD.pgls(Y~size*group, data=gdf, phy = plethtree)
 means<-array(NA,dim=c(2,ncol(res$pgls.fitted)))
 for (i in 1:ncol(res$pgls.fitted)){
   means[,i]<-tapply(res$pgls.fitted[,i],group,mean)
 }
 colnames(means)<-colnames(res$pgls.fitted)
 rownames(means)<-c("Glut","Cin")
 means   #Gluts proportionally larger, especially limbs
 
### PLOTS
plot(plethtree)
plotGMPhyloMorphoSpace(phy=plethtree,A=Y,node.labels = FALSE,tip.labels = FALSE)


###################################
#PhyloMorphoSpace By Hand (to adapt for this dataset)
N <- length(plethtree$tip.label)
x<-Y
x <- x[plethtree$tip.label, ]
gp.plot<-group[plethtree$tip.label]
sz.plot<-size[plethtree$tip.label,]
names <- row.names(x)
anc.states <- NULL
for (i in 1:ncol(x)) {
  options(warn = -1)
  tmp <- as.vector(ace(x[, i], compute.brlen(plethtree, 1), type = "continuous", 
                       method = "ML")$ace)
  anc.states <- cbind(anc.states, tmp)
}
colnames(anc.states) <- NULL
all.data <- rbind(x, anc.states)
phylo.PCA <- prcomp(all.data)$x

limits = function(x, s) {
  r = range(x)
  rc = scale(r, scale = F)
  l = mean(r) + s * rc }

# Plotting regular phylomorphospace
library(calibrate)
plot(phylo.PCA, type = "n", xlim = limits(phylo.PCA[, 1], 1.1), 
     ylim = limits(phylo.PCA[, 2], 1.1), asp = 1)
for (i in 1:nrow(plethtree$edge)) {
  lines(phylo.PCA[(plethtree$edge[i, ]), 1], phylo.PCA[(plethtree$edge[i, ]), 2],
        type = "l", pch = 21, col = "black", lwd = 2)}
N <- length(plethtree$tip.label)

points(phylo.PCA[(N + 1):nrow(phylo.PCA), ], pch = 21, bg = "white", cex= 1.00) #nodes
points(phylo.PCA[which(gp.plot=="Large"), ], pch = 22, 
       bg = "gray", cex= 2) # tips
points(phylo.PCA[which(gp.plot=="Small"), ], pch = 21, 
       bg = "gray", cex= 2) # tips
#textxy(phylo.PCA[1:N, 1], phylo.PCA[1:N, 2], rownames(phylo.PCA), cex= 0.5) # tip labels

