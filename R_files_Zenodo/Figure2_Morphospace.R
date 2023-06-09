# Figure 2
# interspecific morphospace

load("ruminant_data.Rdata")

# first, interspecific allometry... to measure intraspecific PC1s against
# PGLS - axis of evolutionary allometry of cranial shape in ruminants
R.size.pgls <- procD.pgls(f1 = coords~csize, phy = tree, data = R.df)
summary(R.size.pgls)

R.size.pgls.sf <- procD.pgls(f1 = coords~csize*subfamily, phy = tree, data = R.df)
summary(R.size.pgls.sf)



as.proc <- gpagen(ruminants)
allspecs <- as.proc$coords
allspecs.size <- as.proc$Csize

allspecs.scores <- matrix(nrow=dim(allspecs)[3],ncol=dim(R.pca$rotation)[2])
for(i in 1:dim(allspecs)[3]){
  x <- rotonto(mshape(R),allspecs[,,i],scale=F)$yrot - mshape(R)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  allspecs.scores[i,] <- a.preds[1,] %*% R.pca$rotation
}

family.cols <- c()
family.pch <- c()
 for(i in 1:NROW(allspecs.scores)){
   name <- name.id[which(name.id[,1]==si.FV[i]),3]#specimenInfo$species.FV[i]),3]
   nf <- full.inventory1[which(rownames(full.inventory1)==name),2]
   if(nf=="Cervidae"){family.cols[[i]] <- addTrans(family.colors[1],255*0.15) ; family.pch[[i]] <- 21}else{
   if(nf=="Bovidae"){family.cols[[i]] <- addTrans(family.colors[2],255*0.15); family.pch[[i]] <- 22}else{
   if(nf=="Tragulidae"){family.cols[[i]] <- addTrans(family.colors[3],255*0.15); family.pch[[i]] <- 24}else{
   if(nf=="Moschidae"){family.cols[[i]] <- addTrans(family.colors[4],255*0.15); family.pch[[i]] <- 25}}}}}
 family.cols <- unlist(family.cols)
  family.pch <- unlist(family.pch)

 cervid.scores <- allspecs.scores[which(family.cols == addTrans(family.colors[1],255*0.15)),]
 cervid.chull <- chull(cervid.scores[,1:2])
 cervid.chull <- c(cervid.chull,cervid.chull[1])
 bovid.scores <- allspecs.scores[which(family.cols == addTrans(family.colors[2],255*0.15)),]
 bovid.chull <- chull(bovid.scores[,1:2])
 bovid.chull <- c(bovid.chull,bovid.chull[1])
 moschid.scores <- allspecs.scores[which(family.cols == addTrans(family.colors[4],255*0.15)),]
 moschid.chull <- chull(moschid.scores)
 moschid.chull <- c(moschid.chull,moschid.chull[1])
 tragulid.scores <- allspecs.scores[which(family.cols == addTrans(family.colors[3],255*0.15)),]
 tragulid.chull <- chull(tragulid.scores[,1:2])
 tragulid.chull <- c(tragulid.chull,tragulid.chull[1])


# R.bm <- mvBM(R.df$tree, data= R.pca$x[,1:5])

 R.d <- distRoot(R.df$tree)[1]



# Interspecific morphospace - PC1 & 2

# with densities
par(fig=c(0,.9,0,.9),mar=c(4,4,1,1))
# ADD ALLOMETRY VECTOR
   plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1.2,
       xlab = paste("PC1: ",toString(round(R.pca$d[1]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("PC2: ",toString(round(R.pca$d[2]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""))
grid(lty=1,col = "lightgray")

points(allspecs.scores[,1],allspecs.scores[,2],cex = rescale(log(allspecs.size),c(0.5,2.5)),col=family.cols,bg=family.cols,pch=family.pch)
  lines(cervid.scores[cervid.chull,1:2],col=family.colors[1],lwd=3)
  lines(bovid.scores[bovid.chull,1:2],col=family.colors[2],lwd=3)
  lines(moschid.scores[moschid.chull,1:2],col=family.colors[4],lwd=3)
  lines(tragulid.scores[tragulid.chull,1:2],col=family.colors[3],lwd=3)

  mixtools::ellipse(mu = as.numeric(R.bm$theta)[1:2], sigma = R.bm$sigma[1:2,1:2]*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
  mixtools::ellipse(mu = as.numeric(R.bm$theta)[1:2], sigma = R.bm$sigma[1:2,1:2]*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)

  addtree2(R.pca$x[,1:2],tree = R.df$tree)
  points(R.pca$x[,1],R.pca$x[,2],pch=pts.pch,bg = pts.col,lwd = 3, cex = rescale(log(R.size),c(0.5,2.25)))
  points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 2, lwd = 2)
  legend("topleft", c("Tragulidae","Cervidae","Moschidae","Bovidae","Common Ancestor"), col = 'black',pt.lwd = 2, pt.cex = 2,pt.bg = c(family.colors[c(3,1,4,2)],"goldenrod1"), pch = c(24,21,25,22,23),title.adj = -1, cex = 1.25)


# adding the density plots:
w <- table(si.FV)

# make vector of the relative weights of each specimen, based on the amount of specimens representing the species in this dataset
wts <- 1:length(si.FV)
names(wts) <- si.FV
for(i in 1:length(si.FV)){
  wts[i] <- 1/w[which(names(w)==si.FV[i])]
}

all.fams <- c()
for(i in 1:length(family.cols)){
  a <- name.id[which(name.id[,1] == si.FV[i]),3]
  all.fams[[i]] <- R.df$family[which(R.df.n == a)]
}
all.fams <- unlist(all.fams)


library(sm)

uaf <- unique(all.fams) ; uaf <- uaf[c(2,3,4,1)]

wts2 <- wts
wts2 <- (wts2/sum(wts2))


par(fig=c(0,0.9,0.78,1),new=T)
plot(NA,xlim=c(-0.2,0.2),ylim=c(0,6),xlab=NA,ylab=NA,
     frame = F, axes = F)
#axis(4,labels=FALSE)
for(i in 1:4){
famid <- which(all.fams == uaf[i])
wtsf <- wts2[famid]
dens <- density(allspecs.scores[famid,1],weights=wtsf/sum(wtsf))
dens$y <- dens$y * sum(wtsf)
polygon(dens,col = addTrans(family.colors[i],255*0.1), border = family.colors[i], lwd =2)
}


par(fig=c(0.785,1,0,0.9),new=T)
plot(NA,ylim=c(-0.2,0.2),xlim=c(0,6),xlab=NA,ylab=NA,
     frame = F, axes = F)
#axis(3, labels = FALSE)
for(i in 1:4){
famid <- which(all.fams == uaf[i])
wtsf <- wts2[famid]
dens <- density(allspecs.scores[famid,2],weights=wtsf/sum(wtsf))
dens$y <- dens$y * sum(wtsf)
dens2 <- dens
dens$x <- dens2$y
dens$y <- dens2$x
polygon(dens,col = addTrans(family.colors[i],255*0.1), border = family.colors[i], lwd =2)
}


par(mfrow = c(1,1), mar = c(5,5,2,2))
# THE PHYLOGENETIC DISTRIBUTION SCHEMATIC
###
library(mvtnorm)
x.points <- seq(-0.5,0.5,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
mu <- R.bm$theta[1:2]
sigma <- R.bm$sigma*R.d
sigma <- sigma[1:2,1:2]
for (i in 1:100) {
for (j in 1:100) {
  xyfoo <- c(x.points[i],y.points[j])
  z[i,j] <- mvtnorm::dmvnorm(xyfoo, mean=mu, sigma=sigma)
}
}


par(mar=c(5,5,2,2))
  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1.2,
       xlab = "PC1",ylab="PC2")
.filled.contour(x.points,y.points,z, levels = seq(0,max(z),length.out = 50), col = colorRampPalette(c("gray10","gray85"))(50))
    mixtools::ellipse(mu = mu, sigma = sigma, alpha = 0.05, col = 'Gray', lwd = 4)
   mixtools::ellipse(mu = mu, sigma = sigma, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)
   addtree2(R.pca$x[,1:2], tree = R.df$tree, col = 'white')
  points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 2, lwd = 2)
shape::colorlegend(col= colorRampPalette(c("gray10","gray85"))(50), zlim=c(0,max(z)), zlevels = NULL,
            posx = c(0.81,0.91),posy=c(0.06,0.26), main = "expected density\\nof species", main.col = "white", main.cex = 1.2)



# SHAPE VARIATION ALONG THE MORPHOSPACE
####
density.2d <- kde2d.weighted(allspecs.scores[,1], allspecs.scores[,2], w = wts,
                             lims = c(range(allspecs.scores[,1]),range(allspecs.scores[,2])))

par(mar=c(5,5,2,2))
  plot(NA,xlim=c(-0.2,0.2),ylim=c(-.2,.21),asp=T,cex.lab=1.2,
       xlab = "PC1",ylab="PC2")
grid(lty=1)
#.filled.contour(density.2d$x,density.2d$y,density.2d$z,
#                col = colorRampPalette(c("gray95","gray5"))(24),
#                levels=seq(1,max(density.2d$z),length.out=24))
points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 2, lwd = 2)
    mixtools::ellipse(mu = mu, sigma = sigma, alpha = 0.05, col = 'Gray', lwd = 4)
   mixtools::ellipse(mu = mu, sigma = sigma, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)

xseg <- c(-0.15,0,0.15)   ; yseg <- c(-0.15,0,0.15)
scores <- expand.grid(xseg,yseg)
points(scores, pch = 16, lwd = 2, cex = 1)


shapes <- matrix(NA, ncol = dim(R)[1]*dim(R)[2], nrow = nrow(scores))
PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(R)))+R.pca$rotation[,1]*scores[i,1]+R.pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  if(i==5){SHADE<-F}else{SHADE<-T}
  scrs <- as.matrix(scores[i,])
  up <- c(scrs[1],scrs[2]+0.0375)
  down <- c(scrs[1]-0.0071,scrs[2]-0.025)
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(R)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.165, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = up,
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)
  tlm <- shapes1[dor.drop,-2,i] ; tlm[,2] <- -tlm[,2]
  rlm <- mshape(R)[dor.drop,-2] ; rlm[,2] <- -rlm[,2]
  tps(target.lm = tlm,
      reference.lm = rlm,
      n.grid.col = 24,
      scale.lm = 0.16, mag = 1,
      shade = SHADE, shade.trans = 0.7,
      add = T, at = down,
      links = links.dor,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)
}

text(c(-0.025,-0.025), y=c(0.069,-0.004),labels=c("lateral","dorsal"))


# DIET SURFACE
# BROWSING-GRAZING LANDSCAPE
####
plot(NA, xlim = range(R.pca$x[,1]), ylim = range(R.pca$x[,2]), asp = 1,
     xlab = paste("PC1: ",toString(round(R.pca$d[1]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""),
       ylab = paste("PC2: ",toString(round(R.pca$d[2]/sum(R.pca$d),digits=3)*100),"% of shape variation",sep=""),
     main = "Percent Grass in Diet Data")
grid(lty=1,col="lightgray")
points(R.pca$x[,1:2],cex = rescale.numeric(R.size,c(1.5,3)))
points(R.pca$x[R2.df.n,1:2],cex = rescale.numeric(R.size,c(1.5,3))[R2.df.n],
       col = make.color.vec(100-R2.df$percentgrass,11,brewer_pal(palette = "BrBG")), pch = 19)
text(x = c(-0.2,-0.2),
     y = c(0.08, 0.18),
     labels = c("0%", "100%"))
shape::colorlegend(col = brewer_pal(palette = "BrBG",direction=-1)(11), zlim = c(0,100), zlevels = seq(0,100,length.out = 11),
            posx = c(0.1,0.2), posy = c(0.7,0.9), main = "", lab.col = "transparent")


pgrass.surf <- surf.ls(2,R.pca$x[R2.df.n,1],R.pca$x[R2.df.n,2],100-R2.df$percentgrass)
pgrass.trmat <- trmat(pgrass.surf,-.3,.3,-.3,.3,n = 24)
anova(pgrass.surf)

#plot(NA,xlab = "PC1", ylab = "PC2", xlim = range(R.pca$x[,1]), ylim = range(R.pca$x[,2]), asp = T, cex.lab = 1.5)
plot(NA,xlab = "", ylab = "", xlim = range(R.pca$x[,1]), ylim = range(R.pca$x[,2]), asp = 1, axes = F)
grid(lty=1) ;
.filled.contour(x=pgrass.trmat$x,y=pgrass.trmat$y,z = pgrass.trmat$z,col = brewer_pal(palette = "BrBG")(11),levels = seq(0,100,length.out = 11))
box()
#contour(x=pgrass.trmat$x,y=pgrass.trmat$y,z = pgrass.trmat$z,col = brewer_pal(palette = "BrBG")(11),levels = seq(100,0,length.out = 11),lwd = 5, labels = NULL, add = T)
#mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.05, col = 'Gray', lwd = 4)
#mixtools::ellipse(mu = as.numeric(R.bm$theta), sigma = R.bm$sigma*R.d, alpha = 0.01, col = 'Gray', lwd = 4, lty = 3)
 addtree2(R.pca$x[,1:2],tree = R.df$tree, col = "black", lwd = 1)
points(rbind(R.pca$anc.x[1,1:2],R.pca$anc.x[1,1:2]), pch = 23, bg = 'goldenrod1', cex = 3, lwd = 2)
#points(R.pca$x[R2.df.n,1:2],cex = rescale.numeric(R.size,c(0.75,2.25))[R2.df.n],
 #      col = make.color.vec(100-R2.df$percentgrass,11,brewer_pal(palette = "BrBG")), pch = 19)
shape::colorlegend(col = brewer_pal(palette = "BrBG",direction=-1)(11), zlim = c(0,100), zlevels = seq(0,100,length.out = 11),
            posx = c(0.1,0.3), posy = c(0.7,0.9), main = "", lab.col = "transparent")
