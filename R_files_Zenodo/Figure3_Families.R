# Figure 3
# Family-specific morphospaces


#########################


plot.contour <- FALSE
f2.CL <- 1.25

cervid.size <- R.size[which(fam=='Cervidae')]
cervid.size <- cervid.size[cervid.tree$tip.label] # making sure things are in the correct order

bovid.size <- R.size[which(fam=='Bovidae')]
bovid.size <- bovid.size[bovid.tree$tip.label] # making sure things are in the correct order


# cervid phylomorphospace

# predicting PC scores of minimum and maximum shapes in cervidae
allo.Cervidae <- procD.pgls(f1 = cervid.lm~cervid.size, phy = cervid.tree, data = NULL); anova(allo.Cervidae)
#             Df       SS        MS     Rsq      F      Z Pr(>F)
# cervid.size  1 0.016565 0.0165655 0.35361 17.506 4.2427  0.001 **
# Residuals   32 0.030281 0.0009463 0.64639
# Total       33 0.046847
#write.csv(summary(allo.Cervidae)$table,file="results/CervidAlloPGLS.csv")
allo.Cervidae1 <- procD.lm(f1 = cervid.lm~cervid.size, data = NULL); anova(allo.Cervidae1)


# without Alces, just to see
allo.Cervidae2 <- procD.pgls(f1 = cervid.lm[,,-c(31,32)]~cervid.size[-c(31,32)], phy = keep.tip(cervid.tree,cervid.tree$tip.label[-c(31,32)]), data = NULL); anova(allo.Cervidae2)


cervid.fit <- mvBM(tree = cervid.tree, data = cervid.pca$x[,1:2], model = "BM1")
cervid.dist <- as.numeric(distRoot(cervid.tree)[1])
buf <- 1.33 ; buf2 <- buf*1.5

cervid.x.points <- seq((range(cervid.pca$x[,1])*buf2)[1],(range(cervid.pca$x[,1])*buf2)[2],length.out=100)
cervid.y.points <- seq((range(cervid.pca$x[,2])*buf2)[1],(range(cervid.pca$x[,2])*buf2)[2],length.out=100)
cervid.z <- matrix(0,nrow=100,ncol=100)
cervid.mu <-  as.numeric(cervid.fit$theta)

for (i in 1:100) {
for (j in 1:100) {
cervid.z[i,j] <- mixtools::dmvnorm(c(cervid.x.points[i],cervid.y.points[j]), mu=cervid.mu,sigma=cervid.fit$sigma * cervid.dist)
}
}
#contour(cervid.x.points,cervid.y.points,cervid.z )
#plot(expand.grid(cervid.x.points,cervid.y.points),cex = rescale.numeric(cervid.z ), pch=19)
NLEVEL <- 24
cervid.levels <- seq(0,max(range(cervid.z)), length.out = NLEVEL)


bovid.mean.scores <- matrix(nrow=dim(bovid.lm)[3],ncol=4)
for(i in 1:dim(bovid.lm)[3]){
  x <- rotonto(mshape(cervid.lm),bovid.lm[,,i],scale=FALSE)$yrot - mshape(cervid.lm)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  bovid.mean.scores[i,] <- a.preds[1,] %*% (cervid.pca$rotation)[,1:4]
}
rownames(bovid.mean.scores) <- dimnames(bovid.lm)[[3]]


cervidPC1.d <- round(cervid.pca$d[1]/sum(cervid.pca$d),3)*100
cervidPC2.d <- round(cervid.pca$d[2]/sum(cervid.pca$d),3)*100

cervid.allo.plot <- plot(allo.Cervidae, type = "regression", predictor = cervid.size,reg.type ="RegScore", pch = 21, lwd = 2,
                         bg = family.colors[1], cex = rescale.numeric(cervid.size,to = c(1,3.5)), xlab = "log-transformed centroid size")
cervid.preds <- shape.predictor(allo.Cervidae$GM$pgls.fitted, x= cervid.allo.plot$RegScore, Intercept = FALSE,
                        predmin = max(cervid.allo.plot$RegScore),
                        predmax = min(cervid.allo.plot$RegScore))
cervid.allo.scores2 <- matrix(data=NA,nrow=2,ncol=27)
x <- bindArr(rotonto(mshape(cervid.lm),cervid.preds$predmin,scale=TRUE)$yrot,rotonto(mshape(cervid.lm),cervid.preds$predmax,scale=TRUE)$yrot,along=3)
x <- sweep(x, 1:2, mshape(cervid.lm))
c.preds <- two.d.array(x)
cervid.allo.scores2 <- c.preds %*% (cervid.pca$rotation)[,1:27]


par(mar=c(5,5,3,2))
# plotting the phylomorphospace
plot(NA,xlim = range(cervid.pca$x[,1])*buf, ylim =range(cervid.pca$x[,2])*buf, asp = T,
     xlab = paste("PC1: ",cervidPC1.d,"% of variation",sep=""),
     ylab = paste("PC2: ",cervidPC2.d,"% of variation",sep=""), cex.lab = 1.25)
title("Cervidae",adj=0,line=0.1,cex.main=4,col.main=family.colors[5],font=2)
grid(lty=1,col=addTrans('lightgray',255*0.5))
plot.contour = FALSE
if(plot.contour){.filled.contour(cervid.x.points,cervid.y.points,cervid.z,
          levels = cervid.levels, col = addTrans(colorRampPalette(c('white',family.colors[3],family.colors[1],family.colors[5]))(NLEVEL), 255*0.33))}
cervid.ellipse95 <- mixtools::ellipse(mu = as.numeric(cervid.fit$theta), sigma = cervid.fit$sigma*cervid.dist,
        alpha = 0.05, col = family.colors[1], lwd = 3)
cervid.ellipse99 <- mixtools::ellipse(mu = as.numeric(cervid.fit$theta), sigma = cervid.fit$sigma*cervid.dist,
        alpha = 0.01, col = family.colors[1], lwd = 3, lty = 2)
addtree2(bovid.mean.scores[,1:2],bovid.tree,"darkgray",1, lwd = 1)
addtree2(cervid.pca$x[,1:2],cervid.tree,"black",1, lwd = 1)
points(x=bovid.mean.scores[,1:2],pch= 15, col = addTrans(family.colors[2],255*0.5), cex = rescale.numeric(bovid.size,c(1,4)))
points(x = cervid.pca$x[,1:2], pch=21,bg=family.colors[1],cex = rescale.numeric(cervid.size,c(1,4)))
arrows(cervid.allo.scores2[1,1],cervid.allo.scores2[1,2],cervid.allo.scores2[2,1],cervid.allo.scores2[2,2],
        lwd = 7)
arrows(cervid.allo.scores2[1,1],cervid.allo.scores2[1,2],cervid.allo.scores2[2,1],cervid.allo.scores2[2,2],
       col="white", lwd = 3)
points(cervid.fit$theta, bg = "goldenrod1", pch = 23, cex =2.5, lwd = 3)
text(x = c(0.16,0.2025), y = c(0.05,0.05), labels = c('95%','99%'), cex = , col = family.colors[1])

cervid.pts95.1 <- matrix( c(0,.14,.165,-.105,0.0965,.09,0,-0.111), nrow = 4, ncol=2)
    points(x = cervid.pts95.1,lwd=3,pch=21,bg='white',cex=2)
cervid.pts95 <- cervid.pts95.1
cervid.pts95[1,1] <-  -0.05; cervid.pts95[2,1] <- 0.21
cervid.pts95[3,] <- c(.18,-.025) ; cervid.pts95[4,] <- c(-0.05,-0.0875)
scores <- cervid.pts95.1
shapes <- matrix(NA, ncol = dim(cervid.lm)[1]*dim(cervid.lm)[2], nrow = nrow(scores))
#PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(cervid.lm)))+cervid.pca$rotation[,1]*scores[i,1]+cervid.pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)

# to add the TPS figs
for(i in 1:NROW(scores)){
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(cervid.lm)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = .11, mag = 1,
      shade = T, shade.trans = 0.7,
      add = T, at = as.matrix(cervid.pts95)[i,],
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)
}
tps(target.lm = R[lat.drop,-3,"Alces_alces"],
      reference.lm = mshape(cervid.lm)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = .11, mag = 1,
      shade = T, shade.trans = 0.7,
      add = T, at = c(.2,-.12),
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)




###




# Bovidae phylomorphospace
allo.Bovidae <- procD.pgls(f1 = bovid.lm~bovid.size, phy = bovid.tree, data = NULL); anova(allo.Bovidae)$table
#write.csv(anova(allo.Bovidae)$table,file="results/BovidAlloPGLS.csv")
allo.Bovidae1 <- procD.lm(f1 = bovid.lm~bovid.size, data = NULL); anova(allo.Bovidae1)$table

bovid.fit <- mvBM(tree = bovid.tree, data = bovid.pca$x[,1:2], model = "BM1")
bovid.dist <- as.numeric(distRoot(bovid.tree)[1])
buf <- 1.2 ; buf2 <- buf*1.5

bovid.x.points <- seq((range(bovid.pca$x[,1])*buf2)[1],(range(bovid.pca$x[,1])*buf2)[2],length.out=100)
bovid.y.points <- seq((range(bovid.pca$x[,2])*buf2)[1],(range(bovid.pca$x[,2])*buf2)[2],length.out=100)
bovid.z <- matrix(0,nrow=100,ncol=100)
bovid.mu <-  as.numeric(bovid.fit$theta)

for (i in 1:100) {
for (j in 1:100) {
bovid.z[i,j] <- mixtools::dmvnorm(c(bovid.x.points[i],bovid.y.points[j]), mu=bovid.mu,sigma=bovid.fit$sigma * bovid.dist)
}
}
#contour(bovid.x.points,bovid.y.points,bovid.z )
#plot(expand.grid(bovid.x.points,bovid.y.points),cex = rescale.numeric(bovid.z ), pch=19)

bovid.levels <- seq(0,max(range(bovid.z)), length.out = NLEVEL)










cervid.mean.scores <- matrix(nrow=dim(cervid.lm)[3],ncol=4)
for(i in 1:dim(cervid.lm)[3]){
  x <- rotonto(mshape(bovid.lm),cervid.lm[,,i],scale=FALSE)$yrot - mshape(bovid.lm)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  cervid.mean.scores[i,] <- a.preds[1,] %*% (bovid.pca$rotation)[,1:4]
}
rownames(cervid.mean.scores) <- dimnames(cervid.lm)[[3]]

bovidPC1.d <- round(bovid.pca$d[1]/sum(bovid.pca$d),3)*100
bovidPC2.d <- round(bovid.pca$d[2]/sum(bovid.pca$d),3)*100

bovid.allo.plot <- plot(allo.Bovidae, type = "regression", predictor = bovid.size ,reg.type ="RegScore", xlab = "log-transformed centroid size", pch =21, lwd = 2, bg = family.colors[2], cex = rescale.numeric(bovid.size,to=c(1,4)))
bovid.preds <- shape.predictor(allo.Bovidae$GM$pgls.fitted, x= bovid.allo.plot$RegScore, Intercept = FALSE,
                        predmin = min(bovid.allo.plot$RegScore),
                        predmax = max(bovid.allo.plot$RegScore))

bovid.allo.scores2 <- matrix(data=NA,nrow=2,ncol=27)
x <- bindArr(rotonto(mshape(bovid.lm),bovid.preds$predmin,scale=TRUE)$yrot,rotonto(mshape(bovid.lm),bovid.preds$predmax,scale=TRUE)$yrot,along=3)
x <- sweep(x, 1:2, mshape(bovid.lm))
c.preds <- two.d.array(x)
bovid.allo.scores2 <- c.preds %*% (bovid.pca$rotation)[,1:27]




plot(NA,xlim = range(bovid.pca$x[,1])*buf, ylim =range(bovid.pca$x[,2])*buf, asp = T,
     xlab = paste("PC1: ",bovidPC1.d,"% of variation",sep=""),
     ylab = paste("PC2: ",bovidPC2.d,"% of variation",sep=""), cex.lab = f2.CL)
grid(lty=1,col=addTrans('lightgray',255*0.5))
title("Bovidae",adj=0,line=0.1,cex.main=4,col.main=family.colors[6],font=2)
if(plot.contour){.filled.contour(bovid.x.points,bovid.y.points,bovid.z,
          levels = bovid.levels, col = addTrans(colorRampPalette(c('white',family.colors[4],family.colors[2],family.colors[6]))(NLEVEL), 255*0.33))}
bovid.ellipse95 <- mixtools::ellipse(mu = as.numeric(bovid.fit$theta), sigma = bovid.fit$sigma*bovid.dist,
        alpha = 0.05, col = family.colors[2], lwd = 3)
bovid.ellipse99 <- mixtools::ellipse(mu = as.numeric(bovid.fit$theta), sigma = bovid.fit$sigma*bovid.dist,
        alpha = 0.01, col = family.colors[2], lwd = 3, lty = 2)
addtree2(cervid.mean.scores[,1:2],cervid.tree,"darkgray",1, lwd = 1)
addtree2(bovid.pca$x[,1:2],bovid.tree,"black",1, lwd = 1)

points(x=cervid.mean.scores[,1:2],pch= 16, col = addTrans(family.colors[1],255*0.5), cex = rescale.numeric(cervid.size,c(1,4)))
points(x = bovid.pca$x[,1:2], pch=22,bg=family.colors[2],cex = rescale.numeric(bovid.size,c(1,4)))
arrows(bovid.allo.scores2[1,1],bovid.allo.scores2[1,2],bovid.allo.scores2[2,1],bovid.allo.scores2[2,2],
        lwd = 7)
arrows(bovid.allo.scores2[1,1],bovid.allo.scores2[1,2],bovid.allo.scores2[2,1],bovid.allo.scores2[2,2],
       col= "white", lwd = 3)
points(bovid.fit$theta, bg = "goldenrod1", pch = 23, cex =2.5, lwd = 3)

text(x = c(0.11,0.15), y = c(0.05,0.05), labels = c('95%','99%'), cex = f2.CL, col = family.colors[2])

bovid.pts95.1 <- rbind(matrix( c(0,0,0.1185,-.0975), nrow = 2, ncol = 2), matrix( c(0.137,-0.148,0,0), nrow = 2, ncol = 2))
  points(x = bovid.pts95.1,lwd=3,pch=21,bg='white',cex=2)


bovid.pts95 <- bovid.pts95.1
bovid.pts95[1,2] <- 0.14 ; bovid.pts95[2,] <- c(0.03,-.13)
bovid.pts95[3,1] <- .2 ; bovid.pts95[4,] <- c(-0.155,0.025)
scores <- bovid.pts95.1
shapes <- matrix(NA, ncol = dim(bovid.lm)[1]*dim(bovid.lm)[2], nrow = nrow(scores))
PC <- bovid.pca$x[,1:2]
for(i in 1:nrow(shapes)){
  shapes[i,] <- c(t(mshape(bovid.lm)))+bovid.pca$rotation[,1]*scores[i,1]+bovid.pca$rotation[,2]*scores[i,2]
}
shapes1 <- arrayspecs(A = shapes, k = 3, p = 25)


for(i in 1:NROW(scores)){
  tps(target.lm = shapes1[lat.drop,-3,i],
      reference.lm = mshape(bovid.lm)[lat.drop,-3],
      n.grid.col = 24,
      scale.lm = 0.115, mag = 1,
      shade = T, shade.trans = 0.7,
      add = T, at = as.matrix(bovid.pts95)[i,],
      links = links.lat,plot.ref.links = T,
      ref.link.aes = list(lwd = 2, col = "darkgray"),
      link.aes = list(lwd = 2, col = "black"),
      grid.aes = list(col = 'gray', trans = 0.8),
      plot.lm = FALSE)
}





# ruminant populations
####

# specimen scores

# all cervid specimens
cervid.spec.n <- name.id[name.id[,3] %in% names(f.fam)[f.fam=="Cervidae"],1] ; names(cervid.spec.n) <- name.id[name.id[,3] %in% names(f.fam)[f.fam=="Cervidae"],3]

cervid.specs <- specimenInfo$species.FV %in% name.id[name.id[,3] %in% names(f.fam)[f.fam=="Cervidae"],1]
cervid.specs2 <- specimenInfo$species.FV %in% c("Odo_he_he1","Odo_he_ca2","Odo_he_co4","Odo_vi_bo6","Odo_vi_le3","Odo_vi_co7")
cervid.spec <- ruminant.lm[1:25,,cervid.specs|cervid.specs2]
cervid.size <- log(apply(cervid.spec,3,cSize))

bovid.spec.n <- name.id[name.id[,3] %in% names(f.fam)[f.fam=="Bovidae"],1] ; names(bovid.spec.n) <- name.id[name.id[,3] %in% names(f.fam)[f.fam=="Bovidae"],3]

bovid.specs <- specimenInfo$species.FV %in% name.id[name.id[,3] %in% names(f.fam)[f.fam=="Bovidae"],1]
bovid.spec <- ruminant.lm[1:25,,bovid.specs]
bovid.size <- log(apply(bovid.spec,3,cSize))


cervid.spec.n2 <- name.id[name.id[,3] %in% names(fam)[fam=="Cervidae"],1] ; names(cervid.spec.n2) <- name.id[name.id[,3] %in% names(fam)[fam=="Cervidae"],3]
bovid.spec.n2 <- name.id[name.id[,3] %in% names(fam)[fam=="Bovidae"],1] ; names(bovid.spec.n2) <- name.id[name.id[,3] %in% names(fam)[fam=="Bovidae"],3]
# loop that calculates specimen PC scores for each specimen (i) for each species (j)

ALSPSIZE <- rescale.numeric(log(allspecs.size),c(1,4))

cervidspecs.scores <- c()
cervid.SIZES <- c()
for(j in 1:length(cervid.spec.n)){
  SPECS <- allspecs[,,which(dimnames(allspecs)[[3]] == cervid.spec.n[[j]])]
  cervid.SIZES[[j]] <- ALSPSIZE[which(dimnames(allspecs)[[3]] == cervid.spec.n[[j]])]
  SPECSscores <- matrix(nrow=dim(SPECS)[3],ncol=dim(cervid.pca$rotation)[2])

for(i in 1:dim(SPECS)[3]){
  x <- rotonto(mshape(cervid.lm),SPECS[,,i],scale=FALSE, reflection = T)$yrot - mshape(cervid.lm)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  SPECSscores[i,] <- a.preds[1,] %*% cervid.pca$rotation
}
  cervidspecs.scores[[j]] <- SPECSscores
  names(cervidspecs.scores)[[j]] <- names(cervid.spec.n)[[j]]
}

bovidspecs.scores <- c()
bovid.SIZES <- c()
for(j in 1:length(bovid.spec.n)){
  SPECS <- allspecs[,,which(dimnames(allspecs)[[3]] == bovid.spec.n[[j]])]
  SPECSscores <- matrix(nrow=dim(SPECS)[3],ncol=dim(bovid.pca$rotation)[2])
  bovid.SIZES[[j]] <- ALSPSIZE[which(dimnames(allspecs)[[3]] == bovid.spec.n[[j]])]

for(i in 1:dim(SPECS)[3]){
  x <- rotonto(mshape(bovid.lm),SPECS[,,i],scale=F)$yrot - mshape(bovid.lm)
  x <- bindArr(x,x,along=3)
  a.preds <- two.d.array(x)
  SPECSscores[i,] <- a.preds[1,] %*% bovid.pca$rotation
}
  bovidspecs.scores[[j]] <- SPECSscores
  names(bovidspecs.scores)[[j]] <- names(bovid.spec.n)[[j]]
}


cervid.pop.col <- make.color.vec(c(0,na.omit(Vproj[cervid.tree$tip.label])), 100, colorRampPalette(c('white',family.colors[1])))
cervid.pop.col <- cervid.pop.col[-1]
bovid.pop.col <- make.color.vec(c(0,na.omit(Vproj[bovid.tree$tip.label])), 100, colorRampPalette(c('white',family.colors[2])))
bovid.pop.col <- bovid.pop.col[-1]

cervid.specs.size <- log(allspecs.size[which(names(allspecs.size) %in% cervid.spec.n)])
bovid.specs.size <- log(allspecs.size[which(names(allspecs.size) %in% bovid.spec.n)])


# the plots:
plot(NA,xlim = range(cervid.pca$x[,1])*buf, ylim =range(cervid.pca$x[,2])*buf, asp = T,
     xlab = paste("PC1: ",cervidPC1.d,"% of variation",sep=""),
     ylab = paste("PC2: ",cervidPC2.d,"% of variation",sep=""), cex.lab = f2.CL)
grid(lty=1,col=addTrans('lightgray',255*0.5))
addtree2(x = cervid.pca$x[,1:2],tree = cervid.tree, col = "lightgray")

for(i in 1:length(cervidspecs.scores)){
  xy <- cervidspecs.scores[[i]][,1:2]
  points(xy, pch = 19, col = addTrans(family.colors[1],255*0.1), cex=cervid.SIZES[[i]])
  mixtools::ellipse(mu = colMeans(xy), sigma = cov(xy), alpha = 0.2, col = 'black', lwd = 4, lty = 1)
  mixtools::ellipse(mu = colMeans(xy), sigma = cov(xy), alpha = 0.2, col = family.colors[1], lwd = 2, lty = 1)
}
arrows(cervid.allo.scores2[1,1],cervid.allo.scores2[1,2],cervid.allo.scores2[2,1],cervid.allo.scores2[2,2],
        lwd = 7)
arrows(cervid.allo.scores2[1,1],cervid.allo.scores2[1,2],cervid.allo.scores2[2,1],cervid.allo.scores2[2,2],
       col= "white", lwd = 3)
points(cervid.fit$theta, bg = "goldenrod1", pch = 23, cex =2.5, lwd = 3)


plot(NA,xlim = range(bovid.pca$x[,1])*buf, ylim =range(bovid.pca$x[,2])*buf, asp = T,
     xlab = paste("PC1: ",bovidPC1.d,"% of variation",sep=""),
     ylab = paste("PC2: ",bovidPC2.d,"% of variation",sep=""), cex.lab = f2.CL)
grid(lty=1,col=addTrans('lightgray',255*0.5))
addtree2(x = bovid.pca$x[,1:2],tree = bovid.tree, col = "lightgray")

for(i in 1:length(bovidspecs.scores)){
  xy <- bovidspecs.scores[[i]][,1:2]
  points(xy, pch = 15, col = addTrans(family.colors[2],255*0.1), cex = bovid.SIZES[[i]])
  mixtools::ellipse(mu = colMeans(xy), sigma = cov(xy), alpha = 0.2, col = 'black', lwd = 4, lty = 1)
  mixtools::ellipse(mu = colMeans(xy), sigma = cov(xy), alpha = 0.2, col = family.colors[2], lwd = 2, lty = 1)
}
arrows(bovid.allo.scores2[1,1],bovid.allo.scores2[1,2],bovid.allo.scores2[2,1],bovid.allo.scores2[2,2],
        lwd = 7)
arrows(bovid.allo.scores2[1,1],bovid.allo.scores2[1,2],bovid.allo.scores2[2,1],bovid.allo.scores2[2,2],
       col='white', lwd = 3)
points(bovid.fit$theta, bg = "goldenrod1", pch = 23, cex =2.5, lwd = 3)



# contours of face length and nasal retraction
cervid.grid <- pca.grid(cervid.pca, PC1.segments = 10, PC2.segments = 10)

#cervid.grid.x.points <- unlist(lapply(cervid.grid$score.grid,FUN = function(x)x[1]))
#cervid.grid.y.points <- unlist(lapply(cervid.grid$score.grid,FUN = function(x)x[2]))
cervid.z.nasalretraction <- apply(cervid.grid$shape.grid,3,nasal.retraction)
cervid.z.facelength <- apply(cervid.grid$shape.grid,3,face.length)

CNRsurf <- surf.ls(np = 2,cervid.grid$score.grid[,1], cervid.grid$score.grid[,2],cervid.z.nasalretraction)
CNRmat <- trmat(CNRsurf,range(cervid.grid$score.grid[,1])[1],range(cervid.grid$score.grid[,1])[2],
                range(cervid.grid$score.grid[,2])[1],range(cervid.grid$score.grid[,2])[2], n = 100)

# inset nasal retraction surface
plot(NA, xlim = range(cervid.grid$score.grid[,1])*0.875,
        ylim = range(cervid.grid$score.grid[,2])*0.875, xlab = "", ylab = "", axes = F, frame = T, asp = 1,
     main = "Nasal Retraction", cex.main = 2.5)
.filled.contour(x=CNRmat$x,y=CNRmat$y,z=CNRmat$z, levels = seq(min(CNRmat$z), max(CNRmat$z), length.out = 24),
                col = colorRampPalette(c('white', family.colors[5], "black"))(24)) ; box()


CFLsurf <- surf.ls(np = 2,cervid.grid$score.grid[,1], cervid.grid$score.grid[,2],cervid.z.facelength)
CFLmat <- trmat(CFLsurf,range(cervid.grid$score.grid[,1])[1],range(cervid.grid$score.grid[,1])[2],
                range(cervid.grid$score.grid[,2])[1],range(cervid.grid$score.grid[,2])[2], n = 100)

# inset face space
plot(NA, xlim = range(cervid.grid$score.grid[,1])*0.875,
        ylim = range(cervid.grid$score.grid[,2])*0.875, xlab = "", ylab = "", axes = F, frame = T, asp = 1,
     main = "Face Length", cex.main = 2.5)
.filled.contour(x=CFLmat$x,y=CFLmat$y,z=CFLmat$z, levels = seq(min(CFLmat$z), max(CFLmat$z), length.out = 24),
                col = colorRampPalette(c('white', family.colors[5], "black"))(24)) ; box()



bovid.grid <- pca.grid(bovid.pca, PC1.segments = 10, PC2.segments = 10)

#cervid.grid.x.points <- unlist(lapply(cervid.grid$score.grid,FUN = function(x)x[1]))
#cervid.grid.y.points <- unlist(lapply(cervid.grid$score.grid,FUN = function(x)x[2]))
bovid.z.nasalretraction <- apply(bovid.grid$shape.grid,3,nasal.retraction)
bovid.z.facelength <- apply(bovid.grid$shape.grid,3,face.length)

BNRsurf <- surf.ls(np = 2,bovid.grid$score.grid[,1], bovid.grid$score.grid[,2],bovid.z.nasalretraction)
BNRmat <- trmat(BNRsurf,range(bovid.grid$score.grid[,1])[1],range(bovid.grid$score.grid[,1])[2],
                range(bovid.grid$score.grid[,2])[1],range(bovid.grid$score.grid[,2])[2], n = 100)

# inset face space
plot(NA, xlim = range(bovid.grid$score.grid[,1])*0.875,
        ylim = range(bovid.grid$score.grid[,2])*0.875, xlab = "", ylab = "", axes = F, frame = T, asp = 1,
     main = "Nasal Retraction", cex.main = 2.5)
.filled.contour(x=BNRmat$x,y=BNRmat$y,z=BNRmat$z, levels = seq(min(BNRmat$z), max(BNRmat$z), length.out = 24),
                col = colorRampPalette(c('white', family.colors[6], "black"))(24)) ; box()

BFLsurf <- surf.ls(np = 2,bovid.grid$score.grid[,1], bovid.grid$score.grid[,2],bovid.z.facelength)
BFLmat <- trmat(BFLsurf,range(bovid.grid$score.grid[,1])[1],range(bovid.grid$score.grid[,1])[2],
                range(bovid.grid$score.grid[,2])[1],range(bovid.grid$score.grid[,2])[2], n = 100)

plot(NA, xlim = range(bovid.grid$score.grid[,1])*0.875,
        ylim = range(bovid.grid$score.grid[,2])*0.875, xlab = "", ylab = "", axes = F, frame = T, asp = 1,
     main = "Face Length", cex.main = 2.5)
.filled.contour(x=BFLmat$x,y=BFLmat$y,z=BFLmat$z, levels = seq(min(BFLmat$z), max(BFLmat$z), length.out = 24),
                col = colorRampPalette(c('white', family.colors[6], "black"))(24)) ; box()


#########################
