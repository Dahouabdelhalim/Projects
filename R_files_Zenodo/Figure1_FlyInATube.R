# Figure 1:
# evolutionary allometry of relative face length in Ruminantia, and how it impacts morphospace occupation


load("ruminant_data.Rdata")

colpal <- c("firebrick", "dodgerblue3", "#F2BC1B")

fl <- R.df$facelength
size <- R.df$csize

sfl.fit <- mvBM(tree = R.df$tree, data = cbind(size,fl), model = "BM1")

nsim <- 1000

corrs <- array(NA, dim = c(length(size),2,nsim))
for(i in 1:nsim){
  corrs[,,i] <- sim.corrs(tree = R.df$tree, vcv = sfl.fit$sigma, anc = sfl.fit$theta)
}

 size.anc <- fastAnc(tree = R.df$tree, x = size)[1]
 fl.anc <- fastAnc(tree = R.df$tree, x = fl)[1]

size.bm <- fastBM(tree = R.df$tree, a = size.anc, sig2 = sfl.fit$sigma[1,1], nsim = nsim)
fl.bm <- fastBM(tree = R.df$tree, a = fl.anc, sig2 = sfl.fit$sigma[2,2], nsim = nsim)

# mvSIM(tree = R.df$tree, nsim = nsim, model = "BMM", param = )

bm <- array(dim = dim(corrs))
for(i in 1:nsim){
bm[,1,i] <- size.bm[,i]
bm[,2,i] <- fl.bm[,i]
}

bm1 <- bm[,,1]
corrs1 <- corrs[,,1]
for(i in 1:(nsim-1)){
  bm1 <- rbind(bm1,bm[,,(i+1)])
  corrs1 <- rbind(corrs1,corrs[,,(i+1)])
  }
# or , use do.call() ^



# the plot
par(mar=c(5,5,4,4))
plot(NA, asp = diff(range(size))/diff(range(fl)),
          xlim = c(4,7), ylim = c(0.25,0.85),
     xlab = "Cranium Size", ylab = "Relative Face Length", axes = T,
     cex.lab = 2, cex.main = 3)
grid(lty=1,col=addTrans('lightgray',255*0.5))
polygon(ellipse(mu = colMeans(bm1), sigma = cov(bm1), alpha = 0.05, lwd = 8, newplot = FALSE, draw = FALSE), col = addTrans("lightgray",255*0.1))
ellipse(mu = colMeans(bm1), sigma = cov(bm1), alpha = 0.05, col = "gray", lwd = 6)
polygon(ellipse(mu = colMeans(corrs1), sigma = cov(corrs1), alpha = 0.05, lwd = 8, newplot = FALSE, draw = FALSE), col = addTrans(colpal[2],255*0.1))
ellipse(mu = colMeans(corrs1), sigma = cov(corrs1), alpha = 0.05, col = colpal[2], lwd = 6)
addtree2(x = cbind(size,fl), tree = R.df$tree, col = "gray", lwd = 2)
points( cbind(size,fl), bg = pts.col, pch = pts.pch, cex = rescale.numeric(size, to = c(1.25,3.4)),lwd= 2)
points(x = size.anc, y = fl.anc, bg = colpal[3], pch = 23, cex = 3,lwd= 3)
legend(x = 'topright',legend = c("Expected Distributions:","No Correlation", "Observered Correlation"), pch = c(NA,21,21),pt.bg = c('transparent',addTrans("lightgray",255*0.5),addTrans(colpal[2],255*0.5)) ,col = c('transparent','gray',colpal[2]), pt.cex = c(0,4,4), cex = 1.75,pt.lwd=3)
legend("topleft", c("Tragulidae","Cervidae","Moschidae","Bovidae","Common Ancestor"),pt.lwd = 2, pt.cex = 3.5,pt.bg = c(family.colors[c(3,1,4,2)],colpal[3]), pch = c(24,21,25,22,23), cex = 1.75)

