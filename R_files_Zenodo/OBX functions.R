#-----------------------------------------------------------------------
# Tyler, C. and Kowalewski, M., 2016, PRSB
#
# R Functions
#
# Last updated: December 6, 2016
#
# written by: M. Kowalewski and C. Tyler
#-----------------------------------------------------------------------

# Beta variance function based on Bray-Curtis distance
beta.var <- function(x) {
    SS <- sum(vegdist(x, method='bray')^2)/nrow(x)
    BDT <-ncol(x)/mean(specnumber(x)) - 1
    return(c('S'=sum(colSums(x)>0), 'N'=sum(x), 'TSumSq'=SS,'BDTotal'=BDT, 
           'Beta Variance'=round(SS/(sum(rowSums(x)>0)-1),4)))
  }

# Subsampling functions (beta variance and beta shannon)
myrar1 <- function(x, min, times) sapply(replicate(times, list(rrarefy(x,min))), beta.var)
myrar2 <- function(x, min, times) sapply(replicate(times, list(rrarefy(x,min))), H, lev='beta')

# Subsampling functions for estimating number of species retained in a subsampled dataset
find.sp <- function(x, loc.min) {
    y <- rrarefy(x,loc.min)
    return(sum(colSums(y)>0))
}

# function assembling pairwise distances, depth differences, and habitat type
bgrad <- function(data, depth, habitat) {
 n.comp <- (nrow(data)-1)*nrow(data)/2
 k <- 0
 out1<- vector(mode='numeric', length=n.comp*4)
 for (i in 1:nrow(data)) {
   for (j in 1:nrow(data)) {
     if (j>i) {
       k <- k + 1
       BC<- sum(abs(data[i,]-data[j,]))/(sum(data[i,])+sum(data[j,]))
       depthd<- abs(depth[i]-depth[j])
       out1[(4*k-3):(4*k)] <- c(i,j,(1-BC),depthd)
      }
     }
    }
   return(matrix(out1, n.comp, 4, byrow=T))
  }

# function for decay models
decay<-function(x){
  fitH<- x[which(x[,1]*x[,2]>0),] # get rid of 0 values (not allowed for log operations)
  expH<- lm(log(fitH[,2])~ fitH[,1]) # estimate exponent
  xlH<- seq(0.1,15,0.1) # set a series of points for range of depth difference values (for plotting over a given range of x values)
  ylH<- exp(expH$coefficients[2]*xlH+expH$coefficients[1]) # compute predicted y values (bray-curtis values)
  adjr2 <- summary(expH)$adj.r.squared # compute adjusted r-square
  decay<-cbind(xlH,ylH,expH$coefficients[2])
  return(decay)
}

# function to plot beta gradient figure
b.g.plot <- function(x, y, bin, my.yl=c(0,1), my.xl=c(0,15), ct=mean, mcex=0.5, 
                     mxlab='', mylab='', mave=F, bin2=10) {
   plot(x, y, type='n', ylim=my.yl, xlim=my.xl,
        xlab=mxlab, ylab=mylab, las=1)
    bins <- round(1/bin*x)/(1/bin)
    l.y <- tapply(y,bins,ct)
    l.x <- as.numeric(names(l.y))
    if (mave) {
     n.y <- length(y) - length(y)%%bin2
     x1 <- x[order(x)]; y1 <- y[order(x)]
     l.y <- apply(matrix(y1[1:n.y],n.y/bin2,bin2,byrow=T),1,ct) 
     l.x <- apply(matrix(x1[1:n.y],n.y/bin2,bin2,byrow=T),1,ct)
     } 
    d.l <- decay(cbind(x,y))
    points(x, y, pch='.', col='black', cex=mcex) 
#   points(d.l[,1], d.l[,2], type='l', lty=1, lwd=2, col='skyblue') # decay line (disabled)
    points(l.x, l.y, type='l', col='white', lwd=6)
    points(l.x, l.y, type='l', col='red', lwd=2)
    text(my.xl[2], my.yl[2], bquote(lambda==.(round(d.l[1,3],3))), cex=0.8, pos=2)
}


# confidence intervals for within-habitat mutlivariate dispersion
ci.d<-function(x, habitat){
  betad <-vegdist(x,method='bray')
  all.hmd <- betadisper(betad,habitat,type='centroid')
  out <- tapply(all.hmd$distances,all.hmd$group,mean)
  return(out)
  }

# habitat-randomized confidence intervals for multivariate dispersion
ci.rep<-function(x, habitat){
  y <- x[sample(1:nrow(x),replace=F),]
  betad <- vegdist(y, method='bray')
  all.hmd <- betadisper(betad, hab2, type='centroid')
  out <- tapply(all.hmd$distances, all.hmd$group, mean)
  return(out)
}
