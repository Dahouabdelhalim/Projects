#load the required library 
library(GenABEL)
library(tidyverse)
library(ggplot2)
library(lattice)
library(latticeExtra)

#read in the results file for each trait (these were obtained using the "LOCO_GWAS.sh" and the "asaMap_GWAS.sh" scripts)
#as an example, results from the LOCO GWAS for SVL are used here
results<-read.csv(file="for_plotting.wald.lambda_adjusted.logp.txt",header=TRUE,sep=' ')

#check that inflation factor is at or close to 1 for adjusted P values
pvals<-results$Pc
pvals_uncorrected<-results$P
estlambda(pvals, plot=FALSE)
estlambda(pvals_uncorrected, plot=FALSE)

###plot results with ggplot (R code adapted from https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/)
nCHR <- length(unique(results$CHR))
results$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(results$CHR)){
  nbp[i] <- max(results[results$CHR == i,]$BP)
  results[results$CHR == i,"BPcum"] <- results[results$CHR == i,"BP"] + s
  s <- s + nbp[i]
}
axis.set <- results %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

#significance thresholds (change the denominator depending on the number of SNPs used in the analysis)
1/113755
0.05/113755

#manhattan plots for GWAS results (as per Fig. 3A and Fig. S23)
ggplot(results, aes(x=BPcum, y=logP)) +
  geom_hline(yintercept = -log10(0.05/113755), linetype="dashed", colour="firebrick", size=0.4)+
  geom_hline(yintercept = -log10(1/113755), linetype="dashed", colour="black", size=0.4)+
  geom_point(aes(color=as.factor(CHR)), size = 1)+
    scale_color_manual(values = rep(c("#576cadff","#6ab7c1ff"), nCHR))+
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0),limits = c(0,7)) +
  labs(x = NULL, y = NULL) + 
  theme_bw() +
  theme(legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank())


###code for QQplots (as per Fig. S23) obtained from https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}

corrected<-qqunif.plot(pvals,draw.conf = FALSE,col="black")
uncorrected<-qqunif.plot(pvals_uncorrected,draw.conf = FALSE, col="orange")
corrected + as.layer(uncorrected)
