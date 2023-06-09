library(svd)
library(dplyr)
library(plyr)
library(tsne)
library(Rtsne)
library(RColorBrewer)

# --------------------------------------------------
# get variable genes from normalized UMI counts
# --------------------------------------------------
# m: matrix normalized by UMI counts
.get_variable_gene<-function(m) {
  
  df<-data.frame(mean=colMeans(m),cv=apply(m,2,sd)/colMeans(m),var=apply(m,2,var))
  df$dispersion<-with(df,var/mean)
  df$mean_bin<-with(df,cut(mean,breaks=c(-Inf,quantile(mean,seq(0.1,1,0.05)),Inf)))
  var_by_bin<-ddply(df,"mean_bin",function(x) {
    data.frame(bin_median=median(x$dispersion),
               bin_mad=mad(x$dispersion))
  })
  df$bin_disp_median<-var_by_bin$bin_median[match(df$mean_bin,var_by_bin$mean_bin)]
  df$bin_disp_mad<-var_by_bin$bin_mad[match(df$mean_bin,var_by_bin$mean_bin)]
  df$dispersion_norm<-with(df,abs(dispersion-bin_disp_median)/bin_disp_mad)
  df
}
# -----------------------------------------------------
# compute top n principle components with propack 
# -----------------------------------------------------
.do_propack <- function(x,n) {
  use_genes <- which(colSums(x) > 1)
  m <- x[,use_genes]
  bc_tot <- rowSums(m)
  median_tot <- median(bc_tot)
  m <- sweep(m, 1, median_tot/bc_tot, '*')
  m <- sweep(m, 2, colMeans(m), '-')
  m <- sweep(m, 2, apply(m, 2, sd), '/')
  ppk<-propack.svd(as.matrix(m),neig=n)
  pca<-t(ppk$d*t(ppk$u))
  list(ppk=ppk,pca=pca, m=m,use_genes=use_genes)
}
# ---------------------------------------------
# normalize the gene barcode matrix by umi
# ---------------------------------------------
.normalize_by_umi <-function(x) {
  cs <- colSums(x)
  x_use_genes <- which(cs >= 1)
  x_filt<-x[,x_use_genes]
  rs<-rowSums(x_filt)
  rs_med<-median(rs)
  x_norm<-x_filt/(rs/rs_med)
  list(m=x_norm,use_genes=x_use_genes)
}



doTSNE <- function(mat, perplexity){
  matt <- t(mat)
  l<-.normalize_by_umi(matt)
  m_n<-l$m
  df<-.get_variable_gene(m_n) 
  disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
  df$used<-df$dispersion_norm >= disp_cut_off
  set.seed(0)
  m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]
  #m_n_1000 <- t(t(m_n_1000)-hkmean)
  tsne_n_1000<-Rtsne(1-cor(t(m_n_1000),method="spearman"),pca=T,check_duplicates = FALSE, perplexity= perplexity)
  #tsne_n_1000<-Rtsne(m_n_1000,pca=T)
  tsn <- tsne_n_1000$Y
  rownames(tsn)<-colnames(mat)
  return(tsn)
}

plotGene <- function(mat, tsne, genes){
  index <- mat[which(rownames(mat) %in% genes),]
  if(length(genes)>1){
    index <- apply(index, 2, mean)
  }
  
  #zscore normalize
  index <- (index - mean(index))/sd(index)
  threshold <- 2
  index[which(index>threshold)] <- threshold
  index[which(index < -threshold)] <- -threshold
  
  #color it up
  rbPal <-  colorRampPalette(c("#F8FA0D", "#F6DA23", "#F8BA43","#A5BE6A","#2DB7A3","#1389D2", "#352A86"))
  Col <- rev(rbPal(100))[as.numeric(cut(index,breaks = 100))]
  par(mar=c(5.1, 5.1, 4.1, 11.1), xpd=TRUE)
  plot(tsne,pch=16,cex=1.25,col=Col,las=1,xlab="tSNE1",ylab="tSNE2",cex.axis=1.5,cex.lab=1.5)
  
}

plotMetric<- function(mat, tsne, metric){
  index <- metric
  
  #zscore normalize
  index <- (index - mean(index))/sd(index)
  threshold <- 2
  index[which(index>threshold)] <- threshold
  index[which(index < -threshold)] <- -threshold
  
  #color it up
  rbPal <-  colorRampPalette(c("#F8FA0D", "#F6DA23", "#F8BA43","#A5BE6A","#2DB7A3","#1389D2", "#352A86"))
  Col <- rev(rbPal(100))[as.numeric(cut(index,breaks = 100))]
  par(mar=c(5.1, 5.1, 4.1, 11.1), xpd=TRUE)
  plot(tsne,pch=16,cex=1.25,col=Col,las=1,xlab="tSNE1",ylab="tSNE2",cex.axis=1.5,cex.lab=1.5)
  
}

library(RColorBrewer)
plotCategory<- function(tsne, category, cols = NULL, inset = -0.6, labels = NULL){
  cat <- category
  assign <- sort(unique(cat))
  if(!is.null(labels)){
    assign <- labels
  }
  pal <- colorRampPalette(brewer.pal(length(assign), "Set1"))(length(assign))
  if(!is.null(cols)){
    pal <- adjustcolor(cols, alpha.f = 0.75)
  }
  names(pal) <- assign
  cols2 <- pal[cat]
  
  par(mar=c(5.1, 6.1, 4.1, 11.1), xpd=TRUE)
  plot(tsne,pch=16,cex=0.75,col=cols2,las=1,xlab="",ylab="",cex.axis=1.5,cex.lab=1.5)
 
   legend("topright",inset=c(inset,0), y.intersp = 0.9, legend= assign, ncol = 1, pch = 19, pt.cex = 2.5, 
         col = adjustcolor(pal, alpha.f =1),bty = "n", cex = 1.25)
  
}

#zscore normalize
#index <- (index - mean(index))/sd(index)
#threshold <- 0.5
#index[which(index>threshold)] <- threshold
#index[which(index < -threshold)] <- -threshold
#
##color it up
#rbPal <-  colorRampPalette(c("#F8FA0D", "#F6DA23", "#F8BA43","#A5BE6A","#2DB7A3","#1389D2", "#352A86"))
#Col <- rev(rbPal(100))[as.numeric(cut(index,breaks = 100))]
#par(mar=c(5.1, 5.1, 4.1, 11.1), xpd=TRUE)
#plot(m,pch=16,cex=1.25,col=Col,las=1,xlab="tSNE1",ylab="tSNE2",cex.axis=1.5,cex.lab=1.5)
#
