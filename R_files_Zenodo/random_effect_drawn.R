##Random Effect Plotting
setwd("~/Dropbox/Things we generated during Shuxi's visit/R code/Random effects plot code")
library(ggplot2)
rm(list=ls())
load("plot_random_effect_draw.RData")

df=data.frame(name=cluster.name,
              mean=apply(cluster_rd_m,2,mean),
              min=apply(cluster_rd_m,2,FUN=function(x){quantile(x,0.025)}),
              max=apply(cluster_rd_m,2,FUN=function(x){quantile(x,0.975)}))

#df$name=factor(df$name,levels=df$name[order(df$mean)])
df$name=factor(df$name,levels=df$name)
p <- ggplot(df,aes(name,mean,shape=name))

#Added horizontal line at y=0, error bars to points and points with size two
p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=min, ymax=max), width=0,color="black") + geom_point(aes(size=2)) 

#Removed legends and with scale_shape_manual point shapes set to 1 and 16
p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=rep(1,12))

#Changed appearance of plot (black and white theme) and x and y axis labels
p <- p + theme_bw() + xlab("Group identification") + ylab("Mediator Effect size")
p <- p+ coord_flip()
pdf("../../Figures/Random_Effect/grp_rd_on_dsi_f.pdf",width=6,height=6)
p
dev.off()

df=data.frame(name=cluster.name,
              mean=apply(cluster_rd_y,2,mean),
              min=apply(cluster_rd_y,2,FUN=function(x){quantile(x,0.025)}),
              max=apply(cluster_rd_y,2,FUN=function(x){quantile(x,0.975)}))

#df$name=factor(df$name,levels=df$name[order(df$mean)])
df$name=factor(df$name,levels=df$name)
p <- ggplot(df,aes(name,mean,shape=name))

#Added horizontal line at y=0, error bars to points and points with size two
p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=min, ymax=max), width=0,color="black") + geom_point(aes(size=2)) 

#Removed legends and with scale_shape_manual point shapes set to 1 and 16
p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=rep(1,12))

#Changed appearance of plot (black and white theme) and x and y axis labels
p <- p + theme_bw() + xlab("Group identification") + ylab("Outcome Effect size")
p <- p+ coord_flip()
pdf("../../Figures/Random_Effect/grp_rd_on_gc.pdf",width=6,height=6)
p
dev.off()

df=data.frame(name=hydro.cluster.name,
              mean=apply(hydro_cluster_rd_m,2,mean),
              min=apply(hydro_cluster_rd_m,2,FUN=function(x){quantile(x,0.025)}),
              max=apply(hydro_cluster_rd_m,2,FUN=function(x){quantile(x,0.975)}))

#df$name=factor(df$name,levels=df$name[order(df$mean)])
df$name=factor(df$name,levels=df$name)
p <- ggplot(df,aes(name,mean,shape=name))

#Added horizontal line at y=0, error bars to points and points with size two
p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=min, ymax=max), width=0,color="black") + geom_point(aes(size=2)) 

#Removed legends and with scale_shape_manual point shapes set to 1 and 16
p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=rep(1,19))

#Changed appearance of plot (black and white theme) and x and y axis labels
p <- p + theme_bw() + xlab("Hydrological year") + ylab("Mediator Effect size")
p <- p+ coord_flip()
pdf("../../Figures/Random_Effect/hydro_rd_on_dsi_f.pdf",width=6,height=6)
p
dev.off()


df=data.frame(name=hydro.cluster.name,
              mean=apply(hydro_cluster_rd_y,2,mean),
              min=apply(hydro_cluster_rd_y,2,FUN=function(x){quantile(x,0.025)}),
              max=apply(hydro_cluster_rd_y,2,FUN=function(x){quantile(x,0.975)}))

#df$name=factor(df$name,levels=df$name[order(df$mean)])
df$name=factor(df$name,levels=df$name)
p <- ggplot(df,aes(name,mean,shape=name))

#Added horizontal line at y=0, error bars to points and points with size two
p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=min, ymax=max), width=0,color="black") + geom_point(aes(size=2)) 

#Removed legends and with scale_shape_manual point shapes set to 1 and 16
p <- p + guides(size=FALSE,shape=FALSE) + scale_shape_manual(values=rep(1,19))

#Changed appearance of plot (black and white theme) and x and y axis labels
p <- p + theme_bw() + xlab("Hydrological year") + ylab("Outcome Effect size")
p <- p+ coord_flip()
pdf("../../Figures/Random_Effect/hydro_rd_on_gc.pdf",width=6,height=6)
p
dev.off()



