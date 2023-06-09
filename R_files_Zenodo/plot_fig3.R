rm(list=ls()) # clear all

source("common_functions.R") # loads common functions

calculate.xdot <- function(x,dists,muN,nmin,nmax,c,b,w){
	n <- nmin:nmax # group size support
	# gradient of selection (xdot)
	xdot <- matrix(data=NA,ncol=length(dists),nrow=length(x))	
	for(i in 1:length(dists)){
		if(dists[i] == "constant"){ # calculate xdot for const
			xdot[,i] <- x*(1-x)*((b/muN)*(1-x+w*x)^(muN-1) - c)
		} else { # calculate xdot for pois, geom and waring
			param <- find.param(muN,nmin,nmax,dists[i])
			pn <- truncated.dist(param,nmin,nmax,dists[i])
			for(j in 1:length(x)){
				xdot[j,i] <- x[j]*(1-x[j])*(F_pggds_b_c(x[j],w,b,c,muN,n,pn))
			}	
		}
	}
	xdot
}

plot.xdot <- function(x,xdot,ylim,dists,ltys,cols,lwd,cex,cex.axis,cex.lab){
	plot(x,xdot[,1],type="l",col=cols[1],lty=ltys[1],ylim=ylim,lwd=lwd,xlab=NA,ylab=NA,las=1,xaxt="n",cex.axis=cex.axis)
	for(i in 2:length(dists)){
		matplot(x,xdot[,i],type="l",lty=ltys[i],col=cols[i],lwd=lwd,add=TRUE)
	}
	axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c(0,0.2,0.4,0.6,0.8,1.0),cex.axis=cex.axis)
	title(ylab=expression(paste("gradient of selection (",dot(x),")")),cex.lab=cex.lab,line=5)
	abline(h=0)
}

# group-size distribution parameters
muN <- 5 # mean group size
nmin <- 2 # min group size
nmax <- 100 # max group size

dists <- c("constant","poisson","geometric","waring") # group-size distributions
ltys <- c(3,4,2,1) # line types
cols <- c("black","blue","orange","red") # line colors

# proportion of Cs (x)
x <- seq(0,1,len=1000)

# plot figure
pdf("fig3.pdf",bg="white",width=7,height=2*7)
par(oma=c(0,0,0,0),mar=c(5,7,1,1),mgp=c(4,1,0))
layout(matrix(c(1,2), 2, 1, byrow = TRUE))

# top panel
xdot <- calculate.xdot(x,dists,muN,nmin,nmax,c=1,b=15,w=0.5)
plot.xdot(x,xdot,ylim=c(-0.1,0.22),dists,ltys,cols,lwd=5,cex=2,cex.axis=1.5,cex.lab=2)
text(0.875,0.02,expression(x[F]),cex=2.5,col="red")
text(0.425,-0.015,expression(x[f]),cex=2.5,col="black")
legend(x=0.6,0.2,lty=c(3,4,2,1),lwd=3,legend=dists,col=cols,cex=1.5,pt.lwd=5)


# bottom panel
xdot <- calculate.xdot(x,dists,muN,nmin,nmax,c=1,b=1.5,w=1.5)
plot.xdot(x,xdot,ylim=c(-0.1,0.1),dists,ltys,cols,lwd=5,cex=2,cex.axis=1.5,cex.lab=2)
legend(x=0.6,0.2,lty=ltys,lwd=3,legend=dists,col=cols,cex=1.5,pt.lwd=5)
title(xlab=expression(paste("fraction of cooperators (",x,")")),cex.lab=2,line=3)
text(0.25,0.01,expression(x[F]),cex=2.5,col="red")
text(0.75,-0.01,expression(x[f]),cex=2.5,col="black")

dev.off()


