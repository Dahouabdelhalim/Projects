rm(list=ls()) # clear all

source("common_functions.R") # loads common functions

calculate.xF.vd <- function(gamma,dists,muN,nmin,nmax){
	n <- nmin:nmax
	xF <- matrix(data=NA,ncol=length(dists),nrow=length(gamma))
	for(i in 1:length(dists)){
		if(dists[i] == "constant"){
			# calculate xF for const
			xF[,1] <- xF_vd(gamma=gamma,muN=muN,n=NA,pn=NA,dist="constant")
		} else { # calculate xF for pois, geom or waring
			param <- find.param(muN=muN,nmin=nmin,nmax=nmax,dist=dists[i])
			pn <- truncated.dist(param=param,nmin=nmin,nmax=nmax,dist=dists[i])
			for(j in 1:length(gamma)){
				xF[j,i] <- xF_vd(gamma=gamma[j],muN=muN,n=n,pn=pn,dist=dists[i])
			}	
		}
	}
	xF
}

calculate.xdot <- function(x,dists,muN,nmin,nmax,c,b){
	n <- nmin:nmax # group size support
	# gradient of selection (xdot)
	xdot <- matrix(data=NA,ncol=length(dists),nrow=length(x))	
	for(i in 1:length(dists)){
		if(dists[i] == "constant"){ # calculate xdot for const
			xdot[,i] <- x*(1-x)*(b*(1-x)^(muN-1) - c)
		} else { 	# calculate xdot for pois, geom and waring
			param <- find.param(muN,nmin,nmax,dists[i])
			pn <- truncated.dist(param,nmin,nmax,dists[i])
			for(j in 1:length(x)){
				xdot[j,i] <- x[j]*(1-x[j])*(F_vd_b_c(x[j],b,c,muN,n,pn))
			}	
		}
	}
	xdot
}

plot.xF.cb <- function(gamma,xF,dists,ltys,cols,lwd){
	plot(gamma,xF[,1],ylim=c(0,1),type="l",col=cols[1],xlab=expression(paste("cost-to-benefit ratio (",gamma,")")),ylab=NA,las=1,lty=ltys[1],lwd=lwd,cex.axis=2.5,cex.lab=3.0,xaxt="n")
	axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c(0,0.2,0.4,0.6,0.8,1.0),cex.axis=2.5)
	title(ylab=expression(paste("stable equilibrium (",x[F],")")),cex.lab=3,line=5)
	for(i in 2:length(dists)){
		matplot(gamma,xF[,i],type="l",col=cols[i],lty=ltys[i],lwd=lwd,add=TRUE)
	}
	legend(x=0.35,y=0.9,lty=ltys,legend=dists,col=cols,pt.bg=cols,cex=2.5,lwd=5)
}

plot.xF.bc <- function(gamma_inv,dists,ltys,cols,lwd){
	plot(gamma_inv,xF[,1],ylim=c(0,1),type="l",col=cols[1],xlab=expression(paste("benefit-to-cost ratio (",1/gamma,")")),ylab=NA,las=1,lty=ltys[1],lwd=lwd,cex.axis=2.5,cex.lab=3.0,xaxt="n")
	axis(side=1,at=c(0,20,40,60,80,100),labels=c(0,20,40,60,80,100),cex.axis=2.5)
	for(i in 2:length(dists)){
		matplot(gamma_inv,xF[,i],type="l",col=cols[i],lty=ltys[i],lwd=lwd,add=TRUE)
	}
}

plot.xdot <- function(x,xdot,ylim,dists,ltys,cols,lwd,cex,cex.axis,cex.lab){
	plot(x,xdot[,1],type="l",col=cols[1],lty=ltys[1],ylim=ylim,lwd=lwd,xlab=NA,ylab=NA,las=1,xaxt="n",cex.axis=cex.axis)
	for(i in 2:length(dists)){
		matplot(x,xdot[,i],type="l",lty=ltys[i],col=cols[i],lwd=lwd,add=TRUE)
	}
	abline(h=0)
}


# group-size distribution parameters
muN <- 5 # mean group size
nmin <- 2 # min group size
nmax <- 100 # max group size

dists <- c("constant","poisson","geometric","waring") # group-size distributions
ltys <- c(3,4,2,1) # line types
cols <- c("black","blue","orange","red") # line colors

# cost-to-benefit ratios (gamma) and benefit-to-cost ratios (gamma_inv)
gamma <- seq(0,1,len=1000)[2:999]
gamma_inv <- seq(1,100,len=1000)[2:999]

# proportion of Cs (x)
x <- seq(0,1,len=1000)

pdf("fig5.pdf",bg="white",width=14,height=14)
par(oma=c(0,0,0,0),mar=c(7,8,1,1),mgp=c(4,1,0))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

# top left plot
xF <- calculate.xF.vd(gamma=gamma,dists=dists,muN=muN,nmin=nmin,nmax=nmax)
plot.xF.cb(gamma,xF,dists,ltys,cols,lwd=5)

# top right plot
xF <- calculate.xF.vd(gamma=1/gamma_inv,dists=dists,muN=muN,nmin=nmin,nmax=nmax)
plot.xF.bc(gamma_inv,dists,ltys,cols,lwd=5)

# bottom left plot
xdot <- calculate.xdot(x,dists,muN,nmin,nmax,c=1,b=20)
plot.xdot(x,xdot,ylim=c(-0.2,1.2),dists,ltys,cols,lwd=5,cex=2,cex.axis=2.5,cex.lab=2)
title(xlab=expression(paste("fraction of cooperators (",x,")")),cex.lab=3,line=4)
axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c(0,0.2,0.4,0.6,0.8,1.0),cex.axis=2.5)
title(ylab=expression(paste("gradient of selection (",dot(x),")")),cex.lab=3,line=5)
text(0.8,0.1,expression(x[F]),cex=2.5,col="red")
text(0.475,-0.045,expression(x[f]),cex=2.5,col="black")


# bottom right plot
xdot <- calculate.xdot(x,dists,muN,nmin,nmax,c=1,b=5/2)
plot.xdot(x,xdot,ylim=c(-0.25,0.05),dists,ltys,cols,lwd=5,cex=2,cex.axis=2.5,cex.lab=2)
title(xlab=expression(paste("fraction of cooperators (",x,")")),cex.lab=3,line=4)
axis(side=1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c(0,0.2,0.4,0.6,0.8,1.0),cex.axis=2.5)
text(0.1,-0.015,expression(x[F]),cex=2.5,col="red")
text(0.25,0.015,expression(x[f]),cex=2.5,col="black")


dev.off()
