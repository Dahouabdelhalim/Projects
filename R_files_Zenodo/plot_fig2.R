rm(list=ls()) # clears all

### left panel ###

source("common_functions.R") # loads common functions

# calculate gamma2 = E[w^(N-1)]/muN
calculate.gamma2 <- function(dists,muN,nmin,nmax,w){
	n <- nmin:nmax
	gamma2 <- matrix(data=NA,ncol=length(dists),nrow=length(w))
	for(i in 1:length(dists)){
		if(dists[i] == "constant"){ # constant
			gamma2[,1] <- w^(muN-1)/muN
		} else { # poisson, geometric or waring
			param <- find.param(muN,nmin,nmax,dists[i])
			pn <- truncated.dist(param,nmin,nmax,dists[i])
			for(j in 1:length(w)){
				gamma2[j,i] <- sum(pn*w[j]^n)/(w[j]*muN)
			}
		}
	}
	gamma2
}

# plot gamma2
plot.gamma2 <- function(w,gamma2,dists,ltys,cols,lwd=5,cex=2,cex.axis=1.5,cex.lab=1.5){
	plot(w,gamma2[,1]*muN,log="y",ylim=c(0.01,100),type="l",lty=ltys[1],col=cols[1],
xlab=NA,ylab=NA,xaxs="i",yaxs="i",las=1,lwd=lwd,yaxt="n",cex.axis=cex.axis,cex.lab=cex.lab)
	for(i in 2:length(dists)){
		matplot(w,gamma2[,i]*muN,log="y",ylim=c(0.01,100),type="l",lty=ltys[i],col=cols[i],lwd=lwd,add=TRUE)
	}
	abline(v=1)
	abline(h=1,lwd=lwd,lty=1,col="black")
	axis(2, at=c(0.01,0.1,1,10,100),labels=c(expression(10^-2),expression(10^-1),expression(10^0),expression(10^1),expression(10^2)), las=1, pos=0, cex.axis=cex.axis)
	legend(x=1.05,y=0.25,lty=c(1,3,4,2,1),lwd=2,legend=c(expression(gamma[1]),expression(paste(gamma[2]," (constant)")),expression(paste(gamma[2]," (poisson)")),expression(paste(gamma[2]," (geometric)")),expression(paste(gamma[2]," (waring)"))),col=c("black","black","blue","orange","red"),cex=cex.axis)
	text(x=0.5,y=10,cex=cex.axis,"defection")
	text(x=1.7,y=2,cex=cex.axis,"bi-stability")
	text(x=1.5,y=0.45,cex=cex.axis,"cooperation")
	text(x=0.4,y=0.7,cex=cex.axis,"co-existence")
	title(xlab=expression(paste("synergy/discounting (",w,")")),cex.lab=cex.lab,line=3)
	title(ylab=expression(paste("normalized cost-to-benefit ratio (",gamma,mu[N],")")),cex.lab=cex.lab,line=4)
}

### right panel ###

library(lattice) # (for function levelplot)
library(reshape) # (for function melt)

# define palette of colors
rgb.palette.1 <- colorRampPalette(c("red", "orange", "yellow", "green", "light blue", "blue"), space = "rgb")

# load data
load("results_pggds.RData")

# melt results
results <- melt(results, id=c("w","gammamuN")) 

# params
cex.axis <- 1.5
cex.lab <- 1.5

# create lattice object
lplot <- levelplot(value~w*log10(gammamuN)|variable,
data=results,
col.regions=heat.colors(20),
# aspect=1,
colorkey=list(labels=list(cex=1.5)),
as.table=TRUE,
layout=c(2,2),
par.strip.text=list(cex=cex.axis),
xlab=list(label=expression(paste("synergy/discounting (",w,")")), cex=cex.lab),
ylab=list(label=expression(paste("normalized cost-to-benefit ratio (",gamma,mu[N],")")), cex=cex.lab),
scales=list(axs='i',x=list(at=c(0,0.5,1,1.5,2),limits=c(0,2)),y=list(at=c(-2,-1,0,1,2),limits=c(-2,2),labels=c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2}))),cex=cex.axis),
	panel=function(...){
		panel.levelplot.raster(...)
		panel.lines(x=c(0,2),y=c(0,0),type="l",lty="solid",col="black")
		panel.lines(x=c(1,1),y=c(-2,2),type="l",lty="dotdash",col="black")
	}
)

### plot both left and right panel together ###

library(gridBase) # for functions pushViewport, popViewport, gridOMI, etc

pdf("fig2.pdf",bg="white",width=2*7,height=7)
par(mar=c(4,6,3,3))

plot.new()

pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(rep(1,2), c("null", "null")))))

pushViewport(viewport(layout.pos.col = 1))

par(omi = gridOMI(), new = TRUE)

muN <- 5 # mean group size
nmin <- 2 # min group size
nmax <- 100 # max group size
w <- seq(0,2,len=200) # synergy/discounting w
dists <- c("constant","poisson","geometric","waring") # group-size distributions
ltys <- c(3,4,2,1) # line types
cols <- c("black","blue","orange","red") # line colors
gamma2 <- calculate.gamma2(dists,muN,nmin,nmax,w) # calculate gamma2
plot.gamma2(w,gamma2,dists,ltys,cols,lwd=5,cex=2,cex.axis=1.5,cex.lab=1.5) # plot gamma2

popViewport()
pushViewport(viewport(layout.pos.col = 2))

print(lplot, newpage = FALSE)

popViewport(2)

dev.off()



