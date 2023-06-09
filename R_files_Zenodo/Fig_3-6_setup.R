# Figure setup

# Add in some columns

dat$TRT <- -1
dat[dat$Treatment=="LN",]$TRT <- 0
dat[dat$Treatment=="MN",]$TRT <- 1
dat[dat$Treatment=="HN",]$TRT <- 2
dat[dat$Treatment=="PHN",]$TRT <- 3

TRT <- rep(-1,nrow(dat))
TRT[dat$Treatment=="LN"] <- 0
TRT[dat$Treatment=="MN"] <- 1
TRT[dat$Treatment=="HN"] <- 2
TRT[dat$Treatment=="PHN"] <- 3
jTRT <- jitter(TRT)

xtick <- 0:3

plot_u_se <- function(xoff,ct,cl,pc,logy){
	if(logy){
		points(0+xoff,exp(ct[1,1]+ct[2,1]),col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS LN
		segments(0+xoff,exp(ct[1,1]+ct[2,1]-ct[2,2]),0+xoff,exp(ct[1,1]+ct[2,1]+ct[2,2]),col=cl)
		points(1+xoff,exp(ct[1,1]+ct[3,1]),col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS MN
		segments(1+xoff,exp(ct[1,1]+ct[3,1]-ct[3,2]),1+xoff,exp(ct[1,1]+ct[3,1]+ct[3,2]),col=cl)
		points(2+xoff,exp(ct[1,1]+0),col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS HN
		segments(2+xoff,exp(ct[1,1]+0-ct[1,2]),2+xoff,exp(ct[1,1]+0+ct[1,2]),col=cl)
		points(3+xoff,exp(ct[1,1]+ct[4,1]),col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS PHN
		segments(3+xoff,exp(ct[1,1]+ct[4,1]-ct[4,2]),3+xoff,exp(ct[1,1]+ct[4,1]+ct[4,2]),col=cl)
	}else{
		points(0+xoff,ct[1,1]+ct[2,1],col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS LN
		segments(0+xoff,ct[1,1]+ct[2,1]-ct[2,2],0+xoff,ct[1,1]+ct[2,1]+ct[2,2],col=cl)
		points(1+xoff,ct[1,1]+ct[3,1],col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS MN
		segments(1+xoff,ct[1,1]+ct[3,1]-ct[3,2],1+xoff,ct[1,1]+ct[3,1]+ct[3,2],col=cl)
		points(2+xoff,ct[1,1]+0,col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS HN
		segments(2+xoff,ct[1,1]+0-ct[1,2],2+xoff,ct[1,1]+0+ct[1,2],col=cl)
		points(3+xoff,ct[1,1]+ct[4,1],col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS PHN
		segments(3+xoff,ct[1,1]+ct[4,1]-ct[4,2],3+xoff,ct[1,1]+ct[4,1]+ct[4,2],col=cl)
	}
}
plot_u_se_lme <- function(xoff,ct,bm,cl,pc,logy){
	if(logy){
		points(0+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[2,1]),col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS LN
		segments(0+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[2,1]-ct[2,2]),0+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[2,1]+ct[2,2]),col=cl)
		points(1+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[3,1]),col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS MN
		segments(1+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[3,1]-ct[3,2]),1+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[3,1]+ct[3,2]),col=cl)
		points(2+xoff,exp(ct[1,1]+bm*ct[5,1]+0),col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS HN
		segments(2+xoff,exp(ct[1,1]+bm*ct[5,1]+0-ct[1,2]),2+xoff,exp(ct[1,1]+bm*ct[5,1]+0+ct[1,2]),col=cl)
		points(3+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[4,1]),col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS PHN
		segments(3+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[4,1]-ct[4,2]),3+xoff,exp(ct[1,1]+bm*ct[5,1]+ct[4,1]+ct[4,2]),col=cl)
	}else{
		points(0+xoff,ct[1,1]+bm*ct[5,1]+ct[2,1],col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS LN
		segments(0+xoff,ct[1,1]+bm*ct[5,1]+ct[2,1]-ct[2,2],0+xoff,ct[1,1]+bm*ct[5,1]+ct[2,1]+ct[2,2],col=cl)
		points(1+xoff,ct[1,1]+bm*ct[5,1]+ct[3,1],col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS MN
		segments(1+xoff,ct[1,1]+bm*ct[5,1]+ct[3,1]-ct[3,2],1+xoff,ct[1,1]+bm*ct[5,1]+ct[3,1]+ct[3,2],col=cl)
		points(2+xoff,ct[1,1]+bm*ct[5,1]+0,col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS HN
		segments(2+xoff,ct[1,1]+bm*ct[5,1]+0-ct[1,2],2+xoff,ct[1,1]+bm*ct[5,1]+0+ct[1,2],col=cl)
		points(3+xoff,ct[1,1]+bm*ct[5,1]+ct[4,1],col=cl,pch=pc,cex=1.2) # Mean estimate for ROPS PHN
		segments(3+xoff,ct[1,1]+bm*ct[5,1]+ct[4,1]-ct[4,2],3+xoff,ct[1,1]+bm*ct[5,1]+ct[4,1]+ct[4,2],col=cl)
	}
}
plot_siglets <- function(xoff,yv,siglets,cl){
	text((0:3)+xoff,c(yv,yv,yv,yv),siglets,col=cl)
}