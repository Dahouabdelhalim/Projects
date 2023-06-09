### This code produces Figure 2 in the manuscript: 
### Burgess SC, Powell J, Bueno M. Dispersal, kin aggregation, and the fitness consequences of not spreading sibling larvae.
# Code finalized Jan 2022
# Any comments or error reporting, please contact Scott Burgess. sburgess@bio.fsu.edu

# R version 4.0.5 (2021-03-31) -- "Shake and Throw"

# Load required libraries
library(plotrix) # version 3.8-2; used for polar.plot()

# Import data
# set.wd() # set working directory, or file path, before importing data
dat <- read.csv("Figure 2 data.csv")

# Change 'deployment' to a factor
dat$deployment <- as.factor(dat$deployment)
# summary(dat)

# Summarize dat
# Table 1; number of settlement plates per deployment
with(dat,table(deployment,treatment))

# Figure 1 caption; Number of poles per distance 
with(dat[dat$deployment==2,],table(distance.m))

# Number of settlement plates per distance per deployment
with(dat,table(deployment,distance.m))



# Settlement area of each settlement plate
settlement.area <- 21.5 * 7 # length cm x width cm; cm^2



# Calculate the average settlement per plate per day (settlers.day) over all 'Control' deployments
datControl <- dat[dat$treatment=="Control",]
mean.settlers <- mean(datControl$settlers.day) # over all three control deployments, there were on average 0.39 settlers cm-2 day-1.
# with(datControl, aggregate(settlers.day,by=list(Distance=distance.m),mean))
# with(datControl, aggregate(settlers.day,by=list(Deployment=deployment),mean))

qpois(p=0.95,lambda=mean.settlers)  # 95% probability of observing 2 or fewer settlers per cm2 per day from background

threshold <- 2 # upper limit to the expected average 'background' level of settlement



# Prepare Figure 2
max.settlement <- 330 # define maximum value, for plotting
ymax <- max.settlement / settlement.area / 3 / 7 # per cm2 per day (=3) per colony (=7) 

radius = c(0,
	0.25,
	0.5,
	1,
	2,
	4,
	8,
	12)
pole.number.at.each.distance <- c(1,
								4,
								8-(5-1),
								16-(9-1),
								28-(17-1),
								51-(29-1),
								97-(52-1),
								170-(98-1))
circ <- function(r){2*pi*r} # function to calculate circumference
arc.lengths = circ(r=radius) / pole.number.at.each.distance  # distance between poles, given pole.number.at.each.distance
# pole.number.at.each.distance / circ(radius) # poles per meter of circumference at each distance = sampling effort per circumference
# sum(pole.number.at.each.distance) # total number of poles
# theta = arc length / radius

degree.func <- function(L,R){
	ifelse(R>0, (180*L) / (pi*R), 0)
}


###### FIGURE 2
# Set colors
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
cols <- c("#999999","#D55E00")

dev.new(width=8,height=3.5)
par(oma=c(1,3,1,0),mfrow=c(1,3))

# a)
# Add lines
par(mar=c(0,2,0,1))
polar.plot(	lengths=NA,
			polar.pos=NA,
			radial.lim=radius,
			show.radial.grid=F,
			show.grid.labels=F,
			grid.unit="meters",
			labels="",
			point.symbols=19,
			start=270,clockwise=T,rp.type="s")
for(i in 1:length(radius)){
	theta = degree.func(arc.lengths[i],radius[i])
	angles.of.the.poles <- (theta * 0:pole.number.at.each.distance[i])
	radius.vec <- rep(radius[i],length(angles.of.the.poles)) # the radius of the poles
polar.plot(	lengths=radius.vec,
			polar.pos=angles.of.the.poles,
			radial.lim=radius,
			cex=1, 
			point.col='black',
			show.radial.grid=F,
			show.grid.labels=F,
			labels="",
			point.symbols=19,
			start=270,clockwise=T,rp.type="s",add=T)
}
# text(0,-13.5,"12m",cex=1.5)					
# text(0,-9.5,"8m",cex=1.5)					
# text(0,-5.5,"4m",cex=1.5)					
# text(0,-2.9,"2m",cex=1.2)					
axis(side=1,at=seq(-12,12,2),cex.axis=1.3)
mtext(side=1,cex=1.2,line=3,text="X distance (meters)")
axis(side=2,at=seq(-12,12,2),cex.axis=1.3,las=1)
mtext(side=2,cex=1.2,line=3,text="Y distance (meters)")
mtext(side=3,adj=0,cex=1.2,text="a) Settlement plate array")

# b)
par(mar=c(0,6,0,0))
plot(c(0,12),c(0,ymax),type="n",bty="l",ylab="",xlab="",xaxt="n",yaxt="n")
with(dat[dat$deployment =="2" & dat$settlers.day<=threshold,], 
	points(distance.m,settlers/settlement.area/3/7,
	cex=2,
	pch=19,
	# col=adjustcolor(cols[1],alpha.f=0.6)
	col=cols[1]
	))
with(dat[dat$deployment =="2" & dat$settlers.day>threshold,], 
	points(distance.m,settlers/settlement.area/3/7,
	cex=2,
	pch=19,
	# col=adjustcolor(cols[2],alpha.f=0.6)
	col=cols[2]
	))

axis(side=1,at=radius,cex.axis=1.5)
axis(side=2,at=seq(0,ymax,0.02),cex.axis=1.5,las=1)
mtext(side=2,line=5.5,cex=1.2,adj=0.5,text="Number of settlers")
mtext(side=2,line=3.5,cex=1.2,adj=0.5,expression(paste("(",cm^-2,d^-1,colony^-1,")",sep="")))
mtext(side=3,adj=0,line=0,cex=1.2,"b) Trial 1")

# c)
par(mar=c(0,4,0,2))
plot(c(0,12),c(0, 0.04),type="n",bty="l",ylab="",xlab="",xaxt="n",yaxt="n")
with(dat[dat$deployment =="4" & dat$settlers.day<=threshold,], 
	points(distance.m,settlers/settlement.area/3/7,
	cex=2,
	pch=19,
	# col=adjustcolor(cols[1],alpha.f=0.6)
	col=cols[1]
	))
with(dat[dat$deployment =="4" & dat$settlers.day>threshold,], 
	points(distance.m,settlers/settlement.area/3/7,
	cex=2,
	pch=19,
	# col=adjustcolor(cols[2],alpha.f=0.6)
	col=cols[2]
))

axis(side=1,at=radius,cex.axis=1.5)
axis(side=2,at=seq(0,0.04,0.01),cex.axis=1.5,las=1)
mtext(side=3,adj=0,line=0,cex=1.2,"c) Trial 2")

mtext(side=1,line=-2,cex=1.2,outer=T,adj=0.9,"Dispersal distance (meters from maternal colony)")

	