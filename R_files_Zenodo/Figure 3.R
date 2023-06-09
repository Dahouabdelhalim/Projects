### This code produces Figure 3 in the manuscript: 
### Burgess SC, Powell J, Bueno M. Dispersal, kin aggregation, and the fitness consequences of not spreading sibling larvae.
# Code finalized Jan 2022
# Any comments or error reporting, please contact Scott Burgess. sburgess@bio.fsu.edu

# R version 4.0.5 (2021-03-31) -- "Shake and Throw"

# Load required libraries
library("vegan") # version 2.5-7. Needed for Clark Evans Corrected function.R. 
library("spatstat") # version 2.3-0. Needed for clarkevans() and clarkevans.test()

source("Clark Evans Corrected function.R") # Uses 'vegan' and spatstat', but produces more summary information
# ClarkEvans() comes from 'Clark Evans Corrected function.R', which uses clarkevans() from library(spatstat)


# Import data
# set.wd() # set working directory, or file path, before importing data
dat.all <- read.csv("Figure 3 data.csv")

# Change 'Plate.ID' to a factor
dat.all$Plate.ID <- factor(dat.all$Plate.ID)
# summary(dat.all)

# Prepare indices
X.col <- which(names(dat.all) %in% c("True.X"))
Y.col <- which(names(dat.all) %in% c("True.Y"))
Plate.ID.col <- which(names(dat.all) %in% c("Plate.ID"))
IDcols.to.use <- c("D.Date","R.Date","D.Group","Plate.ID", "Location")
IDcols.to.use12 <- c(IDcols.to.use,"Deployment","Plate.ID.Deployment")
Lquad <- 153 # length of one side of quadrat in mm, to calculate realized settlement area 




# Subset data for Deployment = 1 (first three days in the field)
dat1 <- dat.all[dat.all$Deployment==1,]

# Prepare data frames to add data to.
tmp1 <- dat1[!duplicated(dat1$Plate.ID),which(names(dat1) %in% IDcols.to.use)]
dat.ag <- cbind(tmp1,matrix(NA,nrow=nrow(tmp1),ncol=13)) # ncol=13 because there are 13 values outputted by ClarkEvans() 
plate.ind <- unique(dat1$Plate.ID) # Vector of Plate.ID's


# Calculate Clark Evans index within plates
# R = 1 = random. R < 1 = aggregated. R > 1 = regular
# The for loop calculates indices for Day 1 - 3
for(i in 1:nrow(tmp1)){
	tmp.dat <- dat1[dat1$Plate.ID==plate.ind[i],]
	hold <- ClarkEvans(coords=tmp.dat[,X.col:Y.col],quadratsize=Lquad,boundary.strip=FALSE) 
	dat.ag[i,] <- cbind(tmp1[i,],hold)		
}
names(dat.ag) <- c(names(tmp1),names(hold))




# Subset data for only those plates with 2 Deployments 
plates12 <- unique(dat.all$Plate.ID[which(dat.all$Deployment==2)]) # Plate.ID's that had second deployment
dat12 <- dat.all[dat.all$Plate.ID%in%plates12,]
dat12$Plate.ID.Deployment <- interaction(dat12$Plate.ID, dat12$Deployment)

# Prepare data frames to add data to.
tmp1 <- dat12[!duplicated(dat12$Plate.ID.Deployment),which(names(dat12) %in% IDcols.to.use12)]
dat.ag12 <- cbind(tmp1,matrix(NA,nrow=nrow(tmp1),ncol=13)) # ncol=13 because there are 13 values outputted by ClarkEvans() 
plate.ind <- unique(tmp1$Plate.ID.Deployment)

# Calculate Clark Evans index within plates
# R = 1 = random. R < 1 = aggregated. R > 1 = regular
# The for loop calculates indices for settlement over Day 1 - 3 and Day 4 - 6
for(i in 1:nrow(tmp1)){
	tmp.dat <- dat12[dat12$Plate.ID.Deployment==plate.ind[i],]
	hold <- ClarkEvans(coords=tmp.dat[,X.col:Y.col],quadratsize=Lquad,boundary.strip=FALSE) 
	dat.ag12[i,] <- cbind(tmp1[i,],hold)	
}
names(dat.ag12) <- c(names(tmp1),names(hold)) 





# Calculate Clark Evans index within plates, settlement over day 1 - 6.
tmp1 <- dat12[!duplicated(dat12$Plate.ID),which(names(dat12) %in% IDcols.to.use)]
dat.ag1plus2 <- cbind(tmp1,matrix(NA,nrow=nrow(tmp1),ncol=13)) # ncol=13 because there are 13 values outputted by ClarkEvans()
plate.ind <- unique(tmp1$Plate.ID)

# Calculate Clark Evans index within plates
# R = 1 = random. R < 1 = aggregated. R > 1 = regular
# The for loop calculates indices for Day 1 - 6
for(i in 1:nrow(tmp1)){
	tmp.dat <- dat12[dat12$Plate.ID==plate.ind[i],]
	hold <- ClarkEvans(coords=tmp.dat[,X.col:Y.col],quadratsize=Lquad,boundary.strip=FALSE) 
	dat.ag1plus2[i,] <- cbind(tmp1[i,],hold)
}
names(dat.ag1plus2) <- c(names(tmp1),names(hold))



# Statements in Results section
m1 <- lm(settler.density*100/3 ~ D.Group-1, data=dat.ag)
m1data <- data.frame(estimate=round(coef(m1),3),round(confint(m1),3))

m1data[which(m1data$estimate==min(m1data$estimate)),] # 0.024 (0.011 –  0.037, 95% CI)
m1data[which(m1data$estimate==max(m1data$estimate)),] # 0.058 (0.048 –  0.069, 95% CI)
# The average settlement density varied ~2.5-fold among the seven deployments, from 0.024 (0.011 –  0.037, 95% CI) in one deployment to 0.058 (0.048 –  0.069, 95% CI) in another deployment,


m1 <- lm(settler.number ~ D.Group-1, data=dat.ag)
m1data <- data.frame(estimate=round(coef(m1),3),round(confint(m1),3))
m1data[which(m1data$estimate==min(m1data$estimate)),] # 17 (7.788 –  26.212, 95% CI)
m1data[which(m1data$estimate==max(m1data$estimate)),] # 41 (33.476 –  48.522, 95% CI)
# which corresponds to 17 (7.788 – 26.212, 95% CI) to 41 (33.478 – 48.522, 95% CI) settlers per 15 × 15cm settlement plate.




# Make Figure 3 
cex.size <- 2.5
cex.lab <- 1.5

cols <- c("#009E73","#0072B2","#D55E00") # green, blue, red  (color-blind friendly)

quartz(width=12,height=4)
par(mar=c(2,3,4,1),oma=c(3,5,1,1),mfrow=c(1,3))

#a) 
plot(c(0,0.09),c(0.6,1.3),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=seq(0,0.1,0.01),cex.axis=2)
axis(side=2,at=seq(0,2,0.1),cex.axis=2,las=1)
abline(h=1,lty=2)
with(dat.ag,points(settler.density*100/3, R.clarkevans, 
	pch=19,
	# col=ifelse(p.value.agg.clarkevans<0.05,adjustcolor(cols[3],alpha.f=0.6),
			# ifelse(p.value.reg.clarkevans<0.05,	adjustcolor(cols[1],alpha.f=0.6),	
			# adjustcolor(cols[2],alpha.f=0.6))
	col=ifelse(p.value.agg.clarkevans<0.05,cols[3],
			ifelse(p.value.reg.clarkevans<0.05,	cols[1],	
			cols[2])
			),
	cex=cex.size))
mtext(expression(paste("Number of settlers (",cm^-2,d^-1,")",sep="")),
	side=1,line=4,cex=cex.lab)
mtext(expression(paste("Index of aggregation (",italic("R"),")",sep="")),side=2,line=2,cex=cex.lab,outer=T)
mtext("a)",adj=0,line=1,side=3,cex=cex.lab)

# b)
plot(c(0.5,2.5),c(0.6,1.3),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=c(1,2),labels=c("Day 1-3", "Day 1-6"),cex.axis=2)
axis(side=2,at=seq(0,2,0.1),cex.axis=2,las=1)
lines(c(0,2.1),c(1,1),lty=2)

foo <- cbind.data.frame(dat.ag12[dat.ag12$Deployment==1,], 
	R16=dat.ag1plus2$R.clarkevans,
	p.value.agg16=dat.ag1plus2$p.value.agg.clarkevans,
	p.value.reg16=dat.ag1plus2$p.value.reg.clarkevans)

with(foo,points(rep(1,length(R.clarkevans)), R.clarkevans, 
	pch=19,
	col=ifelse(p.value.agg.clarkevans<0.05,cols[3],
			ifelse(p.value.reg.clarkevans<0.05,	cols[1],	
			cols[2])),
	cex=cex.size))

with(foo,points(rep(2,length(R16)), R16,
	pch=19,
	col=ifelse(p.value.agg16 <0.05,cols[3],
			ifelse(p.value.agg16 <0.05,	cols[1],	
			cols[2])),
	cex=cex.size))

plate.ind <- unique(dat.ag12$Plate.ID)
for(i in 1:length(plate.ind)){
with(foo[foo$Plate.ID==plate.ind[i],],lines(c(1,2), c(R.clarkevans,R16)))}
# lines(c(0,2.2),c(1,1),lty=2)
# abline(h=1,lty=2)
mtext("Deployment duration",side=1,line=3.5,cex=cex.lab)
mtext("b)",adj=0,line=1,side=3,cex=cex.lab)

# c)
plot(c(0,0.09),c(0,2),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=seq(0,0.1,0.01),cex.axis=2)
axis(side=2,at=seq(0,2,0.5),cex.axis=2,las=1)
with(dat.ag,points(settler.density*100/3,min.NN/10,
	pch=19,
	# col=ifelse(p.value.agg.clarkevans<0.05,adjustcolor(cols[3],alpha.f=0.6),
			# ifelse(p.value.reg.clarkevans<0.05,	adjustcolor(cols[1],alpha.f=0.6),	
			# adjustcolor(cols[2],alpha.f=0.6))
	col=ifelse(p.value.agg.clarkevans<0.05,cols[3],
			ifelse(p.value.reg.clarkevans<0.05,	cols[1],	
			cols[2])
			),
	cex=cex.size))
mtext(expression(paste("Number of settlers (",cm^-2,d^-1,")",sep="")),
	side=1,line=4,cex=cex.lab)
mtext("Average minimum distance\\n to nearest neighbor (cm)",side=2,line=4,adj=-0.1,cex=cex.lab)
mtext("c)",adj=0,line=1,side=3,cex=cex.lab)

legend('topright',legend=c("Regular","Random","Aggregrated"),
		bty="n",
		pt.cex=2,
		cex=1.5,
		pch=19,
		col=c(cols[1],
			cols[2],
			cols[3]),		
		)
