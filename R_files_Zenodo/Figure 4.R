### This code produces Figure 4 in the manuscript: 
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
dat.all <- read.csv("Figure 4 data.csv")

# Change 'Sheet.ID' to a factor
dat.all$Sheet.ID <- factor(dat.all$Sheet.ID)
# summary(dat.all)


# # Summarize structure of the data
# with(dat.all,table(Relatedness,Sheet.ID)) 
# datstr <- with(dat.all, aggregate(True.X,by=list(Sheet.ID=Sheet.ID,Relatedness=Relatedness),mean))
# datstr <- datstr[,which(names(datstr)!="x")]
# datstr <- datstr[order(datstr$Sheet.ID),]

# dat.ag[dat.ag$settler.number<=3,] # Don't use Sheet.ID 18 because settlement was too low
# range(dat.ag[dat.ag$settler.number>3,which(names(dat.ag)%in%c("settler.number"))]) # 4 to 20 settlers per sheet
# with(datstr[datstr$Sheet.ID!=18,],table(Relatedness)) 


# Prepare indices
X.col <- which(names(dat.all) %in% c("True.X"))
Y.col <- which(names(dat.all) %in% c("True.Y"))
Sheet.ID.col <- which(names(dat.all) %in% c("Sheet.ID"))
IDcols.to.use <- c("Sheet.ID","Relatedness")
Lquad <- 61 # length of one side of quadrat in mm, to calculate realized settlement area

# Prepare data frames
tmp1 <- dat.all[!duplicated(dat.all$Sheet.ID),which(names(dat.all) %in% IDcols.to.use)]
dat.ag <- cbind.data.frame(tmp1,matrix(NA,nrow=nrow(tmp1),ncol=13)) # ncol=13 because there are 13 values outputted by ClarkEvans()  
sheet.ind <- unique(dat.all$Sheet.ID)

# Calculate Clark Evans index within plates 
# R = 1 = random. R < 1 = aggregated. R > 1 = regular
for(i in 1:nrow(dat.ag)){
	tmp.dat <- dat.all[dat.all$Sheet.ID==unique(tmp1$Sheet.ID)[i],]
	hold <- ClarkEvans(coords=tmp.dat[,X.col:Y.col],quadratsize=Lquad,boundary.strip=FALSE) 
	dat.ag[i,] <- cbind(tmp1[i,],hold)	
}
names(dat.ag) <- c(names(tmp1),names(hold))



# Set the minimum number of settlers to use for downstream analyses (greater than this number)
# From paper: "We only analyzed spatial patterns on settlement plates with seven or more settlers to avoid bias from low n (Donnelly 1978)."
minsett <- 6


# Statements in Results section
m1 <- lm(R.clarkevans ~ Relatedness, data=dat.ag[dat.ag$settler.number>minsett,])
# par(mfrow=c(2,2));plot(m1)
anova(m1)
# there was no difference in the spatial arrangement of settlers between groups of sibling larvae and groups of non-sibling larvae (F1,21=5.5e-6, p=0.99)

m1 <- lm(min.NN ~ 1, data=dat.ag[dat.ag$settler.number> minsett,])
round(coef(m1),2);round(confint(m1),2)
# The average minimum distance to the nearest neighbor on each plate averaged 4.49 mm (95% CI: 3.4 â€“ 5.58), 

m1 <- lm(min.NN ~ Relatedness, data=dat.ag[dat.ag$settler.number> minsett,])
anova(m1)
# and did not vary between the sibling and non-sibling treatments (F1,21=0.37, p=0.55) (Figure 4c,d).




# Make Figure 4 
cex.size <- 2
cex.lab <- 1.2

cols <- c("#009E73","#0072B2","#D55E00") # green, blue, red  (color-blind friendly)


quartz(width=6,height=5)
nf <- layout(matrix(1:8,2,4,byrow=T),c(3,1,3,1),c(4,4,4,4))
# layout.show(nf)

# a)
par(mar=c(4,1,1,0),oma=c(1,6,2,1))
plot(c(0,0.6),c(0.6,1.7),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=seq(0,0.6,0.1),cex.axis=1.8)
axis(side=2,at=seq(0,2,0.2),cex.axis=1.8,las=1)
abline(h=1,lty=2)
with(dat.ag[dat.ag$Relatedness=="Sib" & dat.ag$settler.number>minsett,],points(settler.density*100, R.clarkevans, 
	pch=19,
	# col=ifelse(p.value.agg.clarkevans<0.05,adjustcolor(cols[3],alpha.f=0.6),
			# ifelse(p.value.reg.clarkevans<0.05,	adjustcolor(cols[1],alpha.f=0.6),	
			# adjustcolor(cols[2],alpha.f=0.6))
			# ),
	col=ifelse(p.value.agg.clarkevans<0.05,cols[3],
			ifelse(p.value.reg.clarkevans<0.05,	cols[1],	
			cols[2])
			),
	cex=cex.size))
mtext("a) Sibs",adj=0,line=0,side=3,cex=cex.lab)

mtext(expression(paste("Index of aggregation (",italic("R"),")",sep="")),adj=0.4,side=2,line=4,cex=cex.lab)

par(mar=c(4,0,1,0))
with(dat.ag[dat.ag$Relatedness=="Sib" & dat.ag$settler.number>minsett,],boxplot(R.clarkevans,axes=F,ylim=c(0.6,1.7)))


# b)
par(mar=c(4,1,1,0))
plot(c(0,0.6),c(0.6,1.7),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=seq(0,0.6,0.1),cex.axis=1.8)
axis(side=2,at=seq(0,2,0.2),labels=NA)
abline(h=1,lty=2)
with(dat.ag[dat.ag$Relatedness=="NonSibs" & dat.ag$settler.number>minsett,],points(settler.density*100, R.clarkevans, 
	pch=19,
	# col=ifelse(p.value.agg.clarkevans<0.05,adjustcolor(cols[3],alpha.f=0.6),
			# ifelse(p.value.reg.clarkevans<0.05,	adjustcolor(cols[1],alpha.f=0.6),	
			# adjustcolor(cols[2],alpha.f=0.6))
			# ),
	col=ifelse(p.value.agg.clarkevans<0.05,cols[3],
			ifelse(p.value.reg.clarkevans<0.05,	cols[1],	
			cols[2])
			),
	cex=cex.size))
mtext("b) Non-sibs",adj=0,line=0,side=3,cex=cex.lab)

legend(0.4,1.8,legend=c("Regular","Random","Aggregrated"),
		bty="n",
		xpd=NA,
		pt.cex=1.7,
		cex=1.5,
		pch=19,
		# col=c(adjustcolor(cols[1],alpha.f=0.6),
			# adjustcolor(cols[2],alpha.f=0.6),
			# adjustcolor(cols[3],alpha.f=0.6)),
		col=c(cols[1],
			cols[2],
			cols[3]),
			)

par(mar=c(4,0,1,0))
with(dat.ag[dat.ag$Relatedness=="NonSibs" & dat.ag$settler.number>minsett,],boxplot(R.clarkevans,axes=F,ylim=c(0.6,1.7)))

# c)
par(mar=c(4,1,1,0))
plot(c(0,0.6),c(0,15),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=seq(0,0.6,0.1),cex.axis=1.8)
axis(side=2,at=seq(0,15,2),cex.axis=1.8,las=1)
with(dat.ag[dat.ag$Relatedness=="Sib" & dat.ag$settler.number>minsett,],points(settler.density*100, min.NN, 
	pch=19,
	# col=ifelse(p.value.agg.clarkevans<0.05,adjustcolor(cols[3],alpha.f=0.6),
			# ifelse(p.value.reg.clarkevans<0.05,	adjustcolor(cols[1],alpha.f=0.6),	
			# adjustcolor(cols[2],alpha.f=0.6))
			# ),
	col=ifelse(p.value.agg.clarkevans<0.05,cols[3],
			ifelse(p.value.reg.clarkevans<0.05,	cols[1],	
			cols[2])
			),
	cex=cex.size))
mtext("c) Sibs",adj=0,line=0,side=3,cex=cex.lab)

mtext("Average minimum distance\\n to nearest neighbor (cm)",side=2,line=3,cex=cex.lab)

par(mar=c(4,0,1,0))
with(dat.ag[dat.ag$Relatedness=="Sib" & dat.ag$settler.number>minsett,],boxplot(min.NN,axes=F,ylim=c(0,15)))

# d)
par(mar=c(4,1,1,0))
plot(c(0,0.6),c(0,15),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="l")
axis(side=1,at=seq(0,0.6,0.1),cex.axis=1.8)
axis(side=2,at=seq(0,15,2),,labels=NA)
with(dat.ag[dat.ag$Relatedness=="NonSibs" & dat.ag$settler.number>minsett,],points(settler.density*100, min.NN, 
	pch=19,
	# col=ifelse(p.value.agg.clarkevans<0.05,adjustcolor(cols[3],alpha.f=0.6),
			# ifelse(p.value.reg.clarkevans<0.05,	adjustcolor(cols[1],alpha.f=0.6),	
			# adjustcolor(cols[2],alpha.f=0.6))
			# ),
	col=ifelse(p.value.agg.clarkevans<0.05,cols[3],
			ifelse(p.value.reg.clarkevans<0.05,	cols[1],	
			cols[2])
			),
	cex=cex.size))
mtext("d) Non-sibs",adj=0,line=0,side=3,cex=cex.lab)

par(mar=c(4,0,1,0))
with(dat.ag[dat.ag$Relatedness=="NonSibs" & dat.ag$settler.number>minsett,],boxplot(min.NN,axes=F,ylim=c(0,15)))

mtext(expression(paste("Number of settlers (",cm^-2,")",sep="")),
	side=1,line=0,cex=cex.lab,outer=T)

