## R code to produce Figure 1 in: 
## Johnston EC, Wyatt ASJ, Leichter JJ, Burgess SC. Niche differences in co-occuring cryptic coral species (Pocillopora spp). Coral Reefs.
## Code written by Scott Burgess. March 2021. Send comments or corrections to sburgess@bio.fsu.edu
## R version 3.6.3 (2020-02-29)

# Import data on temperature and light
PO_data <- read.csv("Figure 1 and 4 Physical Data.csv")

sind <- c(1,2,4,5)
cols <- c("grey80","grey55","grey30","black")
SCnames <- c(
	"Max Daily Temp",
	"Min Daily Temp",
	"Temp Variance\\n(Mean)",
	"Temp Variance\\n(Max)",
	"Mean Light",
	"Min Light")


dev.new(width=6,height=5)
par(mfrow=c(2,3),mar=c(2,4,2,1),oma=c(3,2,0,0))

# A
plot(c(4,21),c(28.5,30.3),type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="l")
for(i in 1:length(sind)){
	with(PO_data[PO_data $Site==sind[i],],lines(as.numeric(as.character(Depth)), Temp.MaxDailyMean.JunNov,lwd=3,lty=1,col=cols[i]))
}
for(i in 1:length(sind)){
	with(PO_data[PO_data $Site==sind[i],],lines(as.numeric(as.character(Depth)), Temp.MaxDailyMean.DecMay,lwd=3,lty=1,col=cols[i]))
}
axis(side=1,at=c(5,10,20),cex.axis=1.2)
axis(side=2,at=seq(25,32,0.1),cex.axis=1.2,las=1)
mtext(side=2,line=3.5,expression(paste(degree,"C")))
mtext(side=3,eval(paste("b) ",SCnames[1],sep="")),cex=1)
text("Dec-May",x=15,y=29.7,cex=1.2)
text("Jun-Nov",x=15,y=29,cex=1.2)

# B
plot(c(4,21),c(26.2,27.3),type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="l")
for(i in 1:length(sind)){
	with(PO_data[PO_data$Site==sind[i],],lines(as.numeric(as.character(Depth)), Temp.MinDailyMean.DecMay,lwd=3,lty=1,col=cols[i]))
}
axis(side=1,at=c(5,10,20),cex.axis=1.2)
axis(side=2,at=seq(25,32,0.2),cex.axis=1.2,las=1)
mtext(side=2,line=3.5,expression(paste(degree,"C")))
mtext(side=3,eval(paste("c) ",SCnames[2],sep="")),bty="n",cex=1)

# C
plot(c(4,21),c(0.007,0.04),type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="l")
for(i in 1:length(sind)){
	with(PO_data[PO_data$Site==sind[i],],lines(as.numeric(as.character(Depth)), Temp.MeanDailyVariance.DecMay,lwd=3,lty=1,col=cols[i]))
}
axis(side=1,at=c(5,10,20),cex.axis=1.2)
axis(side=2,at=seq(0,0.1,0.005),cex.axis=1.2,las=1)
mtext(side=2,line=3.5,expression(paste(degree,C^2)))
mtext(side=3,adj=0.1,padj=0.6,eval(paste("d) ",SCnames[3],sep="")),bty="n",cex=1)

# D
plot(c(4,21),c(0.1,1.1),type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="l")
for(i in 1:length(sind)){
	with(PO_data[PO_data$Site==sind[i],],lines(as.numeric(as.character(Depth)), Temp.MaxDailyVariance.DecMay,lwd=3,lty=1,col=cols[i]))
}
axis(side=1,at=c(5,10,20),cex.axis=1.2)
axis(side=2,at=seq(0,1.5,0.2),cex.axis=1.2,las=1)
mtext(side=2,line=3.5,expression(paste(degree,C^2)))
mtext(side=3,adj=0.1,padj=0.6,eval(paste("e) ",SCnames[4],sep="")),bty="n",cex=1)

# E
plot(c(4,21),c(15,40),type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="l")
for(i in 1:length(sind)){
	with(PO_data[PO_data$Site==sind[i],],lines(as.numeric(as.character(Depth)), Light.Mean.Em2d.DecFeb,lwd=3,lty=1,col=cols[i]))
}
axis(side=1,at=c(5,10,20),cex.axis=1.2)
axis(side=2,at=seq(15,40,5),cex.axis=1.2,las=1)
mtext(side=2,line=3,expression(paste("E",m^2,d^-1)))
mtext(side=3,adj=0,eval(paste("f) ",SCnames[5],sep="")),bty="n",cex=1)

# F
plot(c(4,21),c(0.2,1.7),type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="l")
for(i in 1:length(sind)){
	with(PO_data[PO_data$Site==sind[i],],lines(as.numeric(as.character(Depth)), Light.Min.Em2d.DecFeb,lwd=3,lty=1,col=cols[i]))
}
axis(side=1,at=c(5,10,20),cex.axis=1.2)
axis(side=2,at=seq(0,2,0.1),cex.axis=1.2,las=1)
mtext(side=2,line=3,expression(paste("E",m^2,d^-1)))
mtext(side=3,adj=0,eval(paste("g) ",SCnames[6],sep="")),bty="n",cex=1)

legend('topright',legend=c("Site 1", "Site 2", "Site 4", "Site 5"), col=cols,bty="n",lty=1,lwd=3,cex=1.1)

mtext(side=1,line=1,"Depth (m)",outer=T)
