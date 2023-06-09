#Compare fitness of mutations occurring in EV5 (=rtsgp) replay populations in 
#the ancestor and EV5 genotypes
#Fig3

data = read.table("~/Dropbox/R/Wunsche/fig3_anc-ev.txt", header=T)

#convert relative fitness to mutation effects (s)
data[,c(2,4)]=data[,c(2,4)]-1

pdf("~/Dropbox/R/Wunsche/fig3_anc-ev.pdf",width=2.5,height=2.5)

par(mar=c(2,2,0.2,0.2))
plot(0,
	xlim=c(-0.05, 0.07),
	ylim=c(-0.05, 0.07),
	pch=1, col="white",
	xlab="", ylab="", axes=FALSE)
box()	
axis(1,labels = FALSE, at = c(-0.05, 0, 0.05), tcl = -0.25 )
axis(2, labels = FALSE, at = c(-0.05, 0, 0.05), tcl = -0.25)
mtext(c(-0.05, 0, 0.05), side = 1, line = 0.2, at = c(-0.05, 0, 0.05), cex = 0.7)	
mtext(c(-0.05, 0, 0.05), side = 2, line = 0.4, at = c(-0.05, 0, 0.05), cex = 0.7)		
mtext("Mutation fitness in ancestor", 1, line=1, cex=1)	
mtext(expression('Mutation fitness in Ev'^5), 2, line=1, cex=1)

#isocline
abline(0,1, lty=2, col="grey")

#Error bars
for(i in 1:nrow(data)){
	#Y error bars
	lines(c(data[i,2],data[i,2]), 
	c(data[i,4]-data[i,5]/2, data[i,4]+data[i,5]/2),
	col="gray")
	#X error bars
	lines(c(data[i,2]-data[i,3]/2, data[i,2]+data[i,3]/2), 
	c(data[i,4],data[i,4]),
	col="gray")

}

points(data[,2],data[,4], col="gold", pch=ifelse(data[,6]<=0.05,16,1), cex=1.5, lwd=2)

text(-0.01, 0.06, "Positive epistasis", cex=0.75)
text(0.025, -0.045, "Negative epistasis", cex=0.75)

dev.off()