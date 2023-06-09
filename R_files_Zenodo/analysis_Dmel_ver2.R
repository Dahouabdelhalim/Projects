####------------------------------------------------------#####
#Analysis of Drosophila melanogaster choice trial data 
#Coded in R
####------------------------------------------------------#####

####------RELATEDNESS TRIALS------####

relatednessDat <- read.csv("relatednessChoiceTrials_Dmel.csv", header=TRUE)

related <- relatednessDat$related
unrelated <- relatednessDat$unrelated

wilcox.test(related,unrelated,paired=TRUE, correct=FALSE)

relatedErr <- sd(related)/sqrt(length(related))*1.96
unrelatedErr <- sd(unrelated)/sqrt(length(unrelated))*1.96

####------FAMILIARITY TRIALS------####

familiarityDat <- read.csv("familiarityChoiceTrials_Dmel.csv", header=TRUE)

familiar <- familiarityDat$familiar
unfamiliar <- familiarityDat$unfamiliar

wilcox.test(familiar,unfamiliar,paired=TRUE)

familiarErr <- sd(familiar)/sqrt(length(familiar))*1.96
unfamiliarErr <- sd(unfamiliar)/sqrt(length(unfamiliar))*1.96

############Relatedness and familiarity plot############

#tiff("kinAndFamiliarityPlot.tiff",height=6,width=7, units='in', res=600)

split.screen(rbind(c(0.08,0.54,0.15,0.9),c(0.54,1,0.15,0.9)))

screen(1)

par(mar=c(0.1,1.5,0.5,1))

plot(c(0.15,0.52), c(-0.2, 6.2), type='n', ylab="Cannibal count", xlab='', xaxt='n', las=1)
	grid(lty=1, lwd=0.6, col='grey')
	axis(1, at = c(0.25, 0.45), labels = c("Kin", "Non-Kin"), cex.axis=1.25)
	points(0.25, mean(related), pch=15, cex=2, col=rgb(30/252, 30/252, 30/252))
	lines(c(0.25, 0.25), c(mean(related)+relatedErr, mean(related)-relatedErr), col=rgb(30/252, 30/252, 30/252), lwd=2.5)
	points(jitter(rep(0.2, length(related))), jitter(related), pch=16, cex=1, col=rgb(30/252, 30/252, 30/252))
	points(0.45, mean(unrelated), pch=15, cex=2, col=rgb(175/252, 175/252, 175/252))
	lines(c(0.45, 0.45), c(mean(unrelated)+unrelatedErr, mean(unrelated)-unrelatedErr), col=rgb(175/252, 175/252, 175/252), lwd=2.5)
	points(jitter(rep(0.4, length(unrelated))), jitter(unrelated), pch=16, cex=1, col=rgb(175/252, 175/252, 175/252))
	mtext("A",side=3,adj=0,line=0.5,cex=1.5)
	lines(c(0.25, 0.45), c(5.5, 5.5))
	#text(0.35, 5.65, labels="*", cex=2)
	text(0.35, 5.67, label=expression(paste(italic("p"), "=0.014")))
	mtext("Cannibal count", side=2, line=2, cex=1.5)
	
screen(2)
par(mar=c(0.1,1.5,0.5,1))

plot(c(0.15,0.52), c(-0.2, 6.2), type='n', ylab="Cannibal count", xlab='', xaxt='n', las=1)
	grid(lty=1, lwd=0.6, col='grey')
	axis(1, at = c(0.25, 0.45), labels = c("Familiar", "Unfamiliar"), cex.axis=1.25)
	points(0.25, mean(familiar), pch=15, cex=2, col=rgb(30/252, 30/252, 30/252))
	lines(c(0.25, 0.25), c(mean(familiar)+familiarErr, mean(familiar)-familiarErr), col=rgb(30/252, 30/252, 30/252), lwd=2.5)
	points(jitter(rep(0.2, length(familiar))), jitter(familiar), pch=16, cex=1, col=rgb(30/252, 30/252, 30/252))
	points(0.45, mean(unfamiliar), pch=15, cex=2, col=rgb(175/252, 175/252, 175/252))
	lines(c(0.45, 0.45), c(mean(unfamiliar)+unfamiliarErr, mean(unfamiliar)-unfamiliarErr), col=rgb(175/252, 175/252, 175/252), lwd=2.5)
	points(jitter(rep(0.4, length(unfamiliar))), jitter(unfamiliar), pch=16, cex=1, col=rgb(175/252, 175/252, 175/252))
	mtext("B",side=3,adj=0,line=0.5,cex=1.5)
	lines(c(0.25, 0.45), c(5.5, 5.5))
	#text(0.35, 5.65, labels="NS", cex=1)
	text(0.35, 5.67, label=expression(paste(italic("p"), "=0.101")))


close.screen(all=T)

#dev.off()

####-----CO-CANNIBALISM TRIALS-----####

library('lmerTest')
library('lme4')

dat <- read.csv("coCannibalismSkew_Dmel.csv", header=TRUE)

oneHourDat <- dat[dat$Time==1, ]

#tendency to cannibalise is not different between the 2 genotypes
wilcox.test(oneHourDat$Total.wt, oneHourDat$Total.GFP, paired=TRUE)

#number of wt cannibals is not significantly different from the number of GFP on a given victim
wtLong <- c(oneHourDat$wt.left, oneHourDat$wt.right)
GFPlong <- c(oneHourDat$GFP.left, oneHourDat$wt.right)

wilcox.test(wtLong, GFPlong, paired=TRUE)

#Victim skew vs time

times <- c(1, 2, 3, 6)
timePoints <- c()
skew <- c()
ID <- c()

for(i in 1:length(times)){
	tempDat <- dat[dat$Time==times[i], ]
	skew <- c(skew, abs(tempDat$Total.left-tempDat$Total.right))
	timePoints <- c(timePoints, tempDat$Time)
	ID <- c(ID, tempDat$ID)
}
	
expectedDist <- rpois(length(skew), mean(skew))

chisq.test(expectedDist, skew) #data follows expected distribution, negative binomial not necessary 

#Random effect for ID to account for pseudoreplication produce by repeated measures	
model <- glmer(skew~poly(timePoints,2)+(1|ID), family='poisson')
model2 <- glmer(skew~timePoints+(1|ID), family='poisson')

anova(model,model2) #Model not significantly better fit than model 2 - use model 2

summary(model2)

#Model used to generate stdError
model <- glm(skew~timePoints, family='poisson')

predFrame <- data.frame(timePoints=seq(1, 6, length=1000))

predictedVals <- predict(model, predFrame, type = 'response', level = 0.95, se.fit = TRUE)


upperSE <- predictedVals$fit + predictedVals$se.fit*1.96
lowerSE <- predictedVals$fit - predictedVals$se.fit*1.96

#tiff("coCannSkew.tiff",width = 4, height = 5.25, units = 'in', res = 600)

par(mar=c(4,4,1.4,1.4))

plot(seq(1,6, length=1000), predictedVals$fit, type='n', ylab="Skew (|side 1 - side 2|)", xlab="Time (hours)", ylim=c(1,3.5), las=1, cex.lab=1.25)
	grid(lty=1, lwd=0.4, col='grey')
	polygon(c(seq(1,6,length=1000),rev(seq(1,6,length=1000))),c(upperSE,rev(lowerSE)),col=rgb(181/250, 171/250, 170/250, 0.4),border=F)
	lines(seq(1,6, length=1000), predictedVals$fit, lwd=3, lty=1)
	lines(seq(1,6, length=1000), upperSE, lty=2)
	lines(seq(1,6, length=1000), lowerSE, lty=2)
	
#dev.off()

