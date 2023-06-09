#Produce initial fitness vs. change in fitness figure
#compare different fits for relationship
#Also draw as sections along a fitness trajectory

library(nlme)

data = read.table("~/Dropbox/R/Wunsche/fitnesschange.tab")

#Mean Founder fitness
data.mean = tapply(data$fitness, data$genotype, mean)-1

#Founder fitness change mean weighted by 1/sd of each descendent population's fitness estimate
data.Wmean = NULL
genotypes = unique(data$genotype)
for(i in 1:length(genotypes)){
	temp = data[data[,1]==genotypes[i],4:5]
	temp2 = weighted.mean(temp[,1],1/temp[,2])
	data.Wmean = c(data.Wmean,temp2)
}

data.Wmean=data.Wmean-1

#initial fitness of founder genotypes
initial.fitness = c(1,1.012,1.103,1.197,1.228,1.328,1.397,1.551)-1

############
#regressions
md1 = lm(data.Wmean ~ initial.fitness)
md1a = lm(data.Wmean ~ 1)
md2 = lm(data.Wmean ~ initial.fitness + I(initial.fitness^2))
anova(md1,md2)

#negative exponential
md3 = nls(data.Wmean ~ C * exp(-k*initial.fitness),   start=c(C=10, k=-1))

#power law -- note offset to ancestor fitness change
md4 = nls(data.Wmean ~ (b*initial.fitness+ (1/data.Wmean[1]))^a,   start=c(b = .15, a=-1))

AIC(md3, md1a, md1, md4)

#Calculate variance explained by regression models (per Bates at: https://stat.ethz.ch/pipermail/r-help/2006-November/118263.html)
#exponential
Rsq.E = 1 - sum(resid(md3)^2)/sum((data.Wmean - mean(data.Wmean))^2)
#power law
Rsq.PL = 1 - sum(resid(md4)^2)/sum((data.Wmean - mean(data.Wmean))^2)

##################################
#Fig. 1b: fitness change of replay populaitons started from different founders
#founder names
f = unique(data$genotype)
color = c("black","purple","blue","dark green","green","gold","orange","red")
pdf("~/Dropbox/R/Wunsche/fitnesschange_all-points_pwr-law.pdf",width=2.5,height=2.5)

par(mar=c(2,2,0.2,0.2))
plot(0,
	xlim=c(1,1.6),ylim=c(-0.05,0.17), pch=1, col="white",
	xlab="", ylab="", axes=FALSE)
box()	
axis(1,labels = FALSE, at = c(1,1.5),tcl = -0.25 )
axis(2, labels = FALSE, at = c(0,0.1), tcl = -0.25)
mtext(c(1,1.5), side = 1, line = 0.2, at = c(1,1.5), cex = 0.7)	
mtext(c(0,0.1), side = 2, line = 0.4, at = c(0,0.1), cex = 0.7)		
mtext("Founder genotype fitness", 1, line=1, cex=1)	
mtext("Change in fitness", 2, line=1, cex=1)

#raw data each founder
for(i in 1:length(f)){
	temp=data[data$genotype==f[i],4]-1
	points(rep(initial.fitness[i]+1,length(temp)),temp,col=color[i])
	
}

#mean fitness each founder
points(initial.fitness+1,data.Wmean, col="black",cex=1.5, pch=16)

#draw regression lines
#exponential
x=seq(0,0.55, length.out=100)
#p=predict(md3,list(initial.fitness=x))
#pred2=data.frame(x,p)[order(x),]
#lines(pred2[,1]+1,pred2[,2],lty=2)

#power law
p4=predict(md4,list(initial.fitness=x))
pred4=data.frame(x,p4)[order(x),]
lines(pred4[,1]+1,pred4[,2],lty=2)

dev.off()

##################################
#Plot initial vs. final fitness
par(mar=c(2,2,0.2,0.2))
plot(0,
	xlim=c(1,1.6),ylim=c(1,1.7), pch=1, col="white",
	xlab="", ylab="", axes=FALSE)
box()	
axis(1,labels = FALSE, at = c(1,1.5),tcl = -0.25 )
axis(2, labels = FALSE, at = c(0,0.1), tcl = -0.25)
mtext(c(1,1.5), side = 1, line = 0.2, at = c(1,1.5), cex = 0.7)	
mtext(c(0,0.1), side = 2, line = 0.4, at = c(0,0.1), cex = 0.7)		
mtext("Founder genotype fitness", 1, line=1, cex=1)	
mtext("Change in fitness", 2, line=1, cex=1)

#raw data each founder
for(i in 1:length(f)){
	temp=data[data$genotype==f[i],4]-1
	points(rep(initial.fitness[i]+1,length(temp)),temp+initial.fitness[i]+1,col=color[i])
	
}
#isocline
abline(0,1, col="black", lty=2)

#mean fitness each founder
points(initial.fitness+1,data.Wmean-1, col="black",cex=1.5, pch=16)

#draw regression lines
#exponential
x=seq(0,0.55, length.out=100)
p=predict(md3,list(initial.fitness=x))
pred2=data.frame(x,p)[order(x),]
lines(pred2[,1]+1,pred2[,2],lty=2)

