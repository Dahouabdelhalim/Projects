#pykF fitness competiitons -- fitness effects across backgrounds
library(lme4)
library(nlme)
library(lmerTest) #p values for fixed effects in lmer
library(afex) #p-values for all effects in lmer
library(plyr)
library(car)

#input data
folder = "~/Dropbox/R/pykF_background/"
d = read.table(paste(folder,"pykF_ev_allele_transfer.txt", sep=""), header=T, skip=4)

#estimate fitness (r)
m.ref = log(d$ref.final*100^d$competition.days/d$ref.initial)
m.test = log(d$test.final*100^d$competition.days/d$test.initial)
r = m.test/m.ref
data = cbind(d,r)

#account for marker cost in 20K A-2 and A-4 competitions
marker.control.A2 = mean(data[data$reference.ID=="TC1512",ncol(data)]) - 1
data[which(data$reference.ID=="TC2024"),ncol(data)] <- data[which(data$reference.ID=="TC2024"),ncol(data)] + marker.control.A2
marker.control.A4 = mean(data[data$reference.ID=="TC1514",ncol(data)]) - 1
data[which(data$reference.ID=="TC2026"),ncol(data)] <- data[which(data$reference.ID=="TC2026"),ncol(data)] + marker.control.A4
#and remove marker control data
data=data[-which(data$reference.ID=="TC1512"),]
data=data[-which(data$reference.ID=="TC1514"),]

#SEM function
SEM=function(d){
	return(sd(d)/sqrt(length(d)))
}

#95CI function
CI95=function(d){
	SE=sd(d)/sqrt(length(d))
	return(SE*qt(0.975,df=length(d)-1))
}

#error bar function
error.bars=function(t,bg){
	offset=ifelse(bg=="Original",0,ifelse(bg=="Del",1,2))
	#return mean-CI and mean+CI
	return(c(t[1+offset]-t[4+offset]/2,t[1+offset]+t[4+offset]/2))	
}
	
#summary table of fitness values per genotype and time
backgrounds = unique(data$reference.ID)
populations = unique(data$population)

alleles = c("Original","Del","A301S")
fit.summary=matrix(NA, nrow=length(backgrounds), ncol=6)
for(j in 1:length(alleles)){
	for(i in 1: length(backgrounds)){
		fit.summary[i,j] = mean(data[data$reference.ID==backgrounds[i] & data$test.genotype==alleles[j],ncol(data)])
		fit.summary[i,j+3] = CI95(d=data[data$reference.ID==backgrounds[i] & data$test.genotype==alleles[j],ncol(data)])
	}
}
rownames(fit.summary)=populations
colnames(fit.summary)=c("original.r","Del.r","A301S.r","original.CI","Del.CI","A301S.CI")
#drop 'original' from ancestor (control competition)
fit.summary[1,c(1,4)]<-NA

#########################################
#ANOVA
#remove ancestor competition
data.NoAnc = data[-which(data$reference.ID=="REL607"),]
#remove 'original' test.genotype (which is not balanced)
data.NoOriginal = data[-which(data$test.genotype=="Original"),]
#removed 'original' mutation and ancestral genotype
data.NoOrAnc = data.NoOriginal[-which(data.NoOriginal$reference.ID=="REL607"),]
#without original or ancestor
m=aov(r ~ population*test.genotype,data=data.NoOrAnc)
summary(m)
#type II makes no difference on balanced data
mII = Anova(m, type=2)
summary(mII)
#without original or ancestor
m1=aov(r ~ population + test.genotype,data=data.NoOrAnc)
summary(m1)
#variance explained by each term:
temp=round(100*(summary(m1)[[1]][2][[1]] / sum(summary(m1)[[1]][2][[1]])),2)
writeLines(paste("Percent variance explained:\\nPopulation: ",temp[1],"\\nAllele: ",temp[2],"\\nPopulation*Allele: ",temp[3],sep=""))

#same, but treating population as a random effect
m1.r.full = lmer(r ~ test.genotype + (1|population) + (1|population:test.genotype), data = data.NoOrAnc, REML=F)
m1.r.noInt = lmer(r ~ test.genotype + (1 | population), data = data.NoOrAnc, REML=F)
m1.r.min = lm(r ~ test.genotype, data = data.NoOrAnc)
m1.r.min2 = lmer(r ~ 1+ (1|population), data = data.NoOrAnc)
summary(m1.r.full)
summary(m1.r.noInt)
anova(m1.r.noInt)

anova(m1.r.full,m1.r.noInt)
anova(m1.r.full,m1.r.min)
anova(m1.r.full,m1.r.min2)

anova(m1.r.noInt,m1.r.min)
anova(m1.r.noInt,m1.r.min2)

#using afex to give p-values
mixed(r ~ test.genotype + (1|population), data = data.NoOrAnc)

#or, paired t-test comparing background effect across two mutations
means=ddply(data.NoOrAnc, c("test.genotype","reference.ID"),function(df) mean(df$r))
t.test(means[1:11,3], means[12:22,3],paired=T)

#loop through each background and test for an allele effect (here up to 3 alleles can be considered)
strains=unique(data$reference.genotype)
genotype.tests=NULL
for(i in strains){
	d.sub = data[data$reference.genotype==i,]
	m = summary(aov(r~test.genotype,d.sub))
	p = m[[1]][[5]][1]
	genotype.tests = rbind(genotype.tests,c(i,p))
}
#########################################
#plots
#Comparison of pykF mutation effect in each strain
xlabels = c("Ancestor", "Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara-5", "Ara-6", "Ara+1", "Ara+2","Ara+3","Ara+5","Ara+6")
fig<-c("pykF_ev_allele_transfer.pdf")
pdf(paste(folder,fig,sep=""),onefile=F,width=4,height=3.5)

par(oma=c(0,0,0,0))
par(mar=c(4,4,0.5,0.5))
offset=0.15
plot(fit.summary[,1],
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.95,1.4),xlim=c(0.5,12.5),
	col="White")
#error bars
for(j in 1:length(alleles)){	
	for(i in 1:nrow(fit.summary)){	
		lines(c(i+(j*offset),i+(j*offset)),error.bars(t=fit.summary[i,], bg=alleles[j]))	
	}
}
#points
points(1:nrow(fit.summary)+1*offset,fit.summary[,1], pch=21, 
	col="Red", lwd=2,
	bg="White"
)
points(1:nrow(fit.summary)+2*offset,fit.summary[,2], pch=22, 
	col="Blue", lwd=2,
	bg=c("White","Red","White","White","White","White","White","White","White","White","White","White")
)
points(1:nrow(fit.summary)+3*offset,fit.summary[,3], pch=23, 
	col="Dark Green", lwd=2,
	bg=c("White","White","White","White","White","Red","White","Red","White","White","Red","White")
)

axis(1,at=(1:nrow(fit.summary)+offset),label=F)
axis(2,at=c(1,1.15,1.3),label=F)
text((1:nrow(fit.summary)+offset), par("usr")[3]-0.03,xlabels, xpd=T, srt=45, cex=0.8, adj=c(1,0))	
mtext(c(0,0.15,0.3),2, at=c(1,1.15,1.3), line=0.65, cex=0.85)
mtext("Population",1, line=2.5)
mtext(expression(Effect~of~italic(pykF)~mutation~(s)),2, line=1.8)

legend("topleft",pch=c(21,22,23), col=c("Red","Blue","Dark green"), bg="White",bty='n',title=expression(italic(pykF)~allele),legend=c("Original","Indel","A301S"), cex=0.75)

dev.off()
