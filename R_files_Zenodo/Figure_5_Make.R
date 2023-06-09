## R code to produce Figure 5, and the accompanying analyses presented in the text, in: 
## Burgess SC, Johnston EC, Wyatt ASJ, Leichter JJ, Edmunds PJ (in press) Response diversity in corals: hidden differences in bleaching mortality among cryptic Pocillopora species. Ecology.
## Code written by Scott Burgess. Feb 2021. Send comments or corrections to sburgess@bio.fsu.edu
## R version 3.6.3 (2020-02-29)

library(plyr)
library(dplyr)
library(tidyr)

# Import data
# setwd("")
dat <- read.csv("Data on Haplotypes 2019.csv")

# Prepare data
dat$Trip.Site <- interaction(dat$Trip,dat$Site)
dat1 <- dat[dat$Trip=="1",] # Data from February 2019
dat2 <- dat[dat$Trip=="2",] # Data from August 2019


# How many corals ID'd genetically in Feb at 10m
length(dat1$Species) # 68 colonies haplotyped
length(dat1[!is.na(dat1$Photos_max.length.cm), which(names(dat1)=="Species")]) # 67 colonies haplotyped and size measured
dplyr::count(dat1,vars=Site) # Site 1 n=7; Site 2 n=30; Site 3 n=19; Site 5 n=12.
length(dat1[!is.na(dat1$Mortality), which(names(dat1)=="Species")]) # 51 colonies haplotyped, size measured, and survival recorded in Aug.

# How many corals ID'd genetically in Aug at 10m
length(dat2$Species) # 394 colonies haplotyped
dplyr::count(dat2,vars=Site) # Site 1 n=65; Site 2 n=69; Site 3 n=42; Site 4 n=83; Site 5 n=68; Site 6 n=67.


# How many corals ID'd genetically in Feb and Aug at 10m?
length(dat$Species) # 462 colonies haplotyped

# How many haplotypes found in Feb and Aug?
dplyr::count(dat,vars=Species.haplotype) # 12 haplotypes total
dplyr::count(dat,vars=Species) # At least 5 nominal species. 

feb01 <- with(dat1,table(Species,Mortality))
feb01 <- feb01[c(5,6,4,7,3,1,2),]; feb01 # 14 colonies died, 12 were haplotype 11

# How many haplotypes found in Feb where size and mortality was measured?
tmp <- dplyr::count(dat1[!is.na(dat1$Mortality),],vars=Species) 
tmp[c(5,6,4,7,3,1,2),]
sum(tmp$n)

# How many haplotypes found in Feb?
feb <- dplyr::count(dat1,vars=Species) # 7 species in Feb
feb <- feb[c(5,6,4,7,3,1,2),]
feb$prop <- feb$n / sum(feb$n)

# How many haplotypes found in Aug?
aug <- dplyr::count(dat2,vars=Species) # 7 species in Aug
aug <- aug[c(5,6,4,7,3,1,2),]
aug <- aug[match(feb$vars,aug$vars),]
aug$prop <- aug$n / sum(aug$n)

# Size and mortality measured
dplyr::count(dat1[!is.na(dat1$Mortality),],vars=Site) # Site 1 n=6; Site 2 n=27; Site 3 n=9; Site 5 n=9.


# Mortality rates for all haplotypes
m <- glm(Mortality ~ Photos_max.length.cm, data=dat1,family="binomial")
m0 <- glm(Mortality ~ 1, data=dat1,family="binomial")
anova(m,m0,test="Chisq")
slope1 <- c(as.numeric(round(coef(m)[2],2)), as.numeric(round(confint(m)[2,],2)))
mean1 <- c(as.numeric(coef(m0)),as.numeric(confint(m0)))
round(plogis(mean1),2) # 
round(exp(slope1),2)
n1 <- length(dat1$Mortality[!is.na(dat1$Mortality)])
svec <- seq(0,max(dat1$Photos_max.length.cm,na.rm=T),0.1)
p <- predict(m,list(Photos_max.length.cm=svec),se.fit=T)
p1 <- data.frame(size=svec,
	mean=p$fit,
	upr = p$fit + 2*p$se.fit,
	lwr = p$fit - 2*p$se.fit)


# Mortality rates for haplotype 11
sp <- "Haplotype 11"
foo <- dat1[dat1$Species.haplotype==sp,]
n11 <- length(foo$Mortality[!is.na(foo$Mortality)])
mb <- glm(Mortality ~ Photos_max.length.cm, data=foo,family="binomial")
mb0 <- glm(Mortality ~ 1, data=foo,family="binomial")
summary(mb0)

# Mortality did not increase with colony size for haplotype 11
anova(mb,mb0,test="Chisq")
slope11 <- c(as.numeric(round(coef(mb)[2],2)), as.numeric(round(confint(mb)[2,],2)))
mean11 <- c(as.numeric(coef(mb0)),as.numeric(confint(mb0)))

# Mortality rate for haplotype 11
round(plogis(mean11),2)

svec <- seq(0,max(dat1$Photos_max.length.cm,na.rm=T),0.1)
p <- predict(mb,list(Photos_max.length.cm=svec),se.fit=T)
p11 <- data.frame(size=svec,
	mean=p$fit,
	upr = p$fit + 2*p$se.fit,
	lwr = p$fit - 2*p$se.fit)



# Species names for plotting
species.vec <- c("Haplotype 1 (P. eydouxi)",
	"Haplotype 1 (P. meandrina)",
	"Haplotype 2 (P. cf. effusus)",
	"Haplotype 3 ('P. verrucosa')",
	"Haplotype 8a",
	"Haplotype 10",
	"Haplotype 11")

# Species names for indexing dataframe
haplotype.vec <- c("P. eydouxi",
	"P. meandrina",
	"P. cf. effusus",
	"P. verrucosa",
	"Haplotype 8a",
	"Haplotype 10",
	"Haplotype 11")

# Colorblind-friendly colours
cols <- c("#56B4E9", # P.e
	"#0072B2", # P.m 
	"#B2DF8A", # P. cf. e
	"#E69F00", # P.v
	"#CC79A7", # 8a
	"#D55E00", # 10
	"#009E73" # 11
)
	
Haplotype.Colors <- as.data.frame(cbind(haplotype.vec, cols))

hap.vec <- unique(dat1[!is.na(dat1$Mortality),which(names(dat1)=="Species")])
hap.vec <- hap.vec[c(5,3,4,7,6,1,2)]  # 





##### Make Figure 5 
quartz(width=11,height=8)
nf <- layout(matrix(1:15,5,3,byrow=T),c(3,3,3),c(4,1,1,1,1))
par(oma=c(4,3,3,0))
# layout.show(nf)

size.vec <- seq(0,60,10)
prob.vec <- seq(0,1,0.2)
xmax <- 65
xmin <- -5
ymax=0.1
brks <- seq(0,xmax,by=5)

par(mar=c(1,3,0,0))
plot(c(xmin,xmax),c(0,1),type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="l")
with(p1, polygon(c(rev(size),size),c(rev(plogis(upr)),plogis(lwr)),col="grey",border=NA))
with(p1, lines(size,plogis(mean),lwd=2))
axis(side=1,at=size.vec,labels=size.vec,cex.axis=1.5)
axis(side=2,at=prob.vec,las=1,cex.axis=1.5)

segments(xmin,plogis(mean1[2]), xmin,plogis(mean1[3]),lwd=2)
points(xmin,plogis(mean1[1]),pch=19,cex=2,col="grey")

for(i in 1:length(hap.vec)){
	foo <- dat1[dat1$Species==hap.vec[i],]
	cols <- as.character(Haplotype.Colors[Haplotype.Colors$haplotype.vec==as.character(hap.vec[i]),2])
	with(foo[foo$Mortality==0,],rug(Photos_max.length.cm,side=1,ticksize=0.05,lwd=4,col=cols))
	with(foo[foo$Mortality==1,],rug(Photos_max.length.cm,side=3,ticksize=0.05,lwd=4,col=cols))
}
dat1[dat1$Photos_max.length.cm>30,]
mtext(eval(paste("a) All mtORF haplotypes (n = ",n1,")",sep="")),side=3,line=1.5,cex=1.2,adj=0)

plot(c(xmin,xmax),c(0,1),type="n",xaxt="n",yaxt="n",ylab="",xlab="",bty="l")
with(p11, polygon(c(rev(size),size),c(rev(plogis(upr)),plogis(lwr)),col="grey",border=NA))
with(p11, lines(size,plogis(mean),lwd=2))
axis(side=1,at=size.vec,labels=size.vec,cex.axis=1.5)
axis(side=2,at=prob.vec,las=1,cex.axis=1.5)

segments(xmin,plogis(mean11[2]), xmin,plogis(mean11[3]),lwd=2)
points(xmin,plogis(mean11[1]),pch=19,cex=2,col="grey")

cols <- as.character(Haplotype.Colors[Haplotype.Colors$haplotype.vec=="Haplotype 11",2])
foo <- dat1[dat1$Species.haplotype=="Haplotype 11",]
with(foo[foo$Mortality==0,],rug(Photos_max.length.cm,side=1,ticksize=0.05,lwd=4,col=cols))
with(foo[foo$Mortality==1,],rug(Photos_max.length.cm,side=3,ticksize=0.05,lwd=4,col=cols))

mtext("Colony size (cm, longest axis)",side=1,line=3,cex=1.2,adj=0.3)
mtext(eval(paste("b) mtORF haplotype 11 (n = ",n11,")",sep="")),side=3,line=1.5,cex=1.2,adj=0)

# Barplot 
par(mar=c(3,0,3,0),xpd=T)
pie.col<-as.character(Haplotype.Colors[match(feb$vars,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(feb$prop),xlim=c(0,10),las=1,cex.axis=1.5,col=pie.col)
barplot(as.matrix(aug$prop),xlim=c(0,10),las=1,cex.axis=1.5,col=pie.col,add=T,space=2)


mtext("Feb\\n2019",side=3,at=0.5,line=1,cex=1.2)
mtext("Aug\\n2019",side=3,at=2.5,line=1,cex=1.2)
mtext("Proportion",side=2,line=4,cex=1.2)
mtext("c) Relative abundance of haplotypes",side=3,line=4.5,cex=1.2,adj=0)
axis(side=1,at=0.7,eval(paste("n=\\n",sum(feb$n),sep="")),cex.axis=1.5,padj=0.8)
axis(side=1,at=2.5,eval(paste("n=\\n",sum(aug$n),sep="")),cex.axis=1.5,padj=0.8)

text(3.3,0.03,eval(paste(species.vec[1])),cex=1.5,adj=0)
segments(3.0,0.03,3.3,0.03)
text(3.3,0.3,eval(paste(species.vec[2])),cex=1.5,adj=0)
segments(3.0,0.31,3.3,0.31)
text(3.3,0.52,eval(paste(species.vec[3])),cex=1.5,adj=0)
segments(3.0,0.55,3.3,0.53)
text(3.3,0.57,eval(paste(species.vec[4])),cex=1.5,adj=0)
segments(3.0,0.57,3.3,0.58)
text(3.3,0.63,eval(paste(species.vec[5])),cex=1.5,adj=0)
segments(3.0,0.63,3.3,0.63)
text(3.3,0.8,eval(paste(species.vec[6])),cex=1.5,adj=0)
segments(3.0,0.8,3.3,0.8)
text(3.3,0.94,eval(paste(species.vec[7])),cex=1.5,adj=0)
segments(3.0,0.95,3.3,0.95)


# Histograms
par(mar=c(2,3,1,0),xpd=T)
foo <- dat1[dat1$Species=="P. meandrina",which(names(dat1)=="Photos_max.length.cm")]
histlim <- round_any(max(foo,na.rm=T),5,f=ceiling)
hist.col<-as.character(Haplotype.Colors[Haplotype.Colors$haplotype.vec=="P. meandrina",2])
hist(foo,breaks=seq(0,histlim,5),freq=F,main="",ylab="",xlab="",col=hist.col,cex.axis=1.5,las=1,ylim=c(0,ymax),xlim=c(xmin,xmax))
legend(60,0.12,legend="Haplotype 1 (P. meandrina)",bty="n",cex=1.5,adj=1)

plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")
plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")

foo <- dat1[dat1$Species=="Haplotype 10",which(names(dat1)=="Photos_max.length.cm")]
histlim <- round_any(max(foo,na.rm=T),5,f=ceiling)
hist.col<-as.character(Haplotype.Colors[Haplotype.Colors$haplotype.vec=="Haplotype 10",2])
hist(foo,breaks=seq(0,histlim,5),freq=F,main="",ylab="",xlab="",col=hist.col,cex.axis=1.5,las=1,ylim=c(0,ymax),xlim=c(xmin,xmax))
legend(60,0.14,legend="Haplotype 10",bty="n",cex=1.5,adj=1)

plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")
plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")

foo <- dat1[dat1$Species=="Haplotype 11",which(names(dat1)=="Photos_max.length.cm")]
histlim <- round_any(max(foo,na.rm=T),5,f=ceiling)
hist.col<-as.character(Haplotype.Colors[Haplotype.Colors$haplotype.vec=="Haplotype 11",2])
hist(foo,breaks=seq(0,histlim,5),freq=F,main="",ylab="",xlab="",col=hist.col,cex.axis=1.5,las=1,ylim=c(0,ymax-0.05),xlim=c(xmin,xmax))
legend(60,0.06,legend="Haplotype 11",bty="n",cex=1.5,adj=1)

plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")
plot(c(0,1),c(0,1),type="n",ylab="",xlab="",xaxt="n",yaxt="n",bty="n")

foo <- dat1[!dat1$Species %in% c("P. meandrina", "Haplotype 10", "Haplotype 11"),which(names(dat1)=="Photos_max.length.cm")]
histlim <- round_any(max(foo,na.rm=T),5,f=ceiling)
hist.col<-"grey"
hist(foo,breaks=seq(0,histlim,5),freq=F,main="",ylab="",xlab="",col=hist.col,cex.axis=1.5,las=1,ylim=c(0,ymax),xlim=c(xmin,xmax))
legend(60,0.1,legend="All other haplotypes ",bty="n",cex=1.5,adj=1)
mtext("Colony size (cm, longest axis)",side=1,line=3,cex=1.2,adj=0.3)

mtext("Probability of bleaching mortality",side=2,line=1,outer=T,cex=1.2,adj=1)
mtext("Proportion of colonies",side=2,line=1,outer=T,cex=1.2,adj=0.2)

