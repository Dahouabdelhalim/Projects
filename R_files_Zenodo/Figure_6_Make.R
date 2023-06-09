## R code to produce Figure 6, and the accompanying analyses presented in the text, in: 
## Burgess SC, Johnston EC, Wyatt ASJ, Leichter JJ, Edmunds PJ (in press) Response diversity in corals: hidden differences in bleaching mortality among cryptic Pocillopora species. Ecology.
## Code written by Scott Burgess. Feb 2021. Send comments or corrections to sburgess@bio.fsu.edu
## R version 3.6.3 (2020-02-29)


library(tidyr)
library(vegan)
library(BiodiversityR)

# Import data
# setwd("")
dat <- read.csv("Data on Haplotypes 2019.csv")

# Prepare data
dat$Trip.Site <- interaction(dat$Trip,dat$Site)

# Function to calculate proportion of haplotypes or species in each Trip and Site
Prop.function <- function(x,ind,type){
	foo <- x[x$Trip.Site==ind,which(names(x)==type)]
	foo1 <- as.data.frame(table(foo))
	foo2 <- as.data.frame(prop.table(table(foo)))
	foo1$Prop <- foo2$Freq
	foo1$Total <- rep(sum(foo1$Freq),(dim(foo1)[1]))
	foo1$Trip.Site <- rep(ind,dim(foo1)[1])
	return(foo1)
}

ind.vec <- unique(dat$Trip.Site)
haplotype.vec <- unique(dat$Species)

dat.prop.haplotype <- list()
for(i in 1:length(ind.vec)){
	foo <- Prop.function(x=dat,ind=ind.vec[i],type="Species")
	dat.prop.haplotype <- rbind(dat.prop.haplotype,foo)		
}
ss <- strsplit(as.character(dat.prop.haplotype$Trip.Site),".",fixed=T)
dat.prop.haplotype$Trip <- sapply(ss,function(x) x[1])
dat.prop.haplotype$Site <- sapply(ss,function(x) x[2])
names(dat.prop.haplotype)[1] <- "Haplotype"
dat.prop.haplotype$Haplotype <- factor(dat.prop.haplotype$Haplotype, levels=haplotype.vec)
# dat.prop.haplotype[1:14,]

dmulti <- spread(dat.prop.haplotype[,c(1,3,5:7)],Haplotype,Prop)
dmulti <- dmulti[order(dmulti$Trip),]
dmat <- dmulti[,-1:-3]

# Get Bray-Curtis distance matrix
dismat <- vegdist(dmat,method="bray")
# Try PCoA
pcoa1 <- cmdscale(dismat,eig=T)
pcoa1$points[,2] <- pcoa1$points[,2]*-1
pcoa1 <- add.spec.scores(pcoa1,dmat,method="pcoa.scores",Rscale=T,scaling=1,multi=1)

# Plot - the analysis uses ALL haplotypes, but here we're only plotting species scores for the common haplotypes
# Extract Site scores (weighted sums of species scores)
xy <- as.matrix(pcoa1$points)
# Extract species scores
haps <- pcoa1$cproj
haps <- haps[c(5,1,4,7,6,3,2),]
hap.names <- c(
	"P. eydouxi",
	"P. meandrina",
	"P. cf. effusus",
	"P. verrucosa",
	"Haplotype 8a",
	"Haplotype 10",
	"Haplotype 11")

haplotype.vec <- c("P. eydouxi",
	"P. meandrina",
	"P. cf. effusus",
	"P. verrucosa",
	"Haplotype 8a",
	"Haplotype 10",
	"Haplotype 11")

cols <- c("#56B4E9", # P.e
	"#0072B2", # P.m 
	"#B2DF8A", # P. cf. e
	"#E69F00", # P.v
	"#CC79A7", # 8a
	"#D55E00", # 10
	"#009E73" # 11
)
	
Haplotype.Colors <- as.data.frame(cbind(haplotype.vec, cols))


name.cols<-as.character(Haplotype.Colors[match(hap.names,Haplotype.Colors$haplotype.vec),2])

# Custom dataframe sort
dat.prop.haplotype$Haplotype <- factor(dat.prop.haplotype$Haplotype, levels=hap.names)
dat.prop.haplotype <- dat.prop.haplotype[order(dat.prop.haplotype$Haplotype),]


## Make Figure 6
quartz(width=7,height=5)
# nf <- layout(matrix(1:15,5,3,byrow=T),c(3,3,3),c(4,1,1,1,1))
par(oma=c(1,1,1.5,1), mar=c(4,4,0,0),fig=c(0,0.7,0,1))

plot(xy,xlim=c(-1,0.9),ylim=c(-0.6,0.8),xaxt="n",yaxt="n",type="n",cex.lab=1.2,ylab="",xlab="")
abline(v=0,h=0,lty=2,col="grey")
text(haps, hap.names,cex=1,col=name.cols)
# Add location of each Site 
pchs <- c(0,1,5,2,15,19,23,8,17,4)
points(xy,pch=pchs,cex=2,col=adjustcolor("black",alpha.f=0.6),bg=adjustcolor("black",alpha.f=0.6),lwd=2)
axis(side=1,at=seq(-1,1,0.2),cex.axis=1.2)
axis(side=2,at=seq(-1,1,0.2),las=2,cex.axis=1.2)
mtext("PCoA 1",side=1,line=2.5,cex=1.2)
mtext("PCoA 2",side=2,line=3.5,cex=1.2)
legend(-0.9,0.8,legend=c(1,2,3,4,5,6),bty="n")
points(rep(-0.85,4),c(0.72,0.64,0.56,0.4),pch=pchs[1:4],cex=1.2,col=adjustcolor("black",alpha.f=0.6),bg=adjustcolor("black",alpha.f=0.6),lwd=2)
points(rep(-0.6,6),c(0.72,0.64,0.56,0.48,0.4,0.32),pch=pchs[5:10],cex=1.2,col=adjustcolor("black",alpha.f=0.6),bg=adjustcolor("black",alpha.f=0.6),lwd=2)
text(-0.71,0.8,paste("Feb","Site","Aug"))

mtext("a)",side=3,adj=0,cex=1.5)
mtext("b)",side=3,adj=1.1,cex=1.5)

par(mar=c(0,4,1,0),fig=c(0.7,1,0.6,1),new=T)
foo1 <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="1.1",which(names(dat.prop.haplotype)=="Prop")]
foo1names <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="1.1",which(names(dat.prop.haplotype)=="Haplotype")]
pie.col<-as.character(Haplotype.Colors[match(foo1names,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(foo1),xlim=c(0,10),las=1,cex.axis=1,col=pie.col,width=1.5)
axis(side=1,at=0.8,line=-1,lty=0,label="Feb",hadj=0.6)
mtext(side=3,"Site 1",adj=0.1)
mtext(side=2,"Proportion of colonies",line=2.5,adj=70)

foo1 <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="2.1",which(names(dat.prop.haplotype)=="Prop")]
foo1names <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="2.1",which(names(dat.prop.haplotype)=="Haplotype")]
pie.col<-as.character(Haplotype.Colors[match(foo1names,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(foo1),xlim=c(0,10),las=1,yaxt="n",col=pie.col,add=T,space=1.5,width=1.5)
axis(side=1,at=3.5,line=-1,lty=0,label="Aug")

foo1 <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="1.2",which(names(dat.prop.haplotype)=="Prop")]
foo1names <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="1.2",which(names(dat.prop.haplotype)=="Haplotype")]
pie.col<-as.character(Haplotype.Colors[match(foo1names,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(foo1),xlim=c(0,10),las=1,yaxt="n",col=pie.col,add=T,space=4,width=1.5)
axis(side=1,at=7,line=-1,lty=0,label="Feb",hadj=0.6)
mtext(side=3,"Site 2",adj=0.95)

foo1 <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="2.2",which(names(dat.prop.haplotype)=="Prop")]
foo1names <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="2.2",which(names(dat.prop.haplotype)=="Haplotype")]
pie.col<-as.character(Haplotype.Colors[match(foo1names,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(foo1),xlim=c(0,10),las=1,yaxt="n",col=pie.col,add=T,space=5.5,width=1.5)
axis(side=1,at=9.5,line=-1,lty=0,label="Aug")

par(mar=c(0,4,1,0),fig=c(0.7,1,0.1,0.5),new=T)
foo1 <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="1.3",which(names(dat.prop.haplotype)=="Prop")]
foo1names <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="1.3",which(names(dat.prop.haplotype)=="Haplotype")]
pie.col<-as.character(Haplotype.Colors[match(foo1names,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(foo1),xlim=c(0,10),las=1,cex.axis=1,col=pie.col,width=1.5)
axis(side=1,at=0.8,line=-1,lty=0,label="Feb",hadj=0.6)
mtext(side=3,"Site 3",adj=0.1)

foo1 <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="2.3",which(names(dat.prop.haplotype)=="Prop")]
foo1names <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="2.3",which(names(dat.prop.haplotype)=="Haplotype")]
pie.col<-as.character(Haplotype.Colors[match(foo1names,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(foo1),xlim=c(0,10),las=1,yaxt="n",col=pie.col,add=T,space=1.5,width=1.5)
axis(side=1,at=3.5,line=-1,lty=0,label="Aug")

foo1 <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="1.5",which(names(dat.prop.haplotype)=="Prop")]
foo1names <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="1.5",which(names(dat.prop.haplotype)=="Haplotype")]
pie.col<-as.character(Haplotype.Colors[match(foo1names,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(foo1),xlim=c(0,10),las=1,yaxt="n",col=pie.col,add=T,space=4,width=1.5)
axis(side=1,at=7,line=-1,lty=0,label="Feb",hadj=0.6)
mtext(side=3,"Site 5",adj=0.95)

foo1 <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="2.5",which(names(dat.prop.haplotype)=="Prop")]
foo1names <- dat.prop.haplotype[dat.prop.haplotype$Trip.Site=="2.5",which(names(dat.prop.haplotype)=="Haplotype")]
pie.col<-as.character(Haplotype.Colors[match(foo1names,Haplotype.Colors$haplotype.vec),2])
barplot(as.matrix(foo1),xlim=c(0,10),las=1,yaxt="n",col=pie.col,add=T,space=5.5,width=1.5)
axis(side=1,at=9.5,line=-1,lty=0,label="Aug")





