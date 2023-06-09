require(adegenet)

###Set working directory
setwd("")

###Define pop as locality, other as region (Alps, Tohoku, Hokkaido)
dataManub<-read.structure("ManubPaperDip.str", n.ind=72, n.loc=99, onerowperind=FALSE, 
col.lab=1, col.pop=2, col.others=3, row.marknames=NULL, NA.char="-9", ask=FALSE, quiet=FALSE)

##Population Stats with Hierfstat, all Manubriatum samples
require(hierfstat)
require(pegas)

mlocalityhier<-genind2hierfstat(dataManub)
mlh<-mlocalityhier[order(mlocalityhier$pop), ]
localstats<-basic.stats(mlocalityhier)

pop<-rep(letters[1:14],table(mlh$pop))
ManubFST<-pairwise.fst(dataManub, res.type=c("matrix"))
ManubFSTCI<-boot.ppfst(mlh)
ManubFSTCI<-as.matrix(ManubFSTCI)

ManubFisCI<-boot.ppfis(mlh)

############Manubriatum nucleotide diversity
############
ManubSeq<-read.dna("ManubPaperDip.phy")
manubdiv<-nuc.div(ManubSeq, pairwise.deletion = TRUE)

ManubAlps<-read.dna("ManubNaga.phy")
ManubTohoku<-read.dna("ManubTohoku.phy")
ManubHokk<-read.dna("ManubHokk.phy")

nuc.div(ManubTohoku, variance=TRUE, pairwise.deletion = FALSE)
nuc.div(ManubHokk, variance=TRUE, pairwise.deletion = FALSE)

###ManubAlps (got NaNs otherwise)
R.b=5000 #number of replicates
pi.values.b <- numeric(R.b) #empty vector to store results
for (i in 1:R.b) {
  subsample <- ManubAlps[sample(nrow(ManubAlps), size=nrow(ManubAlps), replace=T),] #bootstrap does subsampling with replacement
  pi <- nuc.div(subsample,pairwise.deletion=TRUE)
  pi.values.b[i] <- pi
}

mean.pi.b <- mean(pi.values.b,na.rm=TRUE) #mean value for pi
mean.pi.b
hist(pi.values.b) #plot histogram of values
sd(pi.values.b, na.rm=TRUE) #standard deviation
quantile(pi.values.b, c(.025, .975), na.rm=TRUE) #95% confidence interval

mean.pi.b/99
(sd(pi.values.b, na.rm=TRUE))/99

###Manubriatum DAPC
dataManubMt<-read.structure("ManubPaperDipmtDNA.str", n.ind=65, n.loc=99, onerowperind=FALSE, 
col.lab=1, col.pop=2, col.others=NULL, row.marknames=NULL, NA.char="-9", 
pop=NULL, ask=FALSE, quiet=FALSE)

data_scaledManub <- scaleGen(dataManubMt, center=FALSE, scale=FALSE, NA.method=c("mean"))

##Priors for DAPC
ManubGroupings<-read.table("ManubPaperGroupings.txt", row.names=1)
ManubRegion<-ManubGroupings$V2
ManubElevation<-ManubGroupings$V3
ManubMtNet<-ManubGroupings$V4
ManubElevONLY<-ManubGroupings$V5
ManubMtPhylo<-ManubGroupings$V6

things<-attributes(dataManubMt@tab)
ManubNames<-things[[2]][[1]]
ManubLoc<-dataManubMt@pop
ManubLoc<-as.matrix(ManubLoc)
rownames(ManubLoc)<-ManubNames

require(poppr)

###Manubriatum DAPC
gcManub<-as.genclone(dataManubMt)
Mregionx <- xvalDapc(tab(gcManub, NA.method = "mean"),ManubRegion)
Mlocx <- xvalDapc(tab(gcManub, NA.method = "mean"),ManubLoc)
Mbothx <- xvalDapc(tab(gcManub, NA.method = "mean"),ManubElevation)
MElevONLYx <- xvalDapc(tab(gcManub, NA.method = "mean"),ManubElevONLY)
MMtPhylox<-xvalDapc(tab(gcManub, NA.method = "mean"),ManubMtPhylo)

dataforMt<-read.structure("ManubPaperDipmtDNAONLY2.str", n.ind=59, n.loc=99, onerowperind=FALSE, 
col.lab=1, col.pop=2, col.others=NULL, row.marknames=NULL, NA.char="-9", 
pop=NULL, ask=FALSE, quiet=FALSE)
gcManubMts<-as.genclone(dataforMt)
ManubMtNoNAs <- ManubMtNet[!is.na(ManubMtNet)]
MMtxNet <- xvalDapc(tab(gcManubMts, NA.method = "mean"),ManubMtNoNAs)

(summary(Mlocx$DAPC)$assign.prop)
(summary(Mregionx$DAPC)$assign.prop)
(summary(MElevONLYx$DAPC)$assign.prop)
(summary(Mbothx$DAPC)$assign.prop)
(summary(MMtxNet$DAPC)$assign.prop)
(summary(MMtPhylox$DAPC)$assign.prop)

####Globosum population statistics
###Define pop as locality
dataGlobo<-read.structure("GloboPaperDip2.str", n.ind=41, n.loc=626, onerowperind=FALSE, 
col.lab=1, col.pop=2, col.others=3, row.marknames=NULL, NA.char="-9", ask=FALSE, quiet=FALSE)

require(hierfstat)
require(pegas)

glocalityhier<-genind2hierfstat(dataGlobo)
glocalstats<-basic.stats(glocalityhier)
GloboFST<-pairwise.fst(dataGlobo, res.type=c("matrix"))
glocalityhier<-genind2hierfstat(dataGlobo,pop=NULL)
gregionhier<-genind2hierfstat(dataGlobo,pop=GloboRegion)
glocstats<-basic.stats(glocalityhier)
pairwise.fst(dataGlobo, pop=GloboRegion, res.type=c("matrix"))
boot.ppfst(dat=gregionhier,quant=c(0.025,0.975),diploid=TRUE)
vcglobo<-boot.vc(levels=levels2,loci=gregionhier,diploid=TRUE,nboot=1000,quant=c(0.025,0.5,0.975))
varcomp.glob(levels=levels2, loci=mregionhier, diploid=TRUE)
names(vcmanub)
boot.ppfis(mlocalityhier, nboot=5000) 

###Nucleotide Diversity for Globosum

GloboSeq<-read.dna("GloboPaperDip2.phy")
globodiv<-nuc.div(GloboSeq,pairwise.deletion = TRUE)
globodiv/626

GloboTohoku<-read.dna("GloboTohoku.phy")
GloboHokk<-read.dna("GloboHokk.phy")

###Globosum from Tohoku
R.gt=5000 #number of replicates
pi.values.gt <- numeric(R.gt) #empty vector to store results
for (i in 1:R.gt) {
  subsample <- GloboTohoku[sample(nrow(GloboTohoku), size=nrow(GloboTohoku), replace=T),] #bootstrap does subsampling with replacement
  pi <- nuc.div(subsample,pairwise.deletion=TRUE)
  pi.values.gt[i] <- pi
}

mean.pi.gt <- mean(pi.values.gt, na.rm=TRUE) #mean value for pi
mean.pi.gt
hist(pi.values.gt) #plot histogram of values
sd(pi.values.gt, na.rm=TRUE) #standard deviation
quantile(pi.values.gt, c(.025, .975),na.rm=TRUE) #95% confidence interval

mean.pi.gt/626
(sd(pi.values.gt, na.rm=TRUE))/626

###Globosum from Hokkaido
R.gh=5000 #number of replicates
pi.values.gh <- numeric(R.gh) #empty vector to store results
for (i in 1:R.gh) {
  subsample <- GloboHokk[sample(nrow(GloboHokk), size=nrow(GloboHokk), replace=T),] #bootstrap does subsampling with replacement
  pi <- nuc.div(subsample,pairwise.deletion=TRUE)
  pi.values.gh[i] <- pi
}

mean.pi.gh <- mean(pi.values.gh, na.rm=TRUE) #mean value for pi
mean.pi.gh
hist(pi.values.gh) #plot histogram of values
sd(pi.values.gh, na.rm=TRUE) #standard deviation
quantile(pi.values.gh, c(.025, .975),na.rm=TRUE) #95% confidence interval

mean.pi.gh/626
(sd(pi.values.gh, na.rm=TRUE))/626

####Globosum DAPC
setwd("")
dataGloboMt<-read.structure("GloboPaperDipMtOnly.str", n.ind=40, n.loc=626, onerowperind=FALSE, 
col.lab=1, col.pop=2, col.others=NULL, row.marknames=NULL, NA.char="-9", ask=FALSE, quiet=FALSE)

data_scaledGlobo <- scaleGen(dataGloboMt, center=FALSE, scale=FALSE, NA.method=c("mean"))

##Priors for DAPC
GloboGroupings<-read.table("GloboPaperGroupings.txt", row.names=1)
GloboRegion<-GloboGroupings$V2
GloboElevation<-GloboGroupings$V3
GloboMtNet<-GloboGroupings$V4
GloboElevONLY<-GloboGroupings$V5
GloboMtPhylo<-GloboGroupings$V6

things2<-attributes(dataGloboMt@tab)
GloboNames<-things2[[2]][[1]]
GloboLoc<-dataGloboMt@pop
GloboLoc<-as.matrix(GloboLoc)
rownames(GloboLoc)<-GloboNames

##Globosum DAPC
gcGlobo<-as.genclone(dataGloboMt)
Gregionx <- xvalDapc(tab(gcGlobo, NA.method = "mean"),GloboRegion)
Glocx <- xvalDapc(tab(gcGlobo, NA.method = "mean"),GloboLoc)
Gbothx <- xvalDapc(tab(gcGlobo, NA.method = "mean"),GloboElevation)
GMtxNet <- xvalDapc(tab(gcGlobo, NA.method = "mean"),GloboMtNet)
GElevONLYx <- xvalDapc(tab(gcGlobo, NA.method = "mean"),GloboElevONLY)
GMtxPhylo<- xvalDapc(tab(gcGlobo, NA.method = "mean"),GloboMtPhylo)

(summary(Glocx$DAPC)$assign.prop)
(summary(Gregionx$DAPC)$assign.prop)
(summary(GElevONLYx$DAPC)$assign.prop)
(summary(Gbothx$DAPC)$assign.prop)
(summary(GMtxNet$DAPC)$assign.prop)
(summary(GMtxPhylo$DAPC)$assign.prop)





