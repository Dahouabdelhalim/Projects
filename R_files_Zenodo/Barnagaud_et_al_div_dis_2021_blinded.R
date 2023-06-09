### R script for Barnagaud et al, Diversity & Distributions, 2021 ###
# Trait-habitat associations explain novel bird assemblages mixing native and alien species across New-Zealand landscapes
# last update : 04 / 10 / 2021
# for any question or comment please contact Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu
###---------------------------------------------------------------###

graph=1
#---------------------#
### useful packages ###
#---------------------#

library(ade4) ; library(adiv)  # for multivariate analyses
library(factoextra) ; library(adegraphics) # graphical representations of dudi objects
library(geiger)  ; library(ape) ; library(picante) ; library(phytools) ; library(phylolm) # for phylogenetic analyses
library(raster)  ; library(spdep) # processing of spatial data
library(RColorBrewer) ; library(visreg) ; library(heplots) ; library(candisc) # for graphics
library(lmtest) # for testing the independence of manova residuals
library(reshape2)  ; library(beepr) # misc

#----------#
### data ###
#----------#
# NOTE: tables L.txt and R.txt are subject to the agreement of Scion Research.  Given the contractual terms imposed by Scion, a simple citation of the article is not a sufficient credit.
# Hence, all species names have been blinded. Please contact the authors before any use.

L = read.table("L_blinded.txt",header=T,sep="\\t",row.names=1) # site * species matrix (counts)
Q = read.table("Q_blinded.txt",header=T,sep="\\t",as.is=F,row.names=1) # traits
R = read.table("R.txt",header=T,sep="\\t") # environmental variables
TBC2 = read.table("introduction_effort_blinded.txt",sep="\\t",header=T,row.names=1) # introduction effort variables
NZtree_P = read.tree("NZBirds_Prum_blinded.tre") # phylogenetic tree
xy0 = read.table("geographic_coordinates.txt",sep="\\t",header=T) # geographic coordinates of points
status=read.table("species_status_blinded.txt",header=T,sep="\\t")
hab = read.table("hab.txt",header=T,sep="\\t") # habitat table
traitbin = read.table("traits_bird_NZ_bin_blinded.txt",
											header=T,sep="\\t",as.is=F) # table Q reformatted
Q_bin = traitbin[,c(2:28)]

#--------------------#
#### process data ####
#--------------------#

# replace NA in TBC2 by the average value of the column
mrel = mean(TBC2$Released_tot,na.rm = T) 
TBC2[is.na(TBC2$Released_tot),"Released_tot"] = mrel

mintro = mean(TBC2$Nb_intro,na.rm = T) # compute mean
TBC2[is.na(TBC2$Nb_intro),"Nb_intro"] = mintro # replace NA

Pdist = cophenetic.phylo(NZtree_P)

# process spatial data
xy = xy0[,c(2,3)] # centroids of sites (order is the same as in R and L)
gneigh = gabrielneigh(as.matrix(xy)) # Gabriel neighbour matrix
nb = graph2nb(gneigh)
neig = nb2neig(nb)
mat_S = neig2mat(neig)
vec_S = scores.neig(neig) # scores used in the PCA (cf. Pavoine et al., 2011, appendix)

#------------------------------------------------------------------#
#### preliminary tests as in Pavoine et al., Journ. Ecol., 2011 ####
#------------------------------------------------------------------#

# check that "functions_Barnagaud_DD_2021.R" are sourced
# ! some functions take a rather long time to run
source("functions_Barnagaud_DD_2021.R")

# test spatial autocorrelation in environmental variables
gm = gearymoran(mat_S, R) 
plot(gm)

# test phylogenetic signal in traits
traitsQ = Q[,c(2,6:8,10)]
traitsN = Q[,c(1,3:5,9)]
distT = dist.ktab(ktab.list.df(list(traitsQ,traitsN)),c("Q","N"))
phy = NZtree_P
phystot = rtest.decdiv(phy, rep(1, 48), as.dist(as.matrix(distT)[names(phy$leaves),
																																 names(phy$leaves)]), nrep = 999, vranking =
											 	"droot", optiontest = "less", ties.method = "average", option = 3)
listdis = ldist.ktab(ktab.list.df(list(traitsQ, traitsN)), c("Q", "N"), scan = TRUE)
1
physTest = rtest.decdiv(phy, rep(1, 48), as.dist(as.matrix(listdis$Social)[names(phy$leaves),
																																					 names(phy$leaves)]), nrep = 999, vranking = "droot",
												optiontest = "less", ties.method = "average", option = 3)

# test phylogenetic clustering vs dispersion
for (i in 1){
	TQE <- TPQE(as.data.frame(t(L)), distT) #; plot(TQE)
	PQE <- TPQE(as.data.frame(t(L)),
							as.dist(as.matrix(phy$Wdist)[names(L), names(L)])) #; plot(PQE)
	beep(sound = "treasure")
}

#--------------------#
#### RLQ analysis ####
#--------------------#

afc = dudi.coa(L, scannf = F, nf = 2) # correspondence analysis on L
hs = dudi.hillsmith(Q, row.w = afc$cw, scannf = F, nf = 2) # Hill & Smith analysis on Q
acp = dudi.pca(R, row.w = afc$lw, scannf = F, nf = 2) # PCA on R

# classical (non corrected) RLQ
RLQ = rlq(acp,afc,hs, scannf = F, nf = 2)


pco = dudi.pco(as.dist(as.matrix(NZtree_Pphylog$Wdist)[names(L), names(L)]), afc$cw, 
							 scannf = F, nf = 2) # principal coordinates analysis on phylogenetic distances

acp_S = dudi.pca(vec_S, row.w = afc$lw, scannf = FALSE, nf = 2) # PCA on spatial matrix

#-------------------------------------------------------------------------------------------------------#
#### RLQ analysis incorporating phylogeny and spatial autocorrelation - ESLTP in Pavoine et al. 2011 ####
#-------------------------------------------------------------------------------------------------------#

ESLTP = rlqESLTP(dudiE = acp, dudiS = acp_S, dudiL = afc, dudiT = hs, dudiP = pco,
								 scannf = FALSE, nf = 2)

(ESLTP$eig)/sum(ESLTP$eig) # % variance explained

# test phylogenetic signal on axes
axesRLQsp = ESLTP$mQ
distT_axes = dist.ktab(ktab.list.df(list(axesRLQsp)),c("Q"))
phy = NZtree_Pphylog
signal = rtest.decdiv(phy, rep(1, 48), as.dist(as.matrix(distT_axes)[names(phy$leaves),
																																		 names(phy$leaves)]), nrep = 999, vranking = "droot",
											optiontest = "less", ties.method = "average", option = 3)
listdis_axes <- ldist.ktab(ktab.list.df(list(axesRLQsp)), c("Q"), scan = TRUE)
1
physTest_axes <- rtest.decdiv(phy, rep(1, 48), as.dist(as.matrix(listdis_axes$NorS1)[names(phy$leaves),
																																										 names(phy$leaves)]), nrep = 999, vranking = "droot",
															optiontest = "less", ties.method = "average", option = 3)

#----------------------------------#
### graphical display : Figure 2 ###
#----------------------------------#

# environmental variables (Fig 2a)
cor = round(cor(ESLTP$mR,R),2)
rownames(cor) = c("Ax1","Ax2")
matrixE = as.matrix(cor)
melted_matrixE = melt(matrixE)
	x11();ggplot(data = melted_matrixE, aes(Var2, Var1, fill = value))+
		geom_tile(color = "white")+
		scale_fill_gradient2(low = "blue", high = "red", mid = 'white',
												 midpoint = 0, limit = c(-1,1), space = "Lab",
												 name="Pearson\\nCorrelation") +
		theme_minimal()+ 
		theme(axis.text.x = element_text(angle = 45, vjust = 1, 
																		 size = 12, hjust = 1))+
		coord_fixed()+
		geom_text(aes(Var2, Var1, label = value), color = "black", size = 4)

# traits (Fig 2b and 2c)
	

	w2 = ESLTP$col.w
	tab2 = Q[,c(2,6:8,10)]
	tab2 = scalewt(tab2,w2)
	corP1 <- (t(tab2)%*%diag(w2)%*%ESLTP$mQ[,1])[,1]
	corP2 <- (t(tab2)%*%diag(w2)%*%ESLTP$mQ[,2])[,1]
	corP = data.frame(corP1)
	corP[,2] = corP2
	
	#Spearman's correlation coefficient for binary variables
	hs_bin = dudi.hillsmith(Q_bin, row.w = afc$cw, scannf = F, nf = 2) #Hill & Smith
	ESLTP_bin = rlqESLTP(dudiE = acp, dudiS = acp_S, dudiL = afc, dudiT = hs_bin, 
											 dudiP = pco, scannf = FALSE, nf = 2)
	w2b = ESLTP_bin$col.w
	tab2b = Q_bin[,c(1:7,9:16,18:20,24:26)]
	tab2b = as.data.frame(apply(tab2b, 2, rank))
	tab2b = scalewt(tab2b,w2b)
	corS1 <-  t(tab2b)%*%diag(w2b)%*%scalewt(rank(ESLTP_bin$mQ[,1]), w2b)
	corS2 <-  t(tab2b)%*%diag(w2b)%*%scalewt(rank(ESLTP_bin$mQ[,2]), w2b)
	corS = data.frame(corS1)
	corS[,2] = corS2
	
	#graphical display (2b)
	colnames(corS) = c("Ax1","Ax2")
	matrixT = as.matrix(round(corS,2))
	melted_matrixT = melt(matrixT)
	if (graph > 0) {
		x11();ggplot(data = melted_matrixT, aes(Var1, Var2, fill = value))+
			geom_tile(color = "white")+
			scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
													 midpoint = 0, limit = c(-1,1), space = "Lab",
													 name="Spearman\\nCorrelation") +
			theme_minimal()+ 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, 
																			 size = 12, hjust = 1))+
			coord_fixed()+
			geom_text(aes(Var1, Var2, label = value), color = "black", size = 4)
	}
	
	# graphical display (2c)
	colnames(corP) = c("Ax1","Ax2")
	matrixT = as.matrix(round(corP,2))
	melted_matrixT = melt(matrixT)
	if (graph > 0) {
		x11();ggplot(data = melted_matrixT, aes(Var1, Var2, fill = value))+
			geom_tile(color = "white")+
			scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
													 midpoint = 0, limit = c(-1,1), space = "Lab",
													 name="Pearson\\nCorrelation") +
			theme_minimal()+ 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, 
																			 size = 24, hjust = 1))+
			coord_fixed()+
			geom_text(aes(Var1, Var2, label = value), color = "black", size = 7)
	}
	
#------------------------------------#
#### graphical display : Figure 3 ####
#------------------------------------#

# species
status$Status=factor(status$Status)
gspecies = s.class(ESLTP$lQ, fac=status$Status, paxes.draw = T, col = c("#D55E00","#0099CC","#000000",
																																				"purple"), ellipseSize = 1, pbackground.col = "white",
									 plabels.boxes = list(col = "black", alpha = 0,lwd = 2),
									 plabels.cex = 0, starSize = 1, pellipses.alpha = 0.3, porigin.lty = 'dotted',
									 plegend.size = 1.2, pgrid.draw = F, xlab = 'RLQ 1 (74%)', ylab = 'RLQ 2 (25%)')
glab = s.label(ESLTP$lQ, pbackground.col = "white", plabels = list(cex = 0.8, optim = T,
																																	 col = "black"), plabels.boxes = list(col = "black", border = "black",
																																	 																		 alpha = 0.4, lwd = 2), ppoints.alpha = 0, add = F, plabels.boxes.draw = F)

x11(); gespeces = superpose(gspecies,glab, plot = T)

# sites
FTyp=as.character(hab$forest_type)
FTyp[which(FTyp=="Exotic_forest")]="plantation forest"
FTyp[which(FTyp=="Not_forest")]="non forest"
FTyp[which(FTyp=="Indigeneous_forest")]="native forest"
FTyp=factor(FTyp)

x11();gsite = s.class(ESLTP$lR[c(1:404,406:629,631,632,634:917),],paxes.draw = T,
											FTyp[c(1:404,406:629,631,632,634:917)],
											ellipseSize = 1.5, col = c("#333333","#CCCCCC","#996633"),
											porigin.lty = 'dotted',
											plabels.optim = T, pbackground.col = "white", plabels.cex = 0,
											starSize = 0, pellipses.alpha = 0.8, plegend.drawKey = T, pgrid.draw = F,
											xlab = 'RLQ 1 (74%)', ylab = 'RLQ 2 (25%)')

### Figure 4 was entirely generated under qGIS with data from the ESLTP object ###

#--------------#
### Figure 5 ###
#--------------#

# first axis
phy = NZtree_P ; ax = 1
CB = cbind.data.frame(ESLTP$lQ_givenT[phy$tip.label, ax], ESLTP$lQ_givenP[phy$tip.label, 
																																					ax], ESLTP$lQ[phy$tip.label, ax])
colnames(CB) = c("trait-based", "phylogeny-based","global") 
if (graph > 0) {
	x11();dotchart.phylog(NZtree_Pphylog,CB)
}

# second axis
phy = NZtree_P ; ax = 2
CB = cbind.data.frame(ESLTP$lQ_givenT[phy$tip.label, ax], ESLTP$lQ_givenP[phy$tip.label, 
																																					ax], ESLTP$lQ[phy$tip.label, ax])
colnames(CB) = c("trait-based", "phylogeny-based","global") 
if (graph > 0) {
	x11();dotchart.phylog(NZtree_Pphylog,CB)
}

#--------------#
### Figure 6 ###
#--------------#

### signal phylogénétique sur les covariables


# phylogenetic signal
introsp=rownames(TBC2)
NZtree_alien=drop.tip(NZtree_P,NZtree_P$tip.label[-match(introsp, NZtree_P$tip.label)])
variables.for.test=TBC2[NZtree_alien$tip.label,]

date=variables.for.test$Date ; names(date)=rownames(variables.for.test)
test_signal=phylosig(NZtree_alien,date,method="lambda",test=T)

nbintro=variables.for.test$Nb_intro ; names(nbintro)=rownames(variables.for.test)
test_signal=phylosig(NZtree_alien,nbintro,method="lambda",test=T)

released=variables.for.test$Released_tot ; names(released)=rownames(variables.for.test)
test_signal=phylosig(NZtree_alien,released,method="lambda",test=T)

# non phylogenetic manova
TBC2$Ax1 = ESLTP$lQ[c(3,5,6,8,11,12,16,20,22,23,26,31,33,36,40:43,48),1]
TBC2$Ax2 = ESLTP$lQ[c(3,5,6,8,11,12,16,20,22,23,26,31,33,36,40:43,48),2]

# phylogenetically constrainted analysis
vari=as.matrix(TBC2[,c("Ax1","Ax2")])
Date=TBC2$Date ; names(Date)=rownames(TBC2)
Nb_intro=TBC2$Nb_intro ; names(Nb_intro)=rownames(Nb_intro)
Released_tot=TBC2$Released_tot ; names(Released_tot)=rownames(Released_tot)

# phylogenetic regression
phy.man=phylolm(vari~ Date+Nb_intro+Released_tot,phy=NZtree_alien,model="BM")

# partial residuals
p1=visreg(res.man,xvar="Date")
p2=visreg(res.man,xvar="Nb_intro")
p3=visreg(res.man,xvar="Released_tot")

# plots
par(mfrow=c(3,2))
plot(p1[[1]],bty="n",ylab="RLQ1",xlab="",points.par=list(cex=1,pch=21,bg="black",col="black"),line.par=list(col="black"),cex.lab=1.5)
mtext(text="(a)",side=3,at=1840,font=2,cex=1)
mtext(text="Date of first introduction",side=1,at=1930,line=3)
plot(p1[[2]],bty="n",ylab="RLQ2",xlab="",points.par=list(cex=1,pch=21,bg="black",col="black"),line.par=list(col="black"),cex.lab=1.5)
mtext(text="(b)",side=3,at=1840,font=2,cex=1)
plot(p2[[1]],bty="n",ylab="RLQ1",xlab="",points.par=list(cex=1,pch=21,bg="black",col="black"),line.par=list(col="black"),cex.lab=1.5)
mtext(text="(c)",side=3,at=5,font=2,cex=1)
mtext(text="Nb of introduction events",side=1,at=30,line=3)
plot(p2[[2]],bty="n",ylab="RLQ2",xlab="",points.par=list(cex=1,pch=21,bg="black",col="black"),line.par=list(col="black"),cex.lab=1.5)
mtext(text="(d)",side=3,at=5,font=2,cex=1)
plot(p3[[1]],bty="n",ylab="RLQ1",xlab="",points.par=list(cex=1,pch=21,bg="black",col="black"),line.par=list(col="black"),cex.lab=1.5)
mtext(text="(e)",side=3,at=100,font=2,cex=1)
mtext(text="Nb of individuals released",side=1,at=1000,line=3)
plot(p3[[2]],bty="n",ylab="RLQ2",xlab="",points.par=list(cex=1,pch=21,bg="black",col="black"),line.par=list(col="black"),cex.lab=1.5)
mtext(text="(f)",side=3,at=100,font=2,cex=1)
