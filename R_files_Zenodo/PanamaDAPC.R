##load libraries
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("devtools")
#install_github("thibautjombart/adegenet")
library("adegenet")
##load data##
Allpop<- read.genepop("PanamaDAPC.gen", ncode=2L, quiet=FALSE)
##look the data
Allpop

X1 <- tab(Allpop, freq=TRUE, NA.method="mean")
D1 <- dist(X1)

#call Populations
obj1 <- Allpop[sample(1:54, 54)] 
pop(obj1)


temp1 <- table(pop(Allpop))
temp2 <- summary(Allpop)

## Alleles per population 
barplot(temp2$pop.n.all, xlab="Population", ylab="Number of SNPS", col=funky(4))
##Heterozigozity Observed
barplot(temp2$Hobs, xlab="Locus", ylab="Ho", ylim =c(0,1))
##Heterozigosity Expected
barplot(temp2$Hexp, xlab = "Locus", ylab = "He", ylim =c(0,1))
###Expected heterozygozity is 50%
plot(temp2$Hexp, temp2$Hobs, pch=20, cex=3, xlim=c(.1,1), ylim=c(.1,1))
abline(0,1,lty=2)

n.pop <- seppop(Allpop)
mean.hobs <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs))) 

mean.hobs[is.nan(mean.hobs)] <- NA 

barplot(mean.hobs) 

mean.exp <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hexp))) 

mean.exp[is.nan(mean.exp)] <- NA 

barplot(mean.exp) 

plot(mean.exp, mean.hobs, pch=20, cex=3, xlim=c(.1,.9), ylim=c(.1,.9), text(mean.exp, mean.hobs, row.names(temp1), cex=0.6, pos=4, col="green"))
abline(0,1,lty=2)

plot(mean.exp, mean.hobs, pch=20, cex=3, xlim=c(.1,.4), ylim=c(.1,.4), text(mean.exp, mean.hobs, row.names(mean.hobs), cex=0.6, pos=4, col="green"))
abline(0,1,lty=2)


##Population sample size
barplot(table(pop(Allpop)), col=funky(4), las=3,
        xlab="Population", ylab="Sample size")



## Alleles per population 
barplot(temp2$pop.n.all, xlab="Population", ylab="Number of SNPS", col=funky(4))

##Temp2 info
str (temp2)

##load hierfstat library
library(hierfstat)

#Estadistico F
fstat(Allpop)
fstat(Allpop, fstonly=TRUE)

#Paired Fst
AllpopFst <- pairwise.fst(Allpop,res.type= c("matrix"))
AllpopFst[1:4,1:4]


##NJ
crocs.tree <- nj(AllpopFst)
plot(crocs.tree, type="unr", tip.col=funky(nPop(Allpop)), font=0.2)
annot <- round(crocs.tree$edge.length,2)
edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()
table.paint(AllpopFst, col.labels=1:4)

temp2 <- AllpopFst
diag(temp2) <- NA
boxplot(temp2, col=funky(nPop(Allpop)), las=3,
        xlab= NULL, ylab="Fst")
boxplot(temp2, col=funky(nPop(Allpop)), notch=TRUE , las=3,
        xlab= NULL, ylab="Fst")
##PCA
x.crocs <- tab(Allpop, freq=TRUE, NA.method="mean")
pca.crocs <- dudi.pca(x.crocs, center=TRUE, scale=FALSE)

 
pca.crocs

s.label(pca.crocs$li, clabel=0.7, pch=0.5, sub = "PCA")

s.class(pca.crocs$li, fac=pop(Allpop), clabel=0.5, col=funky(16), axesel = FALSE, cstar = 0, cgrid = 0, cpoint=2, cellipse = 2)

library(wordcloud)
## Loading required package: RColorBrewer
textplot(pca.crocs$li[,1], pca.crocs$li[,2], words=rownames(x.crocs), cex=0.8, new=TRUE, xlim=c(-35,35), ylim=c(-35,35))


##PCA 1 vs 2
s.class(pca.crocs$li, fac=pop(Allpop),
        col=transp(funky(4),.6),
        axesel=FALSE, cstar=0, cpoint=3, clabel=0.7)
add.scatter.eig(pca.crocs$eig[1:7],3,1,2, ratio=.2, posi="bottomleft")

##PCA 2 vs 3
s.class(pca.crocs$li, fac=pop(Allpop),
        xax=2, yax=3, col=transp(funky(15),.6),
        axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.crocs$eig[1:7],3,2,3, ratio=.20, posi="bottomleft")

##PCA 1 vs 3
s.class(pca.crocs$li, fac=pop(Allpop),
        xax=1, yax=3, col=transp(funky(15),.6),
        axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.crocs$eig[1:7],3,2,3, ratio=.20, posi="bottomleft")

pca.crocs$eig[1]

pc1 <- pca.crocs$li[,1]
var(pc1)
var(pc1)*54/54
mean(pc1^2)
n <- length(pc1)
0.5*mean(dist(pc1)^2)*((n-1)/n)

eig.perc <- 100*pca.crocs$eig/sum(pca.crocs$eig)
head(eig.perc)

loadingplot(pca.crocs$c1^2)

contrib <- loadingplot(pca.crocs$c1^2, axis=2, factor=pop, byfac=TRUE,
                       thres=quantile(pca.crocs$c1[,2]^2, .98) )


contrib <- loadingplot(pca.crocs$c1^2, axis=2,
                       thres=quantile(pca.crocs$c1[,2]^2, .98) )
contrib$var.names

##Kmeans
set.seed(1)
crocgrp<-find.clusters(Allpop, max.n.clust =20, stat= "BIC", scale=FALSE)

##crocgrp<-find.clusters(Batch69NM1BCHpop, max.n.clust = 20, stat= "AIC")##

plot(crocgrp$Kstat, type="o", xlab="Number of clusters (K)", ylab="BIC",
     col="blue", main="Detection based on BIC")

names(crocgrp)
head(crocgrp$Kstat,30)
head(crocgrp$grp, 54)

table.value(table(pop(Allpop), crocgrp$grp), col.lab=paste("K", 1:3))


#DAPC 
dapc1 <- dapc(Allpop, pop=crocgrp$grp, scale=FALSE)

dapc1 
scatter(dapc1, col=funky(3), scree.pca=FALSE)

scatter(dapc1, col=funky(3),legend=TRUE, mstree=TRUE,
        cstar=0, axesell=FALSE, clab=0, cex=15, bg=grey(.2),
        scree.pca=TRUE,segcol="white")


scatter(dapc1, posi.da="bottomleft", grp=dapc1$grp, pch=30, cex.lab=1, col=funky(3), cstar = 5, cellipse = 12)

scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=funky(3), solid=.5,
        cex=2,clab=0.75, leg=TRUE, txt.leg=paste("K",1:3), cellipse = 5)

scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=funky(14), solid=.5,
        cex=2,clab=0, leg=TRUE, txt.leg=paste("K",1:14))


scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=funky(14), solid=.4, cex=2.5, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, txt.leg=paste("K",1:11))

par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=3,
       cex=3, lwd=8, col="black")

points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=funky(14))


scatter(dapc1, col=funky(3), scree.pca=TRUE, posi.da="bottomleft",
        posi.pca="topleft", cex=5)

contrib <- loadingplot(dapc1$var.contr, axis=2,
                       thres=quantile(pca.crocs$c1[,2]^2, .99) )
contrib$var.names

temp <- seploc(Allpop)
snp146667 <- tab(temp("146667_105.03"))

set.seed(4)
contrib <- loadingplot(dapc1$var.contr, axis=2,
                       thres=.07, lab.jitter=1)


assignplot(dapc1, subset=1:54)
compoplot(dapc1, posi="bottomleft",
          txt.leg=paste("K", 1:3),
          ncol=0.5, xlab="individuals", col=funky(3))

dapc3 <- dapc(Allpop, n.da=2, n.pca=45)
compoplot(dapc3, subset=1:54, posi="bottomright",
          txt.leg=paste("Cluster", 1:3),
          ncol=1, col=funky(13))
