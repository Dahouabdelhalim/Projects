# This R code includes code for all the R analyses in the Guo et al. paper,
#   as well as for other analyses not included in that paper.
#
# A few of the code chunks take a LONG time to execute!

## ----libraries----------------------------------------------------------------
 library(car)
 library(reshape)
 library(ggplot2)
 library(xtable)
 library(MASS)
 library(fields)
 library(vegan)
 library(cluster)
 library(magrittr)
 library(mmod)
 library(ape)
 library(ade4)
 library(adegenet)
 library(stringr)
 library(ggrepel)
 library(viridis)
 library(gplots)
 library(ggfortify)
 library(reshape2)
 library(gridExtra)
 library(Rfast)
 library(Matrix)
 library(cowplot)
 library(poppr)


## ----subset-------------------------------------------------------------------
 Subset <- function (x, ...) {
     droplevels(subset(x, ...))
}


## ----Read---------------------------------------------------------------------
gendat <- read.csv("AllPatch-v3.csv", sep=",", strip.white=TRUE)
areas <- read.csv("PatchSizes.csv")
gendat$Patch <- factor(gendat$Patch, levels = c("1", "2", "3",
    "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"))


## ----dim----------------------------------------------------------------------
d.gendat <- dim(gendat)


## ----basic--------------------------------------------------------------------
start.gendat <- 7
Nloci <- dim(gendat)[2] - (start.gendat -1)
Nindivs <- dim(gendat)[1]
Npatch <- length(levels(gendat$Patch))
Patches <- as.character(levels(gendat$Patch))


## ----Inf-loci-fn--------------------------------------------------------------
fPolymorphicLoci <- function(DataFrame){
  Nloci - sum(duplicated(t(DataFrame[,start.gendat:d.gendat[2]])))
}


## ----patches------------------------------------------------------------------
Npolymorphic.loci.patch <- vector()
Richness <- vector()
Nindivs.marker <- vector()
count <- 0
for (i in 1:length(Patches)){
  count <- count + 1
  data.name <- paste("P",Patches[i],".dat", sep="")
  tmp <- subset(gendat, gendat$Patch == levels(gendat$Patch)[count])
  assign(data.name, tmp)
  Npolymorphic.loci.patch[count] <- fPolymorphicLoci(tmp)
  Richness[count] <- sum(as.logical(rowSums(t(tmp[,start.gendat:d.gendat[2]]))))
  newrow <- as.numeric(rowSums(t(tmp[,start.gendat:d.gendat[2]])))
  Nindivs.marker <- rbind(Nindivs.marker,newrow)
}
Nindiv.patch <- table(gendat$Patch)
List.of.Data.Frames <- paste(Patches,".dat",sep="")
rownames(Nindivs.marker) <- Patches


## ----results="asis", echo=FALSE-----------------------------------------------

    t.rich <- rbind(Patches,Richness)
    rownames(t.rich) <- c("Patch","Marker richness")
    tab.rich <- xtable(t.rich, caption="Marker richness, by patch.", label="tab:rich")
    tab.rich
    #print(tab.summary, sanitize.text.function = function(x) {x})


## ----figheat, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----
#my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
heat<- heatmap.2(Nindivs.marker, col=viridis, density.info="none", trace="none",
   dendrogram="none", labRow=Patches, xlab="Locus", ylab="Patch", cexRow=1.3)


## ----figfreqmarker------------------------------------------------------------
 freq.marker <-   sweep(Nindivs.marker, 1, Nindiv.patch, "/")
 heatfreq<- heatmap.2(freq.marker,col=viridis, density.info="none", trace="none",
   dendrogram="none", labRow=Patches, xlab="Locus", ylab="Patch", cexRow=1.3)


## ----figheatfreq, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----
heatfreq<- heatmap.2(freq.marker, col=viridis, density.info="none", trace="none",
   dendrogram="none", labRow=Patches, xlab="Locus", ylab="Patch", cexRow=1.3)


## ----figdens, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----
par(mfrow = c(3,5))
for (i in 1:13){
    hist(subset(freq.marker[i,], freq.marker[i,] > 0), col="black",
      xlab="",ylab="",
      main = paste("Patch",Patches[i]))
}
mtext("Marker frequency", outer=TRUE,side=1, line = -2)
mtext("# of markers", outer=TRUE,side=2, line = -2)



## ----U.genotypes--------------------------------------------------------------
Nduplicated.genotypes <- sum(duplicated(gendat[,start.gendat:d.gendat[2]]))
Nunique.genotypes <- dim(gendat)[1] - Nduplicated.genotypes
Nunique.genotypes.patch <- Nindiv.patch
Nunique.genotypes
Nunique.genotypes.patch


## ----Inf.loci-----------------------------------------------------------------
N.monomorphic.loci <- sum(duplicated(t(gendat[,start.gendat:d.gendat[2]])))
Npolymorphic.loci <- Nloci - sum(duplicated(t(gendat[,start.gendat:d.gendat[2]])))
N.monomorphic.loci
Npolymorphic.loci


## ----Inf.loci.patch-----------------------------------------------------------
Npolymorphic.loci.patch


## ----Rich.patch---------------------------------------------------------------
Richness


## ----List.inf.loci------------------------------------------------------------
List.monomorphic.loci <- which (duplicated(t(gendat[,start.gendat:d.gendat[2]]))==TRUE)
List.polymorphic.loci <- which (duplicated(t(gendat[,start.gendat:d.gendat[2]]))==FALSE)


## ----Summary-table------------------------------------------------------------
Sample <- c(Patches,"Total")
N.indivs <- c(Nindiv.patch, Nindivs)
N.unique <- c(Nunique.genotypes.patch, Nunique.genotypes)
N.Polymorphic <- as.integer(c(Npolymorphic.loci.patch, Npolymorphic.loci))
N.Monomorphic <- as.integer(c(Nloci - Npolymorphic.loci.patch, Nloci - Npolymorphic.loci))
Sum <- data.frame(Sample, N.indivs, N.unique, N.Polymorphic, N.Monomorphic)


## ----results="asis", echo=FALSE-----------------------------------------------
    tab.summary <- xtable(Sum, caption="Summary of sample size, unique genotypes, polymorphic loci.", label="tab:summary")
    tab.summary
    #print(tab.summary, sanitize.text.function = function(x) {x})


## ----tmpfig-------------------------------------------------------------------
regexp <- "[[:digit:]]+"
 patchnum <- as.numeric(Patches)
 tmpfig <- data.frame(Patch = patchnum, N.indivs = Sum$N.indivs[1:13], N.Polymorphic =
   Sum$N.Polymorphic[1:13], Area = areas$area)


## ----figsampleinf, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----
 p <- ggplot(tmpfig, aes(x = N.indivs, y = N.Polymorphic))  + labs(y = "Number of polymorphic loci", x = "Number of individuals sampled") + geom_point()
 p + geom_text_repel(aes(label=Patch), box.padding = unit(0.45, "lines"))


## ----genalex.format1----------------------------------------------------------
tmp <- gendat
tmp$Location <- NULL
tmp$Zone <- NULL

tmp <- cbind(Ramet=tmp$Sample.No., tmp)
tmp$Sample.No. <- NULL


## ----genalex.format2----------------------------------------------------------
tmp$Ramet <- paste(tmp$Patch, tmp$Ramet, sep = "-")
tmp$Ramet <- gsub("Patch","",tmp$Ramet)


## ----genalex.format3----------------------------------------------------------
X <- tmp$x
Y <- tmp$y
tmp$x <- NULL
tmp$y <- NULL
fileheader <-  vector(mode="character")
fileheader <- (names(tmp))
fileheader <- c(fileheader, "")
fileheader <- c(fileheader, c("X","Y"))
fileheader <- gsub("Locus.","",fileheader)
line1 <- vector(mode="character",length=191)
line2 <- vector(mode="character",length=191)
l1 <- c(Nloci,Nindivs,Npatch, Nindiv.patch)
line1[1:length(l1)] <- l1
#alpha.patches <- c("Patch1", "Patch2", "Patch3", "Patch5", "Patch6",
#  "Patch7", "Patch8", "Patch9", "Patch10", "Patch11", "Patch12",
#  "Patch13", "Patch14")
l2 <- c("Leymus data", "", "13 patches", "", Patches)
line2[1:length(l2)] <- l2


## ----genalex.out, echo=FALSE--------------------------------------------------
sink("out.csv", type="output")

write.table(t(line1), file="out.csv", sep=",", col.names=FALSE,row.names=FALSE, quote=FALSE, append=FALSE)
write.table(t(line2), file="out.csv",  sep=",", col.names=FALSE,row.names=FALSE, quote=FALSE, append=TRUE)
write.table(t(fileheader), file="out.csv",  sep=",", col.names=FALSE,row.names=FALSE, quote=FALSE, append=TRUE)
write.table(
  cbind(tmp,vector(mode="character",length=dim(tmp)[[1]]),X,Y),
    col.names=FALSE,row.names=FALSE, quote=FALSE, sep=",", file="out.csv", append=TRUE
  )
sink()


## ----read.genalex-------------------------------------------------------------
 (genalex2 <- read.genalex("out.csv", ploidy=4, geo=TRUE))
 (poppdata<- poppr::poppr(genalex2))


## ----results="asis", echo=FALSE-----------------------------------------------
    tab.poppdata <- xtable(poppdata[,1:12], caption = "Diversity
        estimates from naive data set, assuming each individual
        is unique.", label="tab:naive-diverse")

    tab.poppdata
    #print(tab.summary, sanitize.text.function = function(x) {x})


## ----occurrence-sums----------------------------------------------------------
tmpdat <- gendat[,7:192]
occurrence_sums <- colSums(tmpdat)
o.occur <- sort(occurrence_sums)
plot(o.occur, xlab="Rank order of occurrences", ylab="Occurrences", main="")
plot(o.occur, xlab="Rank order of occurrences", ylab="Occurrences", main="", log="y")


## ----nmarker-diffs,cache=TRUE-------------------------------------------------
mydata <- gendat[,7:192]

nr <- 357
difmat <- matrix(nrow=nr, ncol=nr+1,0)
for(i in 1:nr){
    for (j in (i+1):nr){
        difmat[i,j] <- nnzero (mydata[i,] - mydata[j,])
        }
}
diftab <- table(difmat[difmat > 0])
plot(diftab/sum(diftab), xlab="Number of marker differences", ylab="Density", main="")


## ----diff10-------------------------------------------------------------------
 1 - sum(diftab[1:10])/sum(diftab)


## ----pc-----------------------------------------------------------------------
 pc <- as.genclone(genalex2, threads = 1L)


## ----ngen---------------------------------------------------------------------
pvdist.tot <- prevosti.dist(pc)
ngen <- vector()
cuts <-       seq(0.04,0.08, 0.001)
for (i in 1:length(cuts)){
    mlg.filter(pc, distance = pvdist.tot) <- cuts[i]
    ngen[i] <- length(unique(pc$mlg))[[1]]
}
ngen.df <- data.frame(cuts,ngen)


## ----figngen, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----
  ggplot(ngen.df, aes(x=cuts, y=ngen)) +
    geom_point(colour = "red", size = 3) +
    labs(x="Threshold", y = "Number of MLLs")



## ----pc2----------------------------------------------------------------------
mlg.filter(pc, distance = pvdist.tot) <- 0.064
pc


## ----diversity-pc-------------------------------------------------------------
 (poppc <- poppr::poppr(pc))


## ----results="asis", echo=FALSE-----------------------------------------------
    tab.poppc <- xtable(poppc[,1:12], caption = "Diversity
        estimates from collapsed data set, using 13 clones.",
        label="tab:collapsed-diverse")
    tab.poppc
    #print(tab.summary, sanitize.text.function = function(x) {x})


## ----wss----------------------------------------------------------------------

wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))

for (i in 2:15) wss[i] <- sum(kmeans(mydata,
  centers=i, nstart=50)$withinss)


## ----n-pca--------------------------------------------------------------------
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")


## ----fit-kmeans---------------------------------------------------------------
fit <- kmeans(mydata, 13)

## ----agg, echo=FALSE----------------------------------------------------------
agg <- aggregate(mydata,by=list(fit$cluster),FUN=mean)


## ----clusplot2----------------------------------------------------------------
zz <- as.integer(gendat$Patch)
pca1 <- autoplot(pam(mydata,13), colour="Patch", data=gendat)

pca2 <- autoplot(pam(mydata,13), data=gendat)


## ----clusplot2a---------------------------------------------------------------
grid.arrange(ggplotGrob(pca1),ggplotGrob(pca2), heights=c(2,1.6))


## ----subsets------------------------------------------------------------------
 p1 <- popsub(genalex2,"1")
 p2 <- popsub(genalex2,"2")
 p3 <- popsub(genalex2,"3")
 p5 <- popsub(genalex2,"5")
 p6 <- popsub(genalex2,"6")
 p7 <- popsub(genalex2,"7")
 p8 <- popsub(genalex2,"8")
 p9 <- popsub(genalex2,"9")
 p10 <- popsub(genalex2,"10")
 p11 <- popsub(genalex2,"11")
 p12 <- popsub(genalex2,"12")
 p13 <- popsub(genalex2,"13")
 p14 <- popsub(genalex2,"14")

## ----test.diseq---------------------------------------------------------------
 diseq <- matrix(data=NA, nrow=Npatch, ncol=4)
 diseq[1,] <- ia(p1, sample=999, plot=FALSE, quiet=TRUE)
 diseq[2,] <- ia(p2, sample=999, plot=FALSE, quiet=TRUE)
 diseq[3,] <- ia(p3, sample=999, plot=FALSE, quiet=TRUE)
 diseq[4,] <- ia(p5, sample=999, plot=FALSE, quiet=TRUE)
 diseq[5,] <- ia(p6, sample=999, plot=FALSE, quiet=TRUE)
 diseq[6,] <- ia(p7, sample=999, plot=FALSE, quiet=TRUE)
 diseq[7,] <- ia(p8, sample=999, plot=FALSE, quiet=TRUE)
 diseq[8,] <- ia(p9, sample=999, plot=FALSE, quiet=TRUE)
 diseq[9,] <- ia(p10, sample=999, plot=FALSE, quiet=TRUE)
 diseq[10,] <- ia(p11, sample=999, plot=FALSE, quiet=TRUE)
 diseq[11,] <- ia(p12, sample=999, plot=FALSE, quiet=TRUE)
 diseq[12,] <- ia(p13, sample=999, plot=FALSE, quiet=TRUE)
 diseq[13,] <- ia(p14, sample=999, plot=FALSE, quiet=TRUE)

disq <- as.data.frame(diseq)
disq <- cbind(l2[5:17], disq)
colnames(disq) <- c("Patch","Ia","P.Ia","rbarD","P.rbarD")


## ----results="asis", echo=FALSE-----------------------------------------------
    tab.summary <- xtable(disq, caption="Index of association ($I_a$) for the naive data, normalized index of association ($\\\\bar{r}_{D}$), and $P$ values for each one, testing (by resampling) the hypothesis that the value equals zero.", label="tab:disq")
    colnames(tab.summary) <- c("Patch","$I_{a}$", "$P$", "$\\\\bar{r}_{D}$","$P$")
    print(tab.summary, include.rownames=FALSE, sanitize.text.function = function(x) {x})


## ----subsets.c----------------------------------------------------------------
 p1.c <- popsub(pc,"1")
 p2.c <- popsub(pc,"2")
 p3.c <- popsub(pc,"3")
 p5.c <- popsub(pc,"5")
 p6.c <- popsub(pc,"6")
 p7.c <- popsub(pc,"7")
 p8.c <- popsub(pc,"8")
 p9.c <- popsub(pc,"9")
 p10.c <- popsub(pc,"10")
 p11.c <- popsub(pc,"11")
 p12.c <- popsub(pc,"12")
 p13.c <- popsub(pc,"13")
 p14.c <- popsub(pc,"14")


## ----test.diseq.c-------------------------------------------------------------
 diseq.c <- matrix(data=NA, nrow=Npatch, ncol=4)
 diseq.c[1,] <- ia(p1.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[2,] <- ia(p2.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[3,] <- ia(p3.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[4,] <- ia(p5.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[5,] <- ia(p6.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[6,] <- ia(p7.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[7,] <- ia(p8.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[8,] <- ia(p9.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[9,] <- ia(p10.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[10,] <- ia(p11.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[11,] <- ia(p12.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[12,] <- ia(p13.c, sample=999, plot=FALSE, quiet=TRUE)
 diseq.c[13,] <- ia(p14.c, sample=999, plot=FALSE, quiet=TRUE)

 disq.c <- as.data.frame(diseq.c)
 disq.c <- cbind(l2[5:17],disq.c)
 colnames(disq.c) <- c("Patch","Ia","P.Ia","rbarD","P.rbarD")


## ----results="asis", echo=FALSE-----------------------------------------------
    tab.summary <- xtable(disq.c, caption="Index of association ($I_a$) for the MLL data, normalized index of association ($\\\\bar{r}_{D}$), and $P$ values for each one, testing (by resampling) the hypothesis that the value equals zero.", label="tab:disq.c")
    colnames(tab.summary) <- c("Patch","$I_{a}$", "$P$", "$\\\\bar{r}_{D}$","$P$")
    print(tab.summary, include.rownames=FALSE, sanitize.text.function = function(x) {x})


## ----make.tree----------------------------------------------------------------
 ddist <- diss.dist(genalex2)
 theTree <- ddist %>%
   nj() %>%
     ladderize()


## ----make.clusters------------------------------------------------------------
clusters <- find.clusters(genalex2,n.clust=13, n.pca=150)


## ----figplottree, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE----
set.seed(2342)
cols<-rainbow(13)
plot.phylo(theTree, cex=0.4, tip.color=cols[clusters$grp],
  label.offset= 0.006125)
#nodelabels(theTree$node.label, adj = c(1.3, -0.5),
#  frame = "n", cex = 0.4,font=1,xpd=TRUE)
axisPhylo(3)


## ----clusptch-----------------------------------------------------------------
clus <- pam(pc,13)$clustering
ptch <- sub('\\\\-.*', '', names(clus))
#table(ptch,clus)


## ----results="asis", echo=FALSE-----------------------------------------------
    tab.ptchclus <- xtable(table(ptch,clus), caption="Assignment of ramets to clusters, by patch (clusters
        across top, patches down left side).", label="tab:ptchclus")
    print(tab.ptchclus, include.rownames=TRUE, sanitize.text.function = function(x) {x})


## ----genetic.distances, echo=FALSE--------------------------------------------

P14<- popsub(genalex2, sublist="14")
pvdist.14 <- provesti.dist(P14)

P13<- popsub(genalex2, sublist="13")
pvdist.13 <- provesti.dist(P13)

P5<- popsub(genalex2, sublist="5")
pvdist.5 <- provesti.dist(P5)

P7<- popsub(genalex2, sublist="7")
pvdist.7 <- provesti.dist(P7)

P11<- popsub(genalex2, sublist="11")
pvdist.11 <- provesti.dist(P11)

P1<- popsub(genalex2, sublist="1")
pvdist.1 <- provesti.dist(P1)

P12<- popsub(genalex2, sublist="12")
pvdist.12 <- provesti.dist(P12)

P10<- popsub(genalex2, sublist="10")
pvdist.10 <- provesti.dist(P10)

P9<- popsub(genalex2, sublist="9")
pvdist.9 <- provesti.dist(P9)

P2<- popsub(genalex2, sublist="2")
pvdist.2 <- provesti.dist(P2)

P3<- popsub(genalex2, sublist="3")
pvdist.3 <- provesti.dist(P3)

P6<- popsub(genalex2, sublist="6")
pvdist.6 <- provesti.dist(P6)

P8<- popsub(genalex2, sublist="8")
pvdist.8 <- provesti.dist(P8)


## ----physical.distances,echo=FALSE--------------------------------------------
P1.loc <- P1.dat[, 4:5]
P1.dist <- rdist(P1.loc, P1.loc)
P1.test <- mantel.rtest (as.dist(pvdist.1),as.dist(P1.dist), nrepet=1000)

P2.loc <- P2.dat[, 4:5]
P2.dist <- rdist(P2.loc, P2.loc)
P2.test <- mantel.rtest (as.dist(pvdist.2),as.dist(P2.dist), nrepet=1000)

P3.loc <- P3.dat[, 4:5]
P3.dist <- rdist(P3.loc, P3.loc)
P3.test <- mantel.rtest (as.dist(pvdist.3),as.dist(P3.dist), nrepet=1000)

P5.loc <- P5.dat[, 4:5]
P5.dist <- rdist(P5.loc, P5.loc)
P5.test <- mantel.rtest (as.dist(pvdist.5),as.dist(P5.dist), nrepet=1000)

P6.loc <- P6.dat[, 4:5]
P6.dist <- rdist(P6.loc, P6.loc)
P6.test <- mantel.rtest (as.dist(pvdist.6),as.dist(P6.dist), nrepet=1000)

P7.loc <- P7.dat[, 4:5]
P7.dist <- rdist(P7.loc, P7.loc)
P7.test <- mantel.rtest (as.dist(pvdist.7),as.dist(P7.dist), nrepet=1000)

P8.loc <- P8.dat[, 4:5]
P8.dist <- rdist(P8.loc, P8.loc)
P8.test <- mantel.rtest (as.dist(pvdist.8),as.dist(P8.dist), nrepet=1000)

P9.loc <- P9.dat[, 4:5]
P9.dist <- rdist(P9.loc, P9.loc)
P9.test <- mantel.rtest (as.dist(pvdist.9),as.dist(P9.dist), nrepet=1000)

P10.loc <- P10.dat[, 4:5]
P10.dist <- rdist(P10.loc, P10.loc)
P10.test <- mantel.rtest (as.dist(pvdist.10),as.dist(P10.dist), nrepet=1000)

P11.loc <- P11.dat[, 4:5]
P11.dist <- rdist(P11.loc, P11.loc)
P11.test <- mantel.rtest (as.dist(pvdist.11),as.dist(P11.dist), nrepet=1000)

P12.loc <- P12.dat[, 4:5]
P12.dist <- rdist(P12.loc, P12.loc)
P12.test <- mantel.rtest (as.dist(pvdist.12),as.dist(P12.dist), nrepet=1000)

P13.loc <- P13.dat[, 4:5]
P13.dist <- rdist(P13.loc, P13.loc)
P13.test <- mantel.rtest (as.dist(pvdist.13),as.dist(P13.dist), nrepet=1000)

P14.loc <- P14.dat[, 4:5]
P14.dist <- rdist(P14.loc, P14.loc)
P14.test <- mantel.rtest (as.dist(pvdist.14),as.dist(P14.dist), nrepet=1000)


## ----mantel.tests-------------------------------------------------------------
obs <- c(P1.test$obs, P2.test$obs, P3.test$obs, P5.test$obs, P6.test$obs, P7.test$obs, P8.test$obs, P9.test$obs, P10.test$obs, P11.test$obs, P12.test$obs, P13.test$obs, P14.test$obs)

pvalues <- c(P1.test$pvalue, P2.test$pvalue, P3.test$pvalue, P5.test$pvalue, P6.test$pvalue, P7.test$pvalue, P8.test$pvalue, P9.test$pvalue, P10.test$pvalue, P11.test$pvalue, P12.test$pvalue, P13.test$pvalue, P14.test$pvalue)

mantel.results <- cbind(Patch=as.numeric(Patches), obs,pvalues)
mantel.results


## ----mantel.tests.MLL---------------------------------------------------------
for (i in Patches){
  xt <- genclone2genind(subset(pc, pc$strata == i))
  locations <- xt$other$xy
  dist <- rdist(locations,locations)
  test <- mantel.rtest(as.dist(dist),provesti.dist(xt), nrepet=1000)
  testname <- paste("P",i,".test.MLL",sep="")
  assign(testname,test)
 }


## ----mantel.tests.comb.MLL----------------------------------------------------
obs.MLL <- c(P1.test.MLL$obs, P2.test$obs, P3.test$obs, P5.test.MLL$obs, P6.test.MLL$obs, P7.test.MLL$obs, P8.test.MLL$obs, P9.test.MLL$obs, P10.test.MLL$obs, P11.test.MLL$obs, P12.test.MLL$obs, P13.test.MLL$obs, P14.test.MLL$obs)

pvalues.MLL <- c(P1.test.MLL$pvalue, P2.test.MLL$pvalue, P3.test.MLL$pvalue, P5.test.MLL$pvalue, P6.test.MLL$pvalue, P7.test.MLL$pvalue, P8.test.MLL$pvalue, P9.test.MLL$pvalue, P10.test.MLL$pvalue, P11.test.MLL$pvalue, P12.test.MLL$pvalue, P13.test.MLL$pvalue, P14.test.MLL$pvalue)

mantel.results.MLL <- cbind(Patch=as.numeric(Patches), obs.MLL, pvalues.MLL)
mantel.results.MLL


## ----p.values.mantel----------------------------------------------------------
cor(mantel.results[,3],mantel.results.MLL[,3])


## ----mantel.corr--------------------------------------------------------------
P1.correlog <- mantel.correlog(pvdist.1, D.geo=P1.dist,nperm=99)
P2.correlog <- mantel.correlog(pvdist.2, D.geo=P2.dist,nperm=99)
P3.correlog <- mantel.correlog(pvdist.3, D.geo=P3.dist,nperm=99)
P5.correlog <- mantel.correlog(pvdist.5, D.geo=P5.dist,nperm=99)
P6.correlog <- mantel.correlog(pvdist.6, D.geo=P6.dist,nperm=99)
P7.correlog <- mantel.correlog(pvdist.7, D.geo=P7.dist,nperm=99)
P8.correlog <- mantel.correlog(pvdist.8, D.geo=P8.dist,nperm=99)
P9.correlog <- mantel.correlog(pvdist.9, D.geo=P9.dist,nperm=99)
P10.correlog <- mantel.correlog(pvdist.10, D.geo=P10.dist,nperm=99)
P11.correlog <- mantel.correlog(pvdist.11, D.geo=P11.dist,nperm=99)
P12.correlog <- mantel.correlog(pvdist.12, D.geo=P12.dist,nperm=99)
P13.correlog <- mantel.correlog(pvdist.13, D.geo=P13.dist,nperm=99)
P14.correlog <- mantel.correlog(pvdist.14, D.geo=P14.dist,nperm=99)


## ----figmantelcorr, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE----
par(mfrow = c(3,5))
plot(P1.correlog)
mtext(text="1", side=3, line = 1)
plot(P2.correlog)
mtext(text="2", side=3, line = 1)
plot(P3.correlog)
mtext(text="3", side=3, line = 1)
plot(P5.correlog)
mtext(text="5", side=3, line = 1)
plot(P6.correlog)
mtext(text="6", side=3, line = 1)
plot(P7.correlog)
mtext(text="7", side=3, line = 1)
plot(P8.correlog)
mtext(text="8", side=3, line = 1)
plot(P9.correlog)
mtext(text="9", side=3, line = 1)
plot(P10.correlog)
mtext(text="10", side=3, line = 1)
plot(P11.correlog)
mtext(text="11", side=3, line = 1)
plot(P12.correlog)
mtext(text="12", side=3, line = 1)
plot(P13.correlog)
mtext(text="13", side=3, line = 1)
plot(P14.correlog)
mtext(text="14", side=3, line = 1)


## ----read.areas---------------------------------------------------------------
mantel.results.d <- as.data.frame(mantel.results)
mantel.results.d$area <- areas$area


## ----figmantel, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE------
p <- ggplot(mantel.results.d, aes(x = area, y = obs)) + geom_point() + labs(y = expression(paste("Observed ", rho)), x = expression(paste("Area in ",m^2)))
p + geom_text_repel(aes(label=Patches), box.padding = unit(0.45, "lines"))


## ----nei----------------------------------------------------------------------
y <- mantel.results.d
y$Patch  <- as.character(y$Patch)
y$Hexp <- poppdata$Hexp[1:13]


## ----regnei-------------------------------------------------------------------
lmnei <- lm(y$Hexp ~ log(y$area,10))
lmnei


## ----fignei, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE---------
p <- ggplot(y, aes(x = area, y = Hexp)) + geom_point() +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_log10() +
    labs(y = expression(H[exp]), x = expression(paste("Area in ",m^2)))
p + geom_text_repel(aes(label=Patches), box.padding = unit(0.45, "lines"))


## ----regdiv-------------------------------------------------------------------
y$ShannonW <- poppdata$H[1:13]
y$Simpson <- poppdata$lambda[1:13]
lmShannonW <- lm(y$ShannonW ~ log(y$area,10))
lmSimpson <- lm(y$Simpson ~ log(y$area,10))
summary(lmnei)
summary(lmSimpson)
summary(lmShannonW)


## ----figShannonW, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE----
p <- ggplot(y, aes(x = area, y = ShannonW)) + geom_point() +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_log10() +
    labs(y = "H", x = expression(paste("Area in ",m^2)))
p + geom_text_repel(aes(label=Patches), box.padding = unit(0.45, "lines"))


## ----figSimpson, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE-----
p <- ggplot(y, aes(x = area, y = Simpson)) + geom_point() +
    geom_smooth(method='lm', se=FALSE) +
    scale_x_log10() +
    labs(y = expression(paste(lambda)), x = expression(paste("Area in ",m^2)))
p + geom_text_repel(aes(label=Patches), box.padding = unit(0.45, "lines"))


## ----adjust.y-----------------------------------------------------------------
 y$p <- as.numeric(y$Patch)



## ----cor.ia-------------------------------------------------------------------
 cor.test(y$area,disq$Ia)


## ----cor.rbard----------------------------------------------------------------
 cor.test(y$area,disq$rbarD)


## ----N.inf--------------------------------------------------------------------
Sum$Sample <- as.character(Sum$Sample)
z <- Sum
z <- z[order(z$Sample),]
y$N.polymorphic <- z$N.Polymorphic[1:13]


## ----figninf, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE--------
p <- ggplot(y, aes(x = area, y = N.polymorphic)) +geom_point() + scale_x_log10() + labs(y = "N polymorphic loci", x = expression(paste("Area in ",m^2)))
p + geom_text_repel(aes(label=Patches), box.padding = unit(0.45, "lines"))


## ----spatial-genetic, tidy=TRUE, echo=FALSE-----------------------------------
t1 <- data.frame(SpatialDist=as.vector(as.dist(P1.dist)), GeneticDist=as.vector(as.dist(pvdist.1)))
t2 <- data.frame(SpatialDist=as.vector(as.dist(P2.dist)), GeneticDist=as.vector(as.dist(pvdist.2)))
t3 <- data.frame(SpatialDist=as.vector(as.dist(P3.dist)), GeneticDist=as.vector(as.dist(pvdist.3)))
t5 <- data.frame(SpatialDist=as.vector(as.dist(P5.dist)), GeneticDist=as.vector(as.dist(pvdist.5)))
t6 <- data.frame(SpatialDist=as.vector(as.dist(P6.dist)), GeneticDist=as.vector(as.dist(pvdist.6)))
t7 <- data.frame(SpatialDist=as.vector(as.dist(P7.dist)), GeneticDist=as.vector(as.dist(pvdist.7)))
t8 <- data.frame(SpatialDist=as.vector(as.dist(P8.dist)), GeneticDist=as.vector(as.dist(pvdist.8)))
t9 <- data.frame(SpatialDist=as.vector(as.dist(P9.dist)), GeneticDist=as.vector(as.dist(pvdist.9)))
t10 <- data.frame(SpatialDist=as.vector(as.dist(P10.dist)), GeneticDist=as.vector(as.dist(pvdist.10)))
t11 <- data.frame(SpatialDist=as.vector(as.dist(P11.dist)), GeneticDist=as.vector(as.dist(pvdist.11)))
t12 <- data.frame(SpatialDist=as.vector(as.dist(P12.dist)), GeneticDist=as.vector(as.dist(pvdist.12)))
t13 <- data.frame(SpatialDist=as.vector(as.dist(P13.dist)), GeneticDist=as.vector(as.dist(pvdist.13)))
t14 <- data.frame(SpatialDist=as.vector(as.dist(P14.dist)), GeneticDist=as.vector(as.dist(pvdist.14)))


## ----figspatgen, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE-----
par(mfrow=c(4,4))
plot(t1, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="1", side=3, line = 1)
plot(t2, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="2", side=3, line = 1)
plot(t3, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="3", side=3, line = 1)
plot(t5, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="5", side=3, line = 1)
plot(t6, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="6", side=3, line = 1)
plot(t7, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="7", side=3, line = 1)
plot(t8, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="8", side=3, line = 1)
plot(t9, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="9", side=3, line = 1)
plot(t10, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="10", side=3, line = 1)
plot(t11, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="11", side=3, line = 1)
plot(t12, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="12", side=3, line = 1)
plot(t13, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="13", side=3, line = 1)
plot(t14, xlim=c(0,20), ylim = c(0, 0.18))
mtext(text="14", side=3, line = 1)



## ----set-params---------------------------------------------------------------
Nmarkers <- seq(10,150,20)
nsets.pick <- 30


## ----create.matrices----------------------------------------------------------
samplist <-vector()
for(i in 1:length(Nmarkers)){
    namev <- paste("samp.",as.character(Nmarkers[i]),sep="")
    assign(namev,t(replicate(nsets.pick, sort(sample(186,Nmarkers[i],replace=FALSE)))))
    samplist[i] <- namev
}


## ----rowVar-fn----------------------------------------------------------------
rowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}


## ----poly-fn------------------------------------------------------------------
fPolymorphic.subs <- function(DataFrame, indices,nmark){
    df.tmp <- gendat[,c(seq(1:6),indices)]
    d.tmp <- dim(df.tmp)[2]
    return(nmark - sum(duplicated(t(DataFrame[,start.gendat:d.tmp]))))}


## ----Rarefaction1,cache=TRUE--------------------------------------------------
 #Effect on diversity statistics and estimates of diversity stats and of N genotypes for pop
 mtmp <-matrix(NA,nrow=length(samplist),ncol = 30 )
 vtmp <-vector()
 npolymorph.tmp <- vector()
 npoly <- matrix(NA,nrow=length(samplist),ncol = 30)
 ppr.sample <- list()
 count <- 0
 for(j in 1:length(samplist)){
    for(i in 1:30){
        count <- count + 1
        indices <- get(samplist[j])[i,] # +6
        vtmp[i] <- dim(unique(gendat[,indices]))[[1]]
        npolymorph.tmp[i] <- fPolymorphic.subs(gendat,indices,Nmarkers[j])
#        print(poppr(genalex2[,indices]))
        ppr.sample[[count]] <- poppr::poppr(genalex2[,indices])[,1:12]
    }
    mtmp[j,] <- vtmp
    npoly[j,] <- npolymorph.tmp
}


## ----DiversityStats-and-rarefaction,cache=TRUE--------------------------------
 imin <- vector();
 for (i in 1:8){
  imin[i] <- (30 + 1) + ((i -2) *30)
 }
 imax <- imin + 29

 # H
 H.mean <- matrix(nrow=8, ncol=14)
 H.var <- matrix(nrow=8, ncol=14)
 for(j in 1:length(samplist)){
    H <- lapply(ppr.sample[imin[j]:imax[j]], `[`, c('H'))
    A.H <- matrix(nrow=14, ncol=30, unlist(H[1:30]), byrow=FALSE)
    H.mean[j,] <- rowMeans(A.H)
    H.var[j,] <- rowVar(A.H)
}
 colnames(H.mean) <- c(Patches,"Total")
 colnames(H.var) <- c(Patches,"Total")
 rownames(H.mean) <- Nmarkers
 rownames(H.var) <- Nmarkers
 H.sd <- sqrt(H.var)
 H.sd[is.infinite(H.sd)] <- NA

 # G - should just be exp(H) but we'll calculate anyway in case there's some issue
 G.mean <- matrix(nrow=8, ncol=14)
 G.var <- matrix(nrow=8, ncol=14)
 for(j in 1:length(samplist)){
    G <- lapply(ppr.sample[imin[j]:imax[j]], `[`, c('G'))
    A.G <- matrix(nrow=14, ncol=30, unlist(G[1:30]), byrow=FALSE)
    G.mean[j,] <- rowMeans(A.G)
    G.var[j,] <- rowVar(A.G)
}
 colnames(G.mean) <- c(Patches,"Total")
 colnames(G.var) <- c(Patches,"Total")
 rownames(G.mean) <- Nmarkers
 rownames(G.var) <- Nmarkers
 G.sd <- sqrt(G.var)
 G.sd[is.infinite(G.sd)] <- NA

 # lambda
 lambda.mean <- matrix(nrow=8, ncol=14)
 lambda.var <- matrix(nrow=8, ncol=14)
 for(j in 1:length(samplist)){
    lambda <- lapply(ppr.sample[imin[j]:imax[j]], `[`, c('lambda'))
    A.lambda <- matrix(nrow=14, ncol=30, unlist(lambda[1:30]), byrow=FALSE)
    lambda.mean[j,] <- rowMeans(A.lambda)
    lambda.var[j,] <- rowVar(A.lambda)
}
 colnames(lambda.mean) <- c(Patches,"Total")
 colnames(lambda.var) <- c(Patches,"Total")
 rownames(lambda.mean) <- Nmarkers
 rownames(lambda.var) <- Nmarkers
 lambda.sd <- sqrt(lambda.var)
 lambda.sd[is.infinite(lambda.sd)] <- NA

 # Hexp
 Hexp.mean <- matrix(nrow=8, ncol=14)
 Hexp.var <- matrix(nrow=8, ncol=14)
 for(j in 1:length(samplist)){
    Hexp <- lapply(ppr.sample[imin[j]:imax[j]], `[`, c('Hexp'))
    A.Hexp <- matrix(nrow=14, ncol=30, unlist(Hexp[1:30]), byrow=FALSE)
    Hexp.mean[j,] <- rowMeans(A.Hexp)
    Hexp.var[j,] <- rowVar(A.Hexp)
}
 colnames(Hexp.mean) <- c(Patches,"Total")
 colnames(Hexp.var) <- c(Patches,"Total")
 rownames(Hexp.mean) <- Nmarkers
 rownames(Hexp.var) <- Nmarkers
 Hexp.sd <- sqrt(Hexp.var)
 Hexp.sd[is.infinite(Hexp.sd)] <- NA

 #Ia
 Ia.mean <- matrix(nrow=8, ncol=14)
 Ia.var <- matrix(nrow=8, ncol=14)
 A.Ia.missing <- matrix(nrow=8, ncol=14)
 for(j in 1:length(samplist)){
    Ia <- lapply(ppr.sample[imin[j]:imax[j]], `[`, c('Ia'))
    A.Ia <- matrix(nrow=14, ncol=30, unlist(Ia[1:30]), byrow=FALSE)
    A.Ia.missing[j,] <-  rowSums(is.nan(A.Ia))
    A.Ia[is.nan(A.Ia)] <- NA
    Ia.mean[j,] <- rowMeans(A.Ia, na.rm=TRUE)
    Ia.var[j,] <- rowVars(A.Ia, na.rm=TRUE)
}
 colnames(Ia.mean) <- c(Patches,"Total")
 colnames(Ia.var) <- c(Patches,"Total")
 rownames(Ia.mean) <- Nmarkers
 rownames(Ia.var) <- Nmarkers
 colnames(A.Ia.missing) <- c(Patches,"Total")
 rownames(A.Ia.missing) <- Nmarkers
 Ia.sd <- sqrt(Ia.var)
 Ia.sd[is.infinite(Ia.sd)] <- NA

 # rbarD
 rbarD.mean <- matrix(nrow=8, ncol=14)
 rbarD.var <- matrix(nrow=8, ncol=14)
 A.rbarD.missing <- matrix(nrow=8, ncol=14)
 for(j in 1:length(samplist)){
    rbarD <- lapply(ppr.sample[imin[j]:imax[j]], `[`, c('rbarD'))
    A.rbarD <- matrix(nrow=14, ncol=30, unlist(rbarD[1:30]), byrow=FALSE)
    A.rbarD.missing[j,] <-  rowSums(is.nan(A.rbarD))
    A.rbarD[is.nan(A.rbarD)] <- NA
    rbarD.mean[j,] <- rowMeans(A.rbarD, na.rm=TRUE)
    rbarD.var[j,] <- rowVars(A.rbarD, na.rm=TRUE)
}
 colnames(rbarD.mean) <- c(Patches,"Total")
 colnames(rbarD.var) <- c(Patches,"Total")
 rownames(rbarD.mean) <- Nmarkers
 rownames(rbarD.var) <- Nmarkers
 colnames(A.rbarD.missing) <- c(Patches,"Total")
 rownames(A.rbarD.missing) <- Nmarkers
 rbarD.sd <- sqrt(rbarD.var)
 rbarD.sd[is.infinite(rbarD.sd)] <- NA



## ----graphs.rbar--------------------------------------------------------------
 melt.rbar <- melt(rbarD.mean, id=Nmarkers)
 colnames(melt.rbar) <- c("Nmarkers","Patch","rbarD")
 p1 <- ggplot(melt.rbar,aes(x=Nmarkers,y=rbarD, group=Patch,
        colour=Patch)) + geom_point() + geom_line() +
        xlab("N markers") + ylab(expression(bar(r)[d]))
 melt.rbar.sd <- melt(rbarD.sd, id=Nmarkers)
 colnames(melt.rbar.sd) <- c("Nmarkers","Patch","rbarSD")
 p2 <- ggplot(melt.rbar.sd,aes(x=Nmarkers,y=rbarSD, group=Patch,
        colour=Patch)) + geom_point() + geom_line() +
        scale_y_log10() +
        theme(legend.position="none") + xlab("N markers") +
        ylab(expression(SD(bar(r)[d])))


## ----figrbar-rarefaction, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----
grid.arrange(ggplotGrob(p1), ggplotGrob(p2), widths=c(2,1.6))


## ----lmlnsd-------------------------------------------------------------------
    sd.mtmp <- sqrt(rowVar(mtmp))
    lnsd.mtmp <- log(sd.mtmp)
    lnsd.mtmp[is.infinite(lnsd.mtmp)] <- NA

    lm.lnsd <- lm(lnsd.mtmp ~ Nmarkers, na.action=na.exclude)
    summary(lm.lnsd)


## ----figmeansd-marker, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----
    par(mfrow=c(1,2))
    plot(rowMeans(mtmp) ~ Nmarkers, pch = 16, font.lab = 2, xlab="Number of markers", ylab="Mean N of unique genotypes")
    plot(log(sqrt(rowVar(mtmp))) ~ Nmarkers, pch = 16, font.lab = 2, xlab="Number of markers", ylab="log[SD(Mean N of unique genotypes)]")
    abline(lm.lnsd)


## ----lm.se--------------------------------------------------------------------
    lm.se <- lm(sqrt(rowVar(mtmp))/rowMeans(mtmp) ~ rowMeans(mtmp))
    summary(lm.se)


## ----figmeanse-marker, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----
    plot(sqrt(rowVar(mtmp))/rowMeans(mtmp) ~ rowMeans(mtmp), pch = 16, font.lab = 2, xlab="Mean N of unique genotypes", ylab="SE(Mean N of unique genotypes)")
    abline(lm.se)


## ----nmarkers.tot,results="asis", echo=FALSE----------------------------------
    Nmarkers.total <- cbind(as.integer(Nmarkers), rowMeans(mtmp),rowVar(mtmp))
    colnames(Nmarkers.total) <- c("N markers", "Mean","Sample variance")
    tab.nm.tot <- xtable(Nmarkers.total, digits = c(0,0,2,2), caption="Mean and sample variance of number of unique genotypes identified, by number of markers used.", label="tab:nm.tot")
    tab.nm.tot


## ----ByPatchStats-------------------------------------------------------------
patch.samp <- list()
mtmp <- matrix(NA,nrow=length(samplist),ncol = 30 )
vtmp <- vector()
for(k in Patches){
    tmp1 <- subset(gendat, gendat$Patch==k)
    for(j in 1:length(samplist)){
        for(i in 1:30){
            indices <- get(samplist[j])[i,] +6
            vtmp[i] <- dim(unique(tmp1[,indices]))[[1]]
            npolymorph.tmp[i] <- fPolymorphic.subs(tmp1,indices,Nmarkers[j])
        }
        mtmp[j,] <- vtmp
    patch.samp[[k]] <- mtmp
    }
}


## ----meanvar------------------------------------------------------------------
lmean <- lapply(1:13, function(i) rowMeans(patch.samp[[i]]))
lvar <-  lapply(1:13, function(i) rowVar(patch.samp[[i]]))

lmat <- unlist(lmean)
lmvar <- unlist(lvar)

lmat <- matrix(lmat,nrow=13, byrow=TRUE)
lvar <- matrix(lmvar,nrow=13, byrow=TRUE)


## ----meanvar-manip------------------------------------------------------------
lmat.t <- as.data.frame(lmat)
lmat.t$Patch <- Patches
rownames(lmat.t) <- c()
lmat.t$Patch <- as.factor(c(14,13,5,7,11,1,12,10,9,2,3,6,8))
colnames(lmat.t) <- c(Nmarkers,"Patch")
plot_data <- melt(lmat.t,id.var="Patch")
colnames(plot_data) <- c("Patch","markers","mean")

p1 <- ggplot(plot_data, aes(x = markers,
    y = mean,
    group = Patch,
    colour = Patch))   +
    geom_point() +
    geom_line(aes(linetype=Patch, color=Patch)) +
    labs(x ="Number of markers", y = "Mean N unique genotypes")


lvar.t <- as.data.frame(lvar)
lvar.t$Patch <- as.factor(c(14,13,5,7,11,1,12,10,9,2,3,6,8))
colnames(lvar.t) <- c(Nmarkers, "Patch")

plot_var <- melt(lvar.t,id.var="Patch")
colnames(plot_var) <- c("Patch","markers","var")

p2 <- ggplot(plot_var, aes(x = markers,
    y = var,
    group = Patch,
    colour = Patch))   +
    geom_point() +
    geom_line(aes(linetype=Patch, color=Patch)) +
    labs(x ="Number of markers", y = "Var(N unique genotypes)") +
    theme(legend.position = "none")


## ----figmeanvar-patch, fig.keep="high", echo=FALSE, include=TRUE, warning=FALSE, exec=FALSE----

grid.arrange(ggplotGrob(p1),ggplotGrob(p2), widths=c(2,1.6))
