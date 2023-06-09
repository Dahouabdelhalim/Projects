#--------------------------------------------------------------------------------
#
# R SCRIPT for # Tyler, C. and Kowalewski, M., 2017, PRSB
#
# Last updated: February 21, 2017
#
# written by: M. Kowalewski and C. Tyler
#
# questions and comments shold be directed to kowalewski@ufl.edu
#--------------------------------------------------------------------------------

#------------ NOTE TO USERS -----------------------------------------------------
# This script and associated functions ('OBX Functions.R')
# are highly idiosyncratic and specifically designed for the datasets
# used in this study. The script may be adaptable to other datasets, but not
# painlessly. The documentation has been provided here so the interested reader
# can reproduce the results and explore other parametrizations 
# not reported in ths study.
#--------------------------------------------------------------------------------

rm(list=ls(all=TRUE))  	# remove all previously created objects (optional)
options(scipen=20)	# optional

#-------- libraries -------------------------------------------------------------
library(vegan)
library(vegetarian)
library(data.table)

#-------- custom functions, and datasets ----------------------------------------
source(file='OBX Functions.R')	# upload custom functions
all <- read.csv('Final_All_Taxa_BetaDiv.csv',head=TRUE)
mol <- read.csv('Final_Moll_BetaDiv.csv',head=TRUE)
nmol <- read.csv('Final_v2_All_ButMoll_BetaDiv.csv',head=TRUE)
rob <- read.csv('Final_RobustMoll_BetaDiv.csv',head=TRUE)
dead <- read.csv('Final_Dead.csv',head=TRUE)
dmol <- read.csv('Final_Dead_Moll_BetaDiv.csv',head=TRUE)
refs <- read.csv('Refs_local.csv')

# -------- color schemes for figures --------------------------------------------
# Figure 1
histcol <- c('black', 'grey30', 'gray40', 'gray50', 'gray60', 'gray70') 
exccol <- 'black'
astcol <- 'black'
pluscol <- 'black'
# Figure 2
ordcol <- c('grey40','orangered','wheat2','mediumseagreen','royalblue') # habitat colors
# Figure 4 
mycolline <- c('orangered3', 'dodgerblue3', 'dodgerblue4', 'chartreuse4', 'black', 'chocolate2') # dataset colors for Fig. 4G

#------ Assemble datasets into a single list and extract key factors ------------
dts <- list(LA=all, MLA=mol, NMLA=nmol, RMLA=rob, DA=dead, MDA=dmol)
dts <- sapply(dts, function(x) x[-39,-1]) 	#remove sample 39 and label column
refs <- refs[-39,]					# locality info (locality 39 excluded)
depth<- refs[,2] 						# water depth values
habitat <- as.factor(refs[,3])			# habitat type

#----- 1. Basic data summary -----------------------------------------------------
data.info <- function(x) c('number of taxa'=sum(colSums(x)>0),'number of specimens'=sum(x))
table.LA <- sapply(dts[1:4], data.info)
table.LA <- rbind(table.LA, 'percent taxa'=round(100*table.LA[1,]/table.LA[1,1],1),
                  'percent specimens'=round(100*table.LA[2,]/table.LA[2,1],1))
table.DA <- sapply(dts[5:6], data.info)
table.DA <- rbind(table.DA, 'percent taxa'=round(100*table.DA[1,]/table.DA[1,1],1),
                  'percent specimens'=round(100*table.DA[2,]/table.DA[2,1],1))
# table.S1 (not used - internal to this script)
(table.S1 <- cbind(table.LA, table.DA)[c(1,3,2,4),]) # reshuffle rows to sort by species and specimen

#--- Remove empty species column (also can remove singleton species)-------------
sapply(dts, function(x) sum(colSums(x)==0)) 			# any empty species columns?
dts2 <- sapply(dts, function(x) x[,-which(colSums(x)==0)]) 	# remove empty species columns
sapply(dts2, function(x) sum(colSums(x)==0)) 			# any empty species columns?
sapply(dts2, dim) 							# should match species in table.S1
dts3 <- sapply(dts, function(x) x[,which(colSums(x>0)>=2)]) # remove singleton species
sapply(dts3, dim) 							# likley fewer species than in table.S1

#--- 2. Beta variance (Table 1, Table 3, and Figure 1) ---------------------------

#### Literature Data
hs.ae<-c(2.5,2.6,1.4,3.7) 									#impacted ecosystems literature beta values
n.ae<-c(37,44,105,19) 										#total number of species
hs.nat<-c(3.9,2.8,2.4,2.1,1.8,2.1,1.4,4.6,3.9,3.9,3,1.3,1.2,3,2.4,2.2,1.9)	#natural ecosystems litertature beta values
n.nat<-c(43,125,45,35,38,41,81,65,86,85,85,54,87,66,65,78,64)			#total number of species
hab.ae<-c('S','S','R','C')#habitat S-sandy, C-coral, R-rocky reef, G-grass
hab.nat<-c('S','R','S','S','S','C','G','C','S','S','S','R','R','R','R','R','R')
class.ae<-c('P','P','M','P')#p=polychaete, M-mollusc
class.nat<-c('P','M','C','C','C','P','M','M','N','N','N','M','M','M','P','P','P')
blit <- c(hs.ae,hs.nat)
slit <- c(n.ae,n.nat)

# Check if the custom-written function {beta.var} is consistent with {beta.div} function of Legendre
# Legendre function downloaded from http://adn.biol.umontreal.ca/~numericalecology/Rcode/
source(file='Legendre_beta.div.function.R')	
beta.div(dts2[[1]],method='%difference', nperm=1)[[1]]		# raw counts
beta.var(dts2[[1]])[c(3,5)]							# raw counts
beta.div(decostand(dts2[[1]], 'total'),method='%difference', nperm=1)[[1]]	# rel. abundance
beta.var(decostand(dts2[[1]], 'total'))[c(3,5)] 			# rel. abundance
beta.div(wisconsin(dts2[[1]]),method='%difference', nperm=1)[[1]] # wisconsin
beta.var(wisconsin(dts2[[1]]))[c(3,5)] 					# wisconsin
beta.div(wisconsin(sqrt(dts2[[1]])),method='%difference', nperm=1)[[1]]	# sqrt-wisconsin
beta.var(wisconsin(sqrt(dts2[[1]])))[c(3,5)]				# sqrt-wisconsin

# Compute various beta diversity metrics
times <- 1000 								# set to 1000 for final version
loc.min <- apply(sapply(dts2, rowSums),1,min)   		# find minimum n for each locality across datasets
print(paste('NOTE: standardized datasets were all reduced to N =', sum(loc.min)))
B.V <- t(sapply(dts2, beta.var))					# Beta variance of the datasets
BVs1 <- sapply(dts2, myrar1, loc.min, times)			# Beta variance subsampled for non-standardized datasets
B.V.rs <- BVs1[seq(5,5*times,5),]					# Select subsampled Beta Variance values
B.V.std <- apply(BVs1[seq(5,5*times,5),],2,mean)		# Compute dataset means of Beta variance values
S.Beta <- sapply(dts2, H, lev='beta')				# Beta Shannon of the datasets
BSs1 <- sapply(dts2, myrar2, loc.min, times)			# Beta Shannon subsamped for standardized datasets
B.S.std <- apply(BSs1, 2, mean)					# Compute dataset means of Beta Shannon values
adj.sp <- replicate(times, sapply(dts2, find.sp, loc.min))	# Estimate sample sizes for standardized data

# Beta dispersion homogeneity test
dts2w <- sapply(dts2, decostand, 'total')				# rel-abundance dts2
dts4 <-rbindlist(dts2w,fill=TRUE)					# assemble datasets into a single matrix
dts4[is.na(dts4)]<-0							# fill in 0 values
dts2r <- sapply(dts2, rrarefy, loc.min)				# sample standardized dataset
dts2wr <- sapply(dts2r, decostand, 'total')			# rel-abundance dts2r
dts2wr <- sapply(dts2wr, data.frame)				# reformat
dts4r <- rbindlist(dts2wr, fill=T)  				# assemble datasets into a single matrix
dts4r[is.na(dts4r)]<-0							# fill in 0 values
dgroup <- factor(rep(names(dts2), each=51, times=1), levels=names(dts))
mod <- betadisper(vegdist(dts4, method='bray'), group=dgroup, type="median", bias.adjust=F)
mod.HSD <- TukeyHSD(mod)
table3A <- mod.HSD$group
mod2 <- betadisper(vegdist(dts4r, method='bray'), group=dgroup, type="median", bias.adjust=F)
mod.HSD2 <- TukeyHSD(mod2)
table3B <- mod.HSD2$group
beta.disp <- tapply(mod$distances, sort(dgroup), mean)
beta.disp2 <- tapply(mod2$distances, sort(dgroup), mean)
table1 <- cbind(B.V[,1:2], 'Shared S' = table.S1[2,], 'Shared N' = table.S1[4,],
                B.V[,c(3,5)], 'Beta Variance (subsampled)'=B.V.std, 'Shannon Beta'=S.Beta,
                'Shannon Beta (subsampled)'=B.S.std, 'Beta Dispersion'=beta.disp,
                'Beta Dispersion (subsampled)'=beta.disp2)
table1[,1:2] <- round(table1[,1:2],)
table1[,3:4] <- round(table1[,3:4],1)
table1[,3:4] <- round(table1[,3:4],1)
table1[,5:11] <- round(table1[,5:11],3)
table3 <- t(rbind('all data'=table3A[c(4,6,9),4],
                   'sample-standardized data'=table3B[c(4,6,9),4]))
#******************
# TABLES 1 AND 3
#******************
table1
table3
write.csv(t(table1), 'table1_12_12_16.csv')
write.csv(table3, 'table3_12_12_16.csv')


# Comparison of Beta Diversity Metrics (Table S6) (Fig. S6x)
cols2 <- c('red', 'orange', 'black', 'darkgray', 'blue', 'skyblue')
beta.mtr <- table1[,c(6:11)]
#******************
# TABLE S7
#******************
(tableS7 <- round(cor(beta.mtr),3))
write.csv(tableS7, 'tableS7_12_12_16.csv')
norm.mtr <- decostand(beta.mtr, 'standardize', MARGIN=2)
pool.mtr <- apply(norm.mtr, 1, mean)
#******************
# FIGURE S4
#******************
plot(pool.mtr, norm.mtr[,1], type='n', ylim=c(-2,2), las=1, 
     xlab='consensus beta diversity metrics',
     ylab='individual beta diversity metrics')
 abline(h=0, v=0, lty=3, col='gray')
 for (i in 1:6) {
     points(pool.mtr, norm.mtr[,i], pch=c(21), col=cols2[i])
     }
 text(pool.mtr, apply(norm.mtr,1,max), names(dts2), pos=3, cex=0.8)
 legend(-1.32, 2.2, pch=21, col=cols2, bty='n', cex=0.8,
       paste(colnames(beta.mtr),'r =', round(cor(pool.mtr, norm.mtr),3)))


#******************
# FIGURE 1
#******************
distr.d <- BSs1
start.bin <- 1
split.screen(c(5,5))
  for (i in 1:25) {
  sc <- screen(i)
  par(mar=c(0,0,0,0))
  if (i < 22) {
   box(col='black')
   mtext(side=3, line=-0.9, col='black', LETTERS[i], cex=0.8, adj=0.025)
   subsc <- split.screen(c(6,1), screen=sc)
   max.v <- max(blit[i],distr.d)+0.2
   for (j in 1:6) {
   screen(subsc[j])
   hist(distr.d[,j], breaks=seq(start.bin, max.v, 0.02), border=NA, col=histcol[j], main='', axes=F)
   abline(v=quantile(distr.d[,1], prob=0.975), col=histcol[1], lty=3, lwd=0.5)
   abline(v=quantile(distr.d[,1], prob=0.025), col=histcol[1], lty=3, lwd=0.5)
    if (j>1) {
     ifelse(quantile(distr.d[,1], prob=0.975) < blit[i], a1 <- 0, a1 <- 1)
     ifelse(quantile(distr.d[,j], prob=0.975) < blit[i], a2 <- 0, a2 <- 1)
     ifelse(quantile(distr.d[,1], prob=0.025) < blit[i], a3 <- 0, a3 <- 1)
     ifelse(quantile(distr.d[,j], prob=0.025) < blit[i], a4 <- 0, a4 <- 1)
     if (abs(a1-a2) + abs(a3-a4) == 0) mtext(side=1, '+', col=pluscol, cex=0.8, line=-0.9, adj=0.025)
     if (abs(a1-a2) + abs(a3-a4) != 0) mtext(side=1, '!', col=exccol, cex=0.8, line=-0.9, adj=0.025)
    }
   if (j==3) {points(cbind(blit[i], 50), cex=0.7, pch=8, col=astcol)}
   }
  }
 if (i==22) {
 par(mar=c(1,1,1,1))
 hist(distr.d[,1], breaks=seq(start.bin, max.v, 0.02), border=NA, col=histcol[1], main='', axes=F)
 mtext(side=1, 'Beta Shannon', cex=0.6, line=-0.2)
 mtext(side=4, 'number of iterations', cex=0.6, line=-1.2)
 par(mar=c(0,0,0,0))
 box(col='black')
 mtext(side=3, line=-0.9, col='black', LETTERS[i], cex=0.8, adj=-0.22)
 }

 if (i==23) {
 plot(1:10, 1:10, type='n', axes=F)
   mtext(side=3, line=-2, 'Beta Shannon', cex=0.6)
   mtext(side=3, line=-2.6, 'estimates in', cex=0.6)
   mtext(side=3, line=-3.2, 'other case studies', cex=0.6)
   points(5,5, pch=8, cex=0.7, col=astcol)
   box(col='black')
 mtext(side=3, line=-0.9, col='black', LETTERS[i], cex=0.8, adj=0.025)
 }

 if (i==24) {
  sc2 <- screen(i)
  box(col='black')
  mtext(side=3, line=-0.9, col='black', LETTERS[i], cex=0.8, adj=0.025)
  subsc <- split.screen(c(6,1), screen=sc2)
  max.v <- max(distr.d)+0.2
  for (j in 1:6) {
   screen(subsc[j])
   info.hist <- hist(distr.d[,j], breaks=seq(start.bin, max.v, 0.02), plot=F)
   hist(distr.d[,j], breaks=seq(start.bin, max.v, 0.02), border=NA, col=histcol[j], main='', axes=F)
   text(max(distr.d[,j])+0.1, 0.7*max(info.hist$counts), names(dts2[j]), col=histcol[j], cex=0.6) 
   }
 }
  if (i==25) {
 plot(1:10, 1:10, type='n', axes=F)
   text(5,8, 'MLA concordant 90.4%', col=histcol[2], cex=0.6)
   text(5,6.5, 'NMLA concordant 90.4%', col=histcol[3], cex=0.6)
   text(5,5, 'RMLA concordant 76.1%', col=histcol[4], cex=0.6)
   text(5,3.5, 'DA concordant 81.0%', col=histcol[5], cex=0.6)
   text(5,2, 'MDA concordant 90.4%', col=histcol[6], cex=0.6)
   box(col='black')
   mtext(side=3, line=-0.9, col='black', LETTERS[i], cex=0.8, adj=0.025)
 }
}
close.screen(all.screens = TRUE)


#------- 3. NMDS ORDINATIONS (FIGURE 2) ------------------------------------------
nmds.d <- dts2 # choose dataset to be analyzed (dts2 - with singletons, dts3 - without singletons)
kk <- 3 # change kk to 2 to perform NMDS in three dimensions
my.ord <- function(x, habitat, label, label2) {
 in.o <- metaMDS(x, k=kk)	# autotransform = (wisconsin(sqrt(data)))
 fin.o <- metaMDS(x, ,k=kk, previous.best=in.o)
 ordiplot(fin.o, type="n", xlim=c(-1.5,1.5), cex.lab=1, las=1, axes=F, xlab='', ylab='')
   points(fin.o,pch=c(25,22,23,24,21)[habitat],
         bg=ordcol[habitat],
                cex=1.4)
  ordihull(fin.o,groups=habitat,draw="polygon",label=F)
  axis(1, tck=-.02, padj=-0.7, cex.axis=0.8)
  axis(2, tck=-.02, hadj=0.7, las=1, cex.axis=0.8)
  box()
  mtext(label, side=3, adj=0.5, line=-1, cex=0.8)
  mtext(label2, side=3, line=-1.1, adj=0.01,cex=0.9)
  mtext(side=2, 'nMDS2', line=2, cex=0.8)
  mtext(side=1, 'nMDS1', line=2, cex=0.8)
  mtext(side=3, adj=0.95, line=-1, cex=0.7, paste('stress=',round(fin.o$stress,3),sep=''))
}

#******************
# FIGURE 2
#******************
op <- par(mfrow=c(2,2), pty='square', mar=c(4,4,0.2,0.2))
    my.ord(nmds.d[[2]], habitat, names(nmds.d)[2], LETTERS[1])
    my.ord(nmds.d[[3]], habitat, names(nmds.d)[3], LETTERS[2])
    my.ord(nmds.d[[1]], habitat, names(nmds.d)[1], LETTERS[3])
    my.ord(nmds.d[[5]], habitat, names(nmds.d)[5], LETTERS[4])
par(op)

#******************
# FIGURE S3
#******************
# live mollusks versus dead mollusks
op <- par(mfrow=c(1,2), pty='square', mar=c(1.5,1.5,0.2,0.2))
    my.ord(nmds.d[[2]], habitat, names(nmds.d)[2], LETTERS[5])
    my.ord(nmds.d[[6]], habitat, names(nmds.d)[6], LETTERS[6])
par(op)

#------- 4. PERMANOVA BY HABITAT(Table 2) ----------------------------------------
my.pm <- function(x) {
  adn <- adonis(x ~ habitat, permutations=9999, method='bray')
  return(c(adn[[1]][[4]][1], adn[[1]][[6]][1]))
}
out.pm <- sapply(dts2, my.pm)
#******************
# TABLE 2
#******************
(table2 <- data.frame('pseudo-F'=round(t(out.pm)[,1],2), 'p-value'=t(out.pm)[,2]))
write.csv(table2, 'table2_12_12_16.csv')

#--------- 5. HOMOGENEITY OF MULTIVARIATE DISPERSIONS ACROSS HABITATS (Figure 3) --
hab2 <- habitat[which(habitat!=1)]
hab2 <- droplevels(hab2)
dth <- sapply(dts2, function(x) x[which(habitat!=1),])
(observed.dist<-sapply(dth,ci.d, hab2))
rep.data <- replicate(1000,sapply(dth,ci.rep, hab2))
ci.data <- apply(rep.data,c(1,2),quantile,prob=c(0.025,0.975))
ci.data1 <- apply(rep.data,c(1,2),quantile,prob=c(0.005,0.995))
ci.data2 <- apply(rep.data,c(1,2),quantile,prob=c(0.05,0.95))
ci.data3 <- apply(rep.data,c(1,2),quantile,prob=c(0.25,0.75))
r.mean <- apply(rep.data,c(1,2),mean)
corr.lower<-as.matrix(ci.data[1,,]); corr.upper<-as.matrix(ci.data[2,,])
corr.lower1<-as.matrix(ci.data1[1,,]); corr.upper1<-as.matrix(ci.data1[2,,])
corr.lower2<-as.matrix(ci.data2[1,,]); corr.upper2<-as.matrix(ci.data2[2,,])
corr.lower3<-as.matrix(ci.data3[1,,]); corr.upper3<-as.matrix(ci.data3[2,,])
y.mean<-as.vector(observed.dist)
x.cor<-c(c(1:4),seq(1.1,4.1,1),seq(1.2,4.2,1),seq(1.3,4.3,1),seq(1.4,4.4,1), seq(1.5,4.5,1))
mycolsym <- rep(c('white','black','black','gray','white','black'),4)
y.mean2 <- y.mean[order(x.cor)]
x.cor2 <- sort(x.cor)

#******************
# FIGURE 3
#******************
op <- par(mgp=c(2.5,0.5,0.2),mar=c(4,4,1,0.3),mfrow=c(1,1))
plot(x.cor2,y.mean2,ylim=c(0.1,0.7),axes=F,ylab="average distance to centroid",xlab="",cex=0)
  a<-0
  for(i in 1:4){
   for(j in 1:6){
    a<-a+1  
    upper<-corr.upper[i,j]; lower<-corr.lower[i,j]
    upper1<-corr.upper1[i,j]; lower1<-corr.lower1[i,j]
    upper2<-corr.upper2[i,j]; lower2<-corr.lower2[i,j]
    upper3<-corr.upper3[i,j]; lower3<-corr.lower3[i,j]
    points(rep(x.cor2[a],2),c(lower1,upper1),type='l', lwd=1, col='lightsteelblue4', lend=2)
    points(rep(x.cor2[a],2),c(lower,upper),type='l', lwd=3, col='skyblue4', lend=2)
    points(rep(x.cor2[a],2),c(lower2,upper2),type='l', lwd=6, col='skyblue3', lend=2)
    points(rep(x.cor2[a],2),c(lower3,upper3),type='l', lwd=10, col='skyblue2', lend=2)
    points(x.cor2[a],y.mean2[a],pch=c(22,22,8,22,25,25)[j],bg=mycolsym[j])
    if (i==2 & j==5) {
    points(rep(4.2,2),c(lower1,upper1)-0.2,type='l', lwd=1, col='lightsteelblue4', lend=2)
    points(rep(4.2,2),c(lower,upper)-0.2,type='l', lwd=3, col='skyblue4', lend=2)
    points(rep(4.2,2),c(lower2,upper2)-0.2,type='l', lwd=6, col='skyblue3', lend=2)
    points(rep(4.2,2),c(lower3,upper3)-0.2,type='l', lwd=10, col='skyblue2', lend=2)
    text(4.2, 0.29, '50% CI', col='skyblue2', cex=0.7, pos=4)
    text(4.2, 0.22, '90% CI', col='skyblue3', cex=0.7, pos=4)
    text(4.2, 0.17, '95% CI', col='skyblue4', cex=0.7, pos=4)
    text(4.2, 0.12, '99% CI', col='lightsteelblue', cex=0.7, pos=4)
    }
   }
  }
  axis(1, tck=0,labels=F)
  axis(2, tck=-0.02, cex.axis=0.8,las=1, hadj=1.3)
  mtext(side=1,at=1.2,"harbor",line=1)
  mtext(side=1,at=2.2,"backsound",line=1)
  mtext(side=1,at=3.2,"nearshore",line=1)
  mtext(side=1,at=4.2,"offshore",line=1)
  legend(3.6,0.3,names(dth),pch=c(22,22,8,22,25,25), pt.bg=mycolsym,bty='n',cex=0.8)
par(op)

#------- 6. BETA GRADIENT ANALYSIS (FIGURE 4) ------------------------------------
#### BETA GRADIENT (FIGURE 4) (uses datasets in dts2 - singletons included)
# This loop may take ~5 minutes or more
out.bg <- list(mode='numeric', length=length(dts))
for (i in 1:length(dts2)) {
 out.bg[[i]] <- bgrad(dts2[[i]], depth, habitat)
 setTxtProgressBar(txtProgressBar(min = 0, max = length(dts2), style = 3), i)
 }
cat("\\n")	# run together with the for loop

mylwd = c(2,1,1,1,2,1)
mylty = c(1,1,2,3,1,1)
mycolline <- c('orange', 'green1', 'dodgerblue4', 'green1', 'black', 'darkgray') # dataset colors for Fig. 4G

##########
# FIGURE S6
##########
# proportion of variance in faunal similarity accounted for by difference in depth
# using various binning intervals
outvar <- vector(mode='numeric')
j <- 0
for (i in 1:6) {
 for (k in seq(0.1, 1, 0.1)) {
 j <- j + 1
 mybin <- k*round(out.bg[[i]][,4]/k)
 binBG <- tapply(out.bg[[i]][,3], mybin, mean)
 outvar[j] <- cor(as.numeric(names(binBG)), binBG)^2
 }
}
outvar2 <- matrix(outvar, 6, j/6, byrow=T)
plot(seq(0.1, 1, 0.1), outvar2[1,], type='n', ylim=c(0,1), xlab='bin width [m]',
     ylab='coefficient of determination', las=1)
    for (i in 1:6) {
    points(seq(0.1, 1, 0.1), outvar2[i,], type='o', pch=16, 
           lwd=mylwd[i], col=mycolline[i], lty=mylty[i])
    }
legend(0.8, 0.3, legend=names(dts2), col=mycolline, pch=16, lty=mylty, lwd=mylwd, bty='n')
 
#******************
# FIGURE 4
#******************
mybin <- 1		# choose bininng resolution [meters]
myct <- mean	# choose metric of central tendency
op <- par(mfrow=c(2,4), mar=c(2,3,1,0), oma=c(4,2.5,0,1))
for (i in (1:8)) {
  if (i<7) {
  b.g.plot(out.bg[[i]][,4], out.bg[[i]][,3], bin=mybin, ct=myct, mcex=2)
  mtext(side=3, adj=0.05, line=-1.5, LETTERS[i], cex=0.7)
  }
  if (i==7) {
  plot(out.bg[[1]][,4], out.bg[[1]][,3], type='n', xlim=c(0,15), ylim=c(0,0.4), xlab='', ylab='')
    for (k in 1:length(dts)) {
    x <- out.bg[[k]][,4]; y <- out.bg[[k]][,3] 
    bins <- round(1/mybin*x)/(1/mybin)
    l.y <- tapply(y,bins,myct)
    l.x <- as.numeric(names(l.y))
    points(l.x, l.y, type='l', col=mycolline[k], lty=mylty[k], lwd=mylwd[k])
    }
    legend(7.4,0.42, names(dts), col=mycolline, lty=mylty, lwd=mylwd, cex=0.8, bty='n')
    mtext(side=3, adj=0.05, line=-1.5, LETTERS[i], cex=0.7)
  }
  if (i==8) {
  plot(out.bg[[1]][,4], out.bg[[1]][,3], type='n', xlim=c(0,15), ylim=c(0,0.4), xlab='', ylab='')
    for (k in 1:length(dts)) {
    x <- out.bg[[k]][,4]; y <- out.bg[[k]][,3] 
    d.l <- decay(cbind(x,y))
    points(d.l[,1], d.l[,2], type='l', lty=mylty[k], lwd=mylwd[k], col=mycolline[k])
    }
    legend(7.4,0.42, names(dts), col=mycolline, lty=mylty, cex=0.8, bty='n', lwd=mylwd)
    mtext(side=3, adj=0.05, line=-1.5, LETTERS[i], cex=0.7)
 }
 if (i %in% c(1,5)) mtext(side=2, line=3, adj=0.5, 'Bray-Curtis similarity', cex=0.8)
 if (i %in% c(5:8)) mtext(side=1, line=2.5, adj=0.5, 'depth difference [m]', cex=0.8)
}
par(op)

#------- 7. EVENNESS - PIE (Figure S2, Table S4 and Table S5) ----------------------
#******************
# FIGURE S2
#******************
op <- par(mar=c(2.5,2.5,1.2,0.7), mgp=c(1.4,0.25,0), mfrow=c(3,1))
PIE <- function(x) (length(x)/(length(x)-1))*(1-sum((x/sum(x))**2))
tablealpha1 <- NULL;  tablealpha2 <- NULL; tablealpha3 <- NULL
myfilter <- 30                        # set filter for the final analysis
for (comp in 1:3) {
 if (comp==1) {i <- 2; j <- 3}
 if (comp==2) {i <- 1; j <- 5}
 if (comp==3) {i <- 2; j <- 6}
 for (filter in seq(30,100,5)) {                        # replace 'seq(30,100,5))' with '30)' to do analyses just for one filter
  (badsam <- unique(c(which(rowSums(dts2[[i]])<filter),which(rowSums(dts2[[j]])<filter))))
  LAdata <- dts2[[i]][-badsam,]
  DAdata <- dts2[[j]][-badsam,] # change to alld3[[2]] to evaluate live mollusks, 5 for DA
  smalln <- apply(cbind(rowSums(LAdata),rowSums(DAdata)),1, min)
  print(unique(c(which(rowSums(LAdata)<filter),which(rowSums(DAdata)<filter)))) # should return 'integer(0)'
  DPIE <- apply(DAdata,1,PIE)-apply(LAdata,1,PIE)
  DDIV <- log(rarefy(DAdata,filter,MARGIN=1))-log(rarefy(LAdata,filter,MARGIN=1))
  ymax <- round(max(max(DPIE),abs(min(DPIE))),1)
  xmax <- round(max(max(DDIV),abs(min(DDIV))),1)
  mysample <- function (x) mean(sample(x,replace=T))
  bootlist <- rep(list(DDIV),10000)
  bootlist2 <- rep(list(DPIE),10000)
  outDIV <- sapply(bootlist,mysample) 
  outPIE <- sapply(bootlist2,mysample)
  if (filter==myfilter) {DPIEo <- DPIE; DDIVo <- DDIV; outDIVo <- outDIV; outPIEo <- outPIE}
   if (comp==1) {
     tablealpha1 <- rbind(tablealpha1, c('filter'=filter, 'pDIV'=(sum(outDIV[-1]<=0)+1)/10000, 'pPIE'=(sum(outPIE[-1]<=0)+1)/10000,
                                    '% samples retained'=round(100*length(DPIE)/nrow(dts2[[i]]),1), 'mean DeltaDIV'=mean(DDIV), 'mean DeltaPIE'=mean(DPIE)))
     }
   if (comp==2) {
     tablealpha2 <- rbind(tablealpha2, c('filter'=filter, 'pDIV'=(sum(outDIV[-1]<=0)+1)/10000, 'pPIE'=(sum(outPIE[-1]<=0)+1)/10000,
                                    '% samples retained'=round(100*length(DPIE)/nrow(dts2[[i]]),1), 'mean DeltaDIV'=mean(DDIV), 'mean DeltaPIE'=mean(DPIE)))
     }
   if (comp==3) {
     tablealpha3 <- rbind(tablealpha3, c('filter'=filter, 'pDIV'=(sum(outDIV[-1]<=0)+1)/10000, 'pPIE'=(sum(outPIE[-1]<=0)+1)/10000,
                                    '% samples retained'=round(100*length(DPIE)/nrow(dts2[[i]]),1), 'mean DeltaDIV'=mean(DDIV), 'mean DeltaPIE'=mean(DPIE)))
     }
  if (filter==myfilter) {
    plot(DDIV, DPIE, cex=0, xlim=c(-1.5,1.5),ylim=c(-1,1),cex.axis=0.7,tck=-0.015,
         ylab=expression(Delta * ~ "PIE"),xlab=expression(Delta*~"ln"~"S"), las=1)
    abline(h=0,v=0, col='black', lwd=1, lty=3)
    points(DDIV, DPIE, cex=log(smalln,10), col='black', pch=21)
    }
  if (filter!=myfilter) {
    points(cbind(c(min(outDIV),max(outDIV)),c(mean(DPIE),mean(DPIE))),type='l',lwd=0.2, col='blue')
    points(cbind(c(mean(DDIV),mean(DDIV)),c(min(outPIE),max(outPIE))),type='l',lwd=0.2, col='blue')
    }
  if (filter==100) {
    points(cbind(c(min(outDIVo),max(outDIVo)),c(mean(DPIEo),mean(DPIEo))),type='l',lwd=3, col='white')
    points(cbind(c(min(outDIVo),max(outDIVo)),c(mean(DPIEo),mean(DPIEo))),type='l',lwd=1, col='red')
    points(cbind(c(mean(DDIVo),mean(DDIVo)),c(min(outPIEo),max(outPIEo))),type='l',lwd=3, col='white')
    points(cbind(c(mean(DDIVo),mean(DDIVo)),c(min(outPIEo),max(outPIEo))),type='l',lwd=1, col='red')
    points(mean(DDIVo), mean(DPIEo), cex=1, col='red', pch=16)
  }
 }
 #  mtext(paste('Comparison of', names(dts2)[i], 'and',  names(dts2)[j]),cex=0.6, adj=0.95, line=-1.2)
    mtext(LETTERS[comp], adj=0.02, cex=0.8, line=-1.2)
    mtext(paste('PIE', names(dts2)[i], '<', names(dts2)[j]),cex=0.6, adj=0.75,line=-1.2)   
    mtext(paste('PIE', names(dts2)[i], '<', names(dts2)[j]),cex=0.6, adj=0.25,line=-1.2)
    mtext(paste('   S', names(dts2)[i], '<', names(dts2)[j]),cex=0.6, adj=0.75,line=-2.2)
    mtext(paste('   S', names(dts2)[i], '>', names(dts2)[j]),cex=0.6, adj=0.25,line=-2.2)
    mtext(paste('PIE', names(dts2)[i], '>', names(dts2)[j]),cex=0.6, adj=0.75,line=-11.7)
    mtext(paste('PIE', names(dts2)[i], '>', names(dts2)[j]),cex=0.6, adj=0.25, line=-11.7)
    mtext(paste('   S', names(dts2)[i], '<', names(dts2)[j]),cex=0.6, adj=0.75, line=-12.7)
    mtext(paste('   S', names(dts2)[i], '>', names(dts2)[j]),cex=0.6, adj=0.25, line=-12.7)
}
par(op)
#******************
# TABLES S3-S5
#******************
(tableS3 <- tablealpha1)
(tableS4 <- tablealpha2)
(tableS5 <- tablealpha3)
write.csv(tableS3, 'tableS3_12_12_16.csv')
write.csv(tableS4, 'tableS4_12_12_16.csv')
write.csv(tableS5, 'tableS5_12_12_16.csv')

#------------ 8. PAIRWISE COMPARISON OF SIMILARITY (Supplementary Figure 5) ------
#******************
# FIG S5
#******************
myplot <- function(x,y,z1,z2,r) {
  plot(x, y, cex=1, main='', xlab=z1, ylab=z2,cex.axis=1, las=1, cex.lab=1.2)
  mtext(side=3, adj=0.10, line=-1.2, paste('rho =',r), cex=0.7)
  mtext(side=3, adj=0.10, line=-2.2, 'p < 0.0001', cex=0.7)
  mtext(side=3, line=-1.35, col='black', LETTERS[i], cex=1, adj=0.018)
}
op <- par(mfrow=c(1,3), mar=c(3,3,1.5,0.65), mgp=c(1.85,0.5,0), tck=-0.01)
ax.lab <- paste('ln(Bray-Curtis)',names(dts2))
for (i in 1:3) {
  if (i==1) {a=3; b=2}	# NMLA vs. MLA
  if (i==2) {a=1; b=5}	# LA vs DA
  if (i==3) {a=2; b=6}	# MLA cs. MDA
 x1 <- out.bg[[a]][,3]
 y1 <- out.bg[[b]][,3]
 rem0 <- unique(which(x1==0), which(y1==0))
 x <- log(x1[-rem0])
 y <- log(y1[-rem0])
 z1 <- ax.lab[a]
 z2 <- ax.lab[b]
 r <- round(cor(x,y,method='spearman'),3)
 myplot(x,y,z1,z2,r)
}
par(op)

#--------- SUMMARY STATISTIC FOR PAIRWISE COMPARISONS (Supplemental Table 6) ----
tableS6 <- t(sapply(out.bg, function(x) c(minimum=min(x[,3]), maximum=max(x[,3]), mean=mean(x[,3]), 'standard deviation'=sd(x[,3]))))
rownames(tableS6) <- names(dts2)
tableS6 <- round(tableS6, 2)
#******************
# TABLE S6
#******************
tableS6
write.csv(tableS6, 'tableS6_12_12_16.csv')

#----------------------------------- END MAIN SCRIPT ----------------------------

##### Additional figures and tables (not published)

# relative differences between dataset estimates
diffMLNL <- round(abs(table1.final[3,6:10]-table1.final[2,6:10])/apply(table1.final[2:3,6:10],2,max),2)
diffLD <- round(abs(table1.final[1,6:10]-table1.final[5,6:10])/apply(table1.final[c(1,5),6:10],2,max),2)
diffMLMD <- round(abs(table1.final[2,6:10]-table1.final[6,6:10])/apply(table1.final[c(2,6),6:10],2,max),2)
(table.rel.diff <- rbind('Mollusk LA vs. Non-mollusk LA'=diffMLNL, 'LA vs. DA'=diffLD, 'MLA vs. MDA'=diffMLMD))

# additional beta diversity plots (not included in the paper)
plot(slit, blit, type='n', ylim=c(0,5), xlim=c(20,205), ylab='Beta Shannon', xlab='Total Number of Species')
     abline(h=c(max(table1[,8]),min(table1[,8])), lwd=0.8, lty=3, col='blue')
     abline(h=c(max(table1[,9]),min(table1[,9])), lwd=0.8, lty=3, col='red')
     abline(h=c(max(blit), min(blit)), lwd=0.8, lty=3, col='black')
     points(slit, blit, pch=16, col='black', cex=0.7)
     text(sapply(dts2,dim)[2,], table1[,8], names(dts2), cex=0.7, col='blue')
     text(apply(adj.sp,1,mean), table1[,9], names(dts2), cex=0.7, col='red')
     points(cbind(c(205,205),c(min(table1[,8]),max(table1[,8]))),type='l', lwd='3', col='blue', lend=2)
     points(cbind(c(200,200),c(min(table1[,9]),max(table1[,9]))),type='l', lwd='3', col='red', lend=2)
     points(cbind(c(195,195),c(min(blit),max(blit))),type='l', lwd='3',col='black', lend=2)
     text(205,1.6, '(c)', srt=90, cex=0.75, adj=0.13, col='blue')  #range of published estimates
     text(200,1.9, '(b)', srt=90, cex=0.75, adj=0.06, col='red')  #observed range (adjusted)
     text(195,4.7, '(a)', srt=90, cex=0.75, adj=0, col='black')  #observed range
plot(slit, blit, type='n', ylim=c(0,5), xlim=c(20,205), ylab='Beta Shannon', xlab='Total Number of Species')
     abline(h=c(max(table1[,9]),min(table1[,9])), lwd=0.8, lty=3, col='red')
     abline(h=c(max(blit), min(blit)), lwd=0.8, lty=3, col='black')
     points(slit, blit, pch=16, col='black', cex=0.7)
     points(apply(adj.sp,1,mean), table1[,9], pch=21, cex=1, col='red')
     text(apply(adj.sp,1,mean), table1[,9], names(dts2), cex=0.7, col='red', pos=c(3,3,3,1,3), offset=0.4)
     points(cbind(c(200,200),c(min(table1[,9]),max(table1[,9]))),type='l', lwd='3', col='red', lend=2)
     points(cbind(c(195,195),c(min(blit),max(blit))),type='l', lwd='3',col='black', lend=2)
     text(200,1.9, '(b)', srt=90, cex=0.75, adj=0.06, col='red')  #observed range (adjusted)
     text(195,4.7, '(a)', srt=90, cex=0.75, adj=0, col='black')  #observed range

op <- par(oma=c(0,2,0,0))
plot(mod.HSD, las=1)
par(op)

#----------------------------------- END ----------------------------------------
