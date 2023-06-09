#Script to identify loci with excess heterozygotes to remove
# Script created by Monia Hasselhort and C.A. Buerkle, modified by M.E.F. LaCava


#### Identify loci with heterozygosity above a given percentage ####
genest <- read.table("pntest_mean_398_sorted.txt", header=F)
tmp.het.count <- apply(genest, 1, function(x) sum(abs(x-1)< 0.05, na.rm=T))
tmp.af<-apply(genest, 1, mean, na.rm=T)/2
#plot heterozygote count against allele frequency 
# - peak at center are loci with 0.5 allele frequencing in almost all individuals (unlikely to be biologically real)
plot(tmp.af,tmp.het.count)
#plot 
tmp<-apply(genest, 1, function(x) length(which(x==1)))
plot(tmp.het.count, tmp, xlim=c(0,400), ylim=c(0,400))

#count number of loci with heterozygosity in more than X% of individuals
length(which(tmp.het.count > 0.95*398))

#remove improbable loci from dataset
genest.subset <- genest[-which(tmp.het.count > 0.95*398),] #32 loci that are hetero in >95% of samples

#correlation between genest and genest.subset covariance matrices
g.covar<-getcovarmat(genest)
gsubset.covar<-getcovarmat(genest.subset)
cor(g.covar[lower.tri(g.covar)], gsubset.covar[lower.tri(gsubset.covar)])






#### Probabilistically identify loci with excess heterozygosity ####

#import data:
myvcf<-read.table("ph398_maf0.01_miss0.5_ind0.6_dp3_BayesNoOutliers_GT_DP.txt", stringsAsFactors=F)
## this looks like this:
## scaffold_11:53 APH_EE_001 0/0 99,0
## scaffold_11:53 APH_EE_002 0/0 49,0
## scaffold_11:53 APH_EE_003 0/0 99,0
## scaffold_11:53 APH_EE_005 0/0 99,0
## scaffold_11:53 APH_EE_007 0/0 99,0
## scaffold_11:53 APH_EE_009 0/0 59,0
## scaffold_11:53 APH_EE_010 0/0 100,0
## scaffold_11:53 APH_EE_013 0/0 99,0
## scaffold_11:53 APH_EE_015 0/0 63,0
## scaffold_11:53 APH_EE_017 0/0 96,0


AD<-matrix(as.numeric(unlist(strsplit(myvcf$V4, ","))), ncol=2, byrow=T)
ADhetprob<-apply(AD[myvcf$V3 == "0/1",],1, function(x) pbinom(x[1], prob=0.5, size=sum(x)))

cor(ADhetprob, rowSums(AD[myvcf$V3=="0/1",]))
#theirs: [1] 0.002717052
#mine: [1] -0.03793837

extreme_ADhetprob_by_contig<-by(ADhetprob, myvcf$V1[myvcf$V3 == "0/1"], function(x)sum(x>0.95 | x<0.05)/length(x))
mean_ADhetprob_by_contig<-by(ADhetprob, myvcf$V1[myvcf$V3 == "0/1"], mean)
N_ADhetprob_by_contig<-by(ADhetprob, myvcf$V1[myvcf$V3 == "0/1"], length)

plot(mean_ADhetprob_by_contig ~ by(rowSums(AD[myvcf$V3=="0/1",]), myvcf$V1[myvcf$V3 == "0/1"], mean))


##  length(ADhetprob)/nrow(myvcf)
## [1] 0.2529876

#subset genotype matrix according to above calc
nind <- 398
nloci<- 4949
g20<-matrix(scan("pntest_mean_ph398_maf0.01_miss0.5_ind0.6_dp3_BayesNoOutliers.recode.txt",n=nind*nloci,sep=" "),nrow=nloci,ncol=nind,byrow=T)
g20.subset<-g20[ mean_ADhetprob_by_contig < 0.9 & mean_ADhetprob_by_contig > 0.1,]
# 3947 loci retained

pc20<-do.pca(g20.subset)
pc20Summary<-summary(pc20)

par(mfrow=c(1,1), pty='s')
plot(pc20$x[,'PC1'], pc20$x[,'PC2'],col=inds$HRColor ,  pch=19 , cex=0.8,
     xlab =  paste("PC1 (", round(pc20Summary$importance[2,1]*100, 1), "%)", sep=""),
     ylab =  paste("PC2 (", round(pc20Summary$importance[2,2]*100, 1), "%)", sep=""), 
     main="PCA 3947loci in 398 inds")
#this PCA looks very similar to my original PCA - actually looks more clumped, when I would have expected it to 
# look more spread out. but PC importance increased slightly for both PC1 and PC2
plot(pc20$x[,'PC2'], pc20$x[,'PC3'],col=inds$HRColor ,  pch=19 , cex=0.8,
     xlab =  paste("PC2 (", round(pc20Summary$importance[2,2]*100, 1), "%)", sep=""),
     ylab =  paste("PC3 (", round(pc20Summary$importance[2,3]*100, 1), "%)", sep=""))



#check correlation between original covariance matrix and subset matrix
getcovarmat<-function(gmat, write.gcov=FALSE, inds=""){
  gmn<-apply(gmat,1,mean, na.rm=T)
  gmnmat<-matrix(gmn,nrow=nrow(gmat),ncol=ncol(gmat))
  gprime<-gmat-gmnmat ## remove mean
  
  gcovarmat<-matrix(NA,nrow=ncol(gmat),ncol=ncol(gmat))
  for(i in 1:ncol(gmat)){
    for(j in i:ncol(gmat)){
      if (i==j){
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
      }
      else{
        gcovarmat[i,j]<-cov(gprime[,i],gprime[,j], use="pairwise.complete.obs")
        gcovarmat[j,i]<-gcovarmat[i,j]
      }
    }
  }
  gcovarmat	
}

g20.covar<-getcovarmat(g20)
g20subset.covar<-getcovarmat(g20.subset)
cor( g20.covar[lower.tri(g20.covar)], g20subset.covar[lower.tri(g20subset.covar)])

