##################
## Script for analyzing Arctic char data
## for Klobucar et al. (2021)
## Script & analysis by Jessica Rick, jrick@uwyo.edu
##################

library(vcfR)
library(dartR)
library(dplyr)
library(ggplot2)
library(adegenet)
library(dartR)

colors<-rainbow(13,start=0, end=0.9,alpha = 0.6)

### ----- function to calculate genetic covariance matrix and pca --------- 
do.pca<-function(gmat){  ### inds in columns, loci in rows
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
  prcomp(~.,data=as.data.frame(gcovarmat),center=TRUE,scale=FALSE,na.action=na.exclude)
}

### ------ function to calculate DAPC threshold ------------------- 

randomize.dapc <- function(gl.obj,pop,npca,niter=10,verbose=TRUE){
  gl.obj$pop <- pop
  orig.mat <- as.matrix(gl.obj)
  rand.loadings <- c()
  
  for (i in 1:niter){
    if(verbose==TRUE){
      print(paste("running iteration ",i,sep=""))
    }
    rand.mat <- apply(orig.mat, 2, function(x) sample(x, replace = T))
    rand.gl <- as.genlight(rand.mat)
    rm(rand.mat)
    indNames(rand.gl) <- indNames(gl.obj)
    ploidy(rand.gl) <- 2
    rand.gl$pop <- gl.obj$pop
    
    # remove any NA loci
    toRemove <- is.na(glMean(rand.gl, alleleAsUnit = FALSE)) # TRUE where NA
    which(toRemove) # position of entirely non-typed loci
    glNoNA <- rand.gl[, !toRemove]
    
    dapc <- dapc(glNoNA,pop=pop(glNoNA),
                 n.pca=npca, n.da=3)
    loadings <- dapc$var.contr[,1]
    rm(rand.gl)
    #hist(loadings)
    #print(summary(loadings))
    
    rand.loadings <- append(rand.loadings,loadings)
  }
  #rand.hist <- hist(rand.loadings)
  #quantile(rand.loadings,c(0.95,0.975,0.99,1),na.rm=TRUE)
  return(quantile(rand.loadings,c(0.99),na.rm=TRUE))
}

###------------------


#######################################################
## Import and clean up data, select only focal lakes ##
#######################################################
char_vcfR<-read.vcfR("variants_nolowcov10k_miss0.5_maf0.01.recode.vcf")

chargen <- vcfR2genlight(char_vcfR)

fishinfo <- read.csv('AC_GBSdata_Master_Clean.csv', 
                   header=TRUE,stringsAsFactors = FALSE)
pairedinfoChar <- left_join(data.frame(New_ID=indNames(chargen)),
                            fishinfo,by="New_ID")

pop(chargen) <- pairedinfoChar$LAKE
chargen <- gl.compliance.check(chargen)

# select only individuals from focal lakes
lakes <- c("GTH60","GTH59","GTH58","GTH57",
           "FOG1","FOG2","FOG3","FOG5")
chargen_fog_lter <- gl.keep.pop(chargen, lakes)

## extracting genotype matrix
char_alleles <- t(as.matrix(chargen_fog_lter))

dim(char_alleles)
head(colnames(char_alleles)) ## to make sure that the names look good


#######################################################
### calculate missing data for these snps, per indv ###
#######################################################
missingness<-data.frame(New_ID = colnames(char_alleles),
  Missing = rep(0,times=length(colnames(char_alleles))),
  Heterozygosity = rep(0,times=length(colnames(char_alleles))))

for (i in 1:length(colnames(char_alleles))){
  missingness$Missing[i] <- round((length(which(is.na(char_alleles[,i]))))/length(rownames(char_alleles)),3)
  missingness$Heterozygosity[i] <- round((length(which((char_alleles[,i] == 1))))/(length((char_alleles[,i]))-length(which(is.na(char_alleles[,i])))),3)
}

par(mfrow=c(1,2))
hist(as.numeric(missingness$Heterozygosity), 
  main=paste("Heterozygosity, mean ",round(mean(as.numeric(missingness[,3])),3),sep=""),
  xlab="Heterozygosity")
hist(as.numeric(missingness$Missing), 
  main=paste("Missingness, mean ",round(mean(as.numeric(missingness[,2])),3),sep=""),
  xlab="% missing data")
plot(as.numeric(missingness$Missing),as.numeric(missingness$Heterozygosity),xlab="missingness",ylab="heterozygosity")
summary(as.numeric(missingness$Missing))
summary(as.numeric(missingness$Heterozygosity))

### check to see whether any indv have high levels of missing data ###
lowcov<-missingness[which(missingness[,2] > 0.8),1]
lowcov
highhet <- missingness[missingness[,3] > 0.3,1]
highhet 

to.remove <- match(c(lowcov,highhet),
                   colnames(char_alleles))
char_alleles_nolowcov <- char_alleles[,-c(to.remove)]
dim(char_alleles_nolowcov) 

missingness.nolowcov <- missingness[-c(to.remove),]
plot(missingness.nolowcov[,2], missingness.nolowcov[,3],
  xlab="missingness",ylab="heterozyg")


####################
## Now, PCA time! ##
####################

char_pca <- do.pca(char_alleles_nolowcov)
pcSummary <- summary(char_pca)
scree <- plot(char_pca, type="lines") 

### combining fish info with PC results ##

pairedinfoChar<-left_join(as.data.frame(missingness.nolowcov),
                          fishinfo,by="New_ID",all.x=TRUE) 

pcaAll <- data.frame(sample.id = pairedinfoChar$New_ID,
                     year = factor(pairedinfoChar$YEAR),
                     lake = factor(pairedinfoChar$LAKE),
                     length = factor(pairedinfoChar$LENGTH_mm),
                     weight = as.numeric(pairedinfoChar$WEIGHT_g),
                     EV1 = char_pca$x[,1],    # the first eigenvector
                     EV2 = char_pca$x[,2],    # the second eigenvector
                     EV3 = char_pca$x[,3],    # the third eigenvector
                     EV4 = char_pca$x[,4],
                     EV5 = char_pca$x[,5],
                     stringsAsFactors = FALSE)


### plotting pca, colored by lake ###

par(mfrow=c(1,4),mar=c(4,4,1,6),oma=c(1,1,1,6))
plot(pcaAll$EV1, pcaAll$EV2, pch=19, cex=1, lwd=1, 
     col=scales::alpha(colors[pcaAll$lake],0.5),
     xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""))

plot(pcaAll$EV2, pcaAll$EV3, pch=19,  cex=1, lwd=1, 
     col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""))

plot(pcaAll$EV3, pcaAll$EV4, pch=19,  cex=1, lwd=1, 
     col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""))

plot(pcaAll$EV4, pcaAll$EV5, pch=19,  cex=1, lwd=1, col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""),
     ylab=paste("PC5 (", round(pcSummary$importance[2,5]*100, 1), "%)", sep=""))

legend(x=0.4,y=0.20,legend=levels(pcaAll$lake),
       col=(colors),border=NULL,pch=19,
       bty="n", cex=1, pt.cex=2, pt.lwd=2, 
       xpd=TRUE, horiz=FALSE)


## plotting with phenotypically sexed individuals

sexed <- read.csv("../data/char_sex_data.csv", header=T)
sexed_lter_fog <- merge(pcaAll_fog_lter, sexed, 
  by.x="sample.id", by.y="NEW_ID",
  all.x=TRUE, all.y=FALSE)

par(mfrow=c(1,3))

par(mar=c(6,6,1,2),oma=c(1,1,1,6))
plot(sexed_lter_fog$EV1, sexed_lter_fog$EV2, pch=21, cex=2, lwd=1.5,bg=scales::alpha(lake_cols[sexed_lter_fog$lake],0.3),col=lake_cols[sexed_lter_fog$lake],
     xlab=paste("PC1 (", round(pcSummary_fog_lter$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_fog_lter$importance[2,2]*100, 1), "%)", sep=""),
     cex.lab=2)
points(sexed_lter_fog$EV1, sexed_lter_fog$EV2, col="black",  cex=2, lwd=2, pch=c(2,3)[sexed_lter_fog$SEX])

plot(sexed_lter_fog$EV2, sexed_lter_fog$EV3, pch=21,  cex=2, lwd=1, bg=scales::alpha(colors[sexed_lter_fog$lake],0.6),
     xlab=paste("PC2 (", round(pcSummary_fog_lter$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary_fog_lter$importance[2,3]*100, 1), "%)", sep=""),
     cex.lab=2)
points(sexed_lter_fog$EV2, sexed_lter_fog$EV3, col="black",  cex=2, lwd=2, pch=c(2,3)[sexed_lter_fog$SEX])

plot(sexed_lter_fog$EV3, sexed_lter_fog$EV4, pch=21,  cex=2, lwd=1, bg=scales::alpha(colors[sexed_lter_fog$lake],0.6),
     xlab=paste("PC3 (", round(pcSummary_fog_lter$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary_fog_lter$importance[2,4]*100, 1), "%)", sep=""),
     cex.lab=2)
legend(x=-0.001,y=0.003,legend=levels(sexed_lter_fog$lake),col=(colors),border=NULL,pch=19,bty="n", cex=1, pt.cex=2, pt.lwd=2, xpd=TRUE, horiz=FALSE)
points(sexed_lter_fog$EV3, sexed_lter_fog$EV4, col="black",  cex=2, lwd=2, pch=c(2,3)[sexed_lter_fog$SEX])
legend(x=0,y=0,legend=levels(sexed_lter_fog$SEX),col="black",pch=c(2,3)[sexed_lter_fog$SEX],bty="n", cex=1, pt.cex=2, pt.lwd=2)

plot(sexed_lter_fog$EV4, sexed_lter_fog$EV5, pch=21,  cex=2, lwd=1, bg=scales::alpha(colors[sexed_lter_fog$lake],0.6),
     xlab=paste("PC4 (", round(pcSummary_fog_lter$importance[2,4]*100, 1), "%)", sep=""),
     ylab=paste("PC5 (", round(pcSummary_fog_lter$importance[2,5]*100, 1), "%)", sep=""),
     cex.lab=2)
points(sexed_lter_fog$EV4, sexed_lter_fog$EV5, col="black",  cex=2, lwd=2, pch=c(2,3)[sexed_lter_fog$SEX])


sexed_lter_fog$genetic_sex <- case_when(sexed_lter_fog$EV5 > -0.005 ~ "male",
                                        sexed_lter_fog$EV5 < -0.005 ~ "female")


##################################
## pulling vcf data for sex DAPC
## removing low coverage and ambiguous individuals
##################################

chargen_lter_fog <- chargen_fog_lter[indNames(chargen_fog_lter) %in% sexed_lter_fog$sample.id]

ploidy(chargen_lter_fog) <- 2
pop(chargen_lter_fog) <- sexed_lter_fog$genetic_sex

### DAPC optimization for number of axes to retain
char.mat <- as.matrix(chargen_lter_fog)
char.mat.noNA <- gtools::na.replace(char.mat, mean, na.rm=T)
char.xval <- xvalDapc(char.mat.noNA,
                      grp=pop(chargen_lter_fog),
                      n.da=2,xval.plot=FALSE,training.set=0.9)
char.xval[2:6] # number of PCs to achieve lowest MSE and highest success

dapc.char <- dapc(chargen_lter_fog,
                  n.pca = as.integer(char.xval[6][[1]]), 
                  n.da = 2)

### DAPC randomization to determine loadings threshold
## it will throw an error if there are any NAs in pop
## (this will take a long while to run)

sex.threshold <- randomize.dapc(chargen_lter_fog,
                                 pop=pop(chargen_lter_fog),
                                 npca=as.integer(char.xval[6][[1]]),
                                 niter=100,
                                 verbose=TRUE)

scatter(dapc.char, scree.da=TRUE, bg="white", pch=20, cell=0, cstar=0, solid=.4, cex=2,clab=0, leg=TRUE)

loadingplot(dapc.char$var.contr[,1],cex.lab=0.5,srt=90,
            byfac=TRUE, threshold = sex.threshold)

sex.loci <- data.frame(dapc.char$var.contr[dapc1$var.contr[,1]>sex.threshold,])

sex.loci$SNP <- row.names(sex.loci)
sex.loci$locus <- chargen_lter_fog$loc.names[as.numeric(sex.loci$SNP)]
sex.loci$chrom <- as.character(chargen_lter_fog$chromosome[as.numeric(sex.loci$SNP)])

## plot of chromosomes with sex loci

par(pin=c(8,2),mar=c(4,2,1,1))
barplot(sort(table(sex.loci$chrom),decreasing=TRUE),
        las=2,
        xlab="",
        ylab="number of significant sex loci")

sex.chr <- unique(sex.loci$chrom)

write.csv(sex.loci,"char_sex_loci.csv",quote=FALSE,row.names=FALSE)

## heterozygosity for sex loci
## first, need to pull only sex loci

chargen_sex <- gl.keep.loc(chargen_lter_fog,sex.loci$locus)

chardiv<-gl.Ho(chargen_lter_fog)
chardiv_sex <- gl.Ho(chargen_sex)

plot(chardiv, xlab='Locus number',ylab='Observed heterozygosity')

## plotting heterozygosity at sex loci

plot(chardiv_sex, col="blue", pch=19)

charSexA <- gl.keep.pop(chargen_lter_fog, "male")
charSexB <- gl.keep.pop(chargen_lter_fog, "female")

chargenSexA_sex <- charSexA[,sex.loci$locus]
chargenSexB_sex <- charSexB[,sex.loci$locus]

chardiv_sexA <- gl.Ho(charSexA)
chardiv_sexB <- gl.Ho(charSexB)
chardiv_sexA_sexloci <- gl.Ho(chargenSexA_sex)
chardiv_sexB_sexloci <- gl.Ho(chargenSexB_sex)

par(mfrow=c(1,2))
plot(chardiv_sexA, xlab='Locus',ylab='Observed heterozygosity', main="All loci, colored by group",xaxt='n')
points(chardiv_sexB, xlab='Locus',ylab='Observed heterozygosity',col="turquoise")

plot(chardiv_sexA_sexloci, xlab='Locus',ylab='Observed heterozygosity',ylim=c(0,1), main="Sex Loci, colored by group",xaxt='n')
axis(1, at=1:length(chardiv_sexA_sexloci), labels=row.names(as.data.frame(chardiv_sexA_sexloci)), las=2)
points(chardiv_sexB_sexloci, xlab='Locus number',ylab='Observed heterozygosity',col="turquoise")
## pretty cool-- looks like most loci are heterozyg in one group and homozyg in the other!

# rm(chargen)
# rm(chargen_inv)
# rm(chargen_noambig) 
# rm(chargen_nolowcov) 
# rm(chargen_nolowcov_nounk) 
# rm(chargen_nosex) 
# rm(chargen_sex) 
# rm(chargenSexA_sex) 
# rm(chargenSexB_sex)