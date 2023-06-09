##################
## Script for analyzing Arctic char data for LTER lakes
## for Klobucar et al. (2021)
## Script & analysis by Jessica Rick, jrick@uwyo.edu
##################

library(vcfR)
library(ggplot2)
library(dartR)
library(adegenet)
library(cowplot)

colors<-rainbow(13,start=0, end=0.9,alpha = 0.6)
colors2<-c("#8DD3C7","#BEBADA","#FB8072","#80B1D3",
           "#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
colors.mt<-c("#FFFFFF","#8DD3C7","#BEBADA","#FB8072","#80B1D3",
          "#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")

shapes<-c(0,1,2,3)

### ----- function to calc pca from genetic covariances --------- 
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



### ------------------------- 

#######################################
## Import data from VCF and metadata ##
#######################################

char_vcfR<-read.vcfR("variants_LTER_nolowcov10k_miss0.5_maf0.01_noSex.vcf")
chargen <- vcfR2genlight(char_vcfR)

fishinfo<-read.csv('AC_GBSdata_Master_Clean.csv', 
                   header=TRUE,stringsAsFactors = FALSE)

pairedinfoChar <- left_join(data.frame(New_ID=indNames(chargen)),
                            fishinfo,by="New_ID")

pop(chargen) <- pairedinfoChar$LAKE

chargen_LTER <- gl.keep.pop(chargen,
                            c("GTH57", "GTH58", "GTH59", "GTH60"))

## extracting genotype matrix
char_alleles <- t(as.matrix(chargen_LTER))

dim(char_alleles)
head(colnames(char_alleles)) ## to make sure that the names look good

#######################################################
### calculate missing data for these snps, per indv ###
### and filter indv with lots of missing data #########
#######################################################
missingness <- data.frame(New_ID = colnames(char_alleles),
                Missing = numeric(length(colnames(char_alleles))),
                Heterozygosity = numeric(length(colnames(char_alleles))))

for (i in 1:length(colnames(char_alleles))){
  missingness$Missing[i] <- round((length(which(is.na(char_alleles[,i]))))/length(rownames(char_alleles)),3)
  missingness$Heterozygosity[i] <- round((length(which((char_alleles[,i] == 1))))/(length((char_alleles[,i]))-length(which(is.na(char_alleles[,i])))),3)
}

par(mfrow=c(1,2))
hist(as.numeric(missingness$Heterozygosity), 
  main=paste("Heterozygosity, mean ",round(mean(as.numeric(missingness$Heterozygosity)),3),sep=""),
  xlab="Heterozygosity")
hist(as.numeric(missingness$Missing), 
  main=paste("Missingness, mean ",round(mean(as.numeric(missingness$Missing)),3),sep=""),
  xlab="% missing data")
plot(as.numeric(missingness$Missing),as.numeric(missingness$Heterozygosity),
  xlab="missingness", ylab="heterozygosity")
summary(as.numeric(missingness$Missing))
summary(as.numeric(missingness$Heterozygosity))

### check to see whether any indv have high levels of missing data, and remove them ###
lowcov<-missingness[which(missingness$Missing > 0.8),1]
lowcov
highhet <- missingness[missingness$Heterozygosity > 0.3,1]

to.remove <- match(c(lowcov,highhet),colnames(char_alleles))
char_alleles_nolowcov <- char_alleles[,-c(to.remove)]
attributes(char_alleles_nolowcov)$dim ## now have 129 indv

missingness.nolowcov <- missingness[-c(to.remove),]
plot(missingness.nolowcov$Missing, missingness.nolowcov$Heterozygosity,
  xlab="missingness",ylab="heterozyg")

chargen_nolowcov <- gl.filter.callrate(chargen,method="ind",threshold=0.2,recalc=F,mono.rm=F)

##################################
## Now, do PCA and plot results ##
##################################

char_pca <- do.pca(char_alleles_nolowcov)
pcSummary <- summary(char_pca)
scree <- plot(char_pca, type="lines") ## looks like first 5ish are interesting

pairedinfoChar<-left_join(as.data.frame(missingness.nolowcov),
                          fishinfo,by="New_ID",all.x=TRUE) 
head(pairedinfoChar)

### combining fish info with PC results ##

pcaAll <- data.frame(sample.id = pairedinfoChar$New_ID,
                     year = factor(pairedinfoChar$YEAR),
                     lake = factor(pairedinfoChar$LAKE),
                     length = as.numeric(pairedinfoChar$LENGTH_mm),
                     weight = as.numeric(pairedinfoChar$WEIGHT_g),
                     EV1 = char_pca$x[,1],    # the first eigenvector
                     EV2 = char_pca$x[,2],    # the second eigenvector
                     EV3 = char_pca$x[,3],    # the third eigenvector
                     EV4 = char_pca$x[,4],
                     EV5 = char_pca$x[,5],
                     stringsAsFactors = FALSE)

### plotting by lake ###

par(mfrow=c(1,4))

plot(pcaAll$EV1, pcaAll$EV2, pch=19, cex=1, lwd=1, col=scales::alpha(colors[pcaAll$lake],0.5),
     xlab=paste("PC1 (", round(pcSummary$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     main="LTER Lakes")
plot(pcaAll$EV2, pcaAll$EV3, pch=19,  cex=1, lwd=1, col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""))
plot(pcaAll$EV3, pcaAll$EV4, pch=19,  cex=1, lwd=1, col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""))
legend(x=0,y=0,legend=levels(pcaAll$lake),col=(colors),border=NULL,pch=19,
  bty="n", cex=1, pt.cex=2, pt.lwd=2, xpd=TRUE, horiz=FALSE)
plot(pcaAll$EV4, pcaAll$EV5, pch=19,  cex=1, lwd=1, col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""),
     ylab=paste("PC5 (", round(pcSummary$importance[2,5]*100, 1), "%)", sep=""))

### PC axis 1 vs. length
pcaAll$ecolength[pairedinfoChar$] <- 2
pcaAll$ecolength[pcaAll$length < 450] <- 1
pcaAll$ecolength <- as.factor(pcaAll$ecolength)

par(mfrow=c(1,3))
hist(pcaAll$length[pcaAll$lake == "GTH60" & pcaAll$EV1 < 0],breaks=30,
     xlab="length (mm)", main="Lake LTER348")
abline(v=450,lty=2,col="red")

boxplot(EV1 ~ ecolength, data=pcaAll[pcaAll$lake == "GTH60" & pcaAll$EV1 < 0,],
        ylim=c(min(pcaAll$EV1),0),main="Lake LTER348")
points(pcaAll$ecolength[pcaAll$lake == "GTH60"],pcaAll$EV1[pcaAll$lake == "GTH60"],col="red")
#t.test(EV1 ~ ecolength, data=pcaAll[pcaAll$EV1 > 0 & pcaAll$EV2 > 0,])
#leveneTest(EV1 ~ ecolength, data=pcaAll[pcaAll$EV1 > 0 & pcaAll$EV2 > 0,])

plot(pcaAll$EV1[pcaAll$lake == "GTH60" & pcaAll$EV1 < 0],
     pcaAll$length[pcaAll$lake == "GTH60" & pcaAll$EV1 < 0],
     ylab="length (mm)",
     xlab="PC Axis 1",
     main="Lake LTER348")

mod1<-lm(pcaAll$EV1[pcaAll$lake == "GTH60" & pcaAll$EV1 < 0] ~ pcaAll$length[pcaAll$lake == "GTH60" & pcaAll$EV1 < 0])
summary(mod1)

###################
## LTER 348 only ##
###################

LTER348 <- as.character(pcaAll$sample.id[pcaAll$lake == "GTH60"])

to.keep <- match(c(LTER348),colnames(char_alleles))
chargen_LTER348 <- gl.keep.ind(chargen,LTER348)
char_alleles_LTER348 <- char_alleles[,to.keep]

missingness.LTER348 <- missingness[c(to.keep),]
plot(missingness.LTER348[,2],missingness.LTER348[,3])

char_pca_LTER348 <- do.pca(char_alleles_LTER348)
pcSummary_LTER348 <- summary(char_pca_LTER348)
scree_LTER348 <- plot(char_pca_LTER348, type="lines") ## looks like first 4 are interesting

pairedinfoLTER348 <-merge(missingness.LTER348,fishinfo,
  by.x="New_ID",all.x=TRUE) 
head(pairedinfoLTER348)

### combining fish info with PC results ##

pcaAll_LTER348 <- data.frame(sample.id = pairedinfoLTER348$New_ID,
                     year = factor(pairedinfoLTER348$YEAR),
                     lake = factor(pairedinfoLTER348$LAKE),
                     length = as.numeric(pairedinfoLTER348$LENGTH_mm),
                     weight = as.numeric(pairedinfoLTER348$WEIGHT_g),
                     sizeclass = factor(pairedinfoLTER348$LENGTH_CLASS),
                     EV1 = char_pca_LTER348$x[,1],    # the first eigenvector
                     EV2 = char_pca_LTER348$x[,2],    # the second eigenvector
                     EV3 = char_pca_LTER348$x[,3],    # the third eigenvector
                     EV4 = char_pca_LTER348$x[,4],
                     EV5 = char_pca_LTER348$x[,5],
                     stringsAsFactors = FALSE)

pcaAll_LTER348$ecolength[pcaAll_LTER348$length < 440] <- 1
pcaAll_LTER348$ecolength[pcaAll_LTER348$length > 440] <- 2

colors3 <- c("olivedrab","darkmagenta")

#############################
par(mfrow=c(1,3))
hist(pcaAll_LTER348$length,breaks=20,main="Lake LTER348",xlab="length (mm)")
abline(v=440, col="red",lty=2)

boxplot(EV1 ~ sizeclass, data=pcaAll_LTER348,
        main="Lake LTER348",ylab="PC1")
points(pcaAll_LTER348$ecolength,
  pcaAll_LTER348$EV1,
  col=colors3[pcaAll_LTER348$ecolength])
#t.test(EV1 ~ ecolength, data=pcaAll_LTER)
#leveneTest(EV1 ~ ecolength, data=pcaAll[pcaAll$EV1 > 0 & pcaAll$EV2 > 0,])

plot(pcaAll_LTER348$EV1,
     pcaAll_LTER348$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_LTER348$importance[2,1]*100, 1), "%)", sep=""),
     main="Lake LTER348",
     col=colors3[pcaAll_LTER348$sizeclass],pch=19)
#text(x=pcaAll_LTER$EV1,y=pcaAll_LTER$length,pcaAll_LTER$sample.id,pos=4)

##############################

par(mfrow=c(1,3))
plot(pcaAll_LTER348$EV1, pcaAll_LTER348$EV2, pch=19, cex=1, lwd=1, col=(scales::alpha(colors3[pcaAll_LTER348$sizeclass],1)),
     xlab=paste("PC1 (", round(pcSummary_LTER348$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_LTER348$importance[2,2]*100, 1), "%)", sep=""))
legend(x="bottomleft",c("large","medium"),col=colors3,pch=19)

plot(pcaAll_LTER348$EV2, pcaAll_LTER348$EV3, pch=19, cex=1, lwd=1, col=(scales::alpha(colors3[pcaAll_LTER348$sizeclass],1)),
     xlab=paste("PC2 (", round(pcSummary_LTER348$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary_LTER348$importance[2,3]*100, 1), "%)", sep=""))
legend(x="bottomleft",c("large","medium"),col=colors3,pch=19)

plot(pcaAll_LTER348$EV3, pcaAll_LTER348$EV4, pch=19, cex=1, lwd=1, col=(scales::alpha(colors3[pcaAll_LTER348$sizeclass],1)),
     xlab=paste("PC3 (", round(pcSummary_LTER348$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary_LTER348$importance[2,4]*100, 1), "%)", sep=""))
legend(x="bottomleft",c("large","medium"),col=colors3,pch=19)

plot(pcaAll_LTER348$EV1,
     pcaAll_LTER348$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_LTER348$importance[2,1]*100, 1), "%)", sep=""),
     col=colors3[pcaAll_LTER348$sizeclass],pch=19)

#### ggplot
library(cowplot)

pca.lter348 <- ggplot(data = pcaAll_LTER348[!is.na(pcaAll_LTER348$sizeclass),], 
  aes(x=EV1,y=EV2,color=as.factor(sizeclass)))
pca.lter348.plot <- pca.lter348 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#006bb3","#4db8ff"),labels=c("large", "medium")) +
  scale_fill_discrete(guide=FALSE)+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        #legend.position=c(0.85,0.1),
        legend.position="none",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",fill=NA),
        text=element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.border=element_rect(color="black",fill=NA))+
  xlab(paste("\\nPC1 (", round(pcSummary_LTER348$importance[2,1]*100, 1), "%)", sep=""))+
  ylab(paste("PC2 (", round(pcSummary_LTER348$importance[2,2]*100, 1), "%)\\n", sep=""))
  #ylim(0,0.00025)

length.lter348 <- ggplot(data = pcaAll_LTER348[!is.na(pcaAll_LTER348$sizeclass),], 
  aes(x=EV1,y=length,color=as.factor(sizeclass)))
length.lter348.plot <- length.lter348 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#006bb3","#4db8ff","lightgray"),labels=c("large", "medium"),guide=FALSE) +
  scale_fill_discrete(guide=FALSE)+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        #legend.position=c(0.9,0.1),
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black",fill=NA),
        text=element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.border=element_rect(color="black",fill=NA))+
  xlab(paste("\\nPC1 (", round(pcSummary_LTER348$importance[2,1]*100, 1), "%)", sep=""))+
  ylab("Length (mm)\\n")
  #ylim(0,0.00025)

plot_grid(pca.lter348.plot, length.lter348.plot, labels = NULL)

###################
## LTER 347 only ##
###################

LTER347 <- as.character(pcaAll$sample.id[pcaAll$lake == "GTH59"])

to.keep <- match(c(LTER347),colnames(char_alleles))
char_alleles_LTER347 <- char_alleles[,c(to.keep)]

missingness.LTER347 <- missingness[c(to.keep),]
plot(missingness.LTER347[,2],missingness.LTER347[,3])

char_pca_LTER347 <- do.pca(char_alleles_LTER347)
pcSummary_LTER347 <- summary(char_pca_LTER347)
scree_LTER347 <- plot(char_pca_LTER347, type="lines") ## looks like first 4 are interesting

pairedinfoLTER347 <- merge(missingness.LTER347,fishinfo,by.x="New_ID",all.x=TRUE) 
head(pairedinfoLTER347)

### combining fish info with PC results ##

pcaAll_LTER347 <- data.frame(sample.id = pairedinfoLTER347$New_ID,
                     year = factor(pairedinfoLTER347$YEAR),
                     lake = factor(pairedinfoLTER347$LAKE),
                     length = as.numeric(pairedinfoLTER347$LENGTH_mm),
                     weight = as.numeric(pairedinfoLTER347$WEIGHT_g),
                     sizeclass = factor(pairedinfoLTER347$LENGTH_CLASS),
                     EV1 = char_pca_LTER347$x[,1],    # the first eigenvector
                     EV2 = char_pca_LTER347$x[,2],    # the second eigenvector
                     EV3 = char_pca_LTER347$x[,3],    # the third eigenvector
                     EV4 = char_pca_LTER347$x[,4],
                     EV5 = char_pca_LTER347$x[,5],
                     stringsAsFactors = FALSE)

pcaAll_LTER347$ecolength[pcaAll_LTER347$length < 440] <- 1
pcaAll_LTER347$ecolength[pcaAll_LTER347$length > 440] <- 2

#############################
par(mfrow=c(1,3))
hist(pcaAll_LTER347$length,breaks=20,
  main="Lake LTER347",xlab="length (mm)")
abline(v=440, col="red",lty=2)

plot(pcaAll_LTER347$EV1,
     pcaAll_LTER347$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_LTER$importance[2,1]*100, 1), "%)", sep=""),
     main="Lake LTER347",
     col=colors3[pcaAll_LTER347$sizeclass],pch=19)
#text(x=pcaAll_LTER$EV1,y=pcaAll_LTER$length,pcaAll_LTER$sample.id,pos=4)

##############################


#### ggplot
pca.lter347 <- ggplot(data = pcaAll_LTER347[!is.na(pcaAll_LTER347$sizeclass),], aes(x=EV1,y=EV2,color=as.factor(sizeclass)))
pca.lter347.plot <- pca.lter347 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#006bb3","#4db8ff"),labels=c("large", "medium")) +
  scale_fill_discrete(guide=FALSE)+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        #legend.position=c(0.85,0.15),
        legend.position="none",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",fill=NA),
        text=element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.border=element_rect(color="black",fill=NA))+
  xlab(paste("\\nPC1 (", round(pcSummary_LTER347$importance[2,1]*100, 1), "%)", sep=""))+
  ylab(paste("PC2 (", round(pcSummary_LTER347$importance[2,2]*100, 1), "%)\\n", sep=""))
  #ylim(0,0.00025)

pca.lter347.plot

###################
## LTER 345 only ##
###################

LTER345 <- as.character(pcaAll$sample.id[pcaAll$lake == "GTH57"])

to.keep <- match(c(LTER345), colnames(char_alleles))
char_alleles_LTER345 <- char_alleles[,c(to.keep)]

missingness.LTER345 <- missingness[c(to.keep),]
plot(missingness.LTER345[,2],missingness.LTER345[,3])

char_pca_LTER345 <- do.pca(char_alleles_LTER345)
pcSummary_LTER345 <- summary(char_pca_LTER345)
scree_LTER345 <- plot(char_pca_LTER345, type="lines") ## looks like first 4 are interesting

pairedinfoLTER345 <- merge(missingness.LTER345,fishinfo,by.x="New_ID",all.x=TRUE) 
head(pairedinfoLTER345)

### combining fish info with PC results ##

pcaAll_LTER <- data.frame(sample.id = pairedinfoLTER345$New_ID,
                     year = factor(pairedinfoLTER345$YEAR),
                     lake = factor(pairedinfoLTER345$LAKE),
                     length = as.numeric(pairedinfoLTER345$LENGTH_mm),
                     weight = as.numeric(pairedinfoLTER345$WEIGHT_g),
                     sizeclass = factor(pairedinfoLTER345$LENGTH_CLASS),
                     EV1 = char_pca_LTER345$x[,1],    # the first eigenvector
                     EV2 = char_pca_LTER345$x[,2],    # the second eigenvector
                     EV3 = char_pca_LTER345$x[,3],    # the third eigenvector
                     EV4 = char_pca_LTER345$x[,4],
                     EV5 = char_pca_LTER345$x[,5],
                     stringsAsFactors = FALSE)

pcaAll_LTER345$ecolength[pcaAll_LTER345$length < 440] <- 1
pcaAll_LTER345$ecolength[pcaAll_LTER345$length > 440] <- 2

#############################
par(mfrow=c(1,2))
hist(pcaAll_LTER345$length,breaks=20,
  main="Lake LTER345",xlab="length (mm)")
abline(v=440, col="red",lty=2)

plot(pcaAll_LTER345$EV1,
     pcaAll_LTER345$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_LTER345$importance[2,1]*100, 1), "%)", sep=""),
     main="Lake LTER345",
     col=colors3[pcaAll_LTER345$sizeclass],pch=19)
#text(x=pcaAll_LTER$EV1,y=pcaAll_LTER$length,pcaAll_LTER$sample.id,pos=4)

##############################


#### ggplot
pca.lter345 <- ggplot(data = pcaAll_LTER345[!is.na(pcaAll_LTER345$sizeclass),], 
  aes(x=EV1,y=EV2,color=as.factor(sizeclass)))
pca.lter345.plot <- pca.lter345 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#006bb3","#4db8ff"),labels=c("large", "medium")) +
  scale_fill_discrete(guide=FALSE)+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        #legend.position=c(0.15,0.9),
        legend.position="none",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",fill=NA),
        text=element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.border=element_rect(color="black",fill=NA))+
  xlab(paste("\\nPC1 (", round(pcSummary_LTER345$importance[2,1]*100, 1), "%)", sep=""))+
  ylab(paste("PC2 (", round(pcSummary_LTER345$importance[2,2]*100, 1), "%)\\n", sep=""))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  #ylim(0,0.00025)

pca.lter345.plot
