##################
## Script for analyzing Arctic char data for FOG lakes
## for Klobucar et al. (2021)
## Script & analysis by Jessica Rick, jrick@uwyo.edu
##################

library(vcfR)
library(dartR)
library(tidyverse)
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

##################################
## Import VCF data and metadata ##
##################################
char_vcfR <- read.vcfR("variants_FOG_nolowcov10k_miss0.5_maf0.01_noSex.recode.vcf")
char_gt <- extract.gt(char_vcfR, element = "GT", mask = FALSE, as.numeric = FALSE,
           return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE,
           convertNA = FALSE)
summary_gt <- summary(char_gt)

chargen <- vcfR2genlight(char_vcfR)

fishinfo<-read.csv('AC_GBSdata_Master_Clean.csv', 
                   header=TRUE,stringsAsFactors = FALSE)
pairedinfoChar <- left_join(data.frame(New_ID=indNames(chargen)),
  fishinfo,by="New_ID")

pop(chargen) <- pairedinfoChar$LAKE

## extracting genotype matrix
char_alleles <- t(as.matrix(chargen))

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
scree <- plot(char_pca, type="lines") ## looks like first 1 are interesting

pairedinfoChar<-left_join(as.data.frame(missingness.nolowcov),
                          fishinfo,
                          by="New_ID",all.x=TRUE) 
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
     ylab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""))

plot(pcaAll$EV2, pcaAll$EV3, pch=19,  cex=1, lwd=1, col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC2 (", round(pcSummary$importance[2,2]*100, 1), "%)", sep=""),
     ylab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""))

plot(pcaAll$EV3, pcaAll$EV4, pch=19,  cex=1, lwd=1, col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC3 (", round(pcSummary$importance[2,3]*100, 1), "%)", sep=""),
     ylab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""))
legend(x=0.4,y=0.20,legend=levels(pcaAll$lake),col=(colors),border=NULL,pch=19,bty="n", cex=1, pt.cex=2, pt.lwd=2, xpd=TRUE, horiz=FALSE)

plot(pcaAll$EV4, pcaAll$EV5, pch=19,  cex=1, lwd=1, col=scales::alpha(colors[pcaAll$lake],0.6),
     xlab=paste("PC4 (", round(pcSummary$importance[2,4]*100, 1), "%)", sep=""),
     ylab=paste("PC5 (", round(pcSummary$importance[2,5]*100, 1), "%)", sep=""))

### PC axis 1 vs. length
pcaAll$ecolength[pcaAll$length > 225] <- 2
pcaAll$ecolength[pcaAll$length < 225] <- 1
pcaAll$ecolength <- as.factor(pcaAll$ecolength)

boxplot(EV1 ~ ecolength, data=pcaAll[pcaAll$lake == "FOG3",],
        main="Lake FOG3",ylim=c(0,max(pcaAll$EV1)))
points(pcaAll$ecolength[pcaAll$lake == "FOG3"],pcaAll$EV1[pcaAll$lake == "FOG3"],col="red")
t.test(EV1 ~ ecolength, data=pcaAll[pcaAll$lake == "FOG3",])
#leveneTest(EV1 ~ ecolength, data=pcaAll[pcaAll$EV1 > 0 & pcaAll$EV2 > 0,])

plot(pcaAll$EV1[pcaAll$lake == "FOG3"],
     pcaAll$length[pcaAll$lake == "FOG3"],
     ylab="length (mm)",
     xlab="PC Axis 1",
     main="Lake FOG3")
#abline(h=225,col="red")

mod1<-lm(pcaAll$EV1[pcaAll$lake == "FOG3"] ~ pcaAll$length[pcaAll$lake == "FOG3"])
summary(mod1)

#####################
## Lake FOG 3 only ##
#####################

FOG3 <- as.character(pcaAll$sample.id[pcaAll$lake == "FOG3"])

to.keep <- match(c(FOG3),colnames(char_alleles))
char_alleles_FOG3 <- char_alleles[,c(to.keep)]

missingness.FOG3 <- missingness[c(to.keep),]

char_pca_FOG3 <- do.pca(char_alleles_FOG3)
pcSummary_FOG3 <- summary(char_pca_FOG3)
scree_FOG3 <- plot(char_pca_FOG3, type="lines") ## looks like first 1 are interesting

pairedinfoFOG3 <- merge(missingness.FOG3,fishinfo,
  by.x="New_ID",all.x=TRUE,sort=FALSE) 
head(pairedinfoFOG3)

### combining fish info with PC results ##

pcaAll_FOG3 <- data.frame(sample.id = pairedinfoFOG3$New_ID,
                     year = factor(pairedinfoFOG3$YEAR),
                     lake = factor(pairedinfoFOG3$LAKE),
                     length = as.numeric(pairedinfoFOG3$LENGTH_mm),
                     weight = as.numeric(pairedinfoFOG3$WEIGHT_g),
                     lib = factor(pairedinfoFOG3$Library),
                     missing = as.numeric(pairedinfoFOG3$Missing),
                     heterozyg = as.numeric(as.character(pairedinfoFOG3$Heterozygosity)),
                     sizeclass = as.factor(pairedinfoFOG3$LENGTH_CLASS),
                     EV1 = char_pca_FOG3$x[,1],    # the first eigenvector
                     EV2 = char_pca_FOG3$x[,2],    # the second eigenvector
                     EV3 = char_pca_FOG3$x[,3],    # the third eigenvector
                     EV4 = char_pca_FOG3$x[,4],
                     EV5 = char_pca_FOG3$x[,5],
                     stringsAsFactors = FALSE)

pcaAll_FOG3$ecolength[pcaAll_FOG3$length < 230] <- 1
pcaAll_FOG3$ecolength[pcaAll_FOG3$length > 230] <- 2

colors3 <- c("olivedrab","darkmagenta")

#############################
par(mfrow=c(1,3))
hist(pcaAll_FOG3$length,breaks=20,
  main="Lake FOG3",xlab="length (mm)")
abline(v=225, col="red",lty=2)

boxplot(EV1 ~ ecolength, data=pcaAll_FOG3,
        main="Lake FOG3",ylab="PC1")
points(pcaAll_FOG3$ecolength,pcaAll_FOG3$EV1,col=colors3[pcaAll_FOG3$ecolength])
t.test(EV1 ~ ecolength, data=pcaAll_FOG3)
#leveneTest(EV1 ~ ecolength, data=pcaAll[pcaAll$EV1 > 0 & pcaAll$EV2 > 0,])

##############################

par(mfrow=c(1,2))
plot(pcaAll_FOG3$EV1, pcaAll_FOG3$EV2, pch=19, cex=1, lwd=1, col=(scales::alpha(colors3[pcaAll_FOG3$length],1)),
     xlab=paste("PC1 (", round(pcSummary_FOG3$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_FOG3$importance[2,2]*100, 1), "%)", sep=""))
legend(x="bottomright",c("<230mm",">=230mm"),
  col=colors3,pch=19)

plot(pcaAll_FOG3$EV1,
     pcaAll_FOG3$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_FOG3$importance[2,1]*100, 1), "%)", sep=""),
     col=colors3[pcaAll_FOG3$ecolength],pch=19)

#### ggplot plots ####

pca.fog3 <- ggplot(data = pcaAll_FOG3, aes(x=EV1,y=EV2,color=as.factor(ecolength)))
pca.fog3.plot <- pca.fog3 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#ff3333","#ff8080"),labels=c("<230mm", ">=230mm")) +
  scale_fill_discrete(guide=FALSE)+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position=c(0.83,0.9),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",fill=NA),
        text=element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.border=element_rect(color="black",fill=NA))+
  xlab(paste("\\nPC1 (", round(pcSummary_FOG3$importance[2,1]*100, 1), "%)", sep=""))+
  ylab(paste("PC2 (", round(pcSummary_FOG3$importance[2,2]*100, 1), "%)\\n", sep=""))
  #ylim(0,0.00025)

length.fog3 <- ggplot(data = pcaAll_FOG3, aes(x=EV1,y=length,color=as.factor(ecolength)))
length.fog3.plot <- length.fog + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#ff3333","#ff8080"),labels=c("<230mm", ">=230mm"),guide=FALSE) +
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
  xlab(paste("\\nPC1 (", round(pcSummary_FOG3$importance[2,1]*100, 1), "%)", sep=""))+
  ylab("Length (mm)\\n")
  #ylim(0,0.00025)

plot_grid(pca.fog3.plot, length.fog3.plot, labels = NULL)

hist(pcaAll_FOG3$EV1)
abline(v = -0.0003, col="red")

group1 <- pcaAll_FOG3[pcaAll_FOG3$EV1 < -0.0003,]
group2 <- pcaAll_FOG3[pcaAll_FOG3$EV1 > -0.0003,]

table(group1$size_class)
table(group2$size_class)

plot(pcaAll_FOG3$EV1,pcaAll_FOG3$EV2,
  pch=19,col=pcaAll_FOG3$size_class)
## still doesn't correspond with these groups!

#### ggplot plots ####
library(cowplot)

pca.fog3 <- ggplot(data = subset(pcaAll_FOG3,!is.na(size_class)),
                  aes(x=EV1,y=EV2,color=as.factor(size_class)),na.rm=T)
pca.fog3.plot <- pca.fog3 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#ff3333","#ff8080","darkred"),labels=c("large", "medium", "small")) +
  scale_fill_discrete(guide=FALSE)+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position=c(0.85,0.87),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",fill=NA),
        text=element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.border=element_rect(color="black",fill=NA))+
  xlab(paste("\\nPC1 (", round(pcSummary_FOG3$importance[2,1]*100, 1), "%)", sep=""))+
  ylab(paste("PC2 (", round(pcSummary_FOG3$importance[2,2]*100, 1), "%)\\n", sep=""))
  #ylim(0,0.00025)

length.fog3 <- ggplot(data = subset(pcaAll_FOG3,!is.na(size_class)),
                     aes(x=EV1,y=length,color=as.factor(size_class)),na.rm=T)
length.fog3.plot <- length.fog3 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#ff3333","#ff8080","darkred"),labels=c("large", "medium", "small"),guide=FALSE) +
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
  xlab(paste("\\nPC1 (", round(pcSummary_FOG3$importance[2,1]*100, 1), "%)", sep=""))+
  ylab("Length (mm)\\n")
  #ylim(0,0.00025)

plot_grid(pca.fog3.plot, length.fog3.plot, labels = NULL)

#####################
## Lake FOG 5 only ##
#####################

FOG5 <- as.character(pcaAll$sample.id[pcaAll$lake == "FOG5"])

to.keep <- match(c(FOG5),colnames(char_alleles))
char_alleles_FOG5 <- char_alleles[,c(to.keep)]

missingness.FOG5 <- missingness[c(to.keep),]

char_pca_FOG5 <- do.pca(char_alleles_FOG5)
pcSummary_FOG5 <- summary(char_pca_FOG5)
scree_FOG5 <- plot(char_pca_FOG5, type="lines") 

pairedinfoFOG5<-merge(missingness.FOG5,fishinfo,
  by.x="New_ID",all.x=TRUE,sort=FALSE)

### combining fish info with PC results ##

pcaAll_FOG <- data.frame(sample.id = pairedinfoFOG5$New_ID,
                     year = factor(pairedinfoFOG5$YEAR),
                     lake = factor(pairedinfoFOG5$LAKE),
                     length = as.numeric(pairedinfoFOG5$LENGTH_mm),
                     weight = as.numeric(pairedinfoFOG5$WEIGHT_g),
                     lib = factor(pairedinfoFOG5$Library),
                     plates = factor(pairedinfoFOG5$Plate),
                     missing = as.numeric(pairedinfoFOG5$Missing),
                     heterozyg = pairedinfoFOG5$Heterozygosity,
                     date = as.factor(pairedinfoFOG5$DATE),
                     recap = as.factor(pairedinfoFOG5$RECAP),
                     sizeclass = as.factor(pairedinfoFOG5$LENGTH_CLASS),
                     EV1 = char_pca_FOG5$x[,1],    # the first eigenvector
                     EV2 = char_pca_FOG5$x[,2],    # the second eigenvector
                     EV3 = char_pca_FOG5$x[,3],    # the third eigenvector
                     #EV4 = char_pca_FOG5$x[,4],
                     #EV5 = char_pca_FOG5$x[,5],
                     stringsAsFactors = FALSE)


#############################
par(mfrow=c(1,2))
hist(pcaAll_FOG5$length,breaks=20,
  main="Lake FOG2",xlab="length (mm)")
abline(v=225, col="red",lty=2)

plot(pcaAll_FOG5$EV1,
     pcaAll_FOG5$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_FOG5$importance[2,1]*100, 1), "%)", sep=""),
     main="Lake FOG5",
     col=colors3[pcaAll_FOG5$sizeclass],pch=19)

##############################

par(mfrow=c(1,2))
plot(pcaAll_FOG5$EV1, pcaAll_FOG5$EV2, pch=19, cex=1, lwd=1, col=(scales::alpha(colors3[pcaAll_FOG5$sizeclass],1)),
     xlab=paste("PC1 (", round(pcSummary_FOG5$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_FOG5$importance[2,2]*100, 1), "%)", sep=""))
legend(x="bottomright",c("medium","large"),col=colors3,pch=19)

plot(pcaAll_FOG5$EV1,
     pcaAll_FOG5$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_FOG5$importance[2,1]*100, 1), "%)", sep=""),
     col=colors3[pcaAll_FOG5$sizeclass],pch=19)

#### ggplot plots ####
pca.fog5 <- ggplot(data = pcaAll_FOG5, 
                  aes(x=EV1,y=EV2,color=as.factor(sizeclass)),na.rm=T)
pca.fog5.plot <- pca.fog5 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#ff3333","#ff8080","darkred"),labels=c("large", "medium", "small")) +
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position=c(0.15,0.12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",fill=NA),
        text=element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.border=element_rect(color="black",fill=NA))+
  xlab(paste("\\nPC1 (", round(pcSummary_FOG5$importance[2,1]*100, 1), "%)", sep=""))+
  ylab(paste("PC2 (", round(pcSummary_FOG5$importance[2,2]*100, 1), "%)\\n", sep=""))
  #ylim(0,0.00025)

length.fog5 <- ggplot(data = pcaAll_FOG5, 
                     aes(x=EV1,y=length,color=as.factor(sizeclass)),
                     na.rm=T)
length.fog5.plot <- length.fog5 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#ff3333","#ff8080","darkred"),labels=c("large", "medium", "small"),guide=FALSE) +
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
  xlab(paste("\\nPC1 (", round(pcSummary_FOG5$importance[2,1]*100, 1), "%)", sep=""))+
  ylab("Length (mm)\\n")
  #ylim(0,0.00025)

plot_grid(pca.fog5.plot, length.fog5.plot, labels = NULL)

#####################
## Lake FOG 1 only ##
#####################

FOG1 <- as.character(pcaAll$sample.id[pcaAll$lake == "FOG1"])

to.keep <- match(c(FOG1),colnames(char_alleles))
char_alleles_FOG1 <- char_alleles[,c(to.keep)]

missingness.FOG1 <- missingness[c(to.keep),]

char_pca_FOG1 <- do.pca(char_alleles_FOG1)
pcSummary_FOG1 <- summary(char_pca_FOG1)
scree_FOG1 <- plot(char_pca_FOG1, type="lines") 

pairedinfoFOG1<-merge(missingness.FOG1,fishinfo,
  by.x="New_ID",all.x=TRUE,sort=FALSE)

### combining fish info with PC results ##

pcaAll_FOG <- data.frame(sample.id = pairedinfoFOG1$New_ID,
                     year = factor(pairedinfoFOG1$YEAR),
                     lake = factor(pairedinfoFOG1$LAKE),
                     length = as.numeric(pairedinfoFOG1$LENGTH_mm),
                     weight = as.numeric(pairedinfoFOG1$WEIGHT_g),
                     lib = factor(pairedinfoFOG1$Library),
                     plates = factor(pairedinfoFOG1$Plate),
                     missing = as.numeric(pairedinfoFOG1$Missing),
                     heterozyg = pairedinfoFOG1$Heterozygosity,
                     date = as.factor(pairedinfoFOG1$DATE),
                     recap = as.factor(pairedinfoFOG1$RECAP),
                     sizeclass = as.factor(pairedinfoFOG1$LENGTH_CLASS),
                     EV1 = char_pca_FOG1$x[,1],    # the first eigenvector
                     EV2 = char_pca_FOG1$x[,2],    # the second eigenvector
                     EV3 = char_pca_FOG1$x[,3],    # the third eigenvector
                     #EV4 = char_pca_FOG1$x[,4],
                     #EV5 = char_pca_FOG1$x[,5],
                     stringsAsFactors = FALSE)


#############################
par(mfrow=c(1,2))
hist(pcaAll_FOG1$length,breaks=20,
  main="Lake FOG2",xlab="length (mm)")
abline(v=225, col="red",lty=2)

plot(pcaAll_FOG1$EV1,
     pcaAll_FOG1$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_FOG1$importance[2,1]*100, 1), "%)", sep=""),
     main="Lake FOG1",
     col=colors3[pcaAll_FOG1$sizeclass],pch=19)

##############################

par(mfrow=c(1,2))
plot(pcaAll_FOG1$EV1, pcaAll_FOG1$EV2, pch=19, cex=1, lwd=1, col=(scales::alpha(colors3[pcaAll_FOG1$sizeclass],1)),
     xlab=paste("PC1 (", round(pcSummary_FOG1$importance[2,1]*100, 1), "%)", sep=""),
     ylab=paste("PC2 (", round(pcSummary_FOG1$importance[2,2]*100, 1), "%)", sep=""))
legend(x="bottomright",c("medium","large"),col=colors3,pch=19)

plot(pcaAll_FOG1$EV1,
     pcaAll_FOG1$length,
     ylab="Length (mm)",
     xlab=paste("PC1 (", round(pcSummary_FOG1$importance[2,1]*100, 1), "%)", sep=""),
     col=colors3[pcaAll_FOG1$sizeclass],pch=19)

#### ggplot plots ####
pca.fog1 <- ggplot(data = pcaAll_FOG1, 
                  aes(x=EV1,y=EV2,color=as.factor(sizeclass)),na.rm=T)
pca.fog1.plot <- pca.fog1 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#ff3333","#ff8080","darkred"),labels=c("large", "medium", "small")) +
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position=c(0.15,0.12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black",fill=NA),
        text=element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        panel.border=element_rect(color="black",fill=NA))+
  xlab(paste("\\nPC1 (", round(pcSummary_FOG1$importance[2,1]*100, 1), "%)", sep=""))+
  ylab(paste("PC2 (", round(pcSummary_FOG1$importance[2,2]*100, 1), "%)\\n", sep=""))
  #ylim(0,0.00025)

length.fog1 <- ggplot(data = pcaAll_FOG1, 
                     aes(x=EV1,y=length,color=as.factor(sizeclass)),
                     na.rm=T)
length.fog1.plot <- length.fog1 + geom_point(size=4) +
  #stat_ellipse(type = "norm", linetype = 2,geom="polygon",alpha=0.1) +
  scale_color_manual(values=c("#ff3333","#ff8080","darkred"),labels=c("large", "medium", "small"),guide=FALSE) +
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
  xlab(paste("\\nPC1 (", round(pcSummary_FOG1$importance[2,1]*100, 1), "%)", sep=""))+
  ylab("Length (mm)\\n")
  #ylim(0,0.00025)

plot_grid(pca.fog1.plot, length.fog1.plot, labels = NULL)