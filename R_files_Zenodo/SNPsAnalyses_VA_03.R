##############################################################################################
# COWPEA - ILLUMINA 50,000 SNP DATA ANALYSES - 20 GENOTYPES
############################################################################################
rm(list = ls())
#Libraries
install.packages("influence.ME")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
install.packages("rgl")
install.packages("snpReady")
library(rgl)
library(snpReady)
install.packages("SNPassoc")
library(SNPassoc)
install.packages("gdmp")
install.packages("adegenet")
library(adegenet)
library(glmmADMB)
library(data.table)
library(tidyr)
library(lme4)
library(MuMIn)
library(multcomp)
library(car)
library(lsmeans)
library(multcompView)
library(ggplot2)
library(Rmisc)
library(emmeans)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(RVAideMemoire)
library(mlbench)
library(MASS)
library(pROC)
library(Rfast)
library(dplyr)
library(matrixStats)
library(raster)
library(ggmap)
library(maptools)
library(maps)
library(ggplot2)
library(urbnmapr)
library(fitdistrplus)
library(logspline)
library(car)
library(Rmisc)
library(effects)
library(tidyverse)
library(urbnmapr)
library(lmerTest)
library(dfoptim)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(adephylo)
library(poppr)
library(ape)
library(RColorBrewer)
library(viridis)
library(influence.ME)

#****************************************************
#READING FILES
#****************************************************

#SNP data set contains  51, 128 SNPs and 20 cowpea lines
SNPS <- read.csv("Cowpea.51128SNPs.20lines.csv", header = TRUE, sep=",", strip.white = TRUE)
Genotype.Info <- read.csv("GenotypeInformation.csv", header = TRUE, sep=",", strip.white = TRUE)


#Reading Trait data:
Traits <- read.csv("CSEdata12Mar2020.csv", header = TRUE, sep=",", strip.white = TRUE)

#*********************************************************************
#BUILDING NEIHBOR-JOINING TREE WITH 51, 128 SNPs & all 20 lines
#*********************************************************************

# We are going to filter two columns from the dataset that have information about the chromosome
SNP.simple <- SNPS[,c(1,4:23)]
length(colnames(SNP.simple))
rownames(SNP.simple) <- SNP.simple$SNP  #The first column contains the SNPs, or what we want it to be the rownames
SNP.simple$SNP <- NULL
#head(SNP.simple) ---#checking order of things
#We need to transpose the matrix so it can match the format used by snpReady
SNP.simple.t <- as.data.frame(t(as.matrix(SNP.simple)))
#head(SNP.simple.t) #Now each SNP is a column
SNP.simple.t2 <- as.matrix(SNP.simple.t)
SNP.simple.t2[1:6,1:5]
SNP.simple.t2[ SNP.simple.t2 == "--" ] <- NA  #Need to remove the "--" and convert into NAs

#Recoding can also be done with the library snpReady, the good thing here is that it help us filter the SNPs prior to analyses
#Now we need to filter some of the SNPs that have more than 50% missing data
library(snpReady)
#if using cluster:
#colnames(SNP.simple.t2) <- SNP.simple.t2[1,]
#SNP.simple.t2 <- SNP.simple.t2[2:21,]
#class(SNP.simple.t2)

geno.ready <- raw.data(data = as.matrix(SNP.simple.t2), frame = "wide", base=TRUE, sweep.sample = 0.5, call.rate = 0.95, maf = 0.10, imput = FALSE)
geno.ready$report #This let us know which SNPs were removed
length(colnames(geno.ready$M.clean))  #34,762 left after filtering for missing data and minimum allele frequencies
x <- geno.ready$M.clean
x[10:5,34:5]  #checking matrix
snpfiltered.names <- colnames(geno.ready$M.clean)
length(snpfiltered.names)
rownames(geno.ready$M.clean)

#Estimating heterozygosity and inbreeding of each line:
pop.gen <- popgen(M = x)
he.fis <- pop.gen$whole$Genotypes
write.csv(he.fis, "he.fis.csv")

#extracting chromosome number and position for these filtered snps:
snpsselected <- dplyr::filter(SNPS, SNP %in% snpfiltered.names)
length(snpsselected$SNP) #checking that the number of snps matches

#extracting the line names to be able to assign the genepools
Genotype.Info.filtered <- dplyr::filter(Genotype.Info,Genotype %in% rownames(geno.ready$M.clean))
Genotype.Info.filtered$Genotype2 <- factor(Genotype.Info.filtered$Genotype, levels=rownames(geno.ready$M.clean))
Genotype.Info.filtered2 <- Genotype.Info.filtered[order(Genotype.Info.filtered$Genotype2),]
Genotype.Info.filtered2$Genepool
#Genotype.Info.filtered2$Genepool4


#Using a genlight object from adegenet to do genetic analyses
library(adegenet)
x1 <- new("genlight", x, parallel=FALSE) #if using cluster
x1@chromosome <- as.factor(snpsselected$Chr) #alternatives ways
x1@position <- snpsselected$Pos #alternative way

#assigning pop:
pop(x1) <- Genotype.Info.filtered2$Genepool
indNames(x1) 
glPlot(x1)


#Converting into fastSTRUCTURE format
library(dartR)
ploidy(x1) <- c("2")  #ploidy needs to be set for the transformation to work
gl2faststructure(x1,outfile = "Cowpea_faststr.str",outpath = ".",probar = FALSE,verbose = NULL)

#converting genlight to genind and then to STRUCTURE format
x1 #genlight object, needs to be transformed into a gind object
gind <- gl2gi(x1)  #converting genlight object to genind
genind2structure(gind, file="cowpeas.str", pops=FALSE)


#Making PCA
glPlot(x1, posi="topleft")
pca1 <- glPca(x1)
scatter(pca1, posi="bottomright")


pop80.pca <- glPca(x1, nf = 3)
barplot(100*pop80.pca$eig/sum(pop80.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance explained", line = 2)
title(xlab="Eigenvalues", line = 1)

#Calculating percent variation for each PC1
eigen.values <- data.frame(pop80.pca$eig)
total <- sum(eigen.values$pop80.pca.eig)
pc1.varpercent <- ((eigen.values$pop80.pca.eig[1])/total)*100
pc2.varpercent <- ((eigen.values$pop80.pca.eig[2])/total)*100
pc3.varpercent <- ((eigen.values$pop80.pca.eig[3])/total)*100


#Obtaining scores for each principal component to plot the PCA
pop80.pca.scores <- as.data.frame(pop80.pca$scores)
pop80.pca.scores$pop <- pop(x1)
pop80.pca.scores$ind <- indNames(x1)

library(ggplot2)
set.seed(9)
cols <- viridis(n = 3)
#cols <- viridis(n = nPop(x1))
p <- ggplot(pop80.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
#p <- p + geom_text(aes(label=ind), hjust=-0.5, vjust=0, size=3, position=position_jitter(width=0,height=3),check_overlap = F) 
p <- p + ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = ind))
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 

pca.cowpea <- p + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                                 axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                                 axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                                 axis.text.x = element_text(size=16, colour ="black"),
                                 axis.text.y = element_text(size=16, colour ="black"),
                                 panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                                 panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())

pca.cowpea <- p + theme_bw(base_size = 16)
pca.cowpea



#Making Neighbor-Joining Tree
library(ape)
library(phangorn)

x1 #genlight object, needs to be transformed into a gind object
gind <- gl2gi(x1)  #converting genlight object to genind

tree.upgma <- poppr::aboot(gind, tree = "upgma", distance = nei.dist, sample = 100, ploidy=2 , showtree = F, cutoff = 50, quiet = T)
tree.NJ <- poppr::aboot(gind, tree = "nj", distance = nei.dist, sample = 100, ploidy=2 , showtree = F, cutoff = 50, quiet = T)
#root_treeNJ <- midpoint(tree.NJ)
#root_PI632876 <- root(tree.NJ, outgroup = "PI632876", resolve.root = TRUE)

library("viridis")           
cols <- c("#440154FF","#21908CFF","#FDE725FF")

upgmatree.nei <- ggtree(tree.upgma, layout='rectangular') %<+% Genotype.Info.filtered2 +
  geom_tiplab(aes(label=label), align=FALSE, size=4) +
  #geom_text2(aes(subset=!isTip, label=Isolate.name),size=2, align=TRUE) +
  geom_tippoint(aes(color=Genepool), alpha=.9, size=2,show.legend = TRUE) + scale_colour_manual(values=cols) +
  #scale_color_manual(values=species.col, name="Host species") +
  geom_treescale(fontsize=3, offset=2, y=20) +
  #geom_label_repel(aes(label=bootstrap, fill=bootstrap)) #option when using the branchlenght tree
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label))), hjust=-.2, size=2.5) # & as.numeric(label) > 80))  #option if filtering bootstraps and using newick format

njtree.nei <- ggtree(tree.NJ, layout='rectangular') %<+% Genotype.Info.filtered2 +
  geom_tiplab(aes(label=label), align=FALSE, size=4, offset = 0.01) +
  #geom_text2(aes(subset=!isTip, label=Isolate.name),size=2, align=TRUE) +
  geom_tippoint(aes(color=Genepool), alpha=.9, size=3,show.legend = TRUE) + scale_colour_manual(values=cols) +
  #scale_color_manual(values=species.col, name="Host species") +
  geom_treescale(fontsize=3, offset=2, y=20) +
  #geom_label_repel(aes(label=bootstrap, fill=bootstrap)) #option when using the branchlenght tree
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label))), hjust=-.2, size=2.5) # & as.numeric(label) > 80))  #option if filtering bootstraps and using newick format


ggsave("upgmatree.nei.pdf", plot=upgmatree.nei, device="pdf", scale = 1, width = 30, height = 20, units = "cm", dpi = 600)
ggsave("njtree.nei.pdf", plot=njtree.nei, device="pdf", scale = 1, width = 25, height = 20, units = "cm", dpi = 600)

#ploting midpoint rooted tree:
njtree.nei_midroot <- ggtree(root_treeNJ, layout='rectangular') %<+% Genotype.Info.filtered2 +
  geom_tiplab(aes(label=label), align=FALSE, size=2) +
  #geom_text2(aes(subset=!isTip, label=Isolate.name),size=2, align=TRUE) +
  geom_tippoint(aes(color=Genepool), alpha=.9, size=1,show.legend = TRUE) + scale_colour_manual(values=cols) +
  #scale_color_manual(values=species.col, name="Host species") +
  geom_treescale(fontsize=3, offset=2, y=20) +
  #geom_label_repel(aes(label=bootstrap, fill=bootstrap)) #option when using the branchlenght tree
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label))), hjust=-.2, size=2.5) # & as.numeric(label) > 80))  #option if filtering bootstraps and using newick format

ggsave("njtree.nei_Midroot.pdf", plot=njtree.nei_midroot, device="pdf", scale = 1, width = 30, height = 20, units = "cm", dpi = 600)

njtree.nei_Rwild <- ggtree(root_PI632876, layout='rectangular') %<+% Genotype.Info.filtered2 +
  geom_tiplab(aes(label=label), align=FALSE, size=2) +
  #geom_text2(aes(subset=!isTip, label=Isolate.name),size=2, align=TRUE) +
  geom_tippoint(aes(color=Genepool), alpha=.9, size=1,show.legend = TRUE) + scale_colour_manual(values=cols) +
  #scale_color_manual(values=species.col, name="Host species") +
  geom_treescale(fontsize=3, offset=2, y=20) +
  #geom_label_repel(aes(label=bootstrap, fill=bootstrap)) #option when using the branchlenght tree
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label))), hjust=-.2, size=2.5) # & as.numeric(label) > 80))  #option if filtering bootstraps and using newick format

ggsave("njtree.nei_RootPI632876.pdf", plot=njtree.nei_Rwild, device="pdf", scale = 1, width = 30, height = 20, units = "cm", dpi = 600)


ggsave("upgmatree.nei.pdf", plot=upgmatree.nei, device="pdf", scale = 1, width = 30, height = 20, units = "cm", dpi = 600)
ggsave("njtree.nei.pdf", plot=njtree.nei, device="pdf", scale = 1, width = 30, height = 20, units = "cm", dpi = 600)




#################################################################
# ESTIMATING GENETIC DIVERSITY
#################################################################
x.mat <- as.matrix(x1) # gl.pops80 is a genlight object
x.mat[3:4,5:6] #checking the matrix
ind <- x1@ind.names  #Including the individual names
population <- as.character(x1@pop) #including the genepool names

x.gind <- df2genind(x.mat, sep = "/", ploidy = 2, ind.names=ind, pop=population) #making an genind object
library(hierfstat)

x.gind2 <- genind2hierfstat(x.gind)  #this is coverting the gind object into a hierfstat object
x.stats <- basic.stats(x.gind2) #Fst following Nei (1987), contains all stats
#x.stats <-he.fis.csv

Ho <- as.data.frame(x.stats$Ho)
mean(Ho$`Landrace genepool 2`)
se(Ho$`Landrace genepool 2`)

mean(Ho$`Landrace genepool 1`)
se(Ho$`Landrace genepool 1`)

mean(Ho$Wild)
se(Ho$Wild)

Hs <- as.data.frame(x.stats$Hs)
mean(Hs$`Landrace genepool 2`)
se(Hs$`Landrace genepool 2`)

mean(Hs$`Landrace genepool 1`)
se(Hs$`Landrace genepool 1`)

mean(Hs$Wild)
se(Hs$Wild)


#Gene diversity (Hs): -- which is is the chance that if you take two random alleles from your sample they are different
library(tidyr)
Hs <- as.data.frame(x.stats$Hs) 
Hs.short <- gather(Hs, Genepool, Hs)
Hs.short$Genepool <- as.character(Hs.short$Genepool)
loci <- rownames(Hs)
Hs.short$loci <- loci

Hs.plot <- ggplot(Hs.short, aes(x=Genepool, y=Hs, fill=Genepool)) + 
  geom_boxplot(outlier.colour="blue", outlier.shape=8,
               outlier.size=4, notch=TRUE) + #The notch displays a confidence interval around the median which is normally based on the median +/- 1.58*IQR/sqrt(n). Notches are used to compare groups; if the notches of two boxes do not overlap, this is a strong evidence that the medians differ
  scale_fill_manual(values = cols) +
  stat_summary(fun=mean, geom="point", shape=19, size=4) +
  labs(y= "Gene diversity (Hs)") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=16, colour ="black", angle=95),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())

#Observed heterozygosity (Ho):
Ho <- as.data.frame(x.stats$Ho)
Ho.short <- gather(Ho, Genepool, Ho)
Ho.short$Genepool <- as.character(Ho.short$Genepool)
loci <- rownames(Ho)
Ho.short$loci <- loci

cols <- c("#440154FF","#21908CFF","#FDE725FF")

Ho.plot <- ggplot(Ho.short, aes(x=Genepool, y=Ho, fill=Genepool)) + 
  geom_boxplot(outlier.colour="blue", outlier.shape=8,
               outlier.size=4, notch=TRUE) + #The notch displays a confidence interval around the median which is normally based on the median +/- 1.58*IQR/sqrt(n). Notches are used to compare groups; if the notches of two boxes do not overlap, this is a strong evidence that the medians differ
  scale_fill_manual(values = cols) +
  labs(y= "Observed Heterozygosity (Ho)") +
  stat_summary(fun=mean, geom="point", shape=19, size=4) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=16, colour ="black", angle=95),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())

library(Rmisc)
Ho_sum <- summarySE(Ho.short, measurevar="Ho", groupvars = "Genepool", na.rm=TRUE, conf.interval = 0.95)
                     
Ho.plot <- ggplot(Ho_sum, aes(x=Genepool, y=Ho, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=Ho - se, ymax=Ho+se), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual(values = cols) +
  labs(y= "Observed Heterozygosity (Ho)") +
  #stat_summary(fun=mean, geom="point", shape=19, size=4) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=16, colour ="black", angle=95),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())

library(tidyr)
Hs <- as.data.frame(x.stats$Hs)
Hs.short <- gather(Hs, Genepool,Hs)
Hs.short$Genepool <- as.character(Hs.short$Genepool)
loci <- rownames(Hs)
Hs.short$loci <- loci

Hs_sum <- summarySE(Hs.short, measurevar="Hs", groupvars = "Genepool", na.rm=TRUE, conf.interval = 0.95)

Hs.plot <- ggplot(Hs_sum, aes(x=Genepool, y=Hs, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=Hs - se, ymax=Hs+se), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual(values = cols) +
  labs(y= "Gene diversity (Hs)") +
  #stat_summary(fun=mean, geom="point", shape=19, size=4) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=16, colour ="black", angle=95),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())




diver.stats.plot <- ggarrange(Ho.plot,Hs.plot, nrow = 1, ncol=2, common.legend = TRUE, legend = "top")
ggsave("diver.stats.plot.pdf", plot=diver.stats.plot, scale = 1, width = 20, height = 15, units = "cm", dpi = 600)


#Testing whether there are significance differences in genetic diversity among genepools:
#* Trying linear mixed models:
library(lme4)
descdist(Hs.short$Hs, discrete=FALSE) #checking distribution of Hs -- looks like Beta distributed
hist(Hs.short$Hs) #not looking normal but let's try different transformations:
hist(subset(Hs.short, Genepool == "Landrace genepool 1")$Hs)
hist(subset(Hs.short, Genepool == "Landrace genepool 2")$Hs)
hist(subset(Hs.short, Genepool == "Wild")$Hs)


lmerhs <- lmer(Hs ~ Genepool + (1|loci), data=Hs.short)
#lmerhs <- lm(Hs ~ loci + Genepool, data=Hs.short)
hist(resid(lmerhs,type="deviance"))
qqnorm(resid(lmerhs,type="deviance"))
qqline(resid(lmerhs,type="deviance"))
plot(lmerhs) # not normal & weird heteroscesdacity

lmerho <- lmer(log(Ho + 1) ~ Genepool + (1|loci), data=Ho.short)
hist(resid(lmerho,type="deviance"))
qqnorm(resid(lmerho,type="deviance"))
qqline(resid(lmerho,type="deviance"))
plot(lmerho) #not normal & weird heteroscesdacity

#* Trying generalized linear mixed models with a Beta distribution:
#* Following AppendixS3 from Douman and Weedon 2019: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13234
install.packages("betareg")
install.packages("TMB")
install.packages("glmmTMB", type="source")
install.packages("lmtest")
install.packages("boot")
install.packages("emmeans")
install.packages("rstan", dependencies = TRUE)
install.packages("brms")


library(betareg)
library(plyr)
library(lmtest)
library(glmmTMB)
library(boot)
library(emmeans)
library(brms)

##* For HS:
#starting with a basic model without loci:
bm1 <- betareg(Hs ~ Genepool, data = Hs.short)
#rescaling the data to remove 0s and 1s
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

Hs.short$Hs.t <- transform01(Hs.short$Hs)

bm1 <- betareg(Hs.t ~ Genepool, data = Hs.short)
summary(bm1)

#now including loci as a random factor
glm.1 <- glmmTMB(Hs.t ~ Genepool, data = Hs.short, family = list(family = "beta", link = "logit"))
summary(glm.1)
glm.2 <- glmmTMB(Hs.t ~ Genepool + (1 | loci), data = Hs.short, family = list(family = "beta", link = "logit"))
summary(glm.2)
#The dispersion (ϕ) for the beta model is 1.88 and in this model specification is assumed to be the same for all treatments.
#to relax this assumption  we can incorporate a covariate model for the dispersion parameter. This allows difference variances in Hs within lineages
glm.3 <- update(glm.2, dispformula = ~Genepool)
summary(glm.3)


bbmle::AICtab(glm.2, glm.3) #model allowing non-constant precision for the Genepools is preferred 
lrtest(glm.3,glm.1)
lrtest(glm.3,glm.2)


#posthoc-tests:
lsmeans(glm.3, pairwise ~ Genepool)
f.lsm <- lsmeans(glm.3, "Genepool")


##* For Ho:
#starting with a basic model without loci:
bm1 <- betareg(Ho ~ Genepool, data = Ho.short)
#rescaling the data to remove 0s and 1s
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

Ho.short$Ho.t <- transform01(Ho.short$Ho)

bm1 <- betareg(Ho.t ~ Genepool, data = Ho.short)
summary(bm1)

#now including loci as a random factor
glm.1 <- glmmTMB(Ho.t ~ Genepool, data = Ho.short, family = list(family = "beta", link = "logit"))
summary(glm.1)
glm.2 <- glmmTMB(Ho.t ~ Genepool + (1 | loci), data = Ho.short, family = list(family = "beta", link = "logit"))
summary(glm.2)
#The dispersion (ϕ) for the beta model is 1.88 and in this model specification is assumed to be the same for all treatments.
#to relax this assumption  we can incorporate a covariate model for the dispersion parameter. This allows difference variances in Hs within lineages
glm.3.ho <- update(glm.2, dispformula = ~Genepool)
summary(glm.3.ho)


bbmle::AICtab(glm.2, glm.3.ho) #model allowing non-constant precision for the Genepools is preferred 
lrtest(glm.3.ho,glm.1)
lrtest(glm.3.ho,glm.2)

#posthoc-tests:
lsmeans(glm.3.ho, pairwise ~ Genepool)
f.lsm <- lsmeans(glm.3.ho, "Genepool")


#*************when including the two wild-lineages as separate -- four lineage analyses:
pop(x1) <- Genotype.Info.filtered2$Genepool4  #re-assigning genepools
x.mat <- as.matrix(x1) # gl.pops80 is a genlight object
x.mat[3:4,5:6] #checking the matrix
ind <- x1@ind.names  #Including the individual names
population <- as.character(x1@pop) #including the genepool names

x.gind <- df2genind(x.mat, sep = "/", ploidy = 2, ind.names=ind, pop=population) #making an genind object
library(hierfstat)
x.gind2 <- genind2hierfstat(x.gind)  #this is coverting the gind object into a hierfstat object
x.stats <- basic.stats(x.gind2) #Fst following Nei (1987), contains all stats


library(tidyr)
Hs <- as.data.frame(x.stats$Hs)
Hs.short <- gather(Hs, Genepool,Hs)
Hs.short$Genepool <- as.character(Hs.short$Genepool)
loci <- rownames(Hs)
Hs.short$loci <- loci


cols <- c("#440154FF","#21908CFF","#FDE725FF", "gold")

Hs.plot <- ggplot(Hs.short, aes(x=Genepool, y=Hs, fill=Genepool)) + 
  geom_boxplot(outlier.colour="blue", outlier.shape=8,
               outlier.size=4, notch=TRUE) + #The notch displays a confidence interval around the median which is normally based on the median +/- 1.58*IQR/sqrt(n). Notches are used to compare groups; if the notches of two boxes do not overlap, this is a strong evidence that the medians differ
  scale_fill_manual(values = cols) +
  stat_summary(fun=mean, geom="point", shape=19, size=4) +
  labs(y= "Gene diversity (Hs)") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=16, colour ="black", angle=95),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())

#Observed heterozygosity (Ho):
Ho <- as.data.frame(x.stats$Ho)
Ho.short <- gather(Ho, Genepool, Ho)
Ho.short$Genepool <- as.character(Ho.short$Genepool)
loci <- rownames(Ho)
Ho.short$loci <- loci


Ho.plot <- ggplot(Ho.short, aes(x=Genepool, y=Ho, fill=Genepool)) + 
  geom_boxplot(outlier.colour="blue", outlier.shape=8,
               outlier.size=4, notch=TRUE) + #The notch displays a confidence interval around the median which is normally based on the median +/- 1.58*IQR/sqrt(n). Notches are used to compare groups; if the notches of two boxes do not overlap, this is a strong evidence that the medians differ
  scale_fill_manual(values = cols) +
  labs(y= "Observed Heterozygosity (Ho)") +
  stat_summary(fun=mean, geom="point", shape=19, size=4) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=16, colour ="black", angle=95),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())

diver.stats.plot.4l <- ggarrange(Ho.plot,Hs.plot, nrow = 1, ncol=2, common.legend = TRUE, legend = "top")
ggsave("diver.stats.plot.4L.pdf", plot=diver.stats.plot.4l, scale = 1, width = 20, height = 15, units = "cm", dpi = 600)


Hs.short$Hs.t <- transform01(Hs.short$Hs)
bm1 <- betareg(Hs.t ~ Genepool, data = Hs.short)
summary(bm1)

#now including loci as a random factor
glm.1 <- glmmTMB(Hs.t ~ Genepool, data = Hs.short, family = list(family = "beta", link = "logit"))
summary(glm.1)
glm.2 <- glmmTMB(Hs.t ~ Genepool + (1 | loci), data = Hs.short, family = list(family = "beta", link = "logit"))
summary(glm.2)
#The dispersion (ϕ) for the beta model is 1.88 and in this model specification is assumed to be the same for all treatments.
#to relax this assumption  we can incorporate a covariate model for the dispersion parameter. This allows difference variances in Hs within lineages
glm.3 <- update(glm.2, dispformula = ~Genepool)
summary(glm.3)


bbmle::AICtab(glm.2, glm.3) #model allowing non-constant precision for the Genepools is preferred 
lrtest(glm.3,glm.1)
lrtest(glm.3,glm.2)


#posthoc-tests:
lsmeans(glm.3, pairwise ~ Genepool)
f.lsm <- lsmeans(glm.3, "Genepool")


##* For Ho:
#starting with a basic model without loci:
bm1 <- betareg(Ho ~ Genepool, data = Ho.short)
#rescaling the data to remove 0s and 1s
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

Ho.short$Ho.t <- transform01(Ho.short$Ho)

bm1 <- betareg(Ho.t ~ Genepool, data = Ho.short)
summary(bm1)

#now including loci as a random factor
glm.1 <- glmmTMB(Ho.t ~ Genepool, data = Ho.short, family = list(family = "beta", link = "logit"))
summary(glm.1)
glm.2 <- glmmTMB(Ho.t ~ Genepool + (1 | loci), data = Ho.short, family = list(family = "beta", link = "logit"))
summary(glm.2)
#The dispersion (ϕ) for the beta model is 1.88 and in this model specification is assumed to be the same for all treatments.
#to relax this assumption  we can incorporate a covariate model for the dispersion parameter. This allows difference variances in Hs within lineages
glm.3.ho <- update(glm.2, dispformula = ~Genepool)
summary(glm.3.ho)


bbmle::AICtab(glm.2, glm.3.ho) #model allowing non-constant precision for the Genepools is preferred 
lrtest(glm.3.ho,glm.1)
lrtest(glm.3.ho,glm.2)

#posthoc-tests:
lsmeans(glm.3.ho, pairwise ~ Genepool)
f.lsm <- lsmeans(glm.3.ho, "Genepool")


####################################################
# Testing genetic differentiation among lineages
####################################################

FST <- genet.dist(x.gind2, method = "WC84")
fst.pairwise <- pairwise.WCfst(x.gind2,diploid=TRUE) #Estimates pairwise FSTs according to Weir and Cockerham (1984
write.csv(fst.pairwise, "fst.pairwise.csv")
boot.fst <- boot.ppfst(dat=x.gind2,nboot=100,quant=c(0.025,0.975),diploid=TRUE)
write.csv(boot.fst, "boot.fst.csv")





#*********************************************************************
#ANALYZING TRAIT DATA
#*********************************************************************
#Because Two lines were discarded for the symbiosis trait analysis (PI632891 & PI632890) we will need to process a bit our SNP data set:

#Removing two of the cowpea lines that didn't formed nodules to construct G-matrix
removedlines <- c("PI632891","PI632890")
SNP_18lines <- SNPS[,!(names(SNPS) %in% removedlines)]
length(colnames(SNP_18lines))
head(SNP_18lines)
colnames(SNP_18lines)[4:15] <- c("TVu1280", "Tvu13305", "MuinanaLaw", "NamuesseD", "TVu14346", "Nhacoongo3",
                                 "Tvu9848", "TVu14971", "TVu15591", "TVu8834", "TVu9492", "TVu3804")

# We are going to filter two columns from the dataset that have information about the chromosome
SNP_18lines.simple <- SNP_18lines[,c(1,4:21)]
length(colnames(SNP_18lines.simple))
rownames(SNP_18lines.simple) <- SNP_18lines.simple$SNP  #The first column contains the SNPs, or what we want it to be the rownames
SNP_18lines.simple$SNP <- NULL
head(SNP_18lines.simple)
#We need to transpose the matrix so it can match the format used by snpReady
SNP_18lines.simple.t <- as.data.frame(t(as.matrix(SNP_18lines.simple)))
head(SNP_18lines.simple.t) #Now each SNP is a column
SNP.18lines.simple.t2 <- as.matrix(SNP_18lines.simple.t)
SNP.18lines.simple.t2[20,5]
SNP.18lines.simple.t2[ SNP.18lines.simple.t2 == "--" ] <- NA  #Need to remove the "--" and convert into NAs
#Recoding with the library snpReady, and filtering SNPs
library(snpReady)
geno.ready2 <- raw.data(data = as.matrix(SNP.18lines.simple.t2), frame = "wide", base=TRUE, sweep.sample = 0.5, call.rate = 0.95, maf = 0.10, imput = FALSE)
geno.ready2$report #This let us know which SNPs were removed
length(colnames(geno.ready2$M.clean))  #34,762 left after filtering for missing data and minimum allele frequencies
x2 <- geno.ready2$M.clean

#Obtaining G-matrix using the kindship relationship described in VanRanden 2008
#two matrices additive and dominance are generated
G <- G.matrix(M = x2, method = "VanRaden", format = "wide")
Ga <- G$Ga  #Additive genetic matrix
  
#checking if names match:
 Lines.names.traits <- unique(Traits$Line) #list of the lines in the Trait data
 Lines.names.Gu <- rownames(Ga) #List of the lines in the matrix data
  ltraits <- unlist(Lines.names.traits) #need to remove format
  lgenetic <- unlist(Lines.names.Gu) #need to remove format
  difg <- setdiff(ltraits,lgenetic) #what is different is ltraits? -- some lines had different names so I need it to change that in the original trait data

#Obtaining additive, dominance and epistatic relationship matrix:
library(sommer)
  A <- A.mat(x2)
  D <- D.mat (x2)
  E <- E.mat(x2)
  

#******************************************************************************
  #DEFINING FUNCTIONS TO BE USED FOR EACH TRAIT ANALYSIS
#****************************************************************************** 
  
  #*Boxplot to visualize outliers
  library(ggplot2)
  bp <- function(mydata, x.var1, y.var, factor2){
    dummy.data <- mydata[,c(x.var1, y.var, factor2)]
    p <- colnames(dummy.data)
    x <- dummy.data[,c(x.var1)]
    y <- dummy.data [, c(y.var)]
    c <- dummy.data [, c(factor2)]
    outplot <- ggplot(data = dummy.data, aes(x= x, y= y, fill= c)) + 
      geom_boxplot(aes(fill=c), outlier.colour="black") +
      xlab(p[1]) + ylab(p[2])
      print(outplot)
      }
        #usage of this function:
        #bp(Traits.nocontrol, "Treatment", "NumberOfNodules", "Line") #data, factorx, factoryy, fill factor
        
  #*Boxplot to remove outliers: creates a list of outliers to be removed from main data
  outliers.fun <- function (mydata, y.var, x.var, c.var,d.var) {
    dummy.data <- mydata[,c(x.var, y.var, c.var, d.var)]
    p <- colnames(dummy.data)
    x <- dummy.data[,c(x.var)]
    y <- dummy.data [, c(y.var)]
    c <- dummy.data [, c(c.var)]
    d <- dummy.data[,c(d.var)]
    
    b =boxplot(y ~ x*c*d, data=dummy.data)$out  #This tells the specific values that are outlier
  o <- mydata[y %in% b,] #This let me see which values are considered outliers
  return(o)
  }
  #function usage:
  #outliers.fun(Traits.nocontrol, "NumberOfNodules", "Treatment", "Line", "Genepool")
 
  
  
  
  #*We will be using different models of variance-covariance structure using the R package sommer:
  
  
  
  
  
  
#******************************************************************************
#USING CLASSICAL QUANT DESIGNS BASED ON THE INBREED DESIGN THAT WE HAVE
#******************************************************************************

#we need to subset the control treatment from our data since we don't want to compare this treatment:
Traits.nocontrol <- subset(Traits, Treatment != "C")
Traits.nocontrol <- subset(Traits.nocontrol, NumberOfNodules > 0)
Traits.nocontrol$Treatment <- factor(Traits.nocontrol$Treatment)
Traits.nocontrol$Line <- factor(Traits.nocontrol$Line) 

Traits.nocontrol$Block <- as.character(Traits.nocontrol$Block) #Making sure these are read as character
Traits.nocontrol$HarvestWeek <- as.character(Traits.nocontrol$HarvestWeek)

 
#----------------------------------------------------------------------
# -------------   TRAIT 1: Host growth  ----------------
#----------------------------------------------------------------------

#let's check how the data looks like
qqPlot(Traits.nocontrol$Host.Growth.Difference)
#Let's have data log-transformed in case of heteroscesdacity in models & not normality
Traits.nocontrol$log.Host.Growth <- log((Traits.nocontrol$Host.Growth.Difference) + 1 - min(Traits.nocontrol$Host.Growth.Difference))
qqPlot(Traits.nocontrol$log.Host.Growth)

#looking the distribution of the data:
bp(Traits.nocontrol, "Treatment", "Host.Growth.Difference", "Genepool")
bp(Traits.nocontrol, "Treatment", "Host.Growth.Difference", "Line")


#Let's find out which ones are the outliers
o <- outliers.fun(Traits.nocontrol, "Host.Growth.Difference", "Treatment", "Line", "Genepool")

#Let's see what happens when we remove the most extrem outlier: Plant with 480 nodules
length(o$log.Host.Growth)
length(Traits.nocontrol$log.Host.Growth)
Traits.nocontrol.noout <- Traits.nocontrol[!Traits.nocontrol$log.Host.Growth %in% o$log.Host.Growth,]
length(Traits.nocontrol.noout$log.Host.Growth)

#NOW LET'S MODEL IF THAT GENETIC VARIATION DIFFERS AMONG TREATMENTS
# Using sommer package - which uses the GBLUPs
#####################################################   

#Model A: Model assuming homogeneous variances across strain treatments (e.g. across environments). In other words, genotypes have the same variance component across treatments (do not vary across treatments) 
  #*Compound simmetry model (CS): this is our null model - equal among-line variance (VA) in all treatments
  #* Here are using all lines independent of their genepool
install.packages("sommer")
library(sommer)
dev.off()
head(Traits.nocontrol)
#let's try to remove one of the outliers:

HG.CSmodel <- mmer(log.Host.Growth ~ Treatment + Genepool + DaysSinceInoculation,tolparinv = 1e-03,
                              random = ~ Line + Block + Treatment:Line,
                              rcov = ~ units,
                              data = Traits.nocontrol)


Traits.nocontrol.noout
summary(HG.CSmodel)
plot(HG.CSmodel) #testing residuals to asses the model
plot(HG.CSmodel,id=0.05,idLabels=~.obs)
HG.CSmodel$residuals
HG.CSmodel$AIC

    #Representation of this model in lme4:
    HG.CSmodel.lmer <- lmer(log.Host.Growth ~ DaysSinceInoculation + Genepool + Treatment + (1|Block) + (1|Line) + (1|Treatment:Line), data = Traits.nocontrol)
    Anova(HG.CSmodel.lmer)
    plot(HG.CSmodel.lmer) #looking at heteroscesdacity
    qqnorm(residuals(HG.CSmodel.lmer))
    library(influence.ME)
    infl <- influence( HG.CSmodel.lmer, obs = TRUE)
    #Calculate Cook's distance:
    cooks.distance(infl)
    #Plot Cook's distance:
    plot(infl, which = "cook")
    
    #finding random effects (e.g. blups)
    ranef(HG.CSmodel.lmer)
    plot(allEffects(HG.CSmodel.lmer, partial.residuals = TRUE))
    dotplot(ranef(HG.CSmodel.lmer, condVar = TRUE)) #Visualizing the random effects
    qqmath(ranef(HG.CSmodel.lmer, condVar = TRUE)) #This also does the the same but it gets the Q-Qplots
    

#Model B:  assuming a heterogeneous variance across strain treatments (e.g. across environments). 
#In other words this model allows the among-line variance (genetic variance) to differ among treatments
#This is done using ds function (it requires CS+DIAG)
                      #Assumes that each environment has a different GxE variance

HG.HTmodel <-mmer(log.Host.Growth~Treatment + Genepool + DaysSinceInoculation,tolparinv = 1e-03,
                             random=~Line+ Block + vs(ds(Treatment),Line),
                             rcov=~ vs(ds(Treatment),units),
                             data=Traits.nocontrol)
summary(HG.HTmodel)
plot(HG.HTmodel)
HG.HTmodel$AIC

            #In lmer4 this would be equivalent to:
            HG.HSmodel.lmer <- lmer(log.Host.Growth ~ DaysSinceInoculation + HarvestWeek + Treatment + (Treatment|Line), data = Traits.nocontrol.noout)
            plot(HG.HSmodel.lmer)
            ranef(HG.HSmodel.lmer)
            plot(allEffects(HG.HSmodel.lmer, partial.residuals = TRUE))
            dotplot(ranef(HG.HSmodel.lmer, condVar = TRUE)) #Visualizing the random effects
            qqmath(ranef(HG.HSmodel.lmer, condVar = TRUE)) #This also does the the same but it gets the Q-Qplots


#Model B.1:Using our relationship matrix
                        HG.HTmodel.pedigree <-mmer(log.Host.Growth~Treatment + Genepool + DaysSinceInoculation + HarvestWeek,tolparinv = 1e-03,
                                                   random=~ Line+ vs(ds(Treatment),Line, Gu=A),    #G <- A.mat(x2) could also be used instead
                                                   rcov=~ vs(ds(Treatment),units),
                                                   data=Traits.nocontrol)
                        
                    
                        plot(HG.HTmodel.pedigree)
                        HG.HTmodel.pedigree$U$Line
                        HG.HTmodel.pedigree$AIC
                        
                       

#Model C : assuming an unstructure variance, meaning that among the levels of certain factor (i.e. Environments) there’s a covariance struture of a second random effect(i.e. Genotypes)
#Estimates different GxE at each treatment
          HG.USmodel <-mmer(log.Host.Growth~Treatment + Genepool + DaysSinceInoculation + HarvestWeek,tolparinv = 1e-03,
                  random=~ vs(us(Treatment),Line),
                  rcov=~ vs(us(Treatment),units),
                  data=Traits.nocontrol)

plot(HG.USmodel)
summary(HG.USmodel)
HG.USmodel$AIC


                  #Model C.1. using pedigree information
                  HG.USmodel.pedigree <-mmer(log.Host.Growth ~ Treatment + Genepool + HarvestWeek,tolparinv = 1e-0,
                                    random=~ vs(us(Treatment),Line, Gu=A),
                                    rcov=~ vs(us(Treatment),units),
                                    data=Traits.nocontrol) #This model says that sistem is singular
                  plot(HG.USmodel.pedigree) #singular model
                  

                  

#######  Which model one is the best model?
AIC <- c(HG.CSmodel$AIC,HG.HTmodel$AIC, HG.HTmodel.pedigree$AIC,HG.USmodel$AIC, "NA") 
Model <- c("Constrained", "Heterogeneous", "Heterogeneous-Amatrix", "Unstructured", "Unstructured-Amatrix")
VarModels <- data.frame(cbind(Model,AIC))
VarModels[order(VarModels$AIC),] #model with the lowest AIC is the best approximation   
               

#####  Performing like-lihood ratio tests among models:

#let's estimate the number of covariate parameters that were estimated:
summary(HG.CSmodel) 
summary(HG.HTmodel)
summary(HG.HTmodel.pedigree)
summary(HG.USmodel)
HG.USmodel$PevU

df.a = 3
df.b = 4
df.b1 = 4
df.c = 6

#* Model A vs B
chi.square= (-2*(-98.95224) + 2*(-97.30872)) #difference between the log-likelihood of both models
DF= abs(df.a - df.b)
pchisq(chi.square, df=DF, lower.tail = F) #Based on this there is among-genotype variance within treatments

#* Model A vs B.1
chi.square= (-2*(-98.95224) + 2*(-96.96055)) #difference between the log-likelihood of both models
DF = abs(df.a - df.b1)
pschi.A.B.1 <- pchisq(chi.square, df=DF, lower.tail = F) 


#* Model A vs C
chi.square= (-2*(-98.95224) + 2*(-94.10389)) #difference between the log-likelihood of both models
DF=  abs(df.a - df.c) #number of var comp estimated
pchisq(chi.square, df=DF, lower.tail = F) 

#* Model B vs C
chi.square= (-2*(-97.30872) + 2*(-94.10389)) #difference between the log-likelihood of both models
DF= abs(df.b - df.c)
pchisq(chi.square, df=DF, lower.tail = F) 


#Getting BLUPS for Host growth based on the unstructured variance model without the kinship matrix
BLUP.A <- as.data.frame(HG.USmodel$U$`A:Line`) #obtaining BLUPs
A.SE <- as.data.frame(HG.USmodel$PevU$`AL:Line`) #obtaining standard errors of BLUPs

#Plotting reaction norms with just the raw BLUPS
BLUP.US.A <- as.data.frame(HG.USmodel$U$`A:Line`)   
BLUP.US.AL <- as.data.frame(HG.USmodel$U$`AL:Line`)
BLUP.US.L <- as.data.frame(HG.USmodel$U$`L:Line`)

BLUPS.Hostgrowth.US <- BLUP.US.A
BLUPS.Hostgrowth.US$A <- BLUP.US.A$log.Host.Growth
BLUPS.Hostgrowth.US$AL <- BLUP.US.AL$log.Host.Growth
BLUPS.Hostgrowth.US$L <- BLUP.US.L$log.Host.Growth
BLUPS.Hostgrowth.US$Line <- rownames(BLUP.US.A) 

BLUPS.Hostgrowth.US <-BLUPS.Hostgrowth.US[, c(2:5)] #subsetting for only the raw BLUPs
BLUPS.Hostgrowth.US.short <- gather(BLUPS.Hostgrowth.US, Treatment, BLUPs, -Line)
BLUPS.Hostgrowth.US.short$Genepool <- c("2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2",
                                             "2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2",
                                             "2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2")


BLUPS.Hostgrowth.US.short$bt.BLUP <- exp(BLUPS.Hostgrowth.US.short$BLUPs)

library("viridis")           
cols <- viridis(n = 3) #selecting colors for each genepool
HostGrowht.USVar.BLUPs.plot <- ggplot(BLUPS.Hostgrowth.US.short, aes(x=Treatment, y=BLUPs, group=Line, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  #geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=0.05, colour="black") +
  geom_point(shape=21, colour = "black", size=3) +
  #geom_text(data=BLUPS.Hostgrowth.US.short %>% filter(Treatment == "A"), 
          #  aes(label=Line), hjust=1.5, vjust=0, size=3, check_overlap = F) +
  ggrepel::geom_text_repel(data=BLUPS.Hostgrowth.US.short %>% filter(Treatment == "A"), aes(x = Treatment, y = BLUPs, label = Line)) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("Unstructured Variance Model-GBLUPS") +
  labs(y="log Host growth difference (BLUPs)", x= "Strain Treatment") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 


HostGrowht.USVar.btBLUPs.plot <- ggplot(BLUPS.Hostgrowth.US.short, aes(x=Treatment, y=bt.BLUP, group=Line, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  geom_point(shape=21, colour = "black", size=3) +
  #geom_text(data=BLUPS.Hostgrowth.US.short %>% filter(Treatment == "A"), 
            #aes(label=Line), hjust=1.5, vjust=0, size=3, check_overlap = F) +
  #scale_fill_viridis(discrete = TRUE) +
  #ggtitle("Unstructured Variance Model-BLUPS") +
  labs(y="Host growth difference (mg, back-transformed BLUPs)", x= "Strain Treatment") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 


#Let's compare with raw host growth:

#let's compare the predicted trait values with the observed values:
HostGrowth.raw.means = summarySE(data=Traits.nocontrol.noout, measurevar= "log.Host.Growth", groupvars = c("Genepool", "Treatment","Line"), na.rm=TRUE)
HostGrowth.raw.means$bt.HG <- exp(HostGrowth.raw.means$log.Host.Growth)

HostGrowth.RN.raw <- ggplot(HostGrowth.raw.means, aes(x=Treatment, y=log.Host.Growth, group=Line, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  #geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=0.05, colour="black") +
  geom_point(shape=21, colour = "black", size=3) +
  #geom_text(aes(label=Line),hjust=0, vjust=0) +
  geom_text(data=HostGrowth.raw.means %>% filter(Treatment== "A"), 
            aes(label=Line), hjust=1.2, vjust=0, size=3, check_overlap = F) +
  geom_text(data=HostGrowth.raw.means %>% filter(Treatment== "L"), 
            aes(label=Line), hjust=-0.5, vjust=0, size=3, check_overlap = F) +
  scale_fill_viridis(discrete = TRUE) +
  labs(y="log Host growth (mg)", x= "Strain Treatment") +
  ggtitle("Raw Mean Host growth") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 



###########################Calculating heritability:
#sommer produces the output for (V1) the additive genetic variance, (V2) the variance between A and AL, (V10) the residuals (termed "units").
#heritability estimates based on the Unstructured variance model:

sum <- summary(HG.USmodel)
varcomp <- data.frame(sum$varcomp[1])
length(varcomp$VarComp)  #we have 12 variance components, let's add numbers to then set up the variances used for heritabilities
varcomp$number <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#now using the pin() function we can access these different variances:
VP.A = pin(HG.USmodel, Vp.HG~(V1+V7)) #phenotypic variance within treatment A equals the Vg (among lines) + residual (among-replicates)
VP.AL = pin(HG.USmodel, Vp.HG~(V3+V9)) #phenotypic variance within treatment AL
VP.L = pin(HG.USmodel, Vp.HG~(V6+V12)) #phenotypic variance within treatment L

Hb.A = pin(HG.USmodel, Vp.HG~V1/(V1+V7)) #Heritability
Hb.AL = pin(HG.USmodel, Vp.HG~ V3/(V3+V9)) #Heritability
Hb.L = pin(HG.USmodel, Vp.HG~V6/(V6+V12)) #Heritability

#overall VP (sum of all genetic variances)
VP = pin(HG.USmodel, Vp.HG~(V1+V7+V3+V9+V6+V12)) #phenotypic variance conditional on fixed effects in the model
VG = pin(HG.USmodel, Vg.HG~(V1+V3+V6))
Hb= VG/VP


#*******************************************************************
# --- Complex models to estimate variance components for genepools
#******************************************************************
#Models are intercept only models

#Model A: Constrained variance model:
#*assumes that the random genotype and genotype by treatment interaction effects are constant. It involved only two genetic parameters; a variance and a correlation
#*residuals and blocks nested within sites have scaled identity variance structures
ModelA <- mmer(log.Host.Growth ~ Treatment + Genepool,
                   random = ~ Line + Block + Treatment:Line + Genepool:Line,
                   rcov = ~ units,
                   data = Traits.nocontrol)

summary(ModelA)
plot(ModelA)


#Model B: heterogeneous variance model - residual error and random effect(line) allowed to vary among treatments
#allows among-genotype variances to differ among treatments
ModelB <- mmer(log.Host.Growth ~ Treatment + Genepool + DaysSinceInoculation,
                   random = ~ Line + Block + vs(ds(Treatment), Line),
                   rcov = ~ vs(ds(Treatment),units),
                   data = Traits.nocontrol)

summary(ModelB)
plot(ModelB)

#Model C: heterogeneous variance model - testing if the among-genotype variances differed among treatments, while keeping
#*allows among-genotype variances to differ among genepools
ModelC <- mmer(log.Host.Growth ~ Treatment + DaysSinceInoculation,
                   random = ~ Line + Genepool + vs(ds(Genepool), Line),
                   rcov = ~ vs(ds(Genepool),units),
                   data = Traits.nocontrol)

summary(ModelC)
plot(ModelC)

#Model D: heterogeneous variance model - testing if the among-genotype variances differed among treatments, while keeping
#*
ModelD <- mmer(log.Host.Growth ~ Treatment + Genepool + DaysSinceInoculation,
                   random = ~ Line + vs(ds(Treatment),ds(Genepool), Line),
                   rcov = ~ vs(ds(Treatment),ds(Genepool),units),
                   data = Traits.nocontrol)


summary(ModelD)   
plot(ModelD)



anova(lm(Block ~ DaysSinceInoculation, Traits.nocontrol))


#Model D.1: heterogeneous variance model including relationship matrix
HS.HSmodel.pedigree <- mmer(log.Host.Growth ~ Treatment + Genepool + DaysSinceInoculation,
                            random = ~ Line + Block + vs(ds(Treatment), ds(Genepool), Line, Gu=A),
                            rcov = ~ vs(ds(Treatment), ds(Genepool),units),
                            data = Traits.nocontrol)

summary(HS.HSmodel.pedigree)
plot(HS.HSmodel.pedigree)


#----Performing like-lihood ratio tests among models:

#let's estimate the number of covariate parameters that were estimated:
summary(HS.CSmodel) 
summary(HS.HSmodel)
summary(HS.HSmodel.pedigree)
summary(HS.USmodel)
summary(HS.USmodel.pedigree)
HG.USmodel$PevU

df.a = 5
df.D = 20 # four VG components are now estimated - also 3 residual variance components
df.D1 = 20


#* Model A vs D
chi.square= (-2*(-106.51) + 2*(79.69)) #difference between the log-likelihood of both models
DF= abs(df.a - df.b)
pchisq(chi.square, df=DF, lower.tail = F) #Based on this there is among-genotype variance within treatments
#*conclusion: among-genotype variances differed among treatments
#* Model A vs D.1
chi.square= (-2*(-106.51) + 2*(84.6886)) #difference between the log-likelihood of both models
DF = abs(df.a - df.b1)
pchisq(chi.square, df=DF, lower.tail = F) 
#*conclusion: among-genotype variances differed among treatments, while controling for line relationships


###########################Calculating heritability:
#sommer produces the output for (V1) the additive genetic variance, (V2) the variance between A and AL, (V10) the residuals (termed "units").
#heritability estimates based on the Heterogeneous variance model, that includes the relationship matrix:

sum <- summary(NN.HSmodel.pedigree)
varcomp <- data.frame(sum$varcomp[1])
length(varcomp$VarComp)  #we have 12 variance components, let's add numbers to then set up the variances used for heritabilities
varcomp$number <- c(1:20)
#now using the pin() function we can access these different variances:
VP.A.genepool1 = pin(NN.HSmodel.pedigree, Vp.NN~(V3+V12)) #phenotypic variance within treatment A equals the Vg (among lines) + residual (among-replicates)
VP.A.genepool2 = pin(NN.HSmodel.pedigree, Vp.NN~(V4+V13)) #phenotypic variance within treatment A equals the Vg (among lines) + residual (among-replicates)
VP.A.wild = pin(NN.HSmodel.pedigree, Vp.NN~(V5+V14)) #phenotypic variance within treatment A equals the Vg (among lines) + residual (among-replicates)

VP.AL.genepool1 = pin(NN.HSmodel.pedigree, Vp.NN~(V6+V15)) #phenotypic variance within treatment AL
VP.AL.genepool2 = pin(NN.HSmodel.pedigree, Vp.NN~(V7+V16)) #phenotypic variance within treatment AL
VP.AL.wild = pin(NN.HSmodel.pedigree, Vp.NN~(V8+V17)) #phenotypic variance within treatment AL

VP.L.genepool1 = pin(NN.HSmodel.pedigree, Vp.NN~(V9+V18)) #phenotypic variance within treatment L
VP.L.genepool2 = pin(NN.HSmodel.pedigree, Vp.NN~(V10+V19)) #phenotypic variance within treatment L
VP.L.wild = pin(NN.HSmodel.pedigree, Vp.NN~(V11+V20)) #phenotypic variance within treatment L


VG.A.genepool1 = pin(NN.HSmodel.pedigree, Vg.NN~(V3))
VG.A.genepool2 = pin(NN.HSmodel.pedigree, Vg.NN~(V4))
VG.A.wild = pin(NN.HSmodel.pedigree, Vg.NN~(V5))

VG.AL.genepool1 = pin(NN.USmodel.pedigree, Vg.NN~(V6))
VG.AL.genepool2 = pin(NN.USmodel.pedigree, Vg.NN~(V7))
VG.AL.wild = pin(NN.USmodel.pedigree, Vg.NN~(V8))

VG.L.genepool1 = pin(NN.USmodel.pedigree, Vg.NN~(V9))
VG.L.genepool2 = pin(NN.USmodel.pedigree, Vg.NN~(V10))
VG.L.wild = pin(NN.USmodel.pedigree, Vg.NN~(V11))

Hb.A.genepool1 = pin(NN.HSmodel.pedigree, Hb.NN~V3/(V3+V12)) #Heritability
Hb.A.genepool2 = pin(NN.HSmodel.pedigree, Hb.NN~V4/(V4+V13)) #Heritability
Hb.A.genepool3 = pin(NN.HSmodel.pedigree, Hb.NN~V5/(V5+V14)) #Heritability

Hb.AL.genepool1 = pin(NN.HSmodel.pedigree, Hb.NN~V6/(V6+V15)) #Heritability
Hb.AL.genepool2 = pin(NN.HSmodel.pedigree, Hb.NN~V7/(V7+V16)) #Heritability
Hb.AL.genepool3 = pin(NN.HSmodel.pedigree, Hb.NN~V8/(V8+V17)) #Heritability

Hb.L.genepool1 = pin(NN.HSmodel.pedigree, Hb.NN~V9/(V9+V18)) #Heritability
Hb.L.genepool2 = pin(NN.HSmodel.pedigree, Hb.NN~V10/(V10+V19)) #Heritability
Hb.L.genepool3 = pin(NN.HSmodel.pedigree, Hb.NN~V11/(V11+V20)) #Heritability


#overall VP (sum of all genetic variances)
VP = pin(NN.USmodel.pedigree, Vp.NN~(V1+V7+V3+V9+V6+V12)) #phenotypic variance conditional on fixed effects in the model
VG = pin(NN.USmodel.pedigree, Vg.NN~(V1+V3+V6))
Hb= VG/VP


#Obtaining BLUPS from heterogeneous variance model with A-matrix
c <- data.frame(Traits.nocontrol$Line)
c$Genepool <- Traits.nocontrol$Genepool
c <- unique.data.frame(c)
c.genepool1 <- subset(c, Genepool == "One")
c.genepool2 <- subset(c, Genepool == "Two")
c.wild <- subset(c, Genepool == "Wild")

blups.A.one <- data.frame(NN.HSmodel.pedigree$U$`A:One:Line`)
blups.A.one$Lines <- rownames(blups.A.one)
blups.A.one2 <-blups.A.one[blups.A.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.A.one2$Genepool <- c("One")

blups.A.Two <- data.frame(NN.HSmodel.pedigree$U$`A:Two:Line`)
blups.A.Two$Lines <- rownames(blups.A.Two)
blups.A.Two2 <-blups.A.Two[blups.A.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.A.Two2$Genepool <- c("Two")

blups.A.Wild <- data.frame(NN.HSmodel.pedigree$U$`A:Wild:Line`)
blups.A.Wild$Lines <- rownames(blups.A.Wild)
blups.A.Wild2 <-blups.A.Wild[blups.A.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.A.Wild2$Genepool <- c("Wild")

blups.A <- rbind(blups.A.one2,blups.A.Two2,blups.A.Wild2)
blups.A$Treatment <- c("A")

blups.AL.one <- data.frame(NN.HSmodel.pedigree$U$`AL:One:Line`)
blups.AL.one$Lines <- rownames(blups.AL.one)
blups.AL.one2 <-blups.AL.one[blups.AL.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.AL.one2$Genepool <- c("One")

blups.AL.Two <- data.frame(NN.HSmodel.pedigree$U$`AL:Two:Line`)
blups.AL.Two$Lines <- rownames(blups.AL.Two)
blups.AL.Two2 <-blups.AL.Two[blups.AL.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.AL.Two2$Genepool <- c("Two")

blups.AL.Wild <- data.frame(NN.HSmodel.pedigree$U$`AL:Wild:Line`)
blups.AL.Wild$Lines <- rownames(blups.AL.Wild)
blups.AL.Wild2 <-blups.AL.Wild[blups.AL.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.AL.Wild2$Genepool <- c("Wild")

blups.AL <- rbind(blups.AL.one2,blups.AL.Two2,blups.AL.Wild2)
blups.AL$Treatment <- c("AL")

blups.L.one <- data.frame(NN.HSmodel.pedigree$U$`L:One:Line`)
blups.L.one$Lines <- rownames(blups.L.one)
blups.L.one2 <-blups.L.one[blups.L.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.L.one2$Genepool <- c("One")

blups.L.Two <- data.frame(NN.HSmodel.pedigree$U$`L:Two:Line`)
blups.L.Two$Lines <- rownames(blups.L.Two)
blups.L.Two2 <-blups.L.Two[blups.L.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.L.Two2$Genepool <- c("Two")

blups.L.Wild <- data.frame(NN.HSmodel.pedigree$U$`L:Wild:Line`)
blups.L.Wild$Lines <- rownames(blups.L.Wild)
blups.L.Wild2 <-blups.L.Wild[blups.L.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.L.Wild2$Genepool <- c("Wild")

blups.L <- rbind(blups.L.one2,blups.L.Two2,blups.L.Wild2)
blups.L$Treatment <- c("L")

blups.NN <- rbind(blups.L,blups.AL, blups.A)
blups.NN$btBLUPs <- (blups.NN$sqr.NN)^2

#Plotting reaction norms with just the raw BLUPS
NodNum.USVar.GBLUPs.plot <- ggplot(blups.NN, aes(x=Treatment, y=btBLUPs, group=Lines, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  geom_point(shape=21, colour = "black", size=3) +
  geom_text(data=blups.NN %>% dplyr::filter(Treatment == "L"), 
            aes(label=Lines), hjust=-0.5, vjust=0, size=3, check_overlap = F) +
  scale_fill_viridis(discrete = TRUE) +
  labs(y="Number of Nodules (GBLUPs)", x= "Strain Treatment") +
  ggtitle("Unstructured Model -A matrix") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 












#***********************************************************************************************
#-----------------------------------------------------------------------
#-----------------TRAIT 2: Number of Nodules  ----------------
#----------------------------------------------------------------------

#Basic linear-mixed model with line as a random factor
qqPlot(Traits.nocontrol$NumberOfNodules)
Traits.nocontrol$sqr.NN <- sqrt(Traits.nocontrol$NumberOfNodules)
Traits.nocontrol$arcsqr.NN <- asin(sqrt(Traits.nocontrol$NumberOfNodules))
qqPlot(Traits.nocontrol$sqr.NN)

#looking the distribution of the data:
bp(Traits.nocontrol, "Treatment", "NumberOfNodules", "Genepool")
bp(Traits.nocontrol, "Treatment", "sqr.NN", "Line")


#Let's find out which ones are the outliers
o <- outliers.fun(Traits.nocontrol, "sqr.NN", "Treatment", "Line", "Genepool")

Traits.nocontrol.noout <- Traits.nocontrol[!Traits.nocontrol$sqr.NN %in% o$sqr.NN,]
length(Traits.nocontrol.noout$sqr.NN)



#####################################################
# Using sommer package - which uses the GBLUPs
#####################################################      

   
#Model A: constraining genetic variation to be same across treatments
NN.CSmodel.1 <- mmer(sqr.NN ~ Treatment  + Genepool + DaysSinceInoculation + HarvestWeek,
                   random = ~ Line + Treatment:Line + Treatment:Block,
                   rcov = ~ units,
                   data = Traits.nocontrol)

NN.CSmodel.0 <- mmer(sqr.NN ~ Treatment  + Genepool + DaysSinceInoculation + HarvestWeek,
                     random = ~ Line + Treatment:Line,
                     rcov = ~ units,
                     data = Traits.nocontrol)

NN.CSmodel.0.0 <- mmer(sqr.NN ~ Treatment  + Genepool + DaysSinceInoculation,
                     random = ~ Line + Treatment:Line,
                     rcov = ~ units,
                     data = Traits.nocontrol)

anova(NN.CSmodel.0,NN.CSmodel.1) #Not significant Block effect
anova(NN.CSmodel.0.0,NN.CSmodel.0)

summary(NN.CSmodel.1)
summary(NN.CSmodel.0)
summary(NN.CSmodel.0.0)
plot(NN.CSmodel.1)
plot(NN.CSmodel.0)
plot(NN.CSmodel.0.0)


#final model:
NN.CSmodel <- mmer(sqr.NN ~ Treatment  + Genepool + DaysSinceInoculation,
                     random = ~ Line + Treatment:Line,
                     rcov = ~ units,
                     data = Traits.nocontrol)

summary(NN.CSmodel)
plot(NN.CSmodel)

                    #Representation of this model in lme4:
                    NN.CSmodel.lmer <- lmer(sqr.NN ~ DaysSinceInoculation + HarvestWeek + Treatment + (1|Line) + (1|Treatment:Line), data = Traits.nocontrol.noout)
                    plot(NN.CSmodel.lmer) #looking at heteroscesdacity
                    qqnorm(residuals(NN.CSmodel.lmer))
                    library(influence.ME)
                    infl <- influence(NN.CSmodel.lmer, obs = TRUE)
                    #Calculate Cook's distance:
                    cooks.distance(infl)
                    #Plot Cook's distance:
                    plot(infl, which = "cook")
                    
                    #finding random effects (e.g. blups)
                    ranef(NN.CSmodel.lmer)
                    plot(allEffects(NN.CSmodel.lmer, partial.residuals = TRUE))
                    dotplot(ranef(NN.CSmodel.lmer, condVar = TRUE)) #Visualizing the random effects
                    qqmath(ranef(NN.CSmodel.lmer, condVar = TRUE)) #This also does the the same but it gets the Q-Qplots
                    



#Model B: allowing genetic variation to differ across treatments
NN.HTmodel.1 <-mmer(sqr.NN~Treatment + Genepool + DaysSinceInoculation + HarvestWeek,tolparinv = 1e-03,
                  random=~Line+ vs(ds(Treatment),Line),
                  rcov=~ vs(ds(Treatment),units),
                  data=Traits.nocontrol)


NN.HTmodel.0 <-mmer(sqr.NN~Treatment + Genepool + DaysSinceInoculation,tolparinv = 1e-03,
                  random=~Line+ vs(ds(Treatment),Line),
                  rcov=~ vs(ds(Treatment),units),
                  data=Traits.nocontrol)

anova(NN.HTmodel.0, NN.HTmodel.1) #Better model without HarvestWeek

#final HS model:
NN.HTmodel <-mmer(sqr.NN~Treatment + Genepool + DaysSinceInoculation,tolparinv = 1e-03,
                    random=~Line+ vs(ds(Treatment),Line),
                    rcov=~ vs(ds(Treatment),units),
                    data=Traits.nocontrol)
summary(NN.HTmodel)
plot(NN.HTmodel)



                    #In lmer4 this would be equivalent to:
                    NN.HSmodel.lmer <- lmer(sqr.NN ~ DaysSinceInoculation + Treatment + (Treatment|Line), data = Traits.nocontrol.noout)
                    plot(NN.HSmodel.lmer)
                    ranef(NN.HSmodel.lmer)
                    plot(allEffects(NN.HSmodel.lmer, partial.residuals = TRUE))
                    dotplot(ranef(NN.HSmodel.lmer, condVar = TRUE)) #Visualizing the random effects
                    qqmath(ranef(NN.HSmodel.lmer, condVar = TRUE)) #This also does the the same but it gets the Q-Qplots


#Model B.1:
NN.HTmodel.pedigree <-mmer(sqr.NN~Treatment + Genepool + DaysSinceInoculation,tolparinv = 1e-03,
                                   random=~ Line + vs(ds(Treatment),Line, Gu=A),    #G <- A.mat(x2) could also be used instead
                                   rcov=~ vs(ds(Treatment),units),
                                   data=Traits.nocontrol)
        plot(NN.HTmodel.pedigree)
        
        
#Model C : assuming an unstructure variance, meaning that among the levels of certain factor (i.e. Environments) there’s a covariance struture of a second random effect(i.e. Genotypes)
#Estimates different GxE at each treatment
NN.USmodel <-mmer(sqr.NN ~ Treatment + Genepool + DaysSinceInoculation,tolparinv = 1e-00,
                  random=~ vs(us(Treatment), Line),
                  rcov=~ vs(us(Treatment),units),
                  data=Traits.nocontrol)
summary(NN.USmodel)
plot(NN.USmodel)

#Model C.1. using pedigree information:
NN.USmodel.pedigree <-mmer(sqr.NN ~ Treatment + Genepool + DaysSinceInoculation,tolparinv = 1e-1,
                           random=~ vs(us(Treatment),Line, Gu=A),
                           rcov=~ vs(us(Treatment),units),
                           data=Traits.nocontrol)
summary(NN.USmodel.pedigree)
plot(NN.USmodel.pedigree)



##################  What is the best model?
AIC <- c(NN.CSmodel$AIC,NN.HTmodel$AIC, NN.HTmodel.pedigree$AIC,NN.USmodel$AIC,NN.USmodel.pedigree$AIC) 
Model <- c("Constrained", "Heterogeneous", "Heterogeneous-Amatrix", "Unstructured", "Unstructured-Amatrix")
VarModels <- data.frame(cbind(Model,AIC))
VarModels[order(VarModels$AIC),] #model with the lowest AIC is the best approximation   


#----Performing like-lihood ratio tests among models:

#let's estimate the number of covariate parameters that were estimated:
summary(NN.CSmodel) 
summary(NN.HTmodel)
summary(NN.HTmodel.pedigree)
summary(NN.USmodel)
summary(NN.USmodel.pedigree)
HG.USmodel$PevU


df.a = 3
df.b = 4 # four VG components are now estimated - also 3 residual variance components
df.b1 = 4
df.c = 7

#* Model A vs B
chi.square= (-2*(42.55257) + 2*(46.17664)) #difference between the log-likelihood of both models
DF= abs(df.a - df.b)
pchisq(chi.square, df=DF, lower.tail = F) #Based on this there is among-genotype variance within treatments
          #*conclusion: among-genotype variances differed among treatments
#* Model A vs B.1
chi.square= (-2*(42.55257) + 2*(47.51)) #difference between the log-likelihood of both models
DF = abs(df.a - df.b1)
pchisq(chi.square, df=DF, lower.tail = F) 
          #*conclusion: among-genotype variances differed among treatments, while controling for line relationships
#* Model A vs C
chi.square= (-2*(42.55257) + 2*(54.65)) #difference between the log-likelihood of both models
DF=  abs(df.a - df.c) #number of var comp estimated
sqrt(pchisq(chi.square, df=DF, lower.tail = F)) #because of the lower value, and since this is a marginal test we did a square-root

#* Model B vs C
chi.square= (-2*(46.17664) + 2*(54.65)) #difference between the log-likelihood of both models
DF= abs(df.b - df.c)
pchisq(chi.square, df=DF, lower.tail = F)
          #*conlusion: genotype effects have different variances at each treatment and different correlations between pairs of treatments

#* Model B vs C.1
chi.square= (-2*(46.17664) + 2*(54.6014)) #difference between the log-likelihood of both models
DF= abs(df.b - df.c)
sqrt(pchisq(chi.square, df=DF, lower.tail = F)) 



###########################Calculating heritability:
#sommer produces the output for (V1) the additive genetic variance, (V2) the variance between A and AL, (V10) the residuals (termed "units").
#heritability estimates based on the Unstructured variance model, that includes the relationship matrix:

sum <- summary(NN.USmodel.pedigree)
varcomp <- data.frame(sum$varcomp[1])
length(varcomp$VarComp)  #we have 12 variance components, let's add numbers to then set up the variances used for heritabilities
varcomp$number <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
#now using the pin() function we can access these different variances:
VP.A = pin(NN.USmodel.pedigree, Vp.NN~(V1+V7)) #phenotypic variance within treatment A equals the Vg (among lines) + residual (among-replicates)
VP.AL = pin(NN.USmodel.pedigree, Vp.NN~(V3+V9)) #phenotypic variance within treatment AL
VP.L = pin(NN.USmodel.pedigree, Vp.NN~(V6+V12)) #phenotypic variance within treatment L

VG.A = pin(NN.USmodel.pedigree, Vg.NN~(V1))
VG.AL = pin(NN.USmodel.pedigree, Vg.NN~(V3))
VP.L = pin(NN.USmodel.pedigree, Vg.NN~(V6))

Hb.A = pin(NN.USmodel.pedigree, Hb.NN~V1/(V1+V7)) #Heritability
Hb.AL = pin(NN.USmodel.pedigree, Hb.NN~ V3/(V3+V9)) #Heritability
Hb.L = pin(NN.USmodel.pedigree, Hb.NN~V6/(V6+V12)) #Heritability

#overall VP (sum of all genetic variances)
VP = pin(NN.USmodel.pedigree, Vp.NN~(V1+V7+V3+V9+V6+V12)) #phenotypic variance conditional on fixed effects in the model
VG = pin(NN.USmodel.pedigree, Vg.NN~(V1+V3+V6))
Hb= VG/VP


#Plotting reaction norms with just the raw BLUPS
BLUP.US.A <- as.data.frame(NN.USmodel$U$`A:Line`)   
BLUP.US.AL <- as.data.frame(NN.USmodel$U$`AL:Line`)
BLUP.US.L <- as.data.frame(NN.USmodel$U$`L:Line`)

BLUPS.NumberOfNodules.US <- BLUP.US.A
BLUPS.NumberOfNodules.US$A <- BLUP.US.A$sqr.NN
BLUPS.NumberOfNodules.US$AL <- BLUP.US.AL$sqr.NN
BLUPS.NumberOfNodules.US$L <- BLUP.US.L$sqr.NN
BLUPS.NumberOfNodules.US$Line <- rownames(BLUP.US.A) 

BLUPS.NumberOfNodules.US <-BLUPS.NumberOfNodules.US[, c(2:5)] #subsetting for only the raw BLUPs
BLUPS.HostGrowth.US.short.0
BLUPS.NumberOfNodules.US.short <- gather(BLUPS.NumberOfNodules.US, Treatment, BLUPs, -Line)
BLUPS.NumberOfNodules.US.short$Genepool <- c("2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2",
                                             "2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2",
                                             "2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2")
                                             

BLUPS.NumberOfNodules.US.short$bt.BLUP <- (BLUPS.NumberOfNodules.US.short$BLUPs)^2


NumberOfNodules.BLUPs.plot <- ggplot(BLUPS.NumberOfNodules.US.short, aes(x=Treatment, y=bt.BLUP, group=Line, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  geom_point(shape=21, colour = "black", size=3) +
  geom_text(data=BLUPS.NumberOfNodules.US.short %>% filter(Treatment == "A"), 
  aes(label=Line), hjust=1.2, vjust=0, size=3, check_overlap = F) +
  scale_fill_viridis(discrete = TRUE) +
  labs(y="Number of Nodules (back-transformed BLUPs)", x= "Strain Treatment") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 



#Plotting reaction norms with just the raw BLUPS
BLUP.US.A <- as.data.frame(NN.USmodel.pedigree$U$`A:Line`)   
BLUP.US.AL <- as.data.frame(NN.USmodel.pedigree$U$`AL:Line`)
BLUP.US.L <- as.data.frame(NN.USmodel.pedigree$U$`L:Line`)

BLUPS.NumberOfNodules.US <- BLUP.US.A
BLUPS.NumberOfNodules.US$A <- BLUP.US.A$sqr.NN
BLUPS.NumberOfNodules.US$AL <- BLUP.US.AL$sqr.NN
BLUPS.NumberOfNodules.US$L <- BLUP.US.L$sqr.NN
BLUPS.NumberOfNodules.US$Line <- rownames(BLUP.US.A) 

BLUPS.NumberOfNodules.US <-BLUPS.NumberOfNodules.US[, c(2:5)] #subsetting for only the raw BLUPs

BLUPS.NumberOfNodules.US.short <- gather(BLUPS.NumberOfNodules.US, Treatment, BLUPs, -Line)
BLUPS.NumberOfNodules.US.short$Genepool <- c("2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2",
                                             "2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2",
                                             "2","2","2", "Wild", "Wild", "Wild", "Wild", "Wild", "Wild",
                                             "2","2","1","1", "1", "1", "1", "1", "2")


BLUPS.NumberOfNodules.US.short$bt.BLUP <- (BLUPS.NumberOfNodules.US.short$BLUPs)^2


NodNum.USVar.btGBLUPs.plot <- ggplot(BLUPS.NumberOfNodules.US.short, aes(x=Treatment, y=bt.BLUP, group=Line, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  geom_point(shape=21, colour = "black", size=3) +
  geom_text(data=BLUPS.NumberOfNodules.US.short %>% filter(Treatment == "A"), 
            aes(label=Line), hjust=1.2, vjust=0, size=3, check_overlap = F) +
  scale_fill_viridis(discrete = TRUE) +
  labs(y="Number of Nodules (back-transformed GBLUPs)", x= "Strain Treatment") +
  ggtitle("Unstructured Model -A matrix (Back-transformed GBLUPS)") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 



#let's compare the predicted trait values with the observed values:
Nod.num.raw.means = summarySE(data=Traits.nocontrol, measurevar= "sqr.NN", groupvars = c("Genepool", "Treatment","Line"), na.rm=TRUE)
Nod.num.raw.means$bt.NodNum <- (Nod.num.raw.means$sqr.NN)^2

NodNum.RN.raw <- ggplot(Nod.num.raw.means, aes(x=Treatment, y=bt.NodNum, group=Line, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  #geom_errorbar(aes(ymin=Value-SE, ymax=Value+SE), width=0.05, colour="black") +
  geom_point(shape=21, colour = "black", size=3) +
  #geom_text(aes(label=Line),hjust=0, vjust=0) +
  geom_text(data=Nod.num.raw.means %>% filter(Treatment== "A"), 
  aes(label=Line), hjust=1.2, vjust=0, size=3, check_overlap = F) +
  geom_text(data=Nod.num.raw.means %>% filter(Treatment== "L"), 
            aes(label=Line), hjust=-0.5, vjust=0, size=3, check_overlap = F) +
  labs(y="Number of Nodules (back-transfomed means)", x= "Strain Treatment") +
  ggtitle("Raw Mean Nodule Number - Back-transformed square-root") +
  scale_fill_viridis(discrete = TRUE) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 


FIGURE2 <- ggarrange(NodNum.RN.raw, NodNum.USVar.btGBLUPs.plot,ncol=2, nrow=1, common.legend = TRUE, labels=c("a", "b"))
ggsave("NodNum_ReactNorms.pdf", plot=FIGURE2, device="pdf", scale = 1, width = 40, height = 20, units = "cm", dpi = 300)


#*******************************************************************
# --- Complex models to estimate variance components for genepools
#******************************************************************
#Models are intercept only models

#Model A: Constrained variance model:
#*assumes that the random genotype and genotype by treatment interaction effects are constant. It involved only two genetic parameters; a variance and a correlation
#*residuals and blocks nested within sites have scaled identity variance structures
NN.CSmodel <- mmer(sqr.NN ~ Treatment + DaysSinceInoculation,
                   random = ~ Line + Genepool + Treatment:Genepool + Genepool:Line,
                   rcov = ~ units,
                   data = Traits.nocontrol)

summary(NN.CSmodel)
plot(NN.CSmodel)


#Model B: heterogeneous variance model - testing if the among-genotype variances differed among treatments, while keeping
#*
NN.HSmodel <- mmer(sqr.NN ~ Treatment + DaysSinceInoculation,
                   random = ~ Line + Genepool + vs(ds(Treatment),ds(Genepool), Line),
                   rcov = ~ vs(ds(Treatment),ds(Genepool),units),
                   data = Traits.nocontrol)


summary(NN.HSmodel)   
plot(NN.HSmodel)


#Model B.1: heterogeneous variance model including relationship matrix
NN.HSmodel.pedigree <- mmer(sqr.NN ~ Treatment + DaysSinceInoculation,
                   random = ~ Line + Genepool + vs(ds(Treatment), ds(Genepool), Line, Gu=A),
                   rcov = ~ vs(ds(Treatment), ds(Genepool),units),
                   data = Traits.nocontrol)

summary(NN.HSmodel.pedigree)
plot(NN.HSmodel.pedigree)

#Model C: Unstructured variance model
#*We can relax a uniform G structure by allowing different genetic (line) variances at each treatment and genepool
NN.USmodel <-mmer(sqr.NN ~ Treatment + DaysSinceInoculation,tolparinv = 1e-03,
                           random=~ vs(us(Treatment),us(Genepool),Line),
                           rcov=~ vs(us(Treatment),units),
                           data=Traits.nocontrol)
summary(NN.USmodel)
plot(NN.USmodel)

#Model C.1: Unstructured-with Pedigree 
NN.USmodel.pedigree <-mmer(sqr.NN ~ Treatment + DaysSinceInoculation,tolparinv = 1e-00,
                           random=~ vs(us(Treatment),us(Genepool),Line, Gu=A),
                           rcov=~ vs(us(Treatment),units),
                           data=Traits.nocontrol)
summary(NN.USmodel.pedigree)
plot(NN.USmodel.pedigree) #(Singular system)


#----Performing like-lihood ratio tests among models:

#let's estimate the number of covariate parameters that were estimated:
summary(NN.CSmodel) 
summary(NN.HSmodel)
summary(NN.HSmodel.pedigree)
summary(NN.USmodel)
summary(NN.USmodel.pedigree)
HG.USmodel$PevU

df.a = 5
df.b = 20 # four VG components are now estimated - also 3 residual variance components
df.b1 = 20
df.c = 54

#* Model A vs B
chi.square= (-2*(42.28196) + 2*(79.68761)) #difference between the log-likelihood of both models
DF= abs(df.a - df.b)
pchisq(chi.square, df=DF, lower.tail = F) #Based on this there is among-genotype variance within treatments
#*conclusion: among-genotype variances differed among treatments
#* Model A vs B.1
chi.square= (-2*(42.55257) + 2*(84.6886)) #difference between the log-likelihood of both models
DF = abs(df.a - df.b1)
pchisq(chi.square, df=DF, lower.tail = F) 
#*conclusion: among-genotype variances differed among treatments, while controling for line relationships

#* Model A vs C
chi.square= (-2*(42.55257) + 2*(64.15616)) #difference between the log-likelihood of both models
DF=  abs(df.a - df.c) #number of var comp estimated
sqrt(pchisq(chi.square, df=DF, lower.tail = F)) #because of the lower value, and since this is a marginal test we did a square-root


#* Model B vs C
chi.square= (-2*(79.68761) + 2*(64.15616)) #difference between the log-likelihood of both models
DF= abs(df.b - df.c)
pchisq(chi.square, df=DF, lower.tail = F)

anova(NN.HSmodel, NN.USmodel)


#*conlusion: genotype effects does not have different variances at each treatment and different correlations between pairs of treatments

#* Model B vs C.1 -- singular fit in C.1


###########################Calculating heritability:
#sommer produces the output for (V1) the additive genetic variance, (V2) the variance between A and AL, (V10) the residuals (termed "units").
#heritability estimates based on the Heterogeneous variance model, that includes the relationship matrix:

sum <- summary(NN.HSmodel.pedigree)
varcomp <- data.frame(sum$varcomp[1])
length(varcomp$VarComp)  #we have 12 variance components, let's add numbers to then set up the variances used for heritabilities
varcomp$number <- c(1:20)
#now using the pin() function we can access these different variances:
VP.A.genepool1 = pin(NN.HSmodel.pedigree, Vp.NN~(V3+V12)) #phenotypic variance within treatment A equals the Vg (among lines) + residual (among-replicates)
VP.A.genepool2 = pin(NN.HSmodel.pedigree, Vp.NN~(V4+V13)) #phenotypic variance within treatment A equals the Vg (among lines) + residual (among-replicates)
VP.A.wild = pin(NN.HSmodel.pedigree, Vp.NN~(V5+V14)) #phenotypic variance within treatment A equals the Vg (among lines) + residual (among-replicates)

VP.AL.genepool1 = pin(NN.HSmodel.pedigree, Vp.NN~(V6+V15)) #phenotypic variance within treatment AL
VP.AL.genepool2 = pin(NN.HSmodel.pedigree, Vp.NN~(V7+V16)) #phenotypic variance within treatment AL
VP.AL.wild = pin(NN.HSmodel.pedigree, Vp.NN~(V8+V17)) #phenotypic variance within treatment AL

VP.L.genepool1 = pin(NN.HSmodel.pedigree, Vp.NN~(V9+V18)) #phenotypic variance within treatment L
VP.L.genepool2 = pin(NN.HSmodel.pedigree, Vp.NN~(V10+V19)) #phenotypic variance within treatment L
VP.L.wild = pin(NN.HSmodel.pedigree, Vp.NN~(V11+V20)) #phenotypic variance within treatment L


VG.A.genepool1 = pin(NN.HSmodel.pedigree, Vg.NN~(V3))
VG.A.genepool2 = pin(NN.HSmodel.pedigree, Vg.NN~(V4))
VG.A.wild = pin(NN.HSmodel.pedigree, Vg.NN~(V5))

VG.AL.genepool1 = pin(NN.USmodel.pedigree, Vg.NN~(V6))
VG.AL.genepool2 = pin(NN.USmodel.pedigree, Vg.NN~(V7))
VG.AL.wild = pin(NN.USmodel.pedigree, Vg.NN~(V8))

VG.L.genepool1 = pin(NN.USmodel.pedigree, Vg.NN~(V9))
VG.L.genepool2 = pin(NN.USmodel.pedigree, Vg.NN~(V10))
VG.L.wild = pin(NN.USmodel.pedigree, Vg.NN~(V11))

Hb.A.genepool1 = pin(NN.HSmodel.pedigree, Hb.NN~V3/(V3+V12)) #Heritability
Hb.A.genepool2 = pin(NN.HSmodel.pedigree, Hb.NN~V4/(V4+V13)) #Heritability
Hb.A.genepool3 = pin(NN.HSmodel.pedigree, Hb.NN~V5/(V5+V14)) #Heritability

Hb.AL.genepool1 = pin(NN.HSmodel.pedigree, Hb.NN~V6/(V6+V15)) #Heritability
Hb.AL.genepool2 = pin(NN.HSmodel.pedigree, Hb.NN~V7/(V7+V16)) #Heritability
Hb.AL.genepool3 = pin(NN.HSmodel.pedigree, Hb.NN~V8/(V8+V17)) #Heritability

Hb.L.genepool1 = pin(NN.HSmodel.pedigree, Hb.NN~V9/(V9+V18)) #Heritability
Hb.L.genepool2 = pin(NN.HSmodel.pedigree, Hb.NN~V10/(V10+V19)) #Heritability
Hb.L.genepool3 = pin(NN.HSmodel.pedigree, Hb.NN~V11/(V11+V20)) #Heritability


#overall VP (sum of all genetic variances)
VP = pin(NN.USmodel.pedigree, Vp.NN~(V1+V7+V3+V9+V6+V12)) #phenotypic variance conditional on fixed effects in the model
VG = pin(NN.USmodel.pedigree, Vg.NN~(V1+V3+V6))
Hb= VG/VP


#Obtaining BLUPS from heterogeneous variance model with A-matrix
c <- data.frame(Traits.nocontrol$Line)
c$Genepool <- Traits.nocontrol$Genepool
c <- unique.data.frame(c)
c.genepool1 <- subset(c, Genepool == "One")
c.genepool2 <- subset(c, Genepool == "Two")
c.wild <- subset(c, Genepool == "Wild")

blups.A.one <- data.frame(NN.HSmodel.pedigree$U$`A:One:Line`)
blups.A.one$Lines <- rownames(blups.A.one)
blups.A.one2 <-blups.A.one[blups.A.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.A.one2$Genepool <- c("One")

blups.A.Two <- data.frame(NN.HSmodel.pedigree$U$`A:Two:Line`)
blups.A.Two$Lines <- rownames(blups.A.Two)
blups.A.Two2 <-blups.A.Two[blups.A.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.A.Two2$Genepool <- c("Two")

blups.A.Wild <- data.frame(NN.HSmodel.pedigree$U$`A:Wild:Line`)
blups.A.Wild$Lines <- rownames(blups.A.Wild)
blups.A.Wild2 <-blups.A.Wild[blups.A.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.A.Wild2$Genepool <- c("Wild")

blups.A <- rbind(blups.A.one2,blups.A.Two2,blups.A.Wild2)
blups.A$Treatment <- c("A")

blups.AL.one <- data.frame(NN.HSmodel.pedigree$U$`AL:One:Line`)
blups.AL.one$Lines <- rownames(blups.AL.one)
blups.AL.one2 <-blups.AL.one[blups.AL.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.AL.one2$Genepool <- c("One")

blups.AL.Two <- data.frame(NN.HSmodel.pedigree$U$`AL:Two:Line`)
blups.AL.Two$Lines <- rownames(blups.AL.Two)
blups.AL.Two2 <-blups.AL.Two[blups.AL.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.AL.Two2$Genepool <- c("Two")

blups.AL.Wild <- data.frame(NN.HSmodel.pedigree$U$`AL:Wild:Line`)
blups.AL.Wild$Lines <- rownames(blups.AL.Wild)
blups.AL.Wild2 <-blups.AL.Wild[blups.AL.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.AL.Wild2$Genepool <- c("Wild")

blups.AL <- rbind(blups.AL.one2,blups.AL.Two2,blups.AL.Wild2)
blups.AL$Treatment <- c("AL")

blups.L.one <- data.frame(NN.HSmodel.pedigree$U$`L:One:Line`)
blups.L.one$Lines <- rownames(blups.L.one)
blups.L.one2 <-blups.L.one[blups.L.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.L.one2$Genepool <- c("One")

blups.L.Two <- data.frame(NN.HSmodel.pedigree$U$`L:Two:Line`)
blups.L.Two$Lines <- rownames(blups.L.Two)
blups.L.Two2 <-blups.L.Two[blups.L.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.L.Two2$Genepool <- c("Two")

blups.L.Wild <- data.frame(NN.HSmodel.pedigree$U$`L:Wild:Line`)
blups.L.Wild$Lines <- rownames(blups.L.Wild)
blups.L.Wild2 <-blups.L.Wild[blups.L.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.L.Wild2$Genepool <- c("Wild")

blups.L <- rbind(blups.L.one2,blups.L.Two2,blups.L.Wild2)
blups.L$Treatment <- c("L")

blups.NN <- rbind(blups.L,blups.AL, blups.A)
blups.NN$btBLUPs <- (blups.NN$sqr.NN)^2

#Plotting reaction norms with just the raw BLUPS
NodNum.USVar.GBLUPs.plot <- ggplot(blups.NN, aes(x=Treatment, y=btBLUPs, group=Lines, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  geom_point(shape=21, colour = "black", size=3) +
  geom_text(data=blups.NN %>% dplyr::filter(Treatment == "L"), 
          aes(label=Lines), hjust=-0.5, vjust=0, size=3, check_overlap = F) +
  scale_fill_viridis(discrete = TRUE) +
  labs(y="Number of Nodules (GBLUPs)", x= "Strain Treatment") +
  ggtitle("Unstructured Model -A matrix") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 


#Obtaining BLUPS from heterogeneous variance model without A-matrix
c <- data.frame(Traits.nocontrol$Line)
c$Genepool <- Traits.nocontrol$Genepool
c <- unique.data.frame(c)
c.genepool1 <- subset(c, Genepool == "One")
c.genepool2 <- subset(c, Genepool == "Two")
c.wild <- subset(c, Genepool == "Wild")

blups.A.one <- data.frame(NN.HSmodel$U$`A:One:Line`)
blups.A.one$Lines <- rownames(blups.A.one)
blups.A.one2 <-blups.A.one[blups.A.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.A.one2$Genepool <- c("One")

blups.A.Two <- data.frame(NN.HSmodel$U$`A:Two:Line`)
blups.A.Two$Lines <- rownames(blups.A.Two)
blups.A.Two2 <-blups.A.Two[blups.A.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.A.Two2$Genepool <- c("Two")

blups.A.Wild <- data.frame(NN.HSmodel$U$`A:Wild:Line`)
blups.A.Wild$Lines <- rownames(blups.A.Wild)
blups.A.Wild2 <-blups.A.Wild[blups.A.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.A.Wild2$Genepool <- c("Wild")

blups.A <- rbind(blups.A.one2,blups.A.Two2,blups.A.Wild2)
blups.A$Treatment <- c("A")

blups.AL.one <- data.frame(NN.HSmodel$U$`AL:One:Line`)
blups.AL.one$Lines <- rownames(blups.AL.one)
blups.AL.one2 <-blups.AL.one[blups.AL.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.AL.one2$Genepool <- c("One")

blups.AL.Two <- data.frame(NN.HSmodel$U$`AL:Two:Line`)
blups.AL.Two$Lines <- rownames(blups.AL.Two)
blups.AL.Two2 <-blups.AL.Two[blups.AL.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.AL.Two2$Genepool <- c("Two")

blups.AL.Wild <- data.frame(NN.HSmodel$U$`AL:Wild:Line`)
blups.AL.Wild$Lines <- rownames(blups.AL.Wild)
blups.AL.Wild2 <-blups.AL.Wild[blups.AL.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.AL.Wild2$Genepool <- c("Wild")

blups.AL <- rbind(blups.AL.one2,blups.AL.Two2,blups.AL.Wild2)
blups.AL$Treatment <- c("AL")

blups.L.one <- data.frame(NN.HSmodel$U$`L:One:Line`)
blups.L.one$Lines <- rownames(blups.L.one)
blups.L.one2 <-blups.L.one[blups.L.one$Lines %in% c.genepool1$Traits.nocontrol.Line,]
blups.L.one2$Genepool <- c("One")

blups.L.Two <- data.frame(NN.HSmodel$U$`L:Two:Line`)
blups.L.Two$Lines <- rownames(blups.L.Two)
blups.L.Two2 <-blups.L.Two[blups.L.Two$Lines %in% c.genepool2$Traits.nocontrol.Line,]
blups.L.Two2$Genepool <- c("Two")

blups.L.Wild <- data.frame(NN.HSmodel$U$`L:Wild:Line`)
blups.L.Wild$Lines <- rownames(blups.L.Wild)
blups.L.Wild2 <-blups.L.Wild[blups.L.Wild$Lines %in% c.wild$Traits.nocontrol.Line,]
blups.L.Wild2$Genepool <- c("Wild")

blups.L <- rbind(blups.L.one2,blups.L.Two2,blups.L.Wild2)
blups.L$Treatment <- c("L")

blups.NN <- rbind(blups.L,blups.AL, blups.A)
blups.NN$btBLUPs <- (blups.NN$sqr.NN)^2

#Plotting reaction norms with just the raw BLUPS
NN.HSVar.BLUPs.plot <- ggplot(blups.NN, aes(x=Treatment, y=btBLUPs, group=Lines, fill=Genepool)) +
  geom_line(linetype = "solid", size =0.5) +
  geom_point(shape=21, colour = "black", size=3) +
  geom_text(data=blups.NN %>% dplyr::filter(Treatment == "L"), 
            aes(label=Lines), hjust=-0.5, vjust=0, size=3, check_overlap = F) +
  scale_fill_viridis(discrete = TRUE) +
  labs(y="Number of Nodules (BLUPs)", x= "Strain Treatment") +
  ggtitle("Heterogeneous Model") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 





#----------------------------------------------------------------------------------------
#-------------------------------TRAIT 3: Investment  ----------------
#-----------------------------------------------------------------------------------------

#Basic linear-mixed model with line as a random factor
qqPlot(Traits.nocontrol$Investment)
Traits.nocontrol$logInvestment <- log(Traits.nocontrol$Investment + 1)
qqPlot(Traits.nocontrol$logInvestment) # same as without transformation
length(Traits.nocontrol$Investment) 

#Boxplot to visualize outliers:
boxp.Invest <- bp(Traits.nocontrol, "Treatment", "Investment", "Line")
outliers.fun(Traits.nocontrol, "NumberOfNodules", "Treatment", "Line", "Genepool")

library(sommer)
#Model A: constraining genetic variation to be same across treatments
        Invest.CSmodel <- mmer(Investment ~ Treatment + Genepool,
                           random = ~ Line + Block + Treatment:Line + Genepool:Line,
                           rcov = ~ units,
                           data = Traits.nocontrol)
        
        summary(Invest.CSmodel)
        plot(Invest.CSmodel) #lots of heteroscesdacity
        
        #let's try transforming the data:
        Traits.nocontrol$log.invest <- log(Traits.nocontrol$Investment)
        Traits.nocontrol$sqrt.invest <- sqrt(Traits.nocontrol$Investment)
        Traits.nocontrol$asinsqrt.invest <- asin(sqrt(Traits.nocontrol$Investment))
        
        Invest.CSmodel <- mmer(asinsqrt.invest  ~ Treatment + Genepool + DaysSinceInoculation + HarvestWeek,
                               random = ~ Line + Treatment:Line + Block,
                               rcov = ~ units,
                               data = Traits.nocontrol)
        
        summary(Invest.CSmodel)
        plot(Invest.CSmodel) #sqrt improved the Heteroscesdacity
        
        
        
#Model B: allowing genetic variation to differ across treatments
        
        Invest.HTmodel <-mmer(asinsqrt.invest~Treatment + Genepool + DaysSinceInoculation + HarvestWeek,tolparinv = 1e-03,
                          random=~ Block + Line + vs(ds(Treatment),Line),
                          rcov=~ vs(ds(Treatment),units),
                          data=Traits.nocontrol)
        summary(Invest.HTmodel)
        plot(Invest.HTmodel)
        
        
        
        
        #Model B.1:
        Invest.HTmodel.pedigree <-mmer(asinsqrt.invest~Treatment + Genepool +  DaysSinceInoculation + HarvestWeek ,tolparinv = 1e-03,
                                   random=~Line+ Block + vs(ds(Treatment),Line, Gu=A),    #G <- A.mat(x2) could also be used instead
                                   rcov=~ vs(ds(Treatment),units),
                                   data=Traits.nocontrol)
        anova(Invest.HTmodel.pedigree,Invest.CSmodel) 
        anova(Invest.HTmodel.pedigree, Invest.HTmodel) 
        plot(Invest.HTmodel.pedigree)
        summary(Invest.HTmodel.pedigree)

#Model C : assuming an unstructure variance, meaning that among the levels of certain factor (i.e. Environments) there’s a covariance struture of a second random effect(i.e. Genotypes)
        #Estimates different GxE at each treatment
        Invest.USmodel <-mmer(asinsqrt.invest ~ Treatment + Genepool + DaysSinceInoculation + HarvestWeek,tolparinv = 1e-03,
                          random=~ Block + vs(us(Treatment),Line),
                          rcov=~ vs(us(Treatment),units),
                          data=Traits.nocontrol)
        summary(Invest.USmodel)
        plot(Invest.USmodel)        
        anova(Invest.HTmodel,Invest.USmodel)
        

        
