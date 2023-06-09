##################################################################
# Analyses of 1536 SNPs from Hyunh et al. 2013 
##################################################################

library(adegenet)
library(ggpubr)
library(vcfR)
library(phangorn)
library(ape)
library(poppr)
library(adegenet) 
library(adegraphics) 
library(pegas) 
library(StAMPP) 
library(lattice) 
library(gplots) 
library(ape) 
library(ggmap) 
library(ape)
library(adegenet)
library(poppr)

# Reading data:
SNPS <- read.csv("../Hyunh_GoldenGate_SNPs.csv", header = TRUE, sep=",", strip.white = FALSE, stringsAsFactors=F)
str(SNPS)
metadata <- read.csv("../Hyunh_GoldenGate_Metadata.csv", header = TRUE, sep=",", strip.white = TRUE, stringsAsFactors=F)
str(metadata)
factormetadata$Genepool

# Removing SNPs that failed genotyping & were monomorphic in the wild lines:
#SNPS.f <- subset(SNPS, Status_W != "Failed" & Status_W != "Mono")

SNPS.f <- SNPS[, c(1:3, 8:446)] #removing some extra columns, but keeping chromosome positions
str(SNPS.f)
SNPS.c <- SNPS.f[, c(1, 4:442)]
str(SNPS.c)

#Transposing matrix:
rownames(SNPS.c ) <- SNPS.c$SNP  #The first column contains the SNPs, or what we want it to be the rownames
SNPS.c$SNP <- NULL
SNPS.c.t <- as.data.frame(t(as.matrix(SNPS.c)))

SNP.simple.t2 <- as.matrix(SNPS.c.t)
SNP.simple.t2[6:20,3:56] #checking matrix
SNP.simple.t2[ SNP.simple.t2 == "--" ] <- NA  #Need to remove the "--" and convert into NAs
SNP.simple.t2[1:20,3:10] #checking matrix

#Recoding can also be done with the library snpReady, the good thing here is that it help us filter the SNPs prior to analyses
#Now we need to filter some of the SNPs that have more than 50% missing data
library(snpReady)
geno.ready <- raw.data(data = as.matrix(SNP.simple.t2), frame = "wide", base=TRUE, sweep.sample = 0.8, call.rate = 0.90, maf = 0.05, imput = FALSE)
geno.ready$report #This let us know which SNPs were removed
length(colnames(geno.ready$M.clean))  #753 left after filtering for missing data and minimum allele frequencies
x <- geno.ready$M.clean
x[10:5,34:5]  #checking matrix

#Using a genlight object from adegenet to do genetic analyses
library(sf)
library(adegenet)
x1 <- new("genlight", x, parallel=FALSE)

#Assign genepools and line names
pop(x1) <- metadata$Genepool
indNames(x1) <- metadata$Line


#Exporting genlight to faststructure; analyses run in cluster see faststructure.sh
library(dartR)
ploidy(x1) <- c("2")  #ploidy needs to be set for the transformation to work
gl2faststructure(x1,outfile = "Cowpea_faststr.str",outpath = ".",probar = FALSE,verbose = NULL)

gl2structure(x1, indNames = NULL, addcolumns = NULL,
  ploidy = 2,
  exportMarkerNames = TRUE,
  outfile = "gl.str",
  outpath = ".",
  verbose = NULL
)

#converting genlight to genind and then to STRUCTURE format
x1 #genlight object, needs to be transformed into a gind object
gind <- gl2gi(x1)  #converting genlight object to genind

genind2structure(gind, file="cowpeas.str", pops=FALSE)



#Making NEIGHBOR-JOINING AND UPGMA  TREES

library(dartR)
x1 #genlight object, needs to be transformed into a gind object
gind <- gl2gi(x1)  #converting genlight object to genind

tree.upgma <- poppr::aboot(gind, tree = "upgma", distance = nei.dist, sample = 100, ploidy=2 , showtree = F, cutoff = 50, quiet = T)
tree.NJ <- poppr::aboot(gind, tree = "nj", distance = nei.dist, sample = 100, ploidy=2 , showtree = F, cutoff = 50, quiet = T)
root_treeNJ <- midpoint(tree.NJ)

BiocManager::install("ggtree")
library(dplyr)
library(ggtree)

install.packages("viridis")  # Install
library("viridis")           # Load
cols <- c("gray","#440154FF","#21908CFF","#FDE725FF")

rownames(metadata_sorted) <- metadata_sorted$Line

upgmatree.nei <- ggtree(tree.upgma, layout='rectangular') %<+% metadata +
  geom_tiplab(aes(label=label), align=FALSE, size=4) +
  #geom_text2(aes(subset=!isTip, label=Isolate.name),size=2, align=TRUE) +
  geom_tippoint(aes(color=Genepool), alpha=.9, size=2,show.legend = TRUE) + scale_colour_manual(values=cols) +
  #scale_color_manual(values=species.col, name="Host species") +
  geom_treescale(fontsize=3, offset=2, y=20) +
  #geom_label_repel(aes(label=bootstrap, fill=bootstrap)) #option when using the branchlenght tree
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label))), hjust=-.2, size=2.5) # & as.numeric(label) > 80))  #option if filtering bootstraps and using newick format

njtree.nei <- ggtree(tree.NJ, layout='rectangular') %<+% metadata +
  geom_tiplab(aes(label=label), align=FALSE, size=2) +
  #geom_text2(aes(subset=!isTip, label=Isolate.name),size=2, align=TRUE) +
  geom_tippoint(aes(color=Genepool), alpha=.9, size=1,show.legend = TRUE) + scale_colour_manual(values=cols) +
  #scale_color_manual(values=species.col, name="Host species") +
  geom_treescale(fontsize=3, offset=2, y=20) +
  #geom_label_repel(aes(label=bootstrap, fill=bootstrap)) #option when using the branchlenght tree
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label))), hjust=-.2, size=2.5) # & as.numeric(label) > 80))  #option if filtering bootstraps and using newick format


ggsave("upgmatree.nei.pdf", plot=upgmatree.nei, device="pdf", scale = 1, width = 30, height = 60, units = "cm", dpi = 600)
ggsave("njtree.nei.pdf", plot=njtree.nei, device="pdf", scale = 1, width = 30, height = 60, units = "cm", dpi = 600)

#ploting midpoint rooted tree:

njtree_rooted.nei <- ggtree(root_treeNJ, layout='fan') %<+% metadata +
  geom_tiplab(aes(label=label), align=FALSE, size=1, offset = 0.02 ) +
  #geom_text2(aes(subset=!isTip, label=Isolate.name),size=2, align=TRUE) +
  geom_tippoint(aes(color=Genepool), alpha=.9, size=1,show.legend = TRUE) + scale_colour_manual(values=cols) +
  #scale_color_manual(values=species.col, name="Host species") +
  geom_treescale(fontsize=3, offset=2) 
  #geom_label_repel(aes(label=bootstrap, fill=bootstrap)) #option when using the branchlenght tree
  #geom_text2(aes(label=label, subset = !is.na(as.numeric(label))), hjust=-.2, size=2.5) # & as.numeric(label) > 80))  #option if filtering bootstraps and using newick format


ggsave("njtree.nei_root_circle.pdf", plot=njtree_rooted.nei, device="pdf", scale = 1, width = 30, height = 30, units = "cm", dpi = 600)


#Making Principal Component Analyses -  PCA
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
pop80.pca.scores$Sachs_Lab <- metadata$Sachs_Lab

library(ggplot2)
library(viridis)
set.seed(9)
#cols <- viridis(n = 2)
#cols <- viridis(n = nPop(x1))
cols <- c("gray","#440154FF","#21908CFF","#FDE725FF")

p <- ggplot(pop80.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
#p <- p + geom_text(aes(label=ind), hjust=-0.5, vjust=0, size=3, position=position_jitter(width=0,height=3),check_overlap = F) 
#p <- p + ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = ind))
#p <- p + geom_text(data=subset(pop80.pca.scores, Sachs_Lab == "YES"), label=ind)
p <- p + ggrepel::geom_text_repel(aes(label=ifelse(Sachs_Lab == "YES",as.character(ind),'')),position = "identity", max.overlaps=100, size=3, color="black")
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p
pca.cowpea <- p + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                        axis.text.x = element_text(size=16, colour ="black"),
                        axis.text.y = element_text(size=16, colour ="black"),
                        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
                        legend.position = "top")
                        


ggsave("pca.cowpeav.pdf", plot=pca.cowpea, device="pdf", scale = 1, width = 20, height = 20, units = "cm", dpi = 600)


p <- ggplot(pop80.pca.scores, aes(x=PC1, y=PC3, colour=pop)) 
p <- p + geom_point(size=2)
#p <- p + geom_text(aes(label=ind), hjust=-0.5, vjust=0, size=3, position=position_jitter(width=0,height=3),check_overlap = F) 
#p <- p + ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = ind))
#p <- p + geom_text(data=subset(pop80.pca.scores, Sachs_Lab == "YES"), label=ind)
p <- p + ggrepel::geom_text_repel(aes(label=ifelse(Sachs_Lab == "YES",as.character(ind),'')),position = "identity", max.overlaps=100, size=3, color="black")
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p
pca.cowpea.2 <- p + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                        axis.title.x = element_text(face="bold", vjust=1.0, size = 16),
                        axis.title.y = element_text(face="bold", vjust=1.0, size = 16),
                        axis.text.x = element_text(size=16, colour ="black"),
                        axis.text.y = element_text(size=16, colour ="black"),
                        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
                        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
                        legend.position = "top")


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
mean(Ho$`Landrace-Genepool2`)
se(Ho$`Landrace-Genepool2`)

mean(Ho$`Landrace-Genepool1`)
se(Ho$`Landrace-Genepool1`)

mean(Ho$Wild)
se(Ho$Wild)

Hs <- as.data.frame(x.stats$Hs)
mean(Hs$`Landrace-Genepool2`)
se(Hs$`Landrace-Genepool2`)

mean(Hs$`Landrace-Genepool1`)
se(Hs$`Landrace-Genepool1`)

mean(Hs$Wild)
se(Hs$Wild)



#Gene diversity (Hs): -- which is is the chance that if you take two random alleles from your sample they are different
library(tidyr)
Hs <- x.stats$Hs
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
Ho <- x.stats$Ho
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


####################################################
# Testing genetic differentiation among lineages
####################################################
FST <- genet.dist(x.gind2, method = "WC84")
fst.pairwise <- pairwise.WCfst(x.gind2,diploid=TRUE) #Estimates pairwise FSTs according to Weir and Cockerham (1984
write.csv(fst.pairwise, "fst.pairwise.csv")
boot.fst <- boot.ppfst(dat=x.gind2,nboot=100,quant=c(0.025,0.975),diploid=TRUE)
write.csv(boot.fst, "boot.fst.csv")


#########################################
# Running LEA
##########################################

library(LEA)
#Creating a .geno file to work with LEA:
#note: the file was .ped format and was created using the program Link_Impute

input = "./cowpeas.str"
output = "./"
output = struct2geno(input, ploidy=2, FORMAT = 2, extra.column = 1, extra.row = 1 )


#TESTING ADMIXTURE & GENETIC CLUSTERS
# main options
# K = number of ancestral populations, here 15
# entropy = TRUE: computes the cross-entropy criterion,
# CPU = 4 the number of CPUs.
#Note: run this in the cluster not in your own pc
set.seed=800
project = NULL
project = snmf("cowpeas.str.geno",K = 1:10,entropy = TRUE,repetitions = 100,project = "new", CPU=4)

#project = load.snmfProject("cowpeas.str.snmfProject") #if wanting to read the project after it was created

#plot cross-entropy criterion for all runs in the snmf project
#Naming the file
pdf(file = "crossentropy.pdf")
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off() #Saving the file


#select the best run for K = 3
best = which.min(cross.entropy(project, K = 3))


best = which.min(cross.entropy(project, K = 3))
qmatrix = Q(project, K = 3, run = best)
qmatrix
cluster <- apply(qmatrix, 1, which.max)  #getting cluster assignments for each individual


#metadata <- read.csv("metadata_sorted.csv", header = TRUE, sep=",", strip.white = TRUE, stringsAsFactors=F)
#str(metadata)
metadata <- read.csv("Hyunh_GoldenGate_Metadata.csv", header = TRUE, sep=",", strip.white = TRUE, stringsAsFactors=F)
str(metadata)
names_ind <- metadata$X...Line
new_qmatrix <- cbind(cluster, qmatrix) #add these assignments to sort ind by this
new_qmatrix <- cbind(names_ind, new_qmatrix)
new_qmatrix_sorted <- new_qmatrix[order(new_qmatrix[,2], decreasing=TRUE),] #sorting by cluster assignment
#new_qmatrix_sorted <- new_qmatrix_sorted[order(new_qmatrix_sorted[,3], decreasing=FALSE),]
new_qmatrix_sorted <- new_qmatrix_sorted[, c(3:5)] #removing the added columns used to sort
write.csv(new_qmatrix_sorted, "new_qmatrix_sorted_K3.csv")

metadata$Assignment.K3 <- cluster
metadata2 <-metadata[order(metadata$Assignment.K3),]
names_sample <- metadata2$X...Line



cols <- c("#440154FF","#FDE725FF", "#21908CFF")
#cols <- c("#FDE725FF", "#21908CFF","#440154FF")
pdf(file = "LEA_clusters_QmatK3.pdf", width = 60, height=10)
barplot(t(new_qmatrix_sorted), 
        col = cols,
        border=NA,
        space=0.0, #to reduce space between barplots
        names.arg = names_sample,
        cex.names=0.6, 
        horiz = FALSE,
        las = 2,
        #xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix") 
dev.off() 



#Let's sort by the Genepool Assignment in Hyunh:
names_ind <- metadata$X...Line
genepools <- metadata$Genepool
new_qmatrix <- cbind(cluster, qmatrix) #add these assignments to sort ind by this
new_qmatrix <- cbind(names_ind, genepools,new_qmatrix)
new_qmatrix_sorted <- new_qmatrix[order(new_qmatrix[,2], decreasing=FALSE),] #sorting by cluster assignment
#new_qmatrix_sorted <- new_qmatrix_sorted[order(new_qmatrix_sorted[,3], decreasing=FALSE),]
new_qmatrix_sorted <- new_qmatrix_sorted[, c(4:6)] #removing the added columns used to sort
write.csv(new_qmatrix_sorted, "new_qmatrix_sorted_K3.csv")
 
metadata$Assignment.K3 <- cluster
write.csv(metadata, "metadata_k3.csv")
metadata2 <-metadata[order(metadata$Genepool),]
names_sample <- metadata2$X...Line

cols <- c("#440154FF","#FDE725FF", "#21908CFF")
#cols <- c("#FDE725FF", "#21908CFF","#440154FF")
pdf(file = "LEA_clusters_QmatK3_bygenepool.pdf", width = 60, height=10)
barplot(t(new_qmatrix_sorted), 
        col = cols,
        border=NA,
        space=0.0, #to reduce space between barplots
        names.arg = names_sample,
        cex.names=0.6, 
        horiz = FALSE,
        las = 2,
        #xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix") 
dev.off() 


#Let's sort by the Q score and Genepool Assignment in Hyunh:
names_ind <- metadata$X...Line
genepools <- metadata$Genepool
new_qmatrix <- cbind(cluster, qmatrix) #add these assignments to sort ind by this
new_qmatrix <- cbind(names_ind, genepools,new_qmatrix)

new_qmatrix_sorted <- new_qmatrix[order(new_qmatrix[,3],new_qmatrix[,2],qmatrix[,2],decreasing=FALSE),] #sorting by cluster assignment
#new_qmatrix_sorted <- new_qmatrix_sorted[order(new_qmatrix_sorted[,3], decreasing=FALSE),]
new_qmatrix_sorted <- new_qmatrix_sorted[, c(4:6)] #removing the added columns used to sort
#write.csv(new_qmatrix_sorted, "new_qmatrix_sorted_K3.csv")

metadata$Assignment.K3 <- cluster
metadata$V1 <-  qmatrix[,1]
metadata$V2 <-  qmatrix[,2]
metadata$V3 <-  qmatrix[,3]
#write.csv(metadata, "metadata_k3.csv")
metadata2 <-metadata[order(metadata$Assignment.K3,metadata$Genepool,metadata$V2,decreasing = FALSE),]
names_sample <- metadata2$X...Line
#names_sample <- metadata2$Genepool
cols <- c("#440154FF","#FDE725FF", "#21908CFF")
#cols <- c("#FDE725FF", "#21908CFF","#440154FF")
pdf(file = "LEA_clusters_QmatK3_bygenepoolandQ3.pdf", width = 60, height=10)
barplot(t(new_qmatrix_sorted), 
        col = cols,
        border=NA,
        space=0.0, #to reduce space between barplots
        names.arg = names_sample,
        cex.names=0.6, 
        horiz = FALSE,
        las = 2,
        #xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix") 
dev.off() 






#ploting if K=2, sorting by genepools
best = which.min(cross.entropy(project, K = 2))
qmatrix = Q(project, K = 2, run = best)
qmatrix
cluster <- apply(qmatrix, 1, which.max)  #getting cluster assignments for each individual

names_ind <- metadata$X...Line
genepools <- metadata$Genepool

new_qmatrix <- cbind(cluster, qmatrix) #add these assignments to sort ind by this
new_qmatrix <- cbind(names_ind, genepools,new_qmatrix)
new_qmatrix_sorted <- new_qmatrix[order(new_qmatrix[,2], decreasing=FALSE),] #sorting by cluster assignment
#new_qmatrix_sorted <- new_qmatrix_sorted[order(new_qmatrix_sorted[,3], decreasing=FALSE),]
new_qmatrix_sorted <- new_qmatrix_sorted[, c(4:5)] #removing the added columns used to sort

metadata$Assignment.K2 <- cluster
metadata2 <-metadata[order(metadata$Genepool),]
names_sample <- metadata2$X...Line



cols <- viridis(n = 3)
cols <- c("#440154FF","#FDE725FF")
pdf(file = "LEA_clusters_QmatK2_genepools.pdf", width = 60, height=10)
barplot(t(new_qmatrix_sorted), 
        col = cols,
        border=NA,
        space=0.0, #to reduce space between barplots
        names.arg = names_sample,
        cex.names=0.6, 
        horiz = FALSE,
        las = 2,
        #xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix") 
dev.off() 

#Ploting K=4
cols <- viridis(n = 4)
best = which.min(cross.entropy(project, K = 4))
qmatrix = Q(project, K = 4, run = best)
qmatrix
cluster <- apply(qmatrix, 1, which.max)  #getting cluster assignments for each individual

names_ind <- metadata$X...Line
genepools <- metadata$Genepool

new_qmatrix <- cbind(cluster, qmatrix) #add these assignments to sort ind by this
new_qmatrix <- cbind(names_ind, genepools,new_qmatrix)
new_qmatrix_sorted <- new_qmatrix[order(new_qmatrix[,2], decreasing=FALSE),] #sorting by cluster assignment
#new_qmatrix_sorted <- new_qmatrix_sorted[order(new_qmatrix_sorted[,3], decreasing=FALSE),]
new_qmatrix_sorted <- new_qmatrix_sorted[, c(4:7)] #removing the added columns used to sort

metadata$Assignment.K4 <- cluster
metadata2 <-metadata[order(metadata$Genepool),]
names_sample <- metadata2$X...Line



cols <- viridis(n = 4)
#cols <- c("#21908CFF","#440154FF", "#FDE725FF","gray")
#cols <- c("#440154FF","#31688EFF","#FDE725FF","#35B779FF")
pdf(file = "LEA_clusters_QmatK4_genepools.pdf", width = 60, height=10)
barplot(t(new_qmatrix_sorted), 
        col = cols,
        border=NA,
        space=0.0, #to reduce space between barplots
        names.arg = names_sample,
        cex.names=0.6, 
        horiz = FALSE,
        las = 2,
        #xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix") 
dev.off() 


#Ploting K=5
best = which.min(cross.entropy(project, K = 5))
qmatrix = Q(project, K = 5, run = best)
qmatrix
cluster <- apply(qmatrix, 1, which.max)  #getting cluster assignments for each individual

names_ind <- metadata$X...Line
genepools <- metadata$Genepool

new_qmatrix <- cbind(cluster, qmatrix) #add these assignments to sort ind by this
new_qmatrix <- cbind(names_ind, genepools,new_qmatrix)
new_qmatrix_sorted <- new_qmatrix[order(new_qmatrix[,2], decreasing=FALSE),] #sorting by cluster assignment
#new_qmatrix_sorted <- new_qmatrix_sorted[order(new_qmatrix_sorted[,3], decreasing=FALSE),]
new_qmatrix_sorted <- new_qmatrix_sorted[, c(4:8)] #removing the added columns used to sort

metadata$Assignment.K5 <- cluster
metadata2 <-metadata[order(metadata$Genepool),]
names_sample <- metadata2$X...Line



cols <- viridis(n = 5)
#cols <- c("#21908CFF","#440154FF", "#FDE725FF","gray")
#cols <- c("#440154FF","#31688EFF","#FDE725FF","#35B779FF")
pdf(file = "LEA_clusters_QmatK5_genepools.pdf", width = 60, height=10)
barplot(t(new_qmatrix_sorted), 
        col = cols,
        border=NA,
        space=0.0, #to reduce space between barplots
        names.arg = names_sample,
        cex.names=0.6, 
        horiz = FALSE,
        las = 2,
        #xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix") 
dev.off() 

#Ploting K=6
best = which.min(cross.entropy(project, K = 6))
qmatrix = Q(project, K = 6, run = best)
qmatrix
cluster <- apply(qmatrix, 1, which.max)  #getting cluster assignments for each individual

names_ind <- metadata$X...Line
genepools <- metadata$Genepool

new_qmatrix <- cbind(cluster, qmatrix) #add these assignments to sort ind by this
new_qmatrix <- cbind(names_ind, genepools,new_qmatrix)
new_qmatrix_sorted <- new_qmatrix[order(new_qmatrix[,2], decreasing=FALSE),] #sorting by cluster assignment
#new_qmatrix_sorted <- new_qmatrix_sorted[order(new_qmatrix_sorted[,3], decreasing=FALSE),]
new_qmatrix_sorted <- new_qmatrix_sorted[, c(4:9)] #removing the added columns used to sort

metadata$Assignment.K6 <- cluster
metadata2 <-metadata[order(metadata$Genepool),]
names_sample <- metadata2$X...Line



cols <- viridis(n = 6)
#cols <- c("#21908CFF","#440154FF", "#FDE725FF","gray")
#cols <- c("#440154FF","#31688EFF","#FDE725FF","#35B779FF")
pdf(file = "LEA_clusters_QmatK6_genepools.pdf", width = 60, height=10)
barplot(t(new_qmatrix_sorted), 
        col = cols,
        border=NA,
        space=0.0, #to reduce space between barplots
        names.arg = names_sample,
        cex.names=0.6, 
        horiz = FALSE,
        las = 2,
        #xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix") 
dev.off() 


#SUBSETING 20 COWPEAS

best = which.min(cross.entropy(project, K = 3))
qmatrix = Q(project, K = 3, run = best)
qmatrix
cluster <- apply(qmatrix, 1, which.max)  #getting cluster assignments for each individual

metadata <- read.csv("Hyunh_GoldenGate_Metadata.csv", header = TRUE, sep=",", strip.white = TRUE, stringsAsFactors=F)
str(metadata)
names_ind <- metadata$X...Line
new_qmatrix <- cbind(cluster, qmatrix) #add these assignments to sort ind by this
new_qmatrix <- cbind(names_ind, new_qmatrix)
new_qmatrix_sorted <- new_qmatrix[order(new_qmatrix[,2], decreasing=FALSE),] #sorting by cluster assignment
#new_qmatrix_sorted <- new_qmatrix_sorted[order(new_qmatrix_sorted[,3], decreasing=FALSE),]
sachscowpeas <- subset(metadata, Sachs_Lab == "YES")
colnames(new_qmatrix_sorted) <- c("X...Line", "cluster", "V1", "V2", "V3")
new_qmatrix_sachscow <- merge(sachscowpeas,new_qmatrix_sorted, by="X...Line")
new_qmatrix_sachscow2 <- new_qmatrix_sachscow[order(new_qmatrix_sachscow[,5], decreasing=FALSE),]
new_qmatrix_sachscow2 <- new_qmatrix_sachscow[, c(6:8)]

cols <- c("#440154FF","#FDE725FF", "#21908CFF")
#cols <- c("#FDE725FF", "#21908CFF","#440154FF")
names_sample <- new_qmatrix_sachscow$X...Line
pdf(file = "LEA_QmatK3_bygenepool_SachsCowpeas.pdf", width = 18, height=5)
barplot(t(new_qmatrix_sachscow2), 
        col = cols,
        border=NA,
        space=0.0, #to reduce space between barplots
        names.arg = names_sample,
        cex.names=1.0, 
        horiz = FALSE,
        las = 2,
        #xlab = "Individuals",
        ylab = "Ancestry proportions",
        main = "Ancestry matrix") 
dev.off() 
