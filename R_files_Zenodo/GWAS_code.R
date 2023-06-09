library("tidyverse")
library("vcfR")
library("multtest")
library("gplots")
library("LDheatmap")
library("genetics")
library("ape")
library("EMMREML")
library("scatterplot3d")
library("compiler")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")
library("stringr")
library("lme4")
library("car")
library("emmeans")
library("qqman")



######################################
### --- Change vcf into hapmap --- ###
######################################

vcf <- read.vcfR("data/CdiffSNPs_GWAS_10_filt.recode.vcf.gz", verbose = TRUE)

info <- getFIX(vcf)
GT <- extract.gt(vcf, return.alleles=TRUE, IDtoRowNames = FALSE, convertNA = FALSE)

GT <- gsub("A.A", "A", GT)
GT <- gsub("T.T", "T", GT)
GT <- gsub("C.C", "C", GT)
GT <- gsub("G.G", "G", GT)
GT <- gsub("A.G", "R", GT)
GT <- gsub("G.A", "R", GT)
GT <- gsub("C.T", "Y", GT)
GT <- gsub("T.C", "Y", GT)
GT <- gsub("G.C", "S", GT)
GT <- gsub("C.G", "S", GT)
GT <- gsub("A.T", "W", GT)
GT <- gsub("T.A", "W", GT)
GT <- gsub("T.G", "K", GT)
GT <- gsub("G.T", "K", GT)
GT <- gsub("C.A", "M", GT)
GT <- gsub("A.C", "M", GT)
GT <- gsub("\\\\W+", "N", GT)

colnames(GT) <- gsub("\\\\.", "\\\\_", colnames(GT))

rs <- paste(info[,1], info[,2], sep = "_")
alleles <- paste(info[,4], info[,5], sep = "/")
chrom <- info[,1]
pos <- info[,2]
strand <- "+"
assembly <- NA
center <- NA
protLSID <- NA
assayLSID <- NA
panelLSID <- NA
Qcode <- NA

hapmap <- cbind(rs, alleles, chrom, pos, strand, assembly, center, protLSID, assayLSID, panelLSID, Qcode, GT)
write.table(hapmap, file="data/Cdiff_hapmap.txt", sep="\\t", row.names = FALSE, quote=FALSE)



################################################
### --- Filter hapmap for use with GAPIT --- ###
################################################

myG <- read.table("data/Cdiff_hapmap.txt", head = TRUE)

# Remove triallelic SNPs
myG <- myG[grep(",", myG$alleles, invert=TRUE),]

# Make chromosome names numeric and sequential for GAPIT

# Pull out chromosomes
list <- c("TR001.Ccrd1", "TR001.Ccrd2", "TR001.Ccrd3", "TR001.Ccrd4", "TR001.Ccrd5", "TR001.Ccrd6", "TR001.Ccrd7", "TR001.Ccrd8", "TR001.Ccrd9", "TR001.Ccrd10", "TR001.Ccrd11", "TR001.Ccrd12", "TR001.Ccrd13",  "TR001.Ccrd14", "TR001.Ccrd15", "TR001.Ccrd16", "TR001.Ccrd17")
myG1 <- myG[myG$chrom %in% list,]
myG1$chrom <- gsub("TR001.Ccrd", "", myG1$chrom)

# Collapse other contigs onto another chromosome
myG2 <- myG[!myG$chrom %in% list,]
myG2$chrom <- rep(18, length(nrow(myG2)))
myG2$pos <- seq(1, nrow(myG2)*10000, by=10000)

# Merge tables 
myG3 <- rbind(myG1, myG2)
write.table(myG3, file="data/Cdiff_hapmap.txt", sep="\\t", row.names = FALSE, quote=FALSE)



###################################
### --- Do GWAS using GAPIT --- ###
###################################

# load phenotypic data
myY <- read.table("data/Cdiff_pheno.txt", head = TRUE)
myY <- myY[,c(1,3:8)]

# load genotypic data
myG <- read.table("data/Cdiff_hapmap.txt", head = FALSE)
myG[1:10,1:15]
# calculate missing data
missing <- apply(myG, 1, function(x){ sum(as.numeric(x[12:ncol(myG)]=="N"), na.rm = TRUE) })
hist(missing)
# calculate heterozygosity
n.homozygous <- apply(myG, 1, function(x){ sum(as.numeric(x[12:ncol(myG)]=="A" | x[12:ncol(myG)]=="T" | x[12:ncol(myG)]=="C" | x[12:ncol(myG)]=="G"), na.rm = TRUE) })
heterozygosity <- 1 - (n.homozygous/(ncol(myG)-11-missing))
hist(heterozygosity)

# count number of individuals in each population
table(gsub("_.*", "", intersect(unlist(myG[1,12:383]), myY[,1])))

# Run GAPIT
setwd(paste0(getwd(), "/GAPIT_output_final"))
# use Major.allele.zero = TRUE to make sign of allelic effect will be relative to the minor allele
myGAPIT <- GAPIT(Y = myY, G = myG, PCA.total = 30, Model.selection = TRUE, kinship.algorithm = "VanRaden")
# according to the model selection, we should not include PCs for any trait

# test other ways methods to make sure the overall patterns do not change
setwd(paste0(getwd(), "/GAPIT_output_temp"))

# missing SNPs imputed with major allele
myGAPIT <- GAPIT(Y = myY, G = myG, PCA.total = 0, kinship.algorithm = "VanRaden", SNP.impute = "Major") # no difference

# try with no missing data
myG_nomissing <- myG[missing < 1,]
myGAPIT <- GAPIT(Y = myY, G = myG_nomissing, PCA.total = 0, kinship.algorithm = "VanRaden") # no difference

# try without CHR 18
myG_CHRonly <- myG[myG$V3!=18,]
myGAPIT <- GAPIT(Y = myY, G = myG_CHRonly, PCA.total = 0, kinship.algorithm = "VanRaden") # no difference

# force GAPIT to include some PCs
myGAPIT <- GAPIT(Y = myY, G = myG, PCA.total = 3, kinship.algorithm = "VanRaden")  # no difference
myGAPIT <- GAPIT(Y = myY, G = myG, PCA.total = 6, kinship.algorithm = "VanRaden") # just barely not significant

# try different kinship matrices
myGAPIT <- GAPIT(Y = myY, G = myG, PCA.total = 0, kinship.algorithm = "Loiselle") # no difference
myGAPIT <- GAPIT(Y = myY, G = myG, PCA.total = 0, kinship.algorithm = "EMMA") # no difference

# try the old version of GAPIT
myGAPIT <- GAPIT2(Y = myY, G = myG, PCA.total = 0, kinship.algorithm = "VanRaden")  # no difference



########################################################################################
### --- Look at the association between associated marker, phenotype, and origin --- ###
########################################################################################

# load phenotype and genotype files
phenos <- read.table("data/Cdiff_pheno.txt", head = TRUE) 
genos <- read.table("data/Cdiff_hapmap.txt", head = FALSE)

# pull out the associated locus
loci <- t(rbind(genos[1,], 
                genos[genos$V1=="TR001.Ccrd12_11392699",]))
colnames(loci) <- c("ind", "LWchr12")

# merge genotypes at the associated locus with phenotypic info
d <- merge(phenos, loci, by.x="Taxa", by.y="ind") %>% 
  mutate(Pop = str_split_fixed(Taxa, "_", 3)[,1], Mom = str_split_fixed(Taxa, "_", 3)[,2]) 
d$Origin <- factor(d$Origin, levels = c("inv", "nat"), labels = c("invasive", "native"))
d$LWchr12 <- factor(d$LWchr12, exclude = "N", levels = c("C", "Y", "T"))
d$Mom <- factor(d$Mom)

# look at plots
boxplot(d$MaxLfWdth ~ d$Origin, ylab = "Maximum leaf width (mm)", xlab = "Origin", staplewex = 0, range = 0, boxwex = 0.5)
Anova(lmer(MaxLfWdth ~ d$Origin + (1 | Mom/Pop), data = d), test.statistic = "F")
plot(lmer(MaxLfWdth ~ d$Origin + (1 | Mom/Pop), data = d))

boxplot(d$MaxLfWdth ~ d$LWchr12, ylab = "Maximum leaf width (mm)", xlab = "Genotype", staplewex = 0, range = 0, boxwex = 0.5)
Anova(lmer(MaxLfWdth ~ LWchr12 + (1 | Mom/Pop), data = d), test.statistic = "F")
plot(lmer(MaxLfWdth ~ LWchr12 + (1 | Mom/Pop), data = d))
emmeans(lmer(MaxLfWdth ~ LWchr12 + (1 | Mom/Pop), data = d), list(pairwise ~ LWchr12), adjust = "tukey")

mosaicplot(table(d$Origin, d$LWchr12), xlab = "Origin", ylab = "Genotype", main = NULL)
table(d$Origin, d$LWchr12)
fisher.test(matrix(c(547,114,81,2), ncol=2)) # p-value = 8.357e-05



##############################################
## --- Make final manhattan and qqplots --- ##
##############################################

# load data
MLW <- read.csv("GAPIT_output_final/GAPIT.MLM.MaxLfWdth.GWAS.Results.csv", header=TRUE) %>% 
  select(SNP, CHR=Chromosome, BP=Position, P=P.value) %>% filter(complete.cases(.))
p.adjust(MLW$P, method="fdr")
CD <- read.csv("GAPIT_output_final/GAPIT.MLM.CrownDiam.mm.GWAS.Results.csv", header=TRUE) %>% 
  select(SNP, CHR=Chromosome, BP=Position, P=P.value) %>% filter(complete.cases(.))
p.adjust(CD$P, method="fdr")
MLC <- read.csv("GAPIT_output_final/GAPIT.MLM.MaxLfCount.GWAS.Results.csv", header=TRUE) %>% 
  select(SNP, CHR=Chromosome, BP=Position, P=P.value) %>% filter(complete.cases(.))
p.adjust(MLC$P, method="fdr")
MLL <- read.csv("GAPIT_output_final/GAPIT.MLM.MaxLfLgth.GWAS.Results.csv", header=TRUE) %>% 
  select(SNP, CHR=Chromosome, BP=Position, P=P.value) %>% filter(complete.cases(.))
p.adjust(MLL$P, method="fdr")
SM <- read.csv("GAPIT_output_final/GAPIT.MLM.ShootMass.g.GWAS.Results.csv", header=TRUE) %>% 
  select(SNP, CHR=Chromosome, BP=Position, P=P.value) %>% filter(complete.cases(.))
p.adjust(SM$P, method="fdr")
SLA <- read.csv("GAPIT_output_final/GAPIT.MLM.SLA.GWAS.Results.csv", header=TRUE) %>% 
  select(SNP, CHR=Chromosome, BP=Position, P=P.value) %>% filter(complete.cases(.))
p.adjust(SLA$P, method="fdr")

layout(matrix(c(1,1,2), 1, 3, byrow = TRUE))
manhattan(MLW, col = c("black", "grey40"), highlight = c("TR001.Ccrd12_11392699", "TR001.Ccrd12_11392552"), main = "Maximum leaf width", suggestiveline = -log10(0.05/3508))
qq(MLW$P)

layout(matrix(c(1,1,2,3,3,4,5,5,6,7,7,8,9,9,10), 5, 3, byrow = TRUE))
par(mar = c(4,4,2,2))
manhattan(SM, col = c("black", "grey40"), main = "Shoot mass", ylim = c(0,5), suggestiveline = -log10(0.05/3508))
qq(SM$P, ylim = c(0,5))
manhattan(CD, col = c("black", "grey40"), main = "Root crown diameter", ylim = c(0,5), suggestiveline = -log10(0.05/3508))
qq(CD$P, ylim = c(0,5))
manhattan(MLC, col = c("black", "grey40"), main = "Maximum leaf count", ylim = c(0,5), suggestiveline = -log10(0.05/3508))
qq(MLC$P, ylim = c(0,5))
manhattan(MLL, col = c("black", "grey40"), main = "Maximum leaf length", ylim = c(0,5), suggestiveline = -log10(0.05/3508))
qq(MLL$P, ylim = c(0,5))
manhattan(SLA, col = c("black", "grey40"), main = "SLA", ylim = c(0,5), suggestiveline = -log10(0.05/3508))
qq(SLA$P, ylim = c(0,5))



############################################################
## --- Look at associated markers and flanking regions --- ##
############################################################

dat <- read.csv("GAPIT_output_final/GAPIT.MLM.MaxLfWdth.GWAS.Results.csv", header=TRUE) %>%
  filter(Chromosome==12) %>%
  select(BP=Position, P=P.value, MAF=maf) %>%
  mutate(logP = -log10(P))

plot(dat$BP, dat$logP, xlab="genomic position (bp)", ylab="-log(p-value)", main = "Chromosome 12 - Maximum leaf width", pch = 16)

plot(dat$BP, dat$logP, xlab="genomic position (bp)", ylab="-log(p-value)", main = "Chromosome 12 - Maximum leaf width", pch = 16, xlim = c(7e6, 1.5e7))
abline(v = 10895280, col = "blue")
abline(v = 11723974, col = "blue")

dat %>% filter(BP > 1e7) %>% filter(BP < 1.3e7) %>% arrange(BP)

# plot a 1Mb region around the SNPs 
abline(v = 10892552, col = "red")
abline(v = 11892699, col = "red")