setwd("/Users/morganslevin/Google Drive File Stream/My Drive/FAU/Dissertation/ZEFI/R")

#install packages
install.packages("data.table") #need to run in R, not R.Studio
install.packages(c("tidyverse",
                   "FSA",
                   "dplyr",
                   "reshape",
                   "factoextra",
                   "arsenal",
                   "ggplot2",
                   "Rmisc",
                   "devtools",
                   "picante")) #these are fine to run in RStudio
devtools::install_github("bryandmartin/corncob")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
BiocManager::install("DESeq2")


#Load libraries#
library(metagenomeSeq)
library(phyloseq)
library(vegan)
library(ggplot2)
library(plyr)
library(tidyverse)
library(FSA)
library(dplyr)
library(reshape)
library(factoextra)
library(cluster)
library(arsenal)
library(Rmisc)
library(corncob)
library(magrittr)
library(picante)
library(decontam)
library(plyr)
library(DESeq2)

#order of operations
#1 load data and subset to exclude pre-study samples 
#2 decontam
#4 filter
#5 normalize
#6 analyses


#File Upload#
asv_table = read.delim(file.choose(), row.names=1, header = T) #file = "feature-table.tsv    #open file in Excel/TextEdit and remove 1st row in file before import
taxa_table = read.delim(file.choose(), row.names = 1) #file= "taxonomy -> taxonomy2.tsv"   #need to delete Confidence column and separate Taxonomy column using "Text to Columns" (Data tab) in Excel by "semicolon" and populate column headers with taxonomic ranks
taxa_table = as.matrix(taxa_table)
DottedMETA = read.delim(file.choose(), row.names=1) #mapping file
dotted_meta = sample_data(DottedMETA)
ASV = otu_table(asv_table, taxa_are_rows = TRUE)
TAX = tax_table(taxa_table)
TREE =  read_tree(file.choose()) #file = tree -> tree.nwk

data <- merge_phyloseq(ASV, TAX, dotted_meta, TREE)

#Change taxa names if needed#
colnames(tax_table(data)) #check the current column names
colnames(tax_table(data))=c("Kingdom","Phylum","Class","Order","Family","Genus", "Species") #rename them
colnames(tax_table(data)) #check the renaming worked

#Check number of samples, separate into studies, recheck samples#
nsamples(data) #Should be 98#
ntaxa(data)#1701

#remove "single" practice cloacal swab samples that are unrelated to the study
studydata <- subset_samples(data, Timepoint!="pre-study") 

#Remove contaminants with the combined method of prevalence and frequency (based on blank negative controls & associated DNA quantifications ["conc"])
sample_data(studydata)$is.neg <- sample_data(studydata)$Sex == "blank"
contamdf.prev <- isContaminant(studydata, conc = studydata@sam_data$DNA.Quant, method = "combined", neg = "is.neg")
table(contamdf.prev$contaminant) #removes 20 "contaminants"

#Creating the decontaminated phyloseq object (decontam_data) #
decontam_data <- prune_taxa(!contamdf.prev$contaminant,studydata) #removes all ASVs that appear in blanks
decontam_data #inspect object; use for downstream subsetting and analysis

#FILTER decontaminated data
ntaxa(decontam_data)#check number of taxa: 1681 (1701 taxa - 20 contaminant taxa)
Fdecontam_data <- filter_taxa(decontam_data, function(x) sum(x) >10, TRUE)
ntaxa(Fdecontam_data) #315#

#remove blank negative controls from filtered and decontamninated data
Fnoblanksdata <- subset_samples(Fdecontam_data, Sex!="blank") 


#mean read depth before decontam/filtering
mean(sample_sums(studydata)); sd(sample_sums(studydata))/sqrt(length(sample_sums(studydata))) # 18758.49 ± 1234.903 reads

#mean read depth after decontam/filtering
mean(sample_sums(Fnoblanksdata)); sd(sample_sums(Fnoblanksdata))/sqrt(length(sample_sums(Fnoblanksdata))) # 17710.97 ± 1331.352 reads



#### ZEROED VARIANCE STABILIZING TRANSFORMATION WITH DESEQ2 ###
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
dds = phyloseq_to_deseq2(Fnoblanksdata, ~1)

#Calculate geometric means prior to estimate size factors#
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

#Conduct DESEQ2 test#
dds = DESeq(dds, fitType="local")

DESeqeddata = Fnoblanksdata
DESeqeddata

otu_table(DESeqeddata) <- otu_table(getVarianceStabilizedData(dds), taxa_are_rows = TRUE)

DESeqeddata

#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
#https://github.com/joey711/phyloseq/issues/445
#Make any negative values equal to 0 since the more negatie they are the more likely they were to be zero over very small and unlikely to affect your results

zeroedDESeqeddata <- DESeqeddata
zeroedDESeqeddata

DESeq2_otu_table <- as.data.frame(otu_table(zeroedDESeqeddata))

DESeq2_otu_table[DESeq2_otu_table < 0] <- 0

otu_table(zeroedDESeqeddata) <-otu_table(DESeq2_otu_table, taxa_are_rows = TRUE)

zeroedDESeqeddata
trialsdata <- zeroedDESeqeddata #rename because existing downstream script uses "trialsdata"

#Check it properly normalized
shapiro.test(sample_sums(zeroedDESeqeddata)) # p>0.05 --> normal
hist(sample_sums(zeroedDESeqeddata), xlab = "Sequences per Sample", main = "Zeroed Variance Stabilizing Transformation")

#Prove to yourself that you have the same percentage of positive numbers throughout the transformation thus making things that are negative after transformation zero is fine. 
z <- otu_table(Fnoblanksdata)
table(as.vector(z) > 0) / prod(dim(z))

#FALSE       TRUE 
#0.92614638 0.07385362  

z <- otu_table(DESeqeddata)
table(as.vector(z) > 0) / prod(dim(z))

#FALSE       TRUE 
#0.92614638 0.07385362 

z <- otu_table(zeroedDESeqeddata)
table(as.vector(z) > 0) / prod(dim(z))

#FALSE       TRUE 
#0.92614638 0.07385362  



### COMPARE ZVST TO CSS
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

#Install and load
#BiocManager::install("metagenomeSeq")
library(metagenomeSeq)

#Convert file types to use in metagenomeSeq
MGS <- phyloseq_to_metagenomeSeq(Fnoblanksdata) 

#Perform normalization following: https://bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
CSSp <- cumNormStatFast(MGS)
CSSp
MGS <- cumNorm(MGS, p = CSSp)

#Store post-transformation abundance table
normmybiom <- MRcounts(MGS, norm = T)

#Create copy of filtered data to be modified
MGSeddata = Fnoblanksdata
MGSeddata

#Switch out old ASV table for transformed one
otu_table(MGSeddata) <-otu_table(normmybiom, taxa_are_rows = TRUE)
MGSeddata

### CSS DID NOT PROPERLY NORMALIZE
shapiro.test(sample_sums(MGSeddata)) # p<0.05 --> non-normal
hist(sample_sums(MGSeddata), xlab = "Sequences per Sample", main = "Cumulative Sum Scaling")

#print comparative panel for supplement 2 figure s1
tiff("CSS histogram.tif", width = 1000, height  = 1000, res = 300)
hist(sample_sums(MGSeddata), xlab = "Sequences per Sample", main = "Cumulative Sum Scaling")
dev.off()

tiff("ZVST histogram.tif", width = 1000, height  = 1000, res = 300)
hist(sample_sums(zeroedDESeqeddata), xlab = "Sequences per Sample", main = "Zeroed Variance
     Stabilizing Transformation")
dev.off()







#subset samples based on a metadata category

#noblanks_trialsdata <- subset_samples(trialsdata, Timepoint != "blank") #remove blanks for most analyses
maledata <- subset_samples(trialsdata, Sex=="Male") #just males
femaledata <- subset_samples(trialsdata, Sex=="Female") #just females
pretrialdata <- subset_samples(trialsdata, Timepoint=="pre") #only pre-trial samples
malepretrialsdata <- subset_samples(pretrialdata, Sex == "Male")
femalepretrialsdata <- subset_samples(pretrialdata, Sex == "Female")
posttrialdata <- subset_samples(trialsdata, Timepoint=="post") #only post-trial samples
maleposttrialsdata <- subset_samples(posttrialdata, Sex == "Male")
femaleposttrialsdata <- subset_samples(posttrialdata, Sex == "Female")


#from non-CSS-normalized data for differential abundances
maledata_nodeseq <- subset_samples(Fnoblanksdata, Sex == "Male")
femaledata_nodeseq <- subset_samples(Fnoblanksdata, Sex == "Female")

malepredata_nodeseq <- subset_samples(maledata_nodeseq, Timepoint == "pre")
femalepredata_nodeseq <- subset_samples(femaledata_nodeseq, Timepoint == "pre")

malepostdata_nodeseq <- subset_samples(maledata_nodeseq, Timepoint == "post")
femalepostdata_nodeseq <- subset_samples(femaledata_nodeseq, Timepoint == "post")


###CREATE PCOAS FOR A QUICK LOOK AT THE DATA###

#Make PCOA of with Sex as color and sample Timepoint as shape for all samples#
trialsdataHornPCoAO <- ordinate(trialsdata,"PCoA", distance = "horn")
trialsdataHornPCoAP=plot_ordination(trialsdata, trialsdataHornPCoAO, color="Sex", shape = "Timepoint")+ 
  geom_point(size=4, alpha=0.4) + 
  ggtitle("Zebra Finch Cloacal Microbiome 
PCoA with Morisita-Horn Distance") +
  theme(plot.title = element_text(color = "black", hjust = 0.5)) +
  scale_colour_manual(name="Sex", values=c("black", "red", "dodgerblue"))
trialsdataHornPCoAP

tiff("PCoA.horn by sex & timepoint.tif", width = 1400, height = 1400, res = 300)
trialsdataHornPCoAP
dev.off()

#Run total PERMANOVA by sex on unmerged data#
pretrialssampledf <-data.frame(sample_data(pretrialdata))
posttrialssampledf <- data.frame(sample_data(posttrialdata))

#morisita - pre
trialsMorisita.pre <-phyloseq::distance(pretrialdata, "horn")
adonis(trialsMorisita.pre~Sex, data=pretrialssampledf, permutations=9999)

#morisita - post
trialsMorisita.post <-phyloseq::distance(posttrialdata, "horn")
adonis(trialsMorisita.post~Sex, data=posttrialssampledf, permutations=9999)


#unweighted unifrac - pre
trialsUnifrac.pre <-phyloseq::distance(pretrialdata, "unifrac")
adonis(trialsUnifrac.pre~Sex, data=pretrialssampledf, permutations=9999)

#unweighted unifrac - post
trialsUnifrac.post <-phyloseq::distance(posttrialdata, "unifrac")
adonis(trialsUnifrac.post~Sex, data=posttrialssampledf, permutations=9999)


#weighted unifrac - pre
trialsWunifrac.pre <-phyloseq::distance(pretrialdata, "wunifrac")
adonis(trialsWunifrac.pre~Sex, data=pretrialssampledf, permutations=9999)

#weighted unifrac - post
trialsWunifrac.post <-phyloseq::distance(posttrialdata, "wunifrac")
adonis(trialsWunifrac.post~Sex, data=posttrialssampledf, permutations=9999)

#Sex: multiple corrections for multiple distance metrics, order = morisita-horn, unweighted unifrac, weighted unifrac
p.adjust(c(0.0003,0.1621,0.0002), method = "BH") #pre
p.adjust(c(0.0467,0.1036,0.1442), method = "BH") #post


#Timepoint (p-values from external analysis in PRIMER-e software): multiple corrections for multiple distance metrics, order = morisita-horn, unweighted unifrac, weighted unifrac
p.adjust(c(0.34,0.63,0.38), method = "BH")

########################################
#  MDS PLOTS FOR BETA DIVERSITY BY SEX #
########################################
#how to add ellipses: https://github.com/joey711/phyloseq/issues/323

dist_methods <-unlist(distanceMethodList)
print(dist_methods)

dist_methods <-dist_methods[-c(3:13,15:47)]
print(dist_methods)

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(trialsdata, method=i)
  # Calculate ordination
  iMDS  <- ordinate(trialsdata, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(trialsdata, iMDS, color="Sex", shape = "Timepoint")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color = Sex, shape = Timepoint))
p = p + geom_point(size=2, alpha=1) + scale_color_manual(values = c("black", "red"))
p = p + facet_wrap(~distance, scales="free", nrow = 3)
p = p + theme(plot.title = element_text(hjust = 0.5))
p = p + stat_ellipse(type = "norm")
p

tiff("MDS panel_sex x timepoint.tif", width = 1400, height = 3000, res = 300)
p
dev.off()



#####Estimate Richness
trialsdata.asvtab <- as.data.frame(trialsdata@otu_table)
trialsdata.tree <- trialsdata@phy_tree
df.pd <- pd(t(trialsdata.asvtab), trialsdata.tree,include.root=T) # t(ou_table) transposes the table for use in picante and the tree file comes from the first code chunck we used to read tree file (see making a phyloseq object section).
write.csv(df.pd, "PD.csv") #save as csv... 
PD <- read.csv("PD.csv") #...then re-import...
PD <- PD[order(PD$X),] #...to order by sample ID
shannon_diversity <- as.data.frame(diversity(t(otu_table(trialsdata)), "shannon"))
write.csv(shannon_diversity, "shannon_diversity.csv")
shannon_diversity <- read.csv("shannon_diversity.csv")
shannon_diversity <- shannon_diversity[order(shannon_diversity$X),]
alpha_diversity <- cbind(as(PD, "data.frame"), as(shannon_diversity, "data.frame"))
alpha_diversity <- alpha_diversity[,-4]
names(alpha_diversity)[1:4] <- c("SampleID", "PD", "Observed.ASVs", "Shannon")

#merge alpha_diversity.csv and cognition trial/etc metadata
assay_metadata <- read.csv("cognition assay metadata.csv")
alpha_and_trial_data <- cbind(alpha_diversity, assay_metadata[,2:27])
#View(alpha_and_trial_data)

#subset by timepoint
alpha_and_trial_data_PRE <- subset(alpha_and_trial_data, alpha_and_trial_data$Timepoint=="pre")
alpha_and_trial_data_POST <- subset(alpha_and_trial_data, alpha_and_trial_data$Timepoint=="post")


#pairwise tests of difference between pre- post-assay
trialsdata_alpha_diversity_pairwise <- alpha_and_trial_data[-c(17,36,49,68),] #remove pre- samples with no post- samples

#Shannon
shannon.ttest <- t.test(trialsdata_alpha_diversity_pairwise$Shannon
                        [trialsdata_alpha_diversity_pairwise$Timepoint=="pre"],
                        trialsdata_alpha_diversity_pairwise$Shannon
                        [trialsdata_alpha_diversity_pairwise$Timepoint=="post"],
                        paired = T)$p.value # p = 0.8659

#Observed
observed.ttest <- t.test(trialsdata_alpha_diversity_pairwise$Observed.ASVs
                         [trialsdata_alpha_diversity_pairwise$Timepoint=="pre"],
                         trialsdata_alpha_diversity_pairwise$Observed.ASVs
                         [trialsdata_alpha_diversity_pairwise$Timepoint=="post"],
                         paired = T)$p.value # p = 0.2503

#PD
pd.ttest <- t.test(trialsdata_alpha_diversity_pairwise$PD
                   [trialsdata_alpha_diversity_pairwise$Timepoint=="pre"],
                   trialsdata_alpha_diversity_pairwise$PD
                   [trialsdata_alpha_diversity_pairwise$Timepoint=="post"],
                   paired = T)$p.value # p = 0.9941

p.adjust(c(shannon.ttest, observed.ttest, pd.ttest), method = "BH")
#shannon = 0.38, observed = 0.38, PD = 0.994
#***********all ns --> no difference between pre- and post-trial alpha diversities



#SUMMARY PLOTS

#Make a jitterplot for shannon diversity by timepoint#
shannon_by_timepoint <- ggplot() + 
  geom_jitter(data = trialsdata_alpha_diversity_pairwise, aes(x=Timepoint, y=Shannon, colour = Sex), width = 0.2) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Shannon Diversity") +
  xlab("") + 
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        legend.title = element_text(hjust=0.5)) +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Male"]), 
                 mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Male"]), 
                    mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Male"])),
           color = "red", size = 1.25) +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Female"]), 
                 mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Female"])),
           yend = c(mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Female"]), 
                    mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Female"])),
           color = "black", size = 1.25)
tiff("shannon_by_timepoint.tif", width = 1000, height = 2000, res = 300)
shannon_by_timepoint
dev.off()

#Make a jitterplot for Observed ASVs by timepoint#
Observed_by_timepoint <- ggplot() + 
  geom_jitter(data = alpha_and_trial_data, aes(x=Timepoint, y=Observed.ASVs, colour = Sex), width = 0.2) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Observed ASVs") +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = "none") +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Male"]), 
                 mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Male"]), 
                    mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Male"])),
           color = "red", size = 1.25) +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Female"]), 
                 mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Female"])),
           yend = c(mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Female"]), 
                    mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Female"])),
           color = "black", size = 1.25)
tiff("Observed_by_timepoint.tif", width = 1000, height = 1000, res = 300)
Observed_by_timepoint
dev.off()

#Make a jitterplot for Faith's Phylogenetic diversity by timepoint#
PD_by_timepoint <- ggplot() + 
  geom_jitter(data = alpha_and_trial_data, aes(x=Timepoint, y=PD, colour = Sex), width = 0.2) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Faith's Phylogenetic Diversity") +
  xlab("") +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = "none") +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Male"]), 
                 mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Male"]), 
                    mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Male"])),
           color = "red", size = 1.25) +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Female"]), 
                 mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Female"])),
           yend = c(mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Female"]), 
                    mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Female"])),
           color = "black", size = 1.25)
tiff("PD_by_timepoint.tif", width = 1000, height = 1000, res = 300)
PD_by_timepoint
dev.off()

tiff("alphas by sex and timepoint.tif", width = 1000, height = 2500, res = 300)
require(gridExtra)
grid.arrange(shannon_by_timepoint, Observed_by_timepoint, PD_by_timepoint, ncol = 1)
dev.off()

#Test the statistical difference between Sexes#
p2adjust <- c(kruskal.test(Shannon ~ Sex, data = alpha_and_trial_data_PRE)$p.value, 
              kruskal.test(Observed.ASVs ~ Sex, data = alpha_and_trial_data_PRE)$p.value ,
              kruskal.test(PD ~ Sex, data = alpha_and_trial_data_PRE)$p.value,
              kruskal.test(Shannon ~ Sex, data = alpha_and_trial_data_POST)$p.value ,
              kruskal.test(Observed.ASVs ~ Sex, data = alpha_and_trial_data_POST)$p.value,
              kruskal.test(PD ~ Sex, data = alpha_and_trial_data_POST)$p.value)

p.adjust(p2adjust, method = "BH")    #********* alpha diversity does not differ between sexes
#0.87 0.82 0.73 0.73 0.73 0.73

#summary for Table S1
##males
#male shannon pre
mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Male"]); sd(
  alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Male"])/
  sqrt(length((alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Male"])))

#male shannon post
mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Male"]); sd(
  alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Male"])/
  sqrt(length((alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Male"])))

#male observed pre
mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Male"]); sd(
  alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Male"])/
  sqrt(length((alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Male"])))

#male observed post
mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Male"]); sd(
  alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Male"])/
  sqrt(length((alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Male"])))

#male faith's pre
mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Male"]); sd(
  alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Male"])/
  sqrt(length((alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Male"])))

#male faith's post
mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Male"]); sd(
  alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Male"])/
  sqrt(length((alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Male"])))

##females
#female shannon pre
mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Female"]); sd(
  alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Female"])/
  sqrt(length((alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Female"])))

#female shannon post
mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Female"]); sd(
  alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Female"])/
  sqrt(length((alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Female"])))

#female observed pre
mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Female"]); sd(
  alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Female"])/
  sqrt(length((alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Female"])))

#female observed post
mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Female"]); sd(
  alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Female"])/
  sqrt(length((alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Female"])))

#female faith's pre
mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Female"]); sd(
  alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Female"])/
  sqrt(length((alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Female"])))

#female faith's post
mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Female"]); sd(
  alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Female"])/
  sqrt(length((alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Female"])))

#PLOTS
#Make a jitterplot for shannon diversity by sex#
shannon_by_sex <- ggplot() + 
  geom_jitter(data = alpha_and_trial_data, aes(x=Sex, y=Shannon, colour = Timepoint), width = 0.2) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Shannon Diversity") +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(color = "black")) + 
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Female"]), 
                 mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Female"]), 
                    mean(alpha_and_trial_data_PRE$Shannon[alpha_and_trial_data_PRE$Sex=="Male"])),
           color = "red", size = 1.25) +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Female"]), 
                 mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Female"]), 
                    mean(alpha_and_trial_data_POST$Shannon[alpha_and_trial_data_POST$Sex=="Male"])),
           color = "black", size = 1.25)
tiff("shannon_by_sex.tif", width = 1000, height = 1000, res = 300)
shannon_by_sex
dev.off()

#Make a jitterplot for observed ASVs by sex#
Observed_by_sex <- ggplot() + 
  geom_jitter(data = alpha_and_trial_data, aes(x=Sex, y=Observed.ASVs, colour = Timepoint), width = 0.2) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Observed ASVs") +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(color = "black")) + 
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Female"]), 
                 mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Female"]), 
                    mean(alpha_and_trial_data_PRE$Observed.ASVs[alpha_and_trial_data_PRE$Sex=="Male"])),
           color = "red", size = 1.25) +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Female"]), 
                 mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Female"]), 
                    mean(alpha_and_trial_data_POST$Observed.ASVs[alpha_and_trial_data_POST$Sex=="Male"])),
           color = "black", size = 1.25)
tiff("Observed_by_sex.tif", width = 1000, height = 1000, res = 300)
Observed_by_sex
dev.off()

#Make a jitterplot for Faith's phylogenetic diversity by sex#
PD_by_sex <- ggplot() + 
  geom_jitter(data = alpha_and_trial_data, aes(x=Sex, y=PD, colour = Timepoint), width = 0.2) +
  scale_color_manual(values = c("black", "red")) +
  ylab("Faith's Phylogenetic Diversity") +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.line = element_line(color = "black")) + 
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Female"]), 
                 mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Female"]), 
                    mean(alpha_and_trial_data_PRE$PD[alpha_and_trial_data_PRE$Sex=="Male"])),
           color = "red", size = 1.25) +
  annotate("segment", x = c(.7, 1.7), xend = c(1.3,2.3), 
           y = c(mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Female"]), 
                 mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Male"])),
           yend = c(mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Female"]), 
                    mean(alpha_and_trial_data_POST$PD[alpha_and_trial_data_POST$Sex=="Male"])),
           color = "black", size = 1.25)
tiff("PD_by_sex.tif", width = 1000, height = 1000, res = 300)
PD_by_sex
dev.off()





#############################################
####### TRIAL DATA vs ALPHA DIVERSITY #######
#############################################

####Build and Inspect Linear Models for cognition trial data vs 

shannon_pre <- lm(Shannon ~ Novel.Foraging + Color.Reversal + Color.Association, data = alpha_and_trial_data_PRE)
summary(shannon_pre) #ns
shannon_post <- lm(Shannon ~ Novel.Foraging + Color.Reversal + Color.Association, data = alpha_and_trial_data_POST)
summary(shannon_post) #ns


####graph of NF~Shannon relationship
shannon_v_NF <- ggplot(alpha_and_trial_data, aes(x=Novel.Foraging, y=Shannon, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Novel Foraging Performance
  (number of trials)", y = "Shannon Diversity") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  annotate("rect", xmin = 34, xmax = 52, ymin = 1.8, ymax = 2.3,
           alpha = .1) +
  annotate("text", size = 3, x = 43, y = c(2.2,2.1,2.0,1.9), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.53",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.41",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

#tiff("shannon_v_NF.tif", width = 1200, height = 1200, res = 300)
shannon_v_NF
dev.off()


## Observed ASVs
observed_pre <- lm(Observed.ASVs ~ Novel.Foraging + Color.Reversal + Color.Association, data = alpha_and_trial_data_PRE)
summary(observed_pre) #ns
observed_post <- lm(Observed.ASVs ~ Novel.Foraging + Color.Reversal + Color.Association, data = alpha_and_trial_data_POST)
summary(observed_post) #ns

####graph of NF~Observed relationship
Observed_v_NF <- ggplot(alpha_and_trial_data, aes(x=Novel.Foraging, y=Observed.ASVs, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Novel Foraging Performance
  (number of trials)", y = "Observed ASVs") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  annotate("rect", xmin= 34, xmax = 52, ymin = 57, ymax = 67,
           alpha = .1) +
  annotate("text", size = 3, x = 43, y = c(65,63,61,59), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.59",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.31",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

tiff("Observed_v_NF.tif", width = 1500, height = 1500, res = 300)
Observed_v_NF
dev.off()


#### Faith's Phylogenetic Diversity
PD_pre <- lm(PD ~ Novel.Foraging + Color.Reversal + Color.Association, data = alpha_and_trial_data_PRE)
summary(PD_pre) #ns
PD_post <- lm(PD ~ Novel.Foraging + Color.Reversal + Color.Association, data = alpha_and_trial_data_POST)
summary(PD_post) #ns


####graph of approaching-significance NF~PD relationship
PD_v_NF <- ggplot(alpha_and_trial_data, aes(x=Novel.Foraging, y=PD, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Novel Foraging Performance
  (number of trials)", y = "Faith's Phylogenetic Diversity") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") + 
  annotate("rect", xmin = 31, xmax = 55, ymin = 0.3, ymax = 1.7,
           alpha = .1) +
  annotate("text", size = 3, x = 43, y = c(1.5,1.2,0.9,0.6), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.80",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.79",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

tiff("PD_v_NF.tif", width = 1700, height = 1500, res = 300)
PD_v_NF
dev.off()



####graph of  Shannon~CA relationship
Shannon_v_CA <- ggplot(alpha_and_trial_data, aes(x=Color.Association, y=Shannon, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Color Association Performance
  (number of trials)", y = "") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") + 
  annotate("rect", xmin = 24, xmax = 36, ymin = 1.8, ymax = 2.3,
           alpha = .1) +
  annotate("text", size = 3, x = 30, y = c(2.2,2.1,2.0,1.9), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.72",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.92",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

tiff("Shannon vs CA.tif", width = 1500, height = 1500, res = 300)
Shannon_v_CA
dev.off()

####graph of  Observed~CA relationship
Observed_v_CA <- ggplot(alpha_and_trial_data, aes(x=Color.Association, y=Observed.ASVs, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Color Association Performance
  (number of trials)", y = "") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") + 
  annotate("rect", xmin = 24, xmax = 36, ymin = 49, ymax = 62,
           alpha = .1) +
  annotate("text", size = 3, x = 30, y = c(60,57,54,51), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.82",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.93",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

tiff("Observed vs CA.tif", width = 1500, height = 1500, res = 300)
Observed_v_CA
dev.off()


####graph of  PD~CA relationship
PD_v_CA <- ggplot(alpha_and_trial_data, aes(x=Color.Association, y=PD, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Color Association Performance
  (number of trials)", y = "") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") + 
  annotate("rect", xmin = 24, xmax = 36, ymin = 0.5, ymax = 1.9,
           alpha = .1) +
  annotate("text", size = 3, x = 30, y = c(1.6,1.3,1.0,0.7), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.79",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.86",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

tiff("PD vs CA.tif", width = 1700, height = 1500, res = 300)
PD_v_CA
dev.off()


####graph of  Shannon~CR relationship
Shannon_v_CR <- ggplot(alpha_and_trial_data, aes(x=Color.Reversal, y=Shannon, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Color Reversal Performance
  (number of trials)", y = "") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") + 
  annotate("rect", xmin = 24, xmax = 36, ymin = 1.8, ymax = 2.3,
           alpha = .1) +
  annotate("text", size = 3, x = 30, y = c(2.2,2.1,2.0,1.9), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.54",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.65",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

tiff("Shannon vs CR.tif", width = 1500, height = 1500, res = 300)
Shannon_v_CR
dev.off()

####graph of  Observed~CR relationship
Observed_v_CR <- ggplot(alpha_and_trial_data, aes(x=Color.Reversal, y=Observed.ASVs, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Color Reversal Performance
  (number of trials)", y = "") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") + 
  annotate("rect", xmin = 24, xmax = 36, ymin = 49, ymax = 62,
           alpha = .1) +
  annotate("text", size = 3, x = 30, y = c(60,57,54,51), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.94",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.52",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

tiff("Observed vs CR.tif", width = 1500, height = 1500, res = 300)
Observed_v_CR
dev.off()


####graph of  PD~CR relationship
PD_v_CR <- ggplot(alpha_and_trial_data, aes(x=Color.Reversal, y=PD, color = Timepoint)) + 
  scale_color_manual(values = c("#3300FF", "#FF33FF")) +
  stat_smooth(method = lm, data = alpha_and_trial_data_POST, formula = y~x,
              fullrange = T, color = "#3300FF", fill = "#3300FF", alpha = 0.1) +
  stat_smooth(method = lm, data = alpha_and_trial_data_PRE, formula = y~x,
              fullrange = T, color = "#FF33FF", fill = "#FF33FF", alpha = 0.1) +
  geom_point(size = 1) +
  labs(x = "Color Reversal Performance
  (number of trials)", y = "") +
  theme_classic() +
  theme(panel.background = element_rect(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=13),
        plot.background = element_rect(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") + 
  annotate("rect", xmin = 24, xmax = 36, ymin = 0.5, ymax = 1.9,
           alpha = .1) +
  annotate("text", size = 3, x = 30, y = c(1.6,1.3,1.0,0.7), 
           colour = c("#3300FF", "#3300FF", "#FF33FF", "#FF33FF"),
           label = c("italic(P)(post) == 0.90",
                     "italic(n) == 34",
                     "italic(P)(pre) == 0.52",
                     "italic(n) == 38"), 
           parse = T, color = "black") 

tiff("PD vs CR.tif", width = 1500, height = 1500, res = 300)
PD_v_CR
dev.off()

tiff("alphas v performance.tif", width = 3000, height = 3000, res = 300)
require(gridExtra)
grid.arrange(shannon_v_NF, Shannon_v_CA, Shannon_v_CR, 
             Observed_v_NF, Observed_v_CA, Observed_v_CR,
             PD_v_NF, PD_v_CA, PD_v_CR, ncol = 3, nrow = 3)
dev.off() 


##########################
#    ALPHA DIVERSITY 
# MULTIPLE TEST CORRECTION
##########################

# order below = "c(Shannon, Observed, Faith's)"
p.adjust(c(0.407,0.313908,0.789), method = "BH") #pre NF
p.adjust(c(0.534,0.589369,0.798), method = "BH") #post NF
p.adjust(c(0.924,0.930027,0.856), method = "BH") #pre CA
p.adjust(c(0.715,0.817045,0.782), method = "BH") #post CA
p.adjust(c(0.645,0.523793,0.518), method = "BH") #pre CR
p.adjust(c(0.540,0.941933,0.893), method = "BH") #post CR




#############################################
####### TRIAL DATA vs BETA DIVERSITY #######
#############################################

malepretrialssampledf <- data.frame(sample_data(malepretrialsdata))
femalepretrialssampledf <- data.frame(sample_data(femalepretrialsdata))
maleposttrialssampledf <- data.frame(sample_data(maleposttrialsdata))
femaleposttrialssampledf <- data.frame(sample_data(femaleposttrialsdata))

###################
# HORN (MORISITA) #
###################

#Run total PERMANOVA on unmerged data based on variables#
pretrialsMorisita <-phyloseq::distance(pretrialdata, "horn")
posttrialsMorisita <-phyloseq::distance(posttrialdata, "horn")

#repeat for male-only pre-trial data
malepretrialsMorisita <-phyloseq::distance(malepretrialsdata, "horn")

#repeat for male-only post-trial data
maleposttrialsMorisita <-phyloseq::distance(maleposttrialsdata, "horn")

#repeat for female-only pre-trial data
femalepretrialsMorisita <-phyloseq::distance(femalepretrialsdata, "horn")

#repeat for female-only post-trial data
femaleposttrialsMorisita <-phyloseq::distance(femaleposttrialsdata, "horn")


#PERMANOVA by Novel Foraging - PRE
#male
adonis(malepretrialsMorisita ~ Novel.Foraging, data = malepretrialssampledf, permutations=9999) #significant
#female
adonis(femalepretrialsMorisita ~ Novel.Foraging, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Novel Foraging - POST
#male
adonis(maleposttrialsMorisita ~ Novel.Foraging, data = maleposttrialssampledf, permutations=9999) #significant
#female
adonis(femaleposttrialsMorisita ~ Novel.Foraging, data = femaleposttrialssampledf, permutations=9999) #ns


#PERMANOVA by Color Ass. - PRE
#male
adonis(malepretrialsMorisita ~ Color.Association, data = malepretrialssampledf, permutations=9999) #ns
#female
adonis(femalepretrialsMorisita ~ Color.Association, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Ass. - POST
#male
adonis(maleposttrialsMorisita ~ Color.Association, data = maleposttrialssampledf, permutations=9999) #ns
#female
adonis(femaleposttrialsMorisita ~ Color.Association, data = femaleposttrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Rev. - PRE
#male
adonis(malepretrialsMorisita ~ Color.Reversal, data = malepretrialssampledf, permutations=9999) #ns
#female
adonis(femalepretrialsMorisita ~ Color.Reversal, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Rev. - POST
#male
adonis(maleposttrialsMorisita ~ Color.Reversal, data = maleposttrialssampledf, permutations=9999) #ns
#female
adonis(femaleposttrialsMorisita ~ Color.Reversal, data = femaleposttrialssampledf, permutations=9999) #ns




#################
#      UNIFRAC       #
#################

#Run total PERMANOVA on unmerged data based on variables#
pretrialsUnifrac <-phyloseq::distance(pretrialdata, "unifrac")
posttrialsUnifrac <-phyloseq::distance(posttrialdata, "unifrac")

#repeat for male-only pre-trial data
malepretrialsUnifrac <-phyloseq::distance(malepretrialsdata, "unifrac")

#repeat for male-only post-trial data
maleposttrialsUnifrac <-phyloseq::distance(maleposttrialsdata, "unifrac")

#repeat for female-only pre-trial data
femalepretrialsUnifrac <-phyloseq::distance(femalepretrialsdata, "unifrac")

#repeat for female-only post-trial data
femaleposttrialsUnifrac <-phyloseq::distance(femaleposttrialsdata, "unifrac")

#PERMANOVA by Novel Foraging - PRE
#male
adonis(malepretrialsUnifrac ~ Novel.Foraging, data = malepretrialssampledf, permutations=9999) #significant
#female
adonis(femalepretrialsUnifrac ~ Novel.Foraging, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Novel Foraging - POST
#male
adonis(maleposttrialsUnifrac ~ Novel.Foraging, data = maleposttrialssampledf, permutations=9999) #significant
#female
adonis(femaleposttrialsUnifrac ~ Novel.Foraging, data = femaleposttrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Ass. - PRE
#male
adonis(malepretrialsUnifrac ~ Color.Association, data = malepretrialssampledf, permutations=9999) #ns
#female
adonis(femalepretrialsUnifrac ~ Color.Association, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Ass. - POST
#male
adonis(maleposttrialsUnifrac ~ Color.Association, data = maleposttrialssampledf, permutations=9999) #ns
#female
adonis(femaleposttrialsUnifrac ~ Color.Association, data = femaleposttrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Rev. - PRE
#male
adonis(malepretrialsUnifrac ~ Color.Reversal, data = malepretrialssampledf, permutations=9999) #ns
#female
adonis(femalepretrialsUnifrac ~ Color.Reversal, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Rev. - POST
#male
adonis(maleposttrialsUnifrac ~ Color.Reversal, data = maleposttrialssampledf, permutations=9999) #ns
#female
adonis(femaleposttrialsUnifrac ~ Color.Reversal, data = femaleposttrialssampledf, permutations=9999) #ns



#################
#  WEIGHTED  UNIFRAC #
#################

#Run total PERMANOVA on unmerged data based on variables#
pretrialsWunifrac <-phyloseq::distance(pretrialdata, "wunifrac")
posttrialsWunifrac <-phyloseq::distance(posttrialdata, "wunifrac")

#repeat for male-only pre-trial data
malepretrialsWunifrac <-phyloseq::distance(malepretrialsdata, "wunifrac")

#repeat for male-only post-trial data
maleposttrialsWunifrac <-phyloseq::distance(maleposttrialsdata, "wunifrac")

#repeat for female-only pre-trial data
femalepretrialsWunifrac <-phyloseq::distance(femalepretrialsdata, "wunifrac")

#repeat for female-only post-trial data
femaleposttrialsWunifrac <-phyloseq::distance(femaleposttrialsdata, "wunifrac")

#PERMANOVA by Novel Foraging - PRE
#male
adonis(malepretrialsWunifrac ~ Novel.Foraging, data = malepretrialssampledf, permutations=9999) #significant
#female
adonis(femalepretrialsWunifrac ~ Novel.Foraging, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Novel Foraging - POST
#male
adonis(maleposttrialsWunifrac ~ Novel.Foraging, data = maleposttrialssampledf, permutations=9999) #significant
#female
adonis(femaleposttrialsWunifrac ~ Novel.Foraging, data = femaleposttrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Ass. - PRE
#male
adonis(malepretrialsWunifrac ~ Color.Association, data = malepretrialssampledf, permutations=9999) #ns
#female
adonis(femalepretrialsWunifrac ~ Color.Association, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Ass. - POST
#male
adonis(maleposttrialsWunifrac ~ Color.Association, data = maleposttrialssampledf, permutations=9999) #ns
#female
adonis(femaleposttrialsWunifrac ~ Color.Association, data = femaleposttrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Rev. - PRE
#male
adonis(malepretrialsWunifrac ~ Color.Reversal, data = malepretrialssampledf, permutations=9999) #ns
#female
adonis(femalepretrialsWunifrac ~ Color.Reversal, data = femalepretrialssampledf, permutations=9999) #ns

#PERMANOVA by Color Rev. - POST
#male
adonis(maleposttrialsWunifrac ~ Color.Reversal, data = maleposttrialssampledf, permutations=9999) #ns
#female
adonis(femaleposttrialsWunifrac ~ Color.Reversal, data = femaleposttrialssampledf, permutations=9999) #ns


##########################
#      BETA DIVERSITY 
# MULTIPLE TEST CORRECTION
##########################

# order below = "c(Morisita-Horn, UniFrac, Weighted UniFrac)"
p.adjust(c(0.5234, 0.5577, 0.57), method = "BH") #male pre NF
p.adjust(c(0.4741, 0.4328, 0.0941), method = "BH") #male post NF
p.adjust(c(0.1842, 0.7006, 0.0989), method = "BH") #female pre NF
p.adjust(c(0.0715, 0.2794, 0.0489), method = "BH") #female post NF
p.adjust(c(0.4139, 0.1047, 0.0696), method = "BH") #male pre CA
p.adjust(c(0.9925, 0.9149, 0.9828), method = "BH") #male post CA
p.adjust(c(0.1516, 0.8252, 0.9618), method = "BH") #female pre CA
p.adjust(c(0.4196, 0.1882, 0.3343), method = "BH") #female post CA
p.adjust(c(0.3672, 0.8567, 0.4887), method = "BH") #male pre CR
p.adjust(c(0.2581, 0.6567, 0.8175), method = "BH") #male post CR
p.adjust(c(0.5128, 0.1543, 0.3824), method = "BH") #female pre CR
p.adjust(c(0.5841, 0.6196, 0.6457), method = "BH") #female post CR

###############
#  MDS PLOTS  
###############


#MALE AND FEMALE POST NOVEL FORAGING WEIGHTED UNIFRAC
dist_methods <-unlist(distanceMethodList)
print(dist_methods)

dist_methods <-dist_methods[2]
print(dist_methods)

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(posttrialdata, method=i)
  # Calculate ordination
  iMDS  <- ordinate(posttrialdata, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(posttrialdata, iMDS, color="Novel.Foraging", shape = "Sex")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
MDS_NF_mf_post_wunifrac = ggplot(df, aes(Axis.1, Axis.2, color = Novel.Foraging, shape = Sex)) +
  geom_point(size=2) +
  facet_wrap(~distance, scales="free") +
  scale_shape_discrete(name = "Sex") + 
  scale_color_manual(name = "Cognitive\\nPerformance", values = c("orange", "#9933FF", "#00CC66")) +
  ggtitle("Novel Foraging (both sexes post-trial)")+
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_blank(),
        legend.position = "none")+
  stat_ellipse(type = "t", data = subset(df, Sex == "Male"), linetype = 1) + 
  stat_ellipse(type = "t", data = subset(df, Sex == "Female"), linetype = 2)
MDS_NF_mf_post_wunifrac

tiff("MDS both post NF_wunifrac.tif", width = 1000, height = 1000, res = 300)
MDS_NF_mf_post_wunifrac
dev.off()


#FEMALE POST NOVEL FORAGING HORN
dist_methods <-unlist(distanceMethodList)
print(dist_methods)

dist_methods <-dist_methods[14]
print(dist_methods)

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(femaleposttrialsdata, method=i)
  # Calculate ordination
  iMDS  <- ordinate(femaleposttrialsdata, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(femaleposttrialsdata, iMDS, color="Novel.Foraging")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
MDS_NF_f_post_horn = ggplot(df, aes(Axis.1, Axis.2, color = Novel.Foraging)) +
  geom_point(size=2) +
  scale_color_manual(values = c("orange", "#9933FF", "#00CC66")) +
  facet_wrap(~distance, scales="free") +
  ggtitle("Novel Foraging (female post-trial)")+
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.title = element_blank()) +
  stat_ellipse(type = "t", data = subset(df, Sex == "Female"), linetype = 2)
MDS_NF_f_post_horn

tiff("MDS NF_f_post_horn.tif", width = 1000, height = 1000, res = 300)
MDS_NF_f_post_horn
dev.off()


#FEMALE PRE NOVEL FORAGING WEIGHTED UNIFRAC
dist_methods <-unlist(distanceMethodList)
print(dist_methods)

dist_methods <-dist_methods[2]
print(dist_methods)

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(femalepretrialsdata, method=i)
  # Calculate ordination
  iMDS  <- ordinate(femalepretrialsdata, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(femalepretrialsdata, iMDS, color="Novel.Foraging")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
MDS_NF_f_post_wunifrac = ggplot(df, aes(Axis.1, Axis.2, color = Novel.Foraging)) +
  geom_point(size=2) +
  scale_color_manual(values = c("orange", "#9933FF", "#00CC66")) +
  facet_wrap(~distance, scales="free") +
  ggtitle("Novel Foraging (female pre-trial)")+
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.title = element_blank())+
  stat_ellipse(type = "t", data = subset(df, Sex == "Female"), linetype = 2)
MDS_NF_f_post_wunifrac

tiff("MDS NF_f_pre_wunifrac.tif", width = 1000, height = 1000, res = 300)
MDS_NF_f_post_wunifrac
dev.off()


#MALE PRE COLOR ASSOCIATION WEIGHTED UNIFRAC
dist_methods <-unlist(distanceMethodList)
print(dist_methods)

dist_methods <-dist_methods[2]
print(dist_methods)

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(malepretrialsdata, method=i)
  # Calculate ordination
  iMDS  <- ordinate(malepretrialsdata, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(malepretrialsdata, iMDS, color="Color.Association")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
MDS_CA_m_pre_wunifrac = ggplot(df, aes(Axis.1, Axis.2, color = Color.Association)) +
  geom_point(size=2) +
  scale_color_manual(values = c("orange", "#9933FF", "#00CC66")) +  facet_wrap(~distance, scales="free") +
  ggtitle("Color Association (male pre-trial)")+
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.title = element_blank()) +
  stat_ellipse(type = "t", data = subset(df, Sex == "Male"))
MDS_CA_m_pre_wunifrac

tiff("MDS CA_m_pre_wunifrac.tif", width = 1000, height = 1000, res = 300)
MDS_CA_m_pre_wunifrac
dev.off()



################## macro-level summary graphs of all beta diversity metrics ##################

####NOVEL FORAGING
dist_methods <-unlist(distanceMethodList)
print(dist_methods)

dist_methods <-dist_methods[-c(3:13,15:47)]
print(dist_methods)

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(trialsdata, method=i)
  # Calculate ordination
  iMDS  <- ordinate(trialsdata, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(trialsdata, iMDS, color="Novel.Foraging", shape = "Sex")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
MDS_NF = ggplot(df, aes(Axis.1, Axis.2, color = Novel.Foraging, shape = Sex)) +
  geom_point(size=1.5) +
  scale_shape_discrete(name = "Sex") + 
  scale_color_manual(name = "Cognitive\\nPerformance", values = c("orange", "#9933FF", "#00CC66")) +
  facet_wrap(~distance, scales="free") +
  ggtitle("Novel Foraging")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.title = element_blank()) + 
  stat_ellipse(type = "t", data = subset(df, Sex == "Male"), linetype = 1) + 
  stat_ellipse(type = "t", data = subset(df, Sex == "Female"), linetype = 2)
MDS_NF

tiff("MDS NF.tif", width = 2100, height = 700, res = 300)
MDS_NF
dev.off()

####COLOR ASSOCIATION
dist_methods <-unlist(distanceMethodList)
print(dist_methods)

dist_methods <-dist_methods[-c(3:13,15:47)]
print(dist_methods)

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(trialsdata, method=i)
  # Calculate ordination
  iMDS  <- ordinate(trialsdata, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(trialsdata, iMDS, color="Color.Association", shape = "Sex")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
MDS_CA = ggplot(df, aes(Axis.1, Axis.2, color = Color.Association, shape = Sex)) +
  geom_point(size=1.5) +
  scale_color_manual(name = "Cognitive\\nPerformance", values = c("orange", "#9933FF", "#00CC66")) +
  facet_wrap(~distance, scales="free") +
  ggtitle("Color Association")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.title = element_blank()) +
  stat_ellipse(type = "t", data = subset(df, Sex == "Male"), linetype = 1) + 
  stat_ellipse(type = "t", data = subset(df, Sex == "Female"), linetype = 2)
MDS_CA

tiff("MDS CA.tif", width = 2100, height = 700, res = 300)
MDS_CA
dev.off()

####COLOR REVERSAL
dist_methods <-unlist(distanceMethodList)
print(dist_methods)

dist_methods <-dist_methods[-c(3:13,15:47)]
print(dist_methods)

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(trialsdata, method=i)
  # Calculate ordination
  iMDS  <- ordinate(trialsdata, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(trialsdata, iMDS, color="Color.Reversal", shape = "Sex")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
MDS_CR = ggplot(df, aes(Axis.1, Axis.2, color = Color.Reversal, shape = Sex)) +
  geom_point(size=1.5) +
  scale_color_manual(name = "Cognitive\\nPerformance", values = c("orange", "#9933FF", "#00CC66")) +
  facet_wrap(~distance, scales="free") +
  ggtitle("Color Reversal")+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.title = element_blank())+
  stat_ellipse(type = "t", data = subset(df, Sex == "Male"), linetype = 1) + 
  stat_ellipse(type = "t", data = subset(df, Sex == "Female"), linetype = 2)
MDS_CR

tiff("MDS CR.tif", width = 2100, height = 700, res = 300)
MDS_CR
dev.off()



####################   CORNCOB   #########################

#### install/load lines below if not done already at top of this script

#install.packages("devtools")
#devtools::install_github("bryandmartin/corncob")
#library(phyloseq)
#library(corncob)
#library(magrittr)


### Need to do separately for males and females b/c beta diversity differed by sex

### Male Pre CA
# inspect data array
malepredata_nodeseq
Fmalepredata_nodeseq <- filter_taxa(malepredata_nodeseq, function(x) sum(x) >0, TRUE)
malepredata_nodeseq.clean <- clean_taxa_names(Fmalepredata_nodeseq)

zefi_genus.m_pre <- malepredata_nodeseq.clean %>% 
  tax_glom("Genus")


### ANALYSIS OF MULTIPLE TAXA ###

## at genus level for color association##
# We specify the covariates of our model using formula and phi.formula as before, except we no longer include the response term because we are testing multiple taxa. We also specify which covariates we want to test for by removing them in the formula\\_null and phi.formula\\_null arguments. #The difference between the formulas and the null version of the formulas are the variables that we test.
set.seed(1)
male_genus_CA_da_analysis <- differentialTest(formula = ~ Color.Association,
                                              phi.formula = ~ Color.Association,
                                              formula_null = ~ 1,
                                              phi.formula_null = ~ Color.Association,
                                              test = "Wald", boot = FALSE,
                                              data = zefi_genus.m_pre,
                                              fdr_cutoff = 0.05)

# see a list of differentially-abundant taxa using: 
male_genus_CA_da_analysis$significant_taxa #OTUs 21 (Pseudomonas, P=0.0001), 54 (Pasteurellaceae P=0.005), 56 (Gallibacterium P<0.001), 60 (Stenotrophomonas P=0.01), 71 (Helicobacter P=0.03) 122 (Enterococcus P=0.02)

# now repeat for differential variance:
set.seed(1)
male_genus_CA_dv_analysis <- differentialTest(formula = ~ Color.Association,
                                              phi.formula = ~ Color.Association,
                                              formula_null = ~ Color.Association,
                                              phi.formula_null = ~ 1,
                                              data = zefi_genus.m_pre,
                                              test = "LRT", boot = FALSE,
                                              fdr_cutoff = 0.05)
# see a list of differentially-variable taxa using:
male_genus_CA_dv_analysis$significant_taxa #"OTU21"  "OTU56"  "OTU141"

# switch the OTU labels to taxonomic labels using otu_to_taxonomy. We supply our OTU labels as strings for the OTU argument. We supply the {phyloseq} object for the data argument.
otu_to_taxonomy(OTU = male_genus_CA_da_analysis$significant_taxa, data = zefi_genus.m_pre)
otu_to_taxonomy(OTU = male_genus_CA_dv_analysis$significant_taxa, data = zefi_genus.m_pre)

# We can examine the p-values of our tests using:
male_genus_CA_da_analysis$p

# We can examine the p-values after controlling for the false discovery rate (where the values are now adjusted to control the false discovery rate at 0.05) using:
male_genus_CA_da_analysis$p_fdr ######## NEED TO GET P VALUES FOR EACH GENUS FOR TEXT

# plot the model coefficients of our results:
corncob_male_genus_CA <- plot(male_genus_CA_da_analysis)

tiff("corncob male pre genus CA.tif", width = 4000, height = 1000, res = 300)
corncob_male_genus_CA
dev.off()


# Finally, we can see a list of any taxa for which we were not able to fit a model using:
which(is.na(male_genus_da_analysis$p)) %>% names


### Male POST NF
# inspect data array
malepostdata_nodeseq
Fmalepostdata_nodeseq <- filter_taxa(malepostdata_nodeseq, function(x) sum(x) >0, TRUE)
malepostdata_nodeseq.clean <- clean_taxa_names(Fmalepostdata_nodeseq)

zefi_genus.m_post <- malepostdata_nodeseq.clean %>% 
  tax_glom("Genus")

# We specify the covariates of our model using formula and phi.formula as before, except we no longer include the response term because we are testing multiple taxa. We also specify which covariates we want to test for by removing them in the formula\\_null and phi.formula\\_null arguments. #The difference between the formulas and the null version of the formulas are the variables that we test.
set.seed(1)
male_genus_NF_da_analysis <- differentialTest(formula = ~ Novel.Foraging,
                                              phi.formula = ~ Novel.Foraging,
                                              formula_null = ~ 1,
                                              phi.formula_null = ~ Novel.Foraging,
                                              test = "Wald", boot = FALSE,
                                              data = zefi_genus.m_post,
                                              fdr_cutoff = 0.05)

# see a list of differentially-abundant taxa using: 
male_genus_NF_da_analysis$significant_taxa #107 (Enterococcus P=0.01)

# now repeat for differential variance:
set.seed(1)
male_genus_NF_dv_analysis <- differentialTest(formula = ~ Novel.Foraging,
                                              phi.formula = ~ Novel.Foraging,
                                              formula_null = ~ Novel.Foraging,
                                              phi.formula_null = ~ 1,
                                              data = zefi_genus.m_post,
                                              test = "LRT", boot = FALSE,
                                              fdr_cutoff = 0.05)
# see a list of differentially-variable taxa using:
male_genus_NF_dv_analysis$significant_taxa # none

# switch the OTU labels to taxonomic labels using otu_to_taxonomy. We supply our OTU labels as strings for the OTU argument. We supply the {phyloseq} object for the data argument.
otu_to_taxonomy(OTU = male_genus_NF_da_analysis$significant_taxa, data = zefi_genus.m_post)
#otu_to_taxonomy(OTU = male_genus_NF_dv_analysis$significant_taxa, data = zefi_genus.m_post)

# We can examine the p-values of our tests using:
male_genus_NF_da_analysis$p

# We can examine the p-values after controlling for the false discovery rate (where the values are now adjusted to control the false discovery rate at 0.05) using:
male_genus_NF_da_analysis$p_fdr ######## NEED TO GET P VALUES FOR EACH GENUS FOR TEXT

# plot the model coefficients of our results:
corncob_male_genus_NF <- plot(male_genus_NF_da_analysis)

tiff("corncob male post genus NF.tif", width = 4000, height = 1000, res = 300)
corncob_male_genus_NF
dev.off()


# Finally, we can see a list of any taxa for which we were not able to fit a model using:
which(is.na(male_genus_NF_da_analysis$p)) %>% names



### Females

# inspect data array
femaledata_nodeseq
Ffemaledata_nodeseq <- filter_taxa(femaledata_nodeseq, function(x) sum(x) >0, TRUE)
femaledata_nodeseq.clean <- clean_taxa_names(Ffemaledata_nodeseq)

zefi_genus.f <- femaledata_nodeseq.clean %>% 
  tax_glom("Genus")

## at genus level for novel foraging##
# We specify the covariates of our model using formula and phi.formula as before, except we no longer include the response term because we are testing multiple taxa. We also specify which covariates we want to test for by removing them in the formula\\_null and phi.formula\\_null arguments. #The difference between the formulas and the null version of the formulas are the variables that we test.
set.seed(1)
female_genus_NF_da_analysis <- differentialTest(formula = ~ Novel.Foraging,
                                                phi.formula = ~ Novel.Foraging,
                                                formula_null = ~ 1,
                                                phi.formula_null = ~ Novel.Foraging,
                                                test = "Wald", boot = FALSE,
                                                data = zefi_genus.f,
                                                fdr_cutoff = 0.05)

# see a list of differentially-abundant taxa using: 
female_genus_NF_da_analysis$significant_taxa #OTU39 (Gallibacterium P<0.001), OTU115 (Catellicoccus P=0.003), OTU116 (Rothia P=0.01)

# now repeat for differential variance:
set.seed(1)
female_genus_NF_dv_analysis <- differentialTest(formula = ~ Color.Reversal,
                                                phi.formula = ~ Color.Reversal,
                                                formula_null = ~ Color.Reversal,
                                                phi.formula_null = ~ 1,
                                                data = zefi_genus.f,
                                                test = "LRT", boot = FALSE,
                                                fdr_cutoff = 0.05)
# see a list of differentially-variable taxa using:
female_genus_NF_dv_analysis$significant_taxa #OTU10, OTU33

# switch the OTU labels to taxonomic labels using otu_to_taxonomy. We supply our OTU labels as strings for the OTU argument. We supply the {phyloseq} object for the data argument.
otu_to_taxonomy(OTU = female_genus_NF_da_analysis$significant_taxa, data = zefi_genus.f)
otu_to_taxonomy(OTU = female_genus_NF_dv_analysis$significant_taxa, data = zefi_genus.f)

# We can examine the p-values of our tests using:
female_genus_NF_da_analysis$p

# We can examine the p-values after controlling for the false discovery rate (where the values are now adjusted to control the false discovery rate at 0.05) using:
female_genus_NF_da_analysis$p_fdr ######## NEED TO GET P VALUES FOR EACH GENUS FOR TEXT

# plot the model coefficients of our results:
corncob_female_genus_NF <- plot(female_genus_NF_da_analysis)

tiff("corncob female genus NF.tif", width = 3800, height = 1000, res = 300)
corncob_female_genus_NF
dev.off()


# Finally, we can see a list of any taxa for which we were not able to fit a model using:
which(is.na(male_genus_da_analysis$p)) %>% names




##############CALCULATE TAXA-SPECIFIC RELATIVE ABUNDANCES (to pair with corncob results)#################

##MALES - PRE - COLOR ASSOCIATION
##Create table ready for making stacked bar graph for Genus <1%  by Sample##
malepredataphy <- transform_sample_counts(malepredata_nodeseq, function(x) 100*x/sum(x))
# agglomerate taxa
glom_genus.m <- tax_glom(malepredataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat_genus.m <- psmelt(glom_genus.m)
# convert Genus to a character vector from a factor because R
dat_genus.m$Genus <- as.character(dat_genus.m$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means_genus.m <- ddply(dat_genus.m, ~Genus, function(x) c(mean=mean(x$Abundance)))
#make separate table with taxa-specific abundances
se_genus.m <- ddply(dat_genus.m, ~Genus, function(x) c(se=sd(x$Abundance)/sqrt(length(x$Abundance))))
abundances<- as.data.frame(cbind(means_genus.m, se_genus.m$se))
#remove unnessary columns
dat_genus.m_CA <- subset(dat_genus.m, select=c(Sex, Timepoint, Abundance, Genus, Color.Association))
#Summarize based upon target parameter
dat_genus.m_CA <- summarySE(data=dat_genus.m_CA, measurevar="Abundance", groupvars=c("Color.Association", "Genus"), na.rm = TRUE)
#Arrange by Genus
dat_genus.m_CA <- arrange(dat_genus.m_CA, Genus)

##MALES - POST - NOVEL FORAGING
##Create table ready for making stacked bar graph for Genus <1%  by Sample##
malepostdataphy <- transform_sample_counts(malepostdata_nodeseq, function(x) 100*x/sum(x))
# agglomerate taxa
glom_genus.m <- tax_glom(malepostdataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat_genus.m <- psmelt(glom_genus.m)
# convert Genus to a character vector from a factor because R
dat_genus.m$Genus <- as.character(dat_genus.m$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means_genus.m <- ddply(dat_genus.m, ~Genus, function(x) c(mean=mean(x$Abundance)))
#make separate table with taxa-specific abundances
se_genus.m <- ddply(dat_genus.m, ~Genus, function(x) c(se=sd(x$Abundance)/sqrt(length(x$Abundance))))
abundances<- as.data.frame(cbind(means_genus.m, se_genus.m$se))
#remove unnessary columns
dat_genus.m_NF <- subset(dat_genus.m, select=c(Sex, Timepoint, Abundance, Genus, Novel.Foraging))
#Summarize based upon target parameter
dat_genus.m_NF <- summarySE(data=dat_genus.m_NF, measurevar="Abundance", groupvars=c("Novel.Foraging", "Genus"), na.rm = TRUE)
#Arrange by Genus
dat_genus.m_NF <- arrange(dat_genus.m_NF, Genus)


##FEMALES - NOVEL FORAGING
##Create table ready for making stacked bar graph for Genus <1%  by Sample##
femaledataphy <- transform_sample_counts(femaledata_nodeseq, function(x) 100*x/sum(x))
# agglomerate taxa
glom_genus.f <- tax_glom(femaledataphy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat_genus.f <- psmelt(glom_genus.f)
# convert Genus to a character vector from a factor because R
dat_genus.f$Genus <- as.character(dat_genus.f$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means_genus.f <- ddply(dat_genus.f, ~Genus, function(x) c(mean=mean(x$Abundance)))
#make separate table with taxa-specific abundances
se_genus.f <- ddply(dat_genus.f, ~Genus, function(x) c(se=sd(x$Abundance)/sqrt(length(x$Abundance))))
abundances<- as.data.frame(cbind(means_genus.f, se_genus.f$se))
#remove unnessary columns
dat_genus.f <- subset(dat_genus.f, select=c(Sex, Timepoint, Abundance, Genus, Novel.Foraging))
#Summarize based upon target parameter
dat_genus.f <- summarySE(data=dat_genus.f, measurevar="Abundance", groupvars=c("Novel.Foraging", "Genus"), na.rm = TRUE)
#Arrange by Genus
dat_genus.f <- arrange(dat_genus.f, Genus)

