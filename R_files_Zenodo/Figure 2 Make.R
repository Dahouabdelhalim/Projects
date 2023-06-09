## R code to produce Figure 2 in:
## Johnston EC, Cunning R, Burgess SC. Cophylogeny and specificity between cryptic coral species (Pocillopora spp.) at Moâ€™orea and their symbionts (Symbiodiniaceae).
## Code written by Erika Johnston. February 2022. Send comments or corrections to ejohnston@bio.fsu.edu
## R version 4.1.2 (2022-03-01)

# initialization ----------------------------------------------------------
rm(list=ls()) #clear all variables

# Load libraries
library(vcfR)
library(adegenet)
library(poppr)
library(ape)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

# Load SNP data
Poc_SNPs.VCF <- read.vcfR("Figure 2.vcf")

# Import sample ID data
samps.data <- read.table("Figure 2.txt", sep = "\\t", header = TRUE)

### Make sure sample names are correct in both files -- Samples need to be in same order
# First, convert to genlight object below, then run gl.Poc_SNPs@ind.names to reorder samples in .txt file
all(colnames(Poc_SNPs.VCF@gt) [-1] == samps.data$Sample_ID)

# Convert to genlight object
gl.Poc_SNPs <- vcfR2genlight(Poc_SNPs.VCF)

# Specify ploidy
ploidy(gl.Poc_SNPs) <- 2

# Add haplotype/species info to genlight object
pop(gl.Poc_SNPs) <- samps.data$Hap_Spp

glPca(gl.Poc_SNPs)

## Find clusters -- Choose a value more than expected clusters: 15 chosen because there are 11 haplotypes
Poc_SNPs.BIC <- find.clusters(gl.Poc_SNPs, max.n.clust=15, scale=FALSE) # Choose 40 PCs because explains most data and 4 clusters because of lowest BIC score

#### First run DAPC with 4 clusters and 40 PCs based on the above
Poc_SNPs.dapc <- dapc(gl.Poc_SNPs, n.da=4, n.pca=40, scale=FALSE)
# Then optimize DAPC
temp <- optim.a.score(Poc_SNPs.dapc) # output suggests 7 PCAs so run dapc again with 7 
Poc_SNPs.dapc <- dapc(gl.Poc_SNPs, n.da=4, n.pca=7, scale=FALSE)

# Generate sample assignment predictions
Poc_SNPs.pred <- predict(Poc_SNPs.dapc)

### Extract sample assignments
tmp <- as.data.frame(Poc_SNPs.dapc$posterior)
tmp$Sample_ID <- rownames(tmp)
tmp <- melt(tmp, id = c("Sample_ID"))
names(tmp)[2:3] <- c("Group", "Posterior")
new <- merge(tmp, samps.data, by="Sample_ID") # add Hap_Spp data to tmp
new$Sample_ID <- gsub("COL_", "", new$Sample_ID) # remove "COL_" from sample ids

# Reorder Pocs
new$Hap_Spp = factor(new$Hap_Spp,levels=c("Hap10","Hap3a","Hap3b","Hap3f","Hap3h","Hap1a_Pm","Hap8a","Hap1a_Pgra","Hap11","Hap2","Hap6a"))

# Custom label
Hap.samp.labs <- c("Hap10","Hap3a","Hap3b","Hap3f","Hap3h","P mean","Hap8a","P gra","Hap11","Hap2","P lig")
names(Hap.samp.labs) <- c("Hap10","Hap3a","Hap3b","Hap3f","Hap3h","Hap1a_Pm","Hap8a","Hap1a_Pgra","Hap11","Hap2","Hap6a")

### Plot
ggplot(new, aes(x = Sample_ID, y = Posterior, fill = Group)) +
geom_bar(stat = "identity") +
facet_grid(~Hap_Spp, scales="free", space="free",
           labeller = labeller(Hap_Spp = Hap.samp.labs)) +
theme_bw(base_size = 6) +
ylab("Posterior membership probability") +
theme(legend.position='none') +
scale_fill_manual(name="Species/haplotype",
                  values = c("Hap10"= "#D55E00",
                             "Hap3a"= "#faa7c0",
                             "Hap3b"= "#E69F00",
                             "Hap3f"= "#a36e40",
                             "Hap3h"= "#ed2b2b",
                             "Hap1a_Pm"= "#0072B2",
                             "Hap8a"= "#CC79A7",
                             "Hap1a_Pgra"= "#56B4E9",
                             "Hap11"= "#009E73",
                             "Hap2"= "#B2DF8A",
                             "Hap6a"= "#266637"),
                  labels=c("Haplotype 10","P. verrucosa 3a","P. verrucosa 3b","P. verrucosa 3f","P. verrucosa 3h","P. meandrina","Haplotype 8a","P. grandis","Haplotype 11","Haplotype 2","P. ligulata"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7), axis.title.x=element_blank())


# Determine what percentage of genetic variance is explained by each axis
percent = Poc_SNPs.dapc$eig/sum(Poc_SNPs.dapc$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,70),
        names.arg = round(percent, 1))

# Create a data frame containing individual coordinates
ind_coords = as.data.frame(Poc_SNPs.dapc$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3","Axis4")

# Add a column containing individuals
ind_coords$Ind = samps.data$Sample_ID

# Add a column with Haplotype/species
ind_coords$haps = samps.data$Hap_Spp

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3, Axis4) ~ haps, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords <- merge(ind_coords, centroid, by="haps", suffix = c("",".cen"))

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Order Pocs
ind_coords$haps = factor(ind_coords$haps,levels=c("Hap10","Hap3a","Hap3b","Hap3f","Hap3h","Hap1a_Pm","Hap8a","Hap1a_Pgra","Hap11","Hap2","Hap6a"))

# Plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen), show.legend = F)+
  # points
  geom_point(aes(fill = haps), shape = 21, size = 3, show.legend = T)+
  scale_fill_manual(name="Species/haplotype",
                    values = c("Hap10"= "#D55E00",
                               "Hap3a"= "#faa7c0",
                               "Hap3b"= "#E69F00",
                               "Hap3f"= "#a36e40",
                               "Hap3h"= "#ed2b2b",
                               "Hap1a_Pm"= "#0072B2",
                               "Hap8a"= "#CC79A7",
                               "Hap1a_Pgra"= "#56B4E9",
                               "Hap11"= "#009E73",
                               "Hap2"= "#B2DF8A",
                               "Hap6a"= "#266637"),
                    labels=c("Haplotype 10","P. verrucosa 3a","P. verrucosa 3b","P. verrucosa 3f","P. verrucosa 3h","P. meandrina","Haplotype 8a","P. grandis","Haplotype 11","Haplotype 2","P. ligulata"))+
  labs(x = xlab, y = ylab)+
  theme_bw()


# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 3 (", format(round(percent[3], 1), nsmall=1)," %)", sep="")

# Plot axis 1 vs. 3
ggplot(data = ind_coords, aes(x = Axis1, y = Axis3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis3.cen), show.legend = F)+
  # points
  geom_point(aes(fill = haps), shape = 21, size = 3, show.legend = T)+
  scale_fill_manual(name="Species/haplotype",
                    values = c("Hap10"= "#D55E00",
                               "Hap3a"= "#faa7c0",
                               "Hap3b"= "#E69F00",
                               "Hap3f"= "#a36e40",
                               "Hap3h"= "#ed2b2b",
                               "Hap1a_Pm"= "#0072B2",
                               "Hap8a"= "#CC79A7",
                               "Hap1a_Pgra"= "#56B4E9",
                               "Hap11"= "#009E73",
                               "Hap2"= "#B2DF8A",
                               "Hap6a"= "#266637"),
                    labels=c("Haplotype 10","P. verrucosa 3a","P. verrucosa 3b","P. verrucosa 3f","P. verrucosa 3h","P. meandrina","Haplotype 8a","P. grandis","Haplotype 11","Haplotype 2","P. ligulata"))+
  labs(x = xlab, y = ylab)+
  theme_bw()
