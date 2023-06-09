## R code to produce Figure 2 in: 
## Johnston EC, Wyatt ASJ, Leichter JJ, Burgess SC. Niche differences in co-occuring cryptic coral species (Pocillopora spp). Coral Reefs.
## Code written by Erika Johnston. March 2021. Send comments or corrections to ejohnston@bio.fsu.edu
## R version 3.6.2 (2021-02-29)


#load packages
library(ape)
library(pegas)
library(adegenet)
library(dplyr)
library(RColorBrewer)

### Import haplotype data
haplotypes_r_Aug2019 <- read.csv("Figure 2 Data.csv")

### Import fasta data
Moo_Aug2019 <- read.dna("Figure 2 mtORF.fasta", format="fasta")

Moo_Aug2019_haps <- haplotype(Moo_Aug2019)
Moo_Aug2019_haps

count(haplotypes_r_Aug2019, "Species.haplotype")
Moo_Aug2019_Net <- haploNet(Moo_Aug2019_haps)

hap.hap_Aug2019 <- table(haplotypes_r_Aug2019$hap, haplotypes_r_Aug2019$Species.haplotype)

### Plot haplotype network first by haplotype number
plot(Moo_Aug2019_Net, size = attr(Moo_Aug2019_Net,"freq")*.03, scale.ratio = 0.6,
     cex = 0.7, lty = 1, lwd = 1, show.mutation = 2, pie = NULL)


### Set colors
cols_Aug2019 <- c("Haplotype 10"= "#D55E00",
                  "Haplotype 3b"= "#E69F00",
                  "Haplotype 1a_Pe"= "#56B4E9",
                  "Haplotype 1a_Pm"= "#0072B2",
                  "Haplotype 1c"= "#808080",
                  "Haplotype 1d"= "#808080",
                  "Haplotype 8a"= "#CC79A7",
                  "Haplotype 3h"= "#808080",
                  "Haplotype 3f"= "#808080",
                  "Haplotype 3a"= "#808080",
                  "Haplotype 11"= "#009E73",
                  "Haplotype 2"= "#B2DF8A",
                  "Haplotype X"= "#808080",
                  "Haplotype 9"= "#808080",
                  "Haplotype 1e"= "#808080")

### Plot
quartz(width=8,height=8)
plot(Moo_Aug2019_Net, size = attr(Moo_Aug2019_Net,"freq")*.03, scale.ratio = 0.6,
     cex = 0.7, lty = 1, lwd = 1, show.mutation = 2, pie = hap.hap_Aug2019, legend=FALSE, labels=FALSE, bg=cols_Aug2019)
