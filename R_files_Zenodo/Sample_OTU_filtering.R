#-----------------PACKAGES---------------------
#packages you will need
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("phyloseq", "DECIPHER", "phangorn")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

library(vegan)
library(car)
library (DESeq2); packageVersion("DESeq2")
library(Rfast) #to take row maximums
library (ggthemes)
library("data.table"); packageVersion("data.table")
library(adaptiveGPCA)
library(ggrepel)
library(ggplot2)
library(phyloseq)
# if there is an instolation probelem use this: chooseCRANmirror()
##-- Rad the phyloseq ---------------------------------------------------------------------------------------
data <- readRDS("Israel_microbiome_run2_crane_phyloseq_data_Updated.rds")
#Get rid of zero taxa in the crane data (not carnes)
data = filter_taxa(data, function(x) sum(x) > 0, TRUE)
data

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##+++++++++++++ NEGATIVE CONTROL+++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###=== Decontam method ========================================================================
# we work with a phyloseq containing all cranes (also trapped) and the negative controls

data1 <- subset_samples(data, Species%in% "Crane" | SampleID%in%c("Negative_control_Ammon_7-10-18_crane",
                                                                "Negative_control_Candace_7-21-18_crane", 
                                                                "Negative_control_Ammon_8-7-18_crane", "H2O_blank_Crane_plate"))
a1<-data.frame(data1@sam_data)
#install.packages("decontam")
library("decontam")
#Let's try removing contaminants using decontam package
#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

#create column with "true" in negative control
sample_data(data1)$is.neg <- sample_data(data1)$SampleID%in%c("Negative_control_Ammon_7-10-18_crane",
                                                            "Negative_control_Candace_7-21-18_crane", 
                                                            "Negative_control_Ammon_8-7-18_crane",
                                                            "H2O_blank_Crane_plate")
#PREVALENCE: Contaminants are identified by increased prevalence in negative controls
prev5 <- isContaminant(data1, method="prevalence", neg="is.neg", threshold = 0.5)
table(prev5$contaminant) #summary: 27 true==contaminent: 0.3% of total, 16.7% neg cong (161 ASVs)

## STEPS:
# we DON'T remove rare ASVs anymore
# We will remove th eContaminants later based on this endetificarion
# Remove non bacteria
# create phyloseq of unclsified and remove them
# for the most prevelent unclasified do blast (google it) to show that tehy are host DNA
# https://blast.ncbi.nlm.nih.gov/Blast.cgi --> choose nucleoyide BLAST and and insert the real sequense (OTU table)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++ Remove irelevant samples ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# choose what samples we want
filtered.data.pared.down.Feces <- subset_samples(data, Species%in% "Crane" & 
                                                   Duplicate%in%  "no" & PCR%in% "Good" & Trapped%in% "No" &
                                                   Status4Analysis%in% "High_priority")
# choose what samples we don't want (Soil Sample)
filtered.data.pared.down.Feces <- subset_samples(filtered.data.pared.down.Feces, extraction_id!="Gg235")
# Choose only taxa relevant for our comparisons:
#data_in <- subset_samples(filtered.data.pared.down.Feces, (Year=="2017" | Year=="2018") & (WeekSinceBreeding>2 & WeekSinceBreeding<7 | 
                                                                       # WeekSinceBreeding>10 & WeekSinceBreeding<16 |
                                                                       # WeekSinceBreeding>20 & WeekSinceBreeding<29  |
                                                                       # WeekSinceBreeding>30 & WeekSinceBreeding<35))

data_in <- subset_samples(filtered.data.pared.down.Feces, (Year=="2017" | Year=="2018") & (WeekSinceBreeding>2 & WeekSinceBreeding<7 | 
                                                              WeekSinceBreeding>10 & WeekSinceBreeding<16 |
                                                                WeekSinceBreeding>20 & WeekSinceBreeding<29))
#-------------------------------------------------
# DON'T RUN TBEFORE removing the negative control, otherwise won't work, it's to know the read number
data_in <- filter_taxa(data_in, function(x) sum(x) > 0, TRUE)
# 5123 taxa and 167 samples
#-------------------------------------------------

#NEGATIVE CONTROL: remove reads based on decontam analysis:
data.noncontam.prev5 <- prune_taxa(!prev5$contaminant, data_in)
data.noncontam.prev5 #8527 ASVs
filtered.data.pared.down.Feces <- filter_taxa(data.noncontam.prev5, function(x) sum(x) > 0, TRUE)


filtered.data.pared.down.Feces  # n = 167
#otu_table()   OTU Table:         [5100 taxa and 167 samples]
# the negative control is 27 ASVs, but only 23 of them still present in the samples
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++ Remove non bacteria ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

filtered.data1 <- filtered.data.pared.down.Feces #this is the input for Ammon's code
#Filtering out taxa that should not be in the dataset
#Learned this from: http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#Select just the mitochondrial sequences.  Note that it is case sensitive, so Mitochondria is different than mitochondria.
mitochondrial.seqs <- subset_taxa(filtered.data1, Family == "Mitochondria")
#Get rid of zero taxa 
mitochondrial.seqs.unique = filter_taxa(mitochondrial.seqs, function(x) sum(x) > 0, TRUE)
#86 taxa were mitochondrial
mitochondrial.seqs.unique

#Select cholorplast sequences, but could not find any
chloroplast.seqs <- subset_taxa(filtered.data1, Class == "Chloroplast")
#Error suggests that there just aren't Choloplasts 

nonbacterial.seqs <- subset_taxa(filtered.data1, Kingdom != "Bacteria")
#Get rid of zero taxa 
nonbacterial.seqs.unique = filter_taxa(nonbacterial.seqs, function(x) sum(x) > 0, TRUE)
#found 66 taxa meeting this criterion.  They are all Archaea according to the taxtable
nonbacterial.seqs.unique


##THUS 152 samples must be removed!!!! (Non bacteria and mithondiria)

archaea.and.mitochondria <- subset_taxa(filtered.data1, Kingdom != "Bacteria" | Family == "Mitochondria" | Class == "Chloroplast")
#Get rid of zero taxa 
archaea.and.mitochondria2 = filter_taxa(archaea.and.mitochondria, function(x) sum(x) > 0, TRUE)
#152 taxa
archaea.and.mitochondria2

archaea.and.mitochondria.taxa <- colnames(otu_table(archaea.and.mitochondria2))
length(archaea.and.mitochondria.taxa) #this is what suppose to go (152 OTUs)

#--SOLUTION:---
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

archaea.and.mitochondrial.filtered.data = pop_taxa(filtered.data1,archaea.and.mitochondria.taxa)
#archaea.and.mitochondrial.filtered.data = filter_taxa(archaea.and.mitochondrial.filtered.data, function(x) sum(x) > 0, TRUE)
#Has 4948 taxa.  5100 taxa originally - 152 taxa for Archaea and mitochondria 
archaea.and.mitochondrial.filtered.data

#Assign this new filtered file back to filtered.data1 so that the code below works.
filtered.data.pared.down.Feces <- archaea.and.mitochondrial.filtered.data
#otu_table()   OTU Table:         [4948 taxa and 167 samples]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++ Remove Unassigned Phlya ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# See the prevelne of the different fhylum inclding the NA
prev0 = apply(X = otu_table(filtered.data.pared.down.Feces), #data with C- taxa
              MARGIN = ifelse(taxa_are_rows(filtered.data.pared.down.Feces), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(filtered.data.pared.down.Feces),
                    tax_table(filtered.data.pared.down.Feces))

#Are there phyla that are comprised of mostly low-prevalence features? Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})                      

# Find the taxa with not assigned a phylum
NoPhylum <- subset_taxa(filtered.data.pared.down.Feces, is.na(Phylum) | Phylum %in% c("", "uncharacterized"))
# There are 74 taxa like that 
prev0NP = apply(X = otu_table(NoPhylum), #data with C- taxa
              MARGIN = ifelse(taxa_are_rows(NoPhylum), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdfNP = data.frame(Prevalence = prev0NP,
                    TotalAbundance = taxa_sums(NoPhylum),
                    tax_table(NoPhylum))
# we will now inspect (using BLAST) 5 most abundant and 5 most prevelant ones:
#https://blast.ncbi.nlm.nih.gov/Blast.cgi --> choose nucleoyide BLAST and and insert the real sequense (OTU table)
prevdfNP$revPrevalence = prevdfNP$Prevalence/199
prevdfNP=prevdfNP[order(prevdfNP$Prevalence, decreasing = TRUE),]
mostPrev=prevdfNP[1:5,]
prevdfNP=prevdfNP[order(prevdfNP$TotalAbundance, decreasing = TRUE),]
MostAbund=prevdfNP[1:5,]
# prevelence # Abundance
#       16       173            # Hannaella oryzae (yeast) mitochondrion
#       13       127            # Rhizopus oryzae mitochondrion (microfungus that occurs as a saprotroph in soil, dung, and rotting vegetation)
#       1        103            # Chlorella sorokiniana mitochondrion (freshwater green microalga)
#       14       60             # Uncultured bacterium clone 773 16S
#       14       35             # Rhizopus oryzae mitochondrion
#       7        29             # Hannaella oryzae (yeast) mitochondrion

#Get rid of any taxa that are not assigned a phylum
filtered.data.pared.down.Feces.1 <- subset_taxa(filtered.data.pared.down.Feces, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
#74 OTUs removed : 4948-4874

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++ Remove outliers ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# remove outliers (dry feces in Emek Israel collected with Gabe)
filtered.data.pared.down.Feces.1<- subset_samples(filtered.data.pared.down.Feces.1, SampleID!="Gg247")
filtered.data.pared.down.Feces.1<- subset_samples(filtered.data.pared.down.Feces.1, SampleID!="Gg248")
filtered.data.pared.down.Feces.1<- subset_samples(filtered.data.pared.down.Feces.1, SampleID!="Gg249")
filtered.data.pared.down.Feces.1<- subset_samples(filtered.data.pared.down.Feces.1, SampleID!="Gg250")
filtered.data.pared.down.Feces.1<- subset_samples(filtered.data.pared.down.Feces.1, SampleID!="Gg251")
#Get rid of zero taxa
filtered.data.pared.down.Feces<- filter_taxa(filtered.data.pared.down.Feces.1, function(x) sum(x) > 0, TRUE)
# 4821 taxa and 162 samples


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++ See what we are left with and give unique OTU ID ++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#See how many reads are left
titles_reads=c("sum","mean","st.dev", "median", "min","max")
values_reads=print(c(sum(sample_sums(filtered.data.pared.down.Feces.1)),
                     mean(sample_sums(filtered.data.pared.down.Feces.1)),
                     sd(sample_sums(filtered.data.pared.down.Feces.1)), 
                     median(sample_sums(filtered.data.pared.down.Feces.1)),
                     min(sample_sums(filtered.data.pared.down.Feces.1)),
                     max(sample_sums(filtered.data.pared.down.Feces.1))))
reads_data=data.frame(titles_reads,values_reads)
reads_data

## THIS WAS WITH NEGATIVE CONTROL REMOVED (The old wrong way)
#titles_reads values_reads
#1          sum   767548.000
#2         mean     3857.025
#3       st.dev     2922.300
#4       median     3267.000
#5          min      428.000
#6          max    17021.000

## THIS IS THE CORRECT ONE DONE WITH decontam

#   titles_reads values_reads
#1          sum   2261851.00
#2         mean     13962.04
#3       st.dev      4377.14
#4       median     13844.50
#5          min      5547.00
#6          max     25670.00
#====Add OTU ID
CleanData=filtered.data.pared.down.Feces
OTU<-otu_table(CleanData)
TAX<-data.frame(tax_table(CleanData))
x <- 1:nrow(TAX)
mynames <- as.character(x)
tax_table(CleanData) <- cbind(tax_table(CleanData), OTUID=mynames)

# SAVE The clean data
saveRDS (CleanData, file="Crane_CleanDataNewAproach.rds") 
