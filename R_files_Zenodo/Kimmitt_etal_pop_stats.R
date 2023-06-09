# title: Kimmitt et al. "Genetic evidence for widespread population size expansion in North American boreal birds prior to the Last Glacial Maximum" Population Statistics
# author: abby kimmitt, ben winger
#this script is evaluating Tajima's D and Fst

#load packages 
library(ape)
library(PopGenome)
library(adegenet)
library(hierfstat)
library(pegas)

#clean environment
rm(list=ls())

#set working directories 
snippet=#snippet is the 5 letter code for the species 

#define wd1 as where the fasta files are located
#define wd2 as where the metadata on the files are located


#read in the fasta alignment 
faspath <- paste0(wd1,"/",snippet,".fasta")
alignment <- read.dna(faspath, format="fasta")

#Tajima's D
d<-tajima.test(alignment)

# Import metadata
mpath<- paste0(wd2,"/","Kimmitt_etal_sample_data.csv")
meta <- read.csv(mpath, header = T)

# check to make sure names of all samples match meta file 
name_check <- rownames(alignment)[which(rownames(alignment) %in% meta$fasta_label == F)]
print(name_check)

names<-rownames(alignment)
# append meta data to rownames
append.meta <- function(Seq){
  for(i in 1:nrow(Seq)){
    if(any(meta$fasta_label %in% rownames(Seq)[i]) == T){
      data <- meta[meta$fasta_label == rownames(Seq)[i],]
      rownames(Seq)[i] <- paste(data$fasta_label, data$population, sep = "_") ### add desired metadata here
    }
  }
  return(Seq)
}

all<- append.meta(alignment) # all populations

# create population labels 
Labels <- rownames(all)
CN <- grep("CN", Labels, ignore.case = F) #MI and MN
NE <- grep("NE", Labels, ignore.case = F) #NY and VT (Northeast)
MB <- grep("MB", Labels, ignore.case = F) #Manitoba
AB <- grep("AB", Labels, ignore.case = F) #Alberta

# Convert alignment to Genind object, preserving populations defined above
pop.factor <- vector("character", nrow(all)) 
pop.factor[CN] <- "CN"
pop.factor[NE] <- "NE"
pop.factor[MB] <- "MB"
pop.factor[AB] <- "AB"

pop.factor <- as.factor(pop.factor)

a.genind <- DNAbin2genind(all, pop = pop.factor)
a.genpop <- genind2genpop(a.genind, pop = pop.factor)
a.genpop

#using hierfstat and adgenet to look at Fst 
#can't use Fst function in pegas on haploid data
a.hier<-genind2hierfstat(a.genind)
fst<-pairwise.WCfst(a.hier, diploid=F)
