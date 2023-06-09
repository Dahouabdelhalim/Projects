library(vcfR)
library(raster)
library(foreign)
library(stringr)
setwd("c:/Users/jsons/Dropbox/R_projects/Phenig_reanalysis/final_data/Uploaded_to_dryad/")
all <- read.vcfR("./Phenig_final_biallele_SNP_maf005_annotated.recode.vcf")
geno <- extract.gt(all) # Character matrix containing the genotypes

position <- getPOS(all) # Positions in bp
chromosome <- getCHROM(all) # Chromosome information
G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2

table(as.vector(G))
dim(G)

pop <- read.table("./Phenig.pop.txt", header = F) #Load the individuals with the population
colnames(pop)<- c("Ind", "pop")
env<- read.csv("Environmental_variation.csv", header = T)#Load the environmental variation
rownames(pop)<- pop$Ind
rownames(env)<- env$ind
env <- env[rownames(pop),]

identical(env$ind, pop$Ind) #check the order of the environmental matrix and the pop matrix
identical(colnames(geno), pop$Ind)  #check the order of the genotype matrix and the pop matrix

chromosome<-as.numeric(str_split_fixed(chromosome, "_",n = 2)[,2])#Remove the text "scaffold" from the chromosome and make the data numeric 
#We put all the data into a list
Phenig <- list("SNP"=paste(chromosome,"_",position, sep = ""),"position"=position, "chromosome"=chromosome, "G"=G, "pop"= pop$pop, "id"=pop$Ind,"env" = env)
