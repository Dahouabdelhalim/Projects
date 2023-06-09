####Script for constructing SNP-based pedigree for Helgeland house sparrows.
#Alina Niskanen
#Modified May 2020

require(sequoia)

#Bring the genotype data in RAW format
#Convert the genotype data into right format for Sequoia
Genotypes_temp <- GenoConvert(InFile = "Genotypes_sequoia.raw")
Individuals_temp <- rownames(Genotypes_temp)

#Import lifehistory data
LifeHistData <- read.table("Sequoia_lifehistory.txt", header = T, stringsAsFactors = F)


#Include only individuals that are in the lifehistory data 
Genotypes <- Genotypes_temp[rownames(Genotypes_temp) %in% LifeHistData[,1],]
Individuals <- rownames(Genotypes)


#Perform pedigree reconstruction based on SNP data, include also sibship clustering
PedOUT <- sequoia(GenoM = Genotypes, LifeHistData = LifeHistData, SeqList = NULL, MaxSibIter = 10, Err = 0.002, MaxMismatch = 3, Tfilter = -2, Tassign = 0.5, MaxSibshipSize = 100, DummyPrefix = c("F", "M"), Complex = "full", FindMaybeRel = TRUE, CalcLLR = TRUE, quiet = FALSE) 

names(PedOUT)
head(PedOUT$Specs)

#Output the pedigree
write.table(PedOUT$Pedigree, "pedigree.txt", col.names = T, row.names = F, quote = F)





