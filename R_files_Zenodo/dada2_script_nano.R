rm(list=ls())
gc()

library(ggplot2)
library(tidyr)
library(Rmisc)
library(scales)
# library(plotrix)
library(Hmisc)
# library(corrplot)
# library(MuMIn)
library(stringr)
# library(reshape2)
library(dplyr)
library(MASS)
library(Rcpp)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.9")
library(dada2)
library(reshape2)

path <- "C:/Users/chammoud/Desktop/nano2/no_primer/"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1_cut.fastq"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])



filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), rm.phix=TRUE, minLen = 150,compress=TRUE, multithread=FALSE) 
head(out)
# 
# path2 <- "C:/Users/chammoud/Desktop/Miseqv3.2/indexes/"
# indFs <- sort(list.files(path2, pattern=".I1.fastq", full.names = TRUE))
# indRs <- sort(list.files(path2, pattern=".I2.fastq", full.names = TRUE))
# filtindFs <- file.path(path2, "filtered", paste0(sample.names, "_IF_filt.fastq.gz"))
# filtindRs <- file.path(path2, "filtered", paste0(sample.names, "_IR_filt.fastq.gz"))
# 
# out_index <- filterAndTrim(indFs, filtindFs, indRs, filtindRs, maxN=0, maxEE=c(1,1), truncQ=30,rm.phix=TRUE, minLen = 8,compress=TRUE, multithread=FALSE) 
# head(out_index)

# plotQualityProfile(filtindFs[1:2])
# plotQualityProfile(filtRs[1:2])


errF <- learnErrors(filtFs, multithread=TRUE,randomize=TRUE,nbases = 2e8)
errR <- learnErrors(filtRs, multithread=TRUE,randomize=TRUE,nbases = 2e8)

# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, pool = FALSE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = FALSE, multithread=TRUE)


dadaFs[[1]]
dadaRs[[1]]

# mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,minOverlap = 10, maxMismatch = 1,  verbose=TRUE, returnRejects= TRUE)
# head(mergers[[1]])

seqtabFWD <- makeSequenceTable(dadaFs)
dim(seqtabFWD)
table(nchar(getSequences(seqtabFWD)))

seqtabFWD.nochim <- removeBimeraDenovo(seqtabFWD, method="consensus", verbose=TRUE)

dim(seqtabFWD.nochim)
sum(seqtabFWD.nochim)/sum(seqtabFWD)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), rowSums(seqtabFWD.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "nonchim")
rownames(track) <- sample.names
head(track)

seqnum <- paste0("Seq", seq(ncol(seqtabFWD.nochim)))
uniqueSeqs <- as.list(colnames(seqtabFWD.nochim))
library("seqinr")
write.fasta(uniqueSeqs, seqnum, "seqtabFWD.nochim_nano2.fasta")


seqtabFWD_nochim=data.frame(seqtabFWD.nochim)
seqtabFWD_nochim$sample=row.names(seqtabFWD_nochim)
rownames(seqtabFWD_nochim) <- c()
seqtabFWD_nochim_long=seqtabFWD_nochim%>%gather(sequence,count,1:237)
seqtabFWD_nochim_long=seqtabFWD_nochim_long%>%filter(count!=0)



seqtabREV <- makeSequenceTable(dadaRs)
dim(seqtabREV)
table(nchar(getSequences(seqtabREV)))

seqtabREV.nochim <- removeBimeraDenovo(seqtabREV, method="consensus", verbose=TRUE)

dim(seqtabREV.nochim)
sum(seqtabREV.nochim)/sum(seqtabREV)

seqtabREV_nochim=data.frame(seqtabREV.nochim)
seqtabREV_nochim$sample=row.names(seqtabREV_nochim)
rownames(seqtabREV_nochim) <- c()
seqtabREV_nochim_long=seqtabREV_nochim%>%gather(sequence,count,1:206)
seqtabREV_nochim_long=seqtabREV_nochim_long%>%filter(count!=0)

seqtabFWD_nochim_long_A1 = seqtabFWD_nochim_long%>%filter(grepl("_A1.18",sample))
seqnum <- paste0("Seq", seq(length(seqtabFWD_nochim_long_A1)))
uniqueSeqs <- as.list(seqtabFWD_nochim_long_A1$sequence)
library("seqinr")
write.fasta(uniqueSeqs, seqnum, "seqtabFWD_nochim_long_A1.fasta")



mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,minOverlap = 15, maxMismatch = 1,  verbose=TRUE, returnRejects= TRUE)
head(mergers[[1]])

seqta <- makeSequenceTable(mergers)
dim(seqta)
table(nchar(getSequences(seqta)))

seqtab.nochim <- removeBimeraDenovo(seqta, method="consensus", verbose=TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqta)
seqtab_nochim=data.frame(seqtab.nochim)
seqtab_nochim$sample=row.names(seqtab_nochim)
rownames(seqtab_nochim) <- c()
seqtab_nochim_long=seqtab_nochim%>%gather(sequence,count,1:219)
seqtab_nochim_long=seqtab_nochim_long%>%filter(count!=0)


seqnum <- paste0("Seq", seq(ncol(seqtab.nochim)))
uniqueSeqs <- as.list(colnames(seqtab.nochim))
library("seqinr")
write.fasta(uniqueSeqs, seqnum, "uniqueSeqs_nano2.fasta")

BLASTED <-  read.csv("C:/Users/chammoud/Desktop/Miseqv3.2/input/BLASTED.csv", sep = ";", dec = ",")
