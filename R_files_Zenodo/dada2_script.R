rm(list=ls())
gc()

## Load packages

library(ggplot2)
library(tidyr)
library(Rmisc)
library(scales)
library(Hmisc)
library(stringr)
library(dplyr)
library(MASS)
library(Rcpp)
library(dada2)
library(reshape2)

## Connect to source folder

path <- "C:/Users/chammoud/Desktop/Miseqv3.2/no_primer/"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_R1_cut.fastq"), `[`, 1)

## Plot quality before trimming

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

## Trim and filter reads

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE=c(2,2), rm.phix=TRUE, minLen = 150,compress=TRUE, multithread=FALSE) 
head(out)

## Plot quality after trimming

plotQualityProfile(filtindFs[1:2])
plotQualityProfile(filtRs[1:2])

## Learn error rates

errF <- learnErrors(filtFs, multithread=TRUE,randomize=TRUE,nbases = 5e8)
errR <- learnErrors(filtRs, multithread=TRUE,randomize=TRUE,nbases = 5e8)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## Dereplicate reads

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Run dada2 algorithm to call variants

dadaFs <- dada(derepFs, err=errF, pool = FALSE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = FALSE, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

## Merge reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,minOverlap = 20, maxMismatch = 1,  verbose=TRUE, returnRejects= FALSE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

## Remove chimera

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track=data.frame(track)
track=track%>%mutate(prop_used=nonchim/input)
sum(track[1:447,]$nonchim)/sum(track[1:447,]$input)

## Write unique ASVs in fasta

seqnum <- paste0("Seq", seq(ncol(seqtab.nochim)))
uniqueSeqs <- as.list(colnames(seqtab.nochim))
library("seqinr")
write.fasta(uniqueSeqs, seqnum, "uniqueSeqs.fasta")

## Reshape ASVs table in long format

seqtab_nochim=seqtab.nochim
seqtab_nochim$sample=row.names(seqtab_nochim)
rownames(seqtab_nochim) <- c()
seqtab_nochim_long=seqtab_nochim%>%gather(sequence,count,1:477)
seqtab_nochim_long$shortID=sub("\\\\..*", '', seqtab_nochim_long$sample)
seqtab_nochim_long$shortID=sub("GC07942", '', seqtab_nochim_long$shortID)

