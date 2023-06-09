
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(phyloseq); packageVersion("phyloseq")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path<-"PacBio data/McMunnFullITS/"
fns <- list.files(path,  pattern="fastq")

ITS1FKYO2 <- "TAGAGGAAGTAAAAGTCGTAA"
ITS4KYO1 <- "TCCTCCGCTTWTTGWTWTGC"

rc <- dada2:::rc
theme_set(theme_bw())

path.out <- "figures"
path.rds <- "RDS"

nops <- file.path(path, "noprimers", basename(fns))
prim <- removePrimers(paste(path,fns,sep=""), nops, primer.fwd=ITS1FKYO2, primer.rev=dada2:::rc(ITS4KYO1), orient=TRUE)

filts <- file.path(path, "noprimers", "filtered", basename(fns))
track <- filterAndTrim(nops, filts, minQ=3, minLen=200, maxLen=850, maxN=0, rm.phix=FALSE, maxEE=2)

lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

drp <- derepFastq(filts, verbose=TRUE)


err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)

saveRDS(err, file.path(path.rds, "SCC_BB_SF_ITS_errors.rds"))

plotErrors(err)

dd <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)
cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2], denoised=sapply(dd, function(x) sum(x$denoised)))

st <- makeSequenceTable(dd)
st

tax <- assignTaxonomy(st, "PacBio data/McMunnFullITS/training/sh_general_release_dynamic_s_all_04.02.2020.fasta", multithread=TRUE) # Slowest part
head(unname(tax))

bim <- isBimeraDenovo(st, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim)
sum(st[,bim])/sum(st)

dim(st)
sample.names <- sapply(strsplit(fns, "-"), function(x) paste(x[1], x[2], sep="-"))
rownames(st) <- sample.names

saveRDS(st, file.path("RDS", "SCCBBSF_PacBio_st_ITS.rds"))
saveRDS(tax, file.path("RDS", "SCCBBSF_PacBio_UNITE_ITS.rds"))

write.csv(tax, "taxCheckITS.csv")
write.csv(st, "sampCheckITS.csv")
its<-phyloseq(otu_table(st, taxa_are_rows = FALSE), tax_table(tax))

saveRDS(its , "ITS_DADA2output")

