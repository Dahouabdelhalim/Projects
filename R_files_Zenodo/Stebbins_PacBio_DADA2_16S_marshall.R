
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
path<-"PacBio data/McMunnFull16S/"
fns <- list.files(path,  pattern="fastq")

F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"

rc <- dada2:::rc
theme_set(theme_bw())

path.out <- "figures"
path.rds <- "RDS"

nops <- file.path(path, "noprimers", basename(fns))
prim <- removePrimers(paste(path,fns,sep=""), nops, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE)

filts <- file.path(path, "noprimers", "filtered", basename(fns))
track <- filterAndTrim(nops, filts, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)

lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

ilts <- file.path(path2, "noprimers", "filtered", basename(fns))
track <- filterAndTrim(nops, filts, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)

drp <- derepFastq(filts, verbose=TRUE)


err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)

saveRDS(err, file.path(path.rds, "SCC_BB_SF_errors.rds"))

plotErrors(err)

dd <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)
cbind(ccs=prim[,1], primers=prim[,2], filtered=track[,2], denoised=sapply(dd, function(x) sum(x$denoised)))

st <- makeSequenceTable(dd)

tax <- assignTaxonomy(st, "PacBio data/McMunnFull16S/Training/silva_nr_v138_train_set.fa.gz", multithread=TRUE) # Slowest part
head(unname(tax2))

bim <- isBimeraDenovo(st, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim)
sum(st[,bim])/sum(st)


sample.names <- sapply(strsplit(fns, "-"), function(x) paste(x[1], x[2], sep="-"))
rownames(st) <- sample.names

saveRDS(st, file.path("RDS", "SCCBBSF_PacBio_st.rds"))
saveRDS(tax, file.path("RDS", "SCCBBSF_PacBio_Silva128.rds"))

write.csv(tax, "taxCheck.csv")
write.csv(st, "sampCheck.csv")
bact16s<-phyloseq(otu_table(st, taxa_are_rows = FALSE), tax_table(tax))

saveRDS(bact16s , "16S_DADA2output")
