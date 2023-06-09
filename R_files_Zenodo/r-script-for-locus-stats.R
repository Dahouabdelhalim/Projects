# R script to calculate mean bootstrap support and GC variance across taxa
# for a set of gene trees and associated alignments.

# Script modified from Borowiec et al. (2015).
# Borowiec, M.L., Lee, E.K., Chiu, J.C. & Plachetzki, D.C. (2015) Extracting phylogenetic 
# signal and accounting for bias in whole-genome data sets supports the Ctenophora as 
# sister to remaining Metazoa. BMC Genomics, 16, 987.

# Alignment files should be in fasta format and deposited within a single folder 
# labelled "Single_gene_alignments". Files should be named "uce-#.fasta". 

# Gene-tree files should be in Newick format, have support values, and be deposited 
# within a single directory labelled "RAxML_single_genes". Files should be named 
# "RAxML_bipartitions.uce-#".

# Install required libraries.
install.packages("ape")
install.packages("seqinr")
install.packages("data.table")

# Load needed libraries.
library("ape")
library("seqinr")
library("data.table")

# Set working directories and load all tree files with support.
setwd("/working/directory/")
getwd()

# For the relative paths to work correctly.
trees_dir <- file.path("RAxML_single_genes/")
trees_bip <- dir(path=trees_dir, pattern="*bipartitions.uce-*")

# Load all alignment files in FASTA format
fasta_dir <- file.path("Single_gene_alignments/")
alignments <- dir(path=fasta_dir, pattern="*.fasta")

### MEAN BOOTSTRAP SUPPORT ###

# This code should work with any Newick tree files
# that have some measure of support at nodes, 
# including PP from Bayesian analysis.

# define a function to calculate average support
Avg_support <- function(file) {
  
  # get locus number from filename
  locus_no <- sub("RAxML_bipartitions.(uce-[0-9]+)", "\\\\1", perl=TRUE, x=file)
  # read in tree
  tree <- read.tree(paste(trees_dir,file, sep=""))
  # store support values in a vector
  support <- c(as.numeric(tree$node.label))
  # calculate average support
  avg_supp <- mean(support, na.rm=T)
  return(c(locus_no,avg_supp))
  
}

# loop over all files
average_bootstrap <- lapply(trees_bip, Avg_support)
average_bootstrap <- data.frame(matrix(unlist(average_bootstrap), nrow=(length(average_bootstrap)), byrow=T))
colnames(average_bootstrap) <- c("locus", "Average_bootstrap")
write.table(average_bootstrap,"average-bootstrap.txt",sep="\\t")

#### BASE COMPOSITION & HETEROGENEITY ###

# Define a function to calculate base composition, gc content, and gc variance

gc_content <- function(file) {

	alignment <- read.alignment(file=paste(fasta_dir, file, sep=""), format="fasta")
	locus_no <- sub("(uce-[0-9]+).fasta", "\\\\1", perl=TRUE, x=file)
	alignbin <- as.DNAbin(alignment)
	gc <- GC.content(alignbin)
	base_comp <- base.freq(alignbin)
	gc_var <- var(sapply(alignbin, GC.content))
	return(list(locus_no,gc,base_comp, gc_var))

}

# loop over all files
bases <- lapply(alignments, gc_content)
bases <- data.frame(matrix(unlist(bases), nrow=(length(bases)), byrow=T))
colnames(bases) <- c("locus", "gc_content", "a","c","g","t", "gc_var")
write.table(bases,"base-comp.txt",sep="\\t")