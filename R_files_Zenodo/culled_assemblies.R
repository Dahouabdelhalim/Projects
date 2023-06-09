#These assemblies were culled using different thresholds for contig size 
#generating fasta files composed of 1Kb contigs, 1.5kb, 2Kb, 2.5Kb:
#written with input from Dr. Matthew Settles at University of California Davis
  library(Biostrings)
fa<- readDNAStringSet("assembled_contigs.fasta")
fa
#lists contents of fasta (output from assembly).
#Starts with "A DNAStringSet instance of length 269287"
#Note that there is a column for width. We'll sort by this column.
fa2<- fa[width(fa)>=1000]
fa2
#Now starts with "A DNAStringSet instance of length 34669"
#only those contigs 1000 bp or greater
writeXStringSet(fa2,"/home/mlatvis/folder/culled_assembly.fasta")
