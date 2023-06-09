# This is code to replicate the analysis from Alleman et al. "Tandem-Running and Scouting Behavior is Characterized by Up-Regulation of Learning and Memory formation genes within the Ant Brain"
# GO enrichment analysis of T. longispinosus
# Transcriptome constructed using Trinity
# Differentially expressed contigs obtained using DESeq2
# Modules of contigs obtained using WGCNA
# Script by M. Stoldt

####################################################################################################################
## Loading required packages
####################################################################################################################

library(topGO)
library(evaluate)

####################################################################################################################
## Establishment of directories
####################################################################################################################

# Specify working directory for project.
setwd("/My/Working/Directory/")

####################################################################################################################
## Run TopGO 
####################################################################################################################

# Read filtered countsmatrix
load("matrix_cook_Tlongi.RData")
contigs<-as.data.frame(rownames(matrix_cook))

# Read in the mapping file with Contigs associated to GO terms produced by Interproscan
geneID2GO<-readMappings("Tlongi_Interproscan_Output.txt")
geneID2GO<-geneID2GO[match(contigs$contigs,names(geneID2GO))]
Universe_genelist = names(geneID2GO)

# choose class of GO terms (e.g. "BP"= Biological Process)
ontology_type = "BP"

####################################################################################################################
## For DEGs 
####################################################################################################################

# Set directory where your input files (DEGs) are
direc<-"/Directory/DEGs/"
# List all files in this directory
resFiles <- list.files(direc, include.dirs = FALSE)

# Do for every file in that folder
for (i in 1:length(resFiles)) {
  
  # Ceate a nice output name
  fileName <- resFiles[i]
  split_x<-strsplit(resFiles[i], "[.]")
  split_y<-strsplit(split_x[[1]][1], "_")
  output_name<-paste(split_y[[1]][2], "_", split_y[[1]][3],"_",split_y[[1]][4], sep="")
  
  # Read in the file
  help<-read.table(paste(direc,fileName, sep=""),header=TRUE, sep=",")
  # Only use the gene IDs for further analysis
  gene_IDs <- as.character(help$ID)
  # See where transcriptome and gene list overlap
  gene_list <- factor(as.integer(Universe_genelist %in% gene_IDs))
  names(gene_list)<-Universe_genelist
  
  # Create a topGO object
  GO_data <- new("topGOdata", description=paste("GO_Enrichment_",output_name,sep = ""), ontology=ontology_type, allGenes= gene_list, annot = annFUN.gene2GO, gene2GO=geneID2GO)
  # Run fishers exact test using weight01 algorithm 
  result_topGO <- runTest(GO_data, algorithm = "weight01", statistic = "fisher",topNodes = numSigGenes(GO_data))
  # Write results into gene table
  result_table <- try(GenTable(GO_data, Fisher = result_topGO, orderBy = "Fisher", ranksOf = "Fisher", topNodes=length(result_topGO@score)))
  if(inherits(result_table, "try-error")){
    next
  }
  # Write results into file
  write.csv(result_table, file=paste(getwd(),"/TopGO_Results_",output_name,".csv", sep=""))
}

####################################################################################################################
## For WGCNA 
####################################################################################################################

# Do same as described above using the directory where the module contig lists are in 