BiocManager::install("topGO")
BiocManager::install("tagcloud")
library("topGO")
library("RColorBrewer")
library("tagcloud")


GO_universe <- readMappings("Mnat_go_corr") # use reduced input file here containing only the orthologs in the counts table.
interesting_genes <- #list of genes of interest.
comp <- # give name of comparison for labeling the output file, e.g. SminVT0

  
geneID2GO <- GO_universe

allRes_BP <- NULL
allRes_MF <- NULL
allRes_CC <- NULL
genesOfInterest <- NULL

geneUniverse <- names(geneID2GO)
genesOfInterest <- as.character(interesting_genes)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

node = 5
algorithm = "classic"

####BP

myGOdataBP <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize = node)

# run the Fisher's exact tests
#test <- runTest(myGOdataBP, algorithm="parentchild", statistic="fisher")
test <- runTest(myGOdataBP,  algorithm = algorithm, statistic="fisher")

mysummary <- summary(attributes(test)$score <= 0.05)

numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.05

if(numsignif > 0){
allRes_BP <- GenTable(myGOdataBP, pvalue=test, topNodes = numsignif)
}else{allRes_BP <- NULL}


#####print into files
out_BP <- paste0(comp, "_BP_table")

write.table(allRes_BP, file = out_BP ,quote=F,row.names = F, sep = "\\t")

if(length(allRes_BP$Term) > 1){
  for(i in 1:length(allRes_BP$pvalue)){
    if(allRes_BP$pvalue[i] == "< 1e-30"){
      allRes_BP$pvalue[i] <- 1e-30 
    }
  }
}  

if(length(allRes_BP$Term) > 20){
  allRes_BP <- allRes_BP[1:20,]
}
  

tag_BP <- paste0(comp, "_tagcloud.pdf")

pdf(file = tag_BP,width=30, height = 5)
colors <- colorRampPalette( brewer.pal( 12, "Paired" ) )( length(allRes_BP$Term) )
tagcloud(strmultline(allRes_BP$Term), weights = -log(as.numeric(allRes_BP$pvalue)),col=colors)
dev.off()

