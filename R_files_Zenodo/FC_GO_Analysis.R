# The genetic basis of variation in sexual aggression: evolution versus social plasticity

# FC RNAseq Gene Ontology Analysis
# Andrew M. Scott

# https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

#### NOTE before running:
#### Script depends on objects generated from running FC_Figures, FC_Vector_Corr, FC_Plasticity_DE and FC_Artificial_Selection_DE Scripts!

library(topGO)
library(Rgraphviz)

gene.assoc <- read_csv("Gene_Ontology/gene_association.csv")
go.list <- read_csv("Gene_Ontology/go-basic.csv")

##### Sig Gene Lists

# Plasticity
plast.go.hits <- plast.hits.vc # 375 hits
plast.go.hits <- plast.go.hits[order(row.names(plast.go.hits)), ]

# Artificial selection
AS.trt.MD.hits
apeglm.AS.hits_unique <- apeglm.AS.hits[!(row.names(apeglm.AS.hits) %in% row.names(AS.trt.MD.hits)), ]
apeglm.AS.hits_unique$gene <- row.names(apeglm.AS.hits_unique)
apeglm.AS.hits_unique <- apeglm.AS.hits_unique[ ,c(6,7,9)]

AS.go.hits <- AS.trt.MD.hits[ ,c(1:3)]
AS.go.hits <- rbind(AS.go.hits, apeglm.AS.hits_unique) # 1030 hits, including unique hits from apeglm 

# Overlap
overlap.go.hits <- overlapping.hits.expr

### Make csvs for using with WebGestalt
#write.csv(plast.go.hits, file = "plastgohits.csv")
#write.csv(AS.go.hits, file = "asgohits.csv")
#write.csv(overlap.go.hits, file = "overlapgohits.csv")

# annotations (gene iD -> GO ID)
FBgn.2.GO.long <- gene.assoc[, c(1, 4)]
FBgn.2.GO <- vector(mode = "list", length = length(unique(FBgn.2.GO.long$FBgn)))
names(FBgn.2.GO) <- unique(FBgn.2.GO.long$FBgn)

for (g in names(FBgn.2.GO)) {
  v <- FBgn.2.GO.long[FBgn.2.GO.long$FBgn == g, ]
  FBgn.2.GO[[g]] <- v$GO
}

str(head(FBgn.2.GO))


##### Plasticity GO analysis

# list of gene identifiers
plast.go.geneList <- factor(as.integer(row.names(plast.trt.MD) %in% row.names(plast.go.hits)))

# create vector of sig (1), nonsig (0) genes in full plasticity gene list
names(plast.go.geneList) <- row.names(plast.trt.MD)

# create topGOdata object
GO.data.plast <- new("topGOdata", ontology = "MF", allGenes = plast.go.geneList, 
                     annot = annFUN.gene2GO, gene2GO = FBgn.2.GO, nodeSize = 5) # minimum 5 annotated genes per GO term

# Fisher's exact test to test enrichment
test.stat.Fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher Test")
plast.resultFisher <- getSigGroups(GO.data.plast, test.stat = test.stat.Fisher)

plast.res.table <- GenTable(GO.data.plast, classic = plast.resultFisher, orderBy = "classic", ranksOf = "classic", topNodes = 35)

# write.csv(plast.res.table, file = "Plasticity_GO.csv")

### FOR JAN 2022 REVISION - FDR adjusted p-values

allGO.plast <- usedGO(GO.data.plast)
plast.res.FULL.table <- GenTable(GO.data.plast, classic = plast.resultFisher, orderBy = "classic", ranksOf = "classic", topNodes = length(allGO.plast))
plast.res.FULL.table$fdr.p <- p.adjust(plast.res.FULL.table$classic, method = "fdr")

# write.csv(plast.res.FULL.table, file = "Plasticity_GO_withFDR.csv")

# GO graph
par(cex = 0.85)
showSigOfNodes(GO.data.plast, score(plast.resultFisher), firstSigNodes = 20, useInfo = 'all')


##### AS GO analysis

# list of gene identifiers
AS.go.geneList <- factor(as.integer(row.names(AS.trt.MD) %in% row.names(AS.go.hits)))

# create vector of sig (1), nonsig (0) genes in full plasticity gene list
names(AS.go.geneList) <- row.names(AS.trt.MD)

# create topGOdata object
GO.data.AS <- new("topGOdata", ontology = "MF", allGenes = AS.go.geneList, 
                  annot = annFUN.gene2GO, gene2GO = FBgn.2.GO, nodeSize = 5)

# Fisher's exact test to test enrichment
test.stat.Fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher Test")
AS.resultFisher <- getSigGroups(GO.data.AS, test.stat = test.stat.Fisher)

AS.res.table <- GenTable(GO.data.AS, classic = AS.resultFisher, orderBy = "classic", ranksOf = "classic", topNodes = 14)

# GO graph
par(cex = 0.85)
showSigOfNodes(GO.data.AS, score(AS.resultFisher), firstSigNodes = 14, useInfo = 'all')

# write.csv(AS.res.table, file = "AS_GO.csv")

### FOR JAN 2022 REVISION - FDR adjusted p-values

allGO.AS <- usedGO(GO.data.AS)
AS.res.FULL.table <- GenTable(GO.data.AS, classic = AS.resultFisher, orderBy = "classic", ranksOf = "classic", topNodes = length(allGO.AS))
AS.res.FULL.table$fdr.p <- p.adjust(AS.res.FULL.table$classic, method = "fdr")

# write.csv(AS.res.FULL.table, file = "AS_GO_withFDR.csv")


##### Overlap GO analysis (using AS list as total list for sampling)

# list of gene identifiers
overlap.go.geneList <- factor(as.integer(row.names(AS.trt.MD) %in% row.names(overlap.go.hits)))

# create vector of sig (1), nonsig (0) genes in full plasticity gene list
names(overlap.go.geneList) <- row.names(AS.trt.MD)

# create topGOdata object
GO.data.overlap <- new("topGOdata", ontology = "MF", allGenes = overlap.go.geneList, 
                       annot = annFUN.gene2GO, gene2GO = FBgn.2.GO, nodeSize = 5)

# Fisher's exact test to test enrichment
test.stat.Fisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher Test")
overlap.resultFisher <- getSigGroups(GO.data.overlap, test.stat = test.stat.Fisher)

overlap.res.table <- GenTable(GO.data.overlap, classic = overlap.resultFisher, orderBy = "classic", ranksOf = "classic", topNodes = 5)

# GO graph
par(cex = 0.85)
showSigOfNodes(GO.data.overlap, score(overlap.resultFisher), firstSigNodes = 5, useInfo = 'all')

# write.csv(overlap.res.table, file = "Overlap_GO.csv")

### FOR JAN 2022 REVISION - FDR adjusted p-values

allGO.overlap <- usedGO(GO.data.overlap)
overlap.res.FULL.table <- GenTable(GO.data.overlap, classic = overlap.resultFisher, orderBy = "classic", ranksOf = "classic", topNodes = length(allGO.overlap))
overlap.res.FULL.table$fdr.p <- p.adjust(overlap.res.FULL.table$classic, method = "fdr")

#write.csv(overlap.res.FULL.table, file = "overlap_GO_withFDR.csv")