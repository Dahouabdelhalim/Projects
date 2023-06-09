# External packages required to run the code below

library(Biobase)
library(GOexpress)
library(biomaRt)
library(ggplot2)
library(RColorBrewer)
library(tools)

# Move to the folder containing the "MDM.eSet.rds" file
#setwd('./path_to_file')

# Load the ExpressionSet in the appropriate variable name for the code below
MDM.eSet = readRDS(file = 'MDM.eSet.rds')

# Download biomart annotations from Ensembl 77 (for microarray!) ----------

# View the available datasets in the Ensembl release 77 archive
listMarts(host = 'oct2014.archive.ensembl.org')

# Prepare the connection to the Bos taurus Ensembl annotation dataset
ensembl77 = useMart(
  host = 'oct2014.archive.ensembl.org',
  biomart = 'ENSEMBL_MART_ENSEMBL',dataset = 'btaurus_gene_ensembl')

# View the available gene annotations
listAttributes(mart = ensembl77, page = 'feature_page')

# Download the gene annotations
allgenes = getBM(
  attributes = c('affy_bovine', 'external_gene_name', 'description'),
  mart = ensembl77)
# Rename the microarray-specific header to the header supported by GOexpress
colnames(allgenes)[1] = 'gene_id'

# Download the gene ontology annotations
allGO = getBM(
  attributes = c('go_id', 'name_1006', 'namespace_1003'),
  mart = ensembl77)

# Download the table mapping genes to gene ontologies
GO_genes = getBM(attributes = c('affy_bovine', 'go_id'), mart = ensembl77)
# Rename the microarray-specific header to the header supported by GOexpress
colnames(GO_genes)[1] = 'gene_id'
# Extra step for sanity check:
# some gene ontologies have no annotated genes
# some genes have no annotated gene ontologies
# Remove those blank fields
sum(GO_genes$go_id == '') # number of empty GO identifiers in the query
GO_genes = GO_genes[GO_genes$go_id != '',]
sum(GO_genes$gene_id == '') # number of empty GO identifiers in the query
GO_genes = GO_genes[GO_genes$gene_id != '',]

# Save the annotations in local files
save(allgenes, file='allgenes.rdata')
save(allGO, file='allGO.rdata')
save(GO_genes, file='GO_genes.rdata')

# GOexpress (Infection) ---------------------------------------------------

# Set the random seed to allow reproducible results
set.seed(4598)

# Run the GO_analyse function with desired parameters and local annotations
GOx.Infection = GO_analyse(
  eSet = MDM.eSet, f = 'Infection',
  subset = list(Time=c('2HR','6HR','24HR')),
  GO_genes = GO_genes, all_GO = allGO, all_genes = allgenes)

# Save the result of the analysis to a R session file
save(GOx.Infection, file='GOx.Infection.rda')
# Compress the file to save space
resaveRdaFiles('GOx.Infection.rda')

# Look at the top-ranked genes prior to filtering
head(GOx.Infection$genes, n=20)
# Export the entire unfiltered gene scoring table to a TAB-delimited file
write.table(
  x = GOx.Infection$genes, file = 'GOx.Infection.genes.txt',
  sep = '\\t')

# Expression plot and profile ---------------------------------------------

# Expressin profile of some top-ranking genes (prior to filtering)
expression_profiles_symbol(
  gene_symbol = 'CCL5', result = GOx.Infection, eSet = MDM.eSet,
  x_var = 'Hours.post.infection', seriesF = 'Animal.Infection', line.size = 3,
  ylab = expression('log'[2]*' intensity'),
  xlab = 'Hours post-infection')
expression_plot_symbol(
  gene_symbol = 'CCL5', result = GOx.Infection, eSet = MDM.eSet,
  x_var = 'Hours post-infection',
  ylab = expression('log'[2]*' intensity'))

# Pvalue ------------------------------------------------------------------

# Compute and append p-values to the gene ontology scoring table
# Duration ~ 30 min
GOx.Infection.pval = pValue_GO(
  GOx.Infection, N = 1000,
  ranked.by = 'Rank', rank.by = 'p.val')

# Export the entire unfiltered table to a TAB-delimited file
write.table(
  x = GOx.Infection.pval$GO, file = 'GOx.Infection.pValue.GO.txt',
  sep = '\\t', row.names = F)

# Filtering of ontologies -------------------------------------------------

# Filtering for all ontologies with >= 15 genes and P-value <= 0.05
GOx.infection.filtered = subset_scores(
  result = GOx.Infection.pval, total=15, p.val = 0.05)

head(GOx.infection.filtered$GO, n=20)

write.table(
  x = GOx.infection.filtered$GO,
  file = 'GOx.Infection.filtered.GO.txt', sep = '\\t', row.names = F)

# List the top genes in the filtered data 
# Genes without annotations were discarded, leaving gaps in the apparent ranks
head(GOx.infection.filtered$genes, n=20)

# Heatmap chemokines ------------------------------------------------------

# Expression heatmap of genes annotated with chemokine activity
heatmap_GO(
  go_id = 'GO:0008009', result = GOx.infection.filtered, eSet = MDM.eSet,
  cexRow = 0.9, labRow = "Infection.Time",
  margins = c(12, 7))

# Table_genes for chemokine activity --------------------------------------

write.table(
  x = table_genes(
    go_id = 'GO:0008009', result = GOx.infection.filtered, data.only = TRUE),
  file = 'table-genes_chemokines_dataset.txt', sep = '\\t')
