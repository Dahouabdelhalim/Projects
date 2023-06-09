if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ALDEx2")
library(ALDEx2)

#Carbon_Nitrogen_Phosphorus_Sulfur metabolism pathways from KEGG
CNPS_metabolism <- read.table("CNPS_metabolism.txt")

#Predicted KO functional abundances from PICRUSt2
pred_metagenome_unstrat <- read.delim("pred_metagenome_unstrat.tsv")
#Predicted MetaCyc Pathway functional abundances from PICRUSt2
path_abun_unstrat <- read.delim("path_abun_unstrat.tsv")

### KOs 
metagenome <- pred_metagenome_unstrat
rownames(metagenome) <- metagenome$X.function
metagenome <- metagenome[-c(1)]
conds <- c(rep("Sub", 3), rep('Int', 3), rep("Sec", 3))
ko <- aldex.clr(round(metagenome), conds, mc.samples=333, denom="all", verbose=F)
ko.kw <- aldex.kw(ko)
sig_kos <- subset(ko.kw, (ko.kw$glm.ep <= 0.05))
write.csv(sig_kos, "aldex_significant_all_KO.csv")

### Metabolism associated KOs 
metagenome <- pred_metagenome_unstrat

#note we are subsetting the set of all KOs to only be those associated with Carbon_Nitrogen_Phosphorus_Sulfur metabolism
metagenome <-merge(metagenome, CNPS_metabolism, by.x='X.function', by.y='V1')
rownames(metagenome) <- metagenome$X.function
metagenome <- metagenome[-c(1)]
conds <- c(rep("Sub", 3), rep('Int', 3), rep("Sec", 3))
ko <- aldex.clr(round(metagenome), conds, mc.samples=333, denom="all", verbose=F)
ko.kw <- aldex.kw(ko)
sig_cnps <- subset(ko.kw, (ko.kw$glm.ep <= 0.05))
write.csv(sig_cnps, "aldex_significant_CNPS_KO.csv")

### MetaCyc pathways
pathway <- path_abun_unstrat
rownames(pathway) <- pathway$X.pathway
pathway <- pathway[-c(1)]
conds <- c(rep("Sub", 3), rep('Int', 3), rep("Sec", 3))
x <- aldex.clr(round(pathway), conds, mc.samples=333, denom="all", verbose=F)
x.kw <- aldex.kw(x)
sig_pathway <- subset(x.kw, x.kw$glm.ep <= 0.05)
write.csv(sig_pathway, "aldex_significant_pathway.csv")
