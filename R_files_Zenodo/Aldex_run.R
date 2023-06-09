if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")

library(ALDEx2)
setwd('C:/Gresham/Project_Gravimondo/BGS_final/Aldex_results_2/')

#Carbon_Nitrogen_Phosphorus_Sulfur metabolism pathways from KEGG
CNPS_metabolism <- read.table("C:/Gresham/Project_Gravimondo/BGS_final/sigilo_results/_ko_list.log")
#Predicted KO functional abundances from PICRUSt2

pred_metagenome_unstrat <- read.delim("C:/Gresham/Project_Gravimondo/BGS_final/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv")
#Predicted MetaCyc Pathway functional abundances from PICRUSt2
#path_abun_unstrat <- read.delim("path_abun_unstrat.tsv")

### KOs 
metagenome <- pred_metagenome_unstrat
rownames(metagenome) <- metagenome$function.
metagenome <- metagenome[-c(1)]

#CT1	CT2	CT3	G1	G2	G3	OL1	OL2	OL3	
conds <- c(rep("sub", 3), rep('in', 3), rep("sup", 3))
ko <- aldex.clr(round(metagenome), conds, mc.samples=333, denom="all", verbose=F)
ko.kw <- aldex.kw(ko)

### Metabolism associated KOs 
metagenome <- pred_metagenome_unstrat
#note we are subsetting the set of all KOs to only be those associated with Carbon_Nitrogen_Phosphorus_Sulfur metabolism
metagenome <-merge(metagenome, CNPS_metabolism, by.x='function.', by.y='V1')
rownames(metagenome) <- metagenome$function.
metagenome <- metagenome[-c(1)]

#CT1	CT2	CT3	G1	G2	G3	OL1	OL2	OL3	
conds <- c(rep("sub", 3), rep('int', 3), rep("sup", 3))
ko <- aldex.clr(round(metagenome), conds, mc.samples=333, denom="all", verbose=F)
ko.kw <- aldex.kw(ko)
sig_cnps <- subset(ko.kw, (ko.kw$glm.ep <= 0.05))
write.csv(sig_cnps, "aldex_significant_all_glm.ep_KO.csv")
sig_cnps <- subset(ko.kw, (ko.kw$kw.ep <= 0.05))
write.csv(sig_cnps, "aldex_significant_all_kw.ep_KO.csv")
sig_cnps <- subset(ko.kw, (ko.kw$glm.BH <= 0.05))
write.csv(sig_cnps, "aldex_significant_all_glm.BH_KO.csv")
sig_cnps <- subset(ko.kw, (ko.kw$kw.BH <= 0.05))
write.csv(sig_cnps, "aldex_significant_all_kw.BH_KO.csv")


sig_cnps$row_names <- row.names(sig_cnps)

#sig_cnps <- subset(ko.kw, (ko.kw$kw.ep <= 0.05))
sig_cnps <- subset(ko.kw, (ko.kw$glm.ep <= 0.05))
write.csv(sig_cnps, "aldex_significant_CT_G_OL_CNPS_KO.csv")

### MetaCyc pathways
#pathway <- path_abun_unstrat
#rownames(pathway) <- pathway$X.pathway
#pathway <- pathway[-c(1)]

#conds <- c(rep("Sub", 3), rep('Int', 3), rep("Sec", 3))
#x <- aldex.clr(round(pathway), conds, mc.samples=333, denom="all", verbose=F)
#x.kw <- aldex.kw(x)

#sig_pathway <- subset(x.kw, x.kw$glm.ep <= 0.05)
#write.csv(sig_pathway, "aldex_significant_pathway.csv")