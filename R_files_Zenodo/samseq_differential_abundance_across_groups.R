#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# Assess differential abundance across social groups using nonparametric samseq method

rm(list=ls())
graphics.off()

library("phyloseq")
library("ggplot2")
library("plyr")
library("RColorBrewer")
library("samr")

setwd("~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/")
getwd()

# input: 

# pre-processed phyloseq object with raw reads (with infrequent otus taken out)
load("data/sifaka_trimmed_phyloseq_raw_openref.RData") #sifaka0 (same as kmnp_data)
sample_data(sifaka0)

# pre-processed phyloseq object with normalized reads
load("data/sifaka_trimmed_phyloseq_normalized_openref.RData") #kmnp_trimmed

# Differential abundance of genus OTUS across social groups
# 	- samseq method
#		- can use raw reads (sifaka0 object)

# output: 
# - significant bacteria genera: sam1 object
#######################################################################
## Differential Abundance of Phyla OTUS
#######################################################################
# Convert the otu and taxa tables to data frame
summary(tax_table(sifaka0))
otum = data.frame(OTUID = rownames(otu_table(sifaka0)),otu_table(sifaka0))
head(otum)
taxm = data.frame(OTUID = rownames(tax_table(sifaka0)),tax_table(sifaka0)) 
taxm <- as.matrix(taxm)
taxm <- data.frame(taxm)

# Merge the otu table and genus column
taxOTU = merge(otum,taxm,by="OTUID")
head(taxOTU)
str(taxOTU)
taxOTUphyla = taxOTU[,-c(1,49,51:55)] #remove other ranks than phyla
tagg = aggregate(.~Phylum,taxOTUphyla,FUN = sum) #bacteria phylum x sifaka individual contingency table
tax.family = tagg[,1] #bacteria phyla names
tagg1 = tagg[,-1] #otu table according to phyla
group = sample_data(sifaka0)$Group

gensums = apply(tagg1,1,sum) #sum up counts according to bacteria phyla
gen.100 = as.character(tax.family[gensums>50]) #bacteria phyla with over 50 counts
gensagg.100 = cbind(tagg1,gensums)
gensagg.100 = subset(gensagg.100, gensums>50) 
gensagg.100 <- gensagg.100[,-48]

sam1 = SAMseq(gensagg.100,group,resp.type="Multiclass",genename=gen.100, nperms=1000, nresamp= 100, geneid = gen.100,random.seed=6, fdr.output=0.05)
siggenes = data.frame(sam1$siggenes.table$genes.up)
siggenes = siggenes[,-2]
names(siggenes) = c("Phylum","Score","I","II","III","IV", "V","VI","Camp","q")
siggenes
siggenes$adj.p <- as.numeric(as.character(siggenes$q))
siggenes$adj.p <- siggenes$adj.p / 100
siggenes
siggenes2 <- subset(siggenes, adj.p < 0.05)
siggenes2

#######################################################################
## Differential Abundance of Family OTUS
#######################################################################
# Convert the otu and taxa tables to data frame
summary(tax_table(sifaka0))
otum = data.frame(OTUID = rownames(otu_table(sifaka0)),otu_table(sifaka0))
head(otum)
taxm = data.frame(OTUID = rownames(tax_table(sifaka0)),tax_table(sifaka0)) 
taxm <- as.matrix(taxm)
taxm[is.na(taxm)] <- "Unclassified"
taxm <- data.frame(taxm)

# Merge the otu table and genus column
phyla_order = paste(taxm$Phylum, sep="-", taxm$Order)
P_O_F = paste(phyla_order, sep="-", taxm$Family)
taxOTU = merge(otum,taxm,by="OTUID")
taxOTU = cbind(taxOTU,P_O_F)

## Fix this so you grab it
taxOTUfamily = taxOTU[,-c(1,49:55)] #remove other ranks
tagg = aggregate(.~P_O_F,taxOTUfamily,FUN = sum) #bacteria genus x sifaka individual contingency table
tagg$P_O_F
tagg = tagg[-c(27,35,60,62,93,100),] #remove unclassified families
tax.family = tagg[,1] #bacteria family names
tagg1 = tagg[,-1] #otu table according to family
group = sample_data(sifaka0)$Group

rownames(tagg1) <- tax.family

gensums = apply(tagg1,1,sum) #sum up counts according to bacteria family
gen.100 = as.character(tax.family[gensums>100]) #bacteria families with over 100 counts
gensagg.100 = cbind(tagg1,gensums)
gensagg.100 = subset(gensagg.100, gensums>100) 
gensagg.100 = gensagg.100[,-48]

sam1 = SAMseq(gensagg.100,group,resp.type="Multiclass",genename=gen.100, nperms=1000, nresamp= 100, geneid = gen.100,random.seed=6, fdr.output=0.05)
str(sam1)
siggenes = data.frame(sam1$siggenes.table$genes.up)
siggenes = siggenes[,-2]
names(siggenes) = c("Family","Score","I","II","III","IV", "V","VI","Camp","q")
siggenes$adj.p <- as.numeric(as.character(siggenes$q))
siggenes$adj.p <- siggenes$adj.p / 100
siggenes2 <- subset(siggenes, adj.p <= 0.05)
siggenes2

write.csv(siggenes2, file="data/diff_bacteria_families.csv", row.names=FALSE)

save(sam1,file="data/siggenes_families_sifaka_openref.RData")
load("data/siggenes_families_sifaka_openref.RData")

#######################################################################
## Differential Abundance of Genus OTUS
#######################################################################
# just genus
# Convert the otu and taxa tables to data frame
otum = data.frame(OTUID = rownames(otu_table(sifaka0)),otu_table(sifaka0))
taxm = data.frame(OTUID = rownames(tax_table(sifaka0)),tax_table(sifaka0)) 
taxm <- as.matrix(taxm)
taxm[is.na(taxm)] <- "Unclassified"
taxm <- data.frame(taxm)

# select the genus column
genus = taxm$Genus
head(genus)
head(taxm)

taxOTU = merge(otum,taxm,by="OTUID")
head(taxOTU)
str(taxOTU)

taxOTUgenus = taxOTU[,-c(1,49:53,55)] #remove other ranks
tagg = aggregate(.~Genus,taxOTUgenus,FUN = sum) #bacteria genus x sifaka individual contingency table

tax.genus = tagg[,1] #bacteria genus names
tagg1 = tagg[,-1] #otu table according to genus
colsums = colSums(tagg1)
group = sample_data(sifaka0)$Group

rownames(tagg1) <- tax.genus

gensums = apply(tagg1,1,sum) #sum up counts according to bacteria genus
gen.100 = as.character(tax.genus[gensums>100]) #bacteria genera with over 100 counts
gensagg.100 = cbind(tagg1,gensums)
gensagg.100 = subset(gensagg.100, gensums>100) 
rownames(gensagg.100)
colnames(gensagg.100)
gensagg.100 <- gensagg.100[-50,-48]
gen.100 <- gen.100[-50]

sam1 = SAMseq(gensagg.100,group,resp.type="Multiclass",genename=gen.100, nperms=1000, nresamp= 100, geneid = gen.100,random.seed=6)
siggenes = data.frame(sam1$siggenes.table$genes.up)
siggenes = siggenes[,-2]
names(siggenes) = c("Genus","Score","I","II","III","IV", "V","VI","Camp","q")
siggenes$adj.p <- as.numeric(as.character(siggenes$q))
siggenes$adj.p <- siggenes$adj.p / 100
siggenes2 <- subset(siggenes, adj.p < 0.05)
siggenes2

save(sam1,file="data/siggenes_genera_sifaka_openref.RData")
load("data/siggenes_genera_sifaka_openref.RData")

write.csv(siggenes2, file="data/diff_bacteria_genera.csv", row.names=FALSE)