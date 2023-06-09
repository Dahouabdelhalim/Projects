# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

rm(list=ls())
graphics.off()

library(phyloseq)  #bioconductor 
library(vegan)
library(DESeq2) #bioconductor 
library(samr)
library(dplyr)
library(rstatix)
library(rcompanion)
library(stringr)

library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
setwd(dir)

## import phyloseq objects

## phyloseq object (kmnp_trimmed)
# normalized sequence counts, limited to ASVs shared across at least two samples 
## marked individuals only with sample CST assignments
load("Rdata/sifaka_trimmed_phyloseq_normalized_dada2_CSTs.RData")
head(sample_data(kmnp_trimmed))

#phyloseq object (kmnp2) 
# raw ASV counts including singletons; required for samr analyses
# marked individuals only
load("Rdata/sifaka_allotus_phyloseq_openref_2016.RData")
head(sample_data(kmnp2))

## merge so that samples with raw ASV counts have CST assignments
kmnp_cst=kmnp_trimmed
cst_df = sample_data(kmnp_cst) %>% as_tibble() %>% dplyr::select(X.SampleID,Group,Year,CST) %>%arrange(Group,CST)

cst_ids = cst_df$X.SampleID

reduced_kmnp2 = kmnp2 %>% subset_samples(X.SampleID %in% cst_ids)
new_df = left_join(as_tibble(sample_data(reduced_kmnp2)),cst_df %>% dplyr::select(X.SampleID,CST),by="X.SampleID") 
rownames(new_df)<- new_df$X.SampleID
sample_data(reduced_kmnp2) = data.frame(new_df)

#######################################################################
## Figure S13: Differential Abundance of Phyla across CSTs
#######################################################################
reduced_kmnp2

TT = tax_table(reduced_kmnp2)
TT[is.na(TT)] <- "Unclassified"
tax_table(reduced_kmnp2) <- TT

otum = data.frame(OTUID = rownames(t(otu_table(reduced_kmnp2))),t(otu_table(reduced_kmnp2)))
taxm = data.frame(OTUID = rownames(tax_table(reduced_kmnp2)),tax_table(reduced_kmnp2))
taxm <- as.matrix(taxm)
taxm <- data.frame(taxm)

# Merge the otu table and genus column
taxOTU = merge(otum,taxm,by="OTUID")
#remove other ranks than phyla
taxOTUphyla = taxOTU[,!(colnames(taxOTU)%in%c("OTUID","Kingdom","Class","Order","Family","Genus","Species"))]
tagg = aggregate(.~Phylum,taxOTUphyla,FUN = sum) #bacteria phyla x sifaka individual contingency table
tax.family = tagg[,1] #bacteria phyla names
tagg1 = tagg[,-1] #asv table according to phyla
CST = sample_data(reduced_kmnp2)$CST

gensums = rowSums(tagg1)
names(gensums) = tax.family
total = sum(gensums)
rel_abund = gensums/total
sort(rel_abund,decreasing = T)
sort(rel_abund[rel_abund>0.01],decreasing = T)

gen.100 = as.character(names(rel_abund[rel_abund>0.01])) #bacteria genera with > 1% relative abundance
gensagg.100 = cbind(tagg1,gensums,rel_abund)
gensagg.100 = subset(gensagg.100, rel_abund>0.01) 

rownames(gensagg.100)
colnames(gensagg.100)
length(colnames(gensagg.100))
colnames(gensagg.100)

gensagg.100 = gensagg.100[,!(colnames(gensagg.100)%in% c("gensums","rel_abund"))]
sam1 = SAMseq(gensagg.100,CST,resp.type="Multiclass",genename=gen.100, nperms=1000, nresamp= 100, geneid = gen.100,random.seed=6, fdr.output=0.05)
siggenes = data.frame(sam1$siggenes.table$genes.up)
siggenes = siggenes[,-2]
siggenes = siggenes %>% filter(q.value...<=0.05)

phylum_long = bind_cols(gensagg.100,Phylum=gen.100)%>% tidyr::pivot_longer(cols=AB.15:ZN7)

phylum_long = left_join(phylum_long,data.frame(sample_data(reduced_kmnp2)),by=c("name"="X.SampleID"))

# annotation table with adjusted pvals and y-position of the labels
reduced=phylum_long %>% filter(Phylum %in% siggenes$Gene.ID)
# Pairwise comparisons
pwc <- reduced %>% 
  group_by(Phylum) %>%
  rstatix::dunn_test(value ~ CST, p.adjust.method = "BH") 
pwc
pwc$group_name = paste(pwc$group1,pwc$group2,sep="-")

pwc_group = group_by(pwc,Phylum)

cld_fun <- function(sub_df){
  ret_df = cldList(p.adj ~ group_name,data = sub_df,threshold = 0.05)
  return(ret_df)
}
dunn_df = pwc %>%
  group_by(Phylum)%>%
  do(as.data.frame(cld_fun(.)))%>%
  ungroup()
names(dunn_df)[2]<-"CST"

max_values = reduced %>%
  group_by(Phylum)%>%
  dplyr::summarize(max_y = max(value)+200)

dunn_df = left_join(dunn_df,max_values,by="Phylum")
brewer = brewer.pal(n = 8, name = "RdBu")
reduced$CST = as.factor(reduced$CST)
reduced$CST  = factor(reduced$CST ,levels=c("1","2","3","4","5","6","7"))

length(unique(reduced$name))
phylum_abund = ggplot(reduced,aes(x=as.factor(CST),y=value,fill=CST))+
  geom_boxplot()+
  facet_wrap(~Phylum,scales="free_y")+
  geom_text(data = dunn_df, aes(y = max_y, label = Letter), 
            position = position_dodge(width = .75),color="black")+
  scale_fill_manual(values = c(brewer[1],brewer[2],brewer[3],brewer[4],brewer[5],brewer[6],brewer[7]),
                    breaks=c("1","2","3","4","5","6","7"),
                    labels=c("1","2","3","4","5","6","7"))+
  ylab(expression(paste(log[10], " (","Phylotype abundance",")")))+
  xlab("CST")+
  theme_bw(base_size = 14)+
  theme(legend.position = "none")+
  scale_y_log10()
phylum_abund
save_plot(phylum_abund,filename = "figures/SFig13_enriched_phyla_sifaka_2012-2016.pdf",base_width = 12,base_height = 6)
#######################################################################
## Figure S14: Differential Abundance of Genera across CSTs
#######################################################################
# Convert the otu and taxa tables to data frame
red <- subset_taxa(reduced_kmnp2, Genus!="Unclassified")
red

summary(tax_table(red))
otum = data.frame(OTUID = rownames(t(otu_table(red))),t(otu_table(red)))
head(otum)
taxm = data.frame(OTUID = rownames(tax_table(red)),tax_table(red))
taxm <- as.matrix(taxm)
taxm <- data.frame(taxm)

P_F_G = paste(taxm$Phylum,taxm$Family,taxm$Genus,sep="-")

taxOTU = merge(otum,taxm,by="OTUID")
taxOTU = cbind(taxOTU,P_F_G)
taxOTUgenus= taxOTU[,!(colnames(taxOTU)%in%c("OTUID","Kingdom","Phylum", "Class","Order","Family","Genus","Species"))]

tagg = aggregate(.~P_F_G,taxOTUgenus,FUN = sum) #bacteria genus x sifaka individual contingency table
tax.genus = tagg[,1] #bacteria genus names
tagg1 = tagg[,-1] 
CST = sample_data(red)$CST

rownames(tagg1) <- tax.genus

gensums = apply(tagg1,1,sum) #sum up counts according to bacteria genus
total = sum(gensums)
rel_abund = gensums/total
sort(rel_abund[rel_abund>0.01])

gen.100 = as.character(names(rel_abund[rel_abund>0.01])) #bacteria genera with > 1% relative abundance
gensagg.100 = cbind(tagg1,gensums,rel_abund)
gensagg.100 = subset(gensagg.100, rel_abund>0.01) 

gensagg.100 = gensagg.100[,!(colnames(gensagg.100)%in% c("gensums","rel_abund"))]
sam1 = SAMseq(gensagg.100,CST,resp.type="Multiclass",genename=gen.100, nperms=1000, nresamp= 100, geneid = gen.100,random.seed=6, fdr.output=0.05)
siggenes = data.frame(sam1$siggenes.table$genes.up)
siggenes = siggenes[,-2]
names(siggenes)[10]<-"q"
head(siggenes)
siggenes2 <- subset(siggenes, q <= 0.05)
siggenes2$Score.d. = as.numeric(siggenes2$Score.d.)
siggenes2 %>% arrange(desc(Score.d.))

family_long = bind_cols(gensagg.100,family=gen.100)%>% tidyr::pivot_longer(cols=AB.15:ZN7)

family_long = left_join(family_long,data.frame(sample_data(reduced_kmnp2)),by=c("name"="X.SampleID"))

# annotation table with adjusted pvals and y-position of the labels
reduced=family_long %>% filter(family %in% siggenes2$Gene.ID)

# Pairwise comparisons
pwc <- reduced %>% 
  group_by(family) %>%
  rstatix::dunn_test(value ~ CST, p.adjust.method = "BH") 
pwc
pwc$group_name = paste(pwc$group1,pwc$group2,sep="-")
unique(pwc$family)

split_df = split(pwc, pwc$family)
signif_list = lapply(split_df, function(x) filter(x, p.adj.signif != "ns"))
signif_list = signif_list[sapply(signif_list, nrow)>0]
family_keep = names(signif_list)
length(family_keep)

cld_fun <- function(sub_df){
  ret_df = cldList(p.adj ~ group_name,data = sub_df,threshold = 0.05)
  return(ret_df)
}
dunn_df = pwc %>%
  filter(family %in% family_keep)%>%
  group_by(family)%>%
  do(as.data.frame(cld_fun(.)))%>%
  ungroup()
names(dunn_df)[2]<-"CST"

max_values = reduced %>%
  group_by(family)%>%
  dplyr::summarize(max_y = max(value)+50)

dunn_df = left_join(dunn_df,max_values,by="family")

brewer = brewer.pal(n = 8, name = "RdBu")
reduced$CST = as.factor(reduced$CST)
reduced$CST  = factor(reduced$CST ,levels=c("1","2","3","4","5","6","7"))

name_mat = str_split_fixed(reduced$family, "-", n = 4)
phy_fam = paste(paste0("P.",name_mat[,1]),paste0("F.",name_mat[,2]),sep=" ")
genus = paste(paste0("G.",name_mat[,3]),name_mat[,4],sep=" ")
genus = gsub("_"," ",genus)
genus = str_trim(genus,"right")
whole_name = paste0(phy_fam,"\\n",genus)
reduced$label = whole_name


name_mat = str_split_fixed(dunn_df$family, "-", n = 4)
phy_fam = paste(paste0("P.",name_mat[,1]),paste0("F.",name_mat[,2]),sep=" ")
genus = paste(paste0("G.",name_mat[,3]),name_mat[,4],sep=" ")
genus = gsub("_"," ",genus)
genus = str_trim(genus,"right")
whole_name = paste0(phy_fam,"\\n",genus)
dunn_df$label = whole_name


genus_abund = ggplot(reduced,aes(x=as.factor(CST),y=value,fill=CST))+
  geom_boxplot()+
  facet_wrap(~label,scales="free_y")+
  geom_text(data = dunn_df, aes(y = max_y, label = Letter), 
            position = position_dodge(width = .75),color="black")+
  scale_fill_manual(values = c(brewer[1],brewer[2],brewer[3],brewer[4],brewer[5],brewer[6],brewer[7]),
                    breaks=c("1","2","3","4","5","6","7"),
                    labels=c("1","2","3","4","5","6","7"))+
  ylab(expression(paste(log[10], " (","Phylotype abundance",")")))+
  xlab("CST")+
  theme_bw(base_size = 14)+
  theme(legend.position = "none",strip.text.x = element_text(size = 8))+
  scale_y_log10()
genus_abund
save_plot(genus_abund,filename = "figures/SFig14_enriched_genera_sifaka_2012-2016.pdf",base_width = 16,base_height = 10)
