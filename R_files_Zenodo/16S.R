library(mctoolsr) # https://github.com/leffj/mctoolsr
library(vegan)
library(tidyverse)
library(nlme)
library(patchwork)

############## Load data
##############
##############

otu_table_fp <- "seqtab_wTax_v2.txt"
mapping_fp <- "sample_metadata_v2.txt"
all <- load_taxa_table(otu_table_fp, mapping_fp)

############## Misc initial filtering/processing/tabulating
##############
##############

#replace sequence IDs with sample IDs
colnames(all$data_loaded) <- all$map_loaded$SampleID
rownames(all$map_loaded) <- all$map_loaded$SampleID

all$map_loaded$num_reads_raw <- colSums(all$data_loaded) #add raw seq depth to map file

# Remove any ASVs with <100 total sequences (likely spurious)
low_abun_ASVs <- rownames(all$data_loaded)[rowSums(all$data_loaded)<100]
length(low_abun_ASVs); length(low_abun_ASVs)/length(all$data_loaded[,1])  #the total number and prop. of these
all_filt <- filter_taxa_from_input(input = all, taxa_IDs_to_remove = low_abun_ASVs)

## taxonomy curation
# Note Euk ASV in 'dysbiotic' male bee from AEC dataset
bee_AEC35 <- filter_data(all_filt, filter_cat = "SampleID", keep_vals = "AEC-35")
return_top_taxa(bee_AEC35,n=10) # ASV_77 classified as Eukaryota
# ASV_77 BLASTs 100% seq id to Aspergillus niger and Mucor sp.
# => replace NA genus with family level classification, and add Fungi as genus classification for ASV_77
levels_to_add <- as.character(unique(all_filt$taxonomy_loaded$taxonomy5[all_filt$taxonomy_loaded$taxonomy6=="NA"]))
levels(all_filt$taxonomy_loaded$taxonomy6) <- c(levels(all_filt$taxonomy_loaded$taxonomy6),levels_to_add, "Fungi")
all_filt$taxonomy_loaded$taxonomy6[all_filt$taxonomy_loaded$taxonomy6=="NA"] <- all_filt$taxonomy_loaded$taxonomy5[all_filt$taxonomy_loaded$taxonomy6=="NA"]
all_filt$taxonomy_loaded$taxonomy6[all_filt$taxonomy_loaded$taxonomy8=="ASV_77"] <- "Fungi"
#for remaining NA at L=6, label 'unidentified'
levels(all_filt$taxonomy_loaded$taxonomy6) <- c(levels(all_filt$taxonomy_loaded$taxonomy6),"unidentified")
all_filt$taxonomy_loaded$taxonomy6[all_filt$taxonomy_loaded$taxonomy6=="NA"] <- "unidentified"

#remove mitochondria and Cplasts
all_filt_noMitoCplast <- filter_taxa_from_input(all_filt, taxa_to_remove = c("Mitochondria","Chloroplast"))
all_filt_noMitoCplast$map_loaded$prop_MitoCplast <- (1-colSums(all_filt_noMitoCplast$data_loaded) / colSums(all_filt$data_loaded))
return_top_taxa(all_filt_noMitoCplast, n=25)
# all typical bee gut bacteria, except Sodalis (x-contamination control) and Corynebacterium and Klebsiella

#examine negative controls:
neg_controls <- all_filt_noMitoCplast$map_loaded$SampleID_dup[all_filt_noMitoCplast$map_loaded$SampleType =="ExtractionBlank"]
neg_controls; length(neg_controls)
#note that NTC_iSeqRun2 (a PCR NTC) has very high seq depth, probably honeybee gut microbiome DNA
neg_control_data <- filter_data(all_filt_noMitoCplast, filter_cat = "SampleType", keep_vals = "ExtractionBlank")
sort(colSums(neg_control_data$data_loaded)) # seq depth
mean(sort(colSums(neg_control_data$data_loaded))) # mean seq depth
return_top_taxa(neg_control_data, n=10) #Burkholderia, Pseudomonas, Nitrospirillum, Sphingomonas, Rhizobiaceae; also some bee taxa (Schmidhempelia)
neg_control_data.Burk <- filter_taxa_from_input(neg_control_data, taxa_to_keep = "Burkholderia")
sum(colSums(neg_control_data.Burk$data_loaded)) / sum(colSums(neg_control_data$data_loaded)) # total prop. of all seqs
neg_control_data.Pseud <- filter_taxa_from_input(neg_control_data, taxa_to_keep = "Pseudomonas")
sum(colSums(neg_control_data.Pseud$data_loaded)) / sum(colSums(neg_control_data$data_loaded)) # total prop. of all seqs

#final sample set
bees <- filter_data(all_filt_noMitoCplast, "SampleType", keep_vals = "Bee")

#check for honey bee contamination w/ Frischella and Bartonella as indicator taxa
HB_check <- filter_taxa_from_input(bees, taxa_to_keep = c("Bartonella","Frischella"))
prop_BartFrisch_in_bumble_samples <- colSums(HB_check$data_loaded) / colSums(bees$data_loaded)
#median and max proportion of honey bee bacteria in bee samples
mean(prop_BartFrisch_in_bumble_samples); median(prop_BartFrisch_in_bumble_samples); max(prop_BartFrisch_in_bumble_samples) 
#very low

#check for cross-contamination using Sodalis + control
sod_check <- filter_taxa_from_input(bees, taxa_to_keep = c("Sodalis"))
prop_Sod_in_bees <- colSums(sod_check$data_loaded) / colSums(bees$data_loaded)
mean(prop_Sod_in_bees); median(prop_Sod_in_bees); max(prop_Sod_in_bees) #median and max proportion of Sodalis in bee samples

# rarefy
rev(sort(colSums(bees$data_loaded))) #distribution of library sizes
rar_depth <- 2200 #3 samples discarded

set.seed(100)
bees_rar <- single_rarefy(bees, depth = rar_depth)
sort(colSums(bees_rar$data_loaded)) #confirm






############## qPCR-corrected taxon abundance estimates
############## hindgut samples only (commercial dataset)
############## 

bees_mod <- filter_data(bees_rar, filter_cat = "GutType", keep_vals = "HG")

#add qPCR data to metadata
qpcr <- read.table("submission_files/qpcr_data.txt", header=T) %>%
  dplyr::select(SampleID,Amount_SYBR_Copies) %>%
  #see qpcr code for more details
  mutate(Amount_SYBR_Copies_gut = Amount_SYBR_Copies * 80) %>%
  group_by(SampleID) %>%
  dplyr::summarize(qPCR_mean_16Scopies_gut = mean(Amount_SYBR_Copies_gut))
bees_mod$map_loaded <- inner_join(bees_mod$map_loaded, qpcr, by="SampleID")

#convert counts to proportions (note: already rarefied to <rar_depth>)
colSums(bees_mod$data_loaded)
temp_df <- bees_mod$data_loaded / rar_depth
colSums(temp_df) #confirm
# now, multiply ASV proportions by total 16S counts to get absolute counts
identical(bees_mod$map_loaded$SampleID_dup,colnames(temp_df)) #confirm that samples in the same order btwn ASV table and map file
# with transposing, multiplies proportions by total 16S counts for each sample
asv_counts <- t(t(temp_df) * bees_mod$map_loaded$qPCR_mean_16Scopies_gut)
bees_mod$data_loaded <- as.data.frame(asv_counts) #replace ASV table
colSums(bees_mod$data_loaded) #confirm

# genus-level plot  **** note relative set to F ****
# code otherwise copied from below taxon summary

ts_all_mod <- summarize_taxonomy(bees_mod, level = 6, report_higher_tax = T, relative=F) #level 6 = genus level
long_ts_all_mod <- ts_all_mod %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_dup") %>%  #this moves the sample ID's from being row names to being their own column
  pivot_longer(!SampleID_dup, names_to="Micro_Genus", values_to="count") %>%  #"pivot" from wide to long format, everything but the Sample ID column
  full_join(bees_mod$map_loaded, by="SampleID_dup") #merge with the metadata

#keep most abundant genera
#(sames as ones defined as dominant, below)
num_taxa <- 8
topX_taxa_mod <- long_ts_all_mod %>%
  group_by(Micro_Genus) %>%
  dplyr::summarize(mean_count=mean(count)) %>%
  arrange(desc(mean_count)) %>%
  head(n=num_taxa)
topX_taxa_mod
long_ts_all_top_mod <- filter(long_ts_all_mod, Micro_Genus %in% topX_taxa_mod$Micro_Genus)

#split the taxonomy string into multiple columns
long_ts_all_top_mod <- long_ts_all_top_mod %>%
  mutate(Micro_Genus = gsub(" ","",Micro_Genus)) %>% #remove space in between taxa names
  separate(col=Micro_Genus,sep=";",remove = T, into=c("Domain", #may need to add "Species" if you get that classification with SILVA
                                                      "Phylum",
                                                      "Class",
                                                      "Order",
                                                      "Family",
                                                      "Micro_Genus")) %>%
  #manual curation of genera names
  mutate(Micro_Genus = gsub("Bifidobacteriaceae","Unclassified Bifidobacteriaceae",Micro_Genus),
         Micro_Genus = gsub("CandidatusSchmidhempelia","Schmidhempelia",Micro_Genus))

#order genera by their overall abundance (mean across all samples)
top_genera_ordered_mod <- long_ts_all_top_mod %>%
  group_by(Micro_Genus) %>%
  dplyr::summarize(mean_count=mean(count)) %>%
  arrange(desc(mean_count))

#reorder genera for plot, from highest->lowest abundance
long_ts_all_top_mod$Micro_Genus <- factor(long_ts_all_top_mod$Micro_Genus,top_genera_ordered_mod$Micro_Genus)

#calculate log10, reset (0) values
long_ts_all_top_mod$log10_count <- log10(long_ts_all_top_mod$count)
long_ts_all_top_mod$log10_count[long_ts_all_top_mod$log10_count=="-Inf"] <- 0

#filter, to only show data after the colonization phase (day 0 -> ~4)
long_ts_all_top_mod_after_col <- filter(long_ts_all_top_mod, Age_days > 4)

#plot

dotplot.by.age.abs.abun <- ggplot(long_ts_all_top_mod_after_col, aes(x=Age_days, y=log10_count)) +
  geom_point(size=3, alpha=0.4, color="blue") +
  theme_bw() +
  scale_x_continuous(breaks=seq(from=5, to=max(long_ts_all_top_mod_after_col$Age_days), by=10)) +
  xlab("Age (days)") + ylab ("log10(# 16S rRNA gene copies)") +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14)) + 
  facet_wrap(~ Micro_Genus, nrow=2, ncol=4, scales = "free") +
  theme(strip.text = element_text(face = "italic", size=9))
dotplot.by.age.abs.abun

# test whether Bombiscardovia indeed changes in est. absolute abundance w/ age 
bombis_abs <- filter(long_ts_all_top_mod_after_col, Micro_Genus=="Bombiscardovia" )
m_b_abs <- lme(log10_count ~ Age_days, random=~1|Colony,data=bombis_abs) 
summary(m_b_abs)


############## Taxon summary
############## and stacked bar plots
##############

#genus-level summary for stacked bar chart
ts_all_rar <- summarize_taxonomy(bees_rar, level = 6, report_higher_tax = T) #level 6 = genus level
long_ts_all <- ts_all_rar %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_dup") %>%  #this moves the sample ID's from being row names to being their own column
  pivot_longer(!SampleID_dup, names_to="Micro_Genus", values_to="prop") %>%  #"pivot" from wide to long format, everything but the Sample ID column
  full_join(bees_rar$map_loaded, by="SampleID_dup") #merge with the metadata

#keep only the dominant genera,
#defined as â‰¥1% [rounded] mean proportion across samples
#note this is ALL bee gut samples (HG and MG, TJH and AEC samples)
topX_taxa <- long_ts_all %>%
  group_by(Micro_Genus) %>%
  dplyr::summarize(mean_prop=mean(prop)) %>%
  arrange(desc(mean_prop)) %>%
  filter(mean_prop>0.009)
topX_taxa 

long_ts_all_top <- filter(long_ts_all, Micro_Genus %in% topX_taxa$Micro_Genus)

#split the taxonomy string into multiple columns
long_ts_all_top <- long_ts_all_top %>%
  mutate(Micro_Genus = gsub(" ","",Micro_Genus)) %>% #remove space in between taxa names
  separate(col=Micro_Genus,sep=";",remove = T, into=c("Domain",
                                                      "Phylum",
                                                      "Class",
                                                      "Order",
                                                      "Family",
                                                      "Micro_Genus")) %>%
  #manual curation of genus names
  mutate(Micro_Genus = gsub("-Caballeronia-Paraburkholderia","",Micro_Genus)) %>%
  mutate(Micro_Genus = gsub("CandidatusSchmidhempelia","Candidatus Schmidhempelia",Micro_Genus)) %>%
  mutate(Micro_Genus = gsub("Bifidobacteriaceae","Unclassified Bifidobacteriaceae",Micro_Genus))

#order genera by their overall abundance (mean across all samples)
#note this is for ALL bee gut samples (HG and MG, TJH and AEC samples)
#hence, ranks will vary in subgroups of samples
top_genera_ordered <- long_ts_all_top %>%
  group_by(Micro_Genus) %>%
  dplyr::summarize(mean_prop=mean(prop)) %>%
  arrange(desc(mean_prop))

#reorder genera for plot, from highest->lowest overall abundance
long_ts_all_top$Micro_Genus <- factor(long_ts_all_top$Micro_Genus,top_genera_ordered$Micro_Genus)



### stacked bar plot

#TJH sample set

long_ts_all_top$Colony <- factor(long_ts_all_top$Colony, c("B","Y","W"))
long_ts_all_top$Age_days_discrete <- as.factor(long_ts_all_top$Age_days)

stackbar_dat_TJH_HG <- long_ts_all_top %>%  filter(Project=="TJH_BAM", GutType=="HG")
stackbar_dat_TJH_MG <- long_ts_all_top %>%  filter(Project=="TJH_BAM", GutType=="MG")

# slight modification to AEC dataset top taxa, to add Fungi and Klebsiella abundant in one bee
# also to add Burkholderia which is only abundant in the day-zero bees
long_ts_all_top_AEC <- filter(long_ts_all, Micro_Genus %in% c("Eukaryota; NA; NA; NA; NA; Fungi", #ASV_77 string
                                                              "Bacteria; Proteobacteria; Gammaproteobacteria; Burkholderiales; Burkholderiaceae; Burkholderia-Caballeronia-Paraburkholderia",
                                                              "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Klebsiella",
                                                              topX_taxa$Micro_Genus))
long_ts_all_top_AEC <- long_ts_all_top_AEC %>%
  mutate(Micro_Genus = gsub(" ","",Micro_Genus)) %>% #remove space in between taxa names
  separate(col=Micro_Genus,sep=";",remove = T, into=c("Domain", #may need to add "Species" if you get that classification with SILVA
                                                      "Phylum",
                                                      "Class",
                                                      "Order",
                                                      "Family",
                                                      "Micro_Genus")) %>%
  mutate(Micro_Genus = gsub("-Caballeronia-Paraburkholderia","",Micro_Genus)) %>%
  mutate(Micro_Genus = gsub("CandidatusSchmidhempelia","Candidatus Schmidhempelia",Micro_Genus)) %>%
  mutate(Micro_Genus = gsub("Bifidobacteriaceae","Unclassified Bifidobacteriaceae",Micro_Genus))
top_genera_ordered_AEC <- long_ts_all_top_AEC %>%
  group_by(Micro_Genus) %>%
  dplyr::summarize(mean_prop=mean(prop)) %>%
  arrange(desc(mean_prop))
long_ts_all_top_AEC$Micro_Genus <- factor(long_ts_all_top_AEC$Micro_Genus,top_genera_ordered_AEC$Micro_Genus)
long_ts_all_top_AEC$Colony <- factor(long_ts_all_top_AEC$Colony, c("B","Y","W"))
long_ts_all_top_AEC$Age_days_discrete <- as.factor(long_ts_all_top_AEC$Age_days)



#custom color scheme for top taxa
genera <- c(levels(long_ts_all_top$Micro_Genus),"Burkholderia", "Klebsiella","Fungi")
genera
colors <- c("#d7191c","#ffffbf","#fdae61", #schmid, snod, gill
             "#3182bd","#deebf7", #lacto, apilacto
            "#bcbddc","#9ecae1","#756bb1",  #Bifidobac, Bombilacto, Bombisc
            "black", # Burk, only abundant in day zero bees
            "green", "darkgrey") ## klebsiella and Fungi; only abundant in AEC data
colors_df <- data.frame(genera,colors,stringsAsFactors = F)


## TJH samples, hindgut, day 1+
stackbar_dat_TJH_HG_day1plus <- filter(stackbar_dat_TJH_HG, Age_days>0)
colors_df_NoCand <- colors_df %>%
  mutate(genera = gsub("Candidatus ","",genera))

bars.plot.TJH_HG <- ggplot(stackbar_dat_TJH_HG_day1plus, aes(Age_days_discrete, prop, fill=Micro_Genus)) +
  geom_bar(stat="identity") + ylim(0,1) + theme_bw() +
  theme(axis.text = element_text(size=10), axis.title.y = element_text(size=12)) +
  xlab("Age (days)") + ylab ("Proportion of sequences") + 
  scale_fill_manual(values=colors_df_NoCand$colors,
                    labels=colors_df_NoCand$genera,
                    name="Genus") +
  guides(fill=guide_legend(label.theme = element_text(face = "italic",size=10),
                           title.theme = element_text(size=12))) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(size=8, angle=45,hjust=1),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=12),
        strip.text = element_text(size=12)) + 
  facet_wrap(~Colony,nrow=3,ncol=1, labeller = label_both)
bars.plot.TJH_HG


# summary version of adult community composition, to include with larval data

means.stackbar_dat_TJH_HG <- stackbar_dat_TJH_HG %>%
  group_by(Micro_Genus) %>%
  dplyr::summarise(mean_prop = mean(prop)) %>%
  mutate(sample=as.factor("1"))

means.bars.plot.TJH_HG <- ggplot(means.stackbar_dat_TJH_HG, aes(sample, mean_prop, fill=Micro_Genus)) +
  geom_bar(stat="identity") + ylim(0,1) + theme_bw() +
  theme(axis.text = element_text(size=10), axis.title.y = element_text(size=12)) +
  xlab("Adult hindguts") + ylab ("Mean proportion of sequences") + 
  scale_fill_manual(values=colors_df_NoCand$colors,
                    labels=colors_df_NoCand$genera,
                    name="Genus") +
  guides(fill=guide_legend(label.theme = element_text(face = "italic",size=10),
                           title.theme = element_text(size=12))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=12),
        strip.text = element_text(size=12))
means.bars.plot.TJH_HG



## TJH samples, midgut, day 1+
stackbar_dat_TJH_MG_day0 <- filter(stackbar_dat_TJH_MG,
                                        Age_days_discrete != "0") #remove newly emerged bees for this plot; see below
bars.plot.TJH_MG <- ggplot(stackbar_dat_TJH_MG_day0, aes(Age_days_discrete, prop, fill=Micro_Genus)) +
  geom_bar(stat="identity") + ylim(0,1) + theme_bw() +
  theme(axis.text = element_text(size=10), axis.title.y = element_text(size=12)) +
  xlab("Age (days)") + ylab ("Proportion of sequences") + 
  scale_fill_manual(values=colors_df_NoCand$colors,
                    labels=colors_df_NoCand$genera,
                    name="Genus") +
  guides(fill=guide_legend(label.theme = element_text(face = "italic",size=10),
                           title.theme = element_text(size=12))) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(size=8, angle=45,hjust=1),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=12),
        strip.text = element_text(size=12)) + 
  facet_wrap(~Colony,nrow=3,ncol=1, labeller = label_both)
bars.plot.TJH_MG



## TJH samples, midgut + hindgut, day 0 only
stackbar_dat_TJH_day0 <- long_ts_all_top_AEC %>%
  filter(Project=="TJH_BAM",Age_days_discrete == "0") %>%
  filter(Micro_Genus != "Fungi") #remove fungi bc not abundant in these samples
stackbar_dat_TJH_day0$GutType <- factor(stackbar_dat_TJH_day0$GutType, c("MG","HG")) #reorder for plot

bars.plot.TJH_day0 <- ggplot(stackbar_dat_TJH_day0, aes(BeeID, prop, fill=Micro_Genus)) +
  geom_bar(stat="identity") + ylim(0,1) + theme_bw() +
  theme(axis.text = element_text(size=10), axis.title.y = element_text(size=12)) +
  xlab("Bee ID") + ylab ("Proportion of sequences") + 
  scale_fill_manual(values=colors_df_NoCand$colors,
                    labels=colors_df_NoCand$genera,
                    name="Genus") +
  guides(fill=guide_legend(label.theme = element_text(face = "italic",size=10),
                           title.theme = element_text(size=12))) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(size=8, angle=45,hjust=1),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=12),
        strip.text = element_text(size=12)) + 
  facet_grid(GutType~Colony, scales="free_x",labeller = label_both)
bars.plot.TJH_day0


## AEC samples (no queens)
stackbar_dat_AEC  <- long_ts_all_top_AEC %>%  filter(Project=="AEC")
stackbar_dat_AEC_noQ <- filter(stackbar_dat_AEC, Caste != "queen", #no age data
                               Micro_Genus !="Burkholderia")  #almost nonexistent in this dataset
                                         
stackbar_dat_AEC_noQ <- stackbar_dat_AEC_noQ %>%
  mutate(Age_SampleID = paste(Age_days_discrete, SampleID_dup, sep="_")) %>%
  arrange(Age_days)
sampleorder <- unique(stackbar_dat_AEC_noQ$Age_SampleID)

stackbar_dat_AEC_noQ$Age_SampleID <- factor(stackbar_dat_AEC_noQ$Age_SampleID, sampleorder)

#edit colors to exclude Burk.
colors_df_AEC <- filter(colors_df_NoCand, genera != "Burkholderia")

bars.plot.AEC <- ggplot(stackbar_dat_AEC_noQ, aes(Age_SampleID, prop, fill=Micro_Genus)) +
  geom_bar(stat="identity") + ylim(0,1) + theme_bw() +
  theme(axis.text = element_text(size=10), axis.title.y = element_text(size=12)) +
  xlab("Age_BeeID") + ylab ("Proportion of sequences") + 
  scale_fill_manual(values=colors_df_AEC$colors,
                    labels=colors_df_AEC$genera,
                    name="Genus") +
  guides(fill=guide_legend(label.theme = element_text(face = "italic",size=10),
                           title.theme = element_text(size=12))) +
  theme(panel.grid.major.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(axis.text.x = element_text(size=8, angle=45,hjust=1),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=12),
        strip.text = element_text(size=12)) + 
  facet_grid(Species~Caste, scales="free_x",labeller = label_both)
bars.plot.AEC






############## Taxa associated with age POST-DAY 4
############## and dot plots
##############

# note, this is using only the top taxa (genera) identified above, not low-abundance taxa

dotplot_df_TJH <- stackbar_dat_TJH_HG %>%
  filter(Age_days >4) %>%   # filter out the youngest bees
  dplyr::select(BeeID,Age_days,Colony,Micro_Genus,prop) %>%
  mutate(Micro_Genus = gsub(" ","",Micro_Genus)) 

dotplot_df_AEC <- stackbar_dat_AEC_noQ %>%
  dplyr::select(BeeID,Age_days,Colony,Micro_Genus,prop) %>%
  mutate(Micro_Genus = gsub(" ","",Micro_Genus))


### stats for TJH HG samples:
### age X prop. nonparametric correlations (spearman) for all dominant taxa,
## and correct for mult comp
# code adapted from https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

#Define functions to pass to map
cor_model <- function(df){
  dfdf=data.frame(df)
  cor.test(x=dfdf$prop, y=dfdf$Age_days, method = "spearman")
}
cor_pval <- function(df){
  dfdf=data.frame(df)
  cor.test(x=dfdf$prop, y=dfdf$Age_days, method = "spearman")$p.value
}
cor_rho <- function(df){
  dfdf=data.frame(df)
  cor.test(x=dfdf$prop, y=dfdf$Age_days, method = "spearman")$estimate
}
#Create nested data frames by genus and loop over each using map 
cor_results <- dotplot_df_TJH %>%
  dplyr::select(-BeeID,-Colony) %>%
  group_by(Micro_Genus) %>%
  nest() %>%
  mutate(cor_test = map(data, cor_model),
         p_value = map(data, cor_pval),
         rho = map(data, cor_rho))                       

#unnest results
cor_results <- cor_results %>%
  dplyr::select(Micro_Genus, p_value, rho) %>%
  unnest(cols=c(p_value,rho))
cor_results$p_fdr <- p.adjust(cor_results$p_value, method="fdr") #wasnt able to get working with mutate in-place
cor_results <- filter(cor_results, p_fdr < 0.05) %>% arrange(p_fdr)
cor_results


# additional stats for taxa IDd as sig. assoc. with age:
#linear mixed-effects model
#https://ourcodingclub.github.io/tutorials/mixed-models/#FERE
#colony as random effect (random intercept model)

dotplot_df_TJH_Schmid <- filter(dotplot_df_TJH, Micro_Genus=="CandidatusSchmidhempelia")
dotplot_df_TJH_Gill <- filter(dotplot_df_TJH, Micro_Genus=="Gilliamella")
dotplot_df_TJH_Bombis <- filter(dotplot_df_TJH, Micro_Genus=="Bombiscardovia")

m_s <- lme(prop ~ Age_days, random=~1|Colony,data=dotplot_df_TJH_Schmid) 
summary(m_s)
plot(m_s)
qqnorm(resid(m_s)); qqline(resid(m_s))

m_g <- lme(prop ~ Age_days, random=~1|Colony,data=dotplot_df_TJH_Gill) 
summary(m_g)
plot(m_g)
qqnorm(resid(m_g)); qqline(resid(m_g))

m_b <- lme(prop ~ Age_days, random=~1|Colony,data=dotplot_df_TJH_Bombis) 
summary(m_b)
plot(m_b)
qqnorm(resid(m_b)); qqline(resid(m_b))



## now dot plot, for only Schmid and Gill
# note, now including ALL AGES EXCEPT NEWs (day 1+)

dotplot_df_TJH_Schmid_Gill <- stackbar_dat_TJH_HG %>%
  filter(Micro_Genus %in% c("Candidatus Schmidhempelia","Gilliamella","Bombiscardovia"),
         Age_days >0) %>%
         mutate(Micro_Genus = gsub("Candidatus ","",Micro_Genus))

dotplot_df_TJH_Schmid_Gill$Micro_Genus <- factor(dotplot_df_TJH_Schmid_Gill$Micro_Genus,
                                                        c("Schmidhempelia","Gilliamella"))

dotplot_Schmid_Gill <- ggplot(dotplot_df_TJH_Schmid_Gill, aes(x=Age_days, y=prop, fill=Colony)) +
  geom_point(size=3, alpha=0.7, shape=21) +
  theme_bw() +
  xlab("Age (days)") + ylab ("Proportion of sequences") + 
  scale_fill_manual(values=c("Y"="#ffeda0",
                              "B"="#9ecae1",
                              "W"="#f0f0f0")) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14),
        strip.text = element_text(face = "italic", size=11)) + 
  facet_wrap(~ Micro_Genus, scales = "free", ncol=2, nrow=2) +
  geom_smooth(method="lm", fill="darkgrey",color="black")
dotplot_Schmid_Gill

## now dot plot, for only schmid and gill 
# compare w/ AEC samples, impatiens only, workers only

dotplot_df_AEC_Schmid_Gill <- stackbar_dat_AEC_noQ %>%
  filter(Micro_Genus %in% c("Candidatus Schmidhempelia","Gilliamella"),
         Species =="impatiens",
         Caste =="worker",
         Age_days >0) %>% #filter out NEWs
        mutate(Micro_Genus = gsub("Candidatus ","",Micro_Genus))

dotplot_df_AEC_Schmid_Gill$Micro_Genus <- factor(dotplot_df_AEC_Schmid_Gill$Micro_Genus,
                                                        c("Schmidhempelia","Gilliamella"))

dotplot_Schmid_Gill_AEC <- ggplot(dotplot_df_AEC_Schmid_Gill, aes(x=Age_days, y=prop)) +
  geom_point(size=3, alpha=0.7, shape=21, fill="grey") +
  theme_bw() +
  xlab("Age (days)") + ylab ("Proportion of sequences") + 
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title = element_text(size=10),
        strip.text = element_text(face = "italic", size=8)) + 
  facet_wrap(~ Micro_Genus, scales = "free", ncol=1, nrow=2)
dotplot_Schmid_Gill_AEC






############## ASV-level plots
############## (for most abundant bacterial genera)
##############

## TJH samples, hindgut, day 1+
sampleIDs_bees_rar_TJH_HG_day1 <- stackbar_dat_TJH_HG %>%
  filter(Age_days_discrete != "0") %>%
  dplyr::select(SampleID_dup) %>%
  distinct()
bees_rar_TJH_HG_day1 <- filter_data(bees_rar, filter_cat = "SampleID_dup", keep_vals = sampleIDs_bees_rar_TJH_HG_day1$SampleID_dup)

snod_data <- filter_taxa_from_input(bees_rar_TJH_HG_day1, taxa_to_keep = "Snodgrassella", at_spec_level = 6)
gill_data <- filter_taxa_from_input(bees_rar_TJH_HG_day1, taxa_to_keep = "Gilliamella", at_spec_level = 6)
schmid_data <- filter_taxa_from_input(bees_rar_TJH_HG_day1, taxa_to_keep = "Candidatus Schmidhempelia", at_spec_level = 6)

#for ASV level plot, filter the taxon summary (with ALL ASVs) at level=8
ts_all_rar_L8 <- summarize_taxonomy(bees_rar_TJH_HG_day1, level = 8, report_higher_tax = T) #level 8 = ASV
long_ts_all_L8 <- ts_all_rar_L8 %>%
  t() %>%  #transpose
  as.data.frame() %>%
  rownames_to_column(var = "SampleID_dup") %>%  #this moves the sample ID's from being row names to being their own column
  pivot_longer(!SampleID_dup, names_to="Micro_Genus", values_to="prop") %>%  #"pivot" from wide to long format, everything but the Sample ID column
  full_join(bees_rar$map_loaded, by="SampleID_dup") %>% #merge with the metadata
  #split the taxonomy string into multiple columns
  mutate(Micro_Genus = gsub(" ","",Micro_Genus)) %>% #remove space in between taxa names
  separate(col=Micro_Genus,sep=";",remove = T, into=c("Domain", #may need to add "Species" if you get that classification with SILVA
                                                      "Phylum",
                                                      "Class",
                                                      "Order",
                                                      "Family",
                                                      "Micro_Genus",
                                                      "Species",
                                                      "ASV"))
long_ts_all_L8_snod <- long_ts_all_L8 %>% filter(ASV %in% unique(snod_data$taxonomy_loaded$taxonomy8))
long_ts_all_L8_gill <- long_ts_all_L8 %>% filter(ASV %in% unique(gill_data$taxonomy_loaded$taxonomy8))
long_ts_all_L8_schmid <- long_ts_all_L8 %>% filter(ASV %in% unique(schmid_data$taxonomy_loaded$taxonomy8))


#sort ASVs by mean proportion, keep only top 10 for visualization
asvs_sorted_top_snod <- long_ts_all_L8_snod %>% group_by(ASV) %>% dplyr::summarize(mean_prop=mean(prop)) %>% arrange(desc(mean_prop)) %>% head(n=5)
long_ts_all_snod_top <- filter(long_ts_all_L8_snod, ASV %in% asvs_sorted_top_snod$ASV)
long_ts_all_snod_top$ASV <- factor(long_ts_all_snod_top$ASV, rev(asvs_sorted_top_snod$ASV))

asvs_sorted_top_gill <- long_ts_all_L8_gill %>% group_by(ASV) %>% dplyr::summarize(mean_prop=mean(prop)) %>% arrange(desc(mean_prop)) %>% head(n=5)
long_ts_all_gill_top <- filter(long_ts_all_L8_gill, ASV %in% asvs_sorted_top_gill$ASV)
long_ts_all_gill_top$ASV <- factor(long_ts_all_gill_top$ASV, rev(asvs_sorted_top_gill$ASV))

asvs_sorted_top_schmid <- long_ts_all_L8_schmid %>% group_by(ASV) %>% dplyr::summarize(mean_prop=mean(prop)) %>% arrange(desc(mean_prop)) %>% head(n=5)
long_ts_all_schmid_top <- filter(long_ts_all_L8_schmid, ASV %in% asvs_sorted_top_schmid$ASV)
long_ts_all_schmid_top$ASV <- factor(long_ts_all_schmid_top$ASV, rev(asvs_sorted_top_schmid$ASV))

#create x axis color mapping
axis_colors <- long_ts_all_snod_top %>% dplyr::select(SampleID, Colony) %>% distinct() %>% group_by(Colony) %>% dplyr::count() %>%
  mutate(color = case_when(Colony == "B" ~ "#9ecae1",
                           Colony == "W" ~ "#f0f0f0",
                           Colony == "Y" ~ "#ffeda0"))
axis_colors_vector <- c(rep(axis_colors$color[1],axis_colors$n[1]),
                        rep(axis_colors$color[2],axis_colors$n[2]),
                        rep(axis_colors$color[3],axis_colors$n[3]))

heatmap_plot_snod <- ggplot(long_ts_all_snod_top, aes(SampleID, ASV)) +
  geom_tile(aes(fill=prop), color="white") +
  scale_fill_gradient(name="proportion",
                      low = "white",high = "black") +
  scale_x_discrete(labels=rep("-",94)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, angle=0,
                                   color=axis_colors_vector),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        legend.position = "none") +
  xlab("") +
  ylab("Snodgrassella") +
  ggtitle("Top 5 ASVs in each genus")

heatmap_plot_gill<- ggplot(long_ts_all_gill_top, aes(SampleID, ASV)) +
  geom_tile(aes(fill=prop), color="white") +
  scale_fill_gradient(name="proportion",
                      low = "white",high = "black") +
  scale_x_discrete(labels=rep("-",94)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, angle=0,
                                   color=axis_colors_vector),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12)) +
  xlab("") +
  ylab("Gilliamella")

heatmap_plot_schmid <- ggplot(long_ts_all_schmid_top, aes(SampleID, ASV)) +
  geom_tile(aes(fill=prop), color="white") +
  scale_fill_gradient(name="proportion",
                      low = "white",high = "black") +
  scale_x_discrete(labels=rep("-",94)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, angle=0,
                                   color=axis_colors_vector),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size=12),
        legend.position = "none") +
  xlab("Bees sorted by colony") +
  ylab("Schmidhempelia")


heatmap_plot_snod /
heatmap_plot_gill /
heatmap_plot_schmid







############## alpha diversity
##############
############## 

## as above, use  TJH samples, hindgut, day 1+
sampleIDs_bees_rar_TJH_HG_day1 <- stackbar_dat_TJH_HG %>%
  filter(Age_days_discrete != "0") %>%
  dplyr::select(SampleID_dup) %>%
  distinct()
bees_rar_TJH_HG_day1 <- filter_data(bees_rar, filter_cat = "SampleID_dup", keep_vals = sampleIDs_bees_rar_TJH_HG_day1$SampleID_dup)

bees_rar_TJH_HG_day1$map_loaded$shannon <- diversity(t(bees_rar_TJH_HG_day1$data_loaded), index = "shannon") #Shannon per individual

plot_shan.by.age <- ggplot(bees_rar_TJH_HG_day1$map_loaded, aes(x=Age_days, y=shannon, fill=Colony)) +
  geom_point(size=3, alpha=0.7, shape=21) +
  theme_bw() +
  xlab("Age (days)") + ylab ("Shannon diversity") + 
  scale_fill_manual(values=c("Y"="#ffeda0",
                             "B"="#9ecae1",
                             "W"="#f0f0f0")) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14),
        legend.position = "none") + 
  facet_wrap(~ Colony, scales = "fixed", labeller = label_both, nrow=3, ncol=1)
plot_shan.by.age

#stats for shannon changes AFTER COLONIZATION (not included in paper)
# see above dotplot stats for same age division
shan <- bees_rar_TJH_HG_day1$map_loaded %>%
  filter(Age_days >4)   # filter out the youngest bees

#LME model, colony as random effect
m_shan <- lme(shannon ~ Age_days, random=~1|Colony,data=shan) 
summary(m_shan) 
hist(resid(m_shan))
qqnorm(resid(m_shan)); qqline(resid(m_shan))






############## beta diversity
############## (also, TJH HG samples, day 1+)
############## 

#ordination
bc.dm <- calc_dm(bees_rar_TJH_HG_day1$data_loaded, method = "bray_sq_trans")
ord <- calc_ordination(bc.dm, 'nmds')
nmds_df <- cbind.data.frame(ord, bees_rar_TJH_HG_day1$map_loaded)

new_to_old_colorscale <- c("#762A83","#E7D4E8","#D9F0D3","#1B7837") #matches colors from strain plot, rnaseq plot

ord.age.plot <- ggplot(nmds_df, aes(x=MDS1, y=MDS2,fill=Age_days)) +
  theme_bw() +
  geom_point(size=3, alpha=0.85, pch=21) +
  scale_fill_gradientn(name="Age (days)",
                       colors = new_to_old_colorscale) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size=12),
        legend.text = element_text(size=11),
        legend.title = element_text(size=12)) +
  xlab("MDS Axis 1") + ylab("MDS Axis 2") 
ord.age.plot


### db-RDA

bee.dat <- bees_rar_TJH_HG_day1$data_loaded
bee.dat.t <- t(bee.dat)
bee.meta <- bees_rar_TJH_HG_day1$map_loaded

#here include both colony and age as predictors in the model
dbRDA <- capscale(bee.dat.t ~ Age_days+Colony, bee.meta, dist="bray")

plot(dbRDA) #base R biplot; red + are ASVs, black O are bees, X are colony centroids
anova(dbRDA) #of model
anova(dbRDA, by="terms", permu=999)
