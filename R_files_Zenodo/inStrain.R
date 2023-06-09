library(tidyverse)
library(patchwork)
library(vroom)

################ inStrain compare
################ pairwise popANI-based analysis
################ Fig. S9

inStrain_popANI <- read.table("resubmission_files/inStrain_popANI.txt", sep='\\t',header=T)

#merge with bin IDs from GTDB-tk
bin.tax <- read.table("resubmission_files/gtdbtk.bac120.summary.tsv", sep='\\t', header=T) %>%
  dplyr::rename(genome = user_genome) %>%
  dplyr::select(genome,classification) %>%
  separate(col=classification,sep=";",remove = F, into=c("Domain",
                                                         "Phylum",
                                                         "Class",
                                                         "Order",
                                                         "Family",
                                                         "Genus",
                                                         "Species"))

inStrain_popANI_wTax <- bin.tax %>%
  inner_join(inStrain_popANI,by="genome") %>%
  dplyr::select(genome,name1,name2,popANI,Genus)

#merge with sample metadata (for each sample in the comparison pair)
metadata <- read.table("resubmission_files/sample_metadata_v2.txt",header=T,sep='\\t') %>%
  dplyr::select(novogene_id_shotgun,Age_days,Colony)

inStrain_popANI_wTax %<>% dplyr::rename(novogene_id_shotgun=name1) #rename to join for 1st sample

inStrain_popANI_wTax_merge1 <- inStrain_popANI_wTax %>%
  inner_join(metadata, by="novogene_id_shotgun") %>%
  dplyr::rename(name.x = novogene_id_shotgun,
                novogene_id_shotgun = name2)
inStrain_popANI_wTax_merge2 <- inStrain_popANI_wTax_merge1 %>%
  inner_join(metadata, by="novogene_id_shotgun") %>%
  dplyr::rename(name.y = novogene_id_shotgun)

head(inStrain_popANI_wTax_merge2) # x corresponds to "name1" (the 1st sample), y corresponds to "name2"

#if/else same vs different colony comparison
df <- inStrain_popANI_wTax_merge2 %>%
  #use if else to classify comparisons by whether the pair is from the same colony
  mutate(colony_comparison = ifelse(Colony.x == Colony.y, "same","different"))

###example plot for %popANI ~ same-vs-different colony
filter(df, Genus=="g__Schmidhempelia") %>% ggplot(aes(x=colony_comparison, y=popANI)) + geom_boxplot()

#convert to same/different strain based on different popANI thresholds
df_strains_shared <- df %>%
  mutate(strains_99.999 = ifelse(popANI > 0.99999, "same_strain","different_strain")) %>%
  mutate(strains_99.99 = ifelse(popANI > 0.9999, "same_strain","different_strain")) %>%
  mutate(strains_99.9 = ifelse(popANI > 0.999, "same_strain","different_strain"))

#plot strain sharing X colony-sharing for different threshold levels
plot_99.999 <- ggplot(df_strains_shared, aes(x=colony_comparison, fill=strains_99.999)) +
  geom_bar(position="fill") +
  facet_wrap(~ genome, scales="free_y") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("Prop. comparisons") + ggtitle("99.999% popANI")
plot_99.99 <- ggplot(df_strains_shared, aes(x=colony_comparison, fill=strains_99.99)) +
  geom_bar(position="fill") +
  facet_wrap(~ genome, scales="free_y") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("Prop. comparisons") + ggtitle("99.99% popANI")
plot_99.9 <- ggplot(df_strains_shared, aes(x=colony_comparison, fill=strains_99.9)) +
  geom_bar(position="fill") +
  facet_wrap(~ genome, scales="free_y") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("Prop. comparisons") + ggtitle("99.9% popANI")

plot_99.999 / plot_99.99 / plot_99.9
# 99.99 seems to have the best resolution



### chi sq test for strain-sharing (99.99% only) ~ colony vs age class, by genome

#add discrete-age-class-comparison to df_strains_shared
#use classes defined in the RNAseq analysis
df_strains_shared %<>%
  mutate(
  discrete_ages.x = cut(Age_days.x, breaks=c(-Inf,1,19,43, Inf),
                           labels=c("new","young","middle","old")),
  discrete_ages.y = cut(Age_days.y, breaks=c(-Inf,1,19,43, Inf),
                        labels=c("new","young","middle","old")),
  age_comparison = ifelse(discrete_ages.x == discrete_ages.y, "same_age_group","different_age_group"))
#confirm
df_strains_shared %>% arrange(desc(Age_days.x)) %>% dplyr::select(age_comparison,Age_days.x,Age_days.y,discrete_ages.x,discrete_ages.y) %>% 
  head(n=10)

#all genomes ~ colony and ~ age class
genomes <- unique(df_strains_shared$genome)
num_genomes <- length(genomes)
chisqs.df <- data.frame(character(num_genomes),numeric(num_genomes),numeric(num_genomes))
names(chisqs.df) <- c("genome","col_chisq.pval","age_chisq.pval")
for (Q in 1:num_genomes) {
  genome.Q <- filter(df_strains_shared, genome==genomes[Q])
  chisqs.df$genome[Q] <- genomes[Q]
  #chisq test for colony membership
  col.test.Q <- chisq.test(table(genome.Q$colony_comparison, genome.Q$strains_99.99))
  print(col.test.Q$expected)
  chisqs.df$col_chisq.pval[Q] <- col.test.Q$p.value
  #chisq test for age group membership
  age.test.Q <- chisq.test(table(genome.Q$age_comparison, genome.Q$strains_99.99))
  print(age.test.Q$expected)
  chisqs.df$age_chisq.pval[Q] <- age.test.Q$p.value
}
# note that D_601_12_bin.6 has expected values < 5, hence flag
# p vals are >>> 0.05; visualizing, proportions are similar

chisqs.df <- chisqs.df %>%
  mutate(col_chisq.pval.FDR.adj = p.adjust(col_chisq.pval, method = "fdr"),
         col.sig = ifelse(col_chisq.pval.FDR.adj < 0.05, "sig","ns"),
         age_chisq.pval.FDR.adj = p.adjust(age_chisq.pval, method = "fdr"),
         age.sig = ifelse(age_chisq.pval.FDR.adj < 0.05, "sig","ns")) %>%
  arrange(col_chisq.pval.FDR.adj)
head(chisqs.df)



#merge with MAG metadata for plotting
# (colony test only)

chisqs.df.w.tax <- inner_join(chisqs.df, bin.tax, by="genome") %>%
  dplyr::select(Genus,col.sig)

#sort genomes by colony test p val
chisqs.df.for_plot <- chisqs.df %>%
  inner_join(bin.tax, by="genome") %>%
  mutate(Genus = gsub("g__","",Genus),
         Genus = substr(Genus,start=1,stop=6),
         Genus_Genome = paste(Genus,genome,sep="-")) 

df_strains_shared.for_plot <- df_strains_shared %>%
  mutate(Genus = gsub("g__","",Genus),
         Genus = substr(Genus,start=1,stop=6),
         Genus_Genome = paste(Genus,genome,sep="-"))

df_strains_shared.for_plot <- left_join(df_strains_shared.for_plot,chisqs.df.for_plot,by="Genus_Genome")
#sorted by p val
df_strains_shared.for_plot$Genus_Genome <- factor(df_strains_shared.for_plot$Genus_Genome, chisqs.df.for_plot$Genus_Genome)
#sort by sig then nonsig
df_strains_shared.for_plot$col.sig <- factor(df_strains_shared.for_plot$col.sig, c("sig","ns"))

strains_by_colony_plot <- ggplot(df_strains_shared.for_plot,
                                 aes(x=colony_comparison,
                                     fill=strains_99.99,
                                     alpha=col.sig)) +
  geom_bar(position="fill") +
  facet_wrap(~ Genus_Genome) +
  scale_fill_manual(name="Strain",
                    labels=c("different","same"),
                    values=c("#d8b365","#5ab4ac")) +
  scale_alpha_manual(name="Ï‡2 test",
                     labels=c("FDR-adj. p < 0.05",
                              "FDR-adj. p > 0.05"),
                     values=c(1.0,0.3)) +
  theme_bw() +
  theme(strip.text = element_text(size=8),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Proportion of comparisons") + xlab("Colony")
strains_by_colony_plot
