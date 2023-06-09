library(vroom)
library(tidyverse)
library(nlme)

#load scaffold-to-bin (STB) file
stb <- read.table("resubmission_files/derep_98_genomes.stb", sep='\\t', header=F)
colnames(stb) <- c("scaffold","genome")
stb$genome <- gsub(".fa","",stb$genome)

#load  metadata
metadata <- read.table("resubmission_files/sample_metadata_v2.txt",header=T,sep='\\t') %>%
  dplyr::select(novogene_id_shotgun,Age_days,Colony)

#load bin metadata: bin IDs from GTDB-tk
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

## analyze RPKM from pileup output

all_mag_mapping <- read.table("resubmission_files/MAG_mapping_data.txt",header=T,sep='\\t') %>%
  #merge with STB file
  inner_join(stb,by="scaffold") %>%
  mutate(novogene_id_shotgun = as.factor(novogene_id_shotgun))

RPM_samples <- all_mag_mapping %>%
  #calculate total mapped reads-per-million for each sample
  group_by(novogene_id_shotgun) %>%
  dplyr::summarise(num_reads_mapped_per_M = (sum(Reads)/1000000))

# calculate RPKM for each bin for each sample
# (reads per kb per million mapped reads)
all_mag_mapping_RPKM <- all_mag_mapping %>%
  inner_join(RPM_samples, by="novogene_id_shotgun") %>%
  group_by(novogene_id_shotgun,genome) %>%
  dplyr::summarise(
    total_reads_bin = sum(Reads),
    total_length_bin_kb = (sum(Length)/1000),
    RPKM_bin = total_reads_bin / (num_reads_mapped_per_M*total_length_bin_kb)) %>%
  dplyr::select(novogene_id_shotgun,genome,RPKM_bin) %>%
  distinct()

# merge with sample and bin metadata
all_mag_mapping_RPKM %<>% full_join(metadata, by="novogene_id_shotgun") %>%
  inner_join(bin.tax, by="genome") %>%
  mutate(Genus = gsub("g__","",Genus))


## plot by age

#order genera by abundance
all_mag_mapping_RPKM$Genus <- factor(all_mag_mapping_RPKM$Genus,
                                     c("Snodgrassella","Schmidhempelia","Gilliamella",
                                       "Lactobacillus","Bifidobacterium","Bombilactobacillus",
                                       "Apilactobacillus"))

plot_mag_mapping_RPKM <- ggplot(all_mag_mapping_RPKM, aes(x=Age_days, y=RPKM_bin, fill=genome)) +
  geom_point(size=3, alpha=0.8, shape=21) +
  facet_wrap(~Genus,ncol=4,nrow=3, scales="free") +
  scale_fill_manual(values=c(
    "D_601_12_bin.3"="#a6611a",
    "D_601_12_bin.6"="#80cdc1",
    "D_530_05_bin.3"="#018571",
    "D_601_04_bin.12"="#a6611a",
    "D_602_12_bin.2"="#dfc27d",
    "D_604_07_bin.12"="#80cdc1",
    "D_601_12_bin.13"="#a6611a",
    "maxbin_D_603_02.003"="#a6611a",
    "D_601_05_bin.2"="#018571",
    "D_602_12_bin.10"="#dfc27d",
    "maxbin_D_605_04.005"="#a6611a",
    "D_601_01_bin.1"="#a6611a",
    "maxbin_D_530_08.002"="#dfc27d",
    "maxbin_D_604_09.001"="#018571",
    "maxbin_D_606_07.003"="#a6611a")) +
  geom_smooth(method="lm",se=T, color="black", size=0.5)+
  ylab("Genome abundance in sample (RPKM)") +
  xlab("Age (days)") +
  theme_bw()  +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic"))
plot_mag_mapping_RPKM





## RPKM stats for schmid and gill

rpkm_stats <- all_mag_mapping_RPKM %>% ungroup() %>% dplyr::select(novogene_id_shotgun,Age_days,Colony,genome,RPKM_bin)
rpkm_stats_maxbin_D_603_02.003 <- filter(rpkm_stats, genome=="maxbin_D_603_02.003") #gill
rpkm_stats_D_601_01_bin.1 <- filter(rpkm_stats, genome=="D_601_01_bin.1") #schmid


m1 <- lme(RPKM_bin ~ Age_days, random=~1|Colony,data=rpkm_stats_maxbin_D_603_02.003) 
summary(m1)
hist(resid(m1))
qqnorm(resid(m1)); qqline(resid(m1))

m2 <- lme(RPKM_bin ~ Age_days, random=~1|Colony,data=rpkm_stats_D_601_01_bin.1) 
summary(m2)
hist(resid(m2))
qqnorm(resid(m2)); qqline(resid(m2))

