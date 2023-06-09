#phyloFlash analysis
library(tidyverse)

#load phyloflash output
pf_raw <- read.table("resubmission_files/phyloFlash_output.txt",header=T,sep='\\t')

#load metadata
metadata <- read.table("resubmission_files/sample_metadata_v2.txt",header=T,sep='\\t') %>%
  dplyr::rename(HostSpecies = Species)


#normalize counts to proportions using total SSU reads identified by phyloflash
pf_counts <- pf_raw %>%
  full_join(metadata, by="novogene_id_shotgun") %>% #merge with metadata
#create new taxon column
  mutate(taxon = case_when(Species == "Metazoa" ~ "Host",
                           Domain == "Bacteria" ~ "Bacteria",
                           Genus == "Fungi" ~ "Fungi",
                           Class == "Chloroplastida" ~ "Plant"
                           )) %>% #non-matching rows as labeled NA for taxon
  drop_na(taxon) %>%
  #filter out host here, to calculate taxa as prop. of identified*, non-host sequences
  # *note v. minor things not included in the 'taxon' column are dropped
  filter(taxon != "Host")

pf_total_nonhost_reads <- pf_counts  %>%
  group_by(novogene_id_shotgun) %>%
  dplyr::summarize(total_nonhost_reads=sum(num_reads))

pf_proportions <- pf_counts %>%
  inner_join(pf_total_nonhost_reads, by="novogene_id_shotgun") %>%
  mutate(prop_nonhost_reads=(num_reads/total_nonhost_reads)) %>%
#calculate total proportions for each non-host taxon
  group_by(SampleID, Age_days, Colony, taxon) %>%
  dplyr::summarize(total_prop_nonhost_reads = sum(prop_nonhost_reads))

#fill in missing 0's
pf_proportions_fill <- pf_proportions %>%
  ungroup() %>%
  dplyr::select(SampleID,taxon,total_prop_nonhost_reads) %>%
  tidyr::complete(SampleID,taxon, fill=list(total_prop_nonhost_reads=0)) %>%
  arrange(SampleID) %>%
  inner_join(metadata, by="SampleID")
  
# median proportions of nonhost reads
plant <- pf_proportions_fill %>% filter(taxon=="Plant")
median(plant$total_prop_nonhost_reads)
bacteria <- pf_proportions_fill %>% filter(taxon=="Bacteria")
median(bacteria$total_prop_nonhost_reads)
fungi <- pf_proportions_fill %>% filter(taxon=="Fungi")
median(fungi$total_prop_nonhost_reads)

#plot plant + fungi relative abundances by age
# (Fig. S7)

pf_proportions_fill %>%
  filter(taxon %in% c("Plant","Fungi"))  %>%
  ggplot(aes(x=Age_days,y=total_prop_nonhost_reads)) +
    geom_point(alpha=0.7, size=3) +
    facet_grid(rows = vars(Colony),
             cols = vars(taxon),
             scales = "free_y",
             labeller = label_both) +
    theme_bw() +
    ylab("Prop. in non-host SSU rRNA sequences") + xlab("Age (days)")
