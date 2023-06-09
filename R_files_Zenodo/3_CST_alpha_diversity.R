# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

rm(list=ls())

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
setwd(dir)

library(phyloseq)
library(dplyr)
library(tidyr)

library(rcompanion)
library(rstatix) ## tidy stat tests

## plotting
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

## phyloseq object (kmnp_trimmed)
# normalized sequence counts, limited to ASVs shared across at least two samples 
## marked individuals only with sample CST assignments
load("Rdata/sifaka_trimmed_phyloseq_normalized_dada2_CSTs.RData")
head(sample_data(kmnp_trimmed))
#phyloseq object (kmnp2) 
# raw ASV counts including singletons; required for alpha diversity 
# marked individuals only
load("Rdata/sifaka_allotus_phyloseq_openref_2016.RData")
save(kmnp2,file="Rdata/sifaka_allotus_phyloseq_openref_2016.RData")
head(sample_data(kmnp2))
## estimating alpha diversity requires raw counts, including singletons

###################################################################
# Figure S15: CST alpha diversity
###################################################################
kmnp2 ### raw ASV counts for alpha diversity metrics
sample_data(kmnp2)
kmnp_trimmed

kmnp_cst = kmnp_trimmed ## CST assignments for each sample
kmnp_cst

cst_df = sample_data(kmnp_cst) %>% as_tibble() %>% dplyr::select(X.SampleID,Group,Year,CST) %>%arrange(Group,CST)

cst_ids = cst_df$X.SampleID

reduced_kmnp2 = kmnp2 %>% subset_samples(X.SampleID %in% cst_ids)
new_df = left_join(as_tibble(sample_data(reduced_kmnp2)),cst_df %>% dplyr::select(X.SampleID,CST),by="X.SampleID") 
rownames(new_df)<- new_df$X.SampleID
sample_data(reduced_kmnp2) = data.frame(new_df)

sample_data(reduced_kmnp2)$CST = as.factor(sample_data(reduced_kmnp2)$CST)
sample_data(reduced_kmnp2)$CST = factor(sample_data(reduced_kmnp2)$CST,levels=c("1","2","3","4","5","6","7"))


brewer = brewer.pal(n = 8, name = "RdBu")

richness_df = estimate_richness(reduced_kmnp2)
richness_df$X.SampleID = rownames(richness_df)

richness_df = left_join(richness_df,data.frame(sample_data(reduced_kmnp2)),by="X.SampleID")
richness_df_long = richness_df %>% pivot_longer(cols=Observed:Fisher,names_to="Measure",values_to="value")


# Pairwise comparisons
pwc <- richness_df_long %>% 
  group_by(Measure) %>%
  rstatix::dunn_test(value ~ CST, p.adjust.method = "BH") 
pwc
pwc$group_name = paste(pwc$group1,pwc$group2,sep="-")

pwc_group = group_by(pwc,Measure)
unique(pwc$Measure)
cld_fun <- function(sub_df){
  ret_df = cldList(p.adj ~ group_name,data = sub_df,threshold = 0.05)
  return(ret_df)
}
dunn_df = pwc %>%
  filter(Measure %in% c("Chao1","Simpson"))%>%
  group_by(Measure)%>%
  do(as.data.frame(cld_fun(.)))%>%
  ungroup()
names(dunn_df)[2]<-"CST"

max_values = richness_df_long %>%
  dplyr::group_by(Measure)%>%
  dplyr::summarize(max_y = max(value)+ (max(value)/50))

dunn_df = left_join(dunn_df,max_values,by="Measure")

brewer = brewer.pal(n = 8, name = "RdBu")
richness_df_long$CST = as.factor(richness_df_long$CST)
richness_df_long$CST  = factor(richness_df_long$CST ,levels=c("1","2","3","4","5","6","7"))

unique(richness_df_long$Measure)

richness_plot = ggplot(richness_df_long %>% filter(Measure %in% c("Chao1","Simpson")),aes(x=as.factor(CST),y=value,fill=CST))+
  geom_boxplot()+
  facet_wrap(~Measure,scales="free_y",labeller = as_labeller(c(Chao1 = "Chao1 Richness", Simpson= "Simpson's Index") ),
             strip.position = "left")+
  geom_text(data = dunn_df %>% filter(Measure %in% c("Chao1","Simpson")), aes(y = max_y, label = Letter), 
            position = position_dodge(width = .75),color="black")+
  scale_fill_manual(values = c(brewer[1],brewer[2],brewer[3],brewer[4],brewer[5],brewer[6],brewer[7]),
                    breaks=c("1","2","3","4","5","6","7"),
                    labels=c("1","2","3","4","5","6","7"))+
  xlab("CST")+
  theme_bw(base_size = 14)+
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.placement = "outside")+
  ylab(NULL) 

richness_plot
save_plot(richness_plot,filename = "figures/SFig15_CST_richness_plot.pdf",base_width = 8,base_height = 4)
