# Perofsky et al. "Social groups constrain the spatiotemporal dynamics of wild sifaka gut microbiomes"
# Contact: Amanda Perofsky, amanda.perofsky@nih.gov

rm(list=ls())
graphics.off()

library(phyloseq)
library(vegan)
library(dplyr)

## set your working directory to where files are located
## input files are in "Rdata" folder
## figures are saved in "figures" folder
dir <- "~/OneDrive - National Institutes of Health/NIH_Laptop_Updates_Post_Damage/Documents/Sifaka_KMNP_2016/Final_Code/Cleaned_Scripts/" ## set your working directory to where files are located
## set your working directory to where files are located
setwd(dir)

## import phyloseq object
## normalized sequence counts, limited to ASVs shared across at least two samples 
## marked individuals only 
load("Rdata/sifaka_trimmed_phyloseq_normalized_dada2.RData") #kmnp_trimmed (beta diversity, shared OTUs)
kmnp_trimmed
head(sample_data(kmnp_trimmed))
#######################################################################
## PERMANOVA
#######################################################################
sample_data(kmnp_trimmed)
unique(sample_data(kmnp_trimmed)$Group)

## limit to marked individuals in Groups I to VI
kmnp_trimmed %>%
  subset_samples(Group %in% c("I","II","III","IV","V","VI") & Age != "unknown" & Sex != "unknown") -> subset_kmnp
phyloseq::nsamples(subset_kmnp)#305 samples
df_beta = data.frame(sample_data(subset_kmnp))

#######################################################################
## 2012 PERMANOVA: two extraction methods
#######################################################################
#35 individuals
set.seed(7*11*13) 
R2scores <- list()
Pvalues <- list()
for(i in 1:100){
  subset_kmnp_df = df_beta %>% filter(Year=="2012" & Age != "unknown" & Sex != "unknown") %>% group_by(Name) %>% sample_n(1) %>% ungroup()
  subset_kmnp_red = subset_kmnp %>% subset_samples(X.SampleID %in% subset_kmnp_df$X.SampleID)
  subset_kmnp_red  <- prune_taxa(taxa_sums(subset_kmnp_red )>0, subset_kmnp_red )
  bdist <- phyloseq::distance(subset_kmnp_red,"bray")
  ad1 <- adonis2(bdist~ Age + Sex + Group + Extraction, by="margin",as(sample_data(subset_kmnp_red),"data.frame"), permutations=1000)
  R2scores[[i]] <- ad1$R2[1:4]
  Pvalues[[i]] <- ad1$`Pr(>F)`[1:4]
}
adonis_output_2012 = as.data.frame(do.call(rbind, R2scores))
names(adonis_output_2012) <- c("Age","Sex","Group","Extraction")
apply(adonis_output_2012,2,mean)
#     Age        Sex      Group   Extraction 
#0.03466951 0.02593882 0.40224210 0.04131905 

## p-values
adonis_output_2012 = as.data.frame(do.call(rbind, Pvalues))
names(adonis_output_2012) <- c("Age","Sex","Group","Extraction")
apply(adonis_output_2012,2,mean)
#       Age         Sex       Group  Extraction 
# 0.546903097 0.212697303 0.000999001 0.040059940 
#######################################################################
## Multi-year PERMANOVA
#######################################################################
## limit to marked individuals sampled in at least 2 years
## stratify by extraction method
#######################################################################
# Year + Group + Age + Sex 
#######################################################################
set.seed(7*11*13) 
R2scores <- list()
Pvalues <- list()
for(i in 1:100){
  subset_kmnp_df = df_beta %>% 
    group_by(Name) %>% filter(n()>1) %>% ungroup %>% group_by(Year,Name) %>% 
    sample_n(1) %>% ungroup() %>% group_by(Name) %>% filter(n()>1)
  subset_kmnp_red = subset_kmnp %>% subset_samples(X.SampleID %in% subset_kmnp_df$X.SampleID)
  subset_kmnp_red  <- prune_taxa(taxa_sums(subset_kmnp_red )>0, subset_kmnp_red )
  bdist <- phyloseq::distance(subset_kmnp_red,"bray")
  perm <- how(nperm = 1000)
  dat = as(sample_data(subset_kmnp_red),"data.frame")
  setBlocks(perm) <- with(dat, Extraction)
  ad1 <- adonis2(bdist~ Year + Group + Age + Sex, by="margin",dat, permutations=perm)
  R2scores[[i]] <- ad1$R2[1:4]
  Pvalues[[i]] <- ad1$`Pr(>F)`[1:4]
}
adonis_output = as.data.frame(do.call(rbind, R2scores))
names(adonis_output) <- c("Year","Group","Age","Sex")
apply(adonis_output,2,mean)
#       Year      Group        Age        Sex 
# 0.06006647 0.23961967 0.02429700 0.02834390 
# p values
adonis_output2 = as.data.frame(do.call(rbind, Pvalues))
names(adonis_output2) <-  c("Year","Group","Age","Sex")
apply(adonis_output2,2,mean,na.rm=T)
# Year       Group         Age         Sex 
#0.001508492 0.000999001 0.354835165 0.013486513 

#######################################################################
# Year + ID
#######################################################################
set.seed(7*11*13) 
R2scores <- list()
Pvalues <- list()
for(i in 1:100){
  subset_kmnp_df = df_beta %>% 
    group_by(Name) %>% filter(n()>1) %>% ungroup %>% group_by(Year,Name) %>% 
    sample_n(1) %>% ungroup() %>% group_by(Name) %>% filter(n()>1)
  subset_kmnp_red = subset_kmnp %>% subset_samples(X.SampleID %in% subset_kmnp_df$X.SampleID)
  subset_kmnp_red  <- prune_taxa(taxa_sums(subset_kmnp_red )>0, subset_kmnp_red )
  bdist <- phyloseq::distance(subset_kmnp_red,"bray")
  perm <- how(nperm = 1000)
  dat = as(sample_data(subset_kmnp_red),"data.frame")
  setBlocks(perm) <- with(dat, Extraction)
  ad1 <- adonis2(bdist~ Name + Year, by="margin",dat, permutations=perm)
  R2scores[[i]] <- ad1$R2[1:2]
  Pvalues[[i]] <- ad1$`Pr(>F)`[1:2]
}
adonis_output = as.data.frame(do.call(rbind, R2scores))
names(adonis_output)<- c("Name","Year")
apply(adonis_output,2,mean)
# Name       Year 
#0.52885396 0.05139409

#p-value
adonis_output = as.data.frame(do.call(rbind, Pvalues))
names(adonis_output) <-  c("Name","Year")
apply(adonis_output,2,mean,na.rm=T)
# Name        Year         
#0.000999001 0.002097902 
#######################################################################
# Year + Age + Sex 
#######################################################################
set.seed(7*11*13) 
R2scores <- list()
Pvalues <- list()
for(i in 1:100){
  subset_kmnp_df = df_beta %>% 
    group_by(Name) %>% filter(n()>1) %>% ungroup %>% group_by(Year,Name) %>% 
    sample_n(1) %>% ungroup() %>% group_by(Name) %>% filter(n()>1)
  subset_kmnp_red = subset_kmnp %>% subset_samples(X.SampleID %in% subset_kmnp_df$X.SampleID)
  subset_kmnp_red  <- prune_taxa(taxa_sums(subset_kmnp_red )>0, subset_kmnp_red )
  bdist <- phyloseq::distance(subset_kmnp_red,"bray")
  perm <- how(nperm = 1000)
  dat = as(sample_data(subset_kmnp_red),"data.frame")
  setBlocks(perm) <- with(dat, Extraction)
  ad1 <- adonis2(bdist~ Age + Sex + Year, by="margin",dat, permutations=perm)
  R2scores[[i]] <- ad1$R2[1:3]
  Pvalues[[i]] <- ad1$`Pr(>F)`[1:3]
}
adonis_output = as.data.frame(do.call(rbind, R2scores))
names(adonis_output)<- c("Age","Sex","Year")
apply(adonis_output,2,mean)
# Age        Sex       Year        
#0.03518276 0.02368082 0.06626752 
#p-value
adonis_output = as.data.frame(do.call(rbind, Pvalues))
names(adonis_output) <-  c("Age","Sex", "Year")
apply(adonis_output,2,mean,na.rm=T)
#       Age         Sex        Year  
#0.213516484 0.096093906 0.004355644 