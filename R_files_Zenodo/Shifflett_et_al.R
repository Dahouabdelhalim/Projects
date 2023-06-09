#### Supplementary R code for Shifflett et al. 2023

###################
## Load packages ##
###################
library(dplyr)
library(readr) # reading samtools files
library(readxl) # use to read in excel files 
library(writexl) # use to save data frames as excel spreedsheets on computer
library(tibble)
library(vegan)
library(ggplot2)
library(fields) # use for rdist.earth for geographic distance calculation
library(otuSummary) #use to turn distance matrix into a data frame
library(tidyr)
library(forcats)
library(plotrix)
library(data.table)
library(ggrepel)

####################################################
## processing ospC mapping stats (from SAMtools) ##
####################################################
## set working directory
setwd("~/Desktop/Shifflett et al. 2023 Dryad/ospC mapping statistics")

## see files
list.files() # all files in the folder
(dat_files <- list.files()[grep("_coverage_samtools.txt", list.files())]) # files for analysis

## see sample names
(ids <- sub("_coverage_samtools.txt", "", dat_files))

## data processing for all files using a for-loop
dat <- list() # empty list to save data

for(i in seq_along(dat_files)){
  
  ## load data
  file.i <- dat_files[i] # select file name
  
  dat.i <- read_tsv(file.i, col_names = FALSE, show_col_types = FALSE) # read file
  
  ## status update
  cat(paste0("Working on file ", file.i))
  
  ## process data
  dat.i.n <- dat.i %>%
    filter(X2 >200 & X2<550) %>%
    group_by(X1) %>%
    summarise(mean_cov_middle = mean(X3), sd_cov_middle = sd(X3)) %>%
    mutate(rel_cov_middle = (mean_cov_middle/sum(mean_cov_middle))) %>%
    arrange(desc(rel_cov_middle))
  
  (dat.i.nn <-  dat.i %>%
      group_by(X1) %>%
      summarize(depth_mean = mean(X3),
                depth_sd = sd(X3),
                n_5X = sum(ifelse(X3 >= 5, 1, 0)),
                ref_len = length(unique(X2))) %>%
      mutate(prop_5X = n_5X / ref_len) %>%
      arrange(desc(depth_mean)) %>%
      left_join(dat.i.n, by = "X1") %>%
      rename("Ref_allele" = "X1") %>%
      mutate(Sample_ID = ids[i]))
  
  ## save to list
  dat[[i]] <- dat.i.nn
  
}

## convert to data frame
dat.n <- do.call("rbind", dat)


################################################
## calculate and test prevalence across sites ##
################################################
## set wd
setwd("~/Desktop/Shifflett et al. 2023 Dryad/")

## read prevalence PCR and qPCR data
prev <- read_excel("Shifflett et al. 2023 Dryad Metadata.xlsx", sheet = "PCR vs qPCR") 

prev[prev$Sample_ID == "IS-49",which(names(prev) %in% "PCR")] <- 0 #changed Is-49 from infected to uninfected to account for B. Myamotoi

## calculate infections and totals by site
prev.n <- prev %>%
  group_by(Site) %>%
  summarize(n = length(unique(Sample_ID)),
            inf_PCR = sum(PCR),
            inf_qPCR = sum(qPCR))

## PCR (all sites) - Fisher's exact test
prev.n %>%
  mutate(uninf_PCR = n - inf_PCR) %>%
  select(inf_PCR, uninf_PCR) %>%
  fisher.test() # P = 0.07022

## PCR (remove S2 and S3) - Fisher's exact test
prev.n %>%
  filter(! n < 10) %>%
  mutate(uninf_PCR = n - inf_PCR) %>%
  select(inf_PCR, uninf_PCR) %>%
  fisher.test() # P = 0.05969

## PCR (All sites) sum total number of infections 
sum(prev.n$inf_PCR) #28 B. burgdorferi infections 

## qPCR (all sites) - Fisher's exact test
prev.n %>%
  mutate(uninf_qPCR = n - inf_qPCR) %>%
  select(inf_qPCR, uninf_qPCR) %>%
  fisher.test() # P =0.1438

## qPCR (remove S2 and S3) - Fisher's exact test
prev.n %>%
  filter(! n < 10) %>%
  mutate(uninf_qPCR = n - inf_qPCR) %>%
  select(inf_qPCR, uninf_qPCR) %>%
  fisher.test() # 0.1361


## correlation PCR vs. qPCR
cor(prev$PCR, prev$qPCR) # r = 0.9123 


##############################
## choosing true infections ##
##############################

## read ospC allele data
al <- read_excel("Shifflett et al. 2023 Dryad Metadata.xlsx", sheet = "ospC allele table")

## select DBI lab id and true id for joining
al.ids <- al %>%
  select("DBI label", "ID")

## join true id
dat.nn <- dat.n %>%
  left_join(al.ids, by = c("Sample_ID" = "DBI label"))

## Filter to only include prop_5X =1, but retain the samples that didn't map
dat.n3 <- dat.nn %>%
  filter(prop_5X == 1) %>%
  right_join(al.ids, by = "ID") 


## recreate original "al" data frame 
al.n <- dat.n3 %>%
  select(Ref_allele, prop_5X, ID) %>%
  pivot_wider(id_cols = ID, 
              names_from = Ref_allele, 
              values_from = prop_5X,
              values_fill = 0) %>%
  select(-"NA") %>%
  mutate(allelic_multiplicity = rowSums(select_if(., is.numeric))) %>%
  rename("T" = "T_F128", "A" = "A_B31",
         "D" = "D_CA-11-2A", "K" = "K_297", "E" = "E_N40",
         "B" = "B_64b", "M" = "M_29805", "G" = "G_72a",
         "H" = "H_156a", "N" = "N_F004", "C" = "C_JD1",
         "F" = "F_F084", "I" = "I_WI91-23", "vsp" = "vsp_N030",
         "O" = "O_N045", "U" = "U_94a")

## Is-49 has vsp, but should get a zero for allelic multiplicity; fix here
al.n[al.n$ID == "Is-49", which(names(al.n) %in% "allelic_multiplicity")] <- 0

## add human infectious info to table
al.nn <- dat.n3 %>%
  select(Ref_allele, ID) %>%
  mutate(human_infectious_allele = ifelse(Ref_allele %in% c("A_B31",
                                                            "B",
                                                            "I_WI91-23",
                                                            "K_297"),
                                          1, 0)) %>%
  group_by(ID) %>%
  summarize(human_inf_allele_total = sum(human_infectious_allele)) %>%
  right_join(al.n, by = "ID") %>%
  relocate(human_inf_allele_total, .after = last_col()) %>%
  mutate(ordering = as.numeric(sub("[^-]*-", "", ID))) %>%
  mutate(order_group = sub("-[0-9]*", "", ID)) %>%
  arrange(order_group, ordering) %>%
  select(-order_group, -ordering)


##################################################################
## Bray-Curtis dissimilarity between sites relative allele freq ##
##################################################################

## calculate frequency of each allele at each site
al_site_freq <- al %>%
  select(ID, Locations) %>%
  right_join(dat.n3, by = "ID") %>%
  select(ID, Locations, Ref_allele) %>%
  filter(!is.na(Ref_allele)) %>%
  filter(!Ref_allele == "vsp_N030") %>%
  group_by(Ref_allele, Locations) %>%
  summarize(n_allele = length(Ref_allele))

## calculate relative frequency of each allele at each site
al_site_freq.rel <- al_site_freq %>%
  group_by(Locations) %>%
  summarize(total_alleles = sum(n_allele)) %>%
  right_join(al_site_freq, by = "Locations") %>%
  mutate(rel_allele = n_allele / total_alleles)

## Sites with >3 alleles
al_site_freq %>%
  group_by(Locations) %>%
  summarize(total_alleles = sum(n_allele)) # K1, K2, NC1, NC3


## Bray-Curtis dissimilarities relative allele freq > 3 alleles
(sites.bc.rel <- al_site_freq.rel %>%
    filter(total_alleles > 3) %>%
    pivot_wider(id_cols = Locations, 
                names_from = Ref_allele, 
                values_from = rel_allele,
                values_fill = 0) %>%
    column_to_rownames("Locations") %>%
    vegdist(method = "bray")) ## Supplementary table 6

## mean +/- se of previous dissimilarity matrix
mean(sites.bc.rel) # mean
sd(sites.bc.rel)/sqrt(length(sites.bc.rel)) # se


#######################################################################
## Geographic distances between sites vs. Bray Curtis dissimilarites ##
#######################################################################

                ## Geographic coordinates not included ##

## read geographic data
geog <- read_excel()

## select sites from Bray Curtis dissimilarity (3 or more alleles present)
geog.dist <- geog %>%
  filter(Code %in% attr(sites.bc.rel, "Labels")) %>%
  select(-Site) %>%
  select(Lon, Lat, Code) %>%
  column_to_rownames("Code") %>%
  rdist.earth(miles = FALSE) %>%
  as.dist(upper = FALSE)

## add attributes labels back in
attr(geog.dist, "Labels") <- attr(sites.bc.rel, "Labels")

## plot the two matrices
plot(geog.dist, sites.bc.rel) 

## correlation
mantel(geog.dist, sites.bc.rel) # Mantel r = 0.75, P= 0.042

##turn matrix into data.frame 
geog.dist.df <- matrixConvert(geog.dist, colname = c("site1","site2", "dist"))

sites.bc.rel.df <- matrixConvert(sites.bc.rel, colname = c("site1","site2", "dist"))

## merged values of both data frames 
dat.nn <- data.frame(geog.dist.df$site1,geog.dist.df$site2, 
                     "sites"= c("K1 & K2","K1 & NC1", "K2 & NC1","K1 & NC3",
                                "K2 & NC3", "NC1 & NC3"),
                     geog.dist.df$dist,sites.bc.rel.df$dist)


## FIGURE 1: Graph matrices using ggplot 
## Graph 
ggplot(dat.nn, aes(x=geog.dist.df.dist, y=sites.bc.rel.df.dist)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Geographic distance between sites (km)") +
  ylab("Bray-Curtis dissimilarity between sites") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(20, 90)) +
  geom_text_repel(size=3,aes(label= sites)) +
  theme_classic() +
  theme(axis.text.x = element_text(size =10, color="black"),
        axis.text.y = element_text(size =10, color="black"))


######################
## Mixed infections ##
######################

## mean number of mixed infections and standard error
MixedInfections <- al.nn %>%
  filter(allelic_multiplicity > 1) 

mean(MixedInfections$allelic_multiplicity) # 3

stderror <- function(x) sd(x)/sqrt(length(x))
stderror(MixedInfections$allelic_multiplicity) #0.3205726

## histogram of allelic multiplicity
al.nn %>%
  filter(allelic_multiplicity > 1) %>%
  ggplot(aes(x = allelic_multiplicity)) +
  geom_histogram(binwidth = 1, color = "white") +
  scale_x_continuous("Number of" ~italic(ospC)~ "alleles in a mixed infection",
                     breaks = 2:7) +
  ylab("Number of mixed infections") +
  theme_bw() #Figure 2 5x5 dimensions 


## calculate allele frequency in single infections only
(rel_sing <- al %>%
    select(ID, Locations) %>%
    right_join(al.nn, by = "ID") %>%
    filter(allelic_multiplicity == 1) %>%
    select(ID, Locations) %>%
    left_join(dat.n3, by = "ID") %>%
    group_by(Ref_allele) %>%
    summarize(n_allele_sing = length(Ref_allele)) %>%
    mutate(rel_n_allele_sing = n_allele_sing/sum(n_allele_sing)))


## calculate allele frequency in mixed infections only
(rel_mix <- al %>%
    select(ID, Locations) %>%
    right_join(al.nn, by = "ID") %>%
    filter(allelic_multiplicity > 1) %>%
    select(ID, Locations) %>%
    left_join(dat.n3, by = "ID") %>%
    group_by(Ref_allele) %>%
    summarize(n_allele_mix = length(Ref_allele)) %>%
    mutate(rel_n_allele_mix = n_allele_mix/sum(n_allele_mix)))


## join single and mixed data and rename alleles to abbreviated terms 
# Rename alleles 
rename_ospC <- data.frame("Ref_allele"= c("A_B31","B_64b","C_JD1","D_CA-11-2A",
                                          "E_N40","I_WI91-23","K_297","O_N045",
                                          "T_F128","U_94a","F_F084","G_72a",
                                          "H_156a","N_F004","M_29805"),
                          allele= c("A","B","C","D", "E","I","K","O","T","U",
                                    "F","G","H","N","M"))

(rel_freqs <- rel_sing %>%
    full_join(rel_mix, by = "Ref_allele") %>%
    replace_na(list(n_allele_sing = 0,
                    rel_n_allele_sing = 0,
                    n_allele_mix = 0,
                    rel_n_allele_mix = 0)) %>%
    left_join(rename_ospC, by="Ref_allele"))


## FIGURE 3
ggplot(rel_freqs, aes(x=rel_n_allele_sing, y=rel_n_allele_mix)) +
  geom_point() +
  geom_smooth(method=lm,se = FALSE) +
  xlab("Relative" ~italic(ospC)~ "allele frequency in single infections") +
  ylab("Relative" ~italic(ospC)~ "allele frequency in mixed infections") +
  geom_text_repel(aes(label=allele)) +
  theme_classic() +
  theme(axis.text.x = element_text(size =10, color="black"),
        axis.text.y = element_text(size =10, color="black"))# export PDF 5 by 6 in


## correlation
cor.test(rel_freqs$rel_n_allele_sing, rel_freqs$rel_n_allele_mix, 
         method = "spearman") # rho = 0.537, P = 0.039

cor.test(rel_freqs$rel_n_allele_sing, rel_freqs$rel_n_allele_mix, 
         method = "pearson") # r0.737, P = 0.0017, t=3.94


## mean number of mixed infections and standard error
MixedInfections <- al.nn %>%
  filter(allelic_multiplicity > 1) 

mean(MixedInfections$allelic_multiplicity) # 3

stderror <- function(x) sd(x)/sqrt(length(x))
stderror(MixedInfections$allelic_multiplicity) #0.3205726

######################
## Chao 1 estimator ##
######################

## check estimator at each site first
al %>%
  select(ID, Locations) %>%
  right_join(al.nn, by = "ID") %>%
  select(-ID, -vsp, -allelic_multiplicity, -human_inf_allele_total) %>%
  group_by(Locations) %>%
  summarize(across(U:O, sum)) %>%
  column_to_rownames("Locations") %>%
  estimateR()

## check estimator overall
al %>%
  select(ID, Locations) %>%
  right_join(al.nn, by = "ID") %>%
  select(-ID, -vsp, -allelic_multiplicity, -human_inf_allele_total, -Locations) %>%
  summarize(across(U:O, sum)) %>%
  estimateR() # Chao1 = 16 +/- 2.29, we observed 15


####################
## Create Table 2 ## 
####################
## Calculate unique number of alleles at each site
UniqueAlleles <- al %>%
  select(ID,Locations) %>%
  right_join(dat.n3, by = "ID") %>%
  group_by(Locations) %>%
  summarize(Unique_alleles = length(unique(Ref_allele)))

## Create table and sum allele freq across sites
I.scap<- al %>%
  select(ID,Locations) %>%
  right_join(al.nn, by = "ID") %>%
  select(-vsp,-ID,-human_inf_allele_total, -allelic_multiplicity) %>%
  group_by(Locations) %>%
  summarize(across(U:O, sum)) %>%
  mutate(Total_Alleles = rowSums(select_if(., is.numeric))) %>%
  left_join(UniqueAlleles, by ="Locations")

## Calculate allele frequencies across all sites and rename data to better join to table above 
AlleleFreq <- I.scap %>%
  summarize(across(U:O, sum))

AlleleFreq$Locations <- NA

AlleleFreq$Total_Alleles <- NA

AlleleFreq$Unique_alleles <- NA

AlleleFreq.n <- AlleleFreq%>%
  select(Locations,U:O,Total_Alleles,Unique_alleles)

## Table 2: Bind tables 
I.scap.n <- I.scap %>%
  bind_rows(AlleleFreq.n)


#################################
## create supplementary table 5##
#################################
SuppTable5 <- al %>%
  select(ID, Locations) %>%
  right_join(al.nn, by = "ID") %>%
  filter(allelic_multiplicity > 1) %>%
  select(-vsp,-human_inf_allele_total) %>%
  rename('Allele multiplicity' = allelic_multiplicity)


## Calculate allele frequencies across all sites and rename data to better join to table above 
AlleleFreq5 <- SuppTable5 %>%
  summarize(across(U:O, sum))

AlleleFreq5$Locations <- NA

AlleleFreq5$'Allele multiplicity' <- NA

AlleleFreq.5n <- AlleleFreq5%>%
  select(Locations,U:O,'Allele multiplicity')

## Supp Table 5: Bind tables 
SuppTable5.n <- SuppTable5 %>%
  bind_rows(AlleleFreq.5n)


######################################
## Supplementary Table 2 using VE2 ##
#####################################
Supp.Table2 <- dat.n %>%
  filter(Sample_ID== "VE2") %>%
  select(-depth_mean,-depth_sd,-n_5X,-ref_len,-Sample_ID) %>%
  rename('Mean depth of coverage' = mean_cov_middle) %>%
  rename('Standard deviation of depth of coverage' = sd_cov_middle) %>%
  rename('% of nucleotides with at least 5X depth of coverage' = prop_5X) %>%
  rename("ospC allele" = Ref_allele) %>%
  rename('Relative coverage(%)'=rel_cov_middle) %>%
  mutate(across(where(is.numeric), round, 2))

###############################################
## Use dat.n3 to create supplementary table 3##
###############################################
SuppTable3 <- dat.n3 %>%
  select(-depth_mean,-depth_sd,-n_5X,-ref_len,-Sample_ID,-`DBI label`) %>%
  rename('Mean depth of coverage' = mean_cov_middle) %>%
  rename('Standard deviation of depth of coverage' = sd_cov_middle) %>%
  rename('% of nucleotides with at least 5X depth of coverage' = prop_5X) %>%
  rename("ospC allele" = Ref_allele) %>%
  rename('Relative coverage'=rel_cov_middle) %>%
  mutate(across(where(is.numeric), round, 2)) 

