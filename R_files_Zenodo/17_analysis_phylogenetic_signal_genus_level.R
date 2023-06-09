library(ape)
library(tidyverse)
library(LCVP)
library(picante)
library(phytools)
library(phylowood)

#load(file = "output/phylowood_sp_phyl_sig")

# get IW/genus data
# load species tree
tr <- read.tree("input/smith_and_brown/GBMB.tre")

# extract species/genus dataframe for replication
sp <- data.frame(species = tr$tip.label)
sp$genus <- str_split_fixed(sp$species, n = 2, pattern = "_")[,1]

## the number of IW species per genus
iw_species <- phylowood_sp %>% 
  filter(geo_insularwoody) %>% 
  group_by(tax_genus) %>% 
  count() %>% 
  rename(genus = tax_genus, iw_species = n)

# use the smith and brown total tree to get an approximation of the number of species per genus

splist <- data.frame(species = read.tree("input/smith_and_brown/ALLMB.tre")$tip.label)
splist$genus <- str_split_fixed(splist$species, n = 2, pattern = "_")[,1]

total_species <- splist %>% 
  group_by(genus) %>% 
  count() %>% 
  rename(total_species = n)

dat <- total_species %>% 
  left_join(iw_species)%>% 
  replace_na(list(iw_species = 0)) %>% 
  mutate(fraction_iw = round(iw_species / total_species, 5)) %>% 
  mutate(binary_iw = ifelse(iw_species > 0, 1, 0))

# filter to only those genera that are in the genetic tree
dat <- dat %>% filter(genus %in% sp$genus)

# replicate the analyses over thousand trees, for each tree one species is randomly sampled as the genus tip
# parallelize on cluster, takes a lot of time
for(i in 1:1000){
 
   print(i)
  # select the sub tree
  sub <- sp %>% group_by(genus) %>% sample_n(1)
  tr_sub <- keep.tip(tr, sub$species)
  tr_sub$tip.label <- str_split_fixed(tr_sub$tip.label, pattern = "_", n =2 )[,1]
 
  # Categorical phylogenetic distance
  comm <- dat %>% 
    mutate(iw = binary_iw) %>% 
    mutate(no_iw = as.numeric(!as.logical(binary_iw))) %>% 
    dplyr::select(genus, iw, no_iw) %>% 
    t()
  
  nam <- comm[1,]
  comm <- as.matrix(comm[-1,])
  comm <- apply(comm,2, "as.numeric")
  comm <- as.data.frame(comm)
  names(comm) <-  nam
  # 
  # phydist <- cophenetic(tr_sub)
  # ses.mpd <- ses.mpd(samp = comm, 
  #                    dis = phydist,
  #                    null.model = "taxa.labels",
  #                    abundance.weighted=FALSE, 
  #                    runs=99)
  # 
  # ses.mntd <- ses.mntd(samp = comm, 
  #                      dis = phydist,
  #                      null.model = "taxa.labels",
  #                      abundance.weighted=FALSE, 
  #                      runs=99)
  # 
  # 
  # # Pagels lmabda and blombergs K
  test_dat <- dat$fraction_iw
  names(test_dat) <- dat$genus
  test_dat <- test_dat[tr_sub$tip.label]
  
  ## Blombergs K
  BBK <- phylosig(tree = tr_sub,
                  x = test_dat,
                  test = TRUE)
  
  if(i == 1){
    write_csv(data.frame(BBK$K, BBK$P), file = "output/genus_level_blomberg_K.csv")
  }else{
    write_csv(data.frame(BBK$K, BBK$P), file = "output/genus_level_blomberg_K.csv", append = TRUE)
  }
  
  # Pagels lambda
  lambda <- phylosig(tree = tr,
                     x = test_dat,
                     method = "lambda",
                     test = TRUE)
  
  if(i == 1){
    write_csv(data.frame(lambda$lambda, lambda$logL, lambda$logL0, lambda$P), file = "output/genus_level_pagels_lambda.csv")
  }else{
    write_csv(data.frame(lambda$lambda, lambda$logL, lambda$logL0, lambda$P), file = "output/genus_level_pagels_lambda.csv", append = TRUE)
  }
}
