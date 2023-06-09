library(ape)
library(tidyverse)
library(LCVP)
library(picante)
library(phytools)

# set seed
set.seed(2505)

# load data
load(file = "output/number_of_IW_species_per_family.rda")
load(file = "output/eudicot_phylogeny.rda")

# get estimate of species number per family
dat <- LCVP::tab_lcvp %>% 
  # accepted only
  filter(Status == "accepted") %>% 
  # remove infra-specific ranks
  mutate(species = paste(str_split_fixed(Input.Taxon, pattern = " ", n = 3)[, 1],
                         str_split_fixed(Input.Taxon, pattern = " ", n = 3)[, 2])) %>% 
  dplyr::select(species, Family) %>%
  distinct() %>% 
  # count species per family
  group_by(Family) %>% 
  count() %>% 
  #merge in information on insular woodiness
  mutate(family = tolower(Family)) %>% 
  left_join(fams, by = "family") %>% 
  mutate(iw_species = ifelse(is.na(iw_species), 0, iw_species)) %>% 
  # filter to species occurring in the tree
  filter(family %in% tr$tip.label) %>% 
  # calcualte the fracction of IW species in each family
  ungroup() %>% 
  dplyr::select(family, total_species = n, iw_species) %>% 
  mutate(fraction_iw = round(iw_species / total_species, 5)) %>% 
  mutate(binary_iw = ifelse(iw_species > 0, 1, 0))

# Three families from the tree are missing. Why - NOt accepted according to LCVP
tr$tip.label[!tr$tip.label %in% dat$family]

# Exclude the missing species from the tree as well
tr <- keep.tip(tr, dat$family)

# Visualize the proportion the tree
plo <- dat$fraction_iw
names(plo) <- dat$family
plo <- data.frame(plo)
plo[plo == 0] <-  NA


p11 <- ggtree(tr, layout="circular", size = 0.1) +
  geom_tiplab(offset=35, size=2)


p11 <-  gheatmap(p11,
                 plo, 
                 width=0.2, 
                 low="blue",
                 high="red",
                 colnames = FALSE, 
                 offset = -12, 
                 font.size = 12)+
  #  scale_fill_viridis(option="D", name="Number of\\nDW\\nspecies", direction =-1)+
  scale_fill_viridis(option="D",
                     name = "Fraction of\\nIW\\nspecies",
                     direction =1,
                     na.value = "black")+
  theme(legend.key.height = unit(1, "cm"),
        legend.key.width = unit(2.5, "cm"),
        legend.text = element_text(size = 12),
        legend.position = "bottom")

ggsave(p11, file = "output/supplementary_figure_proportion_of_iw_phylogeny.jpg", width = 8)
  
# Categorical phylogenetic distance
comm <- dat %>% 
  mutate(iw = binary_iw) %>% 
  mutate(no_iw = as.numeric(!as.logical(binary_iw))) %>% 
  dplyr::select(family, iw, no_iw) %>% 
  t()

nam <- comm[1,]
comm <- as.matrix(comm[-1,])
comm <- apply(comm,2, "as.numeric")
comm <- as.data.frame(comm)
names(comm) <-  nam

phydist <- cophenetic(tr)
ses.mpd <- ses.mpd(samp = comm, 
                   dis = phydist,
                   null.model = "taxa.labels",
                   abundance.weighted=FALSE, 
                   runs=99)

ses.mntd <- ses.mntd(samp = comm, 
                   dis = phydist,
                   null.model = "taxa.labels",
                   abundance.weighted=FALSE, 
                   runs=99)

# 
# # Pagels lmabda and blombergs K
# test_dat <- dat$fraction_iw
# names(test_dat) <- dat$family
# test_dat <- test_dat[tr$tip.label]
# 
# ## Blombergs K
# BBK <- phylosig(tree = tr,
#          x = test_dat,
#          test = TRUE)
# 
# # PAgels lambda
# lambda <- phylosig(tree = tr,
#          x = test_dat,
#          method = "lambda",
#          test = TRUE)
# 
