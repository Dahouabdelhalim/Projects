#### Isotopic Sex Differences Meta-Analysis - Phylogenetic Signal
#### Joshua Bauld


## read in dataframe

remove(list = ls())


## clear console and plots

library(metafor)
library(broom)
library(ggtext)
library(ggbeeswarm)
library(ggpubr)
library(gt)
library(gtsummary)
library(flextable)
library(scales)
library(ape)
library(rotl)
library(phylosignal)
library(adephylo)
library(phylobase)
library(tidyverse)


## set working directory


## read in data

df_effects <- 

#### Nitrogen Mean Differences ####

df_meanN <- df_effects %>% drop_na(Mean_DiffN)

tree_meanN <- read.nexus("phylogeny_tree_meanN.nex")

## create df for tree tips

tips_meanN <- as_tibble(tree_meanN$tip.label) %>% rename(Species = value)

traits_meanN <- df_meanN %>% select(Species, Mean_DiffN) %>% distinct(Species, .keep_all =TRUE)

tip_traits_meanN <- left_join(tips_meanN, traits_meanN) %>% select(Mean_DiffN)

## add traits to the tree in a phylo4d object

p4d_meanN <- phylo4d(tree_meanN, tip_traits_meanN)

## plot tree alongside traits

barplot(p4d_meanN, tree.type = "phylo", tree.ladderize = TRUE)

## make correlogram for trait

crlg_meanN <- phyloCorrelogram(p4d_meanN, trait = "Mean_DiffN")

## plot both

plot(crlg_meanN,
     main="Phylogenetic Correlogram - Nitrogen Mean Sex Difference, All Species")

## very little global phylogenetic signal for nitrogen mean differences across all species in dataset


#### Carbon Mean Differences ####

df_meanC <- df_effects %>% drop_na(Mean_DiffC)

tree_meanC <- read.nexus("phylogeny_tree_meanC.nex")

## create df for tree tips

tips_meanC <- as_tibble(tree_meanC$tip.label) %>% rename(Species = value)

traits_meanC <- df_meanC %>% select(Species, Mean_DiffC) %>% distinct(Species, .keep_all =TRUE)

tip_traits_meanC <- left_join(tips_meanC, traits_meanC) %>% select(Mean_DiffC)

## add traits to the tree in a phylo4d object

p4d_meanC <- phylo4d(tree_meanC, tip_traits_meanC)

## plot tree alongside traits

barplot(p4d_meanC, tree.type = "phylo", tree.ladderize = TRUE)

## make correlogram for each trait

crlg_meanC <- phyloCorrelogram(p4d_meanC, trait = "Mean_DiffC")

## plot both

plot(crlg_meanC,
     main="Phylogenetic Correlogram - Carbon Mean Sex Difference, All Species")

## very little global phylogenetic signal for carbon mean differences across all species in dataset - very small positive phylogenetic signal at medium lags


#### Nitrogen Variation Differences ####

df_varN <- df_effects %>% drop_na(lnVRN)

tree_varN <- read.nexus("phylogeny_tree_varN.nex")

## create df for tree tips

tips_varN <- as_tibble(tree_varN$tip.label) %>% rename(Species = value)

traits_varN <- df_varN %>% select(Species, lnVRN) %>% distinct(Species, .keep_all =TRUE)

tip_traits_varN <- left_join(tips_varN, traits_varN) %>% select(lnVRN)

## add traits to the tree in a phylo4d object

p4d_varN <- phylo4d(tree_varN, tip_traits_varN)

## plot tree alongside traits

barplot(p4d_varN, tree.type = "phylo", tree.ladderize = TRUE)

## make correlogram for each trait

crlg_varN <- phyloCorrelogram(p4d_varN, trait = "lnVRN")

## plot both

plot(crlg_varN,
     main="Phylogenetic Correlogram - Nitrogen Variation Sex Difference, All Species")

## again very little autocorrelation - a negligible amount at medium lags 


#### Carbon Variation Differences ####

df_varC <- df_effects %>% drop_na(lnVRC)

tree_varC <- read.nexus("phylogeny_tree_varC.nex")

## create df for tree tips

tips_varC <- as_tibble(tree_varC$tip.label) %>% rename(Species = value)

traits_varC <- df_varC %>% select(Species, lnVRC) %>% distinct(Species, .keep_all =TRUE)

tip_traits_varC <- left_join(tips_varC, traits_varC) %>% select(lnVRC)

## add traits to the tree in a phylo4d object

p4d_varC <- phylo4d(tree_varC, tip_traits_varC)

## plot tree alongside traits

barplot(p4d_varC, tree.type = "phylo", tree.ladderize = TRUE)

## make correlogram for each trait

crlg_varC <- phyloCorrelogram(p4d_varC, trait = "lnVRC")

## plot both

plot(crlg_varC,
     main="Phylogenetic Correlogram - Carbon Variation Sex Difference, All Species")

## no significant autocorrelation at any lag distance


#### Non-Gape_Limited Carnivores, Nitrogen Mean Differences ####

df_carn <- df_effects %>% filter(Diet == "Carnivore") %>% 
  filter(Gape_Lim == "No") %>% drop_na(Mean_DiffN)

tree_carn <- read.nexus("phylogeny_tree_carn.nex")

## create df for tree tips

tips_carn <- as_tibble(tree_carn$tip.label) %>% rename(Species = value)

traits_carn <- df_carn %>% select(Species, Mean_DiffN) %>% distinct(Species, .keep_all =TRUE)

tip_traits_carn <- left_join(tips_carn, traits_carn) %>% select(Mean_DiffN)

## add traits to the tree in a phylo4d object

p4d_carn <- phylo4d(tree_carn, tip_traits_carn)

## plot tree alongside traits

barplot(p4d_carn, tree.type = "phylo", tree.ladderize = TRUE)

## make correlogram for each trait

crlg_carn <- phyloCorrelogram(p4d_carn, trait = "Mean_DiffN")

## plot both

plot(crlg_varC,
     main="Phylogenetic Correlogram - Nitrogen Mean Sex Difference, Non-Gape-Limited Carnivores")

## no significant autocorrelation at any lag distance


#### Gape-Limited Predators


df_gape <- df_effects %>% 
  filter(Diet == "Carnivore") %>% 
  filter(Gape_Lim == "Yes")

tree_gape <-  read.nexus("phylogeny_tree_gape.nex")


## create df for tree tips

tips <- as_tibble(tree_gape$tip.label) %>% rename(Species = value)

traits <- df_gape %>% select(Species, Mean_DiffN) %>% distinct(Species, .keep_all =TRUE)

tip_traits <- left_join(tips, traits) %>% select(Mean_DiffN)

## add traits to the tree in a phylo4d object

p4d_gape <- phylo4d(tree_gape, tip_traits)

## plot tree alongside traits

barplot(p4d_gape, tree.type = "phylo", tree.ladderize = TRUE)

## make correlogram for each trait

crlh_gape_meanN <- phyloCorrelogram(p4d_gape, trait = "Mean_DiffN")

## plot both

plot(crlh_gape_meanN,
     main="Phylogenetic Correlogram - Nitrogen Mean Sex Difference, Gape-Limited Carnivores")

## non-sig phylogenetic signal across the range of phylogenetic distances, but meta-regression still affected by accounting for phylogeny - check for local phylogenetic signal

## first test global phylogenetic signal

phyloSignal(p4d = p4d_gape, method = "all")

## no significant phylogenetic signal, for either trait, using any method - but these are a global measures

## now caluclate LIPA - local indicators of phylogenetic association

lipa_gape <- lipaMoran(p4d_gape)
lipa_gape_p4d <- lipaMoran(p4d_gape, as.p4d = TRUE)


## convert to dataframe to plot with ggplot, as easier to integrate with other plots later


df_lipa_gape <- as_tibble(as.data.frame(lipa_gape) %>% rownames_to_column() %>% 
  rename(species = rowname,
         lipa_score = Mean_DiffN,
         rank = Mean_DiffN.1,
         p_value = Mean_DiffN.2))

## add group column

df_lipa_gape <- df_lipa_gape %>% 
  mutate(group = 
           if_else(species == "Laticauda saintgironsi", "Snakes",
                   if_else(species == "Nerodia rhombifer", "Snakes",
                           if_else(species == "Nerodia erythrogaster", "Snakes",
                                   if_else(species == "Nerodia sipedon", "Snakes", "Fish")))))

## order factor levels for plotting

df_lipa_gape$group <- factor(df_lipa_gape$group, levels = c("Fish", "Snakes"))


df_lipa_gape %>% 
  ggplot(aes(x = fct_reorder(species, desc(lipa_score)), y = lipa_score)) +
  geom_hline(yintercept = 0, size = 1) +
  geom_hline(yintercept = -1, linetype = "dotted", size = 1) +
  geom_hline(yintercept = -0.5, linetype = "dotted", size = 1) +
  geom_hline(yintercept = 0.5, linetype = "dotted", size = 1) +
  theme_classic() +
  geom_point(aes(col = group, fill = group), 
             size = 7, shape = 21) +
  scale_color_manual(values = c("darkblue", "firebrick3"),
                     labels = c("Fish", "Snakes")) +
  scale_fill_manual(values = c("darkblue", "firebrick3"),
                     labels = c("Fish", "Snakes")) +
  coord_flip() +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 15, face = "italic"),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20, face = "bold"),
        legend.position = "left") +
  guides(col = "none") +
  labs(y = "Phylogenetic Signal",
       x = "Species",
       fill = "Group") 

## save df_lipa_gape, to load df and re-run this plot code in main data-analysis script, to combine this plot with the gape-limitation model

write_csv(df_lipa_gape, "df_lipa_gape.csv")
