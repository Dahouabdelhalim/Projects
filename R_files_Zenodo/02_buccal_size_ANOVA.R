# phylogenetic analysis of whether buccal cavity size depends on 
# winnowing, mouthbrooding, and their interaction

# bookkeeping ####
# libraries
library(phytools)
library(geomorph)
library(tidyverse)

# read in data:
buccal_df <- read.csv("data/buccal_areas.csv")
buccal_df <- buccal_df[which(buccal_df$region == "neotropics"), ]

# read in tree:
tree <- read.tree("data/McGee2020_tree.tre")

# log transform the numeric variables:
# index of numeric variables
num_idx <- which(sapply(buccal_df, is.numeric))
buccal_log <- as.data.frame(sapply(buccal_df[ , num_idx],
                                   log10))
buccal_log <- cbind(buccal_log, buccal_df[ , -num_idx])

# individual body size correction
# we're squaring standard length because area is a squared value
buccal_log$buccal_residuals <- lm(buccal_log$buccal_size ~ buccal_log$standard_length^2)$residuals
buccal_log$inxn <- paste(buccal_log$feeding, 
                         buccal_log$reproduction, sep = " & ")

# get it at the species level
buccal_sp <- buccal_log %>%
  group_by(species) %>%
  summarise(feeding = unique(feeding),
            mcgee_label = unique(mcgee_label),
            reproduction = unique(reproduction),
            buccal_proportion = mean(buccal_proportion),
            buccal_size = mean(buccal_size),
            buccal_residuals = mean(buccal_residuals),
            standard_length = mean(standard_length),
            head_size = mean(head_size),
            inxn = unique(inxn))

# tree trimming
namecheck <- geiger::name.check(tree, buccal_sp, data.names = buccal_sp$mcgee_label)
trimmed_tree <- drop.tip(tree, tip = namecheck$tree_not_data)

# reorder
buccal_sp <- buccal_sp[match(trimmed_tree$tip.label,
                             buccal_sp$mcgee_label), ]

# phylogenetic body size correction ####

# get indices for size correction
n_idx <- which(unlist(lapply(buccal_log, is.numeric)))
x_idx <- n_idx[grep("standard_length", names(n_idx))]
y_idx <- n_idx[-x_idx]

# make matrices/vectors
x <- setNames(buccal_log[ , x_idx], buccal_log$mcgee_label)
Y <- as.matrix(buccal_log[ , y_idx])
rownames(Y) <- names(x)

# get phylogenetic residuals
# note -- this averages individuals together!
resid_df <- phyl.resid(trimmed_tree, x = x^2, Y = Y)

# reorder
tip_order <- match(trimmed_tree$tip.label, rownames(resid_df$resid))
resid_df <- resid_df$resid[tip_order, ]

# and add in discrete traits again:
sp_idx <- match(rownames(resid_df), buccal_log$mcgee_label)
resid_df <- cbind(resid_df, buccal_log[sp_idx, -n_idx])

# phylogenetic ANOVA (Adams & Collyer 2018) ####

# make the geomorph dataframe with raw values:
gmdf_log <- geomorph.data.frame(buccal_sp, phy = trimmed_tree)

# and with phylogenetic residuals:
gmdf_phyl_resid <- geomorph.data.frame(resid_df, phy = trimmed_tree)

# size correction method 1: use SL as a variable
aov_1 <- procD.pgls(data = gmdf_log,
                    buccal_size ~ 
                      standard_length^2 + 
                      reproduction * feeding,
                    phy = phy, iter = 9999)
summary(aov_1)

# method 2: use the buccal residuals (already size-corrected)
aov_2 <- procD.pgls(data = gmdf_log,
                    buccal_residuals ~ 
                      reproduction * feeding,
                    phy = phy, iter = 9999)
summary(aov_2)

# not really qualitatively different; I prefer including SL as a covariate
summary(aov_1) # size, length as a covariate
summary(aov_2) # non-phylo-size-corrected residuals

# reproduction & the interaction term are significant when we 
# fail to account for relatedness, because we essentially
# quadruple the number of "independent" mouthbrooders...
aov_3 <- procD.lm(data = gmdf_log,
                  buccal_size ~ standard_length^2 + 
                    reproduction * feeding,
                  iter = 9999)
summary(aov_3) 

# the effect is not present for head size:
aov_ctrl <- procD.pgls(data = gmdf_log,
                       head_size ~ standard_length^2 + reproduction * feeding,
                       phy = phy, iter = 9999)
summary(aov_ctrl)

# even without phylogenetic weighting:
summary(procD.lm(data = gmdf_log,
                 head_size ~ standard_length^2 + reproduction * feeding,
                 iter = 9999))
