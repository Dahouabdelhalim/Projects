## run Medusa method of estimating diversification rates

rm(list=ls())
setwd("~/Documents/HYDIV")
library(ape)
library(geiger)
library(dplyr)

set.seed(679)

# read in tree
v <- read.tree("data/full_tree.tre")

# read in family sizes
fam_use <- read.csv("data/FamSizes_use.csv", header = TRUE)

# run Medusa with estimated threshold
res <- medusa(phy = v, richness = fam_use, ncores = 32) # Appropriate  aicc-threshold for a tree of 459 tips is: 8.089111
zsum <- data.frame(res$zSummary)

# run Medusa with threshold = 2
res_2 <- medusa(phy = v, richness = fam_use, ncores = 32, threshold = 2)
zsum2 <- data.frame(res_2$zSummary)

# run Medusa with threshold = 4
res_4 <- medusa(phy = v, richness = fam_use, ncores = 32, threshold = 4)
zsum4 <- data.frame(res_4$zSummary)

# save summary infos, join them in microsoft excel sheets, keep only tip data, label according to tip order
write.csv(zsum, "data/medusa_zsum.csv", row.names = FALSE, na = ".")
write.csv(zsum_2, "data/medusa_zsum2.csv", row.names = FALSE, na = ".")
write.csv(zsum_4, "data/medusa_zsum4.csv", row.names = FALSE, na = ".")
