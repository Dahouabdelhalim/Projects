# Delta introgression analysis for Labridae Phylogenomics Project
# Code modified from https://github.com/tauanajc/Cunha_Reimer_Giribet_2021_SystBio/tree/main/Supplementary%20Code%20S4
# See also Vanderpool et al. 2020 for more on the Delta statistic

setwd("~/Dropbox/FISHLIFE/Westneat_Lab/Labridae_ExonCap/introgression/introgression_cf/")

# Load required libraries

library(ggplot2)
library(tidyverse)
library(data.table)



# Import observed concordance factors. Legend to column names are at the bottom of this file.

d_NT = read.delim("observed_NT/observed.cf.stat",
                   header = T, comment.char = '#')
d_AA = read.delim("observed_AA/observed.cf.stat",
                    header = T, comment.char = '#')

d_NT$tree = "Topology 1" # M1-IQTREE-NT
d_AA$tree = "Topology 2" # M1-IQTREE-AA

d = rbind(d_NT, d_AA)

# Import concordance factors calculated from resampling gene trees (with replacement) 2000 times:

read_plus <- function(flnm){ # function to save filename in a new column
  fread(flnm) %>% 
    mutate(filename = basename(flnm))}

resampled_NT = list.files(path = "resample/resample_NT/",
                           pattern = "*.cf.stat", full.names = T) %>% 
  map_df(~read_plus(.)) # use new function instead of fread

resampled_AA = list.files(path = "resample/resample_AA/",
                          pattern = "*.cf.stat", full.names = T) %>% 
  map_df(~read_plus(.)) # use new function instead of fread

resampled_NT$tree = "Topology 1" # All-IQTREE-NT
resampled_AA$tree = "Topology 2" # All-IQTREE-AA

resampled = rbind(resampled_NT, resampled_AA)

# Calculate delta statistic on observed data

observed_delta = d %>%
  rowwise() %>%
  mutate(prop_discord_trees = sum(gDF1, gDF2),
         delta = if_else((gDF1_N + gDF2_N)==0, 0,
                         abs((gDF1_N - gDF2_N)/(gDF1_N + gDF2_N))))

# Define nodes of interest with more than 5% discordant gene trees
# A lot of discordant gene trees can indicate a number of things, ILS, introgression, error in the gene trees etc

most_discord = observed_delta %>%
  filter(prop_discord_trees > 5)
most_discord

# Z-score resampling

# Calculate delta score on resampled concordance factors:

resampled_delta = resampled %>%
  rowwise() %>%
  mutate(delta = if_else((gDF1_N + gDF2_N)==0, 0,
                         abs((gDF1_N - gDF2_N)/(gDF1_N + gDF2_N))))

# Calculate mean and sd of delta for each node based on the 2000 resamples (a null distribution):

nulldist = resampled_delta %>%
  group_by(tree, ID) %>% # group by node/resample
  summarize(mean = mean(delta), sd = sd(delta)) # mean and sd for all resamples of that node
nulldist %>%
  group_split

# For each branch, calculate a standardized z-score of observed deltas relative to the null 
# distribution created with the resampling, then check probability of such a score and calculate 
# the pvalue. One-tailed test looking for CDF_P>0.95 (pvalue<0.05) as evidence of introgression

test = most_discord %>%
  left_join(nulldist) %>%
  mutate(zscore = (delta - mean)/sd,
         CDF_prob = pnorm(zscore),
         pvalue = 1-CDF_prob)

test %>%
  select(-starts_with(c("g","s","m","L"))) %>%
  group_by(tree) %>%
  arrange(desc(CDF_prob)) %>%
  group_split()



