#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################

### Predictors of pairwise similarity 
## 1. Dispersal History
## 2. Habitat Overlap/Proximity
## 3. Relatedness

rm(list=ls())
graphics.off()

library(igraph)
library(vegan)
library(agricolae)
library(coin)

focal_dir <- "~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/"
setwd(focal_dir)
getwd()

load(file="data/sifaka_GLMM_covariates_duration.RData") #beta_div
mapping_file <- read.table("data/Perofsky_mapping_file_sifaka.txt", header=T)
head(mapping_file)
mapping_file <- mapping_file[,-c(21,22)]
write.table(mapping_file, "data/Perofsky_mapping_file_sifaka.txt", sep="\\t")

##############################################
## Dispersal History
##############################################
## limit analysis to individuals within the same group
beta_div_res <- beta_div[beta_div$Group =="Same"&beta_div$Pairwise_Res != "Both-Immigrants",]
beta_div_res <- droplevels(beta_div_res)
table(beta_div_res$Pairwise_Res)
nrow(beta_div_res)#97


## Two-sided approximative Wilcoxon-Mann-Whitney test
wilcox_test(Bray_Curtis ~ Pairwise_Res, data = beta_div_res,
            alternative = "two.sided",
            distribution = approximate(B = 10000))

wilcox_test(W_Unifrac ~ Pairwise_Res, data = beta_div_res,
            alternative = "two.sided",
            distribution = approximate(B = 10000))

##############################################
## Habitat Overlap/Proximity
##############################################
prox2 <- independence_test(Bray_Curtis ~ Proximity,
                           data = beta_div,
                           distribution = approximate(B = 10000),
                           ytrafo = function(data)
                             trafo(data, numeric_trafo = rank_trafo)) 
print(pvalue(prox2, method="single-step"))
pvalue(prox2)

##############################################
## Relatedness/Vertical Inheritance
##############################################

## population-level (without looking at social group affiliation)
mat2 <- independence_test(Bray_Curtis ~ Maternal_Line,
                          data = beta_div,
                          distribution = approximate(B = 10000),
                          ytrafo = function(data)
                            trafo(data, numeric_trafo = rank_trafo)) 
print(pvalue(mat2, method="single-step"))
pvalue(mat2)

## within social groups
beta_div2 <- beta_div[beta_div$Group=="Same",]
mat_within <- independence_test(Bray_Curtis ~ Maternal_Line,
                                data = beta_div2,
                                distribution = approximate(B = 10000),
                                ytrafo = function(data)
                                  trafo(data, numeric_trafo = rank_trafo)) 
print(pvalue(mat_within, method="single-step"))
pvalue(mat_within)