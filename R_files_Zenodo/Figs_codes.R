
##############################################################
# Author: 
# Antica Culina
# NIOO-KNAW & University of Oxford 
# Email: a.culina@yahoo.com


##############################################################
# Description of script and Instructions
##############################################################

# This code is used to Plot the Figures of the main MS and Supplementary Figure S1 (time-lag bias)

#############################################################################################

# Set working directory to wherever you have downloaded the data and code
#session info  can be found at the end of this code


rm(list=ls())
require(ggplot2)
library(MCMCglmm)

sessionInfo()

before <- read.table(file.choose(), header = T, sep = ";")  # breeding success before divorce and occurrence of divorce MA

after<- read.table(file.choose(), header = T, sep = ";") # occurrence of divorce and breeding success after divorce and MA

FvsM <- read.table(file.choose(), header = T, sep = ";") # between sexes benefits and costs of partner change MA

#####################
# Figure 2, panel A)
#####################


library(psych)
fisherz2r() # for trnsforming Zr into r

# create table of effect sizes and estimates of the main random-effect models for each MA
x <- c("before", "after", "MvsF")
effect_size <- c(-0.079, -0.118, 0.033)
L <- c(-0.366, -0.299, -0.131)
U <- c(0.144, 0.105, 0.195)


set.seed(0815)
df <- data.frame(x, F, L, U)


ggplot(df, aes(x = x, y = effect_size)) + geom_point(size = 5) +
  geom_errorbar(aes(ymax = U, ymin = L), width = 0.3, size = 1) +
  geom_point(data = before, 
             mapping = aes(x=2, y = Z), color = "darkgreen", pch = 1, size=log(before$N), position = position_jitter(w = 0.06, h = 0)) +
  geom_point(data = after, 
                  mapping = aes(x = 1, y = Z), color = "darkblue", pch = 1, size=log(after$N), position = position_jitter(w = 0.06, h = 0)) +
  geom_point (data = FvsM,
                  mapping = aes(x=3, y = Z), color = "darkred", pch = 1, size=log(FvsM$N), position = position_jitter(w = 0.06, h = 0)) + coord_flip() 


# without color

ggplot(df, aes(x = x, y = effect_size)) + geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L), width = 0.3) +
  geom_point(data = before, 
             mapping = aes(x = 2, y = Z), pch = 1, size=log(before$N)) + 
  geom_point(data = after, 
             mapping = aes(x = 1, y = Z), pch = 1, size=log(after$N)) +
  geom_point (data = FvsM,
              mapping = aes(x=3, y = Z), pch = 1, size=log(FvsM$N)) +  coord_flip()



#################
### Fig 2 panel B
#################

## plotting 'before' dichotomisation


x <- c("non dichotomised", "dichotomised")
effect_size <- c(-0.285, 0.014)
L <- c(-0.573, -0.313)
U <- c(0.063, 0.314)

set.seed(0815)
df <- data.frame(x, F, L, U)

ggplot(df, aes(x = x, y = effect_size)) + geom_point(size = 4) + geom_errorbar(aes(ymax = U, ymin = L), width = 0.3)+ coord_flip()

# with data points for effects
before_1 <- subset(before, type_of_study=='obs')

ggplot(df, aes(x = x, y = effect_size)) + geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L), width = 0.3) +
  geom_point(data = subset(before, dichotomisation == "no"), 
             mapping = aes(x = 2, y = Z), color = "darkgreen", pch = 1, size=log(subset(before, dichotomisation == "no")$N)) +
  geom_point(data = subset(before, dichotomisation == "yes"), 
             mapping = aes(x = 1, y = Z), color = "darkgreen", pch = 1, size=log(subset(before, dichotomisation == "yes")$N)) +  coord_flip()

#################
### Fig 2 panel C
#################

## plotting 'after' after failure (note, this analyises was done an a subset of ES that had enough levels of the factor)

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


afterf <- completeFun(after, "is_renesting_after_failure")


x <- c("no", "yes")
effect_size <- c(-0.017, -0.155)
L <- c(-0.217, -0.315)
U <- c(0.198, 0.025)

set.seed(0815)
df <- data.frame(x, F, L, U)

ggplot(df, aes(x = x, y = effect_size)) + geom_point(size = 4) + geom_errorbar(aes(ymax = U, ymin = L), width = 0.3)+ coord_flip()

ggplot(df, aes(x = x, y = effect_size)) + geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L), width = 0.3) +
  geom_point(data = subset(after4, is_renesting_after_failure=='yes'), 
             mapping = aes(x = 2, y = Z), color = "darkblue", pch = 1, size=log(subset(after, is_renesting_after_failure=='yes')$N)) +
  geom_point(data = subset(after4, is_renesting_after_failure=='sometimes'), 
             mapping = aes(x = 1, y = Z), color = "darkblue", pch = 1, size=log(subset(after4, is_renesting_after_failure=='sometimes')$N)) +  coord_flip()


#===============================================================
# session info
#===============================================================

#R version 4.0.4 (2021-02-15)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 14393)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_Europe.1252  LC_CTYPE=English_Europe.1252    LC_MONETARY=English_Europe.1252
#[4] LC_NUMERIC=C                    LC_TIME=English_Europe.1252    

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggplot2_3.3.5  rmeta_3.0      metafor_2.4-0  phangorn_2.7.0 MCMCglmm_2.32  ape_5.5        coda_0.19-4   
#[8] Matrix_1.3-2  

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.6       compiler_4.0.4   pillar_1.6.1     tools_4.0.4      lifecycle_1.0.0  tibble_3.1.2    
#[7] gtable_0.3.0     nlme_3.1-152     lattice_0.20-41  pkgconfig_2.0.3  rlang_0.4.11     fastmatch_1.1-0 
#[13] igraph_1.2.6     DBI_1.1.1        parallel_4.0.4   withr_2.4.2      dplyr_1.0.6      generics_0.1.0  3[19] vctrs_0.3.8      tidyselect_1.1.1 grid_4.0.4       glue_1.4.2       R6_2.5.0         fansi_0.5.0     
#[25] tensorA_0.36.2   purrr_0.3.4      corpcor_1.6.9    magrittr_2.0.1   scales_1.1.1     codetools_0.2-18
#[31] ellipsis_0.3.2   assertthat_0.2.1 cubature_2.0.4.2 colorspace_2.0-1 quadprog_1.5-8   utf8_1.2.1      
#[37] munsell_0.5.0    crayon_1.4.1    
> 