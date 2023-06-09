# This script re-runs the comparative analyses on two variations on the original
# trait classifications to verify the robustness of the results:
# 
# 1. INCLUSIVE: species with ambiguous reports of mouthbrooding are
# classified as mouthbrooders.
#
# 2: EXCLUSIVE: Ambiguous species are excluded from the analysis entirely.

library(phytools)

# read in the tree:
tree <- read.tree("data/McGee2020_tree.tre")
tree <- rotateNodes(tree, c(1078, 1061, 939)) # for visual clarity

# Analyses are identical to those in the 01_trait_evolution.R script, but using
# different CSVs; see that script for code comments. The results of the analyses
# as run on the author's computer are included as RDS files already and loaded
# and summarized at the end of the script, but readers can recreate these analyses
# by uncommenting the below scripts.

# Inclusive analyses ####

# # read in data
# data <- read.csv("data/McGee2020_tips_inclusive.csv")
# reproduction <- setNames(data$reproduction, data$tip.label)
# feeding <- setNames(data$feeding, data$tip.label)
# inxn <- setNames(paste(reproduction, feeding, sep = " & "), data$tip.label)
# {inxn <- gsub("substrate brooding & other", "Neither behavior", inxn)
#   inxn <- gsub("substrate brooding & winnowing", "Winnowing only", inxn)
#   inxn <- gsub("mouthbrooding & other", "Mouthbrooding only", inxn)
#   inxn <- gsub("mouthbrooding & winnowing", "Both behaviors", inxn)}

# # pagel
# fit_dep <- fitPagel(tree, x =  reproduction, y = feeding)
# fit_ind <- fitPagel(tree, x = reproduction, y = feeding, dep.var = "x")
# fit_w_dep_m <- fitPagel(tree, x = reproduction, y = feeding, dep.var = "y")
# pagel_aic <- c(fit_dep$independent.AIC, # neither trait affects the other
#                fit_dep$dependent.AIC, # both traits depend on each other
#                fit_ind$dependent.AIC, # mouthbrooding depends on winnowing, but not the other way around
#                fit_w_dep_m$dependent.AIC) # winnowing depends on mouthbrooding, but not the other way around
# aic.w(pagel_aic)
# saveRDS("data/pagel_inclusive.rds")

# # Mk models
# mk_sym <- fitMk(tree, inxn)
# mk_ard <- fitMk(tree, inxn, model = "ARD")
# mk_er <- fitMk(tree, inxn, model = "ER")
# saveRDS(setNames(list(mk_ard, mk_er, mk_sym),
#                  c("ARD", "ER", "SYM")),
#                  file = "data/inclusive_mk_models.rds")

# # stochastic character maps
# nsim <- 1000
# simmap_ard <- make.simmap(tree, inxn, model = "ARD", nsim = nsim)
# simmap_er <- make.simmap(tree, inxn, model = "ER", nsim = nsim)
# simmap_sym <- make.simmap(tree, inxn, model = "SYM", nsim = nsim)
# # saveRDS(simmap_ard, file = "data/inclusive_simmap_ard.RDS")
# # saveRDS(simmap_er, file = "data/inclusive_simmap_er.RDS")
# # saveRDS(simmap_sym, file = "data/inclusive_simmap_sym.RDS")

# # weighted simmaps
# inclusive_mk_models <- readRDS("data/inclusive_mk_models.rds")
# aic_weights <- aic.w(unlist(lapply(mk_models, function(i) AIC(i))))
# nsum <- 1000
# sample_freq <- round(aic_weights * nsum)
# simmap_files <- c("data/inclusive_simmap_ard.RDS",
#                   "data/inclusve_simmap_er.RDS",
#                   "data/inclusive_simmap_sym.RDS")
# for (i in 1:length(simmap_files)) { 
#   simmap_obj <- readRDS(simmap_files[i])
#   if (i == 1) {
#     simmap_subsample <- simmap_obj[sample(1:length(simmap_obj),
#                                           sample_freq[i])]
#   } else {
#     simmap_subsample <- append(simmap_subsample,
#                                simmap_obj[sample(1:length(simmap_obj),
#                                                  sample_freq[i])])
#   }
#   rm(simmap_obj)
# }

# class(simmap_subsample) <- c("multiSimmap", "multiPhylo")
# saveRDS(simmap_subsample, file = "inclusive_simmap_weighted_avg.RDS")

# simmap_subsample <- readRDS("data/simmap_weighted_avg.RDS")
# sum_dmap <- density.multiSimmap(simmap_subsample)
# saveRDS(sum_dmap, "data/inclusive_sum_dmap.RDS)

# # q matrix
# inclusive_mk_models <- readRDS("data/inclusive_mk_models.rds")
# qmatrices <- lapply(inclusive_mk_models,
#                     as.Qmatrix)
# qmatrix_array <- abind::abind(qmatrices[[1]],
#                               qmatrices[[2]],
#                               qmatrices[[3]], 
#                               along = 3)
# qmatrix_avg <- qmatrices[[1]]
# for (i in 1:nrow(qmatrix_array)) {
#   for (j in 1:ncol(qmatrix_array)) {
#     qmatrix_avg[i, j] <- weighted.mean(qmatrix_array[i, j, ], aic_weights)
#   }
# }
# qmatrix_na <- qmatrix_avg
# for (i in 1:nrow(qmatrix_na)) {
#   qmatrix_na[i, i] <- NA
# }
# saveRDS(qmatrix_na, file = "data/inclusive_qmatrix.rds")

# Exclusive analyses ####

# # read in data
# data <- read.csv("data/McGee2020_tips_exclusive.csv")
# reproduction <- setNames(data$reproduction, data$tip.label)
# feeding <- setNames(data$feeding, data$tip.label)
# inxn <- setNames(paste(reproduction, feeding, sep = " & "), data$tip.label)
# {inxn <- gsub("substrate brooding & other", "Neither behavior", inxn)
#   inxn <- gsub("substrate brooding & winnowing", "Winnowing only", inxn)
#   inxn <- gsub("mouthbrooding & other", "Mouthbrooding only", inxn)
#   inxn <- gsub("mouthbrooding & winnowing", "Both behaviors", inxn)}

# # pagel
# fit_dep <- fitPagel(tree, x =  reproduction, y = feeding)
# fit_ind <- fitPagel(tree, x = reproduction, y = feeding, dep.var = "x")
# fit_w_dep_m <- fitPagel(tree, x = reproduction, y = feeding, dep.var = "y")
# pagel_aic <- c(fit_dep$independent.AIC, # neither trait affects the other
#                fit_dep$dependent.AIC, # both traits depend on each other
#                fit_ind$dependent.AIC, # mouthbrooding depends on winnowing, but not the other way around
#                fit_w_dep_m$dependent.AIC) # winnowing depends on mouthbrooding, but not the other way around
# aic.w(pagel_aic)
# saveRDS("data/pagel_exclusive.rds")

# # Mk models
# mk_sym <- fitMk(tree, inxn)
# mk_ard <- fitMk(tree, inxn, model = "ARD")
# mk_er <- fitMk(tree, inxn, model = "ER")
# saveRDS(setNames(list(mk_ard, mk_er, mk_sym),
#                  c("ARD", "ER", "SYM")),
#                  file = "data/exclusive_mk_models.rds")

# # stochastic character maps
# nsim <- 1000
# simmap_ard <- make.simmap(tree, inxn, model = "ARD", nsim = nsim)
# simmap_er <- make.simmap(tree, inxn, model = "ER", nsim = nsim)
# simmap_sym <- make.simmap(tree, inxn, model = "SYM", nsim = nsim)
# # saveRDS(simmap_ard, file = "data/exclusive_simmap_ard.RDS")
# # saveRDS(simmap_er, file = "data/exclusive_simmap_er.RDS")
# # saveRDS(simmap_sym, file = "data/exclusive_simmap_sym.RDS")

# # weighted simmaps
# exclusive_mk_models <- readRDS("data/exclusive_mk_models.rds")
# aic_weights <- aic.w(unlist(lapply(mk_models, function(i) AIC(i))))
# nsum <- 1000
# sample_freq <- round(aic_weights * nsum)
# simmap_files <- c("data/exclusive_simmap_ard.RDS",
#                   "data/exclusive_simmap_er.RDS",
#                   "data/exclusive_simmap_sym.RDS")
# for (i in 1:length(simmap_files)) { 
#   simmap_obj <- readRDS(simmap_files[i])
#   if (i == 1) {
#     simmap_subsample <- simmap_obj[sample(1:length(simmap_obj),
#                                           sample_freq[i])]
#   } else {
#     simmap_subsample <- append(simmap_subsample,
#                                simmap_obj[sample(1:length(simmap_obj),
#                                                  sample_freq[i])])
#   }
#   rm(simmap_obj)
# }

# class(simmap_subsample) <- c("multiSimmap", "multiPhylo")
# saveRDS(simmap_subsample, file = "exclusive_simmap_weighted_avg.RDS")

# simmap_subsample <- readRDS("data/simmap_weighted_avg.RDS")
# sum_dmap <- density.multiSimmap(simmap_subsample)
# saveRDS(sum_dmap, "data/exclusive_sum_dmap.RDS)

# # q matrix
# exclusive_mk_models <- readRDS("data/exclusive_mk_models.rds")
# qmatrices <- lapply(exclusive_mk_models,
#                     as.Qmatrix)
# qmatrix_array <- abind::abind(qmatrices[[1]],
#                               qmatrices[[2]],
#                               qmatrices[[3]], 
#                               along = 3)
# qmatrix_avg <- qmatrices[[1]]
# for (i in 1:nrow(qmatrix_array)) {
#   for (j in 1:ncol(qmatrix_array)) {
#     qmatrix_avg[i, j] <- weighted.mean(qmatrix_array[i, j, ], aic_weights)
#   }
# }
# qmatrix_na <- qmatrix_avg
# for (i in 1:nrow(qmatrix_na)) {
#   qmatrix_na[i, i] <- NA
# }
# saveRDS(qmatrix_na, file = "data/exclusive_qmatrix.rds")

# Sensitivity results ####

# Pagel:

# Read in tests:
pagel_original <- readRDS("data/pagel_original.rds")
pagel_inclusive <- readRDS("data/pagel_inclusive.rds")
pagel_exclusive <- readRDS("data/pagel_exclusive.rds")

# calculate AIC weights:
pagel_aic <- cbind(aic.w(pagel_original$pagel_aic) * 100,
                   aic.w(pagel_inclusive$pagel_aic) * 100,
                   aic.w(pagel_exclusive$pagel_aic) * 100)
colnames(pagel_aic) <- c("original", "inclusive", "exclusive")
rownames(pagel_aic) <- c("independent", "both_dependent", "MB_dep_WN", "WN_dep_MB")

# calculate AIC differences between models:
pagel_delta_aic <- cbind(pagel_original$pagel_aic - min(pagel_original$pagel_aic),
                         pagel_inclusive$pagel_aic - min(pagel_inclusive$pagel_aic),
                         pagel_exclusive$pagel_aic - min(pagel_exclusive$pagel_aic))
rownames(pagel_delta_aic) <- rownames(pagel_aic)
colnames(pagel_delta_aic) <- colnames(pagel_aic)

# the model with the most support will have:
# 1. highest AIC weights; 2. delta AIC of 0
print(round(pagel_aic, 2)); print(round(pagel_delta_aic, 2))





# Mk models (AIC):

original_mk_models <- readRDS("data/original_mk_models.rds")
inclusive_mk_models <- readRDS("data/inclusive_mk_models.rds")
exclusive_mk_models <- readRDS("data/exclusive_mk_models.rds")

# calculate delta AIC for each model for each case
mk_delta_aic <- rbind(unlist(lapply(original_mk_models, function(i) AIC(i))),
                      unlist(lapply(inclusive_mk_models, function(i) AIC(i))),
                      unlist(lapply(exclusive_mk_models, function(i) AIC(i))))
rownames(mk_delta_aic) <- c("original", "inclusive", "exclusive")
mk_delta_aic <- apply(mk_delta_aic, 1, function(i) i - min(i))
# symmetrical rates is best supported for all cases again, although the 
# difference between SYM and ER is < 2 for the inclusive case




# Averaged simmaps: 
original_dmap <- readRDS("data/original_sum_dmap.rds")
inclusive_dmap <- readRDS("data/inclusive_sum_dmap.rds")
exclusive_dmap <- readRDS("data/exclusive_sum_dmap.rds")

# median number of transitions:
cbind(original_dmap$meds,
      inclusive_dmap$meds,
      exclusive_dmap$meds)

# mean number of transitions:
round(cbind(original_dmap$means,
            inclusive_dmap$means,
            exclusive_dmap$means), digits = 1)




# Transition rate (Q) matrices:
weighted_qmatrix <- function(mk_models_list, aic_weights) {
  qmatrices <- lapply(mk_models_list,
                      as.Qmatrix)
  qmatrix_array <- array(as.numeric(unlist(qmatrices)), 
                         dim = c(nrow(qmatrices[[1]]),
                                 ncol(qmatrices[[1]]),
                                 3))
  
  qmatrix_avg <- qmatrices[[1]]
  for (i in 1:nrow(qmatrix_array)) {
    for (j in 1:ncol(qmatrix_array)) {
      qmatrix_avg[i, j] <- weighted.mean(qmatrix_array[i, j, ], aic_weights)
    }
  }
  
  return(qmatrix_avg)
}
qmatrix_original <- weighted_qmatrix(original_mk_models, 
                        unlist(lapply(original_mk_models, function(i) AIC(i))))
qmatrix_inclusive <- weighted_qmatrix(inclusive_mk_models,
                        unlist(lapply(inclusive_mk_models, function(i) AIC(i))))
qmatrix_exclusive <- weighted_qmatrix(exclusive_mk_models,
                        unlist(lapply(exclusive_mk_models, function(i) AIC(i))))
qmatrix_list <- list(round(qmatrix_original, digits = 5),
                     round(qmatrix_inclusive, digits = 5),
                     round(qmatrix_exclusive, digits = 5))
for (i in 1:length(qmatrix_list)) {
  qmatrix_list[[i]] <- qmatrix_list[[i]][c(3, 2, 4, 1), c(3, 2, 4, 1)]
  for (j in 1:nrow(qmatrix_list[[i]])) {
    qmatrix_list[[i]][j, j] <- NA
  }
}
names(qmatrix_list) <- c("original", "inclusive", "exclusive")
lapply(qmatrix_list, function(i) i * 10^3)
