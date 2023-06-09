# analyses for trait evolution

# bookkeeping ####
# libraries
library(phytools)

# Read in the tree (trimmed phylogeny from McGee et al. 2020)
tree <- read.tree("data/McGee2020_tree.tre")
tree <- rotateNodes(tree, c(1078, 1061, 939)) # for visual clarity

# Read in the tip states:
data <- read.csv("data/McGee2020_tips.csv")

# Make named variables
reproduction <- setNames(data$reproduction, data$tip.label)
feeding <- setNames(data$feeding, data$tip.label)
inxn <- setNames(paste(reproduction, feeding, sep = " & "), data$tip.label)

# shorten the interaction names
{inxn <- gsub("substrate brooding & other", "Neither behavior", inxn)
  inxn <- gsub("substrate brooding & winnowing", "Winnowing only", inxn)
  inxn <- gsub("mouthbrooding & other", "Mouthbrooding only", inxn)
  inxn <- gsub("mouthbrooding & winnowing", "Both behaviors", inxn)}


# part 1: pagel 1994 correlation test ####

# this is a model where the evolution of each trait depends on the other:
fit_dep <- fitPagel(tree, x =  reproduction, y = feeding)

# this is a model where the evolution of mouthbrooding depends on winnowing
fit_ind <- fitPagel(tree, x = reproduction, y = feeding, dep.var = "x")

# winnowing depends on mouthbrooding, but not the other way around
fit_w_dep_m <- fitPagel(tree, x = reproduction, y = feeding, dep.var = "y")

# AIC support:
pagel_aic <- c(fit_dep$independent.AIC, # neither trait affects the other
               fit_dep$dependent.AIC, # both traits depend on each other
               fit_ind$dependent.AIC, # mouthbrooding depends on winnowing, but not the other way around
               fit_w_dep_m$dependent.AIC) # winnowing depends on mouthbrooding, but not the other way around
aic.w(pagel_aic)

# alternatively, load previously run test results:
pagel_test <- readRDS("data/pagel_original.rds")
aic.w(pagel_test$pagel_aic) 

# part 2: Mk model weights ####
mk_sym <- fitMk(tree, inxn)
mk_ard <- fitMk(tree, inxn, model = "ARD")
mk_er <- fitMk(tree, inxn, model = "ER")

# to save these:
# saveRDS(setNames(list(mk_ard, mk_er, mk_sym),
#                  c("ARD", "ER", "SYM")), 
#                  file = "data/original_mk_models.rds")

# to load previously run tests instead:
mk_models <- readRDS("data/original_mk_models.rds")
aic.w(unlist(lapply(mk_models, function(i) AIC(i))))
# greatest support for symmetrical rates

# part 3: stochastic character maps ####
# the original code was run with nsim = 10,000 maps and takes several hours to
# run;  I've changed it to 1000 here to avoid doing that to anyone's computer
# unintentionally (ideally you should run this on a cluster)
nsim <- 1000

# symmetrical rates model
simmap_sym <- make.simmap(tree, inxn,
                          model = "SYM",
                          nsim = nsim)
# saveRDS(simmap_sym,
#         file = "data/simmap_sym.RDS")

# equal rates model
simmap_er <- make.simmap(tree, inxn,
                         model = "ER",
                         nsim = nsim)
# saveRDS(simmap_er,
#         file = "data/simmap_er.RDS")

# all-rates-different model
simmap_ard <- make.simmap(tree, inxn,
                          model = "ARD",
                          nsim = nsim)
# saveRDS(simmap_ard,
#         file = "data/simmap_ard.RDS")


# part 4: weighted simmaps and q matrices ####
aic_weights <- aic.w(unlist(lapply(mk_models, function(i) AIC(i))))

# so first, we sample a number of stochastic character maps
# from each set of simmaps proportional to the AIC weight 
# of that model

# again, nsum was originally 10,000 (for the paper);
# changed it to 1000 here
nsum <- 1000
sample_freq <- round(aic_weights * nsum)

# read in simmap files and sample them according to their AIC weights
simmap_files <- c("data/simmap_ard.RDS",
                  "data/simmap_er.RDS",
                  "data/simmap_sym.RDS")

for (i in 1:length(simmap_files)) { 
  simmap_obj <- readRDS(simmap_files[i])
  if (i == 1) {
    simmap_subsample <- simmap_obj[sample(1:length(simmap_obj),
                                          sample_freq[i])]
  } else {
    simmap_subsample <- append(simmap_subsample,
                               simmap_obj[sample(1:length(simmap_obj),
                                                 sample_freq[i])])
  }
  rm(simmap_obj)
}

# classify & save
class(simmap_subsample) <- c("multiSimmap", "multiPhylo")
# saveRDS(simmap_subsample, file = "simmap_weighted_avg.RDS")

# summaries ####

# to generate the density map:
simmap_subsample <- readRDS("data/simmap_weighted_avg.RDS")
sum_dmap <- density.multiSimmap(simmap_subsample)


# alternatively:
sum_dmap <- readRDS("data/original_sum_dmap.rds")
data.frame(means = sum_dmap$means,
           medians = sum_dmap$meds,
           mins = sum_dmap$mins,
           maxs = sum_dmap$maxs)

# q matrices
qmatrices <- lapply(list(mk_ard, mk_er, mk_sym),
                    as.Qmatrix)
qmatrix_array <- abind::abind(qmatrices[[1]],
                              qmatrices[[2]],
                              qmatrices[[3]], 
                              along = 3)
qmatrix_avg <- qmatrices[[1]]
for (i in 1:nrow(qmatrix_array)) {
  for (j in 1:ncol(qmatrix_array)) {
    qmatrix_avg[i, j] <- weighted.mean(qmatrix_array[i, j, ], aic_weights)
  }
}
qmatrix_na <- qmatrix_avg
for (i in 1:nrow(qmatrix_na)) {
  qmatrix_na[i, i] <- NA
}

print(round(qmatrix_na * 10^3, digits = 2))
