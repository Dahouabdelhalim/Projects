################################################################################ 

# SCRIPT 1: SIMULATION OF CHRONOGRAMS/TRANSFORMATION INTO PHYLOGRAMS

# AUTHOR: Jeremy D. Wilson

# ARTICLE: Wilson et al. 2022 - Chronogram or phylogram for ancestral state
# estimation? Model-fit statistics indicate the branch lengths underlying a
# binary characterâ€™s evolution


################################################################################

# LOAD PACKAGES-----------------------------------------------------------------


require(ape)
require(TreeSimGM)
require(geiger)
require(tidyverse)
require(FossilSim)



# ENTERING IN KEY DATA----------------------------------------------------------


num_trees <- 5 # This was set to 5000 for article.

taxon_number_range <- 10:100 # Set to 10:1000 for article.

data <- tibble(.rows = num_trees, replicate = 1:num_trees)


# SIMULATING ULTRAMETRIC TREES WITH AGE-DEPENDANT SPECIATION & EXTINCTION-------


sim.taxa.mod <- function(NumTrees, TaxonNoRange){
  sims <- replicate(NumTrees,
                    sim.taxa(numbsim = 1,
                             n = sample(TaxonNoRange, 1),
                             waitsp = "rweibull(0.4, 1)", #Age dep. speciation.
                             waitext = "rweibull(0.4, 2)", #Age dep. extinction.
                             symmetric = T,
                             complete = F)
  )
  names(sims) <- paste("tree", 1:NumTrees, sep = "_")
  sims <- purrr::map2(sims,
                      sample(1:100, num_trees, replace = T),
                      function(x, y){ rescale(x, "depth", y ) } # Rescaling to random depth.
  )
  sims <- purrr::map(sims, ladderize)
  return(sims) 
}

data$chronogram <- sim.taxa.mod(num_trees, taxon_number_range)
data$chronogram_depth <- purrr::map_dbl(data$chronogram, tree.max)


# PHYLOGRAM TRANSFORMATIONS-----------------------------------------------------

# Transformation 1.

autocorr.transform <- function(Phylo){
  zero_length_branches <- TRUE
  while(zero_length_branches == TRUE){
    for(i in Phylo$edge[,2]){
      xvalue <- 0
      while(xvalue <= 0){
        xvalue <- rnorm(1, mean = 1, sd = 0.2)
      }
      Phylo <- rescale(Phylo,
                       model = "lrate",
                       node = i,
                       rate = xvalue)
    }
    zero_length_branches <- any(Phylo$edge.length <= 0)
  }
  
  return(Phylo)
}

# Transformation 2.

uncorr.transform <- function(Phylo){
  zero_length_branches <- TRUE
  while(zero_length_branches == TRUE){
    for(i in 1:length(Phylo$edge.length)){
      xvalue <- 0
      while(xvalue <= 0){
        xvalue <- rnorm(1, mean = 1, sd = 0.4)
      }
      Phylo$edge.length[i] <- 
        Phylo$edge.length[i] * xvalue
    }
    zero_length_branches <- any(Phylo$edge.length <= 0)
  }
  return(Phylo)
}

data$phylogram_type <- sample(c("autocorrelated", "uncorrelated"),
                              num_trees,
                              replace = T)

data$phylogram <- purrr::map2(data$chronogram,
                              data$phylogram_type,
                              function(chrono, phylo_type){
                                if(phylo_type == "autocorrelated"){
                                  autocorr.transform(chrono)
                                } else {
                                  uncorr.transform(chrono)
                                }
                              }
)

data$phylogram_depth <- purrr::map_dbl(data$phylogram, tree.max)
data$tree_size <- purrr::map_dbl(data$chronogram,
                                 function(x) length(x$tip.label)
)


# SAVE TREE SETS----------------------------------------------------------------


save(data, file = "trees.RData")