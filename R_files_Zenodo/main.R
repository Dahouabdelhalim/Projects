# Main R script to run the analysis
rm(list = ls())
source("R/mechanistic_mods.R")
source("R/likelihoods.R")
source("R/fit_model.R")
source("R/helpers.R")

experiments <- c("behavior_volume2", 
               "density", 
               "host_number", 
               "karvonen1", 
               "karvonen2", 
               "paller1", 
               "paller2", 
               "time", 
               "volume")

# loop over each experiment, estimating parameters and saving output
for (i in seq_along(experiments)) {
  filename <- file.path("R", paste0(experiments[i], ".R"))
  source(filename)
  rm(list = c("d", "model_fits", "aic_tab"))
}
