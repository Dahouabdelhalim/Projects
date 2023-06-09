#this script was written by the Senisphy author Gustavo Paterno, who kindly provided
#it in order to run these analyses for models with >1 predictor. He has given his 
#permission for the code to be provided - please ensure you cite him accordingly:
#Paterno, G.B., C. Penone, and G.D.A. Werner, sensiPhy: An r-package for sensitivity 
#analysis in phylogenetic comparative methods. Methods in Ecology and Evolution, 2018. 
#9(6): p. 1461-1467.

# Usage example for influ_phylm2 
# Load sensiPhy
library(sensiPhy)


# Load new functions (replace out/emma_functions to your working directory)
#use this to get the working directory's pathway should you need it:
getwd()

source("emma_functions/influ_phylm2_Paterno_et_al.R")
source("emma_functions/plot_influ_phylm2_Paterno_et_al.R")
source("emma_functions/summary_influ_phylm2_Paterno_et_al.R")

# Load data
data("alien")

# Example 1 (only 1 predictor)
# run analysis:
influ1 <- influ_phylm2(log(gestaLen) ~ log(adultMass), phy = alien$phy[[1]], 
                     data = alien$data)

# Check full data model (all comparisons are in relation to this summary)
influ1$full.model.estimates$coef

# To check summary results:
# The output is a list with estimates for the influential species for each
# model predictor (in this case only two (see full model above))

# The estimates table include:
# Removed species: the species that was removed
# term: which model term is being evaluated (all from the full model summary)
# estimate: the estimate for the given term
# DIFestimate: the absolute difference in estimate in comparison to the full model
# sDIFestimate: the standardized difference in estimate
# Chage(%): the percentage of change in relation to the full model estimate
# Pvalue: the term pvaleu after removing the species (check if the significance of
# the terms changed after removing each species)
# influential: TRUE is sDIF > cutoff
# cutoff: the cutoff level defined in the function call
# This table is repeated for each model term (from the full model)
# Aditionally the summary also prints only the name of the species with 
# Standardize Difference in Estimate > than the cutoff

# Check estimates
summary_influ2(influ1)$estimates

# Most influential species (for each model term)
summary_influ2(influ1)$influential

# Visual diagnostics
# This plot shows the distribution of estimates for each model term after removing
# species one by one. Helps to the sensitivity of model estimates to individual species removal
# The vertical lines represent the full model estimate
sensi_plot2(influ1)[[1]]
# So The variation in the log(adulMass) estimate is quite low (0.140-0.155)

# This plot shows the distribution of pvalues for each model term after removing
# species one by one. Helps to check if the significance test from the full model
# summary is sensitivity to individual species removal
# Vertical lines represent (pvalue = 0.05) for reference.
sensi_plot2(influ1)[[2]]
# In all simulations, the pvalue for log(adultMass) remains fay below 0.05.

# Example 2 (only 2 predictors)
# run analysis:
influ2 <- influ_phylm2(log(gestaLen) ~ log(adultMass) + homeRange, phy = alien$phy[[1]], 
                       data = alien$data)

# Check full data model (all comparisons are in relation to this summary)
# Each term will have its own output on the sensitivity analysis
influ2$full.model.estimates$coef

# Check estimates
summary_influ2(influ2)$estimates

# Check estimates for a given term
summary_influ2(influ2)$estimates$`log(adultMass)`
summary_influ2(influ2)$estimates$homeRange

# Most influential species (for each model term)
summary_influ2(influ2)$influential

# Visual diagnostics
# Distribution of estimates 
sensi_plot2(influ2)[[1]]

# Pvalues
sensi_plot2(influ)[[2]]
# homeRange is never significant (pvalue > 0.05)
