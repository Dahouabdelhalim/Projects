#this script was written by the Senisphy author Gustavo Paterno, who kindly provided
#it in order to run these analyses for models with >1 predictor. He has given his 
#permission for the code to be provided - please ensure you cite him accordingly:
#Paterno, G.B., C. Penone, and G.D.A. Werner, sensiPhy: An r-package for sensitivity 
#analysis in phylogenetic comparative methods. Methods in Ecology and Evolution, 2018. 
#9(6): p. 1461-1467.

sensi_plot2 <- function(x){
  estim <- x$sensi.estimates
  coef0 <- data.frame(x$full.model.estimates$coef)
  coef0$term <- x$terms
  
  #Estimate plots
  g1 <-
    ggplot(estim, aes(x = estimate, fill = term, group = term, color = term)) +
    geom_histogram(show.legend = FALSE) +
    geom_vline(data = coef0, aes(xintercept = Estimate)) +
    facet_grid(~term, scales = "free") +
    theme_bw(base_size = 12) 
  #Pvalues plots
  
  g2 <-
    ggplot(estim, aes(x = pval, fill = term, group = term, color = term)) +
    geom_histogram(show.legend = FALSE, color = "white") +
    geom_vline(xintercept = 0.05) +
    facet_grid(~term, scales = "free") +
    theme_bw(base_size = 12) 
  return(list(estimates = g1,
              pvalue = g2))
}