#this script was written by the Senisphy author Gustavo Paterno, who kindly provided
#it in order to run these analyses for models with >1 predictor. He has given his 
#permission for the code to be provided - please ensure you cite him accordingly:
#Paterno, G.B., C. Penone, and G.D.A. Werner, sensiPhy: An r-package for sensitivity 
#analysis in phylogenetic comparative methods. Methods in Ecology and Evolution, 2018. 
#9(6): p. 1461-1467.

# summary function for influ_phylolm2
summary_influ2 <- function(x){
  estim <- x$sensi.estimates
  estim_influ <- estim[abs(estim$sDIFestimate) > x$cutoff, ]
  ord <- order(estim_influ$term, abs(estim_influ$sDIFestimate), decreasing = TRUE)
  estim_influ <- estim_influ[ord, ]
  estim_influ$cutoff <- paste("|sDIF| >", x$cutoff)
  
  colnames(estim_influ) <- c(
    "removed species", "term", "estimate", "DIFestimate", "sDIFestimate", 
    "change (%)", "Pvalue", "influential", "cutoff"
  )
  estim_out <- split(estim_influ, estim_influ$term)
  sp_out <- split(estim_influ[c("removed species")], estim_influ$term)
  sp_out <- lapply(sp_out, setNames, c("influential species"))
  
  res <- list(estimates = estim_out,
              influential = sp_out)
  return(res)
}
