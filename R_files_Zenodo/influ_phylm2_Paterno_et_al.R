#this script was written by the Senisphy author Gustavo Paterno, who kindly provided
#it in order to run these analyses for models with >1 predictor. He has given his 
#permission for the code to be provided - please ensure you cite him accordingly:
#Paterno, G.B., C. Penone, and G.D.A. Werner, sensiPhy: An r-package for sensitivity 
#analysis in phylogenetic comparative methods. Methods in Ecology and Evolution, 2018. 
#9(6): p. 1461-1467.

# Function to perform influ sensitivity analysis for multiple regression models
influ_phylm2 <- function(formula, data, phy, model = "lambda", cutoff = 2,
           track = TRUE, ...) {
    if (class(formula) != "formula")
      stop("formula must be class 'formula'")
    if (class(data) != "data.frame")
      stop("data must be class 'data.frame'")
    if (class(phy) != "phylo")
      stop("phy must be class 'phylo'")
    if ((model == "trend") && (ape::is.ultrametric(phy)))
      stop("Trend is unidentifiable for ultrametric trees., see ?phylolm for details")
    else
    
    # Check match between data and phy
    data_phy <- match_dataphy(formula, data, phy)
    #Calculates the full model, extracts model parameters
    full.data <- data_phy$data
    phy <- data_phy$phy
    N               <- nrow(full.data)
    mod.0           <- phylolm::phylolm(formula,
                                        data = full.data,
                                        model = model,
                                        phy = phy)
    coef.0 <- coef(mod.0)
    terms  <- names(coef.0)
    sp     <- rownames(full.data)
    
    # Fit models--------------------------------------------------------------------
    total <- N
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    
    mods <- lapply(seq_len(N), function(x) {
      setTxtProgressBar(pb, x)
      crop.data <- full.data[c(1:N)[-x], ]
      crop.phy  <-  ape::drop.tip(phy, phy$tip.label[x])
      phylolm::phylolm(formula,
                       data = crop.data,
                       model = model,
                       phy = crop.phy)
     
    })
    close(pb)
    
    modc <- do.call("rbind", lapply(mods, coef))
    modc <- lapply(seq_len(N),
                   function(x) {
                     t(summary(mods[[x]])$coefficients)[1, ]
                   })
    modc <- do.call("rbind", modc)
    
    modp <- lapply(seq_len(N),
                   function(x) {
                     t(summary(mods[[x]])$coefficients)[4, ]
                   })
    modp <- do.call("rbind", modp)
    
    # Calculate absolute difference in estimate
    sensi.estimates <- data.frame()
    
    for (i in terms) {
      estimate   <- modc[, i]
      DIF        <- modc[, i] - coef.0[i]
      sDIF       <- DIF / sd(DIF)
      pval       <- modp[, i]
      perc       <- round((abs(DIF / coef.0[i])) * 100, digits = 1)
      estim      <-
        data.frame(
          species       = sp,
          term          = i,
          estimate      = estimate,
          DIFestimate   = DIF,
          sDIFestimate   = sDIF,
          estimate.perc = perc,
          pval
        )
      sensi.estimates <- rbind(sensi.estimates, estim)
    }

    # Influential species-----------------------------------------------------------
    sensi.estimates$influential <- abs(sensi.estimates$sDIFestimate) > 2
    
    #Creates a list with full model estimates:
    param0 <- list(
      coef = phylolm::summary.phylolm(mod.0)$coefficients,
      aic = phylolm::summary.phylolm(mod.0)$aic,
      optpar = mod.0$optpar
    )
    
    
    #Generates output:
    res <- list(
      call = match.call(),
      cutoff = cutoff,
      formula = formula,
      terms   = terms,
      full.model.estimates = param0,
      sensi.estimates = sensi.estimates,
      data = full.data
    )
    return(res)
  }
