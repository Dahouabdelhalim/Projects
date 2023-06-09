
## MVN random deviates function
rmvn <- function(n,mu,sig) {
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

### Function to compute the credible interval  for the photosynthesis endpoint
compute_CI_photo <- function(model, dat, nsims = 1000){
  
  Xp <- predict(model,dat,type="lpmatrix") # Get linear predictor matrix minus the offset
  
  which_re <- grepl("block|sample", names(coef(model)))
  mean_block_effect <- coef(model)[grepl("block", names(coef(model)))] %>% mean()
  #mean_sample_id_effect <- coef(model)[grepl("sample", names(coef(model)))] %>% mean()
  
  # 1000 replicate param. vectors
  br <- rmvn(nsims,coef(model),model$Vp) 
  
  # Remove random effects columns to make condition on average random effect levels
  Xp <- Xp[,!which_re]
  br <- br[, !which_re]
  
  
  # Intialise the matrix for storage of celld estimates
  celldmult <- matrix(NA, nrow = nrow(dat), ncol = nsims)
  
  for (i in 1:nsims){ 
    pr <- Xp %*% br[i,] + mean_block_effect #+ mean_sample_id_effect ## replicate predictions
    celldmult[,i] <- pr 
  }
  
  # Convert predictions to the natural scale (added offset here) and bind to the data
  df_bayes <- cbind(dat, celldmult + dat$t0) 
  
  return(df_bayes)
  
}


### Function to compute the credible interval data for the growth endpoint
compute_CI_growth <- function(model, dat, nsims = 1000){
  
  Xp <- predict(model,dat,type="lpmatrix") # Get linear predictor matrix minus the offset
  
  which_re <- grepl("block|sample", names(coef(model)))
  mean_block_effect <- coef(model)[grepl("block", names(coef(model)))] %>% mean()
  #mean_sample_id_effect <- coef(model)[grepl("sample", names(coef(model)))] %>% mean()
  
  # 1000 replicate param. vectors
  br <- rmvn(nsims,coef(model),model$Vp) 
  
  # Remove random effects columns to make condition on average random effect levels
  Xp <- Xp[,!which_re]
  br <- br[, !which_re]
  
  
  # Intialise the matrix for storage of celld estimates
  celldmult <- matrix(NA, nrow = nrow(dat), ncol = nsims)
  
  for (i in 1:nsims){ 
    pr <- Xp %*% br[i,] + mean_block_effect #+ mean_sample_id_effect ## replicate predictions
    celldmult[,i] <- pr 
  }
  
  # Convert predictions to the natural scale (added offset here) and bind to the data
  df_bayes <- cbind(dat, exp(celldmult + log(dat$t0))) 
  
  return(df_bayes)
  
}
