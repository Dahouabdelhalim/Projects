
library(rlang)

## MVN random deviates function
rmvn <- function(n,mu,sig) {
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

### Function to compute the credible interval data

compute_CI_data_exp <- function(model, dat, nsims = 1000){
  
  Xp <- predict(model,dat,type="lpmatrix") # Get linear predictor matrix minus the offset
  
  which_re <- grepl("block|sample", names(coef(model)))
  #mean_block_effect <- coef(model)[grepl("block", names(coef(model)))] %>% mean()
  #mean_sample_id_effect <- coef(model)[grepl("sample", names(coef(model)))] %>% mean()
  
  # 1000 replicate param. vectors
  br <- rmvn(nsims,coef(model),model$Vp) 
  
  # Remove random effects columns to make condition on average random effect levels
  Xp <- Xp[,!which_re]
  br <- br[, !which_re]
  
  
  # Intialise the matrix for storage of celld estimates
  celldmult <- matrix(NA, nrow = nrow(dat), ncol = nsims)
  
  for (i in 1:nsims){ 
    pr <- Xp %*% br[i,] #+ mean_block_effect #+ mean_sample_id_effect ## replicate predictions
    celldmult[,i] <- pr 
  }
  
  # Convert predictions to the natural scale (added offset here) and bind to the data
  df_bayes <- cbind(dat, exp(celldmult + log(dat$t0))) 
  
  return(df_bayes)
  
}


compute_CI_data_linear <- function(model, dat, nsims = 1000){
  
  Xp <- predict(model,dat,type="lpmatrix") # Get linear predictor matrix minus the offset
  
  which_re <- grepl("block|sample", names(coef(model)))
  #mean_block_effect <- coef(model)[grepl("block", names(coef(model)))] %>% mean()
  #mean_sample_id_effect <- coef(model)[grepl("sample", names(coef(model)))] %>% mean()
  
  # 1000 replicate param. vectors
  br <- rmvn(nsims,coef(model),model$Vp) 
  
  # Remove random effects columns to make condition on average random effect levels
  Xp <- Xp[,!which_re]
  br <- br[, !which_re]
  
  
  # Intialise the matrix for storage of celld estimates
  celldmult <- matrix(NA, nrow = nrow(dat), ncol = nsims)
  
  for (i in 1:nsims){ 
    pr <- Xp %*% br[i,] #+ mean_block_effect #+ mean_sample_id_effect ## replicate predictions
    celldmult[,i] <- pr 
  }
  
  # Convert predictions to the natural scale (added offset here) and bind to the data
  df_bayes <- cbind(dat, celldmult + dat$t0) 
  
  return(df_bayes)
  
}



compute_IR_and_CIs <- function(dat, model, response = c("growth", "photosynthesis")[[1]], variable = c("Diuron", "DIN")[[1]]){

  # ================================================
  # Computes the IR and CIs for a model, given some prediction data, the endpoint and the stressor (Diuron or DIN) 
  # ================================================
  
  # Select the model we need based on the endpoint
  compute_CI_data <- switch(response, "growth" = compute_CI_data_exp, "photosynthesis" = compute_CI_data_linear)
  
  # What is the variable used in the experiment
  variable_e <- parse_expr(paste0(variable, "_num"))
  
  # Substitute the variable of interest's name into the code
  body_expr <- expr({
    
    # Make predictions based
    all_bayes <- compute_CI_data(model = model, dat = dat, nsims = 1000) %>% within({
      Treatment <- case_when(Light_num == 80 & !!variable_e == 0 ~ "Control",
                             Light_num == 80  ~ paste0(variable,"_only"),
                             !!variable_e == 0  ~ "Light_only",
                             TRUE ~ "Light_Stress")
    }) 
    
    # # Generate the control matrix
    control_bayes <- all_bayes %>% filter(Treatment == "Control")
    control_mat <- control_bayes[ match(all_bayes$hours, control_bayes$hours) ,]
    
    # Generate light only matrix
    light_bayes <- all_bayes %>% filter(!!variable_e == 0)
    light_mat <- light_bayes[ 
      match(paste(all_bayes$hours, all_bayes$Light_num), 
            paste(light_bayes$hours, light_bayes$Light_num)) ,]
    
    # Generate Stress only matrix (Diuron or DIN)
    stress_bayes <- all_bayes %>% filter(Light_num == 80 )
    stress_mat <- stress_bayes[ 
      match(paste(all_bayes$hours, with(all_bayes, !!variable_e)), 
            paste(stress_bayes$hours, with(stress_bayes, !!variable_e))) ,]
    
    
    # Compute the IR matrix
    IR_mat2 <- log((select(all_bayes, `1`:`1000`)/select(control_mat, `1`:`1000`))/
                     ((select(light_mat, `1`:`1000`)/select(control_mat, `1`:`1000`))*
                        (select(stress_mat, `1`:`1000`)/select(control_mat, `1`:`1000`))))
    
    # Compute the credible intervals
    CItrend <- cbind(data.frame(t(apply(IR_mat2, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))), dat) %>% within({
      Treatment <- case_when(Light_num == 80 & !!variable_e == 0 ~ "Control",
                             Light_num == 80  ~ paste0(variable,"_only"),
                             !!variable_e == 0  ~ "Light_only",
                             TRUE ~ paste0("Light = ", Light_num, ", ", variable," = ", !!variable_e))}) %>% 
      filter(!(Treatment %in% c("Light_only", paste0(variable,"_only") ,  "Control")))
    
  })
  
  eval(body_expr)
  
  return(CItrend)
  
}
