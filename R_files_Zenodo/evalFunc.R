# Eval function
evalFunc <- function(x) {
  # MPA_amount <- sum(x==0)
  # if (MPA_amount>(fraction_reserves*(length(x)-2*extra_patches)) | any(x[c(1:extra_patches,((length(x)-extra_patches):length(x)))] == 0)){ # Stops edges from being chosen
  #   return(0)
  # }else{
    MPA_metric = vector('numeric',length = num_trials)

    #set.seed(123)
    for (trial in 1:length(MPA_metric)){
      # source('../scripts/generate_disturbances.R')
      MPA_metric[trial] <- simulate_model(r=r_vector,K=K,omega=omega_value, connectivity_matrix=connectivity_matrix , fishing_vector = x*fishing_vector, length(x), max_time=max_time, disturbance_regime=disturbance_regime_trial[[trial]],objective = objective)
    }
    return(-(sum(MPA_metric))) 
 # }
}
