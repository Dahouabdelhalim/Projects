




# Parameter values

num_patches=16
#specify_fraction_reserves = TRUE
#if (specify_fraction_reserves == FALSE){fraction_reserves = 0.05}

#num_reserves = 4#floor(fraction_reserves*num_patches)
extra_patches=10
total_patches=num_patches+2*extra_patches

max_time= 100
num_trials = 100


fishing_values_vec = 1#c(0.5,0.75,1)#seq(0.3,0.9,by=0.25)
r_value_vec = 3#2.2#seq(2.5,3,by=1)
omega_value_vec = 1.2#seq(1.2,1.4,by=0.2)
dispersal_value_vec = c(0.7)#seq(0.02,1,by=0.3)
K = 1

disturbance_prob_vec = seq(0,0.03,0.005)#seq(0,0.015,0.005)
disturbance_size_vec = 2#c(1,2)
disturbance_magnitude_vec = c(0.95)#seq(0.7,0.9,0.2)

objective_vec = c('persistence')#,'fishing','population_size')

# fishing_vector = rep(fishing_value,total_patches)
# r = rep(r_value,total_patches)
# omega=omega_value


params_list <- expand.grid(fishing_values_vec,r_value_vec,omega_value_vec,dispersal_value_vec,disturbance_prob_vec,disturbance_size_vec,disturbance_magnitude_vec,objective_vec)
optimal_reserves_matrix = matrix(0,nrow=nrow(params_list),ncol=num_reserves) 


names(params_list) <- c('fishing_values_vec','r_value_vec','omega_value_vec','dispersal_value_vec','disturbance_prob_vec','disturbance_size_vec','disturbance_magnitude_vec','objective_vec')

params_list$optimal_spacing = 0 
params_list$mean_objective_value = 0 


  
for (param_index in 1:nrow(params_list)){
  fishing_value = params_list$fishing_values_vec[param_index]
  fishing_vector = rep(fishing_value,total_patches)
  r_value = params_list$r_value_vec[param_index]
  r_vector = rep(r_value, total_patches)
  omega_value = params_list$omega_value_vec[param_index]
  
  dispersal_probability_param = params_list$dispersal_value_vec[param_index]
  source('../scripts/n_patch/set_up_linear_landscape.R')
  total_patches=nrow(connectivity_matrix)
  
  disturbance_prob_param =params_list$disturbance_prob_vec[param_index] 
  disturbance_size_param =  params_list$disturbance_size_vec[param_index]
  disturbance_magnitude =  params_list$disturbance_magnitude_vec[param_index]
  source('../scripts/n_patch/linear_disturbances_setup.R')
  
  objective = as.character(params_list$objective_vec[param_index])
  
 # Run code for one combination of parameters and return optimal configuration... 
  source('../scripts/n_patch/run_code_for_N_patch_model.R')
  params_list$optimal_spacing[param_index] <- mean(diff(optimal_reserves)-1)
  params_list$mean_objective_value[param_index] <- -ranked_configs[1,1+num_reserves]
  optimal_reserves_matrix[param_index,] <- optimal_reserves

 # print(paste(param_index,':',params_list$optimal_spacing[param_index],sep = ''))
  #print(head(ranked_configs))
}


#params_list <- cbind(params_list,optimal_reserves_matrix)




