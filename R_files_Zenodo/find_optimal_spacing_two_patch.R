

#Function to find optimal spacing in two-patch model
find_optimal_spacing = function(r1,K1,delta,epsilon,omega,prob_disturb, severity_disturb, gamma,objective,output_time_series=FALSE){
  
  # Objective choice
  #objective = 'persistence'
  #objective = 'yield'
  
  #parameters
  
  num_runs= 100
  
  #distance = 30
  distance_vector = seq(1,50,1)
  total_pop = vector('numeric',length(distance_vector))
  persistence = NA*vector('numeric',length(distance_vector))
  total_var = vector('numeric',length(distance_vector))
  total_yield = vector('numeric',length(distance_vector))
  best_run = NULL
  
  #prob_disturb = 0.15;severity_disturb=0.1
  #delta=0.03
  #source('../scripts/two_patch/simulate_model.R')
  for (i in 1:length(distance_vector)){
    #set.seed(123)
    distance = distance_vector[i]
    model_output= matrix(0,nrow=num_runs,ncol=4)
    for (trial in 1:nrow(model_output)){
      model_output[trial,] <- simulate_model(r1,K1,delta,epsilon,omega,prob_disturb, severity_disturb, gamma, distance)
      
      
      
    }
    
    
    
    #total_pop[distance] = mean(model_output[,1])
    #if (mean(model_output[,1])>0.01){ # The 0.01 is a threshold to call the population extinct
    persistence[i] = sum(model_output[,2])
    total_pop[i] = mean(model_output[,1])
    total_var[i] = mean(model_output[,3])
    total_yield[i] = mean(model_output[,4])
    
    ### option to break loop at speed up code
    best_run = c(best_run,distance_vector[which(persistence==max(persistence,na.rm=T))[1]])
    #print(best_run)
    if (i>10 & all(tail(best_run,10)==tail(best_run,1))){
      #print(length(best_run))
      break}
    
    ###
    #}else{
    #  persistence[distance] = NA
    #  total_pop[distance] = NA
  }
  
 
  #print(min(persistence))
  if (objective=='persistence'){
    return(c(distance_vector[which(persistence==max(persistence,na.rm=T))[1]],persistence[which(persistence==max(persistence,na.rm=T))[1]]/num_runs))
    
  }else{
  return(c(distance_vector[which(total_yield==max(total_yield,na.rm=T))[1]],total_yield[which(total_yield==max(total_yield,na.rm=T))[1]]/num_runs))
  }

  
  #return(c(which(total_pop==max(total_pop))[1], max(total_pop)))
  #my_loess <- loess(near_extinct~distance_vector)
  #return(c(distance_vector[which(predict(my_loess)==min(predict(my_loess)))[1]],max(total_pop[which(predict(my_loess)==min(predict(my_loess)))[1]])))
}

