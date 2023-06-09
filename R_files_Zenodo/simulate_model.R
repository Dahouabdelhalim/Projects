


# This runs a single simulation of the two-patch-model

simulate_model = function(r1,K1,delta,epsilon,omega,prob_disturb, severity_disturb, gamma, distance, output_time_series=FALSE){
  
  #set.seed(123)
  num_years=100
  N1 = 0.6 + vector('numeric',num_years)
  N2 = 0.6 + vector('numeric',num_years)
  Yield =  vector('numeric',num_years)
  for (t in 1:num_years){
    
    #r1s = rnorm(1,mean=r1,sd=0.3)
    # Beverton-Holt model
    N1_temp = rnorm(1,mean=r1,sd=0.25)*(N1[t]^omega)/(1 + (N1[t]^omega)/K1)
    N2_temp = rnorm(1,mean=r1,sd=0.25)*(N2[t]^omega)/(1 +(N2[t]^omega)/K1)
    ## Dispersal
    N1[t+1] =  N1_temp - epsilon*N1_temp + exp(-delta*distance)*epsilon*N2_temp
    N2[t+1] =  N2_temp - epsilon*N2_temp + exp(-delta*distance)*epsilon*N1_temp
    
    ## Fisheries yield calculation
    Yield[t+1] = (1 - exp(-delta*distance))*(epsilon*N2_temp) + (1 - exp(-delta*distance))*(epsilon*N1_temp)
    
    
    ########Alternative with dispersal first ######
    # N1_temp =  N1[t] - epsilon*N1[t] + exp(-delta*distance)*epsilon*N2[t]
    # N2_temp =  N2[t] - epsilon*N2[t] + exp(-delta*distance)*epsilon*N1[t]
    # 
    # N1[t+1] = r1*(N1_temp^omega)/(1 + (N1_temp^omega)/K1)
    # N2[t+1] = (r1-1)*(N2_temp^omega)/(1 + (N2_temp^omega)/K1)
    # 
    ##########
    # Disturbance
    if (rbinom(n = 1,size=1,prob = prob_disturb)){
      N1[t+1] = severity_disturb*N1[t+1]
      N2[t+1] = ifelse(rbinom(1,1,prob=exp(-gamma*distance)),severity_disturb,1)*N2[t+1]
    }
    if (rbinom(n = 1,size=1,prob = prob_disturb)){
      N1[t+1] = ifelse(rbinom(1,1,prob=exp(-gamma*distance)),severity_disturb,1)*N1[t+1]
      N2[t+1] = severity_disturb*N2[t+1]
    }

  }
  

  #persist = ifelse((N1[num_years]+N2[num_years])>0.01,1,0)
# Return either the time series or the measures of persistence   
  if (output_time_series == FALSE){
    return(c(mean(N1[11:num_years]+N2[11:num_years]),ifelse((N1[num_years]+N2[num_years])>0.01,1,0),var(N1[11:num_years]+N2[11:num_years]),sum(Yield[11:num_years])))
  }else{
    return(rbind(N1,N2))
  }
  
   
}
