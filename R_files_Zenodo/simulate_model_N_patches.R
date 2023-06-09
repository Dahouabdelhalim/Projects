


# N-patch model


simulate_model = function(r,K,omega, connectivity_matrix, fishing_vector,  total_patches, max_time, disturbance_regime,objective){
  #N = matrix(0.60,nrow = total_patches,ncol=max_time)
  N = matrix(0,nrow = total_patches,ncol=max_time)
  N[which(fishing_vector==0),]=0.1 # Set up initial patches
  
  total_fished = matrix(0,nrow = total_patches,ncol=max_time)
  #disturbance = rbinom(n = 100,size = num_patches,prob = 0.01)
  #r = r
  for (t in 1:(max_time-1)){
    
    # Beverton-Holt model rnorm(num_patches,r,0.2)
    N_temp = rnorm(total_patches,r,0.5)*(N[1:total_patches,t]^omega)/(1 + (N[1:total_patches,t]^omega)/K) #+ rnorm(num_patches,0,0.01)
    
    # Dispersal
    N[1:total_patches,t+1] =   1*connectivity_matrix%*%N_temp#exp(-d*distance)*dp*N2_temp
    
    # Remove biomass from patches that are only there for boundary conditions (this also means dispersers sent outside the landscape cannot be caught be fishers)
    N[1:extra_patches,t+1]=0
    N[(total_patches-extra_patches+1):total_patches,t+1]=0
    
    # Fishing 
    total_fished[1:total_patches,t+1] = fishing_vector*N[1:total_patches,t+1] # Amount fished
    N[1:total_patches,t+1] = N[1:total_patches,t+1] - total_fished[1:total_patches,t+1]
    
    # Disturbance (this disturbance happens after recruitment, dispersal, and fishing)
    N[1:total_patches,t+1] = disturbance_regime[t,1:total_patches]*N[1:total_patches,t+1]
    
  }
  
  
  # Should I be summing up all of N here? 
   #matplot(t(N),type = 'l',ylim=c(0,1))
  # sum(colSums(N[,1:max_time])>num_patches*0.1)
  switch(objective,
  persistence     = {return(ifelse(sum(N[1:total_patches,max_time])>0.1,1,0))},
  population_size = {return(sum(colSums(N[,11:ncol(N)])))},
  fishing         = {return(sum(colSums(total_fished)[11:ncol(N)]))}
  
  )
}