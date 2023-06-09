#This script creates functions to simulate resource dynamics
#in a population of moving animals
#Claire Teitelbaum
#Date created: 29 April 2020

#uniform distribution, random at each point in time
#bounded by 0 and 1
sim_A_unif = function(tmax,min=0,max=1){
  A = runif(tmax , min = min , max = max)
  return(A)
}

#normal distribution, random at each point in time
sim_A_norm = function(tmax,mean=0.5,sd=0.5){
  A = rnorm(tmax , mean = mean , sd = sd)
  A[A<0] = 0
  A[A>1] = 1
  return(A)
}

#random walk
#range is the amount that A can change in a single timestep
sim_A_rw = function(tmax,range=0.1){
  A = vector(mod = "numeric" , length = tmax)
  A[1] = runif(1 , 0 , 1)
  for(i in 2:length(A)){
    #if A is very close to 1 or 0, implement reflecting boundary
    if(A[i-1] < range){
      min = -1*A[i]
      max = range + A[i]
    } else if(A[i-1] > (1-range)){
      max = 1 - A[i-1]
      min = -1*range - (1 - A[i-1])
    } else{
      min = -1*range
      max = range
    }
    A[i] = A[i-1] + runif(1,min = min , max = max)
    if(A[i]<0) A[i] = 0
    if(A[i]>1) A[i] = 1
  }
  return(A)
}

#functions for linear interpolation of resource availability 
#to create some autocorrelation at natural sites
#b is the time frame at which resources fluctuate randomly
#linear interpolation between these points
sim_A_linear_interp = function(tmax,b){
  vals = runif(ceiling(tmax/b)+1 , 0 , 1)
  A = c(sapply(1:(length(vals)-1) , function(x){
    seq(vals[x],vals[x+1],length.out = pars$b)
  }))
  return(A)
}