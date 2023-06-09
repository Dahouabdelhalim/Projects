surv <- function(pop, phi, days){
  # function to calculate fate of POP individuals with daily survival PHI after DAYS
  # input: pop, phi, days
  # processing
  pop <- rep(1,pop)
  #tmp <- rbinom(pop, 1, phi)
  pop
  for(i in 1:days){
    pop <- pop*(rbinom(pop,1,phi))
    #print(i)
    #print(pop)
  }
  return(pop)
}

