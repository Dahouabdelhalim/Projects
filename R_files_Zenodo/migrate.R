migrate <- function(num, night_prob, rmort, num_guards, guard_limit){   #,fec, bmort){
#   input: number of crossing adults susceptible to road mortality, proportion of adults migrating each night, fecunedity, road mortality, number xing guards each night,
#       maximum number of animals saved per guard per night # drop breeding mortality
#   processing:
#       mortality of adults crossing road per night
#       #drop: fecundity of surviving adults
#       # drop: mortality of adults that survive inbound migration
#   output:
#       adults that make it into the pond to breed on each night
#       # drop: number of adults that will make outbound migration
  #print(flag) debugging parameter, set unique integer to flag different function calls
  if(length(num)>1) print("Check input - pop vector rather than pop total")

  # create object to store migrants each night, should fix to get whole numbers
  num * night_prob -> nxing
  round(nxing) -> nxing

  # calculate potential mortality
  pmort <- rep(NA, length(nxing))
  for(i in 1:length(nxing)){
     pmort[i] <- sum(rbinom(nxing[i], 1, rmort))
  }
  # number of salamanders saved from cars
  num_guards*guard_limit -> numsave
  # 'save' salamanders
  pmort - numsave -> mort
  # set mortality to zero if guards can save more salamanders than those that could be killed
  mort[mort<0] <- 0
  # subtract killed individuals from individuals that crossed road
  nxing - mort -> numcrossed
  return(numcrossed)

}


