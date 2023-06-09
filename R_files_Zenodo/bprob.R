# bprob.r: function that uses breeding probability to split population into breeders and nonbreeders at start of breeding season
# input: population size, breeding probability (1 - temporary emigration; estimated from multistate mark-recapture model)
# processing: determine whether each adult breeds or doesn't breed
# output: return vector with 1/0 for breeding/nonbreeding
bprob <- function(pop, prob){
 out <- rbinom(pop, 1, prob)
 return(out)
# could modify function to deal with Markovian temporary emigration
}


# not run:
# adults <- 100
# bprob(adults,.5) -> breeders

