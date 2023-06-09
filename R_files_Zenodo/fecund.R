fecund <- function(pop, fec){
  # function that creates number of individuals produced by adults with Poisson distribution
  # could focus on overall number of metamorphs relative to number of adults, skipping larval dynamics
  return(rpois(pop, fec))
}
