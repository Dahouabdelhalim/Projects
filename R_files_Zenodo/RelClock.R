
rm( list = ls() ) 

# Create gene trees with relaxed clock.
# Step 2.


library( ape )
set.seed( 12345 )


# Number of gene trees
L <- 50

# Number of replicates
R <- 100

# s2 for the lognormal distribution of rates
s2 <- c( 0.01, 0.25 )

# Where to save gene trees
GTr.S1 <- as.list( numeric( R ) )
GTr.S2 <- as.list( numeric( R ) )


for( rep in 1:R ){

  ###################
  #                 #
  #     STEP  1     #
  #                 #
  ###################


  # Create L identical gene trees to the species tree
  gtree.s1 <- list()    # gene trees with  low rate variation among branches
  gtree.s2 <- list()    # gene trees with high rate variation among branches
  for( i in 1 : L ){
    gtree.s1[[ i ]] <- read.tree( "speciesPhylogram.tre" )
    gtree.s2[[ i ]] <- read.tree( "speciesPhylogram.tre" )
  }


  ###################
  #                 #
  #     STEP  2     #
  #                 #
  ###################


  # Number of branches
  nb <- length( gtree.s1[[ 1 ]]$edge.length )

  # Mean substitution rate over all loci and branches (in 100My)
  m0 <- 0.5

  # Mean substitution rate of each gene over branches
  mg.s1 <- rgamma( L, 10, 10 / m0 )
  mg.s2 <- rgamma( L, 10, 10 / m0 )

  # Simulate substitution rates in branches for all genes
  # assuming the independent rates model
  W.s1 <- W.s2 <- matrix( numeric( L * nb ), L, nb )
  for( i in 1:L ){
    W.s1[ i , ] <- rnorm( nb, log( mg.s1[ i ] ) - s2[1] / 2, sqrt( s2[1] ) )
    W.s2[ i , ] <- rnorm( nb, log( mg.s2[ i ] ) - s2[2] / 2, sqrt( s2[2] ) )
  }
  expW.s1 <- exp( W.s1 ); expW.s2 <- exp( W.s2 )


  # Multiply branch time duration with the branch rates
  for( i in 1:L ){
    gtree.s1[[ i ]]$edge.length <- gtree.s1[[ i ]]$edge.length * expW.s1[ i, ]
    gtree.s2[[ i ]]$edge.length <- gtree.s2[[ i ]]$edge.length * expW.s2[ i, ]
  }


  ###################
  #                 #
  #     STEP  3     #
  #                 #
  ###################

  # Multiply branch lengths with 3 to turn distances in 
  # subst/codon for use in codon models
  for( i in 1:L ){
    gtree.s1[[ i ]]$edge.length <- gtree.s1[[ i ]]$edge.length * 3
    gtree.s2[[ i ]]$edge.length <- gtree.s2[[ i ]]$edge.length * 3
  }

  
  # Save trees
  GTr.S1[[ rep ]] <- gtree.s1
  GTr.S2[[ rep ]] <- gtree.s2
  
  
  # Export trees
  if( rep == 1 ){ dir.create( "geneTrees_S1" ); dir.create( "geneTrees_S2" ) }
  
  lapply( gtree.s1, write.tree, paste( "geneTrees_S1//geneTrees_S1_rep", rep, ".tree", sep="" ), append = TRUE )
  lapply( gtree.s2, write.tree, paste( "geneTrees_S2//geneTrees_S2_rep", rep, ".tree", sep="" ), append = TRUE )

}
