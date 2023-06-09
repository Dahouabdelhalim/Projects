
rm( list = ls() )

# Create a species tree with branch time 
# durations to be used for simulation.
# Here the time unit is 100My.


library( ape )

tree <- read.tree( "speciesCladogram.tre" )

# Use the following ages for the nodes (time unit = 100My)
tAB <- 0.10
tAC <- 0.15
tAD <- 0.25
tAE <- 0.40
tAF <- 0.55
tAG <- 0.95
tHI <- 0.50
tAI <- 1.00

# Make all branches equal in time
tree$edge.length[ 1:16 ] = 1

# Define branch lengths
tree$edge.length[ 1 ] = tAI - tAG
tree$edge.length[ 2 ] = tAG - tAF
tree$edge.length[ 3 ] = tAF - tAE
tree$edge.length[ 4 ] = tAE - tAD
tree$edge.length[ 5 ] = tAD - tAC
tree$edge.length[ 6 ] = tAC - tAB
tree$edge.length[ c( 7, 8 ) ] = tAB
tree$edge.length[ 9 ] = tAC
tree$edge.length[ 10 ] = tAD
tree$edge.length[ 11 ] = tAE
tree$edge.length[ 12 ] = tAF
tree$edge.length[ 13 ] = tAG
tree$edge.length[ 14 ] = tAI - tHI
tree$edge.length[ c( 15, 16 ) ] = tHI


plot( tree )

tree$node.label <- NULL

write.tree( tree, "speciesPhylogram.tre" )

