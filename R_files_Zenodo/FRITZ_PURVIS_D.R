## load packages
library(caper)
library(ape)

#The distributions dr and db can also be used to assign p-values to dobs, i.e., if dobs is larger than 95% of dr values then the 
#distribution of the trait is significantly more overdispersed than the random expectation, 
#if dobs is less than 95% of db values, the character is significantly more clumped than the Brownian expectation.

## read tree
tree <- read.tree ("C_tree.nwk")

## read data
dat <- read.csv ("Discrete_13.csv", header = TRUE)
HEDYCHIUM <- comparative.data(phy = tree, data = dat, names.col = Species, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

## calculating phylogenetic signals for characters 1 to 13
result1 = phylo.d(data=HEDYCHIUM, binvar = slender, permut = 1000)
result1

result2 = phylo.d(data=HEDYCHIUM, binvar = appear, permut = 1000)
result2

result3 = phylo.d(data=HEDYCHIUM, binvar = rachis, permut = 1000)
result3

result4 = phylo.d(data=HEDYCHIUM, binvar = bract, permut = 1000)
result4

result5 = phylo.d(data=HEDYCHIUM, binvar = density, permut = 1000)
result5

result6 = phylo.d(data=HEDYCHIUM, binvar = cincinnus, permut = 1000)
result6

result7 = phylo.d(data=HEDYCHIUM, binvar = flowers, permut = 1000)
result7

result8 = phylo.d(data=HEDYCHIUM, binvar = calyx, permut = 1000)
result8

result9 = phylo.d(data=HEDYCHIUM, binvar = bending, permut = 1000)
result9

result10 = phylo.d(data=HEDYCHIUM, binvar = claw, permut = 1000)
result10

result11 = phylo.d(data=HEDYCHIUM, binvar = labellum, permut = 1000)
result11

result12= phylo.d(data=HEDYCHIUM, binvar = blotch, permut = 1000)
result12

result13= phylo.d(data=HEDYCHIUM, binvar = stigma, permut = 1000)
result13