# graph
# from div/clustering
# from lime.ps-graph

# version 20-4-16

# installs graph in PsiMinL and R-functions

# given: 
# years
# n.sou		= number of sources



# name of network
network <- paste(paste("lili", n.sou, min(years), max(years), "v2", sep="-"), "csv", sep=".")
network <- paste("data/rec.sou", network, sep="/")
print(network, quote=FALSE)

# nr.places = number of initialised places
nr.places <- 16
# nr of places in the PsiMinL container:
print(nr.places, quote=FALSE)

places <- 1:nr.places

# load network
# make container
p <- NULL
# Garbage Collection:
gc() 
# package for link-wise minimisation od Psi
require(PsiMinL)
# Loading required package: PsiMinL
# Loading required package: Rcpp
print(date(), quote=FALSE)
print("loads network", quote=FALSE)
p <- psiMinFileL(network, func = "psi*", numObjects = nr.places)
print(date(), quote=FALSE)
# created PsiMinL Container with nr.places PsiMinL objects

# functions
###############

# to initialise Tree of Change (ToC) for offspring we ned the symmetric difference to ram:
# (ToC is used in confined local search to keep track of including/excluding same links.)
symm.diff.ram <- function(L, M=linkIdsL(pmoL = ram)) {c(setdiff(L, M), setdiff(M, L))}

# to cross ewes with ram we need their union and intersection with ram:
union.ram <- function(L, M=linkIdsL(pmoL = ram)) {union(M, L)}
intersect.ram <- function(L, M=linkIdsL(pmoL = ram)) {intersect(M, L)}