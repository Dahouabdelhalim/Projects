# adapt Ward seed
# version 20-5-1

# given 
# n.clusters = number of clusters
# seed.nr <=  n.cluster

# data
# lili
# n.clusters

# path
path <- paste("data/Ward", n.clusters, sep="-")

# link list
#load("data/lili-300-2011-2015.RObj")
load("data/rec.sou/lili.2011-2015-v2.300.RObj")

# the seed from Ward as set of sources
# or from intersection of results in v.Ls
if(WARD)
{
	n.clusters <- 2
	path <- paste("data/Ward", n.clusters, sep="-")
	load(paste("data/co-cit-all-top-300/seed-c", seed.nr, ".RObj", sep=""))
	print(length(seed))
	n.sou <- length(seed)
# all links of seed nodes
	L.0 <- lili$n.1 %in% seed | lili$n.2 %in% seed 
	m <- nrow(lili)/2
} 
if(INTER)
{
	load(file="data/Ward-15/v.Ls-v1.RObj")
	m <- nrow(v.Ls)
	L.0 <- v.Ls[,L1] & v.Ls[,L2]
	L.0 <- c(L.0, L.0)
	rm(v.Ls)
}
if(CPLC)
{
	print(TOWN, quot=FALSE)
	load(file=paste("data/CPLC/seven/L", TOWN, "v1.RObj", sep="-"))
	L.0 <- L
	rm(L)
}


# we have to adapt the seed with highest resolution
rel.R <- 1/20
rel.res <- TRUE
inclusion.first <- FALSE

source("scripts/Psi-local-L-search-man-seed-v4.R")

rm(lili)


n.sou <- 300
years <- 2011:2015
source("scripts/graph-v2.R")
# mutation variance and renewal variance (in percent):
mv <- 1; rv <- 6
R.levels <- c(20, 10, 5, 4, 3)


rm(L)
rm(L.0)


sizes.costs <- cbind(size = best.sizes[c(1, length(best.sizes))], cost = best.costs[c(1, length(best.costs))])
sc.file <- paste("sizes.costs", seed.nr, "RObj", sep=".")
save(sizes.costs, file=paste(path, sc.file, sep="/"))

