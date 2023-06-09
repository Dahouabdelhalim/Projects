# parameters
# from div/clustering
# from lime.ps-parameters-v3

# version 18-8-19

# all combinations of parameters are set after giving a vector with three main parameters

print(parameters)

# m
m <- as.vector(parameters[1, "n.links"])

# year
#year <- as.vector(parameters[1, "year"])

# mutation variance
mut.var <- as.vector(parameters[1, "mut.var"])
m.var <- mut.var/100

# mutation variance used in renewal of flock
ren.var <- as.vector(parameters[1, "ren.var"])
r.var <- ren.var/100

# source of seeds
seed.source <- as.vector(parameters[1, "seed.source"])

# seed.nr
seed.nr <- as.vector(parameters[1, "seed"])

# phase
phase <- as.vector(parameters[1, "phase"])

# constants:

# maximal age of ram in dependence of ram size
# cf. selection:
k.in.0 <- m
k.in.1 <- 200000
mar.1 <- 40
k.in.2 <- 100000
mar.2 <- 80
mar.3 <- 100
max.age.ram <- mar.3


# population size:
print(pop.size <- 8, quote=FALSE)
nr.places <- 2*pop.size 

# nr. crossovers: 
nr.crossovers <- pop.size/2

# nr mutants:
mut.rate <- 1/2
nr.mutants <- floor(pop.size*mut.rate) 

# minimal relative range of valid communities:
#rel.R <- .1

# minimal relative size of main component
# which with we proceed
# (if main component is smaller we say the subgraph is split)
# (only after dynamic exclusion; confined exclusion does not need it
#  because we test whether main component is not to far from ram)
min.mc <- 1 - rel.R

# to avoid endless evolution
max.nr.gen <- 5000

# to avoid too short evolution
min.nr.gen <- 200

# to trigger renewal
min.mean.innov.rate <- 1/32

# minimal age of ram for renewal
ren.age.ram <- 10

# parameters for decision whether we have a substantial new ram
# whose genes should be stored
# ratio of new and old Psi less than
ram.impr <- .95
# ratio of distance of new from old ram to size of old larger than 
ram.diff <- .05


# protocols are stored in subdirectory protocols in directory of n.sou sources with rel resolution rel.R
#path.proto <- paste("data/PsiMinL", n.sou, min(years), max(years), "R", round(100*rel.R), sep="-")
path.proto <- paste(path, "protocols", sep="/")
# file names have this structure:
# paste("seed", seed.nr, "phase", phase, "mut.var", mut.var, "ren.var", ren.var, "exp", experiment, "protocol.txt", sep="-")
# we determine number of experiment:

files <- dir(path.proto)
files <- files[substr(files, 6, 9)==as.character(seed.nr)]

# phase
if (length(files) > 0)
{
	files <- files[substr(files, 17, 17)==phase]
}

# mut.var
if (length(files) > 0)
{
	dig.mv <- 1
	if(mut.var < 10)
	{
		files <- files[substr(files, 27, 27)==mut.var]
	} else {
		files <- files[substr(files, 27, 28)==mut.var]
		dig.mv <- 2
	}
	
}

# ren.var
if (length(files) > 0)
{
	dig.rv <- 1
	if(ren.var < 10)
	{
		files <- files[substr(files, 36+dig.mv, 36+dig.mv)==ren.var]
	} else {
		files <- files[substr(files, 36+dig.mv, 37+dig.mv)==ren.var]
		dig.rv <- 2
	}
}

# we set the experiment id number to last id + 1:
exp.nrs <- 0
experiment <- 1
if (length(files) > 0)
{
# exp.nrs can have more than one digit:
	if (length(files) > 9) 	
	{
		experiment <- max(as.numeric(substr(files[1:length(files)], 41+dig.rv+dig.mv, 42+dig.rv+dig.mv)), na.rm=TRUE) + 1
	} else {
		exp.nrs <- as.numeric(substr(files, 41+dig.rv+dig.mv, 41+dig.rv+dig.mv))
		experiment <- max(exp.nrs) + 1
	}
}
print(paste("experiment", experiment), quote=FALSE)	
#