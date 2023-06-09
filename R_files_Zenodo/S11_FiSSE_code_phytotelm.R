### Running FiSSE for Tonini et al. 2020 Evolution dataset for Evolution submission, 
#	"Improving inference and avoiding over-interpretation of hidden-state diversification 
#	models: specialized plant breeding has no effect on diversification in frogs"

# Load relevant libraries
	library(ape)
	library(phangorn)
	library(diversitree)
	library(phytools)

### Load data and functions
	tree <- read.tree("Tonini.tree.1579sp.tre")
	dat <- read.csv("Tonini.phyto.dat.1579sp.csv")
	traits <- dat[,2]
	names(traits) <- dat[,1]
	traits <- traits[tree$tip.label]
	# The tree is ultrametric, but fails the 'is.ultrametric' check (due precision of branch-length estimates). Force:
		tree <- force.ultrametric(tree, method = 'nnls') 

### FiSSE analysis
	set.seed(1852) # Set random-number seed to ensure reproducibility of simulation-based P-value
	res <- FISSE.binary(tree, traits, reps = 999)
	res
	#	$lambda0 (rate for non-phytotelm breeding)
	#	[1] 0.0701576

	#	$lambda1 (rate for phytotelm breeding)
	#	[1] 0.06830072 # Phytotelm rate is even lower under this approach

	#	$pval (proportion of simulations where observed test statistic (lambda1 - lambda0) is greater than the simulated value)
	#	[1] 0.487 # Two-tailed result is 0.974 (see below)

### Two-tailed pvalue is obtained as
	pval_2tailed   <- min(res$pval, 1 - res$pval) * 2 # P-value = 0.974