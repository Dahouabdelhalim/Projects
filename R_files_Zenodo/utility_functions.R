# Functions used in analysis

# Plumage complexity
getcomplexity <- function(states, graph) {
	adj <- as_adj(graph) * (1-as.matrix(dist(states)))
	newgraph <- graph_from_adjacency_matrix(adj, mode="undirected")
	# filter out only 1 states
	newgraph <- delete_vertices(newgraph, which(states!=1))
	# get complexity score (number of contiguous plumage regions a given color)
	components(newgraph)$no
}
# getcomplexity(anc[1,], graph=g1)
# getcomplexity(anc[1,], graph=g2)


# 
# hpd <- function(x) {
# 	if (is.null(x)) {
# 		NULL
# 	} else {
# 		mu <- mean(x)
# 		ci <- HPDinterval(as.mcmc(x))
# 		c(mean=mu, lower=ci[1], upper=ci[2])	
# 	}
# }

# Compute 95% credible intervals
hpd <- function(x) {
	int <- HPDinterval(as.mcmc(x))
	c('y'=mean(x), 'ymin'=int[1], 'ymax'=int[2])
}

# Print out 95% credible intervals
hpdprint <- function(x, round=3) {
	int <- HPDinterval(as.mcmc(x))
	res <- c('y'=mean(x), 'ymin'=int[1], 'ymax'=int[2])
	res <- round(res, round)
	paste0(res[1], " [", res[2], ", ", res[3], "]")
	
}


# Calculate differences among parameters from MCMC output
# Note: ggplot has functions with '_' after it - e.g., select_() that allows for quoted variable names
makeletterplot <- function(df, param) {
	source("~/R/mult_ci_test.R")
	require(multcompView)
	xx <- df %>% group_by(colormech) %>% select_(param, "colormech") %>% mutate(id=1:n()) %>% spread_("colormech", param) %>% as.matrix()
	xx <- as.matrix(na.omit(xx))
	N <- nrow(xx)
	set.seed(1980)
	comb <- combn(colnames(xx), m=2)
	diffs <- xx[sample(1:N, N), comb[1, ]] - xx[sample(1:N, N), comb[2, ]]
	colnames(diffs) <- paste0(comb[1,], "-", comb[2,])
	res1 <- mult_ci_test(data.frame(xx[,-1]), hpd=TRUE, seq=TRUE)
	dif1 <- setNames(res1$signif, row.names(res1))
	let1 <- multcompLetters(dif1)$Letters  # get significance test letters
	f <- formula(paste0(param, "~colormech"))
	df2 <- aggregate(f, df, median)
	df2$signif <- let1[match(df2$colormech, names(let1))]
	ggplot(df, aes_string(x="colormech", y=param, fill="colormech")) +
		geom_violin(trim=F) +
		scale_fill_manual(values=pal) +
		theme_minimal() +
		geom_text(data=df2, aes(label=signif), col="white") +
		coord_flip() +
		theme(legend.position="none")
}


# Calculate changes to/from all plumage regions at each time step
rev_changes <- function(x) {
	# subset events, nodes (only want cases where changes occur)
	events <- x[[2]] # all events
	all.nodes <- x[[1]]
	keeps <- sapply(events, nrow) != 1 & all.nodes[,'node']!='225'
	events <- events[keeps]
	all.nodes <- x[[1]][keeps, ]
	# don't run, unless you're ready!! takes a long time ~5 minutes
	res <- pblapply(1:length(events), function(ii) {
		# setup start/end points
		branch.start <- as.numeric(strsplit(all.nodes[ii, 'start'], "")[[1]])
		branch.end <- as.numeric(strsplit(all.nodes[ii, 'end'], "")[[1]])
		switches <- events[[ii]]
		# generate matrix for all changes along branches
		branch.changes <- matrix(NA, nrow = nrow(switches) + 2, ncol = 22)  # matrix to hold results
		branch.changes[1, ] <- branch.start  # set starting state (first row of matrix)
		branch.changes[nrow(switches) + 2, ] <- branch.end  # set ending state (last row of matrix)
		for (i in 1:nrow(switches)) {
			branch.changes[i + 1, ] <- branch.changes[i, ]
			branch.changes[i + 1, switches[i, 'charnum'] + 1] <- switches[i, 'tostate']
		}
		# sanity checks
		all.equal(branch.changes[1, ], branch.start)
		all.equal(branch.changes[nrow(branch.changes), ], branch.end)
		# record "dispersal" events to/from different areas
		branch.trans <- array(NA, dim = c(ncol(branch.changes), ncol(branch.changes), nrow(branch.changes)-1))
		for (i in 1:(nrow(branch.changes) - 1)) {
			for (j in 1:ncol(branch.changes)) {
				for (k in 1:ncol(branch.changes)) {
					branch.trans[j, k, i] <- ifelse(branch.changes[i, j] == 1 & branch.changes[i, k] == 0 & branch.changes[i+1, k] == 1, 1, 0)
				}
			}
		}
		branch.trans
	})
	res
}

# Calculate harmonic mean from likelihood
harmmean <- function(x, burnin=0) {
	1/mean(1/(x[ceiling(burnin*length(x)):length(x)]))
}



# Color functions

# convert data.frame to vismodel object
df2vismodel <- function(x) {
	rownames(x) <- x[,1]
	x <- x[,-1]
	x <- x[,c('u','s','m','l')]
	x <- x[match(dat$illustrator_code, rownames(x)), ]  # reorder to match anatomical distance matrix
	class(x) <- class(teal.vm)
	attr(x, "visualsystem.achromatic") <- "none"
	attr(x, "qcatch") <- "Qi"
	attr(x, "relative") <- FALSE
	attr(x, "conenumb") <- 4
	attr(x, "data.visualsystem.chromatic") <- attr(teal.vm, "data.visualsystem.chromatic")
	x
}

# tree= phylogeny
# X = trait matrix
picmulti <- function(tree, X) {
	rownames(X) <- gsub(" ", "_", rownames(X))
	if (is.null(tree$node.label)) {
		tree <- makeNodeLabel(tree)
	}
	treep <- drop.tip(tree, which(!tree$tip %in% rownames(X)))	
	X <- X[treep$tip, ]
	pics <- sapply(1:ncol(X), function(x) {
		pic(x=X[,x], phy=treep)
	})
	rates <- apply(pics, 1, function(x) {mean(abs(x)^2)})
	subtree <- treep
	list(rates=rates, phy=subtree)
}

# function to compute complexity (patch = contiguous plumage regions with all JND < 1, complexity = N of patches)
getcomplexity_jnd <- function(jndmat, graph, plot=FALSE) {
	require(igraph)
	adj <- as_adj(graph) * (1-jndmat)
	newgraph <- graph_from_adjacency_matrix(adj, mode="undirected")
	V(newgraph)$name <- dat$illustrator_code
	if (plot) {
		plot(newgraph)
	}
	components(newgraph)$no
}




