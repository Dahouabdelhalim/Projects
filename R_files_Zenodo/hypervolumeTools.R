## HYPERVOLUME TOOLS: Set of tools that expands functionality of hypervolume package
# Written by: Jun Ying Lim and Susan Kennedy
# From publication:
# Kennedy, S. R., Lim, J. Y., Clavel, J., Krehenwinkel, H., and Gillespie, R. G. 2019. Spider webs, stable isotopes and molecular gut content analysis: Multiple lines of evidence support trophic niche differentiation in a community of Hawaiian spiders. Functional Ecology.

generateHypervolume <- function(data, id.col){
	# For a set of groups (e.g., species), create hypervolumes for each group
	#
	# Args:
	# 	data = data.frame containing continuous variables from which to create hypervolumes
	# 	id.col = column in `data` specifying groupings
	#
	# Returns:
	# 	a list of hypervolume objects for each group (= species)
	
	idlist <- as.vector(unique(data[,id.col]))
	
	# Estimate bandwidth for kernel density
	#est_bw <- estimate_bandwidth(data[,! names(data) %in% id.col])
	
	# Generate hypervolume for each unique id (i.e., species)
	res <- list()
	for(i in 1:length(idlist)){
	  res[[idlist[i]]] <- hypervolume(data = data[data[id.col] == idlist[i], ! names(data) %in% id.col])
		#res[[idlist[i]]] <- hypervolume(data = data[data[id.col] == idlist[i], ! names(data) %in% id.col],
	  #								bandwidth = est_bw)
										
	}
	
	# Return results
	return(res)
}

pairwiseHypervolumeOverlap <- function(list, method){
	# Generate a pairwise matrix of hypervolume overlaps given hypervolume objects in a list
	#
	# Args:
	# 	list = list of hypervolume objects
  #   method = "sorensen", "jaccard"
	#
	# Returns:
	#	diagonal matrix of pairwise sorenson overlap
	hypervol_dim <- length(list)
	res <- matrix(NA, nrow = hypervol_dim, ncol = hypervol_dim)
	colnames(res) <- names(list)
	rownames(res) <- names(list)
	
	# Create pairs
	pairs <- t(combn(1:hypervol_dim, 2))
	
	# Calculate overlap for pairs
	overlap <- vector()
	for(i in 1:nrow(pairs)){
		hv1 <- list[[pairs[i,][1]]]
		hv2 <- list[[pairs[i,][2]]]
		set <- hypervolume_set(hv1, hv2, check.memory = FALSE, num.points.max = 10^6)
		overlap[i] <- hypervolume_overlap_statistics(set)[method]
	}
	res[lower.tri(res)] <- overlap
	return(res)
}