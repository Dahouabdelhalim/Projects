install.packages("phangorn")
install.packages("phytools")
install.packages("phylobase")
library(phangorn)
library(phytools)
library(phylobase)



###################################################################################################

#Function to compute the branch length distance and binary character difference for each pair of taxa within a tree
#tree = an object of class phylo
#morph = a binary character matrix of class phydat

coherence <- function(tree, morph){
	char_n <- length(morph[[1]])
	TN <- Ntip(tree)-1
	listOne <- tree$tip.label[c(1:TN)]
	taxon_distances <- data.frame(Taxon_A = 0, Taxon_B = 0, seq_distance = 0, BL_distance = 0)
	for(x in 1:length(listOne)){
		listTwo <-  tree$tip.label[c((x+1):length( tree$tip.label))]
		for(i in 1:length(listTwo)){
			seq_diff <- as.numeric(unlist( morph[(listOne[x])])) - as.numeric(unlist( morph[(listTwo[i])]))
			seq_distance <- length(seq_diff) - sum(seq_diff == 0)
			MRC_ancestor <- MRCA( tree, c(listOne[x], listTwo[i]))
			BL_distance <- (nodeheight( tree, match(listOne[[x]],  tree$tip.label)) - nodeheight( tree, MRC_ancestor)) + (nodeheight( tree, match(listTwo[[i]],  tree$tip.label)) - nodeheight( tree, MRC_ancestor))
			newrow <- data.frame(Taxon_A = listOne[x], Taxon_B = listTwo[i], seq_distance = seq_distance, BL_distance = BL_distance)
			taxon_distances <- rbind(taxon_distances, newrow)}}		
	taxon_distances <- taxon_distances[-c(1), ]
	return(taxon_distances)}
	
###################################################################################################
	
#example script to demonstarte compution of raw and adjusted morphological coherence
	#create a random tree
		tree1 <- rtree(30)
	#create a random binary character matrix comprising 200 characters
		morph.m <- matrix(0, Ntip(tree1), 200)
		rownames(morph.m) <- tree1$tip.label
		for(a in 1:Ntip(tree1)){
			for(b in 1:200){
				morph.m[[a, b]] <- sample(0:1, 1)
			}
		}
		phydat_m <- as.phyDat(morph.m, type="USER", levels = c(0, 1))
	#compute taxon pair distances
		distances <- coherence(tree1, phydat_m)
	#compute Spearmans Rank correlation (raw morphological coherence)
		cor(distances$seq_distance, distances$BL_distance, method = "spearman")
	#strip taxon pairs with long branches
		distances2 <- subset(distances, BL_distance <= (max(distances$BL_distance))/2, select = )
	#compute Spearmans Rank correlation (adjusted morphological coherence)
		cor(distances2$seq_distance, distances2$BL_distance, method = "spearman")
		
###################################################################################################

#Function to compute stemminess of a tree (Fiala & Sokal, 1985)

#tree = an object of class phylo

stemminess <- function(tree){
	Nodes <- Nnode(tree)
	Tips <- Ntip(tree)
	Edges <- tree$edge.length
	tree2 <- phylo4(tree)
	sim_df <- head(tree2, (Tips+Nodes+1))
	stemminess <- data.frame(node = 0, values = 0)
	for(i in (Tips + 2):(Tips + Nodes)){
		dec <- descendants(tree2, i, "all")
		values <- sim_df$edge.length[[i]]
		for(j in dec){
			values <- c(values, sim_df$edge.length[[j]])
			}
		newrow <- data.frame(node = i, values = values[[1]]/sum(values))
		stemminess <- rbind(stemminess, newrow)}
	stemminess <- stemminess[-c(1), ] 
	mean_s <- mean(stemminess$values)	
	return(mean_s)}
		
###################################################################################################

#function to calculate the number of bipartitions shared between 2 unrooted trees. Symmetric metric. For unrooted trees, shared nodes = shared bibartitions + 1. 

SB <- function(tree1, tree2){
	tree1 <- unroot(tree1)
	tree2 <- unroot(tree2)
	RF <- RF.dist(tree1, tree2)
	part_all <- (Nnode(tree1) -1) + (Nnode(tree2) -1)
	return((part_all - RF)/2)}

###################################################################################################

#function to calculate the number of bipartitions found in unrooted tree 1 that are not found in unrooted tree 2. Asymmetric metric. For unrooted trees, unique nodes = unique bibartitions. 
   
UB <- function(tree1, tree2){
	tree1 <- unroot(tree1)
	tree2 <- unroot(tree2)
	RF <- RF.dist(tree1, tree2)
	part1 <- (Nnode(tree1) -1)
	part2 <- (Nnode(tree2) -1)
	part_all <- (Nnode(tree1) -1) + (Nnode(tree2) -1)
	SB <- ((part_all - RF)/2)
	return(part1 - SB)}
