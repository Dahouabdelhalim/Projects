
#Ape and pavo are required to run these functions.
require(ape)
require(pavo)


#Functions to replicate the analyses presented by Price and Eaton (2014).
#Citation:
#Price, J. J., and M. D. Eaton. 2014. Reconstructing the evolution of sexual dichromatism: current color diversity does not reflect past rates of male and female change. Evolution 68:2026â€“2037.

#Create a function to perform ASR on all quantum catch variables, and calculate deltaS for tips and nodes. Input is a dataframe of male quantum catch values, a dataframe of female quantum catch values, and a phylogeny. The output is a list containing 1) node ancestral states - node that the male states are the first set of species names and the female states are the second set of species names 2) The node deltaS values 3) The species deltaS values.

price_eaton_asr <- function(tree,m_qcatch,f_qcatch) {
	
	mtree <- tree
	muir_region_m <- m_qcatch
	muir_region_f <- f_qcatch
	muir_region_mf <- as.data.frame(rbind(muir_region_m,muir_region_f)[,1:4])
	
	muir_u_asr_m <- ace(phy=mtree,x=as.numeric(as.character(muir_region_m$u)),type="continuous",method="REML",model="BM")
	muir_u_asr_f <- ace(phy=mtree,x=as.numeric(as.character(muir_region_f$u)),type="continuous",method="REML",model="BM")
	muir_s_asr_m <- ace(phy=mtree,x=as.numeric(as.character(muir_region_m$s)),type="continuous",method="REML",model="BM")
	muir_s_asr_f <- ace(phy=mtree,x=as.numeric(as.character(muir_region_f$s)),type="continuous",method="REML",model="BM")
	muir_m_asr_m <- ace(phy=mtree,x=as.numeric(as.character(muir_region_m$m)),type="continuous",method="REML",model="BM")
	muir_m_asr_f <- ace(phy=mtree,x=as.numeric(as.character(muir_region_f$m)),type="continuous",method="REML",model="BM")
	muir_l_asr_m <- ace(phy=mtree,x=as.numeric(as.character(muir_region_m$l)),type="continuous",method="REML",model="BM")
	muir_l_asr_f <- ace(phy=mtree,x=as.numeric(as.character(muir_region_f$l)),type="continuous",method="REML",model="BM")
	
	#Combine all ancestral states into a dataframe so we can calculate color distances
	muir_node_qc_m <- data.frame(muir_u_asr_m$ace,muir_s_asr_m$ace,muir_m_asr_m$ace,muir_l_asr_m$ace)
	rownames(muir_node_qc_m) <- paste(names(muir_u_asr_m$ace),"_m",sep="")
	colnames(muir_node_qc_m) <- c("u","s","m","l")
	muir_node_qc_f <- data.frame(muir_u_asr_f$ace,muir_s_asr_f$ace,muir_m_asr_f$ace,muir_l_asr_f$ace)
	rownames(muir_node_qc_f) <- paste(names(muir_u_asr_m$ace),"_f",sep="")
	colnames(muir_node_qc_f) <- c("u","s","m","l")
	muir_node_qc_mf <- as.data.frame(rbind(muir_node_qc_m,muir_node_qc_f))
	
	names_muir <- rownames(muir_region_mf)
	muir_region_mf <- apply(muir_region_mf,2,as.numeric)
	rownames(muir_region_mf) <- names_muir
	
	#Calculate color distances
	muir_node_coldist <- coldist(muir_node_qc_mf,qcatch="Qi",vis="tetra",noise="neural",subset=NULL,achro=T,n1=1,n2=2,n3=2,n4=4,v=0.05)
	muir_sp_coldist <- coldist(muir_region_mf,qcatch="Qi",vis="tetra",noise="neural",subset=NULL,achro=T,n1=1,n2=2,n3=2,n4=4,v=0.05)
	
	#Extract species and node values.
	muir_node_dS <- vector()
	muir_sp_dS <- vector()
	
	for (i in 1:nrow(muir_node_qc_m)){
			muir_node_dS[i] <- muir_node_coldist[intersect(grep(names(muir_u_asr_m$ace)[i],muir_node_coldist$patch1), grep(names(muir_u_asr_f$ace)[i],muir_node_coldist$patch2)),"dS"]
	}
	names(muir_node_dS) <- names(muir_u_asr_m$ace)
	
	for (i in 1:nrow(muir_region_m)){
		muir_sp_dS[i] <- muir_sp_coldist[intersect(grep(rownames(muir_region_m)[i],muir_sp_coldist$patch1), grep(rownames(muir_region_m)[i],muir_sp_coldist$patch2)),"dS"]
	}
	names(muir_sp_dS) <- rownames(muir_region_m)
	
	muir_tip_dSlabels <- muir_sp_dS[mtree$tip.label]
	names(muir_tip_dSlabels) <- names(mtree$tip.label)
	
	#Gather results
	results <- list(muir_node_qc_mf, muir_node_dS, muir_sp_dS)
	names(results) <- c("node_qc_values","node_dS","species_dS")
	
	return(results)
}


################################################################################
#This function will calculate delta S for species. It requires a data frame of quantum catch values, with species designations as row names. It can be used to calculate deltaS among species within each sex, but would require being run for males and females separately. The result is an upper triangluar data frame with all pairwise comparsions.

price_eaton_sp_dS <- function(qcvals){
	sp_names <- rownames(qcvals)
	
	#Ensure that the quantum cone catches are of type numeric
	qcvals <- as.data.frame(apply(qcvals,2,FUN= function(x) as.numeric(as.character(x))))
	rownames(qcvals) <- sp_names
	
	#Calculate the color distances
	sp_dists <- coldist(qcvals,qcatch="Qi",vis="tetra",noise="neural",subset=NULL,achro=T,n1=1,n2=2,n3=2,n4=4,v=0.05)
	
	#Aggregate results into a matrix, note that the matrix will be upper triangular, as each comparison is only available once - we use the try() function to silence the comparisons not conducted.
	res_matrix <- matrix(nrow=length(sp_names),ncol=length(sp_names))
	for (i in 1:nrow(res_matrix)){
		for (j in 1:ncol(res_matrix)){
			try(res_matrix[i,j] <- sp_dists[intersect(grep(sp_names[i],sp_dists$patch1),grep(sp_names[j],sp_dists$patch2)),"dS"],silent=TRUE)
		}
	}
	
	rownames(res_matrix) <- colnames(res_matrix) <- sp_names
	
	return(res_matrix)
}



# ################################################################################
# #The following function will calculate dS between each ancestral and daughter node pair on the phylogeny. It requires a phylo object, set of quantum cone catches for internal nodes and set of quantum catches for species. Quantum cone catches can be obtained by ancestral state reconstruction.



price_eaton_node_dS <- function(tree,node_qc,sp_qc){
	
	#First, ensure all node quantum catch vectors are numeric.
	nodenums <- rownames(node_qc)
	node_qc <- as.data.frame(apply(node_qc,2,FUN= function(x) as.numeric(as.character(x))))
	rownames(node_qc) <- nodenums
	
	#Ensure all species quantum catch vectors are numeric, and translate names into numbers associated with tree structure.
	tree_names <- tree$tip.label
	names(tree_names) <- seq(1:length(tree_names))
	tree_nums <- names(tree_names)
	names(tree_nums) <- tree_names
	sp_names <- rownames(sp_qc)
	sp_nums <- tree_nums[sp_names]
	
	sp_qc <- as.data.frame(apply(sp_qc,2,FUN= function(x) as.numeric(as.character(x))))
	rownames(sp_qc) <- sp_nums	
	
	all_qc <- as.data.frame(rbind(sp_qc,node_qc))
	
	all_dists <- coldist(all_qc,qcatch="Qi",vis="tetra",noise="neural",subset=NULL,achro=T,n1=1,n2=2,n3=2,n4=4,v=0.05)
	all_dists$patch1 <- as.character(all_dists$patch1)
	all_dists$patch2 <- as.character(all_dists$patch2)
		
	trans_res <- matrix(nrow=nrow(tree$edge),ncol=3)
	
	for (i in 1:nrow(tree$edge)){
		trans_res[i,1:2] <- tree$edge[i,]
		trans_dS <- all_dists[intersect(which(all_dists$patch1 == tree$edge[i,1]),which(all_dists$patch2 == tree$edge[i,2])),"dS"]
		if (length(trans_dS) == 0){
			trans_dS <- all_dists[intersect(which(all_dists$patch2 == tree$edge[i,1]),which(all_dists$patch1 == tree$edge[i,2])),"dS"]		}

		trans_res[i,3] <- trans_dS
		}
	colnames(trans_res) <- c("anc","dec","dS")
	trans_res <- as.data.frame(trans_res)
	
	return(trans_res)
}

