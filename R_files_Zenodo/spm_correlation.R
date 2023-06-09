spm_correlation <- function(g, p_H = 0.1, k = NA, iterationNumber = 30){
  
  SumCor = 0
  for (icount in 1:iterationNumber)
  {
  # Make weights 1 to ensure functionality
  E(g)$weight <- 1
  
  # Network stats and adjacency matrix
  N <- vcount(g)
  L <- ecount(g)
  A <- as_adjacency_matrix(g, type = "both")
  
  # Size of the probe set E^p
  links_to_remove <- round(L*p_H)
  
  # Determine the link perturbation set dlt_E at random
  to_remove1 <- sample(L, links_to_remove)
  g_perturbed_mid <- delete_edges(g, to_remove1)
  
  A_R_mid <- as_adjacency_matrix(g_perturbed_mid, type = "both")
  dlt_A1 <- A - A_R_mid
  
  L_mid <- ecount(g_perturbed_mid)
  to_remove2 <- sample(L_mid, links_to_remove)
  g_perturbed <- delete_edges(g_perturbed_mid, to_remove2)
  
  A_R <- as_adjacency_matrix(g_perturbed, type = "both")
  dlt_A2 <- A_R_mid - A_R
  
  # Compute the adjacency matrix of the graph after removal of dlt_E
  A_R <- as_adjacency_matrix(g_perturbed, type = "both")
  
  # Diagonalise A_R
  if(is.na(k) || k == N){
    XLX <- eigen(A_R, symmetric = T)
  } else{
    XLX <- eigs_sym(A_R, k = k, which = "LA")
  }
  
  # Compute the adjacency matrix of the perturbation set
  #dlt_A <- A - A_R
  #
  # dlt_A1 <- as_adjacency_matrix(subgraph.edges(g, E(g)[to_remove1]), type = "both")
  # paste(dlt_A1)
  # dlt_A2 <- as_adjacency_matrix(subgraph.edges(g, E(g)[to_remove2]), type = "both")
  # paste(dlt_A2)
  
  # Compute correction terms for eigenvalues based on dlt_A
  dlt_l1 <- rowSums((t(XLX$vectors) %*% dlt_A1) *
                      t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
  
  dlt_l2 <- rowSums((t(XLX$vectors) %*% dlt_A2) *
                      t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
  SumCor = SumCor + cor(dlt_l1, dlt_l2, method = "pearson")
  }
  return(SumCor/iterationNumber)
}