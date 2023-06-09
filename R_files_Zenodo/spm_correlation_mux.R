spm_correlation_mux <- function(g, p_H = 0.25, k = NA){
  # Make weights 1 to ensure functionality
  E(g)$weight <- 1
  
  # Network stats and adjacency matrix
  N <- vcount(g)
  L <- ecount(g)
  L_A <- ecount(global_gA)
  A <- as_adjacency_matrix(g, type = "both")
  AA <- as_adjacency_matrix(global_gA, type = "both")
  
  # Size of the probe set E^p
  links_to_remove <- round(L*p_H)
  links_to_remove_A <- round(L_A*global_mux_pH)
  
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
  
  
  to_remove_A1 <- sample(L_A, links_to_remove_A)
  gA_perturbed_mid <- delete_edges(global_gA, to_remove_A1)
  AA_R_mid <- as_adjacency_matrix(gA_perturbed_mid, type = "both")
  dlt_AA1 <- AA - AA_R_mid
  L_A_mid <- ecount(gA_perturbed_mid)
  
  to_remove_A2 <- sample(L_A_mid, links_to_remove_A)
  gA_perturbed <- delete_edges(gA_perturbed_mid, to_remove_A2)
  
  AA_R <- as_adjacency_matrix(gA_perturbed, type = "both")
  dlt_AA2 <- AA_R_mid - AA_R
  
  
  # Diagonalise A_R
  if(is.na(k) || k == N){
    XLX <- eigen(A_R, symmetric = T)
  } else{
    XLX <- eigs_sym(A_R, k = k, which = "LA")
  }
  
  # Diagonalise AA
  if(is.na(k) || k == N){
    A_XLX <- eigen(AA, symmetric = T)
  }else{
    A_XLX <- eigs_sym(AA, k = k, which = "LA")
  }
  
  A_dlt_l1 <- rowSums((t(XLX$vectors) %*% dlt_AA1) * 
                       t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
  
  A_dlt_l2 <- rowSums((t(XLX$vectors) %*% dlt_AA2) * 
                        t(XLX$vectors)) / colSums(XLX$vectors * XLX$vectors)
  
  return(cor(A_dlt_l1, A_dlt_l2, method = "pearson"))
}