library(igraph)
form_net <- function(edgelist_filepath) {
  graph_edgetable <- read.table(file=edgelist_filepath, header=FALSE, sep="")
  max_node_number <- max(graph_edgetable)
  adj_mat <- matrix(0, max_node_number, max_node_number)
  adj_mat[cbind(source=graph_edgetable[, 1], target=graph_edgetable[, 2])] <- 1
  # adj_mat[cbind(graph_edgetable[, 2], graph_edgetable[, 1])] <- 1
  return(adj_mat)
}