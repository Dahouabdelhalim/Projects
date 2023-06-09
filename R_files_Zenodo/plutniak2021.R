# install and load required packages ####
if (! requireNamespace("igraph", quietly = TRUE)){
  install.packages("igraph")
}
if (! requireNamespace("sna", quietly = TRUE)){
  install.packages("sna")
  }
if (! requireNamespace("RBGL", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("RBGL")
}
if (! requireNamespace("qpgraph", quietly = TRUE)){
  BiocManager::install("qpgraph")
}

library(igraph)
library(sna)
library(RBGL)
library(qpgraph)

# Harary and Ross 1951 implementation ####
.getC <- function(M, unicliquals, Mrowsums){
  Cgroup.list <- list()
  CgroupIndex.list <- list()
  CgroupPrime.list <- list()
  for(i in 1:nrow(M) ){
    # skip if the point already belongs to a clique:
    if(i %in% unlist(CgroupIndex.list)) next  
    # skip if no unicliqual points:
    if(! i %in% unicliquals ) next
    CgroupIndex <- sort(c(i, which(M[, i] != 0)))
    Cgroup <- sort(colnames(M)[c(i, which(M[, i] != 0))])
    CgroupPrime <- which( Mrowsums == Mrowsums[i] ) 
    CgroupPrime <- CgroupPrime[CgroupPrime %in% CgroupIndex]
    # add to results 
    Cgroup.list <- append(Cgroup.list, list(Cgroup))
    CgroupIndex.list <- append(CgroupIndex.list, list(CgroupIndex))
    CgroupPrime.list <- append(CgroupPrime.list, list(CgroupPrime))
  }
  list(Cgroup.list = Cgroup.list,
       CgroupPrime.list = CgroupPrime.list)
}

.substract.matrix <- function(mat, CgroupPrime.list){
  i <- 1:nrow(mat)
  i <- i[ ! 1:nrow(mat) %in% unlist(CgroupPrime.list) ]
  mat[i, i]
}

.extract.cliques <- function(mat){
  cliques.list <- list() 
  M <- mat * mat %*% t(mat)
  np <- apply(mat, 1, sum) # similar to vertices' degrees
  Mrowsums <- rowSums(M)
  unicliquals <- which( Mrowsums == np * (np - 1) )
  
  if(length(unicliquals) > 0){
    res <- .getC(M, unicliquals, Mrowsums)
    cliques <- res$Cgroup.list
    # get only cliques with at least 3 vertices:
    cliques <- cliques[sapply(cliques, function(x) length(x) > 2 )]
    cliques.list <- append(cliques.list, cliques)
    mat <- .substract.matrix(mat, res$CgroupPrime.list)
  }
  else{
    CgroupIndex <- sort(c(1, which(M[, 1] != 0)))
    Cgroup <- c(1, which(mat[, 1] != 0))
    submat1 <- mat[Cgroup, Cgroup]
    submat2 <- mat[-1, -1]
    mat <- list(submat1, submat2)
  }
  list(mat, cliques.list)
}

haross.cliques <- function(mat){
  # initial tests:
  if( ! is.matrix(mat) ) stop("The argument is not a matrix.") 
  if( ncol(mat) != nrow(mat) ) stop("A square matrix is required.") 
  if( is.null(colnames(mat)) & is.null(rownames(mat)) ){
    colnames(mat) <- 1:ncol(mat)
    rownames(mat) <- 1:nrow(mat)
  }
  # set variables:
  cliques.list.final <- list() 
  mat.list <- list(mat)
  
  repeat{ # repeat while the sum of the matrix values > 0
    # run the main function:
    res <- lapply(mat.list, .extract.cliques)
    # sort results:
    #   1) extract and add the cliques to the list:
    cliques.list.final <- append(cliques.list.final, 
                                 lapply(res, function(x) x[[2]])
    )
    #   2) extract the list of matrices:
    mat.list <- lapply(res, function(x) x[[1]] )
    # if the list is too nested, unnest:
    if( is.list(mat.list[[1]]) & length(mat.list[[1]]) > 1 ) {
      mat.list <- unlist(mat.list, recursive = F)
    }
    # keep only the matrices with more than 2 points:
    mat.list <- mat.list[ sapply(mat.list, function(x) sum(x) > 2 ) ]
    # if there is no more matrices, break
    if( length(mat.list) == 0 ) break
  }
  
  unlist(cliques.list.final, recursive = F)
}

# Make table from Gardin & Garelli 1961, fig. 10 p. 868: ####
GardinGarelli.df <- as.matrix(rbind(
  c("amur-ishtar", "laqipun"),
  c("hina", "amur-ishtar"),
  c("hina", "im(i)d-ilum"),
  c("hina", "laqipun"),
  c("im(i)d-ilum", "amur-ishtar"),
  c("im(i)d-ilum", "laqipun"),
  c("pushu-kin", "hina"),
  c("pushu-kin", "laqipun"),
  c("amur-ishtar", "assur-nada"),
  c("assur-imitti", "amur-ishtar"),
  c("assur-imitti", "assur-nada"),
  c("assur-taklaku", "amur-ishtar"),
  c("assur-taklaku", "assur-nada"),
  c("assur-taklaku", "pushu-kin"),
  c("pushu-kin", "assur-nada"),
  c("assur-tab", "im(i)d-ilum"),
  c("assur-tab", "pushu-kin"),
  c("buzazu", "enna-sin"),
  c("buzazu", "shu-belim"),
  c("mannum-balum-assur", "enna-sin"),
  c("mannum-balum-assur", "shu-belim"),
  c("amur-ishtar", "pushu-kin"),
  c("im(i)d-ilum", "pushu-kin"),
  c("enna-sin", "shu-belim")
))
colnames(GardinGarelli.df) <- c("from", "to")

# cliques given by Gardin and Garelli 1961:
GardinGarelli.res <- list(
     c("amur-ishtar", "pushu-kin", "hina", "im(i)d-ilum", "laqipun"),
     c("amur-ishtar", "pushu-kin", "assur-nada", "assur-taklaku"),
     c("amur-ishtar", "assur-imitti", "assur-nada"),
     c("pushu-kin", "assur-tab", "im(i)d-ilum"),
     c("buzazu", "enna-sin", "shu-belim"),
     c("enna-sin", "mannum-balum-assur", "shu-belim")
)

# Generate different formats: ####
# matrix:
GardinGarelli.g <- igraph::graph_from_data_frame(GardinGarelli.df, directed=F)
GardinGarelli.df <- igraph::as_adjacency_matrix(GardinGarelli.g, sparse=F)
# graphNEL:
GardinGarelli.nel <- igraph::igraph.to.graphNEL(GardinGarelli.g)
# sna network:
GardinGarelli.net <- network::network(GardinGarelli.df, directed=F)


# Clique algorithms: ####

## .... Harary and Ross 1957 ####
harary.res <- haross.cliques(GardinGarelli.df)
harary.res <- harary.res[ order(sapply(harary.res, length), decreasing=T) ]

## .... Bron and Kerbosch 1973 (RBGL) ####
bron.res <- RBGL::maxClique(GardinGarelli.nel)$maxCliques

## .... Makino and Uno 2004 (sna) ####
makino.res <- sna::clique.census(GardinGarelli.net, mode="graph",
                                 tabulate.by.vertex=F)$cliques
makino.res <- lapply(makino.res, function(x)
  lapply(x, function(y) network.vertex.names(GardinGarelli.net)[y]) )
makino.res <- unlist(makino.res, recursive = F)

## .... Ostergard 2001 (qpgraph) ####
ostergard.res <- qpgraph::qpGetCliques(GardinGarelli.df)
ostergard.res <- lapply(ostergard.res,
                        function(i) rownames(GardinGarelli.df)[i]) 

## .... Eppstein et al. 2010 (igraph)  ####
eppstein.res <- igraph::max_cliques(GardinGarelli.g, min=3)
eppstein.res <- lapply(eppstein.res, names)

# Compare results: ####

res.list <- list(
  "Gardin 1961" = GardinGarelli.res,
  "Harary 1957" = harary.res,
  "Bron 1973" = bron.res,
  "Makino 2004" = makino.res,
  "Osertgard 2001" = ostergard.res,
  "Eppstein 2010" = eppstein.res
)

res.list <- lapply(res.list, function(x) lapply(x,  sort))
res.list <- lapply(res.list, function(x) lapply(x,  paste, collapse="/"))

res.tab <- lapply(res.list, function(x)  res.list[[2]] %in% x)
res.tab <- do.call("rbind", res.tab)
res.tab <- t(res.tab)

res.tab <- cbind(id = c(1:6, NA, NA),
      size = sapply(harary.res, length), 
      res.tab)
res.tab

