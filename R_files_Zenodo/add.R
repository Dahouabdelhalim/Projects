

add.tip.to.stem <- function(tip.name, clade.spp, phy){
  require(phytools)
  require(dplyr)
  
  if(length(clade.spp)<2){
    stop("`clade.spp` must be more than two species")
  }
  
  mono.check <- is.monophyletic(phy, clade.spp)
  if(mono.check==FALSE){
    stop("`clade.spp` are not monophyletic" )
  }
  
  clade.mrca <- findMRCA(phy, tips=clade.spp)
  # branch length leading to mrca
  parent.branch.length <- phy$edge.length[phy$edge[,2]==clade.mrca]
  # sample along that branch length
  pos <- sample(seq(from = 0.01,
                    to= parent.branch.length,
                    length.out = 100),
                1)
  out.phy <- bind.tip(tree = phy, 
                      tip.label = tip.name, 
                      where = clade.mrca,
                      position = pos)
  
  return(out.phy)
}


add.tip.to.clade <- function(tip.name, clade.spp, phy){
  require(phytools)
  require(dplyr)
  
  # is clade a single tip
  if(length(clade.spp)==1){
    target.node <- grep(TRUE, phy$tip.label==clade.spp)
    bind.information <- .generate.bind.location(target.node,
                                                phy)
    output.phy <- bind.tip(tree = phy,
                           tip.label = tip.name,
                           where = bind.information$bind.node,
                           position = bind.information$bind.location)
  } else {
  
  # The MRCA of the clade species
  clade.mrca <- findMRCA(phy, clade.spp, "node")
  # all nodes within this part of the tree (including tips)
  all.nodes <- getDescendants(phy, clade.mrca)
  # where the new tip will go
  target.node <- sample(all.nodes, 1)
  # If bind location is a tip then can only be bound above, if not, then
  # can be above or below
  bind.information <- .generate.bind.location(target.node, phy)
  # add the tip
  output.phy <- bind.tip(tree = phy, 
                         tip.label = tip.name, 
                         where = bind.information$bind.node, 
                         position = bind.information$bind.location)
  }
  
  return(output.phy)
}


.generate.bind.location <- function(node, phy){
  is.node.tip <- node <= Ntip(phy)
  if(is.node.tip==TRUE){
    edge.length <- phy$edge.length[phy$edge[,2]==node]
    bind.location <- sample(seq(from = 0.01,
                                to = edge.length,
                                length.out = 500), 1)
    output.node <- node
  } else {
    root.tip.side <- sample(c("rootward", "tipward"), 1)
    if(root.tip.side=="rootward"){
      edge.length <- phy$edge.length[phy$edge[,2]==node]
      bind.location <- sample(seq(from = 0.01,
                                  to = edge.length,
                                  length.out = 500), 1)
      output.node <- node
    } else {
      which.daughter <- .sample.daughters(phy, node)
      edge.length <- phy$edge.length[phy$edge[,2]==which.daughter]
      bind.location <- sample(seq(from = 0.01,
                                  to = edge.length,
                                  length.out = 500), 1)
      output.node <- which.daughter
    }
    
  }
  output <- list(bind.location = bind.location,
                 bind.node = output.node)
  return(output)
}

.sample.daughters <- function(phy, node){
  both.daughters <- phy$edge[phy$edge[,1]==node,2]
  sampled.daughter <- sample(both.daughters, 1)
  return(sampled.daughter)
}



