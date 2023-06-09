

generate.treePL.config <- function(phy, dat, diversification="budding", input.tree.name,
                                   numsites, output.name){
  require(phytools)
  
  ## Information for the start of the configuration file 
  config.file.line.text <- vector(mode="list", length=4)
  
  # name of the phylogeny file being scaled
  config.file.line.text[[1]] <- paste("treefile =", input.tree.name)
  
  # smoothing parameter
  config.file.line.text[[2]] <- "smooth = 100"
  
  # number of sites in their alignment
  config.file.line.text[[3]] <- paste("numsites =", numsites)
  
  
  ## Now iterate over every node in the tree, deriving the age from either the older
  ## or younger of the descendants dependent on the diversification process
  store.location <- 4
  
  all.nodes <- seq(from = Ntip(phy)+1, length.out = Nnode(phy), by = 1)
  
  # Can be used to plot node ages on the phylogeny
  node.text <- vector()
  
  # Node conflict check data.frame
  conflict.check.df <- data.frame(node = all.nodes, min = NA, max = NA)
  
  for(i in 1:length(all.nodes)){
    
    node <- all.nodes[i]
    # Used to name each calibration
    node.name <- paste0("node.",node) 
    
    num.daughters <- .num.daughters(node, phy)
    if(num.daughters==2){
      
      # The tips associated with the node
      node.desc <- getDescendants(tree = phy, node = node)
      tips.in.calibration <- node.desc[node.desc<=Ntip(phy)]
      tip.names <- phy$tip.label[tips.in.calibration]
      tip.names.single.vec <- paste0(tip.names, collapse = " ")
      
      # Which daughter defines the age
      daughter.ages <- .generate.daughter.ages(node, phy, dat)
      
      if(nrow(daughter.ages[[1]])>0 & nrow(daughter.ages[[2]])>0){
        
        if(diversification=="budding"){
          target.ages <- daughter.ages[[.choose.process(daughter.ages, "younger")]][1,]
        } else {
          target.ages <- daughter.ages[[.choose.process(daughter.ages, "older")]][1,]
        }
        
        node.text[i] <- paste(as.character(target.ages$base_stage_upper_Ma), as.character(target.ages$base_stage_lower_Ma))
        conflict.check.df[i,"min"] <- target.ages[,"base_stage_upper_Ma"]
        conflict.check.df[i,"max"] <- target.ages[,"base_stage_lower_Ma"]
        
        
        # Define the calibration
        config.file.line.text[[store.location]] <- paste("mrca = ",node.name, tip.names.single.vec)
        store.location <- store.location+1
        # min age for calibration
        config.file.line.text[[store.location]] <- paste("min = ", node.name, target.ages[,"base_stage_upper_Ma"])
        store.location <- store.location+1
        # max age for calibration
        config.file.line.text[[store.location]] <- paste("max = ", node.name, target.ages[,"base_stage_lower_Ma"])
        store.location <- store.location+1
      }
    }
  }
  
  ## Remove conflicts
  conflict.check.df$checkResults <- .conflict.checks(conflict.check.df, phy)
  
  nodes.to.remove <- paste0("node.", conflict.check.df$node[conflict.check.df$checkResults=="conflict"])
  if(length(nodes.to.remove)>0){
    lines.to.drop <- vector()
    for(jj in 1:length(nodes.to.remove)){
      lines.to.drop <- c(lines.to.drop, grep(nodes.to.remove[jj], config.file.line.text))
    }
    config.file.line.text <- config.file.line.text[-lines.to.drop]
  }
  
  
  ## Now we just need to add in the name of the output file from the scaling
  length.file <- length(config.file.line.text)
  final.entry <- length.file + 1
  
  config.file.line.text[[final.entry]] <- paste("outfile =", output.name)
  
  return(config.file.line.text)
}


.conflict.checks <- function(conflict.check.df, phy){
  test.results <- lapply(1:nrow(conflict.check.df), function(ii){
    if(.is.only.node(conflict.check.df$node[ii], phy)==TRUE){
      return("safe")
    } else {
      
      
      daughters <- getDescendants(phy, conflict.check.df$node[ii])
      daughter.nodes <- daughters[daughters>Ntip(phy)]
      
      daughter.mins <- conflict.check.df[conflict.check.df$node %in% daughter.nodes, "min"]
      daughter.maxs <- conflict.check.df[conflict.check.df$node %in% daughter.nodes, "max"]
      
      test1 <- conflict.check.df[ii,"min"] > daughter.maxs # needs to be true
      
      if(FALSE %in% test1){
        return("conflict")
      } else{
        return("safe")
      }
    }
  }) %>% unlist
}


.generate.daughter.ages <- function(node, phy, dat){
  both.daughters <- phy$edge[phy$edge[,1]==node,2]
  res <- lapply(1:2, function(ee){
    desc <- getDescendants(phy, both.daughters[ee])
    tips <- desc[desc<=Ntip(phy)]
    tip.names <- phy$tip.label[tips]
    sub.dat <- dat[dat$family %in% tip.names, c("base_stage_upper_Ma","base_stage_middle_Ma", "base_stage_lower_Ma")]
    sub.dat <- sub.dat[sub.dat$base_stage_middle_Ma == max(sub.dat$base_stage_middle_Ma),]
  })
  return(res)
}


.choose.process <- function(daughter.ages, choice){
  age.test <- daughter.ages[[1]]$base_stage_middle_Ma > daughter.ages[[2]]$base_stage_middle_Ma
 if(choice=="younger"& age.test[1]==TRUE){
   return(2)
 } else {
   return(1)
 }
}

#.conflict.checks <- function(conflict.check.df, phy){
#  test.results <- lapply(1:nrow(conflict.check.df), function(ii){
#    daughters <- phy$edge[phy$edge[,1]==conflict.check.df$node[ii],2]
#    daughter.mins <- conflict.check.df[conflict.check.df$node %in% daughters, "min"]
#    test <- conflict.check.df[ii,"max"] < daughter.mins
#    if(TRUE %in% test){
#      return("conflict")
#    } else{
#      return("safe")
#    }
#  }) %>% unlist
#}


.is.only.node <- function(node, phy){
  desc <- getDescendants(phy, node)
  nodes <- desc[desc>Ntip(phy)]
  if(length(nodes)==0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}





.num.daughters <- function(node, phy){
  tt <- table(phy$edge[,1]==node)
  num.daughters <- tt[names(tt)==TRUE]
  names(num.daughters) <- NULL
  return(num.daughters)
}

