
# source functions at the bottom of this file

library(paleotree)
library(dplyr)
library(phytools)

r.vals <- c(0.2, 0.4, 0.6, 0.8)
num.repeat <- 20

result.list<- vector(mode="list", length = length(r.vals)*num.repeat)
store.location <- 1

for(k in 1:length(r.vals)){
  print(paste(k, "of", length(r.vals)))
  working.r <- r.vals[k]
  
  for(j in 1:num.repeat){
    print(paste("iteration", j))
    # Simulate a fossil record
    record <- simFossilRecord(p = 0.25,
                              q = 0.15,
                              nruns = 1,
                              nTotalTaxa = 20,
                              nExtant = 20) 
    taxa <- fossilRecord2fossilTaxa(record)
    
    # simulate a fossil record with imperfect sampling with sampleRanges
    rangesCont <- sampleRanges(taxa,
                               r = working.r) %>% as.data.frame
    # run also with rTimeRatio = 5 to simulate decreasing probability of sampling 
    # through time
    
    # update sampling in rangesCont where FAD == 0.0 Myr
    for(ee in 1:nrow(rangesCont)){
      if(rangesCont$FAD[ee]==0){
        start.val <- taxa[,"orig.time"][ee]
        seq.vals <- seq(from = 0.001, to = start.val, length.out=1000)
        rangesCont$FAD[ee] <- sample(seq.vals, 1)
      }
    }

    # let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
    cladogram <- taxa2cladogram(taxa,
                                plot = FALSE)
    
    # resolve polytomies
    cladogram.resolved <- multi2di(cladogram)
    cladogram.resolved$edge.length <- rep(1, nrow(cladogram.resolved$edge))
    write.tree(cladogram.resolved, file="backbone_tree.tre")
    
    # scale the topology using what data is available in rangesCont
    bud.config <- generate.treePL.config.cal3(phy = cladogram.resolved, 
                                dat = rangesCont[!is.na(rangesCont$FAD),], 
                                diversification = "budding", 
                                input.tree.name = "backbone_tree.tre", 
                                numsites = 1000, 
                                output.name = "scaled_tree_budding.tre")
    
    fileConn <- file("config_budding")
    writeLines(as.character(bud.config), fileConn)
    close(fileConn)
    
    bif.config <- generate.treePL.config.cal3(phy = cladogram.resolved, 
                                dat = rangesCont[!is.na(rangesCont$FAD),], 
                                diversification = "bifurcating", 
                                input.tree.name = "backbone_tree.tre", 
                                numsites = 1000, 
                                output.name = "scaled_tree_bifurcating.tre")
    
    fileConn <- file("config_bifurcating")
    writeLines(as.character(bif.config), fileConn)
    close(fileConn)
    
    # run treePL
    system(paste("treePL", "config_budding"))
    system(paste("treePL", "config_bifurcating"))
    
    timescaled.phy.budding <- read.tree("scaled_tree_budding.tre")
    timescaled.phy.bifurcating <- read.tree("scaled_tree_bifurcating.tre")
    
    ghost.lineage.calc.budding <- calculate.node.ghost.lineage(timescaled.phy = timescaled.phy.budding, 
                                                               taxa = taxa)
    ghost.lineage.calc.bifurcating <- calculate.node.ghost.lineage(timescaled.phy = timescaled.phy.bifurcating,
                                                                   taxa = taxa)
    
    itt.res <- tibble(bud.gl = sum(ghost.lineage.calc.budding$ghost.lineage),
                      bif.gl = sum(ghost.lineage.calc.bifurcating$ghost.lineage),
                      r = working.r)
    
    result.list[[store.location]] <- itt.res
    store.location <- store.location + 1
  }
  
}

saveRDS(result.list, file="result.list.RDS")

###

library(ggplot2)
library(ggthemes)

result.list <- readRDS("result.list.RDS")

plot.df <- do.call(rbind, result.list)
plot.df$diff <- plot.df$bif.gl - plot.df$bud.gl

pp <- ggplot(plot.df, aes(x = r, y = diff)) +
geom_point() +
theme_bw() +
labs(x = "Sampling rate",
  y = "Bifurcating ghost lineage - budding ghost lineage")


mean.vals <- by(plot.df$diff, plot.df$r, mean)
sd.vals <- by(plot.df$diff, plot.df$r, sd)

plot.df.alt <- data.frame(r = names(mean.vals),
  mean.vals = as.numeric(mean.vals),
  sd.vals = as.numeric(sd.vals))
plot.df.alt$lower <- plot.df.alt$mean.vals - plot.df.alt$sd.vals
plot.df.alt$upper <- plot.df.alt$mean.vals + plot.df.alt$sd.vals

pp <- ggplot(plot.df.alt, aes(x = r, y = mean.vals)) +
#geom_point() +
geom_hline(yintercept = 0, linetype="dashed") + 
geom_pointrange(aes(x=r, y=mean.vals, ymin=lower, ymax=upper)) + 
theme_bw() +
labs(x = "Sampling rate",
  y = "Bifurcating ghost lineage - budding ghost lineage")

##### functions #####


calculate.node.ghost.lineage <- function(timescaled.phy, taxa){
  all.nodes <- seq(from = Ntip(timescaled.phy)+1,
                   by = 1,
                   length.out = Nnode(timescaled.phy))
  gl.results <- lapply(all.nodes, .single.node.gl, timescaled.phy = timescaled.phy, taxa=taxa) %>% bind_rows
  return(gl.results)
}


.single.node.gl <- function(node, timescaled.phy, taxa){
  node.age <- max(nodeHeights(timescaled.phy)) - nodeheight(timescaled.phy, node)
  node.daughters <- .get.daughters(node, timescaled.phy)
  
  # First daughter
  first.daughter.gl <- .get.ghost.lineage(daughter.node = node.daughters[1],
                                          timescaled.phy = timescaled.phy,
                                          taxa = taxa,
                                          target.node.age = node.age)
  # Second daughter
  second.daughter.gl <- .get.ghost.lineage(daughter.node = node.daughters[2],
                                           timescaled.phy = timescaled.phy,
                                           taxa = taxa,
                                           target.node.age = node.age)
  
  total.gl <- first.daughter.gl + second.daughter.gl
  output <- tibble(node = node,
                   ghost.lineage = total.gl)
  return(output)
}


.get.daughters <- function(node, timescaled.phy){
  timescaled.phy$edge[timescaled.phy$edge[,1]==node,2] %>% return
}

.get.ghost.lineage <- function(daughter.node, timescaled.phy, taxa, target.node.age){
  taxa <- as.data.frame(taxa)
  is.daughter.node.a.tip <- daughter.node <= Ntip(timescaled.phy)
  if(is.daughter.node.a.tip==TRUE){
    daughter.orig.time <- taxa$orig.time[rownames(taxa)==timescaled.phy$tip.label[daughter.node]]
  } else {
    node.desc <- getDescendants(timescaled.phy, daughter.node)
    tip.desc <- node.desc[node.desc<=Ntip(timescaled.phy)]
    daughter.orig.time <- taxa$orig.time[rownames(taxa)%in%timescaled.phy$tip.label[tip.desc]] %>% max
  }
  ghost.lineage <- target.node.age - daughter.orig.time
  # ghost lineage can't be negative
  if(ghost.lineage <0 ){
    ghost.lineage <- 0
  }
  return(ghost.lineage)
}



generate.treePL.config.cal3 <- function(phy, 
                                        dat, 
                                        diversification="budding", 
                                        input.tree.name,
                                        numsites, 
                                        output.name){
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
       
        target.ages <- as.data.frame(target.ages)
        node.text[i] <- paste(as.character(round(target.ages$LAD, 2)), as.character(round(target.ages$FAD, 2)))
        conflict.check.df[i,"min"] <- target.ages[,"LAD"]
        conflict.check.df[i,"max"] <- target.ages[,"FAD"]
        
        
        # Define the calibration
        config.file.line.text[[store.location]] <- paste("mrca = ",node.name, tip.names.single.vec)
        store.location <- store.location+1
        # min age for calibration
        #config.file.line.text[[store.location]] <- paste("min = ", node.name, target.ages[,"LAD"])
        config.file.line.text[[store.location]] <- paste("min = ", node.name, target.ages[,"FAD"]-0.01)
        store.location <- store.location+1
        # max age for calibration
        config.file.line.text[[store.location]] <- paste("max = ", node.name, target.ages[,"FAD"])
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
      
      test1 <- conflict.check.df[ii,"max"] > daughter.maxs # needs to be true
      
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
    sub.dat <- dat[rownames(dat) %in% tip.names, ] %>% as.data.frame
    sub.dat <- sub.dat[!is.na(sub.dat$FAD),]
    #sub.dat <- sub.dat[sub.dat$base_stage_middle_Ma == max(sub.dat$base_stage_middle_Ma),]
    sub.dat <- sub.dat[sub.dat$FAD == max(sub.dat$FAD),]
  })
  return(res)
}


.choose.process <- function(daughter.ages, choice){
  #age.test <- daughter.ages[[1]]$base_stage_middle_Ma > daughter.ages[[2]]$base_stage_middle_Ma
  age.test <- daughter.ages[[1]]$FAD > daughter.ages[[2]]$FAD
  vec <- c(1, 2)
  if(age.test==TRUE){
    names(vec) <- c("older", "younger")
  } else {
    names(vec) <- c("younger", "older")
  }
  
  if(choice=="younger"){
    vec[names(vec)=="younger"] %>% return
  } else {
    vec[names(vec)=="older"] %>% return
  }
}


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




