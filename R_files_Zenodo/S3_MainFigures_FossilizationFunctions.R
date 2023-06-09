####PRESERVATION FUNCTIONS

# Load all the packages we need
if(!require(magrittr)){install.packages("magrittr")}
if(!require(reshape2)){install.packages("reshape2")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(dplyr)){install.packages("dplyr")}
if(!require(tidyr)){install.packages("tidyr")}
if(!require(stringr)){install.packages("stringr")}



# set the p values
set.ps <- function(p.eql, p.hard, p.med, p.soft){
  p <- vector()
  
  p["eql"]  <- p.eql # all species have equal probability of fossilizing
  
  p["hard"] <- p.hard # probability that something hard bodied fossilizes
  p["med"]  <- p.med # probability that something midway fossilizes
  p["soft"] <- p.soft # probability that something soft bodied fossilizes
  
  
  return(p)
}

# replace spaces in species names with .
fixnames.fw <- function(fw){
  rownames(fw) <-gsub("[[:space:]]", ".", rownames(fw))
  rownames(fw) <- str_replace_all(rownames(fw), "^[:punct:]", "X.")
  rownames(fw) <- str_replace_all(rownames(fw), "[:punct:]", ".")
  #Trailing punctuation
  rownames(fw)<-gsub('^\\\\.|\\\\.$', '',rownames(fw))
  
  colnames(fw) <-gsub("[[:space:]]", ".", colnames(fw))
  colnames(fw) <- str_replace_all(colnames(fw), "^[:punct:]", "X.")
  colnames(fw) <- str_replace_all(colnames(fw), "[:punct:]", ".")
  #Trailing punctuation
  colnames(fw)<-gsub('^\\\\.|\\\\.$', '',colnames(fw))
  
  # Check which names don't match
  #colnames(fw)[which(!colnames(fw) %in% rownames(fw))]
  #rownames(fw)[which(!rownames(fw) %in% colnames(fw))]
  
  if(length(colnames(fw)[which(colnames(fw) %in% rownames(fw))]) != nrow(fw)){
    stop("Colnames don't match rownames")
  }
  return(fw)
}

# replace spaces in species names with .
fixnames.nodes <- function(nodes){
  nodes$Species <- gsub("[[:space:]]", ".", nodes$Species)
  nodes$Species <- gsub(",", "X.", nodes$Species)
  nodes$Species <- str_replace_all(nodes$Species, "[:punct:]", ".")
  nodes$Species<-gsub('^\\\\.|\\\\.$', '',nodes$Species)
  
  return(nodes)
}


# assign traits
library(stringr)
traits <- function(nodes,name){
  hs_meta<-read.csv("S3_PresGrAssignments_updated_200902.csv",stringsAsFactors=F)
  names(hs_meta)[names(hs_meta) == 'BodyType'] <- 'HardSoft'
  hs_meta$node_name<-gsub("[[:space:]]", ".", hs_meta$node_name)
  hs_meta$node_name <- str_replace_all(hs_meta$node_name, "^[:punct:]", "X.")
  hs_meta$node_name <- str_replace_all(hs_meta$node_name, "[:punct:]", ".")
  
  nodes$Food_web<-name
  nodes$node_name<-nodes$Species
  # add hard/soft data
  nodes <- left_join(nodes, hs_meta[,c("node_name","HardSoft","Food_web")], by =c("node_name","Food_web"))
  nodes[is.na(nodes$HardSoft),"HardSoft"] <- 2
  
  return(nodes)
}


# assign probabilities to each species based on their traits
assign.p <- function(nodes, p.eql = 0.5, p.hard = 0.75, p.med = 0.5, p.soft = 0.25, vary=FALSE){
  
  
  p <- set.ps(p.eql, p.hard, p.med, p.soft)
  # assign fossilization probability (assuming 1 = soft, 2 = med, 3 = hard)
  nodes<-as.data.frame(nodes)
  nodes[,"p.eql"] <- p["eql"]
  
  if(vary==TRUE){
    #Vary fossilization probabilities
    nodes[,"p.hs"] <- ifelse(nodes[,"HardSoft"] == 1, sample(x=c(1:33333),size=nrow(nodes),replace=T)/100000, 
                         ifelse(nodes["HardSoft"] == 2, sample(x=c(33333:66666),size=nrow(nodes),replace=T)/100000,
                                sample(x=c(66666:99999),size=nrow(nodes),replace=T)/100000))
  }else{
    nodes[,"p.hs"] <- ifelse(nodes[,"HardSoft"] == 1, p["soft"], ifelse(nodes[,"HardSoft"] == 2, p["med"], p["hard"]))
  }
  
  nodes[,"p.mean"] <- mean(nodes[,"p.hs"])
  
  return(nodes)
  
}


# make a "fossil" network
fossilize <- function(data, runs, trmnt, E.shape1 = 0, E.shape2 = 1, E.distr = "uniform", null = FALSE) {
  fossils <- data.frame(species = data$Species)
  
  for (i in 1:runs) {
    keep1 <- data
    
    keep1$samples <- rbeta(nrow(keep1), shape1 = E.shape1, shape2 = E.shape2)
    
    if (null == "null.shuffle") {
      keep1[, trmnt] <- sample(keep1[, trmnt])
    } else {
      
    }
    
    fossils[, paste0("run", i)] <- ifelse(keep1[, "samples"] < keep1[, trmnt], 1, 0)
  }
  
  
  return(fossils)
}



# recreate the food web based on the species present
recreate.web <- function(fossilnodes, orig.web, runs) {
  fossilwebs <- list()
  #fossilwebs[[1]] <- make.trophic.species.web(orig.web)
  fossilwebs[[1]] <- orig.web
  
  # Make names into indices and find primary producers
  orig.web2 <- orig.web
  rownames(orig.web2) <- colnames(orig.web2)
  
  # Find primary producers
  primprods <- orig.web2[, colSums(orig.web2) == 0]
  primprods <- (as.character(colnames(primprods)))
  
  
  for (run in 1:runs) {
    fossilwebs[[run + 1]] <- orig.web2
    
    ## Should be able to streamline this if we're only using square adjacency matrices
    sp.rmv <- fossilnodes[which(fossilnodes[, paste0("run", run)] == 0), "species"]
    row.rmv <- which(row.names(fossilwebs[[run + 1]]) %in% sp.rmv)
    col.rmv <- which(colnames(fossilwebs[[run + 1]]) %in% sp.rmv)
    
    if (length(col.rmv) > 0) {
      fossilwebs[[run + 1]] <- fossilwebs[[run + 1]][-row.rmv, -col.rmv]
    } else if (length(col.rmv) == 0 & length(row.rmv) > 0) {
      fossilwebs[[run + 1]] <- fossilwebs[[run + 1]][-row.rmv, ]
    } else if (length(row.rmv) == 0 & length(col.rmv) == 0) {
      fossilwebs[[run + 1]] <- fossilwebs[[run + 1]]
    } else {
      print("Something's wrong in recreate.web()")
    }
    
    # Knock out isolated nodes
    if (length(fossilwebs[[run + 1]]) > 1) {
      
      # Make names into indices
      fosweb <- fossilwebs[[run + 1]]
      fosweb2 <- fosweb
      
      # Identify primary producers in fossil web
      primprods2 <- intersect(primprods, rownames(fosweb2))
      fosweb2 <- fosweb2[primprods2, ]
      # Identify PPs with no consumers
      #primprods3 <- fosweb2[rowSums(fosweb2) == 0, ]
      primprods3 <- ifelse(nrow(fosweb2)==NULL, ,fosweb2[rowSums(fosweb2) == 0, ])
      # Delete primprods3 from fossil web
      primprods3 <- (as.character(rownames(primprods3)))
      
      # Identify consumers without resources in fossil web
      cons <- rownames(fosweb)
      cons2 <- cons[!cons %in% primprods]
      cons3 <- fosweb[,cons2]
      cons4 <- ifelse(nrow(cons3)==NULL, ,cons3[, colSums(cons3) == 0])
      cons4 <- as.character(colnames(cons4))
      
      # Only drop PPs if there is something to drop
      if (length(primprods3) > 0) {
        fosweb3 <- fossilwebs[[run + 1]]
        fossilwebs[[run + 1]] <- fosweb3[!(rownames(fosweb3) %in% primprods3), !(colnames(fosweb3) %in% primprods3)]
      } else {
        fossilwebs[[run + 1]] <- fosweb
      }
      
      # Only drop consumers if there is something to drop
      if (length(cons4) > 0) {
        fosweb3 <- fossilwebs[[run + 1]]
        fossilwebs[[run + 1]] <- fosweb3[!(rownames(fosweb3) %in% cons4), !(colnames(fosweb3) %in% cons4)]
      } else {
        fossilwebs[[run + 1]] <- fossilwebs[[run + 1]]
      }
    } else {
    }
    
    #Make trophic web & drop isolated nodes
    if (length(fossilwebs[[run + 1]]) > 1){
      Isolated = which(degree(graph_from_adjacency_matrix(as.matrix(fossilwebs[[run+1]])))==0)
      G2 = delete.vertices(graph_from_adjacency_matrix(as.matrix(fossilwebs[[run+1]])), Isolated)
      try<-as.data.frame(as.matrix(get.adjacency(G2)))
      fossilwebs[[run+1]]<-try
      #fossilwebs[[run + 1]]<-make.trophic.species.web(fossilwebs[[run + 1]])
    }else{
    }
    
  }
  return(fossilwebs)
}

# Calculate metrics on the fossilized food webs

calc.metrics <- function(fw.list,runs){
  

  swtl<-sapply(fw.list,calc.mean.swtl)

  SOI<-sapply(fw.list,calc.SOI)
  diam<-sapply(fw.list,calc.diam)
  C<-sapply(fw.list,calc.C)

  S<-sapply(fw.list,calc.S)
  degree<-sapply(fw.list,calc.norm.mean.degree)

  cc<-sapply(fw.list,calc.cc)

  avpath<-sapply(fw.list,calc.mean.path.length)

  
  metrics <- as.data.frame(cbind(
    swtl,

    SOI,
    diam,
    C,
  
    S,
    degree,
    cc,

    avpath

  ))
  
  nonOG<-suppressWarnings(mutate_all(metrics[2:(runs+1),], function(x) as.numeric(as.character(x))))
  nonOGmean <- colMeans(nonOG, na.rm = TRUE)
  nonOGsd <- apply(nonOG,2,sd, na.rm = TRUE)
  metrics[(nrow(metrics)+1),] <- nonOGmean
  metrics[(nrow(metrics)+1),] <- nonOGsd
  metrics$run <- 0:(nrow(metrics)-1)   
  metrics$type<-"Fossilized"
  metrics[1,"type"] <- "Original"
  metrics[nrow(metrics)-1,"type"] <- "Mean"
  metrics[nrow(metrics),"type"] <- "SD"
  
  return(metrics)
}



#Set node and fw names to be the same
set.sp.names<-function(fwlist) {
  fw_fw<-fwlist[[1]]
  fw_nodes<-fwlist[[2]]
  
  fw_nodes$Species<-rownames(fw_fw)
  colnames(fw_fw)<-rownames(fw_fw)
  
  new_fw_list<-list(FW=fw_fw,NODE=fw_nodes)
  
  return(new_fw_list)
}


# Create niche model webs - adapted from trophic package
niche2<-function(S,C) {
  
  if(S>1 & C>0 & is.finite(S) & is.finite(C)) {
    
    n.i <- sort(runif(S), decreasing = F)
    r.i <- suppressWarnings(rbeta(S, 1, ((1/(2 * C)) - 1)) * n.i)
    c.i <- suppressWarnings(runif(S, r.i/2, n.i))
    
    if(any(is.na(c.i))==FALSE & any(is.na(r.i))==FALSE & any(is.na(n.i))==FALSE) {
      
      a <- matrix(0, nrow = S, ncol = S)
      for (i in 2:S) {
        for (j in 1:S) {
          if (n.i[j] > (c.i[i] - (0.5 * r.i[i])) & n.i[j] < 
              (c.i[i] + 0.5 * r.i[i])) {
            a[j, i] <- 1
          }
        }
      }
      
    }else{ 
      
      a<-NULL
      
    }
    
    return(a)
  }else{return(NULL)}
}  



# Calculate number of hard/soft taxa in a fossilized web
calc.att.numbers <- function(fw.node.list,node.attributes,runs){
  
  att_fw.list<-merge(node.attributes,fw.node.list,by.x="Species",by.y="species",all.x=TRUE)
  att_fw.list$run0<-1
  
  numbas_list<-c()
  for(row in 0:(ncol(fw.node.list)-1)){
    col_run<-paste("run",row,sep="")
    att_sub<-att_fw.list[,c(col_run,"HardSoft")]
    att_sub$selected_run<-att_sub[,c(col_run)]
    att_sub<-subset(att_sub, selected_run==1)
    
    att_sub$HardSoft<-factor(att_sub$HardSoft, levels = 1:3)
    
    if(nrow(att_sub)>0){
      numbas<-as.matrix(rbind(cbind(as.data.frame(table(att_sub$HardSoft)),att="HardSoft")
      ))
      numbas_list<-rbind(numbas_list,cbind(numbas,run=row))
      
    }else{
      numbas<-as.matrix(rbind(c(1,0,"HardSoft"),
                              c(2,0,"HardSoft"),
                              c(3,0,"HardSoft")))
      numbas_list<-rbind(numbas_list,cbind(numbas,run=row))
      
    }
    
    
  }
  return(numbas_list)
}









