require(igraph)
require(NetIndices)
require(broom)

# Trophic level

# Diameter
calc.diam <- function(fw) {
  if(length(unique(c(rownames(fw),colnames(fw)))) == 0){ diam<-NA }else{
  diam <- diameter(graph_from_adjacency_matrix(as.matrix(fw)))
  }
  return(diam)
}

# Connectance
calc.C <- function(fw){
  if(length(unique(c(rownames(fw),colnames(fw)))) == 0){
    C = NA
  } else C <- length(fw[fw > 0]) / (length(unique(c(rownames(fw),colnames(fw))))^2)
  return(C)
}

#E.dens
calc.edge.dens <- function(fw) {
  ed <- edge_density(graph_from_adjacency_matrix(as.matrix(fw)))
  return(ed)
}

# Size
calc.S <- function(fw){
  if(length(unique(c(rownames(fw),colnames(fw)))) == 0){ S<-NA 
  }else{
  S <- vcount(graph_from_adjacency_matrix(as.matrix(fw)))
  }
  return(S)
}


# Isolated species

# Degree (in vs out)
calc.mean.degree <- function(fw){
  degree <- length(fw[fw > 0])/ length(unique(c(rownames(fw),colnames(fw))))
  return(degree)
}
calc.norm.mean.degree <- function(fw) {
  if(length(unique(c(rownames(fw),colnames(fw)))) == 0){ dg<-NA}else{
  dg <- degree(graph_from_adjacency_matrix(as.matrix(fw)), normalized = TRUE)
  }
  return(mean(dg))
}

# In degree
calc.in.degree <- function(fw){
  if(length(unique(c(rownames(fw),colnames(fw)))) == 0){ in_degree<-NA}else{
  in_degree <- mean(degree(graph_from_adjacency_matrix(as.matrix(fw),mode=c("directed")),mode=c("in")))
  }
  return(in_degree)
}

# Out degree
calc.out.degree <- function(fw){
  if(length(unique(c(rownames(fw),colnames(fw)))) == 0){ out_degree<-NA}else{
  out_degree <- mean(degree(graph_from_adjacency_matrix(as.matrix(fw),mode=c("directed")),mode=c("out")))
  }
  return(out_degree)
}



# Degree distribution
get.degree.dist <- function(fw) {
  dd <- degree_distribution(graph_from_adjacency_matrix(as.matrix(fw)))
  return(dd)
}

# Transitivity/Clustering Coefficient
calc.cc <- function(fw) {
  if(length(unique(c(rownames(fw),colnames(fw)))) == 0){ cc <-NA}else{
  cc <- transitivity(graph_from_adjacency_matrix(as.matrix(fw)), "global")
  }
  return(cc)
}

# Modularity
calc.mod <- function(fw) {
  mod <- modularity(graph_from_adjacency_matrix(as.matrix(fw)), "global")
  return(mod)
}

# Centrality 
calc.mean.close <- function(fw) {
  close <- closeness(graph_from_adjacency_matrix(as.matrix(fw)))
  return(mean(close))
}
calc.mean.between <- function(fw) {
  btwn <- betweenness(graph_from_adjacency_matrix(as.matrix(fw)),normalized=TRUE)
  return(mean(btwn))
}
calc.mean.eigen <- function(fw) {
  eigen <- eigen_centrality(graph_from_adjacency_matrix(as.matrix(fw)))
  return(mean(eigen$vector))
}


#Node-based centrality measures
calc.eigen <- function(fw) {
  eigen <- data.frame(id = c(1:26), 
                      eigen = eigen_centrality(graph_from_adjacency_matrix(as.matrix(fw)))$vector)
  return(eigen)
}

calc.between <- function(fw) {
  btwn <- data.frame(id = c(1:26), 
                      btwn = betweenness(graph_from_adjacency_matrix(as.matrix(fw))))
  return(btwn)
}

#Node-based degree
calc.degree <- function(fw) {
  dg <- data.frame(id = (rownames(fw)), 
                   degree = degree(graph_from_adjacency_matrix(as.matrix(fw))))
  return(dg)
}

#Node-based clustering coefficient
calc.node.cc <- function(fw) {
  cc <- data.frame(id = (rownames(fw)), 
                   cc = transitivity(graph_from_adjacency_matrix(as.matrix(fw)), 
                                     "local"))
  return(cc)
  
}

calc.mean.path.length <- function(fw) {
  if(length(unique(c(rownames(fw),colnames(fw)))) == 0){ avpath<-NA}else{
  avpath <- igraph::average.path.length(graph_from_adjacency_matrix(as.matrix(fw)))
  }
  return(avpath)
}

# network trophic level (see pg 116 Roopnarine and Dineen 2018)
calc.nw.troph.level <- function(fw){
  # get all the primary producers (they don't eat anything)
  pp <- colnames(fw)[which(colSums(fw) == 0)]
  # turn the fw into a graph object
  gFW <- graph_from_adjacency_matrix(as.matrix(fw))
  # get shortest paths of all species to all species
  sp <- shortest.paths(gFW, v = V(gFW), to = V(gFW))
  
  # create a vector to hold network trophic levels and set primary producers to 1
  ntl <- vector(mode = "numeric", length = nrow(fw))
  names(ntl) <- colnames(fw)
  ntl[pp] <- 1
  
  # run though all the species that aren't primary producers
  for(i in colnames(fw)[! colnames(fw) %in% pp]){
    # figure out who their prey are
    prey <- colnames(fw)[which(fw[,i] != 0)]
    
    sum <- 0
    # and get the mean shortest path of all their prey to a primary producer
    for(j in prey){
      sum <- min(sp[j,pp])
    }
    # calculate network trophic level
    ntl[i] <- 2 + 1/length(prey) * sum
  }
  return(ntl)
}

calc.nw.troph.level.mean<-function(fw) {
  ntl.dat<-calc.nw.troph.level(fw)
  mean(ntl.dat,na.rm=TRUE)
}

#Trophic level
calc.trophlev <- function(fw) {
  if(length(fw)>1){
    tl <- NetIndices::TrophInd(as.matrix(fw))
    trophlev <-mean(tl$TL,na.rm=TRUE)
  }else{
    trophlev<-0
  }
  return(trophlev)
}

#Omnivory index level
calc.oindex <- function(fw) {
  if(length(fw)>1){
    tl <- NetIndices::TrophInd(as.matrix(fw))
    oindex <-mean(tl$OI,na.rm=TRUE)
  }else{
    oindex<-0
  }
  return(oindex)
}

#Calculate trophic species
#Overlap is the degree by which resources and consumers must overlap for X taxa to be considered trophic species
#The "perfect" definition of a trophic species (where X taxa have identical consumers and resources) is an overlap of 1
#An overlap of 0.1 just means that they share a very small number of consumers/resources
#Overlap must be greater than 0
calc.troph.species <- function(fw, overlap=0.000000001) {
  
  if(overlap==0){
    print("Need to use 0.000000001")
  }else{
  }
  
  #Calculate common consumers
  df<-as.data.frame(fw)
  
  if(sum(fw)>0){
    df2<-as.data.frame(df)
    #List consumers of all taxa
    df2$true_cols <- apply(df2, 1, function(data)
      names(which(data == 1)))
    l<-  as.list(df2$true_cols)
    nms <- combn( names(l) , 2 , FUN = paste0 , collapse = "___" , simplify = FALSE )
    # Make the combinations of list elements
    ll <- combn( l , 2 , simplify = FALSE )
    # Intersect the list elements
    out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
    # Output with names
    com_cons<-t(as.data.frame(setNames( out , nms )))
    com_cons <- data.frame(do.call('rbind', strsplit(as.character(rownames(com_cons)),'___',fixed=TRUE)),overlap=com_cons[,1],row.names = NULL)
    com_cons<-unique(com_cons[,])
    #Combine with number of taxa
    com_cons_totals<-as.data.frame(rowSums(df),group=0)
    com_cons<-merge(com_cons,com_cons_totals,by.x="X1",by.y=0,all.x=TRUE)
    com_cons<-merge(com_cons,com_cons_totals,by.x="X2",by.y=0,all.x=TRUE)
    
    #Common resources
    df<-t(as.data.frame(fw))
    df2<-as.data.frame(df)
    df2$true_cols <- apply(df2, 1, function(data)
      names(which(data == 1)))
    l<-  as.list(df2$true_cols)
    nms <- combn( names(l) , 2 , FUN = paste0 , collapse = "___" , simplify = FALSE )
    # Make the combinations of list elements
    ll <- combn( l , 2 , simplify = FALSE )
    # Intersect the list elements
    out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
    # Output with names
    com_res<-t(as.data.frame(setNames( out , nms )))
    com_res <- data.frame(do.call('rbind', strsplit(as.character(rownames(com_res)),'___',fixed=TRUE)),overlap=com_res[,1],row.names = NULL)
    com_res<-unique(com_res[,])
    #com_cons<-tidyr::spread(foo, X1, X2)
    com_res_totals<-as.data.frame(rowSums(df),group=0)
    com_res<-merge(com_res,com_res_totals,by.x="X1",by.y=0,all.x=TRUE)
    com_res<-merge(com_res,com_res_totals,by.x="X2",by.y=0,all.x=TRUE)
    
    common_list<-as.data.frame(rbind(cbind(com_cons,taxatype="common_consumers"),cbind(com_res,taxatype="common_resources")))
    colnames(common_list)<-c("X2","X1","overlap","total_X1","total_X2","taxatype")
    common_list$overlap_ratio<-common_list$overlap/(apply(common_list[,c("total_X1","total_X2")],1,min))
    common_list$overlap_ratio[is.nan(common_list$overlap_ratio)] <- 0
    
    common_list_reshaped<-reshape::cast(as.data.frame(common_list[,c("X2","X1","taxatype","overlap_ratio")]),X1+X2~taxatype,value="overlap_ratio")
    common_list_reshaped$overlapped_combined<-(as.numeric(as.character(common_list_reshaped$common_consumers))*as.numeric(as.character(common_list_reshaped$common_resources)))
    
    graphed_common_list<-subset(common_list_reshaped,overlapped_combined>=overlap)
    
    #In order to define trophic species, cluster by modularity
    #Generate network
    gg<-graph_from_edgelist(as.matrix((graphed_common_list[,c("X1","X2")])),directed=F)
    
    #Add isolates to edgelist (present, but not connected to anything) - important as we are looking at effects of node removal, which would increase isolates
    all_names<-unique(rownames(df),colnames(df))
    gg_names<-as.character(unique(graphed_common_list[,c("X1")]))
    gg_names<-append(gg_names,as.character(unique(graphed_common_list[,c("X2")])))
    gg_names<-unique(gg_names)
    missing_names<-all_names[!(all_names %in% gg_names)]
    
    
    #With isolates
    graphed_common_list2<-rbind(graphed_common_list[,c("X1","X2","overlapped_combined")],cbind(X1=missing_names,X2=missing_names,overlapped_combined=""))
    graphed_common_list2$overlapped_combined <- sub("^$", 0, graphed_common_list2$overlapped_combined)
    graphed_common_list2$overlapped_combined<-as.numeric(as.character(graphed_common_list2$overlapped_combined))
    gg_iso<-graph_from_edgelist(as.matrix(graphed_common_list2[,c("X1","X2")]),directed=F)
    #If two taxa don't exceed the overlap threshold, then their overlap is set to 0 (i.e. no trophic species relationship)
    graphed_common_list2[,c("overlapped_combined")]<-ifelse(graphed_common_list2[,c("overlapped_combined")]>=overlap,graphed_common_list2[,c("overlapped_combined")],0)
    E(gg_iso)$weight<-as.numeric(as.character(graphed_common_list2[,c("overlapped_combined")]))
    gg_iso<-igraph::simplify(gg_iso, remove.loops=TRUE)
    lou<-cluster_louvain(gg_iso,weights=E(gg_iso)$weight)
    
    #Generate list of trophic species
    troph_spec_list<-as.data.frame(cbind(V(gg_iso)$name,lou$membership))
    print(paste("Number of trophic species: ",length(unique(troph_spec_list$V2)),"overlap level: ", overlap, "total number of taxa: ",nrow(fw)))
    row<-c(NumberOfTrophicSpecies=length(unique(troph_spec_list$V2)),OverlapLevel=overlap,NumberOfTotalSpecies=nrow(fw))
    return(row)
  }else{
    #If there are species in the web
    if(sum(fw)>0){
      troph_spec_list<-as.data.frame(cbind(rownames(fw),1:nrow(fw)))
      print(paste("Number of trophic species: ",length(unique(troph_spec_list$V2)),"overlap level: ", overlap, "total number of taxa: ",nrow(fw)))
      row<-c(NumberOfTrophicSpecies=length(unique(troph_spec_list$V2)),OverlapLevel=overlap,NumberOfTotalSpecies=nrow(fw))
      return(row)
    }else{
      #If there are no species left
      print(paste("Number of trophic species: ",0,"overlap level: ", overlap, "total number of taxa: ",0))
      return(c(0,overlap,0))
    }
  }
}



#Compare numbers of HardSoft and BenthicPelagic
calc.att.numbers <- function(fw.node.list,node.attributes,runs){
  
  att_fw.list<-merge(node.attributes,fw.node.list,by.x="Species",by.y="species",all.x=TRUE)
  att_fw.list$run0<-1
  
  numbas_list<-c()
  for(row in 0:(ncol(fw.node.list)-1)){
    col_run<-paste("run",row,sep="")
    att_sub<-att_fw.list[,c(col_run,"HardSoft","BenthicPelagic")]
    att_sub$selected_run<-att_sub[,c(col_run)]
    att_sub<-subset(att_sub, selected_run==1)
    
    att_sub$HardSoft<-factor(att_sub$HardSoft, levels = 1:3)
    att_sub$BenthicPelagic<-factor(att_sub$BenthicPelagic, levels = 1:3)
    
    if(nrow(att_sub)>0){
      numbas<-as.matrix(rbind(cbind(as.data.frame(table(att_sub$HardSoft)),att="HardSoft"),
                              cbind(as.data.frame(table(att_sub$BenthicPelagic)),att="BenthicPelagic")))
      numbas_list<-rbind(numbas_list,cbind(numbas,run=row))
      
    }else{
      numbas<-as.matrix(rbind(c(1,0,"HardSoft"),
                              c(2,0,"HardSoft"),
                              c(3,0,"HardSoft"),
                              c(1,0,"BenthicPelagic"),
                              c(2,0,"BenthicPelagic"),
                              c(3,0,"BenthicPelagic")))
      numbas_list<-rbind(numbas_list,cbind(numbas,run=row))
      
    }
    
    
  }
  return(numbas_list)
}




#Compare in and out degree distributions, calculate AIC & BIC
degree_distribution_in <- function(fw, name){
  AIC <- Cumulative <- Exp <- fit <- model <- truncated <- NULL
  totaldegree<- degree(fw, mode="in")
  K <- 0:max(totaldegree)
  For.Graph<- data.frame(K = K, Cumulative = NA, Scenario = name)
  for(i in 1:length(K)){
    For.Graph$Cumulative[i] <- sum(totaldegree>K[i])/length(totaldegree)
  }
  
  exp.model <- nls(Cumulative~exp(-K/y),start= list(y=0.1), data = For.Graph)
  
  For.Graph$Exp <- predict(exp.model)
  Summs.exp <- glance(exp.model)
  Summs.exp$model <- "Exponential"
  
  power <- filter(For.Graph, K != 0 & Cumulative != 0)
  powerlaw.model <- nls(Cumulative~K^-y, start= list(y=0), data = power)
  power$power <- predict(powerlaw.model)
  For.Graph <- full_join(For.Graph, power)
  Summs.power <- glance(powerlaw.model)
  Summs.power$model <- "Power"
  
  truncated.powerlaw.model <- nls(Cumulative~(K^-y)*(exp(-K/y)), start = list(y=1), data = power)
  power$truncated <- predict(truncated.powerlaw.model)
  For.Graph <- full_join(For.Graph, power)
  Summs.truncated <- glance(truncated.powerlaw.model)
  Summs.truncated$model <- "truncated"
  
  Summs <- full_join(Summs.exp, Summs.power)
  Summs <- full_join(Summs, Summs.truncated)
  Summs <- arrange(Summs, AIC)
  
  DF2 <- For.Graph %>% filter(K != 0 & Cumulative != 0) %>% gather(key = model, value = fit, Exp, power, truncated)
  
  g <- ggplot(DF2, aes_string(x = "K", y = "Cumulative")) + geom_line() + geom_point()+ scale_x_log10() + scale_y_log10() + theme_classic() + geom_line(aes_string(y ="fit", color = "model"))
  
  g
  
  return(list(DDvalues = For.Graph, models = Summs, graph = g))
}



degree_distribution_out <- function(fw, name){
  AIC <- Cumulative <- Exp <- fit <- model <- truncated <- NULL
  totaldegree<- degree(fw, mode="out")
  K <- 0:max(totaldegree)
  For.Graph<- data.frame(K = K, Cumulative = NA, Scenario = name)
  for(i in 1:length(K)){
    For.Graph$Cumulative[i] <- sum(totaldegree>K[i])/length(totaldegree)
  }
  
  exp.model <- nls(Cumulative~exp(-K/y),start= list(y=0.1), data = For.Graph)
  
  For.Graph$Exp <- predict(exp.model)
  Summs.exp <- glance(exp.model)
  Summs.exp$model <- "Exponential"
  
  power <- filter(For.Graph, K != 0 & Cumulative != 0)
  powerlaw.model <- nls(Cumulative~K^-y, start= list(y=0), data = power)
  power$power <- predict(powerlaw.model)
  For.Graph <- full_join(For.Graph, power)
  Summs.power <- glance(powerlaw.model)
  Summs.power$model <- "Power"
  
  truncated.powerlaw.model <- nls(Cumulative~(K^-y)*(exp(-K/y)), start = list(y=1), data = power)
  power$truncated <- predict(truncated.powerlaw.model)
  For.Graph <- full_join(For.Graph, power)
  Summs.truncated <- glance(truncated.powerlaw.model)
  Summs.truncated$model <- "truncated"
  
  Summs <- full_join(Summs.exp, Summs.power)
  Summs <- full_join(Summs, Summs.truncated)
  Summs <- arrange(Summs, AIC)
  
  DF2 <- For.Graph %>% filter(K != 0 & Cumulative != 0) %>% gather(key = model, value = fit, Exp, power, truncated)
  
  g <- ggplot(DF2, aes_string(x = "K", y = "Cumulative")) + geom_line() + geom_point()+ scale_x_log10() + scale_y_log10() + theme_classic() + geom_line(aes_string(y ="fit", color = "model"))
  
  g
  
  return(list(DDvalues = For.Graph, models = Summs, graph = g))
}



##calculate prey-averaged trophic level for each node in the foodweb
calc.prey.avg.TL = function(fw) {
  network = graph_from_adjacency_matrix(as.matrix(fw), mode = "directed")
  preynum = degree(network, mode = "in")
  tl = NetIndices::TrophInd(as.matrix(fw))
  pa.tl = vector()
  for(j in 1:ncol(fw)) {
    name = colnames(fw)[j]
    tlj = 0
    prey = rownames(fw)[which(fw[,j] != 0)]
    if(length(prey) != 0) {
    for(i in 1:length(prey)) {
      tli = tl[prey[i],]$TL
      nj = preynum[[name]]
      tlj = tlj + (1 + (tli/nj))
    }
      pa.tl = c(pa.tl, tlj)
    }else {pa.tl = c(pa.tl, NA)}
    
  }
  names(pa.tl) = colnames(fw)
  return(pa.tl)
}


##Calculate short weighted trophic level - via Susanne Kortsch
calc.swtl<-function(web){
  
  if(length(web)>1){
    
    library(NetIndices)
    library(igraph)  
    
    web<-as.matrix(web)
    rownames(web)<-colnames(web)
    
    basal <- rownames(subset(web, apply(web, 2, sum)==0) & apply(web, 1, sum)!= 0)
    edge.list_web <- graph.adjacency(web, mode = "directed");
    paths_prey <- shortest.paths(graph = edge.list_web, v= V(edge.list_web),
                                 to = V(edge.list_web)[basal], mode = "in", weights = NULL, algorithm = "unweighted")
    paths_prey[is.infinite(paths_prey)] <- NA
    shortest_paths <- suppressWarnings(as.matrix(apply(paths_prey, 1, min, na.rm=TRUE)))
    #longest_paths <- as.matrix(apply(paths_prey, 1, max, na.rm=TRUE))
    in_deg <- apply(web, 2, sum) 				## Take each columns, sum  the rows
    out_deg <- apply(web, 1, sum) 			## Take each rows, sum the columns
    
    
    # Shortest TL
    sTL <- 1 + shortest_paths  # Commonly, detritus have a TL value of 1. (Shortest path to basal = 0)
    # Longest TL
    #lTL <- 1 + longest_paths
    
    S<-dim(web)[1]
    
    # Creating the matrix  
    short_TL_matrix <- matrix(NA, S, S)
    #long_TL_matrix <- matrix(NA, S, S)
    prey_ave <- rep(NA, S) # averaged prey
    chain_ave<-rep(NA, S)
    #
    for(j in 1:S){
      for(i in 1:S){
        lij <- web[i, j] 					# is the interaction matrix
        prey_ave[j] <- 1 / in_deg[j] 			# if Bas species in, no in-degrees ; re-attribute a value of 0
        short_TL_matrix[i,j] <- lij * sTL[i]   	# Shortest path in matrix
        #long_TL_matrix[i,j] <- lij * lTL[i]
      }  
    }
    
    prey_ave[which(prey_ave == Inf)] <- 0
    prey_ave[which(prey_ave == -Inf)] <- 0
    
    short_TL_matrix[is.nan(short_TL_matrix)] <- 0

    
    short_TL_matrix[which(short_TL_matrix == Inf)] <- 0
    short_TL_matrix[which(short_TL_matrix == -Inf)] <- 0
    

    
    sumShortTL <- as.matrix(apply(short_TL_matrix, 2, sum)) # sum all shortest path

    
    # Short-weighted TL
    # weigth by the number of prey
    
    SWTL <- matrix(data = NA, nrow = S, ncol = 1)
    for(i in 1:S){
      SWTL[i] <- 1 + (prey_ave[i] * sumShortTL[i])  
    }
    

    
    
    TL_NI<-TrophInd(web)
    TL_NetInd<- TL_NI[,1]

    TLS<-cbind(SWTL, TL_NetInd)
    colnames(TLS)<-c("SWTL","NetInd_TL")
    rownames(TLS)<-rownames(web)
    
    return(TLS)
  }else{
    return(NA)
  }
  
}


##Calc mean SWTL
calc.mean.swtl<-function(web) {
  if(length(web)>1){
    dat<-calc.swtl(web)
    dat2<-mean(dat[,c("SWTL")],na.rm=TRUE)
      }else{
    return(NA)
  }
}






##Make perfect trophic species web
make.trophic.species.web <- function(fw) {
  if(length(fw)>1){
    
    dupl_cols <- duplicated(fw)
    dupl_rows <- duplicated(t(fw))
    dupl <- which(dupl_rows+dupl_cols == 2)
    
    newnw <- as.data.frame(fw[-dupl,-dupl])
    
  }else{
    newnw<-NULL
  }
  return(newnw)
}


####SWOI----
diet_get <- function (Tij) {
  IntFlowTo <- rowSums(Tij)
  p <- Tij
  for (i in 1:ncol(Tij)){
    p[i, ] <- Tij[i, ]/IntFlowTo[i]
  }
  p[is.na(p)] <- 0
  return(p)
}

calc.swoi<-function(web){
  
  if(length(web)>1){
  Tij <- t(web)
  p <- diet_get(Tij)
  
  
  ncomp <- ncol(Tij)
  A <- -p
  diag(A) <- 1
  B <- rep(1, ncomp)
  
  
  TL<-calc.swtl(web)[,1]
  
  OI <- vector(length = nrow(web))
  for (i in 1:ncomp) OI[i] <- sum((TL - (TL[i] - 1))^2 * p[i, ])
  
  return(OI)
  }else{
    return(NA)
  }
  
}




####SW SOI----
calc.SOI<-function (Flow = NULL, Tij = t(Flow), Import = NULL, Export = NULL, 
                    Dead = NULL) 
{
  if(length(Flow)>1){
    OI<-calc.swoi(Flow)
    ncomp <- ncol(Flow)
    Q <- vector(length = ncomp)
    OlogQ<-vector(length = ncomp)
    
    for (i in 1:ncomp) Q[i] <- sum(Flow[,i])
    
    for (i in 1:ncomp) OlogQ[i] <-OI[i]*log(Q[i])
    ##Don't know if this is ok, to replace NaN with 0
    OlogQ[is.nan(OlogQ)] = 0
    
    SOI<-sum(OlogQ)/log(sum(Q))
  }else{
    SOI<-NA
  }
  return(SOI)
  
}

calc.troph.species.LIST <- function(fw, overlap=0.000000001) {
  
  if(overlap==0){
    print("Need to use 0.000000001")
  }else{
  }
  
  #Calculate common consumers
  df<-as.data.frame(fw)
  
  if(sum(fw)>0){
    df2<-as.data.frame(df)
    #List consumers of all taxa
    df2$true_cols <- apply(df2, 1, function(data)
      names(which(data == 1)))
    l<-  as.list(df2$true_cols)
    nms <- combn( names(l) , 2 , FUN = paste0 , collapse = "___" , simplify = FALSE )
    # Make the combinations of list elements
    ll <- combn( l , 2 , simplify = FALSE )
    # Intersect the list elements
    out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
    # Output with names
    com_cons<-t(as.data.frame(setNames( out , nms )))
    com_cons <- data.frame(do.call('rbind', strsplit(as.character(rownames(com_cons)),'___',fixed=TRUE)),overlap=com_cons[,1],row.names = NULL)
    com_cons<-unique(com_cons[,])
    #Combine with number of taxa
    com_cons_totals<-as.data.frame(rowSums(df),group=0)
    com_cons<-merge(com_cons,com_cons_totals,by.x="X1",by.y=0,all.x=TRUE)
    com_cons<-merge(com_cons,com_cons_totals,by.x="X2",by.y=0,all.x=TRUE)
    
    #Common resources
    df<-t(as.data.frame(fw))
    df2<-as.data.frame(df)
    df2$true_cols <- apply(df2, 1, function(data)
      names(which(data == 1)))
    l<-  as.list(df2$true_cols)
    nms <- combn( names(l) , 2 , FUN = paste0 , collapse = "___" , simplify = FALSE )
    # Make the combinations of list elements
    ll <- combn( l , 2 , simplify = FALSE )
    # Intersect the list elements
    out <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
    # Output with names
    com_res<-t(as.data.frame(setNames( out , nms )))
    com_res <- data.frame(do.call('rbind', strsplit(as.character(rownames(com_res)),'___',fixed=TRUE)),overlap=com_res[,1],row.names = NULL)
    com_res<-unique(com_res[,])
    #com_cons<-tidyr::spread(foo, X1, X2)
    com_res_totals<-as.data.frame(rowSums(df),group=0)
    com_res<-merge(com_res,com_res_totals,by.x="X1",by.y=0,all.x=TRUE)
    com_res<-merge(com_res,com_res_totals,by.x="X2",by.y=0,all.x=TRUE)
    
    common_list<-as.data.frame(rbind(cbind(com_cons,taxatype="common_consumers"),cbind(com_res,taxatype="common_resources")))
    colnames(common_list)<-c("X2","X1","overlap","total_X1","total_X2","taxatype")
    common_list$overlap_ratio<-common_list$overlap/(apply(common_list[,c("total_X1","total_X2")],1,min))
    common_list$overlap_ratio[is.nan(common_list$overlap_ratio)] <- 0
    
    common_list_reshaped<-reshape::cast(as.data.frame(common_list[,c("X2","X1","taxatype","overlap_ratio")]),X1+X2~taxatype,value="overlap_ratio")
    common_list_reshaped$overlapped_combined<-(as.numeric(as.character(common_list_reshaped$common_consumers))*as.numeric(as.character(common_list_reshaped$common_resources)))
    
    return(common_list_reshaped)
    
  }else{
    #If there are species in the web
    if(sum(fw)>0){
      troph_spec_list<-as.data.frame(cbind(rownames(fw),1:nrow(fw)))
      print(paste("Number of trophic species: ",length(unique(troph_spec_list$V2)),"overlap level: ", overlap, "total number of taxa: ",nrow(fw)))
      row<-c(NumberOfTrophicSpecies=length(unique(troph_spec_list$V2)),OverlapLevel=overlap,NumberOfTotalSpecies=nrow(fw))
      return(row)
    }else{
      #If there are no species left
      print(paste("Number of trophic species: ",0,"overlap level: ", overlap, "total number of taxa: ",0))
      return(c(0,overlap,0))
    }
  }
}






####FUNCTIONS REQUIRED FOR SYSTEM OMNIVORY INDEX


##==============================================================================
## Implementation of Network indices, as in
## Latham, LG II 2006 Ecol. Modeling 192: 586-600
##
## Implemented: Julius Kones      - University Nairobi
##              Karline Soetaert  - Netherlands Institute of Ecology
##
## Two local functions:
##     InternalNetwork
##     Diet
##==============================================================================

##==============================================================================
InternalNetwork <- function (Tij  ,          # to-from
                             Import,         # flow from external (colNr Tij)
                             Export)         # flow to external (colNr Tij)
  
{                      
  
  ##------------------------------------------------------------------------
  ## Tij[i,j] is a matrix with Tij[i,j]  flow from j to i
  ## note: component position in rows and columns must be the same - not checked
  ##------------------------------------------------------------------------
  
  if (is.character(Import))
    import <- which(colnames(Tij)%in%Import) else
      import <- Import
    if (length(import) != length(Import))
      stop("Import not recognized")
    if (is.character(Export))
      export <- which(rownames(Tij)%in%Export) else
        export <- Export
      if (length(import) != length(Import))
        stop("Import not recognized")
      if (length(export) != length(Export))
        stop("Export not recognized")
      
      ##
      ## CHECK THE INPUT
      ##
      
      # Flow or Tij should be inputted
      if (is.null(Tij))
        stop ("cannot calculate indices - Flow or Tij should be inputted")
      
      #_________________________________________________________________________________
      # NUMBER OF COMPARTMENTS, export, import, internal flows,..
      #_________________________________________________________________________________
      
      # Size of the matrices; without the externals, the matrix has to be square
      
      ncomp     <- ncol(Tij)-length(import)
      if (ncomp != nrow(Tij)-length(export))
        stop ("cannot calculate indices - internal flow input matrix not square ")
      
      #_________________________________________________________________________________
      # ARRAYS DECLARATION         
      #_________________________________________________________________________________
      
      # indices to elements of T that are internal 
      iN  <- setdiff(1:nrow(Tij),export)  # internal rows    of Tij
      jN  <- setdiff(1:ncol(Tij),import)  # internal columns of Tij
      
      
      # Total internal flows, externals removed.
      Tint        <- Tij
      if (! is.null(export))
        Tint <- Tint[-export,]
      if (! is.null(import))
        Tint <- Tint[,-import]
      
      # Total flows, including flow to/from externals
      FlowFrom   <- colSums(Tij)
      FlowTo     <- rowSums(Tij)
      FlowFromC  <- FlowFrom[jN]    # just the total from internal compartments
      FlowToC    <- FlowTo  [iN]
      
      return(list(Tint=Tint,iN=iN,jN=jN,
                  import=import,export=export,
                  FlowFrom=FlowFrom,
                  FlowTo  = FlowTo,
                  FlowFromC=FlowFromC,
                  FlowToC  =FlowToC))
      
} # END InternalNetwork

##==============================================================================
##
## Internal function: estimates the diet composition
##
##==============================================================================
Diet <- function (Tint,                   # Calculates diet composition
                  Dead=NULL,              # index from Dead to Tint
                  iN=1:nrow(Tint))      
  
{
  
  ## p matrix contains the diet composition of predator i
  IntFlowTo   <- rowSums(Tint)  # Total food ingested
  
  p           <- Tint
  
  for (i in 1:ncol(Tint))
    p[i,] <- Tint[i,]/IntFlowTo[i]
  
  p[is.na(p)] <- 0
  
  ## take into account dead matter; Dead refers to column/row in Tij
  ## N$iN maps from Tint to Tij
  
  if (! is.null(Dead))
    p[which(iN%in% Dead),]<-0
  return(p)
  
} # END Diet   
