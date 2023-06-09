#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################

## The following function generates and plots networks and calculate edge weights based on sampling effort (minutes observed). 
## Edge weights are based on the duration of contacts per minutes of observation. 
## Observation times are calculated using the scan data and durations of contacts are based on focal data.

####looking at the network on a population level, contact duration
net.foc.population <- function(date.ref=as.Date("6/1/2008", "%m/%d/%Y"), duration=90, after=T, beh="Approach_<_1M", edgewidth=100000, direction="directed",arrowsize=.5, labelcex=1, labeldist=.5, labels=T, margin=2, lay=99)
  #### inputs: reference/start date, amount of time from the reference date (duration), before or after start date, behavior type, edge width, directed or undirected graph
  
{ #beginning of net.foc function
  
  ###first, we want to only look at focal individuals and behavioral data for the indicated period
  ##list of individuals for this group and period	
  if (after) #if want to look at time period after reference date (after=T)
  {
    foc1 <- dat2[dat2$Date > date.ref & dat2$Date < date.ref+duration,] #focal data after reference date but before duration ends																															
    census1 <- census[census$Date > date.ref & census$Date < date.ref+duration,]  #census data for after reference date but before duration ends
    scan1 <- scan[scan$Date > date.ref & scan$Date < date.ref+duration,]	#scan data after reference daste but before duration ends
  } 
  else #if after=F (want to look at time period before )
  {
    foc1 <- dat2[dat2$Date < date.ref & dat2$Date > date.ref - duration,] #focal data before reference date but greater than reference date minus duration
    census1 <- census[census$Date < date.ref & census$Date > date.ref - duration,] #census data before reference date but greater than reference date minus duration
    scan1 <- scan[scan$Date < date.ref & scan$Date > date.ref - duration,]# #scan data before reference date but greater than reference date minus duration
    
  }
  
  foc1$Focal <- droplevels(foc1$Focal)
  foc1$Initiator <- droplevels(foc1$Initiator)
  foc1$Receiver <- droplevels(foc1$Receiver)
  scan1$Focal <- droplevels(scan1$Focal)
  
  scan.indiv <- levels(scan1$Focal)
  foc.indiv <- union(levels(foc1$Initiator),levels(foc1$Receiver))
  lsind2 <- sort(as.character(na.omit(unique(c(scan.indiv,foc.indiv)))))
  
  foc2 <- foc1[foc1$Focal%in%lsind2 & foc1$Initiator%in%lsind2 & foc1$Receiver%in%lsind2, ] #create new focal set that includes focal individuals, initators, and receivers in "lsind" vector
  
  scan2 <- scan1[scan1$Focal%in%lsind2, ] ##create new scan dataset that includes focal individuals in lsind
  
  ##create an adjacency matrix and social network for the behavioral data 
  if(nrow(foc2)==0 | nrow(scan2)==0) return(NA) #if there is no scan data or no focal data for the indicated time period, return NA
  else {
    ##observation times
    scan2$Focal  <-  factor(scan2$Focal, levels=lsind2) #factor() returns an object of class "factor" that has a set of integer codes the length of x with a "levels" attribute of mode 				character and unique entries 
    hours <- apply(table(scan2$Focal, factor(scan2$Date):factor(scan2$ScanTime)), 1, function(v) sum(v>0))*10
    
    #calculate the number of minutes each dyad is observed by building a contingency table with the times scan focal individuals  were observed
    ##unique scan timepoints for each focal individual = each unit corresponds to one 10 minute increment 
    ##multipy by 10 to get number of minutes observed    
    foc2 <- foc2[foc2$Behavior%in%beh,] #foc2 contains only behavior type indicated by user in the function 

    foc2$Start <- strptime(foc2$Start, "%H:%M:%S") #reformat time; updated focal data txt file (up to 2012)
    foc2$Finish <- strptime(foc2$Finish, "%H:%M:%S") #updated focal data txt file (up to 2012)
    contact.duration <- as.numeric(difftime(foc2$Finish,foc2$Start, units="secs"))/60 #calculate contact duration by difference b/w start and end times (number of minutes)
    foc2$Initiator <- factor(foc2$Initiator, levels=lsind2)
    foc2$Receiver <- factor(foc2$Receiver, levels=lsind2)
    m <- matrix(tapply(contact.duration, foc2$Initiator:foc2$Receiver, sum), length(lsind2), length(lsind2)) #sum of contact durations for both receivers and initiators
    m[is.na(m)] <- 0 #turn NAs into 0s
    diag(m) <- 0 #no self loops
    temp <- matrix(hours, length(lsind2), length(lsind2))
    weight <- temp + t(temp) #create a matrix of hours each focal individual is observed (sampling effort includes amount of time each individual within a dyad was observed)
    m2 <- m/(weight)
    m2[is.na(m2)] <- 0 #m2 is weighted adjacency matrix   
    m2[m2==Inf] <- 0
    
    g <- graph.adjacency(m2, mode=direction,weighted=TRUE, diag=FALSE) #graph the adjacency matrix, this becomes a network object "g" 
    edge_weights <- E(g)$weight
    inv_edge_weights <- (1/(E(g)$weight)) #these will be used for calculating path length 
    V(g)$name <- lsind2
    V(g)$obs_hours <- as.vector(hours)
    edge.w <- t(m2)[t(m2)>0] #show edges for interactions with weights > 0
    hours <- as.vector(hours) #amount of time observed
    g <- set.vertex.attribute(graph=g, name="obs_hours", index=V(g), value=hours)
    g <- delete.vertices(g,which(degree(g)<1)) #get rid of individuals with degrees less than 1 (ie, no contacts)
    sifaka.edgelist <- cbind(get.edgelist(g, names=TRUE), E(g)$weight)#list of all the edges of the network 
    g_adj <- get.adjacency(g, attr="weight", sparse=FALSE)
    g_adj2 <- get.adjacency(g, edges=TRUE, attr="weight", sparse=FALSE)
    centrality <- apply(g_adj+t(g_adj),1,sum) #weighted degree = sum of all edge weights
    ebc <- edge.betweenness.community(g)
    
    vcol <- ifelse(sex.dat[match(V(g)$name, sex.dat[,1]),2]=="M", colors()[430], colors()[367]) #separate colors for males and females
    V(g)$sex <- as.character(sex.dat[match(V(g)$name, sex.dat[,1]),2])
    V(g)$size <- 10
    E(g)$color <- "black"
    nodesize <- 10
    if(class(lay)!="function") lay1 <- layout.fruchterman.reingold(g) else lay1 <- lay 
    
    if(labels) labs <- V(g)$name else labs <- "" #use individuals in lsind to label nodes
    
    if (after) #if want to look at duration after reference date
    {
      begin.date <- date.ref
      end.date <- date.ref+duration
    } 
    else #if want to look at duration before reference date
    {
      end.date <- date.ref
      begin.date <- date.ref-duration
    }
    plot(g, layout=lay1, vertex.color=vcol, edge.width=edge.w*edgewidth, edge.arrow.size=arrowsize, vertex.label=labs, vertex.label.dist=labeldist, vertex.label.cex=labelcex, main=paste(format(begin.date,"%b %d, %Y"), "-", format(end.date, "%b %d, %Y")), margin=c(0,0,.3,0)) #plot the network object        
   
    #####NETWORK STATISTICS
    sifaka.density <- graph.density(g, loops=FALSE) #edge density of entire network
    network.size <- vcount(g)
    female.number <- length(V(g)[V(g)[sex=="F"]]) #number of females in network
    male.number <- length(V(g)[V(g)[sex=="M"]]) #number of males in network
    paths <- shortest.paths(g, v=V(g), to=V(g), mode="all", weights=inv_edge_weights)
    
    ###Degree
    sifaka.degree <- degree(g, mode="all", loops=FALSE) #list of individual degrees
    degree.distribution <- degree.distribution(g, cumulative=FALSE) #degree distribution
    
    ###Degree corrected for sampling effort
    sifaka.degree.corrected <- sifaka.degree/V(g)$obs_hours
    sifaka.degree.corrected[is.infinite(sifaka.degree.corrected)] <- 0
    sifaka.degree.corrected[is.na(sifaka.degree.corrected)] <- 0
    
    ### In-Degree
    sifaka.in.degree <- degree(g, loops=FALSE, mode="in") #in degree 

    ###In-Degree corrected for sampling effort
    sifaka.in.degree.corrected <- sifaka.in.degree/V(g)$obs_hours
    sifaka.in.degree.corrected[is.infinite(sifaka.in.degree.corrected)] <- 0
    sifaka.in.degree.corrected[is.na(sifaka.in.degree.corrected)] <- 0
    
    ##Out-degree
    sifaka.out.degree <- degree(g, loops=FALSE, mode ="out")
    
    ###Out-Degree corrected for sampling effort
    sifaka.out.degree.corrected <- sifaka.out.degree/V(g)$obs_hours
    sifaka.out.degree.corrected[is.infinite(sifaka.out.degree.corrected)] <- 0
    sifaka.out.degree.corrected[is.na(sifaka.out.degree.corrected)] <- 0

    ###Centrality metrics
    
    ##closeness
    sifaka.closeness <- closeness(g, mode="all", weights=edge_weights) #closeness centrality
    sifaka.in.closeness <- closeness(g, mode="in")
    sifaka.out.closeness <- closeness(g, mode="out")
    
    ##betweenness 
    sifaka.betweenness <- betweenness(g, directed=FALSE, weights=edge_weights) #betweenness centrality
    
    ###Weighted degree 
    sifaka.strength <- graph.strength(g,vids=V(g), mode="all", loops=FALSE, weights=edge_weights) #combined weighted degree in and out degree
    
    sifaka.strength.in <- graph.strength(g,vids=V(g), mode="in", loops=FALSE, weights=edge_weights) #weighted in degree
    
    sifaka.strength.out <- graph.strength(g,vids=V(g), mode="out", loops=FALSE) #weighted out degree  
    
    return(list(ind=V(g)$name, 
                g = g,
                size=network.size,
                mat=m2, 
                edgelist = sifaka.edgelist,
                ebc=ebc,
                hours = V(g)$obs_hours,
                mat2=g_adj,
                mat3=g_adj2,
                weight_matrix=weight,
                edge_weights = edge_weights,
                sifaka_centrality = centrality, #same as strength
                sifaka_edgelist=sifaka.edgelist, 
                male_number =male.number,
                female_number=female.number,
                
                ####degree
                degree=sifaka.degree, 
                degree_distribution = degree.distribution, 
                
                ####general stats for whole network
                density = sifaka.density, 
                paths = paths,

                ###in degree
                in_degree= sifaka.in.degree, 
                
                ###out degree
                out_degree=sifaka.out.degree, 
                
                ###Degree corrected for sampling effort
                sifaka_degree_corrected = sifaka.degree.corrected, 
                
                ###In-Degree corrected for sampling effort
                sifaka_in_degree_corrected = sifaka.in.degree.corrected, 
                
                ###Out-Degree corrected for sampling effort
                sifaka_out_degree_corrected = sifaka.out.degree.corrected,      
                
                ###betweenness
                sifaka_betweenness = sifaka.betweenness,
                
                ####weighted degree
                weighted_degree= sifaka.strength, #same as centrality
                weighted_in_degree= sifaka.strength.in, 
                weighted_out_degree=sifaka.strength.out, 
                sex1=V(g)$sex,
                sex=sex.dat[match(V(g)$name, sex.dat[,1]),2])) #output the individuals, weighted adjacency matrix, males, and females
  }
} #end of function