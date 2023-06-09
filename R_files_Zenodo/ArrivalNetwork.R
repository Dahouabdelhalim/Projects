#######Arrival Time Networks#####

### Last Updated: 14/12/2022
### Alex Chan
## Added multi-feeder support

# rm(list=ls())

##Required packages
library('dplyr')

BuildNetworkArrive <- function(df,time.col="coded.time",ID.col="birdID",feeder.col=FALSE, interval = 300, t=30){
  #browser()
  ## df: dataframe input
  ## ID.col: Column name of BirdID
  ## time.col: column name of the time, rounded to seconds
  ## interval: threshold for estimating flocks
  ## t: time threshold to define arrival time associations
  if(sum(colnames(df)==time.col)==0){
    stop(paste("Column Name:", time.col, "not found!"))
  }
  
  if(sum(colnames(df)==ID.col)==0){
    stop(paste("Column Name:", ID.Col, "not found!"))
  }
  
  if(feeder.col == FALSE){
    #do nothing
    FeederBool =FALSE
    FeederNum = 1
  }else{
    if(sum(colnames(df)==feeder.col)==0){#there is multiple feeder data
      stop(paste("Column Name:", feeder.col, "not found!"))
    }
    names(df)[which(names(df)==feeder.col)] <- "Feeder"
    FeederBool =TRUE 
    FeederVec = unique(df$Feeder)
    FeederNum = length(FeederVec)
  }
  
  
  #Rename time column based on user specification
  names(df)[which(names(df)==time.col)] <- "coded.time"
  names(df)[which(names(df)==ID.col)] <- "birdID"
  
  #initialize interaction matrix
  ID <- unique(df$birdID)
  InteractionMat <- matrix(data=0, nrow=length(ID), ncol=length(ID),
                           dimnames= list(ID,ID))
  
  ## 1. get arrival times of all birds, defined by set interval of down time
  for(f in 1:FeederNum){
    if(FeederBool==TRUE){
      FeedSubset <- subset(df,df$Feeder==FeederVec[f])
      print(paste("Analyzing location", FeederVec[f]))
    }else{
      FeedSubset <- df
    }
  
    
    
    BirdArrival <- data.frame(matrix(ncol=2, nrow=0))
    names(BirdArrival) <- c("birdID", "ArriveTime")
    
    for(i in 1:length(ID)){
      BirdDF <- filter(FeedSubset, FeedSubset$birdID ==ID[i])
      if(nrow(BirdDF)==0){ ##If id doesnt exist in this subset
        next
      }
      
      BirdDF <- BirdDF %>% arrange(coded.time) %>% mutate(time.diff = as.numeric(coded.time - lag(coded.time))) #calculate time elapsed before next entry
      BirdDF$time.diff[1] <- 0 #replace first NA with 0 to assign flock
      BirdDF$flock <- cumsum(ifelse(BirdDF$time.diff> interval,1,0)) #assign flock id based on time elapsed (600s)
      BirdSum <- BirdDF %>% group_by(flock) %>% 
        summarise(ArriveTime = min(coded.time)) %>%
        mutate(birdID = ID[i])%>%
        select(-flock)
      BirdArrival <- rbind(BirdArrival, BirdSum)
      
    }
    
    ##2. Find arrival times within t seconds
    BirdArrival <- BirdArrival %>% arrange(ArriveTime) %>% mutate(time.diff = as.numeric(ArriveTime - lag(ArriveTime))) #calculate time elapsed before next entry
    BirdArrival$time.diff[1]<- 0 #replace first NA with 0 to assign flock
    BirdArrival$flock <- cumsum(ifelse(BirdArrival$time.diff> t,1,0)) +1 #assign flock id based on time elapsed (600s)
    #Added plus 1 to make sure the first flock is also counted
    
    
    
    for(j in (1:max(BirdArrival$flock)+1)){
      #browser()
      flockSub <- filter(BirdArrival, BirdArrival$flock==(j-1))
      
      if(nrow(flockSub) >1){
        #multiple birds
        BirdIndex <- which(colnames(InteractionMat) %in% flockSub$birdID)
        InteractionMat[BirdIndex,BirdIndex] <- InteractionMat[BirdIndex,BirdIndex]+ 1
        diag(InteractionMat)[BirdIndex] <- diag(InteractionMat)[BirdIndex]-1 ##remove this from diagonal, remove over counting
        
        
      }else{
        #only one bird
        BirdID <- as.character(flockSub$birdID[1])
        InteractionMat[BirdID,BirdID] <- InteractionMat[BirdID, BirdID]+ 1
      }
      
    }
  }
  
  ##3. build SRI edgelist
  ##Get edge list with SRI
  TotalOcc <- apply(InteractionMat, 1, sum) #Total number of occurences each ind showed up
  UpperMat <- InteractionMat*upper.tri(InteractionMat,diag=T)
  diag(UpperMat) <- 0
  IndexAssoc <- which(UpperMat >0) #matrix index with associations
  
  ##Build edge list
  edgelist <- data.frame(matrix(ncol=3, nrow=length(IndexAssoc)))
  names(edgelist)<-c("to","from", "weight")
  
  for(i in 1:length(IndexAssoc)){
    MatIndex <- arrayInd(IndexAssoc[i], dim(UpperMat))
    edgelist[i,"to"] <- rownames(UpperMat)[MatIndex[,1]]
    edgelist[i, "from"] <- colnames(UpperMat)[MatIndex[,2]]
    SRI <- UpperMat[IndexAssoc[i]]/sum(TotalOcc[MatIndex[,1]], TotalOcc[MatIndex[,2]])
    edgelist[i, "weight"] <- SRI
    
  }
  return(edgelist)
  
}