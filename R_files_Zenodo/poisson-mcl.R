# script implemented by Molly Liu with reference to code written by 
#   Connor Chato (https://github.com/PoonLab/clustuneR)
# reviewed and revised by Art Poon

library(igraph)
library(MCL)
library(netdiffuseR)
library(lhs)
library(parallel)
library(fields)
library(data.table)


#### 1. data processing ####
#read file
# file generated as CSV output from TN93
oel <- read.csv(gzfile("tn93.csv.gz"), stringsAsFactors = F)

# separate id and time elements, e.g., MH355537_2002
temp1 <- sapply(oel$ID1, function(x) strsplit(x,'_')[[1]])
temp2 <- sapply(oel$ID2, function(x) strsplit(x,'_')[[1]])

oel2 <- data.frame(
  ID1=as.factor(temp1[1,]),
  t1=as.numeric(temp1[2,]),
  ID2=as.factor(temp2[1,]), 
  t2=as.numeric(temp2[2,]),
  Distance=as.numeric(oel$Distance)
)

# Create a data frame from the imported edge list.
original_edge_list <- data.frame(
  ID1=as.factor(temp1[1,]), #as.character(temp1[1,]), 
  t1=as.numeric(temp1[2,]),
  ID2=as.factor(temp2[1,]), 
  t2=as.numeric(temp2[2,]),
  Distance = as.numeric(oel$Distance)
)


# find node_list
a1<- unique(oel$ID1)
a2 <- unique(oel$ID2)
train_node_list <- unique(append(a1,a2))
train_node_list <- sapply(train_node_list, function(x) strsplit(x,'_')[[1]])
node_list <- data.frame(
  ID1=as.factor(train_node_list[1,]), #as.character(temp1[1,]), 
  t1=as.numeric(train_node_list[2,]))

# delete nodes earlier than 2003
new_node_list <- data.table(subset(node_list, node_list$t1 > 2002))
new_node_list1 <- new_node_list[, .SD[1:400], by=t1]
new_node_list2 <- subset(new_node_list1, new_node_list1$ID1 != "NA")
nodes <- new_node_list2$ID1

new_edge_list <- original_edge_list[original_edge_list$ID1 %in% nodes & 
                                      original_edge_list$ID2 %in% nodes,]

#### end ####

#### 2. apply MCL and connected components method on known cases ####
# 2.1 apply connected components method
ccomponents <- function(matrix) {
  g <- graph_from_adjacency_matrix(matrix)
  clu <- components(g)
  groups_clu <- igraph::groups(clu)
  connectcomp_result <- sapply(1:length( groups_clu), function(i) {list(groups_clu[[i]])})
  return(connectcomp_result)
  }


# 2.2 Applying MCL method
mcls <- function(matrix, whole_data, tMax, exps, infl) {
  mcl_result <- mcl(x = matrix, expansion = exps, inflation = infl, 
                    allow1=TRUE, addLoops=TRUE)
  # use to set 100, delete max.iter to see the result
  # return NA if doesn't mcl function return error
  if(length(mcl_result) != 3) { return(NA) }
  else {
    # get cluster list for each cluster
    a <- mcl_result$Cluster
    cluster_index_list <- lapply(1:max(a), function(i) {which(a==i)})
    a1<- levels(factor(whole_data$ID1))
    a2 <- levels(factor(whole_data$ID2))
    train_node_list <- levels(factor(c(a1,a2)))
    
    #mape node list to id list
    cluster_id_list <- list()
    for (i in 1:length(cluster_index_list))
    {
      cluster_id_list[[i]] <- train_node_list[cluster_index_list[[i]]]
    }
    mcl_cluster <- list()
    mcl_cluster$train_node_list <- train_node_list
    mcl_cluster$cluster_id_list <- cluster_id_list
    return(mcl_cluster)
  }
  
}


#### end ####

#### 3. add new cases to cluster result ####
# there exist a node doesn't connected to any previous node
# cross out the new nodes if thet only connected to each other 

cluster_to_pdata <- function(cluster_result, train_node_list, test_set, tMax,
                             train_data) {
  test_set2 <- subset(test_set, test_set$t1 != test_set$t2)
  #swap the most recent year to ID1
  test_set2a <- subset(test_set2,test_set2$t1 ==  tMax)
  test_set2b <- subset(test_set2,test_set2$t2 ==  tMax)
  step1 <- test_set2b %>% relocate(ID1, .after = t2)
  step2 <- step1 %>% relocate(t1, .after = ID1)
  rename_step2 <- data_frame(ID1 = step2$ID2, t1 = step2$t2,ID2 = step2$ID1,t2 = step2$t1,Distance= step2$Distance)
  test_set3 <- rbind(test_set2a,rename_step2)
  
  # train_node_list <- mcl_cluster$train_node_list``
  node_year <- list()
  for (i in train_node_list){
    if (i %in% train_data$ID1){
      node_year[[i]] <-  original_edge_list[original_edge_list$ID1 == i,]$t1[1]
    }
    if (i %in% train_data$ID2){
      node_year[[i]] <- original_edge_list[original_edge_list$ID2 == i,]$t2[1]
    }
  }
  
  # sort the new cases list and select the one with shortest distance
  test_set_list <- ddply(test_set3, .(ID1), function(x) x[which.min(x$Distance),])
  
  test_set4 = subset(test_set3, test_set3$ID2%in%train_node_list)
  least_distance_list <- ddply(test_set4, .(ID1), function(x) x[which.min(x$Distance),])
  l1 <- as.numeric(lengths(test_set_list[1]))

  # then find where is that node belong to, form p_data result
  # p_data contain 3 columns: # of known cases, # of new cases, how recent
  col1 <- lengths(cluster_result)
  
  colname = c(1:length(cluster_result))
  cluster_mapping = melt(setNames(cluster_result, colname))
  colnames(cluster_mapping) = c("ID2", "cluster")
  cluster_mapping$cluster = as.factor(cluster_mapping$cluster)
  jointable = merge(cluster_mapping,least_distance_list)
  jointable = as.data.table(jointable)
  
  df = as.data.table(table(jointable$cluster))
  colnames(df) = c("cluster","count")
  df$cluster = as.integer(df$cluster)
  
  df =  df[order(df$cluster),]
  col2 = df$count
  
  colname = c(1:length(cluster_result))
  cluster_mapping = melt(setNames(cluster_result, colname))
  colnames(cluster_mapping) = c("node","cluster")
  
  yearofnode = melt(node_year)
  colnames(yearofnode) = c("year","node")
  join_year_cluster = merge(yearofnode, cluster_mapping)
  join_year_cluster$year = tMax-join_year_cluster$year 
  join_year_cluster = as.data.frame(join_year_cluster)
  join_year_cluster = as.data.table(join_year_cluster)
  
  nodessum = join_year_cluster[,.(sum(year)), by = "cluster"]
  nodessum$cluster = as.numeric(nodessum$cluster)
  nodessum = setorder(nodessum,cluster)
  col3 = nodessum$V1
  pdata_result <- data.frame(known_cases = col1,new_cases = col2,how_recent = col3)
  l2 <- sum(pdata_result$new_cases)
  pdata_result$total_case <- l2/l1 
  print(pdata_result$total_case)
  return (pdata_result)
  
}
#### end ####


#### possion regression ####
# each cluster is an observation, the distribution of new cases among those 
# clusters is the outcome
aicfunc <- function(pdata_result) {
  
  # change the one have 0 in know cases
  pdata_remove0 <-subset(pdata_result, pdata_result$known_cases > 0)
  #first poisson regression
  #rhs= number of known cases in each cluster
  #lhs= number of new cases in each cluster
  pmodel1 <- glm(pdata_remove0$new_cases ~ pdata_remove0$known_cases, poisson)
  
  #second poisson regression:
  #rhs= how recent the known cases:recent year-current
  #lhs= number of new cases in each cluster
  # pmodel2 <- glm(pdata_remove0$new_cases ~ pdata_remove0$how_recent, poisson)
  pdata_remove0$how_recent <- pdata_remove0$how_recent / 
    pdata_remove0$known_cases
  pmodel2 <- glm(pdata_remove0$new_cases ~ pdata_remove0$how_recent + 
                   pdata_remove0$known_cases, poisson)
  
  #compare AIC
  aic_value <- pmodel2$aic - pmodel1$aic
  
  return(aic_value)
}

#### end ####

####hypercube sampling ####
#choose number of points
set.seed(2)
h <- 500
#simulatex
lhs <- maximinLHS(h,3)
#minimum and maximum values for each parameter
threshold.min <- 0.0
threshold.max <- 0.06
expansion.min <- 2
expansion.max <- 25
inflation.min <- 2
inflation.max <- 25
# generate a "parameter set" by rescaling our simulated latin hypercube sample
params.set <- cbind(
  threshold = lhs[ ,1] * (threshold.max-threshold.min) + threshold.min, 
  expansion = lhs[ ,2] * (expansion.max-expansion.min) + expansion.min,
  inflation = lhs[ ,3]*(inflation.max-inflation.min) + inflation.min
)
# cycle through different simulated parameter sets and their dAIC 

#### end ####

main <- function(i) {
  # Params at the first colum (params[,1]) is a set of distance thresholds
  # We filter the edge list by whichever param[i,1] is (Row i, column 1)
  params <- as.list((params.set[i,]))
  edge_list <- new_edge_list[which(new_edge_list$Distance < params[[1]]),]
  # separate known cases(train set) and new cases(test set)
  tMax <- max(edge_list["t1"], edge_list["t2"])
  train_set <- subset(edge_list, edge_list["t1"] < tMax & edge_list["t2"] < tMax)
  test_set <- subset(edge_list, edge_list["t1"] > tMax-1 | edge_list["t2"] > tMax-1)
  
  train_data <- data.frame(ID1 = train_set["ID1"],ID2= train_set["ID2"])
  #change the matrix to adjacency matrix(for traning data)
  train_matrix <- edgelist_to_adjmat(train_data)
  
  mcl_cluster <- mcls(train_matrix,train_data,tMax,params[[2]], params[[3]])
  if (any(is.na(mcl_cluster))) { return(c(rep(NA, 3))) }
  else {
    cc_cluster <- ccomponents(train_matrix)
    mcl_pdata <- cluster_to_pdata(mcl_cluster$cluster_id_list, 
                                  mcl_cluster$train_node_list, test_set, tMax,
                                  train_data)
    cc_pdata <- cluster_to_pdata(cc_cluster, mcl_cluster$train_node_list, 
                                 test_set, tMax, train_data)
    mcl_total <- mcl_pdata$total_case
    mcl_aic <- aicfunc(mcl_pdata[,1:3])
    cc_aic <- aicfunc(cc_pdata[,1:3])
    list1 <- list(mcl_aic, cc_aic, mcl_total)
    return(list1)
    }
}


# run sapply from own computer for testing
#res <- sapply(1:2, main)
#run mclapply if using lab computer
t0 <- Sys.time()
#res <- mclapply(1:h, main, mc.cores = 16)  # requires about 5 hours!
load("poisson-mcl.RData")
tf <- difftime(Sys.time(), t0, units = "mins")


aic_result <- data.frame(params.set)
aic_result$mcl_daic <-  as.numeric(lapply(res, function(x) as.numeric(x[1])))
aic_result$cc_daic <- as.numeric(lapply(res, function(x) as.numeric(x[2])))


#remove the runs what doesn't have mcl result
aic_removeNA <- subset(aic_result,(!is.na(aic_result$mcl_daic) | 
                                     !is.na(aic_result$cc_daic) ))

# workspace now contains all objects needed for plotting.R
