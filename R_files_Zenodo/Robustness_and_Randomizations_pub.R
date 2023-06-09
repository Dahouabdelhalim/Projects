#Robustness calculations and randomizations

#Use this script to calculate robustness estimates of groups and communities, as well as calculate whether 
#community structure in observed networks is more structured than expected by random chance. 

#Robustness code based off of Shizuka and Farine paper
#Randomization code based off Farine's asnipe package


####*Robustness of Groups noK####

###Files required:
#GBInoK for each year
#BirdlistnoK for each year

library(here)

#Load GBI
gbi16noK = read.csv(here::here("Output files","gbi16noK.csv"))
gbi17noK = read.csv(here::here("Output files","gbi17noK.csv"))
gbi18noK = read.csv(here::here("Output files","gbi18noK.csv"))
gbi19noK = read.csv(here::here("Output files","gbi19noK.csv"))

#Load Birdlists
birdlist16noKsub = read.csv(here::here("Output files","birdlist16noKsub.csv"))
birdlist17noKsub = read.csv(here::here("Output files","birdlist17noKsub.csv"))
birdlist18noKsub = read.csv(here::here("Output files","birdlist18noKsub.csv"))
birdlist19noKsub = read.csv(here::here("Output files","birdlist19noKsub.csv"))


library(igraph)
library(asnipe)
library(assortnet)
library(fpc)

#Function to calculate r_c, with default number of bootstraps = 100, and
#default option to not plot result. Plot will be saved as pdf file in R output folder as
#"rc_result.pdf". 
#Function modified from Shizuka and Farine paper

calc_rc=function(data, year, n.bootstraps=100, plot.result=F){
  
  # 1. Calculate network
  network <- get_network(data,data_format="GBI", association_index="SRI")
  network = network[order(rownames(network)),order(colnames(network))]  
  
  # 2. Subset the network to individuals seen enough times in original network
  if(year==2016) {network = network[birdlist16noKsub$Bird,birdlist16noKsub$Bird]}
  if(year==2017) {network = network[birdlist17noKsub$Bird,birdlist17noKsub$Bird]}
  if(year==2018) {network = network[birdlist18noKsub$Bird,birdlist18noKsub$Bird]}
  if(year==2019) {network = network[birdlist19noKsub$Bird,birdlist19noKsub$Bird]}
  
  
  # 3. Create space to store results from bootstraps
  network.community <- matrix(0,ncol(network),ncol(network)) #columns is # of individuals
  network.present <- matrix(0,ncol(network),ncol(network))
  
  # 4. Calculate community membership of the observed network
  networkA = sweep(network,1,1) #Subtract one from each association index to get association distance
  networkB = sweep(networkA,1,-1,FUN = "*") #Make it positive 
  dm = as.dist(networkB)
  fit <- hclust(dm, method="average") 
  
  #Get list of # groups by height
  groupdata = data.frame(seq(1,0,by=-0.01))
  colnames(groupdata) = "Height"
  groupdata$NumGroups = NA
  for (i in 1:nrow(groupdata)) {
    groups = cutree(fit,h=groupdata$Height[i])
    grouplist = data.frame(table(groups))
    groupdata$NumGroups[i] = nrow(grouplist)
  }
  groupdata = groupdata[!duplicated(groupdata$NumGroups),]
  
  #Figure out where to cut each dendrogram using silhouette distance
  groupdata$average.within = NA
  groupdata$avg.silwidth = NA
  #Get an error for first and last cutpoints, so that's why 2: and nrow -1
  for (i in 2:(nrow(groupdata)-1)) { 
    groups = data.frame(cutree(fit, h=groupdata$Height[i]))
    colnames(groups) = "groupID"
    groupdata$average.within[i] = cluster.stats(networkB,groups$groupID)$average.within
    groupdata$avg.silwidth[i] = cluster.stats(networkB,groups$groupID)$avg.silwidth
  }
  groupdata = groupdata[order(groupdata$avg.silwidth, decreasing =T),]
  h_obs = groupdata$Height[1]
  group.obs = cutree(fit,h=h_obs)
  
  # 3. Main bootstrapping method: i) Bootstrap the observed data, ii) recalculate the network, 
  #    iii) recalculate group membership, iv) check if both individuals are observed
  
  for (i in 1:n.bootstraps) {
    # This step bootstraps the sampling periods
    gbi.boot <- data[sample(1:nrow(data),nrow(data),replace=TRUE),]
    network.boot <- get_network(gbi.boot,data_format="GBI", association_index="SRI")
    network.boot = network.boot[order(rownames(network.boot)),order(colnames(network.boot))]  
    
    if(year==2016) {network.boot = network.boot[birdlist16noKsub$Bird,birdlist16noKsub$Bird]}
    if(year==2017) {network.boot = network.boot[birdlist17noKsub$Bird,birdlist17noKsub$Bird]}
    if(year==2018) {network.boot = network.boot[birdlist18noKsub$Bird,birdlist18noKsub$Bird]}
    if(year==2019) {network.boot = network.boot[birdlist19noKsub$Bird,birdlist19noKsub$Bird]}
    
    # This step calculates the community membership from the bootstrapped network
    network.bootA = sweep(network.boot,1,1) #Subtract one from each association index to get association distance
    network.bootB = sweep(network.bootA,1,-1,FUN = "*") #Make it positive 
    dm.boot = as.dist(network.bootB)
    fit.boot <- hclust(dm.boot, method="average") 
    
    #Get list of # groups by height
    groupdata.boot = data.frame(seq(1,0,by=-0.01))
    colnames(groupdata.boot) = "Height"
    groupdata.boot$NumGroups = NA
    for (i in 1:nrow(groupdata.boot)) {
      groups = cutree(fit.boot,h=groupdata.boot$Height[i])
      grouplist = data.frame(table(groups))
      groupdata.boot$NumGroups[i] = nrow(grouplist)
    }
    groupdata.boot = groupdata.boot[!duplicated(groupdata.boot$NumGroups),]
    
    #Figure out where to cut each dendrogram using silhouette distance
    groupdata.boot$average.within = NA
    groupdata.boot$avg.silwidth = NA
    #Get an error for first and last cutpoints, so that's why 2: and nrow -1
    for (i in 2:(nrow(groupdata.boot)-1)) { 
      groups = data.frame(cutree(fit.boot, h=groupdata.boot$Height[i]))
      colnames(groups) = "groupID"
      groupdata.boot$average.within[i] = cluster.stats(network.bootB,groups$groupID)$average.within
      groupdata.boot$avg.silwidth[i] = cluster.stats(network.bootB,groups$groupID)$avg.silwidth
    }
    groupdata.boot = groupdata.boot[order(groupdata.boot$avg.silwidth, decreasing =T),]
    h_boot = groupdata.boot$Height[1]
    group.boot = cutree(fit.boot, h=h_boot)
    
    # This step adds 1 to any dyads in the same community
    # Joe: Just says, are they in the same community, true or false. If true, add 1 to that pair of birds
    network.community <- network.community + outer(group.boot, group.boot,"==")
    #Joe: can use any list here, outer puts one list across top, another across side, and makes a matrix 
    #to compare them
    
    # This step adds 1 to any dyads that are both present (in this case if they have at least 1 edge)
    #Joe: This is a test to make sure that both birds in each dyad were actually present in the 
    #sampled gbi. If they were not both present, then you can't compare whether they're in the same
    #group or not. Dai and Damien's notation is confusing, it's not saying that every dyad has an edge,
    #It's saying that every bird has an edge. 
    network.present <- network.present + outer((rowSums(network.boot)>0),(rowSums(network.boot)>0),"*")
    
  }
  # End bootstrap
  
  # Calculate proportion of times observed in the same community
  #Joe: for each time two birds were in the sampled data, how often were they in the same community? 
  P <- network.community/network.present
  P[!is.finite(P)] <- 0 #make the NAs = 0 
  
  
  # Calculate assortment from known community membership
  rc <- assortment.discrete(P,group.obs)$r
  
  #if the argument plot.result=T, then generate plot network of probabilities that 
  #nodes are assigned to the same community in bootstraps. It will be saved as pdf 
  #file called "rc_result.pdf" in your R output folder
  if(plot.result) {
    pdf("rc_result.pdf")
    diag(P)=0
    g=graph.adjacency(P, "undirected", weighted=T)
    plot(g, edge.width=E(g)$weight, vertex.label="", vertex.size=5, 
         vertex.color=(group.obs))
    dev.off()
  }
  
  return(rc)
}
#end function


#calc_rc(gbi16noK,year=2016, n.bootstraps=1000, plot.result=F) #0.969 
#calc_rc(gbi17noK,year=2017, n.bootstraps=1000, plot.result=F) #0.9412019
#calc_rc(gbi18noK,year=2018, n.bootstraps=1000, plot.result=F) #0.9714841
#calc_rc(gbi19noK,year=2019, n.bootstraps=1000, plot.result=F) #0.9812271


#__________________________________________________________________________#



####Group structure vs Random Networks####

###Files required:
#GBInoK for each year
#BirdlistnoK for each year
#Observed silhouette width for each year


#Load GBI
gbi16noK = read.csv(here::here("Output files","gbi16noK.csv"))
gbi17noK = read.csv(here::here("Output files","gbi17noK.csv"))
gbi18noK = read.csv(here::here("Output files","gbi18noK.csv"))
gbi19noK = read.csv(here::here("Output files","gbi19noK.csv"))

#Load Birdlists
birdlist16noKsub = read.csv(here::here("Output files","birdlist16noKsub.csv"))
birdlist17noKsub = read.csv(here::here("Output files","birdlist17noKsub.csv"))
birdlist18noKsub = read.csv(here::here("Output files","birdlist18noKsub.csv"))
birdlist19noKsub = read.csv(here::here("Output files","birdlist19noKsub.csv"))

#Load average silhoutte widths
avg.sw = read.csv(here::here("Output files","AverageSilhouetteWidths.csv"))
avg.s16 = avg.sw[which(avg.sw$year==2016),]$avg.s
avg.s17 = avg.sw[which(avg.sw$year==2017),]$avg.s
avg.s18 = avg.sw[which(avg.sw$year==2018),]$avg.s
avg.s19 = avg.sw[which(avg.sw$year==2019),]$avg.s


library(igraph)
library(asnipe)
library(assortnet)
library(fpc)


#Function to compare average silhouette width of random networks to average sil width of observed network
observed_vs_random <- function(gbi,perm.number,birdlistnames,observed.sil.width,year) {
  #Get names of birds in GBI
  names = colnames(gbi)
  #Make random networks. Don't worry that it says no association matrix provided - it calculates 
  #it from the gbi
  random_networks <- network_permutation(gbi,data_format="GBI", identities = names, 
                                         permutations=perm.number, association_index = "SRI")
  #Then subset the random networks based on the list of birds seen enough times
  #This method swaps observations of two birds - so the number of times a bird was seen does not change, 
  #only the observations it was seen in change.
  rand.net = random_networks[,which(names %in% birdlistnames),which(names %in% birdlistnames)]
  #Then create dataframe to store avg silhouette widths of random networks
  sil.df = data.frame(matrix(NA, nrow = perm.number, ncol = 1))
  colnames(sil.df) = "Avg.S"
  #Then calculate the silhouette widths for each random network
  for (h in c(1:perm.number)) {
    randn = rand.net[h,,] #For each random network
    network.bootA = sweep(randn,1,1) #Subtract one from each association index to get association distance
    network.bootB = sweep(network.bootA,1,-1,FUN = "*") #Make it positive 
    dm.boot = as.dist(network.bootB) #make it a distance matrix
    fit.boot <- hclust(dm.boot, method="average") #get dendrogram
    
    #Get list of # groups by height
    groupdata.boot = data.frame(seq(1,0,by=-0.01))
    colnames(groupdata.boot) = "Height"
    groupdata.boot$NumGroups = NA
    for (i in 1:nrow(groupdata.boot)) {
      groups = cutree(fit.boot,h=groupdata.boot$Height[i])
      grouplist = data.frame(table(groups))
      groupdata.boot$NumGroups[i] = nrow(grouplist)
    }
    groupdata.boot = groupdata.boot[!duplicated(groupdata.boot$NumGroups),]
    
    #Figure out where to cut each dendrogram using silhouette distance
    groupdata.boot$avg.silwidth = NA
    #Get an error for first and last cutpoints, so that's why 2: and nrow -1
    for (i in 2:(nrow(groupdata.boot)-1)) { 
      groups = data.frame(cutree(fit.boot, h=groupdata.boot$Height[i]))
      colnames(groups) = "groupID"
      groupdata.boot$avg.silwidth[i] = cluster.stats(network.bootB,groups$groupID)$avg.silwidth
    }
    #plot(groupdata.boot$Height,groupdata.boot$avg.silwidth, pch=19)
    #lines(groupdata.boot$Height,groupdata.boot$avg.silwidth)
    groupdata.boot = groupdata.boot[order(groupdata.boot$avg.silwidth, decreasing =T),]
    sil.df$Avg.S[h] = groupdata.boot$avg.silwidth[1]
  }
  #Then plot histogram of observed compared to random networks
  hist(sil.df$Avg.S, xlim=c(min(sil.df$Avg.S)-0.001,observed.sil.width+0.001), xlab="Silhouette Width",
       main=year,breaks = 25)
  abline(v=observed.sil.width, lwd=3,col="red")
  P = 1 - sum(observed.sil.width > sil.df$Avg.S)/perm.number
  cat("\\n")
  print(P)
}

par(xpd=FALSE)
par(mfrow=c(2,2))
# #2016 p<0.001
# observed_vs_random(gbi=gbi16noK,perm.number = 1000,birdlistnames = birdlist16noKsub$Bird,
#                      observed.sil.width = avg.s16,year="2016")
# # #2017 p=0.001
# observed_vs_random(gbi=gbi17noK,perm.number = 1000,birdlistnames = birdlist17noKsub$Bird,
#                     observed.sil.width = avg.s17,year="2017")
# # #2018 p<0.001
# observed_vs_random(gbi=gbi18noK,perm.number = 1000,birdlistnames = birdlist18noKsub$Bird,
#                     observed.sil.width = avg.s18,year="2018")
# # #2019 p<0.001
# observed_vs_random(gbi=gbi19noK,perm.number = 1000,birdlistnames = birdlist19noKsub$Bird,
#                     observed.sil.width = avg.s19,year="2019")
par(mfrow=c(1,1))




#__________________________________________________________________________#




####*Robustness of Communities####

#This analysis looks to see if communities are robust to sampling error. It assumes that social group structure is 
#robust (which we've shown) so when bootstrap is done to gbi, individuals are still assumed to belong to their
#same group as in the observed data. Then the communities are calculated to see if the groups form the same communities
#in the bootstrapped data. Initially thought that maybe re-calculating social groups before calculating robustness of
#communities would be a good idea, but if birds are in different social groups from before, how do you tell if the 
#social communities are robust if different birds are in the groups making up the communities? Since social groups are 
#extremely robust to sampling error, can assume birds are in same social groups and see how robust communities are
#to sampling error. 


###Files required:
#gbi for each year (kerfuffle sightings included)
#reduced association matrix
#group membership
#community membership
#birdlistnoKc for each year - birdlist of birds not in single-bird groups


#Load GBI with kerfuffle sightings included (all sightings)
gbi16 = read.csv(here::here("Output files","gbi16.csv"))
gbi17 = read.csv(here::here("Output files","gbi17.csv"))
gbi18 = read.csv(here::here("Output files","gbi18.csv"))
gbi19 = read.csv(here::here("Output files","gbi19.csv"))

##Load reduced association matrix 
#check.names=F means it will not add "X" to the column names
rnet16 = read.csv(here::here("Output files","rnet16.csv"),row.names = 1, check.names=F)
rnet17 = read.csv(here::here("Output files","rnet17.csv"),row.names = 1, check.names=F)
rnet18 = read.csv(here::here("Output files","rnet18.csv"),row.names = 1, check.names=F)
rnet19 = read.csv(here::here("Output files","rnet19.csv"),row.names = 1, check.names=F)

#Load group membership
group16 = read.csv(here::here("Output files","group16.csv"))
group17 = read.csv(here::here("Output files","group17.csv"))
group18 = read.csv(here::here("Output files","group18.csv"))
group19 = read.csv(here::here("Output files","group19.csv"))

#Load community membership
comm16 = read.csv(here::here("Output files","communities16.csv"))
comm17 = read.csv(here::here("Output files","communities17.csv"))
comm18 = read.csv(here::here("Output files","communities18.csv"))
comm19 = read.csv(here::here("Output files","communities19.csv"))

#Load birdlists with birds in communities - birds in a group on their own have been removed
birdlist16noKsubc = read.csv(here::here("Output files","birdlist16noKsubc.csv"))
birdlist17noKsubc = read.csv(here::here("Output files","birdlist17noKsubc.csv"))
birdlist18noKsubc = read.csv(here::here("Output files","birdlist18noKsubc.csv"))
birdlist19noKsubc = read.csv(here::here("Output files","birdlist19noKsubc.csv"))

library(igraph)
library(asnipe)
library(assortnet)
library(fpc)
library(reshape2)
library(tibble)

#Function to calculate r_c, with default number of bootstraps = 100, and
#default option to not plot result. Plot will be saved as pdf file in R output folder as
#"rc_result.pdf".

calc_rc_com=function(fullgbi, reduced.matrix, grp.membership, com.membership, birdlistnames, 
                     n.bootstraps=100, plot.result=F){
  
  # 1. Get network
  network = reduced.matrix  
  
  # 2. Create space to store results from bootstraps
  network.community <- matrix(0,ncol(network),ncol(network)) #columns is # of social groups
  network.present <- matrix(0,ncol(network),ncol(network))
  
  # 3. Main bootstrapping method: i) Bootstrap the observed data, ii) recalculate the network, 
  #    iii) recalculate group membership, iv) check if both individuals are observed
  
  for (i in 1:n.bootstraps) {
    # This step bootstraps the sampling periods
    gbi.boot <- fullgbi[sample(1:nrow(fullgbi),nrow(fullgbi),replace=TRUE),]
    network.boot <- get_network(gbi.boot,data_format="GBI", association_index="SRI")
    network.boot = network.boot[order(rownames(network.boot)),order(colnames(network.boot))]
    #Subset network to birds that went into the original communities analysis - all of these birds 
    #were seen enough times and are not in groups on their own. 
    names = colnames(network.boot) 
    network.boot = network.boot[which(names %in% birdlistnames),which(names %in% birdlistnames)]
    
    #Melt and add groups
    rnetmelt = melt(network.boot)
    
    #Match each variable (bird) to its group
    rnetmelt = merge(rnetmelt,grp.membership,by.x="Var1",by.y="Bird")
    rnetmelt = merge(rnetmelt,grp.membership,by.x="Var2",by.y="Bird")
    
    #cast into reduced network
    rnet = as.matrix(acast(rnetmelt,Social.Group.x~Social.Group.y,value.var = "value",fun.aggregate = sum))
    diag(rnet) = 0
    
    #calculate modularity community membership
    rneti = graph.adjacency(rnet, mode="undirected", diag=FALSE, weighted=TRUE)
    rfg = fastgreedy.community(rneti)
    rfgm = data.frame(print(membership(rfg)))
    rfgm = rownames_to_column(rfgm,var="Group")
    colnames(rfgm)[2]="Community"
    
    # This step adds 1 to any dyads in the same community
    # Joe: Just says, are they in the same community, true or false. If true, add 1 to that pair of birds
    network.community <- network.community + outer(rfgm$Community, rfgm$Community,"==")
    #Joe: can use any list here, outer puts one list across top, another across side, and makes a matrix 
    #to compare them
    
    # This step adds 1 to any dyads that are both present (in this case if they have at least 1 edge)
    #Joe: This is a test to make sure that both birds in each dyad were actually present in the 
    #sampled gbi. If they were not both present, then you can't compare whether they're in the same
    #group or not. Dai and Damien's notation is confusing, it's not saying that every dyad has an edge,
    #It's saying that every bird has an edge. 
    network.present <- network.present + outer((rowSums(rnet)>0),(rowSums(rnet)>0),"*")
    
  }
  # End bootstrap
  
  # Calculate proportion of times observed in the same community
  #Joe: for each time two birds were in the sampled data, how often were they in the same community? 
  P <- network.community/network.present
  P[!is.finite(P)] <- 0 #make the NAs = 0 
  
  # Calculate assortment from known community membership
  rc <- assortment.discrete(P,com.membership)$r
  
  #if the argument plot.result=T, then generate plot network of probabilities that 
  #nodes are assigned to the same community in bootstraps. It will be saved as pdf 
  #file called "rc_result.pdf" in your R output folder
  if(plot.result) {
    pdf("rc_result.pdf")
    diag(P)=0
    g=graph.adjacency(P, "undirected", weighted=T)
    plot(g, edge.width=E(g)$weight, vertex.label="", vertex.size=5, 
         vertex.color=(group.obs))
    dev.off()
  }
  
  return(rc)
}
#end function

# #2016 = 0.83
# calc_rc_com(fullgbi = gbi16, reduced.matrix = rnet16, grp.membership = group16,
# com.membership = comm16$Community, birdlistnames = birdlist16noKsubc$Bird,
# n.bootstraps=1000, plot.result=F)
# 
# #2017 = 0.81
# calc_rc_com(fullgbi = gbi17, reduced.matrix = rnet17, grp.membership = group17,
# com.membership = comm17$Community, birdlistnames = birdlist17noKsubc$Bird,
# n.bootstraps=1000, plot.result=F)
# 
# #2018 = 0.90
# calc_rc_com(fullgbi = gbi18, reduced.matrix = rnet18, grp.membership = group18,
# com.membership = comm18$Community, birdlistnames = birdlist18noKsubc$Bird,
# n.bootstraps=1000, plot.result=F)
# 
# #2019 = 0.87
# calc_rc_com(fullgbi = gbi19, reduced.matrix = rnet19, grp.membership = group19,
# com.membership = comm19$Community, birdlistnames = birdlist19noKsubc$Bird,
# n.bootstraps=1000, plot.result=F)




#__________________________________________________________________________#




####Community structure vs Random Networks####

###Files required
#bird files for each year - files that lead into calculating gbi
#gbi for each year (kerfuffle sightings included)
#birdlistnoK files
#observed modularity values


#Load bird files for each year
bird16 = read.csv(here::here("Input files","bird16.csv"))
bird17 = read.csv(here::here("Input files","bird17R.csv"))
bird18 = read.csv(here::here("Input files","bird18.csv"))
bird19 = read.csv(here::here("Input files","bird19.csv"))

#Load GBI with kerfuffle sightings included (all sightings)
gbi16 = read.csv(here::here("Output files","gbi16.csv"))
gbi17 = read.csv(here::here("Output files","gbi17.csv"))
gbi18 = read.csv(here::here("Output files","gbi18.csv"))
gbi19 = read.csv(here::here("Output files","gbi19.csv"))

#Load Birdlists
birdlist16noKsub = read.csv(here::here("Output files","birdlist16noKsub.csv"))
birdlist17noKsub = read.csv(here::here("Output files","birdlist17noKsub.csv"))
birdlist18noKsub = read.csv(here::here("Output files","birdlist18noKsub.csv"))
birdlist19noKsub = read.csv(here::here("Output files","birdlist19noKsub.csv"))

#Load modularity scores from observed networks
mod.all = read.csv(here::here("Output files","modularity_scores.csv"))
modularity16 = mod.all[which(mod.all$year==2016),]$modularity
modularity17 = mod.all[which(mod.all$year==2017),]$modularity
modularity18 = mod.all[which(mod.all$year==2018),]$modularity
modularity19 = mod.all[which(mod.all$year==2019),]$modularity

##Get kerfuffle sampling points for each year 
#First remove duplicate sightings which will give same number of observations as in GBI
#Then get a numbered column that will match column numbers in GBI
bird16noD = bird16[!duplicated(bird16$Sighting),]
bird16noD$rowname = seq(1:nrow(bird16noD))
bird17noD = bird17[!duplicated(bird17$Sighting),]
bird17noD$rowname = seq(1:nrow(bird17noD))
bird18noD = bird18[!duplicated(bird18$Sighting),]
bird18noD$rowname = seq(1:nrow(bird18noD))
bird19noD = bird19[!duplicated(bird19$Sighting),]
bird19noD$rowname = seq(1:nrow(bird19noD))

#Then get kerfuffle points
ksp16 = bird16noD[which(bird16noD$Kerfuffle=="Y"),]$rowname
ksp17 = bird17noD[which(bird17noD$Kerfuffle=="Y"),]$rowname
ksp18 = bird18noD[which(bird18noD$Kerfuffle=="Y"),]$rowname
ksp19 = bird19noD[which(bird19noD$Kerfuffle=="Y"),]$rowname


###Run the function
#This is a function to compare modularity of random network communities to modularity of observed network.
#I built it off of Damien's network_permutation function, heavily borrowed a lot of code from that.
#Basic idea is randomize full network, then subset kerfuffle sightings from randomized GBI, then calculate
#social groups from subsetted network, then reduce full network using social groups as nodes, finally 
#calculate modularity of reduced random network. Kerfuffle sightings are those associated with courtship 
#and aggression.

observed_vs_random_com <- function(gbi,birdlistnames,observed.modularity,permutations,ksp,year) {
  #Dataframe to store modularity values of random networks 
  randommods = data.frame(matrix(NA,nrow=permutations))
  colnames(randommods) = "randommod"
  #Get names of birds in GBI
  names = colnames(gbi)
  #Get association matrix
  net = get_network(gbi,data_format = "GBI")
  #Get subsetted gbi - take out kerfuffle points and get network from that gbi
  gbinoK = gbi[!rownames(gbi) %in% ksp,]
  netnoK = get_network(gbinoK,data_format = "GBI")
  
  ###Permutations of gbi
  gbi_perm = gbi
  net2 = net
  net3 = netnoK
  
  #Network function
  #This recalculates the SRI for the two birds that get swapped in the GBI observations below and any birds 
  #connected to them.
  #input is gbi_perm, GroupBy is the column number. First line adds one to the rows where the bird that got 
  #swapped was present in. Now any bird associating with the focal bird will have a value of 2. 
  #Any bird that was not present will have a value of 1. 
  do.SR_perm <- function(GroupBy,input) {
    tmp <- input[ ,GroupBy] + input
    x <- colSums(tmp==2)
    yab <- colSums(tmp==1)
    out <- (x / (x + yab))
    out
  }
  
  #Do permutations of GBI
  #I think "break" breaks out of the repeat loop if the conditions below are met. Code from network_permutation
  #Permute the full GBI. net2 and net3 are going to change with each swap. Each swap gets added onto the network
  #on top of the others that have come before it. 
  for (n in 1:permutations) {
    repeat {
      a <- which(gbi_perm>0,arr.ind=TRUE)
      s <- sample(1:nrow(a),1)
      first <- a[s,]
      second <- a[sample(1:nrow(a),1),]
      if (first[1]!=second[1]&first[2]!=second[2]&(sum(gbi_perm[first[1],first[2]]) > 0 & sum(gbi_perm[second[1],second[2]]) > 0) & (sum(gbi_perm[second[1],first[2]]) == 0 & sum(gbi_perm[first[1],second[2]]) == 0)) { break; }
    }
    
    #Do the permutations/make the swaps and set the old values to 0. 
    gbi_perm[second[1],first[2]] <- gbi_perm[first[1],first[2]]
    gbi_perm[first[1],second[2]] <- gbi_perm[second[1],second[2]]
    gbi_perm[first[1],first[2]] <- 0
    gbi_perm[second[1],second[2]] <- 0
    
    #Subset kerfuffle sightings from gbi permutation
    #Using list of sample points (Group from BirdI) associated with kerfuffles, take out same 
    #sample points from this GBI. Can then use this GBI to calculate social groups without the kerffufle pts
    gbi_perm_sub = gbi_perm[!rownames(gbi_perm) %in% ksp,]
    
    #Calculate new SRI values for birds associating with swapped birds
    tmp1 <- do.SR_perm(first[2],gbi_perm)
    tmp2 <- do.SR_perm(second[2],gbi_perm)
    tmp3 <- do.SR_perm(first[2],gbi_perm_sub)
    tmp4 <- do.SR_perm(second[2],gbi_perm_sub)
    
    #Add in new SRI values to association matrix, in both columns and rows
    net2[,first[2]] <- tmp1 
    net2[first[2],] <- tmp1
    net2[,second[2]] <- tmp2
    net2[second[2],] <- tmp2
    
    #Do the same for subsetted association matrix
    net3[,first[2]] <- tmp3
    net3[first[2],] <- tmp3
    net3[,second[2]] <- tmp4
    net3[second[2],] <- tmp4
    
    #Make the diagonals of the association matrices equal to 0.
    diag(net2) <- 0
    diag(net3) <- 0
    
    #Subset both association matrices to birds seen enough times
    #These are birds that were seen enough times in the noK network. Because they had been seen enough 
    #times we could assign social groups to them. 
    net2.s = net2[which(names %in% birdlistnames),which(names %in% birdlistnames)] 
    net3.s = net3[which(names %in% birdlistnames),which(names %in% birdlistnames)] 
    
    #Calculate social groups for net3 - no kerfuffle sightings
    network.bootA = sweep(net3.s,1,1) #Subtract one from each association index to get association distance
    network.bootB = sweep(network.bootA,1,-1,FUN = "*") #Make it positive 
    dm.boot = as.dist(network.bootB) #make it a distance matrix
    fit.boot <- hclust(dm.boot, method="average") #get dendrogram
    
    #Get list of # groups by height
    groupdata.boot = data.frame(seq(1,0,by=-0.01))
    colnames(groupdata.boot) = "Height"
    groupdata.boot$NumGroups = NA
    for (i in 1:nrow(groupdata.boot)) {
      groups = cutree(fit.boot,h=groupdata.boot$Height[i])
      grouplist = data.frame(table(groups))
      groupdata.boot$NumGroups[i] = nrow(grouplist)
    }
    groupdata.boot = groupdata.boot[!duplicated(groupdata.boot$NumGroups),]
    
    #Figure out where to cut each dendrogram using silhouette distance
    groupdata.boot$avg.silwidth = NA
    #Get an error for first and last cutpoints, so that's why 2: and nrow -1
    for (i in 2:(nrow(groupdata.boot)-1)) { 
      groups = data.frame(cutree(fit.boot, h=groupdata.boot$Height[i]))
      colnames(groups) = "groupID"
      groupdata.boot$avg.silwidth[i] = cluster.stats(network.bootB,groups$groupID)$avg.silwidth
    }
    groupdata.boot = groupdata.boot[order(groupdata.boot$avg.silwidth, decreasing =T),]
    cutH = groupdata.boot$Height[1]
    
    #Get who's in each group
    groupID = data.frame(cutree(fit.boot, h=cutH))
    groupID = rownames_to_column(groupID)
    colnames(groupID)=c("Bird","Social.Group")
    
    #Remove individuals in their own group from group ID list and network
    grouptable = data.frame(table(groupID$Social.Group))
    grouptable = grouptable[which(grouptable$Freq!=1),]
    groupID = groupID[groupID$Social.Group %in% grouptable$Var1,]
    net2.s = net2.s[groupID$Bird,groupID$Bird]
    
    ###Reduce full network with groups
    #Melt association matrix and match up group identities to each bird
    net2.smelt = melt(net2.s)
    
    #Match each variable (bird) to its group
    net2.smelt = merge(net2.smelt,groupID,by.x="Var1",by.y="Bird")
    net2.smelt = merge(net2.smelt,groupID,by.x="Var2",by.y="Bird")
    
    #Cast it back into a matrix based on Social Group and sum values for like associations
    net2.sr = as.matrix(acast(net2.smelt,Social.Group.x~Social.Group.y,value.var = "value",fun.aggregate = sum))
    diag(net2.sr) = 0 
    
    #Setup for igraph
    net2.sr.ig = graph.adjacency(net2.sr, mode="undirected", diag=FALSE, weighted=TRUE)
    fg.net2 = fastgreedy.community(net2.sr.ig)
    
    #Get modularity value of reduced network and place into dataframe
    randommods$randommod[n] = modularity(fg.net2)
    #print(c("round",n))
  }
  par(xpd=FALSE)
  hist(randommods$randommod,xlim=c(min(randommods$randommod)-0.01,observed.modularity+0.01),
       main=year,xlab="Modularity", breaks=50)
  abline(v=observed.modularity,lwd=3,col="red")
  P = 1 - sum(observed.modularity > randommods$randommod)/permutations
  cat("\\n")
  print(P)
}

par(mfrow=c(2,2))
#Make sure you use the correct birdlist names
# #2016
# observed_vs_random_com(gbi=gbi16, birdlistnames = birdlist16noKsub$Bird, observed.modularity = modularity16,
#                        permutations = 1000, ksp=ksp16,year="2016")
# #2017
# observed_vs_random_com(gbi=gbi17, birdlistnames = birdlist17noKsub$Bird, observed.modularity = modularity17,
#                        permutations = 1000, ksp=ksp17,year="2017")
# #2018
# observed_vs_random_com(gbi=gbi18, birdlistnames = birdlist18noKsub$Bird, observed.modularity = modularity18,
#                        permutations = 1000, ksp=ksp18,year="2018")
# #2019
# observed_vs_random_com(gbi=gbi19, birdlistnames = birdlist19noKsub$Bird, observed.modularity = modularity19,
#                        permutations = 1000, ksp=ksp19,year="2019")
par(mfrow=c(1,1))





#__________________________________________________________________________#





####Mean association index within versus across groups - difference greater than random?####

#Want to know if the difference between the mean association index of within social groups connections
#and the mean association index of across-group connections is greater than expected by random chance. 

#Method: 
#Generate random networks using a datastream permutation using the entire dataset
#Then assume all birds are in the same social group as in the observed network
#Calculate within and across-group association indices and compare difference to observed

#The silhoutte width and social group methods sort of show this already, that they associate more than expected by chance
#but this is another way to back that up. Idea is that all birds are in the same groups as in observed, if you randomized
#the network connections, would there still be a large difference between within versus across group connections. 
#Answer is almost certainly no given how structured the network is. 

#Load in mean SRI's and calculate difference - from Networks_Dendrograms_Groups_Communities.R script
obsSRIwithin = 0.619
obsSRIacross = 0.053
obsdiff = obsSRIwithin - obsSRIacross


#Function to compare average silhouette width of random networks to average sil width of observed network
within_across_vs_random_SRIs <- function(gbi,perm.number,birdlistnoKsubc) {
  
  #First order gbi so random networks will be in same order as group list later
  gbi = gbi[order(colnames(gbi))]
  
  #Get names of birds in GBI
  names = colnames(gbi) 
  
  #Get groups from birdlistnoKsubc because that has the same birds in it as the random networks will. 
  #And get rownames of groups dataframe - random networks do not include rownames and column names - bird IDs. but it is in the 
  #same order as groups because of above
  groups = data.frame(birdlistnoKsubc$Bird,birdlistnoKsubc$Social.Group)
  colnames(groups) = c("Bird","Social.Group")
  groups = rownames_to_column(groups)
  
  #Make random networks. Don't worry that it says no association matrix provided - it calculates 
  #it from the gbi
  random_networks <- network_permutation(gbi,data_format="GBI", identities = names, 
                                         permutations=perm.number, association_index = "SRI")

  #Subset random networks by birdslistnoKsubc - birds seen enough times and not in their own group. 
  #Do this after random networks created so all individuals go into network creation like normal. 
  rand.net = random_networks[,which(names %in% birdlistnoKsubc$Bird),which(names %in% birdlistnoKsubc$Bird)]
  
  #Create empty lists to hold 
  within_list = list()
  across_list = list()
  
  #For each random network, calculate the difference in within vs across group associations 
  for (i in 1:nrow(rand.net)) {
    #Melt random network and merge with group IDs. 
    rnetmelt = melt(rand.net[i,,])
    rnetmelt = merge(rnetmelt,groups,by.x="Var1",by.y="rowname")
    rnetmelt = merge(rnetmelt,groups,by.x="Var2",by.y="rowname")
    
    #Get same social group and remove rows where bird is the same
    rnetmelt.same = rnetmelt[which(rnetmelt$Social.Group.x==rnetmelt$Social.Group.y),] 
    rnetmelt.same = rnetmelt.same[which(rnetmelt.same$Var2!=rnetmelt.same$Var1),] 
    
    #Order within rows then remove duplicate pairs
    rnetmelt.same.birds = data.frame(rnetmelt.same$Var2,rnetmelt.same$Var1)
    rnetmelt.same.birds = data.frame(t(apply(rnetmelt.same.birds,1,sort)))
    rnetmelt.same = data.frame(rnetmelt.same.birds$X1,rnetmelt.same.birds$X2,rnetmelt.same$value,
                               rnetmelt.same$Social.Group.x,rnetmelt.same$Social.Group.y)
    colnames(rnetmelt.same) = c("Bird1","Bird2","SRI","SG1","SG2")
    rnetmelt.same = rnetmelt.same %>% distinct(Bird1,Bird2,.keep_all = T)
    
    within_list[[i]] = rnetmelt.same
    
    #Get different social group and remove rows where bird is the same
    #Make sure connections are greater than 0 - so there needs to be a connection, just not in the same social group
    rnetmelt.different = rnetmelt[which(rnetmelt$Social.Group.x!=rnetmelt$Social.Group.y),] 
    rnetmelt.different = rnetmelt.different[which(rnetmelt.different$value>0),]
    rnetmelt.different = rnetmelt.different[which(rnetmelt.different$Var2!=rnetmelt.different$Var1),] 
    
    #Order within rows then remove duplicate pairs
    rnetmelt.different.birds = data.frame(rnetmelt.different$Var2,rnetmelt.different$Var1)
    rnetmelt.different.birds = data.frame(t(apply(rnetmelt.different.birds,1,sort)))
    rnetmelt.different = data.frame(rnetmelt.different.birds$X1,rnetmelt.different.birds$X2,rnetmelt.different$value,
                                    rnetmelt.different$Social.Group.x,rnetmelt.different$Social.Group.y)
    colnames(rnetmelt.different) = c("Bird1","Bird2","SRI","SG1","SG2")
    rnetmelt.different = rnetmelt.different %>% distinct(Bird1,Bird2,.keep_all = T)
    
    across_list[[i]] = rnetmelt.different
  }
  
  return(list(within_list,across_list))
  
}

#Generate random data for within and across lists for all years
# within_across16 = within_across_vs_random_SRIs(gbi=gbi16,perm.number = 1000,birdlistnoKsubc = birdlist16noKsubc)
# within_across17 = within_across_vs_random_SRIs(gbi=gbi17,perm.number = 1000,birdlistnoKsubc = birdlist17noKsubc)
# within_across18 = within_across_vs_random_SRIs(gbi=gbi18,perm.number = 1000,birdlistnoKsubc = birdlist18noKsubc)
# within_across19 = within_across_vs_random_SRIs(gbi=gbi19,perm.number = 1000,birdlistnoKsubc = birdlist19noKsubc)

#Get random within lists for all years
within16 = within_across16[[1]]
within17 = within_across17[[1]]
within18 = within_across18[[1]]
within19 = within_across19[[1]]

#Calculate mean values for each of these random lists
within_all = matrix(nrow=1000,ncol=1)

for (i in 1:1000) {
  within.i = rbind(within16[[i]],within17[[i]],within18[[i]],within19[[i]])
  within_all[i] = mean(within.i$SRI)
}

#Get random across lists for all years
across16 = within_across16[[2]]
across17 = within_across17[[2]]
across18 = within_across18[[2]]
across19 = within_across19[[2]]

#Calculate mean values for each of these random lists
across_all = matrix(nrow=1000,ncol=1)

for (i in 1:1000) {
  across.i = rbind(across16[[i]],across17[[i]],across18[[i]],across19[[i]])
  across_all[i] = mean(across.i$SRI)
}

#Put together into dataframe and calculate p-value
within_across_all = data.frame(within_all,across_all)
within_across_all$difference = within_across_all$within_all - within_across_all$across_all
hist(within_across_all$difference,main="",breaks=25,xlab="Difference")
abline(v=obsdiff,col="red",lwd=3)

sum(within_across_all$difference > obsdiff)/1000 #p=0.02




