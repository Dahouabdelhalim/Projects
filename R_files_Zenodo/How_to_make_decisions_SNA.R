#R script to replicate all the analyses from the manuscript: How to make methodological decisions when inferring social networks. (Ecology and Evolution)

#Convert RFID data to continuous time
RFID_data_to_continuous(RFID_data, tags, tags$Last_colony_captured, same_tags_as_col=TRUE, keeptesttags= TRUE, keepfeeders=TRUE)
rm(RFID_data)
#add time bins to the RIFD
RFID_data_to_birds_fit_time_bin(RFID_continuous_time, 60, tags, duplicates=FALSE, removetest=TRUE)

#Remove test tags and birds not captured
Tag_time_bin_remove_tags(Tag_time_bin, tags, "2017-01-01")

#remove df not needed
rm(RFID_continuous_time)

#save to gmmevents
path="PATH_TO_SAVE_GMMM_EVENTS"

#Change here for the colony that you want to use (colonies: 11,20,27,43,71)
#colonyID
colonyID=71

Gmmevents_within_Gmmevents(Tag_time_bin, colonyID, path, all=TRUE)

########## co-occurrence network ############
#get the network of the co-occurrences
library(asnipe)
net_occorrence <-get_network(gbi2)
cv_occorrence<-cv_network(net_occorrence)

#test if cv is generated randomly

#create a storage matrix for the results
randomcv_occorrence<- matrix(NA, nrow=1000, ncol=1)
networks.rand <- network_permutation(gbi2, association_matrix=net_occorrence, identities = as.character(colnames(net_occorrence)), within_location=FALSE, within_day=FALSE)
#loop
str(networks.rand)
for(i in 1:1000) {
  #extract the association matrix from the networks.rand variable
  net.tmp <- networks.rand[i,,]
  #calculate individuals' eigenvector centralities in the network
  #take this code from the above (on the original netowrk)
  cv_occorrence_rand<-cv_network(net.tmp)
  #now do the same statistical test as above
  randomcv_occorrence[i,1]<-cv_occorrence_rand
  print(i)
}

hist(randomcv_occorrence[,1])
abline(v=cv_occorrence, col="red")
P <- 2*sum(randomcv_occorrence[,1] >= cv_occorrence)/1000
P
plot(randomcv_occorrence)
abline(h=cv_occorrence, col ="red")



##Test assortment
groups<-groups[!is.na(groups$tag),]

#erase groups not duplicated (just one individual)
groups<-groups[ groups$group %in% groups$group[duplicated(groups$group)],]


#get the group ID
groups<-groups[groups$tag %in% colnames(net_occorrence),]
net_occorrence<-net_occorrence[(rownames(net_occorrence) %in% groups$tag),]
net_occorrence<-net_occorrence[,(colnames(net_occorrence) %in% groups$tag)]
net_occorrence<-net_occorrence[order(colnames(net_occorrence)),]
net_occorrence<-net_occorrence[,order(colnames(net_occorrence))]
groups<-groups[order(groups$tag),]

library(asnipe)
library(assortnet)
assort_occorrence <- assortment.discrete(net_occorrence, groups$group, weighted=TRUE, SE=TRUE)
assort_occorrencer<-assort_occorrence$r
#erase individuals that are not in the group ID
gbi2<-gbi2[,colnames(gbi2) %in% groups$tag]
#order gbi2
gbi2<-gbi2[,order(colnames(gbi2))]
networks.rand <- network_permutation(gbi2, association_matrix=net_occorrence, identities = as.character(colnames(net_occorrence)), within_location=FALSE, within_day=FALSE)

#create a storage matrix for the results
randomassort_occorrencer<- matrix(NA, nrow=1000, ncol=1)

#loop
str(networks.rand)
for(i in 1:1000) {
  #extract the association matrix from the networks.rand variable
  net.tmp <- networks.rand[i,,]
  #calculate individuals' eigenvector centralities in the network
  #take this code from the above (on the original netowrk)
  assorti <- assortment.discrete(net.tmp, groups$group, weighted=TRUE, SE=TRUE)
  #now do the same statistical test as above
  randomassort_occorrencer[i,1]<- assorti$r
  print(i)
}

hist(randomassort_occorrencer[,1])
abline(v=assort_occorrencer, col="red")
P <- 2*sum(randomassort_occorrencer[,1] >=assort_occorrencer )/1000
P
plot(randomassort_occorrencer)
abline(h=assort_occorrencer, col ="red")

#########Overlap of time network#############

#convert to start end events
Continuous_time_to_start_end(Tag_time_bin)

#get network of overlapping times
network_overlap_time<-Start_end_to_edgelist(Start_end_events, 5, directed=FALSE)$network_overlap_time

Start_end_events$Start<-as.numeric(Start_end_events$Start)

Start_end_events$End<-as.numeric(Start_end_events$End)
events$Start<-as.numeric(events$Start)
events$End<-as.numeric(events$End)

###########CV  and assortment of the overlap networkcv
cv_overlap_time<-cv_network(network_overlap_time)

#erase groups not duplicated (i.e. it is not a group because it has only one individual)
groups<-groups[ groups$group %in% groups$group[duplicated(groups$group)],]

#set the groups in same order
groups<-groups[(groups$tag %in% rownames(network_overlap_time)),]
network_overlap_time<-network_overlap_time[(rownames(network_overlap_time) %in% groups$tag),]
network_overlap_time<-network_overlap_time[,(colnames(network_overlap_time) %in% groups$tag)]
groups<-groups[order(groups$tag),]
network_overlap_time<-network_overlap_time[order(rownames(network_overlap_time)),]
network_overlap_time<-network_overlap_time[,order(colnames(network_overlap_time))]


#assortment r
network_overlap_time_assortr<-assortment.discrete(network_overlap_time, groups$group, weighted=TRUE, SE=FALSE)$r

#add event ID to the events from the gmm
events$Event<-NA

events<-events[order(events$Start),]

events$Event<-rep(seq_len(nrow(events)%/%1+1L),each=1,len=nrow(events))
rownames(events)<-NULL
events$Start<-events$Start*60
events$End<-events$End*60

#add events to the Start_end_events
Start_end_events$Event<-NA
Start_end_events$notbird<-NA
Start_end_events<-Start_end_events[order(Start_end_events$Start),]
events<-events[order(events$Start),]
rownames(Start_end_events)<-NULL
rownames(events)<-NULL
#add event
for(i in 1:nrow(Start_end_events)){
  if(nrow(subset(events, events$Start<=Start_end_events$Start[i] & events$End>=Start_end_events$Start[i]))>0){
    Start_end_events$Event[i]<-subset(events, events$Start<=Start_end_events$Start[i] & events$End>=Start_end_events$Start[i])$Event
  }
  print(i/nrow(Start_end_events)*100)
}
#last row is na for some reason
Start_end_events<-Start_end_events[!is.na(Start_end_events$Event),]

#plot the events
for(i in 1:nrow(Tag_time_bin1Bin)){
  Tag_time_bin1Bin$countbirds[i]<-nrow(Tag_time_bin1Bin[which(Tag_time_bin1Bin$TimeBin[i]==Tag_time_bin1Bin$TimeBin),])
}
birdrateplot<-seq(min(Tag_time_bin1Bin$TimeBin), max(Tag_time_bin1Bin$TimeBin))
birdrateplot<-data.frame(birdrateplot)
birdrateplot$countbirds<-0
for(i in 1:nrow(birdrateplot)){
  if(nrow(Tag_time_bin1Bin[which(birdrateplot$birdrateplot[i]==Tag_time_bin1Bin$TimeBin),])>0){
    birdrateplot$countbirds[i]<-Tag_time_bin1Bin[which(birdrateplot$birdrateplot[i]==Tag_time_bin1Bin$TimeBin),]$countbirds[1]
  }
}
plot(birdrateplot$birdrateplot, birdrateplot$countbirds, type="l", lwd =1, ylim=c(0, 20))
for(i in 1:nrow(events)){
  segments(events$Start[i]/60, 10, events$End[i]/60, 10,lwd =2, col= 'black')
}

#randomization of the data
random_overlap_cv<- matrix(NA, nrow=1000, ncol=1)
random_overlap_assort<- matrix(NA, nrow=1000, ncol=1)

for(i in 1:1000){
Tagi<-sapply(by(Start_end_events$Tag, Start_end_events$Event, function(x){sample(x, replace = FALSE)}), identity)
Tagi<-as.character(unlist(Tagi))
Start_end_eventsi<-Start_end_events
Start_end_eventsi$Tag<-as.character(Tagi)
net_overlapi<-Start_end_to_edgelist(Start_end_eventsi,5,directed=FALSE)$network_overlap_time
cvi<-cv_network(net_overlapi)
random_overlap_cv[i,1]<-cvi
groupsi<-groups[(groups$tag %in% rownames(net_overlapi)),]
net_overlapi<-net_overlapi[(rownames(net_overlapi) %in% groupsi$tag),]
net_overlapi<-net_overlapi[,(colnames(net_overlapi) %in% groupsi$tag)]
groupsi<-groupsi[order(groupsi$tag),]
net_overlapi<-net_overlapi[order(rownames(net_overlapi)),]
net_overlapi<-net_overlapi[,order(colnames(net_overlapi))]
random_overlap_assort[i,1]<-assortment.discrete(net_overlapi, groupsi$group, weighted=TRUE, SE=FALSE)$r
print(i)
}

hist(random_overlap_cv[,1])
abline(v=cv_overlap_time, col="red")
P <- 2*sum(random_overlap_cv[,1] <=cv_overlap_time )/1000
P


hist(random_overlap_assort[,1])
abline(v=network_overlap_time_assortr, col="red")
P <- 2*sum(random_overlap_assort[,1] >=network_overlap_time_assortr )/1000
P