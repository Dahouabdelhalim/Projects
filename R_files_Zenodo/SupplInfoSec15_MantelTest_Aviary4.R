

##### Load necessary packages
library(asnipe)
library(igraph)
library(sna)
library(assortnet)
library(qgraph)



#### Adult-recipientMother ####

##### Load the data
Event_All_Count <- read.csv("Event_All_Count.csv",stringsAsFactors=FALSE, sep=",")

# chick info
dataChickID <- read.csv("dataChickID_toPublish20200223.csv",stringsAsFactors=FALSE,sep=",")

# parent info
dataParentsID <- read.csv("dataParentsID.csv",stringsAsFactors=FALSE,sep=",")


Event_All_Count_A1<-Event_All_Count[which(Event_All_Count$Aviary=="4"),]
Event_All_Count_A1 <- Event_All_Count_A1[which(Event_All_Count_A1$Category=="Unrelated"),]


## Which adults and offspring of adults (herein, recipient adults) were observed together?
SubjectAdult<-NULL
AllAdult<-unique(Event_All_Count_A1$AdultID)
Present_RecipientAdult_temp2<-NULL
for(i in 1:length(AllAdult)){
  temp_Event_All_Count_A1<-Event_All_Count_A1[which(Event_All_Count_A1$AdultID==AllAdult[i]),]
  SubjectAdult_temp<-NULL
  Present_RecipientAdult_temp<-unique(unlist(strsplit(temp_Event_All_Count_A1$MotherOfPresentNonOffspring," ")))
  length_old<-length(Present_RecipientAdult_temp2)
  for(j in 1:length(Present_RecipientAdult_temp)){
    temp_RecipientAdult<-Present_RecipientAdult_temp[j]
    if(AllAdult[i]!=temp_RecipientAdult){
      Present_RecipientAdult_temp2<-c(Present_RecipientAdult_temp2,temp_RecipientAdult)
    }
  }
  length_new<-length(Present_RecipientAdult_temp2)
  length_to_add<-length_new-length_old
  if(length(Present_RecipientAdult_temp2)>0){
    SubjectAdult_temp<-rep(AllAdult[i],length_to_add)
  }
  SubjectAdult<-c(SubjectAdult,SubjectAdult_temp)
}
length(SubjectAdult)==length(Present_RecipientAdult_temp2)#must be TRUE
ObservedCombination_FedNonFed_female<-data.frame(SubjectAdult=SubjectAdult,Present_RecipientAdult=Present_RecipientAdult_temp2)







##### [Pre-breeding Social Network *include only adults in the provision network]
# Load the data
load("2018_A4_prebreedingMtx.RData") #prebreedmtx
# Change the column names from ring number to backpack
individuals_temp<-NULL
for(i in 1:ncol(prebreedmtx)){
  backpack_temp<-dataParentsID$qr_code[which(colnames(prebreedmtx)[i]==dataParentsID$ind_ID)]
  individuals_temp<-c(individuals_temp,backpack_temp)
}
rownames(prebreedmtx)<-individuals_temp
colnames(prebreedmtx)<-rownames(prebreedmtx)
# Plot the network(igraph)
net<-graph.adjacency(prebreedmtx,mode="undirected", weighted=TRUE,diag=FALSE)
e <- get.edgelist(net, names=F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net), area=8*(vcount(net)^2), repulse.rad = vcount(net)^3.6)
#add colours to adults
A1_individuals<-data.frame(Individuals=individuals_temp, Sex=NA)
A1_individuals$Sex_Colour<-"green"
for(i in 1:nrow(A1_individuals)){
  temp<-A1_individuals[i,] 
  AdultsinA1<-dataParentsID[dataParentsID$aviary=="4",]
  temp$Sex<-AdultsinA1$sex[which(AdultsinA1$qr_code==temp$Individuals)]
  temp$Sex_Colour[which(temp$Sex=="f")]<-"red"
  temp$Sex_Colour[which(temp$Sex=="m")]<-"lightblue"
  A1_individuals[i,]<-temp
}
V(net)$color <- "grey"
V(net)$color <- A1_individuals$Sex_Colour[match(V(net)$name,as.character(A1_individuals$Individuals))]



##### Provisioning Network (Subject Adult -> Mother of recipient juvenile) #####
#### Prepare the dataframe for network 
### subset the dataframe, Event, to exclude the events in which the parents of recipient chicks had unbackpacked offspring. This is because we cannot identify how much exactly the subject adults fed to the offspring of the parents of recipient chicks

Aviary1_ProvisionNetwork_Mother_Event<-data.frame(Event=Event_All_Count_A1$EventID,AdultID=Event_All_Count_A1$AdultID,FedChickID=Event_All_Count_A1$Who,Categoty=Event_All_Count_A1$Category,NoNQPeriod_FatherOfPresentNonOffspring=as.character(Event_All_Count_A1$No_NQOffspring_Period_FatherOfPresentNonOffspring),NoNQPeriod_MotherOfPresentNonOffspring=as.character(Event_All_Count_A1$No_NQOffspring_Period_MotherOfPresentNonOffspring),NoNQPeriod_FathersOFRecipients=as.character(Event_All_Count_A1$No_NQOffspring_Period_FathersOFRecipients),NoNQPeriod_MotherOFRecipients=as.character(Event_All_Count_A1$No_NQOffspring_Period_MothersOFRecipients),Sum_Count=Event_All_Count_A1$Sum_Count,Day=Event_All_Count_A1$Day,MotherOfRecipientChick=Event_All_Count_A1$MotherOfRecipientChick,FatherOfRecipientChick=Event_All_Count_A1$FatherOfRecipientChick,FatherOfPresentNonOffspring=Event_All_Count_A1$FatherOfPresentNonOffspring,MotherOfPresentNonOffspring=Event_All_Count_A1$MotherOfPresentNonOffspring)


AllMale_NoNQPeriod<-NULL
AllFemale_NoNQPeriod<-NULL
AllAdult_NoNQPeriod<-NULL
Aviary1_ProvisionNetwork_Mother_Event$AllMale_NoNQPeriod<-NA
Aviary1_ProvisionNetwork_Mother_Event$AllFemale_NoNQPeriod<-NA
Aviary1_ProvisionNetwork_Mother_Event$AllAdult_NoNQPeriod<-NA
Aviary1_ProvisionNetwork_Mother_Event$NoNQPeriod_FatherOfPresentNonOffspring<-as.character(Aviary1_ProvisionNetwork_Mother_Event$NoNQPeriod_FatherOfPresentNonOffspring)
Aviary1_ProvisionNetwork_Mother_Event$NoNQPeriod_MotherOfPresentNonOffspring<-as.character(Aviary1_ProvisionNetwork_Mother_Event$NoNQPeriod_MotherOfPresentNonOffspring)
Aviary1_ProvisionNetwork_Mother_Event$NoNQPeriod_FathersOFRecipients<-as.character(Aviary1_ProvisionNetwork_Mother_Event$NoNQPeriod_FathersOFRecipients)
Aviary1_ProvisionNetwork_Mother_Event$NoNQPeriod_MotherOFRecipients<-as.character(Aviary1_ProvisionNetwork_Mother_Event$NoNQPeriod_MotherOFRecipients)

for(i in 1:nrow(Aviary1_ProvisionNetwork_Mother_Event)){
  temp<-Aviary1_ProvisionNetwork_Mother_Event[i,]
  AllMale_NoNQPeriod_temp<-NULL
  AllFemale_NoNQPeriod_temp<-NULL
  AllAdult_NoNQPeriod_temp<-NULL
  if(temp$NoNQPeriod_FatherOfPresentNonOffspring!=""){
    AllMale_NoNQPeriod_temp<-c(unlist(strsplit(temp$NoNQPeriod_FatherOfPresentNonOffspring," ")),unlist(strsplit(temp$NoNQPeriod_FathersOFRecipients, " ")))
    temp$AllMale_NoNQPeriod<-paste(AllMale_NoNQPeriod_temp,collapse = " ")
  }else if(temp$NoNQPeriod_FatherOfPresentNonOffspring==""){
    AllMale_NoNQPeriod_temp<-c(unlist(strsplit(temp$NoNQPeriod_FathersOFRecipients, " ")))
    temp$AllMale_NoNQPeriod<-paste(AllMale_NoNQPeriod_temp,collapse = " ")
  }
  if(temp$NoNQPeriod_MotherOfPresentNonOffspring!=""){
    AllFemale_NoNQPeriod_temp<-c(unlist(strsplit(temp$NoNQPeriod_MotherOfPresentNonOffspring," ")),unlist(strsplit(temp$NoNQPeriod_MotherOFRecipients, " ")))
    temp$AllFemale_NoNQPeriod<-paste(AllFemale_NoNQPeriod_temp,collapse = " ")
  }else if(temp$NoNQPeriod_FatherOfPresentNonOffspring==""){
    AllFemale_NoNQPeriod_temp<-c(unlist(strsplit(temp$NoNQPeriod_MotherOFRecipients, " ")))
    temp$AllFemale_NoNQPeriod<-paste(AllFemale_NoNQPeriod_temp,collapse = " ")
  }
  AllAdult_NoNQPeriod_temp<-c(AllMale_NoNQPeriod_temp,AllFemale_NoNQPeriod_temp)
  temp$AllAdult_NoNQPeriod<-paste(AllAdult_NoNQPeriod_temp,collapse = " ")
  Aviary1_ProvisionNetwork_Mother_Event[i,]<-temp
}
table(Aviary1_ProvisionNetwork_Mother_Event$AdultID==Aviary1_ProvisionNetwork_Mother_Event$FedChickID)#should be all FALSE



##### [Provision Network: Adults -> MotherOFrecipient] #####
library(dplyr)
library(tidyr)
###### Prepare the dataframe for network
### Summarize the fed events by unique combination of Adults-MotherOFrecipient
#Create the Adult-MotherOFrecipient combination column
Aviary1_ProvisionNetwork_Mother_Event <- Aviary1_ProvisionNetwork_Mother_Event%>%
  unite(AM_Combinations, AdultID, MotherOfRecipientChick, sep = "_", remove = FALSE)

#Weight the edge based on opportunity to feed
Table_num_fedADULTtoMOTHER<-table(Aviary1_ProvisionNetwork_Mother_Event$AM_Combinations)
AM_Combinations<-names(Table_num_fedADULTtoMOTHER)
Opportunity<-NULL
RecipientMother<-NULL
SubjectAdults<-unique(Aviary1_ProvisionNetwork_Mother_Event$AdultID)
RecipientAdults<-unique(Aviary1_ProvisionNetwork_Mother_Event$MotherOfRecipientChick)
AllAdults<-c(SubjectAdults,RecipientAdults)
AllAdults<-unique(AllAdults)
Aviary1_RecipientMother_ProvisionNetwork_sum<-data.frame(matrix(NA,0,3))
colnames(Aviary1_RecipientMother_ProvisionNetwork_sum)<-c("SubjectAdult","RecipientMother","EdgeWeight")
Opportunity<-NULL
for(i in 1:length(AM_Combinations)){
  tempAM<-AM_Combinations[i]
  Opportunity_temp<-NULL
  AllInfo<-Aviary1_ProvisionNetwork_Mother_Event[which(Aviary1_ProvisionNetwork_Mother_Event$AM_Combinations==tempAM),]
  All_Events<-unique(AllInfo$Event)
  
  SumCount<-NULL
  SumPresenceChickofMother<-NULL
  for(j in 1:length(All_Events)){
    SumCount_temp<-NULL
    SumPresenceChickofMother_temp<-NULL
    temp_All_Events<-All_Events[j]
    SumCount_temp<-unique(AllInfo$Sum_Count[which(AllInfo$Event==temp_All_Events)])
    SumCount<-c(SumCount,SumCount_temp)
    
    RecipientMom<-AllInfo$MotherOfRecipientChick[which(AllInfo$Event==temp_All_Events)]
    
    AnyMom<-unlist(strsplit(as.character(AllInfo$MotherOfPresentNonOffspring[which(AllInfo$Event==temp_All_Events)]), " "))
    
    SumPresenceChickofMother_temp<-sum(RecipientMom %in% AnyMom)
    SumPresenceChickofMother<-c(SumPresenceChickofMother,SumPresenceChickofMother_temp)
  }
  SumCount<-sum(SumCount)
  SumPresenceChickofMother<-sum(SumPresenceChickofMother)
  
  Opportunity_temp<-nrow(AllInfo)/(SumCount*SumPresenceChickofMother)
  Opportunity<-c(Opportunity,Opportunity_temp)
}
AM_Combinations
Opportunity
AM_EdgeWeight_Table<-data.frame(AM_Combinations=AM_Combinations,EdgeWeight=Opportunity)#Regard "Opportunity" as "Edgeweight"
A1ProvNet_MotherRec<-data.frame(AM_Combinations=AM_EdgeWeight_Table$AM_Combinations,SubjectAdult=NA,RecipientMother=NA,EdgeWeight=AM_EdgeWeight_Table$EdgeWeight)
colnames(A1ProvNet_MotherRec)<-c("AM_Combinations","SubjectAdult","RecipientMother","EdgeWeight")
for(i in 1:nrow(A1ProvNet_MotherRec)){
  A1ProvNet_MotherRec$SubjectAdult[i]<-unlist(strsplit(as.character(A1ProvNet_MotherRec$AM_Combinations[i]),"_"))[1]
}
for(i in 1:nrow(A1ProvNet_MotherRec)){
  A1ProvNet_MotherRec$RecipientMother[i]<-unlist(strsplit(as.character(A1ProvNet_MotherRec$AM_Combinations[i]),"_"))[2]
}
A1ProvNet_MotherRec<-A1ProvNet_MotherRec[,-1]
A1ProvNet_MotherRec<-data.frame(A1ProvNet_MotherRec)
A1ProvNet_MotherRec_2<-A1ProvNet_MotherRec
A1ProvNet_MotherRec<-A1ProvNet_MotherRec[!is.infinite(A1ProvNet_MotherRec_2$EdgeWeight),]
A1ProvNet_MotherRec_2$EdgeWeight[is.infinite(A1ProvNet_MotherRec_2$EdgeWeight)]<-NA #set edge weight NA for provisioning events towards offspring
# identify breeding pairs 
dataChickID_A1<-dataChickID[which(dataChickID$fosteraviary=="4"),]
Pairs_A1<-data.frame(RearingPapa=dataChickID_A1$RearingFather,RearingMama=dataChickID_A1$RearingMother)
Pairs_A1<-unique(Pairs_A1)





##### Plot network
##### http://www.shizukalab.com/toolkits/sna/weighted-edgelists
##### weight = Count
A1ProvNet_MotherRec<-as.data.frame(A1ProvNet_MotherRec) 
A1ProvNet_MotherRec$SubjectAdult<-as.character(A1ProvNet_MotherRec$SubjectAdult)
A1ProvNet_MotherRec$RecipientMother<-as.character(A1ProvNet_MotherRec$RecipientMother)#Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
A1ProvNet_MotherRec<-A1ProvNet_MotherRec[which(A1ProvNet_MotherRec$SubjectAdult!=A1ProvNet_MotherRec$RecipientMother),]
A1ProvNet_MotherRec=as.matrix(A1ProvNet_MotherRec) #igraph needs the edgelist to be in matrix format
A1ProvNet_MotherRec_network=graph.edgelist(A1ProvNet_MotherRec[,1:2]) #We first greate a network from the first two columns, which has the list of vertices
E(A1ProvNet_MotherRec_network)$Count=as.numeric(A1ProvNet_MotherRec[,3]) #We then add the edge weights to this network by assigning an edge attribute called 'weight'. 
### want to add colour
SubjectAdult<-unique(A1ProvNet_MotherRec_2$SubjectAdult)
RecipientMother<-unique(as.vector(A1ProvNet_MotherRec_2$RecipientMother))
A1_individuals<-data.frame(Individuals=c(SubjectAdult,RecipientMother), Sex=NA)
A1AllMales<-dataParentsID$qr_code[which(dataParentsID$sex=="m"&dataParentsID$aviary=="4")]
A1AllFemales<-dataParentsID$qr_code[which(dataParentsID$sex=="f"&dataParentsID$aviary=="4")]
A1_individuals$Sex[which(A1_individuals$Individuals%in%A1AllMales)]<-"Male"
A1_individuals$Sex[which(A1_individuals$Individuals%in%A1AllFemales)]<-"Female"
A1_individuals$Colour<-"grey" #defalt=grey
A1_individuals$Colour[which(A1_individuals$Sex=="Male")]<-"lightblue"
A1_individuals$Colour[which(A1_individuals$Sex=="Female")]<-"tomato"

V(A1ProvNet_MotherRec_network)$color <- "grey"
V(A1ProvNet_MotherRec_network)$color <- A1_individuals$Colour[match(V(A1ProvNet_MotherRec_network)$name,as.character(A1_individuals$Individuals))]
### plot 
pdf("A4ProvNet_MotherRec_network.pdf")
coords <- layout_with_fr(A1ProvNet_MotherRec_network)
plot(A1ProvNet_MotherRec_network,layout=coords,edge.width=(5*E(A1ProvNet_MotherRec_network)$Count^1.2))
dev.off()





### Prepare the provisioning network metrics
A1ProvNet_MotherRec_2<-as.data.frame(A1ProvNet_MotherRec_2) 
SubjectAdult<-unique(A1ProvNet_MotherRec_2$SubjectAdult)
RecipientMother<-unique(as.vector(A1ProvNet_MotherRec_2$RecipientMother))
A1ProvNet_MotherRec_2$SubjectAdult<-as.character(A1ProvNet_MotherRec_2$SubjectAdult)
A1ProvNet_MotherRec_2$RecipientMother<-as.character(A1ProvNet_MotherRec_2$RecipientMother)#Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
A1ProvNet_MotherRec_2=as.matrix(A1ProvNet_MotherRec_2) #igraph needs the edgelist to be in matrix format
A1ProvNet_MotherRec_network=graph.edgelist(A1ProvNet_MotherRec_2[,1:2]) #We first greate a network from the first two columns, which has the list of vertices
E(A1ProvNet_MotherRec_network)$Count=as.numeric(A1ProvNet_MotherRec_2[,3]) #We then add the edge weights to this network by assigning an edge attribute called 'weight'. 





##### Change dataframe to draw provisioning network into matrix
as_adjacency_matrix(A1ProvNet_MotherRec_network,sparse = F, attr = "Count") 
as.matrix(as_adjacency_matrix(A1ProvNet_MotherRec_network,sparse = F, attr = "Count"))

prov.female <- as.matrix(as_adjacency_matrix(A1ProvNet_MotherRec_network,sparse = F, attr = "Count"))
#put NA into the column of themselves->themselves
diag(prov.female)<-NA
#put NA into the column of males (recipient)
maleAdults<-dataParentsID$qr_code[which(dataParentsID$sex=="m"&dataParentsID$aviary=="4")]
maleAdults<-intersect(AllAdults,maleAdults)
maleAdults<-as.character(maleAdults)
prov.female[,maleAdults]<-NA
#put NA into the cells of breeding partners
femaleAdults<-dataParentsID$qr_code[which(dataParentsID$sex=="f"&dataParentsID$aviary=="4")]
femaleAdults<-intersect(AllAdults,femaleAdults)
femaleAdults<-as.character(femaleAdults)
for(i in 1:length(femaleAdults)){
  RMama<-femaleAdults[i]
  RPapa<-as.character(Pairs_A1$RearingPapa[which(Pairs_A1$RearingMama==RMama)])
  RPapa<-intersect(RPapa,maleAdults)
  if(length(RPapa)>0){
    for(j in 1:length(RPapa)){
      RPapa_temp<-RPapa[j]
      prov.female[RPapa_temp,RMama]<-NA
      prov.female[RMama,RPapa_temp]<-NA
    }
  }
}
#put NA into the cells of combination of adults who were never observed together
AllProv<-rownames(prov.female)
AllSubAd<-unique(ObservedCombination_FedNonFed_female$SubjectAdult)

for(k in 1:length(AllProv)){
  temp_AllProv<-AllProv[k]
  Sub_or_not<-length(intersect(temp_AllProv,AllSubAd))
  if(Sub_or_not>0){
    for(i in 1:length(AllSubAd)){
      SubAd<-as.character(AllSubAd[i])
      RecipientAd_temp<-as.character(ObservedCombination_FedNonFed_female$Present_RecipientAdult[which(ObservedCombination_FedNonFed_female$SubjectAdult==SubAd)])
      RecipientAd_NA_temp<-setdiff(AllProv,RecipientAd_temp)
      if(length(RecipientAd_NA_temp)>0){
        for(j in 1:length(RecipientAd_NA_temp)){
          RecipientAd_NA_temp2<-RecipientAd_NA_temp[j]
          prov.female[SubAd,RecipientAd_NA_temp2]<-NA
        }
      }
    }
  }else if(Sub_or_not==0){
    prov.female[temp_AllProv,]<-NA
  }
}







### Extract the individuals in the provision network, and delete other individuals from social network
## Match the order of individuals between provision network and social network
soc <- prebreedmtx[which(rownames(prebreedmtx) %in% rownames(prov.female)), which(colnames(prebreedmtx) %in% colnames(prov.female))]
soc <- soc[order(rownames(soc)),order(colnames(soc))]
########## Exclude breeding pairs from social network
# identify breeding pairs 
dataChickID_A1<-dataChickID[which(dataChickID$fosteraviary=="4"),]
Pairs_A1<-data.frame(RearingPapa=dataChickID_A1$RearingFather,RearingMama=dataChickID_A1$RearingMother)
Pairs_A1<-unique(Pairs_A1)
for(i in 1:nrow(Pairs_A1)){
  Papa<-Pairs_A1$RearingPapa[i]
  Mama<-Pairs_A1$RearingMama[i]
  temp_soc<-soc[rownames(soc)%in%Papa,colnames(soc)%in%Mama]
  temp_soc2<-soc[rownames(soc)%in%Mama,colnames(soc)%in%Papa]
  if(length(temp_soc)>0){
    soc[rownames(soc)%in%Papa,colnames(soc)%in%Mama]<-0
    soc[rownames(soc)%in%Mama,colnames(soc)%in%Papa]<-0
  }
}



##### Social Network (pre-breeding) excluding individuals which are not in provisioning network
# Plot the network(igraph)
net<-graph.adjacency(soc,mode="undirected", weighted=TRUE,diag=FALSE)
e <- get.edgelist(net, names=F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net), area=8*(vcount(net)^2), repulse.rad = vcount(net)^3.6)
individuals_temp<-rownames(soc)
#add colours to adults
A1_individuals<-data.frame(Individuals=individuals_temp, Sex=NA)
A1_individuals$Sex_Colour<-"green"
for(i in 1:nrow(A1_individuals)){
  temp<-A1_individuals[i,] 
  AdultsinA1<-dataParentsID[dataParentsID$aviary=="4",]
  temp$Sex<-AdultsinA1$sex[which(AdultsinA1$qr_code==temp$Individuals)]
  temp$Sex_Colour[which(temp$Sex=="f")]<-"tomato"
  temp$Sex_Colour[which(temp$Sex=="m")]<-"lightblue"
  A1_individuals[i,]<-temp
}
V(net)$color <- "grey"
V(net)$color <- A1_individuals$Sex_Colour[match(V(net)$name,as.character(A1_individuals$Individuals))]

pdf("Aviary4_PrebreedingSocialNetwork_RecipientMother.pdf")
coords <- layout_with_fr(net)
plot(net,layout=coords,edge.width=(200*E(net)$weight^1.2))
dev.off()










#### Adult-recipientFather ####

##### Load the data
Event_All_Count <- read.csv("Event_All_Count.csv",stringsAsFactors=FALSE, sep=",")

# chick info
dataChickID <- read.csv("dataChickID_toPublish20200223.csv",stringsAsFactors=FALSE,sep=",")

# parent info
dataParentsID <- read.csv("dataParentsID.csv",stringsAsFactors=FALSE,sep=",")


Event_All_Count_A1<-Event_All_Count[which(Event_All_Count$Aviary=="4"),]
Event_All_Count_A1 <- Event_All_Count_A1[which(Event_All_Count_A1$Category=="Unrelated"),]


## Which adults and offspring of adults (herein, recipient adults) were observed together?
SubjectAdult<-NULL
AllAdult<-unique(Event_All_Count_A1$AdultID)
Present_RecipientAdult_temp2<-NULL
for(i in 1:length(AllAdult)){
  temp_Event_All_Count_A1<-Event_All_Count_A1[which(Event_All_Count_A1$AdultID==AllAdult[i]),]
  SubjectAdult_temp<-NULL
  Present_RecipientAdult_temp<-unique(unlist(strsplit(temp_Event_All_Count_A1$FatherOfPresentNonOffspring," ")))
  length_old<-length(Present_RecipientAdult_temp2)
  for(j in 1:length(Present_RecipientAdult_temp)){
    temp_RecipientAdult<-Present_RecipientAdult_temp[j]
    if(AllAdult[i]!=temp_RecipientAdult){
      Present_RecipientAdult_temp2<-c(Present_RecipientAdult_temp2,temp_RecipientAdult)
    }
  }
  length_new<-length(Present_RecipientAdult_temp2)
  length_to_add<-length_new-length_old
  if(length(Present_RecipientAdult_temp2)>0){
    SubjectAdult_temp<-rep(AllAdult[i],length_to_add)
  }
  SubjectAdult<-c(SubjectAdult,SubjectAdult_temp)
}
length(SubjectAdult)==length(Present_RecipientAdult_temp2)
ObservedCombination_FedNonFed_male<-data.frame(SubjectAdult=SubjectAdult,Present_RecipientAdult=Present_RecipientAdult_temp2)


##When were the observations done?
ObservationDays_List<-unique(as.character(Event_All_Count_A1$Day))
##Which adults attended as which sides (e.g. giver or reciever)?
SubjectAdults<-unique(Event_All_Count_A1$AdultID)
RecipientAdults<-unique(Event_All_Count_A1$FatherOfRecipientChick)
AllAdults<-c(SubjectAdults,RecipientAdults)
AllAdults<-unique(AllAdults)







##### [Pre-breeding Social Network *include only adults in the provision network]
# Load the data
load("2018_A4_prebreedingMtx.RData") #prebreedmtx
# Change the column names from ring number to backpack
individuals_temp<-NULL
for(i in 1:ncol(prebreedmtx)){
  backpack_temp<-dataParentsID$qr_code[which(colnames(prebreedmtx)[i]==dataParentsID$ind_ID)]
  individuals_temp<-c(individuals_temp,backpack_temp)
}
rownames(prebreedmtx)<-individuals_temp
colnames(prebreedmtx)<-rownames(prebreedmtx)
# Plot the network(igraph)
net<-graph.adjacency(prebreedmtx,mode="undirected", weighted=TRUE,diag=FALSE)
e <- get.edgelist(net, names=F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net), area=8*(vcount(net)^2), repulse.rad = vcount(net)^3.6)
#add colours to adults
A1_individuals<-data.frame(Individuals=individuals_temp, Sex=NA)
A1_individuals$Sex_Colour<-"green"
for(i in 1:nrow(A1_individuals)){
  temp<-A1_individuals[i,] 
  AdultsinA1<-dataParentsID[dataParentsID$aviary=="4",]
  temp$Sex<-AdultsinA1$sex[which(AdultsinA1$qr_code==temp$Individuals)]
  temp$Sex_Colour[which(temp$Sex=="f")]<-"red"
  temp$Sex_Colour[which(temp$Sex=="m")]<-"lightblue"
  A1_individuals[i,]<-temp
}
V(net)$color <- "grey"
V(net)$color <- A1_individuals$Sex_Colour[match(V(net)$name,as.character(A1_individuals$Individuals))]



##### Provisioning Network (Subject Adult -> Mother of recipient juvenile) #####
#### Prepare the dataframe for network 
### subset the dataframe, Event, to exclude the events in which the parents of recipient chicks had unbackpacked offspring. This is because we cannot identify how much exactly the subject adults fed to the offspring of the parents of recipient chicks

Aviary1_ProvisionNetwork_Father_Event<-data.frame(Event=Event_All_Count_A1$EventID,AdultID=Event_All_Count_A1$AdultID,FedChickID=Event_All_Count_A1$Who,Categoty=Event_All_Count_A1$Category,NoNQPeriod_FatherOfPresentNonOffspring=as.character(Event_All_Count_A1$No_NQOffspring_Period_FatherOfPresentNonOffspring),NoNQPeriod_MotherOfPresentNonOffspring=as.character(Event_All_Count_A1$No_NQOffspring_Period_MotherOfPresentNonOffspring),NoNQPeriod_FathersOFRecipients=as.character(Event_All_Count_A1$No_NQOffspring_Period_FathersOFRecipients),NoNQPeriod_MotherOFRecipients=as.character(Event_All_Count_A1$No_NQOffspring_Period_MothersOFRecipients),Sum_Count=Event_All_Count_A1$Sum_Count,Day=Event_All_Count_A1$Day,MotherOfRecipientChick=Event_All_Count_A1$MotherOfRecipientChick,FatherOfRecipientChick=Event_All_Count_A1$FatherOfRecipientChick,FatherOfPresentNonOffspring=Event_All_Count_A1$FatherOfPresentNonOffspring,MotherOfPresentNonOffspring=Event_All_Count_A1$MotherOfPresentNonOffspring)


AllMale_NoNQPeriod<-NULL
AllFemale_NoNQPeriod<-NULL
AllAdult_NoNQPeriod<-NULL
Aviary1_ProvisionNetwork_Father_Event$AllMale_NoNQPeriod<-NA
Aviary1_ProvisionNetwork_Father_Event$AllFemale_NoNQPeriod<-NA
Aviary1_ProvisionNetwork_Father_Event$AllAdult_NoNQPeriod<-NA
Aviary1_ProvisionNetwork_Father_Event$NoNQPeriod_FatherOfPresentNonOffspring<-as.character(Aviary1_ProvisionNetwork_Father_Event$NoNQPeriod_FatherOfPresentNonOffspring)
Aviary1_ProvisionNetwork_Father_Event$NoNQPeriod_MotherOfPresentNonOffspring<-as.character(Aviary1_ProvisionNetwork_Father_Event$NoNQPeriod_MotherOfPresentNonOffspring)
Aviary1_ProvisionNetwork_Father_Event$NoNQPeriod_FathersOFRecipients<-as.character(Aviary1_ProvisionNetwork_Father_Event$NoNQPeriod_FathersOFRecipients)

for(i in 1:nrow(Aviary1_ProvisionNetwork_Father_Event)){
  temp<-Aviary1_ProvisionNetwork_Father_Event[i,]
  AllMale_NoNQPeriod_temp<-NULL
  AllFemale_NoNQPeriod_temp<-NULL
  AllAdult_NoNQPeriod_temp<-NULL
  if(temp$NoNQPeriod_FatherOfPresentNonOffspring!=""){
    AllMale_NoNQPeriod_temp<-c(unlist(strsplit(temp$NoNQPeriod_FatherOfPresentNonOffspring," ")),unlist(strsplit(temp$NoNQPeriod_FathersOFRecipients, " ")))
    temp$AllMale_NoNQPeriod<-paste(AllMale_NoNQPeriod_temp,collapse = " ")
  }else if(temp$NoNQPeriod_FatherOfPresentNonOffspring==""){
    AllMale_NoNQPeriod_temp<-c(unlist(strsplit(temp$NoNQPeriod_FathersOFRecipients, " ")))
    temp$AllMale_NoNQPeriod<-paste(AllMale_NoNQPeriod_temp,collapse = " ")
  }
  if(temp$NoNQPeriod_MotherOfPresentNonOffspring!=""){
    AllFemale_NoNQPeriod_temp<-c(unlist(strsplit(temp$NoNQPeriod_MotherOfPresentNonOffspring," ")),unlist(strsplit(as.character(temp$NoNQPeriod_MotherOFRecipients), " ")))
    temp$AllFemale_NoNQPeriod<-paste(AllFemale_NoNQPeriod_temp,collapse = " ")
  }else if(temp$NoNQPeriod_FatherOfPresentNonOffspring==""){
    AllFemale_NoNQPeriod_temp<-c(unlist(strsplit(as.character(temp$NoNQPeriod_MotherOFRecipients), " ")))
    temp$AllFemale_NoNQPeriod<-paste(AllFemale_NoNQPeriod_temp,collapse = " ")
  }
  AllAdult_NoNQPeriod_temp<-c(AllMale_NoNQPeriod_temp,AllFemale_NoNQPeriod_temp)
  temp$AllAdult_NoNQPeriod<-paste(AllAdult_NoNQPeriod_temp,collapse = " ")
  Aviary1_ProvisionNetwork_Father_Event[i,]<-temp
}
table(Aviary1_ProvisionNetwork_Father_Event$AdultID==Aviary1_ProvisionNetwork_Father_Event$FedChickID)#should be all FALSE



##### [Provision Network: Adults -> MotherOFrecipient] #####
library(dplyr)
library(tidyr)
###### Prepare the dataframe for network
### Summarize the fed events by unique combination of Adults-MotherOFrecipient
#Create the Adult-MotherOFrecipient combination column
Aviary1_ProvisionNetwork_Father_Event <- Aviary1_ProvisionNetwork_Father_Event%>%
  unite(AF_Combinations, AdultID, FatherOfRecipientChick, sep = "_", remove = FALSE)

#Weight the edge based on opportunity to feed
Table_num_fedADULTtoFATHER<-table(Aviary1_ProvisionNetwork_Father_Event$AF_Combinations)
AF_Combinations<-names(Table_num_fedADULTtoFATHER)
Opportunity<-NULL
RecipientFahter<-NULL
SubjectAdults<-unique(Aviary1_ProvisionNetwork_Father_Event$AdultID)
RecipientAdults<-unique(Aviary1_ProvisionNetwork_Father_Event$FatherOfRecipientChick)
AllAdults<-c(SubjectAdults,RecipientAdults)
AllAdults<-unique(AllAdults)

Aviary1_RecipientFather_ProvisionNetwork_sum<-data.frame(matrix(NA,0,3))
colnames(Aviary1_RecipientFather_ProvisionNetwork_sum)<-c("SubjectAdult","RecipientFather","EdgeWeight")
Opportunity<-NULL
for(i in 1:length(AF_Combinations)){
  tempAF<-AF_Combinations[i]
  Opportunity_temp<-NULL
  AllInfo<-Aviary1_ProvisionNetwork_Father_Event[which(Aviary1_ProvisionNetwork_Father_Event$AF_Combinations==tempAF),]
  All_Events<-unique(AllInfo$Event)
  
  SumCount<-NULL
  SumPresenceChickofFather<-NULL
  for(j in 1:length(All_Events)){
    SumCount_temp<-NULL
    SumPresenceChickofFather_temp<-NULL
    temp_All_Events<-All_Events[j]
    SumCount_temp<-unique(AllInfo$Sum_Count[which(AllInfo$Event==temp_All_Events)])
    SumCount<-c(SumCount,SumCount_temp)
    
    RecipientDad<-AllInfo$FatherOfRecipientChick[which(AllInfo$Event==temp_All_Events)]
    
    AnyDad<-unlist(strsplit(as.character(AllInfo$FatherOfPresentNonOffspring[which(AllInfo$Event==temp_All_Events)]), " "))
    
    SumPresenceChickofFather_temp<-sum(RecipientDad %in% AnyDad)
    SumPresenceChickofFather<-c(SumPresenceChickofFather,SumPresenceChickofFather_temp)
  }
  SumCount<-sum(SumCount)
  SumPresenceChickofFather<-sum(SumPresenceChickofFather)
  
  Opportunity_temp<-nrow(AllInfo)/(SumCount*SumPresenceChickofFather)
  Opportunity<-c(Opportunity,Opportunity_temp)
}
AF_Combinations
Opportunity
AF_EdgeWeight_Table<-data.frame(AF_Combinations=AF_Combinations,EdgeWeight=Opportunity)#Regard "Opportunity" as "Edgeweight"
A1ProvNet_FatherRec<-data.frame(AF_Combinations=AF_EdgeWeight_Table$AF_Combinations,SubjectAdult=NA,RecipientFather=NA,EdgeWeight=AF_EdgeWeight_Table$EdgeWeight)
colnames(A1ProvNet_FatherRec)<-c("AF_Combinations","SubjectAdult","RecipientFather","EdgeWeight")
for(i in 1:nrow(A1ProvNet_FatherRec)){
  A1ProvNet_FatherRec$SubjectAdult[i]<-unlist(strsplit(as.character(A1ProvNet_FatherRec$AF_Combinations[i]),"_"))[1]
}
for(i in 1:nrow(A1ProvNet_FatherRec)){
  A1ProvNet_FatherRec$RecipientFather[i]<-unlist(strsplit(as.character(A1ProvNet_FatherRec$AF_Combinations[i]),"_"))[2]
}
A1ProvNet_FatherRec<-A1ProvNet_FatherRec[,-1]
A1ProvNet_FatherRec<-data.frame(A1ProvNet_FatherRec)
A1ProvNet_FatherRec_2<-A1ProvNet_FatherRec
A1ProvNet_FatherRec<-A1ProvNet_FatherRec[!is.infinite(A1ProvNet_FatherRec_2$EdgeWeight),]
A1ProvNet_FatherRec_2$EdgeWeight[is.infinite(A1ProvNet_FatherRec_2$EdgeWeight)]<-NA #set edge weight NA for provisioning events towards offspring
## need to set edge weight NA for provisioning events towards offspring of the adults who were not able to fulfill criteria on the observation days
#set edge weight NA for provisioning events towards offspring
# identify breeding pairs 
dataChickID_A1<-dataChickID[which(dataChickID$fosteraviary=="4"),]
Pairs_A1<-data.frame(RearingPapa=dataChickID_A1$RearingFather,RearingMama=dataChickID_A1$RearingMother)
Pairs_A1<-unique(Pairs_A1)



##### Plot network
##### http://www.shizukalab.com/toolkits/sna/weighted-edgelists
##### weight = Count
A1ProvNet_FatherRec<-as.data.frame(A1ProvNet_FatherRec) 
A1ProvNet_FatherRec$SubjectAdult<-as.character(A1ProvNet_FatherRec$SubjectAdult)
A1ProvNet_FatherRec$RecipientFather<-as.character(A1ProvNet_FatherRec$RecipientFather)#Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
A1ProvNet_FatherRec<-A1ProvNet_FatherRec[which(A1ProvNet_FatherRec$SubjectAdult!=A1ProvNet_FatherRec$RecipientFather),]
A1ProvNet_FatherRec=as.matrix(A1ProvNet_FatherRec) #igraph needs the edgelist to be in matrix format
A1ProvNet_FatherRec_network=graph.edgelist(A1ProvNet_FatherRec[,1:2]) #We first greate a network from the first two columns, which has the list of vertices
E(A1ProvNet_FatherRec_network)$Count=as.numeric(A1ProvNet_FatherRec[,3]) #We then add the edge weights to this network by assigning an edge attribute called 'weight'. 
### want to add colour
SubjectAdult<-unique(A1ProvNet_FatherRec_2$SubjectAdult)
RecipientFather<-unique(as.vector(A1ProvNet_FatherRec_2$RecipientFather))
A1_individuals<-data.frame(Individuals=c(SubjectAdult,RecipientFather), Sex=NA)
A1AllMales<-dataParentsID$qr_code[which(dataParentsID$sex=="m"&dataParentsID$aviary=="4")]
A1AllFemales<-dataParentsID$qr_code[which(dataParentsID$sex=="f"&dataParentsID$aviary=="4")]
A1_individuals$Sex[which(A1_individuals$Individuals%in%A1AllMales)]<-"Male"
A1_individuals$Sex[which(A1_individuals$Individuals%in%A1AllFemales)]<-"Female"
A1_individuals$Colour<-"grey" #defalt=grey
A1_individuals$Colour[which(A1_individuals$Sex=="Male")]<-"lightblue"
A1_individuals$Colour[which(A1_individuals$Sex=="Female")]<-"tomato"

V(A1ProvNet_FatherRec_network)$color <- "grey"
V(A1ProvNet_FatherRec_network)$color <- A1_individuals$Colour[match(V(A1ProvNet_FatherRec_network)$name,as.character(A1_individuals$Individuals))]
### plot 
pdf("A4ProvNet_FatherRec_network.pdf")
coords <- layout_with_fr(A1ProvNet_FatherRec_network)
plot(A1ProvNet_FatherRec_network,layout=coords,edge.width=(5*E(A1ProvNet_FatherRec_network)$Count^1.2))
dev.off()





### Prepare the provisioning network metrics
A1ProvNet_FatherRec_2<-as.data.frame(A1ProvNet_FatherRec_2) 
SubjectAdult<-unique(A1ProvNet_FatherRec_2$SubjectAdult)
RecipientFather<-unique(as.vector(A1ProvNet_FatherRec_2$RecipientFather))
A1ProvNet_FatherRec_2$SubjectAdult<-as.character(A1ProvNet_FatherRec_2$SubjectAdult)
A1ProvNet_FatherRec_2$RecipientFather<-as.character(A1ProvNet_FatherRec_2$RecipientFather)#Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
A1ProvNet_FatherRec_2=as.matrix(A1ProvNet_FatherRec_2) #igraph needs the edgelist to be in matrix format
A1ProvNet_FatherRec_network=graph.edgelist(A1ProvNet_FatherRec_2[,1:2]) #We first greate a network from the first two columns, which has the list of vertices
E(A1ProvNet_FatherRec_network)$Count=as.numeric(A1ProvNet_FatherRec_2[,3]) #We then add the edge weights to this network by assigning an edge attribute called 'weight'. 





##### Change dataframe to draw provisioning network into matrix
as_adjacency_matrix(A1ProvNet_FatherRec_network,sparse = F, attr = "Count") 
as.matrix(as_adjacency_matrix(A1ProvNet_FatherRec_network,sparse = F, attr = "Count"))

prov.male <- as.matrix(as_adjacency_matrix(A1ProvNet_FatherRec_network,sparse = F, attr = "Count"))
#put NA into the column of themselves->themselves
diag(prov.male)<-NA
#put NA into the column of females (recipient)
femaleAdults<-dataParentsID$qr_code[which(dataParentsID$sex=="f"&dataParentsID$aviary=="4")]
femaleAdults<-intersect(AllAdults,femaleAdults)
femaleAdults<-as.character(femaleAdults)
prov.male[,femaleAdults]<-NA
#put NA into the cells of breeding partners
maleAdults<-dataParentsID$qr_code[which(dataParentsID$sex=="m"&dataParentsID$aviary=="4")]
maleAdults<-intersect(AllAdults,maleAdults)
maleAdults<-as.character(maleAdults)
for(i in 1:length(maleAdults)){
  RPapa<-maleAdults[i]
  RMama<-as.character(Pairs_A1$RearingMama[which(Pairs_A1$RearingPapa==RPapa)])
  RMama<-intersect(RMama,femaleAdults)
  if(length(RMama)>0){
    for(j in 1:length(RMama)){
      RMama_temp<-RMama[j]
      prov.male[RMama_temp,RPapa]<-NA
      prov.male[RPapa,RMama_temp]<-NA
    }
  }
}
#put NA into the cells of combination of adults who were never observed together
AllProv<-rownames(prov.male)
AllSubAd<-unique(ObservedCombination_FedNonFed_male$SubjectAdult)

for(k in 1:length(AllProv)){
  temp_AllProv<-AllProv[k]
  Sub_or_not<-length(intersect(temp_AllProv,AllSubAd))
  if(Sub_or_not>0){
    for(i in 1:length(AllSubAd)){
      SubAd<-as.character(AllSubAd[i])
      RecipientAd_temp<-as.character(ObservedCombination_FedNonFed_male$Present_RecipientAdult[which(ObservedCombination_FedNonFed_male$SubjectAdult==SubAd)])
      RecipientAd_NA_temp<-setdiff(AllProv,RecipientAd_temp)
      if(length(RecipientAd_NA_temp)>0){
        for(j in 1:length(RecipientAd_NA_temp)){
          RecipientAd_NA_temp2<-RecipientAd_NA_temp[j]
          prov.male[SubAd,RecipientAd_NA_temp2]<-NA
        }
      }
    }
  }else if(Sub_or_not==0){
    prov.male[temp_AllProv,]<-NA
  }
}



### Extract the individuals in the provision network, and delete other individuals from social network
## Match the order of individuals between provision network and social network
soc <- prebreedmtx[which(rownames(prebreedmtx) %in% rownames(prov.male)), which(colnames(prebreedmtx) %in% colnames(prov.male))]
soc <- soc[order(rownames(soc)),order(colnames(soc))]
########## Exclude breeding pairs from social network
# identify breeding pairs 
dataChickID_A1<-dataChickID[which(dataChickID$fosteraviary=="4"),]
Pairs_A1<-data.frame(RearingPapa=dataChickID_A1$RearingFather,RearingMama=dataChickID_A1$RearingMother)
Pairs_A1<-unique(Pairs_A1)
for(i in 1:nrow(Pairs_A1)){
  Papa<-Pairs_A1$RearingPapa[i]
  Mama<-Pairs_A1$RearingMama[i]
  temp_soc<-soc[rownames(soc)%in%Papa,colnames(soc)%in%Mama]
  temp_soc2<-soc[rownames(soc)%in%Mama,colnames(soc)%in%Papa]
  if(length(temp_soc)>0){
    soc[rownames(soc)%in%Papa,colnames(soc)%in%Mama]<-0
    soc[rownames(soc)%in%Mama,colnames(soc)%in%Papa]<-0
  }
}



##### Social Network (pre-breeding) excluding individuals which are not in provisioning network
# Plot the network(igraph)
net<-graph.adjacency(soc,mode="undirected", weighted=TRUE,diag=FALSE)
e <- get.edgelist(net, names=F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net), area=8*(vcount(net)^2), repulse.rad = vcount(net)^3.6)
individuals_temp<-rownames(soc)
#add colours to adults
A1_individuals<-data.frame(Individuals=individuals_temp, Sex=NA)
A1_individuals$Sex_Colour<-"green"
for(i in 1:nrow(A1_individuals)){
  temp<-A1_individuals[i,] 
  AdultsinA1<-dataParentsID[dataParentsID$aviary=="4",]
  temp$Sex<-AdultsinA1$sex[which(AdultsinA1$qr_code==temp$Individuals)]
  temp$Sex_Colour[which(temp$Sex=="f")]<-"tomato"
  temp$Sex_Colour[which(temp$Sex=="m")]<-"lightblue"
  A1_individuals[i,]<-temp
}
V(net)$color <- "grey"
V(net)$color <- A1_individuals$Sex_Colour[match(V(net)$name,as.character(A1_individuals$Individuals))]

pdf("Aviary4_PrebreedingSocialNetwork_RecipientFather.pdf")
coords <- layout_with_fr(net)
plot(net,layout=coords,edge.width=(200*E(net)$weight^1.2))
dev.off()





##########
##### Prepare network matrics for Mantel test
### need to sum up the prov.male and prov.female
inds<-unique(c(rownames(prov.male),rownames(prov.female)))#prov.femaleに含まれてなくてprov.maleにある
inds <- inds[order(inds)]
prov <- matrix(0, nrow=length(inds), ncol=length(inds))
colnames(prov) <- rownames(prov)  <- inds
prov.female[is.na(prov.female)]<-0
prov.male[is.na(prov.male)]<-0
prov[rownames(prov.female),rownames(prov.female)] <- prov[rownames(prov.female),rownames(prov.female)] + prov.female
prov[rownames(prov.male),rownames(prov.male)] <- prov[rownames(prov.male),rownames(prov.male)] + prov.male
#put NA into the column of themselves->themselves
diag(prov)<-NA
#put NA into the cells of breeding partners
maleAdults<-dataParentsID$qr_code[which(dataParentsID$sex=="m"&dataParentsID$aviary=="4")]
maleAdults<-intersect(AllAdults,maleAdults)
maleAdults<-as.character(maleAdults)
for(i in 1:length(maleAdults)){
  RPapa<-maleAdults[i]
  RMama<-as.character(Pairs_A1$RearingMama[which(Pairs_A1$RearingPapa==RPapa)])
  RMama<-intersect(RMama,femaleAdults)
  if(length(RMama)>0){
    for(j in 1:length(RMama)){
      RMama_temp<-RMama[j]
      prov[RMama_temp,RPapa]<-NA
      prov[RPapa,RMama_temp]<-NA
    }
  }
}
#put NA into the cells of combination of adults who were never observed together
ObservedCombination_FedNonFed<-rbind(ObservedCombination_FedNonFed_female,ObservedCombination_FedNonFed_male)
ObservedCombination_FedNonFed<-unique(ObservedCombination_FedNonFed)

AllProv<-rownames(prov)
AllSubAd<-unique(ObservedCombination_FedNonFed$SubjectAdult)

for(k in 1:length(AllProv)){
  temp_AllProv<-AllProv[k]
  Sub_or_not<-length(intersect(temp_AllProv,AllSubAd))
  if(Sub_or_not>0){
    for(i in 1:length(AllSubAd)){
      SubAd<-as.character(AllSubAd[i])
      RecipientAd_temp<-as.character(ObservedCombination_FedNonFed$Present_RecipientAdult[which(ObservedCombination_FedNonFed$SubjectAdult==SubAd)])
      RecipientAd_NA_temp<-setdiff(AllProv,RecipientAd_temp)
      if(length(RecipientAd_NA_temp)>0){
        for(j in 1:length(RecipientAd_NA_temp)){
          RecipientAd_NA_temp2<-RecipientAd_NA_temp[j]
          prov[SubAd,RecipientAd_NA_temp2]<-NA
        }
      }
    }
  }else if(Sub_or_not==0){
    prov[temp_AllProv,]<-NA
  }
}

soc <- prebreedmtx[rownames(prebreedmtx) %in% inds, colnames(prebreedmtx) %in% inds]
soc <- soc[order(rownames(soc)),order(colnames(soc))]



library(vegan)
mantel(prov, soc, na.rm = T) #r=-0.1551, Significance 0.89



AllAdults<-rownames(prov)
femaleAdults<-dataParentsID$qr_code[which(dataParentsID$sex=="f"&dataParentsID$aviary=="4")]
femaleAdults<-intersect(AllAdults,femaleAdults)
femaleAdults<-as.character(femaleAdults)
maleAdults<-dataParentsID$qr_code[which(dataParentsID$sex=="m"&dataParentsID$aviary=="4")]
maleAdults<-intersect(AllAdults,maleAdults)
maleAdults<-as.character(maleAdults)
length(AllAdults)
length(femaleAdults)
length(maleAdults)



##### provisioning network
prov[is.na(prov)]<-0
to_delete<-NULL#delete individuals whose columns and rows (both) are all 0.
for(i in 1:nrow(prov)){
  testprov_temp<-rownames(prov)[i]
  to_delete_temp<-NULL
  if(sum(prov[testprov_temp,]==0)==nrow(prov)){
    if(sum(prov[,testprov_temp]==0)==nrow(prov)){
      to_delete_temp<-testprov_temp
      to_delete<-c(to_delete,to_delete_temp)
    }
  }
}
to_keep<-setdiff(rownames(prov),to_delete)
prov<-prov[to_keep,to_keep]


net2<-graph.adjacency(prov,mode="directed", weighted=TRUE,diag=FALSE)
e <- get.edgelist(net2, names=F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net2), area=8*(vcount(net2)^2), repulse.rad = vcount(net2)^3.6)
individuals_temp<-rownames(prov)
#add colours to adults
inds<-data.frame(Individuals=individuals_temp, Sex=NA)
inds$Sex_Colour<-"green"
for(i in 1:nrow(inds)){
  temp<-inds[i,] 
  AdultsinA1<-dataParentsID[dataParentsID$aviary=="4",]
  temp$Sex<-AdultsinA1$sex[which(AdultsinA1$qr_code==temp$Individuals)]
  temp$Sex_Colour[which(temp$Sex=="f")]<-"tomato"
  temp$Sex_Colour[which(temp$Sex=="m")]<-"lightblue"
  inds[i,]<-temp
}
V(net2)$color <- "grey"
V(net2)$color <- inds$Sex_Colour[match(V(net2)$name,as.character(inds$Individuals))]

pdf("Paper_Aviary4_ProvisioningNetwork.pdf")
coords <- layout_with_fr(net2)
plot(net2,layout=coords, edge.width=(10*E(net2)$weight^1.2))
dev.off()



##### Social Network (pre-breeding) excluding individuals which are not in provisioning network
soc<-soc[which(rownames(soc) %in% rownames(prov)), which(colnames(soc) %in% colnames(prov))]
soc <- soc[order(rownames(soc)),order(colnames(soc))]

rownames(soc)==rownames(prov)#must be all TRUE

# Plot the network(igraph)
net<-graph.adjacency(soc,mode="undirected", weighted=TRUE,diag=FALSE)
e <- get.edgelist(net, names=F)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net), area=8*(vcount(net)^2), repulse.rad = vcount(net)^3.6)
individuals_temp<-rownames(soc)
#add colours to adults
A1_individuals<-data.frame(Individuals=individuals_temp, Sex=NA)
A1_individuals$Sex_Colour<-"green"
for(i in 1:nrow(A1_individuals)){
  temp<-A1_individuals[i,] 
  AdultsinA1<-dataParentsID[dataParentsID$aviary=="4",]
  temp$Sex<-AdultsinA1$sex[which(AdultsinA1$qr_code==temp$Individuals)]
  temp$Sex_Colour[which(temp$Sex=="f")]<-"tomato"
  temp$Sex_Colour[which(temp$Sex=="m")]<-"lightblue"
  A1_individuals[i,]<-temp
}
V(net)$color <- "grey"
V(net)$color <- A1_individuals$Sex_Colour[match(V(net)$name,as.character(A1_individuals$Individuals))]

pdf("Paper_Aviary4_PrebreedingSocialNetwork.pdf")
plot(net,layout=coords, edge.width=(200*E(net)$weight^1.2))
dev.off()

