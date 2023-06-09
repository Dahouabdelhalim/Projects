#Create networks, dendrograms, calculate social groups and social communities

#Use this script to make networks, dendrograms, calculate group and community membership, export files used 
#in other scripts. Also calculate sexes and classes per group and community. 

#In this code, I use Resident, Auxilliary, and 1 year old terms. 
#In paper changed these terms but still refer to same individuals
#Resident = Previous breeder
#Auxilliary = Helper (other)
#1 Year old = hatched previous breeding season (other)


####Load Data####
library(here)

#bird17 has an R because it includes birds in removal communities as part of an experiment performed that year
#Including them in the initial network to calculate association indices but then removing them before calculating
#group membership

#Load in bird dataframes
bird16 = read.csv(here::here("Input files","bird16.csv"))
bird17R = read.csv(here::here("Input files","bird17R.csv"))
bird18 = read.csv(here::here("Input files","bird18.csv"))
bird19 = read.csv(here::here("Input files","bird19.csv"))

#Load in agesexstatus files
agesexstatus16 = read.csv(here::here("Input files","agesexstatus16.csv"))
agesexstatus17 = read.csv(here::here("Input files","agesexstatus17.csv"))
agesexstatus18 = read.csv(here::here("Input files","agesexstatus18.csv"))
agesexstatus19 = read.csv(here::here("Input files","agesexstatus19.csv"))

#Get dataframes without kerfuffle interactions (within 8 min of courtship or aggression)
bird16noK = bird16[which(bird16$Kerfuffle=="N"),]
bird17noK = bird17R[which(bird17R$Kerfuffle=="N"),]
bird18noK = bird18[which(bird18$Kerfuffle=="N"),]
bird19noK = bird19[which(bird19$Kerfuffle=="N"),]

#Get individuals list for network
individuals16noK = data.frame(bird16noK$Bird, bird16noK$Sighting)
colnames(individuals16noK) = c("ID","Sighting")
individuals17noK = data.frame(bird17noK$Bird, bird17noK$Sighting)
colnames(individuals17noK) = c("ID","Sighting")
individuals18noK = data.frame(bird18noK$Bird, bird18noK$Sighting)
colnames(individuals18noK) = c("ID","Sighting")
individuals19noK = data.frame(bird19noK$Bird, bird19noK$Sighting)
colnames(individuals19noK) = c("ID","Sighting")


#Get Group by Individual Matrix 
library(asnipe)
gbi16noK = get_group_by_individual(individuals16noK, data_format = "individuals")
gbi17noK = get_group_by_individual(individuals17noK, data_format = "individuals")
gbi18noK = get_group_by_individual(individuals18noK, data_format = "individuals")
gbi19noK = get_group_by_individual(individuals19noK, data_format = "individuals")


##Get filtered network so I can calculate degree
#2016
network16noK = get_network(gbi16noK, data_format= "GBI", association_index = "SRI") 
network16noK = network16noK[order(rownames(network16noK)), order(colnames(network16noK))]
birdorder16noK = data.frame(rownames(network16noK))
colnames(birdorder16noK)[1] = "Bird"
birdlist16noK = birdorder16noK

#2017
network17noK = get_network(gbi17noK, data_format= "GBI", association_index = "SRI") 
network17noK = network17noK[order(rownames(network17noK)), order(colnames(network17noK))]
birdorder17noK = data.frame(rownames(network17noK))
colnames(birdorder17noK)[1] = "Bird"
birdlist17noK = birdorder17noK

#2018
network18noK = get_network(gbi18noK, data_format= "GBI", association_index = "SRI") 
network18noK = network18noK[order(rownames(network18noK)), order(colnames(network18noK))]
birdorder18noK = data.frame(rownames(network18noK))
colnames(birdorder18noK)[1] = "Bird"
birdlist18noK = birdorder18noK

#2019
network19noK = get_network(gbi19noK, data_format= "GBI", association_index = "SRI") 
network19noK = network19noK[order(rownames(network19noK)), order(colnames(network19noK))]
birdorder19noK = data.frame(rownames(network19noK))
colnames(birdorder19noK)[1] = "Bird"
birdlist19noK = birdorder19noK

####Sight Freq for filtered network####
#2016
sightfreqnoK16 = data.frame(colSums(gbi16noK))
sightfreqnoK16$Bird = rownames(sightfreqnoK16)
colnames(sightfreqnoK16)[1] = "SightFreq"
birdlist16noK = merge(birdlist16noK,sightfreqnoK16,by="Bird")
library(igraph)
net16noK = graph.adjacency(network16noK, mode="undirected", diag=FALSE, weighted=TRUE)
deg16noK = degree(net16noK)
birdlist16noK = data.frame(birdlist16noK,deg16noK)

plot(birdlist16noK$SightFreq,birdlist16noK$deg16noK,pch=19, xlab = "Sighting Frequency",
     ylab = "Degree")
cor.test(birdlist16noK$SightFreq,birdlist16noK$deg16noK)
birdlist16noKsub = birdlist16noK[which(birdlist16noK$SightFreq>32),]
plot(birdlist16noKsub$SightFreq,birdlist16noKsub$deg16noK,pch=19,xlab = "Sighting Frequency",
     ylab = "Degree")
cor.test(birdlist16noKsub$SightFreq,birdlist16noKsub$deg16noK)


#2017
sightfreqnoK17 = data.frame(colSums(gbi17noK))
sightfreqnoK17$Bird = rownames(sightfreqnoK17)
colnames(sightfreqnoK17)[1] = "SightFreq"
birdlist17noK = merge(birdlist17noK,sightfreqnoK17,by="Bird")
library(igraph)
net17noK = graph.adjacency(network17noK, mode="undirected", diag=FALSE, weighted=TRUE)
deg17noK = degree(net17noK)
birdlist17noK = data.frame(birdlist17noK,deg17noK)

plot(birdlist17noK$SightFreq,birdlist17noK$deg17noK,pch=19, xlab = "Sighting Frequency",
     ylab = "Degree")
cor.test(birdlist17noK$SightFreq,birdlist17noK$deg17noK)
birdlist17noKsub = birdlist17noK[which(birdlist17noK$SightFreq>54),] #54 for obs17 at 8min noK
plot(birdlist17noKsub$SightFreq,birdlist17noKsub$deg17noK,pch=19,xlab = "Sighting Frequency",
     ylab = "Degree")
cor.test(birdlist17noKsub$SightFreq,birdlist17noKsub$deg17noK)



#2018
sightfreqnoK18 = data.frame(colSums(gbi18noK))
sightfreqnoK18$Bird = rownames(sightfreqnoK18)
colnames(sightfreqnoK18)[1] = "SightFreq"
birdlist18noK = merge(birdlist18noK,sightfreqnoK18,by="Bird")
library(igraph)
net18noK = graph.adjacency(network18noK, mode="undirected", diag=FALSE, weighted=TRUE)
deg18noK = degree(net18noK)
birdlist18noK = data.frame(birdlist18noK,deg18noK)

plot(birdlist18noK$SightFreq,birdlist18noK$deg18noK,pch=19, xlab = "Sighting Frequency",
     ylab = "Degree")
cor.test(birdlist18noK$SightFreq,birdlist18noK$deg18noK)
birdlist18noKsub = birdlist18noK[which(birdlist18noK$SightFreq>9),]
plot(birdlist18noKsub$SightFreq,birdlist18noKsub$deg18noK,pch=19,xlab = "Sighting Frequency",
     ylab = "Degree")
cor.test(birdlist18noKsub$SightFreq,birdlist18noKsub$deg18noK)


#2019
sightfreqnoK19 = data.frame(colSums(gbi19noK))
sightfreqnoK19$Bird = rownames(sightfreqnoK19)
colnames(sightfreqnoK19)[1] = "SightFreq"
birdlist19noK = merge(birdlist19noK,sightfreqnoK19,by="Bird")
library(igraph)
net19noK = graph.adjacency(network19noK, mode="undirected", diag=FALSE, weighted=TRUE)
deg19noK = degree(net19noK)
birdlist19noK = data.frame(birdlist19noK,deg19noK)

plot(birdlist19noK$SightFreq,birdlist19noK$deg19noK,pch=19, xlab = "Sighting Frequency",
     ylab = "Degree")
cor.test(birdlist19noK$SightFreq,birdlist19noK$deg19noK)
birdlist19noKsub = birdlist19noK[which(birdlist19noK$SightFreq>16),]
plot(birdlist19noKsub$SightFreq,birdlist19noKsub$deg19noK,pch=19,xlab = "Sighting Frequency",
     ylab = "Degree")
cor.test(birdlist19noKsub$SightFreq,birdlist19noKsub$deg19noK)

#Put all on same graph
par(mfrow=c(2,2))
par(xpd=FALSE)
plot(birdlist16noK$SightFreq,birdlist16noK$deg16noK,pch=19, xlab = "Sighting Frequency",
     ylab = "Degree", main="2016", lwd=0.5)
abline(v=32,lwd=3, lty=2, col="red")
plot(birdlist17noK$SightFreq,birdlist17noK$deg17noK,pch=19, xlab = "Sighting Frequency",
     ylab = "Degree", main="2017", lwd=0.5)
abline(v=54,lwd=3, lty=2, col="red")
plot(birdlist18noK$SightFreq,birdlist18noK$deg18noK,pch=19, xlab = "Sighting Frequency",
     ylab = "Degree", main="2018", lwd=0.5)
abline(v=9,lwd=3, lty=2, col="red")
plot(birdlist19noK$SightFreq,birdlist19noK$deg19noK,pch=19, xlab = "Sighting Frequency",
     ylab = "Degree", main="2019", lwd=0.5)
abline(v=17,lwd=3, lty=2, col="red")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(birdlist16noKsub$SightFreq,birdlist16noKsub$deg16noK,pch=19,xlab = "Sighting Frequency",
     ylab = "Degree")
plot(birdlist17noKsub$SightFreq,birdlist17noKsub$deg17noK,pch=19,xlab = "Sighting Frequency",
     ylab = "Degree")
plot(birdlist18noKsub$SightFreq,birdlist18noKsub$deg18noK,pch=19,xlab = "Sighting Frequency",
     ylab = "Degree")
plot(birdlist19noKsub$SightFreq,birdlist19noKsub$deg19noK,pch=19,xlab = "Sighting Frequency",
     ylab = "Degree")
par(mfrow=c(1,1))




####Distance Matrix for filtered network####
#Get distances based on association indices
network16noKA = sweep(network16noK,1,1) #Subtract one from each association index to get association distance
network16noKB = sweep(network16noKA,1,-1,FUN = "*") #Make it positive 
network17noKA = sweep(network17noK,1,1)
network17noKB = sweep(network17noKA,1,-1,FUN = "*")
network18noKA = sweep(network18noK,1,1)
network18noKB = sweep(network18noKA,1,-1,FUN = "*")
network19noKA = sweep(network19noK,1,1)
network19noKB = sweep(network19noKA,1,-1,FUN = "*")


#Take out removal birds before plotting dendrogram
RorC = read.csv(here::here("Input files","RemovalorControl 2017.csv"), header=TRUE)
removalbirds = RorC[which(RorC$TreatmentGroup=="Removal"),]
birdlist17noKsub = birdlist17noKsub[!birdlist17noKsub$Bird %in% removalbirds$Bird,]


#Subset for sight freq - this is after network has been calculated using all birds
#But only assigning groups to birds seen over XX times
#This works because birdlist16noKsub still has all of it's levels in birdlist16noKsub$Bird
network16noKB = network16noKB[birdlist16noKsub$Bird,birdlist16noKsub$Bird]
network17noKB = network17noKB[birdlist17noKsub$Bird,birdlist17noKsub$Bird]
network18noKB = network18noKB[birdlist18noKsub$Bird,birdlist18noKsub$Bird]
network19noKB = network19noKB[birdlist19noKsub$Bird,birdlist19noKsub$Bird]

#Get distance matrix
#Average below plots the average of a pairs connection to the next bird as the distance of 
#the branches. UPGMA method.
d16noK = as.dist(network16noKB)
fit16noK <- hclust(d16noK, method="average") 
d17noK = as.dist(network17noKB)
fit17noK <- hclust(d17noK, method="average") 
d18noK = as.dist(network18noKB)
fit18noK <- hclust(d18noK, method="average") 
d19noK = as.dist(network19noKB)
fit19noK <- hclust(d19noK, method="average") 

#Plot basic dendrograms
plot(fit16noK) # display dendogram
plot(fit17noK)
plot(fit18noK)
plot(fit19noK)


library(dendextend)

fit16noKd = as.dendrogram(fit16noK)
dlabels16noK = data.frame(labels(fit16noKd))
fit17noKd = as.dendrogram(fit17noK)
dlabels17noK = data.frame(labels(fit17noKd))
fit18noKd = as.dendrogram(fit18noK)
dlabels18noK = data.frame(labels(fit18noKd))
fit19noKd = as.dendrogram(fit19noK)
dlabels19noK = data.frame(labels(fit19noKd))

#colors = data.frame(c("Resident Male","Resident Female","1 Year-old Male","Dispersing Female",
#"Unknown","Auxilliary Male"),c("#D55E00","#4DA1D1","#D55E00","#4DA1D1","#999999","#D55E00"))
colors = data.frame(c("M","F","U"),c("#D55E00","#4DA1D1","#999999"))
colnames(colors)=c("Sex","colors")
colors16 = merge(colors,agesexstatus16,by="Sex",all.y=T)
colors17 = merge(colors,agesexstatus17,by="Sex",all.y=T)
colors18 = merge(colors,agesexstatus18,by="Sex",all.y=T)
colors19 = merge(colors,agesexstatus19,by="Sex",all.y=T)

shapes = data.frame(c("Resident Male","Resident Female","1 Year-old Male","Dispersing Female",
                      "Unknown","Auxilliary Male"),c(19,19,15,15,18,15))
colnames(shapes)=c("class","shape")
colors16 = merge(colors16,shapes,by="class",all.x=T)
colors17 = merge(colors17,shapes,by="class",all.x=T)
colors18 = merge(colors18,shapes,by="class",all.x=T)
colors19 = merge(colors19,shapes,by="class",all.x=T)

#Not sorting below keeps the dataframe in the same order as the dendrogram so colors get assigned correctly
dlabels16noK = merge(dlabels16noK,colors16,by.x="labels.fit16noKd.",by.y="Bird", sort=F, all.x = T) #don't sort
dlabels17noK = merge(dlabels17noK,colors17,by.x="labels.fit17noKd.",by.y="Bird", sort=F) #don't sort
dlabels18noK = merge(dlabels18noK,colors18,by.x="labels.fit18noKd.",by.y="Bird", sort=F) #don't sort
dlabels19noK = merge(dlabels19noK,colors19,by.x="labels.fit19noKd.",by.y="Bird", sort=F) #don't sort


#Get list of # groups by height
groupdatanoK16 = data.frame(seq(1,0,by=-0.01))
colnames(groupdatanoK16) = "Height"
groupdatanoK16$NumGroups = NA
for (i in 1:nrow(groupdatanoK16)) {
  groups = cutree(fit16noK,h=groupdatanoK16$Height[i])
  grouplist = data.frame(table(groups))
  groupdatanoK16$NumGroups[i] = nrow(grouplist)
}
groupdatanoK16 = groupdatanoK16[!duplicated(groupdatanoK16$NumGroups),]

groupdatanoK17 = data.frame(seq(1,0,by=-0.01))
colnames(groupdatanoK17) = "Height"
groupdatanoK17$NumGroups = NA
for (i in 1:nrow(groupdatanoK17)) {
  groups = cutree(fit17noK,h=groupdatanoK17$Height[i])
  grouplist = data.frame(table(groups))
  groupdatanoK17$NumGroups[i] = nrow(grouplist)
}
groupdatanoK17 = groupdatanoK17[!duplicated(groupdatanoK17$NumGroups),]


groupdatanoK18 = data.frame(seq(1,0,by=-0.01))
colnames(groupdatanoK18) = "Height"
groupdatanoK18$NumGroups = NA
for (i in 1:nrow(groupdatanoK18)) {
  groups = cutree(fit18noK,h=groupdatanoK18$Height[i])
  grouplist = data.frame(table(groups))
  groupdatanoK18$NumGroups[i] = nrow(grouplist)
}
groupdatanoK18 = groupdatanoK18[!duplicated(groupdatanoK18$NumGroups),]


groupdatanoK19 = data.frame(seq(1,0,by=-0.01))
colnames(groupdatanoK19) = "Height"
groupdatanoK19$NumGroups = NA
for (i in 1:nrow(groupdatanoK19)) {
  groups = cutree(fit19noK,h=groupdatanoK19$Height[i])
  grouplist = data.frame(table(groups))
  groupdatanoK19$NumGroups[i] = nrow(grouplist)
}
groupdatanoK19 = groupdatanoK19[!duplicated(groupdatanoK19$NumGroups),]



#Figure out where to cut each dendrogram using silhouette distance
library(fpc)
#2016
groupdatanoK16$average.within = NA
groupdatanoK16$avg.silwidth = NA
#Get an error for first and last cutpoints, so that's why 2: and nrow -1
for (i in 2:(nrow(groupdatanoK16)-1)) { 
  groups = data.frame(cutree(fit16noK, h=groupdatanoK16$Height[i]))
  colnames(groups) = "groupID"
  groupdatanoK16$average.within[i] = cluster.stats(network16noKB,groups$groupID)$average.within
  groupdatanoK16$avg.silwidth[i] = cluster.stats(network16noKB,groups$groupID)$avg.silwidth
}

plot(groupdatanoK16$Height,groupdatanoK16$average.within, pch=19)
lines(groupdatanoK16$Height,groupdatanoK16$average.within)
plot(groupdatanoK16$Height,groupdatanoK16$avg.silwidth, pch=19, xlab="Height",
     ylab = "Average Silhouette Distance")
lines(groupdatanoK16$Height,groupdatanoK16$avg.silwidth)


#2017
groupdatanoK17$average.within = NA
groupdatanoK17$avg.silwidth = NA
#Get an error for first and last cutpoints, so that's why 2: and nrow -1
for (i in 2:(nrow(groupdatanoK17)-1)) { 
  groups = data.frame(cutree(fit17noK, h=groupdatanoK17$Height[i]))
  colnames(groups) = "groupID"
  groupdatanoK17$average.within[i] = cluster.stats(network17noKB,groups$groupID)$average.within
  groupdatanoK17$avg.silwidth[i] = cluster.stats(network17noKB,groups$groupID)$avg.silwidth
}

plot(groupdatanoK17$Height,groupdatanoK17$average.within, pch=19)
lines(groupdatanoK17$Height,groupdatanoK17$average.within)
plot(groupdatanoK17$Height,groupdatanoK17$avg.silwidth, pch=19)
lines(groupdatanoK17$Height,groupdatanoK17$avg.silwidth)


#2018
groupdatanoK18$average.within = NA
groupdatanoK18$avg.silwidth = NA
#Get an error for first and last cutpoints, so that's why 2: and nrow -1
for (i in 2:(nrow(groupdatanoK18)-1)) { 
  groups = data.frame(cutree(fit18noK, h=groupdatanoK18$Height[i]))
  colnames(groups) = "groupID"
  groupdatanoK18$average.within[i] = cluster.stats(network18noKB,groups$groupID)$average.within
  groupdatanoK18$avg.silwidth[i] = cluster.stats(network18noKB,groups$groupID)$avg.silwidth
}

plot(groupdatanoK18$Height,groupdatanoK18$average.within, pch=19)
lines(groupdatanoK18$Height,groupdatanoK18$average.within)
plot(groupdatanoK18$Height,groupdatanoK18$avg.silwidth, pch=19)
lines(groupdatanoK18$Height,groupdatanoK18$avg.silwidth)


#2019
groupdatanoK19$average.within = NA
groupdatanoK19$avg.silwidth = NA
#Get an error for first and last cutpoints, so that's why 2: and nrow -1
for (i in 2:(nrow(groupdatanoK19)-1)) { 
  groups = data.frame(cutree(fit19noK, h=groupdatanoK19$Height[i]))
  colnames(groups) = "groupID"
  groupdatanoK19$average.within[i] = cluster.stats(network19noKB,groups$groupID)$average.within
  groupdatanoK19$avg.silwidth[i] = cluster.stats(network19noKB,groups$groupID)$avg.silwidth
}

plot(groupdatanoK19$Height,groupdatanoK19$average.within, pch=19)
lines(groupdatanoK19$Height,groupdatanoK19$average.within)
plot(groupdatanoK19$Height,groupdatanoK19$avg.silwidth, pch=19)
lines(groupdatanoK19$Height,groupdatanoK19$avg.silwidth)


#Get cut heights for dendrograms
groupdatanoK16order = groupdatanoK16[order(groupdatanoK16$avg.silwidth, decreasing =T),]
h16 = groupdatanoK16order$Height[1]
avg.s16 = groupdatanoK16order$avg.silwidth[1]

groupdatanoK17order = groupdatanoK17[order(groupdatanoK17$avg.silwidth, decreasing =T),]
h17 = groupdatanoK17order$Height[1]
avg.s17 = groupdatanoK17order$avg.silwidth[1]

groupdatanoK18order = groupdatanoK18[order(groupdatanoK18$avg.silwidth, decreasing =T),]
h18 = groupdatanoK18order$Height[1]
avg.s18 = groupdatanoK18order$avg.silwidth[1]

groupdatanoK19order = groupdatanoK19[order(groupdatanoK19$avg.silwidth, decreasing =T),]
h19 = groupdatanoK19order$Height[1]
avg.s19 = groupdatanoK19order$avg.silwidth[1]


#Plot all silhouette plots together: 
par(mfrow=c(2,2))
plot(groupdatanoK16$Height,groupdatanoK16$avg.silwidth, pch=19, xlab="Association Distance (1-SRI)",
     ylab = "Mean Silhouette Width",main="2016")
lines(groupdatanoK16$Height,groupdatanoK16$avg.silwidth)
abline(v=h16,lty=2,lwd=3,col="red")
plot(groupdatanoK17$Height,groupdatanoK17$avg.silwidth, pch=19, xlab="Association Distance (1-SRI)",
     ylab = "Mean Silhouette Width",main="2017")
lines(groupdatanoK17$Height,groupdatanoK17$avg.silwidth)
abline(v=h17,lty=2,lwd=3,col="red")
plot(groupdatanoK18$Height,groupdatanoK18$avg.silwidth, pch=19, xlab="Association Distance (1-SRI)",
     ylab = "Mean Silhouette Width",main="2018")
lines(groupdatanoK18$Height,groupdatanoK18$avg.silwidth)
abline(v=h18,lty=2,lwd=3,col="red")
plot(groupdatanoK19$Height,groupdatanoK19$avg.silwidth, pch=19, xlab="Association Distance (1-SRI)",
     ylab = "Mean Silhouette Width",main="2019")
lines(groupdatanoK19$Height,groupdatanoK19$avg.silwidth)
abline(v=h19,lty=2,lwd=3,col="red")
par(mfrow=c(1,1))



####Plot dendrograms noK####
fit16noKL = as.dendrogram(fit16noK,hang=0.07) %>% 
  #set("branches_k_color", value = c("gray")) %>% 
  #color_branches(k=1,col="grey") %>%
  set("branches_lwd",2.5) %>%
  #set("labels_colors",h=0.9) %>% 
  set("leaves_pch", dlabels16noK$shape) %>% 
  set("leaves_cex",1.5) %>%
  set("leaves_col", as.character(dlabels16noK$colors))
plot(fit16noKL,ylab="Association Distance (1-SRI)",leaflab = "none")#, main="2016 Season")
abline(h=h16,lty=2, lwd=2.5,col="red")


fit17noKL = as.dendrogram(fit17noK,hang=0.07) %>% 
  #set("branches_k_color", value = c("skyblue", "orange", "grey"), h=0.9) %>% 
  #color_branches(k=1,col="grey") %>%
  set("branches_lwd",2.5) %>%
  #set("labels_colors",h=0.9) %>% 
  set("leaves_pch", dlabels17noK$shape) %>% 
  set("leaves_cex",1.5) %>%
  set("leaves_col", as.character(dlabels17noK$colors))
plot(fit17noKL,ylab="Association Distance (1-SRI)",leaflab="none")#,main="2017 Season")
abline(h=h17,lty=2,lwd=2.5, col="red")

fit18noKL = as.dendrogram(fit18noK,hang=0.07) %>% 
  #set("branches_k_color", value = c("skyblue", "orange", "grey"), h=0.9) %>% 
  #color_branches(k=1,col="grey") %>%
  set("branches_lwd",2.5) %>%
  #set("labels_colors",h=0.9) %>% 
  set("leaves_pch", dlabels18noK$shape) %>% 
  set("leaves_cex",1.5) %>%
  set("leaves_col", as.character(dlabels18noK$colors))
plot(fit18noKL,ylab="Association Distance (1-SRI)",leaflab="none")#, main="2018 Season")#,leaflab = "none")
abline(h=h18,lty=2,lwd=2.5,col="red")

fit19noKL = as.dendrogram(fit19noK,hang=0.07) %>% 
  #set("branches_k_color", value = c("skyblue", "orange", "grey"), h=0.9) %>% 
  #color_branches(k=1,col="grey") %>%
  set("branches_lwd",2.5) %>%
  #set("labels_colors",h=0.9) %>% 
  set("leaves_pch", dlabels19noK$shape) %>% 
  set("leaves_cex",1.5) %>%
  set("leaves_col", as.character(dlabels19noK$colors))
plot(fit19noKL,ylab="Association Distance (1-SRI)",leaflab="none")#,main="2019 Season")#,leaflab = "none")
abline(h=h19,lty=2,lwd=2.5,col="red")




#__________________________________________________________________________#




####*Group Structure Analyses####
#Who's in each group?
library(tibble)
#2016
group16 = data.frame(cutree(fit16noK, h=h16))
group16 = rownames_to_column(group16)
colnames(group16)=c("Bird","Social.Group")

#Average group size for 2016
group16t = data.frame(table(group16$Social.Group))
colnames(group16t)=c("Group","Freq")
#remove groups that only have one bird bc that means we didn't get good group info for those individuals
group16t = group16t[which(group16t$Freq!=1),] 
groupsize16 = mean(group16t$Freq)

#2017
group17 = data.frame(cutree(fit17noK, h=h17))
group17 = rownames_to_column(group17)
colnames(group17)=c("Bird","Social.Group")

#Average group size for 2017
group17t = data.frame(table(group17$Social.Group))
colnames(group17t)=c("Group","Freq")
#remove groups that only have one bird bc that means we didn't get good group info for those individuals
group17t = group17t[which(group17t$Freq!=1),] 
groupsize17 = mean(group17t$Freq)

#2018
group18 = data.frame(cutree(fit18noK, h=h18))
group18 = rownames_to_column(group18)
colnames(group18)=c("Bird","Social.Group")

#Average group size for 2018
group18t = data.frame(table(group18$Social.Group))
colnames(group18t)=c("Group","Freq")
#remove groups that only have one bird bc that means we didn't get good group info for those individuals
group18t = group18t[which(group18t$Freq!=1),] 
groupsize18 = mean(group18t$Freq)

#2019
group19 = data.frame(cutree(fit19noK, h=h19))
group19 = rownames_to_column(group19)
colnames(group19)=c("Bird","Social.Group")

#Average group size for 2019
group19t = data.frame(table(group19$Social.Group))
colnames(group19t)=c("Group","Freq")
#remove groups that only have one bird bc that means we didn't get good group info for those individuals
group19t = group19t[which(group19t$Freq!=1),] 
groupsize19 = mean(group19t$Freq)



####Social Classes per Group####
library(cowplot)
library(dplyr)
library(readr)

setwd("/Users/Joe/Documents/Research/R Packages/RainCloudPlots-master/tutorial_R")
source("R_rainclouds.R")

#Put group identity into birdlist
birdlist16noKsubg = merge(birdlist16noKsub,group16,by="Bird")
birdlist17noKsubg = merge(birdlist17noKsub,group17,by="Bird")
birdlist18noKsubg = merge(birdlist18noKsub,group18,by="Bird")
birdlist19noKsubg = merge(birdlist19noKsub,group19,by="Bird")

#Bring in age, sex, and status - has previous year's breeding status and this year's non-breeding status
birdlist16noKsubg = merge(birdlist16noKsubg,agesexstatus16,by="Bird",all.x=T)
birdlist17noKsubg = merge(birdlist17noKsubg,agesexstatus17,by="Bird",all.x=T)
birdlist18noKsubg = merge(birdlist18noKsubg,agesexstatus18,by="Bird",all.x=T)
birdlist19noKsubg = merge(birdlist19noKsubg,agesexstatus19,by="Bird",all.x=T)

#Remove single individual groups
birdlist16noKsubgu = birdlist16noKsubg[birdlist16noKsubg$Social.Group %in% group16t$Group,]
birdlist17noKsubgu = birdlist17noKsubg[birdlist17noKsubg$Social.Group %in% group17t$Group,]
birdlist18noKsubgu = birdlist18noKsubg[birdlist18noKsubg$Social.Group %in% group18t$Group,]
birdlist19noKsubgu = birdlist19noKsubg[birdlist19noKsubg$Social.Group %in% group19t$Group,]

#Get number of social classes per group
library(plyr)
library(dplyr)
library(tidyr)

groupclass16 = count(birdlist16noKsubgu,c("Social.Group","class"))
groupclass16 = complete(groupclass16,Social.Group,class,fill=list(freq=0))
groupclass16$Year = "2016"
groupclass17 = count(birdlist17noKsubgu,c("Social.Group","class"))
groupclass17 = complete(groupclass17,Social.Group,class,fill=list(freq=0))
groupclass17$Year = "2017"
groupclass18 = count(birdlist18noKsubgu,c("Social.Group","class"))
groupclass18 = complete(groupclass18,Social.Group,class,fill=list(freq=0))
groupclass18$Year = "2018"
groupclass19 = count(birdlist19noKsubgu,c("Social.Group","class"))
groupclass19 = complete(groupclass19,Social.Group,class,fill=list(freq=0))
groupclass19$Year = "2019"


#Combine years
groupclassall = rbind(groupclass16,groupclass17,groupclass18,groupclass19)

###Resident Males
#Get Resident males
groupclassallrm = groupclassall[which(groupclassall$class=="Resident Male"),]

#Get mean values
groupclassallrmsum = data.frame(Year=c("2016","2017","2018","2019"),
                                Mean=c(mean(groupclass16[which(groupclass16$class=="Resident Male"),]$freq),
                                       mean(groupclass17[which(groupclass17$class=="Resident Male"),]$freq),
                                       mean(groupclass18[which(groupclass18$class=="Resident Male"),]$freq),
                                       mean(groupclass19[which(groupclass19$class=="Resident Male"),]$freq)))

#Plot number of resident males per group
library(ggplot2)
A = ggplot(data=groupclassallrm,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.2,h=0),size=2)  +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  ylab("Individuals per Group") + ylim(-0.1,4.1) + geom_boxplot(fill=NA) +
  geom_point(position=position_jitter(width=0.2,h=0),size=2) +
  ggtitle("Resident Males") + geom_point(data=groupclassallrmsum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5)) 

###Resident Females
#Get Resident females
groupclassallrf = groupclassall[which(groupclassall$class=="Resident Female"),]

#Get mean values
groupclassallrfsum = data.frame(Year=c("2016","2017","2018","2019"),
                                Mean=c(mean(groupclass16[which(groupclass16$class=="Resident Female"),]$freq),
                                       mean(groupclass17[which(groupclass17$class=="Resident Female"),]$freq),
                                       mean(groupclass18[which(groupclass18$class=="Resident Female"),]$freq),
                                       mean(groupclass19[which(groupclass19$class=="Resident Female"),]$freq)))

#Plot number of resident females per group
B = ggplot(data=groupclassallrf,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.2,h=0),size=2)  +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  ylab("Individuals per Group") + ylim(-0.1,4.1) + geom_boxplot(fill=NA) +
  geom_point(position=position_jitter(width=0.2,h=0),size=2) +
  ggtitle("Resident Females") + geom_point(data=groupclassallrfsum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5)) 

###Auxilliary Males
groupclassallax = groupclassall[which(groupclassall$class=="Auxilliary Male"),]

#Get mean values
groupclassallaxsum = data.frame(Year=c("2016","2017","2018","2019"),
                                Mean=c(mean(groupclass16[which(groupclass16$class=="Auxilliary Male"),]$freq),
                                       mean(groupclass17[which(groupclass17$class=="Auxilliary Male"),]$freq),
                                       mean(groupclass18[which(groupclass18$class=="Auxilliary Male"),]$freq),
                                       mean(groupclass19[which(groupclass19$class=="Auxilliary Male"),]$freq)))

#Plot number of Auxilliary males per group
C = ggplot(data=groupclassallax,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.2,h=0),size=2)  +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  ylab("Individuals per Group") + ylim(-0.1,4.1) + geom_boxplot(fill=NA) +
  geom_point(position=position_jitter(width=0.2,h=0),size=2) +
  ggtitle("Auxilliary Males") + geom_point(data=groupclassallaxsum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5)) 

###1 Year Old Males
#Get One Year Old Males
groupclassall1m = groupclassall[which(groupclassall$class=="1 Year-old Male"),]

#Get mean values
groupclassall1msum = data.frame(Year=c("2016","2017","2018","2019"),
                                Mean=c(mean(groupclass16[which(groupclass16$class=="1 Year-old Male"),]$freq),
                                       mean(groupclass17[which(groupclass17$class=="1 Year-old Male"),]$freq),
                                       mean(groupclass18[which(groupclass18$class=="1 Year-old Male"),]$freq),
                                       mean(groupclass19[which(groupclass19$class=="1 Year-old Male"),]$freq)))

#Plot number of one year old males per group
D = ggplot(data=groupclassall1m,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.2,h=0),size=2)  +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  ylab("Individuals per Group") + ylim(-0.1,4.1) + geom_boxplot(fill=NA) +
  geom_point(position=position_jitter(width=0.2,h=0),size=2) +
  ggtitle("One-year Old Males") + geom_point(data=groupclassall1msum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5))  

###Dispersing Females
#Get Dispersing Females
groupclassalldf = groupclassall[which(groupclassall$class=="Dispersing Female"),]

#Get mean values
groupclassalldfsum = data.frame(Year=c("2016","2017","2018","2019"),
                                Mean=c(mean(groupclass16[which(groupclass16$class=="Dispersing Female"),]$freq),
                                       mean(groupclass17[which(groupclass17$class=="Dispersing Female"),]$freq),
                                       mean(groupclass18[which(groupclass18$class=="Dispersing Female"),]$freq),
                                       mean(groupclass19[which(groupclass19$class=="Dispersing Female"),]$freq)))

#Plot number of Dispersing Females per group
E = ggplot(data=groupclassalldf,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.2,h=0),size=2)  +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  ylab("Individuals per Group") + ylim(-0.1,4.1) + geom_boxplot(fill=NA) +
  geom_point(position=position_jitter(width=0.2,h=0),size=2) +
  ggtitle("Dispersing Females") + geom_point(data=groupclassalldfsum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5))  


###Unknowns
groupclassallu = groupclassall[which(groupclassall$class=="Unknown"),]

#Get mean values
groupclassallusum = data.frame(Year=c("2016","2017","2018","2019"),
                                Mean=c(mean(groupclass16[which(groupclass16$class=="Unknown"),]$freq),
                                       mean(groupclass17[which(groupclass17$class=="Unknown"),]$freq),
                                       mean(groupclass18[which(groupclass18$class=="Unknown"),]$freq),
                                       mean(groupclass19[which(groupclass19$class=="Unknown"),]$freq)))

#Plot number of Unknowns per group
G = ggplot(data=groupclassallu,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.2,h=0),size=2)  +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  ylab("Individuals per Group") + ylim(-0.1,4.1) + geom_boxplot(fill=NA) +
  geom_point(position=position_jitter(width=0.2,h=0),size=2) +
  ggtitle("Unknowns") + geom_point(data=groupclassallusum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5))  



##Plot all six social classes
#Don't need a legend here because x axis specifies categories.
library(ggpubr)
ggarrange(A,B,C,D,E,G) 


#Average group size across all years
groupall = rbind(group16t,group17t,group18t,group19t)
mean(groupall$Freq)


####Social Classes per group - all years in one plot####
groupclassallsum = data.frame(class=c("Resident Male","Resident Female","Auxilliary Male",
                                      "1 Year-old Male","Dispersing Female","Unknown"),
                              Mean=c(mean(groupclassall[which(groupclassall$class=="Resident Male"),]$freq),
                                     mean(groupclassall[which(groupclassall$class=="Resident Female"),]$freq),
                                     mean(groupclassall[which(groupclassall$class=="Auxilliary Male"),]$freq),
                                     mean(groupclassall[which(groupclassall$class=="1 Year-old Male"),]$freq),
                                     mean(groupclassall[which(groupclassall$class=="Dispersing Female"),]$freq),
                                     mean(groupclassall[which(groupclassall$class=="Unknown"),]$freq)))

#Order factors
groupclassall$class = factor(groupclassall$class,levels=c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male",
                                                          "Dispersing Female","Unknown"))

#Get multi-level labels for each factor
groupclassall$class = as.factor(groupclassall$class)
levels(groupclassall$class) <- gsub(" ", "\\n", levels(groupclassall$class))
levels(groupclassallsum$class) <- gsub(" ", "\\n", levels(groupclassallsum$class))

ggplot(data=groupclassall,aes(x=class,y=freq)) + 
  geom_flat_violin(aes(group=class),position = position_nudge(x=0.27,y=0),trim=F) +
  #geom_point(position=position_jitter(width=0.2,h=0.07),size=1.5, alpha=0.2) +
  geom_boxplot(data=groupclassall,aes(x=class,y=freq),colour="black",width=0.3) +
  #geom_point(data=groupclassallsum,aes(x=class,y=Mean),colour="red",size=2) +
  #geom_errorbar(data=groupclassallsum,aes(x=class,y=Mean,ymin=Mean-Stderrormean,ymax=Mean+Stderrormean),
                #width=0.15,colour="red") + 
  ylab("Individuals per Group") + ylim(-0.7,5) + xlab("Social Class")  +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
  panel.grid.major.y = element_line( size=.2, color="gray" )) # explicitly set the horizontal lines

ggplot(data=groupclassall,aes(x=class,y=freq,colour=Year)) + 
  geom_boxplot(data=groupclassall,aes(x=class,y=freq),colour="black",width=0.3) +
  geom_flat_violin(aes(group=class),position = position_nudge(x=0.27,y=0),trim=F) +
  geom_point(position=position_jitter(width=0.2,h=0.07),size=1.5) +
  geom_boxplot(data=groupclassall,aes(x=class,y=freq),colour="black",width=0.3,fill="transparent") +
  ylab("Individuals per Group")  + ylim(-0.7,5) + 
  xlab("Social Class") + scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
    panel.grid.major.y = element_line( size=.2, color="gray" )) # explicitly set the horizontal lines




#__________________________________________________________________________#



####*Communities####

#Include all sightings, including sightings associated with kerfuffle behaviors in community detection. 
#But only calculate communities for birds that were assigned to a group and were not in a group on their own
#First have to calculate the network with all sightings. 


#Get individuals list for network
individuals16 = data.frame(bird16$Bird, bird16$Sighting)
colnames(individuals16) = c("ID","Sighting")
individuals17 = data.frame(bird17R$Bird, bird17R$Sighting)
colnames(individuals17) = c("ID","Sighting")
individuals18 = data.frame(bird18$Bird, bird18$Sighting)
colnames(individuals18) = c("ID","Sighting")
individuals19 = data.frame(bird19$Bird, bird19$Sighting)
colnames(individuals19) = c("ID","Sighting")


#Get Group by Individual Matrix 
library(asnipe)
gbi16 = get_group_by_individual(individuals16, data_format = "individuals")
gbi17 = get_group_by_individual(individuals17, data_format = "individuals")
gbi18 = get_group_by_individual(individuals18, data_format = "individuals")
gbi19 = get_group_by_individual(individuals19, data_format = "individuals")


##Get network
#2016
network16 = get_network(gbi16, data_format= "GBI", association_index = "SRI") 
network16 = network16[order(rownames(network16)), order(colnames(network16))]
birdorder16 = data.frame(rownames(network16))
colnames(birdorder16)[1] = "Bird"

#2017
network17 = get_network(gbi17, data_format= "GBI", association_index = "SRI") 
network17 = network17[order(rownames(network17)), order(colnames(network17))]
birdorder17 = data.frame(rownames(network17))
colnames(birdorder17)[1] = "Bird"

#2018
network18 = get_network(gbi18, data_format= "GBI", association_index = "SRI") 
network18 = network18[order(rownames(network18)), order(colnames(network18))]
birdorder18 = data.frame(rownames(network18))
colnames(birdorder18)[1] = "Bird"

#2019
network19 = get_network(gbi19, data_format= "GBI", association_index = "SRI") 
network19 = network19[order(rownames(network19)), order(colnames(network19))]
birdorder19 = data.frame(rownames(network19))
colnames(birdorder19)[1] = "Bird"


#Match group ID to birdlist and take out individuals in a group on their own - we didn't get good 
#group/social data on them - groupt has been filtered to remove single bird groups above
birdlist16noKsubc = merge(birdlist16noKsub,group16,by="Bird")
birdlist17noKsubc = merge(birdlist17noKsub,group17,by="Bird")
birdlist18noKsubc = merge(birdlist18noKsub,group18,by="Bird")
birdlist19noKsubc = merge(birdlist19noKsub,group19,by="Bird")
birdlist16noKsubc = birdlist16noKsubc[birdlist16noKsubc$Social.Group %in% group16t$Group,]
birdlist17noKsubc = birdlist17noKsubc[birdlist17noKsubc$Social.Group %in% group17t$Group,]
birdlist18noKsubc = birdlist18noKsubc[birdlist18noKsubc$Social.Group %in% group18t$Group,]
birdlist19noKsubc = birdlist19noKsubc[birdlist19noKsubc$Social.Group %in% group19t$Group,]

#Match list of birds from network with all sightings to list of birds from subsetted network. 
#Usually don't need this step but different number of birds (different levels) 
#in the network and the birdlist dataframe - came from different datasets. In order to subset
#the kerfuffle network correctly I need the same number of levels the network has
birdorder16c = birdorder16[birdorder16$Bird %in% birdlist16noKsubc$Bird,]
birdorder17c = birdorder17[birdorder17$Bird %in% birdlist17noKsubc$Bird,]
birdorder18c = birdorder18[birdorder18$Bird %in% birdlist18noKsubc$Bird,]
birdorder19c = birdorder19[birdorder19$Bird %in% birdlist19noKsubc$Bird,]

#Subset network to match birdlistnoKsubc
network16subc = network16[birdorder16c,birdorder16c]
network17subc = network17[birdorder17c,birdorder17c]
network18subc = network18[birdorder18c,birdorder18c]
network19subc = network19[birdorder19c,birdorder19c]


####Communities of Groups####

#Question: How  do groups interact to form larger social communities? 
#Use size reduction of complex networks method from Arenas et al. 2007 paper to create a reduced network.
#The reduced network G in which each of these groups is replaced by a single node may be easily 
#defined in the following way: the weight W`rs between the nodes which represent
#groups r and s is the sum of all the weights connecting vertices in these groups. Internal connectivity of 
#nodes within groups is left in, meaning the overall strength of the reduced network is the same as the 
#original network. Leaving internal connections in vs taking them out does not influence community assignment
#using modularity on the reduced network

#In this reduced network, each group becomes a node, then we see how groups are structured into larger
#social communities
library(reshape2)

#Melt association matrix and match up group identities to each bird
rnetmelt16 = melt(network16subc)
rnetmelt17 = melt(network17subc)
rnetmelt18 = melt(network18subc)
rnetmelt19 = melt(network19subc)

#Match each variable (bird) to its group
rnetmelt16 = merge(rnetmelt16,group16,by.x="Var1",by.y="Bird")
rnetmelt16 = merge(rnetmelt16,group16,by.x="Var2",by.y="Bird")
rnetmelt17 = merge(rnetmelt17,group17,by.x="Var1",by.y="Bird")
rnetmelt17 = merge(rnetmelt17,group17,by.x="Var2",by.y="Bird")
rnetmelt18 = merge(rnetmelt18,group18,by.x="Var1",by.y="Bird")
rnetmelt18 = merge(rnetmelt18,group18,by.x="Var2",by.y="Bird")
rnetmelt19 = merge(rnetmelt19,group19,by.x="Var1",by.y="Bird")
rnetmelt19 = merge(rnetmelt19,group19,by.x="Var2",by.y="Bird")

##Cast it back into a matrix based on Social Group and sum values for like associations
#By using sum here you keep all of the interactions in the network. The actual network properties do not
#change, just the way you're showing them. See Arenas et al. 2007. 
rnet16 = as.matrix(acast(rnetmelt16,Social.Group.x~Social.Group.y,value.var = "value",fun.aggregate = sum))
rnet17 = as.matrix(acast(rnetmelt17,Social.Group.x~Social.Group.y,value.var = "value",fun.aggregate = sum))
rnet18 = as.matrix(acast(rnetmelt18,Social.Group.x~Social.Group.y,value.var = "value",fun.aggregate = sum))
rnet19 = as.matrix(acast(rnetmelt19,Social.Group.x~Social.Group.y,value.var = "value",fun.aggregate = sum))

#Set diagonals = 0
diag(rnet16) = 0
diag(rnet17) = 0
diag(rnet18) = 0
diag(rnet19) = 0

#Setup for igraph
rnet16i = graph.adjacency(rnet16, mode="undirected", diag=FALSE, weighted=TRUE)
rfg16 = fastgreedy.community(rnet16i)
rl16 <- layout.fruchterman.reingold(rnet16i)

rnet17i = graph.adjacency(rnet17, mode="undirected", diag=FALSE, weighted=TRUE)
rfg17 = fastgreedy.community(rnet17i)
rl17 <- layout.fruchterman.reingold(rnet17i)

rnet18i = graph.adjacency(rnet18, mode="undirected", diag=FALSE, weighted=TRUE)
rfg18 = fastgreedy.community(rnet18i)
rl18 <- layout.fruchterman.reingold(rnet18i)

rnet19i = graph.adjacency(rnet19, mode="undirected", diag=FALSE, weighted=TRUE)
rfg19 = fastgreedy.community(rnet19i)
rl19 <- layout.fruchterman.reingold(rnet19i)

#Get modularity values
modularity16 = modularity(rfg16)
modularity17 = modularity(rfg17)
modularity18 = modularity(rfg18)
modularity19 = modularity(rfg19)


#Plot social communities of groups where each group is a node
plot(rnet16i, vertex.size = 7, layout = layout_with_graphopt,
     vertex.label.cex = 0.5,
     edge.color = "black", vertex.label.color = "black",
     mark.groups=communities(rfg16))

plot(rnet17i, vertex.size = 7, layout = layout_with_graphopt,
     vertex.label.cex = 0.5,
     edge.color = "black", vertex.label.color = "black",
     mark.groups=communities(rfg17))

plot(rnet18i, vertex.size = 7, layout = layout_with_graphopt,
     vertex.label.cex = 0.5,
     edge.color = "black", vertex.label.color = "black",
     mark.groups=communities(rfg18))

plot(rnet19i, vertex.size = 7, layout = layout_with_graphopt,
     vertex.label.cex = 0.5,
     edge.color = "black", vertex.label.color = "black",
     mark.groups=communities(rfg19))



####Social classes per community####
library(tibble)

#Match up community identity to bird list for each group
rfg16m = data.frame(print(membership(rfg16)))
rfg16m = rownames_to_column(rfg16m,var="Group")
colnames(rfg16m)[2]="Community"
birdlist16noKsubc = merge(birdlist16noKsubc,rfg16m,by.x="Social.Group",by.y="Group")

rfg17m = data.frame(print(membership(rfg17)))
rfg17m = rownames_to_column(rfg17m,var="Group")
colnames(rfg17m)[2]="Community"
birdlist17noKsubc = merge(birdlist17noKsubc,rfg17m,by.x="Social.Group",by.y="Group")

rfg18m = data.frame(print(membership(rfg18)))
rfg18m = rownames_to_column(rfg18m,var="Group")
colnames(rfg18m)[2]="Community"
birdlist18noKsubc = merge(birdlist18noKsubc,rfg18m,by.x="Social.Group",by.y="Group")

rfg19m = data.frame(print(membership(rfg19)))
rfg19m = rownames_to_column(rfg19m,var="Group")
colnames(rfg19m)[2]="Community"
birdlist19noKsubc = merge(birdlist19noKsubc,rfg19m,by.x="Social.Group",by.y="Group")

#Groups per community for table
rfg16mt = table(rfg16m$Community)
min(rfg16mt)
max(rfg16mt)
mean(rfg16mt)

rfg17mt = table(rfg17m$Community)
min(rfg17mt)
max(rfg17mt)
mean(rfg17mt)

rfg18mt = table(rfg18m$Community)
min(rfg18mt)
max(rfg18mt)
mean(rfg18mt)

rfg19mt = table(rfg19m$Community)
min(rfg19mt)
max(rfg19mt)
mean(rfg19mt)

#Number of communities in each year 
length(rfg16mt)
length(rfg17mt)
length(rfg18mt)
length(rfg19mt)

#Birds per community for table
bpc16 = data.frame(table(birdlist16noKsubc$Community))
mean(bpc16$Freq)
min(bpc16$Freq)
max(bpc16$Freq)

bpc17 = data.frame(table(birdlist17noKsubc$Community))
mean(bpc17$Freq)
min(bpc17$Freq)
max(bpc17$Freq)

bpc18 = data.frame(table(birdlist18noKsubc$Community))
mean(bpc18$Freq)
min(bpc18$Freq)
max(bpc18$Freq)

bpc19 = data.frame(table(birdlist19noKsubc$Community))
mean(bpc19$Freq)
min(bpc19$Freq)
max(bpc19$Freq)

#Bring in age and sex status
birdlist16noKsubc = merge(birdlist16noKsubc,agesexstatus16,by="Bird",all.x=T)
birdlist17noKsubc = merge(birdlist17noKsubc,agesexstatus17,by="Bird",all.x=T)
birdlist18noKsubc = merge(birdlist18noKsubc,agesexstatus18,by="Bird",all.x=T)
birdlist19noKsubc = merge(birdlist19noKsubc,agesexstatus19,by="Bird",all.x=T)

#Get number of social classes per community
library(dplyr)
library(tidyr)
library(plyr)

comclass16 = count(birdlist16noKsubc,c("Community","class"))
comclass16 = complete(comclass16,Community,class,fill=list(freq=0))
comclass16$Year = "2016"
comclass17 = count(birdlist17noKsubc,c("Community","class"))
comclass17 = complete(comclass17,Community,class,fill=list(freq=0))
comclass17$Year = "2017"
comclass18 = count(birdlist18noKsubc,c("Community","class"))
comclass18 = complete(comclass18,Community,class,fill=list(freq=0))
comclass18$Year = "2018"
comclass19 = count(birdlist19noKsubc,c("Community","class"))
comclass19 = complete(comclass19,Community,class,fill=list(freq=0))
comclass19$Year = "2019"


#Combine years
comclassall = rbind(comclass16,comclass17,comclass18,comclass19)

###Resident Males
#Get Resident males
comclassallrm = comclassall[which(comclassall$class=="Resident Male"),]

#Get mean values
comclassallrmsum = data.frame(Year=c("2016","2017","2018","2019"),
                              Mean=c(mean(comclass16[which(comclass16$class=="Resident Male"),]$freq),
                                     mean(comclass17[which(comclass17$class=="Resident Male"),]$freq),
                                     mean(comclass18[which(comclass18$class=="Resident Male"),]$freq),
                                     mean(comclass19[which(comclass19$class=="Resident Male"),]$freq)))

#Plot number of resident males per community
A = ggplot(data=comclassallrm,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.1,h=0),size=2)  + geom_boxplot(fill=NA) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  geom_point(data=comclassallrmsum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  ggtitle("Resident Males") + ylab("Individuals per Community") +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none",
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5)) 


###Resident Females
#Get Resident females
comclassallrf = comclassall[which(comclassall$class=="Resident Female"),]

#Get mean values
comclassallrfsum = data.frame(Year=c("2016","2017","2018","2019"),
                              Mean=c(mean(comclass16[which(comclass16$class=="Resident Female"),]$freq),
                                     mean(comclass17[which(comclass17$class=="Resident Female"),]$freq),
                                     mean(comclass18[which(comclass18$class=="Resident Female"),]$freq),
                                     mean(comclass19[which(comclass19$class=="Resident Female"),]$freq)))

#Plot number of resident females per community
B = ggplot(data=comclassallrf,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.1,h=0),size=2)  + geom_boxplot(fill=NA) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  geom_point(data=comclassallrfsum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  ggtitle("Resident Females") + ylab("Individuals per Community") +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none",
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5))


###Auxilliary Males
comclassallax = comclassall[which(comclassall$class=="Auxilliary Male"),]

#Get mean values
comclassallaxsum = data.frame(Year=c("2016","2017","2018","2019"),
                              Mean=c(mean(comclass16[which(comclass16$class=="Auxilliary Male"),]$freq),
                                     mean(comclass17[which(comclass17$class=="Auxilliary Male"),]$freq),
                                     mean(comclass18[which(comclass18$class=="Auxilliary Male"),]$freq),
                                     mean(comclass19[which(comclass19$class=="Auxilliary Male"),]$freq)))

#Plot number of Auxilliary males per community
C = ggplot(data=comclassallax,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.1,h=0),size=2)  + geom_boxplot(fill=NA) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  geom_point(data=comclassallaxsum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  ggtitle("Auxilliary Males") + ylab("Individuals per Community") +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none",
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5))


###1 Year Old Males
#Get One Year Old Males
comclassall1m = comclassall[which(comclassall$class=="1 Year-old Male"),]

#Get mean values
comclassall1msum = data.frame(Year=c("2016","2017","2018","2019"),
                              Mean=c(mean(comclass16[which(comclass16$class=="1 Year-old Male"),]$freq),
                                     mean(comclass17[which(comclass17$class=="1 Year-old Male"),]$freq),
                                     mean(comclass18[which(comclass18$class=="1 Year-old Male"),]$freq),
                                     mean(comclass19[which(comclass19$class=="1 Year-old Male"),]$freq)))

#Plot number of one year old males per com
D = ggplot(data=comclassall1m,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.1,h=0),size=2)  + geom_boxplot(fill=NA) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  geom_point(data=comclassall1msum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  ggtitle("One Year-old Males") + ylab("Individuals per Community") +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none",
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5))


###Dispersing Females
#Get Dispersing Females
comclassalldf = comclassall[which(comclassall$class=="Dispersing Female"),]

#Get mean values
comclassalldfsum = data.frame(Year=c("2016","2017","2018","2019"),
                              Mean=c(mean(comclass16[which(comclass16$class=="Dispersing Female"),]$freq),
                                     mean(comclass17[which(comclass17$class=="Dispersing Female"),]$freq),
                                     mean(comclass18[which(comclass18$class=="Dispersing Female"),]$freq),
                                     mean(comclass19[which(comclass19$class=="Dispersing Female"),]$freq)))

#Plot number of Dispersing Females per community
E = ggplot(data=comclassalldf,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.1,h=0),size=2)  + geom_boxplot(fill=NA) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  geom_point(data=comclassalldfsum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  ggtitle("Dispersing Females") + ylab("Individuals per Community") +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none",
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5))

#Plot all graphs together
ggarrange(A,B,C,D)


#Get average number of birds per community for all years
bpc.all = rbind(bpc16,bpc17,bpc18,bpc19)
mean(bpc.all$Freq)

###Unknowns
comclassallunk = comclassall[which(comclassall$class=="Unknown"),]

#Get mean values
comclassallunksum = data.frame(Year=c("2016","2017","2018","2019"),
                               Mean=c(mean(comclass16[which(comclass16$class=="Unknown"),]$freq),
                                      mean(comclass17[which(comclass17$class=="Unknown"),]$freq),
                                      mean(comclass18[which(comclass18$class=="Unknown"),]$freq),
                                      mean(comclass19[which(comclass19$class=="Unknown"),]$freq)))

#Plot number of Unknowns per community
G = ggplot(data=comclassallunk,aes(x=Year,y=freq,fill=Year,colour=Year)) + 
  #geom_flat_violin(position = position_nudge(x=0.25,y=0),adjust=2,trim=F)+
  geom_point(position=position_jitter(width=0.1,h=0),size=2)  + geom_boxplot(fill=NA) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  geom_point(data=comclassallunksum,aes(x=Year,y=Mean),colour="black",size=2.5) +
  ggtitle("Unknowns") + ylab("Individuals per Community") +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none",
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5))


##Plot all six social classes
#Don't need a legend here because x axis specifies categories.
ggarrange(A,B,C,D,E,G) 



####Social Classes per community - all years in one plot####
comclassallsum = data.frame(class=c("Resident Male","Resident Female","Auxilliary Male",
                                      "1 Year-old Male","Dispersing Female","Unknown"),
                              Mean=c(mean(comclassall[which(comclassall$class=="Resident Male"),]$freq),
                                     mean(comclassall[which(comclassall$class=="Resident Female"),]$freq),
                                     mean(comclassall[which(comclassall$class=="Auxilliary Male"),]$freq),
                                     mean(comclassall[which(comclassall$class=="1 Year-old Male"),]$freq),
                                     mean(comclassall[which(comclassall$class=="Dispersing Female"),]$freq),
                                     mean(comclassall[which(comclassall$class=="Unknown"),]$freq)))

#Order factors
comclassall$class = factor(comclassall$class,levels=c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male",
                                                          "Dispersing Female","Unknown"))

#Get multi-level labels for each factor
comclassall$class = as.factor(comclassall$class)
levels(comclassall$class) <- gsub(" ", "\\n", levels(comclassall$class))
levels(comclassallsum$class) <- gsub(" ", "\\n", levels(comclassallsum$class))


ggplot(data=comclassall,aes(x=class,y=freq)) + 
  geom_boxplot(data=comclassall,aes(x=class,y=freq),colour="black",width=0.3) +
  geom_flat_violin(aes(group=class),position = position_nudge(x=0.27,y=0),trim=F) +
  #geom_point(position=position_jitter(width=0.15,h=0.07),size=1.5,color="dark gray") +
  ylab("Individuals per Social Community") + ylim(-0.7,8) + xlab("Social Class")  +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" )) # explicitly set the horizontal lines

ggplot(data=comclassall,aes(x=class,y=freq,colour=Year)) + 
  geom_flat_violin(aes(group=class),position = position_nudge(x=0.27,y=0),trim=F) +
  geom_point(position=position_jitter(width=0.15,h=0.07),size=1.5) +
  geom_boxplot(data=comclassall,aes(x=class,y=freq),colour="black",width=0.3,fill="transparent") +
  ylab("Individuals per Social Community")  + ylim(-0.7,8) + 
  xlab("Social Class") + scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  theme_cowplot() + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" )) # explicitly set the horizontal lines





#__________________________________________________________________________#




####Look at average association index within versus across groups####

#Putting this here because I want to use birdlistnoksubc - birds that that were listed in a group and not on their own
#With this analysis want to know what the average association index of within group connecitons is compared to across group connections
#rnetmelt is already the melted network with group IDs merged and still has association index

#Function to measure associations in vs across groups

within_vs_across = function(rnetmelt) {
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
  
  return(list(rnetmelt.same,rnetmelt.different))
}

within_across16 = within_vs_across(rnetmelt = rnetmelt16)
within_across17 = within_vs_across(rnetmelt = rnetmelt17)
within_across18 = within_vs_across(rnetmelt = rnetmelt18)
within_across19 = within_vs_across(rnetmelt = rnetmelt19)

within16 = within_across16[[1]]
within17 = within_across17[[1]]
within18 = within_across18[[1]]
within19 = within_across19[[1]]

withinall = rbind(within16,within17,within18,within19)
mean(withinall$SRI)
min(withinall$SRI)
max(withinall$SRI)

across16 = within_across16[[2]]
across17 = within_across17[[2]]
across18 = within_across18[[2]]
across19 = within_across19[[2]]

acrossall = rbind(across16,across17,across18,across19)
mean(acrossall$SRI)
min(acrossall$SRI)
max(acrossall$SRI)




#__________________________________________________________________________#


####Write files needed for other scripts####


#GBInoK files - Need to delete first column in .csv file
# write.csv(gbi16noK,file=here::here("Output files","gbi16noK.csv"))
# write.csv(gbi17noK,file=here::here("Output files","gbi17noK.csv"))
# write.csv(gbi18noK,file=here::here("Output files","gbi18noK.csv"))
# write.csv(gbi19noK,file=here::here("Output files","gbi19noK.csv"))

#BirdlistnoK files - delete first column in .csv file
# write.csv(birdlist16noKsub,file=here::here("Output files","birdlist16noKsub.csv"))
# write.csv(birdlist17noKsub,file=here::here("Output files","birdlist17noKsub.csv"))
# write.csv(birdlist18noKsub,file=here::here("Output files","birdlist18noKsub.csv"))
# write.csv(birdlist19noKsub,file=here::here("Output files","birdlist19noKsub.csv"))

#Average silhoutte distance
# avg.sw.all = data.frame(year= c(2016,2017,2018,2019),
# avg.s = c(avg.s16,avg.s17,avg.s18,avg.s19))
# write.csv(avg.sw.all,file=here::here("Output files","Average Silhouette Widths.csv"))

#GBI all sightings files 
# write.csv(gbi16,file=here::here("Output files","gbi16.csv"))
# write.csv(gbi17,file=here::here("Output files","gbi17.csv"))
# write.csv(gbi18,file=here::here("Output files","gbi18.csv"))
# write.csv(gbi19,file=here::here("Output files","gbi19.csv"))

#Reduced matrix files
# write.csv(rnet16,file=here::here("Output files","rnet16.csv"))
# write.csv(rnet17,file=here::here("Output files","rnet17.csv"))
# write.csv(rnet18,file=here::here("Output files","rnet18.csv"))
# write.csv(rnet19,file=here::here("Output files","rnet19.csv"))

#Group membership
# write.csv(group16,file=here::here("Output files","group16.csv"))
# write.csv(group17,file=here::here("Output files","group17.csv"))
# write.csv(group18,file=here::here("Output files","group18.csv"))
# write.csv(group19,file=here::here("Output files","group19.csv"))

#Community membership
# write.csv(rfg16m,file=here::here("Output files","communities16.csv"))
# write.csv(rfg17m,file=here::here("Output files","communities17.csv"))
# write.csv(rfg18m,file=here::here("Output files","communities18.csv"))
# write.csv(rfg19m,file=here::here("Output files","communities19.csv"))

#Birdlists with birds in communities
# write.csv(birdlist16noKsubc,file=here::here("Output files","birdlist16noKsubc.csv"))
# write.csv(birdlist17noKsubc,file=here::here("Output files","birdlist17noKsubc.csv"))
# write.csv(birdlist18noKsubc,file=here::here("Output files","birdlist18noKsubc.csv"))
# write.csv(birdlist19noKsubc,file=here::here("Output files","birdlist19noKsubc.csv"))


#Write modularity scores to a dataframe
# mod.all = data.frame(year= c(2016,2017,2018,2019),
# modularity = c(modularity16,modularity17,modularity18,modularity19))
# write.csv(mod.all,file=here::here("Output files","modularity_scores.csv"))


####Counts for Supplemental Table 1####

#Count sampling points (sightings) for table 1 group network
length(table(bird16noK$Sighting))
length(table(bird18noK$Sighting))
length(table(bird19noK$Sighting))

#For 2017, remove removal birds - want number of sampling points and observations that the control birds
#were involved in for table 1 - number of observations that went towards calculation of groups
bird17noKnoR = bird17noK[!bird17noK$Bird %in% removalbirds$Bird,]
length(table(bird17noKnoR$Sighting))

#Count observations for table 1 group network
length(table(bird16noK$Observation))
length(table(bird17noKnoR$Observation))
length(table(bird18noK$Observation))
length(table(bird19noK$Observation))

#Count sampling points (sightings) for table 1 community network
length(table(bird16$Sighting))
length(table(bird18$Sighting))
length(table(bird19$Sighting))

#For 2017 remove removal birds
bird17noR = bird17R[!bird17R$Bird %in% removalbirds$Bird,]
length(table(bird17noR$Sighting))

#Count observations for table 1 community network
length(table(bird16$Observation))
length(table(bird17noR$Observation))
length(table(bird18$Observation))
length(table(bird19$Observation))

#Mean number of sampling points per observation
sppo16 = bird16 %>% distinct(Sighting,Observation)
sppo16 = data.frame(table(sppo16$Observation))
mean(sppo16$Freq)
sppo17 = bird17R %>% distinct(Sighting,Observation)
sppo17 = data.frame(table(sppo17$Observation))
mean(sppo17$Freq)
sppo18 = bird18 %>% distinct(Sighting,Observation)
sppo18 = data.frame(table(sppo18$Observation))
mean(sppo18$Freq)
sppo19 = bird19 %>% distinct(Sighting,Observation)
sppo19 = data.frame(table(sppo19$Observation))
mean(sppo19$Freq)

