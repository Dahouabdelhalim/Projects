##Assortment of groups across seasons 2

#Use this file to look at the relationships between non-breeding social structure and breeding groups. 
#Two tests

#1. What is the assortment of breeding group members during the non-breeding season? This compares breeding group
#membership to network structure - association indices. Are previous and upcoming breeding groups more assorted
#duing the non-breeding season than expected by chance? And if so which is more assorted? Although end up having
#different birds to compare so not able to make many conclusions about if NBS more reflects previous or upcoming 
#BS. 

#2. Given that previous breeding groups are more assorted during the non-breeding season than expected by chance, 
#what percentage of dyads from previous breeding groups are in the same non-breeding group?
#What percentage of non-breeding dyads are in the same upcoming breeding group? 



#Load birdlists from network

library(here)

#Load birdlists with birds in communities and breeding groups - birds in a non-breeding group on their own have been removed
birdlist16noKsubcb = read.csv(here::here("Output files","birdlist16noKsubcb.csv"))
birdlist17noKsubcb = read.csv(here::here("Output files","birdlist17noKsubcb.csv"))
birdlist18noKsubcb = read.csv(here::here("Output files","birdlist18noKsubcb.csv"))
birdlist19noKsubcb = read.csv(here::here("Output files","birdlist19noKsubcb.csv"))




####Assortment code####

#Function to compare assortment of previous breeding groups and upcoming breeding groups within the non-breeding network
#and to see if these assortment values are different from expected by chance. 
#Since we already know the network is very structured, doesn't make sense to do a datastream permutation here
#like had originally. I want ask, given the network structure, are the connections within the network more driven by 
#previous breeding season group structure or upcoming breeding season group structure? 
#and are these values different from chance? 
#So it seems like a node-permutation will be best here. First version of this script used a datastream permutation
#But within the node-permutation, only swap individuals of the same class to keep the network structured in a similar way

#Tried to write this function so that only birds with previous breeding groups, seen enough times in non-breeding network, and 
#had upcoming breeding groups would be run, so that comparisons to previous and upcoming seasons would use the same individuals.
#However sample size was extremely small when tried that - 11 to 33 individuals per year. So went now comparisons to previous 
#and upcoming breeding groups do not use the same individuals, some overlap but some differences to get a larger sample size. 
#Now birds just need to have been seen enough times in non-breeding network and have either a prev or upcoming breeding group. 

#Function for assortment of breeding and non-breeding groups
#Structure function like groups compared to random. 
#1. Generate observed network from GBI 
#2. Subset for birds that have previous and upcoming breeding groups
#3. Get assortment of previous breeding groups in network and upcoming breeding groups in network
#4. Randomize 1000 times within class which prev breeding or upcoming breeding group each bird belongs to
#5. Get assortment for the randomized networks
#6. Compare observed to random for both prev and upcoming
#7. Compare observed assortment values of prev and upcoming breeding groups

library(asnipe)
library(assortnet)
library(dplyr)
library(reshape2)


#Load data 
gbi16 = read.csv(here::here("Output files","gbi16.csv"))
gbi17 = read.csv(here::here("Output files","gbi17.csv"))
gbi18 = read.csv(here::here("Output files","gbi18.csv"))
gbi19 = read.csv(here::here("Output files","gbi19.csv"))


#Function
assortment.brgroups = function(gbi,perm.number,birdlist) {
  #Get network and order by bird
  network <- get_network(gbi,data_format="GBI", association_index="SRI")
  network = network[order(rownames(network)),order(colnames(network))]
  names = rownames(network)
  
  #Get list of birds with previous breeding groups and make sure all breeding groups occur at least twice
  #They have to occur at least twice to measure assortment between them
  birdlist.pr = birdlist[which(birdlist$PrevBrGrp!=""),]
  birdlist.pr = birdlist.pr[!is.na(birdlist.pr$Social.Group),]
  dups.pr = birdlist.pr[duplicated(birdlist.pr$PrevBrGrp),]$PrevBrGrp
  birdlist.pr = birdlist.pr[birdlist.pr$PrevBrGrp %in% dups.pr,]
  
  #Get list of birds with upcoming breeding groups and make sure all breeding groups occur at least twice
  #They have to occur at least twice to measure assortment between them
  birdlist.up = birdlist[which(birdlist$UpcomingBrGrp!=""),]
  birdlist.up = birdlist.up[!is.na(birdlist.up$Social.Group),]
  dups.up = birdlist.up[duplicated(birdlist.up$UpcomingBrGrp),]$UpcomingBrGrp
  birdlist.up = birdlist.up[birdlist.up$UpcomingBrGrp %in% dups.up,]
  
  
  #Subset network based on birds with previous breeding group and seen enough times (in birdlist)
  network.pr = network[which(names %in% birdlist.pr$Bird),which(names %in% birdlist.pr$Bird)]
  network.up = network[which(names %in% birdlist.up$Bird),which(names %in% birdlist.up$Bird)]
  
  #Get assortment of previous and upcoming breeding groups in network
  assort.obs.pr = assortment.discrete(network.pr,birdlist.pr$PrevBrGrp,weighted = T,SE=FALSE)$r
  assort.obs.up = assortment.discrete(network.up,birdlist.up$UpcomingBrGrp,weighted = T,SE=FALSE)$r
  
  #Randomize breeding groups but only swap within classes
  assort.rand.pr = rep(0,perm.number)
  assort.rand.up = rep(0,perm.number)
  
  for (i in c(1:perm.number)) {
    #Order by class
    birdlist.pr.r = birdlist.pr[order(birdlist.pr$class),]
    #Randomly sample within group without replacement - this randomly orders within class, but keeps class order the same
    birdlist.pr.r2 = birdlist.pr %>% group_by(class) %>% sample_frac(replace=F) %>% select(class,PrevBrGrp) 
    #Then add the random prev breeding gorup to the birdlist dataframe
    birdlist.pr.r$PrevBrGrp.r = birdlist.pr.r2$PrevBrGrp
    #Order the birdlist again by bird because assortment just matches the lists
    birdlist.pr.r = birdlist.pr.r[order(birdlist.pr.r$Bird),]
    #Then run assortment for the random previous breeding groups
    assort.rand.pr[i] = assortment.discrete(network.pr,birdlist.pr.r$PrevBrGrp.r,weighted = T,SE=FALSE)$r
  }
  
  
  for (i in c(1:perm.number)) {
    #Order by class
    birdlist.up.r = birdlist.up[order(birdlist.up$class),]
    #Randomly sample within group without replacement - this randomly orders within class, but keeps class order the same
    birdlist.up.r2 = birdlist.up %>% group_by(class) %>% sample_frac(replace=F) %>% select(class,UpcomingBrGrp) 
    #Then add the random prev breeding gorup to the birdlist dataframe
    birdlist.up.r$UpcomingBrGrp.r = birdlist.up.r2$UpcomingBrGrp
    #Order the birdlist again by bird because assortment just matches the lists
    birdlist.up.r = birdlist.up.r[order(birdlist.up.r$Bird),]
    #Then run assortment for the random previous breeding groups
    assort.rand.up[i] = assortment.discrete(network.up,birdlist.up.r$UpcomingBrGrp.r,weighted = T,SE=FALSE)$r
  }
  
  par(mfrow=c(1,2))
  #Compare observed to random
  hist(assort.rand.pr,breaks=25,xlim=c(-0.5,1),main="Previous vs Random",xlab="Assortment")
  abline(v=assort.obs.pr,col="red",lwd=3)
  hist(assort.rand.up,breaks=25,xlim=c(-0.5,1),main="Upcoming vs Random",xlab="Assortment")
  abline(v=assort.obs.up,col="red",lwd=3)
  par(mfrow=c(1,1))
  
  #calculate p-value for both
  P.pr = 1 - sum(assort.obs.pr > assort.rand.pr)/perm.number
  cat("\\n")
  print(P.pr)
  P.up = 1 - sum(assort.obs.up > assort.rand.up)/perm.number
  print(P.up)
  return(c(assort.obs.pr,assort.obs.up))
}

#Run the function to look at assortment of breeding and non-breeding groups
# assortment.brgroups(gbi=gbi16,perm.number = 1000,birdlist=birdlist16noKsubcb) #p<0.001 p<0.001
# assortment.brgroups(gbi=gbi17,perm.number = 1000,birdlist=birdlist17noKsubcb) #p<0.001 p<0.001
# assortment.brgroups(gbi=gbi18,perm.number = 1000,birdlist=birdlist18noKsubcb) #p<0.001 p<0.001
# assortment.brgroups(gbi=gbi19,perm.number = 1000,birdlist=birdlist19noKsubcb) #p<0.001 p<0.001

#In all years non-breeding groups are more assorted with the past breeding season than the upcoming
#and more assorted than expected by chance. But since we're using different combinations of individuals
#to calculate assortment between non-breeding and previous and non-breeding and upcoming breeding groups, 
#cannot really say for sure which one is more assorted. 

#Get average assortment of the four years
assort16 = assortment.brgroups(gbi=gbi16,perm.number = 1,birdlist=birdlist16noKsubcb)
assort16.pr = assort16[1]
assort16.up = assort16[2]
assort17 = assortment.brgroups(gbi=gbi17,perm.number = 1,birdlist=birdlist17noKsubcb)
assort17.pr = assort17[1]
assort17.up = assort17[2]
assort18 = assortment.brgroups(gbi=gbi18,perm.number = 1,birdlist=birdlist18noKsubcb)
assort18.pr = assort18[1]
assort18.up = assort18[2]
assort19 = assortment.brgroups(gbi=gbi19,perm.number = 1,birdlist=birdlist19noKsubcb)
assort19.pr = assort19[1]
assort19.up = assort19[2]

mean(assort16.pr,assort17.pr,assort18.pr,assort19.pr)
mean(assort16.up,assort17.up,assort18.up,assort19.up)

#Report average assortment in results - previous: 0.84, upcoming: 0.73

#Now that we know previous and upcoming breeding groups are more assorted during the non-breeding season than expected by chance, see
#what relationships are respsonsible for that assortment.








####Percentage of dyads same across seasons code####

#Tried to run previous breeding groups such that only one bird of a dyad had to be recorded in a non-breeding group and 
#the other bird could not be recorded in a non-breeding group, concluding that if they were truly in the same group we would've seen them both enough
#times to include in the network, but that ended up saying some known male-female pairs with decent connections during the 
#non-breeding season and were together in both breeding seasons on either side were not in the same non-breeding group. 
#So requiring that both birds were identified to a non-breeding group, sample size will be lower but can be more confident and 
#should be more accurate. 


#Load in each individual's upcoming breeding group - the first group it was associated with in the upcoming breeding season
upstatus16 = read.csv(here::here("Input files","Upcoming_status16.csv")) %>% select(Bird,Status)
colnames(upstatus16)[2] = "Upstatus"
upstatus17 = read.csv(here::here("Input files","Upcoming_status17.csv")) %>% select(Bird,Status)
colnames(upstatus17)[2] = "Upstatus"
upstatus18 = read.csv(here::here("Input files","Upcoming_status18.csv")) %>% select(Bird,Status)
colnames(upstatus18)[2] = "Upstatus"
upstatus19 = read.csv(here::here("Input files","Upcoming_status19.csv")) %>% select(Bird,Status)
colnames(upstatus19)[2] = "Upstatus"

#Load in agesexstatus files
agesexstatus16 = read.csv(here::here("Input files","agesexstatus16.csv"))
agesexstatus17 = read.csv(here::here("Input files","agesexstatus17.csv"))
agesexstatus18 = read.csv(here::here("Input files","agesexstatus18.csv"))
agesexstatus19 = read.csv(here::here("Input files","agesexstatus19.csv"))

#Merge upcoming status and birdlist
birdlist16noKsubcb = merge(birdlist16noKsubcb,upstatus16,by="Bird",all.x=T)
birdlist17noKsubcb = merge(birdlist17noKsubcb,upstatus17,by="Bird",all.x=T)
birdlist18noKsubcb = merge(birdlist18noKsubcb,upstatus18,by="Bird",all.x=T)
birdlist19noKsubcb = merge(birdlist19noKsubcb,upstatus19,by="Bird",all.x=T)


#Function
compare_dyads = function(birdlist,agesexstatus) {
  
  #Get agesexstatus columns needed
  agesexstatus = agesexstatus %>% select(Bird,class)
  
  #Get list of birds with previous breeding groups and a non-breeding group
  birdlist.pr = birdlist[which(birdlist$PrevBrGrp!=""),]
  birdlist.pr = birdlist.pr %>% select(Bird,Social.Group,PrevBrGrp)
  birdlist.pr = birdlist.pr[!is.na(birdlist.pr$Social.Group),]
  
  
  #For upcoming statuses, need to know if a bird is a first time breeder or a previous breeder or a helper in the upcoming season (first breeding group) 
  #So first take out Unknown class birds from status - we don't know what they were the previous season - so they could've bred before
  #Also take out upcoming status = other - don't know really what class they were at the beginning of the season
  birdlist.up = merge(birdlist,agesexstatus,by="Bird",all=T)
  birdlist.up = birdlist.up[which(birdlist.up$class.y!="Unknown"),]
  birdlist.up = birdlist.up[which(birdlist.up$Upstatus!="Other"),]

  #Get list of birds with upcoming breeding groups and non-breeding social groups
  birdlist.up = birdlist.up[which(birdlist.up$UpcomingBrGrp!=""),]
  birdlist.up = birdlist.up[!is.na(birdlist.up$Social.Group),]
  

  #Now get an updated Upstatus based on first time breeding or not
  birdlist.up$Upstatus2 = NA
  for (i in 1:nrow(birdlist.up)) {
    if(birdlist.up$class.y[i] == "Resident Male" | birdlist.up$class.y[i]=="Resident Female") 
      {birdlist.up$Upstatus2[i]=paste("Old",birdlist.up$Upstatus[i])} else 
      {if(birdlist.up$Upstatus[i]=="Helper") {birdlist.up$Upstatus2[i]="Helper"} else 
      {birdlist.up$Upstatus2[i]=paste("First",birdlist.up$Upstatus[i])}}
  }
  birdlist.up = birdlist.up %>% select(Bird,Social.Group,UpcomingBrGrp,class.y,Upstatus2)

  #Get all combinations of birds
  dyads.pr = expand.grid(birdlist.pr$Bird,birdlist.pr$Bird)
  dyads.pr = data.frame(t(apply(dyads.pr, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(dyads.pr) = c("Bird1","Bird2")
  dyads.pr = dyads.pr[dyads.pr$Bird1!=dyads.pr$Bird2,] #Remove pairs where bird is the same in both columns
  dyads.pr = dyads.pr %>% distinct(Bird1,Bird2) #Remove duplicate pairs

  #Get all combinations of birds
  dyads.up = expand.grid(birdlist.up$Bird,birdlist.up$Bird)
  dyads.up = data.frame(t(apply(dyads.up, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(dyads.up) = c("Bird1","Bird2")
  dyads.up = dyads.up[dyads.up$Bird1!=dyads.up$Bird2,] #Remove pairs where bird is the same in both columns
  dyads.up = dyads.up %>% distinct(Bird1,Bird2) #Remove duplicate pairs
  
  #Merge variables for both birds
  dyads.pr = merge(dyads.pr,birdlist.pr,by.x="Bird1",by.y="Bird")
  dyads.pr = merge(dyads.pr,birdlist.pr,by.x="Bird2",by.y="Bird")
  
  #Merge variables for both birds
  dyads.up = merge(dyads.up,birdlist.up,by.x="Bird1",by.y="Bird")
  dyads.up = merge(dyads.up,birdlist.up,by.x="Bird2",by.y="Bird")
  
  
  #For previous breeding groups, get birds that are in the same previous breeding group 
  dyads.pr.same = dyads.pr[dyads.pr$PrevBrGrp.x==dyads.pr$PrevBrGrp.y,] #Same previous breeding group
  dyads.pr.same = merge(dyads.pr.same,agesexstatus,by.x="Bird2",by.y="Bird")
  dyads.pr.same = merge(dyads.pr.same,agesexstatus,by.x="Bird1",by.y="Bird")
  
  #For upcoming breeding groups, get birds that are in the same upcoming breeding group
  dyads.up.same = dyads.up[dyads.up$UpcomingBrGrp.x==dyads.up$UpcomingBrGrp.y,] 
  
  #See if in same group
  dyads.pr.same$NBgroup.same = NA
  for (i in 1:nrow(dyads.pr.same)) {
    if(dyads.pr.same$Social.Group.x[i]==dyads.pr.same$Social.Group.y[i]) {dyads.pr.same$NBgroup.same[i]="Yes"} else {dyads.pr.same$NBgroup.same[i]="No"}
  }
  
  dyads.up.same$NBgroup.same = NA
  for (i in 1:nrow(dyads.up.same)) {
    if(dyads.up.same$Social.Group.x[i]==dyads.up.same$Social.Group.y[i]) {dyads.up.same$NBgroup.same[i]="Yes"} else {dyads.up.same$NBgroup.same[i]="No"}
  }
  
  #Order classes within rows so they match the class list order within rows later - need this for calculating percentages
  #Otherwise would have duplicates for each pair of classes and would need to adjust for that
  #This now mismatches to the bird IDs but don't need them anymore, so just grab whether previous breeding group was the same and add in
  dyads.pr.same.classes = dyads.pr.same %>% select(class.x,class.y)
  dyads.pr.same.classes = data.frame(t(apply(dyads.pr.same.classes, 1, sort)),dyads.pr.same$NBgroup.same)
  colnames(dyads.pr.same.classes) = c("class.x","class.y","NBgroup.same") 
  
  dyads.up.same.classes = dyads.up.same %>% select(Upstatus2.x, Upstatus2.y)
  dyads.up.same.classes = data.frame(t(apply(dyads.up.same.classes, 1, sort)),dyads.up.same$NBgroup.same)
  colnames(dyads.up.same.classes) = c("class.x","class.y","UpcomingBrGrp.same")
  
  ##Previous breeding groups
  #Get a list of previous class combinations
  classes.pr = c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male","Dispersing Female")
  classes.pr = expand.grid(classes.pr,classes.pr)
  classes.pr = data.frame(t(apply(classes.pr, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(classes.pr) = c("class.x","class.y")
  classes.pr = classes.pr %>% distinct(class.x,class.y)
  classes.pr$percent = NA
  
  #Calculate the percentages for each social class - how many non-breeding within-group dyads were in the 
  #same previous breeding group?

  for (i in 1:nrow(classes.pr)) {
    cl = dyads.pr.same.classes[which(dyads.pr.same.classes$class.x %in% classes.pr$class.x[i] & 
                              dyads.pr.same.classes$class.y %in% classes.pr$class.y[i]),]
    count = nrow(cl[which(cl$NBgroup.same=="Yes"),])
    if (count==0) {classes.pr$percent[i]=0} else {
      classes.pr$percent[i]=count/nrow(cl)
    }
  }
  
  #Put it into an array
  pr.percentage.array = as.matrix(acast(classes.pr,class.x~class.y,value.var="percent",fun.aggregate = sum,fill=NA_real_))
  
  #Now do the same thing, but for counts
  classes.pr$count = NA
  classes.pr$total = NA
  
  for (h in 1:nrow(classes.pr)) {
    cl = dyads.pr.same.classes[which(dyads.pr.same.classes$class.x %in% classes.pr$class.x[h] & 
                                       dyads.pr.same.classes$class.y %in% classes.pr$class.y[h]),]
    count = nrow(cl[which(cl$NBgroup.same=="Yes"),])
    classes.pr$count[h]=count
    classes.pr$total[h]=nrow(cl)
  }
  pr.count.array = as.matrix(acast(classes.pr,class.x~class.y,value.var ="count",fun.aggregate = sum,fill=NA_real_))
  pr.total.array = as.matrix(acast(classes.pr,class.x~class.y,value.var ="total",fun.aggregate = sum,fill=NA_real_))
  
  
  ##Upcoming breeding groups
  #Get a list of upcoming class combinations
  classes.up = c("Old Male","Old Female","Helper","First Male","First Female")
  classes.up = expand.grid(classes.up,classes.up)
  classes.up = data.frame(t(apply(classes.up, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(classes.up) = c("class.x","class.y")
  classes.up = classes.up %>% distinct(class.x,class.y)
  classes.up$percent = NA
  
  #Calculate the percentages for each social class - how many non-breeding within-group dyads were in the 
  #same previous breeding group?
  
  for (i in 1:nrow(classes.up)) {
    cl = dyads.up.same.classes[which(dyads.up.same.classes$class.x %in% classes.up$class.x[i] & 
                                       dyads.up.same.classes$class.y %in% classes.up$class.y[i]),]
    count = nrow(cl[which(cl$UpcomingBrGrp.same=="Yes"),])
    if (count==0) {classes.up$percent[i]=0} else {
      classes.up$percent[i]=count/nrow(cl)
    }
  }
  
  #Put it into an array
  up.percentage.array = as.matrix(acast(classes.up,class.x~class.y,value.var="percent",fun.aggregate = sum,fill=NA_real_))
  
  #Now do the same thing, but for counts
  classes.up$count = NA
  classes.up$total = NA
  
  for (h in 1:nrow(classes.up)) {
    cl = dyads.up.same.classes[which(dyads.up.same.classes$class.x %in% classes.up$class.x[h] & 
                                       dyads.up.same.classes$class.y %in% classes.up$class.y[h]),]
    count = nrow(cl[which(cl$UpcomingBrGrp.same=="Yes"),])
    if (count==0) {classes.up$count[h]=0} else {
      classes.up$count[h]=count}
    classes.up$total[h]=nrow(cl)
  }
  up.count.array = as.matrix(acast(classes.up,class.x~class.y,value.var ="count",fun.aggregate = sum,fill=NA_real_))
  up.total.array = as.matrix(acast(classes.up,class.x~class.y,value.var ="total",fun.aggregate = sum,fill=NA_real_))
  
  
  
  toreturn = list(pr.percentage.array,pr.count.array,pr.total.array,up.percentage.array,up.count.array,up.total.array,dyads.pr.same,dyads.up.same)
  print(pr.percentage.array)
  print(pr.count.array)
  print(pr.total.array)
  print(up.percentage.array)
  print(up.count.array)
  print(up.total.array)
  return(toreturn)
  
}

#Run for all years
dyads16 = compare_dyads(birdlist = birdlist16noKsubcb,agesexstatus = agesexstatus16)
dyads17 = compare_dyads(birdlist = birdlist17noKsubcb,agesexstatus = agesexstatus17)
dyads18 = compare_dyads(birdlist = birdlist18noKsubcb,agesexstatus = agesexstatus18)
dyads19 = compare_dyads(birdlist = birdlist19noKsubcb,agesexstatus = agesexstatus19)

##Previous 
#Get counts for each year and combine
dyads16.count.pr = dyads16[[2]]
dyads17.count.pr = dyads17[[2]]
dyads18.count.pr = dyads18[[2]]
dyads19.count.pr = dyads19[[2]]
dyads.count.pr = dyads16.count.pr + dyads17.count.pr + dyads18.count.pr + dyads19.count.pr

#Get totals for each year and combine
dyads16.total.pr = dyads16[[3]]
dyads17.total.pr = dyads17[[3]]
dyads18.total.pr = dyads18[[3]]
dyads19.total.pr = dyads19[[3]]
dyads.total.pr = dyads16.total.pr + dyads17.total.pr + dyads18.total.pr + dyads19.total.pr

#Get percentage of previous group members in the same non-breeding group, all years combined
dyads.percent.pr = dyads.count.pr/dyads.total.pr
print(dyads.percent.pr)
print(dyads.count.pr)
print(dyads.total.pr)


#Two resident females in group together in 2016 bs and 2017 bs
rf2.17 = dyads17[[7]] #2016 - B636 looks like male had two females, this is WWV's group - got a second female partway through season
rf2.18 = dyads18[[7]] #2017 - B643 - second female looked like a helper according to techs. She was with her parents. May have just been
#a really late disperser. 


##Upcoming
#Get counts for each year and combine
dyads16.count.up = dyads16[[5]]
dyads17.count.up = dyads17[[5]]
dyads18.count.up = dyads18[[5]]
dyads19.count.up = dyads19[[5]]
dyads.count.up = dyads16.count.up + dyads17.count.up + dyads18.count.up + dyads19.count.up

#Get totals for each year and combine
dyads16.total.up = dyads16[[6]]
dyads17.total.up = dyads17[[6]]
dyads18.total.up = dyads18[[6]]
dyads19.total.up = dyads19[[6]]
dyads.total.up = dyads16.total.up + dyads17.total.up + dyads18.total.up + dyads19.total.up

#Get percentage of previous group members in the same non-breeding group, all years combined
dyads.percent.up = dyads.count.up/dyads.total.up
print(dyads.percent.up)
print(dyads.count.up)
print(dyads.total.up)



####How many of the resident male-resident female pairs (previous breeders) are new pairs?####
dyads16.pr = dyads16[[7]]
dyads16.up = dyads16[[8]]
dyads17.pr = dyads17[[7]]
dyads17.up = dyads17[[8]]
dyads18.pr = dyads18[[7]]
dyads18.up = dyads18[[8]]
dyads19.pr = dyads19[[7]]
dyads19.up = dyads19[[8]]

#Get old birds - previous breeders
dyads16.up.old = dyads16.up %>% filter(Upstatus2.x=="Old Male" | Upstatus2.x=="Old Female")
dyads16.up.old = dyads16.up.old %>% filter(Upstatus2.y=="Old Male" | Upstatus2.y=="Old Female")
dyads17.up.old = dyads17.up %>% filter(Upstatus2.x=="Old Male" | Upstatus2.x=="Old Female")
dyads17.up.old = dyads17.up.old %>% filter(Upstatus2.y=="Old Male" | Upstatus2.y=="Old Female")
dyads18.up.old = dyads18.up %>% filter(Upstatus2.x=="Old Male" | Upstatus2.x=="Old Female")
dyads18.up.old = dyads18.up.old %>% filter(Upstatus2.y=="Old Male" | Upstatus2.y=="Old Female")
dyads19.up.old = dyads19.up %>% filter(Upstatus2.x=="Old Male" | Upstatus2.x=="Old Female")
dyads19.up.old = dyads19.up.old %>% filter(Upstatus2.y=="Old Male" | Upstatus2.y=="Old Female")

#Merge previous breeding group
birdlist16noKsubcb.prev = birdlist16noKsubcb %>% select(Bird,PrevBrGrp)
birdlist17noKsubcb.prev = birdlist17noKsubcb %>% select(Bird,PrevBrGrp)
birdlist18noKsubcb.prev = birdlist18noKsubcb %>% select(Bird,PrevBrGrp)
birdlist19noKsubcb.prev = birdlist19noKsubcb %>% select(Bird,PrevBrGrp)

dyads16.up.old = merge(dyads16.up.old,birdlist16noKsubcb.prev,by.x="Bird2",by.y="Bird")
dyads16.up.old = merge(dyads16.up.old,birdlist16noKsubcb.prev,by.x="Bird1",by.y="Bird")
dyads17.up.old = merge(dyads17.up.old,birdlist17noKsubcb.prev,by.x="Bird2",by.y="Bird")
dyads17.up.old = merge(dyads17.up.old,birdlist17noKsubcb.prev,by.x="Bird1",by.y="Bird")
dyads18.up.old = merge(dyads18.up.old,birdlist18noKsubcb.prev,by.x="Bird2",by.y="Bird")
dyads18.up.old = merge(dyads18.up.old,birdlist18noKsubcb.prev,by.x="Bird1",by.y="Bird")
dyads19.up.old = merge(dyads19.up.old,birdlist19noKsubcb.prev,by.x="Bird2",by.y="Bird")
dyads19.up.old = merge(dyads19.up.old,birdlist19noKsubcb.prev,by.x="Bird1",by.y="Bird")

dyads.up.old.all = rbind(dyads16.up.old,dyads17.up.old,dyads18.up.old,dyads19.up.old)
dyads.up.old.all %>% count(PrevBrGrp.x==PrevBrGrp.y)
#28 were in the same group, 7 were definitely not, and 1 pair was probably not - male doesn't have a breeding group 
#but female does. Go with 28 and 8. 













