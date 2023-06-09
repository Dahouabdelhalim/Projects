###Calculate 1yo non-breeding social connections script

#Using this script to calculate observed non-breeding social connections of young males. 
#Files created this script will be used to determine whether and how non-breeding social environment influences
#whether and when males molt into nuptial plumage

library(here)

####Load data####

##Load bird files - non-breeding social network data - come from non-breeding social structure paper
bird16 = read.csv(here::here("Input files","bird16.csv"))
bird17R = read.csv(here::here("Input files","bird17R.csv"))
bird18 = read.csv(here::here("Input files","bird18.csv"))
bird19 = read.csv(here::here("Input files","bird19.csv"))

##Load individual files - need to know when males turned bright if they did
ind16 = read.csv(here::here("Input files","Individuals2016 7_13_20.csv"))
ind17 = read.csv(here::here("Input files","Individuals2017 7_13_20.csv"))
ind18 = read.csv(here::here("Input files","Individuals2018 7_13_20.csv"))
ind19 = read.csv(here::here("Input files","Individuals2019 7_13_20.csv"))

##Load in birdlist files - need to know what non-breeding social group each male was in
#Birdlist files have summary information for non-breeding networks, also has FW# data
#birdlistnoKsubc is the list of individuals that were placed into social groups and were not in a social group on their own
birdlist16 = read.csv(here::here("Input files","birdlist16noKsubc.csv"))
birdlist17 = read.csv(here::here("Input files","birdlist17noKsubc.csv"))
birdlist18 = read.csv(here::here("Input files","birdlist18noKsubc.csv"))
birdlist19 = read.csv(here::here("Input files","birdlist19noKsubc.csv"))

##Load age and sex data
agesex16 = read.csv(here::here("Input files","Ages and Sexes 2016.csv"))
agesex17 = read.csv(here::here("Input files","Ages and Sexes 2017.csv"))
agesex18 = read.csv(here::here("Input files","Ages and Sexes 2018.csv"))
agesex19 = read.csv(here::here("Input files","Ages and Sexes 2019.csv"))

##Load in breeding groups for upcoming breeding seasons
b_groups16 = read.csv(here::here("Input files","Breeding Groups 2016_2016colors_short.csv"),stringsAsFactors = F)
b_groups16$Date.Created = as.Date(b_groups16$Date.Created,"%m/%d/%y")
b_groups17 = read.csv(here::here("Input files","Breeding Groups 2017_2017colors_short.csv"),stringsAsFactors = F)
b_groups17$Date.Created = as.Date(b_groups17$Date.Created,"%m/%d/%y")
b_groups18 = read.csv(here::here("Input files","Breeding Groups 2018_2018colors_short.csv"),stringsAsFactors = F)
b_groups18$Date.Created = as.Date(b_groups18$Date.Created,"%m/%d/%y")
b_groups19 = read.csv(here::here("Input files","Breeding Groups 2019_2019colors_short.csv"),stringsAsFactors = F)
b_groups19$Date.Created = as.Date(b_groups19$Date.Created,"%m/%d/%y")

##Load in all nests for nestling breeding groups
#This file uses FWnumbers instead of colors - so use FWnumbers to match
all_nests = read.csv(here::here("Input files","All RBFW Nests_social_14_19.csv"))

##Load in list of implant males
implantms = read.csv(here::here("Input files","Implant males.csv"))

##Load list of birds in removal communities in 2017 season and their FWnumbers
removals = read.csv(here::here("Input files","RemovalorControl 2017.csv"))


##Get a list of young males in each year
library(dplyr)
blym16 = birdlist16 %>% filter(Sex=="M",Age.Exact=="exact",(Current.Age==1 | Current.Age==2))
blym17 = birdlist17 %>% filter(Sex=="M",Age.Exact=="exact",(Current.Age==1 | Current.Age==2))
blym18 = birdlist18 %>% filter(Sex=="M",Age.Exact=="exact",(Current.Age==1 | Current.Age==2))
blym19 = birdlist19 %>% filter(Sex=="M",Age.Exact=="exact",(Current.Age==1 | Current.Age==2))

#Make molt numeric
ind16$Molt = as.numeric(as.character(ind16$Molt))
ind17$Molt = as.numeric(as.character(ind17$Molt))
ind18$Molt = as.numeric(as.character(ind18$Molt))
ind19$Molt = as.numeric(as.character(ind19$Molt))

#Format dates
ind16$Date = as.Date(ind16$Date,"%m/%d/%y")
ind17$Date = as.Date(ind17$Date,"%m/%d/%y")
ind18$Date = as.Date(ind18$Date,"%m/%d/%y")
ind19$Date = as.Date(ind19$Date,"%m/%d/%y")
bird16$Date = as.Date(bird16$Date,"%m/%d/%y")
bird17R$Date = as.Date(bird17R$Date,"%m/%d/%y")
bird18$Date = as.Date(bird18$Date,"%m/%d/%y")
bird19$Date = as.Date(bird19$Date,"%m/%d/%y")

#Remove blank molt scores
ind16 = ind16 %>% filter(!is.na(Molt))
ind17 = ind17 %>% filter(!is.na(Molt))
ind18 = ind18 %>% filter(!is.na(Molt))
ind19 = ind19 %>% filter(!is.na(Molt))


##Remove T-implant males from blym lists - Males from removal communities already removed because not in birdlist17
#T-implants occurred in 2018 and 2019
implant18 = implantms %>% filter(Year==2018)
implant19 = implantms %>% filter(Year==2019)
blym18 = blym18 %>% filter(!FWnumber %in% implant18$Fwnumber)
blym19 = blym19 %>% filter(!FWnumber %in% implant19$Fwnumber)




####Look at sight freq and degree relationship before young males molted####


#Run a function separately for each non-breeding network, but for each network, for 1 and 2 year old males, 
#want to know what their social connections were before they got to intermediate plumage. 
#So for each network need to first check that we saw the young males enough before they got to intermediate plumage.

#First function - get a list of young males. For each young male, get intermediate plumage date if has one. If has one, cut 
#bird file at that date - get dates before that date then generate network and calculate number of times seen 
#and degree before molt. Save to dataframe. Once I have those, I can calculate the relationship between 
#the number of times a young male was seen before molt and its degree. If young males were only seen a few times
#before they molted doesn't make sense to include them in analyses of the importance of social environment on molt date.



library(asnipe)
library(igraph)

#testing
blym = blym16
ind = ind16
bird = bird16
i=3

sight_freq_deg_before_molt = function(blym,ind,bird) {
  blym$SightFreq.bfm = NA
  blym$deg.bfm = NA
  #Generate initial network matrix to use for birds that did not molt
  bird.i = bird %>% select(Bird,Sighting) #Select columns to calculate gbi
  gbi.i = get_group_by_individual(bird.i, data_format = "individuals") #calculate gbi
  network.i = get_network(gbi.i, data_format= "GBI", association_index = "SRI") #calculate network
  sight.freq.i = data.frame(colSums(gbi.i)) #calculate sighting frequency 
  sight.freq.i$Bird = rownames(sight.freq.i) #calculate sighting frequency 
  net.i = graph.adjacency(network.i, mode="undirected", diag=FALSE, weighted=TRUE) #convert to igraph object to get degree
  sight.freq.i$deg = degree(net.i) #calculate degree
  
  for (i in 1:nrow(blym)) {
    ind.male = ind[ind$Bird %in% blym$Bird[i],] #Get sightings of male of interest
    ind.male.b = ind.male %>% filter(ind.male$Molt>33) #Get sightings where that male was bright
    if(nrow(ind.male.b)!=0) { #If the male was ever seen as bright
    bright.date = ind.male.b[order(ind.male.b$Date),]$Date[1] #Order by date and get first bright date
    bird.m = bird %>% filter(Bird %in% blym$Bird[i]) #Get social observations of bird of interest
    first.obs.date = bird.m[order(bird.m$Date),]$Date[1] #Get first date of social observations that included that bird
    
    if(bright.date<=first.obs.date) { #If molt date was before the first observation date...
      blym$SightFreq.bfm[i]=NA #Make sight freq before molt = NA
      blym$deg.bfm[i]=NA} #and deg before molt = NA
   
     else { #if there were social observations of the male before he was bright, then...
      bird.d = bird %>% filter(Date<bright.date) #Subset bird dataframe to those dates
      bird.d = bird.d %>% select(Bird,Sighting) #Select columns to calculate gbi
      gbi = get_group_by_individual(bird.d, data_format = "individuals") #calculate gbi
      network = get_network(gbi, data_format= "GBI", association_index = "SRI") #calculate network
      sight.freq = data.frame(colSums(gbi)) #calculate sighting frequency 
      sight.freq$Bird = rownames(sight.freq) #calculate sighting frequency 
      net = graph.adjacency(network, mode="undirected", diag=FALSE, weighted=TRUE) #convert to igraph object to get degree
      sight.freq$deg = degree(net) #calculate degree
      #first() turns it into a number rather than a dataframe - only one value here so this works. In later function first brought over
      #the entire column, not just the first value...
      blym$SightFreq.bfm[i] = sight.freq %>% filter(Bird==blym$Bird[i]) %>% select(colSums.gbi.) %>% first() #Add sight freq to dataframe
      blym$deg.bfm[i] = sight.freq %>% filter(Bird==blym$Bird[i]) %>% select(deg) %>% first() #Add degree to dataframe
    }} 
    
    else { #if the male was never seen as bright, use the iniital calculations 
      blym$SightFreq.bfm[i] = sight.freq.i %>% filter(Bird==blym$Bird[i]) %>% select(colSums.gbi.i.) %>% first() #Add sight freq to dataframe
      blym$deg.bfm[i] = sight.freq.i %>% filter(Bird==blym$Bird[i]) %>% select(deg) %>% first() #Add degree to dataframe
    }
  }
  return(blym)
}

##Run the function to get sight frequency and degree for young males
blym16.sfd = sight_freq_deg_before_molt(blym=blym16,ind=ind16,bird=bird16)
blym17.sfd = sight_freq_deg_before_molt(blym=blym17,ind=ind17,bird=bird17R)
blym18.sfd = sight_freq_deg_before_molt(blym=blym18,ind=ind18,bird=bird18)
blym19.sfd = sight_freq_deg_before_molt(blym=blym19,ind=ind19,bird=bird19)


##For each year, check the relationship between sight freq and degree before a male made it to intermediate plumage.
#2016
plot(blym16.sfd$SightFreq.bfm,blym16.sfd$deg.bfm,xlim=c(0,80),ylim=c(0,30),pch=19,col=blym16.sfd$Current.Age)
cor.test(blym16.sfd$SightFreq.bfm,blym16.sfd$deg.bfm)
#There is a slight correlation but these males were all seen quite a few times 

#2017
plot(blym17.sfd$SightFreq.bfm,blym17.sfd$deg.bfm,xlim=c(0,220),ylim=c(0,55),pch=19,col=blym17.sfd$Current.Age)
cor.test(blym17.sfd$SightFreq.bfm,blym17.sfd$deg.bfm)
#These males were all seen many times

#2018
plot(blym18.sfd$SightFreq.bfm,blym18.sfd$deg.bfm,xlim=c(0,180),ylim=c(0,30),pch=19,col=blym17.sfd$Current.Age)
cor.test(blym18.sfd$SightFreq.bfm,blym18.sfd$deg.bfm)
#Lowest-sighted male is RRY - only seen 12 times before he made it to 33% - remove him from analyses, that's not seen very many times
#Second lowest-sighted male is LZI who we don't have breeding data for, so he will be removed from analyses, 
#he actually disappeared (maybe dispersed) during the non-breeding season
blym18.sfd = blym18.sfd %>% filter(Bird!="RRY")
plot(blym18.sfd$SightFreq.bfm,blym18.sfd$deg.bfm,xlim=c(0,180),ylim=c(0,30),pch=19,col=blym17.sfd$Current.Age)
cor.test(blym18.sfd$SightFreq.bfm,blym18.sfd$deg.bfm)

#2019
plot(blym19.sfd$SightFreq.bfm,blym19.sfd$deg.bfm,xlim=c(0,240),ylim=c(0,45),pch=19,,col=blym17.sfd$Current.Age)
cor.test(blym19.sfd$SightFreq.bfm,blym19.sfd$deg.bfm)
#Slight correlation but most of these birds seen quite a few times, no obvious relationship






####Measure each young male's observed non-breeding social connections to different social classes before he molted#### 

##Generate yearly networks and create an agesex file that matches the birds in the network
#Function to go from gbi to network
gbi_to_network = function(bird) {
  bird = bird %>% select(Bird,Sighting) #Select columns to calculate gbi
  gbi = get_group_by_individual(bird, data_format = "individuals") #calculate gbi
  network = get_network(gbi, data_format= "GBI", association_index = "SRI") #calculate network
  network = network[order(rownames(network)), order(colnames(network))]
  return(network)
}

#Get networks
network16 = gbi_to_network(bird16)
network16 = network16[order(rownames(network16)), order(colnames(network16))]
network17 = gbi_to_network(bird17R)
network17 = network17[order(rownames(network17)), order(colnames(network17))]
network18 = gbi_to_network(bird18)
network18 = network18[order(rownames(network18)), order(colnames(network18))]
network19 = gbi_to_network(bird19)
network19 = network19[order(rownames(network19)), order(colnames(network19))]


#For 2016 - need to add BN1 and BN2 to agesex 16 - two old unbanded bright males (molted early) we never caught
#Call them three year olds since they were bright early
agesex16.add = as.data.frame(matrix(nrow=2,ncol=ncol(agesex16)))
colnames(agesex16.add) <- colnames(agesex16)
agesex16.add$Bird = c("BN1","BN2")
agesex16.add$Sex = c("M","M")
agesex16.add$Current.Age = c(3,3)
agesex16.add$Age.Exact. = c("min","min")
agesex16 = rbind(agesex16,agesex16.add)

#Get agesex that matches network order
agesex16_net = agesex16 %>% filter(Bird %in% rownames(network16))
agesex16_net$Bird = as.character(agesex16_net$Bird)
agesex16_net = agesex16_net[order(agesex16_net$Bird),]
agesex17_net = agesex17 %>% filter(Bird %in% rownames(network17))
agesex17_net$Bird = as.character(agesex17_net$Bird)
agesex17_net = agesex17_net[order(agesex17_net$Bird),]
agesex18_net = agesex18 %>% filter(Bird %in% rownames(network18))
agesex18_net$Bird = as.character(agesex18_net$Bird)
agesex18_net = agesex18_net[order(agesex18_net$Bird),]
agesex19_net = agesex19 %>% filter(Bird %in% rownames(network19))
agesex19_net$Bird = as.character(agesex19_net$Bird)
agesex19_net = agesex19_net[order(agesex19_net$Bird),]

#Merge agesex and birdlist (groups) to add non-breeding social groups to agesex file
groups16 = birdlist16 %>% select(Bird,Social.Group)
groups17 = birdlist17 %>% select(Bird,Social.Group)
groups18 = birdlist18 %>% select(Bird,Social.Group)
groups19 = birdlist19 %>% select(Bird,Social.Group)
agesex16_net = merge(agesex16_net,groups16,by="Bird",all.x=T)
agesex17_net = merge(agesex17_net,groups17,by="Bird",all.x=T)
agesex18_net = merge(agesex18_net,groups18,by="Bird",all.x=T)
agesex19_net = merge(agesex19_net,groups19,by="Bird",all.x=T)

#Convert birds without social groups to social group 0 to prevent warnings in function
#Social group 0 does not exist in the non-breeding network and all of the males of interest in blym.sfd have a social group 
agesex16_net[is.na(agesex16_net$Social.Group),]$Social.Group <- 0
agesex17_net[is.na(agesex17_net$Social.Group),]$Social.Group <- 0
agesex18_net[is.na(agesex18_net$Social.Group),]$Social.Group <- 0
agesex19_net[is.na(agesex19_net$Social.Group),]$Social.Group <- 0

#Check that agesex_net order = network order
agesex16_net$Bird == rownames(network16)
agesex17_net$Bird == rownames(network17)
agesex18_net$Bird == rownames(network18)
agesex19_net$Bird == rownames(network19)





####One-year-old male weighted degree function####
#Function to get weighted degree connections to social classes before molt
#testing
blym.sfd = blym16.sfd
bird = bird16
ind = ind16
agesex_net = agesex16_net
b_groups = b_groups16
network_full = network16
i=4

connections_before_bright_1yo_wdeg = function(blym.sfd,ind,bird,agesex_net,b_groups,network_full) {
  blym.sfd$wdeg_1yo_males = NA #Total weighted degree to all 1yo males 
  blym.sfd$wdeg_1yo_males_samegrp = NA #Total weighted degree to 1yo males in same social group
  blym.sfd$wdeg_1yo_males_diffgrps = NA #Total weighted degree to 1yo males in different social groups
  blym.sfd$wdeg_old_males = NA #Total weighted degree to all old males(2+)
  blym.sfd$wdeg_old_males_samegrp = NA #Total weighted degree to old males (2+) in same social group (father and brother)
  blym.sfd$wdeg_oldest_male_samegrp = NA #Total weighted degree to oldest male(s) in the same group - probably the social father, if multiple males of same old age, use both
  blym.sfd$wdeg_old_males_diffgrps = NA #Total weighted degree to old males (2+) in different social groups
  blym.sfd$wdeg_females = NA #Total weighted degree to all females
  blym.sfd$wdeg_females_diffgrps = NA #Total weighted degree to females in different social groups
  blym.sfd$wdeg_female_paired = NA #Weighted degree of connection to female the male eventually paired with
  blym.sfd$wdeg_all = NA #Total weighted degree
  
  #Get just one year olds
  blym.sfd = blym.sfd %>% filter(Current.Age==1)
  
  #Start loops
  for (i in 1:nrow(blym.sfd)) {
    if(is.na(blym.sfd$SightFreq.bfm[i])) {} #Do nothing if is NA - means he molted before social observations and can't measure connections before molt - from prev function
    else { #otherwise...
    male = blym.sfd$Bird[i]
    ind.male = ind[ind$Bird %in% male,] #Get sightings of male of interest
    ind.male.b = ind.male %>% filter(ind.male$Molt>33) #Get sightings where that male was intermediate
    
    #Get network before male of interest was bright
    if(nrow(ind.male.b)!=0) { #If the male was ever seen as bright
      bright.date = ind.male.b[order(ind.male.b$Date),]$Date[1] #Order by date and get first bright date
      bird.m = bird %>% filter(Bird %in% male) #Get social observations of bird of interest
      bird.d = bird %>% filter(Date<bright.date) #Subset bird dataframe to dates before first seen as bright
      bird.d = bird.d %>% select(Bird,Sighting) #Select columns to calculate gbi
      gbi = get_group_by_individual(bird.d, data_format = "individuals") #calculate gbi
      network = get_network(gbi, data_format= "GBI", association_index = "SRI") #calculate network
      network = network[order(rownames(network)), order(colnames(network))] #Order network after it's been calculated with subsetted gbi
      
      #Get agesex_net that matches network order
      agesex_net_no = agesex_net %>% filter(Bird %in% rownames(network))
      agesex_net_no = agesex_net_no[order(agesex_net_no$Bird),]
      
      #Check that agesex_net_no matches network order
      if(all(agesex_net_no$Bird==rownames(network))) {} else {stop("agesex_net_no and network not equal")} 
      
      #Get row number of male of interest in the network/agesex_net_no so I can pull out the male's value in the rowsums step below
      h = as.numeric(rownames(agesex_net_no[agesex_net_no$Bird==male,]))
      
      #Get male's social group
      male.social.group = agesex_net_no[h,]$Social.Group
      
      #Calculate connections
      #rowsum calculates the sum of the false [1] and true [2] conditions of the columns that match the characteristics
      #Need agesex_net_no and network to be in same order
      wdeg_1yo_males = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1))[,h] #[,h] gets the row for the male of interest
      wdeg_1yo_males_samegrp = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1 & agesex_net_no$Social.Group==male.social.group))[,h]
      wdeg_1yo_males_diffgrps = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1 & agesex_net_no$Social.Group!=male.social.group))[,h]
      wdeg_old_males = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2))[,h]
      wdeg_old_males_samegrp = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2 & agesex_net_no$Social.Group==male.social.group))[,h]
      oldest_male_in_group = agesex_net_no %>% filter(Sex=="M",Current.Age>=2,Social.Group==male.social.group,Bird!=male) #Get old males in group that are not the focal male
      oldest_male_in_group = oldest_male_in_group %>% arrange(desc(Current.Age)) %>% filter(Current.Age==oldest_male_in_group$Current.Age[1]) #Get the oldest male in the group or if multiple males of same old age, get both males
      if(nrow(oldest_male_in_group)>0) {wdeg_oldest_male_samegrp = rowsum(network,(agesex_net_no$Bird %in% oldest_male_in_group$Bird))[,h]} else 
      {wdeg_oldest_male_samegrp = NA} #If no oldest male other than focal male, make = NA
      wdeg_old_males_diffgrps = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2 & agesex_net_no$Social.Group!=male.social.group))[,h]
      wdeg_females = rowsum(network,(agesex_net_no$Sex=="F"))[,h]
      wdeg_females_diffgrps = rowsum(network,(agesex_net_no$Sex=="F" & agesex_net_no$Social.Group!=male.social.group))[,h]
      female_paired_groups = b_groups %>% filter(b_groups$Male %in% male) #Get breeding groups male was breeding male in
      female_paired = female_paired_groups %>% arrange(Date.Created) %>% select(Female) #Take first group. Adding first() onto here did not work, don't know why...
      if(nrow(female_paired)>0) {wdeg_female_paired = rowsum(network,(agesex_net_no$Bird==female_paired[1,]))[,h]} else #rowsum for values that equal that first female
      {wdeg_female_paired = NA} #If no paired female in group, make = NA
      wdeg_all = rowSums(network)[h]

      #Put values into dataframe
      blym.sfd$wdeg_1yo_males[i] = wdeg_1yo_males[2] #2 is true if it's there. There should always be false. 2 = NA if not there
      blym.sfd$wdeg_1yo_males_samegrp[i] = wdeg_1yo_males_samegrp[2]
      blym.sfd$wdeg_1yo_males_diffgrps[i] = wdeg_1yo_males_diffgrps[2]
      blym.sfd$wdeg_old_males[i] = wdeg_old_males[2]
      blym.sfd$wdeg_old_males_samegrp[i]  = wdeg_old_males_samegrp[2]
      blym.sfd$wdeg_oldest_male_samegrp[i] = wdeg_oldest_male_samegrp[2]
      blym.sfd$wdeg_old_males_diffgrps[i] = wdeg_old_males_diffgrps[2]
      blym.sfd$wdeg_females[i] = wdeg_females[2]
      blym.sfd$wdeg_females_diffgrps[i] = wdeg_females_diffgrps[2]
      blym.sfd$wdeg_female_paired[i] = wdeg_female_paired[2]
      blym.sfd$wdeg_all[i] = wdeg_all

    } else {  #Do same for males that did not turn bright during non-breeding social observations - don't have to recalculate network here
      
      #Get agesex_net that matches network_full order
      agesex_net_no = agesex_net %>% filter(Bird %in% rownames(network_full))
      agesex_net_no = agesex_net_no[order(agesex_net_no$Bird),]
      
      #Check that agesex_net_no matches network_full order
      if(all(agesex_net_no$Bird!=rownames(network_full))) {stop("agesex_net_no and network_full not equal")} 
      
      #Get row number of male of interest in the network_full/agesex_net_no so I can pull out the male's value in the rowsums step below
      h = as.numeric(rownames(agesex_net_no[agesex_net_no$Bird==male,]))
      
      #Get male's social group
      male.social.group = agesex_net_no[h,]$Social.Group
      
      #Calculate connections
      #rowsum calculates the sum of the false [1] and true [2] conditions of the columns that match the characteristics
      #Need agesex_net_no and network_full to be in same order
      wdeg_1yo_males = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1))[,h] #[,h] gets the row for the male of interest
      wdeg_1yo_males_samegrp = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1 & agesex_net_no$Social.Group==male.social.group))[,h]
      wdeg_1yo_males_diffgrps = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1 & agesex_net_no$Social.Group!=male.social.group))[,h]
      wdeg_old_males = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2))[,h]
      wdeg_old_males_samegrp = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2 & agesex_net_no$Social.Group==male.social.group))[,h]
      oldest_male_in_group = agesex_net_no %>% filter(Sex=="M",Current.Age>=2,Social.Group==male.social.group,Bird!=male) #Get old males in group that are not the focal male
      oldest_male_in_group = oldest_male_in_group %>% arrange(desc(Current.Age)) %>% filter(Current.Age==oldest_male_in_group$Current.Age[1]) #Get the oldest male in the group or if multiple males of same old age, get both males
      if(nrow(oldest_male_in_group)>0) {wdeg_oldest_male_samegrp = rowsum(network_full,(agesex_net_no$Bird %in% oldest_male_in_group$Bird))[,h]} else 
      {wdeg_oldest_male_samegrp = NA} #If no oldest male other than focal male, make = NA
      wdeg_old_males_diffgrps = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2 & agesex_net_no$Social.Group!=male.social.group))[,h]
      wdeg_females = rowsum(network_full,(agesex_net_no$Sex=="F"))[,h]
      wdeg_females_diffgrps = rowsum(network_full,(agesex_net_no$Sex=="F" & agesex_net_no$Social.Group!=male.social.group))[,h]
      female_paired_groups = b_groups %>% filter(b_groups$Male %in% male) #Get breeding groups male was breeding male in
      female_paired = female_paired_groups %>% arrange(Date.Created) %>% select(Female) #Take first group. Adding first() onto here did not work, don't know why...
      if(nrow(female_paired)>0) {wdeg_female_paired = rowsum(network_full,(agesex_net_no$Bird==female_paired[1,]))[,h]} else #rowsum for values that equal that first female
      {wdeg_female_paired = NA} #If no paired female in group, make = NA
      wdeg_all = rowSums(network_full)[h]
      
      #Put values into dataframe
      blym.sfd$wdeg_1yo_males[i] = wdeg_1yo_males[2] #2 is true if it's there. There should always be false. 2 = NA if not there
      blym.sfd$wdeg_1yo_males_samegrp[i] = wdeg_1yo_males_samegrp[2]
      blym.sfd$wdeg_1yo_males_diffgrps[i] = wdeg_1yo_males_diffgrps[2]
      blym.sfd$wdeg_old_males[i] = wdeg_old_males[2]
      blym.sfd$wdeg_old_males_samegrp[i]  = wdeg_old_males_samegrp[2]
      blym.sfd$wdeg_oldest_male_samegrp[i] = wdeg_oldest_male_samegrp[2]
      blym.sfd$wdeg_old_males_diffgrps[i] = wdeg_old_males_diffgrps[2]
      blym.sfd$wdeg_females[i] = wdeg_females[2]
      blym.sfd$wdeg_females_diffgrps[i] = wdeg_females_diffgrps[2]
      blym.sfd$wdeg_female_paired[i] = wdeg_female_paired[2]
      blym.sfd$wdeg_all[i] = wdeg_all
    }
    } 
  }
  return(blym.sfd)
}

#Calculate one year old male connections using weighted degree
blym.sfd.1yo.wdeg16=connections_before_bright_1yo_wdeg(blym.sfd=blym16.sfd,ind=ind16,bird=bird16,agesex_net=agesex16_net,
                                                     b_groups=b_groups16,network_full=network16)
blym.sfd.1yo.wdeg17=connections_before_bright_1yo_wdeg(blym.sfd=blym17.sfd,ind=ind17,bird=bird17,agesex_net=agesex17_net,
                                                     b_groups=b_groups17,network_full=network17)
blym.sfd.1yo.wdeg18=connections_before_bright_1yo_wdeg(blym.sfd=blym18.sfd,ind=ind18,bird=bird18,agesex_net=agesex18_net,
                                                     b_groups=b_groups18,network_full=network18)
blym.sfd.1yo.wdeg19=connections_before_bright_1yo_wdeg(blym.sfd=blym19.sfd,ind=ind19,bird=bird19,agesex_net=agesex19_net,
                                                     b_groups=b_groups19,network_full=network19)
#In 2018 VVI was not in any social observations before he made it to intermediate plumage




####One-year-old male binary degree function####
#Function to get binary degree connections to social classes before molt
#Same as above, just converting network to binary degree early in function

connections_before_bright_1yo_bdeg = function(blym.sfd,ind,bird,agesex_net,b_groups,network_full) {
  blym.sfd$bdeg_1yo_males = NA #Total binary degree to all 1yo males 
  blym.sfd$bdeg_1yo_males_samegrp = NA #Total binary degree to 1yo males in same social group
  blym.sfd$bdeg_1yo_males_diffgrps = NA #Total binary degree to 1yo males in different social groups
  blym.sfd$bdeg_old_males = NA #Total binary degree to all old males(2+)
  blym.sfd$bdeg_old_males_samegrp = NA #Total binary degree to old males (2+) in same social group (father and brother)
  blym.sfd$bdeg_oldest_male_samegrp = NA #Total binary degree to oldest male(s) in the same group - probably the social father, if multiple males of same old age, use both
  blym.sfd$bdeg_old_males_diffgrps = NA #Total binary degree to old males (2+) in different social groups
  blym.sfd$bdeg_females = NA #Total binary degree to all females
  blym.sfd$bdeg_females_diffgrps = NA #Total binary degree to females in different social groups
  blym.sfd$bdeg_female_paired = NA #binary degree of connection to female the male eventually paired with
  
  #Get just one year olds
  blym.sfd = blym.sfd %>% filter(Current.Age==1)
  
  #Make the network_full binary
  network_full[network_full>0] <- 1 #Make the network binary
  
  #Start loops
  for (i in 1:nrow(blym.sfd)) {
    if(is.na(blym.sfd$SightFreq.bfm[i])) {} #Do nothing if is NA - means he molted before social observations and can't measure connections before molt - from prev function
    else { #otherwise...
      male = blym.sfd$Bird[i]
      ind.male = ind[ind$Bird %in% male,] #Get sightings of male of interest
      ind.male.b = ind.male %>% filter(ind.male$Molt>33) #Get sightings where that male was intermediate
      
      #Get network before male of interest was bright
      if(nrow(ind.male.b)!=0) { #If the male was ever seen as bright
        bright.date = ind.male.b[order(ind.male.b$Date),]$Date[1] #Order by date and get first bright date
        bird.m = bird %>% filter(Bird %in% male) #Get social observations of bird of interest
        bird.d = bird %>% filter(Date<bright.date) #Subset bird dataframe to dates before first seen as bright
        bird.d = bird.d %>% select(Bird,Sighting) #Select columns to calculate gbi
        gbi = get_group_by_individual(bird.d, data_format = "individuals") #calculate gbi
        network = get_network(gbi, data_format= "GBI", association_index = "SRI") #calculate network
        network[network>0] <- 1 #Make the network binary
        network = network[order(rownames(network)), order(colnames(network))] #Order network after it's been calculated with subsetted gbi
        
        
        #Get agesex_net that matches network order
        agesex_net_no = agesex_net %>% filter(Bird %in% rownames(network))
        agesex_net_no = agesex_net_no[order(agesex_net_no$Bird),]
        
        #Check that agesex_net_no matches network order
        if(all(agesex_net_no$Bird!=rownames(network))) {stop("agesex_net_no and network not equal")} 
        
        #Get row number of male of interest in the network/agesex_net_no so I can pull out the male's value in the rowsums step below
        h = as.numeric(rownames(agesex_net_no[agesex_net_no$Bird==male,]))
        
        #Get male's social group
        male.social.group = agesex_net_no[h,]$Social.Group
        
        #Calculate connections
        #rowsum calculates the sum of the false [1] and true [2] conditions of the columns that match the characteristics
        #Need agesex_net_no and network to be in same order
        bdeg_1yo_males = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1))[,h] #[,h] gets the row for the male of interest
        bdeg_1yo_males_samegrp = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1 & agesex_net_no$Social.Group==male.social.group))[,h]
        bdeg_1yo_males_diffgrps = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1 & agesex_net_no$Social.Group!=male.social.group))[,h]
        bdeg_old_males = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2))[,h]
        bdeg_old_males_samegrp = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2 & agesex_net_no$Social.Group==male.social.group))[,h]
        oldest_male_in_group = agesex_net_no %>% filter(Sex=="M",Current.Age>=2,Social.Group==male.social.group,Bird!=male) #Get old males in group that are not the focal male
        oldest_male_in_group = oldest_male_in_group %>% arrange(desc(Current.Age)) %>% filter(Current.Age==oldest_male_in_group$Current.Age[1]) #Get the oldest male in the group or if multiple males of same old age, get both males
        if(nrow(oldest_male_in_group)>0) {bdeg_oldest_male_samegrp = rowsum(network,(agesex_net_no$Bird %in% oldest_male_in_group$Bird))[,h]} else 
        {bdeg_oldest_male_samegrp = NA} #If no oldest male other than focal male, make = NA
        bdeg_old_males_diffgrps = rowsum(network,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2 & agesex_net_no$Social.Group!=male.social.group))[,h]
        bdeg_females = rowsum(network,(agesex_net_no$Sex=="F"))[,h]
        bdeg_females_diffgrps = rowsum(network,(agesex_net_no$Sex=="F" & agesex_net_no$Social.Group!=male.social.group))[,h]
        female_paired_groups = b_groups %>% filter(b_groups$Male %in% male) #Get breeding groups male was breeding male in
        female_paired = female_paired_groups %>% arrange(Date.Created) %>% select(Female) #Take first group. Adding first() onto here did not work, don't know why...
        if(nrow(female_paired)>0) {bdeg_female_paired = rowsum(network,(agesex_net_no$Bird==female_paired[1,]))[,h]} else #rowsum for values that equal that first female
        {bdeg_female_paired = NA} #If no paired female in group, make = NA
        
        #Put values into dataframe
        blym.sfd$bdeg_1yo_males[i] = bdeg_1yo_males[2] #2 is true if it's there. There should always be false. 2 = NA if not there
        blym.sfd$bdeg_1yo_males_samegrp[i] = bdeg_1yo_males_samegrp[2]
        blym.sfd$bdeg_1yo_males_diffgrps[i] = bdeg_1yo_males_diffgrps[2]
        blym.sfd$bdeg_old_males[i] = bdeg_old_males[2]
        blym.sfd$bdeg_old_males_samegrp[i]  = bdeg_old_males_samegrp[2]
        blym.sfd$bdeg_oldest_male_samegrp[i] = bdeg_oldest_male_samegrp[2]
        blym.sfd$bdeg_old_males_diffgrps[i] = bdeg_old_males_diffgrps[2]
        blym.sfd$bdeg_females[i] = bdeg_females[2]
        blym.sfd$bdeg_females_diffgrps[i] = bdeg_females_diffgrps[2]
        blym.sfd$bdeg_female_paired[i] = bdeg_female_paired[2]
        
        
      } else {  #Do same for males that did not turn bright during non-breeding social observations - don't have to recalculate network here
        
        #Get agesex_net that matches network_full order
        agesex_net_no = agesex_net %>% filter(Bird %in% rownames(network_full))
        agesex_net_no = agesex_net_no[order(agesex_net_no$Bird),]
        
        #Check that agesex_net_no matches network_full order
        if(all(agesex_net_no$Bird==rownames(network_full))) {} else {stop("agesex_net_no and network_full not equal")} 
        
        #Get row number of male of interest in the network_full/agesex_net_no so I can pull out the male's value in the rowsums step below
        h = as.numeric(rownames(agesex_net_no[agesex_net_no$Bird==male,]))
        
        #Get male's social group
        male.social.group = agesex_net_no[h,]$Social.Group
        
        #Calculate connections
        #rowsum calculates the sum of the false [1] and true [2] conditions of the columns that match the characteristics
        #Need agesex_net_no and network_full to be in same order
        bdeg_1yo_males = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1))[,h] #[,h] gets the row for the male of interest
        bdeg_1yo_males_samegrp = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1 & agesex_net_no$Social.Group==male.social.group))[,h]
        bdeg_1yo_males_diffgrps = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age==1 & agesex_net_no$Social.Group!=male.social.group))[,h]
        bdeg_old_males = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2))[,h]
        bdeg_old_males_samegrp = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2 & agesex_net_no$Social.Group==male.social.group))[,h]
        oldest_male_in_group = agesex_net_no %>% filter(Sex=="M",Current.Age>=2,Social.Group==male.social.group,Bird!=male) #Get old males in group that are not the focal male
        oldest_male_in_group = oldest_male_in_group %>% arrange(desc(Current.Age)) %>% filter(Current.Age==oldest_male_in_group$Current.Age[1]) #Get the oldest male in the group or if multiple males of same old age, get both males
        if(nrow(oldest_male_in_group)>0) {bdeg_oldest_male_samegrp = rowsum(network_full,(agesex_net_no$Bird %in% oldest_male_in_group$Bird))[,h]} else 
        {bdeg_oldest_male_samegrp = NA} #If no oldest male other than focal male, make = NA
        bdeg_old_males_diffgrps = rowsum(network_full,(agesex_net_no$Sex=="M" & agesex_net_no$Current.Age>=2 & agesex_net_no$Social.Group!=male.social.group))[,h]
        bdeg_females = rowsum(network_full,(agesex_net_no$Sex=="F"))[,h]
        bdeg_females_diffgrps = rowsum(network_full,(agesex_net_no$Sex=="F" & agesex_net_no$Social.Group!=male.social.group))[,h]
        female_paired_groups = b_groups %>% filter(b_groups$Male %in% male) #Get breeding groups male was breeding male in
        female_paired = female_paired_groups %>% arrange(Date.Created) %>% select(Female) #Take first group. Adding first() onto here did not work, don't know why...
        if(nrow(female_paired)>0) {bdeg_female_paired = rowsum(network_full,(agesex_net_no$Bird==female_paired[1,]))[,h]} else #rowsum for values that equal that first female
        {bdeg_female_paired = NA} #If no paired female in group, make = NA
        
        #Put values into dataframe
        blym.sfd$bdeg_1yo_males[i] = bdeg_1yo_males[2] #2 is true if it's there. There should always be false. 2 = NA if not there
        blym.sfd$bdeg_1yo_males_samegrp[i] = bdeg_1yo_males_samegrp[2]
        blym.sfd$bdeg_1yo_males_diffgrps[i] = bdeg_1yo_males_diffgrps[2]
        blym.sfd$bdeg_old_males[i] = bdeg_old_males[2]
        blym.sfd$bdeg_old_males_samegrp[i]  = bdeg_old_males_samegrp[2]
        blym.sfd$bdeg_oldest_male_samegrp[i] = bdeg_oldest_male_samegrp[2]
        blym.sfd$bdeg_old_males_diffgrps[i] = bdeg_old_males_diffgrps[2]
        blym.sfd$bdeg_females[i] = bdeg_females[2]
        blym.sfd$bdeg_females_diffgrps[i] = bdeg_females_diffgrps[2]
        blym.sfd$bdeg_female_paired[i] = bdeg_female_paired[2]
      }
    } 
  }
  return(blym.sfd)
}

#Calculate one year old male connections using binary degree
blym.sfd.1yo.bdeg16=connections_before_bright_1yo_bdeg(blym.sfd=blym16.sfd,ind=ind16,bird=bird16,agesex_net=agesex16_net,
                                                       b_groups=b_groups16,network_full=network16)
blym.sfd.1yo.bdeg17=connections_before_bright_1yo_bdeg(blym.sfd=blym17.sfd,ind=ind17,bird=bird17,agesex_net=agesex17_net,
                                                       b_groups=b_groups17,network_full=network17)
blym.sfd.1yo.bdeg18=connections_before_bright_1yo_bdeg(blym.sfd=blym18.sfd,ind=ind18,bird=bird18,agesex_net=agesex18_net,
                                                       b_groups=b_groups18,network_full=network18)
blym.sfd.1yo.bdeg19=connections_before_bright_1yo_bdeg(blym.sfd=blym19.sfd,ind=ind19,bird=bird19,agesex_net=agesex19_net,
                                                       b_groups=b_groups19,network_full=network19)


##Combine wdeg and bdeg files
blym.sfd.1yo.16 = cbind(blym.sfd.1yo.wdeg16, blym.sfd.1yo.bdeg16 %>% select(contains("bdeg")))
blym.sfd.1yo.17 = cbind(blym.sfd.1yo.wdeg17, blym.sfd.1yo.bdeg17 %>% select(contains("bdeg")))
blym.sfd.1yo.18 = cbind(blym.sfd.1yo.wdeg18, blym.sfd.1yo.bdeg18 %>% select(contains("bdeg")))
blym.sfd.1yo.19 = cbind(blym.sfd.1yo.wdeg19, blym.sfd.1yo.bdeg19 %>% select(contains("bdeg")))


##Remove some unnecessary columns
blym.sfd.1yo.16 = blym.sfd.1yo.16 %>% select(-Year.variable,-PrevBreedingYear,-PrevBreedingStatus,-Status.variable,-SightFreq,-deg16noK)
blym.sfd.1yo.17 = blym.sfd.1yo.17 %>% select(-Year.variable,-PrevBreedingYear,-PrevBreedingStatus,-Status.variable,-SightFreq,-deg17noK)
blym.sfd.1yo.18 = blym.sfd.1yo.18 %>% select(-Year.variable,-PrevBreedingYear,-PrevBreedingStatus,-Status.variable,-SightFreq,-deg18noK)
blym.sfd.1yo.19 = blym.sfd.1yo.19 %>% select(-Year.variable,-PrevBreedingYear,-PrevBreedingStatus,-Status.variable,-SightFreq,-deg19noK)


##Make all NA values = 0
blym.sfd.1yo.16[is.na(blym.sfd.1yo.16)]<-0
blym.sfd.1yo.17[is.na(blym.sfd.1yo.17)]<-0
blym.sfd.1yo.18[is.na(blym.sfd.1yo.18)]<-0
blym.sfd.1yo.19[is.na(blym.sfd.1yo.19)]<-0





####Determine if one-year-old males were in a non-breeding group with the male they were hatched to####

library(reshape2)

##Function 
#testing
blym.sfd.1yo = blym.sfd.1yo.17
i=3
birdlist = birdlist17

social_group_equals_nestling = function(blym.sfd.1yo,all_nests,birdlist) {
  all_nestsm = melt(all_nests,id.vars=1:6) #melt all_nests
  colnames(all_nestsm)[8] = "nestling.fwno"
  blym.sfd.1yo$Group_with_social_father = NA #Get blank column
  
  #Start loop
  for (i in 1:nrow(blym.sfd.1yo)) {
    all_nestsm_male = all_nestsm %>% filter(nestling.fwno==blym.sfd.1yo$FWnumber[i]) 
    if(nrow(all_nestsm_male)!=0) {
      soc.father = all_nestsm_male$Male.fwno
      father.soc.group = birdlist %>% filter(FWnumber==soc.father) %>% select(Social.Group)
      if(nrow(father.soc.group)!=0) {if(father.soc.group==blym.sfd.1yo$Social.Group[i]) {blym.sfd.1yo$Group_with_social_father[i]="Yes"} else 
      {blym.sfd.1yo$Group_with_social_father[i]="No"}} else {blym.sfd.1yo$Group_with_social_father[i]="No"}
    } else {
      blym.sfd.1yo$Group_with_social_father[i] = "no.data"
    }
  }
  return(blym.sfd.1yo)
}

#Run the function and add to same dataframe
blym.sfd.1yo.16 = social_group_equals_nestling(blym.sfd.1yo=blym.sfd.1yo.16,all_nests=all_nests,birdlist=birdlist16)
blym.sfd.1yo.17 = social_group_equals_nestling(blym.sfd.1yo=blym.sfd.1yo.17,all_nests=all_nests,birdlist=birdlist17)
blym.sfd.1yo.18 = social_group_equals_nestling(blym.sfd.1yo=blym.sfd.1yo.18,all_nests=all_nests,birdlist=birdlist18)
blym.sfd.1yo.19 = social_group_equals_nestling(blym.sfd.1yo=blym.sfd.1yo.19,all_nests=all_nests,birdlist=birdlist19)
#This function could generate incorrect data if any of the young males that did have natal nests had a male that was not
#banded. But very few nests in all_nests have no father and the offspring from those nests are not in these lists of 1yo males. 



####Record whether for paired 1-year-old males we knew the identify of their female####

#This will be important in looking at wdeg and bdeg to their paired females - if the female was UNB or Unknown cannot include those males
#in those tests

##Function
#testing
blym.sfd.1yo = blym.sfd.1yo.18
b_groups = b_groups18
agesex = agesex18
i=1

unknown_female = function(blym.sfd.1yo,b_groups,agesex) {
  blym.sfd.1yo$unknown_female = NA #Get a blank column
  agesex = agesex %>% filter(Bird!="") #Remove birds without band combos (nestlings) from agesex so males with blank females get caught as having an unknown female
  #Start loop
  for (i in 1:nrow(blym.sfd.1yo)) {
    Female = b_groups %>% filter(Male==blym.sfd.1yo$Bird[i]) %>% arrange(Date.Created) %>% select(Female) #Get female of male's first breeding group
    if(nrow(Female)!=0) {
    if(Female[1,] %in% agesex$Bird) {blym.sfd.1yo$unknown_female[i]="No"} else {blym.sfd.1yo$unknown_female[i]="Yes"}
    } else {blym.sfd.1yo$unknown_female[i]=NA}
    }
  return(blym.sfd.1yo)
}

#Run function and add to same dataframe
blym.sfd.1yo.16 = unknown_female(blym.sfd.1yo=blym.sfd.1yo.16,b_groups=b_groups16,agesex=agesex16)
blym.sfd.1yo.17 = unknown_female(blym.sfd.1yo=blym.sfd.1yo.17,b_groups=b_groups17,agesex=agesex17)
blym.sfd.1yo.18 = unknown_female(blym.sfd.1yo=blym.sfd.1yo.18,b_groups=b_groups18,agesex=agesex18)
blym.sfd.1yo.19 = unknown_female(blym.sfd.1yo=blym.sfd.1yo.19,b_groups=b_groups19,agesex=agesex19)



####Calculate relatedness among one-year-old males and their social father if in the same non-breeding social group####

##Load relatedness data
rel = read.csv(here::here("Input files","RelatednessEstimates_8xdenovo80_14to19.csv"))

#Get pairs in both columns so order of birds in function does not matter
rel2 = rel %>% select(Pair,ind2,ind1,wang)
colnames(rel2)[2:3] = c("ind1","ind2")
rel = rbind(rel,rel2)


##Function
#Testing
blym.sfd.1yo = blym.sfd.1yo.17
i=4

rel_to_social_father = function(blym.sfd.1yo,rel,all_nests) {
  blym.sfd.1yo$rel_to_social_father_samegrp = NA #Get a blank column
  all_nestsm = melt(all_nests,id.vars=1:6) #melt all_nests
  colnames(all_nestsm)[8] = "nestling.fwno"
  
  #Start Loop
  for (i in 1:nrow(blym.sfd.1yo)) { 
    if(blym.sfd.1yo$Group_with_social_father[i]=="no.data" | blym.sfd.1yo$Group_with_social_father[i]=="No") {
      blym.sfd.1yo$rel_to_social_father_samegrp[i] = NA} else {
        all_nestsm_male = all_nestsm %>% filter(nestling.fwno==blym.sfd.1yo$FWnumber[i]) 
        rel.males = rel %>% filter(ind1==all_nestsm_male$Male.fwno & ind2==blym.sfd.1yo$FWnumber[i])
        blym.sfd.1yo$rel_to_social_father_samegrp[i] = rel.males$wang}
  }
  return(blym.sfd.1yo)
}

#Run function
blym.sfd.1yo.16 = rel_to_social_father(blym.sfd.1yo=blym.sfd.1yo.16,rel=rel,all_nests=all_nests)
blym.sfd.1yo.17 = rel_to_social_father(blym.sfd.1yo=blym.sfd.1yo.17,rel=rel,all_nests=all_nests)
blym.sfd.1yo.18 = rel_to_social_father(blym.sfd.1yo=blym.sfd.1yo.18,rel=rel,all_nests=all_nests)
blym.sfd.1yo.19 = rel_to_social_father(blym.sfd.1yo=blym.sfd.1yo.19,rel=rel,all_nests=all_nests)




####Determine if a one-year-old male and his paired female were in the same non-breeding social group####

#Testing
blym.sfd.1yo = blym.sfd.1yo.17
b_groups = b_groups17
birdlist = birdlist17
i=2

group_with_female = function(blym.sfd.1yo,b_groups,birdlist) {
  blym.sfd.1yo$Group_with_paired_female = NA
  for (i in 1:nrow(blym.sfd.1yo)) {
    Female = b_groups %>% filter(Male==blym.sfd.1yo$Bird[i]) %>% arrange(Date.Created) %>% select(Female) #Get female of male's first breeding group
    if(nrow(Female)!=0) {Female.sg = birdlist %>% filter(Bird==Female[1,]) %>% select(Social.Group)
    if(nrow(Female.sg)!=0) {
      if(Female.sg[1,]==blym.sfd.1yo$Social.Group[i]) {blym.sfd.1yo$Group_with_paired_female[i]="Yes"} else
      {blym.sfd.1yo$Group_with_paired_female[i]="No"}} 
    else {blym.sfd.1yo$Group_with_paired_female[i]="No"}} 
    else {blym.sfd.1yo$Group_with_paired_female[i] = NA}} 
  return(blym.sfd.1yo)
}

#Run function
blym.sfd.1yo.16 = group_with_female(blym.sfd.1yo=blym.sfd.1yo.16,b_groups=b_groups16,birdlist=birdlist16)
blym.sfd.1yo.17 = group_with_female(blym.sfd.1yo=blym.sfd.1yo.17,b_groups=b_groups17,birdlist=birdlist17)
blym.sfd.1yo.18 = group_with_female(blym.sfd.1yo=blym.sfd.1yo.18,b_groups=b_groups18,birdlist=birdlist18)
blym.sfd.1yo.19 = group_with_female(blym.sfd.1yo=blym.sfd.1yo.19,b_groups=b_groups19,birdlist=birdlist19)



####Get non-breeding group size for one-year old males####

#Testing
blym.sfd.1yo = blym.sfd.1yo.16
birdlist = birdlist16
i=1

NB_group_size = function(blym.sfd.1yo,birdlist) {
  blym.sfd.1yo$NB_group_size = NA
  gs = data.frame(table(birdlist$Social.Group))
  colnames(gs)[1] = "Social.Group"
  for (i in 1:nrow(blym.sfd.1yo)) {
    group.size = gs %>% filter(Social.Group==blym.sfd.1yo$Social.Group[i]) %>% select(Freq)
    blym.sfd.1yo$NB_group_size[i] = group.size[1,]
  }
  return(blym.sfd.1yo)
}

#Run function
blym.sfd.1yo.16 = NB_group_size(blym.sfd.1yo=blym.sfd.1yo.16,birdlist=birdlist16)
blym.sfd.1yo.17 = NB_group_size(blym.sfd.1yo=blym.sfd.1yo.17,birdlist=birdlist17)
blym.sfd.1yo.18 = NB_group_size(blym.sfd.1yo=blym.sfd.1yo.18,birdlist=birdlist18)
blym.sfd.1yo.19 = NB_group_size(blym.sfd.1yo=blym.sfd.1yo.19,birdlist=birdlist19)





####Calculate wdeg and bdeg to unrelated females####

#Testing
blym.sfd.1yo = blym.sfd.1yo.18
i=8
agesex_net = agesex18_net
network_full = network18
bird = bird18
ind = ind18

wdeg_to_unrelated_females = function(blym.sfd.1yo,bird,agesex_net,network_full,ind) {
  blym.sfd.1yo$wdeg_unrelated_females = NA #Get blank column
  
  #Start loops
  for (i in 1:nrow(blym.sfd.1yo)) {
    if(is.na(blym.sfd.1yo$SightFreq.bfm[i])) {} #Do nothing if is NA - means he molted before social observations and can't measure connections before molt - from prev function
    else { #otherwise...
      male = blym.sfd.1yo$Bird[i]
      ind.male = ind[ind$Bird %in% male,] #Get sightings of male of interest
      ind.male.b = ind.male %>% filter(ind.male$Molt>33) #Get sightings where that male was intermediate
      
      #Get network before male of interest was bright
      if(nrow(ind.male.b)!=0) { #If the male was ever seen as bright
        bright.date = ind.male.b[order(ind.male.b$Date),]$Date[1] #Order by date and get first bright date
        bird.m = bird %>% filter(Bird %in% male) #Get social observations of bird of interest
        bird.d = bird %>% filter(Date<bright.date) #Subset bird dataframe to dates before first seen as bright
        bird.d = bird.d %>% select(Bird,Sighting) #Select columns to calculate gbi
        gbi = get_group_by_individual(bird.d, data_format = "individuals") #calculate gbi
        network = get_network(gbi, data_format= "GBI", association_index = "SRI") #calculate network
        network = network[order(rownames(network)), order(colnames(network))] #Order network after it's been calculated with subsetted gbi
        
        #Get agesex_net in order of new network
        agesex_net_no = agesex_net %>% filter(Bird %in% rownames(network))
        agesex_net_no = agesex_net_no[order(agesex_net_no$Bird),]
        
        fwn = blym.sfd.1yo$FWnumber[i] #Get Fwnumber of male
        rel.fwn = rel %>% filter(ind1==fwn) #Filter relatedness values for pairs with him
        agesex_net_no = agesex_net_no %>% select(Bird,FWNo,Sex) #Select only columns needed from agesex_net which is same birds in network
        rel.fwn = merge(rel.fwn,agesex_net_no,by.x="ind2",by.y="FWNo",all.y=T) #Merge but keep all of agesex_net so can put in order of network
        rel.fwn = rel.fwn %>% arrange(Bird) #Get in network order
        if(all(rel.fwn$Bird==rownames(network))) {} else {stop("agesex_net_no and network not equal")} #Check order
        
        #Get related or unrelated - rowsum seems to need "==" rather than a "<" term. 
        rel.fwn$Related = NA
        for (b in 1:nrow(rel.fwn)) {if(rel.fwn$wang[b]<=0.15 | is.na(rel.fwn$wang[b])) {rel.fwn$Related[b]="Unrelated"} else {rel.fwn$Related[b]="Related"}}
        
        #Get wdeg values
        h = as.numeric(rownames(agesex_net_no[agesex_net_no$Bird==blym.sfd.1yo$Bird[i],])) #Get row of bird of interest
        wdeg_URF = rowsum(network,(rel.fwn$Sex=="F" & rel.fwn$Related=="Unrelated"))[,h] #[,h] gets the row for the male of interest
        
        #Put values into dataframe
        blym.sfd.1yo$wdeg_unrelated_females[i] = wdeg_URF[2] #2 is true if it's there. There should always be false. 2 = NA if not there
    
        
      } else {  #Do same for males that did not turn bright during non-breeding social observations - don't have to recalculate network here
        
        fwn = blym.sfd.1yo$FWnumber[i] #Get Fwnumber of male
        rel.fwn = rel %>% filter(ind1==fwn) #Filter relatedness values for pairs with him
        agesex_net = agesex_net %>% select(Bird,FWNo,Sex) #Select only columns needed from agesex_net which is same birds in network
        rel.fwn = merge(rel.fwn,agesex_net,by.x="ind2",by.y="FWNo",all.y=T) #Merge but keep all of agesex_net so can put in order of network
        rel.fwn = rel.fwn %>% arrange(Bird) #Get in network order
        if(all(rel.fwn$Bird==rownames(network_full))) {} else {stop("agesex_net_no and network not equal")} #Check order
        
        #Get related or unrelated - rowsum seems to need "==" rather than a "<" term. 
        rel.fwn$Related = NA
        for (b in 1:nrow(rel.fwn)) {if(rel.fwn$wang[b]<=0.15 | is.na(rel.fwn$wang[b])) {rel.fwn$Related[b]="Unrelated"} else {rel.fwn$Related[b]="Related"}}
        
        #Get wdeg values
        h = as.numeric(rownames(agesex_net[agesex_net$Bird==blym.sfd.1yo$Bird[i],])) #Get row of bird of interest
        wdeg_URF = rowsum(network_full,(rel.fwn$Sex=="F" & rel.fwn$Related=="Unrelated"))[,h] #[,h] gets the row for the male of interest
        
        #Put values into dataframe
        blym.sfd.1yo$wdeg_unrelated_females[i] = wdeg_URF[2] #2 is true if it's there. There should always be false. 2 = NA if not there
      }
    } 
  }
  return(blym.sfd.1yo)
}

#Run function
blym.sfd.1yo.16 = wdeg_to_unrelated_females(blym.sfd.1yo = blym.sfd.1yo.16,bird=bird16,agesex_net = agesex16_net,network_full = network16,ind=ind16)
blym.sfd.1yo.17 = wdeg_to_unrelated_females(blym.sfd.1yo = blym.sfd.1yo.17,bird=bird17,agesex_net = agesex17_net,network_full = network17,ind=ind17)
blym.sfd.1yo.18 = wdeg_to_unrelated_females(blym.sfd.1yo = blym.sfd.1yo.18,bird=bird18,agesex_net = agesex18_net,network_full = network18,ind=ind18)
blym.sfd.1yo.19 = wdeg_to_unrelated_females(blym.sfd.1yo = blym.sfd.1yo.19,bird=bird19,agesex_net = agesex19_net,network_full = network19,ind=ind19)




bdeg_to_unrelated_females = function(blym.sfd.1yo,bird,agesex_net,network_full,ind) {
  blym.sfd.1yo$bdeg_unrelated_females = NA #Get blank column
  
  #Make the network_full binary
  network_full[network_full>0] <- 1 #Make the network binary
  
  #Start loops
  for (i in 1:nrow(blym.sfd.1yo)) {
    if(is.na(blym.sfd.1yo$SightFreq.bfm[i])) {} #Do nothing if is NA - means he molted before social observations and can't measure connections before molt - from prev function
    else { #otherwise...
      male = blym.sfd.1yo$Bird[i]
      ind.male = ind[ind$Bird %in% male,] #Get sightings of male of interest
      ind.male.b = ind.male %>% filter(ind.male$Molt>33) #Get sightings where that male was intermediate
      
      #Get network before male of interest was bright
      if(nrow(ind.male.b)!=0) { #If the male was ever seen as bright
        bright.date = ind.male.b[order(ind.male.b$Date),]$Date[1] #Order by date and get first bright date
        bird.m = bird %>% filter(Bird %in% male) #Get social observations of bird of interest
        bird.d = bird %>% filter(Date<bright.date) #Subset bird dataframe to dates before first seen as bright
        bird.d = bird.d %>% select(Bird,Sighting) #Select columns to calculate gbi
        gbi = get_group_by_individual(bird.d, data_format = "individuals") #calculate gbi
        network = get_network(gbi, data_format= "GBI", association_index = "SRI") #calculate network
        network = network[order(rownames(network)), order(colnames(network))] #Order network after it's been calculated with subsetted gbi
        network[network>0] <- 1 #Make the network binary
        
        #Get agesex_net in order of new network
        agesex_net_no = agesex_net %>% filter(Bird %in% rownames(network))
        agesex_net_no = agesex_net_no[order(agesex_net_no$Bird),]
        
        fwn = blym.sfd.1yo$FWnumber[i] #Get Fwnumber of male
        rel.fwn = rel %>% filter(ind1==fwn) #Filter relatedness values for pairs with him
        agesex_net_no = agesex_net_no %>% select(Bird,FWNo,Sex) #Select only columns needed from agesex_net which is same birds in network
        rel.fwn = merge(rel.fwn,agesex_net_no,by.x="ind2",by.y="FWNo",all.y=T) #Merge but keep all of agesex_net so can put in order of network
        rel.fwn = rel.fwn %>% arrange(Bird) #Get in network order
        if(all(rel.fwn$Bird==rownames(network))) {} else {stop("agesex_net_no and network not equal")} #Check order
        
        #Get related or unrelated - rowsum seems to need "==" rather than a "<" term. 
        rel.fwn$Related = NA
        for (b in 1:nrow(rel.fwn)) {if(rel.fwn$wang[b]<=0.15 | is.na(rel.fwn$wang[b])) {rel.fwn$Related[b]="Unrelated"} else {rel.fwn$Related[b]="Related"}}
        
        #Get bdeg values
        h = as.numeric(rownames(agesex_net_no[agesex_net_no$Bird==blym.sfd.1yo$Bird[i],])) #Get row of bird of interest
        bdeg_URF = rowsum(network,(rel.fwn$Sex=="F" & rel.fwn$Related=="Unrelated"))[,h] #[,h] gets the row for the male of interest
        
        #Put values into dataframe
        blym.sfd.1yo$bdeg_unrelated_females[i] = bdeg_URF[2] #2 is true if it's there. There should always be false. 2 = NA if not there
        
        
      } else {  #Do same for males that did not turn bright during non-breeding social observations - don't have to recalculate network here
        
        fwn = blym.sfd.1yo$FWnumber[i] #Get Fwnumber of male
        rel.fwn = rel %>% filter(ind1==fwn) #Filter relatedness values for pairs with him
        agesex_net = agesex_net %>% select(Bird,FWNo,Sex) #Select only columns needed from agesex_net which is same birds in network
        rel.fwn = merge(rel.fwn,agesex_net,by.x="ind2",by.y="FWNo",all.y=T) #Merge but keep all of agesex_net so can put in order of network
        rel.fwn = rel.fwn %>% arrange(Bird) #Get in network order
        if(all(rel.fwn$Bird==rownames(network_full))) {} else {stop("agesex_net_no and network not equal")} #Check order
        
        #Get related or unrelated - rowsum seems to need "==" rather than a "<" term. 
        rel.fwn$Related = NA
        for (b in 1:nrow(rel.fwn)) {if(rel.fwn$wang[b]<=0.15 | is.na(rel.fwn$wang[b])) {rel.fwn$Related[b]="Unrelated"} else {rel.fwn$Related[b]="Related"}}
        
        #Get bdeg values
        h = as.numeric(rownames(agesex_net[agesex_net$Bird==blym.sfd.1yo$Bird[i],])) #Get row of bird of interest
        bdeg_URF = rowsum(network_full,(rel.fwn$Sex=="F" & rel.fwn$Related=="Unrelated"))[,h] #[,h] gets the row for the male of interest
        
        #Put values into dataframe
        blym.sfd.1yo$bdeg_unrelated_females[i] = bdeg_URF[2] #2 is true if it's there. There should always be false. 2 = NA if not there
      }
    } 
  }
  return(blym.sfd.1yo)
}

#Run function
blym.sfd.1yo.16 = bdeg_to_unrelated_females(blym.sfd.1yo = blym.sfd.1yo.16,bird=bird16,agesex_net = agesex16_net,network_full = network16,ind=ind16)
blym.sfd.1yo.17 = bdeg_to_unrelated_females(blym.sfd.1yo = blym.sfd.1yo.17,bird=bird17,agesex_net = agesex17_net,network_full = network17,ind=ind17)
blym.sfd.1yo.18 = bdeg_to_unrelated_females(blym.sfd.1yo = blym.sfd.1yo.18,bird=bird18,agesex_net = agesex18_net,network_full = network18,ind=ind18)
blym.sfd.1yo.19 = bdeg_to_unrelated_females(blym.sfd.1yo = blym.sfd.1yo.19,bird=bird19,agesex_net = agesex19_net,network_full = network19,ind=ind19)


##Make NA's zero
blym.sfd.1yo.16$wdeg_unrelated_females[is.na(blym.sfd.1yo.16$wdeg_unrelated_females)]<-0
blym.sfd.1yo.17$wdeg_unrelated_females[is.na(blym.sfd.1yo.17$wdeg_unrelated_females)]<-0
blym.sfd.1yo.18$wdeg_unrelated_females[is.na(blym.sfd.1yo.18$wdeg_unrelated_females)]<-0
blym.sfd.1yo.19$wdeg_unrelated_females[is.na(blym.sfd.1yo.19$wdeg_unrelated_females)]<-0

blym.sfd.1yo.16$bdeg_unrelated_females[is.na(blym.sfd.1yo.16$bdeg_unrelated_females)]<-0
blym.sfd.1yo.17$bdeg_unrelated_females[is.na(blym.sfd.1yo.17$bdeg_unrelated_females)]<-0
blym.sfd.1yo.18$bdeg_unrelated_females[is.na(blym.sfd.1yo.18$bdeg_unrelated_females)]<-0
blym.sfd.1yo.19$bdeg_unrelated_females[is.na(blym.sfd.1yo.19$bdeg_unrelated_females)]<-0




####Make sure that for each 1yo male we collected social data throughout the non-breeding season for him####
ind16.blym = ind16 %>% filter(Bird %in% blym.sfd.1yo.16$Bird)
ind17.blym = ind17 %>% filter(Bird %in% blym.sfd.1yo.17$Bird)
ind18.blym = ind18 %>% filter(Bird %in% blym.sfd.1yo.18$Bird)
ind19.blym = ind19 %>% filter(Bird %in% blym.sfd.1yo.19$Bird)

library(ggplot2)
ggplot(data=ind16.blym,aes(x=Date,y=Bird)) + geom_point() + geom_line() + xlim(as.Date("2015-06-15"),as.Date("2015-09-01"))
ggplot(data=ind17.blym,aes(x=Date,y=Bird)) + geom_point() + geom_line() + xlim(as.Date("2016-06-15"),as.Date("2016-09-01"))
ggplot(data=ind18.blym,aes(x=Date,y=Bird)) + geom_point() + geom_line() + xlim(as.Date("2017-06-15"),as.Date("2017-09-01"))
ggplot(data=ind19.blym,aes(x=Date,y=Bird)) + geom_point() + geom_line() + xlim(as.Date("2018-06-15"),as.Date("2018-09-01"))

#All look good except for LZI who gets removed in whether a male molted script because he wasn't seen in the breeding season




####Combine years and write 1-year-old male social data####

#Add year first
blym.sfd.1yo.16$Year = "2016"
blym.sfd.1yo.17$Year = "2017"
blym.sfd.1yo.18$Year = "2018"
blym.sfd.1yo.19$Year = "2019"

#Combine
oneyo_social_connections = rbind(blym.sfd.1yo.16,blym.sfd.1yo.17,blym.sfd.1yo.18,blym.sfd.1yo.19)

#Remove individuals that were not in social observations before they made it to intermediate 
oneyo_social_connections = oneyo_social_connections %>% filter(SightFreq.bfm!=0)

##Write
#write.csv(oneyo_social_connections,row.names = F,here::here("Output files","oneyo_social_connections_observed.csv"))

#Can then bring this file into whether male molted or not and look at models using social connections



