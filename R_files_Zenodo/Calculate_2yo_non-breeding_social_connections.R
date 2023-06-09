###Calculate 2yo non-breeding social connections script

#Using this script to calculate observed non-breeding social connections of 2yo males. 
#Files created this script will be used to determine whether and how non-breeding social environment influences
#whether and when males molt into nuptial plumage

#Running this for only 2017 since using data from that year to look at social values and timing of molt in 2yo males

library(here)
library(dplyr)
library(reshape2)
library(lubridate)
library(ggplot2)
library(cowplot)
library(asnipe)


##Load bird17 files - non-breeding social interaction data - comes from non-breeding social structure paper
#Includes control and removal communities
bird17 = read.csv(here::here("Input files","bird17R.csv"))
bird17$Date = as.Date(bird17$Date,"%m/%d/%y") #Convert date 
bird17 = bird17 %>% arrange(Date) #Order by date
bird17$jdate = yday(bird17$Date)

##Load individual file - need to know when males turned bright if they did
ind17 = read.csv(here::here("Input files","Individuals2017 7_13_20.csv"))
ind17$Molt = as.numeric(as.character(ind17$Molt))
ind17$Date = as.Date(ind17$Date,"%m/%d/%y")
ind17$jdate = yday(ind17$Date)

#Agesex file
agesex17 = read.csv(here::here("Input files","Ages and Sexes 2017.csv"))

##Load in breeding groups for upcoming breeding seasons
b_groups17 = read.csv(here::here("Input files","Breeding Groups 2017_2017colors_short.csv"),stringsAsFactors = F)
b_groups17$Date.Created = as.Date(b_groups17$Date.Created,"%m/%d/%y")

##Load in list of 2yo males in 2017 with good molt data
males = read.csv(here::here("Output files","2yo_male_list_2017.csv"))

##Load in relatedness data 
rel = read.csv(here::here("Input files","RelatednessEstimates_8xdenovo80_14to19.csv"))





####Measure each 2yo male's observed non-breeding social connections to different social classes before he molted#### 

##Get individual set up for network
individuals17 = data.frame(bird17$Bird, bird17$Sighting)
colnames(individuals17) = c("ID","Sighting")

##Get Group by Individual Matrix 
gbi17 = get_group_by_individual(individuals17, data_format = "individuals")

##Get filtered network so I can calculate degree
network17 = get_network(gbi17, data_format= "GBI", association_index = "SRI") 
network17 = network17[order(rownames(network17)), order(colnames(network17))]
birdorder17 = data.frame(rownames(network17))
colnames(birdorder17)[1] = "Bird"
birdlist17 = birdorder17

##Add age and sex data
agesex17.s = agesex17 %>% select(FWNo,Bird,Sex,Current.Age,Ageexact.)
birdlist17 = merge(birdlist17,agesex17.s,by="Bird") %>% arrange(Bird)

##Make sure birdlist and network in same order
birdlist17$Bird == rownames(network17)

##Get relatedness dataframe set up
rel2 = rel %>% select(Pair,ind2,ind1,wang) #Get columns in different order 
colnames(rel2)[2:3] = c("ind1","ind2")
rel = rbind(rel,rel2) #Combine so all individuals in both columns

#Get relatedness values for birds present in 2017 network
rel = rel %>% filter(ind1 %in% birdlist17$FWNo) %>% filter(ind2 %in% birdlist17$FWNo)




####Two-year-old male weighted degree function
#Function to get weighted degree connections to social classes before molt

connections_before_bright_2yo_wdeg = function(males,bird17,network17,birdlist17) {
  males$wdeg_12yo_males = NA #Total weighted degree to all 1 and 2yo males 
  males$wdeg_old_males = NA #Total weighted degree to all older males(3+)
  males$wdeg_females = NA #Total weighted degree to all females
  males$wdeg_female_paired = NA #Weighted degree of connection to female the male eventually paired with
  males$wdeg_all = NA #Total weighted degree
  males$wdeg_females_unrelated = NA #Weighted degree to unrelated females
  
  #Start loops
  for (i in 1:nrow(males)) {
    male = males$Bird[i]

    #Get network before male of interest was bright
    bird.d = bird17 %>% filter(jdate<males$molt.date[i]) #Subset bird dataframe to dates before first seen as bright
    bird.d = bird.d %>% select(Bird,Sighting) #Select columns to calculate gbi
    gbi = get_group_by_individual(bird.d, data_format = "individuals") #calculate gbi
    network = get_network(gbi, data_format= "GBI", association_index = "SRI") #calculate network
    network = network[order(rownames(network)), order(colnames(network))] #Order network after it's been calculated with subsetted gbi
    
    #Get agesex_net that matches network order
    birdlist17.no = birdlist17 %>% filter(Bird %in% rownames(network)) %>% arrange(Bird)

    #Check that birdlist17.no matches network order
    if(all(birdlist17.no$Bird==rownames(network))) {} else {stop("birdlist17.no and network not equal")} 
    
    #Get row number of male of interest in the network/birdlist17.no so I can pull out the male's value in the rowsums step below
    h = as.numeric(rownames(birdlist17.no[birdlist17.no$Bird==as.character(male),]))

    #Calculate connections
    #rowsum calculates the sum of the false [1] and true [2] conditions of the columns that match the characteristics
    #Need birdlist17.no and network to be in same order
    wdeg_12yo_males = rowsum(network,(birdlist17.no$Sex=="M" & birdlist17.no$Current.Age<=2))[,h] #[,h] gets the row for the male of interest
    wdeg_old_males = rowsum(network,(birdlist17.no$Sex=="M" & birdlist17.no$Current.Age>=3))[,h]
    wdeg_females = rowsum(network,(birdlist17.no$Sex=="F"))[,h]
    wdeg_female_paired = rowsum(network,(birdlist17.no$Bird==as.character(males$Female[i])))[,h] 
    wdeg_all = rowSums(network)[h]
    
    #Get relatedness between male and females
    male.fwn = birdlist17.no %>% filter(Bird==as.character(male)) %>% select(FWNo) #Get male's fwnumber
    rel.male = rel %>% filter(ind1 == male.fwn[1,]) #Get related dataframe for only pairs with male
    birdlist17.no = merge(birdlist17.no,rel.male,by.x="FWNo",by.y="ind2",all.x=T) #Combine with birdlist
    birdlist17.no[is.na(birdlist17.no$wang),]$wang <-0 #Make NA's zero
    wdeg_females_unrelated = rowsum(network,(birdlist17.no$Sex=="F" & birdlist17.no$wang<=0.15))[,h]
    
    #Put values into dataframe
    males$wdeg_12yo_males[i] = wdeg_12yo_males[2] #2 is true if it's there. There should always be false. 2 = NA if not there
    males$wdeg_old_males[i] = wdeg_old_males[2]
    males$wdeg_females[i] = wdeg_females[2]
    males$wdeg_female_paired[i] = wdeg_female_paired[2]
    males$wdeg_all[i] = wdeg_all
    males$wdeg_females_unrelated[i] = wdeg_females_unrelated[2]
  }
  return(males)
}

#Calculate two year old male connections using weighted degree
males.wdeg=connections_before_bright_2yo_wdeg(males=males,bird17=bird17,birdlist17=birdlist17,network17=network17)


###Write to .csv
#write.csv(males.wdeg,here::here("Output files","twoyo_social_connections2.csv"),row.names = F)





