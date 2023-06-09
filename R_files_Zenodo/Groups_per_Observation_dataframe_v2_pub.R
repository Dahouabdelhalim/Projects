#Groups per observation dataframe and graphs

#Use this script to create the groups per observation dataframe used in the climwin sliding window analysis


####Load Data####

library(here)

#Load in bird dataframes
bird16 = read.csv(here::here("Input files","bird16.csv"))
bird17 = read.csv(here::here("Input files","bird17R.csv"))
bird18 = read.csv(here::here("Input files","bird18.csv"))
bird19 = read.csv(here::here("Input files","bird19.csv"))

#Take out removal community sightings from 2017 
RorC = read.csv(here::here("Input files","RemovalorControl 2017.csv"), header=TRUE)
bird17R = merge(bird17,RorC,by="Bird")
bird17R = bird17R[which(bird17R$TreatmentGroup=="Removal"),]
bird17 = bird17[!bird17$Sighting %in% bird17R$Sighting,]

#Load in social group IDs
group16 = read.csv(here::here("Output files","group16.csv"))
group17 = read.csv(here::here("Output files","group17.csv"))
group18 = read.csv(here::here("Output files","group18.csv"))
group19 = read.csv(here::here("Output files","group19.csv"))

#Get social group sizes
group16t = data.frame(table(group16$Social.Group))
colnames(group16t) = c("Social.Group","Social.Group.Size")
group17t = data.frame(table(group17$Social.Group))
colnames(group17t) = c("Social.Group","Social.Group.Size")
group18t = data.frame(table(group18$Social.Group))
colnames(group18t) = c("Social.Group","Social.Group.Size")
group19t = data.frame(table(group19$Social.Group))
colnames(group19t) = c("Social.Group","Social.Group.Size")

#Remove birds from group list if they're in a group on their own. These are birds that were above the 
#cutoff in the dendrogram, meaning we didn't measure their social environment well and they shouldn't 
#be assigned a group
group16tsingle = group16t[which(group16t$Social.Group.Size==1),]
group17tsingle = group17t[which(group17t$Social.Group.Size==1),]
group18tsingle = group18t[which(group18t$Social.Group.Size==1),]
group19tsingle = group19t[which(group19t$Social.Group.Size==1),]

group16 = group16[!group16$Social.Group %in% group16tsingle$Social.Group,]
group17 = group17[!group17$Social.Group %in% group17tsingle$Social.Group,]
group18 = group18[!group18$Social.Group %in% group18tsingle$Social.Group,]
group19 = group19[!group19$Social.Group %in% group19tsingle$Social.Group,]


#Load in birds dataframes - shows number of banded birds seen and number of unknowns (banded or 
#unbanded birds that we didn't get the band combination of). 
birds16 = read.csv(here::here("Input files","birds16_banded_unks.csv"))
birds17 = read.csv(here::here("Input files","birds17_banded_unks.csv"))
birds18 = read.csv(here::here("Input files","birds18_banded_unks.csv"))
birds19 = read.csv(here::here("Input files","birds19_banded_unks.csv"))

#Combine bird and group dataframes
bird16 = merge(bird16,group16,by="Bird",all.x=T)
bird17 = merge(bird17,group17,by="Bird",all.x=T)
bird18 = merge(bird18,group18,by="Bird",all.x=T)
bird19 = merge(bird19,group19,by="Bird",all.x=T)

#Get banded birds with no social group
bird16NA = bird16[is.na(bird16$Social.Group),]
bird17NA = bird17[is.na(bird17$Social.Group),]
bird18NA = bird18[is.na(bird18$Social.Group),]
bird19NA = bird19[is.na(bird19$Social.Group),]

#Remove NAs so not counting them as separate groups yet
bird16s = bird16[!is.na(bird16$Social.Group),]
bird17s = bird17[!is.na(bird17$Social.Group),]
bird18s = bird18[!is.na(bird18$Social.Group),]
bird19s = bird19[!is.na(bird19$Social.Group),]

#Get number of NAs per sighting (Birds with no social group per sighting)
NAbysight16 = data.frame(table(bird16NA$Sighting))
colnames(NAbysight16)=c("Sighting","NAs")
NAbysight17 = data.frame(table(bird17NA$Sighting))
colnames(NAbysight17)=c("Sighting","NAs")
NAbysight18 = data.frame(table(bird18NA$Sighting))
colnames(NAbysight18)=c("Sighting","NAs")
NAbysight19 = data.frame(table(bird19NA$Sighting))
colnames(NAbysight19)=c("Sighting","NAs")


#Get number of groups per sighting for birds that have groups
groupbysight16 = aggregate(bird16s["Social.Group"],by=bird16s["Sighting"],function(x) length(unique(x)))
colnames(groupbysight16)[2] = "Social.Groups"
groupbysight17 = aggregate(bird17s["Social.Group"],by=bird17s["Sighting"],function(x) length(unique(x)))
colnames(groupbysight17)[2] = "Social.Groups"
groupbysight18 = aggregate(bird18s["Social.Group"],by=bird18s["Sighting"],function(x) length(unique(x)))
colnames(groupbysight18)[2] = "Social.Groups"
groupbysight19 = aggregate(bird19s["Social.Group"],by=bird19s["Sighting"],function(x) length(unique(x)))
colnames(groupbysight19)[2] = "Social.Groups"


#Merge groups by sighting with NAs by sighting
groupbysight16 = merge(groupbysight16,NAbysight16,by="Sighting",all.x=T)
groupbysight17 = merge(groupbysight17,NAbysight17,by="Sighting",all.x=T)
groupbysight18 = merge(groupbysight18,NAbysight18,by="Sighting",all.x=T)
groupbysight19 = merge(groupbysight19,NAbysight19,by="Sighting",all.x=T)

#Set NAs in NA column to zero
groupbysight16[3][is.na(groupbysight16$NAs),] <- 0 
groupbysight17[3][is.na(groupbysight17$NAs),] <- 0 
groupbysight18[3][is.na(groupbysight18$NAs),] <- 0 
groupbysight19[3][is.na(groupbysight19$NAs),] <- 0 


#Bring in number of birds in each sighting, observation ID, and date
numberbirds16 = data.frame(bird16$Sighting,bird16$Number,bird16$Observation,bird16$Date)
numberbirds16 = numberbirds16[!duplicated(numberbirds16$bird16.Sighting),]
colnames(numberbirds16) = c("Sighting","Number","Observation","Date")
numberbirds17 = data.frame(bird17$Sighting,bird17$Number,bird17$Observation,bird17$Date)
numberbirds17 = numberbirds17[!duplicated(numberbirds17$bird17.Sighting),]
colnames(numberbirds17) = c("Sighting","Number","Observation","Date")
numberbirds18 = data.frame(bird18$Sighting,bird18$Number,bird18$Observation,bird18$Date)
numberbirds18 = numberbirds18[!duplicated(numberbirds18$bird18.Sighting),]
colnames(numberbirds18) = c("Sighting","Number","Observation","Date")
numberbirds19 = data.frame(bird19$Sighting,bird19$Number,bird19$Observation,bird19$Date)
numberbirds19 = numberbirds19[!duplicated(numberbirds19$bird19.Sighting),]
colnames(numberbirds19) = c("Sighting","Number","Observation","Date")

groupbysight16 = merge(groupbysight16,numberbirds16,by="Sighting")
groupbysight17 = merge(groupbysight17,numberbirds17,by="Sighting")
groupbysight18 = merge(groupbysight18,numberbirds18,by="Sighting")
groupbysight19 = merge(groupbysight19,numberbirds19,by="Sighting")

#Format date as date
groupbysight16$Date = as.Date(groupbysight16$Date,"%m/%d/%y")
groupbysight17$Date = as.Date(groupbysight17$Date,"%m/%d/%y")
groupbysight18$Date = as.Date(groupbysight18$Date,"%m/%d/%y")
groupbysight19$Date = as.Date(groupbysight19$Date,"%m/%d/%y")


#Bring in number of unknown birds per sighting (includes Unknowns (didn't get band combo) and UNBs)
unknowns16 = data.frame(birds16$Sighting,birds16$Banded,birds16$Unknowns)
colnames(unknowns16) = c("Sighting","Banded","Unknowns")
unknowns17 = data.frame(birds17$Sighting,birds17$Banded,birds17$Unknowns)
colnames(unknowns17) = c("Sighting","Banded","Unknowns")
unknowns18 = data.frame(birds18$Sighting,birds18$Banded,birds18$Unknowns)
colnames(unknowns18) = c("Sighting","Banded","Unknowns")
unknowns19 = data.frame(birds19$Sighting,birds19$Banded,birds19$Unknowns)
colnames(unknowns19) = c("Sighting","Banded","Unknowns")

groupbysight16 = merge(groupbysight16,unknowns16,by="Sighting")
groupbysight17 = merge(groupbysight17,unknowns17,by="Sighting")
groupbysight18 = merge(groupbysight18,unknowns18,by="Sighting")
groupbysight19 = merge(groupbysight19,unknowns19,by="Sighting")

#Calculate Percent Unknowns
groupbysight16$PercentUnk = groupbysight16$Unknowns/groupbysight16$Number
groupbysight17$PercentUnk = groupbysight17$Unknowns/groupbysight17$Number
groupbysight18$PercentUnk = groupbysight18$Unknowns/groupbysight18$Number
groupbysight19$PercentUnk = groupbysight19$Unknowns/groupbysight19$Number


#Get number of grouped birds - birds whose band combos were recorded and are grouped
groupbysight16$Grouped.Birds = groupbysight16$Number-groupbysight16$NAs-groupbysight16$Unknowns
groupbysight16$PercentGrouped = groupbysight16$Grouped.Birds/groupbysight16$Number
groupbysight17$Grouped.Birds = groupbysight17$Number-groupbysight17$NAs-groupbysight17$Unknowns
groupbysight17$PercentGrouped = groupbysight17$Grouped.Birds/groupbysight17$Number
groupbysight18$Grouped.Birds = groupbysight18$Number-groupbysight18$NAs-groupbysight18$Unknowns
groupbysight18$PercentGrouped = groupbysight18$Grouped.Birds/groupbysight18$Number
groupbysight19$Grouped.Birds = groupbysight19$Number-groupbysight19$NAs-groupbysight19$Unknowns
groupbysight19$PercentGrouped = groupbysight19$Grouped.Birds/groupbysight19$Number

#Remove sightings where more than 50% of birds were Unknown combos or Unbanded
groupbysight16 = groupbysight16[which(groupbysight16$PercentUnk<=0.5),]
groupbysight17 = groupbysight17[which(groupbysight17$PercentUnk<=0.5),]
groupbysight18 = groupbysight18[which(groupbysight18$PercentUnk<=0.5),]
groupbysight19 = groupbysight19[which(groupbysight19$PercentUnk<=0.5),]

#Remove sightings where more than 50% of birds were not banded and did not have groups
groupbysight16sub = groupbysight16[which(groupbysight16$PercentGrouped>=0.5),]
groupbysight17sub = groupbysight17[which(groupbysight17$PercentGrouped>=0.5),]
groupbysight18sub = groupbysight18[which(groupbysight18$PercentGrouped>=0.5),]
groupbysight19sub = groupbysight19[which(groupbysight19$PercentGrouped>=0.5),]


###Add an additional group to number of social groups if many unknowns present or if birds without 
#a group are present
#First calculate group size totals of groups that are in each sighting if all groups were present.
#If the total number of birds is higher than the number of grouped birds, then any unknowns are 
#probably just grouped birds that didn't get seen
#If the group total is less than the number of birds, then unknowns can be added to the pool for calculating 
#other groups
##Merge bird dataframes with group sizes
bird16g = merge(bird16,group16t,by="Social.Group",all.x=T)
bird17g = merge(bird17,group17t,by="Social.Group",all.x=T)
bird18g = merge(bird18,group18t,by="Social.Group",all.x=T)
bird19g = merge(bird19,group19t,by="Social.Group",all.x=T)

#Make NAs = 0 so they can be added in aggregation step below
bird16g[11][is.na(bird16g$Social.Group.Size),] <- 0 
bird17g[11][is.na(bird17g$Social.Group.Size),] <- 0 
bird18g[11][is.na(bird18g$Social.Group.Size),] <- 0 
bird19g[11][is.na(bird19g$Social.Group.Size),] <- 0 

#Get one line of each group per sighting
library(dplyr)
bird16g = bird16g %>% distinct(Sighting,Social.Group,.keep_all=T)
bird17g = bird17g %>% distinct(Sighting,Social.Group,.keep_all=T)
bird18g = bird18g %>% distinct(Sighting,Social.Group,.keep_all=T)
bird19g = bird19g %>% distinct(Sighting,Social.Group,.keep_all=T)

#Aggregate number of individuals in each group by sighting
#Get number of groups per sighting for birds that have groups
grouppossible16 = aggregate(bird16g["Social.Group.Size"],by=bird16g["Sighting"],FUN=sum)
colnames(grouppossible16)[2] = "Possible.Group.Size"
grouppossible17 = aggregate(bird17g["Social.Group.Size"],by=bird17g["Sighting"],FUN=sum)
colnames(grouppossible17)[2] = "Possible.Group.Size"
grouppossible18 = aggregate(bird18g["Social.Group.Size"],by=bird18g["Sighting"],FUN=sum)
colnames(grouppossible18)[2] = "Possible.Group.Size"
grouppossible19 = aggregate(bird19g["Social.Group.Size"],by=bird19g["Sighting"],FUN=sum)
colnames(grouppossible19)[2] = "Possible.Group.Size"

#Add into groupbysight dataframe
groupbysight16sub = merge(groupbysight16sub,grouppossible16,by="Sighting")
groupbysight17sub = merge(groupbysight17sub,grouppossible17,by="Sighting")
groupbysight18sub = merge(groupbysight18sub,grouppossible18,by="Sighting")
groupbysight19sub = merge(groupbysight19sub,grouppossible19,by="Sighting")


###Calculate whether an extra group is needed
#NAs are truly other groups (band combo seen but not in an observed group), so for sightings with an NA, 
#add an extra social group. Also if number of unknowns/UNBs outnumber total possible group size - 
#the number of birds in all social groups present combined, know there should be another social group
#That should give me a new social groups present that is more reflective of what really happened. 

#2016
groupbysight16sub$Social.Groups.Total = NA
for (i in 1:nrow(groupbysight16sub)) {
  if(groupbysight16sub$NAs[i] > 1) {
    groupbysight16sub$Social.Groups.Total[i] = groupbysight16sub$Social.Groups[i] + 1} else
      if((groupbysight16sub$Unknowns[i] + groupbysight16sub$Grouped.Birds[i]) > groupbysight16sub$Possible.Group.Size[i]) {
        groupbysight16sub$Social.Groups.Total[i] = groupbysight16sub$Social.Groups[i] + 1} else {
          groupbysight16sub$Social.Groups.Total[i] = groupbysight16sub$Social.Groups[i]
        }}
hist(groupbysight16sub$Social.Groups.Total - groupbysight16sub$Social.Groups)

#2017
groupbysight17sub$Social.Groups.Total = NA
for (i in 1:nrow(groupbysight17sub)) {
  if(groupbysight17sub$NAs[i] > 1) {
    groupbysight17sub$Social.Groups.Total[i] = groupbysight17sub$Social.Groups[i] + 1} else
      if((groupbysight17sub$Unknowns[i] + groupbysight17sub$Grouped.Birds[i]) > groupbysight17sub$Possible.Group.Size[i]) {
        groupbysight17sub$Social.Groups.Total[i] = groupbysight17sub$Social.Groups[i] + 1} else {
          groupbysight17sub$Social.Groups.Total[i] = groupbysight17sub$Social.Groups[i]
        }}
hist(groupbysight17sub$Social.Groups.Total - groupbysight17sub$Social.Groups)

#2018
groupbysight18sub$Social.Groups.Total = NA
for (i in 1:nrow(groupbysight18sub)) {
  if(groupbysight18sub$NAs[i] > 1) {
    groupbysight18sub$Social.Groups.Total[i] = groupbysight18sub$Social.Groups[i] + 1} else
      if((groupbysight18sub$Unknowns[i] + groupbysight18sub$Grouped.Birds[i]) > groupbysight18sub$Possible.Group.Size[i]) {
        groupbysight18sub$Social.Groups.Total[i] = groupbysight18sub$Social.Groups[i] + 1} else {
          groupbysight18sub$Social.Groups.Total[i] = groupbysight18sub$Social.Groups[i]
        }}
hist(groupbysight18sub$Social.Groups.Total - groupbysight18sub$Social.Groups)

#2019
groupbysight19sub$Social.Groups.Total = NA
for (i in 1:nrow(groupbysight19sub)) {
  if(groupbysight19sub$NAs[i] > 1) {
    groupbysight19sub$Social.Groups.Total[i] = groupbysight19sub$Social.Groups[i] + 1} else
      if((groupbysight19sub$Unknowns[i] + groupbysight19sub$Grouped.Birds[i]) > groupbysight19sub$Possible.Group.Size[i]) {
        groupbysight19sub$Social.Groups.Total[i] = groupbysight19sub$Social.Groups[i] + 1} else {
          groupbysight19sub$Social.Groups.Total[i] = groupbysight19sub$Social.Groups[i]
        }}
hist(groupbysight19sub$Social.Groups.Total - groupbysight19sub$Social.Groups)


#Control for observation - observations were anywhere from 5min to an hour long, control for observation
#by only using one sighting per observation. Use the sighting that had the max number of social groups present
gbs16subMax = groupbysight16sub[order(groupbysight16sub$Social.Groups.Total,decreasing=T),] 
gbs16subMax = gbs16subMax[!duplicated(gbs16subMax$Observation),] 
gbs17subMax = groupbysight17sub[order(groupbysight17sub$Social.Groups.Total,decreasing=T),] 
gbs17subMax = gbs17subMax[!duplicated(gbs17subMax$Observation),] 
gbs18subMax = groupbysight18sub[order(groupbysight18sub$Social.Groups.Total,decreasing=T),] 
gbs18subMax = gbs18subMax[!duplicated(gbs18subMax$Observation),] 
gbs19subMax = groupbysight19sub[order(groupbysight19sub$Social.Groups.Total,decreasing=T),] 
gbs19subMax = gbs19subMax[!duplicated(gbs19subMax$Observation),] 

hist(gbs16subMax$NAs + gbs16subMax$Unknowns)
hist(gbs17subMax$NAs + gbs17subMax$Unknowns)
hist(gbs18subMax$NAs + gbs18subMax$Unknowns)
hist(gbs19subMax$NAs + gbs19subMax$Unknowns)

hist(gbs16subMax$Social.Groups.Total - gbs16subMax$Social.Groups)
hist(gbs17subMax$Social.Groups.Total - gbs17subMax$Social.Groups)
hist(gbs18subMax$Social.Groups.Total - gbs18subMax$Social.Groups)
hist(gbs19subMax$Social.Groups.Total - gbs19subMax$Social.Groups)

#Plot social groups per observation by date
library(cowplot)
library(ggpubr)
A = ggplot() + geom_point(data=gbs16subMax, aes(x=Date,y=Social.Groups.Total)) + 
  geom_smooth(data=gbs16subMax, aes(x=Date,y=Social.Groups.Total),color="#D55E00") + 
  ylab("Max number of groups/observation") + ggtitle("2016") +
  scale_x_date(date_breaks = "1 month",limits=as.Date(c('2015-06-15','2015-09-01')),date_labels="%b") +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10)) + 
  theme_cowplot(font_size=12) + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5)) 
B = ggplot() + geom_point(data=gbs17subMax, aes(x=Date,y=Social.Groups.Total)) + 
  geom_smooth(data=gbs17subMax, aes(x=Date,y=Social.Groups.Total),color="#0072B2") + 
  ylab("Max number of groups/observation") + ggtitle("2017") +
  scale_x_date(date_breaks = "1 month",limits=as.Date(c('2016-06-15','2016-09-01')),date_labels="%b") +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10)) + 
  theme_cowplot(font_size=12) + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5)) 
C = ggplot() + geom_point(data=gbs18subMax, aes(x=Date,y=Social.Groups.Total)) + 
  geom_smooth(data=gbs18subMax, aes(x=Date,y=Social.Groups.Total),color="#009E73") +
  ylab("Max number of groups/observation") + ggtitle("2018") +
  scale_x_date(date_breaks = "1 month",limits=as.Date(c('2017-06-15','2017-09-01')),date_labels="%b") +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10)) + 
  theme_cowplot(font_size=12) + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5)) 
D = ggplot() + geom_point(data=gbs19subMax, aes(x=Date,y=Social.Groups.Total)) + 
  geom_smooth(data=gbs19subMax, aes(x=Date,y=Social.Groups.Total),color="#F0E442") + 
  ylab("Max number of groups/observation") + ggtitle("2019") +
  scale_x_date(date_breaks = "1 month",limits=as.Date(c('2018-06-15','2018-09-01')),date_labels="%b") +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10)) + 
  theme_cowplot(font_size=12) + theme(panel.grid.major.x = element_blank(), # remove the vertical grid lines 
                          panel.grid.major.y = element_line( size=.2, color="gray" ), legend.position = "none", # explicitly set the horizontal lines
                          plot.title = element_text(size=14, family="Helvetica",hjust=0.5)) 
ggarrange(A,B,C,D)



library(cowplot)
library(dplyr)



####Create final dataframe for running models####

##Build model dataframe

##Determine the focal group for each observation
#First remove duplicate sighting/group combinations - only want each group represented once per sighting
#This ends up including the number of sightings the focal group was seen in each observation which 
#I will use to test for an effect of observation length on groups per obs. That test will ask whether
#following the focal group for a long time could have caused us to push the focal group into other groups. 
library(reshape2)
bird16d = bird16 %>% distinct(Sighting,Social.Group,.keep_all=T)
bird16dc = dcast(bird16d,Observation~Social.Group,fun.aggregate = length) #Cast to get #times each social group appears in an obs
bird16dcm = melt(bird16dc,id.vars=1) #Melt that cast to put into a list
colnames(bird16dcm) = c("Observation","Focal.Group","Freq")
bird16dcm = bird16dcm[which(bird16dcm$Freq!=0),] #Remove zeroes
bird16dcm = bird16dcm[which(bird16dcm$Focal.Group!="NA"),] #Remove NA social groups
bird16dcm = bird16dcm[order(bird16dcm$Freq, decreasing=T),] #Order so highest numbers at top
bird16dcm = bird16dcm[!duplicated(bird16dcm$Observation),] #Take most represented social group in each observation
colnames(bird16dcm)[3] = "Focal.Sightings" #Number of sightings the focal group was in for each obs

bird17d = bird17 %>% distinct(Sighting,Social.Group,.keep_all=T)
bird17dc = dcast(bird17d,Observation~Social.Group,fun.aggregate = length) #Cast to get #times each social group appears in an obs
bird17dcm = melt(bird17dc,id.vars=1) #Melt that cast to put into a list
colnames(bird17dcm) = c("Observation","Focal.Group","Freq")
bird17dcm = bird17dcm[which(bird17dcm$Freq!=0),] #Remove zeroes
bird17dcm = bird17dcm[which(bird17dcm$Focal.Group!="NA"),] #Remove NA social groups
bird17dcm = bird17dcm[order(bird17dcm$Freq, decreasing=T),] #Order so highest numbers at top
bird17dcm = bird17dcm[!duplicated(bird17dcm$Observation),] #Take most represented social group in each observation
colnames(bird17dcm)[3] = "Focal.Sightings" #Number of sightings the focal group was in for each obs

bird18d = bird18 %>% distinct(Sighting,Social.Group,.keep_all=T)
bird18dc = dcast(bird18d,Observation~Social.Group,fun.aggregate = length) #Cast to get #times each social group appears in an obs
bird18dcm = melt(bird18dc,id.vars=1) #Melt that cast to put into a list
colnames(bird18dcm) = c("Observation","Focal.Group","Freq")
bird18dcm = bird18dcm[which(bird18dcm$Freq!=0),] #Remove zeroes
bird18dcm = bird18dcm[which(bird18dcm$Focal.Group!="NA"),] #Remove NA social groups
bird18dcm = bird18dcm[order(bird18dcm$Freq, decreasing=T),] #Order so highest numbers at top
bird18dcm = bird18dcm[!duplicated(bird18dcm$Observation),] #Take most represented social group in each observation
colnames(bird18dcm)[3] = "Focal.Sightings" #Number of sightings the focal group was in for each obs

bird19d = bird19 %>% distinct(Sighting,Social.Group,.keep_all=T)
bird19dc = dcast(bird19d,Observation~Social.Group,fun.aggregate = length) #Cast to get #times each social group appears in an obs
bird19dcm = melt(bird19dc,id.vars=1) #Melt that cast to put into a list
colnames(bird19dcm) = c("Observation","Focal.Group","Freq")
bird19dcm = bird19dcm[which(bird19dcm$Freq!=0),] #Remove zeroes
bird19dcm = bird19dcm[which(bird19dcm$Focal.Group!="NA"),] #Remove NA social groups
bird19dcm = bird19dcm[order(bird19dcm$Freq, decreasing=T),] #Order so highest numbers at top
bird19dcm = bird19dcm[!duplicated(bird19dcm$Observation),] #Take most represented social group in each observation
colnames(bird19dcm)[3] = "Focal.Sightings" #Number of sightings the focal group was in for each obs

#Match to gbs dataframes
gbs16subMax = merge(gbs16subMax,bird16dcm,by="Observation")
gbs17subMax = merge(gbs17subMax,bird17dcm,by="Observation")
gbs18subMax = merge(gbs18subMax,bird18dcm,by="Observation")
gbs19subMax = merge(gbs19subMax,bird19dcm,by="Observation")


#Get focal group size
#Match to gbs dataframes
gbs16subMax = merge(gbs16subMax,group16t,by.x="Focal.Group",by.y="Social.Group")
gbs17subMax = merge(gbs17subMax,group17t,by.x="Focal.Group",by.y="Social.Group")
gbs18subMax = merge(gbs18subMax,group18t,by.x="Focal.Group",by.y="Social.Group")
gbs19subMax = merge(gbs19subMax,group19t,by.x="Focal.Group",by.y="Social.Group")


#Get community info
comm16 = read.csv(here::here("Output files","communities16.csv"))
comm17 = read.csv(here::here("Output files","communities17.csv"))
comm18 = read.csv(here::here("Output files","communities18.csv"))
comm19 = read.csv(here::here("Output files","communities19.csv"))

#Match up community to focal group
gbs16subMax = merge(gbs16subMax,comm16,by.x="Focal.Group",by.y="Group",all.x=T)
gbs17subMax = merge(gbs17subMax,comm17,by.x="Focal.Group",by.y="Group",all.x=T)
gbs18subMax = merge(gbs18subMax,comm18,by.x="Focal.Group",by.y="Group",all.x=T)
gbs19subMax = merge(gbs19subMax,comm19,by.x="Focal.Group",by.y="Group",all.x=T)

#Get community size
comm16t = data.frame(table(comm16$Community))
colnames(comm16t) = c("Community","Groups.in.Community")
comm17t = data.frame(table(comm17$Community))
colnames(comm17t) = c("Community","Groups.in.Community")
comm18t = data.frame(table(comm18$Community))
colnames(comm18t) = c("Community","Groups.in.Community")
comm19t = data.frame(table(comm19$Community))
colnames(comm19t) = c("Community","Groups.in.Community")

#Add community sizes to dataframe
gbs16subMax = merge(gbs16subMax,comm16t,by="Community")
gbs17subMax = merge(gbs17subMax,comm17t,by="Community")
gbs18subMax = merge(gbs18subMax,comm18t,by="Community")
gbs19subMax = merge(gbs19subMax,comm19t,by="Community")

###Add in waypoint IDs to use in NDVI and groups per sighting analysis 
#Can use (all.x = T) in the merge commands below to check that all waypoint numbers in sightings match up
#to an actual waypoint. In 2016 a couple JFW waypoints missing - not in waypoint file either. Must've 
#forgotten to take them. 
#In 2017 missing 20 SAD GPS points from her last two days in the season - must not have downloaded them

#2016 didn't include initials in waypoint in sightings/bird dataframes
library(tidyr)
bird16.gps = bird16
bird16.gps=unite(bird16.gps, WP, c(Initials, WP), remove=FALSE, sep="") #Match initials to waypoint
bird16.gps = bird16.gps[c(2,4)] #Select sighting and WP columns
bird16.gps = bird16.gps[!duplicated(bird16.gps$Sighting),] #Remove duplicate sightings
gbs16subMax = merge(gbs16subMax,bird16.gps,by="Sighting") #add into gbssubMax dataframe

#2017
bird17.gps = bird17
bird17.gps = bird17.gps[c(2,8)] #Select sighting and WP columns
bird17.gps = bird17.gps[!duplicated(bird17.gps$Sighting),] #Remove duplicate sightings
gbs17subMax = merge(gbs17subMax,bird17.gps,by="Sighting") #add into gbssubMax dataframe

#2018
bird18.gps = bird18
bird18.gps = bird18.gps[c(2,8)] #Select sighting and WP columns
bird18.gps = bird18.gps[!duplicated(bird18.gps$Sighting),] #Remove duplicate sightings
gbs18subMax = merge(gbs18subMax,bird18.gps,by="Sighting") #add into gbssubMax dataframe

#2019
bird19.gps = bird19
bird19.gps = bird19.gps[c(2,8)] #Select sighting and WP columns
bird19.gps = bird19.gps[!duplicated(bird19.gps$Sighting),] #Remove duplicate sightings
gbs19subMax = merge(gbs19subMax,bird19.gps,by="Sighting") #add into gbssubMax dataframe


#Add a year column
gbs16subMax$Year = "2016"
gbs17subMax$Year = "2017"
gbs18subMax$Year = "2018"
gbs19subMax$Year = "2019"


#Combine gbssubmax dataframes into single dataframe for glmms and sliding window analyses
gbssubMaxall = rbind(gbs16subMax,gbs17subMax,gbs18subMax,gbs19subMax)

#Make focal group a unique factor - currently different years use the same focal group numbers
gbssubMaxall$Focal.Group.Year = with(gbssubMaxall,interaction(gbssubMaxall$Focal.Group,gbssubMaxall$Year))

#Make community a unique factor - currently different years use the same community numbers
gbssubMaxall$Community.Year = with(gbssubMaxall,interaction(gbssubMaxall$Community,gbssubMaxall$Year))




#Get day of year
library(lubridate)
gbssubMaxall$jdate = yday(gbssubMaxall$Date)

#Select columns of interest for model
gbssubMaxall.print = gbssubMaxall[c(1,4,7,8,15,16,17,18,19,20,21,22,23)]


#Write to .csv
#write.csv(gbssubMaxall.print,here::here("Output files","gbssubmaxall.csv"),row.names = F)



