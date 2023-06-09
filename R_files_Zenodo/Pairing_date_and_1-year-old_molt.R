####1yo male pairing date script

#Use this script to see how pairing date influences how far 1-year-old males moult into nuptial plumage


####Load data####

library(here)
library(dplyr)
library(reshape2)
library(lubridate)
library(ggplot2)
library(cowplot)
library(lme4)
library(ggpubr)

##Load bird files - non-breeding social network data - come from non-breeding social structure paper
bird16 = read.csv(here::here("Input files","bird16.csv"))
bird16$Date = as.Date(bird16$Date,"%m/%d/%y") #Convert date 
bird16 = bird16 %>% arrange(Date) #Order by date
bird17 = read.csv(here::here("Input files","bird17R.csv"))
bird17$Date = as.Date(bird17$Date,"%m/%d/%y") #Convert date 
bird17 = bird17 %>% arrange(Date) #Order by date
bird18 = read.csv(here::here("Input files","bird18.csv"))
bird18$Date = as.Date(bird18$Date,"%m/%d/%y") #Convert date 
bird18 = bird18 %>% arrange(Date) #Order by date
bird19 = read.csv(here::here("Input files","bird19.csv"))
bird19$Date = as.Date(bird19$Date,"%m/%d/%y") #Convert date 
bird19 = bird19 %>% arrange(Date) #Order by date

##Load bird files from breeding - social network data collected during beginning of breeding season
bird18r = read.csv(here::here("Input files","bird18_breeding.csv"))
bird18r$Date = as.Date(bird18r$Date)#Convert date 
bird18r = bird18r %>% arrange(Date) #Order by date
bird19r = read.csv(here::here("Input files","bird19_breeding.csv"))
bird19r$Date = as.Date(bird19r$Date) #Convert date 
bird19r = bird19r %>% arrange(Date) #Order by date

##Combine bird files
bird18 = bird18 %>% select(-Kerfuffle)
bird18 = rbind(bird18,bird18r)
bird18 = bird18 %>% arrange(Date)
bird19 = bird19 %>% select(-Kerfuffle)
bird19 = rbind(bird19,bird19r)
bird19 = bird19 %>% arrange(Date)

##Load individual files - need to know when males turned bright if they did
ind16 = read.csv(here::here("Input files","Individuals2016 7_13_20.csv"))
ind16$Molt = as.numeric(as.character(ind16$Molt))
ind16$Date = as.Date(ind16$Date,"%m/%d/%y")
ind17 = read.csv(here::here("Input files","Individuals2017 7_13_20.csv"))
ind17$Molt = as.numeric(as.character(ind17$Molt))
ind17$Date = as.Date(ind17$Date,"%m/%d/%y")
ind18 = read.csv(here::here("Input files","Individuals2018 7_13_20.csv"))
ind18$Molt = as.numeric(as.character(ind18$Molt))
ind18$Date = as.Date(ind18$Date,"%m/%d/%y")
ind19 = read.csv(here::here("Input files","Individuals2019 7_13_20.csv"))
ind19$Molt = as.numeric(as.character(ind19$Molt))
ind19$Date = as.Date(ind19$Date,"%m/%d/%y")

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

##Load in birdlist files to get Fwnumbers of 1yo males with social data
birdlist16 = read.csv(here::here("Input files","birdlist16noKsubc.csv"))
birdlist17 = read.csv(here::here("Input files","birdlist17noKsubc.csv"))
birdlist18 = read.csv(here::here("Input files","birdlist18noKsubc.csv"))
birdlist19 = read.csv(here::here("Input files","birdlist19noKsubc.csv"))

##Load in capture data to make sure females were all banded before observations
captures = read.csv(here::here("Input files","All RBFW Captures.csv"))
captures$Date = as.Date(captures$Date,"%m/%d/%y")

##Load in 1yo paired male data
pag1p = read.csv(here::here("Output files","oneyo_paired.csv"))
#Implant males and males in removal communities have already been taken out. 
#These are all males that paired

##Get years when we have good social data 
pag1p = pag1p %>% filter(Year %in% c("2016","2017","2018","2019"))



####Setup dataframe####

##Subset to 1yo males with social data
birdlist16.s = birdlist16 %>% filter(Current.Age==1)
birdlist17.s = birdlist17 %>% filter(Current.Age==1)
birdlist18.s = birdlist18 %>% filter(Current.Age==1)
birdlist19.s = birdlist19 %>% filter(Current.Age==1)
Fwnumbers.social = rbind(data.frame(Fwnumber = birdlist16.s$FWnumber),data.frame(Fwnumber = birdlist17.s$FWnumber),
                         data.frame(Fwnumber = birdlist18.s$FWnumber),data.frame(Fwnumber = birdlist19.s$FWnumber))
pag1p = pag1p %>% filter(Fwnumber %in% Fwnumbers.social$Fwnumber)


##Drop colors column - comes from database which shows each birds last color combo - and add colors that year my matching to agesex file
pag1p = pag1p %>% select(-Colors)
agesex16.c = agesex16 %>% select(FWNo,Bird,Year)
agesex17.c = agesex17 %>% select(FWNo,Bird,Year)
agesex18.c = agesex18 %>% select(FWNo,Bird,Year)
agesex19.c = agesex19 %>% select(FWNo,Bird,Year)
agesex.c = rbind(agesex16.c,agesex17.c,agesex18.c,agesex19.c)
pag1p = merge(pag1p,agesex.c,by.x=c("Fwnumber","Year"),by.y=c("FWNo","Year"))


##Get the ID of the male's female - first female paired to in that season - already have Group ID for first breeding group where male is breeder in pag1p
#Also get group size and IDs of any other birds in breeding group
bgroups = rbind(b_groups16,b_groups17,b_groups18,b_groups19)
pag1p$Female = NA
pag1p$Bgroupsize = NA
pag1p$Other1 = NA
pag1p$Other2 = NA
for (i in 1:nrow(pag1p)) {
  bgroup = bgroups %>% filter(Group.Number==pag1p$Group[i])
  bgroupF = bgroup %>% select(Female)
  pag1p$Female[i] = bgroupF[[1]]
  bgroup.m = melt(id.vars=1:3,bgroup)
  bgroup.m = bgroup.m %>% filter(value!="")
  pag1p$Bgroupsize[i] = nrow(bgroup.m)
  pag1p$Other1[i] = bgroup.m$value[3]
  pag1p$Other2[i] = bgroup.m$value[4]
}
#RYR - techs were not confident in assigning his female, whether they really were paired or not - remove him in next step


##Remove UNK or UNB females
pag1p = pag1p %>% filter(!Female %in% c("UNB","UNK",""))


##Change GZR's mate to HIB from BWR - NOTES:
##For GZR, the database lists him as paired to his mother (BWR, age 5), but behaviorally it appeared he was paired with HIB a 1+min female
#before eventually pairing with RVR older female next door. GZR, BWR, and HIB were in non-breeding group together, in sightings have 
#instance written down that GZR and HIB were chasing RVR and duetting, then GZR switched to duetting with RVR and during the whole thing BWR
#did not take part. Think that shows that GZR was paired with HIB and not his mother. His mother eventually paired with LZB when LZB appeared
GZR = pag1p %>% filter(Bird=="GZR")
pag1p = pag1p %>% filter(Bird!="GZR")
GZR$Female = "HIB"
GZR$Other1 = "BWR"
pag1p = rbind(pag1p,GZR)
#Checked other groups in database to make sure they were accurate.


####Make sure all of the females were banded before observations took place, or like HIB were very easy to assign once banded
#Birds that were already banded will be in individual files as having been seen before or early in observations
captures.colors = captures %>% filter(Colors!="") #Get capture records with color bands
female.first = data.frame(female=pag1p$Female,year=pag1p$Year,first.colors=as.Date(NA),first.seen=as.Date(NA))

for (i in 1:nrow(female.first)) {
  female.fwn = agesex.c %>% filter(Bird==as.character(female.first$female[i]),Year==female.first$year[i]) %>% select(FWNo)
  captures.colors.f = captures.colors %>% filter(Fwnumber==female.fwn[1,]) %>% arrange(Date) %>% slice(1) 
  female.first$first.colors[i]=captures.colors.f$Date
  if(female.first$year[i]=="2016") {ind = ind16}
  if(female.first$year[i]=="2017") {ind = ind17}
  if(female.first$year[i]=="2018") {ind = ind18}
  if(female.first$year[i]=="2019") {ind = ind19}
  ind.female = ind %>% filter(Bird==as.character(female.first$female[i])) %>% arrange(Date) %>% slice(1)
  female.first$first.seen[i]=ind.female$Date
}

##Some birds we were able to backdate due to them being the only unbanded in an area when we banded them
#So get earliest date - banded or seen for each female, and if not before present/banded before week 28 of year
female.firstm = melt(female.first,id.vars=1:2)
female.firstm = female.firstm %>% arrange(value) %>% distinct(female,.keep_all = T)

##Flag females with first dates that might cause observation problems - week 28 is beginning of July
female.firstm$flag = NA
for (i in 1:nrow(female.firstm)) {
  value.year = year(female.firstm$value[i])
  if((female.firstm$year[i]-value.year)>1) {female.firstm$flag[i]="ok"} else {
    if(month(female.firstm$value[i])<7) {female.firstm$flag[i]="ok"} else {female.firstm$flag[i]="check"}
  }
}

##Check sightings to see if young male was ever seen consistently with an UNB during observations, if not, then ok, if he was, 
#then remove the male since we can't be sure when he started associating with the female. 
female.firstm.check = female.firstm %>% filter(flag=="check")

#2016:
#GRL of GRY and GRL - GRL banded/first seen 7/26 - GRY not with UNBs hardly at all before then, spent his time with his group
#RHB of WLB and RHB - RHB banded/first seen 8/1 - WLB spent his time with his groups, very few interactions with UNBs in early season in large kerfuffles
#VBY of YIV and VBY - VBY banded/first seen 10/23 - YIV with his group throughout July and August, did have 2 sampling points with group and UNB in middle, but only with group after that
#YVV of WIZ and YVV - YVV banded/first seen 12/13 - WIZ with his group and very few others throughout July and August

#2017:
#YRI of IWW and YRI - YRI banded/first seen 8/19 - YRI was a metal banded bird before capture and IWW consistently seen with metal banded birds during non-breeding season - TAKE OUT
#BRI of GGG and BRI - BRI banded/first seen 10/9 - GGG with his group throughout July and August, no UNBs

#2018:
#HHB of WGL and HHB - HHB banded/first seen 9/11 - WGL did have one observation with multiple groups present where there was an UNB in early July but from there on out it was just him and his group until he paired with HHB in Sept

#2019:
#ZBB of VLR and ZBB - ZBB banded/first seen 8/8 - VLR only with his group before pairing with ZBB
#YZI of LVB and YZI - YZI banded/first seen 9/11 - LVB only with his group throughout July and August
#ZRY of RYL and ZRY - ZRY banded/first seen 9/14 - RYL with group in late July/August but did interact in large kerfuffles with an UNB included in early July quite often, TAKE OUT to be safe
#ZWR of LIB and ZWR - ZWR banded/first seen 9/22 - LIB in same group therefore same boat as ZRY, so TAKE OUT to be safe
#BBV of WWR and BBV - BBV banded/first seen 9/28 - WWR with group throughout July/August

#Females to REMOVE:
remove.females = c("YRI","ZRY","ZWR")
pag1p = pag1p %>% filter(!Female %in% remove.females)


###Get relatedness between male and female
rel = read.csv(here::here("Input files","RelatednessEstimates_8xdenovo80_14to19.csv"))
rel2 = rel %>% select(Pair,ind2,ind1,wang) #Get columns in different order 
colnames(rel2)[2:3] = c("ind1","ind2")
rel = rbind(rel,rel2) #Combine so all individuals in both columns

pag1p$MFrel = NA
for (i in 1:nrow(pag1p)) {
  female.fwn = agesex.c %>% filter(Bird==pag1p$Female[i],Year==pag1p$Year[i]) %>% select(FWNo)
  rel.sub = rel %>% filter(ind1==pag1p$Fwnumber[i] & ind2==female.fwn$FWNo) %>% select(wang)
  if(nrow(rel.sub)>0) {pag1p$MFrel[i] = rel.sub[[1]]} else {pag1p$MFrel[i]=NA}
}

#Probably only one mother-offspring pair - RWI and BLL. Although BRG and YRY were at 0.32. Looking at database no obvious connection
#other than location. YRY is much older than BRG. Could be half siblings, not obvious mother-son relationship. They were in the same
#non-breeding group and were together from the beginning of the season. BRG is more related to a different older female - actually BLL
#at 0.43. She is probably his mother and he shared a parent with YRY. 

##Mike recommended removing RWI from dataset since he was paired to his mother and will be hard to tell when they "paired"
#pag1p = pag1p %>% filter(Bird != "RWI")
##Tried taking this bird out but the model ended up failing when I did that - residuals were very poor. Think it makes sense to 
#leave him in. He was paired eventually, we're just uncertain as to when that happened since he was with his mom the whole time. 
#He still had the opportunity to take on bright plumage. 





####Reclassify ornamented or not based on whether it came from adventitious molt or not####
pag1p.orn = pag1p %>% filter(Ornamented=="Ornamented")
ind.orn = data.frame(Bird=NA,jdate=NA,Molt=NA)

for (i in 1:nrow(pag1p.orn)) {
  if(pag1p.orn$Year[i]=="2016") {ind=ind16}
  if(pag1p.orn$Year[i]=="2017") {ind=ind17}
  if(pag1p.orn$Year[i]=="2018") {ind=ind18}
  if(pag1p.orn$Year[i]=="2019") {ind=ind19}
  ind.male = ind %>% filter(Bird==as.character(pag1p.orn$Bird[i]))
  ind.male$jdate = yday(ind.male$Date)
  ind.male = ind.male %>% select(Bird,jdate,Molt)
  ind.orn = rbind(ind.orn,ind.male)
}
ind.orn = ind.orn[-1,]
ind.orn = ind.orn %>% filter(Molt!="")

##Plot to see difference in true molt vs adventitious molt
ggplot(data=ind.orn,aes(x=jdate,y=Molt,group=Bird,color=Bird)) + geom_point() + geom_line()
#GRY's sighting at 40% came in November from one observer in one sighting. Could be true or could be a mistake, but either way it clearly
#looks like a different molt trajectory (adventitious molt) than the others. Classify him as brown in ornamented.early column. 
pag1p$Ornamented.early = pag1p$Ornamented
GRY = pag1p %>% filter(Bird=="GRY")
pag1p = pag1p %>% filter(Bird!="GRY")
GRY$Ornamented.early = "Brown"
pag1p = rbind(pag1p,GRY)





####Plot weighted degree among pairs by week####

#Get week data in bird files
bird16$week = week(bird16$Date)
bird17$week = week(bird17$Date)
bird18$week = week(bird18$Date)
bird19$week = week(bird19$Date)

#For loop for plots
wdeg.weeks = data.frame(week=NA,bird=NA,wdeg=NA,seen=NA) #Create empty dataframe to rbind to in for loop
pag1p = pag1p %>% arrange(as.character(Bird))

par(mfrow=c(3,3))
for (i in 1:nrow(pag1p)) {
  year = pag1p$Year[i]
  male = pag1p$Bird[i]
  female = pag1p$Female[i]
  if(year=="2016") {bird = bird16} 
  if(year=="2017") {bird = bird17} 
  if(year=="2018") {bird = bird18} 
  if(year=="2019") {bird = bird19} 
  weeks = data.frame(seq(24,48))
  colnames(weeks) = "week"
  weeks$bird = pag1p$Bird[i]
  weeks$wdeg = NA
  weeks$seen = NA
  for (h in 1:nrow(weeks)) {
    bird.week = bird %>% filter(week==weeks$week[h])
    bird.male = bird.week %>% filter(bird.week$Bird==as.character(male))
    bird.female = bird.week %>% filter(bird.week$Bird==as.character(female))
    X = nrow(bird.male %>% filter(bird.male$Sighting %in% bird.female$Sighting))
    Ya = nrow(bird.male %>% filter(!bird.male$Sighting %in% bird.female$Sighting))
    Yb = nrow(bird.female %>% filter(!bird.female$Sighting %in% bird.male$Sighting))
    weeks$wdeg[h] = X/(X+Ya+Yb)
    weeks$seen[h] = X+Ya+Yb
  }
  weeks.print = weeks %>% filter(!is.na(wdeg)) %>% filter(seen>=3)
  plot(weeks.print$week,weeks.print$wdeg,ylim=c(0,1),xlim=c(24,48),pch=19,main=paste(male,pag1p$Ornamented.early[i]),
       xlab="Week of year",ylab="W-degree")
  lines(weeks.print$week,weeks.print$wdeg)
  wdeg.weeks = rbind(wdeg.weeks,weeks)
}
par(mfrow=c(1,1))

#Should be noted when looking at these graphs that GZR paired with HIB, then quickly switched mates, hence the decrease in wdeg to his paired female

#Remove the first row from wdeg.weeks - it was blank from setup
wdeg.weeks = wdeg.weeks[-1,]

###Look at what time of year I have enough social data to compare among weeks
wdeg.weeks.noNA = wdeg.weeks %>% filter(!is.na(wdeg)) %>% filter(seen>=3)
ggplot(data=wdeg.weeks.noNA,aes(x=week,y=bird)) + geom_point() + geom_line() 

#Almost all males have enough data to use weeks 28-32 and used linear interpolation to fill in missing weeks






####Sliding window approach for pairing and molt#### 

#Want to see if wdeg connections to female at certain dates influence whether a 1yo male molts into nuptial plumage

##Get first and last data points for each individual
m.first = wdeg.weeks.noNA %>% arrange(week) %>% distinct(bird,.keep_all = T) %>% select(-wdeg,-seen)
colnames(m.first)[1] = "first.week"
m.last = wdeg.weeks.noNA %>% arrange(desc(week)) %>% distinct(bird,.keep_all = T) %>% select(-wdeg,-seen)
colnames(m.last)[1] = "last.week"

#Merge with dataframe
pag1p = merge(pag1p,m.first,by.x="Bird",by.y="bird")
pag1p = merge(pag1p,m.last,by.x="Bird",by.y="bird")


##Interpolate missing data when individuals missing points
wdeg.weeks.int = data.frame(week=NA,bird=NA,wdeg=NA,wdeg.int=NA)
for (i in 1:nrow(pag1p)) {
  wdeg.weeks.ind = wdeg.weeks %>% filter(bird==as.character(pag1p$Bird[i])) %>% select(-seen)
  wdeg.weeks.ind = wdeg.weeks.ind %>% filter(week>=pag1p$first.week[i])
  wdeg.weeks.ind = wdeg.weeks.ind %>% filter(week<=pag1p$last.week[i])
  wdeg.weeks.ind$wdeg.int = NA
  if (nrow(wdeg.weeks.ind %>% filter(is.na(wdeg)))==0) {wdeg.weeks.int = rbind(wdeg.weeks.int,wdeg.weeks.ind)} else {
    wdeg.weeks.ind.int = data.frame(approx(x=wdeg.weeks.ind$week,y=wdeg.weeks.ind$wdeg,n=nrow(wdeg.weeks.ind)))
    wdeg.weeks.ind$wdeg.int = wdeg.weeks.ind.int$y
    wdeg.weeks.int = rbind(wdeg.weeks.int,wdeg.weeks.ind)
  }
}
wdeg.weeks.int = wdeg.weeks.int[-1,]

##Did interpolation work? - only wdeg = NA or wdeg.int = NA rows should not be test=TRUE - looks good
interpolation.test = wdeg.weeks.int
interpolation.test$test = wdeg.weeks.int$wdeg==wdeg.weeks.int$wdeg.int
nrow(interpolation.test %>% filter(!is.na(wdeg),!is.na(wdeg.int),test==FALSE))

##For birds without missing data, move values to int column
for (i in 1:nrow(wdeg.weeks.int)) {if(is.na(wdeg.weeks.int$wdeg.int[i])) {wdeg.weeks.int$wdeg.int[i]=wdeg.weeks.int$wdeg[i]}}

##Get Ornamented.early as factor
pag1p$Ornamented.early = as.factor(pag1p$Ornamented.early)




####Function for sliding windows
pairing_window = function(range.start,range.end,wdeg.weeks.int,pag1p) {
  
  ##For the specified range, get the individuals that have data within the range 
  #Need to have wdeg measures for on or before and on or after the range limits
  pag1pw = pag1p %>% filter(first.week<=range.start)
  pag1pw = pag1pw %>% filter(last.week>=range.end)
  pag1pw$Year = as.factor(pag1pw$Year)
  
  #Get all possible week windows
  windows = expand.grid(start=seq(range.start,range.end),end=seq(range.start,range.end))
  windows$diff = windows$start-windows$end #Get difference between windows
  windows = windows %>% filter(diff<=0) %>% select(-diff)
  
  #Get aic value of base model
  aic.base = AIC(glm(Ornamented.early ~ Year,data=pag1pw,family="binomial"))
  
  #Make output dataframe
  output = data.frame(matrix(ncol=6,nrow=nrow(windows)))
  colnames(output) = c("start.week","end.week","ModelBeta","deltaAIC","Range.start","Range.end")
  
  #For each possible window
  for (i in 1:nrow(windows)) {
    #Get a mean wdeg variable for each individual in pag1pw
    pag1pw$mean.wdeg = NA
    #For each row in pag1pw
    for (h in 1:nrow(pag1pw)) {
      #Get wdeg data for that male
      wdeg.weeks.intm = wdeg.weeks.int %>% filter(bird==pag1pw$Bird[h]) 
      #Filter for weeks in window
      wdeg.weeks.intm = wdeg.weeks.intm %>% filter(week>=windows$start[i],week<=windows$end[i])
      #Put the mean value in pag1pw
      pag1pw$mean.wdeg[h]=mean(wdeg.weeks.intm$wdeg.int)
    }
    
    #Get the model with wdeg included
    model.wdeg = glm(Ornamented.early ~ mean.wdeg + Year,data=pag1pw,family="binomial")
    
    #Get AIC
    aic.wdeg = AIC(model.wdeg)
    
    #Add values to output dataframe
    output$start.week[i] = windows$start[i]
    output$end.week[i] = windows$end[i]
    output$ModelBeta[i] = coef(model.wdeg)[2]
    output$deltaAIC[i] = aic.wdeg - aic.base
    output$Range.start[i] = range.start
    output$Range.end[i] = range.end
  }
  
  #Order dataframe by lowest AIC values compared to baseline
  output = output %>% arrange(deltaAIC)

  #Output dataframes
  return(output)
}



###Run the function at a couple different ranges 
pairing.output.28.32 = pairing_window(range.start = 28,range.end = 32,wdeg.weeks.int = wdeg.weeks.int,pag1p = pag1p)
pairing.output.27.32 = pairing_window(range.start = 27,range.end = 32,wdeg.weeks.int = wdeg.weeks.int,pag1p = pag1p)
pairing.output.27.33 = pairing_window(range.start = 27,range.end = 33,wdeg.weeks.int = wdeg.weeks.int,pag1p = pag1p)
pairing.output.27.35 = pairing_window(range.start = 27,range.end = 35,wdeg.weeks.int = wdeg.weeks.int,pag1p = pag1p)

#No matter the range, 30-30 is the best window. Cannot compare deltaAIC values across runs since datasets are different but can look
#at what window comes out on top. 
library(climwin)
pairing.output.28.32.climwin = pairing.output.28.32
colnames(pairing.output.28.32.climwin) = c("WindowOpen","WindowClose","ModelBeta","deltaAICc")
plotdelta(dataset=pairing.output.28.32.climwin)




###Function for randomizations to make sure best window is not a false poisitive

#Since using a sliding window approach and testing multiple fixed effects versions, likelihood of getting a false positive is high,
#so need to make sure that observed deltaAIC value is lower than expected by chance. Making sure that a meaningful window 
#could not be discovered by random chance. 
pairing_window_rand = function(range.start,range.end,wdeg.weeks.int,pag1p,randomizations) {
  
  ##For the specified range, get the individuals that have data within the range 
  #Need to have wdeg measures for on or before and on or after the range limits
  pag1pw = pag1p %>% filter(first.week<=range.start)
  pag1pw = pag1pw %>% filter(last.week>=range.end)
  pag1pw$Year = as.factor(pag1pw$Year)
  
  #Get all possible week windows
  windows = expand.grid(start=seq(range.start,range.end),end=seq(range.start,range.end))
  windows$diff = windows$start-windows$end #Get difference between windows
  windows = windows %>% filter(diff<=0) %>% select(-diff)
  
  #Get aic value of base model
  aic.base = AIC(glm(Ornamented.early ~ Year,data=pag1pw,family="binomial"))
  
  #Make output dataframe
  output = data.frame(matrix(ncol=1,nrow=randomizations))
  colnames(output) = "aic.rand"
  
  for (k in 1:randomizations) {
    #Randomize the bird ID to mix up wdeg and ornamented.early pairings - keeps the same pairings between year and ornamented.early
    pag1pw$Bird.rand = sample(pag1pw$Bird,replace=F)
    #Get a deltaAIC variable for windows
    windows$deltaAIC = NA
    
    #For each possible window
    for (i in 1:nrow(windows)) {
      #Get a mean wdeg variable for each individual in pag1pw
      pag1pw$mean.wdeg = NA
      #For each row in pag1pw
      for (h in 1:nrow(pag1pw)) {
        #Get wdeg data for that male
        wdeg.weeks.intm = wdeg.weeks.int %>% filter(bird==pag1pw$Bird.rand[h]) 
        #Filter for weeks in window
        wdeg.weeks.intm = wdeg.weeks.intm %>% filter(week>=windows$start[i],week<=windows$end[i])
        #Put the mean value in pag1pw
        pag1pw$mean.wdeg[h]=mean(wdeg.weeks.intm$wdeg.int)
      }
      
      #Get the model with wdeg included
      model.wdeg = glm(Ornamented.early ~ mean.wdeg + Year,data=pag1pw,family="binomial")
      
      #Get AIC
      windows$deltaAIC[i] = AIC(model.wdeg) - aic.base
    }
    #Get best deltaAIC from all windows and put into output for that randomization
    windows = windows %>% arrange(deltaAIC)
    output$aic.rand[k] = windows$deltaAIC[1]
  }
  #Output dataframes
  return(output)
}


##Run 100 randomizations and compare observed to random AIC values
#pairing.output.rand.28.32 = pairing_window_rand(range.start = 28,range.end = 32,wdeg.weeks.int,pag1p,randomizations=100)
ggplot(data=pairing.output.rand.28.32,aes(x=aic.rand)) + geom_histogram(color="black",fill="gray") + theme_cowplot() +
  geom_vline(aes(xintercept = pairing.output.28.32$deltaAIC[1]),color="red",size=1) + xlab("deltaAIC")

##Get p-value
obs.AIC = pairing.output.28.32$deltaAIC[1]
sum(obs.AIC>pairing.output.rand.28.32$aic.rand)/100
#p=0.01

#So finding such an important observed window by random is quite unlikely, so ok to move forward with analyses looking at importance of that window
#in predicting molt into nuptial plumage

#Top window actually produces a categorical model, where in this dataset, all individuals that were ornamented early were always seen 
#with their females in week 30

#Apply this same window to the full dataset, using all individuals that we had data on at week 30. Sort of like training the model then applying
#it to a larger dataset. Found the sliding window with individuals that had data within the entire range of interest then can apply that window
#to a larger dataset. 



####Function for a single window
pairing_window_single = function(range.start,range.end,start.week,end.week,wdeg.weeks.int,pag1p) {
  
  ##For the specified start-end range, get the individuals that have data within the range of the single window
  wwint.start = wdeg.weeks.int %>% filter(week==start.week)
  wwint.end = wdeg.weeks.int %>% filter(week==end.week)
  pag1pw = pag1p %>% filter(Bird %in% wwint.start$bird)
  pag1pw = pag1p %>% filter(Bird %in% wwint.end$bird)
  pag1pw$Year = as.factor(pag1pw$Year)
  
  #Get a mean wdeg variable for each individual in pag1pw
  pag1pw$mean.wdeg = NA
  
  #For each row in pag1pw
  for (h in 1:nrow(pag1pw)) {
    #Get wdeg data for that male
    wdeg.weeks.intm = wdeg.weeks.int %>% filter(bird==pag1pw$Bird[h]) 
    #Filter for weeks in window
    wdeg.weeks.intm = wdeg.weeks.intm %>% filter(week>=start.week,week<=end.week)
    #Put the mean value in pag1pw
    pag1pw$mean.wdeg[h]=mean(wdeg.weeks.intm$wdeg.int)
  }
  
  #Get the model with wdeg included
  model.wdeg = glm(Ornamented.early ~ mean.wdeg + Year,data=pag1pw,family="binomial")
  
  #Output dataframes
  return(list(model.wdeg,pag1pw))
}

#Get single window model - standard errors are so large here due to reference level 2016 having all zeros for response variable 
pag1pw.30.30m = pairing_window_single(range.start=28,range.end=32,start.week = 30,end.week = 30,wdeg.weeks.int = wdeg.weeks.int,pag1p = pag1p)[[1]]
summary(pag1pw.30.30m)
pag1pw.30.30m.coef = coef(pag1pw.30.30m)[2]

#Get single window dataframe
pag1pw.30.30 = pairing_window_single(range.start=28,range.end=32,start.week = 30,end.week = 30,wdeg.weeks.int = wdeg.weeks.int,pag1p = pag1p)[[2]]

#See how you can relevel to get less extreme standard errors in the model
pag1pw.30.30b = pag1pw.30.30
pag1pw.30.30b$Year = relevel(pag1pw.30.30b$Year,ref="2018")
pag1pw.30.30mb = glm(Ornamented.early ~ mean.wdeg + Year,data=pag1pw.30.30,family="binomial")
summary(pag1pw.30.30mb)

#Run model without year and see that extreme standard errors are gone
pag1pw.30.30mc = glm(Ornamented.early ~ mean.wdeg,data=pag1pw.30.30,family="binomial")
summary(pag1pw.30.30mc)
pag1pw.30.30mc.coef = coef(pag1pw.30.30mc)[2]


###Check residuals of model
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = pag1pw.30.30m)
plot(simulationOutput, asFactor = F) #Plot residual plots 
plotResiduals(simulationOutput,pag1pw.30.30$mean.wdeg, quantreg = T) #Plot residuals for 
testDispersion(simulationOutput)
#Look good



###Randomizations to test importance of mean.wdeg to female in models at 30 week window
##Keeps same pairings between year and ornamented early
pag1pw.30.30.rand = function(coef.obs, perm.number) {
  rands = data.frame(matrix(ncol=1,nrow=perm.number))
  colnames(rands) = "rands"
  for (i in 1:perm.number) {
    pag1pw.30.30$mean.wdeg.rand = sample(pag1pw.30.30$mean.wdeg,replace=F)
    model.rand = glm(Ornamented.early ~ mean.wdeg.rand + Year,data=pag1pw.30.30,family="binomial")
    rands$rands[i] = coef(model.rand)[2]
  }
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30,fill="gray",color="black") + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficient values"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
}

##Run randomization
#pag1pw.30.30.rand(coef.obs = pag1pw.30.30m.coef,perm.number = 1000)
#p<0.001 or p=0.001


####Paper: #Find that the best association to female window to explain likelihood of molt was 30 weeks, and this window was not a false positive.
#Males that were closely associated with their female by the 30th week of the year were more likely to molt into nuptial plumage than those that 
#were not (coef=7.33,p=0.001).


####Plots

#Get ornamented early data in wdeg.weeks dataframe
orn.early = pag1p %>% select(Bird,Ornamented.early)
wdeg.weeks.plot = merge(wdeg.weeks.noNA,orn.early,by.x="bird",by.y="Bird")

##On 8/26 of 2018, GZR switched females from HIB to RVR - I saw it happen. RVR old female's mate had disappeared and GZR and HIB were 
#duetting a bunch on the territory boundary against her, then after a couple minutes of GZR and HIB chasing RVR, GZR and RVR began chasing
#HIB and then GZR and RVR were duetting. The next day GZR and RVR were together and not with HIB. 
#So 8/26 is end of week 34, 8/27 is beginning of week 35, remove GZR-HIB data from week 35 for this plot since he was no longer paired to 
#her at that point. 
wdeg.weeks.plot = wdeg.weeks.plot %>% filter(!(bird=="GZR" & week==35)) 

##Get means and standard errors for each group
wdeg.weeks.plotm = wdeg.weeks.plot %>% group_by(Ornamented.early,week) %>% summarise(mean(wdeg)) #Get mean wdeg by week for groups
colnames(wdeg.weeks.plotm)[3] = "wdeg" 
library(plotrix) #std.error package
wdeg.weeks.plots = wdeg.weeks.plot %>% group_by(Ornamented.early,week) %>% summarise(std.error(wdeg)) #Get std.error by week for groups
colnames(wdeg.weeks.plots)[3] = "std.error" #Rename column
wdeg.weeks.plotms = data.frame(wdeg.weeks.plotm,wdeg.weeks.plots$std.error) #Combine mean and std.error dataframes
colnames(wdeg.weeks.plotms)[4] = "std.error" #Rename column
wdeg.weeks.plotms$Ornamented.early = factor(wdeg.weeks.plotms$Ornamented.early,levels=c("Ornamented","Brown")) #Order factor for plotting

##Plot mean weighted degree to paired female and std error by week of year for ornamented and brown males
ggplot(data=wdeg.weeks.plotms,aes(x=week,y=wdeg)) + theme_cowplot() + geom_errorbar(aes(ymin=wdeg-std.error,ymax=wdeg+std.error),color="black") +
  xlab("Week in year") + ylab("Mean weighted degree to paired female") + scale_x_continuous(breaks=c(25,30,35,40,45)) +
  geom_point(size=3.5,color="black") + geom_point(size=2.5,aes(color=Ornamented.early)) + scale_color_manual(values=c("black","gray")) +
  labs(color="Breeding\\nphenotype") + guides(color = guide_legend(override.aes = list(size=3.5,linetype=0)))

##Plot again with colors
ggplot(data=wdeg.weeks.plotms,aes(x=week,y=wdeg)) + theme_cowplot(font_size = 12) + geom_errorbar(aes(ymin=wdeg-std.error,ymax=wdeg+std.error),color="black") +
  xlab("Week in year") + ylab("Mean association score to\\nfuture mate") + scale_x_continuous(breaks=c(25,30,35,40,45)) +
  geom_point(size=3,color="black") + geom_point(size=2,aes(color=Ornamented.early)) + scale_color_manual(values=c("#D82F2B","#C6A375")) +
  labs(color="Breeding\\nphenotype") + guides(color = guide_legend(override.aes = list(size=3.5,linetype=0)))





###Prediction plot - normal ben bolker script didn't work here. Both get same prediction values, couldn't test if got 
#same CI values. 

#Get ornamented.early in binary form for plotting
pag1pw.30.30$Ornamented.early.binary = NA
for (i in 1:nrow(pag1pw.30.30)) {if(pag1pw.30.30$Ornamented.early[i]=="Brown") {pag1pw.30.30$Ornamented.early.binary[i]=0} else {
  pag1pw.30.30$Ornamented.early.binary[i]=1}
}

#Get prediction dataframe 
pag1pw.30.30m.df = data.frame(expand.grid(mean.wdeg=seq(0,1,by=0.05),Year=as.factor(2019)))

#Predict values 
ilink <- family(pag1pw.30.30m)$linkinv
pd <- cbind(pag1pw.30.30m.df, predict(pag1pw.30.30m, pag1pw.30.30m.df, type = "link", se.fit = TRUE)[1:2])
pd <- transform(pd, Fitted = ilink(fit), Upper = ilink(fit + (2 * se.fit)),
                Lower = ilink(fit - (2 * se.fit)))

ggplot(data=pag1pw.30.30,aes(x=mean.wdeg,y=Ornamented.early.binary)) + geom_point(position=position_jitter(width = 0,height=0.03),size=2,alpha=0.75) +
  theme_cowplot(font_size = 12) + geom_ribbon(data = pd, aes(ymin = Lower, ymax = Upper, x = mean.wdeg), alpha = 0.1, inherit.aes = FALSE) +
  geom_line(data = pd, aes(y = Fitted, x = mean.wdeg),size=1) + xlim(-0.1,1.1) + 
  scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1)) + scale_y_continuous(breaks=c(0,1),labels=c("0","1")) + 
  ylab("Likelihood of breeding in\\nornamented plumage") + xlab("Association score to future mate at week 30")




####Check condition likelihood of ornamentation relationship within these males####

##Read in condition data
cond = read.csv(here::here("Output files","Condition scores 15_19.csv"))
cond.s = cond %>% select(Year,Fwnumber,jdate,mtresid)
cond.s$Year = as.factor(cond.s$Year)

##Merge captures for the same individual
cond.s = cond.s %>% group_by(Year,Fwnumber) %>% summarise(mtresid = mean(mtresid),mean.jdate = mean(jdate))

##Merge 
pag1pw.cond = merge(pag1pw.30.30,cond.s,by=c("Year","Fwnumber"))
pag1pw.cond$Fwnumber = as.factor(pag1pw.cond$Fwnumber)

###First question - did the males that became ornamented differ in condition from those that did not?
ggplot(data=pag1pw.cond,aes(x=Ornamented.early,y=mtresid)) + geom_boxplot() + geom_point()
condm = lm(mtresid~Ornamented.early + Year,data=pag1pw.cond)
summary(condm)

#Check residuals
library(car)
hist(resid(condm))
qqPlot(condm)
#Look good

##Likelihood ratio test #Need to get X2 value from sum of squares
condm.oe = lm(mtresid~Year,data=pag1pw.cond)
anova(condm.oe,condm)

##Paper: Males in this group that acquired ornamented plumage did not differ in condition from those that did not
#acquire ornamented plumage (F=2.82,p=0.105)

###Second question - does condition explain why some males acquired ornamented plumage and others did not? Add to connection to mate model
##Model - adding mean.jdate didn't change anything, year was singular
condm2 = glm(Ornamented.early ~ mean.wdeg + mtresid + Year,data=pag1pw.cond,family="binomial")
summary(condm2)

###Check residuals of model 
simulationOutput <- simulateResiduals(fittedModel = condm2)
plot(simulationOutput, asFactor = F) #Plot residual plots 
testDispersion(simulationOutput)
#Look ok

##Likelihood ratio test
condm2.mt = glm(Ornamented.early ~ mean.wdeg + Year,data=pag1pw.cond,family="binomial")
anova(condm2.mt,condm2,test="LRT")

##Paper: Within these males, residual condition did not influence the likelihood of ornamented plumage acquisition (X2=0.05,p=0.821)





####Plot timing of pairing/association and timing of molt together####

###Need a plot to show that pairing typically comes before molt. Have 
#6 males that were ornamented and had good social data, need to calculate approximate
#timing of pairing for those six birds. 

##Select ornamented 1yos
pair.molt = pag1pw.30.30 %>% filter(Ornamented.early=="Ornamented")

##Get first dates seen as ornamented. For 2018 birds, I don't trust the plumage scores by LGH. I assigned very different
#scores on the same day or soon after for two of these males. So take out those sightings for this graph. 
#Note - these daily plumage scores for 1-year-olds are not used in any other analyses so this problem doesn't affect anything else.
pair.molt18 = pair.molt %>% filter(Year==2018)
pair.molt19 = pair.molt %>% filter(Year==2019)
ind18.pm = ind18 %>% filter(!Initials=="LGH") %>% filter(Bird %in% pair.molt18$Bird) %>% arrange(Date) %>% 
  filter(Molt>33) %>% distinct(Bird,.keep_all = T)
ind19.pm = ind19 %>% filter(Bird %in% pair.molt19$Bird) %>% arrange(Date) %>% 
  filter(Molt>33) %>% distinct(Bird,.keep_all = T) 
ind18.pmd = ind18.pm %>% select(Bird,Date) %>% mutate(pair.molt="First date\\n male seen in\\nornamented plumage")
ind19.pmd = ind19.pm %>% select(Bird,Date) %>% mutate(pair.molt="First date\\n male seen in\\nornamented plumage")
ind.pmd = rbind(ind18.pmd,ind19.pmd)

##Get first date breeding pair seen alone or in breeding group. Similar for loop as above to get w-deg by week
pair.alone = matrix(nrow=0,ncol=2) #Empty matrix
colnames(pair.alone) = c("Bird","alone.date") 
for (i in 1:nrow(pair.molt)) {
  year = pair.molt$Year[i]
  male = pair.molt$Bird[i]
  female = pair.molt$Female[i]
  if(year=="2018") {bird = bird18} 
  if(year=="2019") {bird = bird19} 
  bird.male = bird %>% filter(bird$Bird==as.character(male))
  bird.female = bird %>% filter(bird$Bird==as.character(female))
  bird.both = bird.male %>% filter(bird.male$Sighting %in% bird.female$Sighting) %>% filter(Number==2) %>% arrange(Date) %>%
    slice_head() %>% select(Bird,Date) %>% mutate(pair.molt="First date\\npair seen\\nalone")
  pair.alone = rbind(pair.alone,bird.both)
  }

##Combine dataframes
pair.molt.plot = rbind(pair.alone,ind.pmd)
pair.molt.plot$jdate = yday(pair.molt.plot$Date)

##Reorder factors to plot
pair.molt.plot$pair.molt = factor(pair.molt.plot$pair.molt,levels=c("First date\\npair seen\\nalone","First date\\n male seen in\\nornamented plumage"))

##For ploting break up into multiple dataframes
pair.molt.plot2 = pair.molt.plot %>% filter(Bird %in% c("RRY","BRG"))
pair.molt.plot = pair.molt.plot %>% filter(!Bird %in% c("RRY","BRG"))

##Plot
ggplot(data=pair.molt.plot,aes(x=pair.molt,y=jdate)) + #geom_boxplot(aes(fill=pair.molt)) + 
  geom_line(aes(x=pair.molt,y=jdate,group=Bird),size=0.8) + geom_point(size=3,color="black") + geom_point(size=2,aes(color=pair.molt)) + 
  geom_line(data=pair.molt.plot2, aes(x=pair.molt,y=jdate,group=Bird),size=0.8) + geom_point(data=pair.molt.plot2,size=3,color="black") + 
  geom_point(data=pair.molt.plot2,size=2,aes(color=pair.molt)) + ylim(165,260) +
  theme_cowplot() + theme(legend.position = "") + xlab("") + ylab("Day of year") + scale_color_manual(values=c("#C6A375","#D82F2B")) 

##Plot alongside other plot
A = ggplot(data=pair.molt.plot,aes(x=pair.molt,y=jdate)) + #geom_boxplot(aes(fill=pair.molt)) + 
  geom_line(aes(x=pair.molt,y=jdate,group=Bird),size=0.8) + geom_point(size=3,color="black") + geom_point(size=2,aes(color=pair.molt)) + 
  geom_line(data=pair.molt.plot2, aes(x=pair.molt,y=jdate,group=Bird),size=0.8) + geom_point(data=pair.molt.plot2,size=3,color="black") + 
  geom_point(data=pair.molt.plot2,size=2,aes(color=pair.molt)) + ylim(165,260) +
  theme_cowplot() + theme(legend.position = "") + xlab("") + ylab("Day of year") + scale_color_manual(values=c("#C6A375","#D82F2B"))
B = ggplot(data=wdeg.weeks.plotms,aes(x=week,y=wdeg)) + theme_cowplot(font_size = 12) + geom_errorbar(aes(ymin=wdeg-std.error,ymax=wdeg+std.error),color="black") +
  xlab("Week in year") + ylab("Mean association score to\\nfuture mate") + scale_x_continuous(breaks=c(25,30,35,40,45)) +
  geom_point(size=3,color="black") + geom_point(size=2,aes(color=Ornamented.early)) + scale_color_manual(values=c("#D82F2B","#C6A375")) +
  labs(color="Breeding\\nphenotype") + guides(color = guide_legend(override.aes = list(size=3.5,linetype=0)))

ggarrange(A,B,align="h")



####Social data####

#Want to test whether non-breeding social connections influence a 1-year-old male's plumage
#Use network node permutations to test significance. Using a datastream (pre-network) permutation would not work 
#well here because we already know the network is very structured, more structured than random. 
#Double permutation method is for measuring a trait's influence on a network measure, not a network measure on a trait. 


###Load in non-breeding social connections data
one_social = read.csv(here::here("Output files","oneyo_social_connections_observed.csv"))

##Drop unnecessary columns from one_social 
one_social = one_social %>% select(-Sex,-Year,-Current.Age,-Age.Exact,-Bird,)

##Merge social data with previous dataset
pag1pw.social = merge(pag1pw.30.30,one_social,by.x="Fwnumber",by.y="FWnumber")

##Are all males with high wdeg to female in a non-breeding group with their female?
ggplot(data=pag1pw.social,aes(x=Group_with_paired_female,y=mean.wdeg)) + geom_boxplot() + theme_cowplot() + geom_point(size=2)
#One male that was not - VHR - was seen with his female in a large interaction in week 30, then he didn't pair with her until a bit later

##Get males that were with in the same non-breeding group as their female
pag1pw.sg = pag1pw.social %>% filter(Group_with_paired_female=="Yes")
#Make sure these are the same individuals on the far right of the previous plot
ggplot(data=pag1pw.sg,aes(x=mean.wdeg,y=Ornamented.early)) + geom_point(position=position_jitter(width = 0.05,height=0),size=3,alpha=0.75) + 
  xlim(0,1.1)

##Did these males differ in condition? - Does not appear so
pag1pw.sgcond = merge(pag1pw.sg,cond.s,by=c("Year","Fwnumber"))
ggplot(data=pag1pw.sgcond,aes(x=Ornamented.early,y=mtresid)) + geom_boxplot() + geom_point()



####See if social relationships can explain why some of these males that are in a group with their non-breeding female molted into nuptial plumage
##and others did not



####In a group with another male or not - use larger dataset for this - here male can have molted before or after social data
pag1pw.30.30$another_male_in_group = NA
pag1pw.30.30$Female_in_nbgroup = NA
for (i in 1:nrow(pag1pw.30.30)) {
  year=pag1pw.30.30$Year[i]
  if(year=="2016") {birdlist = birdlist16} 
  if(year=="2017") {birdlist = birdlist17} 
  if(year=="2018") {birdlist = birdlist18} 
  if(year=="2019") {birdlist = birdlist19} 
  sg = birdlist %>% filter(Bird==as.character(pag1pw.30.30$Bird[i])) %>% select(Social.Group)
  sgs = birdlist %>% filter(Social.Group==sg[1,])
  sgf = sgs %>% filter(Bird==pag1pw.30.30$Female[i])
  if(nrow(sgf)==1) {pag1pw.30.30$Female_in_nbgroup[i]="Yes"} else 
  {pag1pw.30.30$Female_in_nbgroup[i]="No"}
  sgm = birdlist %>% filter(Social.Group==sg[1,]) %>% filter(Sex=="M")
  if(nrow(sgm)>1) {pag1pw.30.30$another_male_in_group[i]="Yes"} else 
  {pag1pw.30.30$another_male_in_group[i]="No"}
  }

##Filter for males in a non-breeding group with their female
pag1pw.30.30.nbg = pag1pw.30.30 %>% filter(Female_in_nbgroup=="Yes")

##Look at a table of results
table(pag1pw.30.30.nbg$Ornamented.early,pag1pw.30.30.nbg$another_male_in_group)
#Males that took on ornamented plumage were in a non-breeding social group with their female and never in a non-breeding group with another male

##Fisher exact test to look at these ratios
fisher.test(x=pag1pw.30.30.nbg$another_male_in_group,y=pag1pw.30.30.nbg$Ornamented.early)

##Paper: This could have been due to chance (p=0.182)



####Now for males that we have social data before they molted, look to see if their social interactions differed before molt
##Sample sizes here are probably too small to do binomial models, just run linear models to see if groups were different


####Weighted degree models####


####Weighted degree to other birds

#Get weighted degree to other birds other than the male's female
pag1pw.sg$wdeg_notfemale = pag1pw.sg$wdeg_all - pag1pw.sg$wdeg_female_paired
#Plot
ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_notfemale)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to individuals\\nother than male's female") + theme(legend.position = "")
#Model
pag1pw.sgm1 = lm(wdeg_notfemale~Ornamented.early + Year,data=pag1pw.sg)
summary(pag1pw.sgm1) 
pag1pw.sgm1.coef = coef(pag1pw.sgm1)[2]

#Check residuals
library(car)
hist(resid(pag1pw.sgm1))
qqPlot(pag1pw.sgm1)
#Look good


##Randomization to test for significance
pag1pw.sgm1.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag1pw.sg$Ornamented.early.rand = sample(pag1pw.sg$Ornamented.early,replace=F)
    rands[i] = coef(lm(wdeg_notfemale~Ornamented.early.rand + Year,data=pag1pw.sg))[2]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = 1 - sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  return(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
           geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
}

#Run permutation - 10000 permutations gets pretty reliable (unchanging) p-value
#B = pag1pw.sgm1.rand(perm.number = 10000,coef.obs=pag1pw.sgm1.coef)

##Paper: Males that took on ornamented plumage had lower wdeg to birds other than their female (coef=-1.38, p=0.019)

##Plot together
A = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_notfemale)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to individuals\\nother than male's female") + theme(legend.position = "")
ggarrange(A,B,nrow=1)






####Weighted degree to other males

#Get weighted degree to other birds other than the male's female
pag1pw.sg$wdeg_males = pag1pw.sg$wdeg_1yo_males + pag1pw.sg$wdeg_old_males
#Plot
ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_males)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to other males") + theme(legend.position = "")
#Model
pag1pw.sgm2 = lm(wdeg_males~Ornamented.early + Year,data=pag1pw.sg)
summary(pag1pw.sgm2) 
pag1pw.sgm2.coef = coef(pag1pw.sgm2)[2]

#Check residuals
library(car)
hist(resid(pag1pw.sgm2))
qqPlot(pag1pw.sgm2)
#Look ok


##Randomization to test for significance
pag1pw.sgm2.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag1pw.sg$Ornamented.early.rand = sample(pag1pw.sg$Ornamented.early,replace=F)
    rands[i] = coef(lm(wdeg_males~Ornamented.early.rand + Year,data=pag1pw.sg))[2]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = 1 - sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  return(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
           geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
}

#Run permutation - 10000 permutations gets pretty reliable (unchanging) p-value
#B = pag1pw.sgm2.rand(perm.number = 10000,coef.obs=pag1pw.sgm2.coef)

##Paper: Males that took on ornamented plumage had lower wdeg to other males (coef=-1.12, p=0.013)

##Plot together
A = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_males)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to other males") + theme(legend.position = "")
ggarrange(A,B,nrow=1)





####Weighted degree to old males

#Plot
ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_old_males)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to old males") + theme(legend.position = "")
#Model
pag1pw.sgm3 = lm(wdeg_old_males~Ornamented.early + Year,data=pag1pw.sg)
summary(pag1pw.sgm3) 
pag1pw.sgm3.coef = coef(pag1pw.sgm3)[2]

#Check residuals
library(car)
hist(resid(pag1pw.sgm3))
qqPlot(pag1pw.sgm3)
#Look ok


##Randomization to test for significance
pag1pw.sgm3.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag1pw.sg$Ornamented.early.rand = sample(pag1pw.sg$Ornamented.early,replace=F)
    rands[i] = coef(lm(wdeg_old_males~Ornamented.early.rand + Year,data=pag1pw.sg))[2]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = 1 - sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  return(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
           geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
}

#Run permutation - 10000 permutations gets pretty reliable (unchanging) p-value
#B = pag1pw.sgm3.rand(perm.number = 10000,coef.obs=pag1pw.sgm3.coef)

##Paper: Males that took on ornamented plumage had a lower wdeg to old males (coef=-0.29, p=0.042)

##Plot together
A = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_old_males)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to old males") + theme(legend.position = "")
ggarrange(A,B,nrow=1)






####Weighted degree to 1yo males

#Plot
ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_1yo_males)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to 1-year-old males") + theme(legend.position = "")
#Model
pag1pw.sgm4 = lm(wdeg_1yo_males~Ornamented.early + Year,data=pag1pw.sg)
summary(pag1pw.sgm4) 
pag1pw.sgm4.coef = coef(pag1pw.sgm4)[2]

#Check residuals
library(car)
hist(resid(pag1pw.sgm4))
qqPlot(pag1pw.sgm4)
#Look ok


##Randomization to test for significance
pag1pw.sgm4.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag1pw.sg$Ornamented.early.rand = sample(pag1pw.sg$Ornamented.early,replace=F)
    rands[i] = coef(lm(wdeg_1yo_males~Ornamented.early.rand + Year,data=pag1pw.sg))[2]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = 1 - sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  return(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
           geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency") )
}

#Run permutation - 10000 permutations gets pretty reliable (unchanging) p-value
#B = pag1pw.sgm4.rand(perm.number = 10000,coef.obs=pag1pw.sgm4.coef)

##Paper: Males that took on ornamented plumage had fewer connections to 1yo males (coef=-0.83, p=0.048)

##Plot together
A = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_1yo_males)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to 1-year-old males") + theme(legend.position = "")
ggarrange(A,B,nrow=1)






####Weighted degree to females in different groups

#Plot
ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_females_diffgrps)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to females in different groups") + theme(legend.position = "")
#Model
pag1pw.sgm5 = lm(wdeg_females_diffgrps~Ornamented.early + Year,data=pag1pw.sg)
summary(pag1pw.sgm5) 
pag1pw.sgm5.coef = coef(pag1pw.sgm5)[2]

#Check residuals
library(car)
hist(resid(pag1pw.sgm5))
qqPlot(pag1pw.sgm5)
#Look ok


##Randomization to test for significance
pag1pw.sgm5.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag1pw.sg$Ornamented.early.rand = sample(pag1pw.sg$Ornamented.early,replace=F)
    rands[i] = coef(lm(wdeg_females_diffgrps~Ornamented.early.rand + Year,data=pag1pw.sg))[2]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = 1 - sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  return(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
           geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
}

#Run permutation - 10000 permutations gets pretty reliable (unchanging) p-value
#B = pag1pw.sgm5.rand(perm.number = 10000,coef.obs=pag1pw.sgm5.coef)

##Paper: Males that took on ornamented plumage had a trend towards fewer connections to females in different groups - so not molting necessarily based on 
#number of interactions with females (coef=-0.26, p=0.057)

##Plot together
A = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_females_diffgrps)) + geom_boxplot(aes(fill=Ornamented.early)) + 
  geom_point() + theme_cowplot() + ylim(0,7) + scale_fill_manual(values=c("#C6A375","#D82F2B")) + xlab("Plumage") +
  ylab("Weighted degree to females in\\ndifferent groups") + theme(legend.position = "")
ggarrange(A,B,nrow=1)





####Plot all boxplots and randomization plots for supplemental at once####

A = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_notfemale)) + geom_boxplot(aes(fill=Ornamented.early)) + geom_point() + 
  theme_cowplot() + ylim(0,7) + xlab("Plumage") + ylab("Weighted degree to all individuals\\nbut the male's female") + 
  scale_fill_manual(values=c("light gray","#cc7032")) + theme(legend.position = "None") + xlab("")

B = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_males)) + geom_boxplot(aes(fill=Ornamented.early)) + geom_point() + 
  theme_cowplot() + ylim(0,7) + xlab("Plumage") + ylab("Weighted degree to males") + 
  scale_fill_manual(values=c("light gray","#cc7032")) + theme(legend.position = "None") + xlab("")

C = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_old_males)) + geom_boxplot(aes(fill=Ornamented.early)) + geom_point() + 
  theme_cowplot() + ylim(0,7) + xlab("Plumage") + ylab("Weighted degree to old males") + 
  scale_fill_manual(values=c("light gray","#cc7032")) + theme(legend.position = "None") + xlab("")

D = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_1yo_males)) + geom_boxplot(aes(fill=Ornamented.early)) + geom_point() + 
  theme_cowplot() + ylim(0,7) + xlab("Plumage") + ylab("Weighted degree to young males") + 
  scale_fill_manual(values=c("light gray","#cc7032")) + theme(legend.position = "None") + xlab("")

E = ggplot(data=pag1pw.sg,aes(x=Ornamented.early,y=wdeg_females_diffgrps)) + geom_boxplot(aes(fill=Ornamented.early)) + geom_point() + 
  theme_cowplot() + ylim(0,7) + xlab("Plumage") + ylab("Weighted degree to females\\n in different groups") + 
  scale_fill_manual(values=c("light gray","#cc7032")) + theme(legend.position = "None") + xlab("")
ggarrange(A,B,C,D,E)


##Will have to do randomization plots manually in AD

