####When do 2-year-olds molt script


##Note to self - in ic_sp plots cannot calculate CIs - see survCIs documentation in icenReg: 
#This function is not compatible with ic_np or ic_sp objects, as the distribution of the baseline distribution of these estimators 
#is still an open question.

##Using interval censoring in this script - Interval censoring occurs when a response is known only up to an interval. So in 
#icenReg dataframe, males that were first seen in bright plumage are time 1 = 0, time 2 = seen bright. Know that event happened 
#sometime in that interval but not sure exactly when. 


####Load data####

library(here)
library(dplyr)
library(reshape2)
library(lubridate)
library(ggplot2)
library(cowplot)
library(viridis)

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
ind15 = read.csv(here::here("Input files","Individuals2015 7_13_20.csv"))
ind15$Molt = as.numeric(as.character(ind15$Molt))
ind15$Date = as.Date(ind15$Date,"%m/%d/%y")
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
agesex15 = read.csv(here::here("Input files","Ages and Sexes 2015.csv"))
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

##Load in capture data 
captures = read.csv(here::here("Input files","All RBFW Captures.csv"))
captures$Date = as.Date(captures$Date,"%m/%d/%y")

##Load in capture data with plumage scores
captures.p = read.csv(here::here("Input files","All RBFW Captures plumage scores.csv"))
captures.p$Date = as.Date(captures.p$Date,"%m/%d/%y")

##Load in 2yo male ornamented data
pag2 = read.csv(here::here("Output files","twoyo_ornamented.csv"))
#Implant males and males in removal communities have already been taken out. 
#These are all males that molted into intermediate or nuptial plumage

##Get 2015-2019 - years when we were there during the non-breeding season to get molt dates
pag2 = pag2 %>% filter(Year %in% c("2015","2016","2017","2018","2019"))

##Match up colors to bird using agesex files - this will be that individuals colors in the year of interest pag2$Year and will match to ind file
agesex15.s = agesex15 %>% select(Year,FWNo,Bird)
agesex16.s = agesex16 %>% select(Year,FWNo,Bird)
agesex17.s = agesex17 %>% select(Year,FWNo,Bird)
agesex18.s = agesex18 %>% select(Year,FWNo,Bird)
agesex19.s = agesex19 %>% select(Year,FWNo,Bird)
agesexall = rbind(agesex15.s,agesex16.s,agesex17.s,agesex18.s,agesex19.s)
pag2 = merge(pag2,agesexall,by.x=c("Fwnumber","Year"),by.y=c("FWNo","Year"))

##Make sure no birds changed color combos within a year
captures.p.test = captures.p %>% filter(Fwnumber %in% pag2$Fwnumber)
captures.p.test = captures.p.test %>% distinct(Year,Fwnumber,Colors)
check = captures.p.test %>% group_by(Year,Fwnumber) %>% filter(n()>1)
#So in captures, need to change those birds to the combo they are in agesex files for that year 8730 = BRV, 7776 = RGG
check1 = captures.p %>% filter(Year==2017 & Colors=="RGX") 
check1$Colors = "RGG"
check2 = captures.p %>% filter(Year==2019 & Colors=="BRG")
check2$Colors = "BRV"
captures.p = rbind(check1,check2,captures.p)
captures.p = captures.p %>% distinct(Fwnumber,Date,.keep_all = T) %>% filter(Fwnumber %in% pag2$Fwnumber) %>% arrange(Date) #Remove later values 

##Make sure capture colors are the same as colors in agesex files for that year
capture.colorcheck = pag2 %>% select(Fwnumber, Year, Bird)
capture.colorcheck$colors.same = NA
for (i in 1:nrow(capture.colorcheck)) {
  capture.colors = captures.p %>% filter(Year==capture.colorcheck$Year[i],Fwnumber==capture.colorcheck$Fwnumber[i])
  if(nrow(capture.colors)>0) {if(capture.colors$Colors[1]==as.character(capture.colorcheck$Bird[i])) {capture.colorcheck$colors.same[i]="Yes"}
  }
} #All either yes or NA so captures are good. NA is bird was not captured that year. 




####Set up molt date dataframe####

##Get columns of interest
ind15.s = ind15 %>% select("Bird","Date","Molt") %>% mutate(method="sighting")
ind16.s = ind16 %>% select("Bird","Date","Molt") %>% mutate(method="sighting")
ind17.s = ind17 %>% select("Bird","Date","Molt") %>% mutate(method="sighting")
ind18.s = ind18 %>% select("Bird","Date","Molt") %>% mutate(method="sighting")
ind19.s = ind19 %>% select("Bird","Date","Molt") %>% mutate(method="sighting")

##Combine ind list with capture plumages
captures.p.15 = captures.p %>% filter(Year==2015) %>% mutate(Bird=Colors,Molt = Brightness,method="capture") %>% select(Bird,Date,Molt,method)
ind15.s = rbind(ind15.s,captures.p.15)
captures.p.16 = captures.p %>% filter(Year==2016) %>% mutate(Bird=Colors,Molt = Brightness,method="capture") %>% select(Bird,Date,Molt,method)
ind16.s = rbind(ind16.s,captures.p.16)
captures.p.17 = captures.p %>% filter(Year==2017) %>% mutate(Bird=Colors,Molt = Brightness,method="capture") %>% select(Bird,Date,Molt,method)
ind17.s = rbind(ind17.s,captures.p.17)
captures.p.18 = captures.p %>% filter(Year==2018) %>% mutate(Bird=Colors,Molt = Brightness,method="capture") %>% select(Bird,Date,Molt,method)
ind18.s = rbind(ind18.s,captures.p.18)
captures.p.19 = captures.p %>% filter(Year==2019) %>% mutate(Bird=Colors,Molt = Brightness,method="capture") %>% select(Bird,Date,Molt,method)
ind19.s = rbind(ind19.s,captures.p.19)


##Get birds of interest - birds in pag2
pag2.15 = pag2 %>% filter(Year=="2015")
ind15.s = ind15.s %>% filter(Bird %in% pag2.15$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2015")
pag2.16 = pag2 %>% filter(Year=="2016")
ind16.s = ind16.s %>% filter(Bird %in% pag2.16$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2016")
pag2.17 = pag2 %>% filter(Year=="2017")
ind17.s = ind17.s %>% filter(Bird %in% pag2.17$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2017")
pag2.18 = pag2 %>% filter(Year=="2018")
ind18.s = ind18.s %>% filter(Bird %in% pag2.18$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2018")
pag2.19 = pag2 %>% filter(Year=="2019")
ind19.s = ind19.s %>% filter(Bird %in% pag2.19$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2019")

##Combine, order, and get julian date
ind.all = rbind(ind15.s,ind16.s,ind17.s,ind18.s,ind19.s)
ind.all = ind.all %>% arrange(Year,Bird,Date)
ind.all$jdate = yday(ind.all$Date)

##If date is in January, make Julian date +365
for (i in 1:nrow(ind.all)) {
  if(month(ind.all$Date[i])==1) {ind.all$jdate[i] = ind.all$jdate[i] + 365}
}






####Look at timing of molt across years####

##Plot molt by date
ggplot(data=ind.all,aes(x=jdate,y=Molt)) + geom_point(aes(color=Year),size=1) + geom_line(aes(group=Bird,color=Year),size=0.7) + 
  theme_cowplot() + scale_color_viridis(discrete = T)
#2015 and 2017 might be noticeably later than other years

##Take out the two intermediate males from 2017 
ind.alls = ind.all %>% filter(!Bird %in% (pag2 %>% filter(Plumage=="Interm") %>% select(Bird))$Bird)

##Plot again
ggplot(data=ind.alls,aes(x=jdate,y=Molt)) + geom_point(aes(color=Year),size=1) + geom_line(aes(group=Bird,color=Year),size=0.7) + 
  theme_cowplot() + scale_color_viridis(discrete = T)

##Get week data - add one if in January - one week and two week intervals
ind.alls$week = week(ind.alls$Date)
for (i in 1:nrow(ind.alls)) {
  if(ind.alls$week[i]==1) {ind.alls$week[i]=54}}
ind.alls$twoweek = ceiling(ind.alls$week/2) *2 #two week intervals - divide week by 2 and round up, then multiply by two to get later week of 2-week interval

##For these plots - use sighting data only, capture data muddy the averages a little bit since they don't match sightings exactly
ind.alls.nocap = ind.alls %>% filter(method=="sighting")

##Get one sighting per week for each individual - best pluamge score of that week
duplicated(pag2$Bird) #Below code works - distinct for bird and year because no duplicate color combos in this dataset
ind.alls.d = ind.alls.nocap %>% arrange(desc(Molt)) %>% distinct(Bird,twoweek,.keep_all = T) #Get one sighting of each bird per two week interval

##Group by year and week and get average plumage score and get count of individuals seen in that week
ind.alls.week = ind.alls.d %>% group_by(Year,twoweek) %>% summarise(Molt=mean(Molt))
ind.alls.weekn = ind.alls.d %>% group_by(Year,twoweek) %>% count(twoweek)
ind.alls.week$count = ind.alls.weekn$n

##Require at least 5 individuals present to keep the week in 
ind.alls.week = ind.alls.week %>% filter(count>=5)

##Plot year means
ggplot(data=ind.alls.week,aes(x=twoweek,y=Molt)) + geom_line(aes(color=Year),size=1.8) + 
  theme_cowplot() + scale_color_viridis(discrete = T) + xlab("Week of year") + ylab("Mean 2-year-old male plumage score") + 
  geom_point(size=2,aes(color=Year))
#Some major differences across years in when 2yo males molted







####Look at timing of molt within males####

###How regular is timing of molt data?
##Plot for individual males - jdate
ggplot(data=ind.alls,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=Year),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")

##Plot for individual males - jdate but color points by method
ggplot(data=ind.alls,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=method),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_manual(values=c("red","black")) +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2)
#Including capture data does help bridge the gap between some individual's sightings, but once plumages are stable, it can muddy the trajectory a little

##Plot for individual males - week
ggplot(data=ind.alls,aes(x=week,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=Year),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T)

##Plot for individual males - twoweeks
ggplot(data=ind.alls,aes(x=twoweek,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=Year),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T)


###Measure difference in time between last dull observation and first intermediate/first bright
##Filter pag2 for only birds with molt data - and adventitious molt birds removed
pag2s = pag2 %>% filter(Bird %in% ind.alls$Bird) 

##For loop to get difference in dates
pag2s$first.seen = NA
pag2s$times.seen = NA
pag2s$last.dull = NA
pag2s$first.int = NA
pag2s$last.int = NA
pag2s$first.bright = NA
pag2s$first.plum.afterdull = NA
pag2s$first.bright.plum = NA
pag2s$first.plum
for (i in 1:nrow(pag2s)) {
  bird = pag2s$Bird[i]
  ind.bird = ind.alls %>% filter(Bird==as.character(bird)) %>% arrange(jdate)
  pag2s$first.plum[i] = ind.bird$Molt[1]
  pag2s$first.seen[i] = ind.bird$jdate[1] 
  pag2s$times.seen[i] = nrow(ind.bird)
  ind.dull = ind.bird %>% filter(Molt<33) %>% arrange(desc(jdate))
  pag2s$last.dull[i] = ind.dull$jdate[1]
  ind.int = ind.bird %>% filter(Molt>=33,Molt<=66) %>% arrange(jdate)
  pag2s$first.int[i] = ind.int$jdate[1]
  ind.int2 = ind.int %>% arrange(desc(jdate))
  pag2s$last.int[i] = ind.int2$jdate[1]
  ind.bright= ind.bird %>% filter(Molt>66) %>% arrange(jdate)
  pag2s$first.bright[i] = ind.bright$jdate[1]
  pag2s$first.bright.plum[i] = ind.bright$Molt[1]
  ind.firstplumAD = ind.bird %>% filter(Molt>=33) %>% arrange(jdate)
  pag2s$first.plum.afterdull[i] = ind.firstplumAD$Molt[1]
}

##Get difference in dates
pag2s$dull.to.bright = pag2s$first.bright-pag2s$last.dull
pag2s$dull.to.int = pag2s$first.int-pag2s$last.dull
pag2s$int.to.bright = pag2s$first.bright-pag2s$last.int
pag2s$fint.to.bright = pag2s$first.bright-pag2s$first.int
pag2s$dull.to.intorbright = NA
for (i in 1:nrow(pag2s)) {if(!is.na(pag2s$dull.to.int[i])) {
  pag2s$dull.to.intorbright[i] = pag2s$dull.to.int[i]} else {pag2s$dull.to.intorbright[i] = pag2s$dull.to.bright[i]}
  }


##For a male that makes it to intermediate, how long does it take to get to bright? 
#Get individuals that have an intermediate plumage sighting, then whose first bright plumage sighting was not fully bright
pag2s.test = pag2s %>% filter(first.plum.afterdull<=66,first.bright.plum<=90)
ind.alls.test = ind.alls %>% filter(Bird %in% pag2s.test$Bird) #%>% filter(Molt>=33,Molt<=90)
ggplot(data=ind.alls.test,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=Year),size=1.2) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T,name="Season") + xlim(135,300) + 
  geom_hline(aes(yintercept=66),lty=2,color="gray") + xlab("Day of year") + ylab("Plumage score") + ylim(0,100)
ggplot(pag2s.test,aes(x=fint.to.bright)) + geom_histogram(fill="gray",color="black",binwidth = 2) + theme_cowplot() + 
  xlab("Days between first intermediate plumage (33-66%)\\nand first red-black plumage (67-90%)") + ylab("Count") + xlim(0,30)
#Appears for the most part that once a male made it to intermediate, he made it to bright very quickly afterwards 
#So base analyses around the first date seen ornamented (intermediate or bright). Should have those data for more individuals 
#than just first date seen bright. 
mean(pag2s.test$fint.to.bright)

##Since molt from intermediate to bright so fast, and all of these males eventually made it to bright, if a male has
#an intermediate plumage score, then leave him in. 


##Since molt from intermediate to lower bright (70-85) was so fast, if a male does not have an intermediate plumage score
#but does have a lower bright score within 20 days of his last dull sighting, leave him in. This also based on 
#no or very few individuals maintained a bright plumage of 70-85? They all made it to 90 or greater and stayed there? 
#Correct - could look up most common brightest plumage for each individual 


##Look at most common bright plumage for each individual
pag2s$bright.mostseen = NA
for (i in 1:nrow(pag2s)) {
  bird = pag2s$Bird[i]
  ind.bright = ind.alls %>% filter(Bird==as.character(bird)) %>% filter(Molt>66) %>% arrange(jdate) %>% distinct(jdate,.keep_all = T)
  ind.bright.table = data.frame(table(ind.bright$Molt)) %>% arrange(desc(Freq))
  pag2s$bright.mostseen[i] = as.numeric(as.character(ind.bright.table$Var1[1]))
}
ggplot(pag2s,aes(x=bright.mostseen)) + geom_histogram(color="black",fill="gray",binwidth = 1) + theme_cowplot() +
  xlab("Most common bright plumage for each male")
#Some males plateaued at 90 and 95, most at 100, but none less than 90% red-black. So that means any observation between 70 and 85 should
#be in within the quick molting trajectory timeline. 



####Filter based on above findings
##Select birds with an intermediate plumage or early bright plumage less than 90 - we would have seen these males during molt and molt is usually very fast
pag2s.fpd = pag2s %>% filter(first.plum.afterdull<90)

##For males that went from dull to very bright quickly, restrict to individuals seen less than 20 days between last dull and first bright
#Look
pag2s.dvb = pag2s %>% filter(!first.plum.afterdull<90) %>% filter(!is.na(dull.to.bright))
ggplot(pag2s.dvb,aes(x=dull.to.bright)) + geom_histogram(color="black",fill="gray",binwidth = 5) + theme_cowplot()
ind.alls.test2 = ind.alls %>% filter(Bird %in% pag2s.dvb$Bird)
ggplot(data=ind.alls.test2,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=Year),size=1.2) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) + xlim(135,300) + 
  geom_hline(aes(yintercept=66),lty=2,color="gray") + xlab("Day of year") + ylab("Plumage score") + ylim(0,100) +
  geom_hline(aes(yintercept=33),lty=2,color="gray")
#Restrict
pag2s.dvb = pag2s.dvb %>% filter(dull.to.bright<=20)
ind.alls.test3 = ind.alls %>% filter(Bird %in% pag2s.dvb$Bird)
ggplot(data=ind.alls.test3,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=Year),size=1.2) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) + xlim(135,300) + 
  geom_hline(aes(yintercept=66),lty=2,color="gray") + xlab("Day of year") + ylab("Plumage score") + ylim(0,100) +
  geom_hline(aes(yintercept=33),lty=2,color="gray")
#Looks better

##For males that were first seen as bright - 90 or above
#Look
pag2s.fsb = pag2s %>% filter(first.plum>66)
ind.alls.test4 = ind.alls %>% filter(Bird %in% pag2s.fsb$Bird)
ggplot(data=ind.alls.test4,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=Year),size=1.2) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) + xlim(135,300) + 
  geom_hline(aes(yintercept=66),lty=2,color="gray") + xlab("Day of year") + ylab("Plumage score") + ylim(0,100) +
  geom_hline(aes(yintercept=33),lty=2,color="gray")
ggplot(data=pag2s,aes(x=first.seen,y=Year)) + geom_line(aes(group=Year)) + geom_point() + 
  geom_point(data=pag2s.fsb,aes(x=first.seen,y=Year),color="red")
#Those all look decent, can have them be left-censored in the survival analysis 

###Combine
pag2ss = rbind(pag2s.fpd,pag2s.dvb,pag2s.fsb)
pag2ss = pag2ss %>% distinct(Fwnumber,.keep_all = T) #Remove duplicates - found RLY from 2016 was in two of the categories above
ind.allss = ind.alls %>% filter(Bird %in% pag2ss$Bird)

###Plot trajectories
##Plot for individual males - jdate
ggplot(data=ind.allss,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=Year),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")

##Plot for individual males - jdate but color points by method
ggplot(data=ind.allss,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=method),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_manual(values=c("red","black")) +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2)

##Get first int or bright plumage - molt date
pag2ss$molt.date = NA
for (i in 1:nrow(pag2ss)) {
  if(!is.na(pag2ss$first.int[i])) {pag2ss$molt.date[i] = pag2ss$first.int[i]} else {pag2ss$molt.date[i] = pag2ss$first.bright[i]}
}

##Plot molt date colors
pag2ss.moltdate = pag2ss %>% select(Bird,molt.date)
ind.allss.plot = ind.allss
ind.allss.plot = merge(ind.allss.plot,pag2ss.moltdate,by="Bird")
ind.allss.plot$molt.date.bird = paste(ind.allss.plot$molt.date,ind.allss.plot$Bird)
ind.allss.plot$Bird.year = paste(ind.allss.plot$Bird,ind.allss.plot$Year)
birdorder = ind.allss.plot %>% arrange(molt.date) %>% select(Bird.year) %>% distinct()
ind.allss.plot$Bird.year = factor(ind.allss.plot$Bird.year,levels=birdorder$Bird.year)

#Plot with bird and year labels
ggplot(data=ind.allss.plot,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=molt.date),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(Bird.year)) + 
  scale_color_viridis(discrete = F, begin=1,end=0,trans="reverse","Molt date") +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")

#Plot with molt date and bird labels
ggplot(data=ind.allss.plot,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=molt.date),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(molt.date.bird)) + 
  scale_color_viridis(discrete = F, begin=1,end=0,trans="reverse","Molt date") +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")

#Plot on one graph, color by molt date
ggplot(data=ind.allss.plot,aes(x=jdate,y=Molt))  + geom_line(size=1,aes(group=Bird,color=molt.date)) +  geom_point(aes(color=molt.date),size=1) +
  theme_cowplot() + scale_color_viridis(discrete = F) +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")

#Plot on one graph, color by year
ggplot(data=ind.allss.plot,aes(x=jdate,y=Molt))  + geom_line(size=1,aes(group=Bird,color=Year)) +  geom_point(aes(color=Year),size=1) +
  theme_cowplot() + scale_color_viridis(discrete = T) +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")





####Survival analyses to look at timing of molt####

library(survival)
library(survminer)
library(icenReg)

##Set up columns for survival analysis - if left censored (first seen as bright so don't know true molt date), then time1 = infinity (NA) and
#time2 = date first seen bright. If we know the true date (or have an approximation), then time1 and time 2 are the same
#Survival package will not run coxph on interval-censored data. icenReg will - make left-censored data time1 = 0 according to paper - Anderson-Bergman
pag2ss$molt.time1 = NA
pag2ss$molt.time2 = NA
for (i in 1:nrow(pag2ss)) {
  if(pag2ss$first.plum[i] >= 90) {
    pag2ss$molt.time1[i]=0
    pag2ss$molt.time2[i]=pag2ss$molt.date[i]} else {
      pag2ss$molt.time1[i]=pag2ss$molt.date[i]
      pag2ss$molt.time2[i]=pag2ss$molt.date[i]
    }
}
pag2ss.op = pag2ss %>% filter(!is.na(Ornamented.prev)) 
pag2ss.op$Year = as.factor(pag2ss.op$Year)


##Survival models to use - Anderson-Bergman paper says that semi-parametric is better for interval censored data than parametric where
#results can be heavily influenced by model choice and interval censored results are hard to diagnose. Semi-parametric works with coxph
#and includes some of the metrics of non-parametric models that make them better for interval-censored data. 

###Write pag2ss.op to .csv for use in climate sliding window analyses
#write.csv(pag2ss.op,here::here("Output files","twoyo_moltdates_dataframe.csv"),row.names = F)

##Write pag2ss for plotting rainfall and timing of molt plots
#write.csv(pag2ss,here::here("Output files","twoyo_moltdates_all_dataframe.csv"),row.names = F)


##Plot individuals that end up getting used in climate sliding windows - those that were not first seen in full nuptial plumage
pag2.sdm = pag2ss.op %>% filter(molt.time1>0)
ind.allss.plot.sdm = ind.allss.plot %>% filter(Bird %in% pag2.sdm$Bird)

ggplot(data=ind.allss.plot.sdm,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=molt.date),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(Bird.year)) + 
  scale_color_viridis(discrete = F, begin=1,end=0,trans="reverse","Molt date") +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")


####Whether the male was ornamented previously and paired previously####

##Does being ornamented previously influence timing of molt for 2yo males? 
#1000 or more bootstraps are best if using Wald p-values
#Not sure if can do a likelihood ratio test here. survival does implement one for their coxph models, but not sure if anything might be different
#for semi-parametric models. Best to use permutation tests - better than wald p-values. Do not need bootstrap for that, bootstrap only 
#used for standard error and p-values. Permutations compare coefficients of variable of interest. 

##Run both proportional hazards and proportional odds models without the Paired.prev and Year variables to first check model diagnostics
cph.ph <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.op,model="ph") #,bs_samples=1000)
cph.po <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.op,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph) #Same for both models
diag_covar(cph.ph,varName = "Ornamented.prev")
diag_covar(cph.po,varName = "Ornamented.prev")
#About the same

##Run both proportional hazards and proportional odds models without the Ornamented.prev and Year variables to first check model diagnostics
cph.ph <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev, data = pag2ss.op,model="ph") #,bs_samples=1000)
cph.po <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev, data = pag2ss.op,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph) #Same for both models
diag_covar(cph.ph,varName = "Paired.prev")
diag_covar(cph.po,varName = "Paired.prev")
#About the same


###Full PH model 
cph.full <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev + Paired.prev + Year, data = pag2ss.op,model="ph") #,bs_samples=1000)
summary(cph.full)
cph.op.coef = coef(cph.full)[1]
cph.pp.coef = coef(cph.full)[2]

##Plot
plot.data = data.frame(Ornamented.prev=c("Ornamented","Brown"),Year="2019",Paired.prev="Yes")
rownames(plot.data) = c("Ornamented","Brown")
plot(cph.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")

plot.data = data.frame(Paired.prev=c("Yes","No"),Year="2019",Ornamented.prev="Ornamented")
rownames(plot.data) = c("Paired","Not paired")
plot(cph.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")


##Permutation test for whether a male was ornamented previously
cph.op.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.op$Ornamented.prev.rand = sample(pag2ss.op$Ornamented.prev,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev.rand + Paired.prev + Year, data = pag2ss.op,model="ph"))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#A = cph.op.rand(perm.number = 1000,coef.obs = cph.op.coef)

##Paper: When looking at all years, males that were ornamented previously were maybe slightly more likely to molt into nuptial plumage earlier than 
#males that were not ornamented as 1-year-olds (coef=0.55,p=0.068).



##Permutation test for whether a male was paired previously
cph.pp.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.op$Paired.prev.rand = sample(pag2ss.op$Paired.prev,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev + Paired.prev.rand + Year, data = pag2ss.op,model="ph"))[2]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#B = cph.pp.rand(perm.number = 1000,coef.obs = cph.pp.coef)

##Paper: When looking at all years, males that were paired previously were slightly more likely to molt into nuptial plumage earlier than 
#males that were not paired as 1-year-olds (coef=0.56,p=0.093).

##Plot both randomization plots
ggarrange(A,B,labels=c("A","B"))


###Plot all years ornamented previously
#Set paired.prev = yes, because all ornamented males were paired, so does not make sense to have a not-paired ornamented male in prediction

#Get ornamented data from model
plot.orn = expand.grid(Ornamented.prev="Ornamented",Paired.prev="Yes",Year=c(2015:2019),Day=130:365) 
plot.orn$Year = as.factor(plot.orn$Year)
est.orn = data.frame(getFitEsts(cph.full,newdata = plot.orn,q=plot.orn$Day),Year=plot.orn$Year,Day=plot.orn$Day)
colnames(est.orn)[1] = "est"

#Get brown data from model
plot.brown = expand.grid(Ornamented.prev="Brown",Paired.prev="Yes",Year=c(2015:2019),Day=130:365) 
plot.brown$Year = as.factor(plot.brown$Year)
est.brown =  data.frame(getFitEsts(cph.full,newdata = plot.brown,q=plot.brown$Day),Year=plot.brown$Year,Day=plot.brown$Day)
colnames(est.brown)[1] = "est"

#Combine
est.orn$Ornamented.prev="Ornamented"
est.brown$Ornamented.prev="Brown"
est.orn.brown = rbind(est.orn,est.brown)

#Plot
ggplot() + geom_line(data=est.orn,aes(x=Day,y=est,group=Year,color=Year),size=1) + 
  geom_line(data=est.brown,aes(x=Day,y=est,group=Year,color=Year),size=1) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_viridis(discrete = T)

ggplot(data=est.orn.brown,aes(x=Day,y=est)) + geom_line(aes(group=Ornamented.prev,color=Ornamented.prev),size=1.5) + 
  facet_wrap(facets=vars(Year)) + scale_color_manual(values=c("gray55","darkorange2"),name="Plumage in\\nprevious year") + theme_cowplot() +
  ylab("Likelihood of molt") + xlab("Day of year")



###Plot all years paired previously
#Set plumage = brown, because all ornamented males were paired, so does not make sense to have a not-paired ornamented male in prediction

#Get paired data from model
plot.pairy = expand.grid(Ornamented.prev="Brown",Paired.prev="Yes",Year=c(2015:2019),Day=130:365) 
plot.pairy$Year = as.factor(plot.pairy$Year)
est.pairy = data.frame(getFitEsts(cph.full,newdata = plot.pairy,q=plot.pairy$Day),Year=plot.pairy$Year,Day=plot.pairy$Day)
colnames(est.pairy)[1] = "est"

#Get brown data from model
plot.pairn = expand.grid(Ornamented.prev="Brown",Paired.prev="No",Year=c(2015:2019),Day=130:365) 
plot.pairn$Year = as.factor(plot.pairn$Year)
est.pairn =  data.frame(getFitEsts(cph.full,newdata = plot.pairn,q=plot.pairn$Day),Year=plot.pairn$Year,Day=plot.pairn$Day)
colnames(est.pairn)[1] = "est"

#Combine
est.pairy$Paired.prev="Yes"
est.pairn$Paired.prev="No"
est.pairy.pairn = rbind(est.pairy,est.pairn)

#Plot
ggplot() + geom_line(data=est.pairy,aes(x=Day,y=est,group=Year,color=Year),size=1) + 
  geom_line(data=est.pairn,aes(x=Day,y=est,group=Year,color=Year),size=1) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_viridis(discrete = T)

ggplot(data=est.pairy.pairn,aes(x=Day,y=est)) + geom_line(aes(group=Paired.prev,color=Paired.prev),size=1.5) + 
  facet_wrap(facets=vars(Year)) + scale_color_manual(values=c("gray55","darkorange2"),name="Paired in\\nprevious year") + theme_cowplot() +
  ylab("Likelihood of molt") + xlab("Day of year")


###Might be worth breaking these analyses up by wet and dry years - groups do seem different in wet years but not dry years, 
#similar to cockburn findings.







####Look at wet years separately####

pag2ss.op.wet = pag2ss.op %>% filter(Year %in% c("2016","2018","2019"))

##Run both proportional hazards and proportional odds models without the Paired.prev and Year variables to first check model diagnostics
cph.ph.wet1 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.op.wet,model="ph") #,bs_samples=1000)
cph.po.wet1 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.op.wet,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph.wet1) #Same for both models
diag_covar(cph.ph.wet1,varName = "Ornamented.prev")
diag_covar(cph.po.wet1,varName = "Ornamented.prev")
#About the same

##Run both proportional hazards and proportional odds models without the Ornamented.prev and Year variables to first check model diagnostics
cph.ph.wet2 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev, data = pag2ss.op.wet,model="ph") #,bs_samples=1000)
cph.po.wet2 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev, data = pag2ss.op.wet,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph.wet2) #Same for both models
diag_covar(cph.ph.wet2,varName = "Paired.prev")
diag_covar(cph.po.wet2,varName = "Paired.prev")
#About the same


###Full PH model - cannot add Year to this model - is singular with year
cph.wet.full <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev + Paired.prev, data = pag2ss.op.wet,model="ph") #,bs_samples=1000)
summary(cph.wet.full)
cph.op.wet.coef = coef(cph.wet.full)[1]
cph.pp.wet.coef = coef(cph.wet.full)[2]

##Plot
plot.data = data.frame(Ornamented.prev=c("Ornamented","Brown"),Year="2019",Paired.prev="Yes")
rownames(plot.data) = c("Ornamented","Brown")
plot(cph.wet.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")

plot.data = data.frame(Paired.prev=c("Yes","No"),Year="2019",Ornamented.prev="Ornamented")
rownames(plot.data) = c("Paired","Not paired")
plot(cph.wet.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")


##Permutation test for whether a male was ornamented previously
cph.op.wet.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.op.wet$Ornamented.prev.rand = sample(pag2ss.op.wet$Ornamented.prev,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev.rand + Paired.prev, data = pag2ss.op.wet,model="ph"))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#A = cph.op.wet.rand(perm.number = 1000,coef.obs = cph.op.wet.coef)

##Paper: In wet years, males that were ornamented as 1-year-olds molted earlier as 2-year-olds than males that were not ornamented as 1-year-olds
#(coef=1.26,p=0.014). 



##Permutation test for whether a male was paired previously
cph.pp.wet.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.op.wet$Paired.prev.rand = sample(pag2ss.op.wet$Paired.prev,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev + Paired.prev.rand, data = pag2ss.op.wet,model="ph"))[2]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#B = cph.pp.wet.rand(perm.number = 1000,coef.obs = cph.pp.wet.coef)

##Paper: In wet years, males that were paired previously were not more likely to molt into ornamented plumage earlier than males that 
#were not paired previously (coef=0.62,p=0.151). 

##Plot both:
ggarrange(A,B,labels=c("A","B"))




###Plot above results for wet years:

##I didn't put year in the model so doesn't make sense to plot yearly trends in this section

##Ornamented previously for wet years
#Get ornamented data from model
plot.orn = expand.grid(Ornamented.prev="Ornamented",Paired.prev="Yes",Day=130:365) 
est.orn = data.frame(getFitEsts(cph.wet.full,newdata = plot.orn,q=plot.orn$Day),Day=plot.orn$Day)
colnames(est.orn)[1] = "est"

#Get brown data from model
plot.brown = expand.grid(Ornamented.prev="Brown",Paired.prev="Yes",Day=130:365) 
est.brown =  data.frame(getFitEsts(cph.wet.full,newdata = plot.brown,q=plot.brown$Day),Day=plot.brown$Day)
colnames(est.brown)[1] = "est"

#Combine
est.orn$Ornamented.prev="Ornamented"
est.brown$Ornamented.prev="Brown"
est.orn.brown = rbind(est.orn,est.brown)

#Plot
ggplot() + geom_line(data=est.brown,aes(x=Day,y=est,color=Ornamented.prev),size=1.5) +
  geom_line(data=est.orn,aes(x=Day,y=est,color=Ornamented.prev),size=1.5) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_manual(values=c("#C6A375","#D82F2B"),name="Plumage in\\nprevious year") +
  ylab("Likelihood of molt") + xlab("Day of year")




###Paired previously for wet years
#Get paired data from model
plot.pairy = expand.grid(Ornamented.prev="Brown",Paired.prev="Yes",Day=130:365) 
est.pairy = data.frame(getFitEsts(cph.wet.full,newdata = plot.pairy,q=plot.pairy$Day),Day=plot.pairy$Day)
colnames(est.pairy)[1] = "est"

#Get brown data from model
plot.pairn = expand.grid(Ornamented.prev="Brown",Paired.prev="No",Day=130:365) 
est.pairn =  data.frame(getFitEsts(cph.wet.full,newdata = plot.pairn,q=plot.pairn$Day),Day=plot.pairn$Day)
colnames(est.pairn)[1] = "est"

#Combine
est.pairy$Paired.prev="Yes"
est.pairn$Paired.prev="No"
est.pairy.pairn = rbind(est.pairy,est.pairn)

#Plot
ggplot() + geom_line(data=est.pairn,aes(x=Day,y=est,color=Paired.prev),size=1.5) +
  geom_line(data=est.pairy,aes(x=Day,y=est,color=Paired.prev),size=1.5) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_manual(values=c("#C6A375","#D82F2B"),name="Paired in\\nprevious year") +
  ylab("Likelihood of molt") + xlab("Day of year")








####Look at dry years separately####

pag2ss.op.dry = pag2ss.op %>% filter(Year %in% c("2015","2017"))

##Run both proportional hazards and proportional odds models without the Paired.prev and Year variables to first check model diagnostics
cph.ph.dry1 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.op.dry,model="ph") #,bs_samples=1000)
cph.po.dry1 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.op.dry,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph.dry1) #Same for both models
diag_covar(cph.ph.dry1,varName = "Ornamented.prev")
diag_covar(cph.po.dry1,varName = "Ornamented.prev")
#About the same

##Run both proportional hazards and proportional odds models without the Ornamented.prev and Year variables to first check model diagnostics
cph.ph.dry2 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev, data = pag2ss.op.dry,model="ph") #,bs_samples=1000)
cph.po.dry2 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev, data = pag2ss.op.dry,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph.dry2) #Same for both models
diag_covar(cph.ph.dry2,varName = "Paired.prev")
diag_covar(cph.po.dry2,varName = "Paired.prev")
#About the same


###Full PH model - cannot add Year to this model - is singular with year
cph.dry.full <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev + Paired.prev, data = pag2ss.op.dry,model="ph") #,bs_samples=1000)
summary(cph.dry.full)
cph.op.dry.coef = coef(cph.dry.full)[1]
cph.pp.dry.coef = coef(cph.dry.full)[2]

##Plot
plot.data = data.frame(Ornamented.prev=c("Ornamented","Brown"),Year="2019",Paired.prev="Yes")
rownames(plot.data) = c("Ornamented","Brown")
plot(cph.dry.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")

plot.data = data.frame(Paired.prev=c("Yes","No"),Year="2019",Ornamented.prev="Ornamented")
rownames(plot.data) = c("Paired","Not paired")
plot(cph.dry.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")


##Permutation test for whether a male was ornamented previously
cph.op.dry.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.op.dry$Ornamented.prev.rand = sample(pag2ss.op.dry$Ornamented.prev,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev.rand + Paired.prev, data = pag2ss.op.dry,model="ph"))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#A = cph.op.dry.rand(perm.number = 1000,coef.obs = cph.op.dry.coef)

##Paper: In dry years, no difference in molt dates in relation to previous ornmentation (coef=0.21,p=0.332).


##Permutation test for whether a male was paired previously
cph.pp.dry.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.op.dry$Paired.prev.rand = sample(pag2ss.op.dry$Paired.prev,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev + Paired.prev.rand, data = pag2ss.op.dry,model="ph"))[2]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#B = cph.pp.dry.rand(perm.number = 1000,coef.obs = cph.pp.dry.coef)

##Paper: In dry years, no difference in molt dates among males that were paired as 1-year-olds and those that were not (coef=0.32,p=0.304)

##Plot both
ggarrange(A,B,labels=c("A","B"))



###Plot above results for dry years:

##I didn't put year in the model so doesn't make sense to plot yearly trends in this section

##Ornamented previously for dry years
#Get ornamented data from model
plot.orn = expand.grid(Ornamented.prev="Ornamented",Paired.prev="Yes",Day=130:365) 
est.orn = data.frame(getFitEsts(cph.dry.full,newdata = plot.orn,q=plot.orn$Day),Day=plot.orn$Day)
colnames(est.orn)[1] = "est"

#Get brown data from model
plot.brown = expand.grid(Ornamented.prev="Brown",Paired.prev="Yes",Day=130:365) 
est.brown =  data.frame(getFitEsts(cph.dry.full,newdata = plot.brown,q=plot.brown$Day),Day=plot.brown$Day)
colnames(est.brown)[1] = "est"

#Combine
est.orn$Ornamented.prev="Ornamented"
est.brown$Ornamented.prev="Brown"
est.orn.brown = rbind(est.orn,est.brown)

#Plot
ggplot() + geom_line(data=est.brown,aes(x=Day,y=est,color=Ornamented.prev),size=1.5) +
  geom_line(data=est.orn,aes(x=Day,y=est,color=Ornamented.prev),size=1.5) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_manual(values=c("#C6A375","#D82F2B"),name="Plumage in\\nprevious year") +
  ylab("Likelihood of molt") + xlab("Day of year")




###Paired previously for dry years
#Get paired data from model
plot.pairy = expand.grid(Ornamented.prev="Brown",Paired.prev="Yes",Day=130:365) 
est.pairy = data.frame(getFitEsts(cph.dry.full,newdata = plot.pairy,q=plot.pairy$Day),Day=plot.pairy$Day)
colnames(est.pairy)[1] = "est"

#Get brown data from model
plot.pairn = expand.grid(Ornamented.prev="Brown",Paired.prev="No",Day=130:365) 
est.pairn =  data.frame(getFitEsts(cph.dry.full,newdata = plot.pairn,q=plot.pairn$Day),Day=plot.pairn$Day)
colnames(est.pairn)[1] = "est"

#Combine
est.pairy$Paired.prev="Yes"
est.pairn$Paired.prev="No"
est.pairy.pairn = rbind(est.pairy,est.pairn)

#Plot
ggplot() + geom_line(data=est.pairn,aes(x=Day,y=est,color=Paired.prev),size=1.5) +
  geom_line(data=est.pairy,aes(x=Day,y=est,color=Paired.prev),size=1.5) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_manual(values=c("#C6A375","#D82F2B"),name="Paired in\\nprevious year") +
  ylab("Likelihood of molt") + xlab("Day of year")









####Within males that were paired previously, does staying paired to the same female influence molt date?####

##Get males that were paired previously
pag2ss.sf = pag2ss.op %>% filter(Paired.prev=="Yes")

##Run both proportional hazards and proportional odds models without the Paired.prev and Year variables to first check model diagnostics
cph.ph3 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.sf,model="ph") #,bs_samples=1000)
cph.po3 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.sf,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph3) #Same for both models
diag_covar(cph.ph3,varName = "Ornamented.prev")
diag_covar(cph.po3,varName = "Ornamented.prev")
#The ornamented previously covariates cross a lot here, so much that should maybe not be included in the model? 
#Removing it from the model did not change anything

##Run both proportional hazards and proportional odds models without the Ornamented.prev and Year variables to first check model diagnostics
cph.ph4 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF, data = pag2ss.sf,model="ph") #,bs_samples=1000)
cph.po4 <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF, data = pag2ss.sf,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph4) #Same for both models
diag_covar(cph.ph4,varName = "Paired.prev.sameF")
diag_covar(cph.po4,varName = "Paired.prev.sameF")
#About the same but PO slightly larger difference. Use ph - coxph - to keep consistent. 

###Full PH model 
cph.sf.full <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF + Ornamented.prev + Year, data = pag2ss.sf,model="ph") #,bs_samples=1000)
summary(cph.sf.full)
cph.sf.coef = coef(cph.sf.full)[1]

##Plot
plot.data = data.frame(Paired.prev.sameF=c("Yes","No"),Year="2019",Ornamented.prev="Ornamented")
rownames(plot.data) = c("Same Female","New Female")
plot(cph.sf.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")


##Permutation test for whether a male was paired to the same female previously
cph.pp.sf.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.sf$Paired.prev.sameF.rand = sample(pag2ss.sf$Paired.prev.sameF,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF.rand + Ornamented.prev + Year, data = pag2ss.sf,model="ph"))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
}

##Run permutation
#cph.pp.sf.rand(perm.number = 1000,coef.obs = cph.sf.coef)

##Paper: Within males that paired previously, no significant difference in timing of molt between males that were paired to the same 
#female and those that were not (coef=0.41,p=0.167).




###Plot all years paired to the same female previously
#Set plumage = brown, because all ornamented males were paired, so does not make sense to have a not-paired ornamented male in prediction

#Get sameF yes data from model
plot.pairsfy = expand.grid(Paired.prev.sameF="Yes",Year=c(2015:2019),Day=130:365,Ornamented.prev="Brown") 
plot.pairsfy$Year = as.factor(plot.pairsfy$Year)
est.pairsfy = data.frame(getFitEsts(cph.sf.full,newdata = plot.pairsfy,q=plot.pairsfy$Day),Year=plot.pairsfy$Year,Day=plot.pairsfy$Day)
colnames(est.pairsfy)[1] = "est"

#Get sameF no data from model
plot.pairsfn = expand.grid(Paired.prev.sameF="No",Year=c(2015:2019),Day=130:365,Ornamented.prev="Brown") 
plot.pairsfn$Year = as.factor(plot.pairsfn$Year)
est.pairsfn =  data.frame(getFitEsts(cph.sf.full,newdata = plot.pairsfn,q=plot.pairsfn$Day),Year=plot.pairsfn$Year,Day=plot.pairsfn$Day)
colnames(est.pairsfn)[1] = "est"

#Combine
est.pairsfy$Paired.prev.sameF="Yes"
est.pairsfn$Paired.prev.sameF="No"
est.pairsfy.pairsfn = rbind(est.pairsfy,est.pairsfn)

#Plot
ggplot() + geom_line(data=est.pairsfy,aes(x=Day,y=est,group=Year,color=Year),size=1) + 
  geom_line(data=est.pairsfn,aes(x=Day,y=est,group=Year,color=Year),size=1) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_viridis(discrete = T)

ggplot(data=est.pairsfy.pairsfn,aes(x=Day,y=est)) + geom_line(aes(group=Paired.prev.sameF,color=Paired.prev.sameF),size=1.5) + 
  facet_wrap(facets=vars(Year)) + scale_color_manual(values=c("gray55","darkorange2"),name="Paired to same\\nfemale in previous\\nyear") + theme_cowplot() +
  ylab("Likelihood of molt") + xlab("Day of year")







####Within males that were paired previously, does staying paired to the same female influence molt date? - wet years####

pag2ss.sf.wet = pag2ss.sf %>% filter(Year %in% c("2016","2018","2019"))


##Run both proportional hazards and proportional odds models without the Paired.prev and Year variables to first check model diagnostics
cph.ph3.wet <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.sf.wet,model="ph") #,bs_samples=1000)
cph.po3.wet <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.sf.wet,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph3.wet) #Same for both models
diag_covar(cph.ph3.wet,varName = "Ornamented.prev")
diag_covar(cph.po3.wet,varName = "Ornamented.prev")
#The ornamented previously covariates cross a lot here, so much that should maybe not be included in the model. 

##Run both proportional hazards and proportional odds models without the Ornamented.prev and Year variables to first check model diagnostics
cph.ph4.wet <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF, data = pag2ss.sf.wet,model="ph") #,bs_samples=1000)
cph.po4.wet <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF, data = pag2ss.sf.wet,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph4.wet) #Same for both models
diag_covar(cph.ph4.wet,varName = "Paired.prev.sameF")
diag_covar(cph.po4.wet,varName = "Paired.prev.sameF")
#About the same 

###Full PH model 
cph.sf.wet.full <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF + Ornamented.prev, data = pag2ss.sf.wet,model="ph") #,bs_samples=1000)
summary(cph.sf.wet.full)
cph.sf.wet.coef = coef(cph.sf.wet.full)[1]

##Plot
plot.data = data.frame(Paired.prev.sameF=c("Yes","No"),Ornamented.prev="Ornamented")
rownames(plot.data) = c("Same Female","New Female")
plot(cph.sf.wet.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")


##Permutation test for whether a male was paired to the same female previously
cph.pp.sf.wet.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.sf.wet$Paired.prev.sameF.rand = sample(pag2ss.sf.wet$Paired.prev.sameF,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF.rand + Ornamented.prev, data = pag2ss.sf.wet,model="ph"))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
}

##Run permutation
#cph.pp.sf.wet.rand(perm.number = 1000,coef.obs = cph.sf.wet.coef)

##Paper: Within wet years, within males that paired previously, no significant difference in timing of molt between males that were paired to the same 
#female and those that were not (coef=0.50,p=0.157).




###Plot all years paired to the same female previously
#Set plumage = brown, because all ornamented males were paired, so does not make sense to have a not-paired ornamented male in prediction

#Get sameF yes data from model
plot.pairsfy = expand.grid(Paired.prev.sameF="Yes",Day=130:365,Ornamented.prev="Brown") 
est.pairsfy = data.frame(getFitEsts(cph.sf.wet.full,newdata = plot.pairsfy,q=plot.pairsfy$Day),Day=plot.pairsfy$Day)
colnames(est.pairsfy)[1] = "est"

#Get sameF no data from model
plot.pairsfn = expand.grid(Paired.prev.sameF="No",Day=130:365,Ornamented.prev="Brown") 
est.pairsfn =  data.frame(getFitEsts(cph.sf.wet.full,newdata = plot.pairsfn,q=plot.pairsfn$Day),Day=plot.pairsfn$Day)
colnames(est.pairsfn)[1] = "est"

#Combine
est.pairsfy$Paired.prev.sameF="Yes"
est.pairsfn$Paired.prev.sameF="No"
est.pairsfy.pairsfn = rbind(est.pairsfy,est.pairsfn)

#Plot
ggplot() +  geom_line(data=est.pairsfn,aes(x=Day,y=est,color=Paired.prev.sameF),size=1.5) +
  geom_line(data=est.pairsfy,aes(x=Day,y=est,color=Paired.prev.sameF),size=1.5) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_manual(values=c("gray55","darkorange2"),name="Paired to same\\nfemale in previous\\nyear") + 
  ylab("Likelihood of molt") + xlab("Day of year")






####Within males that were paired previously, does staying paired to the same female influence molt date? - dry years####

pag2ss.sf.dry = pag2ss.sf %>% filter(Year %in% c("2015","2017"))


##Run both proportional hazards and proportional odds models without the Paired.prev and Year variables to first check model diagnostics
cph.ph3.dry <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.sf.dry,model="ph") #,bs_samples=1000)
cph.po3.dry <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Ornamented.prev, data = pag2ss.sf.dry,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph3.dry) #Same for both models
diag_covar(cph.ph3.dry,varName = "Ornamented.prev")
diag_covar(cph.po3.dry,varName = "Ornamented.prev")

##Run both proportional hazards and proportional odds models without the Ornamented.prev and Year variables to first check model diagnostics
cph.ph4.dry <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF, data = pag2ss.sf.dry,model="ph") #,bs_samples=1000)
cph.po4.dry <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF, data = pag2ss.sf.dry,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph4.dry) #Same for both models
diag_covar(cph.ph4.dry,varName = "Paired.prev.sameF")
diag_covar(cph.po4.dry,varName = "Paired.prev.sameF")
#About the same 

###Full PH model 
cph.sf.dry.full <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF + Ornamented.prev, data = pag2ss.sf.dry,model="ph") #,bs_samples=1000)
summary(cph.sf.dry.full)
cph.sf.dry.coef = coef(cph.sf.dry.full)[1]

##Plot
plot.data = data.frame(Paired.prev.sameF=c("Yes","No"),Ornamented.prev="Ornamented")
rownames(plot.data) = c("Same Female","New Female")
plot(cph.sf.dry.full,plot.data,xlab="Day of year",ylab="Likelihood of molting",fun="cdf")


##Permutation test for whether a male was paired to the same female previously
cph.pp.sf.dry.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pag2ss.sf.dry$Paired.prev.sameF.rand = sample(pag2ss.sf.dry$Paired.prev.sameF,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Paired.prev.sameF.rand + Ornamented.prev, data = pag2ss.sf.dry,model="ph"))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
}

##Run permutation
#cph.pp.sf.dry.rand(perm.number = 1000,coef.obs = cph.sf.dry.coef)

##Paper: Within dry years, within males that paired previously, no significant difference in timing of molt between males that were paired to the same 
#female and those that were not (coef=0.58,p=0.164).



###Plot all years paired to the same female previously
#Set plumage = brown, because all ornamented males were paired, so does not make sense to have a not-paired ornamented male in prediction

#Get sameF yes data from model
plot.pairsfy = expand.grid(Paired.prev.sameF="Yes",Day=130:365,Ornamented.prev="Brown") 
est.pairsfy = data.frame(getFitEsts(cph.sf.dry.full,newdata = plot.pairsfy,q=plot.pairsfy$Day),Day=plot.pairsfy$Day)
colnames(est.pairsfy)[1] = "est"

#Get sameF no data from model
plot.pairsfn = expand.grid(Paired.prev.sameF="No",Day=130:365,Ornamented.prev="Brown") 
est.pairsfn =  data.frame(getFitEsts(cph.sf.dry.full,newdata = plot.pairsfn,q=plot.pairsfn$Day),Day=plot.pairsfn$Day)
colnames(est.pairsfn)[1] = "est"

#Combine
est.pairsfy$Paired.prev.sameF="Yes"
est.pairsfn$Paired.prev.sameF="No"
est.pairsfy.pairsfn = rbind(est.pairsfy,est.pairsfn)

#Plot
ggplot() +  geom_line(data=est.pairsfn,aes(x=Day,y=est,color=Paired.prev.sameF),size=1.5) +
  geom_line(data=est.pairsfy,aes(x=Day,y=est,color=Paired.prev.sameF),size=1.5) +
  ylim(0,1) + xlim(0,365) + theme_cowplot() + scale_color_manual(values=c("gray55","darkorange2"),name="Paired to same\\nfemale in previous\\nyear") + 
  ylab("Likelihood of molt") + xlab("Day of year")



