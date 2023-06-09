####When do old males molt script 

##Looking at males that are 3 years old and older


####Load Data####

##Load plumage data from database
library(here)
plumages = read.csv(here::here("Input files","All RBFW ind histories plumages.csv"))

##Load in age data from database
ages = read.csv(here::here("Input files","All RBFW ind histories for ages.csv"))

##Load in breeding group data from database
groups = read.csv(here::here("Input files","All RBFW Groups.csv"))
groups$Date.Created = as.Date(groups$Date.Created,"%m/%d/%y")


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

##Load in capture data with plumage scores
captures.p = read.csv(here::here("Input files","All RBFW Captures plumage scores.csv"))
captures.p$Date = as.Date(captures.p$Date,"%m/%d/%y")



####Get plumage and age files into long format so I can merge####
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(lme4)
library(DHARMa)
library(lubridate)
library(viridis)


#Melt plumages twice
plumages.years = plumages %>% select(contains("Year")) #get year data
plumages.plums = plumages %>% select(1:4,contains("Plumage")) #get plumage data
plumages.years.m = melt(plumages.years,id.vars = NULL) #melt both
plumages.plums.m = melt(plumages.plums,id.vars = 1:4) #melt both
plumages.plums.m = plumages.plums.m %>% select(-variable) 
plumages.py = data.frame(plumages.plums.m,plumages.years.m$value) #combine
colnames(plumages.py)[5:6] = c("Plumage","Year") #rename columns
plumages.py = plumages.py[!is.na(plumages.py$Year),] #remove rows with no data

#Order plumages.py by plumage and take top plumage for each year
plumages.py$Plumage = as.factor(plumages.py$Plumage)
plumages.py = plumages.py[order(factor(plumages.py$Plumage,levels=c("Bright","Interm","Dull",""))),]
plumages.pyd = plumages.py %>% distinct(Fwnumber,Year,.keep_all = T)

#Melt ages twice
ages.years = ages %>% select(contains("Year")) #get year data
ages.ages = ages %>% select(1:4,contains("Age")) #get age data
ages.years.m = melt(ages.years,id.vars = NULL) #melt both
ages.ages.m = melt(ages.ages,id.vars = 1:4) #melt both
ages.ages.m = ages.ages.m %>% select(-variable) 
ages.py = data.frame(ages.ages.m,ages.years.m$value) #combine
colnames(ages.py)[5:6] = c("Age","Year") #rename columns
ages.py = ages.py[!is.na(ages.py$Year),] #remove rows with no data
ages.py = ages.py %>% select(FWno,Age,Age.exact,Year)

#Get a single age for each year
ages.pyd = ages.py %>% distinct(FWno,Year,.keep_all = T)

#Merge files
p.a = merge(plumages.pyd,ages.pyd,by.x=c("Fwnumber","Year"),by.y =c("FWno","Year"))





####Set up 3+ male dataframe####

##Select 3+ males 
p.a = p.a %>% filter(Sex=="M",Age>=3)

##Match up colors to bird using agesex files - this will be that individuals colors in the year of interest pag2$Year and will match to ind file
agesex15.s = agesex15 %>% select(Year,FWNo,Bird)
agesex16.s = agesex16 %>% select(Year,FWNo,Bird)
agesex17.s = agesex17 %>% select(Year,FWNo,Bird)
agesex18.s = agesex18 %>% select(Year,FWNo,Bird)
agesex19.s = agesex19 %>% select(Year,FWNo,Bird)
agesexall = rbind(agesex15.s,agesex16.s,agesex17.s,agesex18.s,agesex19.s)
pa3 = merge(p.a,agesexall,by.x=c("Fwnumber","Year"),by.y=c("FWNo","Year"))


##Make sure no birds changed color combos within a year
captures.p.test = captures.p %>% filter(Fwnumber %in% pa3$Fwnumber)
captures.p.test = captures.p.test %>% distinct(Year,Fwnumber,Colors)
check = captures.p.test %>% group_by(Year,Fwnumber) %>% filter(n()>1)
#So in captures, need to change those birds to the combo they are in agesex files for that year 7776 = RGG
check1 = captures.p %>% filter(Year==2017 & Colors=="RGX") 
check1$Colors = "RGG"
captures.p = rbind(check1,captures.p)
captures.p = captures.p %>% distinct(Fwnumber,Date,.keep_all = T) %>% filter(Fwnumber %in% pa3$Fwnumber) %>% arrange(Date) #Remove later values 

##Make sure capture colors are the same as colors in agesex files for that year
capture.colorcheck = pa3 %>% select(Fwnumber, Year, Bird)
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


##Get birds of interest - birds in pa3
pa3.15 = pa3 %>% filter(Year=="2015")
ind15.s = ind15.s %>% filter(Bird %in% pa3.15$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2015")
pa3.16 = pa3 %>% filter(Year=="2016")
ind16.s = ind16.s %>% filter(Bird %in% pa3.16$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2016")
pa3.17 = pa3 %>% filter(Year=="2017")
ind17.s = ind17.s %>% filter(Bird %in% pa3.17$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2017")
pa3.18 = pa3 %>% filter(Year=="2018")
ind18.s = ind18.s %>% filter(Bird %in% pa3.18$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2018")
pa3.19 = pa3 %>% filter(Year=="2019")
ind19.s = ind19.s %>% filter(Bird %in% pa3.19$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2019")

##Combine, order, and get julian date
ind.all = rbind(ind15.s,ind16.s,ind17.s,ind18.s,ind19.s)
ind.all = ind.all %>% arrange(Year,Bird,Date)
ind.all$jdate = yday(ind.all$Date)

##If date is in January, make Julian date +365
for (i in 1:nrow(ind.all)) {
  if(month(ind.all$Date[i])==1) {ind.all$jdate[i] = ind.all$jdate[i] + 365}
}




####Look at timing of molt across years####

##Get a bird.year factor for grouping
ind.all$Bird.year = paste(ind.all$Bird,ind.all$Year,sep="")

##Plot molt by date
ggplot(data=ind.all,aes(x=jdate,y=Molt)) + geom_point(aes(color=Year),size=1) + geom_line(aes(group=Bird.year,color=Year),size=0.7) + 
  theme_cowplot() + scale_color_viridis(discrete = T)


##Get week data - add one if in January - one week and two week intervals
ind.all$week = week(ind.all$Date)
for (i in 1:nrow(ind.all)) {
  if(ind.all$week[i]==1) {ind.all$week[i]=54}}
ind.all$twoweek = ceiling(ind.all$week/2) *2 #two week intervals - divide week by 2 and round up, then multiply by two to get later week of 2-week interval

##For these plots - use sighting data only, capture data muddy the averages a little bit since they don't match sightings exactly
ind.all.nocap = ind.all %>% filter(method=="sighting")

##Get one sighting per week for each individual - best pluamge score of that week
duplicated(pa3$Bird) #Below code works - distinct for bird and year because no duplicate color combos in this dataset
ind.all.d = ind.all.nocap %>% arrange(desc(Molt)) %>% distinct(Bird,Year,twoweek,.keep_all = T) #Get one sighting of each bird per two week interval

##Group by year and week and get average plumage score and get count of individuals seen in that week
ind.all.week = ind.all.d %>% group_by(Year,twoweek) %>% summarise(Molt=mean(Molt))
ind.all.weekn = ind.all.d %>% group_by(Year,twoweek) %>% count(twoweek)
ind.all.week$count = ind.all.weekn$n

##Require at least 5 individuals present to keep the week in 
ind.all.week = ind.all.week %>% filter(count>=5)

##Plot year means
ggplot(data=ind.all.week,aes(x=twoweek,y=Molt)) + geom_line(aes(color=Year),size=1.8) + 
  theme_cowplot() + scale_color_viridis(discrete = T) + xlab("Week of year") + ylab("Mean old male (3+) plumage score") + 
  geom_point(size=2,aes(color=Year)) + ylim(0,100)
#In 2017 3yo males were much later in molting than any other years 





####Look at timing of molt within males####

###Measure difference in time between last dull observation and first intermediate/first bright
##Filter pa3 for only birds with molt data
pa3s = pa3 %>% filter(Bird %in% ind.all$Bird) 

##For loop to get difference in dates
pa3s$first.seen = NA
pa3s$times.seen = NA
pa3s$last.dull = NA
pa3s$first.int = NA
pa3s$last.int = NA
pa3s$first.bright = NA
pa3s$first.plum.afterdull = NA
pa3s$first.bright.plum = NA
pa3s$first.plum = NA
for (i in 1:nrow(pa3s)) {
  bird = pa3s$Bird[i]
  year = pa3s$Year[i]
  ind.bird = ind.all %>% filter(Bird==as.character(bird),Year==year) %>% arrange(jdate)
  pa3s$first.plum[i] = ind.bird$Molt[1]
  pa3s$first.seen[i] = ind.bird$jdate[1] 
  pa3s$times.seen[i] = nrow(ind.bird)
  ind.dull = ind.bird %>% filter(Molt<33) %>% arrange(desc(jdate))
  pa3s$last.dull[i] = ind.dull$jdate[1]
  ind.int = ind.bird %>% filter(Molt>=33,Molt<=66) %>% arrange(jdate)
  pa3s$first.int[i] = ind.int$jdate[1]
  ind.int2 = ind.int %>% arrange(desc(jdate))
  pa3s$last.int[i] = ind.int2$jdate[1]
  ind.bright= ind.bird %>% filter(Molt>66) %>% arrange(jdate)
  pa3s$first.bright[i] = ind.bright$jdate[1]
  pa3s$first.bright.plum[i] = ind.bright$Molt[1]
  ind.firstplumAD = ind.bird %>% filter(Molt>=33) %>% arrange(jdate)
  pa3s$first.plum.afterdull[i] = ind.firstplumAD$Molt[1]
}

##Get difference in dates
pa3s$dull.to.bright = pa3s$first.bright-pa3s$last.dull
pa3s$dull.to.int = pa3s$first.int-pa3s$last.dull
pa3s$int.to.bright = pa3s$first.bright-pa3s$last.int
pa3s$fint.to.bright = pa3s$first.bright-pa3s$first.int
pa3s$dull.to.intorbright = NA
for (i in 1:nrow(pa3s)) {if(!is.na(pa3s$dull.to.int[i])) {
  pa3s$dull.to.intorbright[i] = pa3s$dull.to.int[i]} else {pa3s$dull.to.intorbright[i] = pa3s$dull.to.bright[i]}
}



##Select birds with an intermediate plumage or early bright plumage less than 90 - we would have seen these males during molt and molt is usually very fast
pa3s.fpd = pa3s %>% filter(first.plum.afterdull<90)

##For males that went from dull to very bright quickly, restrict to individuals seen less than 20 days between last dull and first bright
#Look
pa3s.dvb = pa3s %>% filter(!first.plum.afterdull<90) %>% filter(!is.na(dull.to.bright))
ggplot(pa3s.dvb,aes(x=dull.to.bright)) + geom_histogram(color="black",fill="gray",binwidth = 5) + theme_cowplot()
ind.all.test2 = ind.all %>% filter(Bird %in% pa3s.dvb$Bird)
ggplot(data=ind.all.test2,aes(x=jdate,y=Molt))  + geom_line(size=1,aes(group=Bird.year)) +  geom_point(aes(color=Year),size=1.2) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T)  + 
  geom_hline(aes(yintercept=66),lty=2,color="gray") + xlab("Day of year") + ylab("Plumage score") + ylim(0,100) +
  geom_hline(aes(yintercept=33),lty=2,color="gray")
#Restrict
pa3s.dvb = pa3s.dvb %>% filter(dull.to.bright<=20)
ind.all.test3 = ind.all %>% filter(Bird %in% pa3s.dvb$Bird)
ggplot(data=ind.all.test3,aes(x=jdate,y=Molt))  + geom_line(size=1,aes(group=Bird.year)) +  geom_point(aes(color=Year),size=1.2) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) + xlim(135,300) + 
  geom_hline(aes(yintercept=66),lty=2,color="gray") + xlab("Day of year") + ylab("Plumage score") + ylim(0,100) +
  geom_hline(aes(yintercept=33),lty=2,color="gray")
#Looks better

##For males that were first seen as bright - 90 or above
#Look
pa3s.fsb = pa3s %>% filter(first.plum>=90)
ind.all.test4 = merge(ind.all,pa3s.fsb,by=c("Bird","Year"))
ggplot(data=ind.all.test4,aes(x=jdate,y=Molt))  + geom_line(size=1,aes(group=Bird.year)) +  geom_point(aes(color=Year),size=1.2) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) + xlim(135,300) + 
  geom_hline(aes(yintercept=66),lty=2,color="gray") + xlab("Day of year") + ylab("Plumage score") + ylim(0,100) +
  geom_hline(aes(yintercept=33),lty=2,color="gray")
ggplot(data=pa3s,aes(x=first.seen,y=Year)) + geom_line(aes(group=Year)) + geom_point() + 
  geom_point(data=pa3s.fsb,aes(x=first.seen,y=Year),color="red")
#Some males seen really late. Use only males seen before day 200 - will match with 2yo analysis
pa3s.fsb = pa3s.fsb %>% filter(first.seen<=200)



###Combine
pa3ss = rbind(pa3s.fpd,pa3s.dvb,pa3s.fsb)
pa3ss = pa3ss %>% distinct(Fwnumber,Year,.keep_all = T) #Remove duplicates - none here but good check
pa3ss$Bird.year = paste(pa3ss$Bird,pa3ss$Year,sep="")
ind.alls = ind.all %>% filter(Bird.year %in% pa3ss$Bird.year)



###Plot trajectories
##Plot for individual males - jdate
ggplot(data=ind.alls,aes(x=jdate,y=Molt))  + geom_line(size=1,aes(group=Bird.year)) +  geom_point(aes(color=Year),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")

##Plot for individual males - jdate but color points by method
# ggplot(data=ind.alls,aes(x=jdate,y=Molt))  + geom_line(size=1) +  geom_point(aes(color=method),size=1) +
#   theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_manual(values=c("red","black")) +
#   geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2)



##Look at how many 3yo males we saw during molt and across how many years
pa3ss.sdm = pa3ss %>% filter(first.plum<90)
table(pa3ss.sdm$Year)
length(unique(pa3ss.sdm$Fwnumber))

##Plot
ind.alls.sdm = ind.alls %>% filter(Bird.year %in% pa3ss.sdm$Bird.year)
ggplot(data=ind.alls.sdm,aes(x=jdate,y=Molt))  + geom_line(size=1,aes(group=Bird.year)) +  geom_point(aes(color=Year),size=1) +
  theme_cowplot() + facet_wrap(facets=vars(as.character(Bird))) + scale_color_viridis(discrete = T) +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  xlab("Day of year") + ylab("Plumage score")

##That's 37 males seen during molt, one seen in two years. Can run same sliding window climat analyses for 3yos and see if same 
#rainfall period is responsible for predicting when they're likely to molt. Would be interesting to know if they're 
#responding to the same rainfall period as 2yos just maybe able to respond to it more quickly and molt faster? 


##Get molt dates
pa3ss.sdm$molt.date = NA
for (i in 1:nrow(pa3ss.sdm)) {
  if(is.na(pa3ss.sdm$first.int[i])) {pa3ss.sdm$molt.date[i]=pa3ss.sdm$first.bright[i]} else {
    if(pa3ss.sdm$first.int[i]<=pa3ss.sdm$first.bright[i] | is.na(pa3ss.sdm$first.bright[i])) {pa3ss.sdm$molt.date[i]=pa3ss.sdm$first.int[i]} else {
      pa3ss.sdm$molt.date[i] = pa3ss.sdm$first.bright[i]}
  }
}





####Look at whether paired previously and ornamented previously influence molt date####




##Determine if what their plumage was in the previous year
pa3ss.sdm$Plumage.prev = NA
for (i in 1:nrow(pa3ss.sdm)) {
  pag.male = plumages.pyd %>% filter(Fwnumber==pa3ss.sdm$Fwnumber[i]) %>% filter(Year==(as.numeric(as.character(pa3ss.sdm$Year[i]))-1))
  if(nrow(pag.male)==0) {pa3ss.sdm$Plumage.prev[i]=NA} else {
    pa3ss.sdm$Plumage.prev[i] = as.character(pag.male$Plumage)
  }
}
pa3ss.sdm$Ornamented.prev = NA
for (i in 1:nrow(pa3ss.sdm)) {
  if(pa3ss.sdm$Plumage.prev[i]=="") {pa3ss.sdm$Ornamented.prev[i]=NA} else {
    if(pa3ss.sdm$Plumage.prev[i]=="Dull") {pa3ss.sdm$Ornamented.prev[i]="Brown"} else {pa3ss.sdm$Ornamented.prev[i]="Ornamented"}}
}
pa3ss.sdm$Ornamented.prev = factor(pa3ss.sdm$Ornamented.prev,levels=c("Brown","Ornamented")) 


##Determine if they were paired in the previous year
pa3ss.sdm$Paired.prev = NA
for (i in 1:nrow(pa3ss.sdm)) {
  pag.male = groups %>% filter(Male.fwno==pa3ss.sdm$Fwnumber[i]) %>% filter(Year==(as.numeric(as.character(pa3ss.sdm$Year[i]))-1))
  if(nrow(pag.male)==0) {pa3ss.sdm$Paired.prev[i]="No"} else {
    pa3ss.sdm$Paired.prev[i]="Yes"
  }
}


##In the previous year, only two males were not paired, and only one male was dull, with one other male not seen in the previous year


##Write three-year-old dataframe to use for climate analyses
#write.csv(pa3ss.sdm,here::here("Output files","threeyo_moltdates_sdm_dataframe.csv"),row.names = F)

##Write three-year old dataframe with all males - even those first seen as bright for climate and molt date plots
#write.csv(pa3ss,here::here("Output files","threeyo_moltdates_dataframe.csv"),row.names = F)

##Write three-year-olds in 2017 dataframe for social connections analyses
pa3ss.sdm17 = pa3ss.sdm %>% filter(Year=="2017")
#write.csv(pa3ss.sdm17,here::here("Output files","3yo_male_list_2017.csv"),row.names=F)






####See if age influences when 3yo+ males molt####

##Get molt dates
pa3ss$molt.date = NA
for (i in 1:nrow(pa3ss)) {
  if(is.na(pa3ss$first.int[i])) {pa3ss$molt.date[i]=pa3ss$first.bright[i]} else {
    if(pa3ss$first.int[i]<=pa3ss$first.bright[i] | is.na(pa3ss$first.bright[i])) {pa3ss$molt.date[i]=pa3ss$first.int[i]} else {
      pa3ss$molt.date[i] = pa3ss$first.bright[i]}
  }
}


##Survival packages
library(survival)
library(survminer)
library(icenReg)

##Set up columns for survival analysis - if left censored (first seen as bright so don't know true molt date), then time1 = infinity (NA) and
#time2 = date first seen bright. If we know the true date (or have an approximation), then time1 and time 2 are the same
#Survival package will not run coxph on interval-censored data. icenReg will - make left-censored data time1 = 0 according to paper - Anderson-Bergman
pa3ss$molt.time1 = NA
pa3ss$molt.time2 = NA
for (i in 1:nrow(pa3ss)) {
  if(pa3ss$first.plum[i] >= 90) {
    pa3ss$molt.time1[i]=0
    pa3ss$molt.time2[i]=pa3ss$molt.date[i]} else {
      pa3ss$molt.time1[i]=pa3ss$molt.date[i]
      pa3ss$molt.time2[i]=pa3ss$molt.date[i]
    }
}
pa3ss$Year = as.factor(pa3ss$Year)



####Age and molt date for 3+ males

##Select only males with exact age
pa3ss.age = pa3ss %>% filter(Age.exact == "exact")

##Run both proportional hazards and proportional odds models without the Paired.prev and Year variables to first check model diagnostics
cph.ph <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Age, data = pa3ss.age,model="ph") #,bs_samples=1000)
cph.po <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Age, data = pa3ss.age,model="po") #,bs_samples=1000)

##Check diagnostics and covariates - both assume that Kaplan Meier should not cross for a covariate, whatever model shows greater difference is
#typically better
diag_baseline(cph.ph) #Same for both models
diag_covar(cph.ph,varName = "Age")
diag_covar(cph.po,varName = "Age")
#About the same


###Full PH model 
cph.full <- ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Age + Year, data = pa3ss.age,model="ph") #,bs_samples=1000)
summary(cph.full)
cph.age.coef = coef(cph.full)[1]

##Plot
plot.data = data.frame(Age=c(3,4,5,6,7),Year="2015")
rownames(plot.data) = c(3,4,5,6,7)
plot(cph.full,plot.data,xlab="Day of year",ylab="Likelihood of acquiring ornamented plumage",fun="cdf")


##Permutation test for Age
cph.age.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    pa3ss.age$Age.rand = sample(pa3ss.age$Age,replace=F)
    rands[i] = coef(ic_sp(Surv(molt.time1,molt.time2,type="interval2") ~ Age.rand + Year, data = pa3ss.age,model="ph"))[1]
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
#cph.age.rand(perm.number = 1000,coef.obs = cph.age.coef)

##Paper: within 3+ males, age is still an important factor for determining molt date, with older males molting earlier (coef=0.54,p=0.006)


