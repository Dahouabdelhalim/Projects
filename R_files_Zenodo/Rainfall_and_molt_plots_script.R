####Rainfall and molt plots script

##Using this script to plot rainfall, NDVI, and timing of molt across age classes

library(here)
library(dplyr)
library(reshape2)
library(lubridate)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)

####Load data####

##Load rainfall data
rain = read.csv(here::here("Input files","rain_2004-2018.csv"))
rain$Date = as.Date(rain$Date,"%m/%d/%y")

##Load NDVI data
ndvi = read.csv(here::here("Input files","Landsat7_all_years.csv"))
ndvi$Date = as.Date(ndvi$Date,"%m/%d/%y")


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

##Load molt date dataframes - these are not all birds in these age classes, but are birds that we had good molt dates for or first saw as bright
m2 = read.csv(here::here("Output files","twoyo_moltdates_all_dataframe.csv"))
m3 = read.csv(here::here("Output files","threeyo_moltdates_dataframe.csv"))

##Load whether males molted or not dataframe
pag = read.csv(here::here("Output files","all_males_moltornot.csv"))

##Load in list of implant males
implantms = read.csv(here::here("Input files","Implant males.csv"))




####Calculate molt dates for 1yo males####

#Haven't done this yet since very few 1yo males made it to bright

##Get 1yo males in each year
agesex15.1 = agesex15 %>% filter(Sex=="M",Current.Age==1,Age.exact.=="exact")
agesex16.1 = agesex16 %>% filter(Sex=="M",Current.Age==1,Age.Exact.=="exact")
agesex17.1 = agesex17 %>% filter(Sex=="M",Current.Age==1,Ageexact.=="exact")
agesex18.1 = agesex18 %>% filter(Sex=="M",Current.Age==1,Ageexact.=="exact")
agesex19.1 = agesex19 %>% filter(Sex=="M",Current.Age==1,Ageexact.=="exact")

##Get individual list for those males
ind15.1 = ind15 %>% filter(Bird %in% agesex15.1$Bird)
ind16.1 = ind16 %>% filter(Bird %in% agesex16.1$Bird)
ind17.1 = ind17 %>% filter(Bird %in% agesex17.1$Bird)
ind18.1 = ind18 %>% filter(Bird %in% agesex18.1$Bird)
ind19.1 = ind19 %>% filter(Bird %in% agesex19.1$Bird)

##Get males that made it to ornamented plumage
ind15.1b = ind15.1 %>% filter(Molt>=33) %>% arrange(Date) %>% distinct(Bird,.keep_all = T)
ind16.1b = ind16.1 %>% filter(Molt>=33) %>% arrange(Date) %>% distinct(Bird,.keep_all = T)
ind17.1b = ind17.1 %>% filter(Molt>=33) %>% arrange(Date) %>% distinct(Bird,.keep_all = T)
ind18.1b = ind18.1 %>% filter(Molt>=33) %>% arrange(Date) %>% distinct(Bird,.keep_all = T)
ind19.1b = ind19.1 %>% filter(Molt>=33) %>% arrange(Date) %>% distinct(Bird,.keep_all = T)

##Get julian dates
ind15.1b$jdate = yday(ind15.1b$Date)
ind16.1b$jdate = yday(ind16.1b$Date)
ind17.1b$jdate = yday(ind17.1b$Date)
ind18.1b$jdate = yday(ind18.1b$Date)
ind19.1b$jdate = yday(ind19.1b$Date)

##Get year and select important columns
ind15.1b = ind15.1b %>% mutate(Year="2015") %>% select(Year,Bird,Date,jdate,Molt)
ind16.1b = ind16.1b %>% mutate(Year="2016") %>% select(Year,Bird,Date,jdate,Molt)
ind17.1b = ind17.1b %>% mutate(Year="2017") %>% select(Year,Bird,Date,jdate,Molt)
ind18.1b = ind18.1b %>% mutate(Year="2018") %>% select(Year,Bird,Date,jdate,Molt)
ind19.1b = ind19.1b %>% mutate(Year="2019") %>% select(Year,Bird,Date,jdate,Molt)

##Combine
ind1 = rbind(ind15.1b,ind16.1b,ind17.1b,ind18.1b,ind19.1b)


###Remove T-implant males and males from removal communities in 2017 season from dataset

##Remove T-implant males - no 1yo males made it to ornamented before 10/1 in removal communities 
implantms.t = implantms %>% filter(Implant=="T")
ind1 = ind1 %>% filter(!Bird %in% implantms.t$Colors)


##Plot molt by day 
ggplot(ind1,aes(x=jdate,y=Molt)) + geom_point(aes(color=Year),size=2) + theme_cowplot() + xlim(0,365) + ylim(0,100)
#Lots of 1yo males were not seen in oranmented plumage until late. Know that many of these birds were experiencing adventitious molt
#Last molt date of 2yo male that was not adventitious was Oct 1, jdate 274
ggplot(ind1,aes(x=jdate,y=Molt)) + geom_point(aes(color=Year),size=2) + theme_cowplot() + xlim(0,365) + ylim(0,100) + 
  geom_vline(aes(xintercept=274))
#Consider molt after 274 to probably be adventitious. For one year olds, some of these first ornamented dates might not be 
#extremely accurate since we don't have nearly as many sightings from the breeding season. But it should show the natural 
#trends well enough for these graphs - there are a decent number of 1yo males that obtain nuptial plumage late in the year
#probably through adventitious molt. 



####Quality check 1yo molt dates####

##Just looking at indYEAR.1 files - males that only have a month or so between dull and ornamented sightings leaving in. 
#Males that have greater than that, taking out. Many 1yo males are getting to ornamented plumage in late months. 

##2015:
remove15 = data.frame(birds=c("YIY"))
##2016:
#Looks good
##2017: 
remove17 = data.frame(birds=c("BLV"))
##2018:
remove18 = data.frame(birds=c("LHG"))
##2019:
remove19 = data.frame(birds=c("IZR"))

##Remove those males from ind1 list
remove1s = rbind(remove15,remove17,remove18,remove19)
ind1 = ind1 %>% filter(!Bird %in% remove1s$birds)




####Set up 3yo dataframe####

##Get molt dates
m3$molt.date = NA
for (i in 1:nrow(m3)) {
  if(is.na(m3$first.int[i])) {m3$molt.date[i]=m3$first.bright[i]} else {
    if(m3$first.int[i]<=m3$first.bright[i] | is.na(m3$first.bright[i])) {m3$molt.date[i]=m3$first.int[i]} else {
      m3$molt.date[i] = m3$first.bright[i]}
  }
}

##Get molt.time1 to match m2
m3$molt.time1 = NA
for (i in 1:nrow(m3)) {if(m3$first.plum[i]>=90) {m3$molt.time1[i]=0} else {m3$molt.time1[i]=m3$molt.date[i]}
}

##Going to start plot at 182, July 1st, so if any 3yos were first seen as bright after July 1st, remove them
#No 2yo males with this problem
m3 = m3 %>% filter(!(molt.time1==0 & molt.date>182))




####Get total males in each year####

##Get males that did not molt in each year
pag.d = pag %>% filter(Plumage=="Dull") %>% filter(Year %in% c("2015","2016","2017","2018","2019"))

##Summarize by years
pag.ds = pag.d %>% group_by(Year,Age) %>% summarise(Brown = n())

##Summarize number ornamented by year for each age class
ind1s = ind1 %>% group_by(Year) %>% summarise(Ornamented = n())
m2s = m2 %>% group_by(Year) %>% summarise(Ornamented = n())
m3s = m3 %>% group_by(Year) %>% summarise(Ornamented = n())

##Get total birds for each year 1yos
pag.ds1 = pag.ds %>% filter(Age==1) 
pag.ds1 = merge(pag.ds1,ind1s,by="Year")
pag.ds1$Total = pag.ds1$Brown + pag.ds1$Ornamented

##Get total birds for each year 2yos
pag.ds2 = pag.ds %>% filter(Age==2)
pag.ds2 = merge(pag.ds2,m2s,by="Year",all.y=T)
pag.ds2[is.na(pag.ds2$Brown),]$Brown <- 0
pag.ds2$Total = pag.ds2$Brown + pag.ds2$Ornamented





####Get percentage bright for each date####

##Start plot at July 1st - otherwise going to run into the problem of having some males seen for the first time after the date
#in ornamented plumage. By starting at July 1st I already cut off the few males in that category. So now % bright on July 1st and following
#will be accurate. 

##Three-year-olds
pb3yos = expand.grid(Age=3,jdate=c(182:365),percentB=NA,Year=c(2015:2019))

for (i in 1:nrow(pb3yos)) {
  #For each day, filter m3 for that year
  m3.year = m3 %>% filter(Year==pb3yos$Year[i])
  #For each day, get birds molted by that day
  m3.day = m3.year %>% filter(molt.date<=pb3yos$jdate[i])
  #Add in percentage molted by that day 
  pb3yos$percentB[i] = nrow(m3.day)/nrow(m3.year)
}


##Two-year-olds
pb2yos = expand.grid(Age=2,jdate=c(182:365),percentB=NA,Year=c(2015:2019))

for (i in 1:nrow(pb2yos)) {
  #For each row, get the total
  m2.total = (pag.ds2 %>% filter(Year==pb2yos$Year[i]))$Total
  #For each day, filter m2 for that year
  m2.year = m2 %>% filter(Year==pb2yos$Year[i])
  #For each day, get birds molted by that day
  m2.day = m2.year %>% filter(molt.date<=pb2yos$jdate[i])
  #Add in percentage molted by that day 
  pb2yos$percentB[i] = nrow(m2.day)/m2.total
}


##One-year-olds
pb1yos = expand.grid(Age=1,jdate=c(182:365),percentB=NA,Year=c(2015:2019))

for (i in 1:nrow(pb1yos)) {
  #For each row, get the total
  m1.total = (pag.ds1 %>% filter(Year==pb1yos$Year[i]))$Total
  #For each day, filter m1 for that year
  m1.year = ind1 %>% filter(Year==pb1yos$Year[i])
  #For each day, get birds molted by that day
  m1.day = m1.year %>% filter(jdate<=pb1yos$jdate[i])
  #Add in percentage molted by that day 
  pb1yos$percentB[i] = nrow(m1.day)/m1.total
}

##Get earlier dates for each age class
pb3yose = expand.grid(Age=3,jdate=c(1:181),percentB=0,Year=c(2015:2019))
pb2yose = expand.grid(Age=2,jdate=c(1:181),percentB=0,Year=c(2015:2019))
pb1yose = expand.grid(Age=1,jdate=c(1:181),percentB=0,Year=c(2015:2019))

##Combine all ages
pball = rbind(pb1yos,pb2yos,pb3yos,pb3yose,pb2yose,pb1yose)
pball$Age.Year = paste(pball$Age,pball$Year,sep="")
ggplot(pball,aes(x=jdate,y=percentB)) + geom_point() + geom_line(aes(group=Age.Year)) + facet_wrap(facets="Year")

##Get weeks
pball$date = as.Date("2015-01-01")+pball$jdate -1
pball$week = week(pball$date)

##Drop week 53 - no %bright numbers change with week 53 and week 53 is really part of week 1 of next year
pball = pball %>% filter(week!=53)

##Summarize by week - take max percentB value from each week
pball.week = pball %>% group_by(Year,week,Age) %>% summarise(percentB = max(percentB))
pball.week$Age.Year = paste(pball.week$Age,pball.week$Year)
pball.week$Age = factor(pball.week$Age,levels=c("3","2","1"))
pball.week$percentB = pball.week$percentB * 100

##Plot by week with facet plot
ggplot(pball.week,aes(x=week,y=percentB)) + geom_line(aes(group=Age.Year)) + geom_point(size=2,aes(fill=Age),color="black",pch=21) +
  facet_wrap(facets="Year",nrow=1) + theme_cowplot() + scale_fill_manual(values=c("black","gray","white")) + xlim(0,52) + 
  ylab("Percent\\nornamented") + xlab("Week of year")

##Set up weeks for plotting continuously
pball.week$week2 = NA
for (i in 1:nrow(pball.week)) {
  if(pball.week$Year[i]==2015) {pball.week$week2[i] = pball.week$week[i]}
  if(pball.week$Year[i]==2016) {pball.week$week2[i] = pball.week$week[i] + 52}
  if(pball.week$Year[i]==2017) {pball.week$week2[i] = pball.week$week[i] + 104}
  if(pball.week$Year[i]==2018) {pball.week$week2[i] = pball.week$week[i] + 156}
  if(pball.week$Year[i]==2019) {pball.week$week2[i] = pball.week$week[i] + 208}
}

##Plot by week without a good x-axis scale
ggplot(pball.week,aes(x=week2,y=percentB)) + geom_line(aes(group=Age.Year)) + geom_point(size=2,aes(fill=Age),color="black",pch=21) +
  theme_cowplot() + scale_fill_manual(values=c("black","gray","white")) + 
  ylab("Percent\\nornamented") + xlab("Week of year") + scale_x_continuous(breaks=seq(0,260,by=10))

##Get scale for weeks - this gets 10,20, etc week of each year and the corresponding week2 continuous variable
week.scale = pball.week %>% filter(week %in% c(10,20,30,40,50)) %>% distinct(week2)
week0 = data.frame(week2=0,Year=2015,week=0)
week.scale2 = rbind(week0,data.frame(week.scale))

##Plot by week - good x-axis scale
ggplot(pball.week,aes(x=week2,y=percentB)) + geom_line(aes(group=Age.Year)) + geom_point(size=2,aes(fill=Age),color="black",pch=21) +
  theme_cowplot() + scale_fill_manual(values=c("black","gray","white")) + 
  ylab("Percent\\nornamented") + xlab("Week of year") + scale_x_continuous(breaks=week.scale2$week2,labels=week.scale2$week)






####Rainfall####


##Restrict rain to 2015-2019 seasons
rain = rain %>% filter(Year %in% c("2015","2016","2017","2018","2019"))

##Get cumulative rainfall
rain$Year.cuml = NA
rain$Year.cuml[1] = rain$amountcm[1]
for (i in 2:nrow(rain)) {
  if(rain$Year[i]==rain$Year[i-1]) {rain$Year.cuml[i] = rain$Year.cuml[i-1] + rain$amountcm[i]} else {
    rain$Year.cuml[i] = rain$amountcm[i]}
}

##Get rainfall by week
rain$week = week(rain$Date)

##Add week 53 dates to week 52
for (i in 1:nrow(rain)) {
  if(rain$week[i]==53) {rain$week[i]=52}
}

##Summarize rainfall by week
rain.week = rain %>% group_by(Year,week) %>% summarise(rain.cuml = max(Year.cuml),rain.total=sum(amountcm))

##Set up weeks for plotting continuously
rain.week$week2 = NA
for (i in 1:nrow(rain.week)) {
  if(rain.week$Year[i]==2015) {rain.week$week2[i] = rain.week$week[i]}
  if(rain.week$Year[i]==2016) {rain.week$week2[i] = rain.week$week[i] + 52}
  if(rain.week$Year[i]==2017) {rain.week$week2[i] = rain.week$week[i] + 104}
  if(rain.week$Year[i]==2018) {rain.week$week2[i] = rain.week$week[i] + 156}
  if(rain.week$Year[i]==2019) {rain.week$week2[i] = rain.week$week[i] + 208}
}

##Separate rainfall into yearly dataframes to plot cumulative lines separately
rain.week.15 = rain.week %>% filter(Year=="2015")
rain.week.16 = rain.week %>% filter(Year=="2016")
rain.week.17 = rain.week %>% filter(Year=="2017")
rain.week.18 = rain.week %>% filter(Year=="2018")
rain.week.19 = rain.week %>% filter(Year=="2019")

##Plot rainfall with facet plot
ggplot(rain.week,aes(x=week,y=rain.cuml)) + geom_line(size=1,color="dark blue") + 
  facet_wrap(facets="Year",nrow=1) + theme_cowplot() +  xlim(0,52) + 
  ylab("Rainfall in cm") + xlab("Week of year") + geom_line(aes(x=week,y=rain.total)) 
ggplot(rain.week,aes(x=week,y=rain.total)) + geom_line(size=0.75) + 
  facet_wrap(facets="Year",nrow=1) + theme_cowplot() +  xlim(0,52) + ylim(0,30) +
  ylab("") + xlab("Week of year") 

##Plot rainfall continuous
ggplot(rain.week,aes(x=week2,y=rain.total)) + geom_line(size=0.75) + 
  theme_cowplot() + ylim(0,30) + ylab("") + xlab("Week of year") + 
  scale_x_continuous(breaks=week.scale2$week2,labels=week.scale2$week)

ggplot(rain.week,aes(x=week2,y=rain.total)) + geom_line(size=0.75) + 
  theme_cowplot() + ylab("") + xlab("Week of year") + 
  scale_x_continuous(breaks=week.scale2$week2,labels=week.scale2$week) +
  geom_line(data=rain.week,aes(x=week2,y=rain.cuml),size=1,color="dark blue")

##Plot rainfall but cumulative lines separately
ggplot(rain.week,aes(x=week2,y=rain.total)) + geom_line(size=0.75) + 
  theme_cowplot() + ylab("") + xlab("Week of year") + 
  scale_x_continuous(breaks=week.scale2$week2,labels=week.scale2$week) +
  geom_line(data=rain.week.15,aes(x=week2,y=rain.cuml),size=1,color="dark blue") + 
  geom_line(data=rain.week.16,aes(x=week2,y=rain.cuml),size=1,color="dark blue") + 
  geom_line(data=rain.week.17,aes(x=week2,y=rain.cuml),size=1,color="dark blue") + 
  geom_line(data=rain.week.18,aes(x=week2,y=rain.cuml),size=1,color="dark blue") + 
  geom_line(data=rain.week.19,aes(x=week2,y=rain.cuml),size=1,color="dark blue")
  

####Plot both
A = ggplot(pball.week,aes(x=week2,y=percentB)) + geom_line(aes(group=Age.Year)) + geom_point(size=2,aes(fill=Age),color="black",pch=21) +
  theme_cowplot() + scale_fill_manual(values=c("black","gray","white")) + theme(legend.position = "none") +
  ylab("Percent\\nornamented") + xlab("Week of year") + scale_x_continuous(breaks=week.scale2$week2,labels=week.scale2$week)
B = ggplot(rain.week,aes(x=week2,y=rain.total)) + geom_line(size=0.75) + 
  theme_cowplot() + ylim(0,30) + ylab("Rainfall (cm)") + xlab("Week of year") + 
  scale_x_continuous(breaks=week.scale2$week2,labels=week.scale2$week)
ggarrange(A,B,ncol=1,align="v")




####NDVI####

##Get julian date for NDVI
ndvi$jdate = yday(ndvi$Date)

##Get season (Year for NDVI)
ndvi$Year = year(ndvi$Date) + 1

##Restrict NDVI to 2015-2019 seasons
ndvi = ndvi %>% filter(Year %in% c("2015","2016","2017","2018","2019"))

##Get NDVI by week
ndvi$week = week(ndvi$Date)

##Add week 53 dates to week 52
for (i in 1:nrow(ndvi)) {
  if(ndvi$week[i]==53) {ndvi$week[i]=52}
}

##Summarize ndvi by week
ndvi.week = ndvi %>% group_by(Year,week) %>% summarise(NDVI = mean(NDVI))

##Set up weeks for plotting continuously
ndvi.week$week2 = NA
for (i in 1:nrow(ndvi.week)) {
  if(ndvi.week$Year[i]==2015) {ndvi.week$week2[i] = ndvi.week$week[i]}
  if(ndvi.week$Year[i]==2016) {ndvi.week$week2[i] = ndvi.week$week[i] + 52}
  if(ndvi.week$Year[i]==2017) {ndvi.week$week2[i] = ndvi.week$week[i] + 104}
  if(ndvi.week$Year[i]==2018) {ndvi.week$week2[i] = ndvi.week$week[i] + 156}
  if(ndvi.week$Year[i]==2019) {ndvi.week$week2[i] = ndvi.week$week[i] + 208}
}

##Plot NDVI with facet wrap
ggplot(ndvi.week,aes(x=week,y=NDVI)) + geom_area(fill="light green",alpha=0.8) +
  geom_line(size=0.7) + xlim(0,52) + ylim(0,1) + xlab("Week of year") + ylab("NDVI") +
  facet_wrap(facets="Year",nrow=1) + theme_cowplot()

##Plot ndvi continuous
ggplot(ndvi.week,aes(x=week2,y=NDVI)) + geom_line(size=0.75) + 
  theme_cowplot() + ylim(0,1) + ylab("") + xlab("Week of year") + 
  scale_x_continuous(breaks=week.scale2$week2,labels=week.scale2$week)

##Match NDVI scale to rainfall scale 
ndvi.week$rain.scale = ndvi.week$NDVI * 30

##Plot NDVI with rain scale
ggplot(ndvi.week,aes(x=week2,y=rain.scale)) + geom_area(fill="light green",alpha=0.8) +
  geom_line(size=0.7) + ylim(0,30) + xlab("Week of year") + ylab("NDVI") + theme_cowplot()

##Plot rain and NDVI together
ggplot(rain.week,aes(x=week2,y=rain.total)) + 
  geom_area(data=ndvi.week,aes(x=week2,y=rain.scale),fill="#9bca72",alpha=0.8) +
  geom_line(data=ndvi.week,aes(x=week2,y=rain.scale),size=0.7,color="#9bca72") +
  geom_line(size=0.7,color="#91bae8") + 
  geom_area(fill="#91bae8",alpha=0.8) +
  theme_cowplot() +  ylim(0,30) +
  ylab("") + xlab("Week of year") + ylab("Rainfall in cm\\n") + 
  scale_x_continuous(breaks=week.scale2$week2,labels=week.scale2$week)





####Plot plumage, NDVI, and rainfall####

##Make a dataframe to show season extents
seasons = data.frame(season=c(2015,2016,2017,2018,2019),start=c(22,74,126,178,230),end=c(55,107,159,211,263))
seasons = melt(seasons,id.vars = 1)

##Make a dataframe to show approximate non-breeding and breeding periods
#Breeding mid-August (week 34 = Aug 20) to mid-January (week 3 = Jan 15)
breed = data.frame(breed=c("breed1","breed2","breed3","breed4","breed5","breed6"),start=c(1,34,86,138,190,242),end=c(3,55,107,159,211,263))
breed = melt(breed,id.vars = 1)
nbreed = data.frame(nb=c("nb1","nb2","nb3","nb4","nb5"),start=c(4,56,108,160,212),end=c(33,85,137,189,241))
nbreed = melt(nbreed,id.vars = 1)

##Add in blanks to week scale
week.scale3 = data.frame(week2 = seq(0,263,by=2),Year="",week="")
week.scale3 = rbind(week.scale2,week.scale3)
week.scale3 = week.scale3 %>% distinct(week2,.keep_all = T) %>% arrange(week2)

##Plot them all - these don't align perfectly on the vertical axis
#Rectangles here block off misleading areas - when the area fields do not drop off exactly or the 3yos start at zero and race up 
#up to 100. More accurate to have the 3yo line start at where they were on July 1st - usually above 80%
A = ggplot(pball.week,aes(x=week2,y=percentB)) + 
  geom_area(aes(fill=Age,group=Age),position="identity",alpha=0.8) +  geom_line(aes(group=Age.Year,color=Age),size=1) + 
  theme_cowplot() + scale_color_manual(values=c("#2f6680","#94473b","#cc7032")) + 
  scale_fill_manual(values=c("#2f6680","#94473b","#cc7032")) +
  ylab("Percent\\nornamented") + xlab("Week of year") + theme(legend.position = "none") + 
  scale_x_continuous(breaks=week.scale3$week2,labels=week.scale3$week) +
  geom_rect(mapping=aes(xmin=0, xmax=26.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=52, xmax=78.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=104, xmax=130.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=156, xmax=182.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=208, xmax=234.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=260, xmax=261, ymin=-2, ymax=101),fill="white") +
  geom_line(data=seasons,aes(x=value,y=116,group=season),size=1) +
  geom_line(data=breed,aes(x=value,y=108,group=breed),size=1.5,color="dark gray") +
  geom_line(data=nbreed,aes(x=value,y=108,group=nb),size=0.5,color="dark gray",lty=2) +
  scale_y_continuous(breaks=c(0,25,50,75,100))

B = ggplot(rain.week,aes(x=week2,y=rain.total)) + 
  geom_area(data=ndvi.week,aes(x=week2,y=rain.scale),fill="#87c757",alpha=0.7) +
  geom_line(data=ndvi.week,aes(x=week2,y=rain.scale),size=1,color="#87c757") +
  geom_line(size=0.7,color="#0272b1") + 
  geom_area(fill="#0272b1",alpha=0.8) +
  theme_cowplot() + ylim(0,30) +
  ylab("") + xlab("Week of year") + ylab("Rainfall (cm)\\nand NDVI") +
  scale_x_continuous(breaks=week.scale3$week2,labels=week.scale3$week)

##Plot together - don't use this figure - USE VERSION BELOW! 
ggarrange(A,B,ncol=1,align="v")

##Going to make my own key in AD for this figure. 




###This is absurd, but to get plot x axes to line up exactly, they have to have all parts of each other. 
##So plot with all parts, using alpha to hide unwanted parts. This actually works to line up the x axes. 

A = ggplot(pball.week,aes(x=week2,y=percentB)) + 
  geom_area(aes(fill=Age,group=Age),position="identity",alpha=0.8) +  geom_line(aes(group=Age.Year,color=Age),size=1) + 
  theme_cowplot() + scale_color_manual(values=c("#2f6680","#94473b","#cc7032")) + 
  scale_fill_manual(values=c("#2f6680","#94473b","#cc7032")) +
  ylab("Percent\\nornamented") + xlab("Week of year") + theme(legend.position = "none") + 
  scale_x_continuous(breaks=week.scale3$week2,labels=week.scale3$week) +
  geom_rect(mapping=aes(xmin=0, xmax=26.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=52, xmax=78.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=104, xmax=130.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=156, xmax=182.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=208, xmax=234.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=260, xmax=261, ymin=-2, ymax=101),fill="white") +
  geom_line(data=seasons,aes(x=value,y=116,group=season),size=1) +
  geom_line(data=breed,aes(x=value,y=108,group=breed),size=1.5,color="dark gray") +
  geom_line(data=nbreed,aes(x=value,y=108,group=nb),size=0.5,color="dark gray",lty=2) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) + 
  geom_area(data=ndvi.week,aes(x=week2,y=rain.scale),fill="#87c757",alpha=0) +
  geom_line(data=ndvi.week,aes(x=week2,y=rain.scale),size=1,color="#87c757",alpha=0) +
  geom_line(data=rain.week,aes(x=week2,y=rain.total),size=0.7,color="#0272b1",alpha=0) + 
  geom_area(data=rain.week,aes(x=week2,y=rain.total),fill="#0272b1",alpha=0)

B = ggplot(pball.week,aes(x=week2,y=percentB)) + 
  geom_area(aes(fill=Age,group=Age),position="identity",alpha=0) +  
  geom_line(aes(group=Age.Year,color=Age),size=1,alpha=0) + 
  theme_cowplot() + scale_color_manual(values=c("#2f6680","#94473b","#cc7032")) + 
  scale_fill_manual(values=c("#2f6680","#94473b","#cc7032")) +
  xlab("Week of year") + theme(legend.position = "none") + 
  scale_x_continuous(breaks=week.scale3$week2,labels=week.scale3$week) +
  geom_rect(mapping=aes(xmin=0, xmax=26.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=52, xmax=78.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=104, xmax=130.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=156, xmax=182.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=208, xmax=234.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=260, xmax=261, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_line(data=seasons,aes(x=value,y=116,group=season),size=1,alpha=0) +
  geom_line(data=breed,aes(x=value,y=108,group=breed),size=1.5,color="dark gray",alpha=0) +
  geom_line(data=nbreed,aes(x=value,y=108,group=nb),size=0.5,color="dark gray",lty=2,alpha=0) +
  geom_area(data=ndvi.week,aes(x=week2,y=rain.scale),fill="#87c757",alpha=0.7) +
  geom_line(data=ndvi.week,aes(x=week2,y=rain.scale),size=1,color="#87c757") +
  geom_line(data=rain.week,aes(x=week2,y=rain.total),size=0.7,color="#0272b1") + 
  geom_area(data=rain.week,aes(x=week2,y=rain.total),fill="#0272b1",alpha=0.8) +
  ylim(0,30) + ylab("Rainfall (cm)\\nand NDVI")

ggarrange(A,B,ncol=1,align="v")



##Get NDVI y-axis scale - do same plotting of everything just in case
A = ggplot(pball.week,aes(x=week2,y=percentB)) + 
  geom_area(aes(fill=Age,group=Age),position="identity",alpha=0.8) +  geom_line(aes(group=Age.Year,color=Age),size=1) + 
  theme_cowplot() + scale_color_manual(values=c("#2f6680","#94473b","#cc7032")) + 
  scale_fill_manual(values=c("#2f6680","#94473b","#cc7032")) +
  ylab("Percent\\nornamented") + xlab("Week of year") + theme(legend.position = "none") + 
  scale_x_continuous(breaks=week.scale3$week2,labels=week.scale3$week) +
  geom_rect(mapping=aes(xmin=0, xmax=26.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=52, xmax=78.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=104, xmax=130.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=156, xmax=182.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=208, xmax=234.5, ymin=-2, ymax=101),fill="white") +
  geom_rect(mapping=aes(xmin=260, xmax=261, ymin=-2, ymax=101),fill="white") +
  geom_line(data=seasons,aes(x=value,y=116,group=season),size=1) +
  geom_line(data=breed,aes(x=value,y=108,group=breed),size=1.5,color="dark gray") +
  geom_line(data=nbreed,aes(x=value,y=108,group=nb),size=0.5,color="dark gray",lty=2) +
  scale_y_continuous(breaks=c(0,25,50,75,100)) + 
  geom_area(data=ndvi.week,aes(x=week2,y=rain.scale),fill="#87c757",alpha=0) +
  geom_line(data=ndvi.week,aes(x=week2,y=rain.scale),size=1,color="#87c757",alpha=0) +
  geom_line(data=rain.week,aes(x=week2,y=rain.total),size=0.7,color="#0272b1",alpha=0) + 
  geom_area(data=rain.week,aes(x=week2,y=rain.total),fill="#0272b1",alpha=0)

B = ggplot(pball.week,aes(x=week2,y=percentB)) + 
  geom_area(aes(fill=Age,group=Age),position="identity",alpha=0) +  
  geom_line(aes(group=Age.Year,color=Age),size=1,alpha=0) + 
  theme_cowplot() + scale_color_manual(values=c("#2f6680","#94473b","#cc7032")) + 
  scale_fill_manual(values=c("#2f6680","#94473b","#cc7032")) +
  xlab("Week of year") + theme(legend.position = "none") + 
  scale_x_continuous(breaks=week.scale3$week2,labels=week.scale3$week) +
  geom_rect(mapping=aes(xmin=0, xmax=26.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=52, xmax=78.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=104, xmax=130.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=156, xmax=182.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=208, xmax=234.5, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_rect(mapping=aes(xmin=260, xmax=261, ymin=-2, ymax=101),fill="white",alpha=0) +
  geom_line(data=seasons,aes(x=value,y=116,group=season),size=1,alpha=0) +
  geom_line(data=breed,aes(x=value,y=108,group=breed),size=1.5,color="dark gray",alpha=0) +
  geom_line(data=nbreed,aes(x=value,y=108,group=nb),size=0.5,color="dark gray",lty=2,alpha=0) +
  geom_area(data=ndvi.week,aes(x=week2,y=NDVI),fill="#87c757",alpha=0.7) +
  geom_line(data=ndvi.week,aes(x=week2,y=NDVI),size=1,color="#87c757") +
  geom_line(data=rain.week,aes(x=week2,y=rain.total),size=0.7,color="#0272b1",alpha=0) + 
  geom_area(data=rain.week,aes(x=week2,y=rain.total),fill="#0272b1",alpha=0) +
  ylim(0,1) + ylab("Rainfall (cm)\\nand NDVI")

ggarrange(A,B,ncol=1,align="v")



