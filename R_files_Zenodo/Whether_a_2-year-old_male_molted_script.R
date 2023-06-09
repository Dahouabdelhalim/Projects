###Whether a two-year-old male molted script

#Use this script for analyses looking at why some 2yo males molt and others do not



####Load Data####

##Load plumage data from database
library(here)
plumages = read.csv(here::here("Input files","All RBFW ind histories plumages.csv"))

##Load in age data from database
ages = read.csv(here::here("Input files","All RBFW ind histories for ages.csv"))

##Load in breeding group data from database
groups = read.csv(here::here("Input files","All RBFW Groups.csv"))
groups$Date.Created = as.Date(groups$Date.Created,"%m/%d/%y")

##Load in list of implant males
implantms = read.csv(here::here("Input files","Implant males.csv"))

##Load list of birds in removal communities in 2017 season and their FWnumbers
removals = read.csv(here::here("Input files","RemovalorControl 2017.csv"))
ages17 = read.csv(here::here("Input files","Ages and Sexes 2017.csv"))





####Get plumage and age files into long format so I can merge####
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(lme4)
library(DHARMa)


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


####Remove T-implant males and males from removal communities in 2017 season from dataset

##Remove T-implant males
implantms.t = implantms %>% filter(Implant=="T")
p.a$Drop = NA
for (i in 1:nrow(p.a)) {
  if(p.a$Fwnumber[i] %in% implantms.t$Fwnumber & p.a$Year[i] %in% implantms.t$Year) {p.a$Drop[i]="drop"} else {p.a$Drop[i]="no"}
}

##Remove Removal males
removals = merge(removals,ages17,by="Bird",all.x=T) #Merge community ID with ages to get FWnumber. This ages file has colorbands specific to what they were in 2017 season
removals = removals %>% filter(TreatmentGroup=="Removal") %>% select(TreatmentGroup,FWNo)
#Removals happened in 2017 season
for (i in 1:nrow(p.a)) {
  if(p.a$Fwnumber[i] %in% removals$FWNo & p.a$Year[i]==2017) {p.a$Drop[i]="drop"}
}

##Remove both
p.a = p.a %>% filter(Drop!="drop") %>% select(-Drop)






####Set up two year old dataframe####

###Select males and remove individuals without a plumage score
#A few individuals would disappear for a year and reappear and we don't know what plumage they were that year
#Also many nestlings don't have plumage scores entered - get birds one year old and older
#Select only individuals with a known age (age = exact)
p.a = p.a %>% filter(Sex=="M",Plumage!="",Age>=1,Age.exact=="exact")



###Make sure all individuals have a breeding group for that season - shows that we followed them into the breeding season
#thereby giving them time to molt if they were going to. A male that was only seen during the non-breeding season may not 
#have had time to molt. Then order by status with male.fwno first because want to know if they were ever paired that season.
#In some instances a male paired, then lost his female, then was a helper. Also order by date.created to get first male.fwno status
#of each season so later can look up his female. 
groups = groups %>% select(Year,Group,Date.Created,contains("fwno")) #Get columns that have Fwnumbers
groups.m = melt(groups,id.vars = 1:3) #Melt
groups.m = groups.m %>% filter(value!="") #Remove blanks
colnames(groups.m)[4:5] = c("Status","Fwnumber") #Rename columns
#Set male.fwno as top level in status factor so can order by that - want to know if male was ever a breeder in a year
groups.m$Status = factor(groups.m$Status,levels=c("Male.fwno","Female.fwno","Female.2.fwno","Helper1.fwno",
                                                  "Helper2.fwno","Helper3.fwno","Helper4.fwno","Helper5.fwno","Other1.fwno",
                                                  "Other2.fwno","Other3.fwno"))
groups.m = groups.m %>% arrange(Status,Date.Created) #Order by status and date created - status first
#Remove duplicates based on ID and year, getting highest breeding status group for each individual in each year. 
groups.md = groups.m %>% distinct(Year,Fwnumber,.keep_all = T) 

#Merge groups.md and p.a to make sure all males have a breeding group
p.a.g = merge(p.a, groups.md,by=c("Fwnumber","Year"))



####Get two-year-olds
pag2 = p.a.g %>% filter(Age==2)
pag2$Year = as.factor(pag2$Year)


##Assign whether a male was Ornamented or not (binary variable)
pag2$Ornamented = NA
for (i in 1:nrow(pag2)) {
  if(pag2$Plumage[i]=="Dull") {pag2$Ornamented[i]="Brown"} else {pag2$Ornamented[i]="Ornamented"}
}
pag2$Ornamented = factor(pag2$Ornamented,levels=c("Brown","Ornamented")) 


##Determine if two year olds were ever paired in the season they were that age
pag2$Paired = NA
for (i in 1:nrow(pag2)) {
  if(pag2$Status[i]=="Male.fwno") {pag2$Paired[i]="Yes"} else {pag2$Paired[i]="No"}
}
pag2$Paired = as.factor(pag2$Paired)


##Determine if what their plumage was in the previous year
pag2$Plumage.prev = NA
for (i in 1:nrow(pag2)) {
  pag.male = p.a.g %>% filter(Fwnumber==pag2$Fwnumber[i]) %>% filter(Year==(as.numeric(as.character(pag2$Year[i]))-1))
  if(nrow(pag.male)==0) {pag2$Plumage.prev[i]=NA} else {
    pag2$Plumage.prev[i] = as.character(pag.male$Plumage)
  }
}
pag2$Ornamented.prev = NA
for (i in 1:nrow(pag2)) {
  if(is.na(pag2$Plumage.prev[i])) {pag2$Ornamented.prev[i]=NA} else {
  if(pag2$Plumage.prev[i]=="Dull") {pag2$Ornamented.prev[i]="Brown"} else {pag2$Ornamented.prev[i]="Ornamented"}}
  }
pag2$Ornamented.prev = factor(pag2$Ornamented.prev,levels=c("Brown","Ornamented")) 


##Determine if they were paired in the previous year
pag2$Paired.prev = NA
for (i in 1:nrow(pag2)) {
  pag.male = p.a.g %>% filter(Fwnumber==pag2$Fwnumber[i]) %>% filter(Year==(as.numeric(as.character(pag2$Year[i]))-1))
  if(nrow(pag.male)==0) {pag2$Paired.prev[i]=NA} else {
    if(pag.male$Status=="Male.fwno") {pag2$Paired.prev[i]="Yes"} else {pag2$Paired.prev[i]="No"}
  }
}


##Determine if males were paired to the same female as the previous year
#First get groups.m in descending order - want male's last female, whereas before wanted his first female of the season
groups.m2 = groups.m %>% arrange(Status,desc(Date.Created)) #Order by status and date created - status first
#Remove duplicates based on ID and year, getting highest breeding status group for each individual in each year. 
groups.m2d = groups.m %>% distinct(Year,Fwnumber,.keep_all = T) 

pag2$Paired.prev.sameF = NA
for (i in 1:nrow(pag2)) {
  if(is.na(pag2$Paired.prev[i]) | pag2$Paired.prev[i]=="No") {pag2$Paired.prev.sameF[i]=NA} else {
    brg = groups.m2d %>% filter(Fwnumber==pag2$Fwnumber[i],Year==(as.numeric(as.character(pag2$Year[i]))-1))
    female.prev = groups %>% filter(Group==brg$Group) %>% select(Female.fwno)
    female.current = groups %>% filter(Group==pag2$Group[i]) %>% select(Female.fwno)
    if(female.prev[1,]==female.current[1,]) {pag2$Paired.prev.sameF[i]="Yes"} else {pag2$Paired.prev.sameF[i]="No"}
    }
  }



####Why do some 2yo's molt and others do not?####

#Not tons of data to work with here, only 5 males did not molt into nuptial plumage at 2 years of age and only one male from 
#years we have social data did this - GBI
table(pag2$Ornamented)


###How does previous ornamentation influence molt to bright? First variable in table on side, second on top:
table(pag2$Ornamented,pag2$Ornamented.prev)
#The five males that were brown as 2yos were brown in the previous year. All males that were ornamented as 1-year-olds 
#were ornamented in their second year (25 males). 74 males that were brown as 1-year-olds were ornamented as 2-year-olds. 

##Fisher exact test:
fisher.test(pag2$Ornamented,pag2$Ornamented.prev)
#Shows that the difference in proportions is not significant. So previously ornamented males being ornamented as 2yos could be due to 
#random chance, mainly because very few males were brown in their first year. 


###How does pairing status as a 2-year-old influence molt to bright? 
table(pag2$Ornamented,pag2$Paired)
#There were two bright helper 2yo males, probably these males paired up then their mate disappeared and they re-joined their parents. 
#We've seen that happen with a few other birds. All other birds were paired, including brown and ornamented birds.




###How does pairing status as a 1yo influence molt to bright? 

#Because fisher test showed no effect of ornamented prev on ornamented as a 2yo, use all birds. 
#Get only birds with paired.prev
pag2op = pag2 %>% filter(!is.na(Ornamented.prev))

#Model
pag2op.m = glmer(Ornamented~Paired.prev + (1|Year),data=pag2op,family="binomial")
summary(pag2op.m)

#Residuals
simulationOutput <- simulateResiduals(fittedModel = pag2op.m)
plot(simulationOutput, asFactor = T) #Plot residual plots 
plotResiduals(simulationOutput,pag2op$Paired.prev, quantreg = T) 
testDispersion(simulationOutput)

#Likelihood ratio test
pag2op.m2 = glmer(Ornamented~1 + (1|Year),data=pag2op,family="binomial")
anova(pag2op.m,pag2op.m2)

#Paper: 2yo males that were paired as 1yos were maybe slightly more likely to molt into nuptial plumage than 
#males that were helpers as 1-year-olds (X2=2.72,p=0.099). 



###How does being paired to the same female influence molt into bright? 
pag2opsf = pag2op %>% filter(!is.na(Paired.prev.sameF))
table(pag2opsf$Ornamented,pag2opsf$Paired.prev.sameF)
pag2opsf.m = glmer(Ornamented~Paired.prev.sameF + (1|Year),data=pag2opsf,family="binomial")
summary(pag2opsf.m)

#Residuals
simulationOutput <- simulateResiduals(fittedModel = pag2opsf.m)
plot(simulationOutput, asFactor = T) #Plot residual plots 
plotResiduals(simulationOutput,pag2opsf$Paired.prev.sameF, quantreg = T) 
testDispersion(simulationOutput)

#Likelihood ratio test
pag2opsf.m2 = glmer(Ornamented~1 + (1|Year),data=pag2opsf,family="binomial")
anova(pag2opsf.m,pag2opsf.m2)

#Paper: Within paired males that were dull in their first year, there was no effect of whether or not they were with their same female or 
#with a new female in their second year on their likelihood of molting into nuptail plumage (X2=0.01,p=0.926)



###Write pag2 ornamented males to file for when males molted analysis
##Before writing, get rid of colors column - colors comes from ind history which is the male's last colors, not necessarily his colors
#when he was 2 years old
pag2.ornamented = pag2 %>% filter(Ornamented=="Ornamented") %>% select(-Colors)
#write.csv(pag2.ornamented,here::here("Output files","twoyo_ornamented.csv"),row.names = F)

