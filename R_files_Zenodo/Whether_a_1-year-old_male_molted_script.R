###Whether a one-year-old male molted script

#Use this script for analyses looking at which males molt into nuptial plumage and why some males molt and others do not


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

##Load in nestling hatch dates
hatch = read.csv(here::here("Input files","All RBFW hatch dates.csv"))

##Load list of birds in removal communities in 2017 season and their FWnumbers
removals = read.csv(here::here("Input files","RemovalorControl 2017.csv"))
ages17 = read.csv(here::here("Input files","Ages and Sexes 2017.csv"))





####Get plumage and age files into long format so I can merge####
library(dplyr)
library(reshape2)

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

##Remove Removal males - males in removal communities
removals = merge(removals,ages17,by="Bird",all.x=T) #Merge community ID with ages to get FWnumber. This ages file has colorbands specific to what they were in 2017 season
removals = removals %>% filter(TreatmentGroup=="Removal") %>% select(TreatmentGroup,FWNo)
#Removals happened in 2017 season
for (i in 1:nrow(p.a)) {
  if(p.a$Fwnumber[i] %in% removals$FWNo & p.a$Year[i]==2017) {p.a$Drop[i]="drop"}
}

##Remove both
p.a = p.a %>% filter(Drop!="drop") %>% select(-Drop)




####How does age influence whether males molt into nuptial plumage?####

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



###Histogram of ages represented
library(ggplot2)
library(cowplot)
ggplot(data=p.a.g,aes(x=Age)) + geom_histogram(binwidth = 0.5) + theme_cowplot() +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7))



###Aggregate by plumage
pag.percent = p.a.g %>% group_by(Age,Plumage) %>% summarise(n=n()) %>% mutate(freq=n/sum(n))

#Add levels for dull and interm at other age classes
ages.levels = c(3:7)  
plum.levels = c("Dull","Interm")
pag.levels = expand.grid(ages.levels,plum.levels) #Get all combinations
colnames(pag.levels) = c("Age","Plumage")
pag.levels = pag.levels %>% mutate(n=0,freq=0) #Add other columns
pag.percent = bind_rows(pag.percent,pag.levels) #Combine dataframes
pag.percent$Plumage = factor(pag.percent$Plumage,levels=c("Dull","Interm","Bright"),ordered = T) #Order factors for plotting



###Plot plumage by age
ggplot(data=pag.percent,aes(x=Age,y=freq,fill=Plumage)) + geom_bar(stat='identity',position = 'dodge',color="black") + 
  theme_cowplot() + scale_x_continuous(breaks=c(1,2,3,4,5,6,7)) + 
  scale_fill_manual(values=c("white","grey","black"),labels=c("Dull","Intermediate","Bright")) + ylab("Percentage") 



###Plot plumage by age but group 3+ birds since all the same
#Assign 3 year olds and older to Age 3
p.a.g$Age.threeplus = NA
for (i in 1:nrow(p.a.g)) {
  if(p.a.g$Age[i]>=3) {p.a.g$Age.threeplus[i]=3} else {p.a.g$Age.threeplus[i]=p.a.g$Age[i]}
}

##Get ornamented or not in p.a.g 
p.a.g$Ornamented = NA
for (i in 1:nrow(p.a.g)) {
  if(p.a.g$Plumage[i]=="Dull") {p.a.g$Ornamented[i]="Brown"} else {p.a.g$Ornamented[i]="Ornamented"}
}
p.a.g$Ornamented = factor(p.a.g$Ornamented,levels=c("Brown","Ornamented")) 


#Aggregate by plumage again
pag.percent.3plus = p.a.g %>% group_by(Age.threeplus,Ornamented) %>% summarise(n=n()) %>% mutate(freq=n/sum(n))
threeplus.extra = data.frame(matrix(nrow=1,ncol=0))
threeplus.extra = threeplus.extra %>% mutate(Age.threeplus=3,Ornamented=c("Brown"),n=0,freq=0)
pag.percent.3plus = bind_rows(pag.percent.3plus,threeplus.extra)
pag.percent.3plus$Ornamented = factor(pag.percent.3plus$Ornamented,levels=c("Brown","Ornamented"),ordered = T) #Order factors for plotting
pag.percent.3plus$Age.threeplus = as.factor(pag.percent.3plus$Age.threeplus)

#Get ornamented by age column for fill
pag.percent.3plus$Age.orn = paste(pag.percent.3plus$Age.threeplus,pag.percent.3plus$Ornamented,sep="")

#Plot
ggplot(data=pag.percent.3plus,aes(x=Age.threeplus,y=freq,fill=Ornamented)) + 
  geom_bar(stat='identity',position = 'dodge',color="black",size=1) + 
  theme_cowplot() + #scale_fill_manual(values=c("white","gray"),labels=c("Brown","Ornamented")) + 
  scale_fill_manual(values=c("#C6A375","#D82F2B")) +
  ylab("Frequency") + xlab("Age") + scale_x_discrete(labels=c("1","2","3+"))


##Write p.a.g for rain and molt date plots - need to know how many males of each age class did not make it to bright plumage in each year
#write.csv(p.a.g,here::here("Output files","all_males_moltornot.csv"),row.names = F)




####Set up one year old dataframe####

pag1 = p.a.g %>% filter(Age==1)
pag1$Year = as.factor(pag1$Year)


##Assign whether a male was Ornamented or not (binary variable)
pag1$Ornamented = NA
for (i in 1:nrow(pag1)) {
  if(pag1$Plumage[i]=="Dull") {pag1$Ornamented[i]="Brown"} else {pag1$Ornamented[i]="Ornamented"}
}
pag1$Ornamented = factor(pag1$Ornamented,levels=c("Brown","Ornamented")) 


##Determine if one year olds were ever paired in the season they were that age
pag1$Paired = NA
for (i in 1:nrow(pag1)) {
  if(pag1$Status[i]=="Male.fwno") {pag1$Paired[i]="Yes"} else {pag1$Paired[i]="No"}
}
pag1$Paired = as.factor(pag1$Paired)


##Get hatch dates for one-year-old males
hatch.m = melt(hatch,id.vars = 1:5) #melt
hatch.m = hatch.m %>% filter(value!="") #Remove blanks
colnames(hatch.m)[7] = "Fwnumber" 
hatch.m = hatch.m %>% select(Fwnumber,Hatch.date,Group) #Select only columns needed
hatch.m = hatch.m %>% filter(Hatch.date!="") #Remove nestlings without hatch dates
hatch.m$Hatch.date = as.Date(hatch.m$Hatch.date,"%m/%d/%y") #convert to date
colnames(hatch.m)[3] = "Natal.group"
pag1 = merge(pag1,hatch.m, by="Fwnumber",all.x=T) #Merge with pag1

#Get hatch date in julian date form
library(lubridate)
pag1$Hatch.jdate = yday(pag1$Hatch.date)

#If the hatch date is in January, add 365
for (i in 1:nrow(pag1)) { if(!is.na(pag1$Hatch.date[i])) {
  if(month(pag1$Hatch.date[i])==1) {pag1$Hatch.jdate[i]=pag1$Hatch.jdate[i]+365}}
}


###Get number of male helpers at natal nest for one-year-old males
groups.h = groups %>% select(Group,contains("Helper")) #select helper columns
groups.hm = melt(groups.h,id.vars=1) #melt
groups.hm$value[groups.hm$value==""] <- NA #make sure blank cells are NA

#Make sure the helpers are male or unknown sex (assume male) - so remove female helpers
sexes = plumages %>% select(Fwnumber,Sex) #Get a list of sexes
groups.hm = merge(groups.hm,sexes,by.x="value", by.y="Fwnumber",all.x=T) #Merge by helper ID
colnames(groups.hm)[1] = "Fwnumber"

#For female helpers, replace Fwnumber with NA so it doesn't get counted in next step
for (i in 1:nrow(groups.hm)) { if(!is.na(groups.hm$Sex[i])) {
  if(groups.hm$Sex[i]=="F") {groups.hm$Fwnumber[i]=NA}}
}

#Get counts of helper males associated with each nest
groups.hm.counts = groups.hm %>% group_by(Group) %>% summarise(Natal.group.helpers = sum(!is.na(Fwnumber))) #Get counts of number of helpers in each group
pag1 = merge(pag1,groups.hm.counts,by.x="Natal.group",by.y="Group",all.x=T) #Merge with pag1

#Get binary version of helpers at nest or not 
pag1$Natal.group.helpers.binary = NA
for (i in 1:nrow(pag1)) { if(!is.na(pag1$Natal.group.helpers[i])) {
  if(pag1$Natal.group.helpers[i]>=1) {pag1$Natal.group.helpers.binary[i]="Yes"} else
  {pag1$Natal.group.helpers.binary[i]="No"} }
}
pag1$Natal.group.helpers.binary = as.factor(pag1$Natal.group.helpers.binary)




####For paired males, get information on the female they were first paired to

##For GZR, the database lists him as paired to his mother (BWR, age 5), but behaviorally it appeared he was paired with HIB a 1+min female
#before eventually pairing with RVR older female next door. GZR, BWR, and HIB were in non-breeding group together, in sightings have 
#instance written down that GZR and HIB were chasing RVR and duetting, then GZR switched to duetting with RVR and during the whole thing BWR
#did not take part. Think that shows that GZR was paired with HIB and not his mother. His mother eventually paired with LZB when LZB appeared

###So for now, change GZR's female to HIB - age 1, age exact = min. 
b771 = groups %>% filter(Group=="B771")
groups = groups %>% filter(Group!="B771")
b771$Female.fwno = b771$Other1.fwno
groups = rbind(groups,b771)
groups = groups %>% arrange(Date.Created)

#Now get female data
pag1$Female.age = NA  
pag1$Female.age.exact = NA 
pag1$Female.prev = NA
for (i in 1:nrow(pag1)) {if(pag1$Paired[i]=="No") {
  pag1$Female.age[i]=NA 
  pag1$Female.age.exact[i] = NA} else {
    brg = groups %>% filter(Group==pag1$Group[i]) #Get male's breeding group
    female = brg %>% select(Female.fwno)
    female.age = ages.pyd %>% filter(Year==brg$Year,FWno==brg$Female.fwno) #Get female's data
    if(nrow(female.age)==0) { #if UNB or not in ages file, then = NA
      pag1$Female.age[i]=NA 
      pag1$Female.age.exact[i] = NA} else { 
        pag1$Female.age[i] = female.age$Age
        pag1$Female.age.exact[i] = as.character(female.age$Age.exact) #Have to use character here because it turns this column to an integer
        female.prev = groups.md %>% filter(Fwnumber==female[1,],Year==(brg$Year-1)) %>% select(Status) #Get female status in previous year if we have it
        pag1$Female.prev[i] = as.character(female.prev[1,])
        }} 
}
#Convert back to factor
pag1$Female.age.exact = as.factor(pag1$Female.age.exact)
pag1$Female.prev = as.factor(pag1$Female.prev)
#Checked the males that were paired but female data = NA, all had unknown or UNB females. 







####How does pairing status influence whether one year olds molt into bright?####

#Realized that all bright one-year-old males were paired, but not all paired males were bright
#No helpers were bright
#Plot this

##First aggregate data by pairing status
pag1.percent.paired = pag1 %>% group_by(Paired,Ornamented) %>% summarise(n=n()) %>% mutate(freq=n/sum(n)) #summarize
pag1.pp.extra = data.frame(matrix(nrow=1,ncol=0)) #get other levels
pag1.pp.extra = pag1.pp.extra %>% mutate(Paired="No",Ornamented=c("Ornamented"),n=0,freq=0) #get other levels
pag1.percent.paired = bind_rows(pag1.percent.paired,pag1.pp.extra) #combine
pag1.percent.paired$Ornamented = factor(pag1.percent.paired$Ornamented, levels=c("Brown","Ornamented"),ordered = T) #set order of levels

#Plot
ggplot(data=pag1.percent.paired,aes(x=Paired,y=freq,fill=Ornamented)) + 
  geom_bar(stat='identity',position = 'dodge',color="black",size=1) + 
  theme_cowplot() + scale_fill_manual(values=c("#C6A375","#D82F2B"),labels=c("Brown","Ornamented"),name="Plumage") + 
  ylab("Frequency") + xlab("Breeding status") + scale_x_discrete(labels=c("Helper","Paired"))

##Are these ratios different from chance? 
fisher.test(pag1$Paired,pag1$Ornamented)
#Yes different from chance

##Paper: Only paired males were ornamented and this difference was significant (fisher exact test p<0.001)



#So further analyses should look within paired males - why are some paired males bright and other not? 
pag1.p = pag1 %>% filter(Paired=="Yes")

#Write pag1.p to .csv to use in 1yo male pairing date script
#write.csv(pag1.p,here::here("Output files","oneyo_paired.csv"),row.names = F)

#But first look at why some males are paired and others are not





####Why do some males pair but others do not?####

#Get paired in binary form, (binary variable; 1 = "success", 0 = "failure")
pag1$paired.binary = NA
for (i in 1:nrow(pag1)) {
  if(pag1$Paired[i]=="Yes") {pag1$paired.binary[i]="1"} else {pag1$paired.binary[i]="0"}
}
pag1$paired.binary = as.numeric(pag1$paired.binary) #Make it numeric for binomial models search "binomial" to see details


####Look at effects of hatch date and helpers on whether a male paired or not
#Get birds with natal data
pag1.n = pag1 %>% filter(Hatch.jdate!="")

#Plot Hatch date
ggplot(data=pag1.n,aes(x=Hatch.jdate,y=paired.binary)) + geom_point() + geom_smooth() + theme_cowplot()

#Plot number of helpers at natal nest
ggplot(data=pag1.n,aes(x=Natal.group.helpers.binary,y=paired.binary)) + geom_point() + geom_boxplot() + theme_cowplot()


###Model 
library(glmmTMB)
library(lme4)
pag1.nm1 = glmer(paired.binary~scale(Hatch.jdate) + scale(Natal.group.helpers) + (1|Natal.group),data=pag1.n,family="binomial")
summary(pag1.nm1) #Year as a random effect was singular - explained no variation

#Model with interaction
pag1.nm2 = glmer(paired.binary~scale(Hatch.jdate) + scale(Natal.group.helpers) + scale(Hatch.jdate):scale(Natal.group.helpers) + (1|Natal.group),data=pag1.n,family="binomial")
summary(pag1.nm2) #Year as a random effect was singular - explained no variation

#Test if interaction improved model
anova(pag1.nm1,pag1.nm2)
#Interaction did not improve the model, go with pag1.nm1


###Check residuals of model
#These steps from DHARMa vignette 
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = pag1.nm1)
plot(simulationOutput, asFactor = F) #Plot residual plots
plotResiduals(simulationOutput,pag1.n$Hatch.jdate, quantreg = T) #Plot again
#This plot looks decent
plotResiduals(simulationOutput,pag1.n$Natal.group.helpers, quantreg = T,asFactor = F)
#Looks good
testDispersion(simulationOutput)
#Looks good
simulationOutput = recalculateResiduals(simulationOutput , group = pag1.n$Year)
testDispersion(simulationOutput)
#Looks good
simulationOutput = recalculateResiduals(simulationOutput , group = pag1.n$Natal.group)
testDispersion(simulationOutput)
#Looks good


###Likelihood ratio tests to get significance of individual fixed effects
pag1.n.noHJD = glmer(paired.binary~scale(Natal.group.helpers) + (1|Natal.group),data=pag1.n,family="binomial")
pag1.n.noNGH = glmer(paired.binary~scale(Hatch.jdate) + (1|Natal.group),data=pag1.n,family="binomial")

#LRT
anova(pag1.nm1,pag1.n.noHJD)
anova(pag1.nm1,pag1.n.noNGH)
#No effect of hatch date or the number of helpers at the natal group on whether or not a 1-year-old male paired
#Given the above graphs, no surprise that these fixed effects are not important

#Plot effects for fun
library(effects)
library(MASS)
plot(Effect(focal.predictors = "Hatch.jdate",pag1.nm1))
plot(Effect(focal.predictors = "Natal.group.helpers",pag1.nm1))
detach("package:MASS") #detach MASS package so "select" function in dplyr still works


###Paper: 
#Neither hatch date nor the number of male helpers at the 1-year-old male's natal nest were predictive of whether a 
#1-year-old male paired during the breeding season or not (hatch date: X2<0.01, p=1.00; natal helpers: X2<0.01, p=0.941). 

#Don't need plots of model effects

#Methods - standardized and centered all variables to assist with model convergence. Binomial model, logit. 








####How do hatch date, number of helpers influence whether one-year-olds molt into bright?####

##Get birds with nest data - hatch date and helper data
pag1.pn = pag1.p %>% filter(!is.na(Hatch.date))
pag1.pn$Plumage = factor(pag1.pn$Plumage,levels=c("Dull","Interm","Bright"),ordered = T)

##Get Ornamented.binary
pag1.pn$Ornamented.binary = NA
for (i in 1:nrow(pag1.pn)) {
  if(pag1.pn$Ornamented[i]=="Ornamented") {pag1.pn$Ornamented.binary[i]=1} else {pag1.pn$Ornamented.binary[i]=0}
}

#Plot Hatch date
ggplot(data=pag1.pn,aes(x=Ornamented,y=Hatch.jdate)) + geom_boxplot() + geom_point()  + theme_cowplot()
ggplot(data=pag1.pn,aes(x=Hatch.jdate,y=Ornamented.binary)) + geom_point() + geom_smooth() + theme_cowplot()

#Plot number of helpers
ggplot(data=pag1.pn,aes(x=Ornamented,y=Natal.group.helpers)) + geom_boxplot() + geom_point() + theme_cowplot()
ggplot(data=pag1.pn,aes(x=Natal.group.helpers,y=Ornamented.binary)) + geom_point() + theme_cowplot()

##See if hatch date and helpers are related
pag1.pn$Natal.group.helpers.cat = factor(pag1.pn$Natal.group.helpers,levels=c("0","1","2"))
ggplot(data=pag1.pn,aes(x=Natal.group.helpers.cat,y=Hatch.jdate)) + geom_boxplot() + geom_point() + theme_cowplot()
#Could be something going on? 



####Model hatch date and number of helpers 

##Scale variables here so their importance can be compared using their coefficients

##For random effects, - nest natal group ID within year, since natal groups are specific to specific years - doesn't actually 
#change anything though since groups are labeled uniquely. Could have it not nested and would get same results.

#Model
pag1.pn.m1 = glmer(Ornamented~scale(Hatch.jdate) + scale(Natal.group.helpers) + (1|Fwnumber) + (1|Year/Natal.group),data=pag1.pn,family="binomial")
summary(pag1.pn.m1)

####Check residuals
#Observation level random effect helped the residuals.
simulationOutput <- simulateResiduals(fittedModel = pag1.pn.m1)
plot(simulationOutput, asFactor = F) #Plot residual plots 
plotResiduals(simulationOutput,pag1.pn$Hatch.jdate, quantreg = T) 
testDispersion(simulationOutput)
#Look decent


####Likelhood ratio tests to get p-values of fixed effects

##Results to report in paper comes from here - these are a bit better than Wald Z scores in normal model summary output
#First get models with and without a variable each time
pag1.pn.m1.HJD = glmer(Ornamented~scale(Natal.group.helpers) + (1|Fwnumber) + (1|Year/Natal.group),data=pag1.pn,family="binomial")
pag1.pn.m1.NGH = glmer(Ornamented~scale(Hatch.jdate) + (1|Fwnumber) + (1|Year/Natal.group),data=pag1.pn,family="binomial")

#Likelihood ratio tests
anova(pag1.pn.m1,pag1.pn.m1.HJD)
anova(pag1.pn.m1,pag1.pn.m1.NGH)


####Plot odds ratio plots
library(MASS) #MASS is needed to run Effect function
plot(Effect(focal.predictors = "Hatch.jdate",pag1.pn.m1))
plot(Effect(focal.predictors = "Natal.group.helpers",pag1.pn.m1))
detach("package:MASS") #detach MASS package so "select" function in dplyr still works



####Paper 
#Include likelihood ratio tests in text of results, something like this: 
#Hatch date was an important predictor of whether a one-year old male molted into intermediate or bright plumage,
#(Hatch date: Wald's X2=7.87, p=0.005), while the number of helpers alone was nonsignificant (X2=1.20, p=0.158).

#Include odds ratio plots in supplemental materials?
#Or plot better prediction plots? 

#Methods - standardized and centered all variables to assist with model convergence and to allow us to assess the importance 
#of fixed effects and interactions simultaneously (Schielzeth, 2010).


####Plot prediction plot for hatch date
#Get prediction model format
hatch.predict = glmer(Ornamented~Hatch.jdate + Natal.group.helpers + (1|Fwnumber) + (1|Year/Natal.group),data=pag1.pn,family="binomial")
summary(hatch.predict)


#Function to get confidence intervals from Ben Bolker
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-crit*pred.se,
                conf.high=pred0+crit*pred.se))
}


#Get prediction dataframe - need to have columns for each of the fixed effects here
hatch.predict.df = expand.grid(Hatch.jdate=seq(from=230,to=370,by=1),Natal.group.helpers=mean(pag1.pn$Natal.group.helpers))

#Predict values 
hatch.predict.p = predict(hatch.predict,newdata=hatch.predict.df,type="response",re.form=NA)

#Get confidence intervals
hatch.predict.ci = easyPredCI(hatch.predict,newdata=hatch.predict.df)

#Combine prediction dataframes
hatch.predict.df2 = cbind(hatch.predict.df,molt.likelihood = hatch.predict.p,hatch.predict.ci)

#Plot
ggplot(data=hatch.predict.df2,aes(x=Hatch.jdate,y=molt.likelihood)) + 
  geom_point(data=pag1.pn,aes(x=Hatch.jdate,y=jitter(Ornamented.binary,amount=0.01)),
             alpha=0.7,size=2) +
  geom_ribbon(aes(x=Hatch.jdate,y=molt.likelihood,ymin=conf.low,ymax=conf.high),alpha=0.1) +
  geom_line(color="black",size=1.2) + theme_cowplot() + xlab("Hatch date (day of year)") +
  ylab("Likelihood of acquiring\\nornamented plumage")





####Plot prediction plot for number of helpers
#Get prediction model format
helpers.predict = glmer(Ornamented~Hatch.jdate + Natal.group.helpers + (1|Fwnumber) + (1|Year/Natal.group),data=pag1.pn,family="binomial")
summary(helpers.predict)

#Get prediction dataframe - need to have columns for each of the fixed effects here
helpers.predict.df = expand.grid(Hatch.jdate=mean(pag1.pn$Hatch.jdate),Natal.group.helpers=seq(from=0,to=2,by=1))

#Predict values 
helpers.predict.p = predict(helpers.predict,newdata=helpers.predict.df,type="response",re.form=NA)

#Get confidence intervals
helpers.predict.ci = easyPredCI(helpers.predict,newdata=helpers.predict.df)

#Combine prediction dataframes
helpers.predict.df2 = cbind(helpers.predict.df,molt.likelihood = helpers.predict.p,helpers.predict.ci)

#Plot
ggplot(data=helpers.predict.df2,aes(x=Natal.group.helpers,y=molt.likelihood)) + 
  geom_point(data=pag1.pn,aes(x=jitter(Natal.group.helpers,amount=0),y=jitter(Ornamented.binary,amount=0.02)),
             alpha=0.7,size=2) +
  geom_ribbon(aes(x=Natal.group.helpers,y=molt.likelihood,ymin=conf.low,ymax=conf.high),alpha=0.1) +
  geom_line(color="black",size=1.2) + theme_cowplot() + xlab("Number of male helpers at natal nest") +
  ylab("Likelihood of molt into ornamented plumage")









####How does the female they paired to influence whether one-year-olds molt into bright?####

##This section sets up the baseline model to be used in the sliding window analyses

####Look at female age data

##First aggregate data by female age using only exact ages
pag1p.female.age = pag1.p %>% filter(Female.age!="")
pag1p.female.age = pag1p.female.age %>% filter(Female.age.exact=="exact")
pag1p.female.age$Female.age.three.plus = NA
for (i in 1:nrow(pag1p.female.age)) {if(pag1p.female.age$Female.age[i]>=3) {pag1p.female.age$Female.age.three.plus[i]="3+"} else 
{pag1p.female.age$Female.age.three.plus[i]=pag1p.female.age$Female.age[i]}
}
pag1p.female.age.p = pag1p.female.age  %>% group_by(Female.age.three.plus,Plumage) %>% summarise(n=n()) %>% mutate(freq=n/sum(n)) #summarize
pag1p.female.age.p$Plumage = factor(pag1p.female.age.p$Plumage, levels=c("Dull","Interm","Bright"),ordered = T) #set order of levels


##What range of ages did one-year old males pair with?
ggplot(data=pag1p.female.age,aes(x=Female.age)) + geom_histogram(binwidth = 0.5) + theme_cowplot()

#Plot relationship between female age and molt to bright
ggplot(data=pag1p.female.age.p,aes(x=Female.age.three.plus,y=freq,fill=Plumage)) + 
  geom_bar(stat='identity',position = 'dodge',color="black") + 
  theme_cowplot() + scale_fill_manual(values=c("white","grey","black"),labels=c("Dull","Intermediate","Bright")) + 
  ylab("Frequency") + xlab("Age of Paired Female") + ylim(0,1)


##Aggregate based on whether female has bred before - using age and female previous breeding status
#If a female is age 1 min, remove. If a female is 2 or older but was a helper or other breeding status at best last year, then say not bred before
#All of the 2 year old and older females have previous statuses - see pag1p.female.bredb4
pag1p.female.bredb4 = pag1.p %>% filter(Female.age!="")
Fage.one.min.remove = pag1p.female.bredb4 %>% filter(Female.age==1,Female.age.exact=="min") #Remove males with females that are age=1, age.exact=min 
pag1p.female.bredb4 = pag1p.female.bredb4 %>% filter(!Fwnumber %in% Fage.one.min.remove$Fwnumber) #Don't know if these females have bred before 
pag1p.female.bredb4$Bredb4 = NA
for (i in 1:nrow(pag1p.female.bredb4)) {if(pag1p.female.bredb4$Female.age[i]==1) {pag1p.female.bredb4$Bredb4[i]="No"} else {
  if(pag1p.female.bredb4$Female.prev[i]=="Female.fwno" | is.na(pag1p.female.bredb4$Female.prev[i])) {pag1p.female.bredb4$Bredb4[i]="Yes"} else {
    pag1p.female.bredb4$Bredb4[i]="No"}
  }
}
pag1p.female.bredb4$Bredb4 = as.factor(pag1p.female.bredb4$Bredb4)
pag1p.female.bredb4.p = pag1p.female.bredb4 %>% group_by(Bredb4,Ornamented) %>% summarise(n=n()) %>% mutate(freq=n/sum(n)) #summarize
pag1p.female.bredb4.p$Ornamented = factor(pag1p.female.bredb4.p$Ornamented, levels=c("Brown","Ornamented"),ordered = T) #set order of levels

#Plot
ggplot(data=pag1p.female.bredb4.p,aes(x=Bredb4,y=freq,fill=Ornamented)) + 
  geom_bar(stat='identity',position = 'dodge',color="black",size=1) + 
  theme_cowplot() + scale_fill_manual(values=c("#C6A375","#D82F2B"),labels=c("Brown","Ornamented"),name="Plumage") + 
  ylab("Frequency") + xlab("Female breeding experience") + scale_x_discrete(labels=c("Inexperienced","Experienced")) + ylim(0,1)




####Model whether female breeding experience/pairing method influences plumage type

#Model this in two ways - one without nest data (hatch date and helper data) and one with those data (smaller dataset)


###Model without nest data - large dataset
pag1.fbb4.m1 = glmer(Ornamented~Bredb4 + (1|Year),data=pag1p.female.bredb4,family="binomial")
summary(pag1.fbb4.m1)

###Model with nest data - reduced dataset
pag1p.female.bredb4.h = pag1p.female.bredb4 %>% filter(Hatch.jdate!="")
pag1.fbb4.m2 = glmer(Ornamented~Bredb4 + scale(Hatch.jdate) + scale(Natal.group.helpers) + (1|Year/Natal.group),data=pag1p.female.bredb4.h,family="binomial")
summary(pag1.fbb4.m2)


##Result: Pairing method is important when using the larger dataset and not including nest data, but when hatch date and natal helpers
#are included in the model it's not important anymore. 
#Explore below:


##Model with reduced dataset but without nest fixed effects, see if reduced sample size or inclusion of nest fixed 
#effects is driving nonsignificance of pairing method
pag1.fbb4.m3 = glmer(Plumage~Bredb4 + (1|Year/Natal.group),data=pag1p.female.bredb4.h,family="binomial")
summary(pag1.fbb4.m3)
#When you take nest data out of the model, pairing method is still not important in the reduced dataset

##Is hatch date related to pairing method? 
ggplot(data=pag1p.female.bredb4.h,aes(x=Bredb4,y=Hatch.jdate)) + geom_boxplot() + geom_point() + theme_cowplot()
t.test(pag1p.female.bredb4.h$Hatch.jdate~pag1p.female.bredb4.h$Bredb4)
#No relationship between hatch date and pairing method
#So why is the bredb4 variable not important anymore? Try plotting percentage version using pag1p.female.bredb4.h




####Compare percentage graphs across the two datasets
#Get percentages for nest data dataset
pag1p.female.bredb4.hp = pag1p.female.bredb4.h %>% group_by(Bredb4,Plumage) %>% summarise(n=n()) %>% mutate(freq=n/sum(n)) #summarize
pag1p.female.bredb4.hp$Plumage = factor(pag1p.female.bredb4.hp$Plumage, levels=c("Dull","Interm","Bright"),ordered = T) #set order of levels

##Plot
library(ggpubr)
#Nest data excluded large dataset
A = ggplot(data=pag1p.female.bredb4.p,aes(x=Bredb4,y=freq,fill=Plumage)) + 
  geom_bar(stat='identity',position = 'dodge',color="black") + 
  theme_cowplot() + scale_fill_manual(values=c("white","grey","black"),labels=c("Dull","Intermediate","Bright")) + 
  ylab("Percentage") + xlab("Pairing method") + scale_x_discrete(labels=c("New pair","Filled vacancy")) + ylim(0,1) +
  ggtitle("Nest data excluded,large dataset")
#Nest data included (small dataset)
B = ggplot(data=pag1p.female.bredb4.hp,aes(x=Bredb4,y=freq,fill=Plumage)) + 
  geom_bar(stat='identity',position = 'dodge',color="black") + 
  theme_cowplot() + scale_fill_manual(values=c("white","grey","black"),labels=c("Dull","Intermediate","Bright")) + 
  ylab("Percentage") + xlab("Pairing method") + scale_x_discrete(labels=c("New pair","Filled vacancy")) + ylim(0,1) +
  ggtitle("Nest data included, small dataset")

##Compare
ggarrange(A,B)

##Result:
#When you subset to 1yo males with nest data, the percentage of bright males in new pairs goes up compared to 
#when dataframe was not subsetted to include males with nest data. 
#Ratios of plumage levels for filled vacancy stay almost exactly the same across the datasets
#By selecting for birds with hatch dates could be removing individuals that hatched after our field season. 

#See how numbers of males in each group compare across dataframe
table(pag1p.female.bredb4$Bredb4)
table(pag1p.female.bredb4.h$Bredb4)
#Groups decrease by 26 and 18 males, 44 total lost when nest data included in model. 
#Losing lots of new pair dull birds when restricting to males that have hatch dates. 




####Look to see if males missing nest data were hatched late or immigrated

##An individual could have no nest data if they hatched after the breeding field season ended or if they immigrated into our population.  
#If they were hatched late, they should be included in analyses, but if they immigrated, they could be less likely to molt for reasons
#other than what we can test here - maybe they took longer to find a mate? Were not socially accepted into the population? 

##Can test if the birds missing nest data hatched late by looking to see if their parents are in our population using relatedness scores.
#Males that hatched after the breeding field season should be highly related to a male and female in the population that were not 1 year old
#the season that the male was a 1-year-old (not the same age as the male).

##Load relatedness data and merge with birds missing nest data in fbb4 dataframe
rel = read.csv(here::here("Input files","RelatednessEstimates_8xdenovo80_14to19.csv"))
pag1p.female.bredb4.nonest = pag1p.female.bredb4 %>% filter(is.na(Hatch.jdate)) %>% select(Fwnumber,Year) #Get individuals without nest data in fbb4 dataframe
rel.nonest = rel %>% select(Pair,ind2,ind1,wang) #Get columns in different order 
colnames(rel.nonest)[2:3] = c("ind1","ind2")
rel.nonest = rbind(rel,rel.nonest) #Combine so all individuals in both columns

##See if I should be including 1-year-old males from all years in this relatedness analysis. I only ran paternity for 2014-2019 so 
#SNP data does not go back all the way to 2011. Do have some individuals from later years if they were active in any years between 2014-2019
fwnos = read.csv(here::here("Input files","FWno Ages by Year.csv")) #Load in fairywren numbers and ages
rel.nn.yeartest = rel.nonest %>% distinct(ind1) #Get a single individual ID for each individual
rel.nn.yeartest = merge(rel.nn.yeartest,fwnos,by.x="ind1",by.y="FWno") #Merge with fwnos to get years active for each individual with SNP data
rel.nn.yeartest$Year = factor(rel.nn.yeartest$Year,levels=c("2011","2012","2013","2014","2015","2016","2017","2018","2019"))
ggplot(data=rel.nn.yeartest,aes(x=Year)) + geom_bar() #Plot
#Looks like I don't have much SNP data from 2011-2013, not surprising since I didn't run individuals active in those years unless they were also active
#in later years. So do not include 1-yo males in relatedness testing from 2011-2013, good chance they may have hatched in our population but we 
#don't have SNP data for their parents. 

##Need to subset for individuals the males are related to that were at least 2 years old the year the male was a 1yo - want to get birds that could've
#been the male's parent
rel.nonest = merge(rel.nonest,pag1p.female.bredb4.nonest,by.x="ind1",by.y="Fwnumber") #Merge to filter for 1yo males without nest data
rel.nonest = rel.nonest %>% filter(!Year %in% c("2011","2012","2013"))
length(unique(rel.nonest$ind1)) #Have relatedness data for 45 1yo males out of 49 total without nest data
rel.nonest.ages = merge(rel.nonest,fwnos,by.x=c("ind2","Year"),by.y=c("FWno","Year")) #Merge ind2 by fairywren number and year to get that individuals age when the 1yo male was 1yo
rel.nonest.ages = rel.nonest.ages %>% filter(Age>=2) #Year corresponds year when 1yo male was 1yo, so select only individuals that are older

##Select top related male and top related female for year 1yo male and look at distributions of relatedness to see if parents in population
rel.nonest.agesF = rel.nonest.ages %>% filter(Sex=="F")
rel.nonest.agesF = rel.nonest.agesF %>% arrange(-wang) %>% distinct(ind1,.keep_all = T)
rel.nonest.agesM = rel.nonest.ages %>% filter(Sex=="M")
rel.nonest.agesM = rel.nonest.agesM %>% arrange(-wang) %>% distinct(ind1,.keep_all = T)

#Plot each 1-year-old male's top related female score
ggplot(data=rel.nonest.agesF,aes(x=wang)) + geom_histogram(binwidth = 0.03,fill="light blue",color="black") + theme_cowplot() + xlim(0,1) +
  ggtitle("1-year-old males missing nest data: \\n relatedness to top older female") + theme(plot.title = element_text(hjust = 0.5))
ggplot(data=rel.nonest.agesM,aes(x=wang)) + geom_histogram(binwidth = 0.03,fill="light blue",color="black") + theme_cowplot() + xlim(0,1) +
  ggtitle("1-year-old males missing nest data: \\n relatedness to top older Male") + theme(plot.title = element_text(hjust = 0.5))

#These results indicate that the birds with missing data were from our population, meaning they likely hatched after the breeding season 
#field work ended, as we don't miss many nests during the breeding season according to Mike. 
#Going to run both version of models below though just in case and to have if I ever need. Will report model without nest data in paper. 



####Check residuals of model with large dataset - hatch data not included
simulationOutput <- simulateResiduals(fittedModel = pag1.fbb4.m1)
plot(simulationOutput) #Plot residual plots 
plotResiduals(simulationOutput,pag1p.female.bredb4$Bredb4, quantreg = T) 
testDispersion(simulationOutput)



####Likelihood ratio tests

###Model without nest data
##Results to report in paper comes from here - these are a bit better than Wald Z scores in normal model summary output
#Get models without a variable each time - comparing to pag1.fbb4.ordm1 
pag1.fbb4.BB4 = glmer(Ornamented~1 + (1|Year),data=pag1p.female.bredb4,family="binomial")

#Likelihood ratio test
anova(pag1.fbb4.m1,pag1.fbb4.BB4)


####Plot odds ratio plots
library(MASS)
plot(Effect(focal.predictors = "Bredb4",pag1.fbb4.m1))
detach("package:MASS") #detach MASS package so "select" function in dplyr still works


####Paper

##No climate variables were important in predicting plumage - all best windows were false positives. So below results stand. 

#Results might look like: One-year-old males that paired with an older female by filling a vacany were more likely to molt
#into intermediate or red-black plumage than males that formed a new pair with 1-year-old females (X2=6.43, p=0.011). 

#Methods - standardized and centered all variables to assist with model convergence.

#Include odds ratio plot for model without nest data in supplemental materials? Need new axes limits




####Which dataset to use for sliding window climate analyses? 

#Could use full dataset - pag1.p or fbb4 dataset and would be able to include that fixed effect in the model. 

##Check that fbb4 and full dataset have similar ratio of plumage types
pag1p.percent = pag1.p %>% group_by(Plumage) %>% summarise(n=n()) %>% mutate(freq=n/sum(n))
pag1p.percent$Plumage = factor(pag1p.percent$Plumage,levels=c("Dull","Interm","Bright"),ordered = T)
pag1p.female.bredb4.percent = pag1p.female.bredb4 %>% group_by(Plumage) %>% summarise(n=n()) %>% mutate(freq=n/sum(n))
pag1p.female.bredb4.percent$Plumage = factor(pag1p.female.bredb4.percent$Plumage,levels=c("Dull","Interm","Bright"),ordered = T)

#Plot
A = ggplot(data=pag1p.percent,aes(x=Plumage,y=freq,fill=Plumage)) + geom_bar(stat='identity',position = 'dodge',color="black")+ 
  theme_cowplot() + scale_fill_manual(values=c("white","grey","black"),labels=c("Dull","Intermediate","Bright")) + ylab("Frequency") + ylim(0,1) +
  ggtitle("Full dataset")
B = ggplot(data=pag1p.female.bredb4.percent,aes(x=Plumage,y=freq,fill=Plumage)) + geom_bar(stat='identity',position = 'dodge',color="black")+ 
  theme_cowplot() + scale_fill_manual(values=c("white","grey","black"),labels=c("Dull","Intermediate","Bright")) + ylab("Frequency") + ylim(0,1) +
  ggtitle("Pairing method dataset")
ggarrange(A,B)

##Ratios of plumage types are nearly the same across the two datasets. Use the pairing method dataset so that the baseline model can be a bit more
#defined and hopefully the likelihood of getting a false positive will be lower than a baseline model with no fixed effects. 

###Write fbb4 csv file to use in sliding window analyses

#write.csv(pag1p.female.bredb4,row.names=F,here::here("Output Files","paired_1yo_males_fbb4.csv"))





####For 1yo males that filled a vacancy, how does timing of that event influence molt into nuptial plumage?####

##Get 1yo males that filled vacancy
pag1p.fbredb4.fv = pag1p.female.bredb4 %>% filter(Bredb4=="Yes")
table(pag1p.fbredb4.fv$Ornamented) #About half ornamented and half brown

##Get female's last group from previous breeding season and see if her male from the previous season was active in the current year
#Also see if the female's male from last year had a breeding group in the current season - would expect most males not to have one
#if the 1yo male and female's group is the female's first group
pag1p.fbredb4.fv$Female.lastmaleactive = NA
pag1p.fbredb4.fv$Female.lastmaleactive.brg = NA
for (i in 1:nrow(pag1p.fbredb4.fv)) {
  female = groups %>% filter(Group==pag1p.fbredb4.fv$Group[i]) %>% select(Female.fwno)
  brgs = groups %>% filter(Female.fwno==female[1,]) %>% filter(Year==(as.numeric(as.character(pag1p.fbredb4.fv$Year[i]))-1)) 
  female.prevmale = brgs %>% arrange(desc(Date.Created)) %>% slice(1) %>% select(Male.fwno)
  active = ages.pyd %>% filter(FWno==female.prevmale[1,],Year==pag1p.fbredb4.fv$Year[i])
  if(nrow(active)>0) {pag1p.fbredb4.fv$Female.lastmaleactive[i] = "Yes"} else {pag1p.fbredb4.fv$Female.lastmaleactive[i] = "No"}
  brgs.m = groups %>% filter(Male.fwno==female.prevmale[1,],Year==pag1p.fbredb4.fv$Year[i])
  if(nrow(brgs.m)>0) {pag1p.fbredb4.fv$Female.lastmaleactive.brg[i] = "Yes"} else {pag1p.fbredb4.fv$Female.lastmaleactive.brg[i] = "No"}
}

##See if the group that the 1yo male and female were in was the female's first breeding group from that season
pag1p.fbredb4.fv$Female.firstgroup = NA
for (i in 1:nrow(pag1p.fbredb4.fv)) {
  female = groups %>% filter(Group==pag1p.fbredb4.fv$Group[i]) %>% select(Female.fwno)
  brgs = groups %>% filter(Female.fwno==female[1,]) %>% filter(Year==pag1p.fbredb4.fv$Year[i]) 
  first = brgs %>% arrange(Date.Created) %>% slice(1) %>% select(Group)
  if(pag1p.fbredb4.fv$Group[i]==first[1,]) {pag1p.fbredb4.fv$Female.firstgroup[i]="Yes"} else {pag1p.fbredb4.fv$Female.firstgroup[i]="No"}
}

##Whether the 1yo's breeding group was the female's first breeding group is the better measure here, because divorces do happen, but
#interesting to look at whether the previous male is still around or not. 

####Plot
pag1p.fbredb4.fv.plot = pag1p.fbredb4.fv %>% group_by(Female.firstgroup,Ornamented) %>% summarise(n=n()) %>% mutate(freq=n/sum(n)) #summarize
pag1p.fbredb4.fv.plot$Ornamented = factor(pag1p.fbredb4.fv.plot$Ornamented, levels=c("Brown","Ornamented"),ordered = T) #set order of levels
pag1p.fbredb4.fv.plot$Female.firstgroup = factor(pag1p.fbredb4.fv.plot$Female.firstgroup, levels=c("Yes","No"),ordered = T)

ggplot(data=pag1p.fbredb4.fv.plot,aes(x=Female.firstgroup,y=freq,fill=Ornamented)) + 
  geom_bar(stat='identity',position = 'dodge',color="black",size=1) + 
  theme_cowplot() + scale_fill_manual(values=c("#C6A375","#D82F2B"),labels=c("Brown","Ornamented"),name="Plumage") + 
  ylab("Percentage") + xlab("Females first mate of season") + ylim(0,1)


####Model to see if filling a breeding vacancy earlier or later influences nuptial molt

#Table 
table(pag1p.fbredb4.fv$Female.firstgroup,pag1p.fbredb4.fv$Ornamented)

#Model
pag1.fbb4.fv = glmer(Ornamented~Female.firstgroup + (1|Year),data=pag1p.fbredb4.fv,family="binomial")
summary(pag1.fbb4.fv)


##Check residuals of model 
simulationOutput <- simulateResiduals(fittedModel = pag1.fbb4.fv)
plot(simulationOutput) #Plot residual plots 
plotResiduals(simulationOutput,pag1p.female.bredb4$Female.firstgroup, quantreg = T) 
testDispersion(simulationOutput)


###Likelihood ratio test
#Model without fixed effect
pag1.fbb4.nofv = glmer(Ornamented~1 + (1|Year),data=pag1p.fbredb4.fv,family="binomial")

anova(pag1.fbb4.fv,pag1.fbb4.nofv)

####Plot odds ratio plots
library(MASS)
plot(Effect(focal.predictors = "Female.firstgroup",pag1.fbb4.fv))
detach("package:MASS") #detach MASS package so "select" function in dplyr still works

###Paper: there is a trend towards males that filled a vacancy that was open during the non-breeding season or early breeding season
#being more likely to molt into nuptial plumage (X2 = 3.12, p=0.077). But this analysis is limited by broad time intervals, some of these
#1yo males could have filled the vacancy during the non-breeding season, others could have filled the vacancy during the breeding season.



##Now ask in more detail - does when males pair with a female influence whether they molt into nuptial plumage or not

##Do that in separate script that focuses on social network data.



####Do males that hatch earlier pair earlier?####
pag1p.fbredb4.fv.plot2 = pag1p.fbredb4.fv %>% filter(!is.na(Hatch.jdate))
ggplot(data=pag1p.fbredb4.fv.plot2,aes(x=Female.firstgroup,y=Hatch.jdate)) + geom_boxplot(fill="light gray") + geom_point() + theme_cowplot() +
  xlab("Filled vacancy timing") + scale_x_discrete(labels=c("Late","Early")) + ylab("Hatch date") + ylim(200,380)
pag1p.fbredb4.fv.plot2$Female.firstgroup = as.factor(pag1p.fbredb4.fv.plot2$Female.firstgroup)
t.test(pag1p.fbredb4.fv.plot2$Hatch.jdate~pag1p.fbredb4.fv.plot2$Female.firstgroup)


###Paper: Males that hatched earlier were more likely to be their female's first mate of the breeding season (df=8.8, t=2.76, p=0.022)




####Condition and likelihood of molt in 1-year-old males####

#This is non-breeding condition - still only using captures at day 200 and before (July 19th)

##Read in condition data
cond = read.csv(here::here("Output files","Condition scores 15_19.csv"))

##See how much data if merge with full dataset vs fbb4
pag1.pc = merge(pag1.p,cond,by=c("Year","Fwnumber")) #89 observations
pag1p.female.bredb4c = merge(pag1p.female.bredb4,cond,by=c("Year","Fwnumber")) #65 observations

##Merge captures for the same individual
pag1.pc2 = pag1.pc %>% group_by(Year,Fwnumber) %>% summarise(mtresid = mean(mtresid),mean.jdate = mean(jdate))
pag1.pc3 = merge(pag1.pc,pag1.pc2,by=c("Year","Fwnumber")) #merge with original to get other columns
pag1.pc3 = pag1.pc3 %>% distinct(Fwnumber,.keep_all = T) #Get one observation per individual

##Model with full dataset - Year was singular. Adding mean.jdate to the model did not change anything
pag1p.cond.m1 = glm(Ornamented~mtresid.y,data=pag1.pc3,family="binomial")
summary(pag1p.cond.m1)

##Check residuals
simulationOutput <- simulateResiduals(fittedModel = pag1p.cond.m1)
plot(simulationOutput) #Plot residual plots 
plotResiduals(simulationOutput,pag1.pc3$mtresid.y, quantreg = T) 
testDispersion(simulationOutput)
#Look good

##Likelihood ratio test
pag1p.cond.m1.noc = glm(Ornamented~1,data=pag1.pc3,family="binomial")
anova(pag1p.cond.m1.noc,pag1p.cond.m1,test="LRT")

##Paper: Residual condition was not important for determining whether or not a 1-year-old male molted into ornamented plumage
#(X2=-0.58,p=0.45)



####Is condition correlated with hatch date? 
pag1.pc3.hd = pag1.pc3 %>% filter(!is.na(Hatch.jdate))
ggplot(data=pag1.pc3.hd,aes(x=Hatch.jdate,y=mtresid.y)) + geom_point() + geom_smooth(method="lm") + theme_cowplot() + 
  xlab("Hatch date") + ylab("Residual condition (mass/tarsus)")
cor.test(pag1.pc3.hd$mtresid.y,pag1.pc3.hd$Hatch.jdate)

####Paper: Hatch date was not correlated with early non-breeding season condition score (n=36, r=0.11, p=0.51)



####Does hatch date influence residual condition within a year? 
library(climwin)
ggplot(data=pag1.pc3.hd,aes(x=Hatch.jdate,y=mtresid.y)) + geom_point() + geom_smooth(method="lm",aes(group=Year)) + theme_cowplot() + 
  xlab("Hatch date") + ylab("Residual condition (mass/tarsus)")
##Center residual condition on year to take out across-year variation 
pag1.pc3.hd = droplevels(pag1.pc3.hd)
pag1.pc3.hd$mtresid.yc = wgdev(pag1.pc3.hd$mtresid.y,pag1.pc3.hd$Year)
ggplot(data=pag1.pc3.hd,aes(x=Hatch.jdate,y=mtresid.yc)) + geom_point() + geom_smooth(method="lm",aes(group=Year)) + theme_cowplot() + 
  xlab("Hatch date") + ylab("Residual condition (mass/tarsus)")
pag1p.hc.m1 = lm(mtresid.yc~Hatch.jdate,data=pag1.pc3.hd) #Year was singular. Adding mean.jdate didn't change anything
summary(pag1p.hc.m1)

##Check residuals
library(car)
hist(resid(pag1p.hc.m1))
qqPlot(pag1p.hc.m1)
#Look good

##LRT
pag1p.hc.m1.noh = lm(mtresid.yc~1,data=pag1.pc3.hd) 
anova(pag1p.hc.m1,pag1p.hc.m1.noh,test="LRT")

####Paper: Hatch date did not influence residual condition within a year (n=36, x2<0.01, p=0.926)




