####Molt by July script####

##Using this script to look how rainfall influences when males molt across years
#Might be able to compare btw 3+ and 2yo males. 

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
library(climwin)

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






####Set up molt by July dataframe####

##Select 2yo and 3+ males
p.a = p.a %>% filter(Age>=2)

##Select only males that made it to ornamented plumage
p.a = p.a %>% filter(Plumage=="Interm" | Plumage=="Bright")

##Match up colors to bird using agesex files - this will be that individuals colors in the year of interest pag2$Year and will match to ind file
agesex15.s = agesex15 %>% select(Year,FWNo,Bird)
agesex16.s = agesex16 %>% select(Year,FWNo,Bird)
agesex17.s = agesex17 %>% select(Year,FWNo,Bird)
agesex18.s = agesex18 %>% select(Year,FWNo,Bird)
agesex19.s = agesex19 %>% select(Year,FWNo,Bird)
agesexall = rbind(agesex15.s,agesex16.s,agesex17.s,agesex18.s,agesex19.s)
paa = merge(p.a,agesexall,by.x=c("Fwnumber","Year"),by.y=c("FWNo","Year"))

##Make sure no birds changed color combos within a year
captures.p.test = captures.p %>% filter(Fwnumber %in% paa$Fwnumber)
captures.p.test = captures.p.test %>% distinct(Year,Fwnumber,Colors)
check = captures.p.test %>% group_by(Year,Fwnumber) %>% filter(n()>1)
#So in captures, need to change those birds to the combo they are in agesex files for that year 7776 = RGG, 8730 = BRV
check1 = captures.p %>% filter(Year==2017 & Colors=="RGX") 
check1$Colors = "RGG"
check2 = captures.p %>% filter(Year==2019 & Colors=="BRG")
check2$Colors = "BRV"
captures.p = rbind(check1,check2,captures.p)
captures.p = captures.p %>% distinct(Fwnumber,Date,.keep_all = T) %>% filter(Fwnumber %in% paa$Fwnumber) %>% arrange(Date) #Remove later values 

##Make sure capture colors are the same as colors in agesex files for that year
capture.colorcheck = paa %>% select(Fwnumber, Year, Bird)
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


##Get birds of interest - birds in paa
paa.15 = paa %>% filter(Year=="2015")
ind15.s = ind15.s %>% filter(Bird %in% paa.15$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2015")
paa.16 = paa %>% filter(Year=="2016")
ind16.s = ind16.s %>% filter(Bird %in% paa.16$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2016")
paa.17 = paa %>% filter(Year=="2017")
ind17.s = ind17.s %>% filter(Bird %in% paa.17$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2017")
paa.18 = paa %>% filter(Year=="2018")
ind18.s = ind18.s %>% filter(Bird %in% paa.18$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2018")
paa.19 = paa %>% filter(Year=="2019")
ind19.s = ind19.s %>% filter(Bird %in% paa.19$Bird) %>% filter(!is.na(Molt)) %>% mutate(Year="2019")

##Combine, order, and get julian date
ind.all = rbind(ind15.s,ind16.s,ind17.s,ind18.s,ind19.s)
ind.all = ind.all %>% arrange(Year,Bird,Date)
ind.all$jdate = yday(ind.all$Date)

##If date is in January, make Julian date +365
for (i in 1:nrow(ind.all)) {
  if(month(ind.all$Date[i])==1) {ind.all$jdate[i] = ind.all$jdate[i] + 365}
}



###Measure difference in time between last dull observation and first intermediate/first bright
##Filter paa for only birds with molt data
paas = paa %>% filter(Bird %in% ind.all$Bird) 

##For loop to get difference in dates
paas$first.seen = NA
paas$times.seen = NA
paas$last.dull = NA
paas$first.int = NA
paas$last.int = NA
paas$first.bright = NA
paas$first.plum.afterdull = NA
paas$first.bright.plum = NA
paas$first.plum = NA
for (i in 1:nrow(paas)) {
  bird = paas$Bird[i]
  year = paas$Year[i]
  ind.bird = ind.all %>% filter(Bird==as.character(bird),Year==year) %>% arrange(jdate)
  paas$first.plum[i] = ind.bird$Molt[1]
  paas$first.seen[i] = ind.bird$jdate[1] 
  paas$times.seen[i] = nrow(ind.bird)
  ind.dull = ind.bird %>% filter(Molt<33) %>% arrange(desc(jdate))
  paas$last.dull[i] = ind.dull$jdate[1]
  ind.int = ind.bird %>% filter(Molt>=33,Molt<=66) %>% arrange(jdate)
  paas$first.int[i] = ind.int$jdate[1]
  ind.int2 = ind.int %>% arrange(desc(jdate))
  paas$last.int[i] = ind.int2$jdate[1]
  ind.bright= ind.bird %>% filter(Molt>66) %>% arrange(jdate)
  paas$first.bright[i] = ind.bright$jdate[1]
  paas$first.bright.plum[i] = ind.bright$Molt[1]
  ind.firstplumAD = ind.bird %>% filter(Molt>=33) %>% arrange(jdate)
  paas$first.plum.afterdull[i] = ind.firstplumAD$Molt[1]
}

#Remove records where bird not seen that year
paas = paas %>% filter(first.plum!="")

##Birds either need to be seen as bright by July 1, or they need to be seen as dull after July 1.
paas.july = paas %>% filter(last.dull>182 | (first.int<=182 | first.bright<=182))

##Determine if a male was bright by July 1st - July 1st is usually the date when we've seen most of the males
paas.july$firstbi = NA
for (i in 1:nrow(paas.july)) {if(is.na(paas.july$first.int[i])) {paas.july$firstbi[i]=paas.july$first.bright[i]} else {
  paas.july$firstbi[i]=paas.july$first.int[i]}
  }
paas.july$brightbyjuly = NA
for (i in 1:nrow(paas.july)) {if(paas.july$firstbi[i]<=182) {paas.july$brightbyjuly[i]="Yes"} else {paas.july$brightbyjuly[i]="No"}}
paas.july$brightbyjuly = as.factor(paas.july$brightbyjuly)
paas.july$Age = as.factor(paas.july$Age)


##Look at how many males of each age class were ornamented by July
table(paas.july$brightbyjuly,paas.july$Age)
paas.july.sum = paas.july %>% group_by(Age,brightbyjuly,.drop=F) %>% tally()
ggplot(data=paas.july.sum,aes(x=Age,y=n,fill=brightbyjuly)) + geom_bar(stat="identity",position=position_dodge(),color="black") + 
  scale_fill_manual(values=c("#C6A375","#D82F2B"),name="Ornamented\\nby July") + theme_cowplot() + ylab("Count")






####Sliding windows for molt by July####

##Load rainfall data
rain = read.csv(here::here("Input files","rain_2004-2018.csv"))
rain$Date = as.Date(rain$Date,"%m/%d/%y")

##Baseline model 
bm = glmer(brightbyjuly~1 + (1|Fwnumber),family="binomial",data=paas.july)
summary(bm)

##Check residuals
simulationOutput <- simulateResiduals(fittedModel = bm)
plot(simulationOutput, asFactor = F) #Plot residual plots 
testDispersion(simulationOutput)

##Add a date column
paas.july$date = as.Date(paste((paas.july$Year-1),"1","10",sep="-"))


####Sliding window - nothing beyond 200 - tried up to 400.

##Adding Fwnumber as a random effect led to major convergence errors and gave same window

##Get response variable as 0/1 to work for k-fold cross validation
paas.july$brightbyjuly2 = NA
for(i in 1:nrow(paas.july)) {if(paas.july$brightbyjuly[i]=="Yes") {paas.july$brightbyjuly2[i] = 1} else {paas.july$brightbyjuly2[i] = 0}}

##k-fold cross validation with k=10 returned the same best window and same model statistics. So not using k-fold in final model 
#to speed up randomization process. 

#Get age as integer
paas.july$Age.int = as.integer(as.character(paas.july$Age))

# rain.july = slidingwin(xvar=list(Rain = rain$amountcm),
#                        #k=10,
#                        cdate = rain$Date,
#                        bdate = paas.july$date, #This will get the correct year
#                        baseline = glm(brightbyjuly2~Age.int,family="binomial",data=paas.july),
#                        type = "absolute",
#                        cinterval = "day",
#                        range = c(200,0),
#                        refday = c(1,7),
#                        exclude = c(20,-1),
#                        stat = "sum",
#                        fun = "lin")

##Save model file
#saveRDS(rain.july,here::here("Output files","moltbyjuly_200_0_rainfall_norandom_age.csv"))

##Read in model file
rain.july = readRDS(here::here("Output files","moltbyjuly_200_0_rainfall_norandom_age.csv"))
#Look to see what windows are of best model 
head(rain.july[[1]]$Dataset,n=20)
#Get dataset
rain.july.dataset = rain.july[[1]]$Dataset
#Check summary
summary(rain.july[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=rain.july.dataset)
#Plotall
#plotall(rain.july.dataset)
#Get best model data
rain.july.bestdata = rain.july[[1]]$BestModelData
rain.july.bestdata$Fwnumber = paas.july$Fwnumber
rain.july.bestdata$Year = paas.july$Year
#Confirm that above model summary is correct by running model outside of climwin with climate data it found. 
rain.july.m = glmer(yvar~scale(climate) + scale(Age.int) + (1|Fwnumber),data=rain.july.bestdata,family="binomial")
summary(rain.july.m)

##Plot effects of the best model
library(effects)
library(MASS)
plot(Effect(focal.predictors = "climate",rain.july.m))
detach("package:MASS") #interferes with dplyr

##Get percentages
rain.july.p = rain.july.bestdata %>% group_by(Year) %>% mutate(count = n()) %>% group_by(Year,yvar,count) %>% summarise(p=n()) %>% 
  mutate(per=(p/count)*100) %>% filter(yvar=="1")
rain.july.bestdata.short = rain.july.bestdata %>% distinct(Year,.keep_all = T) %>% select(Year,climate)
rain.july.p = merge(rain.july.p,rain.july.bestdata.short,by="Year")

##Plot 
ggplot(data=rain.july.p,aes(x=climate,y=per)) + geom_point(size=3) + geom_smooth(method="lm",color="black") + 
  theme_cowplot() + xlim(0,60) + ylim(0,100) +
  xlab("Total rainfall (cm) from Feb 24th - May 3rd") + ylab("Percentage of males age 2 and older\\nin ornamented plumage by July 1st")

##Window is days 55-123 - February 24th - May 3rd

##LRT for rainfall variable
rain.july.m.norain = glmer(yvar~scale(Age.int) + (1|Fwnumber),data=rain.july.bestdata,family="binomial")
anova(rain.july.m,rain.july.m.norain) 

##Paper: X2=107.25, p<0.001

##Save above plot for plotting next to 2yo rainfall hazard ratio plot
A = ggplot(data=rain.july.p,aes(x=climate,y=per)) + geom_point(size=3) + geom_smooth(method="lm",color="black") + 
  theme_cowplot() + xlim(0,60) + ylim(0,100) +
  xlab("Total rainfall (cm) from Feb 24th - May 3rd") + ylab("Percentage of males age 2 and older\\nin ornamented plumage by July 1st")
#saveRDS(A,file=here::here("Plots","bright_by_july.rds"))





####Randomization test

##Randomizing the bdate keeps age and brightbyjuly together. 

##Randwin function didn't work here, not sure why. Here's custom function for randomization
randlist = seq(1,100,by=1)

rand.july.func = function() {
  paas.july.rand = paas.july
  output = data.frame(matrix(ncol=1,nrow=length(randlist)))
  colnames(output)="deltaAICc"
  
  for (i in 1:length(randlist)) {
    #Randomize dates across years
    paas.july.rand$date.rand = sample(paas.july.rand$date,replace = F)
    #Run climate window with random dates
    test.rand = slidingwin(xvar=list(Rain = rain$amountcm),
                           cdate = rain$Date,
                           bdate = paas.july.rand$date.rand, #This will get the correct year
                           baseline = glm(brightbyjuly~Age.int,family="binomial",data=paas.july.rand),
                           type = "absolute",
                           cinterval = "day",
                           range = c(200,0),
                           refday = c(1,7),
                           exclude = c(20,-1),
                           stat = "sum",
                           fun = "lin")
    #Get top AIC value
    rand = test.rand[[1]]$Dataset[1,1]
    #Add top AIC value to output
    output$deltaAICc[i] = rand
  }
  return(output)
}

##Run randomization
#rand.july = rand.july.func()

##Save randomization output 
#write.csv(rand.july,here::here("Output files","Molt_by_july_randomization_200_1_lin.csv"),row.names=F)

##Read in randomization
rand.july = read.csv(here::here("Output files","Molt_by_july_randomization_200_1_lin.csv"))

##Get p-value
pvalue(dataset=rain.july.dataset,datasetrand=rand.july,metric="AIC",sample.size=5)

##Plot randomizations
ggplot(data=rand.july,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=rain.july.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")


