####Pairing date and 2-year-old male molt timing


#Going to try using 2017 data for asking how pairing date relates to when a 2yo male molts. 2017 was a dry year where all the males molted
#late, so we have molt data on many 2yo males compared to other years. 2017 was also the year we did the removal experiment, so will 
#need to take that into account or maybe even include that in analysis steps


####Load data####

library(here)
library(dplyr)
library(reshape2)
library(lubridate)
library(ggplot2)
library(cowplot)
library(viridis)
library(asnipe)


##Load bird17 files - non-breeding social interaction data - comes from non-breeding social structure paper
#Includes control and removal communities
bird17 = read.csv(here::here("Input files","bird17R.csv"))
bird17$Date = as.Date(bird17$Date,"%m/%d/%y") #Convert date 
bird17 = bird17 %>% arrange(Date) #Order by date

##Load individual file - need to know when males turned bright if they did
ind17 = read.csv(here::here("Input files","Individuals2017 7_13_20.csv"))
ind17$Molt = as.numeric(as.character(ind17$Molt))
ind17$Date = as.Date(ind17$Date,"%m/%d/%y")

#Agesex file
agesex17 = read.csv(here::here("Input files","Ages and Sexes 2017.csv"))

##Load in breeding groups for upcoming breeding seasons
b_groups17 = read.csv(here::here("Input files","Breeding Groups 2017_2017colors_short.csv"),stringsAsFactors = F)
b_groups17$Date.Created = as.Date(b_groups17$Date.Created,"%m/%d/%y")

##Load in all breeding groups for groups from previous breeding season
groups = read.csv(here::here("Input files","All RBFW Groups.csv"))
groups$Date.Created = as.Date(groups$Date.Created,"%m/%d/%y")

##Load in capture data to make sure females were all banded before observations
captures = read.csv(here::here("Input files","All RBFW Captures.csv"))
captures$Date = as.Date(captures$Date,"%m/%d/%y")

##Load list of birds in removal communities in 2017 season and their FWnumbers
removals = read.csv(here::here("Input files","RemovalorControl 2017.csv"))

##Load plumages list to get pluamge in previous year
plumages = read.csv(here::here("Input files","All RBFW ind histories plumages.csv"))



####Look at relationship between times seen and degree####

#bird17 includes both kerfuffle and non-keruffle sampling points

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
birdlist17 = merge(birdlist17,agesex17.s,by="Bird")

##Sighting frequency
sightfreq17 = data.frame(colSums(gbi17))
sightfreq17$Bird = rownames(sightfreq17)
colnames(sightfreq17)[1] = "SightFreq"
birdlist17 = merge(birdlist17,sightfreq17,by="Bird")
library(igraph)
net17 = graph.adjacency(network17, mode="undirected", diag=FALSE, weighted=TRUE)
deg17 = degree(net17)
birdlist17 = data.frame(birdlist17,deg17)
detach("package:igraph")

##Determine if male of interest - 2yo male
birdlist17$moi = NA
for (i in 1:nrow(birdlist17)) {if(birdlist17$Current.Age[i]==2 & birdlist17$Sex[i]=="M" & birdlist17$Ageexact.[i]=="exact") 
  {birdlist17$moi[i]="Yes"} else {birdlist17$moi[i]="No"}}

##Plot degree by sighting frequency and color by agesex
ggplot(data=birdlist17,aes(x=SightFreq,y=deg17)) + geom_point(aes(color=moi),size=2) + theme_cowplot() + 
  scale_color_manual(values=c("black","red"))
#Definitely should drop any birds seen below 20 times, question is did we see BBL (54 times) enough? 
#Maybe? Overall degree doesn't matter much, more want to know we saw him enough that we saw him with his female
#if he ever was with her. Leave in for now and keep an eye on him. 

##Get list of 2yo males of interest
birdlist.moi = birdlist17 %>% filter(moi=="Yes",SightFreq>20)



####Look at molt trajectories for 2yo males####
ind.moi = ind17 %>% filter(Bird %in% birdlist.moi$Bird)
ind.moi$jdate = yday(ind.moi$Date)

##Plot all males
ggplot(data=ind.moi,aes(x=jdate,y=Molt)) + geom_line(aes(group=Bird,color=Bird),size=1) + theme_cowplot() +
  scale_color_viridis(discrete = T) + geom_point(aes(color=Bird))

##Plot individual plots
ggplot(data=ind.moi,aes(x=jdate,y=Molt)) + geom_line(aes(group=Bird,color=Bird),size=1) + theme_cowplot() +
  scale_color_viridis(discrete = T) + geom_point(aes(color=Bird)) + facet_wrap(facets="Bird")

##Clearly two males disappeared before making it to intermediate/ornamented plumage. Take out VRG and WLB. 
birdlist.moi = birdlist.moi %>% filter(Bird!="VRG",Bird!="WLB")
##Also VLL had quite a large gap between 0 and 100 - 40 days, take him out. All other males' gaps are less than 20 days
birdlist.moi = birdlist.moi %>% filter(Bird!="VLL")
ind.moi = ind.moi %>% filter(Bird %in% birdlist.moi$Bird)

#Plot revised individuals
ggplot(data=ind.moi,aes(x=jdate,y=Molt)) + geom_line(aes(group=Bird),size=1) + theme_cowplot() +
  scale_color_viridis(discrete = T) + geom_point(aes(color=Bird)) + facet_wrap(facets="Bird") +
  geom_hline(aes(yintercept=66),size=0.5,color="gray",lty=2) + geom_hline(aes(yintercept=33),size=0.5,color="gray",lty=2) +
  theme(legend.position = "")
#RHV only made it to 40% but he did make it to ornamented plumage before disappearing





####Get breeding group/female data for males of interest####

##Get upcoming breeding group
birdlist.moi$Bgroup = NA
birdlist.moi$Female = NA
for (i in 1:nrow(birdlist.moi)) {
  b_groups17.moi = b_groups17 %>% filter(Male==birdlist.moi$Bird[i]) %>% arrange(Date.Created)
  birdlist.moi$Bgroup[i] = b_groups17.moi$Group.Number[1]
  birdlist.moi$Female[i] = b_groups17.moi$Female[1]
}
#All paired except for RHV who disappeared before breeding began


##Get female's FWnumber
birdlist17.s = birdlist17 %>% select(Bird,FWNo) %>% rename(Female.fwno = FWNo)
birdlist.moi = merge(birdlist.moi,birdlist17.s,by.x="Female",by.y="Bird",all.x=T)


##Determine if males were paired in the previous season
birdlist.moi$Paired.prev = NA
birdlist.moi$SameF = NA
for(i in 1:nrow(birdlist.moi)) {
  groups.moi = groups %>% filter(Male.fwno==birdlist.moi$FWNo[i]) %>% filter(Year==2016) %>% arrange(desc(Date.Created))
  if(nrow(groups.moi)>0) {birdlist.moi$Paired.prev[i] = "Yes"
  if(is.na(birdlist.moi$Female.fwno[i])) {birdlist.moi$SameF[i] = "No"} else {
    if(groups.moi$Female.fwno[1]==birdlist.moi$Female.fwno[i]) {birdlist.moi$SameF[i] = "Yes"} else {birdlist.moi$SameF[i] = "No"}}
  } else {birdlist.moi$Paired.prev[i] = "No"}
}
#RHV was paired in the previous year, but doesn't look like that female was active in 2017 when he was a 2yo


##For RHV, see if he was well connected to any females
rhv = data.frame(network17["RHV",]) #Get RHV column
colnames(rhv) = "SRI"
ggplot(rhv,aes(x=SRI)) + geom_histogram() + theme_cowplot() #One bird stands out - VBV 
#VBV is a 2yo female that was paired last year with an older male who was not a breeder in 2017. So looks like her and RHV paired up

#Make sure VBV was not more connected to a different male
vbv = data.frame(network17["VBV",]) #Get VBV column - #RHV was top by far
colnames(vbv) = "SRI"
ggplot(vbv,aes(x=SRI)) + geom_histogram() + theme_cowplot() 

#Add that info into birdlist
birdlist.rhv = birdlist.moi %>% filter(Bird=="RHV")
birdlist.rhv$Female = "VBV"
birdlist.rhv$Female.fwno = 7778
birdlist.moi = birdlist.moi %>% filter(Bird!="RHV")
birdlist.moi = rbind(birdlist.moi,birdlist.rhv)





####Check relatedness of male-female pairings####
rel = read.csv(here::here("Input files","RelatednessEstimates_8xdenovo80_14to19.csv"))
rel2 = rel %>% select(Pair,ind2,ind1,wang) #Get columns in different order 
colnames(rel2)[2:3] = c("ind1","ind2")
rel = rbind(rel,rel2) #Combine so all individuals in both columns

#Get relatedness values
birdlist.moi$MFrel = NA
for (i in 1:nrow(birdlist.moi)) {
  rel.sub = rel %>% filter(ind1==birdlist.moi$FWNo[i] & ind2==birdlist.moi$Female.fwno[i]) %>% select(wang)
  if(nrow(rel.sub)>0) {birdlist.moi$MFrel[i] = rel.sub[[1]]} else {birdlist.moi$MFrel[i]=NA}
}

#Plot
ggplot(birdlist.moi,aes(x=MFrel)) + geom_histogram(binwidth = 0.02) + theme_cowplot() + xlim(-0.2,1)
#Looks like mostly all unrelated, could have one distant sibling pairing. 




####Make sure all females were banded during the non-breeding season/before or during social observations####

##First see if females were in a breeding group in the previous year
females = birdlist.moi %>% select(Female,Female.fwno)
females$Bred.prev = NA
for (i in 1:nrow(females)) {
  groups.f = groups %>% filter(Female.fwno==females$Female.fwno[i],Year=="2016")
  if(nrow(groups.f)>0) {females$Bred.prev[i]="Yes"} else {females$Bred.prev[i]="No"}
}
#Only three females were not known breeders in 2016. Two were color banded at beginning of 2017 season in early June - 6/11 and 6/12
#Other was listed as a helper in 2016 breeding season. 
#So all good on females being banded before social observations 



####Determine whether a male was ornamented previously or not####

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

#Get ornamented.prev variable
birdlist.moi$ornamented.prev = NA
for (i in 1:nrow(birdlist.moi)) {
  plum.m = plumages.pyd %>% filter(Fwnumber==birdlist.moi$FWNo[i]) %>% filter(Year=="2016") 
  if(plum.m$Plumage=="Interm" | plum.m$Plumage=="Bright") {birdlist.moi$ornamented.prev[i]="Yes"} else {birdlist.moi$ornamented.prev[i]="No"}
}
birdlist.moi$ornamented.prev = as.factor(birdlist.moi$ornamented.prev)



####Calculate molt dates####

##Get first date seen in ornmanented plumage for each male
ind.moi$week = week(ind.moi$Date)
birdlist.moi$molt.date = NA
birdlist.moi$molt.week = NA
for(i in 1:nrow(birdlist.moi)) {
  ind.moi.m = ind.moi %>% filter(Bird==as.character(birdlist.moi$Bird[i])) %>% arrange(jdate) %>% filter(Molt>=33)
  birdlist.moi$molt.date[i] = ind.moi.m$jdate[1]
  birdlist.moi$molt.week[i] = ind.moi.m$week[1]
}



####Plot weighted degree among pairs by week####

#Get week data in bird files
bird17$week = week(bird17$Date)

#For loop for plots
wdeg.weeks = data.frame(week=NA,bird=NA,wdeg=NA,seen=NA) #Create empty dataframe to rbind to in for loop
birdlist.moi = birdlist.moi %>% arrange(as.character(Bird))

par(mfrow=c(4,4))
for (i in 1:nrow(birdlist.moi)) {
  male = birdlist.moi$Bird[i]
  female = birdlist.moi$Female[i]
  bird = bird17
  weeks = data.frame(seq(25,35))
  colnames(weeks) = "week"
  weeks$bird = birdlist.moi$Bird[i]
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
  plot(weeks.print$week,weeks.print$wdeg,ylim=c(0,1),xlim=c(24,48),pch=19,main=paste(male,birdlist.moi$molt.week[i]),
       xlab="Week of year",ylab="Wdegree to female")
  lines(weeks.print$week,weeks.print$wdeg)
  abline(v=birdlist.moi$molt.week[i],col="red")
  wdeg.weeks = rbind(wdeg.weeks,weeks)
}
par(mfrow=c(1,1))

#Remove the first row from wdeg.weeks - it was blank from setup
wdeg.weeks = wdeg.weeks[-1,]

###Look at what time of year I have enough social data to compare among weeks
wdeg.weeks.noNA = wdeg.weeks %>% filter(!is.na(wdeg)) %>% filter(seen>=3)
ggplot(data=wdeg.weeks.noNA,aes(x=week,y=bird)) + geom_point() + geom_line() 

#Almost all males have enough data to use weeks 27-33 and used linear interpolation to fill in missing weeks
#Most males pretty connected to their females early on, so not expecting much in terms of when they're connected to their female
#but after looking at that, can look at connections to other birds - subtract wdeg to female from total wdeg. Or maybe wdeg to 
#group members from total wdeg. 




####For sliding window looking at pairing date, linear model or survival analysis make more sense?####

##Add in removal/control data
removals = removals %>% filter(TreatmentGroup!="X") #Remove X's and then drop levels - otherwise survival model thinks X is a factor level
removals = droplevels(removals) #No 2yo males from community X
birdlist.moi = merge(birdlist.moi,removals,by="Bird")

##Use week 30 as a test - just need an arbitrary value for now
wdeg.weeks.test = wdeg.weeks %>% filter(week==30)
birdlist.moi.test = merge(birdlist.moi,wdeg.weeks.test,by.x="Bird",by.y="bird")
birdlist.moi.test$Paired.prev = as.factor(birdlist.moi.test$Paired.prev)

##Try a linear model looking at relationship between wdeg to female at week 30 and molt date
md.test = lm(molt.date~wdeg + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi.test)
summary(md.test)
hist(resid(md.test))
library(car)
qqPlot(md.test)
plot(molt.date~wdeg,data=birdlist.moi.test,pch=19)
#Linear model works fine

##Try a survival analysis looking at relationship between wdeg to female at week 30 and molt date
library(survival)
library(survminer)
birdlist.moi.test$Status = 1
md.test2 = coxph(Surv(molt.date,Status)~wdeg + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi.test)
summary(md.test2)
plot(survfit(md.test2))
AIC(md.test2)

##Looks like either linear model or survival analysis would work. If trying to predict the occurence of an event then 
#survival analysis probably better. 





####Check assumptions of baseline survival analysis model####

##Get the baseline model 
sm = coxph(Surv(molt.date,Status)~Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi.test)
summary(sm)

#Does ornamented.prev add anything to this model? 
sm.no.op = coxph(Surv(molt.date,Status)~Paired.prev + TreatmentGroup,data=birdlist.moi.test)
AIC(sm.no.op) - AIC(sm) #Leave in

##Test proportional hazards assumption
test.ph = cox.zph(sm)
test.ph
ggcoxzph(test.ph)
#Looks good

##Test influential observations
ggcoxdiagnostics(sm, type = "deviance", linear.predictions = TRUE)
#ok, not great






####Sliding window approach for pairing and molt#### 

#Want to see if wdeg connections to female at certain dates influence whether a 1yo male molts into nuptial plumage

##Interpolate missing data when individuals missing points
wdeg.weeks.int = data.frame(week=NA,bird=NA,wdeg=NA,wdeg.int=NA)
for (i in 1:nrow(birdlist.moi)) {
  wdeg.weeks.ind = wdeg.weeks %>% filter(bird==as.character(birdlist.moi$Bird[i])) %>% select(-seen)
  wdeg.weeks.ind = wdeg.weeks.ind %>% filter(week>=27,week<=33)
  wdeg.weeks.ind$wdeg.int = NA
  if (nrow(wdeg.weeks.ind %>% filter(is.na(wdeg)))==0) {wdeg.weeks.int = rbind(wdeg.weeks.int,wdeg.weeks.ind)} else {
    wdeg.weeks.ind.int = data.frame(approx(x=wdeg.weeks.ind$week,y=wdeg.weeks.ind$wdeg,n=nrow(wdeg.weeks.ind)))
    wdeg.weeks.ind$wdeg.int = wdeg.weeks.ind.int$y
    wdeg.weeks.int = rbind(wdeg.weeks.int,wdeg.weeks.ind)
  }
}
wdeg.weeks.int = wdeg.weeks.int[-1,]

##For birds without missing data, move values to int column
for (i in 1:nrow(wdeg.weeks.int)) {if(is.na(wdeg.weeks.int$wdeg.int[i])) {wdeg.weeks.int$wdeg.int[i]=wdeg.weeks.int$wdeg[i]}}




####Function for sliding windows
pairing_window = function(range.start,range.end,wdeg.weeks.int,birdlist.moi) {
  
  #Get all possible week windows
  windows = expand.grid(start=seq(range.start,range.end),end=seq(range.start,range.end))
  windows$diff = windows$start-windows$end #Get difference between windows
  windows = windows %>% filter(diff<=0) %>% select(-diff)
  
  #Get aic value of base model
  birdlist.moi$Status = 1
  birdlist.moi$Paired.prev = as.factor(birdlist.moi$Paired.prev)
  aic.base = AIC(coxph(Surv(molt.date,Status)~Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi))
  
  #Make output dataframe
  output = data.frame(matrix(ncol=6,nrow=nrow(windows)))
  colnames(output) = c("start.week","end.week","ModelBeta","deltaAIC","Range.start","Range.end")
  
  #For each possible window
  for (i in 1:nrow(windows)) {
    #Get a mean wdeg variable for each individual in birdlist.moiw
    birdlist.moi$mean.wdeg = NA
    #For each row in birdlist.moi
    for (h in 1:nrow(birdlist.moi)) {
      #Get wdeg data for that male
      wdeg.weeks.intm = wdeg.weeks.int %>% filter(bird==birdlist.moi$Bird[h]) 
      #Filter for weeks in window
      wdeg.weeks.intm = wdeg.weeks.intm %>% filter(week>=windows$start[i],week<=windows$end[i])
      #Put the mean value in birdlist.moiw
      birdlist.moi$mean.wdeg[h]=mean(wdeg.weeks.intm$wdeg.int)
    }
    
    #Get the model with wdeg included
    model.wdeg = coxph(Surv(molt.date,Status)~mean.wdeg + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi)
    
    #Get AIC
    aic.wdeg = AIC(model.wdeg)
    
    #Add values to output dataframe
    output$start.week[i] = windows$start[i]
    output$end.week[i] = windows$end[i]
    output$ModelBeta[i] = coef(model.wdeg)[1]
    output$deltaAIC[i] = aic.wdeg - aic.base
    output$Range.start[i] = range.start
    output$Range.end[i] = range.end
  }
  
  #Order dataframe by lowest AIC values compared to baseline
  output = output %>% arrange(deltaAIC)
  
  #Output dataframes
  return(output)
}

###Run the function
pairing.output.27.33 = pairing_window(range.start = 27,range.end = 33,wdeg.weeks.int = wdeg.weeks.int,birdlist.moi = birdlist.moi)

#Plot delta plot
library(climwin)
pairing.output.27.33.climwin = pairing.output.27.33
colnames(pairing.output.27.33.climwin) = c("WindowOpen","WindowClose","ModelBeta","deltaAICc")
plotdelta(dataset=pairing.output.27.33.climwin)

##Result: The timing of when a 2-year-old male became associated with his female did not influence timing of molt for 2-year-olds 
#in the 2017 season. No windows improved upon the baseline model. All males were pretty connected to their females during the 
#non-breeding season. Get same results when using linear model instead of survival analysis




####Social connections and timing of molt####

##First export basic list to know who to calculate social connections for
bmoi.export = birdlist.moi %>% select(Bird,FWNo,Sex,Current.Age,Ageexact.,Female,Female.fwno,molt.date)
#write.csv(bmoi.export,here::here("Output files","2yo_male_list_2017.csv"),row.names = F)

#Get each male's connections to older males, to unrelated females, etc. before molt - do in separate script? 
#Load those into survival analysis and see if they influence anything. 

##Load in social connections
bmoi.social = read.csv(here::here("Output files","twoyo_social_connections.csv"))
bmoi.social = bmoi.social %>% select(Bird,wdeg_12yo_males,wdeg_old_males,wdeg_females,wdeg_female_paired,wdeg_all,wdeg_females_unrelated)
birdlist.moi = merge(birdlist.moi,bmoi.social,by="Bird")
birdlist.moi$Status = 1


####Wdeg_all before molt survival analysis
sm.wdeg = coxph(Surv(molt.date,Status)~wdeg_all + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi)
summary(sm.wdeg)
sm.wdeg.coef = coef(sm.wdeg)[1]

##Test proportional hazards assumption
test.ph1 = cox.zph(sm.wdeg)
test.ph1
ggcoxzph(test.ph1)
#Looks good

##Test influential observations
ggcoxdiagnostics(sm.wdeg, type = "deviance", linear.predictions = TRUE)
#Could be better but using permutations

##Test linearirty of continuous variable - http://www.sthda.com/english/wiki/cox-model-assumptions
ggcoxfunctional(Surv(molt.date, Status) ~ wdeg_all + log(wdeg_all) + sqrt(wdeg_all), data = birdlist.moi)
#Not very linear but transformations do not help. Sample size is quite small.

##Permutation test for significance of wdeg variable
sm.wdeg.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    birdlist.moi$wdeg_all.rand = sample(birdlist.moi$wdeg_all,replace=F)
    rands[i] = coef(coxph(Surv(molt.date,Status)~wdeg_all.rand + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = 1-sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#A = sm.wdeg.rand(perm.number = 1000,coef.obs = sm.wdeg.coef)

##Result: total weighted degree not a predictor of when a 2yo male molted (coef=-0.27,p=0.287)




####Wdeg_12yo_males before molt survival analysis
sm.wdeg12yo = coxph(Surv(molt.date,Status)~wdeg_12yo_males + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi)
summary(sm.wdeg12yo)
sm.wdeg12yo.coef = coef(sm.wdeg12yo)[1]

##Test proportional hazards assumption
test.ph1 = cox.zph(sm.wdeg12yo)
test.ph1
ggcoxzph(test.ph1)
#Looks good

##Test influential observations
ggcoxdiagnostics(sm.wdeg12yo, type = "deviance", linear.predictions = TRUE)
#Looks mostly good

##Test linearirty of continuous variable - http://www.sthda.com/english/wiki/cox-model-assumptions
ggcoxfunctional(Surv(molt.date, Status) ~ wdeg_12yo_males + log(wdeg_12yo_males) + sqrt(wdeg_12yo_males), data = birdlist.moi)
#Not very linear but transformations do not help. Sample size is quite small.

##Permutation test for significance of wdeg variable
sm.wdeg12yo.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    birdlist.moi$wdeg_12yo_males.rand = sample(birdlist.moi$wdeg_12yo_males,replace=F)
    rands[i] = coef(coxph(Surv(molt.date,Status)~wdeg_12yo_males.rand + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = 1-sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#B = sm.wdeg12yo.rand(perm.number = 1000,coef.obs = sm.wdeg12yo.coef)

##Result: Weighted degree to other young males (1 and 2 year olds) was not a predictor of when a 2yo male molted (coef=-0.49,p=0.325)






####Wdeg_old_males before molt survival analysis
sm.wdegoldm = coxph(Surv(molt.date,Status)~wdeg_old_males + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi)
summary(sm.wdegoldm)
sm.wdegoldm.coef = coef(sm.wdegoldm)[1]

##Test proportional hazards assumption
test.ph1 = cox.zph(sm.wdegoldm)
test.ph1
ggcoxzph(test.ph1)
#Looks good

##Test influential observations
ggcoxdiagnostics(sm.wdegoldm, type = "deviance", linear.predictions = TRUE)
#Could be better but using permutations so shouldn't matter too much

##Test linearirty of continuous variable - http://www.sthda.com/english/wiki/cox-model-assumptions
ggcoxfunctional(Surv(molt.date, Status) ~ wdeg_old_males + log(wdeg_old_males) + sqrt(wdeg_old_males), data = birdlist.moi)
#Not very linear but transformations do not help. Sample size is quite small.

##Permutation test for significance of wdeg variable
sm.wdegoldm.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    birdlist.moi$wdeg_old_males.rand = sample(birdlist.moi$wdeg_old_males,replace=F)
    rands[i] = coef(coxph(Surv(molt.date,Status)~wdeg_old_males.rand + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi))[1]
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
#C = sm.wdegoldm.rand(perm.number = 1000,coef.obs = sm.wdegoldm.coef)

##Result: Weighted degree to older males was not a predictor of when a 2yo male molted (coef=0.29,p=0.4)





####Wdeg_females before molt survival analysis
sm.wdegf = coxph(Surv(molt.date,Status)~wdeg_females + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi)
summary(sm.wdegf)
sm.wdegf.coef = coef(sm.wdegf)[1]

##Test proportional hazards assumption
test.ph1 = cox.zph(sm.wdegf)
test.ph1
ggcoxzph(test.ph1)
#Looks good

##Test influential observations
ggcoxdiagnostics(sm.wdegf, type = "deviance", linear.predictions = TRUE)
#Could be better but using permutations so shouldn't matter too much

##Test linearirty of continuous variable - http://www.sthda.com/english/wiki/cox-model-assumptions
ggcoxfunctional(Surv(molt.date, Status) ~ wdeg_females + log(wdeg_females) + sqrt(wdeg_females), data = birdlist.moi)
#Not very linear but transformations do not help. Sample size is quite small.

##Permutation test for significance of wdeg variable
sm.wdegf.rand = function(perm.number,coef.obs) {
  rands = matrix(ncol=1,nrow=perm.number)
  for (i in 1:perm.number) { 
    birdlist.moi$wdeg_females.rand = sample(birdlist.moi$wdeg_females,replace=F)
    rands[i] = coef(coxph(Surv(molt.date,Status)~wdeg_females.rand + Paired.prev + TreatmentGroup + ornamented.prev,data=birdlist.moi))[1]
  }
  rands = data.frame(rands)
  print(ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
          geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency"))
  P = 1-sum(coef.obs < rands$rands)/perm.number
  print(P)
  coef.rand = coef.obs - mean(rands$rands)
  print(coef.rand)
  graph = ggplot(data=rands,aes(x=rands)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = coef.obs,size=1,color="red") + theme_cowplot() + xlab("Coefficients") + ylab("Frequency")
  return(graph)
}

##Run permutation
#D = sm.wdegf.rand(perm.number = 1000,coef.obs = sm.wdegf.coef)

##Result: Weighted degree to females was not a predictor of when a 2yo male molted (coef=-0.26,p=0.345)


##Plot all four randomizations above
ggarrange(A,B,C,D,labels=c("A","B","C","D"))


