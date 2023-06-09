#LOAD DATA----
Choice<-read.csv("....2-Choice-Rdata.csv")

#2-CHOICE RESULTS----
#Is small or large chosen first often in 2-choice experiment? 
#_____________________________________________________________
Firstchoice<-table(Choice$X1stRaid)
Firstchoice
ct_first<-chisq.test(c(8,8)) #and no video data on 3 trials
ct_first #X-squared = 0, df = 1, p-value = 1
#They are chosen equally


# Was the first nest found the first to be raided?
#_________________________________________________
NoChoice<-Choice[Choice$X1stRaid==Choice$X1stFind_Nest,]
NoChoice # 4Large and 7Small were raided when they were the first nest found

#How many times were both nests raided? 
#______________________________________
table(Choice$RaidSumm)
# both large small 
# 12     4     3
#63.17% of the time

#Percent of the time raid the first thing you find
#in 11 trials, the first nest found was the first to be raided
#in 4 trails, there is no data
#there are 19 trials total
# 19 - 4 = 15 viable trials
11/15 #0.7333 trials raided the first nest they found

#How often were both nests found before raiding? 
#______________________________________________
FoundBefore<-table(Choice$BothFoundB4)
FoundBefore
#When both found, 3 chose Large, 1 chose small
#4/14 trials, data missing from 5 trials

#LOGISTIC REGRESSION----
#____________________
library(lme4)
#Does the outcome depend on colony size? Nope
Choice$Logic_1stRaid<-as.logical(Choice$X1stRaid=="large")#Conver to logical variable
ChoiceLog<-Choice[-c(4,7,10,14,20,21),] #Remove empty cells
logistic_1stRaid<-lm(Logic_1stRaid~Pw,family=binomial(link = "logit"),data=ChoiceLog)
summary(logistic_1stRaid)

#Does time spent searching depend on colony size? Nope
lm_SearchTime<-lm(SearchTime_1stFind~Pw,na.action="na.omit",data=Choice)
lm_SearchTime
summary(lm_SearchTime)

#CALCULATING ENCOUTNER RATES----
#Mean Search time = 7476 sec = 2h 4m 36s
#_______________________________________
summary(Choice$SearchTime_1stFind)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  861    3088    3594    7476   12330   23620       6 
n<-length(Choice$SearchTime_1stFind[!is.na(Choice$SearchTime_1stFind)])
sd<-sd(Choice$SearchTime_1stFind,na.rm = TRUE)
se<-sd(Choice$SearchTime_1stFind,na.rm = TRUE)/sqrt(length(Choice$SearchTime_1stFind[!is.na(Choice$SearchTime_1stFind)])) 
#SD = 6780.13 = 1h 53m 0s
#SE = 1880.47 - 31m 20s

error<-qt(0.975,df=n-1)*sd/sqrt(n)
#error = 4097.191

#Confidence Intervals
upci<-7476+error # upper CI = 11573.19
loci<-7476-error # lower CI = 3378.809

library(lubridate)
seconds_to_period(7476)
seconds_to_period(6780)
seconds_to_period(1880)

#LATENCY TO RAID----
hist(Choice$Latency)
x<-Choice[Choice$X1stRaid=="small",]
y<-Choice[Choice$X1stRaid=="large",]
wilcox.test(x$Latency, y$Latency, alternative = "greater")
