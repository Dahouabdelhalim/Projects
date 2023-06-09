#Analysis file for "Decoupling cooperation and punishment in humans shows that punishment is not an altruistic trait." by Burton-Chellew, M.N. & Guerin, C. 2021

#note the experiment/project had an official codename of "PIN".


# Step 1: 
print("Step 1: import the z-Tree data")
# using package zTree (from Kirchkamp_2019_JBEF) to upload all the data

if(!require(zTree)){
  install.packages("zTree")
  library(zTree)
}

#uploading data according to their instructions
(myFilenames <- list.files(pattern = "*.xls"))
allData <- zTreeTables(myFilenames)
summary(allData)

#Here is the dimensions and a subset of the globals tables:
dim(allData[["globals"]])


#40 variables, here is a description of them all
# Date = the date and time of the session
# Treatment = does not mean experimental treatment, but means z-Tree file order in experiment (1 = instructions, 2 = strategy method stage, 3 = public good game)
# Period = period/round of that stage ('treatment')
# NumPeriods = how many that stage had (1, 1, or 5 for the repeated pgg)
# RepeatTreatment = ignore
# ExchangeRate = the exchange rate we used between Monetary Units and Swiss CHF, note this was reduced from 0.05 to 0.04 after session 3
# Ascending/Descending = coding of the words for coding purposes in the subjects table
# Random1-3 = random number from 0-20 to simulate computer players in strategy method
# Computer1-3,..Sum,...Average,...AverageInteger = integer version of Random1-3 to represent computer players and summaries
# AuctionStop / noStop = ignore
# N = group size for public good (4)
# Endowment = endowment for public good (20)
# MPCR = marginal per capita return for public good (0.4)
# EfficientyFactor = N*MPCR (1.6)
# PunishmentFactor = how much deductions are multiplied, so punisher spends 1 MU to deduct 3 MU from punishee
# PunishmentLimit = the maximum deduction a player could spend on another individual (spend 6 to deduct 18), NOTE: THIS CHANGED FROM 6 TO 10 after session 2
# PunishmentBudget = The fresh endowment given to players to enable them to punish others (3*limit because 3 other group members), NOTE: THIS CHANGED FROM 18 TO 30 AFTER SESSION 2, this was so a punisher could feasibly equalize payoffs between themself and a free rider who had 20 MU more than them (each unit spent on punishment reduces the inequality by 2)
# Yes/ No = coding of the words for subjects coding purposes
# Tyrant / Diplomat / Solider / Peasant = coding of those words numerically as 1-4 for subjects table 
# AllSoldiers, TyrantvSoldiers, TyrantvPeasants, DiplomatvSoldiers = coding of those names for 'Worlds'
# World = the World (aka 'Scenario') that session modeled 
#(each session was either AllSoldiders(1); TyrantvSoldiers(2); or TyrantvPeasants(3), and DiplomatvSoldiers was not conducted in the end)
# 'AllSoldiers', i.e. same as Fehr & Gachter 2002 design, aka 'Mutual Punishers')
# TyrantsvSoldiers (1 immune punisher, aka "Immune Punisher') 
# TyrantsvPeasants (1 punisher, aka 'Sole Punisher')
# DiplomatvSoldiers (1 immune player who cannot punish, 3 punishers) - not conducted
# Session = Session number (1-20)
# Special = Tyrant (not important, can ignore)
# Normal = either Soldier or Peasant (not important, can ignore)

#And here is a subset of the subjects tables: 
allData[["subjects"]][100:105,1:5]
dim(allData[["subjects"]])
#$subjects:'data.frame': 2964 obs. of  205 variables:

#Importing questionnaires: z-Tree stores information from the questionnaire at the end of the experiment in a table with the extension
#.sbj.
( myQuestFiles <- list.files(pattern="*.sbj") )
## [1] "160215_0810.sbj" "160215_0949.sbj"

#Now myQuestFiles is a vector of filenames. Each name refers to the subjects from one session of the experiment (here we have only two sessions). While the .xls table stores subject specific information as one line per subject, the .sbj table uses one column per subject. The function zTreeSbj reads a vector of files and transposes them to one line per subject. We read subjects from all sessions with the following command:
allQuest <- zTreeSbj(myQuestFiles) ## reading 160215_0810.sbj ... ## reading 160215_0949.sbj ...

#The data frame allQuest contains all subjects from all experiments. Here is a subset:
allQuest[,5:6]
levels(as.factor(allQuest$Gender))
levels(as.factor(allQuest$Age))
summary(as.factor(allQuest$Gender))
summary(as.factor(allQuest$Age))

#We can easily merge data from the questionnaire with the subjects table:
subjectsJoined <- merge(allData[["subjects"]],allQuest)

#Here is a subset:
subjectsJoined[990:1000,c("Date","Subject", "Group","Profit","Age","Gender", "Period", "Contribution")]
dim(subjectsJoined)
#each subject (participant) has 7 rows of data, 1 for instructions phase, 1 for strategy method phase, and 1-5 for rPGG phase.
#I do not know why the periods do not seem to be in order


###################### Step 2: data wrangling #####################################
print("Step 2 - renaming and creating variables")
###########################################################

if(!require(data.table)){
  install.packages("data.table")
  library(data.table)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}


# joining some key variables from the Globals table into the subjects table
GlobalsData <- allData[["globals"]]
FromGlobals = GlobalsData  %>% filter(Treatment==3 & Period ==1) %>% select(Date, c(World, Session, PunishmentBudget))
SbjEdit = inner_join(subjectsJoined, FromGlobals, by = "Date", suffix = c("raw", "joined")) 
nrow(SbjEdit) == nrow(subjectsJoined) #check still the same rows but different columns(+2) as SubjectsData
class(SbjEdit)
SbjEdit <- data.table(SbjEdit)
class(SbjEdit)
dim(SbjEdit)

#overwriting endowment column to always be 20 (replacing NAs) in subjects data for glmms later on
SbjEdit[,Endowment := 20]

#creating Session Number and Character variables, can be useful
SbjEdit <- rename(SbjEdit, SessionNumber = Session)
SbjEdit$SessionCh <- as.factor((paste("S",SbjEdit$SessionNumber, sep = "")))
levels(SbjEdit$SessionCh)
class(SbjEdit$SessionCh)
SbjEdit[Period ==3 & Subject ==11, SessionCh] #check all 20 sessions entered, as every session had a Period ==3 and Subject ==1.

#arrange session character factor levels in numerical order
SbjEdit$SessionCh <- factor(SbjEdit$SessionCh, levels = c(
  "S1",
  "S2",
  "S3",
  "S4",
  "S5",
  "S6",
  "S7",
  "S8",
  "S9",
  "S10",
  "S11",
  "S12",
  "S13",
  "S14",
  "S15",
  "S16",
  "S17",
  "S18",
  "S19",
  "S20"
))
levels(SbjEdit$SessionCh)


#Grouping sessions by Scenario type rather than chronologically
SbjEdit[SessionNumber == 1, SessionType := "IP.1"]
SbjEdit[SessionNumber == 3, SessionType := "IP.2"]
SbjEdit[SessionNumber == 6, SessionType := "IP.3"]
SbjEdit[SessionNumber == 8, SessionType := "IP.4"]
SbjEdit[SessionNumber == 9, SessionType := "IP.5"]
SbjEdit[SessionNumber == 10, SessionType := "IP.6"]
SbjEdit[SessionNumber == 15, SessionType := "IP.7"]
SbjEdit[SessionNumber == 16, SessionType := "IP.8"]

SbjEdit[SessionNumber == 2, SessionType := "SP.1"]
SbjEdit[SessionNumber == 4, SessionType := "SP.2"]
SbjEdit[SessionNumber == 5, SessionType := "SP.3"]
SbjEdit[SessionNumber == 7, SessionType := "SP.4"]
SbjEdit[SessionNumber == 11, SessionType := "SP.5"]
SbjEdit[SessionNumber == 12, SessionType := "SP.6"]
SbjEdit[SessionNumber == 13, SessionType := "SP.7"]
SbjEdit[SessionNumber == 14, SessionType := "SP.8"]

SbjEdit[SessionNumber == 17, SessionType := "MP.1"]
SbjEdit[SessionNumber == 18, SessionType := "MP.2"]
SbjEdit[SessionNumber == 19, SessionType := "MP.3"]
SbjEdit[SessionNumber == 20, SessionType := "MP.4"]

SbjEdit$SessionType <- factor(SbjEdit$SessionType, levels = c(
  "MP.1",
  "MP.2",
  "MP.3",
  "MP.4",
  "IP.1",
  "IP.2",
  "IP.3",
  "IP.4",
  "IP.5",
  "IP.6",
  "IP.7",
  "IP.8",
  "SP.1",
  "SP.2",
  "SP.3",
  "SP.4",
  "SP.5",
  "SP.6",
  "SP.7",
  "SP.8"
))
levels(SbjEdit$SessionType)


#need to create Unique Participant IDentity (UPID)
SbjEdit$UPID <- (paste(SbjEdit$Date,SbjEdit$Subject, sep = "-"))
SbjEdit$UPID
uniqueN(SbjEdit$UPID) #should be 420

# create Unique Group IDentity (UGID), note unique groups only ever exist for one round
SbjEdit$UGID <- (paste(SbjEdit$Date,SbjEdit$Group,SbjEdit$Period, sep = "-"))
SbjEdit$UGID
uniqueN(SbjEdit$UGID) #should be 420/4*5 = 525

#setting relevant variables as factors and renaming as necessary
SbjEdit[World ==1, Scenario := "Mutual Punishers" ]
SbjEdit[World ==2, Scenario := "Immune Punisher" ]
SbjEdit[World ==3, Scenario := "Sole Punisher" ]
#arrange world/treatment factor levels how I want
SbjEdit$Scenario <- factor(SbjEdit$Scenario, levels = c("Mutual Punishers", "Immune Punisher", "Sole Punisher"))

SbjEdit <- rename(SbjEdit, Stage = Treatment)  #renaming 'Treatment' as 'Stage' of the experiment. Treatment is a zTree misnomer. The real treatments in this experiment are the 'World' (aka Scenario) and player subtype within each World.
nrow(SbjEdit) / 420 #this should equal 7, as 7 rows of data per participant (N=420). However, there was a problem in that in Session 3 we accidentally restarted the public good game instead of going to the questionnaire phase. We realized immediately and cancelled the game. This means there are 24 extra rows of data as the session had 24 participants, each enter 1 new round (they were all instructed to enter 0 as their contribution, the participant payments were somewhat inflated in this session).
SbjEdit <- SbjEdit[Stage != 4]  #this is to fix the above problem.
nrow(SbjEdit) /420  #as you can see, there are now 24 fewer rows, and 7 rows per participant, so the modification worked correctly.

#changing some variables from numerical to categorical
SbjEdit$Type <- as.factor(SbjEdit$Type)
SbjEdit$Immune <- as.factor(SbjEdit$Immune)
SbjEdit$Punisher <- as.factor(SbjEdit$Punisher)
dim(SbjEdit)
class(SbjEdit)

#key variables
#Contribution: this is the individual input per round to the public good (0-20 MU)
#MySumDeductionofOthers: this is the total spend that individual made in that round 'punishing' others, I change the name to 'PunishmentSpending' below
#GroupPunOfMe: this is how many MU were deducted from focal player (in multiples of 3 because 1 MU spent punishing focal player causes focal player to lose 3 MU). Renamed to 'PunishmentReceived'.

#duplicating and renaming the punishment variables to something simpler
SbjEdit[,PunishmentSpending := MySumDeductionofOthers] #How much the UPID spent on punishing in that round (multiples of x1)
SbjEdit[,Punished := ifelse(PunishmentSpending == 0, 0, 1)] #binary variable, did the UPID punish anyone that round?
SbjEdit$Punished <- as.factor(SbjEdit$Punished)
SbjEdit[,PunishmentReceived := GroupPunOfMe] #How much the UPID was punished that round (multiples of x3, so the fine received, not what it cost the punishers)
SbjEdit[,GotPunished := ifelse(PunishmentReceived == 0, 0,1)] #binary variable, was the UPID punished that round or not?
SbjEdit$GotPunished <- as.factor(SbjEdit$GotPunished)


##### Step 3: classifying Strategy Method responses ########
print("Step 3 - classifying Strategy Method")
###########################################################

# Please note. Sessions 1-3 were problematic for the strategy method, a simple solution is to ignore all strategy method responses from sessions 1-3 (although that is not strictly necessary).
# For the purposes of another paper, investigating the Strategy Method in detail, the participants played the SM with both computers and then with humans, or vice versa. They also faced either Ascending (0-20) or Descending (20-0) input boxes. These comparisons are primarily for another paper. However, due to some coding mistakes in the first 3 sessions, many responses are invalid here. Full details below.

# variable SMorder refers to whether the entry boxes were ascending (0-20, as usual) or descending (20-0), not this is not about order of humans then computers or vice versa
#SMorder ==1, players saw on their screen a contribution schedule that ran from 0 to 20 in ascending order (so top of the screen was 0, then as they moved down the screen the value increased)
# in contrast, SMorder ==0 means players saw a descending/decreasing contribution schedule.
#in session 1-3, i randomly assigned the order separately for the with human and with computer phase, upon reflection this was not a good idea, because participants maybe didn't notice the second time that the numbers were in different order (if they did not get the same order). Then they were not sure if this was a coding error or not. I was alerted to this by a participant who told me the numbers were in the 'wrong order'. This also goes to show a difficulty of overriding participant's intuition in incentivsed experiments.
# anyway, after session 3 we changed the coding so they got the same order for both, to do this we used a variable called SMorder, instead of SMwCorder (i.e. SMwithComputersOrder) and SMwHorder (although these both still appear in the code they were not of any impact)

#another problem in Sessions 1 and 3 was that some participants saw an error, in which case all the descending numbers did not go from 0 to 20, but all just said 20! Although they had been told to expect 0-20 so some still probably answered as if that was the case, but of course this data is worthless now.

# Specifically, in Session 1 (humans first) only participants that saw an ascending order in the with humans phase, and then an ascending order again with the computers phase, are valid. This is 7 of 20 participants (subjects 5,6, 9, 10, 11, 12, 18)

# Specifically, in Session 2 (humans first again), the descending 20s problem had been fixed for this session but the inconsistency between with humans and computers remained. Therefore only participants that saw ascending:ascending, or descending:descending, are valid, 12 in total (subjects 3, 4, 6, 7, 8, 10, 11, 12, 14, 16, 19, 20)

# Specifically, in Session 3 (computers first), both problems occur (because I failed to realize I had to fix the descending 20s for both session orders, human or computer first). In this session, the players played with computers first, and all of these are valid. Then they played with humans, only those that has ascending numbers and also had had ascending numbers with computers are valid for this with human data, which is only 4 participants (8, 9, 10, 24)

#making a new variable that records withe the Strategy Method with Humans (SMH) was valid or not (SMHvalid).
SMHvalid <- SbjEdit %>%
  filter(Stage ==2) %>%
  select(SessionNumber, Subject, UPID, SMwHorder, SMwCorder) %>%
  mutate(SMHvalid := if_else(SessionNumber < 3 & SMwHorder ==0, 0, if_else(SessionNumber == 3 & SMwCorder !=1 | SessionNumber ==3 & SMwHorder != 1, 0, 1))) %>%
  select(UPID, SMHvalid)
SMHvalid
sum(SMHvalid$SMHvalid) #ok this seems to work and gives the expected value of 380 valid and 40 invalid (invalids are 10 from S1, 10 from S2, 20 from S3). Need to join to SbjEdit data

#inner join of SMValid to SbjEdit
nrow(SbjEdit) / 420 #should be 7 each
SbjEdit = inner_join(SbjEdit, SMHvalid, by="UPID", suffix = c("raw", "joined"))
SbjEdit <- data.table(SbjEdit)
dim(SbjEdit)
class(SbjEdit$SMHvalid)
summary((as.factor(SbjEdit$SMHvalid)))/7 #ok this returns correct values (divided by 7 because each participant has 7 rows of data)
nrow(SbjEdit) /420 #should still be 7 each

#making a new variable that records whether the Strategy Method with Computers (SMC) was valid or not (SMCvalid).
SMCvalid <- SbjEdit %>%
  filter(Stage ==2) %>%
  select(SessionNumber, Subject, UPID, SMwHorder, SMwCorder) %>%
  mutate(SMCvalid := if_else(SessionNumber == 1 & SMwHorder ==1 & SMwCorder ==1, 1, if_else(SessionNumber == 2 & SMwCorder == SMwHorder, 1, if_else(SessionNumber > 2, 1, 0)))) %>%
  select(UPID, SMCvalid)
SMCvalid
sum(SMCvalid$SMCvalid) #ok this seems to work and gives the expected value of 399 valid and 21 invalid (13 invalid from S1, 8 invalid from S2, rest all valid). Need to join to SbjEdit data

#inner join again
nrow(SbjEdit) / 420 #should be 7 each
SbjEdit = inner_join(SbjEdit, SMCvalid, by="UPID", suffix = c("raw", "joined"))
SbjEdit <- data.table(SbjEdit)
dim(SbjEdit)
class(SbjEdit$SMCvalid)
summary((as.factor(SbjEdit$SMCvalid)))/7 #ok this returns correct values (divided by 7 because each participant has 7 rows of data)
nrow(SbjEdit) /420 #should still be 7 each


#################
print("Classifying SM (with humans and with computers in Stage 2 - before the rPGG")
SbjEdit[1,54:74] #list of the SM variables
SbjEdit[1,98:118] #same, but for the SM with computers
SMvariables <- colnames(SbjEdit[0,c(54:74, 98:118)]) #now putting all the relevant col names into a vector


if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
} #needed to change Predictor variable to numeric by dropping text

if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
} #need to gather variables



SMdataCombi <- SbjEdit %>%
  filter(Stage==2) %>% #filter SbjEdit data to stage 2 only (the Strategy method stage)
  gather(key="Predictor", value="Response", all_of(SMvariables)) %>% #if variable starts with CCC is computer, if HCC is human
  arrange(UPID) %>%
  mutate(PredictorN = str_remove(Predictor, "HCContribution")) %>%
  mutate(PredictorN = str_remove(PredictorN, "CCContribution"))
SMdataCombi$PredictorN <- as.numeric(SMdataCombi$PredictorN)
SMdataCombi <- data.table(SMdataCombi)
SMdataCombi$Groupmates <- ifelse(grepl("^HC", SMdataCombi$Predictor), "Human", "Computer")
SMdataCombi[Predictor=="HCContribution5"] #checking
SMdataCombi[Predictor=="CCContribution5"] #checking
nrow(SMdataCombi)/21/2 # should equal the number of participants, 420. 21 rows for each SM, x2 SMs
#remember to filter out the invalid responses before analyses


SMcors <- SMdataCombi %>%
  group_by(UPID, Groupmates, Stage) %>%
  summarise(Session = mean(SessionNumber),
            SMHvalid = mean(SMHvalid),
            SMCvalid = mean(SMCvalid),
            meanR = mean(Response),
            R0 = Response[PredictorN==0],
            R1 = Response[PredictorN==1],
            R2 = Response[PredictorN==2],
            R3 = Response[PredictorN==3],
            R4 = Response[PredictorN==4],
            R5 = Response[PredictorN==5],
            R6 = Response[PredictorN==6],
            R7 = Response[PredictorN==7],
            R8 = Response[PredictorN==8],
            R9 = Response[PredictorN==9],
            R10 = Response[PredictorN==10],
            R11 = Response[PredictorN==11],
            R12 = Response[PredictorN==12],
            R13 = Response[PredictorN==13],
            R14 = Response[PredictorN==14],
            R15 = Response[PredictorN==15],
            R16 = Response[PredictorN==16],
            R17 = Response[PredictorN==17],
            R18 = Response[PredictorN==18],
            R19 = Response[PredictorN==19],
            R20 = Response[PredictorN==20],
            maxR = max(Response), #don't know where this occurred
            minR = min(Response), #don't know where this occured, but if maxR == minR then flat profile
            Monot_i = ifelse(minR != maxR, all(Response == cummax(Response)),0), #this calculates if monotonically increasing
            Monot_d = ifelse(minR != maxR, all(Response == cummin(Response)),0), #this calculates if monotonically decreasing
            Pearson = cor(PredictorN,Response, method="pearson")) %>% #this calculates the Pearson correlation
  mutate( CC = ifelse(Monot_i ==1 | (Pearson > 0.5 & R20 > meanR), 1, 0), #Conditional cooperator a la Thoni, has a pearson cor > 0.5, ensuring the r20 is above the meanR eliminates strong triangle cooperators a la Thoni
          PCC = ifelse(Pearson == 1, 1, 0), #Perfect Conditional Cooperators match PredictorN perfectly, so Pearson correlation has to be 1.0 (and this is only case)
          iPCC = ifelse(CC ==1 & PCC ==0, 1, 0), #imperfect, perfect and imperfect are mutually exclusive
          NCC = ifelse(Monot_d ==1 | (Pearson < -0.5 & R0 > meanR), 1, 0), #Negative conditional cooperators have a Pearson cor < -0.5
          PNCC = ifelse(Pearson == -1, 1, 0), #Perfectly Negative Conditional cooperators inverse match (20-0, 19-1, 18-2,..), cor = -1
          iPNCC = ifelse(NCC ==1 & PNCC ==0, 1, 0), #imperfect, perfect and imperfect are mutually exclusive
          Angel = ifelse(meanR == 20, 1, 0), #Angels are Max or Perfect Cooperators give 100% (20) in all cases, so meanR has to be 20
          UC = ifelse((maxR == minR) & meanR > 0 & meanR < 20, 1, 0), #flat profile, not 0 and not 20
          FR = ifelse(meanR == 0, 1, 0),
          Other = ifelse((CC + NCC + Angel + UC + FR)==0, 1, 0),
          Overlaps = ifelse((CC + NCC + Angel + UC + FR)>1, 1,0),
          Stype = ifelse(CC ==1, "CC", ifelse(FR==1, "FR", "Other"))) #if change the CC==1 to PCC==1 then can use just perfect CCs, change line 588 also
SMcors
SMcors[is.na(SMcors)] <- 0
SMcors
class(SMcors)
SMcors <- data.table(SMcors)

print("Summary statistics from SM")
summary(as.factor(SMcors[Stage==2 & Groupmates=="Human" & SMHvalid==1 ]$Stype)) #274 of 380 valid responses were Conditional Cooperators, 40 were Free Riders, 66 other, if use all, it is 306 cc, 46 FR, 68 other
summary(as.factor(SMcors[Stage==2 & Groupmates=="Computer" & SMCvalid==1]$Stype)) #274 of 399 valid responses were Conditional Cooperators with Computers, 60 were Free Riders, 65 other, if use all, it is 65 free riders (287 CC and 68)


summary(SMcors[Stage==2 & Groupmates=="Human" & SMHvalid==1]$Pearson) #mean Pearson correlation with Humans was 0.6388, median = 0.9230
summary(SMcors[Stage==2 & Groupmates=="Computer" & SMCvalid==1]$Pearson) #mean with computers was 0.6243, median = 0.8938
summary(SMcors[Stage==2 & Groupmates=="Human" & SMHvalid==1]$meanR) #mean response (contribution) with humans across all 21 options was 6.556, median = 7.381
summary(SMcors[Stage==2 & Groupmates=="Computer" & SMCvalid==1]$meanR) #and with computers the mean was 5.915, median = 6.286

# ok now need to add these Stypes by stage and groupmate to SbjEdit
# un-gather (spread aka pivot_wider) variables to avoid multiplying dataset
SMcorskey <- select(SMcors, c(UPID, Stage, Groupmates, Stype, meanR, Pearson))
SMtypes <- pivot_wider(SMcorskey, UPID, names_from = Groupmates, values_from = c(Stype,meanR, Pearson))

dim(SbjEdit)
SbjEditSM = inner_join(SbjEdit, SMtypes, by="UPID", suffix = c("raw", "joined"))
dim(SbjEditSM)

SbjEdit <- SbjEditSM %>%
  mutate(Consistent = ifelse(Stype_Computer == Stype_Human, 1, 0),
         DoubleType = ifelse(Stype_Human == "CC" & Stype_Computer == "FR", "True CC", #i.e. Rational CC
                             ifelse(Stype_Human == "FR" & Stype_Computer == "FR", "True FR", #i.e. Homo Eoconomicus
                                    ifelse(Stype_Human == "CC" & Stype_Computer == "CC", "Confused CC","Other"))))
SbjEdit <- data.table(SbjEdit)
SbjEditdoublevalid <- SbjEdit[SMHvalid ==1 & SMCvalid==1]
summary(as.factor(SbjEditdoublevalid[Stage==2]$DoubleType)) #this is repeated for all periods/stages
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "CC" & Stype_Computer == "CC"]$DoubleType))
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "CC" & Stype_Computer == "FR"]$DoubleType))
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "CC" & Stype_Computer == "Other"]$DoubleType))
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "FR" & Stype_Computer == "CC"]$DoubleType))
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "FR" & Stype_Computer == "FR"]$DoubleType))
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "FR" & Stype_Computer == "Other"]$DoubleType))
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "Other" & Stype_Computer == "CC"]$DoubleType))
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "Other" & Stype_Computer == "FR"]$DoubleType))
summary(as.factor(SbjEditdoublevalid[Stage==2 & Stype_Human == "Other" & Stype_Computer == "Other"]$DoubleType))
#ok have 3-4 variables for SM in SbjEdit now.
#Stype_Human - a la Thoni & Volk, CC, FR, Other (triangles are other)
#Stype_Computer - same but with computerized groupmates, within participant design
#Consistent - if got same classification with humans and with computers
#DoubleType - a classification based on their combined response, Homo Economicus is double FR, true CC is CC with humans and FR with computers, CC with both is confused CC, take care not to include participants that did not get to do a valid response for BOTH with humans and compputers.
#these are probably for separate paper on Strategy Method




################# Step 4: Analyzing data, Results ####################
print("Step 4 - Analyzing Data")
print("Step 4.1 - Analyzing Contributions")
#####################################################################
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(lme4)){
  install.packages("lme4")
  library(lme4)
}

#restricting data to just rows of public good game (PGG)
SbjEditPGG <- SbjEdit[!is.na(Contribution) & Stage ==3] #this is redundant, stage 3 is PGG and Contribution should only be na in other stages, no harm.
nrow(SbjEditPGG) / 420 / 5 #nrow / nparticipants / nrounds should equal 1 if done correctly

#borrow color blind colours from icompara black box figure
#colourN12M4 = "#56B4E9"
#colourN3M8 = "#CC79A7"
#colourN3M4 = "#009E73"

#making figure of contributions over time by scenario and immunity
MyShapes = c(21,21)
DefaultColoursGreen = c("#F8766D","#7CAE00")
BlackRed = c("Black", "Red")
ColorBlindFriendly = c("#CC79A7","#56B4E9")
ColourChoice = ColorBlindFriendly

#making a categorical subgroup for plotting, session plus type of player (immune or not)
SbjEditPGG[Immune == 1, SubGroup := SessionNumber+0.1]
SbjEditPGG[Immune == 0, SubGroup := SessionNumber]

#Figure showing what happened with contributions over time
FigCon <- ggplot(SbjEditPGG, aes(Period, Contribution, fill = Immune)) +
  geom_line(aes(color = Immune, group = SubGroup), stat="smooth", method = "lm", alpha = 0.75, show.legend = FALSE)+ #adds linear model regression for each session (separately for immune and non-immune players)
  stat_summary(aes(), fun = "mean", geom = "line", size = 1.25)+ #non statistical, shows the time dynaics from round to round, visual aid only
  stat_summary(aes(), geom="linerange", fun.data=mean_cl_boot, fun.args=list(conf.int=0.95), linetype = 1, alpha = 1, color = "grey20")+ #adds 95% bootstrapped confidence intervals for each round
  #scale_linetype_manual(values = c(1,2), guide = FALSE)+ #did have different linetypes for immune not-immune but decided to use color aesthetic instaed
  coord_cartesian(xlim = c(0.5,5.5), ylim = c(-0.2,20), expand = FALSE)+
  stat_summary(fun = "mean", size=1.1, shape = 21)+ #shows a marker for the mean contribution in each round
  stat_summary(aes(colour = Immune), fun = "median", size=0.8, shape = 23, fill="white", show.legend = FALSE)+ #shows median, which is useful especially for non-normal distributions, but makes figure too crowded i think?
  scale_shape_manual(values=MyShapes, labels = c("Can be punished", "Immune"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  scale_color_manual(values = ColourChoice,labels = c("Can be punished", "Immune"))+
  labs(x="Game round (1-5)", y="Mean contribution (0-20 MU)", shape = "Player role", fill = "Player role", color = "Player role")+
  theme_classic()+
  theme(legend.position = "top")+
  facet_wrap(as.factor(SbjEditPGG$Scenario))
FigCon +  
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14, angle = 0))
#save as MS_Figure_Contributions
ggsave("MS_Figure_2_contributions.pdf",height = 185, width = 180, units = c("mm"))
#what are these warnings about missing values for? ah it's about error bars, not a problem as make own error bars for CIs.

#Figure restricting Immune data to just Conditional Cooperators
CCfigdata <- SbjEditPGG[Immune == 0 | (Immune == 1 & Stype_Human =="CC" & SMHvalid ==1)]
dim(CCfigdata)
dim(SbjEditPGG) # difference is 165 rows, divided by 5 means 33 immune participants excluded
FigConCC <- ggplot(CCfigdata, aes(Period, Contribution, fill = Immune)) +
  geom_line(aes(color = Immune, group = SubGroup), stat="smooth", method = "lm", alpha = 0.75, show.legend = FALSE)+ #adds linear model regression for each session (separately for immune and non-immune players)
  stat_summary(aes(), fun = "mean", geom = "line", size = 1.25)+ #non statistical, shows the time dynaics from round to round, visual aid only
  stat_summary(aes(), geom="linerange", fun.data=mean_cl_boot, fun.args=list(conf.int=0.95), linetype = 1, alpha = 1, color = "grey20")+ #adds 95% confidence intervals for each round
  #scale_linetype_manual(values = c(1,2), guide = FALSE)+ #did have different linetypes for immune not-immune but decided to use color aesthetic instaed
  coord_cartesian(xlim = c(0.5,5.5), ylim = c(-0.2,20), expand = FALSE)+
  stat_summary(fun = "mean", size=1.1, shape = 21)+ #shows a marker for the mean contribution in each round
  stat_summary(aes(colour = Immune), fun = "median", size=0.8, shape = 23, fill="white", show.legend = FALSE)+ #shows median, which is useful especially for non-normal distributions, but makes figure too crowded i think?
  scale_shape_manual(values=MyShapes, labels = c("Can be punished", "Immune (CC)"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished", "Immune (CC)"))+
  scale_color_manual(values = ColourChoice,labels = c("Can be punished", "Immune (CC)"))+
  labs(x="Game round (1-5)", y="Mean contribution (0-20 MU)", shape = "Player role", fill = "Player role", color = "Player role")+
  theme_classic()+
  theme(legend.position = "top")+
  facet_wrap(as.factor(CCfigdata$Scenario))
FigConCC +theme(axis.text=element_text(size=14),
                axis.title=element_text(size=16),
                legend.title=element_text(size=16),
                legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14, angle = 0))
#save as MS_Figure_Contributions
#ggsave("MS_Figure_ContributionsCC.pdf",height = 185, width = 180, units = c("mm")) #not used in ms


#Figure restricting all data to just Conditional Cooperators
CCfigdata <- SbjEditPGG[Stype_Human =="CC" & SMHvalid ==1]
dim(CCfigdata)
dim(SbjEditPGG) # difference is 165 rows, divided by 5 means 33 immune participants excluded
FigConCC <- ggplot(CCfigdata, aes(Period, Contribution, fill = Immune)) +
  geom_line(aes(color = Immune, group = SubGroup), stat="smooth", method = "lm", alpha = 0.75, show.legend = FALSE)+ #adds linear model regression for each session (separately for immune and non-immune players)
  stat_summary(aes(), fun = "mean", geom = "line", size = 1.25)+ #non statistical, shows the time dynaics from round to round, visual aid only
  stat_summary(aes(), geom="linerange", fun.data=mean_cl_boot, fun.args=list(conf.int=0.95), linetype = 1, alpha = 1, color = "grey20")+ #adds 95% confidence intervals for each round
  #scale_linetype_manual(values = c(1,2), guide = FALSE)+ #did have different linetypes for immune not-immune but decided to use color aesthetic instaed
  coord_cartesian(xlim = c(0.5,5.5), ylim = c(-0.2,20), expand = FALSE)+
  stat_summary(fun = "mean", size=1.1, shape = 21)+ #shows a marker for the mean contribution in each round
  stat_summary(aes(colour = Immune), fun = "median", size=0.8, shape = 23, fill="white", show.legend = FALSE)+ #shows median, which is useful especially for non-normal distributions, but makes figure too crowded i think?
  scale_shape_manual(values=MyShapes, labels = c("Can be punished (CC)", "Immune (CC)"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished (CC)", "Immune (CC)"))+
  scale_color_manual(values = ColourChoice,labels = c("Can be punished (CC)", "Immune (CC)"))+
  labs(x="Game round (1-5)", y="Mean contribution (0-20 MU)", shape = "Player role", fill = "Player role", color = "Player role")+
  theme_classic()+
  theme(legend.position = "top")+
  facet_wrap(as.factor(CCfigdata$Scenario))
FigConCC +theme(axis.text=element_text(size=14),
                axis.title=element_text(size=16),
                legend.title=element_text(size=16),
                legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14, angle = 0))
#save as MS_Figure_Contributions
ggsave("MS_Supplementary_Figure_4_conditional_cooperators.pdf",height = 185, width = 180, units = c("mm"))


#Supplementary figure on aggregate contribution patterns by session
SessionConts <- ggplot(SbjEditPGG, aes(Period, Contribution, color = Scenario, fill = Scenario, linetype = Immune, shape = Immune)) +
  stat_summary(aes(), fun = "mean", geom = "line", size = 0.5)+
  stat_summary(fun = "mean", size = 0.2)+
  stat_summary(aes(color = Scenario, linetype = Immune), geom="linerange", fun.data=mean_cl_boot, fun.args=list(conf.int=0.95), alpha = 1)+ #adds 95% confidence intervals for each round
  scale_shape_manual(values =c(21,1), labels = c("Can be punished", "Immune"))+
  scale_y_continuous(limits = c(0, 20), breaks = seq(0,20,5), expand = c(0,0)) +
  scale_linetype_manual(values = c("solid","dotted"), labels = c("Can be punished", "Immune"))+
  labs(x="Game round (1-5)", y="Mean contribution (0-20 MU)", linetype = "Player role", shape = "Player role", color = "Scenario", fill = "Scenario")+
  theme_bw()+
  #theme(legend.position = "top", legend.box = "vertical")+ #this moved the legend to the top, but was harder to read
  facet_wrap(~SessionType, ncol =4)
SessionConts + theme(axis.text=element_text(size=12),
                     axis.title=element_text(size=16),
                     legend.title=element_text(size=12),
                     legend.text=element_text(size=11),
                     strip.text.x = element_text(size = 12, angle = 0),
                     strip.text.y = element_text(size = 14))
#save as MS_Supplementary_Figure_Session_Contributions
ggsave("MS_Supplementary_Figure_1_Session_Contributions.pdf",height = 185, width = 180, units = c("mm"))


# Making a histogram of contributions
print("making histograms")
#unfortunately the sample sizes are not idenctical across scenarios, so one has to adjust for that to make relative frequencies, which are more comparable than raw numbers.

SbjEditPGGSample <- SbjEditPGG %>%
  group_by(Scenario, Immune, Period) %>%
  summarise(Sample = n_distinct(UPID))
SbjEditPGGSample

SbjEditPGGHIST <- SbjEditPGG %>%
  group_by(Contribution, Scenario, Immune, Period) %>%
  summarise(CAmount = sum(Contribution+1)) %>% #need to +1 so that when sum the zeroes it is not still zero
  mutate(Amount = CAmount / (Contribution+1)) #divide sum contributions+1 ==X by X+1 to get frequency
SbjEditPGGHIST
SbjEditPGGHIST <- data.table(SbjEditPGGHIST)
SbjEditPGGHIST[Scenario == "Mutual Punishers", Sample := 80]
SbjEditPGGHIST[Scenario == "Immune Punisher" & Immune == 0, Sample := 126]
SbjEditPGGHIST[Scenario == "Immune Punisher" & Immune == 1, Sample := 42]
SbjEditPGGHIST[Scenario == "Sole Punisher" & Immune == 0, Sample := 129]
SbjEditPGGHIST[Scenario == "Sole Punisher" & Immune == 1, Sample := 43]
SbjEditPGGHIST
SbjEditPGGHIST[,Perc := Amount*100/(Sample)]
SbjEditPGGHIST

SbjEditPGGHIST[Scenario == "Immune Punisher" & Period ==1 & Immune ==0]

datacheck <- SbjEditPGGHIST %>%
  group_by(Period, Scenario, Immune) %>%
  summarise(totalP = sum(Perc))
summary(datacheck$totalP) #totalP should equal 100 percent in all cases

# New facet label names for Period variable
Period.labs <- c("Round 1", "Round 2", "Round 3", "Round 4", "Round 5")
names(Period.labs) <- c("1","2","3","4","5")

# New facet label names for World variable
Scenario.labs <- c("Mutual Punishers (N=80)", "Immune Punisher (N=42+126)", "Sole Punisher (N=43+129)") #changed my mind and removed sample sizes as font has to be too small
names(Scenario.labs) <- c("Mutual Punishers", "Immune Punisher", "Sole Punisher")

FigHist <- 
  ggplot(SbjEditPGGHIST, aes(y=Perc/100, x=Contribution, colour = Immune, fill = Immune))+
  geom_col(position="dodge")+
  scale_colour_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  labs(x="Contribution (0-20 MU)", y="Relative frequency", fill = "Player role", color = "Player role")+
  theme_bw() +
  coord_cartesian(ylim = c(0,0.6))+
  theme(legend.position = "top")+
  scale_x_continuous(breaks = seq(0,20,10))+
  facet_grid(Scenario ~ Period, labeller = labeller(Period = Period.labs))
FigHist +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        strip.text.x = element_text(size = 12, angle = 0),
        strip.text.y = element_text(size = 14))
ggsave("MS_Supplementary_Figure_2_Histograms.pdf",height = 185, width = 180, units = c("mm"))

# Repeating making a histogram of contributions for Conditional cooperators only
print("making histograms CC")
#unfortunately the sample sizes are not idenctical across scenarios, so one has to adjust for that to make relative frequencies, which are more comparable than raw numbers.

#restricting data to CCs
SbjEditPGGcc <- SbjEditPGG[Stype_Human =="CC" & SMHvalid ==1]

SbjEditPGGSampleCC <- SbjEditPGGcc %>%
  group_by(Scenario, Immune, Period) %>%
  summarise(Sample = n_distinct(UPID))
SbjEditPGGSampleCC

SbjEditPGGHISTcc <- SbjEditPGGcc %>%
  group_by(Contribution, Scenario, Immune, Period) %>%
  summarise(CAmount = sum(Contribution+1)) %>% #need to +1 so that when sum the zeroes it is not still zero
  mutate(Amount = CAmount / (Contribution+1)) #divide sum contributions+1 ==X by X+1 to get frequency
SbjEditPGGHISTcc
SbjEditPGGHISTcc <- data.table(SbjEditPGGHISTcc)
SbjEditPGGHISTcc[Scenario == "Mutual Punishers", Sample := 54]
SbjEditPGGHISTcc[Scenario == "Immune Punisher" & Immune == 0, Sample := 72]
SbjEditPGGHISTcc[Scenario == "Immune Punisher" & Immune == 1, Sample := 24]
SbjEditPGGHISTcc[Scenario == "Sole Punisher" & Immune == 0, Sample := 96]
SbjEditPGGHISTcc[Scenario == "Sole Punisher" & Immune == 1, Sample := 28]
SbjEditPGGHISTcc
SbjEditPGGHISTcc[,Perc := Amount*100/(Sample)]
SbjEditPGGHISTcc

SbjEditPGGHISTcc[Scenario == "Immune Punisher" & Period ==1 & Immune ==0]

datacheck <- SbjEditPGGHISTcc %>%
  group_by(Period, Scenario, Immune) %>%
  summarise(totalP = sum(Perc))
summary(datacheck$totalP) #totalP should equal 100 percent in all cases

FigHistcc <- 
  ggplot(SbjEditPGGHISTcc, aes(y=Perc/100, x=Contribution, colour = Immune, fill = Immune))+
  geom_col(position="dodge")+
  scale_colour_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  labs(x="Contribution (0-20 MU)", y="Relative frequency among Conditional Cooperators.", fill = "Player role", color = "Player role")+
  coord_cartesian(ylim = c(0,0.6))+
  theme_bw() +
  theme(legend.position = "top")+
  scale_x_continuous(breaks = seq(0,20,10))+
  facet_grid(Scenario ~ Period, labeller = labeller(Period = Period.labs))
FigHistcc +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16),
        strip.text.x = element_text(size = 12, angle = 0),
        strip.text.y = element_text(size = 14))
ggsave("MS_Supplementary_Figure_3_HistogramsCC.pdf",height = 185, width = 180, units = c("mm"))


# Model 1: was cooperation stable in Mutual Punishers scenario?
MPdata <- SbjEditPGG[Scenario == "Mutual Punishers"] #restricting to one scenario (Mutual Punishers Scenario, MPS)
ModelMPS <- glmer(Contribution/Endowment ~ Period +(1|UPID) +(1|Date), data = MPdata, family = "binomial", weight = Endowment)
summary(ModelMPS) #no significant decline with time, controlling for session and UPID with random intercepts
#Period estimate = -0.02084 ± 0.0189, z value = -1.157, P = 0.2471

MPmeans <- MPdata %>%
  group_by(Immune, Period) %>%
  summarise(mean(Contribution))
MPmeans*5

MPmedians <- MPdata %>%
  group_by(Immune, Period) %>%
  summarise(median(Contribution))
MPmedians

# Model 2a: was cooperation stable in Immune Punishers scenario?
IPdata <- SbjEditPGG[Scenario == "Immune Punisher"]
IPdata <- within(IPdata, Immune <- relevel(Immune, ref = "0")) #change this from 0 to 1 if want immune players to be reference coefficient

ModelIPSa <- glmer(Contribution/Endowment ~ Period +(1|UPID) +(1|Date), data = IPdata, family = "binomial", weight = Endowment)
summary(ModelIPSa) #overall, significant decline with time, 
#Period estimate = -0.12019 ± 0.01325, z value = -9.074, P < 0.001

# Model 2b: was decline dependent on immunity in Immune Punisher scenario?
ModelIPSb <- glmer(Contribution/Endowment ~ Period * Immune +(1|UPID) +(1|Date), data = IPdata, family = "binomial", weight = Endowment)
summary(ModelIPSb) #significant interaction between immunity and period.
#Period estimate for non immune = -0.07770 ± 0.01455, z value = -5.338, P < 0.001, suggesting both immune and non-immune players declined
#Period estimate for  immune = -0.07770-0.24834 = -0.32611
#estimate of difference = -0.24834 ± 0.03582, z value = -6.933, P < 0.001, but decline for immune players was greater (faster)
#the decline was faster among the immune players

# Model 2c: were R1 contributions different depending on immunity?
ModelIPSc <- glmer(Contribution/Endowment ~ Immune+(1|Date), data = IPdata[Period==1], family = "binomial", weight = Endowment)
summary(ModelIPSc) #significant difference in round 1
#R1 estimate for non immune = -0.37660 ±0.10170
#R1 estimate for immune = -0.37660-0.44388
#estimate of difference = -0.44388 ±0.08543, z-value = -5.196, P < 0.001

IPmeans <- IPdata %>%
  group_by(Immune, Period) %>%
  summarise(mean(Contribution))
IPmeans
IPmeans*5 #for percentages

IPmedians <- IPdata %>%
  group_by(Immune, Period) %>%
  summarise(median(Contribution))
IPmedians

# Model 3a: was cooperation stable in Sole Punisher scenario?
SPdata <- SbjEditPGG[Scenario == "Sole Punisher"]
ModelSPSa <- glmer(Contribution/Endowment ~ Period +(1|UPID) +(1|Date), data = SPdata, family = "binomial", weight = Endowment)
summary(ModelSPSa) #overall, significant decline with time, 
#Period estimate = -0.32774 ± 0.01348, z value = -24.32, P < 0.001

# Model 3b: was decline dependent on immunity in Sole Punishers scenario?
ModelSPSb <- glmer(Contribution/Endowment ~ Period * Immune +(1|UPID) +(1|Date), data = SPdata, family = "binomial", weight = Endowment)
summary(ModelSPSb) #significant interaction between immunity and period.
#Period estimate for non immune = -0.28163 ± 0.01499, z value = -18.789, P < 0.001, suggesting both immune and non-immune players declined
#Period estimate for  immune = -0.28163-0.23198 = -0.51361
#estimate of difference = -0.23198 ± 0.03469, z value = -6.688, P < 0.001, but decline for immune players was greater (faster)
#the decline was faster among the immune players

# Model 3c: were R1 contributions different depending on immunity?
ModelSPSc <- glmer(Contribution/Endowment ~ Immune +(1|Date), data = SPdata[Period==1], family = "binomial", weight = Endowment)
summary(ModelSPSc) #significant difference in round 1, it's small, but significant because large sample size i guess
#R1 estimate for non immune = -0.15935 ±0.05936
#R1 estimate for immune = -0.15935-0.24167
#estimate of difference = -0.24167 ±0.08019, z-value = -3.014, P < 0.00258

SPmeans <- SPdata %>%
  group_by(Immune, Period) %>%
  summarise(mean(Contribution))
SPmeans
SPmeans*5

SPmedians <- SPdata %>%
  group_by(Immune, Period) %>%
  summarise(median(Contribution))
SPmedians

###repeating above analyses on contributions but restricting data for immune players to just those that categorize as Conditional Cooperators in prior strategy method with humans.
SMfreqs <- SbjEdit[Stage==3 & Period==1 & SMHvalid==1 ] %>% #important to restrict to SMHvalid cases.
  group_by(Scenario, Immune, Stype_Human) %>%
  summarise(N = n()) #this gives the breakdown of social types per scenario and immunity, 24 valid immune conditional cooperators in immune punisher scenario
SMfreqs

#can restrict to just CC for only Immune participants (1st row), or all participants (2nd row), results basically the same
#1
SbjEditPGGcc <- SbjEditPGG[Immune == 0 | (Immune == 1 & Stype_Human =="CC" & SMHvalid ==1)] #this more useful for seeing how immune CCs react to rest of group, if including all non immune individuals, regardless of Strategy Method, don't need to be valid.
#2
#SbjEditPGGcc <- SbjEditPGG[Stype_Human =="CC" & SMHvalid ==1] #this more useful for comparing directly among CCs if immunity matters.

# Model 2a: was cooperation stable in Immune Punishers scenario?
IPdata <- SbjEditPGGcc[Scenario == "Immune Punisher"]
ModelIPSa <- glmer(Contribution/Endowment ~ Period +(1|UPID) +(1|Date), data = IPdata, family = "binomial", weight = Endowment)
summary(ModelIPSa) #overall, still significant decline with time

# Model 2b: was decline dependent on immunity in Immune Punisher scenario?
ModelIPSb <- glmer(Contribution/Endowment ~ Period * Immune +(1|UPID) +(1|Date), data = IPdata, family = "binomial", weight = Endowment)
summary(ModelIPSb) #still significant interaction between immunity and period.
#estimate of difference = -0.23153 ± 0.04153, z value = -5.575, P < 0.001, decline for immune players was still greater (faster) even if just CCs


# Model 2c: were R1 contributions different depending on immunity?
ModelIPSc <- glmer(Contribution/Endowment ~ Immune+(1|Date), data = IPdata[Period==1], family = "binomial", weight = Endowment)
summary(ModelIPSc) #No significant difference in round 1, which fits confused learners, as CC start cooperative but decline, even though we control for stable social norm
#estimate of difference = 0.08074 ±0.10514, z-value = 0.768, P = 0.442560

IPmeansCC <- IPdata %>%
  group_by(Immune, Period) %>%
  summarise(mean(Contribution))
IPmeansCC
IPmeansCC*5 #for percentages

IPmediansCC <- IPdata %>%
  group_by(Immune, Period) %>%
  summarise(median(Contribution))
IPmediansCC

# Model 3a: was cooperation stable in Sole Punisher scenario?
SPdata <- SbjEditPGGcc[Scenario == "Sole Punisher"]

# Model 3b: was decline dependent on immunity in Sole Punishers scenario?
ModelSPSb <- glmer(Contribution/Endowment ~ Period * Immune +(1|UPID) +(1|Date), data = SPdata, family = "binomial", weight = Endowment)
summary(ModelSPSb) # Still significant interaction between immunity and period.
#estimate of difference = -0.23327 ± 0.03780, z value = -6.171, P < 0.001, decline for immune players was still greater (faster) even if just CCs

# Model 3c: were R1 contributions different depending on immunity?
ModelSPSc <- glmer(Contribution/Endowment ~ Immune +(1|Date), data = SPdata[Period==1], family = "binomial", weight = Endowment)
summary(ModelSPSc) #Just significant difference in round 1, but in wrong direction
#estimate of difference = 0.18876 ±0.09453, z-value = 1.997, P < 0.0459

SPmeansCC <- SPdata %>%
  group_by(Immune, Period) %>%
  summarise(mean(Contribution))
SPmeansCC*5

SPmediansCC <- SPdata %>%
  group_by(Immune, Period) %>%
  summarise(median(Contribution))
SPmediansCC

#Making figure of R1 and R5, by FR and CC and Immune or not
ColourChoice2 = c("#CC79A7", "black")
PrefsPlot <- ggplot(SbjEditPGG[Stype_Human != "Other" & (Period == 1 | Period ==5)  & SMHvalid==1  & Scenario != "Sole Punisher"], aes(Period, Contribution, colour = Immune, shape = Immune))+
  geom_point(aes(alpha = Immune, shape = Immune, fill = Immune))+
  #geom_line(aes(color = Immune, group = UPID), stat="smooth", method = "lm", alpha = 0.25, show.legend = FALSE)+ #adds linear model regression for each session (separately for immune and non-immune players)
  stat_summary(aes(color = Immune, group = UPID, alpha = Immune, size = Immune), fun = "mean", geom = "line")+
  stat_summary(fun = "mean", geom = "line", size =1.8)+
  stat_summary(aes(), fun="mean", fill = "white")+
  stat_summary(aes(), geom="linerange", fun.data=mean_cl_boot, fun.args=list(conf.int=0.95), linetype = 1, alpha = 1)+ #adds 95% bootstrapped confidence intervals for each round
  scale_size_manual(values = c(0.33, 0.66), labels = c("Can be punished", "Immune"))+
  scale_shape_manual(values = c(21,21), labels = c("Can be punished", "Immune"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  scale_color_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  scale_alpha_manual(values = c(0.25, 0.5),labels = c("Can be punished", "Immune"))+
  labs(x="Game round (1-5)", y="Contribution (0-20 MU)", size = "Player role", shape = "Player role", fill = "Player role", color = "Player role", alpha = "Player role")+
  facet_grid(~Stype_Human)+
  theme_classic()+
  theme(legend.position = "top")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14, angle = 0))
xlabs <- c("1", "", "", "", "5")
PrefsPlot  + scale_x_continuous(labels= xlabs)
ggsave("MS_Figure_3_Preferences_and_deterrence.pdf",height = 185, width = 180, units = c("mm"))
#The threat of punishment converts Free Riders into cooperators. Experience with immunity turns Conditional Cooperators into Free Riders.


Ns <- SbjEditPGG[ (Period == 1 | Period ==5)  & SMHvalid & Scenario != "Sole Punisher"] %>%
  group_by(Stype_Human, Immune) %>%
  summarise(N = n_distinct(UPID))
Ns #others have been excluded

##################### analyzing punishment ##################
print("Step 4.2 - Analyzing Punishment")
#############################################################

#Frequency of punishment
#step one, use dplyr to work out how many individuals punished at one point 
#step two, and work out how many individuals punished only one round, or multiple rounds

dim(SbjEditPGG) /420 #This variable we calculated earlier, scores a 1 if the player punished anyone in that round
#now to use Dplyr to see how many 1s a player gets across the 5 rounds


summary1 <- SbjEditPGG %>%
  group_by(UPID, Scenario, Immune) %>% #should be 5 rows per UPID
  summarise(timesP = sum(as.numeric(Punished))-5) %>% #this value tells you how many rounds the UPID punished in
  group_by(Scenario, Immune, timesP) %>%
  summarise(N = n()) #now this counts how many upids (rows) there are for each value of timesP
summary1 # Supplementary Table 2

#step three, work out how many times someone punished 0, 1, 2, or 3 groupmates
#all by scenario and immune status

#how many times an individual punished 0,1,2 or 3 groupmates, Supplementary Table 1
SbjEditPGG[,PunMP :="blank"]
SbjEditPGG[Scenario=="Mutual Punishers", PunMP := if_else(MyDeductionofP1 > 0 & MyDeductionofP2 > 0 & MyDeductionofP3 > 0,"pun3", if_else(MyDeductionofP1 > 0 & MyDeductionofP2 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP1 > 0 & MyDeductionofP3 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP2 > 0 & MyDeductionofP3 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP1 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP2 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP3 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP4 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(PunishmentSpending ==0,"pun0", "pun2")))))))))]
summary(as.factor(SbjEditPGG$PunMP))

SbjEditPGG[,PunIPni :="blank"]
SbjEditPGG[Scenario=="Immune Punisher" & Immune == 0, PunIPni := if_else(MyDeductionofP1 > 0 & MyDeductionofP2 > 0 & MyDeductionofP3 > 0,"pun3", if_else(MyDeductionofP1 > 0 & MyDeductionofP2 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP1 > 0 & MyDeductionofP3 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP2 > 0 & MyDeductionofP3 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP1 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP2 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP3 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP4 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(PunishmentSpending ==0,"pun0", "pun2")))))))))]
summary(as.factor(SbjEditPGG$PunIPni))

SbjEditPGG[,PunIPi :="blank"]
SbjEditPGG[Scenario=="Immune Punisher" & Immune == 1, PunIPi := if_else(MyDeductionofP1 > 0 & MyDeductionofP2 > 0 & MyDeductionofP3 > 0,"pun3", if_else(MyDeductionofP1 > 0 & MyDeductionofP2 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP1 > 0 & MyDeductionofP3 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP2 > 0 & MyDeductionofP3 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP1 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP2 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP3 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP4 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(PunishmentSpending ==0,"pun0", "pun2")))))))))]
summary(as.factor(SbjEditPGG$PunIPi))

SbjEditPGG[,PunSP :="blank"]
SbjEditPGG[Scenario=="Sole Punisher" & Immune == 1, PunSP := if_else(MyDeductionofP1 > 0 & MyDeductionofP2 > 0 & MyDeductionofP3 > 0,"pun3", if_else(MyDeductionofP1 > 0 & MyDeductionofP2 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP1 > 0 & MyDeductionofP3 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP2 > 0 & MyDeductionofP3 > 0 & MyDeductionofP4 > 0,"pun3", if_else(MyDeductionofP1 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP2 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP3 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(MyDeductionofP4 == PunishmentSpending & PunishmentSpending > 0,"pun1", if_else(PunishmentSpending ==0,"pun0", "pun2")))))))))]
summary(as.factor(SbjEditPGG$PunSP))


GroupPunCosts <- SbjEditPGG %>%
  filter(SessionNumber > 2) %>% #sessions 1 and 2 had a smaller punishment budget of 18 instead of 30, limit of -6 spending per victim, so filter when looking at earnings
  group_by(Scenario, Immune) %>% #do with or without immune accordingly
  summarise(meanSpend = mean(MySumDeductionofOthers),
            meanLoss = mean(GroupPunOfMe),
            meanCont = mean(Contribution),
            meanSuruplus = mean((Contribution*0.6)),
            meanEarnings = mean(Profit))
GroupPunCosts

#doing differently, to get bootstrapped confidence intervals also
if(!require(Hmisc)){
  install.packages("Hmisc")
  library(Hmisc)
}

print("testing if Immune players punish less")

#summarizing Mutual Punishers scenario, punishment and earnings
mean_cl_boot(SbjEditPGG[Scenario=="Mutual Punishers"]$PunishmentSpending)
mean_cl_boot(SbjEditPGG[Scenario=="Mutual Punishers"]$Contribution*0.6)
mean_cl_boot(SbjEditPGG[Scenario=="Mutual Punishers"]$Profit)


#summarizing Immune Punisher scenario
mean_cl_boot(SbjEditPGG[Scenario=="Immune Punisher" & Immune ==0 & SessionNumber >0]$PunishmentSpending) #Sessions 1 and 2 had a smaller punishment budget/limit of 18 MU and 6 MU max spend per target, instead of 30 and 10, but doesn't seem to matter.
mean_cl_boot(SbjEditPGG[Scenario=="Immune Punisher" & Immune ==1]$PunishmentSpending)
mean_cl_boot(SbjEditPGG[Scenario=="Immune Punisher"]$Contribution*0.6)
mean_cl_boot(SbjEditPGG[Scenario=="Immune Punisher" & Immune ==0]$Profit)
mean_cl_boot(SbjEditPGG[Scenario=="Immune Punisher" & Immune ==1]$Profit)

PunIP <- glmer(PunishmentSpending/PunishmentBudget ~ Immune+as.factor(Period) +(1|UPID) +(1|Date), data = SbjEditPGG[Scenario=="Immune Punisher"], family = "binomial", weight = PunishmentBudget)
summary(PunIP) #Controlling for period within the Immune Punisher scenario, the immune individuals spent significantly less on punishment, glmer controlling for UPID and session, estimate = -0.8587 ±0.3277, z value = -2.621, P = 0.00877

SbjEditPGG[, PerTargetSpend := if_else(Scenario == "Immune Punisher" & Immune ==0, PunishmentSpending/2, PunishmentSpending/3)]

NImodel <- lmer(PerTargetSpend ~ Scenario +(1|UPID), data = SbjEditPGG[Immune ==0 & Scenario != "Sole Punisher"])
summary(NImodel) #did spending per potential target differ for non-immune individuals between mutual punishers and immune punisher scenarios?
basemodel <- lmer(PerTargetSpend ~ 1 +(1|UPID), data = SbjEditPGG[Immune ==0 & Scenario != "Sole Punisher"])
summary(basemodel) #basal comparison
anova(NImodel, basemodel) #LRT

#Summarizing Sole Punisher scenario
mean_cl_boot(SbjEditPGG[Scenario=="Sole Punisher" & Immune ==1]$PunishmentSpending)
mean_cl_boot(SbjEditPGG[Scenario=="Sole Punisher" & Immune ==0]$PunishmentSpending)
mean_cl_boot(SbjEditPGG[Scenario=="Sole Punisher"]$Contribution*0.6)
mean_cl_boot(SbjEditPGG[Scenario=="Sole Punisher" & Immune ==1]$Profit)
mean_cl_boot(SbjEditPGG[Scenario=="Sole Punisher" & Immune ==0]$Profit)

PunImmune <- glmer(PunishmentSpending/PunishmentBudget ~ Scenario+ as.factor(Period) +(1|UPID), data = SbjEditPGG[Immune ==1], family = "binomial", weight = PunishmentBudget)
summary(PunImmune) #Controlling for period within the Immune Punisher scenario, the immune individuals spent significantly less on punishment, glmer controlling for UPID and session, estimate = -0.85855 ±0.32760, z value = -2.621, P = 0.00877


#making figure of punishment rates (per capita spending)
FigPunZoom <- ggplot(SbjEditPGG[Immune ==1 | Scenario != "Sole Punisher"], aes(Period, PunishmentSpending, colour = Immune, fill = Immune)) +
  geom_line(aes(color = Immune, group = SubGroup), stat="smooth", method = "lm", alpha = 0.75, show.legend = FALSE)+
  stat_summary(aes(), fun = "mean", geom = "line", size = 1.25, colour = "black")+
  stat_summary(aes(), geom="linerange", fun.data=mean_cl_boot, fun.args=list(conf.int=0.95), linetype = 1, alpha = 1, color = "grey20")+
  #scale_linetype_manual(values = c(1,2), guide = FALSE)+
  #scale_y_continuous(limits=c(0,30), breaks = seq(0,30,5))+
  coord_cartesian(ylim = c(0,10))+ #individuals could spend up to 30 MU, but mean is around 2-3, so axis is truncated (zoomed in) at 10 to enable comparisons
  #stat_summary(aes(color = Immune),fun = "max", size=0.55, shape =17)+
  stat_summary(fun = "mean", size=1.1, shape = 21, colour = "black")+
  stat_summary(aes(colour = Immune), fun = "median", size=0.8, shape = 23, fill="white", show.legend = FALSE)+
  #scale_shape_manual(values=MyShapes, labels = c("Non-immune", "Immune"))+
  scale_shape_manual(values=MyShapes, labels = c("Can be punished", "Immune"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  scale_color_manual(values = ColourChoice,labels = c("Can be punished", "Immune"))+
  labs(x="Game round (1-5)", y="Mean punishment spending (0-30 MU)", shape = "Player role", fill = "Player role", color = "Player role")+
  theme_classic()+
  theme(legend.position = "top")+
  facet_wrap(~Scenario)
FigPunZoom +theme(axis.text=element_text(size=14),
                  axis.title=element_text(size=16),
                  legend.title=element_text(size=16),
                  legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14, angle = 0))
ggsave("MS_Figure_4_PunishmentSpending.pdf",height = 185, width = 180, units = c("mm"))

#make figure of number of punishers (punisher frequency) per round
FigPunFreq <- ggplot(SbjEditPGG[Immune ==1 | Scenario != "Sole Punisher"], aes(x=Period, y=as.numeric(Punished)-1, color = Immune, fill = Immune))+
  stat_smooth(aes(),method = "lm", size = 1.25, se=TRUE, alpha=0.5, show.legend = FALSE)+
  stat_summary(aes(), fun = "mean", geom = "line", size = 1.25, colour = "black")+
  stat_summary(fun = "mean", size=1.1, shape = 21, colour = "black")+
  scale_shape_manual(values=MyShapes, labels = c("Can be punished", "Immune"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  scale_color_manual(values = ColourChoice,labels = c("Can be punished", "Immune"))+
  coord_cartesian(ylim = c(0,1))+
  #stat_summary(aes(), geom="linerange", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), linetype = 1, alpha = 1, color = "grey20")+ #adds 95% confidence intervals for each round
  labs(x="Game round (1-5)", y="Frequency of punishers (0-1)", shape = "Player role", fill = "Player role", color = "Player role")+
  theme_classic()+
  theme(legend.position = "top")+
  facet_wrap(~Scenario)
FigPunFreq +theme(axis.text=element_text(size=14),
                  axis.title=element_text(size=16),
                  legend.title=element_text(size=16),
                  legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14, angle = 0))
ggsave("MS_Supplementary_Figure_6_PunishmentFrequency.pdf",height = 185, width = 180, units = c("mm"))

#calculate raw values of pun frequencies for reporting
PunPerRound <- SbjEditPGG %>%
  group_by(Scenario, Immune, Period) %>%
  summarise(sumPun = sum(as.numeric(Punished)-1),
            sampleN = n()) %>%
  mutate(Freq = sumPun / sampleN)
PunPerRound
#Mutual punishers, frequency goes from 53% to 64% then 60% until end, marginal increase or constant
#Immune Punishers, non immune, freq goes from 52% and ends with 52%, a high of 56%, very constant
#Immune punishers immune, freq never reaches over 50%, goes 45 to 36%, with a low of 31%, significant decline?
#Sole punishers, freq goes from 58 to 47% steady decline, significant?

#statistical models of frequency over time
#requires logistic regression 0/1
FreqData <- SbjEditPGG[Immune ==1 | Scenario != "Sole Punisher"] #excluding non-immune individuals in Sole Punisher scenario who could not punish
ModelFreqs <- glmer(Punished ~ Period*Immune +(1|UPID), data = FreqData, family = "binomial")
summary(ModelFreqs) #significant interaction, faster decline in immune together and control for scenario than in non-immune

FreqData <- SbjEditPGG[Immune ==1] #excluding non-immune individuals in Sole Punisher scenario who could not punish
ModelFreqs <- glmer(Punished ~ Period+Scenario +(1|UPID), data = FreqData, family = "binomial")
summary(ModelFreqs) #significant decline in immune individuals, controlling for scenario 


############################################################
print("ANALYSING EARNINGS FIGURE")
############################################################

FigEarnings <- ggplot(SbjEditPGG[SessionNumber > 2], aes(Period, MyFinalEarnings, shape = Immune, fill = Immune)) +
  annotate(geom="text", x=3, y=52, label="Endowment (50 MU)",
           color="orange")+
  annotate(geom="text", x=3, y = 64, label="Cooperative maximum (62 MU)", colour = "dark green")+
  annotate(geom="text", x=3, y = 75.5, label="Selfish maximum (74 MU)", colour = "black")+
  geom_line(aes(color = Immune, group = SubGroup), stat="smooth", method = "lm", alpha = 0.75, show.legend = FALSE)+
  stat_summary(aes(), fun = "mean", geom = "line", size = 1.25)+
  stat_summary(aes(), geom="linerange", fun.data=mean_cl_boot, fun.args=list(conf.int=0.95), linetype = 1, alpha = 1, color = "grey20")+
  #scale_linetype_manual(values = c(1,2), guide = FALSE)+
  scale_y_continuous(limits=c(0,76), breaks = seq(0,75,15))+
  #coord_cartesian(ylim = c(0,75)))+
  #stat_summary(aes(color = Immune),fun = "max", size=0.55, shape =17)+
  geom_line(y=50, colour = "orange", size = 1.5, linetype =1.25)+
  geom_line(y=62, colour = "darkgreen", size = 1.5, linetype =1.25)+
  geom_line(y=74, colour = "black", size = 1.5, linetype =1.25)+
  stat_summary(fun = "mean", size=1)+
  stat_summary(aes(colour = Immune), fun = "max", size=0.55, shape = 2, alpha = 0.95, show.legend = FALSE, fill = "white")+
  stat_summary(aes(colour = Immune), fun = "min", size=0.55, shape = 6, alpha = 0.95, show.legend = FALSE)+
  scale_shape_manual(values=MyShapes, labels = c("Can be punished", "Immune"))+
  scale_fill_manual(values = ColourChoice, labels = c("Can be punished", "Immune"))+
  scale_color_manual(values = ColourChoice,labels = c("Can be punished", "Immune"))+
  labs(x="Game round", y="Final earnings (0-74 MU)")+
  labs(x="Game round (1-5)", y="Final earnings (0-74 MU)", shape = "Player role", fill = "Player role", color = "Player role")+
  theme_classic()+
  theme(legend.position = "top")+
  facet_wrap(~Scenario)
FigEarnings +theme(axis.text=element_text(size=14),
                   axis.title=element_text(size=16),
                   legend.title=element_text(size=16),
                   legend.text=element_text(size=16))+
  theme(strip.text.x = element_text(size = 14, angle = 0))
#save as MS_Figure_Earnings
ggsave("MS_Supplementary_Figure_5_Earnings.pdf",height = 185, width = 180, units = c("mm"))

#endowment = 20 MU for stage 1, and 30 MU for stage 2 (except was 18 in sessions 1 and 2), so starting point is 50 MU
#cooperative equilibrium is full cooperate and no punishment, =(4x20)*0.4 = 32 + 30 = 62 MU
#free rider maximum is full defect with full cooperators, and no punishment spending, = 20 + (3*(20*0.4)) = 44+30 = 74 MU

#Statistically testing if earnings were different to endowment.
TotalEarnings <- SbjEditPGG %>% #summarizing earnings per individual across all 5 rounds
  filter (SessionNumber > 2) %>% #exclude sessions with smaller punishment budget endowment
  group_by(UPID, Scenario, Immune, Date) %>%
  summarise(sumProfit = sum(Profit))
TotalEarnings
TotalEarnings <- data.table(TotalEarnings)

TP <- TotalEarnings$sumProfit
shapiro.test(TP)


MPe <- TotalEarnings[Scenario=="Mutual Punishers"]$sumProfit
ks.test(MPe, "pnorm", mean=mean(MPe), sd=sd(MPe))
shapiro.test(MPe)
mean_cl_boot(MPe) #reference is 5*50 = 250, upper bound is below
250/215.2 #average earnings would have been 16% higher if no one had punished nor contributed

IPnie <- TotalEarnings[Scenario=="Immune Punisher" & Immune == "0"]$sumProfit
ks.test(IPnie, "pnorm", mean=mean(IPnie), sd=sd(IPnie))
shapiro.test(IPnie)
mean_cl_boot(IPnie) #upper bound is below 250, the decimal point values on these changed after a reinstalling of r-studio after some problems. The change is super minor, but I do not know what causes it. If you run this on your machine, it may give slightly differnt values I guess.

IPie <- TotalEarnings[Scenario=="Immune Punisher" & Immune == "1"]$sumProfit
ks.test(IPie, "pnorm", mean=mean(IPie), sd=sd(IPie))
shapiro.test(IPie)
mean_cl_boot(IPie) #lower bound is above 250

SPnie <- TotalEarnings[Scenario=="Sole Punisher" & Immune == "0"]$sumProfit
ks.test(SPnie, "pnorm", mean=mean(SPnie), sd=sd(SPnie))
shapiro.test(SPnie)
mean_cl_boot(SPnie) # 

SPie <- TotalEarnings[Scenario=="Sole Punisher" & Immune == "1"]$sumProfit
ks.test(SPie, "pnorm", mean=mean(SPie), sd=sd(SPie))
shapiro.test(SPie)
mean_cl_boot(SPie) #CIs straddle 250

#Q: did earnings differ by scenario?
ModelBig<- lmer(sumProfit ~ Scenario +(1|Date), data = TotalEarnings)
summary(ModelBig)
anova(ModelBig)

ModelSmall <- lmer(sumProfit ~ 1 +(1|Date), data = TotalEarnings)
summary(ModelSmall)
anova(ModelSmall)
anova(ModelBig, ModelSmall) #Earnings significantly varied by Scenario. X2 = 20.0, df=2, P < 0.001

#Q: did earnings differ by immunity?
ModelBig<- lmer(sumProfit ~ Immune+Scenario +(1|Date), data = TotalEarnings)
summary(ModelBig)
anova(ModelBig)

ModelSmall <- lmer(sumProfit ~ Scenario +(1|Date), data = TotalEarnings)
summary(ModelSmall)
anova(ModelSmall)
anova(ModelBig, ModelSmall) #

#Q: but did this depend on scenario? 
ModelBig<- lmer(sumProfit ~ Immune*Scenario +(1|Date), data = TotalEarnings)
summary(ModelBig)
anova(ModelBig)

ModelSmall <- lmer(sumProfit ~ Immune+Scenario +(1|Date), data = TotalEarnings)
summary(ModelSmall)
anova(ModelSmall)
anova(ModelBig, ModelSmall) # Yes, effect of immunity depended on scenario, X2 = 79.7, df = 1, P < 0.0001

#Q: within Immune Punisher, did immunity affect earnings?
ModelBig<- lmer(sumProfit ~ Immune +(1|Date), data = TotalEarnings[Scenario=="Immune Punisher"])
summary(ModelBig)
anova(ModelBig)

ModelSmall <- lmer(sumProfit ~ 1 +(1|Date), data = TotalEarnings[Scenario=="Immune Punisher"])
summary(ModelSmall)
anova(ModelSmall)
anova(ModelBig, ModelSmall) #Immune Punishers made significantly more than non-immune in IP scenario, X2 = 132.7, df = 1, P < 0.001, estimate was +58.2 ±3.93 MU

#Q: and within Sole Punisher, did immunity affect earnings?
ModelBig<- lmer(sumProfit ~ Immune +(1|Date), data = TotalEarnings[Scenario=="Sole Punisher"])
summary(ModelBig)
anova(ModelBig)

ModelSmall <- lmer(sumProfit ~ 1 +(1|Date), data = TotalEarnings[Scenario=="Sole Punisher"])
summary(ModelSmall)
anova(ModelSmall)
anova(ModelBig, ModelSmall) #Sole Punisher made significantly more than non immune in SP scenario, X2 = 13.9, df = 1, P < 0.001, estimate was +11.3 MU ±2.99


#did any group reach full contributions ('paradise', SumC == 80/80)?
paradise <- SbjEditPGG %>%
  group_by(Scenario, UGID) %>%
  summarise(SumC = mean(SumC),
            Scenario = first(Scenario)) 
paradise
paradise <- data.table(paradise)

ggplot(paradise[Scenario == "Mutual Punishers"], aes(SumC))+
  geom_histogram(binwidth = 1) +
  facet_wrap(~Scenario)


paradiselimits <- paradise %>%
  group_by(Scenario) %>%
  summarise(maxC = max(SumC),
            minC = min(SumC))
paradiselimits #no the highest contributing group was 63/80, in the MP scenario. 63/80 is 15.75 per groupmember, 79%.

############################### Social, anti-social, and hypocritical punishment ################

# step one
#calculating for each player if PlayerX is a possible target for either Hypocritical (HT), anti-social (antiT), or social-punishment (socT)
#can then look to see if they tartgeted PlayerX, and then code the punishment accordingly

#calculate social punishment, calculate anti-social punishment

summary(as.factor(SbjEditPGG$P1socT))

summary(as.factor(SbjEditPGG$P1antiT))
nrow(SbjEditPGG[RandomPlayerID==4])

SbjEditPGG[,P1socT := if_else(Contribution > ContPlayer1, 1, 0)][,P2socT := if_else(Contribution > ContPlayer2, 1, 0)][,P3socT := if_else(Contribution > ContPlayer3, 1, 0)][,P4socT := if_else(Contribution > ContPlayer4, 1, 0)]
# This caluculates if PlayerX(1-4), PX, a potential target for social punishment, i.e. if PX contributed less than the focal player, in which case the punishment could be considered 'social'

SbjEditPGG[, P1antiT := if_else(Contribution <= ContPlayer1 & RandomPlayerID != 1, 1, 0)][, P2antiT := if_else(Contribution <= ContPlayer2 & RandomPlayerID != 2, 1, 0)][, P3antiT := if_else(Contribution <= ContPlayer3 & RandomPlayerID != 3, 1, 0)][, P4antiT := if_else(Contribution <= ContPlayer4 & RandomPlayerID != 4, 1, 0)]

SbjEditPGG[, P1antiT := if_else(P1socT ==0 & RandomPlayerID !=1, 1, 0)][, P2antiT := if_else(P2socT ==0 & RandomPlayerID !=2, 1, 0)][, P3antiT := if_else(P3socT ==0 & RandomPlayerID !=3, 1, 0)][, P4antiT := if_else(P4socT ==0 & RandomPlayerID !=4, 1, 0)] #have to exclude those cases where PX is actually just the focal player (RandomPlayerID)

summary(as.factor(SbjEditPGG$P1socT)) #1-4
summary(as.factor(SbjEditPGG$P1antiT)) #bit more targets for anti because of ties (both contribute 20, 10, or 0 for example)
# total should equal 1575 each time (2,100 rows, but 1/4, 525, are self, leaving 1,575 to be either social or anti-social)

#have now calculated targets as either social or anti-social.
#now calcuLated if only punished social, anti or both.

SbjEditPGG[, ididsocP := if_else((MyDeductionofP1 > 0 & P1socT ==1 | MyDeductionofP2 > 0 & P2socT ==1| MyDeductionofP3 > 0 & P3socT ==1 | MyDeductionofP4 > 0 & P4socT ==1), 1,0)] #if punished socially
summary(as.factor(SbjEditPGG$ididsocP)) #566 times someone punished at least one person socially. But may have done anti-social too.

SbjEditPGG[, ididantiP := if_else(MyDeductionofP1 > 0 & P1antiT ==1 | MyDeductionofP2 > 0 & P2antiT ==1| MyDeductionofP3 > 0 & P3antiT ==1 | MyDeductionofP4 > 0 & P4antiT ==1, 1,0)] #if punished 'anti-socially', considered an anti-social punisher [note at this stage can be multiple things]
summary(as.factor(SbjEditPGG$ididantiP)) #250 times, someone punished at least one person anti-socially. But may have done social too.

SbjEditPGG[, myPtype := if_else(ididantiP ==0 & ididsocP ==0, "none", if_else(ididantiP ==0 & ididsocP ==1, "social", if_else(ididantiP ==1 & ididsocP ==0, "anti", if_else(ididantiP ==1 & ididsocP ==1, "both", "99"))))]#99 is just to check no mistakes made, see a 99 and something has not been classified
summary(as.factor(SbjEditPGG$myPtype)) #need to summarise by scenario and immunity
#bear in mind that an individual cannot punish below (or both) if they are the lowest contributor. A top contributor can still punish 'anti-socially' if they tie with someone and punish them.        

#using dplyr to summarise by scenario and immunity
PTFs <- SbjEditPGG %>%
  group_by(Scenario, Immune) %>%
  summarise(Nsocial = sum(myPtype == "social"),
            Nanti = sum(myPtype == "anti"),
            Nboth = sum(myPtype =="both"),
            Nnon = sum(myPtype =="none"),
            N = n()) %>% 
  mutate(PercSoc = Nsocial/(N-Nnon))
PTFs #Punishment Type FrequencieS

Ftesting <-
  matrix(c((188+219), (32+17+85+27), (38+60), (32+7+40+10)), #these are the values in PTFs summary
         nrow = 2,
         dimnames = list(Normal = c("Social", "NonSoc"),
                         Immune = c("Social", "NonSoc")))
fisher.test(Ftesting, alternative = "two.sided")



#we differentiate between anti-social punishment, against top contributor(s), and hypocritical punishment, which is when an individual punishes those that contribute more than them but not the top contributor(s). Worth nothing that all anti-social punishment could be hypocritical, but here we conservatively estimate levels of hypocrisy by ruling out cases where they punish the top. If a player only punishes intermediate contributors, and they are a lower conrtributor, this is undeniably hypocritical.  

#first we calculate if a groupmember was a potential target of anti-social punishment under this new framework
SbjEditPGG[,P1antiT := if_else(Contribution < ContPlayer1 & ContPlayer1 == maxC, 1, 0)][,P2antiT := if_else(Contribution < ContPlayer2 & ContPlayer2 == maxC, 1, 0)][,P3antiT := if_else(Contribution < ContPlayer3 & ContPlayer3 == maxC, 1, 0)][,P4antiT := if_else(Contribution < ContPlayer4 & ContPlayer4 == maxC, 1, 0)]
# if P1antiT ==1 means Player 1 (P1) was a potential target (T) for antisocial punishment (anti) from focal player's perspective. This calculates if PX was the top or equal top contributor in the group, in which case punishing them could be considered anti-social punishment

#then we calculate if a groupmember was a potential target for hypocritical punishment.
SbjEditPGG[,P1HT := if_else(Contribution < ContPlayer1 & ContPlayer1 < maxC | Contribution == ContPlayer1 & RandomPlayerID !=1, 1, 0)][,P2HT := if_else(Contribution < ContPlayer2 & ContPlayer2 < maxC | Contribution == ContPlayer2 & RandomPlayerID !=2, 1, 0)][,P3HT := if_else(Contribution < ContPlayer3 & ContPlayer3 < maxC | Contribution == ContPlayer3 & RandomPlayerID !=3, 1, 0)][,P4HT := if_else(Contribution < ContPlayer4 & ContPlayer4 < maxC | Contribution == ContPlayer4 & RandomPlayerID !=4, 1, 0)]
# Was PlayerX a potential target for Hypocritical Punishment? this calculates if PX contributed the  more than the focal player, but was not the highest contributor, or either the same as focal player, in which case punishing them could be considered 'hypocritical' (it's not antisoc pun if they did not pun the top contributor).

#Then we work out if players did any acts consistent with social, anti-social, and/or hypocritical punishment
SbjEditPGG[, ididantiP := if_else(MyDeductionofP1 > 0 & P1antiT ==1 | MyDeductionofP2 > 0 & P2antiT ==1| MyDeductionofP3 > 0 & P3antiT ==1 | MyDeductionofP4 > 0 & P4antiT ==1, 1,0)] 
summary(as.factor(SbjEditPGG$ididantiP)) #123 times someone did anti-social

SbjEditPGG[, ididHP := if_else((MyDeductionofP1 > 0 & P1HT ==1 | MyDeductionofP2 > 0 & P2HT ==1| MyDeductionofP3 > 0 & P3HT ==1 | MyDeductionofP4 > 0 & P4HT ==1), 1,0)] 
summary(as.factor(SbjEditPGG$ididHP)) #204 times someone did hypocritical

SbjEditPGG[, ididsocP := if_else((MyDeductionofP1 > 0 & P1socT ==1 | MyDeductionofP2 > 0 & P2socT ==1| MyDeductionofP3 > 0 & P3socT ==1 | MyDeductionofP4 > 0 & P4socT ==1), 1,0)] 
summary(as.factor(SbjEditPGG$ididsocP)) #566 times someone did social

#Then we assign type, to social, anti-social, hypocritical, other, and non
#we then have to reassign type, with the new possibility of being hypocritical.
SbjEditPGG[, myPtype3 := if_else(ididantiP ==1 & ididsocP ==0, "anti", if_else(ididHP ==1 & ididantiP ==0 ,"Hypo", if_else(ididsocP ==1 & ididantiP ==0 & ididHP ==0, "social", if_else(PunishmentSpending==0, "none", "other"))))]
summary(as.factor(SbjEditPGG$myPtype3))

#using dplyr to summarise by scenario and immunity
PTFs3 <- SbjEditPGG %>%
  group_by(Scenario, Immune) %>%
  summarise(Nsocial = sum(myPtype3 == "social"),
            Nanti = sum(myPtype3 == "anti"),
            NHypo = sum(myPtype3 =="Hypo"),
            Nnon = sum(myPtype3 =="none"),
            Nother = sum(myPtype3 =="other"),
            N = n()) %>% 
  mutate(PercHypo = NHypo/(N-Nnon))
PTFs3 #Punishment Type FrequencieS 3 types (social, anti, hypocritical)

others <-SbjEditPGG[myPtype3=="other"] %>%
  select(Scenario, Immune,Contribution, ididHP, ididantiP, ididsocP, SumC, ContPlayer1, ContPlayer2, ContPlayer3, ContPlayer4, MyDeductionofP1, MyDeductionofP2, MyDeductionofP3, MyDeductionofP4) #shows what the unclassified did

#Was hypocritical punishment more common among immune individuals?
Ftesting3 <-
  matrix(c((19+44), (188+18+12+219+47+21), (26+38), (38+11+2+60+10+2)), #these are the values in PTFs summary
         nrow = 2,
         dimnames = list(Normal = c("Hypo", "NonHypo"),
                         Immune = c("Hypo", "NonHypo")))
fisher.test(Ftesting3, alternative = "two.sided") #Significantly more hypocritical punishment among immune


SbjEditPGG[, AltPun := if_else(myPtype3 == "social", 1, 0)][, AntiPun := if_else(myPtype3 == "anti", 1, 0)][, HypoPun := if_else(myPtype3 == "Hypo", 1, 0)][, OtherPun := if_else(myPtype3 == "other", 1, 0)][, NoPun := if_else(myPtype == "none", 1, 0)]

#measuring punishment consistency
PunConsistency <- SbjEditPGG %>%
  group_by(UPID, Scenario, Immune) %>%
  summarise(sumAlt = sum(AltPun),
            sumAnti = sum(AntiPun),
            sumHypo = sum(HypoPun),
            sumOther = sum(OtherPun),
            sumNoPun = sum(NoPun)) %>%
  mutate(onlyAlt = if_else(sumAlt > 0 & sumAnti ==0 & sumHypo==0, 1,0),
         onlyAnti = if_else(sumAnti > 0 & sumAlt ==0 & sumHypo==0,1,0),
         onlyHypo = if_else(sumHypo > 0 & sumAlt ==0 & sumAnti ==0, 1,0),
         neverAlt = if_else(sumAlt ==0 & (sumAnti > 0 | sumHypo > 0), 1, 0),
         Altplus = if_else(sumAlt > 0 & (sumAnti >0 | sumHypo > 0), 1, 0),
         onlyOther = if_else(sumOther > 0 & sumAlt ==0 & sumAnti ==0 & sumHypo ==0, 1,0),
         NeverPun = if_else(sumNoPun ==5, 1, 0)) %>%
  group_by(Scenario, Immune) %>%
  summarise(sumOnlyAlt = sum(onlyAlt),
            sumOnlyAnti = sum(onlyAnti),
            sumOnlyHypo = sum(onlyHypo),
            sumNeverAlt = sum(neverAlt),
            sumAltplus = sum(Altplus),
            sumOnlyOther = sum(onlyOther),
            sumNeverPun = sum(NeverPun)) %>%
  mutate(Npun = sumOnlyAlt + sumNeverAlt + sumAltplus + sumOnlyOther,
         N = Npun + sumNeverPun,
         percentPOnlyAlt = sumOnlyAlt / Npun)
PunConsistency

TotalConsistency <- PunConsistency %>%
  ungroup() %>%
  summarise(NonlyAlt = sum(sumOnlyAlt),
            Npun = sum(Npun))  %>%
  mutate(PercentAlt = NonlyAlt / Npun)
TotalConsistency
#114 of 220 punishers were consistent 'social' punishers, everyone else either did not punish or sometimes punished anti-socially / hypocritically / other (all 3 of these are anti-social or hypocritical in the broad sense)
#114/220
#so only 52% of those that punished were consistently social punishers (they may have only punished once, point is they never punished anti-socially/hypocritically/other).

print("finished here")


#post review, spending on punishment depending on budget
#did anyone punish 6/6 in sessions 1-2?
#1 was Immune Punisher (N = 20, 20 punishers possible)
#2 was Sole Punisher (N = 20, 5 punishers possible)
SbjEdit[SessionNumber <3 & (MyDeductionofP1==6 | MyDeductionofP2==6 | MyDeductionofP3==6 | MyDeductionofP4==6), UPID]
#Yes 8 UPIDs (out of 25 possible punishers) did, a total of 14 times.

#Did spending on punishment increase after increasing budget?
SbjEdit[SessionNumber < 3 & Stage ==3 & Type ==1,  mean(MySumDeductionofOthers)]
SbjEdit[SessionNumber > 2 & Stage ==3 & Type ==1 & Scenario != "Mutual Punishers",  mean(MySumDeductionofOthers)]
# No, it was 2.1 before, and 2.0 afterwards.

#post review, request about differences in variance
Vari <- SbjEdit %>%
  filter(Stage==3 & Period < 3) %>%
  group_by(Scenario, Immune, Period) %>%
  summarise(meanC = mean(Contribution),
            Vari = var(Contribution, na.rm=TRUE))
Vari

VariMP <- SbjEdit %>%
  filter(Stage==3 & Scenario == "Mutual Punishers") %>%
  summarise(meanC = mean(Contribution),
            Vari = var(Contribution, na.rm=TRUE))
VariMP

VariNIP <- SbjEdit %>%
  filter(Stage==3 & Scenario == "Immune Punisher" & Immune ==0) %>%
  summarise(meanC = mean(Contribution),
            Vari = var(Contribution, na.rm=TRUE))
VariNIP

varData <- SbjEdit[Stage==3]
var.test(varData[Immune==0 & Scenario == "Immune Punisher"]$Contribution, varData[Immune==1 & Scenario == "Immune Punisher"]$Contribution)
var.test(varData[Immune==0 & Scenario == "Sole Punisher"]$Contribution, varData[Immune==1 & Scenario == "Sole Punisher"]$Contribution)
var.test(varData[Immune==0 & Scenario != "Mutual Punishers"]$Contribution, varData[Immune==1 & Scenario != "Mutual Punishers"]$Contribution)
var.test(varData[Scenario == "Mutual Punishers"]$Contribution, varData[Immune==0 & Scenario != "Mutual Punishers"]$Contribution)

var.test(varData[Period ==1 & Immune==0 & Scenario == "Immune Punisher"]$Contribution, varData[Period ==1 & Immune==1 & Scenario == "Immune Punisher"]$Contribution)
var.test(varData[Period ==1 & Immune==0 & Scenario == "Sole Punisher"]$Contribution, varData[Period ==1 & Immune==1 & Scenario == "Sole Punisher"]$Contribution)
var.test(varData[Period ==1 & Immune==0 & Scenario != "Mutual Punishers"]$Contribution, varData[Period ==1 & Immune==1 & Scenario != "Mutual Punishers"]$Contribution)
var.test(varData[Period ==1 & Scenario == "Mutual Punishers"]$Contribution, varData[Period ==1 & Immune==0 & Scenario != "Mutual Punishers"]$Contribution)
?var.test()
#No significant or noteworthy differences



