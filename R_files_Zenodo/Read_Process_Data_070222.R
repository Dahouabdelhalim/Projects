#Code to read model predictions from COVID-hub
#Version: July 2, 2022
#This code was run on the FASRC Cluster
#It took 3 hours and 14 minutes on one core, and used 32 GB of memory
#Most (2h 18 minutes) due to the initial data read step
#The code was run on:
#R version 3.5.1
#Packages viridis version 0.5.1 and viridisLite version 0.3.0
#(these two packages are only needed for plotting figure 1)

###### 0. Required Data Downloads ######
#Running this code requires the following data downloads:
#1) COVID_hub_data deposited on Zenodo
#Available at: 
#https://zenodo.org/record/6301718/files/reichlab/covid19-forecast-hub-zenodo-20220227.zip?download=1
#Deposited on February 27, 2022
#Downloaded on: July 01, 2022
#Size: 3.31 GB (zipped); ~30 GB (unzipped)

#2) State population from the U.S. Census Bureau
#Note: this is only required to plot figure 1
#U.S. Census Bureau, 2021
#Annual Population Estimates, Estimated Components of Resident Population Change, and Rates of the Components of Resident Population Change for the United States, States, District of Columbia, and Puerto Rico: April 1, 2020 to July 1, 2021 (NST-EST2021-ALLDATA)
#Available at:
#https://www2.census.gov/programs-surveys/popest/datasets/2020-2021/state/totals/NST-EST2021-alldata.csv
#Accessed on: January 28, 2022
#Last modified: December 21, 2022
#Size: 12 kB

#Loading viridis package 
#This is the only package required, and is only required to plot Figure 1
library("viridis")

###### 1. Changing Data Locations #######
#Running the code requires changing the following data input and output locations to the appropriate folders
setwd("~/Spengler_echoma/Ernani/COVID_Project/RSOS/Revision_July_2022/")

#a) COVID-hub data location:
input.data.location <- "~/Spengler_echoma/Ernani/COVID_Project/RSOS/Revision_July_2022/COVID_hub_data/reichlab-covid19-forecast-hub-17fb2b6/data-processed"
truth.data.location <- "~/Spengler_echoma/Ernani/COVID_Project/RSOS/Revision_July_2022/COVID_hub_data/reichlab-covid19-forecast-hub-17fb2b6/data-truth/truth-Incident Deaths.csv"

#b) USCB State population data location
state.population.location <- "~/Spengler_echoma/Ernani/COVID_Project/RSOS/Revision_July_2022/State_Populations_USCB/NST-EST2021-alldata.csv"

#c) Output data location, where final processed data and Figure 1 will be stored
location.write.final.data <- "~/Spengler_echoma/Ernani/COVID_Project/RSOS/Revision_July_2022/Final_Expert_Prediction_Data"

###### 2. Reading COVID_hub data #######
#Note: COVID hub data is spread throughout thousands of files across 114 models
#Models were numbered alphabetically
#Some files for models numbered 28, 58, and 77 generated warnings when reading in the initial run
#Aux.i check.files has the expert numbers (folders)
#The warnings are just due to location (coded as factor vs. character) 
#That occurred -- when at a later date models changed locations -- e.g. moved from state to counties
#New county codes don't match old factors for states. R will then assign NA.

#Not all files are "problematic" -- the problematic files are:
#Files numbers 38 to 43 for model 28 
#Files numbers 46 and 47 for model 58
#File number 19 for model 77
#Aux.j.check.files has the respective files for each expert where warnings occurred

#For these three models, we will read these problematic files separately from the rest
#And then merge them

######    2.1 Initial data read ######
aux.s.t <- Sys.time()

setwd(input.data.location)
options(warn=1)
Model.List <- list.dirs()
Model.List <- Model.List[-1]
File.List <- list()
File.List.csv <- list()
aux.variables <- list()

for (i in 1:length(Model.List))
{
  setwd(Model.List[[i]])
  File.List[[i]] <- list.files()
  File.List.csv[[i]] <- list.files(pattern = ".csv")
  if(length(File.List.csv[[i]])>0)
  {
    for(j in 1:length(File.List.csv[[i]]))
    {
      aux.read <- read.csv(File.List.csv[[i]][j])
      if (j==1)
      {
        aux.var.names <- names(aux.read)
        assign(paste0("Exp",i),aux.read)
      }  else 
      {
        aux.bind <- rbind(get(paste0("Exp",i)),aux.read)
        assign(paste0("Exp",i),aux.bind)
        aux.var.names <- c(aux.var.names,names(aux.read))
      }
    }
  } else {
    print(paste0("Model ",i, " has no files"))
  }
  
  aux.variables[[i]] <- aux.var.names
  setwd(input.data.location)
  print(i)
  
}
aux.e.t <- Sys.time()
aux.e.t-aux.s.t

######    2.2. Re-Reading files that created warnings ######
######       2.2.1 Reading only those files with issues ######
aux.j.check.files <- list(c(38:43),c(46:47),c(19))
aux.i.check.files <- c(28,58,77)
# Check.Files <- rep("a",9)
# cf <-1
for (i in aux.i.check.files)
{
  aux.j <- which(i==aux.i.check.files)
  setwd(Model.List[[i]])
  File.List[[i]] <- list.files()
  File.List.csv[[i]] <- list.files(pattern = ".csv")
  if(length(File.List.csv[[i]])>0)
  {
    for(j in aux.j.check.files[[aux.j]])
    {
      print(paste0("aux.j.checkfiles[[aux.j]] = ",aux.j.check.files[[aux.j]][1]))
      aux.read <- read.csv(File.List.csv[[i]][j])
      if (j==aux.j.check.files[[aux.j]][1])
      {
        aux.var.names <- names(aux.read)
        assign(paste0("ExpB",i),aux.read)
      }  else
      {
        aux.bind <- rbind(get(paste0("ExpB",i)),aux.read)
        assign(paste0("ExpB",i),aux.bind)
        aux.var.names <- c(aux.var.names,names(aux.read))
      }
      print(Model.List[[i]])
      print(File.List.csv[[i]][j])
      # Check.Files[cf] <- as.character(paste0(Model.List[i],"/",File.List.csv[[i]][j]))
      print(paste0("j = ",j))
      # cf <- cf+1
    }
  } else {
    print(paste0("Model ",i, " has no files"))
  }
  
  aux.variables[[i]] <- aux.var.names
  setwd(input.data.location)
  print(paste0("i = ",i))
}

######       2.2.2 Reading the remainder files for models with files that generated warnings ######
aux.variables.r <- list()
for (i in aux.i.check.files)
{
  aux.j <- which(i==aux.i.check.files)
  setwd(Model.List[[i]])
  File.List[[i]] <- list.files()
  File.List.csv[[i]] <- list.files(pattern = ".csv")
  
  aux.j2 <- aux.j.check.files[[aux.j]]
  aux.j3 <- 1:length(File.List.csv[[i]])
  aux.j3 <- aux.j3[!aux.j3 %in% aux.j2]
  
  if(length(aux.j3)>0)
  {
    for(j in aux.j3)
    {
      aux.read <- read.csv(File.List.csv[[i]][j])
      if (j==aux.j3[1])
      {
        aux.var.names <- names(aux.read)
        assign(paste0("ExpR",i),aux.read)
      }  else 
      {
        aux.bind <- rbind(get(paste0("ExpR",i)),aux.read)
        assign(paste0("ExpR",i),aux.bind)
        aux.var.names <- c(aux.var.names,names(aux.read))
      }
      print(paste0("j = ",j))
    }
  } else {
    print(paste0("Model ",i, " has no other files"))
  }
  
  aux.variables.r[[i]] <- aux.var.names
  setwd(input.data.location)
  print(paste0("i = ",i))
}

######       2.2.3. Merging files with and without issues ######
#Making location character
ExpB28$location <- as.character(ExpB28$location)
ExpR28$location <- as.character(ExpR28$location)
ExpB58$location <- as.character(ExpB58$location)
ExpR58$location <- as.character(ExpR58$location)
ExpB77$location <- as.character(ExpB77$location)
ExpR77$location <- as.character(ExpR77$location)

#Merging
Exp28_Corrected <- rbind(ExpB28,ExpR28)
Exp58_Corrected <- rbind(ExpB58,ExpR58)
Exp77_Corrected <- rbind(ExpB77,ExpR77)

#Replacing old with new reads/files for the 3 experts with problematic files
Exp28 <- Exp28_Corrected
Exp58 <- Exp58_Corrected
Exp77 <- Exp77_Corrected

aux.e.t.2 <- Sys.time()
aux.e.t.2-aux.s.t

###### 3. Adjusting some variable types ######
#A few models had variables assigned inconsistent types during the read
#i)Model 46 had some NA's coded as character/string "NA"
#ii) Models 12, 63, 68, 70, and 105 had values (predictions) as "character"/"string"
#For this second case, it was due to retracted forecasts. These are coded as "NULL" (character)
#We transform values into numeric
#in order to do so, we assign the value ("code") -987.789 for these "NULL"/retracted forecasts
#(This value does not appear anywhere else in the dataset)

#i) Adjusting model 46
aux.filter <- which((!substr(Exp46$quantile,1,1) %in% 0:9) & !is.na(Exp46$quantile))
Exp46$quantile[aux.filter] <- NA

#ii) Adjusting Models 12, 63, 68, 70, and 105
Exp12$value[Exp12$value=="NULL"] <- as.character(-987.789)
Exp63$value[Exp63$value=="NULL"] <- as.character(-987.789)
Exp68$value[Exp68$value=="NULL"] <- as.character(-987.789)
Exp70$value[Exp70$value=="NULL"] <- as.character(-987.789)
Exp105$value[Exp105$value=="NULL"] <- as.character(-987.789)

Exp12$value <- as.numeric(as.character(Exp12$value))
Exp63$value <- as.numeric(as.character(Exp63$value))
Exp68$value <- as.numeric(as.character(Exp68$value))
Exp70$value <- as.numeric(as.character(Exp70$value))
Exp105$value <- as.numeric(as.character(Exp105$value))

###### 4. Filtering and merging data ######

######    4.1. Filtering to incidental deaths 4 weeks ahead and merging files for all models ######
aux.dims <- matrix(rep(0,2*114),ncol=2)
aux.dims[1,1] <- dim(Exp1)[1]
aux.dims[1,2] <- dim(Exp1)[2]+1

Exp1.4wk <- Exp1[Exp1$target=="4 wk ahead inc death",]

All.4wk.data <- Exp1.4wk
All.4wk.data$ExpID <- as.numeric(1)

#merging
for (aux.e in c(2:40,42:114))
{
  print(aux.e)
  aux.var <- get(paste0("Exp",aux.e))
  test.4wk <- length(aux.var$target[aux.var$target=="4 wk ahead inc death"])
  if (test.4wk > 0)
  {
    aux.var <- aux.var[aux.var$target=="4 wk ahead inc death",]
    aux.var$forecast_date <- as.character(aux.var$forecast_date)
    aux.var$target <- as.character(aux.var$target)
    aux.var$target_end_date <- as.character(aux.var$target_end_date)
    aux.var$quantile <- as.numeric(as.character(aux.var$quantile))
    aux.var$type <- as.character(aux.var$type)
    aux.var$location <- as.character(aux.var$location)
    aux.var$value <- as.numeric(aux.var$value)
    aux.var$ExpID <- as.numeric(aux.e)
    All.4wk.data <- rbind(All.4wk.data,aux.var)
    aux.dims[aux.e,1] <- dim(aux.var)[1]
    aux.dims[aux.e,2] <- dim(aux.var)[2]
  } else {
    aux.dims[aux.e,1] <- 0
    aux.dims[aux.e,2] <- dim(aux.var)[2]
  }
  rm(list=(paste0("Exp",aux.e)))
}

######    4.2. Assigning final model IDs ######
#Changing some model IDs
#We initially assigned model IDs in alphabetical order
#Two models (currently numbered 65 and 73), however, were added very recently. 
#We'll move their assigned numbers to 113 and 114 ('end' of dataset)
#These models' forecasts are not ultimately included since they do not have data for 2021

All.4wk.data$ExpID[All.4wk.data$ExpID==65] <- 115
All.4wk.data$ExpID[All.4wk.data$ExpID==73] <- 116
aux.which <- which(All.4wk.data$ExpID %in% 66:72)
All.4wk.data$ExpID[aux.which] <- All.4wk.data$ExpID[aux.which] - 1
aux.which <- which(All.4wk.data$ExpID %in% 74:114)
All.4wk.data$ExpID[aux.which] <- All.4wk.data$ExpID[aux.which] - 2
All.4wk.data$ExpID[All.4wk.data$ExpID==115] <- 113
All.4wk.data$ExpID[All.4wk.data$ExpID==116] <- 114

Changed.Model.List <- Model.List[c(1:64,66:72,74:114,65,73)]

######    4.3. Filtering data to select experts who predicted in all months of 2021 ######
All.4wk.data$forecast_date <- as.character(All.4wk.data$forecast_date)
All.4wk.data$target <- as.character(All.4wk.data$target)
All.4wk.data$target_end_date <- as.character(All.4wk.data$target_end_date)
All.4wk.data$type <- as.character(All.4wk.data$type)
All.4wk.data$location <- as.character(All.4wk.data$location)

All.4wk.data$target_end_year <- substr(All.4wk.data$target_end_date,1,4)
All.4wk.data$target_end_month <- substr(All.4wk.data$target_end_date,6,7)
All.4wk.data$target_end_day <- substr(All.4wk.data$target_end_date,9,10)

All.4wk.data <- All.4wk.data[All.4wk.data$target_end_year=="2021",]

All.4wk.data$location[All.4wk.data$location=="1"] <- "01"
All.4wk.data$location[All.4wk.data$location=="2"] <- "02"
All.4wk.data$location[All.4wk.data$location=="3"] <- "03"
All.4wk.data$location[All.4wk.data$location=="4"] <- "04"
All.4wk.data$location[All.4wk.data$location=="5"] <- "05"
All.4wk.data$location[All.4wk.data$location=="6"] <- "06"
All.4wk.data$location[All.4wk.data$location=="7"] <- "07"
All.4wk.data$location[All.4wk.data$location=="8"] <- "08"
All.4wk.data$location[All.4wk.data$location=="9"] <- "09"

aux.quantile <- which(All.4wk.data$quantile %in% c(0.05,0.25,0.5,0.75,0.95))
Five.Quantile <- All.4wk.data[aux.quantile,]

#Filtering experts with forecasts for all months of 2021
temporal.coverage <- as.matrix(table(Five.Quantile$ExpID,Five.Quantile$target_end_month))
Exp.coverage <- rowSums(temporal.coverage>0)
Exp.entire2021 <- rownames(temporal.coverage)[which(Exp.coverage==12)] #Total of 27

Five.Quantile <- Five.Quantile[Five.Quantile$ExpID %in% Exp.entire2021,]

######    4.4. Addressing duplicate predictions ######
#4,160 predictions in the dataset are duplicates (two predictions for the same model & target_end_date)
#They are all for 4 versions of the CU model ('CU-nochange', 'CU-scenario_low', 'CU-scenario_mid', and 'CU-select')
#All these predictions are from December 2020, with a target end date of January 2021
#In all cases, the predictions were done twice, 3 days apart
#We will use the later prediction and discard the earlier one

Five.Quantile$ForecastDate <- as.Date(Five.Quantile$forecast_date)
Latest.forecast <- aggregate(Five.Quantile$ForecastDate, by=list(ExpID=Five.Quantile$ExpID,
                                                                 target_end_date=Five.Quantile$target_end_date,
                                                                 location=Five.Quantile$location,
                                                                 quantile=Five.Quantile$quantile),
                             FUN=max)
names(Latest.forecast)[which(names(Latest.forecast)=="x")] <- "Latest_Pred"

Five.Quantile <- merge(x=Five.Quantile,y=Latest.forecast,by=c("ExpID",
                                                              "target_end_date",
                                                              "location",
                                                              "quantile"),
                       all=TRUE)

aux.dup <- which(Five.Quantile$forecast_date != as.character(Five.Quantile$Latest_Pred))

Five.Quantile <- Five.Quantile[-aux.dup,]

#Also removing the 5 retracted forecasts
#These were just 5 forecasts by Model 103 ('USC-SI_kJalpha'), For Florida, in July 2021 
#i.e., forecast date = 06/13/21, target end date = 07/10/21
Five.Quantile <- Five.Quantile[-which(Five.Quantile$value=="-987.789"),] 

######    4.5. Restricting data to the 50 U.S. states, D.C., and National-level forecasts ######
#We will remove the following locations:
#FIPS code 60: American Samoa
#FIPS code 66: Guam
#FIPS code 69: Northern Mariana Islands
#FIPS code 72: Puerto Rico
#FIPS code 78: Virgin Islands
Five.Quantile <- Five.Quantile[!Five.Quantile$location %in% c("60","66","69","72","78"),]

###### 5. Checking overlap of model forecasts -- finding common subset of predictions ######
#We need a common set of questions/forecasts (i.e., a combination of location and date) for all models
#e.g. we can only include forecasts for the state Florida in week 10, if all models forecasted for that state and date (week)
#There is a tradeoff between models included and questions included
#More questions means less common models with predictions, more models means smaller subset of common questions

######    5.1. Creating Question List ######
#Each combination of week and state is a question
#i.e., predicting incident deaths in Massachusetts in Week 1 of our dataset (2021)
#note that we analyze predictions of indicental deaths carried out 4 weeks ahead of time
Question.List <- subset(Five.Quantile, select=c("target_end_date","location"))
Question.List <- unique(Question.List)
Question.List <- as.data.frame(Question.List)
Question.List <- Question.List[order(Question.List$target_end_date,Question.List$location,decreasing=FALSE),]
Question.List$LabelQuestion <- 1:length(Question.List$location)

######    5.2. Creating Final Dataset ######
#Long version
Five.Quantile.Long <- merge(x=Five.Quantile,y=Question.List,by=c("location","target_end_date"),all=TRUE)

#Wide version
aux.wide <- subset(Five.Quantile.Long, select=c("LabelQuestion","ExpID","quantile","value"))

#Reshaping to create Wide dataset
Five.Quantile.Wide <- reshape(aux.wide,
                              timevar="quantile",
                              idvar=c("ExpID","LabelQuestion"),
                              direction="wide")

######    5.3. Filtering forecasts by expert ######
Exp.List <- unique(Five.Quantile.Wide$ExpID)
Exp.List <- Exp.List[order(Exp.List,decreasing=FALSE)]

Miss.Q <- list() #Questions that are missing for each expert (i.e., no predictions)
Have.Q <- list() #Questions for which experts made predictions

for (i in 1:length(Exp.List))
{
  aux.var <- Five.Quantile.Wide[Five.Quantile.Wide$ExpID==Exp.List[i],]
  aux.qu <- unique(aux.var$LabelQuestion)
  aux.qt <- aux.var$LabelQuestion
  aux.qu <- aux.qu[order(aux.qu,decreasing=FALSE)]
  aux.qt <- aux.qt[order(aux.qt,decreasing=FALSE)]
  if(sum(aux.qu!=aux.qt)>0)
  {
    print(paste0("Error -- Model ",i," has duplicated predictions"))
  }
  aux.qa <- Question.List$LabelQuestion
  aux.qm <- aux.qa[!aux.qa %in% aux.qu]
  if(length(aux.qm)>0)
  {
    Miss.Q[[i]] <- aux.qm
  } else {
    print(paste0("Model ",i, "has complete data -- predictions for all questions"))
    Miss.Q[[i]] <- 0
  }
  if(length(aux.qu)>0)
  {
    Have.Q[[i]] <- aux.qu
  } else {
    print(paste0("Error -- Model ",i, "has no data, i.e., no predictions"))
    Have.Q[[i]] <- 0
  }
}

lengths(Miss.Q) #Missing questions by expert -- 1 indicates complete data as I add a "0" there in this case (length=1 element)
lengths(Have.Q) #Predictions by expert (2,704 = complete)

######    5.4. Finding common subsets of questions ######
#Removing the complete experts from this part of the analysis
#These experts have predictions for all questions so they will always be included
#i.e., regardless of what the subset of questions selected is
Exp.List.Analysis <- Exp.List[-which(Miss.Q=="0")]
Exp.List.Complete <- Exp.List[which(Miss.Q=="0")]

aux.df.exp <- as.data.frame(t(Exp.List.Analysis))
names(aux.df.exp) <- Exp.List.Analysis
aux.df.exp[1,] <- 0
aux.df.exp[2,] <- 1

Expert.Combinations <- expand.grid(aux.df.exp)
names(Expert.Combinations)
sum(duplicated(Expert.Combinations))

NonCommon.Q <- list() #This is the list of Non-common questions for each combination of models
# Because these questions are not common, they would be excluded from the analysis
# If the given combination/subset of experts is selected
# With 19 experts/models, there are 2^19 or 524,288 possible combinations of models to check

aux.nlines <- dim(Expert.Combinations)[1] 

Have.Q.Analysis <- Have.Q[-which(Miss.Q=="0")]
aux.s.t <- Sys.time()
#Going through the 524,288 combinations
for (i in 2:aux.nlines)
{
  aux.row <- which(Expert.Combinations[i,]>0)
  aux.c <- Reduce(intersect,Have.Q.Analysis[aux.row])
  aux.t <- Question.List$LabelQuestion
  aux.i <- aux.t[!aux.t %in% aux.c]
  NonCommon.Q[[i]] <- aux.i
  if(i%%10000==0)
  {
    print(i)
  }
}
aux.e.t <- Sys.time()
aux.e.t-aux.s.t
#This step took less than 10 minutes 

###### 6. Reading Truth Data (i.e., realizations) ######
#We take as week 1 the epidemiological week ending on Saturday Jan 02, 2021.
#This week therefore starts on 12/27/2020, and we call this week 1
#Note that the CDC uses a slightly different Epidemiological week number
#Our week 1 might be included in year 2020 (i.e., the last week of 2020), and not 2021

#List of states
location.list <- unique(Question.List$location)
location.list <- location.list[order(location.list,decreasing=FALSE)]

Truth_IncDeath <- read.csv(truth.data.location)

#The "truth" data are provided as incident deaths per day, so we aggregate to obtain weekly counts
Truth_IncDeath <- Truth_IncDeath[Truth_IncDeath$location %in% location.list,]
Truth_IncDeath$newdate <- as.Date(Truth_IncDeath$date)
Truth_IncDeath <- Truth_IncDeath[Truth_IncDeath$newdate>=as.Date("2020-12-27"),]
Truth_IncDeath$diffdate <- Truth_IncDeath$newdate-as.Date("2020-12-27")
Truth_IncDeath$diffdate <- as.numeric(Truth_IncDeath$diffdate)
Truth_IncDeath$Week <- (Truth_IncDeath$diffdate %/% 7)+1
summary(Truth_IncDeath$Week)
Truth_IncDeath <- Truth_IncDeath[Truth_IncDeath$Week<=52,]

Truth_Weekly_IncDeaths <- aggregate(Truth_IncDeath$value,by=list(Week=Truth_IncDeath$Week,
                                                                 location=Truth_IncDeath$location),
                                    FUN=sum)
names(Truth_Weekly_IncDeaths)[which(names(Truth_Weekly_IncDeaths)=="x")] <- "Weekly_Inc_Deaths"
names(Truth_IncDeath)

Truth_Weekly_IncDeaths$target_end_date <- as.Date("2021-01-02")
Truth_Weekly_IncDeaths$target_end_date <- as.Date("2021-01-02") + 7*(Truth_Weekly_IncDeaths$Week-1)

Truth_Weekly_IncDeaths$location <- as.character(Truth_Weekly_IncDeaths$location)

###### 7. Excluding Questions with non-positive realizations ######
#We excluding questions where non-positive deaths were reported
#i.e., either 0 deaths or negative deaths reported as truth
#From 2,704 questions, i.e. 52 states x 52 weeks, we have 2,666 with positive realizations
#i.e., we exclude 38 of the 2,704 questions
Question.List$target_end_date <- as.Date(Question.List$target_end_date)
Question.List <- merge(x=Question.List, y=Truth_Weekly_IncDeaths, by=c("location","target_end_date"), all=TRUE)

Question.List.Complete <- Question.List
Question.List <- Question.List[Question.List$Weekly_Inc_Deaths>0,]


###### 8. Choosing the common subset of questions ######
#8 Models (Experts) provide complete data, i.e., forecasts for all states and weeks of interest
#That leaves 19 experts with non-complete data
#Therefore:
#Including all questions requires restricting the dataset to 8 experts
#Including all experts requires excluding data
#There is a tradeoff between more experts or more questions included

#The total number of questions is 2,666, after non-positive realizations were reported
#Checking the maximum number of questions that can be included (i.e., that is common), 
#as a function of the number of models included
#If models <= 8 we can choose the 8 models with complete data
#and work with a complete dataset (2,666 questions)
#If we want to include more than 8 models, we check
#what combination of models has the maximum common set of questions

NonCommon.Q.Filtered <- list()
for (i in 1:length(NonCommon.Q))
{
  NonCommon.Q.Filtered[[i]] <- NonCommon.Q[[i]][NonCommon.Q[[i]] %in% Question.List$LabelQuestion]
  if(i%%10000==0)
  {
    print(i)
  }
}
Common.Experts <- list()
for (i in 1:19)
{
  Common.Experts[[i]] <- which(rowSums(Expert.Combinations)==i)
  Excluded.Questions <- min(lengths(NonCommon.Q[Common.Experts[[i]]]))
  print(paste0("Maximum number of common questions if number of Experts is ",i+8, " = ",2666-Excluded.Questions))
  
}

######    8.1. Checking spatial and temporal coverage ######
#We also want to check what is the spatial and temporal coverage of the commons sets of questions
Coverage.by.State <- matrix(rep(-1,524288*52),ncol=52)
Coverage.by.Week <- matrix(rep(-1,524288*52),ncol=52)
colnames(Coverage.by.State) <- unique(Question.List$location)

aux.mg <- as.data.frame(0:51)
names(aux.mg) <- "aux.loc"

Qcount <- as.data.frame(Question.List$LabelQuestion)
names(Qcount) <- "aux.q"
Qcount$aux.loc <- Qcount$aux.q%%52
Qcount$aux.temp <- ((Qcount$aux.q-1) %/% 52)
Full.Q.State.count <- aggregate(Qcount$aux.q,by=list(Qcount$aux.loc),FUN=length)
names(Full.Q.State.count) <- c("aux.loc","count")
Full.Q.State.count <- Full.Q.State.count[order(Full.Q.State.count$aux.loc,decreasing=FALSE),]
Full.Q.Week.count <- aggregate(Qcount$aux.q,by=list(Qcount$aux.temp),FUN=length)
names(Full.Q.Week.count) <- c("aux.temp","count")
Full.Q.Week.count <- Full.Q.Week.count[order(Full.Q.Week.count$aux.temp,decreasing=FALSE),]

for (i in 1:19)
{
  aux.list <- Common.Experts[[i]]
  for(j in 1:length(aux.list))
  {
    aux.row <- aux.list[j]
    aux.qs <- as.data.frame(NonCommon.Q.Filtered[[aux.row]])
    names(aux.qs) <- "aux.var"
    aux.qs$aux.loc <- aux.qs$aux.var %% 52
    aux.qs$aux.temp <- ((aux.qs$aux.var-1) %/% 52)
    
    aux.count.state <- aggregate(aux.qs$aux.var,by=list(aux.qs$aux.loc),FUN=length)
    names(aux.count.state) <- c("aux.loc","count.state")
    aux.count.state <- merge(x=aux.mg,y=aux.count.state,by="aux.loc",all=TRUE)
    names(aux.count.state)
    aux.count.state <- aux.count.state[order(aux.count.state$aux.loc,decreasing=FALSE),]
    aux.count.state$count.state[is.na(aux.count.state$count.state)] <- 0
    
    aux.count.week <- aggregate(aux.qs$aux.var,by=list(aux.qs$aux.temp),FUN=length)
    names(aux.count.week) <- c("aux.loc","count.week")
    aux.count.week <- merge(x=aux.mg,y=aux.count.week,by="aux.loc",all=TRUE)
    names(aux.count.week)
    aux.count.week <- aux.count.week[order(aux.count.week$aux.loc,decreasing=FALSE),]
    aux.count.week$count.week[is.na(aux.count.week$count.week)] <- 0
    
    Coverage.by.State[aux.row,] <- Full.Q.State.count$count - aux.count.state$count
    Coverage.by.Week[aux.row,] <- Full.Q.Week.count$count - aux.count.week$count
  }
  print(i)
}

#Common weeks by each state (column) for 19+8 = 27 models included
#By state
Coverage.by.State[Common.Experts[[19]],]
#Total:
sum(Coverage.by.State[Common.Experts[[19]],])

#Common weeks by each state (column) for 18+8 = 26 models included
#that is, for each possible combination of 26 models
#19 choose 18 = 19 possible combinations
#(we are always including the 8 models with complete data)
#By state
Coverage.by.State[Common.Experts[[18]],]
#Total:
rowSums(Coverage.by.State[Common.Experts[[18]],])

#Common weeks by each state (column) for 17+8 = 25 models included
#that is, for each possible combination of 25 models
#19 choose 17 = 19*18/2 = 171 possible combinations
#(we are always including the 8 models with complete data)
#By state
Coverage.by.State[Common.Experts[[17]],]
#Total:
rowSums(Coverage.by.State[Common.Experts[[17]],])

#Common states by each week (column) for 19+8 = 27 models included
#By week:
Coverage.by.Week[Common.Experts[[19]],]
#Total:
sum(Coverage.by.Week[Common.Experts[[19]],])

#Common weeks by each week (column) for 18+8 = 26 models included
#that is, for each possible combination of 26 models
#19 choose 18 = 19 possible combinations
#(we are always including the 8 models with complete data)
#By week:
Coverage.by.Week[Common.Experts[[18]],]
#Total
rowSums(Coverage.by.Week[Common.Experts[[18]],])

#Common weeks by each week (column) for 17+8 = 25 models included
#that is, for each possible combination of 25 models
#19 choose 17 = 19*18/2 = 171 possible combinations
#(we are always including the 8 models with complete data)
#By week:
Coverage.by.Week[Common.Experts[[17]],]
#Total
rowSums(Coverage.by.Week[Common.Experts[[17]],])

#In order to improve spatial and temporal coverage, we will remove two infrequently cited models
#These models are "Microsoft-DeepSTIA" and "UMich-RidgeTfReg" 
#We will therefore keep 25 of the 27 models
#Our chosen set is 
Common.Experts[[17]][25]
Expert.Combinations[Common.Experts[[17]][25],]
#plus the 8 experts with forecasts for all questions (i.e., complete data). Those were:
Exp.List.Complete

#The final list of questions is
Excluded.Questions <- NonCommon.Q.Filtered[Common.Experts[[17]][25]]
Excluded.Questions <- as.vector(Excluded.Questions[[1]])
Final.Question.List <- Question.List[!(Question.List$LabelQuestion %in% Excluded.Questions),]

#Our set of 25 includes three ensemble models:
substr(Changed.Model.List[c(12,13,15)],3,nchar(Changed.Model.List[c(12,13,15)]))
#Since they had identical forecasts through November 2021, we will keep only one: 
#Later in the analysis we choose to keep expert 12 only, i.e. "COVIDhub_CDC-ensemble", 
#That is, we removed experts 13 and 15
#This does not affect spatial and temporal coverage since all three ensembles had complete data
#That is, they forecasted all weeks and locations/states included.

#Experts (Models) Included:
Experts.Included <- as.numeric(colnames(Expert.Combinations)[which(Expert.Combinations[Common.Experts[[17]][25],]==1)])
Experts.Included <- c(Experts.Included,Exp.List.Complete)
Experts.Included <- Experts.Included[order(Experts.Included,decreasing=FALSE)]

Names.Experts.Included <- Changed.Model.List[Experts.Included]
Names.Experts.Included <- substr(Names.Experts.Included,3,nchar(Names.Experts.Included))

#Experts (Models) Excluded
Experts.Excluded <- as.numeric(colnames(Expert.Combinations)[which(Expert.Combinations[Common.Experts[[17]][25],]==0)])
Experts.Excluded <- Experts.Excluded[order(Experts.Excluded,decreasing=FALSE)]

Names.Experts.Excluded <- Changed.Model.List[Experts.Excluded]
Names.Experts.Excluded <- substr(Names.Experts.Excluded,3,nchar(Names.Experts.Excluded))

###### 9. Exporting final dataset ######
Final.Expert.List <- as.data.frame(Experts.Included)
names(Final.Expert.List) <- "Expert_ID"
Final.Expert.List$Expert_Name <- substr(Changed.Model.List[Final.Expert.List$Expert_ID],3,nchar(Changed.Model.List[Final.Expert.List$Expert_ID]))

Location.List <- c("Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","District of Columbia","Florida",
                   "Georgia","Hawaii","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine",
                   "Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire",
                   "New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island",
                   "South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virginia","Washington","West Virginia","Wisconsin",
                   "Wyoming","United States")
Location.FIPS <- unique(Question.List$location)

Locations <- as.data.frame(Location.List)
names(Locations) <- "location_name"
Locations$FIPS_Code <- Location.FIPS

Final.Question.List <- merge(x=Final.Question.List, y=Locations, by.x="location", by.y="FIPS_Code", all.x=TRUE, all.y=FALSE)
Final.Question.List <- subset(Final.Question.List, select=c("LabelQuestion",
                                                            "target_end_date",
                                                            "location",
                                                            "location_name",
                                                            "Weekly_Inc_Deaths"))
names(Final.Question.List)[5] <- "Realization"

Final.Question.List <- Final.Question.List[order(Final.Question.List$LabelQuestion, decreasing=FALSE),]
Final.Realizations <- subset(Final.Question.List, select=c("LabelQuestion","Realization"))
Final.Expert.Predictions <- Five.Quantile.Wide[(Five.Quantile.Wide$LabelQuestion %in% Final.Question.List$LabelQuestion) & (Five.Quantile.Wide$ExpID %in% Final.Expert.List$Expert_ID),]

names(Final.Expert.Predictions) <- c("LabelQuestion","Expert_ID",
                                     paste0("Percent_",c(5,25,50,75,95)))

#Writing Question List, Expert List, and Realizations
write.csv(Final.Question.List, file=paste0(location.write.final.data,"/Question_List.csv"), row.names=FALSE)
write.csv(Final.Expert.List, file=paste0(location.write.final.data,"/Expert_List.csv"), row.names=FALSE)
write.csv(Final.Realizations, file=paste0(location.write.final.data,"/realizations.csv"), row.names=FALSE)

#Writing Expert predictions (one file for each expert)
for (i in 1:length(Final.Expert.List$Expert_ID))
{
  expert.selection <- Final.Expert.List$Expert_ID[i]
  predictions.selection <- Final.Expert.Predictions[Final.Expert.Predictions$Expert_ID==expert.selection,]
  filtered.question.list <- predictions.selection$LabelQuestion
  filtered.question.list <- filtered.question.list[order(filtered.question.list, decreasing=FALSE)]
  correct.question.list <- Final.Question.List$LabelQuestion
  correct.question.list <- correct.question.list[order(correct.question.list, decreasing=FALSE)]
  check.different <- sum(filtered.question.list != correct.question.list)
  check.same <- sum(filtered.question.list == correct.question.list)
  
  if(check.different != 0)
  {
    print(paste0("Error -- Question List for Expert ", expert.selection, " is not correct"))
  }
  if(check.same == length(Final.Question.List$LabelQuestion))
  {
    print(paste0("Question List for Expert ", expert.selection, " is right"))
  }
  predictions.selection <- subset(predictions.selection, select=c("Expert_ID","LabelQuestion",
                                                                  paste0("Percent_",c(5,25,50,75,95))))
  predictions.selection <- predictions.selection[order(predictions.selection$LabelQuestion, decreasing=FALSE),]
  write.csv(predictions.selection, file=paste0(location.write.final.data,"/Exp_",expert.selection,".csv"),row.names=FALSE)
}

###### 10. Plotting Figure 1 ######
######    10.1 reading and filtering State population ######
State.Pop <- read.csv(state.population.location)
State.Pop <- subset(State.Pop, select=c("STATE","NAME","POPESTIMATE2021"))
State.Pop <- State.Pop[!State.Pop$NAME %in% c("Northeast Region",
                                              "Midwest Region",
                                              "South Region",
                                              "West Region",
                                              "Puerto Rico"),]
State.Pop$STATE <- as.character(State.Pop$STATE)
State.Pop$NAME <- as.character(State.Pop$NAME)
State.Pop$STATE[State.Pop$STATE=="0"] <- "US"
State.Pop$STATE[State.Pop$STATE=="1"] <- "01"
State.Pop$STATE[State.Pop$STATE=="2"] <- "02"
State.Pop$STATE[State.Pop$STATE=="4"] <- "04"
State.Pop$STATE[State.Pop$STATE=="5"] <- "05"
State.Pop$STATE[State.Pop$STATE=="6"] <- "06"
State.Pop$STATE[State.Pop$STATE=="8"] <- "08"
State.Pop$STATE[State.Pop$STATE=="9"] <- "09"

Locations$Code <- c("AL","AK","AZ","AR","CA","CO","CT","DE","DC","FL","GA","HI","ID","IL","IN","IA","KS","KY","LA","ME",
                    "MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI",
                    "SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY","US")

names(Locations)
State.Pop <- merge(x=State.Pop, y=subset(Locations, select=c("FIPS_Code","Code")),
                   by.x="STATE", by.y="FIPS_Code", all=TRUE)

Order.States <- State.Pop$STATE[order(State.Pop$POPESTIMATE2021,decreasing=TRUE)]

######    10.2. Creating plotting dataset ######
Question.List.Complete <- merge(x=Question.List.Complete, y=State.Pop, 
                                by.x="location", by.y="STATE", all=TRUE)

Question.List.Complete.Ordered <- Question.List.Complete[order(Question.List.Complete$POPESTIMATE2021,Question.List.Complete$LabelQuestion,decreasing=FALSE),]
Question.List.Complete.Ordered$Weekly_Death_Rate <- Question.List.Complete.Ordered$Weekly_Inc_Deaths/Question.List.Complete.Ordered$POPESTIMATE2021

Question.List.Complete.Ordered$Log_Daily_Death_Rate <- log10(Question.List.Complete.Ordered$Weekly_Death_Rate/7)
Question.List.Complete.Ordered$Log_Daily_Death_Rate[is.na(Question.List.Complete.Ordered$Log_Daily_Death_Rate)]

aux.cut <- which(Question.List.Complete.Ordered$Weekly_Inc_Deaths<=0)
range(Question.List.Complete.Ordered$Log_Daily_Death_Rate[-aux.cut])
Min.Rate <- min(Question.List.Complete.Ordered$Log_Daily_Death_Rate[-aux.cut])
Max.Rate <- max(Question.List.Complete.Ordered$Log_Daily_Death_Rate[-aux.cut])

#Creating color gradient
N.Color.Cuts <- 10000
Gray.Shade <- 0.55
Color.Scale <- viridis(N.Color.Cuts,alpha=1,begin=0,end=1,direction=1)

Color.Gap <- (Max.Rate-Min.Rate)/(N.Color.Cuts)
Question.List.Complete.Ordered$Color.Shade <- ceiling(((Question.List.Complete.Ordered$Log_Daily_Death_Rate-Min.Rate)/Color.Gap))
Question.List.Complete.Ordered$Color.Shade[aux.cut] <- -1
Question.List.Complete.Ordered$Color.Shade[Question.List.Complete.Ordered$Color.Shade==0] <- 1
range(Question.List.Complete.Ordered$Color.Shade[-aux.cut])

Question.List.Complete.Ordered$Plot.Color <- 0
Question.List.Complete.Ordered$Plot.Color[-aux.cut] <- Color.Scale[Question.List.Complete.Ordered$Color.Shade[-aux.cut]]
Question.List.Complete.Ordered$Plot.Color[aux.cut] <- gray(Gray.Shade)

Question.List.Complete.Ordered$Order.Index <- 1:2704
rownames(Question.List.Complete.Ordered) <- 1:2704

#Creating plot points, i.e. rectangles
aux.rec.width <- 1
aux.rec.height <- aux.rec.width

Question.List.Complete.Ordered$Rec.xleft <- (Question.List.Complete.Ordered$Order.Index-1)%%52
Question.List.Complete.Ordered$Rec.xright <- Question.List.Complete.Ordered$Rec.xleft + aux.rec.width
Question.List.Complete.Ordered$Rec.ybottom <- (Question.List.Complete.Ordered$Order.Index-1)%/%52
Question.List.Complete.Ordered$Rec.ytop <- Question.List.Complete.Ordered$Rec.ybottom + aux.rec.height

Location.Legend <- State.Pop$Code[order(State.Pop$POPESTIMATE2021,decreasing=FALSE)]

Data.Select <- which(!Question.List.Complete.Ordered$LabelQuestion %in% Excluded.Questions)

######    10.3. Plotting ######
#Plotting options
Week.Legend <- 1:52
x.tick.width <- 0.05
y.tick.width <- 0.5
col.ticks <- "black"
col.labels <- "black"
main.y.offset <- 0
main.x.offset <- 0
main.x.title.offset <- -3
main.y.title.offset <- 2
plot.ylim <- c(-5,main.y.offset+52+10) 
plot.xlim <- c(-8,70) 
aux.ysize <- plot.ylim[2]-plot.ylim[1]
aux.xsize <- plot.xlim[2]-plot.xlim[1]
adj.ratio <- aux.ysize/aux.xsize

#Plot legend options
Leg.height <- 52
Leg.width <- 5
Leg.bottom <- 0
Leg.top <- Leg.bottom+Leg.height
Leg.xleft  <- rep((main.x.offset+52+3),10000)
Leg.xright <- Leg.xleft+Leg.width
Leg.ybottom <- Leg.bottom+(Leg.height/10000)*0:9999
Leg.ytop <- Leg.bottom+(Leg.height/10000)*1:10000
Leg.Lab <- c("0.07","0.1","0.3","1","3","10","30","61")
Leg.Lab.M <- "Daily Death Rate [per million]"
aux.pos.2 <- c(10^(Min.Rate),1e-6*c(0.1,0.3,1,3,10,30),10^(Max.Rate))
Leg.Lab.Pos <- Leg.bottom + ((log10(aux.pos.2)-Min.Rate)/(Max.Rate-Min.Rate))*Leg.height

Plot.Titles <- c("Chosen Dataset (N=23)")
Color.Contour.Select <- which(names(Question.List.Complete.Ordered)=="Plot.Color")
pdf.file.path <- paste0(location.write.final.data,"/Figure_1_070122.pdf")

plot(0)
pdf(file=pdf.file.path,width=5,height=5*adj.ratio)
par(mar=c(0.1,0.1,0.1,0.1),xaxs="i",yaxs="i")
plot(y=-100,x=-100,
     ylim=plot.ylim,
     xlim=plot.xlim,
     xaxt="n",
     yaxt="n",
     ylab="",
     xlab="",
     bty="n")

#Contour Plot
rect(xleft=Question.List.Complete.Ordered[Data.Select,]$Rec.xleft+main.x.offset,
     xright=Question.List.Complete.Ordered[Data.Select,]$Rec.xright+main.x.offset,
     ybottom=Question.List.Complete.Ordered[Data.Select,]$Rec.ybottom+main.y.offset,
     ytop=Question.List.Complete.Ordered[Data.Select,]$Rec.ytop+main.y.offset,
     density=NA,
     border=NA,
     col=Question.List.Complete.Ordered[Data.Select,Color.Contour.Select])

#x-Axis for Question plot
#ticks
rect(xleft=(1:52)-0.5-x.tick.width+main.x.offset,
     xright=(1:52)-0.5+x.tick.width+main.x.offset,
     ybottom=52+main.y.offset,
     ytop=52+y.tick.width+main.y.offset,
     density=NA,
     border=NA,
     col=col.ticks)
#labels
text(x=(1:52)-0.5-main.x.offset,
     y=52+y.tick.width+main.y.offset,
     labels=Week.Legend,
     cex=0.5,col=col.labels,
     srt=90, font=2, adj=c(0,0.5))

#y-Axis for Question Plot
#ticks
rect(ybottom=(1:52)-0.5-x.tick.width+main.y.offset,
     ytop=(1:52)-0.5+x.tick.width+main.y.offset,
     xright=0+main.x.offset,
     xleft=0-y.tick.width+main.x.offset,
     density=NA,
     border=NA,
     col=col.ticks)

#labels
text(y=(1:52)-0.5+main.y.offset,
     x=0-y.tick.width+main.x.offset,
     labels=Location.Legend,
     cex=0.5,col=col.labels,
     font=2, adj=c(1,0.5))

#x title for question plot
text(y=26+main.y.offset,
     x=main.x.title.offset+main.x.offset-1.5,
     "States (ordered by population)",
     srt=90,adj=c(0.5,0),
     font=2, cex=0.8, col="black")

#y title for question plot
text(x=26+main.x.offset,
     y=52+main.y.title.offset+main.y.offset+1,
     "Weeks (January-December)",
     adj=c(0.5,0),
     font=2, cex=0.8, col="black")

#Plot titles
text(x=26+main.x.offset,
     y=52+main.y.title.offset+main.y.offset+7.5,
     labels=Plot.Titles,
     adj=c(0.5,1),
     font=2, cex=1, col="black")

#Plot borders
rect(xleft=main.x.offset,
     xright=main.x.offset+52,
     ybottom=main.y.offset,
     ytop=main.y.offset+52,
     border="black",lwd=1)

#Plot.legend
rect(xleft=Leg.xleft,
     xright=Leg.xright,
     ybottom=Leg.ybottom,
     ytop=Leg.ytop,
     density=NA,
     border=NA,
     col=Color.Scale)

rect(xleft=Leg.xleft[1],
     xright=Leg.xright[1],
     ybottom=Leg.ybottom[1],
     ytop=Leg.ytop[10000],
     border="black",
     lwd=1)

rect(xleft=Leg.xright[1],
     xright=Leg.xright[1]+2*y.tick.width,
     ybottom=Leg.Lab.Pos-2*x.tick.width,
     ytop=Leg.Lab.Pos+2*x.tick.width,
     density=NA, border=NA,
     col=col.ticks)

text(x=Leg.xright[1]+2*y.tick.width,
     y=Leg.Lab.Pos,
     labels=Leg.Lab,
     cex=0.8, font=2, col="black",
     adj=c(0,0.5))

text(x=Leg.xright[1]+2*y.tick.width+5,
     y=Leg.bottom+((Leg.top-Leg.bottom)/2),
     labels=Leg.Lab.M,
     cex=1, font=2, col="black",
     srt=90,
     adj=c(0.5,1))

rect(xleft=main.x.offset,
     xright=main.x.offset+1,
     ytop=main.y.offset-2,
     ybottom=main.y.offset-3,
     density=NA,border="black",
     col=gray(Gray.Shade))

text(x=main.x.offset+2,
     y=main.y.offset-2.5,
     labels=bquote("Realizations"<="0 : excluded"),
     cex=0.8,adj=c(0,0.5))
dev.off() 
