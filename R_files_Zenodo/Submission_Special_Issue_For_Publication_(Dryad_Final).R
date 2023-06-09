##############################################################
#                                                            #
#     Title: Developmental Sex Differences in Hyenas         #
#     Author: S. Kevin McCormick                             #
#     Date: 12/05/2022                                       #
#     Purpose: Comparisons of Unsolicited Submissive Acts    #
#                                                            #
#                                                            #
##############################################################


#set working directory
rm(list = ls())
setwd("C:/Users/skmcc/Documents/R/R_wd")
options(stringsAsFactors = FALSE)
library(tidyverse)
library(sjstats)
library(glmmTMB)
set.seed(101)
library(broom)
library(broom.mixed)
library(tidyverse)



#Import unsolicited submissive data
Data_File_Sub1 <- read.csv("SubmissiveCountsSimplifiedDateMom_Dryad.csv")
str(Data_File_Sub1)
tblRanks <- read.csv("tblRanks.csv")

#Restrict it to the clans that have good data
Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$clan == "talek" |
                                   Data_File_Sub1$clan == "serena.n"|
                                   Data_File_Sub1$clan == "serena.s",]

#Since we are looking at rates, we should restrict the duration
#For now we will start with greater than 10 minutes
#Since very long sessions are few and far between,
#and because I'm seeing some entries that look too long,
#I'm going to restrict it to under 2 hours as well (or 120 minutes)

Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$Duration > 10,]
Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$Duration < 120,]

#For now we will also remove natal den location contexts
Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$location != "n",]
Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$location != "g",]

#Need reliable sexes
Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$Actor_Sex != "u",]
Data_File_Sub1$Actor_Sex[Data_File_Sub1$Actor_Sex == "f"] <- "Female"
Data_File_Sub1$Actor_Sex[Data_File_Sub1$Actor_Sex == "m"] <- "Male"

#Something odd in the translation from ACCESS to R in the unknown age class
unique(Data_File_Sub1[Data_File_Sub1$Actor_Age_Class == "Unknown" &
                        Data_File_Sub1$Actor_Sex == "Female",]$Actor_ID)

Data_File_Sub1$Actor_Age_Class <- Data_File_Sub1$Actor_Age_Class2
Data_File_Sub1$Actor_Age_Class[Data_File_Sub1$Actor_Age <= 1] <- "Cub"
Data_File_Sub1$Actor_Age_Class[Data_File_Sub1$Actor_Age >= 1] <- "Sub"
Data_File_Sub1$Actor_Age_Class[Data_File_Sub1$Actor_Age > 2] <- "Adult"
Data_File_Sub1$Actor_Age_Class[is.na(Data_File_Sub1$Actor_Age)] <- "Unknown"

unique(Data_File_Sub1[Data_File_Sub1$Actor_Age_Class == "Unknown" &
                        Data_File_Sub1$Actor_Sex == "Female",]$Actor_ID)

#Need to remove aliens that had missing info
unique(Data_File_Sub1[Data_File_Sub1$Actor_Status == "u" &
                        Data_File_Sub1$Actor_Sex == "Female",]$Actor_ID)
unique(Data_File_Sub1[Data_File_Sub1$Actor_Status == "t" &
                        Data_File_Sub1$Actor_Sex == "Female",]$Actor_ID)
Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$Actor_Status != "u",]
Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$Actor_Status != "t",]

unique(Data_File_Sub1[Data_File_Sub1$Actor_Age_Class == "Unknown" &
                        Data_File_Sub1$Actor_Sex == "Female",]$Actor_ID)

#Make adult females from the beginning of the study that don't have ages adults
Data_File_Sub1$Actor_Age_Class[Data_File_Sub1$Actor_Age_Class == "Unknown" &
                                 Data_File_Sub1$Actor_Sex == "Female"] <- "Adult"


#Make a Male Adult Age class
Data_File_Sub1$Actor_Sex[Data_File_Sub1$Actor_Age > 2 &
                           Data_File_Sub1$Actor_Sex == "Male" &
                           Data_File_Sub1$Actor_Status == "r"] <- "Natal_Male"


Data_File_Sub1$Actor_Age_Class[Data_File_Sub1$Actor_Sex == "Natal_Male"] <- "Adult"

Data_File_Sub1$Actor_Sex[Data_File_Sub1$Actor_Status == "i" &
                           Data_File_Sub1$Actor_Sex == "Male"] <- "Immigrant_Male"

Data_File_Sub1$Actor_Age_Class[Data_File_Sub1$Actor_Sex == "Immigrant_Male"] <- "Adult"

#Still have some unknowns in our actors
#These are relatively few, so I am just removing the entries or blanks

Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$Actor_Age_Class != "Unknown",]
Data_File_Sub1 <- Data_File_Sub1[Data_File_Sub1$Actor_Sex != "",]
str(Data_File_Sub1)

#Subset the data on Actor Age Class
#Reminder: We are interested in general spontaneous submissive behavior of the Actor
#based on the sex of the acting aggressor.

AdultData <- Data_File_Sub1[Data_File_Sub1$Actor_Age_Class == "Adult",]
SubData <- Data_File_Sub1[Data_File_Sub1$Actor_Age_Class == "Sub",]
CubData <- Data_File_Sub1[Data_File_Sub1$Actor_Age_Class == "Cub",]


CubData$Duration <- CubData$Duration/60
SubData$Duration <- SubData$Duration/60
AdultData$Duration <- AdultData$Duration/60


#Adults

AdultData$year <- as.numeric(AdultData$year)
tblRanks$year <- as.numeric(tblRanks$year)
AdultData$id <- AdultData$Actor_ID
str(AdultData)

#Assign known ranks
AdultData <- left_join(AdultData, tblRanks,
                       by = c("id", "year", "clan"), na_matches = "never")
summary(AdultData$stan_rank)     #13583 remaining
str(AdultData)

#Assign mom ranks if needed
AdultData <- left_join(AdultData, tblRanks,
                       by = c("mom" = "id", "year", "clan"), na_matches = "never")
colnames(AdultData)[16] <- "stan_rank"
colnames(AdultData)[17] <- "mom_stan_rank"
str(AdultData)

##Females assigned yearly ranks for each year where they are at least 18 months old at the start of the year
##Add mom's rank for anyone under age 2.5 without a rank yet
AdultData$stan_rank <- ifelse(is.na(AdultData$stan_rank) & !is.na(AdultData$Actor_Age) & (AdultData$Actor_Age < 2.5),
                              AdultData$mom_stan_rank, AdultData$stan_rank)
summary(AdultData$stan_rank)     #11726 remaining

##Males assigned yearly ranks for each year where they are at least 5 years old at the start of the year [or when they are immigrants]
##Add mom's rank for males under age 5 without a rank yet
AdultData$stan_rank <- ifelse(is.na(AdultData$stan_rank) & !is.na(AdultData$Actor_Age) & (AdultData$Actor_Age < 6) & AdultData$Actor_Sex == "Natal_Male",
                              AdultData$mom_stan_rank, AdultData$stan_rank)
summary(AdultData$stan_rank)     #9944 remaining
AdultData$last_year <- NA
for(year.step in 1:4){
  AdultData$last_year <- AdultData$year - year.step
  AdultData <- left_join(AdultData, tblRanks,
                         by = c("mom" = "id", "last_year" = "year", "clan"), na_matches = "never")
  colnames(AdultData)[16] <- "stan_rank"
  colnames(AdultData)[19] <- "mom_stan_rank_last_year"

  if(year.step < 2){
    ##Females assigned yearly ranks for each year where they are at least 18 months old at the start of the year
    ##Add mom's rank for anyone under age 2.5 without a rank yet
    AdultData$stan_rank <- ifelse(is.na(AdultData$stan_rank) & !is.na(AdultData$Actor_Age) & (AdultData$Actor_Age < 2.5),
                                  AdultData$mom_stan_rank_last_year, AdultData$stan_rank)
  }
  if(year.step < 5){
    ##Males assigned yearly ranks for each year where they are at least 5 years old at the start of the year [or when they are immigrants]
    ##Add mom's rank for males under age 5 without a rank yet
    AdultData$stan_rank <- ifelse(is.na(AdultData$stan_rank) & !is.na(AdultData$Actor_Age) & (AdultData$Actor_Age < 6) & AdultData$Actor_Sex == "Natal_Male",
                                  AdultData$mom_stan_rank_last_year, AdultData$stan_rank)
  }
  AdultData <- AdultData[,c(1:19)]
  print(summary(AdultData$stan_rank))
}

#Look at remaining animals
remaining <- filter(AdultData, is.na(stan_rank))
remaining <- unique(remaining[,c("clan", "year", "id")])
remaining <- filter(remaining, !(clan == "talek" & year > 2013))
check <- remaining %>% group_by(clan, id) %>% summarise(num_years = length(year))
AdultData <- AdultData[!is.na(AdultData$stan_rank),]

adult.unsol.appease.zipoiss.RANK <- glmmTMB(Count ~ Actor_Sex + Group_Size + stan_rank +
                                           offset(log(Duration)) +
                                             (1|session),
                                         data = AdultData,
                                         ziformula = ~1,
                                         family = poisson())

summary(adult.unsol.appease.zipoiss.RANK)

#Cubs

CubData$year <- as.numeric(CubData$year)
CubData$id <- CubData$Actor_ID


#Assign known ranks
CubData <- left_join(CubData, tblRanks,
                     by = c("id", "year", "clan"), na_matches = "never")
summary(CubData$stan_rank)     #13583 remaining

#Assign mom ranks
CubData <- left_join(CubData, tblRanks,
                     by = c("mom" = "id", "year", "clan"), na_matches = "never")
colnames(CubData)[16] <- "stan_rank"
colnames(CubData)[17] <- "mom_stan_rank"

##Females assigned yearly ranks for each year where they are at least 18 months old at the start of the year
##Add mom's rank for anyone under age 2.5 without a rank yet
CubData$stan_rank <- ifelse(is.na(CubData$stan_rank) & !is.na(CubData$Actor_Age) & (CubData$Actor_Age < 2.5),
                            CubData$mom_stan_rank, CubData$stan_rank)
summary(CubData$stan_rank)     #11726 remaining

##Males assigned yearly ranks for each year where they are at least 5 years old at the start of the year [or when they are immigrants]
##Add mom's rank for males under age 5 without a rank yet
CubData$stan_rank <- ifelse(is.na(CubData$stan_rank) & !is.na(CubData$Actor_Age) & (CubData$Actor_Age < 6) & CubData$Actor_Sex == "Natal_Male",
                            CubData$mom_stan_rank, CubData$stan_rank)
summary(CubData$stan_rank)     #9944 remaining
CubData$last_year <- NA
for(year.step in 1:4){
  CubData$last_year <- CubData$year - year.step
  CubData <- left_join(CubData, tblRanks,
                       by = c("mom" = "id", "last_year" = "year", "clan"), na_matches = "never")
  colnames(CubData)[16] <- "stan_rank"
  colnames(CubData)[19] <- "mom_stan_rank_last_year"

  if(year.step < 2){
    ##Females assigned yearly ranks for each year where they are at least 18 months old at the start of the year
    ##Add mom's rank for anyone under age 2.5 without a rank yet
    CubData$stan_rank <- ifelse(is.na(CubData$stan_rank) & !is.na(CubData$Actor_Age) & (CubData$Actor_Age < 2.5),
                                CubData$mom_stan_rank_last_year, CubData$stan_rank)
  }
  if(year.step < 5){
    ##Males assigned yearly ranks for each year where they are at least 5 years old at the start of the year [or when they are immigrants]
    ##Add mom's rank for males under age 5 without a rank yet
    CubData$stan_rank <- ifelse(is.na(CubData$stan_rank) & !is.na(CubData$Actor_Age) & (CubData$Actor_Age < 6) & CubData$Actor_Sex == "Natal_Male",
                                CubData$mom_stan_rank_last_year, CubData$stan_rank)
  }
  CubData <- CubData[,c(1:19)]
  print(summary(CubData$stan_rank))
}

#Look at remaining animals
remaining <- filter(CubData, is.na(stan_rank))
remaining <- unique(remaining[,c("clan", "year", "id")])
remaining <- filter(remaining, !(clan == "talek" & year > 2013))
check <- remaining %>% group_by(clan, id) %>% summarise(num_years = length(year))
CubData <- CubData[!is.na(CubData$stan_rank),]

cub.unsol.appease.zipoiss.RANK <- glmmTMB(Count ~ Actor_Sex + Group_Size + stan_rank +
                                            offset(log(Duration)) +
                                            (1|session),
                                          data = CubData,
                                          ziformula = ~1,
                                          family = poisson())

summary(cub.unsol.appease.zipoiss.RANK)

#Subs

SubData$year <- as.numeric(SubData$year)
SubData$id <- SubData$Actor_ID


#Assign known ranks
SubData <- left_join(SubData, tblRanks,
                     by = c("id", "year", "clan"), na_matches = "never")
summary(SubData$stan_rank)

#Assign mom ranks
SubData <- left_join(SubData, tblRanks,
                     by = c("mom" = "id", "year", "clan"), na_matches = "never")
colnames(SubData)[16] <- "stan_rank"
colnames(SubData)[17] <- "mom_stan_rank"

##Females assigned yearly ranks for each year where they are at least 18 months old at the start of the year
##Add mom's rank for anyone under age 2.5 without a rank yet
SubData$stan_rank <- ifelse(is.na(SubData$stan_rank) & !is.na(SubData$Actor_Age) & (SubData$Actor_Age < 2.5),
                            SubData$mom_stan_rank, SubData$stan_rank)
summary(SubData$stan_rank)

##Males assigned yearly ranks for each year where they are at least 5 years old at the start of the year [or when they are immigrants]
##Add mom's rank for males under age 5 without a rank yet
SubData$stan_rank <- ifelse(is.na(SubData$stan_rank) & !is.na(SubData$Actor_Age) & (SubData$Actor_Age < 6) & SubData$Actor_Sex == "Natal_Male",
                            SubData$mom_stan_rank, SubData$stan_rank)
summary(SubData$stan_rank)
SubData$last_year <- NA
for(year.step in 1:4){
  SubData$last_year <- SubData$year - year.step
  SubData <- left_join(SubData, tblRanks,
                       by = c("mom" = "id", "last_year" = "year", "clan"), na_matches = "never")
  colnames(SubData)[16] <- "stan_rank"
  colnames(SubData)[19] <- "mom_stan_rank_last_year"

  if(year.step < 2){
    ##Females assigned yearly ranks for each year where they are at least 18 months old at the start of the year
    ##Add mom's rank for anyone under age 2.5 without a rank yet
    SubData$stan_rank <- ifelse(is.na(SubData$stan_rank) & !is.na(SubData$Actor_Age) & (SubData$Actor_Age < 2.5),
                                SubData$mom_stan_rank_last_year, SubData$stan_rank)
  }
  if(year.step < 5){
    ##Males assigned yearly ranks for each year where they are at least 5 years old at the start of the year [or when they are immigrants]
    ##Add mom's rank for males under age 5 without a rank yet
    SubData$stan_rank <- ifelse(is.na(SubData$stan_rank) & !is.na(SubData$Actor_Age) & (SubData$Actor_Age < 6) & SubData$Actor_Sex == "Natal_Male",
                                SubData$mom_stan_rank_last_year, SubData$stan_rank)
  }
  SubData <- SubData[,c(1:19)]
  print(summary(SubData$stan_rank))
}

#Look at remaining animals
remaining <- filter(SubData, is.na(stan_rank))
remaining <- unique(remaining[,c("clan", "year", "id")])
remaining <- filter(remaining, !(clan == "talek" & year > 2013))
check <- remaining %>% group_by(clan, id) %>% summarise(num_years = length(year))
SubData <- SubData[!is.na(SubData$stan_rank),]


sub.unsol.appease.zipoiss.RANK <- glmmTMB(Count ~ Actor_Sex + Group_Size + stan_rank +
                                            offset(log(Duration)) +
                                            (1|session),
                                          data = SubData,
                                          ziformula = ~1,
                                          family = poisson())

summary(sub.unsol.appease.zipoiss.RANK)

# just a bit of house keeping / data checking
class(Data_File_Sub1$Count)
Data_File_Sub1$Actor_Sex <- as.factor(Data_File_Sub1$Actor_Sex)
Data_File_Sub1$Actor_Age_Class <- as.factor(Data_File_Sub1$Actor_Age_Class)
levels(Data_File_Sub1$Actor_Age_Class)

# descriptive statistics table
bivar_agg <- Data_File_Sub1 %>%
  group_by(Actor_Age_Class,Actor_Sex) %>%
  summarise (n.agg.cnt = sum(!is.na(Count)),
             mean.agg.cnt = round (mean(Count,
                                        na.rm = T),2),
             stdev.agg.cnt = round (sd(Count,
                                       na.rm = T), 2),
             med.agg.cnt = round(median(Count,
                                        na.rm = T), 2),
             min.agg.cnt = round(min(Count,
                                     na.rm = T), 2),
             max.agg.cnt = round(max(Count,
                                     na.rm = T), 2)) %>%
  ungroup()


#Cub model

cub_sub_est <- tidy(cub.unsol.appease.zipoiss.RANK,  effects = "fixed", conf.int=TRUE,
                    exponentiate = T)

# Remove intercepts and covariates and unencessary columns
cub_sub_est <- cub_sub_est %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# add age class variable
cub_sub_est$age.class <- c('cub')


# Sub model

# use broom to create a tidy data frame of the sub aggression model output
sub_sub_est <- tidy(sub.unsol.appease.zipoiss.RANK,  effects = "fixed", conf.int=TRUE,
                    exponentiate = T)

# Remove intercepts and covariates and unencessary columns
sub_sub_est <- sub_sub_est %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# add age class variable
sub_sub_est$age.class <- c('subadult')

# Adult model

# use broom to create a tidy data frame of the adult aggression model output
adult_sub_est <- tidy(adult.unsol.appease.zipoiss.RANK,  effects = "fixed", conf.int=TRUE,
                      exponentiate = T)


# Remove intercepts and covariates and unencessary columns
adult_sub_est <- adult_sub_est %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# add age class variable
adult_sub_est$age.class <- c('adult', 'adult')


# combine tidy cub, sub and adult estmiates dataframes for graphing
sub_est <- rbind(cub_sub_est, sub_sub_est, adult_sub_est)
sub_est$age.class[sub_est$age.class == "cub"] <- "Male Cubs"
sub_est$age.class[sub_est$age.class == "subadult"] <- "Male Sub-Adults"
sub_est$age.class[sub_est$age.class == "adult" &
                    sub_est$term == "Actor_SexImmigrant_Male"] <- "Adult Immigant Males"
sub_est$age.class[sub_est$age.class == "adult" &
                    sub_est$term == "Actor_SexNatal_Male"] <- "Adult Natal Males"
# set the levels and order of the age.class factor
sub_est <-
  transform(sub_est,
            age.class = factor(age.class,
                               levels = c('Male Cubs','Male Sub-Adults','Adult Immigant Males', 'Adult Natal Males')))

#Graph results using dotwhisker, broom, dplyr, and ggplot2 packages
sub_est.plot <-
  ggplot(sub_est, aes(x = age.class, y = estimate)) +
  geom_hline(yintercept = 1, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_vline(xintercept = c(1.5, 2.5), linetype =1, size = 1)+
  geom_point(aes(color = age.class), size = 2, alpha = 1,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(color = age.class, ymin=conf.low, ymax=conf.high), width=.1,
                position=position_dodge(.5)) +
  scale_color_manual(breaks = c('Male Cubs','Male Sub-Adults','Adult Immigant Males', 'Adult Natal Males'),
                     values = c('darkgreen', 'darkblue', 'brown','brown')) +
  #coord_flip() + # flip x and y axes
  theme(text = element_text(size=14)) +
  theme(legend.position = 'none') +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = 'black',
                                 size = 1, linetype = 'solid')) +
  scale_x_discrete(labels = c('Male\\nCubs', 'Male\\nSub-Adults', 'Adult\\nImmigrant\\nMales', 'Adult\\nNatal\\nMales'), expand = expansion(add = c(1, 0)))+
  xlab('Actor') +
  ylab("Incidence Rate Ratio of Submissive Acts")+
  annotate('text', x = 0.1, y = 1.03, label = expression(italic('Female submitters')), color = 'red', hjust = 0, size = 10/.pt)


jpeg('Figure 2.pdf', width = 7.50, height = 5.00)
sub_est.plot
dev.off()


print(sub_est.plot)

# Summary statistics for the data tables

AdultData$Rate <- AdultData$Count/AdultData$Group_Size/AdultData$Duration

Adult_Summary <- summarize(group_by(AdultData,Actor_Sex),
                           lengthID = length(unique(Actor_ID)),
                           lengthSession = length(session),
                           sumDuration = sum(Duration),
                           sumCount = sum(Count),
                           meanRate = mean(Rate),
                           SDRate = sd(Rate))
Adult_Summary

SubData$Rate <- SubData$Count/SubData$Group_Size/SubData$Duration

Sub_Summary <- summarize(group_by(SubData,Actor_Sex),
                         lengthID = length(unique(Actor_ID)),
                         lengthSession = length(session),
                         sumDuration = sum(Duration),
                         sumCount = sum(Count),
                         meanRate = mean(Rate),
                         SDRate = sd(Rate))
Sub_Summary

CubData$Rate <- CubData$Count/CubData$Group_Size/CubData$Duration

Cub_Summary <- summarize(group_by(CubData,Actor_Sex),
                         lengthID = length(unique(Actor_ID)),
                         lengthSession = length(session),
                         sumDuration = sum(Duration),
                         sumCount = sum(Count),
                         meanRate = mean(Rate),
                         SDRate = sd(Rate))
Cub_Summary

Submission_Summary <- rbind(Adult_Summary,Sub_Summary,Cub_Summary)
Submission_Summary
write.csv(Submission_Summary,"C:/Users/skmcc/Documents/R/R_wd/Submission_Summary.csv" )


AdultDataIMales <- subset(AdultData, Actor_Sex == "Immigrant_Male")
AdultDataNMales <- subset(AdultData, Actor_Sex == "Natal_Male")
AdultDataIMalesYear <- summarize(group_by(AdultDataIMales,year),
                                maleCount = length(unique(Actor_ID)))
AdultDataNMalesYear <- summarize(group_by(AdultDataNMales,year),
                                 maleCount = length(unique(Actor_ID)))
AdultDataIMalesYear$Status <- "Immigrant"
AdultDataNMalesYear$Status <- "Natal"
AllMalesYear <- rbind(AdultDataIMalesYear,AdultDataNMalesYear)
AllMalesYear
summarize(group_by(AllMalesYear,Status),
          meanMales = mean(maleCount))

