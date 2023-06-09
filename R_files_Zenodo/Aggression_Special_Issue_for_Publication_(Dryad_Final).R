##############################################################
#                                                            #
#     Title: Developmental Sex Differences in Hyenas         #
#     Author: S. Kevin McCormick                             #
#     Date: 12/05/2022                                       #
#     Purpose: Comparisons of Unsolicited Aggressive Acts    #
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


#Import data
Data_File_Agg1 <- read.csv("AggressiveCountsSimplifiedClanFixMom_Dryad.csv")
tblRanks <- read.csv("tblRanks.csv")
Data_File_Agg1$year <- format(as.Date(Data_File_Agg1$date, format="%m/%d/%Y"),"%Y")


#Restrict it to the clans that have good data
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$clan == "talek" |
                                   Data_File_Agg1$clan == "serena.n"|
                                   Data_File_Agg1$clan == "serena.s",]

#Since we are looking at rates, we should restrict the duration
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$Duration > 10,]
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$Duration < 120,]

#For now we will also remove natal den location contexts
#Counter Attacks Removed in ACCESS prior to exporting the file
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$location != "n",]
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$location != "g",]

#Need to reliable sexes in this file since it is difficult do it in ACCESS
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$Actor_Sex != "u",]
Data_File_Agg1$Actor_Sex[Data_File_Agg1$Actor_Sex == "f"] <- "Female"
Data_File_Agg1$Actor_Sex[Data_File_Agg1$Actor_Sex == "m"] <- "Male"

#Something odd in the translation from ACCESS to R in the unknown age class
unique(Data_File_Agg1[Data_File_Agg1$Actor_Age_Class == "Unknown" &
                        Data_File_Agg1$Actor_Sex == "Female",]$Actor_ID)

Data_File_Agg1$Actor_Age_Class <- Data_File_Agg1$Actor_Age_Class2
Data_File_Agg1$Actor_Age_Class[Data_File_Agg1$Actor_Age <= 1] <- "Cub"
Data_File_Agg1$Actor_Age_Class[Data_File_Agg1$Actor_Age > 1] <- "Sub"
Data_File_Agg1$Actor_Age_Class[Data_File_Agg1$Actor_Age > 2] <- "Adult"
Data_File_Agg1$Actor_Age_Class[is.na(Data_File_Agg1$Actor_Age)] <- "Unknown"

unique(Data_File_Agg1[Data_File_Agg1$Actor_Age_Class == "Unknown" &
                        Data_File_Agg1$Actor_Sex == "Female",]$Actor_ID)

#Need to remove aliens
unique(Data_File_Agg1[Data_File_Agg1$Actor_Status == "u" &
                        Data_File_Agg1$Actor_Sex == "Female",]$Actor_ID)
unique(Data_File_Agg1[Data_File_Agg1$Actor_Status == "t" &
                        Data_File_Agg1$Actor_Sex == "Female",]$Actor_ID)
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$Actor_Status != "u",]
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$Actor_Status != "t",]

unique(Data_File_Agg1[Data_File_Agg1$Actor_Age_Class == "Unknown" &
                        Data_File_Agg1$Actor_Sex == "Female",]$Actor_ID)

#Make females from the beginning of the study that don't have ages adults
Data_File_Agg1$Actor_Age_Class[Data_File_Agg1$Actor_Age_Class == "Unknown" &
                                 Data_File_Agg1$Actor_Sex == "Female"] <- "Adult"


#Make a Male Adult Age class

Data_File_Agg1$Actor_Sex[Data_File_Agg1$Actor_Age > 2 &
                           Data_File_Agg1$Actor_Sex == "Male" &
                           Data_File_Agg1$Actor_Status == "r"] <- "Natal_Male"

###Need to check dates, and the most up to date would be tblclan membership
###Link on date range in tblClan membership

Data_File_Agg1$Actor_Age_Class[Data_File_Agg1$Actor_Sex == "Natal_Male"] <- "Adult"

Data_File_Agg1$Actor_Sex[Data_File_Agg1$Actor_Status == "i" &
                           Data_File_Agg1$Actor_Sex == "Male"] <- "Immigrant_Male"

Data_File_Agg1$Actor_Age_Class[Data_File_Agg1$Actor_Sex == "Immigrant_Male"] <- "Adult"

#Still have some unknowns in our actors
#Many of them are missing sex too
#These are relatively few, so I am just removing the entries or blanks

Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$Actor_Age_Class != "Unknown",]
Data_File_Agg1 <- Data_File_Agg1[Data_File_Agg1$Actor_Sex != "",]
str(Data_File_Agg1)

#Subset the data on Actor Age Class
#Reminder: We are interested in general spontaneous aggressive behavior of the Actor
#based on the sex of the acting aggressor.

AdultData <- Data_File_Agg1[Data_File_Agg1$Actor_Age_Class == "Adult",]
SubData <- Data_File_Agg1[Data_File_Agg1$Actor_Age_Class == "Sub",]
CubData <- Data_File_Agg1[Data_File_Agg1$Actor_Age_Class == "Cub",]

CubData$Duration <- CubData$Duration/60
SubData$Duration <- SubData$Duration/60
AdultData$Duration <- AdultData$Duration/60

#####Adding Rank to Adults
AdultData$year <- as.numeric(AdultData$year)
tblRanks$year <- as.numeric(tblRanks$year)
AdultData$id <- AdultData$Actor_ID
str(AdultData)

#Assign known ranks
AdultData <- left_join(AdultData, tblRanks,
                       by = c("id", "year", "clan"), na_matches = "never")
summary(AdultData$stan_rank)
str(AdultData)

#Assign mom ranks
AdultData <- left_join(AdultData, tblRanks,
                       by = c("mom" = "id", "year", "clan"), na_matches = "never")
colnames(AdultData)[16] <- "stan_rank"
colnames(AdultData)[17] <- "mom_stan_rank"
str(AdultData)

##Females assigned yearly ranks for each year where they are at least 18 months old at the start of the year
##Add mom's rank for anyone under age 2.5 without a rank yet
AdultData$stan_rank <- ifelse(is.na(AdultData$stan_rank) & !is.na(AdultData$Actor_Age) & (AdultData$Actor_Age < 2.5),
                              AdultData$mom_stan_rank, AdultData$stan_rank)
summary(AdultData$stan_rank)     

##Males assigned yearly ranks for each year where they are at least 5 years old at the start of the year [or when they are immigrants]
##Add mom's rank for males under age 5 without a rank yet
AdultData$stan_rank <- ifelse(is.na(AdultData$stan_rank) & !is.na(AdultData$Actor_Age) & (AdultData$Actor_Age < 6) & AdultData$Actor_Sex == "Natal_Male",
                              AdultData$mom_stan_rank, AdultData$stan_rank)
summary(AdultData$stan_rank)     
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
  AdultData <- AdultData[,c(1:18)]
  print(summary(AdultData$stan_rank))
}

#Look at remaining animals
remaining <- filter(AdultData, is.na(stan_rank))
remaining <- unique(remaining[,c("clan", "year", "id")])
remaining <- filter(remaining, !(clan == "talek" & year > 2013))
check <- remaining %>% group_by(clan, id) %>% summarise(num_years = length(year))
AdultData <- AdultData[!is.na(AdultData$stan_rank),] 

adult.aggression.zipoiss.RANK <- glmmTMB(Count ~ Actor_Sex + Group_Size + stan_rank +
                                           offset(log(Duration)),
                                         data = AdultData,
                                         ziformula = ~1,
                                         family = list(family = 'poisson',
                                                       link = 'log'))

summary(adult.aggression.zipoiss.RANK)


#Cubs adding rank
CubData$year <- as.numeric(CubData$year)
CubData$id <- CubData$Actor_ID

#Assign known ranks
CubData <- left_join(CubData, tblRanks,
                     by = c("id", "year", "clan"), na_matches = "never")

summary(CubData$stan_rank)     

#Assign mom ranks
CubData <- left_join(CubData, tblRanks,
                     by = c("mom" = "id", "year", "clan"), na_matches = "never")
colnames(CubData)[16] <- "stan_rank"
colnames(CubData)[17] <- "mom_stan_rank"

##Females assigned yearly ranks for each year where they are at least 18 months old at the start of the year
##Add mom's rank for anyone under age 2.5 without a rank yet
CubData$stan_rank <- ifelse(is.na(CubData$stan_rank) & !is.na(CubData$Actor_Age) & (CubData$Actor_Age < 2.5),
                            CubData$mom_stan_rank, CubData$stan_rank)
summary(CubData$stan_rank)     

##Males assigned yearly ranks for each year where they are at least 5 years old at the start of the year [or when they are immigrants]
##Add mom's rank for males under age 5 without a rank yet
CubData$stan_rank <- ifelse(is.na(CubData$stan_rank) & !is.na(CubData$Actor_Age) & (CubData$Actor_Age < 6) & CubData$Actor_Sex == "Natal_Male",
                            CubData$mom_stan_rank, CubData$stan_rank)
summary(CubData$stan_rank)     
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
  CubData <- CubData[,c(1:18)]
  print(summary(CubData$stan_rank))
}

#Look at remaining animals
remaining <- filter(CubData, is.na(stan_rank))
remaining <- unique(remaining[,c("clan", "year", "id")])
remaining <- filter(remaining, !(clan == "talek" & year > 2013))
check <- remaining %>% group_by(clan, id) %>% summarise(num_years = length(year))
CubData <- CubData[!is.na(CubData$stan_rank),] 

cub.aggression.zipoiss.RANK <- glmmTMB(Count ~ Actor_Sex + Group_Size + stan_rank +
                                         offset(log(Duration)) +
                                         (1|session),
                                       data = CubData,
                                       ziformula = ~1,
                                       family = list(family = 'poisson',
                                                     link = 'log'))

summary(cub.aggression.zipoiss.RANK)


#Subs adding rank
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
  SubData <- SubData[,c(1:18)]
  print(summary(SubData$stan_rank))
}

#Look at remaining animals
remaining <- filter(SubData, is.na(stan_rank))
remaining <- unique(remaining[,c("clan", "year", "id")])
remaining <- filter(remaining, !(clan == "talek" & year > 2013))
check <- remaining %>% group_by(clan, id) %>% summarise(num_years = length(year))
SubData <- SubData[!is.na(SubData$stan_rank),] 

sub.aggression.zipoiss.RANK <- glmmTMB(Count ~ Actor_Sex + Group_Size + stan_rank +
                                         offset(log(Duration)),
                                       data = SubData,
                                       ziformula = ~1,
                                       family = list(family = 'poisson',
                                                     link = 'log'))

summary(sub.aggression.zipoiss.RANK)

# just a bit of house keeping / data checking before prepping the figure
class(Data_File_Agg1$Count)
Data_File_Agg1$Actor_Sex <- as.factor(Data_File_Agg1$Actor_Sex)
Data_File_Agg1$Actor_Age_Class <- as.factor(Data_File_Agg1$Actor_Age_Class)
levels(Data_File_Agg1$Actor_Age_Class)

#cub model

cub_agg_est <- tidy(cub.aggression.zipoiss.RANK,  effects = "fixed", conf.int=TRUE,
                    exponentiate = T)

# Remove intercepts and covariates and unencessary columns
cub_agg_est <- cub_agg_est %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# add age class variable
cub_agg_est$age.class <- c('cub')

# sub model

# use broom to create a tidy data frame of the sub aggression model output
sub_agg_est <- tidy(sub.aggression.zipoiss.RANK,  effects = "fixed", conf.int=TRUE,
                    exponentiate = T)

# Remove intercepts and covariates and unencessary columns
sub_agg_est <- sub_agg_est %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# add age class variable
sub_agg_est$age.class <- c('subadult')

# Kevin's adult model

# use broom to create a tidy data frame of the adult aggression model output
adult_agg_est <- tidy(adult.aggression.zipoiss.RANK,  effects = "fixed", conf.int=TRUE,
                      exponentiate = T)



# Remove intercepts and covariates and unencessary columns
adult_agg_est <- adult_agg_est %>%
  filter(grepl('Actor', term)) %>%
  select(term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# add age class variable
adult_agg_est$age.class <- c('adult', 'adult')


# combine tidy cub, sub and adult estmiates dataframes for graphing
agg_est <- rbind(cub_agg_est, sub_agg_est, adult_agg_est)
agg_est$age.class[agg_est$age.class == "cub"] <- "Male Cubs"
agg_est$age.class[agg_est$age.class == "subadult"] <- "Male Sub-Adults"
agg_est$age.class[agg_est$age.class == "adult" &
                    agg_est$term == "Actor_SexImmigrant_Male"] <- "Adult Immigant Males"
agg_est$age.class[agg_est$age.class == "adult" &
                    agg_est$term == "Actor_SexNatal_Male"] <- "Adult Natal Males"

# set the levels and order of the age.class factor
agg_est <-
  transform(agg_est,
            age.class = factor(age.class,
                               levels = c('Male Cubs','Male Sub-Adults', 'Adult Natal Males', 'Adult Immigant Males')))
#Graph results using dotwhisker, broom, dplyr, and ggplot2 packages
agg_est.plot <-
  ggplot(agg_est, aes(x = age.class, y = estimate)) +
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
  ylab("Incidence Rate Ratio of Aggressive Acts")+
  annotate('text', x = 0.1, y = 1.02, label = expression(italic('Female aggressors')), color = 'red', hjust = 0, size = 10/.pt)


pdf('Figure 1.pdf', width = 7.50, height = 5.00)
agg_est.plot
dev.off()


print(agg_est.plot)

## f) Save Plot
# use ggsave to save the linearization plot
ggsave('agg_est.plot.pdf', plot = agg_est.plot,
       device = NULL,
       path = here(),
       scale = 1, width = 12,
       height = 9,
       units = c('in'), dpi = 300, limitsize = TRUE)

save(agg_est, adult.aggression.zipoiss.RANK, cub.aggression.zipoiss.RANK, sub.aggression.zipoiss.RANK,  file = 'data_Figure_1.Rdata')

#Summary Stats

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

Aggression_Summary <- rbind(Adult_Summary,Sub_Summary,Cub_Summary)
Aggression_Summary
write.csv(Aggression_Summary,"C:/Users/skmcc/Documents/R/R_wd/Aggression_Summary.csv" )


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

