###Analysis code for OAJournos_PublicData.csv dataset. 

library(tidyverse)
library(here)
library(skimr) # install.packages('skimr')
library(kableExtra) # install.packages('kableExtra')
library(janitor) #install.packages('janitor')
library(readxl) #install.packages('readxl')
library(mice) #install.packages('mice')
library(VIM) #install.packages('VIM')
library(glmnet) #install.packages('glmnet')
library(mlogit) #install.packages('mlogit')
library(readr) #install.packages('readr')
library(naniar) #install.packages('naniar')
library(caret)   #install.packages('caret')
library(foreign)      #install.packages('foreign')
library(nnet)      #install.packages('nnet')
library(reshape2)      #install.packages('reshape2')
library(car)   #install.packages('car')
#install.packages('scales')


#Read in public data
survey_results <- read_csv(here("data", "OAJournos_PublicData.csv"))

#Change values to shorter names and also collapse to help with analyzing - if do not want to collapse, see next section
survey_results[survey_results == "Yes, I am employed by a news organization"] <- "Employed"
survey_results[survey_results == "Yes, I work as a freelancer for a news organization"] <- "Freelance"
survey_results[survey_results == "Master's degree"] <- "Grad degree"
survey_results[survey_results == "PhD degree"] <- "Grad degree"
survey_results[survey_results == "I did not receive a master's or PhD degree"] <- "NA"
survey_results[survey_results == "Extremely uncomfortable"] <-  "Uncomfortable"
survey_results[survey_results == "Somewhat comfortable" | survey_results == "Extremely comfortable"] <- "Comfortable"
survey_results[survey_results == "Science-focused magazine or other outlet"] <- "ScienceOutlet"
survey_results[survey_results == "More than 75 percent" | survey_results == "More than 75 percent of my news stories cite a scientific study"] <- "51+"
survey_results[survey_results == "26 to 50 percent" | survey_results == "26 to 50 percent of my news stories cite a scientific study"] <- "-50"
survey_results[survey_results == "51 to 75 percent" | survey_results == "51 to 75 percent of my news stories cite a scientific study"] <- "51+"
survey_results[survey_results == "25 percent or less" | survey_results == "25 percent or fewer of my news stories cite a scientific study"] <- "-50"
survey_results[survey_results == "Very important - I will not reference a scholarly study if I cannot access the full-text version"] <- "Important"
survey_results[survey_results == "Pretty important - I will take as many steps as I can to find a copy of the scientific study before relying on another source"] <- "Important"
survey_results[survey_results == "Somewhat important - I try to find the scientific study, but if I can't get a copy of it or it's too expensive, I will use another source"] <- "Somewhat important"
survey_results[survey_results == "It hasn't changed" | survey_results == "My impressions and views have stayed the same"] <- "Same"
survey_results[survey_results == "It's gotten harder"] <- "Harder"
survey_results[survey_results == "Yes, before COVID-19 research"] <- "Yes pre Covid"
survey_results[survey_results == "Yes, but only in the past year during coverage of COVID-19 research"] <- "Yes post Covid"
survey_results[survey_results == "I had heard of the term but did not know what it meant" | survey_results == "I had heard the term but did not know what it meant"] <- "Heard of"
survey_results[survey_results == "I had never heard of the term"] <- "No"
survey_results[survey_results == "I will cite or refer to them without concern to whether they have been peer reviewed or published"] <- "Cite"
survey_results[survey_results == "I will cite or refer to them only if I know they have been peer reviewed"] <- "Cautious/don't cite"
survey_results[survey_results == "I do not pay attention to if a study came from an open database"] <- "Cite"
survey_results[survey_results == "I will cite or refer to them only if I know they have been peer reviewed and published"] <- "Cautious/don't cite"
survey_results[survey_results == "I will not cite or refer to them"] <- "Cautious/don't cite"
survey_results[survey_results == "I have no concern about using them"] <- "Little/no concern"
survey_results[survey_results == "I am hesitant but will use them if I am familiar with the journal" | survey_results == "I am hesitant  but will use them from journals I trust"] <- "Hesitant"
survey_results[survey_results == "I do not pay attention to the open access status of a journal" | survey_results == "I do not pay attention to the open access status of articles in paywalled journals"] <- "Little/no concern"
survey_results[survey_results == "I have some concerns but will use them unless there are red flags"] <- "Hesitant"
survey_results[survey_results == "My impressions and views have improved"] <- "Improved"
survey_results[survey_results == "Not at all important"] <- "Not important"
survey_results[survey_results == "Pretty important"] <- "Important"
survey_results[survey_results == "Somewhat important"] <- "Not important"
survey_results[survey_results == "Very important"] <- "Important"


####Change values to shorter names without collapsing values.
survey_results[survey_results == "Yes, I am employed by a news organization"] <- "Employed"
survey_results[survey_results == "Yes, I work as a freelancer for a news organization"] <- "Freelance"
survey_results[survey_results == "Master's degree"] <- "Masters"
survey_results[survey_results == "PhD degree"] <- "PhD"
survey_results[survey_results == "I did not receive a master's or PhD degree"] <- "NA"
survey_results[survey_results == "Science-focused magazine or other outlet"] <- "ScienceOutlet"
survey_results[survey_results == "More than 75 percent" | survey_results == "More than 75 percent of my news stories cite a scientific study"] <- "75+"
survey_results[survey_results == "26 to 50 percent" | survey_results == "26 to 50 percent of my news stories cite a scientific study"] <- "26-50"
survey_results[survey_results == "51 to 75 percent" | survey_results == "51 to 75 percent of my news stories cite a scientific study"] <- "51-75"
survey_results[survey_results == "25 percent or less" | survey_results == "25 percent or fewer of my news stories cite a scientific study"] <- "-25"
survey_results[survey_results == "Very important - I will not reference a scholarly study if I cannot access the full-text version"] <- "Very"
survey_results[survey_results == "Pretty important - I will take as many steps as I can to find a copy of the scientific study before relying on another source"] <- "Pretty"
survey_results[survey_results == "Somewhat important - I try to find the scientific study, but if I can't get a copy of it or it's too expensive, I will use another source"] <- "Somewhat"
survey_results[survey_results == "It hasn't changed" | survey_results == "My impressions and views have stayed the same"] <- "Same"
survey_results[survey_results == "It's gotten harder"] <- "Harder"
survey_results[survey_results == "Yes, before COVID-19 research"] <- "Yes pre Covid"
survey_results[survey_results == "Yes, but only in the past year during coverage of COVID-19 research"] <- "Yes post Covid"
survey_results[survey_results == "I had heard of the term but did not know what it meant" | survey_results == "I had heard the term but did not know what it meant"] <- "Heard of"
survey_results[survey_results == "I had never heard of the term"] <- "No"
survey_results[survey_results == "I will cite or refer to them without concern to whether they have been peer reviewed or published"] <- "Cite"
survey_results[survey_results == "I will cite or refer to them only if I know they have been peer reviewed"] <- "Cite if peer reviewed"
survey_results[survey_results == "I do not pay attention to if a study came from an open database"] <- "Don't pay attention"
survey_results[survey_results == "I will cite or refer to them only if I know they have been peer reviewed and published"] <- "Cite if peer reviewed, published"
survey_results[survey_results == "I will not cite or refer to them"] <- "Don't cite"
survey_results[survey_results == "I have no concern about using them"] <- "No concern"
survey_results[survey_results == "I am hesitant but will use them if I am familiar with the journal" | survey_results == "I am hesitant  but will use them from journals I trust"] <- "Hesitant"
survey_results[survey_results == "I do not pay attention to the open access status of a journal" | survey_results == "I do not pay attention to the open access status of articles in paywalled journals"] <- "Don't pay attention"
survey_results[survey_results == "I have some concerns but will use them unless there are red flags"] <- "Some concern"
survey_results[survey_results == "My impressions and views have improved"] <- "Improved"
survey_results[survey_results == "Not at all important"] <- "Not at all"
survey_results[survey_results == "Pretty important"] <- "Pretty"
survey_results[survey_results == "Somewhat important"] <- "Somewhat"
survey_results[survey_results == "Very important"] <- "Very"


#Setting factors and levels - based on collapsed values
survey_results$audSize <- factor(survey_results$audSize, order = TRUE, levels =c("Regional/local", "National"))
survey_results$education <- factor(survey_results$education, order = TRUE, levels =c("BS or lower", "Grad degree"))
survey_results$comfortArts <- factor(survey_results$comfortArts, order = TRUE, levels = c("Uncomfortable", "Comfortable"))
survey_results$needFullText <- factor(survey_results$needFullText)
survey_results$needFullText <- relevel(survey_results$needFullText, "Not important")
survey_results$needFreeText <- factor(survey_results$needFreeText, order = TRUE, levels =c("Not important", "Important"))
survey_results$prepKnow <- factor(survey_results$prepKnow, order = TRUE, levels = c("No", "Heard of", "Yes pre Covid", "Yes post Covid"))
survey_results$postKnow <- factor(survey_results$postKnow, order = TRUE, levels = c("No", "Heard of", "Yes pre Covid", "Yes post Covid"))
survey_results$OAKnow <- factor(survey_results$OAKnow, order = TRUE, levels = c("No", "Heard of", "Yes pre Covid", "Yes post Covid"))
survey_results$greenOAKnow <- factor(survey_results$greenOAKnow, order = TRUE, levels = c("Cautious/don't cite", "Cite"))
survey_results$goldOAKnow <- factor(survey_results$goldOAKnow, order = TRUE, levels = c("Hesitant", "Little/no concern"))
survey_results$hybridOAKnow <- factor(survey_results$hybridOAKnow, order = TRUE, levels = c("Hesitant", "Little/no concern"))

####Set factors and levels - based on distinct values
survey_results$audSize <- factor(survey_results$audSize, order = TRUE, levels =c("Regional/local", "National"))
survey_results$education <- factor(survey_results$education, order = TRUE, levels =c("BS or lower", "Masters", "PhD"))
survey_results$comfortArts <- factor(survey_results$comfortArts, order = TRUE, levels = c("Extremely uncomfortable", "Somewhat comfortable", "Extremely comfortable"))
survey_results$needFullText <- factor(survey_results$needFullText)
survey_results$needFullText <- relevel(survey_results$needFullText, "Somewhat")
survey_results$needFreeText <- factor(survey_results$needFreeText, order = TRUE, levels =c("Not at all", "Somewhat", "Pretty", "Very"))
survey_results$prepKnow <- factor(survey_results$prepKnow, order = TRUE, levels = c("No", "Heard of", "Yes pre Covid", "Yes post Covid"))
survey_results$postKnow <- factor(survey_results$postKnow, order = TRUE, levels = c("No", "Heard of", "Yes pre Covid", "Yes post Covid"))
survey_results$OAKnow <- factor(survey_results$OAKnow, order = TRUE, levels = c("No", "Heard of", "Yes pre Covid", "Yes post Covid"))
survey_results$greenOAKnow <- factor(survey_results$greenOAKnow, order = TRUE, levels = c("Don't pay attention", "Don't cite", "Cite if peer reviewed, published", "Cite if peer reviewed", "Cite"))
survey_results$goldOAKnow <- factor(survey_results$goldOAKnow, order = TRUE, levels = c("Don't pay attention", "Hesitant", "Some concern", "No concern"))
survey_results$hybridOAKnow <- factor(survey_results$hybridOAKnow, order = TRUE, levels = c("Don't pay attention", "Hesitant", "Some concern", "No concern"))

#Test that factor setting worked - should return TRUE. If FALSE, not factors
print(is.factor(survey_results$audSize))

#Test that the levels were set correctly - should be Some BS or lower, Masters, PhD
levels(survey_results$education)

#One way to get overview of my data
skim(survey_results)

#Look for missing values by column. Returns the numerical row for each missing value. 
colSums(is.na(survey_results))




summary(survey_results)

table(survey_results$USJourno)
table(survey_results$education)
table(survey_results$degreeSubj)
table(survey_results$outletType)
table(survey_results$audSize)
table(survey_results$timeSciJ)
table(survey_results$percentCite)
table(survey_results$comfortArts)
table(survey_results$peerReview)
table(survey_results$needFreeText)
table(survey_results$needFullText)
table(survey_results$ideaSources)
table(survey_results$prepKnow)
table(survey_results$postKnow)
table(survey_results$OAKnow)
table(survey_results$greenOAKnow)
table(survey_results$goldOAKnow)
table(survey_results$hybridOAKnow)
table(survey_results$knowPredatory)
table(survey_results$predConcern)


#Cross tab to compare IVs to DVs - plug and chug in the IVs and DVs as want to compare
survey_results %>%
  group_by(percentCite, comfortArts) %>% 
  summarize(count_by_needGold = n())

#Run chi square test for knowledge of preprints and postprints with all different IVs
chisq.test(table(survey_results$prepKnow, survey_results$education))
chisq.test(table(survey_results$postKnow, survey_results$education))
chisq.test(table(survey_results$prepKnow, survey_results$percentCite))
chisq.test(table(survey_results$postKnow, survey_results$percentCite))
chisq.test(table(survey_results$prepKnow, survey_results$audSize))
chisq.test(table(survey_results$postKnow, survey_results$audSize))
chisq.test(table(survey_results$prepKnow, survey_results$USJourno))
chisq.test(table(survey_results$postKnow, survey_results$USJourno))
chisq.test(table(survey_results$prepKnow, survey_results$timeSciJ))
chisq.test(table(survey_results$postKnow, survey_results$timeSciJ))
chisq.test(table(survey_results$prepKnow, survey_results$needFullText))
chisq.test(table(survey_results$postKnow, survey_results$needFullText))
chisq.test(table(survey_results$prepKnow, survey_results$needFreeText))
chisq.test(table(survey_results$postKnow, survey_results$needFreeText))
chisq.test(table(survey_results$prepKnow, survey_results$yearsWorkedGrouped))
chisq.test(table(survey_results$postKnow, survey_results$yearsWorkedGrouped))
#They all got errors - aka had at least one cell with too low frequency count. 
#so means need to use fisher's test on all of them

#Fisher's exact test for preprint and postprint knowledge vs. all IVs
fisher.test(table(survey_results$prepKnow, survey_results$education))
fisher.test(table(survey_results$postKnow, survey_results$education))
fisher.test(table(survey_results$prepKnow, survey_results$percentCite))
fisher.test(table(survey_results$postKnow, survey_results$percentCite))
fisher.test(table(survey_results$prepKnow, survey_results$audSize))
fisher.test(table(survey_results$postKnow, survey_results$audSize))
fisher.test(table(survey_results$prepKnow, survey_results$USJourno))
fisher.test(table(survey_results$postKnow, survey_results$USJourno))
fisher.test(table(survey_results$prepKnow, survey_results$timeSciJ))
fisher.test(table(survey_results$postKnow, survey_results$timeSciJ))
fisher.test(table(survey_results$prepKnow, survey_results$needFullText))
fisher.test(table(survey_results$postKnow, survey_results$needFullText))
fisher.test(table(survey_results$prepKnow, survey_results$needFreeText))
fisher.test(table(survey_results$postKnow, survey_results$needFreeText))

#IVs with p value < 0.05: 
  #For preprints, education, percentCite, audSize
  #For postprints, percentCite, audSize, 


chisq.test(table(survey_results$greenOAKnow, survey_results$education))
chisq.test(table(survey_results$goldOAKnow, survey_results$education))
chisq.test(table(survey_results$hybridOAKnow, survey_results$education))
chisq.test(table(survey_results$greenOAKnow, survey_results$percentCite))
chisq.test(table(survey_results$goldOAKnow, survey_results$percentCite))
chisq.test(table(survey_results$hybridOAKnow, survey_results$percentCite))
chisq.test(table(survey_results$greenOAKnow, survey_results$audSize))
chisq.test(table(survey_results$goldOAKnow, survey_results$audSize))
chisq.test(table(survey_results$hybridOAKnow, survey_results$audSize))
chisq.test(table(survey_results$greenOAKnow, survey_results$USJourno))
chisq.test(table(survey_results$goldOAKnow, survey_results$USJourno))
chisq.test(table(survey_results$hybridOAKnow, survey_results$USJourno))
chisq.test(table(survey_results$greenOAKnow, survey_results$timeSciJ))
chisq.test(table(survey_results$goldOAKnow, survey_results$timeSciJ))
chisq.test(table(survey_results$hybridOAKnow, survey_results$timeSciJ))
chisq.test(table(survey_results$greenOAKnow, survey_results$needFullText))
chisq.test(table(survey_results$goldOAKnow, survey_results$needFullText))
chisq.test(table(survey_results$hybridOAKnow, survey_results$needFullText))
chisq.test(table(survey_results$greenOAKnow, survey_results$needFreeText))
chisq.test(table(survey_results$goldOAKnow, survey_results$needFreeText))
chisq.test(table(survey_results$hybridOAKnow, survey_results$needFreeText))
chisq.test(table(survey_results$greenOAKnow, survey_results$outletType))
chisq.test(table(survey_results$goldOAKnow, survey_results$outletType))
chisq.test(table(survey_results$hybridOAKnow, survey_results$outletType))
#IVs with too small frequency counts:
  #For green - 
  #For gold - percentCite
  #For hybrid - education, percent
#IVs with significance:
  #For green - percentCite, 
  #For gold - education, 
  #For hybrid - 

fisher.test(table(survey_results$greenOAKnow, survey_results$education))
fisher.test(table(survey_results$goldOAKnow, survey_results$education))
fisher.test(table(survey_results$hybridOAKnow, survey_results$education))
fisher.test(table(survey_results$greenOAKnow, survey_results$percentCite))
fisher.test(table(survey_results$goldOAKnow, survey_results$percentCite))
fisher.test(table(survey_results$hybridOAKnow, survey_results$percentCite))
fisher.test(table(survey_results$greenOAKnow, survey_results$audSize))
fisher.test(table(survey_results$goldOAKnow, survey_results$audSize))
fisher.test(table(survey_results$hybridOAKnow, survey_results$audSize))
fisher.test(table(survey_results$greenOAKnow, survey_results$USJourno))
fisher.test(table(survey_results$goldOAKnow, survey_results$USJourno))
fisher.test(table(survey_results$hybridOAKnow, survey_results$USJourno))
fisher.test(table(survey_results$greenOAKnow, survey_results$timeSciJ))
fisher.test(table(survey_results$goldOAKnow, survey_results$timeSciJ))
fisher.test(table(survey_results$hybridOAKnow, survey_results$timeSciJ))
fisher.test(table(survey_results$greenOAKnow, survey_results$needFullText))
fisher.test(table(survey_results$goldOAKnow, survey_results$needFullText))
fisher.test(table(survey_results$hybridOAKnow, survey_results$needFullText))
fisher.test(table(survey_results$greenOAKnow, survey_results$needFreeText))
fisher.test(table(survey_results$goldOAKnow, survey_results$needFreeText))
fisher.test(table(survey_results$hybridOAKnow, survey_results$needFreeText))
fisher.test(table(survey_results$greenOAKnow, survey_results$outletType))
fisher.test(table(survey_results$goldOAKnow, survey_results$outletType))
fisher.test(table(survey_results$hybridOAKnow, survey_results$outletType))
fisher.test(table(survey_results$greenOAKnow, survey_results$prepKnow))
fisher.test(table(survey_results$goldOAKnow, survey_results$prepKnow))
fisher.test(table(survey_results$hybridOAKnow, survey_results$prepKnow))
fisher.test(table(survey_results$greenOAKnow, survey_results$postKnow))
fisher.test(table(survey_results$goldOAKnow, survey_results$postKnow))
fisher.test(table(survey_results$hybridOAKnow, survey_results$postKnow))

#Determine how many people actually answered the question of how they check if journal is predatory
predKnow <- subset(survey_results, knowPredatory == "Yes")
predConcern <- subset(predKnow, predConcern != "No")
table(survey_results$knowPredatory)
table(predKnow$predConcern)
table(predConcern$predConcern)
which(is.na(predConcern$predCheck))
sum(predConcern$predCheck == "")

#Same as above but for checking the websites
which(is.na(predConcern$websiteCheck))
sum(predConcern$websiteCheck == "")


  table(survey_results$knowPredatory)


###Data visuals - create a bunch of bar charts to visually explore data

ggplot(data = survey_results, aes(x = education)) +
  geom_bar()

education_graph <- survey_results %>%
  filter(!is.na(education)) %>%
  ggplot(aes(x= education)) + 
  geom_bar()

education_graph

yearsWorked_graph <- survey_results  %>%
  filter(!is.na(yearsWorked)) %>%
  ggplot(aes(x = yearsWorked)) +
  geom_bar()

yearsWorked_graph

ggplot(data = survey_results, aes(x = outletType)) +
  geom_bar()

ggplot(data = survey_results, aes(x = audSize)) +
  geom_bar()

ggplot(data = survey_results, aes(x = timeSciJ)) +
  geom_bar()

ggplot(data = survey_results, aes(x = comfortArts)) +
  geom_bar()

ggplot(data = survey_results, aes(x = peerReview)) +
  geom_bar()

ggplot(data = survey_results, aes(x = percentCite)) +
  geom_bar()  

ggplot(data = survey_results, aes(x = needFullText)) +
  geom_bar()

ggplot(data = survey_results, aes(x = needFreeText)) +
  geom_bar()

ggplot(data = survey_results, aes(x = easyFreeText)) +
  geom_bar()

easyFree_graph <- survey_results %>%
  filter(!is.na(easyFreeText)) %>%
  ggplot(aes(x= easyFreeText)) + 
  geom_bar()

easyFree_graph

ggplot(data = survey_results, aes(x = prepKnow)) +
  geom_bar()

ggplot(data = survey_results, aes(x = postKnow)) +
  geom_bar()

ggplot(data = survey_results, aes(x = OAKnow)) +
  geom_bar()

ggplot(data = survey_results, aes(x = greenOAKnow)) +
  geom_bar() +
  coord_flip()

goldOAKnow_graph <- survey_results %>%
  filter(!is.na(goldOAKnow)) %>%
  ggplot(aes(x= goldOAKnow)) + 
  geom_bar()

goldOAKnow_graph

hybridOAKnow_graph <- survey_results %>%
  filter(!is.na(hybridOAKnow)) %>%
  ggplot(aes(x= hybridOAKnow)) + 
  geom_bar()

hybridOAKnow_graph

ggplot(data = survey_results, aes(x = hybridOAKnow)) +
  geom_bar() +
  coord_flip()

ggplot(data = survey_results, aes(x = viewsChanged)) +
  geom_bar()

ggplot(data = survey_results, aes(x = knowPredatory)) +
  geom_bar()

predConcern_graph <- survey_results %>%
  filter(!is.na(predConcern)) %>%
  ggplot(aes(x= predConcern)) + 
  geom_bar()

predConcern_graph

ggplot(data = survey_results, aes(x = predConcern)) +
  geom_bar()

greenKnowByEd_graph <- survey_results %>%
  filter(!is.na(greenOAKnow)) %>%
  filter(!is.na(education)) %>%
  ggplot(aes(x= greenOAKnow, fill = education)) + 
  geom_bar()
coord_flip()

greenKnowByEd_graph

ggplot(survey_results, aes(x = greenOAKnow, fill = education)) + 
  geom_bar() +
  coord_flip()


table(survey_results$needFreeText)