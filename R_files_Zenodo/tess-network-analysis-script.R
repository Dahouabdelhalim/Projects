setwd("tess-motivation-survey")

#----------- Preparation 

# Read data extracted from Coney
raw.data <- read.csv("tess-network-results.csv", sep =";")
raw.data <- cbind(raw.data, tag.question = paste0(raw.data$tag,".",raw.data$question_id))
# Remove test users
raw.data.filtered<-raw.data[!raw.data$user %in% c("testCefriel", "testCefriel2"),]

# Select unfinished survey completions
raw.data.unfinished<-raw.data.filtered[raw.data.filtered$totalDuration == 'unfinished', ]
# Remove unfinished survey completions
raw.data.finished<-raw.data.filtered[!raw.data.filtered$totalDuration == 'unfinished', ]
rm(raw.data.filtered)


#----------- Analyse finished survey completions

# Select answers to closed questions with values between 1 and 5

# Find users completing the survey multiple times, keep only the most recent survey completion
library(dplyr)
res <- raw.data.finished %>% group_by(user,session) %>% summarise(Freq=n())
res$duplicated.user <- duplicated(res$user)
duplicated.user<- res[res$duplicated.user == TRUE, "user"]

for(i in duplicated.user$user){
  b<-raw.data.finished[raw.data.finished$user == i, ]
  b.first.question<- b[b$question_id == unique(b$question_id[b$question_id == min(b$question_id)]),]
  b.first.question$date.time <- as.POSIXlt(paste(b.first.question$date, b.first.question$time), format="%d-%m-%Y %H:%M:%S")
  # Extract the session
  session.to.keep<-as.character(b.first.question$session[b.first.question$date.time == min(b.first.question$date.time)])
  session.to.delete<-as.character(b.first.question$session[!b.first.question$session == session.to.keep])
  # Remove the duplicated survey completion
  raw.data.finished<-raw.data.finished[! (raw.data.finished$user == i & raw.data.finished$session %in% session.to.delete), ]
}

write.csv(raw.data.finished, file =  "completed-survey-all-answers.csv", row.names = F)

# Select answers to closed questions with values between 1 and 5
raw.data.finished <- raw.data.finished[raw.data.finished$minValue == '1' & raw.data.finished$maxValue == '5' ,  ]

rm(i, session.to.delete, session.to.keep, res, b, b.first.question, duplicated.user, raw.data)

# Pivoting table: one row for each user and one column for each question, as value the numerical value given as answer
library(reshape)
pivot.for.question<-cast(raw.data.finished, user ~ tag.question, fun.aggregate = mean)
names(pivot.for.question) <- c("user" , "achievement.1" , "achievement.2" , "belongingness.1" , "belongingness.2",     
                               "benevolence.1" ,  "benevolence.2" , "conformity.1" , "conformity.2" ,  "data usage.1" , 
                               "data usage.2" ,  "engagement.1" , "engagement.2" ,   "engagement.3" ,   "global motivation",
                               "hedonism.1" ,    "hedonism.2" ,   "power.1" , "power.2" , "routine.1" ,    
                               "routine.2" ,     "self-direction.1" , "self-direction.2" , "stimulation.1" , "stimulation.2" , 
                               "universalism.1" ,  "universalism.2")

# Save result as a csv
write.csv(pivot.for.question, file =  "completed-survey-closed-answers.csv", row.names = F)


# Correlation analysis

completed<-read.csv(file =  "completed-survey-closed-answers.csv", header = T)

# Average of questions for each latent variable
completed$ach <- rowMeans(subset(completed, select = c(achievement.1, achievement.2)), na.rm = TRUE)
completed$bel <- rowMeans(subset(completed, select = c(belongingness.1, belongingness.2)), na.rm = TRUE)
completed$ben <- rowMeans(subset(completed, select = c(benevolence.1, benevolence.2)), na.rm = TRUE)
completed$conf <- rowMeans(subset(completed, select = c(conformity.1, conformity.2)), na.rm = TRUE)
completed$hed <- rowMeans(subset(completed, select = c(hedonism.1, hedonism.2)), na.rm = TRUE)
completed$pwr <- rowMeans(subset(completed, select = c(power.1, power.2)), na.rm = TRUE)
completed$rout <- rowMeans(subset(completed, select = c(routine.1, routine.2)), na.rm = TRUE)
completed$self <- rowMeans(subset(completed, select = c(self.direction.1, self.direction.2)), na.rm = TRUE)
completed$stim <- rowMeans(subset(completed, select = c(stimulation.1, stimulation.2)), na.rm = TRUE)
completed$univ <- rowMeans(subset(completed, select = c(universalism.1, universalism.2)), na.rm = TRUE)
completed$engag <- rowMeans(subset(completed, select = c(engagement.1, engagement.2, engagement.3)), na.rm = TRUE)
completed$datause <- rowMeans(subset(completed, select = c(data.usage.1, data.usage.2)), na.rm = TRUE)

# Correlation between global motivation and each latent variable
d.cor<-completed
global.corr <- data.frame(var=c(),corr=c(),pv=c(),sign=c(), global.motivation= c())
for(i in 28:39){
  t <- cor.test(d.cor[,i], d.cor$global.motivation)
  nome <- names(d.cor)[i]
  corr <- t$estimate
  pv <- t$p.value
  sign <- ifelse(pv<0.001, "***", ifelse(pv<0.01, "**", ifelse(pv<0.05, "*", "")))
  global.motivation<- paste0(round(corr,3),sign)
  global.corr <- rbind(global.corr, data.frame(nome, corr, pv, sign, global.motivation))
}


# Creation of the final csv summarising the obtained results
library(dplyr)
all.results<-data.frame(variable = character(), mean =  double(), stdev = double())

ach<-data.frame(variable = "achievement", mean =  mean(completed$ach), stdev = sd(completed$ach))
bel<-data.frame(variable = "belongingness",  mean = mean(completed$bel), stdev = sd(completed$bel))
ben<-data.frame(variable = "benevolence",  mean = mean(completed$ben), stdev = sd(completed$ben))
conf<-data.frame(variable = "conformity",  mean = mean(completed$conf), stdev = sd(completed$conf))
hed<-data.frame(variable = "hedonism",  mean = mean(completed$hed), stdev = sd(completed$hed))
pwr<-data.frame(variable = "power",  mean = mean(completed$pwr), stdev = sd(completed$pwr))
rout<-data.frame(variable = "routine",  mean = mean(completed$rout), stdev = sd(completed$rout))
self<-data.frame(variable = "self-direction",  mean = mean(completed$self), stdev = sd(completed$self))
stim<-data.frame(variable = "stimulation",  mean = mean(completed$stim), stdev = sd(completed$stim))
univ<-data.frame(variable = "universalism",  mean = mean(completed$univ), stdev = sd(completed$univ))

all.results<- rbind(all.results, ach, bel, ben, conf, hed, pwr, rout, self, stim, univ)
correlations <- global.corr

library(plyr)

correlations$nome<-mapvalues(correlations$nome, from = c("ach", "bel", "ben", "conf", "hed", "pwr", "rout", "self", "stim", "univ"), to = c("achievement", "belongingness", "benevolence", "conformity", "hedonism", "power", "routine", "self-direction", "stimulation", "universalism"))
correlations<-correlations[,c(1,2,3,4)]
names(correlations)<-c("variable", "correlation", "pvalue", "significance")

all.results<-merge(all.results, correlations, by="variable")
write.csv(all.results, "tess-network-analysis-results.csv", row.names= F)