setwd("open-soil-atlas")

## DATA PREPARATION
raw.data<- read.csv("open-soil-atlas-results.csv", header=T)

# Keep only completed surveys
raw.data.unfinished<-raw.data[raw.data$totalDuration == 'unfinished', ]
# length(unique(raw.data.unfinished$user))

raw.data<-raw.data[!raw.data$totalDuration == 'unfinished', ]
# length(unique(raw.data$user))

# Keep only questions related to motivation 
motivation.questions<-raw.data[raw.data$tag %in% c("achievement", "conformity", "self-direction", "stimulation", "routine", "hedonism", "power", "belongingness", "benevolence", "universalism", "global motivation"),]
unique(motivation.questions$tag)

# Create unique id made usign tag+id.question to perform pivoting of the table  
motivation.questions<-cbind(motivation.questions, tag.question = paste0(motivation.questions$tag, motivation.questions$questionId))


## --  pivoting table (one row for each user and one column for each question. As value the numerical value given as answer)
library(reshape2) 
library(reshape) 

match.tag.question<-unique(motivation.questions[, c("tag.question", "question")])
# Rows now represent users, each column is related to a tag.question. Value of each cell is taken from the "value" column of the original file (default behaviour of the function, if a "value" column doesn't exist in the original table you need to specify which column should be used).
pivot.motivation.questions<-cast(motivation.questions, user ~ tag.question , fun.aggregate = mean)
pivot.motivation.questions<-pivot.motivation.questions[ , !(names(pivot.motivation.questions) %in% c("global motivation1188"))]


# Modify question names to uniform analysis according to the survey-motivation-template
names(pivot.motivation.questions) <- c("user", "achievement.2", "achievement.1", "belongingness.1", "benevolence.1", "benevolence.2", "conformity.1", 
"global motivation", "hedonism.1", "power.2", "routine.1", "self-direction.1", "self-direction.2", "stimulation.1", 
"stimulation.2", "universalism.1", "universalism.2")

write.csv(pivot.motivation.questions, "open-soil-atlas-pivot-questions.csv", row.names = F)

#--------------------------- ANALYSIS RESULTS BY QUESTION

# Compute average and variance for each question 
ans<- read.csv("open-soil-atlas-pivot-questions.csv", header=TRUE)

mean.question<-round(sapply(ans[,c(2:17)], mean),2)
var.questions<-round(sapply(ans[,c(2:17)], var), 2)

df.questions<- data.frame(mean = mean.question, var = var.questions)


#--------------------------- ANALYSIS RESULTS BY MOTIVATING FACTOR

# Compute average and variance for each motivating factor 
ans<- read.csv("open-soil-atlas-pivot-questions.csv", header=TRUE)

library(dplyr)

factors.summary= data.frame(factor= character(), mean=numeric(), var=numeric())

all.factors<-c("Achievement", "Belongingness", "Benevolence", "Conformity", "Hedonism", "Power", "Routine", "Self.direction", "Stimulation", "Universalism", "Global.motivation")

for (k in all.factors){
  f<-ans %>% select(starts_with(k)) 
  num.col<-ncol(f)
  v.final=vector()
  for(i in 1:num.col){
    v<-as.vector(f[,i])
    v.final<-c(v.final, v)
  }
  
  f.mean<-round(mean(v.final), 2)
  f.var<-round(var(v.final),2)
  new.row<-data.frame(factor= k, mean=f.mean, var=f.var)
  factors.summary<-rbind(factors.summary, new.row) # Add row to the dataframe
  
  rm(f, f.mean, f.var, new.row, v.final, v, num.col)
  
}


#--------------------------- CORRELATION ANALYSIS

completed<- read.csv("open-soil-atlas-pivot-questions.csv", header=TRUE)

library(dplyr)
# average questions for each tag
# average results for each tag

#achievement
ach.subset<-completed %>% select(starts_with("ach"))
completed$ach <- rowMeans(ach.subset, na.rm = TRUE)

#belongingness
bel.subset<-completed %>% select(starts_with("bel"))
completed$bel <- rowMeans(bel.subset, na.rm = TRUE)

#benevolence
ben.subset<-completed %>% select(starts_with("ben"))
completed$ben <- rowMeans(ben.subset, na.rm = TRUE)

#conformity
conf.subset<-completed %>% select(starts_with("conf"))
completed$conf <- rowMeans(conf.subset, na.rm = TRUE)

#hedonism
hed.subset<-completed %>% select(starts_with("hed"))
completed$hed <- rowMeans(hed.subset, na.rm = TRUE)

#power
pwr.subset<-completed %>% select(starts_with("pow"))
completed$pwr <- rowMeans(pwr.subset, na.rm = TRUE)

#routine
rout.subset<-completed %>% select(starts_with("rout"))
completed$rout <- rowMeans(rout.subset, na.rm = TRUE)

#self-direction
self.subset<-completed %>% select(starts_with("self"))
completed$self <- rowMeans(self.subset, na.rm = TRUE)

#stimulation
stim.subset<-completed %>% select(starts_with("stim"))
completed$stim <- rowMeans(stim.subset, na.rm = TRUE)

#universalism
univ.subset<-completed %>% select(starts_with("univ"))
completed$univ <- rowMeans(univ.subset, na.rm = TRUE)


# Correlation between motivating factors and global motivation
d.cor<-completed
global.corr <- data.frame(var=c(),corr=c(),pv=c(),sign=c(), global.motivation= c())
for(i in 18:27){ # Change considering dataframe size
  t <- cor.test(d.cor[,i], d.cor$global.motivation)
  nome <- names(d.cor)[i]
  corr <- t$estimate
  pv <- t$p.value
  sign <- ifelse(pv<0.001, "***", ifelse(pv<0.01, "**", ifelse(pv<0.05, "*", "")))
  global.motivation<- paste0(round(corr,3),sign)
  global.corr <- rbind(global.corr, data.frame(nome, corr, pv, sign, global.motivation))
}
