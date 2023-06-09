setwd("water-sentinels")


##---------- PREPARAZIONE DATI
raw.data<- read.csv("water-sentinels-results.csv", header=T)
raw.data<-raw.data[1:4,8:28]

raw.data.backup<-raw.data

library(stringr)
# Keep only tag as header for the column
names(raw.data)<-str_replace(names(raw.data), "\\\\.\\\\.\\\\..+","")

# Substitute answers with their numeric value
library(plyr)
#Stimulation.1
levels(raw.data$Stimulation.1)
raw.data$Stimulation.1<- revalue(raw.data$Stimulation.1, c("em parte"= 4, "exactamente"= 5, "n達o influenciou"=3) )
raw.data$Stimulation.1<-as.numeric(as.character(raw.data$Stimulation.1))

#Stimulation.2
levels(raw.data$Stimulation.2)
raw.data$Stimulation.2<- revalue(raw.data$Stimulation.2, c("em parte"= 4, "exactamente"= 5) )
raw.data$Stimulation.2<-as.numeric(as.character(raw.data$Stimulation.2))

#Routine.1
levels(raw.data$Routine.1)
raw.data$Routine.1<- revalue(raw.data$Routine.1, c("nunca"= 1) )
raw.data$Routine.1<-as.numeric(as.character(raw.data$Routine.1))

#Achievement.1
levels(raw.data$Achievement.1)
raw.data$Achievement.1<- revalue(raw.data$Achievement.1, c("sim" = 5, "sim, um pouco"=4) )
raw.data$Achievement.1<-as.numeric(as.character(raw.data$Achievement.1))

#Achievement.2
levels(raw.data$Achievement.2)
raw.data$Achievement.2<- revalue(raw.data$Achievement.2, c("sim" = 5) )
raw.data$Achievement.2<-as.numeric(as.character(raw.data$Achievement.2))

#Power.1
levels(raw.data$Power.1)
raw.data$Power.1<- revalue(raw.data$Power.1, c("n達o muito" = 2,"n達o, de todo" = 1, "sim"= 5, "sim, um pouco" = 4 ) )
raw.data$Power.1<-as.numeric(as.character(raw.data$Power.1))

#Power.2
levels(raw.data$Power.2)
raw.data$Power.2<- revalue(raw.data$Power.2, c("nada" = 1) )
raw.data$Power.2<-as.numeric(as.character(raw.data$Power.2))

#Belongingness.1
levels(raw.data$Belongingness.1)
raw.data$Belongingness.1<- revalue(raw.data$Belongingness.1, c("muito influenciada" = 5 , "n達o, de todo" = 1, "neutro" = 3))
raw.data$Belongingness.1<-as.numeric(as.character(raw.data$Belongingness.1))


#Belongingness.2
levels(raw.data$Belongingness.2)
raw.data$Belongingness.2<- revalue(raw.data$Belongingness.2, c("sim" = 5))
raw.data$Belongingness.2<-as.numeric(as.character(raw.data$Belongingness.2))

#Conformity.1
levels(raw.data$Conformity.1)
raw.data$Conformity.1<- revalue(raw.data$Conformity.1, c("algumas pessoas" = 4, "poucos participantes" = 3))
raw.data$Conformity.1<-as.numeric(as.character(raw.data$Conformity.1))

#Benevolence.2
levels(raw.data$Benevolence.2)
raw.data$Benevolence.2<- revalue(raw.data$Benevolence.2, c("Sim, definitivamente" = 5, "Sim, principalmente por isso" =4))
raw.data$Benevolence.1<-as.numeric(as.character(raw.data$Benevolence.1))

#Universalism.1
levels(raw.data$Universalism.1)
raw.data$Universalism.1<- revalue(raw.data$Universalism.1, c("Sim, definitivamente" = 5))
raw.data$Universalism.1<-as.numeric(as.character(raw.data$Universalism.1))

write.csv(raw.data, "water-sentinels-pivot-questions.csv", row.names = F)


#--------------------------- ANALYSIS RESULTS BY QUESTION

# Compute average and variance for each question 
ans<- read.csv("water-sentinels-pivot-questions.csv", header=TRUE)

mean.question<-round(sapply(ans[,c(1:21)], mean),2)
var.questions<- round(sapply(ans[,c(1:21)], var), 2)

df.questions<- data.frame(mean= mean.question, var = var.questions)


#--------------------------- ANALYSIS RESULTS BY MOTIVATING FACTOR

# Compute average and variance for each motivating factor 
ans<- read.csv("water-sentinels-pivot-questions.csv", header=TRUE)

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
  factors.summary<-rbind(factors.summary, new.row) #aggiunta di riga al dataframe
  
  rm(f, f.mean, f.var, new.row, v.final, v, num.col)
  
}


#--------------------------- CORRELATION ANALYSIS

completed<- read.csv("water-sentinels-pivot-questions.csv", header=TRUE)

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
for(i in 22:31){ # Change considering dataframe size
  t <- cor.test(d.cor[,i], d.cor$Global.motivation)
  nome <- names(d.cor)[i]
  corr <- t$estimate
  pv <- t$p.value
  sign <- ifelse(pv<0.001, "***", ifelse(pv<0.01, "**", ifelse(pv<0.05, "*", "")))
  global.motivation<- paste0(round(corr,3),sign)
  global.corr <- rbind(global.corr, data.frame(nome, corr, pv, sign, global.motivation))
}
### Warning! Standard Dev is zero because global motivation is 5 for each survey completion