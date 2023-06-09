instance_nov2018<- read.table("//IITSCHRIWS001/S4HRI_localrepo/Experiments/InStance_Questionnaire/instance_Nov2018.txt", header=TRUE, sep="\\t",dec=",", na.strings = " ",quote = "")
View(instance_nov2018)

mean(instance_nov2018$IS37_03)#mean age
sd(instance_nov2018$IS37_03)
min(instance_nov2018$IS37_03)
max(instance_nov2018$IS37_03)
length(which(instance_nov2018$Gender==1))
mean(instance_nov2018$IS41_01)#mean Education
sd(instance_nov2018$IS41_01)
min(instance_nov2018$IS41_01)
max(instance_nov2018$IS41_01)
### Analysis (sample =106)
shapiro.test(instance_nov2018$Avg_Score)
mean(instance_nov2018$Avg_Score)
sd(instance_nov2018$Avg_Score)
ISS_106<-instance_nov2018$Avg_Score
t.test(ISS_106,mu=0)
#Multiple regression for Gender, Age, Education, Children, Siblings
fit2 <- lm(Avg_Score~Gender + IS37_03 + IS41_01 + IS37_04 + IS37_05, data=instance_nov2018)
summary(fit2)
  #Correlation with MindEyes and Vollm
  
cor.test(instance_nov2018$Avg_Score, instance_nov2018$MindEyes_Correct,method="spearman")
cor.test(instance_nov2018$Avg_Score,instance_nov2018$Vollm_Correct,method="spearman")


#linear regression for familiarity
length(which(instance_nov2018$Fam==1))
length(which(instance_nov2018$Fam==2))
fit <- lm(Avg_Score~Fam, data=instance_nov2018)
summary(fit)


### Create a subset of participants who are not familiar with robots
NonFam_nov2018 <- subset(instance_nov2018, Fam==1)
Fam_nov2018 <- subset(instance_nov2018, Fam==2)
mean(NonFam_nov2018$Avg_Score)
mean(Fam_nov2018$Avg_Score)
sd(NonFam_nov2018$Avg_Score)
sd(Fam_nov2018$Avg_Score)
ISS_Fam<-Fam_nov2018$Avg_Score
t.test(ISS_Fam,mu=0)
ISS_NonFam<-NonFam_nov2018$Avg_Score
t.test(ISS_NonFam, mu=0)
### Standard error of means

perc_NonFam<-length(which(instance_nov2018$Fam==1))/106
perc_Fam<-length(which(instance_nov2018$Fam==2))/106


### Analysis on Non Familiar (sample =106)

shapiro.test(NonFam_nov2018$Avg_Score)

#Multiple regression for Gender, Age, Education, Children, Siblings
NonFam_nov2018_fit2 <- lm(Avg_Score~Gender + IS37_03 + IS41_01 + IS37_04 + IS37_05, data=NonFam_nov2018)
summary(NonFam_nov2018_fit2)

cor.test(NonFam_nov2018$Avg_Score,NonFam_nov2018$MindEyes_Correct,method="spearman")
cor.test(NonFam_nov2018$Avg_Score,NonFam_nov2018$Vollm_Correct,method="spearman")

## Gender
fit <- lm(Avg_Score~Gender, data=NonFam_cl)
summary(fit)

## Age
fit <- lm(Avg_Score~IS37_03, data=NonFam_cl)
summary(fit)

## Education
fit <- lm(Avg_Score~IS41_01, data=NonFam_cl)
summary(fit)

### Children
fit <- lm(Avg_Score~IS37_04, data=NonFam_cl)
summary(fit)

### Siblings
fit <- lm(Avg_Score~IS37_05, data=NonFam_cl)
summary(fit)

### Chi-Square of the distribution of mechanistic vs intentional answers

NonFam_nov2018$Group[NonFam_nov2018$Avg_Score > 50] <- 2
NonFam_nov2018$Group[NonFam_nov2018$Avg_Score < 50] <- 1

cs <- table(NonFam_nov2018$Group)
chisq.test(cs)



### Compare score of different groups to 50 (0 of our scale)
#Not familiar group
aggregate(x_notfam$Avg_Score, by=list(x_notfam$Group), FUN=mean)
G1 <- NonFam_nov2018$Avg_Score[NonFam_nov2018$Group==1]
G2 <- NonFam_nov2018$Avg_Score[NonFam_nov2018$Group==2]
length(which(NonFam_nov2018$Group==1))
length(which(NonFam_nov2018$Group==2))
perc_NonFam_mec<-length(which(NonFam_nov2018$Group==1))/89
perc_NonFam_ment<-length(which(NonFam_nov2018$Group==2))/89


t.test(G1,mu=50)
t.test(G2,mu=50)
#Frequencies Familar group
Fam_nov2018$Group[Fam_nov2018$Avg_Score > 50] <- 2
Fam_nov2018$Group[Fam_nov2018$Avg_Score < 50] <- 1
aggregate(x_fam$Avg_Score, by=list(x_fam$Group), FUN=mean)
G1_Fam_nov_2018 <- Fam_nov2018$Avg_Score[Fam_nov2018$Group==1]
G2_Fam_nov2018 <- Fam_nov2018$Avg_Score[Fam_nov2018$Group==2]
length(which(Fam_nov2018$Group==1))
length(which(Fam_nov2018$Group==2))

perc_Fam_mec<-length(which(Fam_nov2018$Group==1))/17
perc_Fam_ment<-length(which(Fam_nov2018$Group==2))/17

### Standard error of means
sd(G1, na.rm=TRUE) / sqrt(length(G1[!is.na(G1)])) 
sd(G2, na.rm=TRUE) / sqrt(length(G2[!is.na(G2)])) 

#plots
hist(instance_nov2018$Avg_Score, main="InStance Scores", xlab="Scores", border="black", col="grey",xlim=c(0,100),ylim=c(0,0.04),las=1, 
     breaks=25, prob = TRUE)
lines(density(instance_nov2018$Avg_Score))
hist(NonFam_nov2018$Avg_Score, main="InStance Scores: Not familiar group", xlab="Scores", border="black", col="grey",xlim=c(0,100),ylim=c(0,0.04),las=1, 
     breaks=25, prob = TRUE)
lines(density(NonFam_nov2018$Avg_Score))

NonFamMec_nov2018 <- subset(NonFam_nov2018, NonFam_nov2018$Group==1)
NonFamMent_nov2018 <- subset(NonFam_nov2018, NonFam_nov2018$Group==2)

hist(NonFamMec_nov2018$Avg_Score, main="InStance Scores:Non Familiar InStance Mechanistic Group", xlab="Scores", border="black", col="grey",xlim=c(0,100),ylim=c(0,0.06),las=1, breaks=25, prob = TRUE)
lines(density(NonFamMec_nov2018$Avg_Score))
hist(NonFamMent_nov2018$Avg_Score, main="InStance Scores:Not Familiar InStance Mentalistic Gruop", xlab="Scores", border="black", col="grey",xlim=c(0,100),ylim=c(0,0.25),las=1, 
     breaks=25, prob = TRUE)
lines(density(NonFamMent_nov2018$Avg_Score))

#test per Bimodal distribution
library(diptest)
dip(NonFam_nov2018$Avg_Score)#result=0.03133407
dip(NonFam_nov2018_mec$Avg_Score) #result=0.03430837
dip(NonFamMent_nov2018$Avg_Score)#result=0.07412467
dip(instance_nov2018$Avg_Score) #result=0.0265179
shapiro.test(NonFamMec_nov2018$Avg_Score)#W = 0.96418, p-value = 0.06735
shapiro.test(NonFamMent_nov2018$Avg_Score)#W = 0.93106, p-value = 0.07339