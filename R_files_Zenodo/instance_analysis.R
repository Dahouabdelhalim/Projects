### Script to analyze InStance data

#instance <- read.table("C:/Users/ebaykara/Documents/InStance_Questionnaire/data_s4hri_2018-05-17_14-57.txt", header=TRUE, sep="\\t",dec=".", na.strings = " ",quote = "")

### This code reads the file with the new data, name it differently than the main data file
instq <- read.csv("C:/Users/ebaykara/Documents/InStance_Questionnaire/data_s4hri_2018-05-17_14-57.csv", header=TRUE, sep=",",dec=".", na.strings = " ")
new_test <- read.csv("C:/Users/ebaykara/Documents/InStance_Questionnaire/data_s4hri_2018-05-21_10-23.csv", header=TRUE, sep=",",dec=".", na.strings = " ")

## Run the function files to include them in the Global Environment
subt <- function(x){  ##subt is to subtract 1 from the scoring of each question (1-101->0-100)
  
  x <- x-1
  return (x)

}

rev <- function(x){   ##Reverse the scoring for questions 19-35, where intentional and mechanistic are reversed in the questionnaire.
  
  x <- 100-x
  return (x)
  
}
  

ff <- dput(grep("^IS.._..$", names(new_test), value = TRUE)) # two letter variable names

qq <- data.frame(mapply(subt, new_test[,ff[1:35]]))
rr <- mapply(rev, qq[,ff[19:35]])
qdata <- cbind(qq[,1:18],rr)

new_test <- new_test[ , -which(names(new_test) %in% ff[1:35])]  ## Remove the original scoring data (1-101, reverse ordered)
new_test <- cbind(new_test,qdata)  ## Add the new calculated scores for each question

instance$Avg_Score <- rowMeans(instance[,ff[2:35]])  ### Calculate average score for each participant


instance$Fam <- ifelse(instance$IS43==1,1,2)  ## Code Familiarity as 1=No, 2=yes


### Coding of the gender 1=F, 2=M
instance$Gender[instance$IS37_02 == "f"] <- 1
instance$Gender[instance$IS37_02 == "F"] <- 1
instance$Gender[instance$IS37_02 == "m"] <- 2
instance$Gender[instance$IS37_02 == "M"] <- 2


### Analysis 

fit <- lm(Avg_Score~Fam, data=instance)
summary(fit)

aggregate(test$Avg_Score, by=list(test$Fam), FUN=mean)
length(which(instance$Fam==1))
length(which(instance$Fam==2))

fit <- lm(Avg_Score~Gender + IS37_03 + IS41_01 + IS37_04 + IS37_05, data=NonFam_cl)
summary(fit)

### Create a subset of participants who are not familiar with robots
NonFam <- subset(test,Fam==1)

cor.test(NonFam$Avg_Score,NonFam$MindEyes_Correct,method="spearman")
cor.test(NonFam$Avg_Score,NonFam$Vollm_Correct,method="spearman")

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

NonFam$Group[NonFam$Avg_Score > 50] <- 2
NonFam$Group[NonFam$Avg_Score < 50] <- 1

cs <- table(NonFam$Group)
chisq.test(cs)



### Compare score of different groups to 50 (0 of our scale)
aggregate(x_notfam$Avg_Score, by=list(x_notfam$Group), FUN=mean)
G1 <- NonFam$Avg_Score[NonFam$Group==1]
G2 <- NonFam$Avg_Score[NonFam$Group==2]

t.test(G1,mu=50)
t.test(G2,mu=50)

### Standard error of means
sd(G1, na.rm=TRUE) / sqrt(length(G1[!is.na(G1)])) 
sd(G2, na.rm=TRUE) / sqrt(length(G2[!is.na(G2)])) 
