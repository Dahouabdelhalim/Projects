#Analyses run with R4.1.0
library(tidyverse)
library(dplyr)
library(afex)
library(psych)
library(emmeans)
library(effectsize)

#Experiment 1 Analysis
Exp1_Stroop <- read.table("Exp1_with_filter.csv", head=TRUE, sep=",")
str(Exp1_Stroop)
Exp1_Stroop$Participant <- factor(Exp1_Stroop$Participant)
Exp1_Stroop$Congruency <- factor(Exp1_Stroop$Congruency)
#Descriptive Statistics
#mean and sd for RT Trimmed
describeBy(Exp1_Stroop$RTTrimmed, Exp1_Stroop$Congruency)
#count of Errors **note the number of errors does not include those removed by trimming**
group_by(Exp1_Stroop, Congruency) %>% 
  dplyr::summarise(sumNA = sum(is.na(RT))) #this shows errors before RTs excluded

#RT Analysis
Exp1 <-group_by(Exp1_Stroop, Participant, Congruency) %>%
  dplyr::summarise(n = n(),
                   mean = mean(RTTrimmed, na.rm = TRUE),
                   sd = sd(RTTrimmed, na.rm = TRUE)
  )
Exp1.Comp <- t.test(mean ~ Congruency, data = Exp1, paired = TRUE, var.equal = TRUE)
Exp1.Comp
cohens_d(mean ~ Congruency, data = Exp1, paired = TRUE)

#Accuracy Analysis
Exp1_Err <-group_by(Exp1_Stroop, Participant, Congruency) %>%
  dplyr::summarise(n = n(),
                   mean = mean(Correct, na.rm = TRUE),
                   sum = sum(Correct),
  )
Exp1.ErrComp <- t.test(sum() ~ Congruency, data = Exp1_Err, paired = TRUE)
Exp1.ErrComp
cohens_d(sum ~ Congruency, data = Exp1_Err, paired = TRUE)

#Experiment 2 Analysis
CombExp <- read.csv("Exp2+3_Stroop.csv")
Exp2_nofb <- filter(CombExp, Experiment =="2")
str(Exp2_nofb)
Exp2_nofb$Participant <- factor(Exp2_nofb$Participant)
Exp2_nofb$PairCondition_reduced <- factor(Exp2_nofb$PairCondition_reduced)
Exp2_nofb$Congruency <- factor(Exp2_nofb$Congruency)

#Descriptive Statistics
#mean and sd for RT Trimmed
describeBy(Exp2_nofb$TrimmedRT, list(Exp2_nofb$Congruency, Exp2_nofb$PairCondition_reduced))
#descriptive statistics for errors
sum(is.na(Exp2_nofb$RT))
group_by(Exp2_nofb, PairCondition_reduced, Congruency) %>% 
  dplyr::summarise(sumNA = sum(is.na(RT))) #this shows errors before RTs excluded
group_by(Exp2_nofb, PairCondition_reduced, Congruency) %>% 
  dplyr::summarise(sumNA = sum(is.na(TrimmedRT))) #this is overall and includes exclusions
#graph
nofb.BarGraph<-ggplot(Exp2_nofb, aes(Congruency, TrimmedRT, fill=Congruency)) +
  geom_bar(stat="summary", fun="mean", position = "dodge") + 
  facet_grid(.~PairCondition_reduced) +
  xlab("Condition") + ylab("Reaction Time") +
  scale_fill_brewer(palette="Dark2") +
  theme(legend.position="none")
nofb.BarGraph
#with error bars
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
data2 <- data_summary(nofbstroop, varname="TrimmedRT", 
                    groupnames=c("PairCondition_reduced", "Condition"))
head(data2)
ggplot(data2, aes(x=Condition, y=TrimmedRT, fill=Condition)) + 
  geom_bar(stat="summary", fun = "mean", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=TrimmedRT-sd, ymax=TrimmedRT+sd), width=.2,
                position=position_dodge(.9))+
  facet_grid(.~PairCondition_reduced) +
  xlab("Condition") + ylab("Reaction Time") +
  theme(legend.position="none") +
  theme_minimal()

#ANOVA
#RT Analysis
# Exclude the missing observations
Exp2_nofb_drop <-Exp2_nofb %>%
  na.omit()		
dim(Exp2_nofb_drop)
# IV (between): PairCondition_reduced
# IV (within): Condition
# DV:          TrimmedRT
Exp2_nofb_drop$Participant <- factor(Exp2_nofb_drop$Participant)
Exp2_nofb_drop$PairCondition_reduced <- factor(Exp2_nofb_drop$PairCondition_reduced)
Exp2_nofb_drop$Congruency <- factor(Exp2_nofb_drop$Congruency)

aov_nofb <- aov_car(TrimmedRT ~ PairCondition_reduced*Congruency + Error(Participant/Congruency), data=Exp2_nofb_drop, anova_table = list(es = "pes"))
nice(aov_nofb)
emmeans(aov_nofb, specs = pairwise ~ Congruency|PairCondition_reduced, type = "response")

#Accuracy Analysis - uses the data with NAs included!
aov_nofb_err <- aov_car(Errors ~ PairCondition_reduced*Congruency + Error(Participant/Congruency), data=Exp2_nofb, anova_table = list(es = "pes"))
nice(aov_nofb_err)
emmeans(aov_nofb_err, specs = pairwise ~ Congruency|PairCondition_reduced, type = "response")

#Experiment 3 Analysis
Exp3_fb <- filter(CombExp, Experiment =="3")
str(Exp3_fb)
Exp3_fb$Participant <- factor(Exp3_fb$Participant)
Exp3_fb$PairCondition_reduced <- factor(Exp3_fb$PairCondition_reduced)
Exp3_fb$Congruency <- factor(Exp3_fb$Congruency)

#Descriptive Statistics
#mean and sd for RT Trimmed
describeBy(Exp3_fb$TrimmedRT, list(Exp3_fb$Congruency, Exp3_fb$PairCondition_reduced))
#descriptive statistics for errors
group_by(Exp3_fb, PairCondition_reduced, Congruency) %>% 
  dplyr::summarise(sumNA = sum(is.na(RT))) #this shows errors before RTs excluded
group_by(Exp3_fb, PairCondition_reduced, Congruency) %>% 
  dplyr::summarise(sumNA = sum(is.na(TrimmedRT))) #this is overall, and includes exclusions
#graph
fb.BarGraph<-ggplot(Exp3_fb, aes(Congruency, TrimmedRT, fill=Congruency)) +
  geom_bar(stat="summary", fun="mean", position = "dodge") + 
  facet_grid(.~PairCondition_reduced) +
  xlab("Condition") + ylab("Reaction Time") +
  scale_fill_brewer(palette="Dark2") +
  theme(legend.position="none")
fb.BarGraph
#with error bars
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
data2 <- data_summary(Exp3_fb, varname="TrimmedRT", 
                      groupnames=c("PairCondition_reduced", "Congruency"))
head(data2)
ggplot(data2, aes(x=Congruency, y=TrimmedRT, fill=Congruency)) + 
  geom_bar(stat="summary", fun = "mean", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=TrimmedRT-sd, ymax=TrimmedRT+sd), width=.2,
                position=position_dodge(.9))+
  facet_grid(.~PairCondition_reduced) +
  xlab("Condition") + ylab("Reaction Time") +
  theme(legend.position="none") +
  theme_minimal()

#ANOVA
#RT Analysis
# Exclude the missing observations
Exp3_fb_drop <-Exp3_fb %>%
  na.omit()		
dim(Exp3_fb_drop)
# IV (between): PairCondition_reduced
# IV (within): Condition
# DV:          TrimmedRT
Exp3_fb_drop$Participant <- factor(Exp3_fb_drop$Participant)
Exp3_fb_drop$PairCondition_reduced <- factor(Exp3_fb_drop$PairCondition_reduced)
Exp3_fb_drop$Congruency <- factor(Exp3_fb_drop$Congruency)

aov_fb <- aov_car(TrimmedRT ~ PairCondition_reduced*Congruency + Error(Participant/Congruency), data=Exp3_fb_drop, anova_table = list(es = "pes"))
nice(aov_fb)
emmeans(aov_fb, specs = pairwise ~ Congruency|PairCondition_reduced, type = "response")

#Accuracy Analysis - uses the data with NAs included!
aov_fb_err <- aov_car(Errors ~ PairCondition_reduced*Congruency + Error(Participant/Congruency), data=Exp3_fb, anova_table = list(es = "pes"))
nice(aov_fb_err)
emmeans(aov_fb_err, specs = pairwise ~ Congruency|PairCondition_reduced, type = "response")

#CombinedAnalysis
#RT Analysis
# Exclude the missing observations
CombExp_drop <-CombExp %>%
  na.omit()		
dim(CombExp_drop)
# IV (between): Experiment
# IV (between): PairCondition_reduced
# IV (within): Condition
# DV:          TrimmedRT
CombExp_drop$Experiment <- factor(CombExp_drop$Experiment)
CombExp_drop$Participant <- factor(CombExp_drop$Participant)
CombExp_drop$PairCondition_reduced <- factor(CombExp_drop$PairCondition_reduced)
CombExp_drop$Congruency <- factor(CombExp_drop$Congruency)

aov_both <- aov_car(TrimmedRT ~ Experiment*PairCondition_reduced*Congruency + Error(Participant/Congruency), data=CombExp_drop, anova_table = list(es = "pes"))
nice(aov_both)

#Accuracy Analysis - uses the data with NAs included!
aov_both_err <- aov_car(Errors ~ Experiment*PairCondition_reduced*Congruency + Error(Participant/Congruency), data=CombExp, anova_table = list(es = "pes"))
nice(aov_both_err)
