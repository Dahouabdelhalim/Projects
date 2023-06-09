setwd("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis")
library(lme4); library(lmerTest); library(stats); library(ggplot2); library(readr); library(readxl)
library(Rmisc); library(pastecs); library(psych); library(multcomp)


####------------------------ AUs selection process DO NOT USE!!!!!------------------------####
## PRODUCTION
FactorProduction <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/0_AU selection/AU selected_Production_PCA.xlsx", 
                               sheet = "Aggregated")
View(FactorProduction)
SelectionGuilt<-  glm(GuiltSelf_P1_T3 ~ FAC1_1 + FAC2_1 + FAC3_1 + FAC4_1 + FAC5_1 + FAC6_1 + FAC7_1 + FAC8_1, dat=FactorProduction)
summary(SelectionGuilt)

## PERCEPTION
FactorPerception <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/0_AU selection/AU selected_Perception_PCA.xlsx", 
                               sheet = "Aggregated2")
View(FactorPerception)
SelectionGuilt<-  glm(GuiltJudged_P2P1 ~ FAC1_1 + FAC2_1 + FAC3_1, dat=FactorPerception)
summary(SelectionGuilt)

####------------------------------ UK stay ----------------------------####
UK <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/Personality_summary.xlsx")
t.test(UK$Gp, UK$UK_stay)

####---------------------------STUDY A---------------------------------####
###---------------- AFFECT -------------------------------###
AffectA <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyA/Data_StudyA_Affect.xlsx")
View(AffectA)

t.test(AffectA$`1PositiveAffect`,AffectA$`3PositiveAffectAfter`, paired = TRUE) # Pos Panas 1 vs 3
t.test(AffectA$`1NegativeAffect`,AffectA$`3NegativeAffectAfter`, paired = TRUE) # Neg Panas 1 vs 3
t.test(AffectA$`1Guilty`,AffectA$`3Guilty`, paired = TRUE) # Guilty Panas 1 vs 3
t.test(AffectA$`1Ashamed`,AffectA$`3Ashamed`, paired = TRUE)
t.test(AffectA$`1Distressed`,AffectA$`3Distressed`, paired = TRUE)
t.test(AffectA$`1Proud`,AffectA$`3Proud`)

t.test(AffectA$`2PositiveAffect`,AffectA$`3PositiveAffectAfter`, paired = TRUE) # Pos Panas 2 vs 3
t.test(AffectA$`2NegativeAffect`,AffectA$`3NegativeAffectAfter`, paired = TRUE) # Neg Panas 2 vs 3
t.test(AffectA$`2Guilty`,AffectA$`3Guilty`, paired = TRUE) # Guilty Panas 2 vs 3
t.test(AffectA$`2Ashamed`,AffectA$`3Ashamed`, paired = TRUE)
t.test(AffectA$`2Distressed`,AffectA$`3Distressed`, paired = TRUE)
t.test(AffectA$`2Proud`,AffectA$`3Proud`, paired = TRUE)

t.test(AffectA$`1PositiveAffect`,AffectA$`2PositiveAffect`, paired = TRUE) # Pos Panas 1 vs 2
t.test(AffectA$`1NegativeAffect`,AffectA$`2NegativeAffect`, paired = TRUE) # Neg Panas 1 vs 2
t.test(AffectA$`1Guilty`,AffectA$`2Guilty`, paired = TRUE) # Guilty Panas 1 vs 2
t.test(AffectA$`1Ashamed`,AffectA$`2Ashamed`, paired = TRUE)
t.test(AffectA$`1Distressed`,AffectA$`2Distressed`, paired = TRUE)
t.test(AffectA$`1Proud`,AffectA$`2Proud`, paired = TRUE)

t.test(AffectA$`3Ashamed`, AffectA$`3Guilty`, paired = TRUE)

# Affect plot
Affect_Plot_12 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyA/Data_StudyA_Affect_Plot_12.xlsx")
Sum <- summarySE(Affect_Plot_12, measurevar="Score", groupvars=c("Affect")) # PANAS 1 vs 2
ggplot(Sum, aes(x=factor(Affect), y=Score, group = 1))  + 
  geom_point(size = 0.7) + 
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_errorbar(aes(ymin=Score-se, ymax=Score+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white'))

Affect_Plot_23 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyA/Data_StudyA_Affect_Plot_23.xlsx")
Sum <- summarySE(Affect_Plot_23, measurevar="Score", groupvars=c("Affect")) # PANAS 2 vs 3
ggplot(Sum, aes(x=factor(Affect), y=Score, group = 1))  + 
  geom_point(size = 0.7) + 
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_errorbar(aes(ymin=Score-se, ymax=Score+se), size = .6, width=.2, alpha=0.8) + 
  xlab("Self-reported affect") + ylab("Change before/after feedback") +
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white'))

###----------------- GENERAL RESULTS Study A -----------------------####
StudyA <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyA/Data_StudyA.xlsx")
View(StudyA)

h <- hist(StudyA$Friend_index, xlab = "Friendship index", main = "Histogram with normal curve")
xfit <- seq (min(StudyA$Friend_index), max(StudyA$Friend_index), length = 40)
yfit <- dnorm(xfit, mean = mean(StudyA$Friend_index), sd = sd(StudyA$Friend_index))
yfit <- yfit*diff(h$mids[1:2])*length(StudyA$Friend_index)
lines(xfit, yfit, col = "blue", lwd = 2)
qplot(sample = StudyA$Friend_index, stat = "qq")

t.test(StudyA$`YOU`, StudyA$`OTHER`, paired = TRUE)
describe(StudyA)

summary(glm(Percent_Other ~ GuiltSelf_P1_T3, dat=StudyA))
money_split <- glm(Percent_Other ~ GuiltSelf_P1_T3 + Friend_index + FAC4_2 + Congruence2 + ConditionP1, dat= StudyA)
summary(money_split)

#summary(glm(Percent_Other ~ FAC1_1 + FAC2_1 + FAC3_1 + FAC4_1 + FAC5_1 + FAC6_1 + FAC7_1 + FAC8_1,  dat=StudyA))
#summary(glm(GuiltSelf_P1_T3 ~ FAC1_1 + FAC2_1 + FAC3_1 + FAC4_1 + FAC5_1 + FAC6_1 + FAC7_1 + FAC8_1,  dat=StudyA))
summary(glm(GuiltSelf_P1_T3 ~ Friend_index, dat=StudyA))
summary(glm(GuiltSelf_P1_T3 ~ ConditionP1 + EthnicGp, dat=StudyA))
summary(glm(JudgedGuilt_P2P1 ~ GuiltSelf_P1_T3, dat=StudyA))
summary(glm(Friend_index ~ Congruence + EthnicGp, dat=StudyA))
cor.test(StudyA$GuiltSelf_P1_T3,StudyA$JudgedGuilt_P2P1) #ouch...

Sum <- summarySE(StudyA, measurevar="Percent_Other", groupvars=c("GuiltSelf_P1_BI"))
ggplot(Sum, aes(x=factor(GuiltSelf_P1_BI), y=Percent_Other, group = 1))  + 
  geom_point(size = 0.7) + 
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_smooth(method=lm, se=T, fullrange=F, size=.3, colour="black", linetype="dashed") + 
  geom_errorbar(aes(ymin=Percent_Other-se, ymax=Percent_Other+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white')) + 
  scale_y_continuous(limits=c(0, 1)) + 
  xlab("Guilt Felt") + ylab("Reward split") +
  geom_hline(yintercept = .5, linetype="dotted", alpha=0.7) +
  geom_hline(yintercept = .466, linetype-"dotted", alpha=0.7, color ="red") +
  geom_hline(yintercept = .533, linetype="dashed", alpha=0.7, color = "red")
  

Sum <- summarySE(StudyA, measurevar="You_Other", groupvars=c("GuiltSelf_P1_BI"))
ggplot(Sum, aes(x=factor(GuiltSelf_P1_BI), y=You_Other, group = 1))  + 
  geom_point(size = 0.7) + 
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_smooth(method=lm, se=T, fullrange=F, size=.3, colour="black", linetype="dashed") + 
  geom_errorbar(aes(ymin=You_Other -se, ymax=You_Other +se), size = .6, width=.2, alpha=0.8) +
  xlab("Self-reported Guilt") + ylab("Number of extra coins kept for self") +
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white')) 

Sum <- summarySE(StudyA, measurevar="Percent_Other", groupvars=c("Congruence2"))
ggplot(Sum, aes(x=factor(Congruence2), y=Percent_Other, group = 1))  + 
  geom_point(size = 0.7) + 
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  #geom_smooth(method=lm, se=T, fullrange=F, size=.3, colour="black", linetype="dashed") + 
  geom_errorbar(aes(ymin=Percent_Other-se, ymax=Percent_Other+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white')) + 
  scale_y_continuous(limits=c(0, 1)) + 
  xlab("Pair - Culture") + ylab("Reward split") +
  geom_hline(yintercept = .5, linetype="dotted", alpha=0.7)+
  geom_hline(yintercept = .466, linetype-"dotted", alpha=0.7, color ="red") +
  geom_hline(yintercept = .533, linetype="dashed", alpha=0.7, color = "red")

#--------------Impact of GASP on Money Split-------------#
summary(lm(Percent_Other ~ G_R, dat= StudyA))
summary(lm(Percent_Other ~ G_NBE, dat= StudyA))
summary(lm(Percent_Other ~ S_W, dat= StudyA))
summary(lm(Percent_Other ~ S_NSE, dat= StudyA))

#--------------Impact of GASP on Felt Guilt-------------#
summary(lm(GuiltSelf_P1_T3 ~ G_R,  dat=StudyA))
summary(lm(GuiltSelf_P1_T3 ~ G_NBE, dat=StudyA))
summary(lm(GuiltSelf_P1_T3 ~ S_W, dat=StudyA))
summary(lm(GuiltSelf_P1_T3 ~ S_NSE, dat=StudyA))
cor.test(Personality$G_NBE,Personality$S_NSE)


###----------------- PERSONALITY DATA ---------------------####
PersonalityA <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyA/Data_StudyA_Personality.xlsx")
View(PersonalityA)

#Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism

Money <- glm(Money_split ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityA)
summary(Money)
summary(glm(Money_split ~ AnnoyanceSelf, dat=PersonalityA))
summary(glm(Money_split ~ AnnoyancePartner, dat=PersonalityA))
summary(glm(Money_split ~ WhoseFault, dat=PersonalityA))

Guilt <- glm(Guilt3 ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityA)
summary(Guilt)
summary(glm(AnnoyanceSelf ~ Guilt3, dat = PersonalityA))

Annoyance <- glm(AnnoyanceSelf ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityA)
summary(Annoyance)


####---------------------------STUDY B---------------------------------####
###---------------- AFFECT -------------------------------###
AffectB <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB_Affect.xlsx")
View(AffectB)

t.test(AffectB$`1PositiveAffect`,AffectB$`3PositiveAffectAfter`, paired = TRUE) # Pos Panas 1 vs 3
t.test(AffectB$`1NegativeAffect`,AffectB$`3NegativeAffectAfter`, paired = TRUE) # Neg Panas 1 vs 3
t.test(AffectB$`1Guilty`,AffectB$`3Guilty`, paired = TRUE) # Guilty Panas 1 vs 3
t.test(AffectB$`1Ashamed`,AffectB$`3Ashamed`, paired = TRUE)
t.test(AffectB$`1Distressed`,AffectB$`3Distressed`, paired = TRUE)
t.test(AffectB$`1Proud`,AffectB$`3Proud`, paired = TRUE)

t.test(AffectB$`2PositiveAffect`,AffectB$`3PositiveAffectAfter`, paired = TRUE) # Pos Panas 2 vs 3
t.test(AffectB$`2NegativeAffect`,AffectB$`3NegativeAffectAfter`, paired = TRUE) # Neg Panas 2 vs 3
t.test(AffectB$`2Guilty`,AffectB$`3Guilty`, paired = TRUE) # Guilty Panas 2 vs 3
t.test(AffectB$`2Ashamed`,AffectB$`3Ashamed`, paired = TRUE)
t.test(AffectB$`2Distressed`,AffectB$`3Distressed`, paired = TRUE)
t.test(AffectB$`2Proud`,AffectB$`3Proud`, paired = TRUE)

t.test(AffectB$`1PositiveAffect`,AffectB$`2PositiveAffect`, paired = TRUE) # Pos Panas 1 vs 2
t.test(AffectB$`1NegativeAffect`,AffectB$`2NegativeAffect`, paired = TRUE) # Neg Panas 1 vs 2
t.test(AffectB$`1Guilty`,AffectB$`2Guilty`, paired = TRUE) # Guilty Panas 1 vs 2
t.test(AffectB$`1Ashamed`,AffectB$`2Ashamed`, paired = TRUE)
t.test(AffectB$`1Distressed`,AffectB$`2Distressed`, paired = TRUE)
t.test(AffectB$`1Proud`,AffectB$`2Proud`, paired = TRUE)

t.test(AffectB$`3Ashamed`, AffectB$`3Guilty`, paired = TRUE)

# Affect plot
Affect_Plot_12 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_Studyb_Affect_Plot_12.xlsx")
Sum <- summarySE(Affect_Plot_12, measurevar="Score", groupvars=c("Affect")) # PANAS 1 vs 2
ggplot(Sum, aes(x=factor(Affect), y=Score, group = 1))  + geom_point(size = 0.7) + geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + geom_errorbar(aes(ymin=Score-se, ymax=Score+se), size = .6, width=.2, alpha=0.8) + theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white'))

Affect_Plot_23 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB_Affect_Plot_23.xlsx")
Sum <- summarySE(Affect_Plot_23, measurevar="Score", groupvars=c("Affect")) # PANAS 2 vs 3
ggplot(Sum, aes(x=factor(Affect), y=Score, group = 1))  + 
  geom_point(size = 0.7) + geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_errorbar(aes(ymin=Score-se, ymax=Score+se), size = .6, width=.2, alpha=0.8) + 
  xlab("Self-reported affect") + ylab("Change before/after feedback") +
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white'))

###----------------- GENERAL RESULTS Study B -----------------------####
StudyB <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB.xlsx", sheet = "All_Factors")
View(StudyB)

StudyBall <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB.xlsx", sheet = "All_condition")

h2 <- hist(StudyB$Friend_index, xlab = "Friendship index", main = "Histogram with normal curve")
xfit <- seq (min(StudyB$Friend_index), max(StudyB$Friend_index), length = 40)
yfit <- dnorm(xfit, mean = mean(StudyB$Friend_index), sd = sd(StudyB$Friend_index))
yfit <- yfit*diff(h2$mids[1:2])*length(StudyB$Friend_index)
lines(xfit, yfit, col = "blue", lwd = 2)
qplot(sample = StudyB$Friend_index, stat = "qq")

t.test(StudyB$`YOU`, StudyB$`OTHER`, paired = TRUE)
t.test(StudyB$`OriginalSplitOther`, StudyB$`OTHER`, paired = TRUE)
describe(StudyB)

money_split <- glm(MagnitudeSplit ~ FAC4_4 + Congruence2 + ConditionP2 + GuiltJudged_P2P1bi:FriendBI, dat= StudyB)
summary(money_split)

summary(glm(GuiltJudged_P2P1 ~ Friend_index + ConditionP2 + Congruence2, dat=StudyBall))

summary(glm(MagnitudeSplit ~ GuiltJudged_P2P1bi:FriendBI, data=StudyB))

#summary(glm(MagnitudeSplit ~ FAC1_1 + FAC2_1 + FAC3_1, data=StudyB))
#summary(glm(GuiltJudged_P2P1 ~ FAC1_1 + FAC2_1 + FAC3_1, dat=StudyB))
summary(glm(GuiltJudged_P2P1 ~ Friend_index, dat=StudyB))
summary(glm(GuiltJudged_P2P1 ~ ConditionP2, dat=StudyB))

StudyB2 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB.xlsx", sheet = "figure")
labels<- c("Weak friendship", "Strong friendship")
names(labels)<-c("0Weak", "Strong")
Sum <- summarySE(StudyB2, measurevar="MagnitudeSplit", groupvars=c("FriendBI","GuiltJudged_P2P1bi"))
ggplot(Sum, aes(x=factor(GuiltJudged_P2P1bi), y=MagnitudeSplit, group = 1))  + 
  geom_point(size = 0.7) + 
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_smooth(method=lm, se=T, fullrange=F, size=.3, colour="black", linetype="dashed") + 
  geom_errorbar(aes(ymin=MagnitudeSplit-se, ymax=MagnitudeSplit+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white')) + 
  facet_grid(.~ FriendBI, labeller = labeller(FriendBI = labels)) + 
  scale_y_continuous(limits=c(0, 0.3)) + 
  xlab("Guilt Judged") + ylab("Change in reward split") +
  geom_hline(yintercept = .13333, linetype-"dotted", alpha=0.7, color ="red") +
  geom_hline(yintercept = .2, linetype="dashed", alpha=0.7, color = "red") +
  geom_hline(yintercept = .16666, linetype="dotted", alpha=0.7)


StudyB2 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB.xlsx", sheet = "figure")
labels<- c("Weak friendship", "Strong friendship")
names(labels)<-c("0Weak", "Strong")
Sum <- summarySE(StudyB2, measurevar="ChangeSplit", groupvars=c("FriendBI","GuiltJudged_P2P1bi"))
ggplot(Sum, aes(x=factor(GuiltJudged_P2P1bi), y=ChangeSplit, group = 1))  + 
  geom_point(size = 0.7) + 
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_smooth(method=lm, se=T, fullrange=F, size=.3, colour="black", linetype="dashed") + 
  geom_errorbar(aes(ymin=ChangeSplit-se, ymax=ChangeSplit+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white')) + 
  facet_grid(.~ FriendBI, labeller = labeller(FriendBI = labels)) + 
  xlab("Judgement of Guilt in partner") + ylab("Number of coins moved")

Sum <- summarySE(StudyB2, measurevar="You_Other", groupvars=c("FriendBI","GuiltJudged_P2P1bi"))
ggplot(Sum, aes(x=factor(GuiltJudged_P2P1bi), y=You_Other, group = 1))  + 
  geom_point(size = 0.7) + 
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_smooth(method=lm, se=T, fullrange=F, size=.3, colour="black", linetype="dashed") + 
  geom_errorbar(aes(ymin=You_Other-se, ymax=You_Other+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white')) + 
  facet_grid(.~ FriendBI, labeller = labeller(FriendBI = labels)) + 
  xlab("Judgement of Guilt in partner") + ylab("Number of extra coins kept for self")

Sum <- summarySE(StudyB, measurevar="MagnitudeSplit", groupvars=c("FriendBI","ConditionP2"))
ggplot(Sum, aes(x=factor(FriendBI), y=MagnitudeSplit, group = 1))  + 
  geom_point(size = 0.7) + geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_smooth(method=lm, se=T, fullrange=F, size=.3, colour="black", linetype="dashed") +
  geom_errorbar(aes(ymin=MagnitudeSplit-se, ymax=MagnitudeSplit+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white')) +
  facet_grid(.~ ConditionP2) + 
  scale_y_continuous(limits=c(0, 0.3)) + geom_hline(yintercept = .16666, linetype="dotted", alpha=0.7)

Sum <- summarySE(StudyB2, measurevar="You_Other", groupvars=c("GuiltJudged_P2P1bi"))
ggplot(Sum, aes(x=factor(GuiltJudged_P2P1bi), y=You_Other, group = 1))  + 
  geom_point(size = 0.7) + geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_smooth(method=lm, se=T, fullrange=F, size=.3, colour="black", linetype="dashed") +
  geom_errorbar(aes(ymin=You_Other-se, ymax=You_Other+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white'))

#--------------Correlation Felt - Judged Guilt-------------#
cor.test(StudyB$GuiltSelf_P1,StudyB$GuiltJudged_P2P1) #ouch...

StudyB2 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB.xlsx", sheet = "Condition1")
cor.test(StudyB2$Guilt3_P1,StudyB2$GuiltJudged_P2P1)

StudyB3 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB.xlsx", sheet = "Condition2")
cor.test(StudyB3$Guilt3_P1,StudyB3$GuiltJudged_P2P1)

#--------------Impact of GASP on Money Split-------------#
summary(lm(MagnitudeSplit ~ G_R, dat= StudyB))
summary(lm(MagnitudeSplit ~ G_NBE, dat= StudyB))
summary(lm(MagnitudeSplit ~ S_W, dat= StudyB))
summary(lm(MagnitudeSplit ~ S_NSE, dat= StudyB))

summary(lm(Money_split ~ G_R, dat= StudyBall))
summary(lm(Money_split ~ G_NBE, dat= StudyBall))
summary(lm(Money_split ~ S_W, dat= StudyBall))
summary(lm(Money_split ~ S_NSE, dat= StudyBall))

#--------------Impact of GASP on Judged Guilt-------------#
summary(lm(GuiltJudged_P2P1 ~ G_R,  dat=StudyBall))
summary(lm(GuiltJudged_P2P1 ~ G_NBE, dat=StudyBall))
summary(lm(GuiltJudged_P2P1 ~ S_W, dat=StudyBall))
summary(lm(GuiltJudged_P2P1 ~ S_NSE, dat=StudyBall))

###----------------- Personality Data ---------------------####
PersonalityBall <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB_Personality.xlsx", sheet = "All")
View(PersonalityBall)

#Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism

Money <- glm(Money_split ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityBall)
summary(Money)
summary(glm(Money_split ~ WhoseFault, dat=PersonalityBall))
summary(glm(Money_split ~ AnnoyancePartner, dat=PersonalityBall))
summary(glm(Money_split ~ AnnoyanceSelf, dat=PersonalityBall))

Guilt <- glm(GuiltJudged_P2P1 ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityBall)
summary(Guilt)
summary(glm(GuiltJudged_P2P1 ~ AnnoyanceSelf + AnnoyancePartner, dat = PersonalityBall))

summary(lm(Guilt3 ~ G_NBE, dat = PersonalityBall))
summary(lm(Guilt3 ~ G_R, dat = PersonalityBall))
summary(lm(Guilt3 ~ S_NSE, dat = PersonalityBall))
summary(lm(Guilt3 ~ S_W, dat = PersonalityBall))

cor.test(PersonalityBall$G_NBE,PersonalityBall$S_NSE)

GuiltP2 <- glm(Guilt3 ~Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityBall)
summary(GuiltP2)
summary(glm(Guilt3 ~ AnnoyanceSelf + AnnoyancePartner, dat = PersonalityBall))


Annoyance <- glm(AnnoyanceSelf ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityBall)
summary(Annoyance)
AnnoyanceP <- glm(AnnoyancePartner ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityBall)
summary(AnnoyanceP)


####---------------------STUDY B- CONDITION 3--------------------------####
###---------------- Affect -------------------------------###
AffectB3 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB3_Affect.xlsx")
View(AffectB3)

t.test(AffectB3$`1PositiveAffect`,AffectB3$`3PositiveAffectAfter`, paired = TRUE) # Pos Panas 1 vs 3
t.test(AffectB3$`1NegativeAffect`,AffectB3$`3NegativeAffectAfter`, paired = TRUE) # Neg Panas 1 vs 3
t.test(AffectB3$`1Guilty`,AffectB3$`3Guilty`, paired = TRUE) # Guilty Panas 1 vs 3

t.test(AffectB3$`2PositiveAffect`,AffectB3$`3PositiveAffectAfter`, paired = TRUE) # Pos Panas 2 vs 3
t.test(AffectB3$`2NegativeAffect`,AffectB3$`3NegativeAffectAfter`, paired = TRUE) # Neg Panas 2 vs 3
t.test(AffectB3$`2Guilty`,AffectB3$`3Guilty`, paired = TRUE) # Guilty Panas 2 vs 3
t.test(AffectB3$`2Ashamed`,AffectB3$`3Ashamed`, paired = TRUE)
t.test(AffectB3$`2Distressed`,AffectB3$`3Distressed`, paired = TRUE)
t.test(AffectB3$`2Proud`,AffectB3$`3Proud`, paired = TRUE)

t.test(AffectB3$`3Ashamed`, AffectB3$`3Guilty`, paired = TRUE)

# Affect plot
Affect_Plot_13 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB3_Affect_Plot_13.xlsx")
Sum <- summarySE(Affect_Plot_13, measurevar="Score", groupvars=c("Affect")) # PANAS 1 vs 3
ggplot(Sum, aes(x=factor(Affect), y=Score, group = 1))  + geom_point(size = 0.7) + geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + geom_errorbar(aes(ymin=Score-se, ymax=Score+se), size = .6, width=.2, alpha=0.8) + theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white'))

Affect_Plot_23 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB3_Affect_Plot_23.xlsx")
Sum <- summarySE(Affect_Plot_23, measurevar="Score", groupvars=c("Affect")) # PANAS 2 vs 3
ggplot(Sum, aes(x=factor(Affect), y=Score, group = 1))  + geom_point(size = 0.7) + 
  xlab("Self-reported affect") + ylab("Change before/after feedback") +
  geom_bar(stat="identity", fill="lightgrey", colour="black", width=0.9) + 
  geom_errorbar(aes(ymin=Score-se, ymax=Score+se), size = .6, width=.2, alpha=0.8) + 
  theme(strip.background = element_blank(), panel.border = element_rect(fill=NA, colour="black", size=0.7), panel.background = element_rect(fill = 'white'))


###----------------- GENERAL RESULTS StudyB3 -----------------------####
StudyB3 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB3.xlsx")
View(StudyB3)

t.test(StudyB3$`YOU`, StudyB3$`OTHER`, paired = TRUE)
describe(StudyB3)

money_split <- glm(Percent_Other ~ GuiltJudged_P2P1 + Friend_index + Congruence2, dat= StudyB3)
summary(money_split)

summary(lm(P1_Guilty ~ Friend_index, dat=StudyB3))

###----------------- Personality Data ---------------------####
PersonalityB3 <- read_excel("G:/1_STUDY3-PairStudy/DATA_Study3-PairStudy/Data&Analysis/1_StudyB/Data_StudyB3_Personality.xlsx")
View(PersonalityB3)

Money <- glm(Money_split ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat= PersonalityB3)
summary(Money)
summary(lm(Money_split ~ G_R, dat= PersonalityB3))
summary(lm(Money_split ~ G_NBE, dat= PersonalityB3))
summary(lm(Money_split ~ S_W, dat= PersonalityB3))
summary(lm(Money_split ~ S_NSE, dat= PersonalityB3))

Guilt <- glm(Guilt_P1 ~ Extraversion + Agreebleness + Conscientiousness + Emotional_Stability + Openness + Machiavelism + Psychopathy + Narcissism, dat = PersonalityB3)
summary(Guilt)
summary(lm(Guilt_P1 ~ G_R, dat= PersonalityB3))
summary(lm(Guilt_P1 ~ G_NBE, dat= PersonalityB3))
summary(lm(Guilt_P1 ~ S_W, dat= PersonalityB3))
summary(lm(Guilt_P1 ~ S_NSE, dat= PersonalityB3))
