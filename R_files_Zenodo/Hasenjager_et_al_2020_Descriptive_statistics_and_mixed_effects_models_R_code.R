###############################################################################
#Code for descriptive statistics and mixed effects models presented in:
#Hasenjager, M. J., et al. (2020). Proc. R. Soc. B. doi: 10.1098/rspb.2020.1871
###############################################################################

#Set working directory to where the data files are located on your computer
setwd()

#Load data sets
#Individual-level data
indivData<-read.csv("Hasenjager et al 2020_Guppy individual-level data.csv", header=TRUE)

#Group-level data
groupData<-read.csv("Hasenjager et al 2020_Guppy group-level data.csv", header=TRUE)

#To run analyses without group B7, which only contained 9 group members
#indivData<-subset(indivData, GroupID != "B7")
#groupData<-subset(groupData, GroupID != "B7")

#Load packages
library(lme4)
library(emmeans)
library(ggplot2)
library(ggbeeswarm)
library(DHARMa)
library(glmmTMB)
library(MuMIn)

#######################
#Descriptive statistics
#######################

#Total number of solvers
sum(indivData$Solved)

#Total number of individuals across all groups
nrow(indivData)

#Bold-dominated groups:
#Number of solvers
sum(indivData[indivData$GroupComposition=="Bold", ]$Solved)
#Mean number of solvers per group
mean(groupData[groupData$GroupComposition=="Bold", ]$NumberSolves)
#Minimum and maximum number of solvers per group
min(groupData[groupData$GroupComposition=="Bold", ]$NumberSolves)
max(groupData[groupData$GroupComposition=="Bold", ]$NumberSolves)

#Mixed groups:
#Number of solvers
sum(indivData[indivData$GroupComposition=="Mix", ]$Solved)
#Mean number of solvers per group
mean(groupData[groupData$GroupComposition=="Mix", ]$NumberSolves)
#Minimum and maximum number of solvers per group
min(groupData[groupData$GroupComposition=="Mix", ]$NumberSolves)
max(groupData[groupData$GroupComposition=="Mix", ]$NumberSolves)

#Shy-dominated groups:
#Number of solvers
sum(indivData[indivData$GroupComposition=="Shy", ]$Solved)
#Mean number of solvers per group
mean(groupData[groupData$GroupComposition=="Shy", ]$NumberSolves)
#Minimum and maximum number of solvers per group
min(groupData[groupData$GroupComposition=="Shy", ]$NumberSolves)
max(groupData[groupData$GroupComposition=="Shy", ]$NumberSolves)

#Number of groups in which the first individual to solve the task was bold
#Bold-dominated groups (out of 11 groups)
sum(lengths(regmatches(groupData[groupData$GroupComposition=="Bold", ]$FirstSolverPers, 
	gregexpr("Bold", groupData[groupData$GroupComposition=="Bold", ]$FirstSolverPers))))

#Mixed groups (out of 11 groups)
sum(lengths(regmatches(groupData[groupData$GroupComposition=="Mix", ]$FirstSolverPers, 
	gregexpr("Bold", groupData[groupData$GroupComposition=="Mix", ]$FirstSolverPers))))

#Shy-dominated groups (out of 12 groups)
sum(lengths(regmatches(groupData[groupData$GroupComposition=="Shy", ]$FirstSolverPers, 
	gregexpr("Bold", groupData[groupData$GroupComposition=="Shy", ]$FirstSolverPers))))

#Median and interquartile range for the duration (sec) of gaps between consecutive solving events within a group
#Bold-dominated groups
quantile(indivData[indivData$GroupComposition=="Bold", ]$LatencySincePreviousSolve, na.rm=TRUE)

#Mixed groups
quantile(indivData[indivData$GroupComposition=="Mix", ]$LatencySincePreviousSolve, na.rm=TRUE)

#Shy-dominated groups
quantile(indivData[indivData$GroupComposition=="Shy", ]$LatencySincePreviousSolve, na.rm=TRUE)


#Number of individuals for which at least 1 feeding strike was observed
#Bold individuals
nrow(indivData[which(indivData$Personality=="Bold" & indivData$Fed=="1"), ])

#Shy individuals
nrow(indivData[which(indivData$Personality=="Shy" & indivData$Fed=="1"), ])


#Percentage of individuals that fed during their first entry into the device
(sum(indivData$FirstFeedingVisit==1, na.rm=TRUE)/sum(indivData$Fed, na.rm=TRUE))*100


#Mean (SD) latency to begin feeding upon solving the task
#Bold individuals
mean(indivData[indivData$Personality=="Bold", ]$FeedingLatency, na.rm=TRUE)
sd(indivData[indivData$Personality=="Bold", ]$FeedingLatency, na.rm=TRUE)

#Shy individuals
mean(indivData[indivData$Personality=="Shy", ]$FeedingLatency, na.rm=TRUE)
sd(indivData[indivData$Personality=="Shy", ]$FeedingLatency, na.rm=TRUE)

##################################
#To reproduce Figure 2 (main text)
##################################

ggplot(data = groupData, aes(y = NumberSolves/GroupSize, x = GroupComposition, colour = GroupComposition, fill = GroupComposition)) + 
	ylab("Proportion of solvers per group") + xlab("Group personality composition") +
	geom_boxplot(show.legend = FALSE, lwd = 0.8, width=0.5) +
	geom_beeswarm(show.legend = FALSE, size = 3, cex = 2) +
	scale_colour_manual(values=c("#CC79A7", "#E69F00", "#56B4E9")) + 
	scale_fill_manual(values=c("thistle", "lightgoldenrod", "lightcyan")) +
	scale_x_discrete(labels=c("Bold-dominated", "Mixed", "Shy-dominated")) +  
	scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),limits=c(0,1))+
	theme(axis.text=element_text(size=18),axis.title=element_text(size=18)) +
	theme(axis.title.x=element_text(margin=margin(t=12)), axis.title.y=element_text(margin=margin(r=12))) +
	theme(panel.background=element_rect(fill="white"), axis.line.x.bottom=element_line(size=1), axis.line.y.left=element_line(size=1))


#########################
#Mean elective group size
#########################

#Standardise the categorical input variables
mixGPC<-stdize(indivData$MixComposition)
shyGPC<-stdize(indivData$ShyComposition)
Pers<-stdize(indivData$Personality)

#Fit the global model
mEGS1<-lmer(MeanElectiveGroupSize ~ mixGPC * Pers + shyGPC * Pers + BodyLengthCentred + (1|GroupID), REML=TRUE, data=indivData)

#Inspection of model residuals reveals no issues
R1<-resid(mEGS1)
qqnorm(R1)
qqline(R1)
shapiro.test(R1)
plot(R1~indivData$GroupComposition)
plot(R1~indivData$Personality)
plot(R1~indivData$BodyLengthCentred)
F1<-fitted(mEGS1)
plot(R1~F1)

#Get parameter estimates (presented in Table S3)
summary(mEGS1)

#Likelihood ratio tests
#Models are first refitted with maximum likelihood estimation

#Interaction term between group composition and individual personality
mEGS1.ML<-lmer(MeanElectiveGroupSize ~ mixGPC * Pers + shyGPC * Pers + BodyLengthCentred + (1|GroupID), REML=FALSE, data=indivData)
mEGS2.ML<-lmer(MeanElectiveGroupSize ~ mixGPC + shyGPC + Pers + BodyLengthCentred + (1|GroupID), REML=FALSE, data=indivData)
anova(mEGS1.ML, mEGS2.ML)

#Main effect of group composition
mEGS3.ML<-lmer(MeanElectiveGroupSize ~ Pers + BodyLengthCentred +(1|GroupID), REML=FALSE, data=indivData)
anova(mEGS2.ML, mEGS3.ML)

#Main effect of personality
mEGS4.ML<-lmer(MeanElectiveGroupSize ~ mixGPC + shyGPC + BodyLengthCentred +(1|GroupID), REML=FALSE, data=indivData)
anova(mEGS2.ML, mEGS4.ML)

#Main effect of body length
mEGS5.ML<-lmer(MeanElectiveGroupSize ~ mixGPC * Pers + shyGPC * Pers +(1|GroupID), REML=FALSE, data=indivData)
anova(mEGS1.ML, mEGS5.ML)

#Estimated marginal means for mean elective group size for bold and shy individuals
#Fit model without interaction term using REML estimation
mEGS2<-lmer(MeanElectiveGroupSize ~ mixGPC + shyGPC + Pers + BodyLengthCentred + (1|GroupID), REML=TRUE, data=indivData)

#Get EMMs
emmeans(mEGS2, pairwise~Pers)

#########################################
#Probability of solving the foraging task
#########################################

#Standardise the categorical input variables
mixGPC<-stdize(indivData$MixComposition)
shyGPC<-stdize(indivData$ShyComposition)
Pers<-stdize(indivData$Personality)

S1<-glmer(Solved ~ mixGPC * Pers + shyGPC * Pers + BodyLengthCentred +(1|GroupID), family=binomial, data=indivData)

#Model validation carried out with the DHARMa package using simulated residuals
resSim<-simulateResiduals(S1, n=10000)
plot(resSim)
testQuantiles(resSim)
testUniformity(resSim)
testDispersion(resSim)
testZeroInflation(resSim)
plotResiduals(resSim, indivData$GroupComposition)
plotResiduals(resSim, indivData$Personality)
plotResiduals(resSim, indivData$BodyLengthCentred)

#Get parameter estimates (presented in Table S4)
summary(S1)

#Likelihood ratio tests

#Interaction term between group composition and individual personality
S2<-glmer(Solved ~ mixGPC + shyGPC + Pers + BodyLengthCentred +(1|GroupID), family=binomial, data=indivData)
anova(S1, S2)

#Main effect of group composition
S3<-glmer(Solved ~ Pers + BodyLengthCentred +(1|GroupID), family=binomial, data=indivData)
anova(S2, S3)

#Main effect of personality
S4<-glmer(Solved ~ mixGPC + shyGPC + BodyLengthCentred +(1|GroupID), family=binomial, data=indivData)
anova(S2, S4)

#Main effect of body length
S5<-glmer(Solved ~ mixGPC * Pers + shyGPC * Pers +(1|GroupID), family=binomial, data=indivData)
anova(S1, S5)

#Estimated marginal means for solving probability by group composition
S6<-glmer(Solved ~ GroupComposition + Personality + BodyLengthCentred +(1|GroupID), family=binomial, data=indivData)
emmeans(S6, pairwise~GroupComposition, type="response")

###################################
#Entry rate of informed individuals
###################################

#Restrict data to only those individuals that solved the task
indivData_Solvers<-indivData[indivData$Solved==1, ]

#Create variable indicating the duration of time each individual was informed during the trial
informedTime<-(1200-indivData_Solvers$SolveLatency)/60

#Mean (SD) entries into the device per minute by informed individuals
mean(indivData_Solvers$Entries/informedTime)
sd(indivData_Solvers$Entries/informedTime)

#Standardise the categorical input variables
mixGPC_S<-stdize(indivData_Solvers$MixComposition)
shyGPC_S<-stdize(indivData_Solvers$ShyComposition)
Pers_S<-stdize(indivData_Solvers$Personality)

#Fit global count model using a truncated Poisson distribution (by definition, informed individuals entred the device at least once)
ER1<-glmmTMB(Entries ~ mixGPC_S * Pers_S + shyGPC_S * Pers_S + BodyLengthCentred + offset(log(informedTime)) + (1|GroupID), family = truncated_poisson, data=indivData_Solvers)

#Model validation carried out with the DHARMa package using simulated residuals
resSim<-simulateResiduals(ER1, n=10000)
plot(resSim)
testQuantiles(resSim)
testUniformity(resSim)
testDispersion(resSim)
plotResiduals(resSim, indivData_Solvers$GroupComposition)
plotResiduals(resSim, indivData_Solvers$Personality)
plotResiduals(resSim, indivData_Solvers$BodyLengthCentred)

#Get parameter estimates (presented in Table S5)
summary(ER1)

#Likelihood ratio tests

#Interaction term between group composition and individual personality
ER2<-glmmTMB(Entries ~ mixGPC_S + shyGPC_S + Pers_S + BodyLengthCentred + offset(log(informedTime)) + (1|GroupID), family = truncated_poisson, data=indivData_Solvers)
anova(ER1, ER2)

#Main effect of group composition
ER3<-glmmTMB(Entries ~ Pers_S + BodyLengthCentred + offset(log(informedTime)) + (1|GroupID), family = truncated_poisson, data=indivData_Solvers)
anova(ER2, ER3)

#Main effect of personality
ER4<-glmmTMB(Entries ~ mixGPC_S + shyGPC_S + BodyLengthCentred + offset(log(informedTime)) + (1|GroupID), family = truncated_poisson, data=indivData_Solvers)
anova(ER2, ER4)

#Main effect of body length
ER5<-glmmTMB(Entries ~ mixGPC_S * Pers_S + shyGPC_S * Pers_S + offset(log(informedTime)) + (1|GroupID), family = truncated_poisson, data=indivData_Solvers)
anova(ER1, ER5)

#######################
#To reproduce Figure S4
#######################

#Remove individual with an estimated entry rate of 18.2 times per minute due to solving the task 3 sec before the end of the trial
S4plotData<-indivData_Solvers[indivData_Solvers$Entries/informedTime < 3, ]

#Create variable indicating the duration of time each individual was informed during the trial
S4plotInformedTime<-(1200-S4plotData$SolveLatency)/60

#Fix order of groups from earliest to latest within group compositions
S4plotData$GroupID <- factor(S4plotData$GroupID, levels = c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B12", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M11", "M12", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12"))

ggplot(data = S4plotData, aes(y = S4plotData$Entries/S4plotInformedTime, x = GroupID, colour = GroupComposition, fill = GroupComposition)) + 
	ylab("Entries per min") + xlab("Group ID") +
	geom_boxplot(show.legend = FALSE, lwd = 0.8) +
	scale_colour_manual(values=c("#CC79A7", "#E69F00", "#56B4E9")) + 
	scale_fill_manual(values=c("thistle", "lightgoldenrod", "lightcyan")) +
	scale_x_discrete(labels=c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B12", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M11", "M12", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12")) +  
	theme(axis.text=element_text(size=12, angle=90),axis.title=element_text(size=12))

