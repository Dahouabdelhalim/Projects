#Changelog

#7-8-2017   File created
#7-10-2017 added ANOVA
#7-11-2017 added chi-square
#7-12-2017 replaced ANOVA with GLMM analysis
#8-26-2017 completed microbe pref data set
#9-3-2017 added prelim supernatant pref data set
#9-8-2017 ,PRE supernatant pref data set
#9-19-2017, complete supernatant pref data set
#11-27-2017, added init pref vs test pref for diff condit data set
#12-8-2017, added full diff cond data set, dharma and effects packages
#12-16-2017, added frankenflower data. All data sets now included.
#1-4-2018, added number of visits to reach learn criterion. For diff cond exp



#LOAD EVERYTHING AFTER THIS

#Name of data file:
FILENAME1 <- "GLMER_Preference_Control-Microbe_Exacum.csv"
FILENAME2 <- "ANOVA_Percent-Rejection_Exacum.csv"
FILENAME3 <- "TTest_Preference_Control-Yeast-Bacteria-Super.csv"
FILENAME4 <- "GLMER_ControlPropAppReject_Control-Microbe_Exacum.csv"
FILENAME5 <- "GLMER_MicrobePropAppReject_Control-Microbe_Exacum.csv"
FILENAME6 <- "LM_WildBees_MicrobeCount-vs-HeadWidth.csv"
FILENAME7 <- "WilcoxTest_WildBees_SizeVsPollenCollection.csv"
FILENAME8 <- "LM_UrbanFlowers_MicrobeCount-vs-Width.csv"
#FILENAME9 <- "GLMER_ControlPropLandReject_Control-Microbe_Exacum.csv" #Not useful
#FILENAME10 <- "GLMER_MicrobePropLandReject_Control-Microbe_Exacum.csv" #Not useful
FILENAME11 <- "GLMER_Preference_Control-Super_Exacum.csv"
FILENAME12 <- "GLMER_ControlPropAppReject_Control-Super_Exacum.csv"
FILENAME13 <- "GLMER_SuperPropAppReject_Control-Super_Exacum.csv"
FILENAME14 <- "TTest_Preference_ArtifFlwr-Microbe.csv"
FILENAME15 <- "GLMER_ArtiFlwr-FullPreference.csv"
FILENAME16 <- "GLMER_ArtiFlwr_Microbe.csv"
FILENAME17 <- "GLMER_ArtiFlwr_Control.csv"
FILENAME18 <- "TTest_Preference_Frankenflower.csv"
FILENAME19 <- "TTest_NoVisits-Training_ArtiFlwer.csv"

#Set the directory your data are in:
DATADIR <- "C:\\\\Users\\\\Avery\\\\Dropbox\\\\ALR PEEP Experiments\\\\P1 - Floral Microbiome Bee Preferences\\\\#Stats" #Lenovo computer


#Load data:
setwd(DATADIR); dataset1 <- read.table(FILENAME1, sep=",", header=TRUE)
setwd(DATADIR); dataset2 <- read.table(FILENAME2, sep=",", header=TRUE)
setwd(DATADIR); dataset3 <- read.table(FILENAME3, sep=",", header=TRUE)
setwd(DATADIR); dataset4 <- read.table(FILENAME4, sep=",", header=TRUE)
setwd(DATADIR); dataset5 <- read.table(FILENAME5, sep=",", header=TRUE)
setwd(DATADIR); dataset6 <- read.table(FILENAME6, sep=",", header=TRUE)
setwd(DATADIR); dataset7 <- read.table(FILENAME7, sep=",", header=TRUE)
setwd(DATADIR); dataset8 <- read.table(FILENAME8, sep=",", header=TRUE)
#setwd(DATADIR); dataset9 <- read.table(FILENAME9, sep=",", header=TRUE) #Not useful
#setwd(DATADIR); dataset10 <- read.table(FILENAME10, sep=",", header=TRUE) #Not useful
setwd(DATADIR); dataset11 <- read.table(FILENAME11, sep=",", header=TRUE)
setwd(DATADIR); dataset12 <- read.table(FILENAME12, sep=",", header=TRUE)
setwd(DATADIR); dataset13 <- read.table(FILENAME13, sep=",", header=TRUE)
setwd(DATADIR); dataset14 <- read.table(FILENAME14, sep=",", header=TRUE)
setwd(DATADIR); dataset15 <- read.table(FILENAME15, sep=",", header=TRUE)
setwd(DATADIR); dataset16 <- read.table(FILENAME16, sep=",", header=TRUE)
setwd(DATADIR); dataset17 <- read.table(FILENAME17, sep=",", header=TRUE)
setwd(DATADIR); dataset18 <- read.table(FILENAME18, sep=",", header=TRUE)
setwd(DATADIR); dataset19 <- read.table(FILENAME19, sep=",", header=TRUE)

#install.packages("XNomial")

#Load the libraries:
library("mgcv")
library("lme4")
library("lsmeans")
library("multcomp")
library("lmerTest")
library("car")
library("effects")
library("DHARMa")
library("ggplot2")
library("DescTools")



#
#G Test
#
#Determine whether first visit (innate) choice significantly different
#Differential Conditioning experiment
#

InitialPref = c(19,18) #Sums to the total N
GTest(InitialPref) #G = 0.02703, X-squared df = 1, p-value = 0.8694, N=37





#
#Paired T-tests to determine whether Overall Initial control vs microbe Preference sig different
#Differential Conditioning experiment
#

#Inputs:
IP_Control <- dataset14$InitialPref_overall1
IP_Microbe <- dataset14$InitialPref_overall2

#Compare Microbe vs Control, paired t-test
#t.test(DCM1,DCM2,paired=TRUE)
t.test(IP_Control,IP_Microbe,paired=TRUE) #t = -0.4221, df = 36, p-value = 0.6755, N=37
shapiro.test(IP_Control) #Normal
shapiro.test(IP_Microbe) #Normal
var.test(IP_Control,IP_Microbe) #Normal






#
#Paired T-tests to determine whether Initial vs Test Preference sig different
#Differential Conditioning experiment
#

#Inputs:
DCM1 <- dataset14$InitialPref_Microbe
DCM2 <- dataset14$TestPref_Microbe
DCC1 <- dataset14$InitialPref_Control
DCC2 <- dataset14$TestPref_Control

#Compare Initial vs Trained to Microbe, paired t-test
#t.test(DCM1,DCM2,paired=TRUE)
wilcox.test(DCM1,DCM2,paired=TRUE) #V = 77, p-value = 0.0009766, N=12
shapiro.test(DCM1) #Normal
shapiro.test(DCM2) #Not Normal
var.test(DCC1,DCM2) #Normal

#Compare Initial vs Trained to Control, paired t-test
#t.test(DCC1,DCC2,paired=TRUE)
t.test(DCC1,DCC2,paired=TRUE) #t = -6.3516, df = 11, p-value = 5.433e-05, N=12
shapiro.test(DCC1) #Normal
shapiro.test(DCC2) #Normal
var.test(DCC1,DCC2) #Normal



#
#T-tests to determine whether num visits to reach learning criterion differs with Control and Microbe tmt
#Differential Conditioning experiment
#

#Inputs:
CNV <- dataset19$Control_VisitsToLearn
MNV <- dataset19$Microbe_VisitsToLearn


#Compare number of visits to learning for Control vs Microbe, t-test
#t.test(DCM1,DCM2,paired=TRUE)
wilcox.test(CNV,MNV,paired=FALSE) #W = 0.90681, p-value = 0.1942, N=24
shapiro.test(CNV) #Normal
shapiro.test(MNV) #Not Normal
var.test(CNV,MNV) #Not Normal







#Inputs:
Correct <- dataset15$Correct #Did bee visit the rewarding flower
Treatment <- dataset15$Treatment #How many males in tmt
Colony <- dataset15$Colony
BeeID <- dataset15$BeeID
NumVisits <- dataset15$NumVisits


#Order of BeeID and Colony matters
#VisitNo treated as a repeated measure within BeeID within colony
DC_Correct <- glmer(Correct ~ Treatment * NumVisits + (1 | NumVisits / BeeID / Colony), family = binomial, data = dataset15)
summary(DC_Correct)

#Type 2 sums of squares, for Wald Test reporting. As Z(<p are just approximations; Not useful with multiple factors)
#Note the capital Anova.
Anova(DC_Correct)

#Post-hoc test to see significance of the pairs (Treatment 1 vs Treatment 2)
summary(glht(DC_Correct, linfct=mcp(Treatment="Tukey")))

#RESULTS

#Response: Correct
#Chisq Df Pr(>Chisq)    
#Treatment            0.5378  1   0.463328    
#NumVisits           20.8969  1  4.847e-06 ***
#Treatment:NumVisits  8.4497  1   0.003651 ** 



#
#Plot taking random effects into consideration
#with interaction between a factor and a continuous variable--just specify one predictor to just see a plot of that
#
plot(Effect(focal.predictors=c("Treatment","NumVisits"),DC_Correct)) #Treatment interaction with NumVisits

#Nicer plot, standardized visit number
effect <- Effect(focal.predictors=c("Treatment","NumVisits"),DC_Correct)
effect <- as.data.frame(effect)

#Plot model
Plot1 <- ggplot(effect, aes(NumVisits, fit, col=Treatment, linetype=Treatment)) + 
  geom_line(lwd=.5) + geom_ribbon(aes(ymin=lower, ymax=upper, fill=Treatment), alpha=0.3, colour=NA) + 
  scale_fill_manual(values = c('gray65','gray65'), guide=FALSE) + 
  scale_linetype_manual(values=c("solid","dashed")) + 
  scale_color_manual(values = c('royalblue3','deeppink3')) + 
  ylab("Proportion correct visit") + xlab("Standardized visit number") + theme_classic() + 
  theme(text=element_text(size=14)) + xlim(0,100) + ylim(0,1)

plot(Plot1) #Plot it


#
#Check GLMM model assumptions via DHARMa
#
DC_Correct.1 <- glmer(Correct ~ Treatment * NumVisits + (1 | NumVisits / BeeID / Colony), family = binomial, data = dataset15) 

simulationOutput1 <- simulateResiduals(fittedModel = DC_Correct.1, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks ok and the residual vs. predicted quantile lines are ok
testUniformity(simulationOutput = simulationOutput1) #GOOD: D = 0.020943, p-value = 0.55
testOverdispersionParametric(mod10) #Dispersion test doesnt run
testZeroInflation(simulationOutput1) #Zero inflation test yields ratioObsExp = 1.0011, p-value = 0.462; zero inflation is okay









#
#Plot logistic regression of all data. Need to split up
#Run both together to plot on same graph
#

#Microbe
fit1 = glm(correct ~ visit, data=dataset16, family=binomial)
newdat <- data.frame(visit=seq(min(dataset16$visit), max(dataset16$visit),len=100)) #Predictions
newdat$correct = predict(fit1, newdata=newdat, type="response")
plot(dataset16$visit, dataset16$correct, xlab="Visit Number", ylab="Proportion Choosing Reward", col="darkgrey")
lines(correct ~ visit, newdat, col="orange", lwd=2)

#Control
fit2 = glm(correct ~ visit, data=dataset17, family=binomial)
newdat <- data.frame(visit=seq(min(dataset17$visit), max(dataset17$visit),len=100))
newdat$correct = predict(fit2, newdata=newdat, type="response")

#To plot the two graphs separately
#plot(dataset17$visit, dataset17$correct, xlab="Visit Number", ylab="Proportion Choosing Reward", 
#     col="darkgrey", xlim=range(c(dataset16$visit, dataset17$visit))) 
lines(correct ~ visit, newdat, col="blue", lwd=2)

legend("bottomright", inset=.08, title="",
       c("Microbe","Control"), fill=terrain.colors(3), horiz=FALSE, bty = "n")




#
#Paired T-tests to determine whether Overall Initial control vs microbe Preference sig different
#Frankenflower experiment
#

#Inputs:
RAAC_Control <- dataset18$RAAC_Control
RAAC_Microbe <- dataset18$RAAC_Microbe

#Compare Microbe vs Control, paired t-test
#t.test(DCM1,DCM2,paired=TRUE)
t.test(RAAC_Control,RAAC_Microbe,paired=TRUE) #t = 2.4551, df = 13, p-value = 0.02893, N=14
shapiro.test(RAAC_Control) #Normal
shapiro.test(RAAC_Microbe) #Normal
var.test(RAAC_Control,RAAC_Microbe) #Normal

#Inputs:
AARC_Control <- dataset18$AARC_Control
AARC_Microbe <- dataset18$AARC_Microbe

#Compare Microbe vs Control, paired t-test
#t.test(DCM1,DCM2,paired=TRUE)
wilcox.test(AARC_Control,AARC_Microbe,paired=TRUE) #V = 56, p-value = 0.4846, N=15
shapiro.test(AARC_Control) #Not Normal
shapiro.test(AARC_Microbe) #Not Normal
var.test(AARC_Control,AARC_Microbe) #Normal


#
#G Test
#
#Determine whether first visit (innate) choice significantly different
#Frankenflower experiment
#

InitialPref_RAAC = c(10,4) #Sums to the total N
GTest(InitialPref_RAAC) #G = 2.6566, X-squared df = 1, p-value = 0.1031

InitialPref_AARC = c(9,6) #Sums to the total N
GTest(InitialPref_AARC) #G = 0.60407, X-squared df = 1, p-value = 0.437





#
#Logistic regression via GLMM
#
#Test preference for preference for control flowers vs microbe flowers, for each MICROBE treatment
#This is for the learning experiment with Exacum
#

#Inputs:
TypeChosen <- dataset1$TypeChosen #What flower type the bee visited
Treatment <- dataset1$Treatment #What flower type bee has been trained on
BeeID <- dataset1$BeeID
VisitNo <- dataset1$VisitNo
Colony <- dataset1$Colony

#test <- glmer(Species ~ Experience + (1 | Experience:BeeID), family = binomial)  #nesting 4-26-2015
#Run Logistic Regression, BeeID is a factor, Species is 'binary'; only two species options
microbeBL <- glmer(TypeChosen ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial) #VisitNo treated as a repeated measure within BeeID within colony
summary(microbeBL)

#Type 2 sums of squares, for Wal Test reporting. As Z(<p are just approximations; Not useful with multiple factors)
#Note the capital Anova.
Anova(microbeBL)

#Treatment isn't a random factor, just BeeID; so no nesting!
#Replace 1 with the Treatment if afraid random effect, if bees affecting that; interaction

#Post-hoc test to see significance of the pairs
summary(glht(microbeBL, linfct=mcp(Treatment="Tukey")))


#
#Check GLMM model assumptions via DHARMa
#
microbeBL.1 <- glmer(TypeChosen ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial, data = dataset1) 

simulationOutput1 <- simulateResiduals(fittedModel = microbeBL.1, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks ok and the residual vs. predicted quantile lines are ok
testUniformity(simulationOutput = simulationOutput1) #GOOD: D = 0.022765, p-value = 0.5446
testOverdispersionParametric(mod10) #Dispersion test doesnt run
testZeroInflation(simulationOutput1) #Zero inflation test yields ratioObsExp = 1.0012, p-value = 0.454; zero inflation is okay






#Test proportion of APPROACH rejection for CONTROL flowers only, for each MICROBE treatment
#This is for the learning experiment with Exacum
#

#Inputs:
Behavior <- dataset4$Behavior #Whether the bee rejected or accepted the flower (landed vs approached)
Treatment <- dataset4$Treatment #What flower type bee has been trained on
BeeID <- dataset4$BeeID
VisitNo <- dataset4$VisitNo
Colony <- dataset4$Colony


#Run Logistic Regression, BeeID is a factor, Species is 'binary'; only two species options
CReject <- glmer(Behavior ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial) #VisitNo treated as a repeated measure within BeeID within colony
summary(CReject)

#Type 2 sums of squares, for Wal Test reporting. As Z(<p are just approximations; Not useful with multiple factors)
#Note the capital Anova.
Anova(CReject)

#Treatment isn't a random factor, just BeeID; so no nesting!
#Replace 1 with the Treatment if afraid random effect, if bees affecting that; interaction

#Post-hoc test to see significance of the pairs
summary(glht(CReject, linfct=mcp(Treatment="Tukey")))


#
#Check GLMM model assumptions via DHARMa
#
CReject.1 <- glmer(Behavior ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial, data = dataset4) 

simulationOutput1 <- simulateResiduals(fittedModel = CReject.1, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks ok and the residual vs. predicted quantile lines are ok
testUniformity(simulationOutput = simulationOutput1) #GOOD: D = 0.019737, p-value = 0.747
testOverdispersionParametric(mod10) #Dispersion test doesnt run
testZeroInflation(simulationOutput1) #Zero inflation test yields ratioObsExp = 1.0063, p-value = 0.392; zero inflation is okay



#Test proportion of APPROACH rejection for MICROBE flowers only, for each MICROBE treatment
#This is for the learning experiment with Exacum
#

#Inputs:
Behavior <- dataset5$Behavior #Whether the bee rejected or accepted the flower (landed vs approached)
Treatment <- dataset5$Treatment #What flower type bee has been trained on
BeeID <- dataset5$BeeID
VisitNo <- dataset5$VisitNo
Colony <- dataset5$Colony


#Run Logistic Regression, BeeID is a factor, Species is 'binary'; only two species options
BReject <- glmer(Behavior ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial) #VisitNo treated as a repeated measure within BeeID within colony
summary(BReject)

#Type 2 sums of squares, for Wal Test reporting. As Z(<p are just approximations; Not useful with multiple factors)
#Note the capital Anova.
Anova(BReject)

#Treatment isn't a random factor, just BeeID; so no nesting!
#Replace 1 with the Treatment if afraid random effect, if bees affecting that; interaction

#Post-hoc test to see significance of the pairs
summary(glht(BReject, linfct=mcp(Treatment="Tukey")))


#
#Check GLMM model assumptions via DHARMa
#
BReject.1 <- glmer(Behavior ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial, data = dataset5) 

simulationOutput1 <- simulateResiduals(fittedModel = BReject.1, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks ok and the residual vs. predicted quantile lines are ok
testUniformity(simulationOutput = simulationOutput1) #GOOD: D = 0.025462, p-value = 0.482
testOverdispersionParametric(mod10) #Dispersion test doesnt run
testZeroInflation(simulationOutput1) #Zero inflation test yields ratioObsExp = 1.0174, p-value = 0.327; zero inflation is okay





#
#Logistic regression via GLMM
#
#Test preference for preference for control flowers vs supernatant flowers, for each SUPERNATANT treatment
#This is for the learning experiment with Exacum
#

#Inputs:
TypeChosen <- dataset11$TypeChosen #What flower type the bee visited
Treatment <- dataset11$Treatment #What flower type bee has been trained on
BeeID <- dataset11$BeeID
VisitNo <- dataset11$VisitNo
Colony <- dataset11$Colony

#Run Logistic Regression, BeeID is a factor, Species is 'binary'; only two species options
superBL <- glmer(TypeChosen ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial) #VisitNo treated as a repeated measure within BeeID within colony
summary(superBL)

#Type 2 sums of squares, for Wal Test reporting. As Z(<p are just approximations; Not useful with multiple factors)
#Note the capital Anova.
Anova(superBL)

#Treatment isn't a random factor, just BeeID; so no nesting!
#Replace 1 with the Treatment if afraid random effect, if bees affecting that; interaction

#Post-hoc test to see significance of the pairs
summary(glht(superBL, linfct=mcp(Treatment="Tukey")))


#
#Check GLMM model assumptions via DHARMa
#
superBL.1 <- glmer(TypeChosen ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial, data = dataset11) 

simulationOutput1 <- simulateResiduals(fittedModel = superBL.1, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks ok and the residual vs. predicted quantile lines are ok
testUniformity(simulationOutput = simulationOutput1) #GOOD: D = 0.025462, p-value = 0.482
testOverdispersionParametric(mod10) #Dispersion test doesnt run
testZeroInflation(simulationOutput1) #Zero inflation test yields ratioObsExp = 1.0174, p-value = 0.327; zero inflation is okay





#Test proportion of APPROACH rejection for CONTROL flowers only, for each SUPERNATANT treatment
#This is for the learning experiment with Exacum
#

#Inputs:
Behavior <- dataset12$Behavior #Whether the bee rejected or accepted the flower (landed vs approached)
Treatment <- dataset12$Treatment #What flower type bee has been trained on
BeeID <- dataset12$BeeID
VisitNo <- dataset12$VisitNo
Colony <- dataset12$Colony

#test <- glmer(Species ~ Experience + (1 | Experience:BeeID), family = binomial)  #nesting 4-26-2015
#Run Logistic Regression, BeeID is a factor, Species is 'binary'; only two species options
CReject <- glmer(Behavior ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial) #VisitNo treated as a repeated measure within BeeID
summary(CReject)

#Type 2 sums of squares, for Wal Test reporting. As Z(<p are just approximations; Not useful with multiple factors)
#Note the capital Anova.
Anova(CReject)

#Treatment isn't a random factor, just BeeID; so no nesting!
#Replace 1 with the Treatment if afraid random effect, if bees affecting that; interaction

#Post-hoc test to see significance of the pairs
summary(glht(CReject, linfct=mcp(Treatment="Tukey")))



#
#Check GLMM model assumptions via DHARMa
#
CReject.1 <- glmer(Behavior ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial, data = dataset12) 

simulationOutput1 <- simulateResiduals(fittedModel = CReject.1, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks ok and the residual vs. predicted quantile lines are ok
testUniformity(simulationOutput = simulationOutput1) #GOOD: D = 0.024473, p-value = 0.3957
testOverdispersionParametric(mod10) #Dispersion test doesnt run
testZeroInflation(simulationOutput1) #Zero inflation test yields ratioObsExp = 1.0121, p-value = 0.325; zero inflation is okay



#Test proportion of APPROACH rejection for SUPERNATANT flowers only, for each SUPERNATANT treatment
#This is for the learning experiment with Exacum
#

#Inputs:
Behavior <- dataset13$Behavior #Whether the bee rejected or accepted the flower (landed vs approached)
Treatment <- dataset13$Treatment #What flower type bee has been trained on
BeeID <- dataset13$BeeID
VisitNo <- dataset13$VisitNo
Colony <- dataset13$Colony

#test <- glmer(Species ~ Experience + (1 | Experience:BeeID), family = binomial)  #nesting 4-26-2015
#Run Logistic Regression, BeeID is a factor, Species is 'binary'; only two species options
BReject <- glmer(Behavior ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial) #VisitNo treated as a repeated measure within BeeID within colony
summary(BReject)

#Type 2 sums of squares, for Wal Test reporting. As Z(<p are just approximations; Not useful with multiple factors)
#Note the capital Anova.
Anova(BReject)

#Treatment isn't a random factor, just BeeID; so no nesting!
#Replace 1 with the Treatment if afraid random effect, if bees affecting that; interaction

#Post-hoc test to see significance of the pairs
summary(glht(BReject, linfct=mcp(Treatment="Tukey")))


#
#Check GLMM model assumptions via DHARMa
#
BReject.1 <- glmer(Behavior ~ Treatment + (1 | VisitNo / BeeID / Colony), family = binomial, data = dataset13) 

simulationOutput1 <- simulateResiduals(fittedModel = BReject.1, n = 1000)
plotSimulatedResiduals(simulationOutput = simulationOutput1) #Plot of QQ plot looks ok and the residual vs. predicted quantile lines are ok
testUniformity(simulationOutput = simulationOutput1) #GOOD: D = 0.014777, p-value = 0.927
testOverdispersionParametric(mod10) #Dispersion test doesnt run
testZeroInflation(simulationOutput1) #Zero inflation test yields ratioObsExp = 1.0169, p-value = 0.291; zero inflation is okay




#
#Paired T-tests to determine whether Yeast or Bacteria vs Control sig different
#


#Inputs:
M1 <- dataset3$ControlMLandBuzz
M2 <- dataset3$MicrobeLandBuzz
Y1 <- dataset3$ControlYLandBuzz
Y2 <- dataset3$YeastLandBuzz
B1 <- dataset3$ControlBLandBuzz
B2 <- dataset3$BacteriaLandBuzz
PM1 <- dataset3$Pooled_ControlMLandBuzz
PM2 <- dataset3$Pooled_MicrobeLandBuzz


#Compare Microbe vs Control, paired t-test
t.test(M1,M2,paired=TRUE) #t = 3.2077, df = 16, p-value = 0.005488
shapiro.test(M1) #Normal
shapiro.test(M2) #Normal
var.test(M1,M2) #Normal

#Compare Microbe vs Control, paired t-test
#Pooled across both experiments
wilcox.test(PM1,PM2,paired=TRUE) #V = 473, p-value = 9.633e-05
shapiro.test(PM1) #Not Normal
shapiro.test(PM2) #Not Normal
var.test(PM1,PM2) #Normal



#Compare Bacteria vs Control, paired t-test
t.test(B1,B2,paired=TRUE) #t = 2.5032, df = 13, p-value = 0.02643
shapiro.test(B1) #Normal
shapiro.test(B2) #Normal
var.test(B1,B2) #Normal


#Compare Yeast vs Control, paired t-test
t.test(Y1,Y2,paired=TRUE) #t = 0.16779, df = 13, p-value = 0.8693
shapiro.test(Y1) #Normal
shapiro.test(Y2) #Normal
var.test(Y1,Y2) #Normal




#
#G Test
#
#Determine whether first visit (innate) choice significantly different
#Bacteria vs Yeast vs Naive (not supernatant)

Naive = c(26,10) #Sums to the total N. 4+6=M. 28=C
GTest(Naive) #G = 7.366, X-squared df = 1, p-value = 0.006647, N = 36

Bacteria = c(12,2) #Sums to the total N
GTest(Bacteria) #G = 7.9249, X-squared df = 1, p-value = 0.004876

Yeast = c(8,6) #Sums to the total N
GTest(Yeast) #G = 0.28669, X-squared df = 1, p-value = 0.5923





#
#Linear Model
#
#Wild caught Bombus body size versus haemocytometer microbe cell count
#


#Inputs:
MicrobeCount <- dataset6$MicrobeCount #Number of microbes on bee
HeadWidth <- dataset6$HeadWidth #widest part of head in mm
CollectedPollen <- dataset6$CollectedPollen
#YeastPresent <- dataset6$YeastLevel

#Model
CountByBeeSize.m <-lm(log(MicrobeCount) ~ log(HeadWidth) * CollectedPollen, data=dataset6) #Add in pollen collection presence/absence

#Summary stats
summary(CountByBeeSize.m)

#See if the residuals are distributed normally for each group. 
#Must be done when a continuous distribution (like with lmer, anova etc)
CountByBeeSize.resid = resid(CountByBeeSize.m)

#Plot the residuals
hist(CountByBeeSize.resid)

#qq plot to see if residuals normally distributed
qqPlot(CountByBeeSize.m, main="QQ Plot")

#Check to see the effects of each factor
#Anova(specialvsbouts.m, type = "3") #To look specifically at the interaction
Anova(CountByBeeSize.m, type = "2")  #Multiple R-squared:  0.1221,	Adjusted R-squared:  0.1034 , F-statistic: 6.536 on 1 and 47 DF,  p-value: 0.01386
  #log(HeadWidth)                  0.738  1  0.9415    0.3371    
  #CollectedPollen                44.986  1 57.3545 1.451e-09 ***
  #log(HeadWidth):CollectedPollen  0.107  1  0.1363    0.7137   









#
#Wilcoxon signed-rank test
#
#Are wild-caught nectar foragers are different in size than pollen foragers
#


#Inputs:
P1 <- dataset7$PollenHeadWidth
N1 <- dataset7$NectarHeadWidth

#Compare Microbe vs Control, paired t-test
t.test(P1,N1,paired=FALSE) #t = 3.2077, df = 16, p-value = 0.005488
shapiro.test(P1) #Not normal
shapiro.test(N1) #Not normal
var.test(P1,N1) #Not normal

#Inputs:
HW <- dataset7$HeadWidth
PC <- dataset7$PollenCollected


#Model: wilcoxon-signed rank test
wilcox.test(HW ~ PC, paired=FALSE, data=dataset7) #W = 150, p-value = 0.005008





#
#Linear Model
#
#Urban flowers width versus haemocytometer microbe cell count
#

#Inputs:
MicrobeCount <- dataset8$MicrobeCount #Number of microbes on flower
FlowerWidth <- dataset8$FlowerWidth #widest part of flower in mm
FlowerSpecies <- dataset8$FlowerSpecies
#FlowerID <- dataset6$FlowerID #Repeated measure

#Model
CountByFLWRSize.m <-lm(log(MicrobeCount) ~ log(FlowerWidth) * FlowerSpecies, data=dataset8) #Need to put in FlowerID as repeated measure?

#Summary stats
summary(CountByFLWRSize.m)

#See if the residuals are distributed normally for each group. 
#Must be done when a continuous distribution (like with lmer, anova etc)
CountByFLWRSize.resid = resid(CountByFLWRSize.m)

#Plot the residuals
hist(CountByFLWRSize.resid)

#qq plot to see if residuals normally distributed
qqPlot(CountByFLWRSize.m, main="QQ Plot")

#Check to see the effects of each factor
#Anova(specialvsbouts.m, type = "3") #To look specifically at the interaction
Anova(CountByFLWRSize.m, type = "2")  #Multiple R-squared:  0.1221,	Adjusted R-squared:  0.1034 , F-statistic: 6.536 on 1 and 47 DF,  p-value: 0.01386
#log(HeadWidth)                  0.738  1  0.9415    0.3371    
#CollectedPollen                44.986  1 57.3545 1.451e-09 ***
#log(HeadWidth):CollectedPollen  0.107  1  0.1363    0.7137   