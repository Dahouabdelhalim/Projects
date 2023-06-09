setwd("/Users/Mickey/Documents/Acorn Woodpecker Research/TPA_Manuscript/PRSB/Revision_1/")
df1 <- read.csv("20180405_2016TriadicAwareness_Data_for_pub.csv")

#loading necessary packages
library(survival)
library(coxme)
library(MASS)
library(corrplot)
library(ggplot2)
library(survminer)
library(lme4)
library(lmerTest)
library(emmeans)
library(CircStats)


########### SETTING FACTORS AS FACTORS ##############
df1$Treatment <- as.factor(df1$Treatment)
df1$Female <- as.factor(df1$Female)
df1$Stimulus <- as.factor(df1$Stimulus)

############# CREATING HELMERT-CODED VERSIONS OF COVARIATES #############
#Helmert coding = take each value, subtract the mid-point of the data, then divide
#by 1/2 of the range
#This is done so all the covariates will be on the same scale, since when I tried
#running some of the models with the original covariates I got a warning message
#saying that they were on very different scales

helmert <- function(x) {
  midpoint <- (max(x)+min(x))/2
  range <- max(x)-min(x)
  (x-midpoint)/(0.5*range)
}

df1$ACSD.helm <- helmert(df1$Avg.Caller.Subj.Dist)
df1$Overlap.helm <- helmert(df1$PropOverlap)
df1$Lag.helm <- helmert(df1$StimLagTime)
df1$sdA.helm <- helmert(df1$sdA)
df1$Tot_Dur.helm <- helmert(df1$Tot_Dur)

############### CORRELATION MATRIX ####################

#Create correlation matrix of response variables

#Create a df of the relevant response variables
response <- cbind(df1$F.Time.First.Flight,
                       df1$F.Time.First.Positive.Flight,
                       df1$F.Time.React,
                       df1$F.Time.Closest.Approach,
                       df1$F.Closest.Distance,
                       df1$F.First.Flight.Direction,
                       df1$G.Num.Birds.Rally, 
                       df1$G.Num.Birds.Approaching, 
                  df1$G.Time.First.Call,
                       df1$G.WK.Diff,
                       df1$O.Time.Leave.Tree,
                       df1$O.Time.Closest.Approach,
                       df1$O.Closest.Distance)
colnames(response) <- c("F.Lat.Flight", "F.Lat.Pos.Fl", "F.Lat.React", "F.Lat.Closest", 
                             "F.Dist", "F.Direction",
                             "G.Rally", "G.Num.App", "G.Lat.Call", "G.WK.Rate",
                             "O.Lat.Leave.Tr", "O.Lat.Closest", "O.Dist")

#create the correlation matrix, deleting NA values on case by case basis
CorMat <- cor(response, use="complete.obs")
corrplot(CorMat, method="number", tl.pos="lt", tl.srt=90, number.cex=0.6)


############## SURVIVAL ANALYSES (COX REG) ##############

#Cox regression for F.Time.First.Positive.Flight
F.LPF.obj <- Surv(time=df1$F.Time.First.Positive.Flight, 
                  event=df1$F.Positive.Flight.Censor)
F.LPF.coxme1 <- coxme(F.LPF.obj ~ Treatment + (1|Female) + (1|Stimulus), data=df1)
summary(F.LPF.coxme1)
anova(F.LPF.coxme1)

#Now with all covariates
F.LPF.coxme2 <- coxme(F.LPF.obj ~ Treatment + (1|Female) + (1|Stimulus) + 
                        ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                        Tot_Dur.helm, data=df1)
summary(F.LPF.coxme2)
anova(F.LPF.coxme2) #Treatment still significant. Only sig covar is Tot_Dur.helm


#Cox regression for G.Time.First.Call
G.LCall.obj <- Surv(time=df1$G.Time.First.Call, 
                    event=df1$G.First.Call.Censor)
G.LCall.coxme1 <- coxme(G.LCall.obj ~ Treatment + 
                          (1|Female) + (1|Stimulus), data=df1)
summary(G.LCall.coxme1)
anova(G.LCall.coxme1)

#Now with all covariates
G.LCall.coxme2 <- coxme(G.LCall.obj ~ Treatment + 
                         (1|Female) + (1|Stimulus) + 
                         ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                         Tot_Dur.helm, data=df1)
summary(G.LCall.coxme2)
anova(G.LCall.coxme2) #Treatment p=0.572. Overlap.helm p=0.03784

############## Linear Mixed Model ##################

#LMM for O.Closest.Distance
#rank-transforming response variable
##tied ranks are randomly assigned an order relative to each other
### so there are no ties
LMM.O.Dist1 <- lmer(rank(O.Closest.Distance, na.last=FALSE, 
                         ties.method="random") ~ 
                      Treatment + (1|Female) + 
                      (1|Stimulus), data=df1)
summary(LMM.O.Dist1)
anova(LMM.O.Dist1)

#Checking assumptions for LMM.O.Dist1
hist(resid(LMM.O.Dist1))
plot(predict(LMM.O.Dist1), resid(LMM.O.Dist1))

#Now with all covariates
LMM.O.Dist2 <- lmer(rank(O.Closest.Distance, na.last=FALSE, 
                        ties.method="random") ~ 
                     Treatment + (1|Female) + 
                     (1|Stimulus) + 
                     ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                     Tot_Dur.helm, data=df1)
summary(LMM.O.Dist2)
anova(LMM.O.Dist2) #Can't compute p-value in lmerTest


########LMM for WKDiff

#Use linear mixed model with response variable sqrt transformed
LMM.WKDiff1 <- lmer(sqrt(G.WK.Diff) ~ Treatment + (1|Female) + 
                      (1|Stimulus), data=df1)

summary(LMM.WKDiff1)
anova(LMM.WKDiff1)

#Checking model assumptions for LMM.WKDiff
hist(resid(LMM.WKDiff))  #checking to see if residuals are normally distrib
plot(predict(LMM.WKDiff), resid(LMM.WKDiff)) #plotting resid vs. predicted values

#Now with all covariates
LMM.WKDiff2 <- lmer(sqrt(G.WK.Diff) ~ Treatment + (1|Female) + 
                       (1|Stimulus) + 
                     ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                     Tot_Dur.helm, data=df1)

summary(LMM.WKDiff2)
anova(LMM.WKDiff2) #nothing sig


############# GLMMs ################

#Binomial regression for G.Num.Birds.Approaching
Binom.NumApp1 <- glmer(cbind(G.Num.Birds.Approaching, 
                             G.Group.Size-G.Num.Birds.Approaching) ~ 
                         Treatment + 
                         (1|Female) + 
                         (1|Stimulus), 
                       family=binomial, data=df1)

Binom.NumApp1.reduc <- glmer(cbind(G.Num.Birds.Approaching, 
                                   G.Group.Size-G.Num.Birds.Approaching) ~
                               (1|Female) + 
                               (1|Stimulus), 
                             family=binomial, data=df1)

summary(Binom.NumApp1)

anova(Binom.NumApp1, Binom.NumApp1.reduc)

emmeans(Binom.NumApp1, ~Treatment, type="response")

#Checking assumptions for Binom.NumApp1
table(df1$G.Num.Birds.Approaching > 0, df1$Treatment)


#Now with all covariates
Binom.NumApp2 <- glmer(cbind(G.Num.Birds.Approaching, 
                            G.Group.Size-G.Num.Birds.Approaching) ~ 
                        Treatment + 
                       (1|Female) + 
                       (1|Stimulus) + 
                        ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                        Tot_Dur.helm, 
                     family=binomial, data=df1)

Binom.NumApp2.reduc <- glmer(cbind(G.Num.Birds.Approaching, 
                            G.Group.Size-G.Num.Birds.Approaching) ~
                        (1|Female) + 
                        (1|Stimulus) + 
                          ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                          Tot_Dur.helm, 
                      family=binomial, data=df1)

summary(Binom.NumApp2)

anova(Binom.NumApp2, Binom.NumApp2.reduc)

emmeans(Binom.NumApp2, ~Treatment, type="response")

##################
# Calculating the p-values for the covariates for G.Num.Birds.Approaching

#ACSD.helm
NBA.ACSD <- glmer(cbind(G.Num.Birds.Approaching, 
                            G.Group.Size-G.Num.Birds.Approaching) ~ 
                        Treatment + 
                        (1|Female) + 
                        (1|Stimulus) + 
                        ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                        Tot_Dur.helm, 
                      family=binomial, data=df1)

NBA.ACSD.reduc <- glmer(cbind(G.Num.Birds.Approaching, 
                              G.Group.Size-G.Num.Birds.Approaching) ~ 
                          Treatment + 
                          (1|Female) + 
                          (1|Stimulus) +
                          Overlap.helm + Lag.helm + sdA.helm + 
                          Tot_Dur.helm, 
                        family=binomial, data=df1)

anova(NBA.ACSD, NBA.ACSD.reduc) #p=0.944

#Overlap.helm
NBA.overlap <- glmer(cbind(G.Num.Birds.Approaching, 
                        G.Group.Size-G.Num.Birds.Approaching) ~ 
                    Treatment + 
                    (1|Female) + 
                    (1|Stimulus) + 
                    ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                    Tot_Dur.helm, 
                  family=binomial, data=df1)

NBA.overlap.reduc <- glmer(cbind(G.Num.Birds.Approaching, 
                           G.Group.Size-G.Num.Birds.Approaching) ~ 
                       Treatment + 
                       (1|Female) + 
                       (1|Stimulus) + 
                       ACSD.helm + Lag.helm + sdA.helm + 
                       Tot_Dur.helm, 
                     family=binomial, data=df1)

anova(NBA.overlap, NBA.overlap.reduc) #p=0.4997

#Lag.helm
NBA.Lag <- glmer(cbind(G.Num.Birds.Approaching, 
                           G.Group.Size-G.Num.Birds.Approaching) ~ 
                       Treatment + 
                       (1|Female) + 
                       (1|Stimulus) + 
                       ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                       Tot_Dur.helm, 
                     family=binomial, data=df1)

NBA.Lag.reduc <- glmer(cbind(G.Num.Birds.Approaching, 
                       G.Group.Size-G.Num.Birds.Approaching) ~ 
                   Treatment + 
                   (1|Female) + 
                   (1|Stimulus) + 
                   ACSD.helm + Overlap.helm + sdA.helm + 
                   Tot_Dur.helm, 
                 family=binomial, data=df1)

anova(NBA.Lag, NBA.Lag.reduc) #p=0.00114

#sdA.helm
NBA.sdA <- glmer(cbind(G.Num.Birds.Approaching, 
                       G.Group.Size-G.Num.Birds.Approaching) ~ 
                   Treatment + 
                   (1|Female) + 
                   (1|Stimulus) + 
                   ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                   Tot_Dur.helm, 
                 family=binomial, data=df1)

NBA.sdA.reduc <- glmer(cbind(G.Num.Birds.Approaching, 
                       G.Group.Size-G.Num.Birds.Approaching) ~ 
                   Treatment + 
                   (1|Female) + 
                   (1|Stimulus) + 
                   ACSD.helm + Overlap.helm + Lag.helm + 
                   Tot_Dur.helm, 
                 family=binomial, data=df1)

anova(NBA.sdA, NBA.sdA.reduc) #p=0.3741

#Tot_Dur.helm
NBA.Tot_Dur <- glmer(cbind(G.Num.Birds.Approaching, 
                           G.Group.Size-G.Num.Birds.Approaching) ~ 
                       Treatment + 
                       (1|Female) + 
                       (1|Stimulus) + 
                       ACSD.helm + Overlap.helm + Lag.helm + sdA.helm + 
                       Tot_Dur.helm, 
                     family=binomial, data=df1)

NBA.Tot_Dur.reduc <- glmer(cbind(G.Num.Birds.Approaching, 
                           G.Group.Size-G.Num.Birds.Approaching) ~ 
                       Treatment + 
                       (1|Female) + 
                       (1|Stimulus) + 
                       ACSD.helm + Overlap.helm + Lag.helm + sdA.helm, 
                     family=binomial, data=df1)

anova(NBA.Tot_Dur, NBA.Tot_Dur.reduc) #p=0.0083

###### Testing to see if mean angular moment differed b/t T and C stimuli ########
########## angular moment = metric of call synchrony
x <- df1[df1$Treatment=="T",]
y <- df1[df1$Treatment=="C",]

#in package "CircStats"
watson.two(x$CircMean, y$CircMean, plot=TRUE) #p>0.10


############ Code for creating survival curves for Figure 1 ###########

#Create survival fit object, which has survival probability and 95CIs
LPFfit <- survfit(F.LPF.obj ~ Treatment, data=df1)

#### Plot with survminer
survcurve <- ggsurvplot(LPFfit, fun="event", conf.int=T, 
                        font.x=c(13),
                        font.y=c(13),
                        legend="none",
                        legend.title="",
                        legend.labs=c("Control", "Test"),
                        linetype=c("twodash","solid"))
survcurve$plot <- survcurve$plot + 
  scale_x_continuous(breaks=c(0,30,60,90,120,150,180)) + 
  xlab("Time Since Start of Playback (s)") + 
  ylab("Cumulative Probability of a Positive Flight") + 
  ggplot2::annotate("text", x=188, y=0.935, label="Test", size=4.5) + 
  ggplot2::annotate("text", x=192, y=0.60, label="Control", size=4.5)
survcurve


