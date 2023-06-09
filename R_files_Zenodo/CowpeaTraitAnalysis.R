###########################################################################################################################
# COWPEA SANCTIONS EXPERIMENT - COMPLETE R CODE
###########################################################################################################################
###########################################################################################################################

#############################################
#PACKAGES TO DOWNLOAD AND INSTALL
#############################################
install.packages("Rmisc")
install.packages("ggplot2", dependencies=TRUE)
install.packages ("grid")
install.packages ("gridExtra")
install.packages ("grid")
install.packages("MASS")
install.packages("lme4")
install.packages("agridat")
install.packages("car")
install.packages("lsmeans")
install.packages("lmerTest")
install.packages("lme4")
install.packages("agricolae")
install.packages( "multcomp")
install.packages ("emmeans")
install.packages("pls")
install.packages("multcompView")
install.packages("dplyr")
install.packages("ggthemes")

library(ggplot2)
library(gridExtra)
library(grid)
install.packages("ggpubr")
library(ggpubr)
install.packages("Rmisc")
library(Rmisc)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library (car) ## Needed for testing homogeneity of variances
library(agricolae) # Needed to obtain Tukey groupings
library(MASS)
library(lme4)
library(lsmeans)
library(lmerTest)
library(foreign)
library(multcomp)
library(multcompView)
library(pls) ## Needed for partial least-square analyses
library(dplyr)
library(ggthemes)

##############################################
#TESTING MODELS FOR TABLES
##############################################
#Number of Nodules
##############################################
#Reading data into R
CowpeasAll <- read.csv(file.choose(), header = TRUE, sep=",", strip.white = TRUE)
#Subset data for specific analysis
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules))
#Transforming variable to achieve normality
CowpeasNoControlsNumberOfNodules$SqrtNON <-sqrt(CowpeasNoControlsNumberOfNodules$NumberOfNodules)
#Checking for Normality
ggqqplot(sqrt(CowpeasNoControlsNumberOfNodules$SqrtNON))
ggqqplot(CowpeasNoControlsNumberOfNodules$SqrtNON)
hist(CowpeasNoControlsNumberOfNodules$SqrtNON)
qqnorm(CowpeasNoControlsNumberOfNodules$SqrtNON)
qqline(CowpeasNoControlsNumberOfNodules$SqrtNON)

#Linear Mixed Model
CowpeasAll$sqrtNON <-sqrt(CowpeasAll$NumberOfNodules)

CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, sqrtNON))
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, Treatment == "A" | Treatment == "L" | Treatment == "AL", 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, sqrtNON))
#Model without line as random factor
lm.sqrtNON1 <- lm(sqrtNON ~ DaysSinceInoculation + Treatment*Genepool, CowpeasNoControlsNumberOfNodules)
#Model with line as random factor
lmer.sqrtNON <- lmer(sqrtNON ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON)
anova(lmer.sqrtNON,lm.sqrtNON1)
#Model withouth line and Treatment random factors
lmer.sqrtNON2 <- lmer(sqrtNON ~ DaysSinceInoculation + Treatment*Genepool + (Treatment|Line), CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON)
anova(lmer.sqrtNON2,lmer.sqrtNON)
#Model two is significant we will use model two

#Post hoc tests

lsmeans(lmer.sqrtNON2, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.sqrtNON2, pairwise ~ Treatment|Genepool, adjust ="bonferroni")
lsmeans(lmer.sqrtNON2, pairwise ~ Genepool, adjust ="bonferroni")
lsmeans(lmer.sqrtNON2, pairwise ~ Genepool*Treatment, adjust ="bonferroni")

#################################################
#Soil Community Number of Nodules
################################################

#Linear Mixed Model
CowpeasAll$sqrtNON <-sqrt(CowpeasAll$NumberOfNodules)

CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, sqrtNON))
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, Treatment == "S", 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, sqrtNON))

#Model withouth line as random factor
lmer.sqrtNON <- lmer(sqrtNON ~ DaysSinceInoculation + Genepool + (1|Line), CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON)

#Model without line as random factor
lm.sqrtNON1 <- lm(sqrtNON ~ DaysSinceInoculation + Genepool, CowpeasNoControlsNumberOfNodules)
anova(lmer.sqrtNON,lm.sqrtNON1)
#Post hoc tests

lsmeans(lmer.sqrtNON, pairwise ~ Genepool, adjust ="bonferroni")

#################################################
# Dry Nodule Weight
################################################

CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0, 
                                             select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L", 
                                             select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected <- CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected+1

CowpeasAllDryNoduleWeightCorrected$LogDry.Nodule.Weight.Corrected <-log10(CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected)


#New
lmer.DNB <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB)
lmer.DNB2 <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Treatment*Genepool + (Treatment|Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB)
anova(lmer.DNB,lmer.DNB2)

#Model without line as random factor
lm.DNB <- lm(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Treatment*Genepool, CowpeasAllDryNoduleWeightCorrected)
anova(lmer.DNB, lm.DNB)


#Lsmeans
lsmeans(lmer.DNB, pairwise ~ Genepool*Treatment, adjust ="bonferroni")
lsmeans(lmer.DNB, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.DNB, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.DNB, pairwise ~ Genepool, adjust ="bonferroni")
##############################################
# Soil Community Dry Nodule Weight
###############################################

CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0, 
                                             select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Treatment == "S", 
                                             select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected <- CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected+1

CowpeasAllDryNoduleWeightCorrected$LogDry.Nodule.Weight.Corrected <-log10(CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected)


#New
lmer.DNB <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Genepool + (1|Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB)

#Model without line as random factor
lm.DNB <- lm(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Genepool, CowpeasAllDryNoduleWeightCorrected)
anova(lmer.DNB, lm.DNB)

#Lsmeans
lsmeans(lmer.DNB, pairwise ~ Genepool, adjust ="bonferroni")

##############################################
# Investment
##############################################

#Subset of data for Investment analysis
CowpeasInvestment <- subset(CowpeasAll, Investment > 0,
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

CowpeasInvestment <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L",
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

#Model using interaction between treatment and line:
lmer.I <- lmer(Investment ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), data = CowpeasInvestment)
Anova(lmer.I)
lmer.I2 <- lmer(Investment ~ DaysSinceInoculation + Treatment*Genepool + (Treatment|Line), data = CowpeasInvestment)
Anova(lmer.I2)
anova(lmer.I,lmer.I2)

lm.I <- lm(Investment ~ DaysSinceInoculation + Treatment*Genepool, CowpeasInvestment)
anova(lmer.I, lm.I)

#Lsmeans
lsmeans(lmer.I, pairwise ~ Genepool*Treatment, adjust ="bonferroni")
lsmeans(lmer.I, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.I, pairwise ~ Genepool, adjust ="bonferroni")
lsmeans(lmer.I, pairwise ~ Treatment, adjust ="bonferroni")

##############################################
# Investment SOil COmmunity
##############################################

#Subset of data for Investment analysis
CowpeasInvestment <- subset(CowpeasAll, Investment > 0,
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

CowpeasInvestment <- subset(CowpeasAll, Treatment == "S",
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

#Model using interaction between treatment and line:
lmer.I <- lmer(Investment ~ DaysSinceInoculation + Genepool + (1|Line), data = CowpeasInvestment)
Anova(lmer.I)

lm.I <- lm(Investment ~ DaysSinceInoculation + Genepool, CowpeasInvestment)
anova(lmer.I, lm.I)


#Lsmeans
lsmeans(lmer.I, pairwise ~ Genepool, adjust ="bonferroni")


###############################################
# HOST GROWTH RESPONSE
###############################################

#Subset data for specific analysis
CowpeasHostGrowth <- subset(CowpeasAll, HGR > 0 | HGR < 0, 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, HGR))


CowpeasHostGrowth <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L", 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, HGR))

CowpeasHostGrowth$log.HGR <- log((CowpeasHostGrowth$HGR) + 1 - min(CowpeasHostGrowth$HGR))

#Checking for Normality

ggqqplot(CowpeasHostGrowth$HGR)
ggqqplot(CowpeasHostGrowth$HGR)
hist(CowpeasHostGrowth$HGR)
qqnorm(CowpeasHostGrowth$HGR)
qqline(CowpeasHostGrowth$HGR)

shapiro.test(CowpeasHostGrowth$HGR)

ggqqplot(CowpeasHostGrowth$log.HGR)
ggqqplot(CowpeasHostGrowth$log.HGR)
hist(CowpeasHostGrowth$log.HGR)
qqnorm(CowpeasHostGrowth$log.HGR)
qqline(CowpeasHostGrowth$log.HGR)

shapiro.test(CowpeasHostGrowth$log.HGR)

#Model using interaction between treatment and line:
lmer.HG <- lmer(log.HGR ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), data = CowpeasHostGrowth)
Anova(lmer.HG)
lmer.HG2 <- lmer(log.HGR ~ DaysSinceInoculation + Treatment*Genepool + (Treatment|Line), data = CowpeasHostGrowth)
Anova(lmer.HG2)
anova(lmer.HG,lmer.HG2)

lm.HGR <- lm(log.HGR ~ DaysSinceInoculation + Treatment*Genepool, CowpeasHostGrowth)
anova(lmer.HG, lm.HGR)

#Lsmeans
lsmeans(lmer.HG, pairwise ~ Genepool*Treatment, adjust ="bonferroni")
lsmeans(lmer.HG, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.HG, pairwise ~ Genepool, adjust ="bonferroni")
lsmeans(lmer.HG, pairwise ~ Treatment, adjust ="bonferroni")

##############################################
# Host Growth Response - Soil Community
##############################################
#Subset data for specific analysis
CowpeasHostGrowth <- subset(CowpeasAll, HGR > 0 | HGR < 0, 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, HGR))


CowpeasHostGrowth <- subset(CowpeasAll, Treatment == "S", 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, HGR))

CowpeasHostGrowth$log.HGR <- log((CowpeasHostGrowth$HGR) + 1 - min(CowpeasHostGrowth$HGR))

#Checking for Normality

ggqqplot(CowpeasHostGrowth$HGR)
ggqqplot(CowpeasHostGrowth$HGR)
hist(CowpeasHostGrowth$HGR)
qqnorm(CowpeasHostGrowth$HGR)
qqline(CowpeasHostGrowth$HGR)

shapiro.test(CowpeasHostGrowth$HGR)

ggqqplot(CowpeasHostGrowth$log.HGR)
ggqqplot(CowpeasHostGrowth$log.HGR)
hist(CowpeasHostGrowth$log.HGR)
qqnorm(CowpeasHostGrowth$log.HGR)
qqline(CowpeasHostGrowth$log.HGR)

shapiro.test(CowpeasHostGrowth$log.HGR)

#Model using interaction between treatment and line:
lmer.HG <- lmer(log.HGR ~ DaysSinceInoculation + Genepool + (1|Line), data = CowpeasHostGrowth)
Anova(lmer.HG)
lm.HGR <- lm(log.HGR ~ DaysSinceInoculation + Genepool, CowpeasHostGrowth)
anova(lmer.HG, lm.HGR)

#Lsmeans
lsmeans(lmer.HG, pairwise ~ Genepool, adjust ="bonferroni")

###############################################
#MEAN NODULE BIOMASS
##############################################


CowpeasAllMNW <- subset(CowpeasAll, DryWeightPerNodule > 0, 
                        select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))

#subset for Fix + and Fix -
MNWplusminus <- subset(CowpeasAllMNW, Treatment == "A" | Treatment == "L", 
                       select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))

#Transforming data
MNWplusminus$DryWeightPerNodule <- MNWplusminus$DryWeightPerNodule+1
MNWplusminus$logMNW <-log(MNWplusminus$DryWeightPerNodule)

#Normality  tests
ggqqplot(MNWplusminus$DryWeightPerNodule)
ggqqplot(MNWplusminus$DryWeightPerNodule)
hist(MNWplusminus$DryWeightPerNodule)
qqnorm(MNWplusminus$DryWeightPerNodule)
qqline(MNWplusminus$DryWeightPerNodule)

shapiro.test(MNWplusminus$DryWeightPerNodule)

ggqqplot(MNWplusminus$logMNW)
ggqqplot(MNWplusminus$logMNW)
hist(MNWplusminus$logMNW)
qqnorm(MNWplusminus$logMNW)
qqline(MNWplusminus$logMNW)

shapiro.test(MNWplusminus$logMNW)

#New
#Model withouth line and treatment as random factors
lmer.MNW1 <- lmer(logMNW ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), MNWplusminus)
Anova(lmer.MNW1)

lmer.MNW2 <- lmer(logMNW ~ DaysSinceInoculation + Treatment*Genepool + (Treatment|Line), data = MNWplusminus)
Anova(lmer.MNW2)

anova(lmer.MNW2,lmer.MNW1) #Interaction between line and treatment is significant

lm.MNB <- lm(logMNW ~ DaysSinceInoculation + Treatment*Genepool, MNWplusminus)
anova(lmer.MNW1, lm.MNB)

#Lsmeans
lsmeans(lmer.MNW1, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.MNW1, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.MNW1, pairwise ~ Genepool, adjust ="bonferroni")
lsmeans(lmer.MNW1, pairwise ~ Genepool*Treatment, adjust ="bonferroni")
###############################################
#MEAN NODULE BIOMASS - SOIL COMMUNITY
##############################################

CowpeasAllMNW <- subset(CowpeasAll, DryWeightPerNodule > 0, 
                        select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))

#subset for Fix + and Fix -
MNWplusminus <- subset(CowpeasAllMNW, Treatment == "S", 
                       select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))

#Transforming data
MNWplusminus$DryWeightPerNodule <- MNWplusminus$DryWeightPerNodule+1
MNWplusminus$logMNW <-log(MNWplusminus$DryWeightPerNodule)

#New
#Model withouth line and treatment as random factors
lmer.MNW1 <- lmer(logMNW ~ DaysSinceInoculation + Genepool + (1|Line), MNWplusminus)
Anova(lmer.MNW1)


lm.MNB <- lm(logMNW ~ DaysSinceInoculation + Genepool, MNWplusminus)
anova(lmer.MNW1, lm.MNB)

#Lsmeans
lsmeans(lmer.MNW1, pairwise ~ Genepool, adjust ="bonferroni")

################################################################
# Delta 15N
################################################################

################################################
# Delta fifthteen N
##############################################
# since we are using the delta 15N values we make a subset to avoid using the controls.
#Subset for analysis
fifthteenND <- subset(CowpeasAll, fifthteenN > 0 | fifthteenN < 0,
                      select=c(Line, Treatment, fifthteenN, Genepool, DaysSinceInoculation))

fifthteenNC <- subset(fifthteenND, Treatment == "A" | Treatment == "L",
                      select=c(Line, Treatment, fifthteenN, Genepool, DaysSinceInoculation))
#Model using interaction between treatment and line:
lmer.fifthteenN <- lmer(fifthteenN ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), data = fifthteenNC)
Anova(lmer.fifthteenN)

lm.fifthteenN <- lm(fifthteenN ~ DaysSinceInoculation + Treatment*Genepool, fifthteenNC)
anova(lmer.fifthteenN, lm.fifthteenN)
#Lsmeans
lsmeans(lmer.fifthteenN, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.fifthteenN, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.fifthteenN, pairwise ~ Genepool, adjust ="bonferroni")

################################################
# 4 Genepool Analysis
################################################
#Number of Nodules
##############################################
#Reading data into R
CowpeasAll <- read.csv(file.choose(), header = TRUE, sep=",", strip.white = TRUE)
#Subset data for specific analysis
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules))
#Transforming variable to achieve normality
CowpeasNoControlsNumberOfNodules$SqrtNON <-sqrt(CowpeasNoControlsNumberOfNodules$NumberOfNodules)
#Checking for Normality
ggqqplot(sqrt(CowpeasNoControlsNumberOfNodules$SqrtNON))
ggqqplot(CowpeasNoControlsNumberOfNodules$SqrtNON)
hist(CowpeasNoControlsNumberOfNodules$SqrtNON)
qqnorm(CowpeasNoControlsNumberOfNodules$SqrtNON)
qqline(CowpeasNoControlsNumberOfNodules$SqrtNON)

#Linear Mixed Model
CowpeasAll$sqrtNON <-sqrt(CowpeasAll$NumberOfNodules)

CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules, sqrtNON))
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, Treatment == "A" | Treatment == "L" | Treatment == "AL", 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules, sqrtNON))
#Model without line as random factor
lm.sqrtNON1 <- lm(sqrtNON ~ DaysSinceInoculation + Treatment*Groups, CowpeasNoControlsNumberOfNodules)
#Model with line as random factor
lmer.sqrtNON <- lmer(sqrtNON ~ DaysSinceInoculation + Treatment*Groups + (1|Line), CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON)
anova(lmer.sqrtNON,lm.sqrtNON1)
#Model withouth line and Treatment random factors
lmer.sqrtNON2 <- lmer(sqrtNON ~ DaysSinceInoculation + Treatment*Groups + (Treatment|Line), CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON2)
anova(lmer.sqrtNON2,lmer.sqrtNON)
#Model two is significant we will use model two

#Post hoc tests

lsmeans(lmer.sqrtNON2, pairwise ~ Groups|Treatment, adjust ="bonferroni")
lsmeans(lmer.sqrtNON2, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.sqrtNON2, pairwise ~ Groups, adjust ="bonferroni")
lsmeans(lmer.sqrtNON2, pairwise ~ Groups*Treatment, adjust ="bonferroni")

#################################################
#Soil Community Number of Nodules
################################################

#Linear Mixed Model
CowpeasAll$sqrtNON <-sqrt(CowpeasAll$NumberOfNodules)

CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules, sqrtNON))
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, Treatment == "S", 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules, sqrtNON))

#Model withouth line as random factor
lmer.sqrtNON <- lmer(sqrtNON ~ DaysSinceInoculation + Groups + (1|Line), CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON)

#Model without line as random factor
lm.sqrtNON1 <- lm(sqrtNON ~ DaysSinceInoculation + Groups, CowpeasNoControlsNumberOfNodules)
anova(lmer.sqrtNON,lm.sqrtNON1)
#Post hoc tests

lsmeans(lmer.sqrtNON, pairwise ~ Groups, adjust ="bonferroni")


#################################################
# Dry Nodule Weight
################################################

CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0, 
                                             select=c(Line, Treatment, Groups, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L", 
                                             select=c(Line, Treatment, Groups, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected <- CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected+1

CowpeasAllDryNoduleWeightCorrected$LogDry.Nodule.Weight.Corrected <-log10(CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected)


#New
lmer.DNB <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Treatment*Groups + (1|Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB)
lmer.DNB2 <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Treatment*Groups + (Treatment|Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB)
anova(lmer.DNB,lmer.DNB2)

#Model without line as random factor
lm.DNB <- lm(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Treatment*Groups, CowpeasAllDryNoduleWeightCorrected)
anova(lmer.DNB, lm.DNB)


#Lsmeans
lsmeans(lmer.DNB, pairwise ~ Groups*Treatment, adjust ="bonferroni")
lsmeans(lmer.DNB, pairwise ~ Groups|Treatment, adjust ="bonferroni")
lsmeans(lmer.DNB, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.DNB, pairwise ~ Groups, adjust ="bonferroni")
##############################################
# Soil Community Dry Nodule Weight
###############################################

CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0, 
                                             select=c(Line, Treatment, Groups, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Treatment == "S", 
                                             select=c(Line, Treatment, Groups, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected <- CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected+1

CowpeasAllDryNoduleWeightCorrected$LogDry.Nodule.Weight.Corrected <-log10(CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected)


#New
lmer.DNB <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Groups + (1|Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB)

#Model without line as random factor
lm.DNB <- lm(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Groups, CowpeasAllDryNoduleWeightCorrected)
anova(lmer.DNB, lm.DNB)

#Lsmeans
lsmeans(lmer.DNB, pairwise ~ Groups, adjust ="bonferroni")

##############################################
# Investment
##############################################

#Subset of data for Investment analysis
CowpeasInvestment <- subset(CowpeasAll, Investment > 0,
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, Investment))

CowpeasInvestment <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L",
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, Investment))

#Model using interaction between treatment and line:
lmer.I <- lmer(Investment ~ DaysSinceInoculation + Treatment*Groups + (1|Line), data = CowpeasInvestment)
Anova(lmer.I)
lmer.I2 <- lmer(Investment ~ DaysSinceInoculation + Treatment*Groups + (Treatment|Line), data = CowpeasInvestment)
Anova(lmer.I2)
anova(lmer.I,lmer.I2)

lm.I <- lm(Investment ~ DaysSinceInoculation + Treatment*Groups, CowpeasInvestment)
anova(lmer.I, lm.I)

#Lsmeans
lsmeans(lmer.I, pairwise ~ Groups*Treatment, adjust ="bonferroni")
lsmeans(lmer.I, pairwise ~ Groups|Treatment, adjust ="bonferroni")
lsmeans(lmer.I, pairwise ~ Groups, adjust ="bonferroni")
lsmeans(lmer.I, pairwise ~ Treatment, adjust ="bonferroni")

##############################################
# Investment SOil COmmunity
##############################################

#Subset of data for Investment analysis
CowpeasInvestment <- subset(CowpeasAll, Investment > 0,
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, Investment))

CowpeasInvestment <- subset(CowpeasAll, Treatment == "S",
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, Investment))

#Model using interaction between treatment and line:
lmer.I <- lmer(Investment ~ DaysSinceInoculation + Groups + (1|Line), data = CowpeasInvestment)
Anova(lmer.I)

lm.I <- lm(Investment ~ DaysSinceInoculation + Groups, CowpeasInvestment)
anova(lmer.I, lm.I)


#Lsmeans
lsmeans(lmer.I, pairwise ~ Groups, adjust ="bonferroni")


###############################################
# HOST GROWTH RESPONSE
###############################################

#Subset data for specific analysis
CowpeasHostGrowth <- subset(CowpeasAll, HGR > 0 | HGR < 0, 
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, HGR))


CowpeasHostGrowth <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L", 
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, HGR))

CowpeasHostGrowth$log.HGR <- log((CowpeasHostGrowth$HGR) + 1 - min(CowpeasHostGrowth$HGR))

#Checking for Normality

ggqqplot(CowpeasHostGrowth$HGR)
ggqqplot(CowpeasHostGrowth$HGR)
hist(CowpeasHostGrowth$HGR)
qqnorm(CowpeasHostGrowth$HGR)
qqline(CowpeasHostGrowth$HGR)

shapiro.test(CowpeasHostGrowth$HGR)

ggqqplot(CowpeasHostGrowth$log.HGR)
ggqqplot(CowpeasHostGrowth$log.HGR)
hist(CowpeasHostGrowth$log.HGR)
qqnorm(CowpeasHostGrowth$log.HGR)
qqline(CowpeasHostGrowth$log.HGR)

shapiro.test(CowpeasHostGrowth$log.HGR)

#Model using interaction between treatment and line:
lmer.HG <- lmer(log.HGR ~ DaysSinceInoculation + Treatment*Groups + (1|Line), data = CowpeasHostGrowth)
Anova(lmer.HG)
lmer.HG2 <- lmer(log.HGR ~ DaysSinceInoculation + Treatment*Groups + (Treatment|Line), data = CowpeasHostGrowth)
Anova(lmer.HG2)
anova(lmer.HG,lmer.HG2)

lm.HGR <- lm(log.HGR ~ DaysSinceInoculation + Treatment*Groups, CowpeasHostGrowth)
anova(lmer.HG, lm.HGR)

#Lsmeans
lsmeans(lmer.HG, pairwise ~ Groups*Treatment, adjust ="bonferroni")
lsmeans(lmer.HG, pairwise ~ Groups|Treatment, adjust ="bonferroni")
lsmeans(lmer.HG, pairwise ~ Groups, adjust ="bonferroni")
lsmeans(lmer.HG, pairwise ~ Treatment, adjust ="bonferroni")

##############################################
# Host Growth Response - Soil Community
##############################################
#Subset data for specific analysis
CowpeasHostGrowth <- subset(CowpeasAll, HGR > 0 | HGR < 0, 
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, HGR))


CowpeasHostGrowth <- subset(CowpeasAll, Treatment == "S", 
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, HGR))

CowpeasHostGrowth$log.HGR <- log((CowpeasHostGrowth$HGR) + 1 - min(CowpeasHostGrowth$HGR))

#Checking for Normality

ggqqplot(CowpeasHostGrowth$HGR)
ggqqplot(CowpeasHostGrowth$HGR)
hist(CowpeasHostGrowth$HGR)
qqnorm(CowpeasHostGrowth$HGR)
qqline(CowpeasHostGrowth$HGR)

shapiro.test(CowpeasHostGrowth$HGR)

ggqqplot(CowpeasHostGrowth$log.HGR)
ggqqplot(CowpeasHostGrowth$log.HGR)
hist(CowpeasHostGrowth$log.HGR)
qqnorm(CowpeasHostGrowth$log.HGR)
qqline(CowpeasHostGrowth$log.HGR)

shapiro.test(CowpeasHostGrowth$log.HGR)

#Model using interaction between treatment and line:
lmer.HG <- lmer(log.HGR ~ DaysSinceInoculation + Groups + (1|Line), data = CowpeasHostGrowth)
Anova(lmer.HG)
lm.HGR <- lm(log.HGR ~ DaysSinceInoculation + Groups, CowpeasHostGrowth)
anova(lmer.HG, lm.HGR)

#Lsmeans
lsmeans(lmer.HG, pairwise ~ Groups, adjust ="bonferroni")

###############################################
#MEAN NODULE BIOMASS
##############################################


CowpeasAllMNW <- subset(CowpeasAll, DryWeightPerNodule > 0, 
                        select=c(Line, Treatment, Groups, DaysSinceInoculation, DryWeightPerNodule))

#subset for Fix + and Fix -
MNWplusminus <- subset(CowpeasAllMNW, Treatment == "A" | Treatment == "L", 
                       select=c(Line, Treatment, Groups, DaysSinceInoculation, DryWeightPerNodule))

#Transforming data
MNWplusminus$DryWeightPerNodule <- MNWplusminus$DryWeightPerNodule+1
MNWplusminus$logMNW <-log(MNWplusminus$DryWeightPerNodule)

#Normality  tests
ggqqplot(MNWplusminus$DryWeightPerNodule)
ggqqplot(MNWplusminus$DryWeightPerNodule)
hist(MNWplusminus$DryWeightPerNodule)
qqnorm(MNWplusminus$DryWeightPerNodule)
qqline(MNWplusminus$DryWeightPerNodule)

shapiro.test(MNWplusminus$DryWeightPerNodule)

ggqqplot(MNWplusminus$logMNW)
ggqqplot(MNWplusminus$logMNW)
hist(MNWplusminus$logMNW)
qqnorm(MNWplusminus$logMNW)
qqline(MNWplusminus$logMNW)

shapiro.test(MNWplusminus$logMNW)

#New
#Model withouth line and treatment as random factors
lmer.MNW1 <- lmer(logMNW ~ DaysSinceInoculation + Treatment*Groups + (1|Line), MNWplusminus)
Anova(lmer.MNW1)

lmer.MNW2 <- lmer(logMNW ~ DaysSinceInoculation + Treatment*Groups + (Treatment|Line), data = MNWplusminus)
Anova(lmer.MNW2)

anova(lmer.MNW2,lmer.MNW1) #Interaction between line and treatment is not significant we can use the first model

lm.MNB <- lm(logMNW ~ DaysSinceInoculation + Treatment*Groups, MNWplusminus)
anova(lmer.MNW1, lm.MNB)

#Lsmeans
lsmeans(lmer.MNW1, pairwise ~ Groups|Treatment, adjust ="bonferroni")
lsmeans(lmer.MNW1, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.MNW1, pairwise ~ Groups, adjust ="bonferroni")
lsmeans(lmer.MNW1, pairwise ~ Groups*Treatment, adjust ="bonferroni")
###############################################
#MEAN NODULE BIOMASS - SOIL COMMUNITY
##############################################

CowpeasAllMNW <- subset(CowpeasAll, DryWeightPerNodule > 0, 
                        select=c(Line, Treatment, Groups, DaysSinceInoculation, DryWeightPerNodule))

#subset for Fix + and Fix -
MNWplusminus <- subset(CowpeasAllMNW, Treatment == "S", 
                       select=c(Line, Treatment, Groups, DaysSinceInoculation, DryWeightPerNodule))

#Transforming data
MNWplusminus$DryWeightPerNodule <- MNWplusminus$DryWeightPerNodule+1
MNWplusminus$logMNW <-log(MNWplusminus$DryWeightPerNodule)

#New
#Model withouth line and treatment as random factors
lmer.MNW1 <- lmer(logMNW ~ DaysSinceInoculation + Groups + (1|Line), MNWplusminus)
Anova(lmer.MNW1)


lm.MNB <- lm(logMNW ~ DaysSinceInoculation + Groups, MNWplusminus)
anova(lmer.MNW1, lm.MNB)

#Lsmeans
lsmeans(lmer.MNW1, pairwise ~ Groups, adjust ="bonferroni")

################################################
# Delta fifthteen N
##############################################
# since we are using the delta 15N values we make a subset to avoid using the controls.
#Subset for analysis
fifthteenND <- subset(CowpeasAll, fifthteenN > 0 | fifthteenN < 0,
                      select=c(Line, Treatment, fifthteenN, Genepool, DaysSinceInoculation))

fifthteenNC <- subset(fifthteenND, Treatment == "A" | Treatment == "L",
                      select=c(Line, Treatment, fifthteenN, Genepool, DaysSinceInoculation))

#Model using interaction between treatment and line:
lmer.fifthteenN <- lmer(fifthteenN ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), data = fifthteenNC)
Anova(lmer.fifthteenN)

lm.fifthteenN <- lm(fifthteenN ~ DaysSinceInoculation + Treatment*Genepool, fifthteenNC)
anova(lmer.fifthteenN, lm.fifthteenN)
#Lsmeans
lsmeans(lmer.fifthteenN, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.fifthteenN, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.fifthteenN, pairwise ~ Genepool, adjust ="bonferroni")


###################################################################
##################################################################
#     MEANS
###################################################################

###################################################################
# Number of Nodules
###################################################################

#Subsets for raw means

Numberofnodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                          select=c(Line, Treatment, Genepool, NumberOfNodules))

NONW <- subset(Numberofnodules, Genepool == "Wild", 
               select=c(Line, Treatment, Genepool, NumberOfNodules))
NONC1 <- subset(Numberofnodules, Genepool == "One", 
                select=c(Line, Treatment, Genepool, NumberOfNodules))
NONC2 <- subset(Numberofnodules, Genepool == "Two", 
                select=c(Line, Treatment, Genepool, NumberOfNodules))
#Mean & SE

mean(Numberofnodules$NumberOfNodules)
sem<-sd(Numberofnodules$NumberOfNodules)/sqrt(length(Numberofnodules$NumberOfNodules))
sem
mean(NONW$NumberOfNodules)
sem<-sd(NONW$NumberOfNodules)/sqrt(length(NONW$NumberOfNodules))
sem
mean(NONC1$NumberOfNodules)
sem<-sd(NONC1$NumberOfNodules)/sqrt(length(NONC1$NumberOfNodules))
sem
mean(NONC2$NumberOfNodules)
sem<-sd(NONC2$NumberOfNodules)/sqrt(length(NONC2$NumberOfNodules))
sem

#All lines by treatment Total

NON <- subset(Numberofnodules, Treatment == "A", 
              select=c(Line, Treatment, Genepool, NumberOfNodules))
mean(NON$NumberOfNodules)
sem<-sd(NON$NumberOfNodules)/sqrt(length(NON$NumberOfNodules))
sem

NON <- subset(Numberofnodules, Treatment == "L", 
              select=c(Line, Treatment, Genepool,NumberOfNodules))
mean(NON$NumberOfNodules)
sem<-sd(NON$NumberOfNodules)/sqrt(length(NON$NumberOfNodules))
sem

NON <- subset(Numberofnodules, Treatment == "AL", 
              select=c(Line, Treatment, Genepool,NumberOfNodules))
mean(NON$NumberOfNodules)
sem<-sd(NON$NumberOfNodules)/sqrt(length(NON$NumberOfNodules))
sem

NON <- subset(Numberofnodules, Treatment == "S", 
              select=c(Line, Treatment, Genepool,NumberOfNodules))
mean(NON$NumberOfNodules)
sem<-sd(NON$NumberOfNodules)/sqrt(length(NON$NumberOfNodules))
sem


#Genepool x Treatment 

NONWA <- subset(NONW, Treatment == "A", 
                select=c(Line, Treatment, NumberOfNodules))
mean(NONWA$NumberOfNodules)
sem<-sd(NONWA$NumberOfNodules)/sqrt(length(NONWA$NumberOfNodules))
sem
NONWL <- subset(NONW, Treatment == "L", 
                select=c(Line, Treatment, NumberOfNodules))
mean(NONWL$NumberOfNodules)
sem<-sd(NONWL$NumberOfNodules)/sqrt(length(NONWL$NumberOfNodules))
sem

NONWAL <- subset(NONW, Treatment == "AL", 
                 select=c(Line, Treatment, NumberOfNodules))
mean(NONWAL$NumberOfNodules)
sem<-sd(NONWAL$NumberOfNodules)/sqrt(length(NONWAL$NumberOfNodules))
sem

NONC1A <- subset(NONC1, Treatment == "A", 
                 select=c(Line, Treatment, NumberOfNodules))
mean(NONC1A$NumberOfNodules)
sem<-sd(NONC1A$NumberOfNodules)/sqrt(length(NONC1A$NumberOfNodules))
sem
NONC1L <- subset(NONC1, Treatment == "L", 
                 select=c(Line, Treatment, NumberOfNodules))
mean(NONC1L$NumberOfNodules)
sem<-sd(NONC1L$NumberOfNodules)/sqrt(length(NONC1L$NumberOfNodules))
sem

NONC1AL <- subset(NONC1, Treatment == "AL", 
                  select=c(Line, Treatment, NumberOfNodules))
mean(NONC1AL$NumberOfNodules)
sem<-sd(NONC1AL$NumberOfNodules)/sqrt(length(NONC1AL$NumberOfNodules))
sem

NONC2A <- subset(NONC2, Treatment == "A", 
                 select=c(Line, Treatment, NumberOfNodules))
mean(NONC2A$NumberOfNodules)
sem<-sd(NONC2A$NumberOfNodules)/sqrt(length(NONC2A$NumberOfNodules))
sem
NONC2L <- subset(NONC2, Treatment == "L", 
                 select=c(Line, Treatment, NumberOfNodules))
mean(NONC2L$NumberOfNodules)
sem<-sd(NONC2L$NumberOfNodules)/sqrt(length(NONC2L$NumberOfNodules))
sem
NONC2AL <- subset(NONC2, Treatment == "AL", 
                  select=c(Line, Treatment, NumberOfNodules))
mean(NONC2AL$NumberOfNodules)
sem<-sd(NONC2AL$NumberOfNodules)/sqrt(length(NONC2AL$NumberOfNodules))
sem


NONWS <- subset(NONW, Treatment == "S", 
                select=c(Line, Treatment, NumberOfNodules))
mean(NONWS$NumberOfNodules)
sem<-sd(NONWS$NumberOfNodules)/sqrt(length(NONWS$NumberOfNodules))
sem


NONC1S <- subset(NONC1, Treatment == "S", 
                 select=c(Line, Treatment, NumberOfNodules))
mean(NONC1S$NumberOfNodules)
sem<-sd(NONC1S$NumberOfNodules)/sqrt(length(NONC1S$NumberOfNodules))
sem

NONC2S <- subset(NONC2, Treatment == "S", 
                  select=c(Line, Treatment, NumberOfNodules))
mean(NONC2S$NumberOfNodules)
sem<-sd(NONC2S$NumberOfNodules)/sqrt(length(NONC2S$NumberOfNodules))
sem


####################################################################
# Investment
####################################################################

#Investment
#Subsets 

CowpeasInvestment <- subset(CowpeasAll, Investment > 0, 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

CowpeasMeansInvestmentW <- subset(CowpeasInvestment, Genepool == "Wild", 
                                  select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))
CowpeasMeansInvestmentC1 <- subset(CowpeasInvestment, Genepool == "One", 
                                   select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))
CowpeasMeansInvestmentC2 <- subset(CowpeasInvestment, Genepool == "Two", 
                                   select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

#Mean & SE

mean(CowpeasInvestment$Investment)
sem<-sd(CowpeasInvestment$Investment)/sqrt(length(CowpeasInvestment$Investment))
sem
mean(CowpeasMeansInvestmentW$Investment)
sem<-sd(CowpeasMeansInvestmentW$Investment)/sqrt(length(CowpeasMeansInvestmentW$Investment))
sem
mean(CowpeasMeansInvestmentC1$Investment)
sem<-sd(CowpeasMeansInvestmentC1$Investment)/sqrt(length(CowpeasMeansInvestmentC1$Investment))
sem
mean(CowpeasMeansInvestmentC2$Investment)
sem<-sd(CowpeasMeansInvestmentC2$Investment)/sqrt(length(CowpeasMeansInvestmentC2$Investment))
sem

#Genepool x Treatment 

CowpeasMeansInvestmentWCo <- subset(CowpeasMeansInvestmentW, Treatment == "AL", 
                                    select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentWCo$Investment)
sem<-sd(CowpeasMeansInvestmentWCo$Investment)/sqrt(length(CowpeasMeansInvestmentWCo$Investment))
sem
CowpeasMeansInvestmentWA <- subset(CowpeasMeansInvestmentW, Treatment == "A", 
                                   select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentWA$Investment)
sem<-sd(CowpeasMeansInvestmentWA$Investment)/sqrt(length(CowpeasMeansInvestmentWA$Investment))
sem
CowpeasMeansInvestmentWL <- subset(CowpeasMeansInvestmentW, Treatment == "L", 
                                   select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentWL$Investment)
sem<-sd(CowpeasMeansInvestmentWL$Investment)/sqrt(length(CowpeasMeansInvestmentWL$Investment))
sem

CowpeasMeansInvestmentC1Co <- subset(CowpeasMeansInvestmentC1, Treatment == "AL", 
                                     select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentC1Co$Investment)
sem<-sd(CowpeasMeansInvestmentC1Co$Investment)/sqrt(length(CowpeasMeansInvestmentC1Co$Investment))
sem
CowpeasMeansInvestmentC1A <- subset(CowpeasMeansInvestmentC1, Treatment == "A", 
                                    select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentC1A$Investment)
sem<-sd(CowpeasMeansInvestmentC1A$Investment)/sqrt(length(CowpeasMeansInvestmentC1A$Investment))
sem
CowpeasMeansInvestmentC1L <- subset(CowpeasMeansInvestmentC1, Treatment == "L", 
                                    select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentC1L$Investment)
sem<-sd(CowpeasMeansInvestmentC1L$Investment)/sqrt(length(CowpeasMeansInvestmentC1L$Investment))
sem

CowpeasMeansInvestmentC2Co <- subset(CowpeasMeansInvestmentC2, Treatment == "AL", 
                                     select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentC2Co$Investment)
sem<-sd(CowpeasMeansInvestmentC2Co$Investment)/sqrt(length(CowpeasMeansInvestmentC2Co$Investment))
sem
CowpeasMeansInvestmentC2A <- subset(CowpeasMeansInvestmentC2, Treatment == "A", 
                                    select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentC2A$Investment)
sem<-sd(CowpeasMeansInvestmentC2A$Investment)/sqrt(length(CowpeasMeansInvestmentC2A$Investment))
sem
CowpeasMeansInvestmentC2L <- subset(CowpeasMeansInvestmentC2, Treatment == "L", 
                                    select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentC2L$Investment)
sem<-sd(CowpeasMeansInvestmentC2L$Investment)/sqrt(length(CowpeasMeansInvestmentC2L$Investment))
sem



CowpeasMeansInvestmentWS <- subset(CowpeasMeansInvestmentW, Treatment == "S", 
                                    select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentWS$Investment)
sem<-sd(CowpeasMeansInvestmentWS$Investment)/sqrt(length(CowpeasMeansInvestmentWS$Investment))
sem

CowpeasMeansInvestmentC1S <- subset(CowpeasMeansInvestmentC1, Treatment == "S", 
                                    select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentC1S$Investment)
sem<-sd(CowpeasMeansInvestmentC1S$Investment)/sqrt(length(CowpeasMeansInvestmentC1S$Investment))
sem

CowpeasMeansInvestmentC2S <- subset(CowpeasMeansInvestmentC2, Treatment == "S", 
                                    select=c(Line, Treatment, Investment))
mean(CowpeasMeansInvestmentC2S$Investment)
sem<-sd(CowpeasMeansInvestmentC2S$Investment)/sqrt(length(CowpeasMeansInvestmentC2S$Investment))
sem

#All lines by treatment Total

CowpeasInvestment <- subset(CowpeasAll, Treatment == "AL", 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))
mean(CowpeasInvestment$Investment)
sem<-sd(CowpeasInvestment$Investment)/sqrt(length(CowpeasInvestment$Investment))
sem

CowpeasInvestment <- subset(CowpeasAll, Treatment == "A", 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))
mean(CowpeasInvestment$Investment)
sem<-sd(CowpeasInvestment$Investment)/sqrt(length(CowpeasInvestment$Investment))
sem

CowpeasInvestment <- subset(CowpeasAll, Treatment == "L", 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))
mean(CowpeasInvestment$Investment)
sem<-sd(CowpeasInvestment$Investment)/sqrt(length(CowpeasInvestment$Investment))
sem

CowpeasInvestment <- subset(CowpeasAll, Treatment == "S", 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))
mean(CowpeasInvestment$Investment)
sem<-sd(CowpeasInvestment$Investment)/sqrt(length(CowpeasInvestment$Investment))
sem

###################################################################
# Dry Nodule Biomass
####################################################################


CowpeasDNW <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0,
                     select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

#Subsets 

CowpeasNoControlsDryNoduleWeight <- subset(CowpeasDNW, 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasMeansDryNoduleWeightW <- subset(CowpeasNoControlsDryNoduleWeight, Genepool == "Wild", 
                                       select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasMeansDryNoduleWeightC1 <- subset(CowpeasNoControlsDryNoduleWeight, Genepool == "One", 
                                        select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasMeansDryNoduleWeightC2 <- subset(CowpeasNoControlsDryNoduleWeight, Genepool == "Two", 
                                        select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

#Mean & SE

mean(CowpeasNoControlsDryNoduleWeight$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasNoControlsDryNoduleWeight$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasNoControlsDryNoduleWeight$Dry.Nodule.Weight.Corrected))
sem
mean(CowpeasMeansDryNoduleWeightW$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightW$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightW$Dry.Nodule.Weight.Corrected))
sem
mean(CowpeasMeansDryNoduleWeightC1$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC1$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC1$Dry.Nodule.Weight.Corrected))
sem
mean(CowpeasMeansDryNoduleWeightC2$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC2$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC2$Dry.Nodule.Weight.Corrected))
sem

#All lines by treatment Total

CowpeasNoControlsDryNoduleWeightAdjusted <- subset(CowpeasDNW, Treatment == "AL", 
                                                   select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
mean(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected))
sem

CowpeasNoControlsDryNoduleWeightAdjusted <- subset(CowpeasDNW, Treatment == "A", 
                                                   select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
mean(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected))
sem

CowpeasNoControlsDryNoduleWeightAdjusted <- subset(CowpeasDNW, Treatment == "L", 
                                                   select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
mean(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected))
sem

CowpeasNoControlsDryNoduleWeightAdjusted <- subset(CowpeasDNW, Treatment == "S", 
                                                   select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
mean(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasNoControlsDryNoduleWeightAdjusted$Dry.Nodule.Weight.Corrected))
sem

#Genepool x Treatment 

CowpeasMeansDryNoduleWeightWCo <- subset(CowpeasMeansDryNoduleWeightW, Treatment == "AL", 
                                         select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightWCo$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightWCo$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightWCo$Dry.Nodule.Weight.Corrected))
sem
CowpeasMeansDryNoduleWeightWA <- subset(CowpeasMeansDryNoduleWeightW, Treatment == "A", 
                                        select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightWA$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightWA$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightWA$Dry.Nodule.Weight.Corrected))
sem
CowpeasMeansDryNoduleWeightWL <- subset(CowpeasMeansDryNoduleWeightW, Treatment == "L", 
                                        select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightWL$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightWL$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightWL$Dry.Nodule.Weight.Corrected))
sem

CowpeasMeansDryNoduleWeightC1Co <- subset(CowpeasMeansDryNoduleWeightC1, Treatment == "AL", 
                                          select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightC1Co$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC1Co$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC1Co$Dry.Nodule.Weight.Corrected))
sem
CowpeasMeansDryNoduleWeightC1A <- subset(CowpeasMeansDryNoduleWeightC1, Treatment == "A", 
                                         select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightC1A$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC1A$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC1A$Dry.Nodule.Weight.Corrected))
sem
CowpeasMeansDryNoduleWeightC1L <- subset(CowpeasMeansDryNoduleWeightC1, Treatment == "L", 
                                         select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightC1L$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC1L$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC1L$Dry.Nodule.Weight.Corrected))
sem

CowpeasMeansDryNoduleWeightC2Co <- subset(CowpeasMeansDryNoduleWeightC2, Treatment == "AL", 
                                          select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightC2Co$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC2Co$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC2Co$Dry.Nodule.Weight.Corrected))
sem
CowpeasMeansDryNoduleWeightC2A <- subset(CowpeasMeansDryNoduleWeightC2, Treatment == "A", 
                                         select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightC2A$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC2A$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC2A$Dry.Nodule.Weight.Corrected))
sem
CowpeasMeansDryNoduleWeightC2L <- subset(CowpeasMeansDryNoduleWeightC2, Treatment == "L", 
                                         select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightC2L$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC2L$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC2L$Dry.Nodule.Weight.Corrected))
sem



CowpeasMeansDryNoduleWeightWS <- subset(CowpeasMeansDryNoduleWeightW, Treatment == "S", 
                                         select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightWS$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightWS$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightWS$Dry.Nodule.Weight.Corrected))
sem

CowpeasMeansDryNoduleWeightC1S <- subset(CowpeasMeansDryNoduleWeightC1, Treatment == "S", 
                                          select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightC1S$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC1S$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC1S$Dry.Nodule.Weight.Corrected))
sem

CowpeasMeansDryNoduleWeightC2S <- subset(CowpeasMeansDryNoduleWeightC2, Treatment == "S", 
                                         select=c(Line, Treatment, Dry.Nodule.Weight.Corrected))
mean(CowpeasMeansDryNoduleWeightC2S$Dry.Nodule.Weight.Corrected)
sem<-sd(CowpeasMeansDryNoduleWeightC2S$Dry.Nodule.Weight.Corrected)/sqrt(length(CowpeasMeansDryNoduleWeightC2S$Dry.Nodule.Weight.Corrected))
sem



###################################################################
# Mean Nodule Biomass
###################################################################

#Means and SE for genepool under each treatment
#Subsets 

CowpeasAllMNW <- subset(CowpeasAll, DryWeightPerNodule > 0, 
                        select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))

#subset for Fix + and Fix -

MNWplusminus <- subset(CowpeasAllMNW, Treatment == "A" | Treatment == "L", 
                       select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))


MNWW <- subset(MNWplusminus, Genepool == "Wild", 
               select=c(Line, Treatment, Genepool, DryWeightPerNodule))
MNWC1 <- subset(MNWplusminus, Genepool == "One", 
                select=c(Line, Treatment, Genepool, DryWeightPerNodule))
MNWC2 <- subset(MNWplusminus, Genepool == "Two", 
                select=c(Line, Treatment, Genepool, DryWeightPerNodule))
#Mean & SE

mean(MNWplusminus$DryWeightPerNodule)
sem<-sd(MNWplusminus$DryWeightPerNodule)/sqrt(length(MNWplusminus$DryWeightPerNodule))
sem
mean(MNWW$DryWeightPerNodule)
sem<-sd(MNWW$DryWeightPerNodule)/sqrt(length(MNWW$DryWeightPerNodule))
sem
mean(MNWC1$DryWeightPerNodule)
sem<-sd(MNWC1$DryWeightPerNodule)/sqrt(length(MNWC1$DryWeightPerNodule))
sem
mean(MNWC2$DryWeightPerNodule)
sem<-sd(MNWC2$DryWeightPerNodule)/sqrt(length(MNWC2$DryWeightPerNodule))
sem

#All lines by treatment Total

MNWARS <- subset(MNWplusminus, Treatment == "A", 
                 select=c(Line, Treatment, Genepool, DryWeightPerNodule))
mean(MNWARS$DryWeightPerNodule)
sem<-sd(MNWARS$DryWeightPerNodule)/sqrt(length(MNWARS$DryWeightPerNodule))
sem

MNWLI <- subset(MNWplusminus, Treatment == "L", 
                select=c(Line, Treatment, Genepool,DryWeightPerNodule))
mean(MNWLI$DryWeightPerNodule)
sem<-sd(MNWLI$DryWeightPerNodule)/sqrt(length(MNWLI$DryWeightPerNodule))
sem

#Genepool x Treatment 

MNWWARS <- subset(MNWW, Treatment == "A", 
                  select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWWARS$DryWeightPerNodule)
sem<-sd(MNWWARS$DryWeightPerNodule)/sqrt(length(MNWWARS$DryWeightPerNodule))
sem
MNWC1ARS <- subset(MNWC1, Treatment == "A", 
                   select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWC1ARS$DryWeightPerNodule)
sem<-sd(MNWC1ARS$DryWeightPerNodule)/sqrt(length(MNWC1ARS$DryWeightPerNodule))
sem

MNWC2ARS <- subset(MNWC2, Treatment == "A", 
                   select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWC2ARS$DryWeightPerNodule)
sem<-sd(MNWC2ARS$DryWeightPerNodule)/sqrt(length(MNWC2ARS$DryWeightPerNodule))
sem


MNWWLI <- subset(MNWW, Treatment == "L", 
                 select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWWLI$DryWeightPerNodule)
sem<-sd(MNWWLI$DryWeightPerNodule)/sqrt(length(MNWWLI$DryWeightPerNodule))
sem
MNWC1LI <- subset(MNWC1, Treatment == "L", 
                  select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWC1LI$DryWeightPerNodule)
sem<-sd(MNWC1LI$DryWeightPerNodule)/sqrt(length(MNWC1LI$DryWeightPerNodule))
sem

MNWC2LI <- subset(MNWC2, Treatment == "L", 
                  select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWC2LI$DryWeightPerNodule)
sem<-sd(MNWC2LI$DryWeightPerNodule)/sqrt(length(MNWC2LI$DryWeightPerNodule))
sem

MNWplusminus <- subset(CowpeasAllMNW, Treatment == "S", 
                       select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))


MNWWS <- subset(MNWplusminus, Genepool == "Wild", 
                 select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWWS$DryWeightPerNodule)
sem<-sd(MNWWS$DryWeightPerNodule)/sqrt(length(MNWWS$DryWeightPerNodule))
sem
MNWC1S <- subset(MNWplusminus, Genepool == "One", 
                  select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWC1S$DryWeightPerNodule)
sem<-sd(MNWC1S$DryWeightPerNodule)/sqrt(length(MNWC1S$DryWeightPerNodule))
sem

MNWC2S <- subset(MNWplusminus, Genepool == "Two", 
                  select=c(Line, Treatment, DryWeightPerNodule))
mean(MNWC2S$DryWeightPerNodule)
sem<-sd(MNWC2S$DryWeightPerNodule)/sqrt(length(MNWC2S$DryWeightPerNodule))
sem


###################################################################
#   Host Growth Response
##################################################################

## Subset for means

CowpeasHGR <- subset(CowpeasAll, HGR < 0 | HGR > 0,  
              select=c(Line, Treatment, Genepool, HGR, Groups, DaysSinceInoculation))

CowpeasHGR$log.HGR <- log((CowpeasHGR$HGR) + 1 - min(CowpeasHGR$HGR))
#Means and SE for genepool under each treatment
#Subsets 

HGW <- subset(CowpeasHGR, Genepool == "Wild", 
              select=c(Line, Treatment, Genepool, log.HGR))
HGC1 <- subset(CowpeasHGR, Genepool == "One", 
               select=c(Line, Treatment, Genepool, log.HGR))
HGC2 <- subset(CowpeasHGR, Genepool == "Two", 
               select=c(Line, Treatment, Genepool, log.HGR))
#Mean & SE

mean(data.cnoc$RHG)
sem<-sd(data.cnoc$RHG)/sqrt(length(data.cnoc$RHG))
sem
mean(HGW$RHG)
sem<-sd(HGW$RHG)/sqrt(length(HGW$RHG))
sem
mean(HGC1$RHG)
sem<-sd(HGC1$RHG)/sqrt(length(HGC1$RHG))
sem
mean(HGC2$RHG)
sem<-sd(HGC2$RHG)/sqrt(length(HGC2$RHG))
sem
##################### means for transformed data
mean(CowpeasHGR$log.HGR)
sem<-sd(CowpeasHGR$log.HGR)/sqrt(length(CowpeasHGR$log.HGR))
sem
mean(HGW$RHG)
sem<-sd(HGW$RHG)/sqrt(length(HGW$RHG))
sem
mean(HGC1$RHG)
sem<-sd(HGC1$RHG)/sqrt(length(HGC1$RHG))
sem
mean(HGC2$RHG)
sem<-sd(HGC2$RHG)/sqrt(length(HGC2$RHG))
sem

#All lines by treatment Total

HGARS <- subset(CowpeasHGR, Treatment == "A", 
                select=c(Line, Treatment, Genepool, log.HGR))
mean(HGARS$log.HGR)
sem<-sd(HGARS$log.HGR)/sqrt(length(HGARS$log.HGR))
sem

HGLI <- subset(CowpeasHGR, Treatment == "L", 
               select=c(Line, Treatment, Genepool,log.HGR))
mean(HGLI$log.HGR)
sem<-sd(HGLI$log.HGR)/sqrt(length(HGLI$log.HGR))
sem

HGCo <- subset(CowpeasHGR, Treatment == "AL", 
               select=c(Line, Treatment, Genepool,log.HGR))
mean(HGCo$log.HGR)
sem<-sd(HGCo$log.HGR)/sqrt(length(HGCo$log.HGR))
sem

HGS <- subset(CowpeasHGR, Treatment == "S", 
              select=c(Line, Treatment, Genepool,log.HGR))
mean(HGS$log.HGR)
sem<-sd(HGS$log.HGR)/sqrt(length(HGS$log.HGR))
sem

#Genepool x Treatment 


HGWARS <- subset(HGW, Treatment == "A", 
                 select=c(Line, Treatment, RHG))
mean(HGWARS$RHG)
sem<-sd(HGWARS$RHG)/sqrt(length(HGWARS$RHG))
sem
HGC1ARS <- subset(HGC1, Treatment == "A", 
                  select=c(Line, Treatment, RHG))
mean(HGC1ARS$RHG)
sem<-sd(HGC1ARS$RHG)/sqrt(length(HGC1ARS$RHG))
sem

HGC2ARS <- subset(HGC2, Treatment == "A", 
                  select=c(Line, Treatment, RHG))
mean(HGC2ARS$RHG)
sem<-sd(HGC2ARS$RHG)/sqrt(length(HGC2ARS$RHG))
sem


HGWLI <- subset(HGW, Treatment == "L", 
                select=c(Line, Treatment, RHG))
mean(HGWLI$RHG)
sem<-sd(HGWLI$RHG)/sqrt(length(HGWLI$RHG))
sem
HGC1LI <- subset(HGC1, Treatment == "L", 
                 select=c(Line, Treatment, RHG))
mean(HGC1LI$RHG)
sem<-sd(HGC1LI$RHG)/sqrt(length(HGC1LI$RHG))
sem

HGC2LI <- subset(HGC2, Treatment == "L", 
                 select=c(Line, Treatment, RHG))
mean(HGC2LI$RHG)
sem<-sd(HGC2LI$RHG)/sqrt(length(HGC2LI$RHG))
sem

HGWCo <- subset(HGW, Treatment == "AL", 
                select=c(Line, Treatment, RHG))
mean(HGWCo$RHG)
sem<-sd(HGWCo$RHG)/sqrt(length(HGWCo$RHG))
sem
HGC1Co <- subset(HGC1, Treatment == "AL", 
                 select=c(Line, Treatment, RHG))
mean(HGC1Co$RHG)
sem<-sd(HGC1Co$RHG)/sqrt(length(HGC1Co$RHG))
sem

HGC2Co <- subset(HGC2, Treatment == "AL", 
                 select=c(Line, Treatment, RHG))
mean(HGC2Co$RHG)
sem<-sd(HGC2Co$RHG)/sqrt(length(HGC2Co$RHG))
sem

#Genepool x Treatment (Means of transformed data)


HGWARS <- subset(HGW, Treatment == "A", 
                 select=c(Line, Treatment, log.HGR))
mean(HGWARS$log.HGR)
sem<-sd(HGWARS$log.HGR)/sqrt(length(HGWARS$log.HGR))
sem
HGC1ARS <- subset(HGC1, Treatment == "A", 
                  select=c(Line, Treatment, log.HGR))
mean(HGC1ARS$log.HGR)
sem<-sd(HGC1ARS$log.HGR)/sqrt(length(HGC1ARS$log.HGR))
sem

HGC2ARS <- subset(HGC2, Treatment == "A", 
                  select=c(Line, Treatment, log.HGR))
mean(HGC2ARS$log.HGR)
sem<-sd(HGC2ARS$log.HGR)/sqrt(length(HGC2ARS$log.HGR))
sem


HGWLI <- subset(HGW, Treatment == "L", 
                select=c(Line, Treatment, log.HGR))
mean(HGWLI$log.HGR)
sem<-sd(HGWLI$log.HGR)/sqrt(length(HGWLI$log.HGR))
sem
HGC1LI <- subset(HGC1, Treatment == "L", 
                 select=c(Line, Treatment, log.HGR))
mean(HGC1LI$log.HGR)
sem<-sd(HGC1LI$log.HGR)/sqrt(length(HGC1LI$log.HGR))
sem

HGC2LI <- subset(HGC2, Treatment == "L", 
                 select=c(Line, Treatment, log.HGR))
mean(HGC2LI$log.HGR)
sem<-sd(HGC2LI$log.HGR)/sqrt(length(HGC2LI$log.HGR))
sem

HGWCo <- subset(HGW, Treatment == "AL", 
                select=c(Line, Treatment, log.HGR))
mean(HGWCo$log.HGR)
sem<-sd(HGWCo$log.HGR)/sqrt(length(HGWCo$log.HGR))
sem
HGC1Co <- subset(HGC1, Treatment == "AL", 
                 select=c(Line, Treatment, log.HGR))
mean(HGC1Co$log.HGR)
sem<-sd(HGC1Co$log.HGR)/sqrt(length(HGC1Co$log.HGR))
sem

HGC2Co <- subset(HGC2, Treatment == "AL", 
                 select=c(Line, Treatment, log.HGR))
mean(HGC2Co$log.HGR)
sem<-sd(HGC2Co$log.HGR)/sqrt(length(HGC2Co$log.HGR))
sem


HGWS <- subset(HGW, Treatment == "S", 
               select=c(Line, Treatment, log.HGR))
mean(HGWS$log.HGR)
sem<-sd(HGWS$log.HGR)/sqrt(length(HGWS$log.HGR))
sem
HGC1S <- subset(HGC1, Treatment == "S", 
                select=c(Line, Treatment, log.HGR))
mean(HGC1S$log.HGR)
sem<-sd(HGC1S$log.HGR)/sqrt(length(HGC1S$log.HGR))
sem

HGC2S <- subset(HGC2, Treatment == "S", 
                select=c(Line, Treatment, log.HGR))
mean(HGC2S$log.HGR)
sem<-sd(HGC2S$log.HGR)/sqrt(length(HGC2S$log.HGR))
sem


###############################################################
# Delta 15N
##############################################################

fifthteenND <- subset(CowpeasAll, fifthteenN > 0 | fifthteenN < 0,
                      select=c(Line, Treatment, fifthteenN, Genepool, DaysSinceInoculation))

#All lines by treatment Total

fifthteenNFixplus <- subset(fifthteenND, Treatment == "A", 
                            select=c(Line, Treatment, Genepool, fifthteenN))
mean(fifthteenNFixplus $fifthteenN)
sem<-sd(fifthteenNFixplus $fifthteenN)/sqrt(length(fifthteenNFixplus $fifthteenN))
sem

fifthteenNFixminus <- subset(fifthteenND, Treatment == "L", 
                             select=c(Line, Treatment, Genepool,fifthteenN))
mean(fifthteenNFixminus$fifthteenN)
sem<-sd(fifthteenNFixminus$fifthteenN)/sqrt(length(fifthteenNFixminus$fifthteenN))
sem

#Genepool x Treatment 
fifthteenNWA <- subset(fifthteenNW, Treatment == "A", 
                       select=c(Line, Treatment, fifthteenN))
mean(fifthteenNWA$fifthteenN)
sem<-sd(fifthteenNWA$fifthteenN)/sqrt(length(fifthteenNWA$fifthteenN))
sem
fifthteenNWL <- subset(fifthteenNW, Treatment == "L", 
                       select=c(Line, Treatment, fifthteenN))
mean(fifthteenNWL$fifthteenN)
sem<-sd(fifthteenNWL$fifthteenN)/sqrt(length(fifthteenNWL$fifthteenN))
sem

fifthteenNC1A <- subset(fifthteenNC1, Treatment == "A", 
                        select=c(Line, Treatment, fifthteenN))
mean(fifthteenNC1A$fifthteenN)
sem<-sd(fifthteenNC1A$fifthteenN)/sqrt(length(fifthteenNC1A$fifthteenN))
sem
fifthteenNC1L <- subset(fifthteenNC1, Treatment == "L", 
                        select=c(Line, Treatment, fifthteenN))
mean(fifthteenNC1L$fifthteenN)
sem<-sd(fifthteenNC1L$fifthteenN)/sqrt(length(fifthteenNC1L$fifthteenN))
sem

fifthteenNC2A <- subset(fifthteenNC2, Treatment == "A", 
                        select=c(Line, Treatment, fifthteenN))
mean(fifthteenNC2A$fifthteenN)
sem<-sd(fifthteenNC2A$fifthteenN)/sqrt(length(fifthteenNC2A$fifthteenN))
sem
fifthteenNC2L <- subset(fifthteenNC2, Treatment == "L", 
                        select=c(Line, Treatment, fifthteenN))
mean(fifthteenNC2L$fifthteenN)
sem<-sd(fifthteenNC2L$fifthteenN)/sqrt(length(fifthteenNC2L$fifthteenN))
sem


#####################################################################################
#### Graphs new format
###################################################################################

#Reading data into R
CowpeasAll <- read.csv(file.choose(), header = TRUE, sep=",", strip.white = TRUE)
#Subset data for specific analysis
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, GT))

#Transforming variable to achieve normality
CowpeasNoControlsNumberOfNodules$SqrtNON <-sqrt(CowpeasNoControlsNumberOfNodules$NumberOfNodules)

#Linear Mixed Model
CowpeasAll$sqrtNON <-sqrt(CowpeasAll$NumberOfNodules)

CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, sqrtNON, GT))
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, Treatment == "A" | Treatment == "L" | Treatment == "AL", 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, sqrtNON, GT))

#New
#Model withouth line and treatment as random factors
lmer.sqrtNON <- lmer(sqrtNON ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON)
#Model without line as random factor to test for significance:
lm.sqrtNON <- lm(sqrtNON ~ DaysSinceInoculation + Treatment*Genepool, CowpeasNoControlsNumberOfNodules)
#Log-likelihood test for random factor
anova(lmer.sqrtNON,lm.sqrtNON) #random factor is significant
#Model using interaction between treatment and line:
lmer.sqrtNON2 <- lmer(sqrtNON ~ DaysSinceInoculation + Treatment*Genepool + (1|Line) + (1|Treatment:Line), data = CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON2)
anova(lmer.sqrtNON2)
anova(lmer.sqrtNON2,lmer.sqrtNON) #Interaction between line and treatment is significant

#Lsmeans
lsmeans(lmer.sqrtNON2, pairwise ~ Treatment|Genepool, adjust ="bonferroni")
lsmeans(lmer.sqrtNON2, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.sqrtNON2, pairwise ~ Genepool, adjust ="bonferroni")

lsmeans.NON <- emmeans(lmer.sqrtNON2, pairwise ~ Treatment|Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
NON.cld <- cld(lsmeans.NON$emmeans, Letters = letters)
NON.cld$.group=gsub(" ", "", NON.cld$.group)  ###  Remove spaces in .group 
NON.cld


#Plotting overall nodule number by species and treatment

non.plot <- ggplot(data = NON.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Treatment", values = c("A" = "#440154FF","AL" = "#FDE725FF", "L" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot

NON.cld$Genepool=factor(NON.cld$Genepool, levels=c("Wild",
                                                   "One",
                                                   "Two"))

non.plot <- ggplot(data = NON.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Treatment", values = c("A" = "#440154FF","AL" = "#FDE725FF", "L" = "#21908CFF")) +
  labs(x="Cowpea Genepools", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 

non.plot

######changing the colors of the plots so that they are not equal to the genepool referents used in Fig 1 and S1.

install.packages("viridis")
library(viridis)

non.plot <- ggplot(data = NON.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  scale_fill_viridis_d(option = "mako") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  cleanup +
  labs(x="Cowpea Genepools", y= "Sqrt Nodule Number")


non.plot

#Soil Community

#Subset data for specific analysis
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules))

#Transforming variable to achieve normality
CowpeasNoControlsNumberOfNodules$SqrtNON <-sqrt(CowpeasNoControlsNumberOfNodules$NumberOfNodules)

#Linear Mixed Model
CowpeasAll$sqrtNON <-sqrt(CowpeasAll$NumberOfNodules)

CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, sqrtNON))
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, Treatment == "S", 
                                           select=c(Line, Treatment, Genepool, DaysSinceInoculation, NumberOfNodules, sqrtNON))

#New
lmer.sqrtNON2 <- lmer(sqrtNON ~ DaysSinceInoculation + Genepool + (1|Line), data = CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON2)

#Lsmeans

lsmeans(lmer.sqrtNON2, pairwise ~ Genepool, adjust ="bonferroni")

lsmeans.NON <- emmeans(lmer.sqrtNON2, pairwise ~ Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
lsmeans.NON
NON.cld <- cld(lsmeans.NON$emmeans, Letters = letters)
NON.cld
NON.cld$.group=gsub(" ", "", NON.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
non.plot <- ggplot(data = NON.cld, mapping = aes(x=Genepool, y=emmean, fill=NON.cld$Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "gold", "One" = "darkslateblue", "Two" = "aquamarine4")) +
  labs(x="Inoculation Treatments", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot

NON.cld$Genepool=factor(NON.cld$Genepool, levels=c("Wild",
                                                   "One",
                                                   "Two"))

non.plot <- ggplot(data = NON.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4")) +
  labs(x="Inoculation Treatments", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot

################
# Investment
################

#Subset of data for Investment analysis
CowpeasInvestment <- subset(CowpeasAll, Investment > 0,
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

CowpeasInvestment <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L",
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

#Checking for Normality
ggqqplot(sqrt(CowpeasInvestment$Investment))
ggqqplot(CowpeasInvestment$Investment)
hist(CowpeasInvestment$Investment)
qqnorm(CowpeasInvestment$Investment)
qqline(CowpeasInvestment$Investment)

#Model using interaction between treatment and line:
lmer.I2 <- lmer(Investment ~ DaysSinceInoculation + Treatment*Genepool + (1|Line) + (1|Treatment:Line), data = CowpeasInvestment)
Anova(lmer.I2)

#Lsmeans

lsmeans.I <- emmeans(lmer.I2, pairwise ~ Treatment|Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
I.cld <- cld(lsmeans.I$emmeans, Letters = letters)
I.cld$.group=gsub(" ", "", I.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
I.plot <- ggplot(data = I.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Inoculation Treatments", y= "Investment") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
I.plot

I.cld$Genepool=factor(I.cld$Genepool, levels=c("Wild",
                                               "One",
                                               "Two"))

I.plot <- ggplot(data = I.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Inoculation Treatments", y= "Investment") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

I.plot
###############
# Soil Community
###############


#Subset of data for Investment analysis
CowpeasInvestment <- subset(CowpeasAll, Investment > 0,
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

CowpeasInvestment <- subset(CowpeasAll, Treatment == "S",
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, Investment))

#Model using interaction between treatment and line:
lmer.I2 <- lmer(Investment ~ DaysSinceInoculation + Genepool + (1|Line), data = CowpeasInvestment)
Anova(lmer.I2)

#Lsmeans

lsmeans.I <- emmeans(lmer.I2, pairwise ~ Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
I.cld <- cld(lsmeans.I$emmeans, Letters = letters)
I.cld$.group=gsub(" ", "", I.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
I.plot <- ggplot(data = I.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Investment") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
I.plot

I.cld$Genepool=factor(I.cld$Genepool, levels=c("Wild",
                                               "One",
                                               "Two"))

I.plot <- ggplot(data = I.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4")) +
  labs(x="Inoculation Treatments", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

I.plot

##############################
# Host Growth Response
##############################

#Subset data for specific analysis
CowpeasHostGrowth <- subset(CowpeasAll, HGR > 0 | HGR < 0, 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, HGR))


CowpeasHostGrowth <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L", 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, HGR))

CowpeasHostGrowth$log.HGR <- log((CowpeasHostGrowth$HGR) + 1 - min(CowpeasHostGrowth$HGR))
qqPlot(CowpeasHostGrowth$log.HGR)


#Model using interaction between treatment and line:
lmer.HG <- lmer(log.HGR ~ DaysSinceInoculation + Treatment*Genepool + (1|Line) + (1|Treatment:Line), data = CowpeasHostGrowth)
Anova(lmer.HG)

#Lsmeans
lsmeans(lmer.HG, pairwise ~ Genepool*Treatment, adjust ="bonferroni")

lsmeans.HG <- emmeans(lmer.HG, pairwise ~ Treatment|Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
HG.cld <- cld(lsmeans.HG$emmeans, Letters = letters)
HG.cld$.group=gsub(" ", "", HG.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 

non.plot <- ggplot(data = HG.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Host Growth Response (%)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot

HG.cld$Genepool=factor(HG.cld$Genepool, levels=c("Wild",
                                                 "One",
                                                 "Two"))

non.plot <- ggplot(data = HG.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Inoculation Treatments", y= "Investment") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot

########################################## Soil Community #################################

#Subset data for specific analysis
CowpeasHostGrowth <- subset(CowpeasAll, HGR > 0 | HGR < 0, 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, HGR))
CowpeasHostGrowth <- subset(CowpeasAll, Treatment == "S", 
                            select=c(Line, Treatment, Genepool, DaysSinceInoculation, HGR))
CowpeasHostGrowth$log.HGR <- log((CowpeasHostGrowth$HGR) + 1 - min(CowpeasHostGrowth$HGR))
qqPlot(CowpeasHostGrowth$log.HGR)


#Model using interaction between treatment and line:
lmer.HG <- lmer(log.HGR ~ DaysSinceInoculation + Genepool + (1|Line), data = CowpeasHostGrowth)
Anova(lmer.HG)

#Lsmeans
lsmeans(lmer.HG, pairwise ~ Genepool, adjust ="bonferroni")

lsmeans.HG <- emmeans(lmer.HG, pairwise ~ Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
HG.cld <- cld(lsmeans.HG$emmeans, Letters = letters)
HG.cld$.group=gsub(" ", "", HG.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 

non.plot <- ggplot(data = HG.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Host Growth Response (%)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot

non.plot <- ggplot(data = HG.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Host Growth Response (%)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) 

non.plot

HG.cld$Genepool=factor(HG.cld$Genepool, levels=c("Wild",
                                                 "One",
                                                 "Two"))

non.plot <- ggplot(data = HG.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4")) +
  labs(x="Inoculation Treatments", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot

######################################################
# Dry Nodule Weight
#####################################################

CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0, 
                                             select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L", 
                                             select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected <- CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected+1

CowpeasAllDryNoduleWeightCorrected$LogDry.Nodule.Weight.Corrected <-log10(CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected)


#New

lmer.DNB2 <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Treatment*Genepool + (1|Line) + (1|Treatment:Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB2)

#Lsmeans
lsmeans(lmer.DNB2, pairwise ~ Genepool*Treatment, adjust ="bonferroni")

lsmeans.DNB <- emmeans(lmer.DNB2, pairwise ~ Treatment|Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
DNB.cld <- cld(lsmeans.DNB$emmeans, Letters = letters)
DNB.cld$.group=gsub(" ", "", DNB.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
DNB.plot <- ggplot(data = DNB.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Dry Nodule Biomass") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
DNB.plot

DNB.plot <- ggplot(data = DNB.cld, mapping = aes(x=Treatment, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Dry Nodule Biomass") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
DNB.plot

DNB.cld$Genepool=factor(DNB.cld$Genepool, levels=c("Wild",
                                                   "One",
                                                   "Two"))

non.plot <- ggplot(data = DNB.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Inoculation Treatments", y= "Dry Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot

#############################
# Soil Community
############################


CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0, 
                                             select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Treatment == "S", 
                                             select=c(Line, Treatment, Genepool, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected <- CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected+1

CowpeasAllDryNoduleWeightCorrected$LogDry.Nodule.Weight.Corrected <-log10(CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected)


#New

lmer.DNB2 <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Genepool + (1|Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB2)

#Lsmeans
lsmeans(lmer.DNB2, pairwise ~ Genepool, adjust ="bonferroni")

lsmeans.DNB <- emmeans(lmer.DNB2, pairwise ~ Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
DNB.cld <- cld(lsmeans.DNB$emmeans, Letters = letters)
DNB.cld$.group=gsub(" ", "", DNB.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
DNB.plot <- ggplot(data = DNB.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Dry Nodule Biomass") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
DNB.plot

DNB.cld$Genepool=factor(DNB.cld$Genepool, levels=c("Wild",
                                                   "One",
                                                   "Two"))

non.plot <- ggplot(data = DNB.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4")) +
  labs(x="Inoculation Treatments", y= "Log Dry Nodule Biomass (g)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot
###################################################################
#   Mean Nodule Biomass
##############################################################

#Subset for specific values for analysis
CowpeasAllMNW <- subset(CowpeasAll, DryWeightPerNodule > 0, 
                        select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))

#subset for Fix + and Fix -

MNWplusminus <- subset(CowpeasAllMNW, Treatment == "A" | Treatment == "L", 
                       select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))

#Transforming data
MNWplusminus$DryWeightPerNodule <- MNWplusminus$DryWeightPerNodule+1
MNWplusminus$logMNW <-log(MNWplusminus$DryWeightPerNodule)

#New
#Model withouth line and treatment as random factors
lmer.MNW <- lmer(logMNW ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), MNWplusminus)
Anova(lmer.MNW)
#Model without line as random factor to test for significance:
lm.MNW <- lm(logMNW ~ DaysSinceInoculation + Treatment*Genepool, MNWplusminus)
#Log-likelihood test for random factor
anova(lmer.MNW,lm.MNW) #random factor is significant
#Model using interaction between treatment and line:
lmer.MNW2 <- lmer(logMNW ~ DaysSinceInoculation + Treatment*Genepool + (1|Line) + (1|Treatment:Line), data = MNWplusminus)
Anova(lmer.MNW2)
anova(lmer.MNW2)
anova(lmer.MNW2,lmer.MNW) #Interaction between line and treatment is significant

#Lsmeans
lsmeans(lmer.MNW2, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.MNW2, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.MNW2, pairwise ~ Genepool, adjust ="bonferroni")
lsmeans(lmer.MNW2, pairwise ~ Genepool*Treatment, adjust ="bonferroni")

lsmeans.MNW <- emmeans(lmer.MNW2, pairwise ~ Treatment|Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
MNW.cld <- cld(lsmeans.MNW$emmeans, Letters = letters)
MNW.cld
MNW.cld$.group=gsub(" ", "", MNW.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
MNW.plot <- ggplot(data = MNW.cld, mapping = aes(x=Treatment, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
MNW.plot

MNW.cld$Genepool=factor(MNW.cld$Genepool, levels=c("Wild",
                                                   "One",
                                                   "Two"))

non.plot <- ggplot(data = MNW.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Inoculation Treatments", y= "Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot

#############################################
# Soil Community
#############################################
#Subset for specific values for analysis

CowpeasAllMNW <- subset(CowpeasAll, DryWeightPerNodule > 0, 
                        select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))

#subset for Fix + and Fix -

MNWplusminus <- subset(CowpeasAllMNW, Treatment == "S", 
                       select=c(Line, Treatment, Genepool, DaysSinceInoculation, DryWeightPerNodule))


#Transforming data
MNWplusminus$DryWeightPerNodule <- MNWplusminus$DryWeightPerNodule+1
MNWplusminus$logMNW <-log(MNWplusminus$DryWeightPerNodule)

#New
#Model using interaction between treatment and line:
lmer.MNW2 <- lmer(logMNW ~ DaysSinceInoculation + Genepool + (1|Line), data = MNWplusminus)
Anova(lmer.MNW2)

#Lsmeans
lsmeans(lmer.MNW2, pairwise ~ Genepool, adjust ="bonferroni")
lsmeans.MNW <- emmeans(lmer.MNW2, pairwise ~ Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
MNW.cld <- cld(lsmeans.MNW$emmeans, Letters = letters)
MNW.cld
MNW.cld$.group=gsub(" ", "", MNW.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
MNW.plot <- ggplot(data = MNW.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("Wild" = "#FDE725FF", "One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
MNW.plot

MNW.cld$Genepool=factor(MNW.cld$Genepool, levels=c("Wild",
                                                   "One",
                                                   "Two"))

non.plot <- ggplot(data = MNW.cld, mapping = aes(x=Genepool, y=emmean, fill=Genepool)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4")) +
  labs(x="Inoculation Treatments", y= "Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot

################################################
# Delta fifthteen N
##############################################
# since we are using the delta 15N values we make a subset to avoid using the controls.
#Subset for analysis
fifthteenND <- subset(CowpeasAll, fifthteenN > 0 | fifthteenN < 0,
                      select=c(Line, Treatment, fifthteenN, Genepool, DaysSinceInoculation))

fifthteenNC <- subset(fifthteenND, Treatment == "A" | Treatment == "L",
                      select=c(Line, Treatment, fifthteenN, Genepool, DaysSinceInoculation))

#Model using interaction between treatment and line:
lmer.fifthteenN <- lmer(fifthteenN ~ DaysSinceInoculation + Treatment*Genepool + (1|Line), data = fifthteenNC)
Anova(lmer.fifthteenN)

lm.fifthteenN <- lm(fifthteenN ~ DaysSinceInoculation + Treatment*Genepool, fifthteenNC)
anova(lmer.fifthteenN, lm.fifthteenN)
#Lsmeans
lsmeans(lmer.fifthteenN, pairwise ~ Genepool|Treatment, adjust ="bonferroni")
lsmeans(lmer.fifthteenN, pairwise ~ Treatment, adjust ="bonferroni")
lsmeans(lmer.fifthteenN, pairwise ~ Genepool, adjust ="bonferroni")

lsmeans.fifthteenN <- emmeans(lmer.fifthteenN2, pairwise ~ Treatment|Genepool, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
fifthteenN.cld <- cld(lsmeans.fifthteenN$emmeans, Letters = letters)
fifthteenN.cld$.group=gsub(" ", "", fifthteenN.cld$.group)  ###  Remove spaces in .group 

fifthteenN.cld$Genepool=factor(fifthteenN.cld$Genepool, levels=c("Wild",
                                                                 "One",
                                                                 "Two"))

non.plot <- ggplot(data = fifthteenN.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Inoculation Treatments", y= "Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot

non.plot <- ggplot(data = fifthteenN.cld, mapping = aes(x=Genepool, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Inoculation Treatments", y= "Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup
non.plot

###################################################################
#  4 genepools analysis graphs
###################################################################
#Subset data for specific analysis
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules))

#Transforming variable to achieve normality
CowpeasNoControlsNumberOfNodules$SqrtNON <-sqrt(CowpeasNoControlsNumberOfNodules$NumberOfNodules)

#Linear Mixed Model
CowpeasAll$sqrtNON <-sqrt(CowpeasAll$NumberOfNodules)

CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules, sqrtNON))
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, Treatment == "A" | Treatment == "L" | Treatment == "AL", 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules, sqrtNON))

#New
#Model withouth line and treatment as random factors
lmer.sqrtNON <- lmer(sqrtNON ~ DaysSinceInoculation + Treatment*Groups + (1|Line), CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON)
#Model without line as random factor to test for significance:
lm.sqrtNON <- lm(sqrtNON ~ DaysSinceInoculation + Treatment*Groups, CowpeasNoControlsNumberOfNodules)
#Log-likelihood test for random factor
anova(lmer.sqrtNON,lm.sqrtNON) #random factor is significant
#Model using interaction between treatment and line:
lmer.sqrtNON2 <- lmer(sqrtNON ~ DaysSinceInoculation + Treatment*Groups + (1|Line) + (1|Treatment:Line), data = CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON2)
anova(lmer.sqrtNON2)
anova(lmer.sqrtNON2,lmer.sqrtNON) #Interaction between line and treatment is significant

#Lsmeans
lsmeans.NON <- emmeans(lmer.sqrtNON2, pairwise ~ Treatment|Groups, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
NON.cld <- cld(lsmeans.NON$emmeans, Letters = letters)
NON.cld$.group=gsub(" ", "", NON.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
non.plot <- ggplot(data = NON.cld, mapping = aes(x=Groups, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot


NON.cld$Groups=factor(NON.cld$Groups, levels=c("Wild1","Wild2",
                                               "One",
                                               "Two"))

non.plot <- ggplot(data = NON.cld, mapping = aes(x=Groups, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Cowpea Genepools", y= "Sqrt Number of Nodules") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot
#Angelas Data

#Subset data for specific analysis
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules))

#Transforming variable to achieve normality
CowpeasNoControlsNumberOfNodules$SqrtNON <-sqrt(CowpeasNoControlsNumberOfNodules$NumberOfNodules)

#Linear Mixed Model
CowpeasAll$sqrtNON <-sqrt(CowpeasAll$NumberOfNodules)

CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, NumberOfNodules > 0, 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules, sqrtNON))
CowpeasNoControlsNumberOfNodules <- subset(CowpeasAll, Treatment == "S", 
                                           select=c(Line, Treatment, Groups, DaysSinceInoculation, NumberOfNodules, sqrtNON))

#New
lmer.sqrtNON2 <- lmer(sqrtNON ~ DaysSinceInoculation + Groups + (1|Line), data = CowpeasNoControlsNumberOfNodules)
Anova(lmer.sqrtNON2)

#Lsmeans

lsmeans(lmer.sqrtNON2, pairwise ~ Groups, adjust ="bonferroni")

lsmeans.NON <- emmeans(lmer.sqrtNON2, pairwise ~ Groups, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
NON.cld <- cld(lsmeans.NON$emmeans, Letters = letters)
NON.cld$.group=gsub(" ", "", NON.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
non.plot <- ggplot(data = NON.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot

non.plot <- ggplot(data = NON.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Sqrt Nodule Number") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())

non.plot

NON.cld$Groups=factor(NON.cld$Groups, levels=c("Wild1","Wild2",
                                               "One",
                                               "Two"))

non.plot <- ggplot(data = NON.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4", "royalblue4")) +
  labs(x="Inoculation Treatments", y= "Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot

################
# Investment
################

#Subset of data for Investment analysis
CowpeasInvestment <- subset(CowpeasAll, Investment > 0,
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, Investment))

CowpeasInvestment <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L",
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, Investment))

#Model using interaction between treatment and line:
lmer.I2 <- lmer(Investment ~ DaysSinceInoculation + Treatment*Groups + (1|Line) + (1|Treatment:Line), data = CowpeasInvestment)
Anova(lmer.I2)

#Lsmeans

lsmeans.I <- emmeans(lmer.I2, pairwise ~ Treatment|Groups, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
I.cld <- cld(lsmeans.I$emmeans, Letters = letters)
I.cld$.group=gsub(" ", "", I.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
I.plot <- ggplot(data = I.cld, mapping = aes(x=Groups, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Investment") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
I.plot


I.cld$Groups=factor(I.cld$Groups, levels=c("Wild1","Wild2",
                                           "One",
                                           "Two"))
I.plot <- ggplot(data = I.cld, mapping = aes(x=Groups, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Cowpea Genepools", y= "Investment into symbiosis") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

I.plot


###############
# Angelas Data
###############


#Subset of data for Investment analysis
CowpeasInvestment <- subset(CowpeasAll, Investment > 0,
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, Investment))

CowpeasInvestment <- subset(CowpeasAll, Treatment == "S",
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, Investment))

#Model using interaction between treatment and line:
lmer.I2 <- lmer(Investment ~ DaysSinceInoculation + Groups + (1|Line), data = CowpeasInvestment)
Anova(lmer.I2)

#Lsmeans

lsmeans.I <- emmeans(lmer.I2, pairwise ~ Groups, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
I.cld <- cld(lsmeans.I$emmeans, Letters = letters)
I.cld$.group=gsub(" ", "", I.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
I.plot <- ggplot(data = I.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Investment") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
I.plot

I.cld$Groups=factor(I.cld$Groups, levels=c("Wild1","Wild2",
                                           "One",
                                           "Two"))

I.plot <- ggplot(data = I.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4", "royalblue4")) +
  labs(x="Inoculation Treatments", y= "Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

I.plot

##############################
# Host Growth Response
##############################

#Subset data for specific analysis
CowpeasHostGrowth <- subset(CowpeasAll, HGR > 0 | HGR < 0, 
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, HGR))


CowpeasHostGrowth <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L", 
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, HGR))

CowpeasHostGrowth$log.HGR <- log((CowpeasHostGrowth$HGR) + 1 - min(CowpeasHostGrowth$HGR))
qqPlot(CowpeasHostGrowth$log.HGR)


#Model using interaction between treatment and line:
lmer.HG <- lmer(log.HGR ~ DaysSinceInoculation + Treatment*Groups + (1|Line) + (1|Treatment:Line), data = CowpeasHostGrowth)
Anova(lmer.HG)

#Lsmeans
lsmeans(lmer.HG, pairwise ~ Treatment|Groups, adjust ="bonferroni")

lsmeans.HG <- emmeans(lmer.HG, pairwise ~ Treatment|Groups, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
HG.cld <- cld(lsmeans.HG$emmeans, Letters = letters)
HG.cld$.group=gsub(" ", "", HG.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 

non.plot <- ggplot(data = HG.cld, mapping = aes(x=Groups, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Host Growth Response (%)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot


HG.cld$Groups=factor(HG.cld$Groups, levels=c("Wild1","Wild2",
                                             "One",
                                             "Two"))
non.plot <- ggplot(data = HG.cld, mapping = aes(x=Groups, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Cowpea Genepools", y= "Investment into symbiosis") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

non.plot

########################################## Angela

#Subset data for specific analysis
CowpeasHostGrowth <- subset(CowpeasAll, HGR > 0 | HGR < 0, 
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, HGR))
CowpeasHostGrowth <- subset(CowpeasAll, Treatment == "S", 
                            select=c(Line, Treatment, Groups, DaysSinceInoculation, HGR))
CowpeasHostGrowth$log.HGR <- log((CowpeasHostGrowth$HGR) + 1 - min(CowpeasHostGrowth$HGR))
qqPlot(CowpeasHostGrowth$log.HGR)


#Model using interaction between treatment and line:
lmer.HG <- lmer(log.HGR ~ DaysSinceInoculation + Groups + (1|Line), data = CowpeasHostGrowth)
Anova(lmer.HG)

#Lsmeans
lsmeans(lmer.HG, pairwise ~ Groups, adjust ="bonferroni")

lsmeans.HG <- emmeans(lmer.HG, pairwise ~ Groups, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
HG.cld <- cld(lsmeans.HG$emmeans, Letters = letters)
HG.cld$.group=gsub(" ", "", HG.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 

non.plot <- ggplot(data = HG.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Host Growth Response (%)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
non.plot

HG.cld$Groups=factor(HG.cld$Groups, levels=c("Wild1","Wild2",
                                             "One",
                                             "Two"))

HG.plot <- ggplot(data = HG.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4", "royalblue4")) +
  labs(x="Inoculation Treatments", y= "Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

HG.plot
######################################################
# Dry Nodule Weight
#####################################################

CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0, 
                                             select=c(Line, Treatment, Groups, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Treatment == "A" | Treatment == "AL" | Treatment == "L", 
                                             select=c(Line, Treatment, Groups, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected <- CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected+1

CowpeasAllDryNoduleWeightCorrected$LogDry.Nodule.Weight.Corrected <-log10(CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected)


#New

lmer.DNB2 <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Treatment*Groups + (1|Line) + (1|Treatment:Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB2)

#Lsmeans
lsmeans(lmer.DNB2, pairwise ~ Groups*Treatment, adjust ="bonferroni")

lsmeans.DNB <- emmeans(lmer.DNB2, pairwise ~ Treatment|Groups, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
DNB.cld <- cld(lsmeans.DNB$emmeans, Letters = letters)
DNB.cld$.group=gsub(" ", "", DNB.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
DNB.plot <- ggplot(data = DNB.cld, mapping = aes(x=Groups, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Dry Nodule Biomass") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
DNB.plot

DNB.cld$Groups=factor(DNB.cld$Groups, levels=c("Wild1","Wild2",
                                               "One",
                                               "Two"))
DNB.plot <- ggplot(data = DNB.cld, mapping = aes(x=Groups, y=emmean, fill=Treatment)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_viridis_d(option = "mako") +
  labs(x="Cowpea Genepools", y= "Dry Nodule Biomass") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

DNB.plot
#############################
# Angelas
############################


CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Dry.Nodule.Weight.Corrected > 0, 
                                             select=c(Line, Treatment, Groups, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))
CowpeasAllDryNoduleWeightCorrected <- subset(CowpeasAll, Treatment == "S", 
                                             select=c(Line, Treatment, Groups, DaysSinceInoculation, Dry.Nodule.Weight.Corrected))

CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected <- CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected+1

CowpeasAllDryNoduleWeightCorrected$LogDry.Nodule.Weight.Corrected <-log10(CowpeasAllDryNoduleWeightCorrected$Dry.Nodule.Weight.Corrected)


#New

lmer.DNB2 <- lmer(LogDry.Nodule.Weight.Corrected ~ DaysSinceInoculation + Groups + (1|Line), data = CowpeasAllDryNoduleWeightCorrected)
Anova(lmer.DNB2)

#Lsmeans
lsmeans(lmer.DNB2, pairwise ~ Groups, adjust ="bonferroni")

lsmeans.DNB <- emmeans(lmer.DNB2, pairwise ~ Groups, adjust="Tukey") #Note that by including the type="response" we are back-transforming the data
DNB.cld <- cld(lsmeans.DNB$emmeans, Letters = letters)
DNB.cld$.group=gsub(" ", "", DNB.cld$.group)  ###  Remove spaces in .group 

#Plotting overall nodule number by species and treatment
pd = position_dodge(0.4) 
DNB.plot <- ggplot(data = DNB.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Groups", values = c("Wild1" = "#FDE725FF", "Wild2" = "#CC9900","One" = "#440154FF", "Two" = "#21908CFF")) +
  labs(x="Inoculation Treatments", y= "Log Dry Nodule Biomass") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        axis.title.x = element_text(face="bold", vjust=1.0, size = 14),
        axis.title.y = element_text(face="bold", vjust=1.0, size = 14),
        axis.text.x = element_text(size=16, colour ="black"),
        axis.text.y = element_text(size=16, colour ="black"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
  geom_text(aes(label = .group),
            position = position_dodge(0.6),
            vjust = -4,hjust=0.6,
            color="black",
            size=4)
DNB.plot

DNB.plot <- ggplot(data = DNB.cld, mapping = aes(x=Groups, y=emmean, fill=Groups)) + 
  geom_bar(position = position_dodge (), colour= "black", stat= "identity") +
  geom_errorbar(aes(ymin=emmean + SE, ymax=emmean - SE), width=0.3, colour="black", position=position_dodge(0.9)) +
  scale_fill_manual("Genepool", values = c("royalblue4","royalblue4","royalblue4", "royalblue4")) +
  labs(x="Inoculation Treatments", y= "Mean Nodule Weight") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  cleanup

DNB.plot



