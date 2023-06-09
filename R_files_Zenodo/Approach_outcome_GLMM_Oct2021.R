#Title: Factors determining whether or not an African Wild dog will successfully join a kill/carcass
#Script: Performs GLMMs on "approach" data and tests model assumptions
#Data: Collected from videos of kill sites

#load packages
library(readxl)
library(lme4)
library(MuMIn)
library(ggplot2)


#setwd to source file and read file "Approach outcome - PFQ_Oct2021.xlsx", sheet = "R ready")
#assuming R studio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
Approach_outcome_PFQ <- read_excel("Approach outcome - PFQ_Oct2021.xlsx", sheet = "R ready")

#format data for analysis
View(Approach_outcome_PFQ)
data<-(Approach_outcome_PFQ)
names(data)

# need to remove NA's in data
data<-data[!is.na(data$CarcassCondition), ]


#dependent variable
data$Outcome<-as.factor(data$Outcome)
#possible predictors
data$CarcassCondition<-as.numeric(data$CarcassCondition)
data$IDRelativetoMinPFQ<-as.numeric(data$IDRelativetoMinPFQ)
data$Noonkill_Corrected<-as.numeric(data$Noonkill_Corrected)


#Random terms
data$Encounter<-as.factor(data$Encounter)
data$ID<-as.factor(data$ID)



#check variables for correlation
cor(data$IDRelativetoMinPFQ, data$Noonkill_Corrected)
#not correlated
cor(data$IDRelativetoMinPFQ, data$CarcassCondition)
#not correlated
cor(data$Noonkill_Corrected, data$CarcassCondition)
#not correlated

summary(data)
View(data)

#data exploration
plot(data$Outcome, data$Noonkill_Corrected, xlab= "Outcome 0=Rejection 1= Join", ylab= "Number on kill (corrected)", main = "Number on kill ~ outcome")
plot(data$Outcome, data$IDRelativetoMinPFQ, xlab= "Outcome 0=Rejection 1= Join", ylab= "-ve = approacher has a higher PFQ, 0 = same, +ve approacher has a lower PFQ", main = "Relative to MinPFQ on kill ~ outcome")
plot(data$Outcome, data$CarcassCondition, xlab= "Outcome 0=Rejection 1= Join", ylab= "Carcass Condition; 100 = intact, 0 = skin and bones", main = "Carcass Condition ~ outcome")


#------------------------------------------------------------------------------------------------------------------------------
#Model it 
#Keep both Random terms in, some ID's are not always in every encounter so not nested. Encounters control for prey type, length of video etc.
#Ids, controls for individual variation. These are not points of interest but things to control for.

#Run full model. Models are somewhat similar, and need to control for number on kill and carcass condition
basic<-glmer(Outcome~1+
               (1|Encounter)+(1|ID),family=binomial,data=data)
AICc (basic) #196.4926


full<-glmer(Outcome~(IDRelativetoMinPFQ +CarcassCondition + Noonkill_Corrected)^2 +
             (1|Encounter)+(1|ID),family=binomial,data=data)
AICc(full) #197.0984

#checking overdispersion

install.packages ("blmeco")
library(blmeco) 
dispersion_glmer(full) #0.748 it shouldn't be over ~1.4. all good.

#dredge full model
mdfull<-dredge(full, rank = "AICc")
mdfull


#include competitive models in model averaging
subset(mdfull, delta < 2)
mafull<-model.avg(mdfull, delta<2, rank = "AICc")
summary(mafull)


##############################################################################################################

#The model averaging result seems robust. Noonkill has limited impact, carcass state is minimal, PFQ is key
#checking overdispersion in all candidate models

install.packages ("blmeco")
library(blmeco) 
dispersion_glmer(4) #0.748 it shouldn't be over ~1.4. all good.
m1<-glmer(formula = Outcome ~ CarcassCondition + IDRelativetoMinPFQ + (1 | Encounter) + (1 | ID), data = data, family = binomial)
m2<-glmer(formula = Outcome ~ IDRelativetoMinPFQ +
            (1 | Encounter) + (1 | ID), data = data, family = binomial)
m3<-glmer(formula = Outcome ~ IDRelativetoMinPFQ + Noonkill_Corrected +
            (1 | Encounter) + (1 | ID), data = data, family = binomial)
m4<-glmer(formula = Outcome ~ CarcassCondition + IDRelativetoMinPFQ + Noonkill_Corrected + CarcassCondition:Noonkill_Corrected + 
            (1 | Encounter) + (1 | ID), data = data, family = binomial)
m5<-glmer(formula = Outcome ~ CarcassCondition + IDRelativetoMinPFQ + Noonkill_Corrected +
            (1 | Encounter) + (1 | ID), data = data, family = binomial)
m6<-glmer(formula = Outcome ~ CarcassCondition + IDRelativetoMinPFQ + CarcassCondition:IDRelativetoMinPFQ +
            (1 | Encounter) + (1 | ID), data = data, family = binomial)

dispersion_glmer(m1) #0.7542921 (it shouldn't be over ~1.4. all good.)
dispersion_glmer(m2) #0.747579
dispersion_glmer(m3) #0.7477286
dispersion_glmer(m4) #0.7540024
dispersion_glmer(m5) #0.7520012
dispersion_glmer(m6) #0.7544154



