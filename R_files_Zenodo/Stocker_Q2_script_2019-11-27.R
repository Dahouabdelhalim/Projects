
# Stocker et al. 2019: Cooperation with closely bonded individuals reduces cortisol levels in long-tailed macaques
# Q2: ARE CHANGES IN CORTISOL CONNECTED WITH THE INDIVIDUALS' COOPERATIVE SUCCESS AND/OR WITH OTHER SOCIAL FACTORS?
#--------------------------------------------------------------------


setwd("your_working_directory")

LTM <- read.csv("Stocker_Data_Nov2019.csv")
str(LTM)

# Q2a dataset with cooperation condition only
LTM.Coop <- LTM[!is.na(LTM$Success), ]
LTM.Coop <- LTM.Coop[!is.na(LTM.Coop$cort.1), ]
# ... and no missing delta cort
LTM.Coop.Q2 <- LTM.Coop[!is.na(LTM.Coop$cort.2), ]
LTM.Coop.Q2$Kin.maternal <- droplevels(LTM.Coop.Q2$Kin.maternal)


# Q2b dataset with social control condition only
LTM_noZ <- LTM[!is.na(LTM$delta.2.1), ]
LTM.SocCont <- LTM_noZ[LTM_noZ$Condition == "Control",]
LTM.SocCont$Kin.maternal <- droplevels(LTM.SocCont$Kin.maternal)


#--------------------------------------------------------------------
# full model Q2a

library(lme4)
summary(m.Q2.a <- lmer(delta.2.1 ~ Success + Soc.bond.prop + Sex.indiv + Rank.indiv + Baby + Kin.maternal +  
                         (1|Individual) , data=LTM.Coop.Q2, REML=FALSE))


library(car)
vif(m.Q2.a)    # no mulit-collinearity issue

#standardize model for averaging
library(arm)
m.Q2.a_stdz <- standardize(m.Q2.a)

# Model averaging
library(MuMIn)
options(na.action=na.fail)
all.m.Q2.a <- dredge(m.Q2.a_stdz, rank="AICc") #all model combinations
all.m.Q2.a
sub.m.Q2.a <- subset(all.m.Q2.a, delta<2) # model combinations less than 2 AIC difference to the 'best' models
sub.m.Q2.a
avg.m.Q2.a <- model.avg(sub.m.Q2.a)
summary(avg.m.Q2.a)
importance(avg.m.Q2.a)

confint(avg.m.Q2.a)



#--------------------------------------------------------------------
# full model Q2b

library(lme4)
summary(m.Q2.b <- lmer(delta.2.1 ~ Soc.bond.prop + Sex.indiv + Rank.indiv + Baby + Kin.maternal +
                            (1|Individual), data=LTM.SocCont, REML = F))

library(car)
vif(m.Q2.b) #no mulit-collinearity issue

library(arm)
m.Q2.b_stdz <- standardize(m.Q2.b)

# Model averaging
library(MuMIn)
options(na.action=na.fail)
all.m.Q2.b <- dredge(m.Q2.b_stdz, rank="AICc")
all.m.Q2.b
sub.m.Q2.b <- subset(all.m.Q2.b, delta<2)
sub.m.Q2.b                      # null model is the best!



