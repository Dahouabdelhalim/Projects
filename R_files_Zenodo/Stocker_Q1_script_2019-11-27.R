
# Stocker et al. 2019: Cooperation with closely bonded individuals reduces cortisol levels in long-tailed macaques
# Q1: IS AN INDIVIDUAL'S COOPERATIVE BEHAVIOUR MODULATED BY ITS CORTISOL LEVEL BEFORE COOPERATION, SOCIAL BOND AND/OR RANK?
#--------------------------------------------------------------------


setwd("your_working_directory")

LTM <- read.csv("Stocker_Data_Nov2019.csv")
str(LTM)

# Dataset with cooperation condition only
LTM.Coop <- LTM[!is.na(LTM$Success), ]
LTM.Coop <- LTM.Coop[!is.na(LTM.Coop$cort.1), ]


#---------------------
# check distribution of cooperative success (number of trials)

hist(LTM.Coop$Success, 10)

# check beta distribution
suc_prop<-LTM.Coop$Success/15
LTM.Coop$suc_prop <- suc_prop
LTM.Coop$suc_prop[LTM.Coop$suc_prop==1] <- 0.9999  #cannot be 0 or 1 with beta distribution

library(fitdistrplus)
fit.beta<-fitdist(LTM.Coop$suc_prop, "beta")  
plot(fit.beta)  # better fit than other distributions



#---------------------
# full model

library(glmmADMB)
summary(m.Q1 <- glmmadmb (suc_prop ~ cort.1 + Soc.bond.prop + Rank.indiv + 
                            (1|Individual) + (1|Partner) , data = LTM.Coop, family = "beta"))


library(car)
vif(m.Q1)  #no mulit-collinearity issue


# Model averaging
library(MuMIn)
options(na.action=na.fail)
all.m.Q1 <- dredge(m.Q1, rank="AICc") #all model combinations
all.m.Q1
sub.m.Q1 <- subset(all.m.Q1, delta<2) # model combinations less than 2 AIC difference to the 'best' models
sub.m.Q1 # subset includes null model
#avg.m.Q1 <- model.avg(sub.m.Q1)

