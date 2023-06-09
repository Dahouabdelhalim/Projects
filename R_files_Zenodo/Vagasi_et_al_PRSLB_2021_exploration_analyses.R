####################################################################################
#------------------------------------ R code --------------------------------------- 
# ------------------------------------ for -----------------------------------------
#------------------------ Vágási et al. 2021 Proc.R.Soc.B --------------------------
#------------- "Social groups with diverse personalities mitigate ------------------ 
#--------------------- physiological stress in a songbird" -------------------------
#-----------------------------------------------------------------------------------
#-------------------------- REPEATABILITY OF EXPLORATION ---------------------------
#-----------------------------------------------------------------------------------
####################################################################################

# AIM OF THE ANALYSIS: 
# - to test the repeatability of exploratory behaviour of house sparrows

# GENERAL COMMENTS:
# - code just good to calculate repeatability of exploratory behaviour of house 
#   sparrows for the experimental study on groups personality composition and physiology
# - abbreviations are explained in details in the README file accompanying the data frame  
# - code written by Attila Fülöp: fafeldolgozo@gmail.com


####################################################################################
#------------------------------- Preparatory steps ---------------------------------
####################################################################################

#load libraries
library(lme4)
library(car)
library(rptR)


####################################################################################

#set the working directory
setwd("~/")

#load data frame with the exploration data
datp <- read.csv(file="Vagasi_et_al_PRSLB_2021_PASDOM_exploration_DATA.csv", header=TRUE, sep=",", dec=".", na.strings=c("NA",""," "))
head(datp)
summary(datp)
str(datp)

#set factors and other variables
datp$ring <- factor(datp$ring)
datp$sex <- factor(datp$sex)
datp$comp <- factor(datp$comp)
datp$run <- factor(datp$run)
datp$run_day <- as.numeric(as.character(datp$run_day))
datp$run_day_order <- as.numeric(as.character(datp$run_day_order))

#check summary of data
summary(datp)
str(datp) #OK



####################################################################################
#------------------------------------ ANALYSES -------------------------------------
####################################################################################

################################################
#test the effect of predictors in a LMM - full model with ALL predictors included
m <- lmer(scale(log(exploration + 1)) ~ (1|ring) + (1|run_day) + (1|run_day:run_day_order) + 
            sex + run + comp + sex:run + sex:comp + run:comp,
           data=datp, na.action=na.omit)
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")

m2 <- update(m, .~. -sex:run)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -sex:comp)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sex)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -run:comp)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -comp)
summary(m6)
drop1(m6, test="Chisq")
m7 <- update(m6, .~. -run)
summary(m7)
Anova(m7, type=2) 
drop1(m7, test="Chisq") #MAM reached


################################################
#calculate repeatability based on the full model using the "rptR" R package
set.seed(1234)
fm <- rptGaussian(scale(log(exploration + 1)) ~ (1|ring) + (1|run_day) + (1|run_day:run_day_order) + 
                   sex + run + comp + sex:run + sex:comp + run:comp, 
                 grname=c("ring"), data=datp, nboot=4999, npermut=4999, parallel=TRUE, ncores=2)
print(fm) #r=0.472, SE=0.098, 95%CI=0.305-0.687, P<0.0001, P.perm=0.0004
plot(fm)


#calculate repeatability based on the minimal model
set.seed(1234)
mm <- rptGaussian(scale(log(exploration + 1)) ~ (1|ring) + (1|run_day) + (1|run_day:run_day_order) + 1, 
                 grname=c("ring"), data=datp, nboot=4999, npermut=4999, parallel=TRUE, ncores=2)
print(mm) #r=0.416, SE=0.099, 95%CI=0.200-0.588, P<0.0001, P.perm=0.0002
plot(mm)

