#set file direction 
setwd("<file.directory>")
d<-read.csv("MateChoice-UV-BSA.csv",header=T, sep=",")
#Remove trials when the female have not visit one of the three male
d<-d[d$TrialsToRemoveBecauseFemaleHaveNotSeenThe3males=="N",]

#Adding the boundary length of the male colour pattern
datatransition<-read.csv("transition.csv", header=T,sep=",")

index<-match(d$male,datatransition$male)
index
d<-cbind(d,datatransition[index,2:23])
summary(d)


#------------------------#
#    packagesloading     #
#------------------------#
library(boot)	           #
library(lattice)         #
library(nlme)            #        
library(gplots)          #
library(lme4)            #
library(MASS)            #
library(hglm)            #
library(MCMCglmm)        #
library(graphics)        #
library(car)             #
#to get p on lmer        #
library(lmerTest)        #
#to run plot(allEffects) #
library(effects)         #
library(mgcv)            #
#------------------------#

#######################
#checking data format #
#######################
summary(d)

#trial has to be factor
d$trial<-as.factor(d$trial)

#data in minutes
d$TotalTimefemaleinterested<-d$TotalTimefemaleinterested/60000
d$TotalTimefemalenotinterested<-d$TotalTimefemalenotinterested/60000

d$trial<-as.numeric(d$trial)
#trial 2 from 1 to ....
d$alltrials<-1
d$alltrials[d$groupfish=="B"]<-2
d$alltrials[d$groupfish=="c"]<-3
d$alltrials[d$groupfish=="D"]<-4
d$alltrials[d$groupfish=="E"]<-5
d$alltrials[d$groupfish=="F"]<-6
d$alltrials[d$groupfish=="G"]<-7
d$alltrials[d$groupfish=="H"]<-8
d$alltrials[d$groupfish=="I"]<-9
d$alltrials[d$groupfish=="J"]<-10
d$alltrials[d$groupfish=="K"]<-11
d$alltrials[d$groupfish=="L"]<-12
d$alltrials[d$groupfish=="M"]<-13
d$alltrials[d$groupfish=="N"]<-14
d$alltrials[d$groupfish=="O"]<-15
d$alltrials[d$groupfish=="P"]<-16
d$alltrials[d$groupfish=="Q"]<-17
d$alltrials[d$groupfish=="R"]<-18
d$alltrials[d$groupfish=="S"]<-19
d$alltrials[d$groupfish=="T"]<-20
d$alltrials<-as.numeric(d$alltrials)

# numbering trial order within day
d$trialOrderDay<-1
d$trialOrderDay[d$groupfish=="B"]<-2
d$trialOrderDay[d$groupfish=="c"]<-3
d$trialOrderDay[d$groupfish=="D"]<-4
d$trialOrderDay[d$groupfish=="E"]<-1
d$trialOrderDay[d$groupfish=="F"]<-2
d$trialOrderDay[d$groupfish=="G"]<-3
d$trialOrderDay[d$groupfish=="H"]<-4
d$trialOrderDay[d$groupfish=="I"]<-1
d$trialOrderDay[d$groupfish=="J"]<-2
d$trialOrderDay[d$groupfish=="K"]<-3
d$trialOrderDay[d$groupfish=="L"]<-4
d$trialOrderDay[d$groupfish=="M"]<-1
d$trialOrderDay[d$groupfish=="N"]<-2
d$trialOrderDay[d$groupfish=="O"]<-3
d$trialOrderDay[d$groupfish=="P"]<-4
d$trialOrderDay[d$groupfish=="Q"]<-1
d$trialOrderDay[d$groupfish=="R"]<-2
d$trialOrderDay[d$groupfish=="S"]<-3
d$trialOrderDay[d$groupfish=="T"]<-4
d$trialOrderDay<-as.numeric(d$trialOrderDay)

# day order
d$DayOrder<-1
d$DayOrder[d$date=="24/03/2017"]<-2
d$DayOrder[d$date=="26/03/2017"]<-3
d$DayOrder[d$date=="28/03/2017"]<-4
d$DayOrder[d$date=="30/03/2017"]<-5
d$DayOrder[d$date=="1/04/2017"]<-6
d$trialOrderDay<-as.numeric(d$trialOrderDay)


d$sidefromfemalepov2[d$sidefromfemalepov=="left"]<-"zleft"
d$sidefromfemalepov2[d$sidefromfemalepov=="middle"]<-"middle"
d$sidefromfemalepov2[d$sidefromfemalepov=="right"]<-"right"
d$sidefromfemalepov2<-as.factor(d$sidefromfemalepov2)


###############################################


###########################
#  TEST WITH BSA measures #
###########################
names(d)



#########################
###### MN Delta S  ######
#########################


#selection of random intercept and slope

# WITh Day ORDERED (important to not just use day but the days since beginning of the exp)_ WITH RANDOM INTERCEPT AND RANDOM SLOPE--> THE RANDOM INTERCEPT NEED TO BE ADDED TO THE MODEL


summary(m51<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+(Exp3JLumMN+Exp3JClrMN|female)+(1|DayOrder),data=d)) # aic= 1053.355 # model failed to converge
summary(m52<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+trialOrderDay+(Exp3JLumMN+Exp3JClrMN|female)+(trialOrderDay|DayOrder),data=d))# aic=1058.977 # boundary (singular) fit: see ?isSingular
summary(m53<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+(1|female)+(1|DayOrder),data=d))# aic=1043.519
summary(m54<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+trialOrderDay+(1|female)+(trialOrderDay|DayOrder),data=d))# aic= 1049.144
summary(m55<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+(1|female)+(1|trialOrderDay),data=d))# aic= 1087.439 # boundary (singular) fit: see ?isSingular
summary(m56<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+trialOrderDay+(trialOrderDay|female),data=d))# aic=1073.822
summary(m57<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+DayOrder+trialOrderDay+(DayOrder|female)+(trialOrderDay|DayOrder),data=d))#aic=1045.290
summary(m58<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+DayOrder+(DayOrder|female)+(DayOrder|trialOrderDay),data=d))#aic=1043.356 #boundary (singular) fit: see ?isSingular
summary(m59<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+DayOrder+(DayOrder|female),data=d))#aic=1037.377 # best model 
summary(m60<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+(1|female),data=d))#aic=1085.439

AIC(m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)
# run anova when dropping a random factor
anova (m51,m52,refit=FALSE)  # refit= FALSE because we compare model with different random effects. REML to compare different random
anova (m52,m53,refit=FALSE)
anova (m53,m54,refit=FALSE)
anova (m53,m55,refit=FALSE)
anova (m55,m56,refit=FALSE)
anova (m59,m57,refit=FALSE)
anova (m53,m60,refit=FALSE)
anova (m53,m59,refit=FALSE)
anova (m59,m58,refit=FALSE)
anova (m58,m60,refit=FALSE)
anova (m59,m60,refit=FALSE)

#m59 best model
#REML is a method for estimating variance components in models with random effects. 
#If all effects are fixed, then using REML makes no sense because the first thing REML does, computationally speaking, 
#is removing all fixed effects and evaluating remaining variance that belongs to random effects

#You can compare nested models that only differ in
#the random terms by using the REML likelihood or the ordinary likelihood. If
#you want to compare models that differ in fixed effects terms, then you must use ordinary likelihood. A few words about REML Gary W. Oehlert 2011


#  MODEL selection WITH CONTROL OF TIME (random):  selection on interaction only (not on single variables)
summary(m59<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+DayOrder+(DayOrder|female),data=d))#aic=1037.377 # best model when random intercept is added as fix variable
summary(m1<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+DayOrder+(DayOrder|female),data=d))
summary(m2<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JLumMN+DayOrder+(DayOrder|female),data=d))
summary(m3<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov+DayOrder+(DayOrder|female),data=d))
summary(m4<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov+DayOrder+(DayOrder|female),data=d))
summary(m5<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JClrMN+Exp3JLumMN+sidefromfemalepov+DayOrder+(DayOrder|female),data=d))
summary(m5<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JClrMN+Exp3JLumMN+sidefromfemalepov+DayOrder+(DayOrder|female),data=d,REML=FALSE))

#post hoc test#
library(multcomp)
summary(glht(m5, linfct = mcp(sidefromfemalepov = "Tukey")), test = adjusted("holm"))
#BEST MODEL MEAN DELTA S IS M5

AIC(m1,m2,m3,m4,m5)

#run anova for each fixed variable dropped
anova(m3,m4) 

plot(m5)
hist(resid(m5))
###################################
#DONE july 2020, this refit model to ML
anova(m59,m1) # m59= 1013.5, m1= 1017.4
anova(m1,m2)  #m2= 1015.5
anova(m2,m3)  #m3=1012.2
anova(m3,m4)  #m4=1010.9
anova(m4,m5)  #m5=1008.9 , m5 is still best model, anova show no prob to remove interaction
summary(m5<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JClrMN+Exp3JLumMN+sidefromfemalepov+DayOrder+(DayOrder|female),data=d,REML=FALSE))
allEffects(m5)
ranef(m5)


#testing the interaction between JlumMN and JClrMN
summary(m59<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+DayOrder+(DayOrder|female),data=d))#aic=1037.377 # best model when random intercept is added as fix variable
summary(m591<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrMN+PropMaleResponseOnFemaleInterest*Exp3JLumMN+sidefromfemalepov*Exp3JClrMN+sidefromfemalepov*Exp3JLumMN+DayOrder+Exp3JClrMN*Exp3JLumMN+(DayOrder|female),data=d))
summary(m5<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JClrMN+Exp3JLumMN+sidefromfemalepov+DayOrder+(DayOrder|female),data=d))
summary(m51<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JLumMN+Exp3JClrMN+sidefromfemalepov+DayOrder+Exp3JClrMN*Exp3JLumMN+(DayOrder|female),data=d))#aic=1021.966

anova(m591,m59) # not significant + no effect if interaction dropped
anova(m51,m5) # not significant + no effect if interaction dropped
####################################
plot(allEffects(m5))





#########################
###### cv Delta S  ######
#########################

# WITh Day ORDERED _ WITH RANDOM INTERCEPT AND RANDOM SLOPE--> THE RANDOM INTERCEPT NEED TO BE ADDED TO THE MODEL
summary(m51<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+(Exp3JLumCV+Exp3JClrCV|female)+(1|DayOrder),data=d)) # aic= 1027.817
summary(m52<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+trialOrderDay+(Exp3JLumCV+Exp3JClrCV|female)+(trialOrderDay|DayOrder),data=d))# aic=1033.532
summary(m53<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+(1|female)+(1|DayOrder),data=d))# aic=1019.413
summary(m54<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+trialOrderDay+(1|female)+(trialOrderDay|DayOrder),data=d))# aic= 1025.225
summary(m55<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+(1|female)+(1|trialOrderDay),data=d))# aic= 1061.657
summary(m56<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+trialOrderDay+(trialOrderDay|female),data=d))# aic=1048.411
summary(m57<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+trialOrderDay+(DayOrder|female)+(trialOrderDay|DayOrder),data=d))#aic=1020.599
summary(m58<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+(DayOrder|female)+(DayOrder|trialOrderDay),data=d))#aic=1018.527
summary(m59<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+(DayOrder|female),data=d))#aic=1012.527 # best model when random intercept is added as fix variable
AIC(m51,m52,m53,m54,m55,m56,m57,m58,m59)

anova(m58,m59,refit=FALSE)

plot(m59)

# BEST MODEL WITH CONTROL OF TIME:  selection on interaction not on single variables
summary(m1<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+(DayOrder|female),data=d))#
summary(m2<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+(DayOrder|female),data=d))#
summary(m3<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JLumCV+Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+(DayOrder|female),data=d))#
summary(m4<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JLumCV+Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+(DayOrder|female),data=d))#
summary(m5<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JLumCV+Exp3JClrCV+sidefromfemalepov+DayOrder+(DayOrder|female),data=d))#

##################################################"
#refitting model with ML doing anova
anova(m1,m2)# m1=1017.4
anova(m2,m3)# m2=1015.5
anova(m3,m4)# m3=1012.2
anova(m4,m5)# m4=1010.9
            # m5= 1008.9
#ANOVA show no problem to remove the inteaction-fixed-variable to reach m5, no effect of CV chromatic and achromatic on male attractiveness with m5 same than with m1 

plot(allEffects(m5))
plot(m5)
hist(resid(m5))


###MODEL USED IN THE PAPER  #####
summary(m5<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JLumCV+Exp3JClrCV+sidefromfemalepov+DayOrder+(DayOrder|female),data=d,REML=FALSE))#aic=1008.9
summary(m5<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JLumCV+Exp3JClrCV+sidefromfemalepov2+DayOrder+(DayOrder|female),data=d,REML=FALSE))#aic=1008.9
AIC(m5)
ranef(m5)
###########################################################################


#testing the interaction between JlumCV and JClrCV
summary(m6<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JClrCV+PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+Exp3JClrCV*Exp3JLumCV+(DayOrder|female),data=d,REML=FALSE))
summary(m61<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JLumCV+sidefromfemalepov*Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+Exp3JClrCV*Exp3JLumCV+(DayOrder|female),data=d))#
summary(m62<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Exp3JLumCV+Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+Exp3JClrCV*Exp3JLumCV+(DayOrder|female),data=d))#
summary(m63<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JLumCV+Exp3JClrCV+sidefromfemalepov*Exp3JLumCV+DayOrder+Exp3JClrCV*Exp3JLumCV+(DayOrder|female),data=d))#
summary(m64<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Exp3JLumCV+Exp3JClrCV+sidefromfemalepov+DayOrder+Exp3JClrCV*Exp3JLumCV+(DayOrder|female),data=d,REML=FALSE))#aic=1021.966
step(m6)
anova(m1,m6)
anova(m2,m61)
anova(m3,m62)
anova(m4,m63)
anova(m5,m64) # M5 better AIC and interaction non significant , so can be removed

################################################""


###############################
#TEST transition length effect# 
###############################

#selection random effect


# HERE WITh Day ORDERED _ WITH RANDOM INTERCEPT AND RANDOM SLOPE -->THE RANDOM INTERCEPT NEED TO BE ADDED TO THE MODEL
summary(m53<-lmer(log(TotalTimefemaleinterested)~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi+(1|female)+(1|DayOrder),data=d))# aic=1409.017
summary(m54<-lmer(log(TotalTimefemaleinterested)~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi+trialOrderDay+(1|female)+(trialOrderDay|DayOrder),data=d))# aic= 1409.050
summary(m55<-lmer(log(TotalTimefemaleinterested)~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi+(1|female)+(1|trialOrderDay),data=d))# aic= 1444.141
summary(m56<-lmer(log(TotalTimefemaleinterested)~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi+trialOrderDay+(trialOrderDay|female),data=d))# aic=1431.474
summary(m57<-lmer(log(TotalTimefemaleinterested)~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi+trialOrderDay+DayOrder+(DayOrder|female)+(trialOrderDay|DayOrder),data=d))#aic=1407.875
summary(m58<-lmer(log(TotalTimefemaleinterested)~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi+DayOrder+(DayOrder|female)+(DayOrder|trialOrderDay),data=d))#aic=1409.225
summary(m59<-lmer(log(TotalTimefemaleinterested)~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi+DayOrder+(DayOrder|female),data=d))#aic= 1404.660# best model 
AIC(m53,m54,m55,m56,m57,m58,m59)
anova(m57,m59,refit=FALSE) #
ranef(m59)
#######################################################################################################################################################################
# BEST MODEL m59 no selection on the model because we decided to run selection only on non-significant INTERACTION and here there is no interaction between variables #
#######################################################################################################################################################################
summary(m59<-lmer(log(TotalTimefemaleinterested)~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi+DayOrder+(DayOrder|female),data=d,REML=FALSE))#aic= 1404.660# best model 






#########################################################################
#loading of all transition (boundary lenght and delta S) on mean Delta S#
#########################################################################

###first we need to add the info about the chromatic and achromatic D.S of each transition - for each fish###

#combining transition of Chromatic Delta S data
datatransitionDeltaS<-read.csv("TransitionWithDeltaSinsteadOfLenght.csv", header=T,sep=",") # chromatic Delta S
index<-match(d$male,datatransitionDeltaS$male)
index
d<-cbind(d,datatransitionDeltaS[index,2:22])
summary(d)

#combining transition of Chromatic Delta S data
datatransitionDeltaS<-read.csv("TransitionWithDeltaS-Achrom-insteadOfLenght.csv", header=T,sep=",") # achromatic Delta S
index<-match(d$male,datatransitionDeltaS$male)
index
d<-cbind(d,datatransitionDeltaS[index,2:22])
summary(d)
head(d)



###then we reduce the data set at one male per line###
data.active<-d[!duplicated(d$male),]# we select only 1 male of each because measures are repeted between males

summary(data.active)
head(data.active)
summary(data.active$bl.FzBl)

# we should not have random effect because we don't repeat males and just test how mich weight each of the transition as on the Mean Delta S values
# result part 3 :  factor loading for all transition on the mean delta S #


hist(data.active$Exp3JClrMN)


#----------------#
#effect of lenght
#----------------#
summary(m1<-lm(Exp3JClrMN~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi,data=data.active))# 

#--------------------------------------#
# effect of delta S of each transition
#--------------------------------------#

summary(m2<-lm(Exp3JClrMN~body.bl.DS+body.FzBl.DS+body.or.DS+body.gr.DS+body.sil.DS+body.vi.DS+bl.FzBl.DS+bl.or.DS+bl.gr.DS+bl.sil.DS+bl.vi.DS+FzBl.or.DS+FzBl.gr.DS+FzBl.sil.DS+FzBl.vi.DS+or.gr.DS+or.sil.DS+or.vi.DS+gr.sil.DS+gr.vi.DS+sil.vi.DS,data=data.active))# 
data.active$body.FzBl.DS # all fish possess the transition body.FzBl, so body.FzBl.DS is the same for the 63 fish, this is why =NA ->removed
summary(m2<-lm(Exp3JClrMN~body.bl.DS+body.or.DS+body.gr.DS+body.sil.DS+body.vi.DS+bl.FzBl.DS+bl.or.DS+bl.gr.DS+bl.sil.DS+bl.vi.DS+FzBl.or.DS+FzBl.gr.DS+FzBl.sil.DS+FzBl.vi.DS+or.gr.DS+or.sil.DS+or.vi.DS+gr.sil.DS+gr.vi.DS+sil.vi.DS,data=data.active))# aic=



#########################
#achromatic mean Delta S
#########################
hist(data.active$Exp3JLumMN)

#----------------#
#effect of lenght
#----------------#
summary(m1<-lm(Exp3JLumMN~body.bl+body.FzBl+body.or+body.gr+body.sil+body.vi+bl.FzBl+bl.or+bl.gr+bl.sil+bl.vi+FzBl.or+FzBl.gr+FzBl.sil+FzBl.vi+or.gr+or.sil+or.vi+gr.sil+gr.vi+sil.vi,data=data.active))# aic=1409.017

#--------------------------------------#
# effect of delta S of each transition
#--------------------------------------#
summary(m2<-lm(Exp3JLumMN~body.bl.ADS+body.FzBl.ADS+body.or.ADS+body.gr.ADS+body.sil.ADS+body.vi.ADS+bl.FzBl.ADS+bl.or.ADS+bl.gr.ADS+bl.sil.ADS+bl.vi.ADS+FzBl.or.ADS+FzBl.gr.ADS+FzBl.sil.ADS+FzBl.vi.ADS+or.gr.ADS+or.sil.ADS+or.vi.ADS+gr.sil.ADS+gr.vi.ADS+sil.vi.ADS,data=data.active))# aic=
data.active$body.FzBl.DS # all fish possess the transition body.FzBl, so body.FzBl.DS is the same for the 63 fish ->removed
summary(m2<-lm(Exp3JLumMN~body.bl.ADS+body.or.ADS+body.gr.ADS+body.sil.ADS+body.vi.ADS+bl.FzBl.ADS+bl.or.ADS+bl.gr.ADS+bl.sil.ADS+bl.vi.ADS+FzBl.or.ADS+FzBl.gr.ADS+FzBl.sil.ADS+FzBl.vi.ADS+or.gr.ADS+or.sil.ADS+or.vi.ADS+gr.sil.ADS+gr.vi.ADS+sil.vi.ADS,data=data.active))# aic=



####################################################################
####################################################################
#                test for supplementary material                   #
####################################################################
####################################################################

# variation in BSA measure within each group of males
tapply(d$Exp3JLumMN,d$groupfish,var)

summary(m<-lmer(Exp3JLumMN~groupfish+(1|female),data=d)) 

# time spend with male along the experiment
summary(m<-lmer(log(TotalTimefemaleinterested)~alltrials+(1|female),data=d)) 
plot(allEffects(m))
summary(m<-lmer(log(TotalTimefemaleinterested)~groupfish+(1|female),data=d)) 
#post hoc test#
library(multcomp)
summary(glht(m, linfct = mcp(groupfish = "Tukey")), test = adjusted("holm"))
plot(allEffects(m))

summary(m<-lmer(log(TotalTimefemaleinterested)~date+(1|female),data=d)) 
#post hoc test#
library(multcomp)
summary(glht(m, linfct = mcp(date = "Tukey")), test = adjusted("holm"))
plot(allEffects(m))

summary(m<-lmer(log(TotalTimefemaleinterested)~DayOrder+(1|female),data=d)) 
plot(allEffects(m))


summary(m<-lm(log(TotalTimefemaleinterested)~date*female,data=d)) 
#post hoc test#
library(multcomp)
summary(glht(m, linfct = mcp(date = "Tukey")), test = adjusted("holm"))
plot(allEffects(m))



####################################################################
####################################################################
#                test OPC effect on attractiveness                 #
####################################################################
####################################################################



#################### MEAN ##########################"

#Model selection on the random factor
# WITh Day ORDERED _ WITH RANDOM INTERCEPT AND RANDOM SLOPE --> THE RANDOM INTERCEPT NEED TO BE ADDED TO THE MODEL
summary(m51<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+(Mlum+Mchr|female)+(1|DayOrder),data=d)) # aic= 
summary(m52<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+trialOrderDay+(Mlum+Mchr|female)+(trialOrderDay|DayOrder),data=d))# aic=
summary(m53<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+(1|female)+(1|DayOrder),data=d))# aic=
summary(m54<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+trialOrderDay+(1|female)+(trialOrderDay|DayOrder),data=d))# aic= 
summary(m55<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+(1|female)+(1|trialOrderDay),data=d))# aic= 
summary(m56<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+trialOrderDay+(trialOrderDay|female),data=d))# aic=
summary(m57<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+DayOrder+trialOrderDay+(DayOrder|female)+(trialOrderDay|DayOrder),data=d))#aic=
summary(m58<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+DayOrder+(DayOrder|female)+(DayOrder|trialOrderDay),data=d))#aic
summary(m59<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+DayOrder+(DayOrder|female),data=d))#aic=
summary(m60<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+(1|female),data=d))#aic=

AIC(m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)

#   df       AIC
#m51 20 1002.5217
#m52 23 1009.8248
#m53 15 1001.0080
#m54 18 1008.3153
#m55 15 1037.5368
#m56 17 1023.7441
#m57 21 1005.2324
#m58 20 1001.3949
#m59 17  995.3949
#m60 14 1035.537

anova(m57,m58,refit=FALSE)
anova(m59,m58,refit=FALSE)

#BEST model m59

#NOW selection on fixed variable
summary(m59<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+DayOrder+(DayOrder|female),data=d))#aic=
summary(m1<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Mchr+PropMaleResponseOnFemaleInterest*Mlum+sidefromfemalepov*Mchr+sidefromfemalepov*Mlum+DayOrder+(DayOrder|female),data=d))#aic=
summary(m2<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Mchr+Mlum+sidefromfemalepov2*Mchr+sidefromfemalepov2*Mlum+DayOrder+(DayOrder|female),data=d))#aic=
summary(m3<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Mchr+Mlum+sidefromfemalepov+Mchr+sidefromfemalepov*Mlum+DayOrder+(DayOrder|female),data=d))#aic=
summary(m4<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Mchr+Mlum+sidefromfemalepov2+Mchr+Mlum+DayOrder+(DayOrder|female),data=d))#aic=

#Need to be fitted with ML for comparison of model with different fixed variables, anova does it automatically
anova(m59,m1) #m59=1010.5, m1=1008.5
anova(m1,m2)  #m2=1005.2
anova(m2,m3)  #m3=1005.2
anova(m3,m4)  #m4=1003.8

######MODEL PRESENTED IN PAPER##############
summary(m4<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Mchr+Mlum+sidefromfemalepov2+Mchr+Mlum+DayOrder+(DayOrder|female),data=d, REML=FALSE))#aic=
summary(m4<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+Mchr+Mlum+sidefromfemalepov+Mchr+Mlum+DayOrder+(DayOrder|female),data=d, REML=FALSE))#aic=

#Therefore best model m4 when fitted with ML, and no effect of Mchr and Mlum (when REML, m59, no effect of Mchr and Mlum)


#################### CV ##########################

#Model selection on the random factor
# WITh Day ORDERED _ WITH RANDOM INTERCEPT AND RANDOM SLOPE --> THE RANDOM INTERCEPT NEED TO BE ADDED TO THE MODEL
summary(m51<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+(CVlum+CVchr|female)+(1|DayOrder),data=d)) # aic= 
summary(m52<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+trialOrderDay+(CVlum+CVchr|female)+(trialOrderDay|DayOrder),data=d))# aic=
summary(m53<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+(1|female)+(1|DayOrder),data=d))# aic=
summary(m54<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+trialOrderDay+(1|female)+(trialOrderDay|DayOrder),data=d))# aic= 
summary(m55<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+(1|female)+(1|trialOrderDay),data=d))# aic= 
summary(m56<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+trialOrderDay+(trialOrderDay|female),data=d))# aic=
summary(m57<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+DayOrder+trialOrderDay+(DayOrder|female)+(trialOrderDay|DayOrder),data=d))#aic=
summary(m58<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+DayOrder+(DayOrder|female)+(DayOrder|trialOrderDay),data=d))#aic
summary(m59<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+DayOrder+(DayOrder|female),data=d))#aic=
summary(m60<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+(1|female),data=d))#aic=

AIC(m51,m52,m53,m54,m55,m56,m57,m58,m59,m60)
#    df      AIC
#m51 20 1023.922
#m52 23 1030.168
#m53 15 1012.784
#m54 18 1020.426
#m55 15 1052.152
#m56 17 1038.000
#m57 21 1018.056
#m58 20 1016.085
#m59 17 1010.085
#m60 14 1050.152

anova(m57,m58,refit=FALSE)
anova(m59,m58,refit=FALSE)

#BEST model m59

#NOW selection on fixed variable
summary(m59<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+PropMaleResponseOnFemaleInterest*CVhue+sidefromfemalepov*CVchr+sidefromfemalepov*CVlum+sidefromfemalepov*CVhue+DayOrder+(DayOrder|female),data=d))#aic=
summary(m1<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+PropMaleResponseOnFemaleInterest*CVhue+sidefromfemalepov+CVchr+sidefromfemalepov*CVlum+sidefromfemalepov*CVhue+DayOrder+(DayOrder|female),data=d))#aic=
summary(m2<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+PropMaleResponseOnFemaleInterest*CVhue+sidefromfemalepov+CVchr+CVlum+sidefromfemalepov*CVhue+DayOrder+(DayOrder|female),data=d))#aic=
summary(m3<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVlum+PropMaleResponseOnFemaleInterest*CVhue+sidefromfemalepov+CVchr+CVlum+CVhue+DayOrder+(DayOrder|female),data=d))#aic=
summary(m4<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+PropMaleResponseOnFemaleInterest*CVchr+PropMaleResponseOnFemaleInterest*CVhue+sidefromfemalepov+CVchr+CVlum+CVhue+DayOrder+(DayOrder|female),data=d))#aic=
summary(m5<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+PropMaleResponseOnFemaleInterest*CVchr+sidefromfemalepov+CVchr+CVlum+CVhue+DayOrder+(DayOrder|female),data=d))#aic=
summary(m6<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+sidefromfemalepov+CVchr+CVlum+CVhue+DayOrder+(DayOrder|female),data=d))#aic=

#comparison refit ML (anova)
anova(m59,m1) #m59=1099.0, m1=1008.2
anova(m1,m2)  #m2=1004.3
anova(m2,m3)  #m3=1001.6
anova(m3,m4)  #m4=1000.8
anova(m4,m5)  #m5=998.97, best model m5
anova(m5,m6)  #m6=999.04
#Therefore best model m5, and effect of CVchr (when REML, m59, no effect of CVchr and CVlum)
##############MODEL USED IN PAPER, ML fitted#####################
summary(m5bis<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+PropMaleResponseOnFemaleInterest*CVchr+sidefromfemalepov+CVchr+CVlum+CVhue+DayOrder+(DayOrder|female),data=d,REML=FALSE))#aic=
summary(m5bis<-lmer(log(TotalTimefemaleinterested)~PropMaleResponseOnFemaleInterest+PropMaleResponseOnFemaleInterest*CVchr+sidefromfemalepov2+CVchr+CVlum+CVhue+DayOrder+(DayOrder|female),data=d,REML=FALSE))#aic=

#m5bis
#                                         Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)                             -0.90406    0.36455  38.95789  -2.480  0.01757 *  
#PropMaleResponseOnFemaleInterest         1.89922    0.40044 448.81808   4.743 2.84e-06 ***
#CVchr                                    0.93714    0.33859 448.86241   2.768  0.00588 ** 
#sidefromfemalepovmiddle                 -0.39937    0.07997 446.62266  -4.994 8.50e-07 ***
#sidefromfemalepovright                  -0.17336    0.07813 446.38285  -2.219  0.02701 *  
#CVlum                                    1.23535    0.93636 446.51554   1.319  0.18774    
#CVhue                                   -0.62210    0.70020 448.52526  -0.888  0.37477    
#DayOrder                                -0.07492    0.06835   9.45708  -1.096  0.30013    
#PropMaleResponseOnFemaleInterest:CVchr  -1.70720    1.18471 448.56283  -1.441  0.15027   




#########################################
## CIRCULAR ANOVA FOR OPC ANALYSIS HUE ##
#########################################
library(circular)
require(CircStats)


###For the effect of hue on female choice###
d$huecircular<-circular(d$Mhue, type = c("angles"), units = c("degrees"), zero=pi, rotation="clock")
summary(m1<-lm.circular((d$TotalTimefemaleinterested),d$huecircular))
m1$rho
m1$A.k
m1$kappa
m1$p.values
plot(m1$residuals) # GOOD residual distribution if we DON'T try log d$TotalTimefemale interested
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##








#Zuur et al. (2009; PAGE 122) suggest that "To compare models with nested fixed effects (but with the same random structure), 
#ML estimation must be used and not REML." This indicates to me that I ought to use ML since my random effects are the same in both models,
# but my fixed effects differ. [Zuur et al. 2009. Mixed Effect Models and Extensions in Ecology with R. Springer.]

#REML should not be used when comparing models with different fixed effects. 
#REML, however, often estimates the random effects parameters better and therefore 
#it is sometimes recommended to use ML for comparisons and REML for estimating a single (perhaps final) model.