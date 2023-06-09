#CODE FOR THE ANALYSES FOR MANUSCRIPT:
#Plasticity in social behavior varies with reproductive status 
#in an avian cooperative breeder
#AUTHORS: LITTLE, RUBENSTEIN & GUINDRE-PARKER






#FLEDGING SUCCESS ANALYSES

#Bring in the dataset file
data_fledge<-read.csv("guarding rep successV2.csv")
#Load package needed
library(lme4)

#generalized linear mixed model exploring how nest guarding and environmental factors 
#shape the number of fledglings that successfully leave the nest 
#Corresponds to Table 1 in manuscript
fledge.mod<-glmer(fledglings~scale(total.guarding)+
                    scale(clutch.size)+
                    scale(breed.rain)+
                    scale(pre.breed.rain)+
                    season+
                    (1|social.group.ID),
                  data_fledge,na.action=na.omit,family=poisson)

#Call model output and confidence intervals
summary(fledge.mod)
confint(fledge.mod)










#GUARDING BEHAVIOR ANALYSES

#Bring in the dataset file
data_guard<-read.csv("plasticity superbV2.csv")
#Load package needed
library(lme4)
library(bbmle)


#Guarding behavior is skewed towards 0
#thus Guarding is log transformed for analyses below
hist(data_guard$guarding)
hist(log(data_guard$guarding))





#linear mixed model exploring which intrinsic and extrinsic factors 
#influence nest guarding behavior among individual superb starlings
#Corresponds to Table S1 in ESM
guard_mod<-lmer(log(guarding)~sex*role+
                   sex*relatedness+
                   scale(clutch.size)+
                   scale(age)+
                   scale(breed.rain)+
                   scale(pre.breed.rain)+
                   (res.caretakers)+
                   (1|social.group.ID)+
                   (1|bird.ID),
                 data_guard,na.action=na.omit,REML = F)

#Call model output and confidence intervals
summary(guard_mod)
confint(guard_mod)

#Check model residuals
hist(residuals(guard_mod))
plot(residuals(guard_mod))


#GUARDING RANDOM SLOPE AND INTERCEPT MODEL COMPARISON ANALYSES

#comparing multiple nested models that ranged from 
#no random effect of individual ID (model 1)
#to including random intercepts for each individual (model 2)
#and finally to including random intercepts and random slopes 
#allowing for individual starlings to adjust their nest guarding 
#according to breeding rainfall (model 3).


#BREEDER MODELS
#Model 1 includes only a random effect for social group ID
model1_breeders<-lmer(log(guarding)~sex+
                        scale(clutch.size)+
                        scale(breed.rain)+
                        (1|social.group.ID),
                      subset(data_guard, role=="Breeder"),
                      na.action=na.omit,REML = F)

#Model 2 includes random intercepts for individual bird ID
model2_breeders<-lmer(log(guarding)~sex+
                        scale(clutch.size)+
                        scale(breed.rain)+
                        (1|social.group.ID)+
                        (1|bird.ID),
                      subset(data_guard, role=="Breeder"),
                      na.action=na.omit,REML = F)

#Model 3 includes random intercepts and random slopes for individual bird ID
model3_breeders<-lmer(log(guarding)~sex+
                        scale(clutch.size)+
                        scale(breed.rain)+
                        (1|social.group.ID)+
                        (1|bird.ID)+
                        (0+scale(breed.rain)|bird.ID),
                      subset(data_guard, role=="Breeder"),
                      na.action=na.omit,REML = F)





#HELPER MODELS
#Model 1 includes only a random effect for social group ID
model1_helpers<-lmer(log(guarding)~sex+
                        scale(clutch.size)+
                        scale(breed.rain)+
                        (1|social.group.ID),
                      subset(data_guard, role=="Helper"),
                     na.action=na.omit,REML = F)

#Model 2 includes random intercepts for individual bird ID
model2_helpers<-lmer(log(guarding)~sex+
                        scale(clutch.size)+
                        scale(breed.rain)+
                        (1|social.group.ID)+
                        (1|bird.ID),
                      subset(data_guard, role=="Helper"),
                     na.action=na.omit,REML = F)

#Model 3 includes random intercepts and random slopes for individual bird ID
model3_helpers<-lmer(log(guarding)~sex+
                        scale(clutch.size)+
                        scale(breed.rain)+
                        (1|social.group.ID)+
                        (1|bird.ID)+
                        (0+scale(breed.rain)|bird.ID),
                      subset(data_guard, role=="Helper"),
                     na.action=na.omit,REML = F)



#We identified the best fitting model(s) using AICc, where all models within 2 Î”AICc 
#of the top model are interpreted to fit the dataset equally well
#Corresponds to Table 2 in the manuscript

#Breeder AIC model comparison
bbmle::AICctab(model1_breeders,model2_breeders,model3_breeders,base = T,weights = T)
#Helper AIC model comparison
bbmle::AICctab(model1_helpers,model2_helpers,model3_helpers,base = T,weights = T)







#results of top linear mixed models for BREEDERS
#Corresponds to Table 3 in the manuscript

#Breeder model 2
summary(model2_breeders)
confint(model2_breeders)
#Breeder model 3
summary(model3_breeders)
confint(model3_breeders)




#results of top linear mixed models for HELPERS
#Corresponds to Table 4 in the manuscript

#Helper model 1
summary(model1_helpers)
confint(model1_helpers)
#Helper model 2
summary(model2_helpers)
confint(model2_helpers)












#GUARDING REPEATABILITY ANALYSES

#Load package needed
library(rptR)

#we estimated the repeatability of nest guarding behavior 
#for breeders and for helpers separately
#Corresponds to the "Repeatability in nest guarding behavior" section
#in the RESULTS section of the manuscript

#Breeder repeatability
rep.breeder<-rpt(log(guarding)~sex+
                 scale(clutch.size)+
                 scale(breed.rain)+
                 (1|social.group.ID)+
                 (1|bird.ID), 
               grname = c("bird.ID","social.group.ID"), 
               data =  subset(data_guard, role=="Breeder"), 
               datatype = "Gaussian",  
               nboot = 1000, npermut = 0, adjusted = FALSE)

summary(rep.breeder)



#Helper repeatability
rep.helper<-rpt(log(guarding)~sex+
                   scale(clutch.size)+
                   scale(breed.rain)+
                   (1|social.group.ID)+
                   (1|bird.ID), 
                 grname = c("bird.ID","social.group.ID"), 
                 data =  subset(data_guard, role=="Helper"), 
                 datatype = "Gaussian",  
                 nboot = 1000, npermut = 0, adjusted = FALSE)

summary(rep.helper)




