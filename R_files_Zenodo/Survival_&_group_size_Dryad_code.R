#CODE FOR ARTICLE TITLED: "Survival benefits of group living in a
#fluctuating environment" BY GUINDRE-PARKER & RUBENSTEIN 
#American Naturalist

#LOAD PACKAGES
library(survival)
library(car)

#IMPORT DATA FILES
dataset_males<-read.csv("Survival Groupsize Males Dryad.csv")
dataset_females<-read.csv("Survival Groupsize Females Dryad.csv")




#SURVIVAL MODEL FOR MALES
attach(dataset_males)
surv.mod.male<-Surv(entry.age,exit.age,dead)
mod.male<-coxph(surv.mod.male~cluster(group.id)+cluster(id)+
                  (prerain.zscore+breedrain.zscore+grasscover.zscore)*(groupsize.zscore))
summary(mod.male)

#LOAD 95% CI 
confint(mod.male)

#CHECK PROPORTIONAL HAZARD ASSUMPTION
cox.zph(mod.male)

#CHECK VARIANCE INFLATION FACTORS (see FOX & MONETTE 1992)
vif(mod.male)
detach()




#SURVIVAL MODEL FOR FEMALES
attach(dataset_females)
surv.mod.female<-Surv(entry.age,exit.age,dead)
mod.female<-coxph(surv.mod.female~cluster(group.id)+cluster(id)+
                    (prerain.zscore+breedrain.zscore+grasscover.zscore)*(groupsize.zscore))
summary(mod.female)

#LOAD 95% CI 
confint(mod.female)

#CHECK PROPORTIONAL HAZARD ASSUMPTION
cox.zph(mod.female)

#CHECK VARIANCE INFLATION FACTORS (see FOX & MONETTE 1992)
vif(mod.female)
detach()
