### Interspecific Foraging Behavior Models ####

library(lme4)
library(lmerTest)
library(sjPlot)

##### DATA ####
data.all<-read.csv("Melin_et_al_Interspecific_Foraging_Behavior_data.csv",header=T,stringsAsFactors = T, fileEncoding="UTF-8-BOM")

#### rename primate species to compare against Ateles ####
data.all$Monkey2[data.all$Species=="Howler"] <- "2_Howler"
data.all$Monkey2[data.all$Species=="capuchin"] <- "3_Capuchin"
data.all$Monkey2[data.all$Species=="Ateles"] <- "1_Spider"
data.all$Monkey2<-as.factor(data.all$Monkey2)


##### MODELS #####

####Question 1 Foraging Behavior: Trichromats ONLY #####

trichromats<-subset(data.all,ColorVisionType=="Trichromat")


##sniffing
smell.m1<-glmer(Total.Smell~Monkey2+ offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=trichromats,family="poisson")
smell.m2<-glmer(Total.Smell~1 + offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=trichromats,family="poisson")
anova(smell.m1,smell.m2)
smell.emm.s <- emmeans(smell.m1, "Monkey2")
pairs(smell.emm.s)

##touch
touch.m1<-glmer(Total.Touch~Monkey2+offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=trichromats,family="poisson")
touch.m2<-glmer(Total.Touch~1+offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=trichromats,family="poisson")
anova(touch.m1,touch.m2)
touch.emm.s <- emmeans(touch.m1, "Monkey2")
pairs(touch.emm.s)


##bite/reject
bite.m1<-glmer(Bite.Reject~Monkey2+offset(log(Total))+(1|Animal2)+
                 (1|ScientificName2),data=trichromats,family="poisson")
bite.m2<-glmer(Bite.Reject~1+offset(log(Total))+(1|Animal2)+
                 (1|ScientificName2),data=trichromats,family="poisson")
anova(bite.m1,bite.m2)
bite.emm.s <- emmeans(bite.m1, "Monkey2")
pairs(bite.emm.s)


###### IRRs ######
exp(fixef(accept.m1)) ##incidence rate ratios of fixed effects
exp(confint(smell.m1,method="Wald")) ## confidence interval


####Question 2: Colour Vision Phenotype and Foraging Behavior ######

## Code plant as conspicuous or cryptic 
data.all$Plant.Color[data.all$ScientificName2=="Spondias mombin"]<-"conspicuous"
data.all$Plant.Color[data.all$ScientificName2=="Ficus_red"]<-"conspicuous"
data.all$Plant.Color[data.all$ScientificName2=="Bursera simaruba"]<-"conspicuous"
data.all$Plant.Color[data.all$ScientificName2=="Sciadodendron excelsum"]<-"conspicuous"
data.all$Plant.Color[data.all$ScientificName2=="Simarouba glauca"]<-"conspicuous"
data.all$Plant.Color[data.all$ScientificName2=="Ficus_green"]<-"cryptic"
data.all$Plant.Color<-as.factor(data.all$Plant.Color)

## change reference colour vision phenotype to Trichromat because of unequal sampling (see Methods)
data.all$CT[data.all$ColorVisionType=="Dichromat"] <- "2_Dichromat"
data.all$CT[data.all$ColorVisionType=="Trichromat"] <- "1_Trichromat"
data.all$CT<-as.factor(data.all$CT)

##GLMM Models
#sniffing
color.smell1<-glmer(Total.Smell~CT*Plant.Color+Monkey2+ offset(log(Total))+(1|Animal2)
                    +(1|ScientificName2),data=data.all,family="poisson")
color.smell2<-glmer(Total.Smell~Monkey2+Plant.Color+ offset(log(Total))+(1|Animal2)
                    +(1|ScientificName2),data=data.all,family="poisson")
anova(color.smell1,color.smell2)
color.smell.emm.s <- emmeans(color.smell1, "CT", "Plant.Color")
pairs(color.smell.emm.s, reverse=T)

#manual touch
color.touch1<-glmer(Total.Touch~CT*Plant.Color+ Monkey2 + offset(log(Total))+(1|Animal2)
                    +(1|ScientificName2),data=data.all,family="poisson")
color.touch2<-glmer(Total.Touch~Monkey2+ Plant.Color+offset(log(Total))+(1|Animal2)
                    +(1|ScientificName2),data=data.all,family="poisson")
anova(color.touch1,color.touch2)
color.touch.emm.s <- emmeans(color.touch1, "CT", "Plant.Color")
pairs(color.touch.emm.s, reverse=T)

#bite/reject
color.bite1<-glmer(Bite.Reject~CT*Plant.Color+Monkey2+ offset(log(Total))+(1|Animal2)
                   +(1|ScientificName2),data=data.all,family="poisson")
color.bite2<-glmer(Bite.Reject~Monkey2+Plant.Color+ offset(log(Total))+(1|Animal2)
                   +(1|ScientificName2),data=data.all,family="poisson")
anova(color.bite1,color.bite2)
color.bite.emm.s <- emmeans(color.bite1, "CT", "Plant.Color")
pairs(color.bite.emm.s, reverse=T)

#acceptance
color.accept1<-glmer(Total.Eat~CT*Plant.Color+Monkey2+ offset(log(Total))+(1|Animal2)
                     +(1|ScientificName2),data=data.all,family="poisson")
color.accept2<-glmer(Total.Eat~Monkey2+Plant.Color+ offset(log(Total))+(1|Animal2)
                     +(1|ScientificName2),data=data.all,family="poisson")
anova(color.accept1,color.accept2)
color.accept.emm.s <- emmeans(color.accept1, "CT", "Plant.Color")
pairs(color.accept.emm.s, reverse=T)


#### Supplemental Analysis 1: Alouatta Between Sites and Sexes #####
#Goal: test for effects of Site (SSR in Costa Rica, Isla Agaltepec in Mexico) on sensory behaviours
#Goal 2: test for effects of Sex on sensory behaviours in Alouatta

alouatta<-subset(data.all,Species=="Howler")
alouatta2<-subset(alouatta,Sex!="Unknown")

###For examining site effects, use: data=alouatta, main effect = Site
###For examining effects of Sex, repeat with data=alouatta2, main effect = Sex

##sniffing
smell.m1<-glmer(Total.Smell~Sex+ offset(log(Total))+
                  (1|ScientificName2),data=alouatta2,family="poisson")
smell.m2<-glmer(Total.Smell~1 + offset(log(Total))+
                  (1|ScientificName2),data=alouatta2,family="poisson")
anova(smell.m1,smell.m2)


##touch
touch.m1<-glmer(Total.Touch~Sex +offset(log(Total))+(1|Animal2)
                  ,data=alouatta2,family="poisson")
touch.m2<-glmer(Total.Touch~1+offset(log(Total))+(1|Animal2)
                  ,data=alouatta2,family="poisson")
anova(touch.m1,touch.m2)


##bite/reject
bite.m1<-glmer(Bite.Reject~Sex +offset(log(Total))+(1|Animal2)+
                 (1|ScientificName2),data=alouatta2,family="poisson")
bite.m2<-glmer(Bite.Reject~1+offset(log(Total))+(1|Animal2)+
                 (1|ScientificName2),data=alouatta2,family="poisson")
anova(bite.m1,bite.m2)


#### Alouatta subset with trichromat data (Question 1) ####
## Goal: Test how each Alouatta subset (SSR, Mexico) affects interspecific results

trichromats.SSR<-subset(trichromats,Site=="SSR")
trichromats.Mexico<-subset(trichromats,Site=="Mexico")


##touch ## repeat with data = trichromats.SSR
touch.m1<-glmer(Total.Touch~Monkey2+offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=trichromats.Mexico,family="poisson")
touch.m2<-glmer(Total.Touch~1+offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=trichromats.Mexico,family="poisson")
touch.m3<-glmer(Total.Touch~Monkey3+offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=trichromats.Mexico,family="poisson")
anova(touch.m1,touch.m2)

##bite/reject ## repeat with data = trichromats.SSR
bite.m1<-glmer(Bite.Reject~Monkey2+offset(log(Total))+(1|Animal2)+
                 (1|ScientificName2),data=trichromats.Mexico,family="poisson")
bite.m2<-glmer(Bite.Reject~1+offset(log(Total))+(1|Animal2)+
                 (1|ScientificName2),data=trichromats.Mexico,family="poisson")
bite.m3<-glmer(Bite.Reject~Monkey3+offset(log(Total))+(1|Animal2)+
                 (1|ScientificName2),data=trichromats.Mexico,family="poisson")
anova(bite.m1,bite.m2)


#### Supplemental Analysis 2: Polymorphic Females #####
## Goal: test Question 2 (effect of colour vision phenotype) on sensory behaviours 
##using only polymorphic female dataset

females<-subset(data.all,Sex=="Female")
females<-subset(females,Monkey2!="2_Howler")

##GLMM Models
#sniffing
color.smell3<-glmer(Total.Smell~CT*Plant.Color+Monkey2+ offset(log(Total))+(1|Animal2)
                    +(1|ScientificName2),data=females,family="poisson")
color.smell2<-glmer(Total.Smell~Monkey2+Plant.Color+ offset(log(Total))+(1|Animal2)
                    +(1|ScientificName2),data=females,family="poisson")
anova(color.smell1,color.smell2)
color.smell.emm.s <- emmeans(color.smell3, "CT", "Plant.Color")
pairs(color.smell.emm.s, reverse=T)

#manual touch
color.touch3<-glmer(Total.Touch~CT*Plant.Color+ Monkey2 + offset(log(Total))+(1|Animal2)
                    +(1|ScientificName2),data=females,family="poisson")
color.touch2<-glmer(Total.Touch~Monkey2+ Plant.Color+offset(log(Total))+(1|Animal2)
                    +(1|ScientificName2),data=females,family="poisson")
anova(color.touch3,color.touch2)
color.touch.emm.s <- emmeans(color.touch1, "CT", "Plant.Color")
pairs(color.touch.emm.s, reverse=T)

#bite/reject
color.bite3<-glmer(Bite.Reject~CT*Plant.Color+Monkey2+ offset(log(Total))+(1|Animal2)
                   +(1|ScientificName2),data=females,family="poisson")
color.bite2<-glmer(Bite.Reject~Monkey2+Plant.Color+ offset(log(Total))+(1|Animal2)
                   +(1|ScientificName2),data=females,family="poisson")
anova(color.bite1,color.bite2)
color.bite.emm.s <- emmeans(color.bite1, "CT", "Plant.Color")
pairs(color.bite.emm.s,reverse=T)

#acceptance
color.accept3<-glmer(Total.Eat~CT*Plant.Color+Monkey2+ offset(log(Total))+(1|Animal2)
                     +(1|ScientificName2),data=females,family="poisson")
color.accept2<-glmer(Total.Eat~Monkey2+Plant.Color+ offset(log(Total))+(1|Animal2)
                     +(1|ScientificName2),data=females,family="poisson")
anova(color.accept3,color.accept2)
color.accept.emm.s <- emmeans(color.accept3, "CT", "Plant.Color")
pairs(color.accept.emm.s,reverse=T)


#### Supplementary Analysis SEX
data.all2<-subset(data.all,Sex!="Unknown")
sex.smell.m1<-glmer(Total.Smell~Sex+CT+Monkey2+offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=data.all2,family="poisson")
sex.smell.m2<-glmer(Total.Smell~CT+Monkey2+offset(log(Total))+(1|Animal2)+
                  (1|ScientificName2),data=data.all2,family="poisson")
anova(sex.smell.m1,sex.smell.m2)

sex.touch.m1<-glmer(Total.Touch~Sex+CT+Monkey2+offset(log(Total))+(1|Animal2)+
                      (1|ScientificName2),data=data.all2,family="poisson")
sex.touch.m2<-glmer(Total.Touch~CT+Monkey2+offset(log(Total))+(1|Animal2)+
                      (1|ScientificName2),data=data.all2,family="poisson")
anova(sex.touch.m1,sex.touch.m2)

sex.bite.m1<-glmer(Bite.Reject~Sex+CT+Monkey2+offset(log(Total))+(1|Animal2)+
                      (1|ScientificName2),data=data.all2,family="poisson")
sex.bite.m2<-glmer(Bite.Reject~CT+Monkey2+offset(log(Total))+(1|Animal2)+
                      (1|ScientificName2),data=data.all2,family="poisson")
anova(sex.bite.m1,sex.bite.m2)

sex.accept.m1<-glmer(Total.Eat~Sex+CT+Monkey2+offset(log(Total))+(1|Animal2)+
                     (1|ScientificName2),data=data.all2,family="poisson")
sex.accept.m2<-glmer(Total.Eat~CT+Monkey2+offset(log(Total))+(1|Animal2)+
                     (1|ScientificName2),data=data.all2,family="poisson")
anova(sex.accept.m1,sex.accept.m2)
