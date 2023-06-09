#PELABON ET AL: Effects of population density on static allometry between horn length and body mass in mountain ungulates
#Oikos . DOI: 10.1111/oik.08726. 
#Version August 2021

#Packages#####
library(nlme)
library(lme4)
library(broom)
library(AICcmodavg)
library(lattice)

###########################################################################

# 1. Capra Ibex

# 1.1.  Data ----

data.CI<-read.table("C:\\\\Ibex.txt", header=T)#Data corrected - no duplicated entree. 
head(data.CI)

Hornm<-pmax(data.CI$lengthR, data.CI$lengthL)#selecting the longest of the two horns
year<-data.CI$yearbirth+data.CI$ageatcapture#year of capture
diff.st<-(data.CI$diff-median(data.CI$diff))*(-1)#median date of capture
data.CI<-data.frame(data.CI, year, diff.st, Hornm)


data.CI<-data.CI[data.CI$yearbirth<1999 | data.CI$yearbirth>2002,]#removing transition period
table(data.CI$yearbirth)
#We discard all data for individuals older than 10 years old for males and females
data.CIM<-subset(data.CI, data.CI$sexe=="M"& data.CI$ageatcapture<7)
data.CIF<-subset(data.CI, data.CI$sexe=="F"& data.CI$ageatcapture<7)#selection of the female

# 1.2. Analysis of the effect of density on body mass and horn length-----

    # males mass----
        #Model selection
        mod.CImm1<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):diff.st + factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm2<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):diff.st + factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm3<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm4<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):diff.st + 
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm5<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm6<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm7<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):diff.st+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm8<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm9<-lmer(bodymass~factor(ageatcapture) + factor(birthperiod) + 
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm10<-lmer(bodymass~factor(ageatcapture) + diff.st+
                           (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImm11<-lmer(bodymass~factor(ageatcapture)+
                           (1|year)+(1|ID), data = data.CIM, REML=F)
        
        
        AIC(mod.CImm1, mod.CImm2, mod.CImm3, mod.CImm4, mod.CImm5, mod.CImm6, mod.CImm7, mod.CImm8, mod.CImm9, mod.CImm10, mod.CImm11)
        
        
        #Calculation of the parameter estimates for each peirod separatly
        mod.CImm5a<-lmer(bodymass~-1+ factor(ageatcapture) + diff.st + (1|year)+(1|ID), 
                         data = data.CIM, subset= birthperiod == 1, REML=T)#model period 1
        
        summary(mod.CImm5a)
        
        mod.CImm5b<-lmer(bodymass~-1+ factor(ageatcapture) + diff.st + (1|year)+(1|ID), 
                         data = data.CIM, subset= birthperiod == 2, REML=T)#model period 2
        
        summary(mod.CImm5b)



    # males horn length----
        mod.CImh1<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):diff.st + factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh2<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):diff.st + factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh3<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh4<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):diff.st + 
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh5<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh6<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh7<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):diff.st+
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh8<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh9<-lmer(Hornm~factor(ageatcapture) + factor(birthperiod) + 
                          (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh10<-lmer(Hornm~factor(ageatcapture) + diff.st+
                           (1|year)+(1|ID), data = data.CIM, REML=F)
        
        mod.CImh11<-lmer(Hornm~factor(ageatcapture)+
                           (1|year)+(1|ID), data = data.CIM, REML=F)
        
        
        AIC(mod.CImh1, mod.CImh2, mod.CImh3, mod.CImh4, mod.CImh5, mod.CImh6, mod.CImh7, mod.CImh8, mod.CImh9, mod.CImh10, mod.CImh11)
        
        #estimation of the mean for each age class and each period
        mod.CImh3a<-lmer(Hornm~-1 + diff.st +factor(ageatcapture) + 
                           (1|year), data = data.CIM, subset = birthperiod==1, REML=T)#model period 1
        summary(mod.CImh3a)
        
        
        mod.CImh3b<-lmer(Hornm~-1 + diff.st +factor(ageatcapture) + 
                           (1|year), data = data.CIM, subset = birthperiod==2, REML=T)#model period 2
        summary(mod.CImh3b)
        

    # females mass----
        mod.CIfm1<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):diff.st + factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm2<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):diff.st + factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm3<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm4<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):diff.st + 
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm5<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm6<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF,  REML=F)
        
        mod.CIfm7<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):diff.st+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm8<-lmer(bodymass~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm9<-lmer(bodymass~factor(ageatcapture) + factor(birthperiod) + 
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm10<-lmer(bodymass~factor(ageatcapture) + diff.st+
                           (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfm11<-lmer(bodymass~factor(ageatcapture)+
                           (1|year)+(1|ID), data = data.CIF, REML=F)
        
        
        AIC(mod.CIfm1, mod.CIfm2, mod.CIfm3, mod.CIfm4, mod.CIfm5, mod.CIfm6, mod.CIfm7, mod.CIfm8, mod.CIfm9, mod.CIfm10, mod.CIfm11)
        
        #Parameter estimates 
        mod.CIfm5a<-lmer(bodymass~-1 + factor(ageatcapture) + diff.st +
                           (1|year)+(1|ID), data = data.CIF, subset=birthperiod == 1, REML=T)#model period 1
        
        summary(mod.CIfm5a)
        
        
        mod.CIfm5b<-lmer(bodymass~-1 + factor(ageatcapture) + diff.st +
                           (1|year)+(1|ID), data = data.CIF, subset=birthperiod == 2, REML=T)#model period 2
        
        summary(mod.CIfm5b)
        


    # females horn length-----

        mod.CIfh1<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):diff.st + factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        
        mod.CIfh2<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):diff.st + factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfh3<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfh4<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+ factor(ageatcapture):diff.st + 
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfh5<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfh6<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          diff.st:factor(birthperiod)+
                          (1|year)+(1|ID), data = data.CIF,  REML=F)
        
        mod.CIfh7<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          factor(ageatcapture):diff.st+
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfh8<-lmer(Hornm~factor(ageatcapture) + diff.st + factor(birthperiod) + 
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfh9<-lmer(Hornm~factor(ageatcapture) + factor(birthperiod) + 
                          (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfh10<-lmer(Hornm~factor(ageatcapture) + diff.st+
                           (1|year)+(1|ID), data = data.CIF, REML=F)
        
        mod.CIfh11<-lmer(Hornm~factor(ageatcapture)+
                           (1|year)+(1|ID), data = data.CIF, REML=F)
        
        
        AIC(mod.CIfh1, mod.CIfh2, mod.CIfh3, mod.CIfh4, mod.CIfh5, mod.CIfh6, mod.CIfh7, mod.CIfh8, mod.CIfh9, mod.CIfh10, mod.CIfh11)
        
        
        #Parameter estimates   
        
        mod.CIfh2a<-lmer(Hornm~-1 + factor(ageatcapture) + diff.st + 
                           (1|year), data = data.CIF, subset = birthperiod ==1,  REML=F)  #model period 1
        summary(mod.CIfh2a)
        
        mod.CIfh2b<-lmer(Hornm~-1 + factor(ageatcapture) + diff.st + 
                           (1|year), data = data.CIF, subset = birthperiod ==2,  REML=F)  #mmodel period 2
        summary(mod.CIfh2b)
        



# 1.3. Analysis of the static allometry all data ----       
    # males ---- 
        data.CIM2<-subset(data.CI, data.CI$sexe=="M"& data.CI$ageatcapture<6)
        age<-factor(data.CIM2$ageatcapture) 
        head(data.CIM2)
        
        
        # Model selection
        
        model.ageMa<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            factor(birthperiod):age:log(bodymass)+
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMb<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMc<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMd<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):log(bodymass) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE) 
        
        model.ageMe<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + age:factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)
        
        model.ageMf<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):log(bodymass) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMg<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):age +
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMh<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMi<-lmer(log(Hornm)~log(bodymass) + age + log(bodymass):age +
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMj<-lmer(log(Hornm)~log(bodymass) + age + factor(birthperiod) +
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)
        
        model.ageMk<-lmer(log(Hornm)~log(bodymass) + age + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        AIC(model.ageMa, model.ageMb, model.ageMc, model.ageMd, model.ageMe, model.ageMf, model.ageMg, model.ageMh, model.ageMi, model.ageMj, model.ageMk)
        AICc(model.ageMa)
        AICc(model.ageMb)
        AICc(model.ageMc)
        AICc(model.ageMd)
        AICc(model.ageMe)
        AICc(model.ageMf)
        AICc(model.ageMg)
        AICc(model.ageMh)
        AICc(model.ageMi)
        AICc(model.ageMj)
        AICc(model.ageMk)
        
        
        
        # Parameter estimates
        
        #centering the data
        bodymassC<-residuals(lm(log(bodymass)~age, data=data.CIM2))
        
        
        # Difference in slopes among age classes Period 1
        model.ageMP1a<-lmer(log(Hornm)~log(bodymass) + age +log(bodymass):age + 
                              (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==1, REML=FALSE)  
        
        model.ageMP1b<-lmer(log(Hornm)~log(bodymass) + age +
                              (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==1, REML=FALSE)  
        
        
        AIC(model.ageMP1a, model.ageMP1b)
        
        #Slope and intercept on age centered data
        model.ageMP1a<-lmer(log(Hornm)~-1+ age + bodymassC:age + 
                              (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==1, REML=TRUE)  
        
        summary(model.ageMP1a)#Parameters period 1
        
        # Difference in slopes among age classes periode 2
        model.ageMP2a<-lmer(log(Hornm)~log(bodymass) + age +log(bodymass):age + 
                              (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==2, REML=FALSE)  
        
        model.ageMP2b<-lmer(log(Hornm)~log(bodymass) + age +
                              (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==2, REML=FALSE)  
        
        AIC(model.ageMP2a, model.ageMP2b)
        
        #Best model
        model.ageMP2a<-lmer(log(Hornm)~-1+ age +bodymassC:age + 
                              (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==2, REML=TRUE)  
        
        summary(model.ageMP2a)#Parameters period 2
        
    # females ----
        
        with(data.CI, table(sexe, ageatcapture,birthperiod))
        data.CIF2<-subset(data.CI, data.CI$sexe=="F"& data.CI$ageatcapture<4)
        
        # Model selection
        age<-factor(data.CIF2$ageatcapture)
        
        model.ageFa<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            factor(birthperiod):age:log(bodymass)+
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFb<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFc<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFd<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):log(bodymass) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE) 
        
        model.ageFe<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + age:factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)
        
        model.ageFf<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):log(bodymass) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFg<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):age +
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFh<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFi<-lmer(log(Hornm)~log(bodymass) + age + log(bodymass):age +
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFj<-lmer(log(Hornm)~log(bodymass) + age + factor(birthperiod) +
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)
        
        model.ageFk<-lmer(log(Hornm)~log(bodymass) + age + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        AIC(model.ageFa, model.ageFb, model.ageFc, model.ageFd, model.ageFe, model.ageFf, model.ageFg, model.ageFh, model.ageFi, model.ageFj, model.ageFk)
        
        
        AICc(model.ageFa)
        AICc(model.ageFb)
        AICc(model.ageFc)
        AICc(model.ageFd)
        AICc(model.ageFe)
        AICc(model.ageFf)
        AICc(model.ageFg)
        AICc(model.ageFh)
        AICc(model.ageFi)
        AICc(model.ageFj)
        AICc(model.ageFk)
        
        
        
        # Parameter estimates
        
        bodymassCf<-residuals(lm(log(bodymass)~age, data=data.CIF2))
        
        # Difference in slopes among age classes Period 1
        model.ageFP1a<-lmer(log(Hornm)~log(bodymass) + age +log(bodymass):age + 
                              (1|yearbirth)+(1|ID), data=data.CIF2, subset=birthperiod ==1, REML=FALSE)  
        
        model.ageFP1b<-lmer(log(Hornm)~log(bodymass) + age +
                              (1|yearbirth)+(1|ID), data=data.CIF2, subset=birthperiod ==1, REML=FALSE)  
        
        
        AIC(model.ageFP1a, model.ageFP1b)
        
        #Best model
        model.ageFP1a<-lmer(log(Hornm)~-1+ age + bodymassCf:age + 
                              (1|yearbirth)+(1|ID), data=data.CIF2, subset=birthperiod ==1, REML=TRUE)  
        
        summary(model.ageFP1a)#Parameters period 1
        
        # Difference in slopes among age classes  periode 2 only
        model.ageFP2a<-lmer(log(Hornm)~log(bodymass) + age +log(bodymass):age + 
                              (1|yearbirth)+(1|ID), data=data.CIF2, subset=birthperiod ==2, REML=FALSE)  
        
        model.ageFP2b<-lmer(log(Hornm)~log(bodymass) + age +
                              (1|yearbirth)+(1|ID), data=data.CIF2, subset=birthperiod ==2, REML=FALSE)  
        
        AIC(model.ageFP2a, model.ageFP2b)
        
        #Best model
        
        model.ageFP2a<-lmer(log(Hornm)~ -1 + age + bodymassCf:age +
                              (1|yearbirth)+(1|ID), data=data.CIF2, subset=birthperiod ==2, REML=TRUE)  
        
        summary(model.ageFP2a) #Parameters period 2
  
# 1.4. Analysis of the static allometry 10 years with most extreme density-----
 
        data.CI<-data.CI[data.CI$yearbirth<1994 | data.CI$yearbirth>2005,]#removing 1995 - 2005
        data.CI<-data.CI[data.CI$yearbirth>1984,]#removing before 1985
        data.CI<-data.CI[data.CI$yearbirth<2016,]#removing after 2015
        
      
    # males ----    
        data.CIM2<-subset(data.CI, data.CI$sexe=="M"& data.CI$ageatcapture<6)
        age<-factor(data.CIM2$ageatcapture) 
        head(data.CIM2)
        
        
        # Model selection
        
        model.ageMa<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            factor(birthperiod):age:log(bodymass)+
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMb<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMc<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMd<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):log(bodymass) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE) 
        
        model.ageMe<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + age:factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)
        
        model.ageMf<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):log(bodymass) + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMg<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):age +
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMh<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMi<-lmer(log(Hornm)~log(bodymass) + age + log(bodymass):age +
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        model.ageMj<-lmer(log(Hornm)~log(bodymass) + age + factor(birthperiod) +
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)
        
        model.ageMk<-lmer(log(Hornm)~log(bodymass) + age + 
                            (1|yearbirth)+(1|ID), data=data.CIM2, REML=FALSE)  
        
        AIC(model.ageMa, model.ageMb, model.ageMc, model.ageMd, model.ageMe, model.ageMf, model.ageMg, model.ageMh, model.ageMi, model.ageMj, model.ageMk)
        AICc(model.ageMa)
        AICc(model.ageMb)
        AICc(model.ageMc)
        AICc(model.ageMd)
        AICc(model.ageMe)
        AICc(model.ageMf)
        AICc(model.ageMg)
        AICc(model.ageMh)
        AICc(model.ageMi)
        AICc(model.ageMj)
        AICc(model.ageMk)
        
        
        
        # Parameter estimates
        
        #centering the data
        bodymassC<-residuals(lm(log(bodymass)~age, data=data.CIM2))
      
        
      #Slope and intercept on age centered data
        model.ageMP1a<-lmer(log(Hornm)~-1+ age + bodymassC:age + 
                              (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==1, REML=TRUE)  
        
        summary(model.ageMP1a)#Parameters period 1
        
        model.ageMP2a<-lmer(log(Hornm)~-1+ age +bodymassC:age + 
                              (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==2, REML=TRUE)  
        
        summary(model.ageMP2a)#Parameters period 2
    
    # Females----    
        
        data.CIF2<-subset(data.CI, data.CI$sexe=="F"& data.CI$ageatcapture<4)
        
        # Model selection
        age<-factor(data.CIF2$ageatcapture)
        
        model.ageFa<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            factor(birthperiod):age:log(bodymass)+
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFb<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFc<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):age + log(bodymass):factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFd<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + factor(birthperiod):log(bodymass) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE) 
        
        model.ageFe<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + age:factor(birthperiod) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)
        
        model.ageFf<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):log(bodymass) + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFg<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            factor(birthperiod):age +
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFh<-lmer(log(Hornm)~log(bodymass) + factor(birthperiod) + age + 
                            log(bodymass):age + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFi<-lmer(log(Hornm)~log(bodymass) + age + log(bodymass):age +
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        model.ageFj<-lmer(log(Hornm)~log(bodymass) + age + factor(birthperiod) +
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)
        
        model.ageFk<-lmer(log(Hornm)~log(bodymass) + age + 
                            (1|yearbirth)+(1|ID), data=data.CIF2, REML=FALSE)  
        
        AIC(model.ageFa, model.ageFb, model.ageFc, model.ageFd, model.ageFe, model.ageFf, model.ageFg, model.ageFh, model.ageFi, model.ageFj, model.ageFk)
        
        
        AICc(model.ageFa)
        AICc(model.ageFb)
        AICc(model.ageFc)
        AICc(model.ageFd)
        AICc(model.ageFe)
        AICc(model.ageFf)
        AICc(model.ageFg)
        AICc(model.ageFh)
        AICc(model.ageFi)
        AICc(model.ageFj)
        AICc(model.ageFk)
        
        
        
        # Parameter estimates
        
        bodymassCf<-residuals(lm(log(bodymass)~age, data=data.CIF2))
        
        #model.ageFP1a<-lmer(log(Hornm)~-1+ age + bodymassCf:age + 
        #                      (1|yearbirth)+(1|ID), data=data.CIF2, subset=birthperiod ==1, REML=TRUE)  
        
        #summary(model.ageFP1a)#Parameters period 1
        
        model.ageFP1a<-lm(log(Hornm)~-1+ age + bodymassCf:age, 
                               data=data.CIF2, subset=birthperiod ==1)  
        
        summary(model.ageFP1a)#Parameters period 1
        
        
        model.ageFP2a<-lmer(log(Hornm)~ -1 + age + bodymassCf:age +
                              (1|yearbirth)+(1|ID), data=data.CIF2, subset=birthperiod ==2, REML=TRUE)  
        
        summary(model.ageFP2a) #Parameters period 2        
# 2. Ovis canadensis -----

# 2.1. data----
        #males
        data.BHM<-read.table("C:\\\\Ovis_M.txt", header=T)   #The transition period is removed from the data set 
                                                                                          #We remove the years before 1974 for each sex in the analysis for static allometry
        horn<-pmax(data.BHM$horn.r, data.BHM$horn.l)
        data.BHM<-data.frame(data.BHM, horn)
        
        with(data.BHM, plot(log(horn)~log(bodymass)))
        with(data.BHM, plot(lm(log(horn)~log(bodymass))))
        data.BHM<-data.BHM[-c(30),]  #Removing one outlier
        
        yearcat<-factor(data.BHM$year)
        agecat<-factor(data.BHM$age)
        periodcat<-factor(data.BHM$period)
        data.BHM<-cbind(data.BHM, yearcat, agecat, periodcat)
        data.BHM<-data.frame(data.BHM)
        
        
        
        #females
        data.BHF<-read.table("C:\\\\Ovis_F.txt", header=T)   #The transition period is removed from the data set   
        data.BHF<-subset(data.BHF, data.BHF$age<6)
  
        with(data.BHF, plot(log(horn)~log(bodymass)))
        with(data.BHF, plot(lm(log(horn)~log(bodymass))))
        
        yearcat<-factor(data.BHF$year)
        agecat<-factor(data.BHF$age)
        periodcat<-factor(data.BHF$period)
        data.BHF<-cbind(data.BHF, yearcat, agecat, periodcat)
        data.BHF<-data.frame(data.BHF)
        
# 2.2. Analysis of the effect of density on body mass and horn length----
    # male mass ---- 
        model.mm1<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm2<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm3<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:periodcat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm4<-lmer(bodymass~agecat + periodcat + date +
                          date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm5<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm6<-lmer(bodymass~agecat + periodcat + date +
                          date:periodcat+ (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm7<-lmer(bodymass~agecat + periodcat + date +date:agecat +  (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm8<-lmer(bodymass~agecat + periodcat + date + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm9<-lmer(bodymass~agecat + periodcat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm10<-lmer(bodymass~agecat + date + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mm11<-lmer(bodymass~agecat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        
        AIC(model.mm1, model.mm2, model.mm3, model.mm4, model.mm5, model.mm6, model.mm7, model.mm8, model.mm9, model.mm10, model.mm11)
        
        #Parameters 
        date.st<-data.BHM$date-median(data.BHM$date)#median center measuring date
        model.mm8b<-lmer(bodymass~-1 + date.st + agecat:periodcat + (1|yearcat)+(1|ID), data=data.BHM)
        
        summary(model.mm8b)#Parameter estimates 
        
    # male horn----
        model.mh1<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh2<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh3<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:periodcat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh4<-lmer(horn~agecat + periodcat + date +
                          date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh5<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh6<-lmer(horn~agecat + periodcat + date +
                          date:periodcat+ (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh7<-lmer(horn~agecat + periodcat + date +date:agecat +  (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh8<-lmer(horn~agecat + periodcat + date + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh9<-lmer(horn~agecat + periodcat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh10<-lmer(horn~agecat + date + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        model.mh11<-lmer(horn~agecat + (1|yearcat)+(1|ID), data=data.BHM, REML = F)
        
        
        AIC(model.mh1, model.mh2, model.mh3, model.mh4, model.mh5, model.mh6, model.mh7, model.mh8, model.mh9, model.mh10, model.mh11)
        
        #Parameters
        date.st<-data.BHM$date-median(data.BHM$date)#median center measuring date
        model.mh8b<-lmer(horn~-1 + date.st + agecat:periodcat + (1|yearcat)+(1|ID), data=data.BHM)
        
        summary(model.mh8b)#Parameter estimates 
        
     
        
    # female mass----
       
        model.fm1<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm2<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm3<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:periodcat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm4<-lmer(bodymass~agecat + periodcat + date +
                          date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm5<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm6<-lmer(bodymass~agecat + periodcat + date +
                          date:periodcat+ (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm7<-lmer(bodymass~agecat + periodcat + date +date:agecat +  (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm8<-lmer(bodymass~agecat + periodcat + date + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm9<-lmer(bodymass~agecat + periodcat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm10<-lmer(bodymass~agecat + date + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fm11<-lmer(bodymass~agecat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        
        AIC(model.fm1, model.fm2, model.fm3, model.fm4, model.fm5, model.fm6, model.fm7, model.fm8, model.fm9, model.fm10, model.fm11)
        
        #Parameters 
        date.st<-data.BHF$date-median(data.BHF$date)
        model.fm8b<-lmer(bodymass~-1 + date.st + agecat:periodcat + (1|yearcat)+(1|ID), data=data.BHF)
        
        summary(model.fm8b)#Parameter estimates
        
    # female horn----
        model.fmp1<-lmer(bodymass~-1 + agecat + agecat:date.st + (1|yearcat)+(1|ID), subset = periodcat==1, data=data.BHF)
        summary(model.fmp1)
        
        model.fmp2<-lmer(bodymass~-1 + agecat + agecat:date.st + (1|yearcat)+(1|ID), subset = periodcat==2, data=data.BHF)
        summary(model.fmp2)
        model.fh1<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh2<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh3<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:periodcat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh4<-lmer(horn~agecat + periodcat + date +
                          date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh5<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh6<-lmer(horn~agecat + periodcat + date +
                          date:periodcat+ (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh7<-lmer(horn~agecat + periodcat + date +date:agecat +  (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh8<-lmer(horn~agecat + periodcat + date + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh9<-lmer(horn~agecat + periodcat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh10<-lmer(horn~agecat + date + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        model.fh11<-lmer(horn~agecat + (1|yearcat)+(1|ID), data=data.BHF, REML = F)
        
        
        AIC(model.fh1, model.fh2, model.fh3, model.fh4, model.fh5, model.fh6, model.fh7, model.fh8, model.fh9, model.fh10, model.fh11)
        
        
        #Parameters
        date.st<-data.BHF$date-median(data.BHF$date)#median center measuring date
        
        model.fhp1<-lmer(horn~-1 + agecat + agecat:date.st + (1|yearcat)+(1|ID), subset = periodcat==1, data=data.BHF)
        summary(model.fhp1)#Parameters period 1
        
        model.fhp2<-lmer(horn~-1 + agecat + agecat:date.st + (1|yearcat)+(1|ID), subset = periodcat==2, data=data.BHF)
        summary(model.fhp2)#Parameters period 2
        
        
        
        
        
# 2.3. Analysis of the static allometry all data----
    # males ----- 
        data.BHM<-data.BHM[data.BHM$year>73 ,]#removing data before 1974
        
        #Model selection
        
        model.BHM1<-lmer(log(horn)~log(bodymass) : periodcat : agecat + 
                           log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM2<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM3<-lmer(log(horn)~log(bodymass) : periodcat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM4<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM5<-lmer(log(horn)~ log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM6<-lmer(log(horn)~log(bodymass) : periodcat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM7<-lmer(log(horn)~ log(bodymass) + periodcat : agecat +
                           periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM8<-lmer(log(horn)~log(bodymass) + periodcat + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM9<-lmer(log(horn)~log(bodymass) + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM10<-lmer(log(horn)~ log(bodymass) + periodcat + agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM11<-lmer(log(horn)~ log(bodymass) + agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        AIC(model.BHM1, model.BHM2, model.BHM3, model.BHM4, model.BHM5, model.BHM6, model.BHM7, model.BHM8, model.BHM9, model.BHM10, model.BHM11)
        
        AICc(model.BHM1)
        AICc(model.BHM2)
        AICc(model.BHM3)
        AICc(model.BHM4)
        AICc(model.BHM5)
        AICc(model.BHM6)
        AICc(model.BHM7)
        AICc(model.BHM8)
        AICc(model.BHM9)
        AICc(model.BHM10)
        AICc(model.BHM11)
        
        # Parameter estimates
        
        bodymassCm<-residuals(lm(log(bodymass)~agecat, data=data.BHM))
        
        data.BHM<-data.frame(data.BHM, bodymassCm)
        
        data.BHM1<-subset(data.BHM, data.BHM$period==1)
        data.BHM2<-subset(data.BHM, data.BHM$period==2)
        
        # Difference in slope period 1
        
        model.BHM1a<-lmer(log(horn)~ agecat + log(bodymass)+ log(bodymass): agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHM1, REML=FALSE)  
        
        model.BHM1b<-lmer(log(horn)~ agecat + log(bodymass)+ 
                            (1|yearbirth)+(1|ID), data=data.BHM1, REML=FALSE)  
        
        AIC(model.BHM1a, model.BHM1b)
        
        # Difference in slope period 2
        
        model.BHM2a<-lmer(log(horn)~ agecat + log(bodymass)+ log(bodymass): agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHM2, REML=FALSE)  
        
        model.BHM2b<-lmer(log(horn)~ agecat + log(bodymass)+ 
                            (1|yearbirth)+(1|ID), data=data.BHM2, REML=FALSE)  
        
        AIC(model.BHM2a, model.BHM2b)
        
        #parameters            
        model.BHM1<-lmer(log(horn)~ -1+  agecat + bodymassCm: agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM1, REML=FALSE)  
        summary(model.BHM1)#Parameters period 1
        
        
        model.BHM2<-lmer(log(horn)~ -1+  agecat + bodymassCm: agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM2, REML=FALSE)  
        summary(model.BHM2)#Parameters period 2
        
        
    # females ----- 
        
        data.BHF<-data.BHF[data.BHF$year>1973 ,]#removing data before 1974
        
        
        # Model selection
        model.BHF1<-lmer(log(horn)~log(bodymass) : periodcat : agecat + 
                           log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF2<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF3<-lmer(log(horn)~log(bodymass) : periodcat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF4<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF5<-lmer(log(horn)~ log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF6<-lmer(log(horn)~log(bodymass) : periodcat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF7<-lmer(log(horn)~ log(bodymass) + periodcat : agecat +
                           periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF8<-lmer(log(horn)~log(bodymass) + periodcat + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF9<-lmer(log(horn)~log(bodymass) + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF10<-lmer(log(horn)~ log(bodymass) + periodcat + agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF11<-lmer(log(horn)~ log(bodymass) + agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        
        AICc(model.BHF1)
        AICc(model.BHF2)
        AICc(model.BHF3)
        AICc(model.BHF4)
        AICc(model.BHF5)
        AICc(model.BHF6)
        AICc(model.BHF7)
        AICc(model.BHF8)
        AICc(model.BHF9)
        AICc(model.BHF10)
        AICc(model.BHF11)
        
        
        # Parameter estimates
        
        bodymassCf<-residuals(lm(log(bodymass)~agecat, data=data.BHF))
        
        data.BHF<-data.frame(data.BHF, bodymassCf)
        
        data.BHF1<-subset(data.BHF, data.BHF$period==1)
        data.BHF2<-subset(data.BHF, data.BHF$period==2)
        
        
        # Difference in slope period 1
        
        model.BHF1a<-lmer(log(horn)~ agecat + log(bodymass)+ log(bodymass): agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHF1, REML=FALSE)  
        
        model.BHF1b<-lmer(log(horn)~ agecat + log(bodymass)+ 
                            (1|yearbirth)+(1|ID), data=data.BHF1, REML=FALSE)  
        
        AIC(model.BHF1a, model.BHF1b)
        
        # Difference in slope period 2
        
        model.BHF2a<-lmer(log(horn)~ agecat + log(bodymass)+ log(bodymass): agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHF2, REML=FALSE)  
        
        model.BHF2b<-lmer(log(horn)~ agecat + log(bodymass)+ 
                            (1|yearbirth)+(1|ID), data=data.BHF2, REML=FALSE)  
        
        AIC(model.BHF2a, model.BHF2b)
        
        #parameters            
        model.BHF1<-lmer(log(horn)~ -1+  agecat + bodymassCf: agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF1, REML=FALSE)  
        summary(model.BHF1)#Parameters period 1
        
        
        model.BHF2<-lmer(log(horn)~ -1+  agecat + bodymassCf: agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF2, REML=FALSE)  
        summary(model.BHF2)#Parameters period 2
        
# 2.4. Analysis static allometry 10 years with most extreme density-----
        ##10 years low density 1974 - 1983 -- 10 years with high denisty 1988 - 1997
    # males-----
        data.BHM<-data.BHM[data.BHM$year<85 | data.BHM$year>87,]#removing transition period 84 - 87
        data.BHM<-data.BHM[data.BHM$year<100 ,]#removing after 2000
        
        with(data.BHM, table(year))
        
        #Model selection
        
        model.BHM1<-lmer(log(horn)~log(bodymass) : periodcat : agecat + 
                           log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM2<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM3<-lmer(log(horn)~log(bodymass) : periodcat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM4<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM5<-lmer(log(horn)~ log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM6<-lmer(log(horn)~log(bodymass) : periodcat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM7<-lmer(log(horn)~ log(bodymass) + periodcat : agecat +
                           periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM8<-lmer(log(horn)~log(bodymass) + periodcat + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM9<-lmer(log(horn)~log(bodymass) + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM10<-lmer(log(horn)~ log(bodymass) + periodcat + agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        model.BHM11<-lmer(log(horn)~ log(bodymass) + agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHM, REML=FALSE)  
        
        AIC(model.BHM1, model.BHM2, model.BHM3, model.BHM4, model.BHM5, model.BHM6, model.BHM7, model.BHM8, model.BHM9, model.BHM10, model.BHM11)
        
        AICc(model.BHM1)
        AICc(model.BHM2)
        AICc(model.BHM3)
        AICc(model.BHM4)
        AICc(model.BHM5)
        AICc(model.BHM6)
        AICc(model.BHM7)
        AICc(model.BHM8)
        AICc(model.BHM9)
        AICc(model.BHM10)
        AICc(model.BHM11)
        
        # Parameter estimates
        
        bodymassCm<-residuals(lm(log(bodymass)~agecat, data=data.BHM))
        
        data.BHM<-data.frame(data.BHM, bodymassCm)
        
        data.BHM1<-subset(data.BHM, data.BHM$period==1)
        data.BHM2<-subset(data.BHM, data.BHM$period==2)
        
       
        #parameters            
        model.BHM1<-lmer(log(horn)~ -1+  agecat + bodymassCm: agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM1, REML=FALSE)  
        summary(model.BHM1)#Parameters period 1
        
        
        model.BHM2<-lmer(log(horn)~ -1+  agecat + bodymassCm: agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHM2, REML=FALSE)  
        summary(model.BHM2)#Parameters period 2
        
    # females----
        data.BHF<-data.BHF[data.BHF$year<1985 | data.BHF$year>1987,]#removing transition period 84 - 87
        data.BHF<-data.BHF[data.BHF$year<2000 ,]#removing after 2000
        
        with(data.BHF, table(year))
        
        
        
        # Model selection
        model.BHF1<-lmer(log(horn)~log(bodymass) : periodcat : agecat + 
                           log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF2<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF3<-lmer(log(horn)~log(bodymass) : periodcat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF4<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF5<-lmer(log(horn)~ log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF6<-lmer(log(horn)~log(bodymass) : periodcat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF7<-lmer(log(horn)~ log(bodymass) + periodcat : agecat +
                           periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF8<-lmer(log(horn)~log(bodymass) + periodcat + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF9<-lmer(log(horn)~log(bodymass) + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF10<-lmer(log(horn)~ log(bodymass) + periodcat + agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        model.BHF11<-lmer(log(horn)~ log(bodymass) + agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHF, REML=FALSE)  
        
        
        AICc(model.BHF1)
        AICc(model.BHF2)
        AICc(model.BHF3)
        AICc(model.BHF4)
        AICc(model.BHF5)
        AICc(model.BHF6)
        AICc(model.BHF7)
        AICc(model.BHF8)
        AICc(model.BHF9)
        AICc(model.BHF10)
        AICc(model.BHF11)
        
        
        # Parameter estimates
        
        bodymassCf<-residuals(lm(log(bodymass)~agecat, data=data.BHF))
        
        data.BHF<-data.frame(data.BHF, bodymassCf)
        
        data.BHF1<-subset(data.BHF, data.BHF$period==1)
        data.BHF2<-subset(data.BHF, data.BHF$period==2)
      
        
        #parameters            
        model.BHF1<-lmer(log(horn)~ -1+  agecat + bodymassCf: agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF1, REML=FALSE)  
        summary(model.BHF1)#Parameters period 1
        
        
        model.BHF2<-lmer(log(horn)~ -1+  agecat + bodymassCf: agecat + 
                           (1|yearbirth)+(1|ID), data=data.BHF2, REML=FALSE)  
        summary(model.BHF2)#Parameters period 2
        
        

# 3. Rupicapra rupicapra----
        
# 3.1. data----
        
        data.rr2<-read.table("C:\\\\chamois.txt", header=T)
        head(data.rr2)
        
        yr.birth<-with(data.rr2, yr-age)
        data.rr2<-data.frame(data.rr2, yr.birth)
        
        data.rr2<-data.rr2[data.rr2$yr<1998 | data.rr2$yr>2000,]#removing individuals in the transition period
        data.rr2<-subset(data.rr2, data.rr2$age>0)#removing yearling 
        data.rr2<-subset(data.rr2, data.rr2$age<9)#removing individuals older than 8 year (not enough in each age class)
        data.rr2<-subset(data.rr2, yr.birth>1981)#we remove all individuals born during the episode of Kerato
        
        with(data.rr2, table(yr))
        
        horn<-pmax(data.rr2$hornR, data.rr2$hornL)#selecting the longest of the two horns
        period<-with(data.rr2, ifelse(density=="low", 1,2))#defining the two period until 1998 / after 1998
        dayk<-with(data.rr2, daynr-median(daynr))# Centering the hunting date around the median hunting date
        dayk2<-dayk^2
        data.rr2<-data.frame(data.rr2, horn, period, dayk, dayk2)
        
        data.rrM<-subset(data.rr2, data.rr2$sex=="M")
        head(data.rrM)
        data.rrF<-subset(data.rr2, data.rr2$sex=="F")
        

# 3.2. Analysis of the effect of density on body mass and horn length----
    # males mass----
        mod.rrmm1<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + dayk2:factor(period)+
                          dayk:factor(age):factor(period) + dayk2:factor(age):factor(period)+
                          (1|season), data = data.rrM, REML=F)
        
        
        mod.rrmm2<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + dayk2:factor(period)+
                          (1|season), data = data.rrM, REML=F)
        
        mod.rrmm3<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + 
                          (1|season), data = data.rrM, REML=F)
        
        mod.rrmm4<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          (1|season), data = data.rrM, REML=F)
        
        mod.rrmm5<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + (1|season), data = data.rrM, REML=F)
        
        mod.rrmm6<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + (1|season), data = data.rrM, REML=F)
        
        mod.rrmm7<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + (1|season), data = data.rrM, REML=F)
        
        mod.rrmm8<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          (1|season), data = data.rrM, REML=F)
        
        mod.rrmm9<-lmer(mass~factor(age) + factor(period)+ (1|season), data = data.rrM, REML=F)
        
        mod.rrmm10<-lmer(mass~factor(age) + (1|season), data = data.rrM, REML=F)
        
        
        AIC(mod.rrmm1, mod.rrmm2, mod.rrmm3, mod.rrmm4, mod.rrmm5, mod.rrmm6, mod.rrmm7, mod.rrmm8, mod.rrmm9, mod.rrmm10 )
        
        
        #  Parameter estimates
        
        mod.rrmm5a<-lmer(mass~-1+ factor(age) +  dayk + dayk2 + dayk:factor(age) + 
                           (1|season), data = data.rrM, subset = period ==1, REML=T)
        
        summary(mod.rrmm5a)
        
        mod.rrmm5b<-lmer(mass~-1+ factor(age) +  dayk + dayk2 + dayk:factor(age) + 
                           (1|season), data = data.rrM, subset = period ==2 , REML=T)
        
        summary(mod.rrmm5b)
        
        
        
        
    # males horn length----
        
        mod.rrmh1<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + dayk2:factor(period)+
                          dayk:factor(age):factor(period) + dayk2:factor(age):factor(period)+
                          (1|season), data = data.rrM, REML=F)
        
        
        mod.rrmh2<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + dayk2:factor(period)+
                          (1|season), data = data.rrM, REML=F)
        
        mod.rrmh3<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + 
                          (1|season), data = data.rrM, REML=F)
        
        mod.rrmh4<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          (1|season), data = data.rrM, REML=F)
        
        mod.rrmh5<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + (1|season), data = data.rrM, REML=F)
        
        mod.rrmh6<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + (1|season), data = data.rrM, REML=F)
        
        mod.rrmh7<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + (1|season), data = data.rrM, REML=F)
        
        mod.rrmh8<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          (1|season), data = data.rrM, REML=F)
        
        mod.rrmh9<-lmer(horn~factor(age) + factor(period)+ (1|season), data = data.rrM, REML=F)
        
        mod.rrmh10<-lmer(horn~factor(age) + (1|season), data = data.rrM, REML=F)
        
        
        AIC(mod.rrmh1, mod.rrmh2, mod.rrmh3, mod.rrmh4, mod.rrmh5, mod.rrmh6, mod.rrmh7, mod.rrmh8, mod.rrmh9, mod.rrmh10 )
        
        #parameter estimates 
        
        mod.rrmh7a<-lmer(horn~-1 + factor(age) + dayk + (1|season), data = data.rrM, subset = period ==1, REML=T)
        summary(mod.rrmh7a)
        
        mod.rrmh7b<-lmer(horn~-1 + factor(age) + dayk + (1|season), data = data.rrM, subset = period ==2, REML=T)
        summary(mod.rrmh7b)
        
        
        
        
    # females mass----
        
        
        mod.rrfm1<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + dayk2:factor(period)+
                          dayk:factor(age):factor(period) + dayk2:factor(age):factor(period)+
                          (1|season), data = data.rrF, REML=F)
        
        
        mod.rrfm2<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + dayk2:factor(period)+
                          (1|season), data = data.rrF, REML=F)
        
        mod.rrfm3<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + 
                          (1|season), data = data.rrF, REML=F)
        
        mod.rrfm4<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          (1|season), data = data.rrF, REML=F)
        
        mod.rrfm5<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + (1|season), data = data.rrF, REML=F)
        
        mod.rrfm6<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + (1|season), data = data.rrF, REML=F)
        
        mod.rrfm7<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + (1|season), data = data.rrF, REML=F)
        
        mod.rrfm8<-lmer(mass~factor(age) + factor(period) + factor(age):factor(period)+ 
                          (1|season), data = data.rrF, REML=F)
        
        mod.rrfm9<-lmer(mass~factor(age) + factor(period)+ (1|season), data = data.rrF, REML=F)
        
        mod.rrfm10<-lmer(mass~factor(age) + (1|season), data = data.rrF, REML=F)
        
        
        AIC(mod.rrfm1, mod.rrfm2, mod.rrfm3, mod.rrfm4, mod.rrfm5, mod.rrfm6, mod.rrfm7, mod.rrfm8, mod.rrfm9, mod.rrfm10 )
        
        
        #  Parameter estimates
        
        mod.rrfm6a<-lmer(mass~-1+ factor(age) +dayk + dayk2 + (1|season), data = data.rrF, subset= period ==1, REML=T)
        
        
        summary(mod.rrfm6a)
        
        mod.rrfm6b<-lmer(mass~-1+factor(age) +dayk + dayk2 + (1|season), data = data.rrF, subset= period ==2, REML=T)
        
        summary(mod.rrfm6b)
        
        
        
        
    # females horn length----
        
        mod.rrfh1<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + dayk2:factor(period)+
                          dayk:factor(age):factor(period) + dayk2:factor(age):factor(period)+
                          (1|season), data = data.rrF, REML=F)
        
        
        mod.rrfh2<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + dayk2:factor(period)+
                          (1|season), data = data.rrF, REML=F)
        
        mod.rrfh3<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          dayk:factor(period) + 
                          (1|season), data = data.rrF, REML=F)
        
        mod.rrfh4<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + dayk2:factor(age)+
                          (1|season), data = data.rrF, REML=F)
        
        mod.rrfh5<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + dayk:factor(age) + (1|season), data = data.rrF, REML=F)
        
        mod.rrfh6<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + dayk2 + (1|season), data = data.rrF, REML=F)
        
        mod.rrfh7<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          dayk + (1|season), data = data.rrF, REML=F)
        
        mod.rrfh8<-lmer(horn~factor(age) + factor(period) + factor(age):factor(period)+ 
                          (1|season), data = data.rrF, REML=F)
        
        mod.rrfh9<-lmer(horn~factor(age) + factor(period)+ (1|season), data = data.rrF, REML=F)
        
        mod.rrfh10<-lmer(horn~factor(age) + (1|season), data = data.rrF, REML=F)
        
        
        AIC(mod.rrfh1, mod.rrfh2, mod.rrfh3, mod.rrfh4, mod.rrfh5, mod.rrfh6, mod.rrfh7, mod.rrfh8, mod.rrfh9, mod.rrfh10 )
        
        #parameter estimates 
        
        mod.rrfh5a<-lmer(horn~-1 + factor(age) + dayk +  
                           (1|season), data = data.rrF, subset = period ==1, REML=T)
        
        summary(mod.rrfh5a)
        
        
        mod.rrfh5b<-lmer(horn~-1 + factor(age) + dayk + 
                           (1|season), data = data.rrF, subset = period ==2, REML=T)
        
        summary(mod.rrfh5b)
        
        
        
        

# 3.3. Correcting for the effect of hunting date on mass and horn length----
      
    # male mass----
        
        mod.rr_mm1<-lmer(mass~ factor(age) + dayk + dayk2+ dayk:age + dayk:period + dayk2:age + dayk2:period + 
                           (1|season), data = data.rrM, REML=T)#test for the effect of hunting date 
        summary(mod.rr_mm1)
        
        mod.rr_mm2<-lmer(mass~ factor(age) + dayk + dayk2+ dayk:age + dayk:period + dayk2:period +
                           (1|season), data = data.rrM, REML=T)#test for the effect of hunting date 
        summary(mod.rr_mm2)
        
        mod.rr_mm2b<-lmer(mass~ -1 + factor(age) + dayk + dayk2+ dayk:age + dayk:period + dayk2:period +
                            (1|season), data = data.rrM, REML=T)#test for the effect of hunting date 
        summary(mod.rr_mm2b)
        
        
        masscorr<-data.rrM$mass-((fixef(mod.rr_mm2b)[9]+fixef(mod.rr_mm2b)[12]*data.rrM$period+fixef(mod.rr_mm2b)[11]*data.rrM$age)*data.rrM$dayk+
                                   (fixef(mod.rr_mm2b)[10]+fixef(mod.rr_mm2b)[13]*data.rrM$period)*data.rrM$dayk2)#correcting the mass for hunting date
        
        data.rrM<-data.frame(data.rrM, masscorr)
        with(data.rrM, plot(mass, masscorr))
        with(data.rrM, mean(100*abs(mass-masscorr)/mass))#average correction in % of the mean
        
        
        
       
        
    # male horn----
        
        mod.rr_mh1<-lmer(horn~ factor(age) + dayk + dayk2+ dayk:age + dayk:period + dayk2:age + dayk2:period + 
                           (1|season), data = data.rrM, REML=T)#test for the effect of hunting date 
        summary(mod.rr_mh1)
        
        mod.rr_mh2<-lmer(horn~ factor(age) + dayk + dayk2+ dayk:age + dayk:period + 
                           (1|season), data = data.rrM, REML=T)#test for the effect of hunting date 
        summary(mod.rr_mh2)
        
        mod.rr_mh3<-lmer(horn~ factor(age) + dayk +
                           (1|season), data = data.rrM, REML=T)#test for the effect of hunting date 
        summary(mod.rr_mh3)
        
        
        horncorr<-data.rrM$horn-((fixef(mod.rr_mh3)[9]*data.rrM$dayk))#correcting the mass for hunting date
        data.rrM<-data.frame(data.rrM, horncorr)
        with(data.rrM, plot(horn, horncorr))
        with(data.rrM, mean(100*(horn-horncorr)/horn))
        
        
    # female mass----
        
        
        mod.rr_Fm1<-lmer(mass~ factor(age) + dayk + dayk2+ dayk:age + dayk:period + dayk2:age + dayk2:period + 
                           (1|season), data = data.rrF, REML=T)#test for the effect of hunting date 
        summary(mod.rr_Fm1)
        
        mod.rr_Fm2<-lmer(mass~ factor(age) + dayk + dayk2+ dayk:age + dayk:period + dayk2:period +
                           (1|season), data = data.rrF, REML=T)#test for the effect of hunting date 
        summary(mod.rr_Fm2)
        
        
        mod.rr_Fm3<-lmer(mass~ factor(age) + dayk + dayk2+ dayk:period + dayk2:period +
                           (1|season), data = data.rrF, REML=T)#test for the effect of hunting date 
        summary(mod.rr_Fm3)
        
        
        mod.rr_Fm3b<-lmer(mass~ -1 + factor(age) + dayk + dayk2+ dayk:period + dayk2:period +
                            (1|season), data = data.rrF, REML=T)#test for the effect of hunting date 
        summary(mod.rr_Fm3b)
        
        
        masscorr<-data.rrF$mass-((fixef(mod.rr_Fm3b)[9]+fixef(mod.rr_Fm3b)[11]*data.rrF$period)*data.rrF$dayk+
                                   (fixef(mod.rr_Fm3b)[10]+fixef(mod.rr_Fm3b)[12]*data.rrF$period)*data.rrF$dayk2)#correcting the mass for hunting date
        
        data.rrF<-data.frame(data.rrF, masscorr)
        with(data.rrF, plot(mass, masscorr))
        with(data.rrF, mean(100*abs(mass-masscorr)/mass))#average correction in % of the mean
        
        
        
    # female horn ----
        
        mod.rr_Fh1<-lmer(horn~ factor(age) + dayk + dayk2+ dayk:age + dayk:period + dayk2:age + dayk2:period + 
                           (1|season), data = data.rrF, REML=T)#test for the effect of hunting date 
        summary(mod.rr_Fh1)
        
        mod.rr_Fh2b<-lmer(horn~ -1+ factor(age) + dayk + dayk2+ dayk:age + dayk2:age + 
                            (1|season), data = data.rrF, REML=T)#test for the effect of hunting date 
        summary(mod.rr_Fh2b)
        
        
        mod.rr_Fh2l<-lm(horn~ factor(age) + dayk + dayk2+ dayk:age + dayk2:age, data = data.rrF)#test for the effect of hunting date 
        summary(mod.rr_Fh2l)
        mod.rr_Fh3l<-lm(horn~ factor(age), data = data.rrF)#test for the effect of hunting date 
        summary(mod.rr_Fh3l)
        
        
# 3.4. Analysis of the static allometry all data ----

        data.rrM[data.rrM$age>=5,"age"] <- 5#We include individuals older than 5 into the 5 year age class since these individuals do not differ in mass or horn length
        data.rrF[data.rrF$age>=5,"age"] <- 5# idem females
    
    # males -----
        
       
        table(data.rrM$age)
        head(data.rrM)
        age2<-factor(data.rrM$age)
        data.rrM<-data.frame(data.rrM, age2)
        with(data.rrM, cohort<-factor(cohort))
 
        model.ageMa<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):age2 + log(masscorr):factor(period) + 
                            factor(period):age2:log(masscorr)+
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMb<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):age2 + log(masscorr):factor(period) + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMc<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):age2 + log(masscorr):factor(period) + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMd<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):log(masscorr) + 
                            (1|cohort), data=data.rrM, REML=FALSE) 
        
        model.ageMe<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + age2:factor(period) + 
                            (1|cohort), data=data.rrM, REML=FALSE)
        
        model.ageMf<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):log(masscorr) + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMg<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):age2 +
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMh<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMi<-lmer(log(horn)~log(masscorr) + age2 + log(masscorr):age2 +
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMj<-lmer(log(horn)~log(masscorr) + age2 + factor(period) +
                            (1|cohort), data=data.rrM, REML=FALSE)
        
        model.ageMk<-lmer(log(horn)~log(masscorr) + age2 + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        AIC(model.ageMa, model.ageMb, model.ageMc, model.ageMd, model.ageMe, 
            model.ageMf, model.ageMg, model.ageMh, model.ageMi, model.ageMj, model.ageMk)
        
        # Difference in slope period 1
        
        mod.rrM1a<-lmer(log(horn)~log(masscorr) + age2 +log(masscorr):age2 + 
                          (1|cohort), data=data.rrM, subset=period==1, REML=F) 
        
        mod.rrM1b<-lmer(log(horn)~log(masscorr) + age2 + 
                          (1|cohort), data=data.rrM, subset=period==1, REML=F) 
        
        AIC(mod.rrM1a, mod.rrM1b)
        
        # Difference in slope period 2
        
        mod.rrM2a<-lmer(log(horn)~log(masscorr) + age2 +log(masscorr):age2 + 
                          (1|cohort), data=data.rrM, subset=period==2, REML=F) 
        
        mod.rrM2b<-lmer(log(horn)~log(masscorr) + age2 + 
                          (1|cohort), data=data.rrM, subset=period==2, REML=F) 
        
        AIC(mod.rrM2a, mod.rrM2b)
        
        
        # Parameter estimates
        
        masscorrCm<-residuals(lm(log(masscorr)~age2, data = data.rrM))
        
        mod.rrMg1a<-lmer(log(horn)~-1 +age2 +masscorrCm:age2 + 
                           (1|cohort), data=data.rrM, subset=period==1, REML=TRUE) 
        summary(mod.rrMg1a)
        
        
        mod.rrMg2a<-lmer(log(horn)~-1 +age2 +masscorrCm:age2 + 
                           (1|cohort), data=data.rrM, subset=period==2, REML=TRUE) 
        summary(mod.rrMg2a)
        
    # females-----
  
        
        age2<-factor(data.rrF$age)
       
        data.rrF<-data.frame(data.rrF, age2)
        with(data.rrF, cohort<-factor(cohort))
        
        model.ageFa<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):age2 + log(masscorr):factor(period) + 
                            factor(period):age2:log(masscorr)+
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFb<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):age2 + log(masscorr):factor(period) + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFc<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):age2 + log(masscorr):factor(period) + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFd<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):log(masscorr) + 
                            (1|cohort), data=data.rrF, REML=FALSE) 
        
        model.ageFe<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + age2:factor(period) + 
                            (1|cohort), data=data.rrF, REML=FALSE)
        
        model.ageFf<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):log(masscorr) + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFg<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):age2 +
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFh<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFi<-lmer(log(horn)~log(masscorr) + age2 + log(masscorr):age2 +
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFj<-lmer(log(horn)~log(masscorr) + age2 + factor(period) +
                            (1|cohort), data=data.rrF, REML=FALSE)
        
        model.ageFk<-lmer(log(horn)~log(masscorr) + age2 + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        AIC(model.ageFa, model.ageFb, model.ageFc, model.ageFd, model.ageFe, 
            model.ageFf, model.ageFg, model.ageFh, model.ageFi, model.ageFj, model.ageFk)
        
        
        # Difference in slope period 1
        
        mod.rrF1a<-lmer(log(horn)~log(masscorr) + age2 +log(masscorr):age2 + 
                          (1|cohort), data=data.rrF, subset=period==1, REML=F) 
        
        mod.rrF1b<-lmer(log(horn)~log(masscorr) + age2 + 
                          (1|cohort), data=data.rrF, subset=period==1, REML=F) 
        
        AIC(mod.rrF1a, mod.rrF1b)
        
        # Difference in slope period 2
        
        mod.rrF2a<-lmer(log(horn)~log(masscorr) + age2 +log(masscorr):age2 + 
                          (1|cohort), data=data.rrF, subset=period==2, REML=F) 
        
        mod.rrF2b<-lmer(log(horn)~log(masscorr) + age2 + 
                          (1|cohort), data=data.rrF, subset=period==2, REML=F) 
        
        AIC(mod.rrF2a, mod.rrF2b)
        
        summary(mod.rrF2a)
        
        # Parameter estimates
    
        masscorrCf<-residuals(lm(log(masscorr)~age2, data = data.rrF))
        
        mod.rrF1a<-lmer(log(horn)~-1 +age2 +masscorrCf:age2 + 
                          (1|cohort), data=data.rrF, subset=period==1, REML=TRUE) 
        summary(mod.rrF1a)
        
        
        mod.rrF2a<-lmer(log(horn)~-1 +age2 +masscorrCf:age2 + 
                          (1|cohort), data=data.rrF, subset=period==2, REML=TRUE) 
        summary(mod.rrF2a)
        
# 3.5. Analysis static allometry 10 years with most extreme density-----
        ##10 years low density 1983 - 1995 -- 10 years with high denisty 2001 2010
    # males------
        
        data.rrM<-data.rrM[data.rrM$yr<2011 ,]#removing individuals after 2010
        data.rrM<-data.rrM[data.rrM$yr<1996 | data.rrM$yr>2000,]#removing individuals in the transition period
        table(data.rrM$yr)
        
        model.ageMa<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):age2 + log(masscorr):factor(period) + 
                            factor(period):age2:log(masscorr)+
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMb<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):age2 + log(masscorr):factor(period) + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMc<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):age2 + log(masscorr):factor(period) + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMd<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):log(masscorr) + 
                            (1|cohort), data=data.rrM, REML=FALSE) 
        
        model.ageMe<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + age2:factor(period) + 
                            (1|cohort), data=data.rrM, REML=FALSE)
        
        model.ageMf<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):log(masscorr) + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMg<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):age2 +
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMh<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMi<-lmer(log(horn)~log(masscorr) + age2 + log(masscorr):age2 +
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        model.ageMj<-lmer(log(horn)~log(masscorr) + age2 + factor(period) +
                            (1|cohort), data=data.rrM, REML=FALSE)
        
        model.ageMk<-lmer(log(horn)~log(masscorr) + age2 + 
                            (1|cohort), data=data.rrM, REML=FALSE)  
        
        
        
        AIC(model.ageMa, model.ageMb, model.ageMc, model.ageMd, model.ageMe, 
            model.ageMf, model.ageMg, model.ageMh, model.ageMi, model.ageMj, model.ageMk)
        
        
        # Parameter estimates
        
        masscorrCm<-residuals(lm(log(masscorr)~age2, data = data.rrM))
        
        mod.rrMg1a<-lmer(log(horn)~-1 +age2 +masscorrCm:age2 + 
                           (1|cohort), data=data.rrM, subset=period==1, REML=TRUE) 
        summary(mod.rrMg1a)
        
        
        mod.rrMg2a<-lmer(log(horn)~-1 +age2 +masscorrCm:age2 + 
                           (1|cohort), data=data.rrM, subset=period==2, REML=TRUE) 
        summary(mod.rrMg2a)
        
        
        
    # females-----
        data.rrF<-data.rrF[data.rrF$yr<2011 ,]#removing individuals after 2010
        data.rrF<-data.rrF[data.rrF$yr<1996 | data.rrF$yr>2000,]#removing individuals in the transition period
        table(data.rrF$yr)
        
        age2<-factor(data.rrF$age)
        data.rrF<-data.frame(data.rrF, age2)
        
        
        model.ageFa<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):age2 + log(masscorr):factor(period) + 
                            factor(period):age2:log(masscorr)+
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFb<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):age2 + log(masscorr):factor(period) + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFc<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):age2 + log(masscorr):factor(period) + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFd<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + factor(period):log(masscorr) + 
                            (1|cohort), data=data.rrF, REML=FALSE) 
        
        model.ageFe<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + age2:factor(period) + 
                            (1|cohort), data=data.rrF, REML=FALSE)
        
        model.ageFf<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):log(masscorr) + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFg<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            factor(period):age2 +
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFh<-lmer(log(horn)~log(masscorr) + factor(period) + age2 + 
                            log(masscorr):age2 + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFi<-lmer(log(horn)~log(masscorr) + age2 + log(masscorr):age2 +
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        model.ageFj<-lmer(log(horn)~log(masscorr) + age2 + factor(period) +
                            (1|cohort), data=data.rrF, REML=FALSE)
        
        model.ageFk<-lmer(log(horn)~log(masscorr) + age2 + 
                            (1|cohort), data=data.rrF, REML=FALSE)  
        
        AIC(model.ageFa, model.ageFb, model.ageFc, model.ageFd, model.ageFe, 
            model.ageFf, model.ageFg, model.ageFh, model.ageFi, model.ageFj, model.ageFk)
        
        
        # Parameter estimates
        
        masscorrCf<-residuals(lm(log(masscorr)~age2, data = data.rrF))
        
        mod.rrF1a<-lmer(log(horn)~-1 +age2 +masscorrCf:age2 + 
                          (1|cohort), data=data.rrF, subset=period==1, REML=TRUE) 
        summary(mod.rrF1a)
        
        
        mod.rrF2a<-lmer(log(horn)~-1 +age2 +masscorrCf:age2 + 
                          (1|cohort), data=data.rrF, subset=period==2, REML=TRUE) 
        summary(mod.rrF2a)
        
        
        
# 4. Mountain goat Oreamnos----
# 4.1. data----        
        data.OA<-read.table("C:\\\\Oreamnos.txt", header=T)#Data with only animals of 1 year or more 
        head(data.OA)
        horn<-pmax(data.OA$right_horn, data.OA$left_horn)
        data.OA<-data.frame(data.OA, horn)
        
        data.OA <- data.OA[-c(531, 41),] #two outliers
        data.OA1<-data.OA[data.OA$year<1997| data.OA$year>2000,]#removing observation year 1997 - 2000
        data.OA1<-data.OA1[data.OA1$year<2012| data.OA1$year>2013,]#removing observation year 2012 - 2013
        
        ##Period 1 = low desnity;  2 = High density
        
        yearcat<-factor(data.OA1$year)
        agecat<-factor(data.OA1$age)
        periodcat<-factor(data.OA1$period)
        data.OA1<-cbind(data.OA1, yearcat, agecat, periodcat)
        
        data.OA<-data.frame(data.OA1)
        
        
        data.OAM<-subset(data.OA, data.OA$sexe==2 & data.OA$age<6)#Males 
        data.OAF<-subset(data.OA, data.OA$sexe==1  & data.OA$age<6)#Females
        
# 4.2. Analysis of the effect of density on body mass and horn length----
    # male mass----
        model.mm1<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm2<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm3<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:periodcat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm4<-lmer(bodymass~agecat + periodcat + date +
                          date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm5<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm6<-lmer(bodymass~agecat + periodcat + date +
                          date:periodcat+ (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm7<-lmer(bodymass~agecat + periodcat + date +date:agecat +  (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm8<-lmer(bodymass~agecat + periodcat + date + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm9<-lmer(bodymass~agecat + periodcat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm10<-lmer(bodymass~agecat + date + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mm11<-lmer(bodymass~agecat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        
        AIC(model.mm1, model.mm2, model.mm3, model.mm4, model.mm5, model.mm6, model.mm7, model.mm8, model.mm9, model.mm10, model.mm11)
        
        
        date.st<-data.OAM$date-median(data.OAM$date)
        
        model.mmP1<-lmer(bodymass~-1 + agecat+ date.st:agecat+ (1|yearcat)+(1|ID), data=data.OAM, subset= period ==1)
        summary(model.mmP1)
        
        model.mmP2<-lmer(bodymass~-1 + agecat+ date.st:agecat+ (1|yearcat)+(1|ID), data=data.OAM, subset= period ==2)
        summary(model.mmP2)
        
        
    # male horn----
        model.mh1<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh2<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh3<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:periodcat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh4<-lmer(horn~agecat + periodcat + date +
                          date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh5<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh6<-lmer(horn~agecat + periodcat + date +
                          date:periodcat+ (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh7<-lmer(horn~agecat + periodcat + date +date:agecat +  (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh8<-lmer(horn~agecat + periodcat + date + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh9<-lmer(horn~agecat + periodcat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh10<-lmer(horn~agecat + date + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        model.mh11<-lmer(horn~agecat + (1|yearcat)+(1|ID), data=data.OAM, REML = F)
        
        
        AIC(model.mh1, model.mh2, model.mh3, model.mh4, model.mh5, model.mh6, model.mh7, model.mh8, model.mh9, model.mh10, model.mh11)
        
        
        date.st<-data.OAM$date-median(data.OAM$date)
        
        model.mhP1<-lmer(horn~-1 + agecat+ date.st:agecat+ (1|yearcat)+(1|ID), data=data.OAM, subset= period ==1)
        summary(model.mhP1)
        
        model.mhP2<-lmer(horn~-1 + agecat+ date.st:agecat+ (1|yearcat)+(1|ID), data=data.OAM, subset= period ==2)
        summary(model.mhP2)
        
    # female mass ----

        model.fm1<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4, REML = F)
        
        model.fm2<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + (1|yearcat)+(1|ID), data=data.OAF,subset= agecat!=4, REML = F)
        
        model.fm3<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + date:periodcat + (1|yearcat)+(1|ID), data=data.OAF,subset= agecat!=4, REML = F)
        
        model.fm4<-lmer(bodymass~agecat + periodcat + date +
                          date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.OAF,subset= agecat!=4,REML = F)
        
        model.fm5<-lmer(bodymass~agecat + periodcat + date +
                          agecat:periodcat + (1|yearcat)+(1|ID), data=data.OAF,subset= agecat!=4, REML = F)
        
        model.fm6<-lmer(bodymass~agecat + periodcat + date +
                          date:periodcat+ (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fm7<-lmer(bodymass~agecat + periodcat + date +date:agecat +  (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fm8<-lmer(bodymass~agecat + periodcat + date + (1|yearcat)+(1|ID), data=data.OAF,subset= agecat!=4, REML = F)
        
        model.fm9<-lmer(bodymass~agecat + periodcat + (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fm10<-lmer(bodymass~agecat + date + (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fm11<-lmer(bodymass~agecat + (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        
        AIC(model.fm1, model.fm2, model.fm3, model.fm4, model.fm5, model.fm6, model.fm7, model.fm8, model.fm9, model.fm10, model.fm11)
        
        
        date.st<-data.OAF$date-median(data.OAF$date)
        
        model.fmP1<-lmer(bodymass~-1 + agecat+ date.st:agecat+ (1|yearcat)+(1|ID), data=data.OAF, subset= period ==1 )
        summary(model.fmP1)
        
        model.fmP2<-lmer(bodymass~-1 + agecat+ date.st:agecat+ (1|yearcat)+(1|ID), data=data.OAF, subset= period ==2)
        summary(model.fmP2)
        
        
    # female horn----
        model.fh1<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4, REML = F)
        
        model.fh2<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:agecat + (1|yearcat)+(1|ID), data=data.OAF,subset= agecat!=4, REML = F)
        
        model.fh3<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + date:periodcat + (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fh4<-lmer(horn~agecat + periodcat + date +
                          date:agecat + date:periodcat+ (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fh5<-lmer(horn~agecat + periodcat + date +
                          agecat:periodcat + (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fh6<-lmer(horn~agecat + periodcat + date +
                          date:periodcat+ (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fh7<-lmer(horn~agecat + periodcat + date +date:agecat +  (1|yearcat)+(1|ID), data=data.OAF,subset= agecat!=4,REML = F)
        
        model.fh8<-lmer(horn~agecat + periodcat + date + (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4, REML = F)
        
        model.fh9<-lmer(horn~agecat + periodcat + (1|yearcat)+(1|ID), data=data.OAF,subset= agecat!=4, REML = F)
        
        model.fh10<-lmer(horn~agecat + date + (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        model.fh11<-lmer(horn~agecat + (1|yearcat)+(1|ID), data=data.OAF, subset= agecat!=4,REML = F)
        
        
        AIC(model.fh1, model.fh2, model.fh3, model.fh4, model.fh5, model.fh6, model.fh7, model.fh8, model.fh9, model.fh10, model.fh11)
        
        
        date.st<-data.OAF$date-median(data.OAF$date)
        
        model.fhP1<-lmer(horn~-1 + agecat+ date.st:agecat+ (1|yearcat)+(1|ID), data=data.OAF, subset= period ==1)
        summary(model.fhP1)
        
        model.fhP2<-lmer(horn~-1 + agecat+ date.st:agecat+ (1|yearcat)+(1|ID), data=data.OAF, subset= period ==2)
        summary(model.fhP2)
        
        
        
# 4.3. Analysis of the static allometry all data ----
    # males----
        model.OAM1<-lmer(log(horn)~log(bodymass) : periodcat : agecat + 
                           log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM2<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM3<-lmer(log(horn)~log(bodymass) : periodcat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM4<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM5<-lmer(log(horn)~ log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM6<-lmer(log(horn)~log(bodymass) : periodcat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM7<-lmer(log(horn)~ log(bodymass) + periodcat : agecat +
                           periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE) 
        
        model.OAM8<-lmer(log(horn)~log(bodymass) + periodcat + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE) 
        
        model.OAM9<-lmer(log(horn)~log(bodymass) + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM10<-lmer(log(horn)~ log(bodymass) + periodcat + agecat + 
                            (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM11<-lmer(log(horn)~ log(bodymass) + agecat + 
                            (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        
        AIC(model.OAM1, model.OAM2, model.OAM3, model.OAM4, model.OAM5, model.OAM6, model.OAM7, model.OAM8, model.OAM9, model.OAM10, model.OAM11)
        AICc(model.OAM1)
        AICc(model.OAM2)
        AICc(model.OAM3)
        AICc(model.OAM4)
        AICc(model.OAM5)
        AICc(model.OAM6)
        AICc(model.OAM7)
        AICc(model.OAM8)
        AICc(model.OAM9)
        AICc(model.OAM10)
        AICc(model.OAM11)
        
        
        
        # Difference in slope period 1
        
        mod.OAM1a<-lmer(log(horn)~log(bodymass) + agecat +log(bodymass):agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAM, subset=period==1, REML=F) 
        
        mod.OAM1b<-lmer(log(horn)~log(bodymass) + agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAM, subset=period==1, REML=F) 
        
        AIC(mod.OAM1a, mod.OAM1b)
        
        # Difference in slope period 2
        
        mod.OAM2a<-lmer(log(horn)~log(bodymass) + agecat +log(bodymass):agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAM, subset=period==2, REML=F) 
        
        mod.OAM2b<-lmer(log(horn)~log(bodymass) + agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAM, subset=period==2, REML=F) 
        
        AIC(mod.OAM2a, mod.OAM2b)
        
        # Parameter estimates
        
        massCm<-residuals(lm(log(bodymass)~agecat, data = data.OAM))
        
        mod.OAM1a<-lmer(log(horn)~-1 +agecat +massCm:agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAM, subset=period==1, REML=TRUE) 
        summary(mod.OAM1a)
        
        
        mod.OAM2a<-lmer(log(horn)~-1 +agecat +massCm:agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAM, subset=period==2, REML=TRUE) 
        summary(mod.OAM2a)
        
        
    # female-----
        data.OAF2<-subset(data.OAF, data.OAF$age<3)#There is not enough females older than 2 years in the second period 
        
        
        model.OAF4<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        model.OAF6<-lmer(log(horn)~log(bodymass) : periodcat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        model.OAF8<-lmer(log(horn)~log(bodymass) + periodcat + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE) 
        
        model.OAF9<-lmer(log(horn)~log(bodymass) + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        model.OAF10<-lmer(log(horn)~ log(bodymass) + periodcat + agecat + 
                            (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        model.OAF11<-lmer(log(horn)~ log(bodymass) + agecat + 
                            (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        AIC(model.OAF4, model.OAF6, model.OAF8, model.OAF9, model.OAF10, model.OAF11)
        
        AICc(model.OAF4)
        AICc(model.OAF6)
        AICc(model.OAF8)
        AICc(model.OAF9)
        AICc(model.OAF10)
        AICc(model.OAF11)
        
        # Difference in slope period 1
        
        
        mod.OAF1a<-lmer(log(horn)~log(bodymass) + agecat +log(bodymass):agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAF, subset=period==1, REML=F) 
        
        mod.OAF1b<-lmer(log(horn)~log(bodymass) + agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAF, subset=period==1, REML=F) 
        
        AIC(mod.OAF1a, mod.OAF1b)
        
        # Difference in slope period 2
        
        mod.OAF2a<-lmer(log(horn)~log(bodymass) + agecat +log(bodymass):agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAF2, subset=period==2, REML=F) 
        
        mod.OAF2b<-lmer(log(horn)~log(bodymass) + agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAF2, subset=period==2, REML=F) 
        
        AIC(mod.OAF2a, mod.OAF2b)
        summary(mod.OAF2a)
        
        # Parameter estimates
        
        massCf<-residuals(lm(log(bodymass)~agecat, data = data.OAF))
        
        mod.OAF1a<-lmer(log(horn)~-1 +agecat +massCf:agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAF, subset=period==1, REML=TRUE) 
        summary(mod.OAF1a)
        
        
        mod.OAF2a<-lmer(log(horn)~-1 +agecat +massCf:agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAF, subset=period==2 & age<3, REML=TRUE) 
        summary(mod.OAF2a)
        
        
        
        
        
# 4.4. Analysis static allometry 10 years with most extreme density-----
        ##10 years low density 1988 - 1997 -- 10 years with high denisty 2001 2010
        ##Because the period are slighlty different regarding the transition we reload the original data file
        
    # Data selection---- 
        data.OA<-read.table("C:\\\\Oreamnos.txt", header=T)#Data with only animals of 1 year or more 
        with(data.OA, table(period, yearbirth))
        with(data.OA, table(period, year))
        
        horn<-pmax(data.OA$right_horn, data.OA$left_horn)
        data.OA<-data.frame(data.OA, horn)
        
        data.OA <- data.OA[-c(531, 41),] #two outliers
        
        data.OA1<-data.OA[data.OA$year<1997| data.OA$year>2000,]#removing observation year 1997 - 2000
        data.OA2<-data.OA1[data.OA1$year<2011,]#removing observation after year 2011 
        
        yearcat<-factor(data.OA2$year)
        agecat<-factor(data.OA2$age)
        periodcat<-factor(data.OA2$period)
        data.OA2<-cbind(data.OA2, yearcat, agecat, periodcat)
        
        data.OA<-data.frame(data.OA2)
        
        data.OAM<-subset(data.OA2, data.OA$sexe==2 & data.OA$age<6)#Males 
        data.OAF<-subset(data.OA2, data.OA$sexe==1  & data.OA$age<6)#Females
        
        
    # Males----
        model.OAM1<-lmer(log(horn)~log(bodymass) : periodcat : agecat + 
                           log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM2<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM3<-lmer(log(horn)~log(bodymass) : periodcat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM4<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM5<-lmer(log(horn)~ log(bodymass) : agecat + periodcat : agecat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM6<-lmer(log(horn)~log(bodymass) : periodcat +
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM7<-lmer(log(horn)~ log(bodymass) + periodcat : agecat +
                           periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE) 
        
        model.OAM8<-lmer(log(horn)~log(bodymass) + periodcat + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE) 
        
        model.OAM9<-lmer(log(horn)~log(bodymass) + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM10<-lmer(log(horn)~ log(bodymass) + periodcat + agecat + 
                            (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        model.OAM11<-lmer(log(horn)~ log(bodymass) + agecat + 
                            (1|yearbirth)+(1|ID), data=data.OAM, REML=FALSE)  
        
        
        AIC(model.OAM1, model.OAM2, model.OAM3, model.OAM4, model.OAM5, model.OAM6, model.OAM7, model.OAM8, model.OAM9, model.OAM10, model.OAM11)
        AICc(model.OAM1)
        AICc(model.OAM2)
        AICc(model.OAM3)
        AICc(model.OAM4)
        AICc(model.OAM5)
        AICc(model.OAM6)
        AICc(model.OAM7)
        AICc(model.OAM8)
        AICc(model.OAM9)
        AICc(model.OAM10)
        AICc(model.OAM11)
      
        # Parameter estimates
        
        massCm<-residuals(lm(log(bodymass)~agecat, data = data.OAM))
        
        mod.OAM1a<-lmer(log(horn)~-1 +agecat +massCm:agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAM, subset=period==1, REML=TRUE) 
        summary(mod.OAM1a)
        
        
        mod.OAM2a<-lmer(log(horn)~-1 +agecat +massCm:agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAM, subset=period==2, REML=TRUE) 
        summary(mod.OAM2a)
        
        
    # Females mountain goat----- 
        
        data.OAF2<-subset(data.OAF, data.OAF$age<3)#There is not enough females older than 2 years in the second period 
        
        
        model.OAF4<-lmer(log(horn)~log(bodymass) : periodcat + log(bodymass) : agecat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        model.OAF6<-lmer(log(horn)~log(bodymass) : periodcat + 
                           log(bodymass) + periodcat + agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        model.OAF8<-lmer(log(horn)~log(bodymass) + periodcat + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE) 
        
        model.OAF9<-lmer(log(horn)~log(bodymass) + agecat + 
                           log(bodymass) : agecat + 
                           (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        model.OAF10<-lmer(log(horn)~ log(bodymass) + periodcat + agecat + 
                            (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        model.OAF11<-lmer(log(horn)~ log(bodymass) + agecat + 
                            (1|yearbirth)+(1|ID), data=data.OAF2, REML=FALSE)  
        
        AIC(model.OAF4, model.OAF6, model.OAF8, model.OAF9, model.OAF10, model.OAF11)
        
        AICc(model.OAF4)
        AICc(model.OAF6)
        AICc(model.OAF8)
        AICc(model.OAF9)
        AICc(model.OAF10)
        AICc(model.OAF11)
        
 
        # Parameter estimates
        
        massCf<-residuals(lm(log(bodymass)~agecat, data = data.OAF))
        
        mod.OAF1a<-lmer(log(horn)~-1 +agecat +massCf:agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAF, subset=period==1, REML=TRUE) 
        summary(mod.OAF1a)
        
        
        mod.OAF2a<-lmer(log(horn)~-1 +agecat +massCf:agecat + 
                          (1|yearbirth)+(1|ID), data=data.OAF, subset=period==2 & age<3, REML=TRUE) 
        summary(mod.OAF2a)
        
        
        
###################################################
# 5 figure 2 ----- 
        
        
        library(lme4)
        layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = TRUE))
        
        
#Ibex---- 
    #data-----
        data.CI<-read.table("C:\\\\Ibex.txt", header=T)#Data corrected - no duplicated entree. 
        head(data.CI)
        
        Hornm<-pmax(data.CI$lengthR, data.CI$lengthL)#selecting the longest of the two horns
        year<-data.CI$yearbirth+data.CI$ageatcapture#year of capture
        diff.st<-(data.CI$diff-median(data.CI$diff))*(-1)#median date of capture
        data.CI<-data.frame(data.CI, year, diff.st, Hornm)
        data.CI<-data.CI[data.CI$yearbirth<1999 | data.CI$yearbirth>2002,]#removing transition period
        
    # Male ----
        
        data.CIM2<-subset(data.CI, data.CI$sexe=="M"& data.CI$ageatcapture<6)
        age<-factor(data.CIM2$ageatcapture) 
        
        #Models for the parameters
        model.ageMP1<-lmer(log(Hornm)~-1+ age + log(bodymass):age + 
                             (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==1, REML=TRUE)  
        
        model.ageMP2<-lmer(log(Hornm)~-1+ age +log(bodymass):age + 
                             (1|yearbirth)+(1|ID), data=data.CIM2, subset=birthperiod ==2, REML=TRUE)  
        
        #graph
        par(mar=c(1,5,5.5,4.5))
        
        with(data.CIM2, plot(bodymass, Hornm, log="xy", pch=16, axes = FALSE, 
                             col=c("black", "grey")[as.numeric(birthperiod)],
                             main="Males",cex.main=2.5, xlab="", ylab="Horn length (cm)", cex.lab = 1.8))
        
        axis(side=1, at=seq(0, 120, by=10))
        axis(side=2, at=seq(0, 100, by=10))
        
        
        #period 1
        
        xx1.1<-with(data.CIM2, seq(min(bodymass[ageatcapture==1 & birthperiod == 1]), max(bodymass[ageatcapture==1& birthperiod == 1]), 0.1))
        yy1.1<-exp(fixef(model.ageMP1)[1])*xx1.1^fixef(model.ageMP1)[6]
        lines(xx1.1,yy1.1, col="blue", lwd=1)
        
        text(30, 16, pos=4, label = "1", cex=1.5)
        
        xx1.2<-with(data.CIM2, seq(min(bodymass[ageatcapture==2 & birthperiod == 1]), max(bodymass[ageatcapture==2& birthperiod == 1]), 0.1))
        yy1.2<-exp(fixef(model.ageMP1)[2])*xx1.2^fixef(model.ageMP1)[7]
        lines(xx1.2,yy1.2, col="blue", lwd=1)
        
        text(45,30, pos=4, label = "2", cex=1.5)
        
        xx1.3<-with(data.CIM2, seq(min(bodymass[ageatcapture==3 & birthperiod == 1]), max(bodymass[ageatcapture==3& birthperiod == 1]), 0.1))
        yy1.3<-exp(fixef(model.ageMP1)[3])*xx1.3^fixef(model.ageMP1)[8]
        lines(xx1.3,yy1.3, col="blue", lwd=1)
        
        text(58, 42, pos=4, label = "3", cex=1.5)
        
        
        xx1.4<-with(data.CIM2, seq(min(bodymass[ageatcapture==4 & birthperiod == 1]), max(bodymass[ageatcapture==4& birthperiod == 1]), 0.1))
        yy1.4<-exp(fixef(model.ageMP1)[4])*xx1.4^fixef(model.ageMP1)[9]
        lines(xx1.4,yy1.4, col="blue", lwd=1)
        
        text(72, 50, pos=4, label = "4", cex=1.5)
        
        xx1.5<-with(data.CIM2, seq(min(bodymass[ageatcapture==5 & birthperiod == 1]), max(bodymass[ageatcapture==5& birthperiod == 1]), 0.1))
        yy1.5<-exp(fixef(model.ageMP1)[5])*xx1.5^fixef(model.ageMP1)[10]
        lines(xx1.5,yy1.5, col="blue", lwd=1)
        
        text(89, 68, pos=4, label = "5", cex=1.5)
        
        
        #period 2
        
        xx2.1<-with(data.CIM2, seq(min(bodymass[ageatcapture==1& birthperiod == 2]), max(bodymass[ageatcapture==1& birthperiod == 2]), 0.1))
        yy2.1<-exp(fixef(model.ageMP2)[1])*xx2.1^fixef(model.ageMP2)[6]
        lines(xx2.1,yy2.1, col="red", lwd=1)
        
        xx2.2<-with(data.CIM2, seq(min(bodymass[ageatcapture==2& birthperiod == 2]), max(bodymass[ageatcapture==2& birthperiod == 2]), 0.1))
        yy2.2<-exp(fixef(model.ageMP2)[2])*xx2.2^fixef(model.ageMP2)[7]
        lines(xx2.2,yy2.2, col="red", lwd=1)
        
        xx2.3<-with(data.CIM2, seq(min(bodymass[ageatcapture==3& birthperiod == 2]), max(bodymass[ageatcapture==3& birthperiod == 2]), 0.1))
        yy2.3<-exp(fixef(model.ageMP2)[3])*xx2.3^fixef(model.ageMP2)[8]
        lines(xx2.3,yy2.3, col="red", lwd=1)
        
        xx2.4<-with(data.CIM2, seq(min(bodymass[ageatcapture==4& birthperiod == 2]), max(bodymass[ageatcapture==4& birthperiod == 2]), 0.1))
        yy2.4<-exp(fixef(model.ageMP2)[4])*xx2.4^fixef(model.ageMP2)[9]
        lines(xx2.4,yy2.4, col="red", lwd=1)
        
        xx2.5<-with(data.CIM2, seq(min(bodymass[ageatcapture==5& birthperiod == 2]), max(bodymass[ageatcapture==5& birthperiod == 2]), 0.1))
        yy2.5<-exp(fixef(model.ageMP2)[5])*xx2.5^fixef(model.ageMP2)[10]
        lines(xx2.5,yy2.5, col="red", lwd=1)
        
        
    # female ----
        # Data
        data.CIF2<-subset(data.CI, data.CI$sexe=="F"& data.CI$ageatcapture<4)
        age<-factor(data.CIF2$ageatcapture) 
        
        #model
        model.ageFP1<-lmer(log(Hornm)~-1+ age +log(bodymass):age + 
                             (1|yearbirth), data=data.CIF2, subset=birthperiod ==1, REML=TRUE)  
        
        
        
        model.ageFP2<-lmer(log(Hornm)~-1+log(bodymass) + age +
                             (1|yearbirth), data=data.CIF2, subset=birthperiod ==2, REML=TRUE)  
        
    #Graph
        par(mar=c(1,4,5.5,4.5))
        
        with(data.CIF2, plot(bodymass, Hornm, log="xy", pch=16, axes=FALSE,
                             col=c("black", "grey")[as.numeric(birthperiod)],
                             xlab="", ylab="", cex.lab = 1.8, main ="Females", cex.main=2.5))
        
        axis(side=1, at=seq(0, 120, by=5))
        axis(side=2, at=seq(0, 100, by=5))
        
        # period 1
        
        xx1.1<-with(data.CIF2, seq(min(bodymass[ageatcapture==1& birthperiod == 1]), max(bodymass[ageatcapture==1& birthperiod == 1]), 0.1))
        yy1.1<-exp(fixef(model.ageFP1)[1])*xx1.1^fixef(model.ageFP1)[4]
        lines(xx1.1,yy1.1, col="blue")
        
        text(25, 11, pos=4, label = "1", cex=1.5)
        
        xx1.2<-with(data.CIF2, seq(min(bodymass[ageatcapture==2& birthperiod == 1]), max(bodymass[ageatcapture==2 & birthperiod == 1]), 0.1))
        yy1.2<-exp(fixef(model.ageFP1)[2])*xx1.2^fixef(model.ageFP1)[5]
        lines(xx1.2,yy1.2, col="blue")
        
        text(40, 17, pos=4, label = "2", cex=1.5)
        
        xx1.3<-with(data.CIF2, seq(min(bodymass[ageatcapture==3& birthperiod == 1]), max(bodymass[ageatcapture==3& birthperiod == 1]), 0.1))
        yy1.3<-exp(fixef(model.ageFP1)[3])*xx1.3^fixef(model.ageFP1)[6]
        lines(xx1.3,yy1.3, col="blue")
        
        text(42, 23, pos=4, label = "3", cex=1.5)
        
        # period 2
        
        xx2.1<-with(data.CIF2, seq(min(bodymass[ageatcapture==1& birthperiod == 2]), max(bodymass[ageatcapture==1& birthperiod == 2]), 0.1))
        yy2.1<-exp(fixef(model.ageFP2)[2])*xx2.1^(fixef(model.ageFP2)[1])
        lines(xx2.1,yy2.1, col="red")
        
        xx2.2<-with(data.CIF2, seq(min(bodymass[ageatcapture==2& birthperiod == 2]), max(bodymass[ageatcapture==2 & birthperiod == 2]), 0.1))
        yy2.2<-exp(fixef(model.ageFP2)[3])*xx2.2^fixef(model.ageFP2)[1]
        lines(xx2.2,yy2.2, col="red")
        
        xx2.3<-with(data.CIF2, seq(min(bodymass[ageatcapture==3& birthperiod == 2]), max(bodymass[ageatcapture==3& birthperiod == 2]), 0.1))
        yy2.3<-exp(fixef(model.ageFP2)[4])*xx2.3^fixef(model.ageFP2)[1]
        lines(xx2.3,yy2.3, col="red")
        
        
# Ovis canadensis----
    # data----
        data.BHM<-read.table("C:\\\\Ovis_M.txt", header=T)   
        horn<-pmax(data.BHM$horn.r, data.BHM$horn.l)
        data.BHM<-data.frame(data.BHM, horn)
        
        data.BHM<-data.BHM[-c(30),]  
        
        yearcat<-factor(data.BHM$year)
        agecat<-factor(data.BHM$age)
        periodcat<-factor(data.BHM$period)
        data.BHM<-cbind(data.BHM, yearcat, agecat, periodcat)
        data.BHM<-data.frame(data.BHM)
        
        data.BHM1<-subset(data.BHM, data.BHM$period== 1)
        data.BHM2<-subset(data.BHM, data.BHM$period== 2)
        
        data.BHF<-read.table("C:\\\\Ovis_F.txt", header=T)  
        data.BHF<-subset(data.BHF, data.BHF$age<6)
        
        yearcat<-factor(data.BHF$year)
        agecat<-factor(data.BHF$age)
        periodcat<-factor(data.BHF$period)
        data.BHF<-cbind(data.BHF, yearcat, agecat, periodcat)
        data.BHF<-data.frame(data.BHF) 
        
        data.BHF1<-subset(data.BHF, data.BHF$period== 1)
        data.BHF2<-subset(data.BHF, data.BHF$period== 2)
        
        
    # Male----
        #Models 
        model.BHM1a<-lmer(log(horn)~ -1+  agecat + log(bodymass): agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHM1, REML=FALSE)  
        model.BHM2a<-lmer(log(horn)~ -1+  agecat + log(bodymass): agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHM2, REML=FALSE)  
        #graph
        par(mar=c(1,5,2,4.5))
        with(data.BHM, plot(bodymass, horn, log="xy", pch=16, axes=FALSE,
                            col=c("black", "grey")[as.numeric(periodcat)], 
                            xlab="", ylab="Horn length (cm)", cex.lab=1.8))
        
        axis(side=1, at=seq(0, 120, by=10))
        axis(side=2, at=seq(0, 100, by=10))
        
        #Static allometry period 1
        xx1.2<-with(data.BHM, seq(min(bodymass[agecat==2& periodcat == 1]), max(bodymass[agecat==2& periodcat == 1]), 0.1))
        yy1.2<-exp(fixef(model.BHM1a)[1])*xx1.2^fixef(model.BHM1a)[5]
        lines(xx1.2,yy1.2, col="blue", lwd=1)
        
        text(33, 23, pos=4, label = "1", cex=1.5)
        
        xx1.3<-with(data.BHM, seq(min(bodymass[agecat==3& periodcat == 1]), max(bodymass[agecat==3& periodcat == 1]), 0.1))
        yy1.3<-exp(fixef(model.BHM1a)[2])*xx1.3^fixef(model.BHM1a)[6]
        lines(xx1.3,yy1.3, col="blue", lwd=1)
        
        text(45, 40, pos=4, label = "2", cex=1.5)
        
        xx1.4<-with(data.BHM, seq(min(bodymass[agecat==4& periodcat == 1]), max(bodymass[agecat==4& periodcat == 1]), 0.1))
        yy1.4<-exp(fixef(model.BHM1a)[3])*xx1.4^fixef(model.BHM1a)[7]
        lines(xx1.4,yy1.4, col="blue", lwd=1)
        
        text(55, 54, pos=4, label = "3", cex=1.5)
        
        xx1.5<-with(data.BHM, seq(min(bodymass[agecat==5& periodcat == 1]), max(bodymass[agecat==5& periodcat == 1]), 0.1))
        yy1.5<-exp(fixef(model.BHM1a)[4])*xx1.5^fixef(model.BHM1a)[8]
        lines(xx1.5,yy1.5, col="blue", lwd=1)
        
        text(65, 68, pos=4, label = "4", cex=1.5)
        
        #Static allometry period 2 
        xx2.2<-with(data.BHM2, seq(min(bodymass[agecat==2& periodcat == 2]), max(bodymass[agecat==2& periodcat == 2]), 0.1))
        yy2.2<-exp(fixef(model.BHM2a)[1])*xx2.2^fixef(model.BHM2a)[5]
        lines(xx2.2,yy2.2, col="red", lwd=1)
        
        xx2.3<-with(data.BHM2, seq(min(bodymass[agecat==3& periodcat == 2]), max(bodymass[agecat==3& periodcat == 2]), 0.1))
        yy2.3<-exp(fixef(model.BHM2a)[2])*xx2.3^fixef(model.BHM2a)[6]
        lines(xx2.3,yy2.3, col="red", lwd=1)
        
        xx2.4<-with(data.BHM2, seq(min(bodymass[agecat==4& periodcat == 2]), max(bodymass[agecat==4& periodcat == 2]), 0.1))
        yy2.4<-exp(fixef(model.BHM2a)[3])*xx2.4^fixef(model.BHM2a)[7]
        lines(xx2.4,yy2.4, col="red", lwd=1)
        
        xx2.5<-with(data.BHM2, seq(min(bodymass[agecat==5& periodcat == 2]), max(bodymass[agecat==5& periodcat == 2]), 0.1))
        yy2.5<-exp(fixef(model.BHM2a)[4])*xx2.5^fixef(model.BHM2a)[8]
        lines(xx2.5,yy2.5, col="red", lwd=1)
        
    # Female ----
        
        #Models 
        model.BHF1a<-lmer(log(horn)~ -1+  agecat + log(bodymass): agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHF1, REML=FALSE)  
        model.BHF2a<-lmer(log(horn)~ -1+  agecat + log(bodymass): agecat + 
                            (1|yearbirth)+(1|ID), data=data.BHF2, REML=FALSE)  
        
        par(mar=c(1,4 ,2 , 4.5))
        
        #Graph
        with(data.BHF, plot(bodymass, horn, log="xy", pch=16, axes=FALSE,
                            col=c("black", "grey")[as.numeric(periodcat)],
                            xlab="", ylab="", cex.lab=1.8))
        
        axis(side=1, at=seq(0, 120, by=5))
        axis(side=2, at=seq(0, 100, by=5))
        
        #Static allometry period 1
        xx1.2<-with(data.BHF, seq(min(bodymass[agecat==2& periodcat == 1]), max(bodymass[agecat==2& periodcat == 1]), 0.1))
        yy1.2<-exp(fixef(model.BHF1a)[1])*xx1.2^fixef(model.BHF1a)[5]
        lines(xx1.2,yy1.2, col="blue")
        
        text(25, 11, pos=4, label = "1", cex=1.5)
        
        xx1.3<-with(data.BHF, seq(min(bodymass[agecat==3& periodcat == 1]), max(bodymass[agecat==3& periodcat == 1]), 0.1))
        yy1.3<-exp(fixef(model.BHF1a)[2])*xx1.3^fixef(model.BHF1a)[6]
        lines(xx1.3,yy1.3, col="blue")
        
        text(30, 18, pos=4, label = "2", cex=1.5)
        
        xx1.4<-with(data.BHF, seq(min(bodymass[agecat==4& periodcat == 1]), max(bodymass[agecat==4& periodcat == 1]), 0.1))
        yy1.4<-exp(fixef(model.BHF1a)[3])*xx1.4^fixef(model.BHF1a)[7]
        lines(xx1.4,yy1.4, col="blue")
        
        text(37, 21, pos=4, label = "3", cex=1.5)
        
        xx1.5<-with(data.BHF, seq(min(bodymass[agecat==5& periodcat == 1]), max(bodymass[agecat==5& periodcat == 1]), 0.1))
        yy1.5<-exp(fixef(model.BHF1a)[4])*xx1.5^fixef(model.BHF1a)[8]
        lines(xx1.5,yy1.5, col="blue")
        
        text(40, 23, pos=4, label = "4", cex=1.5)
        
        #Static allometry period 2 
        xx2.2<-with(data.BHF, seq(min(bodymass[agecat==2& periodcat == 2]), max(bodymass[agecat==2& periodcat == 2]), 0.1))
        yy2.2<-exp(fixef(model.BHF2a)[1])*xx2.2^fixef(model.BHF2a)[5]
        lines(xx2.2,yy2.2, col="red")
        
        xx2.3<-with(data.BHF, seq(min(bodymass[agecat==3& periodcat == 2]), max(bodymass[agecat==3& periodcat == 2]), 0.1))
        yy2.3<-exp(fixef(model.BHF2a)[2])*xx2.3^fixef(model.BHF2a)[6]
        lines(xx2.3,yy2.3, col="red")
        
        xx2.4<-with(data.BHF, seq(min(bodymass[agecat==4& periodcat == 2]), max(bodymass[agecat==4& periodcat == 2]), 0.1))
        yy2.4<-exp(fixef(model.BHF2a)[3])*xx2.4^fixef(model.BHF2a)[7]
        lines(xx2.4,yy2.4, col="red")
        
        xx2.5<-with(data.BHF, seq(min(bodymass[agecat==5& periodcat == 2]), max(bodymass[agecat==5& periodcat == 2]), 0.1))
        yy2.5<-exp(fixef(model.BHF2a)[4])*xx2.5^fixef(model.BHF2a)[8]
        lines(xx2.5,yy2.5, col="red")
        
        
        
        
        
        
# male Rupicapra -----
    # data-----     
        data.rr2<-read.table("C:\\\\chamois.txt", header=T)
        head(data.rr2)
        
        yr.birth<-with(data.rr2, yr-age)
        data.rr2<-data.frame(data.rr2, yr.birth)
        
        data.rr2<-data.rr2[data.rr2$yr<1998 | data.rr2$yr>2000,]#removing individuals in the transition period
        data.rr2<-subset(data.rr2, data.rr2$age>0)#removing yearling 
        data.rr2<-subset(data.rr2, data.rr2$age<9)#removing individuals older than 8 year (not enough in each age class)
        data.rr2<-subset(data.rr2, yr.birth>1981)#we remove all individuals born during the episode of Kerato
        
        horn<-pmax(data.rr2$hornR, data.rr2$hornL)/10#selecting the longest of the two horns
        period<-with(data.rr2, ifelse(density=="low", 1,2))#defining the two period until 1998 / after 1998
        dayk<-with(data.rr2, daynr-median(daynr))# Centering the hunting date around the median hunting date
        dayk2<-dayk^2
        data.rr2<-data.frame(data.rr2, horn, period, dayk, dayk2)
        
        data.rrM<-subset(data.rr2, data.rr2$sex=="M")
        
        data.rrF<-subset(data.rr2, data.rr2$sex=="F")
        
        
    # data correction----      
        #male
        mod.rr_mm2b<-lmer(mass~ -1 + factor(age) + dayk + dayk2+ dayk:age + dayk:period + dayk2:period +
                            (1|season), data = data.rrM, REML=T)#test for the effect of hunting date 
        summary(mod.rr_mm2b)
        
        
        masscorr<-data.rrM$mass-((fixef(mod.rr_mm2b)[9]+fixef(mod.rr_mm2b)[12]*data.rrM$period+fixef(mod.rr_mm2b)[11]*data.rrM$age)*data.rrM$dayk+
                                   (fixef(mod.rr_mm2b)[10]+fixef(mod.rr_mm2b)[13]*data.rrM$period)*data.rrM$dayk2)#correcting the mass for hunting date
        
        data.rrM<-data.frame(data.rrM, masscorr)
        #female
        mod.rr_Fm3b<-lmer(mass~ -1 + factor(age) + dayk + dayk2+ dayk:period + dayk2:period +
                            (1|season), data = data.rrF, REML=T)#test for the effect of hunting date 
        summary(mod.rr_Fm3b)
        
        
        masscorr<-data.rrF$mass-((fixef(mod.rr_Fm3b)[9]+fixef(mod.rr_Fm3b)[11]*data.rrF$period)*data.rrF$dayk+
                                   (fixef(mod.rr_Fm3b)[10]+fixef(mod.rr_Fm3b)[12]*data.rrF$period)*data.rrF$dayk2)#correcting the mass for hunting date
        
        data.rrF<-data.frame(data.rrF, masscorr)
        
    # Males-----     
        age2<-factor(data.rrM$age)
        data.rrM<-data.frame(data.rrM, age2)
        head(data.rrM)
        #model period 1-
        
        mod.rrMg1<-lmer(log(horn)~-1 +age2 +log(masscorr):age2 + 
                          (1|cohort), data=data.rrM, subset=period==1, REML=TRUE) 
        
        summary(mod.rrMg1)
        
        #model period 2-
        mod.rrMg2<-lmer(log(horn)~-1 +age2 + log(masscorr):age2 + 
                          (1|cohort), data=data.rrM, subset=period==2, REML=TRUE) 
        
        summary(mod.rrMg2)
        
        #Graph
        par(mar=c(1,5,2,4.5))
        with(data.rrM, plot(masscorr, horn, log="xy", pch=16, axes=FALSE,
                            col=c("black", "grey")[as.numeric(period)],
                            xlab="", ylab="Horn length (cm)", cex.lab=1.8))
        
        axis(side=1, at=seq(0, 120, by=2))
        axis(side=2, at=seq(0, 30, by=2))
        
        #period 1
        xx1.1<-with(data.rrM, seq(min(masscorr[age==1& period == 1]), max(masscorr[age==1& period == 1]), 0.1))
        yy1.1<-exp(fixef(mod.rrMg1)[1])*xx1.1^fixef(mod.rrMg1)[9]
        lines(xx1.1,yy1.1, col="blue")
        
        text(9, 13, pos=4, label = "1", cex=1.5)
        
        xx1.2<-with(data.rrM, seq(min(masscorr[age==2& period == 1]), max(masscorr[age==2& period == 1]), 0.1))
        yy1.2<-exp(fixef(mod.rrMg1)[2])*xx1.2^fixef(mod.rrMg1)[10]
        lines(xx1.2,yy1.2, col="blue")
        
        text(12, 16, pos=4, label = "2", cex=1.5)
        
        xx1.3<-with(data.rrM, seq(min(masscorr[age==3& period == 1]), max(masscorr[age==3& period == 1]), 0.1))
        yy1.3<-exp(fixef(mod.rrMg1)[3])*xx1.3^fixef(mod.rrMg1)[11]
        lines(xx1.3,yy1.3, col="blue")
        
        text(15, 20, pos=4, label = "3", cex=1.5)
        
        #period 2
        xx2.1<-with(data.rrM, seq(min(masscorr[age==1& period == 2]), max(masscorr[age==1& period == 2]), 0.1))
        yy2.1<-exp(fixef(mod.rrMg2)[1])*xx2.1^fixef(mod.rrMg2)[9]
        lines(xx2.1,yy2.1, col="red")
        
        xx2.2<-with(data.rrM, seq(min(masscorr[age==2& period == 2]), max(masscorr[age==2& period == 2]), 0.1))
        yy2.2<-exp(fixef(mod.rrMg2)[2])*xx2.2^fixef(mod.rrMg2)[10]
        lines(xx2.2,yy2.2, col="red")
        
        xx2.3<-with(data.rrM, seq(min(masscorr[age==3& period == 2]), max(masscorr[age==3& period == 2]), 0.1))
        yy2.3<-exp(fixef(mod.rrMg2)[3])*xx2.3^fixef(mod.rrMg2)[11]
        lines(xx2.3,yy2.3, col="red")
        
    # Female ----
        age2<-factor(data.rrF$age)
        data.rrF<-data.frame(data.rrF, age2)
        
        #model period 1-
        mod.rrF1<-lmer(log(horn)~-1 +age2 +log(masscorr):age2 + 
                         (1|cohort), data=data.rrF, subset=period==1, REML=TRUE) 
        
        #model period 2-
        mod.rrF2<-lmer(log(horn)~-1 +age2 + log(masscorr):age2 + 
                         (1|cohort), data=data.rrF, subset=period==2, REML=TRUE) 
        #Graph
        
        par(mar=c(1,4,2, 4.5))
        with(data.rrF, plot(masscorr, horn, log="xy", pch=16, axes=FALSE,
                            col=c("black", "grey")[as.numeric(period)],
                            xlab="", ylab="", cex.lab=1.8))
        
        axis(side=1, at=seq(0, 40, by=5))
        axis(side=2, at=seq(0, 32, by=2))
        
        # period 1
        xx1.1<-with(data.rrF, seq(min(masscorr[age==1& period == 1]), max(masscorr[age==1& period == 1]), 0.1))
        yy1.1<-exp(fixef(mod.rrF1)[1])*xx1.1^fixef(mod.rrF1)[9]
        lines(xx1.1,yy1.1, col="blue")
        
        text(9, 11, pos=4, label = "1", cex=1.5)
        
        xx1.2<-with(data.rrF, seq(min(masscorr[age==2& period == 1]), max(masscorr[age==2& period == 1]), 0.1))
        yy1.2<-exp(fixef(mod.rrF1)[2])*xx1.2^fixef(mod.rrF1)[10]
        lines(xx1.2,yy1.2, col="blue")
        
        text(14, 16, pos=4, label = "2", cex=1.5)
        
        xx1.3<-with(data.rrF, seq(min(masscorr[age==3& period == 1]), max(masscorr[age==3& period == 1]), 0.1))
        yy1.3<-exp(fixef(mod.rrF1)[3])*xx1.3^fixef(mod.rrF1)[11]
        lines(xx1.3,yy1.3, col="blue")
        
        text(30, 18, pos=4, label = "3", cex=1.5)
        
        xx1.4<-with(data.rrF, seq(min(masscorr[age==4& period == 1]), max(masscorr[age==4& period == 1]), 0.1))
        yy1.4<-exp(fixef(mod.rrF1)[4])*xx1.4^fixef(mod.rrF1)[12]
        lines(xx1.4,yy1.4, col="blue")
        
        text(33, 19, pos=4, label = "4", cex=1.5)
        
        
        #period 2
        xx2.1<-with(data.rrF, seq(min(masscorr[age==1& period == 2]), max(masscorr[age==1& period == 2]), 0.1))
        yy2.1<-exp(fixef(mod.rrF2)[1])*xx2.1^fixef(mod.rrF2)[9]
        lines(xx2.1,yy2.1, col="red")
        
        xx2.2<-with(data.rrF, seq(min(masscorr[age==2& period == 2]), max(masscorr[age==2& period == 2]), 0.1))
        yy2.2<-exp(fixef(mod.rrF2)[2])*xx2.2^fixef(mod.rrF2)[10]
        lines(xx2.2,yy2.2, col="red")
        
        xx2.3<-with(data.rrF, seq(min(masscorr[age==3& period == 2]), max(masscorr[age==3& period == 2]), 0.1))
        yy2.3<-exp(fixef(mod.rrF2)[3])*xx2.3^fixef(mod.rrF2)[11]
        lines(xx2.3,yy2.3, col="red")
        
        xx2.4<-with(data.rrF, seq(min(masscorr[age==4& period == 2]), max(masscorr[age==4& period == 2]), 0.1))
        yy2.4<-exp(fixef(mod.rrF2)[4])*xx2.4^fixef(mod.rrF2)[12]
        lines(xx2.4,yy2.4, col="red")
        
        
# male Oreamnos ----
    # data -----
        data.OA<-read.table("C:\\\\Oreamnos.txt", header=T)#Data with only animals of 1 year or more 
        head(data.OA)
        horn<-pmax(data.OA$right_horn, data.OA$left_horn)/10
        data.OA<-data.frame(data.OA, horn)
        
        data.OA <- data.OA[-c(531, 41),] #two outliers
        data.OA1<-data.OA[data.OA$year<1997| data.OA$year>2000,]#removing observation year 1997 - 2000
        
        ##Period 1 = low desnity;  2 = High density
        
        yearcat<-factor(data.OA$year)
        agecat<-factor(data.OA$age)
        periodcat<-factor(data.OA$period)
        data.OA<-cbind(data.OA, yearcat, agecat, periodcat)
        
        data.OA<-data.frame(data.OA)
        
        data.OAM<-subset(data.OA, data.OA$sexe==2 & data.OA$age<6)#Males 
        data.OAF<-subset(data.OA, data.OA$sexe==1  & data.OA$age<6)#Females
        
    # Male----
        mod.OAM1<-lmer(log(horn)~-1 +agecat +log(bodymass):agecat + 
                         (1|yearbirth)+(1|ID), data=data.OAM, subset=period==1, REML=TRUE) 
        
        mod.OAM2<-lmer(log(horn)~-1 +agecat +log(bodymass):agecat + 
                         (1|yearbirth)+(1|ID), data=data.OAM, subset=period==2, REML=TRUE) 
        
        #graph
        par(mar=c(5,5,2,4.5))
        with(data.OAM, plot(bodymass, horn, log="xy", pch=16, axes=FALSE,
                            col=c("black", "grey")[as.numeric(periodcat)],
                            xlab="body mass(kg)", ylab="Horn length (cm)", cex.lab=1.8))
        
        axis(side=1, at=seq(0, 120, by=10))
        axis(side=2, at=seq(0, 28, by=2))
        
        #period 1
        xx1.1<-with(data.OAM, seq(min(bodymass[agecat==1& periodcat == 1]), max(bodymass[agecat==1& periodcat == 1]), 0.1))
        yy1.1<-exp(fixef(mod.OAM1)[1])*xx1.1^fixef(mod.OAM1)[6]
        lines(xx1.1,yy1.1, col="blue")
        
        text(20, 10, pos=4, label = "1", cex=1.5)
        
        xx1.2<-with(data.OAM, seq(min(bodymass[agecat==2& periodcat == 1]), max(bodymass[agecat==2& periodcat == 1]), 0.1))
        yy1.2<-exp(fixef(mod.OAM1)[2])*xx1.2^fixef(mod.OAM1)[7]
        lines(xx1.2,yy1.2, col="blue")
        
        text(29, 18, pos=4, label = "2", cex=1.5)
        
        xx1.3<-with(data.OAM, seq(min(bodymass[agecat==3& periodcat == 1]), max(bodymass[agecat==3& periodcat == 1]), 0.1))
        yy1.3<-exp(fixef(mod.OAM1)[3])*xx1.3^fixef(mod.OAM1)[8]
        lines(xx1.3,yy1.3, col="blue")
        
        text(37, 22, pos=4, label = "3", cex=1.5)
        
        xx1.4<-with(data.OAM, seq(min(bodymass[agecat==4& periodcat == 1]), max(bodymass[agecat==4& periodcat == 1]), 0.1))
        yy1.4<-exp(fixef(mod.OAM1)[4])*xx1.4^fixef(mod.OAM1)[9]
        lines(xx1.4,yy1.4, col="blue")
        
        text(50, 24, pos=4, label = "4", cex=1.5)
        
        xx1.5<-with(data.OAM, seq(min(bodymass[agecat==5& periodcat == 1]), max(bodymass[agecat==5& periodcat == 1]), 0.1))
        yy1.5<-exp(fixef(mod.OAM1)[5])*xx1.5^fixef(mod.OAM1)[10]
        lines(xx1.5,yy1.5, col="blue")
        
        text(96, 25, pos=4, label = "5", cex=1.5)
        
        #period 2
        xx2.1<-with(data.OAM, seq(min(bodymass[agecat==1& periodcat == 2]), max(bodymass[agecat==1& periodcat == 2]), 0.1))
        yy2.1<-exp(fixef(mod.OAM2)[1])*xx2.1^fixef(mod.OAM2)[6]
        lines(xx2.1,yy2.1, col="red")
        
        xx2.2<-with(data.OAM, seq(min(bodymass[agecat==2& periodcat == 2]), max(bodymass[agecat==2& periodcat == 2]), 0.1))
        yy2.2<-exp(fixef(mod.OAM2)[2])*xx2.2^fixef(mod.OAM2)[7]
        lines(xx2.2,yy2.2, col="red")
        
        xx2.3<-with(data.OAM, seq(min(bodymass[agecat==3& periodcat == 2]), max(bodymass[agecat==3& periodcat == 2]), 0.1))
        yy2.3<-exp(fixef(mod.OAM2)[3])*xx2.3^fixef(mod.OAM2)[8]
        lines(xx2.3,yy2.3, col="red")
        
        xx2.4<-with(data.OAM, seq(min(bodymass[agecat==4& periodcat == 2]), max(bodymass[agecat==4& periodcat == 2]), 0.1))
        yy2.4<-exp(fixef(mod.OAM2)[4])*xx2.4^fixef(mod.OAM2)[9]
        lines(xx2.4,yy2.4, col="red")
        
        xx2.5<-with(data.OAM, seq(min(bodymass[agecat==5& periodcat == 2]), max(bodymass[agecat==5& periodcat == 2]), 0.1))
        yy2.5<-exp(fixef(mod.OAM2)[5])*xx2.5^fixef(mod.OAM2)[10]
        lines(xx2.5,yy2.5, col="red")
        
    # Female ----
        
        mod.OAF1<-lmer(log(horn)~-1 +agecat +log(bodymass):agecat + 
                         (1|yearbirth)+(1|ID), data=data.OAF, subset=period==1, REML=TRUE) 
        
        mod.OAF2<-lmer(log(horn)~-1 +agecat +log(bodymass):agecat + 
                         (1|yearbirth)+(1|ID), data=data.OAF, subset=period==2 & age<3, REML=TRUE) 
        
        par(mar=c(5,4,2,4.5))
        
        
        with(data.OAF, plot(bodymass, horn, log="xy", pch=16, axes=FALSE,
                            col=c("black", "grey")[as.numeric(periodcat)],
                            xlab="body mass(kg)", ylab="", cex.lab=1.8))
        
        axis(side=1, at=seq(0, 120, by=5))
        axis(side=2, at=seq(0, 28, by=2))
        
        #period 1
        xx1.1<-with(data.OAF, seq(min(bodymass[agecat==1& periodcat == 1]), max(bodymass[agecat==1& periodcat == 1]), 0.1))
        yy1.1<-exp(fixef(mod.OAF1)[1])*xx1.1^fixef(mod.OAF1)[6]
        lines(xx1.1,yy1.1, col="blue")
        
        text(19, 9, pos=4, label = "1", cex=1.5)
        
        xx1.2<-with(data.OAF, seq(min(bodymass[agecat==2& periodcat == 1]), max(bodymass[agecat==2& periodcat == 1]), 0.1))
        yy1.2<-exp(fixef(mod.OAF1)[2])*xx1.2^fixef(mod.OAF1)[7]
        lines(xx1.2,yy1.2, col="blue")
        
        text(26, 16, pos=4, label = "2", cex=1.5)
        
        xx1.3<-with(data.OAF, seq(min(bodymass[agecat==3& periodcat == 1]), max(bodymass[agecat==3& periodcat == 1]), 0.1))
        yy1.3<-exp(fixef(mod.OAF1)[3])*xx1.3^fixef(mod.OAF1)[8]
        lines(xx1.3,yy1.3, col="blue")
        
        text(45, 21, pos=4, label = "3", cex=1.5)
        
        xx1.4<-with(data.OAF, seq(min(bodymass[agecat==4& periodcat == 1]), max(bodymass[agecat==4& periodcat == 1]), 0.1))
        yy1.4<-exp(fixef(mod.OAF1)[4])*xx1.4^fixef(mod.OAF1)[9]
        lines(xx1.4,yy1.4, col="blue")
        
        text(80, 21.5, pos=4, label = "4", cex=1.5)
        
        xx1.5<-with(data.OAF, seq(min(bodymass[agecat==5& periodcat == 1]), max(bodymass[agecat==5& periodcat == 1]), 0.1))
        yy1.5<-exp(fixef(mod.OAF1)[5])*xx1.5^fixef(mod.OAF1)[10]
        lines(xx1.5,yy1.5, col="blue")
        
        text(77, 24, pos=4, label = "5", cex=1.5)
        
        #period 2 
        xx2.1<-with(data.OAF, seq(min(bodymass[agecat==1& periodcat == 2]), max(bodymass[agecat==1& periodcat == 2]), 0.1))
        yy2.1<-exp(fixef(mod.OAF2)[1])*xx2.1^fixef(mod.OAF2)[3]
        lines(xx2.1,yy2.1, col="red")
        
        xx2.2<-with(data.OAF, seq(min(bodymass[agecat==2& periodcat == 2]), max(bodymass[agecat==2& periodcat == 2]), 0.1))
        yy2.2<-exp(fixef(mod.OAF2)[2])*xx2.2^fixef(mod.OAF2)[4]
        lines(xx2.2,yy2.2, col="red")
        
        
        #Saved 15 x 12 inches 
        
####################################################
        
